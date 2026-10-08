"""Download KBonus-BKG catalogs and light curves from MAST.

The download is resumable: a file counts as done when it exists under its
final name.  Each file is written to ``<name>.part``, checked against the
server's Content-Length and then renamed, so an interrupted run never leaves a
truncated file behind.  Re-run the same command to continue.
"""

import argparse
import os
import queue
import sys
import threading
import time
from pathlib import Path

import requests

from .catalog import CATALOGS, catalog_path, catalog_url, lc_filenames
from .paths import ARCHIVE_BASE, data_dir, lc_filename, lc_url, local_lc_path

CHUNK = 1 << 20
TIMEOUT = (30, 120)  # connect, read (seconds)


def fetch(session, url, dest, retries=5):
    """Download ``url`` to ``dest`` atomically.  Returns the number of bytes."""
    dest = Path(dest)
    dest.parent.mkdir(parents=True, exist_ok=True)
    part = dest.with_name(dest.name + ".part")
    for attempt in range(retries):
        try:
            with session.get(url, stream=True, timeout=TIMEOUT) as response:
                if response.status_code == 404:
                    raise FileNotFoundError(url)
                response.raise_for_status()
                expected = response.headers.get("Content-Length")
                size = 0
                with open(part, "wb") as out:
                    for chunk in response.iter_content(CHUNK):
                        out.write(chunk)
                        size += len(chunk)
                if expected is not None and size != int(expected):
                    raise OSError(f"got {size} of {expected} bytes")
            os.replace(part, dest)
            return size
        except FileNotFoundError:
            raise
        except (requests.RequestException, OSError):
            if attempt == retries - 1:
                part.unlink(missing_ok=True)
                raise
            time.sleep(min(2**attempt * 2, 60))


def downloaded_lc_names(root=None):
    """Return the names of all light curve files present under ``lcs/``."""
    names = set()
    lcs = data_dir(root) / "lcs"
    if not lcs.is_dir():
        return names
    for directory in os.scandir(lcs):
        if directory.is_dir():
            names.update(e.name for e in os.scandir(directory.path) if e.name.endswith("_lc.fits"))
    return names


def download_catalogs(root=None, base=ARCHIVE_BASE):
    session = requests.Session()
    for name in CATALOGS:
        dest = catalog_path(name, root)
        if dest.exists():
            print(f"have    {dest}")
            continue
        size = fetch(session, catalog_url(name, base), dest)
        print(f"fetched {dest} ({size / 1e6:.1f} MB)")


def _rate(nbytes, seconds):
    return nbytes / 1e6 / seconds if seconds > 0 else 0.0


def download_lightcurves(filenames, root=None, workers=8, base=ARCHIVE_BASE, report_every=30):
    """Download the light curves in ``filenames`` that are not already present.

    Returns the list of ``(file name, error)`` that failed; these are also
    written to ``logs/failed_<timestamp>.txt`` under the data directory.
    """
    root = data_dir(root)
    filenames = list(dict.fromkeys(filenames))  # two workers must never write the same file
    have = downloaded_lc_names(root)
    todo = [name for name in filenames if name not in have]
    print(f"{len(filenames)} requested, {len(filenames) - len(todo)} present, {len(todo)} to download")
    if not todo:
        return []

    pending = queue.Queue()
    for name in todo:
        pending.put(name)
    lock = threading.Lock()
    stop = threading.Event()
    state = {"done": 0, "bytes": 0}
    failed = []

    def worker():
        session = requests.Session()
        while not stop.is_set():
            try:
                name = pending.get_nowait()
            except queue.Empty:
                return
            try:
                size = fetch(session, lc_url(name, base=base), local_lc_path(name, root))
                with lock:
                    state["done"] += 1
                    state["bytes"] += size
            except Exception as exc:  # keep going; failures are reported at the end
                with lock:
                    failed.append((name, repr(exc)))

    threads = [threading.Thread(target=worker, daemon=True) for _ in range(workers)]
    start = time.monotonic()
    for thread in threads:
        thread.start()

    def report():
        elapsed = time.monotonic() - start
        with lock:
            done, nbytes, nfailed = state["done"], state["bytes"], len(failed)
        left = len(todo) - done - nfailed
        eta = elapsed / done * left / 3600 if done else float("nan")
        print(
            f"{time.strftime('%Y-%m-%d %H:%M:%S')}  {done}/{len(todo)} files  "
            f"{nbytes / 1e9:.2f} GB  {_rate(nbytes, elapsed):.1f} MB/s  "
            f"{nfailed} failed  ETA {eta:.1f} h",
            flush=True,
        )

    try:
        while any(thread.is_alive() for thread in threads):
            deadline = time.monotonic() + report_every
            for thread in threads:
                thread.join(max(0.0, deadline - time.monotonic()))
            report()
    except KeyboardInterrupt:
        print("interrupted; letting in-flight files finish (Ctrl-C again to abort)", flush=True)
        stop.set()
        for thread in threads:
            thread.join()
        report()

    if failed:
        log = root / "logs" / f"failed_{time.strftime('%Y%m%dT%H%M%S')}.txt"
        log.parent.mkdir(parents=True, exist_ok=True)
        log.write_text("".join(f"{name}\t{error}\n" for name, error in failed))
        print(f"{len(failed)} failures written to {log}; re-run to retry them")
    return failed


def status(root=None):
    """Print how much of each catalog's light curve list is on disk."""
    root = data_dir(root)
    have = downloaded_lc_names(root)
    print(f"data directory: {root}")
    print(f"light curve files on disk: {len(have)}")
    for name in CATALOGS:
        if not catalog_path(name, root).exists():
            print(f"  {name:9s} catalog not downloaded")
            continue
        names = set(lc_filenames(name, root))
        present = len(names & have)
        print(f"  {name:9s} {present:7d} / {len(names):7d}  ({100 * present / len(names):.1f}%)")


def main(argv=None):
    parser = argparse.ArgumentParser(prog="kbonus-data", description=__doc__.split("\n")[0])
    parser.add_argument("--root", help="data directory (default: $KBONUS_DATA_DIR)")
    sub = parser.add_subparsers(dest="command", required=True)

    download = sub.add_parser("download", help="download catalogs or light curves")
    dsub = download.add_subparsers(dest="what", required=True)
    dsub.add_parser("catalogs", help="download the four catalogs")
    lcs = dsub.add_parser("lcs", help="download light curves (resumable)")
    lcs.add_argument(
        "--subset",
        choices=sorted(CATALOGS),
        default="source",
        help="catalog whose light curves to fetch; 'source' is all 606,900 (default)",
    )
    lcs.add_argument("--targets", nargs="+", help="specific targets (KIC numbers or Gaia DR3 ids)")
    lcs.add_argument("--limit", type=int, help="stop after this many missing files")
    lcs.add_argument("--workers", type=int, default=8, help="parallel connections (default 8)")

    sub.add_parser("status", help="show download progress per catalog")

    args = parser.parse_args(argv)
    if args.command == "status":
        status(args.root)
    elif args.what == "catalogs":
        download_catalogs(args.root)
    else:
        if args.targets:
            names = [lc_filename(target) for target in args.targets]
        else:
            if not catalog_path(args.subset, args.root).exists():
                download_catalogs(args.root)
            names = lc_filenames(args.subset, args.root)
        if args.limit is not None:
            have = downloaded_lc_names(args.root)
            names = [n for n in dict.fromkeys(names) if n not in have][: args.limit]
        failed = download_lightcurves(names, args.root, workers=args.workers)
        return 1 if failed else 0
    return 0


if __name__ == "__main__":
    sys.exit(main())
