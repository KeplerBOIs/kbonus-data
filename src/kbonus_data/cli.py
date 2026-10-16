"""Command line interface: ``kbonus-data download|status|index|channel``."""

import argparse
import sys

from .catalog import CATALOGS, catalog_path, lc_filenames
from .download import download_catalogs, download_lightcurves, downloaded_lc_names, status
from .paths import lc_filename


def _download(args):
    if args.what == "catalogs":
        download_catalogs(args.root)
        return 0
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


def _index(args):
    from .index import Index, build_index

    if args.action == "build":
        build_index(
            args.root, args.db, workers=args.workers, prefixes=args.prefixes, limit=args.limit, rebuild=args.rebuild
        )
        return 0
    with Index(args.root, args.db) as idx:
        for key, value in idx.summary().items():
            print(f"{key:13s} {value}")
        if args.quarter is not None:
            counts = idx.channel_counts(args.quarter)
            print(f"sources per channel in Q{args.quarter}:")
            for row in counts:
                print(f"  ch {row['channel']:2d}  {row['n']:6d}")
    return 0


def _channel(args):
    from .channels import DEFAULT_COLUMNS, read_channel
    from .index import Index

    with Index(args.root, args.db) as idx:
        cq = read_channel(
            args.quarter,
            args.channel,
            columns=args.columns or DEFAULT_COLUMNS,
            workers=args.workers,
            index=idx,
            progress=True,
            kepmag=(args.kepmag_min, args.kepmag_max),
            crowdsap=(args.crowdsap_min, None),
            limit=args.limit,
        )
    out = args.out or f"kbonus_q{args.quarter:02d}_ch{args.channel:02d}.npz"
    cq.save(out)
    print(f"{len(cq)} stars x {cq.shape[1]} cadences written to {out}")
    return 0


def main(argv=None):
    parser = argparse.ArgumentParser(prog="kbonus-data", description="KBonus-BKG Kepler light curves")
    parser.add_argument("--root", help="data directory (default: $KBONUS_DATA_DIR, else the current directory)")
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

    index = sub.add_parser("index", help="build or inspect the quarter/channel index")
    isub = index.add_subparsers(dest="action", required=True)
    build = isub.add_parser("build", help="scan light curve headers into the index (resumable)")
    build.add_argument("--workers", type=int, default=12, help="parallel readers (default 12)")
    build.add_argument("--prefixes", nargs="+", help="only these lcs/<prefix> directories")
    build.add_argument("--limit", type=int, help="stop after this many new files")
    build.add_argument("--rebuild", action="store_true", help="rescan files already indexed")
    info = isub.add_parser("status", help="show index counts")
    info.add_argument("--quarter", type=int, help="also list sources per channel in this quarter")
    channel = sub.add_parser("channel", help="save all light curves of one quarter and channel to .npz")
    for p in (build, info, channel):
        p.add_argument("--db", help="index file (default: <root>/index/kbonus_index.sqlite)")

    channel.add_argument("quarter", type=int)
    channel.add_argument("channel", type=int)
    channel.add_argument("--out", help="output file (default kbonus_q<QQ>_ch<CC>.npz)")
    channel.add_argument("--kepmag-min", type=float)
    channel.add_argument("--kepmag-max", type=float)
    channel.add_argument("--crowdsap-min", type=float, help="minimum CROWDSAP (1 = uncrowded)")
    channel.add_argument("--limit", type=int, help="at most this many stars, brightest first")
    channel.add_argument("--columns", nargs="+", help="per-quarter columns to read")
    channel.add_argument("--workers", type=int, default=16)

    args = parser.parse_args(argv)
    if args.command == "status":
        status(args.root)
        return 0
    return {"download": _download, "index": _index, "channel": _channel}[args.command](args)


if __name__ == "__main__":
    sys.exit(main())
