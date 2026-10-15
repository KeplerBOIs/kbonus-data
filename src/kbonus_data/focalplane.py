"""Kepler focal plane numbering: modules, outputs, channels, seasons.

The 21 science modules are numbered 2-4, 6-20 and 22-24 on a 5x5 grid (the
four corners hold the fine guidance sensors).  Each module has four outputs,
giving channels 1-84 in module order.  The spacecraft rolled by 90 degrees
between quarters, so a star falls on one of (up to) four channels depending
on the season.
"""

MODULES = (2, 3, 4, *range(6, 21), 22, 23, 24)
CHANNELS = tuple(range(1, 85))
QUARTERS = tuple(range(18))


def channel(module, output):
    """Return the channel number (1-84) of a module and output (1-4)."""
    if module not in MODULES:
        raise ValueError(f"module must be one of {MODULES}, not {module!r}")
    if output not in (1, 2, 3, 4):
        raise ValueError(f"output must be 1-4, not {output!r}")
    return MODULES.index(module) * 4 + output


def module_output(channel):
    """Return ``(module, output)`` for a channel number (1-84)."""
    if channel not in CHANNELS:
        raise ValueError(f"channel must be 1-84, not {channel!r}")
    index, output = divmod(channel - 1, 4)
    return MODULES[index], output + 1


def season(quarter):
    """Return the season (0-3) of a quarter, as in the SEASON header keyword.

    Q0 and Q1 share season 3 with Q5, Q9, Q13 and Q17.
    """
    if quarter not in QUARTERS:
        raise ValueError(f"quarter must be 0-17, not {quarter!r}")
    return 3 if quarter == 0 else (quarter + 2) % 4
