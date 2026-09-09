"""Show what cropping and interpolation do to a track, before any matching.

uv run python examples/gpx_track_cropping_and_interpolation.py
"""

from __future__ import annotations

import argparse
import logging
import sys

from geopard import Geopard
from geopard.settings import GPX_DATA_DIR

logging.basicConfig(level=logging.INFO, format="%(message)s")


def parse_args() -> argparse.Namespace:
    """Read the command line."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--gold", default=str(GPX_DATA_DIR / "tds_sunnestube_segment.gpx"))
    parser.add_argument("--activity", default=str(GPX_DATA_DIR / "tds_sunnestube_activity_25_25.gpx"))
    parser.add_argument("--radius", type=float, default=20.0)
    parser.add_argument("--save-to", help="Path prefix for PNG output, e.g. docs/images/example")

    return parser.parse_args()


def main() -> int:
    """Draw the crop and the interpolated curves."""
    args = parse_args()
    gp = Geopard()

    gold = gp.gpx_loading(args.gold)
    activity = gp.gpx_loading(args.activity)
    cropped = gp.gpx_track_crop(gold=gold, gpx_data=activity, radius=args.radius)

    logging.info("Gold segment:      %s trackpoints", gold.shape[1])
    logging.info("Activity:          %s trackpoints", activity.shape[1])
    logging.info("Activity, cropped: %s trackpoints", cropped.shape[1])

    gp.plot_track_comparison(
        gold_file_name=args.gold,
        activity_file_name=args.activity,
        radius=args.radius,
        show=args.save_to is None,
        save_to=args.save_to,
    )

    return 0


if __name__ == "__main__":
    sys.exit(main())
