"""Match an activity against a gold-standard segment, and plot the result.

Run with the shipped example tracks::

    uv run python examples/track_matching.py

Or against your own::

    uv run python examples/track_matching.py \
        --gold my_segment.gpx --activity my_activity.gpx --radius 10
"""

from __future__ import annotations

import argparse
import logging
import sys

from geopard import Geopard
from geopard.plotting import plot_regions
from geopard.settings import GPX_DATA_DIR

logging.basicConfig(level=logging.INFO, format="%(message)s")


def parse_args() -> argparse.Namespace:
    """Read the command line."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--gold",
        default=str(GPX_DATA_DIR / "tds_sunnestube_segment.gpx"),
        help="GPX track other activities are compared against.",
    )
    parser.add_argument(
        "--activity",
        default=str(GPX_DATA_DIR / "tds_sunnestube_activity_25_25.gpx"),
        help="GPX track to evaluate.",
    )
    parser.add_argument(
        "--radius",
        type=float,
        default=7.0,
        help="Metres around the gold start/finish. Default: %(default)s",
    )
    parser.add_argument("--start-region", help="CSV polygon to use instead of a radius.")
    parser.add_argument("--finish-region", help="CSV polygon to use instead of a radius.")
    parser.add_argument("--no-plot", action="store_true", help="Skip the figures.")

    return parser.parse_args()


def main() -> int:
    """Match, report, and plot."""
    args = parse_args()
    gp = Geopard()

    start_region = gp.create_polygon(args.start_region) if args.start_region else None
    finish_region = gp.create_polygon(args.finish_region) if args.finish_region else None

    logging.info("Matching %s against %s", args.activity, args.gold)
    response = gp.dtw_match(
        gold_name=args.gold,
        activity_name=args.activity,
        radius=args.radius,
        start_region=start_region,
        finish_region=finish_region,
    )

    if not response.is_success():
        logging.error("No match. Flag %s: %s", response.match_flag, response.error)
        return 1

    gp.parse_response(response)

    if args.no_plot:
        return 0

    gp.plot_track_comparison(
        gold_file_name=args.gold,
        activity_file_name=args.activity,
        radius=args.radius,
    )

    if start_region or finish_region:
        plot_regions(gp.gpx_loading(args.activity), start_region, finish_region)

    return 0


if __name__ == "__main__":
    sys.exit(main())
