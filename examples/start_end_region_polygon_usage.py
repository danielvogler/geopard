"""Draw a track together with the polygons used as start and finish regions.

Regions replace the circular radius when a segment's start is better described
by a shape than by a point -- a finish line across a road, say.

    uv run python examples/start_end_region_polygon_usage.py
"""

from __future__ import annotations

import argparse
import logging
import sys

from geopard import Geopard
from geopard.plotting import plot_regions
from geopard.settings import GPX_DATA_DIR, POLYGON_DATA_DIR

logging.basicConfig(level=logging.INFO, format="%(message)s")


def parse_args() -> argparse.Namespace:
    """Read the command line."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--activity", default=str(GPX_DATA_DIR / "tds_sunnestube_activity_25_25.gpx"))
    parser.add_argument("--start-region", default=str(POLYGON_DATA_DIR / "example_start_region.csv"))
    parser.add_argument("--finish-region", default=str(POLYGON_DATA_DIR / "example_finish_region.csv"))
    parser.add_argument("--save-to", help="Path to write a PNG to.")

    return parser.parse_args()


def main() -> int:
    """Draw the track and its regions."""
    args = parse_args()
    gp = Geopard()

    logging.info("Track:         %s", args.activity)
    logging.info("Start region:  %s", args.start_region)
    logging.info("Finish region: %s", args.finish_region)

    plot_regions(
        track=gp.gpx_loading(args.activity),
        start_region=gp.create_polygon(args.start_region),
        finish_region=gp.create_polygon(args.finish_region),
        show=args.save_to is None,
        save_to=args.save_to,
    )

    return 0


if __name__ == "__main__":
    sys.exit(main())
