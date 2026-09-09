"""Drawing tracks.

A GPX track plotted as a raw latitude/longitude scatter is misleading: away
from the equator a degree of longitude is much shorter than a degree of
latitude, so a matplotlib default squashes the route sideways and a match that
is geometrically fine looks wrong. Everything here corrects for that, and
otherwise stays out of the way.
"""

from __future__ import annotations

import logging
from math import cos, radians
from typing import Any, cast

import numpy as np
from matplotlib import pyplot as plt
from matplotlib.figure import Figure

logger = logging.getLogger(__name__)

# --- palette -------------------------------------------------------------
# Ink on paper, with the segment picked out in amber and the matched endpoints
# in red. Deliberately few colours: a track plot with six of them stops being
# readable at exactly the point it gets interesting.
PAPER = "#FDFDFB"
INK = "#0E0E10"
MUTED = "#6F6F6B"
FAINT = "#A3A39D"
RULE = "#E5E4DE"
GOLD = "#E8820C"
GOLD_FAINT = "#F0A85C"
ACCENT = "#FF000D"

#: Default marker area, in points squared.
MARKER_SIZE = 14
#: Marker area for the handful of points that carry meaning on their own.
HIGHLIGHT_MARKER_SIZE = 90

DEFAULT_FIGSIZE = (12.0, 8.0)
DEFAULT_DPI = 110


def style() -> dict[str, Any]:
    """The rcParams geopard's figures use.

    Returned rather than applied globally, so importing geopard never changes
    how anybody else's plots look. Use it as a context manager::

        with plt.rc_context(geopard.plotting.style()):
            ...
    """
    return {
        "figure.facecolor": PAPER,
        "axes.facecolor": PAPER,
        "savefig.facecolor": PAPER,
        "axes.edgecolor": RULE,
        "axes.labelcolor": MUTED,
        "axes.titlecolor": INK,
        "axes.grid": True,
        "axes.axisbelow": True,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "grid.color": RULE,
        "grid.linewidth": 0.8,
        "text.color": INK,
        "xtick.color": FAINT,
        "ytick.color": FAINT,
        "xtick.labelcolor": MUTED,
        "ytick.labelcolor": MUTED,
        "font.size": 11,
        "axes.titlesize": 15,
        "axes.titleweight": "bold",
        "axes.labelsize": 11,
        "legend.fontsize": 10,
        "legend.frameon": False,
        "figure.dpi": DEFAULT_DPI,
        "scatter.edgecolors": "none",
    }


def gpx_plot(
    fig: Figure,
    gpx_data: np.ndarray,
    plot_info: Any,
    marker_size: int | None = None,
) -> Figure:
    """Scatter one track onto a figure.

    Args:
        fig: The figure to draw on.
        gpx_data: A track whose first two rows are latitude and longitude.
        plot_info: ``[label, marker, colour]``, as matplotlib understands them.
        marker_size: Marker area in points squared.

    Returns:
        The same figure, so calls can be chained.
    """
    logger.debug("Plot gpx data")

    label, marker, colour = plot_info[0], plot_info[1], plot_info[2]
    axes = fig.gca()

    axes.scatter(
        gpx_data[1, :],
        gpx_data[0, :],
        s=MARKER_SIZE if marker_size is None else marker_size,
        marker=marker,
        c=colour,
        label=label,
        linewidths=0,
    )

    finish_axes(axes, gpx_data)

    return fig


def finish_axes(axes: Any, gpx_data: np.ndarray) -> None:
    """Label the axes, fix the aspect ratio, and place the legend."""
    axes.set_xlabel("Longitude [°]")
    axes.set_ylabel("Latitude [°]")
    axes.legend(loc="upper left", bbox_to_anchor=(1.01, 1.0), markerscale=1.8, borderaxespad=0)

    # One degree of longitude is cos(latitude) as long as one of latitude.
    # Scaling by the inverse makes a metre look the same on both axes, which
    # is the difference between a track and a smear.
    mean_latitude = float(np.mean(np.asarray(gpx_data[0, :], dtype=float)))
    axes.set_aspect(1.0 / max(cos(radians(mean_latitude)), 1e-6))


def new_figure(title: str | None = None) -> Figure:
    """A figure carrying geopard's style."""
    fig = plt.figure(figsize=DEFAULT_FIGSIZE, dpi=DEFAULT_DPI, facecolor=PAPER)

    if title:
        fig.suptitle(title, x=0.01, ha="left", color=INK, fontsize=16, weight="bold")

    return fig


def plot_track_comparison(
    gold_file_name: str,
    activity_file_name: str,
    radius: float = 7.0,
    show: bool = True,
    save_to: str | None = None,
) -> list[Figure]:
    """Draw a gold segment, an activity, and what the crop kept.

    Two figures: the tracks as recorded, and the interpolated curves that
    dynamic time warping actually compares. The second is the one worth looking
    at when a match is closer or further than expected.

    Args:
        gold_file_name: Path to the gold segment GPX file.
        activity_file_name: Path to the activity GPX file.
        radius: Radius in metres around the gold start and finish.
        show: Whether to call ``plt.show()``. Turn it off in a script that
            only wants the files.
        save_to: Path prefix to write PNGs to. ``"docs/images/example"`` writes
            ``example_track.png`` and ``example_track_interpolated.png``.

    Returns:
        The two figures, in that order.
    """
    from geopard import gpx as gpx_module
    from geopard import interpolation, selection

    gold = gpx_module.load_track(gold_file_name)
    gold_interpolated = interpolation.interpolate(gold)

    activity = gpx_module.load_track(activity_file_name)
    cropped = selection.crop_track(gold, activity, radius)

    near_start, _ = selection.nearest_neighbours(cropped, gold[:4, 0], radius)
    near_finish, _ = selection.nearest_neighbours(cropped, gold[:4, -1], radius)

    with plt.rc_context(cast(Any, style())):
        tracks = new_figure("Activity cropped to the gold segment")
        gpx_plot(tracks, gold, ["Gold segment", ".", GOLD], MARKER_SIZE * 4)
        gpx_plot(tracks, activity, ["Activity", ".", FAINT], MARKER_SIZE)
        gpx_plot(tracks, cropped, ["Activity, cropped", ".", INK], MARKER_SIZE)
        gpx_plot(tracks, near_start, ["Start candidates", "X", ACCENT], HIGHLIGHT_MARKER_SIZE)
        gpx_plot(tracks, near_finish, ["Finish candidates", "P", ACCENT], HIGHLIGHT_MARKER_SIZE)

        curves = new_figure("What dynamic time warping actually compares")
        gpx_plot(curves, gold, ["Gold, recorded", ".", GOLD_FAINT], MARKER_SIZE * 9)
        gpx_plot(curves, cropped, ["Activity, recorded", ".", FAINT], MARKER_SIZE * 5)
        gpx_plot(curves, gold_interpolated.T, ["Gold, interpolated", ".", GOLD], MARKER_SIZE)
        gpx_plot(
            curves,
            interpolation.interpolate(cropped).T,
            ["Activity, interpolated", ".", INK],
            MARKER_SIZE,
        )

        for fig in (tracks, curves):
            fig.tight_layout()

        if save_to:
            tracks.savefig(f"{save_to}_track.png", bbox_inches="tight", dpi=DEFAULT_DPI)
            curves.savefig(
                f"{save_to}_track_interpolated.png",
                bbox_inches="tight",
                dpi=DEFAULT_DPI,
            )

        if show:
            plt.show()

    return [tracks, curves]


def plot_regions(
    track: np.ndarray,
    start_region: Any | None = None,
    finish_region: Any | None = None,
    show: bool = True,
    save_to: str | None = None,
) -> Figure:
    """Draw a track with its start and finish region polygons.

    Args:
        track: A ``(4, n)`` track.
        start_region: Polygon around the segment start.
        finish_region: Polygon around the segment finish.
        show: Whether to call ``plt.show()``.
        save_to: Path to write a PNG to.

    Returns:
        The figure.
    """
    with plt.rc_context(cast(Any, style())):
        fig = new_figure("Start and finish regions")
        axes = fig.gca()

        for region, colour, label in (
            (start_region, ACCENT, "Start region"),
            (finish_region, INK, "Finish region"),
        ):
            if region is not None:
                axes.fill(*region.exterior.xy, color=colour, alpha=0.12)
                axes.plot(*region.exterior.xy, color=colour, linewidth=2, label=label)

        gpx_plot(fig, track, ["Track", ".", GOLD], MARKER_SIZE)
        fig.tight_layout()

        if save_to:
            fig.savefig(save_to, bbox_inches="tight", dpi=DEFAULT_DPI)

        if show:
            plt.show()

    return fig
