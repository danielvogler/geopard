"""Figures.

Drawing is checked for the things that go wrong silently -- a swapped axis, a
squashed aspect ratio, a file that never gets written -- rather than for how it
looks. Everything runs on the Agg backend, so no window is ever opened.
"""

from __future__ import annotations

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import pytest

from geopard import plotting


@pytest.fixture(autouse=True)
def _close_figures():
    yield
    plt.close("all")


@pytest.fixture(scope="session")
def track(gp, gpx_dir):
    return gp.gpx_loading(str(gpx_dir / "tds_sunnestube_segment.gpx"))


def test_style_is_a_plain_mapping_of_rcparams():
    """Returned rather than applied, so importing geopard changes nobody's plots."""
    style = plotting.style()

    assert style["figure.facecolor"] == plotting.PAPER
    assert style["axes.spines.top"] is False


def test_importing_geopard_does_not_touch_global_rcparams():
    before = plt.rcParams["axes.facecolor"]
    plotting.style()

    assert plt.rcParams["axes.facecolor"] == before


def test_gpx_plot_draws_longitude_on_x_and_latitude_on_y(track):
    """The track rows are latitude-first and the axes are not. Easy to invert."""
    fig = plt.figure()
    plotting.gpx_plot(fig, track, ["Gold", ".", plotting.GOLD])

    axes = fig.gca()
    x_min, x_max = axes.get_xlim()
    y_min, y_max = axes.get_ylim()

    assert 9.0 < x_min < x_max < 10.0
    assert 47.0 < y_min < y_max < 48.0


def test_gpx_plot_labels_the_axes(track):
    fig = plt.figure()
    plotting.gpx_plot(fig, track, ["Gold", ".", plotting.GOLD])

    assert "Longitude" in fig.gca().get_xlabel()
    assert "Latitude" in fig.gca().get_ylabel()


def test_gpx_plot_corrects_the_aspect_ratio_for_latitude(track):
    """At 47°N a degree of longitude is ~0.68 of a degree of latitude.

    Without this the route is drawn a third too wide, and a good match looks
    like a bad one.
    """
    from math import cos, radians

    import numpy as np

    fig = plt.figure()
    plotting.gpx_plot(fig, track, ["Gold", ".", plotting.GOLD])

    mean_latitude = float(np.mean(np.asarray(track[0, :], dtype=float)))
    expected = 1.0 / cos(radians(mean_latitude))

    assert fig.gca().get_aspect() == pytest.approx(expected, rel=1e-9)
    assert 1.4 < expected < 1.5


def test_gpx_plot_returns_the_figure_so_calls_chain(track):
    fig = plt.figure()

    assert plotting.gpx_plot(fig, track, ["Gold", ".", plotting.GOLD]) is fig


def test_gpx_plot_adds_one_legend_entry_per_call(gp, track, gpx_dir):
    activity = gp.gpx_loading(str(gpx_dir / "tds_sunnestube_activity_25_25.gpx"))
    fig = plt.figure()

    plotting.gpx_plot(fig, track, ["Gold", ".", plotting.GOLD])
    plotting.gpx_plot(fig, activity, ["Activity", ".", plotting.INK])

    labels = [text.get_text() for text in fig.gca().get_legend().get_texts()]
    assert labels == ["Gold", "Activity"]


def test_track_comparison_returns_two_figures(gp, gpx_dir):
    figures = gp.plot_track_comparison(
        str(gpx_dir / "tds_sunnestube_segment.gpx"),
        str(gpx_dir / "tds_sunnestube_activity_25_25.gpx"),
        radius=7,
        show=False,
    )

    assert len(figures) == 2


def test_track_comparison_writes_both_files(gp, gpx_dir, tmp_path):
    prefix = tmp_path / "example"
    gp.plot_track_comparison(
        str(gpx_dir / "tds_sunnestube_segment.gpx"),
        str(gpx_dir / "tds_sunnestube_activity_25_25.gpx"),
        radius=7,
        show=False,
        save_to=str(prefix),
    )

    assert (tmp_path / "example_track.png").stat().st_size > 0
    assert (tmp_path / "example_track_interpolated.png").stat().st_size > 0


def test_plot_regions_draws_both_polygons(gp, gpx_dir, polygon_dir, tmp_path):
    figure = plotting.plot_regions(
        track=gp.gpx_loading(str(gpx_dir / "tds_sunnestube_activity_25_25.gpx")),
        start_region=gp.create_polygon(str(polygon_dir / "example_start_region.csv")),
        finish_region=gp.create_polygon(str(polygon_dir / "example_finish_region.csv")),
        show=False,
        save_to=str(tmp_path / "regions.png"),
    )

    labels = [text.get_text() for text in figure.gca().get_legend().get_texts()]

    assert "Start region" in labels
    assert "Finish region" in labels
    assert (tmp_path / "regions.png").stat().st_size > 0


def test_plot_regions_works_with_no_regions_at_all(gp, gpx_dir):
    figure = plotting.plot_regions(
        track=gp.gpx_loading(str(gpx_dir / "tds_sunnestube_segment.gpx")), show=False
    )

    assert figure is not None


def test_the_class_delegates_plotting_to_the_module(gp, track):
    """``gp.gpx_plot`` is part of the documented surface, not just the module."""
    fig = plt.figure()

    assert gp.gpx_plot(fig, track, ["Gold", ".", plotting.GOLD], 20) is fig
    assert fig.gca().get_legend() is not None
