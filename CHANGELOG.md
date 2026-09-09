# Changelog

All notable changes to this project are documented here. The format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and this project
adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

## [0.3.0] - 2026-09-09

A restructuring release. **No matching behaviour changed**: every elapsed time,
match flag and matched trackpoint the library produced before it produces after,
verified against a benchmark suite written before any code was touched.

### Added

- `geopard/__init__.py`, so the `import geopard; geopard.Geopard()` the README
  has always documented actually works. It did not before — the package had no
  `__init__`, and only `from geopard.geopard import Geopard` resolved.
- `MatchFlag`, an `IntEnum` naming the four outcomes (`MATCH`, `GREY_ZONE`,
  `NO_MATCH`, `NOT_EVALUATED`). The integers are unchanged, so existing
  comparisons against `2` and `-1` keep working.
- `interpolate(..., seed=...)`, pinning the coordinate jitter without touching
  the caller's own global RNG.
- `plotting.style()`, `plotting.plot_regions()`, and `save_to` / `show`
  arguments on `plot_track_comparison`, so figures can be produced in a script
  or in CI without opening a window.
- A benchmark test suite: 133 tests, 98% coverage, up from 5 tests. Times, flags
  and trackpoints are asserted exactly; only the warping distance carries a
  tolerance, sized from the measured run-to-run spread.
- `AGENTS.md`, `CLAUDE.md`, `GEMINI.md`, `SECURITY.md`, `CHANGELOG.md`,
  `docs/architecture.md`, and a branded README.

### Changed

- **The 750-line `geopard/geopard.py` is now ten focused modules** under
  `src/geopard/`. `Geopard` keeps every method it had, delegating to them.
  `from geopard.geopard import Geopard` still resolves, and is tested.
- **Packaging moved to `uv` + `hatchling`**, replacing Poetry, `setup.py`,
  `setup.cfg` and `requirements.txt`. The version now has one home,
  `src/geopard/version.py`, read dynamically by `pyproject.toml`.
- **Linting moved to `ruff`**, replacing black, isort, flake8, pylint,
  pydocstyle and autoflake. `mypy` was added and passes.
- **Python floor raised to 3.11**, and **Shapely unpinned to 2.x**. The old
  `Shapely <2` pin had no wheels for any supported Python, so a fresh install
  fell back to compiling from source and usually failed. Verified numerically
  identical across the upgrade.
- `GeopardResponse` is a frozen dataclass. It previously declared six methods
  that were shadowed by instance attributes of the same name and could never be
  called. `is_success()` remains a method, because the shipped examples call it.
- Figures were restyled: a single palette, aspect-ratio correction for latitude
  (a raw scatter squashes a route by a third at 47°N), legends outside the data,
  and readable marker sizes in place of the previous 500-point markers.
- Example scripts take real `--flags` with defaults, and exit with a status code.

### Fixed

- `logging.info("Centroid: %f", centroid[:2])` in `nearest_neighbours` formatted
  a list with `%f`, which failed on every call. It was invisible because the
  logging module swallows formatting errors.
- `parse_response` raised a formatting error on an unmatched response, where
  `dtw` is `None`.
- `requirements.txt` listed `haversine` and `similaritymeasures`, neither of
  which the library imports. Removed with the file.
- The PyPI workflow ran on every push *and every pull request* to `main`, rather
  than on a release. It now triggers on a published release.
- `nearest_neighbours` had `centroid=list` as a default — the type object itself.

### Removed

- `README.rst`, a stale duplicate of `README.md` documenting an API that had
  moved on.
- `setup.py`, `setup.cfg`, `requirements.txt`, `poetry.lock`, `.pylintrc`, and
  the `github-release` workflow, all superseded.

## [0.2.3] - 2023-10-31

Earlier releases predate this changelog. See the
[commit history](https://github.com/danielvogler/geopard/commits/main).
