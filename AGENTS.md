# AGENTS.md

The single source of truth for anyone working in this repository, human or
agent. `CLAUDE.md` and `GEMINI.md` are pointers here and hold no content of
their own.

**There are two jobs in this file.** Matching tracks with the library is §A.
Changing the code is §B. Read the one you are here for.

---

# §A — Using the library

## A1. What you are holding

A Python library that answers one question: **did this activity cover that
segment, and how fast?** Two GPX files in, one `GeopardResponse` out.

It is a library, not a service. It reads local files, holds no credentials,
makes no network calls, and writes nothing except figures you explicitly ask it
to save.

## A2. Setup

```bash
make setup       # venv, dependencies, git hooks
```

Or, as a dependency of something else: `pip install geopard` /
`uv add geopard`.

Check it works without writing anything:

```bash
make test-fast   # the suite minus the end-to-end matches
make match       # matches a shipped example and opens two figures
```

## A3. The one call that matters

```python
import geopard

gp = geopard.Geopard()
response = gp.dtw_match("segment.gpx", "activity.gpx")
```

**Always branch on `response.is_success()` before reading anything else.** Every
field except `match_flag` is `None` when no match was found, and
`response.time.seconds` on a failed match is an `AttributeError` several frames
away from the cause.

```python
if not response.is_success():
    print(response.match_flag, response.error)
    return

print(response.time, response.dtw)
```

**`dtw_match` does not raise when there is no match.** It returns a response
carrying a negative flag. It *does* raise `GeopardException` when the activity
file is missing latitude, longitude or time — a corrupt file is a different
thing from a track that went somewhere else, and only one of them is worth
stopping for.

## A4. Reading the flag

| `match_flag` | `MatchFlag` | What happened | What to do |
|---|---|---|---|
| `2` | `MATCH` | Distance at or below `dtw_threshold` | Use the result |
| `1` | `GREY_ZONE` | Between the threshold and `dtw_threshold × dtw_margin_range` | Look at it. This is a real attempt that nearly qualified — do not report it as a clean match, and do not silently drop it |
| `-1` | `NO_MATCH` | Candidate start/finish pairs existed; none matched | The activity went past both ends but did not follow the route |
| `-2` | `NOT_EVALUATED` | Never ran. `error` says why | Usually nothing near the start or the finish. Check `radius` before concluding the activity is unrelated |

`MatchFlag` is an `IntEnum`, so `response.match_flag == 2` and
`response.match_flag == MatchFlag.MATCH` are the same test. Both are supported
deliberately; callers have compared against the raw integers since the first
release.

## A5. When a match is not what you expected

Work down this list rather than reaching for `dtw_threshold` first. Loosening
the threshold is the one change that makes *every* result less trustworthy.

**Flag `-2`, error "No trackpoints found near centroid".** Nothing came within
`radius` of the segment's start or its finish. Either the activity really did
not go there, or `radius` (default 7 m) is too tight for the reception. Raise it
to 15-20 m before concluding anything. Under tree cover or between buildings,
7 m is genuinely tight.

**Flag `-1`, no error.** Both ends were reached but no pairing warped closely
enough. Plot it before touching any threshold:

```python
gp.plot_track_comparison("segment.gpx", "activity.gpx", radius=7)
```

Usually the answer is visible — a diversion, a shortcut, a wrong turn. If the
two curves genuinely do lie on top of each other, then raise `dtw_threshold`.

**A match, but the time looks too short.** A start candidate was paired with a
finish candidate a few seconds later, on a course whose start and finish regions
overlap. Raise `min_trkps`.

**A match, but the time looks too long.** The fastest pairing failed the
threshold and a slower one passed. Plot it; often the fast attempt was the one
with the GPS glitch.

**Repeated laps.** Nothing special is needed. Every start is paired with every
later finish and the shortest qualifying pair wins, so the result is the best lap
in the activity.

## A6. Reproducibility

Interpolation jitters each coordinate by up to 1e-7 degrees — about a centimetre
— because the spline solver cannot fit through duplicate points. It is the only
stochastic step in the library.

It moves `response.dtw` in the fourth decimal, roughly 0.1% run to run. It does
**not** move `time`, `match_flag`, `start_point` or `end_point`, which are
bit-stable.

**Never report a warping distance to more figures than that.** If you need a
bit-identical number, seed it:

```python
np.random.seed(1234)  # pins a whole pipeline
interpolate(track, seed=1234)  # or one call, leaving the global RNG alone
```

## A7. Costs and limits

- **The accumulated cost matrix is the slow step**, quadratic in the
  interpolation count (1000 × 1000) and written as a Python loop because the
  recurrence is sequential. Budget roughly a second per candidate pair.
- Candidate pairs are tested **shortest first and the loop stops at the first
  match**, so a clean activity costs one warp, and a difficult one costs many.
  Raising `radius` raises the pair count sharply.
- **Elevation is loaded but never used.** Matching is two-dimensional. A
  multi-storey car park matches itself on every level.
- Distances are great-circle on a sphere of mean Earth radius: accurate to about
  0.5%, far below GPS error at these scales.

## A8. What is deliberately absent

**No CLI.** The library is a handful of functions; the examples in `examples/`
are the interface for running it by hand, and they are scripts you can read and
change rather than a command surface to keep stable.

**No batch or database layer.** Matching a season of activities is a `for` loop
over `dtw_match`, and what you want to do with the results is not something this
repository can guess.

**No elevation matching, and no "closest segment out of N" search.** Both are
real features. Neither is a small addition to what is here, and inventing them
before there is a use is how the shape of a library goes wrong.

---

# §B — Changing the code

## B1. Layout

```
src/geopard/
  matching.py        The Geopard class and dtw_match — the pipeline, in order
  selection.py       nearest_neighbours, crop_track — which stretches to test
  dtw.py             Accumulated cost matrix, and one warp
  interpolation.py   Spline resampling to a fixed point count
  gpx.py             GPX file -> (4, n) array; duplicate and corruption handling
  geodesy.py         Haversine distance
  regions.py         CSV -> Shapely polygon
  response.py        GeopardResponse, MatchFlag, log_response
  plotting.py        Palette, style(), and the figures
  exceptions.py      GeopardException
  settings.py        PROJECT_ROOT, for the shipped examples only
  geopard.py         Backwards-compatible import location — do not add to it
tests/               One file per module, plus the end-to-end benchmark
examples/            Runnable scripts, not a stable interface
data/                Two real segments with matching activities
```

Each row of a track is `[latitude, longitude, elevation, time]`, in that order,
and `geopard.gpx` names the offsets. Latitude comes first everywhere in the
library and **longitude comes first in every Shapely geometry**, because Shapely
is `(x, y)`. `regions.create_polygon` and `selection.within_region` are the two
places that swap, and they are the two places to look when a region silently
contains nothing.

`Geopard` holds no state. Its methods delegate to the modules above, and exist so
the whole library is reachable from one object — which is how it has always been
used. Add a function to its module first, then delegate.

## B2. The tests are a benchmark. Treat them as one.

Every elapsed time, flag and trackpoint in `tests/test_matching.py` was measured
against the shipped GPX tracks and is asserted **exactly**. Only `dtw` carries a
tolerance — `DTW_RTOL` in `conftest.py`, 1%, sized from the observed run-to-run
spread of the jitter described in §A6.

**If a change moves one of those numbers, the change is wrong until you can say
why it is right.** That is the whole value of them. They are what made it safe to
split a 750-line module into ten, and to move Shapely across a major version,
without wondering whether the results had drifted.

When adding behaviour:

- Cover the failure path, not just the happy one. The flags that matter most in
  practice are `-1` and `-2`, and they are returned rather than raised, so a test
  that only checks a successful match tells you almost nothing.
- Prefer a property that must hold — distance is symmetric, the crop is
  contiguous in time, interpolation preserves the endpoints — over a value copied
  from a run. Recorded values pin behaviour; properties explain it.
- Mark anything doing a full `dtw_match` over the shipped tracks with
  `@pytest.mark.slow`, so `make test-fast` stays quick.
- No test may need a network or write outside `tmp_path`.

## B3. Checks

```bash
make check      # ruff check, ruff format --check, mypy, pytest — what CI runs
```

Everything must pass before a change is handed over. CI additionally builds the
wheel and imports it, because a packaging mistake is invisible to a test suite
that runs against the source tree.

## B4. Conventions

- **`uv` for everything.** `ruff` formats and lints (line length 110, Google
  docstrings), `mypy` type-checks, `pytest` tests. Commit `uv.lock`.
- **Type hints on every public signature.** Use the `Coordinates` alias from
  `geopard.geodesy` for a point — callers pass Python lists *and* numpy columns,
  and `Sequence[float]` misdescribes half of them.
- **`logging`, never `print`**, outside `examples/`. The library logs a lot at
  INFO, per trackpoint in places, which is why the test suite disables it.
- **`pathlib.Path`, never string concatenation** for paths.
- **No bare `except`, and no `except Exception: pass`.** `dtw_match` catches
  exactly `GeopardException` and turns it into a flag; that is the only place
  swallowing anything, and it is documented as the contract.
- **Seed every stochastic process explicitly**, and say so in the docstring.
- **Module-level helpers carry no leading underscore.** It signals "internal"
  and enforces nothing. Methods on a class may.
- Commits are `<type>(<scope>): <description>` — `feat`, `fix`, `refactor`,
  `docs`, `test`, `chore`, `ci`, `perf`. The scope is the area touched:
  `fix(selection): ...`. This repository does not use gitmoji; do not introduce
  it. A commit-msg hook rejects AI `Co-Authored-By` trailers.
- **Never commit an absolute path** — not in a message, not in code, not in a
  committed agent file. A pre-commit hook refuses `/Users/…` and `/home/…`.

## B5. The things that will bite you

**The version lives in one place**, `src/geopard/version.py`, and `pyproject.toml`
reads it through hatchling's dynamic version. `tests/test_public_api.py` fails if
that wiring is broken. Do not add a second literal.

**`geopard/geopard.py` is a compatibility shim.** `from geopard.geopard import
Geopard` was the only import path before the package had an `__init__`, and
people's scripts use it. It re-exports and nothing else. New code goes in the
module it belongs to.

**`settings.PROJECT_ROOT` is for the examples, not for the library or the
tests.** It resolves the repository root from the source file's location, which
is meaningless in an installed wheel. The tests resolve their own paths in
`conftest.py` on purpose — a change to `settings` must not be able to silently
repoint the suite at different files.

**Gold segments may carry no timestamps at all.** `tds_sunnestube_segment.gpx` is
a drawn route, not a recording. Only the activity's clock is ever read, and
`validate_track` is called on the activity alone. Validating the gold would break
the shipped example.

**`plotting.style()` returns rcParams rather than applying them.** Importing
geopard must never change how somebody else's plots look. Use it with
`plt.rc_context`.

**Latitude/longitude order.** See §B1. It is the bug this codebase is most likely
to grow.
