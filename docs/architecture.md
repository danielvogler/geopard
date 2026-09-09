# Architecture

Why geopard is shaped the way it is, and the decisions that look wrong until you
know why. For how to *use* it, see [AGENTS.md](../AGENTS.md) §A.

---

## The problem, precisely

You have a segment — a stretch of route, as GPX — and an activity that may have
covered it. Answering "did it, and how fast" sounds like a lookup and is not,
for three reasons.

**A start line is not a point.** It is every trackpoint within a few metres of
one, and so is the finish. A GPS recording at one point per second passing
through a 7 m circle contributes several candidates at each end.

**An activity may cover the segment more than once.** Laps, an out-and-back, a
warm-up loop through the start. Each start pairs with each later finish, so
*k* starts and *m* finishes give up to *k×m* plausible attempts inside one file.

**Passing through both ends is not the same as following the route.** A shortcut
that skips a switchback starts and finishes in exactly the right places. Only the
shape between them distinguishes it.

So the pipeline is: narrow down which stretches are even plausible, then compare
their shapes, in an order that lets the search stop early.

---

## The pipeline

```
load ──> validate ──> crop ──> candidate pairs ──> interpolate ──> warp ──> flag
```

| Stage | Module | Decision it makes |
|---|---|---|
| `load_track` | `gpx.py` | Which trackpoints exist at all |
| `validate_track` | `gpx.py` | Whether the activity is usable |
| `nearest_neighbours` | `selection.py` | Which trackpoints count as "at the start" |
| `crop_track` | `selection.py` | Which time window could hold the segment |
| `candidate_segments` | `matching.py` | Which start/finish pairs are worth testing, and in what order |
| `interpolate` | `interpolation.py` | What the two curves look like at a common resolution |
| `accumulated_cost_matrix` | `dtw.py` | How different the two shapes are |
| `GeopardResponse` | `response.py` | What the caller is told |

---

## Decisions worth explaining

### Shortest first, and stop at the first match

Candidates are sorted by elapsed time and warped in that order. The loop exits
the moment one clears `dtw_threshold`.

This is not an optimisation, it is the definition of the answer. Every remaining
pair is slower by construction, so the first qualifying pair *is* the best time
run on the segment. Warping all of them and taking the minimum would give the
same result at many times the cost; warping them in arrival order would give a
different, wrong one.

It is also what makes the cost bearable. A clean activity costs one warp. Only a
difficult one — a wobbly GPS trace, or an activity that passed the start several
times without running the segment — pays for many.

### Dynamic time warping, rather than a point-wise distance

Two people run the same route at different speeds and their trackpoints do not
correspond: at trackpoint 300 one is at the top of the climb and the other is
halfway up. Comparing point *i* to point *i* measures pace, not route.

Warping aligns the two sequences by allowing either to stretch, so the comparison
is of shape alone. The recurrence is Müller (2007) Theorem 4.3:

```
D(n,m) = c(n,m) + min{ D(n-1,m-1), D(n-1,m), D(n,m-1) }
```

with the first row and column as running sums, since with one point on the other
side there is no alignment to choose.

`accumulated_cost_matrix` is a Python double loop over a 1000 × 1000 matrix, and
it is by far the slowest thing here. It stays a loop because the recurrence is
sequential — `D(n,m)` reads `D(n,m-1)`, computed in the same pass — so the inner
dimension does not vectorise. The row-wise dependency, not the array library, is
the constraint.

### Interpolating to a fixed 1000 points

Warping needs both sequences at a comparable resolution. A drawn 61-point route
and a recorded 1616-point activity are not, and the count is not a property of
the route — it is a property of whatever device recorded it.

Resampling both to the same fixed count removes the device from the answer. 1000
is a module constant rather than an argument because *both* curves must agree on
it; exposing it per-call would invite two curves resampled differently, which
compares nothing meaningful.

### The jitter, and why it is not a bug

`splprep` cannot fit a spline through repeated points. Stationary trackpoints are
already dropped at load, but a track can still contain two consecutive fixes
equal to the last decimal, so each coordinate is perturbed by up to 1e-7 degrees
— about a centimetre, two orders of magnitude below GPS error.

It costs reproducibility in the fourth decimal of `dtw`, and nothing else: flags,
elapsed times and matched trackpoints are bit-stable. `interpolate(seed=...)`
pins it when a bit-identical number is needed, and the default path stays on
numpy's global RNG so `np.random.seed(...)` keeps working for callers who already
rely on it.

### Failures are returned, not raised

Matching an activity that never went near the segment is an ordinary outcome, not
an error. A caller sweeping a season of files should not need a `try` around
every call, so `dtw_match` catches `GeopardException` from the selection stage
and returns flag `-2` with the message.

A *corrupt* activity — missing latitude, longitude or time — still raises. That
is a different kind of problem: it means the input is wrong, not that the answer
is no.

### The grey zone

A single threshold turns a continuous measurement into a yes/no, and the
interesting cases sit exactly at the line: a GPS drifting under tree cover, a
diversion around roadworks.

Flag `1` reports the shortest attempt that fell between `dtw_threshold` and
`dtw_threshold × dtw_margin_range`. Crucially it does **not** stop the search — a
later, slower candidate may still be a proper match, and a proper match always
wins over a provisional one. It is only reported if nothing better turns up.

### Elevation is loaded and never used

Every track carries an elevation row, and the distance ignores it. Two-dimensional
matching is a genuine limitation — the same ground plan on a different level
matches — but adding a third dimension means choosing how to weight metres of
climb against metres of ground, and that weight is a modelling decision the
library has no basis to make on the caller's behalf.

### `Geopard` is a facade over free functions

The class holds no state. Every method delegates to a module that can be imported
on its own, which is what makes the modules testable in isolation and the
pipeline legible as a sequence of steps.

The class stays because it is the documented interface and always has been.
Removing it to be tidy would break every existing caller for no gain.

---

## The two coordinate orders

This is the bug the codebase is most likely to grow.

- **Tracks are latitude-first**: every array row is `[lat, lon, ele, time]`, and
  `geopard.gpx` names the offsets.
- **Shapely is `(x, y)`, which means longitude-first.**
- **Region CSVs are latitude-first**, because that is how a person writes a
  coordinate down.

Exactly two functions cross between them — `regions.create_polygon` on the way in,
and `selection.within_region` when testing containment. A region that silently
contains no trackpoints is what a mistake here looks like; there is no error, just
an empty result and a flag `-2`. Both conversions are covered by tests that assert
on WKT and on containment against a real track, which is the only way this fails
loudly.

---

## Layout

```
src/geopard/
  matching.py        Geopard, dtw_match, candidate ordering
  selection.py       nearest_neighbours, crop_track
  dtw.py             accumulated_cost_matrix, dtw_computation
  interpolation.py   interpolate
  gpx.py             load_track, validate_track, duplicate detection
  geodesy.py         spheroid_point_distance, Coordinates
  regions.py         create_polygon
  response.py        GeopardResponse, MatchFlag, log_response
  plotting.py        palette, style, figures
  exceptions.py      GeopardException
  settings.py        PROJECT_ROOT — for the shipped examples only
  geopard.py         backwards-compatible import location
```

Dependencies run one way: `matching` knows about everything, `selection` knows
`geodesy`, `dtw` knows `interpolation`, and the leaves know nothing. Nothing
imports `matching`, which is what keeps the pipeline a sequence rather than a
web.

`settings.PROJECT_ROOT` deliberately sits outside that graph. It resolves the
repository root from a source file's location and is meaningless in an installed
wheel, so it exists for `examples/` and nothing else. The tests resolve their own
paths in `conftest.py` precisely so that a change here cannot repoint the suite
at different data without anybody noticing.
