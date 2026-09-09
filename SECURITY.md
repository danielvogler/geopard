# Security

## Reporting

Report a vulnerability privately through
[GitHub security advisories](https://github.com/danielvogler/geopard/security/advisories/new).
Please do not open a public issue for anything exploitable.

## What this repository is, security-wise

geopard is a library that reads two local files and returns a number. It holds no
credentials, opens no sockets, starts no server, and writes nothing except
figures a caller explicitly asks it to save. There is no deployment to harden.

That leaves a small, real surface:

- **It parses untrusted input.** A GPX file is XML, and a CSV region file is
  read with `csv.DictReader`. Parsing is delegated to `gpxpy` and the standard
  library rather than hand-rolled. If you match files uploaded by other people,
  the XML parser is the component to keep patched — `uv lock --upgrade` moves it.
- **It reads any path it is given.** `gpx_loading` and `create_polygon` open the
  filename passed to them, with no sandboxing, exactly like `open()`. An
  application that takes a path from a user is responsible for validating it
  before handing it over.
- **Memory is proportional to the input.** The accumulated cost matrix is
  1000 × 1000 floats per candidate pair, and the number of pairs grows with
  `radius`. A deliberately pathological activity can make a match slow. Bound
  `radius` and `min_trkps` if you accept files from strangers.
- **No telemetry, no network.** Nothing in the library makes an outbound
  request, and nothing reports usage anywhere.

## Supply chain

- Dependencies are locked in `uv.lock`, which is committed, and CI installs with
  `uv sync --frozen`.
- `pre-commit` runs `detect-private-key` and refuses large added files.
- Releases publish to PyPI through [trusted
  publishing](https://docs.pypi.org/trusted-publishers/) on a published GitHub
  release. No API token exists in this repository or its settings.
