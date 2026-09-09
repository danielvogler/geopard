"""Where the repository lives, for the shipped examples.

Only useful when running from a checkout -- an installed wheel has no example
data next to it. The tests deliberately do not use this; they resolve their own
paths from ``tests/conftest.py``, so a change here cannot silently repoint the
suite at the wrong files.
"""

from __future__ import annotations

from pathlib import Path

#: Repository root, resolved from this file at ``src/geopard/settings.py``.
PROJECT_ROOT_PATH = Path(__file__).resolve().parent.parent.parent

# An editable install inside the project's own virtualenv resolves through
# `.venv/site-packages/...`, so the walk above lands inside the venv rather
# than at the repository that contains it. Cut back to just above `.venv`.
if ".venv" in PROJECT_ROOT_PATH.parts:
    index = PROJECT_ROOT_PATH.parts.index(".venv")
    PROJECT_ROOT_PATH = Path(*PROJECT_ROOT_PATH.parts[:index])

#: The same path as a string, which is how the examples concatenate it.
PROJECT_ROOT = str(PROJECT_ROOT_PATH)

#: Example GPX tracks shipped with the repository.
GPX_DATA_DIR = PROJECT_ROOT_PATH / "data" / "gpx_files"

#: Example start/finish region polygons shipped with the repository.
POLYGON_DATA_DIR = PROJECT_ROOT_PATH / "data" / "csv_polygon_files"
