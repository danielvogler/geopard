"""Generate a variable holding the project location."""
import os

PROJECT_ROOT = os.path.dirname(os.path.abspath(os.path.dirname(__file__)))

if "/.venv/" in PROJECT_ROOT:
    PROJECT_ROOT = PROJECT_ROOT[: PROJECT_ROOT.index("/.venv/")]
