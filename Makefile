.DEFAULT_GOAL := help

# Override on the command line, e.g. `make match GOLD=a.gpx ACTIVITY=b.gpx`.
GOLD     ?= data/gpx_files/tds_sunnestube_segment.gpx
ACTIVITY ?= data/gpx_files/tds_sunnestube_activity_25_25.gpx
RADIUS   ?= 7

.PHONY: help setup hooks test test-fast coverage lint format format-check \
        typecheck check match plots figures clean

help: ## Show this help
	@grep -hE '^[a-zA-Z_-]+:.*?## ' $(MAKEFILE_LIST) \
		| awk 'BEGIN {FS = ":.*?## "}; {printf "  \033[36m%-14s\033[0m %s\n", $$1, $$2}'

setup: ## Create the venv, install dependencies, install git hooks
	uv sync
	uv run pre-commit install --install-hooks

hooks: ## Run every pre-commit hook over the whole repository
	uv run pre-commit run --all-files

test: ## Run the whole test suite
	uv run pytest

test-fast: ## Skip the end-to-end matches over the shipped tracks
	uv run pytest -m "not slow"

coverage: ## Run the suite with a coverage report
	uv run pytest --cov --cov-report=term-missing

lint: ## Run Ruff lint
	uv run ruff check .

format: ## Format with Ruff
	uv run ruff format .

format-check: ## Check formatting without writing
	uv run ruff format --check .

typecheck: ## Run mypy
	uv run mypy

check: lint format-check typecheck test ## Everything CI runs

match: ## Match ACTIVITY against GOLD and plot the result
	uv run python examples/track_matching.py \
		--gold "$(GOLD)" --activity "$(ACTIVITY)" --radius $(RADIUS)

plots: ## Show the cropping and interpolation figures
	uv run python examples/gpx_track_cropping_and_interpolation.py \
		--gold "$(GOLD)" --activity "$(ACTIVITY)" --radius $(RADIUS)

figures: ## Regenerate the figures used in the README
	uv run python examples/gpx_track_cropping_and_interpolation.py \
		--gold "$(GOLD)" --activity "$(ACTIVITY)" --radius $(RADIUS) \
		--save-to docs/images/example
	uv run python examples/start_end_region_polygon_usage.py \
		--save-to docs/images/example_track_start-finish.png

clean: ## Remove the venv and build artefacts
	rm -rf .venv dist build .pytest_cache .mypy_cache .ruff_cache .coverage
	find . -type d -name __pycache__ -prune -exec rm -rf {} +
