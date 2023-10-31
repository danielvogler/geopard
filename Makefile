SHELL := /bin/bash

VENV := .venv

POETRY_VERSION='$(shell poetry --version)'
POETRY_VERSION_SHORT=$(shell echo $(POETRY_VERSION) | sed 's/[^0-9.]*//g')

.DEFAULT_GOAL = help

help:
	@echo "---------------HELP-----------------"
	@echo "Set up the project: make setup"
	@echo "Activate the venv: make venv"
	@echo "Test the project: make test"
	@echo "Clean the project: make clean"
	@echo "------------------------------------"

PYTHON = poetry run python

all: setup

.PHONY: data

venv:
	poetry shell

check_poetry_version:
ifneq "1.3.2" "$(word 1, $(sort 1.3.2 $(POETRY_VERSION_SHORT)))"
	$(warning "Recommended: Update Poetry version to >= 1.3.2")
endif

setup: install check_poetry_version

install:
	@echo "************  Setup poetry ************"
	pip install poetry
	poetry config virtualenvs.in-project true
	poetry install

	@echo "************  Install pre-commit ************"
	pip install pre-commit
	pre-commit install

test:
	pre-commit run -a
	${PYTHON} -m pytest

clean:
	@echo "************  Clean project ************"

	@echo "************  Remove venv ************"
	rm -rf $(VENV)

	@echo "************  Remove *.pyc ************"
	find . -type f -name '*.pyc' -delete
