# ROSALIA Repository Status Report

Date: October 26, 2026
Repository: `Borlaff/ROSALIA`
Scope: packaging, runtime behavior, environment dependencies, and project maturity review

---

## Executive Summary

ROSALIA is a research-focused Python package for modeling stray light and sky background in Roman and Hubble data products. The project has moved toward a conventional Python packaging layout with a `pyproject.toml` configuration, documentation for installation, and a public distribution name (`rosalia-wfi`).

The repository is more organized than a loose research script collection, but several engineering issues remain. The most important are version drift, environment-dependent runtime assumptions, and a few legacy packaging and safety patterns that should be cleaned up before treating the package as production-ready.

---

## Current Package Metadata

The repository currently defines:

- Distribution name: `rosalia-wfi`
- Python requirement: `>=3.9`
- Build backend: `setuptools.build_meta`
- Package version in `pyproject.toml`: `1.3.9.dev0`

The runtime package metadata in `rosalia/__init__.py` does not match this exactly:

- `__version__` fallback: `1.4.0`

This is the main visible metadata inconsistency in the repo. A single canonical version source should be used for both build metadata and runtime import metadata.

---

## Installation and Runtime Notes

The project README includes installation guidance for a clean Conda environment and a pip installation path:

- `pip install rosalia-wfi`
- `ROSALIACACHE` must be set for the cache and calibration files used by the main sky and stray-light workflows

This is a valid operational requirement for the science pipeline, but the runtime behavior should be more explicit and fail gracefully when the cache directory is missing.

---

## Positive Findings

The repo shows several strengths:

- A modern `pyproject.toml` is present
- The project has a clear Python package layout under `rosalia/`
- Scientific dependencies are centralized in one dependency list
- Documentation and installation instructions are present and publicly oriented
- The package is organized around a defined import namespace rather than only standalone scripts

---

## Remaining Issues

### 1. Version drift

The distribution metadata and runtime metadata disagree:

- `pyproject.toml`: `1.3.9.dev0`
- `rosalia/__init__.py`: fallback `1.4.0`

This should be resolved by using a single version source of truth.

### 2. Import-time behavior

`rosalia/__init__.py` eagerly imports many submodules at import time. This improves convenience for users doing `import rosalia as rs`, but it also increases import time and can make the package fragile when optional scientific packages or environment variables are not configured.

### 3. Environment-variable coupling

Several modules rely on `ROSALIACACHE` being set. This is acceptable for the actual pipeline, but it should be validated with explicit error messages rather than causing opaque failures.

### 4. Legacy shell and filesystem patterns

Some code paths still use direct shell calls and file-side effects rather than safer Python APIs. For maintainability and portability, shell calls should be tightened and file handling should be more disciplined.

### 5. Packaging cleanup

The project still has some packaging details worth cleaning up:

- `Cython` is included in the build-system requirements even though the project does not appear to ship compiled extension modules
- package data declarations are broad and may be redundant
- CLI entry points are not expressed as modern `[project.scripts]` metadata

### 6. Quality and test posture

The repository includes a test area and scientific notebooks, but the project would benefit from a clearer separation between:

- runtime dependencies
- documentation dependencies
- development/test tooling
- calibration/cache assets

A stronger CI workflow would improve the long-term maintainability of the package.

---

## Overall Assessment

ROSALIA is in a healthier state than a typical research prototype repository: it has a Python package structure, packaging metadata, installation instructions, and a real project identity. The repo is moving in the right direction.

The main work left is not a rebuild of the project, but a stabilization pass focused on:

- version consistency
- explicit runtime dependency checks
- safer shell and file operations
- cleaner packaging metadata
- stronger quality automation

These are tractable improvements and would significantly raise the reliability of the package for scientific users and contributors.

---

## Recommended Next Actions

1. Standardize the package version across build and runtime metadata.
2. Make `ROSALIACACHE` and other required runtime dependencies fail with clear messages.
3. Replace unsafe shell operations with safer subprocess/pathlib patterns where practical.
4. Clean up package data and define CLI entry points in `pyproject.toml`.
5. Separate test and docs dependencies from runtime dependencies.
6. Add a lightweight CI workflow for import smoke tests and unit tests.

This will bring ROSALIA closer to a stable, reproducible, and professional scientific Python distribution.
