# ROSALIA Codebase Quality & Standards Audit Report

**Date:** September 16, 2026  
**Repository:** `Borlaff/ROSALIA`  
**Scope:** Full repository scan including folder structure, package layout, packaging/pip configuration, test infrastructure, and Python source code across all modules.

---

## 1. Executive Summary

A comprehensive automated and manual code audit of the **ROSALIA** repository revealed numerous deviations from standard Python packaging conventions, PEP standards (PEP 8, PEP 517/518/621), software architecture best practices, security standards, andv resource management rules.

### Key Metrics
- **Total Python Modules Scanned:** 34 files (~13,000+ lines of Python code)
- **Functions with Type Annotations:** 0 / 322 (0.0%)
- **Functions with Docstrings:** 78 / 322 (24.2%)
- **Classes with Docstrings:** 11 / 30 (36.7%)
- **`print()` calls in library code:** 632 instances across 22 modules (instead of standard `logging`)
- **Bare `except:` blocks:** 27 instances silently masking exceptions
- **`fits.open()` / `open()` calls without context managers (`with`):** 52 instances causing file descriptor leaks
- **`os.system()` invocations:** 15 instances executing shell commands via unescaped string concatenation
- **Hardcoded local user paths (`/Users/aborlaff/...`):** 8 instances across multiple modules

---

## 2. Folder Structure and File Locations

### 2.1 Stray and Cluttered Root Files
1. **Redundant Root `__init__.py`:**
   - Location: `/__init__.py`
   - Issue: Having an `__init__.py` at the repository root alongside the `rosalia/` directory treats the top-level repository folder as a Python package, leading to module collision, import ambiguity, and packaging bugs.

2. **Non-Package Assets at Repository Root:**
   - Location: `/PAPER_I/` (Contains LaTeX papers, proceedings, bibtex files, PDFs, EPS, and `.zip` files)
   - Issue: Academic paper drafts and build artifacts are committed directly in the root of the software package repository rather than in a separate repository or dedicated `papers/` / documentation branch.

3. **Heavy Binary Files Tracked Directly in Git:**
   - Location: `/images/` (Contains large uncompressed GIFs and MP4 videos totaling ~86 MB, e.g., `rosalia_star_on_the_edge.gif` [37.8 MB], `earth_moon_sun_from_L2.gif` [26.8 MB], `rosalia_loading.gif` [5.6 MB], `rosalia_loading.mp4` [1.2 MB])
   - Issue: Bloats the repository size for all contributors and CI pipelines. Large assets should be hosted externally, tracked via Git LFS, or compiled dynamically.

4. **Clutter and OS/IDE Artifacts in Source Control:**
   - `.DS_Store` files are present in `/`, `/docs/`, `/images/`, `/rosalia/`.
   - `.vscode/` configuration directory is present in the repository root.
   - Built package artifacts (`/build/`, `/dist/`, `/ROSALIA.egg-info/`, `/ROSALIA_wfi.egg-info/`) are committed or left unignored in the source directory.

5. **Notebooks Directory Organization:**
   - Location: `/notebooks/`
   - Issues:
     - Contains raw downloaded folders: `notebooks/tutorials/raw.githubusercontent.com/...`
     - Contains leftover workspace config: `notebooks/tutorials/ROSALIA.code-workspace`
     - Contains binary outputs: `notebooks/tutorials/output.np`
     - Contains duplicated files and scratch notebooks: `R1_Simulate_Roman_straylight-Copy1.ipynb`, `APT_1023_GhostsStrayLight_CAR171_modified.apt`

### 2.2 Package Subdirectory Structure (`rosalia/`)
1. **Case-Sensitivity Collision:**
   - Location: `rosalia/CORE/` vs `rosalia/core.py`
   - Issue: Directory `rosalia/CORE/` uses ALL_CAPS while `rosalia/core.py` uses lowercase. On case-insensitive file systems (macOS, Windows), this can cause namespace collisions and import shadowing.

2. **Redundant Asset Storage inside Package:**
   - Location: `rosalia/images/raccoon3_logo.png`
   - Issue: Duplicates image assets already present in `/images/`.

3. **Non-Standard Resource Storage:**
   - Location: `rosalia/CORE/detector/RST_WFI_SCA_Subarray_locations.csv`, `rosalia/style/*.mplstyle`
   - Issue: Data files are placed in ad-hoc subdirectories without standard Python packaging resource layout (such as `importlib.resources`).

### 2.3 Test Suite Layout
1. **Split and Misplaced Test Infrastructure:**
   - Location: `rosalia/tests.py` vs `tests/psf.py`
   - Issue: The primary test suite is packaged inside the production library code (`rosalia/tests.py`), while the root `tests/` directory only contains a single file (`tests/psf.py`) which is an unasserted demo script copied from a Jupyter notebook rather than a valid unit test.

### 2.4 CLI Executables (`bin/`)
1. **Legacy `bin/` Directory:**
   - Location: `bin/exposure-inspector`, `bin/rosalia-sky`, `bin/rosalia-stray`, `bin/run-rosalia-tests`, `bin/test-rosalia-parser`
   - Issue: Command-line entry points are placed in a standalone `bin/` directory without Python extensions (`.py`) and installed via legacy `scripts=` rather than standard `[project.scripts]` console entry points.
2. **Developer Scratch Script Committed in `bin/`:**
   - Location: `bin/test-rosalia-parser`
   - Issue: Contains placeholder code described as *"This is a developer canvas to test parsing arguments in Python"*.
3. **Syntax / Header Glitches in Executable:**
   - Location: `bin/exposure-inspector` (Lines 1 & 5)
   - Issue: Contains a duplicate shebang line (`#!/usr/bin/env python`) on line 5 after an import statement.

---

## 3. Pip Installation, Packaging & Configuration

### 3.1 Metadata & Version Inconsistencies
1. **Version Divergence:**
   - `pyproject.toml` defines `version = "1.2.4"`
   - `rosalia/__init__.py` defines `__version__ = "1.2.1"`
   - Result: Users installing the package see version `1.2.4` in `pip list`, but `rosalia.__version__` reports `1.2.1`.

2. **Package Naming Inconsistency:**
   - `pyproject.toml` defines `name = "rosalia-wfi"`
   - Root egg-info contains `ROSALIA.egg-info` and `ROSALIA_wfi.egg-info`
   - `setup.py` packages `rosalia`
   - Result: Discrepancy between PyPI distribution name (`rosalia-wfi`) and import name (`rosalia`).

### 3.2 `setup.py` vs `pyproject.toml` Conflicts
1. **Superfluous Cython Build Dependency:**
   - In `setup.py`: Imports `from Cython.Build import cythonize` and sets `Options.docstrings = True`, `Options.annotate = False`.
   - In `pyproject.toml`: `build-system.requires = ["setuptools>=61.0", "Cython >=0.29.21", "wheel"]`.
   - Issue: The repository contains **zero** `.pyx`, `.pxd`, or C/C++ extensions. Forcing Cython in `build-system` adds unnecessary installation overhead, compile-time friction, and potential build failures on minimal environments.

2. **Broken / Redundant Package Data Declarations:**
   - In `setup.py`: `package_data={'rosalia': ['rosalia/style/*.mplstyle']}` (incorrect syntax: paths inside `package_data['rosalia']` must be relative to the `rosalia` package folder, e.g., `'style/*.mplstyle'`).
   - In `pyproject.toml`: `[tool.setuptools.package-data]` defines `rosalia = ["*.mplstyle", "*.csv", "*.txt", "rosalia/style/*", "**/*.mplstyle", "images/*", "bin/*"]`.
   - Result: Conflicting declarations between `setup.py` and `pyproject.toml`.

3. **Missing `[project.scripts]` in `pyproject.toml`:**
   - CLI tools are wired via `setup(scripts=scripts)` in `setup.py` instead of modern, cross-platform `[project.scripts]` entry points in `pyproject.toml`.

### 3.3 Dependency Management Discrepancies
There is significant divergence across `pyproject.toml`, root `requirements.txt`, and `docs/requirements.txt`:

| Dependency | In `pyproject.toml` | In `requirements.txt` | In `docs/requirements.txt` | Issue / Note |
|---|---|---|---|---|
| `astropy` | ❌ **Missing!** | ✅ Present | ✅ Present | **Critical:** Core dependency missing from `pyproject.toml`! |
| `scipy` | Unpinned | `scipy` + `scipy>=1.18.0` | `scipy` | Invalid duplicate and non-existent version constraint |
| `pytest` | In core `dependencies` | In core | In core | Test dependency in core install |
| `sphinx`, `sphinx-rtd-dark-mode` | In core `dependencies` | In core | In core | Docs dependency in core install |
| `twine`, `build` | In core `dependencies` | In core | ❌ | Packaging tooling in core install |
| `jupyterlab`, `ipython` | In core `dependencies` | `ipython` | `ipython` | Interactive dev tooling in core install |
| `drizzlepac`, `zodipy`, `s3fs`, `crds` | ❌ **Missing** | ✅ Present | Mixed | Missing from `pyproject.toml` |
| `romanisim`, `cartopy`, `asdf`, `galsim` | ✅ Present | ❌ **Missing** | ❌ **Missing** | Missing from `requirements.txt` |

### 3.4 Packaging & Deployment Scripts
1. **`upload_to_pypi.sh`:**
   - Line 2: `rm dist/*` fails if `dist/` is empty or non-existent.
   - Lacks error handling (`set -e`), credential protection, or dry-run checks.
2. **`upgrade_documentation.sh`:**
   - Lines 2-3: Runs `sphinx-build -b html source build/html` twice consecutively without explanation.
3. **`.gitignore` Omissions:**
   - Fails to ignore `.DS_Store`, `.vscode/`, `dist/`, `build/`, `*.egg-info/`, `.apt`, `output.np`, and temporary FITS files.

---

## 4. Code Quality, Anti-Patterns & PEP Violations

### 4.1 Critical Runtime & Crash Bugs on Import
1. **Module-Level Unhandled Environment Variable Access (`KeyError`):**
   - Direct access to `os.environ["ROSALIACACHE"]`, which crashes immediately if the variable is unset:
     - `rosalia/telescopes.py:32`: `ndi_db = pd.read_csv(os.environ["ROSALIACACHE"] + "/CORE/ndi_HST_legacy_bely2003.csv")`
     - `rosalia/telescopes.py:311`: `with sf_api.load.open(os.environ["ROSALIACACHE"] + "/CORE/HST_TLE_history_sat000020580.txt") as f:`
     - `rosalia/psf.py:51`: `psf_archive = os.environ["ROSALIACACHE"] + "/CORE/PSF_ARCHIVE/"`
   - **Impact:** Because `rosalia/__init__.py` imports `rosalia.telescopes` and `rosalia.psf`, executing `import rosalia` immediately crashes with `KeyError: 'ROSALIACACHE'` unless the user has manually configured that specific environment variable.

2. **Hardcoded User-Specific Absolute Paths:**
   - `rosalia/telescopes.py:373`: `exp_name = "/Users/aborlaff/NASA/SPARKLES/notebooks/SATELLITES/HST_ACS_mock_es.fits"`
   - `rosalia/telescopes.py:421`: `psf_hdu = fits.open("/Users/aborlaff/NASA/SPARKLES/notebooks/SATELLITES/psf_tinytim00_psf_scaled_cropped.fits")`
   - `rosalia/telescopes.py:732`: `exp_name = "/Users/aborlaff/NASA/SPARKLES/notebooks/SATELLITES/CSST_mock_es.fits"`
   - `rosalia/telescopes.py:822-824`: References to `/Users/aborlaff/NASA/PAPER_SPARKLES_1/...`
   - `rosalia/telescopes.py:949`: `exp_name = "/Users/aborlaff/NASA/PAPER_SPARKLES_1/notebooks/SATELLITES/ARRAKIHS_mock_es.fits"`
   - `rosalia/tests.py:7`: `#sys.path.append("/Users/aborlaff/NASA/STRAYCOR/")`
   - **Impact:** Functions calling these lines fail on any environment other than the original author's machine.

3. **Overriding Python Special Dunder Attributes:**
   - `rosalia/__init__.py:6`: `__name__ = "ROSALIA"`
   - Overwriting the built-in `__name__` module attribute breaks module introspection and standard tooling.

4. **Circular / Redundant Self-Import:**
   - `rosalia/__init__.py:13`: `import rosalia` inside `rosalia/__init__.py`.

5. **Syntax Warnings:**
   - `rosalia/plots.py:731-733`: Invalid escape sequences `\/`, `\ `, `\_` in non-raw string literals trigger `SyntaxWarning` on Python 3.12+.

---

### 4.2 Security Vulnerabilities & Dangerous Shell Execution
1. **Unsafe `os.system()` with String Concatenation:**
   Instead of using Python's standard `pathlib`, `os.remove`, `shutil`, or `subprocess.run`, raw shell commands are constructed via string concatenation:
   - `rosalia/astrometry.py:158`: `os.system("rm " + rootname + "_crclean.fits")`
   - `rosalia/gaia.py:229-230`: `os.system('mkdir ' + os.environ["ROSALIACACHE"] + '/cache')`
   - `rosalia/hst.py:23, 28`: `os.system("mkdir " + filter_name)`, `os.system("mv " + selected_exposure + " " + filter_name)`
   - `rosalia/skywalker.py:137, 142`: `os.system("mkdir " + filter_name)`, `os.system("mv " + selected_exposure + " " + filter_name)`
   - `rosalia/telescopes.py:560`: `os.system("rm " + outname)`
   - `rosalia/roman.py:756`: `os.system("romanisim-make-image --radec " + str(ra) + " " + str(dec) + ...)`
   - `rosalia/core.py:786`, `rosalia/correct.py:406`, `rosalia/utils.py:550`: `os.system("swarp -dd > swarp.conf")`
   - `rosalia/utils.py:1187`: `os.system("astnoisechisel --tilesize=5,5 -K -h" + str(ext) + " " + file_name)`
   - **Risks:** Shell injection vulnerabilities, failure on paths containing spaces, OS incompatibility (Windows), and unhandled exit codes.

2. **Insecure External Network Endpoints:**
   - `rosalia/skysurf.py:51`: Uses unencrypted `http://skysurf.asu.edu/sky_measurements/` instead of HTTPS.
   - `rosalia/sky.py:637`: Fetches refdata dynamically at runtime from raw GitHub links.

---

### 4.3 Resource Management & File Leaks
1. **Unclosed File Descriptors (`fits.open` & `open`):**
   `astropy.io.fits.open()` is called without `with` context managers across almost every module (52+ instances), including:
   - `rosalia/correct.py`: lines 421, 422, 619
   - `rosalia/gnu.py`: lines 64, 100
   - `rosalia/hst.py`: line 39
   - `rosalia/plots.py`: lines 98, 206, 341, 478, 502
   - `rosalia/psf.py`: lines 97, 163, 282, 381, 402, 405, 437, 520, 837
   - `rosalia/roman.py`: lines 115, 147, 721, 730
   - `rosalia/sky.py`: lines 40, 67, 68
   - `rosalia/telescopes.py`: lines 375, 405, 421, 559, 734, 826-828, 876, 951
   - `rosalia/utils.py`: lines 75, 163, 295, 299, 391, 562, 592, 649, 731, 741, 784, 816, 844, 949, 1003, 1013, 1021, 1100, 1199, 1234, 1236, 1298, 1610, 1735, 1782, 1953, 1967, 1983
   - **Impact:** Memory leaks when opening memory-mapped arrays (`memmap=True`), file handle exhaustion, and file lock contentions.

2. **Temporary Files Dumped in Working Directory:**
   - Multiple functions create hardcoded temporary files directly in the user's current working directory (e.g., `swarp.conf`, `subprocess.out`, `coadd.fits`) rather than using Python's standard `tempfile` module.

---

### 4.4 Global State & Side Effects
1. **Global Warning Suppression on Import:**
   - `rosalia/gaia.py:27`: `warnings.filterwarnings('ignore')`
   - `rosalia/plots.py:23`: `warnings.filterwarnings('ignore')`
   - `rosalia/psf.py:28`: `warnings.filterwarnings('ignore')`
   - **Impact:** Unconditionally silences all Python warnings across the entire interpreter session upon importing `rosalia`.

2. **Global Matplotlib Style Mutation:**
   - `rosalia/plots.py:204`: `plt.style.use(os.path.dirname(rs.__file__) + "/style/nature_style.mplstyle")`
   - **Impact:** Plot functions overwrite global matplotlib stylesheet settings for the entire active session.

3. **Eager Submodule Loading:**
   - `rosalia/__init__.py:14-37`: Imports all 24 submodules at top level on initialization.
   - **Impact:** Extremely slow startup and fragile imports.

---

### 4.5 Exception Handling & Error Suppression
1. **27 Bare `except:` Statements:**
   - Found in: `correct.py` (L201, L210), `horizons.py` (L84), `hst.py` (L46), `plots.py` (L102, L446), `psf.py` (L842), `sky.py` (L81), `skysurf.py` (L25, L34, L43), `skywalker.py` (L92), `telescopes.py` (L122), `utils.py` (L281, L341, L396, L403, L416, L451, L469, L481, L546, L906, L1025, L1110, L1305, L1785).
   - **Anti-Pattern:** Intercepts `KeyboardInterrupt` and `SystemExit`, making scripts impossible to abort with `Ctrl+C` and masking fatal runtime errors.

2. **Absence of Custom Exception Hierarchy:**
   - No custom domain exceptions (e.g., `RosaliaError`, `CalibrationError`, `CoordinateError`) are defined across the codebase.

---

### 4.6 Logging & Print Statements
1. **632 `print()` Calls in Library Code:**
   - Raw `print()` statements are scattered across 22 modules instead of standard `logging.getLogger(__name__)`.

---

### 4.7 Function Design & Code Smells
1. **Mutable Default Arguments:**
   - `rosalia/attitude.py:188`: `def create_subplot_grid(..., projection_list=['hammer', 'hammer', ...])`
   - `rosalia/telescopes.py:393`: `def get_HST_ACS_psf(..., position=[1000, 1000])`
   - `rosalia/utils.py:1246`: `def create_custom_wcs(..., pixel_scale=[-0.11 / 3600.0, 0.11 / 3600.0])`

2. **Prolific Deferred / Function-Level Imports (100+ instances):**
   - Standard modules (`os`, `re`, `multiprocessing`, `copy`, `subprocess`) and third-party packages (`numpy`, `scipy`, `pandas`, `astropy`, `tqdm`, `galsim`) are imported inside function definitions rather than at the module top level.

3. **Non-Standard Aliasing:**
   - `rosalia/utils.py:109, 122`: `import dill as pickle`

4. **Monolithic Modules:**
   - `rosalia/utils.py` (~2,000 lines), `rosalia/psf.py` (1,178 lines), `rosalia/telescopes.py` (1,139 lines), `rosalia/attitude.py` (1,399 lines).

5. **Unresolved TODO Comments:**
   - `rosalia/astrometry.py:14`: `# TODO: We need to find a way to generalize this.`
   - `rosalia/correct.py:600`: `# TODO: Organize the return in a coherent way.`
   - `rosalia/psf.py:33`: `## TODO: Compute the scale factor for the object at (x,y)=(53,69) for...`

---

### 4.8 PEP 8 Naming, Typing & Style Violations
1. **Class Naming (PEP 8 requires `CapWords` / `PascalCase`):**
   - `rosalia/core.py:16`: `class exposure:` (lowercase)
   - `rosalia/ndi.py:100`: `class ndi_euclid:` (snake_case)
   - `rosalia/plots.py:32`: `class style:` (lowercase)
   - `rosalia/plots.py:615`: `class ascii_progress_focal_plane:` (snake_case)
   - `rosalia/sky.py:612`: `class background:` (lowercase)
   - `rosalia/sky.py:857`: `class compute_visibility:` (snake_case)

2. **Function Naming (PEP 8 requires `snake_case`):**
   - `bootima.py:27`: `does_it_fit_in_RAM` (camelCase + uppercase)
   - `psf.py:1069`: `getWCS_galsim_dict_style`
   - `psf.py:1080`: `find_SCA_for_a_target`
   - `sso.py:108`: `get_SS0s_loc_magnitude` (contains typo: number `0` instead of letter `O`)
   - `telescopes.py:489`: `DEPRECATED_get_bestPA`

3. **Constants Naming in `rosalia/constants.py`:**
   - Inconsistent casing: `MJysr_to_Jyarcsec2`, `mapp_sun_V_AB`, `TwoMASS_fnu0_J`, `WISE_W1_delta_mag_Vega_to_AB`.

4. **Type Hints & Docstrings:**
   - **Type hints:** 0 out of 322 functions (0.0%).
   - **Docstrings:** Only 24.2% of functions and 36.7% of classes have docstrings.

---

### 4.9 Test Suite Anti-Patterns
1. **Tests Packaged in Production Library:**
   - `rosalia/tests.py` is shipped inside the production distribution.
2. **Interactive GUI Calls in Automated Tests:**
   - `rosalia/tests.py:142`: `plt.show(block=False)` pops open matplotlib GUI windows during test runs.
3. **Live Unmocked Network Queries:**
   - `test_magnitude_conversion_gaia()` and `test_photometry_superstars()` query remote TAP databases over HTTP without caching or mocks.
4. **Non-Standard Assertions and Return Values:**
   - Test functions use stdout prints and `return(True)` rather than standard pytest assertions.
5. **Disconnected Test Runners:**
   - `pytest.ini_options` is configured to scan for `test_*.py`, which misses `rosalia/tests.py` and `tests/psf.py`.
   - `bin/run-rosalia-tests` uses a custom script printing ASCII art banners.

---

## 5. Summary Table of Scanned Files

| File | Lines | Primary Non-Standard Practices Found |
|---|---|---|
| `__init__.py` (root) | 7 | Stray root package marker; causes root import ambiguity |
| `pyproject.toml` | 93 | Version mismatch; Cython build requirement; missing `astropy`; dev deps in core |
| `setup.py` | 17 | Unused Cython imports; broken `package_data` relative paths; legacy `scripts=` |
| `requirements.txt` | 37 | Invalid `scipy>=1.18.0` duplicate; missing core packages; dev/doc tool pollution |
| `docs/requirements.txt` | 29 | Out of sync with `pyproject.toml` and root `requirements.txt` |
| `upload_to_pypi.sh` | 4 | Fragile `rm dist/*`; lacks safety checks and error handling |
| `upgrade_documentation.sh` | 4 | Duplicate `sphinx-build` invocation; hardcoded relative paths |
| `bin/exposure-inspector` | 33 | Duplicate shebang header; non-standard script packaging |
| `bin/rosalia-sky` | 91 | Script in `bin/` without entry point; uses `os.system` via utils |
| `bin/rosalia-stray` | 63 | Script in `bin/` without entry point |
| `bin/run-rosalia-tests` | 39 | Custom script runner; ASCII art instead of standard test runner |
| `bin/test-rosalia-parser` | 25 | Developer playground scratch script committed in production |
| `rosalia/__init__.py` | 37 | Version mismatch (`1.2.1`); `__name__` overwrite; circular self-import; eager loading |
| `rosalia/albedo.py` | 61 | Missing module docstring; non-standard function naming |
| `rosalia/astrometry.py` | 261 | `os.system("rm ...")`; unresolved TODO; missing module docstrings |
| `rosalia/attitude.py` | 1399 | Mutable default argument; 16 `print()` calls; missing docstrings |
| `rosalia/bootima.py` | 487 | 54 `print()` calls; non-standard function naming (`does_it_fit_in_RAM`) |
| `rosalia/constants.py` | 75 | Inconsistent constant naming; mixed units |
| `rosalia/core.py` | 917 | `os.system("swarp ...")`; `class exposure` (lowercase); 38 `print()` calls |
| `rosalia/correct.py` | 623 | Unclosed FITS files; bare `except:`; `os.system()`; unresolved TODO |
| `rosalia/detectors.py` | 275 | Missing docstrings; non-standard uppercase function naming |
| `rosalia/gaia.py` | 491 | `warnings.filterwarnings('ignore')`; `os.system('mkdir ...')`; unclosed files |
| `rosalia/gnu.py` | 116 | Unclosed FITS files; missing module docstring |
| `rosalia/horizons.py` | 331 | Bare `except:`; default URL parameter pointing to external GitHub repo |
| `rosalia/hst.py` | 96 | `os.system("mkdir/mv")`; bare `except:`; unclosed FITS files |
| `rosalia/irsa.py` | 108 | Missing docstrings; `print()` statements |
| `rosalia/ndi.py` | 360 | `class ndi_euclid` (snake_case); 23 `print()` statements |
| `rosalia/plots.py` | 737 | Global `warnings.filterwarnings('ignore')`; `plt.style.use()` side effect; SyntaxWarnings |
| `rosalia/psf.py` | 1178 | Crash on import (`os.environ["ROSALIACACHE"]`); 74 `print()` calls; unclosed FITS |
| `rosalia/render.py` | 26 | Missing module docstring; undocumented functions |
| `rosalia/roman.py` | 762 | `os.system("romanisim-make-image ...")`; unclosed FITS files; 35 `print()` calls |
| `rosalia/sky.py` | 946 | `class background` (lowercase); bare `except:`; hardcoded remote GitHub URLs |
| `rosalia/skysurf.py` | 113 | Insecure HTTP URLs; 3 bare `except:` statements |
| `rosalia/skywalker.py` | 258 | `os.system("mkdir/mv")`; bare `except:`; unclosed FITS files |
| `rosalia/sso.py` | 158 | Typo in function name `get_SS0s_loc_magnitude`; missing docstrings |
| `rosalia/telescopes.py` | 1139 | Crash on import (`os.environ["ROSALIACACHE"]`); 8 hardcoded `/Users/aborlaff/` paths |
| `rosalia/tests.py` | 235 | In-package tests; live Gaia network queries; `plt.show()` in tests; hardcoded paths |
| `rosalia/utils.py` | 1997 | Monolithic "god" module; 16 bare `except:` blocks; 91 `print()` calls; 25+ unclosed FITS |
| `tests/psf.py` | 28 | Unasserted demonstration script in `tests/`; no test framework |
