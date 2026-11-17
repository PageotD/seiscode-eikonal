# Architecture & Development Guide

This document describes the software design, Rust workspace structure, testing strategy, and development workflows for `seiscode-eikonal`.

---

## 1. Workspace Layout

The repository uses a Cargo workspace with two distinct crates:

```
seiscode-eikonal/
├── crates/
│   ├── seiscode-eikonal-core/       # Pure Rust algorithmic core
│   │   ├── Cargo.toml
│   │   └── src/
│   │       ├── lib.rs              # Public exports and error types
│   │       ├── grid.rs             # 2D Grid structure and index mapping
│   │       ├── stencil.rs          # 1st and 2nd-order Eikonal update formulas
│   │       ├── narrowband.rs       # Min-heap priority queue with lazy deletion
│   │       ├── solver.rs           # Fast Marching Method loop
│   │       └── ray.rs              # RK2 ray backtracking and gradient interpolation
│   └── seiscode-eikonal-py/         # Thin Python extension layer
│       ├── Cargo.toml
│       └── src/
│           ├── lib.rs              # PyO3 bindings and NumPy conversions
│           └── error.rs            # Rust to Python exception mapping
├── python/                         # Python tests and scripts
│   └── tests/
├── examples/                       # Geophysical benchmark scripts
│   └── plot_rays_two_layer.py
├── docs/                           # Documentation
├── Taskfile.yml                    # Local command orchestration
└── pyproject.toml                  # Python package metadata (Maturin backend)
```

---

## 2. Design Principles & Invariants

1. **Separation of Concerns**:
   - `seiscode-eikonal-core` has **zero dependencies on PyO3 or Python**. It can be compiled, benchmarked, and used independently in pure Rust environments.
   - `seiscode-eikonal-py` contains **zero numerical logic**. It is strictly responsible for validating NumPy buffers and calling core functions.

2. **No Unsafe Code in Core**:
   - `seiscode-eikonal-core` is annotated with `#![deny(unsafe_code)]`.
   - Any raw pointer manipulation or unsafe interop is isolated to the Python binding layer where NumPy C buffers are accessed, fully documented with safety invariants.

3. **Zero Panics (`Result<T, SeikonalError>`)**:
   - Public APIs never panic on invalid inputs (e.g. non-positive velocities, negative grid spacings, out-of-bounds sources or receivers).
   - All errors are typed with `SeikonalError` and returned as `Result<T, SeikonalError>`, which PyO3 maps to standard Python exceptions (`ValueError`, `IndexError`, `RuntimeError`).

4. **Double Precision Everywhere (`f64`)**:
   - Geodesic distance tracking, numerical discriminant calculation, and ray backtracking require high precision to prevent numerical drift.
   - Using `f64` across the stack aligns natively with Python standard floats and NumPy `float64` without generic overhead.

---

## 3. Testing Strategy

The repository follows a strict Test-Driven Development (TDD) methodology with layered validation:

### Rust Unit & Integration Tests
- **Stencil Tests (`stencil.rs`)**:
  - Validates 1D and 2D roots against exact pencil-and-paper quadratic solutions.
  - Tests 2nd-order upstream combinations and boundary fallbacks.
- **Narrowband Tests (`narrowband.rs`)**:
  - Tests min-heap insertion, order preservation, and lazy deletion of stale candidate nodes.
- **FMM Integration Tests (`solver.rs`)**:
  - **Homogeneous Model**: Compares computed traveltime against exact analytical Euclidean distance $T(\mathbf{x}) = \frac{\|\mathbf{x} - \mathbf{x}_s\|}{v}$. Checks that max relative error is $< 0.5\%$.
  - **Two-Layer Medium (Snell's Law)**: Tests refractions, critical angles, and head waves (conical waves) along the interface.
- **Ray Backtracking Tests (`ray.rs`)**:
  - Tests straight-line raypaths in homogeneous media.
  - Verifies ray bending across velocity interfaces matching Snell-Descartes refraction law:
    $$\frac{\sin \theta_1}{v_1} = \frac{\sin \theta_2}{v_2}$$

### Python Tests (`python/tests/`)
- Run via `pytest python/tests/`.
- Verifies zero-copy buffer passing, non-contiguous array handling, Fortran-order conversion, and parameter validation exceptions.

---

## 4. Local Development Workflow

We use [Taskfile](https://taskfile.dev) as a single source of truth for all build and verification tasks:

```bash
# Code Formatting & Style
task fmt            # Auto-format Rust code (cargo fmt)
task lint           # Check formatting and Clippy warnings (-D warnings)

# Testing
task test:unit          # Rust unit tests
task test:integration   # Rust integration tests
task python:dev         # Build and install editable wheel (pip install -e .)
task python:test        # Run pytest suite
task test               # Run all tests (Rust + Python)

# Security & Dependency Auditing
task audit:deny     # Check licenses and bans (cargo-deny)
task audit          # Check security vulnerabilities (cargo-audit)

# Full CI Verification
task ci             # Equivalent to the complete GitHub Actions PR pipeline
```

---

## 5. Continuous Integration & Release Pipeline

The CI/CD pipeline on GitHub Actions includes:
- **`pull_request_lint.yml`**: Validates Conventional Commits titles on PRs (e.g. `feat:`, `fix:`, `docs:`).
- **`ci.yml`**: Runs `task ci` on Linux with matrix testing for Python 3.11, 3.12, and 3.13.
- **`release_and_publish.yml`**:
  - Triggered on workflow dispatch or tag creation.
  - Calculates semantic version bumps using `git-cliff` with tags without `v` prefix (e.g. `0.2.0`).
  - Builds Linux (manylinux x86_64) and Windows (x86_64) wheels via `maturin-action`.
  - Publishes binary wheels and sdist to PyPI and GitHub Releases.
