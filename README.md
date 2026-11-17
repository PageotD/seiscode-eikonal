# seiscode-eikonal

[![CI](https://github.com/pageotd/seiscode-eikonal/actions/workflows/ci.yml/badge.svg)](https://github.com/pageotd/seiscode-eikonal/actions/workflows/ci.yml)
[![PyPI version](https://img.shields.io/pypi/v/seiscode-eikonal.svg)](https://pypi.org/project/seiscode-eikonal/)
[![Python versions](https://img.shields.io/pypi/pyversions/seiscode-eikonal.svg)](https://pypi.org/project/seiscode-eikonal/)
[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](https://opensource.org/licenses/MIT)

**seiscode-eikonal** is a high-performance, robust Fast Marching Method (FMM) eikonal solver and seismic ray tracer written in pure Rust with seamless Python bindings.

Part of the **seiscode** ecosystem, it computes first-arrival traveltime fields and traces seismic rays via gradient backtracking with high numerical stability and zero-copy NumPy interoperability.

---

## Key Features

- **Fast & Memory-Safe Numerical Core**: Pure Rust core (`crates/seiscode-eikonal-core`) built with `#![deny(unsafe_code)]` and zero panics.
- **High-Order Upwind Stencils**: Supports both 1st-order and 2nd-order upwind stencils (`order=1` or `order=2`, default `2`) with automatic boundary fallbacks.
- **Singularity Treatment**: Analytical spherical wave initialization in a $3 \times 3$ box around the source point to prevent point-source curvature errors.
- **Seismic Ray Tracing**: Backtracking from receivers to the source using Runge-Kutta 2 (Heun) along $-\nabla T / |\nabla T|$ with bilinear gradient interpolation.
- **Geophysical Conventions**: $z$ positive downwards (depth), C-order arrays (NumPy standard), free and unforced physical units.
- **Zero-Copy Python Bindings**: Direct NumPy array sharing using PyO3 and Maturin.

---

## Installation

### From PyPI (pre-built wheels)

```bash
pip install seiscode-eikonal
```

### From Source (requires Rust >= 1.85)

```bash
git clone https://github.com/pageotd/seiscode-eikonal.git
cd seiscode-eikonal
pip install -e .
```

---

## Quickstart

### 1. Computing 2D Traveltimes (`solve_2d`)

```python
import numpy as np
import seiscode_eikonal as eikonal

# Define grid: 201 nodes in X (width: 2000 m), 101 nodes in Z (depth: 1000 m)
nx, nz = 201, 101
dx, dz = 10.0, 10.0
x0, z0 = 0.0, 0.0  # Top-left grid origin

# Velocity model: 2000 m/s with a vertical gradient
velocities = np.empty((nx, nz), dtype=np.float64)
for iz in range(nz):
    velocities[:, iz] = 2000.0 + iz * 15.0

# Source position (X = 100 m, Z = 50 m depth)
xs, zs = 100.0, 50.0

# Solve Eikonal equation with 2nd-order stencil
traveltimes = eikonal.solve_2d(
    velocities=velocities,
    dx=dx,
    dz=dz,
    xs=xs,
    zs=zs,
    x0=x0,
    z0=z0,
    order=2,
)

print(f"Max traveltime: {traveltimes.max():.3f} s")
```

### 2. Tracing Seismic Rays (`trace_rays`)

```python
# Define receivers along surface and borehole
receivers = np.array([
    [500.0, 0.0],
    [1000.0, 0.0],
    [1500.0, 0.0],
    [1500.0, 500.0],
], dtype=np.float64)

# Trace rays by backtracking along -grad(T)
rays = eikonal.trace_rays(
    traveltimes=traveltimes,
    dx=dx,
    dz=dz,
    xs=xs,
    zs=zs,
    receivers=receivers,
    x0=x0,
    z0=z0,
)

for idx, ray in enumerate(rays):
    print(f"Ray {idx}: {len(ray)} points from receiver {receivers[idx]} to source")
```

### 3. Visualizing Results

An end-to-end example displaying velocity, traveltime contours (isochrones), and raypaths is available in [`examples/plot_rays_two_layer.py`](examples/plot_rays_two_layer.py):

```bash
python examples/plot_rays_two_layer.py
```

![Seismic Ray Tracing Demo](ray_tracing_demo.png)

---

## Conventions & Best Practices

1. **Axes and Coordinates**:
   - $X$ is horizontal distance (index `ix`, range $0 \le ix < nx$).
   - $Z$ is vertical depth, **positive downwards** (index `iz`, range $0 \le iz < nz$).
   - Coordinates map to grid indices as: $x = x_0 + ix \cdot dx$, $z = z_0 + iz \cdot dz$.
2. **Memory Layout**:
   - Arrays passed to `solve_2d` must be 2D C-contiguous arrays of shape `(nx, nz)` with `dtype=np.float64`.
   - If Fortran-ordered arrays are provided, an automatic copy is performed.
3. **Physical Units**:
   - Units are intentionally arbitrary but must be mutually consistent.
   - Example 1: $dx, dz$ in meters ($\text{m}$), $v$ in meters per second ($\text{m/s}$), yields traveltime in seconds ($\text{s}$).
   - Example 2: $dx, dz$ in kilometers ($\text{km}$), $v$ in kilometers per second ($\text{km/s}$), yields traveltime in seconds ($\text{s}$).

---

## Architecture

The project is structured as a two-crate Cargo workspace:

- **`crates/seiscode-eikonal-core`**:
  - Independent, pure Rust library.
  - Implements the grid structure, local upwind stencils, min-heap priority queue with lazy deletion, FMM solver, and Runge-Kutta 2 ray backtracking.
  - Zero PyO3 or Python dependencies; can be tested and benchmarked natively.
- **`crates/seiscode-eikonal-py`**:
  - Thin binding layer using PyO3 and Maturin.
  - Handles NumPy pointer extraction, shape validation, and Rust-to-Python error mapping.

---

## Documentation

Comprehensive documentation is available in the [`docs/`](docs/) directory:

- [**User Guide & Python API**](docs/user-guide.md): Complete parameter descriptions, error handling, and receiver ray tracing.
- [**Numerical & Physical Theory**](docs/theory.md): Eikonal equation formulation, upwind stencil derivations, singularity correction, and ray equations.
- [**Architecture & Development**](docs/architecture.md): Workspace layout, Rust safety invariants, testing strategy, and benchmarking.

---

## Contributing & Testing

Local tasks are orchestrated using [Taskfile](https://taskfile.dev):

```bash
task lint       # Run cargo clippy and fmt checks
task test       # Run all tests (Rust unit + integration + Python pytest)
task audit      # Run cargo-deny and cargo-audit security checks
task ci         # Run the full suite expected by CI
```

---

## License

Distributed under the MIT License. See [LICENSE](LICENSE) for details.
