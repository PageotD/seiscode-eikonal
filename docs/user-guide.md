# User Guide & Python API Reference

`seiscode-eikonal` provides two core capabilities:
1. Solving the 2D Eikonal equation for first-arrival traveltimes via the Fast Marching Method (`solve_2d`).
2. Tracing seismic rays from receivers back to the source point via gradient backtracking (`trace_ray` and `trace_rays`).

---

## 1. Conventions & Coordinate System

### Coordinate Axes
- **$X$-axis**: Horizontal coordinate. Corresponds to the first dimension of the velocity array (length $nx$).
- **$Z$-axis**: Vertical coordinate, **positive downwards** (geophysical depth convention). Corresponds to the second dimension (length $nz$).
- **Physical origin**: `(x0, z0)` corresponds to node index `(0, 0)`.
- **Node coordinates**:
  $$x_i = x_0 + i_x \cdot dx, \quad z_k = z_0 + i_z \cdot dz$$

### Data Types and Memory Layout
- **Precision**: Calculations are performed in double precision (`f64` in Rust, `np.float64` in Python).
- **Memory ordering**: Arrays must be C-contiguous.
  - Shape: `(nx, nz)`.
  - In NumPy: `array.shape == (nx, nz)` and `array.flags.c_contiguous == True`.
  - If a Fortran-ordered array or non-contiguous slice is passed, `seiscode-eikonal` automatically copies it into C-order and logs a warning.

### Units
No physical units are hardcoded. Ensure your units are mutually consistent:
| Length ($dx, dz, x, z$) | Velocity ($v$) | Slowness ($s = 1/v$) | Traveltime ($T$) |
|---|---|---|---|
| meters ($\text{m}$) | $\text{m/s}$ | $\text{s/m}$ | seconds ($\text{s}$) |
| kilometers ($\text{km}$) | $\text{km/s}$ | $\text{s/km}$ | seconds ($\text{s}$) |

---

## 2. API Reference

### `solve_2d`

Solves the 2D Eikonal equation on a Cartesian grid for a single source.

```python
import seiscode_eikonal as eikonal

traveltimes = eikonal.solve_2d(
    velocities: np.ndarray,
    dx: float,
    dz: float,
    xs: float,
    zs: float,
    x0: float = 0.0,
    z0: float = 0.0,
    order: int = 2,
) -> np.ndarray
```

#### Parameters:
- **`velocities`** (*`np.ndarray`, shape `(nx, nz)`, dtype `np.float64`*):
  The 2D velocity model. All values must be strictly positive and finite.
- **`dx`** (*`float`*):
  Grid spacing along the horizontal $X$-axis ($dx > 0$).
- **`dz`** (*`float`*):
  Grid spacing along the vertical $Z$-axis ($dz > 0$).
- **`xs`** (*`float`*):
  Horizontal coordinate of the source point. Must be within $[x_0, x_0 + (nx - 1) \cdot dx]$.
- **`zs`** (*`float`*):
  Vertical depth coordinate of the source point. Must be within $[z_0, z_0 + (nz - 1) \cdot dz]$.
- **`x0`** (*`float`*, optional, default `0.0`):
  Horizontal coordinate of the grid origin (node `[0, 0]`).
- **`z0`** (*`float`*, optional, default `0.0`):
  Vertical coordinate of the grid origin (node `[0, 0]`).
- **`order`** (*`int`*, optional, default `2`):
  Upwind finite-difference order:
  - `1`: First-order upwind scheme (standard FMM).
  - `2`: Second-order upwind scheme (higher accuracy, reduced diagonal grid error). Falls back to 1st-order near boundaries or when 2nd-order upstream nodes are unavailable.

#### Returns:
- **`traveltimes`** (*`np.ndarray`, shape `(nx, nz)`, dtype `np.float64`*):
  2D array of first-arrival traveltimes across the grid.

---

### `trace_ray`

Traces a single seismic ray from a receiver `(xr, zr)` back to the source `(xs, zs)`.

```python
ray = eikonal.trace_ray(
    traveltimes: np.ndarray,
    dx: float,
    dz: float,
    xs: float,
    zs: float,
    xr: float,
    zr: float,
    x0: float = 0.0,
    z0: float = 0.0,
    step_size: float | None = None,
) -> np.ndarray
```

#### Parameters:
- **`traveltimes`** (*`np.ndarray`, shape `(nx, nz)`, dtype `np.float64`*):
  Traveltime field computed by `solve_2d`.
- **`dx`**, **`dz`** (*`float`*):
  Grid spacings along $X$ and $Z$.
- **`xs`**, **`zs`** (*`float`*):
  Source coordinates.
- **`xr`**, **`zr`** (*`float`*):
  Receiver coordinates. Must be within grid bounds.
- **`x0`**, **`z0`** (*`float`*, optional, default `0.0`):
  Grid origin coordinates.
- **`step_size`** (*`float`*, optional, default `None`):
  Integration step length $\Delta s$. If `None` or $\le 0$, defaults to $0.5 \cdot \min(dx, dz)$.

#### Returns:
- **`ray`** (*`np.ndarray`, shape `(N, 2)`, dtype `np.float64`*):
  Discrete trajectory coordinates `[[x_0, z_0], [x_1, z_1], ..., [x_s, z_s]]`.
  The first point `ray[0]` is the receiver, and the last point `ray[-1]` is the exact source position.

---

### `trace_rays`

Traces multiple seismic rays from an array of receivers back to the source.

```python
rays = eikonal.trace_rays(
    traveltimes: np.ndarray,
    dx: float,
    dz: float,
    xs: float,
    zs: float,
    receivers: np.ndarray,
    x0: float = 0.0,
    z0: float = 0.0,
    step_size: float | None = None,
) -> list[np.ndarray]
```

#### Parameters:
- **`receivers`** (*`np.ndarray`, shape `(M, 2)`, dtype `np.float64`*):
  2D array where each row `[xr, zr]` specifies the coordinates of a receiver.
- Other parameters are identical to `trace_ray`.

#### Returns:
- **`rays`** (*`list[np.ndarray]`*):
  A list of $M$ arrays, where each element is an array of shape `(N_i, 2)` containing the trajectory of the $i$-th ray.

---

## 3. Error Handling

`seiscode-eikonal` validates all inputs in the Rust core and maps errors to clear Python exceptions:

```python
import seiscode_eikonal as eikonal
import numpy as np

v = np.ones((50, 50), dtype=np.float64)

# 1. Invalid Grid Spacing
try:
    eikonal.solve_2d(v, dx=-1.0, dz=10.0, xs=250.0, zs=250.0)
except ValueError as err:
    print(f"Caught expected error: {err}")
    # Output: Invalid grid spacing: dx=-1, dz=10 (must be strictly positive)

# 2. Source Outside Bounds
try:
    eikonal.solve_2d(v, dx=10.0, dz=10.0, xs=9999.0, zs=250.0)
except ValueError as err:
    print(f"Caught expected error: {err}")
    # Output: Source position (9999, 250) is out of bounds
```
