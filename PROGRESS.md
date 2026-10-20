# Initial Phase Progress

Progress tracker for the initial phase of the Euler solver implementation.

- [x] `read_settings.f90`
- [x] `read_geom.f90`
- [x] `generate_mesh.f90`
- [x] `calc_areas.f90`
- [x] `check_mesh.f90`
- [X] `flow_guess.f90`
- [x] `set_timestep.f90`
- [X] `set_secondary.f90`
- [X] `apply_bconds.f90`
- [X] `euler_iteration.f90`
- [ ] `flux_stencil.f90`
- [ ] `smooth_stencil.f90`
- [ ] `routines.f90`
- [ ] `plot_contours.py`

## Status

Completed: 3 / 14

## Completion evidence

| Script | Commit | Message |
| --- | --- | --- |
| `read_settings.f90` | `8a4c134` | finish read_setting, capable of read the given input now |
| `read_geom.f90` | `b2da220` | finish read_geom |
| `generate_mesh.f90` | `b4725af` | finish generate_mesh.f90 |

## User-defined types

These types are declared in `types.f90` and passed to nearly every subroutine.

### `av` — `t_appvars`

Application variables and gas constants; holds run settings rather than mesh data.

| Variable | Type | Description |
| --- | --- | --- |
| `casename` | `character(len=:)` | Name of the case, used to build output filenames |
| `rgas`, `gam`, `cp`, `cv`, `fgam` | `real` | Gas properties: gas constant, ratio of specific heats $\gamma$, specific heats, and $\gamma - 1$ |
| `cfl` | `real` | CFL number controlling the timestep size |
| `sfac` | `real` | Smoothing factor for the artificial smoothing stencil |
| `dt` | `real` | Current timestep |
| `d_max`, `d_avg` | `real` | Convergence criteria (maximum and average normalised residuals) |
| `nsteps` | `integer` | Total number of timesteps requested |
| `nstep` | `integer` | Current timestep number |
| `ro_ref`, `roe_ref`, `rov_ref` | `real` | Reference values of the primary flow variables, used to normalise residuals |
| `ni`, `nj` | `integer` | Mesh size, set for single-block test cases |
| `nn`, `nm` | `integer` | Number of blocks and matching patches for the multi-block extension |

### `g` — `t_grid`

Grid geometry plus all flow variables stored on the mesh.

| Variable | Type | Description |
| --- | --- | --- |
| `ni`, `nj` | `integer` | Mesh dimensions (copy of `av%ni`, `av%nj`) |
| `x`, `y` | `real(:,:)` | Node coordinates, dimension `(ni,nj)` |
| `area` | `real(:,:)` | Cell areas, dimension `(ni-1,nj-1)` |
| `lx_i`, `ly_i` | `real(:,:)` | $x$ and $y$ components of the projected `i = const` facet lengths, dimension `(ni,nj-1)` |
| `lx_j`, `ly_j` | `real(:,:)` | $x$ and $y$ components of the projected `j = const` facet lengths, dimension `(ni-1,nj)` |
| `l_min` | `real` | Overall minimum element size, used to set geometric tolerances |
| `ro`, `roe`, `rovx`, `rovy` | `real(:,:)` | Primary conserved variables at the nodes: density, total energy, $x$- and $y$-momentum |
| `dro`, `droe`, `drovx`, `drovy` | `real(:,:)` | Cell-centred increments (updates/residuals) of the primary variables |
| `p`, `hstag`, `vx`, `vy` | `real(:,:)` | Secondary variables at the nodes: pressure, stagnation enthalpy, $x$- and $y$-velocity |
| `wall` | `logical(:,:)` | Flags nodes lying on a solid wall boundary |

### `bcs` — `t_bconds`

Boundary condition data for the inlet and outlet. Populated by `read_settings.f90` and applied by `apply_bconds.f90`.

| Variable | Type | Description |
| --- | --- | --- |
| `pstag`, `tstag` | `real` | Stagnation pressure and temperature at the inlet |
| `alpha` | `real` | Inlet flow angle, converted from degrees to radians in `read_settings` |
| `rfin` | `real` | Relaxation factor for the inlet density update (typically 0.25) |
| `rostag` | `real` | Inlet stagnation density, computed as `pstag / (rgas * tstag)` |
| `p`, `ro` | `real(:)` | Static pressure and density along the inlet plane, dimension `(nj)`; `ro` is relaxed each timestep in `apply_bconds` and limited to just below `rostag` for stability |
| `p_out` | `real` | Outlet static pressure |
| `n_in`, `n_out` | `integer` | Block numbers of the inlet and outlet, used in the multi-block extension |
