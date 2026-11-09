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
- [X] `flux_stencil.f90`
- [X] `smooth_stencil.f90`
- [X] `routines.f90`
- [X] `plot_contours.py`

## Status

Completed: 14 / 14

## Solver configuration (`read_settings`)

- [x] Select subsonic or supersonic boundary conditions from input.
- [x] Select crude/uniform or refined initialization from input.
- [x] Configure deferred correction strength, including disabling it.
- [x] Generate and read feature settings in Python and Fortran.

`Code/generate_case.py` uses `default_settings` and `write_settings` in
`Code/routines.py`. Change the returned settings dictionary in a case generator,
or edit the generated input file directly. The input records are:

```text
casename
gam rgas
cfl sfac d_max
nsteps
ni nj
pstag tstag alpha rfin
p
supersonic guesstype fcorr
```

The last line is optional. Omitting it preserves current `develop` defaults:
`1 2 0.9`. Generated files include it explicitly. Both readers reject malformed
feature settings and values outside the supported ranges.

| Setting | Allowed values | Default | Effect |
| --- | --- | --- | --- |
| `supersonic` | `0`, `1` | `1` | `0`: relaxed subsonic inlet and prescribed outlet pressure from `master` (`73f3384`). `1`: prescribed supersonic inlet and no imposed outlet pressure. |
| `guesstype` | `1`, `2` | `2` | `1`: crude subsonic guess or uniform supersonic inlet state. `2`: refined subsonic guess or supersonic mass-flow/Mach bisection guess. |
| `fcorr` | finite real in `[0, 1]` | `0.9` | Deferred correction strength. `0` disables its effect and clears correction history; `0.9` retains previous behavior. History relaxation remains `0.99 * old + 0.01 * target`. |

`p` is **outlet static pressure** with `supersonic = 0` and **inlet static
pressure** with `supersonic = 1`; use pressure appropriate to the chosen mode.
`alpha` is read in degrees and converted to radians. `rfin` controls inlet
density relaxation only in subsonic mode. Fortran scales `sfac` and `d_max` by
`cfl`, then sets `d_avg = 0.5 * d_max`. Selected feature settings and the pressure
interpretation are printed in the solver log.

Examples of the final line:

```text
1 2 0.9
```

retains current `develop` behavior; `0 2 0` selects restored subsonic behavior
with refined initialization and deferred correction off; `1 1 0` selects a
uniform supersonic guess with deferred correction off.

Flux-stencil fixes, the second-order flux module, and plotting improvements
remain enabled. Fourth-order flux is an unimplemented stub; fourth-order
smoothing, mesh clustering, local timesteps, and multi-block repairs are outside
this configuration change.

### Configuration verification

Validated with `gfortran -fcheck=all -g` on a nonuniform 25-by-9 duct mesh
for 50 steps:

- All eight combinations of boundary mode, guess type and correction off/on
  produced finite flow fields; Python settings round trips passed.
- Legacy input and explicit `1 2 0.9` produced identical final output.
- Default initialization and final output matched original `develop` exactly;
  restored subsonic refined initialization with correction off matched `master`
  exactly for this case.
- Both readers rejected ten malformed/out-of-range records, including nonfinite
  correction strengths. Correction on/off changed the flow in both modes.

These are configuration and regression checks, not full convergence studies.

## Planned improvements for Extension B

Keep each improvement switchable so its effect can be compared with the basic solver and with the other improvements.

- [ ] Fourth-order smoothing: reduce artificial diffusion, retaining second-order smoothing where needed for stability near shocks.
- [ ] Refined mesh density: compare uniform spacing with local mesh clustering to improve accuracy for a given node count. Account for non-uniform spacing in the numerical stencils.
- [ ] Spatially varying timestep: use local timesteps to accelerate convergence towards steady flow. This is a convergence-acceleration method and does not preserve a common physical time across the domain.

Compare accuracy, stability, convergence speed and CPU time with each improvement enabled separately and in combination.

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
| `supersonic`, `guesstype` | `integer` | Boundary mode and initialization selection; see configuration above |
| `fcorr` | `real` | Deferred correction strength (`0` disables correction) |
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
| `corr_ro`, `corr_roe`, `corr_rovx`, `corr_rovy` | `real(:,:)` | Deferred correction histories, initialized to zero |

### `bcs` — `t_bconds`

Boundary condition data for the inlet and outlet. Populated by `read_settings.f90` and applied by `apply_bconds.f90`.

| Variable | Type | Description |
| --- | --- | --- |
| `pstag`, `tstag` | `real` | Stagnation pressure and temperature at the inlet |
| `alpha` | `real` | Inlet flow angle, converted from degrees to radians in `read_settings` |
| `rfin` | `real` | Relaxation factor for the inlet density update (typically 0.25) |
| `rostag` | `real` | Inlet stagnation density, computed as `pstag / (rgas * tstag)` |
| `p` | `real` | Prescribed static pressure: outlet in subsonic mode, inlet in supersonic mode |
| `ro` | `real(:)` | Inlet density, dimension `(nj)`; relaxed and limited below `rostag` in subsonic mode, prescribed in supersonic mode |
| `n_in`, `n_out` | `integer` | Block numbers of the inlet and outlet, used in the multi-block extension |
