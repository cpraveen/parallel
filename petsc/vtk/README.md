# DMDA → VTK output

Writes the exact 2D isentropic vortex solution of the compressible Euler
equations from a PETSc DMDA grid to VTK files: one `.vts` per time step plus
a `.pvd` collection for time-series visualization in ParaView/VisIt.

## What the code does

- Builds a 64×64 DMDA on the domain [-5,5]² (grid size overridable via
  `-da_grid_x` / `-da_grid_y`).
- Solves nothing: it evaluates the exact isentropic vortex
  (freestream (ρ,u,v,p) = (1,1,0,1), γ = 1.4, β = 5) and advects it with the
  freestream velocity u = 1.
- 4 conserved variables (ρ, ρu, ρv, E) are stored in one vector of dof = 4;
  field names appear in the VTK output as `U.Density`, `U.MomentumX`,
  `U.MomentumY`, `U.Energy` (the vector is named `U`).
- Time loop from t = 0 to tfinal = 2.0 with dt = 0.05, saving 41 frames as
  `solution_0000.vts … solution_0040.vts` and appending each to `solution.pvd`.

## Build

Requires PETSc (developed against 3.25.5, installed via Spack at `/opt/spack`).

`PETSC_DIR` must be set to your PETSc installation directory — either as an
environment variable or on the make command line:

    export PETSC_DIR=/opt/spack
    make

or:

    make PETSC_DIR=/opt/spack

The Makefile aborts with an error if `PETSC_DIR` is not set.

## Run

    ./main

Runtime options:

| Option       | Default | Meaning                |
|--------------|---------|------------------------|
| `-dt`        | 0.05    | Time step size         |
| `-tfinal`    | 2.0     | Final time             |
| `-da_grid_x` | 64      | Grid points in x (also y) |
| `-da_grid_y` | 64      | Grid points in y       |

## Output and visualization

- `solution_XXXX.vts` — one structured-grid file per frame (appended raw binary).
- `solution.pvd` — ParaView collection referencing all frames with timesteps.

Open `solution.pvd` in ParaView (or Visit) to animate the advecting vortex.

## Clean up

    make clean

removes the `main` binary, all `solution_*.vts` files, and `solution.pvd`.

## Sanity check

The vortex core (minimum density, ρ ≈ 0.50) starts at the origin and moves
along the x-axis: at t = 2 its center is near x = 2.