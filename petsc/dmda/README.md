# PETSc 2D DMDA Coordinate Example

This program demonstrates how to create a 2D structured grid with **PETSc's DMDA**,
modify its node coordinates by hand, and write the result to a VTK file
(`grid.vts`) for visualization in ParaView.

The grid is a 20×50 (M×N) logical mesh defined on the unit square [0,1]×[0,1].
Each node's coordinates are then remapped from the square onto an annulus
(ring) using a polar transformation, so the output shows a distorted,
non-uniform grid.

## What the code does

The program (`main.c`) works in four steps:

1. **Create the DMDA**
   - `DMDACreate2d` builds a 20×50 global grid with
     - no periodic boundaries (`DM_BOUNDARY_NONE`) in x and y,
     - a star stencil (`DMDA_STENCIL_STAR`) of width 1,
     - one degree of freedom per node.
   - Processor partitioning is left to PETSc (`PETSC_DECIDE`), so the code
     runs with any number of MPI ranks.
   - `DMSetUp` finalizes the DM.

2. **Set uniform coordinates**
   - `DMDASetUniformCoordinates` assigns bounds x ∈ [0,1], y ∈ [0,1]
     (z bounds are required by the API even in 2D; they are set to 0).

3. **Modify the coordinates manually**
   - `DMGetCoordinateDM` / `DMGetCoordinates` retrieve the coordinate
     DM (`cda`) and the coordinates vector.
   - `DMDAVecGetArray` exposes the coordinates as a 2D array of
     `DMDACoor2d` structs, local to the current rank.
   - Each (x, y) point is mapped to polar coordinates
     - `r     = x + 1`        → radius in [1, 2],
     - `theta = 2π · y`       → angle in [0, 2π),
     - x → r·cos(θ), y → r·sin(θ).
   - The result is an annular grid: 50 nodes around the ring, 20 nodes
     in the radial direction.
   - `DMDAVecRestoreArray` gives the array back to PETSc.

4. **Write VTK output**
   - `DMCreateGlobalVector` creates a field vector on the grid
     (here set to the constant 1.0 as a placeholder for a solution field).
   - `PetscViewerVTKOpen` opens a VTK viewer; the `.vts` extension selects
     the VTK Structured Grid format.
   - `VecView` writes the grid topology, the (distorted) coordinates, and
     the field data to `grid.vts`.

## Build and run

Requirements:

- A PETSc installation with the DMDA component (standard).
- An MPI compiler (the `Makefile` uses PETSc's `mpicc` rules).
- `PETSC_DIR` must point at the PETSc installation.

```sh
export PETSC_DIR=/path/to/petsc    # or the equivalent spack module
make
./main
```

The run produces `grid.vts` in the current directory.

> **Running in parallel:** with `mpirun -np N ./main`, PETSc writes one
> `.vts` per rank plus a `.visit` master file, instead of a single
> `grid.vts`. That is expected distributed behavior — point ParaView at the
> `.visit` file in that case.

## Output

`grid.vts` is a VTK `StructuredGrid`:

```
WholeExtent="0 19 0 49 0 0"
```

It contains:

- **Points** — the distorted (annular) node coordinates, and
- **PointData** — the scalar field `"Unnamed"` (constant 1.0 in this example).

Open it with **ParaView**: *File → Open → grid.vts*. Use the "Unstructured/
Structured grid" representation and toggle visibility of the field or of the
mesh surface to inspect the annular grid.

## Key implementation notes / gotchas

- **A VTK file is only written when a `Vec` is viewed.** The `.vts` output
  comes from `VecView(x, viewer)`; the viewer opens the file lazily, so the
  grid is written here only because of the `VecView` call.
- **Access coordinates through the coordinate DM.** When editing node
  coordinates you must extract the array with the coordinate DM (`cda`),
  *not* the base DM (`da`).
- **`DMDAVecGetArray` / `DMDAVecRestoreArray` are a matched pair.** Always
  restore the array before any further PETSc operations that touch the vector.
- **Error handling** uses the modern `PetscCall` macro (which replaced the
  deprecated `CHKERRQ`). Every PETSc call is checked automatically.

## Files

| File       | Purpose                                             |
|------------|-----------------------------------------------------|
| `main.c`   | Source code                                         |
| `Makefile` | Build rules (requires `PETSC_DIR`)                  |
| `main`     | Compiled executable (`make`)                        |
| `grid.vts` | VTK output produced by `./main`                     |