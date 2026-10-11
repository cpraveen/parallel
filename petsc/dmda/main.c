#include <petscdmda.h>

int main(int argc, char **argv) {
    DM           da;
    PetscInt     M = 20, N = 50; /* 10x10 global grid */
    DM           cda;
    Vec          coords;
    DMDACoor2d   **carray;
    PetscInt     xs, ys, xm, ym, i, j;
    PetscViewer  viewer;
    Vec          x;

    /* Initialize PETSc */
    PetscCall(PetscInitialize(&argc, &argv, NULL, "2D DMDA Coordinate Example\n"));

    /* ====================================================================
       1. CREATE THE 2D DMDA
       ==================================================================== */
    PetscCall(DMDACreate2d(PETSC_COMM_WORLD,
                        DM_BOUNDARY_NONE, DM_BOUNDARY_NONE, /* Boundary conditions */
                        DMDA_STENCIL_STAR,                  /* Stencil type */
                        M, N,                               /* Global grid dimensions */
                        PETSC_DECIDE, PETSC_DECIDE,         /* Procs per dimension */
                        1,                                  /* Degrees of freedom per node */
                        1,                                  /* Stencil width */
                        NULL, NULL,                         /* Processor ownership arrays */
                        &da));

    PetscCall(DMSetUp(da)); /* Always call DMSetUp after creation */

    /* ====================================================================
       2. SET UNIFORM COORDINATES
       ==================================================================== */
    /* Set bounds: x from 0.0 to 1.0, y from 0.0 to 1.0.
       Note: zmin/zmax must be provided even for 2D, usually as 0.0 */
    PetscCall(DMDASetUniformCoordinates(da, 0.0, 1.0, 0.0, 1.0, 0.0, 0.0));

    /* ====================================================================
       3. MODIFY COORDINATES MANUALLY (For Non-Uniform/Distorted Grids)
       ==================================================================== */

    /* Get the coordinate DM and the coordinate vector */
    PetscCall(DMGetCoordinateDM(da, &cda));
    PetscCall(DMGetCoordinates(da, &coords));

    /* Get local array bounds for this processor */
    PetscCall(DMDAGetCorners(da, &xs, &ys, NULL, &xm, &ym, NULL));

    /* Extract 2D array of coordinates.
       Crucial: Use the coordinate DM (cda), not the base DM (da) here! */
    PetscCall(DMDAVecGetArray(cda, coords, &carray));

    for (j = ys; j < ys + ym; j++)
        for (i = xs; i < xs + xm; i++)
        {
           double r = carray[j][i].x + 1.0;
           double theta = 2.0 * M_PI * carray[j][i].y;

            /* Example: Perturb the x coordinates slightly based on y */
            carray[j][i].x = r*cos(theta);
            carray[j][i].y = r*sin(theta);

        }

    /* Restore the array to the vector */
    PetscCall(DMDAVecRestoreArray(cda, coords, &carray));

    /* ====================================================================
       4. OUTPUT TO VTK FOR VISUALIZATION
       ==================================================================== */

    /* Create a vector (e.g. a solution field) defined on the DMDA grid. */
    PetscCall(DMCreateGlobalVector(da, &x));
    PetscCall(VecSet(x, 1.0));

    /* Open VTK viewer. Using .vts extension tells PETSc to write VTK Struct Grid */
    PetscCall(PetscViewerVTKOpen(PETSC_COMM_WORLD, "grid.vts", FILE_MODE_WRITE, &viewer));

    /* View the vector. This writes the grid topology and the (distorted)
       coordinates, along with the field data, to grid.vts. */
    PetscCall(VecView(x, viewer));

    /* Destroy the viewer */
    PetscCall(PetscViewerDestroy(&viewer));
    PetscCall(VecDestroy(&x));

    /* Cleanup */
    PetscCall(DMDestroy(&da));
    PetscCall(PetscFinalize());
    return 0;
}