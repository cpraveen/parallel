/*
 * Solves 2D Poisson equation -∆u = f = rhs            in Ω = [0,1]x[0,1]
 *                              u = g = boundary_value in ∂Ω
 *             using Jacobi iterations.
 *
 *  Date: Friday, September 8th, 2023.
 */
#include <iostream>
#include <fstream>
#include <petsc.h>
#include <petscdm.h>
#include <petscsys.h>
#include <petscvec.h>
#include <petscviewerhdf5.h>
#include <stdio.h>

using namespace std;

//------------------------------------------------------------------------------
// Application context to store some data
struct AppCtx
{
  PetscInt uindices[2][2];
  PetscReal dx, dy;
  // DMDACoor2d is just a struct of 2 PetscScalars x and y.
  DMDACoor2d **coord;
};

PetscErrorCode writeVTK(AppCtx *ctx, DM da, Vec u_global, PetscInt c);

//------------------------------------------------------------------------------
PetscReal rhs(PetscInt x, PetscInt y) 
{ 
  return 4.0 * (pow(x, 4) + pow(y, 4)); 
}

//------------------------------------------------------------------------------
PetscReal boundary_value(PetscReal x, PetscReal y) 
{ 
  return x * x + y * y; 
}

//------------------------------------------------------------------------------
PetscErrorCode set_initial_guess(AppCtx *ctx, DM da, Vec u_global)
{
  // Fill boundary values
  PetscInt nx, ny;
  PetscCall(DMDAGetInfo(da, NULL, &nx, &ny, NULL, NULL, NULL, NULL, NULL, NULL,
                        NULL, NULL, NULL, NULL));

  PetscInt ibeg, jbeg, nlocx, nlocy;
  PetscCall(DMDAGetCorners(da, &ibeg, &jbeg, NULL, &nlocx, &nlocy, NULL));

  PetscScalar **u;
  PetscCall(DMDAVecGetArray(da, u_global, &u));

  for (PetscInt j = jbeg; j < jbeg + nlocy; ++j)
    for (PetscInt i = ibeg; i < ibeg + nlocx; ++i)
    {
      if (i == 0 || i == nx - 1 || j == 0 || j == ny - 1)
      {
        u[j][i] = boundary_value(ctx->coord[j][i].x, ctx->coord[j][i].y);
      }
      else
        u[j][i] = 0.0;
    }

  PetscCall(DMDAVecRestoreArray(da, u_global, &u));
  PetscFunctionReturn(PETSC_SUCCESS);
}

//------------------------------------------------------------------------------
PetscErrorCode apply_jacobi(AppCtx *ctx, DM da, Vec u_local,
                            Vec u_global, double &maxdelta)
{
  // Make u_old point to the data in the u_local array
  PetscScalar **u_old;
  PetscCall(DMDAVecGetArrayRead(da, u_local, &u_old));

  // Get a pointer to the portion of memory where we will be putting the new
  // solution.
  PetscScalar **u_new;
  PetscCall(DMDAVecGetArrayWrite(da, u_global, &u_new));

  const PetscReal one_over_dx_squared = 1.0 / (ctx->dx * ctx->dx);
  const PetscReal one_over_dy_squared = 1.0 / (ctx->dy * ctx->dy);
  const PetscReal coefficient = 0.5 / (one_over_dx_squared + one_over_dy_squared);

  // Loop over interior points, dont update boundary points where dirichlet
  // bc is specified.
  maxdelta = 0.0;
  for (PetscInt j = ctx->uindices[1][0]; j < ctx->uindices[1][1]; ++j)
    for (PetscInt i = ctx->uindices[0][0]; i < ctx->uindices[0][1]; ++i)
    {
      // Because u_local is a local array, ghost points will be accessible.
      u_new[j][i] =
          coefficient
          * (rhs(ctx->coord[j][i].x, ctx->coord[j][i].y)
             + (u_old[j - 1][i] + u_old[j + 1][i]) * one_over_dy_squared
             + (u_old[j][i - 1] + u_old[j][i + 1]) * one_over_dx_squared);
      maxdelta = std::max(abs(u_new[j][i] - u_old[j][i]), maxdelta);
    }

  PetscCall(DMDAVecRestoreArrayRead(da, u_local, &u_old));
  PetscCall(DMDAVecRestoreArrayWrite(da, u_global, &u_new));
  PetscFunctionReturn(PETSC_SUCCESS);
}

//------------------------------------------------------------------------------
int main(int argc, char *argv[])
{
  char      help[]  = "Solves -∆u = f\n";
  PetscReal xmin    = -1.0;
  PetscReal xmax    = 1.0;
  PetscReal ymin    = -1.0;
  PetscReal ymax    = 1.0;
  PetscInt  nx      = 100;
  PetscInt  ny      = 100;
  PetscReal eps     = 1e-5;
  PetscInt  itermax = 10000;
  // DM stands for "distributed mesh"
  // da stands for "distributed array"
  DM        da;
  Vec       u_local, u_global;
  AppCtx    ctx;

  // This calls MPI_Init unless we have already.
  PetscCall(PetscInitialize(&argc, &argv, NULL, help));

  // PetscOptionsBegin and ...End are macros. They don't have return
  // values, so they can't be used with PetscCall.
  PetscOptionsBegin(PETSC_COMM_WORLD, "jacobi_",
                    "options for number of grid of points", "");
  // TODO Put defaults here, and declare them as uninit'd consts above.
  PetscCall(PetscOptionsInt("-itermax", "Maximum number of iterations",
                            "main.cc", itermax, &itermax, NULL));
  PetscCall(PetscOptionsReal("-eps", "Step size to halt at", "main.cc", eps,
                             &eps, NULL));
  PetscOptionsEnd();

  // Create a 2D Distributed Mesh of DA type.
  // We choose BOUNDARY_NONE because we don't want a layer of ghosts
  // around the whole domain. We know the values on the boundary.
  PetscCall(DMDACreate2d(PETSC_COMM_WORLD,
                         DM_BOUNDARY_NONE,
                         DM_BOUNDARY_NONE,
                         DMDA_STENCIL_STAR,
                         nx, ny,
                         PETSC_DECIDE, PETSC_DECIDE,
                         1, 1,
                         NULL, NULL,
                         &da));
  // Sets parameters from the options database
  PetscCall(DMSetFromOptions(da));
  // Now actual set_up_. That was just set params from options.
  PetscCall(DMSetUp(da));

  // Sets uniform coordinates for the grid. Z-values are ignored in 2D
  PetscCall(DMDASetUniformCoordinates(da, xmin, xmax, ymin, ymax, 0.0, 0.0));

  // Find h
  PetscCall(DMDAGetInfo(da, NULL, &nx, &ny, NULL, NULL, NULL, NULL, NULL, NULL,
                        NULL, NULL, NULL, NULL));

  ctx.dx = (xmax - xmin) / (PetscReal)(nx);
  ctx.dy = (ymax - ymin) / (PetscReal)(ny);

  PetscCall(DMCreateLocalVector(da, &u_local));
  PetscCall(DMCreateGlobalVector(da, &u_global));
  PetscCall(PetscObjectSetName((PetscObject)u_global, "Poisson solution"));
  PetscCall(PetscObjectSetName((PetscObject)u_local, "Local solution"));

  // Create an array with local coordinates.
  DM           cda;
  Vec          clocal;
  PetscCall(DMGetCoordinateDM(da, &cda));
  // This one has da, not cda
  PetscCall(DMGetCoordinatesLocal(da, &clocal));
  PetscCall(DMDAVecGetArrayRead(cda, clocal, &ctx.coord));

  PetscCall(set_initial_guess(&ctx, da, u_global));

  // Save initial condition.
  writeVTK(&ctx, da, u_global, 0);

  /* ----------------------------
   * END OF SETUP
   * BEGIN ITERATIONS
   * ----------------------------
   */
  // Find the indices to update with Jacobi iterations
  // If the first or last point is on the global boundary then don't update it.
  PetscInt ibeg, jbeg, nlocx, nlocy;
  PetscCall(DMDAGetCorners(da, &ibeg, &jbeg, NULL, &nlocx, &nlocy, NULL));

  ctx.uindices[0][0] = (ibeg == 0) ? 1 : ibeg;
  ctx.uindices[1][0] = (jbeg == 0) ? 1 : jbeg;
  ctx.uindices[0][1] = (ibeg + nlocx == nx) ? ibeg + nlocx - 1 : ibeg + nlocx;
  ctx.uindices[1][1] = (jbeg + nlocy == ny) ? jbeg + nlocy - 1 : jbeg + nlocy;

  PetscReal maxdelta = 2.0 * eps;
  PetscInt  iter     = 0;
  while (iter < itermax && maxdelta > eps)
  {
    // Transfer data from global to local array.
    PetscCall(DMGlobalToLocalBegin(da, u_global, INSERT_VALUES, u_local));
    PetscCall(DMGlobalToLocalEnd(da, u_global, INSERT_VALUES, u_local));

    PetscCall(apply_jacobi(&ctx, da, u_local, u_global, maxdelta));

    // Get the global maxdelta.
    MPI_Allreduce(MPI_IN_PLACE, &maxdelta, 1, MPI_DOUBLE_PRECISION, MPI_MAX,
                  PETSC_COMM_WORLD);

    ++iter;
    PetscCall(PetscPrintf(PETSC_COMM_WORLD, "iter, maxdelta = %6d, %e\n", iter,
                          maxdelta));
  }
  // Save the last one.
  writeVTK(&ctx, da, u_global, 1);

  /* ----------------
   * END OF ITER
   * BEGIN CLEANUP
   * ----------------
   */
  PetscCall(DMDAVecRestoreArrayRead(cda, clocal, &ctx.coord));
  // Always a good idea to destroy everything.
  PetscCall(VecDestroy(&u_global));
  // Do NOT destroy clocal, it is a borrowed vector. It goes when the DM goes.
  // Also do NOT destroy cda. It's causing a segfault somewhere below.
  PetscCall(DMRestoreLocalVector(da, &u_local));
  PetscCall(DMDestroy(&da));

  PetscCall(PetscFinalize());
  PetscFunctionReturn(PETSC_SUCCESS);
}

//------------------------------------------------------------------------------
// from
// https://github.com/aadi-bh/parallel/blob/60bfbfa85302bd9cf97ccc123c99032cfb496173/mpi/poisson3d.cc#L493
/*
 * Writes out a VTK file with the local data in u_global. Corners and
 * coordinates provided by clocal.
 */
PetscErrorCode writeVTK(AppCtx *ctx, DM da, Vec u_global, PetscInt c)
{
  PetscMPIInt rank;
  MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
  char filename[64];
  snprintf(filename, 64, "sol_%d_%d.vtk", c, rank);

  Vec           u_local;
  PetscScalar **sol;
  PetscCall(DMGetLocalVector(da, &u_local));
  // This fills u_local with the latest values, including ghosts.
  PetscCall(DMGlobalToLocal(da, u_global, INSERT_VALUES, u_local));
  PetscCall(DMDAVecGetArrayRead(da, u_local, &sol));

  PetscInt ibeg, jbeg, nlocx, nlocy;
  // We are going to save the latest ghost values too.
  PetscCall(DMDAGetGhostCorners(da, &ibeg, &jbeg, NULL, &nlocx, &nlocy, NULL));

  // Set the data name to the name given to the vector.
  const char *buffer;
  PetscObjectGetName((PetscObject)u_global, &buffer);
  char objName[64];
  strncpy(objName, buffer, 64);
  for (int c = 0, n = strlen(objName); c < n; ++c)
  {
    if (objName[c] == ' ')
      objName[c] = '_';
  }

  ofstream fout;
  fout.open(filename);
  fout << "# vtk DataFile Version 3.0" << endl;
  fout << "Cartesian grid" << endl;
  fout << "ASCII" << endl;
  fout << "DATASET RECTILINEAR_GRID" << endl;

  // 1 for Z because this is 2D
  fout << "DIMENSIONS " << nlocx << " " << nlocy << " " << 1 << endl;

  fout << "X_COORDINATES " << nlocx << " float" << endl;
  for (PetscInt i = ibeg; i < ibeg + nlocx; ++i)
    fout << ctx->coord[jbeg+1][i].x << " ";
  fout << endl;

  fout << "Y_COORDINATES " << nlocy << " float" << endl;
  for (PetscInt j = jbeg; j < jbeg + nlocy; ++j)
    fout << ctx->coord[j][ibeg+1].y << " ";
  fout << endl;

  fout << "Z_COORDINATES " << 1 << " float" << endl;
  fout << 0.0 << endl;

  fout << "POINT_DATA " << nlocx * nlocy << endl;
  fout << "SCALARS " << objName << " double" << endl;
  fout << "LOOKUP_TABLE default" << endl;
  for (PetscInt j = jbeg; j < jbeg + nlocy; ++j)
  {
    for (PetscInt i = ibeg; i < ibeg + nlocx; ++i)
      fout << sol[j][i] << " ";
    fout << endl;
  }
  fout << endl;
  fout.close();

  PetscCall(DMDAVecRestoreArrayRead(da, u_local, &sol));
  PetscCall(DMRestoreLocalVector(da, &u_local));

  if(rank == 0)
  {
    // create file c.visit which contains list of vtk files
    char visit[64];
    snprintf(visit, 64, "%d.visit", c);

    PetscMPIInt nprocs;
    MPI_Comm_size(PETSC_COMM_WORLD, &nprocs);

    ofstream fvis;
    fvis.open(visit);
    fvis << "!NBLOCKS " << nprocs << endl;

    for(PetscMPIInt i=0; i<nprocs; ++i)
    {
      snprintf(filename, 64, "sol_%d_%d.vtk", c, i);
      fvis << filename << endl;
    }
    fvis.close();
    cout << "Wrote file " << visit << endl;
  }

  PetscFunctionReturn(PETSC_SUCCESS);
}
