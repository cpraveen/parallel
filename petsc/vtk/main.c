#include <petscdmda.h>
#include <petscviewer.h>
#include <petscsys.h>
#include <petscmath.h>
#include <stdio.h>

static char help[] = "Saves DMDA solutions to VTK (.vts/.pvd) files.\n\n";

/*
 * Saves the DMDA solution vector to a VTK file.
 *
 * da       : The DMDA context
 * X        : The solution vector
 * filename : The output .vts filename
 */
PetscErrorCode SaveSolutionVTK(DM da, Vec X, const char *filename)
{
   PetscViewer viewer;

   PetscFunctionBeginUser;

   /* Open a VTK viewer in write mode */
   PetscCall(PetscViewerVTKOpen(PetscObjectComm((PetscObject)da), filename,
                                FILE_MODE_WRITE, &viewer));

   /* Write the vector to the viewer */
   PetscCall(VecView(X, viewer));

   /* Clean up the viewer */
   PetscCall(PetscViewerDestroy(&viewer));

   PetscFunctionReturn(0);
}

/*
 * Initializes the PVD file (Run this once at the start of your simulation)
 */
PetscErrorCode InitializePVD(MPI_Comm comm, const char *pvd_filename)
{
   PetscMPIInt rank;
   FILE        *fp;

   PetscFunctionBeginUser;
   MPI_Comm_rank(comm, &rank);

   if (rank == 0)
   {
      fp = fopen(pvd_filename, "w");
      if (fp)
      {
         fprintf(fp, "<?xml version=\"1.0\"?>\n");
         fprintf(fp, "<VTKFile type=\"Collection\" version=\"0.1\" "
                     "byte_order=\"LittleEndian\">\n");
         fprintf(fp, "  <Collection>\n");
         fclose(fp);
      }
   }
   PetscFunctionReturn(0);
}

/*
 * Appends a new timestep to the PVD file (Run this after saving each .vts file)
 */
PetscErrorCode AppendPVD(MPI_Comm comm, const char *pvd_filename,
                         const char *vts_filename, PetscReal time)
{
   PetscMPIInt rank;
   FILE        *fp;

   PetscFunctionBeginUser;
   MPI_Comm_rank(comm, &rank);

   if (rank == 0)
   {
      fp = fopen(pvd_filename, "a");
      if (fp)
      {
         fprintf(fp, "    <DataSet timestep=\"%f\" group=\"\" part=\"0\" "
                     "file=\"%s\"/>\n", (double)time, vts_filename);
         fclose(fp);
      }
   }
   PetscFunctionReturn(0);
}

/*
 * Closes the PVD file tags (Run this once at the end of your simulation)
 */
PetscErrorCode FinalizePVD(MPI_Comm comm, const char *pvd_filename)
{
   PetscMPIInt rank;
   FILE        *fp;

   PetscFunctionBeginUser;
   MPI_Comm_rank(comm, &rank);

   if (rank == 0)
   {
      fp = fopen(pvd_filename, "a");
      if (fp)
      {
         fprintf(fp, "  </Collection>\n");
         fprintf(fp, "</VTKFile>\n");
         fclose(fp);
      }
   }
   PetscFunctionReturn(0);
}

/*
 * Fills the global vector X with the exact 2D isentropic vortex solution of
 * the Euler equations at time `time`. The vortex translates with the
 * freestream, so its core sits at (uinf*time, vinf*time). The conserved
 * variables are ordered as (rho, rho*u, rho*v, E) with freestream
 * (rho, u, v, p) = (1, 1, 0, 1).
 */
PetscErrorCode SetExactVortex(DM da, Vec X, PetscReal time)
{
   DM                   cda;
   Vec                  coords;
   DMDACoor2d         **coors;
   PetscScalar        ***arr;
   PetscInt             i, j, xs, ys, xm, ym;
   PetscReal            gamma = 1.4;
   PetscReal            beta  = 5.0;
   PetscReal            uinf  = 1.0;
   PetscReal            vinf  = 0.0;

   PetscFunctionBeginUser;

   PetscCall(DMGetCoordinateDM(da, &cda));
   PetscCall(DMGetCoordinates(da, &coords));
   PetscCall(DMDAVecGetArray(cda, coords, &coors));
   PetscCall(DMDAVecGetArrayDOF(da, X, &arr));
   PetscCall(DMDAGetCorners(da, &xs, &ys, NULL, &xm, &ym, NULL));

   /* C = (gamma - 1)*beta^2/(8*gamma*pi^2) */
   const PetscReal C   = (gamma - 1.0)*beta*beta/
                         (8.0*gamma*PETSC_PI*PETSC_PI);

   for (j = ys; j < ys + ym; ++j)
   {
      for (i = xs; i < xs + xm; ++i)
      {
         /* Vortex core translated by the freestream: (uinf*t, vinf*t) */
         const PetscReal xp  = PetscRealPart(coors[j][i].x) - uinf*time;
         const PetscReal yp  = PetscRealPart(coors[j][i].y) - vinf*time;
         const PetscReal rsq = xp*xp + yp*yp;
         /* tc = exp((1 - r^2)/2) */
         const PetscReal tc  = PetscExpReal(0.5*(1.0 - rsq));
         /* T = 1 - C*exp(1 - r^2) */
         const PetscReal T   = 1.0 - C*PetscExpReal(1.0 - rsq);
         const PetscReal rho = PetscPowReal(T, 1.0/(gamma - 1.0));
         const PetscReal u   = uinf - beta/(2.0*PETSC_PI)*yp*tc;
         const PetscReal v   = vinf + beta/(2.0*PETSC_PI)*xp*tc;
         /* Isentropic: p = rho^gamma */
         const PetscReal p   = PetscPowReal(rho, gamma);
         const PetscReal E   = p/(gamma - 1.0) + 0.5*rho*(u*u + v*v);

         arr[j][i][0] = rho;
         arr[j][i][1] = rho*u;
         arr[j][i][2] = rho*v;
         arr[j][i][3] = E;
      }
   }

   PetscCall(DMDAVecRestoreArray(cda, coords, &coors));
   PetscCall(DMDAVecRestoreArrayDOF(da, X, &arr));

   PetscFunctionReturn(0);
}

int main(int argc, char **argv)
{
   DM             da;
   Vec            X;
   PetscInt       nvar = 4;
   PetscInt       mx = 64, my = 64;
   PetscReal      xmin = -5.0, xmax = 5.0, ymin = -5.0, ymax = 5.0;
   PetscReal      t, dt = 0.05, tfinal = 2.0;
   PetscInt       step, nsteps;
   char           filename[PETSC_MAX_PATH_LEN];
   const char    *varnames[] = {"Density", "MomentumX", "MomentumY", "Energy"};

   PetscCall(PetscInitialize(&argc, &argv, NULL, help));

   PetscCall(PetscOptionsGetReal(NULL, NULL, "-dt", &dt, NULL));
   PetscCall(PetscOptionsGetReal(NULL, NULL, "-tfinal", &tfinal, NULL));
   nsteps = (PetscInt)(tfinal/dt + 0.5);

   PetscCall(DMDACreate2d(PETSC_COMM_WORLD, DM_BOUNDARY_NONE, DM_BOUNDARY_NONE,
                          DMDA_STENCIL_STAR, mx, my, PETSC_DECIDE, PETSC_DECIDE,
                          nvar, 1, NULL, NULL, &da));
   PetscCall(DMSetFromOptions(da));
   PetscCall(DMSetUp(da));

   /* Physical coordinates, so the VTK output has the correct extents */
   PetscCall(DMDASetUniformCoordinates(da, xmin, xmax, ymin, ymax, 0.0, 0.0));

   for (PetscInt i = 0; i < nvar; ++i)
   {
      PetscCall(DMDASetFieldName(da, i, varnames[i]));
   }

   PetscCall(DMCreateGlobalVector(da, &X));

   /* Name the vector so the VTK arrays are labeled U.<field> */
   PetscCall(PetscObjectSetName((PetscObject)X, "U"));

   PetscCall(InitializePVD(PETSC_COMM_WORLD, "solution.pvd"));
   for (step = 0; step <= nsteps; ++step)
   {
      t = step*dt;
      PetscCall(SetExactVortex(da, X, t));
      PetscCall(PetscSNPrintf(filename, sizeof(filename),
                              "solution_%04d.vts", step));
      PetscCall(SaveSolutionVTK(da, X, filename));
      PetscCall(AppendPVD(PETSC_COMM_WORLD, "solution.pvd", filename, t));
   }
   PetscCall(FinalizePVD(PETSC_COMM_WORLD, "solution.pvd"));

   PetscCall(VecDestroy(&X));
   PetscCall(DMDestroy(&da));
   PetscCall(PetscFinalize());
   return 0;
}
