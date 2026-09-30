/* Solves 2d linear advection equation with periodic bc.
   The grid is cell-centered.
*/
static char help[] = "Solves u_t + u_x + u_y = 0.\n\n";

#include <petscsys.h>
#include <petscdm.h>
#include <petscdmda.h>
#include <petscvec.h>

const double ark[] = {0.0, 3.0/4.0, 1.0/3.0};
const double xmin = 0.0, xmax = 1.0;
const double ymin = 0.0, ymax = 1.0;
double dx, dy;

//------------------------------------------------------------------------------
// Initial condition
//------------------------------------------------------------------------------
double initcond(double x, double y)
{
   return sin(2*M_PI*x) * sin(2*M_PI*y);
}

//------------------------------------------------------------------------------
// Weno reconstruction of Jiang and Shu
// Return left value for face between u0, up1
//------------------------------------------------------------------------------
double weno5(double um2, double um1, double u0, double up1, double up2)
{
   double eps = 1.0e-6;
   double gamma1=1.0/10.0, gamma2=3.0/5.0, gamma3=3.0/10.0;
   double beta1, beta2, beta3;
   double u1, u2, u3;
   double w1, w2, w3;

   beta1 = (13.0/12.0)*pow((um2 - 2.0*um1 + u0),2) +
           (1.0/4.0)*pow((um2 - 4.0*um1 + 3.0*u0),2);
   beta2 = (13.0/12.0)*pow((um1 - 2.0*u0 + up1),2) +
           (1.0/4.0)*pow((um1 - up1),2);
   beta3 = (13.0/12.0)*pow((u0 - 2.0*up1 + up2),2) +
           (1.0/4.0)*pow((3.0*u0 - 4.0*up1 + up2),2);

   w1 = gamma1 / pow(eps+beta1, 2);
   w2 = gamma2 / pow(eps+beta2, 2);
   w3 = gamma3 / pow(eps+beta3, 2);

   u1 = (1.0/3.0)*um2 - (7.0/6.0)*um1 + (11.0/6.0)*u0;
   u2 = -(1.0/6.0)*um1 + (5.0/6.0)*u0 + (1.0/3.0)*up1;
   u3 = (1.0/3.0)*u0 + (5.0/6.0)*up1 - (1.0/6.0)*up2;

   return (w1 * u1 + w2 * u2 + w3 * u3)/(w1 + w2 + w3);
}

//------------------------------------------------------------------------------
// Each rank saves its own solution to a tecplot file
//------------------------------------------------------------------------------
PetscErrorCode savesol(int *c, double t, DM da, Vec ug)
{
   char           filename[32] = "sol";
   PetscMPIInt    rank;
   PetscInt       i, j, nx, ny, ibeg, jbeg, nlocx, nlocy;
   FILE           *fp;
   Vec            ul;
   PetscScalar    **u;

   PetscFunctionBeginUser;
   PetscCall(DMGetLocalVector(da, &ul));
   PetscCall(DMGlobalToLocalBegin(da, ug, INSERT_VALUES, ul));
   PetscCall(DMGlobalToLocalEnd(da, ug, INSERT_VALUES, ul));
   PetscCall(DMDAVecGetArray(da, ul, &u));
   PetscCall(DMDAGetInfo(da,0,&nx,&ny,0,0,0,0,0,0,0,0,0,0));
   PetscCall(DMDAGetCorners(da, &ibeg, &jbeg, 0, &nlocx, &nlocy, 0));

   PetscInt iend = PetscMin(ibeg+nlocx+1, nx);
   PetscInt jend = PetscMin(jbeg+nlocy+1, ny);

   MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
   sprintf(filename, "sol-%03d-%03d.plt", *c, rank);
   fp = fopen(filename,"w");
   fprintf(fp, "TITLE = \"u_t + u_x + u_y = 0\"\n");
   fprintf(fp, "VARIABLES = x, y, sol\n");
   fprintf(fp, "ZONE STRANDID=1, SOLUTIONTIME=%e, I=%d, J=%d, DATAPACKING=POINT\n",
           t, iend-ibeg, jend-jbeg);
   for(j=jbeg; j<jend; ++j)
      for(i=ibeg; i<iend; ++i)
      {
         PetscReal x = xmin + i*dx + 0.5*dx;
         PetscReal y = ymin + j*dy + 0.5*dy;
         fprintf(fp, "%e %e %e\n", x, y, u[j][i]);
      }
   fclose(fp);

   PetscCall(DMDAVecRestoreArray(da, ul, &u));
   PetscCall(DMRestoreLocalVector(da, &ul));

   ++(*c);
   PetscFunctionReturn(PETSC_SUCCESS);
}
//------------------------------------------------------------------------------
int main(int argc, char *argv[])
{
   // some parameters that can overwritten from command line
   PetscReal Tf  = 10.0;
   PetscReal cfl = 0.4;
   PetscInt  si  = 100;
   PetscInt  nx  = 50, ny = 50; // use -da_grid_x, -da_grid_y to override these

   DM       da;
   Vec      ug, ul;
   PetscScalar **u;
   PetscInt i, j, ibeg, jbeg, nlocx, nlocy;
   const PetscInt sw = 3, ndof = 1; // stencil width
   PetscMPIInt rank, size;
   int c = 0; // counter for saving solution files

   PetscFunctionBeginUser;
   PetscCall(PetscInitialize(&argc, &argv, (char*)0, help));

   MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
   MPI_Comm_size(PETSC_COMM_WORLD, &size);

   // Get some command line options
   PetscCall(PetscOptionsGetReal(NULL,NULL,"-Tf",&Tf,NULL));
   PetscCall(PetscOptionsGetReal(NULL,NULL,"-cfl",&cfl,NULL));
   PetscCall(PetscOptionsGetInt(NULL,NULL,"-si",&si,NULL));

   PetscCall(DMDACreate2d(PETSC_COMM_WORLD, DM_BOUNDARY_PERIODIC, DM_BOUNDARY_PERIODIC,
                          DMDA_STENCIL_BOX, nx, ny, PETSC_DECIDE, PETSC_DECIDE, ndof,
                          sw, NULL, NULL, &da));
   PetscCall(DMSetFromOptions(da));
   PetscCall(DMSetUp(da));
   PetscCall(DMDASetUniformCoordinates(da,xmin,xmax,ymin,ymax,0.0,0.0));

   PetscCall(DMDAGetInfo(da,0,&nx,&ny,0,0,0,0,0,0,0,0,0,0));
   dx = (xmax - xmin) / (PetscReal)(nx);
   dy = (ymax - ymin) / (PetscReal)(ny);
   PetscPrintf(PETSC_COMM_WORLD,"nx = %d, dx = %e\n", nx, dx);
   PetscPrintf(PETSC_COMM_WORLD,"ny = %d, dy = %e\n", ny, dy);

   PetscCall(DMCreateGlobalVector(da, &ug));
   PetscCall(PetscObjectSetName((PetscObject) ug, "Solution"));

   PetscCall(DMDAGetCorners(da, &ibeg, &jbeg, 0, &nlocx, &nlocy, 0));
   PetscCall(DMDAVecGetArray(da, ug, &u));
   for(j=jbeg; j<jbeg+nlocy; ++j)
      for(i=ibeg; i<ibeg+nlocx; ++i)
      {
         PetscReal x = xmin + i*dx + 0.5*dx;
         PetscReal y = ymin + j*dy + 0.5*dy;
         u[j][i] = initcond(x,y);
      }
   PetscCall(DMDAVecRestoreArray(da, ug, &u));
   PetscCall(savesol(&c, 0.0, da, ug));

   // Get local view
   PetscCall(DMGetLocalVector(da, &ul));

   PetscInt il, jl, nl, ml;
   PetscCall(DMDAGetGhostCorners(da,&il,&jl,0,&nl,&ml,0));

   // Allocate res[nlocy][nlocx] and uold[nlocy][nlocx]
   PetscReal (*res) [nlocx] = calloc(nlocy, sizeof(*res) );
   PetscReal (*uold)[nlocx] = calloc(nlocy, sizeof(*uold));
   PetscReal umax = sqrt(2.0);
   PetscReal dt = cfl * dx / umax;
   PetscReal lam= dt/(dx*dy);
   PetscReal t = 0.0;
   int it = 0;

   while(t < Tf)
   {
      if(t+dt > Tf) // adjust dt to reach Tf
      {
         dt = Tf - t;
         lam = dt/(dx*dy);
      }
      for(int rk=0; rk<3; ++rk) // loop for rk stages
      {
         PetscCall(DMGlobalToLocalBegin(da, ug, INSERT_VALUES, ul));
         PetscCall(DMGlobalToLocalEnd(da, ug, INSERT_VALUES, ul));

         PetscCall(DMDAVecGetArrayRead(da, ul, &u));

         PetscScalar **unew;
         PetscCall(DMDAVecGetArray(da, ug, &unew));

         if(rk==0)
         {
            for(j=jbeg; j<jbeg+nlocy; ++j)
               for(i=ibeg; i<ibeg+nlocx; ++i)
                  uold[j-jbeg][i-ibeg] = u[j][i];
         }

         for(j=0; j<nlocy; ++j)
            for(i=0; i<nlocx; ++i)
               res[j][i] = 0.0;

         // x fluxes
         for(i=0; i<nlocx+1; ++i)
            for(j=0; j<nlocy; ++j)
            {
               // face between k-1, k
               PetscInt k = il+sw+i;
               PetscInt l = jl+sw+j;
               PetscReal uleft = weno5(u[l][k-3],u[l][k-2],u[l][k-1],u[l][k],u[l][k+1]);
               PetscReal flux = uleft * dy;
               if(i==0)
               {
                  res[j][i] -= flux;
               }
               else if(i==nlocx)
               {
                  res[j][i-1] += flux;
               }
               else
               {
                  res[j][i]   -= flux;
                  res[j][i-1] += flux;
               }
            }

         // y fluxes
         for(j=0; j<nlocy+1; ++j)
            for(i=0; i<nlocx; ++i)
            {
               // face between l-1, l
               PetscInt k = il+sw+i;
               PetscInt l = jl+sw+j;
               PetscReal uleft = weno5(u[l-3][k],u[l-2][k],u[l-1][k],u[l][k],u[l+1][k]);
               PetscReal flux = uleft * dx;
               if(j==0)
               {
                  res[j][i] -= flux;
               }
               else if(j==nlocy)
               {
                  res[j-1][i] += flux;
               }
               else
               {
                  res[j][i]   -= flux;
                  res[j-1][i] += flux;
               }
            }

         // Update solution
         for(j=jbeg; j<jbeg+nlocy; ++j)
            for(i=ibeg; i<ibeg+nlocx; ++i)
               unew[j][i] = ark[rk]*uold[j-jbeg][i-ibeg]
                            + (1.0-ark[rk])*(u[j][i] - lam * res[j-jbeg][i-ibeg]);

         PetscCall(DMDAVecRestoreArrayRead(da, ul, &u));
         PetscCall(DMDAVecRestoreArray(da, ug, &unew));
      }

      t += dt; ++it;
      PetscPrintf(PETSC_COMM_WORLD,"it, t = %d, %f\n", it, t);
      if(it%si == 0 || PetscAbs(t-Tf) < 1.0e-13)
      {
         PetscCall(savesol(&c, t, da, ug));
      }
   }

   // Destroy everything before finishing
   PetscCall(VecDestroy(&ug));
   PetscCall(DMRestoreLocalVector(da, &ul));
   PetscCall(DMDestroy(&da));

   free(res); free(uold);

   PetscCall(PetscFinalize());
   return 0;
}
