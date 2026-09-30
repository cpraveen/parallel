/* Solve Burger's equation in 1d using periodic bc.
   The grid is cell-centered. We use finite volume WENO scheme.
*/
static char help[] = "Solves 1d Burger equation.\n\n";

#include <petscsys.h>
#include <petscdm.h>
#include <petscdmda.h>
#include <petscvec.h>

#define max(a,b) ((a) > (b) ? (a) : (b))
#define min(a,b) ((a) < (b) ? (a) : (b))

const PetscReal ark[] = {0.0, 3.0/4.0, 1.0/3.0};
const PetscReal xmin = 0.0;
const PetscReal xmax = 1.0;

//------------------------------------------------------------------------------
// Initial condition
//------------------------------------------------------------------------------
PetscReal initcond(PetscReal x)
{
   return 1.0 + 0.5 * sin(2*M_PI*x);
}

//------------------------------------------------------------------------------
// Weno reconstruction of Jiang-Shu
// Return left value for face between u0, up1
//------------------------------------------------------------------------------
PetscReal weno5(PetscReal um2, PetscReal um1, PetscReal u0, PetscReal up1, PetscReal up2)
{
   PetscReal eps = 1.0e-6;
   PetscReal gamma1=1.0/10.0, gamma2=3.0/5.0, gamma3=3.0/10.0;
   PetscReal beta1, beta2, beta3;
   PetscReal u1, u2, u3;
   PetscReal w1, w2, w3;

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
// Flux for the Burger equation
//------------------------------------------------------------------------------
PetscReal Flux(PetscReal u)
{
   return 0.5*pow(u,2);
}
//------------------------------------------------------------------------------
// Godunov flux if sonic point is u = 0
//------------------------------------------------------------------------------
PetscReal numflux(PetscReal ul, PetscReal ur)
{
   return max(Flux(max(0,ul)),Flux(min(0,ur)));
}
//------------------------------------------------------------------------------
// Save solution to file
//------------------------------------------------------------------------------
PetscErrorCode savesol(const int nx, const PetscReal dx, Vec ug)
{
   int i, rank;
   static int count = 0;

   PetscFunctionBeginUser;
   MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

   // Gather entire solution on rank=0 process. This is bad thing to do
   // in a real application.
   VecScatter ctx;
   Vec        uall;
   PetscCall(VecScatterCreateToZero(ug,&ctx,&uall));
   // scatter as many times as you need
   PetscCall(VecScatterBegin(ctx,ug,uall,INSERT_VALUES,SCATTER_FORWARD));
   PetscCall(VecScatterEnd(ctx,ug,uall,INSERT_VALUES,SCATTER_FORWARD));
   // destroy scatter context and local vector when no longer needed
   PetscCall(VecScatterDestroy(&ctx));
   if(rank==0)
   {
      PetscScalar *uarray;
      PetscCall(VecGetArray(uall, &uarray));
      FILE *f;
      f = fopen("sol.dat","w");
      for(i=0; i<nx; ++i)
         fprintf(f, "%e %e\n", xmin+i*dx+0.5*dx, uarray[i]);
      fclose(f);
      printf("Wrote solution into sol.dat\n");
      PetscCall(VecRestoreArray(uall, &uarray));
      if(count==0)
      {
         // Initial solution is copied to sol0.dat
         system("cp sol.dat sol0.dat");
         count = 1;
      }
   }
   PetscCall(VecDestroy(&uall));
   PetscFunctionReturn(PETSC_SUCCESS);
}

//------------------------------------------------------------------------------
int main(int argc, char *argv[])
{
   DM       da;
   Vec      ug; // global vector
   Vec      ul; // local vector
   PetscInt i, ibeg, nloc;
   PetscInt nx=100;         // no. of cells, can change via command line
   const PetscInt sw = 3;   // stencil width
   const PetscInt ndof = 1; // no. of dofs per cell
   PetscMPIInt rank, size;
   PetscReal cfl = 0.4;
   PetscReal tfinal = 0.25;

   PetscFunctionBeginUser;
   PetscCall(PetscInitialize(&argc, &argv, (char*)0, help));

   MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
   MPI_Comm_size(PETSC_COMM_WORLD, &size);

   PetscCall(PetscOptionsGetReal(NULL,NULL,"-tfinal",&tfinal,NULL));
   PetscCall(PetscOptionsGetReal(NULL,NULL,"-cfl",&cfl,NULL));

   PetscCall(DMDACreate1d(PETSC_COMM_WORLD, DM_BOUNDARY_PERIODIC, nx, ndof, sw, NULL, &da));
   PetscCall(DMSetFromOptions(da));
   PetscCall(DMSetUp(da));

   PetscCall(DMCreateGlobalVector(da, &ug));

   PetscCall(DMDAGetCorners(da, &ibeg, 0, 0, &nloc, 0, 0));
   PetscCall(DMDAGetInfo(da,0,&nx,0,0,0,0,0,0,0,0,0,0,0));
   PetscReal dx = (xmax - xmin) / (PetscReal)(nx);
   PetscPrintf(PETSC_COMM_WORLD,"nx = %d, dx = %e\n", nx, dx);
   PetscReal umax_loc = 0.0;
   for(i=ibeg; i<ibeg+nloc; ++i)
   {
      PetscReal x = xmin + i*dx + 0.5*dx;
      PetscReal v = initcond(x);
      umax_loc = max(umax_loc,fabs(v)); // max wave speed
      PetscCall(VecSetValues(ug,1,&i,&v,INSERT_VALUES));
   }
   // Find max over all partitions
   PetscReal umax;
   MPI_Allreduce(&umax_loc, &umax, 1, MPI_DOUBLE, MPI_MAX, PETSC_COMM_WORLD);

   PetscCall(VecAssemblyBegin(ug));
   PetscCall(VecAssemblyEnd(ug));

   // save initial condition
   PetscCall(savesol(nx, dx, ug));

   // Get local view
   PetscCall(DMGetLocalVector(da, &ul));

   PetscInt il, nl;
   PetscCall(DMDAGetGhostCorners(da,&il,0,0,&nl,0,0));

   PetscReal res[nloc], uold[nloc];
   PetscReal dt = cfl * dx / (umax + 1.0e-13);
   PetscReal lam= dt/dx;
   PetscReal t = 0.0;

   while(t < tfinal)
   {
      if(t+dt > tfinal) // adjust dt to reach Tf
      {
         dt = tfinal - t;
         lam = dt/dx;
      }
      for(int rk=0; rk<3; ++rk) // loop for rk stages
      {
         PetscCall(DMGlobalToLocalBegin(da, ug, INSERT_VALUES, ul));
         PetscCall(DMGlobalToLocalEnd(da, ug, INSERT_VALUES, ul));

         PetscScalar *u;
         PetscCall(DMDAVecGetArrayRead(da, ul, &u));

         PetscScalar *unew;
         PetscCall(DMDAVecGetArray(da, ug, &unew));

         // First stage, store solution at time level n into uold
         if(rk==0)
            for(i=ibeg; i<ibeg+nloc; ++i) uold[i-ibeg] = u[i];

         for(i=0; i<nloc; ++i)
            res[i] = 0.0;

         // Loop over faces and compute flux
         for(i=0; i<nloc+1; ++i) // local index
         {
            // face between j-1, j
            PetscInt j   = il+sw+i; // global index
            PetscReal ul = weno5(u[j-3],u[j-2],u[j-1],u[j],u[j+1]);
            PetscReal ur = weno5(u[j+2],u[j+1],u[j],u[j-1],u[j-2]);
            PetscReal flux = numflux(ul, ur);
            if(i==0) // first face
            {
               res[i] -= flux;
            }
            else if(i==nloc) // last face
            {
               res[i-1] += flux;
            }
            else // intermediate faces
            {
               res[i]   -= flux;
               res[i-1] += flux;
            }
         }

         // Update solution
         for(i=ibeg; i<ibeg+nloc; ++i)
            unew[i] = ark[rk]*uold[i-ibeg] 
                      + (1.0-ark[rk])*(u[i] - lam * res[i-ibeg]);

         PetscCall(DMDAVecRestoreArrayRead(da, ul, &u));
         PetscCall(DMDAVecRestoreArray(da, ug, &unew));
      }

      t += dt;
      PetscPrintf(PETSC_COMM_WORLD,"t = %f\n", t);
   }

   PetscCall(savesol(nx, dx, ug));

   // Destroy everything before finishing
   PetscCall(DMDestroy(&da));
   PetscCall(VecDestroy(&ug));

   PetscCall(PetscFinalize());
   return 0;
}
