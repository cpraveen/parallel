/*
Laplacian in 2D. Modeled by the partial differential equation

   -div grad u = f,  0 < x,y < 1,

with forcing function

   f = e^{-x^2/\nu} e^{-y^2/\nu}

with Dirichlet boundary conditions

   u = g(x,y) for x = 0, x = 1, y = 0, y = 1

or pure Neumman boundary conditions

This uses multigrid to solve the linear system.
The grid is vertex-centered so there are grid points on the boundary.
*/

static char help[] = "Solves 2D Laplacian using multigrid.\n\n";

#include <petscdm.h>
#include <petscdmda.h>
#include <petscksp.h>

extern PetscErrorCode ComputeMatrix(KSP, Mat, Mat, void*);
extern PetscErrorCode ComputeRHS(KSP, Vec, void*);

typedef enum
{
   DIRICHLET,
   NEUMANN
} BCType;

typedef struct
{
   PetscReal nu;
   BCType    bcType;
} UserContext;

PetscReal BoundaryValue(const PetscReal x, const PetscReal y, UserContext* user)
{
   return PetscExpScalar(-x * x / user->nu) * PetscExpScalar(-y * y / user->nu);
}

PetscReal RHSFun(const PetscReal x, const PetscReal y, UserContext* user)
{
   return PetscExpScalar(-x * x / user->nu) * PetscExpScalar(-y * y / user->nu);
}

int
main(int argc, char** argv)
{
   KSP         ksp;
   DM          da;
   Vec         x;
   UserContext user;
   const char* bcTypes[2] = {"dirichlet", "neumann"};
   PetscInt    bc;

   PetscFunctionBeginUser;
   PetscCall(PetscInitialize(&argc, &argv, NULL, help));
   PetscCall(KSPCreate(PETSC_COMM_WORLD, &ksp));
   PetscCall(DMDACreate2d(PETSC_COMM_WORLD, DM_BOUNDARY_NONE, DM_BOUNDARY_NONE,
                          DMDA_STENCIL_STAR, 3, 3, PETSC_DECIDE, PETSC_DECIDE,
                          1, 1, 0, 0, &da));
   PetscCall(DMDASetInterpolationType(da, DMDA_Q1));
   PetscCall(DMSetFromOptions(da));
   PetscCall(DMSetUp(da));
   PetscCall(DMDASetUniformCoordinates(da, 0, 1, 0, 1, 0, 0));
   PetscCall(DMDASetFieldName(da, 0, "Pressure"));

   user.nu = 0.1;
   bc = (PetscInt)DIRICHLET;
   user.bcType = (BCType)bc;

   PetscOptionsBegin(PETSC_COMM_WORLD, "", "Options for Poisson equation", "DMqq");
   PetscCall(PetscOptionsReal("-nu", "The width of the Gaussian source", "ex29.c", user.nu, &user.nu, NULL));
   PetscCall(PetscOptionsEList("-bc_type", "Type of boundary condition", "ex29.c", bcTypes, 2, bcTypes[0], &bc, NULL));
   PetscOptionsEnd();

   PetscCall(KSPSetComputeRHS(ksp, ComputeRHS, &user));
   PetscCall(KSPSetComputeOperators(ksp, ComputeMatrix, &user));
   PetscCall(KSPSetDM(ksp, da));
   PetscCall(KSPSetFromOptions(ksp));
   PetscCall(KSPSetUp(ksp));
   PetscCall(KSPSolve(ksp, NULL, NULL));

   PetscCall(KSPGetSolution(ksp, &x));
   PetscCall(PetscObjectSetName((PetscObject)x, "Pressure"));
   PetscCall(VecViewFromOptions(x, NULL, "-solution_view"));
   PetscCall(VecDestroy(&x));

   PetscCall(KSPDestroy(&ksp));
   PetscCall(DMDestroy(&da));
   PetscCall(PetscFinalize());
   return 0;
}

PetscErrorCode
ComputeRHS(KSP ksp, Vec b, void* ctx)
{
   UserContext  *user = (UserContext*)ctx;
   PetscInt      i, j, mx, my, xm, ym, xs, ys;
   PetscScalar   Hx, Hy, HydHx, HxdHy, x, y;
   PetscScalar** array;
   DM            da;

   PetscFunctionBeginUser;
   PetscCall(KSPGetDM(ksp, &da));
   PetscCall(DMDAGetInfo(da, 0, &mx, &my, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0));
   Hx    = 1.0 / (PetscReal)(mx - 1);
   Hy    = 1.0 / (PetscReal)(my - 1);
   HxdHy = Hx / Hy;
   HydHx = Hy / Hx;
   PetscCall(DMDAGetCorners(da, &xs, &ys, 0, &xm, &ym, 0));
   PetscCall(DMDAVecGetArray(da, b, &array));
   for (j = ys; j < ys + ym; j++)
   {
      for (i = xs; i < xs + xm; i++)
      {
         x = i * Hx;
         y = j * Hy;
         if (user->bcType == DIRICHLET && (i == 0 || j == 0 || i == mx - 1 || j == my - 1))
         {
            array[j][i] = BoundaryValue(x, y, user)  * 2.0 * (HxdHy + HydHx);
         }
         else
         {
            array[j][i] = RHSFun(x, y, user) * Hx * Hy;
         }
      }
   }
   PetscCall(DMDAVecRestoreArray(da, b, &array));
   PetscCall(VecAssemblyBegin(b));
   PetscCall(VecAssemblyEnd(b));

   /* force right-hand side to be consistent for singular matrix */
   /* note this is really a hack, normally the model would provide you with a consistent right handside */
   if (user->bcType == NEUMANN)
   {
      MatNullSpace nullspace;

      PetscCall(MatNullSpaceCreate(PETSC_COMM_WORLD, PETSC_TRUE, 0, 0, &nullspace));
      PetscCall(MatNullSpaceRemove(nullspace, b));
      PetscCall(MatNullSpaceDestroy(&nullspace));
   }
   PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode
ComputeMatrix(KSP ksp, Mat J, Mat jac, void* ctx)
{
   UserContext *user = (UserContext*)ctx;
   PetscInt     i, j, mx, my, xm, ym, xs, ys;
   PetscScalar  v[5];
   PetscReal    Hx, Hy, HydHx, HxdHy;
   MatStencil   row, col[5];
   DM           da;

   PetscFunctionBeginUser;
   PetscCall(KSPGetDM(ksp, &da));
   PetscCall(DMDAGetInfo(da, 0, &mx, &my, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0));
   Hx    = 1.0 / (PetscReal)(mx - 1);
   Hy    = 1.0 / (PetscReal)(my - 1);
   HxdHy = Hx / Hy;
   HydHx = Hy / Hx;
   PetscCall(DMDAGetCorners(da, &xs, &ys, 0, &xm, &ym, 0));
   for (j = ys; j < ys + ym; j++)
   {
      for (i = xs; i < xs + xm; i++)
      {
         row.i = i;
         row.j = j;
         if (i == 0 || j == 0 || i == mx - 1 || j == my - 1)
         {
            if (user->bcType == DIRICHLET)
            {
               v[0] = 2.0 * (HxdHy + HydHx);
               PetscCall(MatSetValuesStencil(jac, 1, &row, 1, &row, v, INSERT_VALUES));
            }
            else if (user->bcType == NEUMANN)
            {
               PetscInt numx = 0, numy = 0, num = 0;
               if (j != 0)
               {
                  v[num]     = -HxdHy;
                  col[num].i = i;
                  col[num].j = j - 1;
                  numy++;
                  num++;
               }
               if (i != 0)
               {
                  v[num]     = -HydHx;
                  col[num].i = i - 1;
                  col[num].j = j;
                  numx++;
                  num++;
               }
               if (i != mx - 1)
               {
                  v[num]     = -HydHx;
                  col[num].i = i + 1;
                  col[num].j = j;
                  numx++;
                  num++;
               }
               if (j != my - 1)
               {
                  v[num]     = -HxdHy;
                  col[num].i = i;
                  col[num].j = j + 1;
                  numy++;
                  num++;
               }
               v[num]     = numx * HydHx + numy * HxdHy;
               col[num].i = i;
               col[num].j = j;
               num++;
               PetscCall(MatSetValuesStencil(jac, 1, &row, num, col, v, INSERT_VALUES));
            }
         }
         else
         {
            v[0]     = -HxdHy;
            col[0].i = i;
            col[0].j = j - 1;
            v[1]     = -HydHx;
            col[1].i = i - 1;
            col[1].j = j;
            v[2]     = 2.0 * (HxdHy + HydHx);
            col[2].i = i;
            col[2].j = j;
            v[3]     = -HydHx;
            col[3].i = i + 1;
            col[3].j = j;
            v[4]     = -HxdHy;
            col[4].i = i;
            col[4].j = j + 1;
            PetscCall(MatSetValuesStencil(jac, 1, &row, 5, col, v, INSERT_VALUES));
         }
      }
   }
   PetscCall(MatAssemblyBegin(jac, MAT_FINAL_ASSEMBLY));
   PetscCall(MatAssemblyEnd(jac, MAT_FINAL_ASSEMBLY));
   if (user->bcType == NEUMANN)
   {
      MatNullSpace nullspace;
      PetscCall(MatNullSpaceCreate(PETSC_COMM_WORLD, PETSC_TRUE, 0, 0, &nullspace));
      PetscCall(MatSetNullSpace(J, nullspace));
      PetscCall(MatNullSpaceDestroy(&nullspace));
   }
   PetscFunctionReturn(PETSC_SUCCESS);
}

/*TEST

   test:
      args: -pc_type mg -pc_mg_type full -ksp_type fgmres -ksp_monitor_short -da_refine 8 -ksp_rtol 1.e-3

   test:
      suffix: 2
      args: -bc_type neumann -pc_type mg -pc_mg_type full -ksp_type fgmres -ksp_monitor_short -da_refine 8 -mg_coarse_pc_factor_shift_type nonzero
      requires: !single

   test:
      suffix: telescope
      nsize: 4
      args: -ksp_monitor_short -da_grid_x 257 -da_grid_y 257 -pc_type mg -pc_mg_galerkin pmat -pc_mg_levels 4 -ksp_type richardson -mg_levels_ksp_type chebyshev -mg_levels_pc_type jacobi -mg_coarse_pc_type telescope -mg_coarse_pc_telescope_ignore_kspcomputeoperators -mg_coarse_telescope_pc_type mg -mg_coarse_telescope_pc_mg_galerkin pmat -mg_coarse_telescope_pc_mg_levels 3 -mg_coarse_telescope_mg_levels_ksp_type chebyshev -mg_coarse_telescope_mg_levels_pc_type jacobi -mg_coarse_pc_telescope_reduction_factor 4

   test:
      suffix: 3
      args: -ksp_view -da_refine 2 -pc_type mg -pc_mg_distinct_smoothup -mg_levels_up_pc_type jacobi

   test:
      suffix: 4
      args: -ksp_view -da_refine 2 -pc_type mg -pc_mg_distinct_smoothup -mg_levels_up_ksp_max_it 3 -mg_levels_ksp_max_it 4

   testset:
     suffix: aniso
     args: -da_grid_x 10 -da_grid_y 2 -da_refine 2 -pc_type mg -ksp_monitor_short -mg_levels_ksp_max_it 6 -mg_levels_pc_type jacobi
     test:
       suffix: first
       args: -mg_levels_ksp_chebyshev_kind first
     test:
       suffix: fourth
       args: -mg_levels_ksp_chebyshev_kind fourth
     test:
       suffix: opt_fourth
       args: -mg_levels_ksp_chebyshev_kind opt_fourth

   test:
      suffix: 5
      nsize: 2
      requires: hypre !complex
      args: -pc_type mg -da_refine 2 -ksp_monitor -matptap_via hypre -pc_mg_galerkin both

   test:
      suffix: 6
      args: -pc_type svd -pc_svd_monitor ::all

TEST*/
