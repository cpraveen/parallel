# Solve 2d Poisson equation in parallel using Jacobi method

Uses PetSc to implement a parallel Jacobi solver for the Poisson equation in 2 dimensions.

Set PETSC location in variable `PETSC_DIR`.

```
make
mpirun -np 4 ./main
visit -o 1.visit
```
