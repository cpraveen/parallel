# Solve 2d Poisson equation in parallel using Jacobi method

Uses PetSc to implement a parallel Jacobi solver for the Poisson equation in 2 dimensions.

Set PETSC location in variable `PETSC_DIR`.

```
make
mpirun -np 4 ./main
visit -o 1.visit
```

Run with different grid size

```
make
mpirun -np 4 ./main -da_grid_x 200 -da_grid_y 200
visit -o 1.visit
```

Change jacobi options like this

```
make
mpirun -np 4 ./main -jacobi_iter 1000 -jacobi_eps 1e-3
visit -o 1.visit
```
