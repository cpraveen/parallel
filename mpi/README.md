# Parallel programs using MPI

## hello

This code just prints a message from every rank.

Compile the code

```shell
mpif77 -o hello hello.f
mpif90 -o hello hello.f90
mpicc  -o hello hello.c
```

If you have a 4-core cpu, run like this

```shell
mpirun -np 4 ./hello
```

## integrate1, using blocking send/recv

In this example, we compute the integral of cos(x) on [a,b]. We break the interval into as many ranks we have and compute integrals on each sub-interval. Then each rank sends its own result to the root process, where all the sub-integrals are added up to get the full integral.

Compile

```shell
mpif90 -o integrate1 integrate1.f90
mpicc  -o integrate1 integrate1.c
```

Run

```shell
mpirun -np 4 ./integrate1
```

## integrate2, using non-blocking recv

Same as integrate1, but uses non-blocking recv.

Compile

```shell
mpif90 -o integrate2 integrate2.f90
mpicc  -o integrate2 integrate2.c
```

Run

```shell
mpirun -np 4 ./integrate2
```


## integrate3, using reduce

Same as integrate1 and integrate2, but we use `MPI_Reduce` to sum the sub-integrals on the root process.

## linadv1d

Solves linear advection in 1d using Lax-Wendroff and periodic bc.

```shell
mpicc -o linadv1d linadv1d.c
mpirun -np 4 ./linadv1d
python linadv1d.py
```

## poisson3d

Solves Poisson equation in 3d using Jacobi method. Set some parameters in `poisson3d.in`.

```shell
mpif90 -o poisson3d poisson3d.f90
mpirun -np 4 ./poisson3d
```

There is also a C version in `poisson3d.c`.
