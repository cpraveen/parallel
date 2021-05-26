# OpenMP

## hello.f90

```
gfortran -o hello hello.f90 -fopenmp
./hello
OMP_NUM_THREADS=4 ./hello
```

## hello.cc

```
g++ -o hello hello.cc -fopenmp
./hello
OMP_NUM_THREADS=4 ./hello
```

## array_sum.f90
Add two arrays: c = a + b

## integral.f90
Compute integral of a function by quadrature

## sum_reduction.f90
Compute sum of an array

## jacobi.f90
Solve 2d Poisson equation on Cartesian grid using Jacobi iterations
