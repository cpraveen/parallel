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

## integral.f90

## sum_reduction.f90
