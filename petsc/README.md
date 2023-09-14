# PETSc

We need to know where Petsc is installed, set the variable `PETSC_DIR` to point to the location of your Petsc installation. The makefiles will use this variable.

```shell
make hello
./hello -help
./hello -help intro
./hello
mpirun -n 4 ./hello
mpiexec -n 4 ./hello
```
## Makefile

There is a standard makefile included in Petsc which can be used to compile your code. E.g., the `hello.c` can be compiled as

```
make -f $PETSC_DIR/share/petsc/Makefile.user hello
```

This will create the executable file `hello`. You can modify this makefile for more complex settings, see the comments in that file.
