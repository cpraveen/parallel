ifndef PETSC_DIR
$(error PETSC_DIR is not set.)
endif

include ${PETSC_DIR}/lib/petsc/conf/variables
include ${PETSC_DIR}/lib/petsc/conf/rules

.DEFAULT_GOAL := main

main: main.o
	-${CLINKER} -o main main.o ${PETSC_LIB}

include ${PETSC_DIR}/lib/petsc/conf/test

clean::
	-${RM} -f main main.o solution_*.vts solution.pvd
