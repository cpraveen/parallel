TARGETS = hello

all: $(TARGETS)

hello: hello.c
	mpicc -o hello hello.c -I$(PETSC_DIR)/include -L$(PETSC_DIR)/lib -lpetsc

clean:
	rm -rf $(TARGETS) *.o
