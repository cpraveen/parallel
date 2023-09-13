TARGETS = main

all: $(TARGETS)

main: main.cc
	mpic++ -o main main.cc -I$(PETSC_DIR)/include -L$(PETSC_DIR)/lib -lpetsc

clean:
	rm -rf $(TARGETS) *.o
