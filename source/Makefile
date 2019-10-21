#program name
PROGS   = ../SQCF90

#using source code
SRC     = \
def_gate.f90 \
set_gate_list.f90 \
prod_gate.f90 \
main_exp.f90

#source name --> object name
OBJ      = $(SRC:.f90=.o)

#options for compiling
#######################################################
FC      = gfortran
FCFLAGS = -fopenmp -Ofast -march=native
# (for debug ) -Wall -pedantic -fbounds-check -O -Wuninitialized -ffpe-trap=invalid,zero,overflow -fbacktrace -fopenmp 
INCLUDE = -I /opt/OpenBLAS/include/ 
LIB     = -L/opt/OpenBLAS/lib
LINK    = -lopenblas -lpthread

#######################################################

LD      = $(FC)
LDFLAGS = $(FCFLAGS)

#suffixe
.SUFFIXES: .f90

all:$(PROGS)

$(PROGS):$(OBJ)
	$(LD) $(LDFLAGS) -o $@ $(OBJ) $(LIB) $(LINK)

.f90.o: 
	$(FC) $(FCFLAGS) -c $< $(INCLUDE)

clean: 
	rm -f $(OBJ) $(OBJ:.o=.mod) ./$(PROGS)
