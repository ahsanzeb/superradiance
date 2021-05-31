#################################################################################
# makefile: superradiance in organic microcavities
#################################################################################
#	info:
#	    organic molecules are 3 level systems
#	    Hilbert space: 
#	         N-molecules permutation symmetric states
#	              times
#	         nph cavity states
#	       = {N+2 C 2} * nph = (N+2)*(N+1)*nph 
#                (instead of 3^N*nph states in the full Hilbert space)
#################################################################################

FC 	= mpif90  #gfortran
FFLAGS= -O3 -g -fbounds-check -fbacktrace
LOADLIBES = libarpack_mac.a -llapack

# to link statically we might need something like: 
#FFLAGS += -static
#LOADLIBES += -lblas -lgfortran -lm


SRC = main.f basis.f modmain.f maps.f hamiltonian.f diag.f

OBJ = $(SRC:.f=.o)

.SUFFIXES: .o .f

super: $(OBJ)
	$(FC) -o super $(FFLAGS) $(OBJ) $(LOADLIBES)
	cp super ~/bin/.


# Here are the compile steps
# IMPORTANT: remove unused dependencies after testing all modules.
 
modmain.o: modmain.f  
	 $(FC) -c $(FFLAGS) modmain.f 

maps.o: maps.f modmain.o
	 $(FC) -c $(FFLAGS) maps.f 

basis.o: basis.f modmain.o maps.o
	 $(FC) -c $(FFLAGS) basis.f 	

hamiltonian.o: hamiltonian.f modmain.o maps.o basis.o 
	 $(FC) -c $(FFLAGS) hamiltonian.f 	


diag.o: diag.f modmain.o
	 $(FC) -c $(FFLAGS) diag.f 	


#dm.o: dm.f modmain.o
#	 $(FC) -c $(FFLAGS) dm.f 


main.o: main.f modmain.o maps.o basis.o hamiltonian.o diag.o
	 $(FC) -c $(FFLAGS) main.f	

clean:
	rm -rf *.o *.mod super
  
