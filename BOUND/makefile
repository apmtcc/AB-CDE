#Macro definition
FC = mpif90
FFLAGS = -O2

#end of Macro definition

bound.exe: modules.o main.o sub.o bf.o lobatto.o potential.o potdriv.o \
	jacobix.o specfun.o rs.o hamiltonian.o  trans_blas.o \
	fnh3ion_finn.o arpack1.o arpack2.o arpack0.o
 

	$(FC) $(FFLAGS) -o bound.exe modules.o main.o sub.o bf.o lobatto.o \
	potential.o potdriv.o jacobix.o specfun.o rs.o hamiltonian.o \
        fnh3ion_finn.o arpack1.o arpack2.o arpack0.o trans_blas.o  \
        -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core \
        -lmkl_blas95_lp64 -lmkl_lapack95_lp64 -liomp5

clean:
	rm *.o *.mod *.exe

