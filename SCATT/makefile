#Macro definition
FC = mpif90
FFLAGS = -O2

#end of Macro definition

scatt.exe: modules.o main.o sub.o bf.o potential.o potdriv.o jacobix.o specfun.o rs.o \
	hamiltonian.o  trans_blas.o fnh3_finn.o init.o \
	analysis.o progation.o lib_mkl.o 

	$(FC) $(FFLAGS) -o scatt.exe modules.o main.o sub.o bf.o \
	potential.o potdriv.o jacobix.o specfun.o rs.o hamiltonian.o \
        fnh3_finn.o init.o analysis.o progation.o lib_mkl.o \
	trans_blas.o

clean:
	rm *.o *.mod *.exe

