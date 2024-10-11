FORTRAN=ifort
LINK=-lmkl_rt

FLAGS=-mp1 -mieee-fp -fimf-precision=high -O3 -heap-arrays

requisitos_f90 = global_variables.f90 processing_input.f90 source1D1.f90 source1D2.f90 \
		processing_outputs.f90 GEMM3D_EFA.f90 Main_ME.f90 

objetos = global_variables.o processing_input.o source1D1.o source1D2.o \
		processing_outputs.o GEMM3D_EFA.o Main_ME.o 

requisitos_grid3D_f90 = grid3D.f90

objetos_grid3D = grid3D.o

ME_3D_EFA.exec : $(objetos) grid3D.exec
	$(FORTRAN) $(FLAGS) -qopenmp -o ME_3D_EFA.exec $(objetos) $(LINK)

global_variables.o: global_variables.f90
	$(FORTRAN) -c global_variables.f90

processing_input.o: processing_input.f90
	$(FORTRAN) -c processing_input.f90

source1D1.o: source1D1.f90
	$(FORTRAN) -c source1D1.f90

source1D2.o: source1D2.f90
	$(FORTRAN) -qopenmp -c source1D2.f90

processing_outputs.o: processing_outputs.f90
	$(FORTRAN) -qopenmp -c processing_outputs.f90

GEMM3D_EFA.o: GEMM3D_EFA.f90
	$(FORTRAN) -c GEMM3D_EFA.f90 

Main_ME.o: Main_ME.f90
	$(FORTRAN) -c Main_ME.f90 

grid3D.exec :$(objetos_grid3D)
	$(FORTRAN) $(FLAGS) -o grid3D.exec $(objetos_grid3D)

grid3D.o: grid3D.f90
	$(FORTRAN) -c grid3D.f90


clean:	
	rm -rf *.o *.mod *.exec

