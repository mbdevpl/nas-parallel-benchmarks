# Makefile for OMP dummy library routine.

OMP_LIB_DUM = omp_lib_dumf.o omp_lib_dumc.o

libomp_lib_dum.a: $(OMP_LIB_DUM)
	ar r $@ $(OMP_LIB_DUM)

omp_lib_dumf.o: omp_lib_dumf.f
	$(F77) -c $(FFLAGS) omp_lib_dumf.f

omp_lib_dumc.o: omp_lib_dumc.c
	$(CC) -c $(CFLAGS) omp_lib_dumc.c

clean:
	-rm -f *~ *.o *.a
