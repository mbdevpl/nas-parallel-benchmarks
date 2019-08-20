SHELL=/bin/sh
BENCHMARK=bt
BENCHMARKU=BT

include ../config/make.def


OBJS = bt.o  initialize.o exact_solution.o exact_rhs.o \
       set_constants.o adi.o  rhs.o      \
       x_solve.o y_solve.o   \
       z_solve.o add.o error.o verify.o \
       ${COMMON}/print_results.o ${COMMON}/timers.o ${COMMON}/wtime.o

INCLUDES = header.h cnst.h npbparams.h

include ../sys/make.common

#ADDRESSLENGTH =
#PRECISION     =

#FComplier     = HPFC
#CommLayer = -Mmpi
#F77	  = pghpf
#FFLAGS	  = -O3 $(PRECISION) $(ADDRESSLENGTH) $(CommLayer) \
#		-D$(FComplier) -Mautopar -Moverlap=size:2 -O3 -Minfo
#FLINKFLAGS     =  -O3 $(PRECISION) $(ADDRESSLENGTH) $(CommLayer)
#FLINK	= pghpf


# npbparams.h is included by header.h
# The following rule should do the trick but many make programs (not gmake)
# will do the wrong thing and rebuild the world every time (because the
# mod time on header.h is not changed. One solution would be to
# touch header.h but this might cause confusion if someone has
# accidentally deleted it. Instead, make the dependency on npbparams.h
# explicit in all the lines below (even though dependence is indirect).

# header.h: npbparams.h

${PROGRAM}: config ${OBJS}
	${FLINK} ${FLINKFLAGS} -o ${PROGRAM} ${OBJS} ${F_LIB}

.f.o:
	${FCOMPILE} $<


bt.o:             bt.f  $(INCLUDES)
initialize.o:     initialize.f  $(INCLUDES)
exact_solution.o: exact_solution.f  $(INCLUDES)
exact_rhs.o:      exact_rhs.f  $(INCLUDES)
set_constants.o:  set_constants.f  $(INCLUDES)
adi.o:            adi.f  $(INCLUDES)
rhs.o:            rhs.f  $(INCLUDES)
x_solve.o:        x_solve.f  $(INCLUDES)
y_solve.o:        y_solve.f  $(INCLUDES)
z_solve.o:        z_solve.f  $(INCLUDES)
add.o:            add.f $(INCLUDES)
error.o:          error.f $(INCLUDES)
verify.o:         verify.f $(INCLUDES)

clean:
	- rm -f *.o *~ mputil*
	- rm -f  npbparams.h core *.d
