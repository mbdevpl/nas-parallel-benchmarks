SHELL=/bin/sh
BENCHMARK=lu
BENCHMARKU=LU

include ../config/make.def

OBJS = lu.o read_input.o \
       domain.o setcoeff.o setbv.o exact.o setiv.o \
       erhs.o ssor.o calcnp.o indexbounds.o\
       rhs.o l2norm.o jacld.o blts.o jacu.o buts.o error.o \
       pintgr.o verify.o ${COMMON}/print_results.o \
       ${COMMON}/timers.o ${COMMON}/wtime.o

INCLUDES = applu.incl l2norm.h npbparams.h \
		blts.h indexbounds.h loopcontrol.h projarea.h buts.h

include ../sys/make.common


# npbparams.h is included by applu.incl
# The following rule should do the trick but many make programs (not gmake)
# will do the wrong thing and rebuild the world every time (because the
# mod time on header.h is not changed. One solution would be to
# touch header.h but this might cause confusion if someone has
# accidentally deleted it. Instead, make the dependency on npbparams.h
# explicit in all the lines below (even though dependence is indirect).

# applu.incl: npbparams.h

${PROGRAM}: config ${OBJS}
	${FLINK} ${FLINKFLAGS} -o ${PROGRAM} ${OBJS} ${F_LIB}

.f.o :
	${FCOMPILE} $<

lu.o:		lu.f $(INCLUDES)
blts.o:		blts.f $(INCLUDES)
buts.o:		buts.f $(INCLUDES)
calcnp.o:	calcnp.f $(INCLUDES)
erhs.o:		erhs.f $(INCLUDES)
error.o:	error.f $(INCLUDES)
exact.o:	exact.f $(INCLUDES)
jacld.o:	jacld.f $(INCLUDES)
indexbounds.o:	indexbounds.f $(INCLUDES)
jacu.o:		jacu.f $(INCLUDES)
l2norm.o:	l2norm.f $(INCLUDES)
pintgr.o:	pintgr.f $(INCLUDES)
read_input.o:	read_input.f $(INCLUDES)
rhs.o:		rhs.f $(INCLUDES)
setbv.o:	setbv.f $(INCLUDES)
setiv.o:	setiv.f $(INCLUDES)
setcoeff.o:	setcoeff.f $(INCLUDES)
ssor.o:		ssor.f $(INCLUDES)
domain.o:	domain.f $(INCLUDES)
verify.o:	verify.f $(INCLUDES)

clean:
	- /bin/rm -f npbparams.h
	- /bin/rm -f *.o *.d
