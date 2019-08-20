#!-------------------------------------------------------------------------!
#!									   !
#!	  N  A  S     P A R A L L E L	  B E N C H M A R K S  3.0	   !
#!									   !
#!			 J A V A	 V E R S I O N  		   !
#!									   !
#!				 M a k e f i l e			   !
#!									   !
#!-------------------------------------------------------------------------!
#!									   !
#!    Permission to use, copy, distribute and modify this software	   !
#!    for any purpose with or without fee is hereby granted.  We	   !
#!    request, however, that all derived work reference the NAS 	   !
#!    Parallel Benchmarks 3.0. This software is provided "as is"	   !
#!    without express or implied warranty.				   !
#!									   !
#!    Information on NPB 3.0, including the Technical Report NAS-02-008	   !
#!    "Implementation of the NAS Parallel Benchmarks in Java",		   !
#!    original specifications, source code, results and information	   !
#!    on how to submit new results, is available at:			   !
#!									   !
#!	     http://www.nas.nasa.gov/Software/NPB/			   !
#!									   !
#!    Send comments or suggestions to  npb@nas.nasa.gov 		   !
#!									   !
#!	    NAS Parallel Benchmarks Group				   !
#!	    NASA Ames Research Center					   !
#!	    Mail Stop: T27A-1						   !
#!	    Moffett Field, CA	94035-1000				   !
#!									   !
#!	    E-mail:  npb@nas.nasa.gov					   !
#!	    Fax:     (650) 604-3957					   !
#!									   !
#!-------------------------------------------------------------------------!
#!     Translation to Java and to MultiThreaded Code:			   !
#      M. Frumkin							   !
#!     M. Schultz							   !
#!-------------------------------------------------------------------------!
JAVAC = $(JROOT)/bin/javac

bt: BT.class \
	BTThreads/BTBase.class BTThreads/RHSCompute.class  \
	BTThreads/YSolver.class BTThreads/RHSAdder.class  \
	BTThreads/XSolver.class BTThreads/ZSolver.class \
	BMInOut/BMArgs.class BMInOut/BMResults.class \
	Timer.class
	
cg: CG.class \
	CGThreads/CGBase.class CGThreads/CGWorker.class \
	BMInOut/BMArgs.class BMInOut/BMResults.class \
	Random.class  Timer.class

ft: FT.class \
	FTThreads/FTBase.class FTThreads/FFTThread.class \
	FTThreads/EvolveThread.class \
	BMInOut/BMArgs.class BMInOut/BMResults.class \
	Random.class  Timer.class
		
is: IS.class \
	ISThreads/ISBase.class ISThreads/RankThread.class \
	BMInOut/BMArgs.class BMInOut/BMResults.class \
	Random.class  Timer.class
		
lu: LU.class \
	LUThreads/LUBase.class LUThreads/Adder.class \
	LUThreads/RHSCompute.class LUThreads/UpperJac.class \
	LUThreads/LowerJac.class LUThreads/Scale.class \
	BMInOut/BMArgs.class BMInOut/BMResults.class \
	Timer.class

mg: MG.class \
	MGThreads/MGBase.class MGThreads/Interp.class \
	MGThreads/Psinv.class MGThreads/Resid.class MGThreads/Rprj.class \
	BMInOut/BMArgs.class BMInOut/BMResults.class \
	Random.class  Timer.class
		
sp: SP.class \
	SPThreads/SPBase.class SPThreads/RHSAdder.class \
	SPThreads/RHSCompute.class SPThreads/TXInverse.class \
	SPThreads/XSolver.class SPThreads/YSolver.class \
	SPThreads/ZSolver.class \
	BMInOut/BMArgs.class BMInOut/BMResults.class \
	Timer.class
		
all: bt cg ft is mg lu sp
		
BT.class:
	$(JAVAC) $(JAVAFLAGS) BT.java
	
CG.class:
	$(JAVAC) $(JAVAFLAGS) CG.java

FT.class:
	$(JAVAC) $(JAVAFLAGS) FT.java

IS.class:
	$(JAVAC) $(JAVAFLAGS) IS.java

LU.class:
	$(JAVAC) $(JAVAFLAGS) LU.java

MG.class:
	$(JAVAC) $(JAVAFLAGS) MG.java

SP.class:
	$(JAVAC) $(JAVAFLAGS) SP.java

Random.class:
	$(JAVAC) $(JAVAFLAGS) Random.java

Timer.class:
	$(JAVAC) $(JAVAFLAGS) Timer.java

.SUFFIXES : .java.class
.java.class :
	$(JAVAC) $(JAVAFLAGS) $<
	
#delete all class files
clean:
		- rm *.class
		- rm BTThreads/*.class
		- rm SPThreads/*.class
		- rm LUThreads/*.class
		- rm FTThreads/*.class
		- rm ISThreads/*.class
		- rm MGThreads/*.class
		- rm CGThreads/*.class
		- rm BMInOut/*.class
