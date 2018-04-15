SRCFILES	= MPI_functions.cpp \
                  Stor.cpp \
                  GroupFct.cpp \
                  MPQ_Matrix.cpp \
                  DatabaseFct.cpp \
                  PolytopeEquiStab.cpp \
                  PolytopeFct.cpp \
                  MPQ_CDD_LRS.cpp \
                  Heuristic_fct.cpp \
                  DualDescriptionFct.cpp \
                  BankingSystem.cpp \
                  LatticeIsom.cpp \
                  MatrixOperation.cpp \
                  EnumerationConfig.cpp \
                  IterationAlgorithm.cpp \
                  VoronoiRankine.cpp \
                  Skeletton.cpp \
                  Get_md5sum.cpp

PROGRAM		= StandaloneGroupGraph \
                  StandaloneExchGroupPolytope \
                  StandaloneGroupPolytope \
                  StandaloneInitialVertex \
                  StandaloneInverseMatrix \
                  StandaloneFacPolytope \
                  StandaloneWMatIsom \
                  StandaloneFacPolytope_MPI \
                  StandaloneCyclicPolytope \
                  StandaloneOrbitFacePolytopeLP \
                  StandaloneOrbitFacePolytopeLPtree \
                  StandaloneOrbitFacePolytopeLPfast


PROGRAM_CDDDOU	= StandaloneVoronoiRankine

OBJECTS         = $(SRCFILES:.cpp=.o) \
                  $(PROGRAM:%=%.o) \
                  $(PROGRAM_CDDDOU:%=%.o)

LINKOBJECTS     = $(SRCFILES:.cpp=.o)

CC		= mpic++ -Wall -g 

NOLINKFLAG 	= -c
OUTPUTFLAG 	= -o


CFLAGS          = -Wall -ansi -DLRS -O3 -I$(LIBLRSGMP_INCDIR) -I$(LIBBLISS_INCDIR) -I$(LIBCDDGMP_INCDIR) -DGMPRATIONAL -DGMP -I$(PERMLIB_INCDIR) -I$(LIBSHORT_INCDIR)

LIBS            = -lm -lgmp -L$(LIBLRSGMP_LIBDIR) -llrs -L$(LIBBLISS_LIBDIR) $(LIBCDDGMP_LIBDIR)/libcddgmp.a  -lbliss -L$(LIBSHORT_LIBDIR) -lshort

LIBS_CDDDOU     = -lm -lgmp -L$(LIBLRSGMP_LIBDIR) -llrs -L$(LIBBLISS_LIBDIR) $(LIBCDDGMP_LIBDIR)/libcdd.a  -lbliss -L$(LIBSHORT_LIBDIR) -lshort


%.o: %.cpp
	@echo
	@echo compile $<
	$(CC) $(OUTPUTFLAG) $@ $(CFLAGS) $(NOLINKFLAG) $<

all:    $(PROGRAM) $(PROGRAM_CDDDOU)

$(PROGRAM_CDDDOU):     $(OBJECTS)
		@echo
		@echo link $@
		$(CC) $(LINKOBJECTS) $(DEBUGFLAG) $@.o $(LIBS_CDDDOU) -o $@

$(PROGRAM):     $(OBJECTS)
		@echo
		@echo link $@
		$(CC) -DGMP $(LINKOBJECTS) $(DEBUGFLAG) $@.o $(LIBS) -o $@

clean:;
	rm -f *.o
	rm -f $(PROGRAM)

#
# more dependencies
# =================
#

EnumerationConfig.o: 
MatrixOperation.o: 
StandaloneEnumConfig.o: EnumerationConfig.h 
