#pour compiler, taper Taper ./compigvpt2.sh

include OpenBabel232/lib/CONFIG

D0     = $(OBLIBDIR)
D1     = $(OBLIBDIR)/charges
D2     = $(OBLIBDIR)/ops
D3     = $(OBLIBDIR)/descriptors
D4     = $(OBLIBDIR)/fingerprints
D5     = $(OBLIBDIR)/forcefields
D6     = $(OBLIBDIR)/formats
D7     = $(OBLIBDIR)/math
D9     = $(OBLIBDIR)/stereo

O0     = $(OBLIBDIR)/*.o
Oa     = $(OBLIBDIR)/*/*.o

TDIR   = $(OBLIBDIR)/../../src

OBGVPT2     = $(TDIR)/OBGVPT2.o
IGVPT2      = $(TDIR)/iGVPT2.o
JOBIGVPT2   = $(TDIR)/JobiGVPT2.o
N2P2GVPT2      = $(TDIR)/N2P2GVPT2.o
TMGVPT2      = $(TDIR)/TMGVPT2.o


include $(LIBCCHEMIDIR)/CONFIG
MOLDIR    = $(LIBCCHEMIDIR)/src/Molecule
UTILDIR   = $(LIBCCHEMIDIR)/src/Utils
QFFDIR    = $(LIBCCHEMIDIR)/src/QuarticForceField
VPT2DIR    = $(LIBCCHEMIDIR)/src/VPT2
JOBDIR     = $(LIBCCHEMIDIR)/src/JobControl
EMDIR     = $(LIBCCHEMIDIR)/src/EmpriricalCorrections
MMDIR     = $(LIBCCHEMIDIR)/src/MolecularMechanics
QMDIR     = $(LIBCCHEMIDIR)/src/QuantumMechanics
PIMDDIR     = $(LIBCCHEMIDIR)/src/PathIntegral
IN2P2DIR  = $(LIBCCHEMIDIR)/src/InterfaceLibN2P2
TMOLDIR  = $(LIBCCHEMIDIR)/src/InterfaceTensorMol

 

LIBN2P2DIR  = $(LIBCCHEMIDIR)/src/InterfaceLibN2P2/lib
LIBN2P2LIB  = -lnnp -ldftd3 -lgfortran

MOLOBJ    =  $(MOLDIR)/*.o
UTILOBJ   =  $(UTILDIR)/*.o
QFFOBJ   =  $(QFFDIR)/*.o
VPT2OBJ   =  $(VPT2DIR)/*.o

JOBOBJ     =  $(JOBDIR)/*.o
EMOBJ     =  $(EMDIR)/*.o
MMOBJ     =  $(MMDIR)/*.o
QMOBJ     =  $(QMDIR)/*.o
PIMDOBJ   =  $(PIMDDIR)/*.o
IN2P2OBJ  =  $(IN2P2DIR)/*.o
TMOLOBJ  =  $(TMOLDIR)/*.o

SUBDIRS = \
	$(D0) $(D1) $(D2) $(D3) $(D4) $(D5) $(D6) $(D7) $(D9) $(TDIR) \
        $(MOLDIR) $(UTILDIR) $(QFFDIR) $(VPT2DIR)\
	$(MMDIR) $(EMDIR) $(QMDIR) $(JOBDIR) $(PIMDDIR) $(IN2P2DIR) $(TMOLDIR)


OBJECTSOB = \
	$(O0) $(Oa) \

OBJECTSCC = \
        $(MOLOBJ) $(UTILOBJ) $(QFFOBJ) $(VPT2OBJ) \
	$(MMOBJ) $(EMOBJ) $(QMOBJ) $(JOBOBJ) $(PIMDOBJ) $(IN2P2OBJ) $(TMOLOBJ)

OBJECTS = \
	$(O0) $(Oa) \
        $(MOLOBJ) $(UTILOBJ) $(QFFOBJ) $(VPT2OBJ) \
	$(MMOBJ) $(EMOBJ) $(QMOBJ) $(JOBOBJ) $(PIMDOBJ)  $(IN2P2OBJ) $(TMOLOBJ)

TARGETS = all clean dep
 
.PHONY: $(SUBDIRS)

$(TARGETS): $(SUBDIRS)

$(SUBDIRS):
	@$(MAKE) --directory=$@ $(MAKECMDGOALS)

all: $(SUBDIRS)
	ar rcs libcchemi.a $(OBJECTSCC) 
	ar rcs libob.a $(OBJECTSOB) 
	c++ -O -fPIC -o igvpt2 $(IGVPT2) $(JOBIGVPT2) $(IGVPT2UTILS) $(OBGVPT2) $(N2P2GVPT2)  $(WIN32LIB)  $(OMPLIB) $(CLLIB) $(TMGVPT2) -g -lm -L./ -lcchemi -L $(LIBN2P2DIR) $(LIBN2P2LIB) -lob -lcrypto -lgomp $(PYTHONLIB)

clean:
	@$(RM) igvpt2
	@$(RM) *.exe
	@$(RM) p
	@$(RM) pp
	@$(RM) p.c
	@$(RM) *.log
	@$(RM) *.out
	@$(RM) *.aux
	@$(RM) *.cmd
	@$(RM) *.trj
	@$(RM) *.txt
	@$(RM) geom*
	@$(RM) traj*
	@$(RM) *.arc
	@$(RM) PUNCH
	@$(RM) prop*
	@$(RM) res*
	@$(RM) *.a
