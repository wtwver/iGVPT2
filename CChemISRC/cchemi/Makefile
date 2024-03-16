#pour compiler, taper Taper ./compcchemi.sh

include $(LIBCCHEMIDIR)/CONFIG

MAINDIR   = src/Main
JOBDIR     = $(LIBCCHEMIDIR)/src/JobControl
EMDIR     = $(LIBCCHEMIDIR)/src/EmpriricalCorrections
MMDIR     = $(LIBCCHEMIDIR)/src/MolecularMechanics
QMDIR     = $(LIBCCHEMIDIR)/src/QuantumMechanics
PIMDDIR     = $(LIBCCHEMIDIR)/src/PathIntegral
MOLDIR    = $(LIBCCHEMIDIR)/src/Molecule
UTILDIR   = $(LIBCCHEMIDIR)/src/Utils
QFFDIR    = $(LIBCCHEMIDIR)/src/QuarticForceField
VPT2DIR    = $(LIBCCHEMIDIR)/src/VPT2
IHDNNDIR  = $(LIBCCHEMIDIR)/src/InterfaceLibN2P2
TMOLDIR  = $(LIBCCHEMIDIR)/src/InterfaceTensorMol

LIBNNPDIR  = $(LIBCCHEMIDIR)/src/InterfaceLibN2P2/lib
LIBNNPLIB  = -lnnp -ldftd3 -lgfortran


MAINOBJ   = $(MAINDIR)/*.o
JOBOBJ     =  $(JOBDIR)/*.o
EMOBJ     =  $(EMDIR)/*.o
MMOBJ     =  $(MMDIR)/*.o
QMOBJ     =  $(QMDIR)/*.o
PIMDOBJ   =  $(PIMDDIR)/*.o
MOLOBJ    =  $(MOLDIR)/*.o
UTILOBJ   =  $(UTILDIR)/*.o
QFFOBJ   =  $(QFFDIR)/*.o
VPT2OBJ   =  $(VPT2DIR)/*.o
IHDNNOBJ  =  $(IHDNNDIR)/*.o
TMOLOBJ  =  $(TMOLDIR)/*.o
 
SUBDIRS = \
	$(MAINDIR) $(MMDIR) $(EMDIR) $(QMDIR) $(PIMDDIR) $(MOLDIR) $(UTILDIR) $(QFFDIR) $(JOBDIR) $(VPT2DIR) $(IHDNNDIR) $(TMOLDIR)

OBJECTS = \
	$(MAINOBJ) $(MMOBJ) $(EMOBJ) $(QMOBJ) $(PIMDOBJ) $(MOLOBJ) $(UTILOBJ) $(QFFOBJ) $(JOBOBJ) $(VPT2OBJ) $(IHDNNOBJ) $(TMOLOBJ)

TARGETS = all clean dep
 
.PHONY: $(SUBDIRS)

$(TARGETS): $(SUBDIRS)

$(SUBDIRS):
	@$(MAKE) --directory=$@ $(MAKECMDGOALS)

all: cchemi

cchemi: $(SUBDIRS)
	g++ $(LDFLAGS) $(OBJECTS) -o cchemi $(WIN32LIB)  $(OMPLIB) $(CLLIB) -g -lm -L $(LIBNNPDIR) $(LIBNNPLIB)  -lgomp $(PYTHONLIB) 
#	$(CC) $(LDFLAGS) $(OBJECTS) -o cchemi $(WIN32LIB)  $(OMPLIB) $(CLLIB) -g -lm
 
clean:
	@$(RM) cchemi
	@$(RM) cchemi.exe
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
