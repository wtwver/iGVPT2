include CONFIG

JOBDIR     = src/JobControl
EMDIR     = src/EmpriricalCorrections
MMDIR     = src/MolecularMechanics
QMDIR     = src/QuantumMechanics
PIMDDIR   = src/PathIntegral
MOLDIR    = src/Molecule
UTILDIR   = src/Utils
QFFDIR    = src/QuarticForceField
VPT2DIR   = src/VPT2
QFFPOTDIR   = src/QFFPot
IN2P2DIR   =  src/InterfaceLibN2P2
TMOLDIR   =  src/InterfaceTensorMol

JOBOBJ    =  $(JOBDIR)/*.o
MMOBJ     =  $(MMDIR)/*.o
EMOBJ     =  $(EMDIR)/*.o
QMOBJ     =  $(QMDIR)/*.o
PIMDOBJ   =  $(PIMDDIR)/*.o
MOLOBJ    =  $(MOLDIR)/*.o
UTILOBJ   =  $(UTILDIR)/*.o
QFFOBJ    =  $(QFFDIR)/*.o
VPT2OBJ   =  $(VPT2DIR)/*.o
QFFPOTOBJ   =  $(QFFPOTDIR)/*.o
IN2P2OBJ   =  $(IN2P2DIR)/*.o
RMOLOBJ   =  $(TMOLDIR)/*.o
 
SUBDIRS = \
	$(MMDIR) $(EMDIR) $(QMDIR) $(PIMDDIR) $(MOLDIR) $(UTILDIR) $(JOBDIR) $(QFFDIR) $(VPT2DIR) $(QFFPOTDIR) $(IN2P2DIR) $(TMOLDIR)

OBJECTS = \
	$(MMOBJ) $(EMOBJ) $(QMOBJ) $(PIMDOBJ) $(MOLOBJ) $(UTILOBJ) $(JOBOBJ) $(QFFOBJ) $(VPT2OBJ) $(QFFPOTOBJ) $(IN2P2OBJ) $(RMOLOBJ)

TARGETS = all clean dep
 
.PHONY: $(SUBDIRS)

$(TARGETS): $(SUBDIRS)

$(SUBDIRS):
	@$(MAKE) --directory=$@ $(MAKECMDGOALS)

all: cchemilib

cchemilib: $(SUBDIRS)
	$(CC) -fPIC -shared $(LDFLAGS) $(OBJECTS) -o lib/libcchemi.so  $(WIN32LIB)  $(OMPLIB) $(CLLIB) -g -lm
 
clean:
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
