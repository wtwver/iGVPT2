OBJECTS = alias.o atom.o base.o bitvec.o bond.o bondtyper.o builder.o canon.o chains.o chargemodel.o chiral.o data.o descriptor.o doxygen_pages.o fingerprint.o forcefield.o format.o generic.o graphsym.o grid.o griddata.o isomorphism.o kekulize.o locale.o matrix.o mcdlutil.o molchrg.o mol.o obconversion.o oberror.o obiter.o obmolecformat.o obutil.o op.o parsmart.o patty.o phmodel.o pointgroup.o query.o rand.o residue.o ring.o rotamer.o rotor.o spectrophore.o tokenst.o transform.o typer.o plugin.o

include CONFIG

CFLAGS =  $(COMMONCFLAGS) 

all: compiling 

%.o : %.cpp
	$(CCPP) -c $(CFLAGS) $< -o $@

compiling: $(OBJECTS)

clean:
	@$(RM) *.o *~ 

#make dependance
SRCS = $(OBJECTS:%.o=%.cpp)
dep :
	gcc  $(CFLAGS) -MM $(SRCS) > Dep.mk

-include Dep.mk

