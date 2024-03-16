OBJECTS = Utils.o AtomsProp.o HydrogenBond.o CLProp.o QL.o Jacobi.o SList.o CalculTypesAmber.o CreateDefaultPDBTpl.o LoadPDBTemplate.o PDBTemplate.o SavePDBTemplate.o


include ../../CONFIG
CFLAGS =  $(COMMONCFLAGS) 

all: compiling 

%.o : %.c
	$(CC) -c $(CFLAGS) $< -o $@

compiling: $(OBJECTS)

clean:
	@$(RM) *o *.out *.log *.inp *~ 

#make dependance
SRCS = $(OBJECTS:%.o=%.c)
dep :
	gcc  $(CFLAGS) -MM $(SRCS) > Dep.mk

-include Dep.mk

