OBJECTS = addfilename.o addinindex.o addpolarh.o canonical.o fillUC.o forcefield.o gen2D.o gen3d.o largest.o loader.o ophighlight.o opisomorph.o optransform.o partialcharges.o readconformers.o sort.o xout.o

include ../CONFIG

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

