OBJECTS = InterfaceCChemI.o InterfaceCChemIES.o

CFLAGS =  $(COMMONCFLAGS) -DHIGH_DERIVATIVES 

all: compiling 

%.o : %.cpp
	g++ -c -fPIC $(CFLAGS) -I include -L lib -l libnnp.a $< -o $@ 

compiling: $(OBJECTS)

clean:
	@$(RM) *o *.out *.log *.inp *~ 

#make dependance
SRCS = $(OBJECTS:%.o=%.c)
dep :
	g++  $(CFLAGS) -MM $(SRCS) > Dep.mk

-include Dep.mk

