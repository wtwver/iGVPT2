OBJECTS = cansmidescriptor.o cmpdfilter.o filters.o groupcontrib.o inchidescriptor.o smartsdescriptors.o

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

