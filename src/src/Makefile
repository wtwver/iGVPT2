OBJECTS = OBGVPT2.o iGVPT2.o JobiGVPT2.o N2P2GVPT2.o TMGVPT2.o 

include $(OBLIBDIR)/CONFIG
include $(LIBCCHEMIDIR)/CONFIG

D0     = $(OBLIBDIR)

CFLAGS =  $(COMMONCFLAGS) 

all: compiling 

%.o : %.cpp
	$(CCPP) $(N2P2CFLAGS) -I $(LIBCCHEMIDIR)/src -c $(CFLAGS) $< -o $@

compiling: $(OBJECTS)

clean:
	@$(RM) *.o *~ 

#make dependance
SRCS = $(OBJECTS:%.o=%.cpp)
dep :
	gcc  $(CFLAGS) -MM $(SRCS) > Dep.mk

-include Dep.mk

