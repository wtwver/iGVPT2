OBJECTS = cchemi.o


include $(LIBCCHEMIDIR)/CONFIG
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

