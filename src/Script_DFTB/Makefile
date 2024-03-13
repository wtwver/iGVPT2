# Modifier le contenu de ces variables pour adapter ce Makefile

#RM     = rm
#RM     = del

# You must define the LD_LIBRARY_PATH : 
# export LD_LIBRARY_PATH=/home/theochem/allouche/Softwares/boost/lib:$LD_LIBRARY_PATH
# You can also use xcomp

#USEMPI=yes
INC  = -fopenmp
LIB  = -lgomp
CC     =  g++  $(INC)
EXEC     = xRun


SRCS     =  Utils.cpp Matrix.cpp xRun.cpp 
CXXFLAGS = -Wall


# Ne rien changer a partir d'ici 

DEPEND_FILE = .depend
OBJS        = $(SRCS:%.c=%.o)

all: clean depend $(EXEC)

$(EXEC): $(OBJS)
	$(CC) $(OBJS) -o $(EXEC) -lm  $(LIB)

clean:
	$(RM) $(EXEC) *~ *.o \#*

cleanall: clean
	$(RM) $(EXEC) $(DEPEND_FILE)

cleaner: cleanall

depend: $(SRCS)
	g++ -MM $(SRCS) > $(DEPEND_FILE)

-include $(DEPEND_FILE)
