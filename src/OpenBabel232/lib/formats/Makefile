OBJECTS = abinitformat.o acesformat.o acrformat.o adfformat.o alchemyformat.o amberformat.o APIInterface.o balstformat.o bgfformat.o boxformat.o cacaoformat.o cacheformat.o carformat.o castepformat.o cccformat.o chem3dformat.o chemdrawcdx.o chemdrawct.o chemkinformat.o chemtoolformat.o cifformat.o copyformat.o crkformat.o crystal09format.o CSRformat.o cssrformat.o dlpolyformat.o dmolformat.o exampleformat.o fastaformat.o fastsearchformat.o fchkformat.o featformat.o fhformat.o fhiaimsformat.o fingerprintformat.o fpsformat.o freefracformat.o gamessformat.o gamessukformat.o gausscubeformat.o gaussformat.o gausszmatformat.o genbankformat.o ghemicalformat.o groformat.o gromos96format.o gulpformat.o hinformat.o jaguarformat.o lmpdatformat.o MCDLformat.o mdlformat.o mmcifformat.o mmodformat.o MNAformat.o mol2format.o moldenformat.o molproformat.o molreport.o mopacformat.o mpdformat.o mpqcformat.o msiformat.o msmsformat.o nulformat.o nwchemformat.o opendxformat.o outformat.o pcmodelformat.o pdbformat.o pdbqtformat.o posformat.o povrayformat.o pqrformat.o PQSformat.o pwscfformat.o qchemformat.o reportformat.o rsmiformat.o rxnformat.o shelxformat.o smilesformat.o textformat.o thermoformat.o tinkerformat.o titleformat.o turbomoleformat.o unichemformat.o vaspformat.o viewmolformat.o xedformat.o xsfformat.o xtcformat.o xyzformat.o yasaraformat.o zindoformat.o

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

