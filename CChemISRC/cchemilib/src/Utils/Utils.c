/********************************************************************************
 cchemi is an interface to ab initio computational chemistry programs 
 designed for add them many functionalities non available in these packages.
 Copyright (C) 2010 Abdulrahman Allouche (University Lyon 1)

 This program is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.

 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with this program.  If not, see <https://www.gnu.org/licenses/>.
********************************************************************************/

/* Utils.c */
#include <stdlib.h>
#include <ctype.h>
#include <math.h>
#include <time.h>
#include <stdio.h>
#include <string.h>
#include <stdarg.h> 
#include <assert.h> 
#include <sys/stat.h>
#include <sys/types.h>
#include <fcntl.h>
#ifdef ENABLE_MPI
#include <mpi.h>
#endif

#ifdef OS_WIN32
#include <windows.h>
#include <io.h>
#include <direct.h>
#include <io.h>
#define mkdir(p,m) _mkdir(p)
#ifndef S_ISDIR
#define S_ISDIR(mode) ((mode)&_S_IFDIR)
#endif
#else /* OS_WIN32 */
#include <pwd.h>
#include <unistd.h> 
#include <sys/times.h>
#endif /* OS_WIN32 */

#include "../Utils/Types.h"
#include "../Utils/Constants.h"
#include "../Utils/Timer.h"
#include "../Utils/Utils.h"
#include "../Utils/AtomsProp.h"
#include "../Utils/HydrogenBond.h"
#include "../MolecularMechanics/MolecularMechanics.h"
#ifdef ENABLE_CL
#include "../Utils/CLProp.h"
#endif

/********************************************************************************/
char* strdup_vprintf(const char* format, va_list ap)
{
	va_list ap2;
	int size;
	char* buffer;

	va_copy(ap2, ap);
	size = vsnprintf(NULL, 0, format, ap2)+1;
	va_end(ap2);

	buffer = malloc(size+1);
	assert(buffer != NULL);

	vsnprintf(buffer, size, format, ap);
	return buffer;
}
/********************************************************************************/
char* strdup_printf(const char* format, ...)
{
	char* buffer;
	va_list ap;
	va_start(ap, format);
	buffer = strdup_vprintf(format, ap);
	va_end(ap);
	return buffer;
}
/********************************************************************************/
#ifndef OS_WIN32
#define TIMER_TICK      60
static clock_t it;
static struct tms itt;
void timing(double* cpu,double *sys)
{
	it=times(&itt);
	*cpu=(double) itt.tms_utime / (double) TIMER_TICK;
	*sys=(double) itt.tms_stime / (double) TIMER_TICK;
}
#endif
#ifdef OS_WIN32
void addUnitDisk(FILE* file, const char* name)
{
	if(name && strlen(name)>1 && name[1]==':')
		fprintf(file,"%c%c\n", name[0],name[1]);
}
#endif
/********************************************************************************/
char* getTimeStr()
{
	char* str=NULL;
	time_t t;
	struct tm* ts;

	t = time(NULL);
	ts = localtime(&t);
	str = asctime (ts);

	return str;
}
/********************************************************************************/
boolean isABackspace(char *st)
{
        int i;
        for(i=0;i<(int)strlen(st);i++)
        	if(st[i] != ' ' && st[i] !='\n' && st[i] !='\r')
                	return FALSE;
        return TRUE;
}
/********************************************************************************/
void waiting(double tsecond)
{
        TimerType timer;
        double elaps;

        timer_init(timer);
        do{
        	timer_start( timer );
		elaps = exp(1.0);
        	timer_stop(timer);
                elaps = timer_get(timer);
        }while(elaps*1e-6<tsecond);

}
/*************************************************************************************/
void debug(char *fmt,...)
{
#ifdef DEBUG
	va_list ap;
	va_start(ap,fmt);
	vfprintf(stdout, fmt, ap);
	va_end(ap);
#endif

}
/********************************************************************************/
char* getLineChars(char c,int n)
{
	int i;
	char *line = NULL;

	if(n<1) return line;
	line = malloc((n+1)*sizeof(char));
	for(i=0;i<n;i++) line[i] = c;
	line[n] = '\0';

	return line;
	
}
/********************************************************************************/
char* catFile(char* namefile,boolean tabulation)
{
	char *t = NULL;
	char *tsrt = NULL;
	FILE *fd;
	char *dump = NULL;


	t=malloc(BSIZE*sizeof(char));

	fd = fopen(namefile, "r");
	if(fd)
	{
		while(!feof(fd))
  		{
    		if(!fgets(t,BSIZE, fd)) break;
                dump = tsrt;
		if(!tsrt)
		{
			if(tabulation)
				tsrt = strdup_printf("\t%s",t);
			else
				tsrt = strdup_printf("%s",t);
		}
		else
		{
			if(tabulation)
				tsrt = strdup_printf("%s\t%s",tsrt,t);
			else
				tsrt = strdup_printf("%s%s",tsrt,t);
			free(dump);
			dump = NULL;
		}
  		}
 		fclose(fd);
		unlink (namefile);
 	}
	else
	{
		tsrt = NULL;
	}
	free(t);
	t = tsrt;

	return tsrt;
}
/*************************************************************************************/
char *runCommand(char *command)
{
	char *t;
	char *terr = NULL;
	FILE *fd;
	char *temp;
	char *outfile= strdup_printf("%s%stmp%soutfile",cchemiDirectory(), DIR_SEPARATOR_S, DIR_SEPARATOR_S);
	char *errfile= strdup_printf("%s%stmp%serrfile",cchemiDirectory(), DIR_SEPARATOR_S, DIR_SEPARATOR_S);
	char *dump;

	temp = strdup_printf("sh -c '%s >%s 2>%s'",command,outfile,errfile);
	system(temp);

	t=malloc(BSIZE*sizeof(char));

	fd = fopen(errfile, "r");
	if(fd)
	{
  		while(!feof(fd))
  		{
    			if(!fgets(t,BSIZE, fd)) break;
                	dump = terr;
			if(!terr) terr = strdup_printf("%s",t);
			else
			{
				terr = strdup_printf("%s%s",terr,t);
				free(dump);
			}
  		}
 		fclose(fd);
		unlink (errfile);
 	}
 	else terr = NULL;

	fd = fopen(outfile, "r");
	if(fd)
	{
		unlink (outfile);
	}

	free(t);
	free(temp);
	free(outfile);
	free(errfile);

	return terr;
}
/********************************************************************************/
char *getSuffixNameFile(const char* allname)
{
	char* name = strdup(allname);
	int len=strlen(allname);
	int i;
	for(i=len;i>0;i--)
	if(name[i]=='.')
	{
		name[i] = '\0';
		break;
	}
	if(!strstr(name,DIR_SEPARATOR_S))
	{
		char*t=strdup_printf("%s%s%s", getenv("PWD"),DIR_SEPARATOR_S,name);
		free(name);
		name = t;
	}
	return name;
}
/********************************************************************************/
char* getHomeDir()
{
	char *ptr;
	char *copy = NULL;

	ptr = getenv ("HOME");

#ifdef OS_WIN32
	if (ptr == NULL) return strdup("C:\\");
#else
	if (ptr == NULL) return strdup(DIR_SEPARATOR_S);
#endif

/*
	copy = (char *) malloc (strlen (ptr) + 2);
	strcpy (copy, ptr);
	strcat (copy, DIR_SEPARATOR_S);
*/
	copy = strdup_printf("%s%s",ptr,DIR_SEPARATOR_S);

	return copy;
}
/********************************************************************************/
CONST char* cchemiDirectory()
{
	static char *cchemi_dir = NULL;
	char *home_dir;
	char* Version_S = NULL;

	if (cchemi_dir != NULL) return cchemi_dir;

	home_dir = strdup(getHomeDir());

#ifdef OS_WIN32
	Version_S = strdup_printf("%d%d%d",MAJOR_VERSION,MINOR_VERSION,MICRO_VERSION);
	cchemi_dir = strdup_printf("%s%s%s",home_dir,"cchemi",Version_S);
#else
	Version_S = strdup_printf("%d.%d.%d",MAJOR_VERSION,MINOR_VERSION,MICRO_VERSION);
	cchemi_dir = strdup_printf("%s%s%s",home_dir,".cchemi-",Version_S);
#endif

	free(Version_S);
	free(home_dir);
	return cchemi_dir;
}
/*************************************************************************************/
void deleteLastChar(char *str)
{
        str[strlen(str)-1]='\0';
}
/*************************************************************************************/
void createRessourceFile()
{
	saveAtomsProp();
}
/*************************************************************************************/
void addToPath(char* toAdd)
{
#ifdef OS_WIN32
	{
		char t[BSIZE];
		sprintf(t,"%s;%cPATH%c",
		toAdd,
		'%','%');
		if(strlen(t)>1) setenv("PATH",t,TRUE);
	}
#else
	{
		char t[BSIZE];
		sprintf(t,"%s:%cPATH",
		toAdd,
		'$');
		if(strlen(t)>1) setenv("PATH",t,TRUE);
	}
#endif
}
/*************************************************************************************/
void readRessourceFile()
{
	boolean rOK = FALSE;
 
	defineDefaultAtomsProp();
	rOK = readAtomsProp();
	if(!rOK) defineDefaultAtomsProp();
}
/*************************************************************************************/
void uppercase(char *str)
{
	while( *str != '\0')
	{
		if (isalpha((int)*str))
		if (islower((int)*str))
			*str = toupper((int)*str);
		str ++;
	}
}
/*************************************************************************************/
void lowercase(char *str)
{
	while( *str != '\0')
	{
		*str = (char)tolower((int)*str);
	str ++;
	}
}
/*************************************************************************************/
void strfreev (char **str)
{
	int i;
	if (!str) return;

	for (i = 0; str[i] != NULL; i++) free (str[i]);

	free (str);
}
/*************************************************************************************/
char** split(char *str)
{
	char** strsplit= malloc(sizeof(char*));
	int n=0;
	char* t=str;
	char p[BSIZE];
	while(*t!='\n' && *t !='\0')
	{
		if(*t!=' ' && *t !='\t')
		{
			n++;
			strsplit= realloc(strsplit,(n+1)*sizeof(char*));
			sscanf(t,"%s",p);
			strsplit[n-1]= strdup(p);
			while(*t!=' ' && *t !='\t')
			{
				t++;
				if(*t =='\n' || *t =='\0') break;
			}

		}
		else
		{
			while(*t ==' ' || *t == '\t' )
			{
				t++;
				if(*t =='\n' || *t =='\0') break;
			}
		}
	}
	strsplit[n]= NULL;
	return strsplit;
}
/********************************************************************************/
void deleteLastSpaces(char* str)
{
	char *s;

	if(str == NULL)
		return;

	if (!*str)
		return;
	for (s = str + strlen (str) - 1; s >= str && isspace ((unsigned char)*s); s--)
		*s = '\0';
}
/********************************************************************************/
void deleteFirstSpaces(char* str)
{
	char *start;
	int i;
	int lenSpace = 0;

	if(str == NULL)
		return;
	if (!*str)
		return;

	for (start = str; *start && isspace (*start); start++)lenSpace++;

	for(i=0;i<(int)(strlen(str)-lenSpace);i++)
		str[i] = str[i+lenSpace];
	str[strlen(str)-lenSpace] = '\0';
}
/********************************************************************************/
void deleteAllSpaces(char* str)
{
	int i;
	int j;
	boolean Ok = FALSE;

	deleteLastSpaces(str);
	deleteFirstSpaces(str);
	while(!Ok)
	{
		Ok = TRUE;
		for(i=0;i<(int)strlen(str);i++)
		{
			if(isspace(str[i]))
			{
				Ok = FALSE;
				for(j=i;j<(int)strlen(str);j++)
				{
					str[j] = str[j+1];
				}
				break;
			}
		}
	}
}
/**********************************************/
char* getToStr(char* str,char* end)
{
	char* iend = NULL;
	char* res = NULL;
	int len;
	int i;

	if(str == NULL || end == NULL)
		return NULL;

	iend = strstr(str,end);
	if(iend==NULL)
		return strdup(str);
	len = iend - str;
	if(len<1)
		return NULL;

	res = malloc((len+1)*sizeof(char));
	for(i=0;i<len;i++)
		res[i] = str[i];

	res[len] = '\0';
	return res;
	
}
/*************************************************************************************/
static boolean testi(char c)
{
	switch ( c )
	{
	case	'0':
	case	'1':
	case	'2':
	case	'3':
	case	'4':
	case	'5':
	case	'6':
	case	'7':
	case	'8':
	case	'9': return TRUE;
	}
	return FALSE;
}
/*************************************************************************************/
boolean isInteger(CONST char *t)
{
	int i;
	if(!testi(t[0])&& t[0] != '-' ) return FALSE;
	for(i=1;i<strlen(t);i++)
		if(!testi(t[i]) ) return FALSE;
	return TRUE;

}
/*************************************************************************************/
static boolean testascii(char c)
{
	switch ( c )
	{
	case	'0':
	case	'1':
	case	'2':
	case	'3':
	case	'4':
	case	'5':
	case	'6':
	case	'7':
	case	'8':
	case	'9':
	case	'.':
	case	'e':
	case	'E':
	case	'+':
	case	'-':return TRUE;
	}
	return FALSE;
}
/*************************************************************************************/
boolean isFloat(const char *t)
{
	int i;
	for(i=0;i<strlen(t);i++)
		if(!testascii(t[i]) ) return FALSE;
	if(t[0] =='e' || t[0] =='E' ) return FALSE;
	return TRUE;

}
/*************************************************************************************/
void strDeleten(char* str)
{
	char *s;

	if(str == NULL)
		return;

	if (!*str)
		return;
	for (s = str + strlen (str) - 1; s >= str && ((unsigned char)*s)=='\n'; s--)
		*s = '\0';
}
/*************************************************************************************/
char * mystrcasestr(const char *haystack, const char *needle)
{
	char *i, *startn = 0, *j = 0;
	for (i = (char*)haystack; *i; i++)
	{
		if(j)
		{
			if (toupper(*i) == toupper(*j))
			{
				if (!*++j)
				return startn;
			}
			else j = 0;
		}
		else if (toupper(*i) == toupper(*needle))
		{
			j = (char*)needle + 1;
			startn = i;
		}
	}
	return 0;
}
/*************************************************************************************/
boolean readOneReal(FILE* file, char* tag, double*value)
{
	static char *t = NULL; 
	char* TAG = NULL;
	char* pos;
	if(!tag) return FALSE;
	if(!value) return FALSE;
	if(t==NULL) t = malloc(BSIZE*sizeof(char));

	TAG = strdup(tag);
	uppercase(TAG);
	rewind(file);

	while(!feof(file))
  	{
    		if(!fgets(t,BSIZE, file)) break;
		deleteFirstSpaces(t);
		if(t[0]=='#') continue;
		uppercase(t);
		pos = strstr(t,TAG);
		if(!pos) continue;
		if(strstr(pos,"=")) 
		{
			pos = strstr(pos,"=")+1;
		}
		else pos += strlen(TAG)+1;
		free(TAG);
		if(1==sscanf(pos,"%lf",value)) return TRUE;
		return FALSE;
	}
	free(TAG);
	return FALSE;
}
/********************************************************************************/
boolean readOneRealFromAFile(char* namefile, char* tag, double* value)
{
	FILE* file = NULL;
	boolean res;

	if(!namefile) return FALSE;

	file = fopen(namefile, "rb");
	res = readOneReal(file,tag,value);
	fclose(file);
	return res;
}
/*************************************************************************************/
boolean readOneInt(FILE* file, char* tag, int*value)
{
	static char *t = NULL; 
	char* TAG = NULL;
	char* pos;
	if(!tag) return FALSE;
	if(!value) return FALSE;
	if(t==NULL) t = malloc(BSIZE*sizeof(char));

	TAG = strdup(tag);
	uppercase(TAG);
	rewind(file);

	while(!feof(file))
  	{
    		if(!fgets(t,BSIZE, file)) break;
		deleteFirstSpaces(t);
		if(t[0]=='#') continue;
		uppercase(t);
		pos = strstr(t,TAG);
		if(!pos) continue;
		if(strstr(pos,"=")) 
		{
			pos = strstr(pos,"=")+1;
		}
		else pos += strlen(TAG)+1;
		free(TAG);
		if(1==sscanf(pos,"%d",value)) return TRUE;
		return FALSE;
	}
	free(TAG);
	return FALSE;
}
/********************************************************************************/
boolean readOneIntFromAFile(char* namefile, char* tag, int* value)
{
	FILE* file = NULL;
	boolean res;

	if(!namefile) return FALSE;

	file = fopen(namefile, "rb");
	res = readOneInt(file,tag,value);
	fclose(file);
	return res;
}
/*************************************************************************************/
boolean readOneBoolean(FILE* file, char* tag, boolean*value)
{
	static char *t = NULL; 
	char* TAG = NULL;
	char* pos;
	char tmp[100];
	if(!tag) return FALSE;
	if(!value) return FALSE;
	if(t==NULL) t = malloc(BSIZE*sizeof(char));

	TAG = strdup(tag);
	uppercase(TAG);
	rewind(file);

	while(!feof(file))
  	{
    		if(!fgets(t,BSIZE, file)) break;
		deleteFirstSpaces(t);
		if(t[0]=='#') continue;
		uppercase(t);
		pos = strstr(t,TAG);
		if(!pos) continue;
		if(strstr(pos,"=")) 
		{
			pos = strstr(pos,"=")+1;
		}
		else pos += strlen(TAG)+1;
		free(TAG);
		if(1==sscanf(pos,"%s",tmp)) 
		{
			
			if(!strcmp(tmp,"TRUE"))*value = TRUE;
			else *value = FALSE;
			return TRUE;
		}
		return FALSE;
	}
	free(TAG);
	return FALSE;
}
/********************************************************************************/
boolean readOneBooleanFromAFile(char* namefile, char* tag, boolean* value)
{
	FILE* file = NULL;
	boolean res;

	if(!namefile) return FALSE;

	file = fopen(namefile, "rb");
	res = readOneBoolean(file,tag,value);
	fclose(file);
	return res;
}
/********************************************************************************/
boolean readOneStringFromAFile(char* namefile, char* tag, int* value)
{
	FILE* file = NULL;
	boolean res;

	if(!namefile) return FALSE;

	file = fopen(namefile, "rb");
	res = readOneInt(file,tag,value);
	fclose(file);
	return res;
}
/*************************************************************************************/
boolean readOneString(FILE* file, char* tag, char**value)
{
	static char *t = NULL; 
	static char *t2 = NULL; 
	char* TAG = NULL;
	char* pos;
	if(!tag) return FALSE;
	if(!value) return FALSE;
	if(t==NULL) t = malloc(BSIZE*sizeof(char));
	if(t2==NULL) t2 = malloc((BSIZE+2)*sizeof(char));

	TAG = strdup(tag);
	uppercase(TAG);
	rewind(file);

	while(!feof(file))
  	{
    		if(!fgets(t,BSIZE, file)) break;
		deleteFirstSpaces(t);
		if(t[0]=='#') continue;
		sprintf(t2,"%s",t);
		uppercase(t2);
		pos = strstr(t2,TAG);
		if(!pos) continue;
		if(strstr(pos,"=")) 
		{
			pos = strstr(pos,"=")+1;
		}
		else pos += strlen(TAG)+1;
		free(TAG);
		if(strlen(pos)>0) 
		{
			char* p = t+(int)(pos-t2);
			if(*value) free(*value);
			*value = strdup(p);
			strDeleten(*value);
			deleteFirstSpaces(*value);
			deleteLastSpaces(*value);
			return TRUE;
		}
		return FALSE;
	}
	free(TAG);
	return FALSE;
}
/*****************************************************************************/
void setMDOptions(FILE* file, int* updateFrequency, 
double* heatTime, double*equiTime, double* runTime, double* coolTime, 
double* heatTemp,  double*equiTemp, double*runTemp, double*coolTemp, 
double* stepSize, MDIntegratorType* integrator, MDThermostatType* thermostat, double* friction, double* omegaMax, int* Nf, double* collide, double* qNH)
{
	int itmp;
	*updateFrequency = 5;
	readOneInt(file,"updateFrequency",updateFrequency);
	if(*updateFrequency<0) *updateFrequency = 0;

	*heatTime = 1;
	readOneReal(file,"heatTime",heatTime);
	*equiTime = 2;
	readOneReal(file,"equiTime",equiTime);
	*runTime = 10;
	readOneReal(file,"runTime",runTime);
	*coolTime = 10;
	readOneReal(file,"coolTime",coolTime);
	if(*heatTime<0) *heatTime = 1;
	if(*equiTime<0) *equiTime = 1;
	if(*runTime<0) *runTime = 1;
	if(*coolTime<0) *coolTime = 1;

	*heatTemp = 0;
	readOneReal(file,"heatTemp",heatTemp);
	*runTemp = 300;
	readOneReal(file,"runTemp",runTemp);
	*equiTemp = *runTemp;
	*coolTemp = *heatTemp;
	if(*heatTemp<0) *heatTemp = 0;
	if(*equiTemp<0) *runTemp = 300;
	if(*runTemp<0) *runTemp = 300;
	if(*coolTemp<0) *coolTemp = 0;

	*stepSize = 0.5;
	readOneReal(file,"stepSize",stepSize);
	if(*stepSize<0) *stepSize = 1.0;
	if(*stepSize>5) *stepSize = 5.0;


	*integrator = BEEMAN;
	/* *integrator = STOCHASTIC;*/
	if(readOneInt(file,"integrator",&itmp))*integrator = itmp;

	/* *thermostat = ANDERSEN;*/
	*thermostat = BERENDSEN;
	if(readOneInt(file,"thermostat",&itmp))*thermostat = itmp;

	if(*integrator==STOCHASTIC) *thermostat = NONE;
	if(*integrator==QTB) *thermostat = NONE;


	*friction=-1;
	readOneReal(file,"friction",friction);
	*omegaMax=4000;
	readOneReal(file,"omegaMax",omegaMax);
	*Nf=50;
	readOneInt(file,"Nf",Nf);

	*collide = 20;
	readOneReal(file,"collide",collide);
	*qNH = 20;
	readOneReal(file,"qNoseHoover",qNH);
}
/********************************************************************************/
static boolean createUserCChemIDirectory()
{
	if (mkdir (cchemiDirectory(), 0755) < 0) 
	{
		printf(("Installation failed.  Contact system administrator."));
		exit(1);
		return FALSE;
	}
	return TRUE;
}
/*****************************************************************************************************************************************************/
void userInstallVerify()
{
	const char *filename;
	struct stat  stat_buf;

	filename = cchemiDirectory();

	if (stat(filename, &stat_buf) != 0)
	{
		createUserCChemIDirectory();
		saveAtomsProp();
		saveAmberParameters();
		saveHBondsProperties();
		/*system(strdup_printf("ls -l %s",cchemiDirectory()));*/
	}
}
/*****************************************************************************************************************************************************/
void readRessources()
{
	userInstallVerify();
	readAtomsProp();
	loadAmberParameters();
	readHBondsProperties();

}
/*****************************************************************************************************************************************************/
void initSeed(char* inputFileName)
{
	int seed = -1;
	readOneIntFromAFile(inputFileName,"SeedRandom",&seed);
	if(seed>=0) srand((unsigned int)seed);
	else srand((unsigned int)time(NULL));
	printf("seed = %d\n",seed);
}
/*****************************************************************************************************************************************************/
void initAll(int argc, char * argv[])
{
	initSeed(argv[1]);
	defineDefaultAtomsProp();
	initHBondsProperties();
#ifdef ENABLE_CL
	initCLProp();
#endif
#ifdef ENABLE_MPI
	MPI_Init( &argc, &argv );
#endif
}
/*****************************************************************************************************************************************************/
void finalize()
{
#ifdef ENABLE_MPI
	MPI_Finalize();
#endif
}
/*********************************************************************************/
double drandom()
{
	return (rand()/(double)RAND_MAX);
}
/*****************************************************************************************************************************************************/
/*     "normal" generates a random number from a normal Gaussian
     distribution with a mean of zero and a variance of one
*/
double normal()
{
	double v1,v2,rsq;
	double factor;
	static double store;
	static boolean compute = TRUE;

	if (compute)
	{
		do{
         		v1 = 2.0 * drandom()  - 1.0;
         		v2 = 2.0 * drandom () - 1.0;
         		rsq = v1*v1 + v2*v2;
		}while(rsq >= 1.0);
		compute = FALSE;
		factor = sqrt(-2.0*log(rsq)/rsq);
		store = v1 * factor;
		return v2 * factor;
      }
/*     use the second random value computed at the last call */
      else
      {
		compute = TRUE;
		return store;
      }
}
/*****************************************************************************************************************************************************/
Molecule* getFixedNormalModeSampling(char* inputFileName, int nModes, double* frequencies, double** modes, double* reducedMasses, double* quantumNumbers)
{
	Molecule* mol = NULL;
	int i,k,j;
	int nAtoms;
	double* Q = NULL;
	double* dQ = NULL;
	int ntr = 6;

	mol = readMolecule(inputFileName,TRUE);
	if(!mol) return NULL;
	nAtoms = mol->nAtoms;
	if(mol->klass->isLinear(mol)) ntr = 5;

	Q = malloc(nModes*sizeof(double));
	dQ = malloc(nModes*sizeof(double));
	for(i=0;i<nModes;i++) Q[i] = 0.0;
	for(i=0;i<nModes;i++) dQ[i] = 0.0;

	// Q and dQ in AU
	for(i=ntr;i<nModes;i++) 
	{
		double R = (normal()+1)/2.0;
		double fi = fabs(frequencies[i])/AUTOCM1;
		double Ei = (quantumNumbers[i]+0.5)*frequencies[i]/AUTOCM1; 
		double Ai = 0;
		if(Ei<0) Ei = 0;
		if(fi>1e-10) Ai = sqrt(2.0*fabs(Ei))/fi;
		Q[i] = Ai*cos(2*PI*R); 
		dQ[i] = -fi*Ai*sin(2*PI*R); 
	}
	// Q and dQ in AU
	for(i=ntr;i<nModes;i++)  Q[i]  /= sqrt(reducedMasses[i]*AMUTOAU);
	for(i=ntr;i<nModes;i++)  dQ[i] /= sqrt(reducedMasses[i]*AMUTOAU);

	// Q in Ang
	for(i=ntr;i<nModes;i++)  Q[i] *= BOHRTOANG;
	// dQ in Ang/AKMA-time 
	for(i=ntr;i<nModes;i++)  dQ[i] *= BOHRTOANG/(AUTOfs*fsInAKMA);
	
	for(j=0;j<nAtoms;j++) 
	for(k=0;k<3;k++) mol->atoms[j].velocity[k] = 0.0;

	for(i=0;i<nModes;i++) 
	{
		for(j=0;j<nAtoms;j++) 
		for(k=0;k<3;k++) 
		{
			int jd = 3*j+k; 
			 mol->atoms[j].velocity[k] += modes[jd][i]*dQ[i];
			 mol->atoms[j].coordinates[k] += modes[jd][i]*Q[i];
		}
	}
	free(Q);
	free(dQ);
	return mol;
}
/*****************************************************************************************************************************************************/
void printHarmonicVelocities(char* inputFileName, int nModes, double* frequencies, double** modes, double* reducedMasses)
{
	Molecule* mol = NULL;
	int i;
	double T=0;
	double* quantumNumbers = NULL;

	FILE* file = NULL;

	if(nModes<1) return;
	file = fopen(inputFileName,"rb");
	if(!file) return;
	readOneReal(file,"runTemp",&T);
        fclose(file);
	if(T<=0) T = 300;

	quantumNumbers = malloc(nModes*sizeof(double));
	for(i=0;i<nModes;i++) quantumNumbers[i] = 0.0;
	mol = getFixedNormalModeSampling(inputFileName, nModes, frequencies, modes, reducedMasses, quantumNumbers);
	if(!mol) return;

	mol->klass->removeTranslationAndRotation(mol);
	mol->klass->resetConstraints(mol, mol->constraints);
	mol->klass->scaleVelocities(mol, T);
        fprintf(stdout,"========================================================================================================================\n");
        fprintf(stdout,"# Geometry and velocities\n");
        fprintf(stdout,"# Fixed normal modes (ground, n=0) and scaled to T = %f\n",T);
        mol->klass->addGeometry(mol,stdout);
        mol->klass->addVelocities(mol,stdout);
        fprintf(stdout,"========================================================================================================================\n");

	mol->klass->free(mol);
	free(quantumNumbers);
}
/*****************************************************************************************************************************************************/
void addHarmonicVelocities(char* inputFileName, int nModes, double* frequencies, double** modes, double* reducedMasses, double* IRIntensities)
{
	Molecule* mol = NULL;
	FILE* file = NULL;
	int i;
	int dum = -1;
	double T;
	double *quantumNumbers = NULL;
	int ntr = 6;

	if(nModes<1) return;
	file = fopen(inputFileName,"rb");
	readOneInt(file,"HarmonicVelocityModes",&dum);
	readOneReal(file,"runTemp",&T);
	fclose(file);
	if(T<=0) T = 300;
	if(dum>nModes) 
	{
		fprintf(stderr,"Error in HarmonicVelocityModes keyword, the mode number must be <=i %d\n",nModes);
		return;
	}
	if(dum<=0) return;
	quantumNumbers = malloc(nModes*sizeof(double));
	for(i=0;i<nModes;i++) quantumNumbers[i] = 0.0;
	quantumNumbers[dum-1] = 1.0;
	mol = getFixedNormalModeSampling(inputFileName, nModes, frequencies, modes, reducedMasses, quantumNumbers);
	if(!mol) return;
	if(mol->klass->isLinear(mol)) ntr = 5;
	if(dum<=ntr) 
	{
		fprintf(stdout,"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n");
		fprintf(stdout,"Warning : The selected mode is a rotation or translation ; # of trans/rot = %d\n",ntr);
		fprintf(stdout,"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n");
	}

	mol->klass->removeTranslationAndRotation(mol);
	mol->klass->resetConstraints(mol, mol->constraints);
	mol->klass->scaleVelocities(mol, T);

	file = fopen(inputFileName,"ab+");
        mol->klass->addVelocities(mol,file);
	printf("Warning ! The velocities  have been added to the %s file\n",inputFileName);
	fclose(file);
	mol->klass->free(mol);
}
/****************************************************************************************************************************/
double* newVectorDouble(int n)
{
	double* v = malloc(n*sizeof(double));
	return v;
}
/****************************************************************************************************************************/
void initVectorDouble(double* v, int n, double val)
{
	int i;
	if(!v) return;
	for(i = 0;i<n; i++)  v[i] = val;
}
/****************************************************************************************************************************/
void freeVectorDouble(double** v)
{
	if(*v) free(*v);
	*v= NULL;
}
/****************************************************************************************************************************/
double** newMatrixDouble(int nrows, int ncolumns)
{
	double** M  = malloc(nrows*sizeof(double*));
	int i;
	for(i = 0;i<nrows; i++) 
		M[i] = malloc(ncolumns*sizeof(double));
	return M;
}
/****************************************************************************************************************************/
void freeMatrixDouble(double*** M, int nrows)
{
	if(*M) 
	{
		int i;
		for(i = 0;i<nrows; i++) 
			if((*M)[i])free((*M)[i]);
	}
	*M= NULL;
}
/****************************************************************************************************************************/
void initMatrixDouble(double** M, int nrows, int ncolumns, double val)
{
	int i,j;
	if(!M) return;
	for(i = 0;i<nrows; i++) 
		for(j = 0;j<ncolumns; j++)  M[i][j]  = val;
}
/****************************************************************************************************************************/
void printMatrixDouble(double** M, int nrows, int ncolumns)
{
	int i,j;
	for(i = 0;i<nrows; i++) 
	{
		for(j = 0;j<ncolumns; j++) 
      			printf("%f ",M[i][j]);
		printf("\n");
	}
}
/****************************************************************************************************************************/
double*** newCubeDouble(int nrows, int ncolumns, int nslices)
{
	double*** C  = malloc(nrows*sizeof(double**));
	int i,j;
	for(i = 0;i<nrows; i++) 
	{
		C[i] = malloc(ncolumns*sizeof(double*));
		for(j = 0;j<ncolumns; j++) 
			C[i][j] = malloc(nslices*sizeof(double));
	}
	return C;
}
/****************************************************************************************************************************/
void printCubeDouble(double*** C, int nrows, int ncolumns, int nslices)
{
	int i,j,k;
	for(i = 0;i<nrows; i++) 
	{
		for(j = 0;j<ncolumns; j++) 
		{
			for(k = 0;k<nslices; k++) 
      				printf("%f ",C[i][j][k]);
			printf("\n");
		}
		printf("\n");
	}
}
/****************************************************************************************************************************/
void initCubeDouble(double*** C, int nrows, int ncolumns, int nslices, double val)
{
	int i,j,k;
	if(!C) return;
	for(i = 0;i<nrows; i++) 
		for(j = 0;j<ncolumns; j++) 
			for(k = 0;k<nslices; k++) C[i][j][k] = val;
}
/****************************************************************************************************************************/
void freeCubeDouble(double**** C, int nrows, int ncolumns)
{
	if(*C) 
	{
		int i,j;
		for(i = 0;i<nrows; i++) 
		{
			if((*C)[i])
			for(j = 0;j<ncolumns; j++) 
				if((*C)[i][j]) free((*C)[i][j]);
	
			if((*C)[i])free((*C)[i]);
		}
	}
	*C= NULL;
}
/****************************************************************************************************************************/
int* newVectorInt(int n)
{
	int* v = malloc(n*sizeof(int));
	return v;
}
/****************************************************************************************************************************/
void freeVectorInt(int** v)
{
	if(*v) free(*v);
	*v= NULL;
}
/****************************************************************************************************************************/
int** newMatrixInt(int nrows, int ncolumns)
{
	int** M  = malloc(nrows*sizeof(int*));
	int i;
	for(i = 0;i<nrows; i++) 
		M[i] = malloc(ncolumns*sizeof(int));
	return M;
}
/****************************************************************************************************************************/
void freeMatrixInt(int*** M, int nrows)
{
	if(*M) 
	{
		int i;
		for(i = 0;i<nrows; i++) 
			if((*M)[i])free((*M)[i]);
	}
	*M= NULL;
}
/****************************************************************************************************************************/
int*** newCubeInt(int nrows, int ncolumns, int nslices)
{
	int*** C  = malloc(nrows*sizeof(int**));
	int i,j;
	for(i = 0;i<nrows; i++) 
	{
		C[i] = malloc(ncolumns*sizeof(int*));
		for(j = 0;j<ncolumns; j++) 
			C[i][j] = malloc(nslices*sizeof(int));
	}
	return C;
}
/****************************************************************************************************************************/
void freeCubeInt(int**** C, int nrows, int ncolumns)
{
	if(*C) 
	{
		int i,j;
		for(i = 0;i<nrows; i++) 
		{
			if((*C)[i])
			for(j = 0;j<ncolumns; j++) 
				if((*C)[i][j]) free((*C)[i][j]);
	
			if((*C)[i])free((*C)[i]);
		}
	}
	*C= NULL;
}
/****************************************************************************************************************************/
void initMatrixInt(int** M, int nrows, int ncolumns, int val)
{
	int i,j;
	if(!M) return;
	for(i = 0;i<nrows; i++) 
		for(j = 0;j<ncolumns; j++)  M[i][j]  = val;
}
/****************************************************************************************************************************/
void initCubeInt(int*** C, int nrows, int ncolumns, int nslices, int val)
{
	int i,j,k;
	if(!C) return;
	for(i = 0;i<nrows; i++) 
		for(j = 0;j<ncolumns; j++) 
			for(k = 0;k<nslices; k++) C[i][j][k] = val;
}
/****************************************************************************************************************************/
int**** newQuarticInt(int nrows, int ncolumns, int nslices, int nl)
{
	int**** C  = NULL;
	int i,j,k;
	if(nrows<1 || ncolumns<1 || nslices<1) return C;
	C  = malloc(nrows*sizeof(int***));
	for(i = 0;i<nrows; i++) 
	{
		C[i] = malloc(ncolumns*sizeof(int**));
		for(j = 0;j<ncolumns; j++) 
		{
			C[i][j] = malloc(nslices*sizeof(int*));
			for(k = 0;k<nslices; k++) 
				C[i][j][k] = malloc(nl*sizeof(int));
		}
	}
	return C;
}
/****************************************************************************************************************************/
void initQuarticInt(int**** C, int nrows, int ncolumns, int nslices, int nl, int val)
{
	int i,j,k,l;
	if(!C) return;
	for(i = 0;i<nrows; i++) 
		for(j = 0;j<ncolumns; j++) 
			for(k = 0;k<nslices; k++) 
				for(l = 0;l<nl; l++) 
					C[i][j][k][l] = val;
}
/****************************************************************************************************************************/
void freeQuarticInt(int***** C, int nrows, int ncolumns, int nslices)
{
	if(*C) 
	{
		int i,j,k;
		for(i = 0;i<nrows; i++) 
		{
			if((*C)[i])
			for(j = 0;j<ncolumns; j++) 
			{
				if((*C)[i][j])
				{
				for(k = 0;k<nslices; k++) 
					if((*C)[i][j][k]) free((*C)[i][j][k]);
				if((*C)[i][j]) free((*C)[i][j]);
				}
			}
			if((*C)[i])free((*C)[i]);
		}
		free(*C);
	}
	*C= NULL;
}
/****************************************************************************************************************************/
void printVectorDoubleCutOff(double* C, int n, double cutoff)
{
        int i;
        if(!C) return;
        for(i=0;i<n;i++) if(fabs(C[i])>=cutoff) printf("%d %20.10f\n",i+1,C[i]);
}
/****************************************************************************************************************************/
void printMatrixDoubleCutOff(double** M, int nrows, int ncolumns, double cutoff)
{
	int i,j;
	for(i = 0;i<nrows; i++) 
	{
		for(j = 0;j<ncolumns; j++) 
		if(fabs(M[i][j])>=cutoff)
      			printf("%d %d %20.10f\n",i+1,j+1,M[i][j]);
	}
}
/****************************************************************************************************************************/
void printCubeDoubleCutOff(double*** C, int nrows, int ncolumns, int nslices, double cutoff)
{
	int i,j,k;
	for(i = 0;i<nrows; i++) 
	{
		for(j = 0;j<ncolumns; j++) 
		{
			for(k = 0;k<nslices; k++) 
			if(fabs(C[i][j][k])>=cutoff)
      				printf("%d %d %d %20.10f\n",i+1,j+1,k+1,C[i][j][k]);
		}
	}
}
/****************************************************************************************************************************/
void printQuarticDouble(double**** C, int nrows, int ncolumns, int nslices, int nl)
{
	int i,j,k,l;
	for(i = 0;i<nrows; i++) 
	{
		for(j = 0;j<ncolumns; j++) 
		{
			for(k = 0;k<nslices; k++) 
			{
				for(l = 0;l<nl; l++) 
      					printf("%f ",C[i][j][k][l]);
				printf("\n");
			}
		}
		printf("\n");
	}
}
/****************************************************************************************************************************/
void printQuarticDoubleCutOff(double**** C, int nrows, int ncolumns, int nslices, int nl, double cutoff)
{
	int i,j,k,l;
	for(i = 0;i<nrows; i++) 
	{
		for(j = 0;j<ncolumns; j++) 
		{
			for(k = 0;k<nslices; k++) 
			{
				for(l = 0;l<nl; l++) 
				if(fabs(C[i][j][k][l])>=cutoff)
      					printf("%d %d %d %d %20.10f\n",i+1,j+1,k+1,l+1,C[i][j][k][l]);
			}
		}
	}
}
/****************************************************************************************************************************/
/**********************************************************************/
void getRandVect(double len, double V[])
{
	double l = 0;
	int j;
	for(j=0;j<3;j++)
	{
		V [j] = drandom();
		l += V[j]*V[j];
	}
	
	if(l<=0) return;
	l = sqrt(l);
	for(j=0;j<3;j++)
		V [j] *= len/l;
}
/**********************************************************************/
double erfinv( double y )
{
	static double a[] = {0,  0.886226899, -1.645349621,  0.914624893, -0.140543331 };
	static double b[] = {0, -2.118377725,  1.442710462, -0.329097515,  0.012229801 };
	static double c[] = {0, -1.970840454, -1.624906493,  3.429567803,  1.641345311 };
	static double d[] = {0,  3.543889200,  1.637067800 };
	double x=1e100, z;
  
	if ( y < -1. ) return x;
	if ( y >  1. ) return x;
	if ( y >= -.7 )
	{
		if ( y <= .7 )
		{
			z = y*y;
			x = y * (((a[4]*z+a[3])*z+a[2])*z+a[1]) /
			  ((((b[4]*z+b[3])*z+b[2])*z+b[1])*z+1);
		}
		else if ( y < 1 )
		{
			z = sqrt(-log((1-y)/2));
			x = (((c[4]*z+c[3])*z+c[2])*z+c[1]) / ((d[2]*z+d[1])*z+1);
		}
	}
	else
	{
  		z = sqrt(-log((1+y)/2));
  		x = -(((c[4]*z+c[3])*z+c[2])*z+c[1]) / ((d[2]*z+d[1])*z+1);
	}
	return x;
}
/**********************************************************************/
double	maxwel(double mass, double temperature)
{
	/* 
	 *  physical constants in SI units
	 *   ------------------------------
	 *      Kb = 1.380662 E-23 J/K
	 *      Na = 6.022045 E23  1/mol
	 *      e = 1.6021892 E-19 C
	 *      eps = 8.85418782 E-12 F/m
	 *                       
	 *      1 Kcal = 4184.0 J
	 *      1 amu = 1.6605655 E-27 Kg
	 *      1 A = 1.0 E-10 m
	 *                                       
	 *       Internally, AKMA units are used:
	 *       KBOLTZ = Na *Kb  / 1 Kcal
	 */
	/* double Kb = 6.022045e23*1.380662e-23/4184.0;*/
	double beta = sqrt(mass / (2.0*Kb*temperature));
	double rho;
	double xs, ys, zs;
	rho = drandom();
	xs = erfinv(rho)/beta;
	rho = drandom();
	ys = erfinv(rho)/beta;
	rho = drandom();
	zs = erfinv(rho)/beta;

	return sqrt(xs*xs+ys*ys+zs*zs);

}
/**************************************************/
boolean InverseTensor(double mat[3][3],double invmat[3][3])
{
	double t4,t6,t8,t10,t12,t14,t17;
	double d = 0;
	double precision = 1e-12;

	t4 = mat[0][0]*mat[1][1];     
 	t6 = mat[0][0]*mat[1][2];
      	t8 = mat[0][1]*mat[1][0];
      	t10 = mat[0][2]*mat[1][0];
      	t12 = mat[0][1]*mat[2][0];
      	t14 = mat[0][2]*mat[2][0];
      	d =(t4*mat[2][2]-t6*mat[2][1]-t8*mat[2][2]+t10*mat[2][1]+t12*mat[1][2]-t14*mat[1][1]);
	if(fabs(d)<precision) 
	{
      		invmat[0][0] = 0;
      		invmat[0][1] = 0;
      		invmat[0][2] = 0;
      		invmat[1][0] = 0;
      		invmat[1][1] = 0;
      		invmat[1][2] = 0;
      		invmat[2][0] = 0;
      		invmat[2][1] = 0;
      		invmat[2][2] = 0;
		return FALSE;
	}
      	t17 = 1/d;
      	invmat[0][0] = (mat[1][1]*mat[2][2]-mat[1][2]*mat[2][1])*t17;
      	invmat[0][1] = -(mat[0][1]*mat[2][2]-mat[0][2]*mat[2][1])*t17;
      	invmat[0][2] = -(-mat[0][1]*mat[1][2]+mat[0][2]*mat[1][1])*t17;
      	invmat[1][0] = -(mat[1][0]*mat[2][2]-mat[1][2]*mat[2][0])*t17;
      	invmat[1][1] = (mat[0][0]*mat[2][2]-t14)*t17;
      	invmat[1][2] = -(t6-t10)*t17;
      	invmat[2][0] = -(-mat[1][0]*mat[2][1]+mat[1][1]*mat[2][0])*t17;
      	invmat[2][1] = -(mat[0][0]*mat[2][1]-t12)*t17;
      	invmat[2][2] = (t4-t8)*t17;

	return TRUE;
}
/**************************************************/
void computeAngularVelocitiesForALinearMolecule(double* R1, double* R2, double  inert[3][3], double*L, double* vAng)
{
	double rmol[3];/* normal vector parallel to the molecule */
	double ra[3]; /* first vector perpendiculaire to the molecule */
	double rb[3]; /* second vector perpendiculaire to the molecule */
	double dum[3];
	double precision = 1e-12;
	int j,k;
	double d;
	double a00,a01,a11;
	double b0,b1;
	double w0,w1;

	for(j=0;j<3;j++) vAng[j] = 0;

	for(j=0;j<3;j++) rmol[j] = R2[j]- R1[j];

	d = 0;
	for(j=0;j<3;j++) d+=rmol[j]*rmol[j];
	if(d<precision) return;
	d = 1/sqrt(d);
	for(j=0;j<3;j++) rmol[j] *= d;

/*      find two orthogonal vectors to molecule coordinate frame */
	k = 0;
	for(j=1;j<3;j++) if(fabs(rmol[k])>fabs(rmol[j])) k = j;

	for(j=0;j<3;j++) ra[j] = -rmol[k] * rmol[j];
	ra[k] += 1.0;

	d = 0;
	for(j=0;j<3;j++) d+=ra[j]*ra[j];
	if(d<precision) return;
	d = 1/sqrt(d);
	for(j=0;j<3;j++) ra[j] *= d;

	for(j=0;j<3;j++) rb[j] = ra[(j+1)%3]*rmol[(j+2)%3]-rmol[(j+1)%3]*ra[(j+2)%3];

/*      solve the 2 linear system for angular velocity */
	for(k=0;k<3;k++) dum[k] = 0;
	for(k=0;k<3;k++) for(j=0;j<3;j++) dum[k] += inert[j][k]*ra[j];

	a00 = 0;
	for(j=0;j<3;j++) a00 +=ra[j]*dum[j];

	for(k=0;k<3;k++) dum[k] = 0;
	for(k=0;k<3;k++) for(j=0;j<3;j++) dum[k] += inert[j][k]*rb[j];

	a01 = 0;
	for(j=0;j<3;j++) a01 += ra[j]*dum[j];

	a11 = 0;
	for(j=0;j<3;j++) a11 += rb[j]*dum[j];

	b0 = 0;
	for(j=0;j<3;j++) b0 += ra[j]*L[j];

	b1 = 0;
	for(j=0;j<3;j++) b1 += rb[j]*L[j];

	w0 = (a01*b1-a11*b0) / (a01*a01-a00*a11);
	w1 = (b1-a01*w0) / a11;

	for(j=0;j<3;j++) vAng[j] = w0*ra[j] + w1*rb[j];
}
/*************************************************************************************/
boolean readVectorReal(FILE* file, char* tag, int n, double*values)
{
	static char *t = NULL; 
	char* TAG = NULL;
	char* pos = NULL;
	int i;
	int ii;
	double v;
	if(!tag) return FALSE;
	if(!values) return FALSE;
	if(t==NULL) t = malloc(BSIZE*sizeof(char));

	TAG = strdup(tag);
	uppercase(TAG);
	rewind(file);

	while(!feof(file))
  	{
    		if(!fgets(t,BSIZE, file)) break;
		deleteFirstSpaces(t);
		if(t[0]=='#') continue;
		uppercase(t);
		pos = strstr(t,TAG);
		if(!pos) continue;
		break;
	}
	if(!pos) 
	{
		free(TAG);
		return FALSE;
	}
	for(i=0;i<n;i++)
	{
    		if(!fgets(t,BSIZE, file)) break;
		if(2==sscanf(t,"%d %lf",&ii, &v)) 
		{
			if(ii<=n && ii>0) values[ii-1] = v;
			else
			{
				fprintf(stderr,"Erreur dans les donnees de %s: verifie les indices\n",tag);
				exit(1);
			}
		}
		else break;
	}
	free(TAG);
	if(i<=0) return FALSE;
	return TRUE;
}
/********************************************************************************/
int readEnergyAndDipoleFromGabeditFile(char* fileName, double energy[], double dipole[])
{
	FILE* file = NULL;
	static char *t = NULL; 
	int n;
	int i;
	char* pos;
	int nGeoms = 0;
	int nLabels = 0;
	int iEnergy, iDipole;

	if(!fileName) 
	{
		printf("Sorry I cannot read energy and dipole fileName = NULL\n");
		exit(1);
	}
	if(t==NULL) t = malloc(BSIZE*sizeof(char));
	file = fopen(fileName,"rb");
	if(!file)
	{
		printf("Sorry I cannot open %s file\n",fileName);
		exit(1);
	}
	rewind(file);
	iEnergy = -1;
	iDipole = -1;
	while(!feof(file))
  	{
    		if(!fgets(t,BSIZE, file)) break;
		deleteFirstSpaces(t);
		if(t[0]=='#') continue;
		uppercase(t);
		pos = strstr(t,"[GEOMS]");
		if(pos)
		{ 
			while(!feof(file))
			{
    				if(!fgets(t,BSIZE, file)) break;
				deleteFirstSpaces(t);
				if(t[0]=='#') continue;
				break;
			}
			if(2==sscanf(t,"%d%d",&nGeoms,&nLabels))
			{

				for(i=0; i<nLabels; i++)
				{
    					if(!fgets(t,BSIZE, file)) break;
					uppercase(t);
					if(strstr(t,"ENERGY") && strstr(t,"KCAL")) iEnergy = i;
					if(strstr(t,"DIPOLE")) iDipole = i;
				}
			}
			break;
		}
	}
	/* I read the first geometry */
	n = 0;
	for(i=0; i<nLabels; i++)
	{
    		if(!fgets(t,BSIZE, file)) break;
		if(i==iEnergy) { n++; sscanf(t,"%lf",&energy[0]);}
		if(i==iDipole) { n++; sscanf(t,"%lf %lf %lf",&dipole[0], &dipole[1],&dipole[2]);}
	}
	fclose(file);
	return n;
}
/********************************************************************************/
double**** newQuarticDouble(int nrows, int ncolumns, int nslices, int nl)
{
	double**** C  = NULL;
	int i,j,k;
	if(nrows<1 || ncolumns<1 || nslices<1) return C;
	C  = malloc(nrows*sizeof(double***));
	for(i = 0;i<nrows; i++) 
	{
		C[i] = malloc(ncolumns*sizeof(double**));
		for(j = 0;j<ncolumns; j++) 
		{
			C[i][j] = malloc(nslices*sizeof(double*));
			for(k = 0;k<nslices; k++) C[i][j][k] = malloc(nl*sizeof(double));
		}
	}
	return C;
}
/****************************************************************************************************************************/
void initQuarticDouble(double**** C, int nrows, int ncolumns, int nslices, int nl, double val)
{
	int i,j,k,l;
	if(!C) return;
	for(i = 0;i<nrows; i++) 
		for(j = 0;j<ncolumns; j++) 
			for(k = 0;k<nslices; k++) 
				for(l = 0;l<nl; l++) 
					C[i][j][k][l] = val;
}
/****************************************************************************************************************************/
void freeQuarticDouble(double***** C, int nrows, int ncolumns, int nslices)
{
	if(*C) 
	{
		int i,j,k;
		for(i = 0;i<nrows; i++) 
		{
			if((*C)[i])
			for(j = 0;j<ncolumns; j++) 
			{
				if((*C)[i][j])
				{
				for(k = 0;k<nslices; k++) 
					if((*C)[i][j][k]) free((*C)[i][j][k]);
				if((*C)[i][j]) free((*C)[i][j]);
				}
			}
			if((*C)[i])free((*C)[i]);
		}
		free(*C);
	}
	*C= NULL;
}
/********************************************************************************/
double***** newQuinticDouble(int nrows, int ncolumns, int nslices, int nl, int n5)
{
	double***** C  = NULL;
	int i,j,k,l;
	if(nrows<1 || ncolumns<1 || nslices<1 || nl<1 || n5<1) return C;
	C  = malloc(nrows*sizeof(double****));
	for(i = 0;i<nrows; i++) 
	{
		C[i] = malloc(ncolumns*sizeof(double***));
		for(j = 0;j<ncolumns; j++) 
		{
			C[i][j] = malloc(nslices*sizeof(double**));
			for(k = 0;k<nslices; k++) 
			{
				C[i][j][k] = malloc(nl*sizeof(double*));
				for(l = 0;l<nl; l++) 
					C[i][j][k][l] = malloc(n5*sizeof(double));
			}
		}
	}
	return C;
}
/****************************************************************************************************************************/
void initQuinticDouble(double***** C, int nrows, int ncolumns, int nslices, int nl, int n5, double val)
{
	int i,j,k,l,n;
	if(!C) return;
	for(i = 0;i<nrows; i++) 
		for(j = 0;j<ncolumns; j++) 
			for(k = 0;k<nslices; k++) 
				for(l = 0;l<nl; l++) 
				for(n = 0;n<n5; n++) 
					C[i][j][k][l][n] = val;
}
/****************************************************************************************************************************/
void freeQuinticDouble(double****** C, int nrows, int ncolumns, int nslices, int nl)
{
	if(*C) 
	{
		int i,j,k,l;
		for(i = 0;i<nrows; i++) 
		{
			if((*C)[i])
			for(j = 0;j<ncolumns; j++) 
			{
				if((*C)[i][j])
				{
					for(k = 0;k<nslices; k++) 
						if((*C)[i][j][k]) 
						{
							for(l = 0;l<nl; l++) 
								if((*C)[i][j][k][l]) free((*C)[i][j][k][l]);
							free((*C)[i][j][k]);
						}
					if((*C)[i][j]) free((*C)[i][j]);
				}
			}
			if((*C)[i])free((*C)[i]);
		}
		free(*C);
	}
	*C= NULL;
}
/****************************************************************************************************************************/
void symmetrizeMatrixDouble(double** M, int nrows, int ncolumns, double cutOff)
{
	int i,j;
	double x,y;
	if(!M) return;
	for(i = 0;i<nrows; i++) 
		for(j = 0;j<i; j++)  
		{
			if(j>ncolumns-1) continue;
			x =  M[i][j];
			y =  M[j][i];
			if(fabs(x)>cutOff && fabs(y)>cutOff)  M[i][j] = M[j][i] = (x+y)/2;
			else if(fabs(x)>cutOff)  M[i][j] = M[j][i] = x;
			else if(fabs(y)>cutOff)  M[i][j] = M[j][i] = y;
			else M[i][j]  = 0.0;
		}
}
/****************************************************************************************************************************/
static void symmetrizeCubeDoubleIJ(double*** C, int nrows, int ncolumns,  int nslices, double cutOff)
{
	int i,j,k;
	double x,y;
	for(i = 0;i<nrows; i++) 
		for(j = 0;j<i; j++)  
		{
			if(j>ncolumns-1) continue;
			for(k = 0;k<nslices; k++)  
			{
				x =  C[i][j][k];
				y =  C[j][i][k];
				if(fabs(x)>cutOff && fabs(y)>cutOff)  C[i][j][k] = C[j][i][k] = (x+y)/2;
				else if(fabs(x)>cutOff)  C[i][j][k] = C[j][i][k] = x;
				else if(fabs(y)>cutOff)  C[i][j][k] = C[j][i][k] = y;
				else C[i][j][k]  = 0.0;
			}
		}
}
/****************************************************************************************************************************/
static void symmetrizeCubeDoubleJK(double*** C, int nrows, int ncolumns,  int nslices, double cutOff)
{
	int i;
	for(i=0;i<nrows;i++) symmetrizeMatrixDouble(C[i], ncolumns, nslices,  cutOff);
}
/****************************************************************************************************************************/
void symmetrizeCubeDouble(double*** C, int nrows, int ncolumns,  int nslices, double cutOff)
{
	symmetrizeCubeDoubleIJ(C, nrows, ncolumns,  nslices, cutOff);
	symmetrizeCubeDoubleJK(C, nrows, ncolumns,  nslices, cutOff);
	symmetrizeCubeDoubleIJ(C, nrows, ncolumns,  nslices, cutOff);
}
/****************************************************************************************************************************/
/****************************************************************************************************************************/
static void symmetrizeQuarticDoubleIJ(double**** Q, int nrows, int ncolumns,  int nslices, int nq, double cutOff)
{
	int i,j,k,l;
	double x,y;
	for(i = 0;i<nrows; i++) 
		for(j = 0;j<i; j++)  
		{
			if(j>ncolumns-1) continue;
			for(k = 0;k<nslices; k++)  
			for(l = 0;l<nq; l++)  
			{
				x =  Q[i][j][k][l];
				y =  Q[j][i][k][l];
				if(fabs(x)>cutOff && fabs(y)>cutOff)  Q[i][j][k][l] = Q[j][i][k][l] = (x+y)/2;
				else if(fabs(x)>cutOff)  Q[i][j][k][l] = Q[j][i][k][l] = x;
				else if(fabs(y)>cutOff)  Q[i][j][k][l] = Q[j][i][k][l] = y;
				else Q[i][j][k][l]  = 0.0;
			}
		}
}
/****************************************************************************************************************************/
static void symmetrizeQuarticDoubleJKL(double**** Q, int nrows, int ncolumns,  int nslices, int nq, double cutOff)
{
	int i;
	for(i=0;i<nrows;i++) symmetrizeCubeDouble(Q[i], ncolumns, nslices, nq,  cutOff);
}
/****************************************************************************************************************************/
void symmetrizeQuarticDouble(double**** Q, int nrows, int ncolumns,  int nslices, int nq, double cutOff)
{
	symmetrizeQuarticDoubleIJ(Q, nrows, ncolumns,  nslices, nq, cutOff);
	symmetrizeQuarticDoubleJKL(Q, nrows, ncolumns,  nslices, nq, cutOff);
	symmetrizeQuarticDoubleIJ(Q, nrows, ncolumns,  nslices, nq, cutOff);
}
/****************************************************************************************************************************/
boolean goToStr(FILE* file, char* tag)
{
        static char *t = NULL;
        char* TAG = NULL;
        char* pos = NULL;
        if(!tag) return FALSE;
        if(t==NULL) t = malloc(BSIZE*sizeof(char));

        TAG = strdup(tag);
        uppercase(TAG);
        rewind(file);

        while(!feof(file))
        {
                if(!fgets(t,BSIZE, file)) break;
                deleteFirstSpaces(t);
                if(t[0]=='#') continue;
                uppercase(t);
                pos = strstr(t,TAG);
                if(!pos) continue;
                break;
        }
        return pos != NULL;
}

/********************************************************************************/
void computeAndPrintHarmonicVelocitiesCoordinates(Molecule* mol, char* inputFileName, boolean addToInputFile, boolean changeGeom)
{
	int numMode = -1;
	double T = -1;
	FILE* file = fopen(inputFileName,"rb");
	if(!file) return;
	readOneInt(file,"HarmonicVelocityModes",&numMode);
	readOneReal(file,"runTemp",&T);
	fclose(file);
	if(T<=0) T = 300;
	if(numMode>mol->vibration.nModes) 
	{
		fprintf(stderr,"Error in HarmonicVelocityModes keyword, the mode number must be <=i %d\n",mol->vibration.nModes);
		return;
	}
	mol->klass->computeHarmonicVelocitiesCoordinates(mol, T, numMode, changeGeom);
	file = stdout;
	if(addToInputFile) file = fopen(inputFileName,"ab+");
        fprintf(file,"#========================================================================================================================\n");
	if(numMode<1) fprintf(file,"# Fixed normal modes (n=1 for mode %d, n=0 for all other modes) and scaled to T = %f\n",numMode,T);
	else fprintf(file,"# Fixed normal modes (ground, n=0) and scaled to T = %f\n",T);
	if(changeGeom)
	{  
        	fprintf(stdout,"# Geometry and velocities\n");
        	mol->klass->addGeometry(mol,file);
	}
        fprintf(file,"# Harmonic Velocities\n");
        mol->klass->addVelocities(mol,file);
        fprintf(file,"#========================================================================================================================\n");
	if(addToInputFile) fclose(file);
}
/********************************************************************************/
void computeFrequenciesFromFilesDlg(char* inputFileName, boolean oneStep)
{
	Molecule* mol = readMolecule(inputFileName,TRUE);
	double dx = 1e-3;
	char* fileNameOut = strdup_printf("%sFreq.gab",getSuffixNameFile(inputFileName));
	 FILE* file = fopen(inputFileName,"rb");

	readOneReal(file,"dx",&dx);
	fclose(file);

	if(oneStep) mol->klass->computeFrequenciesOneStepFromFiles(mol, inputFileName,  dx);
	else mol->klass->computeFrequenciesFromFiles(mol, inputFileName,  dx);

	printf("Frequencies and modes in the %s file\n",fileNameOut);
	mol->klass->save(mol, fileNameOut);
	computeAndPrintHarmonicVelocitiesCoordinates(mol, inputFileName, FALSE, TRUE);
	computeAndPrintHarmonicVelocitiesCoordinates(mol, inputFileName, TRUE, FALSE);

	mol->klass->free(mol);
	free(fileNameOut);
}
/********************************************************************************/
void computeFrequenciesFromGradFilesDlg(char* inputFileName, boolean oneStep)
{
	Molecule* mol = readMolecule(inputFileName,TRUE);
	double dx = 1e-3;
	char* fileNameOut = strdup_printf("%sFreq.gab",getSuffixNameFile(inputFileName));
	 FILE* file = fopen(inputFileName,"rb");

	readOneReal(file,"dx",&dx);
	fclose(file);

	if(oneStep) mol->klass->computeFrequenciesOneStepFromGradFiles(mol, inputFileName,  dx);
	else mol->klass->computeFrequenciesFromGradFiles(mol, inputFileName,  dx);

	printf("Frequencies and modes in the %s file\n",fileNameOut);
	mol->klass->save(mol, fileNameOut);
	computeAndPrintHarmonicVelocitiesCoordinates(mol, inputFileName, FALSE, TRUE);
	computeAndPrintHarmonicVelocitiesCoordinates(mol, inputFileName, TRUE, FALSE);
	mol->klass->free(mol);
	free(fileNameOut);
}
/********************************************************************************/
void generateCChemIFilesForFrequenciesDlg(char* inputFileName, boolean oneStep)
{
	Molecule* mol = readMolecule(inputFileName,TRUE);
	double dx = 1e-3;
	char* fileNameOut = strdup_printf("%sFreq.gab",getSuffixNameFile(inputFileName));
	 FILE* file = fopen(inputFileName,"rb");

	readOneReal(file,"dx",&dx);
	fclose(file);

	if(oneStep) mol->klass->generateCChemIFilesOneStepForFrequencies(mol, inputFileName,  dx);
	else mol->klass->generateCChemIFilesForFrequencies(mol, inputFileName,  dx);

	mol->klass->free(mol);
	free(fileNameOut);
}
/********************************************************************************/
void generateCChemIGradFilesForFrequenciesDlg(char* inputFileName, boolean oneStep)
{
	Molecule* mol = readMolecule(inputFileName,TRUE);
	double dx = 1e-3;
	char* fileNameOut = strdup_printf("%sFreq.gab",getSuffixNameFile(inputFileName));
	 FILE* file = fopen(inputFileName,"rb");

	readOneReal(file,"dx",&dx);
	fclose(file);

	if(oneStep) mol->klass->generateCChemIGradFilesOneStepForFrequencies(mol, inputFileName,  dx);
	else mol->klass->generateCChemIGradFilesForFrequencies(mol, inputFileName,  dx);

	mol->klass->free(mol);
	free(fileNameOut);
}
/**********************************************************************************************************************************/
int get_num_orbitals_from_aux_mopac_file(FILE* file, char* blockName,  int* begin, int* end)
{
	char t[BSIZE];
	*begin = 0;
	*end = 0;
	 while(!feof(file))
	 {
		if(!fgets(t,BSIZE,file))break;
		if(strstr( t, blockName))
		{
			char* pdest = strstr( t, "=")+1;
			int i = sscanf(pdest,"%d %d",begin,end);
			return i;
		}
	 }
	 return 0;
}
/**********************************************************************************************************************************/
char** get_one_block_from_aux_mopac_file(FILE* file, char* blockName,  int* n)
{
	int nElements = 0;
	char** elements = NULL;
	char t[BSIZE];
	 while(!feof(file))
	 {
		if(!fgets(t,BSIZE,file))break;
		if(strstr( t, blockName))
		{
			char* pdest = strstr( t, "[")+1;
			int i;
			nElements = atoi(pdest);
			if(nElements<1) break;
			else
			{ 
				long int geomposok = ftell(file);
				if(!fgets(t,BSIZE,file))break;
				if(!strstr(t,"# ")) fseek(file, geomposok, SEEK_SET);
			}

			elements = malloc(nElements*sizeof(char*));
			for(i=0;i<nElements;i++)
			{
				int k;
				elements[i] = malloc(100*sizeof(char));
				k = fscanf(file,"%s",elements[i]);
				if(k<1 || strstr(elements[i],"["))
				{
					if(elements)
					{
						for(i=0;i<nElements;i++)
							if(elements[i]) free(elements[i]);
						free(elements);
						elements = NULL;
					}
					break;
				}
				else
				{
					if(!strstr(blockName,"ATOM_EL"))
					for(k=0;k<strlen(elements[i]);k++)
					{
						if(elements[i][k]=='D') elements[i][k]='e';
						if(elements[i][k]=='d') elements[i][k]='e';
					}
				}
			}
			break;
		}
	 }
	 *n = nElements;
	 return elements;
}
/**********************************************************************************************************************************/
char** free_one_string_table(char** table, int n)
{
	if(table)
	{
		int i;
		for(i=0;i<n;i++)
			if(table[i]) free(table[i]);
		free(table);
	}
	return NULL;
}
/********************************************************************************/
void get_dipole_from_gamess_output_file(FILE* file, double dipole[])
{
  	char *t = malloc(BSIZE*sizeof(char));
  	char* t1;

	dipole[0] = dipole[1] = dipole[2] = 0.0;

  	while(!feof(file) )
	{
    		t1 = NULL;
    		if(!fgets(t,BSIZE,file))break;
    		t1 = strstr( t, "ELECTROSTATIC MOMENTS");
		if(t1)
		{
  			while(!feof(file) )
			{
    				if(!fgets(t,BSIZE,file))break;
    				t1 = strstr( t, "DEBYE");
				if(t1)
				{
    					if(!fgets(t,BSIZE,file))break;
					sscanf(t,"%lf %lf %lf",&dipole[0],&dipole[1],&dipole[2]);
					break;
				}
			}
			break;
		}
		else
		{
			if(strstr( t, "END OF PROPERTY" )) break;
		}

	}
	free(t);
}
/********************************************************************************/
void get_dipole_from_turbomole_output_file(FILE* file, double dipole[])
{
  	char *t = malloc(BSIZE*sizeof(char));
  	char dum[100];


	dipole[0] = dipole[1] = dipole[2] = 0.0;

  	while(!feof(file) )
	{
    		if(!fgets(t,BSIZE,file))break;
		if(strstr( t, "electrostatic moments"))
		{
  			while(!feof(file) )
			{
    				if(!fgets(t,BSIZE,file))break;
				if(strstr( t, "dipole moment"))
				{
					double d;
    					if(!fgets(t,BSIZE,file))break;
    					if(!fgets(t,BSIZE,file))break;/* x */
					sscanf(t,"%s %lf %lf %lf",dum, &d, &d, &dipole[0]);
    					if(!fgets(t,BSIZE,file))break;/* y */
					sscanf(t,"%s %lf %lf %lf",dum, &d, &d, &dipole[1]);
    					if(!fgets(t,BSIZE,file))break;/* z */
					sscanf(t,"%s %lf %lf %lf",dum, &d, &d, &dipole[2]);
					break;
				}
			}
			break;
		}
	}
	free(t);
}
/********************************************************************************/
void get_dipole_from_gaussian_output_file(FILE* file, double dipole[])
{
  	char *t = malloc(BSIZE*sizeof(char));
  	char* pdest;
	int ngrad = 0;

	dipole[0] = dipole[1] = dipole[2] = 0.0;

  	while(!feof(file) )
	{
    		pdest = NULL;
		dipole[0] = dipole[1] = dipole[2] = 0.0;
    		if(!fgets(t,BSIZE,file)) break;
    		pdest = strstr( t, "Dipole moment (Debye)");

		if(strstr( t, "Dipole moment") && strstr( t, "Debye")) /* field-independent basis */
		{
    		if(!fgets(t,BSIZE,file)) break;
    		pdest = strstr( t, "X=")+2;
		sscanf(pdest,"%lf",&dipole[0]);
    		pdest = strstr( t, "Y=")+2;
		sscanf(pdest,"%lf",&dipole[1]);
    		pdest = strstr( t, "Z=")+2;
		sscanf(pdest,"%lf",&dipole[2]);
		/*
		Debug("t =%s\n",t);
		*/
		break;
		}
		else
		{
          		pdest = strstr( t, "GradGradGrad" );
			if(pdest)
			{
				ngrad++;
			/*	Debug("ngrad = %d\n",ngrad);*/
			}
			if(ngrad>2)
				break;
		}

	}
	free(t);
}
/********************************************************************************/
void get_dipole_from_molpro_output_file(FILE* file, double dipole[])
{
  	char *t = malloc(BSIZE*sizeof(char));
  	char* t1;
  	char* t2;

	dipole[0] = dipole[1] = dipole[2] = 0.0;

  	while(!feof(file) )
	{
    		t1 = NULL;
    		if(!fgets(t,BSIZE,file)) break;
    		t1 = strstr( t, "DIPOLE MOMENTS:");

		if(t1)
		{
    		t2 = strstr( t1, ":")+2;
		sscanf(t2,"%lf %lf %lf",&dipole[0],&dipole[1],&dipole[2]);
		/*
		Debug("t =%s\n",t);
		*/
		break;
		}
		else
		{
          		t1 = strstr( t, "GEOMETRY OPTIMIZATION STEP" );
			if(t1)
				break;
          		t1 = strstr( t, "SEWARD" );
			if(t1)
				break;
		}

	}
	free(t);
}
/********************************************************************************/
void get_dipole_from_dalton_output_file(FILE* file, double dipole[])
{
  	char *t = malloc(BSIZE*sizeof(char));
  	char* t1;
	char dum[100];

	dipole[0] = dipole[1] = dipole[2] = 0.0;

  	while(!feof(file) )
	{
    		t1 = NULL;
    		if(!fgets(t,BSIZE,file))break;
    		t1 = strstr( t, "Dipole moment components");
		if(t1)
		{
    			if(!fgets(t,BSIZE,file))break;
    			if(!fgets(t,BSIZE,file))break;
    			if(!fgets(t,BSIZE,file))break;
    			if(!fgets(t,BSIZE,file))break;
    			if(!fgets(t,BSIZE,file))break;
			sscanf(t,"%s %lf",dum, &dipole[0]);
    			if(!fgets(t,BSIZE,file))break;
			sscanf(t,"%s %lf",dum, &dipole[1]);
    			if(!fgets(t,BSIZE,file))break;
			sscanf(t,"%s %lf",dum, &dipole[2]);
		/*
			Debug("t =%s\n",t);
		*/
			break;
		}
		else
		{
			if(strstr( t, ">>>>" )) break;
		}

	}
	free(t);
}
/********************************************************************************/
void get_dipole_from_orca_output_file(FILE* file, double dipole[])
{
  	char *t = malloc(BSIZE*sizeof(char));
  	char* pdest;

	dipole[0] = dipole[1] = dipole[2] = 0.0;

  	while(!feof(file) )
	{
    		pdest = NULL;
		dipole[0] = dipole[1] = dipole[2] = 0.0;
    		if(!fgets(t,BSIZE,file)) break;
    		pdest = strstr( t, "Total Dipole Moment");

		if(pdest && strstr( t,":"))
		{
			pdest = strstr( t,":")+1;
			sscanf(pdest,"%lf %lf %lf",&dipole[0],&dipole[1],&dipole[2]);
			break;
		}
	}
	free(t);
}
/********************************************************************************/
void get_dipole_from_vasp_output_file(FILE* file, double dipole[])
{
  	char *t = malloc(BSIZE*sizeof(char));

	dipole[0] = dipole[1] = dipole[2] = 0.0;

	free(t);
}
/********************************************************************************/
void get_dipole_from_nwchem_output_file(FILE* file, double dipole[])
{
  	char *t = malloc(BSIZE*sizeof(char));
  	char* pdest;

	dipole[0] = dipole[1] = dipole[2] = 0.0;

  	while(!feof(file) )
	{
    		pdest = NULL;
		dipole[0] = dipole[1] = dipole[2] = 0.0;
    		if(!fgets(t,BSIZE,file)) break;
    		pdest = strstr( t, "Nuclear Dipole moment");
		if(pdest)
		{
			boolean OK = FALSE;
  			while(!feof(file) )
			{
    				if(!fgets(t,BSIZE,file)) break;
				if(strstr(t,"---------------- ---------------- ----------------"))
				{
					OK = TRUE;
					break;
				}
			}
			if(!OK) break;
    			if(!fgets(t,BSIZE,file)) break;
			sscanf(pdest,"%lf %lf %lf",&dipole[0],&dipole[1],&dipole[2]);
			break;
		}
	}
	free(t);
}
/********************************************************************************/
void get_dipole_from_psicode_output_file(FILE* file, double dipole[])
{
  	char *t = malloc(BSIZE*sizeof(char));
  	char* pdest;

	dipole[0] = dipole[1] = dipole[2] = 0.0;

  	while(!feof(file) )
	{
    		pdest = NULL;
		dipole[0] = dipole[1] = dipole[2] = 0.0;
    		if(!fgets(t,BSIZE,file)) break;
    		pdest = strstr( t, "Nuclear Dipole moment");
		if(pdest)
		{
			boolean OK = FALSE;
  			while(!feof(file) )
			{
    				if(!fgets(t,BSIZE,file)) break;
				if(strstr(t,"---------------- ---------------- ----------------"))
				{
					OK = TRUE;
					break;
				}
			}
			if(!OK) break;
    			if(!fgets(t,BSIZE,file)) break;
			sscanf(pdest,"%lf %lf %lf",&dipole[0],&dipole[1],&dipole[2]);
			break;
		}
	}
	free(t);
}
/********************************************************************************/
void get_dipole_from_qchem_output_file(FILE* file, double dipole[])
{
  	char *t = malloc(BSIZE*sizeof(char));
  	char* pdest;
	int ngrad = 0;

	dipole[0] = dipole[1] = dipole[2] = 0.0;

  	while(!feof(file) )
	{
    		pdest = NULL;
		dipole[0] = dipole[1] = dipole[2] = 0.0;
    		if(!fgets(t,BSIZE,file)) break;
    		pdest = strstr( t, "Dipole Moment (Debye)");

		if(pdest)
		{
    		if(!fgets(t,BSIZE,file)) break;
		else break;
    		pdest = strstr( t, "X")+2;
		if(pdest) sscanf(pdest,"%lf",&dipole[0]);
    		pdest = strstr( t, "Y")+2;
		if(pdest) sscanf(pdest,"%lf",&dipole[1]);
    		pdest = strstr( t, "Z")+2;
		if(pdest) sscanf(pdest,"%lf",&dipole[2]);
		break;
		}
		else
		{
          		pdest = strstr( t, "GradGradGrad" );
			if(pdest)
			{
				ngrad++;
			}
			if(ngrad>2) break;
		}

	}
	free(t);
}
/********************************************************************************/
void get_dipole_from_mopac_output_file(FILE* file, double dipole[])
{
  	char *t = malloc(BSIZE*sizeof(char));
  	char* pdest;
	char dum[100];

	dipole[0] = dipole[1] = dipole[2] = 0.0;

  	while(!feof(file) )
	{
    		pdest = NULL;
		dipole[0] = dipole[1] = dipole[2] = 0.0;
    		if(!fgets(t,BSIZE,file)) break;
    		pdest = strstr( t, "DIPOLE           X         Y         Z");

		if(pdest)
		{
    			if(!fgets(t,BSIZE,file)) break;
    			if(!fgets(t,BSIZE,file)) break;
    			if(!fgets(t,BSIZE,file)) break;
    			pdest = strstr( t, "SUM")+2;
			sscanf(t,"%s %lf %lf %lf",dum,&dipole[0],&dipole[1], &dipole[2]);
			break;
		}
	}
	free(t);
}
/**********************************************************************************************************************************/
void get_dipole_from_mopac_aux_file(FILE* file, double dipole[])
{
  	char t[BSIZE];
  	char* pdest;
	int i;

	rewind(file);
	for(i=0;i<3;i++) dipole[i] = 0;
  	while(!feof(file) )
	{
    		pdest = NULL;
		if(!fgets(t,BSIZE,file))break;
    		pdest = strstr( t, "DIPOLE:DEBYE=");

		if(pdest)
		{
    			pdest = strstr( t, "=")+1;
			for(i=0;i<3;i++) dipole[i] = 0;
			if(pdest) sscanf(pdest,"%lf %lf %lf",&dipole[0], &dipole[1],&dipole[2]);
			break;
		}
	}
}
/********************************************************************************/
void changeDInE(char *st)
{
        int i;
        int l = 0;
        if(!st) return;
        l = strlen(st);
        for(i=0;i<l;i++) if(st[i] == 'D' || st[i] =='d') st[i]='e';
}

/****************************************************************************/
int get_one_int_from_fchk_gaussian_file(FILE* file, char* blockName)
{
	int ipos = 47;
	char t[BSIZE];
	 while(!feof(file))
	 {
		if(!fgets(t,BSIZE,file))break;
		if(strstr( t, blockName))
		{
			if(strlen(t)>ipos+1) return atoi(t+ipos);
			return -1;
		}
	 }
	 return -1;
}
/****************************************************************************/
double get_one_real_from_fchk_gaussian_file(FILE* file, char* blockName)
{
	int ipos = 47;
	char t[BSIZE];
	 while(!feof(file))
	 {
		if(!fgets(t,BSIZE,file))break;
		if(strstr( t, blockName))
		{
			if(strlen(t)>ipos+1) return atof(t+ipos);
			return -1;
		}
	 }
	 return -1;
}
/****************************************************************************/
int* get_array_int_from_fchk_gaussian_file(FILE* file, char* blockName, int* nElements)
{
	int ipos = 43;
	int i;
	char t[BSIZE];
	int* elements = NULL;
	*nElements = 0;
	 while(!feof(file))
	 {
		if(!fgets(t,BSIZE,file))break;
		if(strstr( t, blockName))
		{
			if(!(strstr( t, blockName) && strstr(t,"N=") && strlen(strstr(t,"N="))>2)) return elements;
			if(strlen(t)>ipos+1 && t[ipos]!='I') return elements;
			*nElements = atof(strstr(t,"N=")+2);
			if(*nElements<1) return elements;
			elements = malloc(*nElements*sizeof(int));
			for(i=0;i<*nElements;i++)
			{
				if(1!=fscanf(file,"%d",&elements[i])) break;
			}
			if(i!=*nElements)
			{
				*nElements = 0;
				free(elements);
				return NULL;
			}
			return elements;
		}
	 }
	 return elements;
}
/****************************************************************************/
double* get_array_real_from_fchk_gaussian_file(FILE* file, char* blockName, int* nElements)
{
	int ipos = 43;
	int i;
	char t[BSIZE];
	double* elements = NULL;

	*nElements = 0;
	 while(!feof(file))
	 {
		if(!fgets(t,BSIZE,file))break;
		if(strstr( t, blockName))
		{
			if(!(strstr( t, blockName) && strstr(t,"N=") && strlen(strstr(t,"N="))>2)) return elements;
			if(strlen(t)>ipos+1 && t[ipos]!='R') return elements;
			*nElements = atof(strstr(t,"N=")+2);
			if(*nElements<1) return elements;
			elements = malloc(*nElements*sizeof(double));
			for(i=0;i<*nElements;i++)
				if(1!=fscanf(file,"%lf",&elements[i])) break;
			if(i!=*nElements)
			{
				*nElements = 0;
				free(elements);
				return NULL;
			}
			return elements;
		}
	 }
	 return elements;
}
/****************************************************************************/
char** get_array_string_from_fchk_gaussian_file(FILE* file, char* blockName, int* nElements)
{
	int ipos = 43;
	int i;
	char t[BSIZE];
	char** elements = NULL;
	char type = ' ';

	*nElements = 0;
	 while(!feof(file))
	 {
		if(!fgets(t,BSIZE,file))break;
		if(strstr( t, blockName))
		{
			if(!(strstr( t, blockName) && strstr(t,"N=") && strlen(strstr(t,"N="))>2)) return elements;
			if(strlen(t)>ipos+1 && t[ipos]=='C') type = 'C';
			if(strlen(t)>ipos+1 && t[ipos]=='H') type = 'H';
			if(type!='C' && type!='H') return elements;
			*nElements = atof(strstr(t,"N=")+2);
			if(*nElements<1) return elements;
			elements = malloc(*nElements*sizeof(char*));
			for(i=0;i<*nElements;i++) elements[i] = NULL;
			if(type=='C')
			for(i=0;i<*nElements;i++)
			{
				if(1!=fscanf(file,"%12s",t)) break;
				elements[i] = strdup(t);
			}
			else
			for(i=0;i<*nElements;i++)
			{
				if(1!=fscanf(file,"%8s",t)) break;
				elements[i] = strdup(t);
			}
			
			if(i!=*nElements)
			{
				*nElements = 0;
				free(elements);
				return NULL;
			}
			return elements;
		}
	 }
	 return elements;
}
/****************************************************************************************************************************/
static int getNumAxis(char* s)
{
        if(!s) return 0;
        if(!strcmp(s,"x")) return 1;
        if(!strcmp(s,"X")) return 1;
        if(!strcmp(s,"y")) return 2;
        if(!strcmp(s,"Y")) return 2;
        if(!strcmp(s,"z")) return 3;
        if(!strcmp(s,"Z")) return 3;
        return atoi(s);
}
/****************************************************************************************************************************/
static boolean get2IndexValue(char* t, int* i, int* j, double*value)
{
	char s1[BSIZE];
	char s2[BSIZE];
	*i = 0;
	*j = 0;
	*value = 0;
//	if(mystrcasestr(t,"END")){ printf("%s\n",t);return FALSE;}
	if(mystrcasestr(t,"END"))return FALSE;
	if(3!=sscanf(t,"%s %s %lf",s1, s2,value)) return FALSE;
	*i = getNumAxis(s1);
	*j = getNumAxis(s2);
	return TRUE;
}
/****************************************************************************************************************************/
static boolean get3IndexValue(char* t, int* i, int* j, int* k, double*value)
{
	char s1[BSIZE];
	char s2[BSIZE];
	char s3[BSIZE];
	*i = 0;
	*j = 0;
	*k = 0;
	*value = 0;
	//if(mystrcasestr(t,"END")){ printf("%s\n",t);return FALSE;}
	if(mystrcasestr(t,"END"))return FALSE;
	if(4!=sscanf(t,"%s %s %s %lf",s1, s2,s3, value)) return FALSE;
	*i = getNumAxis(s1);
	*j = getNumAxis(s2);
	*k = getNumAxis(s3);
	return TRUE;
}
/****************************************************************************************************************************/
static boolean get4IndexValue(char* t, int* i, int* j, int* k, int* l, double*value)
{
	char s1[BSIZE];
	char s2[BSIZE];
	char s3[BSIZE];
	char s4[BSIZE];
	*i = 0;
	*j = 0;
	*k = 0;
	*l = 0;
	*value = 0;
	//if(mystrcasestr(t,"END")){ printf("%s\n",t);return FALSE;}
	if(mystrcasestr(t,"END"))return FALSE;
	if(5!=sscanf(t,"%s %s %s %s %lf",s1, s2,s3, s4, value)) return FALSE;
	*i = getNumAxis(s1);
	*j = getNumAxis(s2);
	*k = getNumAxis(s3);
	*l = getNumAxis(s4);
	return TRUE;
}
/****************************************************************************************************************************/
boolean readMatrixReal(FILE* file, char* tag, int nrows, int ncolumns, double**values)
{
	static char *t = NULL; 
	char* TAG = NULL;
	char* pos = NULL;
	int i=0;
	int ii,jj;
	double v;
	int** counter = NULL;
	if(!tag) return FALSE;
	if(!values) return FALSE;
	if(t==NULL) t = malloc(BSIZE*sizeof(char));

	TAG = strdup(tag);
	uppercase(TAG);
	rewind(file);

	while(!feof(file))
  	{
    		if(!fgets(t,BSIZE, file)) break;
		deleteFirstSpaces(t);
		if(t[0]=='#') continue;
		uppercase(t);
		pos = strstr(t,TAG);
		if(!pos) continue;
		break;
	}
	if(!pos) 
	{
		free(TAG);
		return FALSE;
	}
	counter = newMatrixInt(nrows,ncolumns);
	//fprintf(stderr,"nrows=%d ncolums=%d\n",nrows,ncolumns);
	if(!counter) 
	{
		fprintf(stderr,"I am not able to create a %d %d matrix\n",nrows,ncolumns);
		exit(1);
	}
	initMatrixInt(counter, nrows,ncolumns, 0);
	while(!feof(file))
	{
    		if(!fgets(t,BSIZE, file)) break;
		if(get2IndexValue(t,&ii, &jj, &v))
		{
			if(ii<=nrows && ii>0 && jj<=ncolumns && jj>0 ) 
			{ 
				int c=  counter[ii-1][jj-1];  
				counter[ii-1][jj-1]++; 
				i++; 
				//c = 0; // TO DO 
				if(c==0) values[ii-1][jj-1] = v; 
				else values[ii-1][jj-1] = (values[ii-1][jj-1]*c+v)/(c+1);
				//if(c!=0) fprintf(stderr," %d %d v=%f newV=%f c = %d\n",ii,jj,v, values[ii-1][jj-1],c);
			}
			else
			{
				fprintf(stderr,"Erreur dans les donnees de %s: verifie les indices\n",tag);
				fprintf(stderr,"ii=%d jj=%d\n",ii,jj);
				fprintf(stderr,"nr=%d nc=%d\n",nrows,ncolumns);
				exit(1);
			}
		}
		else break;
	}
	free(TAG);
	freeMatrixInt(&counter,nrows);
	if(i<=0) return FALSE;
	return TRUE;
}
/****************************************************************************************************************************/
boolean readCubeReal(FILE* file, char* tag, int nrows, int ncolumns, int nslices, double***values)
{
	static char *t = NULL; 
	char* TAG = NULL;
	char* pos = NULL;
	int i=0;
	int ii,jj,kk;
	double v;
	int*** counter = NULL;
	if(!tag) return FALSE;
	if(!values) return FALSE;
	if(t==NULL) t = malloc(BSIZE*sizeof(char));

	TAG = strdup(tag);
	uppercase(TAG);
	rewind(file);

	while(!feof(file))
  	{
    		if(!fgets(t,BSIZE, file)) break;
		deleteFirstSpaces(t);
		if(t[0]=='#') continue;
		uppercase(t);
		pos = strstr(t,TAG);
		if(!pos) continue;
		break;
	}
	if(!pos) 
	{
		free(TAG);
		return FALSE;
	}
	counter = newCubeInt(nrows,ncolumns,nslices);
	initCubeInt(counter, nrows,ncolumns,nslices, 0);
	while(!feof(file))
	{
    		if(!fgets(t,BSIZE, file)) break;
		if(get3IndexValue(t,&ii, &jj, &kk, &v))
		{
			if(ii<=nrows && ii>0 && jj<=ncolumns && jj>0  &&kk<=nslices && kk>0 ) 
			{ 
				int c=  counter[ii-1][jj-1][kk-1];  
				counter[ii-1][jj-1][kk-1]++; 
				i++; 
				//c = 0;//TO DO
				if(c==0) values[ii-1][jj-1][kk-1] = v; 
				else values[ii-1][jj-1][kk-1] = (values[ii-1][jj-1][kk-1]*c+v)/(c+1);
				//if(c!=0) fprintf(stderr," %d %d %d v=%f newV=%f c = %d\n",ii,jj,kk,v, values[ii-1][jj-1][kk-1],c);
			}
			else
			{
				fprintf(stderr,"Erreur dans les donnees de %s: verifie les indices\n",tag);
				exit(1);
			}
		}
		else break;
	}
	free(TAG);
	freeCubeInt(&counter,nrows,ncolumns);
	if(i<=0) return FALSE;
	return TRUE;
}
/****************************************************************************************************************************/
boolean readQuarticReal(FILE* file, char* tag, int nrows, int ncolumns, int nslices, int nq, double****values)
{
	static char *t = NULL; 
	char* TAG = NULL;
	char* pos = NULL;
	int i=0;
	int ii,jj,kk,ll;
	double v;
	int**** counter = NULL;
	if(!tag) return FALSE;
	if(!values) return FALSE;
	if(t==NULL) t = malloc(BSIZE*sizeof(char));

	TAG = strdup(tag);
	uppercase(TAG);
	rewind(file);

	while(!feof(file))
  	{
    		if(!fgets(t,BSIZE, file)) break;
		deleteFirstSpaces(t);
		if(t[0]=='#') continue;
		uppercase(t);
		pos = strstr(t,TAG);
		if(!pos) continue;
		break;
	}
	if(!pos) 
	{
		free(TAG);
		return FALSE;
	}
	counter = newQuarticInt(nrows,ncolumns,nslices,nq);
	initQuarticInt(counter, nrows,ncolumns,nslices,nq, 0);
	while(!feof(file))
	{
    		if(!fgets(t,BSIZE, file)) break;
		if(get4IndexValue(t,&ii, &jj, &kk, &ll, &v))
		{
			if(ii<=nrows && ii>0 && jj<=ncolumns && jj>0  &&kk<=nslices && kk>0 && ll<=nq && ll>0 )
			{ 
				int c=  counter[ii-1][jj-1][kk-1][ll-1];  
				counter[ii-1][jj-1][kk-1][ll-1]++; 
				i++; 
                                //c = 0;// TO DO
				if(c==0) values[ii-1][jj-1][kk-1][ll-1] = v; 
				else values[ii-1][jj-1][kk-1][ll-1] = (values[ii-1][jj-1][kk-1][ll-1]*c+v)/(c+1);
				//if(c!=0) fprintf(stderr," %d %d %d %d v=%f newV=%f c = %d\n",ii,jj,kk,ll,v, values[ii-1][jj-1][kk-1][ll-1],c);
			}
			else
			{
				fprintf(stderr,"Erreur dans les donnees de %s: verifie les indices\n",tag);
				exit(1);
			}
		}
		else break;
	}
	free(TAG);
	freeQuarticInt(&counter,nrows,ncolumns,nslices);
	if(i<=0) return FALSE;
	return TRUE;
}
/****************************************************************************************************************************/
double getMaxQuarticIJKL(double**** C, int i, int j, int k, int l)
{
	double A[] = { getMaxCubicIJK(C[i],j,k,l), getMaxCubicIJK(C[j],i,k,l), getMaxCubicIJK(C[k],j,i,l), getMaxCubicIJK(C[l],j,k,i)};
	double max=A[0];
	int c;
	for(c=1;c<4;c++) if(fabs(max)<fabs(A[c])) max = A[c];
	return max;
}
/****************************************************************************************************************************/
double getMaxCubicIJK(double*** C, int i, int j, int k)
{
	double A[] = {getMaxMatrixIJ(C[i],j,k), getMaxMatrixIJ(C[j],i,k), getMaxMatrixIJ(C[k],i,j)};
	double max=A[0];
	int c;

	for(c=1;c<3;c++) if(fabs(max)<fabs(A[c])) max = A[c];
	return max;
}
/****************************************************************************************************************************/
double getMaxMatrixIJ(double** C, int i, int j)
{
	double A = C[i][j];
	if(fabs(C[j][i])> fabs(A)) A = C[j][i];
	return A;
}
/****************************************************************************************************************************/
double getMaxCubicIJ(double*** C, int a, int i, int j)
{
	double A = getMaxMatrixIJ(C[a],i,j);
	return A;
}
/****************************************************************************************************************************/
double getMaxQuarticIJK(double**** C, int a, int i, int j, int k)
{
	double A = getMaxCubicIJK(C[a],i,j,k);
	return A;
}
/*************************************************************************************/
void getCrossProduct(double * V1,double * V2, double *PV)
{
	PV[0]= V1[1]*V2[2]-V1[2]*V2[1];
	PV[1]= V1[2]*V2[0]-V1[0]*V2[2];
	PV[2]= V1[0]*V2[1]-V1[1]*V2[0];
} 
/*************************************************************************************/
double getDistancePoints(double* P1,double* P2)
{
	double distance = 0.0;
	int i;
	for(i=0;i<3;i++) distance += (P1[i]- P2[i])*(P1[i]- P2[i]);
	distance = sqrt(distance);
	return distance;
} 
/*************************************************************************************/
double getNorm(double* V)
{
	double norm = 0.0;
	int i;
	for(i=0;i<3;i++) norm += V[i]*V[i];
	return sqrt(norm);
}
/*************************************************************************************/
double getDotProduct(double* V1,double* V2)
{
	double dot = 0.0;
	int i;
	for(i=0;i<3;i++) dot += V1[i]*V2[i];
	return dot;
}
/*************************************************************************************/
double getAngleVectors(double* V1,double* V2)
{
	double angle;
	double modv1v2 = getNorm(V1)*getNorm(V2);
 
	if(fabs(modv1v2)>1e-14 )
	{
		angle = getDotProduct(V1,V2)/modv1v2;
		if(angle<=-1) return 180.0;
		if(angle>=1) return 0.0;
		angle = acos(angle)/DEGTORAD;
	}
	else angle = 0.0;
	return angle;
} 
//  Windows
#ifdef _WIN32
#include <Windows.h>
double get_wall_time(){
    LARGE_INTEGER time,freq;
    if (!QueryPerformanceFrequency(&freq)){
        //  Handle error
        return 0;
    }
    if (!QueryPerformanceCounter(&time)){
        //  Handle error
        return 0;
    }
    return (double)time.QuadPart / freq.QuadPart;
}
double get_cpu_time(){
    FILETIME a,b,c,d;
    if (GetProcessTimes(GetCurrentProcess(),&a,&b,&c,&d) != 0){
        //  Returns total user time.
        //  Can be tweaked to include kernel times as well.
        return
            (double)(d.dwLowDateTime |
            ((unsigned long long)d.dwHighDateTime << 32)) * 0.0000001;
    }else{
        //  Handle error
        return 0;
    }
}

//  Posix/Linux
#else
#include <time.h>
#include <sys/time.h>
double get_wall_time(){
    struct timeval time;
    if (gettimeofday(&time,NULL)){
        //  Handle error
        return 0;
    }
    return (double)time.tv_sec + (double)time.tv_usec * .000001;
}
double get_cpu_time(){
    return (double)clock() / CLOCKS_PER_SEC;
}
#endif

void get3DRandMatrix(double M[3][3])
{
	const double pi = M_PI;
	int i,j;
	double RotMatrix[3][3];
	double x1,x2,x3;
	x1 = drandom();
	x2 = drandom();
	x3 = drandom();

	RotMatrix[0][0] = cos(2.0*pi*x1);
	RotMatrix[0][1] = sin(2.0*pi*x1);
	RotMatrix[0][2] = 0.0;

	RotMatrix[1][0] = -1.0 * sin(2.0*pi*x1);
	RotMatrix[1][1] = cos(2.0*pi*x1);
	RotMatrix[1][2] = 0.0;

	RotMatrix[2][0] = 0.0;
	RotMatrix[2][1] = 0.0;
	RotMatrix[2][2] = 1.0;

	double v[3] = {   
		sqrt(x3)*cos(2.0*pi*x2),
		sqrt(x3)*sin(2.0*pi*x2), 
                sqrt(1.0 - x3) 
		};

	double vv[3][3];
	vv[0][0] = v[0] * v[0];
	vv[0][1] = v[0] * v[1];
	vv[0][2] = v[0] * v[2];

	vv[1][0] = v[1] * v[0];
	vv[1][1] = v[1] * v[1];
	vv[1][2] = v[1] * v[2];

	vv[2][0] = v[2] * v[0];
	vv[2][1] = v[2] * v[1];
	vv[2][2] = v[2] * v[2];

	double H[3][3];

	H[0][0] = 1.0 - 2*vv[0][0];
	H[0][1] = 0.0 - 2*vv[0][1];
	H[0][2] = 0.0 - 2*vv[0][2];
	
	H[1][0] = 0.0 - 2*vv[1][0];
	H[1][1] = 1.0 - 2*vv[1][1];
	H[1][2] = 0.0 - 2*vv[1][2];
	
	H[2][0] = 0.0 - 2*vv[2][0];
	H[2][1] = 0.0 - 2*vv[2][1];
	H[2][2] = 1.0 - 2*vv[2][2];

	for (i = 0; i < 3; i++)
	for (j = 0; j < 3; j++)
	{
		int k;
		M[i][j] = 0.0;
		for(k = 0; k < 3; ++k) M[i][j] += H[i][k] * RotMatrix[k][j];
		M[i][j] = -M[i][j];
	}	
}	
void getRandDirection(double D[])
{
	const double pi = M_PI;
	double x1,x2;
	x1 = drandom();
	x2 = drandom();

	D[2] = cos(pi*x1);
	D[1] = sin(pi*x1)*cos(2.0*pi*x2);
	D[0] = sin(pi*x1)*sin(2.0*pi*x2);
}	
