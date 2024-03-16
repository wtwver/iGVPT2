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

/* MolecularMechanicsDlg.c */
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#ifdef ENABLE_MPI
#include <mpi.h>
#endif

#include "../Utils/Utils.h"
#include "../Utils/AtomsProp.h"
#include "../Utils/Constants.h"
#include "../Utils/Types.h"
#include "../Molecule/Molecule.h"
#include "../MolecularMechanics/ForceField.h"
#include "../MolecularMechanics/MolecularMechanics.h"
#include "../MolecularMechanics/ConjugateGradient.h"
#include "../MolecularMechanics/SteepestDescent.h"
#include "../MolecularMechanics/QuasiNewton.h"
#include "../MolecularMechanics/MolecularDynamics.h"

typedef enum
{
	MMBOND = 0,
	MMBEND = 1,
	MMTORSION = 2,
	MMIMPROPER = 3,
	MMNONBOND = 4,
	MMHBOND  =5 ,
	MMCOULOMB = 6,
	PWCOULOMB = 7,
	PWVANDERWALS = 8
} MMOptions;

typedef enum
{
	GRADQUASINEWTON  = 0,
	GRADSTEEPEST  = 1,
	GRADCONJUGATE = 2,
	GRADHESTENES  = 3,
	GRADFLETCHER  = 4,
	GRADPOLAK     = 5,
	GRADWOLF      = 6
} GradientOptions;

typedef enum
{
	GRADMAXITERATIONS  = 0,
	GRADEPSILON        = 1,
	GRADMAXLINES       = 2,
	GRADINITIALSTEP    = 3,
	GRADFREQUENCY      = 4
} GradientEntrys;

typedef enum
{
	TOLE = 0,
	TOLD = 1
} TOLptions;

#define NGRADENTRYS 5
#define NGRADOPTIONS 7
#define NOPTIONS1 4
#define NOPTIONS2 3
#define NOPTIONS3 2
#define NINTEGOPTIONS 3
#define NTHERMOPTIONS 3
#define NENTRYTOL 2
#define NCONSTRAINTS 3

/*****************************************************************************/
static void checkWallCorrection(FILE* file, ForceField* forceField)
{
	char t[BSIZE];
	char* pos;
	rewind(file);
	forceField->options.addWallCorrection=FALSE;
	while(!feof(file))
  	{
    		if(!fgets(t,BSIZE, file)) break;
		deleteFirstSpaces(t);
		if(t[0]=='#') continue;
		uppercase(t);
		pos = strstr(t,"WALL");
		if(pos && strstr(t,"="))
		{ 
			double E0;
			double rho;
			int nc;
			int n =0;
			pos = strstr(t,"=") + 1;
			n =sscanf(pos,"%lf %lf %d",&E0,&rho,&nc);
			//printf("t=%s\n",t);
			//printf("pos=%s\n",pos);
			if(n==3 && nc%2==0) { forceField->options.addWallCorrection=TRUE; break;}
			else  { 
				fprintf(stderr,"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n");
				fprintf(stderr,"Error during the reading of Wall parameters\n");
				fprintf(stderr,"You must give E0(au), rho (cutoff radius in angstrom) and nc(even integer)\n");
				fprintf(stderr,"Example : Wall=1000.0 10.0 6\n");
				fprintf(stderr,"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n");
				exit(1);
			}
			break;
		}
	}
}
/*****************************************************************************/
static void printEnergyAndGradient(ForceField* forceField)
{
	char* str;
	double gradientNorm = 0;
	int i,j;

	//forceField->klass->calculateEnergy(forceField);
	forceField->klass->calculateGradient(forceField);

	gradientNorm = 0;
	for (  i = 0; i < forceField->molecule.nAtoms; i++ )
		for(j=0;j<3;j++)
			gradientNorm += 
			forceField->molecule.atoms[i].gradient[j]
			*forceField->molecule.atoms[i].gradient[j]; 

	str = strdup_printf(("Gradient Norm  = %0.14f energy = %0.14f(kcal/mol)\n"),
		sqrt(gradientNorm),forceField->molecule.potentialEnergy); 

	printf("%s",str);
	free(str);
	forceField->klass->printEnergies(forceField);
}
/*********************************************************************************/
static void getMultiplicityName(int multiplicity, char* buffer)
{
	if(multiplicity==1) sprintf(buffer,"Singlet");
	else if(multiplicity==2) sprintf(buffer,"Doublet");
	else if(multiplicity==3) sprintf(buffer,"Triplet");
	else if(multiplicity==4) sprintf(buffer,"Quartet");
	else if(multiplicity==5) sprintf(buffer,"Quintet");
	else if(multiplicity==6) sprintf(buffer,"Sextet");
	else sprintf(buffer,"UNKNOWN");
}
/*****************************************************************************/
static boolean createCChemIFiles(int numberOfGeometries, ForceField** geometries, double* energies, char* fileNamePrefix, char* keyWords,char* cchemiCommand)
{
	FILE* file = NULL;
	FILE* fileSH = NULL;
	int i;
	int j;
	int nG = 0;
	char* fileName = NULL;
	char* fileNameSH = NULL;

	if(numberOfGeometries<1) return FALSE;
	if(!geometries) return FALSE;
	if(!energies) return FALSE;
	for(i=0;i<numberOfGeometries;i++) if(geometries[i]) nG++;
	if(nG<1) return FALSE;
	fileNameSH = strdup_printf("%scchemi.sh",fileNamePrefix);
#ifndef OS_WIN32
	fileNameSH = strdup_printf("%scchemi.sh",fileNamePrefix);
#else
	fileNameSH = strdup_printf("%sCChemI.bat",fileNamePrefix);
#endif
 	fileSH = fopen(fileNameSH, "w");
	if(!fileSH) return FALSE;


	for(i=0;i<numberOfGeometries;i++)
	{
		if(!geometries[i]) continue;
 		if(fileName) free(fileName);
		fileName = strdup_printf("%sCCHEMI_%d.inp",fileNamePrefix,i+1);
 		file = fopen(fileName, "w");
		if(!file) return FALSE;
		fprintf(file,"# RunType = Energy, Optimization, MD, MDConfo, REMDConfo\n");
		fprintf(file,"RunType=Optimization\n");
		fprintf(file,"OptimizerType=QuasiNewton\n");
		fprintf(file,"Model=Mopac\n");
		fprintf(file,"SEKeys=AM1\n");
		fprintf(file,"mopacCommand=/home/allouche/Softwares/MOPAC2009/MOPAC2009.exe\n");
		fprintf(file,"orcaCommand=orca\n");
		fprintf(file,"fireflyCommand=firefly\n");
		fprintf(file,"gaussianCommand=g03\n");
		fprintf(file,"#gaussianKeywordsPost=B3LYP/aug-cc-pvdz\n");
		fprintf(file,"#QuasiNewton\n");
		fprintf(file,"useQuasiNewton = TRUE\n");
		fprintf(file,"Geometry\n");
		fprintf(file,"%d %d %d\n",geometries[i]->molecule.nAtoms, geometries[i]->molecule.totalCharge, geometries[i]->molecule.spinMultiplicity);
		for(j=0;j<geometries[i]->molecule.nAtoms;j++)
		{
			int nc = 0;
			int k;
			for(k=0;k<geometries[i]->molecule.nAtoms;k++) 
				if(geometries[i]->molecule.atoms[j].typeConnections&&geometries[i]->molecule.atoms[j].typeConnections[k]>0) nc++;

			fprintf(file," %s %s %s %s %d %f %d %d %f %f %f %d ", 
				geometries[i]->molecule.atoms[j].prop.symbol,
				geometries[i]->molecule.atoms[j].mmType,
				geometries[i]->molecule.atoms[j].pdbType,
				geometries[i]->molecule.atoms[j].residueName,
				geometries[i]->molecule.atoms[j].residueNumber,
				geometries[i]->molecule.atoms[j].charge,
				geometries[i]->molecule.atoms[j].layer,
				geometries[i]->molecule.atoms[j].variable,
				geometries[i]->molecule.atoms[j].coordinates[0],
				geometries[i]->molecule.atoms[j].coordinates[1],
				geometries[i]->molecule.atoms[j].coordinates[2],
				nc
				);
			for(k=0;k< geometries[i]->molecule.nAtoms;k++) 
			{
		 		int nk =  geometries[i]->molecule.atoms[k].N-1;
				if(geometries[i]->molecule.atoms[j].typeConnections && geometries[i]->molecule.atoms[j].typeConnections[nk]>0) 
					fprintf(file," %d %d", nk+1, geometries[i]->molecule.atoms[j].typeConnections[nk]);
			}
			fprintf(file,"\n");
		}
		fprintf(file,"Velocities\n");
		for(j=0;j<geometries[i]->molecule.nAtoms;j++)
		{
			fprintf(file,"%f %f %f", 
				geometries[i]->molecule.atoms[j].velocity[0],
				geometries[i]->molecule.atoms[j].velocity[1],
				geometries[i]->molecule.atoms[j].velocity[2]
				);
			fprintf(file,"\n");
		}
		fprintf(file,"Masses\n");
		for(j=0;j<geometries[i]->molecule.nAtoms;j++)
		{
			fprintf(file,"%f", geometries[i]->molecule.atoms[j].mass);
			fprintf(file,"\n");
		}
		fclose(file);
		fprintf(fileSH,"%s %s\n",cchemiCommand,fileName);
	}
	fclose(fileSH);
#ifndef OS_WIN32
	{
		char buffer[1024];
  		sprintf(buffer,"chmod u+x %s",fileNameSH);
		system(buffer);
	}
#endif
 	if(fileName) free(fileName);
 	if(fileNameSH) free(fileNameSH);
	return TRUE;

}
/*****************************************************************************/
static boolean createMopacFiles(int numberOfGeometries, ForceField** geometries, double* energies, char* fileNamePrefix, char* keyWords, char* mopacCommand)
{
	FILE* file = NULL;
	FILE* fileSH = NULL;
	int i;
	int j;
	int nG = 0;
	char* fileName = NULL;
	char* fileNameSH = NULL;
	char multiplicityStr[100];
#ifdef OS_WIN32
	char c='%';
#endif

	if(numberOfGeometries<1) return FALSE;
	if(!geometries) return FALSE;
	if(!energies) return FALSE;
	for(i=0;i<numberOfGeometries;i++) if(geometries[i]) nG++;
	if(nG<1) return FALSE;
#ifndef OS_WIN32
	fileNameSH = strdup_printf("%sMopac.sh",fileNamePrefix);
#else
	fileNameSH = strdup_printf("%sMopac.bat",fileNamePrefix);
#endif
 	fileSH = fopen(fileNameSH, "w");
	if(!fileSH) return FALSE;
#ifdef OS_WIN32
	fprintf(fileSH,"set PATH=%cPATH%c;\"%s\"\n",c,c,mopacDirectory);
#endif


	for(i=0;i<numberOfGeometries;i++)
	{
		getMultiplicityName( geometries[i]->molecule.spinMultiplicity, multiplicityStr);
		if(!geometries[i]) continue;
 		if(fileName) free(fileName);
		fileName = strdup_printf("%s_%d.mop",fileNamePrefix,i+1);
 		file = fopen(fileName, "w");
		if(!file) return FALSE;
		fprintf(file,"* ===============================\n");
		fprintf(file,"* Input file for Mopac\n");
		fprintf(file,"* MM/SE Energy(kCal/mol) =%f\n",energies[i]);
		fprintf(file,"* ===============================\n");
		fprintf(file,"%s CHARGE=%d %s\n",keyWords, geometries[i]->molecule.totalCharge,multiplicityStr);
		fprintf(file,"\n");
		fprintf(file,"Mopac file generated by CChemI\n");

		for(j=0;j<geometries[i]->molecule.nAtoms;j++)
		{
		fprintf(file," %s %f %d %f %d %f %d\n", 
				geometries[i]->molecule.atoms[j].prop.symbol,
				geometries[i]->molecule.atoms[j].coordinates[0],
				1,
				geometries[i]->molecule.atoms[j].coordinates[1],
				1,
				geometries[i]->molecule.atoms[j].coordinates[2],
				1
				);
		}
		fclose(file);
		fprintf(fileSH,"%s %s\n",mopacCommand,fileName);
	}
	fclose(fileSH);
#ifndef OS_WIN32
	{
		char buffer[1024];
  		sprintf(buffer,"chmod u+x %s",fileNameSH);
		system(buffer);
	}
#endif
 	if(fileName) free(fileName);
 	if(fileNameSH) free(fileNameSH);
	return TRUE;

}
/*****************************************************************************/
static boolean createGaussianFiles(int numberOfGeometries, ForceField** geometries, double* energies, char* fileNamePrefix, char* keyWords,char* gaussianCommand)
{
	FILE* file = NULL;
	FILE* fileSH = NULL;
	int i;
	int j;
	int nG = 0;
	char* fileName = NULL;
	char* fileNameSH = NULL;

	if(numberOfGeometries<1) return FALSE;
	if(!geometries) return FALSE;
	if(!energies) return FALSE;
	for(i=0;i<numberOfGeometries;i++) if(geometries[i]) nG++;
	if(nG<1) return FALSE;
	fileNameSH = strdup_printf("%sGauss.sh",fileNamePrefix);
#ifndef OS_WIN32
	fileNameSH = strdup_printf("%sGauss.sh",fileNamePrefix);
#else
	fileNameSH = strdup_printf("%sGauss.bat",fileNamePrefix);
#endif
 	fileSH = fopen(fileNameSH, "w");
	if(!fileSH) return FALSE;


	for(i=0;i<numberOfGeometries;i++)
	{
		if(!geometries[i]) continue;
 		if(fileName) free(fileName);
		fileName = strdup_printf("%s_%d.com",fileNamePrefix,i+1);
 		file = fopen(fileName, "w");
		if(!file) return FALSE;
		fprintf(file,"#P %s\n",keyWords);
		fprintf(file,"#  Units(Ang,Deg)\n");
		fprintf(file,"\n");
		fprintf(file,"File generated by CChemI\n");
		fprintf(file,"MM/SE Energy(kCal/mol) = %f\n",energies[i]);
		fprintf(file,"\n");
		fprintf(file,"%d %d\n", geometries[i]->molecule.totalCharge, geometries[i]->molecule.spinMultiplicity);
		for(j=0;j<geometries[i]->molecule.nAtoms;j++)
		{
		fprintf(file,"%s %f %f %f\n", 
				geometries[i]->molecule.atoms[j].prop.symbol,
				geometries[i]->molecule.atoms[j].coordinates[0],
				geometries[i]->molecule.atoms[j].coordinates[1],
				geometries[i]->molecule.atoms[j].coordinates[2]
				);
		}
		fprintf(file,"\n");
		fclose(file);
		fprintf(fileSH,"%s %s\n",gaussianCommand,fileName);
	}
	fclose(fileSH);
#ifndef OS_WIN32
	{
		char buffer[1024];
  		sprintf(buffer,"chmod u+x %s",fileNameSH);
		system(buffer);
	}
#endif
 	if(fileName) free(fileName);
 	if(fileNameSH) free(fileNameSH);
	return TRUE;

}
/*****************************************************************************/
static boolean createFireFlyFiles(int numberOfGeometries, ForceField** geometries, double* energies, char* fileNamePrefix, char* keyWords, char* fireflyCommand)
{
	FILE* file = NULL;
	FILE* fileSH = NULL;
	int i;
	int j;
	int nG = 0;
	char* fileName = NULL;
	char* fileNameSH = NULL;
	char buffer[1024];
#ifdef OS_WIN32
	char c='%';
#endif

	if(numberOfGeometries<1) return FALSE;
	if(!geometries) return FALSE;
	if(!energies) return FALSE;
	for(i=0;i<numberOfGeometries;i++) if(geometries[i]) nG++;
	if(nG<1) return FALSE;
#ifndef OS_WIN32
	fileNameSH = strdup_printf("%sPCGam.sh",fileNamePrefix);
#else
	fileNameSH = strdup_printf("%sPCGam.bat",fileNamePrefix);
#endif
 	fileSH = fopen(fileNameSH, "w");
	if(!fileSH) return FALSE;
#ifdef OS_WIN32
	fprintf(fileSH,"set PATH=%cPATH%c;\"%s\"\n",c,c,fireflyDirectory);
#endif


	uppercase(keyWords);
	for(i=0;i<numberOfGeometries;i++)
	{
		if(!geometries[i]) continue;
 		if(fileName) free(fileName);
		fileName = strdup_printf("%sFF_%d.inp",fileNamePrefix,i+1);
 		file = fopen(fileName, "w");
		if(!file) return FALSE;
		fprintf(file,"! ======================================================\n");
		fprintf(file,"!  Input file for FireFly\n"); 
		fprintf(file,"! ======================================================\n");
		if(strstr(keyWords,"RUNTYP"))
		{
			sscanf(strstr(keyWords,"RUNTYP"),"%s",buffer);
			fprintf(file," $CONTRL %s $END\n",buffer);
		}
		if(strstr(keyWords,"SCFTYP"))
		{
			sscanf(strstr(keyWords,"SCFTYP"),"%s",buffer);
			fprintf(file," $CONTRL %s $END\n",buffer);
		}
		else
		{
			if(geometries[i]->molecule.spinMultiplicity==1)
				fprintf(file," $CONTRL SCFTYP=RHF $END\n");
			else
				fprintf(file," $CONTRL SCFTYP=UHF $END\n");
		}

		fprintf(file," $CONTRL ICHARG=%d MULT=%d $END\n", geometries[i]->molecule.totalCharge, geometries[i]->molecule.spinMultiplicity);
		if(strstr(keyWords,"GBASIS"))
		{
			sscanf(strstr(keyWords,"GBASIS"),"%s",buffer);
			fprintf(file," $BASIS %s $END\n",buffer);
		}
		fprintf(file," $DATA\n");
		fprintf(file,"Molecule specification\n");
		fprintf(file,"C1\n");
		for(j=0;j<geometries[i]->molecule.nAtoms;j++)
		{
			char* symbol = geometries[i]->molecule.atoms[j].prop.symbol;
			SAtomsProp prop = propAtomGet(symbol);
			fprintf(file,"%s %f %f %f %f\n", 
				symbol,
				(double)prop.atomicNumber,
				geometries[i]->molecule.atoms[j].coordinates[0],
				geometries[i]->molecule.atoms[j].coordinates[1],
				geometries[i]->molecule.atoms[j].coordinates[2]
				);
		}
		fprintf(file," $END\n");
		fclose(file);

#ifndef OS_WIN32
		if(!strcmp(fireflyCommand,"pcgamess") || !strcmp(fireflyCommand,"nohup pcgamess")||
		!strcmp(fireflyCommand,"firefly") || !strcmp(fireflyCommand,"nohup firefly"))
		{
			fprintf(fileSH,"mkdir %stmp%d\n",fileNamePrefix,i+1);
			fprintf(fileSH,"cd %stmp%d\n",fileNamePrefix,i+1);
			fprintf(fileSH,"cp %s input\n",fileName);
			fprintf(fileSH,"%s -p -o %sFF_%d.log\n",fireflyCommand,fileNamePrefix,i+1);
			fprintf(fileSH,"cd ..\n");
			fprintf(fileSH,"mv PUNCH  %sFF_%d.pun\n",fileNamePrefix,i+1);
			fprintf(fileSH,"/bin/rm -r  %stmp%d\n",fileNamePrefix,i+1);
		}
		else
			fprintf(fileSH,"%s %s",fireflyCommand,fileName);
#else
	 	if(!strcmp(fireflyCommand,"pcgamess") ||
	 	!strcmp(fireflyCommand,"firefly") )
		{
         		fprintf(fileSH,"mkdir %stmp%d\n",fileNamePrefix,i+1);
			addUnitDisk(fileSH, fileNamePrefix);
	 		fprintf(fileSH,"cd %stmp%d\n",fileNamePrefix,i+1);
         		fprintf(fileSH,"copy %s input\n",fileName);
         		fprintf(fileSH,"%s -p -o %sFF_%d.log\n",fireflyCommand,fileNamePrefix,i+1);
	 		fprintf(fileSH,"cd ..\n");
         		fprintf(fileSH,"move PUNCH  %sFF_%d.pun\n",fileNamePrefix,i+1);
         		fprintf(fileSH,"del /Q  %stmp%d\n",fileNamePrefix,i+1);
         		fprintf(fileSH,"rmdir  %stmp%d\n",fileNamePrefix,i+1);
		}
		else
			fprintf(fileSH,"%s %s",fireflyCommand,fileName);
#endif
	}
	fclose(fileSH);
#ifndef OS_WIN32
	{
		char buffer[1024];
  		sprintf(buffer,"chmod u+x %s",fileNameSH);
		system(buffer);
	}
#endif
 	if(fileName) free(fileName);
 	if(fileNameSH) free(fileNameSH);
	return TRUE;

}
/*****************************************************************************/
static boolean createOrcaFiles(int numberOfGeometries, ForceField** geometries, double* energies, char* fileNamePrefix, char* keyWords,char* orcaCommand)
{
	FILE* file = NULL;
	FILE* fileSH = NULL;
	int i,j;
	int nG = 0;
	int nV = 0;
	char* fileName = NULL;
	char* fileOut = NULL;
	char* fileNameSH = NULL;
	Molecule* mol = NULL;

	if(numberOfGeometries<1) return FALSE;
	if(!geometries) return FALSE;
	if(!energies) return FALSE;
	for(i=0;i<numberOfGeometries;i++) if(geometries[i]) nG++;
	if(nG<1) return FALSE;
	fileNameSH = strdup_printf("%sOrca.sh",fileNamePrefix);
#ifndef OS_WIN32
	fileNameSH = strdup_printf("%sOrca.sh",fileNamePrefix);
#else
	fileNameSH = strdup_printf("%sOrca.bat",fileNamePrefix);
#endif
 	fileSH = fopen(fileNameSH, "w");
	if(!fileSH) return FALSE;


	for(j=0;j<numberOfGeometries;j++)
	{
		if(!geometries[j]) continue;
 		if(fileName) free(fileName);
 		if(fileOut) free(fileOut);
		fileOut = strdup_printf("%sORCA_%d.out",fileNamePrefix,j+1);
		fileName = strdup_printf("%sORCA_%d.inp",fileNamePrefix,j+1);
 		file = fopen(fileName, "w");
		if(!file) return FALSE;
		fprintf(file,"! %s\n",keyWords);
		mol = &geometries[j]->molecule;
		fprintf(file,"* xyz %d   %d\n",mol->totalCharge,mol->spinMultiplicity);
		nV = 0;
      		for (i=0;i<mol->nAtoms;i++)
		{
			char X[100];
			char Y[100];
			char Z[100];
			sprintf(X,"%20.14f",mol->atoms[i].coordinates[0]);
			sprintf(Y,"%20.14f",mol->atoms[i].coordinates[1]);
			sprintf(Z,"%20.14f",mol->atoms[i].coordinates[2]);

			fprintf(file," %s  %s %s %s\n",mol->atoms[i].prop.symbol, X,Y,Z);
			if(mol->atoms[i].variable) nV+=3;
		}
		fprintf(file,"*\n");
		if(nV>0&&nV!=3*mol->nAtoms) 
		{
			fprintf(file,"%cgeom Constraints\n",'%');
      			for (i=0;i<mol->nAtoms;i++)
			{
				if(mol->atoms[i].variable)
				{
					fprintf(file,"  {C %d C}\n",i);
				}
			}
			fprintf(file," end #Constraints\n");
			fprintf(file," invertConstraints true\n");
			fprintf(file," end #geom\n");
		}
		fprintf(file,"\n");
		fclose(file);
		fprintf(fileSH,"%s %s > %s\n",orcaCommand,fileName, fileOut);
	}
	fclose(fileSH);
#ifndef OS_WIN32
	{
		char buffer[1024];
  		sprintf(buffer,"chmod u+x %s",fileNameSH);
		system(buffer);
	}
#endif
 	if(fileName) free(fileName);
 	if(fileOut) free(fileOut);
 	if(fileNameSH) free(fileNameSH);
	return TRUE;
}
/*****************************************************************************/
static boolean createConfoGabeditGeometries(int numberOfGeometries, ForceField** geometries, double* energies, char* fileNamePrefix)
{
	FILE* file = NULL;
	int i;
	int j;
	int nG = 0;
	int k;
	int form = 1;
	char* fileName = NULL;

	if(numberOfGeometries<1) return FALSE;
	if(!geometries) return FALSE;
	if(!energies) return FALSE;
	for(i=0;i<numberOfGeometries;i++) if(geometries[i]) nG++;
	if(nG<1) return FALSE;

	k = -1;
	for(i=0;i<numberOfGeometries;i++)
	{
		if(!geometries[i]) continue;

 		if(fileName) free(fileName);
		fileName = strdup_printf("%s_%d.gab",fileNamePrefix,i+1);
 		file = fopen(fileName, "w");
		if(!file) continue;
		fprintf(file,"[Gabedit Format]\n");
		fprintf(file,"\n");
		fprintf(file,"[GEOMS] %d\n",form);
		fprintf(file,"%d 2\n",nG);
		fprintf(file,"energy kcal/mol 1\n");
		fprintf(file,"deltaE K 1\n");
		fprintf(file,"Dipole Debye 1\n");

		geometries[i]->molecule.klass->computeDipole(&geometries[i]->molecule);
		if(k<0) k = i;
		fprintf(file,"%f\n",energies[i]);
		if(k>=0) fprintf(file,"%f\n",(energies[i]-energies[k])*503.21892494);
		else fprintf(file,"0\n");
		fprintf(file,"%0.14f %0.14f %0.14f\n",geometries[i]->molecule.dipole[0],geometries[i]->molecule.dipole[1],geometries[i]->molecule.dipole[2]);
		fprintf(file,"%d %d %d\n",geometries[i]->molecule.nAtoms, geometries[i]->molecule.totalCharge, geometries[i]->molecule.spinMultiplicity);
		for(j=0;j<geometries[i]->molecule.nAtoms;j++)
		{
			int nc = 0;
			int k;
			for(k=0;k<geometries[i]->molecule.nAtoms;k++) 
				if(geometries[i]->molecule.atoms[j].typeConnections&&geometries[i]->molecule.atoms[j].typeConnections[k]>0) nc++;

			fprintf(file," %s %s %s %s %d %0.6f %d %d %0.8f %0.8f %0.8f %d ", 
				geometries[i]->molecule.atoms[j].prop.symbol,
				geometries[i]->molecule.atoms[j].mmType,
				geometries[i]->molecule.atoms[j].pdbType,
				geometries[i]->molecule.atoms[j].residueName,
				geometries[i]->molecule.atoms[j].residueNumber,
				geometries[i]->molecule.atoms[j].charge,
				geometries[i]->molecule.atoms[j].layer,
				geometries[i]->molecule.atoms[j].variable,
				geometries[i]->molecule.atoms[j].coordinates[0],
				geometries[i]->molecule.atoms[j].coordinates[1],
				geometries[i]->molecule.atoms[j].coordinates[2],
				nc
				);
			for(k=0;k< geometries[i]->molecule.nAtoms;k++) 
			{
		 		int nk =  geometries[i]->molecule.atoms[k].N-1;
				if(geometries[i]->molecule.atoms[j].typeConnections && geometries[i]->molecule.atoms[j].typeConnections[nk]>0) 
					fprintf(file," %d %d", nk+1, geometries[i]->molecule.atoms[j].typeConnections[nk]);
			}
			fprintf(file,"\n");
		}
		fclose(file);
	}
	return TRUE;
}
/*****************************************************************************/
static boolean saveConfoGeometries(int numberOfGeometries, ForceField** geometries, double* energies, char* fileNameGeom)
{
	FILE* file = NULL;
	int i;
	int j;
	int nG = 0;
	int k;
	int form = 1;
	char* fileNamePrefix = getSuffixNameFile(fileNameGeom);

	if(numberOfGeometries<1) return FALSE;
	if(!geometries) return FALSE;
	if(!energies) return FALSE;
	for(i=0;i<numberOfGeometries;i++) if(geometries[i]) nG++;
	if(nG<1) return FALSE;

 	file = fopen(fileNameGeom, "w");

	if(!file) return FALSE;

	fprintf(file,"[Gabedit Format]\n");
	fprintf(file,"[GEOCONV]\n");
	fprintf(file,"energy\n");
	for(i=0;i<numberOfGeometries;i++)
		if(geometries[i]) fprintf(file,"%f\n",energies[i]);
	fprintf(file,"max-force\n");
	for(i=0;i<numberOfGeometries;i++)
		if(geometries[i]) fprintf(file,"0.0\n");
	fprintf(file,"rms-force\n");
	for(i=0;i<numberOfGeometries;i++)
		if(geometries[i]) fprintf(file,"0.0\n");

	fprintf(file,"\n");
	fprintf(file,"[GEOMETRIES]\n");
	for(i=0;i<numberOfGeometries;i++)
	{
		if(!geometries[i]) continue;
		fprintf(file,"%d\n",geometries[i]->molecule.nAtoms);
		fprintf(file,"\n");
		for(j=0;j<geometries[i]->molecule.nAtoms;j++)
		fprintf(file," %s %0.8f %0.8f %0.8f\n", 
				geometries[i]->molecule.atoms[j].prop.symbol,
				geometries[i]->molecule.atoms[j].coordinates[0],
				geometries[i]->molecule.atoms[j].coordinates[1],
				geometries[i]->molecule.atoms[j].coordinates[2]
				);
	}
	fprintf(file,"\n");
	fprintf(file,"[GEOMS] %d\n",form);
	fprintf(file,"%d 2\n",nG);
	fprintf(file,"energy kcal/mol 1\n");
	fprintf(file,"deltaE K 1\n");
	fprintf(file,"Dipole Debye 1\n");
	k = -1;
	for(i=0;i<numberOfGeometries;i++)
	{
		if(!geometries[i]) continue;
		geometries[i]->molecule.klass->computeDipole(&geometries[i]->molecule);
		if(k<0) k = i;
		fprintf(file,"%f\n",energies[i]);
		if(k>=0) fprintf(file,"%f\n",(energies[i]-energies[k])*503.21892494);
		else fprintf(file,"0\n");
		fprintf(file,"%0.14f %0.14f %0.14f\n",geometries[i]->molecule.dipole[0],geometries[i]->molecule.dipole[1],geometries[i]->molecule.dipole[2]);
		fprintf(file,"%d %d %d\n",geometries[i]->molecule.nAtoms, geometries[i]->molecule.totalCharge, geometries[i]->molecule.spinMultiplicity);
		for(j=0;j<geometries[i]->molecule.nAtoms;j++)
		{
			int nc = 0;
			int k;
			for(k=0;k<geometries[i]->molecule.nAtoms;k++) 
				if(geometries[i]->molecule.atoms[j].typeConnections&&geometries[i]->molecule.atoms[j].typeConnections[k]>0) nc++;

			fprintf(file," %s %s %s %s %d %0.6f %d %d %0.8f %0.8f %0.8f %d ", 
				geometries[i]->molecule.atoms[j].prop.symbol,
				geometries[i]->molecule.atoms[j].mmType,
				geometries[i]->molecule.atoms[j].pdbType,
				geometries[i]->molecule.atoms[j].residueName,
				geometries[i]->molecule.atoms[j].residueNumber,
				geometries[i]->molecule.atoms[j].charge,
				geometries[i]->molecule.atoms[j].layer,
				geometries[i]->molecule.atoms[j].variable,
				geometries[i]->molecule.atoms[j].coordinates[0],
				geometries[i]->molecule.atoms[j].coordinates[1],
				geometries[i]->molecule.atoms[j].coordinates[2],
				nc
				);
			for(k=0;k< geometries[i]->molecule.nAtoms;k++) 
			{
		 		int nk =  geometries[i]->molecule.atoms[k].N-1;
				if(geometries[i]->molecule.atoms[j].typeConnections && geometries[i]->molecule.atoms[j].typeConnections[nk]>0) 
					fprintf(file," %d %d", nk+1, geometries[i]->molecule.atoms[j].typeConnections[nk]);
			}
			fprintf(file,"\n");
		}
	}
	fclose(file);
	createConfoGabeditGeometries(numberOfGeometries, geometries, energies, fileNamePrefix);
	return TRUE;

}
/*****************************************************************************/
static boolean getEnergyMopac(char* fileNameOut, double* energy)
{
	FILE* file = NULL;
	char buffer[1024];
	char* pdest = NULL;

 	file = fopen(fileNameOut, "r");
	if(!file) return FALSE;
	 while(!feof(file))
	 {
		if(!fgets(buffer,BSIZE,file))break;
		pdest = strstr( buffer, " FINAL HEAT OF FORMATION");
		if(pdest) 
		{
			pdest = strstr( buffer, "=");
			if(pdest)
			{
				if(sscanf(pdest+1,"%lf",energy)==1)
				{
					fclose(file);
					return TRUE;
				}
			}
		}
	 }
	fclose(file);
	return FALSE;
}
/*****************************************************************************/
static boolean runOneMopac(ForceField* geometry, double* energy, char* fileNamePrefix, char* keyWords, char* mopacCommand)
{
	FILE* file = NULL;
	FILE* fileSH = NULL;
	int j;
	char* fileNameIn = NULL;
	char* fileNameOut = NULL;
	char* fileNameSH = NULL;
	char multiplicityStr[100];
	char buffer[1024];
#ifdef OS_WIN32
	char c='%';
#endif

	if(!geometry) return FALSE;
#ifndef OS_WIN32
	fileNameSH = strdup_printf("%sMopacOne.sh",fileNamePrefix);
#else
	fileNameSH = strdup_printf("%sMopacOne.bat",fileNamePrefix);
#endif
 	fileSH = fopen(fileNameSH, "w");
	if(!fileSH) return FALSE;
#ifdef OS_WIN32
	fprintf(fileSH,"set PATH=%cPATH%c;\"%s\"\n",c,c,mopacDirectory);
#endif

	getMultiplicityName( geometry->molecule.spinMultiplicity, multiplicityStr);

	fileNameIn = strdup_printf("%sOne.mop",fileNamePrefix);
 	file = fopen(fileNameIn, "w");
	if(!file) 
	{
 		if(fileNameIn) free(fileNameIn);
 		if(fileNameOut) free(fileNameOut);
 		if(fileNameSH) free(fileNameSH);
		return FALSE;
	}
	fprintf(file,"* ===============================\n");
	fprintf(file,"* Input file for Mopac\n");
	fprintf(file,"* ===============================\n");
	fprintf(file,"%s CHARGE=%d %s\n",keyWords, geometry->molecule.totalCharge,multiplicityStr);
	fprintf(file,"\n");
	fprintf(file,"Mopac file generated by CChemI\n");

	for(j=0;j<geometry->molecule.nAtoms;j++)
	{
	fprintf(file," %s %f %d %f %d %f %d\n", 
			geometry->molecule.atoms[j].prop.symbol,
			geometry->molecule.atoms[j].coordinates[0],
			1,
			geometry->molecule.atoms[j].coordinates[1],
			1,
			geometry->molecule.atoms[j].coordinates[2],
			1
			);
	}
	fclose(file);
#ifndef OS_WIN32
	fprintf(fileSH,"%s %s\n",mopacCommand,fileNameIn);
	fclose(fileSH);
	sprintf(buffer,"chmod u+x %s",fileNameSH);
	system(buffer);
	system(fileNameSH);
#else
	fprintf(fileSH,"\"%s\" \"%s\"\n",mopacCommand,fileNameIn);
	fclose(fileSH);
	sprintf(buffer,"\"%s\"",fileNameSH);
	system(buffer);
#endif

	fileNameOut = strdup_printf("%sOne.out",fileNamePrefix);
	if(getEnergyMopac(fileNameOut,energy))
	{
		char* str = NULL;

		geometry->molecule.klass->readGeomFromMopacOutputFile(&geometry->molecule, fileNameOut, -1);
		if(strstr(keyWords,"AM1")) str = strdup_printf("Energy by AM1/Mopac = %f", *energy);
		else str = strdup_printf("Energy by PM6/Mopac = %f", *energy);
		printf("%s\n",str);
		if(str) free(str);
	}
	else
	{
 		if(fileNameIn) free(fileNameIn);
 		if(fileNameOut) free(fileNameOut);
 		if(fileNameSH) free(fileNameSH);
		return FALSE;
	}

 	if(fileNameIn) free(fileNameIn);
 	if(fileNameOut) free(fileNameOut);
 	if(fileNameSH) free(fileNameSH);
	return TRUE;
}
/*****************************************************************************/
static boolean runMopacFiles(int numberOfGeometries, ForceField** geometries, double* energies, char* fileNamePrefix, char* keyWords, char* mopacCommand)
{
	int i;
	int nG = 0;
	int nM = 0;
	char* str = NULL;
	for(i=0;i<numberOfGeometries;i++)
	{
		if(!geometries[i]) continue;
		nG++;
		if(str) free(str);
		if(strstr(keyWords,"AM1"))
		str = strdup_printf("Minimization by AM1/Mopac of geometry n = %d... Please wait\n", i+1);
		else
		str = strdup_printf("Minimization by PM6/Mopac of geometry n = %d... Please wait\n", i+1);
		printf("%s",str);
		runOneMopac(geometries[i], &energies[i], fileNamePrefix, keyWords, mopacCommand);
		nM++;
	}
	if(str) free(str);
	if(nM==nG) return TRUE;
	return FALSE;

}
/*****************************************************************************/
static boolean getEnergyFireFly(char* fileNameOut, double* energy)
{
	FILE* file = NULL;
	char buffer[1024];
	char* pdest = NULL;
	boolean OK = FALSE;

 	file = fopen(fileNameOut, "r");
	if(!file) return FALSE;
	 while(!feof(file))
	 {
		if(!fgets(buffer,BSIZE,file))break;
		pdest = strstr( buffer, "HEAT OF FORMATION IS");
		if(pdest) 
		{
			pdest = strstr( buffer, "S");
			if(pdest)
			{
				if(sscanf(pdest+1,"%lf",energy)==1)
					OK = TRUE;
			}
		}
	 }
	fclose(file);
	return OK;
}
/*****************************************************************************/
static boolean runOneFireFly(ForceField* geometry, double* energy, char* fileNamePrefix, char* keyWords, char* fireflyCommand)
{
	FILE* file = NULL;
	FILE* fileSH = NULL;
	int j;
	char* fileNameIn = NULL;
	char* fileNameOut = NULL;
	char* fileNameSH = NULL;
	char multiplicityStr[100];
	char buffer[1024];
#ifdef OS_WIN32
	char c='%';
#endif

	if(!geometry) return FALSE;
#ifndef OS_WIN32
	fileNameSH = strdup_printf("%sPCGOne.sh",fileNamePrefix);
#else
	fileNameSH = strdup_printf("%sPCGOne.bat",fileNamePrefix);
#endif
 	fileSH = fopen(fileNameSH, "w");
	if(!fileSH) return FALSE;
#ifdef OS_WIN32
	fprintf(fileSH,"set PATH=%cPATH%c;\"%s\"\n",c,c,fireflyDirectory);
#endif

	getMultiplicityName(geometry->molecule.spinMultiplicity, multiplicityStr);

	fileNameIn = strdup_printf("%sOne.inp",fileNamePrefix);
 	file = fopen(fileNameIn, "w");
	if(!file) 
	{
 		if(fileNameIn) free(fileNameIn);
 		if(fileNameOut) free(fileNameOut);
 		if(fileNameSH) free(fileNameSH);
		return FALSE;
	}
	fprintf(file,"! ======================================================\n");
	fprintf(file,"!  Input file for FireFly\n"); 
	fprintf(file,"! ======================================================\n");
	if(strstr(keyWords,"RUNTYP"))
	{
		sscanf(strstr(keyWords,"RUNTYP"),"%s",buffer);
		fprintf(file," $CONTRL %s $END\n",buffer);
	}
	if(strstr(keyWords,"SCFTYP"))
	{
		sscanf(strstr(keyWords,"SCFTYP"),"%s",buffer);
		fprintf(file," $CONTRL %s $END\n",buffer);
	}
	else
	{
		if(geometry->molecule.spinMultiplicity==1)
			fprintf(file," $CONTRL SCFTYP=RHF $END\n");
		else
			fprintf(file," $CONTRL SCFTYP=UHF $END\n");
	}

	fprintf(file," $CONTRL ICHARG=%d MULT=%d $END\n", geometry->molecule.totalCharge, geometry->molecule.spinMultiplicity);
	if(strstr(keyWords,"GBASIS"))
	{
		sscanf(strstr(keyWords,"GBASIS"),"%s",buffer);
		fprintf(file," $BASIS %s $END\n",buffer);
	}
	fprintf(file," $DATA\n");
	fprintf(file,"Molecule specification\n");
	fprintf(file,"C1\n");
	for(j=0;j<geometry->molecule.nAtoms;j++)
	{
		char* symbol = geometry->molecule.atoms[j].prop.symbol;
		SAtomsProp prop = propAtomGet(symbol);
		fprintf(file,"%s %f %f %f %f\n", 
			symbol,
			(double)prop.atomicNumber,
			geometry->molecule.atoms[j].coordinates[0],
			geometry->molecule.atoms[j].coordinates[1],
			geometry->molecule.atoms[j].coordinates[2]
			);
	}
	fprintf(file," $END\n");
	fclose(file);
	fileNameOut = strdup_printf("%sOne.out",fileNamePrefix);
#ifndef OS_WIN32
	if(!strcmp(fireflyCommand,"pcgamess") || !strcmp(fireflyCommand,"nohup pcgamess")||
	!strcmp(fireflyCommand,"firefly") || !strcmp(fireflyCommand,"nohup firefly"))
	{
		fprintf(fileSH,"mkdir %stmp\n",fileNamePrefix);
		fprintf(fileSH,"cd %stmp\n",fileNamePrefix);
		fprintf(fileSH,"cp %s input\n",fileNameIn);
		fprintf(fileSH,"%s -p -o %s\n",fireflyCommand,fileNameOut);
		fprintf(fileSH,"cd ..\n");
		fprintf(fileSH,"rm PUNCH\n");
		fprintf(fileSH,"/bin/rm -r  %stmp\n",fileNamePrefix);
	}
	else
		fprintf(fileSH,"%s %s",fireflyCommand,fileNameIn);
#else
	if(!strcmp(fireflyCommand,"pcgamess") ||
	!strcmp(fireflyCommand,"firefly") )
	{
        	fprintf(fileSH,"mkdir \"%stmp\"\n",fileNamePrefix);
		addUnitDisk(fileSH, fileNamePrefix);
	 	fprintf(fileSH,"cd \"%stmp\"\n",fileNamePrefix);
         	fprintf(fileSH,"copy \"%s\" input\n",fileNameIn);
         	fprintf(fileSH,"%s -p -o \"%s\"\n",fireflyCommand,fileNameOut);
	 	fprintf(fileSH,"cd ..\n");
         	fprintf(fileSH,"del PUNCH\n");
         	fprintf(fileSH,"del /Q  \"%stmp\"\n",fileNamePrefix);
         	fprintf(fileSH,"rmdir  \"%stmp\"\n",fileNamePrefix);
	}
	else
		fprintf(fileSH,"%s %s",fireflyCommand,fileNameIn);
#endif
	fclose(fileSH);
#ifndef OS_WIN32
	sprintf(buffer,"chmod u+x %s",fileNameSH);
	system(buffer);
	system(fileNameSH);
#else
	sprintf(buffer,"\"%s\"",fileNameSH);
	system(buffer);
#endif
	if(getEnergyFireFly(fileNameOut,energy))
	{
		char* str = NULL;

		geometry->molecule.klass->readGeomFromGamessOutputFile(&geometry->molecule, fileNameOut, -1);
		str = strdup_printf("Energy by FireFly = %f", *energy);
		printf("%s\n",str);
		if(str) free(str);
	}
	else
	{
 		if(fileNameIn) free(fileNameIn);
 		if(fileNameOut) free(fileNameOut);
 		if(fileNameSH) free(fileNameSH);
		return FALSE;
	}

 	if(fileNameIn) free(fileNameIn);
 	if(fileNameOut) free(fileNameOut);
 	if(fileNameSH) free(fileNameSH);
	return TRUE;
}
/*****************************************************************************/
static boolean runFireFlyFiles(int numberOfGeometries, ForceField** geometries, double* energies, char* fileNamePrefix, char* keyWords, char* fireflyCommand)
{
	int i;
	int nG = 0;
	int nM = 0;
	char* str = NULL;
	for(i=0;i<numberOfGeometries;i++)
	{
		if(!geometries[i]) continue;
		nG++;
		if(str) free(str);
		str = strdup_printf("Minimization by FireFly of geometry n = %d... Please wait\n", i+1);
		printf("%s",str);
		runOneFireFly(geometries[i], &energies[i], fileNamePrefix, keyWords, fireflyCommand);
		nM++;
	}
	if(str) free(str);
	if(nM==nG) return TRUE;
	return FALSE;

}
/*************************************************************************************************************/
static boolean getEnergyGaussian(char* fileNameOut, double* energy)
{
	FILE* file = NULL;
	char buffer[1024];
	char* pdest = NULL;
	boolean OK = FALSE;

 	file = fopen(fileNameOut, "r");
	if(!file) return FALSE;
	 while(!feof(file))
	 {
		if(!fgets(buffer,BSIZE,file))break;
		pdest = strstr( buffer, "SCF Done:  E(");
		if(!pdest) if(strstr( buffer, "Energy=") && !strstr( buffer, "hange") ) pdest = strstr( buffer, "Energy=");
		if(pdest) 
		{
			pdest = strstr( buffer, "=");
			if(pdest)
			{
				int l = strlen(pdest);
				int i;
				for(i=0;i<l;i++) if(pdest[i]=='D' || pdest[i]=='E') pdest[i] ='E';
				if(sscanf(pdest+1,"%lf",energy)==1)
				{
					OK = TRUE;
					/* break;*/
				}
			}
		}
	 }
	fclose(file);
	return OK;
}
/*****************************************************************************/
static boolean runOneGaussian(ForceField* geom, double* energy, char* fileNamePrefix, char* keyWords, char* gaussianCommand)
{
	FILE* file = NULL;
	FILE* fileSH = NULL;
	int j;
	char* fileNameIn = NULL;
	char* fileNameOut = NULL;
	char* fileNameSH = NULL;
	char buffer[1024];
	Molecule* mol = &geom->molecule;
	*energy = 0;
#ifdef OS_WIN32
	char c='%';
#endif

	if(!geom) return FALSE;
	if(geom->molecule.nAtoms<1) return FALSE;
#ifndef OS_WIN32
	fileNameSH = strdup_printf("%sGaussOne.sh",fileNamePrefix);
#else
	fileNameSH = strdup_printf("%sGaussOne.bat",fileNamePrefix);
#endif
 	fileSH = fopen(fileNameSH, "w");
	if(!fileSH) return FALSE;
#ifdef OS_WIN32
	fprintf(fileSH,"set PATH=%cPATH%c;\"%s\"\n",c,c,gaussianDirectory);
#endif

	fileNameIn = strdup_printf("%sOne.inp",fileNamePrefix);
 	file = fopen(fileNameIn, "w");
	if(!file) 
	{
 		if(fileNameIn) free(fileNameIn);
 		if(fileNameOut) free(fileNameOut);
 		if(fileNameSH) free(fileNameSH);
		return FALSE;
	}
 	file = fopen(fileNameIn, "w");
	if(!file) 
	{
 		if(fileNameIn) free(fileNameIn);
 		if(fileNameOut) free(fileNameOut);
 		if(fileNameSH) free(fileNameSH);
		return FALSE;
	}
	fprintf(file,"# %s\n",keyWords);
	fprintf(file,"\n");
	fprintf(file,"! ======================================================\n");
	fprintf(file,"!  Input file for Gaussian\n"); 
	fprintf(file,"! ======================================================\n");
	fprintf(file,"\n");

	fprintf(file,"%d %d\n",mol->totalCharge,mol->spinMultiplicity);
	for(j=0;j<mol->nAtoms;j++)
	{
		char* symbol = mol->atoms[j].prop.symbol;
		fprintf(file,"%s %f %f %f\n", 
			symbol,
			mol->atoms[j].coordinates[0],
			mol->atoms[j].coordinates[1],
			mol->atoms[j].coordinates[2]
			);
	}
	fprintf(file,"\n");
	fclose(file);
	fileNameOut = strdup_printf("%sOne.out",fileNamePrefix);
#ifndef OS_WIN32
	fprintf(fileSH,"%s %s",gaussianCommand,fileNameIn);
#else
	fprintf(fileSH,"%s %s",gaussianCommand,fileNameIn);
#endif
	fclose(fileSH);
	{
		char* str = NULL;
		if(strstr(keyWords,"OPT")) str = strdup_printf("Minimization by Gaussian ... Please wait");
		else str = strdup_printf("Computing of energy by Gaussian .... Please wait");
		printf("%s",str);
	}
#ifndef OS_WIN32
	sprintf(buffer,"chmod u+x %s",fileNameSH);
	system(buffer);
	system(fileNameSH);
#else
	sprintf(buffer,"\"%s\"",fileNameSH);
	system(buffer);
#endif
	if(getEnergyGaussian(fileNameOut,energy))
	{
		printf("Energy by Gaussian = %f\n", *energy);
		mol->klass->readGeomFromGaussianOutputFile(mol, fileNameOut, -1);
	}
	else
	{
 		if(fileNameIn) free(fileNameIn);
 		if(fileNameOut) free(fileNameOut);
 		if(fileNameSH) free(fileNameSH);
		return FALSE;
	}

 	if(fileNameIn) free(fileNameIn);
 	if(fileNameOut) free(fileNameOut);
 	if(fileNameSH) free(fileNameSH);
	return TRUE;
}
/*****************************************************************************/
static boolean runGaussianFiles(int numberOfGeometries, ForceField** geometries, double* energies, char* fileNamePrefix, char* keyWords, char* gaussianCommand)
{
	int i;
	int nG = 0;
	int nM = 0;
	char* str = NULL;
	for(i=0;i<numberOfGeometries;i++)
	{
		if(!geometries[i]) continue;
		nG++;
		if(str) free(str);
		printf("Minimization by Gaussian of geometry n = %d... Please wait\n", i+1);
		if(runOneGaussian(geometries[i], &energies[i], fileNamePrefix, keyWords,gaussianCommand)) 
		{
			nM++;
		}
	}
	if(str) free(str);
	if(nM==nG) return TRUE;
	return FALSE;

}
/*****************************************************************************/
static boolean getEnergyOrca(char* fileNameOut, double* energy)
{
	FILE* file = NULL;
	char buffer[1024];
	char* pdest = NULL;
	char* energyTag = "FINAL SINGLE POINT ENERGY";

 	file = fopen(fileNameOut, "r");
	if(!file) return FALSE;
	 while(!feof(file))
	 {
		if(!fgets(buffer,BSIZE,file))break;
		pdest = strstr( buffer, energyTag);
		if(pdest)
		{
			int l = strlen(pdest);
			int i;
			for(i=0;i<l;i++) if(pdest[i]=='D' || pdest[i]=='E') pdest[i] ='E';
		}
		if(pdest &&sscanf(pdest+strlen(energyTag)+1,"%lf",energy)==1)
		{
			fclose(file);
			*energy *=AUTOKCAL;
			return TRUE;
		}
	 }
	fclose(file);
	return FALSE;
}
/*****************************************************************************/
static boolean runOneOrca(ForceField* geom, double* energy, char* fileNamePrefix, char* keyWords, char* orcaCommand)
{
	FILE* file = NULL;
	FILE* fileSH = NULL;
	char* fileNameIn = NULL;
	char* fileNameOut = NULL;
	char* fileNameSH = NULL;
	char multiplicityStr[100];
	char buffer[1024];
	int i;
	int nV;
	Molecule* mol = &geom->molecule;

	if(!mol) return FALSE;
        if(mol->nAtoms<2)return FALSE;
#ifndef OS_WIN32
	fileNameSH = strdup_printf("%sOne.sh",fileNamePrefix);
#else
	fileNameSH = strdup_printf("%sOne.bat",fileNamePrefix);
#endif
 	fileSH = fopen(fileNameSH, "w");
	if(!fileSH) return FALSE;
#ifdef OS_WIN32
	fprintf(fileSH,"@echo off\n");
#endif

	getMultiplicityName(mol->spinMultiplicity, multiplicityStr);

	fileNameIn = strdup_printf("%sOne.inp",fileNamePrefix);
 	file = fopen(fileNameIn, "w");
	if(!file) 
	{
 		if(fileNameIn) free(fileNameIn);
 		if(fileNameOut) free(fileNameOut);
 		if(fileNameSH) free(fileNameSH);
		return FALSE;
	}
	fprintf(file,"# ======================================================\n");
	fprintf(file,"#  Orca input file made in Gabedit\n"); 
	fprintf(file,"# ======================================================\n");
	fprintf(file,"! %s\n",keyWords);
	{
		fprintf(file,"* xyz %d   %d\n",mol->totalCharge,mol->spinMultiplicity);
		nV = 0;
      		for (i=0;i<mol->nAtoms;i++)
		{
			char X[100];
			char Y[100];
			char Z[100];
			sprintf(X,"%20.14f",mol->atoms[i].coordinates[0]);
			sprintf(Y,"%20.14f",mol->atoms[i].coordinates[1]);
			sprintf(Z,"%20.14f",mol->atoms[i].coordinates[2]);

			fprintf(file," %s  %s %s %s\n",mol->atoms[i].prop.symbol, X,Y,Z);
			if(mol->atoms[i].variable) nV+=3;
		}
		fprintf(file,"*\n");
		if(nV>0&&nV!=3*mol->nAtoms) 
		{
			fprintf(file,"%cgeom Constraints\n",'%');
      			for (i=0;i<mol->nAtoms;i++)
			{
				if(mol->atoms[i].variable)
				{
					fprintf(file,"  {C %d C}\n",i);
				}
			}
			fprintf(file," end #Constraints\n");
			fprintf(file," invertConstraints true\n");
			fprintf(file," end #geom\n");
		}
	}

	fclose(file);
	fileNameOut = strdup_printf("%sOne.out",fileNamePrefix);
#ifndef OS_WIN32
	if(!strcmp(orcaCommand,"orca") || !strcmp(orcaCommand,"nohup orca"))
	{
		fprintf(fileSH,"%s %s > %s\n",orcaCommand,fileNameIn,fileNameOut);
		fprintf(fileSH,"exit\n");
	}
	else
		fprintf(fileSH,"%s %s",orcaCommand,fileNameIn);
#else
	 if(!strcmp(orcaCommand,"orca") )
	{
		if(strstr(orcaDirectory,"\"")) fprintf(fileSH,"set PATH=%s;%cPATH%c\n",orcaDirectory,'%','%');
		else fprintf(fileSH,"set PATH=\"%s\";%cPATH%c\n",orcaDirectory,'%','%');
		fprintf(fileSH,"%s %s > %s\n",orcaCommand,fileNameIn,fileNameOut);
		fprintf(fileSH,"exit\n");
	}
	else
		fprintf(fileSH,"%s %s",orcaCommand,fileNameIn);
#endif
	fclose(fileSH);
	{
		char* str = NULL;
		if(strstr(keyWords,"Opt")) str = strdup_printf("Minimization by Orca ... Please wait");
		else str = strdup_printf("Computing of energy by Orca .... Please wait");
		printf("%s\n",str);
		if(str) free(str);
	}
#ifndef OS_WIN32
	sprintf(buffer,"chmod u+x %s",fileNameSH);
	system(buffer);
	system(fileNameSH);
#else
	sprintf(buffer,"\"%s\"",fileNameSH);
	system(buffer);
#endif
	if(getEnergyOrca(fileNameOut,energy))
	{
		printf("Energy by Orca = %f\n", *energy);
		mol->klass->readGeomFromOrcaOutputFile(mol, fileNameOut, -1);
	}
	else
	{
 		if(fileNameIn) free(fileNameIn);
 		if(fileNameOut) free(fileNameOut);
 		if(fileNameSH) free(fileNameSH);
		return FALSE;
	}

 	if(fileNameIn) free(fileNameIn);
 	if(fileNameOut) free(fileNameOut);
 	if(fileNameSH) free(fileNameSH);
	return TRUE;
}
/*****************************************************************************/
static boolean runOrcaFiles(int numberOfGeometries, ForceField** geometries, double* energies, char* fileNamePrefix, char* keyWords, char* orcaCommand)
{
	int i;
	int nG = 0;
	int nM = 0;
	char* str = NULL;
	for(i=0;i<numberOfGeometries;i++)
	{
		if(!geometries[i]) continue;
		nG++;
		if(str) free(str);
		printf("Minimization by Orca of geometry n = %d... Please wait\n", i+1);
		if(runOneOrca(geometries[i], &energies[i], fileNamePrefix, keyWords,orcaCommand)) 
		{
			nM++;
		}
	}
	if(str) free(str);
	if(nM==nG) return TRUE;
	return FALSE;

}
/*****************************************************************************/
static boolean testEqualDistances(double* distancesI, double* distancesJ, int n, double tol)
{
	int k;
	if(!distancesI) return FALSE;
	if(!distancesJ) return FALSE;
	if(n<1) return FALSE;
	for (  k = 0; k < n; k++ )
		if(fabs(distancesI[k]-distancesJ[k])>tol) return FALSE;
	return TRUE;
}
/*****************************************************************************/
static double* getDistancesBetweenAtoms(ForceField* forceField)
{
	double* distances = NULL;
	int i;
	int j;
	int n;
	int k;
	if(forceField->molecule.nAtoms<1) return distances;
	n = forceField->molecule.nAtoms*(forceField->molecule.nAtoms-1)/2;
	distances = malloc(n*sizeof(double));
	n = 0;
	for (  i = 0; i < forceField->molecule.nAtoms-1; i++ )
	for (  j = i+1; j < forceField->molecule.nAtoms; j++ )
	{
		double x = forceField->molecule.atoms[i].coordinates[0]-forceField->molecule.atoms[j].coordinates[0];
		double y = forceField->molecule.atoms[i].coordinates[1]-forceField->molecule.atoms[j].coordinates[1];
		double z = forceField->molecule.atoms[i].coordinates[2]-forceField->molecule.atoms[j].coordinates[2];
		distances[n++] = x*x + y*y + z*z;
	}
	for(i=0;i<n-1;i++)
	{
		k = i;
		for(j=i+1;j<n;j++)
			if(distances[j]<distances[k]) k= j;
		if(k!=i)
		{
			double d = distances[i];
			distances[i] = distances[k];
			distances[k] = d;
		}
	}
	return distances;
}
/*****************************************************************************/
static void removedsToEnd(int numberOfGeometries, ForceField** geometries, double* energies, boolean* removeds)
{
	if(geometries && energies && removeds)
	{
		int i;
		int j;
		int k;
		for(i=0;i<numberOfGeometries-1;i++)
		{
			if(!removeds[i]) continue;
			k = i;
			for(j=i+1;j<numberOfGeometries;j++)
				if(!removeds[j]) { k= j; break;}
			if(k!=i)
			{
				double energy = energies[i];
				boolean r = removeds[i];
				ForceField* g = geometries[i];

				energies[i] = energies[k];
				energies[k] = energy;
				geometries[i] = geometries[k];
				geometries[k] = g;
				removeds[i] = removeds[k];
				removeds[k] = r;
			}
		}
	}
}
/*****************************************************************************/
static void computeRemoveds(int numberOfGeometries, ForceField** geometries, double* energies, boolean *removeds, 
		double tolEnergy, double tolDistance)
{
	int i;
	int j;
	double* distancesI = NULL;
	double* distancesJ = NULL;
	if(tolDistance<=0 && tolEnergy<=0) return;
	if(!geometries || !energies) return;
	if(numberOfGeometries<1) return;
	i = numberOfGeometries-1;
	if(!geometries[i]){ removeds[i]=TRUE;}
	for(i=0;i<numberOfGeometries-1;i++)
	{
		int n;
		if(removeds[i]) continue;
		if(!geometries[i]){ removeds[i]=TRUE; continue;}
		if(tolDistance>0) distancesI =  getDistancesBetweenAtoms(geometries[i]);
		n = geometries[i]->molecule.nAtoms*(geometries[i]->molecule.nAtoms-1)/2;
		for(j=i+1;j<numberOfGeometries;j++)
		{
			if(removeds[j]) continue;
			if(!geometries[j]){ removeds[j]=TRUE; continue;}
			if(tolEnergy>0 && fabs(energies[j]-energies[i])<tolEnergy && geometries[i]->molecule.nAtoms==geometries[j]->molecule.nAtoms)
			{
				if(tolDistance>0) 
				{
					distancesJ =  getDistancesBetweenAtoms(geometries[j]);
					if(testEqualDistances(distancesI, distancesJ, n, tolDistance))
						removeds[j] = TRUE;
					if(distancesJ) free(distancesJ);
					distancesJ = NULL;
				}
				else
					removeds[j] = TRUE;
			}
			if(tolEnergy<0 && tolDistance>0 && geometries[i]->molecule.nAtoms==geometries[j]->molecule.nAtoms)
			{
				distancesJ =  getDistancesBetweenAtoms(geometries[j]);
				if(testEqualDistances(distancesI, distancesJ, n, tolDistance))
					removeds[j] = TRUE;
				if(distancesJ) free(distancesJ);
				distancesJ = NULL;
			}
		}
		if(distancesI) free(distancesI);
		distancesI = NULL;
	}

}
/*****************************************************************************/
static void removeIdenticalGeometries(int* nG, ForceField*** geoms, double** eners, double tolEnergy, double tolDistance)
{
	int i;
	int numberOfGeometries =*nG;
	ForceField** geometries = *geoms; 
	double* energies = *eners;
	boolean* removeds = NULL;
	int newN = 0;
	if(numberOfGeometries<1) return;
	removeds = malloc(numberOfGeometries*sizeof(boolean));
	for(i=0;i<numberOfGeometries;i++) removeds[i] = FALSE;
	computeRemoveds(numberOfGeometries, geometries, energies, removeds, tolEnergy, tolDistance);
	removedsToEnd(numberOfGeometries, geometries, energies, removeds);

	for(i=0;i<numberOfGeometries;i++) 
	{
		if(removeds[i]) 
		{
			if(geometries[i]) freeForceField(geometries[i]);
		}
		else newN++;
	}
	free(removeds);
	if(newN==0) { *nG = newN; return;}
	if(newN==numberOfGeometries) return;
	*nG = newN;
	*eners = realloc(*eners,newN*sizeof(double));
	*geoms = realloc(*geoms,newN*sizeof(ForceField**));

}
/*****************************************************************************/
/*
static int removeIdenticalGeometriesNULL(int numberOfGeometries, ForceField** geometries, double* energies, double tolEnergy, double tolDistance)
{
	int i;
	boolean* removeds = NULL;
	int newN = 0;
	if(numberOfGeometries<1) return 0;
	removeds = malloc(numberOfGeometries*sizeof(boolean));
	for(i=0;i<numberOfGeometries;i++) removeds[i] = FALSE;
	computeRemoveds(numberOfGeometries, geometries, energies, removeds, tolEnergy, tolDistance);
	removedsToEnd(numberOfGeometries, geometries, energies, removeds);

	for(i=0;i<numberOfGeometries;i++) 
	{
		if(removeds[i]) 
		{
			if(geometries[i]) freeForceField(geometries[i]);
			geometries[i] = NULL;
			energies[i] = 1e30;
		}
		else newN++;
	}
	free(removeds);
	return newN;
}
*/
/*****************************************************************************/
static void sortGeometries(int numberOfGeometries, ForceField** geometries, double* energies)
{
	if(geometries && energies)
	{
		int i;
		int j;
		int k;
		for(i=0;i<numberOfGeometries-1;i++)
		{
			k = i;
			for(j=i+1;j<numberOfGeometries;j++)
				if(energies[j]<energies[k]) k= j;
			if(k!=i)
			{
				double energy = energies[i];
				ForceField* g = geometries[i];

				energies[i] = energies[k];
				energies[k] = energy;
				geometries[i] = geometries[k];
				geometries[k] = g;
			}
		}
	}
}
/*****************************************************************************/
static void createPostProcessingFiles(int numberOfGeometries, ForceField** geometries,double* energies,char* fileNameGeom, char* mopacKeywords, char* gaussianKeywords, char* fireflyKeywords, char* orcaKeywords, char* cchemiKeywords, char* message, char* mopacCommand, char* gaussianCommand, char* fireflyCommand, char* orcaCommand, char* cchemiCommand)
{
	if(mopacKeywords)
	{
		char* fileNamePrefix = getSuffixNameFile(fileNameGeom);
		createMopacFiles(numberOfGeometries, geometries, energies, fileNamePrefix, mopacKeywords,mopacCommand);
		strcat(message,fileNamePrefix);
		strcat(message,("_*.mop\n\tFiles for a post processing by Mopac\n\n"));
		if(fileNamePrefix) free(fileNamePrefix);
	}
	if(gaussianKeywords)
	{
		char* fileNamePrefix = getSuffixNameFile(fileNameGeom);
		createGaussianFiles(numberOfGeometries, geometries, energies, fileNamePrefix, gaussianKeywords,gaussianCommand);
		strcat(message,fileNamePrefix);
		strcat(message,("_*.com\n\tFiles for a post processing by Gaussian\n\n"));
		if(fileNamePrefix) free(fileNamePrefix);
	}
	if(orcaKeywords)
	{
		char* fileNamePrefix = getSuffixNameFile(fileNameGeom);
		createOrcaFiles(numberOfGeometries, geometries, energies, fileNamePrefix, orcaKeywords,orcaCommand);
		strcat(message,fileNamePrefix);
		strcat(message,("ORCA_*.inp\n\tFiles for a post processing by Orca\n\n"));
		if(fileNamePrefix) free(fileNamePrefix);
	}
	if(fireflyKeywords)
	{
		char* fileNamePrefix = getSuffixNameFile(fileNameGeom);
		createFireFlyFiles(numberOfGeometries, geometries, energies, fileNamePrefix, fireflyKeywords,fireflyCommand);
		strcat(message,fileNamePrefix);
		strcat(message,("FF_*.inp\n\tFiles for a post processing by FireFly\n\n"));
		if(fileNamePrefix) free(fileNamePrefix);
	}
	if(cchemiKeywords)
	{
		char* fileNamePrefix = getSuffixNameFile(fileNameGeom);
		createCChemIFiles(numberOfGeometries, geometries, energies, fileNamePrefix, cchemiKeywords,cchemiCommand);
		strcat(message,fileNamePrefix);
		strcat(message,("CCHEMI_*.inp\n\tFiles for a post processing by CChemI\n\n"));
		if(fileNamePrefix) free(fileNamePrefix);
	}
}
/*****************************************************************************/
char*  setOptOptions(FILE* file, ConjugateGradientOptions* conjugateGradientOptions, QuasiNewton* quasiNewton)
{
	char* optimizerType = strdup("QuasiNewton");
	readOneString(file,"OptimizerType",&optimizerType);
	if(strstr(optimizerType,"Gradient")) setCGOptions(file, conjugateGradientOptions);
	else if(strstr(optimizerType,"Quasi")) setQNOptions(file, quasiNewton);
	else if(strstr(optimizerType,"Steep")) setCGOptions(file, conjugateGradientOptions);
	else if(strstr(optimizerType,"External"))
	{
		printf("Sorry, ExternalOptimizer with a MM potential is not implemented\n");
		exit(1);
	}
	return optimizerType;
}
/*****************************************************************************/
static int collectGeometriesFromProcessors(int nproc, ForceField** geometriesAll, int numberOfGeometriesMax, double* energiesAll, ForceField** geometries, int numberOfGeometries, double* energies,  double* coords,  double* enerDum, ForceField* forceField, double tolEnergy, double tolDistance)
{
		int numberOfGeometriesAll = 0;
#ifdef ENABLE_MPI
#ifdef DEBUG
		printf("Begin collectGeometriesFromProcessors\n");
#endif
		int j;
		int i;
		int code,tag;
		for(i=0;i<numberOfGeometriesMax;i++) 
		{
			energiesAll[i] = 1e30;
			if(geometriesAll[i]) freeForceField(geometriesAll[i]);
			geometriesAll[i] = NULL;
		}
		for(i=0;i<numberOfGeometries;i++) 
		{
			if(geometries[i])
			{
				geometriesAll[i] = malloc(sizeof(ForceField));
               			*geometriesAll[i] = copyForceField(geometries[i]);
				energiesAll[i] = energies[i];
			}
		}
		numberOfGeometriesAll = numberOfGeometries;
		// get geometries from other proc
		for(j=1;j<nproc;j++)
		{
			int nG = 0;
			tag = 1000;
			MPI_Status status ;
#ifdef DEBUG
			printf("get nGeometries from proc n %d\n", j);
#endif
			code = MPI_Recv(&nG,1,MPI_INT,j,tag,MPI_COMM_WORLD,&status) ;
#ifdef DEBUG
			printf("nGeometries=%d from proc n %d\n",nG, j);
#endif
			if(nG>0) 
			{
				int k;
				int a,b;
				tag = 2000;
				code = MPI_Recv(enerDum,nG,MPI_DOUBLE,j,tag,MPI_COMM_WORLD,&status) ;
				for(k=0;k<nG;k++)
				{
					int nA = forceField->molecule.nAtoms;
					energiesAll[i+k] = enerDum[k];
					tag = 3000+k;
					code = MPI_Recv(coords,nA*3,MPI_DOUBLE,j,tag,MPI_COMM_WORLD,&status) ;
					geometriesAll[i+k] = malloc(sizeof(ForceField));
               				*geometriesAll[i+k] = copyForceField(forceField);
					b = 0;
					for(a=0;a<nA;a++)
					{
#ifdef DEBUG
						printf(" atoms %d C = %f %f %f\n",a,coords[b], coords[b+1], coords[b+2]);
#endif
						geometriesAll[i+k]->molecule.atoms[a].coordinates[0] = coords[b++];
						geometriesAll[i+k]->molecule.atoms[a].coordinates[1] = coords[b++];
						geometriesAll[i+k]->molecule.atoms[a].coordinates[2] = coords[b++];
					}
				}
				i += nG;
				numberOfGeometriesAll += nG;
			}
		}
		sortGeometries(numberOfGeometriesAll, geometriesAll, energiesAll);
		numberOfGeometriesAll = removeIdenticalGeometriesNULL(numberOfGeometriesAll, geometriesAll, energiesAll, tolEnergy, tolDistance);
#ifdef DEBUG
		printf("End collectGeometriesFromProcessors\n");
#endif

#endif /* ENABLE_MPI */
		return numberOfGeometriesAll;
}
/*********************************************************************************************************************************/
static void sendGeometriesToMaster(int rank, ForceField** geometries, int numberOfGeometries, double* energies,  double* coords, int nAtoms)
{
#ifdef ENABLE_MPI
	int j;
	int code,tag;
	j = 0;
	tag = 1000;
#ifdef DEBUG
	printf("Begin sendGeometriesToMaster from proc n %d\n",rank);
#endif
	code = MPI_Send(&numberOfGeometries,1,MPI_INT,j,tag,MPI_COMM_WORLD) ;
	if(numberOfGeometries>0)
	{
		int k;
		int a,b;
		tag = 2000;
		code = MPI_Send(energies,numberOfGeometries,MPI_DOUBLE,j,tag,MPI_COMM_WORLD) ;
		for(k=0;k<numberOfGeometries;k++)
		{
			tag = 3000+k;
			b = 0;
			for(a=0;a<nAtoms;a++)
			{
				coords[b++] = geometries[k]->molecule.atoms[a].coordinates[0];
				coords[b++] = geometries[k]->molecule.atoms[a].coordinates[1];
				coords[b++] = geometries[k]->molecule.atoms[a].coordinates[2];
			}
			code = MPI_Send(coords,nAtoms*3,MPI_DOUBLE,j,tag,MPI_COMM_WORLD);
		}
	}
#ifdef DEBUG
	printf("End sendGeometriesToMaster from proc n %d\n",rank);
#endif
#endif /* ENABLE_MPI */
}
/*****************************************************************************/
void molecularMechanicsDynamicsREMDConfoDlg(char* inputFileName)
{
	ForceField forceField; 
	ForceFieldOptions forceFieldOptions;
	MolecularDynamics molecularDynamics;
	int updateFrequency = 1;
	double heatTime;
	double equiTime;
	double runTime;
	double coolTime;
	double heatTemp; 
	double equiTemp; 
	double runTemp; 
	double coolTemp; 
	double stepSize;
	MDIntegratorType integrator = VERLET;
	char* fileNameGeom = NULL;
	char* fileNameTraj = NULL;
	char* fileNameProp = NULL;
	char* mopacKeywordsPost = NULL;
	char* gaussianKeywordsPost = NULL;
	char* fireflyKeywordsPost = NULL;
	char* orcaKeywordsPost = NULL;
	char* cchemiKeywordsPost = NULL;
	double friction=-1;
	double omegaMax = 4000;
	int Nf = 50;
	double collide = 20;
	double qNH = 20;
	MDThermostatType thermostat = NONE;
	int numberOfGeometries = 2;
	ForceField** geometries = NULL; 
	double* energies = NULL;
	char* optMopacMethod=strdup("PM6");
	char* optGaussianMethod=strdup("AM1");
	char* optFireFlyMethod=strdup("AM1");
	char* optOrcaMethod=strdup("AM1");
	boolean optMopac = FALSE;
	boolean optFireFly = FALSE;
	boolean optGaussian = FALSE;
	boolean optOrca = FALSE;
	boolean optMM = FALSE;
	Molecule mol = *(readMolecule(inputFileName,TRUE));
	char* optimizerType= strdup("QuasiNewton");

	QuasiNewton quasiNewton;
	ConjugateGradientOptions conjugateGradientOptions;
	SteepestDescent steepestDescent;
	ConjugateGradient conjugateGradient;
	int i;
	char message[BSIZE]=" ";
	double tolEnergy = -1;
	double tolDistance = -1;
	char* mopacCommand = strdup("mopac");
	char* gaussianCommand=strdup("g03"); 
	char* orcaCommand=strdup("orca"); 
	char* fireflyCommand=strdup("firefly");
	char* cchemiCommand=strdup("cchemi");
	FILE* file = fopen(inputFileName,"rb");
	double runTempMax = 700;
	int nTemperatures = 10;
	int numberOfExchanges = 10;
	double timeExchange = 1;
	int nproc;
	int rank;
	ForceField** geometriesAll = NULL; 
	int numberOfGeometriesAll = 2;
	int numberOfGeometriesMax = 2;
	double* energiesAll = NULL;
	double* coords = NULL;
	double* enerDum = NULL;
	FILE* logfile = NULL;
#ifdef ENABLE_MPI
	MPI_Comm_rank( MPI_COMM_WORLD,&rank);
	MPI_Comm_size( MPI_COMM_WORLD,&nproc );
#else
	rank = 0;
	nproc = 1;
#endif

	logfile = stdout;
	if(rank!=0)
	{
		char* fileNamePrefix = getSuffixNameFile(inputFileName);
		char* tmp = strdup_printf("%s_%d.log",fileNamePrefix, rank);
		logfile = fopen(tmp,"w");
		free(tmp);
		free(fileNamePrefix);
	}
	fprintf(logfile, "MolecularMechanicsDynamicsREMDConfoDlg Rank#=%d  nproc = %d\n", rank, nproc );
#ifdef DEBUG
	fprintf(logfile, "MolecularMechanicsDynamicsREMDConfoDlg Rank#=%d  nproc = %d\n", rank, nproc );
#endif
	setForceFieldOptions(file, &forceFieldOptions);
	mol.klass->buildMMTypes(&mol, file);

	setMDOptions(file, &updateFrequency, 
		&heatTime, &equiTime, &runTime, &coolTime,
		&heatTemp, &runTemp, &equiTemp, &coolTemp, &stepSize, 
		&integrator, &thermostat, &friction, &omegaMax, &Nf, &collide,&qNH);
	if(thermostat == NONE) 
	{
		fprintf(logfile,"Warning....................\n");
		fprintf(logfile," A thermostat is required for a REMD calculation\n");
		fprintf(logfile," I set it to Berendsen\n");
		thermostat = BERENDSEN;
	}
	readOneReal(file,"runTempMax",&runTempMax);
	if(runTempMax<=runTemp) runTempMax=10*runTemp;
	readOneInt(file,"nTemperatures",&nTemperatures);
	if(nTemperatures<1) nTemperatures=10;

	if(readOneReal(file,"timeExchange",&timeExchange))
	{
		if(timeExchange>=runTime)
		{
			fprintf(logfile, "Error : time of exchange cannot be larger than run time \n");
			fprintf(logfile, "      : check your input file\n");
			exit(1);
		}
	}
	else timeExchange = runTime/10;

	numberOfExchanges = (int)(runTime/timeExchange+0.5);
	fprintf(logfile, "Number of exchanges = %d\n",numberOfExchanges);

	if(numberOfExchanges<1) numberOfExchanges=2;
	readOneString(file,"mopacCommand",&mopacCommand);
	readOneString(file,"gaussianCommand",&gaussianCommand);
	readOneString(file,"fireflyCommand",&fireflyCommand);
	readOneString(file,"orcaCommand",&orcaCommand);
	readOneString(file,"cchemiCommand",&cchemiCommand);
	readOneReal(file,"tolEnergy",&tolEnergy);
	readOneReal(file,"tolDistance",&tolDistance);
	readOneBoolean(file,"ConfoOptMM",&optMM);
	readOneBoolean(file,"ConfoOptMopac",&optMopac);
	readOneString(file,"ConfoOptMopacMethod",&optMopacMethod);
	readOneBoolean(file,"ConfoOptGaussian",&optGaussian);
	readOneString(file,"ConfoOptGaussianMethod",&optGaussianMethod);
	readOneBoolean(file,"ConfoOptFireFly",&optFireFly);
	readOneString(file,"ConfoOptFireFlyMethod",&optFireFlyMethod);
	readOneBoolean(file,"ConfoOptOrac",&optOrca);
	readOneString(file,"ConfoOptOrcaMethod",&optOrcaMethod);

	readOneString(file,"mopacKeywordsPost",&mopacKeywordsPost);
	readOneString(file,"gaussianKeywordsPost",&gaussianKeywordsPost);
	readOneString(file,"fireflyKeywordsPost",&fireflyKeywordsPost);
	readOneString(file,"orcaKeywordsPost",&orcaKeywordsPost);
	readOneString(file,"cchemiKeywordsPost",&cchemiKeywordsPost);
	optimizerType = setOptOptions(file, &conjugateGradientOptions, &quasiNewton);
	if(strstr(optimizerType,"External")  && (optMopac||optFireFly || optGaussian || optOrca))
	{
		printf("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n");
		printf("Sorry, The optimization using the internal optimizer after a MD conformational search is not yes implemented in this software\n");
		printf("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n");
		exit(1);
	}
	/* number for geometries */
	{
		numberOfGeometries = 10;
		readOneInt(file,"numberOfGeometries",&numberOfGeometries);
		if(numberOfGeometries<2) numberOfGeometries = 2;
	}
	/* fileNames  */
	{
		char* suff = getSuffixNameFile(inputFileName);
		fileNameGeom = strdup_printf("%s%s",suff, "Geoms.gab");
		fileNameTraj = strdup_printf("%s%s",suff, "Traj.gab");
		fileNameProp = strdup_printf("%s%s",suff, "Prop.txt");
		free(suff);
	}

/* Optimsation options */ 

	/* Molecule to read */
	if(forceFieldOptions.type==AMBER) forceField = createAmberModel(&mol,forceFieldOptions, logfile);
	else if(forceFieldOptions.type==PAIRWISE) forceField = createPairWiseModel(&mol,forceFieldOptions,logfile);
	setH4CorrectionMM(file, &forceField);
	checkWallCorrection(file, &forceField);

	geometries = runREMD(&molecularDynamics, &forceField,
		updateFrequency, heatTime, equiTime, runTime, heatTemp, runTemp, runTempMax, stepSize, 
		integrator, thermostat, friction, omegaMax, Nf, collide, qNH, numberOfGeometries, nTemperatures, numberOfExchanges, fileNameTraj, fileNameProp);

	if(geometries) 
	{
		int k = 0;
		int i = 0;
#ifdef DEBUG
		fprintf(logfile, "number max of Geometries in rank n %d = %d\n",rank,numberOfGeometries );
#endif
		for(i=0;i<numberOfGeometries;i++) if(geometries[i])k++;
		fprintf(logfile, "number of selected geometries on rank n %d = %d\n",rank, k);
		if(k>0) fprintf(logfile, "fileNameGeom = %s\n",fileNameGeom);
	}
	else fprintf(logfile, "No selected geometriy in rank n %d\n",rank);
	fflush(logfile);

	if(geometries && numberOfGeometries>0) energies = malloc(numberOfGeometries*sizeof(double));
	numberOfGeometriesAll = numberOfGeometries+1;
	numberOfGeometriesMax = numberOfGeometries+1;
	if(nproc>1)
	{
		coords =  malloc(forceField.molecule.nAtoms*3*sizeof(double));
		enerDum =  malloc(numberOfGeometriesMax*sizeof(double));
		if(geometries && numberOfGeometries>0) energiesAll = malloc(numberOfGeometriesMax*sizeof(double));
		if(geometries && numberOfGeometries>0) 
		{
			geometriesAll = malloc(numberOfGeometriesMax*sizeof(ForceField*));
			for(i=0;i<numberOfGeometriesMax;i++) 
			{
				geometriesAll[i] = NULL;
			/*
				geometriesAll[i] = malloc(sizeof(ForceField));
                		*geometries[i] = copyForceField(&forceField);
			*/
			}
		}
	}

	if(geometries && optMM)
	for(i=0;i<numberOfGeometries;i++)
	{
		char* str = NULL;
		energies[i] = 1e30;
		if(!geometries[i]) continue;
		if(str) free(str);
		str = strdup_printf("Minimization of geometry number %d\n", i+1);
		fprintf(logfile,"%s",str);
		fflush(logfile);
		if(str) free(str);


		if(strstr(optimizerType,"Grad"))
		{
			conjugateGradient.logfile= logfile;
			runConjugateGradient(&conjugateGradient, geometries[i], conjugateGradientOptions); 
			energies[i] = conjugateGradient.forceField->klass->calculateEnergyTmp
				(conjugateGradient.forceField, &conjugateGradient.forceField->molecule );
			freeConjugateGradient(&conjugateGradient);
		}
		else if(strstr(optimizerType,"Quasi"))
		{
			QuasiNewton tmpQuasiNewton = quasiNewton;
			tmpQuasiNewton.forceField = geometries[i];
                	tmpQuasiNewton.logfile = logfile;
                	runQuasiNewton(&tmpQuasiNewton);
			energies[i] = tmpQuasiNewton.forceField->klass->calculateEnergyTmp
				(tmpQuasiNewton.forceField, &tmpQuasiNewton.forceField->molecule );
			freeQuasiNewton(&tmpQuasiNewton);

		}
		else
		{
			steepestDescent.logfile= logfile;
			runSteepestDescent(&steepestDescent, geometries[i],
			       	conjugateGradientOptions.updateFrequency,
			       conjugateGradientOptions.maxIterations,
			       conjugateGradientOptions.gradientNorm,
			       conjugateGradientOptions.maxLines);
			energies[i] = steepestDescent.forceField->klass->calculateEnergyTmp
				(steepestDescent.forceField, &steepestDescent.forceField->molecule );
			freeSteepestDescent(&steepestDescent);
		}
		str = strdup_printf("End Minimization of geometry number %d\n", i+1);
		fprintf(logfile,"%s",str);
		fflush(logfile);
		if(str) free(str);
	}
	else 
	{
		for(i=0;i<numberOfGeometries;i++)
		{
			energies[i] = 1e30;
			if(!geometries[i]) continue;
			energies[i] = geometries[i]->klass->calculateEnergyTmp (geometries[i], &geometries[i]->molecule );
		}

	}

	/*  sort by energies */
	{
//#ifdef DEBUG
		fprintf(logfile, "begin sort geometry rank = %d\n",rank);
		fflush(logfile);
//#endif
		sortGeometries(numberOfGeometries, geometries, energies);
//#ifdef DEBUG
		fprintf(logfile, "end sort geometry rank = %d\n",rank);
		fflush(logfile);
//#endif
		removeIdenticalGeometries(&numberOfGeometries, &geometries, &energies,tolEnergy,tolDistance);
//#ifdef DEBUG
		fprintf(logfile, "end removeIdenticalGeometries geometry rank = %d, numberOfGeometries= %d\n",rank, numberOfGeometries);
		fflush(logfile);
//#endif
	}
	if(nproc==1)
	{
		numberOfGeometriesAll = numberOfGeometries;
		geometriesAll = geometries;
		energiesAll = energies;
	}
	else
	{
//#ifdef DEBUG
		fprintf(logfile, "MolecularMechanicsDynamicsREMDConfoDlg Rank#=%d  nproc = %d\n", rank, nproc );
//#endif
		if(rank==0) 
		numberOfGeometriesAll = collectGeometriesFromProcessors(nproc, geometriesAll, numberOfGeometriesMax, energiesAll, geometries, numberOfGeometries, energies,  coords,  enerDum, &forceField,tolEnergy,tolDistance);
		else sendGeometriesToMaster(rank, geometries, numberOfGeometries, energies,  coords, forceField.molecule.nAtoms);
		fprintf(logfile, "End collect&send nGeoms = %d\n", numberOfGeometriesAll );
		fflush(logfile);
	}
	if(rank==0 && saveConfoGeometries(numberOfGeometriesAll, geometriesAll, energiesAll, fileNameGeom))
	{
		sprintf(message,"Created files :\n");
		createPostProcessingFiles(numberOfGeometriesAll, geometriesAll,energiesAll,fileNameGeom, mopacKeywordsPost, gaussianKeywordsPost, fireflyKeywordsPost, orcaKeywordsPost, cchemiKeywordsPost, message, mopacCommand, gaussianCommand, fireflyCommand, orcaCommand, cchemiCommand);
		strcat(message,fileNameGeom);
		strcat(message,("\n\tGeometries selected and optimized using your MM potentials"));
		strcat(message,("\n\tTo read this file through Gabedit: 'Read/CChemI file'\n\n"));
	}
	fprintf(logfile, "End saveConfoGeometries\n" );
	fflush(logfile);
	/* minimazation by mopac*/
	if(optMopac)
	{
		char* fileNamePrefix = getSuffixNameFile(fileNameGeom);
		char* keys = strdup_printf("%s XYZ",optMopacMethod);
		if(runMopacFiles(numberOfGeometries, geometries, energies, fileNamePrefix, keys, mopacCommand) )
		{
			char* fileNameGeomMop = strdup_printf("%sMop.gab",fileNamePrefix);
			sortGeometries(numberOfGeometries, geometries, energies);
			removeIdenticalGeometries(&numberOfGeometries, &geometries, &energies,tolEnergy,tolDistance);
			if(nproc==1)
			{
				numberOfGeometriesAll = numberOfGeometries;
				geometriesAll = geometries;
				energiesAll = energies;
			}
			else
			{
				if(rank==0) 
				numberOfGeometriesAll = collectGeometriesFromProcessors(nproc, geometriesAll, numberOfGeometriesMax, energiesAll, geometries, numberOfGeometries, energies,  coords,  enerDum, &forceField, tolEnergy,tolDistance);
				else sendGeometriesToMaster(rank, geometries, numberOfGeometries, energies,  coords, forceField.molecule.nAtoms);
			}
			if(rank==0 && saveConfoGeometries(numberOfGeometriesAll, geometriesAll, energiesAll, fileNameGeomMop))
			{
				createPostProcessingFiles(numberOfGeometriesAll, geometriesAll,energiesAll,fileNameGeomMop, mopacKeywordsPost, gaussianKeywordsPost, fireflyKeywordsPost, orcaKeywordsPost, cchemiKeywordsPost, message, mopacCommand, gaussianCommand, fireflyCommand, orcaCommand, cchemiCommand);
				strcat(message,fileNameGeomMop);
				strcat(message,("\n\tGeometries after minimization by Mopac/"));
				strcat(message,optMopacMethod);
				strcat(message,("\n\tTo read this file through Gabedit : 'Read/CChemI file'\n\n"));
			}

			free(fileNameGeomMop);
		}
		if(fileNamePrefix) free(fileNamePrefix);
	}
	/* minimazation by FireFly*/
	if(optFireFly)
	{
		char* fileNamePrefix = getSuffixNameFile(fileNameGeom);
		char* keys = strdup_printf("RUNTYP=Optimize GBASIS=%s",optFireFlyMethod);
		if(runFireFlyFiles(numberOfGeometries, geometries, energies, fileNamePrefix, keys,fireflyCommand) )
		{
			char* fileNameGeomFireFly = strdup_printf("%sFireFly.gab",fileNamePrefix);
			sortGeometries(numberOfGeometries, geometries, energies);
			removeIdenticalGeometries(&numberOfGeometries, &geometries, &energies,tolEnergy,tolDistance);
			if(nproc==1)
			{
				numberOfGeometriesAll = numberOfGeometries;
				geometriesAll = geometries;
				energiesAll = energies;
			}
			else
			{
				if(rank==0) numberOfGeometriesAll = collectGeometriesFromProcessors(nproc, geometriesAll, numberOfGeometriesMax, energiesAll, geometries, numberOfGeometries, energies,  coords,  enerDum, &forceField, tolEnergy,tolDistance);
				else sendGeometriesToMaster(rank, geometries, numberOfGeometries, energies,  coords, forceField.molecule.nAtoms);
			}
			if(rank==0 && saveConfoGeometries(numberOfGeometriesAll, geometriesAll, energiesAll, fileNameGeomFireFly))
			{
				createPostProcessingFiles(numberOfGeometriesAll, geometriesAll,energiesAll, fileNameGeomFireFly, mopacKeywordsPost, gaussianKeywordsPost, fireflyKeywordsPost, orcaKeywordsPost, cchemiKeywordsPost, message, mopacCommand, gaussianCommand, fireflyCommand, orcaCommand, cchemiCommand);
				strcat(message,fileNameGeomFireFly);
				strcat(message,("\n\tGeometries after minimization by FireFly"));
				strcat(message,("\n\tTo read this file through Gabedit : 'Read/CChemI file'\n\n"));
			}

			free(fileNameGeomFireFly);
		}
		if(fileNamePrefix) free(fileNamePrefix);
	}
	/* minimazation by Gaussian*/
	if(optGaussian)
	{
		char* fileNamePrefix = getSuffixNameFile(fileNameGeom);
		char* keys = strdup_printf("%s Opt",optGaussianMethod);
		if(runGaussianFiles(numberOfGeometries, geometries, energies, fileNamePrefix, keys,orcaCommand) )
		{
			char* fileNameGeomGaussian = strdup_printf("%sGaussian.gab",fileNamePrefix);
			sortGeometries(numberOfGeometries, geometries, energies);
			removeIdenticalGeometries(&numberOfGeometries, &geometries, &energies,tolEnergy,tolDistance);
			if(saveConfoGeometries(numberOfGeometries, geometries, energies, fileNameGeomGaussian))
			{
				createPostProcessingFiles(numberOfGeometriesAll, geometriesAll,energiesAll, fileNameGeomGaussian, mopacKeywordsPost, gaussianKeywordsPost, fireflyKeywordsPost, orcaKeywordsPost, cchemiKeywordsPost, message, mopacCommand, gaussianCommand, fireflyCommand, orcaCommand, cchemiCommand);
				strcat(message,fileNameGeomGaussian);
				strcat(message,("\n\tGeometries after minimization by Gaussian"));
				strcat(message,("\n\tTo read this file through Gabedit : 'Read/CChemI file'\n\n"));
			}

			free(fileNameGeomGaussian);
		}
		if(fileNamePrefix) free(fileNamePrefix);
	}
	/* minimazation by Orca*/
	if(optOrca)
	{
		char* fileNamePrefix = getSuffixNameFile(fileNameGeom);
		char* keys = strdup_printf("%s Opt",optOrcaMethod);
		if(runOrcaFiles(numberOfGeometries, geometries, energies, fileNamePrefix, keys,orcaCommand) )
		{
			char* fileNameGeomOrca = strdup_printf("%sOrca.gab",fileNamePrefix);
			sortGeometries(numberOfGeometries, geometries, energies);
			removeIdenticalGeometries(&numberOfGeometries, &geometries, &energies,tolEnergy,tolDistance);
			if(saveConfoGeometries(numberOfGeometries, geometries, energies, fileNameGeomOrca))
			{
				createPostProcessingFiles(numberOfGeometriesAll, geometriesAll,energiesAll, fileNameGeomOrca, mopacKeywordsPost, gaussianKeywordsPost, fireflyKeywordsPost, orcaKeywordsPost, cchemiKeywordsPost, message, mopacCommand, gaussianCommand, fireflyCommand, orcaCommand, cchemiCommand);
				strcat(message,fileNameGeomOrca);
				strcat(message,("\n\tGeometries after minimization by Orca"));
				strcat(message,("\n\tTo read this file through Gabedit : 'Read/CChemI file'\n\n"));
			}

			free(fileNameGeomOrca);
		}
		if(fileNamePrefix) free(fileNamePrefix);
	}

	fprintf(logfile, "Begin freegeometriesAll\n" );
	fflush(logfile);
	if(geometriesAll && geometriesAll!=geometries)
	{
		for(i=0;i<numberOfGeometriesAll;i++)
			if(geometriesAll[i]) freeForceField(geometriesAll[i]);
		free(geometriesAll);
	}
	fprintf(logfile, "End freegeometriesAll\n" );
	fflush(logfile);
	if(energiesAll && energiesAll!= energies)
	{
		free(energiesAll);
	}
	fprintf(logfile, "End free energiesAll\n" );
	fflush(logfile);

	if(geometries)
	{
		for(i=0;i<numberOfGeometries;i++)
			if(geometries[i]) freeForceField(geometries[i]);
		free(geometries);
	}
	if(energies) free(energies);
	fprintf(logfile, "%s\n",message);
	if(logfile!=stdout) fclose(logfile);
	fclose(file);
	freeForceField(&forceField);
}
/*****************************************************************************/
void molecularMechanicsDynamicsConfoDlg(char* inputFileName)
{
	ForceField forceField; 
	ForceFieldOptions forceFieldOptions;
	MolecularDynamics molecularDynamics;
	int updateFrequency = 1;
	double heatTime;
	double equiTime;
	double runTime;
	double coolTime;
	double heatTemp; 
	double equiTemp; 
	double runTemp; 
	double coolTemp; 
	double stepSize;
	MDIntegratorType integrator = VERLET;
	char* fileNameGeom = NULL;
	char* fileNameTraj = NULL;
	char* fileNameProp = NULL;
	char* mopacKeywordsPost = NULL;
	char* gaussianKeywordsPost = NULL;
	char* fireflyKeywordsPost = NULL;
	char* orcaKeywordsPost = NULL;
	char* cchemiKeywordsPost = NULL;
	double friction=-1;
	double omegaMax = 4000;
	int Nf = 50;
	double collide = 20;
	double qNH = 20;
	MDThermostatType thermostat = NONE;
	int numberOfGeometries = 2;
	ForceField** geometries = NULL; 
	double* energies = NULL;
	boolean optMM = FALSE;
	boolean optMopac = FALSE;
	boolean optGaussian = FALSE;
	boolean optFireFly = FALSE;
	boolean optOrca = FALSE;
	char* optMopacMethod=strdup("PM6");
	char* optGaussianMethod=strdup("AM1");
	char* optFireFlyMethod=strdup("AM1");
	char* optOrcaMethod=strdup("AM1");
	Molecule mol = *(readMolecule(inputFileName,TRUE));

	QuasiNewton quasiNewton;
	ConjugateGradientOptions conjugateGradientOptions;
	SteepestDescent steepestDescent;
	ConjugateGradient conjugateGradient;
	int i;
	char message[BSIZE]="Created files :\n";
	double tolEnergy = -1;
	double tolDistance = -1;
	char* mopacCommand = strdup("mopac");
	char* gaussianCommand=strdup("g03"); 
	char* orcaCommand=strdup("orca"); 
	char* fireflyCommand=strdup("firefly");
	char* cchemiCommand=strdup("cchemi");
	FILE* file = fopen(inputFileName,"rb");
	char* optimizerType= strdup("QuasiNewton");

	setForceFieldOptions(file, &forceFieldOptions);
	mol.klass->buildMMTypes(&mol, file);

	setMDOptions(file, &updateFrequency, 
		&heatTime, &equiTime, &runTime, &coolTime,
		&heatTemp, &runTemp, &equiTemp, &coolTemp, &stepSize, 
		&integrator, &thermostat, &friction,  &omegaMax, &Nf, &collide,&qNH);
	readOneString(file,"mopacCommand",&mopacCommand);
	readOneString(file,"gaussianCommand",&gaussianCommand);
	readOneString(file,"fireflyCommand",&fireflyCommand);
	readOneString(file,"cchemiCommand",&cchemiCommand);
	readOneString(file,"orcaCommand",&orcaCommand);
	readOneReal(file,"tolEnergy",&tolEnergy);
	readOneReal(file,"tolDistance",&tolDistance);
	readOneBoolean(file,"ConfoOptMM",&optMM);
	readOneBoolean(file,"ConfoOptMopac",&optMopac);
	readOneString(file,"ConfoOptMopacMethod",&optMopacMethod);
	readOneBoolean(file,"ConfoOptGaussian",&optGaussian);
	readOneString(file,"ConfoOptGaussianMethod",&optGaussianMethod);
	readOneBoolean(file,"ConfoOptFireFly",&optFireFly);
	readOneString(file,"ConfoOptFireFlyMethod",&optFireFlyMethod);
	readOneBoolean(file,"ConfoOptOrac",&optOrca);
	readOneString(file,"ConfoOptOrcaMethod",&optOrcaMethod);

	readOneString(file,"mopacKeywordsPost",&mopacKeywordsPost);
	readOneString(file,"gaussianKeywordsPost",&gaussianKeywordsPost);
	readOneString(file,"fireflyKeywordsPost",&fireflyKeywordsPost);
	readOneString(file,"cchemiKeywordsPost",&cchemiKeywordsPost);
	readOneString(file,"orcaKeywordsPost",&orcaKeywordsPost);
	optimizerType = setOptOptions(file, &conjugateGradientOptions, &quasiNewton);
	if(strstr(optimizerType,"External")  && (optMopac||optFireFly || optGaussian || optOrca))
	{
		printf("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n");
		printf("Sorry, The optimization using the internal optimizer after a MD conformational search is not yes implemented in this software\n");
		printf("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n");
		exit(1);
	}
	/* number for geometries */
	{
		numberOfGeometries = 10;
		readOneInt(file,"numberOfGeometries",&numberOfGeometries);
		if(numberOfGeometries<2) numberOfGeometries = 2;
	}
	/* fileName for geometries */
	{
		char* suff = getSuffixNameFile(inputFileName);
		fileNameGeom = strdup_printf("%s%s",suff, "Geoms.gab");
		fileNameTraj = strdup_printf("%s%s",suff, "Traj.gab");
		fileNameProp = strdup_printf("%s%s",suff, "Prop.txt");
		free(suff);
	}

/* Optimsation options */ 

	/* Molecule to read */
	if(forceFieldOptions.type==AMBER) forceField = createAmberModel(&mol,forceFieldOptions, stdout);
	else if(forceFieldOptions.type==PAIRWISE) forceField = createPairWiseModel(&mol,forceFieldOptions, stdout);
	setH4CorrectionMM(file, &forceField);
	checkWallCorrection(file, &forceField);


	geometries = runMolecularDynamicsConfo(&molecularDynamics, &forceField,
		updateFrequency, heatTime, equiTime, runTime, heatTemp, equiTemp, runTemp, stepSize, 
		integrator, thermostat, friction, omegaMax, Nf, collide, qNH, numberOfGeometries, fileNameTraj, fileNameProp);
	printf("End runMolecularDynamicsConfo\n");

	freeForceField(&forceField);
	if(geometries && numberOfGeometries>0) energies = malloc(numberOfGeometries*sizeof(double));

	if(geometries && optMM)
	for(i=0;i<numberOfGeometries;i++)
	{
		char* str = NULL;
		energies[i] = 1e30;
		if(!geometries[i]) continue;
		if(str) free(str);
		str = strdup_printf("Minimization of geometry number %d\n", i+1);
		printf("%s",str);
		if(str) free(str);


		if(strstr(optimizerType,"Grad"))
		{
			conjugateGradient.logfile= stdout;
			runConjugateGradient(&conjugateGradient, geometries[i], conjugateGradientOptions); 
			energies[i] = conjugateGradient.forceField->klass->calculateEnergyTmp
				(conjugateGradient.forceField, &conjugateGradient.forceField->molecule );
			freeConjugateGradient(&conjugateGradient);
		}
		else if(strstr(optimizerType,"Quasi"))
		{
			QuasiNewton tmpQuasiNewton = quasiNewton;
			tmpQuasiNewton.forceField = geometries[i];
                	tmpQuasiNewton.logfile = stdout;
                	runQuasiNewton(&tmpQuasiNewton);
			energies[i] = tmpQuasiNewton.forceField->klass->calculateEnergyTmp
				(tmpQuasiNewton.forceField, &tmpQuasiNewton.forceField->molecule );
			freeQuasiNewton(&tmpQuasiNewton);

		}
		else
		{
			steepestDescent.logfile= stdout;
			runSteepestDescent(&steepestDescent, geometries[i],
			       	conjugateGradientOptions.updateFrequency,
			       conjugateGradientOptions.maxIterations,
			       conjugateGradientOptions.gradientNorm,
			       conjugateGradientOptions.maxLines);
			energies[i] = steepestDescent.forceField->klass->calculateEnergyTmp
				(steepestDescent.forceField, &steepestDescent.forceField->molecule );
			freeSteepestDescent(&steepestDescent);
		}
	}
	else 
	{
		for(i=0;i<numberOfGeometries;i++)
		{
			energies[i] = 1e30;
			if(!geometries[i]) continue;
			energies[i] = geometries[i]->klass->calculateEnergyTmp
				(geometries[i], &geometries[i]->molecule );
		}

	}

	/*  sort by energies */
	{
		printf("sortGeometries\n");
		sortGeometries(numberOfGeometries, geometries, energies);
		printf("removeIdenticalGeometries\n");
		removeIdenticalGeometries(&numberOfGeometries, &geometries, &energies,tolEnergy,tolDistance);
	}
	printf("fileNameGeom = %s\n",fileNameGeom);
	if(saveConfoGeometries(numberOfGeometries, geometries, energies, fileNameGeom))
	{
		createPostProcessingFiles(numberOfGeometries, geometries,energies, fileNameGeom, mopacKeywordsPost, gaussianKeywordsPost, fireflyKeywordsPost, orcaKeywordsPost, cchemiKeywordsPost, message, mopacCommand, gaussianCommand, fireflyCommand, orcaCommand, cchemiCommand);
		strcat(message,fileNameGeom);
		strcat(message,("\n\tGeometries selected and optimized using your MM potentials"));
		strcat(message,("\n\tTo read this file through Gabedit: 'Read/CChemI file'\n\n"));
	}
	/* minimazation by mopac*/
	if(optMopac)
	{
		char* fileNamePrefix = getSuffixNameFile(fileNameGeom);
		char* keys = strdup_printf("%s XYZ",optMopacMethod);
		if(runMopacFiles(numberOfGeometries, geometries, energies, fileNamePrefix, keys, mopacCommand) )
		{
			char* fileNameGeomMop = strdup_printf("%sMop.gab",fileNamePrefix);
			sortGeometries(numberOfGeometries, geometries, energies);
			removeIdenticalGeometries(&numberOfGeometries, &geometries, &energies,tolEnergy,tolDistance);
			if(saveConfoGeometries(numberOfGeometries, geometries, energies, fileNameGeomMop))
			{
				createPostProcessingFiles(numberOfGeometries, geometries,energies, fileNameGeomMop, mopacKeywordsPost, gaussianKeywordsPost, fireflyKeywordsPost, orcaKeywordsPost, cchemiKeywordsPost, message, mopacCommand, gaussianCommand, fireflyCommand, orcaCommand, cchemiCommand);
				strcat(message,fileNameGeomMop);
				strcat(message,("\n\tGeometries after minimization by Mopac/"));
				strcat(message,optMopacMethod);
				strcat(message,("\n\tTo read this file through Gabedit : 'Read/CChemI file'\n\n"));
			}

			free(fileNameGeomMop);
		}
		if(fileNamePrefix) free(fileNamePrefix);
	}
	/* minimazation by FireFly*/
	if(optFireFly)
	{
		char* fileNamePrefix = getSuffixNameFile(fileNameGeom);
		char* keys = strdup_printf("RUNTYP=Optimize GBASIS=%s",optFireFlyMethod);
		if(runFireFlyFiles(numberOfGeometries, geometries, energies, fileNamePrefix, keys,fireflyCommand) )
		{
			char* fileNameGeomFireFly = strdup_printf("%sFireFly.gab",fileNamePrefix);
			sortGeometries(numberOfGeometries, geometries, energies);
			removeIdenticalGeometries(&numberOfGeometries, &geometries, &energies,tolEnergy,tolDistance);
			if(saveConfoGeometries(numberOfGeometries, geometries, energies, fileNameGeomFireFly))
			{
				createPostProcessingFiles(numberOfGeometries, geometries,energies, fileNameGeomFireFly, mopacKeywordsPost, gaussianKeywordsPost, fireflyKeywordsPost, orcaKeywordsPost, cchemiKeywordsPost, message, mopacCommand, gaussianCommand, fireflyCommand, orcaCommand, cchemiCommand);
				strcat(message,fileNameGeomFireFly);
				strcat(message,("\n\tGeometries after minimization by FireFly"));
				strcat(message,("\n\tTo read this file through Gabedit : 'Read/CChemI file'\n\n"));
			}

			free(fileNameGeomFireFly);
		}
		if(fileNamePrefix) free(fileNamePrefix);
	}
	/* minimazation by Orca*/
	if(optOrca)
	{
		char* fileNamePrefix = getSuffixNameFile(fileNameGeom);
		char* keys = strdup_printf("%s Opt",optOrcaMethod);
		if(runOrcaFiles(numberOfGeometries, geometries, energies, fileNamePrefix, keys,fireflyCommand) )
		{
			char* fileNameGeomOrca = strdup_printf("%sOrca.gab",fileNamePrefix);
			sortGeometries(numberOfGeometries, geometries, energies);
			removeIdenticalGeometries(&numberOfGeometries, &geometries, &energies,tolEnergy,tolDistance);
			if(saveConfoGeometries(numberOfGeometries, geometries, energies, fileNameGeomOrca))
			{
				createPostProcessingFiles(numberOfGeometries, geometries,energies, fileNameGeomOrca, mopacKeywordsPost, gaussianKeywordsPost, fireflyKeywordsPost, orcaKeywordsPost, cchemiKeywordsPost, message, mopacCommand, gaussianCommand, fireflyCommand, orcaCommand, cchemiCommand);
				strcat(message,fileNameGeomOrca);
				strcat(message,("\n\tGeometries after minimization by Orca"));
				strcat(message,("\n\tTo read this file through Gabedit : 'Read/CChemI file'\n\n"));
			}

			free(fileNameGeomOrca);
		}
		if(fileNamePrefix) free(fileNamePrefix);
	}
	/* minimazation by Gaussian*/
	if(optGaussian)
	{
		char* fileNamePrefix = getSuffixNameFile(fileNameGeom);
		char* keys = strdup_printf("%s Opt",optGaussianMethod);
		if(runGaussianFiles(numberOfGeometries, geometries, energies, fileNamePrefix, keys,fireflyCommand) )
		{
			char* fileNameGeomGaussian = strdup_printf("%sGaussian.gab",fileNamePrefix);
			sortGeometries(numberOfGeometries, geometries, energies);
			removeIdenticalGeometries(&numberOfGeometries, &geometries, &energies,tolEnergy,tolDistance);
			if(saveConfoGeometries(numberOfGeometries, geometries, energies, fileNameGeomGaussian))
			{
				createPostProcessingFiles(numberOfGeometries, geometries,energies, fileNameGeomGaussian, mopacKeywordsPost, gaussianKeywordsPost, fireflyKeywordsPost, orcaKeywordsPost, cchemiKeywordsPost, message, mopacCommand, gaussianCommand, fireflyCommand, orcaCommand, cchemiCommand);
				strcat(message,fileNameGeomGaussian);
				strcat(message,("\n\tGeometries after minimization by Gaussian"));
				strcat(message,("\n\tTo read this file through Gabedit : 'Read/CChemI file'\n\n"));
			}

			free(fileNameGeomGaussian);
		}
		if(fileNamePrefix) free(fileNamePrefix);
	}

	if(geometries)
	{
		for(i=0;i<numberOfGeometries;i++)
			if(geometries[i]) freeForceField(geometries[i]);
		free(geometries);
	}
	if(energies) free(energies);
	printf("%s\n",message);
	fclose(file);
}
/***********************************************************************************************************************/
static ForceField** getMolecularMechanicsRDConfo(ForceField* forceField, int numberOfGeometries)
{
	int i;
	char* str = NULL;
	ForceField** geometries = NULL;

	if(forceField->molecule.nAtoms<1) return NULL;
	if(numberOfGeometries<2) return NULL;
	geometries = malloc(numberOfGeometries*sizeof(ForceField*));
	for (i = 0; i < numberOfGeometries; i++ )
	{
		geometries[i] = NULL;
		if(i>0) forceField->molecule.klass->setRandomPositions(&forceField->molecule);
		forceField->klass->calculateEnergy(forceField);
		if(str) free(str);
		str = strdup_printf(("Geometry # %d Potential energy =  %0.4f"), i+1, forceField->molecule.potentialEnergy);
		printf("%s\n",str);
		geometries[i] = malloc(sizeof(ForceField));
		*geometries[i] = copyForceField(forceField);
	}
	if(str) free(str);
	return geometries;
}
/*****************************************************************************/
void molecularMechanicsRandomConfoDlg(char* inputFileName)
{
	ForceField forceField; 
	ForceFieldOptions forceFieldOptions;
	int updateFrequency = 1;
	double heatTime;
	double equiTime;
	double runTime;
	double coolTime;
	double heatTemp; 
	double equiTemp; 
	double runTemp; 
	double coolTemp; 
	double stepSize;
	MDIntegratorType integrator = VERLET;
	char* fileNameGeom = NULL;
	//char* fileNameTraj = NULL;
	//char* fileNameProp = NULL;
	char* mopacKeywordsPost = NULL;
	char* gaussianKeywordsPost = NULL;
	char* fireflyKeywordsPost = NULL;
	char* orcaKeywordsPost = NULL;
	char* cchemiKeywordsPost = NULL;
	double friction=-1;
	double omegaMax = 4000;
	int Nf = 50;
	double collide = 20;
	double qNH = 20;
	MDThermostatType thermostat = NONE;
	int numberOfGeometries = 2;
	ForceField** geometries = NULL; 
	double* energies = NULL;
	boolean optMM = FALSE;
	boolean optMopac = FALSE;
	boolean optGaussian = FALSE;
	boolean optFireFly = FALSE;
	boolean optOrca = FALSE;
	char* optMopacMethod=strdup("PM6");
	char* optGaussianMethod=strdup("AM1");
	char* optFireFlyMethod=strdup("AM1");
	char* optOrcaMethod=strdup("AM1");
	Molecule mol = *(readMolecule(inputFileName,TRUE));

	QuasiNewton quasiNewton;
	ConjugateGradientOptions conjugateGradientOptions;
	SteepestDescent steepestDescent;
	ConjugateGradient conjugateGradient;
	int i;
	char message[BSIZE]="Created files :\n";
	double tolEnergy = -1;
	double tolDistance = -1;
	char* mopacCommand = strdup("mopac");
	char* gaussianCommand=strdup("g03"); 
	char* orcaCommand=strdup("orca"); 
	char* fireflyCommand=strdup("firefly");
	char* cchemiCommand=strdup("cchemi");
	FILE* file = fopen(inputFileName,"rb");
	char* optimizerType= strdup("QuasiNewton");

	setForceFieldOptions(file, &forceFieldOptions);
	mol.klass->buildMMTypes(&mol, file);

	setMDOptions(file, &updateFrequency, 
		&heatTime, &equiTime, &runTime, &coolTime,
		&heatTemp, &runTemp, &equiTemp, &coolTemp, &stepSize, 
		&integrator, &thermostat, &friction,  &omegaMax, &Nf, &collide,&qNH);
	readOneString(file,"mopacCommand",&mopacCommand);
	readOneString(file,"gaussianCommand",&gaussianCommand);
	readOneString(file,"fireflyCommand",&fireflyCommand);
	readOneString(file,"cchemiCommand",&cchemiCommand);
	readOneString(file,"orcaCommand",&orcaCommand);
	readOneReal(file,"tolEnergy",&tolEnergy);
	readOneReal(file,"tolDistance",&tolDistance);
	readOneBoolean(file,"ConfoOptMM",&optMM);
	readOneBoolean(file,"ConfoOptMopac",&optMopac);
	readOneString(file,"ConfoOptMopacMethod",&optMopacMethod);
	readOneBoolean(file,"ConfoOptGaussian",&optGaussian);
	readOneString(file,"ConfoOptGaussianMethod",&optGaussianMethod);
	readOneBoolean(file,"ConfoOptFireFly",&optFireFly);
	readOneString(file,"ConfoOptFireFlyMethod",&optFireFlyMethod);
	readOneBoolean(file,"ConfoOptOrac",&optOrca);
	readOneString(file,"ConfoOptOrcaMethod",&optOrcaMethod);

	readOneString(file,"mopacKeywordsPost",&mopacKeywordsPost);
	readOneString(file,"gaussianKeywordsPost",&gaussianKeywordsPost);
	readOneString(file,"fireflyKeywordsPost",&fireflyKeywordsPost);
	readOneString(file,"cchemiKeywordsPost",&cchemiKeywordsPost);
	readOneString(file,"orcaKeywordsPost",&orcaKeywordsPost);
	optimizerType = setOptOptions(file, &conjugateGradientOptions, &quasiNewton);
	if(strstr(optimizerType,"External")  && (optMopac||optFireFly || optGaussian || optOrca))
	{
		printf("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n");
		printf("Sorry, The optimization using the internal optimizer after a MD conformational search is not yes implemented in this software\n");
		printf("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n");
		exit(1);
	}
	/* number for geometries */
	{
		numberOfGeometries = 10;
		readOneInt(file,"numberOfGeometries",&numberOfGeometries);
		if(numberOfGeometries<2) numberOfGeometries = 2;
	}
	/* fileName for geometries */
	{
		char* suff = getSuffixNameFile(inputFileName);
		fileNameGeom = strdup_printf("%s%s",suff, "Geoms.gab");
		//fileNameTraj = strdup_printf("%s%s",suff, "Traj.gab");
		//fileNameProp = strdup_printf("%s%s",suff, "Prop.txt");
		free(suff);
	}

/* Optimsation options */ 

	/* Molecule to read */
	if(forceFieldOptions.type==AMBER) forceField = createAmberModel(&mol,forceFieldOptions, stdout);
	else if(forceFieldOptions.type==PAIRWISE) forceField = createPairWiseModel(&mol,forceFieldOptions, stdout);
	setH4CorrectionMM(file, &forceField);
	checkWallCorrection(file, &forceField);


	geometries = getMolecularMechanicsRDConfo(&forceField, numberOfGeometries);
	printf("End getMolecularMechanicsRDConfo\n");

	freeForceField(&forceField);
	if(geometries && numberOfGeometries>0) energies = malloc(numberOfGeometries*sizeof(double));

	if(geometries && optMM)
	for(i=0;i<numberOfGeometries;i++)
	{
		char* str = NULL;
		energies[i] = 1e30;
		if(!geometries[i]) continue;
		if(str) free(str);
		str = strdup_printf("Minimization of geometry number %d\n", i+1);
		printf("%s",str);
		if(str) free(str);


		if(strstr(optimizerType,"Grad"))
		{
			conjugateGradient.logfile= stdout;
			runConjugateGradient(&conjugateGradient, geometries[i], conjugateGradientOptions); 
			energies[i] = conjugateGradient.forceField->klass->calculateEnergyTmp
				(conjugateGradient.forceField, &conjugateGradient.forceField->molecule );
			freeConjugateGradient(&conjugateGradient);
		}
		else if(strstr(optimizerType,"Quasi"))
		{
			QuasiNewton tmpQuasiNewton = quasiNewton;
			tmpQuasiNewton.forceField = geometries[i];
                	tmpQuasiNewton.logfile = stdout;
                	runQuasiNewton(&tmpQuasiNewton);
			energies[i] = tmpQuasiNewton.forceField->klass->calculateEnergyTmp
				(tmpQuasiNewton.forceField, &tmpQuasiNewton.forceField->molecule );
			freeQuasiNewton(&tmpQuasiNewton);

		}
		else
		{
			steepestDescent.logfile= stdout;
			runSteepestDescent(&steepestDescent, geometries[i],
			       	conjugateGradientOptions.updateFrequency,
			       conjugateGradientOptions.maxIterations,
			       conjugateGradientOptions.gradientNorm,
			       conjugateGradientOptions.maxLines);
			energies[i] = steepestDescent.forceField->klass->calculateEnergyTmp
				(steepestDescent.forceField, &steepestDescent.forceField->molecule );
			freeSteepestDescent(&steepestDescent);
		}
	}
	else 
	{
		for(i=0;i<numberOfGeometries;i++)
		{
			energies[i] = 1e30;
			if(!geometries[i]) continue;
			energies[i] = geometries[i]->klass->calculateEnergyTmp
				(geometries[i], &geometries[i]->molecule );
		}

	}

	/*  sort by energies */
	{
		printf("sortGeometries\n");
		sortGeometries(numberOfGeometries, geometries, energies);
		printf("removeIdenticalGeometries\n");
		removeIdenticalGeometries(&numberOfGeometries, &geometries, &energies,tolEnergy,tolDistance);
	}
	printf("fileNameGeom = %s\n",fileNameGeom);
	if(saveConfoGeometries(numberOfGeometries, geometries, energies, fileNameGeom))
	{
		createPostProcessingFiles(numberOfGeometries, geometries,energies, fileNameGeom, mopacKeywordsPost, gaussianKeywordsPost, fireflyKeywordsPost, orcaKeywordsPost, cchemiKeywordsPost, message, mopacCommand, gaussianCommand, fireflyCommand, orcaCommand, cchemiCommand);
		strcat(message,fileNameGeom);
		strcat(message,("\n\tGeometries selected and optimized using your MM potentials"));
		strcat(message,("\n\tTo read this file through Gabedit: 'Read/CChemI file'\n\n"));
	}
	/* minimazation by mopac*/
	if(optMopac)
	{
		char* fileNamePrefix = getSuffixNameFile(fileNameGeom);
		char* keys = strdup_printf("%s XYZ",optMopacMethod);
		if(runMopacFiles(numberOfGeometries, geometries, energies, fileNamePrefix, keys, mopacCommand) )
		{
			char* fileNameGeomMop = strdup_printf("%sMop.gab",fileNamePrefix);
			sortGeometries(numberOfGeometries, geometries, energies);
			removeIdenticalGeometries(&numberOfGeometries, &geometries, &energies,tolEnergy,tolDistance);
			if(saveConfoGeometries(numberOfGeometries, geometries, energies, fileNameGeomMop))
			{
				createPostProcessingFiles(numberOfGeometries, geometries,energies, fileNameGeomMop, mopacKeywordsPost, gaussianKeywordsPost, fireflyKeywordsPost, orcaKeywordsPost, cchemiKeywordsPost, message, mopacCommand, gaussianCommand, fireflyCommand, orcaCommand, cchemiCommand);
				strcat(message,fileNameGeomMop);
				strcat(message,("\n\tGeometries after minimization by Mopac/"));
				strcat(message,optMopacMethod);
				strcat(message,("\n\tTo read this file through Gabedit : 'Read/CChemI file'\n\n"));
			}

			free(fileNameGeomMop);
		}
		if(fileNamePrefix) free(fileNamePrefix);
	}
	/* minimazation by FireFly*/
	if(optFireFly)
	{
		char* fileNamePrefix = getSuffixNameFile(fileNameGeom);
		char* keys = strdup_printf("RUNTYP=Optimize GBASIS=%s",optFireFlyMethod);
		if(runFireFlyFiles(numberOfGeometries, geometries, energies, fileNamePrefix, keys,fireflyCommand) )
		{
			char* fileNameGeomFireFly = strdup_printf("%sFireFly.gab",fileNamePrefix);
			sortGeometries(numberOfGeometries, geometries, energies);
			removeIdenticalGeometries(&numberOfGeometries, &geometries, &energies,tolEnergy,tolDistance);
			if(saveConfoGeometries(numberOfGeometries, geometries, energies, fileNameGeomFireFly))
			{
				createPostProcessingFiles(numberOfGeometries, geometries,energies, fileNameGeomFireFly, mopacKeywordsPost, gaussianKeywordsPost, fireflyKeywordsPost, orcaKeywordsPost, cchemiKeywordsPost, message, mopacCommand, gaussianCommand, fireflyCommand, orcaCommand, cchemiCommand);
				strcat(message,fileNameGeomFireFly);
				strcat(message,("\n\tGeometries after minimization by FireFly"));
				strcat(message,("\n\tTo read this file through Gabedit : 'Read/CChemI file'\n\n"));
			}

			free(fileNameGeomFireFly);
		}
		if(fileNamePrefix) free(fileNamePrefix);
	}
	/* minimazation by Orca*/
	if(optOrca)
	{
		char* fileNamePrefix = getSuffixNameFile(fileNameGeom);
		char* keys = strdup_printf("%s Opt",optOrcaMethod);
		if(runOrcaFiles(numberOfGeometries, geometries, energies, fileNamePrefix, keys,fireflyCommand) )
		{
			char* fileNameGeomOrca = strdup_printf("%sOrca.gab",fileNamePrefix);
			sortGeometries(numberOfGeometries, geometries, energies);
			removeIdenticalGeometries(&numberOfGeometries, &geometries, &energies,tolEnergy,tolDistance);
			if(saveConfoGeometries(numberOfGeometries, geometries, energies, fileNameGeomOrca))
			{
				createPostProcessingFiles(numberOfGeometries, geometries,energies, fileNameGeomOrca, mopacKeywordsPost, gaussianKeywordsPost, fireflyKeywordsPost, orcaKeywordsPost, cchemiKeywordsPost, message, mopacCommand, gaussianCommand, fireflyCommand, orcaCommand, cchemiCommand);
				strcat(message,fileNameGeomOrca);
				strcat(message,("\n\tGeometries after minimization by Orca"));
				strcat(message,("\n\tTo read this file through Gabedit : 'Read/CChemI file'\n\n"));
			}

			free(fileNameGeomOrca);
		}
		if(fileNamePrefix) free(fileNamePrefix);
	}
	/* minimazation by Gaussian*/
	if(optGaussian)
	{
		char* fileNamePrefix = getSuffixNameFile(fileNameGeom);
		char* keys = strdup_printf("%s Opt",optGaussianMethod);
		if(runGaussianFiles(numberOfGeometries, geometries, energies, fileNamePrefix, keys,fireflyCommand) )
		{
			char* fileNameGeomGaussian = strdup_printf("%sGaussian.gab",fileNamePrefix);
			sortGeometries(numberOfGeometries, geometries, energies);
			removeIdenticalGeometries(&numberOfGeometries, &geometries, &energies,tolEnergy,tolDistance);
			if(saveConfoGeometries(numberOfGeometries, geometries, energies, fileNameGeomGaussian))
			{
				createPostProcessingFiles(numberOfGeometries, geometries,energies, fileNameGeomGaussian, mopacKeywordsPost, gaussianKeywordsPost, fireflyKeywordsPost, orcaKeywordsPost, cchemiKeywordsPost, message, mopacCommand, gaussianCommand, fireflyCommand, orcaCommand, cchemiCommand);
				strcat(message,fileNameGeomGaussian);
				strcat(message,("\n\tGeometries after minimization by Gaussian"));
				strcat(message,("\n\tTo read this file through Gabedit : 'Read/CChemI file'\n\n"));
			}

			free(fileNameGeomGaussian);
		}
		if(fileNamePrefix) free(fileNamePrefix);
	}

	if(geometries)
	{
		for(i=0;i<numberOfGeometries;i++)
			if(geometries[i]) freeForceField(geometries[i]);
		free(geometries);
	}
	if(energies) free(energies);
	printf("%s\n",message);
	fclose(file);
}
/*****************************************************************************/
void molecularMechanicsDynamicsDlg(char* inputFileName)
{
	ForceField forceField; 
	ForceFieldOptions forceFieldOptions;
	MolecularDynamics molecularDynamics;
	int updateFrequency = 1;
	double heatTime;
	double equiTime;
	double runTime;
	double coolTime; 
	double heatTemp; 
	double equiTemp; 
	double runTemp; 
	double coolTemp; 
	double stepSize;
	MDIntegratorType integrator = VERLET;
	char* fileNameTraj = NULL;
	char* fileNameProp = NULL;
	double friction=-1;
	double omegaMax = 4000;
	int Nf = 50;
	double collide = 20;
	double qNH = 20;
	Molecule* mol = readMolecule(inputFileName,TRUE);
	MDThermostatType thermostat = NONE;
	FILE* file = fopen(inputFileName,"rb");

	setForceFieldOptions(file, &forceFieldOptions);
	mol->klass->buildMMTypes(mol, file);

	setMDOptions(file, &updateFrequency, 
		&heatTime, &equiTime, &runTime, &coolTime,
		&heatTemp, &runTemp, &equiTemp, &coolTemp, &stepSize, 
		&integrator, &thermostat, &friction,  &omegaMax, &Nf, &collide,&qNH);

	{
		char* suff = getSuffixNameFile(inputFileName);
		fileNameTraj = strdup_printf("%s%s",suff, "Traj.gab");
		fileNameProp = strdup_printf("%s%s",suff, "Prop.txt");
		free(suff);
	}

	if(forceFieldOptions.type==AMBER) forceField = createAmberModel(mol,forceFieldOptions, stdout);
	else if(forceFieldOptions.type==PAIRWISE) forceField = createPairWiseModel(mol,forceFieldOptions, stdout);
	setH4CorrectionMM(file, &forceField);
	checkWallCorrection(file, &forceField);

	runMolecularDynamics(&molecularDynamics, &forceField,
		updateFrequency, heatTime, equiTime, runTime, coolTime, heatTemp, equiTemp, runTemp, coolTemp, stepSize, 
		integrator, thermostat, friction, omegaMax, Nf, collide, qNH, fileNameTraj, fileNameProp);

	freeForceField(&forceField);
	fclose(file);
}
/********************************************************************************/
void molecularMechanicsMinimizeDlg(char* inputFileName)
{
	ForceField forceField; 
	ForceFieldOptions forceFieldOptions;
	SteepestDescent steepestDescent;
	ConjugateGradient conjugateGradient;
	QuasiNewton quasiNewton;
	ConjugateGradientOptions conjugateGradientOptions;
	char* optimizerType = NULL;
	Molecule* mol = readMolecule(inputFileName,TRUE);
	FILE* file = fopen(inputFileName,"rb");
	char* fileNameOut = strdup_printf("%sOpt.gab",getSuffixNameFile(inputFileName));

	setForceFieldOptions(file, &forceFieldOptions);
	mol->klass->buildMMTypes(mol, file);
	optimizerType = setOptOptions(file, &conjugateGradientOptions, &quasiNewton);


/* Molecule to read */
	if(forceFieldOptions.type==AMBER) forceField = createAmberModel(mol,forceFieldOptions, stdout);
	else if(forceFieldOptions.type==PAIRWISE) forceField = createPairWiseModel(mol,forceFieldOptions,stdout);
	setH4CorrectionMM(file, &forceField);
	checkWallCorrection(file, &forceField);

	if(strstr(optimizerType,"Grad"))
	{
		printf("Minimization by Conjugate Gradient method\n");
		conjugateGradient.logfile= stdout;
		runConjugateGradient(&conjugateGradient, &forceField, conjugateGradientOptions); 
		printf("Optimized geometry saved in %s file\n",fileNameOut);
		forceField.molecule.klass->computeDipole(&forceField.molecule);
		forceField.molecule.klass->save(&forceField.molecule, fileNameOut);
		freeConjugateGradient(&conjugateGradient);
	}
	else if(strstr(optimizerType,"Quasi"))
	{
		printf("Minimization by QuasiNewton method\n");
		quasiNewton.forceField = &forceField; 
               	quasiNewton.logfile = stdout;
                runQuasiNewton(&quasiNewton);
		printf("Optimized geometry saved in %s file\n",fileNameOut);
		forceField.molecule.klass->computeDipole(&forceField.molecule);
		forceField.molecule.klass->save(&forceField.molecule, fileNameOut);
		freeQuasiNewton(&quasiNewton);
	}
	else
	{
		printf("Minimization by steepest descent method\n");
		steepestDescent.logfile= stdout;
		runSteepestDescent(&steepestDescent, &forceField,
			       	conjugateGradientOptions.updateFrequency,
			       conjugateGradientOptions.maxIterations,
			       conjugateGradientOptions.gradientNorm,
			       conjugateGradientOptions.maxLines);
		printf("Optimized geometry saved in %s file\n",fileNameOut);
		forceField.molecule.klass->computeDipole(&forceField.molecule);
		forceField.molecule.klass->save(&forceField.molecule, fileNameOut);
		freeSteepestDescent(&steepestDescent);
	}
	printEnergyAndGradient(&forceField);
	freeForceField(&forceField);
	fclose(file);
	free(fileNameOut);
}
/*****************************************************************************/
void molecularMechanicsEnergyDlg(char* inputFileName)
{
	ForceField forceField;
	ForceFieldOptions forceFieldOptions;
	Molecule* mol = readMolecule(inputFileName,TRUE);
	FILE* file = fopen(inputFileName,"rb");
	char* fileNameOut = strdup_printf("%s.gab",getSuffixNameFile(inputFileName));

	setForceFieldOptions(file, &forceFieldOptions);
	mol->klass->buildMMTypes(mol, file);

/* Molecule to read */
	if(forceFieldOptions.type==AMBER) forceField = createAmberModel(mol,forceFieldOptions, stdout);
	else if(forceFieldOptions.type==PAIRWISE) forceField = createPairWiseModel(mol,forceFieldOptions,stdout);
	setH4CorrectionMM(file, &forceField);
	checkWallCorrection(file, &forceField);

	printEnergyAndGradient(&forceField);
	printf("Geometry saved in %s file\n",fileNameOut);
	forceField.molecule.klass->computeDipole(&forceField.molecule);
	forceField.molecule.klass->save(&forceField.molecule, fileNameOut);
	freeForceField(&forceField);
	fclose(file);
	free(fileNameOut);
}
/*****************************************************************************/
void molecularMechanicsGradientDlg(char* inputFileName)
{
	ForceField forceField;
	ForceFieldOptions forceFieldOptions;
	Molecule* mol = readMolecule(inputFileName,TRUE);
	FILE* file = fopen(inputFileName,"rb");
	char* fileNameOut = strdup_printf("%s.gab",getSuffixNameFile(inputFileName));

	setForceFieldOptions(file, &forceFieldOptions);
	mol->klass->buildMMTypes(mol, file);

/* Molecule to read */
	if(forceFieldOptions.type==AMBER) forceField = createAmberModel(mol,forceFieldOptions, stdout);
	else if(forceFieldOptions.type==PAIRWISE) forceField = createPairWiseModel(mol,forceFieldOptions,stdout);
	setH4CorrectionMM(file, &forceField);
	checkWallCorrection(file, &forceField);

	printEnergyAndGradient(&forceField);
	printf("Geometry saved in %s file\n",fileNameOut);
	forceField.molecule.klass->computeDipole(&forceField.molecule);
	forceField.molecule.klass->save(&forceField.molecule, fileNameOut);
	freeForceField(&forceField);
	fclose(file);
	free(fileNameOut);
}
/*****************************************************************************/
void molecularMechanicsFrequenciesDlg(char* inputFileName)
{
	ForceField forceField;
	ForceFieldOptions forceFieldOptions;
	Molecule* mol = readMolecule(inputFileName,TRUE);
	FILE* file = fopen(inputFileName,"rb");
	char* fileNameOut = strdup_printf("%sFreq.gab",getSuffixNameFile(inputFileName));
	double* frequencies = NULL;
	double* reducedMasses = NULL;
	double* IRIntensities = NULL;
	double** modes = NULL;
	int nModes = 0;
	int i;

	setForceFieldOptions(file, &forceFieldOptions);
	mol->klass->buildMMTypes(mol, file);

/* Molecule to read */
	if(forceFieldOptions.type==AMBER) forceField = createAmberModel(mol,forceFieldOptions, stdout);
	else if(forceFieldOptions.type==PAIRWISE) forceField = createPairWiseModel(mol,forceFieldOptions,stdout);
	printf("Appel de setH4CorrectionMM\n");
	setH4CorrectionMM(file, &forceField);
	printf("Appel de checkWallCorrection\n");
	checkWallCorrection(file, &forceField);
	fclose(file);

	//forceField.klass->calculateEnergy(&forceField);//alrady calculated in calculateGradient
	//forceField.klass->calculateGradient(&forceField); // already in printEnergyAndGrad

	printf("Appelr de computeMMFrequencies\n");
	nModes = computeMMFrequencies(&forceField, &frequencies, & modes, &reducedMasses, &IRIntensities);
	

	printEnergyAndGradient(&forceField);
	printf("Frequencies and modes in the %s file\n",fileNameOut);
	mol->klass->saveFrequencies(mol, fileNameOut, nModes, frequencies, modes, reducedMasses, IRIntensities);
	addHarmonicVelocities(inputFileName, nModes, frequencies, modes, reducedMasses, IRIntensities);
	printHarmonicVelocities(inputFileName, nModes, frequencies, modes, reducedMasses);
	if(frequencies) free(frequencies);
	if(reducedMasses) free(reducedMasses);
	for(i=0;i<nModes;i++) free(modes[i]);
	if(modes) free(modes);
	freeForceField(&forceField);
}
/*****************************************************************************/
void molecularMechanicsOptFrequenciesDlg(char* inputFileName)
{
	ForceField forceField;
	ForceFieldOptions forceFieldOptions;
	Molecule* mol;
	FILE* file;
	char* fileNameFreq = strdup_printf("%sFreq.gab",getSuffixNameFile(inputFileName));
	char* fileNameOpt = strdup_printf("%sOpt.gab",getSuffixNameFile(inputFileName));
	double* frequencies = NULL;
	double* reducedMasses = NULL;
	double* IRIntensities = NULL;
	double** modes = NULL;
	int nModes = 0;
	int i;

	molecularMechanicsMinimizeDlg(inputFileName);

	file = fopen(inputFileName,"rb");
	setForceFieldOptions(file, &forceFieldOptions);
	mol = readMoleculeFromGabeditFile(fileNameOpt);
	mol->klass->setConnections(mol);
	mol->klass->buildMMTypes(mol, file);
	//mol = readMolecule(inputFileName,TRUE);

/* Molecule to read */
	if(forceFieldOptions.type==AMBER) forceField = createAmberModel(mol,forceFieldOptions, stdout);
	else if(forceFieldOptions.type==PAIRWISE) forceField = createPairWiseModel(mol,forceFieldOptions,stdout);
	setH4CorrectionMM(file, &forceField);
	checkWallCorrection(file, &forceField);
	fclose(file);

	nModes = computeMMFrequencies(&forceField, &frequencies, & modes, &reducedMasses, &IRIntensities);

	printEnergyAndGradient(&forceField);
	printf("Frequencies and modes in the %s file\n",fileNameFreq);
	mol->klass->saveFrequencies(mol, fileNameFreq, nModes, frequencies, modes, reducedMasses, IRIntensities);
	addHarmonicVelocities(inputFileName, nModes, frequencies, modes, reducedMasses, IRIntensities);
	printHarmonicVelocities(inputFileName, nModes, frequencies, modes, reducedMasses);
	if(frequencies) free(frequencies);
	if(reducedMasses) free(reducedMasses);
	for(i=0;i<nModes;i++) free(modes[i]);
	if(modes) free(modes);
	freeForceField(&forceField);
	free(fileNameOpt);
	free(fileNameFreq);
}
