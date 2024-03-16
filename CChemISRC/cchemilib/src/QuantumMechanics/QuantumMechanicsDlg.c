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

/* QuantumMechanicsDlg.c */

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
#include "../QuantumMechanics/QuantumMechanicsModel.h"
#include "../QuantumMechanics/QuantumMechanics.h"
#include "../QuantumMechanics/QuantumMechanicsMD.h"
#include "../QuantumMechanics/QuantumMechanicsDlg.h"
#include "../QuantumMechanics/SteepestDescentQM.h"
#include "../QuantumMechanics/ConjugateGradientQM.h"
#include "../QuantumMechanics/QuasiNewtonQM.h"
#include "../EmpriricalCorrections/HydrogenBondCorrection.h"
#include "../EmpriricalCorrections/ShortRangeBasisSetCorrection.h"
#include "../EmpriricalCorrections/DispersionCorrection.h"
#include "../QuantumMechanics/GeneticAlgorithm.h"
#include "../QuantumMechanics/QuantumMechanicsDlg.h"

typedef enum
{
	TOLE = 0,
	TOLD = 1
} TOLptions;

#define NINTEGOPTIONS 3
#define NTHERMOPTIONS 3

#define NENTRYTOL 2

#define NCONSTRAINTS 3

static boolean saveGeometry(Molecule* molecule, double energy, char* fileNameGeom);
static boolean miminizeGeometriesUsingInternalOptimizer(int numberOfGeometries, QuantumMechanicsModel** geometries, double* energies, char* inputFileName);
/*****************************************************************************/
static void saveEnergies(double* energies, int n, char* inputFileName)
{
	FILE* file=NULL;
	int i;
	char* suff ;
	char* fileName ;

	if(!energies) return;

	suff = getSuffixNameFile(inputFileName);
	fileName = strdup_printf("%s%s",suff, "Energies.txt");
	free(suff);

	file=fopen(fileName,"w");
	if(!file) return;

	for(i=0;i<n;i++)
	{
		fprintf(file,"%lf %d\n", energies[i],i+1);
	}
	fclose(file);
	printf("----------------------------------------- \n");
	printf("Energies for all geometries saved in %s file\n",fileName);
	printf("----------------------------------------- \n");
	fflush(stdout);
	free(fileName);
}
static void checkWallCorrection(FILE* file, QuantumMechanicsModel* qmModel)
{
	char t[BSIZE];
	char* pos;
	rewind(file);
	qmModel->addWallCorrection=FALSE;
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
			if(n==3 && nc%2==0) { qmModel->addWallCorrection=TRUE; break;}
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
static void printEnergyAndGradient(QuantumMechanicsModel* qmModel)
{
	char* str;
	double gradientNorm = 0;
	int i,j;

	qmModel->klass->calculateGradient(qmModel);

	gradientNorm = 0;
	for (  i = 0; i < qmModel->molecule.nAtoms; i++ )
		for(j=0;j<3;j++)
			gradientNorm += 
			qmModel->molecule.atoms[i].gradient[j]
			*qmModel->molecule.atoms[i].gradient[j]; 

	str = strdup_printf(("Gradient Norm  = %0.14f energy = %0.14f(kcal/mol) %0.14f(Hartree)\n"),
		sqrt(gradientNorm),qmModel->molecule.potentialEnergy, qmModel->molecule.potentialEnergy/(AUTOKCAL)); 

	printf("%s",str);
	free(str);
}
/*****************************************************************************/
static void printEnergy(QuantumMechanicsModel* qmModel)
{
	char* str;

	qmModel->klass->calculateEnergy(qmModel);
	str = strdup_printf(("energy = %0.14f(kcal/mol) %0.14f(Hartree)\n"),qmModel->molecule.potentialEnergy, qmModel->molecule.potentialEnergy/(AUTOKCAL)); 
	printf("%s",str);
	free(str);
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
static boolean getEnergyMopac(char* fileNameOut, double* energy)
{
	FILE* file = NULL;
	char buffer[1024];
	char* pdest = NULL;

 	file = fopen(fileNameOut, "r");
	*energy=1e10;
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
				int l = strlen(pdest);
				int i;
				for(i=0;i<l;i++) if(pdest[i]=='D' || pdest[i]=='E') pdest[i] ='E';
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
/*************************************************************************************************************/
static boolean putMopacMoleculeInFile(Molecule* mol, FILE* file)
{
        char buffer[BSIZE];
	int i;
	int k1 = 0;
	int k2 = 0;
	int k3 = 0;

	if(mol->nAtoms<1) return FALSE;
      	for (i=0;i<mol->nAtoms;i++)
	{
		k1 = k2 = k3 = 0;
		if(mol->atoms[i].variable) k1 = k2 = k3 = 1;
		sprintf(buffer,"%s  %20.14f %d %20.14f %d %20.14f %d\n",mol->atoms[i].prop.symbol,
				mol->atoms[i].coordinates[0], k1,
				mol->atoms[i].coordinates[1], k2,
				mol->atoms[i].coordinates[2], k3
				);
       		fprintf(file, "%s",buffer);
	}
	return TRUE;
}
/*****************************************************************************/
static boolean runOneMopac(Molecule* mol, char* fileNamePrefix, char* keyWords, char* mopacCommand)
{
	FILE* file = NULL;
	FILE* fileSH = NULL;
	char* fileNameIn = NULL;
	char* fileNameOut = NULL;
	char* fileNameSH = NULL;
	char multiplicityStr[100];
	char buffer[1024];
	double energy = 0;
#ifdef OS_WIN32
	char c='%';
#endif

	if(!mol) return FALSE;
	if(mol->nAtoms<1) return FALSE;
#ifndef OS_WIN32
	fileNameSH = strdup_printf("%sMopacOne.sh",fileNamePrefix);
#else
	fileNameSH = strdup_printf("%sMopacOne.bat",fileNamePrefix);
#endif
 	fileSH = fopen(fileNameSH, "w");
	if(!fileSH) return FALSE;
#ifdef OS_WIN32
	fprintf(fileSH,"@echo off\n");
	fprintf(fileSH,"set PATH=%cPATH%c;\"%s\"\n",c,c,mopacDirectory);
#endif

	getMultiplicityName(mol->spinMultiplicity, multiplicityStr);

	fileNameIn = strdup_printf("%sOne.mop",fileNamePrefix);
 	file = fopen(fileNameIn, "w");
	if(!file) 
	{
 		if(fileNameIn) free(fileNameIn);
 		if(fileNameOut) free(fileNameOut);
 		if(fileNameSH) free(fileNameSH);
		return FALSE;
	}
	if(mol->spinMultiplicity>1)
	fprintf(file,"%s UHF CHARGE=%d %s\n",keyWords,mol->totalCharge,multiplicityStr);
	else
	fprintf(file,"%s CHARGE=%d %s\n",keyWords,mol->totalCharge,multiplicityStr);
	fprintf(file,"\n");
	fprintf(file,"Mopac file generated by Gabedit\n");

  	putMopacMoleculeInFile(mol, file);
	fclose(file);
	{
		char* str = NULL;
		if(!strstr(keyWords,"1SCF") && !strstr(keyWords,"ESP")) str = strdup_printf("Minimization by Mopac/%s ... Please wait",keyWords);
		else if(strstr(keyWords,"ESP")) str = strdup_printf("ESP charges from Mopac/%s ... Please wait",keyWords);
		else str = strdup_printf("Computing of energy by Mopac/%s ... Please wait",keyWords);
		printf("%s\n",str);
		if(str) free(str);
	}
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
	if(getEnergyMopac(fileNameOut,&energy))
	{
		char* str = NULL;
		mol->klass->readGeomFromMopacOutputFile(mol, fileNameOut, -1);
		str = strdup_printf("Energy by Mopac = %f", energy);
		printf("%s\n",str);
		if(str) free(str);
		if(strstr(keyWords,"XYZ"))
		{
			str = strdup_printf("%s.gab",fileNamePrefix);
			saveGeometry(mol, energy, str);
			printf("----------------------------------------- \n");
			printf("Optimized geometry saved in %s file\n",str);
			printf("----------------------------------------- \n");
			if(str) free(str);
		}
		else
		{
			str = strdup_printf("%s.gab",fileNamePrefix);
			saveGeometry(mol, energy, str);
			printf("----------------------------------------- \n");
			printf("Geometry saved in %s file\n",str);
			printf("----------------------------------------- \n");
			if(str) free(str);
		}
	}
	else
	{
		char* str = NULL;
		str = strdup_printf(
				(
				"Sorry, I cannot read the output file : %s "
				" ; Check also the installation of Mopac..."
				),
				fileNameOut
				);
		printf("%s",str);
		if(str) free(str);
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
static boolean getEnergyFireFly(char* fileNameOut, double* energy)
{
	FILE* file = NULL;
	char buffer[1024];
	char* pdest = NULL;
	boolean OK = FALSE;

	*energy=1e10;
 	file = fopen(fileNameOut, "r");
	if(!file) return FALSE;
	 while(!feof(file))
	 {
		if(!fgets(buffer,BSIZE,file))break;
		pdest = strstr( buffer, "HEAT OF FORMATION IS");
		if(!pdest) pdest = strstr( buffer, "FINAL ENERGY IS");
		if(pdest) 
		{
			pdest = strstr( buffer, "S");
			if(pdest)
			{
				int l = strlen(pdest);
				int i;
				for(i=0;i<l;i++) if(pdest[i]=='D' || pdest[i]=='E') pdest[i] ='E';
				if(sscanf(pdest+1,"%lf",energy)==1) OK = TRUE;
				if(OK && strstr( buffer, "FINAL ENERGY IS")) *energy *= AUTOKCAL;
			}
		}
	 }
	fclose(file);
	return OK;
}
/*************************************************************************************************************/
static void putFireFlyMoleculeXYZFixed(Molecule* mol, FILE* file)
{
	int i,k,l;
	int nvar = 0;

        if(mol->nAtoms<2)return;
	nvar = 0;
        for(i=0;i<mol->nAtoms;i++)
		if(mol->atoms[i].variable) nvar+=3;
	/* printf("nvar = %d\n",nvar);*/
	if(nvar==3*mol->nAtoms) return;
	if(nvar==0) return;

        fprintf(file," ");
        fprintf(file, "$STATPT\n");
        fprintf (file,"   IFREEZ(1)=");

	l = 0;
        for(i=0;i<mol->nAtoms;i++)
	{
		if(!mol->atoms[i].variable)
		{
			l++;
			k = i*3+1;
			fprintf(file,"%d, %d, %d ",k,k+1,k+2);
			if(l%10==0) fprintf(file,"\n");
		}
	}
	fprintf(file,"\n ");
        fprintf (file, "$END\n");
}
/*****************************************************************************/
static boolean runOneFireFly(Molecule* mol, char* fileNamePrefix, char* keyWords, char* fireflyCommand)
{
	FILE* file = NULL;
	FILE* fileSH = NULL;
	int j;
	char* fileNameIn = NULL;
	char* fileNameOut = NULL;
	char* fileNameSH = NULL;
	char multiplicityStr[100];
	char buffer[1024];
	double energy = 0;
#ifdef OS_WIN32
	char c='%';
#endif

	if(!mol) return FALSE;
        if(mol->nAtoms<2)return FALSE;
#ifndef OS_WIN32
	fileNameSH = strdup_printf("%sPCGOne.sh",fileNamePrefix);
#else
	fileNameSH = strdup_printf("%sPCGOne.bat",fileNamePrefix);
#endif
 	fileSH = fopen(fileNameSH, "w");
	if(!fileSH) return FALSE;

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
		if(mol->spinMultiplicity==1)
			fprintf(file," $CONTRL SCFTYP=RHF $END\n");
		else
			fprintf(file," $CONTRL SCFTYP=UHF $END\n");
	}

	fprintf(file," $CONTRL ICHARG=%d MULT=%d $END\n",mol->totalCharge,mol->spinMultiplicity);
	if(strstr(keyWords,"GBASIS"))
	{
		sscanf(strstr(keyWords,"GBASIS"),"%s",buffer);
		fprintf(file," $BASIS %s $END\n",buffer);
	}
	if(strstr(keyWords,"Optimize"))
	{
        	fprintf(file, " $STATPT OptTol=1e-4 NStep=500 $END\n");
	}
	if(strstr(keyWords,"Optimize"))
	{
		putFireFlyMoleculeXYZFixed(mol, file);
	}
	fprintf(file," $DATA\n");
	fprintf(file,"Molecule specification\n");
	fprintf(file,"C1\n");
	for(j=0;j<mol->nAtoms;j++)
	{
		char* symbol = mol->atoms[j].prop.symbol;
		SAtomsProp prop = propAtomGet(symbol);
		fprintf(file,"%s %f %f %f %f\n", 
			symbol,
			(double)prop.atomicNumber,
			mol->atoms[j].coordinates[0],
			mol->atoms[j].coordinates[1],
			mol->atoms[j].coordinates[2]
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
         	fprintf(fileSH,"del PUNCH 2> nul\n");
         	fprintf(fileSH,"del /Q  \"%stmp\"\n",fileNamePrefix);
         	fprintf(fileSH,"rmdir  \"%stmp\"\n",fileNamePrefix);
	}
	else
		fprintf(fileSH,"%s %s",fireflyCommand,fileNameIn);
#endif
	fclose(fileSH);
	{
		char* str = NULL;
		if(strstr(keyWords,"Optimiz")) str = strdup_printf("Minimization by FireFly/%s ... Please wait",keyWords);
		else str = strdup_printf("Computing of energy by FireFly/%s .... Please wait",keyWords);
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
	if(getEnergyFireFly(fileNameOut,&energy))
	{
		printf("Energy by FireFly = %f\n", energy);
		if(strstr(keyWords,"Optimiz")) mol->klass->readGeomFromGamessOutputFile(mol, fileNameOut, -1);
		else mol->klass->readGeomFromGamessOutputFile(mol, fileNameOut, 1);
		if(strstr(keyWords,"Optimiz"))
		{
			char* str = strdup_printf("%s.gab",fileNamePrefix);
			saveGeometry(mol, energy, str);
			printf("----------------------------------------- \n");
			printf("Optimized geometry saved in %s file\n",str);
			printf("----------------------------------------- \n");
			if(str) free(str);
		}
		else
		{
			char* str = strdup_printf("%s.gab",fileNamePrefix);
			saveGeometry(mol, energy, str);
			printf("----------------------------------------- \n");
			printf("Geometry saved in %s file\n",str);
			printf("----------------------------------------- \n");
			if(str) free(str);
		}
	}
	else
	{
		printf(
				(
				"Sorry, I cannot read the output file :  %s"
				" ; Check also the installation of FireFly...")
				,
				fileNameOut
				);
		exit(1);
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
static boolean getEnergyOrca(char* fileNameOut, double* energy)
{
	FILE* file = NULL;
	char buffer[1024];
	char* pdest = NULL;
	char* energyTag = "FINAL SINGLE POINT ENERGY";

	*energy=1e10;
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
			*energy *= AUTOKCAL;
			return TRUE;
		}
	 }
	fclose(file);
	return FALSE;
}
/*****************************************************************************/
static boolean runOneOrca(Molecule* mol, char* fileNamePrefix, char* keyWords, char* NameCommandOrca)
{
	FILE* file = NULL;
	FILE* fileSH = NULL;
	char* fileNameIn = NULL;
	char* fileNameOut = NULL;
	char* fileNameSH = NULL;
	char multiplicityStr[100];
	char buffer[1024];
	double energy = 0;
	int i;
	int nV;

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
	if(!strcmp(NameCommandOrca,"orca") || !strcmp(NameCommandOrca,"nohup orca"))
	{
		fprintf(fileSH,"%s %s > %s\n",NameCommandOrca,fileNameIn,fileNameOut);
		fprintf(fileSH,"exit\n");
	}
	else
		fprintf(fileSH,"%s %s",NameCommandOrca,fileNameIn);
#else
	 if(!strcmp(NameCommandOrca,"orca") )
	{
		if(strstr(orcaDirectory,"\"")) fprintf(fileSH,"set PATH=%s;%cPATH%c\n",orcaDirectory,'%','%');
		else fprintf(fileSH,"set PATH=\"%s\";%cPATH%c\n",orcaDirectory,'%','%');
		fprintf(fileSH,"%s %s > %s\n",NameCommandOrca,fileNameIn,fileNameOut);
		fprintf(fileSH,"exit\n");
	}
	else
		fprintf(fileSH,"%s %s",NameCommandOrca,fileNameIn);
#endif
	fclose(fileSH);
	{
		char* str = NULL;
		if(strstr(keyWords,"Opt")) str = strdup_printf("Minimization by Orca/%s ... Please wait",keyWords);
		else str = strdup_printf("Computing of energy by Orca/%s .... Please wait",keyWords);
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
	if(getEnergyOrca(fileNameOut,&energy))
	{
		printf("Energy by Orca = %f\n", energy);
		mol->klass->readGeomFromOrcaOutputFile(mol, fileNameOut, -1);
		if(strstr(keyWords,"Opt"))
		{
			char* str = strdup_printf("%s.gab",fileNamePrefix);
			saveGeometry(mol, energy, str);
			printf("----------------------------------------- \n");
			printf("Optimized geometry saved in %s file\n",str);
			printf("----------------------------------------- \n");
			if(str) free(str);
		}
		else
		{
			char* str = strdup_printf("%s.gab",fileNamePrefix);
			saveGeometry(mol, energy, str);
			printf("----------------------------------------- \n");
			printf("Geometry saved in %s file\n",str);
			printf("----------------------------------------- \n");
			if(str) free(str);
		}
	}
	else
	{
		printf(
				(
				"Sorry, I cannot read the output file :  %s"
				" ; Check also the installation of Orca..."
				),
				fileNameOut
				);
 		if(fileNameIn) free(fileNameIn);
 		if(fileNameOut) free(fileNameOut);
 		if(fileNameSH) free(fileNameSH);
		exit(1);
		return FALSE;
	}

 	if(fileNameIn) free(fileNameIn);
 	if(fileNameOut) free(fileNameOut);
 	if(fileNameSH) free(fileNameSH);
	return TRUE;
}
/*****************************************************************************/
static boolean getEnergyOpenBabel(char* fileNameOut, double* energy)
{
	FILE* file = NULL;
	char buffer[1024];
	char* pdest = NULL;
	//char* energyTag = "TOTAL ENERGY =";
	char* energyTag = "FINAL ENERGY:";

 	file = fopen(fileNameOut, "r");
	if(!file) return FALSE;
	 while(!feof(file))
	 {
		if(!fgets(buffer,BSIZE,file))break;
		pdest = strstr( buffer, energyTag);
		if(pdest &&sscanf(pdest+strlen(energyTag)+1,"%lf",energy)==1)
		{
			fclose(file);
			if(strstr(pdest,"kJ")) *energy /= KCALTOKJ;
			return TRUE;
		}
	 }
	fclose(file);
        return FALSE;
}
/*****************************************************************************/
static boolean runOneOpenBabel(Molecule* mol, char* fileNamePrefix, char* keyWords, char* NameCommandOpenBabel)
{
	FILE* fileSH = NULL;
	char* fileNameIn = NULL;
	char* fileNameOut = NULL;
	char* fileNameSH = NULL;
	char buffer[1024];
	double energy;
#ifdef OS_WIN32
	char c='%';
#endif

	if(mol->nAtoms<1) return FALSE;
#ifndef OS_WIN32
	fileNameSH = strdup_printf("%sOne.sh",fileNamePrefix);
#else
	fileNameSH = strdup_printf("%sOne.bat",fileNamePrefix);
#endif
 	fileSH = fopen(fileNameSH, "w");
	if(!fileSH) return FALSE;
#ifdef OS_WIN32
	fprintf(fileSH,"@echo off\n");
	fprintf(fileSH,"set PATH=%cPATH%c;\"%s\"\n",c,c,openBabelDirectory);
#endif

	fileNameIn = strdup_printf("%sOne.hin",fileNamePrefix);
	fileNameOut = strdup_printf("%sOne.out",fileNamePrefix);

	if(!mol->klass->saveHIN(mol,fileNameIn))
	{
 		if(fileNameIn) free(fileNameIn);
 		if(fileNameOut) free(fileNameOut);
 		if(fileNameSH) free(fileNameSH);
		return FALSE;
	}
#ifndef OS_WIN32
	if(!strcmp(NameCommandOpenBabel,"obgradient") || !strcmp(NameCommandOpenBabel,"obopt"))
	{
		if(strstr(keyWords,"Opt")) fprintf(fileSH,"%s -ff gaff %s > %s\n","obopt",fileNameIn,fileNameOut);
		else fprintf(fileSH,"%s -ff gaff %s > %s\n","obenergy",fileNameIn,fileNameOut);
		fprintf(fileSH,"exit\n");
	}
	else 
	{
		if(strstr(keyWords,"Opt") && !strstr( NameCommandOpenBabel,"obopt")) 
		{
			char** ssplit = NULL;
			int nA = 0;
			int i;
			ssplit = split(NameCommandOpenBabel);
			while(ssplit && ssplit[nA]!=NULL) nA++;
			fprintf(fileSH,"%s ", "obopt");
			for(i=1;i<nA;i++) fprintf(fileSH,"%s ",  ssplit[i]);
			fprintf(fileSH," %s > %s 2>/dev/null", fileNameIn, fileNameOut);
			strfreev(ssplit);
		}
		else fprintf(fileSH,"%s %s > %s 2>/dev/null", NameCommandOpenBabel, fileNameIn, fileNameOut);
	}
#else
	 if(!strcmp(NameCommandOpenBabel,"obabel") )
	{
		if(strstr(openBabelDirectory,"\"")) fprintf(fileSH,"set PATH=%s;%cPATH%c\n",openBabelDirectory,'%','%');
		else fprintf(fileSH,"set PATH=\"%s\";%cPATH%c\n",openBabelDirectory,'%','%');
		if(strstr(keyWords,"Opt")) fprintf(fileSH,"%s -ff gaff %s > %s\n","obminimize",fileNameIn,fileNameOut);
		else fprintf(fileSH,"%s -ff gaff %s > %s\n","obenergy",fileNameIn,fileNameOut);
		fprintf(fileSH,"exit\n");
	}
	else
	{
		if(strstr(keyWords,"Opt") && !strstr( NameCommandOpenBabel,"obopt")) 
		{
			char** ssplit = NULL;
			int nA = 0;
			int i;
			ssplit = split(NameCommandOpenBabel);
			while(ssplit && ssplit[nA]!=NULL) nA++;
			fprintf(fileSH,"%s ", "obopt");
			for(i=1;i<nA;i++) fprintf(fileSH,"%s ",  ssplit[i]);
			fprintf(fileSH," %s > %s 2>/dev/null", fileNameIn, fileNameOut);
			strfreev(ssplit);
		}
		else fprintf(fileSH,"%s %s > %s", NameCommandOpenBabel, fileNameIn, fileNameOut);
	}
#endif
	fclose(fileSH);
#ifndef OS_WIN32
	/*
	sprintf(buffer,"cat %s",fileNameSH);
	system(buffer);
	sprintf(buffer,"cat %s",fileNameIn);
	system(buffer);
	*/



	sprintf(buffer,"chmod u+x %s",fileNameSH);
	system(buffer);
	system(fileNameSH);
#else
	sprintf(buffer,"\"%s\"",fileNameSH);
	system(buffer);
#endif
	if(getEnergyOpenBabel(fileNameOut,&energy))
	{
		printf("Energy by OpenBabel = %f\n", energy);
		mol->klass->readGeomFromOpenBabelOutputFile(mol, fileNameOut, -1);
		if(strstr(keyWords,"Opt"))
		{
			char* str = strdup_printf("%s.gab",fileNamePrefix);
			saveGeometry(mol, energy, str);
			printf("----------------------------------------- \n");
			printf("Optimized geometry saved in %s file\n",str);
			printf("----------------------------------------- \n");
			if(str) free(str);
		}
		else
		{
			char* str = strdup_printf("%s.gab",fileNamePrefix);
			saveGeometry(mol, energy, str);
			printf("----------------------------------------- \n");
			printf("Geometry saved in %s file\n",str);
			printf("----------------------------------------- \n");
			if(str) free(str);
		}
	}
	else
	{
		printf(
				(
				"Sorry, I cannot read the output file :  %s"
				" ; Check also the installation of OpenBabel..."
				),
				fileNameOut
				);
 		if(fileNameIn) free(fileNameIn);
 		if(fileNameOut) free(fileNameOut);
 		if(fileNameSH) free(fileNameSH);
		exit(1);
		return FALSE;
	}

 	if(fileNameIn) free(fileNameIn);
 	if(fileNameOut) free(fileNameOut);
 	if(fileNameSH) free(fileNameSH);
	return TRUE;
}
/*************************************************************************************************************/
/*
static boolean getEnergyGaussian(char* fileNameOut, double* energy)
{
	FILE* file = NULL;
	char buffer[1024];
	char* pdest = NULL;
	boolean OK = FALSE;

	*energy=1e10;
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
					break;
				}
			}
		}
	 }
	fclose(file);
	*energy *= AUTOKCAL;
	return OK;
}
*/
/**********************************************************************/
static boolean getEnergyGaussian(char* fileNameOut, double* energy)
{
	FILE* file = NULL;
	char buffer[1024];
	char* pdest = NULL;
	char* pos = NULL;
	boolean OK = FALSE;

 	file = fopen(fileNameOut, "r");
	if(!file) return FALSE;
	 while(!feof(file))
	 {
		if(!fgets(buffer,BSIZE,file))break;
		pdest = NULL;
		if(strstr( buffer, "SCF Done:  E(") && strstr( buffer, "=")) pdest = strstr( buffer, "=")+1;
		if(!pdest) if(strstr( buffer, "Energy=") && !strstr( buffer, "hange") ) pdest = strstr( buffer, "Energy=")+strlen("Energy=");
		if(!pdest) if(strstr( buffer, "UMP") && strstr( buffer, "=")) { pos = strstr( buffer, "UMP"); pdest = strstr(pos, "=")+1;};
		if(!pdest) if(strstr( buffer, " E(B2PLYP)") && strstr( buffer, "=")) { pos = strstr( buffer, " E(B2PLYP)"); pdest = strstr(pos, "=")+1;};
		if(!pdest) if(strstr( buffer, "E(CI") && strstr( buffer, "=")) { pos = strstr( buffer, "E(CI"); pdest = strstr(pos, "=")+1;};
		if(!pdest) if(strstr( buffer, " E(CORR)") && strstr( buffer, "=")) { pos = strstr( buffer, " E(CORR)"); pdest = strstr(pos, "=")+1;};
		if(!pdest) if(strstr( buffer, "CCSD(T)=")) { pos = strstr( buffer, "CCSD(T)="); pdest = strstr(pos, "=")+1;};
		if(pdest) 
		{
			int l = strlen(pdest);
			int i;
			for(i=0;i<l;i++) if(pdest[i]=='D' || pdest[i]=='E') pdest[i] ='E';
			if(sscanf(pdest+1,"%lf",energy)==1)
			{
				*energy *= AUTOKCAL;
				OK = TRUE;
				/* break;*/
			}
		}
	 }
	fclose(file);
	return OK;
}
/*************************************************************************************************************/
static boolean putGaussianMoleculeInFile(Molecule* mol, FILE* file)
{
        char buffer[BSIZE];
	int i;
	int k = 0;

	if(mol->nAtoms<1) return FALSE;
      	for (i=0;i<mol->nAtoms;i++)
	{
		k = -1;
		if(mol->atoms[i].variable) k = 0;
		sprintf(buffer,"%s %d %20.14f %20.14f %20.14f\n",mol->atoms[i].prop.symbol, k,
				mol->atoms[i].coordinates[0], 
				mol->atoms[i].coordinates[1],
				mol->atoms[i].coordinates[2]
				);
       		fprintf(file, "%s",buffer);
	}
       	fprintf(file, "\n");
	return TRUE;
}
/*****************************************************************************/
static boolean runOneGaussian(Molecule* mol, char* fileNamePrefix, char* keyWords, char* gaussianCommand)
{
	FILE* file = NULL;
	FILE* fileSH = NULL;
	char* fileNameIn = NULL;
	char* fileNameOut = NULL;
	char* fileNameSH = NULL;
	char multiplicityStr[100];
	char buffer[1024];
	double energy = 0;
#ifdef OS_WIN32
	char c='%';
#endif

	if(!mol) return FALSE;
	if(mol->nAtoms<1) return FALSE;
#ifndef OS_WIN32
	fileNameSH = strdup_printf("%sGaussianOne.sh",fileNamePrefix);
#else
	fileNameSH = strdup_printf("%sGaussianOne.bat",fileNamePrefix);
#endif
 	fileSH = fopen(fileNameSH, "w");
	if(!fileSH) return FALSE;
#ifdef OS_WIN32
	fprintf(fileSH,"@echo off\n");
	fprintf(fileSH,"set PATH=%cPATH%c;\"%s\"\n",c,c,gaussianDirectory);
#endif

	getMultiplicityName(mol->spinMultiplicity, multiplicityStr);

	fileNameIn = strdup_printf("%sOne.com",fileNamePrefix);
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
	fprintf(file,"Gaussian file generated by cchemi\n");
	fprintf(file,"\n");
	fprintf(file,"%d %d\n",mol->totalCharge,mol->spinMultiplicity);
  	putGaussianMoleculeInFile(mol, file);
	fclose(file);
	{
		char* str = NULL;
		if(strstr(keyWords,"Opt")) str = strdup_printf("Minimization by Gaussian/%s ... Please wait",keyWords);
		else str = strdup_printf("Energy by Gaussian/%s ... Please wait",keyWords);
		printf("%s\n",str);
		if(str) free(str);
	}
#ifndef OS_WIN32
	fprintf(fileSH,"%s %s\n",gaussianCommand,fileNameIn);
	fclose(fileSH);
	sprintf(buffer,"chmod u+x %s",fileNameSH);
	system(buffer);
	system(fileNameSH);
#else
	fprintf(fileSH,"\"%s\" \"%s\"\n",gaussianCommand,fileNameIn);
	fclose(fileSH);
	sprintf(buffer,"\"%s\"",fileNameSH);
	system(buffer);
#endif

	fileNameOut = strdup_printf("%sOne.log",fileNamePrefix);
	if(getEnergyGaussian(fileNameOut,&energy))
	{
		char* str = NULL;
		mol->klass->readGeomFromGaussianOutputFile(mol, fileNameOut, -1);
		str = strdup_printf("Energy by Gaussian = %f", energy);
		printf("%s\n",str);
		if(str) free(str);
		if(strstr(keyWords,"Opt"))
		{
			str = strdup_printf("%s.gab",fileNamePrefix);
			saveGeometry(mol, energy, str);
			printf("----------------------------------------- \n");
			printf("Optimized geometry saved in %s file\n",str);
			printf("----------------------------------------- \n");
			if(str) free(str);
		}
		else
		{
			str = strdup_printf("%s.gab",fileNamePrefix);
			saveGeometry(mol, energy, str);
			printf("----------------------------------------- \n");
			printf("Geometry saved in %s file\n",str);
			printf("----------------------------------------- \n");
			if(str) free(str);
		}
	}
	else
	{
		char* str = NULL;
		str = strdup_printf(
				(
				"Sorry, I cannot read the output file : %s "
				" ; Check also the installation of Gaussian..."
				),
				fileNameOut
				);
		printf("%s",str);
		if(str) free(str);
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
static boolean getEnergyGeneric(char* fileNameOut, double* energy)
{
	FILE* file = NULL;
	char buffer[1024];
	int i;
 	file = fopen(fileNameOut, "r");
	if(!file) return FALSE;
	if(!fgets(buffer,BSIZE,file)) { fclose(file); return FALSE;}/* first line for energy in Hartree*/

	for(i=0;i<strlen(buffer);i++) if(buffer[i]=='D' || buffer[i]=='d') buffer[i] ='E';
	if(sscanf(buffer,"%lf",energy)==1)
	{
		fclose(file);
		*energy *=AUTOKCAL;
		return TRUE;
	}
	fclose(file);
	return FALSE;
}
/*****************************************************************************/
static boolean runOneGeneric(Molecule* mol, char* fileNamePrefix, char* keyWords, char* genericCommand)
{
	FILE* file = NULL;
	FILE* fileSH = NULL;
	char* fileNameIn = NULL;
	char* fileNameOut = NULL;
	char* fileNameSH = NULL;
	char multiplicityStr[100];
	char buffer[1024];
	double energy = 0;
	int type = 0;
#ifdef OS_WIN32
	char c='%';
#endif

	if(!mol) return FALSE;
	if(mol->nAtoms<1) return FALSE;
#ifndef OS_WIN32
	fileNameSH = strdup_printf("%sGenericOne.sh",fileNamePrefix);
#else
	fileNameSH = strdup_printf("%sGenericOne.bat",fileNamePrefix);
#endif
 	fileSH = fopen(fileNameSH, "w");
	if(!fileSH) return FALSE;
#ifdef OS_WIN32
	fprintf(fileSH,"@echo off\n");
#endif

	getMultiplicityName(mol->spinMultiplicity, multiplicityStr);

	fileNameIn = strdup_printf("%sOne.inp",fileNamePrefix);
	fileNameOut = strdup_printf("%sOne.out",fileNamePrefix);

 	file = fopen(fileNameIn, "w");
	if(!file) 
	{
 		if(fileNameIn) free(fileNameIn);
 		if(fileNameOut) free(fileNameOut);
 		if(fileNameSH) free(fileNameSH);
		return FALSE;
	}
	if(strstr(keyWords,"Opt")) type = 2;
	if(strstr(keyWords,"ENGRAD")) type = 1;
	fprintf(file,"%d\n",type);
	mol->klass->addMolecule(mol,file);
	fclose(file);
	{
		char* str = NULL;
		if(type==2) str = strdup_printf("Minimization by Generic/%s ... Please wait",genericCommand);
		else str = strdup_printf("Energy by Generic/%s ... Please wait",genericCommand);
		printf("%s\n",str);
		if(str) free(str);
	}
#ifndef OS_WIN32
	fprintf(fileSH,"%s %s %s",genericCommand,fileNameIn,fileNameOut);
	fclose(fileSH);
	sprintf(buffer,"chmod u+x %s",fileNameSH);
	system(buffer);
	system(fileNameSH);
#else
	fprintf(fileSH,"\"%s\" \"%s\" \"%s\"",genericCommand,fileNameIn,fileNameOut);
	fclose(fileSH);
	sprintf(buffer,"\"%s\"",fileNameSH);
	system(buffer);
#endif

	if(getEnergyGeneric(fileNameOut,&energy))
	{
		char* str = NULL;
		mol->klass->readGeometry(mol,fileNameOut);
		str = strdup_printf("Energy by %s = %f", genericCommand,energy);
		printf("%s\n",str);
		if(str) free(str);
		if(strstr(keyWords,"Opt"))
		{
			str = strdup_printf("%s.gab",fileNamePrefix);
			saveGeometry(mol, energy, str);
			printf("----------------------------------------- \n");
			printf("Optimized geometry saved in %s file\n",str);
			printf("----------------------------------------- \n");
			if(str) free(str);
		}
		else
		{
			str = strdup_printf("%s.gab",fileNamePrefix);
			saveGeometry(mol, energy, str);
			printf("----------------------------------------- \n");
			printf("Geometry saved in %s file\n",str);
			printf("----------------------------------------- \n");
			if(str) free(str);
		}
	}
	else
	{
		char* str = NULL;
		str = strdup_printf(
				(
				"Sorry, I cannot read the output file : %s "
				" ; Check also the installation of Generic..."
				),
				fileNameOut
				);
		printf("%s",str);
		if(str) free(str);
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
static void runMopac(Molecule* mol, char* fileName, char* keys, char* mopacCommand)
{
	if(mopacCommand&& keys)
	{
		char* fileNamePrefix = getSuffixNameFile(fileName);
		printf("keys = %s\n",keys);
		if(runOneMopac(mol, fileNamePrefix, keys,mopacCommand))
		{
		}
		if(fileNamePrefix) free(fileNamePrefix);
	}
}
/*****************************************************************************/
static void runOrca(Molecule* mol, char* fileName, char* keys, char* NameCommandOrca)
{
	if(NameCommandOrca && keys)
	{
		char* fileNamePrefix = getSuffixNameFile(fileName);
		if(runOneOrca(mol, fileNamePrefix, keys,NameCommandOrca))
		{
			char* fileOut = strdup_printf("%sOne.out",fileNamePrefix);
			mol->klass->readGeomFromOrcaOutputFile(mol, fileOut, -1);
			if(fileOut) free(fileOut);
		}
		if(fileNamePrefix) free(fileNamePrefix);
	}
}
/*****************************************************************************/
static void runOpenBabel(Molecule* mol, char* fileName, char* keys, char* NameCommandOpenBabel)
{
	if(NameCommandOpenBabel && keys)
	{
		char* fileNamePrefix = getSuffixNameFile(fileName);
		if(runOneOpenBabel(mol, fileNamePrefix, keys,NameCommandOpenBabel))
		{
			char* fileOut = strdup_printf("%sOne.out",fileNamePrefix);
			mol->klass->readGeomFromOpenBabelOutputFile(mol, fileOut, -1);
			if(fileOut) free(fileOut);
		}
		if(fileNamePrefix) free(fileNamePrefix);
	}
}
/*****************************************************************************/
static void runFireFly(Molecule* mol, char* fileName, char* keys, char* fireflyCommand)
{
	if(keys && fireflyCommand)
	{
		char* fileNamePrefix = getSuffixNameFile(fileName);
		if(runOneFireFly(mol, fileNamePrefix, keys,fireflyCommand))
		{
		}
		if(fileNamePrefix) free(fileNamePrefix);
	}
}
/*****************************************************************************/
static void runGaussian(Molecule* mol, char* fileName, char* keys, char* gaussianCommand)
{
	if(keys && gaussianCommand)
	{
		char* fileNamePrefix = getSuffixNameFile(fileName);
		if(runOneGaussian(mol, fileNamePrefix, keys,gaussianCommand))
		{
		}
		if(fileNamePrefix) free(fileNamePrefix);
	}
}
/*****************************************************************************/
static void runGeneric(Molecule* mol, char* fileName, char* keys, char* genericCommand)
{
	if(keys && genericCommand)
	{
		char* fileNamePrefix = getSuffixNameFile(fileName);
		if(runOneGeneric(mol, fileNamePrefix, keys, genericCommand))
		{
		}
		if(fileNamePrefix) free(fileNamePrefix);
	}
}
/*****************************************************************************/
void quantumMechanicsDlg(char* inputFileName)
{
	Molecule mol = *(readMolecule(inputFileName, TRUE));
	char* fileName = NULL;
	char* mopacCommand = strdup("/opt/mopac/MOPAC2012");
	char* fireflyCommand = strdup("firefly");
	char* orcaCommand = strdup("orca");
	char* openBabelCommand = strdup("obgradient");
	char* gaussianCommand = strdup("g09");
	char* genericCommand = strdup("runGeneric");
	char* N2P2Dir=strdup(getenv("PWD"));
	char* tmModule=strdup_printf("%s/%s",getenv("PWD"),"CChemITMModule");
	FILE* file = NULL;
	char* runType = NULL;
	char* model = NULL;
	char* QMKeys = NULL;
	if(!inputFileName)
	{
		printf("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n");
		printf("Usage : cchemi inputFileName.inp\n");
		printf("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n");
		exit(1);
	}
	file = fopen(inputFileName, "r");
	if(!file)
	{
		printf("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n");
		printf("Sorry, I cannot open the input file : %s\n",inputFileName);
		printf("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n");
		exit(1);
	}
	if(!readOneString(file,"RunType",&runType)) runType = strdup("ENERGY");
	if(!readOneString(file,"Model",&model)) model = strdup("MOPAC");
	uppercase(model);
	if(!readOneString(file,"QMKeys",&QMKeys)) 
	{
		if(!strcmp(model,"MOPAC")) QMKeys = strdup("PM6-DH+");
		else if(!strcmp(model,"ORCA")) QMKeys = strdup("AM1");
		else if(!strcmp(model,"GAUSSIAN")) QMKeys = strdup("AM1");
		else if(!strcmp(model,"FIREFLY")) QMKeys = strdup("AM1");
		else if(!strcmp(model,"OPENBABEL")) QMKeys = strdup("mmff94");
		else QMKeys = strdup("SCC-DFTB");
	}
	if(!strcmp(model,"OPENBABEL")) mol.klass->buildMMTypes(&mol, file);

	uppercase(runType);
	readOneString(file,"mopacCommand",&mopacCommand);
	readOneString(file,"orcaCommand",&orcaCommand);
	readOneString(file,"openBabelCommand",&openBabelCommand);
	readOneString(file,"HDNNDir",&N2P2Dir);
	readOneString(file,"N2P2Dir",&N2P2Dir);
	readOneString(file,"tmModule",&tmModule);
	readOneString(file,"fireflyCommand",&fireflyCommand);
	readOneString(file,"gaussianCommand",&gaussianCommand);
	readOneString(file,"genericCommand",&genericCommand);
	if(!strcmp(model,"FIREFLY")) QMKeys = strdup("GBASIS=AM1");
	printf("QMKeys = %s\n",QMKeys);
	printf("Model = %s\n",model);
	printf("openBabelCommand = %s\n",openBabelCommand);
	if(!strcmp(model,"MOPAC") && (!strcmp(runType,"ENERGY")|| strstr(QMKeys,"ESP") ||  strstr(QMKeys,"POINT")))
	{
		char* keys = strdup_printf("1SCF %s",QMKeys);
		printf("keys = %s\n",keys);
		printf("mopacCommand = %s\n",mopacCommand);
		fileName = strdup_printf("cchemiMopac.inp");
		runMopac(&mol, fileName, keys, mopacCommand);
		free(keys);
	}
	else if(!strcmp(model,"MOPAC") && strstr(runType,"GRADIENT"))
	{
		char* keys = strdup_printf("1SCF GRAD %s",QMKeys);
		printf("keys = %s\n",keys);
		printf("mopacCommand = %s\n",mopacCommand);
		fileName = strdup_printf("cchemiMopac.inp");
		runMopac(&mol, fileName, keys, mopacCommand);
		free(keys);
	}
	else if(!strcmp(model,"MOPAC") && strstr(runType,"OPTIMIZATION"))
	{
		char* keys = strdup_printf("XYZ %s",QMKeys);
		printf("keys = %s\n",keys);
		printf("mopacCommand = %s\n",mopacCommand);
		fileName = strdup_printf("cchemiMopac.inp");
		runMopac(&mol, fileName, keys, mopacCommand);
		free(keys);
	}
	else if(!strcmp(model,"ORCA") && strstr(runType,"ENERGY"))
	{
		char* keys = strdup_printf("%s",QMKeys);
		fileName = strdup_printf("cchemiOrca.inp");
		runOrca(&mol, fileName, keys, orcaCommand);
		free(keys);
	}
	else if(!strcmp(model,"ORCA") && strstr(runType,"GRADIENT"))
	{
		char* keys = strdup_printf("%s ENGRAD",QMKeys);
		fileName = strdup_printf("cchemiOrca.inp");
		runOrca(&mol, fileName, keys, orcaCommand);
		free(keys);
	}
	else if(!strcmp(model,"ORCA") && strstr(runType,"OPTIMIZATION"))
	{
		char* keys = strdup_printf("%s Opt",QMKeys);
		fileName = strdup_printf("cchemiOrca.inp");
		runOrca(&mol, fileName, keys, orcaCommand);
		free(keys);
	}
	else if(!strcmp(model,"FIREFLY") && strstr(runType,"ENERGY"))
	{
		char* keys = strdup_printf("RUNTYP=Energy %s",QMKeys);
		fileName = strdup_printf("cchemiFF.inp");
		runFireFly(&mol, fileName, keys, fireflyCommand);
		free(keys);
	}
	else if(!strcmp(model,"FIREFLY") && strstr(runType,"GRADIENT"))
	{
		char* keys = strdup_printf("RUNTYP=GRADIENT %s",QMKeys);
		printf("keys=%s\n",QMKeys);
		fileName = strdup_printf("cchemiFF.inp");
		runFireFly(&mol, fileName, keys, fireflyCommand);
		free(keys);
	}
	else if(!strcmp(model,"FIREFLY") && strstr(runType,"OPTIMIZATION"))
	{
		char* keys = strdup_printf("RUNTYP=Optimize %s",QMKeys);
		printf("keys=%s\n",QMKeys);
		fileName = strdup_printf("cchemiFF.inp");
		runFireFly(&mol, fileName, keys, fireflyCommand);
		free(keys);
	}
	else if(!strcmp(model,"GAUSSIAN") && strstr(runType,"ENERGY"))
	{
		char* keys = strdup_printf("%s",QMKeys);
		fileName = strdup_printf("cchemiGauss.com");
		runGaussian(&mol, fileName, keys, gaussianCommand);
		free(keys);
	}
	else if(!strcmp(model,"GAUSSIAN") && strstr(runType,"GRADIENT"))
	{
		char* keys = strdup_printf("%s Force ",QMKeys);
		printf("keys=%s\n",QMKeys);
		fileName = strdup_printf("cchemiGauss.com");
		runGaussian(&mol, fileName, keys, gaussianCommand);
		free(keys);
	}
	else if(!strcmp(model,"GAUSSIAN") && strstr(runType,"OPTIMIZATION"))
	{
		char* keys = strdup_printf("%s Opt",QMKeys);
		printf("keys=%s\n",QMKeys);
		fileName = strdup_printf("cchemiGauss.com");
		runGaussian(&mol, fileName, keys, gaussianCommand);
		free(keys);
	}
	else if(!strcmp(model,"GENERIC") && strstr(runType,"ENERGY"))
	{
		char* keys = strdup_printf("%s",QMKeys);
		fileName = strdup_printf("cchemiGen.com");
		runGeneric(&mol, fileName, keys, genericCommand);
		free(keys);
	}
	else if(!strcmp(model,"GENERIC") && strstr(runType,"GRADIENT"))
	{
		char* keys = strdup_printf("%s ENGRAD",QMKeys);
		fileName = strdup_printf("cchemiGen.com");
		runGeneric(&mol, fileName, keys, genericCommand);
		free(keys);
	}
	else if(!strcmp(model,"GENERIC") && strstr(runType,"OPTIMIZATION"))
	{
		char* keys = strdup_printf("%s Opt",QMKeys);
		fileName = strdup_printf("cchemiGauss.com");
		runGeneric(&mol, fileName, keys, genericCommand);
		free(keys);
	}
	else if(!strcmp(model,"OPENBABEL") && strstr(runType,"ENERGY"))
	{
		char* keys = strdup_printf("%s",QMKeys);
		fileName = strdup_printf("cchemiOB.hin");
		runOpenBabel(&mol, fileName, keys, openBabelCommand);
		free(keys);
	}
	else if(!strcmp(model,"OPENBABEL") && strstr(runType,"GRADIENT"))
	{
		char* keys = strdup_printf("%s ENGRAD",QMKeys);
		fileName = strdup_printf("cchemiOB.hin");
		runOpenBabel(&mol, fileName, keys, openBabelCommand);
		free(keys);
	}
	else if(!strcmp(model,"OPENBABEL") && strstr(runType,"OPTIMIZATION"))
	{
		char* keys = strdup_printf("%s Opt",QMKeys);
		fileName = strdup_printf("cchemiOB.hin");
		runOpenBabel(&mol, fileName, keys, openBabelCommand);
		free(keys);
	}
	else if(!strcmp(model,"N2P2") || !strcmp(model,"HDNN"))
	{
		
		double cfenergy = 1.0/(AUTOKCAL);
        	double cflength = ANGTOBOHR;
		// compute energy using N2P2 potential with data in N2P2Dir
		runN2P2(&mol, N2P2Dir, cflength, cfenergy, 0, 0);
	}
	else if(!strcmp(model,"TENSORMOL") || !strcmp(model,"TM"))
	{
		
		runTM(&mol, tmModule,TRUE);
	}
	printf("fileName = %s\n",fileName);
	fclose(file);
	if(fileName) free(fileName);
}
/*****************************************************************************/
static boolean saveGeometry(Molecule* molecule, double energy, char* fileNameGeom)
{
	boolean Ok = FALSE;
	double oldEnergy = molecule->potentialEnergy;
	molecule->potentialEnergy = energy;
	Ok = molecule->klass->save(molecule,fileNameGeom);
	molecule->potentialEnergy = oldEnergy;
	return Ok;
}
/*****************************************************************************/
static boolean saveConfoGeometries(int numberOfGeometries, QuantumMechanicsModel** geometries, double* energies, char* fileNameGeom)
{
	FILE* file = NULL;
	int i;
	int j;
	int nG = 0;
	int k;
	int form = 1;

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
	fprintf(file,"deltaE eV 1\n");
	k = -1;
	for(i=0;i<numberOfGeometries;i++)
	{
		if(!geometries[i]) continue;
		if(k<0) k = i;
		fprintf(file,"%f\n",energies[i]);
		//if(k>=0) fprintf(file,"%f\n",(energies[i]-energies[k])*503.21892494); // in K
		if(k>=0) fprintf(file,"%f\n",(energies[i]-energies[k])/AUTOKCAL*AUTOEV); // in eV
		else fprintf(file,"0\n");

		Molecule* mol=&geometries[i]->molecule;
		double oldEnergy = mol->potentialEnergy;
		mol->potentialEnergy = energies[i];
		geometries[i]->molecule.klass->addGeometryToGabedit(&geometries[i]->molecule, file);
		mol->potentialEnergy = oldEnergy;
	}
	fclose(file);
	return TRUE;

}
/*****************************************************************************/
static boolean runOneOptMopac(QuantumMechanicsModel* geom, double* energy, char* fileNamePrefix, char* keyWords, char* mopacCommand)
{
	FILE* file = NULL;
	FILE* fileSH = NULL;
	int j;
	char* fileNameIn = NULL;
	char* fileNameOut = NULL;
	char* fileNameSH = NULL;
	char multiplicityStr[100];
	char buffer[1024];
	*energy = 0;
	Molecule* mol = &geom->molecule;
#ifdef OS_WIN32
	char c='%';
#endif

	if(!geom) return FALSE;
	if(geom->molecule.nAtoms<1) return FALSE;
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

	getMultiplicityName(mol->spinMultiplicity, multiplicityStr);

	fileNameIn = strdup_printf("%sOne.mop",fileNamePrefix);
 	file = fopen(fileNameIn, "w");
	if(!file) 
	{
 		if(fileNameIn) free(fileNameIn);
 		if(fileNameOut) free(fileNameOut);
 		if(fileNameSH) free(fileNameSH);
		return FALSE;
	}
	if(mol->spinMultiplicity>1)
	fprintf(file,"%s UHF CHARGE=%d %s\n",keyWords,mol->totalCharge,multiplicityStr);
	else
	fprintf(file,"%s CHARGE=%d %s\n",keyWords,mol->totalCharge,multiplicityStr);
	fprintf(file,"\n");
	fprintf(file,"Mopac file generated by Gabedit\n");

	for(j=0;j<geom->molecule.nAtoms;j++)
	{
	fprintf(file," %s %f %d %f %d %f %d\n", 
			geom->molecule.atoms[j].prop.symbol,
			geom->molecule.atoms[j].coordinates[0],
			1,
			geom->molecule.atoms[j].coordinates[1],
			1,
			geom->molecule.atoms[j].coordinates[2],
			1
			);
	}
	fclose(file);
	{
		char* str = NULL;
		if(!strstr(keyWords,"1SCF") && !strstr(keyWords,"ESP")) str = strdup_printf("Minimization by Mopac/%s ... Please wait",keyWords);
		else if(strstr(keyWords,"ESP")) str = strdup_printf("ESP charges from Mopac/%s ... Please wait",keyWords);
		else str = strdup_printf("Computing of energy by Mopac/%s .... Please wait",keyWords);
		printf("%s\n",str);
		if(str) free(str);
	}
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
		printf("Energy by Mopac = %f\n", *energy);
		mol->klass->readGeomFromMopacOutputFile(mol, fileNameOut, -1);
	}
	else
	{
		printf("I cannot read energy = from %s file\n",fileNameOut);
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
static boolean runMopacFiles(int numberOfGeometries, QuantumMechanicsModel** geometries, double* energies, char* fileNamePrefix, char* keyWords, char* mopacCommand)
{
	int i;
	int nG = 0;
	int nM = 0;
	FILE* logfile = stdout;
	for(i=0;i<numberOfGeometries;i++)
	{
		if(!geometries[i]) continue;
		nG++;
		fprintf(logfile,"Minimization by Mopac of geometry n = %d... Please wait\n", i+1);
		if(runOneOptMopac(geometries[i], &energies[i], fileNamePrefix, keyWords, mopacCommand)) 
		{
			nM++;
		}
		else
		{
			geometries[i]->klass->free(geometries[i]);
			geometries[i] =NULL;
		}
		fflush(logfile);

	}
	/*
	if(nM==nG) return TRUE;
	return FALSE;
	*/
	fprintf(logfile,"Number of Mopac runs with errors = %d\n", nG-nM); fflush(logfile);
	fprintf(logfile,"-------------------------------------------\n"); fflush(logfile);
	return (nM>0);


}
/*****************************************************************************/
static boolean runOneOptGeneric(QuantumMechanicsModel* geom, double* energy, char* fileNamePrefix, char* keyWords, char* genericCommand)
{
	FILE* file = NULL;
	FILE* fileSH = NULL;
	char* fileNameIn = NULL;
	char* fileNameOut = NULL;
	char* fileNameSH = NULL;
	char buffer[1024];
	int type = 0;
	Molecule* mol = &geom->molecule;
	*energy = 0;
#ifdef OS_WIN32
	char c='%';
#endif

	if(!geom) return FALSE;
	if(geom->molecule.nAtoms<1) return FALSE;
#ifndef OS_WIN32
	fileNameSH = strdup_printf("%sGeneOne.sh",fileNamePrefix);
#else
	fileNameSH = strdup_printf("%sGeneOne.bat",fileNamePrefix);
#endif
 	fileSH = fopen(fileNameSH, "w");
	if(!fileSH) return FALSE;

	fileNameIn = strdup_printf("%sOne.inp",fileNamePrefix);
	fileNameOut = strdup_printf("%sOne.out",fileNamePrefix);

 	file = fopen(fileNameIn, "w");
	if(!file) 
	{
 		if(fileNameIn) free(fileNameIn);
 		if(fileNameOut) free(fileNameOut);
 		if(fileNameSH) free(fileNameSH);
		return FALSE;
	}

	if(strstr(keyWords,"Opt")) type = 2;
	if(strstr(keyWords,"ENGRAD")) type = 1;
	fprintf(file,"%d\n",type);
	mol->klass->addMolecule(mol,file);
	fclose(file);
	{
		char* str = NULL;
		if(strstr(keyWords,"OPT")) str = strdup_printf("Minimization by Generic/%s ... Please wait",genericCommand);
		else str = strdup_printf("Computing of energy by Generic/%s .... Please wait ",genericCommand);
		printf("%s",str);
	}
#ifndef OS_WIN32
	fprintf(fileSH,"%s %s %s",genericCommand,fileNameIn,fileNameOut);
	fclose(fileSH);
	sprintf(buffer,"chmod u+x %s",fileNameSH);
	system(buffer);
	system(fileNameSH);
#else
	fprintf(fileSH,"\"%s\" \"%s\" \"%s\"",genericCommand,fileNameIn,fileNameOut);
	fclose(fileSH);
	sprintf(buffer,"\"%s\"",fileNameSH);
	system(buffer);
#endif
	if(getEnergyGeneric(fileNameOut,energy))
	{
		printf("Energy by Generic = %f\n", *energy);
		mol->klass->readGeometry(mol,fileNameOut);
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
static boolean runGenericFiles(int numberOfGeometries, QuantumMechanicsModel** geometries, double* energies, char* fileNamePrefix, char* keyWords, char* genericCommand)
{
	int i;
	int nG = 0;
	int nM = 0;
	FILE* logfile = stdout;
	for(i=0;i<numberOfGeometries;i++)
	{
		if(!geometries[i]) continue;
		nG++;
		fprintf(logfile,"Minimization by Generic of geometry n = %d... Please wait\n", i+1);
		if(runOneOptGeneric(geometries[i], &energies[i], fileNamePrefix, keyWords, genericCommand)) 
		{
			nM++;
		}
		else
		{
			geometries[i]->klass->free(geometries[i]);
			geometries[i] =NULL;
		}
		fflush(logfile);

	}
	/*
	if(nM==nG) return TRUE;
	return FALSE;
	*/
	fprintf(logfile,"Number of generic runs with errors = %d\n", nG-nM); fflush(logfile);
	fprintf(logfile,"-------------------------------------------\n"); fflush(logfile);
	return (nM>0);

}
/*****************************************************************************/
static boolean runOneOptFireFly(QuantumMechanicsModel* geom, double* energy, char* fileNamePrefix, char* keyWords, char* fireflyCommand)
{
	FILE* file = NULL;
	FILE* fileSH = NULL;
	int j;
	char* fileNameIn = NULL;
	char* fileNameOut = NULL;
	char* fileNameSH = NULL;
	char multiplicityStr[100];
	char buffer[1024];
	Molecule* mol = &geom->molecule;
	*energy = 0;
#ifdef OS_WIN32
	char c='%';
#endif

	if(!geom) return FALSE;
	if(geom->molecule.nAtoms<1) return FALSE;
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
		if(mol->spinMultiplicity==1)
			fprintf(file," $CONTRL SCFTYP=RHF $END\n");
		else
			fprintf(file," $CONTRL SCFTYP=UHF $END\n");
	}

	fprintf(file," $CONTRL ICHARG=%d MULT=%d $END\n",mol->totalCharge,mol->spinMultiplicity);
	if(strstr(keyWords,"GBASIS"))
	{
		sscanf(strstr(keyWords,"GBASIS"),"%s",buffer);
		fprintf(file," $BASIS %s $END\n",buffer);
	}
	fprintf(file," $DATA\n");
	fprintf(file,"Molecule specification\n");
	fprintf(file,"C1\n");
	for(j=0;j<geom->molecule.nAtoms;j++)
	{
		char* symbol = geom->molecule.atoms[j].prop.symbol;
		SAtomsProp prop = propAtomGet(symbol);
		fprintf(file,"%s %f %f %f %f\n", 
			symbol,
			(double)prop.atomicNumber,
			geom->molecule.atoms[j].coordinates[0],
			geom->molecule.atoms[j].coordinates[1],
			geom->molecule.atoms[j].coordinates[2]
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
	{
		char* str = NULL;
		if(strstr(keyWords,"OPTIMIZE")) str = strdup_printf("Minimization by FireFly/%s ... Please wait",keyWords);
		else str = strdup_printf("Computing of energy by FireFly/%s .... Please wait",keyWords);
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
	if(getEnergyFireFly(fileNameOut,energy))
	{
		printf("Energy by FireFly = %f\n", *energy);
		mol->klass->readGeomFromGamessOutputFile(mol, fileNameOut, -1);
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
static boolean runFireFlyFiles(int numberOfGeometries, QuantumMechanicsModel** geometries, double* energies, char* fileNamePrefix, char* keyWords, char* fireflyCommand)
{
	int i;
	int nG = 0;
	int nM = 0;
	FILE* logfile = stdout;
	for(i=0;i<numberOfGeometries;i++)
	{
		if(!geometries[i]) continue;
		nG++;
		fprintf(logfile,"Minimization by FireFly of geometry n = %d... Please wait\n", i+1);
		if(runOneOptFireFly(geometries[i], &energies[i], fileNamePrefix, keyWords,fireflyCommand)) 
		{
			nM++;
		}
		else
		{
			geometries[i]->klass->free(geometries[i]);
			geometries[i] =NULL;
		}
		fflush(logfile);

	}
	/*
	if(nM==nG) return TRUE;
	return FALSE;
	*/
	fprintf(logfile,"Number of FireFly generic runs with errors = %d\n", nG-nM); fflush(logfile);
	fprintf(logfile,"-------------------------------------------\n"); fflush(logfile);
	return (nM>0);


}
/*****************************************************************************/
static boolean runOneOptGaussian(QuantumMechanicsModel* geom, double* energy, char* fileNamePrefix, char* keyWords, char* gaussianCommand)
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
		if(strstr(keyWords,"OPT")) str = strdup_printf("Minimization by Gaussian/%s ... Please wait",keyWords);
		else str = strdup_printf("Computing of energy by Gaussian/%s .... Please wait",keyWords);
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
static boolean runGaussianFiles(int numberOfGeometries, QuantumMechanicsModel** geometries, double* energies, char* fileNamePrefix, char* keyWords, char* gaussianCommand)
{
	int i;
	int nG = 0;
	int nM = 0;
	FILE* logfile = stdout;
	for(i=0;i<numberOfGeometries;i++)
	{
		if(!geometries[i]) continue;
		nG++;
		fprintf(logfile,"Minimization by Gaussian of geometry n = %d... Please wait\n", i+1);
		if(runOneOptGaussian(geometries[i], &energies[i], fileNamePrefix, keyWords,gaussianCommand)) 
		{
			nM++;
		}
		else
		{
			geometries[i]->klass->free(geometries[i]);
			geometries[i] =NULL;
		}
		fflush(logfile);

	}
	/*
	if(nM==nG) return TRUE;
	return FALSE;
	*/
	fprintf(logfile,"Number of Gaussian runs with errors = %d\n", nG-nM); fflush(logfile);
	fprintf(logfile,"-------------------------------------------\n"); fflush(logfile);
	return (nM>0);


}
/*****************************************************************************/
static boolean runOneOptOrca(QuantumMechanicsModel* geom, double* energy, char* fileNamePrefix, char* keyWords, char* orcaCommand)
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
		if(strstr(keyWords,"Opt")) str = strdup_printf("Minimization by Orca/%s ... Please wait",keyWords);
		else str = strdup_printf("Computing of energy by Orca/%s .... Please wait",keyWords);
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
static boolean runOrcaFiles(int numberOfGeometries, QuantumMechanicsModel** geometries, double* energies, char* fileNamePrefix, char* keyWords, char* orcaCommand)
{
	int i;
	int nG = 0;
	int nM = 0;
	FILE* logfile = stdout;
	for(i=0;i<numberOfGeometries;i++)
	{
		if(!geometries[i]) continue;
		nG++;
		fprintf(logfile,"Minimization by Orca of geometry n = %d... Please wait\n", i+1);
		if(runOneOptOrca(geometries[i], &energies[i], fileNamePrefix, keyWords,orcaCommand)) 
		{
			nM++;
		}
		else
		{
			geometries[i]->klass->free(geometries[i]);
			geometries[i] =NULL;
		}
		fflush(logfile);

	}
	/*
	if(nM==nG) return TRUE;
	return FALSE;
	*/
	fprintf(logfile,"Number of Orca runs with errors = %d\n", nG-nM); fflush(logfile);
	fprintf(logfile,"-------------------------------------------\n"); fflush(logfile);
	return (nM>0);


}
/*****************************************************************************/
static boolean runOpenBabelFiles(int numberOfGeometries, QuantumMechanicsModel** geometries, double* energies, char* fileNamePrefix, char* keyWords, char* openBabelCommand)
{
	int i;
	int nG = 0;
	int nM = 0;
	FILE* logfile = stdout;
	if(!geometries) return FALSE;
	for(i=0;i<numberOfGeometries;i++)
	{
		if(!geometries[i]) continue;
		nG++;
		fprintf(logfile,"Minimization by OpenBabel of geometry n = %d... Please wait\n", i+1);
		if(runOneOpenBabel(&geometries[i]->molecule, fileNamePrefix, "Opt",openBabelCommand)) 
		{
			energies[i] = geometries[i]->molecule.potentialEnergy;
			nM++;
		}
		else
		{
			geometries[i]->klass->free(geometries[i]);
			geometries[i] =NULL;
		}
		fflush(logfile);

	}
	/*
	if(nM==nG) return TRUE;
	return FALSE;
	*/
	fprintf(logfile,"Number of OpenBabel runs with errors = %d\n", nG-nM); fflush(logfile);
	fprintf(logfile,"-------------------------------------------\n"); fflush(logfile);
	return (nM>0);


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
static double* getDistancesBetweenAtoms(QuantumMechanicsModel* qmModel)
{
	double* distances = NULL;
	int i;
	int j;
	int n;
	int k;
	if(qmModel->molecule.nAtoms<1) return distances;
	n = qmModel->molecule.nAtoms*(qmModel->molecule.nAtoms-1)/2;
	distances = malloc(n*sizeof(double));
	n = 0;
	for (  i = 0; i < qmModel->molecule.nAtoms-1; i++ )
	for (  j = i+1; j < qmModel->molecule.nAtoms; j++ )
	{
		double x = qmModel->molecule.atoms[i].coordinates[0]-qmModel->molecule.atoms[j].coordinates[0];
		double y = qmModel->molecule.atoms[i].coordinates[1]-qmModel->molecule.atoms[j].coordinates[1];
		double z = qmModel->molecule.atoms[i].coordinates[2]-qmModel->molecule.atoms[j].coordinates[2];
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
static void removedsToEnd(int numberOfGeometries, QuantumMechanicsModel** geometries, double* energies, boolean* removeds)
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
				QuantumMechanicsModel* g = geometries[i];

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
static void computeRemoveds(int numberOfGeometries, QuantumMechanicsModel** geometries, double* energies, boolean *removeds, 
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
	if(!geometries[i]) removeds[i] = TRUE;
	for(i=0;i<numberOfGeometries-1;i++)
	{
		int n;
		if(!geometries[i]) removeds[i] = TRUE;
		if(removeds[i]) continue;
		if(tolDistance>0) distancesI =  getDistancesBetweenAtoms(geometries[i]);
		n = geometries[i]->molecule.nAtoms*(geometries[i]->molecule.nAtoms-1)/2;
		for(j=i+1;j<numberOfGeometries;j++)
		{
			if(!geometries[j]) removeds[j] = TRUE;
			if(removeds[j]) continue;
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
static void removeIdenticalGeometries(int* nG, QuantumMechanicsModel*** geoms, double** eners, double tolEnergy, double tolDistance)
{
	int i;
	int numberOfGeometries =*nG;
	QuantumMechanicsModel** geometries = *geoms; 
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
			if(geometries[i]) geometries[i]->klass->free(geometries[i]);
		}
		else newN++;
	}
	free(removeds);
	if(newN==0) { *nG = newN; return;}
	if(newN==numberOfGeometries) return;
	*nG = newN;
	*eners = realloc(*eners,newN*sizeof(double));
	*geoms = realloc(*geoms,newN*sizeof(QuantumMechanicsModel**));

}
/*****************************************************************************/
/*
static int removeIdenticalGeometriesNULL(int numberOfGeometries, QuantumMechanicsModel** geometries, double* energies, double tolEnergy, double tolDistance)
{
	int i;
	boolean* removeds = NULL;
	int newN = 0;
	if(numberOfGeometries<1) return numberOfGeometries;
	removeds = malloc(numberOfGeometries*sizeof(boolean));
	for(i=0;i<numberOfGeometries;i++) removeds[i] = FALSE;
	computeRemoveds(numberOfGeometries, geometries, energies, removeds, tolEnergy, tolDistance);
	removedsToEnd(numberOfGeometries, geometries, energies, removeds);
	for(i=0;i<numberOfGeometries;i++) 
	{
		if(removeds[i]) 
		{
			if(geometries[i]) geometries[i]->klass->free(geometries[i]);
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
static void sortGeometries(int numberOfGeometries, QuantumMechanicsModel** geometries, double* energies)
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
				QuantumMechanicsModel* g = geometries[i];

				energies[i] = energies[k];
				energies[k] = energy;
				geometries[i] = geometries[k];
				geometries[k] = g;
			}
		}
	}
}
/*****************************************************************************/
static boolean createCChemIFiles(int numberOfGeometries, QuantumMechanicsModel** geometries, double* energies, char* fileNamePrefix, char* keyWords,char* cchemiCommand)
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
		fprintf(file,"OptimizerType=External\n");
		fprintf(file,"Model=Mopac\n");
		fprintf(file,"QMKeys=AM1\n");
		fprintf(file,"mopacCommand=/home/allouche/Softwares/MOPAC2012/MOPAC2012.exe\n");
		fprintf(file,"orcaCommand=orca\n");
		fprintf(file,"fireflyCommand=firefly\n");
		fprintf(file,"gaussianCommand=g09\n");
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
static boolean createMopacFiles(int numberOfGeometries, QuantumMechanicsModel** geometries, double* energies, char* fileNamePrefix, char* keyWords, char* mopacCommand)
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
		if(!geometries[i]) continue;
 		if(fileName) free(fileName);
		getMultiplicityName(geometries[i]->molecule.spinMultiplicity, multiplicityStr);
		fileName = strdup_printf("%s_%d.mop",fileNamePrefix,i+1);
 		file = fopen(fileName, "w");
		if(!file) return FALSE;
		if(geometries[i]->molecule.spinMultiplicity>1)
		fprintf(file,"%s UHF CHARGE=%d %s\n",keyWords,geometries[i]->molecule.totalCharge,multiplicityStr);
		else
		fprintf(file,"%s CHARGE=%d %s\n",keyWords,geometries[i]->molecule.totalCharge,multiplicityStr);
		fprintf(file,"\n");
		fprintf(file," Quantum Mechanics Energy(kCal/mol) =%f\n",energies[i]);

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
static boolean createGaussianFiles(int numberOfGeometries, QuantumMechanicsModel** geometries, double* energies, char* fileNamePrefix, char* keyWords, char* gaussianCommand)
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
		fprintf(file,"File generated by Gabedit\n");
		fprintf(file,"Quantum Mechanics Energy(kCal/mol) = %f\n",energies[i]);
		fprintf(file,"\n");
		fprintf(file,"%d %d\n",geometries[i]->molecule.totalCharge,geometries[i]->molecule.spinMultiplicity);
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
static boolean createFireFlyFiles(int numberOfGeometries, QuantumMechanicsModel** geometries, double* energies, char* fileNamePrefix, char* keyWords, char* fireflyCommand)
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

		fprintf(file," $CONTRL ICHARG=%d MULT=%d $END\n",geometries[i]->molecule.totalCharge,geometries[i]->molecule.spinMultiplicity);
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
static boolean createOrcaFiles(int numberOfGeometries, QuantumMechanicsModel** geometries, double* energies, char* fileNamePrefix, char* keyWords,char* orcaCommand)
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
		mol = &geometries[j]->molecule;
		fprintf(file,"! %s\n",keyWords);
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
/*****************************************************************************/
static void createPostProcessingFiles(int numberOfGeometries, QuantumMechanicsModel** geometries,double* energies,char* fileNameGeom, char* mopacKeywords, char* gaussianKeywords, char* fireflyKeywords, char* orcaKeywords,  char* cchemiKeywords, char* message, char* mopacCommand, char* gaussianCommand, char* fireflyCommand, char* orcaCommand, char* cchemiCommand)
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
	if(fireflyKeywords)
	{
		char* fileNamePrefix = getSuffixNameFile(fileNameGeom);
		createFireFlyFiles(numberOfGeometries, geometries, energies, fileNamePrefix, fireflyKeywords,fireflyCommand);
		strcat(message,fileNamePrefix);
		strcat(message,("FF_*.inp\n\tFiles for a post processing by FireFly\n\n"));
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
static int collectGeometriesFromProcessors(int nproc, QuantumMechanicsModel** geometriesAll, int numberOfGeometriesMax, double* energiesAll, QuantumMechanicsModel** geometries, int numberOfGeometries, double* energies,  double* coords,  double* enerDum, QuantumMechanicsModel* qmModel, double tolEnergy, double tolDistance)
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
			if(geometriesAll[i]) geometries[i]->klass->free(geometriesAll[i]);
			geometriesAll[i] = NULL;
		}
		for(i=0;i<numberOfGeometries;i++) 
		{
			if(geometries[i])
			{
				geometriesAll[i] = malloc(sizeof(QuantumMechanicsModel));
               			*geometriesAll[i] = geometries[i]->klass->copy(geometries[i]);
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
//#ifdef DEBUG
			printf("get nGeometries from proc n %d\n", j);
//#endif
			code = MPI_Recv(&nG,1,MPI_INT,j,tag,MPI_COMM_WORLD,&status) ;
//#ifdef DEBUG
			printf("nGeometries=%d from proc n %d\n",nG, j);
//#endif
			if(nG>0) 
			{
				int k;
				int a,b;
				tag = 2000;
				code = MPI_Recv(enerDum,nG,MPI_DOUBLE,j,tag,MPI_COMM_WORLD,&status) ;
				for(k=0;k<nG;k++)
				{
					int nA = qmModel->molecule.nAtoms;
					energiesAll[i+k] = enerDum[k];
					tag = 3000+k;
					code = MPI_Recv(coords,nA*3,MPI_DOUBLE,j,tag,MPI_COMM_WORLD,&status) ;
					geometriesAll[i+k] = malloc(sizeof(QuantumMechanicsModel));
               				*geometriesAll[i+k] = qmModel->klass->copy(qmModel);
					b = 0;
					for(a=0;a<nA;a++)
					{
//#ifdef DEBUG
						printf(" atoms %d C = %f %f %f\n",a,coords[b], coords[b+1], coords[b+2]);
//#endif
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
//#ifdef DEBUG
		printf("End collectGeometriesFromProcessors\n");
//#endif

#endif /* ENABLE_MPI */
		return numberOfGeometriesAll;
}
/*********************************************************************************************************************************/
static void sendGeometriesToMaster(int rank, QuantumMechanicsModel** geometries, int numberOfGeometries, double* energies,  double* coords, int nAtoms)
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
/***********************************************************************************************************************/
static void optAndSortGeometries(
	QuantumMechanicsModel* qmModel, 
	char* inputFileName,
	QuantumMechanicsModel*** geometries,
	int* numberOfGeometries,
	char* mopacCommand,
	char* gaussianCommand,
	char* N2P2Dir,
	char* tmModule,
	char* fireflyCommand,
	char* cchemiCommand,
	char* orcaCommand,
	char* openBabelCommand,
	char* genericCommand,
	char* optMopacMethod,
	char* optGaussianMethod,
	char* optFireFlyMethod,
	char* optOrcaMethod,
	char* optGenericMethod,
	double** energies,
	boolean optMopac,
	boolean optGaussian,
	boolean optGeneric,
	boolean optFireFly,
	boolean optOrca,
	boolean optOpenBabel,
	boolean optN2P2,
	boolean optTM,
	double tolEnergy,
	double tolDistance,
	char* fileNameGeom,
	boolean removeSimilarInertia,
	boolean removeFragmented,
	boolean removeSmallDistance,
	boolean removeSimilarBonds,
	double inertiaTol,
	double sTol,
	double distMaxTol,
	FILE* logfile,
	char* message
)
{
	/* minimazation by mopac PM6*/
	if(optMopac )
	{
		char* fileNamePrefix = getSuffixNameFile(fileNameGeom);
		char* keys = strdup_printf("%s XYZ",optMopacMethod);
		runMopacFiles(*numberOfGeometries, *geometries, *energies, fileNamePrefix, keys, mopacCommand);
		{
			char* fileNameGeomMop = strdup_printf("%sMop.gab",fileNamePrefix);
			sortGeometries(*numberOfGeometries, *geometries, *energies);
			saveEnergies(*energies, *numberOfGeometries, inputFileName);
			removeIdenticalGeometries(numberOfGeometries, geometries, energies, tolEnergy, tolDistance);
			if(removeSimilarInertia) qmModel->klass->removeSimilarInertiaGeometries(*geometries, numberOfGeometries,*energies,logfile,inertiaTol);
			if(removeSimilarBonds) qmModel->klass->removeSimilarBondsGeometries(*geometries, numberOfGeometries, *energies,logfile,sTol, distMaxTol);
			if(removeFragmented) qmModel->klass->removeFragmentedMolecules(*geometries, numberOfGeometries, *energies, logfile);
			if(removeSmallDistance) qmModel->klass->removeSmallDistanceMolecules(*geometries, numberOfGeometries, *energies, logfile);
			if(saveConfoGeometries(*numberOfGeometries, *geometries, *energies, fileNameGeom))
			{
				strcat(message,fileNameGeom);
				strcat(message,("\n\tGeometries after minimization by Mopac/"));
				strcat(message,optMopacMethod);
				strcat(message,("\n\tTo read this file through Gabedit : 'Read/CChemI file'\n\n"));
			}
			free(fileNameGeomMop);
		}
		if(fileNamePrefix) free(fileNamePrefix);
		if(keys)free(keys);
	}
	/* minimazation by gaussian PM6*/
	if(optGaussian )
	{
		char* fileNamePrefix = getSuffixNameFile(fileNameGeom);
		char* keys = strdup_printf("%s Opt",optGaussianMethod);
		if(runGaussianFiles(*numberOfGeometries, *geometries, *energies, fileNamePrefix, keys, gaussianCommand))
		{
			char* fileNameGeomGauss = strdup_printf("%sGauss.gab",fileNamePrefix);
			sortGeometries(*numberOfGeometries, *geometries, *energies);
			saveEnergies(*energies, *numberOfGeometries, inputFileName);
			removeIdenticalGeometries(numberOfGeometries, geometries, energies, tolEnergy, tolDistance);
			if(removeSimilarInertia) qmModel->klass->removeSimilarInertiaGeometries(*geometries, numberOfGeometries,*energies,logfile,inertiaTol);
			if(removeSimilarBonds) qmModel->klass->removeSimilarBondsGeometries(*geometries, numberOfGeometries, *energies,logfile,sTol, distMaxTol);
			if(removeFragmented) qmModel->klass->removeFragmentedMolecules(*geometries, numberOfGeometries, *energies, logfile);
			if(removeSmallDistance) qmModel->klass->removeSmallDistanceMolecules(*geometries, numberOfGeometries, *energies, logfile);
			if(saveConfoGeometries(*numberOfGeometries, *geometries, *energies, fileNameGeom))
			{
				strcat(message,fileNameGeom);
				strcat(message,("\n\tGeometries after minimization by Gaussian/"));
				strcat(message,optGaussianMethod);
				strcat(message,("\n\tTo read this file through Gabedit : 'Read/CChemI file'\n\n"));
			}
			free(fileNameGeomGauss);
		}
		if(fileNamePrefix) free(fileNamePrefix);
		if(keys)free(keys);
	}
	/* minimazation by FireFly*/
	if(optFireFly )
	{
		char* fileNamePrefix = getSuffixNameFile(fileNameGeom);
		char* keys=strdup_printf("RUNTYP=Optimize GBASIS=%s",optFireFlyMethod);
		if(runFireFlyFiles(*numberOfGeometries, *geometries, *energies, fileNamePrefix, keys,fireflyCommand))
		{
			char* fileNameGeomFireFly = strdup_printf("%sFireFly.gab",fileNamePrefix);
			sortGeometries(*numberOfGeometries, *geometries, *energies);
			saveEnergies(*energies, *numberOfGeometries, inputFileName);
			removeIdenticalGeometries(numberOfGeometries, geometries, energies, tolEnergy, tolDistance);
			if(removeSimilarInertia) qmModel->klass->removeSimilarInertiaGeometries(*geometries, numberOfGeometries,*energies,logfile,inertiaTol);
			if(removeSimilarBonds) qmModel->klass->removeSimilarBondsGeometries(*geometries, numberOfGeometries, *energies,logfile,sTol, distMaxTol);
			if(removeFragmented) qmModel->klass->removeFragmentedMolecules(*geometries, numberOfGeometries, *energies, logfile);
			if(removeSmallDistance) qmModel->klass->removeSmallDistanceMolecules(*geometries, numberOfGeometries, *energies, logfile);
			if(saveConfoGeometries(*numberOfGeometries, *geometries, *energies, fileNameGeom))
			{
				strcat(message,fileNameGeom);
				strcat(message,("\n\tGeometries after minimization by FireFly"));
				strcat(message,("\n\tTo read this file through Gabedit : 'Read/Gabedit file'\n\n"));
			}
			free(fileNameGeomFireFly);

		}
		if(fileNamePrefix) free(fileNamePrefix);
		if(keys)free(keys);
	}
	/* minimazation by OpenBabel*/
	if(optOpenBabel )
	{
		char* fileNamePrefix = getSuffixNameFile(fileNameGeom);
		char* keys=strdup_printf("opt");
		if(runOpenBabelFiles(*numberOfGeometries, *geometries, *energies, fileNamePrefix, keys,openBabelCommand))
		{
			char* fileNameGeomOpenBabel = strdup_printf("%sOpenBabel.gab",fileNamePrefix);
			sortGeometries(*numberOfGeometries, *geometries, *energies);
			saveEnergies(*energies, *numberOfGeometries, inputFileName);
			removeIdenticalGeometries(numberOfGeometries, geometries, energies, tolEnergy, tolDistance);
			if(removeSimilarInertia) qmModel->klass->removeSimilarInertiaGeometries(*geometries, numberOfGeometries,*energies,logfile,inertiaTol);
			if(removeSimilarBonds) qmModel->klass->removeSimilarBondsGeometries(*geometries, numberOfGeometries, *energies,logfile,sTol, distMaxTol);
			if(removeFragmented) qmModel->klass->removeFragmentedMolecules(*geometries, numberOfGeometries, *energies, logfile);
			if(removeSmallDistance) qmModel->klass->removeSmallDistanceMolecules(*geometries, numberOfGeometries, *energies, logfile);
			if(saveConfoGeometries(*numberOfGeometries, *geometries, *energies, fileNameGeom))
			{
				strcat(message,fileNameGeom);
				strcat(message,("\n\tGeometries after minimization by OpenBabel"));
				strcat(message,("\n\tTo read this file through Gabedit : 'Read/Gabedit file'\n\n"));
			}
			free(fileNameGeomOpenBabel);

		}
		if(fileNamePrefix) free(fileNamePrefix);
		if(keys)free(keys);
	}
	/* minimazation by Orca*/
	if(optOrca )
	{
		char* fileNamePrefix = getSuffixNameFile(fileNameGeom);
		char* keys=strdup_printf("%s opt",optOrcaMethod);
		if(runOrcaFiles(*numberOfGeometries, *geometries, *energies, fileNamePrefix, keys,orcaCommand))
		{
			char* fileNameGeomOrca = strdup_printf("%sOrca.gab",fileNamePrefix);
			sortGeometries(*numberOfGeometries, *geometries, *energies);
			saveEnergies(*energies, *numberOfGeometries, inputFileName);
			removeIdenticalGeometries(numberOfGeometries, geometries, energies, tolEnergy, tolDistance);
			if(removeSimilarInertia) qmModel->klass->removeSimilarInertiaGeometries(*geometries, numberOfGeometries,*energies,logfile,inertiaTol);
			if(removeSimilarBonds) qmModel->klass->removeSimilarBondsGeometries(*geometries, numberOfGeometries, *energies,logfile,sTol, distMaxTol);
			if(removeFragmented) qmModel->klass->removeFragmentedMolecules(*geometries, numberOfGeometries, *energies, logfile);
			if(removeSmallDistance) qmModel->klass->removeSmallDistanceMolecules(*geometries, numberOfGeometries, *energies, logfile);
			if(saveConfoGeometries(*numberOfGeometries, *geometries, *energies, fileNameGeom))
			{
				strcat(message,fileNameGeom);
				strcat(message,("\n\tGeometries after minimization by Orca"));
				strcat(message,("\n\tTo read this file through Gabedit : 'Read/Gabedit file'\n\n"));
			}
			free(fileNameGeomOrca);

		}
		if(fileNamePrefix) free(fileNamePrefix);
		if(keys)free(keys);
	}
	/* minimazation by generic*/
	if(optGeneric )
	{
		char* fileNamePrefix = getSuffixNameFile(fileNameGeom);
		char* keys = strdup_printf(" Opt");
		if(runGenericFiles(*numberOfGeometries, *geometries, *energies, fileNamePrefix, keys, optGenericMethod))
		{
			char* fileNameGeomGene = strdup_printf("%sGene.gab",fileNamePrefix);
			sortGeometries(*numberOfGeometries, *geometries, *energies);
			saveEnergies(*energies, *numberOfGeometries, inputFileName);
			removeIdenticalGeometries(numberOfGeometries, geometries, energies, tolEnergy, tolDistance);
			if(removeSimilarInertia) qmModel->klass->removeSimilarInertiaGeometries(*geometries, numberOfGeometries,*energies,logfile,inertiaTol);
			if(removeSimilarBonds) qmModel->klass->removeSimilarBondsGeometries(*geometries, numberOfGeometries, *energies,logfile,sTol, distMaxTol);
			if(removeFragmented) qmModel->klass->removeFragmentedMolecules(*geometries, numberOfGeometries, *energies, logfile);
			if(removeSmallDistance) qmModel->klass->removeSmallDistanceMolecules(*geometries, numberOfGeometries, *energies, logfile);
			if(saveConfoGeometries(*numberOfGeometries, *geometries, *energies, fileNameGeom))
			{
				strcat(message,fileNameGeom);
				strcat(message,("\n\tGeometries after minimization by Generic/"));
				strcat(message,optGenericMethod);
				strcat(message,("\n\tTo read this file through Gabedit : 'Read/CChemI file'\n\n"));
			}
			free(fileNameGeomGene);
		}
		if(fileNamePrefix) free(fileNamePrefix);
		if(keys)free(keys);
	}
	/* minimazation with N2P2  */
	//fprintf(stderr,"optN2P2=%d\n",optN2P2);
	if(optN2P2 )
	{
		fprintf(stderr,"optimisation using HDDN potentials\n");
		if(miminizeGeometriesUsingInternalOptimizer(*numberOfGeometries, *geometries, *energies, inputFileName))
		{
			sortGeometries(*numberOfGeometries, *geometries, *energies);
			saveEnergies(*energies, *numberOfGeometries, inputFileName);
			removeIdenticalGeometries(numberOfGeometries, geometries, energies, tolEnergy, tolDistance);
			if(removeSimilarInertia) qmModel->klass->removeSimilarInertiaGeometries(*geometries, numberOfGeometries,*energies,logfile,inertiaTol);
			if(removeSimilarBonds) qmModel->klass->removeSimilarBondsGeometries(*geometries, numberOfGeometries, *energies,logfile,sTol, distMaxTol);
			if(removeFragmented) qmModel->klass->removeFragmentedMolecules(*geometries, numberOfGeometries, *energies, logfile);
			if(removeSmallDistance) qmModel->klass->removeSmallDistanceMolecules(*geometries, numberOfGeometries, *energies, logfile);
			if(saveConfoGeometries(*numberOfGeometries, *geometries, *energies, fileNameGeom))
			{
				strcat(message,fileNameGeom);
				strcat(message,("\n\tGeometries after minimization by N2P2 potential"));
				strcat(message,("\n\tTo read this file through Gabedit : 'Read/Gabedit file'\n\n"));
			}

		}
	}
	/* minimazation with TM  */
	//fprintf(stderr,"optTM=%d\n",optTM);
	if(optTM )
	{
		fprintf(stderr,"optimisation using TM potential\n");
		if(miminizeGeometriesUsingInternalOptimizer(*numberOfGeometries, *geometries, *energies, inputFileName))
		{
			sortGeometries(*numberOfGeometries, *geometries, *energies);
			saveEnergies(*energies, *numberOfGeometries, inputFileName);
			removeIdenticalGeometries(numberOfGeometries, geometries, energies, tolEnergy, tolDistance);
			if(removeSimilarInertia) qmModel->klass->removeSimilarInertiaGeometries(*geometries, numberOfGeometries,*energies,logfile,inertiaTol);
			if(removeSimilarBonds) qmModel->klass->removeSimilarBondsGeometries(*geometries, numberOfGeometries, *energies,logfile,sTol, distMaxTol);
			if(removeFragmented) qmModel->klass->removeFragmentedMolecules(*geometries, numberOfGeometries, *energies, logfile);
			if(removeSmallDistance) qmModel->klass->removeSmallDistanceMolecules(*geometries, numberOfGeometries, *energies, logfile);
			if(saveConfoGeometries(*numberOfGeometries, *geometries, *energies, fileNameGeom))
			{
				strcat(message,fileNameGeom);
				strcat(message,("\n\tGeometries after minimization by TM potential"));
				strcat(message,("\n\tTo read this file through Gabedit : 'Read/Gabedit file'\n\n"));
			}

		}
	}
	if(!optMopac && !optFireFly && !optGaussian && !optGeneric && !optOpenBabel && !optN2P2 && !optTM)
	//if(!optMopac && !optFireFly && !optGaussian && !optGeneric )
	{
		/*  sort by energies */
		sortGeometries(*numberOfGeometries, *geometries, *energies);
		saveEnergies(*energies, *numberOfGeometries, inputFileName);
		removeIdenticalGeometries(numberOfGeometries, geometries, energies, tolEnergy, tolDistance);
		if(removeSimilarInertia) qmModel->klass->removeSimilarInertiaGeometries(*geometries, numberOfGeometries,*energies,logfile,inertiaTol);
		if(removeSimilarBonds) qmModel->klass->removeSimilarBondsGeometries(*geometries, numberOfGeometries, *energies,logfile,sTol, distMaxTol);
		if(removeFragmented) qmModel->klass->removeFragmentedMolecules(*geometries, numberOfGeometries, *energies, logfile);
		if(removeSmallDistance) qmModel->klass->removeSmallDistanceMolecules(*geometries, numberOfGeometries, *energies, logfile);
		/* printf("fileNameGeom = %s\n",fileNameGeom);*/
		if(saveConfoGeometries(*numberOfGeometries, *geometries, *energies, fileNameGeom))
		{
			strcat(message,fileNameGeom);
			strcat(message,("\n\tGeometries selected and optimized using your Quantum Mechanics potentials"));
			strcat(message,("\n\tTo read this file through Gabedit : 'Read/Gabedit file'\n\n"));
		}
	}
}
/***********************************************************************************************************************/
static void quantumMechanicsRemoveSimilarConfo(char* inputFileName)
{
	QuantumMechanicsModel qmModel; 
	QuantumMechanicsModel** geometries = NULL;
	int i;
	char message[BSIZE]="Created files :\n";
	char* dirName = NULL;
	int numberOfGeometries = 0;
	Molecule** mols = readMolecules(inputFileName,FALSE);
	char* mopacCommand = strdup("mopac");
	char* gaussianCommand=strdup("g09"); 
	char* N2P2Dir=strdup(getenv("PWD"));
	char* tmModule=strdup_printf("%s/%s",getenv("PWD"),"CChemITMModule");
	char* fireflyCommand=strdup("firefly");
	char* cchemiCommand=strdup("cchemi");
	char* orcaCommand = strdup("orca");
	char* openBabelCommand = strdup("obgradient");
	char* genericCommand=strdup("runGeneric");
	FILE* file = fopen(inputFileName,"rb");
	char* mopacKeywordsPost = NULL;
	char* gaussianKeywordsPost = NULL;
	char* fireflyKeywordsPost = NULL;
	char* orcaKeywordsPost = NULL;
	char* cchemiKeywordsPost = NULL;
	char* optMopacMethod=strdup("PM6");
	char* optGaussianMethod=strdup("AM1");
	char* optFireFlyMethod=strdup("AM1");
	char* optOrcaMethod=strdup("AM1");
	char* optGenericMethod=strdup("runGeneric");
	char* model = NULL;
	char* QMKeys = NULL;
	int cnt;
	char* optimizerType = NULL;

	double* energies = NULL;
	boolean optMopac = FALSE;
	boolean optGaussian = FALSE;
	boolean optGeneric = FALSE;
	boolean optFireFly = FALSE;
	boolean optOrca = FALSE;
	boolean optOpenBabel = FALSE;
	boolean optN2P2 = FALSE;
	boolean optTM = FALSE;
	double tolEnergy = 0.01;
	double tolDistance = 0.01;
	char* fileNameGeom = NULL;
	Constraints constraints = NOCONSTRAINTS;
	boolean chain=FALSE;
	boolean fragments=FALSE;
	boolean saveFirstGeom=FALSE;
	boolean removeSimilarInertia = FALSE;
	boolean removeFragmented = FALSE;
	boolean removeSmallDistance = FALSE;
	int nTimesGeoms=1;
	double inertiaTol = 0.04; // recomended byjun Zhao et al (2016)
       	//Comprehensive genetic algorithm for ab initio global optimisation of clusters, Molecular Simulation,
	// 42:10, 809-819, DOI: 10.1080/08927022.2015.1121386	     

	boolean removeSimilarBonds = FALSE;
	double sTol=0.02;
	double distMaxTol=0.7;
	// sTol = 0.02 , distMaxTol = 0.7 Ang, recommanded in Jorgensen et al JCTC, 2017
//Mathias S. Jørgensen , Michael N. Groves, and Bjørk Hammer
//J. Chem. Theory Comput., 2017, 13 (3), pp 1486–1493
//DOI: 10.1021/acs.jctc.6b01119
	FILE* logfile = stdout;

	numberOfGeometries = 0;
	while(mols && mols[numberOfGeometries]) numberOfGeometries++;

	readOneString(file,"mopacCommand",&mopacCommand);
	readOneString(file,"gaussianCommand",&gaussianCommand);
	readOneString(file,"genericCommand",&genericCommand);
	readOneString(file,"orcaCommand",&orcaCommand);
	readOneString(file,"openBabelCommand",&openBabelCommand);
	readOneString(file,"HDNNDir",&N2P2Dir);
	readOneString(file,"N2P2Dir",&N2P2Dir);
	readOneString(file,"tmModule",&tmModule);
	readOneString(file,"fireflyCommand",&fireflyCommand);
	readOneString(file,"cchemiCommand",&cchemiCommand);
	readOneReal(file,"tolEnergy",&tolEnergy);
	readOneReal(file,"tolDistance",&tolDistance);
	readOneBoolean(file,"ConfoOptMopac",&optMopac);
	readOneString(file,"ConfoOptMopacMethod",&optMopacMethod);
	readOneBoolean(file,"ConfoOptGaussian",&optGaussian);
	readOneString(file,"ConfoOptGaussianMethod",&optGaussianMethod);
	readOneBoolean(file,"ConfoOptFireFly",&optFireFly);
	readOneString(file,"ConfoOptFireFlyMethod",&optFireFlyMethod);
	readOneBoolean(file,"ConfoOptOrca",&optOrca);
	readOneString(file,"ConfoOptOrcaMethod",&optOrcaMethod);
	readOneBoolean(file,"ConfoOptOpenBabel",&optOpenBabel);
	readOneBoolean(file,"ConfoOptGeneric",&optGeneric);
	readOneString(file,"ConfoOptGenericMethod",&optGenericMethod);
	readOneBoolean(file,"ConfoOptN2P2",&optN2P2);
	readOneBoolean(file,"ConfoOptTM",&optTM);
	readOneBoolean(file,"RDChain",&chain);
	readOneBoolean(file,"RDFragments",&fragments);
	readOneBoolean(file,"RDSaveFirstGeom",&saveFirstGeom);
	readOneInt(file,"nTimesGeoms",&nTimesGeoms);

	if(readOneInt(file,"Constraints",&cnt)) constraints = cnt;

	readOneString(file,"mopacKeywordsPost",&mopacKeywordsPost);
	readOneString(file,"gaussianKeywordsPost",&gaussianKeywordsPost);
	readOneString(file,"fireflyKeywordsPost",&fireflyKeywordsPost);
	readOneString(file,"orcaKeywordsPost",&orcaKeywordsPost);
	readOneString(file,"cchemiKeywordsPost",&cchemiKeywordsPost);
	readOneString(file,"OptimizerType",&optimizerType);
	if(!optimizerType) optimizerType = strdup("External");
	if(!strstr(optimizerType,"External") && (optMopac||optFireFly))
	{
		printf("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n");
		printf("Sorry, The optimization using the internal optimizer with a RD conformational search is not yes implemented in this software\n");
		printf("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n");
		exit(1);
	}
	if(!readOneString(file,"Model",&model)) model = strdup("MOPAC");
	uppercase(model);
	if(!readOneString(file,"QMKeys",&QMKeys)) 
	{
		if(!strcmp(model,"MOPAC")) QMKeys = strdup("PM6-DH+");
		else if(!strcmp(model,"ORCA")) QMKeys = strdup("AM1");
		else if(!strcmp(model,"GAUSSIAN")) QMKeys = strdup("AM1");
		else if(!strcmp(model,"FIREFLY")) QMKeys = strdup("AM1");
		else if(!strcmp(model,"OPENBABEL")) QMKeys = strdup("mmff94");
		else QMKeys = strdup("SCC-DFTB");
	}
	if(!strcmp(model,"OPENBABEL")) 
	{
		int i; for(i=0;i<numberOfGeometries;i++) mols[i]->klass->buildMMTypes(mols[i], file);
	}

	/* fileName for geometries */
	{
		char* suff = getSuffixNameFile(inputFileName);
		dirName = strdup(getenv("PWD"));
		fileNameGeom = strdup_printf("%s%s",suff, "Geoms.gab");
		free(suff);
	}
	if(readOneBoolean(file,"removeSimilarInertia",&removeSimilarInertia) && removeSimilarInertia) readOneReal(file,"InertiaTol",&inertiaTol);
	readOneBoolean(file,"removeFragmented",&removeFragmented);
	readOneBoolean(file,"removeDissociated",&removeFragmented);
	readOneBoolean(file,"removeSmallDistance",&removeSmallDistance);
	if(readOneBoolean(file,"removeSimilarBonds",&removeSimilarBonds) && removeSimilarBonds)
	{
		readOneReal(file,"sTol",&sTol);
		readOneReal(file,"distMaxTol",&distMaxTol);
	}


	printf("Model = %s\n",model);
	printf("QMKeys = %s\n",QMKeys);

	{
		int i; 
		boolean addD3Correction;
		readOneBoolean(file,"addD3Correction",&addD3Correction);
		geometries = malloc(numberOfGeometries*sizeof(QuantumMechanicsModel*));
		for(i=0;i<numberOfGeometries;i++) 
		{
			geometries[i] = malloc(sizeof(QuantumMechanicsModel));
			QuantumMechanicsModel qmModel;
			qmModel = createGenericModel(mols[i], QMKeys, dirName, genericCommand, constraints, stdout);
                	*geometries[i] = qmModel.klass->copy(&qmModel);
			setH4Correction(file,geometries[i]);
			setSRBCorrection(file,geometries[i]);
			geometries[i]->addD3Correction=addD3Correction;
			checkWallCorrection(file, geometries[i]);
		}
	}



	int nOld = numberOfGeometries*nTimesGeoms;
	if(nTimesGeoms>1) numberOfGeometries = nOld;
	qmModel = *geometries[0];

	if(removeSimilarInertia) qmModel.klass->removeSimilarInertiaGeometries(geometries, &numberOfGeometries, energies,stdout,inertiaTol);
	if(removeFragmented) qmModel.klass->removeFragmentedMolecules(geometries, &numberOfGeometries, energies, stdout);
	if(removeSmallDistance) qmModel.klass->removeSmallDistanceMolecules(geometries, &numberOfGeometries, energies, stdout);
	if(removeSimilarBonds) qmModel.klass->removeSimilarBondsGeometries(geometries, &numberOfGeometries, energies,stdout,sTol, distMaxTol);

	if(nTimesGeoms>1) qmModel.klass->cutByInertia(geometries, &numberOfGeometries, energies,nOld/nTimesGeoms,stdout);

	if(geometries && numberOfGeometries>0)
	{
		int i;
		energies = malloc(numberOfGeometries*sizeof(double));
		for(i=0;i<numberOfGeometries;i++)
			energies[i] = geometries[i]->molecule.potentialEnergy;
	}

	optAndSortGeometries(&qmModel, inputFileName, &geometries, &numberOfGeometries, mopacCommand, gaussianCommand, N2P2Dir, tmModule, fireflyCommand, cchemiCommand,
	orcaCommand, openBabelCommand, genericCommand, optMopacMethod, optGaussianMethod, optFireFlyMethod, optOrcaMethod, optGenericMethod, &energies,
	optMopac, optGaussian, optGeneric, optFireFly, optOrca, optOpenBabel, optN2P2, optTM, tolEnergy, tolDistance, fileNameGeom, removeSimilarInertia,
	removeFragmented, removeSmallDistance, removeSimilarBonds, inertiaTol, sTol, distMaxTol, logfile, message);

	if(numberOfGeometries>0 && geometries )
	{
		char* fileNamePrefix = getSuffixNameFile(fileNameGeom);
		createPostProcessingFiles(numberOfGeometries, geometries, energies,fileNamePrefix, mopacKeywordsPost, gaussianKeywordsPost, fireflyKeywordsPost, orcaKeywordsPost, cchemiKeywordsPost, message, mopacCommand, gaussianCommand, fireflyCommand, orcaCommand, cchemiCommand);
		if(fileNamePrefix) free(fileNamePrefix);
	}
	if(geometries)
	{
		for(i=0;i<numberOfGeometries;i++)
			if(geometries[i]) geometries[i]->klass->free(geometries[i]);
		free(geometries);
	}
	if(energies) free(energies);
	if(strlen(message)>20) printf("%s",message);
	if(fileNameGeom)free(fileNameGeom);
	fclose(file);

}
/*****************************************************************************/
static void quantumMechanicsREMDConfo(char* inputFileName)
{
	QuantumMechanicsModel qmModel; 
	QuantumMechanicsMD molecularDynamics;
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
	double friction=-1;
	double omegaMax = 4000;
	int Nf = 50;
	double collide = 20;
	double qNH = 20;
	MDThermostatType thermostat = NONE;
	int numberOfGeometries = 2;
	QuantumMechanicsModel** geometries = NULL; 
	double* energies = NULL;
	boolean optMopac = FALSE;
	boolean optGaussian = FALSE;
	boolean optGeneric = FALSE;
	boolean optFireFly = FALSE;
	boolean optOrca = FALSE;
	boolean optOpenBabel = FALSE;
	boolean optN2P2 = FALSE;
	boolean optTM = FALSE;
	double tolEnergy = 0.01;
	double tolDistance = 0.01;
	Constraints constraints = NOCONSTRAINTS;
	int i;
	char message[BSIZE]="Created files :\n";
	char* dirName = NULL;
	Molecule mol = *(readMolecule(inputFileName,TRUE));
	char* mopacCommand = strdup("mopac");
	char* gaussianCommand=strdup("g09"); 
	char* N2P2Dir=strdup(getenv("PWD"));
	char* tmModule=strdup_printf("%s/%s",getenv("PWD"),"CChemITMModule");
	char* genericCommand=strdup("runGeneric"); 
	char* fireflyCommand=strdup("firefly");
	char* cchemiCommand=strdup("cchemi");
	char* orcaCommand = strdup("orca");
	char* openBabelCommand = strdup("obgradient");
	FILE* file = fopen(inputFileName,"rb");
	char* mopacKeywordsPost = NULL;
	char* gaussianKeywordsPost = NULL;
	char* fireflyKeywordsPost = NULL;
	char* orcaKeywordsPost = NULL;
	char* cchemiKeywordsPost = NULL;
	char* optMopacMethod=strdup("PM6");
	char* optGaussianMethod=strdup("AM1");
	char* optFireFlyMethod=strdup("AM1");
	char* optOrcaMethod=strdup("AM1");
	char* optGenericMethod=strdup("runGeneric");
	char* model = NULL;
	char* QMKeys = NULL;
	double runTempMax = 700;
	int nTemperatures = 10;
	int numberOfExchanges = 10;
	double timeExchange = 1;
	int nproc;
	int rank;
	QuantumMechanicsModel** geometriesAll = NULL; 
	int numberOfGeometriesAll = 2;
	int numberOfGeometriesMax = 2;
	double* energiesAll = NULL;
	double* coords = NULL;
	double* enerDum = NULL;
	FILE* logfile = NULL;
	int cnt = 0;
	char* optimizerType = NULL;
	boolean removeSimilarInertia = FALSE;
	boolean removeFragmented = FALSE;
	boolean removeSmallDistance = FALSE;
	double inertiaTol = 0.04; // recomended byjun Zhao et al (2016)
       	//Comprehensive genetic algorithm for ab initio global optimisation of clusters, Molecular Simulation,
	// 42:10, 809-819, DOI: 10.1080/08927022.2015.1121386	     

	boolean removeSimilarBonds = FALSE;
	double sTol=0.02;
	double distMaxTol=0.7;
	// sTol = 0.02 , distMaxTol = 0.7 Ang, recommanded in Jorgensen et al JCTC, 2017
//Mathias S. Jørgensen , Michael N. Groves, and Bjørk Hammer
//J. Chem. Theory Comput., 2017, 13 (3), pp 1486–1493
//DOI: 10.1021/acs.jctc.6b01119

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
		char* tmp = strdup_printf("cchemi_%d.log",rank);
		logfile = fopen(tmp,"w");
		free(tmp);
	}
	fprintf(logfile, "MolecularMechanicsDynamicsREMDConfoDlg Rank#=%d  nproc = %d\n", rank, nproc );
#ifdef DEBUG
	fprintf(logfile, "MolecularMechanicsDynamicsREMDConfoDlg Rank#=%d  nproc = %d\n", rank, nproc );
#endif

	setMDOptions(file, &updateFrequency, 
		&heatTime, &equiTime, &runTime, &coolTime,
		&heatTemp, &runTemp, &equiTemp, &coolTemp, &stepSize, 
		&integrator, &thermostat, &friction, &omegaMax, &Nf, &collide,&qNH);
	if(thermostat == NONE) 
	{
		fprintf(logfile, "Warning....................\n");
		fprintf(logfile, " A thermostat is required for a REMD calculation\n");
		fprintf(logfile, " I set it to Berendsen\n");
		thermostat = BERENDSEN;
	}
	if(readOneBoolean(file,"removeSimilarInertia",&removeSimilarInertia) && removeSimilarInertia) readOneReal(file,"InertiaTol",&inertiaTol);
	readOneBoolean(file,"removeFragmented",&removeFragmented);
	readOneBoolean(file,"removeDissociated",&removeFragmented);
	readOneBoolean(file,"removeSmallDistance",&removeSmallDistance);
	if(readOneBoolean(file,"removeSimilarBonds",&removeSimilarBonds) && removeSimilarBonds)
	{
		readOneReal(file,"sTol",&sTol);
		readOneReal(file,"distMaxTol",&distMaxTol);
	}


	readOneString(file,"mopacCommand",&mopacCommand);
	readOneString(file,"orcaCommand",&orcaCommand);
	readOneString(file,"openBabelCommand",&openBabelCommand);
	readOneString(file,"HDNNDir",&N2P2Dir);
	readOneString(file,"N2P2Dir",&N2P2Dir);
	readOneString(file,"tmModule",&tmModule);
	readOneString(file,"gaussianCommand",&gaussianCommand);
	readOneString(file,"fireflyCommand",&fireflyCommand);
	readOneString(file,"genericCommand",&genericCommand);
	readOneString(file,"cchemiCommand",&cchemiCommand);
	readOneReal(file,"tolEnergy",&tolEnergy);
	readOneReal(file,"tolDistance",&tolDistance);
	readOneBoolean(file,"ConfoOptMopac",&optMopac);
	readOneString(file,"ConfoOptMopacMethod",&optMopacMethod);
	readOneBoolean(file,"ConfoOptGaussian",&optGaussian);
	readOneString(file,"ConfoOptGaussianMethod",&optGaussianMethod);
	readOneBoolean(file,"ConfoOptGeneric",&optGeneric);
	readOneString(file,"ConfoOptGenericMethod",&optGenericMethod);
	readOneBoolean(file,"ConfoOptFireFly",&optFireFly);
	readOneString(file,"ConfoOptFireFlyMethod",&optFireFlyMethod);
	readOneBoolean(file,"ConfoOptOrca",&optOrca);
	readOneString(file,"ConfoOptOrcaMethod",&optOrcaMethod);
	readOneBoolean(file,"ConfoOptOpenBabel",&optOpenBabel);
	readOneBoolean(file,"ConfoOptN2P2",&optN2P2);
	readOneBoolean(file,"ConfoOptTM",&optTM);

	if(readOneInt(file,"Constraints",&cnt)) constraints = cnt;

	readOneString(file,"mopacKeywordsPost",&mopacKeywordsPost);
	readOneString(file,"gaussianKeywordsPost",&gaussianKeywordsPost);
	readOneString(file,"fireflyKeywordsPost",&fireflyKeywordsPost);
	readOneString(file,"orcaKeywordsPost",&orcaKeywordsPost);
	readOneString(file,"cchemiKeywordsPost",&cchemiKeywordsPost);
	readOneString(file,"OptimizerType",&optimizerType);
	if(!optimizerType) optimizerType = strdup("External");
	if(!strstr(optimizerType,"External")  && (optMopac||optFireFly))
	{
		printf("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n");
		printf("Sorry, The optimization using the internal optimizer after a MD conformational search is not yes implemented in this software\n");
		printf("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n");
		exit(1);
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


	if(numberOfExchanges<1) numberOfExchanges=2;
	/* number for geometries */
	{
		numberOfGeometries = 10;
		readOneInt(file,"numberOfGeometries",&numberOfGeometries);
		if(numberOfGeometries<2) numberOfGeometries = 2;
		fprintf(logfile, "numerOfGeometries = %d\n",numberOfGeometries);
	}
	if(!readOneString(file,"Model",&model)) model = strdup("MOPAC");
	uppercase(model);
	if(!readOneString(file,"QMKeys",&QMKeys)) 
	{
		if(!strcmp(model,"MOPAC")) QMKeys = strdup("PM6-DH+");
		else if(!strcmp(model,"ORCA")) QMKeys = strdup("AM1");
		else if(!strcmp(model,"GAUSSIAN")) QMKeys = strdup("AM1");
		else if(!strcmp(model,"FIREFLY")) QMKeys = strdup("AM1");
		else if(!strcmp(model,"OPENBABEL")) QMKeys = strdup("mmff94");
		else QMKeys = strdup("SCC-DFTB");
	}
	if(!strcmp(model,"OPENBABEL")) mol.klass->buildMMTypes(&mol, file);

	/* fileName for geometries */
	{
		char* suff = getSuffixNameFile(inputFileName);
		dirName = strdup(getenv("PWD"));
		fileNameGeom = strdup_printf("%s%s",suff, "Geoms.gab");
		fileNameTraj = strdup_printf("%s%s",suff, "Traj.gab");
		fileNameProp = strdup_printf("%s%s",suff, "Prop.txt");
		free(suff);
	}


	fprintf(logfile, "Model = %s\n",model);
	fprintf(logfile, "QMKeys = %s\n",QMKeys);
	if(!strcmp(model,"MOPAC")) qmModel = createMopacModel(&mol, QMKeys, dirName, mopacCommand, constraints, logfile);
	else if(!strcmp(model,"FIREFLY")) qmModel = createFireFlyModel(&mol, QMKeys, dirName, fireflyCommand, constraints, logfile);
	else if(!strcmp(model,"ORCA")) qmModel = createOrcaModel(&mol, QMKeys, dirName, orcaCommand, constraints, logfile);
	else if(!strcmp(model,"OPENBABEL")) qmModel = createOpenBabelModel(&mol, QMKeys, dirName, openBabelCommand, constraints, logfile);
	else if(!strcmp(model,"GAUSSIAN")) qmModel = createGaussianModel(&mol, QMKeys, dirName, gaussianCommand, constraints, logfile);
	else if(!strcmp(model,"N2P2") || !strcmp(model,"HDNN")) qmModel = createN2P2Model(&mol, N2P2Dir, constraints, logfile);
	else if(!strcmp(model,"TM") || !strcmp(model,"TensorModel")) qmModel = createTMModel(&mol, tmModule, constraints, logfile);
	else qmModel = createGenericModel(&mol, QMKeys, dirName, genericCommand, constraints, logfile);

	setH4Correction(file,&qmModel);
	setSRBCorrection(file,&qmModel);
	readOneBoolean(file,"addD3Correction",&qmModel.addD3Correction);
	checkWallCorrection(file, &qmModel);

	geometries = runQuantumMechanicsREMDConfo(&molecularDynamics, &qmModel,
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
	fprintf(logfile, "Alloc energiesAll rank= %d\n",rank);
	fflush(logfile);
	if(nproc>1)
	{
		coords =  malloc(qmModel.molecule.nAtoms*3*sizeof(double));
		enerDum =  malloc(numberOfGeometriesMax*sizeof(double));
		if(geometries && numberOfGeometries>0) energiesAll = malloc(numberOfGeometriesMax*sizeof(double));
		if(geometries && numberOfGeometries>0) 
		{
			geometriesAll = malloc(numberOfGeometriesMax*sizeof(QuantumMechanicsModel*));
			for(i=0;i<numberOfGeometriesMax;i++) 
			{
				geometriesAll[i] = NULL;
			}
		}
	}

	fprintf(logfile, "Alloc energies rank= %d\n",rank);
	fflush(logfile);
	if(geometries && numberOfGeometries>0)
	{
		int i;
		energies = malloc(numberOfGeometries*sizeof(double));
		for(i=0;i<numberOfGeometries;i++) energies[i] = 1e30;
		for(i=0;i<numberOfGeometries;i++) 
			if(geometries[i]) energies[i] = geometries[i]->molecule.potentialEnergy;
	}
	fprintf(logfile, "End Alloc energies rank= %d\n",rank);
	fflush(logfile);

	/* minimazation by mopac PM6*/
	if(optMopac )
	{
		char* fileNamePrefix = getSuffixNameFile(fileNameGeom);
		char* keys = strdup_printf("%s XYZ",optMopacMethod);
		runMopacFiles(numberOfGeometries, geometries, energies, fileNamePrefix, keys, mopacCommand);
		{
			char* fileNameGeomMop = strdup_printf("%sMop.gab",fileNamePrefix);
			sortGeometries(numberOfGeometries, geometries, energies);
			removeIdenticalGeometries(&numberOfGeometries, &geometries, &energies, tolEnergy, tolDistance);
			if(removeSimilarInertia) qmModel.klass->removeSimilarInertiaGeometries(geometries, &numberOfGeometries,energies,stdout,inertiaTol);
			if(removeSimilarBonds) qmModel.klass->removeSimilarBondsGeometries(geometries, &numberOfGeometries, energies,stdout,sTol, distMaxTol);
			if(removeFragmented) qmModel.klass->removeFragmentedMolecules(geometries, &numberOfGeometries, energies, stdout);
			if(removeSmallDistance) qmModel.klass->removeSmallDistanceMolecules(geometries, &numberOfGeometries, energies, stdout);
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
					numberOfGeometriesAll = collectGeometriesFromProcessors(nproc, geometriesAll, numberOfGeometriesMax, energiesAll, geometries, numberOfGeometries, energies,  coords,  enerDum, &qmModel,tolEnergy,tolDistance);
				else sendGeometriesToMaster(rank, geometries, numberOfGeometries, energies,  coords, qmModel.molecule.nAtoms);
				fprintf(logfile, "End collect&send nGeoms = %d\n", numberOfGeometriesAll );
				fflush(logfile);
			}
			if(rank==0 && saveConfoGeometries(numberOfGeometriesAll, geometriesAll, energiesAll, fileNameGeom))
			{
				strcat(message,fileNameGeom);
				strcat(message,("\n\tGeometries after minimization by Mopac/"));
				strcat(message,optMopacMethod);
				strcat(message,("\n\tTo read this file through Gabedit : 'Read/CChemI file'\n\n"));
			}
			free(fileNameGeomMop);
		}
		if(fileNamePrefix) free(fileNamePrefix);
		if(keys)free(keys);
	}
	/* minimazation by Gaussian PM6*/
	if(optGaussian )
	{
		char* fileNamePrefix = getSuffixNameFile(fileNameGeom);
		char* keys = strdup_printf("%s Opt",optGaussianMethod);
		if(runGaussianFiles(numberOfGeometries, geometries, energies, fileNamePrefix, keys, gaussianCommand))
		{
			char* fileNameGeomGauss = strdup_printf("%sGauss.gab",fileNamePrefix);
			sortGeometries(numberOfGeometries, geometries, energies);
			removeIdenticalGeometries(&numberOfGeometries, &geometries, &energies, tolEnergy, tolDistance);
			if(removeSimilarInertia) qmModel.klass->removeSimilarInertiaGeometries(geometries, &numberOfGeometries,energies,stdout,inertiaTol);
			if(removeSimilarBonds) qmModel.klass->removeSimilarBondsGeometries(geometries, &numberOfGeometries, energies,stdout,sTol, distMaxTol);
			if(removeFragmented) qmModel.klass->removeFragmentedMolecules(geometries, &numberOfGeometries, energies, stdout);
			if(removeSmallDistance) qmModel.klass->removeSmallDistanceMolecules(geometries, &numberOfGeometries, energies, stdout);
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
					numberOfGeometriesAll = collectGeometriesFromProcessors(nproc, geometriesAll, numberOfGeometriesMax, energiesAll, geometries, numberOfGeometries, energies,  coords,  enerDum, &qmModel,tolEnergy,tolDistance);
				else sendGeometriesToMaster(rank, geometries, numberOfGeometries, energies,  coords, qmModel.molecule.nAtoms);
				fprintf(logfile, "End collect&send nGeoms = %d\n", numberOfGeometriesAll );
				fflush(logfile);
			}
			if(rank==0 && saveConfoGeometries(numberOfGeometriesAll, geometriesAll, energiesAll, fileNameGeom))
			{
				strcat(message,fileNameGeom);
				strcat(message,("\n\tGeometries after minimization by Gaussian/"));
				strcat(message,optGaussianMethod);
				strcat(message,("\n\tTo read this file through Gabedit : 'Read/CChemI file'\n\n"));
			}
			free(fileNameGeomGauss);
		}
		if(fileNamePrefix) free(fileNamePrefix);
		if(keys)free(keys);
	}
	/* minimazation by FireFly*/
	if(optFireFly )
	{
		char* fileNamePrefix = getSuffixNameFile(fileNameGeom);
		char* keys=strdup_printf("RUNTYP=Optimize GBASIS=%s",optFireFlyMethod);
		if(runFireFlyFiles(numberOfGeometries, geometries, energies, fileNamePrefix, keys,fireflyCommand))
		{
			char* fileNameGeomFireFly = strdup_printf("%sFireFly.gab",fileNamePrefix);
			sortGeometries(numberOfGeometries, geometries, energies);
			removeIdenticalGeometries(&numberOfGeometries, &geometries, &energies, tolEnergy, tolDistance);
			if(removeSimilarInertia) qmModel.klass->removeSimilarInertiaGeometries(geometries, &numberOfGeometries,energies,stdout,inertiaTol);
			if(removeSimilarBonds) qmModel.klass->removeSimilarBondsGeometries(geometries, &numberOfGeometries, energies,stdout,sTol, distMaxTol);
			if(removeFragmented) qmModel.klass->removeFragmentedMolecules(geometries, &numberOfGeometries, energies, stdout);
			if(removeSmallDistance) qmModel.klass->removeSmallDistanceMolecules(geometries, &numberOfGeometries, energies, stdout);
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
					numberOfGeometriesAll = collectGeometriesFromProcessors(nproc, geometriesAll, numberOfGeometriesMax, energiesAll, geometries, numberOfGeometries, energies,  coords,  enerDum, &qmModel,tolEnergy,tolDistance);
				else sendGeometriesToMaster(rank, geometries, numberOfGeometries, energies,  coords, qmModel.molecule.nAtoms);
				fprintf(logfile, "End collect&send nGeoms = %d\n", numberOfGeometriesAll );
				fflush(logfile);
			}
			if(rank==0 && saveConfoGeometries(numberOfGeometriesAll, geometriesAll, energiesAll, fileNameGeom))
			{
				strcat(message,fileNameGeom);
				strcat(message,("\n\tGeometries after minimization by FireFly"));
				strcat(message,("\n\tTo read this file through Gabedit : 'Read/Gabedit file'\n\n"));
			}
			free(fileNameGeomFireFly);

		}
		if(fileNamePrefix) free(fileNamePrefix);
		if(keys)free(keys);
	}
	/* minimazation by generic*/
	if(optGeneric )
	{
		char* fileNamePrefix = getSuffixNameFile(fileNameGeom);
		char* keys = strdup_printf(" Opt");
		runGenericFiles(numberOfGeometries, geometries, energies, fileNamePrefix, keys, optGenericMethod);
		{
			char* fileNameGeomGen = strdup_printf("%sGen.gab",fileNamePrefix);
			sortGeometries(numberOfGeometries, geometries, energies);
			removeIdenticalGeometries(&numberOfGeometries, &geometries, &energies, tolEnergy, tolDistance);
			if(removeSimilarInertia) qmModel.klass->removeSimilarInertiaGeometries(geometries, &numberOfGeometries,energies,stdout,inertiaTol);
			if(removeSimilarBonds) qmModel.klass->removeSimilarBondsGeometries(geometries, &numberOfGeometries, energies,stdout,sTol, distMaxTol);
			if(removeFragmented) qmModel.klass->removeFragmentedMolecules(geometries, &numberOfGeometries, energies, stdout);
			if(removeSmallDistance) qmModel.klass->removeSmallDistanceMolecules(geometries, &numberOfGeometries, energies, stdout);
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
					numberOfGeometriesAll = collectGeometriesFromProcessors(nproc, geometriesAll, numberOfGeometriesMax, energiesAll, geometries, numberOfGeometries, energies,  coords,  enerDum, &qmModel,tolEnergy,tolDistance);
				else sendGeometriesToMaster(rank, geometries, numberOfGeometries, energies,  coords, qmModel.molecule.nAtoms);
				fprintf(logfile, "End collect&send nGeoms = %d\n", numberOfGeometriesAll );
				fflush(logfile);
			}
			if(rank==0 && saveConfoGeometries(numberOfGeometriesAll, geometriesAll, energiesAll, fileNameGeom))
			{
				strcat(message,fileNameGeom);
				strcat(message,("\n\tGeometries after minimization by Generic/"));
				strcat(message,optGenericMethod);
				strcat(message,("\n\tTo read this file through Gabedit : 'Read/CChemI file'\n\n"));
			}
			free(fileNameGeomGen);
		}
		if(fileNamePrefix) free(fileNamePrefix);
		if(keys)free(keys);
	}
	/* minimazation by OpenBabel*/
	if(optOpenBabel )
	{
		char* fileNamePrefix = getSuffixNameFile(fileNameGeom);
		char* keys=strdup_printf("opt");
		if(runOpenBabelFiles(numberOfGeometries, geometries, energies, fileNamePrefix, keys,openBabelCommand))
		{
			char* fileNameGeomOpenBabel = strdup_printf("%sOpenBabel.gab",fileNamePrefix);
			sortGeometries(numberOfGeometries, geometries, energies);
			removeIdenticalGeometries(&numberOfGeometries, &geometries, &energies, tolEnergy, tolDistance);
			if(removeSimilarInertia) qmModel.klass->removeSimilarInertiaGeometries(geometries, &numberOfGeometries,energies,stdout,inertiaTol);
			if(removeSimilarBonds) qmModel.klass->removeSimilarBondsGeometries(geometries, &numberOfGeometries, energies,stdout,sTol, distMaxTol);
			if(removeFragmented) qmModel.klass->removeFragmentedMolecules(geometries, &numberOfGeometries, energies, stdout);
			if(removeSmallDistance) qmModel.klass->removeSmallDistanceMolecules(geometries, &numberOfGeometries, energies, stdout);
			if(saveConfoGeometries(numberOfGeometries, geometries, energies, fileNameGeom))
			{
				strcat(message,fileNameGeom);
				strcat(message,("\n\tGeometries after minimization by OpenBabel"));
				strcat(message,("\n\tTo read this file through Gabedit : 'Read/Gabedit file'\n\n"));
			}
			free(fileNameGeomOpenBabel);

		}
		if(fileNamePrefix) free(fileNamePrefix);
		if(keys)free(keys);
	}
	/* minimazation by Orca*/
	if(optOrca )
	{
		char* fileNamePrefix = getSuffixNameFile(fileNameGeom);
		char* keys=strdup_printf("%s opt",optOrcaMethod);
		if(runOrcaFiles(numberOfGeometries, geometries, energies, fileNamePrefix, keys,orcaCommand))
		{
			char* fileNameGeomOrca = strdup_printf("%sOrca.gab",fileNamePrefix);
			sortGeometries(numberOfGeometries, geometries, energies);
			removeIdenticalGeometries(&numberOfGeometries, &geometries, &energies, tolEnergy, tolDistance);
			if(removeSimilarInertia) qmModel.klass->removeSimilarInertiaGeometries(geometries, &numberOfGeometries,energies,stdout,inertiaTol);
			if(removeSimilarBonds) qmModel.klass->removeSimilarBondsGeometries(geometries, &numberOfGeometries, energies,stdout,sTol, distMaxTol);
			if(removeFragmented) qmModel.klass->removeFragmentedMolecules(geometries, &numberOfGeometries, energies, stdout);
			if(removeSmallDistance) qmModel.klass->removeSmallDistanceMolecules(geometries, &numberOfGeometries, energies, stdout);
			if(saveConfoGeometries(numberOfGeometries, geometries, energies, fileNameGeom))
			{
				strcat(message,fileNameGeom);
				strcat(message,("\n\tGeometries after minimization by Orca"));
				strcat(message,("\n\tTo read this file through Gabedit : 'Read/Gabedit file'\n\n"));
			}
			free(fileNameGeomOrca);

		}
		if(fileNamePrefix) free(fileNamePrefix);
		if(keys)free(keys);
	}
	/* minimazation with N2P2  */
	//fprintf(stderr,"optN2P2=%d\n",optN2P2);
	if(optN2P2 )
	{
		fprintf(stderr,"optimisation using HDDN potentials\n");
		if(miminizeGeometriesUsingInternalOptimizer(numberOfGeometries, geometries, energies, inputFileName))
		{
			sortGeometries(numberOfGeometries, geometries, energies);
			saveEnergies(energies, numberOfGeometries, inputFileName);
			removeIdenticalGeometries(&numberOfGeometries, &geometries, &energies, tolEnergy, tolDistance);
			if(removeSimilarInertia) qmModel.klass->removeSimilarInertiaGeometries(geometries, &numberOfGeometries,energies,stdout,inertiaTol);
			if(removeSimilarBonds) qmModel.klass->removeSimilarBondsGeometries(geometries, &numberOfGeometries, energies,stdout,sTol, distMaxTol);
			if(removeFragmented) qmModel.klass->removeFragmentedMolecules(geometries, &numberOfGeometries, energies, stdout);
			if(removeSmallDistance) qmModel.klass->removeSmallDistanceMolecules(geometries, &numberOfGeometries, energies, stdout);
			if(saveConfoGeometries(numberOfGeometries, geometries, energies, fileNameGeom))
			{
				strcat(message,fileNameGeom);
				strcat(message,("\n\tGeometries after minimization by N2P2 potential"));
				strcat(message,("\n\tTo read this file through Gabedit : 'Read/Gabedit file'\n\n"));
			}

		}
	}
	/* minimazation with TM  */
	//fprintf(stderr,"optTM=%d\n",optTM);
	if(optTM )
	{
		fprintf(stderr,"optimisation using TM potential\n");
		if(miminizeGeometriesUsingInternalOptimizer(numberOfGeometries, geometries, energies, inputFileName))
		{
			sortGeometries(numberOfGeometries, geometries, energies);
			saveEnergies(energies, numberOfGeometries, inputFileName);
			removeIdenticalGeometries(&numberOfGeometries, &geometries, &energies, tolEnergy, tolDistance);
			if(removeSimilarInertia) qmModel.klass->removeSimilarInertiaGeometries(geometries, &numberOfGeometries,energies,stdout,inertiaTol);
			if(removeSimilarBonds) qmModel.klass->removeSimilarBondsGeometries(geometries, &numberOfGeometries, energies,stdout,sTol, distMaxTol);
			if(removeFragmented) qmModel.klass->removeFragmentedMolecules(geometries, &numberOfGeometries, energies, stdout);
			if(removeSmallDistance) qmModel.klass->removeSmallDistanceMolecules(geometries, &numberOfGeometries, energies, stdout);
			if(saveConfoGeometries(numberOfGeometries, geometries, energies, fileNameGeom))
			{
				strcat(message,fileNameGeom);
				strcat(message,("\n\tGeometries after minimization by TM potential"));
				strcat(message,("\n\tTo read this file through Gabedit : 'Read/Gabedit file'\n\n"));
			}

		}
	}
	if(!optMopac && !optFireFly && !optGaussian && !optGeneric && !optOpenBabel && !optN2P2 && !optTM)
	{
		/*  sort by energies */
		fprintf(logfile, "Begin sortGeom\n");
		fflush(logfile);
		sortGeometries(numberOfGeometries, geometries, energies);
		fprintf(logfile, "End sortGeom\n");
		fflush(logfile);
		removeIdenticalGeometries(&numberOfGeometries, &geometries, &energies, tolEnergy, tolDistance);
		if(removeSimilarInertia) qmModel.klass->removeSimilarInertiaGeometries(geometries, &numberOfGeometries,energies,stdout,inertiaTol);
		if(removeSimilarBonds) qmModel.klass->removeSimilarBondsGeometries(geometries, &numberOfGeometries, energies,stdout,sTol, distMaxTol);
		if(removeFragmented) qmModel.klass->removeFragmentedMolecules(geometries, &numberOfGeometries, energies, stdout);
		if(removeSmallDistance) qmModel.klass->removeSmallDistanceMolecules(geometries, &numberOfGeometries, energies, stdout);
		fprintf(logfile, "End removeId\n");
		fflush(logfile);
		/* printf("fileNameGeom = %s\n",fileNameGeom);*/
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
			fflush(logfile);
//#endif
			if(rank==0) 
				numberOfGeometriesAll = collectGeometriesFromProcessors(nproc, geometriesAll, numberOfGeometriesMax, energiesAll, geometries, numberOfGeometries, energies,  coords,  enerDum, &qmModel,tolEnergy,tolDistance);
			else sendGeometriesToMaster(rank, geometries, numberOfGeometries, energies,  coords, qmModel.molecule.nAtoms);
			fprintf(logfile, "End collect&send nGeoms = %d\n", numberOfGeometriesAll );
			fflush(logfile);
		}
		if(rank==0 && saveConfoGeometries(numberOfGeometriesAll, geometriesAll, energiesAll, fileNameGeom))
		{
			strcat(message,fileNameGeom);
			strcat(message,("\n\tGeometries selected and optimized using your Quantum Mechanics potentials"));
			strcat(message,("\n\tTo read this file through Gabedit : 'Read/Gabedit file'\n\n"));
		}
	}
	if(rank==0 && numberOfGeometriesAll>0 && geometriesAll )
	{
		createPostProcessingFiles(numberOfGeometriesAll, geometriesAll, energiesAll,fileNameGeom, mopacKeywordsPost, gaussianKeywordsPost, fireflyKeywordsPost, orcaKeywordsPost, cchemiKeywordsPost, message, mopacCommand, gaussianCommand, fireflyCommand, orcaCommand, cchemiCommand);
/*
		if(mopacKeywordsPost)
		{
			char* fileNamePrefix = getSuffixNameFile(fileNameGeom);
			createMopacFiles(numberOfGeometriesAll, geometriesAll, energiesAll, fileNamePrefix, mopacKeywordsPost,mopacCommand);
			strcat(message,fileNamePrefix);
			strcat(message,("_*.mop\n\tFiles for a post processing by Mopac\n\n"));
			if(fileNamePrefix) free(fileNamePrefix);
		}
		if(gaussianKeywordsPost)
		{
			char* fileNamePrefix = getSuffixNameFile(fileNameGeom);
			createGaussianFiles(numberOfGeometriesAll, geometriesAll, energiesAll, fileNamePrefix, gaussianKeywordsPost, gaussianCommand);
			strcat(message,fileNamePrefix);
			strcat(message,("_*.com\n\tFiles for a post processing by Gaussian\n\n"));
			if(fileNamePrefix) free(fileNamePrefix);
		}
		if(fireflyKeywordsPost)
		{
			char* fileNamePrefix = getSuffixNameFile(fileNameGeom);
			createFireFlyFiles(numberOfGeometriesAll, geometriesAll, energiesAll, fileNamePrefix, fireflyKeywordsPost,fireflyCommand);
			strcat(message,fileNamePrefix);
			strcat(message,("FF_*.inp\n\tFiles for a post processing by FireFly\n\n"));
			if(fileNamePrefix) free(fileNamePrefix);
		}
*/
	}
	fprintf(logfile, "Begin freegeometriesAll\n" );
	fflush(logfile);
	if(geometriesAll && geometriesAll!=geometries)
	{
		for(i=0;i<numberOfGeometriesAll;i++)
			if(geometriesAll[i]) geometries[i]->klass->free(geometriesAll[i]);
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
			if(geometries[i]) geometries[i]->klass->free(geometries[i]);
		free(geometries);
	}
	if(energies) free(energies);
	if(strlen(message)>20) printf("%s",message);
	if(fileNameGeom)free(fileNameGeom);
	qmModel.klass->free(&qmModel);
	fclose(file);

}
/*****************************************************************************/
static void quantumMechanicsMDConfo(char* inputFileName)
{
	QuantumMechanicsModel qmModel; 
	QuantumMechanicsMD molecularDynamics;
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
	double friction=-1;
	double collide = 20;
	double qNH = 20;
	MDThermostatType thermostat = NONE;
	int numberOfGeometries = 2;
	QuantumMechanicsModel** geometries = NULL; 
	double* energies = NULL;
	boolean optMopac = FALSE;
	boolean optGaussian = FALSE;
	boolean optGeneric = FALSE;
	boolean optFireFly = FALSE;
	boolean optOrca = FALSE;
	boolean optOpenBabel = FALSE;
	boolean optN2P2 = FALSE;
	boolean optTM = FALSE;
	double tolEnergy = 0.01;
	double tolDistance = 0.01;
	double omegaMax = 4000;
	int Nf = 50;
	Constraints constraints = NOCONSTRAINTS;
	int i;
	char message[BSIZE]="Created files :\n";
	char* dirName = NULL;
	Molecule mol = *(readMolecule(inputFileName,TRUE));
	char* mopacCommand = strdup("mopac");
	char* gaussianCommand=strdup("g09"); 
	char* N2P2Dir=strdup(getenv("PWD"));
	char* tmModule=strdup_printf("%s/%s",getenv("PWD"),"CChemITMModule");
	char* fireflyCommand=strdup("firefly");
	char* cchemiCommand=strdup("cchemi");
	char* orcaCommand = strdup("orca");
	char* openBabelCommand = strdup("obgradient");
	char* genericCommand=strdup("runGeneric");
	FILE* file = fopen(inputFileName,"rb");
	char* mopacKeywordsPost = NULL;
	char* gaussianKeywordsPost = NULL;
	char* fireflyKeywordsPost = NULL;
	char* orcaKeywordsPost = NULL;
	char* cchemiKeywordsPost = NULL;
	char* optMopacMethod=strdup("PM6");
	char* optGaussianMethod=strdup("AM1");
	char* optFireFlyMethod=strdup("AM1");
	char* optOrcaMethod=strdup("AM1");
	char* optGenericMethod=strdup("runGeneric");
	char* model = NULL;
	char* QMKeys = NULL;
	FILE* logfile = stdout;
	int cnt;
	char* optimizerType = NULL;
	boolean removeSimilarInertia = FALSE;
	boolean removeFragmented = FALSE;
	boolean removeSmallDistance = FALSE;
	double inertiaTol = 0.04; // recomended byjun Zhao et al (2016)
       	//Comprehensive genetic algorithm for ab initio global optimisation of clusters, Molecular Simulation,
	// 42:10, 809-819, DOI: 10.1080/08927022.2015.1121386	     

	boolean removeSimilarBonds = FALSE;
	double sTol=0.02;
	double distMaxTol=0.7;
	// sTol = 0.02 , distMaxTol = 0.7 Ang, recommanded in Jorgensen et al JCTC, 2017
//Mathias S. Jørgensen , Michael N. Groves, and Bjørk Hammer
//J. Chem. Theory Comput., 2017, 13 (3), pp 1486–1493
//DOI: 10.1021/acs.jctc.6b01119


	printf("MDConfo begin\n");
	setMDOptions(file, &updateFrequency, 
		&heatTime, &equiTime, &runTime, &coolTime,
		&heatTemp, &runTemp, &equiTemp, &coolTemp, &stepSize, 
		&integrator, &thermostat, &friction, &omegaMax, &Nf, &collide,&qNH);
	if(readOneBoolean(file,"removeSimilarInertia",&removeSimilarInertia) && removeSimilarInertia) readOneReal(file,"InertiaTol",&inertiaTol);
	readOneBoolean(file,"removeFragmented",&removeFragmented);
	readOneBoolean(file,"removeDissociated",&removeFragmented);
	readOneBoolean(file,"removeSmallDistance",&removeSmallDistance);
	if(readOneBoolean(file,"removeSimilarBonds",&removeSimilarBonds) && removeSimilarBonds)
	{
		readOneReal(file,"sTol",&sTol);
		readOneReal(file,"distMaxTol",&distMaxTol);
	}

	readOneString(file,"mopacCommand",&mopacCommand);
	readOneString(file,"gaussianCommand",&gaussianCommand);
	readOneString(file,"genericCommand",&genericCommand);
	readOneString(file,"orcaCommand",&orcaCommand);
	readOneString(file,"openBabelCommand",&openBabelCommand);
	readOneString(file,"HDNNDir",&N2P2Dir);
	readOneString(file,"N2P2Dir",&N2P2Dir);
	readOneString(file,"tmModule",&tmModule);
	readOneString(file,"fireflyCommand",&fireflyCommand);
	readOneString(file,"cchemiCommand",&cchemiCommand);
	readOneReal(file,"tolEnergy",&tolEnergy);
	readOneReal(file,"tolDistance",&tolDistance);
	readOneBoolean(file,"ConfoOptMopac",&optMopac);
	readOneString(file,"ConfoOptMopacMethod",&optMopacMethod);
	readOneBoolean(file,"ConfoOptGaussian",&optGaussian);
	readOneString(file,"ConfoOptGaussianMethod",&optGaussianMethod);
	readOneBoolean(file,"ConfoOptFireFly",&optFireFly);
	readOneString(file,"ConfoOptFireFlyMethod",&optFireFlyMethod);
	readOneBoolean(file,"ConfoOptOrca",&optOrca);
	readOneString(file,"ConfoOptOrcaMethod",&optOrcaMethod);
	readOneBoolean(file,"ConfoOptOpenBabel",&optOpenBabel);
	readOneBoolean(file,"ConfoOptGeneric",&optGeneric);
	readOneString(file,"ConfoOptGenericMethod",&optGenericMethod);
	readOneBoolean(file,"ConfoOptN2P2",&optN2P2);
	readOneBoolean(file,"ConfoOptTM",&optTM);

	if(readOneInt(file,"Constraints",&cnt)) constraints = cnt;

	readOneString(file,"mopacKeywordsPost",&mopacKeywordsPost);
	readOneString(file,"gaussianKeywordsPost",&gaussianKeywordsPost);
	readOneString(file,"fireflyKeywordsPost",&fireflyKeywordsPost);
	readOneString(file,"orcaKeywordsPost",&orcaKeywordsPost);
	readOneString(file,"cchemiKeywordsPost",&cchemiKeywordsPost);
	readOneString(file,"OptimizerType",&optimizerType);
	if(!optimizerType) optimizerType = strdup("External");
	if(!strstr(optimizerType,"External") && (optMopac||optFireFly))
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
	if(!readOneString(file,"Model",&model)) model = strdup("MOPAC");
	uppercase(model);
	if(!readOneString(file,"QMKeys",&QMKeys)) 
	{
		if(!strcmp(model,"MOPAC")) QMKeys = strdup("PM6-DH+");
		else if(!strcmp(model,"ORCA")) QMKeys = strdup("AM1");
		else if(!strcmp(model,"GAUSSIAN")) QMKeys = strdup("AM1");
		else if(!strcmp(model,"FIREFLY")) QMKeys = strdup("AM1");
		else if(!strcmp(model,"OPENBABEL")) QMKeys = strdup("mmff94");
		else QMKeys = strdup("SCC-DFTB");
	}
	if(!strcmp(model,"OPENBABEL")) mol.klass->buildMMTypes(&mol, file);

	/* fileName for geometries */
	{
		char* suff = getSuffixNameFile(inputFileName);
		dirName = strdup(getenv("PWD"));
		fileNameGeom = strdup_printf("%s%s",suff, "Geoms.gab");
		fileNameTraj = strdup_printf("%s%s",suff, "Traj.gab");
		fileNameProp = strdup_printf("%s%s",suff, "Prop.txt");
		free(suff);
	}


	printf("Model = %s\n",model);
	printf("QMKeys = %s\n",QMKeys);
	if(!strcmp(model,"MOPAC")) qmModel = createMopacModel(&mol, QMKeys, dirName, mopacCommand, constraints, stdout);
	else  if(!strcmp(model,"FIREFLY")) qmModel = createFireFlyModel(&mol, QMKeys, dirName, fireflyCommand, constraints, stdout);
	else if(!strcmp(model,"ORCA")) qmModel = createOrcaModel(&mol, QMKeys, dirName, orcaCommand, constraints, stdout);
	else if(!strcmp(model,"OPENBABEL")) qmModel = createOpenBabelModel(&mol, QMKeys, dirName, openBabelCommand, constraints, stdout);
	else if(!strcmp(model,"GAUSSIAN")) qmModel = createGaussianModel(&mol, QMKeys, dirName, gaussianCommand, constraints, stdout);
	else if(!strcmp(model,"N2P2") || !strcmp(model,"HDNN")) qmModel = createN2P2Model(&mol, N2P2Dir, constraints, stdout);
	else if(!strcmp(model,"TM") || !strcmp(model,"TensorModel")) qmModel = createTMModel(&mol, tmModule, constraints, logfile);
	else qmModel = createGenericModel(&mol, QMKeys, dirName, genericCommand, constraints, stdout);
	setH4Correction(file,&qmModel);
	setSRBCorrection(file,&qmModel);
	readOneBoolean(file,"addD3Correction",&qmModel.addD3Correction);
	checkWallCorrection(file, &qmModel);

	geometries = runQuantumMechanicsMDConfo(&molecularDynamics, &qmModel,
		updateFrequency, heatTime, equiTime, runTime, heatTemp, equiTemp, runTemp, stepSize, 
		integrator, thermostat, friction, omegaMax, Nf, collide, qNH, numberOfGeometries, fileNameTraj, fileNameProp);
	if(geometries && numberOfGeometries>0)
	{
		int i;
		energies = malloc(numberOfGeometries*sizeof(double));
		for(i=0;i<numberOfGeometries;i++)
			energies[i] = geometries[i]->molecule.potentialEnergy;
	}
	optAndSortGeometries(&qmModel, inputFileName, &geometries, &numberOfGeometries, mopacCommand, gaussianCommand, N2P2Dir, tmModule, fireflyCommand, cchemiCommand,
	orcaCommand, openBabelCommand, genericCommand, optMopacMethod, optGaussianMethod, optFireFlyMethod, optOrcaMethod, optGenericMethod, &energies,
	optMopac, optGaussian, optGeneric, optFireFly, optOrca, optOpenBabel, optN2P2, optTM, tolEnergy, tolDistance, fileNameGeom, removeSimilarInertia,
	removeFragmented, removeSmallDistance, removeSimilarBonds, inertiaTol, sTol, distMaxTol, logfile, message);

	if(numberOfGeometries>0 && geometries )
	{
		char* fileNamePrefix = getSuffixNameFile(fileNameGeom);
		createPostProcessingFiles(numberOfGeometries, geometries, energies,fileNamePrefix, mopacKeywordsPost, gaussianKeywordsPost, fireflyKeywordsPost, orcaKeywordsPost, cchemiKeywordsPost, message, mopacCommand, gaussianCommand, fireflyCommand, orcaCommand, cchemiCommand);
		if(fileNamePrefix) free(fileNamePrefix);
/*
		if(mopacKeywordsPost)
		{
			char* fileNamePrefix = getSuffixNameFile(fileNameGeom);
			createMopacFiles(numberOfGeometries, geometries, energies, fileNamePrefix, mopacKeywordsPost,mopacCommand);
			strcat(message,fileNamePrefix);
			strcat(message,("_*.mop\n\tFiles for a post processing by Mopac\n\n"));
			if(fileNamePrefix) free(fileNamePrefix);
		}
		if(gaussianKeywordsPost)
		{
			char* fileNamePrefix = getSuffixNameFile(fileNameGeom);
			createGaussianFiles(numberOfGeometries, geometries, energies, fileNamePrefix, gaussianKeywordsPost, gaussianCommand);
			strcat(message,fileNamePrefix);
			strcat(message,("_*.com\n\tFiles for a post processing by Gaussian\n\n"));
			if(fileNamePrefix) free(fileNamePrefix);
		}
		if(fireflyKeywordsPost)
		{
			char* fileNamePrefix = getSuffixNameFile(fileNameGeom);
			createFireFlyFiles(numberOfGeometries, geometries, energies, fileNamePrefix, fireflyKeywordsPost,fireflyCommand);
			strcat(message,fileNamePrefix);
			strcat(message,("FF_*.inp\n\tFiles for a post processing by FireFly\n\n"));
			if(fileNamePrefix) free(fileNamePrefix);
		}
*/
	}
	if(geometries)
	{
		for(i=0;i<numberOfGeometries;i++)
			if(geometries[i]) geometries[i]->klass->free(geometries[i]);
		free(geometries);
	}
	if(energies) free(energies);
	if(strlen(message)>20) printf("%s",message);
	if(fileNameGeom)free(fileNameGeom);
	fclose(file);
	qmModel.klass->free(&qmModel);

}
/*****************************************************************************/
static void quantumMechanicsMD(char* inputFileName)
{
	QuantumMechanicsModel qmModel; 
	QuantumMechanicsMD molecularDynamics;
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
	MDThermostatType thermostat = NONE;
	char* dirName = NULL;
	Constraints constraints = NOCONSTRAINTS;
	Molecule mol = *(readMolecule(inputFileName,TRUE));
	char* mopacCommand = strdup("/opt/mopac/MOPAC2012");
	char* fireflyCommand = strdup("firefly");
	char* orcaCommand = strdup("orca");
	char* N2P2Dir=strdup(getenv("PWD"));
	char* tmModule=strdup_printf("%s/%s",getenv("PWD"),"CChemITMModule");
	char* openBabelCommand = strdup("obgradient");
	char* gaussianCommand = strdup("g09");
	char* genericCommand = strdup("runGeneric");
	FILE* file = fopen(inputFileName,"rb");
	char* model = NULL;
	char* QMKeys = NULL;
	int cnt;
	FILE* logfile = stdout;

	setMDOptions(file, &updateFrequency, 
		&heatTime, &equiTime, &runTime, &coolTime,
		&heatTemp, &runTemp, &equiTemp, &coolTemp, &stepSize, 
		&integrator, &thermostat, &friction, &omegaMax, &Nf, &collide,&qNH);
	if(!readOneString(file,"Model",&model)) model = strdup("MOPAC");
	uppercase(model);
	if(!readOneString(file,"QMKeys",&QMKeys)) 
	{
		if(!strcmp(model,"MOPAC")) QMKeys = strdup("PM6-DH+");
		else if(!strcmp(model,"ORCA")) QMKeys = strdup("AM1");
		else if(!strcmp(model,"GAUSSIAN")) QMKeys = strdup("AM1");
		else if(!strcmp(model,"FIREFLY")) QMKeys = strdup("AM1");
		else if(!strcmp(model,"OPENBABEL")) QMKeys = strdup("mmff94");
		else QMKeys = strdup("SCC-DFTB");
	}
	if(!strcmp(model,"OPENBABEL")) mol.klass->buildMMTypes(&mol, file);

	readOneString(file,"mopacCommand",&mopacCommand);
	readOneString(file,"fireflyCommand",&fireflyCommand);
	readOneString(file,"orcaCommand",&orcaCommand);
	readOneString(file,"openBabelCommand",&openBabelCommand);
	readOneString(file,"HDNNDir",&N2P2Dir);
	readOneString(file,"N2P2Dir",&N2P2Dir);
	readOneString(file,"tmModule",&tmModule);
	readOneString(file,"gaussianCommand",&gaussianCommand);
	readOneString(file,"genericCommand",&genericCommand);

	if(readOneInt(file,"Constraints",&cnt)) constraints = cnt;


	{
		char* suff = getSuffixNameFile(inputFileName);
		dirName = strdup(getenv("PWD"));
		fileNameTraj = strdup_printf("%s%s",suff, "Traj.gab");
		fileNameProp = strdup_printf("%s%s",suff, "Prop.txt");
		free(suff);
	}

	if(!strcmp(model,"MOPAC")) qmModel = createMopacModel(&mol, QMKeys, dirName, mopacCommand, constraints,logfile);
	else if(!strcmp(model,"FIREFLY")) qmModel = createFireFlyModel(&mol, QMKeys,dirName, fireflyCommand, constraints,logfile);
	else if(!strcmp(model,"ORCA")) qmModel = createOrcaModel(&mol, QMKeys, dirName, orcaCommand, constraints, logfile);
	else if(!strcmp(model,"OPENBABEL")) qmModel = createOpenBabelModel(&mol, QMKeys, dirName, openBabelCommand, constraints, logfile);
	else if(!strcmp(model,"GAUSSIAN")) qmModel = createGaussianModel(&mol, QMKeys,dirName, gaussianCommand, constraints,logfile);
	else if(!strcmp(model,"N2P2") || !strcmp(model,"HDNN")) qmModel = createN2P2Model(&mol, N2P2Dir, constraints, logfile);
	else if(!strcmp(model,"TM") || !strcmp(model,"TensorModel")) qmModel = createTMModel(&mol, tmModule, constraints, logfile);
	else qmModel = createGenericModel(&mol, QMKeys,dirName, genericCommand, constraints,logfile);
	setSRBCorrection(file,&qmModel);
	readOneBoolean(file,"addD3Correction",&qmModel.addD3Correction);
	checkWallCorrection(file, &qmModel);

	runQuantumMechanicsMD(&molecularDynamics, &qmModel,
		updateFrequency, heatTime, equiTime, runTime, coolTime, heatTemp, equiTemp, runTemp, coolTemp, stepSize, 
		integrator, thermostat, friction, omegaMax, Nf, collide, qNH, fileNameTraj, fileNameProp);

	qmModel.klass->free(&qmModel);
	free(dirName);
	fclose(file);
}
/***********************************************************************************************************************/
static void quantumMechanicsRDConfo(char* inputFileName)
{
	QuantumMechanicsModel qmModel; 
	QuantumMechanicsModel** geometries = NULL;
	int i;
	char message[BSIZE]="Created files :\n";
	char* dirName = NULL;
	Molecule mol = *(readMolecule(inputFileName,TRUE));
	char* mopacCommand = strdup("mopac");
	char* gaussianCommand=strdup("g09"); 
	char* N2P2Dir=strdup(getenv("PWD"));
	char* tmModule=strdup_printf("%s/%s",getenv("PWD"),"CChemITMModule");
	char* fireflyCommand=strdup("firefly");
	char* cchemiCommand=strdup("cchemi");
	char* orcaCommand = strdup("orca");
	char* openBabelCommand = strdup("obgradient");
	char* genericCommand=strdup("runGeneric");
	FILE* file = fopen(inputFileName,"rb");
	char* mopacKeywordsPost = NULL;
	char* gaussianKeywordsPost = NULL;
	char* fireflyKeywordsPost = NULL;
	char* orcaKeywordsPost = NULL;
	char* cchemiKeywordsPost = NULL;
	char* optMopacMethod=strdup("PM6");
	char* optGaussianMethod=strdup("AM1");
	char* optFireFlyMethod=strdup("AM1");
	char* optOrcaMethod=strdup("AM1");
	char* optGenericMethod=strdup("runGeneric");
	char* model = NULL;
	char* QMKeys = NULL;
	int cnt;
	char* optimizerType = NULL;

	double* energies = NULL;
	boolean optMopac = FALSE;
	boolean optGaussian = FALSE;
	boolean optGeneric = FALSE;
	boolean optFireFly = FALSE;
	boolean optOrca = FALSE;
	boolean optOpenBabel = FALSE;
	boolean optN2P2 = FALSE;
	boolean optTM = FALSE;
	double tolEnergy = 0.01;
	double tolDistance = 0.01;
	int numberOfGeometries = 2;
	char* fileNameGeom = NULL;
	Constraints constraints = NOCONSTRAINTS;
	boolean chain=FALSE;
	boolean fragments=FALSE;
	boolean saveFirstGeom=FALSE;
	boolean removeSimilarInertia = FALSE;
	boolean removeFragmented = FALSE;
	boolean removeSmallDistance = FALSE;
	FILE* logfile = stdout;
	int nTimesGeoms=1;
	double inertiaTol = 0.04; // recomended byjun Zhao et al (2016)
       	//Comprehensive genetic algorithm for ab initio global optimisation of clusters, Molecular Simulation,
	// 42:10, 809-819, DOI: 10.1080/08927022.2015.1121386	     

	boolean removeSimilarBonds = FALSE;
	double sTol=0.02;
	double distMaxTol=0.7;
	// sTol = 0.02 , distMaxTol = 0.7 Ang, recommanded in Jorgensen et al JCTC, 2017
//Mathias S. Jørgensen , Michael N. Groves, and Bjørk Hammer
//J. Chem. Theory Comput., 2017, 13 (3), pp 1486–1493
//DOI: 10.1021/acs.jctc.6b01119

	readOneString(file,"mopacCommand",&mopacCommand);
	readOneString(file,"gaussianCommand",&gaussianCommand);
	readOneString(file,"genericCommand",&genericCommand);
	readOneString(file,"orcaCommand",&orcaCommand);
	readOneString(file,"openBabelCommand",&openBabelCommand);
	readOneString(file,"HDNNDir",&N2P2Dir);
	readOneString(file,"N2P2Dir",&N2P2Dir);
	readOneString(file,"tmModule",&tmModule);
	readOneString(file,"fireflyCommand",&fireflyCommand);
	readOneString(file,"cchemiCommand",&cchemiCommand);
	readOneReal(file,"tolEnergy",&tolEnergy);
	readOneReal(file,"tolDistance",&tolDistance);
	readOneBoolean(file,"ConfoOptMopac",&optMopac);
	readOneString(file,"ConfoOptMopacMethod",&optMopacMethod);
	readOneBoolean(file,"ConfoOptGaussian",&optGaussian);
	readOneString(file,"ConfoOptGaussianMethod",&optGaussianMethod);
	readOneBoolean(file,"ConfoOptFireFly",&optFireFly);
	readOneString(file,"ConfoOptFireFlyMethod",&optFireFlyMethod);
	readOneBoolean(file,"ConfoOptOrca",&optOrca);
	readOneString(file,"ConfoOptOrcaMethod",&optOrcaMethod);
	readOneBoolean(file,"ConfoOptOpenBabel",&optOpenBabel);
	readOneBoolean(file,"ConfoOptGeneric",&optGeneric);
	readOneString(file,"ConfoOptGenericMethod",&optGenericMethod);
	readOneBoolean(file,"ConfoOptN2P2",&optN2P2);
	readOneBoolean(file,"ConfoOptTM",&optTM);
	readOneBoolean(file,"RDChain",&chain);
	readOneBoolean(file,"RDFragments",&fragments);
	readOneBoolean(file,"RDSaveFirstGeom",&saveFirstGeom);
	readOneInt(file,"nTimesGeoms",&nTimesGeoms);

	if(readOneInt(file,"Constraints",&cnt)) constraints = cnt;

	readOneString(file,"mopacKeywordsPost",&mopacKeywordsPost);
	readOneString(file,"gaussianKeywordsPost",&gaussianKeywordsPost);
	readOneString(file,"fireflyKeywordsPost",&fireflyKeywordsPost);
	readOneString(file,"orcaKeywordsPost",&orcaKeywordsPost);
	readOneString(file,"cchemiKeywordsPost",&cchemiKeywordsPost);
	readOneString(file,"OptimizerType",&optimizerType);
	if(!optimizerType) optimizerType = strdup("External");
	if(!strstr(optimizerType,"External") && (optMopac||optFireFly))
	{
		printf("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n");
		printf("Sorry, The optimization using the internal optimizer with a RD conformational search is not yes implemented in this software\n");
		printf("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n");
		exit(1);
	}
	/* number for geometries */
	{
		numberOfGeometries = 10;
		readOneInt(file,"numberOfGeometries",&numberOfGeometries);
		if(numberOfGeometries<2) numberOfGeometries = 2;
	}
	if(!readOneString(file,"Model",&model)) model = strdup("MOPAC");
	uppercase(model);
	if(!readOneString(file,"QMKeys",&QMKeys)) 
	{
		if(!strcmp(model,"MOPAC")) QMKeys = strdup("PM6-DH+");
		else if(!strcmp(model,"ORCA")) QMKeys = strdup("AM1");
		else if(!strcmp(model,"GAUSSIAN")) QMKeys = strdup("AM1");
		else if(!strcmp(model,"FIREFLY")) QMKeys = strdup("AM1");
		else if(!strcmp(model,"OPENBABEL")) QMKeys = strdup("mmff94");
		else QMKeys = strdup("SCC-DFTB");
	}
	if(!strcmp(model,"OPENBABEL")) mol.klass->buildMMTypes(&mol, file);

	/* fileName for geometries */
	{
		char* suff = getSuffixNameFile(inputFileName);
		dirName = strdup(getenv("PWD"));
		fileNameGeom = strdup_printf("%s%s",suff, "Geoms.gab");
		free(suff);
	}
	if(readOneBoolean(file,"removeSimilarInertia",&removeSimilarInertia) && removeSimilarInertia) readOneReal(file,"InertiaTol",&inertiaTol);
	readOneBoolean(file,"removeFragmented",&removeFragmented);
	readOneBoolean(file,"removeDissociated",&removeFragmented);
	readOneBoolean(file,"removeSmallDistance",&removeSmallDistance);
	if(readOneBoolean(file,"removeSimilarBonds",&removeSimilarBonds) && removeSimilarBonds)
	{
		readOneReal(file,"sTol",&sTol);
		readOneReal(file,"distMaxTol",&distMaxTol);
	}


	printf("Model = %s\n",model);
	printf("QMKeys = %s\n",QMKeys);
	if(!strcmp(model,"MOPAC")) qmModel = createMopacModel(&mol, QMKeys, dirName, mopacCommand, constraints, logfile);
	else  if(!strcmp(model,"FIREFLY")) qmModel = createFireFlyModel(&mol, QMKeys, dirName, fireflyCommand, constraints, logfile);
	else if(!strcmp(model,"ORCA")) qmModel = createOrcaModel(&mol, QMKeys, dirName, orcaCommand, constraints, logfile);
	else if(!strcmp(model,"OPENBABEL")) qmModel = createOpenBabelModel(&mol, QMKeys, dirName, openBabelCommand, constraints, logfile);
	else if(!strcmp(model,"GAUSSIAN")) qmModel = createGaussianModel(&mol, QMKeys, dirName, gaussianCommand, constraints, logfile);
	else if(!strcmp(model,"N2P2") || !strcmp(model,"HDNN")) qmModel = createN2P2Model(&mol, N2P2Dir, constraints, logfile);
	else if(!strcmp(model,"TM") || !strcmp(model,"TensorModel")) qmModel = createTMModel(&mol, tmModule, constraints, logfile);
	else qmModel = createGenericModel(&mol, QMKeys, dirName, genericCommand, constraints, logfile);

	setH4Correction(file,&qmModel);
	setSRBCorrection(file,&qmModel);
	readOneBoolean(file,"addD3Correction",&qmModel.addD3Correction);
	checkWallCorrection(file, &qmModel);


	int nOld = numberOfGeometries*nTimesGeoms;
	if(nTimesGeoms>1) numberOfGeometries = nOld;
	if(fragments)
	geometries = qmModel.klass->getQuantumMechanicsRDFConfo(&qmModel, numberOfGeometries, saveFirstGeom);
	else
	geometries = qmModel.klass->getQuantumMechanicsRDConfo(&qmModel, numberOfGeometries, chain, saveFirstGeom);


	if(removeSimilarInertia) qmModel.klass->removeSimilarInertiaGeometries(geometries, &numberOfGeometries, energies,logfile,inertiaTol);
	if(removeFragmented) qmModel.klass->removeFragmentedMolecules(geometries, &numberOfGeometries, energies, logfile);
	if(removeSmallDistance) qmModel.klass->removeSmallDistanceMolecules(geometries, &numberOfGeometries, energies, logfile);
	if(removeSimilarBonds) qmModel.klass->removeSimilarBondsGeometries(geometries, &numberOfGeometries, energies,logfile,sTol, distMaxTol);

	if(nTimesGeoms>1) qmModel.klass->cutByInertia(geometries, &numberOfGeometries, energies,nOld/nTimesGeoms,logfile);

	if(geometries && numberOfGeometries>0)
	{
		int i;
		energies = malloc(numberOfGeometries*sizeof(double));
		for(i=0;i<numberOfGeometries;i++)
			energies[i] = geometries[i]->molecule.potentialEnergy;
	}

	optAndSortGeometries(&qmModel, inputFileName, &geometries, &numberOfGeometries, mopacCommand, gaussianCommand, N2P2Dir, tmModule, fireflyCommand, cchemiCommand,
	orcaCommand, openBabelCommand, genericCommand, optMopacMethod, optGaussianMethod, optFireFlyMethod, optOrcaMethod, optGenericMethod, &energies,
	optMopac, optGaussian, optGeneric, optFireFly, optOrca, optOpenBabel, optN2P2, optTM, tolEnergy, tolDistance, fileNameGeom, removeSimilarInertia,
	removeFragmented, removeSmallDistance, removeSimilarBonds, inertiaTol, sTol, distMaxTol, logfile, message);

	if(numberOfGeometries>0 && geometries )
	{
		char* fileNamePrefix = getSuffixNameFile(fileNameGeom);
		createPostProcessingFiles(numberOfGeometries, geometries, energies,fileNamePrefix, mopacKeywordsPost, gaussianKeywordsPost, fireflyKeywordsPost, orcaKeywordsPost, cchemiKeywordsPost, message, mopacCommand, gaussianCommand, fireflyCommand, orcaCommand, cchemiCommand);
		if(fileNamePrefix) free(fileNamePrefix);
	}
	if(geometries)
	{
		for(i=0;i<numberOfGeometries;i++)
			if(geometries[i]) geometries[i]->klass->free(geometries[i]);
		free(geometries);
	}
	qmModel.klass->free(&qmModel);
	if(energies) free(energies);
	if(strlen(message)>20) printf("%s",message);
	if(fileNameGeom)free(fileNameGeom);
	fclose(file);

}
/***********************************************************************************************************************/
static void quantumMechanicsGAConfo(char* inputFileName)
{
	GeneticAlgorithm* ga; 
	ga = newGeneticAlgorithm(inputFileName);
	if(ga) ga->klass->run(ga);
}
/********************************************************************************/
void quantumMechanicsMolecularDynamicsDlg(char* inputFileName)
{
	quantumMechanicsMD(inputFileName);
}
/***********************************************************************/
void quantumMechanicsMolecularDynamicsConfoDlg(char* inputFileName)
{
	quantumMechanicsMDConfo(inputFileName);
}
/***********************************************************************/
void quantumMechanicsREMDConfoDlg(char* inputFileName)
{
	quantumMechanicsREMDConfo(inputFileName);
}
/***********************************************************************/
void quantumMechanicsRandomConfoDlg(char* inputFileName)
{
	quantumMechanicsRDConfo(inputFileName);
}
/***********************************************************************/
void quantumMechanicsGAConfoDlg(char* inputFileName)
{
	quantumMechanicsGAConfo(inputFileName);
}
/***********************************************************************/
void quantumMechanicsRemoveSimilarConfoDlg(char* inputFileName)
{
	quantumMechanicsRemoveSimilarConfo(inputFileName);
}
/*****************************************************************************/
static char*  setOptOptions(FILE* file, ConjugateGradientQMOptions* conjugateGradientQMOptions, QuasiNewtonQM* quasiNewtonQM)
{
	char* optimizerType = strdup("External");
	readOneString(file,"OptimizerType",&optimizerType);
	if(strstr(optimizerType,"Gradient")) setCGQMOptions(file, conjugateGradientQMOptions);
	else if(strstr(optimizerType,"Quasi")) setQNQMOptions(file, quasiNewtonQM);
	else if(strstr(optimizerType,"Steep")) setCGQMOptions(file, conjugateGradientQMOptions);
	return optimizerType;
}
/*****************************************************************************/
void quantumMechanicsMinimizeExternalDlg(char* inputFileName)
{
	Molecule mol = *(readMolecule(inputFileName, TRUE));
	char* mopacCommand = strdup("/opt/mopac/MOPAC2012");
	char* fireflyCommand = strdup("firefly");
	char* orcaCommand = strdup("orca");
	char* openBabelCommand = strdup("obgradient");
	char* gaussianCommand = strdup("g09");
	char* genericCommand = strdup("runGeneric");
	char* N2P2Dir=strdup(getenv("PWD"));
	FILE* file = NULL;
	char* runType = NULL;
	char* model = NULL;
	char* QMKeys = NULL;
	char* fileName = strdup_printf("%sOpt.gab",getSuffixNameFile(inputFileName));
	if(!inputFileName)
	{
		printf("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n");
		printf("Usage : cchemi inputFileName.inp\n");
		printf("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n");
		exit(1);
	}
	file = fopen(inputFileName, "r");
	if(!file)
	{
		printf("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n");
		printf("Sorry, I cannot open the input file : %s\n",inputFileName);
		printf("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n");
		exit(1);
	}
	if(!readOneString(file,"RunType",&runType)) runType = strdup("ENERGY");
	if(!readOneString(file,"Model",&model)) model = strdup("MOPAC");
	uppercase(model);
	if(!readOneString(file,"QMKeys",&QMKeys)) 
	{
		if(!strcmp(model,"MOPAC")) QMKeys = strdup("PM6-DH+");
		else if(!strcmp(model,"ORCA")) QMKeys = strdup("AM1");
		else if(!strcmp(model,"GAUSSIAN")) QMKeys = strdup("AM1");
		else if(!strcmp(model,"FIREFLY")) QMKeys = strdup("AM1");
		else if(!strcmp(model,"OPENBABEL")) QMKeys = strdup("mmff94");
		else QMKeys = strdup("SCC-DFTB");
	}
	if(!strcmp(model,"OPENBABEL")) mol.klass->buildMMTypes(&mol, file);
	uppercase(runType);
	readOneString(file,"mopacCommand",&mopacCommand);
	readOneString(file,"orcaCommand",&orcaCommand);
	readOneString(file,"openBabelCommand",&openBabelCommand);
	readOneString(file,"HDNNDir",&N2P2Dir);
	readOneString(file,"N2P2Dir",&N2P2Dir);
	readOneString(file,"fireflyCommand",&fireflyCommand);
	readOneString(file,"gaussianCommand",&gaussianCommand);
	readOneString(file,"genericCommand",&genericCommand);
	if(!strcmp(model,"FIREFLY")) QMKeys = strdup("GBASIS=AM1");
	printf("QMKeys = %s\n",QMKeys);
	printf("Model = %s\n",model);
	if(!strcmp(model,"MOPAC"))
	{
		char* keys = strdup_printf("XYZ %s",QMKeys);
		printf("keys = %s\n",keys);
		printf("mopacCommand = %s\n",mopacCommand);
		runMopac(&mol, fileName, keys, mopacCommand);
		free(keys);
	}
	else if(!strcmp(model,"ORCA"))
	{
		char* keys = strdup_printf("%s Opt",QMKeys);
		runOrca(&mol, fileName, keys, orcaCommand);
		free(keys);
	}
	else if(!strcmp(model,"FIREFLY"))
	{
		char* keys = strdup_printf("RUNTYP=Optimize %s",QMKeys);
		printf("keys=%s\n",QMKeys);
		runFireFly(&mol, fileName, keys, fireflyCommand);
		free(keys);
	}
	else if(!strcmp(model,"GAUSSIAN"))
	{
		char* keys = strdup_printf("%s Opt",QMKeys);
		printf("keys=%s\n",QMKeys);
		runGaussian(&mol, fileName, keys, gaussianCommand);
		free(keys);
	}
	else if(!strcmp(model,"GENERIC"))
	{
		char* keys = strdup_printf("%s Opt",QMKeys);
		runGeneric(&mol, fileName, keys, genericCommand);
		free(keys);
	}
	else if(!strcmp(model,"OPENBABEL"))
	{
		char* keys = strdup_printf("%s Opt",QMKeys);
		runOpenBabel(&mol, fileName, keys, openBabelCommand);
		free(keys);
	}
	else if(!strcmp(model,"N2P2") || !strcmp(model,"HDNN"))
	{
		fprintf(stderr,"Optimization by a external minimizer for N2P2 not yet implemented\n");
	}
	else if(!strcmp(model,"TM") || !strcmp(model,"TensorMol"))
	{
		char* tmModule=strdup_printf("%s","CChemITMModule");
		int err = 0;
		readOneString(file,"tmModule",&tmModule);
		readOneString(file,"TensorMolModule",&tmModule);
		if(tmModule) err = runOptTensorMol(&mol, tmModule);
		if(err)
		{
			fprintf(stderr," You must compile cchemi with enable_python = 1 in CONFIG file\n");
		}
		else{
			char* fileNameOut = strdup_printf("%sOpt.gab",getSuffixNameFile(inputFileName));
			printf("Optimized geometry saved in %s file\n",fileNameOut);
			mol.klass->save(&mol,fileNameOut);
			free(fileNameOut);
		}
	}
	fclose(file);
	if(fileName) free(fileName);
}
/*****************************************************************************/
static boolean minimizeGeometry(QuantumMechanicsModel* qmModel, 
	char* optimizerType, SteepestDescentQM steepestDescentQM, 
	QuasiNewtonQM quasiNewtonQM, ConjugateGradientQM conjugateGradientQM,
	ConjugateGradientQMOptions conjugateGradientQMOptions
	)
{
	if(strstr(optimizerType,"Grad"))
	{
		//printf("Minimization by Conjugate Gradient method\n");
		conjugateGradientQM.logfile= stdout;
		runConjugateGradientQM(&conjugateGradientQM, qmModel, conjugateGradientQMOptions); 
	}
	else if(strstr(optimizerType,"Quasi"))
	{
		//printf("Minimization by QuasiNewton method\n");
		quasiNewtonQM.qmModel = qmModel; 
               	quasiNewtonQM.logfile = stdout;
                runQuasiNewtonQM(&quasiNewtonQM);
	}
	else
	{
		//printf("Minimization by steepest descent method\n");
		steepestDescentQM.logfile= stdout;
		runSteepestDescentQM(&steepestDescentQM, qmModel,
			       	conjugateGradientQMOptions.updateFrequency,
			       conjugateGradientQMOptions.maxIterations,
			       conjugateGradientQMOptions.gradientNorm,
			       conjugateGradientQMOptions.maxLines);
	}
	//printEnergyAndGradient(qmModel);
	return TRUE;
}
/*****************************************************************************/
static boolean miminizeGeometriesUsingInternalOptimizer(int numberOfGeometries, QuantumMechanicsModel** geometries, double* energies, char* inputFileName)
{
	SteepestDescentQM steepestDescentQM;
	ConjugateGradientQM conjugateGradientQM;
	QuasiNewtonQM quasiNewtonQM;
	char* optimizerType = NULL;
	ConjugateGradientQMOptions conjugateGradientQMOptions;
	FILE* file = NULL;
	int i;
	int nG = 0;
	int nM = 0;
	if(!geometries) return FALSE;

	file = fopen(inputFileName,"rb");
	optimizerType = setOptOptions(file, &conjugateGradientQMOptions, &quasiNewtonQM);
	fclose(file);

	for(i=0;i<numberOfGeometries;i++)
	{
		if(!geometries[i]) continue;
		nG++;
		fprintf(stdout,"Optimization of geometry # %d / %d\n",i+1, numberOfGeometries);
		if(minimizeGeometry(geometries[i], optimizerType, steepestDescentQM, quasiNewtonQM, conjugateGradientQM,conjugateGradientQMOptions))
		{
			energies[i] = geometries[i]->molecule.potentialEnergy;
			nM++;
		}
		else
		{
			fprintf(stderr,"============>Free geom number %d\n",i+1);
			geometries[i]->klass->free(geometries[i]);
			geometries[i] =NULL;
		}

	}
	return (nM>0);
}
/***********************************************************************/
void quantumMechanicsMinimizeInternalDlg(char* inputFileName)
{
	QuantumMechanicsModel qmModel; 
	SteepestDescentQM steepestDescentQM;
	ConjugateGradientQM conjugateGradientQM;
	QuasiNewtonQM quasiNewtonQM;
	char* optimizerType = NULL;
	ConjugateGradientQMOptions conjugateGradientQMOptions;
	Molecule* mol = readMolecule(inputFileName,TRUE);
	FILE* file = fopen(inputFileName,"rb");
	char* fileNameOut = strdup_printf("%sOpt.gab",getSuffixNameFile(inputFileName));
	char* model = NULL;
	char* QMKeys = strdup("AM1");
	char* dirName = NULL;
	char* mopacCommand = strdup("/opt/mopac/MOPAC2012");
	char* fireflyCommand = strdup("firefly");
	char* orcaCommand = strdup("orca");
	char* N2P2Dir=strdup(getenv("PWD"));
	char* tmModule=strdup_printf("%s/%s",getenv("PWD"),"CChemITMModule");
	char* openBabelCommand = strdup("obgradient");
	char* gaussianCommand = strdup("g09");
	char* genericCommand = strdup("runGeneric");
	Constraints constraints = NOCONSTRAINTS;
	FILE* logfile = stdout;


	optimizerType = setOptOptions(file, &conjugateGradientQMOptions, &quasiNewtonQM);

	if(!readOneString(file,"Model",&model)) model = strdup("MOPAC");
	uppercase(model);
	if(!readOneString(file,"QMKeys",&QMKeys)) 
	{
		if(!strcmp(model,"MOPAC")) QMKeys = strdup("PM6-DH+");
		else if(!strcmp(model,"ORCA")) QMKeys = strdup("AM1");
		else if(!strcmp(model,"GAUSSIAN")) QMKeys = strdup("AM1");
		else if(!strcmp(model,"FIREFLY")) QMKeys = strdup("AM1");
		else if(!strcmp(model,"OPENBABEL")) QMKeys = strdup("mmff94");
		else QMKeys = strdup("SCC-DFTB");
	}
	if(!strcmp(model,"OPENBABEL")) mol->klass->buildMMTypes(mol, file);

	dirName = strdup(getenv("PWD"));
	readOneString(file,"mopacCommand",&mopacCommand);
	readOneString(file,"orcaCommand",&orcaCommand);
	readOneString(file,"openBabelCommand",&openBabelCommand);
	readOneString(file,"HDNNDir",&N2P2Dir);
	readOneString(file,"N2P2Dir",&N2P2Dir);
	readOneString(file,"tmModule",&tmModule);
	readOneString(file,"fireflyCommand",&fireflyCommand);
	readOneString(file,"gaussianCommand",&gaussianCommand);
	readOneString(file,"genericCommand",&genericCommand);

	if(!strcmp(model,"MOPAC")) qmModel = createMopacModel(mol, QMKeys, dirName, mopacCommand, constraints, logfile);
	else if(!strcmp(model,"FIREFLY")) qmModel = createFireFlyModel(mol, QMKeys, dirName, fireflyCommand, constraints, logfile);
	else if(!strcmp(model,"ORCA")) qmModel = createOrcaModel(mol, QMKeys, dirName, orcaCommand, constraints, logfile);
	else if(!strcmp(model,"OPENBABEL")) qmModel = createOpenBabelModel(mol, QMKeys, dirName, openBabelCommand, constraints, logfile);
	else if(!strcmp(model,"GAUSSIAN")) qmModel = createGaussianModel(mol, QMKeys, dirName, gaussianCommand, constraints, logfile);
	else if(!strcmp(model,"N2P2") || !strcmp(model,"HDNN")) qmModel = createN2P2Model(mol, N2P2Dir, constraints, logfile);
	else if(!strcmp(model,"TM") || !strcmp(model,"TensorModel")) qmModel = createTMModel(mol, tmModule, constraints, logfile);
	else qmModel = createGenericModel(mol, QMKeys, dirName, genericCommand, constraints, logfile);
	setH4Correction(file,&qmModel);
	setSRBCorrection(file,&qmModel);
	readOneBoolean(file,"addD3Correction",&qmModel.addD3Correction);
	checkWallCorrection(file, &qmModel);



	if(strstr(optimizerType,"Grad"))
	{
		printf("Minimization by Conjugate Gradient method\n");
		conjugateGradientQM.logfile= stdout;
		runConjugateGradientQM(&conjugateGradientQM, &qmModel, conjugateGradientQMOptions); 
		printf("Optimized geometry saved in %s file\n",fileNameOut);
		qmModel.molecule.klass->save(&qmModel.molecule,fileNameOut);

		freeConjugateGradientQM(&conjugateGradientQM);
	}
	else if(strstr(optimizerType,"Quasi"))
	{
		printf("Minimization by QuasiNewton method\n");
		quasiNewtonQM.qmModel = &qmModel; 
               	quasiNewtonQM.logfile = stdout;
                runQuasiNewtonQM(&quasiNewtonQM);
		printf("Optimized geometry saved in %s file\n",fileNameOut);
		qmModel.molecule.klass->save(&qmModel.molecule,fileNameOut);
		freeQuasiNewtonQM(&quasiNewtonQM);
	}
	else
	{
		printf("Minimization by steepest descent method\n");
		steepestDescentQM.logfile= stdout;
		runSteepestDescentQM(&steepestDescentQM, &qmModel,
			       	conjugateGradientQMOptions.updateFrequency,
			       conjugateGradientQMOptions.maxIterations,
			       conjugateGradientQMOptions.gradientNorm,
			       conjugateGradientQMOptions.maxLines);
		printf("Optimized geometry saved in %s file\n",fileNameOut);
		qmModel.molecule.klass->save(&qmModel.molecule,fileNameOut);
		freeSteepestDescentQM(&steepestDescentQM);
	}
	printEnergyAndGradient(&qmModel);
	qmModel.klass->free(&qmModel);
	fclose(file);
}
/***********************************************************************/
void quantumMechanicsMinimizeDlg(char* inputFileName)
{
	
	char* optimizerType = NULL;
	char* fileName = NULL;
	boolean addH4Correction = FALSE;
	boolean addSRBCorrection = FALSE;
	boolean addD3Correction = FALSE;
	FILE* file = fopen(inputFileName,"rb");
	readOneBoolean(file,"addD3Correction",&addD3Correction);
	readOneString(file,"OptimizerType",&optimizerType);
	if(readOneString(file,"H4FileName",&fileName)) addH4Correction = TRUE;
	if(readOneString(file,"SRBFileName",&fileName)) addSRBCorrection = TRUE;
	fclose(file);
	if(!optimizerType) optimizerType = strdup("External");
	if(strstr(optimizerType,"External") && (addH4Correction || addD3Correction || addSRBCorrection))
	{
		printf("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n");
		printf("I cannot use a External optimizer with H4 or D3 Correction\n");
		printf("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n");
		exit(1);
	}
	if(strstr(optimizerType,"External")) printf("optimization using the Computational chemistry optimizer\n");
	else  printf("We use the internal optimizer\n");
	if(strstr(optimizerType,"External")) quantumMechanicsMinimizeExternalDlg(inputFileName);
	else quantumMechanicsMinimizeInternalDlg(inputFileName);
}
/*****************************************************************************/
void quantumMechanicsEnergyDlg(char* inputFileName)
{
	QuantumMechanicsModel qmModel;
	Molecule* mol = readMolecule(inputFileName,TRUE);
	FILE* file = fopen(inputFileName,"rb");
	char* fileNameOut = strdup_printf("%s.gab",getSuffixNameFile(inputFileName));
	char* model = NULL;
	char* QMKeys = NULL;
	char* dirName = NULL;
	char* mopacCommand = strdup("/opt/mopac/MOPAC2012");
	char* fireflyCommand = strdup("firefly");
	char* orcaCommand = strdup("orca");
	char* N2P2Dir=strdup(getenv("PWD"));
	char* tmModule=strdup_printf("%s/%s",getenv("PWD"),"CChemITMModule");
	char* openBabelCommand = strdup("orca");
	char* gaussianCommand = strdup("g09");
	char* genericCommand = strdup("runGeneric");
	Constraints constraints = NOCONSTRAINTS;
	FILE* logfile = stdout;


	if(!readOneString(file,"Model",&model)) model = strdup("MOPAC");
	uppercase(model);
	if(!readOneString(file,"QMKeys",&QMKeys)) 
	{
		if(!strcmp(model,"MOPAC")) QMKeys = strdup("PM6-DH+");
		else if(!strcmp(model,"ORCA")) QMKeys = strdup("AM1");
		else if(!strcmp(model,"GAUSSIAN")) QMKeys = strdup("AM1");
		else if(!strcmp(model,"FIREFLY")) QMKeys = strdup("AM1");
		else if(!strcmp(model,"OPENBABEL")) QMKeys = strdup("mmff94");
		else QMKeys = strdup("SCC-DFTB");
	}
	if(!strcmp(model,"OPENBABEL")) mol->klass->buildMMTypes(mol, file);

	dirName = strdup(getenv("PWD"));
	readOneString(file,"mopacCommand",&mopacCommand);
	readOneString(file,"orcaCommand",&orcaCommand);
	readOneString(file,"openBabelCommand",&openBabelCommand);
	readOneString(file,"HDNNDir",&N2P2Dir);
	readOneString(file,"N2P2Dir",&N2P2Dir);
	readOneString(file,"tmModule",&tmModule);
	readOneString(file,"fireflyCommand",&fireflyCommand);
	readOneString(file,"gaussianCommand",&gaussianCommand);
	readOneString(file,"genericCommand",&genericCommand);

	if(!strcmp(model,"MOPAC")) qmModel = createMopacModel(mol, QMKeys, dirName, mopacCommand, constraints, logfile);
	else if(!strcmp(model,"FIREFLY")) qmModel = createFireFlyModel(mol, QMKeys, dirName, fireflyCommand, constraints, logfile);
	else if(!strcmp(model,"ORCA")) qmModel = createOrcaModel(mol, QMKeys, dirName, orcaCommand, constraints, logfile);
	else if(!strcmp(model,"OPENBABEL")) qmModel = createOpenBabelModel(mol, QMKeys, dirName, openBabelCommand, constraints, logfile);
	else if(!strcmp(model,"GAUSSIAN")) qmModel = createGaussianModel(mol, QMKeys, dirName, gaussianCommand, constraints, logfile);
	else if(!strcmp(model,"N2P2") || !strcmp(model,"HDNN")) qmModel = createN2P2Model(mol, N2P2Dir, constraints, logfile);
	else if(!strcmp(model,"TM") || !strcmp(model,"TensorModel")) qmModel = createTMModel(mol, tmModule, constraints, logfile);
	else qmModel = createGenericModel(mol, QMKeys, dirName, genericCommand, constraints, logfile);
	setH4Correction(file,&qmModel);
	setSRBCorrection(file,&qmModel);
	readOneBoolean(file,"addD3Correction",&qmModel.addD3Correction);
	checkWallCorrection(file, &qmModel);

	printEnergy(&qmModel);
	printf("Geometry saved in %s file\n",fileNameOut);
	qmModel.molecule.klass->save(&qmModel.molecule,fileNameOut);
	if(qmModel.H4Parameters)
	{
		int i,j;
		printf("H4 Correction\n");
		printf("Gradient(analytic) in kcal/Ang\n");
		for (  i = 0; i < qmModel.molecule.nAtoms; i++ )
			for(j=0;j<3;j++) qmModel.molecule.atoms[i].gradient[j] = 0;
		qmModel.molecule.potentialEnergy = getH4Correction(&qmModel.molecule, qmModel.H4Parameters, TRUE);
		for (  i = 0; i < qmModel.molecule.nAtoms; i++ )
		{
			printf("%-6s ",qmModel.molecule.atoms[i].prop.symbol);
			for(j=0;j<3;j++) printf("%14.8f ",qmModel.molecule.atoms[i].gradient[j]); 
			printf("\n");
		}

		printf("energy = %f(kcal/mol)\n",qmModel.molecule.potentialEnergy);
	}
	if(qmModel.SRBParameters)
	{
		int i,j;
		printf("SRB Correction\n");
		printf("Gradient(analytic) in kcal/Ang\n");
		for (  i = 0; i < qmModel.molecule.nAtoms; i++ )
			for(j=0;j<3;j++) qmModel.molecule.atoms[i].gradient[j] = 0;
		qmModel.molecule.potentialEnergy = getSRBCorrection(&qmModel.molecule, qmModel.SRBParameters, TRUE);
		for (  i = 0; i < qmModel.molecule.nAtoms; i++ )
		{
			printf("%-6s ",qmModel.molecule.atoms[i].prop.symbol);
			for(j=0;j<3;j++) printf("%14.8f ",qmModel.molecule.atoms[i].gradient[j]); 
			printf("\n");
		}

		printf("energy = %f(kcal/mol)\n",qmModel.molecule.potentialEnergy);
	}
	if(qmModel.addD3Correction)
	{
		int i,j;
		printf("Geometry\n");
		for (  i = 0; i < qmModel.molecule.nAtoms; i++ )
		{
			printf("%-6s ",qmModel.molecule.atoms[i].prop.symbol);
			for(j=0;j<3;j++) printf("%14.8f ",qmModel.molecule.atoms[i].coordinates[j]); 
			printf("\n");
		}
		printf("D3 Correction\n");
		printf("Gradient(numeric) in kcal/Ang\n");
		for (  i = 0; i < qmModel.molecule.nAtoms; i++ )
			for(j=0;j<3;j++) qmModel.molecule.atoms[i].gradient[j] = 0;
		qmModel.molecule.potentialEnergy = getD3Correction(&qmModel.molecule, qmModel.method, TRUE);
		for (  i = 0; i < qmModel.molecule.nAtoms; i++ )
		{
			printf("%-6s ",qmModel.molecule.atoms[i].prop.symbol);
			for(j=0;j<3;j++) printf("%14.8f ",qmModel.molecule.atoms[i].gradient[j]); 
			printf("\n");
		}

		printf("energy = %f(kcal/mol)\n",qmModel.molecule.potentialEnergy);
	}

	qmModel.klass->free(&qmModel);
	fclose(file);
}
/*****************************************************************************/
void quantumMechanicsGradientDlg(char* inputFileName)
{
	QuantumMechanicsModel qmModel;
	Molecule* mol = readMolecule(inputFileName,TRUE);
	FILE* file = fopen(inputFileName,"rb");
	char* fileNameOut = strdup_printf("%s.gab",getSuffixNameFile(inputFileName));
	char* model = NULL;
	char* QMKeys = NULL;
	char* dirName = NULL;
	char* mopacCommand = strdup("/opt/mopac/MOPAC2012");
	char* fireflyCommand = strdup("firefly");
	char* orcaCommand = strdup("orca");
	char* N2P2Dir=strdup(getenv("PWD"));
	char* tmModule=strdup_printf("%s/%s",getenv("PWD"),"CChemITMModule");
	char* openBabelCommand = strdup("orca");
	char* gaussianCommand = strdup("g09");
	char* genericCommand = strdup("runGeneric");
	Constraints constraints = NOCONSTRAINTS;
	FILE* logfile = stdout;


	if(!readOneString(file,"Model",&model)) model = strdup("MOPAC");
	uppercase(model);
	if(!readOneString(file,"QMKeys",&QMKeys)) 
	{
		if(!strcmp(model,"MOPAC")) QMKeys = strdup("PM6-DH+");
		else if(!strcmp(model,"ORCA")) QMKeys = strdup("AM1");
		else if(!strcmp(model,"GAUSSIAN")) QMKeys = strdup("AM1");
		else if(!strcmp(model,"FIREFLY")) QMKeys = strdup("AM1");
		else if(!strcmp(model,"OPENBABEL")) QMKeys = strdup("mmff94");
		else QMKeys = strdup("SCC-DFTB");
	}
	if(!strcmp(model,"OPENBABEL")) mol->klass->buildMMTypes(mol, file);

	dirName = strdup(getenv("PWD"));
	readOneString(file,"mopacCommand",&mopacCommand);
	readOneString(file,"orcaCommand",&orcaCommand);
	readOneString(file,"openBabelCommand",&openBabelCommand);
	readOneString(file,"HDNNDir",&N2P2Dir);
	readOneString(file,"N2P2Dir",&N2P2Dir);
	readOneString(file,"tmModule",&tmModule);
	readOneString(file,"fireflyCommand",&fireflyCommand);
	readOneString(file,"gaussianCommand",&gaussianCommand);
	readOneString(file,"genericCommand",&genericCommand);

	if(!strcmp(model,"MOPAC")) qmModel = createMopacModel(mol, QMKeys, dirName, mopacCommand, constraints, logfile);
	else if(!strcmp(model,"FIREFLY")) qmModel = createFireFlyModel(mol, QMKeys, dirName, fireflyCommand, constraints, logfile);
	else if(!strcmp(model,"ORCA")) qmModel = createOrcaModel(mol, QMKeys, dirName, orcaCommand, constraints, logfile);
	else if(!strcmp(model,"OPENBABEL")) qmModel = createOpenBabelModel(mol, QMKeys, dirName, openBabelCommand, constraints, logfile);
	else if(!strcmp(model,"GAUSSIAN")) qmModel = createGaussianModel(mol, QMKeys, dirName, gaussianCommand, constraints, logfile);
	else if(!strcmp(model,"N2P2") || !strcmp(model,"HDNN")) qmModel = createN2P2Model(mol, N2P2Dir, constraints, logfile);
	else if(!strcmp(model,"TM") || !strcmp(model,"TensorModel")) qmModel = createTMModel(mol, tmModule, constraints, logfile);
	else qmModel = createGenericModel(mol, QMKeys, dirName, genericCommand, constraints, logfile);
	setH4Correction(file,&qmModel);
	setSRBCorrection(file,&qmModel);
	readOneBoolean(file,"addD3Correction",&qmModel.addD3Correction);
	checkWallCorrection(file, &qmModel);

	printEnergyAndGradient(&qmModel);
	printf("Geometry saved in %s file\n",fileNameOut);
	qmModel.molecule.klass->save(&qmModel.molecule,fileNameOut);
	if(qmModel.H4Parameters)
	{
		int i,j;
		printf("H4 Correction\n");
		printf("Gradient(analytic) in kcal/Ang\n");
		for (  i = 0; i < qmModel.molecule.nAtoms; i++ )
			for(j=0;j<3;j++) qmModel.molecule.atoms[i].gradient[j] = 0;
		qmModel.molecule.potentialEnergy = getH4Correction(&qmModel.molecule, qmModel.H4Parameters, TRUE);
		for (  i = 0; i < qmModel.molecule.nAtoms; i++ )
		{
			printf("%-6s ",qmModel.molecule.atoms[i].prop.symbol);
			for(j=0;j<3;j++) printf("%14.8f ",qmModel.molecule.atoms[i].gradient[j]); 
			printf("\n");
		}

		printf("energy = %f(kcal/mol)\n",qmModel.molecule.potentialEnergy);
	}
	if(qmModel.SRBParameters)
	{
		int i,j;
		printf("SRB Correction\n");
		printf("Gradient(analytic) in kcal/Ang\n");
		for (  i = 0; i < qmModel.molecule.nAtoms; i++ )
			for(j=0;j<3;j++) qmModel.molecule.atoms[i].gradient[j] = 0;
		qmModel.molecule.potentialEnergy = getSRBCorrection(&qmModel.molecule, qmModel.SRBParameters, TRUE);
		for (  i = 0; i < qmModel.molecule.nAtoms; i++ )
		{
			printf("%-6s ",qmModel.molecule.atoms[i].prop.symbol);
			for(j=0;j<3;j++) printf("%14.8f ",qmModel.molecule.atoms[i].gradient[j]); 
			printf("\n");
		}

		printf("energy = %f(kcal/mol)\n",qmModel.molecule.potentialEnergy);
	}
	if(qmModel.addD3Correction)
	{
		int i,j;
		printf("Geometry\n");
		for (  i = 0; i < qmModel.molecule.nAtoms; i++ )
		{
			printf("%-6s ",qmModel.molecule.atoms[i].prop.symbol);
			for(j=0;j<3;j++) printf("%14.8f ",qmModel.molecule.atoms[i].coordinates[j]); 
			printf("\n");
		}
		printf("D3 Correction\n");
		printf("Gradient(numeric) in kcal/Ang\n");
		for (  i = 0; i < qmModel.molecule.nAtoms; i++ )
			for(j=0;j<3;j++) qmModel.molecule.atoms[i].gradient[j] = 0;
		qmModel.molecule.potentialEnergy = getD3Correction(&qmModel.molecule, qmModel.method, TRUE);
		for (  i = 0; i < qmModel.molecule.nAtoms; i++ )
		{
			printf("%-6s ",qmModel.molecule.atoms[i].prop.symbol);
			for(j=0;j<3;j++) printf("%14.8f ",qmModel.molecule.atoms[i].gradient[j]); 
			printf("\n");
		}

		printf("energy = %f(kcal/mol)\n",qmModel.molecule.potentialEnergy);
	}

	qmModel.klass->free(&qmModel);
	fclose(file);
}
/********************************************************************************/
void quantumMechanicsFrequenciesDlg(char* inputFileName)
{
	double* frequencies = NULL;
	double* reducedMasses = NULL;
	double* IRIntensities = NULL;
	double** modes = NULL;
	int nModes = 0;
	int i;

	QuantumMechanicsModel qmModel;
	Molecule* mol = readMolecule(inputFileName,TRUE);
	FILE* file = fopen(inputFileName,"rb");
	char* fileNameOut = strdup_printf("%sFreq.gab",getSuffixNameFile(inputFileName));
	char* model = NULL;
	char* QMKeys = NULL;
	char* dirName = NULL;
	char* mopacCommand = strdup("/opt/mopac/MOPAC2012");
	char* fireflyCommand = strdup("firefly");
	char* orcaCommand = strdup("orca");
	char* N2P2Dir=strdup(getenv("PWD"));
	char* tmModule=strdup_printf("%s/%s",getenv("PWD"),"CChemITMModule");
	char* openBabelCommand = strdup("orca");
	char* gaussianCommand = strdup("g09");
	char* genericCommand = strdup("runGeneric");
	Constraints constraints = NOCONSTRAINTS;
	FILE* logfile = stdout;
	double dx = -1.0;


	readOneReal(file,"dx",&dx);
	if(!readOneString(file,"Model",&model)) model = strdup("MOPAC");
	uppercase(model);
	if(!readOneString(file,"QMKeys",&QMKeys)) 
	{
		if(!strcmp(model,"MOPAC")) QMKeys = strdup("PM6-DH+");
		else if(!strcmp(model,"ORCA")) QMKeys = strdup("AM1");
		else if(!strcmp(model,"GAUSSIAN")) QMKeys = strdup("AM1");
		else if(!strcmp(model,"FIREFLY")) QMKeys = strdup("AM1");
		else if(!strcmp(model,"OPENBABEL")) QMKeys = strdup("mmff94");
		else QMKeys = strdup("SCC-DFTB");
	}
	if(!strcmp(model,"OPENBABEL")) mol->klass->buildMMTypes(mol, file);

	dirName = strdup(getenv("PWD"));
	readOneString(file,"mopacCommand",&mopacCommand);
	readOneString(file,"orcaCommand",&orcaCommand);
	readOneString(file,"openBabelCommand",&openBabelCommand);
	readOneString(file,"HDNNDir",&N2P2Dir);
	readOneString(file,"N2P2Dir",&N2P2Dir);
	readOneString(file,"tmModule",&tmModule);
	readOneString(file,"fireflyCommand",&fireflyCommand);
	readOneString(file,"gaussianCommand",&gaussianCommand);
	readOneString(file,"genericCommand",&genericCommand);

	if(!strcmp(model,"MOPAC")) qmModel = createMopacModel(mol, QMKeys, dirName, mopacCommand, constraints, logfile);
	else if(!strcmp(model,"FIREFLY")) qmModel = createFireFlyModel(mol, QMKeys, dirName, fireflyCommand, constraints, logfile);
	else if(!strcmp(model,"ORCA")) qmModel = createOrcaModel(mol, QMKeys, dirName, orcaCommand, constraints, logfile);
	else if(!strcmp(model,"OPENBABEL")) qmModel = createOpenBabelModel(mol, QMKeys, dirName, openBabelCommand, constraints, logfile);
	else if(!strcmp(model,"GAUSSIAN")) qmModel = createGaussianModel(mol, QMKeys, dirName, gaussianCommand, constraints, logfile);
	else if(!strcmp(model,"N2P2") || !strcmp(model,"HDNN")) qmModel = createN2P2Model(mol, N2P2Dir, constraints, logfile);
	else if(!strcmp(model,"TM") || !strcmp(model,"TensorModel")) qmModel = createTMModel(mol, tmModule, constraints, logfile);
	else qmModel = createGenericModel(mol, QMKeys, dirName, genericCommand, constraints, logfile);
	setH4Correction(file,&qmModel);
	setSRBCorrection(file,&qmModel);
	readOneBoolean(file,"addD3Correction",&qmModel.addD3Correction);
	checkWallCorrection(file, &qmModel);
	fclose(file);

	printf("Computing of IR spectra... Please wait\n");
	if(dx>0)
	{
		/* use numerical method */
		qmModel.dx = dx;
		nModes = qmModel.klass->computeQMFrequenciesNumeric(&qmModel, &frequencies, &modes, &reducedMasses, &IRIntensities);
	}
	else
	{
		/* use analytical method if available */
		qmModel.dx = 1e-3;
		nModes = qmModel.klass->computeQMFrequencies(&qmModel, &frequencies, &modes, &reducedMasses, &IRIntensities);
	}

	printEnergyAndGradient(&qmModel);
	printf("Frequencies and modes in the %s file\n",fileNameOut);
	mol->klass->saveFrequencies(mol, fileNameOut, nModes, frequencies, modes, reducedMasses, IRIntensities);
	addHarmonicVelocities(inputFileName, nModes, frequencies, modes, reducedMasses, IRIntensities);
	printHarmonicVelocities(inputFileName, nModes, frequencies, modes, reducedMasses);
	if(frequencies) free(frequencies);
	if(reducedMasses) free(reducedMasses);
	for(i=0;i<nModes;i++) free(modes[i]);
	if(modes) free(modes);
	qmModel.klass->free(&qmModel);
	free(fileNameOut);
}
/********************************************************************************/
void quantumMechanicsOptFrequenciesDlg(char* inputFileName)
{
	double* frequencies = NULL;
	double* reducedMasses = NULL;
	double* IRIntensities = NULL;
	double** modes = NULL;
	int nModes = 0;
	int i;

	QuantumMechanicsModel qmModel;
	Molecule* mol = readMolecule(inputFileName,TRUE);
	FILE* file = NULL;
	char* fileNameOut = strdup_printf("%sFreq.gab",getSuffixNameFile(inputFileName));
	char* fileNameOpt = strdup_printf("%sOpt.gab",getSuffixNameFile(inputFileName));
	char* model = NULL;
	char* QMKeys = NULL;
	char* dirName = NULL;
	char* mopacCommand = strdup("/opt/mopac/MOPAC2012");
	char* fireflyCommand = strdup("firefly");
	char* orcaCommand = strdup("orca");
	char* N2P2Dir=strdup(getenv("PWD"));
	char* tmModule=strdup_printf("%s/%s",getenv("PWD"),"CChemITMModule");
	char* openBabelCommand = strdup("obgradient");
	char* gaussianCommand = strdup("g09");
	char* genericCommand = strdup("runGeneric");
	Constraints constraints = NOCONSTRAINTS;
	FILE* logfile = stdout;
	double dx = -1.0;

	quantumMechanicsMinimizeDlg(inputFileName);
	mol = readMoleculeFromGabeditFile(fileNameOpt);
	/* mol->klass->print(mol, logfile);*/

	file = fopen(inputFileName,"rb");

	readOneReal(file,"dx",&dx);

	if(!readOneString(file,"Model",&model)) model = strdup("MOPAC");
	uppercase(model);
	if(!readOneString(file,"QMKeys",&QMKeys)) 
	{
		if(!strcmp(model,"MOPAC")) QMKeys = strdup("PM6-DH+");
		else if(!strcmp(model,"ORCA")) QMKeys = strdup("AM1");
		else if(!strcmp(model,"GAUSSIAN")) QMKeys = strdup("AM1");
		else if(!strcmp(model,"FIREFLY")) QMKeys = strdup("AM1");
		else if(!strcmp(model,"OPENBABEL")) QMKeys = strdup("mmff94");
		else QMKeys = strdup("SCC-DFTB");
	}
	if(!strcmp(model,"OPENBABEL")) mol->klass->buildMMTypes(mol, file);

	dirName = strdup(getenv("PWD"));
	readOneString(file,"mopacCommand",&mopacCommand);
	readOneString(file,"orcaCommand",&orcaCommand);
	readOneString(file,"openBabelCommand",&openBabelCommand);
	readOneString(file,"HDNNDir",&N2P2Dir);
	readOneString(file,"N2P2Dir",&N2P2Dir);
	readOneString(file,"tmModule",&tmModule);
	readOneString(file,"fireflyCommand",&fireflyCommand);
	readOneString(file,"gaussianCommand",&gaussianCommand);
	readOneString(file,"genericCommand",&genericCommand);

	if(!strcmp(model,"MOPAC")) qmModel = createMopacModel(mol, QMKeys, dirName, mopacCommand, constraints, logfile);
	else if(!strcmp(model,"FIREFLY")) qmModel = createFireFlyModel(mol, QMKeys, dirName, fireflyCommand, constraints, logfile);
	else if(!strcmp(model,"ORCA")) qmModel = createOrcaModel(mol, QMKeys, dirName, orcaCommand, constraints, logfile);
	else if(!strcmp(model,"OPENBABEL")) qmModel = createOpenBabelModel(mol, QMKeys, dirName, openBabelCommand, constraints, logfile);
	else if(!strcmp(model,"GAUSSIAN")) qmModel = createGaussianModel(mol, QMKeys, dirName, gaussianCommand, constraints, logfile);
	else if(!strcmp(model,"N2P2") || !strcmp(model,"HDNN")) qmModel = createN2P2Model(mol, N2P2Dir, constraints, logfile);
	else if(!strcmp(model,"TM") || !strcmp(model,"TensorModel")) qmModel = createTMModel(mol, tmModule, constraints, logfile);
	else qmModel = createGenericModel(mol, QMKeys, dirName, genericCommand, constraints, logfile);
	setH4Correction(file,&qmModel);
	setSRBCorrection(file,&qmModel);
	readOneBoolean(file,"addD3Correction",&qmModel.addD3Correction);
	checkWallCorrection(file, &qmModel);
	fclose(file);

	printf("Computing of IR spectra... Please wait\n");
	if(dx>0)
	{
		/* use numerical method */
		qmModel.dx = dx;
		nModes = qmModel.klass->computeQMFrequenciesNumeric(&qmModel, &frequencies, &modes, &reducedMasses, &IRIntensities);
	}
	else
	{
		/* use analytical method if available */
		qmModel.dx = 1e-3;
		nModes = qmModel.klass->computeQMFrequencies(&qmModel, &frequencies, &modes, &reducedMasses, &IRIntensities);
	}

	printEnergyAndGradient(&qmModel);
	printf("Frequencies and modes in the %s file\n",fileNameOut);
	mol->klass->saveFrequencies(mol, fileNameOut, nModes, frequencies, modes, reducedMasses, IRIntensities);
	addHarmonicVelocities(inputFileName, nModes, frequencies, modes, reducedMasses, IRIntensities);
	printHarmonicVelocities(inputFileName, nModes, frequencies, modes, reducedMasses);
	if(frequencies) free(frequencies);
	if(reducedMasses) free(reducedMasses);
	for(i=0;i<nModes;i++) free(modes[i]);
	if(modes) free(modes);
	qmModel.klass->free(&qmModel);
	free(fileNameOut);
	free(fileNameOpt);
}
/********************************************************************************/
