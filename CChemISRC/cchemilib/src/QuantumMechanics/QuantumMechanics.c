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

/* QuantumMechanics.c */

#ifndef OS_WIN32
#include <unistd.h>
#endif

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
#include "../EmpriricalCorrections/HydrogenBondCorrection.h"
#include "../EmpriricalCorrections/ShortRangeBasisSetCorrection.h"
#include "../EmpriricalCorrections/DispersionCorrection.h"
#include "../EmpriricalCorrections/WallCorrection.h"

static void calculateGradientMopac(QuantumMechanicsModel* qmModel);
static void calculateEnergyMopac(QuantumMechanicsModel* qmModel);
static void calculateGradientFireFly(QuantumMechanicsModel* qmModel);
static void calculateEnergyFireFly(QuantumMechanicsModel* qmModel);
static void calculateGradientOrca(QuantumMechanicsModel* qmModel);
static void calculateEnergyOrca(QuantumMechanicsModel* qmModel);
static void calculateGradientGaussian(QuantumMechanicsModel* qmModel);
static void calculateEnergyGaussian(QuantumMechanicsModel* qmModel);
static void calculateGradientGeneric(QuantumMechanicsModel* qmModel);
static void calculateEnergyGeneric(QuantumMechanicsModel* qmModel);
static void calculateGradientOpenBabel(QuantumMechanicsModel* qmModel);
static void calculateEnergyOpenBabel(QuantumMechanicsModel* qmModel);
static void calculateHessianN2P2(QuantumMechanicsModel* qmModel, double **F, double*** dmu);
static void calculateGradientN2P2(QuantumMechanicsModel* qmModel);
static void calculateEnergyN2P2(QuantumMechanicsModel* qmModel);

#ifdef ENABLE_PYTHON
static void calculateGradientTM(QuantumMechanicsModel* qmModel);
static void calculateEnergyTM(QuantumMechanicsModel* qmModel);
#endif

static boolean getDipoleGaussian(char* fileNameOut, double* dipole);

/*****************************************************************************/
void setSRBCorrection(FILE* file,QuantumMechanicsModel* qmModel) 
{
	char* fileName = NULL;

	if(readOneString(file,"SRBCorrection",&fileName) && fileName) 
	{
		char tmp[BSIZE];
		sprintf(tmp,"%s",fileName);
		uppercase(tmp);
		if(!strstr(tmp,"NONE"))
		{
			ShortRangeBasisSetCorrectionParameters parameters;
			if(!strstr(tmp,"DEFAULT")) setShortRangeBasisSetCorrectionParameters(&parameters, fileName, NULL);
			else setShortRangeBasisSetCorrectionParameters(&parameters, NULL, qmModel->method);
			qmModel->SRBParameters = malloc(sizeof(ShortRangeBasisSetCorrectionParameters));
			*(qmModel->SRBParameters) = parameters;
			return;
		}
	}
	qmModel->SRBParameters= NULL;
}
/**********************************************************************/
static void addSRBCorrection(QuantumMechanicsModel* qmModel,boolean addGradient)
{
	if(qmModel->SRBParameters)
		qmModel->molecule.potentialEnergy += getSRBCorrection(&qmModel->molecule, qmModel->SRBParameters, addGradient);
	//if(!qmModel->SRBParameters) { fprintf(stdout,"DEBUG SRB=NULL\n"); exit(1);}
}
/*****************************************************************************/
void setH4Correction(FILE* file,QuantumMechanicsModel* qmModel) 
{
	char* fileName = NULL;

	if(readOneString(file,"H4Correction",&fileName) && fileName) 
	{
		char tmp[BSIZE];
		sprintf(tmp,"%s",fileName);
		uppercase(tmp);
		if(!strstr(tmp,"NONE"))
		{
			HyhrogenBondCorrectionParameters parameters;
			if(!strstr(tmp,"DEFAULT")) setHydrogenBondCorrectionParameters(&parameters, fileName, NULL);
			else setHydrogenBondCorrectionParameters(&parameters, NULL, qmModel->method);
			qmModel->H4Parameters = malloc(sizeof(HyhrogenBondCorrectionParameters));
			*(qmModel->H4Parameters) = parameters;
			return;
		}
	}
	qmModel->H4Parameters= NULL;
}
/**********************************************************************/
static void addH4Correction(QuantumMechanicsModel* qmModel,boolean addGradient)
{
	if(qmModel->H4Parameters)
		qmModel->molecule.potentialEnergy += getH4Correction(&qmModel->molecule, qmModel->H4Parameters, addGradient);
}
/**********************************************************************/
static void addD3Correction(QuantumMechanicsModel* qmModel,boolean addGradient)
{
	if(qmModel->addD3Correction)
		qmModel->molecule.potentialEnergy += getD3Correction(&qmModel->molecule, qmModel->method, addGradient);
}
/**********************************************************************/
static void addWallCorrection(QuantumMechanicsModel* qmModel,boolean addGradient)
{
	if(qmModel->addWallCorrection)
	qmModel->molecule.potentialEnergy += getWallCorrection(&qmModel->molecule, addGradient);
}
/****************************************************************/
static void getMultiplicityName(int multiplicity, char* buffer)
{
//printf("mult = %d\n",multiplicity);
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
/*****************************************************************************/
static boolean getDipoleMopac(char* fileNameOut, double* dipole)
{
	FILE* file = NULL;
	char buffer[1024];
	char* pdest = NULL;

 	file = fopen(fileNameOut, "r");
	if(!file) return FALSE;
	 while(!feof(file))
	 {
		if(!fgets(buffer,BSIZE,file))break;
		if(
		strstr( buffer, "DIPOLE")&&
		strstr( buffer, "X")&&
		strstr( buffer, "Y")&&
		strstr( buffer, "TOTAL"))
		{
			if(!fgets(buffer,BSIZE,file))break;
			if(!fgets(buffer,BSIZE,file))break;
			if(!fgets(buffer,BSIZE,file))break;
			if(strstr( buffer, "SUM"))
			{
				int l,i;
				pdest =strstr(buffer, "SUM")+strlen("SUM")+1;
				l = strlen(pdest);
				for(i=0;i<l;i++) if(pdest[i]=='D' || pdest[i]=='E') pdest[i] ='E';
				if(sscanf(pdest+1,"%lf %lf %lf",&dipole[0],&dipole[1],&dipole[2])==3)
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
/*
static int getNumberOfTv(QuantumMechanicsModel *qmModel)
{
	int i;
	char tmp[100];
	int k = 0;
	for(i=0;i<qmModel->molecule.nAtoms;i++)
	{
		sprintf(tmp,"%s",qmModel->molecule.atoms[i].prop.symbol);
		uppercase(tmp);
		if(!strcmp(tmp,"TV")) k++;
	}
	return k;
}
*/
/*****************************************************************************/
static boolean getGradientMopac(char* fileNameOut, QuantumMechanicsModel *qmModel)
{
	FILE* file = NULL;
	char buffer[1024];
	char* pdest = NULL;
	boolean Ok = FALSE;
	double tmp;
	int i;
	int j;
	int k;
	int dum;
	/* int nTv= 0;*/

 	file = fopen(fileNameOut, "r");
	if(!file) return FALSE;
	/*
	nTv=getNumberOfTv(qmModel);
	printf("nTv=%d\n",nTv);
	*/
	for(i=0;i<qmModel->molecule.nAtoms;i++)
		for(j=0;j<3;j++) 
			qmModel->molecule.atoms[i].gradient[j] =0;

	 while(!feof(file))
	 {
		if(!fgets(buffer,BSIZE,file))break;
		pdest = strstr( buffer, "PARAMETER     ATOM    TYPE            VALUE       GRADIENT");
		if(pdest) 
		{
			for(k=0;k<3*qmModel->molecule.nAtoms;k++)
			{
				if(!fgets(buffer,BSIZE,file))break;
				pdest = strstr( buffer, "CARTESIAN");
				if(!pdest) break;
				j = 0;
				if(strstr( buffer, "CARTESIAN X"))j=0;
				if(strstr( buffer, "CARTESIAN Y"))j=1;
				if(strstr( buffer, "CARTESIAN Z"))j=2;
				if(2==sscanf(buffer,"%d %d",&dum,&i))
				{
					i--;
					if(i<0 || i> qmModel->molecule.nAtoms-1) 
					{
						fclose(file);
						return FALSE;
					}
					pdest = strstr( buffer, "CARTESIAN");
					if(sscanf(pdest+12,"%lf %lf",&tmp,&qmModel->molecule.atoms[i].gradient[j])!=2)
					{
						fclose(file);
						return FALSE;
					}
				}
				
			}
			Ok = TRUE;
			break;
	 	}
		pdest = strstr( buffer, "Cartesian Gradients"); /* MOZYME Keyword */
		if(pdest) 
		{
			char td[100];
			int d;
			if(!fgets(buffer,BSIZE,file))break; /*Atom       X  ....*/
			if(!fgets(buffer,BSIZE,file))break; /* backspace */
			for(k=0;k<qmModel->molecule.nAtoms;k++)
			{
				if(!fgets(buffer,BSIZE,file)) /* 1  O    0.000   -4.566    0.027  */
				{
					fclose(file);
					return FALSE;
				}
				if(1!=sscanf(buffer,"%d",&i)) break;
				i--;
				if(i<0 || i>=qmModel->molecule.nAtoms)
				{
					fclose(file);
					return FALSE;
				}
				if(sscanf(buffer,"%d %s %lf %lf %lf",&d, td, 
						&qmModel->molecule.atoms[i].gradient[0],
						&qmModel->molecule.atoms[i].gradient[1],
						&qmModel->molecule.atoms[i].gradient[2]
					 )
						!=5)
					{
						fclose(file);
						return FALSE;
					}
			}
			Ok = TRUE;
			break;
	 	}
	 }
	fclose(file);
	return Ok;
}
/*****************************************************************************/
static char* runOneMopac(QuantumMechanicsModel* qmModel, char* keyWords)
{
	FILE* file = NULL;
	FILE* fileSH = NULL;
	int j;
	char* fileNameIn = NULL;
	char* fileNameOut = NULL;
	char* fileNameSH = NULL;
	char multiplicityStr[100];
	char buffer[1024];
	Molecule m = qmModel->molecule;
	int rank = 0;
#ifdef ENABLE_MPI
	MPI_Comm_rank( MPI_COMM_WORLD,&rank);
#endif

#ifdef OS_WIN32
	char c='%';
#endif

	if(m.nAtoms<1) return fileNameOut;
#ifndef OS_WIN32
	fileNameSH = strdup_printf("%s%sMopacOne%d.sh",qmModel->workDir,DIR_SEPARATOR_S,rank);
#else
	fileNameSH = strdup_printf("%s%sMopacOne%d.bat",qmModel->workDir,DIR_SEPARATOR_S,rank);
#endif
 	fileSH = fopen(fileNameSH, "w");
	if(!fileSH) return FALSE;

	getMultiplicityName(qmModel->molecule.spinMultiplicity, multiplicityStr);

	fileNameIn = strdup_printf("%s%sOne%d.mop",qmModel->workDir,DIR_SEPARATOR_S,rank);
 	file = fopen(fileNameIn, "w");
	if(!file) 
	{
 		if(fileNameIn) free(fileNameIn);
 		if(fileNameOut) free(fileNameOut);
 		if(fileNameSH) free(fileNameSH);
		return FALSE;
	}
	if(qmModel->molecule.spinMultiplicity>1)
		fprintf(file,"%s UHF CHARGE=%d %s\n",keyWords,qmModel->molecule.totalCharge,multiplicityStr);
	else
		fprintf(file,"%s CHARGE=%d %s\n",keyWords,qmModel->molecule.totalCharge,multiplicityStr);
	fprintf(file,"\n");
	fprintf(file,"Mopac file generated by Gabedit\n");

	for(j=0;j<m.nAtoms;j++)
	{
	fprintf(file," %s %f %d %f %d %f %d\n", 
			m.atoms[j].prop.symbol,
			m.atoms[j].coordinates[0],
			1,
			m.atoms[j].coordinates[1],
			1,
			m.atoms[j].coordinates[2],
			1
			);
	}
	fclose(file);
#ifndef OS_WIN32
	fprintf(fileSH,"%s %s\n",qmModel->nameCommand,fileNameIn);
	fclose(fileSH);
	sprintf(buffer,"chmod u+x %s",fileNameSH);
	system(buffer);
	system(fileNameSH);
#else
	fprintf(fileSH,"\"%s\" \"%s\"\n",qmModel->nameCommand,fileNameIn);
	fclose(fileSH);
	sprintf(buffer,"\"%s\"",fileNameSH);
	system(buffer);
#endif

	unlink(fileNameIn);
	unlink(fileNameSH);
 	if(fileNameIn) free(fileNameIn);
 	if(fileNameSH) free(fileNameSH);
	fileNameOut = strdup_printf("%s%sOne%d.out",qmModel->workDir,DIR_SEPARATOR_S,rank);
	return fileNameOut;
}
/**********************************************************************/
static QuantumMechanicsModel newMopacModel(char* method, char* dirName, char* nameCommand, Constraints constraints, FILE* logfile)
{
	QuantumMechanicsModel qmModel = newQuantumMechanicsModel(method, dirName, nameCommand, NULL, NULL, constraints, logfile);

	qmModel.klass->calculateGradient = calculateGradientMopac;
	qmModel.klass->calculateEnergy = calculateEnergyMopac;

	return qmModel;
}
/**********************************************************************/
static void calculateGradientMopac(QuantumMechanicsModel* qmModel)
{
	int i;
	int j;
	Molecule m = qmModel->molecule;
	char* keyWords = NULL;
	char* fileOut = NULL;
	if(!qmModel) return;
	if(qmModel->molecule.nAtoms<1) return;
	if(!qmModel->method) return;
	keyWords = strdup_printf("%s 1SCF GRAD",qmModel->method);
	fileOut = runOneMopac(qmModel, keyWords);

	if(fileOut)
	{
		for(j=0;j<3;j++)
			for( i=0; i<m.nAtoms;i++)
				m.atoms[i].gradient[j] = 0.0;
		if(!getGradientMopac(fileOut, qmModel))
		{
			printf(("Problem : I cannot compute the Gradient by OpenMopac... "));
			exit(1);
			return;
		}
		getEnergyMopac(fileOut, &qmModel->molecule.potentialEnergy);
		addH4Correction(qmModel,TRUE);
		addSRBCorrection(qmModel,TRUE);
		addD3Correction(qmModel,TRUE);
		addWallCorrection(qmModel,TRUE);
		getDipoleMopac(fileOut, qmModel->molecule.dipole);
		/* printf("Energy = %f\n",qmModel->molecule.potentialEnergy);*/
		free(fileOut);
	}

}
/**********************************************************************/
static void calculateEnergyMopac(QuantumMechanicsModel* qmModel)
{
	char* keyWords = NULL;
	char* fileOut = NULL;
	if(!qmModel) return;
	if(qmModel->molecule.nAtoms<1) return;
	if(!qmModel->method) return;
	keyWords = strdup_printf("%s 1SCF",qmModel->method);
	fileOut = runOneMopac(qmModel, keyWords);
	if(fileOut)
	{
		getEnergyMopac(fileOut, &qmModel->molecule.potentialEnergy);
		addH4Correction(qmModel,FALSE);
		addSRBCorrection(qmModel,FALSE);
		addD3Correction(qmModel,FALSE);
		addWallCorrection(qmModel,FALSE);
		getDipoleMopac(fileOut, qmModel->molecule.dipole);
		free(fileOut);
	}

}
/**********************************************************************/
QuantumMechanicsModel createMopacModel (Molecule* mol, char* method, char* dirName, char* nameCommand, Constraints constraints, FILE* logfile)
{
	QuantumMechanicsModel qmModel = newMopacModel(method, dirName, nameCommand, constraints, logfile);

	qmModel.molecule = *(mol->klass->copy(mol));
	qmModel.molecule.constraints = constraints;
	qmModel.klass->setRattleConstraintsParameters(&qmModel);
	
	return qmModel;
}
/*****************************************************************************/
static boolean getDipoleFireFly(char* fileNameOut, double* dipole)
{
	FILE* file = NULL;
	char buffer[1024];

 	file = fopen(fileNameOut, "r");
	if(!file) return FALSE;

	 while(!feof(file))
	 {
		if(!fgets(buffer,BSIZE,file))break;
		if(
		strstr( buffer, "DX")&&
		strstr( buffer, "DY")&&
		strstr( buffer, "DZ") &&
		strstr( buffer, "(DEBYE)")
		)
		{
			int l,i;
			if(!fgets(buffer,BSIZE,file))break;
			l = strlen(buffer);
			for(i=0;i<l;i++) if(buffer[i]=='D' || buffer[i]=='E') buffer[i] ='E';
			l = sscanf(buffer,"%lf %lf %lf",&dipole[0], &dipole[1], &dipole[2]);
			if(l!=3) for(i=0;i<3;i++) dipole[i] = 0.0;
		}
	 }
	fclose(file);
	return FALSE;
}
/**********************************************************************/
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
		if(!pdest) pdest = strstr( buffer, "FINAL ENERGY IS");
		if(pdest) 
		{
			pdest = strstr( buffer, "S");
			if(pdest)
			{
				int l = strlen(pdest);
				int i;
				for(i=0;i<l;i++) if(pdest[i]=='D' || pdest[i]=='E') pdest[i] ='E';
				if(sscanf(pdest+1,"%lf",energy)==1)
				{
					OK = TRUE;
					if(strstr( buffer, "FINAL ENERGY IS")) *energy *= AUTOKCAL;
					/* break;*/
				}
			}
		}
	 }
	fclose(file);
	return OK;
}
/*****************************************************************************/
static boolean getGradientFireFly(char* fileNameOut, QuantumMechanicsModel *qmModel)
{
	FILE* file = NULL;
	char buffer[1024];
	char stmp[1024];
	char* pdest = NULL;
	boolean Ok = FALSE;
	int itmp;
	int i;
	int j;

 	file = fopen(fileNameOut, "r");
	if(!file) return FALSE;
	 while(!feof(file))
	 {
		if(!fgets(buffer,BSIZE,file))break;
		pdest = strstr( buffer, "ATOM                 E'X               E'Y               E'Z");
		if(pdest) 
		{
			for(i=0;i<qmModel->molecule.nAtoms;i++)
			{
				if(!fgets(buffer,BSIZE,file))break;
				if(sscanf(buffer,"%d %s %lf %lf %lf",&itmp, stmp,
							&qmModel->molecule.atoms[i].gradient[0],
							&qmModel->molecule.atoms[i].gradient[1],
							&qmModel->molecule.atoms[i].gradient[2]
							)!=5)
				{
					fclose(file);
					return FALSE;
				}
				for(j=0;j<3;j++) qmModel->molecule.atoms[i].gradient[j] *= AUTOKCAL/BOHRTOANG;
			}
			Ok = TRUE;
			break;
	 	}
	 }
	fclose(file);
	return Ok;
}
/*****************************************************************************/
static char* runOneFireFly(QuantumMechanicsModel* qmModel, char* keyWords)
{
	FILE* file = NULL;
	FILE* fileSH = NULL;
	int j;
	char* fileNameIn = NULL;
	char* fileNameOut = NULL;
	char* fileNameSH = NULL;
	char multiplicityStr[100];
	char buffer[1024];
	Molecule m = qmModel->molecule;
	char* fileNamePrefix = NULL;
	int rank = 0;
#ifdef ENABLE_MPI
	MPI_Comm_rank( MPI_COMM_WORLD,&rank);
#endif
#ifdef OS_WIN32
	char c='%';
#endif

	if(m.nAtoms<1) return fileNameOut;
#ifndef OS_WIN32
	fileNameSH = strdup_printf("%s%sFireFlyOne%d.sh",qmModel->workDir,DIR_SEPARATOR_S,rank);
#else
	fileNameSH = strdup_printf("%s%sFireFlyOne%d.bat",qmModel->workDir,DIR_SEPARATOR_S,rank);
#endif
 	fileSH = fopen(fileNameSH, "w");
	if(!fileSH) return FALSE;
#ifdef OS_WIN32
	fprintf(fileSH,"@echo off\n");
	fprintf(fileSH,"set PATH=%cPATH%c;\"%s\"\n",c,c,fireflyDirectory);
#endif

	getMultiplicityName(qmModel->molecule.spinMultiplicity, multiplicityStr);

	fileNameIn = strdup_printf("%s%sOne%d.inp",qmModel->workDir,DIR_SEPARATOR_S,rank);
	fileNameOut = strdup_printf("%s%sOne%d.out",qmModel->workDir,DIR_SEPARATOR_S,rank);


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
		if(qmModel->molecule.spinMultiplicity==1)
			fprintf(file," $CONTRL SCFTYP=RHF $END\n");
		else
			fprintf(file," $CONTRL SCFTYP=UHF $END\n");
	}

	fprintf(file," $CONTRL ICHARG=%d MULT=%d $END\n",qmModel->molecule.totalCharge,qmModel->molecule.spinMultiplicity);
	if(strstr(keyWords,"GBASIS"))
	{
		sscanf(strstr(keyWords,"GBASIS"),"%s",buffer);
		fprintf(file," $BASIS %s $END\n",buffer);
	}
	fprintf(file," $DATA\n");
	fprintf(file,"Molecule specification\n");
	fprintf(file,"C1\n");
	for(j=0;j<m.nAtoms;j++)
	{
		char* symbol = m.atoms[j].prop.symbol;
		SAtomsProp prop = propAtomGet(symbol);
		fprintf(file,"%s %f %f %f %f\n", 
			symbol,
			(double)prop.atomicNumber,
			m.atoms[j].coordinates[0],
			m.atoms[j].coordinates[1],
			m.atoms[j].coordinates[2]
			);
	}
	fprintf(file," $END\n");
	fclose(file);
	fileNamePrefix = strdup_printf("%s%sWorkFF",qmModel->workDir,DIR_SEPARATOR_S);
#ifndef OS_WIN32
	if(!strcmp(qmModel->nameCommand,"pcgamess") || !strcmp(qmModel->nameCommand,"nohup pcgamess")
	|| !strcmp(qmModel->nameCommand,"firefly") || !strcmp(qmModel->nameCommand,"nohup firefly"))
	{
		fprintf(fileSH,"mkdir %stmp\n",fileNamePrefix);
		fprintf(fileSH,"cd %stmp\n",fileNamePrefix);
		fprintf(fileSH,"cp %s input\n",fileNameIn);
		fprintf(fileSH,"%s -p -o %s\n",qmModel->nameCommand,fileNameOut);
		fprintf(fileSH,"cd ..\n");
		fprintf(fileSH,"rm PUNCH\n");
		fprintf(fileSH,"/bin/rm -r  %stmp\n",fileNamePrefix);
	}
	else
		fprintf(fileSH,"%s %s",qmModel->nameCommand,fileNameIn);
#else
	 if(!strcmp(qmModel->nameCommand,"pcgamess") ||
	 !strcmp(qmModel->nameCommand,"firefly") )
	{
        	fprintf(fileSH,"mkdir \"%stmp\"\n",fileNamePrefix);
		addUnitDisk(fileSH, fileNamePrefix);
	 	fprintf(fileSH,"cd \"%stmp\"\n",fileNamePrefix);
         	fprintf(fileSH,"copy \"%s\" input\n",fileNameIn);
         	fprintf(fileSH,"%s -p -o \"%s\"\n",qmModel->nameCommand,fileNameOut);
	 	fprintf(fileSH,"cd ..\n");
         	fprintf(fileSH,"del PUNCH 2>nul\n");
         	fprintf(fileSH,"del /Q  \"%stmp\"\n",fileNamePrefix);
         	fprintf(fileSH,"rmdir  \"%stmp\"\n",fileNamePrefix);
	}
	else
		fprintf(fileSH,"%s %s",qmModel->nameCommand,fileNameIn);
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
	unlink(fileNameIn);
	unlink(fileNameSH);
 	if(fileNamePrefix) free(fileNamePrefix);
 	if(fileNameIn) free(fileNameIn);
 	if(fileNameSH) free(fileNameSH);
	return fileNameOut;
}
/**********************************************************************/
static QuantumMechanicsModel newFireFlyModel(char* method, char* dirName, char* nameCommand, Constraints constraints, FILE* logfile)
{
	QuantumMechanicsModel qmModel = newQuantumMechanicsModel(method, dirName, nameCommand, NULL, NULL, constraints, logfile);

	qmModel.klass->calculateGradient = calculateGradientFireFly;
	qmModel.klass->calculateEnergy = calculateEnergyFireFly;

	return qmModel;
}
/**********************************************************************/
static void calculateGradientFireFly(QuantumMechanicsModel* qmModel)
{
	int i;
	int j;
	Molecule m = qmModel->molecule;
	char* keyWords = NULL;
	char* fileOut = NULL;
	if(!qmModel) return;
	if(qmModel->molecule.nAtoms<1) return;
	if(!qmModel->method) return;
	keyWords = strdup_printf("RUNTYP=GRADIENT GBASIS=%s",qmModel->method);
	fileOut = runOneFireFly(qmModel, keyWords);

	if(fileOut)
	{
		for(j=0;j<3;j++)
			for( i=0; i<m.nAtoms;i++)
				m.atoms[i].gradient[j] = 0.0;
		if(!getGradientFireFly(fileOut, qmModel))
		{
#ifdef OS_WIN32
			char* comm = strdup_printf("type %s",fileOut);
#else
			char* comm = strdup_printf("cat %s",fileOut);
#endif
			printf(("Problem : I cannot caculate the Gradient... "));
			printf(("Calculation Stopped "));
			system(comm);
			free(fileOut);
			free(comm);
			exit(1);
			return;
		}
		getEnergyFireFly(fileOut, &qmModel->molecule.potentialEnergy);
		addH4Correction(qmModel,TRUE);
		addSRBCorrection(qmModel,TRUE);
		addD3Correction(qmModel,TRUE);
		addWallCorrection(qmModel,TRUE);
		getDipoleFireFly(fileOut, qmModel->molecule.dipole);
		free(fileOut);
	}

}
/**********************************************************************/
static void calculateEnergyFireFly(QuantumMechanicsModel* qmModel)
{
	char* keyWords = NULL;
	char* fileOut = NULL;
	if(!qmModel) return;
	if(qmModel->molecule.nAtoms<1) return;
	if(!qmModel->method) return;
	keyWords = strdup_printf("RUNTYP=Energy GBASIS=%s",qmModel->method);
	fileOut = runOneFireFly(qmModel, keyWords);
	if(fileOut)
	{
		getEnergyFireFly(fileOut, &qmModel->molecule.potentialEnergy);
		addH4Correction(qmModel,FALSE);
		addSRBCorrection(qmModel,FALSE);
		addD3Correction(qmModel,FALSE);
		addWallCorrection(qmModel,FALSE);
		getDipoleFireFly(fileOut, qmModel->molecule.dipole);
		free(fileOut);
	}

}
/**********************************************************************/
QuantumMechanicsModel createFireFlyModel (Molecule* mol, char* method, char* dirName, char* nameCommand, Constraints constraints, FILE* logfile)
{
	QuantumMechanicsModel qmModel = newFireFlyModel(method,dirName, nameCommand, constraints, logfile);

	qmModel.molecule = *(mol->klass->copy(mol));
	qmModel.molecule.constraints = constraints;
	qmModel.klass->setRattleConstraintsParameters(&qmModel);
	
	return qmModel;
}
/**********************************************************************/
/*
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
					*energy *= AUTOKCAL;
					OK = TRUE;
					// break;
				}
			}
		}
	 }
	fclose(file);
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
/*****************************************************************************/
boolean getGradientGaussian(char* fileNameOut, QuantumMechanicsModel *qmModel)
{
	FILE* file = NULL;
	char buffer[1024];
	char stmp[1024];
	char* pdest = NULL;
	boolean Ok = FALSE;
	int itmp;
	int i;
	int j;

 	file = fopen(fileNameOut, "r");
	if(!file) return FALSE;
	 while(!feof(file))
	 {
		if(!fgets(buffer,BSIZE,file))break;
		pdest = strstr( buffer, " Forces (Hartrees/Bohr)");
		if(pdest) 
		{
			do{
				if(!fgets(buffer,BSIZE,file))break;
			}while(!feof(file)&&!strstr(buffer,"------"));
			if(!strstr(buffer,"------"))break;
			for(i=0;i<qmModel->molecule.nAtoms;i++)
			{
				if(!fgets(buffer,BSIZE,file))break;
				/* printf("%s\n",buffer);*/
				if(sscanf(buffer,"%d %s %lf %lf %lf",&itmp, stmp,
							&qmModel->molecule.atoms[i].gradient[0],
							&qmModel->molecule.atoms[i].gradient[1],
							&qmModel->molecule.atoms[i].gradient[2]
							)!=5)
				{
					fclose(file);
					return FALSE;
				}
				for(j=0;j<3;j++) qmModel->molecule.atoms[i].gradient[j] *= AUTOKCAL/BOHRTOANG;
				for(j=0;j<3;j++) qmModel->molecule.atoms[i].gradient[j] = - qmModel->molecule.atoms[i].gradient[j];
			}
			Ok = TRUE;
			break;
	 	}
	 }
	fclose(file);
	return Ok;
}
/*****************************************************************************/
static char* runOneGaussian(QuantumMechanicsModel* qmModel, char* keyWords)
{
	FILE* file = NULL;
	FILE* fileSH = NULL;
	int j;
	char* fileNameIn = NULL;
	char* fileNameChk = NULL;
	char* fileNameOut = NULL;
	char* fileNameSH = NULL;
	char buffer[1024];
	Molecule m = qmModel->molecule;
	int rank = 0;
#ifdef ENABLE_MPI
	MPI_Comm_rank( MPI_COMM_WORLD,&rank);
#endif
#ifdef OS_WIN32
	char c='%';
#endif

	if(m.nAtoms<1) return fileNameOut;
#ifndef OS_WIN32
	fileNameSH = strdup_printf("%s%sGaussianOne%d.sh",qmModel->workDir,DIR_SEPARATOR_S,rank);
#else
	fileNameSH = strdup_printf("%s%sGaussianOne%d.bat",qmModel->workDir,DIR_SEPARATOR_S,rank);
#endif
 	fileSH = fopen(fileNameSH, "w");
	if(!fileSH) return FALSE;
#ifdef OS_WIN32
	fprintf(fileSH,"@echo off\n");
	fprintf(fileSH,"set PATH=%cPATH%c;\"%s\"\n",c,c,gaussianDirectory);
#endif

	fileNameIn = strdup_printf("%s%sOne%d.com",qmModel->workDir,DIR_SEPARATOR_S,rank);
	fileNameOut = strdup_printf("%s%sOne%d.log",qmModel->workDir,DIR_SEPARATOR_S,rank);
	fileNameChk = strdup_printf("%s%sOne%d",qmModel->workDir,DIR_SEPARATOR_S,rank);

 	file = fopen(fileNameIn, "w");
	if(!file) 
	{
 		if(fileNameIn) free(fileNameIn);
 		if(fileNameOut) free(fileNameOut);
 		if(fileNameSH) free(fileNameSH);
		return FALSE;
	}
	fprintf(file,"%%chk=%s\n",fileNameChk);
	fprintf(file,"# %s\n",keyWords);
	fprintf(file,"\n");
	fprintf(file,"! ======================================================\n");
	fprintf(file,"!  Input file for Gaussian\n"); 
	fprintf(file,"! ======================================================\n");
	fprintf(file,"\n");

	fprintf(file,"%d %d\n",qmModel->molecule.totalCharge,qmModel->molecule.spinMultiplicity);
	for(j=0;j<m.nAtoms;j++)
	{
		char* symbol = m.atoms[j].prop.symbol;
		fprintf(file,"%s %f %f %f\n", 
			symbol,
			m.atoms[j].coordinates[0],
			m.atoms[j].coordinates[1],
			m.atoms[j].coordinates[2]
			);
	}
	fprintf(file,"\n");
	fclose(file);
#ifndef OS_WIN32
	fprintf(fileSH,"%s %s",qmModel->nameCommand,fileNameIn);
#else
	fprintf(fileSH,"%s %s",qmModel->nameCommand,fileNameIn);
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
	unlink(fileNameIn);
	unlink(fileNameSH);
 	if(fileNameIn) free(fileNameIn);
 	if(fileNameSH) free(fileNameSH);
	return fileNameOut;
}
/**********************************************************************/
static QuantumMechanicsModel newGaussianModel(char* method, char* dirName, char* nameCommand, Constraints constraints, FILE* logfile)
{
	QuantumMechanicsModel qmModel = newQuantumMechanicsModel(method, dirName, nameCommand, NULL,  NULL, constraints, logfile);

	qmModel.klass->calculateGradient = calculateGradientGaussian;
	qmModel.klass->calculateEnergy = calculateEnergyGaussian;

	return qmModel;
}
/**********************************************************************/
static void calculateGradientGaussian(QuantumMechanicsModel* qmModel)
{
	int i;
	int j;
	Molecule m = qmModel->molecule;
	char* keyWords = NULL;
	char* fileOut = NULL;
	if(!qmModel) return;
	if(qmModel->molecule.nAtoms<1) return;
	if(!qmModel->method) return;
	if(qmModel->firstRun) keyWords = strdup_printf("%s NoSym Force SCF(XQC) Test Pop=None",qmModel->method);
	else keyWords = strdup_printf("%s Force SCF(XQC) Guess=Read Test Pop=None",qmModel->method);
	fileOut = runOneGaussian(qmModel, keyWords);

	if(fileOut)
	{
		for(j=0;j<3;j++)
			for( i=0; i<m.nAtoms;i++)
				m.atoms[i].gradient[j] = 0.0;
		if(!getGradientGaussian(fileOut, qmModel))
		{
#ifdef OS_WIN32
			char* comm = strdup_printf("type %s",fileOut);
#else
			char* comm = strdup_printf("cat %s",fileOut);
#endif
			printf(("Problem : I cannot caculate the Gradient... "));
			printf(("Calculation Stopped "));
			system(comm);
			free(fileOut);
			free(comm);
			exit(1);
			return;
		}
		getEnergyGaussian(fileOut, &qmModel->molecule.potentialEnergy);
		addH4Correction(qmModel,TRUE);
		addSRBCorrection(qmModel,TRUE);
		addD3Correction(qmModel,TRUE);
		addWallCorrection(qmModel,TRUE);
		getDipoleGaussian(fileOut, qmModel->molecule.dipole);
		free(fileOut);
	}

}
/*****************************************************************************/
static boolean getDipoleGaussian(char* fileNameOut, double* dipole)
{
	FILE* file = NULL;
	char buffer[1024];
	char* pdest = NULL;

 	file = fopen(fileNameOut, "r");
	if(!file) return FALSE;

	 while(!feof(file))
	 {
		if(!fgets(buffer,BSIZE,file))break;
		if(
		strstr( buffer, " X=")&&
		strstr( buffer, " Y=")&&
		strstr( buffer, " Z=")
		)
		{
			int l,i;
			pdest =strstr(buffer, "X=")+strlen("X=")+1;
			l = strlen(pdest);
			for(i=0;i<l;i++) if(pdest[i]=='D' || pdest[i]=='E') pdest[i] ='E';
			l = 0;
			pdest =strstr(buffer, "X=")+strlen("X=")+1;
			l += sscanf(pdest,"%lf",&dipole[0]);
			pdest =strstr(buffer, "Y=")+strlen("Y=")+1;
			l += sscanf(pdest,"%lf",&dipole[1]);
			pdest =strstr(buffer, "Z=")+strlen("Z=")+1;
			l += sscanf(pdest,"%lf",&dipole[2]);
			if(l!=3) for(i=0;i<3;i++) dipole[i] = 0.0;
		}
	 }
	fclose(file);
	return FALSE;
}
/**********************************************************************/
static void calculateEnergyGaussian(QuantumMechanicsModel* qmModel)
{
	char* keyWords = NULL;
	char* fileOut = NULL;
	if(!qmModel) return;
	if(qmModel->molecule.nAtoms<1) return;
	if(!qmModel->method) return;
	keyWords = strdup_printf("%s  SCF(XQC) Test",qmModel->method);
	fileOut = runOneGaussian(qmModel, keyWords);
	if(fileOut)
	{
		getEnergyGaussian(fileOut, &qmModel->molecule.potentialEnergy);
		addH4Correction(qmModel,FALSE);
		addSRBCorrection(qmModel,FALSE);
		addWallCorrection(qmModel,FALSE);
		getDipoleGaussian(fileOut, qmModel->molecule.dipole);
		free(fileOut);
	}

}
/**********************************************************************/
QuantumMechanicsModel createGaussianModel (Molecule* mol, char* method, char* dirName, char* nameCommand, Constraints constraints, FILE* logfile)
{
	QuantumMechanicsModel qmModel = newGaussianModel(method,dirName, nameCommand, constraints, logfile);

	qmModel.molecule = *(mol->klass->copy(mol));
	qmModel.molecule.constraints = constraints;
	qmModel.klass->setRattleConstraintsParameters(&qmModel);
	
	return qmModel;
}
/*****************************************************************************/
static boolean getDipoleOrca(char* fileNameOut, double* dipole)
{
	FILE* file = NULL;
	char buffer[1024];
	char* pdest = NULL;

 	file = fopen(fileNameOut, "r");
	if(!file) return FALSE;

	 while(!feof(file))
	 {
		if(!fgets(buffer,BSIZE,file))break;
		if(
		strstr( buffer, "Total")&&
		strstr( buffer, "Dipole")&&
		strstr( buffer, "Moment")
		)
		{
			int l,i;
			pdest =strstr(buffer, ":")+strlen(":")+1;
			l = strlen(pdest);
			for(i=0;i<l;i++) if(pdest[i]=='D' || pdest[i]=='E') pdest[i] ='E';
			if(sscanf(pdest,"%lf %lf %lf",&dipole[0],&dipole[1],&dipole[2])==3)
			{
				for(i=0;i<3;i++) dipole[i] *= AUTODEB;
			}
		}
	 }
	fclose(file);
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
			*energy *= AUTOKCAL;
			return TRUE;
		}
	 }
	fclose(file);
	return FALSE;
}
/*****************************************************************************/
boolean getGradientOrca(char* fileNameOut, QuantumMechanicsModel *qmModel)
{
	FILE* file = NULL;
	char buffer[1024];
	char* pdest = NULL;
	boolean Ok = FALSE;
	int i;
	int k;
	char* gradTag = "# The current gradient in Eh/bohr";

 	file = fopen(fileNameOut, "r");
	if(!file) return FALSE;
	 while(!feof(file))
	 {
		if(!fgets(buffer,BSIZE,file))break;
		pdest = strstr( buffer, gradTag);
		if(pdest)
		{
			i=0;
			k=0;
	 		while(!feof(file) && i< qmModel->molecule.nAtoms)
	 		{
				if(!fgets(buffer,BSIZE,file))break;
				if(strstr(buffer,"#")) continue;
				/* printf("%s\n",buffer);*/
				if(sscanf(buffer,"%lf",&qmModel->molecule.atoms[i].gradient[k])!=1)
				{
					fclose(file);
					return FALSE;
				}
				k++;
				if(k==3) { k = 0; i++;}
			}
			Ok = TRUE;
		}
	 }
	fclose(file);
	if(Ok)
	{
		for(i=0;i<qmModel->molecule.nAtoms;i++)
		{
			for(k=0;k<3;k++) qmModel->molecule.atoms[i].gradient[k] *= AUTOKCAL/BOHRTOANG;
			//for(k=0;k<3;k++) qmModel->molecule.atoms[i].gradient[k] = - qmModel->molecule.atoms[i].gradient[k];
		}
	}
	return Ok;
}
/*****************************************************************************/
static char* runOneOrca(QuantumMechanicsModel* qmModel, char* keyWords)
{
	FILE* file = NULL;
	FILE* fileSH = NULL;
	int i;
	int nV = 0;
	char* fileNameIn = NULL;
	char* fileNameOut = NULL;
	char* fileNameSH = NULL;
	char buffer[1024];
	Molecule* mol = &qmModel->molecule;
	char* NameCommandOrca = qmModel->nameCommand;
	int rank = 0;
#ifdef ENABLE_MPI
	MPI_Comm_rank( MPI_COMM_WORLD,&rank);
#endif
#ifdef OS_WIN32
	char c='%';
#endif

	if(mol->nAtoms<1) return fileNameOut;
#ifndef OS_WIN32
	fileNameSH = strdup_printf("%s%sOrcaOne%d.sh",qmModel->workDir,DIR_SEPARATOR_S,rank);
#else
	fileNameSH = strdup_printf("%s%sOrcaOne%d.bat",qmModel->workDir,DIR_SEPARATOR_S,rank);
#endif
 	fileSH = fopen(fileNameSH, "w");
	if(!fileSH) return FALSE;
#ifdef OS_WIN32
	fprintf(fileSH,"@echo off\n");
	fprintf(fileSH,"set PATH=%cPATH%c;\"%s\"\n",c,c,orcaDirectory);
#endif

	fileNameIn = strdup_printf("%s%sOne%d.inp",qmModel->workDir,DIR_SEPARATOR_S,rank);
	fileNameOut = strdup_printf("%s%sOne%d.out",qmModel->workDir,DIR_SEPARATOR_S,rank);

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
#ifndef OS_WIN32
	sprintf(buffer,"chmod u+x %s",fileNameSH);
	system(buffer);
	system(fileNameSH);
#else
	sprintf(buffer,"\"%s\"",fileNameSH);
	system(buffer);
#endif
	unlink(fileNameIn);
	unlink(fileNameSH);
 	if(fileNameIn) free(fileNameIn);
 	if(fileNameSH) free(fileNameSH);
	return fileNameOut;
}
/**********************************************************************/
static QuantumMechanicsModel newOrcaModel(char* method, char* dirName, char* nameCommand, Constraints constraints, FILE* logfile)
{
	QuantumMechanicsModel qmModel = newQuantumMechanicsModel(method, dirName, nameCommand, NULL, NULL, constraints, logfile);

	qmModel.klass->calculateGradient = calculateGradientOrca;
	qmModel.klass->calculateEnergy = calculateEnergyOrca;

	return qmModel;
}
/**********************************************************************/
static void calculateGradientOrca(QuantumMechanicsModel* qmModel)
{
	int i;
	int j;
	Molecule m = qmModel->molecule;
	char* keyWords = NULL;
	char* fileOut = NULL;
	if(!qmModel) return;
	if(qmModel->molecule.nAtoms<1) return;
	if(!qmModel->method) return;
	if(qmModel->firstRun) keyWords = strdup_printf("%s ENGRAD ",qmModel->method);
	else keyWords = strdup_printf("%s ENGRAD ",qmModel->method);
	fileOut = runOneOrca(qmModel, keyWords);

	if(fileOut)
	{
		for(j=0;j<3;j++)
			for( i=0; i<m.nAtoms;i++)
				m.atoms[i].gradient[j] = 0.0;
		if(!getGradientOrca(fileOut, qmModel))
		{
#ifdef OS_WIN32
			char* comm = strdup_printf("type %s",fileOut);
#else
			char* comm = strdup_printf("cat %s",fileOut);
#endif
			printf(("Problem : I cannot caculate the Gradient... "));
			printf(("Calculation Stopped "));
			system(comm);
			free(fileOut);
			free(comm);
			exit(1);
			return;
		}
		getEnergyOrca(fileOut, &qmModel->molecule.potentialEnergy);
		addH4Correction(qmModel,TRUE);
		addSRBCorrection(qmModel,TRUE);
		addD3Correction(qmModel,TRUE);
		addWallCorrection(qmModel,TRUE);
		getDipoleOrca(fileOut, qmModel->molecule.dipole);
		free(fileOut);
	}

}
/**********************************************************************/
static void calculateEnergyOrca(QuantumMechanicsModel* qmModel)
{
	char* keyWords = NULL;
	char* fileOut = NULL;
	if(!qmModel) return;
	if(qmModel->molecule.nAtoms<1) return;
	if(!qmModel->method) return;
	keyWords = strdup_printf("%s ",qmModel->method);
	fileOut = runOneOrca(qmModel, keyWords);
	if(fileOut)
	{
		getEnergyOrca(fileOut, &qmModel->molecule.potentialEnergy);
		addH4Correction(qmModel,FALSE);
		addSRBCorrection(qmModel,FALSE);
		addWallCorrection(qmModel,FALSE);
		getDipoleOrca(fileOut, qmModel->molecule.dipole);
		free(fileOut);
	}

}
/**********************************************************************/
QuantumMechanicsModel createOrcaModel (Molecule* mol, char* method, char* dirName, char* nameCommand, Constraints constraints, FILE* logfile)
{
	QuantumMechanicsModel qmModel = newOrcaModel(method,dirName, nameCommand, constraints, logfile);

	qmModel.molecule = *(mol->klass->copy(mol));
	qmModel.molecule.constraints = constraints;
	qmModel.klass->setRattleConstraintsParameters(&qmModel);
	
	return qmModel;
}
/**********************************************************************/
static boolean getDipoleGeneric(char* fileNameOut, double* dipole)
{
	FILE* file = NULL;
	char buffer[1024];
	int i;
 	file = fopen(fileNameOut, "r");
	if(!file) return FALSE;
	if(!fgets(buffer,BSIZE,file)) { fclose(file); return FALSE;}/* first line for energy in Hartree*/
	if(!fgets(buffer,BSIZE,file)) { fclose(file); return FALSE;}/* dipole in au */
	for(i=0;i<strlen(buffer);i++) if(buffer[i]=='D' || buffer[i]=='d') buffer[i] ='E';
	if(sscanf(buffer,"%lf %lf %lf",&dipole[0],&dipole[1],&dipole[2])==3)
	{
		for(i=0;i<3;i++) dipole[i] *= AUTODEB;
		fclose(file);
		return TRUE;
	}
	fclose(file);
	return FALSE;
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
		*energy *= AUTOKCAL;
		return TRUE;
	}
	fclose(file);
	return FALSE;
}
/*****************************************************************************/
boolean getGradientGeneric(char* fileNameOut, QuantumMechanicsModel *qmModel)
{
	FILE* file = NULL;
	char buffer[1024];
	boolean Ok = FALSE;
	int i;
	int j;

 	file = fopen(fileNameOut, "r");
	if(!file) return FALSE;
	if(!fgets(buffer,BSIZE,file)) { fclose(file); return FALSE;}/* first line for energy in Hartree*/
	if(!fgets(buffer,BSIZE,file)) { fclose(file); return FALSE;}/* dipole in au */
	for(i=0;i<qmModel->molecule.nAtoms;i++)
	{
		if(!fgets(buffer,BSIZE,file))break;
		for(j=0;j<strlen(buffer);j++) if(buffer[j]=='D' || buffer[j]=='d') buffer[j] ='E';
		if(sscanf(buffer,"%lf %lf %lf",
					&qmModel->molecule.atoms[i].gradient[0],
					&qmModel->molecule.atoms[i].gradient[1],
					&qmModel->molecule.atoms[i].gradient[2]
					)!=3)
		{
			fclose(file);
			return FALSE;
		}
		for(j=0;j<3;j++) qmModel->molecule.atoms[i].gradient[j] *= AUTOKCAL/BOHRTOANG;
		for(j=0;j<3;j++) qmModel->molecule.atoms[i].gradient[j] = - qmModel->molecule.atoms[i].gradient[j];
	}
	Ok = TRUE;
	fclose(file);
	return Ok;
}
/*****************************************************************************/
char* runOneGeneric(QuantumMechanicsModel* qmModel, char* keyWords)
{
	FILE* file = NULL;
	FILE* fileSH = NULL;
	char* fileNameIn = NULL;
	char* fileNameOut = NULL;
	char* fileNameSH = NULL;
	char buffer[1024];
	Molecule* mol = &qmModel->molecule;
	char* NameCommandGeneric = qmModel->nameCommand;
	int rank = 0;
	int type = 0;
#ifdef ENABLE_MPI
	MPI_Comm_rank( MPI_COMM_WORLD,&rank);
#endif
#ifdef OS_WIN32
	char c='%';
#endif

	if(mol->nAtoms<1) return fileNameOut;
#ifndef OS_WIN32
	fileNameSH = strdup_printf("%s%sGenericOne%d.sh",qmModel->workDir,DIR_SEPARATOR_S,rank);
#else
	fileNameSH = strdup_printf("%s%sGenericOne%d.bat",qmModel->workDir,DIR_SEPARATOR_S,rank);
#endif
 	fileSH = fopen(fileNameSH, "w");
	if(!fileSH) return FALSE;
#ifdef OS_WIN32
	fprintf(fileSH,"@echo off\n");
#endif

	fileNameIn = strdup_printf("%s%sOne%d.inp",qmModel->workDir,DIR_SEPARATOR_S,rank);
	fileNameOut = strdup_printf("%s%sOne%d.out",qmModel->workDir,DIR_SEPARATOR_S,rank);

 	file = fopen(fileNameIn, "w");
	if(!file) 
	{
 		if(fileNameIn) free(fileNameIn);
 		if(fileNameOut) free(fileNameOut);
 		if(fileNameSH) free(fileNameSH);
		return FALSE;
	}
	/*
	fprintf(file,"# ======================================================\n");
	fprintf(file,"#  Generic input file made in CChemI\n"); 
	fprintf(file,"# ======================================================\n");
	*/
	if(strstr(keyWords,"ENGRAD")) type = 1;
	fprintf(file,"%d\n",type);
	mol->klass->addMolecule(mol,file);
	fclose(file);

#ifndef OS_WIN32
	fprintf(fileSH,"%s %s %s",NameCommandGeneric,fileNameIn,fileNameOut);
	fclose(fileSH);
	sprintf(buffer,"chmod u+x %s",fileNameSH);
	system(buffer);
	system(fileNameSH);
#else
	fprintf(fileSH,"\"%s\" \"%s\" \"%s\"",NameCommandGeneric,fileNameIn,fileNameOut);
	fclose(fileSH);
	sprintf(buffer,"\"%s\"",fileNameSH);
	system(buffer);
#endif
	unlink(fileNameIn);
	unlink(fileNameSH);
 	if(fileNameIn) free(fileNameIn);
 	if(fileNameSH) free(fileNameSH);
	return fileNameOut;
}
/**********************************************************************/
static QuantumMechanicsModel newGenericModel(char* method, char* dirName, char* nameCommand, Constraints constraints, FILE* logfile)
{
	QuantumMechanicsModel qmModel = newQuantumMechanicsModel(method, dirName, nameCommand, NULL, NULL, constraints, logfile);

	qmModel.klass->calculateGradient = calculateGradientGeneric;
	qmModel.klass->calculateEnergy = calculateEnergyGeneric;

	return qmModel;
}
/**********************************************************************/
static void calculateGradientGeneric(QuantumMechanicsModel* qmModel)
{
	int i;
	int j;
	Molecule m = qmModel->molecule;
	char* keyWords = NULL;
	char* fileOut = NULL;
	if(!qmModel) return;
	if(qmModel->molecule.nAtoms<1) return;
	if(!qmModel->method) return;
	if(qmModel->firstRun) keyWords = strdup_printf("%s ENGRAD ",qmModel->method);
	else keyWords = strdup_printf("%s ENGRAD ",qmModel->method);
	fileOut = runOneGeneric(qmModel, keyWords);

	if(fileOut)
	{
		for(j=0;j<3;j++)
			for( i=0; i<m.nAtoms;i++)
				m.atoms[i].gradient[j] = 0.0;
		if(!getGradientGeneric(fileOut, qmModel))
		{
#ifdef OS_WIN32
			char* comm = strdup_printf("type %s",fileOut);
#else
			char* comm = strdup_printf("cat %s",fileOut);
#endif
			printf(("Problem : I cannot caculate the Gradient... "));
			printf(("Calculation Stopped "));
			system(comm);
			free(fileOut);
			free(comm);
			exit(1);
			return;
		}
		getEnergyGeneric(fileOut, &qmModel->molecule.potentialEnergy);
		addH4Correction(qmModel,TRUE);
		addSRBCorrection(qmModel,TRUE);
		addD3Correction(qmModel,TRUE);
		addWallCorrection(qmModel,TRUE);
		getDipoleGeneric(fileOut, qmModel->molecule.dipole);
		free(fileOut);
	}

}
/**********************************************************************/
static void calculateEnergyGeneric(QuantumMechanicsModel* qmModel)
{
	char* keyWords = NULL;
	char* fileOut = NULL;
	if(!qmModel) return;
	if(qmModel->molecule.nAtoms<1) return;
	if(!qmModel->method) return;
	keyWords = strdup_printf("%s ",qmModel->method);
	fileOut = runOneGeneric(qmModel, keyWords);
	if(fileOut)
	{
		getEnergyGeneric(fileOut, &qmModel->molecule.potentialEnergy);
		addH4Correction(qmModel,FALSE);
		addSRBCorrection(qmModel,FALSE);
		addD3Correction(qmModel,FALSE);
		addWallCorrection(qmModel,FALSE);
		getDipoleGeneric(fileOut, qmModel->molecule.dipole);
		free(fileOut);
	}

}
/**********************************************************************/
QuantumMechanicsModel createGenericModel (Molecule* mol, char* method, char* dirName, char* nameCommand, Constraints constraints, FILE* logfile)
{
	QuantumMechanicsModel qmModel = newGenericModel(method,dirName, nameCommand, constraints, logfile);

	qmModel.molecule = *(mol->klass->copy(mol));
	qmModel.molecule.constraints = constraints;
	qmModel.klass->setRattleConstraintsParameters(&qmModel);
	
	return qmModel;
}
/**********************************************************************/
/*
static boolean getDipoleOpenBabel(char* fileNameOut, Molecule* mol, double dipole[])
{
	FILE* file = NULL;
	char buffer[1024];
	boolean ok = FALSE;

 	file = fopen(fileNameOut, "r");
	if(!file) return FALSE;

	 while(!feof(file))
	 {
		if(!fgets(buffer,BSIZE,file))break;
		if(
		strstr( buffer, "IDX")&&
		strstr( buffer, "CHARGE")
		)
		{
			int i = 0;
			int k;
			for(k=0;k<3;k++) dipole[k] = 0.0;
			for(i=0;i<mol->nAtoms;i++)
			{
				if(!fgets(buffer,BSIZE,file))break;
				for(k=0;k<3;k++) dipole[k] += mol->atoms[i].coordinates[k]*atof(buffer)*BOHRTOANG;
			}
			if(i!=mol->nAtoms) break;
			for(k=0;k<3;k++) dipole[i] *= AUTODEB;
			ok = TRUE;
		}
	 }
	fclose(file);
	return ok;
}
*/
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
static char* runOneOpenBabel(QuantumMechanicsModel* qmModel, char* keyWords)
{
	FILE* fileSH = NULL;
	char* fileNameIn = NULL;
	char* fileNameOut = NULL;
	char* fileNameSH = NULL;
	char buffer[1024];
	Molecule* mol = &qmModel->molecule;
	char* NameCommandOpenBabel = qmModel->nameCommand;
	int rank = 0;
#ifdef ENABLE_MPI
	MPI_Comm_rank( MPI_COMM_WORLD,&rank);
#endif
#ifdef OS_WIN32
	char c='%';
#endif

	if(mol->nAtoms<1) return fileNameOut;
#ifndef OS_WIN32
	fileNameSH = strdup_printf("%s%sOpenBabelOne%d.sh",qmModel->workDir,DIR_SEPARATOR_S,rank);
#else
	fileNameSH = strdup_printf("%s%sOpenBabelOne%d.bat",qmModel->workDir,DIR_SEPARATOR_S,rank);
#endif
 	fileSH = fopen(fileNameSH, "w");
	if(!fileSH) return FALSE;
#ifdef OS_WIN32
	fprintf(fileSH,"@echo off\n");
	fprintf(fileSH,"set PATH=%cPATH%c;\"%s\"\n",c,c,openBabelDirectory);
#endif

	fileNameIn = strdup_printf("%s%sOne%d.hin",qmModel->workDir,DIR_SEPARATOR_S,rank);
	fileNameOut = strdup_printf("%s%sOne%d.out",qmModel->workDir,DIR_SEPARATOR_S,rank);

	if(!mol->klass->saveHIN(mol,fileNameIn))
	{
 		if(fileNameIn) free(fileNameIn);
 		if(fileNameOut) free(fileNameOut);
 		if(fileNameSH) free(fileNameSH);
		return FALSE;
	}
#ifndef OS_WIN32
	if(!strcmp(NameCommandOpenBabel,"obenergy") || !strcmp(NameCommandOpenBabel,"nohup obenergy"))
	{
		fprintf(fileSH,"%s -ff gaff %s > %s 2>/dev/null\n",NameCommandOpenBabel,fileNameIn,fileNameOut);
		fprintf(fileSH,"exit\n");
	}
	else
		fprintf(fileSH,"%s %s >%s 2>/dev/null",NameCommandOpenBabel,fileNameIn, fileNameOut);
#else
	 if(!strcmp(NameCommandOpenBabel,"obenergy") )
	{
		if(strstr(openBabelDirectory,"\"")) fprintf(fileSH,"set PATH=%s;%cPATH%c\n",openBabelDirectory,'%','%');
		else fprintf(fileSH,"set PATH=\"%s\";%cPATH%c\n",openBabelDirectory,'%','%');
		fprintf(fileSH,"%s -ff gaff %s > %s\n",NameCommandOpenBabel,fileNameIn,fileNameOut);
		fprintf(fileSH,"exit\n");
	}
	else
		fprintf(fileSH,"%s %s > %s",NameCommandOpenBabel,fileNameIn, fileNameOut);
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

	/*
	sprintf(buffer,"cat %s",fileNameOut);
	system(buffer);
	*/
#else
	sprintf(buffer,"\"%s\"",fileNameSH);
	system(buffer);
#endif
	unlink(fileNameIn);
	unlink(fileNameSH);
 	if(fileNameIn) free(fileNameIn);
 	if(fileNameSH) free(fileNameSH);
	return fileNameOut;
}
/**********************************************************************/
static QuantumMechanicsModel newOpenBabelModel(char* method, char* dirName, char* nameCommand, Constraints constraints, FILE* logfile)
{
	QuantumMechanicsModel qmModel = newQuantumMechanicsModel(method, dirName, nameCommand, NULL, NULL, constraints, logfile);

	qmModel.klass->calculateGradient = calculateGradientOpenBabel;
	qmModel.klass->calculateEnergy = calculateEnergyOpenBabel;

	return qmModel;
}
/**********************************************************************/
/*
static void calculateGradientOpenBabelNumeric(QuantumMechanicsModel* qmModel)
{
	double dx = 0.01;
        Molecule* mol = &qmModel->molecule;
        int nAtoms = mol->nAtoms;
	int i,k;
	double Em,Ep;

        for(i=0;i<nAtoms;i++)
        for(k=0;k<3;k++)
                mol->atoms[i].gradient[k] = 0.0;

        for(i=0;i<nAtoms;i++)
        for(k=0;k<3;k++)
        {
                mol->atoms[i].coordinates[k] += dx;
                qmModel->klass->calculateEnergy(qmModel);
                Ep = mol->potentialEnergy;


                mol->atoms[i].coordinates[k] -= 2*dx;
                qmModel->klass->calculateEnergy(qmModel);
                Em = mol->potentialEnergy;

                mol->atoms[i].gradient[k] = (Ep-Em)/dx/2;
                mol->atoms[i].coordinates[k] += dx;
        }
        qmModel->klass->calculateEnergy(qmModel);

}
*/
/*****************************************************************************/
static boolean getGradientsOpenBabel(char* fileNameOut, Molecule* mol)
{
	FILE* file = NULL;
	char buffer[1024];
	char* pdest = NULL;
	//char* energyTag = "TOTAL ENERGY =";
	char* energyTag = "FINAL ENERGY:";
	char* gradTag = "Gradients:";
	boolean kj = FALSE;
	double conv = 1.0;
	int i;

 	file = fopen(fileNameOut, "r");
	if(!file) return FALSE;
	 while(!feof(file))
	 {
		if(!fgets(buffer,BSIZE,file))break;
		pdest = strstr( buffer, energyTag);
		if(pdest &&sscanf(pdest+strlen(energyTag)+1,"%lf",&mol->potentialEnergy)==1)
		{
			if(strstr(pdest,"kJ")) { kj = TRUE;}
			break;
		}
	 }
	 while(!feof(file))
	 {
		if(!fgets(buffer,BSIZE,file))break;
		if(strstr(buffer, gradTag))
		{
			for(i=0;i<mol->nAtoms;i++)
			{
				if(!fgets(buffer,BSIZE,file))break;
				//printf("%s\n",buffer);
				if(sscanf(buffer,"%lf %lf %lf",
					&mol->atoms[i].gradient[0],
					&mol->atoms[i].gradient[1],
					&mol->atoms[i].gradient[2]
					)!=3) break;
			}
			break;
		}
	 }
	if(kj) conv /= KCALTOKJ;
	mol->potentialEnergy *= conv;
	for(i=0;i<mol->nAtoms;i++)
	{
		mol->atoms[i].gradient[0] *= conv;
		mol->atoms[i].gradient[1] *= conv;
		mol->atoms[i].gradient[2] *= conv;
	}

	fclose(file);
	return FALSE;
}
/**********************************************************************/
static void calculateGradientOpenBabelAnalytic(QuantumMechanicsModel* qmModel)
{
	char* keyWords = NULL;
	char* fileOut = NULL;
	if(!qmModel) return;
	if(qmModel->molecule.nAtoms<1) return;
	if(!qmModel->method) return;
	keyWords = strdup_printf("%s ",qmModel->method);
	fileOut = runOneOpenBabel(qmModel, keyWords);
	if(fileOut)
	{
		getGradientsOpenBabel(fileOut, &qmModel->molecule);
		addH4Correction(qmModel,TRUE);
		addSRBCorrection(qmModel,TRUE);
		addWallCorrection(qmModel,TRUE);
		//getDipoleOpenBabel(fileOut, &qmModel->molecule, qmModel->molecule.dipole);
		qmModel->molecule.klass->computeDipole(&qmModel->molecule);
		free(fileOut);
	}
}
/**********************************************************************/
static void calculateGradientOpenBabel(QuantumMechanicsModel* qmModel)
{
	//calculateGradientOpenBabelNumeric(qmModel);
	calculateGradientOpenBabelAnalytic(qmModel);
}
/**********************************************************************/
static void calculateEnergyOpenBabel(QuantumMechanicsModel* qmModel)
{
	char* keyWords = NULL;
	char* fileOut = NULL;
	if(!qmModel) return;
	if(qmModel->molecule.nAtoms<1) return;
	if(!qmModel->method) return;
	keyWords = strdup_printf("%s ",qmModel->method);
	fileOut = runOneOpenBabel(qmModel, keyWords);
	if(fileOut)
	{
		getEnergyOpenBabel(fileOut, &qmModel->molecule.potentialEnergy);
		addSRBCorrection(qmModel,FALSE);
		addWallCorrection(qmModel,FALSE);
		//getDipoleOpenBabel(fileOut, &qmModel->molecule, qmModel->molecule.dipole);
		qmModel->molecule.klass->computeDipole(&qmModel->molecule);
		free(fileOut);
	}

}
/**********************************************************************/
QuantumMechanicsModel createOpenBabelModel (Molecule* mol, char* method, char* dirName, char* nameCommand, Constraints constraints, FILE* logfile)
{
	QuantumMechanicsModel qmModel = newOpenBabelModel(method,dirName, nameCommand, constraints, logfile);

	qmModel.molecule = *(mol->klass->copy(mol));
	qmModel.molecule.constraints = constraints;
	qmModel.klass->setRattleConstraintsParameters(&qmModel);
	
	return qmModel;
}
/**********************************************************************/
static void calculateEnergyN2P2(QuantumMechanicsModel* qmModel)
{
	int err = 0;
	Molecule* mol = &qmModel->molecule;
	if(!qmModel) return;
	if(qmModel->molecule.nAtoms<1) return;
	err = interfaceCChemIComputeEnergy(qmModel->interfaceLibN2P2, mol);
	if(!err)
	{
		addH4Correction(qmModel,TRUE);
		addSRBCorrection(qmModel,TRUE);
		addD3Correction(qmModel,TRUE);
		addWallCorrection(qmModel,TRUE);
		//getDipoleN2P2(fileOut, qmModel->molecule.dipole);
		if(qmModel->interfaceLibN2P2ES)
		{
			double charge = 0;
			err = interfaceCChemIESComputeChargeAndDipole(qmModel->interfaceLibN2P2ES, mol, &charge);
			/*
			if(!err)
			{
				fprintf(stderr,"charge NN = %f\n",charge);
				fprintf(stderr,"dipole NN = %f %f %f\n",mol->dipole[0], mol->dipole[1],mol->dipole[2]);
			}
			else fprintf(stderr,"Error NN charge & dipole calculationn\n");
			mol->dipole[0]= mol->dipole[1]=mol->dipole[2]=0;
			qmModel->molecule.klass->computeDipole(mol);
			fprintf(stderr,"dipole from molecule = %f %f %f\n",mol->dipole[0], mol->dipole[1],mol->dipole[2]);
			*/
		}
		
	}
	else
	{
		fprintf(stderr,"Error : in interfaceCChemIComputeForces, Stopped\n");
		exit(1);
	}

}
/**********************************************************************/
static void calculateGradientN2P2(QuantumMechanicsModel* qmModel)
{
	int err = 0;
	Molecule* mol = &qmModel->molecule;
	if(!qmModel) return;
	if(qmModel->molecule.nAtoms<1) return;

	err = interfaceCChemIComputeGradients(qmModel->interfaceLibN2P2, mol);
	if(!err)
	{
		addH4Correction(qmModel,TRUE);
		addSRBCorrection(qmModel,TRUE);
		addD3Correction(qmModel,TRUE);
		addWallCorrection(qmModel,TRUE);
		//getDipoleN2P2(fileOut, qmModel->molecule.dipole);
		if(qmModel->interfaceLibN2P2ES)
		{
			double charge = 0;
			err = interfaceCChemIESComputeChargeAndDipole(qmModel->interfaceLibN2P2ES, mol, &charge);
			/*
			if(!err)
			{
				fprintf(stderr,"charge NN = %f\n",charge);
				fprintf(stderr,"dipole NN = %f %f %f\n",mol->dipole[0], mol->dipole[1],mol->dipole[2]);
			}
			else fprintf(stderr,"Error NN charge & dipole calculationn\n");
			*/
			
		}
	}
	else
	{
		fprintf(stderr,"Error : in interfaceCChemIComputeForces, Stopped\n");
		exit(1);
	}
}
/**********************************************************************/
static void calculateHessianN2P2(QuantumMechanicsModel* qmModel, double **F, double*** dmu)
{
	int err = 0;
	int nAtoms;
	int i;
	int j,k,c,index,id,jd;
	Molecule* mol = NULL;
	if(!qmModel) return;
	if(qmModel->molecule.nAtoms<1) return;

	if(qmModel->SRBParameters || qmModel->addD3Correction || qmModel->addWallCorrection || qmModel->H4Parameters || !qmModel->interfaceLibN2P2 || !qmModel->interfaceLibN2P2ES)
	{
		fprintf(stderr,"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n");
		fprintf(stderr,"Analytical calculation of Hessian with N2P2\n");
		fprintf(stderr,"AND empirical corrections is not yes implemented\n");
		fprintf(stderr,"Analytical calculation of Hessian with N2P2 is implemented only without corrections\n");
		fprintf(stderr,"Add dx=value to your input file to do that numerically\n");
		fprintf(stderr,"Example : dx=1e-3\n");
		fprintf(stderr,"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n");
		exit(1);
	}
	mol = &qmModel->molecule;
	nAtoms = mol->nAtoms;
	err = interfaceCChemIComputeHessian(qmModel->interfaceLibN2P2, mol, F);

	if(!err)
	{
		for(i=0;i<nAtoms;i++)
        		for(k=0;k<3;k++)
        		{
                		id=3*i+k;
                		for(j=0;j<=i;j++)
                		{
                        		double invm = 1.0/sqrt( mol->atoms[i].mass* mol->atoms[j].mass);
                        		for(c = 0;c<3;c++)
                        		{
                                		jd = 3*j+c;
                                		if(jd>id) continue;
                                		index = jd + id*(id+1)/2;
                                		(*F)[index] *= invm;
                        		}
                		}
        	}
	}

	if(!err)
	{
		if(qmModel->interfaceLibN2P2ES)
		{
			double charge = 0;
			err = interfaceCChemIESComputedDipole(qmModel->interfaceLibN2P2ES, mol, dmu, &charge);
			/*
			if(!err)
			{
				fprintf(stderr,"charge NN = %f\n",charge);
				fprintf(stderr,"dipole NN = %f %f %f\n",mol->dipole[0], mol->dipole[1],mol->dipole[2]);
			}
			else fprintf(stderr,"Error NN charge & dipole calculationn\n");
			*/
			
		}
	}
	else
	{
		fprintf(stderr,"Error : in interfaceCChemIComputeHessian, Stopped\n");
		exit(1);
	}

}
/**********************************************************************/
static QuantumMechanicsModel newN2P2Model(char*  HNNDir, Constraints constraints, FILE* logfile)
{
	QuantumMechanicsModel qmModel = newQuantumMechanicsModel(NULL, NULL, NULL, HNNDir, NULL, constraints, logfile);
	double cfenergy = 1.0/(AUTOKCAL);
	double cflength = ANGTOBOHR;
	double cfdipole = 1.0/(AUTODEB);

	qmModel.interfaceLibN2P2 = newInterfaceCChemI(HNNDir, cflength, cfenergy, 0);
	qmModel.interfaceLibN2P2ES = newInterfaceCChemIES(HNNDir, cflength, cfdipole, 0);
	qmModel.klass->calculateGradient = calculateGradientN2P2;
	qmModel.klass->calculateHessian = calculateHessianN2P2;
	qmModel.klass->calculateEnergy = calculateEnergyN2P2;
	return qmModel;
}
/**********************************************************************/
QuantumMechanicsModel createN2P2Model (Molecule* mol, char* HNNDir, Constraints constraints, FILE* logfile)
{
	QuantumMechanicsModel qmModel = newN2P2Model(HNNDir, constraints, logfile);

	qmModel.molecule = *(mol->klass->copy(mol));
	qmModel.molecule.constraints = constraints;
	qmModel.klass->setRattleConstraintsParameters(&qmModel);
	
	return qmModel;
}
/**********************************************************************/
#ifdef ENABLE_PYTHON
static void calculateEnergyTM(QuantumMechanicsModel* qmModel)
{
	int err = 0;
	Molecule* mol = &qmModel->molecule;
	if(!qmModel) return;
	if(qmModel->molecule.nAtoms<1) return;
	err = interfaceTMComputeEnergy(qmModel->interfaceTM, mol);
	if(!err)
	{
		addH4Correction(qmModel,TRUE);
		addSRBCorrection(qmModel,TRUE);
		addD3Correction(qmModel,TRUE);
		addWallCorrection(qmModel,TRUE);
		//getDipoleTM(fileOut, qmModel->molecule.dipole);
		
	}
	else
	{
		fprintf(stderr,"Error : in interfaceTMComputeForces, Stopped\n");
		exit(1);
	}

}
/**********************************************************************/
static void calculateGradientTM(QuantumMechanicsModel* qmModel)
{
	int err = 0;
	Molecule* mol = &qmModel->molecule;
	if(!qmModel) return;
	if(qmModel->molecule.nAtoms<1) return;

	err = interfaceTMComputeGradients(qmModel->interfaceTM, mol);
	if(!err)
	{
		addH4Correction(qmModel,TRUE);
		addSRBCorrection(qmModel,TRUE);
		addD3Correction(qmModel,TRUE);
		addWallCorrection(qmModel,TRUE);
		//getDipoleTM(fileOut, qmModel->molecule.dipole);
	}
	else
	{
		fprintf(stderr,"Error : in interfaceTMComputeForces, Stopped\n");
		exit(1);
	}

}
/**********************************************************************/
static QuantumMechanicsModel newTMModel(char*  tmModule, Constraints constraints, FILE* logfile)
{
	QuantumMechanicsModel qmModel = newQuantumMechanicsModel(NULL, NULL, NULL, tmModule, NULL, constraints, logfile);

	qmModel.interfaceTM = newInterfaceTM(tmModule);
	qmModel.klass->calculateGradient = calculateGradientTM;
	qmModel.klass->calculateEnergy = calculateEnergyTM;
	return qmModel;
}
/**********************************************************************/
QuantumMechanicsModel createTMModel (Molecule* mol, char* tmModule, Constraints constraints, FILE* logfile)
{
	QuantumMechanicsModel qmModel = newTMModel(tmModule, constraints, logfile);

	qmModel.molecule = *(mol->klass->copy(mol));
	qmModel.molecule.constraints = constraints;
	qmModel.klass->setRattleConstraintsParameters(&qmModel);
	
	return qmModel;
}
#else
/**********************************************************************/
QuantumMechanicsModel createTMModel (Molecule* mol, char* tmModule, Constraints constraints, FILE* logfile)
{
	QuantumMechanicsModel qmModel;
	fprintf(stderr,"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n");
	fprintf(stderr," If you want to use the interface to TensorMol software\n");
	fprintf(stderr," You must compile cchemi with enable_python = 1 in CONFIG file\n");
	fprintf(stderr,"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n");
	exit(1);
	return qmModel;
}
#endif /* ENABLE_PYTHON */
