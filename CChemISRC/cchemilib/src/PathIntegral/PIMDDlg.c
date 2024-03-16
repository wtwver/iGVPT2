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
#include "../Utils/Utils.h"
#include "../Utils/AtomsProp.h"
#include "../Utils/Constants.h"
#include "../Utils/Types.h"
#include "../Molecule/Molecule.h"
#include "../QuantumMechanics/QuantumMechanics.h"
#include "../MolecularMechanics/MolecularMechanics.h"
#include "../EmpriricalCorrections/HydrogenBondCorrection.h"
#include "../EmpriricalCorrections/DispersionCorrection.h"
#include "../PathIntegral/Paths.h"

/*****************************************************************************/
void setPIMDOptions(FILE* file, int* updateFrequency, 
double* heatTime, double*equiTime, double* runTime, double* coolTime, 
double* heatTemp,  double*equiTemp, double*runTemp, double*coolTemp, 
double* stepSize, PIMDThermostatType* thermostat, PIMDTransformationType* transformation, double* friction, double* omegaMax, int* Nf, int* nBeads, int* nNH, int* nNHSteps, int* nSY)
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

	*thermostat = PIMDTHERMOSTATNONE;
	if(readOneInt(file,"PIMDThermostat",&itmp))*thermostat = itmp;

	*transformation = PIMDTRANSFORMATIONNONE;
	if(readOneInt(file,"PIMDTransformation",&itmp))*transformation = itmp;


	*friction=-1;
	readOneReal(file,"friction",friction);

	*nBeads=4;
	readOneInt(file,"nBeads",nBeads);

	*nNH=1;
	readOneInt(file,"lengthOfNHChains",nNH);

	*nNHSteps=5;
	readOneInt(file,"nNHSteps",nNHSteps);

	*nSY=3;
	readOneInt(file,"nSuzukiYoshida",nSY);

	*omegaMax=4000;
	readOneReal(file,"omegaMax",omegaMax);
	*Nf=50;
	readOneInt(file,"Nf",Nf);
}
/****************************************************************************************************************/
void pathIntegralMDDlg(char* inputFileName)
{
	QuantumMechanicsModel qmModel; 
	ForceField forceField; 
	QuantumMechanicsModel* pqmModel = NULL; 
	ForceField* pforceField = NULL; 
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
	char* fileNameTraj = NULL;
	char* fileNameProp = NULL;
	double friction=-1;
	PIMDThermostatType thermostat;
	PIMDTransformationType transformation;
	char* dirName = NULL;
	Constraints constraints = NOCONSTRAINTS;
	Molecule mol = *(readMolecule(inputFileName,TRUE));
	char* mopacCommand = strdup("/opt/mopac/MOPAC2009");
	char* fireflyCommand = strdup("firefly");
	char* orcaCommand=strdup("orca");
	char* gaussianCommand = strdup("g09");
	char* genericCommand = strdup("runGeneric");
	FILE* file = fopen(inputFileName,"rb");
	char* model = NULL;
	char* QMKeys = NULL;
	int cnt;
	int nNH;
	int nNHSteps;
	int nSY;
	int nBeads;
	double omegaMax;
	int Nf;
	
	setPIMDOptions(file, &updateFrequency, 
		&heatTime, &equiTime, &runTime, &coolTime,
		&heatTemp, &runTemp, &equiTemp, &coolTemp, &stepSize, 
		&thermostat, &transformation, &friction, &omegaMax, &Nf, &nBeads, &nNH, &nNHSteps, &nSY);

	if(!readOneString(file,"Model",&model)) model = strdup("MM");
	if(!readOneString(file,"QMKeys",&QMKeys)) 
	{
		if(!strcmp(model,"MOPAC")) QMKeys = strdup("PM6-DH2");
		else QMKeys = strdup("AM1");
	}
	uppercase(model);
	readOneString(file,"mopacCommand",&mopacCommand);
	readOneString(file,"fireflyCommand",&fireflyCommand);
	readOneString(file,"orcaCommand",&orcaCommand);
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

	if(!strcmp(model,"MM"))
	{
		ForceFieldOptions forceFieldOptions;
		setForceFieldOptions(file, &forceFieldOptions);
		if(forceFieldOptions.type==AMBER) forceField = createAmberModel(&mol,forceFieldOptions, stdout);
		else if(forceFieldOptions.type==PAIRWISE) forceField = createPairWiseModel(&mol,forceFieldOptions, stdout);
		pforceField = &forceField;	
	}
	else 
	{
		if(!strcmp(model,"MOPAC")) qmModel = createMopacModel(&mol, QMKeys, dirName, mopacCommand, constraints,stdout);
		else if(!strcmp(model,"FIREFLY")) qmModel = createFireFlyModel(&mol, QMKeys,dirName, fireflyCommand, constraints,stdout);
		else if(!strcmp(model,"ORCA")) qmModel = createOrcaModel(&mol, QMKeys, dirName, orcaCommand, constraints, stdout);
		else if(!strcmp(model,"GAUSSIAN")) qmModel = createGaussianModel(&mol, QMKeys,dirName, gaussianCommand, constraints,stdout);
		else qmModel = createGenericModel(&mol, QMKeys, dirName, genericCommand, constraints, stdout);
		setH4Correction(file,&qmModel);
		readOneBoolean(file,"addD3Correction",&qmModel.addD3Correction);
		pqmModel = &qmModel;	
	}

	
	runPIMD(pforceField, pqmModel, 
		updateFrequency, heatTime, equiTime, runTime, coolTime, heatTemp, equiTemp, runTemp, coolTemp, stepSize, 
		thermostat, transformation, friction, omegaMax, Nf, nBeads, nNH, nNHSteps, nSY, fileNameTraj, fileNameProp);

	if(!strcmp(model,"MM")) freeForceField(&forceField);
	else qmModel.klass->free(&qmModel);

	free(dirName);
	fclose(file);
}
/********************************************************************************/
