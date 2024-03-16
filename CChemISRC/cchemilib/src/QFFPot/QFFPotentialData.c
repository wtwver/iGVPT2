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

/* QFFPotentialData.c */
#include <math.h>
#include "../QFFPot/QFFPotentialData.h"
#include "../Utils/Utils.h"
#include "../Utils/AtomsProp.h"
#include "../Utils/Constants.h"
#include "../Utils/QL.h"

static void readData(QFFPotentialData* qffPotentialData, char* fileName);
static void convertToAU(QFFPotentialData* qffPotentialData);
static void convertToAU2(QFFPotentialData* qffPotentialData);
static void printData(QFFPotentialData* qffPotentialData);
static void freeQFFPotentialData(QFFPotentialData* qffPotentialData);
static void rotModes(QFFPotentialData* qffPotentialData, double u[3][3]);
static void computeQFFParameters(QFFPotentialData* qffPotentialData);
/**********************************************************************/
static void freeQFFPotentialData(QFFPotentialData* qffPotentialData)
{
	int nF = qffPotentialData->nFrequencies;
	int nAtoms = qffPotentialData->molecule.nAtoms;
	freeVectorDouble(&qffPotentialData->effectiveMasses);
	freeVectorDouble(&qffPotentialData->gradients);
	freeVectorDouble(&qffPotentialData->Be);
	freeMatrixDouble(&qffPotentialData->hessian, nF);
	freeCubeDouble(&qffPotentialData->modes, nF, nAtoms);
	freeCubeDouble(&qffPotentialData->MWModes, nF, nAtoms);
	freeCubeDouble(&qffPotentialData->cubic, nF, nF);
	freeCubeDouble(&qffPotentialData->coriolis, 3, nF);
	freeQuarticDouble(&qffPotentialData->quartic, nF, nF, nF);
}
/**********************************************************************/
QFFPotentialData newQFFPotentialData(int n, int nAtoms)
{
	QFFPotentialData qffPotentialData;
	qffPotentialData.klass = malloc(sizeof(QFFPotentialDataClass));
	qffPotentialData.klass->readData = readData;
	qffPotentialData.klass->convertToAU = convertToAU;
	qffPotentialData.klass->convertToAU2 = convertToAU2;
	qffPotentialData.klass->rotModes = rotModes;
	qffPotentialData.klass->free = freeQFFPotentialData;
	qffPotentialData.klass->print = printData;
	qffPotentialData.klass->computeQFFParameters = computeQFFParameters;
	qffPotentialData.nFrequencies = n;
	if(qffPotentialData.nFrequencies<=0) qffPotentialData.nFrequencies = 0;

	qffPotentialData.molecule = *(newMolecule());

	qffPotentialData.gradients = NULL;
	qffPotentialData.gradients = newVectorDouble(qffPotentialData.nFrequencies);
	initVectorDouble(qffPotentialData.gradients, qffPotentialData.nFrequencies, 0.0);

	qffPotentialData.hessian = newMatrixDouble(qffPotentialData.nFrequencies,qffPotentialData.nFrequencies);
	initMatrixDouble(qffPotentialData.hessian, qffPotentialData.nFrequencies, qffPotentialData.nFrequencies, 0.0);

	qffPotentialData.Be = newVectorDouble(3);
	initVectorDouble(qffPotentialData.Be, 3, 1.0);

	qffPotentialData.effectiveMasses = newVectorDouble(qffPotentialData.nFrequencies);
	initVectorDouble(qffPotentialData.effectiveMasses, qffPotentialData.nFrequencies, 1.0);

	qffPotentialData.cubic = newCubeDouble(qffPotentialData.nFrequencies, qffPotentialData.nFrequencies, qffPotentialData.nFrequencies);
	initCubeDouble(qffPotentialData.cubic, qffPotentialData.nFrequencies, qffPotentialData.nFrequencies, qffPotentialData.nFrequencies, 0.0);

	qffPotentialData.coriolis = newCubeDouble(3, qffPotentialData.nFrequencies, qffPotentialData.nFrequencies);
	initCubeDouble(qffPotentialData.coriolis, 3, qffPotentialData.nFrequencies, qffPotentialData.nFrequencies, 0.0);

	qffPotentialData.quartic = newQuarticDouble(qffPotentialData.nFrequencies, qffPotentialData.nFrequencies, qffPotentialData.nFrequencies, qffPotentialData.nFrequencies);
	initQuarticDouble(qffPotentialData.quartic, qffPotentialData.nFrequencies, qffPotentialData.nFrequencies, qffPotentialData.nFrequencies, qffPotentialData.nFrequencies, 0.0);

	qffPotentialData.modes = newCubeDouble(qffPotentialData.nFrequencies, nAtoms, 3);
	initCubeDouble(qffPotentialData.modes, qffPotentialData.nFrequencies, nAtoms, 3, 0.0);

	qffPotentialData.MWModes = newCubeDouble(qffPotentialData.nFrequencies, nAtoms, 3);
	initCubeDouble(qffPotentialData.MWModes, qffPotentialData.nFrequencies, nAtoms, 3, 0.0);
	qffPotentialData.qffPotParameters = NULL;

	return qffPotentialData;

}
/*****************************************************************************/
static void printQFFHessian(QFFPotentialData* qffPotentialData, char* fileName)
{
	FILE* file =fopen(fileName,"w");
	int n = qffPotentialData->nFrequencies;
	double cutoff = 1e-10;
	double** M = qffPotentialData->hessian;
	int i;
	for(i = 0;i<n; i++) 
	{
		if(fabs(M[i][i])>=cutoff)
      			fprintf(file,"%d %d %20.10f\n",i+1,i+1,M[i][i]);
	}
	fclose(file);
}
static void printQFFCubic(QFFPotentialData* qffPotentialData, char* fileName)
{
	FILE* file =fopen(fileName,"w");
	int n = qffPotentialData->nFrequencies;
	double cutoff = 1e-10;
	double*** M = qffPotentialData->cubic;
	int i,j,k;
	for(i = 0;i<n; i++) 
	{
		for(j = 0;j<n; j++) 
		for(k = 0;k<n; k++) 
			if(fabs(M[k][i][j])>=cutoff)
      			fprintf(file,"%d %d %d %20.10f\n",i+1,j+1,k+1,M[k][i][j]);
	}
	fclose(file);
}
static void printQFFCoriolis(QFFPotentialData* qffPotentialData, char* fileName)
{
	FILE* file =fopen(fileName,"w");
	int n = qffPotentialData->nFrequencies;
	int m = 3;
	double cutoff = 1e-10;
	double*** M = qffPotentialData->coriolis;
	int i,j,k;
	for(i = 0;i<n; i++) 
	{
		for(j = 0;j<n; j++) 
		for(k = 0;k<m; k++) 
			if(fabs(M[k][i][j])>=cutoff)
      			fprintf(file,"%d %d %d %20.10f\n",i+1,j+1,k+1,M[k][i][j]);
	}
	fclose(file);
}
static void printQFFRot(QFFPotentialData* qffPotentialData, char* fileName)
{
	FILE* file =fopen(fileName,"w");
	double cutoff = 1e-10;
	double* M = qffPotentialData->Be;
	int i;
	for(i = 0;i<3; i++) 
	{
		if(fabs(M[i])>=cutoff)
      			fprintf(file,"%d %20.10f\n",i+1,M[i]);
	}
	fclose(file);
}
static void printQFFQuartic(QFFPotentialData* qffPotentialData, char* fileName)
{
	FILE* file =fopen(fileName,"w");
	int i,j,k,l;
	double cutoff = 1e-10;
	double**** C = qffPotentialData->quartic;
	int n = qffPotentialData->nFrequencies;
	for(i = 0;i<n; i++) 
	{
		for(j = 0;j<n; j++) 
		{
			for(k = 0;k<n; k++) 
			{
				for(l = 0;l<n; l++) 
				if(fabs(C[i][j][k][l])>=cutoff)
      					fprintf(file,"%d %d %d %d %20.10f\n",i+1,j+1,k+1,l+1,C[i][j][k][l]);
			}
		}
	}
	fclose(file);
}
static void printQFFData(QFFPotentialData* qffPotentialData)
{
	printQFFHessian(qffPotentialData,"data.w");
	printQFFCubic(qffPotentialData,"data.f3");
	printQFFQuartic(qffPotentialData,"data.f4");
	printQFFRot(qffPotentialData,"data.b");
	printQFFCoriolis(qffPotentialData,"data.z");
}
/*****************************************************************************/
static void printModes(QFFPotentialData* qffPotentialData)
{
	printf("\n# imode iAtom xyz value\n");
	printf("Modes\n");
	printCubeDoubleCutOff(qffPotentialData->modes, qffPotentialData->nFrequencies, qffPotentialData->molecule.nAtoms, 3,-1);
	printf("END\n");
}
static void printGradients(QFFPotentialData* qffPotentialData)
{
	printf("\nGradients\n");
	printVectorDoubleCutOff(qffPotentialData->gradients, qffPotentialData->nFrequencies, 1e-10);
	printf("END\n");
}
static void printHessian(QFFPotentialData* qffPotentialData)
{
	printf("\nHessian\n");
	printMatrixDoubleCutOff(qffPotentialData->hessian, qffPotentialData->nFrequencies, qffPotentialData->nFrequencies,1e-10);
	printf("END\n");
}
static void printCubic(QFFPotentialData* qffPotentialData)
{
	printf("\nCubic\n");
	printCubeDoubleCutOff(qffPotentialData->cubic, qffPotentialData->nFrequencies, qffPotentialData->nFrequencies, qffPotentialData->nFrequencies,1e-10);
	printf("END\n");
}
static void printCoriolis(QFFPotentialData* qffPotentialData)
{
	printf("\nCoriolis\n");
	printCubeDoubleCutOff(qffPotentialData->coriolis, 3, qffPotentialData->nFrequencies, qffPotentialData->nFrequencies,1e-10);
	printf("END\n");
}
static void printRot(QFFPotentialData* qffPotentialData)
{
	printf("\nRotation constants\n");
	printVectorDoubleCutOff(qffPotentialData->Be, 3,1e-10);
	printf("END\n");
}
static void printQuartic(QFFPotentialData* qffPotentialData)
{
   printf("\nQuartic\n");
   printQuarticDoubleCutOff(qffPotentialData->quartic, qffPotentialData->nFrequencies, qffPotentialData->nFrequencies, qffPotentialData->nFrequencies, qffPotentialData->nFrequencies,1e-10);
	printf("END\n");
}
static void printData(QFFPotentialData* qffPotentialData)
{
	printModes(qffPotentialData);
	printGradients(qffPotentialData);
	printHessian(qffPotentialData);
	printCubic(qffPotentialData);
	printQuartic(qffPotentialData);
	printCoriolis(qffPotentialData);
	printRot(qffPotentialData);
	//printQFFData(qffPotentialData);
}
/*****************************************************************************/
static void convertToAU(QFFPotentialData* qffPotentialData)
{
	int i,j,k,l;
	int c;
	double conv = 1/AUTOCM1;
	double conv2 = conv*conv;
	double conv3 = conv2*sqrt(conv);
	double conv4 = conv2*conv;
	int nF = qffPotentialData->nFrequencies;
	double convAMUToAU = sqrt(AMUTOAU);
	
	/*
	conv2 *= convAMUToAU*convAMUToAU;
	conv3 *= convAMUToAU*convAMUToAU*convAMUToAU;
	conv4 *= convAMUToAU*convAMUToAU*convAMUToAU*convAMUToAU;
	// To convert in Hartree amu^-n/2 Bohr^-n
	*/

	//fprintf(stderr,"phi111 cm-1=%0.12e\n",qffPotentialData->cubic[1][1][1]);
	for(i=0;i<nF;i++)
	for(j=0;j<nF;j++)
	for(k=0;k<nF;k++)
		qffPotentialData->cubic[i][j][k] *=  conv3*sqrt(qffPotentialData->hessian[i][i]*qffPotentialData->hessian[j][j]*qffPotentialData->hessian[k][k]);

	for(i=0;i<nF;i++)
	for(j=0;j<nF;j++)
	for(k=0;k<nF;k++)
		qffPotentialData->cubic[i][j][k] =  getMaxCubicIJK( qffPotentialData->cubic, i, j, k);
	//fprintf(stderr,"phi111=%0.12e\n",qffPotentialData->cubic[1][1][1]);
	//exit(1);

	for(i=0;i<nF;i++)
	for(j=0;j<nF;j++)
	for(k=0;k<nF;k++)
	for(l=0;l<nF;l++)
		qffPotentialData->quartic[i][j][k][l] *=  conv4*sqrt(qffPotentialData->hessian[i][i]*qffPotentialData->hessian[j][j]*qffPotentialData->hessian[k][k]*qffPotentialData->hessian[l][l]);


	for(i=0;i<nF;i++)
	for(j=0;j<nF;j++)
	for(k=0;k<nF;k++)
	for(l=0;l<nF;l++)
		qffPotentialData->quartic[i][j][k][l] =  getMaxQuarticIJK( qffPotentialData->quartic,i,j,k,l);

	for(i=0;i<qffPotentialData->nFrequencies;i++)
		qffPotentialData->hessian[i][i] =  qffPotentialData->hessian[i][i]* qffPotentialData->hessian[i][i]*conv2;

	for(i=0;i<qffPotentialData->nFrequencies;i++)
	for(k=0;k<qffPotentialData->molecule.nAtoms;k++)
	for(c=0;c<3;c++)
		qffPotentialData->MWModes[i][k][c] *=  convAMUToAU;
}
/*****************************************************************************/
static void convertToAU2(QFFPotentialData* qffPotentialData)
{
	int i,j,k,l;
	int c;
	double conv = 1/AUTOCM1;
	double conv2 = conv*conv;
	double conv3 = conv2*sqrt(conv);
	double conv4 = conv2*conv;
	int nF = qffPotentialData->nFrequencies;
	double convAMUToAU = sqrt(AMUTOAU);
	
	conv2 *= convAMUToAU*convAMUToAU;
	conv3 *= convAMUToAU*convAMUToAU*convAMUToAU;
	conv4 *= convAMUToAU*convAMUToAU*convAMUToAU*convAMUToAU;
	// To convert in Hartree amu^-n/2 Bohr^-n

	//fprintf(stderr,"phi111 cm-1=%0.12e\n",qffPotentialData->cubic[1][1][1]);
	for(i=0;i<nF;i++)
	for(j=0;j<nF;j++)
	for(k=0;k<nF;k++)
		qffPotentialData->cubic[i][j][k] *=  conv3*sqrt(qffPotentialData->hessian[i][i]*qffPotentialData->hessian[j][j]*qffPotentialData->hessian[k][k]);

	for(i=0;i<nF;i++)
	for(j=0;j<nF;j++)
	for(k=0;k<nF;k++)
		qffPotentialData->cubic[i][j][k] =  getMaxCubicIJK( qffPotentialData->cubic, i, j, k);
	//fprintf(stderr,"phi111=%0.12e\n",qffPotentialData->cubic[1][1][1]);
	//exit(1);

	for(i=0;i<nF;i++)
	for(j=0;j<nF;j++)
	for(k=0;k<nF;k++)
	for(l=0;l<nF;l++)
		qffPotentialData->quartic[i][j][k][l] *=  conv4*sqrt(qffPotentialData->hessian[i][i]*qffPotentialData->hessian[j][j]*qffPotentialData->hessian[k][k]*qffPotentialData->hessian[l][l]);


	for(i=0;i<nF;i++)
	for(j=0;j<nF;j++)
	for(k=0;k<nF;k++)
	for(l=0;l<nF;l++)
		qffPotentialData->quartic[i][j][k][l] =  getMaxQuarticIJK( qffPotentialData->quartic,i,j,k,l);

	for(i=0;i<qffPotentialData->nFrequencies;i++)
		qffPotentialData->hessian[i][i] =  qffPotentialData->hessian[i][i]* qffPotentialData->hessian[i][i]*conv2;

	for(i=0;i<qffPotentialData->nFrequencies;i++)
	for(k=0;k<qffPotentialData->molecule.nAtoms;k++)
	for(c=0;c<3;c++)
		qffPotentialData->MWModes[i][k][c] *=  convAMUToAU;
}
/*****************************************************************************/
static void computeMWModes(QFFPotentialData* qffPotentialData)
{
	int i,k,c;
	int nF = qffPotentialData->nFrequencies;
	int nAtoms = qffPotentialData->molecule.nAtoms;
	for(i=0;i<nF;i++)
	for(k=0;k<nAtoms;k++)
	for(c=0;c<3;c++)
		qffPotentialData->MWModes[i][k][c] =  qffPotentialData->modes[i][k][c]*qffPotentialData->molecule.atoms[k].mass/sqrt(qffPotentialData->effectiveMasses[i]);
		//qffPotentialData->MWModes[i][k][c] =  qffPotentialData->modes[i][k][c]*sqrt(qffPotentialData->effectiveMasses[i]);
}
/*****************************************************************************/
static void rotModes(QFFPotentialData* qffPotentialData, double u[3][3])
{
	int i,k;
	int nF = qffPotentialData->nFrequencies;
	int nAtoms = qffPotentialData->molecule.nAtoms;
	double x,y,z;
	for(i=0;i<nF;i++)
	for(k=0;k<nAtoms;k++)
	{
		//qffPotentialData->MWModes[i][k][c] =  qffPotentialData->modes[i][k][c]*qffPotentialData->molecule.atoms[k].mass/sqrt(qffPotentialData->effectiveMasses[i]);
   		x = u[0][0] * qffPotentialData->MWModes[i][k][0] + u[0][1] * qffPotentialData->MWModes[i][k][1]  + u[0][2] * qffPotentialData->MWModes[i][k][2];
   		y = u[1][0] * qffPotentialData->MWModes[i][k][0] + u[1][1] * qffPotentialData->MWModes[i][k][1]  + u[1][2] * qffPotentialData->MWModes[i][k][2];
   		z = u[2][0] * qffPotentialData->MWModes[i][k][0] + u[2][1] * qffPotentialData->MWModes[i][k][1]  + u[2][2] * qffPotentialData->MWModes[i][k][2];

		qffPotentialData->MWModes[i][k][0] = x;
		qffPotentialData->MWModes[i][k][1] = y;
		qffPotentialData->MWModes[i][k][2] = z;
	}
}
/*****************************************************************************/
static void readData(QFFPotentialData* qffPotentialData, char* inputFileName)
{
	FILE* inputFile;
	Molecule mol = *(newMolecule());

        inputFile = fopen(inputFileName,"rb");
	if(!inputFile)
	{
		fprintf(stderr, "==========================================================\n");
		fprintf(stderr, "Sorry, I cannot opent the %s file\n", inputFileName);
		fprintf(stderr, "==========================================================\n");
		exit(1);
	}

	readOneInt(inputFile,"nFrequencies",&qffPotentialData->nFrequencies);
//	fprintf(stdout, "nFrequencies=%d\n", qffPotentialData->nFrequencies);
	fclose(inputFile);
	
	mol.klass->read(&mol, inputFileName);
	*qffPotentialData = newQFFPotentialData(qffPotentialData->nFrequencies,mol.nAtoms);
	mol.klass->free(&mol);

	qffPotentialData->molecule.klass->read(&qffPotentialData->molecule, inputFileName);
	if(qffPotentialData->molecule.nAtoms<1)
	{
		fprintf(stderr, "==========================================================\n");
		fprintf(stderr, "Sorry, I cannot read geometry from the %s file\n", inputFileName);
		fprintf(stderr, "==========================================================\n");
		exit(1);
	}

        inputFile = fopen(inputFileName,"rb");
	readVectorReal(inputFile,"Masses",qffPotentialData->nFrequencies, qffPotentialData->effectiveMasses);
	//fprintf(stderr,"Masse0 = %f\n",qffPotentialData->effectiveMasses[0]);

//	printf("nFrequencies=%d molecule.nAtoms=%d \n",qffPotentialData->nFrequencies, qffPotentialData->molecule.nAtoms);
	readCubeReal(inputFile,"Modes",qffPotentialData->nFrequencies, qffPotentialData->molecule.nAtoms, 3, qffPotentialData->modes);
	
	readVectorReal(inputFile,"Gradients",qffPotentialData->nFrequencies, qffPotentialData->gradients);/* qffPotentialData.gradients presently not used, min surface */
	readMatrixReal(inputFile,"Hessian",qffPotentialData->nFrequencies, qffPotentialData->nFrequencies, qffPotentialData->hessian);
	readCubeReal(inputFile,"Cubic",qffPotentialData->nFrequencies, qffPotentialData->nFrequencies, qffPotentialData->nFrequencies, qffPotentialData->cubic);
	readQuarticReal(inputFile,"Quartic",qffPotentialData->nFrequencies, qffPotentialData->nFrequencies, qffPotentialData->nFrequencies, qffPotentialData->nFrequencies, qffPotentialData->quartic);
	
	readVectorReal(inputFile,"Rotational",3, qffPotentialData->Be);
	readCubeReal(inputFile,"Coriolis",3, qffPotentialData->nFrequencies, qffPotentialData->nFrequencies, qffPotentialData->coriolis);
	//printData(qffPotentialData);
	//printQFFData(qffPotentialData);
	fclose(inputFile);
	computeMWModes(qffPotentialData);
}
/**********************************************************************/
static void computeQFFParameters(QFFPotentialData* qffPotentialData)
{
	int i,j,k,l;
	int nF = qffPotentialData->nFrequencies;
	QFFPotParameters* qffPotParameters = malloc(sizeof(QFFPotParameters));
	double cut = 1e-12;
	int m;
	

	// V1MR
	qffPotParameters->numberOf1MR = 0;
	qffPotParameters->qff1MR = NULL;
	for(i=0;i<nF;i++) if(fabs(qffPotentialData->hessian[i][i])>cut) qffPotParameters->numberOf1MR++;
	for(i=0;i<nF;i++) if(fabs(qffPotentialData->cubic[i][i][i])>cut) qffPotParameters->numberOf1MR++;
	for(i=0;i<nF;i++) if(fabs(qffPotentialData->quartic[i][i][i][i])>cut) qffPotParameters->numberOf1MR++;
	if(qffPotParameters->numberOf1MR>0) qffPotParameters->qff1MR= malloc(qffPotParameters->numberOf1MR*sizeof(QFFPotnMR));

	m = 0;
	for(i=0;i<nF;i++) if(fabs(qffPotentialData->hessian[i][i])>cut)
	{ 
		//fprintf(stderr,"H = %f\n",qffPotentialData->hessian[i][i]);
		qffPotParameters->qff1MR[m].numbers[0] = i;
		qffPotParameters->qff1MR[m].numbers[1] = 2;
		qffPotParameters->qff1MR[m].energy=0.5*qffPotentialData->hessian[i][i];
		qffPotParameters->qff1MR[m].grad=2*qffPotParameters->qff1MR[m].energy;
		m++;
	}
	for(i=0;i<nF;i++) if(fabs(qffPotentialData->cubic[i][i][i])>cut)
	{ 
		qffPotParameters->qff1MR[m].numbers[0] = i;
		qffPotParameters->qff1MR[m].numbers[1] = 3;
		qffPotParameters->qff1MR[m].energy=qffPotentialData->cubic[i][i][i]/6.0;
		qffPotParameters->qff1MR[m].grad=3*qffPotParameters->qff1MR[m].energy;
		m++;
	}
	for(i=0;i<nF;i++) if(fabs(qffPotentialData->quartic[i][i][i][i])>cut)
	{ 
		qffPotParameters->qff1MR[m].numbers[0] = i;
		qffPotParameters->qff1MR[m].numbers[1] = 4;
		qffPotParameters->qff1MR[m].energy=qffPotentialData->quartic[i][i][i][i]/24.0;
		qffPotParameters->qff1MR[m].grad=4*qffPotParameters->qff1MR[m].energy;
		m++;
	}
	

	// free Hessian

	// V2MR
	qffPotParameters->numberOf2MR = 0;
	qffPotParameters->qff2MR = NULL;
	for(i=0;i<nF;i++)
	for(j=0;j<nF;j++)
	{
		if(i==j) continue;
		if(fabs(qffPotentialData->cubic[i][i][j])>cut) qffPotParameters->numberOf2MR++;
		if(fabs(qffPotentialData->quartic[i][i][i][j])>cut) qffPotParameters->numberOf2MR++;
	}
	for(i=0;i<nF;i++)
	for(j=0;j<i;j++)
		if(fabs(qffPotentialData->quartic[i][i][j][j])>cut) qffPotParameters->numberOf2MR++;

	if(qffPotParameters->numberOf2MR>0) qffPotParameters->qff2MR= malloc(qffPotParameters->numberOf2MR*sizeof(QFFPotnMR));
	m = 0;
	for(i=0;i<nF;i++)
	for(j=0;j<nF;j++)
	{
		if(i==j) continue;
		if(fabs(qffPotentialData->cubic[i][i][j])>cut) 
		{
			qffPotParameters->qff2MR[m].numbers[0] = i;
			qffPotParameters->qff2MR[m].numbers[1] = j;
			qffPotParameters->qff2MR[m].numbers[2] = 2;
			qffPotParameters->qff2MR[m].energy=3*qffPotentialData->cubic[i][i][j]/6.0;
			qffPotParameters->qff2MR[m].grad=2*qffPotParameters->qff2MR[m].energy;
			m++;
		}
	}
	for(i=0;i<nF;i++)
	for(j=0;j<nF;j++)
	{
		if(i==j) continue;
		if(fabs(qffPotentialData->quartic[i][i][i][j])>cut)
		{
			qffPotParameters->qff2MR[m].numbers[0] = i;
			qffPotParameters->qff2MR[m].numbers[1] = j;
			qffPotParameters->qff2MR[m].numbers[2] = 3;
			qffPotParameters->qff2MR[m].energy=4*qffPotentialData->quartic[i][i][i][j]/24.0;
			qffPotParameters->qff2MR[m].grad=3*qffPotParameters->qff2MR[m].energy;
			m++;
		}
	}
	for(i=0;i<nF;i++)
	for(j=0;j<i;j++)
	{
		if(fabs(qffPotentialData->quartic[i][i][j][j])>cut)
		{
			qffPotParameters->qff2MR[m].numbers[0] = i;
			qffPotParameters->qff2MR[m].numbers[1] = j;
			qffPotParameters->qff2MR[m].numbers[2] = -2;
			qffPotParameters->qff2MR[m].energy=6*qffPotentialData->quartic[i][i][j][j]/24.0;
			qffPotParameters->qff2MR[m].grad=2*qffPotParameters->qff2MR[m].energy;
			m++;
		}
	}
	// V3MR
	qffPotParameters->numberOf3MR = 0;
	qffPotParameters->qff3MR = NULL;
	for(i=0;i<nF;i++)
	for(j=0;j<i;j++)
	for(k=0;k<j;k++)
	{
		if(fabs(qffPotentialData->cubic[i][j][k])>cut) qffPotParameters->numberOf3MR++;
	}
	for(i=0;i<nF;i++)
	for(j=0;j<nF;j++)
	for(k=0;k<j;k++)
	{
		if(i==j||i==k) continue;
		if(fabs(qffPotentialData->quartic[i][i][j][k])>cut) qffPotParameters->numberOf3MR++;
	}
	if(qffPotParameters->numberOf3MR>0) qffPotParameters->qff3MR= malloc(qffPotParameters->numberOf3MR*sizeof(QFFPotnMR));
	m = 0;
	for(i=0;i<nF;i++)
	for(j=0;j<i;j++)
	for(k=0;k<j;k++)
	{	
		if(fabs(qffPotentialData->cubic[i][j][k])>cut) 
		{
			qffPotParameters->qff3MR[m].numbers[0] = i;
			qffPotParameters->qff3MR[m].numbers[1] = j;
			qffPotParameters->qff3MR[m].numbers[2] = k;
			qffPotParameters->qff3MR[m].numbers[3] = 1;
			qffPotParameters->qff3MR[m].energy=6*qffPotentialData->cubic[i][j][k]/6.0;
			qffPotParameters->qff3MR[m].grad=qffPotParameters->qff3MR[m].energy;
			m++;
		}
	}
	for(i=0;i<nF;i++)
	for(j=0;j<nF;j++)
	for(k=0;k<j;k++)
	{
		if(i==j||i==k) continue;
		if(fabs(qffPotentialData->quartic[i][i][j][k])>cut)
		{
			qffPotParameters->qff3MR[m].numbers[0] = i;
			qffPotParameters->qff3MR[m].numbers[1] = j;
			qffPotParameters->qff3MR[m].numbers[2] = k;
			qffPotParameters->qff3MR[m].numbers[3] = 2;
			qffPotParameters->qff3MR[m].energy=6*2*qffPotentialData->quartic[i][i][j][k]/24.0;
			qffPotParameters->qff3MR[m].grad=2*qffPotParameters->qff3MR[m].energy;
			m++;
		}
	}
	// free Cubic

	qffPotParameters->numberOf4MR = 0;
	qffPotParameters->qff4MR = NULL;
	for(i=0;i<nF;i++)
	for(j=0;j<i;j++)
	for(k=0;k<j;k++)
	for(l=0;l<k;l++)
	{
		if(fabs(qffPotentialData->quartic[i][j][k][l])>cut) qffPotParameters->numberOf4MR++;
	}

	if(qffPotParameters->numberOf4MR>0) qffPotParameters->qff4MR= malloc(qffPotParameters->numberOf4MR*sizeof(QFFPotnMR));
	m = 0;
	for(i=0;i<nF;i++)
	for(j=0;j<i;j++)
	for(k=0;k<j;k++)
	for(l=0;l<k;l++)
	{
		if(fabs(qffPotentialData->quartic[i][j][k][l])>cut)
		{
			qffPotParameters->qff4MR[m].numbers[0] = i;
			qffPotParameters->qff4MR[m].numbers[1] = j;
			qffPotParameters->qff4MR[m].numbers[2] = k;
			qffPotParameters->qff4MR[m].numbers[3] = l;
			qffPotParameters->qff4MR[m].energy=24*qffPotentialData->quartic[i][j][k][l]/24.0;
			qffPotParameters->qff4MR[m].grad=qffPotParameters->qff4MR[m].energy;
			m++;
		}
	}
	qffPotentialData->qffPotParameters = qffPotParameters;
}
