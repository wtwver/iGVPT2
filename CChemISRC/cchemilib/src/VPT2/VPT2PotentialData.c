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

/* VPT2PotentialData.c */
#include <math.h>
#include "VPT2PotentialData.h"

static void readData(VPT2PotentialData* vpt2PotentialData, char* fileName);
/**********************************************************************/
static VPT2VModel newModel(VPT2ModelType type, double alpha, double beta)
{
	VPT2VModel model;
	model.type = type;
	model.alphaHDCPT2 = alpha;
	model.betaHDCPT2 = beta;
	return model;
}
/**********************************************************************/
VPT2PotentialData newVPT2PotentialData(int n)
{
	VPT2PotentialData vpt2PotentialData;
	vpt2PotentialData.klass = malloc(sizeof(VPT2PotentialDataClass));
	vpt2PotentialData.klass->readData = readData;
	vpt2PotentialData.nFrequencies = n;
	if(vpt2PotentialData.nFrequencies<=0) vpt2PotentialData.nFrequencies = 0;
	vpt2PotentialData.gradients = NULL;
	vpt2PotentialData.gradients = newVectorDouble(vpt2PotentialData.nFrequencies);
	initVectorDouble(vpt2PotentialData.gradients, vpt2PotentialData.nFrequencies, 0.0);

	vpt2PotentialData.hessian = newMatrixDouble(vpt2PotentialData.nFrequencies,vpt2PotentialData.nFrequencies);
	initMatrixDouble(vpt2PotentialData.hessian, vpt2PotentialData.nFrequencies, vpt2PotentialData.nFrequencies, 0.0);

	vpt2PotentialData.Be = newVectorDouble(3);
	initVectorDouble(vpt2PotentialData.Be, 3, 1.0);

	vpt2PotentialData.cubic = newCubeDouble(vpt2PotentialData.nFrequencies, vpt2PotentialData.nFrequencies, vpt2PotentialData.nFrequencies);
	initCubeDouble(vpt2PotentialData.cubic, vpt2PotentialData.nFrequencies, vpt2PotentialData.nFrequencies, vpt2PotentialData.nFrequencies, 0.0);

	vpt2PotentialData.coriolis = newCubeDouble(3, vpt2PotentialData.nFrequencies, vpt2PotentialData.nFrequencies);
	initCubeDouble(vpt2PotentialData.coriolis, 3, vpt2PotentialData.nFrequencies, vpt2PotentialData.nFrequencies, 0.0);

	vpt2PotentialData.quartic = newQuarticDouble(vpt2PotentialData.nFrequencies, vpt2PotentialData.nFrequencies, vpt2PotentialData.nFrequencies, vpt2PotentialData.nFrequencies);
	initQuarticDouble(vpt2PotentialData.quartic, vpt2PotentialData.nFrequencies, vpt2PotentialData.nFrequencies, vpt2PotentialData.nFrequencies, vpt2PotentialData.nFrequencies, 0.0);
	vpt2PotentialData.model = newModel(MODEL_VPT2,1.0,5e5);

	vpt2PotentialData.maxFrequencyDifferenceFermi=-1;
	vpt2PotentialData.parametersResonance[0]=-1.0;
	vpt2PotentialData.parametersResonance[1]=-1.0;
	vpt2PotentialData.parametersResonance[2]=-1.0;
	return vpt2PotentialData;

}
/*****************************************************************************/
/*
static void printVPT2KHessian(VPT2PotentialData* vpt2PotentialData, char* fileName)
{
	FILE* file =fopen(fileName,"w");
	int n = vpt2PotentialData->nFrequencies;
	double cutoff = 1e-10;
	double** M = vpt2PotentialData->hessian;
	int i;
	for(i = 0;i<n; i++) 
	{
		if(fabs(M[i][i])>=cutoff)
      			fprintf(file,"%d %d %20.10f\n",i+1,i+1,M[i][i]);
	}
	fclose(file);
}
static void printVPT2KCubic(VPT2PotentialData* vpt2PotentialData, char* fileName)
{
	FILE* file =fopen(fileName,"w");
	int n = vpt2PotentialData->nFrequencies;
	double cutoff = 1e-10;
	double*** M = vpt2PotentialData->cubic;
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
static void printVPT2KCoriolis(VPT2PotentialData* vpt2PotentialData, char* fileName)
{
	FILE* file =fopen(fileName,"w");
	int n = vpt2PotentialData->nFrequencies;
	int m = 3;
	double cutoff = 1e-10;
	double*** M = vpt2PotentialData->coriolis;
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
static void printVPT2KRot(VPT2PotentialData* vpt2PotentialData, char* fileName)
{
	FILE* file =fopen(fileName,"w");
	double cutoff = 1e-10;
	double* M = vpt2PotentialData->Be;
	int i;
	for(i = 0;i<3; i++) 
	{
		if(fabs(M[i])>=cutoff)
      			fprintf(file,"%d %20.10f\n",i+1,M[i]);
	}
	fclose(file);
}
static void printVPT2KQuartic(VPT2PotentialData* vpt2PotentialData, char* fileName)
{
	FILE* file =fopen(fileName,"w");
	int i,j,k,l;
	double cutoff = 1e-10;
	double**** C = vpt2PotentialData->quartic;
	int n = vpt2PotentialData->nFrequencies;
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
*/
/*
static void printVPT2KData(VPT2PotentialData* vpt2PotentialData)
{
	printVPT2KHessian(vpt2PotentialData,"data.w");
	printVPT2KCubic(vpt2PotentialData,"data.f3");
	printVPT2KQuartic(vpt2PotentialData,"data.f4");
	printVPT2KRot(vpt2PotentialData,"data.b");
	printVPT2KCoriolis(vpt2PotentialData,"data.z");
}
*/
/*****************************************************************************/
static void printGradients(VPT2PotentialData* vpt2PotentialData)
{
	printf("\nGradients\n");
	printVectorDoubleCutOff(vpt2PotentialData->gradients, vpt2PotentialData->nFrequencies, 1e-10);
	printf("END\n");
}
static void printHessian(VPT2PotentialData* vpt2PotentialData)
{
	printf("\nHessian\n");
	printMatrixDoubleCutOff(vpt2PotentialData->hessian, vpt2PotentialData->nFrequencies, vpt2PotentialData->nFrequencies,1e-10);
	printf("END\n");
}
static void printCubic(VPT2PotentialData* vpt2PotentialData)
{
	printf("\nCubic\n");
	printCubeDoubleCutOff(vpt2PotentialData->cubic, vpt2PotentialData->nFrequencies, vpt2PotentialData->nFrequencies, vpt2PotentialData->nFrequencies,1e-10);
	printf("END\n");
}
static void printCoriolis(VPT2PotentialData* vpt2PotentialData)
{
	printf("\nCoriolis\n");
	printCubeDoubleCutOff(vpt2PotentialData->coriolis, 3, vpt2PotentialData->nFrequencies, vpt2PotentialData->nFrequencies,1e-10);
	printf("END\n");
}
static void printRot(VPT2PotentialData* vpt2PotentialData)
{
	printf("\nRotation constants\n");
	printVectorDoubleCutOff(vpt2PotentialData->Be, 3,1e-10);
	printf("END\n");
}
static void printQuartic(VPT2PotentialData* vpt2PotentialData)
{
   printf("\nQuartic\n");
   printQuarticDoubleCutOff(vpt2PotentialData->quartic, vpt2PotentialData->nFrequencies, vpt2PotentialData->nFrequencies, vpt2PotentialData->nFrequencies, vpt2PotentialData->nFrequencies,1e-10);
	printf("END\n");
}
static void printData(VPT2PotentialData* vpt2PotentialData)
{
	printGradients(vpt2PotentialData);
	printHessian(vpt2PotentialData);
	printCubic(vpt2PotentialData);
	printQuartic(vpt2PotentialData);
	printCoriolis(vpt2PotentialData);
	printRot(vpt2PotentialData);
}
/*****************************************************************************/
static void readResonancesData(VPT2PotentialData* vpt2PotentialData, char* inputFileName)
{
	FILE* inputFile = fopen(inputFileName,"rb");
	readOneReal(inputFile,"maxFrequencyDifferenceFermi",&vpt2PotentialData->maxFrequencyDifferenceFermi);
	readOneReal(inputFile,"MartinCutOff1",&vpt2PotentialData->parametersResonance[0]);
	readOneReal(inputFile,"MartinCutOff2",&vpt2PotentialData->parametersResonance[1]);
	readOneReal(inputFile,"ZCutOff",&vpt2PotentialData->parametersResonance[2]);
	//printf("Zcutoff=%f\n",vpt2PotentialData->parametersResonance[2]);
	fclose(inputFile);
}
/*****************************************************************************/
static void readData(VPT2PotentialData* vpt2PotentialData, char* inputFileName)
{
	FILE* inputFile;
	char* tmp = NULL;
        inputFile = fopen(inputFileName,"rb");
	if(!inputFile)
	{
		fprintf(stderr, "==========================================================\n");
		fprintf(stderr, "Sorry, I cannot opent the %s file\n", inputFileName);
		fprintf(stderr, "==========================================================\n");
		exit(1);
	}
	readOneInt(inputFile,"nFrequencies",&vpt2PotentialData->nFrequencies);
	fprintf(stdout, "nFrequencies=%d\n", vpt2PotentialData->nFrequencies);
	*vpt2PotentialData = newVPT2PotentialData(vpt2PotentialData->nFrequencies);
	readOneReal(inputFile,"alphaHDCPT2",&vpt2PotentialData->model.alphaHDCPT2);
	readOneReal(inputFile,"betaHDCPT2",&vpt2PotentialData->model.betaHDCPT2);
	
	readOneString(inputFile,"VPT2Model",&tmp);
	if(tmp && mystrcasestr(tmp,"NONE")) vpt2PotentialData->model.type = MODEL_VPT2;
	if(tmp && mystrcasestr(tmp,"DCPT2")) vpt2PotentialData->model.type = MODEL_DCPT2;
	if(tmp && mystrcasestr(tmp,"HDCPT2")) vpt2PotentialData->model.type = MODEL_HDCPT2;
	if(tmp && mystrcasestr(tmp,"GVPT2")) vpt2PotentialData->model.type = MODEL_GVPT2;
	if(tmp && mystrcasestr(tmp,"VPT2+K")) vpt2PotentialData->model.type = MODEL_VPT2K;
	if(tmp) free(tmp);
	if(vpt2PotentialData->model.type==MODEL_GVPT2) printf("VPT2 model : GVPT2\n");
	if(vpt2PotentialData->model.type==MODEL_VPT2K) printf("VPT2 model : VPT2+K\n");
	if(vpt2PotentialData->model.type==MODEL_VPT2) printf("VPT2 model : VPT2\n");
	if(vpt2PotentialData->model.type==MODEL_DCPT2) printf("VPT2 model : DCPT2\n");
	if(vpt2PotentialData->model.type==MODEL_HDCPT2) printf("VPT2 model : HDCPT2\n");

	
	readVectorReal(inputFile,"Gradients",vpt2PotentialData->nFrequencies, vpt2PotentialData->gradients);/* vpt2PotentialData.gradients presently not used, min surface */
	readMatrixReal(inputFile,"Hessian",vpt2PotentialData->nFrequencies, vpt2PotentialData->nFrequencies, vpt2PotentialData->hessian);
	readCubeReal(inputFile,"Cubic",vpt2PotentialData->nFrequencies, vpt2PotentialData->nFrequencies, vpt2PotentialData->nFrequencies, vpt2PotentialData->cubic);
	readQuarticReal(inputFile,"Quartic",vpt2PotentialData->nFrequencies, vpt2PotentialData->nFrequencies, vpt2PotentialData->nFrequencies, vpt2PotentialData->nFrequencies, vpt2PotentialData->quartic);
	
	readVectorReal(inputFile,"Rotational",3, vpt2PotentialData->Be);
	readCubeReal(inputFile,"Coriolis",3, vpt2PotentialData->nFrequencies, vpt2PotentialData->nFrequencies, vpt2PotentialData->coriolis);
	readResonancesData(vpt2PotentialData, inputFileName);
	printData(vpt2PotentialData);
	//printVPT2KData(vpt2PotentialData);
	fclose(inputFile);

	
}
/**********************************************************************/
