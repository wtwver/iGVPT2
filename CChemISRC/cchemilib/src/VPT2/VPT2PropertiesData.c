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

/* VPT2Properties.c */

#include <math.h>
#include "VPT2PropertiesData.h"

//static boolean printMax = FALSE;

static void readData(VPT2PropertiesData* vpt2PropertiesData, char* fileName);
static void printData(VPT2PropertiesData* vpt2PropertiesData);
/**********************************************************************/
static VPT2PropModel newModel(VPT2PropModelType type, double alpha, double beta)
{
	VPT2PropModel model;
	model.type = type;
	model.alphaHDCPT2 = alpha;
	model.betaHDCPT2 = beta;
	return model;
}
/**********************************************************************/
/* PCCP, 2014, 16, 1759-1787, page 1763-4 */
/* CPL, 496 (2010) 157–161 */
/**********************************************************************/
VPT2PropertiesData newVPT2PropertiesData(int nDim, int n)
{
	VPT2PropertiesData vpt2PropertiesData;
	vpt2PropertiesData.klass = malloc(sizeof(VPT2PropertiesDataClass));
	vpt2PropertiesData.klass->readData = readData;
	vpt2PropertiesData.nFrequencies = n;
	vpt2PropertiesData.nDim = nDim;
	if(vpt2PropertiesData.nFrequencies<=0) vpt2PropertiesData.nFrequencies = 0;
	if(vpt2PropertiesData.nDim<=0) vpt2PropertiesData.nDim = 0;
	if(vpt2PropertiesData.nDim>6) vpt2PropertiesData.nDim = 6;
	vpt2PropertiesData.first = NULL;
	vpt2PropertiesData.first = newMatrixDouble(vpt2PropertiesData.nDim,vpt2PropertiesData.nFrequencies);
	initMatrixDouble(vpt2PropertiesData.first, vpt2PropertiesData.nDim, vpt2PropertiesData.nFrequencies, 0.0);

	vpt2PropertiesData.second = newCubeDouble(vpt2PropertiesData.nDim, vpt2PropertiesData.nFrequencies,vpt2PropertiesData.nFrequencies);
	initCubeDouble(vpt2PropertiesData.second, vpt2PropertiesData.nDim, vpt2PropertiesData.nFrequencies, vpt2PropertiesData.nFrequencies, 0.0);

	vpt2PropertiesData.cubic = newQuarticDouble(vpt2PropertiesData.nDim, vpt2PropertiesData.nFrequencies, vpt2PropertiesData.nFrequencies, vpt2PropertiesData.nFrequencies);
	initQuarticDouble(vpt2PropertiesData.cubic, vpt2PropertiesData.nDim, vpt2PropertiesData.nFrequencies, vpt2PropertiesData.nFrequencies, vpt2PropertiesData.nFrequencies, 0.0);
	vpt2PropertiesData.model = newModel(MODEL_PROP_VPT2,1.0,5e5);
	return vpt2PropertiesData;
}
/*****************************************************************************/
static void printFirst(VPT2PropertiesData* vpt2PropertiesData)
{
	printf("\nFirst derivatives\n");
	printMatrixDoubleCutOff(vpt2PropertiesData->first, vpt2PropertiesData->nDim, vpt2PropertiesData->nFrequencies, 1e-10);
	printf("END\n\n");
}
static void printSecond(VPT2PropertiesData* vpt2PropertiesData)
{
	printf("\nSecond derivatives\n");
	printCubeDoubleCutOff(vpt2PropertiesData->second, vpt2PropertiesData->nDim, vpt2PropertiesData->nFrequencies, vpt2PropertiesData->nFrequencies,1e-10);
	printf("END\n\n");
}
static void printCubic(VPT2PropertiesData* vpt2PropertiesData)
{
	printf("\nCubic derivatives\n");
	printQuarticDoubleCutOff(vpt2PropertiesData->cubic, vpt2PropertiesData->nDim, vpt2PropertiesData->nFrequencies, vpt2PropertiesData->nFrequencies, vpt2PropertiesData->nFrequencies,1e-10);
	printf("END\n\n");
}
static void printData(VPT2PropertiesData* vpt2PropertiesData)
{
	printFirst(vpt2PropertiesData);
	printSecond(vpt2PropertiesData);
	printCubic(vpt2PropertiesData);
}
/*****************************************************************************/
static void readData(VPT2PropertiesData* vpt2PropertiesData, char* inputFileName)
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
	readOneInt(inputFile,"nFrequencies",&vpt2PropertiesData->nFrequencies);
	fprintf(stdout, "nFrequencies=%d\n", vpt2PropertiesData->nFrequencies);
	readOneInt(inputFile,"nDim",&vpt2PropertiesData->nDim);
	fprintf(stdout, "nDim=%d\n", vpt2PropertiesData->nDim);
	if(vpt2PropertiesData->nDim<=0)
	{
		fprintf(stderr, "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n");
		fprintf(stderr, "I cannot read properties vpt2PropertiesData from input file : nDim=%d\n", vpt2PropertiesData->nDim);
		fprintf(stderr, "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n\n");
		exit(1);
	}
	*vpt2PropertiesData = newVPT2PropertiesData(vpt2PropertiesData->nDim, vpt2PropertiesData->nFrequencies);

	readMatrixReal(inputFile,"First Derivatives",vpt2PropertiesData->nDim, vpt2PropertiesData->nFrequencies, vpt2PropertiesData->first);
	readCubeReal(inputFile,"Second Derivatives",vpt2PropertiesData->nDim, vpt2PropertiesData->nFrequencies, vpt2PropertiesData->nFrequencies, vpt2PropertiesData->second);
	readQuarticReal(inputFile,"Cubic Derivatives", vpt2PropertiesData->nDim, vpt2PropertiesData->nFrequencies, vpt2PropertiesData->nFrequencies, vpt2PropertiesData->nFrequencies, vpt2PropertiesData->cubic);

	readOneReal(inputFile,"alphaPropHDCPT2",&vpt2PropertiesData->model.alphaHDCPT2);
	readOneReal(inputFile,"betaPropHDCPT2",&vpt2PropertiesData->model.betaHDCPT2);
	
	readOneString(inputFile,"PropModel",&tmp);
	if(tmp && mystrcasestr(tmp,"NONE")) vpt2PropertiesData->model.type = MODEL_PROP_VPT2;
	if(tmp && mystrcasestr(tmp,"VPT2")) vpt2PropertiesData->model.type = MODEL_PROP_VPT2;
	if(tmp && mystrcasestr(tmp,"DCPT2")) vpt2PropertiesData->model.type = MODEL_PROP_DCPT2;
	if(tmp && mystrcasestr(tmp,"HDCPT2")) vpt2PropertiesData->model.type = MODEL_PROP_HDCPT2;
	if(tmp && mystrcasestr(tmp,"VPT2+K")) vpt2PropertiesData->model.type = MODEL_PROP_VPT2K;
	if(tmp && mystrcasestr(tmp,"GVPT2")) vpt2PropertiesData->model.type = MODEL_PROP_GVPT2;
	if(tmp && mystrcasestr(tmp,"GVPT2S")) vpt2PropertiesData->model.type = MODEL_PROP_GVPT2S;
	if(tmp) free(tmp);
	if(vpt2PropertiesData->model.type==MODEL_PROP_VPT2) printf("VPT2 Properties model : VPT2\n");
	if(vpt2PropertiesData->model.type==MODEL_PROP_DCPT2) printf("VPT2 Properties model : DCPT2\n");
	if(vpt2PropertiesData->model.type==MODEL_PROP_HDCPT2) printf("VPT2 Propertiesmodel : HDCPT2\n");
	if(vpt2PropertiesData->model.type==MODEL_PROP_VPT2K) printf("VPT2 Properties model : VPT2+K\n");
	if(vpt2PropertiesData->model.type==MODEL_PROP_GVPT2) printf("VPT2 Properties model : GVPT2 : Gaussian method\n");
	if(vpt2PropertiesData->model.type==MODEL_PROP_GVPT2S) printf("VPT2 Properties model : GVPT2 : Gaussian method but sum on all basis set\n");

	printData(vpt2PropertiesData);

}
