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

/* QFFProperties.c */
#include <math.h>
#include "../QFFPot/QFFModel.h"
#include "../Utils/Utils.h"
#include "../Utils/AtomsProp.h"
#include "../Utils/Constants.h"
#include "../Utils/QL.h"

//static boolean printMax = FALSE;

static void readData(QFFPropertiesData* qffPropertiesData, char* fileName);
static void convertToAU(QFFPropertiesData* qffPropertiesData);
static void convertToAU2(QFFPropertiesData* qffPropertiesData);
static void printData(QFFPropertiesData* qffPropertiesData);
static void freeQFFPropertiesData(QFFPropertiesData* qffPropertiesData);
/**********************************************************************/
static void freeQFFPropertiesData(QFFPropertiesData* qffPropertiesData)
{
	int nF = qffPropertiesData->nFrequencies;
	int nDim = qffPropertiesData->nDim;

	freeVectorDouble(&qffPropertiesData->zero);
	freeMatrixDouble(&qffPropertiesData->first, nDim);
	freeCubeDouble(&qffPropertiesData->second, nDim, nF);
	freeQuarticDouble(&qffPropertiesData->cubic, nDim, nF, nF);
}
/**********************************************************************/
/**********************************************************************/
/* PCCP, 2014, 16, 1759-1787, page 1763-4 */
/* CPL, 496 (2010) 157–161 */
/**********************************************************************/
QFFPropertiesData newQFFPropertiesData(int nDim, int n)
{
	QFFPropertiesData qffPropertiesData;
	qffPropertiesData.klass = malloc(sizeof(QFFPropertiesDataClass));
	qffPropertiesData.klass->readData = readData;
	qffPropertiesData.klass->convertToAU = convertToAU;
	qffPropertiesData.klass->convertToAU2 = convertToAU2;
	qffPropertiesData.klass->free = freeQFFPropertiesData;
	qffPropertiesData.klass->print = printData;
	qffPropertiesData.nFrequencies = n;
	qffPropertiesData.nDim = nDim;
	if(qffPropertiesData.nFrequencies<=0) qffPropertiesData.nFrequencies = 0;
	if(qffPropertiesData.nDim<=0) qffPropertiesData.nDim = 0;
	if(qffPropertiesData.nDim>6) qffPropertiesData.nDim = 6;

	qffPropertiesData.first = NULL;
	qffPropertiesData.first = newMatrixDouble(qffPropertiesData.nDim,qffPropertiesData.nFrequencies);
	initMatrixDouble(qffPropertiesData.first, qffPropertiesData.nDim, qffPropertiesData.nFrequencies, 0.0);

        qffPropertiesData.zero = newVectorDouble(qffPropertiesData.nDim);
        initVectorDouble(qffPropertiesData.zero, qffPropertiesData.nDim, 0.0);

	qffPropertiesData.second = newCubeDouble(qffPropertiesData.nDim, qffPropertiesData.nFrequencies,qffPropertiesData.nFrequencies);
	initCubeDouble(qffPropertiesData.second, qffPropertiesData.nDim, qffPropertiesData.nFrequencies, qffPropertiesData.nFrequencies, 0.0);

	qffPropertiesData.cubic = newQuarticDouble(qffPropertiesData.nDim, qffPropertiesData.nFrequencies, qffPropertiesData.nFrequencies, qffPropertiesData.nFrequencies);
	initQuarticDouble(qffPropertiesData.cubic, qffPropertiesData.nDim, qffPropertiesData.nFrequencies, qffPropertiesData.nFrequencies, qffPropertiesData.nFrequencies, 0.0);
	return qffPropertiesData;
}
/*****************************************************************************/
static void printFirst(QFFPropertiesData* qffPropertiesData)
{
	printf("\nFirst derivatives\n");
	printMatrixDoubleCutOff(qffPropertiesData->first, qffPropertiesData->nDim, qffPropertiesData->nFrequencies, 1e-10);
	printf("END\n\n");
}
static void printSecond(QFFPropertiesData* qffPropertiesData)
{
	printf("\nSecond derivatives\n");
	printCubeDoubleCutOff(qffPropertiesData->second, qffPropertiesData->nDim, qffPropertiesData->nFrequencies, qffPropertiesData->nFrequencies,1e-10);
	printf("END\n\n");
}
static void printCubic(QFFPropertiesData* qffPropertiesData)
{
	printf("\nCubic derivatives\n");
	printQuarticDoubleCutOff(qffPropertiesData->cubic, qffPropertiesData->nDim, qffPropertiesData->nFrequencies, qffPropertiesData->nFrequencies, qffPropertiesData->nFrequencies,1e-10);
	printf("END\n\n");
}
static void printData(QFFPropertiesData* qffPropertiesData)
{
	printFirst(qffPropertiesData);
	printSecond(qffPropertiesData);
	printCubic(qffPropertiesData);
}
/*****************************************************************************/
static void convertToAU(QFFPropertiesData* qffPropertiesData)
{
	double conv = sqrt(1/AUTOCM1);
	double conv1 = conv;
	double conv2 = conv*conv;
	double conv3 = conv*conv*conv;
	int i,j,k,c;
	int nF =  qffPropertiesData->nFrequencies;

	for(i=0;i<nF;i++)
	for(c=0;c<3;c++) qffPropertiesData->first[c][i]*=conv1;

	for(i=0;i<nF;i++)
	for(j=0;j<nF;j++)
	for(c=0;c<3;c++) qffPropertiesData->second[c][i][j]*=conv2;

	for(i=0;i<nF;i++)
	for(j=0;j<nF;j++)
	for(k=0;k<nF;k++)
	for(c=0;c<3;c++) qffPropertiesData->cubic[c][i][j][k]*=conv3;

}
/*****************************************************************************/
static void convertToAU2(QFFPropertiesData* qffPropertiesData)
{
	double conv = sqrt(1/AUTOCM1);
	double conv1 = conv;
	double conv2 = conv*conv;
	double conv3 = conv*conv*conv;
       	double convAMUToAU = sqrt(AMUTOAU);
	int i,j,k,c;
	int nF =  qffPropertiesData->nFrequencies;

        conv1 *= convAMUToAU;
        conv2 *= convAMUToAU*convAMUToAU;
        conv3 *= convAMUToAU*convAMUToAU*convAMUToAU;

	for(i=0;i<nF;i++)
	for(c=0;c<3;c++) qffPropertiesData->first[c][i]*=conv1;

	for(i=0;i<nF;i++)
	for(j=0;j<nF;j++)
	for(c=0;c<3;c++) qffPropertiesData->second[c][i][j]*=conv2;

	for(i=0;i<nF;i++)
	for(j=0;j<nF;j++)
	for(k=0;k<nF;k++)
	for(c=0;c<3;c++) qffPropertiesData->cubic[c][i][j][k]*=conv3;

}
/*****************************************************************************/
static void readData(QFFPropertiesData* qffPropertiesData, char* inputFileName)
{
	FILE* inputFile;
        inputFile = fopen(inputFileName,"rb");
	if(!inputFile)
	{
		fprintf(stderr, "==========================================================\n");
		fprintf(stderr, "Sorry, I cannot opent the %s file\n", inputFileName);
		fprintf(stderr, "==========================================================\n");
		exit(1);
	}
	readOneInt(inputFile,"nFrequencies",&qffPropertiesData->nFrequencies);
//	fprintf(stdout, "nFrequencies=%d\n", qffPropertiesData->nFrequencies);
	readOneInt(inputFile,"nDim",&qffPropertiesData->nDim);
//	fprintf(stdout, "nDim=%d\n", qffPropertiesData->nDim);
	if(qffPropertiesData->nDim<=0)
	{
		fprintf(stderr, "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n");
		fprintf(stderr, "I cannot read properties qffPropertiesData from input file : nDim=%d\n", qffPropertiesData->nDim);
		fprintf(stderr, "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n\n");
		exit(1);
	}
	*qffPropertiesData = newQFFPropertiesData(qffPropertiesData->nDim, qffPropertiesData->nFrequencies);

	readVectorReal(inputFile,"Zero",qffPropertiesData->nDim, qffPropertiesData->zero);
	readMatrixReal(inputFile,"First Derivatives",qffPropertiesData->nDim, qffPropertiesData->nFrequencies, qffPropertiesData->first);
	readCubeReal(inputFile,"Second Derivatives",qffPropertiesData->nDim, qffPropertiesData->nFrequencies, qffPropertiesData->nFrequencies, qffPropertiesData->second);
	readQuarticReal(inputFile,"Cubic Derivatives", qffPropertiesData->nDim, qffPropertiesData->nFrequencies, qffPropertiesData->nFrequencies, qffPropertiesData->nFrequencies, qffPropertiesData->cubic);

	//printData(qffPropertiesData);

}
