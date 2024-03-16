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

/* Job.c */

#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <time.h>
#include <locale.h>
#include <string.h>

#include "Job.h"
#include "../Utils/Types.h"
#include "../Utils/Utils.h"
#include "../MolecularMechanics/MolecularMechanicsDlg.h"
#include "../QuantumMechanics/QuantumMechanicsDlg.h"
#include "../PathIntegral/PIMDDlg.h"
#include "../QuarticForceField/QFFnMR.h"
#include "../VPT2/VPT2.h"

typedef struct _JobLink  JobLink;
struct _JobLink
{
        char* type;
	void (*run)(char* inputFileName, char* model);
};
/******************************************************************************************/
Job newJob(char* inputFileName)
{
	Job job;
	job.inputFileName = strdup(inputFileName);
	return job;
}
/******************************************************************************************/
void freeJob(Job* job)
{
	if(job->inputFileName) free(job->inputFileName);
	job->inputFileName = NULL;
}
/******************************************************************************************/
static void runEnergy(char* inputFileName, char* model) 
{
	if(!strcmp(model,"MM")) molecularMechanicsEnergyDlg(inputFileName);
	else quantumMechanicsEnergyDlg(inputFileName);
}
static void runGradient(char* inputFileName, char* model) 
{
	if(!strcmp(model,"MM")) molecularMechanicsGradientDlg(inputFileName);
	else quantumMechanicsGradientDlg(inputFileName);
}
static void runMinimize(char* inputFileName, char* model) 
{
	if(!strcmp(model,"MM")) molecularMechanicsMinimizeDlg(inputFileName);
	else quantumMechanicsMinimizeDlg(inputFileName);
}
static void runOptFreq(char* inputFileName, char* model) 
{
	if(!strcmp(model,"MM")) molecularMechanicsOptFrequenciesDlg(inputFileName);
	else quantumMechanicsOptFrequenciesDlg(inputFileName);
}
static void runGenerateCChemIFilesForFrequencies(char* inputFileName, char* model) 
{
	generateCChemIFilesForFrequenciesDlg(inputFileName, FALSE);
}
static void runComputeFrequenciesFromFiles(char* inputFileName, char* model) 
{
	computeFrequenciesFromFilesDlg(inputFileName, FALSE);
}
static void runGenerateCChemIFilesForFrequenciesOneStep(char* inputFileName, char* model) 
{
	generateCChemIFilesForFrequenciesDlg(inputFileName, TRUE);
}
static void runComputeFrequenciesFromFilesOneStep(char* inputFileName, char* model) 
{
	computeFrequenciesFromFilesDlg(inputFileName, TRUE);
}
static void runGenerateCChemIGradFilesForFrequencies(char* inputFileName, char* model) 
{
	generateCChemIGradFilesForFrequenciesDlg(inputFileName, FALSE);
}
static void runComputeFrequenciesFromGradFiles(char* inputFileName, char* model) 
{
	computeFrequenciesFromGradFilesDlg(inputFileName, FALSE);
}
static void runGenerateCChemIGradFilesForFrequenciesOneStep(char* inputFileName, char* model) 
{
	generateCChemIGradFilesForFrequenciesDlg(inputFileName, TRUE);
}
static void runComputeFrequenciesFromGradFilesOneStep(char* inputFileName, char* model) 
{
	computeFrequenciesFromGradFilesDlg(inputFileName, TRUE);
}
static void runFrequencies(char* inputFileName, char* model) 
{
	if(!strcmp(model,"MM")) molecularMechanicsFrequenciesDlg(inputFileName);
	else quantumMechanicsFrequenciesDlg(inputFileName);
}
static void runComputeQFFFromEnergiesDipolesFile(char* inputFileName, char* model) 
{
	computeQFFFromEnergiesDipolesFile(inputFileName);
}
static void runComputeQFFFromFiles(char* inputFileName, char* model) 
{
	computeQFFFromFiles(inputFileName);
}
static void runGenerateQFFCChemIFilesForFrequencies(char* inputFileName, char* model) 
{
	generateQFFCChemIFilesForFrequencies(inputFileName);
}
static void runPathIntegralMD(char* inputFileName, char* model) 
{
	pathIntegralMDDlg(inputFileName);
}
static void runREMDConfo(char* inputFileName, char* model) 
{
	if(!strcmp(model,"MM")) molecularMechanicsDynamicsREMDConfoDlg(inputFileName);
	else quantumMechanicsREMDConfoDlg(inputFileName);
}
static void runRemoveSimilarConfo(char* inputFileName, char* model) 
{
	quantumMechanicsRemoveSimilarConfoDlg(inputFileName);
}
static void runGAConfo(char* inputFileName, char* model) 
{
	quantumMechanicsGAConfoDlg(inputFileName);
}
static void runMDConfo(char* inputFileName, char* model) 
{
	if(!strcmp(model,"MM")) molecularMechanicsDynamicsConfoDlg(inputFileName);
	else quantumMechanicsMolecularDynamicsConfoDlg(inputFileName);
}
static void runRandomConfo(char* inputFileName, char* model) 
{
	if(!strcmp(model,"MM")) molecularMechanicsRandomConfoDlg(inputFileName);
	else quantumMechanicsRandomConfoDlg(inputFileName);
}
static void runMD(char* inputFileName, char* model) 
{
	if(!strcmp(model,"MM")) molecularMechanicsDynamicsDlg(inputFileName);
	else quantumMechanicsMolecularDynamicsDlg(inputFileName);
}
static void runVPT2(char* inputFileName, char* model) 
{
	vpt2(inputFileName);
}
static void printListOfRunTypes(JobLink jobsList[], int nJobsList)
{
	int i;
	printf("available runType :\n");
	for(i=0;i<nJobsList;i++) printf("%s\n",jobsList[i].type);
}
/******************************************************************************************/
void runJob(Job* job) 
{
	char* inputFileName = job->inputFileName;
	FILE* file = NULL;
	char* runType = NULL;
	char* model = NULL;
	JobLink jobsList[] = {
		{"ENERGY", runEnergy}, 
		{"GRADIENT", runGradient},
		{"OPTIMIZATION", runMinimize},
		{"OPTFREQ", runOptFreq},
		{"GENERATEFILESFORFREQ", runGenerateCChemIFilesForFrequencies},
		{"COMPUTEFREQUENCIESFROMFILES", runComputeFrequenciesFromFiles},
		{"GENERATEFILESONESTEPFORFREQ", runGenerateCChemIFilesForFrequenciesOneStep},
		{"COMPUTEFREQUENCIESONESTEPFROMFILES", runComputeFrequenciesFromFilesOneStep},
		{"GENERATEGRADFILESFORFREQ", runGenerateCChemIGradFilesForFrequencies},
		{"COMPUTEFREQUENCIESFROMGRADFILES", runComputeFrequenciesFromGradFiles},
		{"GENERATEGRADFILESONESTEPFORFREQ", runGenerateCChemIGradFilesForFrequenciesOneStep},
		{"COMPUTEFREQUENCIESONESTEPFROMGRADFILES", runComputeFrequenciesFromGradFilesOneStep},
		{"FREQ", runFrequencies},
		{"COMPUTEQFFNMRFROMENERG", runComputeQFFFromEnergiesDipolesFile},
		{"COMPUTEQFFNMRFROMFILES", runComputeQFFFromFiles},
		{"GENERATEQFFNMRFILES", runGenerateQFFCChemIFilesForFrequencies},
		{"COMPUTEQFF2MRFROMENERG", runComputeQFFFromEnergiesDipolesFile},
		{"COMPUTEQFF2MRFROMFILES", runComputeQFFFromFiles},
		{"GENERATEQFF2MRFILES", runGenerateQFFCChemIFilesForFrequencies},
		{"PIMD", runPathIntegralMD},
		{"REMDCONFO", runREMDConfo},
		{"REMOVESIMILARCONFO", runRemoveSimilarConfo},
		{"MDCONFO", runMDConfo},
		{"RDCONFO", runRandomConfo},
		{"GACONFO", runGAConfo},
		{"MD",  runMD},
		{"VPT2", runVPT2}
	};
	int i;
	int nJobsList = sizeof(jobsList)/sizeof(JobLink);
	if(!inputFileName) exit(1);
	file = fopen(inputFileName, "r");
	if(!file)
	{
		printf("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n");
		printf("Sorry, I cannot open the input file : %s\n",inputFileName);
		printf("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n");
		exit(1);
	}
	if(!readOneString(file,"RunType",&runType)) runType = strdup("ENERGY");
	if(!readOneString(file,"Model",&model)) model = strdup("MM");
	uppercase(runType);
	uppercase(model);
	fclose(file);
	printf("----------------------------------------------------------\n");
	printf("runType=%s\n",runType);
	printf("model=%s\n",model);
	printf("----------------------------------------------------------\n");
	printf("\n");

	if(strstr(runType,"HELP"))printListOfRunTypes(jobsList,nJobsList);
	else
	for(i=0;i<nJobsList;i++) 
	{
		if(strstr(runType,jobsList[i].type))
		{
			jobsList[i].run(inputFileName, model);
			break;
		}
	}
	printf("\n");
}
