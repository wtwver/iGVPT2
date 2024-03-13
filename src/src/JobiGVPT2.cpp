/*****************************************************************************************
 iGVPT2 is a program for computing anharmonic corrections to vibration frequencies, 
 based on force field expansion of the potential energy surface in normal mode coordinates.
 iGVPT2 supports several computation chemistry packages(CCP).

 Copyright (C) 2020 Abdulrahman Allouche (University Lyon 1)

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
*****************************************************************************************/

/* JobiGVPT2.cpp */

using namespace std;

#include <cstdlib>
#include <cstdio>
#include <unistd.h>
#include <ctime>
#include <locale.h>
#include <string>
#include <cstring>

#include "JobiGVPT2.h"
#include "OBGVPT2.h"
#include "N2P2GVPT2.h"
#include "TMGVPT2.h"

#ifdef __cplusplus
extern "C" {
#include <Utils/Types.h>
#include <Utils/Utils.h>
#include <JobControl/Job.h>
#include <QuarticForceField/QFFnMR.h>
#include <MolecularMechanics/MolecularMechanicsDlg.h>
#include <QuantumMechanics/QuantumMechanicsDlg.h>
#include <PathIntegral/PIMDDlg.h>
#include <VPT2/VPT2.h>
#include <Utils/Constants.h>
}
#endif

struct JobiGVPT2Link
{
        string type;
	void (*run)(char* inputFileName, char* model);
};
/******************************************************************************************/
JobiGVPT2 newJobiGVPT2(char* inputFileName)
{
	JobiGVPT2 job;
	job.inputFileName = strdup(inputFileName);
	return job;
}
/******************************************************************************************/
void freeJobiGVPT2(JobiGVPT2* job)
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
static void runVPT2(char* inputFileName, char* model) 
{
	vpt2(inputFileName);
}
static void runOBGVPT2(char* inputFileName, char* model) 
{
	obgvpt2(inputFileName);
}
static void runOBSelectModes(char* inputFileName, char* model) 
{
	selectModesForTargetModes(inputFileName);
}
static void runN2P2GVPT2(char* inputFileName, char* model) 
{
	n2p2GVPT2(inputFileName);
}
static void runTMGVPT2(char* inputFileName, char* model) 
{
	tmGVPT2(inputFileName);
}
static void printListOfRunTypes(JobiGVPT2Link jobsList[], int nJobiGVPT2sList)
{
	int i;
	printf("available runType :\n");
	for(i=0;i<nJobiGVPT2sList;i++) printf("%s\n",jobsList[i].type.c_str());
}
/******************************************************************************************/
void runJobiGVPT2(JobiGVPT2* job) 
{
	char* inputFileName = job->inputFileName;
	FILE* file = NULL;
	char* runType = NULL;
	char* model = NULL;
	JobiGVPT2Link jobsList[] = {
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
		{"HYBRIDMM", runOBGVPT2},
		{"OBGVPT2", runOBGVPT2},
		{"HDNNGVPT2", runN2P2GVPT2},
		{"N2P2GVPT2", runN2P2GVPT2},
		{"TMGVPT2", runTMGVPT2},
		{"TensorMolGVPT2", runTMGVPT2},
		{"SELECTMODES", runOBSelectModes},
		{"VPT2", runVPT2}
	};
	int i;
	int nJobiGVPT2sList = sizeof(jobsList)/sizeof(JobiGVPT2Link);
	if(!inputFileName) exit(1);
	file = fopen(inputFileName, "r");
	if(!file)
	{
		printf("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n");
		printf("Sorry, I cannot open the input file : %s\n",inputFileName);
		printf("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n");
		exit(1);
	}
	if(!readOneString(file,strdup(string("RunType").c_str()),&runType)) runType = strdup(string("ENERGY").c_str());
	if(!readOneString(file,strdup(string("Model").c_str()),&model)) model = strdup(string("MM").c_str());
	uppercase(runType);
	uppercase(model);
	fclose(file);
	printf("----------------------------------------------------------\n");
	printf("runType=%s\n",runType);
	printf("model=%s\n",model);
	printf("----------------------------------------------------------\n");
	printf("\n");

	if(strstr(runType,"HELP"))printListOfRunTypes(jobsList,nJobiGVPT2sList);
	else
	for(i=0;i<nJobiGVPT2sList;i++) 
	{
		if(strstr(runType,jobsList[i].type.c_str()))
		{
			jobsList[i].run(inputFileName, model);
			break;
		}
	}
	printf("\n");
}
