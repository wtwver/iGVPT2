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

/* MolecularDynamics.c  */
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <time.h>

#ifdef ENABLE_MPI
#include <mpi.h>
#endif

#include "../Utils/Utils.h"
#include "../Utils/AtomsProp.h"
#include "../Utils/Constants.h"
#include "../Utils/Timer.h"
#include "../Utils/Types.h"
#include "../Molecule/Molecule.h"
#include "../MolecularMechanics/MolecularDynamics.h"
#include "../MolecularMechanics/MolecularMechanics.h"
#ifdef ENABLE_CL
#include "../Utils/CLProp.h"
#include "../MolecularMechanics/MolecularDynamicsCL.h"
#endif


/*********************************************************************************/
static void initMD(MolecularDynamics* molecularDynamics, double temperature, double stepSize, MDIntegratorType integratorType, MDThermostatType thermostat, double friction, double omegaMax, int Nf, double collide, double qNH, char* fileNameTraj, char* fileNameProp, int numberOfRunSteps, int index);
static void berendsen(MolecularDynamics* molecularDynamics);
static void scaleV(MolecularDynamics* molecularDynamics);
static void andersen(MolecularDynamics* molecularDynamics);
static void bussi(MolecularDynamics* molecularDynamics);
static void nose_hoover(MolecularDynamics* molecularDynamics);
static void rescaleVelocities(MolecularDynamics* molecularDynamics);
static void computeEnergies(MolecularDynamics* molecularDynamics);
static void applyLangevin(MolecularDynamics* molecularDynamics);
static void updateLangevin(MolecularDynamics* molecularDynamics);
static void resetLangevin(MolecularDynamics* molecularDynamics);
static void applyQTB(MolecularDynamics* molecularDynamics);
static void updateQTB(MolecularDynamics* molecularDynamics);
static void resetQTB(MolecularDynamics* molecularDynamics);
static void applyOneStep(MolecularDynamics* molecularDynamics, int iStep);
static void applyThermostat(MolecularDynamics* molecularDynamics);
static void applyVerlet(MolecularDynamics* molecularDynamics);
static void applyMartynaTuckerman(MolecularDynamics* molecularDynamics);
static void applyBeeman(MolecularDynamics* molecularDynamics);
static void applyStochastic(MolecularDynamics* molecularDynamics);
static void newProperties(MolecularDynamics* molecularDynamics, char* comments);
static void saveProperties(MolecularDynamics* molecularDynamics, int iStep0, int iStep, char* comments);
static void saveTrajectory(MolecularDynamics* molecularDynamics, int iStep);
static double getEKin(MolecularDynamics* molecularDynamics);
static double getKelvin(MolecularDynamics* molecularDynamics);
static void removeTranslationAndRotation(MolecularDynamics* molecularDynamics);
#ifdef ENABLE_CL
static void scaleCLVelocities(MolecularDynamics* molecularDynamics, double scale);
#endif

/**********************************************************************/
#ifdef ENABLE_CL
void resetRandomNumbers(MolecularDynamics* molecularDynamics)
{
	size_t global =  molecularDynamics->nRandomsCL;
	cl_int err;
	CLProp clProp = getCLProp();
	err = clEnqueueNDRangeKernel(clProp.command_queue, molecularDynamics->generateRandomNumbers, 1, NULL, &global, NULL, 0, NULL, NULL);
	if(err != CL_SUCCESS)
	{
		printErrorCLRun(err);
       		printf("I cannot execute generateRandomNumbers\n");
		exit(1);
	}
/*
	clFinish(clProp.command_queue);
	{
		int i;
		cl_float4 r[molecularDynamics->forceField->molecule.nAtoms*2];
		clEnqueueReadBuffer(clProp.command_queue, molecularDynamics->randomsCL, CL_TRUE, 0, sizeof(cl_float4)*molecularDynamics->forceField->molecule.nAtoms*2, r, 0, NULL, NULL);
		for(i=0;i<molecularDynamics->forceField->molecule.nAtoms*2;i++)
		printf("rrrrrrrrrrr=%f %f %f %f\n",r[i].s[0], r[i].s[1], r[i].s[2], r[i].s[3]);
	}
*/
}
void initRandomNumberGenerator(MolecularDynamics* molecularDynamics, unsigned int randomNumberSeed)
{
	int n = molecularDynamics->forceField->molecule.nAtoms*2;
	cl_uint4* seed = malloc(sizeof(cl_uint4)*n);
	unsigned int r = randomNumberSeed;
	cl_int err;
	CLProp clProp = getCLProp();
	int i;
	for (i = 0; i < n; i++)
	{
        	seed[i].x = r = (1664525*r + 1013904223) & 0xFFFFFFFF;
        	seed[i].y = r = (1664525*r + 1013904223) & 0xFFFFFFFF;
        	seed[i].z = r = (1664525*r + 1013904223) & 0xFFFFFFFF;
        	seed[i].w = r = (1664525*r + 1013904223) & 0xFFFFFFFF;
	}
	molecularDynamics->seedCL = clCreateBuffer(clProp.context, CL_MEM_READ_WRITE, sizeof(cl_uint4) * n, NULL, NULL);
	molecularDynamics->randomsCL = clCreateBuffer(clProp.context, CL_MEM_READ_WRITE, sizeof(cl_float4) * n, NULL, NULL);
	molecularDynamics->nRandomsCL = n;
	clEnqueueWriteBuffer(clProp.command_queue, molecularDynamics->seedCL, CL_TRUE, 0, sizeof(cl_uint4) * n, seed, 0, NULL, NULL);
	free(seed);

	molecularDynamics->generateRandomNumbers = clCreateKernel(molecularDynamics->programMD, "generateRandomNumbers", &err);
	clSetKernelArg(molecularDynamics->generateRandomNumbers, 0, sizeof(cl_mem),   &molecularDynamics->randomsCL);
	clSetKernelArg(molecularDynamics->generateRandomNumbers, 1, sizeof(cl_mem),   &molecularDynamics->seedCL);
	clSetKernelArg(molecularDynamics->generateRandomNumbers, 2, sizeof(cl_int), &n);

	resetRandomNumbers(molecularDynamics);
}
#endif
/**********************************************************************************************************/
// return # temperatures
static int getNumExchangeReplica(double* energies, double* temperatures, int nTemperatures, int* numMDForTemperatures)
{
	int it = rand()%(nTemperatures-1);
	int n = numMDForTemperatures[it];
	int next = numMDForTemperatures[it+1];
	double delta = (1.0/(Kb*temperatures[n])- 1.0/(Kb*temperatures[next]))*(energies[n]-energies[next]);
	boolean exchange= FALSE;
	printf("n=%d next=%d delta = %f\n",n,next, delta);
	if(delta<0) exchange= TRUE;
	else
	{
		if(rand()/(double)RAND_MAX>exp(-delta)) exchange = TRUE;
	}
	if(exchange) return it;
	else return -1;
}
/*****************************************************************************************************************/
static void applyExchangeReplica(MolecularDynamics* md, int nTemperatures, int n, int next)
{
	double r;
	int i,j;
	double t;
	FILE* f;
	r = sqrt(md[next].temperature/md[n].temperature);
	/* exchange temperatures */
	t = md[n].temperature;
	md[n].temperature = md[next].temperature;
	md[next].temperature = t;
	/* exchange filesoutput */
	f = md[n].fileProp;
	md[n].fileProp = md[next].fileProp;
	md[next].fileProp = f;

	f = md[n].fileTraj;
	md[n].fileTraj = md[next].fileTraj;
	md[next].fileTraj = f;

	/* scale velocity */
	for ( i = 0; i < md[n].numberOfAtoms; i++)
	{
		for ( j = 0; j < 3; j++)
		{
			 md[n].forceField->molecule.atoms[i].velocity[j] *= r;
			 md[next].forceField->molecule.atoms[i].velocity[j] /= r;
		}
	}
#ifdef ENABLE_CL
	scaleCLVelocities(&md[n], r);
	scaleCLVelocities(&md[next], 1.0/r);
#endif
	removeTranslationAndRotation(&md[n]);
	removeTranslationAndRotation(&md[next]);
	resetQTB(&md[n]);
	resetQTB(&md[next]);
	resetLangevin(&md[n]);
	resetLangevin(&md[next]);
}
/*****************************************************************************************************************/
static void changeOneReplica(MolecularDynamics* md, int n, double newTemperature, char* fileNamePrefixProp,  char* fileNamePrefixTraj)
{
	double r;
	int i,j;
	char*  fileNameProp = NULL;
	char*  fileNameTraj = NULL;
	r = sqrt(newTemperature/md[n].temperature);
	/* change temperature */
//	printf("Old temp = %f newTemp = %f\n", md[n].temperature, newTemperature);
	md[n].temperature = newTemperature;
	/* close old file and open file for new temperture*/
	if(fileNamePrefixProp) fileNameProp = strdup_printf("%s%0.0f.txt", fileNamePrefixProp, newTemperature);
	if(fileNamePrefixTraj)  fileNameTraj = strdup_printf("%s%0.0f.gab", fileNamePrefixTraj, newTemperature);
	if(fileNameProp) md[n].fileProp = fopen(fileNameProp,"a");
	if(fileNameTraj) md[n].fileTraj = fopen(fileNameTraj,"a");

	/* scale velocity */
	for ( i = 0; i < md[n].numberOfAtoms; i++)
	{
		for ( j = 0; j < 3; j++)
			 md[n].forceField->molecule.atoms[i].velocity[j] *= r;
	}
#ifdef ENABLE_CL
	scaleCLVelocities(&md[n], r);
#endif
	removeTranslationAndRotation(&md[n]);
	if(fileNameProp) free(fileNameProp);
	if(fileNameTraj) free(fileNameTraj);
}
/**********************************************************************/
/*
static void exchangeReplicaLocalAll(MolecularDynamics* md, int nTemperatures, int* numMDForTemperatures)
{
	int it = rand()%(nTemperatures-1);
	int n = numMDForTemperatures[it];
	int next =  numMDForTemperatures[it+1];
	double delta = (1.0/(Kb*md[n].temperature)- 1.0/(Kb*md[next].temperature))*(md[n].potentialEnergy-md[next].potentialEnergy);
	boolean exchange= FALSE;
	double t;
	double r = 1.0;
	int i,j;
	FILE* f;
	printf("n=%d delta = %f\n",n,delta);
	if(delta<0) exchange= TRUE;
	else
	{
		if(rand()/(double)RAND_MAX>exp(-delta)) exchange = TRUE;
	}
	if(!exchange)return;
	r = sqrt(md[next].temperature/md[n].temperature);
	t = md[n].temperature;
	md[n].temperature = md[next].temperature;
	md[next].temperature = t;
	numMDForTemperatures[it]=next;
	numMDForTemperatures[it+1]=n;
	f = md[n].fileProp;
	md[n].fileProp = md[next].fileProp;
	md[next].fileProp = f;
	f = md[n].fileTraj;
	md[n].fileTraj = md[next].fileTraj;
	md[next].fileTraj = f;

	for ( i = 0; i < md[n].numberOfAtoms; i++)
	{
		for ( j = 0; j < 3; j++)
		{
			 md[n].forceField->molecule.atoms[i].velocity[j] *= r;
			 md[next].forceField->molecule.atoms[i].velocity[j] /= r;
		}
	}
#ifdef ENABLE_CL
	scaleCLVelocities(&md[n], r);
	scaleCLVelocities(&md[next], 1.0/r);
#endif
	removeTranslationAndRotation(&md[n]);
	removeTranslationAndRotation(&md[next]);
}
*/
/**********************************************************************/
static void exchangeReplica(MolecularDynamics* md,  double* energiesAll, double* energies,  double* temperaturesAll, int* numMDForTemperatures, int* nTemperaturesLocal,
	int nTemperaturesAll, int nTemperatures, int nproc, int rank, int iBegin, char* fileNamePrefixProp, char* fileNamePrefixTraj)
{
	if(nproc ==1) 
	{
		//exchangeReplicaLocalAll(md, nTemperatures, numMDForTemperatures);
		int n=-1,next=-1;
		int it;
		int k;
		for(k=0;k<nTemperatures;k++) energiesAll[k] = md[k].potentialEnergy;
		it = getNumExchangeReplica(energiesAll, temperaturesAll, nTemperaturesAll, numMDForTemperatures);
		if(it>=0)
		{
			n = numMDForTemperatures[it];
			next = numMDForTemperatures[it+1];
			numMDForTemperatures[it]=next;
			numMDForTemperatures[it+1]=n;
		}
		if(n>=0 && n<nTemperatures && next>=0 && next<nTemperatures)
		{
			applyExchangeReplica(md, nTemperatures, n,next);
		}
	}
	else
	{
		if(rank==0)
		{
#ifdef ENABLE_MPI
			MPI_Status status ;
			int dum;
			int code;
			int tag=100;
#endif
			int it;
			int n=-1,next=-1;
			int ii;
			int k;
			int j;
			for(k=0;k<nTemperatures;k++) energiesAll[k] = md[k].potentialEnergy;
			for(j=1;j<nproc;j++) 
			{
#if DEBUG
				fprintf(md[0].forceField->logfile, "Get energies from rank %d\n",j);
#endif
#ifdef ENABLE_MPI
				code = MPI_Recv(energies,nTemperaturesLocal[j],MPI_DOUBLE,j,tag,MPI_COMM_WORLD,&status) ;
#endif
#if DEBUG
				fprintf(md[0].forceField->logfile, "End Get energies from rank %d\n",j);
#endif
				for(ii=0;ii<nTemperaturesLocal[j];ii++) energiesAll[k++] = energies[ii];
			}
			it = getNumExchangeReplica(energiesAll, temperaturesAll, nTemperaturesAll, numMDForTemperatures);
#if DEBUG
			fprintf(md[0].forceField->logfile, "\n\nn=%d\n",it);
#endif
			if(it>=0)
			{
				n = numMDForTemperatures[it];
				next = numMDForTemperatures[it+1];
				numMDForTemperatures[it]=next;
				numMDForTemperatures[it+1]=n;
			}
			// send it to all others 
#ifdef ENABLE_MPI
			tag = 200;
#endif
#if DEBUG
			fprintf(md[0].forceField->logfile, "Send n from 0 to all others rank\n");
#endif
#ifdef ENABLE_MPI
			for(j=1;j<nproc;j++) code = MPI_Send(&it,1,MPI_INT,j,tag,MPI_COMM_WORLD) ;
#endif
#if DEBUG
			fprintf(md[0].forceField->logfile, "End Send n from 0 to all others rank\n");
#endif
			// close file for n and next if nessesary
			if(n>=0 && n<nTemperatures && (next<0 || next>=nTemperatures))
			{
				 if(md[n].fileProp) fclose(md[n].fileProp);
				 if(md[n].fileTraj) fclose(md[n].fileTraj);
			}
			if(next>=0 && next<nTemperatures && (n<0||n>=nTemperatures) )
			{
				 if(md[next].fileProp) fclose(md[next].fileProp);
				 if(md[next].fileTraj) fclose(md[next].fileTraj);
			}
#ifdef ENABLE_MPI
			tag = 300;
			for(j=1;j<nproc;j++) 
			{
#if DEBUG
				fprintf(md[0].forceField->logfile, "Get dum val for check close file at rank %d\n",j);
#endif
				code = MPI_Recv(&dum,1,MPI_INT,j,tag,MPI_COMM_WORLD,&status) ;
#if DEBUG
				fprintf(md[0].forceField->logfile, "End Get dum val for check close file at rank %d\n",j);
#endif
			}
#endif
			/// all files closed on other rank
#if DEBUG
			fprintf(md[0].forceField->logfile, "Send n from 0 to all others rank\n");
#endif
#ifdef ENABLE_MPI
			tag = 400;
			for(j=1;j<nproc;j++) code = MPI_Send(&it,1,MPI_INT,j,tag,MPI_COMM_WORLD) ;
#endif

			if(n>=0 && n<nTemperatures && next>=0 && next<nTemperatures)
			{
#if DEBUG
				fprintf(md[0].forceField->logfile, " Exchange between 2 trajs on one proc rank = 0\n");
#endif
				// 2 trajs on one proc
				applyExchangeReplica(md, nTemperatures, n,next);
			}
			else 
			{
				if(n>=0 && n<nTemperatures) changeOneReplica(md,  n, temperaturesAll[it+1], fileNamePrefixProp, fileNamePrefixTraj);
				if(next>=0 && next<nTemperatures) changeOneReplica(md, next, temperaturesAll[it], fileNamePrefixProp, fileNamePrefixTraj);
			}
		}
#ifdef ENABLE_MPI
		else
		{
			int code, tag=100;
			MPI_Status status ;
			int it;
			int n=-1,next=-1;
			int k;
			for(k=0;k<nTemperatures;k++) energies[k] = md[k].potentialEnergy;
#if DEBUG
			fprintf(md[0].forceField->logfile, "Send energies from %d to 0\n",rank);
#endif
			code = MPI_Send(energies,nTemperatures,MPI_DOUBLE,0,tag,MPI_COMM_WORLD) ;
#if DEBUG
			fprintf(md[0].forceField->logfile, "End Send energies from %d to 0\n",rank);
#endif
			tag = 200;
			code = MPI_Recv(&it,1,MPI_INT,0,tag,MPI_COMM_WORLD,&status) ;
			if(it>=0)
			{
				n = numMDForTemperatures[it];
				next = numMDForTemperatures[it+1];
				numMDForTemperatures[it]=next;
				numMDForTemperatures[it+1]=n;
			}
			n -= iBegin;
			next -= iBegin;
			// close file for n and next if nessesar
			if(n>=0 && n<nTemperatures && (next<0 || next>=nTemperatures))
			{
				 if(md[n].fileProp) fclose(md[n].fileProp);
				 if(md[n].fileTraj) fclose(md[n].fileTraj);
			}
			if(next>=0 && next<nTemperatures && (n<0||n>=nTemperatures) )
			{
				 if(md[next].fileProp) fclose(md[next].fileProp);
				 if(md[next].fileTraj) fclose(md[next].fileTraj);
			}
			// send to 0 for tell it all files closed
			tag = 300;
			code = MPI_Send(&it,1,MPI_INT,0,tag,MPI_COMM_WORLD) ;
			tag = 400;
			code = MPI_Recv(&it,1,MPI_INT,0,tag,MPI_COMM_WORLD,&status) ;
			if(n>=0 && n<nTemperatures && next>=0 && next<nTemperatures)
			{
#if DEBUG
				fprintf(md[0].forceField->logfile, " Exchange between 2 trajs on one proc rank = 0\n");
#endif
				// 2 trajs on one proc
				applyExchangeReplica(md, nTemperatures, n,next);
			}
			else 
			{
				if(n>=0 && n<nTemperatures) changeOneReplica(md,  n, temperaturesAll[it+1], fileNamePrefixProp, fileNamePrefixTraj);
				if(next>=0 && next<nTemperatures) changeOneReplica(md, next, temperaturesAll[it], fileNamePrefixProp, fileNamePrefixTraj);
			}
		}
#endif
	}
}
/**********************************************************************/
ForceField**   runREMD(
		MolecularDynamics* molecularDynamics, ForceField* forceField, 
		int updateFrequency, 
		double heatTime, double equiTime, double runTime, 
		double heatTemperature, double runTemperature, double runTemperatureMax,
		double stepSize,
		MDIntegratorType integratorType,
		MDThermostatType thermostat,
		double friction,
		double omegaMax,
		int Nf,
		double collide,
		double qNH,
		int numberOfGeometries,
		int nTemperaturesAll,
		int numberOfExchanges,
		char* fileNameTraj,
		char* fileNameProp
		)
{
	char* fileNamePrefixProp = NULL;
	char* fileNamePrefixTraj = NULL;
	int i;
	int j;
	int k;
	char* str = NULL;
        double gradientNorm = 0;
	int numberOfHeatSteps = 0;
	int numberOfEquiSteps = 0;
	int numberOfRunSteps = 0;
	double currentTemp;
	int* updateNumber = NULL;
	int n0 = 0;
	ForceField** geometries = NULL;
	int iSel = 0;
	int stepSel = 1;
	int stepExchange = 1;
	MolecularDynamics* md = NULL;
	double* runTemps = NULL;
	double a = 1.0;
	double b = 1.0;
	int nproc = 1;
	int* nTemperaturesLocal = NULL;
	int rank = 0;
	int n = 0;
	int nTemperatures = 0;
	int iBegin = 0;
	double* energies = NULL;
	double* energiesAll = NULL;
	double* temperaturesAll = NULL;
	int* numMDForTemperatures = NULL;
	FILE* logfile = forceField->logfile;

	if(nTemperaturesAll<1) return NULL;
	if(forceField->molecule.nAtoms<1) return NULL;
	if(numberOfGeometries<2) return NULL;

#ifdef ENABLE_MPI
	MPI_Comm_rank( MPI_COMM_WORLD,&rank);
	MPI_Comm_size( MPI_COMM_WORLD,&nproc );
	fprintf(logfile, "Rank#=%d  nproc = %d\n", rank, nproc );
	if(nTemperaturesAll<nproc)
	{
		fprintf(logfile, "nTemperatures<number of procs. I set #proc to nTemperatures=%d\n", nTemperaturesAll);
		nproc = nTemperaturesAll;
	}
	if(rank>nproc-1) return NULL;
#else
	rank = 0;
	nproc = 1;
#endif
	nTemperaturesLocal = malloc(nproc*sizeof(int));
	for(i=0;i<nproc;i++) nTemperaturesLocal[i] = 0;
	n = 0;
	i = nproc-1;
	do{
		nTemperaturesLocal[i]++;
		i--;
		if(i<0) i = nproc-1;
		n++;
		fprintf(logfile,"i=%d n = %d\n",i,n);
	}while(n<nTemperaturesAll);
	/* fprintf(logfile, "basname = %s\n",g_path_get_basename(fileNameTraj));*/
	nTemperatures = nTemperaturesLocal[rank];
	energies = malloc(nTemperaturesAll*sizeof(double));
	energiesAll = malloc(nTemperaturesAll*sizeof(double));
	temperaturesAll = malloc(nTemperaturesAll*sizeof(double));

	iBegin = 0;
	for(i=0;i<rank;i++) iBegin += nTemperaturesLocal[i];
	fprintf(logfile, "Rank#=%d  nTemperatures at this proc  = %d nTemperaturesAll = %d iBegin = %d\n", rank, nTemperatures, nTemperaturesAll,iBegin);

	runTemps = malloc(nTemperatures*sizeof(double));
	if(nTemperaturesAll>1) a = pow(runTemperatureMax/runTemperature,1.0/(nTemperaturesAll-1));
	b = 1.0;
	for(i=0;i<nTemperaturesAll;i++) 
	{
		if(i>=iBegin && i<iBegin+nTemperatures) runTemps[i-iBegin] = runTemperature*b;
		temperaturesAll[i] = runTemperature*b;
		b = b*a;
	}
	fprintf(logfile, "\nTemperatures = ");
	for(k=0;k<nTemperatures;k++) fprintf(logfile, "%f ", runTemps[k]);
	fprintf(logfile, "\n\n");
	fflush(logfile);
	numMDForTemperatures = malloc(nTemperaturesAll*sizeof(int));
	for(i=0;i<nTemperaturesAll;i++)  numMDForTemperatures[i] = i;

	md = malloc(nTemperatures*sizeof(MolecularDynamics));
	for(k=0;k<nTemperatures;k++) 
	{
		md[k] = *molecularDynamics;
		md[k].forceField = malloc(sizeof(ForceField)); 
		*md[k].forceField = copyForceField(forceField);

#ifdef ENABLE_CL
		if(k!=0) initCLForceField (md[k].forceField);
#endif
		md[k].numberOfAtoms = forceField->molecule.nAtoms;
		fprintf(logfile, "nAtoms = %d\n",md[k].numberOfAtoms);
		md[k].updateFrequency = updateFrequency;
	}
	updateNumber = malloc(nTemperatures*sizeof(int));
	for(k=0;k<nTemperatures;k++) updateNumber[k] = 0;

	//if(rank==0)
	{
		geometries = malloc(numberOfGeometries*sizeof(ForceField*));
		for(i=0;i<numberOfGeometries;i++) geometries[i] = NULL;
	}

	currentTemp = heatTemperature/2;
	numberOfHeatSteps = heatTime/stepSize*1000;
	numberOfEquiSteps = equiTime/stepSize*1000;; 
	numberOfRunSteps = runTime/stepSize*1000;; 

	for(k=0;k<nTemperatures;k++) 
	{
		char* fileNamePropk = NULL;
		char* fileNameTrajk = NULL;
		{
			if(fileNameProp) fileNamePrefixProp = getSuffixNameFile(fileNameProp);
			if(fileNamePrefixProp) fileNamePropk = strdup_printf("%s%0.0f.txt",fileNamePrefixProp,runTemps[k]);
			if(fileNameTraj) fileNamePrefixTraj = getSuffixNameFile(fileNameTraj);
			if(fileNamePrefixTraj) fileNameTrajk = strdup_printf("%s%0.0f.gab",fileNamePrefixTraj,runTemps[k]);
		}
		currentTemp = heatTemperature;
		if(numberOfHeatSteps==0) currentTemp = runTemps[k];
		/* printf("Begin initMD\n");*/
		initMD(&md[k],currentTemp,stepSize, integratorType, thermostat, friction, omegaMax, Nf, collide, qNH, fileNameTrajk, fileNamePropk, numberOfRunSteps,k+iBegin);
		/* printf("End initMD\n");*/
	}
	for(k=0;k<nTemperatures;k++) 
	{
		md[k].forceField->klass->calculateGradient(md[k].forceField);
		computeEnergies(&md[k]);
	}
	/* printf("End computeEnergies\n");*/

	iSel = -1;
	if(rank==0)
	{
		if(str) free(str);
		str = strdup_printf("Geometry selected Potential energy =  %0.4f", md[0].potentialEnergy);
		/* redrawMolecule(&md[0].forceField->molecule,str);*/
		fprintf(logfile, "%s\n",str);
		iSel++;
		geometries[iSel] = malloc(sizeof(ForceField));
		*geometries[iSel] = copyForceField(md[0].forceField);
		/* waiting(0.1);*/
	}


	for(k=0;k<nTemperatures;k++) 
	{
		md[k].temperature = heatTemperature;
		rescaleVelocities(&md[k]);
		newProperties(&md[k]," ");
	}

	currentTemp = heatTemperature;
	n0 = 0;
	for (i = 0; i < numberOfHeatSteps; i++ )
	for(k=0;k<nTemperatures;k++) 
	{
		applyOneStep(&md[k],i);
		currentTemp = heatTemperature + ( runTemps[k] - heatTemperature ) * ( ( double )( i + 1 )/ numberOfHeatSteps );
		md[k].temperature = currentTemp;
		rescaleVelocities(&md[k]);
		if (++updateNumber[k] >= md[k].updateFrequency )
		{
			if(str) free(str);
			str = strdup_printf(("MD Heating: %0.2f fs, T = %0.2f K T(t) = %8.2f Kin = %0.4f Pot =  %0.4f Tot =  %0.4f"), 
					i*stepSize,
					md[k].temperature, 
					md[k].kelvin, 
					md[k].kineticEnergy,
					md[k].potentialEnergy,
					md[k].totalEnergy
					);
			/* redrawMolecule(&md[k].forceField->molecule,str);*/
			fprintf(logfile, "%s\n",str);
			updateNumber[k] = 0;
		}
		//saveProperties(&md[k], n0+i+1, i+1," Heating");
	}

	for(k=0;k<nTemperatures;k++) 
	{
		md[k].temperature = runTemps[k];
		rescaleVelocities(&md[k]);
		updateNumber[k] = md[k].updateFrequency;
	}

	n0 += numberOfHeatSteps;
	for (i = 0; i < numberOfEquiSteps; i++ )
	for(k=0;k<nTemperatures;k++) 
	{
		md[k].temperature = runTemps[k];
		applyOneStep(&md[k],i);
		md[k].temperature = runTemps[k];
		/* rescaleVelocities(&md[k]);*/
		applyThermostat(&md[k]);
		if (++updateNumber[k] >= md[k].updateFrequency )
		{
			currentTemp =  runTemps[k];
			if(str) free(str);
			str = strdup_printf(("MD Equilibrium: %0.2f fs, T = %0.2f K  T(t) = %8.2f K Kin = %0.4f Pot =  %0.4f Tot =  %0.4f"), 
					i*stepSize,
					md[k].temperature, 
					md[k].kelvin, 
					md[k].kineticEnergy,
					md[k].potentialEnergy,
					md[k].totalEnergy
					);
			//redrawMolecule(&md[k].forceField->molecule,str);
			fprintf(logfile, "%s\n",str);
			updateNumber[k] = 0;
		}
		//saveProperties(&md[k], n0+i+1, i+1, " Equilibrium");
	}
	for(k=0;k<nTemperatures;k++) 
	{
		md[k].temperature = runTemps[k];
		rescaleVelocities(&md[k]);
		updateNumber[k] = md[k].updateFrequency;
	}

	n0 += numberOfEquiSteps;
/*
	if(str) free(str);
	str = strdup_printf(("Geometry selected Potential energy =  %0.4f"), md[0].potentialEnergy);
	fprintf(logfile, "%s\n",str);
*/
	if(numberOfGeometries>2) stepSel = numberOfRunSteps/(numberOfGeometries-1);
	else stepSel = numberOfRunSteps;
	if(numberOfExchanges>2) stepExchange =  numberOfRunSteps/(numberOfExchanges-1);
	else stepExchange = numberOfRunSteps;

	for(k=0;k<nTemperatures;k++) 
			md[k].temperature = runTemps[k];
	for (i = 0; i < numberOfRunSteps; i++ )
	{
		for(k=0;k<nTemperatures;k++) 
		{
			applyOneStep(&md[k],i);
			applyThermostat(&md[k]);
			if (++updateNumber[k] >= md[k].updateFrequency )
			{
				if(str) free(str);
				str = strdup_printf(("MD Running: %0.2f fs, T = %0.2f K  T(t) = %8.2f K Kin = %0.4f Pot =  %0.4f Tot =  %0.4f"), 
					i*stepSize, 
					md[k].temperature, 
					md[k].kelvin, 
					md[k].kineticEnergy,
					md[k].potentialEnergy,
					md[k].totalEnergy
					);
				//redrawMolecule(&md[k].forceField->molecule,str);
				fprintf(logfile, "%s\n",str);
				updateNumber[k] = 0;
				saveTrajectory(&md[k], i+1);
			}
			saveProperties(&md[k], n0+i+1, i+1," Running");
		}
		for(k=0;k<nTemperatures;k++) 
		if(fabs(temperaturesAll[0]-md[k].temperature)<1e-10 && (i+1)%stepSel==0 && (iSel+1)<numberOfGeometries)
		{
			if(str) free(str);
			str = strdup_printf(("Geometry selected Potential energy =  %0.4f"), md[k].potentialEnergy);
			//redrawMolecule(&md[k].forceField->molecule,str);
			fprintf(logfile, "%s\n",str);
			iSel++;
			geometries[iSel] = malloc(sizeof(ForceField));
			*geometries[iSel] = copyForceField(md[k].forceField);
			/* waiting(0.1);*/
		}
/* Exchange here */
		if((i+1)%stepExchange==0&&nTemperaturesAll>1)
		{
			exchangeReplica(md,  energiesAll, energies,  temperaturesAll, numMDForTemperatures, nTemperaturesLocal,
			nTemperaturesAll, nTemperatures, nproc, rank, iBegin, fileNamePrefixProp, fileNamePrefixTraj);
		}
	}
	for(k=0;k<nTemperatures;k++) 
	if(fabs(temperaturesAll[0]-md[k].temperature)<1e-10 && iSel<numberOfGeometries-1)
	{
		if(str) free(str);
		str = strdup_printf(("Geometry selected Potential energy =  %0.4f"), md[k].potentialEnergy);
		//redrawMolecule(&md[0].forceField->molecule,str);
		fprintf(logfile, "%s\n",str);
		iSel++;
		geometries[iSel] = malloc(sizeof(ForceField));
		*geometries[iSel] = copyForceField(md[k].forceField);
		/* waiting(0.1);*/
	}

	n0 += numberOfRunSteps;

	fprintf(logfile,"End of MD Simulation on rank # %d.\n",rank);
	for(k=0;k<nTemperatures;k++) 
	{
		md[k].forceField->klass->calculateGradient(md[k].forceField);
        	gradientNorm = 0;
		for (i = 0; i < md[k].numberOfAtoms; i++)
			for ( j = 0; j < 3; j++)
                        	gradientNorm += 
					md[k].forceField->molecule.atoms[i].gradient[j] * 
					md[k].forceField->molecule.atoms[i].gradient[j]; 

        	gradientNorm = sqrt( gradientNorm );
		if(str) free(str);
		str = strdup_printf(("T(K)=%0.2f Gradient = %f Ekin = %f (Kcal/mol) EPot =  %0.4f ETot =  %0.4f T(t) = %0.2f"),
			runTemps[k],
			(double)gradientNorm,
			md[k].kineticEnergy,
			md[k].potentialEnergy,
			md[k].totalEnergy,
			md[k].kelvin 
			); 
		//redrawMolecule(&md[0].forceField->molecule,str);
		fprintf(logfile, "%s\n",str);
	}
	if(str) free(str);

	for(k=0;k<nTemperatures;k++) 
	{
		if(md[k].fileTraj)fclose(md[k].fileTraj);
		if(md[k].fileProp)fclose(md[k].fileProp);
	}
	
#if DEBUG
	fprintf(logfile, "Begin freeMolecularDynamics in REMD\n");
	fflush(logfile);
#endif
	for(k=0;k<nTemperatures;k++) 
	{
		freeMolecularDynamics(&md[k]);
	}
#if DEBUG
	fprintf(logfile, "End freeMolecularDynamics in REMD\n");
	fflush(logfile);
#endif

	if(updateNumber) free(updateNumber);
#if DEBUG
	fprintf(logfile, "End free updateNumber\n");
	fflush(logfile);
#endif

	if(runTemps) free(runTemps);

#if DEBUG
	fprintf(logfile, "End free runTemps\n");
	fflush(logfile);
#endif

	fflush(logfile);
	return geometries;
}
/**********************************************************************/
ForceField**    runMolecularDynamicsConfo(
		MolecularDynamics* molecularDynamics, ForceField* forceField, 
		int updateFrequency, 
		double heatTime, double equiTime, double runTime,
		double heatTemperature, double equiTemperature, double runTemperature,
		double stepSize,
		MDIntegratorType integratorType,
		MDThermostatType thermostat,
		double friction,
		double omegaMax,
		int Nf,
		double collide,
		double qNH,
		int numberOfGeometries,
		char* fileNameTraj,
		char* fileNameProp
		)
{
	int i;
	int j;
	char* str = NULL;
        double gradientNorm = 0;
	int numberOfHeatSteps = 0;
	int numberOfEquiSteps = 0;
	int numberOfRunSteps = 0;
	double currentTemp;
	int updateNumber = 0;
	int n0 = 0;
	ForceField** geometries = NULL;
	int iSel = 0;
	int stepSel = 1;
	double e0 = 0;
	double esum = 0;
	double e2sum = 0;
	/* 
	 *  physical constants in SI units
	 *   ------------------------------
	 *      Kb = 1.380662 E-23 J/K
	 *      Na = 6.022045 E23  1/mol
	 *      e = 1.6021892 E-19 C
	 *      eps = 8.85418782 E-12 F/m
	 *                       
	 *      1 Kcal = 4184.0 J
	 *      1 amu = 1.6605655 E-27 Kg
	 *      1 A = 1.0 E-10 m
	 *                                       
	 *       Internally, AKMA units are used:
	 *                                        
	 *       timeFactor = SQRT ( ( 1A )**2 * 1amu * Na  / 1Kcal )
	 *       kBoltzmann = Na *Kb  / 1 Kcal
	*/ 

	/* printf("basname = %s\n",g_path_get_basename(fileNameTraj));*/

	if(forceField->molecule.nAtoms<1) return NULL;
	if(numberOfGeometries<2) return NULL;
	geometries = malloc(numberOfGeometries*sizeof(ForceField*));
	for(i=0;i<numberOfGeometries;i++) geometries[i] = NULL;

	molecularDynamics->forceField = forceField;
	molecularDynamics->numberOfAtoms = forceField->molecule.nAtoms;
	molecularDynamics->updateFrequency = updateFrequency;

	currentTemp = heatTemperature/2;
	
	numberOfHeatSteps = heatTime/stepSize*1000;
	numberOfEquiSteps = equiTime/stepSize*1000;; 
	numberOfRunSteps = runTime/stepSize*1000;; 


	currentTemp = heatTemperature;
	if(numberOfHeatSteps==0) currentTemp = equiTemperature; 
	if(numberOfHeatSteps==0 && numberOfEquiSteps==0 ) currentTemp = runTemperature; 

	initMD(molecularDynamics,currentTemp,stepSize, integratorType, thermostat, friction, omegaMax, Nf, collide, qNH, fileNameTraj, fileNameProp, numberOfRunSteps,0);
	molecularDynamics->forceField->klass->calculateGradient(molecularDynamics->forceField);
	computeEnergies(molecularDynamics);
	e0 = molecularDynamics->potentialEnergy;
	printf("E0 = The first potential energy in kcal = %f\n",e0); 

	iSel = -1;
	{
		if(str) free(str);
		str = strdup_printf(("Geometry #%d selected Potential energy =  %0.4f"), iSel+1,molecularDynamics->potentialEnergy);
		/* redrawMolecule(&molecularDynamics->forceField->molecule,str);*/
		printf("%s\n",str);
		iSel++;
		geometries[iSel] = malloc(sizeof(ForceField));
		*geometries[iSel] = copyForceField(molecularDynamics->forceField);
		/* waiting(0.1);*/
	}

	molecularDynamics->temperature = heatTemperature;
	rescaleVelocities(molecularDynamics);

	currentTemp = heatTemperature;
	n0 = 0;
	newProperties(molecularDynamics," ");
	/*newProperties(molecularDynamics," ----> Heating");*/
	for (i = 0; i < numberOfHeatSteps; i++ )
	{
		molecularDynamics->temperature = currentTemp;
		applyOneStep(molecularDynamics,i);
		currentTemp = heatTemperature + ( runTemperature - heatTemperature ) *
				( ( double )( i + 1 )/ numberOfHeatSteps );
		molecularDynamics->temperature = currentTemp;
		rescaleVelocities(molecularDynamics);
		if (++updateNumber >= molecularDynamics->updateFrequency )
		{
			if(str) free(str);
			str = strdup_printf(("MD Heating: %0.2f fs, T = %0.2f K T(t) = %8.2f Kin = %0.4f Pot =  %0.4f Tot =  %0.4f"), 
					i*stepSize, 
					molecularDynamics->temperature, 
					molecularDynamics->kelvin, 
					molecularDynamics->kineticEnergy,
					molecularDynamics->potentialEnergy,
					molecularDynamics->totalEnergy
					);
			/* redrawMolecule(&molecularDynamics->forceField->molecule,str);*/
			printf("%s\n",str);
			updateNumber = 0;
		}
		//saveProperties(molecularDynamics, n0+i+1, i+1," Heating");
	}

	currentTemp = equiTemperature;
	molecularDynamics->temperature = currentTemp;
	rescaleVelocities(molecularDynamics);
	updateNumber = molecularDynamics->updateFrequency;
	n0 += numberOfHeatSteps;
	/* newProperties(molecularDynamics," ----> Equilibrium");*/
	for (i = 0; i < numberOfEquiSteps; i++ )
	{
		molecularDynamics->temperature = currentTemp;
		applyOneStep(molecularDynamics,i);
		molecularDynamics->temperature = currentTemp;
		/*rescaleVelocities(molecularDynamics);*/
		applyThermostat(molecularDynamics);
		if (++updateNumber >= molecularDynamics->updateFrequency )
		{
			if(str) free(str);
			str = strdup_printf(("MD Equilibrium: %0.2f fs, T = %0.2f K  T(t) = %8.2f K Kin = %0.4f Pot =  %0.4f Tot =  %0.4f"), 
					i*stepSize,
					molecularDynamics->temperature, 
					molecularDynamics->kelvin, 
					molecularDynamics->kineticEnergy,
					molecularDynamics->potentialEnergy,
					molecularDynamics->totalEnergy
					);
			//redrawMolecule(&molecularDynamics->forceField->molecule,str);
			printf("%s\n",str);
			updateNumber = 0;
		}
		//saveProperties(molecularDynamics, n0+i+1, i+1, " Equilibrium");
	}
	updateNumber = molecularDynamics->updateFrequency;

	currentTemp = runTemperature;
	molecularDynamics->temperature = currentTemp;
	rescaleVelocities(molecularDynamics);
	updateNumber = molecularDynamics->updateFrequency;
	n0 += numberOfEquiSteps;
	/* newProperties(molecularDynamics," ----> Runing");*/
	//if(str) free(str);
	//str = strdup_printf(("Geometry selected Potential energy =  %0.4f"), molecularDynamics->potentialEnergy);
	//redrawMolecule(&molecularDynamics->forceField->molecule,str);
	//printf("%s\n",str);
	if(numberOfGeometries>2) stepSel = numberOfRunSteps/(numberOfGeometries-1);
	else stepSel = numberOfRunSteps;
	/* printf("Isel = %d\n",stepSel);*/
	esum = 0;
	e2sum = 0;
	for (i = 0; i < numberOfRunSteps; i++ )
	{
		molecularDynamics->temperature = currentTemp;
		applyOneStep(molecularDynamics,i);
		applyThermostat(molecularDynamics);
		esum  += molecularDynamics->totalEnergy;
		e2sum += molecularDynamics->totalEnergy*molecularDynamics->totalEnergy;
		if (++updateNumber >= molecularDynamics->updateFrequency )
		{
			if(str) free(str);
			str = strdup_printf(("MD Running: %0.2f fs, T = %0.2f K  T(t) = %8.2f K Kin = %0.4f Pot =  %0.4f Tot =  %0.4f Eav = %0.4f sigE = %0.4f Eav-E0(cm^-1) = %0.2f"), 
					i*stepSize, 
					molecularDynamics->temperature, 
					molecularDynamics->kelvin, 
					molecularDynamics->kineticEnergy,
					molecularDynamics->potentialEnergy,
					molecularDynamics->totalEnergy,
					esum/(i+1),
					sqrt(fabs(e2sum/(i+1)-esum/(i+1)*esum/(i+1))),
					(esum/(i+1)-e0)*349.75511054
					);
			//redrawMolecule(&molecularDynamics->forceField->molecule,str);
			printf("%s\n",str);
			updateNumber = 0;
			saveTrajectory(molecularDynamics, i+1);
		}
		if((i+1)%stepSel==0 && (iSel+1)<numberOfGeometries)
		{
			if(str) free(str);
			str = strdup_printf(("Geometry #%d selected Potential energy =  %0.4f"), iSel+1,molecularDynamics->potentialEnergy);
			//redrawMolecule(&molecularDynamics->forceField->molecule,str);
			printf("%s\n",str);
			fflush(stdout);
			iSel++;
			geometries[iSel] = malloc(sizeof(ForceField));
			printf("Begin copy\n");
			fflush(stdout);
			*geometries[iSel] = copyForceField(molecularDynamics->forceField);
			printf("End copy\n");
			fflush(stdout);
			/* waiting(0.1);*/
		}
		saveProperties(molecularDynamics, n0+i+1, i+1," Running");
	}
	if(iSel<numberOfGeometries-1)
	{
		if(str) free(str);
		str = strdup_printf(("Geometry selected Potential energy =  %0.4f"), molecularDynamics->potentialEnergy);
		//redrawMolecule(&molecularDynamics->forceField->molecule,str);
		printf("%s\n",str);
		iSel++;
		geometries[iSel] = malloc(sizeof(ForceField));
		*geometries[iSel] = copyForceField(molecularDynamics->forceField);
		/* waiting(0.1);*/
	}

	updateNumber = molecularDynamics->updateFrequency;
	n0 += numberOfRunSteps;

	molecularDynamics->forceField->klass->calculateGradient(molecularDynamics->forceField);
        gradientNorm = 0;
	for (i = 0; i < molecularDynamics->numberOfAtoms; i++)
		for ( j = 0; j < 3; j++)
                        gradientNorm += 
				molecularDynamics->forceField->molecule.atoms[i].gradient[j] * 
				molecularDynamics->forceField->molecule.atoms[i].gradient[j]; 

        gradientNorm = sqrt( gradientNorm );
	if(str) free(str);
	str = strdup_printf(("End of MD Simulation. Gradient = %f Ekin = %f (Kcal/mol) EPot =  %0.4f ETot =  %0.4f T(t) = %0.2f"),
			(double)gradientNorm,
			molecularDynamics->kineticEnergy,
			molecularDynamics->potentialEnergy,
			molecularDynamics->totalEnergy,
			molecularDynamics->kelvin 
			); 
	//redrawMolecule(&molecularDynamics->forceField->molecule,str);
	printf("%s\n",str);
	free(str);
	if(molecularDynamics->fileTraj)fclose(molecularDynamics->fileTraj);
	if(molecularDynamics->fileProp)fclose(molecularDynamics->fileProp);
	freeMolecularDynamics(molecularDynamics);
	return geometries;
}
/**********************************************************************/
static void printGeometryAndVelocities(MolecularDynamics* molecularDynamics, char* title)
{
	fprintf(stdout,"========================================================================================================================\n");
	fprintf(stdout,"#  Geometry and velocities at %s ; T0(K) = %0.2f\n", title,molecularDynamics->kelvin);
	molecularDynamics->forceField->molecule.klass->addGeometry(& molecularDynamics->forceField->molecule,stdout);
	molecularDynamics->forceField->molecule.klass->addVelocities(& molecularDynamics->forceField->molecule,stdout);
	fprintf(stdout,"========================================================================================================================\n");
}
/**********************************************************************/
void	runMolecularDynamics(
		MolecularDynamics* molecularDynamics, ForceField* forceField, 
		int updateFrequency, 
		double heatTime, double equiTime, double runTime, double coolTime, 
		double heatTemperature, double equiTemperature, double runTemperature, double coolTemperature, 
		double stepSize,
		MDIntegratorType integratorType,
		MDThermostatType thermostat,
		double friction,
		double omegaMax,
		int Nf,
		double collide,
		double qNH,
		char* fileNameTraj,
		char* fileNameProp
		)
{
	int i;
	int j;
	char* str = NULL;
        double gradientNorm = 0;
	int numberOfHeatSteps = 0;
	int numberOfEquiSteps = 0;
	int numberOfRunSteps = 0;
	int numberOfCoolSteps = 0;
	double currentTemp;
	int updateNumber = 0;
	int n0 = 0;
	double e0 = 0;
	double esum = 0;
	double e2sum = 0;
	/* 
	 *  physical constants in SI units
	 *   ------------------------------
	 *      Kb = 1.380662 E-23 J/K
	 *      Na = 6.022045 E23  1/mol
	 *      e = 1.6021892 E-19 C
	 *      eps = 8.85418782 E-12 F/m
	 *                       
	 *      1 Kcal = 4184.0 J
	 *      1 amu = 1.6605655 E-27 Kg
	 *      1 A = 1.0 E-10 m
	 *                                       
	 *       Internally, AKMA units are used:
	 *                                        
	 *       timeFactor = SQRT ( ( 1A )**2 * 1amu * Na  / 1Kcal )
	 *       kBoltzmann = Na *Kb  / 1 Kcal
	*/ 

	/* printf("basname = %s\n",g_path_get_basename(fileNameTraj));*/

	if(forceField->molecule.nAtoms<1) return;

	molecularDynamics->forceField = forceField;
	molecularDynamics->numberOfAtoms = forceField->molecule.nAtoms;
	molecularDynamics->updateFrequency = updateFrequency;

	currentTemp = heatTemperature/2;
	
	numberOfHeatSteps = heatTime/stepSize*1000;
	numberOfEquiSteps = equiTime/stepSize*1000;; 
	numberOfRunSteps = runTime/stepSize*1000;; 
	numberOfCoolSteps = coolTime/stepSize*1000;;


	currentTemp = heatTemperature;
	if(numberOfHeatSteps==0) currentTemp = equiTemperature; 
	if(numberOfHeatSteps==0 && numberOfEquiSteps==0 ) currentTemp = runTemperature; 
	if(numberOfHeatSteps==0 && numberOfEquiSteps==0 && numberOfRunSteps==0 ) currentTemp = coolTemperature; 

	initMD(molecularDynamics,currentTemp,stepSize, integratorType, thermostat, friction, omegaMax, Nf, collide, qNH, fileNameTraj, fileNameProp, numberOfRunSteps,0);
	molecularDynamics->forceField->klass->calculateGradient(molecularDynamics->forceField);

	computeEnergies(molecularDynamics);

	e0 = molecularDynamics->potentialEnergy;
	printf("E0 = The first potential energy in kcal = %f\n",e0); 

	molecularDynamics->temperature = heatTemperature;
	if(numberOfHeatSteps>0) rescaleVelocities(molecularDynamics);

	currentTemp = heatTemperature;
	n0 = 0;
	newProperties(molecularDynamics," ");
	/*newProperties(molecularDynamics," ----> Heating");*/
	if(numberOfHeatSteps>0) printGeometryAndVelocities(molecularDynamics, "the begining of Heating stage");
	for (i = 0; i < numberOfHeatSteps; i++ )
	{
		molecularDynamics->temperature = currentTemp;
		applyOneStep(molecularDynamics,i);
		currentTemp = heatTemperature + ( runTemperature - heatTemperature ) *
				( ( double )( i + 1 )/ numberOfHeatSteps );
		molecularDynamics->temperature = currentTemp;
		rescaleVelocities(molecularDynamics);
		if (++updateNumber >= molecularDynamics->updateFrequency )
		{
			if(str) free(str);
			str = strdup_printf(("MD Heating: %0.2f fs, T = %0.2f K T(t) = %0.2f Kin = %0.4f Pot =  %0.4f Tot =  %0.4f"), 
					i*stepSize, 
					molecularDynamics->temperature, 
					molecularDynamics->kelvin, 
					molecularDynamics->kineticEnergy,
					molecularDynamics->potentialEnergy,
					molecularDynamics->totalEnergy
					);
			//redrawMolecule(&molecularDynamics->forceField->molecule,str);
			printf("%s\n",str);
			updateNumber = 0;
		}
		saveProperties(molecularDynamics, n0+i+1, i+1," Heating");
	}

	currentTemp = equiTemperature;
	molecularDynamics->temperature = currentTemp;
	if(numberOfHeatSteps>0) rescaleVelocities(molecularDynamics);
	updateNumber = molecularDynamics->updateFrequency;
	n0 += numberOfHeatSteps;
	/* newProperties(molecularDynamics," ----> Equilibrium");*/
	if(numberOfEquiSteps) 
	{
		printGeometryAndVelocities(molecularDynamics, "the begining of Equilibrium stage");
	}
	for (i = 0; i < numberOfEquiSteps; i++ )
	{
		molecularDynamics->temperature = currentTemp;
		applyOneStep(molecularDynamics,i);
		molecularDynamics->temperature = currentTemp;
		/* rescaleVelocities(molecularDynamics);*/
		applyThermostat(molecularDynamics);
		if (++updateNumber >= molecularDynamics->updateFrequency )
		{
			if(str) free(str);
			str = strdup_printf(("MD Equilibrium: %0.2f fs, T = %0.2f K  T(t) = %0.2f K Kin = %0.4f Pot =  %0.4f Tot =  %0.4f"), 
					i*stepSize, 
					molecularDynamics->temperature, 
					molecularDynamics->kelvin, 
					molecularDynamics->kineticEnergy,
					molecularDynamics->potentialEnergy,
					molecularDynamics->totalEnergy
					);
			//redrawMolecule(&molecularDynamics->forceField->molecule,str);
			printf("%s\n",str);
			updateNumber = 0;
		}
		saveProperties(molecularDynamics, n0+i+1, i+1, " Equilibrium");
	}
	updateNumber = molecularDynamics->updateFrequency;

	currentTemp = runTemperature;
	molecularDynamics->temperature = currentTemp;
	/* rescaleVelocities(molecularDynamics);*/
	updateNumber = molecularDynamics->updateFrequency;
	n0 += numberOfEquiSteps;
	/* newProperties(molecularDynamics," ----> Runing");*/
	if(numberOfRunSteps>0) printGeometryAndVelocities(molecularDynamics, "the begining of Production stage");
	esum = 0;
	e2sum = 0;
	for (i = 0; i < numberOfRunSteps; i++ )
	{
		molecularDynamics->temperature = currentTemp;
		applyOneStep(molecularDynamics,i);
		applyThermostat(molecularDynamics);
		esum  += molecularDynamics->totalEnergy;
		e2sum += molecularDynamics->totalEnergy*molecularDynamics->totalEnergy;
		if (++updateNumber >= molecularDynamics->updateFrequency )
		{
			if(str) free(str);
			str = strdup_printf(("MD Running: %0.2f fs, T = %0.2f K  T(t) = %8.2f K Kin = %0.4f Pot =  %0.4f Tot =  %0.4f Eav =  %0.4f sigE =  %0.4f Eav-E0(cm^-1) = %0.2f"), 
					i*stepSize, 
					molecularDynamics->temperature, 
					molecularDynamics->kelvin, 
					molecularDynamics->kineticEnergy,
					molecularDynamics->potentialEnergy,
					molecularDynamics->totalEnergy,
					esum/(i+1),
					sqrt(fabs(e2sum/(i+1)-esum/(i+1)*esum/(i+1))),
					(esum/(i+1)-e0)*349.75511054
					);
			//redrawMolecule(&molecularDynamics->forceField->molecule,str);
			printf("%s\n",str);
			updateNumber = 0;
			saveTrajectory(molecularDynamics, i+1);
		}
		saveProperties(molecularDynamics, n0+i+1, i+1," Running");
	}
	if(numberOfCoolSteps>0) printGeometryAndVelocities(molecularDynamics, "the beginning of Cooling stage");
	updateNumber = molecularDynamics->updateFrequency;
	n0 += numberOfRunSteps;
	/* newProperties(molecularDynamics," ----> Cooling");*/
	for (i = 0; i < numberOfCoolSteps; i++ )
	{
		currentTemp = runTemperature - ( runTemperature - coolTemperature ) * 
				( ( double )( i + 1 )/ numberOfCoolSteps );
		molecularDynamics->temperature = currentTemp;
		rescaleVelocities(molecularDynamics);
		molecularDynamics->temperature = currentTemp;
		applyOneStep(molecularDynamics,i);
		if (++updateNumber >= molecularDynamics->updateFrequency )
		{
			if(str) free(str);
			str = strdup_printf(("MD Cooling: %0.2f fs, T = %0.2f K T(t) = %0.2f K Kin = %0.4f Pot =  %0.4f Tot =  %0.4f"), 
					i*stepSize, 
					molecularDynamics->temperature, 
					molecularDynamics->kelvin, 
					molecularDynamics->kineticEnergy,
					molecularDynamics->potentialEnergy,
					molecularDynamics->totalEnergy
					);
			//redrawMolecule(&molecularDynamics->forceField->molecule,str);
			printf("%s\n",str);
			updateNumber = 0;
		}
		saveProperties(molecularDynamics, n0+i+1, i+1," Cooling");
	}
	molecularDynamics->forceField->klass->calculateGradient(molecularDynamics->forceField);
        gradientNorm = 0;
	for (i = 0; i < molecularDynamics->numberOfAtoms; i++)
		for ( j = 0; j < 3; j++)
                        gradientNorm += 
				molecularDynamics->forceField->molecule.atoms[i].gradient[j] * 
				molecularDynamics->forceField->molecule.atoms[i].gradient[j]; 

        gradientNorm = sqrt( gradientNorm );
	if(str) free(str);
	str = strdup_printf(("End of MD Simulation. Gradient = %f Ekin = %f (Kcal/mol) EPot =  %0.4f ETot =  %0.4f T(t) = %0.2f"),
			(double)gradientNorm,
			molecularDynamics->kineticEnergy,
			molecularDynamics->potentialEnergy,
			molecularDynamics->totalEnergy,
			molecularDynamics->kelvin 
			); 
	//redrawMolecule(&molecularDynamics->forceField->molecule,str);
	printf("%s\n",str);
	free(str);
	printGeometryAndVelocities(molecularDynamics, "the end of simulation");
	if(molecularDynamics->fileTraj)fclose(molecularDynamics->fileTraj);
	if(molecularDynamics->fileProp)fclose(molecularDynamics->fileProp);
	freeMolecularDynamics(molecularDynamics);
}
/*********************************************************************************/
static void removeTranslation(MolecularDynamics* molecularDynamics)
{
#ifdef ENABLE_CL
	//size_t global = molecularDynamics->numberOfAtoms;
	size_t global = 1;
	CLProp clProp=getCLProp();
	cl_int err = clEnqueueNDRangeKernel(clProp.command_queue, molecularDynamics->removeTranslation, 1, NULL, &global, NULL, 0, NULL, NULL);
	if(err != CL_SUCCESS)
	{
		printErrorCLRun(err);
       		printf("I cannot execute removeTranslation\n");
		exit(1);
	}
	clFinish(clProp.command_queue);
#else
	Molecule* mol = &molecularDynamics->forceField->molecule;
	mol->klass->removeTranslation(mol); 
#endif
}
/*********************************************************************************/
static void removeRotation(MolecularDynamics* molecularDynamics)
{
#ifdef ENABLE_CL
	//size_t global = molecularDynamics->numberOfAtoms;
	size_t global = 1;
	CLProp clProp=getCLProp();
	cl_int err = clEnqueueNDRangeKernel(clProp.command_queue, molecularDynamics->removeRotation, 1, NULL, &global, NULL, 0, NULL, NULL);
	if(err != CL_SUCCESS)
	{
		printErrorCLRun(err);
       		printf("I cannot execute removeRotation\n");
		exit(1);
	}
	clFinish(clProp.command_queue);
#else
	Molecule* mol = &molecularDynamics->forceField->molecule;
	mol->klass->removeRotation(mol); 
#endif

}
/*********************************************************************************/
static void removeTranslationAndRotation(MolecularDynamics* molecularDynamics)
{
	removeTranslation(molecularDynamics);
	removeRotation(molecularDynamics);
}
/*********************************************************************************/
static void initNH(MolecularDynamics* molecularDynamics, double qNH)
{
	int i;

	if(molecularDynamics->thermostat != NOSEHOOVER) return;
	for(i=0;i<MAXNH;i++)
	{
		molecularDynamics->xNH[i] = 0;
		molecularDynamics->vNH[i] = 0;
		molecularDynamics->qNH[i] = qNH;
		molecularDynamics->gNH[i] = 0;
	}
}
/*********************************************************************************/
static void initSD(MolecularDynamics* molecularDynamics, double friction)
{
	int i;
	if(friction<0) friction = 40;
	molecularDynamics->friction = friction/(fsInAKMA)/1000;

	molecularDynamics->positionFriction = NULL;
	molecularDynamics->velocityFriction = NULL;
	molecularDynamics->accelarationFriction = NULL;
	molecularDynamics->gamma = NULL;
	molecularDynamics->positionRandom = NULL;
	molecularDynamics->velocityRandom = NULL;

	if(molecularDynamics->integratorType != STOCHASTIC) return;

	molecularDynamics->positionFriction = malloc(molecularDynamics->numberOfAtoms *sizeof(double)); 
	molecularDynamics->velocityFriction = malloc(molecularDynamics->numberOfAtoms *sizeof(double)); 
	molecularDynamics->accelarationFriction = malloc(molecularDynamics->numberOfAtoms *sizeof(double)); 
	molecularDynamics->gamma = malloc(molecularDynamics->numberOfAtoms *sizeof(double)); 

	molecularDynamics->positionRandom = malloc(molecularDynamics->numberOfAtoms *sizeof(double*)); 
	for(i=0;i<molecularDynamics->numberOfAtoms;i++)
		molecularDynamics->positionRandom[i] = malloc(3*sizeof(double));

	molecularDynamics->velocityRandom = malloc(molecularDynamics->numberOfAtoms *sizeof(double*)); 
	for(i=0;i<molecularDynamics->numberOfAtoms;i++)
		molecularDynamics->velocityRandom[i] = malloc(3*sizeof(double));

}
/* Martyna&Tuckerman Symplectic reversible integrators: Predictor-corrector methods */
/* Ref : JCP 102, 8071 (1995) */
/*********************************************************************************/
static void initMartynaTuckerman(MolecularDynamics* molecularDynamics)
{
	int i,j;

	//printf("Begin initMartynaTuckerman\n");
	for(j=0;j<3;j++) molecularDynamics->Ftilde[j] = NULL;
	for(j=0;j<3;j++) molecularDynamics->VF[j] = NULL;
	for(j=0;j<3;j++) molecularDynamics->VJ[j] = NULL;

	if(molecularDynamics->integratorType != MARTYNATUCKERMAN) return;

	for(j=0;j<3;j++) molecularDynamics->Ftilde[j] = malloc(molecularDynamics->numberOfAtoms *sizeof(double)); 
	for(j=0;j<3;j++) molecularDynamics->VF[j] = malloc(molecularDynamics->numberOfAtoms *sizeof(double)); 
	for(j=0;j<3;j++) molecularDynamics->VJ[j] = malloc(molecularDynamics->numberOfAtoms *sizeof(double)); 

	for(i=0;i<molecularDynamics->numberOfAtoms;i++)
		for ( j = 0; j < 3; j++)
			molecularDynamics->Ftilde[j][i] =-molecularDynamics->forceField->molecule.atoms[i].gradient[j]/molecularDynamics->forceField->molecule.atoms[i].mass;

	for(i=0;i<molecularDynamics->numberOfAtoms;i++) for ( j = 0; j < 3; j++) molecularDynamics->VF[j][i] = 0.0;
	for(i=0;i<molecularDynamics->numberOfAtoms;i++) for ( j = 0; j < 3; j++) molecularDynamics->VJ[j][i] = 0.0;
	//printf("End initMartynaTuckerman\n");
}
/*********************************************************************************/
static void removeTranslationAndRotationTheta(MolecularDynamics* molecularDynamics)
{
	Molecule* mol = &molecularDynamics->forceField->molecule;
	mol->klass->removeTranslationAndRotationAcceleration(mol,molecularDynamics->theta); 
}
/*********************************************************************************/
static void computeThetaQTB(MolecularDynamics* molecularDynamics)
{
	double sigma=sqrt(2.*molecularDynamics->friction*Kb*molecularDynamics->temperature/molecularDynamics->h);
	int i,j,k;

	/* thetaProg = thetaPaper*sigma/sqrt(m) */
	for(i=0;i<molecularDynamics->numberOfAtoms;i++)
	for(j=0;j<3;j++)
	{
		molecularDynamics->theta[3*i+j] = 0.0;
		for(k=0;k<2*molecularDynamics->Nf;k++)
			molecularDynamics->theta[3*i+j] += molecularDynamics->rnoise[3*i+j][2*molecularDynamics->Nf-1-k]*molecularDynamics->Ht[k];
		molecularDynamics->theta[3*i+j] *= sigma/sqrt(molecularDynamics->forceField->molecule.atoms[i].mass); 
	}
	removeTranslationAndRotationTheta(molecularDynamics);
}
/*********************************************************************************/
static void resetQTB(MolecularDynamics* molecularDynamics)
{
	double* Filter = NULL;
	int k;
	double hbarwOverkT;
	double hbardwOverkT;
	int i,j;
	double T;
	

	if(molecularDynamics->integratorType != QTB) return;
	T = molecularDynamics->temperature;
	if(T<=0) return;

	/* computing of Filter */
	/* Htild/sqrt(kT), sqrt(kT) in sigma */
	Filter = malloc((2*molecularDynamics->Nf)*sizeof(double)); 
	/* h dOmega = pi /Nf */
	hbardwOverkT = 1.0/(molecularDynamics->Nf*molecularDynamics->h*Kb*T);
	for(k=0;k<2*molecularDynamics->Nf;k++)
	{
		int kk= k-molecularDynamics->Nf;
		if(kk==0) Filter[k] = 1.0;
		else
		{
			hbarwOverkT = fabs(kk)*hbardwOverkT;
			Filter[k] = sqrt(hbarwOverkT*(0.5+1./(exp(hbarwOverkT)-1.0)));
			//Filter[k] = 1.0; // to test classic
			Filter[k] *= (kk*M_PI/molecularDynamics->Nf/2)/sin(kk*M_PI/molecularDynamics->Nf/2);
		}
	}
	/* compute Ht */
	for(j=0;j<2*molecularDynamics->Nf;j++)
	{
                molecularDynamics->Ht[j] = 0;
		for(k=0; k<2*molecularDynamics->Nf;k++)
                	molecularDynamics->Ht[j] += Filter[k]*cos(M_PI*(k-molecularDynamics->Nf)*(j-molecularDynamics->Nf)*1.0/molecularDynamics->Nf);
		
               	molecularDynamics->Ht[j] /= 2*molecularDynamics->Nf;
	}
	free(Filter);
	
	for(i=0;i<molecularDynamics->numberOfAtoms;i++)
	for(j=0;j<3;j++)
	for(k=0;k<2*molecularDynamics->Nf;k++)
		 molecularDynamics->rnoise[3*i+j][k] = normal();/* sqrt(h) is in sigma */

	computeThetaQTB(molecularDynamics);

}
/*********************************************************************************/
/* omegaMax in cm-1 */
static void initQTB(MolecularDynamics* molecularDynamics, double omegaMax, double friction, int Nf)
{
/* Refs 
Jean-Louis Barrat , David Rodney
Portable implementation of a quantum thermal bath for molecular dynamics simulations
JOURNAL OF STATISTICAL PHYSICS 670, 144, (2011)
*/
	static double cmM1fsM1 = 2.99792458e-5;
	double Omegafs = omegaMax*cmM1fsM1;/* fs^-1 */ 
	int i,j;
	

	if(Nf<1) Nf = 50;

	molecularDynamics->Ht = NULL;
	molecularDynamics->theta = NULL;
	molecularDynamics->rnoise = NULL;
	molecularDynamics->Nf = 0;
	molecularDynamics->M = 0;

	if(molecularDynamics->integratorType != QTB) return;

	molecularDynamics->Nf = Nf;
	molecularDynamics->h = 1/Omegafs*(fsInAKMA);
	molecularDynamics->M = (int)(molecularDynamics->h/molecularDynamics->dt);
	if(molecularDynamics->M<1) molecularDynamics->M = 1;
	molecularDynamics->h = molecularDynamics->M *molecularDynamics->dt;
	omegaMax = 1.0/molecularDynamics->h*(fsInAKMA)/cmM1fsM1; /* cm-1 */
	if(friction<0) molecularDynamics->friction = (1.0/ molecularDynamics->h)/50;
	else  molecularDynamics->friction = friction/fsInAKMA/1000.0;

	molecularDynamics->Ht = malloc((2*Nf)*sizeof(double)); 
	molecularDynamics->theta = malloc(3*molecularDynamics->numberOfAtoms *sizeof(double)); 
	molecularDynamics->rnoise = malloc(3*molecularDynamics->numberOfAtoms *sizeof(double*)); 
	for(i=0;i<molecularDynamics->numberOfAtoms;i++)
	for(j=0;j<3;j++)
		molecularDynamics->rnoise[3*i+j] = malloc((2*Nf)*sizeof(double)); 
	

	printf("\n");
	printf("*************** QTB Parameters ******************************************************************\n");
	printf("Nf\t\t= %d\n",molecularDynamics->Nf);
	printf("M\t\t= %d\n",molecularDynamics->M);
	printf("dt(fs)\t\t= %f\n",molecularDynamics->dt/fsInAKMA);
	printf("h(fs)\t\t= %f\n",molecularDynamics->h/fsInAKMA);
	printf("gamma(ps^-1)\t= %f\n",molecularDynamics->friction*fsInAKMA*1000);
	printf("omegaMax(cm^-1)\t= %f\n",omegaMax);
	printf("*************************************************************************************************\n");
	printf("\n");

	resetQTB(molecularDynamics);

}
/*********************************************************************************/
static void updateQTB(MolecularDynamics* molecularDynamics)
{
	int i,j,k;

	if(molecularDynamics->temperature<=0) return;
	computeThetaQTB(molecularDynamics);
	/* shift rnoise */
	for(i=0;i<molecularDynamics->numberOfAtoms;i++)
	for(j=0;j<3;j++)
	for(k=0;k<2*molecularDynamics->Nf-1;k++)
		 molecularDynamics->rnoise[3*i+j][k] = molecularDynamics->rnoise[3*i+j][k+1];

	/* add one value to the end */
	for(i=0;i<molecularDynamics->numberOfAtoms;i++)
	for(j=0;j<3;j++)
		 molecularDynamics->rnoise[3*i+j][2*molecularDynamics->Nf-1] = normal(); /* sqrt(h) in sigma */

}
/*********************************************************************************/
static void resetLangevin(MolecularDynamics* molecularDynamics)
{
	if(molecularDynamics->integratorType != LANGEVIN) return;
	updateLangevin(molecularDynamics);
}
/*********************************************************************************/
/* omegaMax in cm-1 */
static void initLangevin(MolecularDynamics* molecularDynamics, double friction)
{
	if(molecularDynamics->integratorType != LANGEVIN) return;

	if(friction<0) friction = 40;
	molecularDynamics->friction = friction/1000.0/fsInAKMA;
	molecularDynamics->theta = malloc(3*molecularDynamics->numberOfAtoms *sizeof(double)); 
	
	printf("\n");
	printf("*************** Langevin Parameters ******************************************************************\n");
	printf("dt(fs)\t\t= %f\n",molecularDynamics->dt/fsInAKMA);
	printf("gamma(ps^-1)\t= %f\n",molecularDynamics->friction*fsInAKMA*1000);
	printf("*************************************************************************************************\n");
	printf("\n");

	resetLangevin(molecularDynamics);

}
/*********************************************************************************/
static void updateLangevin(MolecularDynamics* molecularDynamics)
{
	int i,j;
	double sigma;
	/* update theta */
	if(molecularDynamics->integratorType != LANGEVIN) return;
	sigma=sqrt(2.*molecularDynamics->friction*Kb*molecularDynamics->temperature/molecularDynamics->dt);
	for(i=0;i<molecularDynamics->numberOfAtoms;i++)
	for(j=0;j<3;j++)
	{
		molecularDynamics->theta[3*i+j] = normal();
		molecularDynamics->theta[3*i+j] *= sigma/sqrt(molecularDynamics->forceField->molecule.atoms[i].mass); 
	}
}
/******************************************************************************************************************/
/*********************************************************************************/
static void initMD(MolecularDynamics* molecularDynamics, double temperature, double stepSize, MDIntegratorType integratorType, MDThermostatType thermostat, double friction, double omegaMax, int Nf, double collide, double qNH, char* fileNameTraj, char* fileNameProp, int numberOfRunSteps, int index)
{

#ifdef DEBUG
	printf("Begin initMD\n");
#endif
#ifdef ENABLE_CL
	CLProp clProp=getCLProp();
	cl_int err;
	cl_float dtf;
	cl_float frictionf;
	cl_float temp = molecularDynamics->temperature;
#endif
	int i;
	int j;
	double dt = stepSize * fsInAKMA;

	molecularDynamics->collide = collide;
	molecularDynamics->potentialEnergy = 0;
	molecularDynamics->kineticEnergy = 0;
	molecularDynamics->totalEnergy = 0;
	molecularDynamics->kelvin = 0;
	molecularDynamics->temperature = temperature;
	molecularDynamics->thermostat = NONE;

	molecularDynamics->integratorType = integratorType;
	molecularDynamics->thermostat = thermostat;
	molecularDynamics->fileTraj = NULL;
	molecularDynamics->fileProp = NULL;

	molecularDynamics->index = index;

	molecularDynamics->a = malloc(molecularDynamics->numberOfAtoms *sizeof(double*)); 
	for(i=0;i<molecularDynamics->numberOfAtoms;i++)
		molecularDynamics->a[i] = malloc(3*sizeof(double));

	molecularDynamics->aold = NULL;
	if(molecularDynamics->integratorType==BEEMAN)
	{
		molecularDynamics->aold = malloc(molecularDynamics->numberOfAtoms *sizeof(double*)); 
		for(i=0;i<molecularDynamics->numberOfAtoms;i++)
			molecularDynamics->aold[i] = malloc(3*sizeof(double));
	}
	molecularDynamics->coordinatesOld = NULL;
	molecularDynamics->moved = NULL;
	molecularDynamics->update = NULL;
	if(molecularDynamics->forceField->molecule.constraints!=NOCONSTRAINTS)
	{
		molecularDynamics->coordinatesOld = malloc(molecularDynamics->numberOfAtoms *sizeof(double*)); 
		for(i=0;i<molecularDynamics->numberOfAtoms;i++)
			molecularDynamics->coordinatesOld[i] = malloc(3*sizeof(double));
		molecularDynamics->moved = malloc(molecularDynamics->numberOfAtoms *sizeof(boolean)); 
		molecularDynamics->update = malloc(molecularDynamics->numberOfAtoms *sizeof(boolean)); 

	}

	if(fileNameTraj)
	{
 		molecularDynamics->fileTraj = fopen(fileNameTraj, "w");
		if(molecularDynamics->fileTraj != NULL)
		{
			fprintf(molecularDynamics->fileTraj,"[Gabedit Format]\n");
			fprintf(molecularDynamics->fileTraj,"\n");
			fprintf(molecularDynamics->fileTraj,"[MD]\n");
			if(molecularDynamics->updateFrequency>0) numberOfRunSteps/=molecularDynamics->updateFrequency;
			fprintf(molecularDynamics->fileTraj," %d\n",numberOfRunSteps);
		}
	}
	if(fileNameProp)
	{
 		molecularDynamics->fileProp = fopen(fileNameProp, "w");
	}

	/* srand ( (unsigned)time (NULL));*/
	
	molecularDynamics->dt = dt;
	molecularDynamics->dt_2 = dt/2.0;
	molecularDynamics->dt_4 = dt/4.0;
	molecularDynamics->dt2_2 = dt*dt/2;;
	molecularDynamics->dt_8 = dt/8.0;
	molecularDynamics->dt2_8 = dt*dt/8.0;
	molecularDynamics->dt2 = dt*dt;;
	molecularDynamics->dt3_6 = dt*dt*dt/6;
	molecularDynamics->dt4_48 = dt*dt*dt*dt/48;

#ifdef DEBUG
	printf("Begin initSD\n");
#endif
	initSD(molecularDynamics, friction);
	//printf("End initSD\n");
	initNH(molecularDynamics,qNH);
	//printf("End initNH\n");
	initQTB(molecularDynamics, omegaMax, friction, Nf);
	initLangevin(molecularDynamics, friction);
#ifdef DEBUG
	printf("Begin calcGrad\n");
#endif

#ifndef ENABLE_CL
	molecularDynamics->forceField->klass->calculateGradient(molecularDynamics->forceField);
#endif

#ifdef DEBUG
	//printf("End calcGrad\n");
#endif
	initMartynaTuckerman(molecularDynamics);
	for ( i = 0; i < molecularDynamics->numberOfAtoms; i++)
	{
		double m = molecularDynamics->forceField->molecule.atoms[i].mass;
		for ( j = 0; j < 3; j++) molecularDynamics->a[i][j] = -molecularDynamics->forceField->molecule.atoms[i].gradient[j]/m;
		if(molecularDynamics->aold) 
			for ( j = 0; j < 3; j++) molecularDynamics->aold[i][j]  = molecularDynamics->a[i][j];
		if(molecularDynamics->coordinatesOld)
		       	for ( j = 0; j < 3; j++) molecularDynamics->coordinatesOld[i][j]  = molecularDynamics->forceField->molecule.atoms[i].coordinates[j];
		if(molecularDynamics->moved) molecularDynamics->moved[i] = FALSE; 
		if(molecularDynamics->update) molecularDynamics->update[i] = FALSE; 
	}
	//printf("End initAl\n");
	/* set velocities */
	/* check if velocities are already read from the input file */
	if(molecularDynamics->forceField->molecule.klass->setMaxwellVelocitiesIfNull(&molecularDynamics->forceField->molecule, temperature)) rescaleVelocities(molecularDynamics);
#ifdef DEBUG
	molecularDynamics->forceField->molecule.klass->print(&molecularDynamics->forceField->molecule,stdout);
	printf("Velocity init\n");
	for ( i = 0; i < molecularDynamics->numberOfAtoms; i++)
	{
		for ( j = 0; j < 3; j++) fprintf(stdout, "%f ",molecularDynamics->forceField->molecule.atoms[i].velocity[j]);
		fprintf(stdout,"\n");

	}
#endif
#ifdef ENABLE_CL
#ifdef DEBUG
	printf("Begin CL\n");
#endif
	dtf =  molecularDynamics->dt;
	frictionf = molecularDynamics->friction;
#ifdef DEBUG
	printf("Begin clBuildProgram\n");
#endif
	// create a program from the kernel source code
	printf("Begin molecularDynamics->programMD\n");
	molecularDynamics->programMD = clCreateProgramWithSource(clProp.context,1,(const char **) &mdCLSource, NULL, &err);
	printf("End molecularDynamics->programMD\n");
#ifdef DEBUG
	printf("err = %d\n",err);
#endif
	// compile the program
	if (clBuildProgram(molecularDynamics->programMD, 0, NULL, NULL, NULL, NULL) != CL_SUCCESS)
	{
		char build[2048];
		printf("Error building MD CL program\n");
		clGetProgramBuildInfo(molecularDynamics->programMD, clProp.device_id, CL_PROGRAM_BUILD_LOG, 2048, build, NULL);
		printf("Build Log:\n%s\n",build);
		exit(1);
	}
#ifdef DEBUG
	printf("End clBuildProgram\n");
#endif
	molecularDynamics->aCL = clCreateBuffer(clProp.context, CL_MEM_READ_WRITE, sizeof(cl_float4) * forceField->molecule.nAtoms*2, NULL, NULL);


	// create kernels
	molecularDynamics->copyAccelarations = clCreateKernel(molecularDynamics->programMD, "copyAccelarations", &err);
	clSetKernelArg(molecularDynamics->copyAccelarations, 0, sizeof(cl_mem), &molecularDynamics->aCL);
	clSetKernelArg(molecularDynamics->copyAccelarations, 1, sizeof(cl_int), &forceField->molecule.nAtoms);
#ifdef DEBUG
	printf("End create kernel copyAccelarations err = %d\n",err);
#endif

	molecularDynamics->computeAccelarationsFromGradients = clCreateKernel(molecularDynamics->programMD, "computeAccelarationsFromGradients", &err);
	printf("ok  clCreateKernel(molecularDynamics->programMD  &err)\n");
	clSetKernelArg(molecularDynamics->computeAccelarationsFromGradients, 0, sizeof(cl_mem),   &forceField->gradientBufferCL);
	printf("ok  gradientBufferCL)\n");
	clSetKernelArg(molecularDynamics->computeAccelarationsFromGradients, 1, sizeof(cl_mem),   &molecularDynamics->aCL);
	printf("ok  aCL\n");
	clSetKernelArg(molecularDynamics->computeAccelarationsFromGradients, 2, sizeof(cl_mem),   &forceField->atomsCL);
	printf("ok  atomsCL\n");
	clSetKernelArg(molecularDynamics->computeAccelarationsFromGradients, 3, sizeof(cl_int), &forceField->molecule.nAtoms);
#ifdef DEBUG
	printf("End create kernel computeAccelarationsFromGradients\n");
#endif

	molecularDynamics->initAccelarations = clCreateKernel(molecularDynamics->programMD, "initAccelarations", &err);
	clSetKernelArg(molecularDynamics->initAccelarations, 0, sizeof(cl_mem), &molecularDynamics->aCL);
	clSetKernelArg(molecularDynamics->initAccelarations, 1, sizeof(cl_int), &forceField->molecule.nAtoms);
#ifdef DEBUG
	printf("End create kernel initAccelarations\n");
#endif


	molecularDynamics->computeKineticEnergy = clCreateKernel(molecularDynamics->programMD, "computeKineticEnergy", &err);
	clSetKernelArg(molecularDynamics->computeKineticEnergy, 0, sizeof(cl_mem),   &forceField->energyBufferCL);
	clSetKernelArg(molecularDynamics->computeKineticEnergy, 1, sizeof(cl_mem),   &forceField->atomsCL);
	clSetKernelArg(molecularDynamics->computeKineticEnergy, 2, sizeof(cl_int), &forceField->molecule.nAtoms);
#ifdef DEBUG
	printf("End create kernel computeKineticEnergy\n");
#endif

	molecularDynamics->scaleVelocities = clCreateKernel(molecularDynamics->programMD, "scaleVelocities", &err);
	clSetKernelArg(molecularDynamics->scaleVelocities, 0, sizeof(cl_mem),   &forceField->atomsCL);
	clSetKernelArg(molecularDynamics->scaleVelocities, 1, sizeof(cl_int), &forceField->molecule.nAtoms);

	molecularDynamics->removeTranslation = clCreateKernel(molecularDynamics->programMD, "removeTranslation", &err);
	clSetKernelArg(molecularDynamics->removeTranslation, 0, sizeof(cl_mem), &forceField->atomsCL);
	clSetKernelArg(molecularDynamics->removeTranslation, 1, sizeof(cl_int), &forceField->molecule.nAtoms);
	if(err != CL_SUCCESS) 
	{
		printf("I cannot create removeTranslation\n");
		exit(1);
	}

	molecularDynamics->removeRotation = clCreateKernel(molecularDynamics->programMD, "removeRotation", &err);
	clSetKernelArg(molecularDynamics->removeRotation, 0, sizeof(cl_mem), &forceField->atomsCL);
	clSetKernelArg(molecularDynamics->removeRotation, 1, sizeof(cl_int), &forceField->molecule.nAtoms);
	if(err != CL_SUCCESS) 
	{
		printf("I cannot create removeRotation\n");
		exit(1);
	}

	if(molecularDynamics->integratorType==VERLET) 
	{
		molecularDynamics->applyVerlet1 = clCreateKernel(molecularDynamics->programMD, "applyVerlet1", &err);
		clSetKernelArg(molecularDynamics->applyVerlet1, 0, sizeof(cl_mem),   &molecularDynamics->aCL);
		clSetKernelArg(molecularDynamics->applyVerlet1, 1, sizeof(cl_mem),   &forceField->atomsCL);
		clSetKernelArg(molecularDynamics->applyVerlet1, 2, sizeof(cl_float), &dtf);
		clSetKernelArg(molecularDynamics->applyVerlet1, 3, sizeof(cl_int), &forceField->molecule.nAtoms);

		molecularDynamics->applyVerlet2 = clCreateKernel(molecularDynamics->programMD, "applyVerlet2", &err);
		clSetKernelArg(molecularDynamics->applyVerlet2, 0, sizeof(cl_mem),   &molecularDynamics->aCL);
		clSetKernelArg(molecularDynamics->applyVerlet2, 1, sizeof(cl_mem),   &forceField->atomsCL);
		clSetKernelArg(molecularDynamics->applyVerlet2, 2, sizeof(cl_float), &dtf);
		clSetKernelArg(molecularDynamics->applyVerlet2, 3, sizeof(cl_int), &forceField->molecule.nAtoms);
#ifdef DEBUG
		printf("End create kernel applyVerlet1&2\n");
#endif
	}
	if(molecularDynamics->integratorType==BEEMAN) 
	{
#ifdef DEBUG
		printf("Begin create kernel applyBeeman1&2\n");
#endif
		molecularDynamics->applyBeeman1 = clCreateKernel(molecularDynamics->programMD, "applyBeeman1", &err);
		clSetKernelArg(molecularDynamics->applyBeeman1, 0, sizeof(cl_mem),   &molecularDynamics->aCL);
		clSetKernelArg(molecularDynamics->applyBeeman1, 1, sizeof(cl_mem),   &forceField->atomsCL);
		clSetKernelArg(molecularDynamics->applyBeeman1, 2, sizeof(cl_float), &dtf);
		clSetKernelArg(molecularDynamics->applyBeeman1, 3, sizeof(cl_int), &forceField->molecule.nAtoms);
#ifdef DEBUG
		printf("End create kernel applyBeeman1\n");
#endif

		molecularDynamics->applyBeeman2 = clCreateKernel(molecularDynamics->programMD, "applyBeeman2", &err);
		clSetKernelArg(molecularDynamics->applyBeeman2, 0, sizeof(cl_mem),   &molecularDynamics->aCL);
		clSetKernelArg(molecularDynamics->applyBeeman2, 1, sizeof(cl_mem),   &forceField->atomsCL);
		clSetKernelArg(molecularDynamics->applyBeeman2, 2, sizeof(cl_float), &dtf);
		clSetKernelArg(molecularDynamics->applyBeeman2, 3, sizeof(cl_int), &forceField->molecule.nAtoms);
#ifdef DEBUG
		printf("End create kernel applyBeeman1&2\n");
#endif
	}
	if(molecularDynamics->integratorType==STOCHASTIC) 
	{
		molecularDynamics->frictionBufferCL = clCreateBuffer(clProp.context, CL_MEM_READ_WRITE, sizeof(cl_float4) * forceField->molecule.nAtoms*3, NULL, NULL);
#ifdef DEBUG
		printf("Begin applyStochastic\n");
#endif
		molecularDynamics->applyStochastic1 = clCreateKernel(molecularDynamics->programMD, "applyStochastic1", &err);
		clSetKernelArg(molecularDynamics->applyStochastic1, 0, sizeof(cl_mem),  &molecularDynamics->frictionBufferCL);
		clSetKernelArg(molecularDynamics->applyStochastic1, 1, sizeof(cl_mem),  &molecularDynamics->aCL);
		clSetKernelArg(molecularDynamics->applyStochastic1, 2, sizeof(cl_mem),  &forceField->atomsCL);
		clSetKernelArg(molecularDynamics->applyStochastic1, 3, sizeof(cl_float), &dtf);
		clSetKernelArg(molecularDynamics->applyStochastic1, 4, sizeof(cl_int), &forceField->molecule.nAtoms);
#ifdef DEBUG
		printf("End applyStochastic1\n");
#endif

		molecularDynamics->applyStochastic2 = clCreateKernel(molecularDynamics->programMD, "applyStochastic2", &err);
		clSetKernelArg(molecularDynamics->applyStochastic2, 0, sizeof(cl_mem),   &molecularDynamics->frictionBufferCL);
		clSetKernelArg(molecularDynamics->applyStochastic2, 1, sizeof(cl_mem),   &molecularDynamics->aCL);
		clSetKernelArg(molecularDynamics->applyStochastic2, 2, sizeof(cl_mem),   &forceField->atomsCL);
		clSetKernelArg(molecularDynamics->applyStochastic2, 3, sizeof(cl_float), &dtf);
		clSetKernelArg(molecularDynamics->applyStochastic2, 4, sizeof(cl_int), &forceField->molecule.nAtoms);
#ifdef DEBUG
		printf("End applyStochastic2\n");
#endif
	}
	initRandomNumberGenerator(molecularDynamics, (unsigned int)(time(NULL)));
	if(molecularDynamics->thermostat==ANDERSEN) 
	{
		molecularDynamics->andersen = clCreateKernel(molecularDynamics->programMD, "andersen", &err);
		clSetKernelArg(molecularDynamics->andersen, 0, sizeof(cl_mem),   &molecularDynamics->randomsCL);
		clSetKernelArg(molecularDynamics->andersen, 1, sizeof(cl_mem),   &forceField->atomsCL);
		clSetKernelArg(molecularDynamics->andersen, 2, sizeof(cl_mem),   &molecularDynamics->nRandomsCL);
		clSetKernelArg(molecularDynamics->andersen, 3, sizeof(cl_int), &forceField->molecule.nAtoms);
		clSetKernelArg(molecularDynamics->andersen, 4, sizeof(cl_float), &temp);
#ifdef DEBUG
		printf("End ANDERSEN\n");
#endif
	}

#ifdef DEBUG
	printf("Begin updateGeometryVelocitiesCL--------------->\n");
#endif
	updateGeometryVelocitiesCL(forceField,&forceField->molecule);
	{
		size_t global = molecularDynamics->forceField->molecule.nAtoms;
		cl_int err = clEnqueueNDRangeKernel(clProp.command_queue, molecularDynamics->initAccelarations, 1, NULL, &global, NULL, 0, NULL, NULL);
		if(err != CL_SUCCESS)
		{
			printErrorCLRun(err);
        		printf("I cannot execute initAccelarations\n");
			exit(1);
		}
		clFinish(clProp.command_queue);
	}
	if(molecularDynamics->integratorType==STOCHASTIC) 
	{
		printf("---------->frictionf=%f\n",frictionf);
		molecularDynamics->setFrictionalAndRandomForce = clCreateKernel(molecularDynamics->programMD, "setFrictionalAndRandomForce", &err);
		clSetKernelArg(molecularDynamics->setFrictionalAndRandomForce, 0, sizeof(cl_mem), &molecularDynamics->frictionBufferCL);
		clSetKernelArg(molecularDynamics->setFrictionalAndRandomForce, 1, sizeof(cl_mem), &forceField->atomsCL);
		clSetKernelArg(molecularDynamics->setFrictionalAndRandomForce, 2, sizeof(cl_mem), &molecularDynamics->randomsCL);
		clSetKernelArg(molecularDynamics->setFrictionalAndRandomForce, 3, sizeof(cl_int), &molecularDynamics->nRandomsCL);
		clSetKernelArg(molecularDynamics->setFrictionalAndRandomForce, 4, sizeof(cl_int), &forceField->molecule.nAtoms);
		clSetKernelArg(molecularDynamics->setFrictionalAndRandomForce, 5, sizeof(cl_float), &temp);
		clSetKernelArg(molecularDynamics->setFrictionalAndRandomForce, 6, sizeof(cl_float), &dtf);
		clSetKernelArg(molecularDynamics->setFrictionalAndRandomForce, 7, sizeof(cl_float), &frictionf);
#ifdef DEBUG
		printf("End create kernel setFrictionalAndRandomForce\n");
#endif
	}
#ifdef DEBUG
	printf("Begin calcGrad2\n");
#endif
	molecularDynamics->forceField->klass->calculateGradient(molecularDynamics->forceField);
#ifdef DEBUG
	printf("End calcGrad2\n");
#endif
	{
		size_t global = molecularDynamics->forceField->molecule.nAtoms;
		cl_int err = clEnqueueNDRangeKernel(clProp.command_queue, molecularDynamics->computeAccelarationsFromGradients, 1, NULL, &global, NULL, 0, NULL, NULL);
		if(err != CL_SUCCESS)
		{
			printErrorCLRun(err);
        		printf("I cannot execute computeAccelarationsFromGradients\n");
			exit(1);
		}
		clFinish(clProp.command_queue);
	}
	if(molecularDynamics->aold)
	{
		size_t global = molecularDynamics->forceField->molecule.nAtoms;
		cl_int err = clEnqueueNDRangeKernel(clProp.command_queue, molecularDynamics->copyAccelarations, 1, NULL, &global, NULL, 0, NULL, NULL);
		if(err != CL_SUCCESS)
		{
			printErrorCLRun(err);
        		printf("I cannot execute copyAccelarations\n");
			exit(1);
		}
		clFinish(clProp.command_queue);
	}
	if(molecularDynamics->forceField->molecule.constraints!=NOCONSTRAINTS)
	{
		molecularDynamics->oldPositionBufferCL = clCreateBuffer(clProp.context, CL_MEM_READ_WRITE, sizeof(cl_float4) * forceField->molecule.nAtoms, NULL, NULL);
		molecularDynamics->copyPositions = clCreateKernel(molecularDynamics->programMD, "copyPositions", &err);
		clSetKernelArg(molecularDynamics->copyPositions, 0, sizeof(cl_mem), &molecularDynamics->oldPositionBufferCL);
		clSetKernelArg(molecularDynamics->copyPositions, 1, sizeof(cl_mem), &forceField->atomsCL);
		clSetKernelArg(molecularDynamics->copyPositions, 2, sizeof(cl_int), &forceField->molecule.nAtoms);

		molecularDynamics->initRattles = clCreateKernel(molecularDynamics->programMD, "initRattles", &err);
		clSetKernelArg(molecularDynamics->initRattles, 0, sizeof(cl_mem), &forceField->rattledeltaPositionBufferCL);
		clSetKernelArg(molecularDynamics->initRattles, 1, sizeof(cl_mem), &forceField->rattledeltaVelocityBufferCL);
		clSetKernelArg(molecularDynamics->initRattles, 2, sizeof(cl_mem), &forceField->rattleUpdateBufferCL);
		clSetKernelArg(molecularDynamics->initRattles, 3, sizeof(cl_mem), &forceField->rattleMovedBufferCL);
		clSetKernelArg(molecularDynamics->initRattles, 4, sizeof(cl_int), &forceField->molecule.nAtoms);
		clSetKernelArg(molecularDynamics->initRattles, 5, sizeof(cl_int), &forceField->nBlockRattleBuffer);
#ifdef DEBUG
		printf("End create kernel initRattles\n");
#endif

		molecularDynamics->updateRattle1 = clCreateKernel(molecularDynamics->programMD, "updateRattle1", &err);
		clSetKernelArg(molecularDynamics->updateRattle1, 0, sizeof(cl_mem), &forceField->rattledeltaPositionBufferCL);
		clSetKernelArg(molecularDynamics->updateRattle1, 1, sizeof(cl_mem), &forceField->rattledeltaVelocityBufferCL);
		clSetKernelArg(molecularDynamics->updateRattle1, 2, sizeof(cl_mem), &forceField->rattleUpdateBufferCL);
		clSetKernelArg(molecularDynamics->updateRattle1, 3, sizeof(cl_mem), &forceField->rattleMovedBufferCL);
		clSetKernelArg(molecularDynamics->updateRattle1, 4, sizeof(cl_mem), &forceField->atomsCL);
		clSetKernelArg(molecularDynamics->updateRattle1, 5, sizeof(cl_int), &forceField->molecule.nAtoms);
		clSetKernelArg(molecularDynamics->updateRattle1, 6, sizeof(cl_int), &forceField->nBlockRattleBuffer);
#ifdef DEBUG
		printf("End create kernel updateRattle1\n");
#endif

		molecularDynamics->testDoneRattle = clCreateKernel(molecularDynamics->programMD, "testDoneRattle", &err);
		clSetKernelArg(molecularDynamics->testDoneRattle, 0, sizeof(cl_mem), &forceField->rattleMovedBufferCL);
		clSetKernelArg(molecularDynamics->testDoneRattle, 1, sizeof(cl_mem), &forceField->rattleDoneBufferCL);
		clSetKernelArg(molecularDynamics->testDoneRattle, 2, sizeof(cl_int), &forceField->molecule.nAtoms);
#ifdef DEBUG
		printf("End create kernel testDoneRattle\n");
#endif

		molecularDynamics->applyRattle1 = clCreateKernel(molecularDynamics->programMD, "applyRattle1", &err);
		clSetKernelArg(molecularDynamics->applyRattle1, 0, sizeof(cl_mem), &molecularDynamics->oldPositionBufferCL);
		clSetKernelArg(molecularDynamics->applyRattle1, 1, sizeof(cl_mem), &forceField->rattleConstraintsIndexCL);
		clSetKernelArg(molecularDynamics->applyRattle1, 2, sizeof(cl_mem), &forceField->rattleConstraintsTermsCL);
		clSetKernelArg(molecularDynamics->applyRattle1, 3, sizeof(cl_mem), &forceField->rattledeltaPositionBufferCL);
		clSetKernelArg(molecularDynamics->applyRattle1, 4, sizeof(cl_mem), &forceField->rattledeltaVelocityBufferCL);
		clSetKernelArg(molecularDynamics->applyRattle1, 5, sizeof(cl_mem), &forceField->rattleUpdateBufferCL);
		clSetKernelArg(molecularDynamics->applyRattle1, 6, sizeof(cl_mem), &forceField->rattleMovedBufferCL);
		clSetKernelArg(molecularDynamics->applyRattle1, 7, sizeof(cl_mem), &forceField->atomsCL);
		clSetKernelArg(molecularDynamics->applyRattle1, 8, sizeof(cl_int), &forceField->molecule.nAtoms);
		clSetKernelArg(molecularDynamics->applyRattle1, 9, sizeof(cl_int), &forceField->numberOfRattleConstraintsTerms);
		clSetKernelArg(molecularDynamics->applyRattle1, 10, sizeof(cl_int), &dtf);
#ifdef DEBUG
		printf("End create kernel applyRattle1\n");
#endif


		molecularDynamics->updateRattle2 = clCreateKernel(molecularDynamics->programMD, "updateRattle2", &err);
		clSetKernelArg(molecularDynamics->updateRattle2, 0, sizeof(cl_mem), &forceField->rattledeltaVelocityBufferCL);
		clSetKernelArg(molecularDynamics->updateRattle2, 1, sizeof(cl_mem), &forceField->rattleUpdateBufferCL);
		clSetKernelArg(molecularDynamics->updateRattle2, 2, sizeof(cl_mem), &forceField->rattleMovedBufferCL);
		clSetKernelArg(molecularDynamics->updateRattle2, 3, sizeof(cl_mem), &forceField->atomsCL);
		clSetKernelArg(molecularDynamics->updateRattle2, 4, sizeof(cl_int), &forceField->molecule.nAtoms);
		clSetKernelArg(molecularDynamics->updateRattle2, 5, sizeof(cl_int), &forceField->nBlockRattleBuffer);
#ifdef DEBUG
		printf("End create kernel updateRattle2\n");
#endif

		molecularDynamics->applyRattle2 = clCreateKernel(molecularDynamics->programMD, "applyRattle2", &err);
		clSetKernelArg(molecularDynamics->applyRattle2, 0, sizeof(cl_mem), &forceField->rattleConstraintsIndexCL);
		clSetKernelArg(molecularDynamics->applyRattle2, 1, sizeof(cl_mem), &forceField->rattleConstraintsTermsCL);
		clSetKernelArg(molecularDynamics->applyRattle2, 2, sizeof(cl_mem), &forceField->rattledeltaVelocityBufferCL);
		clSetKernelArg(molecularDynamics->applyRattle2, 3, sizeof(cl_mem), &forceField->rattleUpdateBufferCL);
		clSetKernelArg(molecularDynamics->applyRattle2, 4, sizeof(cl_mem), &forceField->rattleMovedBufferCL);
		clSetKernelArg(molecularDynamics->applyRattle2, 5, sizeof(cl_mem), &forceField->atomsCL);
		clSetKernelArg(molecularDynamics->applyRattle2, 6, sizeof(cl_int), &forceField->molecule.nAtoms);
		clSetKernelArg(molecularDynamics->applyRattle2, 7, sizeof(cl_int), &forceField->numberOfRattleConstraintsTerms);
		clSetKernelArg(molecularDynamics->applyRattle2, 8, sizeof(cl_int), &dtf);
#ifdef DEBUG
		printf("End create kernel applyRattle2\n");
#endif
	}

/* TO DO for rattle
	for ( i = 0; i < molecularDynamics->numberOfAtoms; i++)
	{
		if(molecularDynamics->coordinatesOld)
		       	for ( j = 0; j < 3; j++) molecularDynamics->coordinatesOld[i][j]  = molecularDynamics->forceField->molecule.atoms[i].coordinates[j];
		if(molecularDynamics->moved) molecularDynamics->moved[i] = FALSE; 
		if(molecularDynamics->update) molecularDynamics->update[i] = FALSE; 
	}
*/
#endif

#ifdef DEBUG
	printf("End initMD--------------->\n");
#endif
}
/*********************************************************************************/
static void rescaleVelocities(MolecularDynamics* molecularDynamics)
{
	/* berendsen(molecularDynamics);*/
	if(molecularDynamics->temperature<=0) return;
	scaleV(molecularDynamics);
	resetQTB(molecularDynamics);
	resetLangevin(molecularDynamics);
}
/*********************************************************************************/
#ifdef ENABLE_CL
static void scaleCLVelocities(MolecularDynamics* molecularDynamics, double scale)
{
	{
		size_t global = 1;
		global =  molecularDynamics->numberOfAtoms;
		cl_float s = (cl_float)scale;
		cl_int err;
		CLProp clProp = getCLProp();

		clSetKernelArg(molecularDynamics->scaleVelocities, 2, sizeof(cl_float), &s);
		err = clEnqueueNDRangeKernel(clProp.command_queue, molecularDynamics->scaleVelocities, 1, NULL, &global, NULL, 0, NULL, NULL);
		if(err != CL_SUCCESS)
		{
			printErrorCLRun(err);
        		printf("I cannot execute scaleVelocities\n");
			exit(1);
		}
		clFinish(clProp.command_queue);
	}
//	printf("End scaleCLVelocities\n");
}
#endif
/*********************************************************************************/
static void scaleV(MolecularDynamics* molecularDynamics)
{
	double kelvin = 0;
	int nfree = molecularDynamics->forceField->molecule.nFree;
	double scale = 1.0;
#ifndef ENABLE_CL
	//double mass = 1.0;
	double ekin = 0;
	int i;
	int j;
#endif
	if(molecularDynamics->temperature<=0) return;
	if(nfree<1) return;
#ifdef ENABLE_CL
	{
		size_t global = 1;
		cl_float ekin;
		ekin = getEKin(molecularDynamics);
		kelvin = 2*ekin / ( nfree * Kb);
		scale = sqrt(molecularDynamics->temperature/kelvin);
		global =  molecularDynamics->numberOfAtoms;

#ifdef DEBUG
		printf("temp = %f kelvin = %f scale = %f\n",molecularDynamics->temperature, kelvin, scale);
		printf("===========> Call scaleVelocities\n");
#endif
		scaleCLVelocities(molecularDynamics, scale);
	}
#else 
	ekin = getEKin(molecularDynamics);
	kelvin = 2*ekin / ( nfree * Kb);

	scale = sqrt(molecularDynamics->temperature/kelvin);
#ifdef DEBUG
	printf("temp = %f kelvin = %f scale = %f\n",molecularDynamics->temperature, kelvin, scale);
#endif
	for ( i = 0; i < molecularDynamics->numberOfAtoms; i++)
		if(molecularDynamics->forceField->molecule.atoms[i].variable)
		for ( j = 0; j < 3; j++)
			molecularDynamics->forceField->molecule.atoms[i].velocity[j] *= scale;
#endif
	removeTranslationAndRotation(molecularDynamics);

/*
	ekin = getEKin(molecularDynamics);
	kelvin = 2*ekin / ( nfree * Kb);
	scale = sqrt(molecularDynamics->temperature/kelvin);
	printf("REM temp = %f kelvin = %f scale = %f\n",molecularDynamics->temperature, kelvin, scale);
*/

}
/*********************************************************************************/
static void berendsen(MolecularDynamics* molecularDynamics)
{
	double kelvin = 0;
	int nfree = molecularDynamics->forceField->molecule.nFree;
	double scale = 1.0;
	double dt = molecularDynamics->dt;
	double tautemp = 1.0/(molecularDynamics->collide)*1000*fsInAKMA;
#ifndef ENABLE_CL
	//double mass = 1.0;
	double ekin = 0;
	int i;
	int j;
#endif
	if(molecularDynamics->temperature<=0) return;
	if(nfree<1) return;
#ifdef ENABLE_CL
	{
		size_t global = 1;
		cl_float ekin;
		ekin = getEKin(molecularDynamics);
		kelvin = 2*ekin / ( nfree * Kb);
		scale = sqrt(1.0 + (dt/tautemp)*(molecularDynamics->temperature/kelvin-1.0));
		global =  molecularDynamics->numberOfAtoms;

#ifdef DEBUG
		printf("temp = %f kelvin = %f scale = %f\n",molecularDynamics->temperature, kelvin, scale);
		printf("===========> Call scaleVelocities\n");
#endif
		scaleCLVelocities(molecularDynamics, scale);
	}
#else 
/*
	for ( i = 0; i < molecularDynamics->numberOfAtoms; i++)
	{
		mass = molecularDynamics->forceField->molecule.atoms[i].mass;
		for ( j = 0; j < 3; j++)
			ekin += molecularDynamics->forceField->molecule.atoms[i].velocity[j]*molecularDynamics->forceField->molecule.atoms[i].velocity[j]*
				mass;
	}
	//ekin /= 2;
	//kelvin = 2* ekin / ( nfree * Kb);
	kelvin = ekin / ( nfree * Kb);
*/
	ekin = getEKin(molecularDynamics);
	kelvin = 2*ekin / ( nfree * Kb);

	scale = sqrt(1.0 + (dt/tautemp)*(molecularDynamics->temperature/kelvin-1.0));
#ifdef DEBUG
	printf("temp = %f kelvin = %f scale = %f\n",molecularDynamics->temperature, kelvin, scale);
#endif
	for ( i = 0; i < molecularDynamics->numberOfAtoms; i++)
		if(molecularDynamics->forceField->molecule.atoms[i].variable)
		for ( j = 0; j < 3; j++)
			molecularDynamics->forceField->molecule.atoms[i].velocity[j] *= scale;
#endif
	removeTranslationAndRotation(molecularDynamics);
}
/*********************************************************************************/
static void andersen(MolecularDynamics* molecularDynamics)
{
	double tau = 1.0/molecularDynamics->collide*1000*fsInAKMA; /* in fs */
	double rate;
	if(molecularDynamics->temperature<=0) return;
	if(molecularDynamics->numberOfAtoms<1) return;

	rate = molecularDynamics->dt / tau;
	rate /= pow(molecularDynamics->nvariables,2.0/3.0);

#ifdef ENABLE_CL
	size_t global =  molecularDynamics->forceField->molecule.nAtoms;
	cl_int err;
	CLProp clProp = getCLProp();
	cl_float r = rate;
	cl_float temp = molecularDynamics->temperature;
	resetRandomNumbers(molecularDynamics);
	clSetKernelArg(molecularDynamics->andersen, 4, sizeof(cl_float), &temp);
	clSetKernelArg(molecularDynamics->andersen, 5, sizeof(cl_float), &r);
	err = clEnqueueNDRangeKernel(clProp.command_queue, molecularDynamics->andersen, 1, NULL, &global, NULL, 0, NULL, NULL);
	if(err != CL_SUCCESS)
	{
		printErrorCLRun(err);
       		printf("I cannot execute andersen\n");
		exit(1);
	}
	clFinish(clProp.command_queue);
#else
	int i;
	/* printf("------->rate=%f\n",rate);*/
	for ( i = 0; i < molecularDynamics->numberOfAtoms; i++)
	{
		double trial = drandom();
		double m = molecularDynamics->forceField->molecule.atoms[i].mass;
		if(molecularDynamics->forceField->molecule.atoms[i].variable)
		if(trial<rate)
		{
	/*
			double speed = maxwel(
					molecularDynamics->forceField->molecule.atoms[i].mass,
					molecularDynamics->temperature
					);
			getRandVect(speed, molecularDynamics->forceField->molecule.atoms[i].velocity);
	*/
			double speed = sqrt(Kb* molecularDynamics->temperature/m);
                	double pnorm = normal();
			molecularDynamics->forceField->molecule.atoms[i].velocity[0] = pnorm*speed;
                	pnorm = normal();
			molecularDynamics->forceField->molecule.atoms[i].velocity[1] = pnorm*speed;
                	pnorm = normal();
			molecularDynamics->forceField->molecule.atoms[i].velocity[2] = pnorm*speed;
		}
	}
#endif
}
/*********************************************************************************/
static void bussi(MolecularDynamics* molecularDynamics)
{
	int nfree = molecularDynamics->forceField->molecule.nFree;
	double scale = 1.0;
	double dt = molecularDynamics->dt;
	double tautemp = 1.0/(molecularDynamics->collide)*1000*fsInAKMA;
        double c = exp(-dt/tautemp);
	double ekin = getEKin(molecularDynamics);
	double kelvin = 2*ekin / ( nfree * Kb);
	double d = (1.0-c) * (molecularDynamics->temperature/kelvin) / (nfree);
	double r = normal ();
	double si = 0.0;
	double s = 0.0;
	int i;
#ifndef ENABLE_CL
	int j;
#endif
	if(molecularDynamics->temperature<=0) return;
	if(nfree<1) return;
        for(i=0;i<nfree-1;i++)
	{
            si = normal ();
            s += si*si;
	}
	scale = c + (s+r*r)*d + 2.0*r*sqrt(c*d);
	scale = sqrt(scale);
	if (r+sqrt(c/d)<0)  scale = -scale;
#ifdef ENABLE_CL
	{
		size_t global = 1;
		cl_float s;
		cl_int err;
		CLProp clProp = getCLProp();
		global =  molecularDynamics->numberOfAtoms;
		s = (cl_float)scale;

#ifdef DEBUG
		printf("temp = %f kelvin = %f scale = %f\n",molecularDynamics->temperature, kelvin, scale);
		printf("===========> Call scaleVelocities\n");
#endif

		clSetKernelArg(molecularDynamics->scaleVelocities, 2, sizeof(cl_float), &s);
		err = clEnqueueNDRangeKernel(clProp.command_queue, molecularDynamics->scaleVelocities, 1, NULL, &global, NULL, 0, NULL, NULL);
		if(err != CL_SUCCESS)
		{
			printErrorCLRun(err);
        		printf("I cannot execute scaleVelocities\n");
			exit(1);
		}
		clFinish(clProp.command_queue);
	}
#else 
#ifdef DEBUG
	printf("temp = %f kelvin = %f scale = %f\n",molecularDynamics->temperature, kelvin, scale);
#endif
	for ( i = 0; i < molecularDynamics->numberOfAtoms; i++)
		if(molecularDynamics->forceField->molecule.atoms[i].variable)
		for ( j = 0; j < 3; j++)
			molecularDynamics->forceField->molecule.atoms[i].velocity[j] *= scale;
#endif
	removeTranslationAndRotation(molecularDynamics);
}
/*********************************************************************************/
static void nose_hoover(MolecularDynamics* molecularDynamics)
{
	int nfree = molecularDynamics->forceField->molecule.nFree;
	double scale = 1.0;
	double ekin = getEKin(molecularDynamics);
	double kT = Kb* molecularDynamics->temperature;
	int i,j;
	if(molecularDynamics->temperature<=0) return;
	if(nfree<1) return;
	molecularDynamics->gNH[1] = (molecularDynamics->qNH[0]*molecularDynamics->vNH[0]*molecularDynamics->vNH[0]-kT) / molecularDynamics->qNH[1];
	//printf("gNH = %f\n",molecularDynamics->gNH[1]);
	molecularDynamics->vNH[1] = molecularDynamics->vNH[1] + molecularDynamics->gNH[1]*molecularDynamics->dt_4;
	molecularDynamics->vNH[0] = molecularDynamics->vNH[0] * exp(-molecularDynamics->vNH[1]*molecularDynamics->dt_8);
	molecularDynamics->gNH[0] = (2.0*ekin-nfree*kT) / molecularDynamics->qNH[0];
	molecularDynamics->vNH[0] = molecularDynamics->vNH[0] + molecularDynamics->gNH[0]*molecularDynamics->dt_4;
	molecularDynamics->vNH[0] = molecularDynamics->vNH[0] * exp(-molecularDynamics->vNH[1]*molecularDynamics->dt_8);
	molecularDynamics->xNH[0] = molecularDynamics->xNH[0] + molecularDynamics->vNH[0]*molecularDynamics->dt_2;
	molecularDynamics->xNH[1] = molecularDynamics->xNH[1] + molecularDynamics->vNH[1]*molecularDynamics->dt_2;
	//printf("vnH0 = %f\n",molecularDynamics->vNH[0]);
	scale = exp(-molecularDynamics->vNH[0]*molecularDynamics->dt_2);
	//printf("scale = %f\n",scale);
#ifdef ENABLE_CL
	{
		size_t global = 1;
		cl_float s;
		cl_int err;
		CLProp clProp = getCLProp();
		global =  molecularDynamics->numberOfAtoms;
		s = (cl_float)scale;

#ifdef DEBUG
		printf("===========> Call scaleVelocities\n");
#endif
		clSetKernelArg(molecularDynamics->scaleVelocities, 2, sizeof(cl_float), &s);
		err = clEnqueueNDRangeKernel(clProp.command_queue, molecularDynamics->scaleVelocities, 1, NULL, &global, NULL, 0, NULL, NULL);
		if(err != CL_SUCCESS)
		{
			printErrorCLRun(err);
        		printf("I cannot execute scaleVelocities\n");
			exit(1);
		}
		clFinish(clProp.command_queue);
	}
#else 
	for ( i = 0; i < molecularDynamics->numberOfAtoms; i++)
		if(molecularDynamics->forceField->molecule.atoms[i].variable)
		for ( j = 0; j < 3; j++)
			molecularDynamics->forceField->molecule.atoms[i].velocity[j] *= scale;
#endif
	ekin = ekin * scale * scale;
	molecularDynamics->vNH[0] = molecularDynamics->vNH[0] * exp(-molecularDynamics->vNH[1]*molecularDynamics->dt_8);
	molecularDynamics->gNH[0] = (2.0*ekin-nfree*kT) /  molecularDynamics->qNH[0];
	molecularDynamics->vNH[0] = molecularDynamics->vNH[0] + molecularDynamics->gNH[0]*molecularDynamics->dt_4;
	molecularDynamics->vNH[0] = molecularDynamics->vNH[0] * exp(-molecularDynamics->vNH[1]*molecularDynamics->dt_8);
	molecularDynamics->gNH[1] = ( molecularDynamics->qNH[0]*molecularDynamics->vNH[0]*molecularDynamics->vNH[0]-kT) /  molecularDynamics->qNH[1];
	molecularDynamics->vNH[1] = molecularDynamics->vNH[1] + molecularDynamics->gNH[1]*molecularDynamics->dt_4;

	removeTranslationAndRotation(molecularDynamics);

}
/*********************************************************************************/
static void newAccelaration(MolecularDynamics* molecularDynamics)
{
#ifdef DEBUG
	TimerType timer;
	TimerType timer2;
        timer_init(timer);
       	timer_start( timer );
        timer_init(timer2);
       	timer_start( timer2 );
#endif
#ifdef ENABLE_CL
	CLProp clProp=getCLProp();

#ifdef DEBUG
        printf("Begin newAcc\n");
#endif
	molecularDynamics->forceField->klass->calculateGradient(molecularDynamics->forceField);
#ifdef DEBUG
       	timer_stop(timer2);
        printf("time newAccCalcG (s) = %f\n", timer_get(timer2)*1e-6);
#endif
	if(molecularDynamics->aold)
	{
		size_t global = molecularDynamics->forceField->molecule.nAtoms;
		cl_int err = clEnqueueNDRangeKernel(clProp.command_queue, molecularDynamics->copyAccelarations, 1, NULL, &global, NULL, 0, NULL, NULL);
		if(err != CL_SUCCESS)
		{
			printErrorCLRun(err);
        		printf("I cannot execute copyAccelarations\n");
			exit(1);
		}
		clFinish(clProp.command_queue);
	}
	{
		size_t global = molecularDynamics->forceField->molecule.nAtoms;
		cl_int err = clEnqueueNDRangeKernel(clProp.command_queue, molecularDynamics->computeAccelarationsFromGradients, 1, NULL, &global, NULL, 0, NULL, NULL);
		if(err != CL_SUCCESS)
		{
			printErrorCLRun(err);
        		printf("I cannot execute computeAccelarationsFromGradients\n");
			exit(1);
		}
		clFinish(clProp.command_queue);
	}
#else
	int i;
	int j;
	molecularDynamics->forceField->klass->calculateGradient(molecularDynamics->forceField);
#ifdef DEBUG
       	timer_stop(timer2);
        printf("time newAccCalcG (s) = %f\n", timer_get(timer2)*1e-6);
#endif
	for ( i = 0; i < molecularDynamics->numberOfAtoms; i++)
	{
		double m = molecularDynamics->forceField->molecule.atoms[i].mass;
		if(molecularDynamics->aold)
			for ( j = 0; j < 3; j++)
				molecularDynamics->aold[i][j]  = molecularDynamics->a[i][j];

		for ( j = 0; j < 3; j++)
			molecularDynamics->a[i][j] = -molecularDynamics->forceField->molecule.atoms[i].gradient[j]/m;
	}
#endif
#ifdef DEBUG
       	timer_stop(timer);
        printf("time newAcc (s) = %f\n", timer_get(timer)*1e-6);
#endif
}
/*********************************************************************************/
static void computeEnergies(MolecularDynamics* molecularDynamics)
{
	molecularDynamics->kineticEnergy = getEKin(molecularDynamics);
	//molecularDynamics->potentialEnergy = molecularDynamics->forceField->klass->calculateEnergyTmp(
	//	      molecularDynamics->forceField,&molecularDynamics->forceField->molecule);
	molecularDynamics->potentialEnergy = molecularDynamics->forceField->molecule.potentialEnergy;
	molecularDynamics->totalEnergy = molecularDynamics->kineticEnergy + molecularDynamics->potentialEnergy;
	molecularDynamics->kelvin = getKelvin(molecularDynamics);
}
/*********************************************************************************/
static void applyThermostat(MolecularDynamics* molecularDynamics)
{
	if(molecularDynamics->integratorType == STOCHASTIC) return;
	if(molecularDynamics->integratorType == QTB) return;
	if(molecularDynamics->integratorType == LANGEVIN) return;
	if(molecularDynamics->thermostat == ANDERSEN) andersen(molecularDynamics);
	if(molecularDynamics->thermostat == BERENDSEN) berendsen(molecularDynamics);
	if(molecularDynamics->thermostat == BUSSI) bussi(molecularDynamics);
}
/*********************************************************************************/
static void applyOneStep(MolecularDynamics* molecularDynamics, int iStep)
{
#ifdef DEBUG
	TimerType timer;
        timer_init(timer);
       	timer_start( timer );
#endif
	//printf("molecularDynamics->integratorType=%d Tuckerman=%d\n",molecularDynamics->integratorType,MARTYNATUCKERMAN);
	if(molecularDynamics->integratorType==VERLET) applyVerlet(molecularDynamics);
	else if(molecularDynamics->integratorType==BEEMAN) applyBeeman(molecularDynamics);
	else if(molecularDynamics->integratorType==STOCHASTIC) applyStochastic(molecularDynamics);
	else if(molecularDynamics->integratorType==MARTYNATUCKERMAN) applyMartynaTuckerman(molecularDynamics);
	else if(molecularDynamics->integratorType==LANGEVIN) {
		updateLangevin(molecularDynamics);
		applyLangevin(molecularDynamics);
	}
	else {
		if((iStep+1)%molecularDynamics->M==0) 
				updateQTB(molecularDynamics);
		applyQTB(molecularDynamics);
	}
#ifdef DEBUG
       	timer_stop(timer);
        printf("time applyOneStape (s) = %f\n", timer_get(timer)*1e-6);
#endif
	computeEnergies(molecularDynamics);
}
/*********************************************************************************/
static void applyRattleFirstPortion(MolecularDynamics* molecularDynamics)
{
#ifdef ENABLE_CL
	ForceField* forceField = molecularDynamics->forceField;
	CLProp clProp=getCLProp();
	cl_int done=1;
	int nIter = 0;
	int maxIter = 1000;
	cl_int err;
	size_t global = forceField->molecule.nAtoms*forceField->nBlockRattleBuffer;
	err = clEnqueueNDRangeKernel(clProp.command_queue, molecularDynamics->initRattles, 1, NULL, &global, NULL, 0, NULL, NULL);
	if(err != CL_SUCCESS)
	{
		printErrorCLRun(err);
       		printf("I cannot execute initRattles,\n");
		exit(1);
	}
	clFinish(clProp.command_queue);
#ifdef DEBUG
	printf("applyRattleFirstPortion : end initRattles\n");
#endif
	do{
		nIter++;
		//printf("applyRattleFirstPortion : nIter =%d\n",nIter);
		global =  molecularDynamics->forceField->numberOfRattleConstraintsTerms;
		err = clEnqueueNDRangeKernel(clProp.command_queue, molecularDynamics->applyRattle1, 1, NULL, &global, NULL, 0, NULL, NULL);
		global = forceField->molecule.nAtoms;
		err = clEnqueueNDRangeKernel(clProp.command_queue, molecularDynamics->updateRattle1, 1, NULL, &global, NULL, 0, NULL, NULL);
		global = 1;
		err = clEnqueueNDRangeKernel(clProp.command_queue, molecularDynamics->testDoneRattle, 1, NULL, &global, NULL, 0, NULL, NULL);
		clEnqueueReadBuffer(clProp.command_queue, forceField->rattleDoneBufferCL, CL_TRUE, 0, sizeof(cl_int)*1, &done, 0, NULL, NULL);
		//printf("applyRattleFirstPortion : Done =%d\n",done);
	}while(done<1 && nIter<maxIter);
	if(nIter>=maxIter)
	{
		printf(("Rattle first portion : Warning, distance constraints not satisfied\n"));
		exit(1);
	}
#else
	int i;
	int k;
	int maxIter = 1000;
	double omega = 1.0; 
	double tolerance = 1e-6;
	boolean done = FALSE;
	int nIter = 0;
	int a1 = 0;
	int a2 = 0;
	double r2ij;
	double dot;
	double invMass1;
	double invMass2;
	double delta;
	double term = 0;
	double terms[3];
	double d;
	Molecule* m = &molecularDynamics->forceField->molecule;
	ForceField* forceField = molecularDynamics->forceField;
	double deltaMax = 0;

	if(forceField->molecule.constraints==NOCONSTRAINTS) return;
	for (i = 0; i < molecularDynamics->numberOfAtoms; i++)
	{
			molecularDynamics->moved[i] = molecularDynamics->forceField->molecule.atoms[i].variable;
			//molecularDynamics->moved[i] = TRUE;
			molecularDynamics->update[i] = FALSE;
	}
	/* maxIter *= molecularDynamics->forceField->numberOfRattleConstraintsTerms;*/
	do{
		nIter++;
		done=TRUE;
		deltaMax = 0;
		for (i = 0; i < molecularDynamics->forceField->molecule.numberOfRattleConstraintsTerms; i++)
		{
			a1 = (int)molecularDynamics->forceField->molecule.rattleConstraintsTerms[0][i];
			a2 = (int)molecularDynamics->forceField->molecule.rattleConstraintsTerms[1][i];
			if( !molecularDynamics->moved[a1] && !molecularDynamics->moved[a2] ) continue;
			r2ij = 0;
			for (k=0;k<3;k++)
			{
				d = m->atoms[a2].coordinates[k]-m->atoms[a1].coordinates[k];
				r2ij +=d*d;
			}
			delta = molecularDynamics->forceField->molecule.rattleConstraintsTerms[2][i]-r2ij;
			/* if(fabs(delta)<=tolerance) continue;*/
			if(r2ij>0 && fabs(delta/r2ij)<=tolerance) continue;
			if(deltaMax<fabs(delta)) deltaMax = fabs(delta);
			done = FALSE;
			molecularDynamics->update[a1] = TRUE;
			molecularDynamics->update[a2] = TRUE;
			/* here : rattle image for PBC, not yet implemented */
			dot = 0;
			for (k=0;k<3;k++)
			{
				d = m->atoms[a2].coordinates[k]-m->atoms[a1].coordinates[k];
				dot +=d*(molecularDynamics->coordinatesOld[a2][k]-molecularDynamics->coordinatesOld[a1][k]);
			}
			invMass1 = 1/m->atoms[a1].mass;
			invMass2 = 1/m->atoms[a2].mass;
		        term = omega*delta / (2.0*(invMass1+invMass2)*dot);
			for (k=0;k<3;k++)
			{
				terms[k] = (molecularDynamics->coordinatesOld[a2][k]-molecularDynamics->coordinatesOld[a1][k])*term;
			}
			for (k=0;k<3;k++) m->atoms[a1].coordinates[k] -= terms[k]*invMass1;
			for (k=0;k<3;k++) m->atoms[a2].coordinates[k] += terms[k]*invMass2;

			invMass1 /= molecularDynamics->dt;
			invMass2 /= molecularDynamics->dt;
			for (k=0;k<3;k++) molecularDynamics->forceField->molecule.atoms[a1].velocity[k] -= terms[k]*invMass1;
			for (k=0;k<3;k++) molecularDynamics->forceField->molecule.atoms[a2].velocity[k] += terms[k]*invMass2;
		}
		for (i = 0; i < molecularDynamics->numberOfAtoms; i++)
		{
			molecularDynamics->moved[i] = molecularDynamics->update[i];
			molecularDynamics->update[i] = FALSE;
		}
	}while(!done && nIter<maxIter);
	if(nIter>=maxIter && deltaMax>tolerance*10)
	{
		printf(("Rattle first portion : Warning, distance constraints not satisfied\n"));
		/*
		for (i = 0; i < molecularDynamics->numberOfAtoms; i++)
		{
			printf("atom#%d\n",i);
			printf("Old coord\n");
			for (k=0;k<3;k++) printf("%f ",molecularDynamics->coordinatesOld[i][k]);
			printf("\n");
			printf("New coord\n");
			for (k=0;k<3;k++) printf("%f ",m->atoms[i].coordinates[k]);
			printf("\n");
		}
		exit(1);
		*/
		for (i = 0; i < molecularDynamics->forceField->molecule.numberOfRattleConstraintsTerms; i++)
		{
			a1 = (int)molecularDynamics->forceField->molecule.rattleConstraintsTerms[0][i];
			a2 = (int)molecularDynamics->forceField->molecule.rattleConstraintsTerms[1][i];
			r2ij = 0;
			for (k=0;k<3;k++)
			{
				d = m->atoms[a2].coordinates[k]-m->atoms[a1].coordinates[k];
				r2ij +=d*d;
			}
			delta = molecularDynamics->forceField->molecule.rattleConstraintsTerms[2][i]-r2ij;
			printf("%d %d %s %s r2ij=%f r2Old=%f delta=%f\n",
			a1,a2,
			molecularDynamics->forceField->molecule.atoms[a1].mmType,
			molecularDynamics->forceField->molecule.atoms[a2].mmType,
			r2ij, molecularDynamics->forceField->molecule.rattleConstraintsTerms[2][i],delta);
		}
	}
	for (i = 0; i < molecularDynamics->numberOfAtoms; i++)
	{
		if(!molecularDynamics->forceField->molecule.atoms[i].variable)
		{
			for (k=0;k<3;k++) molecularDynamics->forceField->molecule.atoms[i].velocity[k] = 0.0;
			for (k=0;k<3;k++) m->atoms[i].coordinates[k] =  molecularDynamics->coordinatesOld[i][k];
		}
	}
#endif /* ENABLE_CL*/

}
/*********************************************************************************/
static void applyRattleSecondPortion(MolecularDynamics* molecularDynamics)
{
#ifdef ENABLE_CL
	ForceField* forceField = molecularDynamics->forceField;
	CLProp clProp=getCLProp();
	cl_int done=1;
	int nIter = 0;
	cl_int err;
	int maxIter = 1000;
	size_t global = forceField->molecule.nAtoms*forceField->nBlockRattleBuffer;
	err = clEnqueueNDRangeKernel(clProp.command_queue, molecularDynamics->initRattles, 1, NULL, &global, NULL, 0, NULL, NULL);
	if(err != CL_SUCCESS)
	{
		printErrorCLRun(err);
       		printf("I cannot execute initRattles,\n");
		exit(1);
	}
	clFinish(clProp.command_queue);
#ifdef DEBUG
	printf("applyRattleSecondPortion : End initRattles\n");
#endif
	do{
		nIter++;
		//printf("applyRattleSecondPortion : nIter =%d\n",nIter);
		global =  molecularDynamics->forceField->numberOfRattleConstraintsTerms;
		err = clEnqueueNDRangeKernel(clProp.command_queue, molecularDynamics->applyRattle2, 1, NULL, &global, NULL, 0, NULL, NULL);
		global = forceField->molecule.nAtoms;
		err = clEnqueueNDRangeKernel(clProp.command_queue, molecularDynamics->updateRattle2, 1, NULL, &global, NULL, 0, NULL, NULL);
		global = 1;
		err = clEnqueueNDRangeKernel(clProp.command_queue, molecularDynamics->testDoneRattle, 1, NULL, &global, NULL, 0, NULL, NULL);
		clEnqueueReadBuffer(clProp.command_queue, forceField->rattleDoneBufferCL, CL_TRUE, 0, sizeof(cl_int)*1, &done, 0, NULL, NULL);
	}while(done<1 && nIter<maxIter);
	if(nIter>=maxIter)
	{
		printf(("Rattle second portion : Warning, velocity constraints not satisfied\n"));
		exit(1);
	}
#else
	int i;
	int k;
	int maxIter = 1000;
	double omega = 1.0;
	double tolerance = 1e-6; 
	boolean done = FALSE;
	int nIter = 0;
	int a1 = 0;
	int a2 = 0;
	double r2ij;
	double dot;
	double invMass1;
	double invMass2;
	double term = 0;
	double terms[3];
	double d;
	Molecule* m = &molecularDynamics->forceField->molecule;
	ForceField* forceField = molecularDynamics->forceField;
	double deltaMax = 0;

	if(forceField->molecule.constraints==NOCONSTRAINTS) return;
	tolerance /= molecularDynamics->dt;
	for (i = 0; i < molecularDynamics->numberOfAtoms; i++)
	{
			molecularDynamics->moved[i] = molecularDynamics->forceField->molecule.atoms[i].variable;
			//molecularDynamics->moved[i] = TRUE;
			molecularDynamics->update[i] = FALSE;
	}
	/* maxIter *= molecularDynamics->forceField->numberOfRattleConstraintsTerms;*/
	do{
		nIter++;
		done=TRUE;
		deltaMax = 0;
		for (i = 0; i < molecularDynamics->forceField->molecule.numberOfRattleConstraintsTerms; i++)
		{
			a1 = (int)molecularDynamics->forceField->molecule.rattleConstraintsTerms[0][i];
			a2 = (int)molecularDynamics->forceField->molecule.rattleConstraintsTerms[1][i];
			r2ij = molecularDynamics->forceField->molecule.rattleConstraintsTerms[2][i];
			if( !molecularDynamics->moved[a1] && !molecularDynamics->moved[a2] ) continue;
			/* here : rattle image for PBC, not yet implemented */
			dot = 0;
			for (k=0;k<3;k++)
			{
				d = m->atoms[a2].coordinates[k]-m->atoms[a1].coordinates[k];
				dot +=d*(molecularDynamics->forceField->molecule.atoms[a2].velocity[k]-molecularDynamics->forceField->molecule.atoms[a1].velocity[k]);
			}
			invMass1 = 1/molecularDynamics->forceField->molecule.atoms[a1].mass;
			invMass2 = 1/molecularDynamics->forceField->molecule.atoms[a2].mass;
		        term = -dot / ((invMass1+invMass2)*r2ij);
			if(fabs(term)<=tolerance) continue;
			/* if(fabs(dot/r2ij)<=tolerance) continue;*/
			if(deltaMax<fabs(term)) deltaMax = fabs(term);

			done = FALSE;
			molecularDynamics->update[a1] = TRUE;
			molecularDynamics->update[a2] = TRUE;
		        term *= omega;

			for (k=0;k<3;k++)
			{
				d = m->atoms[a2].coordinates[k]-m->atoms[a1].coordinates[k];
				terms[k] = d*term;
			}
			for (k=0;k<3;k++) molecularDynamics->forceField->molecule.atoms[a1].velocity[k] -= terms[k]*invMass1;
			for (k=0;k<3;k++) molecularDynamics->forceField->molecule.atoms[a2].velocity[k] += terms[k]*invMass2;
		}
		for (i = 0; i < molecularDynamics->numberOfAtoms; i++)
		{
			molecularDynamics->moved[i] = molecularDynamics->update[i];
			molecularDynamics->update[i] = FALSE;
		}
	}while(!done && nIter<maxIter);
	if(nIter>=maxIter && deltaMax>tolerance*10)
	{
		printf(("Rattle second portion : Warning, velocity constraints not satisfied\n"));
		for (i = 0; i < molecularDynamics->forceField->numberOfRattleConstraintsTerms; i++)
		{
			a1 = (int)molecularDynamics->forceField->molecule.rattleConstraintsTerms[0][i];
			a2 = (int)molecularDynamics->forceField->molecule.rattleConstraintsTerms[1][i];
			r2ij = 0;
			for (k=0;k<3;k++)
			{
				d = m->atoms[a2].coordinates[k]-m->atoms[a1].coordinates[k];
				r2ij +=d*d;
			}
			dot = 0;
			for (k=0;k<3;k++)
			{
				d = m->atoms[a2].coordinates[k]-m->atoms[a1].coordinates[k];
				dot +=d*(molecularDynamics->forceField->molecule.atoms[a2].velocity[k]-molecularDynamics->forceField->molecule.atoms[a1].velocity[k]);
			}
			invMass1 = 1/molecularDynamics->forceField->molecule.atoms[a1].mass;
			invMass2 = 1/molecularDynamics->forceField->molecule.atoms[a2].mass;
		        term = -dot / ((invMass1+invMass2)*r2ij);
			printf("%d %d %s %s r2ij=%f r2Old=%f term=%f\n",
			a1,a2,
			molecularDynamics->forceField->molecule.atoms[a1].mmType,
			molecularDynamics->forceField->molecule.atoms[a2].mmType,
			r2ij, molecularDynamics->forceField->molecule.rattleConstraintsTerms[2][i],term);
		}
	}
	for (i = 0; i < molecularDynamics->numberOfAtoms; i++)
			if(!molecularDynamics->forceField->molecule.atoms[i].variable)
			for (k=0;k<3;k++) molecularDynamics->forceField->molecule.atoms[i].velocity[k] = 0.0;
					
#endif /* ENABLE_CL*/
}
/*********************************************************************************/
static void applyVerlet(MolecularDynamics* molecularDynamics)
{
#ifdef ENABLE_CL
	size_t global = molecularDynamics->forceField->molecule.nAtoms;
	cl_int err;
	CLProp clProp=getCLProp();

/*
	{
		size_t global = molecularDynamics->forceField->molecule.nAtoms;
		CLProp clProp = getCLProp();
		cl_int err = clEnqueueNDRangeKernel(clProp.command_queue, molecularDynamics->initAccelarations, 1, NULL, &global, NULL, 0, NULL, NULL);
		if(err != CL_SUCCESS)
		{
			printErrorCLRun(err);
        		printf("I cannot execute computeAccelarationsFromGradients\n");
			exit(1);
		}
		clFinish(clProp.command_queue);
	}
*/
	//newAccelaration(molecularDynamics);
	if(molecularDynamics->forceField->molecule.constraints!=NOCONSTRAINTS)
	{
		size_t global = molecularDynamics->forceField->molecule.nAtoms;
		cl_int err = clEnqueueNDRangeKernel(clProp.command_queue, molecularDynamics->copyPositions, 1, NULL, &global, NULL, 0, NULL, NULL);
		if(err != CL_SUCCESS)
		{
			printErrorCLRun(err);
        		printf("I cannot execute copyPositions\n");
			exit(1);
		}
		clFinish(clProp.command_queue);
	}
	err = clEnqueueNDRangeKernel(clProp.command_queue, molecularDynamics->applyVerlet1, 1, NULL, &global, NULL, 0, NULL, NULL);
	if(err != CL_SUCCESS)
	{
		printErrorCLRun(err);
        	printf("I cannot execute  applyVerlet1\n");
		exit(1);
	}
	clFinish(clProp.command_queue);
	if(molecularDynamics->forceField->molecule.constraints!=NOCONSTRAINTS) applyRattleFirstPortion(molecularDynamics);
	if(molecularDynamics->thermostat == NOSEHOOVER) nose_hoover(molecularDynamics);
	newAccelaration(molecularDynamics);
	err = clEnqueueNDRangeKernel(clProp.command_queue, molecularDynamics->applyVerlet2, 1, NULL, &global, NULL, 0, NULL, NULL);
	if(err != CL_SUCCESS)
	{
		printErrorCLRun(err);
        	printf("I cannot execute  applyVerlet2\n");
		exit(1);
	}
	clFinish(clProp.command_queue);
	if(molecularDynamics->forceField->molecule.constraints!=NOCONSTRAINTS) applyRattleSecondPortion(molecularDynamics);
	if(molecularDynamics->thermostat == NOSEHOOVER) nose_hoover(molecularDynamics);
#else
	int i;
	int j;

	if(molecularDynamics->forceField->molecule.constraints!=NOCONSTRAINTS)
	for (i = 0; i < molecularDynamics->numberOfAtoms; i++)
		for ( j = 0; j < 3; j++)
				molecularDynamics->coordinatesOld[i][j]= molecularDynamics->forceField->molecule.atoms[i].coordinates[j];

	for (i = 0; i < molecularDynamics->numberOfAtoms; i++)
	{
		if(!molecularDynamics->forceField->molecule.atoms[i].variable) continue;
		for ( j = 0; j < 3; j++)
		{
			molecularDynamics->forceField->molecule.atoms[i].coordinates[j] += 
				molecularDynamics->forceField->molecule.atoms[i].velocity[j] * molecularDynamics->dt +
				molecularDynamics->a[i][j]*molecularDynamics->dt2_2;	
		}
		for ( j = 0; j < 3; j++)
			molecularDynamics->forceField->molecule.atoms[i].velocity[j] += molecularDynamics->a[i][j] * molecularDynamics->dt_2;
	}

	if(molecularDynamics->forceField->molecule.constraints!=NOCONSTRAINTS) applyRattleFirstPortion(molecularDynamics);
	if(molecularDynamics->thermostat == NOSEHOOVER) nose_hoover(molecularDynamics);

	newAccelaration(molecularDynamics);

	for (i = 0; i < molecularDynamics->numberOfAtoms; i++)
	{
		if(!molecularDynamics->forceField->molecule.atoms[i].variable) continue;
		for ( j = 0; j < 3; j++)
			molecularDynamics->forceField->molecule.atoms[i].velocity[j] += molecularDynamics->a[i][j] * molecularDynamics->dt_2;
	}
	if(molecularDynamics->forceField->molecule.constraints!=NOCONSTRAINTS) applyRattleSecondPortion(molecularDynamics);
	if(molecularDynamics->thermostat == NOSEHOOVER) nose_hoover(molecularDynamics);
#endif
}
/*********************************************************************************/
static void applyBeeman(MolecularDynamics* molecularDynamics)
{
#ifdef ENABLE_CL
	size_t global = molecularDynamics->forceField->molecule.nAtoms;
	cl_int err;
	CLProp clProp=getCLProp();
	if(molecularDynamics->forceField->molecule.constraints!=NOCONSTRAINTS)
	{
		size_t global = molecularDynamics->forceField->molecule.nAtoms;
		cl_int err = clEnqueueNDRangeKernel(clProp.command_queue, molecularDynamics->copyPositions, 1, NULL, &global, NULL, 0, NULL, NULL);
		if(err != CL_SUCCESS)
		{
			printErrorCLRun(err);
        		printf("I cannot execute copyPositions\n");
			exit(1);
		}
		clFinish(clProp.command_queue);
	}
	err = clEnqueueNDRangeKernel(clProp.command_queue, molecularDynamics->applyBeeman1, 1, NULL, &global, NULL, 0, NULL, NULL);
	if(err != CL_SUCCESS)
	{
		printErrorCLRun(err);
        	printf("I cannot execute applyBeeman1\n");
		exit(1);
	}
	clFinish(clProp.command_queue);
	if(molecularDynamics->forceField->molecule.constraints!=NOCONSTRAINTS) applyRattleFirstPortion(molecularDynamics);
	newAccelaration(molecularDynamics);
	err = clEnqueueNDRangeKernel(clProp.command_queue, molecularDynamics->applyBeeman2, 1, NULL, &global, NULL, 0, NULL, NULL);
	if(err != CL_SUCCESS)
	{
		printErrorCLRun(err);
        	printf("I cannot execute applyBeeman2\n");
		exit(1);
	}
	clFinish(clProp.command_queue);
	if(molecularDynamics->forceField->molecule.constraints!=NOCONSTRAINTS) applyRattleSecondPortion(molecularDynamics);
#else
	int i;
	int j;
	double terms[3];

	if(molecularDynamics->forceField->molecule.constraints!=NOCONSTRAINTS)
	for (i = 0; i < molecularDynamics->numberOfAtoms; i++)
		for ( j = 0; j < 3; j++)
				molecularDynamics->coordinatesOld[i][j]= molecularDynamics->forceField->molecule.atoms[i].coordinates[j];
	for (i = 0; i < molecularDynamics->numberOfAtoms; i++)
	{
		if(!molecularDynamics->forceField->molecule.atoms[i].variable) continue;
		for ( j = 0; j < 3; j++)
			terms[j] = 5.0*molecularDynamics->a[i][j]-molecularDynamics->aold[i][j];

		for ( j = 0; j < 3; j++)
		{
			molecularDynamics->forceField->molecule.atoms[i].coordinates[j] += 
				molecularDynamics->forceField->molecule.atoms[i].velocity[j] * molecularDynamics->dt +
				terms[j]*molecularDynamics->dt2_8;	
		}
		for ( j = 0; j < 3; j++)
			molecularDynamics->forceField->molecule.atoms[i].velocity[j] += terms[j] * molecularDynamics->dt_8;
	}

	if(molecularDynamics->forceField->molecule.constraints!=NOCONSTRAINTS) applyRattleFirstPortion(molecularDynamics);

	newAccelaration(molecularDynamics);

	for (i = 0; i < molecularDynamics->numberOfAtoms; i++)
	{
		if(!molecularDynamics->forceField->molecule.atoms[i].variable) continue;
		for ( j = 0; j < 3; j++)
			molecularDynamics->forceField->molecule.atoms[i].velocity[j] += (3.0*molecularDynamics->a[i][j]+molecularDynamics->aold[i][j]) * molecularDynamics->dt_8;
	}

	if(molecularDynamics->forceField->molecule.constraints!=NOCONSTRAINTS) applyRattleSecondPortion(molecularDynamics);
#endif
}
/*********************************************************************************/
static void applyMartynaTuckerman(MolecularDynamics* molecularDynamics)
{
#ifdef ENABLE_CL
	printf("Sorry, Martyna&Tuckerman  integrator not yet implemented in GPU\n");
	exit(1);
#else
	int i;
	int j;
	Constraints constraints = molecularDynamics->forceField->molecule.constraints;
	Atom* atoms = molecularDynamics->forceField->molecule.atoms;
	
	//printf("Begin applyMartynaTuckerman\n");

	if(constraints!=NOCONSTRAINTS)
	for (i = 0; i < molecularDynamics->numberOfAtoms; i++)
		for ( j = 0; j < 3; j++)
				molecularDynamics->coordinatesOld[i][j]= atoms[i].coordinates[j];

	for (i = 0; i < molecularDynamics->numberOfAtoms; i++)
	{
		if(!atoms[i].variable) continue;
		for ( j = 0; j < 3; j++)
		{
			atoms[i].coordinates[j] += atoms[i].velocity[j] * molecularDynamics->dt +
					(5*molecularDynamics->a[i][j]-molecularDynamics->Ftilde[j][i])*molecularDynamics->dt2_8 +
					molecularDynamics->VF[j][i]*molecularDynamics->dt3_6 +
					molecularDynamics->VJ[j][i]*molecularDynamics->dt4_48;
		}
		for ( j = 0; j < 3; j++)
			atoms[i].velocity[j] += (4*molecularDynamics->a[i][j]+molecularDynamics->Ftilde[j][i]) * molecularDynamics->dt_8+
						molecularDynamics->VF[j][i]*molecularDynamics->dt2_8;
	}
	for (i = 0; i < molecularDynamics->numberOfAtoms; i++)
	for ( j = 0; j < 3; j++)
	{
		double f =  molecularDynamics->Ftilde[j][i];
		molecularDynamics->Ftilde[j][i] = 0.5*(2*molecularDynamics->a[i][j]-molecularDynamics->Ftilde[j][i]) + molecularDynamics->VF[j][i]*molecularDynamics->dt_2;
		molecularDynamics->VJ[j][i] = -molecularDynamics->VJ[j][i]-3.0/molecularDynamics->dt* molecularDynamics->VF[j][i]-1.5/molecularDynamics->dt2*f;
		molecularDynamics->VF[j][i] = -0.5*molecularDynamics->VF[j][i]-1.5/molecularDynamics->dt*f;
	}

	if(constraints!=NOCONSTRAINTS) applyRattleFirstPortion(molecularDynamics);
	if(molecularDynamics->thermostat == NOSEHOOVER) nose_hoover(molecularDynamics);

	//printf("End half applyMartynaTuckerman\n");

	newAccelaration(molecularDynamics);

	for (i = 0; i < molecularDynamics->numberOfAtoms; i++)
	{
		if(!atoms[i].variable) continue;
		for ( j = 0; j < 3; j++)
			atoms[i].velocity[j] += 3*molecularDynamics->a[i][j] * molecularDynamics->dt_8;
	}
	for (i = 0; i < molecularDynamics->numberOfAtoms; i++)
	for ( j = 0; j < 3; j++)
	{
		molecularDynamics->Ftilde[j][i] += 0.5*molecularDynamics->a[i][j];
		molecularDynamics->VF[j][i] += 1.5*molecularDynamics->a[i][j]/molecularDynamics->dt;
		molecularDynamics->VJ[j][i] += 1.5*molecularDynamics->a[i][j]/molecularDynamics->dt2;
	}

	if(constraints!=NOCONSTRAINTS) applyRattleSecondPortion(molecularDynamics);
	if(molecularDynamics->thermostat == NOSEHOOVER) nose_hoover(molecularDynamics);
	//printf("End applyMartynaTuckerman\n");
#endif
}
/*********************************************************************************/
static void newProperties(MolecularDynamics* molecularDynamics, char* comments)
{
	if( molecularDynamics->fileProp == NULL) return;
	fprintf(molecularDynamics->fileProp,"Time0(fs)\tTime(fs)\tTotal Energy(Kcal/mol)\tPotential Energy(kcal/mol)\tKinetic Energy(Kcal/mol)\tT(t) (K)\tTaver(K)\tsigma(T)(K)\tIndex");
	if(molecularDynamics->thermostat==NOSEHOOVER) fprintf(molecularDynamics->fileProp,"\tEtot+Etherm");
	if(comments) fprintf(molecularDynamics->fileProp,"%s\n", comments);
	else fprintf(molecularDynamics->fileProp,"\n");
}
/*********************************************************************************/
static void saveProperties(MolecularDynamics* molecularDynamics, int iStep0, int iStep, char* comments)
{
	double dt = molecularDynamics->dt/(fsInAKMA);
	static double Ttot = 0;
	static double T2tot = 0;
	double Taver = 0;
	double T2aver = 0;

	double totalEnergy =  molecularDynamics->totalEnergy;

	if( molecularDynamics->fileProp == NULL) return;
	if(molecularDynamics->thermostat==NOSEHOOVER)
	{
		int i;
		double kT = Kb* molecularDynamics->temperature;
		double e = molecularDynamics->vNH[0]*molecularDynamics->vNH[0]* molecularDynamics->qNH[0]/2 + (molecularDynamics->forceField->molecule.nFree)*kT* molecularDynamics->xNH[0];
		for(i=1;i<MAXNH;i++) e += molecularDynamics->vNH[i]*molecularDynamics->vNH[i]* molecularDynamics->qNH[i]/2 + kT* molecularDynamics->xNH[i];
		
		totalEnergy += e;
	}
	

	if(iStep==1)
	{
			Ttot = 0;
			T2tot = 0;
	}
	Ttot += molecularDynamics->kelvin;
	T2tot += molecularDynamics->kelvin*molecularDynamics->kelvin;
	Taver = Ttot/iStep;
	T2aver = T2tot/iStep;


	fprintf(molecularDynamics->fileProp,"%f\t%f\t%f\t\t%f\t\t\t%f\t\t\t%f\t%f\t%f\t%d\t", 
			(iStep0)*dt, 
			(iStep)*dt, 
			molecularDynamics->totalEnergy,
			molecularDynamics->potentialEnergy,
			molecularDynamics->kineticEnergy,
			molecularDynamics->kelvin,
			Taver,
			sqrt(fabs(T2aver-Taver*Taver)),
			molecularDynamics->index
			 );
	if( molecularDynamics->thermostat==NOSEHOOVER) fprintf(molecularDynamics->fileProp,"%f\t",totalEnergy);
	if(comments) fprintf(molecularDynamics->fileProp,"%s\n", comments);
	else fprintf(molecularDynamics->fileProp,"\n");
}
/*********************************************************************************/
static void saveTrajectory(MolecularDynamics* molecularDynamics, int iStep)
{
	double dt = molecularDynamics->dt/(fsInAKMA);
	int i;
	if( molecularDynamics->fileTraj == NULL) return;
	// Get geometry from Devices

	getGeometryVelocitiesCL( molecularDynamics->forceField,  &molecularDynamics->forceField->molecule);

	fprintf(molecularDynamics->fileTraj," %d %f %f %f %f nAtoms time(fs) TotalEnery(Kcal/mol) Kinetic Potential\n", 
			molecularDynamics->numberOfAtoms,
			 (iStep)*dt, 
			molecularDynamics->totalEnergy,
			molecularDynamics->kineticEnergy,
			molecularDynamics->potentialEnergy
			 );
	fprintf(molecularDynamics->fileTraj," %s\n", "Coord in Ang, Velocity in AKMA, time in fs");

	for (i = 0; i < molecularDynamics->numberOfAtoms; i++)
	{
		fprintf(molecularDynamics->fileTraj," %s %f %f %f %f %f %f %f %s %s %s %d %d\n", 
				molecularDynamics->forceField->molecule.atoms[i].prop.symbol,
				molecularDynamics->forceField->molecule.atoms[i].coordinates[0],
				molecularDynamics->forceField->molecule.atoms[i].coordinates[1],
				molecularDynamics->forceField->molecule.atoms[i].coordinates[2],
				molecularDynamics->forceField->molecule.atoms[i].velocity[0],
				molecularDynamics->forceField->molecule.atoms[i].velocity[1],
				molecularDynamics->forceField->molecule.atoms[i].velocity[2],
				molecularDynamics->forceField->molecule.atoms[i].charge,
				molecularDynamics->forceField->molecule.atoms[i].mmType,
				molecularDynamics->forceField->molecule.atoms[i].pdbType,
				molecularDynamics->forceField->molecule.atoms[i].residueName,
				molecularDynamics->forceField->molecule.atoms[i].residueNumber,
				molecularDynamics->forceField->molecule.atoms[i].variable
				);
	}
}

/**********************************************************************/
void	freeMolecularDynamics(MolecularDynamics* molecularDynamics)
{

	int i;
	if(molecularDynamics->a)
	{
		int i;
		for(i=0;i<molecularDynamics->numberOfAtoms;i++)
			if(molecularDynamics->a[i]) free(molecularDynamics->a[i]);
		free(molecularDynamics->a);
	}
	if(molecularDynamics->aold)
	{
		int i;
		for(i=0;i<molecularDynamics->numberOfAtoms;i++)
			if(molecularDynamics->aold[i]) free(molecularDynamics->aold[i]);
		free(molecularDynamics->aold);
	}
	if(molecularDynamics->coordinatesOld)
	{
		int i;
		for(i=0;i<molecularDynamics->numberOfAtoms;i++)
			if(molecularDynamics->coordinatesOld[i]) free(molecularDynamics->coordinatesOld[i]);
		free(molecularDynamics->coordinatesOld);
	}
	if(molecularDynamics->moved) free(molecularDynamics->moved);
	if(molecularDynamics->update) free(molecularDynamics->update);
	if(molecularDynamics->positionFriction) free(molecularDynamics->positionFriction);
	if(molecularDynamics->velocityFriction) free(molecularDynamics->velocityFriction);
	if(molecularDynamics->accelarationFriction) free(molecularDynamics->accelarationFriction);
	if(molecularDynamics->gamma) free(molecularDynamics->gamma);
	if(molecularDynamics->positionRandom) free(molecularDynamics->positionRandom);
	if(molecularDynamics->velocityRandom) free(molecularDynamics->velocityRandom);
	if(molecularDynamics->Ht) free(molecularDynamics->Ht);
	if(molecularDynamics->theta) free(molecularDynamics->theta);
	if(molecularDynamics->rnoise)
	{
		for(i=0;i<3*molecularDynamics->numberOfAtoms;i++)
			if(molecularDynamics->rnoise[i]) free(molecularDynamics->rnoise[i]);
		free(molecularDynamics->rnoise);
	}
	for(i=0;i<3;i++) if(molecularDynamics->Ftilde[i]) free(molecularDynamics->Ftilde[i]);
	for(i=0;i<3;i++) if(molecularDynamics->VF[i]) free(molecularDynamics->VF[i]);
	for(i=0;i<3;i++) if(molecularDynamics->VJ[i]) free(molecularDynamics->VJ[i]);
	molecularDynamics->forceField = NULL;
	molecularDynamics->numberOfAtoms = 0;
	molecularDynamics->updateFrequency = 0;
}
/********************************************************************************/
static double getEKin(MolecularDynamics* molecularDynamics)
{
#ifdef ENABLE_CL
	static cl_float ekin;
	CLProp clProp=getCLProp();

	size_t global = 1;
	cl_int err = clEnqueueNDRangeKernel(clProp.command_queue, molecularDynamics->computeKineticEnergy, 1, NULL, &global, NULL, 0, NULL, NULL);
	if(err != CL_SUCCESS)
	{
		printErrorCLRun(err);
       		printf("I cannot execute computeKineticEnergies\n");
		exit(1);
	}
	clFinish(clProp.command_queue);
	clEnqueueReadBuffer(clProp.command_queue, molecularDynamics->forceField->energyBufferCL, CL_TRUE, 0, sizeof(cl_float)*1, &ekin, 0, NULL, NULL);

/*
	int i,j;
	double mass;
	ForceField* forceField = molecularDynamics->forceField;
	clEnqueueReadBuffer(clProp.command_queue, forceField->atomsCL, CL_TRUE, 0, sizeof(cl_float8)*forceField->molecule.nAtoms, forceField->atomsCPU, 0, NULL, NULL);
	ekin = 0;
	for( i=0; i<forceField->molecule.nAtoms;i++)
	{
		mass = molecularDynamics->forceField->molecule.atoms[i].mass;
		for ( j = 0; j < 3; j++)
			ekin += forceField->atomsCPU[i].s[j+5]*forceField->atomsCPU[i].s[j+5]*mass;
	}
	ekin /=2;
*/
#else 
	double ekin = 0;
	int i;
	int j;
	double mass;
	for ( i = 0; i < molecularDynamics->numberOfAtoms; i++)
	{
		mass = molecularDynamics->forceField->molecule.atoms[i].mass;
		for ( j = 0; j < 3; j++)
			ekin += molecularDynamics->forceField->molecule.atoms[i].velocity[j]*molecularDynamics->forceField->molecule.atoms[i].velocity[j]*
				mass;
	}
	ekin /=2;
#endif
	return ekin;
}
/********************************************************************************/
static double getKelvin(MolecularDynamics* molecularDynamics)
{
	int nfree = molecularDynamics->forceField->molecule.nFree;
	if(nfree<1) return 0;
	return 2*getEKin(molecularDynamics) / ( nfree * Kb);
}
/********************************************************************************/
/*
     literature references:

     M. P. Allen, "Brownian Dynamics Simulation of a Chemical
     Reaction in Solution", Molecular Physics, 40, 1073-1087 (1980)

     F. Guarnieri and W. C. Still, "A Rapidly Convergent Simulation
     Method: Mixed Monte Carlo / Stochastic Dynamics", Journal of
     Computational Chemistry, 15, 1302-1310 (1994)
*/
/*********************************************************************************/
#ifndef ENABLE_CL
static void getFrictionalAndRandomForce(MolecularDynamics* molecularDynamics)
{
	double* gamma = molecularDynamics->gamma;
	double* positionFriction = molecularDynamics->positionFriction;
	double* velocityFriction = molecularDynamics->velocityFriction;
	double* accelarationFriction = molecularDynamics->accelarationFriction;
	double** positionRandom = molecularDynamics->positionRandom;
	double** velocityRandom = molecularDynamics->velocityRandom;
	double dt = molecularDynamics->dt;
	
	int n = molecularDynamics->numberOfAtoms;

	int i;
	int j;
	double gdt;
	double egdt;
	double ktm = 0;
	double pterm;
	double vterm;
        double psig;
        double vsig;
        double rho;
        double rhoc;
	double pnorm;
	double vnorm;

	for(i=0;i<n;i++)
        	gamma[i] = molecularDynamics->friction;

	/* printf(" friction = %f\n", molecularDynamics->friction);*/
	for(i=0;i<n;i++)
	{
		gdt = gamma[i] * dt;
		/* printf("gdt = %f\n",gdt);*/
		if (gdt <= 0.0)
		{
               		positionFriction[i] = 1.0;
			velocityFriction[i] = dt;
			accelarationFriction[i] = 0.5 * dt * dt;
			for(j=0;j<3;j++)
			{
                  		positionRandom[i][j] = 0.0;
                  		velocityRandom[i][j] = 0.0;
			}
		}
            	else
		{
			/* analytical expressions when friction coefficient is large */
               		if (gdt>=0.05)
			{
                  		egdt = exp(-gdt);
                  		positionFriction[i] = egdt;
                  		velocityFriction[i] = (1.0-egdt) / gamma[i];
                  		accelarationFriction[i] = (dt-velocityFriction[i]) / gamma[i];
                  		pterm = 2.0*gdt - 3.0 + (4.0-egdt)*egdt;
                  		vterm = 1.0 - egdt*egdt;
                  		rho = (1.0-egdt)*(1.0-egdt) / sqrt(pterm*vterm);
			}
			/* use series expansions when friction coefficient is small */
			else
			{
                  		double gdt2 = gdt * gdt;
                  		double gdt3 = gdt * gdt2;
                  		double gdt4 = gdt2 * gdt2;
                  		double gdt5 = gdt2 * gdt3;
                  		double gdt6 = gdt3 * gdt3;
                  		double gdt7 = gdt3 * gdt4;
                  		double gdt8 = gdt4 * gdt4;
                  		double gdt9 = gdt4 * gdt5;
                  		accelarationFriction[i] = (gdt2/2.0 - gdt3/6.0 + gdt4/24.0
                               	- gdt5/120.0 + gdt6/720.0
                               	- gdt7/5040.0 + gdt8/40320.0
                               	- gdt9/362880.0) / gamma[i]/gamma[i];
                  		velocityFriction[i] = dt - gamma[i]*accelarationFriction[i];
                  		positionFriction[i] = 1.0 - gamma[i]*velocityFriction[i];
                  		pterm = 2.0*gdt3/3.0 - gdt4/2.0
                            	+ 7.0*gdt5/30.0 - gdt6/12.0
                            	+ 31.0*gdt7/1260.0 - gdt8/160.0
                            	+ 127.0*gdt9/90720.0;
                  		vterm = 2.0*gdt - 2.0*gdt2 + 4.0*gdt3/3.0
                            	- 2.0*gdt4/3.0 + 4.0*gdt5/15.0
                            	- 4.0*gdt6/45.0 + 8.0*gdt7/315.0
                            	- 2.0*gdt8/315.0 + 4.0*gdt9/2835.0;
                  		rho = sqrt(3.0) * (0.5 - 3.0*gdt/16.0
                            	- 17.0*gdt2/1280.0
                            	+ 17.0*gdt3/6144.0
                            	+ 40967.0*gdt4/34406400.0
                            	- 57203.0*gdt5/275251200.0
                            	- 1429487.0*gdt6/13212057600.0);
			}
               		ktm = Kb * molecularDynamics->temperature / molecularDynamics->forceField->molecule.atoms[i].mass;
               		psig = sqrt(ktm*pterm) / gamma[i];
               		vsig = sqrt(ktm*vterm);
               		rhoc = sqrt(1.0 - rho*rho);
			for(j=0;j<3;j++)
			{
                		pnorm = normal();
             			vnorm = normal ();
				positionRandom[i][j] = psig * pnorm;
                  		velocityRandom[i][j] = vsig * (rho*pnorm+rhoc*vnorm);
			}
		}
	}
}
#endif
/*********************************************************************************/
static void applyStochastic(MolecularDynamics* molecularDynamics)
{
#ifdef ENABLE_CL
	size_t global = molecularDynamics->forceField->molecule.nAtoms;
	cl_int err;
	CLProp clProp=getCLProp();
	cl_float temp = molecularDynamics->temperature;
	resetRandomNumbers(molecularDynamics);
	clSetKernelArg(molecularDynamics->setFrictionalAndRandomForce, 5, sizeof(cl_float), &temp);
	printf("Call setFrictionalAndRandomForce\n");
	err = clEnqueueNDRangeKernel(clProp.command_queue, molecularDynamics->setFrictionalAndRandomForce, 1, NULL, &global, NULL, 0, NULL, NULL);
	if(err != CL_SUCCESS)
	{
		printErrorCLRun(err);
        	printf("I cannot execute setFrictionalAndRandomForce\n");
		exit(1);
	}
	clFinish(clProp.command_queue);
	printf("Call applyStochastic1\n");
	if(molecularDynamics->forceField->molecule.constraints!=NOCONSTRAINTS)
	{
		size_t global = molecularDynamics->forceField->molecule.nAtoms;
		cl_int err = clEnqueueNDRangeKernel(clProp.command_queue, molecularDynamics->copyPositions, 1, NULL, &global, NULL, 0, NULL, NULL);
		if(err != CL_SUCCESS)
		{
			printErrorCLRun(err);
        		printf("I cannot execute copyPositions\n");
			exit(1);
		}
		clFinish(clProp.command_queue);
	}
	err = clEnqueueNDRangeKernel(clProp.command_queue, molecularDynamics->applyStochastic1, 1, NULL, &global, NULL, 0, NULL, NULL);
	if(err != CL_SUCCESS)
	{
		printErrorCLRun(err);
        	printf("I cannot execute applyStochastic1\n");
		exit(1);
	}
	clFinish(clProp.command_queue);
	if(molecularDynamics->forceField->molecule.constraints!=NOCONSTRAINTS) applyRattleFirstPortion(molecularDynamics);
	printf("Call newAccelaration\n");
	newAccelaration(molecularDynamics);
	printf("Call applyStochastic2\n");

	err = clEnqueueNDRangeKernel(clProp.command_queue, molecularDynamics->applyStochastic2, 1, NULL, &global, NULL, 0, NULL, NULL);
	if(err != CL_SUCCESS)
	{
		printErrorCLRun(err);
        	printf("I cannot execute applyStochastic2\n");
		exit(1);
	}
	clFinish(clProp.command_queue);

	if(molecularDynamics->forceField->molecule.constraints!=NOCONSTRAINTS) applyRattleSecondPortion(molecularDynamics);
	// TO DO remove translation and rotation
	computeEnergies(molecularDynamics);
#else
	double* positionFriction = molecularDynamics->positionFriction;
	double* velocityFriction = molecularDynamics->velocityFriction;
	double* accelarationFriction = molecularDynamics->accelarationFriction;
	double** positionRandom = molecularDynamics->positionRandom;
	double** velocityRandom = molecularDynamics->velocityRandom;
	double**a = molecularDynamics->a;
	int n = molecularDynamics->numberOfAtoms;
	int i;
	int j;
	Atom* atoms = molecularDynamics->forceField->molecule.atoms;

	getFrictionalAndRandomForce(molecularDynamics);

	if(molecularDynamics->forceField->molecule.constraints!=NOCONSTRAINTS)
	for(i=0;i<n;i++)
		for(j=0;j<3;j++)
			molecularDynamics->coordinatesOld[i][j]= molecularDynamics->forceField->molecule.atoms[i].coordinates[j];
	for(i=0;i<n;i++)
	{
		if(!molecularDynamics->forceField->molecule.atoms[i].variable) continue;
		for(j=0;j<3;j++)
			atoms[i].coordinates[j] += molecularDynamics->forceField->molecule.atoms[i].velocity[j]*velocityFriction[i] + a[i][j]*accelarationFriction[i] + positionRandom[i][j];
		for(j=0;j<3;j++)
			molecularDynamics->forceField->molecule.atoms[i].velocity[j] = molecularDynamics->forceField->molecule.atoms[i].velocity[j]*positionFriction[i] + 0.5*a[i][j]*velocityFriction[i];
	}
	if(molecularDynamics->forceField->molecule.constraints!=NOCONSTRAINTS) applyRattleFirstPortion(molecularDynamics);
	newAccelaration(molecularDynamics);

	for (i = 0; i < n; i++)
	{
		if(!molecularDynamics->forceField->molecule.atoms[i].variable) continue;
		for ( j = 0; j < 3; j++)
			molecularDynamics->forceField->molecule.atoms[i].velocity[j] += 0.5*a[i][j]*velocityFriction[i] + velocityRandom[i][j];
	}
	if(molecularDynamics->forceField->molecule.constraints!=NOCONSTRAINTS) applyRattleSecondPortion(molecularDynamics);
	computeEnergies(molecularDynamics);

	removeTranslationAndRotation(molecularDynamics);
#endif
}
/*********************************************************************************/
static void applyQTB(MolecularDynamics* molecularDynamics)
{
	
	int n = molecularDynamics->numberOfAtoms;
	int i;
	int j;
	//double gp = 1.0/(1.0+molecularDynamics->friction*molecularDynamics->dt_2);
	//double gm = (1.0-molecularDynamics->friction*molecularDynamics->dt_2)*gp;
	double gf = (1.0-molecularDynamics->friction*molecularDynamics->dt_2);

#ifdef ENABLE_CL
	printf("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n");
	printf("Sorry, QTB dynamics is not yet implemented on GPU\n");
	printf("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n");
	exit(1);
#endif
	/* printf("gm = %f gp =%f\n",gm,gp);*/
	if(molecularDynamics->forceField->molecule.constraints!=NOCONSTRAINTS)
	for (i = 0; i < n; i++)
		for ( j = 0; j < 3; j++)
				molecularDynamics->coordinatesOld[i][j]= molecularDynamics->forceField->molecule.atoms[i].coordinates[j];

	for (i = 0; i < n; i++)
	if(molecularDynamics->forceField->molecule.atoms[i].variable)
	for ( j = 0; j < 3; j++)
	{
		molecularDynamics->forceField->molecule.atoms[i].coordinates[j] += 
		molecularDynamics->forceField->molecule.atoms[i].velocity[j]*molecularDynamics->dt +
		(molecularDynamics->a[i][j]+molecularDynamics->theta[3*i+j]-molecularDynamics->forceField->molecule.atoms[i].velocity[j]*molecularDynamics->friction)*molecularDynamics->dt2_2;	
	}
	for (i = 0; i < n; i++)
	if(molecularDynamics->forceField->molecule.atoms[i].variable)
	for ( j = 0; j < 3; j++)
			molecularDynamics->forceField->molecule.atoms[i].velocity[j] = 
			gf*molecularDynamics->forceField->molecule.atoms[i].velocity[j] +
			(molecularDynamics->a[i][j]+molecularDynamics->theta[3*i+j])*molecularDynamics->dt_2;
            
	if(molecularDynamics->forceField->molecule.constraints!=NOCONSTRAINTS) applyRattleFirstPortion(molecularDynamics);

	newAccelaration(molecularDynamics);

	for (i = 0; i < n; i++)
	if(molecularDynamics->forceField->molecule.atoms[i].variable)
	for ( j = 0; j < 3; j++)
			molecularDynamics->forceField->molecule.atoms[i].velocity[j] = 
			gf*molecularDynamics->forceField->molecule.atoms[i].velocity[j] +
			(molecularDynamics->a[i][j]+molecularDynamics->theta[3*i+j])*molecularDynamics->dt_2;

	if(molecularDynamics->forceField->molecule.constraints!=NOCONSTRAINTS) applyRattleSecondPortion(molecularDynamics);

	removeTranslationAndRotation(molecularDynamics);

	computeEnergies(molecularDynamics);
}
/*********************************************************************************/
static void applyLangevin(MolecularDynamics* molecularDynamics)
{
	
	int n = molecularDynamics->numberOfAtoms;
	int i;
	int j;
	//double gp = 1.0/(1.0+molecularDynamics->friction*molecularDynamics->dt_2);
	//double gm = (1.0-molecularDynamics->friction*molecularDynamics->dt_2)*gp;
	double gf = (1.0-molecularDynamics->friction*molecularDynamics->dt_2);

#ifdef ENABLE_CL
	printf("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n");
	printf("Sorry, Langevin dynamics is not yet implemented on GPU\n");
	printf("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n");
	exit(1);
#endif
	/* printf("gm = %f gp =%f\n",gm,gp);*/
	if(molecularDynamics->forceField->molecule.constraints!=NOCONSTRAINTS)
	for (i = 0; i < n; i++)
		for ( j = 0; j < 3; j++)
				molecularDynamics->coordinatesOld[i][j]= molecularDynamics->forceField->molecule.atoms[i].coordinates[j];

	for (i = 0; i < n; i++)
	if(molecularDynamics->forceField->molecule.atoms[i].variable)
	for ( j = 0; j < 3; j++)
	{
		molecularDynamics->forceField->molecule.atoms[i].coordinates[j] += 
		molecularDynamics->forceField->molecule.atoms[i].velocity[j]*molecularDynamics->dt +
		(molecularDynamics->a[i][j]+molecularDynamics->theta[3*i+j]-molecularDynamics->forceField->molecule.atoms[i].velocity[j]*molecularDynamics->friction)*molecularDynamics->dt2_2;	
	}
	for (i = 0; i < n; i++)
	if(molecularDynamics->forceField->molecule.atoms[i].variable)
	for ( j = 0; j < 3; j++)
			molecularDynamics->forceField->molecule.atoms[i].velocity[j] = 
			gf*molecularDynamics->forceField->molecule.atoms[i].velocity[j] +
			(molecularDynamics->a[i][j]+molecularDynamics->theta[3*i+j])*molecularDynamics->dt_2;
            
	if(molecularDynamics->forceField->molecule.constraints!=NOCONSTRAINTS) applyRattleFirstPortion(molecularDynamics);

	newAccelaration(molecularDynamics);

	for (i = 0; i < n; i++)
	if(molecularDynamics->forceField->molecule.atoms[i].variable)
	for ( j = 0; j < 3; j++)
			molecularDynamics->forceField->molecule.atoms[i].velocity[j] = 
			gf*molecularDynamics->forceField->molecule.atoms[i].velocity[j] +
			(molecularDynamics->a[i][j]+molecularDynamics->theta[3*i+j])*molecularDynamics->dt_2;

	if(molecularDynamics->forceField->molecule.constraints!=NOCONSTRAINTS) applyRattleSecondPortion(molecularDynamics);

	removeTranslationAndRotation(molecularDynamics);

	computeEnergies(molecularDynamics);
}
