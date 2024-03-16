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

/* QuantumMechanicsMD.c  */

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
#include "../Utils/Types.h"
#include "../Molecule/Molecule.h"
#include "../QuantumMechanics/QuantumMechanicsModel.h"
#include "../QuantumMechanics/QuantumMechanics.h"
#include "../QuantumMechanics/QuantumMechanicsMD.h"


/*********************************************************************************/
static void initMD(QuantumMechanicsMD* molecularDynamics, double temperature, double stepSize, MDIntegratorType integratorType, MDThermostatType thermostat, double friction, double omegaMax, int Nf, double collide, double qNH, char* fileNameTraj, char* fileNameProp, int numberOfRunSteps, int index);
static void berendsen(QuantumMechanicsMD* molecularDynamics);
static void scaleV(QuantumMechanicsMD* molecularDynamics);
static void bussi(QuantumMechanicsMD* molecularDynamics);
static void andersen(QuantumMechanicsMD* molecularDynamics);
static void nose_hoover(QuantumMechanicsMD* molecularDynamics);
static void rescaleVelocities(QuantumMechanicsMD* molecularDynamics);
static void computeEnergies(QuantumMechanicsMD* molecularDynamics);
static void applyOneStep(QuantumMechanicsMD* molecularDynamics, int iStep);
static void applyThermostat(QuantumMechanicsMD* molecularDynamics);
static void applyVerlet(QuantumMechanicsMD* molecularDynamics);
static void applyBeeman(QuantumMechanicsMD* molecularDynamics);
static void applyStochastic(QuantumMechanicsMD* molecularDynamics);
static void applyQTB(QuantumMechanicsMD* molecularDynamics);
static void updateQTB(QuantumMechanicsMD* molecularDynamics);
static void resetQTB(QuantumMechanicsMD* molecularDynamics);
static void applyLangevin(QuantumMechanicsMD* molecularDynamics);
static void updateLangevin(QuantumMechanicsMD* molecularDynamics);
static void resetLangevin(QuantumMechanicsMD* molecularDynamics);
static void newProperties(QuantumMechanicsMD* molecularDynamics, char* comments);
static void saveProperties(QuantumMechanicsMD* molecularDynamics, int iStep0, int iStep, char* comments);
static void saveTrajectory(QuantumMechanicsMD* molecularDynamics, int iStep);
static double getEKin(QuantumMechanicsMD* molecularDynamics);
static double getKelvin(QuantumMechanicsMD* molecularDynamics);
static void removeTranslationAndRotationTheta(QuantumMechanicsMD* molecularDynamics);
static void removeTranslationAndRotation(QuantumMechanicsMD* molecularDynamics);
/*****************************************************************************************************************/
// return # temperatures
static int getNumExchangeReplica(double* energies, double* temperatures, int nTemperatures, int* numMDForTemperatures)
{
	int it = rand()%(nTemperatures-1);
	int n = numMDForTemperatures[it];
	int next = numMDForTemperatures[it+1];
	double delta = (1.0/(Kb*temperatures[n])- 1.0/(Kb*temperatures[next]))*(energies[n]-energies[next]);
	boolean exchange= FALSE;
#ifdef DEBUG
	printf("n=%d next=%d delta = %f\n",n,next, delta);
#endif
	if(delta<0) exchange= TRUE;
	else
	{
		if(rand()/(double)RAND_MAX>exp(-delta)) exchange = TRUE;
	}
	if(exchange) return it;
	else return -1;
}
/*****************************************************************************************************************/
static void applyExchangeReplica(QuantumMechanicsMD* md, int nTemperatures, int n, int next)
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
			 md[n].qmModel->molecule.atoms[i].velocity[j] *= r;
			 md[next].qmModel->molecule.atoms[i].velocity[j] /= r;
		}
	}
	removeTranslationAndRotation(&md[n]);
	removeTranslationAndRotation(&md[next]);
	resetQTB(&md[n]);
	resetQTB(&md[next]);
	resetLangevin(&md[n]);
	resetLangevin(&md[next]);
}
/*****************************************************************************************************************/
static void changeOneReplica(QuantumMechanicsMD* md, int n, double newTemperature, char* fileNamePrefixProp,  char* fileNamePrefixTraj)
{
	double r;
	int i,j;
	char*  fileNameProp = NULL;
	char*  fileNameTraj = NULL;
	r = sqrt(newTemperature/md[n].temperature);
	/* change temperature */
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
			 md[n].qmModel->molecule.atoms[i].velocity[j] *= r;
	}
	removeTranslationAndRotation(&md[n]);
	if(fileNameProp) free(fileNameProp);
	if(fileNameTraj) free(fileNameTraj);
}
/*****************************************************************************************************************/
/*
static void exchangeReplicaLocalAll(QuantumMechanicsMD* md, int nTemperatures, int* numMDForTemperatures)
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
			 md[n].qmModel->molecule.atoms[i].velocity[j] *= r;
			 md[next].qmModel->molecule.atoms[i].velocity[j] /= r;
		}
	}
	removeTranslationAndRotation(&md[n]);
	removeTranslationAndRotation(&md[next]);
}
*/
/*****************************************************************************************************************/
static void exchangeReplica(QuantumMechanicsMD* md,  double* energiesAll, double* energies,  double* temperaturesAll, int* numMDForTemperatures, int* nTemperaturesLocal, int nTemperaturesAll, int nTemperatures, int nproc, int rank, int iBegin, char* fileNamePrefixProp, char* fileNamePrefixTraj)
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
				fprintf(md[0].qmModel->logfile, "Get energies from rank %d\n",j);
#endif
#ifdef ENABLE_MPI
				code = MPI_Recv(energies,nTemperaturesLocal[j],MPI_DOUBLE,j,tag,MPI_COMM_WORLD,&status) ;
#endif
#if DEBUG
				fprintf(md[0].qmModel->logfile, "End Get energies from rank %d\n",j);
#endif
				for(ii=0;ii<nTemperaturesLocal[j];ii++) energiesAll[k++] = energies[ii];
			}
			it = getNumExchangeReplica(energiesAll, temperaturesAll, nTemperaturesAll, numMDForTemperatures);
#if DEBUG
			fprintf(md[0].qmModel->logfile, "\n\nn=%d\n",it);
#endif
			if(it>=0)
			{
				n = numMDForTemperatures[it];
				next = numMDForTemperatures[it+1];
				numMDForTemperatures[it]=next;
				numMDForTemperatures[it+1]=n;
			}
			// send it to all others 
#if DEBUG
			fprintf(md[0].qmModel->logfile, "Send n from 0 to all others rank\n");
#endif
#ifdef ENABLE_MPI
			tag = 200;
			for(j=1;j<nproc;j++) code = MPI_Send(&it,1,MPI_INT,j,tag,MPI_COMM_WORLD) ;
#endif
#if DEBUG
			fprintf(md[0].qmModel->logfile, "End Send n from 0 to all others rank\n");
#endif
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
#ifdef ENABLE_MPI
			tag = 300;
#endif
			for(j=1;j<nproc;j++) 
			{
#if DEBUG
				fprintf(md[0].qmModel->logfile, "Get dum val for check close file at rank %d\n",j);
#endif
#ifdef ENABLE_MPI
				code = MPI_Recv(&dum,1,MPI_INT,j,tag,MPI_COMM_WORLD,&status) ;
#endif
#if DEBUG
				fprintf(md[0].qmModel->logfile, "End Get dum val for check close file at rank %d\n",j);
#endif
			}
			/// all files closed on other rank
#ifdef ENABLE_MPI
			tag = 400;
#endif
#if DEBUG
			fprintf(md[0].qmModel->logfile, "Send n from 0 to all others rank\n");
#endif
#ifdef ENABLE_MPI
			for(j=1;j<nproc;j++) code = MPI_Send(&it,1,MPI_INT,j,tag,MPI_COMM_WORLD) ;
#endif

			if(n>=0 && n<nTemperatures && next>=0 && next<nTemperatures)
			{
#if DEBUG
				fprintf(md[0].qmModel->logfile, " Exchange between 2 trajs on one proc rank = 0\n");
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
			fprintf(md[0].qmModel->logfile, "Send energies from %d to 0\n",rank);
#endif
			code = MPI_Send(energies,nTemperatures,MPI_DOUBLE,0,tag,MPI_COMM_WORLD) ;
#if DEBUG
			fprintf(md[0].qmModel->logfile, "End Send energies from %d to 0\n",rank);
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
				fprintf(md[0].qmModel->logfile, " Exchange between 2 trajs on one proc rank = 0\n");
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
/*****************************************************************************************************************/
QuantumMechanicsModel**    runQuantumMechanicsREMDConfo(
		QuantumMechanicsMD* molecularDynamics, QuantumMechanicsModel* qmModel, 
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
	QuantumMechanicsModel** geometries = NULL;
	int iSel = 0;
	int stepSel = 1;
	int stepExchange = 1;
	QuantumMechanicsMD* md = NULL;
	double* runTemps = NULL;
	double a = 1.0;
	double b = 1.0;
	int nTemperatures = nTemperaturesAll;
	int nproc = 1;
	int* nTemperaturesLocal = NULL;
	int rank = 0;
	int n = 0;
	int iBegin = 0;
	double* energies = NULL;
	double* energiesAll = NULL;
	double* temperaturesAll = NULL;
	int* numMDForTemperatures = NULL;
	FILE* logfile = qmModel->logfile;
	char* fileNamePrefixProp = NULL;
	char* fileNamePrefixTraj = NULL;

	if(nTemperaturesAll<1) return NULL;
	if(qmModel->molecule.nAtoms<1) return NULL;
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


	md = malloc(nTemperatures*sizeof(QuantumMechanicsMD));
	printf("end md alloc\n");
	fflush(logfile);
	for(k=0;k<nTemperatures;k++) 
	{
		md[k] = *molecularDynamics;
		md[k].qmModel = malloc(sizeof(QuantumMechanicsModel)); 
		*md[k].qmModel = qmModel->klass->copy(qmModel);
		md[k].numberOfAtoms = qmModel->molecule.nAtoms;
		fprintf(logfile, "nAtoms = %d\n",md[k].numberOfAtoms);
		md[k].updateFrequency = updateFrequency;
	}
	fprintf(logfile, "end md copy\n");
	fflush(logfile);
	updateNumber = malloc(nTemperatures*sizeof(int));
	for(k=0;k<nTemperatures;k++) updateNumber[k] = 0;

	fprintf(logfile,"\nTemperatures = ");
	for(k=0;k<nTemperatures;k++) fprintf(logfile,"%f ", runTemps[k]);
	fprintf(logfile,"\n\n");
	fflush(logfile);

	geometries = malloc(numberOfGeometries*sizeof(QuantumMechanicsModel*));
	for(i=0;i<numberOfGeometries;i++) geometries[i] = NULL;

	currentTemp = heatTemperature/2;
	numberOfHeatSteps = heatTime/stepSize*1000;
	numberOfEquiSteps = equiTime/stepSize*1000;; 
	numberOfRunSteps = runTime/stepSize*1000;; 

	fprintf(logfile,"begin init\n");
	fflush(logfile);
	for(k=0;k<nTemperatures;k++) 
	{
		char* fileNamePropk = NULL;
		char* fileNameTrajk = NULL;
		if(fileNameProp) fileNamePrefixProp = getSuffixNameFile(fileNameProp);
		if(fileNamePrefixProp) fileNamePropk = strdup_printf("%s%0.0f.txt",fileNamePrefixProp,runTemps[k]);
		if(fileNameTraj) fileNamePrefixTraj = getSuffixNameFile(fileNameTraj);
		if(fileNamePrefixTraj) fileNameTrajk = strdup_printf("%s%0.0f.gab",fileNamePrefixTraj,runTemps[k]);
		currentTemp = heatTemperature;
		if(numberOfHeatSteps==0) currentTemp = runTemps[k];
		initMD(&md[k],currentTemp,stepSize, integratorType, thermostat, friction, omegaMax, Nf, collide, qNH, fileNameTrajk, fileNamePropk, numberOfRunSteps,k+iBegin);
		md[k].qmModel->klass->calculateGradient(md[k].qmModel);
		computeEnergies(&md[k]);
	}
	fprintf(logfile,"end init\n");
	fflush(logfile);

	iSel = -1;
	if(rank==0)
	{
		if(str) free(str);
		str = strdup_printf("Geometry selected Potential energy =  %0.4f", md[0].potentialEnergy);
		/* redrawMolecule(&md[0].qmModel->molecule,str);*/
		fprintf(logfile,"%s\n",str);
		fflush(logfile);
		iSel++;
		geometries[iSel] = malloc(sizeof(QuantumMechanicsModel));
		*geometries[iSel] = md[0].qmModel->klass->copy(md[0].qmModel);
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
		//md[k].temperature = currentTemp;
		applyOneStep(&md[k],i);
		md[k].qmModel->firstRun = FALSE;
		currentTemp = heatTemperature + ( runTemps[k] - heatTemperature ) *
				( ( double )( i + 1 )/ numberOfHeatSteps );
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
			/* redrawMolecule(&md[k].qmModel->molecule,str);*/
			fprintf(logfile,"%s\n",str);
			fflush(logfile);
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
		md[k].qmModel->firstRun = FALSE;
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
			//redrawMolecule(&md[k].qmModel->molecule,str);
			fprintf(logfile, "%s\n",str);
			fflush(logfile);
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
	if(str) free(str);
	str = strdup_printf(("Geometry selected Potential energy =  %0.4f"), md[0].potentialEnergy);
	fprintf(logfile,"%s\n",str);
	fflush(logfile);
	if(numberOfGeometries>2) stepSel = numberOfRunSteps/(numberOfGeometries-1);
	else stepSel = numberOfRunSteps;
	if(numberOfExchanges>2) stepExchange =  numberOfRunSteps/(numberOfExchanges-1);
	else stepExchange = numberOfRunSteps;

	for(k=0;k<nTemperatures;k++) md[k].temperature = runTemps[k];
	for (i = 0; i < numberOfRunSteps; i++ )
	{
		for(k=0;k<nTemperatures;k++) 
		{
			applyOneStep(&md[k],i);
			md[k].qmModel->firstRun = FALSE;
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
				//redrawMolecule(&md[k].qmModel->molecule,str);
				fprintf(logfile,"%s\n",str);
				fflush(logfile);
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
			//redrawMolecule(&md[k].qmModel->molecule,str);
			fprintf(logfile,"%s\n",str);
			fflush(logfile);
			iSel++;
			geometries[iSel] = malloc(sizeof(QuantumMechanicsModel));
			*geometries[iSel] = md[k].qmModel->klass->copy(md[k].qmModel);
			/* waiting(0.1);*/
		}
/* Exchange here */
		if((i+1)%stepExchange==0)
		{
			exchangeReplica(md,  energiesAll, energies,  temperaturesAll, numMDForTemperatures, nTemperaturesLocal, nTemperaturesAll, nTemperatures, nproc, rank,  iBegin, fileNamePrefixProp, fileNamePrefixTraj);
		}
	}
	for(k=0;k<nTemperatures;k++) 
	if(fabs(temperaturesAll[0]-md[k].temperature)<1e-10 && iSel<numberOfGeometries-1)
	{
		if(str) free(str);
		str = strdup_printf(("Geometry selected Potential energy =  %0.4f"), md[k].potentialEnergy);
		//redrawMolecule(&md[0].qmModel->molecule,str);
		fprintf(logfile, "%s\n",str);
		fflush(logfile);
		iSel++;
		geometries[iSel] = malloc(sizeof(QuantumMechanicsModel));
		*geometries[iSel] = md[k].qmModel->klass->copy(md[k].qmModel);
		/* waiting(0.1);*/
	}

	n0 += numberOfRunSteps;

	fprintf(logfile,"End of MD Simulation on rank = %d\n",rank);
	fflush(logfile);
	for(k=0;k<nTemperatures;k++) 
	{
		md[k].qmModel->klass->calculateGradient(md[k].qmModel);
        	gradientNorm = 0;
		for (i = 0; i < md[k].numberOfAtoms; i++)
			for ( j = 0; j < 3; j++)
                        	gradientNorm += 
				md[k].qmModel->molecule.atoms[i].gradient[j] * 
				md[k].qmModel->molecule.atoms[i].gradient[j]; 

        	gradientNorm = sqrt( gradientNorm );
		if(str) free(str);
		str = strdup_printf(("T(K) = %0.2f Gradient = %f Ekin = %f (Kcal/mol) EPot =  %0.4f ETot =  %0.4f T(t) = %0.2f"),
			runTemps[k],
			(double)gradientNorm,
			md[k].kineticEnergy,
			md[k].potentialEnergy,
			md[k].totalEnergy,
			md[k].kelvin 
			); 
		//redrawMolecule(&md[0].qmModel->molecule,str);
		fprintf(logfile, "%s\n",str);
		fflush(logfile);
	}
	if(str) free(str);

	for(k=0;k<nTemperatures;k++) 
	{
		if(md[k].fileTraj)fclose(md[k].fileTraj);
		if(md[k].fileProp)fclose(md[k].fileProp);
		freeQuantumMechanicsMD(&md[k]);
	}
	free(updateNumber);
	free(runTemps);

	return geometries;
}
/**********************************************************************/
QuantumMechanicsModel**    runQuantumMechanicsMDConfo(
		QuantumMechanicsMD* molecularDynamics, QuantumMechanicsModel* qmModel, 
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
	QuantumMechanicsModel** geometries = NULL;
	int iSel = 0;
	int stepSel = 1;
	double e0  = 0;
	double esum  = 0;
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

	if(qmModel->molecule.nAtoms<1) return NULL;
	if(numberOfGeometries<2) return NULL;
	geometries = malloc(numberOfGeometries*sizeof(QuantumMechanicsModel*));
	for(i=0;i<numberOfGeometries;i++) geometries[i] = NULL;

	molecularDynamics->qmModel = qmModel;
	molecularDynamics->numberOfAtoms = qmModel->molecule.nAtoms;
	molecularDynamics->updateFrequency = updateFrequency;

	currentTemp = heatTemperature/2;
	
	numberOfHeatSteps = heatTime/stepSize*1000;
	numberOfEquiSteps = equiTime/stepSize*1000;; 
	numberOfRunSteps = runTime/stepSize*1000;; 


	currentTemp = heatTemperature;
	if(numberOfHeatSteps==0) currentTemp = equiTemperature; 
	if(numberOfHeatSteps==0 && numberOfEquiSteps==0 ) currentTemp = runTemperature; 

	initMD(molecularDynamics,currentTemp,stepSize,integratorType, thermostat, friction, omegaMax, Nf, collide, qNH, fileNameTraj, fileNameProp, numberOfRunSteps,0);
	molecularDynamics->qmModel->klass->calculateGradient(molecularDynamics->qmModel);
	computeEnergies(molecularDynamics);
	e0 = molecularDynamics->potentialEnergy;
	printf("E0 = The first potential energy in kcal = %f\n",e0); 

	iSel =-1;
	if((i+1)%stepSel==0 && (iSel+1)<numberOfGeometries)
	{
		if(str) free(str);
		str = strdup_printf(("Geometry selected Potential energy =  %0.4f"), molecularDynamics->potentialEnergy);
		//redrawMolecule(&molecularDynamics->qmModel->molecule,str);
		printf("%s\n",str);
		iSel++;
		geometries[iSel] = malloc(sizeof(QuantumMechanicsModel));
		*geometries[iSel] = molecularDynamics->qmModel->klass->copy(molecularDynamics->qmModel);
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
		molecularDynamics->qmModel->firstRun = FALSE;
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
			//redrawMolecule(&molecularDynamics->qmModel->molecule,str);
			printf("%s\n",str);
			updateNumber = 0;
		}
		saveProperties(molecularDynamics, n0+i+1, i+1," Heating");
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
		molecularDynamics->qmModel->firstRun = FALSE;
		molecularDynamics->temperature = currentTemp;
		/* rescaleVelocities(molecularDynamics);*/
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
			//redrawMolecule(&molecularDynamics->qmModel->molecule,str);
			printf("%s\n",str);
			updateNumber = 0;
		}
		saveProperties(molecularDynamics, n0+i+1, i+1, " Equilibrium");
	}
	updateNumber = molecularDynamics->updateFrequency;

	currentTemp = runTemperature;
	molecularDynamics->temperature = currentTemp;
	rescaleVelocities(molecularDynamics);
	updateNumber = molecularDynamics->updateFrequency;
	n0 += numberOfEquiSteps;
	/* newProperties(molecularDynamics," ----> Runing");*/
	if(str) free(str);
	str = strdup_printf(("Geometry selected Potential energy =  %0.4f"), molecularDynamics->potentialEnergy);
	//redrawMolecule(&molecularDynamics->qmModel->molecule,str);
	printf("%s\n",str);
	if(numberOfGeometries>2) stepSel = numberOfRunSteps/numberOfGeometries;
	else stepSel = numberOfRunSteps;
	/* printf("Isel = %d\n",stepSel);*/
	esum  = 0;
	e2sum = 0;
	for (i = 0; i < numberOfRunSteps; i++ )
	{
		molecularDynamics->temperature = currentTemp;
		applyOneStep(molecularDynamics,i);
		molecularDynamics->qmModel->firstRun = FALSE;
		applyThermostat(molecularDynamics);
		esum  += molecularDynamics->totalEnergy;
		e2sum += molecularDynamics->totalEnergy*molecularDynamics->totalEnergy;
		if (++updateNumber >= molecularDynamics->updateFrequency )
		{
			if(str) free(str);
			str = strdup_printf(("MD Running: %0.2f fs, T = %0.2f K  T(t) = %8.2f K Kin = %0.4f Pot =  %0.4f Tot =  %0.4f Eav = %0.4f sigE = %0.4f Eva-E0(cm^-1) = %0.2f"), 
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
			//redrawMolecule(&molecularDynamics->qmModel->molecule,str);
			printf("%s\n",str);
			updateNumber = 0;
			saveTrajectory(molecularDynamics, i+1);
		}
		if((i+1)%stepSel==0 && (iSel+1)<numberOfGeometries)
		{
			if(str) free(str);
			str = strdup_printf(("Geometry selected Potential energy =  %0.4f"), molecularDynamics->potentialEnergy);
			//redrawMolecule(&molecularDynamics->qmModel->molecule,str);
			printf("%s\n",str);
			iSel++;
			geometries[iSel] = malloc(sizeof(QuantumMechanicsModel));
			*geometries[iSel] = molecularDynamics->qmModel->klass->copy(molecularDynamics->qmModel);
		}
		saveProperties(molecularDynamics, n0+i+1, i+1," Running");
	}
	if(iSel<numberOfGeometries-1)
	{
		if(str) free(str);
		str = strdup_printf(("Geometry selected Potential energy =  %0.4f"), molecularDynamics->potentialEnergy);
		//redrawMolecule(&molecularDynamics->qmModel->molecule,str);
		printf("%s\n",str);
		iSel++;
		geometries[iSel] = malloc(sizeof(QuantumMechanicsModel));
		*geometries[iSel] = molecularDynamics->qmModel->klass->copy(molecularDynamics->qmModel);
	}

	updateNumber = molecularDynamics->updateFrequency;
	n0 += numberOfRunSteps;

	molecularDynamics->qmModel->klass->calculateGradient(molecularDynamics->qmModel);
        gradientNorm = 0;
	for (i = 0; i < molecularDynamics->numberOfAtoms; i++)
		for ( j = 0; j < 3; j++)
                        gradientNorm += 
				molecularDynamics->qmModel->molecule.atoms[i].gradient[j] * 
				molecularDynamics->qmModel->molecule.atoms[i].gradient[j]; 

        gradientNorm = sqrt( gradientNorm );
	if(str) free(str);
	str = strdup_printf(("End of MD Simulation. Gradient = %f Ekin = %f (Kcal/mol) EPot =  %0.4f ETot =  %0.4f T(t) = %0.2f"),
			(double)gradientNorm,
			molecularDynamics->kineticEnergy,
			molecularDynamics->potentialEnergy,
			molecularDynamics->totalEnergy,
			molecularDynamics->kelvin 
			); 
	//redrawMolecule(&molecularDynamics->qmModel->molecule,str);
	printf("%s\n",str);
	free(str);
	if(molecularDynamics->fileTraj)fclose(molecularDynamics->fileTraj);
	if(molecularDynamics->fileProp)fclose(molecularDynamics->fileProp);
	freeQuantumMechanicsMD(molecularDynamics);
	return geometries;
}
/**********************************************************************/
static void printGeometryAndVelocities(QuantumMechanicsMD* molecularDynamics, char* title)
{
	fprintf(stdout,"========================================================================================================================\n");
	fprintf(stdout,"#  Geometry and velocities at %s ; T0(K) = %0.2f\n", title, molecularDynamics->kelvin);
	molecularDynamics->qmModel->molecule.klass->addGeometry(& molecularDynamics->qmModel->molecule,stdout);
	molecularDynamics->qmModel->molecule.klass->addVelocities(& molecularDynamics->qmModel->molecule,stdout);
	fprintf(stdout,"========================================================================================================================\n");
}
/**********************************************************************/
void	runQuantumMechanicsMD(
		QuantumMechanicsMD* molecularDynamics, QuantumMechanicsModel* qmModel, 
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
	double e0  = 0;
	double esum  = 0;
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

	if(qmModel->molecule.nAtoms<1) return;

	molecularDynamics->qmModel = qmModel;
	molecularDynamics->numberOfAtoms = qmModel->molecule.nAtoms;
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

	initMD(molecularDynamics,currentTemp,stepSize,integratorType, thermostat, friction, omegaMax, Nf, collide, qNH, fileNameTraj, fileNameProp, numberOfRunSteps,0);
	molecularDynamics->qmModel->klass->calculateGradient(molecularDynamics->qmModel);

	molecularDynamics->temperature = heatTemperature;
	if( numberOfHeatSteps>0) rescaleVelocities(molecularDynamics);

	computeEnergies(molecularDynamics);
	e0 = molecularDynamics->potentialEnergy;
	printf("E0 = The first potential energy in kcal = %f\n",e0); 
	if( numberOfHeatSteps>0) printGeometryAndVelocities(molecularDynamics, "beginning of Heating stage");

	currentTemp = heatTemperature;
	n0 = 0;
	newProperties(molecularDynamics," ");
	/*newProperties(molecularDynamics," ----> Heating");*/
	for (i = 0; i < numberOfHeatSteps; i++ )
	{
		molecularDynamics->temperature = currentTemp;
		applyOneStep(molecularDynamics,i);
		molecularDynamics->qmModel->firstRun = FALSE;
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
			//redrawMolecule(&molecularDynamics->qmModel->molecule,str);
			printf("%s\n",str);
			updateNumber = 0;
		}
		saveProperties(molecularDynamics, n0+i+1, i+1," Heating");
	}

	currentTemp = equiTemperature;
	molecularDynamics->temperature = currentTemp;
	if( numberOfHeatSteps>0) rescaleVelocities(molecularDynamics);
	updateNumber = molecularDynamics->updateFrequency;
	n0 += numberOfHeatSteps;
	/* newProperties(molecularDynamics," ----> Equilibrium");*/
	if(numberOfEquiSteps>0) printGeometryAndVelocities(molecularDynamics, "beginning of Equilibrium stage");
	for (i = 0; i < numberOfEquiSteps; i++ )
	{
		molecularDynamics->temperature = currentTemp;
		applyOneStep(molecularDynamics,i);
		molecularDynamics->qmModel->firstRun = FALSE;
		molecularDynamics->temperature = currentTemp;
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
			//redrawMolecule(&molecularDynamics->qmModel->molecule,str);
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
	esum  = 0;
	e2sum = 0;
	if(numberOfRunSteps>0) printGeometryAndVelocities(molecularDynamics, "beginning of Production stage");
	for (i = 0; i < numberOfRunSteps; i++ )
	{
		molecularDynamics->temperature = currentTemp;
		applyOneStep(molecularDynamics,i);
		molecularDynamics->qmModel->firstRun = FALSE;
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
			//redrawMolecule(&molecularDynamics->qmModel->molecule,str);
			printf("%s\n",str);
			updateNumber = 0;
			saveTrajectory(molecularDynamics, i+1);
		}
		saveProperties(molecularDynamics, n0+i+1, i+1," Running");
	}
	if(numberOfCoolSteps>0) printGeometryAndVelocities(molecularDynamics, "the begining of Cooling stage");
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
		molecularDynamics->qmModel->firstRun = FALSE;
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
			//redrawMolecule(&molecularDynamics->qmModel->molecule,str);
			printf("%s\n",str);
			updateNumber = 0;
		}
		saveProperties(molecularDynamics, n0+i+1, i+1," Cooling");
	}
	molecularDynamics->qmModel->klass->calculateGradient(molecularDynamics->qmModel);
        gradientNorm = 0;
	for (i = 0; i < molecularDynamics->numberOfAtoms; i++)
		for ( j = 0; j < 3; j++)
                        gradientNorm += 
				molecularDynamics->qmModel->molecule.atoms[i].gradient[j] * 
				molecularDynamics->qmModel->molecule.atoms[i].gradient[j]; 

        gradientNorm = sqrt( gradientNorm );
	if(str) free(str);
	str = strdup_printf(("End of MD Simulation. Gradient = %f Ekin = %f (Kcal/mol) EPot =  %0.4f ETot =  %0.4f T(t) = %0.2f"),
			(double)gradientNorm,
			molecularDynamics->kineticEnergy,
			molecularDynamics->potentialEnergy,
			molecularDynamics->totalEnergy,
			molecularDynamics->kelvin 
			); 
	//redrawMolecule(&molecularDynamics->qmModel->molecule,str);
	printf("%s\n",str);
	free(str);
	printGeometryAndVelocities(molecularDynamics, "the end of simulation");
	if(molecularDynamics->fileTraj)fclose(molecularDynamics->fileTraj);
	if(molecularDynamics->fileProp)fclose(molecularDynamics->fileProp);
	freeQuantumMechanicsMD(molecularDynamics);
}
/*********************************************************************************/
static void initNH(QuantumMechanicsMD* molecularDynamics, double qNH)
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
static void initSD(QuantumMechanicsMD* molecularDynamics, double friction)
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
/*********************************************************************************/
static void removeTranslationAndRotationTheta(QuantumMechanicsMD* molecularDynamics)
{
	Molecule* molecule = &molecularDynamics->qmModel->molecule;
	double* theta = molecularDynamics->theta;
	molecule->klass->removeTranslationAndRotationAcceleration(molecule,theta);
}
/*********************************************************************************/
static void compteThetaQTB(QuantumMechanicsMD* molecularDynamics)
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
		molecularDynamics->theta[3*i+j] *= sigma/sqrt(molecularDynamics->qmModel->molecule.atoms[i].mass); 
	}
	removeTranslationAndRotationTheta(molecularDynamics);
}
/*********************************************************************************/
static void resetQTB(QuantumMechanicsMD* molecularDynamics)
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

	compteThetaQTB(molecularDynamics);

}
/*********************************************************************************/
/* omegaMax in cm-1 */
static void initQTB(QuantumMechanicsMD* molecularDynamics, double omegaMax, double friction, int Nf)
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
	else molecularDynamics->friction = friction/1000.0/fsInAKMA;

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
static void updateQTB(QuantumMechanicsMD* molecularDynamics)
{
	int i,j,k;

	if(molecularDynamics->temperature<=0) return;
	compteThetaQTB(molecularDynamics);
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
static void resetLangevin(QuantumMechanicsMD* molecularDynamics)
{
	if(molecularDynamics->integratorType != LANGEVIN) return;
	updateLangevin(molecularDynamics);
}
/*********************************************************************************/
/* omegaMax in cm-1 */
static void initLangevin(QuantumMechanicsMD* molecularDynamics, double friction)
{
	if(molecularDynamics->integratorType != LANGEVIN) return;

	if(friction<0) friction = 40;
	molecularDynamics->friction = friction/1000/fsInAKMA;
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
static void updateLangevin(QuantumMechanicsMD* molecularDynamics)
{
	int i,j;
	double sigma;
	/* update theta */
	/* thetaProg = thetaPaper*sigma/sqrt(m) */
	if(molecularDynamics->integratorType != LANGEVIN) return;
	sigma=sqrt(6.*molecularDynamics->friction*Kb*molecularDynamics->temperature/molecularDynamics->dt);
	for(i=0;i<molecularDynamics->numberOfAtoms;i++)
	for(j=0;j<3;j++)
	{
		molecularDynamics->theta[3*i+j] = normal();
		molecularDynamics->theta[3*i+j] *= sigma/sqrt(molecularDynamics->qmModel->molecule.atoms[i].mass); 
	}
}
/*********************************************************************************/
/*
static void printTranslation(QuantumMechanicsMD* molecularDynamics)
{
	double vtot[3] = {0,0,0};
	int i;
	int j;
	double mass = 1.0;
	for ( j = 0; j < 3; j++)
		vtot[j] = 0;
	for ( i = 0; i < molecularDynamics->numberOfAtoms; i++)
	{
		mass = molecularDynamics->qmModel->molecule.atoms[i].mass;
		for ( j = 0; j < 3; j++)
		{
			vtot[j] += mass*molecularDynamics->qmModel->molecule.atoms[i].velocity[j];
		}
	}
	printf("Trans velocity = %f %f %f\n",vtot[0], vtot[1], vtot[2]);
}
*/
/*********************************************************************************/
static void removeTranslationAndRotation(QuantumMechanicsMD* molecularDynamics)
{
	Molecule* molecule = &molecularDynamics->qmModel->molecule;
	molecule->klass->removeTranslationAndRotation(molecule); 
}
/*********************************************************************************/
static void initMD(QuantumMechanicsMD* molecularDynamics, double temperature, double stepSize, MDIntegratorType integratorType, MDThermostatType thermostat, double friction, double omegaMax, int Nf, double collide, double qNH, char* fileNameTraj, char* fileNameProp, int numberOfRunSteps, int index)
{
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
	if(molecularDynamics->qmModel->molecule.constraints!=NOCONSTRAINTS)
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

	initSD(molecularDynamics, friction);
	initNH(molecularDynamics,qNH);
	initQTB(molecularDynamics,omegaMax, friction, Nf);
	initLangevin(molecularDynamics,friction);


	molecularDynamics->qmModel->klass->calculateGradient(molecularDynamics->qmModel);
	for ( i = 0; i < molecularDynamics->numberOfAtoms; i++)
	{
		double m = molecularDynamics->qmModel->molecule.atoms[i].mass;
		for ( j = 0; j < 3; j++)
			molecularDynamics->a[i][j] = -molecularDynamics->qmModel->molecule.atoms[i].gradient[j]/m;
		if(molecularDynamics->aold)
			for ( j = 0; j < 3; j++)
				molecularDynamics->aold[i][j]  = molecularDynamics->a[i][j];
	}
	if(molecularDynamics->qmModel->molecule.klass->setMaxwellVelocitiesIfNull(&molecularDynamics->qmModel->molecule, temperature)) rescaleVelocities(molecularDynamics);
	removeTranslationAndRotation(molecularDynamics);
#ifdef DEBUG
	printf("nfree =%d\n",molecularDynamics->qmModel->molecule.nfree);
#endif
}
/*********************************************************************************/
static void rescaleVelocities(QuantumMechanicsMD* molecularDynamics)
{
	/* berendsen(molecularDynamics);*/
	scaleV(molecularDynamics);
	resetQTB(molecularDynamics);
	resetLangevin(molecularDynamics);
}
/*********************************************************************************/
static void scaleV(QuantumMechanicsMD* molecularDynamics)
{
	int i;
	int j;
	double ekin = 0;
	double kelvin = 0;
	int nfree = molecularDynamics->qmModel->molecule.nFree;
	double scale = 1.0;
	double mass = 1.0;
	if(molecularDynamics->temperature<=0) return;
	if(nfree<1) return;
	for ( i = 0; i < molecularDynamics->numberOfAtoms; i++)
	{
		mass = molecularDynamics->qmModel->molecule.atoms[i].mass;
		for ( j = 0; j < 3; j++)
			ekin += molecularDynamics->qmModel->molecule.atoms[i].velocity[j]*molecularDynamics->qmModel->molecule.atoms[i].velocity[j]*mass;
	}
	/*
	ekin /= 2;
	kelvin = 2* ekin / ( nfree * Kb);
	*/
	kelvin = ekin / ( nfree * Kb);
	scale = sqrt(molecularDynamics->temperature/kelvin);
#ifdef DEBUG
	printf("temp = %f kelvin = %f scale = %f\n",molecularDynamics->temperature, kelvin, scale);
#endif
	for ( i = 0; i < molecularDynamics->numberOfAtoms; i++)
		if(molecularDynamics->qmModel->molecule.atoms[i].variable)
		for ( j = 0; j < 3; j++)
			molecularDynamics->qmModel->molecule.atoms[i].velocity[j] *= scale;
	removeTranslationAndRotation(molecularDynamics);

/*
	ekin = 0;
	for ( i = 0; i < molecularDynamics->numberOfAtoms; i++)
	{
		mass = molecularDynamics->qmModel->molecule.atoms[i].mass;
		for ( j = 0; j < 3; j++)
			ekin += molecularDynamics->qmModel->molecule.atoms[i].velocity[j]*molecularDynamics->qmModel->molecule.atoms[i].velocity[j]*mass;
	}
	kelvin = ekin / ( nfree * Kb);
	scale = sqrt(molecularDynamics->temperature/kelvin);
	printf("Rem temp = %f kelvin = %f scale = %f\n",molecularDynamics->temperature, kelvin, scale);
*/
}
/*********************************************************************************/
static void berendsen(QuantumMechanicsMD* molecularDynamics)
{
	int i;
	int j;
	double ekin = 0;
	double kelvin = 0;
	int nfree = molecularDynamics->qmModel->molecule.nFree;
	double scale = 1.0;
	double dt = molecularDynamics->dt;
	double tautemp = 1.0/(molecularDynamics->collide)*1000*fsInAKMA;
	double mass = 1.0;
	if(molecularDynamics->temperature<=0) return;
	if(nfree<1) return;
	for ( i = 0; i < molecularDynamics->numberOfAtoms; i++)
	{
		mass = molecularDynamics->qmModel->molecule.atoms[i].mass;
		for ( j = 0; j < 3; j++)
			ekin += molecularDynamics->qmModel->molecule.atoms[i].velocity[j]*molecularDynamics->qmModel->molecule.atoms[i].velocity[j]*mass;
	}
	/*
	ekin /= 2;
	kelvin = 2* ekin / ( nfree * Kb);
	*/
	kelvin = ekin / ( nfree * Kb);
	/* if(tautemp>dt) tautemp = dt;*/
	scale = sqrt(1.0 + (dt/tautemp)*(molecularDynamics->temperature/kelvin-1.0));
#ifdef DEBUG
	printf("temp = %f kelvin = %f scale = %f\n",molecularDynamics->temperature, kelvin, scale);
#endif
	for ( i = 0; i < molecularDynamics->numberOfAtoms; i++)
		if(molecularDynamics->qmModel->molecule.atoms[i].variable)
		for ( j = 0; j < 3; j++)
			molecularDynamics->qmModel->molecule.atoms[i].velocity[j] *= scale;
	removeTranslationAndRotation(molecularDynamics);
}
/*********************************************************************************/
static void andersen(QuantumMechanicsMD* molecularDynamics)
{
	int i;
	double tau = 1.0/molecularDynamics->collide*1000*fsInAKMA; /* in fs */
	double rate;
	if(molecularDynamics->temperature<=0) return;
	if(molecularDynamics->numberOfAtoms<1) return;

	rate = molecularDynamics->dt / tau;
	rate /= pow(molecularDynamics->nvariables,2.0/3.0);

	for ( i = 0; i < molecularDynamics->numberOfAtoms; i++)
	{
		double trial = drandom();
		double m = molecularDynamics->qmModel->molecule.atoms[i].mass;
		if(trial<rate)
		{
/*
			double speed = maxwel(
					molecularDynamics->qmModel->molecule.atoms[i].mass,
					molecularDynamics->temperature
					);
			getRandVect(speed, molecularDynamics->qmModel->molecule.atoms[i].velocity);
*/
			double speed = sqrt(Kb* molecularDynamics->temperature/m);
                	double pnorm = normal();
			molecularDynamics->qmModel->molecule.atoms[i].velocity[0] = pnorm*speed;
                	pnorm = normal();
			molecularDynamics->qmModel->molecule.atoms[i].velocity[1] = pnorm*speed;
                	pnorm = normal();
			molecularDynamics->qmModel->molecule.atoms[i].velocity[2] = pnorm*speed;
		}
	}
}
/*********************************************************************************/
static void bussi(QuantumMechanicsMD* molecularDynamics)
{
	int nfree = molecularDynamics->qmModel->molecule.nFree;
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
	int i,j;
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
	for ( i = 0; i < molecularDynamics->numberOfAtoms; i++)
		if(molecularDynamics->qmModel->molecule.atoms[i].variable)
		for ( j = 0; j < 3; j++)
			molecularDynamics->qmModel->molecule.atoms[i].velocity[j] *= scale;
	removeTranslationAndRotation(molecularDynamics);
}
/*********************************************************************************/
static void nose_hoover(QuantumMechanicsMD* molecularDynamics)
{
	int nfree = molecularDynamics->qmModel->molecule.nFree;
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
	molecularDynamics->gNH[0] = (2.0*ekin-molecularDynamics->qmModel->molecule.nFree*kT) / molecularDynamics->qNH[0];
	molecularDynamics->vNH[0] = molecularDynamics->vNH[0] + molecularDynamics->gNH[0]*molecularDynamics->dt_4;
	molecularDynamics->vNH[0] = molecularDynamics->vNH[0] * exp(-molecularDynamics->vNH[1]*molecularDynamics->dt_8);
	molecularDynamics->xNH[0] = molecularDynamics->xNH[0] + molecularDynamics->vNH[0]*molecularDynamics->dt_2;
	molecularDynamics->xNH[1] = molecularDynamics->xNH[1] + molecularDynamics->vNH[1]*molecularDynamics->dt_2;
	//printf("vnH0 = %f\n",molecularDynamics->vNH[0]);
	scale = exp(-molecularDynamics->vNH[0]*molecularDynamics->dt_2);
	//printf("scale = %f\n",scale);
	for ( i = 0; i < molecularDynamics->numberOfAtoms; i++)
		for ( j = 0; j < 3; j++)
			molecularDynamics->qmModel->molecule.atoms[i].velocity[j] *= scale;
	//removeTranslationAndRotation(molecularDynamics);
	ekin = ekin * scale * scale;
	molecularDynamics->vNH[0] = molecularDynamics->vNH[0] * exp(-molecularDynamics->vNH[1]*molecularDynamics->dt_8);
	molecularDynamics->gNH[0] = (2.0*ekin-nfree*kT) /  molecularDynamics->qNH[0];
	molecularDynamics->vNH[0] = molecularDynamics->vNH[0] + molecularDynamics->gNH[0]*molecularDynamics->dt_4;
	molecularDynamics->vNH[0] = molecularDynamics->vNH[0] * exp(-molecularDynamics->vNH[1]*molecularDynamics->dt_8);
	molecularDynamics->gNH[1] = ( molecularDynamics->qNH[0]*molecularDynamics->vNH[0]*molecularDynamics->vNH[0]-kT) /  molecularDynamics->qNH[1];
	molecularDynamics->vNH[1] = molecularDynamics->vNH[1] + molecularDynamics->gNH[1]*molecularDynamics->dt_4;

}
/*********************************************************************************/
static void newAccelaration(QuantumMechanicsMD* molecularDynamics)
{
	int i;
	int j;
	molecularDynamics->qmModel->klass->calculateGradient(molecularDynamics->qmModel);
	for ( i = 0; i < molecularDynamics->numberOfAtoms; i++)
	{
		double m = molecularDynamics->qmModel->molecule.atoms[i].mass;
		if(molecularDynamics->aold)
			for ( j = 0; j < 3; j++)
				molecularDynamics->aold[i][j]  = molecularDynamics->a[i][j];

		for ( j = 0; j < 3; j++)
			molecularDynamics->a[i][j] = -molecularDynamics->qmModel->molecule.atoms[i].gradient[j]/m;
	}
}
/*********************************************************************************/
static void computeEnergies(QuantumMechanicsMD* molecularDynamics)
{
	molecularDynamics->kineticEnergy = getEKin(molecularDynamics);
	molecularDynamics->potentialEnergy = molecularDynamics->qmModel->molecule.potentialEnergy;
	molecularDynamics->totalEnergy = molecularDynamics->kineticEnergy + molecularDynamics->potentialEnergy;
	molecularDynamics->kelvin = getKelvin(molecularDynamics);
}
/*********************************************************************************/
static void applyThermostat(QuantumMechanicsMD* molecularDynamics)
{
	if(molecularDynamics->integratorType == STOCHASTIC) return;
	if(molecularDynamics->integratorType == QTB) return;
	if(molecularDynamics->integratorType == LANGEVIN) return;
	if(molecularDynamics->thermostat == ANDERSEN) andersen(molecularDynamics);
	if(molecularDynamics->thermostat == BERENDSEN) berendsen(molecularDynamics);
	if(molecularDynamics->thermostat == BUSSI) bussi(molecularDynamics);
}
/*********************************************************************************/
static void applyOneStep(QuantumMechanicsMD* molecularDynamics, int iStep)
{
	if(molecularDynamics->integratorType == VERLET) applyVerlet(molecularDynamics);
	else if(molecularDynamics->integratorType == BEEMAN) applyBeeman(molecularDynamics);
	else if(molecularDynamics->integratorType == STOCHASTIC) applyStochastic(molecularDynamics);
	else if(molecularDynamics->integratorType == LANGEVIN) 
	{
		updateLangevin(molecularDynamics);
		applyLangevin(molecularDynamics);
	}
	else {
		if((iStep+1)%molecularDynamics->M==0) updateQTB(molecularDynamics);
		applyQTB(molecularDynamics);
	}
	computeEnergies(molecularDynamics);
	/*
	printTranslation(molecularDynamics);
	printRotation(molecularDynamics);
	*/
	removeTranslationAndRotation(molecularDynamics);

}
/*********************************************************************************/
static void applyRattleFirstPortion(QuantumMechanicsMD* quantumMechanicsMD)
{
	int i;
	int k;
	int maxIter = 100;
	double omega = 1.2; 
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
	Molecule* m = &quantumMechanicsMD->qmModel->molecule;
	QuantumMechanicsModel* qmModel = quantumMechanicsMD->qmModel;
	double deltaMax = 0;

	if(qmModel->molecule.constraints==NOCONSTRAINTS) return;
	for (i = 0; i < quantumMechanicsMD->numberOfAtoms; i++)
	{
			quantumMechanicsMD->moved[i] = quantumMechanicsMD->qmModel->molecule.atoms[i].variable;
			quantumMechanicsMD->update[i] = FALSE;
	}
	maxIter *= quantumMechanicsMD->qmModel->molecule.numberOfRattleConstraintsTerms;
	do{
		nIter++;
		done=TRUE;
		deltaMax = 0;
		for (i = 0; i < quantumMechanicsMD->qmModel->molecule.numberOfRattleConstraintsTerms; i++)
		{
			a1 = (int)quantumMechanicsMD->qmModel->molecule.rattleConstraintsTerms[0][i];
			a2 = (int)quantumMechanicsMD->qmModel->molecule.rattleConstraintsTerms[1][i];
			if( !quantumMechanicsMD->moved[a1] && !quantumMechanicsMD->moved[a2] ) continue;
			r2ij = 0;
			for (k=0;k<3;k++)
			{
				d = m->atoms[a2].coordinates[k]-m->atoms[a1].coordinates[k];
				r2ij +=d*d;
			}
			delta = quantumMechanicsMD->qmModel->molecule.rattleConstraintsTerms[2][i]-r2ij;
			if(deltaMax<fabs(delta)) deltaMax = fabs(delta);
			if(fabs(delta)<=tolerance) continue;
			done = FALSE;
			quantumMechanicsMD->update[a1] = TRUE;
			quantumMechanicsMD->update[a2] = TRUE;
			/* here : rattle image for PBC, not yet implemented */
			dot = 0;
			for (k=0;k<3;k++)
			{
				d = m->atoms[a2].coordinates[k]-m->atoms[a1].coordinates[k];
				dot +=d*(quantumMechanicsMD->coordinatesOld[a2][k]-quantumMechanicsMD->coordinatesOld[a1][k]);
			}
			invMass1 = 1/m->atoms[a1].mass;
			invMass2 = 1/m->atoms[a2].mass;
		        term = omega*delta / (2.0*(invMass1+invMass2)*dot);
			for (k=0;k<3;k++)
			{
				terms[k] = (quantumMechanicsMD->coordinatesOld[a2][k]-quantumMechanicsMD->coordinatesOld[a1][k])*term;
			}
			for (k=0;k<3;k++) m->atoms[a1].coordinates[k] -= terms[k]*invMass1;
			for (k=0;k<3;k++) m->atoms[a2].coordinates[k] += terms[k]*invMass2;

			invMass1 /= quantumMechanicsMD->dt;
			invMass2 /= quantumMechanicsMD->dt;
			for (k=0;k<3;k++) quantumMechanicsMD->qmModel->molecule.atoms[a1].velocity[k] -= terms[k]*invMass1;
			for (k=0;k<3;k++) quantumMechanicsMD->qmModel->molecule.atoms[a2].velocity[k] += terms[k]*invMass2;
		}
		for (i = 0; i < quantumMechanicsMD->numberOfAtoms; i++)
		{
			quantumMechanicsMD->moved[i] = quantumMechanicsMD->update[i];
			quantumMechanicsMD->update[i] = FALSE;
		}
	}while(!done && nIter<maxIter);
	if(nIter>=maxIter && deltaMax>tolerance*10)
	{
		printf(("Rattle first portion : Warning, distance constraints not satisfied\n"));
	}
	for (i = 0; i <  quantumMechanicsMD->numberOfAtoms; i++)
	if(!m->atoms[i].variable)
	{
		for (k=0;k<3;k++)  quantumMechanicsMD->qmModel->molecule.atoms[i].velocity[k] = 0.0;
		for (k=0;k<3;k++)  for (k=0;k<3;k++) m->atoms[i].coordinates[k] =  quantumMechanicsMD->coordinatesOld[i][k];
	}

}
/*********************************************************************************/
static void applyRattleSecondPortion(QuantumMechanicsMD* quantumMechanicsMD)
{
	int i;
	int k;
	int maxIter = 100;
	double omega = 1.2;
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
	Molecule* m = &quantumMechanicsMD->qmModel->molecule;
	QuantumMechanicsModel* qmModel = quantumMechanicsMD->qmModel;
	double deltaMax = 0;

	if(qmModel->molecule.constraints==NOCONSTRAINTS) return;
	tolerance /= quantumMechanicsMD->dt;
	for (i = 0; i < quantumMechanicsMD->numberOfAtoms; i++)
	{
			quantumMechanicsMD->moved[i] = quantumMechanicsMD->qmModel->molecule.atoms[i].variable;
			quantumMechanicsMD->update[i] = FALSE;
	}
	maxIter *= quantumMechanicsMD->qmModel->molecule.numberOfRattleConstraintsTerms;
	do{
		nIter++;
		done=TRUE;
		deltaMax = 0;
		for (i = 0; i < quantumMechanicsMD->qmModel->molecule.numberOfRattleConstraintsTerms; i++)
		{
			a1 = (int)quantumMechanicsMD->qmModel->molecule.rattleConstraintsTerms[0][i];
			a2 = (int)quantumMechanicsMD->qmModel->molecule.rattleConstraintsTerms[1][i];
			r2ij = quantumMechanicsMD->qmModel->molecule.rattleConstraintsTerms[2][i];
			if( !quantumMechanicsMD->moved[a1] && !quantumMechanicsMD->moved[a2] ) continue;
			/* here : rattle image for PBC, not yet implemented */
			dot = 0;
			for (k=0;k<3;k++)
			{
				d = m->atoms[a2].coordinates[k]-m->atoms[a1].coordinates[k];
				dot +=d*(quantumMechanicsMD->qmModel->molecule.atoms[a2].velocity[k]-quantumMechanicsMD->qmModel->molecule.atoms[a1].velocity[k]);
			}
			invMass1 = 1/quantumMechanicsMD->qmModel->molecule.atoms[a1].mass;
			invMass2 = 1/quantumMechanicsMD->qmModel->molecule.atoms[a2].mass;
		        term = -dot / ((invMass1+invMass2)*r2ij);
			if(deltaMax<fabs(term)) deltaMax = fabs(term);
			if(fabs(term)<=tolerance) continue;
			done = FALSE;
			quantumMechanicsMD->update[a1] = TRUE;
			quantumMechanicsMD->update[a2] = TRUE;
		        term *= omega;

			for (k=0;k<3;k++)
			{
				d = m->atoms[a2].coordinates[k]-m->atoms[a1].coordinates[k];
				terms[k] = d*term;
			}
			for (k=0;k<3;k++) quantumMechanicsMD->qmModel->molecule.atoms[a1].velocity[k] -= terms[k]*invMass1;
			for (k=0;k<3;k++) quantumMechanicsMD->qmModel->molecule.atoms[a2].velocity[k] += terms[k]*invMass2;
		}
		for (i = 0; i < quantumMechanicsMD->numberOfAtoms; i++)
		{
			quantumMechanicsMD->moved[i] = quantumMechanicsMD->update[i];
			quantumMechanicsMD->update[i] = FALSE;
		}
	}while(!done && nIter<maxIter);
	if(nIter>=maxIter && deltaMax>tolerance*10)
	{
		printf(("Rattle second portion : Warning, velocity constraints not satisfied\n"));
	}
	for (i = 0; i <  quantumMechanicsMD->numberOfAtoms; i++)
			if(!m->atoms[i].variable)
			for (k=0;k<3;k++)  quantumMechanicsMD->qmModel->molecule.atoms[i].velocity[k] = 0.0;
}
/*********************************************************************************/
static void applyVerlet(QuantumMechanicsMD* molecularDynamics)
{
	int i;
	int j;

	if(molecularDynamics->qmModel->molecule.constraints!=NOCONSTRAINTS)
	for (i = 0; i < molecularDynamics->numberOfAtoms; i++)
		for ( j = 0; j < 3; j++)
				molecularDynamics->coordinatesOld[i][j]= molecularDynamics->qmModel->molecule.atoms[i].coordinates[j];

	for (i = 0; i < molecularDynamics->numberOfAtoms; i++)
	{
		if(!molecularDynamics->qmModel->molecule.atoms[i].variable) continue;

		for ( j = 0; j < 3; j++)
		{
			molecularDynamics->qmModel->molecule.atoms[i].coordinates[j] += 
				molecularDynamics->qmModel->molecule.atoms[i].velocity[j] * molecularDynamics->dt +
				molecularDynamics->a[i][j]*molecularDynamics->dt2_2;	
		}
		for ( j = 0; j < 3; j++)
			molecularDynamics->qmModel->molecule.atoms[i].velocity[j] += molecularDynamics->a[i][j] * molecularDynamics->dt_2;
	}

	if(molecularDynamics->qmModel->molecule.constraints!=NOCONSTRAINTS) applyRattleFirstPortion(molecularDynamics);
	if(molecularDynamics->thermostat==NOSEHOOVER) nose_hoover(molecularDynamics);

	newAccelaration(molecularDynamics);

	for (i = 0; i < molecularDynamics->numberOfAtoms; i++)
		if(molecularDynamics->qmModel->molecule.atoms[i].variable)
		for ( j = 0; j < 3; j++)
			molecularDynamics->qmModel->molecule.atoms[i].velocity[j] += molecularDynamics->a[i][j] * molecularDynamics->dt_2;
	if(molecularDynamics->qmModel->molecule.constraints!=NOCONSTRAINTS) applyRattleSecondPortion(molecularDynamics);
}
/*********************************************************************************/
static void applyBeeman(QuantumMechanicsMD* molecularDynamics)
{
	int i;
	int j;
	double terms[3];
	if(molecularDynamics->qmModel->molecule.constraints!=NOCONSTRAINTS)
	for (i = 0; i < molecularDynamics->numberOfAtoms; i++)
		for ( j = 0; j < 3; j++)
				molecularDynamics->coordinatesOld[i][j]= molecularDynamics->qmModel->molecule.atoms[i].coordinates[j];

	for (i = 0; i < molecularDynamics->numberOfAtoms; i++)
	{
		if(!molecularDynamics->qmModel->molecule.atoms[i].variable) continue;

		for ( j = 0; j < 3; j++)
			terms[j] = 5.0*molecularDynamics->a[i][j]-molecularDynamics->aold[i][j];

		for ( j = 0; j < 3; j++)
		{
			molecularDynamics->qmModel->molecule.atoms[i].coordinates[j] += 
				molecularDynamics->qmModel->molecule.atoms[i].velocity[j] * molecularDynamics->dt +
				terms[j]*molecularDynamics->dt2_8;	
		}
		for ( j = 0; j < 3; j++)
			molecularDynamics->qmModel->molecule.atoms[i].velocity[j] += terms[j] * molecularDynamics->dt_8;
	}

	if(molecularDynamics->qmModel->molecule.constraints!=NOCONSTRAINTS) applyRattleFirstPortion(molecularDynamics);
	if(molecularDynamics->thermostat==NOSEHOOVER) nose_hoover(molecularDynamics);

	newAccelaration(molecularDynamics);

	for (i = 0; i < molecularDynamics->numberOfAtoms; i++)
		if(molecularDynamics->qmModel->molecule.atoms[i].variable)
		for ( j = 0; j < 3; j++)
			molecularDynamics->qmModel->molecule.atoms[i].velocity[j] += (3.0*molecularDynamics->a[i][j]+molecularDynamics->aold[i][j]) * molecularDynamics->dt_8;
	if(molecularDynamics->qmModel->molecule.constraints!=NOCONSTRAINTS) applyRattleSecondPortion(molecularDynamics);
}
/*********************************************************************************/
static void newProperties(QuantumMechanicsMD* molecularDynamics, char* comments)
{
	if( molecularDynamics->fileProp == NULL) return;
	fprintf(molecularDynamics->fileProp,"Time0(fs)\tTime(fs)\tTotal Energy(Kcal/mol)\tPotential Energy(kcal/mol)\tKinetic Energy(Kcal/mol)\tT(t) (K)\tTaver(K)\tsigma(T)(K)\tIndex\tmuX\tmuY\tmuZ");
	if(molecularDynamics->thermostat==NOSEHOOVER) fprintf(molecularDynamics->fileProp,"\tEtot+Etherm");
	if(comments) fprintf(molecularDynamics->fileProp,"%s\n", comments);
	else fprintf(molecularDynamics->fileProp,"\n");
}
/*********************************************************************************/
static void saveProperties(QuantumMechanicsMD* molecularDynamics, int iStep0, int iStep, char* comments)
{
	double dt = molecularDynamics->dt/(fsInAKMA);
	static double Ttot = 0;
	static double T2tot = 0;
	double Taver = 0;
	double T2aver = 0;
	double totalEnergy =  molecularDynamics->totalEnergy;

	if( molecularDynamics->thermostat==NOSEHOOVER)
	{
		int i;
		double kT = Kb* molecularDynamics->temperature;
		double e = molecularDynamics->vNH[0]*molecularDynamics->vNH[0]* molecularDynamics->qNH[0]/2 + (molecularDynamics->qmModel->molecule.nFree)*kT* molecularDynamics->xNH[0];
		for(i=1;i<MAXNH;i++) e += molecularDynamics->vNH[i]*molecularDynamics->vNH[i]* molecularDynamics->qNH[i]/2 + kT* molecularDynamics->xNH[i];
		
		totalEnergy += e;
	}

	if( molecularDynamics->fileProp == NULL) return;
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
	fprintf(molecularDynamics->fileProp,"%f\t%f\t%f\t", 
			molecularDynamics->qmModel->molecule.dipole[0], 
			molecularDynamics->qmModel->molecule.dipole[1], 
			molecularDynamics->qmModel->molecule.dipole[2]);
	if( molecularDynamics->thermostat==NOSEHOOVER) fprintf(molecularDynamics->fileProp,"%f\t",totalEnergy);
	if(comments) fprintf(molecularDynamics->fileProp,"%s\n", comments);
	else fprintf(molecularDynamics->fileProp,"\n");
}
/*********************************************************************************/
static void saveTrajectory(QuantumMechanicsMD* molecularDynamics, int iStep)
{
	double dt = molecularDynamics->dt/(fsInAKMA);
	int i;
	if( molecularDynamics->fileTraj == NULL) return;

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
				molecularDynamics->qmModel->molecule.atoms[i].prop.symbol,
				molecularDynamics->qmModel->molecule.atoms[i].coordinates[0],
				molecularDynamics->qmModel->molecule.atoms[i].coordinates[1],
				molecularDynamics->qmModel->molecule.atoms[i].coordinates[2],
				molecularDynamics->qmModel->molecule.atoms[i].velocity[0],
				molecularDynamics->qmModel->molecule.atoms[i].velocity[1],
				molecularDynamics->qmModel->molecule.atoms[i].velocity[2],
				molecularDynamics->qmModel->molecule.atoms[i].charge,
				molecularDynamics->qmModel->molecule.atoms[i].mmType,
				molecularDynamics->qmModel->molecule.atoms[i].pdbType,
				molecularDynamics->qmModel->molecule.atoms[i].residueName,
				molecularDynamics->qmModel->molecule.atoms[i].residueNumber,
				molecularDynamics->qmModel->molecule.atoms[i].variable
				);
	}
}

/**********************************************************************/
void	freeQuantumMechanicsMD(QuantumMechanicsMD* molecularDynamics)
{

	molecularDynamics->qmModel = NULL;
	molecularDynamics->numberOfAtoms = 0;
	molecularDynamics->updateFrequency = 0;
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
		int i,j;
		for(i=0;i<molecularDynamics->numberOfAtoms;i++)
		for(j=0;j<3;j++)
			if(molecularDynamics->rnoise[3*i+j]) free(molecularDynamics->rnoise[3*i+j]);
		free(molecularDynamics->rnoise);
	}
}
/********************************************************************************/
static double getEKin(QuantumMechanicsMD* molecularDynamics)
{
	double ekin = 0;
	int i;
	int j;
	double mass;
	for ( i = 0; i < molecularDynamics->numberOfAtoms; i++)
	{
		mass = molecularDynamics->qmModel->molecule.atoms[i].mass;
		for ( j = 0; j < 3; j++)
			ekin += molecularDynamics->qmModel->molecule.atoms[i].velocity[j]*molecularDynamics->qmModel->molecule.atoms[i].velocity[j]*
				mass;
	}
	return ekin/2;
}
/********************************************************************************/
static double getKelvin(QuantumMechanicsMD* molecularDynamics)
{
	int nfree = molecularDynamics->qmModel->molecule.nFree;
	/* printf("nfree = %d\n",nfree);*/
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
static void getsFrictionalAndRandomForce(QuantumMechanicsMD* molecularDynamics)
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
               		ktm = Kb * molecularDynamics->temperature / molecularDynamics->qmModel->molecule.atoms[i].mass;
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
/*********************************************************************************/
static void applyStochastic(QuantumMechanicsMD* molecularDynamics)
{
	double* positionFriction = molecularDynamics->positionFriction;
	double* velocityFriction = molecularDynamics->velocityFriction;
	double* accelarationFriction = molecularDynamics->accelarationFriction;
	double** positionRandom = molecularDynamics->positionRandom;
	double** velocityRandom = molecularDynamics->velocityRandom;
	double**a = molecularDynamics->a;
	
	int n = molecularDynamics->numberOfAtoms;
	int i;
	int j;
	Atom* atoms = molecularDynamics->qmModel->molecule.atoms;

	getsFrictionalAndRandomForce(molecularDynamics);

	if(molecularDynamics->qmModel->molecule.constraints!=NOCONSTRAINTS)
	for (i = 0; i < n; i++)
		for ( j = 0; j < 3; j++)
				molecularDynamics->coordinatesOld[i][j]= molecularDynamics->qmModel->molecule.atoms[i].coordinates[j];

	for(i=0;i<n;i++)
	{
		if(!molecularDynamics->qmModel->molecule.atoms[i].variable) continue;
		for(j=0;j<3;j++)
			atoms[i].coordinates[j] += molecularDynamics->qmModel->molecule.atoms[i].velocity[j]*velocityFriction[i] + a[i][j]*accelarationFriction[i] + positionRandom[i][j];
		for(j=0;j<3;j++)
			molecularDynamics->qmModel->molecule.atoms[i].velocity[j] = molecularDynamics->qmModel->molecule.atoms[i].velocity[j]*positionFriction[i] + 0.5*a[i][j]*velocityFriction[i];
	}

	if(molecularDynamics->qmModel->molecule.constraints!=NOCONSTRAINTS) applyRattleFirstPortion(molecularDynamics);

	newAccelaration(molecularDynamics);

	for (i = 0; i < n; i++)
		if(molecularDynamics->qmModel->molecule.atoms[i].variable)
		for ( j = 0; j < 3; j++)
			molecularDynamics->qmModel->molecule.atoms[i].velocity[j] += 0.5*a[i][j]*velocityFriction[i] + velocityRandom[i][j];
	if(molecularDynamics->qmModel->molecule.constraints!=NOCONSTRAINTS) applyRattleSecondPortion(molecularDynamics);
	removeTranslationAndRotation(molecularDynamics);
	computeEnergies(molecularDynamics);
}
/*********************************************************************************/
static void applyQTB(QuantumMechanicsMD* molecularDynamics)
{
	
	int n = molecularDynamics->numberOfAtoms;
	int i;
	int j;
	double gp = 1/(1+molecularDynamics->friction*molecularDynamics->dt_2);
	double gm = (1-molecularDynamics->friction*molecularDynamics->dt_2)*gp;

	/* printf("gm = %f gp =%f\n",gm,gp);*/
	if(molecularDynamics->qmModel->molecule.constraints!=NOCONSTRAINTS)
	for (i = 0; i < n; i++)
		for ( j = 0; j < 3; j++)
				molecularDynamics->coordinatesOld[i][j]= molecularDynamics->qmModel->molecule.atoms[i].coordinates[j];

	for (i = 0; i < n; i++)
		if(molecularDynamics->qmModel->molecule.atoms[i].variable)
		for ( j = 0; j < 3; j++)
		{
			molecularDynamics->qmModel->molecule.atoms[i].coordinates[j] += 
			molecularDynamics->qmModel->molecule.atoms[i].velocity[j]*molecularDynamics->dt +
			(molecularDynamics->a[i][j]+molecularDynamics->theta[3*i+j]-molecularDynamics->qmModel->molecule.atoms[i].velocity[j]*molecularDynamics->friction)*molecularDynamics->dt2_2;	
		}
	for (i = 0; i < n; i++)
		if(molecularDynamics->qmModel->molecule.atoms[i].variable)
		for ( j = 0; j < 3; j++)
			molecularDynamics->qmModel->molecule.atoms[i].velocity[j] = 
			gm*molecularDynamics->qmModel->molecule.atoms[i].velocity[j] +
			gp*(molecularDynamics->a[i][j]+molecularDynamics->theta[3*i+j])*molecularDynamics->dt_2;
            
	if(molecularDynamics->qmModel->molecule.constraints!=NOCONSTRAINTS) applyRattleFirstPortion(molecularDynamics);

	newAccelaration(molecularDynamics);

	for (i = 0; i < n; i++)
		if(molecularDynamics->qmModel->molecule.atoms[i].variable)
		for ( j = 0; j < 3; j++)
			molecularDynamics->qmModel->molecule.atoms[i].velocity[j] += gp*(molecularDynamics->a[i][j]+molecularDynamics->theta[3*i+j])*molecularDynamics->dt_2;

	if(molecularDynamics->qmModel->molecule.constraints!=NOCONSTRAINTS) applyRattleSecondPortion(molecularDynamics);
	removeTranslationAndRotation(molecularDynamics);
	computeEnergies(molecularDynamics);
}
/*********************************************************************************/
static void applyLangevin(QuantumMechanicsMD* molecularDynamics)
{
	
	int n = molecularDynamics->numberOfAtoms;
	int i;
	int j;
	double gp = 1/(1+molecularDynamics->friction*molecularDynamics->dt_2);
	double gm = (1-molecularDynamics->friction*molecularDynamics->dt_2)*gp;

	/* printf("gm = %f gp =%f\n",gm,gp);*/
	if(molecularDynamics->qmModel->molecule.constraints!=NOCONSTRAINTS)
	for (i = 0; i < n; i++)
		for ( j = 0; j < 3; j++)
				molecularDynamics->coordinatesOld[i][j]= molecularDynamics->qmModel->molecule.atoms[i].coordinates[j];

	for (i = 0; i < n; i++)
		if(molecularDynamics->qmModel->molecule.atoms[i].variable)
		for ( j = 0; j < 3; j++)
		{
			molecularDynamics->qmModel->molecule.atoms[i].coordinates[j] += 
			molecularDynamics->qmModel->molecule.atoms[i].velocity[j]*molecularDynamics->dt +
			(molecularDynamics->a[i][j]+molecularDynamics->theta[3*i+j]-molecularDynamics->qmModel->molecule.atoms[i].velocity[j]*molecularDynamics->friction)*molecularDynamics->dt2_2;	
		}
	for (i = 0; i < n; i++)
		if(molecularDynamics->qmModel->molecule.atoms[i].variable)
		for ( j = 0; j < 3; j++)
			molecularDynamics->qmModel->molecule.atoms[i].velocity[j] = 
			gm*molecularDynamics->qmModel->molecule.atoms[i].velocity[j] +
			gp*(molecularDynamics->a[i][j]+molecularDynamics->theta[3*i+j])*molecularDynamics->dt_2;
            
	if(molecularDynamics->qmModel->molecule.constraints!=NOCONSTRAINTS) applyRattleFirstPortion(molecularDynamics);

	newAccelaration(molecularDynamics);

	for (i = 0; i < n; i++)
		if(molecularDynamics->qmModel->molecule.atoms[i].variable)
		for ( j = 0; j < 3; j++)
			molecularDynamics->qmModel->molecule.atoms[i].velocity[j] += gp*(molecularDynamics->a[i][j]+molecularDynamics->theta[3*i+j])*molecularDynamics->dt_2;

	if(molecularDynamics->qmModel->molecule.constraints!=NOCONSTRAINTS) applyRattleSecondPortion(molecularDynamics);
	removeTranslationAndRotation(molecularDynamics);
	computeEnergies(molecularDynamics);
}
/*********************************************************************************/
