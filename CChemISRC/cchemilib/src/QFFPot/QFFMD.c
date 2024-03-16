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

/* QFFMD.c  */

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
#include "../QFFPot/QFFModel.h"
#include "../QFFPot/QFFMD.h"


/*********************************************************************************/
static void initMD(QFFMD* molecularDynamics, double temperature, double stepSize, MDIntegratorType integratorType, MDThermostatType thermostat, double friction, double omegaMax, int Nf, double collide, double qNH, char* fileNameTraj, char* fileNameProp, int numberOfRunSteps, int index);
static void berendsen(QFFMD* molecularDynamics);
static void scaleV(QFFMD* molecularDynamics);
static void bussi(QFFMD* molecularDynamics);
static void andersen(QFFMD* molecularDynamics);
static void nose_hoover(QFFMD* molecularDynamics);
static void rescaleVelocities(QFFMD* molecularDynamics);
static void computeEnergies(QFFMD* molecularDynamics);
static void applyOneStep(QFFMD* molecularDynamics, int iStep);
static void applyThermostat(QFFMD* molecularDynamics);
static void applyVerlet(QFFMD* molecularDynamics);
static void applyBeeman(QFFMD* molecularDynamics);
static void applyStochastic(QFFMD* molecularDynamics);
static void applyQTB(QFFMD* molecularDynamics);
static void updateQTB(QFFMD* molecularDynamics);
static void resetQTB(QFFMD* molecularDynamics);
static void applyLangevin(QFFMD* molecularDynamics);
static void updateLangevin(QFFMD* molecularDynamics);
static void resetLangevin(QFFMD* molecularDynamics);
static void newProperties(QFFMD* molecularDynamics, char* comments);
static void saveProperties(QFFMD* molecularDynamics, int iStep0, int iStep, char* comments);
static void saveTrajectory(QFFMD* molecularDynamics, int iStep);
static double getEKin(QFFMD* molecularDynamics);
static double getKelvin(QFFMD* molecularDynamics);
/*****************************************************************************************************************/
static boolean createNewInputFile(QFFModel* qffModel,char* inputFileName, int iStep)
{
	char t[BSIZE];
	FILE* fileIn = NULL;
	FILE* fileOut = NULL;
	char* fileNameMD = NULL;
        char* suff = getSuffixNameFile(inputFileName);
        fileNameMD = strdup_printf("%s_%d.ici",suff, iStep);
        free(suff);
	fileIn = fopen(inputFileName,"rb");
	fileOut = fopen(fileNameMD,"w");
        free(fileNameMD);
	if(!fileIn) return FALSE;
	if(!fileOut) return FALSE;
	qffModel->klass->printModesAndVelocities(qffModel,fileOut);
	while(!feof(fileIn))
	{
	 	if(!fgets(t,BSIZE,fileIn)) break;
		deleteFirstSpaces(t);
		if(t[0]!='#' && mystrcasestr(t,"integrator"))
		{
			sprintf(t,"integrator = 0\n");
		}
		if(t[0]!='#' && mystrcasestr(t,"equiTime"))
		{
			sprintf(t,"equiTime = 0\n");
		}
		if(t[0]!='#' && mystrcasestr(t,"numberOfGeometries"))
		{
			sprintf(t,"#numberOfGeometries=10\n");
		}
		fprintf(fileOut,"%s",t);
	}
	fclose(fileIn);
	fclose(fileOut);
	return TRUE;
}
/*****************************************************************************************************************/
void	runQFFMD(
		QFFMD* molecularDynamics, QFFModel* qffModel, 
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
		int numberOfGeometries,
		char* fileNameTraj,
		char* fileNameProp,
		char* inputFileName
		)
{
	int i;
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
	int iSel = 0;
	int stepSel = 1;
	FILE* logfile = stdout;
        char* suff = getSuffixNameFile(inputFileName);
	int j;
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
	 *       Internally, Atomic units are used:
	 *                                        
	*/ 

	/* fprintf(logfile,"basname = %s\n",g_path_get_basename(fileNameTraj));*/

	if(qffModel->molecule.nAtoms<1) return;

	molecularDynamics->qffModel = qffModel;
	molecularDynamics->numberOfModes = qffModel->vData.nFrequencies;
	molecularDynamics->updateFrequency = updateFrequency;

	currentTemp = heatTemperature/2;
	
	numberOfHeatSteps = heatTime/stepSize*1000;
	numberOfEquiSteps = equiTime/stepSize*1000;; 
	numberOfRunSteps = runTime/stepSize*1000;; 
	numberOfCoolSteps = coolTime/stepSize*1000;;

	if(numberOfGeometries>2) stepSel = numberOfRunSteps/(numberOfGeometries-1)+1;
        else stepSel = numberOfRunSteps;



	currentTemp = heatTemperature;
	if(numberOfHeatSteps==0) currentTemp = equiTemperature; 
	if(numberOfHeatSteps==0 && numberOfEquiSteps==0 ) currentTemp = runTemperature; 
	if(numberOfHeatSteps==0 && numberOfEquiSteps==0 && numberOfRunSteps==0 ) currentTemp = coolTemperature; 

	initMD(molecularDynamics,currentTemp,stepSize,integratorType, thermostat, friction, omegaMax, Nf, collide, qNH, fileNameTraj, fileNameProp, numberOfRunSteps,0);
	molecularDynamics->qffModel->klass->calculateGradient(molecularDynamics->qffModel);

	molecularDynamics->temperature = heatTemperature;
	if( numberOfHeatSteps>0) rescaleVelocities(molecularDynamics);

	computeEnergies(molecularDynamics);
	//e0 = molecularDynamics->potentialEnergy;
	e0 = 0; // Qi=0;
	fprintf(logfile,"E0 = The first potential energy in AU = %f\n",e0); 
	if( numberOfHeatSteps>0) molecularDynamics->qffModel->klass->printModesAndVelocities(molecularDynamics->qffModel,stdout);



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
			str = strdup_printf(("MD Heating: %0.2f fs, T = %0.2f K T(t) = %0.2f Kin = %0.4f Pot =  %0.4f Tot =  %0.4f"), 
					i*stepSize, 
					molecularDynamics->temperature, 
					molecularDynamics->kelvin, 
					molecularDynamics->kineticEnergy,
					molecularDynamics->potentialEnergy,
					molecularDynamics->totalEnergy
					);
			//redrawMolecule(&molecularDynamics->qffModel->molecule,str);
			fprintf(logfile,"%s\n",str);
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
	if(numberOfEquiSteps>0) molecularDynamics->qffModel->klass->printModesAndVelocities(molecularDynamics->qffModel,stdout);
	for (i = 0; i < numberOfEquiSteps; i++ )
	{
		molecularDynamics->temperature = currentTemp;
		applyOneStep(molecularDynamics,i);
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
			//redrawMolecule(&molecularDynamics->qffModel->molecule,str);
			fprintf(logfile,"%s\n",str);
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
	iSel = -1;
	if(numberOfGeometries>0)
	{
		iSel++;
		for(j=0;j<146;j++) fprintf(logfile,"-"); fprintf(logfile,"\n");
		if(str) free(str);
		str = strdup_printf("Geometry selected Potential energy(kcal) =  %0.4 f; File %s_%d.ici created ", molecularDynamics->potentialEnergy,suff,iSel);
		fprintf(logfile,"%s\n",str);
		fflush(logfile);
		for(j=0;j<146;j++) fprintf(logfile,"-"); fprintf(logfile,"\n");
		createNewInputFile(molecularDynamics->qffModel,inputFileName, iSel);
		// save cchemi file for iSel
	}
	esum  = 0;
	e2sum = 0;
	if(numberOfRunSteps>0) molecularDynamics->qffModel->klass->printModesAndVelocities(molecularDynamics->qffModel,stdout);
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
			//redrawMolecule(&molecularDynamics->qffModel->molecule,str);
			fprintf(logfile,"%s\n",str);
			updateNumber = 0;
			saveTrajectory(molecularDynamics, i+1);
		}
		if(numberOfGeometries>0 && (i+1)%stepSel==0 && (iSel+1)<numberOfGeometries)
                {
			iSel++;
			for(j=0;j<146;j++) fprintf(logfile,"-"); fprintf(logfile,"\n");
			if(str) free(str);
			str = strdup_printf("Geometry selected Potential energy(kcal) =  %0.4 f; File %s_%d.ici created ", molecularDynamics->potentialEnergy,suff,iSel);
			fprintf(logfile,"%s\n",str);
			fflush(logfile);
			for(j=0;j<146;j++) fprintf(logfile,"-"); fprintf(logfile,"\n");
			createNewInputFile(molecularDynamics->qffModel,inputFileName, iSel);
		// save cchemi file for iSel
                }

		saveProperties(molecularDynamics, n0+i+1, i+1," Running");
	}
	if(numberOfGeometries>0 && (iSel+1)<numberOfGeometries)
	{
		iSel++;
		for(j=0;j<146;j++) fprintf(logfile,"-"); fprintf(logfile,"\n");
		if(str) free(str);
		str = strdup_printf("Geometry selected Potential energy(kcal) =  %0.4 f; File %s_%d.ici created ", molecularDynamics->potentialEnergy,suff,iSel);
		fprintf(logfile,"%s\n",str);
		fflush(logfile);
		for(j=0;j<146;j++) fprintf(logfile,"-"); fprintf(logfile,"\n");
		createNewInputFile(molecularDynamics->qffModel,inputFileName, iSel);
	// save cchemi file for iSel
	}
	if(numberOfCoolSteps>0) molecularDynamics->qffModel->klass->printModesAndVelocities(molecularDynamics->qffModel,stdout);
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
			//redrawMolecule(&molecularDynamics->qffModel->molecule,str);
			fprintf(logfile,"%s\n",str);
			updateNumber = 0;
		}
		saveProperties(molecularDynamics, n0+i+1, i+1," Cooling");
	}
	molecularDynamics->qffModel->klass->calculateGradient(molecularDynamics->qffModel);
        gradientNorm = 0;
	for (i = 0; i < molecularDynamics->numberOfModes; i++)
                        gradientNorm += 
				molecularDynamics->qffModel->gradQ[i] * 
				molecularDynamics->qffModel->gradQ[i]; 

        gradientNorm = sqrt( gradientNorm );
	if(str) free(str);
	str = strdup_printf(("End of MD Simulation. Gradient = %f Ekin = %f (Kcal/mol) EPot =  %0.4f ETot =  %0.4f T(t) = %0.2f"),
			(double)gradientNorm*AUTOKCAL,
			molecularDynamics->kineticEnergy,
			molecularDynamics->potentialEnergy,
			molecularDynamics->totalEnergy,
			molecularDynamics->kelvin 
			); 
	//redrawMolecule(&molecularDynamics->qffModel->molecule,str);
	fprintf(logfile,"%s\n",str);
	if(numberOfGeometries>0)
        {
        	char* suff = getSuffixNameFile(inputFileName);
		for(i=0;i<146;i++) fprintf(logfile,"="); fprintf(logfile,"\n");
		if(str) free(str);
		str = strdup_printf("See %s_*.ici created files",suff);
		free(suff);
		fprintf(logfile,"%s\n",str);
		for(i=0;i<146;i++) fprintf(logfile,"="); fprintf(logfile,"\n");
	
	}
	free(str);
	free(suff);
	molecularDynamics->qffModel->klass->printModesAndVelocities(molecularDynamics->qffModel,stdout);
	if(molecularDynamics->fileTraj)fclose(molecularDynamics->fileTraj);
	if(molecularDynamics->fileProp)fclose(molecularDynamics->fileProp);
	freeQFFMD(molecularDynamics);
}
/*********************************************************************************/
static void initNH(QFFMD* molecularDynamics, double qNH)
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
static void initSD(QFFMD* molecularDynamics, double friction)
{
	if(friction<0) friction = 40;
	molecularDynamics->friction = friction/(fsInAU)/1000;

	molecularDynamics->positionFriction = NULL;
	molecularDynamics->velocityFriction = NULL;
	molecularDynamics->accelarationFriction = NULL;
	molecularDynamics->gamma = NULL;
	molecularDynamics->positionRandom = NULL;
	molecularDynamics->velocityRandom = NULL;

	if(molecularDynamics->integratorType != STOCHASTIC) return;

	molecularDynamics->positionFriction = malloc(molecularDynamics->numberOfModes *sizeof(double)); 
	molecularDynamics->velocityFriction = malloc(molecularDynamics->numberOfModes *sizeof(double)); 
	molecularDynamics->accelarationFriction = malloc(molecularDynamics->numberOfModes *sizeof(double)); 
	molecularDynamics->gamma = malloc(molecularDynamics->numberOfModes *sizeof(double)); 

	molecularDynamics->positionRandom = malloc(molecularDynamics->numberOfModes *sizeof(double)); 

	molecularDynamics->velocityRandom = malloc(molecularDynamics->numberOfModes *sizeof(double)); 

}
/*********************************************************************************/
static void compteThetaQTB(QFFMD* molecularDynamics)
{
	//double sigma=sqrt(2.*molecularDynamics->friction*KbInAU*molecularDynamics->temperature/molecularDynamics->h);
	double sigma=sqrt(6.*molecularDynamics->friction*KbInAU*molecularDynamics->temperature/molecularDynamics->h);
	int i,k;

	/* thetaProg = thetaPaper*sigma/sqrt(m) */
	for(i=0;i<molecularDynamics->numberOfModes;i++)
	{
		double mass = 1.0; // Mass-weigted normal coordinates
		//mass = molecularDynamics->qffModel->vData.effectiveMasses[i]*AMUTOAU; // to be use with convertTOAU2
		//mass = molecularDynamics->qffModel->vData.effectiveMasses[i]; 
		molecularDynamics->theta[i] = 0.0;
		for(k=0;k<2*molecularDynamics->Nf;k++)
			molecularDynamics->theta[i] += molecularDynamics->rnoise[i][2*molecularDynamics->Nf-1-k]*molecularDynamics->Ht[k];
		molecularDynamics->theta[i] *= sigma/sqrt(mass);
	}
}
/*********************************************************************************/
static void resetQTB(QFFMD* molecularDynamics)
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
	hbardwOverkT = 1.0/(molecularDynamics->Nf*molecularDynamics->h*KbInAU*T);
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
	
	for(i=0;i<molecularDynamics->numberOfModes;i++)
	for(k=0;k<2*molecularDynamics->Nf;k++)
		 molecularDynamics->rnoise[i][k] = normal();/* sqrt(h) is in sigma */

	compteThetaQTB(molecularDynamics);

}
/*********************************************************************************/
/* omegaMax in cm-1 */
static void initQTB(QFFMD* molecularDynamics, double omegaMax, double friction, int Nf)
{
/* Refs 
Jean-Louis Barrat , David Rodney
Portable implementation of a quantum thermal bath for molecular dynamics simulations
JOURNAL OF STATISTICAL PHYSICS 670, 144, (2011)
*/
	static double cmM1fsM1 = 2.99792458e-5;
	double Omegafs = omegaMax*cmM1fsM1;/* fs^-1 */ 
	int i;
	
	//printf("==========> dt(fs)\t\t= %f\n",molecularDynamics->dt/fsInAU);

	if(Nf<1) Nf = 50;

	molecularDynamics->Ht = NULL;
	molecularDynamics->theta = NULL;
	molecularDynamics->rnoise = NULL;
	molecularDynamics->Nf = 0;
	molecularDynamics->M = 0;

	if(molecularDynamics->integratorType != QTB) return;

	molecularDynamics->Nf = Nf;
	molecularDynamics->h = 1/Omegafs*(fsInAU);
	molecularDynamics->M = (int)(molecularDynamics->h/molecularDynamics->dt);
	if(molecularDynamics->M<1) molecularDynamics->M = 1;
	molecularDynamics->h = molecularDynamics->M *molecularDynamics->dt;
	omegaMax = 1.0/molecularDynamics->h*(fsInAU)/cmM1fsM1; /* cm-1 */
	if(friction<0) molecularDynamics->friction = (1.0/ molecularDynamics->h)/50;
	else molecularDynamics->friction = friction/1000.0/(fsInAU);

	molecularDynamics->Ht = malloc((2*Nf)*sizeof(double)); 
	molecularDynamics->theta = malloc(molecularDynamics->numberOfModes *sizeof(double)); 
	molecularDynamics->rnoise = malloc(molecularDynamics->numberOfModes *sizeof(double*)); 
	for(i=0;i<molecularDynamics->numberOfModes;i++)
		molecularDynamics->rnoise[i] = malloc((2*Nf)*sizeof(double)); 
	

	printf("\n");
	printf("*************** QTB Parameters ******************************************************************\n");
	printf("Nf\t\t= %d\n",molecularDynamics->Nf);
	printf("M\t\t= %d\n",molecularDynamics->M);
	printf("dt(fs)\t\t= %f\n",molecularDynamics->dt/(fsInAU));
	printf("h(fs)\t\t= %f\n",molecularDynamics->h/(fsInAU));
	printf("gamma(ps^-1)\t= %f\n",molecularDynamics->friction*fsInAU*1000);
	printf("omegaMax(cm^-1)\t= %f\n",omegaMax);
	printf("*************************************************************************************************\n");
	printf("\n");

	resetQTB(molecularDynamics);

}
/*********************************************************************************/
static void updateQTB(QFFMD* molecularDynamics)
{
	int i,k;

	if(molecularDynamics->temperature<=0) return;
	compteThetaQTB(molecularDynamics);
	/* shift rnoise */
	for(i=0;i<molecularDynamics->numberOfModes;i++)
	for(k=0;k<2*molecularDynamics->Nf-1;k++)
		 molecularDynamics->rnoise[i][k] = molecularDynamics->rnoise[i][k+1];

	/* add one value to the end */
	for(i=0;i<molecularDynamics->numberOfModes;i++)
		 molecularDynamics->rnoise[i][2*molecularDynamics->Nf-1] = normal(); /* sqrt(h) in sigma */

}
/*********************************************************************************/
static void resetLangevin(QFFMD* molecularDynamics)
{
	if(molecularDynamics->integratorType != LANGEVIN) return;
	updateLangevin(molecularDynamics);
}
/*********************************************************************************/
/* omegaMax in cm-1 */
static void initLangevin(QFFMD* molecularDynamics, double friction)
{
	if(molecularDynamics->integratorType != LANGEVIN) return;

	if(friction<0) friction = 40;
	molecularDynamics->friction = friction/1000/(fsInAU);
	molecularDynamics->theta = malloc(3*molecularDynamics->numberOfModes *sizeof(double)); 
	
	printf("\n");
	printf("*************** Langevin Parameters ******************************************************************\n");
	printf("dt(fs)\t\t= %f\n",molecularDynamics->dt/(fsInAU));
	printf("gamma(ps^-1)\t= %f\n",molecularDynamics->friction*fsInAU*1000);
	printf("*************************************************************************************************\n");
	printf("\n");

	resetLangevin(molecularDynamics);

}
/*********************************************************************************/
static void updateLangevin(QFFMD* molecularDynamics)
{
	int i;
	double sigma;
	/* update theta */
	/* thetaProg = thetaPaper*sigma/sqrt(m) */
	if(molecularDynamics->integratorType != LANGEVIN) return;
	sigma=sqrt(6.*molecularDynamics->friction*KbInAU*molecularDynamics->temperature/molecularDynamics->dt);
	for(i=0;i<molecularDynamics->numberOfModes;i++)
	{
		double mass = 1.0; // Mass-weigted normal coordinates
		//mass = molecularDynamics->qffModel->vData.effectiveMasses[i]*AMUTOAU;
		//mass = molecularDynamics->qffModel->vData.effectiveMasses[i];
		molecularDynamics->theta[i] = normal();
		molecularDynamics->theta[i] *= sigma/sqrt(mass);
	}
}
/*********************************************************************************/
static void initMD(QFFMD* molecularDynamics, double temperature, double stepSize, MDIntegratorType integratorType, MDThermostatType thermostat, double friction, double omegaMax, int Nf, double collide, double qNH, char* fileNameTraj, char* fileNameProp, int numberOfRunSteps, int index)
{
	int i;
	double dt = stepSize * fsInAU;

	//printf("----------> stepSize(fs)\t\t= %f\n",stepSize);
	//printf("----------> dt(fs)\t\t= %f\n",dt/fsInAU);

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

	molecularDynamics->nFree = molecularDynamics->numberOfModes;
	for(i=0;i<molecularDynamics->numberOfModes;i++) if(!molecularDynamics->qffModel->variable[i]) molecularDynamics->nFree--;

	molecularDynamics->a = malloc(molecularDynamics->numberOfModes *sizeof(double)); 

	molecularDynamics->aold = NULL;
	if(molecularDynamics->integratorType==BEEMAN)
	{
		molecularDynamics->aold = malloc(molecularDynamics->numberOfModes *sizeof(double)); 
	}
	molecularDynamics->coordinatesOld = NULL;
	molecularDynamics->moved = NULL;
	molecularDynamics->update = NULL;
	if(molecularDynamics->qffModel->molecule.constraints!=NOCONSTRAINTS)
	{
		molecularDynamics->coordinatesOld = malloc(molecularDynamics->numberOfModes *sizeof(double)); 
		molecularDynamics->moved = malloc(molecularDynamics->numberOfModes *sizeof(boolean)); 
		molecularDynamics->update = malloc(molecularDynamics->numberOfModes *sizeof(boolean)); 

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


	molecularDynamics->qffModel->klass->calculateGradient(molecularDynamics->qffModel);
	for ( i = 0; i < molecularDynamics->numberOfModes; i++)
	{
		double mass = 1.0; // Mass-weigted normal coordinates
		//mass = molecularDynamics->qffModel->vData.effectiveMasses[i]*AMUTOAU;
		//mass = molecularDynamics->qffModel->vData.effectiveMasses[i];
		molecularDynamics->a[i] = -molecularDynamics->qffModel->gradQ[i]/mass;
		if(molecularDynamics->aold) molecularDynamics->aold[i]  = molecularDynamics->a[i];
	}
	if(molecularDynamics->qffModel->klass->setMaxwellVelocitiesIfNull(molecularDynamics->qffModel, temperature)) rescaleVelocities(molecularDynamics);
#ifdef DEBUG
	printf("nfree =%d\n",molecularDynamics->nFree);
#endif
}
/*********************************************************************************/
static void rescaleVelocities(QFFMD* molecularDynamics)
{
	/* berendsen(molecularDynamics);*/
	scaleV(molecularDynamics);
	resetQTB(molecularDynamics);
	resetLangevin(molecularDynamics);
}
/*********************************************************************************/
static void scaleV(QFFMD* molecularDynamics)
{
	int i;
	double ekin = 0;
	double kelvin = 0;
	int nfree = molecularDynamics->nFree;
	double scale = 1.0;
	if(molecularDynamics->temperature<=0) return;
	if(nfree<1) return;
	for ( i = 0; i < molecularDynamics->numberOfModes; i++)
	{
		double mass = 1.0; // Mass-weigted normal coordinates
		//mass = molecularDynamics->qffModel->vData.effectiveMasses[i]*AMUTOAU;
		//mass = molecularDynamics->qffModel->vData.effectiveMasses[i];
		ekin += molecularDynamics->qffModel->velocity[i]*molecularDynamics->qffModel->velocity[i]*mass;
	}
	/*
	ekin /= 2;
	kelvin = 2* ekin / ( nfree * KbInAU);
	*/
	kelvin = ekin / ( nfree * KbInAU);
	scale = sqrt(molecularDynamics->temperature/kelvin);
#ifdef DEBUG
	printf("temp = %f kelvin = %f scale = %f\n",molecularDynamics->temperature, kelvin, scale);
#endif
	for ( i = 0; i < molecularDynamics->numberOfModes; i++)
		if(molecularDynamics->qffModel->variable[i])
			molecularDynamics->qffModel->velocity[i] *= scale;
}
/*********************************************************************************/
static void berendsen(QFFMD* molecularDynamics)
{
	int i;
	double ekin = 0;
	double kelvin = 0;
	int nfree = molecularDynamics->nFree;
	double scale = 1.0;
	double dt = molecularDynamics->dt;
	double tautemp = 1.0/(molecularDynamics->collide)*1000*fsInAU;
	if(molecularDynamics->temperature<=0) return;
	if(nfree<1) return;
	for ( i = 0; i < molecularDynamics->numberOfModes; i++)
	{
		double mass = 1.0; // Mass-weigted normal coordinates
		//mass = molecularDynamics->qffModel->vData.effectiveMasses[i]*AMUTOAU;
		//mass = molecularDynamics->qffModel->vData.effectiveMasses[i];
		ekin += molecularDynamics->qffModel->velocity[i]*molecularDynamics->qffModel->velocity[i]*mass;
	}
	/*
	ekin /= 2;
	kelvin = 2* ekin / ( nfree * KbInAU);
	*/
	kelvin = ekin / ( nfree * KbInAU);
	/* if(tautemp>dt) tautemp = dt;*/
	scale = sqrt(1.0 + (dt/tautemp)*(molecularDynamics->temperature/kelvin-1.0));
#ifdef DEBUG
	printf("temp = %f kelvin = %f scale = %f\n",molecularDynamics->temperature, kelvin, scale);
#endif
	for ( i = 0; i < molecularDynamics->numberOfModes; i++)
		if(molecularDynamics->qffModel->variable[i])
			molecularDynamics->qffModel->velocity[i] *= scale;
}
/*********************************************************************************/
static void andersen(QFFMD* molecularDynamics)
{
	int i;
	double tau = 1.0/molecularDynamics->collide*1000*fsInAU; /* in fs */
	double rate;
	if(molecularDynamics->temperature<=0) return;
	if(molecularDynamics->numberOfModes<1) return;

	rate = molecularDynamics->dt / tau;
	rate /= pow(molecularDynamics->nvariables,2.0/3.0);

	for ( i = 0; i < molecularDynamics->numberOfModes; i++)
	{
		double trial = drandom();
		double mass = 1.0;
		//mass = molecularDynamics->qffModel->vData.effectiveMasses[i]*AMUTOAU;
		//mass = molecularDynamics->qffModel->vData.effectiveMasses[i];
		if(trial<rate)
		{
			double speed = sqrt(KbInAU* molecularDynamics->temperature/mass);
                	double pnorm = normal();
			molecularDynamics->qffModel->velocity[i] = pnorm*speed;
		}
	}
}
/*********************************************************************************/
static void bussi(QFFMD* molecularDynamics)
{
	int nfree = molecularDynamics->nFree;
	double scale = 1.0;
	double dt = molecularDynamics->dt;
	double tautemp = 1.0/(molecularDynamics->collide)*1000*fsInAU;
        double c = exp(-dt/tautemp);
	double ekin = getEKin(molecularDynamics);
	double kelvin = 2*ekin / ( nfree * KbInAU);
	double d = (1.0-c) * (molecularDynamics->temperature/kelvin) / (nfree);
	double r = normal ();
	double si = 0.0;
	double s = 0.0;
	int i;
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
	for ( i = 0; i < molecularDynamics->numberOfModes; i++)
		if(molecularDynamics->qffModel->variable[i])
			molecularDynamics->qffModel->velocity[i] *= scale;
}
/*********************************************************************************/
static void nose_hoover(QFFMD* molecularDynamics)
{
	int nfree = molecularDynamics->nFree;
	double scale = 1.0;
	double ekin = getEKin(molecularDynamics);
	double kT = KbInAU* molecularDynamics->temperature;
	int i;
	if(molecularDynamics->temperature<=0) return;
	if(nfree<1) return;
	molecularDynamics->gNH[1] = (molecularDynamics->qNH[0]*molecularDynamics->vNH[0]*molecularDynamics->vNH[0]-kT) / molecularDynamics->qNH[1];
	//printf("gNH = %f\n",molecularDynamics->gNH[1]);
	molecularDynamics->vNH[1] = molecularDynamics->vNH[1] + molecularDynamics->gNH[1]*molecularDynamics->dt_4;
	molecularDynamics->vNH[0] = molecularDynamics->vNH[0] * exp(-molecularDynamics->vNH[1]*molecularDynamics->dt_8);
	molecularDynamics->gNH[0] = (2.0*ekin-molecularDynamics->nFree*kT) / molecularDynamics->qNH[0];
	molecularDynamics->vNH[0] = molecularDynamics->vNH[0] + molecularDynamics->gNH[0]*molecularDynamics->dt_4;
	molecularDynamics->vNH[0] = molecularDynamics->vNH[0] * exp(-molecularDynamics->vNH[1]*molecularDynamics->dt_8);
	molecularDynamics->xNH[0] = molecularDynamics->xNH[0] + molecularDynamics->vNH[0]*molecularDynamics->dt_2;
	molecularDynamics->xNH[1] = molecularDynamics->xNH[1] + molecularDynamics->vNH[1]*molecularDynamics->dt_2;
	//printf("vnH0 = %f\n",molecularDynamics->vNH[0]);
	scale = exp(-molecularDynamics->vNH[0]*molecularDynamics->dt_2);
	//printf("scale = %f\n",scale);
	for ( i = 0; i < molecularDynamics->numberOfModes; i++)
			molecularDynamics->qffModel->velocity[i] *= scale;
	ekin = ekin * scale * scale;
	molecularDynamics->vNH[0] = molecularDynamics->vNH[0] * exp(-molecularDynamics->vNH[1]*molecularDynamics->dt_8);
	molecularDynamics->gNH[0] = (2.0*ekin-nfree*kT) /  molecularDynamics->qNH[0];
	molecularDynamics->vNH[0] = molecularDynamics->vNH[0] + molecularDynamics->gNH[0]*molecularDynamics->dt_4;
	molecularDynamics->vNH[0] = molecularDynamics->vNH[0] * exp(-molecularDynamics->vNH[1]*molecularDynamics->dt_8);
	molecularDynamics->gNH[1] = ( molecularDynamics->qNH[0]*molecularDynamics->vNH[0]*molecularDynamics->vNH[0]-kT) /  molecularDynamics->qNH[1];
	molecularDynamics->vNH[1] = molecularDynamics->vNH[1] + molecularDynamics->gNH[1]*molecularDynamics->dt_4;

}
/*********************************************************************************/
static void newAccelaration(QFFMD* molecularDynamics)
{
	int i;
	molecularDynamics->qffModel->klass->calculateGradient(molecularDynamics->qffModel);
	for ( i = 0; i < molecularDynamics->numberOfModes; i++)
	{
		double mass =  1.0;
		//mass = molecularDynamics->qffModel->vData.effectiveMasses[i]*AMUTOAU;
		//mass = molecularDynamics->qffModel->vData.effectiveMasses[i];
		if(molecularDynamics->aold)
				molecularDynamics->aold[i]  = molecularDynamics->a[i];

			molecularDynamics->a[i] = -molecularDynamics->qffModel->gradQ[i]/mass;
	}
}
/*********************************************************************************/
static void computeEnergies(QFFMD* molecularDynamics)
{
	molecularDynamics->kineticEnergy = getEKin(molecularDynamics);
	molecularDynamics->potentialEnergy = molecularDynamics->qffModel->molecule.potentialEnergy;
	molecularDynamics->totalEnergy = molecularDynamics->kineticEnergy + molecularDynamics->potentialEnergy;
	molecularDynamics->kelvin = getKelvin(molecularDynamics);

/* convert in kcal/mol*/
	molecularDynamics->kineticEnergy *= AUTOKCAL;
	molecularDynamics->potentialEnergy *= AUTOKCAL; 
	molecularDynamics->totalEnergy *= AUTOKCAL;
}
/*********************************************************************************/
static void applyThermostat(QFFMD* molecularDynamics)
{
	if(molecularDynamics->integratorType == STOCHASTIC) return;
	if(molecularDynamics->integratorType == QTB) return;
	if(molecularDynamics->integratorType == LANGEVIN) return;
	if(molecularDynamics->thermostat == ANDERSEN) andersen(molecularDynamics);
	if(molecularDynamics->thermostat == BERENDSEN) berendsen(molecularDynamics);
	if(molecularDynamics->thermostat == BUSSI) bussi(molecularDynamics);
}
/*********************************************************************************/
static void applyOneStep(QFFMD* molecularDynamics, int iStep)
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

}
/*********************************************************************************/
static void applyConstraintsFirstPortion(QFFMD* quantumMechanicsMD)
{
	int i;
	QFFModel* qffModel = quantumMechanicsMD->qffModel;

	if(qffModel->molecule.constraints==NOCONSTRAINTS) return;
	for (i = 0; i <  quantumMechanicsMD->numberOfModes; i++)
	if(!quantumMechanicsMD->qffModel->variable[i])
	{
		qffModel->velocity[i] = 0.0;
		qffModel->Q[i] =  quantumMechanicsMD->coordinatesOld[i];
	}

}
/*********************************************************************************/
static void applyConstraintsSecondPortion(QFFMD* quantumMechanicsMD)
{
	int i;
	QFFModel* qffModel = quantumMechanicsMD->qffModel;

	if(qffModel->molecule.constraints==NOCONSTRAINTS) return;
	for (i = 0; i <  quantumMechanicsMD->numberOfModes; i++)
	if(!quantumMechanicsMD->qffModel->variable[i])
	{
		qffModel->velocity[i] = 0.0;
		qffModel->Q[i] =  quantumMechanicsMD->coordinatesOld[i];
	}
}
/*********************************************************************************/
static void applyVerlet(QFFMD* molecularDynamics)
{
	int i;

	if(molecularDynamics->qffModel->molecule.constraints!=NOCONSTRAINTS)
	for (i = 0; i < molecularDynamics->numberOfModes; i++)
				molecularDynamics->coordinatesOld[i]= molecularDynamics->qffModel->Q[i];

	for (i = 0; i < molecularDynamics->numberOfModes; i++)
	{
		if(!molecularDynamics->qffModel->variable[i]) continue;

		{
			molecularDynamics->qffModel->Q[i] += 
				molecularDynamics->qffModel->velocity[i] * molecularDynamics->dt +
				molecularDynamics->a[i]*molecularDynamics->dt2_2;	
		}
			molecularDynamics->qffModel->velocity[i] += molecularDynamics->a[i] * molecularDynamics->dt_2;
	}

	applyConstraintsFirstPortion(molecularDynamics);
	if(molecularDynamics->thermostat==NOSEHOOVER) nose_hoover(molecularDynamics);

	newAccelaration(molecularDynamics);

	for (i = 0; i < molecularDynamics->numberOfModes; i++)
		if(molecularDynamics->qffModel->variable[i])
			molecularDynamics->qffModel->velocity[i] += molecularDynamics->a[i] * molecularDynamics->dt_2;
	applyConstraintsSecondPortion(molecularDynamics);
}
/*********************************************************************************/
static void applyBeeman(QFFMD* molecularDynamics)
{
	int i;
	double terms;
	if(molecularDynamics->qffModel->molecule.constraints!=NOCONSTRAINTS)
	for (i = 0; i < molecularDynamics->numberOfModes; i++)
				molecularDynamics->coordinatesOld[i]= molecularDynamics->qffModel->Q[i];

	for (i = 0; i < molecularDynamics->numberOfModes; i++)
	{
		if(!molecularDynamics->qffModel->variable[i]) continue;

			terms = 5.0*molecularDynamics->a[i]-molecularDynamics->aold[i];

		{
			molecularDynamics->qffModel->Q[i] += 
				molecularDynamics->qffModel->velocity[i] * molecularDynamics->dt +
				terms*molecularDynamics->dt2_8;	
		}
			molecularDynamics->qffModel->velocity[i] += terms * molecularDynamics->dt_8;
	}

	applyConstraintsFirstPortion(molecularDynamics);
	if(molecularDynamics->thermostat==NOSEHOOVER) nose_hoover(molecularDynamics);

	newAccelaration(molecularDynamics);

	for (i = 0; i < molecularDynamics->numberOfModes; i++)
		if(molecularDynamics->qffModel->variable[i])
			molecularDynamics->qffModel->velocity[i] += (3.0*molecularDynamics->a[i]+molecularDynamics->aold[i]) * molecularDynamics->dt_8;
	applyConstraintsSecondPortion(molecularDynamics);
}
/*********************************************************************************/
static void newProperties(QFFMD* molecularDynamics, char* comments)
{
	if( molecularDynamics->fileProp == NULL) return;
	fprintf(molecularDynamics->fileProp,"Time0(fs)\tTime(fs)\tTotal Energy(Kcal/mol)\tPotential Energy(Hartree)\tKinetic Energy(Hartree)\tT(t) (K)\tTaver(K)\tsigma(T)(K)\tIndex\tmuX\tmuY\tmuZ");
	if(molecularDynamics->thermostat==NOSEHOOVER) fprintf(molecularDynamics->fileProp,"\tEtot+Etherm");
	if(comments) fprintf(molecularDynamics->fileProp,"%s\n", comments);
	else fprintf(molecularDynamics->fileProp,"\n");
}
/*********************************************************************************/
static void saveProperties(QFFMD* molecularDynamics, int iStep0, int iStep, char* comments)
{
	double dt = molecularDynamics->dt/(fsInAU);
	static double Ttot = 0;
	static double T2tot = 0;
	double Taver = 0;
	double T2aver = 0;
	double totalEnergy =  molecularDynamics->totalEnergy;

	if( molecularDynamics->thermostat==NOSEHOOVER)
	{
		int i;
		double kT = KbInAU* molecularDynamics->temperature;
		double e = molecularDynamics->vNH[0]*molecularDynamics->vNH[0]* molecularDynamics->qNH[0]/2 + (molecularDynamics->nFree)*kT* molecularDynamics->xNH[0];
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
			molecularDynamics->qffModel->molecule.dipole[0]*AUTODEB, 
			molecularDynamics->qffModel->molecule.dipole[1]*AUTODEB, 
			molecularDynamics->qffModel->molecule.dipole[2]*AUTODEB);
	if( molecularDynamics->thermostat==NOSEHOOVER) fprintf(molecularDynamics->fileProp,"%f\t",totalEnergy);
	if(comments) fprintf(molecularDynamics->fileProp,"%s\n", comments);
	else fprintf(molecularDynamics->fileProp,"\n");
}
/*********************************************************************************/
static void saveTrajectory(QFFMD* molecularDynamics, int iStep)
{
	double dt = molecularDynamics->dt/(fsInAU);
	if( molecularDynamics->fileTraj == NULL) return;

	fprintf(molecularDynamics->fileTraj,"%d %f %f %f %f %f %f %f nModes, time(fs) TotalEnery(Hartree) Kinetic Potential, Dipole in AU\n", 
			molecularDynamics->numberOfModes,
			 (iStep)*dt, 
			molecularDynamics->totalEnergy,
			molecularDynamics->kineticEnergy,
			molecularDynamics->potentialEnergy,
			molecularDynamics->qffModel->molecule.dipole[0],
			molecularDynamics->qffModel->molecule.dipole[1],
			molecularDynamics->qffModel->molecule.dipole[2]
			 );
}

/**********************************************************************/
void	freeQFFMD(QFFMD* molecularDynamics)
{

	molecularDynamics->qffModel = NULL;
	molecularDynamics->numberOfModes = 0;
	molecularDynamics->updateFrequency = 0;
	if(molecularDynamics->a)
	{
		free(molecularDynamics->a);
	}
	if(molecularDynamics->aold)
	{
		free(molecularDynamics->aold);
	}
	if(molecularDynamics->coordinatesOld)
	{
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
		int i;
		for(i=0;i<molecularDynamics->numberOfModes;i++)
			if(molecularDynamics->rnoise[i]) free(molecularDynamics->rnoise[i]);
		free(molecularDynamics->rnoise);
	}
}
/********************************************************************************/
static double getEKin(QFFMD* molecularDynamics)
{
	double ekin = 0;
	int i;
	for ( i = 0; i < molecularDynamics->numberOfModes; i++)
	{
		double mass = 1.0; // Mass-weigted normal coordinates
		//mass = molecularDynamics->qffModel->vData.effectiveMasses[i]*AMUTOAU;
		//mass = molecularDynamics->qffModel->vData.effectiveMasses[i];
		ekin += molecularDynamics->qffModel->velocity[i]*molecularDynamics->qffModel->velocity[i]*mass;
	}
	return ekin/2;
}
/********************************************************************************/
static double getKelvin(QFFMD* molecularDynamics)
{
	int nfree = molecularDynamics->nFree;
	/* printf("nfree = %d\n",nfree);*/
	if(nfree<1) return 0;
	return 2*getEKin(molecularDynamics) / ( nfree * KbInAU);
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
static void getsFrictionalAndRandomForce(QFFMD* molecularDynamics)
{
	double* gamma = molecularDynamics->gamma;
	double* positionFriction = molecularDynamics->positionFriction;
	double* velocityFriction = molecularDynamics->velocityFriction;
	double* accelarationFriction = molecularDynamics->accelarationFriction;
	double* positionRandom = molecularDynamics->positionRandom;
	double* velocityRandom = molecularDynamics->velocityRandom;
	double dt = molecularDynamics->dt;
	
	int n = molecularDynamics->numberOfModes;

	int i;
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
		double mass = 1.0; // Mass-weigted normal coordinates
		//mass = molecularDynamics->qffModel->vData.effectiveMasses[i]*AMUTOAU;
		//mass = molecularDynamics->qffModel->vData.effectiveMasses[i];
		gdt = gamma[i] * dt;
		/* printf("gdt = %f\n",gdt);*/
		if (gdt <= 0.0)
		{
               		positionFriction[i] = 1.0;
			velocityFriction[i] = dt;
			accelarationFriction[i] = 0.5 * dt * dt;
			{
                  		positionRandom[i] = 0.0;
                  		velocityRandom[i] = 0.0;
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
		
               		ktm = KbInAU * molecularDynamics->temperature / mass;
               		psig = sqrt(ktm*pterm) / gamma[i];
               		vsig = sqrt(ktm*vterm);
               		rhoc = sqrt(1.0 - rho*rho);
			{
                		pnorm = normal();
             			vnorm = normal ();
				positionRandom[i] = psig * pnorm;
                  		velocityRandom[i] = vsig * (rho*pnorm+rhoc*vnorm);
			}
		}
	}
}
/*********************************************************************************/
static void applyStochastic(QFFMD* molecularDynamics)
{
	double* positionFriction = molecularDynamics->positionFriction;
	double* velocityFriction = molecularDynamics->velocityFriction;
	double* accelarationFriction = molecularDynamics->accelarationFriction;
	double* positionRandom = molecularDynamics->positionRandom;
	double* velocityRandom = molecularDynamics->velocityRandom;
	double*a = molecularDynamics->a;
	
	int n = molecularDynamics->numberOfModes;
	int i;
	double * Q = molecularDynamics->qffModel->Q;

	getsFrictionalAndRandomForce(molecularDynamics);

	if(molecularDynamics->qffModel->molecule.constraints!=NOCONSTRAINTS)
	for (i = 0; i < n; i++)
				molecularDynamics->coordinatesOld[i]= molecularDynamics->qffModel->Q[i];

	for(i=0;i<n;i++)
	{
		if(!molecularDynamics->qffModel->variable[i]) continue;
		Q[i] += molecularDynamics->qffModel->velocity[i]*velocityFriction[i] + a[i]*accelarationFriction[i] + positionRandom[i];
		molecularDynamics->qffModel->velocity[i] = molecularDynamics->qffModel->velocity[i]*positionFriction[i] + 0.5*a[i]*velocityFriction[i];
	}

	applyConstraintsFirstPortion(molecularDynamics);

	newAccelaration(molecularDynamics);

	for (i = 0; i < n; i++)
		if(molecularDynamics->qffModel->molecule.atoms[i].variable)
			molecularDynamics->qffModel->velocity[i] += 0.5*a[i]*velocityFriction[i] + velocityRandom[i];
	applyConstraintsSecondPortion(molecularDynamics);
	computeEnergies(molecularDynamics);
}
/*********************************************************************************/
static void applyQTB(QFFMD* molecularDynamics)
{
	
	int n = molecularDynamics->numberOfModes;
	int i;
	double gp = 1/(1+molecularDynamics->friction*molecularDynamics->dt_2);
	double gm = (1-molecularDynamics->friction*molecularDynamics->dt_2)*gp;

	/* printf("gm = %f gp =%f\n",gm,gp);*/
	if(molecularDynamics->qffModel->molecule.constraints!=NOCONSTRAINTS)
	for (i = 0; i < n; i++)
				molecularDynamics->coordinatesOld[i]= molecularDynamics->qffModel->Q[i];

	for (i = 0; i < n; i++)
		if(molecularDynamics->qffModel->variable[i])
		{
			molecularDynamics->qffModel->Q[i] += 
			molecularDynamics->qffModel->velocity[i]*molecularDynamics->dt +
			(molecularDynamics->a[i]+molecularDynamics->theta[i]-molecularDynamics->qffModel->velocity[i]*molecularDynamics->friction)*molecularDynamics->dt2_2;	
		}
	for (i = 0; i < n; i++)
		if(molecularDynamics->qffModel->variable[i])
			molecularDynamics->qffModel->velocity[i] = 
			gm*molecularDynamics->qffModel->velocity[i] +
			gp*(molecularDynamics->a[i]+molecularDynamics->theta[i])*molecularDynamics->dt_2;
            
	applyConstraintsFirstPortion(molecularDynamics);

	newAccelaration(molecularDynamics);

	for (i = 0; i < n; i++)
		if(molecularDynamics->qffModel->variable[i])
			molecularDynamics->qffModel->velocity[i] += gp*(molecularDynamics->a[i]+molecularDynamics->theta[i])*molecularDynamics->dt_2;

	applyConstraintsSecondPortion(molecularDynamics);
	computeEnergies(molecularDynamics);
}
/*********************************************************************************/
static void applyLangevin(QFFMD* molecularDynamics)
{
	
	int n = molecularDynamics->numberOfModes;
	int i;
	double gp = 1/(1+molecularDynamics->friction*molecularDynamics->dt_2);
	double gm = (1-molecularDynamics->friction*molecularDynamics->dt_2)*gp;

	/* printf("gm = %f gp =%f\n",gm,gp);*/
	if(molecularDynamics->qffModel->molecule.constraints!=NOCONSTRAINTS)
	for (i = 0; i < n; i++)
				molecularDynamics->coordinatesOld[i]= molecularDynamics->qffModel->Q[i];

	for (i = 0; i < n; i++)
		if(molecularDynamics->qffModel->variable[i])
		{
			molecularDynamics->qffModel->Q[i] += 
			molecularDynamics->qffModel->velocity[i]*molecularDynamics->dt +
			(molecularDynamics->a[i]+molecularDynamics->theta[i]-molecularDynamics->qffModel->velocity[i]*molecularDynamics->friction)*molecularDynamics->dt2_2;	
		}
	for (i = 0; i < n; i++)
		if(molecularDynamics->qffModel->variable[i])
			molecularDynamics->qffModel->velocity[i] = 
			gm*molecularDynamics->qffModel->velocity[i] +
			gp*(molecularDynamics->a[i]+molecularDynamics->theta[i])*molecularDynamics->dt_2;
            
	applyConstraintsFirstPortion(molecularDynamics);

	newAccelaration(molecularDynamics);

	for (i = 0; i < n; i++)
		if(molecularDynamics->qffModel->variable[i])
			molecularDynamics->qffModel->velocity[i] += gp*(molecularDynamics->a[i]+molecularDynamics->theta[i])*molecularDynamics->dt_2;

	applyConstraintsSecondPortion(molecularDynamics);
	computeEnergies(molecularDynamics);
}
/*********************************************************************************/
void QFFMDDlg(char* inputFileName)
{
	QFFModel qffModel; 
	QFFMD molecularDynamics;
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
	MDIntegratorType integrator = VERLET;
	char* fileNameTraj = NULL;
	char* fileNameProp = NULL;
	double friction=-1;
	double omegaMax = 4000;
	int Nf = 50;
	double collide = 20;
	double qNH = 20;
	MDThermostatType thermostat = NONE;
	char* dirName = NULL;
	FILE* file = fopen(inputFileName,"rb");
	int numberOfGeometries = 0;

	setMDOptions(file, &updateFrequency, 
		&heatTime, &equiTime, &runTime, &coolTime,
		&heatTemp, &runTemp, &equiTemp, &coolTemp, &stepSize, 
		&integrator, &thermostat, &friction, &omegaMax, &Nf, &collide,&qNH);


	{
		char* suff = getSuffixNameFile(inputFileName);
		dirName = strdup(getenv("PWD"));
		fileNameTraj = strdup_printf("%s%s",suff, "Traj.gab");
		fileNameProp = strdup_printf("%s%s",suff, "Prop.txt");
		free(suff);
	}
	readOneInt(file,"numberOfGeometries",&numberOfGeometries);
	fclose(file);


	qffModel = newQFFModel();
	qffModel.klass->readData(&qffModel, inputFileName, inputFileName);
	//qffModel.klass->convertToAU2(&qffModel);
	qffModel.klass->convertToAU(&qffModel);
	qffModel.klass->computeQFFParameters(&qffModel);
	qffModel.molecule.constraints = NOCONSTRAINTS;

	runQFFMD(&molecularDynamics, &qffModel,
		updateFrequency, heatTime, equiTime, runTime, coolTime, heatTemp, equiTemp, runTemp, coolTemp, stepSize, 
		integrator, thermostat, friction, omegaMax, Nf, collide, qNH, numberOfGeometries, fileNameTraj, fileNameProp,inputFileName);

	qffModel.klass->free(&qffModel);
	free(dirName);
}
