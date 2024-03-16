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

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include "Paths.h"
#include "LangevinThermostat.h"
#include "NoseHooverThermostat.h"
#include "QTBThermostat.h"
#include "NormalModeTrans.h"
#include "StagingTrans.h"
#include "WithoutTrans.h"

static void calculateGradient(Paths* paths);
void initPathsTemperature(Paths* paths, double newTemperature);
void resetPathsTemperature(Paths* paths, double newTemperature);
static void initThermostat(Paths* paths);
static void initTransformation(Paths* paths);
static void computeEnergies(Paths* paths);
static void initIndex(Paths* paths);
static void printPaths(Paths* paths,FILE* file);
static void getGradVPlusHarmonic(Paths* paths, int iBead, int iAtom, double G[]);
static void applyNoThermostat(Paths* paths);
static void initAllPointers(Paths* paths);
static void initVelocities(Paths* paths);
static void removeTranslation(Paths* paths);
static void removeRotation(Paths* paths);
static void removeTranslationAndRotation(Paths* paths);
/*********************************************************************************************************************/
Paths newPaths(
		ForceField* forceField,
		QuantumMechanicsModel* qmModel, 
		int updateFrequency, 
		double temperature, 
		double stepSize,
		PIMDThermostatType thermostat,
		PIMDTransformationType transformation,
		double friction,
		double omegaMax,
		int Nf,
		int nBeads,
		int nNH,
		int nNHSteps,
		int nSY,
		char* fileNameTraj,
		char* fileNameProp
)
{
	int i;
	int nAtoms = 0;
	double dt = stepSize * fsInAKMA;
	Paths opaths;
	Paths* paths = &opaths;
	paths->nBeads = nBeads;
	paths->kT = Kb*temperature;
	paths->temperature = temperature;
	initAllPointers(paths);
	if(forceField) 
	{
		nAtoms = forceField->molecule.nAtoms;
		paths->forceFields = malloc(paths->nBeads*sizeof(ForceField));
		for(i = 0;i<paths->nBeads; i++) 
			paths->forceFields[i] = copyForceField(forceField);
		paths->molecules = malloc(paths->nBeads*sizeof(Molecule*));
		for(i = 0;i<paths->nBeads; i++) 
			paths->molecules[i] = &paths->forceFields[i].molecule;
		paths->qmModels = NULL;
		

	}
	else if(qmModel) 
	{
		nAtoms = qmModel->molecule.nAtoms;
		paths->qmModels = malloc(paths->nBeads*sizeof(QuantumMechanicsModel));
		for(i = 0;i<paths->nBeads; i++) 
			paths->qmModels[i] = qmModel->klass->copy(qmModel);
		paths->molecules = malloc(paths->nBeads*sizeof(Molecule*));
		for(i = 0;i<paths->nBeads; i++) 
			paths->molecules[i] = &paths->qmModels[i].molecule;
		paths->forceFields = NULL;
	}
	else exit(1);
	paths->moleculeCentroid = paths->molecules[0]->klass->copy(paths->molecules[0]);
	initVelocities(paths);
	initIndex(paths);

	paths->friction = friction/fsInAKMA/1000;
	paths->nNH = nNH;
	paths->nNHSteps = nNHSteps;
	paths->nSY = nSY;
	paths->thermostat = thermostat;
	paths->transformation = transformation;
	paths->updateFrequency = updateFrequency;
	paths->dt = dt;
	paths->omegaMax= omegaMax;
	paths->Nf= Nf;

	initPathsTemperature(paths, paths->temperature);

	paths->oneOverNumberOfParticules = 1.0/nAtoms;
	paths->oneOvernBeads = 1.0/paths->nBeads;
	paths->oneOverNumberOfParticulesBeads= 1.0/(1.0*nAtoms*paths->nBeads);

	if(fileNameTraj) paths->fileTraj = fopen(fileNameTraj, "w");
	if(fileNameProp) paths->fileProp = fopen(fileNameProp, "w");

	paths->potentialEnergy = 0;
	paths->kineticEnergy = 0;
	paths->totalEnergy = 0;
	paths->dipole[0] = 0;
	paths->dipole[1] = 0;
	paths->dipole[2] = 0;

		
	paths->klass = malloc(sizeof(PathsClass));
	paths->klass->getGradVPlusHarmonic = getGradVPlusHarmonic;
	paths->klass->calculateGradient = calculateGradient;
	paths->klass->initThermostat = initThermostat;
	paths->klass->initTransformation = initTransformation;
	paths->klass->applyThermostat = applyNoThermostat;
	paths->klass->computeEnergies = computeEnergies;
	paths->klass->removeTranslation = removeTranslation;
	paths->klass->removeRotation = removeRotation;
	paths->klass->removeTranslationAndRotation = removeTranslationAndRotation;

	paths->klass->print = printPaths;
	paths->index = 0;
	paths->F = newCubeDouble(3,paths->nBeads,paths->molecules[0]->nAtoms);
	paths->klass->initTransformation(paths);// transformation before thermostat
	paths->klass->initThermostat(paths);
	paths->klass->calculateForces(paths);

	return opaths;
}
/*********************************************************************************************************************/
void freePaths(Paths* paths)
{
	int i;
	int nAtoms = paths->molecules[0]->nAtoms;
	if(paths->klass) free(paths->klass);
	if(paths->molecules) free(paths->molecules);
	if(paths->forceFields) 
	{
		if(paths->forceFields)
		for(i = 0;i<paths->nBeads; i++) 
			freeForceField(&paths->forceFields[i]);
		if(paths->forceFields) free(paths->forceFields);
		
	}
	if(paths->qmModels)
	{
		if(paths->qmModels)
		for(i = 0;i<paths->nBeads; i++) paths->qmModels[i].klass->free(&paths->qmModels[i]);
		if(paths->qmModels) free(paths->qmModels);
	}
	if(paths->M) freeMatrixDouble(&paths->M,paths->nBeads);
	if(paths->F) freeCubeDouble(&paths->F,3,paths->nBeads);
	if(paths->iNext) freeVectorInt(&paths->iNext);
	if(paths->iPrev) freeVectorInt(&paths->iPrev);
	if(paths->Mprim) freeMatrixDouble(&paths->Mprim,paths->nBeads);
	if(paths->O) freeMatrixDouble(&paths->O,paths->nBeads);
	if(paths->oneOverM) freeMatrixDouble(&paths->oneOverM,paths->nBeads);
	if(paths->LTFMult) freeMatrixDouble(&paths->LTFMult,paths->nBeads);
	if(paths->U) freeCubeDouble(&paths->U,3,paths->nBeads);
	if(paths->P) freeCubeDouble(&paths->P,3,paths->nBeads);
	if(paths->NHP) freeCubeDouble(&paths->NHP,paths->nBeads, 3*nAtoms);
	if(paths->NHF) freeCubeDouble(&paths->NHF,paths->nBeads, 3*nAtoms);
	if(paths->LTPMult) freeVectorDouble(&paths->LTPMult);
	if(paths->oneOverNHM) freeVectorDouble(&paths->oneOverNHM);
	if(paths->NHd) freeVectorDouble(&paths->NHd);
	if(paths->moleculeCentroid) paths->moleculeCentroid->klass->free(paths->moleculeCentroid);
}
/*********************************************************************************************************************/
void initPathsTemperature(Paths* paths, double newTemperature)
{
	int nAtoms = paths->molecules[0]->nAtoms;
	paths->kT = Kb*newTemperature;
	paths->temperature = newTemperature;
	//paths->wp = sqrt(1.0*paths->nBeads)*paths->kT/hbar/M_PI;
	//paths->wp = sqrt(1.0*paths->nBeads)*paths->kT/hbar/M_PI/2;
	paths->wp = sqrt(1.0*paths->nBeads)*paths->kT/hbar;
	paths->wp2 = paths->wp* paths->wp;

	//paths->numberOfPaticulesDimensionOver2Beta = paths->kT*1.5*nAtoms;
	paths->numberOfPaticulesDimensionOver2Beta = paths->kT*0.5*paths->molecules[0]->nFree;
	paths->numberOfPaticulesBeadsDimensionOver2Beta = 1.5*nAtoms*paths->nBeads*paths->kT;
}
/*********************************************************************************************************************/
void resetPathsTemperature(Paths* paths, double newTemperature)
{
	initPathsTemperature(paths, newTemperature);
	paths->klass->initThermostat(paths);
}
/*********************************************************************************************************************/
static void initIndex(Paths* paths)
{
	int i;
	int nBeads = paths->nBeads;

	paths->iNext = newVectorInt(nBeads);
	for(i = 0;i<nBeads; i++) paths->iNext[i]  = (i+1)%nBeads;
	paths->iPrev = newVectorInt(nBeads);
	for(i = 0;i<nBeads; i++) paths->iPrev[i]  = (i-1+nBeads)%nBeads;
}
/*********************************************************************************************************************/
static void calculateGradient(Paths* paths)
{
	int i;
	if(paths->forceFields) 
	for(i = 0;i<paths->nBeads; i++) paths->forceFields[i].klass->calculateGradient(&paths->forceFields[i]);
 
	if(paths->qmModels)
	for(i = 0;i<paths->nBeads; i++) paths->qmModels[i].klass->calculateGradient(&paths->qmModels[i]);
}
/*********************************************************************************************************************/
static 	void getGradVPlusHarmonic(Paths* paths, int iBead, int iAtom, double G[])
{
	int k;
	for(k=0;k<3;k++)
	G[k] = paths->molecules[iBead]->atoms[iAtom].gradient[k] 
	+ paths->M[iBead][iAtom]*
	 paths->wp2*
	paths->molecules[iBead]->atoms[iAtom].coordinates[k];
}
/*********************************************************************************************************************/
/*
static void printAllMolecule(Paths* paths)
{
	Molecule** mols = paths->molecules;
	int i;

	for(i = 0;i<paths->nBeads; i++) 
	{
		printf("Bead #%d\n",i);
		mols[i]->klass->print(mols[i],stdout);
	}
}
*/
/*********************************************************************************************************************/
static void initTransformation(Paths* paths)
{
	if(paths->transformation==PIMDTRANSFORMATIONSTAGING) initStagingTrans(paths);
	else if(paths->transformation==PIMDTRANSFORMATIONNORMALMODE) initNormalModeTrans(paths);
	else initWithoutTrans(paths);
}
/*********************************************************************************************************************/
static void initThermostat(Paths* paths)
{
	if(paths->thermostat==PIMDTHERMOSTATNONE) return;
	if(paths->thermostat==PIMDTHERMOSTATLANGEVIN) initLangevinThermostat(paths);
	if(paths->thermostat==PIMDTHERMOSTATNOSEHOOVER) initNoseHooverThermostat(paths);
	if(paths->thermostat==PIMDTHERMOSTATQTB) initQTBThermostat(paths);
}
/********************************************************************************/
/* Kinetic Virial Energy Estimator*/
static double getEKinVirial(Paths* paths)
{
	double ekin = 0.0;
	int iAtom,i,k;
	Molecule** mols = paths->molecules;
	for(i = 0;i<paths->nBeads; i++) 
	for(iAtom = 0;iAtom<mols[i]->nAtoms; iAtom++) 
	{
		for(k = 0;k<3; k++) ekin += mols[i]->atoms[iAtom].coordinates[k]*mols[i]->atoms[iAtom].gradient[k];
	}
  	ekin *= 0.5*paths->oneOvernBeads;
/*
	double ekin = 0.0;
	int iAtom,i,k;
	double dR[3];
	double RC[3];
	Molecule** mols = paths->molecules;
	int nAtoms = mols[0]->nAtoms;
	for(iAtom = 0;iAtom<nAtoms; iAtom++) 
	{
		for(k = 0;k<3; k++) RC[k] = 0;
		for(i = 0;i<paths->nBeads; i++) 
			for(k = 0;k<3; k++) RC[k] += mols[i]->atoms[iAtom].coordinates[k];

		for(k = 0;k<3; k++) RC[k] *= paths->oneOvernBeads;

		for(i = 0;i<paths->nBeads; i++) 
		{
			for(k = 0;k<3; k++) dR[k] = mols[i]->atoms[iAtom].coordinates[k]-RC[k];
			for(k = 0;k<3; k++) ekin +=dR[k]*mols[i]->atoms[iAtom].gradient[k];

		}
	}
  	ekin *= 0.5*paths->oneOvernBeads;
  	ekin += paths->numberOfPaticulesDimensionOver2Beta;

*/
  	return ekin;
}
/********************************************************************************/
/* Kinetic Primitive Energy Estimator*/
/*
static double getEKinPrimitive(Paths* paths)
{
	int iAtom,i,k;
	double dR[3];
	double ekin = 0.0;
	Molecule** mols = paths->molecules;

	for(i = 0;i<paths->nBeads; i++) 
	for(iAtom = 0;iAtom<paths->molecules[i]->nAtoms; iAtom++) 
	{
		double dot = 0;
		int iNext = paths->iNext[i];
		for(k = 0;k<3; k++) dR[k] = mols[i]->atoms[iAtom].coordinates[k]-mols[iNext]->atoms[iAtom].coordinates[k];
		for(k = 0;k<3; k++) dot += dR[k]*dR[k];
      		ekin += 0.5 * mols[i]->atoms[iAtom].mass*paths->wp2 * dot;
	}
  	ekin = paths->numberOfPaticulesBeadsDimensionOver2Beta - ekin;
	return ekin;
}
*/
/********************************************************************************/
static double getEKin(Paths* paths)
{
	//return paths->klass->getEKinVelocities(paths);
	if(paths->thermostat==PIMDTHERMOSTATQTB) return paths->klass->getEKinVelocities(paths);
	return getEKinVirial(paths);
}
/********************************************************************************/
static double getKelvin(Paths* paths)
{
	int nFree = paths->moleculeCentroid->nFree;
	if(nFree<1) return 0;
	return 2*getEKin(paths) / ( nFree * Kb);
}
/*********************************************************************************************************************/
static double getPotentialEnergy(Paths* paths)
{
	int i;
	double e = 0;
	Molecule** mols = paths->molecules;
	for(i=0;i<paths->nBeads;i++) e += mols[i]->potentialEnergy;
      	e *= paths->oneOvernBeads;
	return e;
}
/*********************************************************************************/
static void computeEnergies(Paths* paths)
{
	paths->kineticEnergy = getEKin(paths);
	paths->potentialEnergy = getPotentialEnergy(paths);
	paths->totalEnergy = paths->kineticEnergy + paths->potentialEnergy;
	paths->kelvin = getKelvin(paths);
}
/*********************************************************************************/
static void rescaleMoments(Paths* paths)
{
	double kelvin = getKelvin(paths);
	double scale = 1.0;
	Molecule** mols = paths->molecules;
	int i,j,iAtom;
	if(paths->temperature<=0) return;

	if(kelvin>0) scale = sqrt(paths->temperature/kelvin);
	if(kelvin<0) scale = -sqrt(-paths->temperature/kelvin);
#ifdef DEBUG
	printf("temp = %f kelvin = %f scale = %f\n",paths->temperature, kelvin, scale);
#endif
	paths->klass->momentToCartezian(paths);
	for(i = 0;i<paths->nBeads; i++) 
	for(iAtom = 0;iAtom<mols[i]->nAtoms; iAtom++) 
		for ( j = 0; j < 3; j++)
			paths->P[j][i][iAtom] *= scale;
}
/*********************************************************************************************************************/
static void newProperties(Paths* paths, char* comments)
{
	if( paths->fileProp == NULL) return;
	fprintf(paths->fileProp,"Time0(fs)\tTime(fs)\tTotal Energy(Kcal/mol)\tPotential Energy(kcal/mol)\tKinetic Energy(Kcal/mol)\tT(t) (K)\tTaver(K)\tsigma(T)(K)\tIndex");
	if(paths->thermostat==PIMDTHERMOSTATNOSEHOOVER) fprintf(paths->fileProp,"\tEtot+Etherm");
	if(comments) fprintf(paths->fileProp,"%s\n", comments);
	else fprintf(paths->fileProp,"\n");
}
/*********************************************************************************/
static void saveProperties(Paths* paths, int iStep0, int iStep, char* comments)
{
	double dt = paths->dt/(fsInAKMA);
	static double Ttot = 0;
	static double T2tot = 0;
	double Taver = 0;
	double T2aver = 0;

	double totalEnergy =  paths->totalEnergy;

	if( paths->fileProp == NULL) return;
	if(paths->thermostat==PIMDTHERMOSTATNOSEHOOVER)
	{
		double e = 0;
		// TO DO
		/*
		int i;
		double kT = paths->kT;
		e += paths->vNH[0]*paths->vNH[0]* paths->qNH[0]/2 + (paths->nFree)*kT* paths->xNH[0];
		for(i=1;i<MAXNH;i++) e += paths->vNH[i]*paths->vNH[i]* paths->qNH[i]/2 + kT* paths->xNH[i];
		*/
		e = getNHEnergy(paths);
		totalEnergy += e;
	}
	

	if(iStep==1)
	{
			Ttot = 0;
			T2tot = 0;
	}
	Ttot += paths->kelvin;
	T2tot += paths->kelvin*paths->kelvin;
	Taver = Ttot/iStep;
	T2aver = T2tot/iStep;


	fprintf(paths->fileProp,"%f\t%f\t%f\t\t%f\t\t\t%f\t\t\t%f\t%f\t%f\t%d\t", 
			(iStep0)*dt, 
			(iStep)*dt, 
			paths->totalEnergy,
			paths->potentialEnergy,
			paths->kineticEnergy,
			paths->kelvin,
			Taver,
			sqrt(fabs(T2aver-Taver*Taver)),
			paths->index
			 );
	if( paths->thermostat==PIMDTHERMOSTATNOSEHOOVER) fprintf(paths->fileProp,"%f\t",totalEnergy);
	if(comments) fprintf(paths->fileProp,"%s\n", comments);
	else fprintf(paths->fileProp,"\n");
}
/*********************************************************************************/
/*
static void resetGeometryCentroid(Paths* paths)
{
	int i,k,iAtom;
	Atom* atoms = paths->moleculeCentroid->atoms;

	for(iAtom = 0;iAtom<paths->moleculeCentroid->nAtoms; iAtom++) 
	for(k = 0;k<3; k++) 
		atoms[iAtom].coordinates[k] = 0.0;

	for(i = 0;i<paths->nBeads; i++) 
	for(iAtom = 0;iAtom<paths->molecules[i]->nAtoms; iAtom++) 
	for(k = 0;k<3; k++) 
		atoms[iAtom].coordinates[k] += paths->molecules[i]->atoms[iAtom].coordinates[k];

	for(iAtom = 0;iAtom<paths->molecules[i]->nAtoms; iAtom++) 
	for(k = 0;k<3; k++) 
		atoms[iAtom].coordinates[k] *= paths->oneOvernBeads;
}
*/
/*********************************************************************************/
/*
static void saveTrajectory(Paths* paths, int iStep)
{
	Molecule* mol = paths->moleculeCentroid;
	double dt = paths->dt/(fsInAKMA);
	int i;
	if( paths->fileTraj == NULL) return;
	// Get geometry from Devices

	resetGeometryCentroid(paths);
	paths->klass->computeEnergies(paths);
	fprintf(paths->fileTraj," %d %f %f %f %f nAtoms time(fs) TotalEnery(Kcal/mol) Kinetic Potential\n", 
			mol->nAtoms,
			 (iStep)*dt, 
			paths->totalEnergy,
			paths->kineticEnergy,
			paths->potentialEnergy
			 );
	fprintf(paths->fileTraj," %s\n", "Coord in Ang, Velocity in AKMA, time in fs");

	paths->klass->momentToCartezian(paths);
	for (i = 0; i < mol->nAtoms; i++)
	{
		fprintf(paths->fileTraj," %s %f %f %f %f %f %f %f %s %s %s %d %d\n", 
				mol->atoms[i].prop.symbol,
				mol->atoms[i].coordinates[0],
				mol->atoms[i].coordinates[1],
				mol->atoms[i].coordinates[2],
				mol->atoms[i].velocity[0],
				mol->atoms[i].velocity[1],
				mol->atoms[i].velocity[2],
				mol->atoms[i].charge,
				mol->atoms[i].mmType,
				mol->atoms[i].pdbType,
				mol->atoms[i].residueName,
				mol->atoms[i].residueNumber,
				mol->atoms[i].variable
				);
	}
}
*/
/*********************************************************************************************************************/
void	runPIMD(
		ForceField* forceField,
		QuantumMechanicsModel* qmModel, 
		int updateFrequency, 
		double heatTime, double equiTime, double runTime, double coolTime, 
		double heatTemperature, double equiTemperature, double runTemperature, double coolTemperature, 
		double stepSize,
		PIMDThermostatType thermostat,
		PIMDTransformationType transformation,
		double friction,
		double omegaMax,
		int Nf,
		int nBeads,
		int nNH,
		int nNHSteps,
		int nSY,
		char* fileNameTraj,
		char* fileNameProp
		)
{
	int i;
	char* str = NULL;
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
	Paths opaths;
	Paths* paths = &opaths;

	if(qmModel && qmModel->molecule.nAtoms<1) return;
	if(forceField && forceField->molecule.nAtoms<1) return;


	currentTemp = heatTemperature/2;
	
	numberOfHeatSteps = heatTime/stepSize*1000;
	numberOfEquiSteps = equiTime/stepSize*1000;; 
	numberOfRunSteps = runTime/stepSize*1000;; 
	numberOfCoolSteps = coolTime/stepSize*1000;;


	currentTemp = heatTemperature;
	if(numberOfHeatSteps==0) currentTemp = equiTemperature; 
	if(numberOfHeatSteps==0 && numberOfEquiSteps==0 ) currentTemp = runTemperature; 
	if(numberOfHeatSteps==0 && numberOfEquiSteps==0 && numberOfRunSteps==0 ) currentTemp = coolTemperature; 

		fprintf(stdout,"gamma(ps^-1)\t\t\t= %f\n",friction*fsInAKMA*1000);
	opaths = newPaths( forceField, qmModel, updateFrequency, currentTemp, stepSize, thermostat, transformation, friction, omegaMax, Nf, nBeads, nNH, nNHSteps, nSY, fileNameTraj, fileNameProp);
		fprintf(stdout,"gamma(ps^-1)\t\t\t= %f\n",paths->friction*fsInAKMA*1000);

	paths->klass->calculateGradient(paths);
	paths->klass->print(paths,stdout);

	resetPathsTemperature(paths, currentTemp);
	rescaleMoments(paths);

	paths->klass->computeEnergies(paths);
	e0 = paths->potentialEnergy;
	printf("E0 = The first potential energy in kcal = %f\n",e0); 

	n0 = 0;
	newProperties(paths," ");
	updateNumber = paths->updateFrequency;
	for (i = 0; i < numberOfHeatSteps; i++ )
	{
		paths->klass->applyOneStep(paths);
		currentTemp = heatTemperature + ( runTemperature - heatTemperature ) * ( ( double )( i + 1 )/ numberOfHeatSteps );
		resetPathsTemperature(paths, currentTemp);
		rescaleMoments(paths);
		if (++updateNumber >= paths->updateFrequency )
		{
			if(str) free(str);
			str = strdup_printf(("MD Heating: %0.2f fs, T = %0.2f K T(t) = %0.2f Kin = %0.4f Pot =  %0.4f Tot =  %0.4f"), 
					i*stepSize, 
					paths->temperature, 
					paths->kelvin, 
					paths->kineticEnergy,
					paths->potentialEnergy,
					paths->totalEnergy
					);
			printf("%s\n",str);
			updateNumber = 0;
		}
		saveProperties(paths, n0+i+1, i+1," Heating");
	}

	updateNumber = paths->updateFrequency;
	n0 += numberOfHeatSteps;
	for (i = 0; i < numberOfEquiSteps; i++ )
	{
		paths->klass->applyOneStep(paths);
		//paths->klass->applyThermostat(paths);
		if (++updateNumber >= paths->updateFrequency )
		{
			if(str) free(str);
			str = strdup_printf(("MD Equilibrium: %0.2f fs, T = %0.2f K  T(t) = %0.2f K Kin = %0.4f Pot =  %0.4f Tot =  %0.4f"), 
					i*stepSize, 
					paths->temperature, 
					paths->kelvin, 
					paths->kineticEnergy,
					paths->potentialEnergy,
					paths->totalEnergy
					);
			printf("%s\n",str);
			updateNumber = 0;
		}
		saveProperties(paths, n0+i+1, i+1, " Equilibrium");
	}
	updateNumber = paths->updateFrequency;
	n0 += numberOfEquiSteps;
	esum  = 0;
	e2sum = 0;
	for (i = 0; i < numberOfRunSteps; i++ )
	{
		paths->klass->applyOneStep(paths);
		//paths->klass->applyThermostat(paths);
		esum  += paths->totalEnergy;
		e2sum += paths->totalEnergy*paths->totalEnergy;
		if (++updateNumber >= paths->updateFrequency )
		{
			if(str) free(str);
			str = strdup_printf(("MD Running: %0.2f fs, T = %0.2f K  T(t) = %8.2f K Kin = %0.4f Pot =  %0.4f Tot =  %0.4f Eav = %0.4f sigE = %0.4f Eav-E0(cm^-1) = %0.2f"), 
					i*stepSize, 
					paths->temperature, 
					paths->kelvin, 
					paths->kineticEnergy,
					paths->potentialEnergy,
					paths->totalEnergy,
					esum/(i+1),
					sqrt(fabs(e2sum/(i+1)-esum/(i+1)*esum/(i+1))),
					(esum/(i+1)-e0)*349.75511054
					);
			printf("%s\n",str);
			updateNumber = 0;
			//saveTrajectory(paths, i+1);
		}
		saveProperties(paths, n0+i+1, i+1," Running");
	}

	updateNumber = paths->updateFrequency;
	n0 += numberOfRunSteps;
	for (i = 0; i < numberOfCoolSteps; i++ )
	{
		currentTemp = runTemperature - ( runTemperature - coolTemperature ) * ( ( double )( i + 1 )/ numberOfCoolSteps );
		resetPathsTemperature(paths, currentTemp);
		rescaleMoments(paths);
		paths->klass->applyOneStep(paths);
		if (++updateNumber >= paths->updateFrequency )
		{
			if(str) free(str);
			str = strdup_printf(("MD Cooling: %0.2f fs, T = %0.2f K T(t) = %0.2f K Kin = %0.4f Pot =  %0.4f Tot =  %0.4f"), 
					i*stepSize, 
					paths->temperature, 
					paths->kelvin, 
					paths->kineticEnergy,
					paths->potentialEnergy,
					paths->totalEnergy
					);
			printf("%s\n",str);
			updateNumber = 0;
		}
		saveProperties(paths, n0+i+1, i+1," Cooling");
	}
	paths->klass->calculateGradient(paths);
	if(str) free(str);
	str = strdup_printf(("End of MD Simulation. Ekin = %f (Kcal/mol) EPot =  %0.4f ETot =  %0.4f T(t) = %0.2f"),
			paths->kineticEnergy,
			paths->potentialEnergy,
			paths->totalEnergy,
			paths->kelvin 
			); 
	printf("%s\n",str);
	free(str);
	if(paths->fileTraj)fclose(paths->fileTraj);
	if(paths->fileProp)fclose(paths->fileProp);
	freePaths(paths);
}
/**************************************************************************************************************************************/
static void printPaths(Paths* paths,FILE* file)
{
	int i;
	int iAtom = 0;
	Molecule* molecule = paths->moleculeCentroid;
	fprintf(file,"\n");
	fprintf(file,"*************** PIMD Parameters ************************************************************************************************\n");
	fprintf(file,"nBeads\t\t\t\t\t= %d\n",paths->nBeads);
	fprintf(file,"dt(fs)\t\t\t\t\t= %f\n",paths->dt/fsInAKMA);
	if(paths->thermostat==PIMDTHERMOSTATLANGEVIN)
	{
		fprintf(file,"Thermostat\t\t\t= %s\n","Langevin thermostat");
		fprintf(file,"Temperature\t\t\t= %f\n",paths->temperature);
		fprintf(file,"gamma(ps^-1)\t\t\t= %f\n",paths->friction*fsInAKMA*1000);
	}
	if(paths->thermostat==PIMDTHERMOSTATNOSEHOOVER)
	{
		fprintf(file,"Thermostat\t\t\t\t= %s\n","Nose-Hoover thermostat");
		fprintf(file,"Temperature\t\t\t\t= %f\n",paths->temperature);
		fprintf(file,"Length of the Nose-Hoover chain\t\t= %d\n",paths->nNH);
		fprintf(file,"Number of Nose-Hoover Steps\t\t= %d\n",paths->nNHSteps);
		fprintf(file,"Number of Suzuki-Yoshida weights\t= %d\n",paths->nSY);
	}
	if(paths->thermostat==PIMDTHERMOSTATQTB)
	{
		fprintf(file,"Thermostat\t\t\t\t= %s\n","Qauntum thermal Bath");
		fprintf(file,"Temperature\t\t\t\t= %f\n",paths->temperature);
		printf("Nf\t\t= %d\n",paths->Nf);
		printf("Number of QTB steps(M)\t\t= %d\n",paths->nQTBSteps);
		printf("dt(fs)\t\t= %f\n",paths->dt/fsInAKMA);
		printf("h(fs)\t\t= %f\n",paths->h/fsInAKMA);
		printf("gamma(ps^-1)\t= %f\n",paths->friction*fsInAKMA*1000);
		printf("omegaMax(cm^-1)\t= %f\n",paths->omegaMax);
	}
	if(paths->thermostat==PIMDTHERMOSTATNONE)
	{
		fprintf(file,"Thermostat\t\t\t\t= %s\n","None");
		fprintf(file,"Temperature\t\t\t\t= %f\n",paths->temperature);
	}
	if(paths->transformation==PIMDTRANSFORMATIONSTAGING)
	{
		fprintf(file,"Transformation\t\t\t\t= %s\n","Staging");
	}
	if(paths->transformation==PIMDTRANSFORMATIONNORMALMODE)
	{
		fprintf(file,"Transformation\t\t\t\t= %s\n","Normal Mode");
	}
	if(paths->transformation==PIMDTRANSFORMATIONNONE)
	{
		fprintf(file,"Transformation\t\t\t\t= %s\n","None");
	}
	fprintf(file,"wp(rad/ps)\t\t\t\t= %f\n",paths->wp*fsInAKMA*1000);
	fprintf(file,"nFree\t\t\t\t\t= %d\n", paths->moleculeCentroid->nFree);
	fprintf(file,"hbar (AKMA)\t\t\t\t= %f\n", hbar);
	for(i = 0;i<paths->nBeads; i++) 
	{
		fprintf(file,"Bead #\t\t= %d\n",i);
		fprintf(file,"%-6s %10s %10s %10s %10s %10s %10s %10s %10s %10s %10s %10s\n","Symbol","M","M'","X","Y","Z","Vx","Vy","Vz","Fx","Fy","Fz");
	for(iAtom = 0;iAtom<paths->molecules[i]->nAtoms; iAtom++) 
		fprintf(file,"%-6s %10.6f %10.6f %10.6f %10.6f %10.6f %10.6f %10.6f %10.6f %10.6f %10.6f %10.6f\n",
			molecule->atoms[iAtom].prop.symbol,
			paths->M[i][iAtom],
			paths->Mprim[i][iAtom],
			molecule->atoms[iAtom].coordinates[0],
			molecule->atoms[iAtom].coordinates[1],
			molecule->atoms[iAtom].coordinates[2],
			molecule->atoms[iAtom].velocity[0],
			molecule->atoms[iAtom].velocity[1],
			molecule->atoms[iAtom].velocity[2],
			paths->F[0][i][iAtom],
			paths->F[1][i][iAtom],
			paths->F[2][i][iAtom]
			);
	}
	fprintf(file,"********************************************************************************************************************************\n");
	fprintf(file,"\n");
}
/*****************************************************************************************************************************/
static void applyNoThermostat(Paths* paths)
{
}
/*********************************************************************************************************************/
static void initVelocities(Paths* paths)
{
	int i;
	//int k,iAtom;
	Molecule* molC = paths->moleculeCentroid;
	Molecule** mols = paths->molecules;
	molC->klass->setMaxwellVelocitiesIfNull(molC, paths->temperature);
	/*
        for ( i = 0; i < paths->nBeads; i++)
        for ( iAtom = 0; iAtom < mols[i]->nAtoms; iAtom++)
	for(k = 0;k<3; k++) 
                	mols[i]->atoms[iAtom].velocity[k] = molC->atoms[iAtom].velocity[k];
	*/
        for ( i = 0; i < paths->nBeads; i++) mols[i]->klass->setMaxwellVelocitiesIfNull(mols[i], paths->temperature);
        for ( i = 0; i < paths->nBeads; i++) mols[i]->klass->removeTranslationAndRotation(mols[i]);
}
/*****************************************************************************************************************************/
static void initAllPointers(Paths* paths)
{
	paths->fileTraj = NULL;
	paths->fileProp = NULL;
	paths->forceFields = NULL;
	paths->qmModels = NULL;
	paths->molecules = NULL;
	paths->moleculeCentroid = NULL;
	paths->M = NULL;
	paths->Mprim = NULL;
	paths->O = NULL;
	paths->F = NULL;
	paths->U = NULL;
	paths->P = NULL;

	paths->oneOverM = NULL;
	paths->iNext = NULL;
	paths->iPrev = NULL;

	paths->LTPMult = NULL;
	paths->LTFMult = NULL;

	paths->NHP = NULL;
	paths->NHF = NULL;
	paths->oneOverNHM = NULL;
	paths->NHd = NULL;
}
/*********************************************************************************/
static void removeTranslation(Paths* paths)
{
	double vtot[3] = {0,0,0};
	int i;
	int iAtom;
	int j;
	double tot = 0.0;
	double*** P = paths->P;

	if(!P) return;
	for ( i = 0; i < paths->nBeads; i++)
	{
		for ( iAtom = 0; iAtom < paths->molecules[i]->nAtoms; iAtom++)
		{
			tot += paths->Mprim[i][iAtom];
			for ( j = 0; j < 3; j++)
			{
				vtot[j] += P[j][i][iAtom];
			}
		}
	}

	for ( j = 0; j < 3; j++) vtot[j] /= tot;


	for ( i = 0; i < paths->nBeads; i++)
		for ( iAtom = 0; iAtom < paths->molecules[i]->nAtoms; iAtom++)
			for ( j = 0; j < 3; j++)
				P[j][i][iAtom] -= paths->Mprim[i][iAtom]*vtot[j];
}
/*********************************************************************************/
static void removeRotation(Paths* paths)
{
	double vtot[3] = {0,0,0};
	double cm[3] = {0,0,0};
	double L[3] = {0,0,0};
	int i;
	int iAtom;
	int j;
	int k;
	double mass = 1.0;
	double totMass = 0.0;
	double cdel[3];
	double vAng[3]={0,0,0};
	double tensor[3][3];
	double invTensor[3][3];
        double xx, xy,xz,yy,yz,zz;
	double*** P = paths->P;
	double*** U = paths->U;
	double n = 0;
	/* find the center of mass coordinates  and total velocity*/
	if(!U)
	{
		paths->molecules[0]->klass->removeRotationMoments(paths->molecules, paths->nBeads, paths->P);
		return;
	}

	for ( i = 0; i < paths->nBeads; i++)
	{
		for ( iAtom = 0; iAtom <  paths->molecules[i]->nAtoms; iAtom++)
		{
			mass = paths->Mprim[i][iAtom];
			totMass += mass;
			n++;
			for ( j = 0; j < 3; j++)
				cm[j] += mass*U[j][i][iAtom];
			for ( j = 0; j < 3; j++)
				vtot[j] += P[j][i][iAtom];
		}
	}


	for ( j = 0; j < 3; j++) cm[j] /= totMass;

	/*   compute the angular momentum  */
	for ( i = 0; i < paths->nBeads; i++)
	{
		for ( iAtom = 0; iAtom <  paths->molecules[i]->nAtoms; iAtom++)
		{
			mass = paths->Mprim[i][iAtom];
			for ( j = 0; j < 3; j++)
			L[j] += (
				U[(j+1)%3][i][iAtom]*P[(j+2)%3][i][iAtom]
			      - U[(j+2)%3][i][iAtom]*P[(j+1)%3][i][iAtom]
			      );
		}
	}
	for ( j = 0; j < 3; j++)
		L[j] -= (
			cm[(j+1)%3]*vtot[(j+2)%3]
		      - cm[(j+2)%3]*vtot[(j+1)%3]
			      );

	/* calculate and invert the inertia tensor */
	for ( k = 0; k < 3; k++)
	for ( j = 0; j < 3; j++)
		tensor[k][j] = 0;
	xx = 0;
	yy = 0;
	zz = 0;
	xy = 0;
	xz = 0;
	yz = 0;
	for ( i = 0; i < paths->nBeads; i++)
	{
		for ( iAtom = 0; iAtom <  paths->molecules[i]->nAtoms; iAtom++)
		{
			mass = paths->Mprim[i][iAtom];
			for ( j = 0; j < 3; j++)
				cdel[j] = U[j][i][iAtom]-cm[j];
			xx +=  cdel[0]*cdel[0]*mass;
			xy +=  cdel[0]*cdel[1]*mass;
			xz +=  cdel[0]*cdel[2]*mass;
			yy +=  cdel[1]*cdel[1]*mass;
			yz +=  cdel[1]*cdel[2]*mass;
			zz +=  cdel[2]*cdel[2]*mass;
		}
	}
	tensor[0][0] = yy+zz;
	tensor[1][0] = -xy;
	tensor[2][0] = -xz;
	tensor[0][1] = -xy;
	tensor[1][1] = xx+zz;
	tensor[2][1] = -yz;
	tensor[0][2] = -xz;
	tensor[1][2] = -yz;
	tensor[2][2] = xx+yy;
	if(InverseTensor(tensor,invTensor))
	{
		for ( j = 0; j < 3; j++)
		{
			vAng[j] = 0;
			for ( k = 0; k < 3; k++)
				vAng[j] += invTensor[j][k]*L[k];
		}
	}
	else
	if( paths->molecules[0]->nAtoms>1)
	{
		double U0[3];
		double U1[3];
		for ( j = 0; j < 3; j++)U0[j] = U[j][0][0];
		for ( j = 0; j < 3; j++)U1[j] = U[j][0][1];
		//printf("!!!!!!!!!!!I cannot invert the rotational Tensor : linear molecule!\n");
		computeAngularVelocitiesForALinearMolecule(U0, U1, tensor, L, vAng);
		//printf("Angular velocityibefore rotation = %f %f %f\n",vAng[0], vAng[1], vAng[2]);
	}
	/*  eliminate any rotation about the system center of mass */
	for ( i = 0; i < paths->nBeads; i++)
	{
		for ( iAtom = 0; iAtom <  paths->molecules[i]->nAtoms; iAtom++)
		{
			for ( j = 0; j < 3; j++)
				cdel[j] = U[j][i][iAtom]-cm[j];
			for ( j = 0; j < 3; j++)
				P[j][i][iAtom] += 
				(cdel[(j+1)%3]*vAng[(j+2)%3]-
				cdel[(j+2)%3]*vAng[(j+1)%3])
				*paths->Mprim[i][iAtom]
				;
		}
	}
}
/*********************************************************************************/
static void removeTranslationAndRotation(Paths* paths)
{
	if(paths->nBeads<1) return;
	paths->klass->removeTranslation(paths);
	paths->klass->removeRotation(paths);
}
