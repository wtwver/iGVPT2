/********************************************************************************
 cchemi is an interface to ab initio computational chemistry programs 
 designed for add them many functionalities non available in these packages.
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
********************************************************************************/

#ifndef __CCHEMILIB_PATHS_H__
#define __CCHEMILIB_PATHS_H__

#include <stdio.h>
#include "../Molecule/Molecule.h"
#include "../Utils/Types.h"
#include "../MolecularMechanics/ForceField.h"
#include "../QuantumMechanics/QuantumMechanicsModel.h"
#include "../Utils/Constants.h"
#include "../Utils/Utils.h"

typedef struct _Paths  Paths;
typedef struct _PathsClass  PathsClass;


struct _PathsClass
{
	void (*getGradVPlusHarmonic)(Paths* paths, int iBead, int iAtom, double G[]);
	void (*calculateGradient)(Paths* paths);
	void (*calculateForces)(Paths* paths);
	void (*applyOneStep)(Paths* paths);
	void (*applyThermostat)(Paths* paths);
	void (*initThermostat)(Paths* paths);
	void (*initTransformation)(Paths* paths);
	void (*momentToCartezian)(Paths* paths);
	void (*computeEnergies)(Paths* paths);
	void (*print)(Paths* paths,FILE* file);
	void (*removeTranslation)(Paths* paths);
	void (*removeRotation)(Paths* paths);
	void (*removeTranslationAndRotation)(Paths* paths);
	double (*getEKinVelocities)(Paths* paths);
};
struct _Paths
{
	PathsClass* klass;
	int nBeads;
	FILE* fileTraj;
	FILE* fileProp;
	ForceField* forceFields;
	QuantumMechanicsModel* qmModels;
	Molecule** molecules;
	Molecule* moleculeCentroid;
	PIMDThermostatType thermostat;
	PIMDTransformationType transformation;
	double** M;
	double** Mprim;
	double** O;
	double*** F;
	double*** U;
	double*** P;

	int updateFrequency;
	double temperature;
	double kT;
	double dt;
	double wp;
	double wp2;
	double oneOverNumberOfParticules;
	double oneOvernBeads;
	double oneOverNumberOfParticulesBeads;
	double numberOfPaticulesBeadsDimensionOver2Beta;
	double numberOfPaticulesDimensionOver2Beta;
	double** oneOverM;
	int* iNext;
	int* iPrev;

	double potentialEnergy;
	double kineticEnergy;
	double totalEnergy;
	double kelvin;
	double dipole[3];
/* Langevin thermostat parameters */
	double friction;// used also in QTB
	double* LTPMult; /* Langevin Momentum Multiplier*/
	double** LTFMult; /* Langevin Force Multiplier*/

/* Nose-Hoover Thermostat*/
	int nNH; /* Length of Nose-Hoover chain*/
	double*** NHP; /* Momenta (Nose-Hoover)*/
	double*** NHF; /* Forces (Nose-Hoover)*/
	double* oneOverNHM; /* Inverse Masses (Nose-Hoover)*/
	int nSY; /* Number of Suzuki-Yoshida weights*/
	int nNHSteps; /* Number of Nose-Hoover Steps*/
	double* NHd; /* Suzuki-Yoshia Weight Factors*/
	double xNHTotal; 

/* QTB */
	double omegaMax;
	double** theta;
	double** Ht;
	double*** rnoise;
	int Nf;
	double h;
	int nQTBSteps;
	double gpQTB;
	double gmQTB;
	int iQTBStep;

	int index; // to be used with RE-PIMD in the future
};
Paths newPaths(ForceField* forceField, QuantumMechanicsModel* qmModel, 
		int updateFrequency, double temperature, double stepSize,
		PIMDThermostatType thermostat, PIMDTransformationType tansformation,
		double friction, double omegaMax, int Nf, int nBeads, int nNH, int nNHSteps, int nSY,
		char* fileNameTraj, char* fileNameProp
);
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
		);
#endif /* __CCHEMILIB_PATHS_H__ */
