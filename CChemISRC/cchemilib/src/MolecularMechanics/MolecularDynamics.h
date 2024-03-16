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

#ifndef __CCHEMILIB_MOLECULARDYNAMICS_H__
#define __CCHEMILIB_MOLECULARDYNAMICS_H__

#include "../Molecule/Molecule.h"
#include "../MolecularMechanics/ForceField.h"

typedef struct _MolecularDynamics  MolecularDynamics;

struct _MolecularDynamics
{
	ForceField* forceField;
	int numberOfAtoms;
	int updateFrequency;
	double** coordinatesOld; 
	double** a; 
	double** aold; 
	boolean* moved;
	boolean* update;
	double dt;
	double dt_2;
	double dt_4;
	double dt2_2;
	double dt_8;
	double dt2_8;
	double dt2;
	double dt3_6;
	double dt4_48;
	MDIntegratorType integratorType;
	double potentialEnergy;
	double kineticEnergy;
	double totalEnergy;
	double kelvin;


	double* positionFriction;
	double* velocityFriction;
	double* accelarationFriction;
	double** positionRandom;
	double** velocityRandom;
	double* gamma;
	double friction;
	double temperature;
	double collide;
	MDThermostatType thermostat;

/* QTB && LANGEVIN */
	double* theta;
/* QTB */
	double** rnoise;
	double* Ht;
	int Nf;
	double h;
	int M;

	FILE* fileTraj;
	FILE* fileProp;
	int nvariables;
	int index;
#define MAXNH 2
	double xNH[MAXNH];
	double vNH[MAXNH];
	double qNH[MAXNH];
	double gNH[MAXNH];

/* Martyna&Tuckerman Symplectic reversible integrators: Predictor-corrector methods */
/* Ref : JCP 102, 8071 (1995) */
	double* Ftilde[3];
	double* VF[3];
	double* VJ[3];

#ifdef ENABLE_CL
	cl_program programMD;

	cl_kernel copyAccelarations;
	cl_kernel computeAccelarationsFromGradients;
	cl_kernel copyPositions;
	cl_kernel initAccelarations;
	cl_kernel applyVerlet1;
	cl_kernel applyVerlet2;
	cl_kernel applyBeeman1;
	cl_kernel applyBeeman2;
	cl_kernel applyStochastic1;
	cl_kernel applyStochastic2;
	cl_kernel setFrictionalAndRandomForce;
	cl_kernel computeKineticEnergy;
	cl_kernel scaleVelocities;
	cl_kernel generateRandomNumbers;
	cl_kernel andersen;
	cl_kernel removeTranslation;
	cl_kernel removeRotation;
	cl_kernel initRattles;
	cl_kernel updateRattle1;
	cl_kernel testDoneRattle;
	cl_kernel applyRattle1;
	cl_kernel updateRattle2;
	cl_kernel applyRattle2;
	cl_mem  aCL;
	cl_mem seedCL;
	cl_mem randomsCL;
	cl_int nRandomsCL;
	cl_mem oldPositionBufferCL;
	cl_mem frictionBufferCL;

#endif
};
void	freeMolecularDynamics(MolecularDynamics* molecularDynamics);
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
		);
ForceField**	runMolecularDynamicsConfo(
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
		);
ForceField**    runREMD(
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
		int nTemperatures,
		int numberOfExchanges,
		char* fileNameTraj,
		char* fileNameProp
		);
#endif /* __CCHEMILIB_MOLECULARDYNAMICS_H__ */

