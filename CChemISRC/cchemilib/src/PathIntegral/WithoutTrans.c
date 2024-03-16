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
#include "../Utils/Constants.h"
#include "../Utils/Utils.h"

static void calculateForcesWithoutTrans(Paths* paths);
static void applyOneStepWithoutTrans(Paths* paths);
static void momentToCartezian(Paths* paths);
static void initMasses(Paths* paths);
static void initMoments(Paths* paths);
static double getEKinVelocities(Paths* paths);
/*********************************************************************************************************************/
void initWithoutTrans(Paths* paths)
{
	initMasses(paths);
	initMoments(paths);

	paths->klass->calculateForces = calculateForcesWithoutTrans;
	paths->klass->applyOneStep = applyOneStepWithoutTrans;
	paths->klass->momentToCartezian = momentToCartezian ;
	paths->klass->getEKinVelocities = getEKinVelocities;
	paths->klass->calculateForces(paths);
	printf("End initWithoutTrans\n");
}
/*********************************************************************************************************************/
static void initMasses(Paths* paths)
{
	int i,iAtom;
	int nAtoms = paths->molecules[0]->nAtoms;
	int nBeads = paths->nBeads;

	paths->M = newMatrixDouble(nBeads,nAtoms);
	for(i = 0;i<paths->nBeads; i++) 
	for(iAtom = 0;iAtom<paths->molecules[i]->nAtoms; iAtom++) 
		paths->M[i][iAtom]  = paths->molecules[i]->atoms[iAtom].mass;

	paths->Mprim = newMatrixDouble(nBeads,nAtoms);
	for(i = 0;i<paths->nBeads; i++) 
	for(iAtom = 0;iAtom<paths->molecules[i]->nAtoms; iAtom++) 
		paths->Mprim[i][iAtom]  = paths->molecules[i]->atoms[iAtom].mass;
}
/*********************************************************************************************************************/
static void initMoments(Paths* paths)
{
	int i,k,iAtom;
	int nAtoms = paths->molecules[0]->nAtoms;
	int nBeads = paths->nBeads;

	paths->P = newCubeDouble(3,nBeads,nAtoms);
	for(k = 0;k<3; k++) 
	for(i = 0;i<paths->nBeads; i++) 
	for(iAtom = 0;iAtom<paths->molecules[i]->nAtoms; iAtom++) 
		paths->P[k][i][iAtom]  = paths->molecules[i]->atoms[iAtom].velocity[k]*paths->Mprim[i][iAtom];
}
/**************************************************************************************************************************/
static void applyOneStepWithoutTrans(Paths* paths)
{
	int iAtom,i,k;
	Molecule** mols;
	mols = paths->molecules;
	double*** P = paths->P;

	//paths->klass->applyThermostat(paths); // to be removed

	paths->klass->applyThermostat(paths); // should be apply 2 times
	//printAllMolecule(paths);
	for(i = 0;i<paths->nBeads; i++) 
	for(iAtom = 0;iAtom<paths->molecules[i]->nAtoms; iAtom++) 
	for(k = 0;k<3; k++) 
	{
		P[k][i][iAtom] += 0.5*paths->F[k][i][iAtom]*paths->dt;
		mols[i]->atoms[iAtom].coordinates[k] +=  P[k][i][iAtom]*paths->dt/paths->Mprim[i][iAtom];
	}
	paths->klass->calculateForces(paths);
	paths->klass->applyThermostat(paths);
	for(i = 0;i<paths->nBeads; i++) 
	for(iAtom = 0;iAtom<paths->molecules[i]->nAtoms; iAtom++) 
	for(k = 0;k<3; k++) 
	{
		P[k][i][iAtom] += 0.5*paths->F[k][i][iAtom]*paths->dt;
	}
	//printf("After the first step\n");
	//printAllMolecule(paths);
  
	//paths->klass->applyThermostat(paths);
	paths->klass->removeTranslationAndRotation(paths);
	paths->klass->computeEnergies(paths);
  
}
/*********************************************************************************************************************/
static void calculateForcesWithoutTrans(Paths* paths)
{
	double dR1[3];
	double dR2[3];
	int iAtom,i,iNext,iPrev;
	int k;
	int nAtoms = paths->molecules[0]->nAtoms;
	Molecule** mols;

	paths->klass->calculateGradient(paths);
	
	mols = paths->molecules;

	for (iAtom = 0; iAtom < nAtoms; iAtom++)
    	for (i = 0; i <paths->nBeads; i++)
	{
		iNext = paths->iNext[i];
		iPrev = paths->iPrev[i];
		//printf("i, IP, IM %d %d %d\n",i,iNext,iPrev);
		for(k=0;k<3;k++) dR1[k] = mols[iNext]->atoms[iAtom].coordinates[k]-mols[i]->atoms[iAtom].coordinates[k];
		for(k=0;k<3;k++) dR2[k] = mols[iPrev]->atoms[iAtom].coordinates[k]-mols[i]->atoms[iAtom].coordinates[k];
		//printf("dR1 = %f %f %f\n",dR1[0],dR1[1],dR1[2]);
		//printf("dR2 = %f %f %f\n",dR2[0],dR2[1],dR2[2]);
		for(k=0;k<3;k++) paths->F[k][i][iAtom] = paths->M[i][iAtom]*paths->wp2*(dR1[k]+dR2[k])
			-paths->oneOvernBeads*mols[i]->atoms[iAtom].gradient[k];
	}
}
/*****************************************************************************************************************************/
static void momentToCartezian(Paths* paths)
{
	int iAtom,i,k;
	Molecule** mols = paths->molecules;

    	for (i = 0; i<paths->nBeads; i++)
		for (iAtom = 0; iAtom < mols[i]->nAtoms; iAtom++)
			for(k = 0;k<3; k++) 
				mols[i]->atoms[iAtom].velocity[k] =  paths->P[k][i][iAtom]/paths->Mprim[i][iAtom];

	return;
}
/********************************************************************************/
/* Kinetic Velocities Energy Estimator*/
static double getEKinVelocities(Paths* paths)
{
	double ekin = 0.0;
	int iAtom,i,j;
	//double RC[3];
	double***P = paths->P;
	Molecule** mols = paths->molecules;
	double dR[3];
	int k;

	double invmass;
	ekin = 0;
	for(i = 0;i<paths->nBeads; i++) 
	for(iAtom = 0;iAtom< mols[i]->nAtoms; iAtom++) 
	{
		invmass =  1.0/paths->Mprim[i][iAtom];
		for ( j = 0; j < 3; j++)
			ekin += P[j][i][iAtom]*P[j][i][iAtom]*invmass;
	}

	for(i = 0;i<paths->nBeads; i++) 
	for(iAtom = 0;iAtom<paths->molecules[i]->nAtoms; iAtom++) 
	{
		double dot = 0;
		int iNext = paths->iNext[i];
		for(k = 0;k<3; k++) dR[k] = mols[i]->atoms[iAtom].coordinates[k]-mols[iNext]->atoms[iAtom].coordinates[k];
		for(k = 0;k<3; k++) dot += dR[k]*dR[k];
      		ekin -=  paths->M[i][iAtom]*paths->wp2 * dot;
	}
  	ekin *= 0.5;
  	return ekin;
}
