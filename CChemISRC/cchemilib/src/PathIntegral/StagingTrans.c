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

static void calculateForcesStaging(Paths* paths);
static void RtoUStaging(Paths* paths);
static void UtoRStaging(Paths* paths);
static void applyOneStepStaging(Paths* paths);
static void momentToCartezianStaging(Paths* paths);
static void initMasses(Paths* paths);
static void initMoments(Paths* paths);
static double getEKinVelocities(Paths* paths);
/*********************************************************************************************************************/
void initStagingTrans(Paths* paths)
{
	int nAtoms = paths->molecules[0]->nAtoms;
	int nBeads = paths->nBeads;

  	/* Positions U */
	paths->U = newCubeDouble(3,nBeads,nAtoms);
	RtoUStaging(paths);

	initMasses(paths);
	initMoments(paths);

	paths->klass->calculateForces = calculateForcesStaging;
	paths->klass->applyOneStep = applyOneStepStaging;
	paths->klass->momentToCartezian = momentToCartezianStaging ;
	paths->klass->getEKinVelocities = getEKinVelocities;
	paths->klass->calculateForces(paths);
	printf("End initStagingTrans\n");
}
/*********************************************************************************************************************/
static void initMasses(Paths* paths)
{
	int i,iAtom;
	int nAtoms = paths->molecules[0]->nAtoms;
	int nBeads = paths->nBeads;

	paths->M  = newMatrixDouble(nBeads,nAtoms);
	for(iAtom = 0;iAtom<nAtoms; iAtom++) paths->M[0][iAtom]  = 0;
	for(i = 1;i<paths->nBeads; i++) 
	for(iAtom = 0;iAtom<nAtoms; iAtom++) 
		paths->M[i][iAtom]  = (i+1.0)/i*paths->molecules[i]->atoms[iAtom].mass;

	paths->Mprim  = newMatrixDouble(nBeads,nAtoms);
	for(iAtom = 0;iAtom<nAtoms; iAtom++) paths->Mprim[0][iAtom]  = paths->molecules[0]->atoms[iAtom].mass;

	for(i = 1;i<paths->nBeads; i++) 
	for(iAtom = 0;iAtom<nAtoms; iAtom++) 
		//paths->Mprim[i][iAtom]  =  paths->molecules[i]->atoms[iAtom].mass;
		paths->Mprim[i][iAtom]  = paths->M[i][iAtom];

}
/*********************************************************************************************************************/
static void initMoments(Paths* paths)
{
	int i,k,iAtom;
	int nAtoms = paths->molecules[0]->nAtoms;
	int nBeads = paths->nBeads;

	paths->P = newCubeDouble(3,nBeads,nAtoms);

	i = 0;
	for(iAtom = 0;iAtom<nAtoms; iAtom++) 
	for(k=0;k<3;k++)
		paths->P[k][i][iAtom] = paths->molecules[i]->atoms[iAtom].velocity[k]*paths->Mprim[i][iAtom];

	for(i = paths->nBeads - 1; i>0 ;  i--) 
	{
		int iN = paths->iNext[i];
		for(iAtom = 0;iAtom<nAtoms; iAtom++) 
		for(k=0;k<3;k++)
		{
			paths->P[k][i][iAtom] =
      			paths->molecules[i]->atoms[iAtom].velocity[k]*paths->Mprim[i][iAtom] 
			- i/(i + 1.0)*paths->molecules[iN]->atoms[iAtom].velocity[k]*paths->Mprim[iN][iAtom]
			- 1.0/(i + 1.0)*paths->P[k][0][iAtom];
		}
	}
}
/**************************************************************************************************************************/
static void applyOneStepStaging(Paths* paths)
{
	int iAtom,i,k;
	double*** P = paths->P;
	double*** U = paths->U;
	paths->klass->applyThermostat(paths); // 2 times

	for(i = 0;i<paths->nBeads; i++) 
	for(iAtom = 0;iAtom<paths->molecules[i]->nAtoms; iAtom++) 
	for(k = 0;k<3; k++) 
	{
		P[k][i][iAtom] += 0.5*paths->F[k][i][iAtom]*paths->dt;
		U[k][i][iAtom] +=  P[k][i][iAtom]*paths->dt/paths->Mprim[i][iAtom];
	}
	UtoRStaging(paths);
	paths->klass->calculateForces(paths);
	paths->klass->applyThermostat(paths);

	for(i = 0;i<paths->nBeads; i++) 
	for(iAtom = 0;iAtom<paths->molecules[i]->nAtoms; iAtom++) 
	for(k = 0;k<3; k++) 
	{
		P[k][i][iAtom] += 0.5*paths->F[k][i][iAtom]*paths->dt;
	}
  
	paths->klass->removeTranslationAndRotation(paths);
	paths->klass->computeEnergies(paths);
  
}

/**************************************************************************************************************************/
static void UtoRStaging(Paths* paths)
{
	int iAtom,i,k;
	Molecule** mols = paths->molecules;
	int nAtoms = paths->molecules[0]->nAtoms;
	for (iAtom = 0; iAtom < nAtoms; iAtom++)
	{
		// first bead
		i = 0;
		for(k = 0;k<3; k++) mols[i]->atoms[iAtom].coordinates[k] = paths->U[k][i][iAtom];
    		for (i = paths->nBeads-1; i > 0; i--)
		{
			int iN = paths->iNext[i];
			for(k = 0;k<3; k++) 
				mols[i]->atoms[iAtom].coordinates[k] = 
      				paths->U[k][i][iAtom] + i/(i + 1.0)* mols[iN]->atoms[iAtom].coordinates[k]
				+ 1.0/(i + 1.0)*paths->U[k][0][iAtom];
		}
	}
}
/**************************************************************************************************************************/
static void RtoUStaging(Paths* paths)
{
	int iAtom,i,k;
	Molecule** mols = paths->molecules;
	int nAtoms = paths->molecules[0]->nAtoms;
	for (iAtom = 0; iAtom < nAtoms; iAtom++)
	{
		i = 0;
		for(k = 0;k<3; k++) 
				paths->U[k][i][iAtom] = mols[i]->atoms[iAtom].coordinates[k];
    		for (i = paths->nBeads-1; i > 0; i--)
		{
			int iN = paths->iNext[i];
			for(k = 0;k<3; k++) 
      				paths->U[k][i][iAtom]  =
				mols[i]->atoms[iAtom].coordinates[k] 
				- i/(i + 1.0)* mols[iN]->atoms[iAtom].coordinates[k]
				- 1.0/(i + 1.0)*paths->U[k][0][iAtom];
		}
	}
}

/*********************************************************************************************************************/
static void calculateForcesStaging(Paths* paths)
{
  	double  gradVA[3];
  	double  gradVB[3];
	int iAtom,i;
	int k;
	int nAtoms = paths->molecules[0]->nAtoms;
	Molecule** mols = paths->molecules;

	paths->klass->calculateGradient(paths);

	for (iAtom = 0; iAtom < nAtoms; iAtom++)
	{
		for(k = 0;k<3; k++) gradVA[k] = 0;

    		for (i = 0; i <paths->nBeads; i++)
		{
		for(k = 0;k<3; k++) gradVA[k] += mols[i]->atoms[iAtom].gradient[k];
		}

		for(k=0;k<3;k++) paths->F[k][0][iAtom] = -paths->oneOvernBeads*gradVA[k];

		for(k = 0;k<3; k++) gradVB[k] = gradVA[k];

    		for (i = 1; i <paths->nBeads; i++)
		{
			for(k = 0;k<3; k++) gradVA[k] =  mols[i]->atoms[iAtom].gradient[k]+(i - 1.0)/i*gradVB[k];
			for(k=0;k<3;k++) paths->F[k][i][iAtom] = -paths->oneOvernBeads*gradVA[k];
			for(k=0;k<3;k++) paths->F[k][i][iAtom] += -paths->M[i][iAtom]*paths->wp2*paths->U[k][i][iAtom];
			for(k = 0;k<3; k++) gradVB[k] = gradVA[k];
		}
	}
}
/**************************************************************************************************************************/
static void momentToCartezianStaging(Paths* paths)
{
	int iAtom,i,k;
	Molecule** mols = paths->molecules;
	int nAtoms = paths->molecules[0]->nAtoms;
	double*** P = paths->P;

	for (iAtom = 0; iAtom < nAtoms; iAtom++)
	{
		// first bead
		i = 0;
		for(k = 0;k<3; k++) mols[i]->atoms[iAtom].velocity[k] = P[k][i][iAtom];
    		for (i = paths->nBeads-1; i > 0; i--)
		{
			int iN = paths->iNext[i];
			for(k = 0;k<3; k++) 
				mols[i]->atoms[iAtom].velocity[k] = 
      				P[k][i][iAtom] + i/(i + 1.0)* mols[iN]->atoms[iAtom].velocity[k]
				+ 1.0/(i + 1.0)*P[k][0][iAtom];
		}
	}
    	for (i = 0; i<paths->nBeads; i++)
	for (iAtom = 0; iAtom < nAtoms; iAtom++)
		for(k = 0;k<3; k++) mols[i]->atoms[iAtom].velocity[k] /= paths->Mprim[i][iAtom];
}
/********************************************************************************/
/* Kinetic Velocities Energy Estimator*/
static double getEKinVelocities(Paths* paths)
{
	double ekin = 0.0;
	int iAtom,i,j;
	//double RC[3];
	double***P = paths->P;
	double***U = paths->U;
	Molecule** mols = paths->molecules;
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
		for(k = 0;k<3; k++) dot += U[k][i][iAtom]*U[k][i][iAtom];
      		ekin -=  paths->M[i][iAtom]*paths->wp2 * dot;
	}
  	ekin *= 0.5;
  	return ekin;
}
