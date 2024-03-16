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
#include "../Utils/QL.h"
#include "../Utils/Utils.h"

/* Refs :
Perez et al : THE JOURNAL OF CHEMICAL PHYSICS 130, 184105 2009
*/

static void calculateForcesNM(Paths* paths);
static void RtoUNM(Paths* paths);
static void UtoRNM(Paths* paths);
static void applyOneStepNormalMode(Paths* paths);
static void sortEigen(int n, double* lambda, double** vectors);
static void momentToCartezianNM(Paths* paths);
static void buildMatrixNormalModeTrans(Paths* paths);
static void initMassesPrim(Paths* paths);
static void initMoments(Paths* paths);
static double getEKinVelocities(Paths* paths);
/*********************************************************************************************************************/
void initNormalModeTrans(Paths* paths)
{
	buildMatrixNormalModeTrans(paths);
	initMassesPrim(paths);
  	/* Positions U */
	paths->U = newCubeDouble(3,paths->nBeads,paths->molecules[0]->nAtoms);
  	RtoUNM(paths);
	/* Moment P*/
	initMoments(paths);
	paths->klass->calculateForces = calculateForcesNM;
	paths->klass->applyOneStep = applyOneStepNormalMode;
	paths->klass->momentToCartezian = momentToCartezianNM;
	paths->klass->getEKinVelocities = getEKinVelocities;
	paths->klass->calculateForces(paths);
}
/*********************************************************************************************************************/
static void initMassesPrim(Paths* paths)
{
	int i,iAtom;
	int nBeads = paths->nBeads;
	int nAtoms = paths->molecules[0]->nAtoms;

	paths->Mprim = newMatrixDouble(nBeads,nAtoms);

	i = 0;
	for(iAtom = 0;iAtom<paths->molecules[i]->nAtoms; iAtom++) 
		paths->Mprim[i][iAtom]  = paths->molecules[i]->atoms[iAtom].mass;

	for(i = 1;i<paths->nBeads; i++) 
	for(iAtom = 0;iAtom<paths->molecules[i]->nAtoms; iAtom++) 
		paths->Mprim[i][iAtom]  = paths->M[i][iAtom];
		//paths->Mprim[i][iAtom]  = paths->molecules[i]->atoms[iAtom].mass;
}
/*****************************************************************************************************************************/
static void applyOneStepNormalMode(Paths* paths)
{
  
	int iAtom,i,k;
	double*** P = paths->P;
	double*** U = paths->U;

	paths->klass->applyThermostat(paths);// must be apply 2 times, time step = dt/2 in the thermostats

	for(i = 0;i<paths->nBeads; i++) 
	for(iAtom = 0;iAtom<paths->molecules[i]->nAtoms; iAtom++) 
	for(k = 0;k<3; k++) 
	{
		P[k][i][iAtom] += 0.5*paths->F[k][i][iAtom]*paths->dt;
		U[k][i][iAtom] +=  P[k][i][iAtom]*paths->dt/paths->Mprim[i][iAtom];
	}

	UtoRNM(paths);
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
static void UtoRNM(Paths* paths)
{
	int i,j,k,iAtom;
	Molecule** mols = paths->molecules;
	int nBeads = paths->nBeads;
	int nAtoms = paths->molecules[0]->nAtoms;
	for (iAtom = 0; iAtom < nAtoms; iAtom++)
	{
    		for (i = 0; i<nBeads; i++)
		{
			for(k = 0;k<3; k++) mols[i]->atoms[iAtom].coordinates[k] = 0.0;
    			for (j = 0; j<nBeads; j++)
			for(k = 0;k<3; k++)
				mols[i]->atoms[iAtom].coordinates[k] += paths->O[i][j]*paths->U[k][j][iAtom];

		}
	}
}
/**************************************************************************************************************************/
static void RtoUNM(Paths* paths)
{
	int i,j,k,iAtom;
	Molecule** mols = paths->molecules;
	int nAtoms = paths->molecules[0]->nAtoms;
	for (iAtom = 0; iAtom < nAtoms; iAtom++)
	{
    		for (i = 0; i<paths->nBeads; i++)
		{
			for(k = 0;k<3; k++) paths->U[k][i][iAtom] = 0.0;

    			for (j = 0; j<paths->nBeads; j++)
			for(k = 0;k<3; k++)
				paths->U[k][i][iAtom] += paths->O[j][i]*mols[j]->atoms[iAtom].coordinates[k] ;

			for(k = 0;k<3; k++) paths->U[k][i][iAtom] *= paths->oneOvernBeads;
		}
	}
/*
	// print O 
	printf("O = \n");
	for(i = 0;i<paths->nBeads; i++) 
	{
	for(j = 0;j<paths->nBeads; j++) 
	{
      		printf("%f ",paths->O[i][j]);
	}
	printf("\n");
	}
	printf("U = \n");
	printCubeDouble( paths->U, 3, paths->nBeads, nAtoms);
*/
}

/*********************************************************************************************************************/
static void calculateForcesNM(Paths* paths)
{
	int i,j,k,iAtom;
	int nAtoms = paths->molecules[0]->nAtoms;
	Molecule** mols = paths->molecules;

	paths->klass->calculateGradient(paths);

	for (iAtom = 0; iAtom < nAtoms; iAtom++)
	{
    		for (i = 0; i <paths->nBeads; i++)
		{
			for(k=0;k<3;k++) paths->F[k][i][iAtom] = 0;

    			for (j = 0; j <paths->nBeads; j++)
			for(k=0;k<3;k++) paths->F[k][i][iAtom] -= mols[j]->atoms[iAtom].gradient[k]*paths->O[j][i];

			for(k=0;k<3;k++) paths->F[k][i][iAtom] *= paths->oneOvernBeads;

			for(k=0;k<3;k++) paths->F[k][i][iAtom] -= paths->M[i][iAtom]*paths->wp2*paths->U[k][i][iAtom];
		}
	}
}
/*****************************************************************************/
static void sortEigen(int n, double* lambda, double** vectors)
{
	int i;
	int j;
	int k;
	double dum;
	if(n<1 || !lambda || !vectors) return;
	for(i=0;i<n;i++)
	{
		k = i;
		for(j=i+1;j<n;j++)
			if(lambda[j]<lambda[k]) k = j;
		if(k==i) continue;
		/* swap i and k vectors */
		dum = lambda[i];
		lambda[i] = lambda[k];
		lambda[k] = dum;
		for(j=0;j<n;j++)
		{
			dum =  vectors[j][i];
			vectors[j][i] = vectors[j][k];
			vectors[j][k] = dum;
		}
	}
}
/**************************************************************************************************************************/
static void momentToCartezianNM(Paths* paths)
{
	int i,j,k,iAtom;
	Molecule** mols = paths->molecules;
	int nBeads = paths->nBeads;
	int nAtoms = paths->molecules[0]->nAtoms;
	double*** P = paths->P;

	for (iAtom = 0; iAtom < nAtoms; iAtom++)
	{
    		for (i = 0; i<nBeads; i++)
		{
			for(k = 0;k<3; k++) mols[i]->atoms[iAtom].velocity[k] = 0.0;
    			for (j = 0; j<nBeads; j++)
			for(k = 0;k<3; k++)
				mols[i]->atoms[iAtom].velocity[k] += paths->O[i][j]*P[k][j][iAtom];
		}
	}

    	for (i = 0; i<paths->nBeads; i++)
	for (iAtom = 0; iAtom < nAtoms; iAtom++)
		for(k = 0;k<3; k++) mols[i]->atoms[iAtom].velocity[k] /= paths->Mprim[i][iAtom];
}
/**************************************************************************************************************************/
static void initMoments(Paths* paths)
{
	int i,j,k,iAtom;
	//Molecule** mols = paths->molecules;
	int nBeads = paths->nBeads;
	int nAtoms = paths->molecules[0]->nAtoms;

	double*** Q = newCubeDouble(3,paths->nBeads,nAtoms);
	paths->P = newCubeDouble(3,paths->nBeads,nAtoms);
	
	//printf("begin carttonm\n");
	for(k=0;k<3;k++)
	{
		for(i = 0;i<paths->nBeads; i++) 
		for(iAtom = 0;iAtom<paths->molecules[i]->nAtoms; iAtom++) 
			Q[k][i][iAtom]  = paths->molecules[i]->atoms[iAtom].velocity[k]*paths->Mprim[i][iAtom];
	}
	//printf("end carttonm 1\n");
	for(i = 0; i<nBeads ;  i++) 
	{
		for(iAtom = 0;iAtom<paths->molecules[i]->nAtoms; iAtom++) 
		for(k=0;k<3;k++) paths->P[k][i][iAtom]=  0.0;

		for(j = 0; j<nBeads ;  j++) 
		for(iAtom = 0;iAtom<paths->molecules[i]->nAtoms; iAtom++) 
		{
			for(k=0;k<3;k++)
			{
				 paths->P[k][i][iAtom] +=
      				paths->O[i][j]* paths->molecules[i]->atoms[iAtom].velocity[k]*paths->Mprim[i][iAtom];
			}
		}
		for(iAtom = 0;iAtom<paths->molecules[i]->nAtoms; iAtom++) 
		for(k=0;k<3;k++)
			paths->P[k][i][iAtom] *= paths->oneOvernBeads;
/*
		for(iAtom = 0;iAtom<paths->molecules[i]->nAtoms; iAtom++) 
		for(k=0;k<3;k++)
				paths->P[k][i][iAtom] =  0.0;
*/
	}
	freeCubeDouble(&Q,3, paths->nBeads);
}
/*********************************************************************************************************************/
static void buildMatrixNormalModeTrans(Paths* paths)
{
	int i,j,k,iAtom;
	int nAtoms = paths->molecules[0]->nAtoms;
	double* A;
	double* lambda;
	int nBeads = paths->nBeads;
	int n2Beads = paths->nBeads*(paths->nBeads+1)/2;
	double sqrtnBeads = sqrt(1.0*nBeads);

	A = newVectorDouble(n2Beads);
	k = 0;
	for(i = 0;i<nBeads; i++) 
	for(j = 0;j<=i; j++) A[k++] = 0;

	k = 0;
	for(i = 0;i<nBeads; i++) 
	for(j = 0;j<=i; j++) 
	{
      		if (i==j) A[k] = 2;
      		else if (i==j+1) A[k] = -1;
		k++;
	}
	k = n2Beads-nBeads;
	if(k!=0) A[k] = -1; // don't change it if nBeads = 1

	// print A 
	printf("A matrix for normal mode transformation\n");
	k = 0;
	for(i = 0;i<nBeads; i++) 
	{
	for(j = 0;j<=i; j++) 
	{
      		printf("%f ",A[k] );
		k++;
	}
	printf("\n");
	}
	lambda = newVectorDouble(nBeads);
	paths->O = newMatrixDouble(nBeads,nBeads);
	eigenQL(nBeads, A, lambda, paths->O);
	sortEigen(nBeads, lambda, paths->O);
	for(i = 0;i<nBeads; i++) lambda[i] *=nBeads;
	k = 0;
	printf("lambda\n");
	for(i = 0;i<nBeads; i++) 
      		printf("%f ",lambda[i] );
	printf("\n");
	printf("End eigenQL\n");
	freeVectorDouble(&A);
	if(paths->O[0][0]<0) sqrtnBeads = -sqrtnBeads;
	for(i = 0;i<nBeads; i++) 
	for(j = 0;j<nBeads; j++) 
		paths->O[i][j] *= sqrtnBeads;

	// print O 
	printf("Transormation matrix in normal mode\n");
	for(i = 0;i<nBeads; i++) 
	{
	for(j = 0;j<nBeads; j++) 
	{
      		printf("%f ",paths->O[i][j]);
	}
	printf("\n");
	}
	paths->M  = newMatrixDouble(nBeads,nAtoms);
	for(i = 0;i<paths->nBeads; i++) 
	for(iAtom = 0;iAtom<paths->molecules[i]->nAtoms; iAtom++) 
		paths->M[i][iAtom]  = fabs(lambda[i]*paths->molecules[i]->atoms[iAtom].mass);
	
	if(paths->nBeads==1)
		for(iAtom = 0;iAtom<paths->molecules[0]->nAtoms; iAtom++) paths->M[0][iAtom]  = 0;
	else
	if(paths->nBeads==2)
	{
		for(iAtom = 0;iAtom<paths->molecules[0]->nAtoms; iAtom++) paths->M[0][iAtom]  = 0;
		for(iAtom = 0;iAtom<paths->molecules[1]->nAtoms; iAtom++) paths->M[1][iAtom]  = 4*2;
		/*
		paths->O[0][0] = 1;
		paths->O[1][0] = 1;
		paths->O[0][1] = 1;
		paths->O[1][1] = -1;
		*/
	}

	freeVectorDouble(&lambda);
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
