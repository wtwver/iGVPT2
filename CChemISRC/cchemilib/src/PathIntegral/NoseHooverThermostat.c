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

static void applyNoseHooverThermostat(Paths* paths);
static void evolveFirstLast(Paths* paths, int iSY);
static void evolveLastFirst(Paths* paths, int iSY);
static void evolveParticleMoments(Paths* paths, int iSY);
static void evolvePositions(Paths* paths, int iSY);
static void initSY(Paths* paths);
static void initMassesNH(Paths* paths);
static void initNHForcesMoments(Paths* paths);
/*********************************************************************************************************************/
void initNoseHooverThermostat(Paths* paths)
{
	initMassesNH(paths);
	initNHForcesMoments(paths);
	initSY(paths);
	paths->klass->applyThermostat = applyNoseHooverThermostat;
	paths->xNHTotal = 0;
	printf("End nosee hoover\n");
}
/*****************************************************************************************************************************/
static void applyNoseHooverThermostat(Paths* paths)
{
	int iNHSteps,iSY;

	for (iNHSteps = 0; iNHSteps < paths->nNHSteps; iNHSteps++)
    	for (iSY = 0; iSY < paths->nSY; iSY += 1)
	{
		evolveLastFirst(paths, iSY);
		evolveParticleMoments(paths, iSY);
		evolvePositions(paths, iSY);
		evolveFirstLast(paths, iSY);
	}
}
/*********************************************************************************************************************/
static void initNHForcesMoments(Paths* paths)
{
	int nAtoms = paths->molecules[0]->nAtoms;
	int nBeads = paths->nBeads;
	paths->NHP = newCubeDouble(nBeads,3*nAtoms,paths->nNH);
	initCubeDouble(paths->NHP, nBeads,3*nAtoms,paths->nNH,0.0);
	paths->NHF = newCubeDouble(nBeads,3*nAtoms,paths->nNH);
	initCubeDouble(paths->NHF, nBeads,3*nAtoms,paths->nNH,0.0);
}
/*********************************************************************************************************************/
static void initMassesNH(Paths* paths)
{
	int i,iAtom;
	int nAtoms = paths->molecules[0]->nAtoms;
	int nBeads = paths->nBeads;

// Eq. 12.6.14, Tuckerman book
	paths->oneOverNHM = newVectorDouble(paths->nNH);
	for (i = 0; i < paths->nNH; i++) paths->oneOverNHM[i] = paths->wp2/paths->kT;

	paths->oneOverM = newMatrixDouble(nBeads,nAtoms);
	for (i = 0; i < paths->nBeads; i++) 
	for (iAtom = 0; iAtom < nAtoms; iAtom++) paths->oneOverM[i][iAtom] = 1.0/paths->Mprim[i][iAtom];
}
/*********************************************************************************************************************/
static void initSY(Paths* paths)
{
// Eqs. 4.11.xx, Tuckerman book
	int i;
	if (paths->nSY == 1)
	{
		paths->NHd = newVectorDouble(paths->nSY);
    		paths->NHd[0] = 1.0;
	} else if (paths->nSY == 3)
	{
		paths->NHd = newVectorDouble(paths->nSY);
    		paths->NHd[0] = 1.0/(2.0 - pow(2.0,1.0/3.0));
    		paths->NHd[2] = paths->NHd[0];
    		paths->NHd[1] = 1.0 - (paths->NHd[0] + paths->NHd[2]);
	} else if (paths->nSY == 5)
	{
		double bi = 1/(2.0*2-1.0);
		double dum = 1.0/(4.0-pow(4.0,bi));
		paths->NHd = newVectorDouble(paths->nSY);
		for(i=0;i<paths->nSY;i++) paths->NHd[i] = dum;
		paths->NHd[2] = 1.0-4.0*dum;
	} else if (paths->nSY == 7)
	{
		int nc;
		nc = paths->nSY/2;
		paths->NHd = newVectorDouble(paths->nSY);
    		paths->NHd[0] = 0.784513610477560;
    		paths->NHd[1] = 0.235573213359357;
    		paths->NHd[2] = -1.17767998417887;
		for(i=0;i<nc;i++) paths->NHd[paths->nSY-1-i] = paths->NHd[i];
    		paths->NHd[nc] = 1;
		for(i=0;i<nc;i++) paths->NHd[nc] -= paths->NHd[i]+paths->NHd[paths->nSY-1-i];
	} else if (paths->nSY == 9)
	{
		int nc;
		nc = paths->nSY/2;
		paths->NHd = newVectorDouble(paths->nSY);
    		paths->NHd[0] = 0.192;
    		paths->NHd[1] = 0.554910818409783619692725006662999;
    		paths->NHd[2] = 0.124659619941888644216504240951585;
    		paths->NHd[3] = -0.843182063596933505315033808282941; 
		for(i=0;i<nc;i++) paths->NHd[paths->nSY-1-i] = paths->NHd[i];
    		paths->NHd[nc] = 1;
		for(i=0;i<nc;i++) paths->NHd[nc] -= paths->NHd[i]+paths->NHd[paths->nSY-1-i];
	}
	else
	{
		printf("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n");
		printf("Suzuki-Yoshida Order must be 1,3,5,7 or 9!\n");
		printf("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n");
		fflush(stdout);
    		exit(1);
	}
  	for (i = 0; i < paths->nSY; i++) paths->NHd[i] *= paths->dt/(paths->nNHSteps*1.0);
}
// Eqs 12.6.14 && 4.11.17 &&  4.11.18 Tuckerman book
/*****************************************************************************************************************************/
static void NHPScale(Paths* paths, int iNH, double f)
{
	int i,j,iAtom;
	double*** NHP = paths->NHP;
      	double* oneOverNHM = paths->oneOverNHM;

	for(i = 0;i<paths->nBeads; i++) 
	for(iAtom = 0;iAtom<paths->molecules[i]->nAtoms; iAtom++) 
	for(j = 0;j<3; j++) 
      		NHP[i][3*iAtom+j][iNH] *=  exp(f * NHP[i][3*iAtom+j][iNH+1] * oneOverNHM[iNH+1]);
}
/*****************************************************************************************************************************/
static void NHPShift(Paths* paths, int iNH, double f)
{
	int i,j,iAtom;
	double*** NHP = paths->NHP;
	double*** NHF = paths->NHF;

	for(i = 0;i<paths->nBeads; i++) 
	for(iAtom = 0;iAtom<paths->molecules[i]->nAtoms; iAtom++) 
	for(j = 0;j<3; j++) 
      		NHP[i][3*iAtom+j][iNH] += f * NHF[i][3*iAtom+j][iNH];
}
/*****************************************************************************************************************************/
static void NHFSet(Paths* paths, int iNH)
{
	int i,j,iAtom;
	double*** NHP = paths->NHP;
	double*** NHF = paths->NHF;
      	double* oneOverNHM = paths->oneOverNHM;
	double*** P = paths->P;

	if(iNH!=0)
	{
		for(i = 0;i<paths->nBeads; i++) 
			for(iAtom = 0;iAtom<paths->molecules[i]->nAtoms; iAtom++) 
			for(j = 0;j<3; j++) 
      				NHF[i][3*iAtom+j][iNH] = NHP[i][3*iAtom+j][iNH-1] * NHP[i][3*iAtom+j][iNH-1] * oneOverNHM[iNH-1] -paths->kT;
	}
	else
	{
		for(i = 0;i<paths->nBeads; i++) 
		for(iAtom = 0;iAtom<paths->molecules[i]->nAtoms; iAtom++) 
		for(j = 0;j<3; j++) 
		{
      			NHF[i][3*iAtom+j][iNH] =  P[j][i][iAtom]*P[j][i][iAtom]*paths->oneOverM[i][iAtom]- paths->kT;
		}
	}
}
/*****************************************************************************************************************************/
static void evolveLastFirst(Paths* paths, int iSY)
{
	int iNH;
	int nNH = paths->nNH;
      	double* NHd = paths->NHd;
	double fscale = -0.125 * NHd[iSY];
	double fshift = 0.25 * NHd[iSY];

	/* Evolve the last therm moment in each chain */
	iNH = nNH-1;
	NHFSet(paths, iNH);
	NHPShift(paths, iNH,  fshift);
      
	/* Evolve the last-1 to the first thermo velocitiy in each chain */
	for (iNH = nNH-2; iNH >= 0; iNH--)
	{
		NHPScale(paths, iNH, fscale);
		NHFSet(paths, iNH);
		NHPShift(paths, iNH,  fshift);
		NHPScale(paths, iNH, fscale);
	}
}
/*****************************************************************************************************************************/
static void evolveFirstLast(Paths* paths, int iSY)
{
	int iNH;
	int nNH = paths->nNH;
      	double* NHd = paths->NHd;
	double fscale = -0.125 * NHd[iSY];
	double fshift = 0.25 * NHd[iSY];

	/* Evolve the 1 to last-1 therm moment in each chain  calculting therm forces as you go along */
	/* 0 => last-1 */
	for (iNH = 0; iNH<paths->nNH-1; iNH++)
	{
		NHPScale(paths, iNH, fscale);
		NHFSet(paths, iNH);
		NHPShift(paths, iNH,  fshift);
		NHPScale(paths, iNH, fscale);
	}
      
	/* last term */
	iNH = nNH-1;
	NHFSet(paths, iNH);
	NHPShift(paths, iNH,  fshift);
}
/*****************************************************************************************************************************/
static void evolveParticleMoments(Paths* paths, int iSY)
{
	int i,j,iAtom;
	double*** NHP = paths->NHP;
      	double* oneOverNHM = paths->oneOverNHM;
      	double* NHd = paths->NHd;
	double*** P = paths->P;
	double f = -0.5 * NHd[iSY]*oneOverNHM[0];

	/* Evolve the particle velocities (by the scaling factor) */
	for(i = 0;i<paths->nBeads; i++) 
	for(iAtom = 0;iAtom<paths->molecules[i]->nAtoms; iAtom++) 
		for(j = 0;j<3; j++) 
			P[j][i][iAtom] *=  exp(f* NHP[i][3*iAtom+j][0]);

}
/*****************************************************************************************************************************/
static void evolvePositions(Paths* paths, int iSY)
{
	int i,j,iNH,iAtom;
	double*** NHP = paths->NHP;
      	double* oneOverNHM = paths->oneOverNHM;
	// double f =  -paths->NHd[iSY]*0.5; // negative from eq 4.11.17 !
	double f =  paths->NHd[iSY]*0.5;

	for (iNH = 0; iNH <paths->nNH; iNH++)
	for(i = 0;i<paths->nBeads; i++) 
	for(iAtom = 0;iAtom<paths->molecules[i]->nAtoms; iAtom++) 
		for(j = 0;j<3; j++) 
		paths->xNHTotal += NHP[i][3*iAtom+j][iNH]*oneOverNHM[iNH]*f;

}
/*****************************************************************************************************************************/
double getNHEnergy(Paths* paths)
{
	int i,j,iAtom;
	double*** NHP = paths->NHP;
      	double* oneOverNHM = paths->oneOverNHM;
	int iNH;

	double energy  = 0;
	for (iNH = 0; iNH <paths->nNH; iNH++)
		for(i = 0;i<paths->nBeads; i++) 
			for(iAtom = 0;iAtom<paths->molecules[i]->nAtoms; iAtom++) 
			for(j = 0;j<3; j++) 
      				energy += NHP[i][3*iAtom+j][iNH-1] * NHP[i][3*iAtom+j][iNH-1] * oneOverNHM[iNH-1];
	energy *=0.5;
	// add potential energy
	energy += paths->xNHTotal* paths->kT;
	return energy;
}
