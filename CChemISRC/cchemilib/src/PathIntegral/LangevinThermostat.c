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
Ceriotti et al., THE JOURNAL OF CHEMICAL PHYSICS 133, 124104 2010
*/
static void applyLangevinThermostat(Paths* paths);
/*********************************************************************************************************************/
void initLangevinThermostat(Paths* paths)
{
	int i,iAtom;
	int nAtoms = paths->molecules[0]->nAtoms;
	int nBeads = paths->nBeads;
	double sqrtnBeads =sqrt(1.0*nBeads);
	double* omegas; /* Normal mode frequencies */
	double* gammas; /* Langevin friction constants */
	if(paths->friction<=0) paths->friction = 1;

	omegas = newVectorDouble(nBeads);
	for(i = 0;i<nBeads; i++) omegas[i] = 2*paths->wp*sqrtnBeads*sin(i*M_PI/nBeads);

	gammas = newVectorDouble(nBeads);
	for(i = 0;i<nBeads; i++) gammas[i] =  paths->friction; /* 2*omegas[i] to be tested */

	paths->LTPMult  = newVectorDouble(nBeads);
	for(i = 0;i<nBeads; i++) paths->LTPMult[i] = exp(-paths->dt/2.0*gammas[i]); 
	//printf("alpha = %f\n",-paths->dt/2.0*gammas[0]); 

	paths->LTFMult  = newMatrixDouble(nBeads,nAtoms);
	for(i = 0;i<nBeads; i++) 
	for(iAtom = 0;iAtom<nAtoms; iAtom++) 
			paths->LTFMult[i][iAtom] = sqrt(paths->Mprim[i][iAtom]*paths->kT)*sqrt(1.0-paths->LTPMult[i]*paths->LTPMult[i]);

	freeVectorDouble(&omegas);
	freeVectorDouble(&gammas);
	paths->klass->applyThermostat = applyLangevinThermostat;
}
/*****************************************************************************************************************************/
static void applyLangevinThermostat(Paths* paths)
{
	int i,k,iAtom;
	int nAtoms = paths->molecules[0]->nAtoms;
	double*** P = paths->P;
	for(i = 0;i<paths->nBeads; i++) 
	{
	for(iAtom = 0;iAtom<nAtoms; iAtom++) 
	for(k = 0;k<3; k++) 
	{
		P[k][i][iAtom] =  P[k][i][iAtom]*paths->LTPMult[i] + paths->LTFMult[i][iAtom]*normal();
	}
//	paths->molecules[i]->klass->removeTranslationAndRotation(paths->molecules[i]);
	}
}

/**************************************************************************************************************************/
