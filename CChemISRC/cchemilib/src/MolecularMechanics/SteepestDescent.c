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

/* SteepestDescent.c  */
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>

#include "../Utils/Utils.h"
#include "../Utils/AtomsProp.h"
#include "../Utils/Constants.h"
#include "../Utils/Types.h"
#include "../Molecule/Molecule.h"
#include "../MolecularMechanics/SteepestDescent.h"

static void Minimize(SteepestDescent* steepestDescent);
/**********************************************************************/
void	runSteepestDescent(
		SteepestDescent* steepestDescent, ForceField* forceField, 
		int updateFrequency, int maxIterations, double epsilon,
		int maxLines
		)
{

	steepestDescent->forceField = forceField;
	steepestDescent->numberOfAtoms = forceField->molecule.nAtoms;
	steepestDescent->updateFrequency = updateFrequency;
	steepestDescent->maxIterations = maxIterations;
	steepestDescent->updateNumber = 0;
	steepestDescent->epsilon = epsilon;
	steepestDescent->rmsDeplacment = 0;
	steepestDescent->maxDeplacment = 0;
	steepestDescent->gradientNorm = 0;
	steepestDescent->maxLines=maxLines;

	Minimize(steepestDescent);
}
/**********************************************************************/
void	freeSteepestDescent(SteepestDescent* steepestDescent)
{

	steepestDescent->forceField = NULL;
	steepestDescent->numberOfAtoms = 0;
	steepestDescent->updateFrequency = 0;
	steepestDescent->maxIterations = 0;
	steepestDescent->updateNumber = 0;
	steepestDescent->epsilon = 0;
	steepestDescent->rmsDeplacment = 0;
	steepestDescent->maxDeplacment = 0;
	steepestDescent->gradientNorm = 0;
	steepestDescent->maxLines=0;

}
/**********************************************************************/
static void Minimize(SteepestDescent* steepestDescent)
{
	double energy;
	int iteration = 0;
	double lastGradientNorm = 1;
	double term = 1;
	char* str = strdup(" ");
	int i;
	int j;
	double f0,f1;
	int ii;
	double fg;

	steepestDescent->updateNumber = 0;

	steepestDescent->forceField->klass->calculateGradient(steepestDescent->forceField);

	steepestDescent->gradientNorm = 0;
	for (  i = 0; i < steepestDescent->numberOfAtoms; i++ )
		for(j=0;j<3;j++)
			steepestDescent->gradientNorm += 
				steepestDescent->forceField->molecule.atoms[i].gradient[j] *
				steepestDescent->forceField->molecule.atoms[i].gradient[j]; 

	lastGradientNorm = sqrt( steepestDescent->gradientNorm );

	while( 
			( lastGradientNorm > steepestDescent->epsilon ) && 
			( iteration++ < steepestDescent->maxIterations )
	     )
	{

		steepestDescent->forceField->klass->calculateGradient(steepestDescent->forceField);
		steepestDescent->gradientNorm = 0;
		for (  i = 0; i < steepestDescent->numberOfAtoms; i++ )
			for(j=0;j<3;j++)
				steepestDescent->gradientNorm += 
					steepestDescent->forceField->molecule.atoms[i].gradient[j] *
					steepestDescent->forceField->molecule.atoms[i].gradient[j]; 

		steepestDescent->gradientNorm = sqrt( steepestDescent->gradientNorm );
		
	
		if(steepestDescent->gradientNorm<1e-12)
			break;

        	f0 = steepestDescent->forceField->klass->calculateEnergyTmp(
		      steepestDescent->forceField,&steepestDescent->forceField->molecule);
		term = 0;
		fg = 1.0;
		if(steepestDescent->gradientNorm>1)
			fg = 1.0/steepestDescent->gradientNorm;

		for(ii=steepestDescent->maxLines;ii>=1;ii--)
		{
			term = ii*0.01;
			for (  i = 0; i < steepestDescent->numberOfAtoms; i++ )
			{
				for(j=0;j<3;j++)
					steepestDescent->forceField->molecule.atoms[i].coordinates[j]-=
					fg*term*steepestDescent->forceField->molecule.atoms[i].gradient[j]; 
			}
			updateGeometryCL(steepestDescent->forceField,NULL);

        		f1 = steepestDescent->forceField->klass->calculateEnergyTmp(
		      	steepestDescent->forceField,&steepestDescent->forceField->molecule);
			if(f1<f0)
				break;
			for (  i = 0; i < steepestDescent->numberOfAtoms; i++ )
			{
				for(j=0;j<3;j++)
					steepestDescent->forceField->molecule.atoms[i].coordinates[ j ] += 
					fg*term*steepestDescent->forceField->molecule.atoms[i].gradient[j]; 
			}
			updateGeometryCL(steepestDescent->forceField,NULL);
		}
		if(ii<=1)
			break;
		lastGradientNorm = steepestDescent->gradientNorm;

		if ( steepestDescent->updateNumber++ >= steepestDescent->updateFrequency )
		{
			free(str);
			str = strdup_printf(("Gradient = %f "),(double)steepestDescent->gradientNorm); 
			/* redrawMolecule(&steepestDescent->forceField->molecule,str);*/
			fprintf(steepestDescent->logfile,"%s\n",str);
			fflush(steepestDescent->logfile);
			steepestDescent->updateNumber = 0;
		}
	}

	steepestDescent->forceField->klass->calculateGradient(steepestDescent->forceField);
	steepestDescent->gradientNorm = 0;
	for (  i = 0; i < steepestDescent->numberOfAtoms; i++ )
		for(j=0;j<3;j++)
			steepestDescent->gradientNorm += 
				steepestDescent->forceField->molecule.atoms[i].gradient[j] *
				steepestDescent->forceField->molecule.atoms[i].gradient[j]; 

	steepestDescent->gradientNorm = sqrt( steepestDescent->gradientNorm );
        energy = steepestDescent->forceField->klass->calculateEnergyTmp(
		      steepestDescent->forceField,&steepestDescent->forceField->molecule);
	free(str);
	str = strdup_printf(("Gradient = %f  Energy = %f (Kcal/mol)"),
			(double)steepestDescent->gradientNorm,(double)energy); 

	/* redrawMolecule(&steepestDescent->forceField->molecule,str);*/
	fprintf(steepestDescent->logfile,"%s\n",str);
	fflush(steepestDescent->logfile);
	free(str);
}
/********************************************************************************/
