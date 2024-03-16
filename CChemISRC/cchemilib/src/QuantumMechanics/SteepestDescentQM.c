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

/* SteepestDescentQM.c  */
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>

#include "../Utils/Utils.h"
#include "../Utils/AtomsProp.h"
#include "../Utils/Constants.h"
#include "../Utils/Types.h"
#include "../Molecule/Molecule.h"
#include "../QuantumMechanics/SteepestDescentQM.h"

static void Minimize(SteepestDescentQM* steepestDescentQM);
/**********************************************************************/
void	runSteepestDescentQM(
		SteepestDescentQM* steepestDescentQM, QuantumMechanicsModel* qmModel, 
		int updateFrequency, int maxIterations, double epsilon,
		int maxLines
		)
{

	steepestDescentQM->qmModel = qmModel;
	steepestDescentQM->numberOfAtoms = qmModel->molecule.nAtoms;
	steepestDescentQM->updateFrequency = updateFrequency;
	steepestDescentQM->maxIterations = maxIterations;
	steepestDescentQM->updateNumber = 0;
	steepestDescentQM->epsilon = epsilon;
	steepestDescentQM->rmsDeplacment = 0;
	steepestDescentQM->maxDeplacment = 0;
	steepestDescentQM->gradientNorm = 0;
	steepestDescentQM->maxLines=maxLines;

	Minimize(steepestDescentQM);
}
/**********************************************************************/
void	freeSteepestDescentQM(SteepestDescentQM* steepestDescentQM)
{

	steepestDescentQM->qmModel = NULL;
	steepestDescentQM->numberOfAtoms = 0;
	steepestDescentQM->updateFrequency = 0;
	steepestDescentQM->maxIterations = 0;
	steepestDescentQM->updateNumber = 0;
	steepestDescentQM->epsilon = 0;
	steepestDescentQM->rmsDeplacment = 0;
	steepestDescentQM->maxDeplacment = 0;
	steepestDescentQM->gradientNorm = 0;
	steepestDescentQM->maxLines=0;

}
/**********************************************************************/
static void Minimize(SteepestDescentQM* steepestDescentQM)
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
	FILE* fileOut = steepestDescentQM->logfile;

	steepestDescentQM->updateNumber = 0;

	steepestDescentQM->qmModel->klass->calculateGradient(steepestDescentQM->qmModel);

	steepestDescentQM->gradientNorm = 0;
	for (  i = 0; i < steepestDescentQM->numberOfAtoms; i++ )
		for(j=0;j<3;j++)
			steepestDescentQM->gradientNorm += 
				steepestDescentQM->qmModel->molecule.atoms[i].gradient[j] *
				steepestDescentQM->qmModel->molecule.atoms[i].gradient[j]; 

	lastGradientNorm = sqrt( steepestDescentQM->gradientNorm );

	fprintf(fileOut,"=============================================================================================================================================\n");
        fprintf(fileOut,"\t\t\tOptimization by Steepest Descent method\n");
        fprintf(fileOut,"---------------------------------------------------------------------------------------------------------------------------------------------\n");
        fprintf(fileOut,"\t\t\tMaxIteration \t\t= %d\n",steepestDescentQM->maxIterations);
        fprintf(fileOut,"\t\t\tEpsilon \t\t= %0.4e\n",steepestDescentQM->epsilon);
        fprintf(fileOut,"\t\t\tMax lines search \t= %d\n",steepestDescentQM->maxLines);
        fprintf(fileOut,"=============================================================================================================================================\n");
        fflush(fileOut); fflush(stderr);


	while( 
			( lastGradientNorm > steepestDescentQM->epsilon*steepestDescentQM->epsilon) && 
			( iteration++ < steepestDescentQM->maxIterations )
	     )
	{

		steepestDescentQM->qmModel->klass->calculateGradient(steepestDescentQM->qmModel);
		steepestDescentQM->gradientNorm = 0;
		for (  i = 0; i < steepestDescentQM->numberOfAtoms; i++ )
			for(j=0;j<3;j++)
				steepestDescentQM->gradientNorm += 
					steepestDescentQM->qmModel->molecule.atoms[i].gradient[j] *
					steepestDescentQM->qmModel->molecule.atoms[i].gradient[j]; 

		steepestDescentQM->gradientNorm = sqrt( steepestDescentQM->gradientNorm );
		
	
		if(steepestDescentQM->gradientNorm<1e-12) break;

		
        	steepestDescentQM->qmModel->klass->calculateEnergy(steepestDescentQM->qmModel);
		f0 = steepestDescentQM->qmModel->molecule.potentialEnergy;
		term = 0;
		fg = 1.0;
		if(steepestDescentQM->gradientNorm>1) fg = 1.0/steepestDescentQM->gradientNorm;

		for(ii=steepestDescentQM->maxLines;ii>=1;ii--)
		{
			term = ii*0.0001;
			for (  i = 0; i < steepestDescentQM->numberOfAtoms; i++ )
			{
				for(j=0;j<3;j++)
					steepestDescentQM->qmModel->molecule.atoms[i].coordinates[j]-=
					fg*term*steepestDescentQM->qmModel->molecule.atoms[i].gradient[j]; 
			}

        		steepestDescentQM->qmModel->klass->calculateEnergy(steepestDescentQM->qmModel);
			f1 = steepestDescentQM->qmModel->molecule.potentialEnergy;
			if(f1<f0)
				break;
			for (  i = 0; i < steepestDescentQM->numberOfAtoms; i++ )
			{
				for(j=0;j<3;j++)
					steepestDescentQM->qmModel->molecule.atoms[i].coordinates[ j ] += 
					fg*term*steepestDescentQM->qmModel->molecule.atoms[i].gradient[j]; 
			}
		}
		if(ii<=1) 
		{
			fprintf(steepestDescentQM->logfile, ("!!!!!!!!!!Steep descent Max lines exessed\n"));
			break;
		}
		lastGradientNorm = steepestDescentQM->gradientNorm;

		if ( steepestDescentQM->updateNumber++ >= steepestDescentQM->updateFrequency )
		{
			free(str);
			str = strdup_printf(("Gradient = %0.10f Energy = %0.10f"),(double)steepestDescentQM->gradientNorm, (double) f0); 
			/* redrawMolecule(&steepestDescentQM->qmModel->molecule,str);*/
			fprintf(steepestDescentQM->logfile,"%s\n",str);
			fflush(steepestDescentQM->logfile);
			steepestDescentQM->updateNumber = 0;
		}
	}

	steepestDescentQM->qmModel->klass->calculateGradient(steepestDescentQM->qmModel);
	steepestDescentQM->gradientNorm = 0;
	for (  i = 0; i < steepestDescentQM->numberOfAtoms; i++ )
		for(j=0;j<3;j++)
			steepestDescentQM->gradientNorm += 
				steepestDescentQM->qmModel->molecule.atoms[i].gradient[j] *
				steepestDescentQM->qmModel->molecule.atoms[i].gradient[j]; 

	steepestDescentQM->gradientNorm = sqrt( steepestDescentQM->gradientNorm );
       	steepestDescentQM->qmModel->klass->calculateEnergy(steepestDescentQM->qmModel);
	energy = steepestDescentQM->qmModel->molecule.potentialEnergy;
	free(str);
	str = strdup_printf(("Gradient = %0.10f  Energy = %0.10f (Kcal/mol)"),
			(double)steepestDescentQM->gradientNorm,(double)energy); 

	/* redrawMolecule(&steepestDescentQM->qmModel->molecule,str);*/
	fprintf(steepestDescentQM->logfile,"%s\n",str);
	fflush(steepestDescentQM->logfile);
	free(str);
}
/********************************************************************************/
