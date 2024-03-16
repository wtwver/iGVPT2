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

/* ConjugateGradientQM.c */
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>

#include "../Utils/Utils.h"
#include "../Utils/AtomsProp.h"
#include "../Utils/Constants.h"
#include "../Utils/Types.h"
#include "../Molecule/Molecule.h"
#include "../QuantumMechanics/ConjugateGradientQM.h"


static double maxarg1,maxarg2;
#define FMAX(a,b) (maxarg1=(a),maxarg2=(b),(maxarg1) > (maxarg2) ?\
        (maxarg1) : (maxarg2))

static void hestenesStiefel(ConjugateGradientQM* conjugateGradientQM);
static void fletcherReeves(ConjugateGradientQM* conjugateGradientQM);
static void polakRibiere(ConjugateGradientQM* conjugateGradientQM);
static void wolfPowell(ConjugateGradientQM* conjugateGradientQM);
static double lineMinimize(ConjugateGradientQM* conjugateGradientQM);
static void bracketMinimum(ConjugateGradientQM* conjugateGradientQM, double pointA[], double pointB[], double pointC[] );
static double oneDimensionalEnergy(ConjugateGradientQM* conjugateGradientQM, double factor );
static double inverseParabolicInterpolation(ConjugateGradientQM* conjugateGradientQM,
		double pointa[], double mid[], double pointb[], 
		double minimum[] );
/**********************************************************************/
void	runConjugateGradientQM(ConjugateGradientQM* conjugateGradientQM, QuantumMechanicsModel* qmModel, 
		ConjugateGradientQMOptions conjugateGradientQMOptions)
{

	int i;
	int minimizerOptions = conjugateGradientQMOptions.method;
	conjugateGradientQM->qmModel = qmModel;
	conjugateGradientQM->numberOfAtoms = qmModel->molecule.nAtoms;
	conjugateGradientQM->updateFrequency = conjugateGradientQMOptions.updateFrequency;
	conjugateGradientQM->maxIterations = conjugateGradientQMOptions.maxIterations;
	conjugateGradientQM->updateNumber = 0;
	conjugateGradientQM->maxLine = conjugateGradientQMOptions.maxLines;
	conjugateGradientQM->epsilon = conjugateGradientQMOptions.gradientNorm;
	conjugateGradientQM->initialStep = conjugateGradientQMOptions.initialStep;
	conjugateGradientQM->rmsDeplacment = 0;
	conjugateGradientQM->maxDeplacment = 0;
	conjugateGradientQM->gradientNorm = 0;
	conjugateGradientQM->initialBracket = conjugateGradientQMOptions.initialStep;
	conjugateGradientQM->lastInitialBracket = 0;
	conjugateGradientQM->term = 0;
	FILE* fileOut = conjugateGradientQM->logfile;
	fprintf(fileOut,"=============================================================================================================================================\n");
        fprintf(fileOut,"\t\t\tOptimization by Conjugate Gradient method\n");
        fprintf(fileOut,"---------------------------------------------------------------------------------------------------------------------------------------------\n");
        fprintf(fileOut,"\t\t\tMaxIteration \t\t= %d\n",conjugateGradientQM->maxIterations);
        fprintf(fileOut,"\t\t\tEpsilon \t\t= %0.4e\n",conjugateGradientQM->epsilon);
        fprintf(fileOut,"\t\t\tMax lines search \t= %d\n",conjugateGradientQM->maxLine);
	switch(minimizerOptions)
	{
		case 1 : fprintf(fileOut,"\t\t\tHestenes Stiefel approach\n");break;
		case 2 : fprintf(fileOut,"\t\t\tFletcher Reeves approach\n");break;
		case 3 : fprintf(fileOut,"\t\t\tPolak Ribiere  approach\n");break;
		case 4 : fprintf(fileOut,"\t\t\tWolf Powell approach\n");break;
	}
        fprintf(fileOut,"=============================================================================================================================================\n");
        fflush(fileOut); fflush(stderr);



	for(i=0;i<3;i++)
	{
		conjugateGradientQM->lastGradient[i] = malloc(conjugateGradientQM->numberOfAtoms*sizeof(double));
		conjugateGradientQM->direction[i] = malloc(conjugateGradientQM->numberOfAtoms*sizeof(double));
	}

	switch(minimizerOptions)
	{
		case 1 : hestenesStiefel(conjugateGradientQM);break;
		case 2 : fletcherReeves(conjugateGradientQM);break;
		case 3 : polakRibiere(conjugateGradientQM);break;
		case 4 : wolfPowell(conjugateGradientQM);break;
	}
}
/**********************************************************************/
void	freeConjugateGradientQM(ConjugateGradientQM* conjugateGradientQM)
{

	int i;
	conjugateGradientQM->qmModel = NULL;
	conjugateGradientQM->numberOfAtoms = 0;
	conjugateGradientQM->updateFrequency = 0;
	conjugateGradientQM->maxIterations = 0;
	conjugateGradientQM->updateNumber = 0;
	conjugateGradientQM->maxLine = 0;
	conjugateGradientQM->epsilon = 0;
	conjugateGradientQM->initialStep = 0;
	conjugateGradientQM->rmsDeplacment = 0;
	conjugateGradientQM->maxDeplacment = 0;
	conjugateGradientQM->gradientNorm = 0;
	conjugateGradientQM->initialBracket = 0;
	conjugateGradientQM->lastInitialBracket = 0;
	conjugateGradientQM->term = 0;

	for(i=0;i<3;i++)
	{
		if(conjugateGradientQM->lastGradient[i] !=NULL)
		{
			free(conjugateGradientQM->lastGradient[i]);
			conjugateGradientQM->lastGradient[i]= NULL;
		}
		if(conjugateGradientQM->direction[i] != NULL)
		{
			free(conjugateGradientQM->direction[i]);
			conjugateGradientQM->direction[i] = NULL; 
		}
	}
}
/**********************************************************************/
static void fletcherReeves(ConjugateGradientQM* conjugateGradientQM)
{
	double lastGradientDotGradient = 0, gradientDotGradient = 0, beta;
	int iterations = 0;
	int i;
	int j;
	double energy;
	char* str = strdup(" ");

	conjugateGradientQM->qmModel->klass->calculateGradient(conjugateGradientQM->qmModel);
	for ( i = 0; i < conjugateGradientQM->numberOfAtoms; i++ )
	{
		for(j=0;j<3;j++)
		{
		conjugateGradientQM->direction[j][ i ] = 
			-conjugateGradientQM->qmModel->molecule.atoms[i].gradient[j];
		lastGradientDotGradient += conjugateGradientQM->qmModel->molecule.atoms[i].gradient[j]
					* conjugateGradientQM->qmModel->molecule.atoms[i].gradient[j]; 
		}
	}
	gradientDotGradient = lastGradientDotGradient;
	while ( 
			( lastGradientDotGradient > conjugateGradientQM->epsilon*conjugateGradientQM->epsilon ) && 
			( iterations++ < conjugateGradientQM->maxIterations ) 
		)
	{

		lineMinimize(conjugateGradientQM);
		conjugateGradientQM->qmModel->klass->calculateGradient(conjugateGradientQM->qmModel);
		lastGradientDotGradient = gradientDotGradient;	
		gradientDotGradient = 0;	
		for (  i = 0; i < conjugateGradientQM->numberOfAtoms; i++ )
			for(j=0;j<3;j++)
			{
				gradientDotGradient += 
					conjugateGradientQM->qmModel->molecule.atoms[i].gradient[j]
					* conjugateGradientQM->qmModel->molecule.atoms[i].gradient[j]; 
			}
		beta = gradientDotGradient / lastGradientDotGradient;
		for (  i = 0; i < conjugateGradientQM->numberOfAtoms; i++ )
		{
			for(j=0;j<3;j++)
			{
				conjugateGradientQM->direction[j][i] = 
					beta *  conjugateGradientQM->direction[j][i] - 
					conjugateGradientQM->qmModel->molecule.atoms[i].gradient[j];
			}
		}
		if ( conjugateGradientQM->updateNumber >= conjugateGradientQM->updateFrequency )
		{
			free(str);
			 str =  strdup_printf(("Iter # %d/%d\t Gradient(kcal/mol/Ang) = %0.14f\tEnergy(kcal/mol) = %0.14f "), iterations, conjugateGradientQM->maxIterations,sqrt(gradientDotGradient),
				conjugateGradientQM->qmModel->molecule.potentialEnergy); 
			/* redrawMolecule(&conjugateGradientQM->qmModel->molecule,str);*/
			fprintf(conjugateGradientQM->logfile, "%s\n",str);
			fflush(conjugateGradientQM->logfile);
			conjugateGradientQM->updateNumber = 0;
		}
		conjugateGradientQM->updateNumber++;
	}	

	conjugateGradientQM->qmModel->klass->calculateGradient(conjugateGradientQM->qmModel);
	conjugateGradientQM->gradientNorm = 0;
	for(  i = 0; i < conjugateGradientQM->numberOfAtoms; i++ )
	{
		for(j=0;j<3;j++)
		{
			conjugateGradientQM->gradientNorm += 
				conjugateGradientQM->qmModel->molecule.atoms[i].gradient[j] *
				conjugateGradientQM->qmModel->molecule.atoms[i].gradient[j]; 
		}
	}
	conjugateGradientQM->gradientNorm = sqrt(conjugateGradientQM->gradientNorm );
	conjugateGradientQM->qmModel->klass->calculateEnergy(conjugateGradientQM->qmModel);
	energy = conjugateGradientQM->qmModel->molecule.potentialEnergy;

	free(str);
	str = strdup_printf(("End optimisation\nGradient(kcal/mol/Ang) = %0.10f  Energy = %0.10f (kcal/mol)\n"),
			(double)conjugateGradientQM->gradientNorm,(double)energy); 

	/* redrawMolecule(&conjugateGradientQM->qmModel->molecule,str);*/
	fprintf(conjugateGradientQM->logfile, "%s\n",str);
	fflush(conjugateGradientQM->logfile);
	free(str);
}
/********************************************************************************/
static void polakRibiere(ConjugateGradientQM* conjugateGradientQM)
{
	double lastGradientDotGradient = 0, gradientDotGradient, beta;
	int iterations = 0;
	double energy;
	int i;
	int j;
	char* str = strdup(" ");

	conjugateGradientQM->qmModel->klass->calculateGradient(conjugateGradientQM->qmModel);

	for ( i = 0; i < conjugateGradientQM->numberOfAtoms; i++ )
	{
		for(j=0;j<3;j++)
		{
		conjugateGradientQM->direction[j][ i ] = 
		conjugateGradientQM->lastGradient[j][i] =
			-conjugateGradientQM->qmModel->molecule.atoms[i].gradient[j];

		lastGradientDotGradient += conjugateGradientQM->qmModel->molecule.atoms[i].gradient[j] 
					* conjugateGradientQM->qmModel->molecule.atoms[i].gradient[j]; 
		}
	}
	gradientDotGradient = lastGradientDotGradient;
	while ( 
			( lastGradientDotGradient > conjugateGradientQM->epsilon*conjugateGradientQM->epsilon ) && 
			( iterations++ < conjugateGradientQM->maxIterations ) 
		)
	{

		lineMinimize(conjugateGradientQM);
		conjugateGradientQM->qmModel->klass->calculateGradient(conjugateGradientQM->qmModel);
		lastGradientDotGradient = gradientDotGradient;
		gradientDotGradient = 0;	
		beta = 0;
		for (  i = 0; i < conjugateGradientQM->numberOfAtoms; i++ )
		{
			for (  j = 0; j < 3; j++ )
			{
				gradientDotGradient += 
					conjugateGradientQM->qmModel->molecule.atoms[i].gradient[j] *
					conjugateGradientQM->qmModel->molecule.atoms[i].gradient[j]; 
				
				beta += ( 
					conjugateGradientQM->qmModel->molecule.atoms[i].gradient[j] - 
					conjugateGradientQM->lastGradient[j][i]
					)* conjugateGradientQM->qmModel->molecule.atoms[i].gradient[j]; 
			}
		}
		beta /= lastGradientDotGradient;
		for (  i = 0; i < conjugateGradientQM->numberOfAtoms; i++ )
		{
			for (  j = 0; j < 3; j++ )
			{
				conjugateGradientQM->direction[j][i] = 
				beta * conjugateGradientQM->direction[j][i] - 
				conjugateGradientQM->qmModel->molecule.atoms[i].gradient[j];
			}
		}
		for ( i = 0; i < conjugateGradientQM->numberOfAtoms; i++ )
		{
			for(j=0;j<3;j++)
			{
				conjugateGradientQM->lastGradient[j][i] = 
				conjugateGradientQM->qmModel->molecule.atoms[i].gradient[j];
			}
		}
		if ( conjugateGradientQM->updateNumber >= conjugateGradientQM->updateFrequency )
		{
			free(str);
			/* str = strdup_printf("gradient = %f ",sqrt(gradientDotGradient)); */
			str = strdup_printf(("Iter # %d/%d\t Gradient(kcal/mol/Ang) = %0.14f\tEnergy(kcal/mol) = %0.14f "), iterations, conjugateGradientQM->maxIterations,sqrt(gradientDotGradient),
				conjugateGradientQM->qmModel->molecule.potentialEnergy); 
			/* redrawMolecule(&conjugateGradientQM->qmModel->molecule,str);*/
			fprintf(conjugateGradientQM->logfile, "%s\n",str);
			fflush(conjugateGradientQM->logfile);
			conjugateGradientQM->updateNumber = 0;
		}
		conjugateGradientQM->updateNumber++;
	}	
	conjugateGradientQM->qmModel->klass->calculateGradient(conjugateGradientQM->qmModel);
	conjugateGradientQM->gradientNorm = 0;
	for(  i = 0; i < conjugateGradientQM->numberOfAtoms; i++ )
	{
		for (  j = 0; j < 3; j++ )
		{
			conjugateGradientQM->gradientNorm += 
				conjugateGradientQM->qmModel->molecule.atoms[i].gradient[j] *
				conjugateGradientQM->qmModel->molecule.atoms[i].gradient[j]; 
		}
	}
	conjugateGradientQM->gradientNorm = sqrt(conjugateGradientQM->gradientNorm );
	conjugateGradientQM->qmModel->klass->calculateEnergy(conjugateGradientQM->qmModel);
	energy = conjugateGradientQM->qmModel->molecule.potentialEnergy;

	free(str);
	str = strdup_printf(("End optimisation\nGradient(kcal/mol/Ang) = %0.10f  Energy = %0.10f (kcal/mol)\n"),
			(double)conjugateGradientQM->gradientNorm,(double)energy); 
	/* redrawMolecule(&conjugateGradientQM->qmModel->molecule,str);*/
	fprintf(conjugateGradientQM->logfile, "%s\n",str);
	fflush(conjugateGradientQM->logfile);
	free(str);
}
/********************************************************************************/
static void hestenesStiefel(ConjugateGradientQM* conjugateGradientQM)
{
	double gradientDotGradient = 1, beta, gradientDiff, denom;
	int iterations = 0;
	double energy;
	int i;
	int j;
	char* str = strdup(" ");

	conjugateGradientQM->qmModel->klass->calculateGradient(conjugateGradientQM->qmModel);
	for ( i = 0; i < conjugateGradientQM->numberOfAtoms; i++ )
	{
		for(j=0;j<3;j++)
		{
		conjugateGradientQM->direction[j][ i ] = 
		conjugateGradientQM->lastGradient[j][i] = 
			-conjugateGradientQM->qmModel->molecule.atoms[i].gradient[j];
		}
	}
	while ( 
			( gradientDotGradient > conjugateGradientQM->epsilon*conjugateGradientQM->epsilon ) && 
			( iterations++ < conjugateGradientQM->maxIterations ) 
		)
	{

		lineMinimize(conjugateGradientQM);
		conjugateGradientQM->qmModel->klass->calculateGradient(conjugateGradientQM->qmModel);
		gradientDotGradient = 0;	
		for (  i = 0; i < conjugateGradientQM->numberOfAtoms; i++ )
		{
			for (  j = 0; j < 3; j++ )
			{
				gradientDotGradient += 
					conjugateGradientQM->qmModel->molecule.atoms[i].gradient[j]*
					conjugateGradientQM->qmModel->molecule.atoms[i].gradient[j];
			}
		}
		beta = 0;
		denom = 0;
		for (  i = 0; i < conjugateGradientQM->numberOfAtoms; i++ )
		{
			for (  j = 0; j < 3; j++ )
			{
				gradientDiff = 
					conjugateGradientQM->qmModel->molecule.atoms[i].gradient[j]-
					conjugateGradientQM->lastGradient[j][i];

				beta += 
					conjugateGradientQM->qmModel->molecule.atoms[i].gradient[j]*
				       	gradientDiff;
				denom += conjugateGradientQM->direction[j][i] * gradientDiff;
			}
		}
		if ( fabs( denom ) > 1.0e-10 )
			beta /= denom;
		else
			beta = 0;

		for (  i = 0; i < conjugateGradientQM->numberOfAtoms; i++ )
		{
			for (  j = 0; j < 3; j++ )
			{
				conjugateGradientQM->direction[j][ i ] = 
					beta * conjugateGradientQM->direction[j][ i ] - 
					conjugateGradientQM->qmModel->molecule.atoms[i].gradient[j];
			}
		}
		for ( i = 0; i < conjugateGradientQM->numberOfAtoms; i++ )
		{
			for(j=0;j<3;j++)
			{
				conjugateGradientQM->lastGradient[j][i] = 
				conjugateGradientQM->qmModel->molecule.atoms[i].gradient[j];
			}
		}
		if ( conjugateGradientQM->updateNumber >= conjugateGradientQM->updateFrequency )
		{
			free(str);
			/* str = strdup_printf(("gradient = %f "),sqrt(gradientDotGradient)); */
			str = strdup_printf(("Iter # %d/%d\t Gradient(kcal/mol/Ang) = %0.14f\tEnergy(kcal/mol) = %0.14f "), iterations, conjugateGradientQM->maxIterations,sqrt(gradientDotGradient),
				conjugateGradientQM->qmModel->molecule.potentialEnergy); 
			/* redrawMolecule(&conjugateGradientQM->qmModel->molecule,str);*/
			fprintf(conjugateGradientQM->logfile, "%s\n",str);
			fflush(conjugateGradientQM->logfile);
			conjugateGradientQM->updateNumber = 0;
		}
		conjugateGradientQM->updateNumber++;
	}	
	conjugateGradientQM->qmModel->klass->calculateGradient(conjugateGradientQM->qmModel);
	conjugateGradientQM->gradientNorm = 0;
	for(  i = 0; i < conjugateGradientQM->numberOfAtoms; i++ )
	{
		for(j=0;j<3;j++)
			conjugateGradientQM->gradientNorm += 
			conjugateGradientQM->direction[j][ i ] *
			conjugateGradientQM->direction[j][ i ];
	}
	conjugateGradientQM->gradientNorm = sqrt( conjugateGradientQM->gradientNorm );
	conjugateGradientQM->qmModel->klass->calculateEnergy(conjugateGradientQM->qmModel);
	energy = conjugateGradientQM->qmModel->molecule.potentialEnergy;

	free(str);
	str = strdup_printf(("End optimisation\nGradient(kcal/mol/Ang) = %0.10f  Energy = %0.10f (kcal/mol)\n"),
			(double)conjugateGradientQM->gradientNorm,(double)energy); 
	/* redrawMolecule(&conjugateGradientQM->qmModel->molecule,str);*/
	fprintf(conjugateGradientQM->logfile, "%s\n",str);
	fflush(conjugateGradientQM->logfile);
	free(str);
}
/********************************************************************************/
static void wolfPowell(ConjugateGradientQM* conjugateGradientQM)
{
	double lastGradientDotGradient = 0, gradientDotGradient, beta;
	int iterations = 0;
	double energy;
	int i;
	int j;
	char* str = strdup(" ");

	conjugateGradientQM->qmModel->klass->calculateGradient(conjugateGradientQM->qmModel);
	for ( i = 0; i < conjugateGradientQM->numberOfAtoms; i++ )
	{
		for(j=0;j<3;j++)
		{
		conjugateGradientQM->direction[j][ i ] = 
		conjugateGradientQM->lastGradient[j][i] =
			-conjugateGradientQM->qmModel->molecule.atoms[i].gradient[j];

		lastGradientDotGradient += conjugateGradientQM->qmModel->molecule.atoms[i].gradient[j]
					* conjugateGradientQM->qmModel->molecule.atoms[i].gradient[j]; 
		}
	}
	while ( 
			( lastGradientDotGradient > conjugateGradientQM->epsilon*conjugateGradientQM->epsilon ) && 
			( iterations++ < conjugateGradientQM->maxIterations ) 
		)
	{

		lineMinimize(conjugateGradientQM);
		conjugateGradientQM->qmModel->klass->calculateGradient(conjugateGradientQM->qmModel);
		gradientDotGradient = 0;	
		for ( i = 0; i < conjugateGradientQM->numberOfAtoms; i++ )
			for ( j = 0; j < 3; j++ )
			{
				gradientDotGradient += 
					conjugateGradientQM->qmModel->molecule.atoms[i].gradient[j]*
					conjugateGradientQM->qmModel->molecule.atoms[i].gradient[j];
			}

		beta = 0;
		for (  i = 0; i < conjugateGradientQM->numberOfAtoms; i++ )
		{
			for ( j = 0; j < 3; j++ )
				beta += ( 
					conjugateGradientQM->qmModel->molecule.atoms[i].gradient[j] - 
					conjugateGradientQM->lastGradient[j][i]
					)* conjugateGradientQM->qmModel->molecule.atoms[i].gradient[j]; 
		}
		beta /= lastGradientDotGradient;
		if ( beta < 0 )
			beta = 0;
		lastGradientDotGradient = gradientDotGradient;
		for (  i = 0; i < conjugateGradientQM->numberOfAtoms; i++ )
			for ( j = 0; j < 3; j++ )
				conjugateGradientQM->direction[j][ i ] = 
					beta * conjugateGradientQM->direction[j][ i ] - 
					conjugateGradientQM->qmModel->molecule.atoms[i].gradient[j];
		
		for ( i = 0; i < conjugateGradientQM->numberOfAtoms; i++ )
		{
			for(j=0;j<3;j++)
			{
				conjugateGradientQM->lastGradient[j][i] = 
				conjugateGradientQM->qmModel->molecule.atoms[i].gradient[j];
			}
		}
		if ( conjugateGradientQM->updateNumber >= conjugateGradientQM->updateFrequency )
		{
			free(str);
			/* str = strdup_printf(("gradient = %f "),sqrt(gradientDotGradient)); */
			str = strdup_printf(("Iter # %d/%d\t Gradient(kcal/mol/Ang) = %0.14f\tEnergy(kcal/mol) = %0.14f "), iterations, conjugateGradientQM->maxIterations,sqrt(gradientDotGradient),
				conjugateGradientQM->qmModel->molecule.potentialEnergy); 
			/* redrawMolecule(&conjugateGradientQM->qmModel->molecule,str);*/
			fprintf(conjugateGradientQM->logfile, "%s\n",str);
			fflush(conjugateGradientQM->logfile);
			conjugateGradientQM->updateNumber = 0;
		}
		conjugateGradientQM->updateNumber++;
	}	
	conjugateGradientQM->qmModel->klass->calculateGradient(conjugateGradientQM->qmModel);
	conjugateGradientQM->gradientNorm = 0;
	for(  i = 0; i < conjugateGradientQM->numberOfAtoms; i++ )
	{
		for(j=0;j<3;j++)
			conjugateGradientQM->gradientNorm += 
			conjugateGradientQM->direction[j][ i ] *
			conjugateGradientQM->direction[j][ i ];
	}
	conjugateGradientQM->gradientNorm = sqrt( conjugateGradientQM->gradientNorm );
	conjugateGradientQM->qmModel->klass->calculateEnergy(conjugateGradientQM->qmModel);
	energy = conjugateGradientQM->qmModel->molecule.potentialEnergy;

	free(str);
	str = strdup_printf(("End optimisation\nGradient(kcal/mol/Ang) = %0.10f  Energy = %0.10f (kcal/mol)\n"),
			(double)conjugateGradientQM->gradientNorm,(double)energy); 
	/* redrawMolecule(&conjugateGradientQM->qmModel->molecule,str);*/
	fprintf(conjugateGradientQM->logfile, "%s\n",str);
	fflush(conjugateGradientQM->logfile);
	free(str);
} 
/**********************************************************************/
static double lineMinimize(ConjugateGradientQM* conjugateGradientQM)
{
	double a;
	double b;
	double c;
	double minimum=0;
	double energy;
        double delta = 1.0e-7;
	int i;
	int j;
		
	a = 0; 
	b = conjugateGradientQM->initialBracket;
	bracketMinimum(conjugateGradientQM, &a, &b, &c );
	energy = inverseParabolicInterpolation(conjugateGradientQM, &a, &b, &c, &minimum );
	conjugateGradientQM->initialBracket = minimum;

	if ( 	( fabs( conjugateGradientQM->initialBracket ) < delta ) || 
		( conjugateGradientQM->initialBracket == conjugateGradientQM->lastInitialBracket ) )
	{
		conjugateGradientQM->initialBracket = 
			rand()/(double)RAND_MAX *conjugateGradientQM->initialStep;

		for (  i = 0; i < conjugateGradientQM->numberOfAtoms; i++ )
		{
			for(j=0;j<3;j++)
			conjugateGradientQM->direction[j][i] = 
				-conjugateGradientQM->qmModel->molecule.atoms[i].gradient[j];
		}
	} 
	conjugateGradientQM->lastInitialBracket = conjugateGradientQM->initialBracket;
	for (  i = 0; i <  conjugateGradientQM->numberOfAtoms; i++ )
	{

			for(j=0;j<3;j++)
				conjugateGradientQM->qmModel->molecule.atoms[i].coordinates[j] += 
				conjugateGradientQM->direction[j][i] * conjugateGradientQM->initialBracket;
	}
	return( energy );
}
/********************************************************************************/
static void bracketMinimum(ConjugateGradientQM* conjugateGradientQM, double pointA[], double pointB[], double pointC[] )
{
	static double GOLDENRATIO = 1.618034;
	double energyA, energyB, energyC, temp;
      	double ulim, u, r, q, fu, denominator;
	int iter = 0;

       	energyA = oneDimensionalEnergy(conjugateGradientQM, pointA[ 0 ] );
        energyB = oneDimensionalEnergy(conjugateGradientQM, pointB[ 0 ] );

        if ( energyB > energyA ) 
	{  
		temp = pointA[ 0 ];
		pointA[ 0 ] = pointB[ 0 ];
		pointB[ 0 ] = temp;
		temp = energyB;
		energyB = energyA;
		energyA = temp;  
        }
        pointC[ 0 ] = pointB[ 0 ] + GOLDENRATIO * ( pointB[ 0 ] - pointA[ 0 ] );
        energyC = oneDimensionalEnergy(conjugateGradientQM, pointC[ 0 ] );
        while ( energyB > energyC )
	{
		iter++;
               	r = ( pointB[ 0 ] - pointA[ 0 ] ) * ( energyB - energyC );
               	q = ( pointB[ 0 ] - pointC[ 0 ] ) * ( energyB - energyA );
		denominator = FMAX( fabs( q - r ), 1.0e-20 );
		if ( ( q - r ) < 0 )
			denominator = -denominator;
               	u = ( pointB[ 0 ] ) - ( ( pointB[ 0 ]-pointC[ 0 ] ) * q - 
			( pointB[ 0 ] - pointA[ 0 ] ) * r ) /
                       	( 2.0 * denominator );
               	ulim = pointB[ 0 ] + 100 * ( pointC[ 0 ] - pointB[ 0 ] );
               	if ( ( pointB[ 0 ] - u ) * ( u - pointC[ 0 ] ) > 0.0 )
		{
                       	fu=oneDimensionalEnergy(conjugateGradientQM, u );
                       	if ( fu < energyC )
			{
                               	pointA[ 0 ] = pointB[ 0 ];
                               	pointB[ 0 ] = u;
                               	energyA = energyB;
                               	energyB = fu;
                               	return;
                       	}
			else if ( fu > energyB )
			{
                               	pointC[ 0 ] = u;
                               	energyC = fu;
                               	return;
			}
                       	u = pointC[ 0 ] + GOLDENRATIO * ( pointC[ 0 ] - pointB[ 0 ] );
                       	fu = oneDimensionalEnergy(conjugateGradientQM, u );
               	}
		else if ( ( pointC[ 0 ] - u ) * ( u - ulim ) > 0.0 )
		{
                       	fu = oneDimensionalEnergy(conjugateGradientQM, u );
                       	if ( fu < energyC )
			{
				pointB[ 0 ] = pointC[ 0 ];
				pointC[ 0 ] = u;
				u = pointC[ 0 ] + GOLDENRATIO * ( pointC[ 0 ] - pointB[ 0 ] );
				energyB = energyC;
				energyC = fu;
				fu = oneDimensionalEnergy(conjugateGradientQM, u );
                       	}
               	}
		else if ( ( u - ulim ) * ( ulim - pointC[ 0 ] ) >= 0.0 )
		{
                       	u = ulim;
                       	fu = oneDimensionalEnergy(conjugateGradientQM, u );
               	}
		else
		{
                       	u = pointC[ 0 ] + GOLDENRATIO * ( pointC[ 0 ] - pointB[ 0 ] );
                       	fu = oneDimensionalEnergy(conjugateGradientQM, u );
               	}
		pointA[ 0 ] = pointB[ 0 ];
		pointB[ 0 ] = pointC[ 0 ];
		pointC[ 0 ] = u;
		energyA = energyB;
		energyB = energyC;
		energyC = fu;
       	}
}
/********************************************************************************/
static double oneDimensionalEnergy(ConjugateGradientQM* conjugateGradientQM, double factor )
{

	int i;
	int j;
	double energy;
	QuantumMechanicsModel copyModel = conjugateGradientQM->qmModel->klass->copy(conjugateGradientQM->qmModel);
	for (  i = 0; i < conjugateGradientQM->numberOfAtoms; i++ )
	{
		for(j=0;j<3;j++)
			copyModel.molecule.atoms[i].coordinates[j] += factor * conjugateGradientQM->direction[j][i];
	}

	conjugateGradientQM->qmModel->klass->calculateEnergy(&copyModel);
	energy = copyModel.molecule.potentialEnergy;
	copyModel.klass->free(&copyModel);
	
	return energy;
}
/********************************************************************************/
static double inverseParabolicInterpolation(ConjugateGradientQM* conjugateGradientQM,
		double pointa[], double mid[], double pointb[], 
		double minimum[] )
{
        int iter;
	int maxIterations = conjugateGradientQM->maxLine;
	static double CGOLD = 0.3819660;
        double a,b,d=0,etemp,fu,fv,fw,fx,p,q,r,tol1,tol2,u,v,w,x,xm;
        double e=0.0, tol=2.0e-4;
	double pointA, pointB;
	double energy;
	
	pointA = pointa[ 0 ];
	pointB = pointb[ 0 ];

        a=(pointA < pointB ? pointA : pointB);
        b=(pointA > pointB ? pointA : pointB);
        x=w=v=mid[ 0 ];
        fw=fv=fx=fu=energy=oneDimensionalEnergy(conjugateGradientQM,x);
        for (iter=1;iter<=maxIterations;iter++)
	{

               	xm=0.5*(a+b);
		tol1 = tol* fabs( x ) + 1.0e-10;
               	tol2=2.0*tol1;
               	if (fabs(x-xm) <= (tol2-0.5*(b-a)))
		{
                       	minimum[0]=x;
                       	return fu;
               	}
               	if (fabs(e) > tol1)
		{
                       	r=(x-w)*(fx-fv);
                       	q=(x-v)*(fx-fw);
                       	p=(x-v)*q-(x-w)*r;
                       	q=2.0*(q-r);
                       	if (q > 0.0) p = -p;
                       	q=fabs(q);
                       	etemp=e;
                       	e=d;
                       	if (fabs(p) >= fabs(0.5*q*etemp) || 
				p <= q*(a-x) || p >= q*(b-x))
                               	d=CGOLD*(e=(x >= xm ? a-x : b-x));
                       	else
			{
                               	d=p/q;
                               	u=x+d;
                               	if (u-a < tol2 || b-u < tol2)
				{
					if ( ( xm - x ) < 0 )
						d = - tol1;
					else
						d = tol1;
				}
                       	}
               	}
		else
		{
                        	d=CGOLD*(e=(x >= xm ? a-x : b-x));
               	}
		if ( fabs( d ) >= tol1 )
		{
			u = x + d;
		}
		else
		{
			if ( d >= 0 )
				u = x + tol1;
			else
				u = x - tol1;
		}
               	fu=oneDimensionalEnergy(conjugateGradientQM,u);
               	if (fu <= fx)
		{
                       	if (u >= x) a=x; else b=x;
			v = w;
			w = x;
			x = u;
			fv = fw;
			fw = fx;
			fx = fu;
               	}
		else
		{
                       	if (u < x) a=u; else b=u;
                       	if (fu <= fw || w == x)
			{
                               	v=w;
                               	w=u;
                               	fv=fw;
                               	fw=fu;
                       	}
			else if (fu <= fv || v == x || v == w)
			{
                               	v=u;
                               	fv=fu;
                       	}
               	}
        }
	return( fu );
}
/*****************************************************************************************************************************************************/
void setCGQMOptions(FILE* file, ConjugateGradientQMOptions* conjugateGradientQMOptions)
{
/* Optimsation options */ 
	conjugateGradientQMOptions->gradientNorm = 1e-3;
	conjugateGradientQMOptions->maxIterations = 100;
	conjugateGradientQMOptions->updateFrequency = 1;
	conjugateGradientQMOptions->maxLines = 25;
	conjugateGradientQMOptions->initialStep = 0.001;
/* 1 : Hestenes Stiefel,  2 : Fletcher Reeves, 3 : Polak Ribiere, 4 : Wolf Powell*/
	conjugateGradientQMOptions->method = 1;

	readOneReal(file,"conjugateGradientGradientNorm",&conjugateGradientQMOptions->gradientNorm);
	readOneInt(file,"conjugateGradientMaxIterations",&conjugateGradientQMOptions->maxIterations);
	readOneInt(file,"conjugateGradientUpdateFrequency",&conjugateGradientQMOptions->updateFrequency);
	readOneInt(file,"conjugateGradientMaxLines",&conjugateGradientQMOptions->maxLines);
	readOneReal(file,"conjugateGradientInitialStep",&conjugateGradientQMOptions->initialStep);
	readOneInt(file,"conjugateGradientMethod",&conjugateGradientQMOptions->method);
	if(conjugateGradientQMOptions->method<1||conjugateGradientQMOptions->method>4) conjugateGradientQMOptions->method=1;
}
