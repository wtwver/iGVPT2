
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

/* ConjugateGradient.c */
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>

#include "../Utils/Utils.h"
#include "../Utils/AtomsProp.h"
#include "../Utils/Constants.h"
#include "../Utils/Types.h"
#include "../Molecule/Molecule.h"
#include "../MolecularMechanics/ConjugateGradient.h"


static double maxarg1,maxarg2;
#define FMAX(a,b) (maxarg1=(a),maxarg2=(b),(maxarg1) > (maxarg2) ?\
        (maxarg1) : (maxarg2))

static void hestenesStiefel(ConjugateGradient* conjugateGradient);
static void fletcherReeves(ConjugateGradient* conjugateGradient);
static void polakRibiere(ConjugateGradient* conjugateGradient);
static void wolfPowell(ConjugateGradient* conjugateGradient);
static double lineMinimize(ConjugateGradient* conjugateGradient);
static void bracketMinimum(ConjugateGradient* conjugateGradient, double pointA[], double pointB[], double pointC[] );
static double oneDimensionalEnergy(ConjugateGradient* conjugateGradient, double factor );
static double inverseParabolicInterpolation(ConjugateGradient* conjugateGradient,
		double pointa[], double mid[], double pointb[], 
		double minimum[] );
/**********************************************************************/
void	runConjugateGradient(ConjugateGradient* conjugateGradient, ForceField* forceField, 
		ConjugateGradientOptions conjugateGradientOptions)
{

	int i;
	int minimizerOptions = conjugateGradientOptions.method;
	conjugateGradient->forceField = forceField;
	conjugateGradient->numberOfAtoms = forceField->molecule.nAtoms;
	conjugateGradient->updateFrequency = conjugateGradientOptions.updateFrequency;
	conjugateGradient->maxIterations = conjugateGradientOptions.maxIterations;
	conjugateGradient->updateNumber = 0;
	conjugateGradient->maxLine = conjugateGradientOptions.maxLines;
	conjugateGradient->epsilon = conjugateGradientOptions.gradientNorm;
	conjugateGradient->initialStep = conjugateGradientOptions.initialStep;
	conjugateGradient->rmsDeplacment = 0;
	conjugateGradient->maxDeplacment = 0;
	conjugateGradient->gradientNorm = 0;
	conjugateGradient->initialBracket = conjugateGradientOptions.initialStep;
	conjugateGradient->lastInitialBracket = 0;
	conjugateGradient->term = 0;

	for(i=0;i<3;i++)
	{
		conjugateGradient->lastGradient[i] = 
			malloc(conjugateGradient->numberOfAtoms*sizeof(double));
		conjugateGradient->direction[i] = 
			malloc(conjugateGradient->numberOfAtoms*sizeof(double));
	}

	conjugateGradient->temporaryMolecule = malloc(sizeof(Molecule));
	conjugateGradient->temporaryMolecule->nAtoms = conjugateGradient->numberOfAtoms;
	conjugateGradient->temporaryMolecule->atoms = 
		malloc(conjugateGradient->numberOfAtoms*sizeof(Atom)); 
	for (  i = 0; i < conjugateGradient->numberOfAtoms; i++ )
	{
		conjugateGradient->temporaryMolecule->atoms[i].charge =  conjugateGradient->forceField->molecule.atoms[i].charge;
		conjugateGradient->temporaryMolecule->atoms[i].mass =  conjugateGradient->forceField->molecule.atoms[i].mass;
	}

	switch(minimizerOptions)
	{
		case 1 : hestenesStiefel(conjugateGradient);break;
		case 2 : fletcherReeves(conjugateGradient);break;
		case 3 : polakRibiere(conjugateGradient);break;
		case 4 : wolfPowell(conjugateGradient);break;
	}
}
/**********************************************************************/
void	freeConjugateGradient(ConjugateGradient* conjugateGradient)
{

	int i;
	conjugateGradient->forceField = NULL;
	conjugateGradient->numberOfAtoms = 0;
	conjugateGradient->updateFrequency = 0;
	conjugateGradient->maxIterations = 0;
	conjugateGradient->updateNumber = 0;
	conjugateGradient->maxLine = 0;
	conjugateGradient->epsilon = 0;
	conjugateGradient->initialStep = 0;
	conjugateGradient->rmsDeplacment = 0;
	conjugateGradient->maxDeplacment = 0;
	conjugateGradient->gradientNorm = 0;
	conjugateGradient->initialBracket = 0;
	conjugateGradient->lastInitialBracket = 0;
	conjugateGradient->term = 0;

	for(i=0;i<3;i++)
	{
		if(conjugateGradient->lastGradient[i] !=NULL)
		{
			free(conjugateGradient->lastGradient[i]);
			conjugateGradient->lastGradient[i]= NULL;
		}
		if(conjugateGradient->direction[i] != NULL)
		{
			free(conjugateGradient->direction[i]);
			conjugateGradient->direction[i] = NULL; 
		}
	}

	if(conjugateGradient->temporaryMolecule != NULL)
	{
		free(conjugateGradient->temporaryMolecule);
		conjugateGradient->temporaryMolecule = NULL;
	}

}
/**********************************************************************/
static void fletcherReeves(ConjugateGradient* conjugateGradient)
{
	double lastGradientDotGradient = 0, gradientDotGradient = 0, beta;
	int iterations = 0;
	int i;
	int j;
	double energy;
	char* str = strdup(" ");

	conjugateGradient->forceField->klass->calculateGradient(conjugateGradient->forceField);
	for ( i = 0; i < conjugateGradient->numberOfAtoms; i++ )
	{
		for(j=0;j<3;j++)
		{
		conjugateGradient->direction[j][ i ] = 
			-conjugateGradient->forceField->molecule.atoms[i].gradient[j];
		lastGradientDotGradient += conjugateGradient->forceField->molecule.atoms[i].gradient[j]
					* conjugateGradient->forceField->molecule.atoms[i].gradient[j]; 
		}
	}
	gradientDotGradient = lastGradientDotGradient;
	while ( 
			( lastGradientDotGradient > conjugateGradient->epsilon ) && 
			( iterations++ < conjugateGradient->maxIterations ) 
		)
	{

		lineMinimize(conjugateGradient);
		conjugateGradient->forceField->klass->calculateGradient(conjugateGradient->forceField);
		lastGradientDotGradient = gradientDotGradient;	
		gradientDotGradient = 0;	
		for (  i = 0; i < conjugateGradient->numberOfAtoms; i++ )
			for(j=0;j<3;j++)
			{
				gradientDotGradient += 
					conjugateGradient->forceField->molecule.atoms[i].gradient[j]
					* conjugateGradient->forceField->molecule.atoms[i].gradient[j]; 
			}
		beta = gradientDotGradient / lastGradientDotGradient;
		for (  i = 0; i < conjugateGradient->numberOfAtoms; i++ )
		{
			for(j=0;j<3;j++)
			{
				conjugateGradient->direction[j][i] = 
					beta *  conjugateGradient->direction[j][i] - 
					conjugateGradient->forceField->molecule.atoms[i].gradient[j];
			}
		}
		if ( conjugateGradient->updateNumber >= conjugateGradient->updateFrequency )
		{
			free(str);
			str = strdup_printf(("gradient = %f "),sqrt(gradientDotGradient)); 
			/* redrawMolecule(&conjugateGradient->forceField->molecule,str);*/
			fprintf(conjugateGradient->logfile, "%s\n",str);
			fflush(conjugateGradient->logfile);
			conjugateGradient->updateNumber = 0;
		}
		conjugateGradient->updateNumber++;
	}	

	conjugateGradient->forceField->klass->calculateGradient(conjugateGradient->forceField);
	conjugateGradient->gradientNorm = 0;
	for(  i = 0; i < conjugateGradient->numberOfAtoms; i++ )
	{
		for(j=0;j<3;j++)
		{
			conjugateGradient->gradientNorm += 
				conjugateGradient->forceField->molecule.atoms[i].gradient[j] *
				conjugateGradient->forceField->molecule.atoms[i].gradient[j]; 
		}
	}
	conjugateGradient->gradientNorm = sqrt(conjugateGradient->gradientNorm );
	updateGeometryCL(conjugateGradient->forceField,  &conjugateGradient->forceField->molecule);
	energy = conjugateGradient->forceField->klass->calculateEnergyTmp
		(conjugateGradient->forceField, &conjugateGradient->forceField->molecule );

	free(str);
	str = strdup_printf(("Gradient = %f  Energy = %f (Kcal/mol) "),
			(double)conjugateGradient->gradientNorm,(double)energy); 

	/* redrawMolecule(&conjugateGradient->forceField->molecule,str);*/
	fprintf(conjugateGradient->logfile, "%s\n",str);
	fflush(conjugateGradient->logfile);
	free(str);
}
/********************************************************************************/
static void polakRibiere(ConjugateGradient* conjugateGradient)
{
	double lastGradientDotGradient = 0, gradientDotGradient, beta;
	int iterations = 0;
	double energy;
	int i;
	int j;
	char* str = strdup(" ");

	conjugateGradient->forceField->klass->calculateGradient(conjugateGradient->forceField);

	for ( i = 0; i < conjugateGradient->numberOfAtoms; i++ )
	{
		for(j=0;j<3;j++)
		{
		conjugateGradient->direction[j][ i ] = 
		conjugateGradient->lastGradient[j][i] =
			-conjugateGradient->forceField->molecule.atoms[i].gradient[j];

		lastGradientDotGradient += conjugateGradient->forceField->molecule.atoms[i].gradient[j] 
					* conjugateGradient->forceField->molecule.atoms[i].gradient[j]; 
		}
	}
	gradientDotGradient = lastGradientDotGradient;
	while ( 
			( lastGradientDotGradient > conjugateGradient->epsilon ) && 
			( iterations++ < conjugateGradient->maxIterations ) 
		)
	{

		lineMinimize(conjugateGradient);
		conjugateGradient->forceField->klass->calculateGradient(conjugateGradient->forceField);
		lastGradientDotGradient = gradientDotGradient;
		gradientDotGradient = 0;	
		beta = 0;
		for (  i = 0; i < conjugateGradient->numberOfAtoms; i++ )
		{
			for (  j = 0; j < 3; j++ )
			{
				gradientDotGradient += 
					conjugateGradient->forceField->molecule.atoms[i].gradient[j] *
					conjugateGradient->forceField->molecule.atoms[i].gradient[j]; 
				
				beta += ( 
					conjugateGradient->forceField->molecule.atoms[i].gradient[j] - 
					conjugateGradient->lastGradient[j][i]
					)* conjugateGradient->forceField->molecule.atoms[i].gradient[j]; 
			}
		}
		beta /= lastGradientDotGradient;
		for (  i = 0; i < conjugateGradient->numberOfAtoms; i++ )
		{
			for (  j = 0; j < 3; j++ )
			{
				conjugateGradient->direction[j][i] = 
				beta * conjugateGradient->direction[j][i] - 
				conjugateGradient->forceField->molecule.atoms[i].gradient[j];
			}
		}
		for ( i = 0; i < conjugateGradient->numberOfAtoms; i++ )
		{
			for(j=0;j<3;j++)
			{
				conjugateGradient->lastGradient[j][i] = 
				conjugateGradient->forceField->molecule.atoms[i].gradient[j];
			}
		}
		if ( conjugateGradient->updateNumber >= conjugateGradient->updateFrequency )
		{
			free(str);
			str = strdup_printf("gradient = %f ",sqrt(gradientDotGradient)); 
			/* redrawMolecule(&conjugateGradient->forceField->molecule,str);*/
			fprintf(conjugateGradient->logfile, "%s\n",str);
			fflush(conjugateGradient->logfile);
			conjugateGradient->updateNumber = 0;
		}
		conjugateGradient->updateNumber++;
	}	
	conjugateGradient->forceField->klass->calculateGradient(conjugateGradient->forceField);
	conjugateGradient->gradientNorm = 0;
	for(  i = 0; i < conjugateGradient->numberOfAtoms; i++ )
	{
		for (  j = 0; j < 3; j++ )
		{
			conjugateGradient->gradientNorm += 
				conjugateGradient->forceField->molecule.atoms[i].gradient[j] *
				conjugateGradient->forceField->molecule.atoms[i].gradient[j]; 
		}
	}
	conjugateGradient->gradientNorm = sqrt(conjugateGradient->gradientNorm );
	updateGeometryCL(conjugateGradient->forceField,  &conjugateGradient->forceField->molecule);
	energy = conjugateGradient->forceField->klass->calculateEnergyTmp
		(conjugateGradient->forceField, &conjugateGradient->forceField->molecule );

	free(str);
	str = strdup_printf(("Gradient = %f  Energy = %f (Kcal/mol) "),
			(double)conjugateGradient->gradientNorm,(double)energy); 
	/* redrawMolecule(&conjugateGradient->forceField->molecule,str);*/
	fprintf(conjugateGradient->logfile, "%s\n",str);
	fflush(conjugateGradient->logfile);
	free(str);
}
/********************************************************************************/
static void hestenesStiefel(ConjugateGradient* conjugateGradient)
{
	double gradientDotGradient = 1, beta, gradientDiff, denom;
	int iterations = 0;
	double energy;
	int i;
	int j;
	char* str = strdup(" ");

	conjugateGradient->forceField->klass->calculateGradient(conjugateGradient->forceField);
	for ( i = 0; i < conjugateGradient->numberOfAtoms; i++ )
	{
		for(j=0;j<3;j++)
		{
		conjugateGradient->direction[j][ i ] = 
		conjugateGradient->lastGradient[j][i] = 
			-conjugateGradient->forceField->molecule.atoms[i].gradient[j];
		}
	}
	while ( 
			( gradientDotGradient > conjugateGradient->epsilon ) && 
			( iterations++ < conjugateGradient->maxIterations ) 
		)
	{

		lineMinimize(conjugateGradient);
		conjugateGradient->forceField->klass->calculateGradient(conjugateGradient->forceField);
		gradientDotGradient = 0;	
		for (  i = 0; i < conjugateGradient->numberOfAtoms; i++ )
		{
			for (  j = 0; j < 3; j++ )
			{
				gradientDotGradient += 
					conjugateGradient->forceField->molecule.atoms[i].gradient[j]*
					conjugateGradient->forceField->molecule.atoms[i].gradient[j];
			}
		}
		beta = 0;
		denom = 0;
		for (  i = 0; i < conjugateGradient->numberOfAtoms; i++ )
		{
			for (  j = 0; j < 3; j++ )
			{
				gradientDiff = 
					conjugateGradient->forceField->molecule.atoms[i].gradient[j]-
					conjugateGradient->lastGradient[j][i];

				beta += 
					conjugateGradient->forceField->molecule.atoms[i].gradient[j]*
				       	gradientDiff;
				denom += conjugateGradient->direction[j][i] * gradientDiff;
			}
		}
		if ( fabs( denom ) > 1.0e-10 )
			beta /= denom;
		else
			beta = 0;

		for (  i = 0; i < conjugateGradient->numberOfAtoms; i++ )
		{
			for (  j = 0; j < 3; j++ )
			{
				conjugateGradient->direction[j][ i ] = 
					beta * conjugateGradient->direction[j][ i ] - 
					conjugateGradient->forceField->molecule.atoms[i].gradient[j];
			}
		}
		for ( i = 0; i < conjugateGradient->numberOfAtoms; i++ )
		{
			for(j=0;j<3;j++)
			{
				conjugateGradient->lastGradient[j][i] = 
				conjugateGradient->forceField->molecule.atoms[i].gradient[j];
			}
		}
		if ( conjugateGradient->updateNumber >= conjugateGradient->updateFrequency )
		{
			free(str);
			str = strdup_printf(("gradient = %f "),sqrt(gradientDotGradient)); 
			/* redrawMolecule(&conjugateGradient->forceField->molecule,str);*/
			fprintf(conjugateGradient->logfile, "%s\n",str);
			fflush(conjugateGradient->logfile);
			conjugateGradient->updateNumber = 0;
		}
		conjugateGradient->updateNumber++;
	}	
	conjugateGradient->forceField->klass->calculateGradient(conjugateGradient->forceField);
	conjugateGradient->gradientNorm = 0;
	for(  i = 0; i < conjugateGradient->numberOfAtoms; i++ )
	{
		for(j=0;j<3;j++)
			conjugateGradient->gradientNorm += 
			conjugateGradient->direction[j][ i ] *
			conjugateGradient->direction[j][ i ];
	}
	conjugateGradient->gradientNorm = sqrt( conjugateGradient->gradientNorm );
	updateGeometryCL(conjugateGradient->forceField,  &conjugateGradient->forceField->molecule);
	energy = conjugateGradient->forceField->klass->calculateEnergyTmp
		(conjugateGradient->forceField, &conjugateGradient->forceField->molecule );

	free(str);
	str = strdup_printf(("Gradient = %f  Energy = %f (Kcal/mol) "),
			(double)conjugateGradient->gradientNorm,(double)energy); 
	/* redrawMolecule(&conjugateGradient->forceField->molecule,str);*/
	fprintf(conjugateGradient->logfile, "%s\n",str);
	fflush(conjugateGradient->logfile);
	free(str);
}
/********************************************************************************/
static void wolfPowell(ConjugateGradient* conjugateGradient)
{
	double lastGradientDotGradient = 0, gradientDotGradient, beta;
	int iterations = 0;
	double energy;
	int i;
	int j;
	char* str = strdup(" ");

	conjugateGradient->forceField->klass->calculateGradient(conjugateGradient->forceField);
	for ( i = 0; i < conjugateGradient->numberOfAtoms; i++ )
	{
		for(j=0;j<3;j++)
		{
		conjugateGradient->direction[j][ i ] = 
		conjugateGradient->lastGradient[j][i] =
			-conjugateGradient->forceField->molecule.atoms[i].gradient[j];

		lastGradientDotGradient += conjugateGradient->forceField->molecule.atoms[i].gradient[j]
					* conjugateGradient->forceField->molecule.atoms[i].gradient[j]; 
		}
	}
	while ( 
			( lastGradientDotGradient > conjugateGradient->epsilon ) && 
			( iterations++ < conjugateGradient->maxIterations ) 
		)
	{

		lineMinimize(conjugateGradient);
		conjugateGradient->forceField->klass->calculateGradient(conjugateGradient->forceField);
		gradientDotGradient = 0;	
		for ( i = 0; i < conjugateGradient->numberOfAtoms; i++ )
			for ( j = 0; j < 3; j++ )
			{
				gradientDotGradient += 
					conjugateGradient->forceField->molecule.atoms[i].gradient[j]*
					conjugateGradient->forceField->molecule.atoms[i].gradient[j];
			}

		beta = 0;
		for (  i = 0; i < conjugateGradient->numberOfAtoms; i++ )
		{
			for ( j = 0; j < 3; j++ )
				beta += ( 
					conjugateGradient->forceField->molecule.atoms[i].gradient[j] - 
					conjugateGradient->lastGradient[j][i]
					)* conjugateGradient->forceField->molecule.atoms[i].gradient[j]; 
		}
		beta /= lastGradientDotGradient;
		if ( beta < 0 )
			beta = 0;
		lastGradientDotGradient = gradientDotGradient;
		for (  i = 0; i < conjugateGradient->numberOfAtoms; i++ )
			for ( j = 0; j < 3; j++ )
				conjugateGradient->direction[j][ i ] = 
					beta * conjugateGradient->direction[j][ i ] - 
					conjugateGradient->forceField->molecule.atoms[i].gradient[j];
		
		for ( i = 0; i < conjugateGradient->numberOfAtoms; i++ )
		{
			for(j=0;j<3;j++)
			{
				conjugateGradient->lastGradient[j][i] = 
				conjugateGradient->forceField->molecule.atoms[i].gradient[j];
			}
		}
		if ( conjugateGradient->updateNumber >= conjugateGradient->updateFrequency )
		{
			free(str);
			str = strdup_printf(("gradient = %f "),sqrt(gradientDotGradient)); 
			/* redrawMolecule(&conjugateGradient->forceField->molecule,str);*/
			fprintf(conjugateGradient->logfile, "%s\n",str);
			fflush(conjugateGradient->logfile);
			conjugateGradient->updateNumber = 0;
		}
		conjugateGradient->updateNumber++;
	}	
	conjugateGradient->forceField->klass->calculateGradient(conjugateGradient->forceField);
	conjugateGradient->gradientNorm = 0;
	for(  i = 0; i < conjugateGradient->numberOfAtoms; i++ )
	{
		for(j=0;j<3;j++)
			conjugateGradient->gradientNorm += 
			conjugateGradient->direction[j][ i ] *
			conjugateGradient->direction[j][ i ];
	}
	conjugateGradient->gradientNorm = sqrt( conjugateGradient->gradientNorm );
	updateGeometryCL(conjugateGradient->forceField,  &conjugateGradient->forceField->molecule);
	energy = conjugateGradient->forceField->klass->calculateEnergyTmp
		(conjugateGradient->forceField, &conjugateGradient->forceField->molecule );

	free(str);
	str = strdup_printf(("Gradient = %f  Energy = %f (Kcal/mol) "),
			(double)conjugateGradient->gradientNorm,(double)energy); 
	/* redrawMolecule(&conjugateGradient->forceField->molecule,str);*/
	fprintf(conjugateGradient->logfile, "%s\n",str);
	fflush(conjugateGradient->logfile);
	free(str);
} 
/**********************************************************************/
static double lineMinimize(ConjugateGradient* conjugateGradient)
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
	b = conjugateGradient->initialBracket;
	bracketMinimum(conjugateGradient, &a, &b, &c );
	energy = inverseParabolicInterpolation(conjugateGradient, &a, &b, &c, &minimum );
	conjugateGradient->initialBracket = minimum;

	if ( 	( fabs( conjugateGradient->initialBracket ) < delta ) || 
		( conjugateGradient->initialBracket == conjugateGradient->lastInitialBracket ) )
	{
		conjugateGradient->initialBracket = 
			rand()/(double)RAND_MAX *conjugateGradient->initialStep;

		for (  i = 0; i < conjugateGradient->numberOfAtoms; i++ )
		{
			for(j=0;j<3;j++)
			conjugateGradient->direction[j][i] = 
				-conjugateGradient->forceField->molecule.atoms[i].gradient[j];
		}
	} 
	conjugateGradient->lastInitialBracket = conjugateGradient->initialBracket;
	for (  i = 0; i <  conjugateGradient->numberOfAtoms; i++ )
	{

			for(j=0;j<3;j++)
				conjugateGradient->forceField->molecule.atoms[i].coordinates[j] += 
				conjugateGradient->direction[j][i] * conjugateGradient->initialBracket;
	}
	updateGeometryCL( conjugateGradient->forceField,NULL);
	return( energy );
}
/********************************************************************************/
static void bracketMinimum(ConjugateGradient* conjugateGradient, double pointA[], double pointB[], double pointC[] )
{
	static double GOLDENRATIO = 1.618034;
	double energyA, energyB, energyC, temp;
      	double ulim, u, r, q, fu, denominator;
	int iter = 0;

       	energyA = oneDimensionalEnergy(conjugateGradient, pointA[ 0 ] );
        energyB = oneDimensionalEnergy(conjugateGradient, pointB[ 0 ] );

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
        energyC = oneDimensionalEnergy(conjugateGradient, pointC[ 0 ] );
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
                       	fu=oneDimensionalEnergy(conjugateGradient, u );
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
                       	fu = oneDimensionalEnergy(conjugateGradient, u );
               	}
		else if ( ( pointC[ 0 ] - u ) * ( u - ulim ) > 0.0 )
		{
                       	fu = oneDimensionalEnergy(conjugateGradient, u );
                       	if ( fu < energyC )
			{
				pointB[ 0 ] = pointC[ 0 ];
				pointC[ 0 ] = u;
				u = pointC[ 0 ] + GOLDENRATIO * ( pointC[ 0 ] - pointB[ 0 ] );
				energyB = energyC;
				energyC = fu;
				fu = oneDimensionalEnergy(conjugateGradient, u );
                       	}
               	}
		else if ( ( u - ulim ) * ( ulim - pointC[ 0 ] ) >= 0.0 )
		{
                       	u = ulim;
                       	fu = oneDimensionalEnergy(conjugateGradient, u );
               	}
		else
		{
                       	u = pointC[ 0 ] + GOLDENRATIO * ( pointC[ 0 ] - pointB[ 0 ] );
                       	fu = oneDimensionalEnergy(conjugateGradient, u );
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
static double oneDimensionalEnergy(ConjugateGradient* conjugateGradient, double factor )
{

	int i;
	int j;
	for (  i = 0; i < conjugateGradient->numberOfAtoms; i++ )
	{
		for(j=0;j<3;j++)
		{
			conjugateGradient->temporaryMolecule->atoms[i].coordinates[j] = 
				conjugateGradient->forceField->molecule.atoms[i].coordinates[j] + 
				factor * conjugateGradient->direction[j][i];
		}
	}
	
	updateGeometryCL(conjugateGradient->forceField,  conjugateGradient->temporaryMolecule);
	return ( conjugateGradient->forceField->klass->calculateEnergyTmp(conjugateGradient->forceField,conjugateGradient->temporaryMolecule) );
}
/********************************************************************************/
static double inverseParabolicInterpolation(ConjugateGradient* conjugateGradient,
		double pointa[], double mid[], double pointb[], 
		double minimum[] )
{
        int iter;
	int maxIterations = conjugateGradient->maxLine;
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
        fw=fv=fx=fu=energy=oneDimensionalEnergy(conjugateGradient,x);
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
               	fu=oneDimensionalEnergy(conjugateGradient,u);
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
void setCGOptions(FILE* file, ConjugateGradientOptions* conjugateGradientOptions)
{
/* Optimsation options */ 
	conjugateGradientOptions->gradientNorm = 1e-3;
	conjugateGradientOptions->maxIterations = 100;
	conjugateGradientOptions->updateFrequency = 1;
	conjugateGradientOptions->maxLines = 25;
	conjugateGradientOptions->initialStep = 0.001;
/* 1 : Hestenes Stiefel,  2 : Fletcher Reeves, 3 : Polak Ribiere, 4 : Wolf Powell*/
	conjugateGradientOptions->method = 1;

	readOneReal(file,"conjugateGradientGradientNorm",&conjugateGradientOptions->gradientNorm);
	readOneInt(file,"conjugateGradientMaxIterations",&conjugateGradientOptions->maxIterations);
	readOneInt(file,"conjugateGradientUpdateFrequency",&conjugateGradientOptions->updateFrequency);
	readOneInt(file,"conjugateGradientMaxLines",&conjugateGradientOptions->maxLines);
	readOneReal(file,"conjugateGradientInitialStep",&conjugateGradientOptions->initialStep);
	readOneInt(file,"conjugateGradientMethod",&conjugateGradientOptions->method);
	if(conjugateGradientOptions->method<1||conjugateGradientOptions->method>4) conjugateGradientOptions->method=1;
}
