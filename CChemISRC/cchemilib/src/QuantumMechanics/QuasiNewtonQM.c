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

/* QuasiNewtonQM.c */
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>

#include "../Utils/Utils.h"
#include "../Utils/AtomsProp.h"
#include "../Utils/Constants.h"
#include "../Utils/Types.h"
#include "../Molecule/Molecule.h"
#include "../QuantumMechanics/QuasiNewtonQM.h"

static double maxarg1,maxarg2;
#define FMIN(a,b) (maxarg1=(a),maxarg2=(b),(maxarg1) > (maxarg2) ?(maxarg2) : (maxarg1))
#define FMAX(a,b) (maxarg1=(a),maxarg2=(b),(maxarg1) > (maxarg2) ?(maxarg1) : (maxarg2))


static int lbfgs( 
		int n , int m , double x[] , double f , double g[] ,
		int diagco , double diag[] , double eps ,
		double xtol , int maxLines,int iflag[]
	);
/**********************************************************************/
void runQuasiNewtonQM(QuasiNewtonQM* quasiNewtonQM)
{

	int j;
	int i,i3;
	int iter;
	int nAtomsX3;
	int diagco = FALSE;
	double* x = NULL;
	double* g = NULL;
	double* diag = NULL;
	int iflag = 0;
	double energy = 0;
	int nAtoms;
	int updateNumber = 0;
	char str[BSIZE];
	double gradientNorm;
	FILE* fileOut = quasiNewtonQM->logfile;

	QuantumMechanicsModel* qmModel = quasiNewtonQM->qmModel;

	if(qmModel->molecule.nAtoms<1) return;
	fprintf(fileOut,"=============================================================================================================================================\n");
        fprintf(fileOut,"\t\t\tOptimization by Quasi-Newton method\n");
        fprintf(fileOut,"---------------------------------------------------------------------------------------------------------------------------------------------\n");
        fprintf(fileOut,"\t\t\tMaxIteration \t\t= %d\n",quasiNewtonQM->maxIterations);
        fprintf(fileOut,"\t\t\tEpsilon \t\t= %0.4e\n",quasiNewtonQM->epsilon);
        fprintf(fileOut,"\t\t\tMax lines search \t= %d\n",quasiNewtonQM->maxLines);
        fprintf(fileOut,"=============================================================================================================================================\n");
        fflush(fileOut); fflush(stderr);


	nAtoms = qmModel->molecule.nAtoms;
	nAtomsX3 = 3*nAtoms;

	diag = malloc(nAtomsX3*sizeof(double));
	x = malloc(nAtomsX3*sizeof(double));
	g = malloc(nAtomsX3*sizeof(double));

	qmModel->klass->calculateEnergy(qmModel);
	qmModel->klass->calculateGradient(qmModel);
	energy = qmModel->molecule.potentialEnergy;

	gradientNorm = 0;
	for (  i = 0; i < nAtoms; i++ )
	for(j=0;j<3;j++) gradientNorm += qmModel->molecule.atoms[i].gradient[j]*qmModel->molecule.atoms[i].gradient[j]; 
	sprintf(str,("Gradient(kcal/mol/Ang) = %0.14f\tEnergy = %0.14f "),sqrt(gradientNorm),energy); 
	fprintf(quasiNewtonQM->logfile,"%s\n",str);
	fflush(quasiNewtonQM->logfile);

	for(iter=0;iter<quasiNewtonQM->maxIterations;iter++)
	{

		//fprintf(stdout,"Begin calculateGradient\n");
		qmModel->klass->calculateGradient(qmModel);
		//fprintf(stdout,"End calculateGradient\n");
		/* energy = qmModel->klass->calculateEnergyTmp(qmModel, &qmModel->molecule );*/
		//qmModel->klass->calculateEnergy(qmModel);
		energy = qmModel->molecule.potentialEnergy;
		//printf("energy = %f energyC = %f\n",energy,qmModel->molecule.potentialEnergy);
		/* set x  and g table from coordinates and gradient */
		for(i=0,i3=0;i<nAtoms;i++)
		{
			if(!qmModel->molecule.atoms[i].variable) continue;
			x[i3  ] = qmModel->molecule.atoms[i].coordinates[0];
			x[i3+1] = qmModel->molecule.atoms[i].coordinates[1];
			x[i3+2] = qmModel->molecule.atoms[i].coordinates[2];

			g[i3  ] = qmModel->molecule.atoms[i].gradient[0];
			g[i3+1] = qmModel->molecule.atoms[i].gradient[1];
			g[i3+2] = qmModel->molecule.atoms[i].gradient[2];
			i3 += 3;
		}
		lbfgs(i3, i3,x, energy,g,diagco,diag,
				quasiNewtonQM->epsilon,quasiNewtonQM->tolerence,
				quasiNewtonQM->maxLines,
				&iflag);
		/*
		lbfgs(nAtomsX3, nAtomsX3,x, energy,g,diagco,diag,
				quasiNewtonQM->epsilon,quasiNewtonQM->tolerence,
				quasiNewtonQM->maxLines,
				&iflag);
				*/
		/* set coordinates from x */
		for(i=0,i3=0;i<nAtoms;i++)
		{
			if(qmModel->molecule.atoms[i].variable) 
			{
				double dmax=quasiNewtonQM->maxStep;
				double dd;
				int c;
				for(c=0;c<3;c++)
				{
					double sign=1.0;
					int ic=i3+c;
					dd = x[ic]-qmModel->molecule.atoms[i].coordinates[c];
					if(dd<0) sign=-1.0;
					if( dmax>0 && fabs(dd)>dmax) x[ic] = qmModel->molecule.atoms[i].coordinates[c]+dmax*sign;
					else qmModel->molecule.atoms[i].coordinates[c] = x[ic];
				}
				i3+=3;
			}
		}

		if ( updateNumber >= quasiNewtonQM->updateFrequency )
		{
			gradientNorm = 0;
			for (  i = 0; i < nAtoms; i++ ) for(j=0;j<3;j++) gradientNorm += qmModel->molecule.atoms[i].gradient[j]*qmModel->molecule.atoms[i].gradient[j]; 


			sprintf(str,("Iter # %d/%d\t Gradient(kcal/mol/Ang) = %0.14f\tEnergy(kcal/mol) = %0.14f "),
			iter,
			quasiNewtonQM->maxIterations,
			sqrt(gradientNorm),energy); 
			fprintf(quasiNewtonQM->logfile,"%s\n",str);
			fflush(quasiNewtonQM->logfile);
			/* redrawMolecule(&qmModel->molecule,str);*/
			updateNumber = 0;
		}
		updateNumber++;
		if(iflag<=0)
			break;
	}
	gradientNorm = 0;
	for (  i = 0; i < nAtoms; i++ )
		for(j=0;j<3;j++)
			gradientNorm += 
			qmModel->molecule.atoms[i].gradient[j]
			*qmModel->molecule.atoms[i].gradient[j]; 

	sprintf(str,("End Optimization\nGradient(kcal/mol/Ang) = %0.14f\tEnergy(kcal/mol) = %0.14f "),sqrt(gradientNorm),energy); 
	/* redrawMolecule(&qmModel->molecule,str);*/
	fprintf(quasiNewtonQM->logfile,"%s\n\n",str);
	fflush(quasiNewtonQM->logfile);
	free(diag);
	free(x);
	free(g);
}
/**********************************************************************/
void	freeQuasiNewtonQM(QuasiNewtonQM* quasiNewtonQM)
{
	quasiNewtonQM->qmModel = NULL;
	quasiNewtonQM->updateFrequency = 0;
	quasiNewtonQM->maxIterations = 0;
	quasiNewtonQM->maxLines = 0;
	quasiNewtonQM->epsilon = 0;
	quasiNewtonQM->tolerence = 0;
	quasiNewtonQM->maxStep = -1;
}
/**********************************************************************/
static double sqr( double x )
{ 
	return x*x;
}
/******************************************************************************/
static double max3( double x, double y, double z )
{
	return x < y ? ( y < z ? z : y ) : ( x < z ? z : x );
}
/** The purpose of this function is to compute a safeguarded step for
  * a linesearch and to update an interval of uncertainty for
  * a minimizer of the function.<p>
*/ 
static void mcstep (   double stx[] , double fx[] , double dx[] ,
		double sty[] , double fy[] , double dy[] ,
		double stp[] , double fp , double dp ,
		int brackt[] , double stpmin , double stpmax , int info[]
   	    )
{
	int bound;
	double gamma, p, q, r, s, sgnd, stpc, stpf, stpq, theta;

	info[0] = 0;

	if ( (
		brackt[0] && 
		( stp[0] <= FMIN ( stx[0] , sty[0] ) || stp[0] >= FMAX ( stx[0] , sty[0] ))
	     ) 
	     || dx[0] * ( stp[0] - stx[0] ) >= 0.0 
	     || stpmax < stpmin
	    )
		return;

	/* Determine if the derivatives have opposite sign.*/

	sgnd = dp * ( dx[0] / fabs ( dx[0] ) );

	if ( fp > fx[0] )
	{
		/* 
		 First case. A higher function value.
		 The minimum is bracketed. If the cubic step is closer
		 to stx than the quadratic step, the cubic step is taken,
		 else the average of the cubic and quadratic steps is taken.
		*/

		info[0] = 1;
		bound = TRUE;
		theta = 3 * ( fx[0] - fp ) / ( stp[0] - stx[0] ) + dx[0] + dp;
		s = max3 ( fabs ( theta ) , fabs ( dx[0] ) , fabs ( dp ) );
		gamma = s * sqrt ( sqr( theta / s ) - ( dx[0] / s ) * ( dp / s ) );
		if ( stp[0] < stx[0] ) gamma = - gamma;
		p = ( gamma - dx[0] ) + theta;
		q = ( ( gamma - dx[0] ) + gamma ) + dp;
		r = p/q;
		stpc = stx[0] + r * ( stp[0] - stx[0] );
		stpq = stx[0] 
			+ ( ( dx[0] / ( ( fx[0] - fp ) / ( stp[0] - stx[0] ) + dx[0] ) ) / 2 )
			*( stp[0] - stx[0] );

		if ( fabs ( stpc - stx[0] ) < fabs ( stpq - stx[0] ) )
			stpf = stpc;
		else
			stpf = stpc + ( stpq - stpc ) / 2;

		brackt[0] = TRUE;
	}
	else if ( sgnd < 0.0 )
	{
		/* Second case. A lower function value and derivatives of
		   opposite sign. The minimum is bracketed. If the cubic
		   step is closer to stx than the quadratic (secant) step,
		   the cubic step is taken, else the quadratic step is taken.
		 */

		info[0] = 2;
		bound = FALSE;
		theta = 3 * ( fx[0] - fp ) / ( stp[0] - stx[0] ) + dx[0] + dp;
		s = max3 ( fabs ( theta ) , fabs ( dx[0] ) , fabs ( dp ) );
		gamma = s * sqrt ( sqr( theta / s ) - ( dx[0] / s ) * ( dp / s ) );
		if ( stp[0] > stx[0] ) gamma = - gamma;
		p = ( gamma - dp ) + theta;
		q = ( ( gamma - dp ) + gamma ) + dx[0];
		r = p/q;
		stpc = stp[0] + r * ( stx[0] - stp[0] );
		stpq = stp[0] + ( dp / ( dp - dx[0] ) ) * ( stx[0] - stp[0] );
		if ( fabs ( stpc - stp[0] ) > fabs ( stpq - stp[0] ) )
			stpf = stpc;
		else
			stpf = stpq;

		brackt[0] = TRUE;
	}
	else if ( fabs ( dp ) < fabs ( dx[0] ) )
	{
		/* Third case. A lower function value, derivatives of the
		   same sign, and the magnitude of the derivative decreases.
		   The cubic step is only used if the cubic tends to infinity
		   in the direction of the step or if the minimum of the cubic
		   is beyond stp. Otherwise the cubic step is defined to be
		   either stpmin or stpmax. The quadratic (secant) step is also
		   computed and if the minimum is bracketed then the the step
		 closest to stx is taken, else the step farthest away is taken.
		 */

		info[0] = 3;
		bound = TRUE;
		theta = 3 * ( fx[0] - fp ) / ( stp[0] - stx[0] ) + dx[0] + dp;
		s = max3 ( fabs ( theta ) , fabs ( dx[0] ) , fabs ( dp ) );
		gamma = s * sqrt ( FMAX ( 0, sqr( theta / s ) - ( dx[0] / s ) * ( dp / s ) ) );
		if ( stp[0] > stx[0] )
		       	gamma = - gamma;
		p = ( gamma - dp ) + theta;
		q = ( gamma + ( dx[0] - dp ) ) + gamma;
		r = p/q;
		if ( r < 0.0 && gamma != 0.0 )
			stpc = stp[0] + r * ( stx[0] - stp[0] );
		else if ( stp[0] > stx[0] )
			stpc = stpmax;
		else
			stpc = stpmin;

		stpq = stp[0] + ( dp / ( dp - dx[0] ) ) * ( stx[0] - stp[0] );
		if ( brackt[0] )
		{
			if ( fabs ( stp[0] - stpc ) < fabs ( stp[0] - stpq ) )
				stpf = stpc;
			else
				stpf = stpq;
		}
		else
		{
			if ( fabs ( stp[0] - stpc ) > fabs ( stp[0] - stpq ) )
				stpf = stpc;
			else
				stpf = stpq;
		}
	}
	else
	{
		/* Fourth case. A lower function value, derivatives of the
		   same sign, and the magnitude of the derivative does
		   not decrease. If the minimum is not bracketed, the step
		   is either stpmin or stpmax, else the cubic step is taken.
		*/

		info[0] = 4;
		bound = FALSE;
		if ( brackt[0] )
		{
			theta = 3 * ( fp - fy[0] ) / ( sty[0] - stp[0] ) + dy[0] + dp;
			s = max3 ( fabs ( theta ) , fabs ( dy[0] ) , fabs ( dp ) );
			gamma = s * sqrt ( sqr( theta / s ) - ( dy[0] / s ) * ( dp / s ) );
			if ( stp[0] > sty[0] ) gamma = - gamma;
			p = ( gamma - dp ) + theta;
			q = ( ( gamma - dp ) + gamma ) + dy[0];
			r = p/q;
			stpc = stp[0] + r * ( sty[0] - stp[0] );
			stpf = stpc;
		}
		else if ( stp[0] > stx[0] )
		{
			stpf = stpmax;
		}
		else
		{
			stpf = stpmin;
		}
	}

	/* Update the interval of uncertainty. This update does not
	   depend on the new step or the case analysis above.
	*/

	if ( fp > fx[0] )
	{
		sty[0] = stp[0];
		fy[0] = fp;
		dy[0] = dp;
	}
	else
	{
		if ( sgnd < 0.0 )
		{
			sty[0] = stx[0];
			fy[0] = fx[0];
			dy[0] = dx[0];
		}
		stx[0] = stp[0];
		fx[0] = fp;
		dx[0] = dp;
	}

	/* Compute the new step and safeguard it.*/

	stpf = FMIN ( stpmax , stpf );
	stpf = FMAX ( stpmin , stpf );
	stp[0] = stpf;

	if ( brackt[0] && bound )
	{
		if ( sty[0] > stx[0] )
		{
			stp[0] = FMIN ( stx[0] + 0.66 * ( sty[0] - stx[0] ) , stp[0] );
		}
		else
		{
			stp[0] = FMAX ( stx[0] + 0.66 * ( sty[0] - stx[0] ) , stp[0] );
		}
	}

	return;
}
/******************************************************************************/
/* Minimize a function along a search direction. */
static void mcsrch ( int n , double x[] , double f , double g[] ,
	      double s[] , int is0 , double stp[] , double ftol , double xtol ,
	      int maxfev , int info[] , int nfev[] , double wa[] )
{

	double LBFGS_gtol = 0.1;
	double LBFGS_stpmin = 1e-16;
	double LBFGS_stpmax = 1e16;
	static int infoc[1];
	int j = 0;
	static double dg = 0, dgm = 0, dginit = 0, dgtest = 0;
	static double dgx[1];
	static double dgxm[1];
	static double dgy[1];
        static double dgym[1];
       	static double finit = 0, ftest1 = 0, fm = 0;
	static double fx[1];
	static double fxm[1];
	static double fy[1];
	static double fym[1];
	static double p5 = 0, p66 = 0;
	static double stx[1];
	static double sty[1];
	static double stmin = 0, stmax = 0, width = 0, width1 = 0, xtrapf = 0;
	static int brackt[1];
	static int stage1 = FALSE;

	p5 = 0.5;
	p66 = 0.66;
	xtrapf = 4;

	if ( info[0] != - 1 )
	{
		infoc[0] = 1;
		if ( 	n <= 0 || stp[0] <= 0 || ftol < 0 || 
			LBFGS_gtol < 0 || xtol < 0 || LBFGS_stpmin < 0 || 
			LBFGS_stpmax < LBFGS_stpmin || maxfev <= 0
		   ) 
			return;

		/* 
		 * Compute the initial gradient in the search direction
		 * and check that s is a descent direction.
		 */

		dginit = 0;

		for ( j = 0 ; j < n ; j++ )
		{
			dginit = dginit + g [j] * s [is0+j];
		}

		if ( dginit >= 0 )
		{
			printf(("The search direction is not a descent direction."));
			return;
		}

		brackt[0] = FALSE;
		stage1 = TRUE;
		nfev[0] = 0;
		finit = f;
		dgtest = ftol*dginit;
		width = LBFGS_stpmax - LBFGS_stpmin;
		width1 = width/p5;

		for ( j = 0 ; j < n ; j++ )
		{
			wa [j] = x [j];
		}

		/*
		 The variables stx, fx, dgx contain the values of the step,
		 function, and directional derivative at the best step.
		 The variables sty, fy, dgy contain the value of the step,
		 function, and derivative at the other endpoint of
		 the interval of uncertainty.
		 The variables stp, f, dg contain the values of the step,
		 function, and derivative at the current step.
		 */

		stx[0] = 0;
		fx[0] = finit;
		dgx[0] = dginit;
		sty[0] = 0;
		fy[0] = finit;
		dgy[0] = dginit;
	}

	while ( TRUE )
	{
		if ( info[0] != -1 )
		{
			/*
			 Set the minimum and maximum steps to correspond
			 to the present interval of uncertainty.
			*/

			if ( brackt[0] )
			{
				stmin = FMIN ( stx[0] , sty[0] );
				stmax = FMAX ( stx[0] , sty[0] );
			}
			else
			{
				stmin = stx[0];
				stmax = stp[0] + xtrapf * ( stp[0] - stx[0] );
			}

			/* Force the step to be within the bounds stpmax and stpmin.*/

			stp[0] = FMAX ( stp[0] , LBFGS_stpmin );
			stp[0] = FMIN ( stp[0] , LBFGS_stpmax );

			/* If an unusual termination is to occur then let
			   stp be the lowest point obtained so far.
			 */

			if ( 	( brackt[0] && ( stp[0] <= stmin || stp[0] >= stmax ) ) ||
			       	nfev[0] >= maxfev - 1 || infoc[0] == 0 || 
				( brackt[0] && stmax - stmin <= xtol * stmax )
			   )
				stp[0] = stx[0];

			/* Evaluate the function and gradient at stp
			   and compute the directional derivative.
			   We return to main program to obtain F and G.
			*/

			for ( j = 0 ; j < n ; j++ )
				x [j] = wa [j] + stp[0] * s [is0+j];

			info[0]=-1;
			return;
		}

		info[0]=0;
		nfev[0] = nfev[0] + 1;
		dg = 0;

		for ( j = 0 ; j < n ; j++ )
		{
			dg = dg + g [j] * s [is0+j];
		}

		ftest1 = finit + stp[0]*dgtest;

		/* Test for convergence.*/

		if ( 	( brackt[0] && ( stp[0] <= stmin || stp[0] >= stmax ) ) || infoc[0] == 0)
		       	info[0] = 6;

		if ( stp[0] == LBFGS_stpmax && f <= ftest1 && dg <= dgtest ) 
			info[0] = 5;

		if ( stp[0] == LBFGS_stpmin && ( f > ftest1 || dg >= dgtest ) ) 
			info[0] = 4;

		if ( nfev[0] >= maxfev )
		       	info[0] = 3;

		if ( brackt[0] && stmax - stmin <= xtol * stmax )
			info[0] = 2;

		if ( f <= ftest1 && fabs ( dg ) <= LBFGS_gtol * ( - dginit ) )
			info[0] = 1;

		/* Check for termination.*/

		if ( info[0] != 0 )
			return;

		/* In the first stage we seek a step for which the modified
		   function has a nonpositive value and nonnegative derivative.
		*/

		if ( stage1 && f <= ftest1 && dg >= FMIN ( ftol , LBFGS_gtol ) * dginit )
			stage1 = FALSE;

		/* 
		 * A modified function is used to predict the step only if
		   we have not obtained a step for which the modified
		   function has a nonpositive function value and nonnegative
		   derivative, and if a lower function value has been
		   obtained but the decrease is not sufficient.
		*/

		if ( stage1 && f <= fx[0] && f > ftest1 )
		{
			/* Define the modified function and derivative values.*/

			fm = f - stp[0]*dgtest;
			fxm[0] = fx[0] - stx[0]*dgtest;
			fym[0] = fy[0] - sty[0]*dgtest;
			dgm = dg - dgtest;
			dgxm[0] = dgx[0] - dgtest;
			dgym[0] = dgy[0] - dgtest;

			/* Call cstep to update the interval of uncertainty
			   and to compute the new step.
			*/

			mcstep ( stx , fxm , dgxm , sty , fym , dgym , stp , fm , dgm , 
					brackt , stmin , stmax , infoc );

			/* Reset the function and gradient values for f.*/

			fx[0] = fxm[0] + stx[0]*dgtest;
			fy[0] = fym[0] + sty[0]*dgtest;
			dgx[0] = dgxm[0] + dgtest;
			dgy[0] = dgym[0] + dgtest;
		}
		else
		{
			/* Call mcstep to update the interval of uncertainty
			   and to compute the new step.
			*/

			mcstep ( stx , fx , dgx , sty , fy , dgy , stp , f , dg ,
					brackt , stmin , stmax , infoc );
		}

		/* Force a sufficient decrease in the size of the
		   interval of uncertainty.
		*/

		if ( brackt[0] )
		{
			if ( fabs ( sty[0] - stx[0] ) >= p66 * width1 )
				stp[0] = stx[0] + p5 * ( sty[0] - stx[0] );
			width1 = width;
			width = fabs ( sty[0] - stx[0] );
		}
	}
}
/**************************************************************************/
static void arrayCopy(double* a,double*b,int n)
{
	int i;
	for(i=0;i<n;i++)
		b[i] = a[i];
}
/************************************************************************************************/
/** Compute the sum of a vector times a scalara plus another vector.
  * Adapted from the subroutine <code>daxpy</code> in <code>lbfgs.f</code>.
  * There could well be faster ways to carry out this operation; this
  * code is a straight translation from the Fortran.
  */ 
static void daxpy ( int n , double da , double dx[] , int ix0, int incx , double dy[] , int iy0, int incy )
{
	int i, ix, iy, m, mp1;

	if ( n <= 0 ) return;

	if ( da == 0 ) return;

	if  ( ! ( incx == 1 && incy == 1 ) )
	{
		ix = 1;
		iy = 1;

		if ( incx < 0 ) ix = ( - n + 1 ) * incx + 1;
		if ( incy < 0 ) iy = ( - n + 1 ) * incy + 1;

		for ( i = 0 ; i < n ; i++ )
		{
			dy [iy0+iy] = dy [iy0+iy] + da * dx [ix0+ix];
			ix = ix + incx;
			iy = iy + incy;
		}

		return;
	}

	m = n % 4;
	if ( m != 0 )
	{
		for ( i = 0 ; i < m ; i++ )
			dy [iy0+i] = dy [iy0+i] + da * dx [ix0+i];

		if ( n < 4 ) return;
	}

	mp1 = m + 1;
	for ( i = mp1-1 ; i < n ; i += 4 )
	{
		dy [iy0+i] = dy [iy0+i] + da * dx [ix0+i];
		dy [iy0+i + 1] = dy [iy0+i + 1] + da * dx [ix0+i + 1];
		dy [iy0+i + 2] = dy [iy0+i + 2] + da * dx [ix0+i + 2];
		dy [iy0+i + 3] = dy [iy0+i + 3] + da * dx [ix0+i + 3];
	}
	return;
}

/** Compute the dot product of two vectors.
  * Adapted from the subroutine <code>ddot</code> in <code>lbfgs.f</code>.
  * There could well be faster ways to carry out this operation; this
  * code is a straight translation from the Fortran.
  */ 
static double ddot ( int n, double dx[], int ix0, int incx, double dy[], int iy0, int incy )
{
	double dtemp;
	int i, ix, iy, m, mp1;

	dtemp = 0;

	if ( n <= 0 ) return 0;

	if ( !( incx == 1 && incy == 1 ) )
	{
		ix = 1;
		iy = 1;
		if ( incx < 0 ) ix = ( - n + 1 ) * incx + 1;
		if ( incy < 0 ) iy = ( - n + 1 ) * incy + 1;
		for ( i = 0 ; i < n ; i++ )
		{
			dtemp = dtemp + dx [ix0+ix] * dy [iy0+iy];
			ix = ix + incx;
			iy = iy + incy;
		}
		return dtemp;
	}

	m = n % 5;
	if ( m != 0 )
	{
		for ( i = 0 ; i < m ; i++ )
			dtemp = dtemp + dx [ix0+i] * dy [iy0+i];
		if ( n < 5 ) return dtemp;
	}

	mp1 = m + 1;
	for ( i = mp1-1 ; i < n ; i += 5 )
	{
		dtemp +=  dx [ix0+i] * dy [ iy0+i] 
			+ dx [ix0+i + 1] * dy [iy0+i + 1] 
			+ dx [ix0+i + 2] * dy [iy0+i + 2] 
			+ dx [ix0+i + 3] * dy [iy0+i + 3] 
			+ dx [ix0+i + 4] * dy [iy0+i + 4];
	}

	return dtemp;
}
/**************************************************************************/
static int lbfgs( 
		int n , int m , double x[] , double f , double g[] ,
		int diagco , double diag[] , double eps ,
		double xtol , int maxLines,int iflag[]
	)
{
	int execute_entire_while_loop = FALSE;
	static double gtol = 0.1;
	static double* solution_cache = NULL;
	static double gnorm = 0, stp1 = 0, ftol = 0;
	static double stp[1];
       	static double ys = 0, yy = 0, sq = 0, yr = 0, beta = 0, xnorm = 0;
	static int iter = 0, nfun = 0, point = 0, ispt = 0, iypt = 0, maxfev = 0;
	static int info[1];
	static int bound = 0, npt = 0, cp = 0, i = 0;
	static int nfev[1];
	static int inmc = 0, iycn = 0, iscn = 0;
	static int finish = FALSE;
	static double* w = NULL;
	static int wlength = 0;
	static int cacheLength = 0;


	if ( w == NULL || wlength != n*(2*m+1)+2*m )
	{
		if(w)
			free(w);

		wlength = n*(2*m+1)+2*m;
		w = malloc(wlength*sizeof(double));
	}
	if ( solution_cache == NULL || cacheLength != n )
	{
		if(solution_cache)
			free(solution_cache);

		cacheLength = n;
		solution_cache = malloc(cacheLength*sizeof(double));
	}

	if ( iflag[0] == 0 )
	{
		/* Initialize.*/

		arrayCopy(x,solution_cache,n);

		iter = 0;

		if ( n <= 0 || m <= 0 )
		{
			iflag[0]= -3;
			printf(("Improper input parameters  (n or m are not positive.)") );
		}

		if ( gtol <= 0.0001 )
		{
			printf(
				(
				"lbfgs: gtol is less than or equal to 0.0001."
				"It has been reset to 0.9."
				)
			      );
			gtol= 0.1;
		}

		nfun= 1;
		point= 0;
		finish= FALSE;
		if ( diagco )
		{
			for ( i = 0 ; i < n ; i++ )
			{
				if ( diag [i] <= 0 )
				{
					iflag[0]=-2;
					printf(
						(
						"The %d-th diagonal element of the inverse"
						" hessian approximation is not positive.")
						,i
					      );
				}
			}
		}
		else
		{
			for ( i = 0 ; i < n ; i++)
				diag [i] = 1;
		}
		ispt= n+2*m;
		iypt= ispt+n*m;

		for ( i = 0 ; i < n ; i++ )
			w [ispt + i] = - g [i] * diag [i];

		gnorm = sqrt ( ddot ( n , g , 0, 1 , g , 0, 1 ) );
		stp1= 1/gnorm;
		ftol= 0.0001; 
		maxfev= maxLines;

		execute_entire_while_loop = TRUE;
	}

	while ( TRUE )
	{
		if ( execute_entire_while_loop )
		{
			iter= iter+1;
			info[0]=0;
			bound=iter-1;
			if ( iter != 1 )
			{
				if ( iter > m ) bound = m;
				ys = ddot ( n , w , iypt + npt , 1 , w , ispt + npt , 1 );
				if ( ! diagco )
				{
					yy = ddot( 
						n , w , iypt + npt , 1 , w , iypt + npt , 1
						);
					for ( i = 0 ; i < n ; i++ )
						diag [i] = ys / yy;
				}
				else
				{
					iflag[0]=2;
					return 1;
				}
			}
		}

		if ( execute_entire_while_loop || iflag[0] == 2 )
		{
			if ( iter != 1 )
			{
				if ( diagco )
				{
					for ( i = 0 ; i < n ; i++ )
					{
						if ( diag [i] <= 0 )
						{
							iflag[0]=-2;
							printf(
							(
							"The %d-th diagonal element"
							" of the inverse hessian approximation"
							" is not positive.")
							, i);
						}
					}
				}
				cp= point;
				if ( point == 0 ) cp = m;
				w [ n + cp -1] = 1 / ys;

				for ( i = 0 ; i < n ; i++ )
					w [i] = - g [i];

				cp= point;

				for ( i = 0 ; i < bound ; i++ )
				{
					cp=cp-1;
					if ( cp == - 1 ) cp = m - 1;
					sq = ddot ( n , w , ispt + cp * n , 1 , w , 0 , 1 );
					inmc=n+m+cp+1;
					iycn=iypt+cp*n;
					w [ inmc -1] = w [ n + cp + 1 -1] * sq;
					daxpy ( n , - w [ inmc -1] , w , iycn , 1 , w , 0 , 1 );
				}

				for ( i = 0 ; i < n ; i++ )
					w [i] = diag [i] * w [i];

				for ( i = 0 ; i < bound ; i++ )
				{
					yr = ddot ( n , w , iypt + cp * n , 1 , w , 0 , 1 );
					beta = w [ n + cp + 1 -1] * yr;
					inmc=n+m+cp+1;
					beta = w [ inmc -1] - beta;
					iscn=ispt+cp*n;
					daxpy ( n , beta , w , iscn , 1 , w , 0 , 1 );
					cp=cp+1;
					if ( cp == m ) cp = 0;
				}

				for ( i = 0 ; i < n ; i++ )
					w [ispt + point * n + i] = w [i];
			}

			nfev[0]=0;
			stp[0]=1;
			if ( iter == 1 ) stp[0] = stp1;

			for ( i = 0 ; i < n ; i++ )
				w [i] = g [i];
		}

		mcsrch(
			n , x , f , g , w , ispt + point * n , stp ,
			ftol , xtol ,maxfev , info , nfev , diag
		      );

		if ( info[0] == - 1 )
		{
			iflag[0]=1;
			return 1;
		}

		if ( info[0] != 1 )
		{
			iflag[0]=-1;
			printf(
			(
			"Line search failed. See documentation of routine mcsrch.\n"
			" Error return of line search: info = %d Possible causes:\n"
			" function or gradient are incorrect, or incorrect tolerances.\n"
			)
			,info[0]);
			return 0;
		}

		nfun= nfun + nfev[0];
		npt=point*n;

		for ( i = 0; i < n ; i++ )
		{
			w [ispt + npt + i] = stp[0] * w [ispt + npt + i];
			w [iypt + npt + i] = g [i] - w [i];
		}

		point=point+1;
		if ( point == m ) point = 0;

		gnorm = sqrt ( ddot ( n , g , 0 , 1 , g , 0 , 1 ) );
		xnorm = sqrt ( ddot ( n , x , 0 , 1 , x , 0 , 1 ) );
		xnorm = FMAX ( 1.0 , xnorm );

		if ( gnorm / xnorm <= eps ) finish = TRUE;

		arrayCopy( x,solution_cache,n);

		if ( finish )
		{
			iflag[0]=0;
			return 0;
		}

		/* from now on, execute whole loop*/
		execute_entire_while_loop = TRUE;
	}
}
/*****************************************************************************************************************************************************/
void setQNQMOptions(FILE* file, QuasiNewtonQM* quasiNewtonQM)
{
	quasiNewtonQM->maxIterations = 100;
	quasiNewtonQM->updateFrequency = 1;
	quasiNewtonQM->epsilon  = 0.001;
	quasiNewtonQM->tolerence = 1e-16;  
	quasiNewtonQM->maxLines =  25;
	quasiNewtonQM->qmModel = NULL;
	quasiNewtonQM->maxStep = -1;// Angstrm, nolimit if maxStep<0
	readOneInt(file,"quasiNewtonMaxIterations",&quasiNewtonQM->maxIterations);
	readOneInt(file,"quasiNewtonUpdateFrequency",&quasiNewtonQM->updateFrequency);
	readOneReal(file,"quasiNewtonEpsilon",&quasiNewtonQM->epsilon);
	readOneReal(file,"quasiNewtonTolerence",&quasiNewtonQM->tolerence);
	readOneInt(file,"quasiNewtonMaxLines",&quasiNewtonQM->maxLines);
	readOneReal(file,"quasiNewtonMaxStep",&quasiNewtonQM->maxStep);
}
