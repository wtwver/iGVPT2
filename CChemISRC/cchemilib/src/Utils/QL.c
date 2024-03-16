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

/* QL.c */
/* To determine the eigenvalues and eigenvectors of a real, symmetric matrix A using the 
 * QL algorithm to determine after tridiagonalisation */

#include <math.h>
#include <stdlib.h>
#include <stdio.h>

/********************************************************************************/
static void reductionToTridiagonal(double **A, int n, double *D, double *E);
static int diagonalisationOfATridiagonalMatrix(double *D, double *E, int n, double **V);
/********************************************************************************/
static void swap2d(double* a, double* b)
{
        double c = *a;
        *a = *b;
        *b = c;
}
/********************************************************************************/
static void sort(int n, double *EVals, double** V)
{
        int i,j,k;

        for(i=0;i<n;i++)
        {
                k = i;
                for(j=i+1;j<n;j++) if(EVals[j]<EVals[k]) k = j;
                if(k==i) continue;
                swap2d(&EVals[i],&EVals[k]);
                for(j=0;j<n;j++) swap2d(&V[j][i],&V[j][k]);
        }
}
/********************************************************************************/
int eigen(double **M, int n, double *EVals, double** V)
{
	double** A;
	double* E;
	int success = 0;
	int i;
	int j;

	if(n<1) return 0;
	A = malloc(n*sizeof(double*));
	for(i=0;i<n;i++) A[i]=malloc(n*sizeof(double));

	for(i=0;i<n;i++)
  	for(j=i;j<n;j++)
    		A[i][j] = A[j][i] = M[i][j];

	E=malloc(n*sizeof(double));
	reductionToTridiagonal(A, n, EVals, E);
	/*
	for(i=0;i<n;i++) printf("EVals[%d]=%f\n",i,EVals[i]);
	*/
	success = diagonalisationOfATridiagonalMatrix(EVals, E, n, A);
	for(i=0;i<n;i++)
	for(j=0;j<n;j++)
		V[i][j] = A[i][j];

	free(E);
	for(i=0;i<n;i++) free(A[i]);
	free(A);
	sort(n, EVals, V);

	return success;
}
/********************************************************************************/
int eigenQL(int n, double *M, double *EVals, double** V)
{
	double** A;
	double* E;
	int ii;
	int success = 0;
	int i;
	int j;

	if(n<1) return 0;
	A = malloc(n*sizeof(double*));
	for(i=0;i<n;i++) A[i]=malloc(n*sizeof(double));

	/* M is an inf symmetric matrix */
	ii = -1;
	for(i=0;i<n;i++)
	for(j=0;j<=i;j++)
	{
		ii++;
		A[i][j] = M[ii];
	}
	for(i=0;i<n;i++)
  	for(j=i+1;j<n;j++)
    		A[i][j] = A[j][i];

	E=malloc(n*sizeof(double));
	reductionToTridiagonal(A, n, EVals, E);
	/*
	for(i=0;i<n;i++) printf("EVals[%d]=%f\n",i,EVals[i]);
	*/
	success = diagonalisationOfATridiagonalMatrix(EVals, E, n, A);
	for(i=0;i<n;i++)
	for(j=0;j<n;j++)
		V[i][j] = A[i][j];

	free(E);
	for(i=0;i<n;i++) free(A[i]);
	free(A);

	return success;
}
/* procedure to reduce a real symmetric matrix to the tridiagonal form that is suitable for input to 
 * diagonalisationOfATridiagonalMatrix.*/
/********************************************************************************/
static void reductionToTridiagonal(double **A, int n, double *D, double *E)
{
	int	l, k, j, i;
	double  scale, hh, h, g, f;
 
	for (i = n-1; i >= 1; i--)
	{
	    l = i - 1;
	    h = scale = 0.0;
	    if (l > 0)
	    {
		   for (k = 0; k <= l; k++) scale += fabs(A[i][k]);
		   if (scale == 0.0) E[i] = A[i][l];
		   else
		   {
			  for (k = 0; k <= l; k++)
			  {
				 A[i][k] /= scale;
				 h += A[i][k] * A[i][k];
			  }
			  f = A[i][l];
			  g = f > 0 ? -sqrt(h) : sqrt(h);
			  E[i] = scale * g;
			  h -= f * g;
			  A[i][l] = f - g;
			  f = 0.0;
			  for (j = 0; j <= l; j++)
			  {
				 A[j][i] = A[i][j] / h;
				 g = 0.0;
				 for (k = 0; k <= j; k++) g += A[j][k] * A[i][k];
				 for (k = j + 1; k <= l; k++) g += A[k][j] * A[i][k];
				 E[j] = g / h;
				 f += E[j] * A[i][j];
			  }
			  hh = f / (h + h);
			  for (j = 0; j <= l; j++)
			  {
				 f = A[i][j];
				 E[j] = g = E[j] - hh * f;
				 for (k = 0; k <= j; k++) A[j][k] -= (f * E[k] + g * A[i][k]);
			  }
		   }
	    } else E[i] = A[i][l];
	    D[i] = h;
	}
	D[0] = 0.0;
	E[0] = 0.0;
	for (i = 0; i < n; i++)
	{
	    l = i - 1;
	    if (D[i])
	    {
		   for (j = 0; j <= l; j++)
		   {
			  g = 0.0;
			  for (k = 0; k <= l; k++) g += A[i][k] * A[k][j];
			  for (k = 0; k <= l; k++) A[k][j] -= g * A[k][i];
		   }
	    }
	    D[i] = A[i][i];
	    A[i][i] = 1.0;
	    for (j = 0; j <= l; j++) A[j][i] = A[i][j] = 0.0;
	}
}
#undef SIGN
#define SIGN(A,B) ((B)<0 ? -fabs(A) : fabs(A))
/* QL algorithm to determine 
 * the eigenvalues and eigenvectors of a real, symmetric, tridiagonal matrix.*/
/********************************************************************************/
static int diagonalisationOfATridiagonalMatrix(double *D, double *E, int n, double **V)
{
	int	m, l, iter, i, k;
	double  s, r, p, g, f, dd, c, b;
 
	for (i = 1; i < n; i++) E[i - 1] = E[i];
	E[n-1] = 0.0;
	for (l = 0; l < n; l++)
	{
	    iter = 0;
	    do
	    {
		   for (m = l; m < n - 1; m++)
		   {
			  dd = fabs(D[m]) + fabs(D[m + 1]);
			  if (fabs(E[m]) + dd == dd) break;
		   }
		   if (m != l)
		   {
			  if (iter++ == 30) return 0;
			  g = (D[l + 1] - D[l]) / (2.0 * E[l]);
			  r = sqrt((g*g) + 1.0);
			  g = D[m] - D[l] + E[l] / (g + SIGN(r, g));
			  s = c = 1.0;
			  p = 0.0;
			  for (i = m - 1; i >= l; i--)
			  {
				 f = s * E[i];
				 b = c * E[i];
				 if (fabs(f) >= fabs(g))
				 {
					c = g / f;
					r = sqrt((c*c) + 1.0);
					E[i + 1] = f * r;
					c *= (s = 1.0 / r);
				 } else
				 {
					s = f / g;
					r = sqrt((s*s) + 1.0);
					E[i + 1] = g * r;
					s *= (c = 1.0 / r);
				 }
				 g = D[i + 1] - p;
				 r = (D[i] - g) * s + 2.0 * c * b;
				 p = s * r;
				 D[i + 1] = g + p;
				 g = c * r - b;
				 for (k = 0; k < n; k++)
				 {
					f = V[k][i + 1];
					V[k][i + 1] = s * V[k][i] + c * f;
					V[k][i] = c * V[k][i] - s * f;
				 }
			  }
			  D[l] = D[l] - p;
			  E[l] = g;
			  E[m] = 0.0;
		   }
	    } while (m != l);
	}
 
	return 1;
}
/*******************************************************************/
/* Given a matrix
a[0..n-1][0..n-1] , this routine replaces it by the LU decomposition of a rowwise
permutation of itself.  a and n are  input.  a is  output,
indx[0..n-1] is  an  output  vector  that  records  the  row  permutation  ected  by  the  partial pivoting;
d is  output  as  1 depending  on  whether  the  number  of  row  interchanges was  even
or odd, respectively.  This routine is used in combination with
lubksb to solve linear equations or  invert  a  matrix
*/
void ludcmp(double **a, int n, int *indx, double *d)
{
	int i,imax=0,j,k;
	double big,dum,sum,temp;
	double *vv;

	vv= malloc(n*sizeof(double));
	*d=1.0;
	for (i=0;i<n;i++) {
		big=0.0;
		for (j=0;j<n;j++)
			if ((temp=fabs(a[i][j])) > big) big=temp;
		if (big == 0.0) fprintf(stderr,"Singular matrix in routine ludcmp");
		vv[i]=1.0/big;
	}
	for (j=0;j<n;j++) {
		for (i=0;i<j;i++) {
			sum=a[i][j];
			for (k=0;k<i;k++) sum -= a[i][k]*a[k][j];
			a[i][j]=sum;
		}
		big=0.0;
		for (i=j;i<n;i++) {
			sum=a[i][j];
			for (k=0;k<j;k++)
				sum -= a[i][k]*a[k][j];
			a[i][j]=sum;
			if ( (dum=vv[i]*fabs(sum)) >= big) {
				big=dum;
				imax=i;
			}
		}
		if (j != imax) {
			for (k=0;k<n;k++) {
				dum=a[imax][k];
				a[imax][k]=a[j][k];
				a[j][k]=dum;
			}
			*d = -(*d);
			vv[imax]=vv[j];
		}
		indx[j]=imax;
		if (a[j][j] == 0.0) a[j][j]=1e-20;
		if (j != n-1) {
			dum=1.0/(a[j][j]);
			for (i=j+1;i<n;i++) a[i][j] *= dum;
		}
	}
	free(vv);

}
/******************************************************/
/* Solves the set of n linear equations A X = B .Here a[0..n-1][0..n-1]
is input, not as the matrix A but rather as its LU decomposition, determined by the routine
ludcmp .  indx[0..n-1] is input as the permutation vector returned by ludcmp.
b[0..n-1] is input as the right-hand side vector B ,  and  returns  with  the  solution  vector
X .  a , n ,and indx are  not  modied  by  this  routine and can be left in place for 
successive calls with different right-hand sides b .  This routine takes into  account the  possibility  that
 b will  begin with  many zero  elements, so it  is  elsecient for use in  matrix  inversion */
void lubksb(double **a, int n, int *indx, double b[])
{
	int i,ii=-1,ip,j;
	double sum;

	for (i=0;i<n;i++) {
		ip=indx[i];
		sum=b[ip];
		b[ip]=b[i];
		if (ii!=-1) for (j=ii;j<=i-1;j++) sum -= a[i][j]*b[j];
		else if (sum) ii=i;
		b[i]=sum;
	}
	for (i=n-1;i>=0;i--) {
		sum=b[i];
		for (j=i+1;j<n;j++) sum -= a[i][j]*b[j];
		b[i]=sum/a[i][i];
	}
}
/********************************************************************************/
int solveSymEqQL(int n, double **A, double* B, double* values)
{
	double** V;
	int *indx;
	int i,j;
	int success = 0;
	double d;
	V = malloc(n*sizeof(double*));
	for(i=0;i<n;i++) V[i]=malloc(n*sizeof(double));

	for(i=0;i<n;i++)
	for(j=0;j<=i;j++)
	{
		V[i][j] = V[j][i] = A[i][j];
	}
	indx = malloc(n*sizeof(int));
	for(i=0;i<n;i++) values[i] = B[i]; 
	/*
	printf("A\n");
	for(i=0;i<n;i++) 
	{
		for(j=0;j<n;j++) printf("%f ",V[i][j]);
		printf("\n");
	}
	printf("B\n");
	for(i=0;i<n;i++) printf("%f ",B[i]);
	printf("\n");
	*/

	ludcmp(V,n,indx,&d);
	lubksb(V,n,indx,values);

	/*
	printf("values\n");
	for(i=0;i<n;i++) printf("%f ",values[i]);
	printf("\n");
	*/
	for(i=0;i<n;i++) free(V[i]);
	free(V);
	free(indx);
	return success;
}
