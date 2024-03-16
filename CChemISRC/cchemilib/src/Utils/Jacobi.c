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

/* Jacobi.c */
#include <math.h>
#include <stdlib.h>
#define ROTATE(a,i,j,k,l) g=a[i][j];h=a[k][l];a[i][j]=g-s*(h+g*tau);\
	a[k][l]=h+s*(g-h*tau);
#define EPS 1e-10

int jacobi(double *M, int n, double d[], double **v, int *nrot)
{
	int j,iq,ip,i;
	double tresh,theta,tau,t,sm,s,h,g,c,*b,*z;
        double **a;
        int k,ki,imin;

	a=malloc(n*sizeof(double*));
        for(i=0;i<n;i++)
           a[i]=malloc(n*sizeof(double));
        iq = -1;
        for(i=0;i<n;i++)
         for(j=i;j<n;j++)
         {
          iq++;
          a[i][j] = M[iq];
         }
        for(i=0;i<n;i++)
         for(j=0;j<i;j++)
          a[i][j] = a[j][i];

	b=malloc(n*sizeof(double));
	z=malloc(n*sizeof(double));
	for (ip=0;ip<n;ip++) {
		for (iq=0;iq<n;iq++) v[ip][iq]=0.0;
		v[ip][ip]=1.0;
	}
	for (ip=0;ip<n;ip++) {
		b[ip]=d[ip]=a[ip][ip];
		z[ip]=0.0;
	}
	*nrot=0;
	for (i=0;i<50;i++) {
		sm=0.0;
		for (ip=0;ip<n-1;ip++) {
			for (iq=ip+1;iq<n;iq++)
				sm += fabs(a[ip][iq]);
		}
		if (fabs(sm)<=EPS) {
			free(z);
			free(b);
		        for(i=0;i<n;i++)
           			free(a[i]);
			free(a);

                        for(k=0;k<n-1;k++)
                        {
			  imin = k;
                          for(ki=k+1;ki<n;ki++)
				if(d[ki]<d[imin])
				   imin = ki;
                          if(imin != k)
			  {
			    sm = d[k];
                            d[k] = d[imin];
			    d[imin] = sm;

                            for(ki=0;ki<n;ki++)
                            {
				sm = v[ki][k];
				v[ki][k] = v[ki][imin];
				v[ki][imin] = sm ;
                            }

			  }
                        }

			return 0;
		}
		if (i < 4)
			tresh=0.2*sm/(n*n);
		else
			tresh=0.0;
		for (ip=0;ip<n-1;ip++) {
			for (iq=ip+1;iq<n;iq++) {
				g=100.0*fabs(a[ip][iq]);
				if (i > 4 && (double)(fabs(d[ip])+g) == (double)fabs(d[ip])
					&& (double)(fabs(d[iq])+g) == (double)fabs(d[iq]))
					a[ip][iq]=0.0;
				else if (fabs(a[ip][iq]) > tresh) {
					h=d[iq]-d[ip];
					if ((double)(fabs(h)+g) == (double)fabs(h))
						t=(a[ip][iq])/h;
					else {
						theta=0.5*h/(a[ip][iq]);
						t=1.0/(fabs(theta)+sqrt(1.0+theta*theta));
						if (theta < 0.0) t = -t;
					}
					c=1.0/sqrt(1+t*t);
					s=t*c;
					tau=s/(1.0+c);
					h=t*a[ip][iq];
					z[ip] -= h;
					z[iq] += h;
					d[ip] -= h;
					d[iq] += h;
					a[ip][iq]=0.0;
					for (j=0;j<=ip-1;j++) {
						ROTATE(a,j,ip,j,iq)
					}
					for (j=ip+1;j<=iq-1;j++) {
						ROTATE(a,ip,j,j,iq)
					}
					for (j=iq+1;j<n;j++) {
						ROTATE(a,ip,j,iq,j)
					}
					for (j=0;j<n;j++) {
						ROTATE(v,j,ip,j,iq)
					}
					++(*nrot);
				}
			}
		}
		for (ip=0;ip<n;ip++) {
			b[ip] += z[ip];
			d[ip]=b[ip];
			z[ip]=0.0;
		}
	}
	/*	Debug("Too many iterations in routine jacobi\n");*/
        free(z);
        free(b);
        for(i=0;i<n;i++)
        	free(a[i]);
	free(a);
        return 1;
}
#undef ROTATE
#undef EPS
