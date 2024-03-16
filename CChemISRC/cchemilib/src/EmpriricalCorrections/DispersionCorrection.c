/* DispersionCorrection.c */
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


/* Grimme at al, JCP, 132, 154104(2010) */

#ifndef OS_WIN32
#include <unistd.h>
#endif

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#ifdef ENABLE_MPI
#include <mpi.h>
#endif

#include "DispersionCorrection.h"
#include "../Utils/Constants.h"
#include "D3Parameters.h"

static double** r0ab = NULL;
static double***** c6ab = NULL;
static double* maxCi = NULL;
/**********************************************************************/
static void setParameters(DispersionParameters* parameters, char* method)
{
/* only the zero version(3) of D3 is implemented. BJ version is not yet implemented */
/* from J. Rezac, P. Hobza J. Chem. Theory Comput. 8, 141-151 (2012)*/
	double s6  =  1.0;
	double s8 = 0.0;
	double rs6 = 1.0;
	double rs8 = 1.0;
	double alp =  14.0;
	sprintf(parameters->method,"%s",method);
	if(strstr(method,"PM6"))
	{
               rs6=1.18;
               alp = 22.0;
               s6=0.88;
               s8=0.0;
	}
	else if(strstr(method,"SCC-DFTB"))
	{
               rs6=1.215;
               alp = 30.0;
               s8=0.0;
	}
	else if(strstr(method,"PM3"))
	{
               rs6=0.90;
               alp = 22.0;
               s8=0.0;
	}
	else if(strstr(method,"AM1"))
	{
               rs6=0.90;
               alp = 15.0;
               s8=0.0;
	}
	else if(strstr(method,"RM1"))
	{
               rs6=1.0;
               alp = 16.0;
               s8=0.0;
	}
	else if(strstr(method,"OM3"))
	{
               rs6=1.14;
               alp = 23.0;
               s8=0.0;
	}
	parameters->s6 = s6;
	parameters->rs6 = rs6;
	parameters->alp6 = alp;
	parameters->s8 = s8;
	parameters->rs8 = rs8;
	parameters->alp8 = alp + 2.0;
	parameters->k1 = 16.0;
	parameters->k2 = 4.0/3.0;
	parameters->k3 = -4.0;
	parameters->rthr = 20000;
}
/****************************************************************************************************/
static void setr0ab()
{
	int i,j,k;
	static const int mm = sizeof(r0ab1D)/sizeof(r0ab1D[0]);
	if(r0ab) return;
	r0ab = malloc(maxElements*sizeof(double*));
	k = 0;
	for(i=0;i<maxElements;i++) r0ab[i] = malloc(maxElements*sizeof(double));

	for(i=0;i<maxElements;i++)
	for(j=0;j<=i;j++)
	{
		r0ab[i][j] = r0ab1D[k]*ANGTOBOHR;
		r0ab[j][i] = r0ab[i][j];
		k++;
		if(k>=mm) break;
	}
}
/****************************************************************************************************/
static void scaleRCovalent(DispersionParameters* parameters)
{
	static int first = 1;
	int i;
	if(!first) return;
	for(i=0;i<maxElements;i++) RCovalent[i] *=parameters->k2*ANGTOBOHR;
	first = 0;
}
/****************************************************************************************************/
static void scaler2r4(DispersionParameters* parameters)
{
	static int first = 1;
	int i;
	if(!first) return;
	for(i=0;i<maxElements;i++) r2r4[i]=sqrt( 0.5*r2r4[i]*sqrt((i+1)*1.0));
	first = 0;
}
/****************************************************************************************************/
static void setc6abtable()
{
	int i,j,ia,ib,k;
	if(c6ab) return;
	maxCi = malloc(maxElements*sizeof(double));
	for(i=0;i<maxElements;i++) maxCi[i] = -1;
	static const int mm = sizeof(c6)/sizeof(c6[0]);

	c6ab = malloc(maxElements*sizeof(double****));
	k = 0;
	for(i=0;i<maxElements;i++) 
	{
		c6ab[i] = malloc(maxElements*sizeof(double***));
		for(j=0;j<maxElements;j++) 
		{
			c6ab[i][j] = malloc(maxCN*sizeof(double**));
			for(ia=0;ia<maxCN;ia++) 
			{
				c6ab[i][j][ia] = malloc(maxCN*sizeof(double*));
				for(ib=0;ib<maxCN;ib++) 
				{
					c6ab[i][j][ia][ib] = malloc(3*sizeof(double));
					for(k=0;k<3;k++)  c6ab[i][j][ia][ib][k] = -1;
				}
			}
		}
	}
	for(i=0;i<mm;i++)
	{
		int iat = iza[i]-1;
		int jat = izb[i]-1;
      		int iadr = 0;
      		int jadr = 0;
		while(iat>=100) {
			iat -=100;
			iadr++;
		}
		while(jat>=100) {
			jat -=100;
			jadr++;
		}
		if(iadr>=maxCN) continue;
		if(jadr>=maxCN) continue;
		if(maxCi[iat]<iadr) maxCi[iat]=iadr;
		if(maxCi[jat]<jadr) maxCi[jat]=jadr;
		c6ab[iat][jat][iadr][jadr][0] = c6[i];
		c6ab[iat][jat][iadr][jadr][1] = cna[i];
		c6ab[iat][jat][iadr][jadr][2] = cnb[i];
		c6ab[jat][iat][jadr][iadr][0] = c6[i];
		c6ab[jat][iat][jadr][iadr][1] = cnb[i];
		c6ab[jat][iat][jadr][iadr][2] = cna[i];
	}
}
/* interpolation of C6 */
static double getC6(DispersionParameters* parameters, int iat, int jat, double cni, double cnj)
{
	double c6mem = -1.0e90;
	double rsum = 0.0;
	double csum = 0.0;
	double c6 = 0.0;
	double cn1,cn2;
	double d;
	double t;
	int i,j;
	for(i=0;i<=maxCi[iat];i++)
	for(j=0;j<=maxCi[jat];j++)
	{
         	c6=c6ab[iat][jat][i][j][0];
         	if(c6>0)
		{
            		c6mem=c6;
            		cn1=c6ab[iat][jat][i][j][1];
            		cn2=c6ab[iat][jat][i][j][2];
            		d=(cn1-cni)*(cn1-cni)+(cn2-cnj)*(cn2-cnj);
            		t=exp(parameters->k3*d);
            		rsum=rsum+t;
            		csum=csum+t*c6;
		}
	}
      if(rsum>0) c6=csum/rsum;
      else { c6=c6mem;}
	return c6;
}
/*  compute coordination numbers by adding an inverse damping function*/
/* store values in cn table */
static void computeCN(Molecule* mol, DispersionParameters* parameters,double* cn)
{
	int iat;
	double d[3],r,damp,xn,rr;
	int i,iz,jz,k;
	double k1 = parameters->k1;
	for(i=0;i<mol->nAtoms;i++)
	{
		xn = 0;
		for(iat=0;iat<mol->nAtoms;iat++)
		{
			if(iat==i) continue;
			for(k=0;k<3;k++) d[k] = mol->atoms[iat].coordinates[k]-mol->atoms[i].coordinates[k];
			r = 0;
			for(k=0;k<3;k++) r+=d[k]*d[k];
			r = sqrt(r);
			r *= ANGTOBOHR;
			iz =  mol->atoms[i].prop.atomicNumber-1;
			jz =  mol->atoms[iat].prop.atomicNumber-1;
			rr = (RCovalent[iz]+RCovalent[jz])/r;
			damp  = 1./(1.+exp(-k1*(rr-1.0)));
            		xn=xn+damp;
		}
      		cn[i] = xn;
	}
}
/*  compute energy */
static double getEnergy(Molecule* mol, DispersionParameters* parameters)
{
	int iat,jat;
	int k;
	int iz,jz;
	double r,r2,r6,r8,tmp,c6,c8;
	double damp6,damp8,rr;
	double* cn = NULL;
	double d[3];
	double e6,e8;
	double rthr;
	double rs6, alp6, rs8, alp8, s6,s8;

	e6 =0;
	e8 =0;
	if(!mol) return 0;
	if(!parameters) return 0;
	if(mol->nAtoms<1) return 0;
	rthr = parameters->rthr;

	rs6 = parameters->rs6;
	rs8 = parameters->rs8;
	alp6 = parameters->alp6;
	alp8 = parameters->alp8;
	s6 = parameters->s6;
	s8 = parameters->s8;
	
/*
	printf("s6 = %f\n",s6);
	printf("rs6 = %f\n",rs6);
	printf("alpha6 = %f\n",alp6);
*/

	cn = malloc(mol->nAtoms*sizeof(double));
	computeCN(mol,parameters,cn);
	for(iat=0;iat<mol->nAtoms-1;iat++)
	for(jat=iat+1;jat<mol->nAtoms;jat++)
	{
		for(k=0;k<3;k++) 
			d[k] = mol->atoms[iat].coordinates[k]-mol->atoms[jat].coordinates[k];
		r2 = 0;
		for(k=0;k<3;k++) r2 +=d[k]*d[k];
		if(r2>rthr) continue;
		r = sqrt(r2);
		r *= ANGTOBOHR;
		r2 = r*r;
		iz =  mol->atoms[iat].prop.atomicNumber-1;
		jz =  mol->atoms[jat].prop.atomicNumber-1;
		rr = r0ab[iz][jz]/r;

		tmp=rs6*rr;
		damp6 =1./( 1.+6.*pow(tmp,alp6) );
		tmp=rs8*rr;
		damp8 =1./( 1.+6.*pow(tmp,alp8) );
		c6 = getC6(parameters, iz, jz, cn[iat], cn[jat]);
		r6=r2*r2*r2;
		r8=r6*r2;
		c8 =3.0*c6*r2r4[iz]*r2r4[jz];
		e6=e6+c6*damp6/r6;
		e8=e8+c8*damp8/r8;
	}
/*
	double x = 0;
	for(iat=0;iat<mol->nAtoms;iat++)
	for(jat=0;jat<mol->nAtoms;jat++)
	{
		iz =  mol->atoms[iat].prop.atomicNumber-1;
		jz =  mol->atoms[jat].prop.atomicNumber-1;
		x += getC6(parameters, iz, jz, cn[iat], cn[jat]);
	}
	printf("x = %f\n",x);
*/
	if(cn) free(cn);
	return -(s6*e6+s8*e8)*AUTOKCAL;
}
/*  compute energy and gradient */
static double computeEnergyAddGradient(Molecule* mol, DispersionParameters* parameters)
{
	double energy = getEnergy(mol,parameters);
	double step=2.e-5*BOHRTOANG;
	int i,j;
	double eplus,emoins;
	double f = 0.5/step;
	if(!mol) return 0;
      	for(i=0;i<mol->nAtoms;i++)
	{
      		for(j=0;j<3;j++) 
		{
			mol->atoms[i].coordinates[j] += step;
			eplus = getEnergy(mol,parameters);
			mol->atoms[i].coordinates[j] -= 2*step;
			emoins = getEnergy(mol,parameters);
      			mol->atoms[i].gradient[j]  += f*(eplus-emoins);
			mol->atoms[i].coordinates[j] += step;
		}
	}
	return energy;
}
/****************************************************************************************************/
double getD3Correction(Molecule* mol, char* method, boolean addGradient)
{
	double e = 0;
	DispersionParameters parameters;
	setr0ab();
	setc6abtable();
	setParameters(&parameters, method);
	scaleRCovalent(&parameters);
	scaler2r4(&parameters);
	if(addGradient) e = computeEnergyAddGradient(mol,&parameters);
	else e = getEnergy(mol,&parameters);
	return e;
}
