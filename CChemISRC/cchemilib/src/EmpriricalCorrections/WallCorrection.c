/* WallCorrection.c */
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

#include "WallCorrection.h"
#include "../Utils/Constants.h"


/**********************************************************************/
/*
static boolean isConnected(Molecule* mol, int j)
{
	int nc = 0;
	int k;
	if(!mol->atoms[j].typeConnections) return FALSE;
        for(k=0;k<mol->nAtoms;k++) if(mol->atoms[j].typeConnections[k]>0) nc++;
	if(nc>0) return TRUE;
	return FALSE;
}
*/
/**********************************************************************/
/*  compute energy */
static double getEnergy(Molecule* mol)
{

	double C[3] = {0,0,0};
	int iat,c;
	int nc = 0;
	double E0 = mol->wall.E0;
	double srho2 =1/(mol->wall.rho*mol->wall.rho);
	int n = mol->wall.n;
	double E = 0;
	for(iat=0;iat<mol->nAtoms;iat++)
	{
		if(mol->atoms[iat].residueNumber==0)
		{
			nc++;
			for(c=0;c<3;c++) C[c]+= mol->atoms[iat].coordinates[c];
		}
	}
	if(nc>0) for(c=0;c<3;c++) C[c] /= nc;
	//printf("nc=%d\n",nc);
	for(iat=0;iat<mol->nAtoms;iat++)
	{
		if(mol->atoms[iat].residueNumber!=0)
		{
			double r2 = 0;
			for(c=0;c<3;c++) r2 +=  (mol->atoms[iat].coordinates[c]-C[c])*(mol->atoms[iat].coordinates[c]-C[c]);
			E += E0*pow(1-exp(-r2*srho2),n);
		}
	}
	//printf("EWall=%f\n",E);

	return E;
}
/*  compute energy and gradient */
/*
static double computeEnergyAddGradientNumeric(Molecule* mol)
{
	double energy = getEnergy(mol);
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
			eplus = getEnergy(mol);
			mol->atoms[i].coordinates[j] -= 2*step;
			emoins = getEnergy(mol);
      			mol->atoms[i].gradient[j]  += f*(eplus-emoins);
			mol->atoms[i].coordinates[j] += step;
		}
	}
	return energy;
}
*/
/*  compute energy and gradient analytically */
static double computeEnergyAddGradientAnalytic(Molecule* mol)
{
	double C[3] = {0,0,0};
	int iat,c;
	int nc = 0;
	double E0 = mol->wall.E0;
	double srho2 =1/(mol->wall.rho*mol->wall.rho);
	int n = mol->wall.n;
	double energy = 0;
	for(iat=0;iat<mol->nAtoms;iat++)
	{
		if(mol->atoms[iat].residueNumber==0)
		{
			nc++;
			for(c=0;c<3;c++) C[c]+= mol->atoms[iat].coordinates[c];
		}
	}
	if(nc>0) for(c=0;c<3;c++) C[c] /= nc;
	for(iat=0;iat<mol->nAtoms;iat++)
	{
		if(mol->atoms[iat].residueNumber!=0)
		{
			double r2 = 0;
			double ex = 0;
			double exn = 0;
			for(c=0;c<3;c++) r2 +=  (mol->atoms[iat].coordinates[c]-C[c])*(mol->atoms[iat].coordinates[c]-C[c]);
			if( r2> mol->wall.rho*mol->wall.rho) printf("r=%f\n",sqrt(r2));
			ex = exp(-r2*srho2);
			exn = pow(1-ex,n-1);
			energy += E0*exn*(1-ex);
			exn  = 2*n*srho2*exn*ex*E0;
			for(c=0;c<3;c++) mol->atoms[iat].gradient[c] += (mol->atoms[iat].coordinates[c]-C[c])*exn;
		}
	}
	//printf("EWall=%f\n",energy);
	return energy;

}
/****************************************************************************************************/
double getWallCorrection(Molecule* mol, boolean addGradient)
{
	double e = 0;
	//printf("E0=%f\n",mol->wall.E0);
	//printf("rho=%f\n",mol->wall.rho);
	//printf("n=%d\n",mol->wall.n);
	if(mol->wall.E0<0 || fabs(mol->wall.E0)<1e-10) return 0;
	//if(addGradient) e = computeEnergyAddGradientNumeric(mol);
	if(addGradient) e = computeEnergyAddGradientAnalytic(mol);
	else e = getEnergy(mol);
	return e;
}
