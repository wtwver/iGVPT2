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

/* Refs 
Jean-Louis Barrat , David Rodney
Portable implementation of a quantum thermal bath for molecular dynamics simulations
JOURNAL OF STATISTICAL PHYSICS 670, 144, (2011)
*/
static void updateQTB(Paths* paths);
static void resetQTB(Paths* paths);
static void computeThetaQTB(Paths* paths);
static void applyQTBThermostat(Paths* paths);
//static double gP(int P, double x);
//static double* compute_gP(int P, double xmax, int n, int ng);
static void removeTranslationAndRotationTheta(Paths* paths);
static void removeTranslationTheta(Paths* paths);
static void removeRotationTheta(Paths* paths);
//static void removeTranslationAndRotationThetaOneBead(Paths* paths, int i);
//static void removeTranslationThetaOneBead(Paths* paths, int i);
//static void removeRotationThetaOneBead(Paths* paths, int i);
static double* compute_gPTable(int P, double dx, int M);
static double interpol(int M, double* g, double dx, double x);
//static double* compute_fPTable(int P, double dx, int M);
/*********************************************************************************/
/* omegaMax in cm-1 */
void initQTBThermostat(Paths* paths)
{
	static double cmM1fsM1 = 2.99792458e-5;
	double Omegafs = paths->omegaMax*cmM1fsM1;/* fs^-1 */ 
	int nAtoms = paths->molecules[0]->nAtoms;
	
		fprintf(stdout,"gamma(ps^-1)\t\t\t= %f\n",paths->friction*fsInAKMA*1000);
	if(paths->Nf<1) paths->Nf = 50;

	//Omegafs *= sqrt(paths->nBeads);
	printf("wp(cm-1)=%f\n",paths->kT*2*paths->nBeads*349.75511054);
	//Omegafs += (paths->kT*2*paths->nBeads*349.75511054*cmM1fsM1);

	paths->h = 1/Omegafs*(fsInAKMA);
	//paths->h = paths->dt;
	paths->nQTBSteps = (int)(paths->h/paths->dt);
	if(paths->nQTBSteps<1) paths->nQTBSteps = 1;
	paths->h = paths->nQTBSteps *paths->dt;
	paths->omegaMax = 1.0/paths->h*(fsInAKMA)/cmM1fsM1; /* cm-1 */
	if(paths->friction<0) paths->friction = (1.0/ paths->h)/50;
		fprintf(stdout,"gamma(ps^-1)\t\t\t= %f\n",paths->friction*fsInAKMA*1000);
	
	paths->Ht = newMatrixDouble(paths->nBeads,2*paths->Nf);
	paths->theta = newMatrixDouble(paths->nBeads,3*nAtoms);
	paths->rnoise = newCubeDouble(paths->nBeads,3*nAtoms,2*paths->Nf);

	paths->gpQTB = 1/(1+paths->friction*paths->dt/2);
	paths->gmQTB = (1-paths->friction*paths->dt/2)*paths->gpQTB;

	printf("gpQTB=%f gmQTB=%f\n",paths->gpQTB,paths->gmQTB);
/*
	paths->gpQTB = 1.0;
	paths->gmQTB = (1-paths->friction*paths->dt/2)*paths->gpQTB;
	printf("gpQTB=%f gmQTB=%f\n",paths->gpQTB,paths->gmQTB);
*/
	//paths->gmQTB = exp(-paths->dt/2.0*paths->friction);
	//paths->gpQTB = 1/(1+paths->friction*paths->dt/2);
	//paths->gmQTB *= paths->gpQTB;
	//printf("gpQTB=%f gmQTB=%f\n",paths->gpQTB,paths->gmQTB);

	paths->iQTBStep = 0;

	paths->klass->applyThermostat = applyQTBThermostat;

	resetQTB(paths);
	printf("End initQTB\n");

}
/*****************************************************************************************************************************/
static void applyQTBThermostat(Paths* paths)
{
	int i,k,iAtom;
	int nAtoms = paths->molecules[0]->nAtoms;
	double*** P = paths->P;
	double*** F = paths->F;

	if((paths->iQTBStep+1)==2*paths->nQTBSteps){ updateQTB(paths); paths->iQTBStep = 0;}
	else paths->iQTBStep++;

	for(i = 0;i<paths->nBeads; i++) 
	{
	for(iAtom = 0;iAtom<nAtoms; iAtom++) 
	for(k = 0;k<3; k++) 
	{
		if((paths->iQTBStep+1)%2==0) F[k][i][iAtom] +=  paths->gpQTB*paths->theta[i][3*iAtom+k]-4*(1-paths->gmQTB)/(1+paths->gmQTB)/paths->dt*P[k][i][iAtom];
		else P[k][i][iAtom] +=  paths->gpQTB*paths->theta[i][3*iAtom+k]*0.5*paths->dt; 

/*
		if((paths->iQTBStep+1)%2==0)
			F[k][i][iAtom] =  paths->gpQTB*(F[k][i][iAtom]+paths->theta[i][3*iAtom+k])+2*(paths->gmQTB-1)/paths->dt*P[k][i][iAtom];
		else
			P[k][i][iAtom] +=  paths->gpQTB*paths->theta[i][3*iAtom+k]*0.5*paths->dt
			+  (paths->gpQTB-1)*F[k][i][iAtom]*0.5*paths->dt; 
*/
	}
	//paths->molecules[i]->klass->removeTranslationAndRotation(paths->molecules[i]);
	}

/*
	if((paths->iQTBStep+1)==paths->nQTBSteps){ updateQTB(paths); paths->iQTBStep = 0;}
	else paths->iQTBStep++;

	Molecule** mols = paths->molecules;
	for(i = 0;i<paths->nBeads; i++) 
	for(iAtom = 0;iAtom<paths->molecules[i]->nAtoms; iAtom++) 
	for(k = 0;k<3; k++) 
	{
		mols[i]->atoms[iAtom].coordinates[k] +=  P[k][i][iAtom]*paths->dt/paths->Mprim[i][iAtom]
		+ (F[k][i][iAtom]+paths->gpQTB*paths->theta[i][3*iAtom+k]-paths->friction*P[k][i][iAtom])
		   *paths->dt*paths->dt/2/paths->Mprim[i][iAtom];

		P[k][i][iAtom]  =  P[k][i][iAtom]*paths->gmQTB;
		P[k][i][iAtom] +=  paths->gpQTB*paths->theta[i][3*iAtom+k]*paths->dt/2;
		P[k][i][iAtom] +=  paths->gpQTB*paths->F[k][i][iAtom]*paths->dt/2;
	}
	paths->klass->calculateForces(paths);

	for(i = 0;i<paths->nBeads; i++) 
	for(iAtom = 0;iAtom<paths->molecules[i]->nAtoms; iAtom++) 
	for(k = 0;k<3; k++) 
	{
		P[k][i][iAtom] +=  paths->gpQTB*paths->theta[i][3*iAtom+k]*paths->dt/2;
		P[k][i][iAtom] +=  paths->gpQTB*paths->F[k][i][iAtom]*paths->dt/2;
	}
	//paths->klass->removeTranslationAndRotation(paths);
	//paths->klass->removeTranslation(paths);
	// used to compute energy
*/
/*
	for(i = 0;i<paths->nBeads; i++) 
	for(iAtom = 0;iAtom<paths->molecules[i]->nAtoms; iAtom++) 
	for(k = 0;k<3; k++) 
	{
		mols[i]->atoms[iAtom].gradient[k] -=  paths->gpQTB*paths->theta[i][3*iAtom+k];
	}
*/
}
/*********************************************************************************/
static void updateQTB(Paths* paths)
{
	int i,iAtom,j,k;
	int nAtoms = paths->molecules[0]->nAtoms;

	if(paths->temperature<=0) return;
	computeThetaQTB(paths);
	/* shift rnoise */

	for(i=0;i<paths->nBeads;i++)
	for(iAtom=0;iAtom<nAtoms;iAtom++)
	for(j=0;j<3;j++)
	for(k=0;k<2*paths->Nf-1;k++)
		 paths->rnoise[i][3*iAtom+j][k] = paths->rnoise[i][3*iAtom+j][k+1];

	/* add one value to the end */
	for(i=0;i<paths->nBeads;i++)
	for(iAtom=0;iAtom<nAtoms;iAtom++)
	for(j=0;j<3;j++)
		 paths->rnoise[i][3*iAtom+j][2*paths->Nf-1] = normal(); /* sqrt(h) in sigma */

}
/*********************************************************************************/
static void computeThetaQTB(Paths* paths)
{
	//double sigma=sqrt(2.*paths->friction*paths->kT/paths->h);
	double sigma=sqrt(2*paths->friction*paths->kT/paths->h);
	int i,iAtom,j,k;
	int nAtoms = paths->molecules[0]->nAtoms;

	/* thetaProg = thetaPaper*sigma/sqrt(m) */
	for(i=0;i<paths->nBeads;i++)
	for(iAtom=0;iAtom<nAtoms;iAtom++)
	for(j=0;j<3;j++)
	{
		paths->theta[i][3*iAtom+j] = 0.0;
		for(k=0;k<2*paths->Nf;k++)
			paths->theta[i][3*iAtom+j] += paths->rnoise[i][3*iAtom+j][2*paths->Nf-1-k]*paths->Ht[i][k];
			//paths->theta[i][3*iAtom+j] += paths->rnoise[0][3*iAtom+j][2*paths->Nf-1-k]*paths->Ht[i][k];
		paths->theta[i][3*iAtom+j] *= sigma*sqrt(paths->Mprim[i][iAtom]);
	}
	//for(i=0;i<paths->nBeads;i++)
	//	paths->molecules[i]->klass->removeTranslationAndRotationVectors(paths->molecules[i],paths->theta[i]);
	//paths->molecules[0]->klass->removeTranslationVectorsCluster(paths->molecules,paths->nBeads,paths->theta);
	removeTranslationAndRotationTheta(paths);
	//removeTranslationTheta(paths);
	//for(i=0;i<paths->nBeads;i++)
		//removeTranslationAndRotationThetaOneBead(paths, i);
		//removeTranslationThetaOneBead(paths, i);
	//removeTranslationTheta(paths);
	//removeRotationTheta(paths);
	//for(i=0;i<paths->nBeads;i++)
		//paths->molecules[i]->klass->removeRotationForce(paths->molecules[i],paths->theta[i]);
		//paths->molecules[i]->klass->removeTranslationForce(paths->molecules[i],paths->theta[i]);
}
/*********************************************************************************/
static void resetQTB(Paths* paths)
{
	double* Filter = NULL;
	int k;
	double hbarwOverkT;
	double hbarwOverkTk;
	double hbardwOverkT;
	int i,iAtom,j;
	double T;
	int nAtoms = paths->molecules[0]->nAtoms;
	double* gp;
	double dx;
	int M;
	double sini2;
	double x2xk2;
	double xk;
	double hwover2;

	paths->iQTBStep = 0;
	T = paths->temperature;
	if(T<=0) return;

	/* computing of Filter */
	/* Htild/sqrt(kT), sqrt(kT) in sigma */
	Filter = malloc((2*paths->Nf)*sizeof(double)); 
	/* h dOmega = pi /Nf */

	//hbardwOverkT = 1.0/(paths->Nf*paths->h*paths->kT);
	//hbardwOverkT = 1.0/(paths->Nf*paths->h*paths->kT*(paths->nBeads));
	hbardwOverkT = 1.0/(paths->Nf*paths->h*paths->kT);
	hbarwOverkTk = sqrt(hbardwOverkT*paths->Nf*hbardwOverkT*paths->Nf+4*paths->nBeads);
	M =  paths->Nf*10*paths->nBeads;
	dx = hbarwOverkTk/2/(M-1);
	gp = compute_gPTable(paths->nBeads, dx, M);
	//gp = compute_fPTable(paths->nBeads, dx, M);
	FILE* f =fopen("gP.txt","w");
	for(k=0;k<M;k++) fprintf(f,"%f %f\n",dx*k, gp[k]);
	fclose(f);
	f =fopen("F.txt","w");
	for(i=0;i<paths->nBeads;i++)
	{
		sini2 = sin(i*M_PI/paths->nBeads);
		sini2 = sini2*sini2;
		sini2 *= 4*paths->nBeads;
		sini2 *= paths->nBeads;
	for(k=0;k<2*paths->Nf;k++)
	{
		int kk= k-paths->Nf;
		int kp= abs(kk);
		printf("k=%d\n",k);
		if(kk==0) Filter[k] = 1.0;
		else
		{
			hbarwOverkT = kp*hbardwOverkT;
			hbarwOverkTk = sqrt(hbarwOverkT*hbarwOverkT+sini2);
			x2xk2 = hbarwOverkT/hbarwOverkTk;
			x2xk2 = x2xk2*x2xk2;
			xk = hbarwOverkTk/2;
			// gpTable
			Filter[k] = sqrt(interpol(M, gp, dx,xk)*x2xk2);
			// fPTable
			//Filter[k] = (interpol(M, gp, dx,xk)*sqrt(x2xk2));
			fprintf(f,"%f %f %f\n", xk, Filter[k], Filter[k]-1);
			//Filter[k] = 1.0;
			hwover2 = xk*paths->kT*paths->h*M_PI;
			//hwover2 = kk*M_PI/paths->Nf/2;
			Filter[k] *= (hwover2)/sin(hwover2);
			//Filter[k] *= (kk*M_PI/paths->Nf/2)/sin(kk*M_PI/paths->Nf/2);
		}
	}
	//for(k=0;k<2*paths->Nf;k++)  Filter[k] *= (sqrt(paths->nBeads));
	/* compute Ht */
	for(j=0;j<2*paths->Nf;j++)
	{
                paths->Ht[i][j] = 0;
		for(k=0; k<2*paths->Nf;k++)
                	paths->Ht[i][j] += Filter[k]*cos(M_PI*(k-paths->Nf)*(j-paths->Nf)*1.0/paths->Nf);
		
               	paths->Ht[i][j] /= 2*paths->Nf;
	}
	}
	fclose(f);
	free(Filter);
	free(gp);
	
	for(i=0;i<paths->nBeads;i++)
	for(iAtom=0;iAtom<nAtoms;iAtom++)
	for(j=0;j<3;j++)
	for(k=0;k<2*paths->Nf;k++)
		 paths->rnoise[i][3*iAtom+j][k] = normal();/* sqrt(h) is in sigma */

	computeThetaQTB(paths);

}
/*********************************************************************************/
static double h(double x)
{
	return x/tanh(x);
}
/*********************************************************************************/
/* Ref ! eq. 17, Ceriotti et al., J. Chem. Phys. 134, 084104 (2011)*/
/*
static double gP(int P, double x)
{
	static double* gold  = NULL;
	static double* g  = NULL;
	static double* xk  = NULL;
	static double* hk  = NULL;
	static int Pold = 0;
	double s;
	double sink;
	int k;
	double alpha = 1.0/P;
	double precision = 1e-8;
	int nMax = 10000;
	int n = 0;
	double norm;
	int j;
	double xx;

	return h(x/P);

	if(Pold ==0 || Pold != P) 
	{
		if(gold) free(gold);
		if(g) free(g);
		if(xk) free(xk);
		if(hk) free(hk);
		gold = malloc(P*sizeof(double));
		g = malloc(P*sizeof(double));
		xk = malloc(P*sizeof(double));
		hk = malloc(P*sizeof(double));
		Pold = P;
	}
	//return h(x/P);
	// i = 0, g = h for x/P, xk/P, k = 1..P-1
	printf("x = %f\n",x);
	xk[0] = x;
	for(k=1;k<P;k++)
	{
		sink  = P*sin(k*M_PI*alpha);
		xk[k] = sqrt(x*x+sink*sink); 
	}
	for(k=0;k<P;k++) gold[k] = h(xk[k]/P);
	for(k=0;k<P;k++) g[k] = gold[k];
	for(k=0;k<P;k++) hk[k] = h(xk[k]);

	n = 0;
	do{
		for(j=0;j<P;j++)
		{
			s = 0;
			for(k=1;k<P;k++)
			{
				xx = xk[j]/xk[k];
				s += gold[k]*xx*xx;
			}
			printf("s = %f\n",s);
			g[j] = alpha*(hk[j]-s)+(1-alpha)*gold[j]; 
		}
		printf("g0 = %f\n",g[0]);
		norm = 0;
		norm = (g[0]-gold[0])*(g[0]-gold[0]);
		norm = sqrt(norm);
		for(j=0;j<P;j++) printf("xk[%d] = %f g[%d]=%f\n",j,xk[j],j,g[j]);
		printf("n = %d norm = %0.12f\n",n,norm);
		for(j=0;j<P;j++) gold[j] = g[j];
		if(norm<precision) break;
		n++;
	}while(n<nMax);
	if(n==nMax)
	{
		
		fprintf(stdout,"Error : divergence in g function of the QTBThermostat.c\n");
		fprintf(stdout,"P = %d x = %f\n",P,x);
		exit(1);
	}
	return gold[0];
}
*/
/*********************************************************************************/
static double interpol(int M, double* g, double dx, double x)
{
	int i = (int)(x/dx)+1;
	double dxi;
	if(i<=0) return 1.0;
	if(i>=M-1) return g[M-1];
	dxi = x-dx*(i-1);
	return (g[i]-g[i-1])/dx*dxi+g[i-1];
}
/*********************************************************************************/
/*
static double* compute_gP(int P, double xmax, int n, int ng)
{
	static double* gold  = NULL;
	static double* g  = NULL;
	static double* x  = NULL;
	static double* hk  = NULL;
	double s;
	int k;
	double alpha = 1.0/P;
	double precision = 1e-8;
	int iterMax = 10000;
	double norm;
	int iter = 0;
	double xx;
	int nn = (int)(n*sqrt(xmax*xmax+P*P)/xmax);
	int M = 2*nn*ng;
	double dx = 2*sqrt(xmax*xmax+P*P)/(M-1);
	double** xk;
	double* gp;
	int i;
	double normold;

	gold = malloc(M*sizeof(double));
	g = malloc(M*sizeof(double));
	x = malloc(M*sizeof(double));
	hk = malloc(M*sizeof(double));
	xk = malloc(M*sizeof(double*));
	for(i=0;i<M;i++) xk[i] = malloc(P*sizeof(double));

	for(i=0;i<M;i++) x[i] = dx*i;

	for(i=0;i<M;i++) 
	for(k=0;k<P;k++)
	{
		// in ceriotti notation : Hceriotti = Htukerman*P
		xk[i][k] = sqrt(P)*sin(k*M_PI/P);
		xk[i][k] *= sqrt(P);
		xk[i][k] = xk[i][k]*xk[i][k];
		xk[i][k] = sqrt(x[i]*x[i]+xk[i][k]);
	}

	gold[0] = 1.0;
	for(i=1;i<M;i++) gold[i] = h(x[i]/P);
	for(i=0;i<M;i++) g[i] = gold[i];

	hk[0] = 1.0;
	for(i=1;i<M;i++) hk[i] = h(x[i]);

	iter = 0;
	normold = 1;
	do{
		for(i=0;i<M;i++)
		{
			s = 0;
			for(k=1;k<P;k++)
			{
				xx = x[i]/xk[i][k];
				s += interpol(M, gold, dx,xk[i][k])*xx*xx;
			}
			//printf("s = %f\n",s);
			g[i] = alpha*(hk[i]-s)+(1-alpha)*gold[i]; 
		}
		printf("g1 = %f\n",g[1]);
		norm = 0;
		for(i=0;i<M/2;i++) norm += (g[i]-gold[i])*(g[i]-gold[i]);
		norm = sqrt(norm/M*2);
		printf("norm = %0.12f\n",norm);
		fflush(stdout);
		if(norm>normold) alpha/=1.1;
		else for(i=0;i<M;i++) gold[i] = g[i];
		if(norm<precision) break;
		normold = norm;
		iter++;
	}while(iter<iterMax);
	if(iter==iterMax)
	{
		
		fprintf(stdout,"Error : divergence in g function of the QTBThermostat.c\n");
		exit(1);
	}
	free(g);
	free(x);
	free(hk);
	for(i=0;i<M;i++) free(xk[i]);
	free(xk);

	gp = malloc(n*sizeof(double));
	for(i=0;i<n;i++)
		gp[i] = interpol(M, gold, dx,i*xmax/(n-1));
	free(gold);

	return gp;
}
*/
/*********************************************************************************/
static double* compute_gPTable(int P, double dx, int M)
{
	static double* gold  = NULL;
	static double* g  = NULL;
	static double* x  = NULL;
	static double* hk  = NULL;
	double s;
	int k;
	double alpha = 1.0/P;
	double precision = 1e-8;
	int iterMax = 10000;
	double norm;
	int iter = 0;
	double xx;
	double** xk;
	int i;
	double normold;

	gold = malloc(M*sizeof(double));
	g = malloc(M*sizeof(double));
	x = malloc(M*sizeof(double));
	hk = malloc(M*sizeof(double));
	xk = malloc(M*sizeof(double*));
	for(i=0;i<M;i++) xk[i] = malloc(P*sizeof(double));

	for(i=0;i<M;i++) x[i] = dx*i;

	for(i=0;i<M;i++) 
	for(k=0;k<P;k++)
	{
		/* in ceriotti notation : Hceriotti = Htukerman*P */
		/* in ceriotti notation P not sqrt(P) */
		xk[i][k] = P*sin(k*M_PI/P);
		//xk[i][k] = sqrt(P)*sin(k*M_PI/P);
		xk[i][k] = xk[i][k]*xk[i][k];
		xk[i][k] = sqrt(fabs(x[i]*x[i]+xk[i][k]));
	}

	gold[0] = 1.0;
	for(i=1;i<M;i++) gold[i] = h(x[i]/P);
	for(i=0;i<M;i++) g[i] = gold[i];

	hk[0] = 1.0;
	for(i=1;i<M;i++) hk[i] = h(x[i]);

	iter = 0;
	normold = 1;
	do{
		for(i=0;i<M;i++)
		{
			s = 0;
			for(k=1;k<P;k++)
			{
				xx = x[i]/xk[i][k];
				s += interpol(M, gold, dx,xk[i][k])*xx*xx;
			}
			//printf("s = %f\n",s);
			g[i] = alpha*(hk[i]-s)+(1-alpha)*gold[i]; 
		}
		printf("g1 = %f\n",g[1]);
		norm = 0;
		//for(i=0;i<M/2;i++) norm += (g[i]-gold[i])*(g[i]-gold[i]);
		//norm = sqrt(norm/M*2);
		for(i=0;i<M;i++) norm += (g[i]-gold[i])*(g[i]-gold[i]);
		//norm = sqrt(norm/M);
		norm = sqrt(norm);
		printf("norm = %0.12f\n",norm);
		fflush(stdout);
		if(norm>normold) alpha/=1.1;
		else for(i=0;i<M;i++) gold[i] = g[i];
		if(norm<precision) break;
		normold = norm;
		iter++;
	}while(iter<iterMax);
	if(iter==iterMax)
	{
		
		fprintf(stdout,"Error : divergence in g function of the QTBThermostat.c\n");
		exit(1);
	}
	free(g);
	free(x);
	free(hk);
	for(i=0;i<M;i++) free(xk[i]);
	free(xk);

	return gold;
}
/*********************************************************************************/
/*
static double* compute_fPTable(int P, double dx, int M)
{
	static double* gold  = NULL;
	static double* g  = NULL;
	static double* x  = NULL;
	static double* hk  = NULL;
	double s;
	int k;
	double alpha = 1.0/P;
	double precision = 1e-8;
	int iterMax = 10000;
	double norm;
	int iter = 0;
	double xx;
	double** xk;
	int i;
	double normold;

	gold = malloc(M*sizeof(double));
	g = malloc(M*sizeof(double));
	x = malloc(M*sizeof(double));
	hk = malloc(M*sizeof(double));
	xk = malloc(M*sizeof(double*));
	for(i=0;i<M;i++) xk[i] = malloc(P*sizeof(double));

	for(i=0;i<M;i++) x[i] = dx*i;

	for(i=0;i<M;i++) 
	for(k=0;k<P;k++)
	{
		// in ceriotti notation : Hceriotti = Htukerman*P 
		// in ceriotti notation P not sqrt(P) 
		xk[i][k] = P*sin(k*M_PI/P);
		xk[i][k] = xk[i][k]*xk[i][k];
		xk[i][k] = sqrt(x[i]*x[i]+xk[i][k]);
	}

	gold[0] = 1.0;
	for(i=1;i<M;i++) gold[i] = sqrt(h(x[i]/P));
	for(i=0;i<M;i++) g[i] = gold[i];

	hk[0] = 1.0;
	for(i=1;i<M;i++) hk[i] = sqrt(h(x[i]));

	iter = 0;
	normold = 1;
	do{
		for(i=0;i<M;i++)
		{
			s = 0;
			for(k=1;k<P;k++)
			{
				xx = x[i]/xk[i][k];
				s += interpol(M, gold, dx,xk[i][k])*xx;
			}
			//printf("s = %f\n",s);
			g[i] = alpha*(hk[i]-s)+(1-alpha)*gold[i]; 
		}
		printf("g1 = %f\n",g[1]);
		norm = 0;
		//for(i=0;i<M/2;i++) norm += (g[i]-gold[i])*(g[i]-gold[i]);
		//norm = sqrt(norm/M*2);
		for(i=0;i<M;i++) norm += (g[i]-gold[i])*(g[i]-gold[i]);
		norm = sqrt(norm/M);
		printf("norm = %0.12f\n",norm);
		fflush(stdout);
		if(norm>normold) alpha/=1.1;
		else for(i=0;i<M;i++) gold[i] = g[i];
		if(norm<precision) break;
		normold = norm;
		iter++;
	}while(iter<iterMax);
	if(iter==iterMax)
	{
		
		fprintf(stdout,"Error : divergence in g function of the QTBThermostat.c\n");
		exit(1);
	}
	free(g);
	free(x);
	free(hk);
	for(i=0;i<M;i++) free(xk[i]);
	free(xk);

	return gold;
}
*/
/*********************************************************************************/
/* Theta is a force */
static void removeTranslationTheta(Paths* paths)
{
	double atot[3] = {0,0,0};
	int i;
	int iAtom;
	int j;
	double totMass = 0.0;
	double** theta = paths->theta;

	if(!theta) return;
	for ( i = 0; i < paths->nBeads; i++)
	{
		for ( iAtom = 0; iAtom < paths->molecules[i]->nAtoms; iAtom++)
		{
			double mass = paths->Mprim[i][iAtom];
			totMass += mass;
			for ( j = 0; j < 3; j++)
			{
				atot[j] += theta[i][3*iAtom+j];
			}
		}
	}

	for ( j = 0; j < 3; j++) atot[j] /= totMass;


	for ( i = 0; i < paths->nBeads; i++)
		for ( iAtom = 0; iAtom < paths->molecules[i]->nAtoms; iAtom++)
			for ( j = 0; j < 3; j++)
				theta[i][3*iAtom+j] -= paths->Mprim[i][iAtom]*atot[j];
}
/*********************************************************************************/
/* Theta is a force */
static void removeRotationTheta(Paths* paths)
{
	double ptot[3] = {0,0,0};
	double cm[3] = {0,0,0};
	double L[3] = {0,0,0};
	int i;
	int iAtom;
	int j;
	int k;
	double mass = 1.0;
	double totMass = 0.0;
	double cdel[3];
	double vAng[3]={0,0,0};
	double tensor[3][3];
	double invTensor[3][3];
        double xx, xy,xz,yy,yz,zz;
	//boolean linear = FALSE;
	double** theta = paths->theta;
	double*** U = paths->U;
	double n = 0;
	/* find the center of mass coordinates  and total velocity*/
	if(!U)
	{
		paths->molecules[0]->klass->removeRotationForceCluster(paths->molecules,paths->nBeads,paths->theta);
		return;
	}

	for ( i = 0; i < paths->nBeads; i++)
	{
		for ( iAtom = 0; iAtom <  paths->molecules[i]->nAtoms; iAtom++)
		{
			mass = paths->Mprim[i][iAtom];
			totMass += mass;
			n++;
			for ( j = 0; j < 3; j++)
				cm[j] += mass*U[j][i][iAtom];
			for ( j = 0; j < 3; j++)
				ptot[j] += theta[i][3*iAtom+j];
		}
	}


	for ( j = 0; j < 3; j++) cm[j] /= totMass;

	/*   compute the angular momentum  */
	for ( i = 0; i < paths->nBeads; i++)
	{
		for ( iAtom = 0; iAtom <  paths->molecules[i]->nAtoms; iAtom++)
		{
			mass = paths->Mprim[i][iAtom];
			for ( j = 0; j < 3; j++)
			L[j] += (
				U[(j+1)%3][i][iAtom]*theta[i][3*iAtom+(j+2)%3]
			      - U[(j+2)%3][i][iAtom]*theta[i][3*iAtom+(j+1)%3]
			      );
		}
	}
	for ( j = 0; j < 3; j++)
		L[j] -= (
			cm[(j+1)%3]*ptot[(j+2)%3]
		      - cm[(j+2)%3]*ptot[(j+1)%3]
			      );

	/* calculate and invert the inertia tensor */
	for ( k = 0; k < 3; k++)
	for ( j = 0; j < 3; j++)
		tensor[k][j] = 0;
	xx = 0;
	yy = 0;
	zz = 0;
	xy = 0;
	xz = 0;
	yz = 0;
	for ( i = 0; i < paths->nBeads; i++)
	{
		for ( iAtom = 0; iAtom <  paths->molecules[i]->nAtoms; iAtom++)
		{
			mass = paths->Mprim[i][iAtom];
			for ( j = 0; j < 3; j++)
				cdel[j] = U[j][i][iAtom]-cm[j];
			xx +=  cdel[0]*cdel[0]*mass;
			xy +=  cdel[0]*cdel[1]*mass;
			xz +=  cdel[0]*cdel[2]*mass;
			yy +=  cdel[1]*cdel[1]*mass;
			yz +=  cdel[1]*cdel[2]*mass;
			zz +=  cdel[2]*cdel[2]*mass;
		}
	}
	tensor[0][0] = yy+zz;
	tensor[1][0] = -xy;
	tensor[2][0] = -xz;
	tensor[0][1] = -xy;
	tensor[1][1] = xx+zz;
	tensor[2][1] = -yz;
	tensor[0][2] = -xz;
	tensor[1][2] = -yz;
	tensor[2][2] = xx+yy;
	if(InverseTensor(tensor,invTensor))
	{
		for ( j = 0; j < 3; j++)
		{
			vAng[j] = 0;
			for ( k = 0; k < 3; k++)
				vAng[j] += invTensor[j][k]*L[k];
		}
	}
	else
	if(paths->molecules[0]->nAtoms>1)
	{
		double U0[3];
		double U1[3];
		for ( j = 0; j < 3; j++)U0[j] = U[j][0][0];
		for ( j = 0; j < 3; j++)U1[j] = U[j][0][1];
		//printf("!!!!!!!!!!!I cannot invert the rotational Tensor : linear molecule!\n");
		computeAngularVelocitiesForALinearMolecule(U0, U1, tensor, L, vAng);
		//printf("Angular velocity before rotation = %f %f %f\n",vAng[0], vAng[1], vAng[2]);
		//linear = TRUE;
	}
	/*  eliminate any rotation about the system center of mass */
	for ( i = 0; i < paths->nBeads; i++)
	{
		for ( iAtom = 0; iAtom <  paths->molecules[i]->nAtoms; iAtom++)
		{
			for ( j = 0; j < 3; j++)
				cdel[j] = U[j][i][iAtom]-cm[j];
			for ( j = 0; j < 3; j++)
				theta[i][3*iAtom+j] += 
				(cdel[(j+1)%3]*vAng[(j+2)%3]-
				cdel[(j+2)%3]*vAng[(j+1)%3])
				*paths->Mprim[i][iAtom]
				;
		}
	}
/* Check */
/*
	for ( j = 0; j < 3; j++) L[j] = 0;
	for ( i = 0; i < paths->nBeads; i++)
	{
		for ( iAtom = 0; iAtom <  paths->molecules[i]->nAtoms; iAtom++)
		{
			mass = paths->Mprim[i][iAtom];
			for ( j = 0; j < 3; j++)
			L[j] += (
				U[(j+1)%3][i][iAtom]*theta[i][3*iAtom+(j+2)%3]
			      - U[(j+2)%3][i][iAtom]*theta[i][3*iAtom+(j+1)%3]
			      );
		}
	}
	for ( j = 0; j < 3; j++)
		L[j] -= (
			cm[(j+1)%3]*ptot[(j+2)%3]
		      - cm[(j+2)%3]*ptot[(j+1)%3]
			      );
	if(linear &&paths->molecules[0]->nAtoms>1)
	{
		double U0[3];
		double U1[3];
		for ( j = 0; j < 3; j++)U0[j] = U[j][0][0];
		for ( j = 0; j < 3; j++)U1[j] = U[j][0][1];
		computeAngularVelocitiesForALinearMolecule(U0, U1, tensor, L, vAng);
		printf("Linear molecule Angular velocity = %f %f %f\n",vAng[0], vAng[1], vAng[2]);
	}
	else
	{
		for ( j = 0; j < 3; j++)
		{
			vAng[j] = 0;
			for ( k = 0; k < 3; k++)
				vAng[j] += invTensor[j][k]*L[k];
		}
		printf("Angular velocity = %f %f %f\n",vAng[0], vAng[1], vAng[2]);
	}
*/
}
/*********************************************************************************/
static void removeTranslationAndRotationTheta(Paths* paths)
{
	if(paths->nBeads<1) return;
	removeTranslationTheta(paths);
	removeRotationTheta(paths);
}
/*********************************************************************************/
/* Theta is a force */
/*
static void removeTranslationThetaOneBead(Paths* paths, int i)
{
	double atot[3] = {0,0,0};
	int iAtom;
	int j;
	double totMass = 0.0;
	double** theta = paths->theta;

	if(!theta) return;
	{
		for ( iAtom = 0; iAtom < paths->molecules[i]->nAtoms; iAtom++)
		{
			totMass += paths->Mprim[i][iAtom];
			for ( j = 0; j < 3; j++)
			{
				atot[j] += theta[i][3*iAtom+j];
			}
		}
	}

	for ( j = 0; j < 3; j++) atot[j] /= totMass;


		for ( iAtom = 0; iAtom < paths->molecules[i]->nAtoms; iAtom++)
			for ( j = 0; j < 3; j++)
				theta[i][3*iAtom+j] -= paths->Mprim[i][iAtom]*atot[j];
}
*/
/*********************************************************************************/
/* Theta is a force */
/*
static void removeRotationThetaOneBead(Paths* paths, int i)
{
	double ptot[3] = {0,0,0};
	double cm[3] = {0,0,0};
	double L[3] = {0,0,0};
	int iAtom;
	int j;
	int k;
	double mass = 1.0;
	double totMass = 0.0;
	double cdel[3];
	double vAng[3]={0,0,0};
	double tensor[3][3];
	double invTensor[3][3];
        double xx, xy,xz,yy,yz,zz;
	//boolean linear = FALSE;
	double** theta = paths->theta;
	double*** U = paths->U;
	// find the center of mass coordinates  and total velocity
	if(!U)
	{
		paths->molecules[i]->klass->removeRotationForce(paths->molecules[i],paths->theta[i]);
		return;
	}

	{
		for ( iAtom = 0; iAtom <  paths->molecules[i]->nAtoms; iAtom++)
		{
			mass = paths->Mprim[i][iAtom];
			totMass += mass;
			for ( j = 0; j < 3; j++)
				cm[j] += mass*U[j][i][iAtom];
			for ( j = 0; j < 3; j++)
				ptot[j] += theta[i][3*iAtom+j];
		}
	}


	for ( j = 0; j < 3; j++) cm[j] /= totMass;

	//   compute the angular momentum 
	{
		for ( iAtom = 0; iAtom <  paths->molecules[i]->nAtoms; iAtom++)
		{
			mass = paths->Mprim[i][iAtom];
			for ( j = 0; j < 3; j++)
			L[j] += (
				U[(j+1)%3][i][iAtom]*theta[i][3*iAtom+(j+2)%3]
			      - U[(j+2)%3][i][iAtom]*theta[i][3*iAtom+(j+1)%3]
			      );
		}
	}
	for ( j = 0; j < 3; j++)
		L[j] -= (
			cm[(j+1)%3]*ptot[(j+2)%3]
		      - cm[(j+2)%3]*ptot[(j+1)%3]
			      );

	// calculate and invert the inertia tensor
	for ( k = 0; k < 3; k++)
	for ( j = 0; j < 3; j++)
		tensor[k][j] = 0;
	xx = 0;
	yy = 0;
	zz = 0;
	xy = 0;
	xz = 0;
	yz = 0;
	{
		for ( iAtom = 0; iAtom <  paths->molecules[i]->nAtoms; iAtom++)
		{
			mass = paths->Mprim[i][iAtom];
			for ( j = 0; j < 3; j++)
				cdel[j] = U[j][i][iAtom]-cm[j];
			xx +=  cdel[0]*cdel[0]*mass;
			xy +=  cdel[0]*cdel[1]*mass;
			xz +=  cdel[0]*cdel[2]*mass;
			yy +=  cdel[1]*cdel[1]*mass;
			yz +=  cdel[1]*cdel[2]*mass;
			zz +=  cdel[2]*cdel[2]*mass;
		}
	}
	tensor[0][0] = yy+zz;
	tensor[1][0] = -xy;
	tensor[2][0] = -xz;
	tensor[0][1] = -xy;
	tensor[1][1] = xx+zz;
	tensor[2][1] = -yz;
	tensor[0][2] = -xz;
	tensor[1][2] = -yz;
	tensor[2][2] = xx+yy;
	if(InverseTensor(tensor,invTensor))
	{
		for ( j = 0; j < 3; j++)
		{
			vAng[j] = 0;
			for ( k = 0; k < 3; k++)
				vAng[j] += invTensor[j][k]*L[k];
		}
	}
	else
	if(paths->molecules[0]->nAtoms>1)
	{
		double U0[3];
		double U1[3];
		for ( j = 0; j < 3; j++)U0[j] = U[j][0][0];
		for ( j = 0; j < 3; j++)U1[j] = U[j][0][1];
		printf("!!!!!!!!!!!I cannot invert the rotational Tensor : linear molecule!\n");
		computeAngularVelocitiesForALinearMolecule(U0, U1, tensor, L, vAng);
		//printf("Angular velocityibefore rotation = %f %f %f\n",vAng[0], vAng[1], vAng[2]);
		//linear = TRUE;
	}
	//  eliminate any rotation about the system center of mass 
	{
		for ( iAtom = 0; iAtom <  paths->molecules[i]->nAtoms; iAtom++)
		{
			for ( j = 0; j < 3; j++)
				cdel[j] = U[j][i][iAtom]-cm[j];
			for ( j = 0; j < 3; j++)
				theta[i][3*iAtom+j] += 
				(cdel[(j+1)%3]*vAng[(j+2)%3]-
				cdel[(j+2)%3]*vAng[(j+1)%3])
				*paths->Mprim[i][iAtom]
				;
		}
	}
}
*/
/*********************************************************************************/
/*
static void removeTranslationAndRotationThetaOneBead(Paths* paths, int i)
{
	if(paths->nBeads<1) return;
	removeTranslationThetaOneBead(paths,i);
	removeRotationThetaOneBead(paths,i);
}
*/
