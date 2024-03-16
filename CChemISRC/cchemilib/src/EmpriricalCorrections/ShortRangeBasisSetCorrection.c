/* ShortRangeBasisSetCorrection.c */
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

/* Reference: J. Rezac, P. Hobza J. Chem. Theory Comput. 8, 141-151 (2012)*/

#ifndef OS_WIN32
#include <unistd.h>
#endif

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <ctype.h>
#ifdef ENABLE_MPI
#include <mpi.h>
#endif

#include "ShortRangeBasisSetCorrection.h"
#include "../Utils/Constants.h"
#include "../Utils/Utils.h"

static int setDefParameters(ShortRangeBasisSetCorrectionParameters* parameters, char* method);
static void printParameters(ShortRangeBasisSetCorrectionParameters* parameters, FILE* file);
/**********************************************************************/
static void printDefParameters(char* method, FILE* file)
{
	ShortRangeBasisSetCorrectionParameters pars;
	setDefParameters(&pars, method);
	printParameters(&pars, file);
}
/**********************************************************************/
static int setDefParameters(ShortRangeBasisSetCorrectionParameters* parameters, char* method)
{
	sprintf(parameters->method,"%s",method);
	if(strstr(method,"SCC-DFTB"))
	{
		parameters->nBonds=1;
        	parameters->sRBBonds = malloc(parameters->nBonds*sizeof(SRBBond));

		sprintf(parameters->sRBBonds[0].symbol1,"O");
		sprintf(parameters->sRBBonds[0].symbol2,"H");
		parameters->sRBBonds[0].A = -0.232106;
		parameters->sRBBonds[0].beta = 0.747067;
		parameters->sRBBonds[0].C = 0.0788854;
		parameters->sRBBonds[0].gamma = 0.0792195;
		parameters->sRBBonds[0].omega = 0.25;  // 1/rcut, rcut = 4 Bohr
	}
	else if(strstr(method,"BLYP"))
	{
		parameters->nBonds=1;
        	parameters->sRBBonds = malloc(parameters->nBonds*sizeof(SRBBond));

		sprintf(parameters->sRBBonds[0].symbol1,"O");
		sprintf(parameters->sRBBonds[0].symbol2,"H");
		parameters->sRBBonds[0].A = -0.232106;
		parameters->sRBBonds[0].beta = 0.747067;
		parameters->sRBBonds[0].C = 0.0788854;
		parameters->sRBBonds[0].gamma = 0.0792195;
		parameters->sRBBonds[0].omega = 0.25;  // 1/rcut, rcut = 4 Bohr
	}
	else
	{
		fprintf(stderr,"I cannot set SRB parameters, method %s is unknown\n"
	 	"The known methods are : SCC-DFTB, BLYP\n"
		,method
		);
		fprintf(stderr,"You can also give your parameters using an input file (value in atomic unit)\n"
		"Format of file \n"
		"Symbol1 Symbol2 AValue betaValue\n"
		"Symbol1 Symbol2 AValue betaValue\n"
		"Symbol1 Symbol2 AValue betaValue\n"
		);
		fprintf(stderr,"Here are the default parameters for several method\n");
		printDefParameters("SCC-DFTB", stderr);
		printDefParameters("BLYP", stderr);
		exit(1);
		return 1;
	}
	printParameters(parameters, stdout);
	return 0;
}
/**********************************************************************************************************/
/* SRB correction calculation*/
static double getSRB(ShortRangeBasisSetCorrectionParameters* parameters, Molecule* molecule, boolean addGradient)
{
	double e_corr_sum = 0;
	double rijxyz[3];
	int i,j,k,c;
	static double convg = AUTOKCAL*ANGTOBOHR;
	static double conve = AUTOKCAL;
	static double twoOverPI05 = 2.0/sqrt(M_PI);

	for(i=0;i<molecule->nAtoms-1;i++)
	{
		for (j = i+1; j < molecule->nAtoms; j++)
		{ 
			double rij = getDistance(&molecule->atoms[i],&molecule->atoms[j]);
			rij = rij*ANGTOBOHR;
			//fprintf(stdout,"DEBUG ij = %d, %d rij=%lf\n", i,j,rij);
			for(k=0;k<parameters->nBonds;k++)
			if(
				( !strcmp(molecule->atoms[i].prop.symbol,parameters->sRBBonds[k].symbol1) && !strcmp(molecule->atoms[j].prop.symbol,parameters->sRBBonds[k].symbol2)) ||
				( !strcmp(molecule->atoms[j].prop.symbol,parameters->sRBBonds[k].symbol1) && !strcmp(molecule->atoms[i].prop.symbol,parameters->sRBBonds[k].symbol2)) 
			)
			{
			// V= A*exp(beta*rij)*erfc(w rij)+C*exp(gamma*rij)*erf(w rij)
			// erfc(x) = 1 - erf(x)
			// derive de erf(x)= 2/sqrt(M_PI)*exp(-x*x)
			// deriv xi [A*exp(beta*rij)*erfc(w rij)] = A*exp(beta*rij) [ beta erfc(w rij)-2 w/sqrt(M_PI)*exp(-(w*rij)^2) ] (xi-xj)/rij
			// deriv xi [C*exp(gamma*rij)*erf(w rij)] = C*exp(gamma*rij) [ gamma erf(w rij)+2 w/sqrt(M_PI)*exp(-(w*rij)^2) ] (xi-xj)/rij

				double A=parameters->sRBBonds[k].A;
				double beta=parameters->sRBBonds[k].beta;
				double gamma=parameters->sRBBonds[k].gamma;
				double C=parameters->sRBBonds[k].C;
				double w=parameters->sRBBonds[k].omega;
				//fprintf(stdout,"DEBUG ij = %d, %d A=%lf beta=%lf, C=%lf gamma=%lf omega=%lf\n", i,j,A,beta,C,gamma,omega);
                                //w = 0;
				//C = 0;
				double wr = rij*w;
				double erfwr = erf(wr);
				double erfcwr = 1-erfwr;
				double Abeta= A*exp(beta*rij);
				double Cgamma= C*exp(gamma*rij);
				e_corr_sum += Abeta*erfcwr+Cgamma*erfwr;
				if (addGradient) 
				{
					double derfwr = twoOverPI05*w*exp(-wr*wr);
					double term1 = Abeta*(beta*erfcwr-derfwr);
					double term2 = Cgamma*(gamma*erfwr+derfwr);
					double term = (term1+term2)/rij;
					term *= convg;
					for(c=0;c<3;c++) rijxyz[c] =  term*(molecule->atoms[i].coordinates[c] - molecule->atoms[j].coordinates[c])*ANGTOBOHR;
					for(c=0;c<3;c++) molecule->atoms[i].gradient[c] += rijxyz[c];
					for(c=0;c<3;c++) molecule->atoms[j].gradient[c] -= rijxyz[c];
                		}
			}
		}
	}
	e_corr_sum *= conve;
	return e_corr_sum;
}
/**********************************************************************/
static int readOneLine(char* buffer, char* symbol1, char* symbol2, double* pA, double* pBeta, double* pOmega, double* pC, double* pGamma)
{
	int k=sscanf(buffer,"%s %s %le %le %le %le %le", symbol1, symbol2, pA, pBeta, pOmega, pC, pGamma);
	//fprintf(stdout,"DEBUG k=%d\n",k);
	if(k>=4)
	{
		int i;
		for(i=1;i<strlen(symbol1);i++) symbol1[i]=tolower(symbol1[i]);
		for(i=1;i<strlen(symbol2);i++) symbol2[i]=tolower(symbol2[i]);
	}

	return k;
}
/**********************************************************************/
static boolean getOneString(char* buffer, char* name, char*pval)
{
        char* st = strstr(buffer,name);
        if(st)
        {
                char* beg = strstr(buffer,"=");
                int k = 0;
                if(beg) k=sscanf(beg+1,"%s",pval);
                else k=sscanf(st+strlen(st),"%s",pval);
                if(k==1) return TRUE;
        }
        return FALSE;
}
/**********************************************************************/
static void printParameters(ShortRangeBasisSetCorrectionParameters* parameters, FILE* file)
{
	int i;
	fprintf(file,"---------------------------------------------------\n");
	fprintf(file,"SRB Correction Parameters:\n");
	fprintf(file,"Method=%s\n",parameters->method);
	for(i=0;i<parameters->nBonds;i++) 
		fprintf(file,"%s %s %le %le %le %le %le\n",
			parameters->sRBBonds[i].symbol1, parameters->sRBBonds[i].symbol2,
			parameters->sRBBonds[i].A,parameters->sRBBonds[i].beta,
			parameters->sRBBonds[i].omega,
			parameters->sRBBonds[i].C,parameters->sRBBonds[i].gamma
			);
	fprintf(file,"End\n");
}
/**********************************************************************/
int readShortRangeBasisSetCorrectionParameters(ShortRangeBasisSetCorrectionParameters* parameters, char* fileName)
{
	char buffer[BSIZE];
	FILE* file = fopen(fileName,"r");
	int nBonds=0;
	char symbol1[10];
	char symbol2[10];
	double A;
	double beta;
	double C;
	double gamma;
	double omega;
	int k;

	sprintf(parameters->method,"%s","Generic");
	while(file && !feof(file))
	{
		if(!fgets(buffer,BSIZE,file))break;
		uppercase(buffer);
		if(strstr(buffer,"END")) break;
		if(!getOneString(buffer, "METHOD", parameters->method))
			if(readOneLine(buffer, symbol1, symbol2, &A, &beta, &omega, &C, &gamma)>=4) nBonds++;
	}
	if(!file)
	{
		fprintf(stderr,"Sorry I cannot opent %s file\n",fileName);
		exit(1);
	}
	if(nBonds<1)
	{
		fprintf(stderr,"Sorry I cannot read the SRB parameters from %s file\n",fileName);
		exit(1);
	}
	rewind(file);
	parameters->nBonds=nBonds;
        parameters->sRBBonds = malloc(parameters->nBonds*sizeof(SRBBond));
	int i=0;
	while(file && !feof(file))
	{
		if(!fgets(buffer,BSIZE,file))break;
		uppercase(buffer);
		if(strstr(buffer,"END")) break;
		if(strstr(buffer, "METHOD")) continue;
		k =readOneLine(buffer, symbol1, symbol2, &A, &beta, &omega, &C, &gamma);
		if(k>=4)
		{
			sprintf(parameters->sRBBonds[i].symbol1,"%s", symbol1);
			sprintf(parameters->sRBBonds[i].symbol2,"%s", symbol2);
			parameters->sRBBonds[i].A = A; 
			parameters->sRBBonds[i].beta = beta;
			parameters->sRBBonds[i].omega = 0;
			parameters->sRBBonds[i].C = 0; 
			parameters->sRBBonds[i].gamma = 0.0;
			if(k>=5) parameters->sRBBonds[i].omega = omega; 
			if(k>=6)
			{
				parameters->sRBBonds[i].C = C; 
				parameters->sRBBonds[i].gamma = gamma;
			}
			i++;
		}
	}

	printParameters(parameters, stdout);
	return 0;
}
/**********************************************************************/
int setShortRangeBasisSetCorrectionParameters(ShortRangeBasisSetCorrectionParameters* parameters, char* fileName, char* method)
{
	if(fileName) return readShortRangeBasisSetCorrectionParameters(parameters, fileName);
	return setDefParameters(parameters, method);
	return 1;
}
/****************************************************************************************************/
double getSRBCorrection(Molecule* molecule, ShortRangeBasisSetCorrectionParameters* parameters, boolean addGradient)
{
	double e;
	//fprintf(stdout," DEBUG je suis dans getSRBCorrection\n");
	//if(!parameters) { fprintf(stdout," DEBUG SRB parameters= NULL\n"); return 0; }
	if(!parameters) return 0;
	e = getSRB(parameters, molecule, addGradient);
	/* printf("SRB Energy=%f\n",e);*/
	return e;
}
