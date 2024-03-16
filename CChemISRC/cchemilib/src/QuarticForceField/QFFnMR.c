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

/* QFFnMR.c */
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <ctype.h>

#include "../Utils/Utils.h"
#include "../Utils/Types.h"
#include "../Utils/Constants.h"
#include "../Molecule/Molecule.h"

/************************************************************************************************************/
typedef struct _NMRQFFDipole
{
	double values[3];
}NMRQFFDipole;

typedef struct _NMRQFF
{
	int numberOfFrequencies;
	double* gradients;
	double* frequencies;
	int numberOfFirstDipolesInput;
	int numberOfEnergies;
	int numberOfDipoles;
	double** firstDipolesInput;
	double* calculatedFrequencies;
	double* mass;
	double* delta;//[i:0..Nfreq]
	double V0;
	double **VI; // [i:0..Nfreq][0..5]
	double ***VIJ; // [i:0..Nfreq][j:0..Nfreq][0..3]
	double ***VI3J; // [i:0..Nfreq][j:0..Nfreq][0..3]
	double ****VIJK; // [i:0..Nfreq][j:0..Nfreq][[j:0..Nfreq][0..7]
	double *****VIJKL; // [i:0..Nfreq][j:0..Nfreq][j:0..Nfreq][j:0..Nfreq][0..16]

	double** cubicEnergies;//[i:0..Nfreq][j:0..Nfreq] for tiii, tiij
	double*** cubicEnergiesIJK;//[i:0..Nfreq][j:0..Nfreq][j:0..Nfreq] for tijk, i # j # k
	double** quarticEnergiesIIJJ;//[i:0..Nfreq][j:0..i-1] for uiijj
	double** quarticEnergiesIIIJ;//[i:0..Nfreq][j:0..Nfreq] for uiiij and uiiii
	double*** quarticEnergiesIIJK;//[i:0..Nfreq][j:0..Nfreq][k:0..Nfreq] for uiijk i # j # k
	double**** quarticEnergiesIJKL;//[i:0..Nfreq][j:0..Nfreq][k:0..Nfreq][l:0..Nfreq] for uijkl i # j # k # l

	double *dipole0; // [0..2] // 0..2 : x,y,z
	double ***dipolesI; // mu[i:0..Nfreq][0..5][0..2] // 0..2 : x,y,z
	double ****dipolesIJ; // mu[i:0..Nfreq][j:0..Nfreq][0..3][0..2] // 0..2 : x,y,z
	double ****dipolesI3J; // mu[i:0..Nfreq][j:0..Nfreq][0..3][0..2] // 0..2 : x,y,z
	double *****dipolesIJK; // mu[i:0..Nfreq][j:0..Nfreq][k:0..Nfreq][0..8][0..2] // 0..2 : x,y,z
	double** firstDipoles;//[i:0..Nfreq][0..2] 
	double*** secondDipoles;//[i:0..Nfreq][i:0..Nfreq][0..2] 
	double*** cubicDipoles;//[i:0..Nfreq][i:0..Nfreq][0..2] // diii, diij
	double**** cubicDipolesIJK;//[i:0..Nfreq][i:0..Nfreq][i:0..Nfreq][0..2] // dijk, i # j # k
	double*** quarticDipolesIIJJ;//[i:0..Nfreq][j:0..i-1][0..2] for diijj
	double*** quarticDipolesIIIJ;//[i:0..Nfreq][j:0..Nfreq][0..2] for diiij and diiii
	double**** quarticDipolesIIJK;//[i:0..Nfreq][j:0..Nfreq][k:0..Nfreq][0..2] for diijk i # j # j

	Molecule mol;
}NMRQFF;
static void computeQFFDerivativesEnergies(NMRQFF* qffConstants, int order, double* energies);
static void computeQFFDerivativesDipoles(NMRQFF* qffConstants, int order, double* dipoles[]);
static char* saveNMRQFFAppend(NMRQFF* qffConstants, char* inputFileName, Molecule* mol);
/************************************************************************************************************/
static void initnMRQFF0(NMRQFF* qffConstants, int nf)
{
	qffConstants->numberOfFrequencies = nf;
	qffConstants->numberOfFirstDipolesInput = 0;
	qffConstants->numberOfEnergies = 0;
	qffConstants->numberOfDipoles = 0;

	//printf("Begin initNMRQFF0 nf = %d\n",nf);

	qffConstants->frequencies = newVectorDouble(qffConstants->numberOfFrequencies);
        initVectorDouble(qffConstants->frequencies, qffConstants->numberOfFrequencies, 0.0);
	qffConstants->gradients = newVectorDouble(qffConstants->numberOfFrequencies);
        initVectorDouble(qffConstants->gradients, qffConstants->numberOfFrequencies, 0.0);
	//printf("End initgrad\n");
	qffConstants->calculatedFrequencies = newVectorDouble(qffConstants->numberOfFrequencies);
        initVectorDouble(qffConstants->calculatedFrequencies, qffConstants->numberOfFrequencies, 0.0);
	qffConstants->mass = newVectorDouble(qffConstants->numberOfFrequencies);
        initVectorDouble(qffConstants->mass, qffConstants->numberOfFrequencies, 0.0);
	qffConstants->delta = newVectorDouble(qffConstants->numberOfFrequencies);
        initVectorDouble(qffConstants->delta, qffConstants->numberOfFrequencies, 0.0);
	//printf("End initNMRQFF0\n");

	qffConstants->mol.nAtoms = 0;
	qffConstants->mol.atoms = NULL;

}
/************************************************************************************************************/
static void initnMRQFFEnergies(NMRQFF* qffConstants, int order)
{
	qffConstants->V0 = 0;
	qffConstants->VI = newMatrixDouble(qffConstants->numberOfFrequencies,6);
        initMatrixDouble(qffConstants->VI, qffConstants->numberOfFrequencies,6, 0.0);
	qffConstants->cubicEnergies = newMatrixDouble(qffConstants->numberOfFrequencies,qffConstants->numberOfFrequencies);
        initMatrixDouble(qffConstants->cubicEnergies, qffConstants->numberOfFrequencies,qffConstants->numberOfFrequencies, 0.0);

	qffConstants->cubicEnergiesIJK = newCubeDouble(qffConstants->numberOfFrequencies,qffConstants->numberOfFrequencies,qffConstants->numberOfFrequencies);
        initCubeDouble(qffConstants->cubicEnergiesIJK, qffConstants->numberOfFrequencies,qffConstants->numberOfFrequencies,qffConstants->numberOfFrequencies, 0.0);
	//printf("End initCubeDouble\n");

	qffConstants->quarticEnergiesIIJJ = newMatrixDouble(qffConstants->numberOfFrequencies,qffConstants->numberOfFrequencies);
        initMatrixDouble(qffConstants->quarticEnergiesIIJJ, qffConstants->numberOfFrequencies,qffConstants->numberOfFrequencies, 0.0);
	qffConstants->quarticEnergiesIIIJ = newMatrixDouble(qffConstants->numberOfFrequencies,qffConstants->numberOfFrequencies);
        initMatrixDouble(qffConstants->quarticEnergiesIIIJ, qffConstants->numberOfFrequencies,qffConstants->numberOfFrequencies, 0.0);
	qffConstants->quarticEnergiesIIJK = newCubeDouble(qffConstants->numberOfFrequencies,qffConstants->numberOfFrequencies,qffConstants->numberOfFrequencies);
        initCubeDouble(qffConstants->quarticEnergiesIIJK, qffConstants->numberOfFrequencies,qffConstants->numberOfFrequencies,qffConstants->numberOfFrequencies, 0.0);

        qffConstants->quarticEnergiesIJKL = newQuarticDouble(qffConstants->numberOfFrequencies,qffConstants->numberOfFrequencies,qffConstants->numberOfFrequencies,qffConstants->numberOfFrequencies);
        initQuarticDouble(qffConstants->quarticEnergiesIJKL, qffConstants->numberOfFrequencies,qffConstants->numberOfFrequencies,qffConstants->numberOfFrequencies,qffConstants->numberOfFrequencies,0.0);
	//printf("End quarticEnergiesIIJK\n");

        qffConstants->VIJ = newCubeDouble(qffConstants->numberOfFrequencies,qffConstants->numberOfFrequencies,4);
        initCubeDouble(qffConstants->VIJ, qffConstants->numberOfFrequencies,qffConstants->numberOfFrequencies,4,0.0);
        qffConstants->VI3J = newCubeDouble(qffConstants->numberOfFrequencies,qffConstants->numberOfFrequencies,4);
        initCubeDouble(qffConstants->VI3J, qffConstants->numberOfFrequencies,qffConstants->numberOfFrequencies,4,0.0);

	//printf("Begin VIJK\n");
        qffConstants->VIJK =  NULL;
	if(order>=3)
	{
		// save VIJK in file
		/*
        	qffConstants->VIJK = newQuarticDouble(qffConstants->numberOfFrequencies,qffConstants->numberOfFrequencies,qffConstants->numberOfFrequencies,8);
        	initQuarticDouble(qffConstants->VIJK, qffConstants->numberOfFrequencies,qffConstants->numberOfFrequencies,qffConstants->numberOfFrequencies,8,0.0);
		*/
	}
	//printf("End VIJK\n");

	//printf("Begin VIJKL\n");
        qffConstants->VIJKL = NULL;
	if(order>=4)
	{
        	qffConstants->VIJKL = newQuinticDouble(qffConstants->numberOfFrequencies,qffConstants->numberOfFrequencies,qffConstants->numberOfFrequencies,qffConstants->numberOfFrequencies,16);
        	initQuinticDouble(qffConstants->VIJKL, qffConstants->numberOfFrequencies,qffConstants->numberOfFrequencies,qffConstants->numberOfFrequencies,qffConstants->numberOfFrequencies,16,0.0);
	}
	//printf("End VIJKL\n");
}
/************************************************************************************************************/
static void initnMRQFFDipoles(NMRQFF* qffConstants, int order)
{
	qffConstants->quarticDipolesIIJJ = newCubeDouble(qffConstants->numberOfFrequencies,qffConstants->numberOfFrequencies,3);
        initCubeDouble(qffConstants->quarticDipolesIIJJ, qffConstants->numberOfFrequencies,qffConstants->numberOfFrequencies,3, 0.0);
	qffConstants->quarticDipolesIIIJ = newCubeDouble(qffConstants->numberOfFrequencies,qffConstants->numberOfFrequencies,3);
        initCubeDouble(qffConstants->quarticDipolesIIIJ, qffConstants->numberOfFrequencies,qffConstants->numberOfFrequencies,3, 0.0);
        qffConstants->quarticDipolesIIJK = newQuarticDouble(qffConstants->numberOfFrequencies,qffConstants->numberOfFrequencies,qffConstants->numberOfFrequencies,3);
        initQuarticDouble(qffConstants->quarticDipolesIIJK, qffConstants->numberOfFrequencies,qffConstants->numberOfFrequencies,qffConstants->numberOfFrequencies,3,0.0);
	//printf("End quarticDipolesIIJK\n");

	qffConstants->firstDipoles = newMatrixDouble(qffConstants->numberOfFrequencies,3);
        initMatrixDouble(qffConstants->firstDipoles, qffConstants->numberOfFrequencies,3, 0.0);

	qffConstants->firstDipolesInput = newMatrixDouble(qffConstants->numberOfFrequencies,3);
        initMatrixDouble(qffConstants->firstDipolesInput, qffConstants->numberOfFrequencies,3, 0.0);

        qffConstants->dipole0 = newVectorDouble(3);
        initVectorDouble(qffConstants->dipole0,3,0.0);

        qffConstants->dipolesI = newCubeDouble(qffConstants->numberOfFrequencies,6,3);
        initCubeDouble(qffConstants->dipolesI, qffConstants->numberOfFrequencies,6,3,0.0);

        qffConstants->secondDipoles = newCubeDouble(qffConstants->numberOfFrequencies,qffConstants->numberOfFrequencies,3);
        initCubeDouble(qffConstants->secondDipoles, qffConstants->numberOfFrequencies,qffConstants->numberOfFrequencies,3,0.0);

        qffConstants->cubicDipoles = newCubeDouble(qffConstants->numberOfFrequencies,qffConstants->numberOfFrequencies,3);
        initCubeDouble(qffConstants->cubicDipoles, qffConstants->numberOfFrequencies,qffConstants->numberOfFrequencies,3,0.0);

        qffConstants->cubicDipolesIJK = newQuarticDouble(qffConstants->numberOfFrequencies,qffConstants->numberOfFrequencies,qffConstants->numberOfFrequencies,3);
        initQuarticDouble(qffConstants->cubicDipolesIJK, qffConstants->numberOfFrequencies,qffConstants->numberOfFrequencies,qffConstants->numberOfFrequencies,3,0.0);

        qffConstants->dipolesIJ = newQuarticDouble(qffConstants->numberOfFrequencies,qffConstants->numberOfFrequencies,4,3);
        initQuarticDouble(qffConstants->dipolesIJ, qffConstants->numberOfFrequencies,qffConstants->numberOfFrequencies,4,3,0.0);

        qffConstants->dipolesI3J = newQuarticDouble(qffConstants->numberOfFrequencies,qffConstants->numberOfFrequencies,4,3);
        initQuarticDouble(qffConstants->dipolesI3J, qffConstants->numberOfFrequencies,qffConstants->numberOfFrequencies,4,3,0.0);

	//printf("Begin dipolesIJK\n");
        qffConstants->dipolesIJK = NULL;
	if(order>=3)
	{
		//save in file
		/*
        	qffConstants->dipolesIJK = newQuinticDouble(qffConstants->numberOfFrequencies,qffConstants->numberOfFrequencies,qffConstants->numberOfFrequencies,8,3);
        	initQuinticDouble(qffConstants->dipolesIJK, qffConstants->numberOfFrequencies, qffConstants->numberOfFrequencies, qffConstants->numberOfFrequencies,8,3,0.0);
		*/
	}
	//printf("End dipolesIJK\n");
}
/************************************************************************************************************/
/*
static void freeNMRQFF(NMRQFF* qffConstants)
{
	freeVectorDouble(&qffConstants->frequencies);
	freeVectorDouble(&qffConstants->gradients);
	freeVectorDouble(&qffConstants->calculatedFrequencies);
	freeVectorDouble(&qffConstants->mass);
	freeVectorDouble(&qffConstants->delta);
	freeMatrixDouble(&qffConstants->VI, qffConstants->numberOfFrequencies);
	freeMatrixDouble(&qffConstants->cubicEnergies, qffConstants->numberOfFrequencies);
	freeMatrixDouble(&qffConstants->quarticEnergiesIIJJ, qffConstants->numberOfFrequencies);
	freeMatrixDouble(&qffConstants->quarticEnergiesIIIJ, qffConstants->numberOfFrequencies);
	freeCubeDouble(&qffConstants->quarticDipolesIIJJ, qffConstants->numberOfFrequencies,qffConstants->numberOfFrequencies);
	freeCubeDouble(&qffConstants->quarticDipolesIIIJ, qffConstants->numberOfFrequencies,qffConstants->numberOfFrequencies);
	freeMatrixDouble(&qffConstants->firstDipoles, qffConstants->numberOfFrequencies);
	freeMatrixDouble(&qffConstants->firstDipolesInput, qffConstants->numberOfFrequencies);
	freeCubeDouble(&qffConstants->VIJ, qffConstants->numberOfFrequencies,qffConstants->numberOfFrequencies);
	freeCubeDouble(&qffConstants->VI3J, qffConstants->numberOfFrequencies,qffConstants->numberOfFrequencies);

	freeVectorDouble(&qffConstants->dipole0);
	freeCubeDouble(&qffConstants->dipolesI, qffConstants->numberOfFrequencies,6);
	freeCubeDouble(&qffConstants->secondDipoles, qffConstants->numberOfFrequencies,qffConstants->numberOfFrequencies);
	freeCubeDouble(&qffConstants->cubicDipoles, qffConstants->numberOfFrequencies,qffConstants->numberOfFrequencies);

	freeQuarticDouble(&qffConstants->dipolesIJ, qffConstants->numberOfFrequencies,qffConstants->numberOfFrequencies,4);
	freeQuarticDouble(&qffConstants->dipolesI3J, qffConstants->numberOfFrequencies,qffConstants->numberOfFrequencies,4);
}
*/
/********************************************************************************/
static void saveNMRQFF(NMRQFF* qffConstants, FILE* file)
{
	char tmp[BSIZE];
	int i,j;
	int nf = qffConstants->numberOfFrequencies;
	char txyz[]={'X','Y','Z'};
	int xyz;

	sprintf(tmp, "#====================================================================================\n");
	fprintf(file,"%s", tmp);   
	sprintf(tmp,"%s",
		"# nMR-QFF constants\n"
		"# See Yagi et al. J. Chem. Phys. 121, 1383 (2004)\n"
		);
	fprintf(file,"%s",tmp);   
	sprintf(tmp, "#====================================================================================\n");
	fprintf(file,"%s",tmp);   

	fprintf(file,"%s","\n");
	fprintf(file,"%s","VPT2Model=GVPT2\n");   
	fprintf(file,"%s","# VPT2Model=DCPT2\n");   
	fprintf(file,"%s","# VPT2Model=HDCPT2\n");   
	fprintf(file,"%s","# alphaHDCPT2=1.0\n");   
	fprintf(file,"%s","# betaHDCPT2=5e5\n");   
	fprintf(file,"%s","\n");
	fprintf(file,"%s","PropModel=GVPT2\n");
	fprintf(file,"%s","# PropModel=HDCPT2\n");
	fprintf(file,"%s","# PropModel=DCPT2\n");
	fprintf(file,"%s","# alphaPropHDCPT2=1.0\n");
	fprintf(file,"%s","# betaPropHDCPT2=5e5\n");
	fprintf(file,"%s","# alphaPropHDCPT2=1.0\n");
	fprintf(file,"%s","# betaPropHDCPT2=5e5\n");
	fprintf(file,"%s","maxFrequencyDifferenceFermi=200\n");
	fprintf(file,"%s","MartinCutOff1=1.0\n");
	fprintf(file,"%s","MartinCutOff2=1.0\n");
	fprintf(file,"%s","# ZCutOff=0.08\n");
	fprintf(file,"%s","\n");
	sprintf(tmp, "#====================================================================================\n");
	fprintf(file,"%s",tmp);   

	fprintf(file,"%s","\n");   
	sprintf(tmp,"nFrequencies=%d\n",nf);
	fprintf(file,"%s",tmp);   
	sprintf(tmp,"nDim=%d\n",3);
	fprintf(file,"%s",tmp);   
	fprintf(file,"%s","\n");   

	sprintf(tmp,"#i Freq(cm-1)  Calc.Freq   dQ(Bohr)  Mass(amu)\tGradient[ H amu^(-1/2) Bohr^(-1)]\n");
	fprintf(file,"%s",tmp);   
	fprintf(file,"%s","Hessian\n");   
	for(i=0;i<nf;i++) 
	{
		sprintf(tmp,"%d %d %0.14f %0.14f %0.14f %0.14f\t%0.14f\n",i+1, i+1, qffConstants->frequencies[i], 
				qffConstants->calculatedFrequencies[i], 
				qffConstants->delta[i], qffConstants->mass[i], qffConstants->gradients[i]);
		fprintf(file,"%s",tmp);   
	}
	fprintf(file,"%s","END\n\n");


	sprintf(tmp,"# i\tj\tk\tReduced values [cm-1]\n");
	fprintf(file,"%s",tmp);   
	fprintf(file,"%s","Cubic\n");   
	for(i=0;i<nf;i++)
        {
                for(j=0;j<nf;j++)
                {
			if(fabs(qffConstants->cubicEnergies[i][j])<1e-12) continue;
			sprintf(tmp,"%d\t%d\t%d\t%14.6f\n",i+1, i+1, j+1, qffConstants->cubicEnergies[i][j]);
			fprintf(file,"%s",tmp);   
                }
        }
	fprintf(file,"%s","END\n\n");


	sprintf(tmp,"# i\tj\tk\tl\tReduced values [cm-1]\n");
	fprintf(file,"%s",tmp);   
	fprintf(file,"%s","Quartic\n");   
	for(i=0;i<nf;i++)
        {
		for(j=0;j<=i;j++)
		{
			if(fabs(qffConstants->quarticEnergiesIIJJ[i][j])<1e-12) continue;
			sprintf(tmp,"%d\t%d\t%d\t%d\t%14.6f\n",i+1, i+1, j+1, j+1, qffConstants->quarticEnergiesIIJJ[i][j]);
			fprintf(file,"%s",tmp);   
		}
		for(j=0;j<nf;j++)
		{
			if(j==i) continue;
			if(fabs(qffConstants->quarticEnergiesIIIJ[i][j])<1e-12) continue;
			sprintf(tmp,"%d\t%d\t%d\t%d\t%14.6f\n",i+1, i+1, i+1, j+1, qffConstants->quarticEnergiesIIIJ[i][j]);
			fprintf(file,"%s",tmp);   
		}
        }
	fprintf(file,"%s","END\n\n");

	if(qffConstants->numberOfFirstDipolesInput==0)
	{
		sprintf(tmp,"#xyz\ti\tValues[au cm^1/2]\n");
		fprintf(file,"%s",tmp);   
		fprintf(file,"%s","First derivatives\n");
		for(i=0;i<nf;i++)
		for(xyz=0;xyz<3;xyz++)
		{
			if(fabs(qffConstants->firstDipoles[i][xyz])<1e-12) continue;
			sprintf(tmp,"%c\t%d\t%14.6f\n",txyz[xyz], i+1,qffConstants->firstDipoles[i][xyz]);
			fprintf(file,"%s",tmp);   
		}
		fprintf(file,"%s","END\n\n");
	}
	else
	{
		sprintf(tmp,"#xyz\ti\tInput values[au cm^1/2]\tCalculated values[au cm^1/2]\n");
		fprintf(file,"%s",tmp);   
		fprintf(file,"%s","First derivatives\n");
		for(i=0;i<nf;i++)
		for(xyz=0;xyz<3;xyz++)
		{
			if(fabs(qffConstants->firstDipolesInput[i][xyz])<1e-12) continue;
			sprintf(tmp,"%c\t%d\t%14.6f\t\t%14.6f\n",txyz[xyz], i+1,qffConstants->firstDipolesInput[i][xyz], qffConstants->firstDipoles[i][xyz]);
			fprintf(file,"%s",tmp);   
		}
		fprintf(file,"%s","END\n\n");
	}

	sprintf(tmp,"#xyz\ti\tj\tValues[au cm]\n");
	fprintf(file,"%s",tmp);   
	fprintf(file,"%s","Second derivatives\n");
	for(i=0;i<nf;i++)
	for(j=0;j<nf;j++)
	for(xyz=0;xyz<3;xyz++)
	{
		if(fabs(qffConstants->secondDipoles[i][j][xyz])<1e-12) continue;
		sprintf(tmp,"%c\t%d\t%d\t%14.6f\n",txyz[xyz], i+1,j+1,qffConstants->secondDipoles[i][j][xyz]);
		fprintf(file,"%s",tmp);   
	}
	fprintf(file,"%s","END\n\n");

	sprintf(tmp,"#xyz\ti\tj\tk\tValues[au cm^3/2]\n");
	fprintf(file,"%s",tmp);   
	fprintf(file,"%s","Cubic derivatives\n");
	for(i=0;i<nf;i++)
	for(j=0;j<nf;j++)
	for(xyz=0;xyz<3;xyz++)
	{
		if(fabs(qffConstants->cubicDipoles[i][j][xyz])<1e-12) continue;
		sprintf(tmp,"%c\t%d\t%d\t%d\t%14.6f\n",txyz[xyz], i+1,i+1,j+1,qffConstants->cubicDipoles[i][j][xyz]);
		fprintf(file,"%s",tmp);   
		if(i!=j)
		{
			sprintf(tmp,"%c\t%d\t%d\t%d\t%14.6f\n",txyz[xyz], i+1,j+1,i+1,qffConstants->cubicDipoles[i][j][xyz]);
			fprintf(file,"%s",tmp);   
			sprintf(tmp,"%c\t%d\t%d\t%d\t%14.6f\n",txyz[xyz], j+1,i+1,i+1,qffConstants->cubicDipoles[i][j][xyz]);
			fprintf(file,"%s",tmp);   
		}
	}
	fprintf(file,"%s","END\n\n");

}
/********************************************************************************/
static char* printnMRQFF(NMRQFF* qffConstants, char* inputFileName)
{
	char* fileNameOut = strdup_printf("%sQFF.txt",getSuffixNameFile(inputFileName));
	FILE* file = fopen(fileNameOut,"w");
	fprintf(stdout,"QFF parameters saved in %s file\n", fileNameOut);
	saveNMRQFF(qffConstants,file);
	fclose(file);
	return (fileNameOut);
}
/************************************************************************************************************/
static boolean readFrequenciesInitNMRQFF(FILE* file, NMRQFF* qffConstants)
{
	char t[BSIZE];
 	int nf = 0;
	boolean Ok = TRUE;
	int nn=1;
	double dum;
	int i;
	int order = 2;
	

	if(!goToStr(file, "Frequencies"))
	{
		fprintf(stderr,"I cannot read the harmonic frequencies\nChech your input file\n");
		return FALSE;
	}
	while(!feof(file))
	{
		if(!fgets(t,BSIZE,file))break;
		nn = sscanf(t,"%lf",&dum);
		if(nn<1) break;
		nf++;
	}
	if(nf==0)
	{
		fprintf(stderr,"I cannot read the harmonic frequencies\nChech your input file\n");
		return FALSE;
	}
	readOneInt(file,"QFFnModes",&order);
	initnMRQFF0(qffConstants,nf);
	rewind(file);
	goToStr(file, "Frequencies");
	for(i=0;i<nf;i++)
	{
		if(!fgets(t,BSIZE,file))break;
		nn = sscanf(t,"%lf",&qffConstants->frequencies[i]);
		if(nn<1) break;
	}
	qffConstants->numberOfFirstDipolesInput = 0;
	if(nn==1)
	{
		rewind(file);
		if(goToStr(file, "First derivatives"))
		for(i=0;i<nf && nn==1 ;i++)
		{
			int xyz;
			for(xyz=0;xyz<3 && nn==1 ;xyz++)
				nn = fscanf(file,"%lf",&qffConstants->firstDipolesInput[i][xyz]);
			qffConstants->numberOfFirstDipolesInput+= nn;
		}
		if(qffConstants->numberOfFirstDipolesInput != nf) qffConstants->numberOfFirstDipolesInput=0;
	}
	if(nn!=1) Ok = FALSE;
	return Ok;
}
/************************************************************************************************************/
static boolean readVectorRealQFF(FILE* file, char* tag, int n, double* values)
{
        char* TAG = NULL;
        int i=0;
	int nn = 0;
        if(!tag) return FALSE;
        if(!values) return FALSE;

        TAG = strdup(tag);
        uppercase(TAG);
        rewind(file);
	if(!goToStr(file, TAG)) 
	{
		fprintf(stderr,"I cannot find %s in our file\n",tag);
		if(TAG) free(TAG);
		return FALSE;
	}
	for(i=0;i<n;i++)
        {
		nn = fscanf(file,"%lf",&values[i]);
                if(nn!=1) break;
        }
	if(i!=n)
	{
		fprintf(stderr,"I cannot read %s\nCheck  the number of values\n",tag);
		return FALSE;
	}
	return TRUE;
}
/************************************************************************************************************/
static int getOrdre(int f, int nGeoms)
{
	int order = -1;
	if(nGeoms == 1 + 6*f) order = 1;
	if(nGeoms == 1 + 6*f+6*f*(f-1) ) order = 2;
	if(nGeoms == 1 + 6*f+6*f*(f-1) + 8*f*(f-1)*(f-2)/6 ) order = 3;
	if(nGeoms == 1 + 6*f+6*f*(f-1) + 8*f*(f-1)*(f-2)/6 + 16*f*(f-1)*(f-2)*(f-3)/24) order = 4;
	if(order<1)
	{
		printf("nModes = %d\n",order);
               	fprintf(stderr,"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n");
               	fprintf(stderr,"Error, the number of files does not correspond to any known nMode (1,2,3,or4)\n");
               	fprintf(stderr,"Check the number of *QFF_*.gab files or the number of energies in your input file\n");
               	fprintf(stderr,"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n");
               	exit(1);
	}
	return order;
}
/************************************************************************************************************/
static int readEnergiesAndDipoles(char* fileName, double* energies[], double* dipoles[])
{
	int n=1;
	int xyz;
        char* TAG = NULL;
	int nGeoms = 0;
	int index = 0;
	int nDipoles = 0;
	FILE* file = fopen(fileName,"rb");
        TAG = strdup("ENERGIES");
        uppercase(TAG);
        rewind(file);
	if(!goToStr(file, TAG)) 
	{
		fprintf(stderr,"I cannot find %s in our file\n",TAG);
		if(TAG) free(TAG);
		return FALSE;
	}

        // computes the number of geometries 
        nGeoms = 0;
        for(index=0; ;index++)
        {
		double e;
		n = fscanf(file,"%lf",&e);
                if(n<1) break;
                nGeoms++;
        }
        if(nGeoms<1)
        {
                fprintf(stderr,"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n");
                fprintf(stderr,"I cannot read energies from %s file\n", fileName);
                fprintf(stderr,"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n");
                exit(1);
        }
	printf("nGeoms = %d\n",nGeoms);
        rewind(file);
	goToStr(file, TAG);
        energies[0] = malloc(nGeoms*sizeof(double));
        for(index=0;index<nGeoms ;index++)
        {
		double e;
		n = fscanf(file,"%lf",&e);
                if(n<1) break;
                energies[0][index] = e;
        }
	free(TAG);
        TAG = strdup("DIPOLES");
        uppercase(TAG);
        rewind(file);
	if(!goToStr(file, TAG)) 
	{
		fprintf(stderr,"I cannot find %s in our file\n",TAG);
		if(TAG) free(TAG);
		return 0;
	}
	nDipoles = 0;
        for(index=0; ;index++)
        {
		double e;
		n = 1;
		for(xyz=0;xyz<3 && n==1;xyz++)
		{
			n = fscanf(file,"%lf",&e);
                	if(n<1) break;
			nDipoles++;
		}
		if(n<1) break;
        }
	nDipoles /= 3;
        if(nGeoms!=nDipoles)
        {
                fprintf(stderr,"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n");
                fprintf(stderr,"Error : number of dipoles != number of geometries in %s file\n",fileName);
                fprintf(stderr,"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n");
                exit(1);
        }
	printf("nDipoles  = %d\n",nDipoles);
        rewind(file);
	goToStr(file, TAG);
	for(xyz=0;xyz<3;xyz++) dipoles[xyz] = malloc(nGeoms*sizeof(double));
        for(index=0;index<nGeoms ;index++)
        {
		double e;
		n = 1;
		for(xyz=0;xyz<3 && n==1;xyz++)
		{
			n = fscanf(file,"%lf",&e);
                	if(n<1) break;
			dipoles[xyz][index] = e;
		}
		if(n<1) break;
        }
	fclose(file);
	return nGeoms;
}
/************************************************************************************************************/
static void computeGradients(NMRQFF* qffConstants)
{
 	int nf = qffConstants->numberOfFrequencies;
	int i;
	for(i=0;i<nf;i++)
	{
		double di = qffConstants->delta[i]*sqrt(qffConstants->mass[i]*AMUTOAU);
		qffConstants->gradients[i]=(
				+qffConstants->VI[i][0]+
				-9*qffConstants->VI[i][1]+
				+45*qffConstants->VI[i][2]+
				-45*qffConstants->VI[i][3]+
				+9*qffConstants->VI[i][4]
				-qffConstants->VI[i][5]
				)/(60.0*di);
	}
}
/************************************************************************************************************/
static void computeFrequencies(NMRQFF* qffConstants)
{
 	int nf = qffConstants->numberOfFrequencies;
	int i;
	for(i=0;i<nf;i++)
	{
		double di = qffConstants->delta[i]*qffConstants->delta[i]*AMUTOAU*qffConstants->mass[i];
		double f = 1/di/180;
		qffConstants->calculatedFrequencies[i]=AUTOCM1*sqrt(f*fabs((
		2*qffConstants->VI[i][0]+
		-27*qffConstants->VI[i][1]+
		+270*qffConstants->VI[i][2]+
		-490*qffConstants->V0+
		+270*qffConstants->VI[i][3]+
		-27*qffConstants->VI[i][4]+
		2*qffConstants->VI[i][5]
		)));
		/*
		qffConstants->calculatedFrequencies[i]=sqrt(fabs((
		  qffConstants->VI[i][2]
		- 2*qffConstants->V0
		+ qffConstants->VI[i][3]
		)
		/di))*AUTOCM1;
		*/
	}
}
/************************************************************************************************************/
static void computeCubicForces(NMRQFF* qffConstants)
{
 	int nf = qffConstants->numberOfFrequencies;
	int i,j,k;
	double f3cm1 = AUTOCM1*sqrt(AUTOCM1)*AUTOCM1;
	if(nf<1) return;
	// tiii

	for(i=0;i<nf;i++)
	{
		double mi = sqrt(qffConstants->mass[i]*AMUTOAU);
		double f = 1.0/(8.0*qffConstants->delta[i]*qffConstants->delta[i]*qffConstants->delta[i]*mi*mi*mi);
		f = f/qffConstants->frequencies[i]/sqrt(qffConstants->frequencies[i])*f3cm1;
		qffConstants->cubicEnergies[i][i] = f*(-qffConstants->VI[i][0]+8*qffConstants->VI[i][1]-13*qffConstants->VI[i][2]+13*qffConstants->VI[i][3]-8*qffConstants->VI[i][4]+qffConstants->VI[i][5]);
		/*
		double mi = sqrt(qffConstants->mass[i]);
		double f = 1.0/(2.0*qffConstants->delta[i]*qffConstants->delta[i]*qffConstants->delta[i]*mi*mi*mi);
		qffConstants->cubicEnergies[i][i] = f*(-1*qffConstants->VI[i][1]+2*qffConstants->VI[i][2]-2*qffConstants->VI[i][3]+1*qffConstants->VI[i][4]);
		*/
	}
	// tiij
	if(qffConstants->numberOfEnergies > 1+6*nf)
	for(i=0;i<nf;i++)
	{
		double fi = 1.0/(2.0*qffConstants->delta[i]*qffConstants->delta[i]*qffConstants->mass[i]*AMUTOAU);
		fi = fi/qffConstants->frequencies[i];
		for(j=0;j<nf;j++)
		{
			double fj,f;
			if(j==i) continue;
			fj = 1.0/(qffConstants->delta[j]*sqrt(qffConstants->mass[j]*AMUTOAU));
			fj = fj/sqrt(qffConstants->frequencies[j]);
			f = fi*fj*f3cm1;
			qffConstants->cubicEnergies[i][j] = f*(
				 (qffConstants->VIJ[i][j][0]+qffConstants->VIJ[i][j][2]-2*qffConstants->VI[j][2])
				-(qffConstants->VIJ[i][j][1]+qffConstants->VIJ[i][j][3]-2*qffConstants->VI[j][3])
			);

		}
	}
	// tijk , i # j # k
	if(qffConstants->numberOfEnergies > 1+6*nf*nf)
	{
	FILE* fileVIJK = fopen("VIJK.bin","rb");
	for(i=0;i<nf;i++)
	{
		double fi;
		fi = 1.0/(qffConstants->delta[i]*sqrt(qffConstants->mass[i]*AMUTOAU));
		fi = fi/sqrt(qffConstants->frequencies[i]);
		for(j=0;j<i;j++)
		{
			double fj;
			fj = 1.0/(qffConstants->delta[j]*sqrt(qffConstants->mass[j]*AMUTOAU));
			fj = fj/sqrt(qffConstants->frequencies[j]);
			for(k=0;k<j;k++)
			{
				double VIJK[8];
				int l;
				double fk,f;
				fk = 1.0/(qffConstants->delta[k]*sqrt(qffConstants->mass[k]*AMUTOAU));
				fk = fk/sqrt(qffConstants->frequencies[k]);
				f = fi*fj*fk/8;
				f =f*f3cm1;

				for(l=0;l<8;l++) fread(&VIJK[l], sizeof(double), 1, fileVIJK);
				/* VIJK */
				/* 0  deltas[j],   deltas[i],  deltas[k] */
				/* 1  deltas[j],   deltas[i], -deltas[k] */
				/* 2  deltas[j],  -deltas[i],  deltas[k] */
				/* 3  deltas[j],  -deltas[i], -deltas[k] */
				/* 4 -deltas[j],  deltas[i],  deltas[k] */
				/* 5 -deltas[j],  deltas[i], -deltas[k] */
				/* 6 -deltas[j], -deltas[i],  deltas[k] */
				/* 7 -deltas[j], -deltas[i], -deltas[k] */
				qffConstants->cubicEnergiesIJK[i][j][k] = f*(
				+VIJK[0]
				-VIJK[1]
				-VIJK[2]
				+VIJK[3]
				-VIJK[4]
				+VIJK[5]
				+VIJK[6]
				-VIJK[7]
				/*
				+qffConstants->VIJK[i][j][k][0]
				-qffConstants->VIJK[i][j][k][1]
				-qffConstants->VIJK[i][j][k][2]
				+qffConstants->VIJK[i][j][k][3]
				-qffConstants->VIJK[i][j][k][4]
				+qffConstants->VIJK[i][j][k][5]
				+qffConstants->VIJK[i][j][k][6]
				-qffConstants->VIJK[i][j][k][7]
				*/
				);
			}
		}
	}
	fclose(fileVIJK);
	}
}
/************************************************************************************************************/
static void computeQuarticForces(NMRQFF* qffConstants)
{
 	int nf = qffConstants->numberOfFrequencies;
	double f4cm1 = AUTOCM1*AUTOCM1*AUTOCM1;
	int i,j,k,l;
	if(nf<1) return;

	// uiiii
	for(i=0;i<nf;i++)
	{
		double mdi = sqrt(qffConstants->mass[i]*AMUTOAU)*qffConstants->delta[i];
		double f = 1.0/(6.0*mdi*mdi*mdi*mdi)*f4cm1;
		f = f/(qffConstants->frequencies[i]*qffConstants->frequencies[i]);
		qffConstants->quarticEnergiesIIJJ[i][i] = f*(
		-qffConstants->VI[i][0]
		+12*qffConstants->VI[i][1]
		-39*qffConstants->VI[i][2]
		+56*qffConstants->V0
		-39*qffConstants->VI[i][3]
		+12*qffConstants->VI[i][4]
		-qffConstants->VI[i][5]
		);
		/*
		double f = 1.0/(mdi*mdi*mdi*mdi)*f4cm1;
		f = f/(qffConstants->frequencies[i]*qffConstants->frequencies[i]);
		qffConstants->quarticEnergiesIIJJ[i][i] = f*(
		qffConstants->VI[i][1]
		-4*qffConstants->VI[i][2]
		+6*qffConstants->V0
		-4*qffConstants->VI[i][3]
		+qffConstants->VI[i][4]
		);
		*/
		qffConstants->quarticEnergiesIIIJ[i][i] = qffConstants->quarticEnergiesIIJJ[i][i];
	}
	// uiiij
	if(qffConstants->numberOfEnergies > 1+6*nf)
	for(i=0;i<nf;i++)
	{
		double mdi = sqrt(qffConstants->mass[i]*AMUTOAU)*qffConstants->delta[i];
		double fi = 1.0/(16.0*mdi*mdi*mdi);
		fi = fi/(qffConstants->frequencies[i]*sqrt(qffConstants->frequencies[i]));
		for(j=0;j<nf;j++)
		{
			double fj,f;
			if(j==i) continue;
			fj = 1.0/(qffConstants->delta[j]*sqrt(qffConstants->mass[j]*AMUTOAU));
			fj = fj/(sqrt(qffConstants->frequencies[j]));
			f = fi*fj*f4cm1;
			qffConstants->quarticEnergiesIIIJ[i][j] = f*(
				 (qffConstants->VI3J[i][j][0]-3*qffConstants->VIJ[i][j][0]+3*qffConstants->VIJ[i][j][2]-qffConstants->VI3J[i][j][2])
				-(qffConstants->VI3J[i][j][1]-3*qffConstants->VIJ[i][j][1]+3*qffConstants->VIJ[i][j][3]-qffConstants->VI3J[i][j][3])
			);
		}
	}
//printf("UIIJJ======================\n");
	// uiijj
	if(qffConstants->numberOfEnergies > 1+6*nf)
	for(i=0;i<nf;i++)
	{
		double mdi = qffConstants->mass[i]*AMUTOAU*qffConstants->delta[i]*qffConstants->delta[i];
		for(j=0;j<i;j++)
		{
			double mdj = qffConstants->mass[j]*AMUTOAU*qffConstants->delta[j]*qffConstants->delta[j];
			double f = 1.0/(mdi*mdj)*f4cm1;
			f = f/(qffConstants->frequencies[i]*qffConstants->frequencies[j]);

			qffConstants->quarticEnergiesIIJJ[i][j] = f*(
				   (qffConstants->VIJ[i][j][0]+qffConstants->VIJ[i][j][2]+qffConstants->VIJ[i][j][1]+qffConstants->VIJ[i][j][3])
				-2*(qffConstants->VI[i][2]+qffConstants->VI[i][3]+qffConstants->VI[j][2]+qffConstants->VI[j][3])
				+4*qffConstants->V0
			);
			qffConstants->quarticEnergiesIIJJ[j][i] = qffConstants->quarticEnergiesIIJJ[i][j];
		}
	}
	// uiijk i # j # k
	if(qffConstants->numberOfEnergies > 1+6*nf*nf)
	{
	FILE* fileVIJK = fopen("VIJK.bin","rb");
	for(i=0;i<nf;i++)
	{
		double mdi = sqrt(qffConstants->mass[i]*AMUTOAU)*qffConstants->delta[i];
		double fi = 1.0/(4.0*mdi*mdi);
		fi = fi/(qffConstants->frequencies[i]*sqrt(qffConstants->frequencies[i]));
		for(j=0;j<i;j++)
		{
			double fj;
			fj = 1.0/(qffConstants->delta[j]*sqrt(qffConstants->mass[j]*AMUTOAU));
			fj = fj/(sqrt(qffConstants->frequencies[j]));
			for(k=0;k<j;k++)
			{
				double VIJK[8];
				int l;
				double fk,f;
				fk = 1.0/(qffConstants->delta[k]*sqrt(qffConstants->mass[k]*AMUTOAU));
				fk = fk/(sqrt(qffConstants->frequencies[k]));
				f = fi*fj*fk*f4cm1;

				for(l=0;l<8;l++) fread(&VIJK[l], sizeof(double), 1, fileVIJK);
				/* VIJK */
				/* 0  deltas[j],   deltas[i],  deltas[k] */
				/* 1  deltas[j],   deltas[i], -deltas[k] */
				/* 2  deltas[j],  -deltas[i],  deltas[k] */
				/* 3  deltas[j],  -deltas[i], -deltas[k] */
				/* 4 -deltas[j],  deltas[i],  deltas[k] */
				/* 5 -deltas[j],  deltas[i], -deltas[k] */
				/* 6 -deltas[j], -deltas[i],  deltas[k] */
				/* 7 -deltas[j], -deltas[i], -deltas[k] */

				/* VIJ */
				/*  0  deltas[j],   deltas[i] */
				/*  1  deltas[j],  -deltas[i] */
				/*  2 -deltas[j],   deltas[i] */
				/*  3 -deltas[j],  -deltas[i] */
				qffConstants->quarticEnergiesIIJK[i][j][k] = 
				f*(
				+VIJK[0]
				-VIJK[2]
				-VIJK[1]
				+VIJK[3]
				+VIJK[4]
				-VIJK[6]
				-VIJK[5]
				+VIJK[7]
				/*
				+qffConstants->VIJK[i][j][k][0]
				-qffConstants->VIJK[i][j][k][2]
				-qffConstants->VIJK[i][j][k][1]
				+qffConstants->VIJK[i][j][k][3]
				+qffConstants->VIJK[i][j][k][4]
				-qffConstants->VIJK[i][j][k][6]
				-qffConstants->VIJK[i][j][k][5]
				+qffConstants->VIJK[i][j][k][7]
				*/
				-2*(qffConstants->VIJ[i][j][0] - qffConstants->VIJ[i][j][2] - qffConstants->VIJ[i][j][1] + qffConstants->VIJ[i][j][3])
				);
			}
		}
	}
	fclose(fileVIJK);
	}
	// uijkl i # j # k # l
	if(qffConstants->numberOfEnergies > 1+6*nf*nf+8*nf*(nf-1)*(nf-2)/6)
	for(i=0;i<nf;i++)
	{
		double fi;
		fi = 1.0/(qffConstants->delta[i]*sqrt(qffConstants->mass[i]*AMUTOAU));
		fi = fi/(sqrt(qffConstants->frequencies[i]));
		for(j=0;j<i;j++)
		{
			double fj;
			fj = 1.0/(qffConstants->delta[j]*sqrt(qffConstants->mass[j]*AMUTOAU));
			fj = fj/(sqrt(qffConstants->frequencies[j]));
			for(k=0;k<j;k++)
			{
				double fk;
				fk = 1.0/(qffConstants->delta[k]*sqrt(qffConstants->mass[k]*AMUTOAU));
				fk = fk/(sqrt(qffConstants->frequencies[k]));
				for(l=0;l<k;l++)
				{
				double fl,f;
				fl = 1.0/(qffConstants->delta[l]*sqrt(qffConstants->mass[l]*AMUTOAU));
				f = fi*fj*fk*fl*f4cm1;
				f = f/16.0;
				/* VIJKL */
				/* 0    deltas[j],   deltas[i],   deltas[k],  deltas[l] */ //+
				/* 1    deltas[j],   deltas[i],   deltas[k], -deltas[l] */ //-
				/* 2    deltas[j],   deltas[i],  -deltas[k],  deltas[l] */ //-
				/* 3    deltas[j],   deltas[i],  -deltas[k], -deltas[l] */ //+
				/* 4    deltas[j],  -deltas[i],   deltas[k],  deltas[l] */ //-
				/* 5    deltas[j],  -deltas[i],   deltas[k], -deltas[l] */ //+
				/* 6    deltas[j],  -deltas[i],  -deltas[k],  deltas[l] */ //+
				/* 7    deltas[j],  -deltas[i],  -deltas[k], -deltas[l] */ //-
				/* 8   -deltas[j],   deltas[i],   deltas[k],  deltas[l] */ //-
				/* 9   -deltas[j],   deltas[i],   deltas[k], -deltas[l] */ //+
				/* 10  -deltas[j],   deltas[i],  -deltas[k],  deltas[l] */ //+
				/* 11  -deltas[j],   deltas[i],  -deltas[k], -deltas[l] */ //-
				/* 12  -deltas[j],  -deltas[i],   deltas[k],  deltas[l] */ //+
				/* 13  -deltas[j],  -deltas[i],   deltas[k], -deltas[l] */ //-
				/* 14  -deltas[j],  -deltas[i],  -deltas[k],  deltas[l] */ //-
				/* 15  -deltas[j],  -deltas[i],  -deltas[k], -deltas[l] */ //+

				qffConstants->quarticEnergiesIJKL[i][j][k][l]= 
				f*(
				+qffConstants->VIJKL[i][j][k][l][0]
				-qffConstants->VIJKL[i][j][k][l][1]
				-qffConstants->VIJKL[i][j][k][l][2]
				+qffConstants->VIJKL[i][j][k][l][3]
				-qffConstants->VIJKL[i][j][k][l][4]
				+qffConstants->VIJKL[i][j][k][l][5]
				+qffConstants->VIJKL[i][j][k][l][6]
				-qffConstants->VIJKL[i][j][k][l][7]
				-qffConstants->VIJKL[i][j][k][l][8]
				+qffConstants->VIJKL[i][j][k][l][9]
				+qffConstants->VIJKL[i][j][k][l][10]
				-qffConstants->VIJKL[i][j][k][l][11]
				+qffConstants->VIJKL[i][j][k][l][12]
				-qffConstants->VIJKL[i][j][k][l][13]
				-qffConstants->VIJKL[i][j][k][l][14]
				+qffConstants->VIJKL[i][j][k][l][15]
				);
				}
			}
		}
	}
}
/************************************************************************************************************/
static void computeFirstDerivativesDipoles(NMRQFF* qffConstants)
{
 	int nf = qffConstants->numberOfFrequencies;
	int i;
	int xyz = 0;
	for(i=0;i<nf;i++)
	{
		double di = 60*qffConstants->delta[i]*sqrt(qffConstants->mass[i]*AMUTOAU);
		double f = 1/di;
		f =f*sqrt(AUTOCM1);
		for(xyz=0;xyz<3;xyz++)
			qffConstants->firstDipoles[i][xyz]=f*(
				+qffConstants->dipolesI[i][0][xyz]+
				-9*qffConstants->dipolesI[i][1][xyz]+
				+45*qffConstants->dipolesI[i][2][xyz]+
				-45*qffConstants->dipolesI[i][3][xyz]+
				+9*qffConstants->dipolesI[i][4][xyz]
				-qffConstants->dipolesI[i][5][xyz]
				);
	}
}
/************************************************************************************************************/
static void changeUnitInputFirstDerivativesDipoles(NMRQFF* qffConstants)
{
        double mu0 = 4*PI*1e-7;
        double eps0 = 1.0/(mu0*slight*slight);
        double   kmmolm1 = 4*PI*PI*PI*NAvogadro/3/hPlank/slight/4/PI/eps0*1e-3*100.0*8.47835267e-30*8.47835267e-30;/* 1e-3 m to km, 100 : cm-1 to m-1 */
	double f = 1.0/sqrt(kmmolm1);

        int nf = qffConstants->numberOfFrequencies;
        int i;
        int xyz = 0;
	if(qffConstants->numberOfFirstDipolesInput==0) return;
        for(i=0;i<nf;i++)
        {
                for(xyz=0;xyz<3;xyz++)
                        qffConstants->firstDipolesInput[i][xyz] *=f;
        }
}
/************************************************************************************************************/
static void computeSecondDerivativesDipoles(NMRQFF* qffConstants)
{
 	int nf = qffConstants->numberOfFrequencies;
	int i;
	int j;
	int xyz = 0;
	// dii
	for(i=0;i<nf;i++)
	{
		double di = qffConstants->delta[i]*qffConstants->delta[i]*AMUTOAU*qffConstants->mass[i];
		double f = 1/(180*di);
		f =f*(AUTOCM1);
		for(xyz=0;xyz<3;xyz++)
			qffConstants->secondDipoles[i][i][xyz]=f*(
				 2*qffConstants->dipolesI[i][0][xyz]
				-27*qffConstants->dipolesI[i][1][xyz]
				+270*qffConstants->dipolesI[i][2][xyz]
				-490*qffConstants->dipole0[xyz]
				+270*qffConstants->dipolesI[i][3][xyz]
				-27*qffConstants->dipolesI[i][4][xyz]
				+2*qffConstants->dipolesI[i][5][xyz]
				);
			/*
			qffConstants->secondDipoles[i][i][xyz]=180*f*(
				+qffConstants->dipolesI[i][2][xyz]
				-2*qffConstants->dipole0[xyz]
				+qffConstants->dipolesI[i][3][xyz]
				);
			*/
	}
	// dij
	if(qffConstants->numberOfEnergies > 1+6*nf)
	for(i=0;i<nf;i++)
	{
		double di = qffConstants->delta[i]*sqrt(qffConstants->mass[i]*AMUTOAU);
		for(j=0;j<nf;j++)
		{
			double dj = qffConstants->delta[j]*sqrt(qffConstants->mass[j]*AMUTOAU);
			double f = 1.0/(4.0*di*dj);
			if(i==j) continue;
			f =f*AUTOCM1;
			for(xyz=0;xyz<3;xyz++)
				qffConstants->secondDipoles[i][j][xyz]=f*(
				  qffConstants->dipolesIJ[i][j][0][xyz]
				- qffConstants->dipolesIJ[i][j][1][xyz]
				- qffConstants->dipolesIJ[i][j][2][xyz]
				+ qffConstants->dipolesIJ[i][j][3][xyz]
				);
			//for(xyz=0;xyz<3;xyz++)
				//qffConstants->secondDipoles[j][i][xyz]=qffConstants->secondDipoles[i][j][xyz];
		}
	}
}
/************************************************************************************************************/
static void computeCubicDerivativesDipoles(NMRQFF* qffConstants)
{
 	int nf = qffConstants->numberOfFrequencies;
	int i,j,k;
	int xyz;
	if(nf<1) return;
	// diii
	for(i=0;i<nf;i++)
	{
		double mi = sqrt(qffConstants->mass[i]*AMUTOAU);
		double f = 1.0/(8.0*qffConstants->delta[i]*qffConstants->delta[i]*qffConstants->delta[i]*mi*mi*mi);
		f =f*AUTOCM1*sqrt(AUTOCM1);
		for(xyz=0;xyz<3;xyz++)
			qffConstants->cubicDipoles[i][i][xyz] = f*
			(
			-   qffConstants->dipolesI[i][0][xyz]
			+8* qffConstants->dipolesI[i][1][xyz]
			-13*qffConstants->dipolesI[i][2][xyz]
			+13*qffConstants->dipolesI[i][3][xyz]
			-8* qffConstants->dipolesI[i][4][xyz]
			+   qffConstants->dipolesI[i][5][xyz]
			);
		//qffConstants->cubicDipoles[i][i][xyz] = f*(qffConstants->dipolesI[i][0][xyz]-3*qffConstants->dipolesI[i][2][xyz]+3*qffConstants->dipolesI[i][3][xyz]-qffConstants->dipolesI[i][5][xyz]);
	}
	// diij
	if(qffConstants->numberOfEnergies > 1+6*nf)
	for(i=0;i<nf;i++)
	{
		double fi = 1.0/(2.0*qffConstants->delta[i]*qffConstants->delta[i]*qffConstants->mass[i]*AMUTOAU);
		for(j=0;j<nf;j++)
		{
			double fj,f;
			if(j==i) continue;
			fj = 1.0/(qffConstants->delta[j]*sqrt(qffConstants->mass[j]*AMUTOAU));
			f = fi*fj;
			f =f*AUTOCM1*sqrt(AUTOCM1);
			for(xyz=0;xyz<3;xyz++)
				qffConstants->cubicDipoles[i][j][xyz] = f*(
				 (qffConstants->dipolesIJ[i][j][0][xyz]+qffConstants->dipolesIJ[i][j][2][xyz]-2*qffConstants->dipolesI[j][2][xyz])
				-(qffConstants->dipolesIJ[i][j][1][xyz]+qffConstants->dipolesIJ[i][j][3][xyz]-2*qffConstants->dipolesI[j][3][xyz])
			);
		}
	}
	// dijk , i # j # k
	if(qffConstants->numberOfEnergies > 1+6*nf*nf)
	{
	FILE* fileDIJK = fopen("DIJK.bin","rb");
	for(i=0;i<nf;i++)
	{
		double fi;
		fi = 1.0/(qffConstants->delta[i]*sqrt(qffConstants->mass[i]*AMUTOAU));
		for(j=0;j<i;j++)
		{
			double fj;
			fj = 1.0/(qffConstants->delta[j]*sqrt(qffConstants->mass[j]*AMUTOAU));
			for(k=0;k<j;k++)
			{
				double fk,f;
				fk = 1.0/(qffConstants->delta[k]*sqrt(qffConstants->mass[k]*AMUTOAU));
				f = fi*fj*fk/8;
				f =f*AUTOCM1*sqrt(AUTOCM1);
				/* VIJK */
				/* 0  deltas[j],   deltas[i],  deltas[k] */
				/* 1  deltas[j],   deltas[i], -deltas[k] */
				/* 2  deltas[j],  -deltas[i],  deltas[k] */
				/* 3  deltas[j],  -deltas[i], -deltas[k] */
				/* 4 -deltas[j],  deltas[i],  deltas[k] */
				/* 5 -deltas[j],  deltas[i], -deltas[k] */
				/* 6 -deltas[j], -deltas[i],  deltas[k] */
				/* 7 -deltas[j], -deltas[i], -deltas[k] */
				double dipolesIJK[8][3];
				int l;
				for(l=0;l<8;l++) 
				for(xyz=0;xyz<3;xyz++)
					fread(&dipolesIJK[l][xyz], sizeof(double), 1, fileDIJK);

				for(xyz=0;xyz<3;xyz++)
				{

				qffConstants->cubicDipolesIJK[i][j][k][xyz] = f*(
				+dipolesIJK[0][xyz]
				-dipolesIJK[1][xyz]
				-dipolesIJK[2][xyz]
				+dipolesIJK[3][xyz]
				-dipolesIJK[4][xyz]
				+dipolesIJK[5][xyz]
				+dipolesIJK[6][xyz]
				-dipolesIJK[7][xyz]
				/*
				+qffConstants->dipolesIJK[i][j][k][0][xyz]
				-qffConstants->dipolesIJK[i][j][k][1][xyz]
				-qffConstants->dipolesIJK[i][j][k][2][xyz]
				+qffConstants->dipolesIJK[i][j][k][3][xyz]
				-qffConstants->dipolesIJK[i][j][k][4][xyz]
				+qffConstants->dipolesIJK[i][j][k][5][xyz]
				+qffConstants->dipolesIJK[i][j][k][6][xyz]
				-qffConstants->dipolesIJK[i][j][k][7][xyz]
				*/
				);
				}
			}
		}
	}
	fclose(fileDIJK);
	for(j=0;j<nf;j++)
	{
		for(i=0;i<j;i++)
		{
			for(k=0;k<i;k++)
			{
				for(xyz=0;xyz<3;xyz++)
                		{
					qffConstants->cubicDipolesIJK[i][j][k][xyz] = qffConstants->cubicDipolesIJK[j][i][k][xyz];
					qffConstants->cubicDipolesIJK[i][k][j][xyz] = qffConstants->cubicDipolesIJK[j][i][k][xyz];
					qffConstants->cubicDipolesIJK[j][k][i][xyz] = qffConstants->cubicDipolesIJK[j][i][k][xyz];
					qffConstants->cubicDipolesIJK[k][i][j][xyz] = qffConstants->cubicDipolesIJK[j][i][k][xyz];
					qffConstants->cubicDipolesIJK[k][j][i][xyz] = qffConstants->cubicDipolesIJK[j][i][k][xyz];
                		}
			}
		}
	}
	}
}
/************************************************************************************************************/
static void computeQuarticDerivativesDipoles(NMRQFF* qffConstants)
{
 	int nf = qffConstants->numberOfFrequencies;
	int i,j,k;
	int xyz;
	if(nf<1) return;

	// diiii
	for(i=0;i<nf;i++)
	{
		double mdi = sqrt(qffConstants->mass[i]*AMUTOAU)*qffConstants->delta[i];
		double f = 1.0/(6*mdi*mdi*mdi*mdi);
		f =f*AUTOCM1*AUTOCM1;
		for(xyz=0;xyz<3;xyz++)
		{
			qffConstants->quarticDipolesIIJJ[i][i][xyz] = f*(
		   	-qffConstants->dipolesI[i][0][xyz]
		   	+12*qffConstants->dipolesI[i][1][xyz]
			-39*qffConstants->dipolesI[i][2][xyz]
			+56*qffConstants->dipole0[xyz]
			-39*qffConstants->dipolesI[i][3][xyz]
			+12*qffConstants->dipolesI[i][4][xyz]
			-qffConstants->dipolesI[i][5][xyz]);
			qffConstants->quarticDipolesIIIJ[i][i][xyz] = qffConstants->quarticDipolesIIJJ[i][i][xyz];
		}
	}
	// diiij
	if(qffConstants->numberOfEnergies > 1+6*nf)
	for(i=0;i<nf;i++)
	{
		double mdi = sqrt(qffConstants->mass[i]*AMUTOAU)*qffConstants->delta[i];
		double fi = 1.0/(16.0*mdi*mdi*mdi);
		for(j=0;j<nf;j++)
		{
			double fj,f;
			if(j==i) continue;
			fj = 1.0/(qffConstants->delta[j]*sqrt(qffConstants->mass[j]*AMUTOAU));
			f = fi*fj;
			f =f*AUTOCM1*AUTOCM1;
			for(xyz=0;xyz<3;xyz++)
				qffConstants->quarticDipolesIIIJ[i][j][xyz] = f*(
				 (	qffConstants->dipolesI3J[i][j][0][xyz]
					-3*qffConstants->dipolesIJ[i][j][0][xyz]
					+3*qffConstants->dipolesIJ[i][j][2][xyz]
					-qffConstants->dipolesI3J[i][j][2][xyz]
				)
				-(
					qffConstants->dipolesI3J[i][j][1][xyz]
					-3*qffConstants->dipolesIJ[i][j][1][xyz]
					+3*qffConstants->dipolesIJ[i][j][3][xyz]
					-qffConstants->dipolesI3J[i][j][3][xyz]
				)
				);
		}
	}
	// diijj
	if(qffConstants->numberOfEnergies > 1+6*nf)
	for(i=0;i<nf;i++)
	{
		double mdi = qffConstants->mass[i]*AMUTOAU*qffConstants->delta[i]*qffConstants->delta[i];
		for(j=0;j<i;j++)
		{
			double mdj = qffConstants->mass[j]*AMUTOAU*qffConstants->delta[j]*qffConstants->delta[j];
			double f = 1.0/(mdi*mdj);
			f =f*AUTOCM1*AUTOCM1;
			for(xyz=0;xyz<3;xyz++)
			{
				qffConstants->quarticDipolesIIJJ[i][j][xyz] = f*(
				   (
					qffConstants->dipolesIJ[i][j][0][xyz]
					+qffConstants->dipolesIJ[i][j][2][xyz]
					+qffConstants->dipolesIJ[i][j][1][xyz]
					+qffConstants->dipolesIJ[i][j][3][xyz]
				   )
				-2*(
					qffConstants->dipolesI[i][2][xyz]
					+qffConstants->dipolesI[i][3][xyz]
					+qffConstants->dipolesI[j][2][xyz]
					+qffConstants->dipolesI[j][3][xyz]
				)
				+4*qffConstants->dipole0[xyz]
				);
				qffConstants->quarticDipolesIIJJ[j][i][xyz] = qffConstants->quarticDipolesIIJJ[i][j][xyz];
			}
		}
	}
	// diijk i # j # k
	if(qffConstants->numberOfEnergies > 1+6*nf*nf)
	{
	FILE* fileDIJK = fopen("DIJK.bin","rb");
	for(i=0;i<nf;i++)
	{
		double mdi = sqrt(qffConstants->mass[i]*AMUTOAU)*qffConstants->delta[i];
		double fi = 1.0/(4.0*mdi*mdi);
		for(j=0;j<i;j++)
		{
			double fj;
			fj = 1.0/(qffConstants->delta[j]*sqrt(qffConstants->mass[j]*AMUTOAU));
			for(k=0;k<j;k++)
			{
				double fk,f;
				fk = 1.0/(qffConstants->delta[k]*sqrt(qffConstants->mass[k]*AMUTOAU));
				f = fi*fj*fk;
				f =f*AUTOCM1*AUTOCM1;
				/* quarticDipolesIIJK */
				/* 0  deltas[j],   deltas[i],  deltas[k] */
				/* 1  deltas[j],   deltas[i], -deltas[k] */
				/* 2  deltas[j],  -deltas[i],  deltas[k] */
				/* 3  deltas[j],  -deltas[i], -deltas[k] */
				/* 4 -deltas[j],  deltas[i],  deltas[k] */
				/* 5 -deltas[j],  deltas[i], -deltas[k] */
				/* 6 -deltas[j], -deltas[i],  deltas[k] */
				/* 7 -deltas[j], -deltas[i], -deltas[k] */

				/* dipolesIJ */
				/*  0  deltas[j],   deltas[i] */
				/*  1  deltas[j],  -deltas[i] */
				/*  2 -deltas[j],   deltas[i] */
				/*  3 -deltas[j],  -deltas[i] */
				for(xyz=0;xyz<3;xyz++)
				{
					double dipolesIJK[8];
					int l;
					for(l=0;l<8;l++) fread(&dipolesIJK[l], sizeof(double), 1, fileDIJK);

					qffConstants->quarticDipolesIIJK[i][j][k][xyz] = f*
					(
				 	+dipolesIJK[0]
				 	-dipolesIJK[2]
				 	-dipolesIJK[1]
				 	+dipolesIJK[3]
				 	+dipolesIJK[4]
				 	-dipolesIJK[6]
				 	-dipolesIJK[5]
				 	+dipolesIJK[7]
					/*
				 	+qffConstants->dipolesIJK[i][j][k][0][xyz]
				 	-qffConstants->dipolesIJK[i][j][k][2][xyz]
				 	-qffConstants->dipolesIJK[i][j][k][1][xyz]
				 	+qffConstants->dipolesIJK[i][j][k][3][xyz]
				 	+qffConstants->dipolesIJK[i][j][k][4][xyz]
				 	-qffConstants->dipolesIJK[i][j][k][6][xyz]
				 	-qffConstants->dipolesIJK[i][j][k][5][xyz]
				 	+qffConstants->dipolesIJK[i][j][k][7][xyz]
					*/
					-2*(
					+qffConstants->dipolesIJ[i][j][0][xyz]
					-qffConstants->dipolesIJ[i][j][2][xyz]
					-qffConstants->dipolesIJ[i][j][1][xyz]
					+qffConstants->dipolesIJ[i][j][3][xyz]
					)
					);
				}
			}
		}
	}
	fclose(fileDIJK);
	for(j=0;j<nf;j++)
	{
		for(i=0;i<j;i++)
		{
			for(k=0;k<i;k++)
			{
				for(xyz=0;xyz<3;xyz++)
                		{
					qffConstants->quarticDipolesIIJK[i][j][k][xyz] = qffConstants->quarticDipolesIIJK[j][i][k][xyz];
					qffConstants->quarticDipolesIIJK[i][k][j][xyz] = qffConstants->quarticDipolesIIJK[j][i][k][xyz];
					qffConstants->quarticDipolesIIJK[j][k][i][xyz] = qffConstants->quarticDipolesIIJK[j][i][k][xyz];
					qffConstants->quarticDipolesIIJK[k][i][j][xyz] = qffConstants->quarticDipolesIIJK[j][i][k][xyz];
					qffConstants->quarticDipolesIIJK[k][j][i][xyz] = qffConstants->quarticDipolesIIJK[j][i][k][xyz];
                		}
			}
		}
	}
	}
}
/************************************************************************************************************/
char* computeQFFFromEnergiesDipolesFile(char* fileName)
{
	NMRQFF* qffConstants;
	boolean Ok = TRUE;
	FILE* file = NULL;
	int nGeoms = 0;
	int order = 0;
	double* energies = NULL;
	double* dipoles[] = { NULL, NULL, NULL};
	char* fileNameOut = NULL;

 	file = fopen(fileName, "rb"); 
        if(!file) 
	{
		fprintf(stderr,"I cannot open the %s file\n",fileName);
		return FALSE;
	}

	qffConstants = (NMRQFF*) malloc(sizeof(NMRQFF));
	qffConstants->numberOfFrequencies = 0;
	Ok =  readFrequenciesInitNMRQFF(file, qffConstants);
	if(Ok) Ok = readVectorRealQFF(file, "Mass", qffConstants->numberOfFrequencies, qffConstants->mass);
	if(Ok) Ok = readVectorRealQFF(file, "Delta", qffConstants->numberOfFrequencies, qffConstants->delta);

	
	printf("Begin readEnergiesAndDipoles\n");
	nGeoms = readEnergiesAndDipoles(fileName, &energies, dipoles);
	printf("End readEnergiesAndDipoles\n");
	order =  getOrdre(qffConstants->numberOfFrequencies, nGeoms);
        qffConstants->numberOfEnergies = nGeoms;
	computeQFFDerivativesEnergies(qffConstants, order, energies);
	printf("Energ computeQFFDerivativesEnergies\n");
	computeQFFDerivativesDipoles(qffConstants, order, dipoles);
	printf("Energ computeQFFDerivativesDipoles\n");

	fileNameOut =  printnMRQFF(qffConstants, fileName);

	return fileNameOut;
}
/********************************************************************************/
static int readGradGabeditFile(Molecule* mol, char* inputFileName, int index, boolean readGrad)
{
        char* fileName = NULL;
        char* prefixName = NULL;
        int i=0;
	int xyz;

        if(!mol || mol->nAtoms<1) return 0;
        if(!inputFileName) return 0;

        prefixName = strdup_printf("%sQFF",getSuffixNameFile(inputFileName));

        fileName = strdup_printf("%s_%d.gab",prefixName,index);
	{
		FILE* file = fopen(fileName,"r");
		if(!file) return 0;
		fclose(file);
	}

	i = readEnergyAndDipoleFromGabeditFile(fileName, &mol->potentialEnergy, mol->dipole);
        //printf("FileName = %s E = %f D = %f %f %f\n", fileName, mol->potentialEnergy, mol->dipole[0],  mol->dipole[1],  mol->dipole[2]);
	if(readGrad) mol->klass->readGradientFromGabeditFile(mol, fileName);
        if(fileName) free(fileName);
	mol->potentialEnergy /= AUTOKCAL;
	for(xyz=0;xyz<3;xyz++)  mol->dipole[xyz] /=  AUTODEB;
        return i;
}
/*****************************************************************************/
static int readEnergiesAndDipolesQFFFromFiles(char* inputFileName, double* energies[], double* dipoles[])
{
	int xyz;
	int index;
	Molecule* mol;
	char* prefixName;
	char* fileName;
	int n;
	int nGeoms = 0;

	prefixName = strdup_printf("%sQFF",getSuffixNameFile(inputFileName));
        fileName = strdup_printf("%s_%d.gab",prefixName,0);

	printf("Reading molecule from %s file\n", fileName);
	mol = readMoleculeFromGabeditFile(fileName);
	if(!mol || mol->nAtoms<1) 
	{
		fprintf(stderr,"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n");
		fprintf(stderr," Error : I cannot read geometry from %s Gabedit file\n",fileName);
		fprintf(stderr,"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n");
		exit(1);
	}
	// computes the number of geometries 
	nGeoms = 0;
	for(index=0; ;index++)
	{
		n = readGradGabeditFile(mol, inputFileName, index, FALSE);
		if(n<1) break;
		nGeoms++;
	}
	if(nGeoms<1)
	{
               	fprintf(stderr,"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n");
               	fprintf(stderr,"I cannot read energy and dipole from %sQFF*.gab gabedit files\n", inputFileName);
               	fprintf(stderr,"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n");
               	exit(1);
	}
	energies[0] = malloc(nGeoms*sizeof(double));
	for(xyz=0;xyz<3;xyz++) dipoles[xyz] = malloc(nGeoms*sizeof(double));
	for(index=0;index<nGeoms ;index++)
	{
		n = readGradGabeditFile(mol, inputFileName, index, FALSE);
		if(n<1) break;
		energies[0][index] = mol->potentialEnergy;
		for(xyz=0;xyz<3;xyz++) dipoles[xyz][index] = mol->dipole[xyz];
	}
	return nGeoms;
}
/*****************************************************************************/
char* computeQFFFromFiles(char* inputFileName)
{
	int i;
	int index;
	Molecule* mol;
	char* prefixName;
	char* fileName;
	FILE* file  = NULL;
	char buffer[BSIZE];
	NMRQFF* qffConstants;
	int order = -1;
	double* energies = NULL;
	double* dipoles[] = {NULL, NULL, NULL};
	int f=0, nGeoms = 0;
	char* fileNameOut = NULL;

	prefixName = strdup_printf("%sQFF",getSuffixNameFile(inputFileName));
        fileName = strdup_printf("%s_%d.gab",prefixName,0);

	//printf("Reading molecule from %s file\n", fileName);
	mol = readMoleculeFromGabeditFile(fileName);
	if(!mol || mol->nAtoms<1) 
	{
		fprintf(stderr,"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n");
		fprintf(stderr," Error : I cannot read geometry from %s Gabedit file\n",fileName);
		fprintf(stderr,"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n");
		exit(1);
	}
	qffConstants = (NMRQFF*) malloc(sizeof(NMRQFF));
	qffConstants->numberOfFrequencies = 0;
	//printf("nAtoms = %d\n", mol->nAtoms);

	// compute number of modes. It is <= 3*mol->nAtoms
	for(i=1;i<=6*3*mol->nAtoms;i++)
	{
        	fileName = strdup_printf("%s_%d.ici",prefixName,i);
		//printf("fileName = %s\n", fileName);
		file = fopen(fileName,"r");
		while(file && !feof(file))
        	{
                	if(!fgets(buffer,BSIZE,file))break;
			//printf("buffer = %s\n", buffer);
                	if(strstr(buffer,"Mode:")) 
			{
				qffConstants->numberOfFrequencies++;
				break;
			}
		}
		if(file)fclose(file);
		else break;
		free(fileName);
	}
	qffConstants->numberOfFrequencies /= 6;
	//printf("numberOfFrequencies = %d\n", qffConstants->numberOfFrequencies);
	initnMRQFF0(qffConstants,qffConstants->numberOfFrequencies);

	index = 1;
	for(i=0;i<qffConstants->numberOfFrequencies;i++)
	{
        	fileName = strdup_printf("%s_%d.ici",prefixName,index);
		file = fopen(fileName,"r");
		while(file && !feof(file))
        	{
                	if(!fgets(buffer,BSIZE,file))break;
                	if(strstr(buffer,"Mode:")) 
			{
				double freq = 0;
				double mass = 0;
				double delta = 0;
				//#Mode: Freq= 625.021755910743 Mass= 1.043203662605 Q= Qeq + 0.644571737151 akI=0.000042
				sscanf(strstr(buffer,"Freq=")+strlen("Freq="),"%lf",&freq);
				sscanf(strstr(buffer,"Mass=")+strlen("Mass="),"%lf",&mass);
				sscanf(strstr(buffer,"Qeq +")+strlen("Qeq +"),"%lf",&delta);
				qffConstants->frequencies[i] = freq;
				qffConstants->mass[i] = mass;
				qffConstants->delta[i] = delta/3.0;
				break;
			}
		}
		fclose(file);
		free(fileName);
		index += 6;
	}
	nGeoms = readEnergiesAndDipolesQFFFromFiles(inputFileName, &energies, dipoles);
	f = qffConstants->numberOfFrequencies;
	order =  getOrdre(f, nGeoms);
        qffConstants->numberOfEnergies = nGeoms;
	computeQFFDerivativesEnergies(qffConstants, order, energies);
	//printf("Energ computeQFFDerivativesEnergies\n");
	computeQFFDerivativesDipoles(qffConstants, order, dipoles);
	//printf("Energ computeQFFDerivativesDipoles\n");

	fileNameOut = saveNMRQFFAppend(qffConstants, inputFileName,NULL);
	
	return fileNameOut;
}
/**********************************************************************************************************************/
void generateQFFCChemIFilesForFrequencies(char* inputFileName)
{
	Molecule* mol = readMolecule(inputFileName,TRUE);
	double delta = 0.5;
	boolean reducedCoordinates = TRUE;
	FILE* file = fopen(inputFileName,"rb");
	int order = 2;

	readOneBoolean(file,"QFFReducedCoordinates",&reducedCoordinates);
	readOneInt(file,"QFFnModes",&order);
	if(!reducedCoordinates) delta = 1e-2;
	readOneReal(file,"QFFDelta",&delta);
	fclose(file);
	printf("delta = %f ",delta);
	if(reducedCoordinates) printf("reducedCoordiantes\n");
	else printf(" Angshtrom\n");
	printf("nAtoms = %d\n",mol->nAtoms);
	printf("nModes = %d\n",mol->vibration.nModes);

	/*
	printf("seet test.gab\n");
	mol->klass->save(mol, "test.gab");
	*/

	mol->klass->generateQFFCChemIFiles(mol, inputFileName, delta,  reducedCoordinates, order);
	printf("order = %d\n",order);
	mol->klass->saveFirstDerivatives(inputFileName,mol);

	mol->klass->free(mol);
}
/**********************************************************************************************************************/
double*** getGeomsQFF(char* inputFileName, Molecule** molecule, int* pOrdre, double** pDeltas, int* nGeoms)
{
	Molecule* mol = readMolecule(inputFileName,TRUE);
	double delta = 0.5;
	boolean reducedCoordinates = TRUE;
	int order = 2;
	FILE* file = fopen(inputFileName,"rb");
	double*** geoms = NULL;

	readOneBoolean(file,"QFFReducedCoordinates",&reducedCoordinates);
	readOneInt(file,"QFFnModes",&order);
	if(!reducedCoordinates) delta = 1e-2;
	readOneReal(file,"QFFDelta",&delta);
	fclose(file);
	printf("delta = %f ",delta);
	if(reducedCoordinates) printf("reducedCoordiantes\n");
	else printf(" Angshtrom\n");
	printf("nAtoms = %d\n",mol->nAtoms);
	printf("nModes = %d\n",mol->vibration.nModes);

	/*
	printf("seet test.gab\n");
	mol->klass->save(mol, "test.gab");
	*/

	geoms = mol->klass->getGeomsQFF(mol, inputFileName, delta,  reducedCoordinates, order, pDeltas, nGeoms);
	*pOrdre = order;

	//mol->klass->free(mol);
	*molecule = mol;
	return geoms;
}
/********************************************************************************/
static char* saveNMRQFFAppend(NMRQFF* qffConstants, char* inputFileName, Molecule* mol)
{

	char tmp[BSIZE];
	int i,j,k,l;
	int nf = qffConstants->numberOfFrequencies;
	char txyz[]={'X','Y','Z'};
	int xyz;
        char* fileNameOut = strdup_printf("%sQFF.txt",getSuffixNameFile(inputFileName));
        FILE* file = fopen(fileNameOut,"w");
        FILE* fileIn = fopen(inputFileName,"rb");
        fprintf(stdout,"QFF parameters saved in %s file\n", fileNameOut);


	while(!feof(fileIn))
      	{
                if(!fgets(tmp,BSIZE,fileIn))break;
		fprintf(file,"%s",tmp);
        }
	fclose(fileIn);

	sprintf(tmp, "#====================================================================================\n");
	fprintf(file,"%s", tmp);   
	sprintf(tmp,"%s",
		"# nMR-QFF constants\n"
		"# See Yagi et al. J. Chem. Phys. 121, 1383 (2004)\n"
		);
	fprintf(file,"%s",tmp);   
	sprintf(tmp, "#====================================================================================\n");
	fprintf(file,"%s",tmp);   

	fprintf(file,"%s","\n");
	fprintf(file,"%s","VPT2Model=GVPT2\n");   
	fprintf(file,"%s","# VPT2Model=DCPT2\n");   
	fprintf(file,"%s","# VPT2Model=HDCPT2\n");   
	fprintf(file,"%s","# alphaHDCPT2=1.0\n");   
	fprintf(file,"%s","# betaHDCPT2=5e5\n");   
	fprintf(file,"%s","\n");
	fprintf(file,"%s","PropModel=GVPT2\n");
	fprintf(file,"%s","# PropModel=HDCPT2\n");
	fprintf(file,"%s","# PropModel=DCPT2\n");
	fprintf(file,"%s","# alphaPropHDCPT2=1.0\n");
	fprintf(file,"%s","# betaPropHDCPT2=5e5\n");
	fprintf(file,"%s","# alphaPropHDCPT2=1.0\n");
	fprintf(file,"%s","# betaPropHDCPT2=5e5\n");
	fprintf(file,"%s","maxFrequencyDifferenceFermi=200\n");
	fprintf(file,"%s","MartinCutOff1=1.0\n");
	fprintf(file,"%s","MartinCutOff2=1.0\n");
	fprintf(file,"%s","# ZCutOff=0.08\n");
	fprintf(file,"%s","\n");
	sprintf(tmp, "#====================================================================================\n");
	fprintf(file,"%s",tmp);   

	fprintf(file,"%s","\n");   
	sprintf(tmp,"nFrequencies=%d\n",nf);
	fprintf(file,"%s",tmp);   
	sprintf(tmp,"nDim=%d\n",3);
	fprintf(file,"%s",tmp);   
	fprintf(file,"%s","\n");   

	sprintf(tmp,"#i Freq(cm-1)  Calc.Freq   dQ(Bohr)  Mass(amu)\tGradient[ H amu^(-1/2) Bohr^(-1)]\n");
	fprintf(file,"%s",tmp);   
	fprintf(file,"%s","Hessian\n");   
	for(i=0;i<nf;i++) 
	{
		sprintf(tmp,"%d %d %0.14f %0.14f %0.14f %0.14f\t%0.14f\n",i+1, i+1, qffConstants->frequencies[i], 
				qffConstants->calculatedFrequencies[i], 
				qffConstants->delta[i], qffConstants->mass[i], qffConstants->gradients[i]);
		fprintf(file,"%s",tmp);   
	}
	fprintf(file,"%s","END\n\n");


	sprintf(tmp,"# i\tj\tk\tReduced values [cm-1]\n");
	fprintf(file,"%s",tmp);   
	fprintf(file,"%s","Cubic\n");   
	for(i=0;i<nf;i++)
        {
                for(j=0;j<nf;j++)
                {
			if(fabs(qffConstants->cubicEnergies[i][j])<1e-12) continue;
			sprintf(tmp,"%d\t%d\t%d\t%14.6f\n",i+1, i+1, j+1, qffConstants->cubicEnergies[i][j]);
			fprintf(file,"%s",tmp);   
                }
                for(j=0;j<i;j++)
                for(k=0;k<j;k++)
                {
			if(fabs(qffConstants->cubicEnergiesIJK[i][j][k])<1e-12) continue;
			sprintf(tmp,"%d\t%d\t%d\t%14.6f\n",i+1, j+1, k+1, qffConstants->cubicEnergiesIJK[i][j][k]);
			fprintf(file,"%s",tmp);   
                }
        }
	fprintf(file,"%s","END\n\n");


	sprintf(tmp,"# i\tj\tk\tl\tReduced values [cm-1]\n");
	fprintf(file,"%s",tmp);   
	fprintf(file,"%s","Quartic\n");   
	for(i=0;i<nf;i++)
        {
		for(j=0;j<=i;j++)
		{
			if(fabs(qffConstants->quarticEnergiesIIJJ[i][j])<1e-12) continue;
			sprintf(tmp,"%d\t%d\t%d\t%d\t%14.6f\n",i+1, i+1, j+1, j+1, qffConstants->quarticEnergiesIIJJ[i][j]);
			fprintf(file,"%s",tmp);   
		}
		for(j=0;j<nf;j++)
		{
			if(j==i) continue;
			if(fabs(qffConstants->quarticEnergiesIIIJ[i][j])<1e-12) continue;
			sprintf(tmp,"%d\t%d\t%d\t%d\t%14.6f\n",i+1, i+1, i+1, j+1, qffConstants->quarticEnergiesIIIJ[i][j]);
			fprintf(file,"%s",tmp);   
		}
                for(j=0;j<i;j++)
                for(k=0;k<j;k++)
		{
			if(fabs(qffConstants->quarticEnergiesIIJK[i][j][k])<1e-12) continue;
			sprintf(tmp,"%d\t%d\t%d\t%d\t%14.6f\n",i+1, i+1, j+1, k+1, qffConstants->quarticEnergiesIIJK[i][j][k]);
			fprintf(file,"%s",tmp);   
		}
                for(j=0;j<i;j++)
                for(k=0;k<j;k++)
                for(l=0;l<k;l++)
		{
			if(fabs(qffConstants->quarticEnergiesIJKL[i][j][k][l])<1e-12) continue;
			sprintf(tmp,"%d\t%d\t%d\t%d\t%14.6f\n",i+1, j+1, k+1, l+1, qffConstants->quarticEnergiesIJKL[i][j][k][l]);
			fprintf(file,"%s",tmp);   
		}
        }
	fprintf(file,"%s","END\n\n");

	if(mol)
	{
		mol->klass->addFirstDerivativeToFile(mol, file);
	}
	if(qffConstants->numberOfFirstDipolesInput==0)
	{
		sprintf(tmp,"#xyz\ti\tValues[au cm^1/2]\n");
		fprintf(file,"%s",tmp);   
		fprintf(file,"%s","First derivatives\n");
		for(i=0;i<nf;i++)
		for(xyz=0;xyz<3;xyz++)
		{
			if(fabs(qffConstants->firstDipoles[i][xyz])<1e-12) continue;
			sprintf(tmp,"%c\t%d\t%14.6f\n",txyz[xyz], i+1,qffConstants->firstDipoles[i][xyz]);
			fprintf(file,"%s",tmp);   
		}
		fprintf(file,"%s","END\n\n");
	}
	else
	{
		sprintf(tmp,"#xyz\ti\tInput values[au cm^1/2]\tCalculated values[au cm^1/2]\n");
		fprintf(file,"%s",tmp);   
		fprintf(file,"%s","First derivatives\n");
		for(i=0;i<nf;i++)
		for(xyz=0;xyz<3;xyz++)
		{
			if(fabs(qffConstants->firstDipolesInput[i][xyz])<1e-12) continue;
			sprintf(tmp,"%c\t%d\t%14.6f\t\t%14.6f\n",txyz[xyz], i+1,qffConstants->firstDipolesInput[i][xyz], qffConstants->firstDipoles[i][xyz]);
			fprintf(file,"%s",tmp);   
		}
		fprintf(file,"%s","END\n\n");
	}

	sprintf(tmp,"#xyz\ti\tj\tValues[au cm]\n");
	fprintf(file,"%s",tmp);   
	fprintf(file,"%s","Second derivatives\n");
	for(i=0;i<nf;i++)
	for(j=0;j<nf;j++)
	for(xyz=0;xyz<3;xyz++)
	{
		if(fabs(qffConstants->secondDipoles[i][j][xyz])<1e-12) continue;
		sprintf(tmp,"%c\t%d\t%d\t%14.6f\n",txyz[xyz], i+1,j+1,qffConstants->secondDipoles[i][j][xyz]);
		fprintf(file,"%s",tmp);   
	}
	fprintf(file,"%s","END\n\n");

	sprintf(tmp,"#xyz\ti\tj\tk\tValues[au cm^3/2]\n");
	fprintf(file,"%s",tmp);   
	fprintf(file,"%s","Cubic derivatives\n");
	for(i=0;i<nf;i++)
	for(j=0;j<nf;j++)
	for(xyz=0;xyz<3;xyz++)
	{
		if(fabs(qffConstants->cubicDipoles[i][j][xyz])<1e-12) continue;
		sprintf(tmp,"%c\t%d\t%d\t%d\t%14.6f\n",txyz[xyz], i+1,i+1,j+1,qffConstants->cubicDipoles[i][j][xyz]);
		fprintf(file,"%s",tmp);   
		if(i!=j)
		{
			sprintf(tmp,"%c\t%d\t%d\t%d\t%14.6f\n",txyz[xyz], i+1,j+1,i+1,qffConstants->cubicDipoles[i][j][xyz]);
			fprintf(file,"%s",tmp);   
			sprintf(tmp,"%c\t%d\t%d\t%d\t%14.6f\n",txyz[xyz], j+1,i+1,i+1,qffConstants->cubicDipoles[i][j][xyz]);
			fprintf(file,"%s",tmp);   
		}
	}
	fprintf(file,"%s","END\n\n");
        fclose(file);
	return fileNameOut;

}
/*****************************************************************************/
static void computeQFFDerivativesEnergies(NMRQFF* qffConstants, int order, double* energies)
{
	int i;
	int j;
	int k;
	int l;
	int index;
	boolean Ok = TRUE;

	initnMRQFFEnergies(qffConstants, order);

	index = 0;
        qffConstants->V0 =energies[index];

	for(j=0;j<qffConstants->numberOfFrequencies;j++)
	{
		int k = 0;
		index++;  /*  3*deltas[j], 0 */
		qffConstants->VI[j][k] = energies[index]; 

		index++; k++; /*  2*deltas[j], 0 */
		qffConstants->VI[j][k] = energies[index]; 

		index++; k++; /*  1*deltas[j], 0 */
		qffConstants->VI[j][k] = energies[index]; 

		index++; k++; /* -1*deltas[j], 0 */
		qffConstants->VI[j][k] = energies[index]; 

		index++; k++; /* -2*deltas[j], 0 */
		qffConstants->VI[j][k] = energies[index]; 
	
		index++; k++; /* -3*deltas[j], 0 */
		qffConstants->VI[j][k] = energies[index]; 
	}

	if(order>=2)
	for(j=0;j<qffConstants->numberOfFrequencies;j++)
	{
		for(i=0;i<j;i++)
		{
			int k;
			index++; /*  deltas[j],   deltas[i] */
                        k = 0;
			qffConstants->VIJ[j][i][k] = energies[index]; 

			index++;/*  deltas[j],  -deltas[i] */
                        k = 1; 
			qffConstants->VIJ[j][i][k] = energies[index]; 

			index++; /* -deltas[j],   deltas[i] */
                        k = 2; 
			qffConstants->VIJ[j][i][k] = energies[index]; 

			index++; /* -deltas[j],  -deltas[i] */
                        k = 3;
			qffConstants->VIJ[j][i][k] = energies[index]; 

                	qffConstants->VIJ[i][j][0] = qffConstants->VIJ[j][i][0];
                	qffConstants->VIJ[i][j][1] = qffConstants->VIJ[j][i][2];
                	qffConstants->VIJ[i][j][2] = qffConstants->VIJ[j][i][1];
                	qffConstants->VIJ[i][j][3] = qffConstants->VIJ[j][i][3];

		}
		for(i=0;i<qffConstants->numberOfFrequencies;i++)
		{
			int k;
			if(i==j) continue;
			index++; /*  deltas[j],  3*deltas[i] */
                        k = 0;
			qffConstants->VI3J[j][i][k] = energies[index]; 

			index++; /*  deltas[j], -3*deltas[i] */
                        k = 1; 
			qffConstants->VI3J[j][i][k] = energies[index]; 

			index++; /* -deltas[j],  3*deltas[i] */
                        k = 2; 
			qffConstants->VI3J[j][i][k] = energies[index]; 

			index++; /* -deltas[j], -3*deltas[i] */
                        k = 3;
			qffConstants->VI3J[j][i][k] = energies[index]; 
		}
	}
	if(order>=3)
	{
	FILE* fileVIJK = fopen("VIJK.bin","wb");
	for(j=0;j<qffConstants->numberOfFrequencies;j++)
	for(i=0;i<j;i++)
	for(k=0;k<i;k++)
	{
			int l = -1;
			index++; /*  deltas[j],   deltas[i],  deltas[k] */
			l++;
			fwrite(&energies[index], sizeof(double), 1, fileVIJK);
			//qffConstants->VIJK[j][i][k][l] = energies[index]; 

			index++; /*  deltas[j],   deltas[i],  -deltas[k] */
			l++;
			//qffConstants->VIJK[j][i][k][l] = energies[index]; 
			fwrite(&energies[index], sizeof(double), 1, fileVIJK);

			index++; /*  deltas[j],   -deltas[i],   deltas[k] */
			l++;
			//qffConstants->VIJK[j][i][k][l] = energies[index]; 
			fwrite(&energies[index], sizeof(double), 1, fileVIJK);

			index++; /*  deltas[j],   -deltas[i],   -deltas[k] */
			l++;
			//qffConstants->VIJK[j][i][k][l] = energies[index]; 
			fwrite(&energies[index], sizeof(double), 1, fileVIJK);

			index++; /*  -deltas[j],   deltas[i],  deltas[k] */
			l++;
			//qffConstants->VIJK[j][i][k][l] = energies[index]; 
			fwrite(&energies[index], sizeof(double), 1, fileVIJK);

			index++; /*  -deltas[j],   deltas[i],  -deltas[k] */
			l++;
			//qffConstants->VIJK[j][i][k][l] = energies[index]; 
			fwrite(&energies[index], sizeof(double), 1, fileVIJK);

			index++; /*  -deltas[j],   -deltas[i],   deltas[k] */
			l++;
			//qffConstants->VIJK[j][i][k][l] = energies[index]; 
			fwrite(&energies[index], sizeof(double), 1, fileVIJK);

			index++; /*  -deltas[j],   -deltas[i],   -deltas[k] */
			l++;
			//qffConstants->VIJK[j][i][k][l] = energies[index]; 
			fwrite(&energies[index], sizeof(double), 1, fileVIJK);

	}
	fclose(fileVIJK);
	}
	if(order>=4)
	for(j=0;j<qffConstants->numberOfFrequencies;j++)
	for(i=0;i<j;i++)
	for(k=0;k<i;k++)
	for(l=0;l<k;l++)
	{
			int n = 0;
			for(n=0;n<16;n++)
			{
				index++;
				qffConstants->VIJKL[j][i][k][l][n]= energies[index]; 
			}
	}

	free(energies);

	if(Ok) computeGradients(qffConstants);
	if(Ok) computeFrequencies(qffConstants);
	if(Ok) computeCubicForces(qffConstants);
	if(Ok) computeQuarticForces(qffConstants);
	//printf("Begin free\n");

	freeMatrixDouble(&qffConstants->VI, qffConstants->numberOfFrequencies);
	freeCubeDouble(&qffConstants->VIJ, qffConstants->numberOfFrequencies,qffConstants->numberOfFrequencies);
	freeCubeDouble(&qffConstants->VI3J, qffConstants->numberOfFrequencies,qffConstants->numberOfFrequencies);
	if(order>=3) freeQuarticDouble(&qffConstants->VIJK, qffConstants->numberOfFrequencies,qffConstants->numberOfFrequencies, qffConstants->numberOfFrequencies);
	if(order>=4) freeQuinticDouble(&qffConstants->VIJKL, qffConstants->numberOfFrequencies,qffConstants->numberOfFrequencies, qffConstants->numberOfFrequencies, qffConstants->numberOfFrequencies);
	//printf("End free\n");
}
/*****************************************************************************/
static void computeQFFDerivativesDipoles(NMRQFF* qffConstants, int order, double* dipoles[])
{
	int i;
	int j;
	int k;
	//int l;
	int index;
	int xyz;
	boolean Ok = TRUE;
	index = 0;

	initnMRQFFDipoles(qffConstants, order);

	for(xyz=0;xyz<3;xyz++) qffConstants->dipole0[xyz] = dipoles[xyz][index];


	for(j=0;j<qffConstants->numberOfFrequencies;j++)
	{
		int k = 0;
		index++;  /*  3*deltas[j], 0 */
		for(xyz=0;xyz<3;xyz++) qffConstants->dipolesI[j][k][xyz] = dipoles[xyz][index]; 

		index++; k++; /*  2*deltas[j], 0 */
		for(xyz=0;xyz<3;xyz++) qffConstants->dipolesI[j][k][xyz] = dipoles[xyz][index]; 

		index++; k++; /*  1*deltas[j], 0 */
		for(xyz=0;xyz<3;xyz++) qffConstants->dipolesI[j][k][xyz] = dipoles[xyz][index]; 

		index++; k++; /* -1*deltas[j], 0 */
		for(xyz=0;xyz<3;xyz++) qffConstants->dipolesI[j][k][xyz] = dipoles[xyz][index]; 

		index++; k++; /* -2*deltas[j], 0 */
		for(xyz=0;xyz<3;xyz++) qffConstants->dipolesI[j][k][xyz] = dipoles[xyz][index]; 
	
		index++; k++; /* -3*deltas[j], 0 */
		for(xyz=0;xyz<3;xyz++) qffConstants->dipolesI[j][k][xyz] = dipoles[xyz][index]; 
	}

	if(order>=2)
	for(j=0;j<qffConstants->numberOfFrequencies;j++)
	{
		for(i=0;i<j;i++)
		{
			int k;
			index++; /*  deltas[j],   deltas[i] */
                        k = 0;
			for(xyz=0;xyz<3;xyz++) qffConstants->dipolesIJ[j][i][k][xyz] = dipoles[xyz][index]; 

			index++;/*  deltas[j],  -deltas[i] */
                        k = 1; 
			for(xyz=0;xyz<3;xyz++) qffConstants->dipolesIJ[j][i][k][xyz] = dipoles[xyz][index]; 

			index++; /* -deltas[j],   deltas[i] */
                        k = 2; 
			for(xyz=0;xyz<3;xyz++) qffConstants->dipolesIJ[j][i][k][xyz] = dipoles[xyz][index]; 

			index++; /* -deltas[j],  -deltas[i] */
                        k = 3;
			for(xyz=0;xyz<3;xyz++) qffConstants->dipolesIJ[j][i][k][xyz] = dipoles[xyz][index]; 

			for(xyz=0;xyz<3;xyz++)
                	{
				qffConstants->dipolesIJ[i][j][0][xyz] = qffConstants->dipolesIJ[j][i][0][xyz];
				qffConstants->dipolesIJ[i][j][1][xyz] = qffConstants->dipolesIJ[j][i][2][xyz];
				qffConstants->dipolesIJ[i][j][2][xyz] = qffConstants->dipolesIJ[j][i][1][xyz];
				qffConstants->dipolesIJ[i][j][3][xyz] = qffConstants->dipolesIJ[j][i][3][xyz];
                	}

		}
		for(i=0;i<qffConstants->numberOfFrequencies;i++)
		{
			int k;
			if(i==j) continue;
			index++; /*  deltas[j],  3*deltas[i] */
                        k = 0;
			for(xyz=0;xyz<3;xyz++) qffConstants->dipolesI3J[j][i][k][xyz] = dipoles[xyz][index]; 

			index++; /*  deltas[j], -3*deltas[i] */
                        k = 1; 
			for(xyz=0;xyz<3;xyz++) qffConstants->dipolesI3J[j][i][k][xyz] = dipoles[xyz][index]; 

			index++; /* -deltas[j],  3*deltas[i] */
                        k = 2; 
			for(xyz=0;xyz<3;xyz++) qffConstants->dipolesI3J[j][i][k][xyz] = dipoles[xyz][index]; 

			index++; /* -deltas[j], -3*deltas[i] */
                        k = 3;
			for(xyz=0;xyz<3;xyz++) qffConstants->dipolesI3J[j][i][k][xyz] = dipoles[xyz][index]; 
		}
	}
	if(order>=3)
	{
	FILE* fileDIJK = fopen("DIJK.bin","wb");
	for(j=0;j<qffConstants->numberOfFrequencies;j++)
	{
		for(i=0;i<j;i++)
		for(k=0;k<i;k++)
		{
			int l = -1;
			index++; /*  deltas[j],   deltas[i],  deltas[k] */
			l++;
			//for(xyz=0;xyz<3;xyz++) qffConstants->dipolesIJK[j][i][k][l][xyz] = dipoles[xyz][index]; 
			for(xyz=0;xyz<3;xyz++) fwrite(&dipoles[xyz][index], sizeof(double), 1, fileDIJK);

			index++; /*  deltas[j],   deltas[i],  -deltas[k] */
			l++;
			//for(xyz=0;xyz<3;xyz++) qffConstants->dipolesIJK[j][i][k][l][xyz] = dipoles[xyz][index]; 
			for(xyz=0;xyz<3;xyz++) fwrite(&dipoles[xyz][index], sizeof(double), 1, fileDIJK);

			index++; /*  deltas[j],   -deltas[i],   deltas[k] */
			l++;
			//for(xyz=0;xyz<3;xyz++) qffConstants->dipolesIJK[j][i][k][l][xyz] = dipoles[xyz][index]; 
			for(xyz=0;xyz<3;xyz++) fwrite(&dipoles[xyz][index], sizeof(double), 1, fileDIJK);

			index++; /*  deltas[j],   -deltas[i],   -deltas[k] */
			l++;
			//for(xyz=0;xyz<3;xyz++) qffConstants->dipolesIJK[j][i][k][l][xyz] = dipoles[xyz][index]; 
			for(xyz=0;xyz<3;xyz++) fwrite(&dipoles[xyz][index], sizeof(double), 1, fileDIJK);

			index++; /*  -deltas[j],   deltas[i],  deltas[k] */
			l++;
			//for(xyz=0;xyz<3;xyz++) qffConstants->dipolesIJK[j][i][k][l][xyz] = dipoles[xyz][index]; 
			for(xyz=0;xyz<3;xyz++) fwrite(&dipoles[xyz][index], sizeof(double), 1, fileDIJK);

			index++; /*  -deltas[j],   deltas[i],  -deltas[k] */
			l++;
			//for(xyz=0;xyz<3;xyz++) qffConstants->dipolesIJK[j][i][k][l][xyz] = dipoles[xyz][index]; 
			for(xyz=0;xyz<3;xyz++) fwrite(&dipoles[xyz][index], sizeof(double), 1, fileDIJK);

			index++; /*  -deltas[j],   -deltas[i],   deltas[k] */
			l++;
			//for(xyz=0;xyz<3;xyz++) qffConstants->dipolesIJK[j][i][k][l][xyz] = dipoles[xyz][index]; 
			for(xyz=0;xyz<3;xyz++) fwrite(&dipoles[xyz][index], sizeof(double), 1, fileDIJK);

			index++; /*  -deltas[j],   -deltas[i],   -deltas[k] */
			l++;
			//for(xyz=0;xyz<3;xyz++) qffConstants->dipolesIJK[j][i][k][l][xyz] = dipoles[xyz][index]; 
			for(xyz=0;xyz<3;xyz++) fwrite(&dipoles[xyz][index], sizeof(double), 1, fileDIJK);

			/*
			for(xyz=0;xyz<3;xyz++)
                	{
				for(l=0;l<8;l++) qffConstants->dipolesIJK[i][j][k][l][xyz] = qffConstants->dipolesIJK[j][i][k][l][xyz];
				for(l=0;l<8;l++) qffConstants->dipolesIJK[i][k][j][l][xyz] = qffConstants->dipolesIJK[j][i][k][l][xyz];
				for(l=0;l<8;l++) qffConstants->dipolesIJK[j][k][i][l][xyz] = qffConstants->dipolesIJK[j][i][k][l][xyz];
				for(l=0;l<8;l++) qffConstants->dipolesIJK[k][i][j][l][xyz] = qffConstants->dipolesIJK[j][i][k][l][xyz];
				for(l=0;l<8;l++) qffConstants->dipolesIJK[k][j][i][l][xyz] = qffConstants->dipolesIJK[j][i][k][l][xyz];
                	}
			*/
		}
	}
	fclose(fileDIJK);
	}
	/*
	if(order>=4)
	for(j=0;j<qffConstants->numberOfFrequencies;j++)
	for(i=0;i<j;i++)
	for(k=0;k<i;k++)
	for(l=0;l<k;l++)
	{
			int n = 0;
			for(n=0;n<16;n++)
			{
				index++;
			}
	}
	*/
	for(xyz=0;xyz<3;xyz++) free(dipoles[xyz]);


	if(Ok) changeUnitInputFirstDerivativesDipoles(qffConstants);
	if(Ok) computeFirstDerivativesDipoles(qffConstants);
	if(Ok) computeSecondDerivativesDipoles(qffConstants);
	if(Ok) computeCubicDerivativesDipoles(qffConstants);
	if(Ok) computeQuarticDerivativesDipoles(qffConstants);

	freeVectorDouble(&qffConstants->dipole0);
	freeCubeDouble(&qffConstants->dipolesI, qffConstants->numberOfFrequencies,6);
	freeQuarticDouble(&qffConstants->dipolesIJ, qffConstants->numberOfFrequencies,qffConstants->numberOfFrequencies,4);
	freeQuarticDouble(&qffConstants->dipolesI3J, qffConstants->numberOfFrequencies,qffConstants->numberOfFrequencies,4);
	if(order>=3) freeQuinticDouble(&qffConstants->dipolesIJK, qffConstants->numberOfFrequencies,qffConstants->numberOfFrequencies,qffConstants->numberOfFrequencies, 8);
}
/*****************************************************************************/
char* computeQFFDerivatives(char* inputFileName, Molecule* mol, int order, double* deltas, int nEnergies, double* energies, double* dipoles[])
{
	int i;
	NMRQFF* qffConstants;
	char* fileNameOut = NULL;

	qffConstants = (NMRQFF*) malloc(sizeof(NMRQFF));
	qffConstants->numberOfFrequencies = 0;
	printf("nAtoms = %d\n", mol->nAtoms);

	qffConstants->numberOfFrequencies = mol->vibration.nModes;
	printf("numberOfFrequencies = %d\n", qffConstants->numberOfFrequencies);
	initnMRQFF0(qffConstants,qffConstants->numberOfFrequencies);
	printf("End initnMRQFF0\n");

	for(i=0;i<qffConstants->numberOfFrequencies;i++)
	{
		qffConstants->frequencies[i] = mol->vibration.modes[i].frequency;
		qffConstants->mass[i] = mol->vibration.modes[i].mass;
		qffConstants->delta[i] = deltas[i];
	}

        qffConstants->numberOfEnergies = nEnergies;
	computeQFFDerivativesEnergies(qffConstants, order, energies);
	printf("Energ computeQFFDerivativesEnergies\n");
	computeQFFDerivativesDipoles(qffConstants, order, dipoles);
	printf("Energ computeQFFDerivativesDipoles\n");

	fileNameOut = saveNMRQFFAppend(qffConstants, inputFileName, mol);
	
	return fileNameOut;
}
