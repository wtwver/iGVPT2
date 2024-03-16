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

/********************************************************************************/
/********************************************************************************/
/**************** compute N2P2QFF using analytical method ***********************/
/********************************************************************************/
#ifdef HIGH_DERIVATIVES
#include "../InterfaceLibN2P2/InterfaceCChemIC.h"
/*****************************************************************************/
// eqs 26-29 : Mackie J. Chem. Phys. 142, 244107 (2015)
static double** getMatrixPartialDerivative(Molecule* mol)
{
// Translation and rotation vector are not needed for trans from cartesian to normal modes
// Omega bar is given by mol->vibration.modes[i].vectors[j][k]
	double** Otild = NULL;
	int i,j,k;
	if(mol->vibration.nModes<1) return Otild;
	if(mol->nAtoms<1) return Otild;
	Otild = malloc(mol->vibration.nModes*sizeof(double*));
	for(i=0;i< mol->vibration.nModes;i++)
		Otild[i] = malloc(3*mol->nAtoms*sizeof(double));

	for(i=0;i< mol->vibration.nModes;i++)
	{
		double mum12 = 1.0/sqrt(mol->vibration.modes[i].mass*AMUTOAU);
		for(k=0;k<mol->nAtoms;k++)  
		for(j=0;j<3;j++) 
			Otild[i][3*k+j] = mol->vibration.modes[i].vectors[j][k]*mum12;
	}
	return Otild;

}
/**************************************************************************************/
static FILE* newQFFAppendPre(char* inputFileName, Molecule* mol)
{

	char tmp[BSIZE];
	int nf = mol->vibration.nModes;
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
		"# QFF constants calculated using analytical method implemented in modified version of N2P2 library\n"
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

	return file;

}
static void QFFN2P2AppendEnergyDerivatives(FILE* file, Molecule* mol, int order, int method, char* HNNDir, double **Otild, int fourIndex)
{
	char tmp[BSIZE];
	double cfenergy = 1.0/(AUTOKCAL);
        double cflength = ANGTOBOHR;
        //double cfdipole = 1.0/(AUTODEB);

	double f1cchemiAu = cfenergy/cflength;
	double f2cchemiAu = f1cchemiAu/cflength;
	double f3cchemiAu = f2cchemiAu/cflength;
	double f4cchemiAu = f3cchemiAu/cflength;

	void* interfaceLibN2P2 = newInterfaceCChemI(HNNDir, cflength, cfenergy, 0);

	DerivativesIC* deriv;

	double f2cm1 = AUTOCM1;
	double f3cm1 = AUTOCM1*sqrt(AUTOCM1)*AUTOCM1;
	double f4cm1 = AUTOCM1*AUTOCM1*AUTOCM1;

	double f1conv = f1cchemiAu;
	double f2conv = f2cm1*f2cchemiAu;
	double f3conv = f3cm1*f3cchemiAu;
	double f4conv = f4cm1*f4cchemiAu;

	int nAtoms = mol->nAtoms;
	int i,j,k,l;
	int a,b,c,d;
	int ca,cb,cc,cd;
	int nf = mol->vibration.nModes;
	
	// Derivatives (in atomic unit) in cartezian coordinates
	deriv = interfaceCChemIComputeHighDerivatives(interfaceLibN2P2, mol, order, method);

	// eqs 26-29 : Mackie J. Chem. Phys. 142, 244107 (2015)
	// (3*nAtoms)(3*nAtoms) matrix

	// Derivatives in normal coordinates
	// Gradients
	if(order<2) return;

	sprintf(tmp,"#i Freq(cm-1)  Calc.Freq   dQ(Bohr)  Mass(amu)\tGradient[ H amu^(-1/2) Bohr^(-1)]\n");
	fprintf(file,"%s",tmp);   
	fprintf(file,"%s","Hessian\n");   
	for(i=0;i<nf;i++) 
	{
		double grad = 0 , f= 0, delta = 0;
		// Gradient
		grad = 0;
        	for(a=0;a<nAtoms;a++) for(ca=0;ca<3;ca++) 
		{
			int ia = 3*a+ca;
			grad +=  Otild[i][ia]*deriv->df[ia];
		}
		grad *= f1conv;
		// Quadratic
		f = 0;
        	for(a=0;a<nAtoms;a++) for(ca=0;ca<3;ca++) 
        	for(b=0;b<nAtoms;b++) for(cb=0;cb<3;cb++) 
		{
			int ia = 3*a+ca;
			int ib = 3*b+cb;
			double o = Otild[i][ia]*Otild[i][ib];
			if(ib>ia) { int t = ia; ia=ib; ib = t; }
			f += o*deriv->d2f[ia][ib];
		}
		f = f2conv*sqrt(fabs(f));
		sprintf(tmp,"%d %d %0.14f %0.14f %0.14f %0.14f\t%0.14f\n",i+1, i+1, mol->vibration.modes[i].frequency, 
				f, 
				delta, mol->vibration.modes[i].mass, grad);
		fprintf(file,"%s",tmp);   
	}
	fprintf(file,"%s","END\n\n");
	// Cubic
	if(order<3) return;
	sprintf(tmp,"# i\tj\tk\tReduced values [cm-1]\n");
	fprintf(file,"%s",tmp);   
	fprintf(file,"%s","Cubic\n");   
	for(i=0;i<nf;i++)
        {
                for(j=0;j<=i;j++)
                for(k=0;k<=j;k++)
                {
			double val = 0;
        		for(a=0;a<nAtoms;a++) for(ca=0;ca<3;ca++) 
        		for(b=0;b<nAtoms;b++) for(cb=0;cb<3;cb++) 
        		for(c=0;c<nAtoms;c++) for(cc=0;cc<3;cc++) 
			{
				int ia = 3*a+ca;
				int ib = 3*b+cb;
				int ic = 3*c+cc;
				double o = Otild[i][ia]*Otild[j][ib]*Otild[k][ic];
				if(ic>ib) { int t = ic; ic=ib; ib = t; }
				if(ib>ia) { int t = ib; ib=ia; ia = t; }
				if(ic>ib) { int t = ic; ic=ib; ib = t; }
				val += o*deriv->d3f[ia][ib][ic];
			}
			val *= f3conv/sqrt(mol->vibration.modes[i].frequency*mol->vibration.modes[j].frequency*mol->vibration.modes[k].frequency);
			if(fabs(val)<1e-12) continue;
			sprintf(tmp,"%d\t%d\t%d\t%14.6f\n",k+1, j+1, i+1, val);
			fprintf(file,"%s",tmp);   
                }
        }
	fprintf(file,"%s","END\n\n");
	//Quartic
	if(order<4) return;
	sprintf(tmp,"# i\tj\tk\tl\tReduced values [cm-1]\n");
	fprintf(file,"%s",tmp);   
	fprintf(file,"%s","Quartic\n");   
	for(i=0;i<nf;i++)
        {
                for(j=0;j<=i;j++)
                for(k=0;k<=j;k++)
                for(l=0;l<=k;l++)
		{
			if(!fourIndex && i!=j && i!=k && i!=l && j!=k && j!=l && k!=l) continue;
			double val = 0;
        		for(a=0;a<nAtoms;a++) for(ca=0;ca<3;ca++) 
        		for(b=0;b<nAtoms;b++) for(cb=0;cb<3;cb++) 
        		for(c=0;c<nAtoms;c++) for(cc=0;cc<3;cc++) 
        		for(d=0;d<nAtoms;d++) for(cd=0;cd<3;cd++) 
			{
				int ia = 3*a+ca;
				int ib = 3*b+cb;
				int ic = 3*c+cc;
				int id = 3*d+cd;
				double o = Otild[i][ia]*Otild[j][ib]*Otild[k][ic]*Otild[l][id];
				if(id>ic) { int t = id; id=ic; ic = t; }
				if(ic>ib) { int t = ic; ic=ib; ib = t; }
				if(ib>ia) { int t = ib; ib=ia; ia = t; }
				if(id>ic) { int t = id; id=ic; ic = t; }
				if(ic>ib) { int t = ic; ic=ib; ib = t; }
				if(id>ic) { int t = id; id=ic; ic = t; }
				val += o*deriv->d4f[ia][ib][ic][id];
			}
			val *= f4conv
			/sqrt(mol->vibration.modes[i].frequency*mol->vibration.modes[j].frequency*mol->vibration.modes[k].frequency*mol->vibration.modes[l].frequency);
			if(fabs(val)<1e-12) continue;
			sprintf(tmp,"%d\t%d\t%d\t%d\t%14.6f\n",l+1,k+1, j+1, i+1, val);
			fprintf(file,"%s",tmp);   
		}
        }
	fprintf(file,"%s","END\n\n");

}
/************************************************************************************************************/
static int readFirstDipolesInput(FILE* file, int nf, double** firstDipolesInput)
{
	int xyz;
	int nn=1;
	int i;
        double mu0 = 4*PI*1e-7;
        double eps0 = 1.0/(mu0*slight*slight);
        double   kmmolm1 = 4*PI*PI*PI*NAvogadro/3/hPlank/slight/4/PI/eps0*1e-3*100.0*8.47835267e-30*8.47835267e-30;/* 1e-3 m to km, 100 : cm-1 to m-1 */
	double f = 1.0/sqrt(kmmolm1);
	int numberOfFirstDipolesInput = 0;
	{
		rewind(file);
		if(goToStr(file, "First derivatives"))
		for(i=0;i<nf && nn==1 ;i++)
		{
			int xyz;
			for(xyz=0;xyz<3 && nn==1 ;xyz++)
				nn = fscanf(file,"%lf",&firstDipolesInput[i][xyz]);
			numberOfFirstDipolesInput += nn;
		}
		if(numberOfFirstDipolesInput != nf) numberOfFirstDipolesInput=0;
	}
	if(numberOfFirstDipolesInput != 0)
	for(i=0;i<nf;i++)
		for(xyz=0;xyz<3;xyz++)
			firstDipolesInput[i][xyz] *= f;

	return numberOfFirstDipolesInput;
}
/******************************************************************************************************************************/
static void QFFN2P2AppendDipoleDerivatives(char* inputFileName, FILE* file, Molecule* mol, int order, int method, char* HNNDir, double** Otild, int threeIndex)
{
	char tmp[BSIZE];

        double cflength = ANGTOBOHR;
        double cfdipole = 1.0/(AUTODEB);

	DerivativesIC** deriv = NULL;
	void* interfaceLibN2P2ES = NULL;
	double f1cm1 = sqrt(AUTOCM1);
	double f2cm1 = AUTOCM1;
	double f3cm1 = AUTOCM1*sqrt(AUTOCM1);

	double f1cchemiAu = cfdipole/cflength;
	double f2cchemiAu = f1cchemiAu/cflength;
	double f3cchemiAu = f2cchemiAu/cflength;

	double f1conv = f1cm1*f1cchemiAu;
	double f2conv = f2cm1*f2cchemiAu;
	double f3conv = f3cm1*f3cchemiAu;

	int nAtoms = mol->nAtoms;
	int i,j,k;
	int a,b,c;
	int ca,cb,cc;
	int nf = mol->vibration.nModes;
	int xyz;
	int numberOfFirstDipolesInput= 0;
	double **firstDipolesInput = NULL;
        FILE* fileIn = fopen(inputFileName,"rb");
	char txyz[]={'X','Y','Z'};

	if(nf>0)
	{
		firstDipolesInput = malloc(nf*sizeof(double*));
		for(i=0;i<nf;i++)
			firstDipolesInput[i] = malloc(3*sizeof(double));
		numberOfFirstDipolesInput =readFirstDipolesInput(fileIn, nf, firstDipolesInput);
	}

        interfaceLibN2P2ES = newInterfaceCChemIES(HNNDir, cflength, cfdipole, 0);
	
	// Derivatives (in atomic unit) in cartezian coordinates
	deriv = interfaceCChemIComputeDipoleHighDerivatives(interfaceLibN2P2ES, mol, order, method);

	// eqs 26-29 : Mackie J. Chem. Phys. 142, 244107 (2015)
	// (3*nAtoms)(3*nAtoms) matrix
	// Derivatives in normal coordinates
	// Gradients
	if(order<1) return;
	if(mol)
	{
		mol->klass->addFirstDerivativeToFile(mol, file);
	}
	if(numberOfFirstDipolesInput==0)
	{
		sprintf(tmp,"#xyz\ti\tValues[au cm^1/2]\n");
		fprintf(file,"%s",tmp);   
		fprintf(file,"%s","First derivatives\n");
		for(i=0;i<nf;i++)
		for(xyz=0;xyz<3;xyz++)
		{
			double val =  0;
        		for(a=0;a<nAtoms;a++) for(ca=0;ca<3;ca++) 
			{
				int ia = 3*a+ca;
				val +=  Otild[i][ia]*deriv[xyz]->df[ia];
			}
			val *=  f1conv;
			if(fabs(val)<1e-12) continue;
			sprintf(tmp,"%c\t%d\t%14.6f\n",txyz[xyz], i+1,val);
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
			double val =  0;
        		for(a=0;a<nAtoms;a++) for(ca=0;ca<3;ca++) 
			{
				int ia = 3*a+ca;
				val +=  Otild[i][ia]*deriv[xyz]->df[ia];
			}
			val *=  f1conv;
			if(fabs(firstDipolesInput[i][xyz])<1e-12) continue;
			sprintf(tmp,"%c\t%d\t%14.6f\t\t%14.6f\n",txyz[xyz], i+1,firstDipolesInput[i][xyz], val);
			fprintf(file,"%s",tmp);   
		}
		fprintf(file,"%s","END\n\n");
	}
	if( firstDipolesInput)
	{
		for(i=0;i<nf;i++) free(firstDipolesInput[i]);
		free(firstDipolesInput);
	}

	// Second derivatives
	if(order<2) return;
	sprintf(tmp,"#xyz\ti\tj\tValues[au cm]\n");
	fprintf(file,"%s",tmp);   
	fprintf(file,"%s","Second derivatives\n");
	for(i=0;i<nf;i++)
	for(j=0;j<nf;j++)
	for(xyz=0;xyz<3;xyz++)
	{
		double val = 0;
        	for(a=0;a<nAtoms;a++) for(ca=0;ca<3;ca++) 
        	for(b=0;b<nAtoms;b++) for(cb=0;cb<3;cb++) 
		{
			int ia = 3*a+ca;
			int ib = 3*b+cb;
			double  o = Otild[i][ia]*Otild[j][ib];
			if(ib>ia) { int t = ia; ia=ib; ib = t; }
			val += o*deriv[xyz]->d2f[ia][ib];
		}
		val *= f2conv;
		if(fabs(val)<1e-12) continue;
		sprintf(tmp,"%c\t%d\t%d\t%14.6f\n",txyz[xyz], i+1,j+1,val);
		fprintf(file,"%s",tmp);   
	}
	fprintf(file,"%s","END\n\n");


	// Cubic derivatives
	if(order<3) return;
	sprintf(tmp,"#xyz\ti\tj\tk\tValues[au cm^3/2]\n");
	fprintf(file,"%s",tmp);   
	fprintf(file,"%s","Cubic derivatives\n");
	for(i=0;i<nf;i++)
	for(j=0;j<nf;j++)
	for(k=0;k<nf;k++)
	for(xyz=0;xyz<3;xyz++)
        {
		double val = 0;
		if(!threeIndex && i!=j && i!=k && j!=k) continue;

        	for(a=0;a<nAtoms;a++) for(ca=0;ca<3;ca++) 
        	for(b=0;b<nAtoms;b++) for(cb=0;cb<3;cb++) 
        	for(c=0;c<nAtoms;c++) for(cc=0;cc<3;cc++) 
		{
			int ia = 3*a+ca;
			int ib = 3*b+cb;
			int ic = 3*c+cc;
			double o = Otild[i][ia]*Otild[j][ib]*Otild[k][ic];
			if(ic>ib) { int t = ic; ic=ib; ib = t; }
			if(ib>ia) { int t = ib; ib=ia; ia = t; }
			if(ic>ib) { int t = ic; ic=ib; ib = t; }
			val += o*deriv[xyz]->d3f[ia][ib][ic];
		}
		val *= f3conv;
		if(fabs(val)<1e-12) continue;
		sprintf(tmp,"%c\t%d\t%d\t%d\t%14.6f\n",txyz[xyz], i+1,j+1,k+1,val);
		fprintf(file,"%s",tmp);
        }
	fprintf(file,"%s","END\n\n");
}
char* computeN2P2QFFAnalyticDerivatives(char* inputFileName, Molecule* mol, int orderEnergy, int orderDipole, int method, char* HNNDir)
{
	char* fileNameOut = NULL;
	double **Otild = NULL;
	FILE* file = NULL;
	int fourIndex = 0;
	int threeIndex = 0;

	
	if(orderEnergy>=5) 
	{
		orderEnergy = 4;
		fourIndex = 1;
	}
	if(orderDipole>=4) 
	{
		orderDipole = 3;
		threeIndex = 1;
	}
	printf("numberOfFrequencies = %d\n", mol->vibration.nModes);
	Otild = getMatrixPartialDerivative(mol);

	printf("newQFFAppendPre\n");
	file =  newQFFAppendPre(inputFileName, mol);
	if(file) fileNameOut = strdup_printf("%sQFF.txt",getSuffixNameFile(inputFileName));
	printf("QFFN2P2AppendEnergyDerivatives\n");
	QFFN2P2AppendEnergyDerivatives(file, mol, orderEnergy, method, HNNDir, Otild, fourIndex);
	printf("QFFN2P2AppendDipoleDerivatives\n");
	QFFN2P2AppendDipoleDerivatives(inputFileName, file, mol, orderDipole, method, HNNDir, Otild, threeIndex);
	fclose(file);
	return fileNameOut;
}
/**********************************************************************************************************************/
#endif
