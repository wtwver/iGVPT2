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

/* VPT2Model.c */
#include <math.h>
#include "VPT2Model.h"

static boolean printMax = FALSE;

static void computeAnharmonic(VPT2Model* vpt2Model);
static void readData(VPT2Model* vpt2Model, char* inputFileName);
static boolean testDegen(VPT2Model* vpt2Model, double w1, double w2);
static boolean testResonanceCubic(VPT2Model* vpt2Model, double w1, double w2, double v);
static boolean testResonanceQuartic(VPT2Model* vpt2Model, double w1, double w2, double v);
static double applyModelPropNew(VPT2Model* vpt2Model, double k2, int ii, int jj, int kk);
/**********************************************************************/
static boolean testFermi(VPT2Model* vpt2Model, int ii, int jj, int kk)
{
	VPT2PotentialData* data = &vpt2Model->vData;
        double maxFrequencyDifferenceFermi = data->maxFrequencyDifferenceFermi;
        double XI = data->parametersResonance[0];
        double XII = data->parametersResonance[1];
        double Z = data->parametersResonance[2];
	int i=abs(ii)-1;
	int j=abs(jj)-1;
	int k=abs(kk)-1;
	double wi;
	double wj;
	double wk;
	double sqrt8=sqrt(1.0/8.0);
	boolean res = FALSE;
	double sijk = 0;

	if(Z<=0 && XI<=0 && XII<=0) 
	{
		XI = 1.0;
		XII = 1.0;
		if(maxFrequencyDifferenceFermi<=0) maxFrequencyDifferenceFermi = 200;
	}

	wi=(ii>0)?data->hessian[i][i]:-data->hessian[i][i];
	wj=(jj>0)?data->hessian[j][j]:-data->hessian[j][j];
	wk=(kk>0)?data->hessian[k][k]:-data->hessian[k][k];
	sijk =wi+wj+wk;
	if(fabs(sijk)<maxFrequencyDifferenceFermi) 
	{
		if(Z>0)
		{
			double fijk=getMaxCubicIJK(vpt2Model->vData.cubic,i,j,k)*sqrt8;
			if(i==j||i==k||j==k) fijk *= 0.5;
			if(fijk/fabs(sijk)>Z) 
			{
					res = TRUE;
					fprintf(stderr,"Warning Fermi resonnance : ");
					fprintf(stderr,"k,l,m = %d %d %d dklm=%f Z=%f Zcut=%f\n",i+1,j+1,k+1,sijk,fijk/fabs(sijk),Z);
			}
		}
		else if(XI>0 && XII>0)
		{
			double X = 0;
			double dijk = fabs(sijk);
			double fijk=getMaxCubicIJK(vpt2Model->vData.cubic,i,j,k);
			X = (fijk*fijk*fijk*fijk)/(dijk*dijk*dijk)/64.0;
			if(i==j||i==k||j==k)
			{
				X *= 0.25;
				if(X>XI) 
				{
					res = TRUE;
					fprintf(stderr,"Warning Fermi resonnance : ");
					fprintf(stderr,"k,l,m = %d %d %d dklm=%f X=%f XI=%f\n",ii,jj,kk,sijk,X,XI);
				}
			}
			else if(X>XII) { 
				res = TRUE;
				fprintf(stderr,"Warning Fermi resonnance : ");
				fprintf(stderr,"k,l,m = %d %d %d dklm=%f X=%f XI=%f\n",ii,jj,kk,sijk,X,XI);
			}
		}
	}
	return res;
}
/**********************************************************************/
static void computeFundamentalsEnergies(VAnharmonic* vAnharmonic)
{
	int i,j;
	if(!vAnharmonic) return;
	if(!vAnharmonic->Xi) return;
	for(i=0;i<vAnharmonic->nFrequencies;i++) 
	{
		vAnharmonic->fundamentals[i] = vAnharmonic->harmonicFrequencies[i]+2*vAnharmonic->Xi[i][i];
		for(j=0;j<vAnharmonic->nFrequencies;j++)
		{
			if(i!=j) vAnharmonic->fundamentals[i] += 0.5*vAnharmonic->Xi[i][j];
		}
	}
}
/**********************************************************************/
static void printFundamentalsEnergies(VAnharmonic* vAnharmonic)
{
	int i;
	if(!vAnharmonic) return;
	if(!vAnharmonic->fundamentals) return;
	printf("Fundamentals\n");
	printf("------------\n");
	printf("%-10s %-10s %-20s %-20s\n","Mode","Quanta","E(harm)","E(anharm)");
	for(i=0;i<vAnharmonic->nFrequencies;i++) 
		printf("%-10d %-10d %-20.8f %-20.8f\n",i+1,1,vAnharmonic->harmonicFrequencies[i], vAnharmonic->fundamentals[i]);
	printf("\n");
}
/**********************************************************************/
static void computeOvertonesEnergies(VAnharmonic* vAnharmonic)
{
	int i,j;
	if(!vAnharmonic) return;
	if(!vAnharmonic->Xi) return;
	for(i=0;i<vAnharmonic->nFrequencies;i++) 
	{
		vAnharmonic->overtones[i] = 2*vAnharmonic->harmonicFrequencies[i]+6*vAnharmonic->Xi[i][i];
		for(j=0;j<vAnharmonic->nFrequencies;j++)
			if(i!=j) vAnharmonic->overtones[i] += vAnharmonic->Xi[i][j];
	}
}
/**********************************************************************/
static void printOvertonesEnergies(VAnharmonic* vAnharmonic)
{
	int i;
	if(!vAnharmonic) return;
	if(!vAnharmonic->Xi) return;
	printf("Overtones\n");
	printf("---------\n");
	printf("%-10s %-10s %-20s %-20s\n","Mode","Quanta","E(harm)","E(anharm)");
	for(i=0;i<vAnharmonic->nFrequencies;i++) 
		printf("%-10d %-10d %-20.8f %-20.8f\n",i+1,2,2*vAnharmonic->harmonicFrequencies[i], vAnharmonic->overtones[i]);
	printf("\n");
}
/**********************************************************************/
static void computeCombinationBandsEnergies(VAnharmonic* vAnharmonic)
{
	int i,j,k;
	if(!vAnharmonic) return;
	if(!vAnharmonic->Xi) return;

	for(i=0;i<vAnharmonic->nFrequencies;i++) 
	for(j=0;j<i;j++) 
	{
		vAnharmonic->combinationBands[i][j] = 
						  vAnharmonic->harmonicFrequencies[i]
						+ vAnharmonic->harmonicFrequencies[j]
						+ 2*vAnharmonic->Xi[i][i]
						+ 2*vAnharmonic->Xi[j][j]
						+ 2*vAnharmonic->Xi[i][j]
						;
		for(k=0;k<vAnharmonic->nFrequencies;k++)
			if(k!=i && k!=j) vAnharmonic->combinationBands[i][j] += 0.5*vAnharmonic->Xi[i][k]+ 0.5*vAnharmonic->Xi[j][k];
	}
}
/**********************************************************************/
static void printCombinationBandsEnergies(VAnharmonic* vAnharmonic)
{
	int i,j;
	if(!vAnharmonic) return;
	if(!vAnharmonic->Xi) return;
	printf("Combination Bands\n");
	printf("-----------------\n");
	printf("%-10s %-10s %-10s %-10s %-20s %-20s\n","Mode","Quanta","Mode","Quanta","E(harm)","E(anharm)");
	for(i=0;i<vAnharmonic->nFrequencies;i++) 
	for(j=0;j<i;j++) 
		printf("%-10d %-10d %-10d %-10d %-20.8f %-20.8f\n",i+1,1,j+1,1,vAnharmonic->harmonicFrequencies[i]+vAnharmonic->harmonicFrequencies[j], vAnharmonic->combinationBands[i][j]);
	printf("\n");
}
/**********************************************************************/
static void printIREnergies(VAnharmonic* vAnharmonic)
{
	printFundamentalsEnergies(vAnharmonic);
	printOvertonesEnergies(vAnharmonic);
	printCombinationBandsEnergies(vAnharmonic);
}
/**********************************************************************/
static double VPT2(double k2, double epsilon)
{
	double vpt2 = 0.0;
	if(fabs(epsilon)>1e-10) vpt2 = k2/epsilon;
	//else printf("warning k2 = %f epslison = %f\n",k2,epsilon);
	return vpt2;
}
/**********************************************************************/
static double VPT2Fermi(VPT2Model* vpt2Model, double k2, int ii, int jj, int kk)
{
	VPT2PotentialData* vData = &vpt2Model->vData;
	double wi;
	double wj;
	double wk;
	double epsilon;
	if(testFermi(vpt2Model, ii, jj, kk)) return 0.0;

	wi=(ii>0)?vData->hessian[ii-1][ii-1]:-vData->hessian[-ii-1][-ii-1];
	wj=(jj>0)?vData->hessian[jj-1][jj-1]:-vData->hessian[-jj-1][-jj-1];
	wk=(kk>0)?vData->hessian[kk-1][kk-1]:-vData->hessian[-kk-1][-kk-1];
	epsilon = wi+wj+wk;
	return k2/epsilon;
}
/**********************************************************************/
/* PT2 pModel: Degeneracy-corrected PT2 (DCPT2)
 Refs: K.M. Kuhler, D.G. Truhlar, A.D. Isaacson,
       J. Chem. Phys. 104, 12, 4664 (1996)
       &
       J. Bloino, M. Biczysko and V. Barone, JCTC, 8, 1015 (2012)
*/
static double DCPT2(double k2, double epsilon)
{
	double sign = (epsilon<=0)?-1.0:1.0;
	double sign2 = (k2<=0)?-1.0:1.0;
	double e = sign*epsilon*0.5;
	double r =  sign*sign2*(sqrt(sign2*k2+e*e)-e);
	//double ex = k2/epsilon;
//	printf("k2=%e\n",k2);
//	printf("e2=%e\n",e*e);
//	printf("r=%e r0=%e diff=%e\n",r,ex,r-ex);
	k2*=1;
	e *=1;
	r =  sign*sign2*(sqrt(sign2*k2+e*e)-e);
//	printf("r=%e r0=%e diff=%e\n",r,ex,r-ex);
	return r;
}
/**********************************************************************/
static double lambda(double k2, double epsilon, double alpha, double beta)
{
	double e2 = epsilon*epsilon/4;
	double x = sqrt(fabs(k2*e2))-beta;
	double l = (tanh(alpha*x)+1.0)/2.0;
//	printf("k2=%e\n",k2);
//	printf("e2=%e\n",e2);
//	printf("k2*e2=%e\n",k2*e2);
//	printf("sqrt(k2*e2)-beta=%e\n",sqrt(fabs(k2*e2))-beta);

//	printf("l%e\n",l);
	return l;
}
/**********************************************************************/
static double HDCPT2(double k2, double epsilon, double alpha, double beta)
{
	double l = lambda(k2,epsilon,alpha,beta);
	double dcpt2 = DCPT2(k2,epsilon);
	double vpt2 = VPT2(k2,epsilon);
//	printf("lambda = %f\n",l);
	return (1-l)*dcpt2+l*vpt2;
}
/**********************************************************************/
static double applyModelEnergies(double k2, double epsilon, VPT2VModel* model)
{
	if(model->type==MODEL_DCPT2) return DCPT2(k2,epsilon);
	if(model->type==MODEL_HDCPT2) return HDCPT2(k2,epsilon,model->alphaHDCPT2,model->betaHDCPT2);
	return VPT2(k2,epsilon);
}
/**********************************************************************/
/* PCCP, 2014, 16, 1759-1787, page 1761 */
static void computeXi(VPT2Model* vpt2Model)
{
	int i,j,k;
	int c;
	double s = 0;
	VAnharmonic* vAnharmonic = &vpt2Model->vAnharmonic;
	VPT2PotentialData* data = &vpt2Model->vData;
	VPT2VModel* model = NULL;
	if(!vAnharmonic) return;
	if(!data) return;
	if(!data->hessian) return;
	if(!vAnharmonic->Xi) return;
	model = &data->model;
	/* compute Xii */
	for(i=0;i<vAnharmonic->nFrequencies;i++)
	{
		double wi=data->hessian[i][i];
		double si = 1/(wi);
		double s16w2 = si*si/16;
		double sqrti = (wi>0)?sqrt(wi):0;
		double kiii=data->cubic[i][i][i];
		double Sii = 0;
		vAnharmonic->Xi[i][i] = 0;
		/* Nii terms  see JCTC 2012, 8, 1015*/
		vAnharmonic->Xi[i][i] += data->quartic[i][i][i][i]/16.0;
		vAnharmonic->Xi[i][i] += -kiii*kiii/wi*5.0/3.0/16.0;
		/* Sii terms */
		Sii = 0;
		for(j=0;j<vAnharmonic->nFrequencies;j++)
		{
			double kiij=getMaxCubicIJK(data->cubic,i,i,j);
			double wj=data->hessian[j][j];
			double sqrtj = (wj>0)?sqrt(wj):0;
			double sj = 1/(wj);
			double sj2 = sj*sj;
			double siij = 1/(wi+wi+wj);
			double Kiij = kiij*sqrti*sqrti*sqrtj;
			double Kiij2 = Kiij*Kiij;
			if(j==i) continue;
		
			Sii += -(2*Kiij2*sj2+Kiij2*0.5*sj*siij)*s16w2;
			Sii += applyModelEnergies(Kiij2*0.5*sj*s16w2,2*wi-wj,model);
		}
		vAnharmonic->Xi[i][i] += Sii;
	}
	/* compute Xij */
	/* Nij terms */
	for(i=0;i<vAnharmonic->nFrequencies;i++)
	{
		double wi=data->hessian[i][i];
		double wi2=wi*wi;
		double kiii=data->cubic[i][i][i];
		for(j=0;j<vAnharmonic->nFrequencies;j++)
		{
			double wj=data->hessian[j][j];
			double wj2=wj*wj;
			double kjjj=getMaxCubicIJK(data->cubic,j,j,j);
			double kiij=getMaxCubicIJK(data->cubic,i,i,j);
			double kijj=getMaxCubicIJK(data->cubic,i,j,j);
			double kiijj=getMaxQuarticIJKL(data->quartic, i, i, j, j);
			double num;
			double denom;
			if(i==j) continue;

			vAnharmonic->Xi[i][j] = 0;
			vAnharmonic->Xi[i][j] += kiijj;

			num = kiii*kijj;
			denom=wi;
			if(fabs(denom)>1e-10) vAnharmonic->Xi[i][j] += -num/denom;

			num = kjjj*kiij;
			denom=wj;
			if(fabs(denom)>1e-10) vAnharmonic->Xi[i][j] += -num/denom;

			s = 0;
			for(k=0;k<vAnharmonic->nFrequencies;k++)
			{
				double wk=data->hessian[k][k];
				double kiik=getMaxCubicIJK(data->cubic,i,i,k);
				double kjjk=getMaxCubicIJK(data->cubic,j,j,k);
				if(k==j || k==i) continue;
				num = kiik*kjjk;
				denom=wk;
				if(fabs(denom)>1e-10) s += -num/denom;
			}
			vAnharmonic->Xi[i][j] += s;
			/* coriolis contributions */
			num = 4*(wi2+wj2);
			denom=wi*wj;
			if(fabs(denom)>1e-10) 
			{
				double s = 0;
				for(c=0;c<3;c++)
				{
					double xetaij=data->coriolis[c][i][j];
					if(fabs(data->coriolis[c][j][i])>fabs(xetaij)) xetaij=data->coriolis[c][j][i];
					s += data->Be[c]*xetaij*xetaij;
				}
				vAnharmonic->Xi[i][j] += num/denom*s;
			}
		}
	}
	for(i=0;i<vAnharmonic->nFrequencies;i++)
		for(j=0;j<vAnharmonic->nFrequencies;j++)
			if(i!=j)vAnharmonic->Xi[i][j] /= 4;
	/* Sij terms */
	/* JCTC, Bloino, 2012 : page 1018 : error in sign, first term in Sij */
	for(i=0;i<vAnharmonic->nFrequencies;i++)
	{
		double wi=data->hessian[i][i];
		double si=1/wi;
		double sqrti = sqrt(wi);

		for(j=0;j<vAnharmonic->nFrequencies;j++)
		{
			double wj=data->hessian[j][j];
			double sqrtj = sqrt(wj);
			double sj=1/wj;
			double s4wiwj=si*sj/4;
			double siij=1/(wi+wi+wj);
			double sjji=1/(wj+wj+wi);
			double kiij=getMaxCubicIJK(data->cubic,i,i,j);
			double kijj=getMaxCubicIJK(data->cubic,i,j,j);
			double Kiij=kiij*sqrti*sqrti*sqrtj;
			double Kijj=kijj*sqrti*sqrtj*sqrtj;
			double Kiij2=Kiij*Kiij;
			double Kijj2=Kijj*Kijj;
			double Sij = 0;
			if(i==j) continue;

			Sij += -Kiij2/8/wi/wi/wj*siij;
			Sij += -applyModelEnergies(Kiij2/8/wi/wi/wj,2*wi-wj,model);
			Sij += -Kijj2/8/wi/wj/wj*sjji;
			Sij += -applyModelEnergies(Kijj2/8/wi/wj/wj,2*wj-wi,model);

			for(k=0;k<vAnharmonic->nFrequencies;k++)
			{
				double wk=data->hessian[k][k];
				double sqrtk = sqrt(wk);
				double sk=1/wk;
				double sijk=1/(wi+wj+wk);
				double kijk=getMaxCubicIJK(data->cubic,i,j,k);
				double Kijk=kijk*sqrti*sqrtj*sqrtk;
				double Kijk2=Kijk*Kijk;
				double k2 = Kijk2*0.5*sk*s4wiwj;
				if(k==j || k==i) continue;
				Sij +=  -k2*sijk;
				/*
				Sij +=-k2/(wk-wi-wj);
				Sij += k2/(wi-wj-wk);
				Sij += k2/(wj-wi-wk);
				*/
				Sij +=  -applyModelEnergies(k2,wk-wi-wj,model);
				Sij += applyModelEnergies(k2,wi-wj-wk,model);
				Sij += applyModelEnergies(k2,wj-wi-wk,model);
			}
			vAnharmonic->Xi[i][j] += Sij;
		}
	}
}
/**********************************************************************/
static VAnharmonic newVAnharmonic(VPT2PotentialData* vData)
{
	VAnharmonic vAnharmonic;
	int n = 0;
	if(vData) n = vData->nFrequencies;
	vAnharmonic.nFrequencies = n;
	if(vAnharmonic.nFrequencies<=0) vAnharmonic.nFrequencies = 0;

	vAnharmonic.harmonicFrequencies = newVectorDouble(vAnharmonic.nFrequencies);
	initVectorDouble(vAnharmonic.harmonicFrequencies, vAnharmonic.nFrequencies, 0.0);

	if(vData && vData->hessian)
	{
		int i=0;
		for(i=0;i<vAnharmonic.nFrequencies;i++) vAnharmonic.harmonicFrequencies[i] = vData->hessian[i][i];
	}
	vAnharmonic.Xi = newMatrixDouble(vAnharmonic.nFrequencies,vAnharmonic.nFrequencies);
	initMatrixDouble(vAnharmonic.Xi, vAnharmonic.nFrequencies, vAnharmonic.nFrequencies, 0.0);
	vAnharmonic.fundamentals = newVectorDouble(vAnharmonic.nFrequencies);
	initVectorDouble(vAnharmonic.fundamentals, vAnharmonic.nFrequencies, 0.0);

	vAnharmonic.overtones = newVectorDouble(vAnharmonic.nFrequencies);
	initVectorDouble(vAnharmonic.overtones, vAnharmonic.nFrequencies, 0.0);

	vAnharmonic.combinationBands = newMatrixDouble(vAnharmonic.nFrequencies,vAnharmonic.nFrequencies);
	initMatrixDouble(vAnharmonic.combinationBands, vAnharmonic.nFrequencies, vAnharmonic.nFrequencies, 0.0);

	return vAnharmonic;
}
/*****************************************************************************/
static void readData(VPT2Model* vpt2Model, char* inputFileName)
{
	FILE* inputFile;
        inputFile = fopen(inputFileName,"rb");
	if(!inputFile)
	{
		fprintf(stderr, "==========================================================\n");
		fprintf(stderr, "Sorry, I cannot opent the %s file\n", inputFileName);
		fprintf(stderr, "==========================================================\n");
		exit(1);
	}
	fclose(inputFile);
	vpt2Model->vData.klass->readData(&vpt2Model->vData, inputFileName);
	vpt2Model->pData.klass->readData(&vpt2Model->pData, inputFileName);
}
/**********************************************************************/
static void computeAnharmonicEnergies(VPT2Model* vpt2Model)
{
	computeXi(vpt2Model);
	printf("\nXi (cm^-1)\n");
	printMatrixDoubleCutOff(vpt2Model->vAnharmonic.Xi, vpt2Model->vAnharmonic.nFrequencies, vpt2Model->vAnharmonic.nFrequencies,1e-10);
	printf("END\n\n");

	computeFundamentalsEnergies(&vpt2Model->vAnharmonic);
	computeOvertonesEnergies(&vpt2Model->vAnharmonic);
	computeCombinationBandsEnergies(&vpt2Model->vAnharmonic);
	printIREnergies(&vpt2Model->vAnharmonic);
}
/**********************************************************************/
/*
static double applyModelProp(double k2, double epsilon, VPT2PropModel* pModel)
{
	if(pModel->type==MODEL_PROP_DCPT2) return DCPT2(k2,epsilon);
	if(pModel->type==MODEL_PROP_HDCPT2) return HDCPT2(k2,epsilon,pModel->alphaHDCPT2,pModel->betaHDCPT2);
	return VPT2(k2,epsilon);
}
*/
/**********************************************************************/
static double applyModelPropNew(VPT2Model* vpt2Model, double k2, int ii, int jj, int kk)
{
	VPT2PropModel* pModel = &vpt2Model->pData.model;
	VPT2PotentialData* vData = &vpt2Model->vData;
	double wi=(ii>0)?vData->hessian[ii-1][ii-1]:-vData->hessian[-ii-1][-ii-1];
	double wj=(jj>0)?vData->hessian[jj-1][jj-1]:-vData->hessian[-jj-1][-jj-1];
	double wk=(kk>0)?vData->hessian[kk-1][kk-1]:-vData->hessian[-kk-1][-kk-1];
	double epsilon = wi+wj+wk;
	if(pModel->type==MODEL_PROP_DCPT2) return DCPT2(k2,epsilon);
	if(pModel->type==MODEL_PROP_HDCPT2) return HDCPT2(k2,epsilon,pModel->alphaHDCPT2,pModel->betaHDCPT2);
	if(pModel->type==MODEL_PROP_VPT2) return VPT2(k2,epsilon);
	return VPT2Fermi(vpt2Model,k2,ii,jj,kk);
}
/**********************************************************************/
/* PCCP, 2014, 16, 1759-1787, page 1763-4 */
/* CPL, 496 (2010) 157–161 */
/**********************************************************************/
static void computeFundamentalsProp(VPT2Model* vpt2Model)
{
	int i,j,k,l;
	int a,c;
	double s = 0;
	double* Pav = NULL;
	double* V = NULL;
	double mu0 = 4*PI*1e-7;
	double eps0 = 1.0/(mu0*slight*slight);
	//double MWQ2q  = hPlank/4/PI/slight;
	double   kmmolm1 = 8*PI*PI*PI*NAvogadro/3/hPlank/slight/4/PI/eps0*1e-3*100.0*8.47835267e-30*8.47835267e-30;/* 1e-3 m to km, 100 : cm-1 to m-1 */
	double s0 = 1.0/sqrt(2.0);
	double s1 = s0/2;
	double s2 = s0/6;
	double S = 1.0;
	//VPT2PropModel* pModel = NULL;
	VPT2PotentialData vData;
	VPT2PropertiesData pData;
	VAnharmonic vAnharmonic;
	PropertiesAnharmonic pAnharmonic;

	if(!vpt2Model) return;
	if(printMax) printf("kmmolm1=%f\n",kmmolm1);
	pData = vpt2Model->pData;
	vData = vpt2Model->vData;
	//pModel = &pData.model;
	pAnharmonic = vpt2Model->pAnharmonic;
	vAnharmonic = vpt2Model->vAnharmonic;

	Pav = malloc(pData.nDim*sizeof(double));
	V = malloc(pData.nDim*sizeof(double));

	for(i=0;i<vData.nFrequencies;i++) 
	{
		double wi=vData.hessian[i][i];
		double si=1/(wi);
		double sqrti=(wi>0)?sqrt(1.0/wi):0;

		pAnharmonic.harmonic[i]     = 0.0;
		pAnharmonic.fundamentals[i] = 0.0;

		// HARMONIC
		for(a=0;a<pData.nDim;a++) Pav[a] = 0;
		for(a=0;a<pData.nDim;a++) V[a] = 0;
		for(a=0;a<pData.nDim;a++)
		{
			double Pi = pData.first[a][i]*sqrti;
			V[a] = s0*S*Pi;
		}
		if(printMax){
		printf("%s %d %s %10.4f %20s","Mode = ",i+1,"Harmonic Frequency = ",wi," Harmonic term ");
		for(a=0;a<pData.nDim;a++) printf("%f ",V[a]);
		printf("\n");
		}

		for(a=0;a<pData.nDim;a++) Pav[a] += V[a];
		for(a=0;a<pData.nDim;a++) V[a] = 0;

		for(a=0;a<pData.nDim;a++)
		{
			s = 0;
			for(j=0;j<vData.nFrequencies;j++) 
			{
				double wj = vData.hessian[j][j];
				double sqrtijj=(wi*wj*wj>0)?sqrt(1.0/(wi*wj*wj)):0;
				//double Pjji = getMaxCubicIJK(pData.cubic[a],j,j,i)*sqrtijj;
				double Pjji = (pData.cubic[a][j][j][i])*sqrtijj;
				double Pijj = (pData.cubic[a][i][j][j])*sqrtijj;
				double Pjij = (pData.cubic[a][j][i][j])*sqrtijj;
				//printf("a=%d P%d%d%d = %f %f %f\n",a,j+1,j+1,i+1,Pjji,Pijj,Pjij);
				s+= (Pjji+Pijj+S*Pjij);
			}
			V[a] += s*s2/2;
		}
		if(printMax){
		printf("%s %d %s %10.4f %20s","Mode = ",i+1,"Harmonic Frequency = ",wi," Electric term ");
		for(a=0;a<pData.nDim;a++) printf("%f ",V[a]);
		printf("\n");
		}

		for(a=0;a<pData.nDim;a++) Pav[a] += V[a];
		for(a=0;a<pData.nDim;a++) V[a] = 0;

		for(a=0;a<pData.nDim;a++)
		{
			s = 0;
			for(j=0;j<vData.nFrequencies;j++) 
			{
				double wj = vData.hessian[j][j];
				double sqrtj=(wj>0)?sqrt(1.0/wj):0;
				double Pj = pData.first[a][j]*sqrtj;
				double simj=(i==j || fabs(wi-wj)<1e-13)?0:1/(wi-wj);
				double A = (1/(wi+wj))*Pj;
				double B = -simj*S*Pj;
				for(k=0;k<vData.nFrequencies;k++) 
				{
					//double wk=vData.hessian[k][k];
					double kijkk = getMaxQuarticIJKL(vData.quartic, i, j, k, k);
					s += kijkk*A;
					if(!testDegen(vpt2Model, wi,wj) && !testResonanceQuartic(vpt2Model,wi,wj,kijkk)) 
						s+= kijkk*B;
				}
			}
			V[a] += -s*s0/8;
		}
		if(printMax){
		printf("%s %d %s %10.4f %20s","Mode = ",i+1,"Harmonic Frequency = ",wi," Quartic term ");
		for(a=0;a<pData.nDim;a++) printf("%f ",V[a]);
		printf("\n");
		}

		for(a=0;a<pData.nDim;a++) Pav[a] += V[a];
		for(a=0;a<pData.nDim;a++) V[a] = 0;

		for(a=0;a<pData.nDim;a++)
		{
			s = 0;
			for(j=0;j<vData.nFrequencies;j++) 
			{
				double wj=vData.hessian[j][j];
				double sqrtij=(wi*wj>0)?sqrt(1.0/(wi*wj)):0;
				//double Pji = getMaxMatrixIJ(pData.second[a],j,i)*sqrtij;
				double Pji = (pData.second[a][j][i])*sqrtij;
				double Pij = (pData.second[a][i][j])*sqrtij;

				for(k=0;k<vData.nFrequencies;k++) 
				{
					double wk=vData.hessian[k][k];
					double sqrtjk=(wj*wk>0)?sqrt(1.0/(wj*wk)):0;
					double kijk = getMaxCubicIJK(vData.cubic,i,j,k);
					double kjkk = getMaxCubicIJK(vData.cubic,j,k,k);
					//double Pjk = getMaxMatrixIJ(pData.second[a],j,k)*sqrtjk;
					//double Pkj = getMaxMatrixIJ(pData.second[a],k,j)*sqrtjk;
					double Pjk = (pData.second[a][j][k])*sqrtjk;
					double Pkj = (pData.second[a][k][j])*sqrtjk;


					s += kijk*(Pjk+Pkj)*(1/(wi+wj+wk));
					//s += -(Pjk+Pkj)*S*applyModelProp(kijk,wi-wj-wk,pModel);
					s += -(Pjk+Pkj)*S*applyModelPropNew(vpt2Model,kijk,i+1,-(j+1),-(k+1));
					s += kjkk/wj*(2*S*Pji+(1+S)*Pij);
				}
			}
			V[a] += -s*s1/8;
		}
		if(printMax){
		printf("%s %d %s %10.4f %20s","Mode = ",i+1,"Harmonic Frequency = ",wi," Mixed term ");
		for(a=0;a<pData.nDim;a++) printf("%f ",V[a]);
		printf("\n");
		}

		for(a=0;a<pData.nDim;a++) Pav[a] += V[a];
		for(a=0;a<pData.nDim;a++) V[a] = 0;

		for(a=0;a<pData.nDim;a++)
		{
			s = 0;
			for(j=0;j<vData.nFrequencies;j++) 
			{
				double wj=vData.hessian[j][j];
				double simj=(i==j)?0:1/(wi-wj);
				double sij=1/(wi+wj);
				double sqrtj=(wj>0)?sqrt(1.0/wj):0;
				double Pj = pData.first[a][j]*sqrtj;
				double A;
				double B;

				if(fabs(Pj)<1e-10) continue;

				for(k=0;k<vData.nFrequencies;k++) 
				{
					double wk=vData.hessian[k][k];
					double sum = 0;
					for(c=0;c<3;c++)
					{
						double cik = getMaxMatrixIJ(vData.coriolis[c],i,k);
						double cjk = getMaxMatrixIJ(vData.coriolis[c],j,k);
						double Be = vData.Be[c];
						sum += Be*cik*cjk;
					}
					A = sij*(sqrt(wi*wj)/wk-wk/sqrt(wi*wj))*Pj;
					s += A*sum;
					if(!testDegen(vpt2Model, wi,wj) && !testResonanceQuartic(vpt2Model,wi,wj,sum)) 
					{
						//printf("xwi = %f wj = %f sum = %f pas de res \n",wi,wj,sum);
						B = S*simj*(sqrt(wi*wj)/wk+wk/sqrt(wi*wj))*Pj;
						s += B*sum;
					}
				}
			}
			V[a] += s*s0/2;
		}
		if(printMax){
		printf("%s %d %s %10.4f %20s","Mode = ",i+1,"Harmonic Frequency = ",wi," Coriolis term ");
		for(a=0;a<pData.nDim;a++) printf("%f ",V[a]);
		printf("\n");
		}
		for(a=0;a<pData.nDim;a++) Pav[a] += V[a];
		for(a=0;a<pData.nDim;a++) V[a] = 0;


		for(a=0;a<pData.nDim;a++)
		{
			s = 0;
			for(j=0;j<vData.nFrequencies;j++) 
			{
				double wj=vData.hessian[j][j];
				double sj = 1/wj;
				double sij=1/(wi+wj);
				double siij=1/(wi+wi+wj);
				double sqrtj=(wj>0)?sqrt(1.0/wj):0;
				double Pj = pData.first[a][j]*sqrtj;
				double dij = (i==j)?1.0:0.0;
				double ddij = (i==j)?0.0:1.0;
				double simj=(i==j)?0:1.0/(wi-wj);

				if(fabs(Pj)<1e-10) continue;
				for(k=0;k<vData.nFrequencies;k++) 
				{
					double wk=vData.hessian[k][k];
					double sk = 1/wk;
					double sijk  = 1/(wi+wj+wk);
					double siik  = 1/(wi+wi+wk);
					double dik = (i==k)?1:0.0;
					double ddik = (i==k)?0:1.0;
					double kijk = getMaxCubicIJK(vData.cubic,i,j,k);

					for(l=0;l<vData.nFrequencies;l++) 
					{
						double wl=vData.hessian[l][l];
						double sjkl=1/(wj+wk+wl);
						double sikl=1/(wi+wk+wl);
						double dil = (i==l)?1:0.0;
						double ddil = (i==l)?0:1.0;
						double dijk = ddij*ddik*ddil;
						double kikl = getMaxCubicIJK(vData.cubic,i,k,l);
						double kjkl = getMaxCubicIJK(vData.cubic,j,k,l);
						double kllk = getMaxCubicIJK(vData.cubic,l,l,k);
						double kABC= kikl*kjkl;
						double kDEF=  kijk*kllk;
						double t;

						double A = 0.0;
						double B = 0.0;
						double C = 0.0;
						double D = 0.0;
						double E = 0.0;
						double F = 0.0;

						A = kABC*dijk*(sij*sjkl+S*sikl*sjkl+sij*sikl);
					        //A  += -applyModelProp(kABC*dijk*sjkl,wi-wk-wl,pModel);
					        A  += -applyModelPropNew(vpt2Model, kABC*dijk*sjkl,i+1,-(k+1),-(l+1));

						B = kABC*dij*(1+dik)*ddil*(0.5*si*sikl+S*0.5*sikl*sikl);
					        //t = applyModelProp(sqrt(kABC*dij*(1+dik)*ddil*0.5),wi-wk-wl,pModel);
					        t = applyModelPropNew(vpt2Model, sqrt(kABC*dij*(1+dik)*ddil*0.5),i+1,-(k+1),-(l+1));
					        B  += -S*t*t;
					        //B  += -applyModelProp(kABC*dij*(1+dik)*ddil*0.5*si,wi-wk-wl,pModel);
					        B  += -applyModelPropNew(vpt2Model, kABC*dij*(1+dik)*ddil*0.5*si,i+1,-(k+1),-(l+1));

						C = kABC*ddij*ddik*dil*(sk*sij+2.0*siik*sij+3.0*sij*sijk+2*S*siik*sijk+3*sk*sijk);
					        //C  += -S*applyModelProp(kABC*ddij*ddik*dil*sk,wi-wj-wk,pModel);
					        C  += -S*applyModelPropNew(vpt2Model, kABC*ddij*ddik*dil*sk,i+1,-(j+1),-(k+1));

						if(!testDegen(vpt2Model, wi,wj) && !testResonanceCubic(vpt2Model,wi,wj,kikl*kjkl))
						{ 
							A += kABC*ddik*ddil*S*simj*(-sjkl);
					        	//A  += S*applyModelProp(kABC*ddik*ddil*simj,wi-wk-wl,pModel);
					        	A  += S*applyModelPropNew(vpt2Model, kABC*ddik*ddil*simj,i+1,-(k+1),-(l+1));

							C += kABC*ddik*dil*S*simj*(-2*sijk-3*sk);
					        	//C  += S*applyModelProp( kABC*ddik*dil*simj,wi-wj-wk,pModel);
					        	C  += S*applyModelPropNew(vpt2Model,  kABC*ddik*dil*simj,i+1,-(j+1),-(k+1));
						}

						D  = kDEF*dij*si*sk*(1.0+(6.0-4.0*S)*dik*dil/9.0);
						E  = kDEF*dijk*(sij*sijk+sk*sij+sk*sijk);
					        //E  += -S*applyModelProp(kDEF*dijk*sk,wi-wj-wk,pModel);
					        E  += -S*applyModelPropNew(vpt2Model, kDEF*dijk*sk,i+1,-(j+1),-(k+1));
						F  = kDEF*dik*ddij*(1.0+dil)*(siij*sij+si*siij);
						F += kDEF*dik*ddij*(dil)*(1.0/3.0*si*sij+S/3.0*si*siij);
						F += kDEF*dik*ddij*(si*sij+S*si*sj);

						if(!testDegen(vpt2Model, wi,wj) && !testResonanceCubic(vpt2Model,wi,wj,kijk*kllk))
						{ 
							E += kDEF*ddik*ddil*S*simj*(-sk);
					        	//E += S*applyModelProp(kDEF*ddik*ddil*simj,wi-wj-wk,pModel);
					        	E += S*applyModelPropNew(vpt2Model, kDEF*ddik*ddil*simj,i+1,-(j+1),-(k+1));
							F += -kDEF*dik*(1.0+dil)*S*simj*si;
							F += -kDEF*dik*(dil)*S*simj*siij;
							F += -kDEF*dik*S*simj*sj;
						}
						s += (A+B+C)*Pj;
						s += (D+E+F)*Pj;
					}
				}
			}
			V[a] += s*s0/16.0;
		}
		if(printMax){
		printf("%s %d %s %10.4f %20s","Mode = ",i+1,"Harmonic Frequency = ",wi," Cubic term ");
		for(a=0;a<pData.nDim;a++) printf("%f ",V[a]);
		printf("\n");
		}

		for(a=0;a<pData.nDim;a++) Pav[a] += V[a];
		for(a=0;a<pData.nDim;a++) V[a] = 0;

		if(printMax){
		printf("%s %d %s %10.4f %20s","Mode = ",i+1,"Harmonic Frequency = ",wi," SUM ");
		for(a=0;a<pData.nDim;a++) printf("%f ",Pav[a]);
		printf("\n");
		}

		pAnharmonic.fundamentals[i] = 0.0;
		for(a=0;a<pData.nDim;a++) pAnharmonic.fundamentals[i] += Pav[a]*Pav[a];
		pAnharmonic.fundamentals[i] *= kmmolm1*vAnharmonic.fundamentals[i];
		pAnharmonic.harmonic[i] = 0.0;
		for(a=0;a<pData.nDim;a++) pAnharmonic.harmonic[i] += pData.first[a][i]*pData.first[a][i]/2;
		pAnharmonic.harmonic[i] *=  kmmolm1;
	}
	
	free(Pav);
	free(V);
}
/**********************************************************************/
static void printFundamentalsProp(VPT2Model* vpt2Model)
{
	int i;
	VAnharmonic vAnharmonic;
	PropertiesAnharmonic pAnharmonic;
	if(!vpt2Model) return;
	vAnharmonic = vpt2Model->vAnharmonic;
	pAnharmonic = vpt2Model->pAnharmonic;

	printf("Fundamentals\n");
	printf("------------\n");
	printf("%-10s %-10s %-20s %-20s %-20s %-20s\n","Mode","Quanta","E(harm)","E(anharm)", "I(harm)","I(anharm)");
	for(i=0;i<vAnharmonic.nFrequencies;i++) 
		printf("%-10d %-10d %-20.8f %-20.8f %-20.8f %-20.8f\n",i+1,1,
		vAnharmonic.harmonicFrequencies[i], vAnharmonic.fundamentals[i],
		pAnharmonic.harmonic[i], pAnharmonic.fundamentals[i]);
	printf("\n");
}
/**********************************************************************/
/* JCP 136, 124108, 2012 Bloino & Barone */
/* JPCA 2015, Bloino DOI: 10.1021/jp509985u */
static void computeOvertonesProp(VPT2Model* vpt2Model)
{
	int i,k;
	int a;
	double* Pav;
	double s;
	double mu0 = 4*PI*1e-7;
	double eps0 = 1.0/(mu0*slight*slight);
	double   kmmolm1 = 8*PI*PI*PI*NAvogadro/3/hPlank/slight/4/PI/eps0*1e-3*100.0*8.47835267e-30*8.47835267e-30;/* 1e-3 m to km, 100 : cm-1 to m-1 */
	double s0 = 1.0/sqrt(2.0);
	double s1 = s0/2;
	//double s2 = s0/6;
	double S = 1;
	//VPT2PropModel* pModel = NULL;
	VPT2PotentialData vData;
	VPT2PropertiesData pData;
	VAnharmonic vAnharmonic;
	PropertiesAnharmonic pAnharmonic;

	if(!vpt2Model) return;

	pData = vpt2Model->pData;
	vData = vpt2Model->vData;
	//pModel = &pData.model;
	pAnharmonic = vpt2Model->pAnharmonic;
	vAnharmonic = vpt2Model->vAnharmonic;

	Pav = malloc(pData.nDim*sizeof(double));
	for(i=0;i<pAnharmonic.nFrequencies;i++) 
	{
		double wi=vData.hessian[i][i];

		pAnharmonic.overtones[i] = 0.0;
		for(a=0;a<pData.nDim;a++) Pav[a] = 0;

		for(a=0;a<pData.nDim;a++)
		{
			double Pii=pData.second[a][i][i]/sqrt(wi*wi);
			Pav[a] += Pii*s1*S;
			s = 0;
			for(k=0;k<vData.nFrequencies;k++) 
			{
				double wk=vData.hessian[k][k];
				double Pk = pData.first[a][k]/sqrt(wk);
				double kiik = getMaxCubicIJK(vData.cubic, i, i, k);
				//s += S*applyModelProp(kiik,wi+wi-wk,pModel)*Pk;
				s += S*applyModelPropNew(vpt2Model, kiik,i+1,i+1,-(k+1))*Pk;
				s += -kiik*Pk/(wi+wi+wk);
			}
			Pav[a] += s*s0/4;
		}
		pAnharmonic.overtones[i] = 0.0;
		for(a=0;a<pData.nDim;a++) pAnharmonic.overtones[i] += Pav[a]*Pav[a];
		pAnharmonic.overtones[i] *= kmmolm1*vAnharmonic.overtones[i];
	}
	free(Pav);
}
/**********************************************************************/
static void printOvertonesProp(VPT2Model* vpt2Model)
{
	int i;
	VAnharmonic vAnharmonic;
	PropertiesAnharmonic pAnharmonic;
	if(!vpt2Model) return;
	vAnharmonic = vpt2Model->vAnharmonic;
	pAnharmonic = vpt2Model->pAnharmonic;
	printf("Overtones\n");
	printf("---------\n");
	printf("%-10s %-10s %-20s %-20s %-20s\n","Mode","Quanta","E(harm)","E(anharm)","I(anharm)");
	for(i=0;i<vAnharmonic.nFrequencies;i++) 
		printf("%-10d %-10d %-20.8f %-20.8f %-20.8f\n",i+1,2,
		2*vAnharmonic.harmonicFrequencies[i], 
		vAnharmonic.overtones[i],pAnharmonic.overtones[i]);
	printf("\n");
}
/**********************************************************************/
static void computeCombinationBandsProp(VPT2Model* vpt2Model)
{
	int i,j,k;
	int a;
	double* Pav;
	double s;
	double mu0 = 4*PI*1e-7;
	double eps0 = 1.0/(mu0*slight*slight);
	double   kmmolm1 = 8*PI*PI*PI*NAvogadro/3/hPlank/slight/4/PI/eps0*1e-3*100.0*8.47835267e-30*8.47835267e-30;/* 1e-3 m to km, 100 : cm-1 to m-1 */
	double s0 = 1.0/sqrt(2.0);
	double s1 = s0/2;
	//double s2 = s0/6;
	double S = 1;
	double sqrt2 = sqrt(2.0);
	//VPT2PropModel* pModel = NULL;
	VPT2PotentialData vData;
	VPT2PropertiesData pData;
	VAnharmonic vAnharmonic;
	PropertiesAnharmonic pAnharmonic;

	if(!vpt2Model) return;
	pAnharmonic = vpt2Model->pAnharmonic;
	vAnharmonic = vpt2Model->vAnharmonic;

	pData = vpt2Model->pData;
	vData = vpt2Model->vData;
	//pModel = &pData.model;

	Pav = malloc(pData.nDim*sizeof(double));
	for(i=0;i<pAnharmonic.nFrequencies;i++) 
	{
		double wi=vData.hessian[i][i];
		for(j=0;j<vData.nFrequencies;j++) 
		{
			double wj=vData.hessian[j][j];
			if(i==j) {vAnharmonic.combinationBands[i][j] = 0.0; continue;}
			for(a=0;a<pData.nDim;a++) Pav[a] = 0;
			for(a=0;a<pData.nDim;a++)
			{
				//double Pij=(pData.second[a][i][j]+pData.second[a][j][i])/sqrt(wi*wj);
				double Pij=getMaxMatrixIJ(pData.second[a],i,j)/sqrt(wi*wj);
				Pav[a] += s1*S*Pij;
				s = 0;
				for(k=0;k<vData.nFrequencies;k++) 
				{
					double wk=vData.hessian[k][k];
					double Pk = pData.first[a][k]/sqrt(wk);
					double kijk = getMaxCubicIJK(vData.cubic, i, j, k);
					//s += S*applyModelProp(kijk,wi+wj-wk,pModel)*Pk;
					s += S*applyModelPropNew(vpt2Model, kijk,i+1,j+1,-(k+1))*Pk;
					s += -kijk*Pk/(wi+wj+wk);
				}
				Pav[a] += s*s0/4;
			}
			for(a=0;a<pData.nDim;a++) Pav[a] *= sqrt2;

			pAnharmonic.combinationBands[i][j] = 0.0;
			for(a=0;a<pData.nDim;a++) pAnharmonic.combinationBands[i][j] += Pav[a]*Pav[a];
			pAnharmonic.combinationBands[i][j] *= kmmolm1*vAnharmonic.combinationBands[i][j];
		}
	}
	free(Pav);
}
/**********************************************************************/
static void printCombinationBandsProp(VPT2Model* vpt2Model)
{
	int i,j;
	VAnharmonic vAnharmonic;
	PropertiesAnharmonic pAnharmonic;

	if(!vpt2Model) return;
	pAnharmonic = vpt2Model->pAnharmonic;
	vAnharmonic = vpt2Model->vAnharmonic;
	printf("Combination Bands\n");
	printf("-----------------\n");
	printf("%-10s %-10s %-10s %-10s %-20s %-20s %-20s\n","Mode","Quanta","Mode","Quanta","E(harm)","E(anharm)","I(anharm)");
	for(i=0;i<vAnharmonic.nFrequencies;i++) 
	for(j=0;j<i;j++) 
		printf("%-10d %-10d %-10d %-10d %-20.8f %-20.8f %-20.8f\n",i+1,1,j+1,1,
		vAnharmonic.harmonicFrequencies[i]+vAnharmonic.harmonicFrequencies[j], 
		vAnharmonic.combinationBands[i][j],pAnharmonic.combinationBands[i][j]);
	printf("\n");
}
/**********************************************************************/
static void printIRProp(VPT2Model* vpt2Model)
{
	printFundamentalsProp(vpt2Model);
	printOvertonesProp(vpt2Model);
	printCombinationBandsProp(vpt2Model);
}
/**********************************************************************/
static boolean testDegen(VPT2Model* vpt2Model, double w1, double w2)
{
        double maxFrequencyDifference11Resonance = vpt2Model->pAnharmonic.maxFrequencyDifference11Resonance;
	if(fabs(w2-w1)<maxFrequencyDifference11Resonance) return TRUE;
	return FALSE;
}
/**********************************************************************/
static boolean testResonanceFreq(VPT2Model* vpt2Model, double w1, double w2)
{
	double* parameters11Resonance = vpt2Model->pAnharmonic.parameters11Resonance;
        double maxFrequencyDifference11Resonance = vpt2Model->pAnharmonic.maxFrequencyDifference11Resonance;
	double min = w1;
	double max = w1;
	if(min>w2) min = w2;
	if(max<w2) max = w2;
	//printf(" w1 = %f w2 = %f thr = %f\n",w1,w2,maxFrequencyDifference11Resonance);
	if(min>parameters11Resonance[0] && max<parameters11Resonance[1] && fabs(w2-w1)<parameters11Resonance[2]) return TRUE;
	if(fabs(w2-w1)<maxFrequencyDifference11Resonance) return TRUE;
	return FALSE;
}
/**********************************************************************/
static boolean testResonanceCubic(VPT2Model* vpt2Model, double w1, double w2, double v)
{
	if(testResonanceFreq(vpt2Model,w1,w2) && fabs(v)>= vpt2Model->pAnharmonic.thresholds11Numerators[0]) return TRUE;
	return FALSE;
}
/**********************************************************************/
static boolean testResonanceQuartic(VPT2Model* vpt2Model, double w1, double w2, double v)
{
	if(testResonanceFreq(vpt2Model,w1,w2) && fabs(v)>= vpt2Model->pAnharmonic.thresholds11Numerators[1]) return TRUE;
	return FALSE;
}
/**********************************************************************/
static PropertiesAnharmonic  newPAnharmonic(VPT2PotentialData*  vData)
{
	PropertiesAnharmonic pAnharmonic;
	int i;
	int n = 0;
	if(vData) n = vData->nFrequencies;
	pAnharmonic.nFrequencies = n;
	if(pAnharmonic.nFrequencies<=0) pAnharmonic.nFrequencies = 0;

	pAnharmonic.harmonic = newVectorDouble(pAnharmonic.nFrequencies);
	initVectorDouble(pAnharmonic.harmonic, pAnharmonic.nFrequencies, 0.0);

	pAnharmonic.fundamentals = newVectorDouble(pAnharmonic.nFrequencies);
	initVectorDouble(pAnharmonic.fundamentals, pAnharmonic.nFrequencies, 0.0);

	pAnharmonic.overtones = newVectorDouble(pAnharmonic.nFrequencies);
	initVectorDouble(pAnharmonic.overtones, pAnharmonic.nFrequencies, 0.0);

	pAnharmonic.combinationBands = newMatrixDouble(pAnharmonic.nFrequencies,pAnharmonic.nFrequencies);
	initMatrixDouble(pAnharmonic.combinationBands, pAnharmonic.nFrequencies, pAnharmonic.nFrequencies, 0.0);

        pAnharmonic.maxFrequencyDifference11Resonance = 2;
        for( i=0;i<2;i++) pAnharmonic.thresholds11Numerators[i] = 10.0;
        for( i=0;i<3;i++) pAnharmonic.parameters11Resonance[i] = 0.0;

	return pAnharmonic;
}
/**********************************************************************/
static void computeAnharmonicProperties(VPT2Model* vpt2Model)
{
	computeFundamentalsProp(vpt2Model);
	computeOvertonesProp(vpt2Model);
	computeCombinationBandsProp(vpt2Model);
	printIRProp(vpt2Model);
}
/**********************************************************************/
VPT2Model newVPT2Model()
{
	VPT2Model vpt2Model;
	vpt2Model.klass = malloc(sizeof(VPT2ModelClass));
	vpt2Model.klass->readData = readData;
	vpt2Model.klass->computeAnharmonic =computeAnharmonic;
	vpt2Model.vData = newVPT2PotentialData(0);
	vpt2Model.pData = newVPT2PropertiesData(0,0);
	vpt2Model.vAnharmonic = newVAnharmonic(&vpt2Model.vData);
	vpt2Model.pAnharmonic = newPAnharmonic(&vpt2Model.vData);
	return vpt2Model;
}
/**********************************************************************/
static void computeAnharmonic(VPT2Model* vpt2Model)
{
	//int nf;
	if(!vpt2Model) return;
	vpt2Model->vAnharmonic = newVAnharmonic(&vpt2Model->vData);
	computeAnharmonicEnergies(vpt2Model);
	//nf = vpt2Model->vAnharmonic.nFrequencies;
	//freeCubeDouble(&vpt2Model->vData.cubic, nf, nf);
	//freeQuarticDouble(&vpt2Model->vData.quartic, nf, nf, nf);
	vpt2Model->pAnharmonic = newPAnharmonic(&vpt2Model->vData);
	//computeAnharmonicEnergies(vpt2Model);
	computeAnharmonicProperties(vpt2Model);
}
/**********************************************************************/
