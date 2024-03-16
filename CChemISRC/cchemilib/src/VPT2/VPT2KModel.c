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

/* VPT2KModel.c */
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "../VPT2/VPT2KModel.h"

static boolean printMax = FALSE;
//static boolean printMax = TRUE;

static void computeAnharmonic(VPT2KModel* vpt2KModel);
static void readData(VPT2KModel* vpt2KModel, char* inputFileName);
static void computeHarmonicProp(VPT2KModel* vpt2KModel, int i, double* Pav);
static void computeFundamentalsProp(VPT2KModel* vpt2KModel, int i, double* Pav);
static void computeOvertonesProp(VPT2KModel* vpt2KModel, int i, double* Pav);
static void computeCombinationBandsProp(VPT2KModel* vpt2KModel, int i, int j, double* Pav);
static double K22(VPT2KModel* vpt2KModel, int k, int l, int m, int n);
static int getDiff(VPT2KModel* vpt2KModel, int* v, int* iDiff);
/**********************************************************************/
static int getGround(VPT2KModel* vpt2KModel)
{
        int i,j;
	int jmax = 0;
        VKAnharmonic* vAnharmonic  = &vpt2KModel->vAnharmonic;
	PropertiesKAnharmonic* pAnharmonic = &vpt2KModel->pAnharmonic;
        int nStates = pAnharmonic->nStates;
	i=0;
	{
		for(j=0;j<nStates;j++)
			if(fabs(vAnharmonic->eigenVectors[i][j])>fabs(vAnharmonic->eigenVectors[i][jmax])) jmax = j;
		if(jmax!=0) 
		{
			fprintf(stderr,"\n*******************************************\n");
			fprintf(stderr,"Warning the ground state is not the 0 state\n");
			fprintf(stderr,"E(ground state)=%f\n",vAnharmonic->eigenValues[0]);
			fprintf(stderr,"E(0)=%f\n",vAnharmonic->eigenValues[jmax]);
			fprintf(stderr,"*******************************************\n\n");
		}
	}
	return jmax;
}
/**********************************************************************/
static void printOneState(VPT2KModel* vpt2KModel, int j, int j0,double cutoff)
{
        VKAnharmonic* vAnharmonic  = &vpt2KModel->vAnharmonic;
	PropertiesKAnharmonic* pAnharmonic = &vpt2KModel->pAnharmonic;
        int nFrequencies = vAnharmonic->nFrequencies;
        int nStates = pAnharmonic->nStates;
        int i;
        int k;
	
	printf("%14.8f ",vAnharmonic->harmonicEnergies[j]); 
	printf("%14.8f ",vAnharmonic->eigenValues[j]-vAnharmonic->eigenValues[j0]);
	printf("%14.8f ",pAnharmonic->harmonic[j]);
	printf("%14.8f ",pAnharmonic->anHarmonic[j]);
	printf("    ");
	for(i=0;i<nStates;i++)
	if(fabs(vAnharmonic->eigenVectors[i][j])>cutoff ) 
	{
		printf("[% 5.3f",vAnharmonic->eigenVectors[i][j]);
		for(k=0;k<nFrequencies;k++) 
		{
			if(vAnharmonic->v[i][k]>0) 
			{
				printf(" %2d(%d)",k+1,vAnharmonic->v[i][k]);
			}
		}
		printf("] ");
	}
	printf("\n");
}
/**********************************************************************/
static void printTitles()
{
	printf("%14s ","E(harm)");
	printf("%14s ","E(anharm)");
	printf("%14s ","I(harm)");
	printf("%14s ","I(anharm)");
	printf("    ");
	printf("%s ","  Coef. Mode(Quanta)");
	printf("\n");
	printf("%14s ","------------");
	printf("%14s ","------------");
	printf("%14s ","------------");
	printf("%14s ","------------");
	printf("    ");
	printf("%s ","--------------------");
	printf("\n");
}
/**********************************************************************/
static void printIREnergiesPropCompact(VPT2KModel* vpt2KModel, double cutoff)
{
        int i,j;
        VKAnharmonic* vAnharmonic  = &vpt2KModel->vAnharmonic;
	PropertiesKAnharmonic* pAnharmonic = &vpt2KModel->pAnharmonic;
        int nFrequencies = vAnharmonic->nFrequencies;
        int nStates = pAnharmonic->nStates;
	int j0 = getGround(vpt2KModel);
	int* iDiff = newVectorInt(nFrequencies);
	int nd = 0;
	printf("\n\n");
	printf("Results With intensities\n");
	printf("========================\n");

	printf("\n");
	printf("%s\n","Fundamental Bands");
	printf("%s\n","-----------------");
	printTitles();
	// Fundamentals 
	for(j=0;j<nStates;j++)
	{
		int imax = 0;
		for(i=0;i<nStates;i++)
			if(fabs(vAnharmonic->eigenVectors[i][j])>fabs(vAnharmonic->eigenVectors[imax][j])) imax = i;
		int* v = vpt2KModel->vAnharmonic.v[imax];
		nd = getDiff(vpt2KModel, v, iDiff);
		if(nd!=1 || (nd==1 && v[iDiff[0]]!=1)) continue;
		printOneState(vpt2KModel, j, j0, cutoff);
	}
	printf("\n");
	printf("%s\n","Overtones");
	printf("%s\n","---------");
	printTitles();
	// Overtones 
	for(j=0;j<nStates;j++)
	{
		int imax = 0;
		for(i=0;i<nStates;i++)
			if(fabs(vAnharmonic->eigenVectors[i][j])>fabs(vAnharmonic->eigenVectors[imax][j])) imax = i;
		int* v = vpt2KModel->vAnharmonic.v[imax];
		nd = getDiff(vpt2KModel, v, iDiff);
		if(nd!=1 || (nd==1 &&  v[iDiff[0]]!=2)) continue;
		printOneState(vpt2KModel, j, j0, cutoff);
	}
	// Combinaition 
	printf("\n");
	printf("%s\n","Combination Bands");
	printf("%s\n","----------------");
	printTitles();
	for(j=0;j<nStates;j++)
	{
		int imax = 0;
		for(i=0;i<nStates;i++)
			if(fabs(vAnharmonic->eigenVectors[i][j])>fabs(vAnharmonic->eigenVectors[imax][j])) imax = i;
		int* v = vpt2KModel->vAnharmonic.v[imax];
		nd = getDiff(vpt2KModel, v, iDiff);
		if(nd!=2) continue;
		printOneState(vpt2KModel, j, j0, cutoff);
	}
	printf("\n");
	free(iDiff);
}
/**********************************************************************/
static void printIREnergiesPropCompact2(VPT2KModel* vpt2KModel, double cutoff)
{
        int i,j;
        int k;
        VKAnharmonic* vAnharmonic  = &vpt2KModel->vAnharmonic;
	PropertiesKAnharmonic* pAnharmonic = &vpt2KModel->pAnharmonic;
        int nFrequencies = vAnharmonic->nFrequencies;
        int nStates = pAnharmonic->nStates;
	int j0 = getGround(vpt2KModel);
	printf("\n\n");
	printf("Results With intensities\n");
	printf("========================\n");
	printf("%14s ","E(harm)");
	printf("%14s ","E(anharm)");
	printf("%14s ","I(harm)");
	printf("%14s ","I(anharm)");
	printf("    ");
	printf("%s ","Coef Mode(Quanta)");
	printf("\n");
	printf("%14s ","-------");
	printf("%14s ","---------");
	printf("%14s ","-------");
	printf("%14s ","---------");
	printf("    ");
	printf("%s ","----------------");
	printf("\n");
	for(j=0;j<nStates;j++)
	{
		printf("%14.8f ",vAnharmonic->harmonicEnergies[j]); 
		printf("%14.8f ",vAnharmonic->eigenValues[j]-vAnharmonic->eigenValues[j0]);
		printf("%14.8f ",pAnharmonic->harmonic[j]);
		printf("%14.8f ",pAnharmonic->anHarmonic[j]);
		printf("    ");
		for(i=0;i<nStates;i++)
		if(fabs(vAnharmonic->eigenVectors[i][j])>cutoff ) 
		{
			printf("[% 5.3f",vAnharmonic->eigenVectors[i][j]);
			for(k=0;k<nFrequencies;k++) 
			{
				if(vAnharmonic->v[i][k]>0) 
				{
					printf(" %2d(%d)",k+1,vAnharmonic->v[i][k]);
				}
			}
			printf("] ");
		}
		printf("\n");
	}
	printf("\n");
}
/**********************************************************************/
/*
static void printIREnergiesCompact(VPT2KModel* vpt2KModel)
{
	int i,j,k;
        VKAnharmonic* vAnharmonic  = &vpt2KModel->vAnharmonic;
	int nFrequencies = vAnharmonic->nFrequencies;
	int nStates = vAnharmonic->nStates;
	int j0 = getGround(vpt2KModel);
	printf("\n\n");
	printf("Results Compact\n");
	printf("===============\n");
	for(i=0;i<nStates;i++)
	{
		int jmax;
		for(k=0;k<nFrequencies;k++)
			printf("%2d ",vAnharmonic->v[i][k]);
		jmax = 0;
		for(j=0;j<nStates;j++)
			if(fabs(vAnharmonic->eigenVectors[i][j])>fabs(vAnharmonic->eigenVectors[i][jmax])) jmax = j;
		printf("%14.8f ",vAnharmonic->eigenValues[jmax]-vAnharmonic->eigenValues[j0]);
		printf("\n");
	}
	printf("\n");
}
*/
/**********************************************************************/
static void printIREnergiesCompact2(VPT2KModel* vpt2KModel, double cutoff)
{
        int i,j;
        int k;
        VKAnharmonic* vAnharmonic  = &vpt2KModel->vAnharmonic;
        int nFrequencies = vAnharmonic->nFrequencies;
	int nStates = vAnharmonic->nStates;
	int j0 = getGround(vpt2KModel);
	printf("\n\n");
	printf("Results\n");
	printf("========================\n");
	printf("%14s ","E(harm)");
	printf("%14s ","E(anharm)");
	printf("    ");
	printf("%s ","Mode(Quanta), Coef");
	printf("\n");
	printf("%14s ","-------");
	printf("%14s ","---------");
	printf("    ");
	printf("%s ","----------------");
	printf("\n");
	for(j=0;j<nStates;j++)
	{
		int imax = 0;
		for(i=0;i<nStates;i++)
		{
			if(fabs(vAnharmonic->eigenVectors[i][j])>fabs(vAnharmonic->eigenVectors[imax][j])) imax = i;
		}
		int* v = vpt2KModel->vAnharmonic.v[imax];
		double E = 0;
		for(k=0;k<nFrequencies;k++)
			E += v[k]*vpt2KModel->vData.hessian[k][k];

		printf("%14.8f ",E);
		printf("%14.8f ",vAnharmonic->eigenValues[j]-vAnharmonic->eigenValues[j0]);
		printf("    ");
		for(i=0;i<nStates;i++)
		if(fabs(vAnharmonic->eigenVectors[i][j])>cutoff ) 
		{
			printf("[% 5.3f",vAnharmonic->eigenVectors[i][j]);
			for(k=0;k<nFrequencies;k++) 
			{
				if(vAnharmonic->v[i][k]>0) 
				{
					printf(" %2d(%d)",k+1,vAnharmonic->v[i][k]);
				}
			}
			printf("] ");
		}
		printf("\n");
	}
	printf("\n");
}
/**********************************************************************/
/*
static void printIREnergies(VKAnharmonic* vAnharmonic)
{
	int i,j,k;
	int nFrequencies = vAnharmonic->nFrequencies;
	int nStates = vAnharmonic->nStates;
	printf("Results\n");
	for(k=0;k<nFrequencies;k++) printf("%2s "," ");
	for(i=0;i<nStates;i++)
		printf("%14.8f ",vAnharmonic->eigenValues[i]);
	printf("\n");
	for(k=0;k<nFrequencies;k++) printf("%2s "," ");
	for(i=0;i<nStates;i++)
		printf("%14.8f ",vAnharmonic->eigenValues[i]-vAnharmonic->eigenValues[0]);
	printf("\n");

	for(i=0;i<nStates;i++)
	{
		for(k=0;k<nFrequencies;k++)
			printf("%2d ",vAnharmonic->v[i][k]);
		for(j=0;j<nStates;j++)
			printf("%14.8f ",vAnharmonic->eigenVectors[i][j]);
		printf("\n");
	}
	printf("\n");
}
*/
/**********************************************************************/
/*
static double D0(VPT2KModel* vpt2KModel, int ii, int jj, int kk)
{
	VPT2PotentialData* data = &vpt2KModel->vData;
	int i=abs(ii)-1;
	int j=abs(jj)-1;
	int k=abs(kk)-1;
	double wi=(ii>0)?data->hessian[i][i]:-data->hessian[i][i];
	double wj=(jj>0)?data->hessian[j][j]:-data->hessian[j][j];
	double wk=(kk>0)?data->hessian[k][k]:-data->hessian[k][k];
	return 1/(wi+wj+wk);
}
*/
/**********************************************************************/
static double D(VPT2KModel* vpt2KModel, int i, int j, int k)
{
	int n = vpt2KModel->vAnharmonic.nFrequencies;
	return vpt2KModel->vAnharmonic.D[i+n][j+n][k+n];
}
/**********************************************************************/
static void setD(VPT2KModel* vpt2KModel, int i, int j, int k, double value, int*** tag, int t)
{
	int n = vpt2KModel->vAnharmonic.nFrequencies;
	vpt2KModel->vAnharmonic.D[i+n][j+n][k+n] = value;
	tag[i+n][j+n][k+n] = t;
}
/**********************************************************************/
static void setDNull(VPT2KModel* vpt2KModel, int ii, int jj, int kk, int***tag)
{
	VPT2PotentialData* data = &vpt2KModel->vData;
	int i=abs(ii);
	int j=abs(jj);
	int k=abs(kk);
	if(data->hessian[j-1][j-1]>data->hessian[i-1][i-1])
	{
		int l=i;
		i=j;
		j=l;
	}
	if(data->hessian[k-1][k-1]>data->hessian[i-1][i-1])
	{
		int l=i;
		i=k;
		k=l;
	}
	setD(vpt2KModel,i,-j,-k,0.0,tag,1);
	setD(vpt2KModel,i,-k,-j,0.0,tag,1);
	setD(vpt2KModel,-j,i,-k,0.0,tag,1);
	setD(vpt2KModel,-j,-k,i,0.0,tag,1);
	setD(vpt2KModel,-k,i,-j,0.0,tag,1);
	setD(vpt2KModel,-k,-j,i,0.0,tag,1);
        setD(vpt2KModel,-i,j,k,0.0,tag,1);
        setD(vpt2KModel,-i,k,j,0.0,tag,1);
        setD(vpt2KModel,j,-i,k,0.0,tag,1);
        setD(vpt2KModel,j,k,-i,0.0,tag,1);
        setD(vpt2KModel,k,-i,j,0.0,tag,1);
        setD(vpt2KModel,k,j,-i,0.0,tag,1);
}
/**********************************************************************/
static boolean testDD(VPT2KModel* vpt2KModel, int ii, int jj, int kk, int ll)
{
	VPT2PotentialData* data = &vpt2KModel->vData;
	int i=abs(ii)-1;
	int k=abs(kk)-1;
	double wi;
	double wk;
	boolean res = FALSE;
	double sik = 0;
	double tolDD = 10.0;

	if(ii!=jj || kk != ll) return FALSE;
	wi=(ii>0)?data->hessian[i][i]:-data->hessian[i][i];
	wk=(kk>0)?data->hessian[k][k]:-data->hessian[k][k];
	sik =wi+wk;
	if(fabs(sik)<tolDD) 
	{
		double Kiikk = K22(vpt2KModel, i, i, k, k);
		if(fabs(Kiikk)>tolDD) res = TRUE;
		else res=FALSE;
		if(res && printMax)
		{
			fprintf(stderr,"Warning Darling-Dennison resonnance : ");
			fprintf(stderr,"i,k = %d %d dw=%f Kiikk=%f tolDD=%f\n",i+1,k+1,sik,Kiikk,tolDD);
		}
	}
	return res;
}
/**********************************************************************/
static double getMaxDiffVPT2Var(VPT2KModel* vpt2KModel, int i, int j, int k, double sijk)
{
	double fijk=getMaxCubicIJK(vpt2KModel->vData.cubic,i,j,k);
	double dijk = fabs(sijk);
	double X = 1e10;
	if(dijk>1e-10)  X = (fijk*fijk*fijk*fijk)/(dijk*dijk*dijk)/64.0;
	if(i==j||i==k||j==k) X *= 0.25;
	return X;
}
/**********************************************************************/
static boolean testFermi(VPT2KModel* vpt2KModel, int ii, int jj, int kk)
{
	VPT2PotentialData* data = &vpt2KModel->vData;
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
		//printf("fabs(sijk)=%f\n",fabs(sijk));
		if(Z>0)
		{
			double fijk=getMaxCubicIJK(vpt2KModel->vData.cubic,i,j,k)*sqrt8;
			if(i==j||i==k||j==k) fijk *= 0.5;
			if(fijk/fabs(sijk)>Z) 
			{
				res = TRUE;
				if(printMax)
				{
					fprintf(stderr,"Warning Fermi resonnance : ");
                                	fprintf(stderr,"k,l,m = %d %d %d dklm=%f Z=%f Zcut=%f\n",ii,jj,kk,sijk,fijk/fabs(sijk),Z);
				}
			}
		}
		else if(XI>0 && XII>0)
		{
		 	double X =  getMaxDiffVPT2Var(vpt2KModel, i, j, k, sijk);
			//printf("X=%f\n",X);
			if(i==j||i==k||j==k)
			{
				if(X>XI) 
				{
					res = TRUE;
					if(printMax)
					{
						fprintf(stderr,"Warning Fermi resonnance : ");
						fprintf(stderr,"k,l,m = %d %d %d dklm=%f X=%f XI=%f\n",ii,jj,kk,sijk,X,XI);
					}
				}
			}
			else if(X>XII) { 
				res = TRUE;
				if(printMax)
				{
					fprintf(stderr,"Warning Fermi resonnance : ");
					fprintf(stderr,"k,l,m = %d %d %d dklm=%f X=%f XI=%f\n",ii,jj,kk,sijk,X,XI);
				}
			}
		}
		//printf("res = %d\n",res);
	}
	return res;
}
/**********************************************************************/
static void initD(VPT2KModel* vpt2KModel)
{
	int k,l,m;
	VKAnharmonic* vAnharmonic = &vpt2KModel->vAnharmonic;
	VPT2PotentialData* data = &vpt2KModel->vData;
	int n = data->nFrequencies;
	int ires = 0;
	double dMax = 0;
	int*** tag = newCubeInt(2*n+1, 2*n+1, 2*n+1);
	initCubeDouble(vAnharmonic->D, 2*n+1, 2*n+1, 2*n+1, -1.0);
	initCubeInt(tag, 2*n+1, 2*n+1, 2*n+1, 0);

	for(k=-n;k<=n;k++)
	{
		int kk = abs(k)-1;
		double s =(k<0)?-1:1;
		double wk=(kk>=0 && kk<n)?data->hessian[kk][kk]*s:0;
		if(k==0) continue;
		for(l=-n;l<=n;l++)
		{
			int ll = abs(l)-1;
			double s =(l<0)?-1:1;
			double wl=(ll>=0 && ll<n)?data->hessian[ll][ll]*s:0;

			if(l==0) continue;
			for(m=-n;m<=n;m++)
			{
				boolean res = FALSE;
				int mm = abs(m)-1;
				double s =(m<0)?-1:1;
				double wm=(mm>=0 && mm<n)?data->hessian[mm][mm]*s:0;
				if(m==0) continue;
				if(tag[k+n][l+n][m+n]>0) continue;
				res = testFermi(vpt2KModel, k,l,m);
				if(res)
				{
					ires++;
					vAnharmonic->D[k+n][l+n][m+n]=0.0;
					setDNull(vpt2KModel, kk+1, ll+1, mm+1,tag);
				}
				else 
				{
					if(fabs(wk+wl+wm)<1e-10) printf("Error !!!!\n");
					setD(vpt2KModel,k,l,m,1.0/(wk+wl+wm),tag,2);
				}
				if(fabs(vAnharmonic->D[k+n][l+n][m+n])>dMax) dMax = vAnharmonic->D[k+n][l+n][m+n];
			}
		}
	}
	if(printMax)
	{
		fprintf(stderr,"Dmax = %f\n",dMax);
		fprintf(stderr,"nresonance = %d\n",ires);
	}
	/*
	{
		FILE* file =fopen("res.in","w");
		for(k=1;k<n;k++)
		for(l=1;l<=k;l++)
		for(m=1;m<=l;m++)
			if(tag[k+n][-l+n][-m+n]==1) fprintf(file,"%d %d %d\n",k,l,m);
		fclose(file);
	}
	*/
	freeCubeInt(&tag,2*n+1,2*n+1);
	//if(ires>0) exit(1);
}
/**********************************************************************/
//1-2 Resonance (Fermi)
static double K12(VPT2KModel* vpt2KModel, int k, int l, int m)
{
	double fklm=getMaxCubicIJK(vpt2KModel->vData.cubic,k,l,m);
// k;l,l 
	if(l==m) return 0.5*fklm;
// General k;lm
	return fklm;
}
/**********************************************************************/
//2-2 Resonance 
static double K22(VPT2KModel* vpt2KModel, int k, int l, int m, int n)
{
	VPT2PotentialData* data = &vpt2KModel->vData;
	double fklmn = getMaxQuarticIJKL(data->quartic,k,l,m,n);
	double cor = 0;
	int i,c;
	double wk=data->hessian[k][k];
	double wl=data->hessian[l][l];
	double wm=data->hessian[m][m];
	double wn=data->hessian[n][n];
	double K = 0;
//	printf("K22 %d %d %d %d \n",k+1,l+1,m+1,n+1);
//      kk;mm
	if(k==l && m==n)
	{
		double cub1 = 0.0;
		double cub2 = 0.0;
		// double wkm = (wk+wm)/2;
		K = fklmn*0.25;
		cor = 0;
		for(c=0;c<3;c++)
		{
			double xeta=data->coriolis[c][k][m];
			if(fabs(data->coriolis[c][m][k])>fabs(xeta)) xeta=data->coriolis[c][m][k];
			cor += data->Be[c]*xeta*xeta;
		}
		K += -cor*(wk+wm)*(wk+wm)/wk/wm; // as VPT2+K
		//K += -cor*4;
		if(printMax) fprintf(stderr,"2-2 Kkkmm %d %d %d %d cor = %f\n",k+1,l+1,m+1,n+1,K-fklmn*0.25);
		
		for(i=0;i<data->nFrequencies;i++)
		{
			double fkki=getMaxCubicIJK(data->cubic,k,k,i);
			double fmmi=getMaxCubicIJK(data->cubic,m,m,i);
			double fkmi=getMaxCubicIJK(data->cubic,k,m,i);
			
			int kk=k+1;
			int mm=m+1;
			int ii=i+1;
                	cub1 += -1.0/16.0*fkki*fmmi*(D(vpt2KModel,kk,kk,ii)+D(vpt2KModel,-kk,-kk,ii)+D(vpt2KModel,mm,mm,ii)+D(vpt2KModel,-mm,-mm,ii));
                	cub2 += -0.25*fkmi*fkmi*(D(vpt2KModel,kk,-mm,ii)+D(vpt2KModel,-kk,mm,ii));
			
			/*
			double a1 = 4*wi/(wi*wi-4*wkm*wkm);
                	cub1 += -1.0/16.0*fkki*fmmi*a1;
                	cub2 += -0.25*fkmi*fkmi*(2.0/wi);
			*/
		}
		K += cub1;
		K += cub2;
		if(printMax) fprintf(stderr,"2-2 Kkkmm %d %d %d %d cub1 = %f\n",k+1,l+1,m+1,n+1,cub1);
		if(printMax) fprintf(stderr,"2-2 Kkkmm %d %d %d %d cub2 = %f\n",k+1,l+1,m+1,n+1,cub2);
		if(printMax) fprintf(stderr,"2-2 Kkkmm %d %d %d %d K = %f\n",k+1,l+1,m+1,n+1,K);
		return K;
	}
//      kl;mm or kk;lm
	if(k==l || m==n)
	{
		if(m==n) { i=m; m=k;n=l;k=l=i;}
                // m=>n l=>m in paper notation
                l=m;
                m=n;

		wk=data->hessian[k][k];
		wl=data->hessian[l][l];
		wm=data->hessian[m][m];

		K = fklmn*0.5;
		cor = 0;
		for(c=0;c<3;c++)
		{
			double xetakl=data->coriolis[c][k][l];
			double xetakm=data->coriolis[c][k][m];
			if(fabs(data->coriolis[c][l][k])>fabs(xetakl)) xetakl=data->coriolis[c][l][k];
			if(fabs(data->coriolis[c][m][k])>fabs(xetakm)) xetakm=data->coriolis[c][m][k];
			cor += data->Be[c]*xetakl*xetakm;
		}
		K += -2*cor*(wk+wl)*(wk+wm)/wk/sqrt(wl*wm);
		for(i=0;i<data->nFrequencies;i++)
		{
			double fkli=getMaxCubicIJK(data->cubic,k,l,i);
			double fkmi=getMaxCubicIJK(data->cubic,k,m,i);
			double fkki=getMaxCubicIJK(data->cubic,k,k,i);
			double flmi=getMaxCubicIJK(data->cubic,l,m,i);
			int kk=k+1;
			int mm=m+1;
			int ll=l+1;
			int ii=i+1;
                	K += -1.0/8.0*fkki*flmi*(D(vpt2KModel,kk,kk,ii)+D(vpt2KModel,-kk,-kk,ii)+D(vpt2KModel,ll,mm,ii)+D(vpt2KModel,-ll,mm,ii));
                	K += -0.25*fkli*fkmi*(D(vpt2KModel,kk,-ll,ii)+D(vpt2KModel,-kk,ll,ii)+D(vpt2KModel,kk,-mm,ii)+D(vpt2KModel,-kk,mm,ii));
		}
//		printf("2-2 Kkklm %d %d %d %d K = %f\n",k+1,k+1,l+1,m+1,K);
		return K;
	}

// General case kl;mn
		K = fklmn;
		cor = 0;
		for(c=0;c<3;c++)
		{
			double xetakl=getMaxMatrixIJ(data->coriolis[c],k,l);
			double xetakm=getMaxMatrixIJ(data->coriolis[c],k,m);
			double xetakn=getMaxMatrixIJ(data->coriolis[c],k,n);
			double xetalm=getMaxMatrixIJ(data->coriolis[c],l,m);
			double xetaln=getMaxMatrixIJ(data->coriolis[c],l,n);
			double xetamn=getMaxMatrixIJ(data->coriolis[c],m,n);

			cor +=  data->Be[c]*xetakl*xetamn*(wk-wl)*(wm-wn);
			cor += -data->Be[c]*xetakm*xetaln*(wk+wm)*(wl+wn);
			cor += -data->Be[c]*xetakn*xetalm*(wk+wn)*(wl+wm);
		}
		K += 2*cor/sqrt(wk*wl*wm*wn);

		for(i=0;i<data->nFrequencies;i++)
		{
			double fkli=getMaxCubicIJK(data->cubic,k,l,i);
			double fmni=getMaxCubicIJK(data->cubic,m,n,i);
			double fkmi=getMaxCubicIJK(data->cubic,k,m,i);
			double flni=getMaxCubicIJK(data->cubic,l,n,i);
			double fkni=getMaxCubicIJK(data->cubic,k,n,i);
			double flmi=getMaxCubicIJK(data->cubic,l,m,i);
			int kk=k+1;
			int mm=m+1;
			int ll=l+1;
			int nn=n+1;
			int ii=i+1;
                	K += -0.25*fkli*fmni*(D(vpt2KModel,kk,ll,ii)+D(vpt2KModel,mm,nn,ii)+D(vpt2KModel,-kk,-ll,ii)+D(vpt2KModel,-mm,-nn,ii));
                	K += -0.25*fkmi*flni*(D(vpt2KModel,kk,-mm,ii)+D(vpt2KModel,ll,-nn,ii)+D(vpt2KModel,-kk,mm,ii)+D(vpt2KModel,-ll,nn,ii));
                	K += -0.25*fkni*flmi*(D(vpt2KModel,kk,-nn,ii)+D(vpt2KModel,ll,-mm,ii)+D(vpt2KModel,-kk,nn,ii)+D(vpt2KModel,-ll,mm,ii));
		}
//		printf("2-2 Kklmn %d %d %d %d K = %f\n",k+1,l+1,m+1,n+1,K);
	return K;
}
/**********************************************************************/
//1-3 Resonance 
/*
static double K13(VPT2KModel* vpt2KModel, int k, int l, int m, int n)
{
	VPT2PotentialData* data = &vpt2KModel->vData;
	double fklmn = getMaxQuarticIJKL(data->quartic,k,l,m,n);
	double cor = 0;
	int i,c;
	double wk=data->hessian[k][k];
	double wl=data->hessian[l][l];
	double wm=data->hessian[m][m];
	double wn=data->hessian[n][n];
	double K = 0;

//	printf("K13 %d %d %d %d \n",k+1,l+1,m+1,n+1);
//      k;lll
	if(l==m && m==n)
	{
		K = fklmn/6.0;
		cor = 0;
		for(i=0;i<data->nFrequencies;i++)
		{
			double fkli=getMaxCubicIJK(data->cubic,k,l,i);
			double flli=getMaxCubicIJK(data->cubic,l,l,i);
			int kk=k+1;
			int ll=l+1;
			int ii=i+1;
                	K += -1.0/8*fkli*flli*(D(vpt2KModel,kk,-ll,ii)+D(vpt2KModel,-kk,ll,ii)+D(vpt2KModel,ll,ll,ii)+D(vpt2KModel,-ll,-ll,ii));
		}
		return K;
	}
//      k;llm
	if(l==m || m==n || l==n)
	{
		if(l==n) { i=m; m=n;n=i;}// permut m and n
		if(m==n) { i=n; n=l;l=i; }//permut n and l
		K = fklmn*0.5;
		cor = 0;
		for(c=0;c<3;c++)
		{
			double xetakl=getMaxMatrixIJ(data->coriolis[c],k,l);
			double xetalm=getMaxMatrixIJ(data->coriolis[c],l,m);
			cor += data->Be[c]*xetakl*xetalm;
		}
		K += 2*cor*(wk+wl)*(wl-wm)/wl/sqrt(wk*wm);
		for(i=0;i<data->nFrequencies;i++)
		{
			double fkli=getMaxCubicIJK(data->cubic,k,l,i);
			double fkmi=getMaxCubicIJK(data->cubic,k,m,i);
			double flli=getMaxCubicIJK(data->cubic,l,l,i);
			double flmi=getMaxCubicIJK(data->cubic,l,m,i);
			int kk=k+1;
			int mm=m+1;
			int ll=l+1;
			int ii=i+1;
                	K += -1.0/8.0*fkmi*flli*(D(vpt2KModel,ll,ll,ii)+D(vpt2KModel,-ll,-ll,ii)+D(vpt2KModel,kk,-mm,ii)+D(vpt2KModel,-kk,mm,ii));
                	K += -0.25*fkli*flmi*(D(vpt2KModel,kk,-ll,ii)+D(vpt2KModel,-kk,ll,ii)+D(vpt2KModel,ll,mm,ii)+D(vpt2KModel,-ll,-mm,ii));
		}
		return K;
	}

// General case k;lmn
		K = fklmn;
		cor = 0;
		for(c=0;c<3;c++)
		{
			double xetakl=getMaxMatrixIJ(data->coriolis[c],k,l);
			double xetakm=getMaxMatrixIJ(data->coriolis[c],k,m);
			double xetakn=getMaxMatrixIJ(data->coriolis[c],k,n);
			double xetalm=getMaxMatrixIJ(data->coriolis[c],l,m);
			double xetaln=getMaxMatrixIJ(data->coriolis[c],l,n);
			double xetamn=getMaxMatrixIJ(data->coriolis[c],m,n);

			cor +=  data->Be[c]*xetakl*xetamn*(wk+wl)*(wm-wn);
			cor +=  data->Be[c]*xetakm*xetaln*(wk+wm)*(wl-wn);
			cor +=  data->Be[c]*xetakn*xetalm*(wk+wn)*(wl-wn);
		}
		K += cor/sqrt(wk*wl*wm*wn);

		for(i=0;i<data->nFrequencies;i++)
		{
			double fkli=getMaxCubicIJK(data->cubic,k,l,i);
			double fmni=getMaxCubicIJK(data->cubic,m,n,i);
			double fkmi=getMaxCubicIJK(data->cubic,k,m,i);
			double flni=getMaxCubicIJK(data->cubic,l,n,i);
			double fkni=getMaxCubicIJK(data->cubic,k,n,i);
			double flmi=getMaxCubicIJK(data->cubic,l,m,i);
			int kk=k+1;
			int mm=m+1;
			int ll=l+1;
			int nn=n+1;
			int ii=i+1;
                	K += -0.25*fkli*fmni*(D(vpt2KModel,-kk,ll,ii)+D(vpt2KModel,-mm,-nn,ii)+D(vpt2KModel,kk,-ll,ii)+D(vpt2KModel,mm,nn,ii));
                	K += -0.25*fkmi*flni*(D(vpt2KModel,kk,-mm,ii)+D(vpt2KModel,ll,nn,ii)+D(vpt2KModel,-kk,mm,ii)+D(vpt2KModel,-ll,-nn,ii));
                	K += -0.25*fkni*flmi*(D(vpt2KModel,kk,-nn,ii)+D(vpt2KModel,ll,mm,ii)+D(vpt2KModel,-kk,nn,ii)+D(vpt2KModel,-ll,-mm,ii));
		}
	return K;
}
*/
/**********************************************************************/
//1-1 Resonance 
/*
static double K11(VPT2KModel* vpt2KModel, int k, int l, int m, int n)
{
	VPT2PotentialData* data = &vpt2KModel->vData;
	double fklmn = getMaxQuarticIJKL(data->quartic,k,l,m,n);
	double cor = 0;
	int i,c;
	double wk=data->hessian[k][k];
	double wl=data->hessian[l][l];
	double wm=data->hessian[m][m];
	double K = 0;
	int ne = 0;
	double s2 = 0.5;
	double s6 = 1.0/6.0;
	double s8 = 1.0/8.0;
	double s24 = 1.0/24.0;
	if(k==l) ne++;
	if(k==m) ne++;
	if(k==n) ne++;
	if(l==m) ne++;
	if(l==n) ne++;
	if(m==n) ne++;


//	printf("K11 %d %d %d %d ne = %d\n",k+1,l+1,m+1,n+1,ne);
//      kl;ml
	if(ne==1)
	{
		if(k==m || k==n) {i=l;l=k;k=i;}//permut kl
		if(l==m ) {i=m;m=n;n=i;}//permut mn

		wk=data->hessian[k][k];
		wl=data->hessian[l][l];
		wm=data->hessian[m][m];

		K = fklmn*s2;
		cor = 0;
		//printf("KQuart %d %d %d %d = %f\n",k+1,l+1,m+1,n+1,K);
		for(c=0;c<3;c++)
		{
			double xetakl=getMaxMatrixIJ(data->coriolis[c],k,l);
			double xetaml=getMaxMatrixIJ(data->coriolis[c],m,l);
			cor += data->Be[c]*xetakl*xetaml;
		}
		//printf("Cori %d %d %d %d = %f\n",k+1,l+1,m+1,n+1,2*cor*(wk*wm+wl*wl)/wl/sqrt(wk*wm));
		K += 2*cor*(wk*wm+wl*wl)/wl/sqrt(wk*wm);
		for(i=0;i<data->nFrequencies;i++)
		{
			double fkli=getMaxCubicIJK(data->cubic,k,l,i);
			double fmli=getMaxCubicIJK(data->cubic,m,l,i);
			double fkmi=getMaxCubicIJK(data->cubic,k,m,i);
			double flli=getMaxCubicIJK(data->cubic,l,l,i);
			double wi=data->hessian[i][i];
			int kk=k+1;
			int ll=l+1;
			int mm=m+1;
			int ii=i+1;
                	K += -s8*fkli*fmli*
				(
					D(vpt2KModel,kk,ll,ii)+D(vpt2KModel,-kk,-ll,ii)+D(vpt2KModel,mm,ll,ii)+D(vpt2KModel,-mm,-ll,ii)+
					D(vpt2KModel,kk,-ll,ii)+D(vpt2KModel,-kk,ll,ii)+D(vpt2KModel,ll,-mm,ii)+D(vpt2KModel,-ll,mm,ii)
				);
                	K += -s8*fkmi*flli* (D(vpt2KModel,kk,-mm,ii)+D(vpt2KModel,-kk,mm,ii)+2/wi);
		}
		//fprintf(stderr,"Kklml %d %d %d %d = %f\n",k+1,l+1,m+1,n+1,K);
		return K;
	}
//      kl;ll
	if(ne==3)
	{
		if(k==l )
		{
			i=k;k=m;m=i;
			i=l;l=n;n=i;
		}//permut kl & mn
		if(k==m || k==n) {i=l;l=k;k=i;}//permut kl

		K = fklmn*s6;
		cor = 0;
		//printf("KQuart %d %d %d %d = %f\n",k+1,l+1,m+1,n+1,K);
		for(i=0;i<data->nFrequencies;i++)
		{
			double fkli=getMaxCubicIJK(data->cubic,k,l,i);
			double flli=getMaxCubicIJK(data->cubic,l,l,i);
			double wi=data->hessian[i][i];
			int kk=k+1;
			int ll=l+1;
			int ii=i+1;
                	K += -s24*fkli*flli*(D(vpt2KModel,ll,ll,ii)+D(vpt2KModel,-ll,-ll,ii)+4/wi);
                	K += -s24*fkli*flli*(2*D(vpt2KModel,kk,-ll,ii)+2*D(vpt2KModel,-kk,ll,ii)+D(vpt2KModel,kk,ll,ii)+D(vpt2KModel,-kk,-ll,ii));
		}
		//fprintf(stderr,"Kklll %d %d %d %d = %f\n",k+1,l+1,m+1,n+1,K);
		return K;
	}

	//printf("K %d %d %d %d = %f\n",k+1,l+1,m+1,n+1,K);
	return K;
}
*/
/**********************************************************************/
/* VPT2-K paper */
static void computeXi(VPT2KModel* vpt2KModel)
{
	int k,l,m;
	int c;
	double s = 0;
	VKAnharmonic* vAnharmonic = &vpt2KModel->vAnharmonic;
	VPT2PotentialData* data = &vpt2KModel->vData;
	//VPT2VModel* model = NULL;
	if(!vAnharmonic) return;
	if(!data) return;
	if(!data->hessian) return;
	if(!vAnharmonic->Xi) return;
	//model = &data->model;
	/* compute Xkk */
	for(k=0;k<vAnharmonic->nFrequencies;k++)
	{
		double Skk = 0;
		vAnharmonic->Xi[k][k] = data->quartic[k][k][k][k]/16.0;
		Skk = 0;
		for(l=0;l<vAnharmonic->nFrequencies;l++)
		{
			double fkkl=getMaxCubicIJK(data->cubic,k,k,l);
			double wl=data->hessian[l][l];
			Skk += fkkl*fkkl*(2/wl+0.5*D(vpt2KModel,l+1,k+1,k+1)-0.5*D(vpt2KModel,-l-1,k+1,k+1));
		}
		vAnharmonic->Xi[k][k] += -Skk/16.0;
	}
	/* compute Xkl */
	for(k=0;k<vAnharmonic->nFrequencies;k++)
	{
		double wk=data->hessian[k][k];
		int kk=k+1;
		for(l=0;l<vAnharmonic->nFrequencies;l++)
		{
			double wl=data->hessian[l][l];
			double skl = wk/wl+wl/wk;
			int ll=l+1;
			if(k==l) continue;
			vAnharmonic->Xi[k][l] = getMaxQuarticIJKL(data->quartic,k,k,l,l)*0.25;
			for(m=0;m<vAnharmonic->nFrequencies;m++)
			{
				double wm=data->hessian[m][m];
				double fklm=getMaxCubicIJK(data->cubic,k,l,m);
				double fkkm=getMaxCubicIJK(data->cubic,k,k,m);
				double fmll=getMaxCubicIJK(data->cubic,m,l,l);
				int mm=m+1;
				vAnharmonic->Xi[k][l] += -0.25*fkkm*fmll/wm;
				vAnharmonic->Xi[k][l] += -0.125*fklm*fklm*
				( D(vpt2KModel,mm,kk, ll) + D(vpt2KModel,mm,-kk,-ll)
				+ D(vpt2KModel,mm,kk,-ll) + D(vpt2KModel,mm,-kk, ll));
			}
			s = 0;
			for(c=0;c<3;c++)
			{
				double xeta=data->coriolis[c][k][l];
				if(fabs(data->coriolis[c][l][k])>fabs(xeta)) xeta=data->coriolis[c][l][k];
				s += data->Be[c]*xeta*xeta;
			}
			vAnharmonic->Xi[k][l] += skl*s;
		}
	}
}
/**********************************************************************/
static void computeC(VPT2KModel* vpt2KModel)
{
	int k,l,m;
	int c;
	double s = 0;
	double C = 0;
	boolean resonance = FALSE;
	VKAnharmonic* vAnharmonic = &vpt2KModel->vAnharmonic;
	VPT2PotentialData* data = &vpt2KModel->vData;
	if(!vAnharmonic) return;
	if(!data) return;
	if(!data->hessian) return;
	if(!vAnharmonic->Xi) return;
	

	s = 0;
	for(k=0;k<vAnharmonic->nFrequencies;k++)
		s += data->quartic[k][k][k][k];
	s  = s/64.0;
	C += s;

	s = 0;
	for(k=0;k<vAnharmonic->nFrequencies;k++)
	{
		int kk=k+1;
		for(l=0;l<vAnharmonic->nFrequencies;l++)
		{
			int ll=l+1;
			double fkkl=getMaxCubicIJK(data->cubic,k,k,l);
			if(k!=l) s += fkkl*fkkl*0.5*(D(vpt2KModel,ll,kk,kk)+D(vpt2KModel,-ll,kk,kk));
		}
	}
	s = 3*s/64.0;
	C += s;

	s = 0;
	for(k=0;k<vAnharmonic->nFrequencies;k++)
	{
		double fkkk=getMaxCubicIJK(data->cubic,k,k,k);
		double wk=data->hessian[k][k];
		s += fkkk*fkkk/wk;
	}
	s  = -7*s/576;
	C += s;

	s = 0;
	for(k=0;k<vAnharmonic->nFrequencies;k++)
	{
		int kk=k+1;
		for(l=0;l<vAnharmonic->nFrequencies;l++)
		{
			int ll=l+1;
			for(m=0;m<vAnharmonic->nFrequencies;m++)
			{
				int mm=m+1;
             			if (
					(!resonance && (k!=l && l!=m && k!=m))
					||
            				( resonance && ((k>l) && l>m && k>m))
				)
				{
					double fklm=getMaxCubicIJK(data->cubic,k,l,m);
  					s += fklm*fklm*(D(vpt2KModel,mm,kk,ll)+D(vpt2KModel,mm,-kk,-ll)-D(vpt2KModel,mm,kk,-ll)-D(vpt2KModel,mm,-kk,ll))/8;
				}
			}
		}
	}
	s = -s /4.0;
	if(!resonance) s = s/6;
	C += s;

	s = 0;
	for(k=0;k<vAnharmonic->nFrequencies;k++)
	{
		for(l=0;l<vAnharmonic->nFrequencies;l++)
		{
                    	if ((!resonance && k!=l) || (resonance && k>l))
			{
				for(c=0;c<3;c++)
				{
					double xetakl = getMaxMatrixIJ(data->coriolis[c],k,l);
					s += -0.5*data->Be[c]*xetakl*xetakl;
				}
			}
		}
	}
	if(!resonance) s = s/2;
	C += s;

	vpt2KModel->vAnharmonic.C = C;
}
/**********************************************************************/
static double Ekl(VPT2KModel* vpt2KModel, int* v)
{
	int k,l;
	int nFrequencies = vpt2KModel->vAnharmonic.nFrequencies;
	double E = 0;
	double shift = 0.5;
	for(k=0;k<nFrequencies;k++)
	{
		double wk=vpt2KModel->vData.hessian[k][k];
		E += wk*(v[k]+shift);
		for(l=0;l<=k;l++)
			E += vpt2KModel->vAnharmonic.Xi[k][l]*(v[k]+shift)*(v[l]+shift);
	}
	return E;
}
/**********************************************************************/
static double Eoff12(VPT2KModel* vpt2KModel, int* vDiff, int* vMax, int noneZero, int negative)
{
	int nFrequencies = vpt2KModel->vAnharmonic.nFrequencies;
	//int i,k,l,m,n;
	int i,k,l,m;
	double E = 0;
//  1-2 Resonances
	k = -1;
	l = -1;
	m = -1;
	//n = -1;
//	printf("Eoff12 noneZero=%d negative=%d \n",noneZero,negative);
	if(noneZero==3) // Resonance type: k;lm
	{
		for(i=0;i<nFrequencies;i++) 
		if(negative==1)
		{
			if(vDiff[i]<0 && k==-1) k=i;
			else if (vDiff[i]>0 && l==-1) l=i;
			else if (vDiff[i]>0 && l!=-1 && m==-1) m=i;
		}
		else if(negative==2)
		{
			if(vDiff[i]>0 && k==-1) k=i;
			else if (vDiff[i]<0 && l==-1) l=i;
			else if (vDiff[i]<0 && l!=-1 && m==-1) m=i;
		}
		if(k==-1 || l==-1 || m==-1) {printf("No 12 reseonance \n");return E;};
		if(testFermi(vpt2KModel,-(k+1),l+1,m+1))
			E=K12(vpt2KModel,k,l,m)*sqrt(vMax[k]*vMax[l]*vMax[m])/2.0/sqrt(2.0);
		return E;
	}
	else if(noneZero==2) // Resonance type: k;ll
	{
		for(i=0;i<nFrequencies;i++) 
		if(negative==1)
		{
			if(abs(vDiff[i])==1 && k==-1) k=i;
			else if (abs(vDiff[i])==2 && l==-1) l=i;
		}
		if(k==-1 || l==-1) {printf("No 12 reseonance \n");return E;};
		if(testFermi(vpt2KModel,-(k+1),l+1,l+1))
			E=K12(vpt2KModel,k,l,l)*sqrt(vMax[k]*(vMax[l]-1)*vMax[l])/2.0/sqrt(2.0);
		return E;
	}
	return E;
}
/**********************************************************************/
static double Eoff22(VPT2KModel* vpt2KModel, int* vDiff, int* vMax, int noneZero, int negative)
{
	int nFrequencies = vpt2KModel->vAnharmonic.nFrequencies;
	int i,k,l,m,n;
	double E = 0;

	k = -1;
	l = -1;
	m = -1;
	n = -1;

	//kl;mn : nonZero = 4;
//	printf("Eoff22 noneZero=%d negative=%d \n",noneZero,negative);
	if(noneZero==4) // Resonance type: kl;mn
	{
		for(i=0;i<nFrequencies;i++) 
		{
			if(vDiff[i]==-1 && k==-1) k=i;
			else if (vDiff[i]==-1 && k!=-1) l=i;
			else if(vDiff[i]==1 && m==-1) m=i;
			else if (vDiff[i]==1 && m!=-1) n=i;
		}
		if(k==-1 || l==-1 || m==-1 || n==-1) {printf("No 22 reseonance noneZero=%d negative=%d \n",noneZero,negative);return E;};
		if(testDD(vpt2KModel,k+1,l+1,-(m+1),-(n+1)))
			E=0.25*K22(vpt2KModel,k,l,m,n)*sqrt(vMax[k]*vMax[l]*vMax[m]*vMax[n]);
		return E;
	}
	else if(noneZero==3) // Resonance type: kl;mm or kk;lm
	{
		for(i=0;i<nFrequencies;i++) 
		if(negative==1)
		{
			if(vDiff[i]<0) k=i;
			else if (vDiff[i]>0 && l==-1) l=i;
			else if (vDiff[i]>0 && l!=-1 && m==-1) m=i;
		}
		else if(negative==2)
		{
			if(vDiff[i]>0) k=i;
			//else if (vDiff[i]<0 && k==-1) l=i;
			else if (vDiff[i]<0 && l==-1) l=i;
			else if (vDiff[i]<0 && l!=-1 && m==-1) m=i;
		}
		if(k==-1 || l==-1 || m==-1) {printf("No 22 reseonance noneZero=%d negative=%d klm = %d %d %d\n",noneZero,negative,k,l,m);return E;};
//		printf(" Resonance type: kl;mm or kk;lm klmn %d %d %d %d\n",k,l,m,n);
		if(testDD(vpt2KModel,k+1,k+1,-(l+1),-(m+1)))
			E=0.25*K22(vpt2KModel,k,k,l,m)*sqrt((vMax[k]-1)*vMax[k]*vMax[l]*vMax[m]);
		return E;
	}
	else if(noneZero==2) // Resonance type: kk;ll
	{
		for(i=0;i<nFrequencies;i++) 
		{
			if(vDiff[i]<0) k=i;
			else if (vDiff[i]>0) l=i;
		}
		if(k==-1 || l==-1) {printf("No 22 reseonance noneZero=%d negative=%d \n",noneZero,negative);return E;};
		if(testDD(vpt2KModel,k+1,k+1,-(l+1),-(l+1)))
			E=0.25*K22(vpt2KModel,k,k,l,l)*sqrt((vMax[k]-1)*vMax[k]*(vMax[l]-1)*vMax[l]);
		// TO DO return E;
		return E/2;
	}
	return E;
}
/**********************************************************************/
/*
static double Eoff13(VPT2KModel* vpt2KModel, int* vDiff, int* vMax, int noneZero, int negative)
{
	int nFrequencies = vpt2KModel->vAnharmonic.nFrequencies;
	int i,k,l,m,n;
	double E = 0;

	k = -1;
	l = -1;
	m = -1;
	n = -1;
//	printf("Eoff13 noneZero=%d negative=%d \n",noneZero,negative);
	//k;lmn : nonZero = 4;
	if(noneZero==4) // Resonance type: k;lmn
	{
		for(i=0;i<nFrequencies;i++) 
		if(negative==1)
		{
			if(vDiff[i]<0) k=i;
			else if (vDiff[i]>0 && l!=-1) l=i;
			else if(vDiff[i]>0 && l!=-1 && m==-1) m=i;
			else if (vDiff[i]>0 && l!=-1 && m!=-1 && n==-1) n=i;
		}
		else if(negative==3)
		{
			if(vDiff[i]>0) k=i;
			else if (vDiff[i]<0 && l!=-1) l=i;
			else if(vDiff[i]<0 && l!=-1 && m==-1) m=i;
			else if (vDiff[i]<0 && l!=-1 && m!=-1 && n==-1) n=i;
		}
		if(k==-1 || l==-1 || m==-1 || n==-1) {printf("No 13 reseonance \n");return E;};
		E=0.25*K13(vpt2KModel,k,l,m,n)*sqrt(vMax[k]*vMax[l]*vMax[m]*vMax[n]);
		return E;
	}
	if(noneZero==3) // Resonance type: k;llm
	{
		for(i=0;i<nFrequencies;i++) 
		if(negative==1)
		{
			if(vDiff[i]<0) k=i;
			else if (abs(vDiff[i])>1 && l!=-1) l=i;
			else if(abs(vDiff[i])==1 && m==-1) m=i;
		}
		else if(negative==2)
		{
			if(vDiff[i]>0) k=i;
			else if (abs(vDiff[i])>1 && l!=-1) l=i;
			else if(abs(vDiff[i])==1 && m==-1) m=i;
		}
		if(k==-1 || l==-1 || m==-1) {printf("No 13 reseonance \n");return E;};
		E=0.25*K13(vpt2KModel,k,l,l,m)*sqrt(vMax[k]*(vMax[l]-1)*vMax[l]*vMax[m]);
		return E;
	}
	if(noneZero==2) // Resonance type: k;lll
	{
		for(i=0;i<nFrequencies;i++) 
		{
			if(abs(vDiff[i])==1) k=i;
			else if (abs(vDiff[i])==3) l=i;
		}
		if(k==-1 || l==-1) {printf("No 13 reseonance \n"); return 0;};
		E=0.25*K13(vpt2KModel,k,l,l,l)*sqrt(vMax[k]*(vMax[l]-1)*(vMax[l]-2)*vMax[l]);
		return E;
	}
	return E;
}
*/
/**********************************************************************/
/*
static double Eoff11(VPT2KModel* vpt2KModel, int* vDiff, int* vMax, int noneZero, int absSum)
{
	int nFrequencies = vpt2KModel->vAnharmonic.nFrequencies;
	int i,k,l,m,n;
	double E = 0;
	double E1,E2;
	double shift = 0.5;

	k = -1;
	l = -1;
	m = -1;
	n = -1;
	if(noneZero==2 && absSum==0) // Resonance type: kl;ll or kl;lm
	{
		for(i=0;i<nFrequencies;i++) 
		{
			if(vDiff[i]<0) k=i;
			else if (vDiff[i]>0) m=i;
		}
		if(k==-1 || m==-1) {printf("No 11 reseonance \n");return E;};
		E1 = 0.75*sqrt(vMax[k]*vMax[m])*vMax[m]*K11(vpt2KModel,k,m,m,m); // kl;ll
		E2 = 0.75*sqrt(vMax[k]*vMax[m])*vMax[k]*K11(vpt2KModel,k,k,k,m); //kk;km
		for(l=0;l<nFrequencies;l++) 
		{
			if(l!=m && l!=k) 
			{
				double E3 = 0.5*sqrt(vMax[k]*vMax[m])*(vMax[l]+shift)*K11(vpt2KModel,k,l,m,l); //kl;lm
				E += E3;
				//if(fabs(E3)>20.0) fprintf(stderr,"k=%d l=%d m=%d Eklml=%f\n",k,l,m,E3);
			}
		}
		E += E1+E2;
		if(fabs(E)>100.0 && printMax) fprintf(stderr,"k=%d m=%d Eklll = %f Ekkkm= %f Eklml=%f\n",k,m,E1,E2,E-E1-E2);
		return E;
	}
	return E;
}
*/
/**********************************************************************/
static double Eoff(VPT2KModel* vpt2KModel, int* vA, int* vB)
{
	int nFrequencies = vpt2KModel->vAnharmonic.nFrequencies;
	int* vDiff = newVectorInt(nFrequencies);
	int* vMax = newVectorInt(nFrequencies);
	double E = 0;
	int negative = 0;
	int noneZero = 0;
	int absSum = 0;
	int sumAbs = 0;
	int i;

	for(i=0;i<nFrequencies;i++) vDiff[i] = vB[i]-vA[i];
	for(i=0;i<nFrequencies;i++) absSum += vDiff[i];
	absSum = abs(absSum);
	for(i=0;i<nFrequencies;i++) sumAbs += abs(vDiff[i]);
	for(i=0;i<nFrequencies;i++) noneZero += (vDiff[i]!=0)?1:0;
	for(i=0;i<nFrequencies;i++) negative += (vDiff[i]<0)?1:0;
	for(i=0;i<nFrequencies;i++) vMax[i] = (vB[i]>vA[i])?vB[i]:vA[i];

	if(sumAbs==3) E = Eoff12(vpt2KModel, vDiff, vMax, noneZero, negative);
	if(sumAbs==4 && absSum==0) E = Eoff22(vpt2KModel, vDiff, vMax, noneZero, negative);
	//if(sumAbs==4 && absSum==2) E = Eoff13(vpt2KModel, vDiff, vMax, noneZero, negative);
	//if(sumAbs==2) E = Eoff11(vpt2KModel, vDiff, vMax, noneZero, absSum);

	freeVectorInt(&vDiff);
	freeVectorInt(&vMax);
	return E;
}
/**********************************************************************/
static void computeHMatrix(VPT2KModel* vpt2KModel)
{
	int i,j;
	VKAnharmonic* vAnharmonic = &vpt2KModel->vAnharmonic;
	int nStates = vAnharmonic->nStates;

//	printf("begin freeMatrixDouble in computeHMatrix\n");
	freeMatrixDouble(&vAnharmonic->H, nStates);
//	printf("End freeMatrixDouble in computeHMatrix\n");
	vAnharmonic->H = newMatrixDouble(nStates, nStates);
	initMatrixDouble(vAnharmonic->H, nStates, nStates, 0.0);
	for(i=0;i<nStates;i++)
	{
		vAnharmonic->H[i][i] = Ekl(vpt2KModel,vAnharmonic->v[i]);
		for(j=0;j<i;j++)
		{
			vAnharmonic->H[i][j] = Eoff(vpt2KModel,vAnharmonic->v[i],vAnharmonic->v[j]);
			//vAnharmonic->H[i][j] = rand()*1.0/RAND_MAX;
			vAnharmonic->H[j][i] = vAnharmonic->H[i][j];
		}
	}
}
/**********************************************************************/
static void diagonalizeHMatrix(VPT2KModel* vpt2KModel)
{
	VKAnharmonic* vAnharmonic = &vpt2KModel->vAnharmonic;
	int nStates = vAnharmonic->nStates;

	if(vAnharmonic->eigenVectors) freeMatrixDouble(&vAnharmonic->eigenVectors, nStates);
	vAnharmonic->eigenVectors = newMatrixDouble(nStates, nStates);
	initMatrixDouble(vAnharmonic->eigenVectors, nStates, nStates, 0.0);

	if(vAnharmonic->eigenValues) freeVectorDouble(&vAnharmonic->eigenValues);
	vAnharmonic->eigenValues = newVectorDouble(nStates);
	initVectorDouble(vAnharmonic->eigenValues, nStates, 0.0);

	eigen(vAnharmonic->H, nStates, vAnharmonic->eigenValues, vAnharmonic->eigenVectors);
}
/**********************************************************************/
static void addGroundState(VKAnharmonic* vAnharmonic)
{
	int nFrequencies = vAnharmonic->nFrequencies;
	if(vAnharmonic->v == NULL)
	{
		vAnharmonic->nStates = 1;
		vAnharmonic->v = newMatrixInt(vAnharmonic->nStates, nFrequencies);
		initMatrixInt(vAnharmonic->v, vAnharmonic->nStates, nFrequencies,0.0);
	}
	else
	{
		
		vAnharmonic->nStates += 1;
		vAnharmonic->v = realloc (vAnharmonic->v, vAnharmonic->nStates*sizeof(int*));
		{
			int ii=vAnharmonic->nStates-1;
			vAnharmonic->v[ii] = malloc(nFrequencies*sizeof(int));
		}
	}
	{
		int ii=vAnharmonic->nStates-1;
		int k;
		for(k=0;k<nFrequencies;k++)
			vAnharmonic->v[ii][k] = 0;
	}
}
/**********************************************************************/
static void addFundamentalsStates(VKAnharmonic* vAnharmonic)
{
	int nFrequencies = vAnharmonic->nFrequencies;
	int k,i;
	if(vAnharmonic->v == NULL)
	{
		vAnharmonic->nStates = nFrequencies;
		vAnharmonic->v = newMatrixInt(vAnharmonic->nStates, nFrequencies);
		initMatrixInt(vAnharmonic->v, vAnharmonic->nStates, nFrequencies,0.0);
	}
	else
	{
		
		vAnharmonic->nStates += nFrequencies;
		vAnharmonic->v = realloc (vAnharmonic->v, vAnharmonic->nStates*sizeof(int*));
		for(i=0;i<nFrequencies;i++)
		{
			int ii=i+vAnharmonic->nStates-nFrequencies;
			vAnharmonic->v[ii] = malloc(nFrequencies*sizeof(int));
		}
	}
	for(i=0;i<nFrequencies;i++)
	{
		int ii=i+vAnharmonic->nStates-nFrequencies;
		for(k=0;k<nFrequencies;k++)
			vAnharmonic->v[ii][k] = (i==k)?1:0;
	}
}
/**********************************************************************/
static void addOvertonesStates(VKAnharmonic* vAnharmonic)
{
	int nFrequencies = vAnharmonic->nFrequencies;
	int k,i;
	if(vAnharmonic->v == NULL)
	{
		vAnharmonic->nStates = nFrequencies;
		vAnharmonic->v = newMatrixInt(vAnharmonic->nStates, nFrequencies);
		initMatrixInt(vAnharmonic->v, vAnharmonic->nStates, nFrequencies,0.0);
	}
	else
	{
		vAnharmonic->nStates += nFrequencies;
		vAnharmonic->v = realloc (vAnharmonic->v, vAnharmonic->nStates*sizeof(int*));
		for(i=0;i<nFrequencies;i++)
		{
			int ii=i+vAnharmonic->nStates-nFrequencies;
			vAnharmonic->v[ii] = malloc(nFrequencies*sizeof(int));
		}
	}
	for(i=0;i<nFrequencies;i++)
	{
		int ii=i+vAnharmonic->nStates-nFrequencies;
		for(k=0;k<nFrequencies;k++)
			vAnharmonic->v[ii][k] = (i==k)?2:0;
	}
}
/**********************************************************************/
static void addCombination11States(VKAnharmonic* vAnharmonic)
{
	int nFrequencies = vAnharmonic->nFrequencies;
	int k,l,i;
	int n = 0;
	int ii;
	for(i=0;i<nFrequencies;i++)
	for(k=i+1;k<nFrequencies;k++) n++;

	if(vAnharmonic->v == NULL)
	{
		vAnharmonic->nStates = nFrequencies;
		vAnharmonic->v = newMatrixInt(vAnharmonic->nStates, nFrequencies);
		initMatrixInt(vAnharmonic->v, vAnharmonic->nStates, nFrequencies,0.0);
	}
	else
	{

		vAnharmonic->nStates += n;
		vAnharmonic->v = realloc (vAnharmonic->v, vAnharmonic->nStates*sizeof(int*));
		for(i=0;i<n;i++)
		{
			int ii=i+vAnharmonic->nStates-n;
			vAnharmonic->v[ii] = malloc(nFrequencies*sizeof(int));
			for(k=0;k<nFrequencies;k++) vAnharmonic->v[ii][k] = 0;
		}
	}
	ii = vAnharmonic->nStates-n;
	for(k=0;k<nFrequencies;k++)
	{
		for(l=k+1;l<nFrequencies;l++) 
		{
			vAnharmonic->v[ii][k] = 1;
			vAnharmonic->v[ii][l] = 1;
			ii++;
		}
	}
}
/**********************************************************************/
static void printStates(VKAnharmonic* vAnharmonic)
{
	int nFrequencies = vAnharmonic->nFrequencies;
	int nStates = vAnharmonic->nStates;
	int k,i;
	//FILE* file = NULL;
	if(nFrequencies>0)
	{
	printf("nFrequencies %d\n",nFrequencies);
	printf("nStates %d\n",nStates);
	for(i=0;i<nStates;i++)
	{
		printf("State n %d\t",i+1);
		for(k=0;k<nFrequencies;k++)
			printf("%d ",vAnharmonic->v[i][k]);
		printf("\n");
	}
	}
	/*
	file =fopen("poly.in","w");
	for(i=0;i<nStates;i++)
	{
		for(k=0;k<nFrequencies;k++)
			fprintf(file,"%d ",vAnharmonic->v[i][k]);
		fprintf(file, "\n");
	}
	fclose(file);
	*/
}
/**********************************************************************/
static VKAnharmonic newVKAnharmonic(VPT2PotentialData* vData)
{
	VKAnharmonic vAnharmonic;
	int n = 0;
	if(vData) n = vData->nFrequencies;
	vAnharmonic.nFrequencies = n;
	if(vAnharmonic.nFrequencies<=0) vAnharmonic.nFrequencies = 0;

	vAnharmonic.Xi = newMatrixDouble(vAnharmonic.nFrequencies,vAnharmonic.nFrequencies);
	initMatrixDouble(vAnharmonic.Xi, vAnharmonic.nFrequencies, vAnharmonic.nFrequencies, 0.0);
	vAnharmonic.C = 0.0;

	vAnharmonic.D = newCubeDouble(2*vAnharmonic.nFrequencies+1, 2*vAnharmonic.nFrequencies+1, 2*vAnharmonic.nFrequencies+1);
	initCubeDouble(vAnharmonic.D, 2*vAnharmonic.nFrequencies+1, 2*vAnharmonic.nFrequencies+1, 2*vAnharmonic.nFrequencies+1, 0.0);

	vAnharmonic.nStates = 0;

	vAnharmonic.v = NULL;
	addGroundState(&vAnharmonic);
	addFundamentalsStates(&vAnharmonic);
	addOvertonesStates(&vAnharmonic);
  	addCombination11States(&vAnharmonic);
	printStates(&vAnharmonic);
	vAnharmonic.H = newMatrixDouble(vAnharmonic.nStates, vAnharmonic.nStates);
	initMatrixDouble(vAnharmonic.H, vAnharmonic.nStates, vAnharmonic.nStates, 0.0);

	vAnharmonic.eigenValues = NULL;
	vAnharmonic.eigenVectors = NULL;

	vAnharmonic.harmonicEnergies = newVectorDouble(vAnharmonic.nStates);
	initVectorDouble(vAnharmonic.harmonicEnergies, vAnharmonic.nStates, 0.0);

	return vAnharmonic;
}
/*****************************************************************************/
static void readData(VPT2KModel* vpt2KModel, char* inputFileName)
{
	FILE* inputFile;
        inputFile = fopen(inputFileName,"rb");
	if(!inputFile)
	{
		fprintf(stderr, "==========================================================\n");
		fprintf(stderr, "Sorry, I cannot open the %s file\n", inputFileName);
		fprintf(stderr, "==========================================================\n");
		exit(1);
	}
	fclose(inputFile);
	vpt2KModel->vData.klass->readData(&vpt2KModel->vData, inputFileName);
	vpt2KModel->pData.klass->readData(&vpt2KModel->pData, inputFileName);
//	printf("End read readData\n");
}
/**********************************************************************/
static void printResonanceMatrix(VPT2KModel* vpt2KModel)
{
	int n = vpt2KModel->vAnharmonic.nStates;
	VKAnharmonic* vAnharmonic  = &vpt2KModel->vAnharmonic;
	int nFrequencies = vAnharmonic->nFrequencies;
        int i,j;
	int k;
	double cutoff = 1e-9;
	double** M = vAnharmonic->H;
	double E0 = M[0][0];
	int nR;
	double df;
	int* iDiff = newVectorInt(nFrequencies);
	int* jDiff = newVectorInt(nFrequencies);
	int ndi,ndj;
	// Fermi Resonances
	nR = 0;
        printf(" Fermi Resonances terms :\n");
        for(i = 0;i<n; i++)
        {
		int ni = 0; 
		for(k=0;k<nFrequencies;k++) ni += vAnharmonic->v[i][k];

                for(j = 0;j<i; j++)
                if(fabs(M[i][j])>=cutoff)
		{
			int nj = 0; 
			double X;
			for(k=0;k<nFrequencies;k++) nj += vAnharmonic->v[j][k];
                        if( !(ni==1 && nj==2) && !(ni==2 && nj==1)) continue;
			df = 0; for(k=0;k<nFrequencies;k++) df += (vAnharmonic->v[j][k]-vAnharmonic->v[i][k])*vpt2KModel->vData.hessian[k][k];
		 	X =  0.0;
			ndi = getDiff(vpt2KModel, vAnharmonic->v[i], iDiff);
			ndj = getDiff(vpt2KModel, vAnharmonic->v[j], jDiff);
			if(ndi==1 && ndj==2) X =  getMaxDiffVPT2Var(vpt2KModel, iDiff[0], jDiff[0], jDiff[1], df);
			if(ndi==2 && ndj==1) X =  getMaxDiffVPT2Var(vpt2KModel, iDiff[0], iDiff[1], jDiff[0], df);
			if(ndi==1 && ndj==1) 
			{
				if(ni==2) X =  getMaxDiffVPT2Var(vpt2KModel, iDiff[0], iDiff[0], jDiff[0], df);
				if(nj==2) X =  getMaxDiffVPT2Var(vpt2KModel, iDiff[0], jDiff[0], jDiff[0], df);
			}

			if(nR==0) printf("%20s %20s %20s %20s %20s\n", "Diff. harm.","PT2-Var","Hii","Hij","Hjj");
                        printf("%20.10f %20.10f ",df,X);
                        printf("%20.10f %20.10f %20.10f ",M[i][i]-E0,M[i][j],M[j][j]-E0);
                	for(k=0;k<nFrequencies;k++) if(vAnharmonic->v[i][k]>0) printf("%2d(%d) ",k+1,vAnharmonic->v[i][k]);
			printf(" --- ");
                	for(k=0;k<nFrequencies;k++) if(vAnharmonic->v[j][k]>0) printf("%2d(%d) ",k+1,vAnharmonic->v[j][k]);
			printf("\n");
			nR++;
		}
        }
        printf(" %d Fermi Resonances \n",nR);
        for(k=0;k<100;k++) if(k<99) printf("-"); else printf("-\n\n");
	// Darling-Dennison Resonances
	nR = 0;
        printf(" Darling-Dennison Resonances terms :\n");
        for(i = 0;i<n; i++)
        {
		int ni = 0; 
		for(k=0;k<nFrequencies;k++) ni += vAnharmonic->v[i][k];

                for(j = 0;j<i; j++)
                if(fabs(M[i][j])>=cutoff)
		{
			int nj = 0; 
			for(k=0;k<nFrequencies;k++) nj += vAnharmonic->v[j][k];
                        if (!(ni==2 && nj==2)) continue;

			df = 0; for(k=0;k<nFrequencies;k++) df += (vAnharmonic->v[j][k]-vAnharmonic->v[i][k])*vpt2KModel->vData.hessian[k][k];
			df /= 2; // 2(wi-wj); df = wi-wj
			if(nR==0) printf("%20s %20s %20s %20s\n", "Diff. harm.","Hii","Hij","Hjj");
                        printf("%20.10f ",df);
                        printf("%20.10f %20.10f %20.10f ",M[i][i]-E0,M[i][j],M[j][j]-E0);
                	for(k=0;k<nFrequencies;k++) if(vAnharmonic->v[i][k]>0) printf("%2d(%d) ",k+1,vAnharmonic->v[i][k]);
			printf(" --- ");
                	for(k=0;k<nFrequencies;k++) if(vAnharmonic->v[j][k]>0) printf("%2d(%d) ",k+1,vAnharmonic->v[j][k]);
			printf("\n");
			nR++;
		}
        }
        printf(" %d Darling-Dennison Resonances \n",nR);
        for(k=0;k<100;k++) if(k<99) printf("-"); else printf("-\n\n");
	nR = 0;
        for(i = 0;i<n; i++)
        {
		int ni = 0; 
		for(k=0;k<nFrequencies;k++) ni += vAnharmonic->v[i][k];

                for(j = 0;j<i; j++)
                if(fabs(M[i][j])>=cutoff)
		{
			int nj = 0; 
			for(k=0;k<nFrequencies;k++) nj += vAnharmonic->v[j][k];
                        if( !(ni==1 && nj==2) && !(ni==2 && nj==1)) continue;
                        if (!(ni==2 && nj==2)) continue;

			df = 0; for(k=0;k<nFrequencies;k++) df += (vAnharmonic->v[j][k]-vAnharmonic->v[i][k])*vpt2KModel->vData.hessian[k][k];
			if(nR==0) printf("%20s %20s %20s %20s\n", "Diff. harm.","Hii","Hij","Hjj");
                        printf("%20.10f ",df);
                        printf("%20.10f %20.10f %20.10f ",M[i][i]-E0,M[i][j],M[j][j]-E0);
                	for(k=0;k<nFrequencies;k++) if(vAnharmonic->v[i][k]>0) printf("%2d(%d) ",k+1,vAnharmonic->v[i][k]);
			printf(" --- ");
                	for(k=0;k<nFrequencies;k++) if(vAnharmonic->v[j][k]>0) printf("%2d(%d) ",k+1,vAnharmonic->v[j][k]);
			printf("\n");
		}
        }
        if(nR>0) printf(" %d others Resonances \n",nR);
        if(nR>0) 
	{
		for(k=0;k<100;k++) 
			if(k<99) printf("-"); 
			else printf("-\n\n");
	}
	free(iDiff);
	free(jDiff);
}
/**********************************************************************/
static void printDiagonalMatrix(VPT2KModel* vpt2KModel)
{
	int n = vpt2KModel->vAnharmonic.nStates;
	VKAnharmonic* vAnharmonic  = &vpt2KModel->vAnharmonic;
	int nFrequencies = vAnharmonic->nFrequencies;
        int i;
	int k;
	double** M = vAnharmonic->H;
	double E0 = M[0][0];
        for(i = 0;i<n; i++)
        {
                printf("%20.10f %20.10f ",M[i][i],M[i][i]-E0);
                for(k=0;k<nFrequencies;k++) if(vAnharmonic->v[i][k]>0) printf("%2d(%d) ",k+1,vAnharmonic->v[i][k]);
		printf("\n");
        }
}
/**********************************************************************/
static void computeAnharmonicEnergies(VPT2KModel* vpt2KModel)
{
	computeXi(vpt2KModel);
	computeC(vpt2KModel);
	printf("\nXkl (cm^-1) by VPT2+K model\n");
	printMatrixDoubleCutOff(vpt2KModel->vAnharmonic.Xi, vpt2KModel->vAnharmonic.nFrequencies, vpt2KModel->vAnharmonic.nFrequencies,1e-10);
	printf("END\n\n");
	printf("C=%f\n",vpt2KModel->vAnharmonic.C);
	computeHMatrix(vpt2KModel);
	/*
	printf("\nH (cm^-1) by VPT2+K model\n");
	printMatrixDoubleCutOff(vpt2KModel->vAnharmonic.H, vpt2KModel->vAnharmonic.nStates, vpt2KModel->vAnharmonic.nStates,100.0);
	printf("END\n\n");
	*/

	printf("\nDiagonal of H (cm^-1) by VPT2+K model\n");
	printDiagonalMatrix(vpt2KModel);
	printf("END\n\n");

	printf("\nResonance Matrix (cm^-1) by VPT2+K model\n");
	printf("========================================\n");
	printResonanceMatrix(vpt2KModel);
	printf("END\n\n");

	diagonalizeHMatrix(vpt2KModel);

	//printIREnergies(&vpt2KModel->vAnharmonic);
	printIREnergiesCompact2(vpt2KModel,0.0);
	//printIREnergiesCompact2(vpt2KModel,0.5);
}

/**********************************************************************/
static PropertiesKAnharmonic newPAnharmonic(VPT2PotentialData* vData, int nStates)
{
	PropertiesKAnharmonic pAnharmonic;
	int i;
	int n = 0;
	if(vData) n = vData->nFrequencies;
	pAnharmonic.nFrequencies = n;
	if(pAnharmonic.nFrequencies<=0) pAnharmonic.nFrequencies = 0;
	pAnharmonic.nStates = nStates;

	pAnharmonic.harmonic = newVectorDouble(pAnharmonic.nStates);
	initVectorDouble(pAnharmonic.harmonic, pAnharmonic.nStates, 0.0);

	pAnharmonic.anHarmonic = newVectorDouble(pAnharmonic.nStates);
	initVectorDouble(pAnharmonic.anHarmonic, pAnharmonic.nStates, 0.0);

        pAnharmonic.maxFrequencyDifference11Resonance = 2.0;
        for( i=0;i<2;i++) pAnharmonic.thresholds11Numerators[i] = 10.0;
        for( i=0;i<3;i++) pAnharmonic.parameters11Resonance[i] = 0.0;

	return pAnharmonic;
}
/**********************************************************************/
static int getDiff(VPT2KModel* vpt2KModel, int* v, int* iDiff)
{
	int nFrequencies = vpt2KModel->vAnharmonic.nFrequencies;
	int nd=0;
	int i=0;

	for(i=0;i<nFrequencies;i++) iDiff[i] = -1;
	for(i=0;i<nFrequencies;i++) {if(v[i]>0) iDiff[nd++] = i;};
	return nd;
}
/**********************************************************************/
static void computeHarmonicValues(VPT2KModel* vpt2KModel)
{
        int i,j;
        int k;
        int a;
        int nd;
	int jmax;
        VKAnharmonic* vAnharmonic  = &vpt2KModel->vAnharmonic;
	PropertiesKAnharmonic* pAnharmonic = &vpt2KModel->pAnharmonic;
        int nFrequencies = vAnharmonic->nFrequencies;
        int nStates = pAnharmonic->nStates;
	int nDim = vpt2KModel->pData.nDim;
	double* V = newVectorDouble(nDim);
	double mu0 = 4*M_PI*1e-7;
	double eps0 = 1.0/(mu0*slight*slight);
	double   kmmolm1 = 8*M_PI*M_PI*M_PI*NAvogadro/3/hPlank/slight/4/M_PI/eps0*1e-3*100.0*8.47835267e-30*8.47835267e-30;/* 1e-3 m to km, 100 : cm-1 to m-1 */
	int* iDiff = newVectorInt(nFrequencies);

	for(i=0;i<nStates;i++)
	{
		double E = 0;
		int* v = vpt2KModel->vAnharmonic.v[i];
		k = 0;
		jmax = 0;
		for(j=0;j<nStates;j++)
			if(fabs(vAnharmonic->eigenVectors[i][j])>fabs(vAnharmonic->eigenVectors[i][jmax])) jmax = j;
		nd = getDiff(vpt2KModel, v, iDiff);
		if(nd==1 && v[iDiff[k]]==1) 
		{
			computeHarmonicProp(vpt2KModel, iDiff[k], V);
			pAnharmonic->harmonic[jmax] = 0.0;
			for(a=0;a<nDim;a++) pAnharmonic->harmonic[jmax] += V[a]*V[a];
			pAnharmonic->harmonic[jmax] *= kmmolm1;
		}
		for(k=0;k<nFrequencies;k++)
			E += v[k]*vpt2KModel->vData.hessian[k][k];
		vAnharmonic->harmonicEnergies[jmax] = E;
		//printf("i=%d jmax = %d E=%f\n",i,jmax,E);
	}
	for(j=0;j<nStates;j++)
	{
		int imax = 0;
		for(i=0;i<nStates;i++)
		{
			if(fabs(vAnharmonic->eigenVectors[i][j])>fabs(vAnharmonic->eigenVectors[imax][j])) imax = i;
		}
		int* v = vpt2KModel->vAnharmonic.v[imax];
		double E = 0;
		for(k=0;k<nFrequencies;k++)
			E += v[k]*vpt2KModel->vData.hessian[k][k];
		vAnharmonic->harmonicEnergies[j] = E;
	}

	free(iDiff);
	free(V);
}
/**********************************************************************/
static void computeAnharmonicProperties(VPT2KModel* vpt2KModel)
{
	int i,j;
	int k;
	int a;
	int nd;
	int nFrequencies = vpt2KModel->vAnharmonic.nFrequencies;
	int nStates = vpt2KModel->pAnharmonic.nStates;
	int* iDiff = newVectorInt(nFrequencies);
	int nDim = vpt2KModel->pData.nDim;
	double** P = newMatrixDouble(nStates, nDim);
	double* V = newVectorDouble(nDim);
	int* v;
	VKAnharmonic* vAnharmonic = &vpt2KModel->vAnharmonic;
	PropertiesKAnharmonic* pAnharmonic = &vpt2KModel->pAnharmonic;
	double mu0 = 4*M_PI*1e-7;
	double eps0 = 1.0/(mu0*slight*slight);
	double   kmmolm1 = 8*M_PI*M_PI*M_PI*NAvogadro/3/hPlank/slight/4/M_PI/eps0*1e-3*100.0*8.47835267e-30*8.47835267e-30;/* 1e-3 m to km, 100 : cm-1 to m-1 */

	initMatrixDouble(P, nStates, nDim, 0.0);
	for(i=0;i<nStates;i++)
	{
		v = vpt2KModel->vAnharmonic.v[i];
		nd = getDiff(vpt2KModel, v, iDiff);
		if(printMax)
		{
			for(k=0;k<nFrequencies;k++) printf("%2d ",v[k]);
			printf(" : ");
			if(nd==0) { printf("\n");continue;}
			for(k=0;k<nd;k++)
				printf(" %2d(%d) ", iDiff[k]+1, v[iDiff[k]]);
		}
		if(nd==1) {
			int k = 0;
			// fund or overtones
			if(v[iDiff[k]]==1) 
			{
				computeFundamentalsProp(vpt2KModel, iDiff[k], P[i]);
				if(printMax)
				{
					printf(" fundamental(%d) ",iDiff[k]+1);
					for(a=0;a<nDim;a++) printf(" %f ",P[i][a]);
					printf("\n");
				}
			}
			else if(v[iDiff[k]]==2) 
			{
				computeOvertonesProp(vpt2KModel, iDiff[k], P[i]);
				if(printMax)
				{
					printf(" overtones(%d) ",iDiff[k]+1);
					for(a=0;a<nDim;a++) printf(" %f ",P[i][a]);
					printf("\n");
				}
			}
			else 
			{
				printf(" not yet implemented(%d)\n",iDiff[k]+1);
			}
		
		}
		if(nd==2) {
		// combination
			int k = 0; 
			int l = 1;
			if(v[iDiff[k]]==1 && v[iDiff[l]]==1)
			{
				computeCombinationBandsProp(vpt2KModel,iDiff[k],iDiff[l], P[i]);
				if(printMax)
				{
					printf(" Combination 1-1(%d,%d) ",iDiff[k]+1,iDiff[l]+1);
					for(a=0;a<nDim;a++) printf(" %f ",P[i][a]);
					printf("\n");
				}
			}
			else 
			{
				printf(" not yet implemented(%d,%d)\n",iDiff[k]+1,iDiff[l]+1);
			}
		}
		if(nd>2)
		printf("nd=%d\n",nd);
	}
	for(j=0;j<nStates;j++)
	{
		double intensity=0;
		for(a=0;a<nDim;a++) V[a] = 0.0;
		if(vpt2KModel->pData.model.type==MODEL_PROP_GVPT2) 
		{
			int imax = 0;
			for(i=1;i<nStates;i++) if(fabs(vAnharmonic->eigenVectors[i][j])>fabs(vAnharmonic->eigenVectors[imax][j])) imax = i;
			//for(a=0;a<nDim;a++) V[a] = vAnharmonic->eigenVectors[imax][j]*P[imax][a];
			for(a=0;a<nDim;a++) V[a] = P[imax][a];
		}
		else
		{
			for(i=0;i<nStates;i++) for(a=0;a<nDim;a++) 
				V[a] += vAnharmonic->eigenVectors[i][j]*P[i][a];
		}
		for(a=0;a<nDim;a++) intensity += V[a]*V[a];
		intensity *= kmmolm1*(vAnharmonic->eigenValues[j]-vAnharmonic->eigenValues[0]);
		pAnharmonic->anHarmonic[j] = intensity;
	}
	freeMatrixDouble(&P,nStates);
	free(V);
	computeHarmonicValues(vpt2KModel);
	printIREnergiesPropCompact2(vpt2KModel,0.0);
	printIREnergiesPropCompact(vpt2KModel,1e-4);
	free(iDiff);
}
/**********************************************************************/
VPT2KModel newVPT2KModel()
{
	VPT2KModel vpt2KModel;
	vpt2KModel.klass = malloc(sizeof(VPT2KModelClass));
	vpt2KModel.klass->readData = readData;
	vpt2KModel.klass->computeAnharmonic =computeAnharmonic;
	vpt2KModel.vData = newVPT2PotentialData(0);
	vpt2KModel.pData = newVPT2PropertiesData(0,0);
	vpt2KModel.vAnharmonic = newVKAnharmonic(&vpt2KModel.vData);
	vpt2KModel.pAnharmonic = newPAnharmonic(&vpt2KModel.vData,vpt2KModel.vAnharmonic.nStates);
	initD(&vpt2KModel);
//	printf("End newVPT2KModel\n");
	return vpt2KModel;
}
/**********************************************************************/
static void computeAnharmonic(VPT2KModel* vpt2KModel)
{
	if(!vpt2KModel) return;
//	printf("Begin newVKAnharmonic\n");
	vpt2KModel->vAnharmonic = newVKAnharmonic(&vpt2KModel->vData);
	vpt2KModel->pAnharmonic = newPAnharmonic(&vpt2KModel->vData, vpt2KModel->vAnharmonic.nStates);
	initD(vpt2KModel);
//	printf("Begin computeAnharmonicEnergies\n");
	computeAnharmonicEnergies(vpt2KModel);
//	printf("End computeAnharmonicEnergies\n");
	computeAnharmonicProperties(vpt2KModel);
}
/**********************************************************************/
static boolean testDegen(VPT2KModel* vpt2KModel, double w1, double w2)
{
        double maxFrequencyDifference11Resonance = vpt2KModel->pAnharmonic.maxFrequencyDifference11Resonance;
	if(fabs(w2-w1)<maxFrequencyDifference11Resonance) return TRUE;
	return FALSE;
}
/**********************************************************************/
static boolean testResonanceFreq(VPT2KModel* vpt2KModel, double w1, double w2)
{
	double* parameters11Resonance = vpt2KModel->pAnharmonic.parameters11Resonance;
        double maxFrequencyDifference11Resonance = vpt2KModel->pAnharmonic.maxFrequencyDifference11Resonance;
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
static boolean testResonanceCubic(VPT2KModel* vpt2KModel, double w1, double w2, double v)
{
	if(testResonanceFreq(vpt2KModel,w1,w2) && fabs(v)>= vpt2KModel->pAnharmonic.thresholds11Numerators[0]) return TRUE;
	return FALSE;
}
/**********************************************************************/
static boolean testResonanceQuartic(VPT2KModel* vpt2KModel, double w1, double w2, double v)
{
	if(testResonanceFreq(vpt2KModel,w1,w2) && fabs(v)>= vpt2KModel->pAnharmonic.thresholds11Numerators[1]) return TRUE;
	return FALSE;
}
/**********************************************************************/
static double VPT2K(VPT2KModel* vpt2KModel, double k2, int ii, int jj, int kk)
{
	double vpt2k = 0.0;
	if(fabs(k2)>1e-10) vpt2k = k2*D(vpt2KModel,ii,jj,kk);
	//printf("ii=%d jj=%d kk=%d D=%f\n",ii,jj,kk,D(vpt2KModel,ii,jj,kk));
	/*
	VPT2PotentialData* vData = &vpt2KModel->vData;
	VPT2PropModel* pModel = &vpt2KModel->pData.model;
	double wi=(ii>0)?vData->hessian[ii-1][ii-1]:-vData->hessian[-ii-1][-ii-1];
	double wj=(jj>0)?vData->hessian[jj-1][jj-1]:-vData->hessian[-jj-1][-jj-1];
	double wk=(kk>0)?vData->hessian[kk-1][kk-1]:-vData->hessian[-kk-1][-kk-1];
	double epsilon = wi+wj+wk;
	double vpt2k = k2/epsilon;
	*/
	return vpt2k;
}
/**********************************************************************/
static double VPT2(double k2, double epsilon)
{
	double vpt2 = 0.0;
	if(fabs(k2)>1e-10) vpt2 = k2/epsilon;
	return vpt2;
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
/*
static double applyModelProp(double k2, double epsilon, VPT2PropModel* pModel)
{
	if(pModel->type==MODEL_PROP_DCPT2) return DCPT2(k2,epsilon);
	if(pModel->type==MODEL_PROP_HDCPT2) return HDCPT2(k2,epsilon,pModel->alphaHDCPT2,pModel->betaHDCPT2);
	return VPT2(k2,epsilon);
}
*/
/**********************************************************************/
static double applyModelPropNew(VPT2KModel* vpt2KModel, double k2, int ii, int jj, int kk)
{
	VPT2PotentialData* vData = &vpt2KModel->vData;
	VPT2PropModel* pModel = &vpt2KModel->pData.model;
	double wi=(ii>0)?vData->hessian[ii-1][ii-1]:-vData->hessian[-ii-1][-ii-1];
	double wj=(jj>0)?vData->hessian[jj-1][jj-1]:-vData->hessian[-jj-1][-jj-1];
	double wk=(kk>0)?vData->hessian[kk-1][kk-1]:-vData->hessian[-kk-1][-kk-1];
	double epsilon = wi+wj+wk;
	if(pModel->type==MODEL_PROP_DCPT2) return DCPT2(k2,epsilon);
	if(pModel->type==MODEL_PROP_HDCPT2) return HDCPT2(k2,epsilon,pModel->alphaHDCPT2,pModel->betaHDCPT2);
	//if(pModel->type==MODEL_PROP_VPT2) {fprintf(stderr,"VPT2 PROP\n");return VPT2(k2,epsilon);}
	if(pModel->type==MODEL_PROP_VPT2) return VPT2(k2,epsilon);
	return VPT2K(vpt2KModel,k2,ii,jj,kk);
}
/**********************************************************************/
/**********************************************************************/
/* PCCP, 2014, 16, 1759-1787, page 1763-4 */
/* CPL, 496 (2010) 157–161 */
/**********************************************************************/
static void computeFundamentalsProp(VPT2KModel* vpt2KModel, int i, double* Pav)
{
	int j,k,l;
	int a,c;
	double s = 0;
	double* V = NULL;
	double mu0 = 4*M_PI*1e-7;
	double eps0 = 1.0/(mu0*slight*slight);
	//double MWQ2q  = hPlank/4/M_PI/slight;
	double   kmmolm1 = 8*M_PI*M_PI*M_PI*NAvogadro/3/hPlank/slight/4/M_PI/eps0*1e-3*100.0*8.47835267e-30*8.47835267e-30;/* 1e-3 m to km, 100 : cm-1 to m-1 */
	double s0 = 1.0/sqrt(2.0);
	double s1 = s0/2;
	double s2 = s0/6;
	double S = 1.0;
	//VPT2PropModel* pModel = NULL;
	VPT2PotentialData vData;
	VPT2PropertiesData pData;
	//VKAnharmonic vAnharmonic;
	//PropertiesKAnharmonic pAnharmonic;

	if(!vpt2KModel) return;
	if(printMax) printf("kmmolm1=%f\n",kmmolm1);
	pData = vpt2KModel->pData;
	vData = vpt2KModel->vData;
	//pModel = &pData.model;
	//pAnharmonic = vpt2KModel->pAnharmonic;
	//vAnharmonic = vpt2KModel->vAnharmonic;

	V = malloc(pData.nDim*sizeof(double));

	{
		double wi=vData.hessian[i][i];
		double si=1/(wi);
		double sqrti=(wi>0)?sqrt(1.0/wi):0;

		// HARMONIC
		for(a=0;a<pData.nDim;a++) Pav[a] = 0;
		for(a=0;a<pData.nDim;a++) V[a] = 0;
		for(a=0;a<pData.nDim;a++)
		{
			double Pi = pData.first[a][i]*sqrti;
			V[a] = s0*S*Pi;
		}
		if(printMax){
		printf("%s %d %s %14.8f %20s","Mode = ",i+1,"Harmonic Frequency = ",wi," Harmonic term ");
		for(a=0;a<pData.nDim;a++) printf("%14.8f ",V[a]);
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
		printf("%s %d %s %14.8f %20s","Mode = ",i+1,"Harmonic Frequency = ",wi," Electric term ");
		for(a=0;a<pData.nDim;a++) printf("%14.8f ",V[a]);
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
				//boolean resDD = FALSE;
				//if(j!=i && testDD(vpt2KModel, i+1, i+1, -(j+1), -(j+1))) { resDD=TRUE;}
				/*
				if(!resDD)
				for(k=0;k<vData.nFrequencies;k++) 
				{
					if(k!=i && k!=j && testDD(vpt2KModel, i+1, i+1, -(k+1), -(k+1)) && testDD(vpt2KModel, j+1, j+1, -(k+1), -(k+1))) { resDD=TRUE;break;}
				}
				*/
				for(k=0;k<vData.nFrequencies;k++) 
				{
					//double wk=vData.hessian[k][k];
					double kijkk = getMaxQuarticIJKL(vData.quartic, i, j, k, k);
					s += kijkk*A;
					// ALLOUCHE Changed remove resDD. gaussian use resDD here
					//if(!resDD && !testDegen(vpt2KModel, wi,wj) && !testResonanceQuartic(vpt2KModel,wi,wj,kijkk)) s+= kijkk*B;
					if(!testDegen(vpt2KModel, wi,wj) && !testResonanceQuartic(vpt2KModel,wi,wj,kijkk)) s+= kijkk*B;
				}
			}
			V[a] += -s*s0/8;
		}
		if(printMax){
		printf("%s %d %s %14.8f %20s","Mode = ",i+1,"Harmonic Frequency = ",wi," Quartic term ");
		for(a=0;a<pData.nDim;a++) printf("%14.8f ",V[a]);
		printf("\n");
		}

		for(a=0;a<pData.nDim;a++) Pav[a] += V[a];
		for(a=0;a<pData.nDim;a++) V[a] = 0;

		int nnRes = 0;
		int nT1 = 0;
		double TT1 = 0.0;
		double TTT1 = 0.0;
		for(a=0;a<pData.nDim;a++)
		{
			s = 0;
			for(j=0;j<vData.nFrequencies;j++) 
			{
				double wj=vData.hessian[j][j];
				double sqrtij=(wi*wj>0)?sqrt(1.0/(wi*wj)):0;
				//double Pji = getMaxMatrixIJ(pData.second[a],j,i)*sqrtij;
				//double Pij = getMaxMatrixIJ(pData.second[a],i,j)*sqrtij;
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

					if(fabs(kijk)>1e-6)
					{
					s += kijk*(Pjk+Pkj)*(1/(wi+wj+wk));
					//s += -(Pjk+Pkj)*S*applyModelProp(kijk,wi-wj-wk,pModel);
					//if(!testFermi(vpt2KModel, i+1,-(j+1),-(k+1))) { s += -(Pjk+Pkj)*kijk/(wi-wj-wk);}
					s += -(Pjk+Pkj)*S*applyModelPropNew(vpt2KModel, kijk, i+1,-(j+1),-(k+1));
					if(testFermi(vpt2KModel, i+1,-(j+1),-(k+1))) nnRes++;
					if(!testFermi(vpt2KModel, i+1,-(j+1),-(k+1))) { nT1++; if(a==2) TT1 += (Pjk+Pkj)*kijk/(wi-wj-wk); TTT1 += kijk/(wi-wj-wk);}
					//printf("%14.8f %14.8f\n",wi, (Pjk+Pkj)*S*applyModelPropNew(vpt2KModel, kijk, i+1,-(j+1),-(k+1))*s1/8);
					//s += -(Pjk+Pkj)*S*kijk/(wi-wj-wk);
					}

					s += kjkk/wj*(2*S*Pji+(1+S)*Pij);
				}
			}
			V[a] += -s*s1/8;
		}
		if(printMax){
		printf("nnRes = %d\n",nnRes/pData.nDim);
		printf("nT1 = %d\n",nT1/pData.nDim);
		printf("TT1 = %14.8e\n",TT1);
		printf("TTT1 = %14.8e\n",TTT1/pData.nDim);
		printf("%s %d %s %14.8f %20s","Mode = ",i+1,"Harmonic Frequency = ",wi," Mixed term ");
		for(a=0;a<pData.nDim;a++) printf("%14.8f ",V[a]);
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
				double simj=(i==j || fabs(wi-wj)<1e-13)?0:1/(wi-wj);
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
					if(!testDegen(vpt2KModel, wi,wj) && !testResonanceQuartic(vpt2KModel,wi,wj,sum)) 
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
		printf("%s %d %s %14.8f %20s","Mode = ",i+1,"Harmonic Frequency = ",wi," Coriolis term ");
		for(a=0;a<pData.nDim;a++) printf("%14.8f ",V[a]);
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
				double simj=(i==j || fabs(wi-wj)<1e-13)?0:1/(wi-wj);

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
					        A  += -applyModelPropNew(vpt2KModel, kABC*dijk*sjkl, i+1,-(k+1),-(l+1));
						B = kABC*dij*(1+dik)*ddil*(0.5*si*sikl+S*0.5*sikl*sikl);
					        //t = applyModelProp(sqrt(kABC*dij*(1+dik)*ddil*0.5),wi-wk-wl,pModel);
					        t = applyModelPropNew(vpt2KModel, sqrt(kABC*dij*(1+dik)*ddil*0.5), i+1,-(k+1),-(l+1));
					        B  += -S*t*t;
					        //B  += -applyModelProp(kABC*dij*(1+dik)*ddil*0.5*si,wi-wk-wl,pModel);
					        B  += -applyModelPropNew(vpt2KModel, kABC*dij*(1+dik)*ddil*0.5*si, i+1,-(k+1),-(l+1));

						C = kABC*ddij*ddik*dil*(sk*sij+2.0*siik*sij+3.0*sij*sijk+2*S*siik*sijk+3*sk*sijk);
					        //C  += -S*applyModelProp(kABC*ddij*ddik*dil*sk,wi-wj-wk,pModel);
					        C  += -S*applyModelPropNew(vpt2KModel, kABC*ddij*ddik*dil*sk,i+1,-(j+1),-(k+1));

						if(!testDegen(vpt2KModel, wi,wj) && !testResonanceCubic(vpt2KModel,wi,wj,kikl*kjkl))
						{ 
							A += kABC*ddik*ddil*S*simj*(-sjkl);
					        	//A  += S*applyModelProp(kABC*ddik*ddil*simj,wi-wk-wl,pModel);
					        	A  += S*applyModelPropNew(vpt2KModel,kABC*ddik*ddil*simj, i+1,-(k+1),-(l+1));
							C += kABC*ddik*dil*S*simj*(-2*sijk-3*sk);
					        	//C  += S*applyModelProp( kABC*ddik*dil*simj,wi-wj-wk,pModel);
					        	C  += S*applyModelPropNew(vpt2KModel, kABC*ddik*dil*simj,i+1,-(j+1),-(k+1));
						}

						D  = kDEF*dij*si*sk*(1.0+(6.0-4.0*S)*dik*dil/9.0);
						E  = kDEF*dijk*(sij*sijk+sk*sij+sk*sijk);
					        //E  += -S*applyModelProp(kDEF*dijk*sk,wi-wj-wk,pModel);
					        E  += -S*applyModelPropNew(vpt2KModel, kDEF*dijk*sk,i+1,-(j+1),-(k+1));
						F  = kDEF*dik*ddij*(1.0+dil)*(siij*sij+si*siij);
						F += kDEF*dik*ddij*(dil)*(1.0/3.0*si*sij+S/3.0*si*siij);
						F += kDEF*dik*ddij*(si*sij+S*si*sj);

						if(!testDegen(vpt2KModel, wi,wj) && !testResonanceCubic(vpt2KModel,wi,wj,kijk*kllk))
						{ 
							E += kDEF*ddik*ddil*S*simj*(-sk);
					        	//E += S*applyModelProp(kDEF*ddik*ddil*simj,wi-wj-wk,pModel);
					        	E += S*applyModelPropNew(vpt2KModel, kDEF*ddik*ddil*simj,i+1,-(j+1),-(k+1));
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
		printf("%s %d %s %14.8f %20s","Mode = ",i+1,"Harmonic Frequency = ",wi," Cubic term ");
		for(a=0;a<pData.nDim;a++) printf("%14.8f ",V[a]);
		printf("\n");
		}

		for(a=0;a<pData.nDim;a++) Pav[a] += V[a];
		for(a=0;a<pData.nDim;a++) V[a] = 0;

		if(printMax){
		printf("%s %d %s %14.8f %20s","Mode = ",i+1,"Harmonic Frequency = ",wi," SUM ");
		for(a=0;a<pData.nDim;a++) printf("%14.8f ",Pav[a]);
		printf("\n");
		}
	}
	free(V);
}
/**********************************************************************/
/* JCP 136, 124108, 2012 Bloino & Barone */
/* JPCA 2015, Bloino DOI: 10.1021/jp509985u */
static void computeOvertonesProp(VPT2KModel* vpt2KModel, int i, double* Pav)
{
	int k;
	int a;
	double s;
	double s0 = 1.0/sqrt(2.0);
	double s1 = s0/2;
	//double s2 = s0/6;
	double S = 1;
	//VPT2PropModel* pModel = NULL;
	VPT2PotentialData vData;
	VPT2PropertiesData pData;
	//VKAnharmonic vAnharmonic;
	//PropertiesKAnharmonic pAnharmonic;

	if(!vpt2KModel) return;

	pData = vpt2KModel->pData;
	vData = vpt2KModel->vData;
	//pModel = &pData.model;
	//pAnharmonic = vpt2KModel->pAnharmonic;
	//vAnharmonic = vpt2KModel->vAnharmonic;
	{
		double wi=vData.hessian[i][i];

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
				s += S*applyModelPropNew(vpt2KModel, kiik,i+1,i+1,-(k+1))*Pk;
				s += -kiik*Pk/(wi+wi+wk);
			}
			Pav[a] += s*s0/4;
		}
	}
}
/**********************************************************************/
static void computeCombinationBandsProp(VPT2KModel* vpt2KModel, int i, int j, double* Pav)
{
	int k;
	int a;
	double s;
	double s0 = 1.0/sqrt(2.0);
	double s1 = s0/2;
	//double s2 = s0/6;
	double S = 1;
	double sqrt2 = sqrt(2.0);
	//VPT2PropModel* pModel = NULL;
	VPT2PotentialData vData;
	VPT2PropertiesData pData;
	//VKAnharmonic vAnharmonic;
	//PropertiesKAnharmonic pAnharmonic;

	if(!vpt2KModel) return;
	//pAnharmonic = vpt2KModel->pAnharmonic;
	//vAnharmonic = vpt2KModel->vAnharmonic;

	pData = vpt2KModel->pData;
	vData = vpt2KModel->vData;
	//pModel = &pData.model;

	{
		double wi=vData.hessian[i][i];
		{
			double wj=vData.hessian[j][j];
			for(a=0;a<pData.nDim;a++) Pav[a] = 0;
			if(i==j) return;
			for(a=0;a<pData.nDim;a++)
			{
				double Pij=(pData.second[a][i][j]+pData.second[a][j][i])/2.0/sqrt(wi*wj);
				//double Pij=getMaxMatrixIJ(pData.second[a],i,j)/sqrt(wi*wj);
				Pav[a] += s1*S*Pij;
				s = 0;
				for(k=0;k<vData.nFrequencies;k++) 
				{
					double wk=vData.hessian[k][k];
					double Pk = pData.first[a][k]/sqrt(wk);
					double kijk = getMaxCubicIJK(vData.cubic, i, j, k);
					//s += S*applyModelProp(kijk,wi+wj-wk,pModel)*Pk;
					s += S*applyModelPropNew(vpt2KModel,kijk,i+1,j+1,-(k+1))*Pk;
					s += -kijk*Pk/(wi+wj+wk);
				}
				Pav[a] += s*s0/4;
			}
			for(a=0;a<pData.nDim;a++) Pav[a] *= sqrt2;
		}
	}
}
/**********************************************************************/
static void computeHarmonicProp(VPT2KModel* vpt2KModel, int i, double* Pav)
{
	VPT2PropertiesData* pData = &vpt2KModel->pData;
	int a;
	for(a=0;a<pData->nDim;a++) Pav[a] = pData->first[a][i]/sqrt(2.0);
}
