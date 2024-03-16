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

/* Molecule.c */
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <ctype.h>

#include "../Utils/Utils.h"
#include "../Utils/AtomsProp.h"
#include "../Utils/Constants.h"
#include "../Utils/Types.h"
#include "../Utils/QL.h"
#include "../Molecule/Molecule.h"
#include "../Utils/CalculTypesAmber.h"
#include "../Utils/PDBTemplate.h"

static boolean** bondedMatrix = NULL;
static void removeTranslation(Molecule* molecule);
static void removeRotation(Molecule* molecule);
static void removeTranslationAndRotation(Molecule* molecule);
static void removeTranslationCluster(Molecule** molecule, int nMols);
static void removeRotationCluster(Molecule** molecules, int nMols);
static void removeTranslationAndRotationCluster(Molecule** molecules, int nMols);

static void removeTranslationForceCluster(Molecule** molecule, int nMols, double** vectors);
static void removeRotationForceCluster(Molecule** molecules, int nMols, double** vectors);
static void removeTranslationAndRotationForceCluster(Molecule** molecules, int nMols, double** vectors);

static void removeTranslationForce(Molecule* molecule, double* f);
static void removeRotationForce(Molecule* molecule, double* f);
static void removeTranslationAndRotationForce(Molecule* molecules, double* f);

static void removeTranslationMoments(Molecule** molecule, int nMols, double*** P);
static void removeRotationMoments(Molecule** molecules, int nMols, double*** P);
static void removeTranslationAndRotationMoments(Molecule** molecules, int nMols, double*** P);

static void removeTranslationAcceleration(Molecule* molecule, double* a);
static void removeRotationAcceleration(Molecule* molecule, double* a);
static void removeTranslationAndRotationAcceleration(Molecule* molecule, double* a);
static boolean printMolecule(Molecule* molecule, FILE* file);
static void scaleVelocities(Molecule* molecule, double temperature);
static double getKelvin(Molecule* molecule);
static void setMaxwellVelocities(Molecule* molecule, double temperature);
static boolean setMaxwellVelocitiesIfNull(Molecule* molecule, double temperature);
static boolean resetConstraints(Molecule* molecule, Constraints constraints);
static void setRattleConstraintsParameters(Molecule* m);
static boolean setRandomPositions(Molecule* molecule);
static boolean setRandomPositionsChain(Molecule* molecule);
static boolean addGeometry(Molecule* molecule,FILE* file);
static boolean addGeometryToGabedit(Molecule* molecule,FILE* file);
static boolean addMolecule(Molecule* molecule,FILE* file);
static boolean addVelocities(Molecule* molecule,FILE* file);
static Molecule* readGeom(char* namefile);
static boolean readGeometry(Molecule* molecule,char* namefile);
static void computeDipole(Molecule* mol);
static void copyChargeInCharge0(Molecule* molecule);
static void setBondHardness(Molecule* molecule, int nBonds, char** atomTypes1, char** atomTypes2, double* hardness);
static void setHardness(Molecule* molecule, int nTypes, char** atomTypes, double* hardness);
static void setElectronegativity(Molecule* molecule, int nTypes, char** atomTypes, double* electronegativity);
static void setWidth(Molecule* molecule, int nTypes, char** atomTypes, double* width);
static void setCharge0(Molecule* molecule, int nTypes, char** atomTypes, double* charge0);
static void setChargesEEM(Molecule* molecule);
static void setChargesACKS2(Molecule* molecule);
static double getEnergyEEM(Molecule* molecule);
static double getEnergyACKS2(Molecule* molecule);
static boolean saveGeometry(Molecule* molecule,char* fileName);
static boolean saveGeometryAndVelocities(Molecule* molecule,char* fileName);
static boolean isLinear(Molecule* molecule);
static boolean saveMolecule(Molecule* molecule,char* fileName);
static boolean saveMol2(Molecule* molecule,char* fileName);
static boolean saveHIN(Molecule* molecule,char* fileName);
static boolean saveFrequencies(Molecule* molecule, char* fileName, int nModes, double* frequencies, double** modes, double* reducedMasses, double* IRIntensities);
static int computeIR(Molecule* mol, double *F, double* dmuX[3]);
static boolean readGradientFromGabeditFile(Molecule* mol, char* namefile);
static int computeFrequenciesFromFiles(Molecule* mol, char* inputFileName, double dx);
static int computeFrequenciesOneStepFromFiles(Molecule* mol, char* inputFileName, double dx);
static int computeFrequenciesFromGradFiles(Molecule* mol, char* inputFileName, double dx);
static int computeFrequenciesOneStepFromGradFiles(Molecule* mol, char* inputFileName, double dx);
static int generateCChemIFilesForFrequencies(Molecule* mol, char* inputFileName, double dx);
static int generateCChemIFilesOneStepForFrequencies(Molecule* mol, char* inputFileName, double dx);
static int generateCChemIGradFilesForFrequencies(Molecule* mol, char* inputFileName, double dx);
static int generateCChemIGradFilesOneStepForFrequencies(Molecule* mol, char* inputFileName, double dx);
static void freeVibrations(Molecule* mol);

static void removeTransRotModes(Molecule* mol);
static void sortFrequencies(Molecule* mol);
static void removeFrequencies(Molecule* mol, double freqMin, double freqMax);
static void removeNotSelectedFrequencies(Molecule* mol, boolean* selected);
static void freeMolecule(Molecule* molecule);
static Molecule* copyMolecule(Molecule* m);
static void computeHarmonicVelocitiesCoordinates(Molecule* mol, double T, int numMode, boolean changeGeom);
static int generateQFFCChemIFiles(Molecule* mol, char* inputFileName, double delta, boolean reducedCoordinates, int ordre);
static double*** getGeomsQFF(Molecule* mol, char* inputFileName, double delta, boolean reducedCoordinates, int ordre, double** pDeltas, int* nGeoms);
static void readGeomFromMopacOutputFile(Molecule* mol, char *fileName, int numgeometry);
static void readGeomFromGamessOutputFile(Molecule* mol, char *fileName, int numgeometry);
static void readGeomFromOrcaOutputFile(Molecule* mol, char* fileName, int numgeometry);
static void readGeomFromOpenBabelOutputFile(Molecule* mol, char* fileName, int numgeometry);
static void readGeomFromGaussianOutputFile(Molecule* mol, char *fileName, int numgeometry);
static void readGeomFromMopacAuxFile(Molecule* mol, char *fileName,int numgeometry);
static void setConnections(Molecule* molecule);
static boolean resetMMTypes(Molecule* mol, char* type);
static boolean buildMMTypes(Molecule* mol, FILE* file);

static boolean addFirstDerivativeToFile(Molecule* molecule, FILE* file);
static Molecule* getMoleculeFromGaussianOutputFile(char *fileName, int numgeometry);
static boolean readVibrationFromGaussianOutputFile(Molecule* mol, char* fileName);
static boolean readDipolesDerivativesFromOrcaHessianFile(Molecule* mol, char* fileName);
static boolean readDipolesDerivativesFromGaussianOutputFile(Molecule* mol, char* fileName);
static boolean readDipolesDerivativesFromGabeditFile(Molecule* mol, char* fileName);
static boolean readDipolesDerivativesFromOrcaOutputFile(Molecule* mol, char* fileName);

static boolean readMolecule0(Molecule* molecule, char* fileName);

static boolean fit2Molecule(Molecule* molToFit, Molecule* molRef, double u[3][3]);
static boolean fitRMSD2Molecule(Molecule* molFit, Molecule* molRef, boolean center);
static void buildStandardOrientation(Molecule* mol,double* centerOfGravity, int* numberOfEquivalentAxes, double* inertialMoment, double axes[3][3]);
static void centerMolecule(Molecule* mol, double C[]);
static char* saveFirstDerivatives(char* inputFileName, Molecule* mol);
static double* getDeltaTable(Molecule* mol, double delta, boolean reducedCoordinates);
static int getQFFOneGeom(Molecule* mol, int mode1, int mode2, double delta1, double delta2, double akOverI1, double** geom);
static int getQFFOneGeom3(Molecule* mol, int mode1, int mode2, int mode3, double delta1, double delta2, double delta3, double** geom);
static int getQFFOneGeom4(Molecule* mol, int mode1, int mode2, int mode3, int mode4, double delta1, double delta2, double delta3, double delta4, double** geom);
static Molecule** readMoleculesFromXYZFile(char *fileName, boolean connections);
static Molecule** readMoleculesFromGabeditFile(char* namefile, boolean connections);
static Molecule** readMoleculesFromCChemIFile(char *fileName, boolean connections);
static void getFileNameToRead(char* fileName, char* fileNameToRead);
static void computePseudoInertia(Molecule* mol, double*pI, double* pI4);
static  boolean similarInertia(Molecule* mol1, Molecule* mol2, double precision);
static  boolean similarBonds(Molecule* mol1, Molecule* mol2, double sTol, double distMaxTol);
static boolean setGeometryToAxes(Molecule* molecule, double axis1[], double axis2[], double axis3[]);
static boolean setRandomOrientation(Molecule* molecule);
static boolean makeLocalSphericalMutation(Molecule* molecule, double rate);
static boolean makeCenterOfMassSphericalMutation(Molecule* molecule);
static  Molecule* makeSphereCutSpliceCrossover(Molecule* mol1, Molecule* mol2, int *err);
static  Molecule* makePlaneCutSpliceCrossover(Molecule* mol1, Molecule* mol2, int *err);
static  void freeTypesOfAtoms(char** types, int nTypes);
static  char** getTypesOfAtoms(Molecule* mol, int* pnTypes);
CCHEMIFileType getTypeFile(char* fileName);
static  boolean smallDistance(Molecule* mol);
static  boolean oneFragment(Molecule* mol);
static double* getDistancesBetweenAtoms(Molecule* mol);
static void swap2Double(double* a, double *b);
static void swap2Int(int* a, int *b);
static  double getSimilatityByBonds(Molecule* mol1, Molecule* mol2, double* pmaxDifference);
static boolean setRandomFragments(Molecule* molecule);
/****************************************************************************************************************************/
static void sortDoubles(double* tab, int n)
{
	int i;
        for (i=0;i<n-1;i++)
	{
		int k = i;
		int j;
		for(j=i+1;j<n;j++) if(tab[j]<tab[k]) k= j;
		if(k!=i) swap2Double(&tab[i],&tab[k]);
	}
}
/****************************************************************************************************************************/
static double sumDoubles(double* tab, int n)
{
	double sum=0;
	int i;
        for (i=0;i<n;i++) sum+= tab[i];
	return sum;
}
/*************************************************************************************************************************/
// sTol = 0.02 , distMaxTol = 0.7 Ang, recommanded in Jorgensen et al JCTC, 2017
static  boolean similarBonds(Molecule* mol1, Molecule* mol2, double sTol, double distMaxTol)
{
	double maxDifference =-1;
	double s =-1;
	s = getSimilatityByBonds(mol1,mol2, &maxDifference);
	fprintf(stderr,"s=%f maxDiffDist=%f\n",s,maxDifference);
	fflush(stderr);
	return (s>=0 && maxDifference>=0 && s<sTol && maxDifference<distMaxTol);
}
/*************************************************************************************************************************/
//Mathias S. Jørgensen , Michael N. Groves, and Bjørk Hammer
//J. Chem. Theory Comput., 2017, 13 (3), pp 1486–1493
//DOI: 10.1021/acs.jctc.6b01119
static  double getSimilatityByBonds(Molecule* mol1, Molecule* mol2, double* pmaxDifference)
{
	int nAtoms = mol1->nAtoms;
	double* distances1 = NULL;
	double* distances2 = NULL;
	double maxDifference = -1;
	double s = -1;
	int nTypes=0;
	int nTypes2=0;
	int nDists=mol2->nAtoms*(mol2->nAtoms-1)/2;

	*pmaxDifference = maxDifference;
	if(mol1->nAtoms<2 || mol2->nAtoms<2 || mol1->nAtoms != mol2->nAtoms) return s;

	char** types = getTypesOfAtoms(mol1, &nTypes);
	if(nTypes<1) return s;

	char** types2 = getTypesOfAtoms(mol2, &nTypes2);
	if(nTypes2 != nTypes)
	{
		freeTypesOfAtoms(types,nTypes);
		freeTypesOfAtoms(types2,nTypes2);
		return s;
	}

	freeTypesOfAtoms(types2,nTypes2);

	nDists=nAtoms*(nAtoms-1)/2;
	distances1 = getDistancesBetweenAtoms(mol1);
	distances2 = getDistancesBetweenAtoms(mol2);
	int* nAtomsByType = malloc(nTypes*sizeof(int));
	int k;
        for (k=0; k<nTypes; k++) nAtomsByType[k] = 0;

	int i;
        for (i=0;i<mol1->nAtoms;i++)
        	for (k=0; k<nTypes; k++) 
			if(!strcmp(mol1->atoms[i].prop.symbol,types[k])) 
				nAtomsByType[k]++;
		
	double* fab1 = malloc(nDists*sizeof(int));
	double* fab2 = malloc(nDists*sizeof(int));

	s=0;
	maxDifference=-1;
	int l;
        for (k=0; k<nTypes; k++) 
        for (l=k; l<nTypes; l++) 
	{
		int n1=0;
		int it1=0;
		int i,j;
        	for (i=0;i<mol1->nAtoms;i++)
        	for (j=i+1;j<mol1->nAtoms;j++)
		{
			if( (!strcmp(mol1->atoms[i].prop.symbol,types[k]) && !strcmp(mol1->atoms[j].prop.symbol,types[l]))
		     || (!strcmp(mol1->atoms[i].prop.symbol,types[l]) && !strcmp(mol1->atoms[j].prop.symbol,types[k]))
			) fab1[it1++] = distances1[n1];
			n1++;
		}
		int n2=0;
		int it2=0;
        	for (i=0;i<mol2->nAtoms;i++)
        	for (j=i+1;j<mol2->nAtoms;j++)
		{
			if( (!strcmp(mol2->atoms[i].prop.symbol,types[k]) && !strcmp(mol2->atoms[j].prop.symbol,types[l]))
		     || (!strcmp(mol2->atoms[i].prop.symbol,types[l]) && !strcmp(mol2->atoms[j].prop.symbol,types[k]))
			) fab2[it2++] = distances2[n2];
			n2++;
		}
		//fprintf(stderr,"n1=%d n2=%d it1=%d it2=%d\n",n1,n2,it1,it2); fflush(stderr);
		if(n1==n2 && it1==it2)
		{
			sortDoubles(fab1,it1);
			sortDoubles(fab2,it2);
			double sfab1=sumDoubles(fab1,it1);
			double sfab2=sumDoubles(fab2,it2);
			//fprintf(stderr,"sfab1=%f sfab2=%f\n",sfab1,sfab2); fflush(stderr);
			s+= fabs(sfab1-sfab2)/(sfab1+sfab2)*(nAtomsByType[k]+nAtomsByType[l])/(double)(nAtoms);
			int i;
        		for (i=0;i<it1;i++)
			{
				double diff=fabs(fab2[i]-fab1[i]);
				if(maxDifference<0 || maxDifference < diff) maxDifference = diff;
			}
		}
	}

	*pmaxDifference = maxDifference;

	free(nAtomsByType);
	free(fab1);
	free(fab2);
	free(distances1);
	free(distances2);
	freeTypesOfAtoms(types,nTypes);
	return s;
}
/****************************************************************************************************************************/
static double* getDistancesBetweenAtoms(Molecule* mol)
{
	double* distances = NULL;
	int n;
	if(mol->nAtoms<1) return distances;
	n = mol->nAtoms*(mol->nAtoms-1)/2;
	distances = malloc(n*sizeof(double));
	n = 0;
	int i,j;
	for (i = 0; i < mol->nAtoms-1; i++ )
	for (j = i+1; j < mol->nAtoms; j++ )
	{
		double d;
		distances[n] = 0;
		int k;
		for (k = 0; k < mol->nAtoms; k++ )
		{
			d = mol->atoms[i].coordinates[k]-mol->atoms[j].coordinates[k];
			distances[n] += d*d;
		}
		distances[n] = sqrt(distances[n]);
		n++;
	}
	return distances;
}
/****************************************************************************************************************************/
static  boolean oneFragment(Molecule* mol)
{
	// list of atoms connected to atom number 1, directly o indeirectly via a string
	int *list = malloc(mol->nAtoms*sizeof(int));
	int *inGroup = malloc(mol->nAtoms*sizeof(int));
	int nA1=1;
	int i;
	for (i = 0; i<mol->nAtoms; i++)  list[i]=-1;
	for (i = 0; i<mol->nAtoms; i++)  inGroup[i]=0;
	Atom* atoms = mol->atoms;
	int n=0;
	list[0]=0;
	inGroup[0]=1;
	int it=0;
	do{
		it++;
		n=0;
		for (i = 0; i<mol->nAtoms; i++) 
		{
			if(inGroup[i]==1) continue;
			int i1;
			for (i1 = 0; i1<nA1 && nA1<mol->nAtoms; i1++) 
			{
				int j = list[i1];
				double distance = 0;
				int k;
				for (k=0;k<3;k++) 
				{ 
					double dij = atoms[i].coordinates[k]-atoms[j].coordinates[k];
					distance +=dij*dij;
				}
				distance = sqrt(distance);
				double bd = atoms[i].prop.covalentRadii + atoms[j].prop.covalentRadii;
				bd *= BOHRTOANG;
				bd *= 1.2;
				if(distance < bd) { list[nA1+n]=i;n++;inGroup[i]=1;break;}
			}
		}
		nA1 += n;
	}while(n!=0 && it<mol->nAtoms);
	free(list);
	free(inGroup);
	return (mol->nAtoms == nA1);
}
/****************************************************************************************************************************/
static  boolean smallDistance(Molecule* mol)
{
	Atom* atoms = mol->atoms;
	boolean small = FALSE;
	int i,j;
	for (i = 0; i<mol->nAtoms-1 && !small; i++) 
	for (j = i+1; j<mol->nAtoms && !small; j++) 
	{
		double distance = 0;
		int k;
		for (k=0;k<3;k++) 
		{ 
			double dij = atoms[i].coordinates[k]-atoms[j].coordinates[k];
			distance +=dij*dij;
		}
		double bd = atoms[i].prop.covalentRadii + atoms[j].prop.covalentRadii;
		bd *= BOHRTOANG;
		bd *=0.6;
		distance = sqrt(distance);
		if(distance < bd) small=TRUE;
		//if(small) { fprintf(stderr,"d=%f bd =%f\n",distance, bd); fflush(stderr);}
	}
	return small;
}
/****************************************************************************************************************************/
static  boolean compare2Molecules(Molecule* mol1, Molecule* mol2, boolean warning)
{
	boolean ok = TRUE;
	int j;
	if(mol1->nAtoms != mol2->nAtoms)
	{
		ok = FALSE;
		fprintf(stderr,"???????????????  Error Error Error ???????????????????????????????????????????????\n");
		fprintf(stderr," The number of atoms is not the same for geometries in the 2 files\n");
		fprintf(stderr,"       Program QFFPot stopped\n");
		fprintf(stderr,"??????????????????????????????????????????????????????????????????????????????????\n");
		exit(1);
	}
	for(j=0;j<mol1->nAtoms;j++)
	{
		if(strcmp(mol1->atoms[j].prop.symbol,mol2->atoms[j].prop.symbol))
		{
			ok = FALSE;
			break;
		}
	}
	if(!ok) 
	{
		fprintf(stderr,"???????????????  Error Error Error ???????????????????????????????????????????????\n");
		fprintf(stderr," Error the order of atoms is not the same for geometries in the 2 files\n");
		fprintf(stderr,"       Program QFFPot stopped\n");
		fprintf(stderr,"??????????????????????????????????????????????????????????????????????????????????\n");
		for(j=0;j<mol1->nAtoms;j++)
		{
			if(strcmp(mol1->atoms[j].prop.symbol,mol2->atoms[j].prop.symbol))
				fprintf(stderr," %s != %s\n", mol1->atoms[j].prop.symbol,mol2->atoms[j].prop.symbol);
		}
		exit(1);
	}
	if(warning)
	{
		for(j=0;j<mol1->nAtoms;j++)
			if(warning && fabs(mol1->atoms[j].mass-mol2->atoms[j].mass)>1e-6)
				fprintf(stderr," %0.12lf != %0.12lf\n", mol1->atoms[j].mass,mol2->atoms[j].mass);
	}
	return ok;

}
/****************************************************************************************************************************/
static boolean readGeom0(Molecule* mol, char* namefile)
{
	Molecule* molNew = readGeom(namefile);
	*mol = *molNew;
	return TRUE;
}
/*******************************************************************************************************************/
static boolean readMolecule0(Molecule* molecule, char* fileName)
{
	return readGeom0(molecule, fileName);
}
/*******************************************************************************************************************/
/****************************************************************************************************************************/
// q2mat Generate a left rotation matrix from a normalized quaternion
static void q2mat (double q[], double u[][3])
{
        u[0][0] = q[0]*q[0] + q[1]*q[1] - q[2]*q[2] - q[3]*q[3];
        u[1][0] = 2.0 * (q[1] * q[2] - q[0] * q[3]);
        u[2][0] = 2.0 * (q[1] * q[3] + q[0] * q[2]);

        u[0][1] = 2.0 * (q[2] * q[1] + q[0] * q[3]);
        u[1][1] = q[0]*q[0] - q[1]*q[1] + q[2]*q[2] - q[3]*q[3];
        u[2][1] = 2.0 * (q[2] * q[3] - q[0] * q[1]);

        u[0][2] = 2.0 *(q[3] * q[1] - q[0] * q[2]);
        u[1][2] = 2.0 * (q[3] * q[2] + q[0] * q[1]);
        u[2][2] = q[0]*q[0] - q[1]*q[1] - q[2]*q[2] + q[3]*q[3];
}
/**********************************************************************/
static void rotMolecule (Molecule* mol, double u[3][3])
{
	double yx, yy, yz;
	int i;
	/*
	fprintf(stderr,"RotMatrix\n");
	for (i = 0; i < 3; i++)
	{
		int j;
		for (j = 0; j < 3; j++) fprintf(stderr,"%f ",u[i][j]);
		fprintf(stderr,"\n");
	}
	fprintf(stderr,"Mol before rot\n");  
	for(i=0;i<mol->nAtoms;i++)
	{
		fprintf(stderr,"%s ",mol->atoms[i].symbol);
		int j;
		for(j=0;j<3;j++) fprintf(stderr,"%f ",mol->atoms[i].coordinates[j]); 
		fprintf(stderr,"\n");
		
	}
	*/

	for (i = 0; i < mol->nAtoms; i++)
	{
   		yx = u[0][0] * mol->atoms[i].coordinates[0] + u[0][1] * mol->atoms[i].coordinates[1] + u[0][2] * mol->atoms[i].coordinates[2];
   		yy = u[1][0] * mol->atoms[i].coordinates[0] + u[1][1] * mol->atoms[i].coordinates[1] + u[1][2] * mol->atoms[i].coordinates[2];
   		yz = u[2][0] * mol->atoms[i].coordinates[0] + u[2][1] * mol->atoms[i].coordinates[1] + u[2][2] * mol->atoms[i].coordinates[2];

   		mol->atoms[i].coordinates[0] = yx;
   		mol->atoms[i].coordinates[1] = yy;
   		mol->atoms[i].coordinates[2] = yz;
	}
	/*
	fprintf(stderr,"Mol after rot\n");  
	for(i=0;i<mol->nAtoms;i++)
	{
		fprintf(stderr,"%s ",mol->atoms[i].symbol);
		int j;
		for(j=0;j<3;j++) fprintf(stderr,"%f ",mol->atoms[i].coordinates[j]); 
		fprintf(stderr,"\n");
		
	}
	*/
}

/**********************************************************************/
void qTransRotFit (Molecule* molFit, Molecule* molRef, double u[][3])
{
	double xxyx, xxyy, xxyz;
	double xyyx, xyyy, xyyz;
	double xzyx, xzyy, xzyz;
	double c[10];
	double d[4];
	double q[4];
	int i;
	int imax;
	double** v = malloc(4*sizeof(double*));
	for(i=0;i<4;i++) v[i] = malloc(4*sizeof(double));

/* generate the lower triangle of the quadratic form matrix */

	xxyx = 0.0;
	xxyy = 0.0;
	xxyz = 0.0;
	xyyx = 0.0;
	xyyy = 0.0;
	xyyz = 0.0;
	xzyx = 0.0;
	xzyy = 0.0;
	xzyz = 0.0;
 
	for (i = 0; i < molRef->nAtoms; i++)
	{
		//double w = molFit->atoms[i].mass*molFit->atoms[i].mass;
		//double w = molFit->atoms[i].mass;
		double w = 1.0;
		xxyx = xxyx + molFit->atoms[i].coordinates[0]*molRef->atoms[i].coordinates[0]*w;
		xxyy = xxyy + molFit->atoms[i].coordinates[0]*molRef->atoms[i].coordinates[1]*w;
		xxyz = xxyz + molFit->atoms[i].coordinates[0]*molRef->atoms[i].coordinates[2]*w;

		xyyx = xyyx + molFit->atoms[i].coordinates[1]*molRef->atoms[i].coordinates[0]*w;
		xyyy = xyyy + molFit->atoms[i].coordinates[1]*molRef->atoms[i].coordinates[1]*w;
		xyyz = xyyz + molFit->atoms[i].coordinates[1]*molRef->atoms[i].coordinates[2]*w;

		xzyx = xzyx + molFit->atoms[i].coordinates[2]*molRef->atoms[i].coordinates[0]*w;
		xzyy = xzyy + molFit->atoms[i].coordinates[2]*molRef->atoms[i].coordinates[1]*w;
		xzyz = xzyz + molFit->atoms[i].coordinates[2]*molRef->atoms[i].coordinates[2]*w;
	}
	for(i = 0; i <10; i++) c[i] = 0.0;

	c[0] = xxyx + xyyy + xzyz;

	c[1] = xzyy - xyyz;
	c[2] = xxyx - xyyy - xzyz;

	c[3] = xxyz - xzyx;
	c[4] = xxyy + xyyx;
	c[5] = xyyy - xzyz - xxyx;

	c[6] = xyyx - xxyy;
	c[7] = xzyx + xxyz;
	c[8] = xyyz + xzyy;
	c[9] = xzyz - xxyx - xyyy;

/* diagonalize c */
	eigenQL(4, c, d, v);
	imax = 0;
	for(i = 0; i <4; i++) if(d[i]>d[imax]) imax = i;

/* extract the desired quaternion */

	q[0] = v[0][imax];
	q[1] = v[1][imax];
	q[2] = v[2][imax];
	q[3] = v[3][imax];

	/*
	fprintf(stderr,"d=%f %f %f %f dimax = %f\n", d[0], d[1], d[2],d[3],d[imax]);
	fprintf(stderr,"Quat=%f %f %f %f norm=%f\n", q[0], q[1], q[2],q[3],sqrt(q[0]*q[0]+q[1]*q[1]+q[2]*q[2]+q[3]*q[3]));
	*/

/* generate the rotation matrix */
	q2mat (q, u);

	for(i=0;i<4;i++) if(v[i]) free(v[i]);
	if(v) free(v);
}

/**********************************************************************/
static void computeMCenterOfMolecule(Molecule* mol, double C[])
{
	int i;
	int j;
	double mt = 0;
	for(j=0;j<3;j++) C[j] = 0;
	for(i=0;i<mol->nAtoms;i++)
        {
	
                //mt += mol->atoms[i].mass;
                mt += 1.0;
		for(j=0;j<3;j++) C[j] += mol->atoms[i].mass*mol->atoms[i].coordinates[j];
		
                //mt += 1.0; for(j=0;j<3;j++) C[j] += mol->atoms[i].coordinates[j];
        }
	if(mt>0) for(j=0;j<3;j++) C[j] /= mt;
	/*
	fprintf(stderr,"Center = "); for(j=0;j<3;j++) fprintf(stderr,"%f ",C[j]); fprintf(stderr,"\n");
	*/
}
/**********************************************************************/
static void shiftMolecule(Molecule* mol, double C[], double scale)
{
	int i;
	int j;
	/*
	fprintf(stderr,"Before Center = "); for(j=0;j<3;j++) fprintf(stderr,"%f ",C[j]); fprintf(stderr,"\n");
	fprintf(stderr,"Mol\n");  
	for(i=0;i<mol->nAtoms;i++)
	{
		fprintf(stderr,"%s ",mol->atoms[i].symbol);
		for(j=0;j<3;j++) fprintf(stderr,"%f ",mol->atoms[i].coordinates[j]); 
		fprintf(stderr,"\n");
		
	}
	*/
	for(i=0;i<mol->nAtoms;i++)
		for(j=0;j<3;j++) mol->atoms[i].coordinates[j] += C[j]*scale;
	/*
	fprintf(stderr,"After Center = "); for(j=0;j<3;j++) fprintf(stderr,"%f ",C[j]); fprintf(stderr,"\n");
	fprintf(stderr,"Mol\n");  
	for(i=0;i<mol->nAtoms;i++)
	{
		fprintf(stderr,"%s ",mol->atoms[i].symbol);
		for(j=0;j<3;j++) fprintf(stderr,"%f ",mol->atoms[i].coordinates[j]); 
		fprintf(stderr,"\n");
		
	}
	*/
}
/**********************************************************************/
static void centerMolecule(Molecule* mol, double C[])
{
	computeMCenterOfMolecule(mol, C);
	shiftMolecule(mol, C, -1.0);
}
/**********************************************************************/
static boolean fit2Molecule(Molecule* molFit, Molecule* molRef, double u[3][3])
{
	double CFit[3];
	double CRef[3];
	centerMolecule(molFit,CFit);
	centerMolecule(molRef,CRef);

	qTransRotFit (molFit, molRef, u);
	rotMolecule (molFit, u);
	/*
	shiftMolecule(molFit,CRef,1.0);
	shiftMolecule(molRef,CRef,1.0);
	*/
	return TRUE;
}
/**********************************************************************/
static boolean fitRMSD2Molecule(Molecule* molFit, Molecule* molRef, boolean center)
{
	const int N = 8;
        double d[8];
        double w[8][3] = {
                {1,1,1},
                {-1,1,1},
                {1,-1,1},
                {1,1,-1},
                {-1,-1,1},
                {-1,1,-1},
                {1,-1,-1},
                {-1,-1,-1}
        };
        double x, y, z;
        double dmin;
	int jmin = 0;
	int i,j;
        double CFit[3];
        double CRef[3];
	int numberOfEquivalentAxes;
	double inertialMoment[3];
	double axes[3][3];
	double C[3];

	fprintf(stderr,"=======> Center molFit\n");
        centerMolecule(molFit,CFit);
	fprintf(stderr,"=======> Center molRef\n");
        centerMolecule(molRef,CRef);

	fprintf(stderr,"=======> std molFit\n");
	buildStandardOrientation(molFit, C, &numberOfEquivalentAxes, inertialMoment, axes);
	molFit->klass->print(molFit, stderr);
	fprintf(stderr,"=======> std molRef\n");
	buildStandardOrientation(molRef, C, &numberOfEquivalentAxes, inertialMoment, axes);
	molRef->klass->print(molRef, stderr);

	fprintf(stderr,"nAtomsRef=%d\n",molRef->nAtoms);
	fprintf(stderr,"nAtomsFit=%d\n",molFit->nAtoms);
	if(molFit->nAtoms != molRef->nAtoms) return FALSE;
        for (j=0;j<N;j++) d[j] = 0;
        for (i=0;i<molRef->nAtoms;i++)
        {
                for (j=0;j<N;j++)
                {
                	x = molRef->atoms[i].coordinates[0]- w[j][0]*molFit->atoms[i].coordinates[0];
                	y = molRef->atoms[i].coordinates[1]- w[j][1]*molFit->atoms[i].coordinates[1];
                	z = molRef->atoms[i].coordinates[2]- w[j][2]*molFit->atoms[i].coordinates[2];
                	d[j] += x*x + y*y + z*z;
                }
        }
        dmin = d[0];
        for (j=1;j<N;j++) if(dmin>d[j]) {dmin=d[j]; jmin = j;}
	fprintf(stderr,"jmin=%d\n",jmin);
        for (i=0;i<molFit->nAtoms;i++)
	{
                molFit->atoms[i].coordinates[0] *=  w[jmin][0];
                molFit->atoms[i].coordinates[1] *=  w[jmin][1];
                molFit->atoms[i].coordinates[2] *=  w[jmin][2];
	}
	if(!center) 
	{
		shiftMolecule(molFit,CFit,1.0);
		shiftMolecule(molRef,CRef,1.0);
	}
	fprintf(stderr,"=======> after minrms molFit\n");
	molFit->klass->print(molFit, stderr);
	fprintf(stderr,"=======> after minrms molRef\n");
	molRef->klass->print(molRef, stderr);
	return TRUE;
}
/********************************************************************************/
/*********************************************************************************/
static boolean buildMMTypes(Molecule* mol, FILE* file)
{
        char* buildMMTypes = NULL;
        if(readOneString(file,"buildMMTypes",&buildMMTypes) && buildMMTypes)
        {
                mol->klass->resetMMTypes(mol, buildMMTypes);
                free(buildMMTypes);
                return TRUE;
        }
        return FALSE;
}

/*********************************************************************************/
static boolean resetMMTypes(Molecule* mol, char* type)
{
	char* tmp = NULL;
	if(!type) return FALSE;
	tmp = strdup(type);
	uppercase(tmp);
	if(!strcmp(tmp,"AMBERFROMPDB")) 
	{
		buildMMTypesFromPDB(mol, FALSE);
		free(tmp);
		return TRUE;
	}
	if(!strcmp(tmp,"AMBERCHARGEFROMPDB")) 
	{
		buildMMTypesFromPDB(mol, TRUE);
		free(tmp);
		return TRUE;
	}
	if(!strcmp(tmp,"AMBER")) 
	{
		calculAmberTypes(mol);
		free(tmp);
		return TRUE;
	}
	fprintf(stderr,"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n");
	fprintf(stderr," Error : sorry the %s is not yet implemented\n",type);
	fprintf(stderr,"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n");
	free(tmp);
	return FALSE;
}
/*********************************************************************************/
static boolean isLinear(Molecule* molecule)
{
	double cm[3] = {0,0,0};
	int i;
	int j;
	int k;
	double mass = 1.0;
	double totMass = 0.0;
	double cdel[3];
	double tensor[3][3];
	double invTensor[3][3];
        double xx, xy,xz,yy,yz,zz;
	Atom* atoms = molecule->atoms;
	int nAtoms = molecule->nAtoms;
	double precision=1e6;
	if(nAtoms<3) return TRUE;

	for ( i = 0; i < nAtoms; i++)
	{
		mass = atoms[i].mass;
		totMass += mass;
		for ( j = 0; j < 3; j++)
			cm[j] += mass*atoms[i].coordinates[j];
	}


	for ( j = 0; j < 3; j++) cm[j] /= totMass;

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
	for ( i = 0; i < nAtoms; i++)
	{
		mass = atoms[i].mass;
		for ( j = 0; j < 3; j++) cdel[j] = atoms[i].coordinates[j]-cm[j];
		xx +=  cdel[0]*cdel[0]*mass;
		xy +=  cdel[0]*cdel[1]*mass;
		xz +=  cdel[0]*cdel[2]*mass;
		yy +=  cdel[1]*cdel[1]*mass;
		yz +=  cdel[1]*cdel[2]*mass;
		zz +=  cdel[2]*cdel[2]*mass;
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
	if(!InverseTensor(tensor,invTensor)) return TRUE;
	/*
	fprintf(stderr,"InvTensor : ");
	for ( j = 0; j < 3; j++) fprintf(stderr,"%lf ",invTensor[j][j]);
	fprintf(stderr,"\n");
	*/
	for ( j = 0; j < 3; j++) if(fabs(invTensor[j][j])>precision) return TRUE;
	return FALSE;
}
/*********************************************************************************/
static void initWallParameters(WallParameters* parameters, double E0, double rho, int n)
{
	parameters->E0 = E0;
	parameters->rho = fabs(rho);
	parameters->n = abs(n);
}
/*********************************************************************************/
static void initVibrations(Molecule* mol, int nModes, int nProps)
{
	int i;
	int j;
	int k;
	if(!mol) return;
	mol->vibration.nModes = 0;
	mol->vibration.nProperties = 0;
	mol->vibration.modes = NULL;
	if(nModes<1) return;
	mol->vibration.nModes = nModes;
	if(nProps<1) nProps = 1;
	mol->vibration.nProperties = nProps;
	mol->vibration.modes = malloc( mol->vibration.nModes*sizeof(VibMode));
	for(i=0;i< mol->vibration.nModes;i++)
	{
		mol->vibration.modes[i].frequency = 0;
		mol->vibration.modes[i].mass = 1;
		for(j=0;j<3;j++) mol->vibration.modes[i].vectors[j] = malloc(mol->nAtoms*sizeof(double));
		for(j=0;j<3;j++) 
			for(k=0;k<mol->nAtoms;k++)  mol->vibration.modes[i].vectors[j][k] = 0.0;

		mol->vibration.modes[i].properties = malloc(mol->vibration.nProperties*sizeof(double));
		for(k=0;k<mol->vibration.nProperties;k++)  mol->vibration.modes[i].properties[k] = 0.0;
	}
}
/*********************************************************************************/
static void freeVibrations(Molecule* mol)
{
	int i;
	int j;
	if(!mol) return;

	for(i=0;i< mol->vibration.nModes;i++)
	{
		mol->vibration.modes[i].frequency = 0;
		mol->vibration.modes[i].mass = 1;
		for(j=0;j<3;j++) if(mol->vibration.modes[i].vectors[j]) free(mol->vibration.modes[i].vectors[j]);
		if(mol->vibration.modes[i].properties) free(mol->vibration.modes[i].properties);
	}
	if(mol->vibration.modes)free(mol->vibration.modes);
	mol->vibration.nProperties = 0;
	mol->vibration.nModes = 0;
	mol->vibration.modes = NULL;
}
/*********************************************************************************/
static Vibration copyVibrations(Molecule* mol)
{
	int i;
	int j;
	int k;
	Vibration vib;
	vib.nModes = 0;
	vib.nProperties = 0;
	vib.modes = NULL;

	if(!mol) return vib;
	if(mol->vibration.nModes<1) return vib;
	vib.nModes = mol->vibration.nModes;
	vib.nProperties = mol->vibration.nProperties;
	vib.modes = malloc(vib.nModes*sizeof(VibMode));
	for(i=0;i<vib.nModes;i++)
	{
		vib.modes[i].frequency = mol->vibration.modes[i].frequency;
		vib.modes[i].mass = mol->vibration.modes[i].mass;
		for(j=0;j<3;j++) vib.modes[i].vectors[j] = malloc(mol->nAtoms*sizeof(double));
		for(j=0;j<3;j++) 
			for(k=0;k<mol->nAtoms;k++)  vib.modes[i].vectors[j][k] = mol->vibration.modes[i].vectors[j][k];

		vib.modes[i].properties = malloc(mol->vibration.nProperties*sizeof(double));
		for(k=0;k<vib.nProperties;k++)  vib.modes[i].properties[k] = mol->vibration.modes[i].properties[k];
	}
	return vib;
}
/*******************************************************************************************************************/
static void initBoxes(Boxes* boxes, double lx, double ly, double lz, double alpha, double beta, double gamma)
{
        double alphaCos = cos(alpha/RADTODEG);
        double betaSin = sin(beta/RADTODEG);
        double betaCos = cos(beta/RADTODEG);
        double gammaSin = sin(gamma/RADTODEG);
        double gammaCos = cos(gamma/RADTODEG);
        double betaTerm = (alphaCos - betaCos*gammaCos) / gammaSin;
        double gammaTerm = sqrt(betaSin*betaSin - betaTerm*betaTerm);
        double volbox = (gammaSin*gammaTerm) * lx * ly * lz;
	double ar1,ar2,ar3,br1,br2,br3,cr1,cr2,cr3;
	boxes->alpha = alpha;
	boxes->beta = beta;
	boxes->gamma = gamma;

//      real space lattice vectors as rows
	ar1 = lx;
      	ar2 = 0.0;
      	ar3 = 0.0;
      	br1 = ly * gammaCos;
      	br2 = ly * gammaSin;
      	br3 = 0.0;
      	cr1 = lz * betaCos;
      	cr2 = lz * betaTerm;
      	cr3 = lz * gammaTerm;
      	boxes->lVectors[0][0] = ar1;
      	boxes->lVectors[0][1] = ar2;
      	boxes->lVectors[0][2] = ar3;
      	boxes->lVectors[1][0] = br1;
      	boxes->lVectors[1][1] = br2;
      	boxes->lVectors[1][2] = br3;
      	boxes->lVectors[2][0] = cr1;
      	boxes->lVectors[2][1] = cr2;
      	boxes->lVectors[2][2] = cr3;

//     reciprocal lattice vectors as columns

	if (volbox >1e-10)
	{
		boxes->reciprocVectors[0][0] = (br2*cr3 - cr2*br3) / volbox;
         	boxes->reciprocVectors[1][0] = (br3*cr1 - cr3*br1) / volbox;
         	boxes->reciprocVectors[2][0] = (br1*cr2 - cr1*br2) / volbox;
         	boxes->reciprocVectors[0][1] = (cr2*ar3 - ar2*cr3) / volbox;
         	boxes->reciprocVectors[1][1] = (cr3*ar1 - ar3*cr1) / volbox;
         	boxes->reciprocVectors[2][1] = (cr1*ar2 - ar1*cr2) / volbox;
         	boxes->reciprocVectors[0][2] = (ar2*br3 - br2*ar3) / volbox;
         	boxes->reciprocVectors[1][2] = (ar3*br1 - br3*ar1) / volbox;
         	boxes->reciprocVectors[2][2] = (ar1*br2 - br1*ar2) / volbox;
	}
}
/**********************************************************************/
static double getGradientNorm(Molecule* molecule)
{
        double gradientNorm = 0;
	int i,j;
	for (i = 0; i < molecule->nAtoms; i++)
		for ( j = 0; j < 3; j++)
                        gradientNorm += 
				molecule->atoms[i].gradient[j] * 
				molecule->atoms[i].gradient[j]; 

        gradientNorm = sqrt( gradientNorm );
	return gradientNorm;
}
/**********************************************************************/
static double getKineticEnergyMass(Molecule* molecule, double* masses)
{
	double ekin = 0;
	int i;
	int j;
	for ( i = 0; i < molecule->nAtoms; i++)
		for ( j = 0; j < 3; j++)
			ekin += molecule->atoms[i].velocity[j]*
			        molecule->atoms[i].velocity[j]*
			        masses[i];
	ekin /=2;
	return ekin;
}
/**********************************************************************/
static double getKineticEnergy(Molecule* molecule)
{
	double ekin = 0;
	int i;
	int j;
	for ( i = 0; i < molecule->nAtoms; i++)
		for ( j = 0; j < 3; j++)
			ekin += molecule->atoms[i].velocity[j]*
			        molecule->atoms[i].velocity[j]*
			        molecule->atoms[i].mass;
	ekin /=2;
	return ekin;
}
/*****************************************************************************/
/* in Debye */
static void computeDipole(Molecule* mol)
{
	int i,k;
	for(k=0;k<3;k++) mol->dipole[k] = 0;
	for(i=0;i<mol->nAtoms;i++)
	for(k=0;k<3;k++)
		mol->dipole[k] += mol->atoms[i].charge*mol->atoms[i].coordinates[k];
	for(k=0;k<3;k++) mol->dipole[k] *= ANGTOBOHR*AUTODEB;
}
/**********************************************************************/
static void copyChargeInCharge0(Molecule* molecule)
{
	int i;
	for ( i = 0; i < molecule->nAtoms; i++) 
		molecule->atoms[i].charge0 = molecule->atoms[i].charge;
}
/**********************************************************************/
Molecule* newMolecule()
{
	int i;
	Molecule*  molecule = malloc(sizeof(Molecule));

	molecule->klass = malloc(sizeof(MoleculeClass));
	molecule->klass->print = printMolecule;
	molecule->klass->scaleVelocities = scaleVelocities;
	molecule->klass->setMaxwellVelocities = setMaxwellVelocities;
	molecule->klass->setMaxwellVelocitiesIfNull = setMaxwellVelocitiesIfNull;
	molecule->klass->getGradientNorm = getGradientNorm;
	molecule->klass->getKelvin = getKelvin;
	molecule->klass->getKineticEnergy = getKineticEnergy;
	molecule->klass->computeDipole = computeDipole;
	molecule->klass->getKineticEnergyMass = getKineticEnergyMass;
	molecule->klass->removeTranslation = removeTranslation ;
	molecule->klass->removeRotation = removeRotation ;
	molecule->klass->removeTranslationAndRotation = removeTranslationAndRotation ;

	molecule->klass->removeTranslationCluster = removeTranslationCluster ;
	molecule->klass->removeRotationCluster = removeRotationCluster ;
	molecule->klass->removeTranslationAndRotationCluster = removeTranslationAndRotationCluster ;

	molecule->klass->removeTranslationForceCluster = removeTranslationForceCluster ;
	molecule->klass->removeRotationForceCluster = removeRotationForceCluster ;
	molecule->klass->removeTranslationAndRotationForceCluster = removeTranslationAndRotationForceCluster ;

	molecule->klass->removeTranslationForce = removeTranslationForce ;
	molecule->klass->removeRotationForce = removeRotationForce ;
	molecule->klass->removeTranslationAndRotationForce = removeTranslationAndRotationForce ;

	molecule->klass->removeTranslationMoments = removeTranslationMoments ;
	molecule->klass->removeRotationMoments = removeRotationMoments ;
	molecule->klass->removeTranslationAndRotationMoments = removeTranslationAndRotationMoments ;

	molecule->klass->removeTranslationAcceleration = removeTranslationAcceleration ;
	molecule->klass->removeRotationAcceleration = removeRotationAcceleration ;
	molecule->klass->removeTranslationAndRotationAcceleration = removeTranslationAndRotationAcceleration ;
	molecule->klass->resetConstraints = resetConstraints;
	molecule->klass->setRattleConstraintsParameters = setRattleConstraintsParameters;
	molecule->klass->setRandomPositions = setRandomPositions;
	molecule->klass->setRandomPositionsChain = setRandomPositionsChain;
	molecule->klass->setRandomFragments = setRandomFragments;
	molecule->klass->addMolecule = addMolecule;
	molecule->klass->addGeometry = addGeometry;
	molecule->klass->addVelocities = addVelocities;
	molecule->klass->addFirstDerivativeToFile=addFirstDerivativeToFile;
	molecule->klass->readGeometry = readGeometry;
	molecule->klass->compare = compare2Molecules;
	molecule->klass->saveGeometry = saveGeometry;
	molecule->klass->save = saveMolecule;
	molecule->klass->saveMol2 = saveMol2;
	molecule->klass->saveHIN = saveHIN;
	molecule->klass->saveFrequencies = saveFrequencies;
	molecule->klass->computeIR = computeIR;
	molecule->klass->removeFrequencies = removeFrequencies;
	molecule->klass->removeNotSelectedFrequencies = removeNotSelectedFrequencies;
	molecule->klass->saveGeometryAndVelocities = saveGeometryAndVelocities;
	molecule->klass->isLinear = isLinear;
	molecule->klass->copyChargeInCharge0 = copyChargeInCharge0;
	molecule->klass->setBondHardness = setBondHardness;
	molecule->klass->setElectronegativity = setElectronegativity;
	molecule->klass->setWidth = setWidth;
	molecule->klass->setHardness = setHardness;
	molecule->klass->setCharge0 = setCharge0;
	molecule->klass->setChargesEEM = setChargesEEM;
	molecule->klass->setChargesACKS2 = setChargesACKS2;
	molecule->klass->getEnergyEEM = getEnergyEEM;
	molecule->klass->getEnergyACKS2 = getEnergyACKS2;
	molecule->klass->addGeometryToGabedit = addGeometryToGabedit;
	molecule->klass->readGradientFromGabeditFile = readGradientFromGabeditFile;
	molecule->klass->computeFrequenciesFromFiles = computeFrequenciesFromFiles;
	molecule->klass->computeFrequenciesOneStepFromFiles = computeFrequenciesOneStepFromFiles;
	molecule->klass->computeFrequenciesFromGradFiles = computeFrequenciesFromGradFiles;
	molecule->klass->computeFrequenciesOneStepFromGradFiles = computeFrequenciesOneStepFromGradFiles;
	molecule->klass->generateCChemIFilesForFrequencies = generateCChemIFilesForFrequencies;
	molecule->klass->generateCChemIFilesOneStepForFrequencies = generateCChemIFilesOneStepForFrequencies;
	molecule->klass->generateCChemIGradFilesForFrequencies = generateCChemIGradFilesForFrequencies;
	molecule->klass->generateCChemIGradFilesOneStepForFrequencies = generateCChemIGradFilesOneStepForFrequencies;


	molecule->klass->free = freeMolecule;
	molecule->klass->copy = copyMolecule;
	molecule->klass->computeHarmonicVelocitiesCoordinates = computeHarmonicVelocitiesCoordinates;
	molecule->klass->generateQFFCChemIFiles = generateQFFCChemIFiles;
	molecule->klass->getGeomsQFF = getGeomsQFF;

	molecule->klass->read = readMolecule0;
	molecule->klass->readGeomFromMopacOutputFile = readGeomFromMopacOutputFile;
	molecule->klass->readGeomFromGamessOutputFile = readGeomFromGamessOutputFile;
	molecule->klass->readGeomFromOrcaOutputFile = readGeomFromOrcaOutputFile;
	molecule->klass->readGeomFromOpenBabelOutputFile = readGeomFromOpenBabelOutputFile;
	molecule->klass->readGeomFromGaussianOutputFile = readGeomFromGaussianOutputFile;
	molecule->klass->readGeomFromMopacAuxFile = readGeomFromMopacAuxFile;
	molecule->klass->removeTransRotModes = removeTransRotModes;
	molecule->klass->setConnections = setConnections;
	molecule->klass->resetMMTypes = resetMMTypes;
	molecule->klass->buildMMTypes = buildMMTypes;
        molecule->klass->fit = fit2Molecule;
        molecule->klass->fitRMSD = fitRMSD2Molecule;
        molecule->klass->buildStandardOrientation = buildStandardOrientation;
        molecule->klass->center = centerMolecule;
        molecule->klass->saveFirstDerivatives = saveFirstDerivatives;
        molecule->klass->getDeltaTable = getDeltaTable;
        molecule->klass->getQFFOneGeom = getQFFOneGeom;
        molecule->klass->getQFFOneGeom3 = getQFFOneGeom3;
        molecule->klass->getQFFOneGeom4 = getQFFOneGeom4;

        molecule->klass-> computePseudoInertia = computePseudoInertia;
        molecule->klass-> similarInertia = similarInertia;
        molecule->klass-> similarBonds = similarBonds;
        molecule->klass-> makeSphereCutSpliceCrossover = makeSphereCutSpliceCrossover;
        molecule->klass-> makePlaneCutSpliceCrossover = makePlaneCutSpliceCrossover;
        molecule->klass-> makeLocalSphericalMutation = makeLocalSphericalMutation;
        molecule->klass-> makeCenterOfMassSphericalMutation = makeCenterOfMassSphericalMutation;
        molecule->klass-> setGeometryToAxes = setGeometryToAxes;
        molecule->klass-> setRandomOrientation = setRandomOrientation;
        molecule->klass-> smallDistance = smallDistance;
        molecule->klass-> oneFragment = oneFragment;
        molecule->klass-> getDistancesBetweenAtoms = getDistancesBetweenAtoms;
        molecule->klass-> getSimilatityByBonds = getSimilatityByBonds;

	molecule->nAtoms = 0;
	molecule->atoms = NULL;
	molecule->potentialEnergy = 0;
	for(i=0;i<3;i++) molecule->dipole[i] = 0;
	molecule->numberOf2Connections = 0;
	for(i=0;i<2;i++)
		molecule->connected2[i] = NULL;
	molecule->numberOf3Connections = 0;
	for(i=0;i<3;i++)
		molecule->connected3[i] = NULL;
	molecule->numberOf4Connections = 0;
	for(i=0;i<4;i++)
		molecule->connected4[i] = NULL;

	molecule->numberOfNonBonded = 0;
	for(i=0;i<2;i++)
		molecule->nonBonded[i] = NULL;

	molecule->bondHardness = NULL;

	molecule->spinMultiplicity = 1;
	molecule->totalCharge = 0;
	molecule->constraints = NOCONSTRAINTS;
	molecule->nFree = 3* molecule->nAtoms;
	molecule->numberOfRattleConstraintsTerms = 0;

	for(i=0;i<RATTLEDIM;i++)
		molecule->rattleConstraintsTerms[i] = NULL;

	initBoxes(&molecule->boxes, -1,-1,-1,0,0,0);
	initWallParameters(&molecule->wall, 0.0,1.0,2);
	initVibrations(molecule, 0, 0);

	return molecule;

}
/**********************************************************************/
static void freeMolecule(Molecule* molecule)
{

	int i;
        if(molecule->nAtoms<=0) return;
	freeVibrations(molecule);

	if(molecule->atoms != NULL)
	{
		for(i=0;i<molecule->nAtoms;i++)
		{
			if(molecule->atoms[i].prop.symbol != NULL)
				free(molecule->atoms[i].prop.symbol);
			if(molecule->atoms[i].mmType !=NULL )
				free(molecule->atoms[i].mmType);
			if(molecule->atoms[i].pdbType !=NULL )
				free(molecule->atoms[i].pdbType);
			if(molecule->atoms[i].typeConnections !=NULL )
				free(molecule->atoms[i].typeConnections);
		}

		free(molecule->atoms);
		molecule->atoms = NULL;
	}
	molecule->nAtoms = 0;
	molecule->potentialEnergy = 0;
	for(i=0;i<3;i++) molecule->dipole[i] = 0;
	molecule->numberOf2Connections = 0;
	for(i=0;i<2;i++)
	{
		if(molecule->connected2[i] != NULL)
			free(molecule->connected2[i]);
		molecule->connected2[i] = NULL;
	}
	molecule->numberOf3Connections = 0;
	for(i=0;i<3;i++)
	{
		if(molecule->connected3[i] != NULL)
			free(molecule->connected3[i]);
		molecule->connected3[i] = NULL;
	}
	molecule->numberOf4Connections = 0;
	for(i=0;i<4;i++)
	{
		if(molecule->connected4[i] != NULL)
			free(molecule->connected4[i]);
		molecule->connected4[i] = NULL;
	}
	if(molecule->klass) free(molecule->klass);
	for(i=0;i<RATTLEDIM;i++) if(molecule->rattleConstraintsTerms[i]) free(molecule->rattleConstraintsTerms[i]); 
}
/*****************************************************************************/
void createBondedMatrix(Molecule* molecule)
{
	int nAtoms = molecule->nAtoms;
	int i;
	int j;

	if(nAtoms<1)
		return;

	bondedMatrix = malloc(nAtoms*sizeof(boolean*));
	for(i=0;i<nAtoms;i++)
		bondedMatrix[i] = malloc(nAtoms*sizeof(boolean));

	for(i=0;i<nAtoms;i++)
	{
		for(j=0;j<nAtoms;j++)
			bondedMatrix[i][j] = FALSE;

		bondedMatrix[i][i] = TRUE;
	}

}
/*****************************************************************************/
void freeBondedMatrix(Molecule* molecule)
{
	int nAtoms = molecule->nAtoms;
	int i;

	if(bondedMatrix == NULL)
	       return;
	for(i=0;i<nAtoms;i++)
		if(bondedMatrix[i] != NULL)
		       	free(bondedMatrix[i]);

	free(bondedMatrix);
	bondedMatrix = NULL;

}
/*****************************************************************************/
void updatebondedMatrix(int a1, int a2)
{
	bondedMatrix[a1][a2] = TRUE;
	bondedMatrix[a2][a1] = TRUE;

}
/*****************************************************************************/
boolean isConnected2(Molecule* molecule,int i,int j)
{
	double distance;
	double dij;
	int k;
	Atom a1 = molecule->atoms[i];
	Atom a2 = molecule->atoms[j];

	if(molecule->atoms[i].typeConnections)
	{
		 	int nj = molecule->atoms[j].N-1;
			if(molecule->atoms[i].typeConnections[nj]>0) return TRUE;
			else return FALSE;
	}
	distance = 0;
	for (k=0;k<3;k++)
	{
		dij = a1.coordinates[k]-a2.coordinates[k];
		distance +=dij*dij;
	}
	distance = sqrt(distance)/BOHRTOANG;

	if(distance<(a1.prop.covalentRadii+a2.prop.covalentRadii)) return TRUE;
  	else return FALSE;
}
/*****************************************************************************/
void set2Connections(Molecule* molecule)
{
	int i;
	int j;
	int k=0;

	k = molecule->nAtoms;
	k = k*(k-1)/2;
	for(i=0;i<2;i++)
		molecule->connected2[i] = malloc(k*sizeof(int));

	k=0;
	for(i=0;i<molecule->nAtoms-1;i++)
		for(j=i+1;j<molecule->nAtoms;j++)
	{
		if(isConnected2(molecule,i,j))
		{
			molecule->connected2[0][k]= i;
			molecule->connected2[1][k]= j;

			updatebondedMatrix(i,j);

			k++;

		}
	}
	molecule->numberOf2Connections = k;
	if(k==0)
		for(i=0;i<2;i++)
		{
			free(molecule->connected2[i]);
			molecule->connected2[i] = NULL;
		}
	else
		for(i=0;i<2;i++)
			molecule->connected2[i] = realloc(molecule->connected2[i],k*sizeof(int));
	/* printing for test*/
	/*
	printf("%d 2 connections : \n",molecule->numberOf2Connections);
	for(k=0;k<molecule->numberOf2Connections;k++)
	{

		i =  molecule->connected2[0][k];
		j =  molecule->connected2[1][k];
		printf("%d-%d ",i,j);
	}
	printf("\n");
	*/


}
/*****************************************************************************/
static void permut(int* a,int *b)
{
	int c = *a;
	*a = *b;
	*b = c;
}
/*****************************************************************************/
boolean  isConnected3(Molecule* molecule,int n,int i,int j, int k)
{
	int c;
	int a1,a2,a3;
	for(c=0;c<n;c++)
	{
		a1 =  molecule->connected3[0][c];
		a2 =  molecule->connected3[1][c];
		a3 =  molecule->connected3[2][c];
		if(a1==i && a2 == j && a3 == k)
			return TRUE;
	}
	return FALSE;

}
/*****************************************************************************/
boolean  connect3(Molecule* molecule,int n,int i,int j, int k)
{
	if(i>k)permut(&i,&k);
	if(!isConnected3(molecule,n,i,j,k))
	{
		molecule->connected3[0][n]= i;
		molecule->connected3[1][n]= j;
		molecule->connected3[2][n]= k;

		updatebondedMatrix(i,j);
		updatebondedMatrix(i,k);
		updatebondedMatrix(j,k);

		return TRUE;
	}
	return FALSE;

}
/*****************************************************************************/
void set3Connections(Molecule* molecule)
{
	int i;
	int j;
	int k=0;
	int l=0;
	int n=0;

	k = molecule->numberOf2Connections*molecule->nAtoms;
	for(i=0;i<3;i++)
		molecule->connected3[i] = malloc(k*sizeof(int));

	n=0;
	for(k=0;k<molecule->numberOf2Connections;k++)
	{
		i = molecule->connected2[0][k];
		j = molecule->connected2[1][k];
		for(l=0;l<molecule->nAtoms;l++)
		{
			if(l!=i && l!=j)
			{
				if( isConnected2(molecule,i,l))
					if( connect3(molecule,n,l,i,j))
						n++;

				if( isConnected2(molecule,j,l))
					if( connect3(molecule,n,i,j,l))
						n++;
			}
		}

	}
	molecule->numberOf3Connections = n;
	if(n==0)
		for(i=0;i<3;i++)
		{
			free(molecule->connected3[i]);
			molecule->connected3[i] = NULL;
		}
	else
		for(i=0;i<3;i++)
			molecule->connected3[i] = realloc(molecule->connected3[i],n*sizeof(int));
	/* printing for test*/
	/*
	printf("%d 3 connections : \n",molecule->numberOf3Connections);
	for(k=0;k<molecule->numberOf3Connections;k++)
	{

		i =  molecule->connected3[0][k];
		j =  molecule->connected3[1][k];
		l =  molecule->connected3[2][k];
		printf("%d-%d-%d ",i,j,l);
	}
	printf("\n");
	*/


}
/*****************************************************************************/
boolean  isConnected4(Molecule* molecule,int n,int i,int j, int k,int l)
{
	int c;
	int a1,a2,a3,a4;
	for(c=0;c<n;c++)
	{
		a1 =  molecule->connected4[0][c];
		a2 =  molecule->connected4[1][c];
		a3 =  molecule->connected4[2][c];
		a4 =  molecule->connected4[3][c];

		if(a1==i && a2 == j && a3 == k && a4 == l)
			return TRUE;
	}
	return FALSE;

}
/*****************************************************************************/
boolean  connect4(Molecule* molecule,int n,int i,int j, int k,int l)
{
	if(i>l)
	{
		permut(&i,&l);
		permut(&j,&k);
	}
	if(!isConnected4(molecule,n,i,j,k,l))
	{
		molecule->connected4[0][n]= i;
		molecule->connected4[1][n]= j;
		molecule->connected4[2][n]= k;
		molecule->connected4[3][n]= l;

		updatebondedMatrix(i,j);
		updatebondedMatrix(i,k);
		updatebondedMatrix(i,l);
		updatebondedMatrix(j,k);
		updatebondedMatrix(j,l);
		updatebondedMatrix(k,l);

		return TRUE;
	}
	return FALSE;

}
/*****************************************************************************/
void set4Connections(Molecule* molecule)
{
	int i;
	int j;
	int k=0;
	int m=0;
	int l=0;
	int n=0;

	k = molecule->numberOf3Connections*molecule->nAtoms;
	for(i=0;i<4;i++)
		molecule->connected4[i] = malloc(k*sizeof(int));

	n=0;
	for(k=0;k<molecule->numberOf3Connections;k++)
	{
		i = molecule->connected3[0][k];
		j = molecule->connected3[1][k];
		m = molecule->connected3[2][k];
		for(l=0;l<molecule->nAtoms;l++)
		{
			/* a refaire voir Set3Co */
			if(l!=i && l!=j && l!= m)
			{
				if( isConnected2(molecule,i,l))
					if(connect4(molecule,n,l,i,j,m))
						n++;
				if( isConnected2(molecule,m,l))
					if(connect4(molecule,n,i,j,m,l))
						n++;
			}
		}

	}
	molecule->numberOf4Connections = n;
	if(n==0)
		for(i=0;i<4;i++)
		{
			free(molecule->connected4[i]);
			molecule->connected4[i] = NULL;
		}
	else
		for(i=0;i<4;i++)
			molecule->connected4[i] = realloc(molecule->connected4[i],n*sizeof(int));
	/* printing for test*/
	/*
	printf("%d 4 connections : \n",molecule->numberOf4Connections);
	for(k=0;k<molecule->numberOf4Connections;k++)
	{

		i =  molecule->connected4[0][k];
		j =  molecule->connected4[1][k];
		l =  molecule->connected4[2][k];
		m =  molecule->connected4[3][k];
		printf("%d-%d-%d-%d ",i,j,l,m);
	}
	printf("\n");
	*/


}
/*****************************************************************************/
void setNonBondedConnections(Molecule* molecule)
{
	int i;
	int j;
	int k;
	int numberOfNonBonded =0;
	int numberOfAtoms = molecule->nAtoms;
	int *nonBonded[2];

	k = numberOfAtoms;
	k = k*(k-1)/2;
	for(i=0;i<2;i++)
		nonBonded[i] = malloc(k*sizeof(int));

	/* list for all nonbonded connections */
	numberOfNonBonded = 0;
	for (  i = 0; i < numberOfAtoms; i++ )
		for (  j = i + 1; j < numberOfAtoms; j++ )
		{
			if ( !bondedMatrix[ i ][ j ] )
			{
				nonBonded[0][numberOfNonBonded] = i;
				nonBonded[1][numberOfNonBonded] = j;
				numberOfNonBonded++;
			}
		}
	if(numberOfNonBonded==0)
		for(i=0;i<2;i++)
		{
			free(nonBonded[i]);
			nonBonded[i] = NULL;
		}
	else
		for(i=0;i<2;i++)
		{
			nonBonded[i] = realloc(nonBonded[i],numberOfNonBonded*sizeof(int));
		}
	molecule->numberOfNonBonded = numberOfNonBonded;
	for(i=0;i<2;i++)
		molecule->nonBonded[i] = nonBonded[i];
	/* printing for test*/
	/*
	printf("%d nonBonded connections : \n",molecule->numberOfNonBonded);
	for(k=0;k<molecule->numberOfNonBonded;k++)
	{

		i =  molecule->nonBonded[0][k];
		j =  molecule->nonBonded[1][k];
		printf("%d-%d ",i,j);
	}
	printf("\n");
	*/
}
/************************************************************************/
static void setMultipleBonds(Molecule* mol)
{
	int* nBonds = NULL;
	int i;
	int j;
	int nAtoms = mol->nAtoms;
	Atom* atoms = NULL;
	if(nAtoms<1) return;
	nBonds = malloc(nAtoms*sizeof(int));
	atoms = mol->atoms;

	for(i=0;i<(int)nAtoms;i++) nBonds[i] = 0;
	for(i=0;i<(int)nAtoms;i++)
		for(j=i+1;j<(int)nAtoms;j++)
			 if(atoms[i].typeConnections && atoms[i].typeConnections[j]!=0) 
			 {
				 nBonds[i] += 1;
				 nBonds[j] += 1;
			 }
	for(i=0;i<nAtoms;i++)
	{
		SAtomsProp Prop_i = atoms[i].prop;
		if(!atoms[i].typeConnections) continue;
		for(j=i+1;j<nAtoms;j++)
		{
			SAtomsProp Prop_j;
			if(atoms[i].typeConnections[j]==0) continue;
			Prop_j = atoms[j].prop;
			if(
		 	nBonds[i] < Prop_i.maximumBondValence &&
		 	nBonds[j] < Prop_j.maximumBondValence 
			)
			{
				atoms[i].typeConnections[j] = 2;
				if(atoms[j].typeConnections) atoms[j].typeConnections[i] = 2;
				nBonds[i] += 1;
				nBonds[j] += 1;
			}
		}
	}
	for(i=0;i<nAtoms;i++)
	{
		SAtomsProp Prop_i = atoms[i].prop;
		if(!atoms[i].typeConnections) continue;
		for(j=i+1;j<nAtoms;j++)
		{
			SAtomsProp Prop_j;
			if(atoms[i].typeConnections[j]==0) continue;
			Prop_j = atoms[j].prop;
			if(
		 	nBonds[i] < Prop_i.maximumBondValence &&
		 	nBonds[j] < Prop_j.maximumBondValence 
			)
			{
				atoms[i].typeConnections[j] = 3;
				if(atoms[j].typeConnections) atoms[j].typeConnections[i] = 3;
				nBonds[i] += 1;
				nBonds[j] += 1;
			}
		}
	}
	free(nBonds);
}
/*****************************************************************************/
static boolean connected(Molecule* mol, int i,int j)
{
	double distance;
	double dif[3];
	int k;
	double d;

	for (k=0;k<3;k++) dif[k] = mol->atoms[i].coordinates[k] -  mol->atoms[j].coordinates[k];
	distance = 0;
	for (k=0;k<3;k++) distance += dif[k]*dif[k];
	distance = sqrt(distance);
  
	d  = mol->atoms[i].width;
	d += mol->atoms[j].width;
	d *= BOHRTOANG;

	//printf("%d %d dis=%f w=%f\n",i+1,j+1,distance,d);
	if(distance<d) return TRUE;
	else return FALSE;
}
/*****************************************************************************/
static void resetTypeConnections(Molecule* mol)
{
	int i,j;
	int nAtoms = mol->nAtoms;
	Atom* atoms = NULL;
	if(nAtoms<1) return;
	atoms = mol->atoms;
	for(i=0;i<nAtoms;i++)
	{
		if(atoms[i].typeConnections) free(atoms[i].typeConnections);
		atoms[i].typeConnections = malloc(nAtoms*sizeof(int));
		for(j=0;j<nAtoms;j++) atoms[i].typeConnections[j] = 0;
		for(j=0;j<nAtoms;j++) 
			if(i!=j && connected(mol,i,j)) atoms[i].typeConnections[j]=1;
	}
	setMultipleBonds(mol);
}
/*****************************************************************************/
static void setConnections(Molecule* molecule)
{
	createBondedMatrix(molecule);

	/* printf("Set Connection\n");*/
	printf(("Establishing connectivity : 2 connections...\n"));
	set2Connections(molecule);
	printf(("Establishing connectivity : 3 connections...\n"));
	set3Connections(molecule);
	printf(("Establishing connectivity : 4 connections...\n"));
	set4Connections(molecule);
	printf(("Establishing connectivity : non bonded ...\n"));
	setNonBondedConnections(molecule);

	freeBondedMatrix(molecule);
}
/*****************************************************************************/
static Molecule* copyMolecule(Molecule* m)
{

	int i;
	int j;
	int k;
	Molecule* molecule = newMolecule();

	molecule->constraints = m->constraints;
	molecule->nFree = m->nFree;
	molecule->klass = malloc(sizeof(MoleculeClass));
	*molecule->klass = *m->klass;
	molecule->potentialEnergy = m->potentialEnergy;
	for(i=0;i<3;i++) molecule->dipole[i] = m->dipole[i];
	molecule->nAtoms = m->nAtoms;
	molecule->spinMultiplicity =  m->spinMultiplicity;
	molecule->totalCharge =  m->totalCharge;
	if( molecule->nAtoms>0) molecule->atoms = malloc(molecule->nAtoms*sizeof(Atom));
	if(m && m->bondHardness)
	{
		if(molecule->nAtoms>0) molecule->bondHardness = malloc(molecule->nAtoms*(molecule->nAtoms+1)/2*sizeof(double));
		for(i=0;i<molecule->nAtoms*(molecule->nAtoms+1)/2;i++) molecule->bondHardness[i] = m->bondHardness[i];
	}
		

	for(i=0;i<molecule->nAtoms;i++)
	{
		molecule->atoms[i].prop = propAtomGet(m->atoms[i].prop.symbol);
		for(j=0;j<3;j++) molecule->atoms[i].coordinates[j] = m->atoms[i].coordinates[j];
		for(j=0;j<3;j++) molecule->atoms[i].gradient[j] = m->atoms[i].gradient[j];
		for(j=0;j<3;j++) molecule->atoms[i].velocity[j] = m->atoms[i].velocity[j];
		molecule->atoms[i].charge = m->atoms[i].charge;
		molecule->atoms[i].charge0 = m->atoms[i].charge0;
		molecule->atoms[i].electronegativity = m->atoms[i].electronegativity;
		molecule->atoms[i].hardness = m->atoms[i].hardness;
		molecule->atoms[i].width = m->atoms[i].width;
		molecule->atoms[i].mass = m->atoms[i].mass;
		molecule->atoms[i].mmType = strdup(m->atoms[i].mmType);
		molecule->atoms[i].pdbType = strdup(m->atoms[i].pdbType);
		molecule->atoms[i].residueName = strdup(m->atoms[i].residueName);
		molecule->atoms[i].residueNumber = m->atoms[i].residueNumber;
		molecule->atoms[i].layer = m->atoms[i].layer;
		molecule->atoms[i].show = m->atoms[i].show;
		molecule->atoms[i].variable = m->atoms[i].variable;
		molecule->atoms[i].N = m->atoms[i].N;
		molecule->atoms[i].rho = m->atoms[i].rho;
		molecule->atoms[i].U = m->atoms[i].U;

		molecule->atoms[i].typeConnections = NULL; 
		if(m->atoms[i].typeConnections)
		{
			int j;
			molecule->atoms[i].typeConnections = malloc(molecule->nAtoms*sizeof(int));
			for(j=0;j<molecule->nAtoms;j++)
				molecule->atoms[i].typeConnections[j] = m->atoms[i].typeConnections[j];
		}
	}
	/*
	printf("End copyCoordinate\n");
	fflush(stdout);
	*/

	molecule->numberOf2Connections = m->numberOf2Connections;
	k = molecule->numberOf2Connections;
	if(k>0)
	for(j=0;j<2;j++)
	{
		molecule->connected2[j] = malloc(k*sizeof(int));
		for(i=0;i<k;i++) molecule->connected2[j][i] = m->connected2[j][i];
	}
	/*
	printf("End copyConnections2\n");
	fflush(stdout);
	*/
	molecule->numberOf3Connections = m->numberOf3Connections;
	k = molecule->numberOf3Connections;
	if(k>0)
	for(j=0;j<3;j++)
	{
		molecule->connected3[j] = malloc(k*sizeof(int));
		for(i=0;i<k;i++) molecule->connected3[j][i] = m->connected3[j][i];
	}
	/*
	printf("End copyConnections3\n");
	fflush(stdout);
	*/
	molecule->numberOf4Connections = m->numberOf4Connections;
	k = molecule->numberOf4Connections;
	if(k>0)
	for(j=0;j<4;j++)
	{
		molecule->connected4[j] = malloc(k*sizeof(int));
		for(i=0;i<k;i++) molecule->connected4[j][i] = m->connected4[j][i];
	}
	/*
	printf("End copyConnections4\n");
	fflush(stdout);
	*/

	molecule->numberOfNonBonded = m->numberOfNonBonded;
	k = molecule->numberOfNonBonded;
	if(k>0)
	for(j=0;j<2;j++)
	{
		molecule->nonBonded[j] = malloc(k*sizeof(int));
		for(i=0;i<k;i++) molecule->nonBonded[j][i] = m->nonBonded[j][i];
	}
	/*
	printf("End copyNonBonded\n");
	fflush(stdout);
	*/

	/*
	printf("End copyGradient\n");
	fflush(stdout);
	*/
	molecule->numberOfRattleConstraintsTerms = m->numberOfRattleConstraintsTerms;
	k = molecule->numberOfRattleConstraintsTerms;
	if(k>0)
	for(i=0;i<RATTLEDIM;i++)
	{
		molecule->rattleConstraintsTerms[i] = malloc(k*sizeof(double));
		for(j=0;j<k;j++) molecule->rattleConstraintsTerms[i][j] = m->rattleConstraintsTerms[i][j];
	}
	/*
	printf("End copyMolecule\n");
	fflush(stdout);
	*/

	molecule->boxes = m->boxes;
	molecule->wall = m->wall;
	molecule->vibration = copyVibrations(m);

	return molecule;
}
/********************************************************************************/
/*
static void getChargesFromGaussianOutputFile(Molecule* mol, FILE* file)
{
  	char t[BSIZE];
  	char dump[BSIZE];
  	char d[BSIZE];
  	char* pdest;
	int i;
	int ngrad=0;


  	while(!feof(file) )
	{
    		pdest = NULL;
    		if(!fgets(t,BSIZE,file)) break;
    		pdest = strstr( t, "Total atomic charges");
		if(!pdest) // Gaussian 03 
    			pdest = strstr( t, "atomic charges");

		if(pdest)
		{
    			if(!fgets(t,BSIZE,file)) break;

			for(i=0;i<mol->nAtoms;i++)
			{
    				if(!fgets(t,BSIZE,file)) break;
				if(sscanf(t,"%s %s %s",dump,dump,d)==3)
				{
					mol->atoms[i].charge = atof(d);
				}
			}
			break;
		}
		else
		{
          		pdest = strstr( t, "GradGradGrad" );
			if(pdest)
			{
				ngrad++;
			}
			if(ngrad>2)
				break;
		}

	}
}
static void getNaturalChargesFromGaussianOutputFile(Molecule* mol, FILE* file)
{
  	char t[BSIZE];
  	char dump[BSIZE];
  	char d[BSIZE];
  	char* pdest;
	int i;
	int ngrad =0;



  	while(!feof(file) )
	{
    		pdest = NULL;
    		if(!fgets(t,BSIZE,file)) break;
    		pdest = strstr( t, "Summary of Natural Population Analysis:");
		if(!pdest) // Gaussian 03 
    			pdest = strstr( t, "Summary of Natural Population Analysis:");

		if(pdest)
		{
    			if(!fgets(t,BSIZE,file)) break;
    			if(!fgets(t,BSIZE,file)) break;


			if(!strstr(t,"Natural Population"))break;
    			if(!fgets(t,BSIZE,file)) break;
			if(!strstr(t,"Natural"))break;
    			if(!fgets(t,BSIZE,file)) break;
			if(!strstr(t,"Charge"))break;
    			if(!fgets(t,BSIZE,file)) break;
			if(!strstr(t,"-------------"))break;

			for(i=0;i<mol->nAtoms;i++)
			{
    				if(!fgets(t,BSIZE,file)) break;
				if(sscanf(t,"%s %s %s",dump,dump,d)==3)
				{
					mol->atoms[i].charge = atof(d);
				}
			}
			break;
		}
		else
		{
          		pdest = strstr( t, "GradGradGrad" );
			if(pdest)
			{
				ngrad++;
			}
			if(ngrad>2)
				break;
		}

	}
}
static void getEpsChargesFromGaussianOutputFile(Molecule* mol, FILE* file)
{
  	char t[BSIZE];
  	char dump[BSIZE];
  	char d[BSIZE];
  	char* pdest;
	int i;
	int ngrad=0;


  	while(!feof(file) )
	{
    		pdest = NULL;
    		if(!fgets(t,BSIZE,file)) break;
    		pdest = strstr( t, "Charges from ESP fit");
		if(!pdest) // Gaussian 03 
    			pdest = strstr( t, "harges from ESP");

		if(pdest)
		{
    			if(!fgets(t,BSIZE,file)) break;
    			if(!fgets(t,BSIZE,file)) break;

			for(i=0;i<mol->nAtoms;i++)
			{
    				if(!fgets(t,BSIZE,file)) break;
				if(sscanf(t,"%s %s %s",dump,dump,d)==3)
				{
					mol->atoms[i].charge = atof(d);

				}
			}
			break;
		}
		else
		{
          		pdest = strstr( t, "GradGradGrad" );
			if(pdest)
			{
				ngrad++;
			}
			if(ngrad>2)
				break;
		}

	}
}
*/
/********************************************************************************/
static void getChargesFromGamessOutputFile(Molecule* mol, FILE* file)
{
  	char t[BSIZE];
  	char dump[BSIZE];
  	char d[BSIZE];
  	char* pdest;
	int i;


  	while(!feof(file) )
	{
    		pdest = NULL;
    		if(!fgets(t,BSIZE,file)) break;
    		pdest = strstr( t, "TOTAL MULLIKEN AND LOWDIN ATOMIC POPULATIONS");

		if(pdest)
		{
    			if(!fgets(t,BSIZE,file)) break;
			for(i=0;i<mol->nAtoms;i++)
			{
    				if(!fgets(t,BSIZE,file)) break;
				if(sscanf(t,"%s %s %s %s %s %s",dump, dump ,dump, dump, dump, d)==6)
				{
					mol->atoms[i].charge = atof(d);
				}
				else break;
			}
			break;
		}
	}
}
/********************************************************************************/
void getChargesFromOrcaOutputFile(Molecule* mol, FILE* file)
{
  	char t[BSIZE];
  	char dump[BSIZE];
  	char d[BSIZE];
  	char* pdest;
	int i;


	for(i=0;i<mol->nAtoms;i++) mol->atoms[i].charge = 0.0;

  	while(!feof(file) )
	{
    		pdest = NULL;
		if(!fgets(t,BSIZE,file)) break;
		//if(strstr(t,"GEOMETRY OPTIMIZATION CYCLE")) break;
    		pdest = strstr( t, "MULLIKEN ATOMIC CHARGES");

		if(pdest)
		{
			boolean OK = FALSE;
  			while(!feof(file) )
			{
    				if(!fgets(t,BSIZE,file)) break;
				if(strstr(t,"----------------"))
				{
					OK = TRUE;
					break;
				}
			}
			if(!OK) break;

			for(i=0;i<mol->nAtoms;i++)
			{
				if(!fgets(t,BSIZE,file)) break;
				if(sscanf(t,"%s %s %s %s",dump,dump,dump,d)==4)
				{
					mol->atoms[i].charge = atof(d);
				}
			}
			break;
		}
	}
}
/********************************************************************************/
void readGeomFromMopacOutputFile(Molecule* mol, char *fileName, int numgeometry)
{
/*      Change only the coordinates of mol */
	char* t;
	boolean OK;
	char *AtomCoord[5];
	FILE *file;
	int idummy;
	int i;
	int j=0;
	int l;
	int numgeom;
	char *pdest;
	long int geomposok = 0;

	for(i=0;i<5;i++) AtomCoord[i]=malloc(BSIZE*sizeof(char));
	t=malloc(BSIZE*sizeof(char));
	 
	file = fopen(fileName, "r");
	if(file ==NULL)
	{
	 	free(t);
	 	printf(("Sorry\nI can not open %s  mopac output file\n"),fileName);
		exit(1);
	 	return;
	}
	numgeom =0;
	OK=FALSE;
	 while(!feof(file))
	 {
		if(!fgets(t,BSIZE,file))break;
		//pdest = strstr( t, "CARTESIAN COORDINATES");
		pdest = strstr( t, " ATOM   CHEMICAL          X               Y               Z");
		if(pdest) 
		{
			//printf(("%s\n"),pdest);
			if(!fgets(t,BSIZE,file)) {pdest=0;break;}
			//printf(("%s\n"),pdest);
			if(!fgets(t,BSIZE,file)) {pdest=0;break;}
			//printf(("%s\n"),pdest);
			//if(!fgets(t,BSIZE,file)) {pdest=0;break;}
			//printf(("%s\n"),pdest);
		}
		if ( pdest )
		{
			numgeom++;
			geomposok = ftell(file);
			if(numgeom == numgeometry )
			{
				OK = TRUE;
				break;
			}
			if(numgeometry<0)
			{
				OK = TRUE;
			}
		}
	 }
	 if(!OK || numgeom == 0)
	 {
		free(t);
	 	printf(("Sorry\nI can not open %s mopac output file\n"),fileName);
		exit(1);
	 	return;
	  }
	j=-1;
	fseek(file, geomposok, SEEK_SET);
	while(!feof(file) )
	{
		if(!fgets(t,BSIZE,file))break;
		if(isABackspace(t))
		{
			break;
		}
		int ii;
		if(sscanf(t,"%d",&ii)==0) {break;}
		j++;
		if(j> mol->nAtoms-1)
		{
	 		free(t);
	 		printf(("Sorry\nnumber of Atoms read > number of atoms in mol. Mopac output file =  %s\n"),fileName);
			exit(1);
	 		return;
		}
		
		for(ii=0;ii<strlen(t);ii++) if (t[ii]=='*') t[ii] = ' ';
		sscanf(t,"%d %s %s %s %s",&idummy,AtomCoord[0],AtomCoord[1],AtomCoord[2],AtomCoord[3]);
		AtomCoord[0][0]=toupper(AtomCoord[0][0]);
		l=strlen(AtomCoord[0]); 
		if(isdigit(AtomCoord[0][1]))l=1;
		if (l==2) AtomCoord[0][1]=tolower(AtomCoord[0][1]);
		if(l==1)sprintf(t,"%c",AtomCoord[0][0]);
		else sprintf(t,"%c%c",AtomCoord[0][0],AtomCoord[0][1]);
		/* test symbol to do */

		mol->atoms[j].coordinates[0]=atof(AtomCoord[1]);
		mol->atoms[j].coordinates[1]=atof(AtomCoord[2]);
		mol->atoms[j].coordinates[2]=atof(AtomCoord[3]);
	  }
	 fclose(file);
	 free(t);
	 for(i=0;i<5;i++) free(AtomCoord[i]);
}
/********************************************************************************/
void readGeomFromGamessOutputFile(Molecule* mol, char *fileName, int numgeometry)
{
	char *t;
	boolean OK;
	char *AtomCoord[5];
	FILE *file;
	int i;
	int j=0;
	int l;
	int numgeom;
	char dum[100];


	for(i=0;i<5;i++) AtomCoord[i]=malloc(BSIZE*sizeof(char));
  
	t=malloc(BSIZE*sizeof(char));
 	file = fopen(fileName, "rb");
	if(file ==NULL)
	{
		free(t);
		printf(("Sorry\nI can not open %s Firefly (gamess) output file\n"),fileName);
		exit(1);
		return;
	}
	numgeom = 0;
	do 
	{
		OK=FALSE;
		while(!feof(file)){
			fgets(t,BSIZE,file);
			if ( numgeometry==1 && strstr(t,"COORDINATES (BOHR)"))
			{
	  			fgets(t,BSIZE,file);
 				numgeom++;
				if((int)numgeom == numgeometry ) { OK = TRUE; break; }
	  		}
			if ( strstr(t,"COORDINATES OF ALL ATOMS ARE (ANGS)"))
			{
	  			fgets(t,BSIZE,file);
	  			fgets(t,BSIZE,file);
 				numgeom++;
				if((int)numgeom == numgeometry ) { OK = TRUE; break; }
				if(numgeometry<0 ) { OK = TRUE; break; }
	  		}
		}
		if(!OK && (numgeom == 0) ){
			free(t);
			printf(("Sorry\nI can not open read geometry from %s Firefly (gamess) output file\n"),fileName);
			exit(1);
			return;
		}
		if(!OK)break;

		j=-1;
		while(!feof(file) )
		{
			fgets(t,BSIZE,file);
			strDeleten(t);
			if (isABackspace(t)) break;
			if ( !strcmp(t,"\n")) break;
			if ( !strcmp(t,"\r\n")) break;
			j++;
			if(j> mol->nAtoms-1)
			{
	 			free(t);
	 			printf(("Sorry\nnumber of Atoms read > number of atoms in mol. Mopac output file =  %s\n"),fileName);
				exit(1);
	 			return;
			}

			sscanf(t,"%s %s %s %s %s",AtomCoord[0],dum, AtomCoord[1], AtomCoord[2],AtomCoord[3]);
			{
				int k;
				for(k=0;k<(int)strlen(AtomCoord[0]);k++) if(isdigit(AtomCoord[0][k])) AtomCoord[0][k] = ' ';
				deleteAllSpaces(AtomCoord[0]);
			}

			AtomCoord[0][0]=toupper(AtomCoord[0][0]);
			l=strlen(AtomCoord[0]);
			if (l==2) AtomCoord[0][1]=tolower(AtomCoord[0][1]);
			mol->atoms[j].coordinates[0]=atof(AtomCoord[1]);
			mol->atoms[j].coordinates[1]=atof(AtomCoord[2]);
			mol->atoms[j].coordinates[2]=atof(AtomCoord[3]);
			if(j==mol->nAtoms-1) break;
		}
		if(OK && numgeometry>=0) break;
	}while(!feof(file));
	fclose(file);
	free(t);
	for(i=0;i<5;i++) free(AtomCoord[i]);
}
/********************************************************************************/
void readGeomFromOrcaOutputFile(Molecule* mol, char* fileName, int numgeometry)
{
	char *t;
	char *AtomCoord[5];
	FILE *file;
	int i;
	int j=0;
	int l;
	int numgeom;
	char *pdest;
	long int geomposok = 0;

	for(i=0;i<5;i++) AtomCoord[i]=malloc(BSIZE*sizeof(char));
	 
	t=malloc(BSIZE*sizeof(char));
	file = fopen(fileName, "r");
	if(file ==NULL)
	{
	 	free(t);
	 	printf(("Sorry\nI can not open the %s  orca output file\n"),fileName);
	 	return;
	}
	numgeom =0;
	 while(!feof(file))
	 {
		if(!fgets(t,BSIZE,file))break;
		pdest = strstr( t, "CARTESIAN COORDINATES (ANGSTROEM)");
		if(pdest) 
		{
			if(!fgets(t,BSIZE,file))break;
			pdest = strstr( t, "---------------------------------");
		}
		if ( pdest )
		{
			numgeom++;
			geomposok = ftell(file);
			if(numgeom == numgeometry )
			{
				break;
			}
			if(numgeometry<0)
			{
			}
		}
	 }
	 if(numgeom == 0)
	 {
		free(t);
		t = strdup_printf(("Sorry\nI can not read geometry from %s  orca output file\n"),fileName);
		return;
	  }
	j=-1;
	fseek(file, geomposok, SEEK_SET);
	while(!feof(file) )
	{
		if(!fgets(t,BSIZE,file))break;
		pdest = strstr( t, "----------------------------------" );
		if (pdest || isABackspace(t))
		{
			break;
		}
		j++;
		if(j> mol->nAtoms-1)
		{
	 		free(t);
	 		printf(("Sorry\nnumber of Atoms read > number of atoms in mol. Mopac output file =  %s\n"),fileName);
			exit(1);
	 		return;
		}
		sscanf(t,"%s %s %s %s",AtomCoord[0],AtomCoord[1],AtomCoord[2],AtomCoord[3]);
		AtomCoord[0][0]=toupper(AtomCoord[0][0]);
		l=strlen(AtomCoord[0]); 
		if(isdigit(AtomCoord[0][1]))l=1;
		if (l==2) AtomCoord[0][1]=tolower(AtomCoord[0][1]);
		if(l==1)sprintf(t,"%c",AtomCoord[0][0]);
		else sprintf(t,"%c%c",AtomCoord[0][0],AtomCoord[0][1]);
		mol->atoms[j].coordinates[0]=atof(AtomCoord[1]);
		mol->atoms[j].coordinates[1]=atof(AtomCoord[2]);
		mol->atoms[j].coordinates[2]=atof(AtomCoord[3]);
	  }
	 fclose(file);
	 free(t);
	 for(i=0;i<5;i++) free(AtomCoord[i]);
}
/********************************************************************************/
void readGeomFromOpenBabelOutputFile(Molecule* mol, char* fileName, int numgeometry)
{
	FILE* file = NULL;
	char buffer[BSIZE];
	char* pdest = NULL;
	//char* energyTag = "TOTAL ENERGY =";
	char* energyTag = "FINAL ENERGY:";
	char* geomTag = "Geometry";
	char* gradTag = "Gradients:";
	double dum;

	printf("Read geom from %s\n",fileName);
 	file = fopen(fileName, "r");
	if(!file) return;
	 while(!feof(file))
	 {
		if(!fgets(buffer,BSIZE,file))break;
		pdest = strstr( buffer, energyTag);
		if(pdest &&sscanf(pdest+strlen(energyTag)+1,"%lf",&mol->potentialEnergy)==1)
		{
			if(strstr(pdest,"kJ")) mol->potentialEnergy /= KCALTOKJ;
			break;
		}
	 }
	 while(!feof(file))
	 {
		if(!fgets(buffer,BSIZE,file))break;
		if(strstr(buffer, geomTag))
		{
			int i;
			for(i=0;i<mol->nAtoms;i++)
			{
				if(!fgets(buffer,BSIZE,file))break;
				//printf("%s\n",buffer);
				if(sscanf(buffer,"%lf %lf %lf %lf",
					&dum,
					&mol->atoms[i].coordinates[0],
					&mol->atoms[i].coordinates[1],
					&mol->atoms[i].coordinates[2]
					)!=4) break;
			}
			break;
		}
	 }
	 while(!feof(file))
	 {
		if(!fgets(buffer,BSIZE,file))break;
		if(strstr(buffer, gradTag))
		{
			int i;
			for(i=0;i<mol->nAtoms;i++)
			{
				if(!fgets(buffer,BSIZE,file))break;
				//printf("%s\n",buffer);
				if(sscanf(buffer,"%lf %lf %lf",
					&mol->atoms[i].gradient[0],
					&mol->atoms[i].gradient[1],
					&mol->atoms[i].gradient[2]
					)!=3) break;
			}
			break;
		}
	 }
	fclose(file);
}
/***********************************************************************************************/
/*
static boolean readOrcaFile_SpinCharge(Molecule* mol, FILE*file)
{
 	char t[BSIZE];
 	char t1[20];
 	char t2[20];
 	char t3[20];
 	char t4[20];
	int dum;
	int n=0;
	mol->spinMultiplicity = 1;
	mol->totalCharge = 0;
 	while(!feof(file))
	{
		if(!fgets(t,BSIZE,file)) break;
 		if (strstr( t,"Total Charge") && strstr( t,"...") )
		{
			if(5==sscanf(t,"%s %s %s %s %d",t1,t2,t3,t4,&dum)) { mol->totalCharge = dum; n++;}
		}
 		if (strstr( t,"Multiplicity") && strstr( t,"...") )
		{
			if(4==sscanf(t,"%s %s %s %d",t1,t2,t3,&dum)) { mol->spinMultiplicity = dum; n++;}
		}
		if(n>=2) break;
	}
	return n>=2;
}
*/
/***********************************************************************************************/
static boolean readOrcaHessianFile_IR(Molecule* mol, FILE*file, char* tag)
{
	int nFreqs = 0;
 	char t[BSIZE];
	double frequency;
	double val[5] = {0,0,0,0,0};
	int nf;
	int i;
	int j;
	int itype = 0;
        double mu0 = 4*PI*1e-7;
        double eps0 = 1.0/(mu0*slight*slight);
        double   kmmolm1 = 4*PI*PI*PI*NAvogadro/3/hPlank/slight/4/PI/eps0*1e-3*100.0*8.47835267e-30*8.47835267e-30;/* 1e-3 m to km, 100 : cm-1 to m-1 */
	double f = 1.0/sqrt(kmmolm1);
 	while(!feof(file))
	{
		if(!fgets(t,BSIZE,file)) break;
 		if (strstr( t,tag) )
		{
			if(!fgets(t,BSIZE,file)) break;
			sscanf(t,"%d",&nFreqs);
			break;
		}
	}
	if(nFreqs<1) return FALSE;
	if(mol->nAtoms*3 != nFreqs)
	{
		fprintf(stderr,"Error : dimension of ir/raman vector is not equal to 3*number of Atoms\n");
		return FALSE;
	}
	
	if(strstr(tag,"raman_")) { itype = 1; f = 1.0;}
	for(i = 0;i<nFreqs;i++)
	{
		if(!fgets(t,BSIZE,file)) break;
		nf = sscanf(t,"%lf %lf %lf %lf %lf", &frequency,&val[itype],&val[2],&val[3],&val[4]);
		if(nf<5) break;
		
		mol->vibration.modes[i].properties[itype] = val[itype];
		for(j=2;j<5;j++) mol->vibration.modes[i].properties[j] = val[j]*f;
		printf("Freq calc from Hessian = %0.5f Freq from orca = %0.5f\n", frequency, mol->vibration.modes[i].frequency);
	}
	return TRUE;
}
/********************************************************************************/
static boolean readOrcaHessianFile_Hessian(Molecule* mol, FILE*file)
{
	int nFreqs = 0;
 	char t[BSIZE];
	int nblock, iblock;
	double** hessian = NULL;
	double* frequencies = NULL;
	double* effectiveMasses = NULL;
	double** modes = NULL;
	double* F = NULL;
	int jh;
	int nf;
	int i,j,k;
	double v[6];
	int ih[6];
	int c;
	rewind(file);
 	while(!feof(file))
	{
    		if(!fgets(t,BSIZE,file)) break;
 		if (strstr( t,"$hessian") )
		{
    			fgets(t,BSIZE,file);
			/* printf("t=%s\n",t);*/
			sscanf(t,"%d",&nFreqs);
			break;
		}
	}
	//printf("nFreqs = %d\n",nFreqs);
	if(nFreqs<1) return FALSE;
	if(mol->nAtoms*3 != nFreqs)
	{
		fprintf(stderr,"Error : dimension of hessian matrix is not equal to 3*number of Atoms\n");
		return FALSE;
	}
	hessian = malloc(nFreqs*sizeof(double*));
	for(j=0;j<nFreqs;j++) hessian[j] = malloc(nFreqs*sizeof(double));
	for(j=0;j<nFreqs;j++) for(i=0;i<nFreqs;i++) hessian[i][j] = 0.0;

	nblock = nFreqs/6; 
	if(nFreqs%6!=0) nblock++;
	for(iblock = 0;iblock<nblock;iblock++)
	{
		if(!fgets(t,BSIZE,file)) break;
		/* printf("t=%s\n",t);*/
		nf = sscanf(t,"%d %d %d %d %d %d", &ih[0],&ih[1],&ih[2], &ih[3],&ih[4],&ih[5]);
		if(iblock==0 && nf>0)
		{
			nblock = nFreqs/nf; 
			if(nFreqs%nf!=0) nblock++;
		}
		for(j=0;j<nFreqs && !feof(file);j++)
		{
			if(!fgets(t,BSIZE,file)) break;
			/* printf("t=%s\n",t);*/
			nf = sscanf(t,"%d %lf %lf %lf %lf %lf %lf",
					&jh,
					&v[0],&v[1],&v[2],
					&v[3],&v[4],&v[5]
					);
			nf--;
			if(jh<nFreqs && jh>-1)
			for(k=0;k<nf;k++)
			{
		
				if(ih[k]<nFreqs && ih[k]>-1) hessian[jh][ih[k]]  = v[k]; 
			}
		}
	}
        for(i=0;i<mol->nAtoms;i++) for(c=0;c<3;c++) 
	for(j=0;j<mol->nAtoms;j++) for(k=0;k<3;k++) 
		hessian[3*i+c][3*j+k] /= sqrt(mol->atoms[i].mass*mol->atoms[j].mass)*AMUTOAU;

	/* symmetrize */
        for(i=0;i<nFreqs;i++)
        for(j=0;j<i;j++)
	{
		hessian[i][j] = (hessian[i][j]+hessian[j][i])/2;
		hessian[j][i] = hessian[i][j];
	}

	/* printf("Ok end div mass\n");*/
	modes = malloc(nFreqs*sizeof(double*));
	for(j=0;j<nFreqs;j++) modes[j] = malloc(nFreqs*sizeof(double));
	frequencies = malloc(nFreqs*sizeof(double));
	effectiveMasses = malloc(nFreqs*sizeof(double));
	/* printf("Ok end read hessian\n");*/
	F = malloc(nFreqs*(nFreqs+1)/2*sizeof(double));
        /* F is an inf symmetric matrix */
        k = -1;
        for(i=0;i<nFreqs;i++)
        for(j=0;j<=i;j++)
        {
                k++;
		F[k] = hessian[i][j];
        }

	eigenQL(nFreqs, F, frequencies, modes);
	free(F);
	/* printf("Ok end eigenQL hessian\n");*/
        /* convert frequencies in cm-1 */
        for(i=0;i<nFreqs;i++)
             if( (frequencies)[i]>0) frequencies[i] = sqrt(frequencies[i])*AUTOCM1;
             else frequencies[i] = -sqrt(-frequencies[i])*AUTOCM1;

        /* for(i=0;i<nFreqs;i++) printf("freq = %f\n",frequencies[i]);*/

        /* compute the effective masses */
        for(i=0;i<nFreqs;i++)
        {
                double m = 0;
                for(j=0;j<mol->nAtoms;j++)
                {
                        double r2 = 0;
                        for(c=0;c<3;c++) r2+= modes[3*j+c][i]*modes[3*j+c][i];
			/* printf("masse = %f\n", GeomOrb[j].Prop.masse);*/
                        m+= r2/(mol->atoms[j].mass);
                }
                if(m<=0) m = 1;
                m = 1.0/m;
                for(j=0;j<mol->nAtoms;j++)
                {
                        double r =sqrt(m)/sqrt(mol->atoms[j].mass);
                        for(c=0;c<3;c++) modes[3*j+c][i] *= r;
                }

                //printf("%f %f\n",(*frequencies)[i],m);
                effectiveMasses[i] = m;
        }
	for(i=0;i<nFreqs;i++)
	{
		mol->vibration.modes[i].frequency = frequencies[i];
		mol->vibration.modes[i].mass = effectiveMasses[i];
		for(c=0;c<3;c++)
		{
			for(j=0;j<mol->nAtoms;j++) mol->vibration.modes[i].vectors[c][j] = modes[3*j+c][i];
		}
	}
	mol->vibration.nModes = nFreqs;
	sortFrequencies(mol);
	/* free tables */
	if(hessian) for(j=0;j<nFreqs;j++) if(hessian[j]) free(hessian[j]);
	if(hessian) free(hessian);
	if(modes) for(j=0;j<nFreqs;j++) if(modes[j]) free(modes[j]);
	if(modes) free(modes);
	if(frequencies) free(frequencies);
	return TRUE;
}
/********************************************************************************/
static boolean readVibrationFromOrcaHessianFile(Molecule* mol, char *fileName)
{
 	FILE *file;


	initVibrations(mol, 3*mol->nAtoms,5);
 	file = fopen(fileName, "rb");
	if(!file) return FALSE;

	if(readOrcaHessianFile_Hessian(mol,file))
	{
		rewind(file);
		readOrcaHessianFile_IR(mol, file, "$ir_spectrum");
		readOrcaHessianFile_IR(mol, file, "$raman_spectrum");
	}
	if(mol->vibration.nModes<1) 
	{
 		char t[BSIZE];
		freeVibrations(mol);
		sprintf(t,"Sorry, I can not read frequencies from '%s' file\n",fileName);
  		fprintf(stderr,"%s\n",t);
		fclose(file);
		return FALSE;
	}
	sortFrequencies(mol);
	removeTransRotModes(mol);
	fclose(file);
	return TRUE;
}
/********************************************************************************/
static Molecule* getMoleculeFromOrcaHessianFile(char *fileName)
{
	Molecule* mol = NULL;
 	char *t;
 	char *AtomCoord[5];
 	FILE *file;
 	int i;
 	int j=0;
 	int l;
	double mass;
	int nAtoms = 0;
	int totalCharge = 0;
	int spinMultiplicity = 1;
 	char t1[20];
 	char t2[20];
 	char t3[20];
 	char t4[20];
	int dum;

 	for(i=0;i<5;i++) AtomCoord[i]=malloc(BSIZE*sizeof(char));
  
 	file = fopen(fileName, "rb");
	t = malloc(BSIZE*sizeof(char));

	while(!feof(file))
	{
		if(!fgets(t,BSIZE,file))break;
		if(strstr( t, "$atoms"))
		{
			if(!fgets(t,BSIZE,file))break;
			sscanf(t,"%d",&nAtoms);
			break;
		}
 		if (strstr( t,"Total Charge") && strstr( t,"...") ) if(5==sscanf(t,"%s %s %s %s %d",t1,t2,t3,t4,&dum)) { totalCharge = dum;}
 		if (strstr( t,"Multiplicity") && strstr( t,"...") ) if(4==sscanf(t,"%s %s %s %d",t1,t2,t3,&dum)) { spinMultiplicity = dum;}
	}
	if(nAtoms == 0)
	{
 		fclose(file);
 		free(t);
 		for(i=0;i<5;i++) free(AtomCoord[i]);
		return NULL;
	}
	mol = newMolecule();
	mol->nAtoms = nAtoms;
	mol->atoms = malloc(mol->nAtoms*sizeof(Atom));
 	for(j=0;j<nAtoms;j++)
  	{
		if(!fgets(t,BSIZE,file))break;
   		sscanf(t,"%s %lf %s %s %s",AtomCoord[0],&mass,AtomCoord[1],AtomCoord[2],AtomCoord[3]);
		/* printf("t=%s\n",t);*/
		AtomCoord[0][0]=toupper(AtomCoord[0][0]);
 		l=strlen(AtomCoord[0]);
       		if (l==2) 
		{
			AtomCoord[0][1]=tolower(AtomCoord[0][1]);
			if(isdigit(AtomCoord[0][1]))l=1;
		}
		if(l==1)sprintf(t,"%c",AtomCoord[0][0]);
	        else sprintf(t,"%c%c",AtomCoord[0][0],AtomCoord[0][1]);

		mol->atoms[j].prop = propAtomGet(t);
                mol->atoms[j].mmType=strdup(t);
                mol->atoms[j].pdbType=strdup(t);
                mol->atoms[j].residueName=strdup(t);
                mol->atoms[j].N=j+1;
                mol->atoms[j].layer=HIGH_LAYER;
                mol->atoms[j].variable=TRUE;
                mol->atoms[j].show=TRUE;
		mol->atoms[j].residueNumber=0;
                mol->atoms[j].charge=0.0;
                mol->atoms[j].charge0=0.0;
                mol->atoms[j].electronegativity=0.0;
                mol->atoms[j].hardness=0.0;
	   	mol->atoms[j].width=mol->atoms[j].prop.covalentRadii;
                mol->atoms[j].mass=mass;
                mol->atoms[j].rho=0.0;
                mol->atoms[j].U = 0.0;

                mol->atoms[j].coordinates[0]=atof(AtomCoord[1])*BOHRTOANG;
                mol->atoms[j].coordinates[1]=atof(AtomCoord[2])*BOHRTOANG;
                mol->atoms[j].coordinates[2]=atof(AtomCoord[3])*BOHRTOANG;

                mol->atoms[j].gradient[0]=0.0;
                mol->atoms[j].gradient[1]=0.0;
                mol->atoms[j].gradient[2]=0.0;

                mol->atoms[j].velocity[0]=0.0;
                mol->atoms[j].velocity[1]=0.0;
                mol->atoms[j].velocity[2]=0.0;
  	}
        mol->dipole[0]=0.0;
        mol->dipole[1]=0.0;
        mol->dipole[2]=0.0;
	rewind(file);
 	get_dipole_from_orca_output_file(file, mol->dipole);
	rewind(file);
	getChargesFromOrcaOutputFile(mol,file);
 	fclose(file);
 	free(t);
 	for(i=0;i<5;i++) free(AtomCoord[i]);
        for(i=0; i<mol->nAtoms; i++) mol->atoms[i].typeConnections = NULL;
	resetTypeConnections(mol);
        int connections = 1;
        if(connections) setConnections(mol);
        else
        {
                createBondedMatrix(mol);
                printf(("Establishing connectivity : non bonded ...\n"));
                setNonBondedConnections(mol);
                freeBondedMatrix(mol);
        }

        /* if all freezed, set all to variable */
        {
                int j = 0;
                for(i=0;i<mol->nAtoms;i++)
                        if(!mol->atoms[i].variable) j++;
                if(j==mol->nAtoms)
                for(i=0;i<mol->nAtoms;i++)
                        mol->atoms[i].variable = TRUE;
        }
	mol->totalCharge = totalCharge;
	mol->spinMultiplicity = spinMultiplicity;
	return mol;
}
/********************************************************************************/
Molecule* readMoleculeFromOrcaHessianFile(char *fileName)
{
	Molecule* mol = getMoleculeFromOrcaHessianFile(fileName);
	readVibrationFromOrcaHessianFile(mol, fileName);
	readDipolesDerivativesFromOrcaHessianFile(mol, fileName);
	return mol;
	
}
/********************************************************************************/
void readGeomFromGaussianOutputFile(Molecule* mol, char *fileName, int numgeometry)
{
	char *t;
	boolean	OK;
	char *AtomCoord[5];
	FILE *file;
	int taille=BSIZE;
	int idummy;
	int i;
	int j=0;
	/* int l;*/
	int numgeom;
	char *pdest;
	int result;
	int itype=0;
	char* strStandard = "Standard orientation:";
	char* strInput = "Input orientation:";
	char* strOther = "orientation:";
	char* strSearch = strOther;

	for(i=0;i<5;i++) AtomCoord[i]=malloc(taille*sizeof(char));
	 
	t=malloc(taille*sizeof(char));
	file = fopen(fileName, "r");
	if(file ==NULL)
	{
		free(t);
	 	t = strdup_printf(("Sorry\nI can not open %s  file "),fileName);
	 	printf("%s\n",t);
	 	free(t);
	 	exit(1);
	}
	while(!feof(file))
	{
		 if(!fgets(t,taille,file))break;
	         if(strstr( t, strStandard))
		 {
			 strSearch = strStandard;
			 break;
		 }
	         if(strstr( t, strInput)) strSearch = strInput;
	}
	fseek(file, 0L, SEEK_SET);
	numgeom =0;
	do 
	{
		OK=FALSE;
		while(!feof(file))
		{
		 	fgets(t,taille,file);
			/*
		 	if(strstr(t,"Charge =") && strstr(t,"Multiplicity ="))
		 	{
				 char* p = strstr(t,"Charge =")+8;
			 	TotalCharges[0] = atoi(p);
			 	p = strstr(t,"Multiplicity =")+14;
			 	SpinMultiplicities[0] = atoi(p);
		 	}
			*/
	         	pdest = strstr( t, strSearch);
	         	result = pdest - t ;
			if ( result >0 )
		 	{
		 		fgets(t,taille,file);
		 		fgets(t,taille,file);
		 		fgets(t,taille,file);
	               		pdest = strstr( t, "Type" );
	               		result = pdest - t ;
	               		if(result>0) itype=1;
	               		else itype=0;
		 		fgets(t,taille,file);
	               		numgeom++;
	               		if(numgeom == numgeometry )
				{
					OK = TRUE;
		 			break;
				}
				OK = TRUE;
				break;
			}
		}
		if(!OK && (numgeom == 0) )
		{
	 		free(t);
	 		t = strdup_printf(("Sorry\nI can not read geometry in  %s  file "),fileName);
	 		printf("%s\n",t);
	 		free(t);
	 		return;
	   	}
	 	if(!OK)break;

	 	j=-1;
	 	while(!feof(file) )
	 	{
	   		fgets(t,taille,file);
	   		pdest = strstr( t, "----------------------------------" );
	   		result = pdest - t ;
	   		if ( result >0 )
	   		{
				/*
				long geomposok = ftell(file);
	     			get_dipole_from_gaussian_output_file(file);
				fseek(file, geomposok, SEEK_SET);
				get_charges_from_gaussian_output_file(file,j+1);
				get_natural_charges_from_gaussian_output_file(file,j+1);
				fseek(file, geomposok, SEEK_SET);
				get_esp_charges_from_gaussian_output_file(file,j+1);
				*/
	     			break;
	   		}
	   		j++;

	   		if(itype==0) sscanf(t,"%d %s %s %s %s",&idummy,AtomCoord[0],AtomCoord[1],AtomCoord[2],AtomCoord[3]);
	   		else sscanf(t,"%d %s %d %s %s %s",&idummy,AtomCoord[0],&idummy,AtomCoord[1],AtomCoord[2],AtomCoord[3]);
			/* to do : test symbol */
			/*
			AtomCoord[0][0]=toupper(AtomCoord[0][0]);
			l=strlen(AtomCoord[0]);
	         	if (l==2) AtomCoord[0][1]=tolower(AtomCoord[0][1]);
			*/


			mol->atoms[j].coordinates[0]=atof(AtomCoord[1]);
			mol->atoms[j].coordinates[1]=atof(AtomCoord[2]);
			mol->atoms[j].coordinates[2]=atof(AtomCoord[3]);

	   		mol->atoms[j].charge=0.0;
	   		mol->atoms[j].charge0=0.0;
	   		mol->atoms[j].electronegativity=0.0;
	   		mol->atoms[j].hardness=0.0;
	   		mol->atoms[j].width=mol->atoms[j].prop.covalentRadii;
	   		mol->atoms[j].rho=0.0;
	   		mol->atoms[j].U = 0.0;
		}
		if(OK && numgeometry>-1) break;
	}while(!feof(file));

	fclose(file);
	free(t);
	for(i=0;i<5;i++) free(AtomCoord[i]);
}
/********************************************************************************/
static boolean readVibrationFromMopacAuxFile(Molecule* mol, char* fileName)
{
	char** freqs = NULL;
	int nFreqs = 0;
	int numberOfFrequencies = 0;
	char** symmetries = NULL;
	int nSymmetries = 0;
	char** modes = NULL;
	int nModes = 0;
	char** intensities = NULL;
	char** effectiveMass = NULL;
	int nIntensities = 0;
	int nEffectiveMass = 0;
	int i,j, c, im;
	FILE* file;
	int nProps = 0;

	if(!mol) return FALSE;
	if(mol->nAtoms<1) return FALSE;
	if(!fileName) return FALSE; 
	file = fopen(fileName,"rb");
	if(!file) { return FALSE; }
	nModes = 0;

	rewind(file);
	freqs = get_one_block_from_aux_mopac_file(file, "VIB._FREQ:CM(-1)[",  &nFreqs);
	/* numberOfFrequencies = nFreqs-6;*/
	numberOfFrequencies = nFreqs;
	if(!freqs || numberOfFrequencies <mol->nAtoms*3)
	{
		char buffer[BSIZE];
		free_one_string_table(freqs, nFreqs);
		sprintf(buffer,"Sorry, I can not read the frequencies from '%s' file\n",fileName);
  		fprintf(stderr,"%s\n",buffer);
		fclose(file);
		return FALSE;
	}
	rewind(file);
	symmetries = get_one_block_from_aux_mopac_file(file, "NORMAL_MODE_SYMMETRY_LABELS[",  &nSymmetries);
	if(!symmetries || nSymmetries<mol->nAtoms*3)
	{
		nSymmetries = nFreqs;
		symmetries = malloc(nSymmetries*sizeof(char*));
		for(i=0;i<nSymmetries;i++)
		{
			symmetries[i] = strdup("UNK");
		}
	}
	rewind(file);
	intensities = get_one_block_from_aux_mopac_file(file, "VIB._T_DIP:ELECTRONS[",  &nIntensities);
	if(!intensities || nIntensities<mol->nAtoms*3)
	{
		nIntensities = nFreqs;
		intensities = malloc(nIntensities*sizeof(char*));
		for(i=0;i<nIntensities;i++)
		{
			intensities[i] = strdup("0.0");
		}
	}
	rewind(file);
	effectiveMass = get_one_block_from_aux_mopac_file(file, "VIB._EFF_MASS:AMU[",  &nEffectiveMass);
	if(!effectiveMass || nEffectiveMass<mol->nAtoms*3)
	{
		nEffectiveMass = nFreqs;
		effectiveMass = malloc(nEffectiveMass*sizeof(char*));
		for(i=0;i<nEffectiveMass;i++)
		{
			effectiveMass[i] = strdup("1.0");
		}
	}

	rewind(file);
	modes = get_one_block_from_aux_mopac_file(file, "NORMAL_MODES[",  &nModes);
	if(!modes || nModes<mol->nAtoms*3*numberOfFrequencies)
	{
		char buffer[BSIZE];
		free_one_string_table(freqs, nFreqs);
		free_one_string_table(symmetries, nSymmetries);
		free_one_string_table(intensities, nIntensities);
		free_one_string_table(effectiveMass, nEffectiveMass);
		free_one_string_table(modes, nModes);
		sprintf(buffer,"Sorry, I can not read the modes of frequencies from '%s' file\n",fileName);
  		fprintf(stderr,"%s\n",buffer);
		fclose(file);
		exit(1);
		return FALSE;
	}
	nModes = numberOfFrequencies;
	nProps = 1;

	initVibrations(mol, nModes, nProps);
	for(i=0;i<mol->vibration.nModes;i++) 
	for(j=0;j<mol->nAtoms;j++) 
	for(c=0;c<3;c++) mol->vibration.modes[i].vectors[c][j] = 0.0; 

	im = 0;
	for(i=0;i<mol->vibration.nModes;i++)
	{
		mol->vibration.modes[i].frequency = atof(freqs[i]);
		mol->vibration.modes[i].properties[0] = atof(intensities[i]);
		mol->vibration.modes[i].mass = atof(effectiveMass[i]);
		for(j=0;j<mol->nAtoms;j++) 
		{
			for(c=0;c<3;c++)
			{
				mol->vibration.modes[i].vectors[c][j] = atof(modes[im]);
				im++;
			}
		}
	}
	sortFrequencies(mol);
	removeTransRotModes(mol);
	free_one_string_table(freqs, nFreqs);
	free_one_string_table(symmetries, nSymmetries);
	free_one_string_table(intensities, nIntensities);
	free_one_string_table(effectiveMass, nEffectiveMass);
	free_one_string_table(modes, nModes);
	return TRUE;
}
/********************************************************************************/
static void readGeomFromMopacAuxFile(Molecule* mol, char *fileName,int numgeometry)
{
	char *t;
	boolean ok;
	char *AtomCoord[5];
	FILE *file;
	int i;
	int j=0;
	int l;
	int numgeom;
	char *pdest;
	long int geomposok = 0;
	char** elements = NULL;
	int nElements = 0;
	char** nuclearCharges = NULL;
	int nNuclearCharges = 0;
	char** partialCharges = NULL;
	int nPartialCharges = 0;

	if(!mol) {exit(1); return;}
	for(i=0;i<5;i++) AtomCoord[i] = malloc(BSIZE*sizeof(char));
	 
	t=malloc(BSIZE*sizeof(char));
	file = fopen(fileName, "rb");
	if(file ==NULL)
	{
	 	free(t);
	 	t = strdup_printf("Sorry\nI can not open %s file ",fileName);
		fprintf(stderr,"Error : %s\n",t);
	 	free(t);
		exit(1);
	 	return;
	}
	ok=FALSE;
	elements = get_one_block_from_aux_mopac_file(file, "ATOM_EL[",  &nElements);
	if(elements) ok = TRUE;
	if(!ok) 
	{
	 	free(t);
	 	t = strdup_printf("Sorry\nI can not read the atom symbols in %s file ",fileName);
		fprintf(stderr,"Error : %s\n",t);
	 	free(t);
		fclose(file);
		exit(1);
	 	return;
	 }
 	geomposok = ftell(file);
	nuclearCharges = get_one_block_from_aux_mopac_file(file, "ATOM_CORE[",  &nNuclearCharges);
	if(!nuclearCharges) fseek(file, geomposok, SEEK_SET);

	numgeom =0;
	 while(!feof(file))
	 {
		if(!fgets(t,BSIZE,file))break;
		if(numgeometry<0) pdest = strstr( t, "ATOM_X_OPT:ANGSTROMS");
		else pdest = strstr( t, "ATOM_X_UPDATED:ANGSTROMS");
		if (pdest)
		{
			numgeom++;
			geomposok = ftell(file);
			if(numgeom == numgeometry )
			{
				ok = TRUE;
				break;
			}
			if(numgeometry<0)
			{
				ok = TRUE;
			}
		}
	 }
	 if(numgeom == 0)
	 {
		free_one_string_table(elements, nElements);
		free(t);
		t = strdup_printf("Sorry\nI can not read geometry in %s file ",fileName);
		fprintf(stderr,"Error : %s\n",t);
		free(t);
		fclose(file);
		exit(1);
		return;
	  }

	mol->klass->free(mol);
	j=-1;
	fseek(file, geomposok, SEEK_SET);
	while(!feof(file) )
	{
		if(!fgets(t,BSIZE,file))break;
		if(strstr( t, "[")
		  || strstr(t,"HEAT_OF_FORM_UPDATED")
		  || strstr( t, "####################################")
		  || isABackspace(t) )
		{
			break;
		}
		if(j+1>nElements)break;
		j++;

		sscanf(t,"%s %s %s",AtomCoord[1],AtomCoord[2],AtomCoord[3]);
		if(j<nElements) sprintf(AtomCoord[0],"%s",elements[j]);
		else sprintf(AtomCoord[0],"X");
		AtomCoord[0][0]=toupper(AtomCoord[0][0]);
		l=strlen(AtomCoord[0]); 
		if (l==2) AtomCoord[0][1]=tolower(AtomCoord[0][1]);
		if(l==1)sprintf(t,"%c",AtomCoord[0][0]);
		else sprintf(t,"%c%c",AtomCoord[0][0],AtomCoord[0][1]);

	        mol->atoms = realloc(mol->atoms,(j+1)*sizeof(Atom));

                mol->atoms[j].prop = propAtomGet(t);
                mol->atoms[j].mmType=strdup(t);
                mol->atoms[j].pdbType=strdup(t);
                mol->atoms[j].residueName=strdup(t);
                mol->atoms[j].N=j+1;
                mol->atoms[j].layer=HIGH_LAYER;
                mol->atoms[j].variable=TRUE;
                mol->atoms[j].show=TRUE;
		mol->atoms[j].residueNumber=0;
                mol->atoms[j].charge=0.0;
                mol->atoms[j].charge0=0.0;
                mol->atoms[j].electronegativity=0.0;
                mol->atoms[j].hardness=0.0;
	   	mol->atoms[j].width=mol->atoms[j].prop.covalentRadii;
                mol->atoms[j].mass=mol->atoms[j].prop.mass;
                mol->atoms[j].rho=0.0;
                mol->atoms[j].U = 0.0;

                mol->atoms[j].coordinates[0]=atof(AtomCoord[1]);
                mol->atoms[j].coordinates[1]=atof(AtomCoord[2]);
                mol->atoms[j].coordinates[2]=atof(AtomCoord[3]);

                mol->atoms[j].gradient[0]=0.0;
                mol->atoms[j].gradient[1]=0.0;
                mol->atoms[j].gradient[2]=0.0;

                mol->atoms[j].velocity[0]=0.0;
                mol->atoms[j].velocity[1]=0.0;
                mol->atoms[j].velocity[2]=0.0;
		if(nuclearCharges && nNuclearCharges>j) mol->atoms[j].charge0 = atof(nuclearCharges[j]);
	}
        mol->nAtoms = j+1;

	if(numgeometry<0)
	{
		fseek(file, geomposok, SEEK_SET);
		partialCharges = get_one_block_from_aux_mopac_file(file, "ATOM_CHARGES[",  &nPartialCharges);
		if(partialCharges)
		{
    			for(j=0;j<mol->nAtoms;j++) 
				if(j<nPartialCharges) mol->atoms[j].charge = atof(partialCharges[j]);
			free_one_string_table(partialCharges, nPartialCharges);
		}
	}
	get_dipole_from_mopac_aux_file(file, mol->dipole);
	fclose(file);
	free(t);

	free_one_string_table(elements, nElements);
        free_one_string_table(nuclearCharges, nNuclearCharges);

	for(i=0;i<5;i++) free(AtomCoord[i]);
        for(i=0; i<mol->nAtoms; i++) mol->atoms[i].typeConnections = NULL;
	resetTypeConnections(mol);
        int connections = 1;
        if(connections) setConnections(mol);
        else
        {
                createBondedMatrix(mol);
                printf(("Establishing connectivity : non bonded ...\n"));
                setNonBondedConnections(mol);
                freeBondedMatrix(mol);
        }

        /* if all freezed, set all to variable */
        {
                int j = 0;
                for(i=0;i<mol->nAtoms;i++) if(!mol->atoms[i].variable) j++;
                if(j==mol->nAtoms) for(i=0;i<mol->nAtoms;i++) mol->atoms[i].variable = TRUE;
        }
}
/********************************************************************************/
int get_connections_one_atom(char* t, int nAtoms, int ibeg, int* connections)
{
	int k;
	int nc;
	int nj;
	char** ssplit = NULL;
	int nA = 0;
	/* int ibeg = 12;*/
	for(k=0;k<nAtoms;k++) connections[k] = 0;
	ssplit = split(t);
	nA = 0;
	while(ssplit && ssplit[nA]!=NULL) nA++;
	if(nA<ibeg)
	{
		strfreev(ssplit);
		return 0;
	}
	nc = atoi(ssplit[ibeg-1]);
	for(k=0;k<2*nc;k+=2) 
	{
		if(!ssplit[ibeg+k]) break;
		if(!ssplit[ibeg+k+1]) break;
		nj = atoi(ssplit[ibeg+k]);
		connections[nj-1] = atoi(ssplit[ibeg+k+1]);
	}

	strfreev(ssplit);

	return 1;
}
/********************************************************************************/
int get_gradients_one_atom(char* t, double G[])
{
	int k;
	for(k=0;k<3;k++) G[k] = 0;
	char tmp[BSIZE];
	char s[BSIZE];
	
	sprintf(tmp,"%s",t);
	uppercase(tmp);
	if(!strstr(tmp,"GRAD")) return 0;
	k = sscanf(strstr(tmp,"GRAD"),"%s %lf %lf %lf",s, &G[0], &G[1], &G[2]);

	if(k!=4) return 0;
	return 1;
}
/********************************************************************************/
int save_gradients_one_atom(FILE* file, double G[])
{
	if(!file) return 0;
	fprintf(file," GRADIENT %0.14f %0.14f %0.14f",G[0], G[1], G[2]);
	return 1;
}
/********************************************************************************/
Molecule* readMoleculeFromMopacAuxFile(char *fileName, int numgeometry)
{
	Molecule* mol = NULL;

	mol = newMolecule();
	readGeomFromMopacAuxFile(mol, fileName,numgeometry);
	readVibrationFromMopacAuxFile(mol, fileName);
	return mol;
}
/********************************************************************************/
Molecule* readMoleculeFromMopacOutputFile(char *fileName, int numgeometry)
{
/*      Change only the coordinates of mol */
	char* t;
	boolean OK;
	char *AtomCoord[5];
	FILE *file;
	int idummy;
	int i;
	int j=0;
	int l;
	int numgeom;
	char *pdest;
	long int geomposok = 0;
	Molecule* mol = NULL;

	for(i=0;i<5;i++) AtomCoord[i]=malloc(BSIZE*sizeof(char));
	t=malloc(BSIZE*sizeof(char));
	 
	file = fopen(fileName, "r");
	if(file ==NULL)
	{
	 	free(t);
	 	printf(("Sorry\nI can not open %s  mopac output file\n"),fileName);
		exit(1);
	 	return NULL;
	}
	numgeom =0;
	OK=FALSE;
	 while(!feof(file))
	 {
		if(!fgets(t,BSIZE,file))break;
		//pdest = strstr( t, "CARTESIAN COORDINATES");
		pdest = strstr( t, " ATOM   CHEMICAL          X               Y               Z");
		if(pdest) 
		{
			if(!fgets(t,BSIZE,file)) {pdest=0;break;}
			if(!fgets(t,BSIZE,file)) {pdest=0;break;}
		}
		if ( pdest )
		{
			numgeom++;
			geomposok = ftell(file);
			if(numgeom == numgeometry )
			{
				OK = TRUE;
				break;
			}
			if(numgeometry<0)
			{
				OK = TRUE;
			}
		}
	 }
	 if(!OK || numgeom == 0)
	 {
		free(t);
	 	printf(("Sorry\nI can not open %s mopac output file\n"),fileName);
		exit(1);
	 	return NULL;
	  }
	j=-1;
	fseek(file, geomposok, SEEK_SET);
	mol = newMolecule();
	while(!feof(file) )
	{
		if(!fgets(t,BSIZE,file))break;
		if(isABackspace(t))
		{
			break;
		}
		j++;
		for(i=0;i<strlen(t);i++) if (t[i]=='*') t[i] = ' ';
		sscanf(t,"%d %s %s %s %s",&idummy,AtomCoord[0],AtomCoord[1],AtomCoord[2],AtomCoord[3]);
		AtomCoord[0][0]=toupper(AtomCoord[0][0]);
		l=strlen(AtomCoord[0]); 
		if(isdigit(AtomCoord[0][1]))l=1;
		if (l==2) AtomCoord[0][1]=tolower(AtomCoord[0][1]);
		if(l==1)sprintf(t,"%c",AtomCoord[0][0]);
		else sprintf(t,"%c%c",AtomCoord[0][0],AtomCoord[0][1]);
		/* test symbol to do */

		mol->atoms = realloc(mol->atoms,(j+1)*sizeof(Atom));

		mol->atoms[j].prop = propAtomGet(t);
		mol->atoms[j].mmType=strdup(t);
		mol->atoms[j].pdbType=strdup(t);
		mol->atoms[j].residueName=strdup(t);
		mol->atoms[j].N=j+1;
		mol->atoms[j].layer=HIGH_LAYER;
		mol->atoms[j].variable=TRUE;
		mol->atoms[j].show=TRUE;
		mol->atoms[j].residueNumber=0;
	   	mol->atoms[j].charge=0.0;
	   	mol->atoms[j].charge0=0.0;
	   	mol->atoms[j].electronegativity=0.0;
	   	mol->atoms[j].hardness=0.0;
	   	mol->atoms[j].width=mol->atoms[j].prop.covalentRadii;
	   	mol->atoms[j].mass=mol->atoms[j].prop.mass;
	   	mol->atoms[j].rho=0.0;
	   	mol->atoms[j].U = 0.0;

		mol->atoms[j].coordinates[0]=atof(AtomCoord[1]);
		mol->atoms[j].coordinates[1]=atof(AtomCoord[2]);
		mol->atoms[j].coordinates[2]=atof(AtomCoord[3]);

		mol->atoms[j].gradient[0]=0.0;
		mol->atoms[j].gradient[1]=0.0;
		mol->atoms[j].gradient[2]=0.0;

		mol->atoms[j].velocity[0]=0.0;
		mol->atoms[j].velocity[1]=0.0;
		mol->atoms[j].velocity[2]=0.0;
	  }
	  mol->nAtoms = j+1;

	 fclose(file);
	 free(t);
	 for(i=0;i<5;i++) free(AtomCoord[i]);
        for(i=0; i<mol->nAtoms; i++) mol->atoms[i].typeConnections = NULL;
	resetTypeConnections(mol);
	int connections = 1;
	if(connections) setConnections(mol);
	else
	{
		createBondedMatrix(mol);
		printf(("Establishing connectivity : non bonded ...\n"));
		setNonBondedConnections(mol);
		freeBondedMatrix(mol);
	}

	/* if all freezed, set all to variable */
	{
		int j = 0;
		for(i=0;i<mol->nAtoms;i++)
			if(!mol->atoms[i].variable) j++;
		if(j==mol->nAtoms)
		for(i=0;i<mol->nAtoms;i++)
			mol->atoms[i].variable = TRUE;
	}
	return mol;
}
/********************************************************************************/
static boolean readVibrationFromGamessOutputFile(Molecule* mol, char* fileName)
{
 	char t[BSIZE];
 	boolean ok;
	int i;
	int j;
	int c;
	int ne;
	int nf;
	int nir;
	int nmass = 0;
	int nfMax = 5;
	double freq[5];
	double ir[5];
	double raman[5];
	double mass[5];
 	char* sdum[5*2];
 	char* tmp;
	int k;
	int nProps;
	FILE* file;

	if(!mol) return FALSE;
        if(mol->nAtoms<1) return FALSE;
        if(!fileName) return FALSE;
        file = fopen(fileName,"rb");
        if(!file) { return FALSE; }

	nProps = 2;

 	ok=FALSE;
 	while(!feof(file))
	{
		if(!fgets(t,BSIZE,file))break;
	 	if ( strstr( t,"FREQUENCY:") )
	  	{
			ok = TRUE;
			break;
	  	}
	}

	if(!ok) return 1;

	initVibrations(mol, 3*mol->nAtoms, nProps);

	for(i=0;i<nfMax*2;i++) sdum[i] = malloc(BSIZE*sizeof(char));

	j = 0;
  	while(!feof(file))
  	{
		int nfi=0;
		if(!strstr( t,"FREQUENCY:")) break;

		tmp = strstr(t,":")+1;
		for(i=0;i<nfMax*2;i++) sprintf(sdum[i]," ");
		nfi = sscanf(tmp,"%s %s %s %s %s %s %s %s %s %s", sdum[0],sdum[1],sdum[2],sdum[3],sdum[4],
				sdum[5],sdum[6],sdum[7],sdum[8],sdum[9]
				);
		if(nfi<1)
		{
			mol->vibration.nModes = j;
                	mol->vibration.modes = realloc(mol->vibration.modes, mol->vibration.nModes*sizeof(VibMode));
			for(i=0;i<nfMax*2;i++) free(sdum[i]);
			return 2;
		}
		nf = 0;
		for(i=0;i<nfi;)
		{
			if(strstr(sdum[i+1],"I"))
			{
				freq[nf] = -atof(sdum[i]);
				i+=2;
			}
			else
			{
				freq[nf] = atof(sdum[i]);
				i+=1;
			}
			nf++;
		}
		nir=-1;
		for(i=0;i<nfMax;i++) ir[i] = 0;
		for(i=0;i<nfMax;i++) raman[i] = 0;
		while(fgets(t,BSIZE,file) && strstr(t,":")) /* REDUCED MASS: IR INTENSITY: RAMAN ACTIVITY: Depol,... backspace */
		{
			if(strstr(t,"MASS:"))
			{
				tmp =  strstr(t,":")+1;
				nmass = sscanf(tmp,"%s %s %s %s %s", sdum[0],sdum[1],sdum[2],sdum[3],sdum[4]);
				if(nf!=nmass)
				{
					mol->vibration.nModes = j;
                			mol->vibration.modes = realloc(mol->vibration.modes, mol->vibration.nModes*sizeof(VibMode));
					for(i=0;i<nfMax*2;i++) free(sdum[i]);
					return 2;
				}
				for(i=0;i<nf;i++) mass[i] = atof(sdum[i]);
			}
			if(strstr(t,"IR"))
			{
				tmp =  strstr(t,":")+1;
				nir = sscanf(tmp,"%s %s %s %s %s", sdum[0],sdum[1],sdum[2],sdum[3],sdum[4]);
				if(nf!=nir)
				{
					mol->vibration.nModes = j;
                			mol->vibration.modes = realloc(mol->vibration.modes, mol->vibration.nModes*sizeof(VibMode));
					for(i=0;i<nfMax*2;i++) free(sdum[i]);
					return 2;
				}
				for(i=0;i<nf;i++) ir[i] = atof(sdum[i]);
			}
			if(strstr(t,"RAMAN"))
			{
				tmp =  strstr(t,":")+1;
				nir = sscanf(tmp,"%s %s %s %s %s", sdum[0],sdum[1],sdum[2],sdum[3],sdum[4]);
				if(nf!=nir)
				{
					mol->vibration.nModes = j;
                			mol->vibration.modes = realloc(mol->vibration.modes, mol->vibration.nModes*sizeof(VibMode));
					for(i=0;i<nfMax*2;i++) free(sdum[i]);
					return 2;
				}
				for(i=0;i<nf;i++) raman[i] = atof(sdum[i]);
			}
		}
		for(i=0;i<nf;i++)
		{
			mol->vibration.modes[j].frequency = freq[i];
			mol->vibration.modes[j].properties[0] = ir[i];
			mol->vibration.modes[j].mass = mass[i];
			mol->vibration.modes[j].properties[1] = raman[i];
			j++;
		}
		for(i=0;i<mol->nAtoms;i++)
		{
			k = j-nf;
    			if(!fgets(t,BSIZE,file)) break;
			c = 0;
			ne = sscanf(t,"%s %s %s %lf %lf %lf %lf %lf %lf", sdum[0],sdum[1],sdum[2],
				&mol->vibration.modes[k  ].vectors[c][i],
				&mol->vibration.modes[k+1].vectors[c][i],
				&mol->vibration.modes[k+2].vectors[c][i],
				&mol->vibration.modes[k+3].vectors[c][i],
				&mol->vibration.modes[k+4].vectors[c][i],
				&mol->vibration.modes[k+5].vectors[c][i]);
			if(ne!=nf+3)return 2;
    			if(!fgets(t,BSIZE,file)) break;
			c = 1;
			ne = sscanf(t," %s %lf %lf %lf %lf %lf %lf",sdum[0],
				&mol->vibration.modes[k  ].vectors[c][i],
				&mol->vibration.modes[k+1].vectors[c][i],
				&mol->vibration.modes[k+2].vectors[c][i],
				&mol->vibration.modes[k+3].vectors[c][i],
				&mol->vibration.modes[k+4].vectors[c][i],
				&mol->vibration.modes[k+5].vectors[c][i]);
			if(ne!=nf+1)return 2;
    			if(!fgets(t,BSIZE,file)) break;
			c = 2;
			ne = sscanf(t," %s %lf %lf %lf %lf %lf %lf",sdum[0],
				&mol->vibration.modes[k  ].vectors[c][i],
				&mol->vibration.modes[k+1].vectors[c][i],
				&mol->vibration.modes[k+2].vectors[c][i],
				&mol->vibration.modes[k+3].vectors[c][i],
				&mol->vibration.modes[k+4].vectors[c][i],
				&mol->vibration.modes[k+5].vectors[c][i]);
			if(ne!=nf+1)return 2;
		}
		for(i=0;i<5*2+2;i++)
		{
			if(!fgets(t,BSIZE,file))break;
		}
		if(i!=5*2+2 || !fgets(t,BSIZE,file))
		{
			mol->vibration.nModes = j;
                	mol->vibration.modes = realloc(mol->vibration.modes, mol->vibration.nModes*sizeof(VibMode));
			for(i=0;i<nfMax*2;i++) free(sdum[i]);
			return 2;
		}
	}
	mol->vibration.nModes = j;
        mol->vibration.modes = realloc(mol->vibration.modes, mol->vibration.nModes*sizeof(VibMode));
	for(j=0;j< mol->vibration.nModes;j++)
	{
			for(i=0;i< mol->nAtoms;i++)
			for(c=0;c<3;c++)
			mol->vibration.modes[j].vectors[c][i] *= sqrt(mol->vibration.modes[j].mass);
	}
	for(i=0;i<nfMax*2;i++) free(sdum[i]);
	sortFrequencies(mol);
	removeTransRotModes(mol);
	return 0;
}
/********************************************************************************/
Molecule* getMoleculeFromGamessOutputFile(char *fileName, int numgeometry)
{
	char *t;
	boolean OK;
	char *AtomCoord[5];
	FILE *file;
	int i;
	int j=0;
	int l;
	int numgeom;
	char dum[100];
	Molecule* mol = NULL;


	for(i=0;i<5;i++) AtomCoord[i]=malloc(BSIZE*sizeof(char));
  
	t=malloc(BSIZE*sizeof(char));
 	file = fopen(fileName, "rb");
	if(file ==NULL)
	{
		free(t);
		printf(("Sorry\nI can not open %s Firefly (gamess) output file\n"),fileName);
		exit(1);
		return NULL;
	}
	numgeom = 0;
	do 
	{
		OK=FALSE;
		while(!feof(file)){
			fgets(t,BSIZE,file);
			if ( numgeometry==1 && strstr(t,"COORDINATES (BOHR)"))
			{
	  			fgets(t,BSIZE,file);
 				numgeom++;
				if((int)numgeom == numgeometry ) { OK = TRUE; break; }
	  		}
			if ( strstr(t,"COORDINATES OF ALL ATOMS ARE (ANGS)"))
			{
	  			fgets(t,BSIZE,file);
	  			fgets(t,BSIZE,file);
 				numgeom++;
				if((int)numgeom == numgeometry ) { OK = TRUE; break; }
				if(numgeometry<0 ) { OK = TRUE; break; }
	  		}
		}
		if(!OK && (numgeom == 0) ){
			free(t);
			printf(("Sorry\nI can not open read geometry from %s Firefly (gamess) output file\n"),fileName);
			exit(1);
			return NULL;
		}
		if(!OK)break;

		//printf("Begin freMole\n");
		if(mol) freeMolecule(mol);
		mol = newMolecule();
		mol->atoms = NULL;
		mol->nAtoms = 0;
		//printf("numgeom=%d numgeometry=%d\n",numgeom,numgeometry);

		j=-1;
		while(!feof(file) )
		{
			fgets(t,BSIZE,file);
			strDeleten(t);
			if (isABackspace(t)) break;
			if ( !strcmp(t,"\n")) break;
			if ( !strcmp(t,"\r\n")) break;
			j++;
			//printf("j=%d t = %s\n",j,t);
			sscanf(t,"%s %s %s %s %s",AtomCoord[0],dum, AtomCoord[1], AtomCoord[2],AtomCoord[3]);
			{
				int k;
				for(k=0;k<(int)strlen(AtomCoord[0]);k++) if(isdigit(AtomCoord[0][k])) AtomCoord[0][k] = ' ';
				deleteAllSpaces(AtomCoord[0]);
			}

			AtomCoord[0][0]=toupper(AtomCoord[0][0]);
			l=strlen(AtomCoord[0]);
			if (l==2) AtomCoord[0][1]=tolower(AtomCoord[0][1]);
			mol->atoms = realloc(mol->atoms,(j+1)*sizeof(Atom));

			mol->atoms[j].prop = propAtomGet(AtomCoord[0]);
			mol->atoms[j].mmType=strdup(AtomCoord[0]);
			mol->atoms[j].pdbType=strdup(AtomCoord[0]);
			mol->atoms[j].residueName=strdup(AtomCoord[0]);
			mol->atoms[j].N=j+1;
			mol->atoms[j].layer=HIGH_LAYER;
			mol->atoms[j].variable=TRUE;
			mol->atoms[j].show=TRUE;
			mol->atoms[j].residueNumber=0;
	   		mol->atoms[j].charge=0.0;
	   		mol->atoms[j].charge0=0.0;
	   		mol->atoms[j].electronegativity=0.0;
	   		mol->atoms[j].hardness=0.0;
	   		mol->atoms[j].width=mol->atoms[j].prop.covalentRadii;
	   		mol->atoms[j].mass= mol->atoms[j].prop.mass;
	   		mol->atoms[j].rho=0.0;
	   		mol->atoms[j].U = 0.0;

			mol->atoms[j].coordinates[0]=atof(AtomCoord[1]);
			mol->atoms[j].coordinates[1]=atof(AtomCoord[2]);
			mol->atoms[j].coordinates[2]=atof(AtomCoord[3]);

			mol->atoms[j].gradient[0]=0.0;
			mol->atoms[j].gradient[1]=0.0;
			mol->atoms[j].gradient[2]=0.0;

			mol->atoms[j].velocity[0]=0.0;
			mol->atoms[j].velocity[1]=0.0;
			mol->atoms[j].velocity[2]=0.0;
        		mol->atoms[j].typeConnections = NULL;
			if(j==mol->nAtoms-1) break;
		}
		if(OK && numgeometry>=0) break;
		mol->nAtoms = j+1;
	}while(!feof(file));
	//printf("end read molecule\n");

	rewind(file);
	getChargesFromGamessOutputFile(mol, file);
	fclose(file);
	free(t);
	for(i=0;i<5;i++) free(AtomCoord[i]);
	//printf("begin resetTypeConnections\n");
        for(i=0; i<mol->nAtoms; i++) mol->atoms[i].typeConnections = NULL;
	resetTypeConnections(mol);
	//printf("end resetTypeConnections\n");
	int connections = 1;
	if(connections) setConnections(mol);
	else
	{
		createBondedMatrix(mol);
		printf(("Establishing connectivity : non bonded ...\n"));
		setNonBondedConnections(mol);
		freeBondedMatrix(mol);
	}

	/* if all freezed, set all to variable */
	{
		int j = 0;
		for(i=0;i<mol->nAtoms;i++)
			if(!mol->atoms[i].variable) j++;
		if(j==mol->nAtoms)
		for(i=0;i<mol->nAtoms;i++)
			mol->atoms[i].variable = TRUE;
	}
	return mol;
}
/********************************************************************************/
Molecule* readMoleculeFromGamessOutputFile(char *fileName, int numgeometry)
{
	//printf("Begin redMolecule from gamess output file\n");
	Molecule* mol = getMoleculeFromGamessOutputFile(fileName, numgeometry);
	//printf("End redMolecule from gamess output file\n");
	readVibrationFromGamessOutputFile(mol, fileName);
	return mol;
	
}
/********************************************************************************/
static boolean readVibrationFromOrcaOutputFile(Molecule* mol, char* fileName)
{
 	char t[BSIZE];
 	char sdum1[BSIZE];
 	boolean ok;
 	FILE *file;
	int nf;
	int n;
	double freq = 0;
	double v[6] ={0,0,0,0,0,0};
	int j;
	int k;
	int nfOld;
	boolean Begin = TRUE;
	int nblock, iblock;
	int ix;
	double dum;
	int numberOfFrequencies = 0;


	initVibrations(mol, 3*mol->nAtoms,2);
 	file = fopen(fileName, "rb");
	if(!file) return FALSE;
 	do 
 	{
 		ok=FALSE;
 		while(!feof(file))
		{
    			if(!fgets(t,BSIZE,file)) break;
	 		if (strstr( t,"VIBRATIONAL FREQUENCIES") ) ok = TRUE;
	 		if (strstr( t,":") && ok ){ ok = TRUE; break;}
		}
		if(!ok) break;
  		numberOfFrequencies = 0;
  		while(!feof(file) )
  		{
			if(!strstr(t,":")) break;
			nf = sscanf(t,"%s %lf",sdum1,&freq);
			if(nf!=2) { ok = FALSE; break;}
  			numberOfFrequencies++;
			if(numberOfFrequencies>mol->nAtoms*3) break;
			k = numberOfFrequencies-1;
			mol->vibration.modes[k].frequency = freq;
			mol->vibration.modes[k].mass = 1.0;
			mol->vibration.modes[k].properties[0] = 0;
			mol->vibration.modes[k].properties[1] = 0;
			if(!fgets(t,BSIZE,file)) break;
		}
		if(!ok) break;
 		ok=FALSE;
 		while(!feof(file))
		{
    			if(!fgets(t,BSIZE,file)) break;
	 		if (strstr( t,"NORMAL MODES") ) ok = TRUE;
	 		if (sscanf(t,"%d",&k)==1 && ok ){ ok = TRUE; break;}
		}
		if(!ok) break;
		nblock = numberOfFrequencies/6; 
		if(numberOfFrequencies%6!=0) nblock++;
		for(iblock = 0;iblock<nblock;iblock++)
		{
			nf = sscanf(t,"%s %lf %lf %lf %lf %lf %lf",
					sdum1,
					&v[0],&v[1],&v[2],
					&v[3],&v[4],&v[5]
					);
			ix = 0;
			for(j=0;j<mol->nAtoms*3 && !feof(file);j++)
			{
				if(!fgets(t,BSIZE,file)) break;
				nf = sscanf(t,"%s %lf %lf %lf %lf %lf %lf",
					sdum1,
					&v[0],&v[1],&v[2],
					&v[3],&v[4],&v[5]
					);
				nf--;
				nfOld = iblock*6;
				for(k=0;k<nf;k++)
				{
					mol->vibration.modes[k+nfOld].vectors[ix][j/3]= v[k]; 
				}
				ix++;
				if(ix>2) ix = 0;
			}
			if(!fgets(t,BSIZE,file)) {ok = FALSE;break;};/* new block */
		}
		if(!ok) break;
		Begin = FALSE;
 		ok=FALSE;
 		while(!feof(file))
		{
			if(!fgets(t,BSIZE,file)) break;
	 		if (strstr( t,"IR SPECTRUM") ) ok = TRUE;
	 		if (strstr( t,"TX")  && strstr( t,"TY") && strstr( t,"TZ") && ok ){ ok = TRUE; break;}
		}
		if(!ok) {continue;}
    		if(fgets(t,BSIZE,file))
  		while(!feof(file) )
  		{
			if(!fgets(t,BSIZE,file)) break;
			n = sscanf(t,"%s %lf %lf", sdum1, &freq,&dum);
			if(n!=3) { break; }
			k = atoi(t);
			mol->vibration.modes[k].properties[0] = dum;
		}
		if(!ok) break;
 		ok=FALSE;
 		while(!feof(file))
		{
			if(!fgets(t,BSIZE,file)) break;
	 		if (strstr( t,"RAMAN SPECTRUM") ) ok = TRUE;
	 		if (strstr( t,"Activity")  && strstr( t,"Depolarization") && ok ){ ok = TRUE; break;}
		}
		if(!ok) {continue;}
    		if(fgets(t,BSIZE,file))
  		while(!feof(file) )
  		{
			if(!fgets(t,BSIZE,file)) break;
			n = sscanf(t,"%s %lf %lf", sdum1, &freq,&dum);
			if(n!=3) { break; }
			k = atoi(t);
			mol->vibration.modes[k].properties[1] = dum;
		}
		if(!ok) break;
		break;
 	}while(!feof(file));
	if((Begin && !ok) || numberOfFrequencies<1) 
	{
		char buffer[BSIZE];
		freeVibrations(mol);
		sprintf(buffer,"Sorry, I can not read frequencies from '%s' file\n",fileName);
  		fprintf(stderr,"%s",buffer);
		return FALSE;
	}
	else
	{
		char buffer[BSIZE];
		sprintf(buffer,
				"\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n"
				"Warning, the effective masses are not available in an orca ourput file\n"
				"It is if you want to generate a geometry along a vibrational mode.\n"
				"These masses are set here to 1.0\n"
				"RECOMMENDED : Use the *.hessian file given by orca\n"
				"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n"
				);
  		fprintf(stdout,"%s\n",buffer);
		mol->vibration.nModes = numberOfFrequencies;
		mol->vibration.modes = realloc(mol->vibration.modes, mol->vibration.nModes*sizeof(VibMode));
		sortFrequencies(mol);
		removeTransRotModes(mol);
	}
	return TRUE;
}
/********************************************************************************/
Molecule* getMoleculeFromOrcaOutputFile(char* fileName, int numgeometry)
{
	char *t;
	char *AtomCoord[5];
	FILE *file;
	int taille=BSIZE;
	int i;
	int j=0;
	int l;
	int numgeom;
	char *pdest;
	long int geomposok = 0;
	Molecule* mol = NULL;
	int totalCharge = 0;
	int spinMultiplicity = 1;
 	char t1[20];
 	char t2[20];
 	char t3[20];
 	char t4[20];
	int dum;

	for(i=0;i<5;i++) AtomCoord[i]=malloc(taille*sizeof(char));
	 
	t=malloc(taille*sizeof(char));
	file = fopen(fileName, "r");
	if(file ==NULL)
	{
	 	free(t);
	 	printf(("Sorry\nI can not open the %s  orca output file\n"),fileName);
	 	return NULL;
	}
	numgeom =0;
	 while(!feof(file))
	 {
		if(!fgets(t,taille,file))break;
 		if (strstr( t,"Total Charge") && strstr( t,"...") ) if(5==sscanf(t,"%s %s %s %s %d",t1,t2,t3,t4,&dum)) { totalCharge = dum;}
 		if (strstr( t,"Multiplicity") && strstr( t,"...") ) if(4==sscanf(t,"%s %s %s %d",t1,t2,t3,&dum)) { spinMultiplicity = dum;}
		pdest = strstr( t, "CARTESIAN COORDINATES (ANGSTROEM)");
		if(pdest) 
		{
			if(!fgets(t,taille,file))break;
			pdest = strstr( t, "---------------------------------");
		}
		if ( pdest )
		{
			numgeom++;
			geomposok = ftell(file);
			if(numgeom == numgeometry )
			{
				break;
			}
			if(numgeometry<0)
			{
			}
		}
	 }
	 if(numgeom == 0)
	 {
		free(t);
		t = strdup_printf(("Sorry\nI can not read geometry from %s  orca output file\n"),fileName);
		return NULL;
	  }
	j=-1;
	fseek(file, geomposok, SEEK_SET);
	mol = newMolecule();
	mol->atoms = NULL;
	mol->nAtoms = 0;


	while(!feof(file) )
	{
		if(!fgets(t,taille,file))break;
		pdest = strstr( t, "----------------------------------" );
		if (pdest || isABackspace(t))
		{
			break;
		}
		j++;
		sscanf(t,"%s %s %s %s",AtomCoord[0],AtomCoord[1],AtomCoord[2],AtomCoord[3]);
		AtomCoord[0][0]=toupper(AtomCoord[0][0]);
		l=strlen(AtomCoord[0]); 
		if(isdigit(AtomCoord[0][1]))l=1;
		if (l==2) AtomCoord[0][1]=tolower(AtomCoord[0][1]);
		if(l==1)sprintf(t,"%c",AtomCoord[0][0]);
		else sprintf(t,"%c%c",AtomCoord[0][0],AtomCoord[0][1]);

		mol->atoms = realloc(mol->atoms,(j+1)*sizeof(Atom));

		mol->atoms[j].prop = propAtomGet(AtomCoord[0]);
		mol->atoms[j].mmType=strdup(AtomCoord[0]);
		mol->atoms[j].pdbType=strdup(AtomCoord[0]);
		mol->atoms[j].residueName=strdup(AtomCoord[0]);
		mol->atoms[j].N=j+1;
		mol->atoms[j].layer=HIGH_LAYER;
		mol->atoms[j].variable=TRUE;
		mol->atoms[j].show=TRUE;
		mol->atoms[j].residueNumber=0;
	   	mol->atoms[j].charge=0.0;
	   	mol->atoms[j].charge0=0.0;
	   	mol->atoms[j].electronegativity=0.0;
	   	mol->atoms[j].hardness=0.0;
	   	mol->atoms[j].width=mol->atoms[j].prop.covalentRadii;
	   	mol->atoms[j].mass= mol->atoms[j].prop.mass;
	   	mol->atoms[j].rho=0.0;
	   	mol->atoms[j].U = 0.0;

		mol->atoms[j].coordinates[0]=atof(AtomCoord[1]);
		mol->atoms[j].coordinates[1]=atof(AtomCoord[2]);
		mol->atoms[j].coordinates[2]=atof(AtomCoord[3]);

		mol->atoms[j].gradient[0]=0;
		mol->atoms[j].gradient[1]=0;
		mol->atoms[j].gradient[2]=0;

		mol->atoms[j].velocity[0]=0.0;
		mol->atoms[j].velocity[1]=0.0;
		mol->atoms[j].velocity[2]=0.0;
	  }
	 mol->nAtoms = j+1;
	 fclose(file);
	 free(t);
	 for(i=0;i<5;i++) free(AtomCoord[i]);
        for(i=0; i<mol->nAtoms; i++) mol->atoms[i].typeConnections = NULL;
	resetTypeConnections(mol);
	int connections = 1;
	if(connections) setConnections(mol);
	else
	{
		createBondedMatrix(mol);
		printf(("Establishing connectivity : non bonded ...\n"));
		setNonBondedConnections(mol);
		freeBondedMatrix(mol);
	}

	/* if all freezed, set all to variable */
	{
		int j = 0;
		for(i=0;i<mol->nAtoms;i++)
			if(!mol->atoms[i].variable) j++;
		if(j==mol->nAtoms)
		for(i=0;i<mol->nAtoms;i++)
			mol->atoms[i].variable = TRUE;
	}
	mol->totalCharge = totalCharge;
	mol->spinMultiplicity = spinMultiplicity;
	return mol;
}
/********************************************************************************/
static boolean readDipolesDerivativesFromOrcaHessianFile(Molecule* mol, char* fileName)
{
	// Here I read derivatives along atomic coordinates and then I compute the derivative along the normal modes
 	char t[BSIZE];
 	boolean ok;
 	FILE *file;
	int nf;
	int i,j,k;
	int xyz;
	double*** T = NULL;
	int nAll = 0;
	//printf("nModes =%d\n",mol->vibration.nModes);
	//printf("fileName =%s\n",fileName);
	if(mol->vibration.nModes<1) return FALSE;
	if(mol->nAtoms<1) return FALSE;
	//printf("fileName =%s\n",fileName);
	T = malloc(3*sizeof(double**));// Dx, Dy, Dz
	for(xyz=0;xyz<3;xyz++) 
	{
		T[xyz] = malloc(3*sizeof(double*)); // 3 coordiantes for each atom
		for(j=0;j<3;j++) T[xyz][j] = malloc(mol->nAtoms*sizeof(double)); 
	}

 	file = fopen(fileName, "rb");
	if(!file) return FALSE;
	//printf("fileName =%s\n",fileName);

/*$dipole_derivatives
9
    -0.786390     0.000000     0.000000
     0.000000    -0.287454     0.000000
     0.000000     0.000000    -0.263611
     0.393191    -0.000000     0.000000
    -0.000000     0.143737     0.147331
     0.000000     0.188756     0.131853
     0.393191    -0.000000    -0.000000
    -0.000000     0.143737    -0.147331
    -0.000000    -0.188756     0.131853
*/

 	ok=FALSE;
	while(!feof(file))
	{
    		if(!fgets(t,BSIZE,file)) break;
 		if (strstr( t,"$dipole_derivatives"))
		{ 
    			if(!fgets(t,BSIZE,file)) break;
			nf=sscanf(t,"%d",&nAll);
			if(nf!=1) break;
			if(nAll!=mol->nAtoms*3) break;
			ok = TRUE; 
			break;
		}
	}
	k = 0;
	if(ok)
	for(i=0;i<mol->nAtoms;i++)
	for(j=0;j<3;j++) 
  	{
    		if(!fgets(t,BSIZE,file)) { ok = FALSE;  break;}
		nf = sscanf(t,"%lf %lf %lf",&T[0][j][i], &T[1][j][i], &T[2][j][i]);
		if(nf!=3) { ok = FALSE; break;}
		k++;
	}
	//printf("k=%d nAll = %d\n",k,nAll);
	if(k!=nAll) ok = FALSE;
	if(ok)
	{
        	double mu0 = 4*PI*1e-7;
        	double eps0 = 1.0/(mu0*slight*slight);
        	double   kmmolm1 = 4*PI*PI*PI*NAvogadro/3/hPlank/slight/4/PI/eps0*1e-3*100.0*8.47835267e-30*8.47835267e-30;/* 1e-3 m to km, 100 : cm-1 to m-1 */
		double f = sqrt(AUTOCM1/AMUTOAU);
		int nPOld = 2;
		double IRI = 0.;
		// allocation already make for dipole derivatives in readOrca
		for(i=0;i<mol->vibration.nModes;i++)
		{
			for(xyz=0;xyz<3;xyz++)
			{
				int jxyz = nPOld+xyz;
				mol->vibration.modes[i].properties[jxyz] = 0;
				for(j=0;j<3;j++) 
				for(k=0;k<mol->nAtoms;k++)
					mol->vibration.modes[i].properties[jxyz] += mol->vibration.modes[i].vectors[j][k]*T[xyz][j][k]/sqrt(mol->vibration.modes[i].mass);
			}
			IRI = 0;
			for(xyz=0;xyz<3;xyz++)
			{
				int jxyz = nPOld+xyz;
				mol->vibration.modes[i].properties[jxyz] *= f;
				IRI += mol->vibration.modes[i].properties[jxyz]*mol->vibration.modes[i].properties[jxyz];
			}
			IRI *= kmmolm1;
			mol->vibration.modes[i].properties[0] = IRI;
			/*
			for(xyz=0;xyz<3;xyz++)
			{
				int jxyz = nPOld+xyz;
				printf("xyz = %d mode = %d mu = %f\n",xyz,i,mol->vibration.modes[i].properties[jxyz]);
			}
			printf("IR (km/mol) = %f\n",mol->vibration.modes[i].properties[0]);
			*/
		}
	}
	for(xyz=0;xyz<3;xyz++) 
	{
		for(j=0;j<3;j++) free(T[xyz][j]);
		free(T[xyz]);
	}
	free(T);
	fclose(file);
	return TRUE;
}
/********************************************************************************/
static boolean readDipolesDerivativesFromGaussianOutputFile(Molecule* mol, char* fileName)
{
 	char t[BSIZE];
 	boolean ok;
 	FILE *file;
	int nf;
	int numberOfFrequencies = 0;
	double* TX = NULL;
	double* TY = NULL;
	double* TZ = NULL;
	//printf("nModes =%d\n",mol->vibration.nModes);
	//printf("fileName =%s\n",fileName);
	if(mol->vibration.nModes<1) return FALSE;
	//printf("fileName =%s\n",fileName);
	TX = malloc(mol->vibration.nModes*sizeof(double));
	TY = malloc(mol->vibration.nModes*sizeof(double));
	TZ = malloc(mol->vibration.nModes*sizeof(double));


 	file = fopen(fileName, "rb");
	if(!file) return FALSE;
	//printf("fileName =%s\n",fileName);

 	ok=TRUE;
  	numberOfFrequencies = 0;
  	while(!feof(file) )
  	{
		int i;
		if(!fgets(t,BSIZE,file)) break;
		if(strstr(t,"Harmonic frequencies") && strstr(t,"IR intensities")) break;
		if(!strstr(t,"Dipole derivatives wrt mode"))  continue;
		if(!strstr(t,":"))  continue;
		for(i=0;i<strlen(t);i++) if(t[i]=='D') t[i]='E';
		nf = sscanf(strstr(t,":")+1,"%lf %lf %lf",&TX[numberOfFrequencies], &TY[numberOfFrequencies], &TZ[numberOfFrequencies]);
		if(nf!=3) { ok = FALSE; break;}
		numberOfFrequencies++;
		if(numberOfFrequencies>mol->vibration.nModes) { ok = FALSE; break;}
	}
	if(numberOfFrequencies!=mol->vibration.nModes) ok = FALSE;
	if(!ok)
	{
			fprintf(stderr,"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n");
			fprintf(stderr," WARNING : I cannot read the dipole derivatives from %s \n",fileName);
			fprintf(stderr,"         : To print these values, add iop(7/33=1) Freq keywords to your  gaussian input file \n");
			fprintf(stderr,"         : For unknown reason iop(7/33=1) does not work with opt keyword! \n");
			fprintf(stderr,"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n");
	}
	if(ok && numberOfFrequencies==mol->vibration.nModes)
	{
        	double mu0 = 4*PI*1e-7;
        	double eps0 = 1.0/(mu0*slight*slight);
        	double   kmmolm1 = 4*PI*PI*PI*NAvogadro/3/hPlank/slight/4/PI/eps0*1e-3*100.0*8.47835267e-30*8.47835267e-30;/* 1e-3 m to km, 100 : cm-1 to m-1 */
		double f = 1.0/sqrt(kmmolm1);
		int i;
		int nPOld = mol->vibration.nProperties;
		mol->vibration.nProperties += 3;
		for(i=0;i<mol->vibration.nModes;i++)
		{
			mol->vibration.modes[i].properties = realloc(mol->vibration.modes[i].properties, mol->vibration.nProperties*sizeof(double));
			mol->vibration.modes[i].properties[nPOld] = TX[i]*f;
			mol->vibration.modes[i].properties[nPOld+1] = TY[i]*f;
			mol->vibration.modes[i].properties[nPOld+2] = TZ[i]*f;
		}
	}
	free(TX);
	free(TY);
	free(TZ);
	fclose(file);
	return TRUE;
}
/********************************************************************************/
static boolean readDipolesDerivativesFromGabeditFile(Molecule* mol, char* fileName)
{
 	FILE *file;
	double** first = NULL;

	if(!mol) return FALSE;
	if(mol->vibration.nModes<1) return FALSE;
	//printf("nModes =%d\n",mol->vibration.nModes);
	//printf("fileName =%s\n",fileName);

 	file = fopen(fileName, "rb");
	if(!file) return FALSE;

	//printf("fileName =%s\n",fileName);
	first = newMatrixDouble(3,mol->vibration.nModes);
	if(!readMatrixReal(file,"First Derivatives",3, mol->vibration.nModes, first))
	{
		printf("I cannot read the dipole derivatives from the %s file\n",fileName);
		freeMatrixDouble(&first, 3);
		return FALSE;
	}
	else
	{
		int i;
		int j;
		int nPOld = mol->vibration.nProperties;
		mol->vibration.nProperties += 3;
		//printf("nprops = %d\n",mol->vibration.nProperties);
		for(i=0;i<mol->vibration.nModes;i++)
		{
			mol->vibration.modes[i].properties = realloc(mol->vibration.modes[i].properties, mol->vibration.nProperties*sizeof(double));
			for(j=0;j<3;j++) mol->vibration.modes[i].properties[nPOld+j] = first[j][i];
			for(j=0;j<3;j++) printf("%f ",mol->vibration.modes[i].properties[nPOld+j]);
			printf("\n");
		}
	}
	freeMatrixDouble(&first, 3);
	fclose(file);
	return TRUE;
}
/********************************************************************************/
static boolean readDipolesDerivativesFromOrcaOutputFile(Molecule* mol, char* fileName)
{
 	char t[BSIZE];
 	char sdum1[BSIZE];
 	char sdum2[BSIZE];
 	boolean ok;
	double freq;
 	FILE *file;
	int nf;
	int numberOfFrequencies = 0;
	double* TX = NULL;
	double* TY = NULL;
	double* TZ = NULL;
	//printf("nModes =%d\n",mol->vibration.nModes);
	//printf("fileName =%s\n",fileName);
	if(mol->vibration.nModes<1) return FALSE;
	//printf("fileName =%s\n",fileName);
	TX = malloc(mol->vibration.nModes*sizeof(double));
	TY = malloc(mol->vibration.nModes*sizeof(double));
	TZ = malloc(mol->vibration.nModes*sizeof(double));

 	file = fopen(fileName, "rb");
	if(!file) return FALSE;
	//printf("fileName =%s\n",fileName);

 	ok=FALSE;
	while(!feof(file))
	{
    		if(!fgets(t,BSIZE,file)) break;
 		if (strstr( t,"TX") && strstr( t,"TY") && strstr( t,"TZ") ) ok = TRUE;
 		if (strstr( t,":") && ok ){ ok = TRUE; break;}
	}
  	numberOfFrequencies = 0;
	if(ok)
  	while(!feof(file) )
  	{
		int i;
		if(!strstr(t,":")) break;
		for(i=0;i<strlen(t);i++) if(t[i]=='(' || t[i]==')') t[i]=' ';
		nf = sscanf(t,"%s %lf %s %lf %lf %lf",sdum1,&freq, sdum2, &TX[numberOfFrequencies], &TY[numberOfFrequencies], &TZ[numberOfFrequencies]);
		if(nf!=6) { ok = FALSE; break;}
		if(fabs(mol->vibration.modes[numberOfFrequencies].frequency-freq)>2.0)
		{
			fprintf(stderr,"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n");
			fprintf(stderr," WARNING : %f # %f \n", mol->vibration.modes[numberOfFrequencies].frequency,freq);
			fprintf(stderr,"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n");
		}
		numberOfFrequencies++;
		if(numberOfFrequencies>mol->vibration.nModes) { ok = FALSE; break;}
		if(!fgets(t,BSIZE,file)) break;
	}
	//printf("nModOld = %d nModnew = %d\n",mol->vibration.nModes,numberOfFrequencies);
	if(ok && numberOfFrequencies==mol->vibration.nModes)
	{
        	double mu0 = 4*PI*1e-7;
        	double eps0 = 1.0/(mu0*slight*slight);
        	double   kmmolm1 = 4*PI*PI*PI*NAvogadro/3/hPlank/slight/4/PI/eps0*1e-3*100.0*8.47835267e-30*8.47835267e-30;/* 1e-3 m to km, 100 : cm-1 to m-1 */
		double f = 1.0/sqrt(kmmolm1);
		int i;
		int nPOld = mol->vibration.nProperties;
		mol->vibration.nProperties += 3;
		for(i=0;i<mol->vibration.nModes;i++)
		{
			mol->vibration.modes[i].properties = realloc(mol->vibration.modes[i].properties, mol->vibration.nProperties*sizeof(double));
			mol->vibration.modes[i].properties[nPOld] = TX[i]*f;
			mol->vibration.modes[i].properties[nPOld+1] = TY[i]*f;
			mol->vibration.modes[i].properties[nPOld+2] = TZ[i]*f;
		}
	}
	free(TX);
	free(TY);
	free(TZ);
	fclose(file);
	return TRUE;
}
/********************************************************************************/
Molecule* readMoleculeFromOrcaOutputFile(char* fileName, int numgeometry)
{
	Molecule* mol = getMoleculeFromOrcaOutputFile(fileName, numgeometry);
	readVibrationFromOrcaOutputFile(mol,fileName);
	readDipolesDerivativesFromOrcaOutputFile(mol,fileName);
	return mol;
}
/********************************************************************************/
static boolean readVibrationFromGaussianOutputFile(Molecule* mol, char* fileName)
{
 	char t[BSIZE];
 	char sdum1[BSIZE];
 	char sdum2[BSIZE];
 	char sdum3[BSIZE];
 	boolean ok;
 	FILE *file;
	int idum,jdum,kdum;
	int nf;
	double freq[3] = {0,0,0};
	double IRIntensity[3] = {0,0,0};
	double mass[3] = {1,1,1};
	double RamanIntensity[3]={ 0,0,0};
	double v[3][3];
	char sym[3][BSIZE];
	int j;
	int k;
	int nProps;
	int nModes = 0;

	if(!mol) return FALSE;
        if(mol->nAtoms<1) return FALSE;
        if(!fileName) return FALSE;
        file = fopen(fileName,"rb");
        if(!file) { return FALSE; }

        nModes = 0;
	nProps = 2;
	initVibrations(mol, 3*mol->nAtoms, nProps);

	for(j=0;j<3;j++)
	{
		sprintf(sym[j]," ");
		for(k=0;k<3;k++) v[j][k] = 0;
	}

 	do 
 	{
 		ok=FALSE;
 		while(!feof(file))
		{
    			if(!fgets(t,BSIZE,file)) break;
	 		/* if ( strstr( t,"reduced masses") )*/
	 		if ( strstr( t,"and normal coordinates:") )
	  		{
				ok = TRUE;
				break;
	  		}
		}
  		while(!feof(file) )
  		{
    			if(!fgets(t,BSIZE,file)) break;
			if(isABackspace(t)) break;
			nf = sscanf(t,"%d %d %d",&idum,&jdum,&kdum);
			if(nf<=0 || nf>3)
			{
				break;
			}
			if(!fgets(t,BSIZE,file)) break;
			sscanf(t,"%s %s %s",sym[0],sym[1],sym[2]);

			if(!fgets(t,BSIZE,file)) break;
			changeDInE(t); 
			sscanf(t,"%s %s %lf %lf %lf", sdum1,sdum2, &freq[0],&freq[1],&freq[2]);
			while(!feof(file))
			{
				if(!fgets(t,BSIZE,file)) break;
				if(strstr(t,"Red."))
				{
					changeDInE(t); 
					sscanf(t,"%s %s %s %lf %lf %lf", sdum1,sdum2, sdum3, &mass[0],&mass[1],&mass[2]);
					break;
				}
			}
			while(!feof(file))
			{
				if(!fgets(t,BSIZE,file)) break;
				if(strstr(t,"IR Inten"))
				{
					changeDInE(t); 
					sscanf(t,"%s %s %s %lf %lf %lf", sdum1,sdum2, sdum3, &IRIntensity[0],&IRIntensity[1],&IRIntensity[2]);
					break;
				}
			}
			while(!feof(file))
			{
				if(!fgets(t,BSIZE,file)) break;
				if(strstr(t,"Raman"))
				{
					changeDInE(t); 
					sscanf(t,"%s %s %s %lf %lf %lf", sdum1,sdum2, sdum3, &RamanIntensity[0],&RamanIntensity[1],&RamanIntensity[2]);
					break;
				}
				if(strstr(t,"Atom ") && strstr(t," AN")) break;
			}

			if(!(strstr(t,"Atom ") && strstr(t," AN")))
			while(!feof(file))
			{
				if(!fgets(t,BSIZE,file)) break;
				if(strstr(t,"Atom ") && strstr(t," AN")) break;
			}
			for(k=0;k<nf;k++)
			{
				int n = nModes + k;
				if(n>3*mol->nAtoms) { ok = FALSE; break;}
				mol->vibration.modes[n].frequency = freq[k];
				mol->vibration.modes[n].properties[0] = IRIntensity[k];
				mol->vibration.modes[n].properties[1] = RamanIntensity[k];
				mol->vibration.modes[n].mass = mass[k];
			}

			for(j=0;j<mol->nAtoms && !feof(file);j++)
			{
				if(!fgets(t,BSIZE,file)) break;
				changeDInE(t); 
				sscanf(t,"%d %d %lf %lf %lf %lf %lf %lf %lf %lf %lf",
					&idum,&jdum,
					&v[0][0],&v[0][1],&v[0][2],
					&v[1][0],&v[1][1],&v[1][2],
					&v[2][0],&v[2][1],&v[2][2]
					);
				for(k=0;k<nf;k++)
				{
					int n = nModes + k;
					if(n>3*mol->nAtoms) { ok = FALSE; break;}
					mol->vibration.modes[n].vectors[0][j]= v[k][0]; 
					mol->vibration.modes[n].vectors[1][j]= v[k][1]; 
					mol->vibration.modes[n].vectors[2][j]= v[k][2]; 
				}
				if(!ok) break;
			}
			if(!ok || j!=mol->nAtoms) { ok = FALSE; freeVibrations(mol); break;}
			nModes += nf;
		}
 	}while(!feof(file));
	if(nModes<1)
	{
		char buffer[BSIZE];
		sprintf(buffer,"Sorry, I can not read frequencies from '%s' file\n",fileName);
  		fprintf(stderr,"%s\n",buffer);
	}
	else
	{
		mol->vibration.nModes = nModes;
		mol->vibration.modes = realloc(mol->vibration.modes, mol->vibration.nModes*sizeof(VibMode));
		sortFrequencies(mol);
		removeTransRotModes(mol);
	}
	return ok;
}
/********************************************************************************/
static Molecule* getMoleculeFromGaussianOutputFile(char *fileName, int numgeometry)
{
	char *t;
	boolean	OK;
	char *AtomCoord[5];
	FILE *file;
	int taille=BSIZE;
	int idummy;
	int i;
	int j=0;
	/* int l;*/
	int numgeom;
	char *pdest;
	int result;
	int itype=0;
	char* strStandard = "Standard orientation:";
	char* strInput = "Input orientation:";
	char* strOther = "orientation:";
	char* strSearch = strOther;
	Molecule* mol = newMolecule();
	char* symbol = NULL;

	for(i=0;i<5;i++) AtomCoord[i]=malloc(taille*sizeof(char));
	 
	t=malloc(taille*sizeof(char));
	file = fopen(fileName, "r");
	if(file ==NULL)
	{
		free(t);
	 	t = strdup_printf(("Sorry\nI can not open %s  file "),fileName);
	 	printf("%s\n",t);
	 	free(t);
	 	exit(1);
		return NULL;
	}
	while(!feof(file))
	{
		 if(!fgets(t,taille,file))break;
	         if(strstr( t, strStandard))
		 {
			 strSearch = strStandard;
			 break;
		 }
	         if(strstr( t, strInput)) strSearch = strInput;
	}
	fseek(file, 0L, SEEK_SET);
	numgeom =0;
	do 
	{
		OK=FALSE;
		while(!feof(file))
		{
		 	fgets(t,taille,file);
		 	if(strstr(t,"Charge =") && strstr(t,"Multiplicity ="))
		 	{
				 char* p = strstr(t,"Charge =")+8;
				mol->totalCharge = atoi(p);
			 	p = strstr(t,"Multiplicity =")+14;
				mol->spinMultiplicity = atoi(p);
		 	}
	         	pdest = strstr( t, strSearch);
	         	result = pdest - t ;
			if ( result >0 )
		 	{
		 		fgets(t,taille,file);
		 		fgets(t,taille,file);
		 		fgets(t,taille,file);
	               		pdest = strstr( t, "Type" );
	               		result = pdest - t ;
	               		if(result>0) itype=1;
	               		else itype=0;
		 		fgets(t,taille,file);
	               		numgeom++;
	               		if(numgeom == numgeometry )
				{
					OK = TRUE;
		 			break;
				}
				OK = TRUE;
				break;
			}
		}
		if(!OK && (numgeom == 0) )
		{
	 		free(t);
	 		t = strdup_printf(("Sorry\nI can not read geometry in  %s  file "),fileName);
	 		printf("%s\n",t);
	 		free(t);
	 		return NULL;
	   	}
	 	if(!OK)break;

	 	j=-1;
		if(mol->atoms) free(mol->atoms);
		mol->atoms = NULL;
		mol->nAtoms = 0;
	 	while(!feof(file) )
	 	{
	   		fgets(t,taille,file);
	   		pdest = strstr( t, "----------------------------------" );
	   		result = pdest - t ;
	   		if ( result >0 )
	   		{
				/*
				long geomposok = ftell(file);
	     			get_dipole_from_gaussian_output_file(file);
				fseek(file, geomposok, SEEK_SET);
				get_charges_from_gaussian_output_file(file,j+1);
				get_natural_charges_from_gaussian_output_file(file,j+1);
				fseek(file, geomposok, SEEK_SET);
				get_esp_charges_from_gaussian_output_file(file,j+1);
				*/
	     			break;
	   		}
	   		j++;

	   		if(itype==0) sscanf(t,"%d %s %s %s %s",&idummy,AtomCoord[0],AtomCoord[1],AtomCoord[2],AtomCoord[3]);
	   		else sscanf(t,"%d %s %d %s %s %s",&idummy,AtomCoord[0],&idummy,AtomCoord[1],AtomCoord[2],AtomCoord[3]);
			/* to do : test symbol */
			/*
			AtomCoord[0][0]=toupper(AtomCoord[0][0]);
			l=strlen(AtomCoord[0]);
	         	if (l==2) AtomCoord[0][1]=tolower(AtomCoord[0][1]);
			*/
			mol->atoms = realloc(mol->atoms,(j+1)*sizeof(Atom));

			symbol = getSymbolUsingZ(atoi(AtomCoord[0]));

			mol->atoms[j].prop = propAtomGet(symbol);
			mol->atoms[j].mmType=strdup(symbol);
			mol->atoms[j].pdbType=strdup(symbol);
			mol->atoms[j].residueName=strdup(symbol);
			mol->atoms[j].N=j+1;
			mol->atoms[j].layer=HIGH_LAYER;
			mol->atoms[j].variable=TRUE;
			mol->atoms[j].show=TRUE;
			mol->atoms[j].residueNumber=0;

			mol->atoms[j].coordinates[0]=atof(AtomCoord[1]);
			mol->atoms[j].coordinates[1]=atof(AtomCoord[2]);
			mol->atoms[j].coordinates[2]=atof(AtomCoord[3]);

			mol->atoms[j].gradient[0]=0;
			mol->atoms[j].gradient[1]=0;
			mol->atoms[j].gradient[2]=0;

			mol->atoms[j].velocity[0]=0.0;
			mol->atoms[j].velocity[1]=0.0;
			mol->atoms[j].velocity[2]=0.0;

	   		mol->atoms[j].charge=0.0;
	   		mol->atoms[j].charge0=0.0;
	   		mol->atoms[j].electronegativity=0.0;
	   		mol->atoms[j].hardness=0.0;
	   		mol->atoms[j].width= mol->atoms[j].prop.covalentRadii;
	   		mol->atoms[j].mass= mol->atoms[j].prop.mass;
	   		mol->atoms[j].rho=0.0;
	   		mol->atoms[j].U = 0.0;
		}
		mol->nAtoms = j+1;
		if(OK && numgeometry>-1) break;
	}while(!feof(file));

	fclose(file);
	free(t);
	for(i=0;i<5;i++) free(AtomCoord[i]);
        for(i=0; i<mol->nAtoms; i++) mol->atoms[i].typeConnections = NULL;
	resetTypeConnections(mol);
	int connections = 1;
	if(connections) setConnections(mol);
	else
	{
		createBondedMatrix(mol);
		printf(("Establishing connectivity : non bonded ...\n"));
		setNonBondedConnections(mol);
		freeBondedMatrix(mol);
	}

	/* if all freezed, set all to variable */
	{
		int j = 0;
		for(i=0;i<mol->nAtoms;i++)
			if(!mol->atoms[i].variable) j++;
		if(j==mol->nAtoms)
		for(i=0;i<mol->nAtoms;i++)
			mol->atoms[i].variable = TRUE;
	}
	printf("nAtoms = %d\n",mol->nAtoms);

	return mol;
}
/********************************************************************************/
Molecule* readMoleculeFromGaussianOutputFile(char *fileName, int numgeometry)
{
	Molecule* mol = getMoleculeFromGaussianOutputFile(fileName, numgeometry);
	readVibrationFromGaussianOutputFile(mol,fileName);
	readDipolesDerivativesFromGaussianOutputFile(mol, fileName);
	return mol;
}
/********************************************************************************/
static boolean readVibrationFromGaussianFChkFile(Molecule* mol, char* fileName)
{
 	FILE *file;
	int nf;
	double* vibE2 = NULL;
	double* vibNM = NULL;
	int idxFreq = 0;
	int idxMass = 0;
	int idxIR = 0;
	int idxRaman = 0;
	int i,j,k;
	int n;
	int nProps;
        int nModes = 0;

        if(!mol) return FALSE;
        if(mol->nAtoms<1) return FALSE;
        if(!fileName) return FALSE;
        file = fopen(fileName,"rb");
        if(!file) { return FALSE; }



	nf = get_one_int_from_fchk_gaussian_file(file,"Number of Normal Modes ");
	if(nf<1)
	{
		fprintf(stderr,"Sorry\nNo normal modes in this file : Use the Freq(SaveNM) option in your input file");
		fclose(file);
		return FALSE;
	}
	rewind(file);
	vibE2 = get_array_real_from_fchk_gaussian_file(file, "Vib-E2 ", &n);
	/* nf frequencies, nf Red. masses , nf Frc consts, nf IR Inten  , nf Raman Activ, nf Depolar (P), nf Depolar (U) */
	if(!vibE2 || n < 5*nf)
	{
		fprintf(stderr,"Sorry\nI can not the frequencies from this file");
		if(vibE2) free(vibE2);
		fclose(file);
		return FALSE;
	}
	rewind(file);
	vibNM = get_array_real_from_fchk_gaussian_file(file, "Vib-Modes ", &n);
	if(!vibNM || n != nf*mol->nAtoms*3)
	{
		fprintf(stderr,"Sorry\nI can not the normal modes from this file");
		printf("n = %d nf*nAtoms*3 = %d\n",n, nf*mol->nAtoms*3);
		if(vibE2) free(vibE2);
		if(vibNM) free(vibNM);
		fclose(file);
		return FALSE;
	}
	fclose(file);
	idxMass = nf;
	idxIR = 3*nf;
	idxRaman = 4*nf;

        nModes = nf;
        nProps = 2;
        initVibrations(mol, nModes, nProps);

	for(k=0;k<nf;k++)
	{
		mol->vibration.modes[k].frequency = vibE2[idxFreq+k];
		mol->vibration.modes[k].properties[0] = vibE2[idxIR+k];
		mol->vibration.modes[k].mass = vibE2[idxMass+k]; 
		mol->vibration.modes[k].properties[1] = vibE2[idxRaman+k];
		for(i=0;i<3;i++)
		{
			for(j=0;j<mol->nAtoms;j++)
			{
				mol->vibration.modes[k].vectors[i][j] = vibNM[k*mol->nAtoms*3+3*j+i];
			}
		}
	}
	if(vibE2) free(vibE2);
	if(vibNM) free(vibNM);
	sortFrequencies(mol);
	removeTransRotModes(mol);

	return TRUE;
}
/********************************************************************************/
static Molecule* getMoleculeFromGaussianFChkFile(char *fileName)
{
 	FILE *file;
	int i,j;
	int n;
	double* coords = NULL;
	double* charges = NULL;
	double* dipole = NULL;
	int* z = NULL;
	double* zn = NULL;
	Molecule* mol;

	file = fopen(fileName, "rb");
	if(file ==NULL)
	{
  		fprintf(stderr,"Sorry\nI can not open this file");
  		return FALSE;
	}

	j = get_one_int_from_fchk_gaussian_file(file,"Number of atoms ");
	if(j<1)
	{
  		fprintf(stderr,"Sorry\nI can not the number of atoms from this file");
  		return FALSE;
	}
	z = get_array_int_from_fchk_gaussian_file(file, "Atomic numbers ", &n);
	if(n!=j)
	{
  		fprintf(stderr,"Sorry\nI can not read the atomic numbers from this file");
  		return FALSE;
	}
	coords = get_array_real_from_fchk_gaussian_file(file, "Current cartesian coordinates  ", &n);
	if(n!=3*j)
	{
  		fprintf(stderr,"Sorry\nI can not read the current cartesian coordinates from this file");
  		return FALSE;
	}
	rewind(file);
	mol = newMolecule();
	mol->nAtoms = j;

    	mol->atoms = malloc(mol->nAtoms*sizeof(Atom));
	for(j=0;j<mol->nAtoms;j++)
	{
    		char* t = strdup(symbAtomGet(z[j]));
               	mol->atoms[j].prop = propAtomGet(t);
                mol->atoms[j].mmType=strdup(t);
                mol->atoms[j].pdbType=strdup(t);
                mol->atoms[j].residueName=strdup(t);
                mol->atoms[j].N=j+1;
                mol->atoms[j].layer=HIGH_LAYER;
                mol->atoms[j].variable=TRUE;
                mol->atoms[j].show=TRUE;
		mol->atoms[j].residueNumber=0;
                mol->atoms[j].charge=0.0;
                mol->atoms[j].charge0=0.0;
                mol->atoms[j].electronegativity=0.0;
                mol->atoms[j].hardness=0.0;
	   	mol->atoms[j].width=mol->atoms[j].prop.covalentRadii;
                mol->atoms[j].mass=mol->atoms[j].prop.mass;
                mol->atoms[j].rho=0.0;
                mol->atoms[j].U = 0.0;

		for(i=0;i<3;i++) mol->atoms[j].coordinates[i]= coords[j*3+i]; 
		for(i=0;i<3;i++) mol->atoms[j].gradient[i] = 0.0;
		for(i=0;i<3;i++) mol->atoms[j].velocity[j]=0.0;
		free(t);
	}

	int connections = 1;
        for(i=0; i<mol->nAtoms; i++) mol->atoms[i].typeConnections = NULL;
	resetTypeConnections(mol);
	if(connections) setConnections(mol);

	if(z) free(z);
	if(coords) free(coords);
	z = NULL;
	coords = NULL;
	rewind(file);
	charges = get_array_real_from_fchk_gaussian_file(file, "NPA Charges ", &n);
	if(n==mol->nAtoms && charges) { for(j=0;j<mol->nAtoms;j++) mol->atoms[j].charge =charges[j]; }
	else
	{
		rewind(file);
		charges = get_array_real_from_fchk_gaussian_file(file, "ESP Charges  ", &n);
		if(n==mol->nAtoms && charges) { for(j=0;j<mol->nAtoms;j++) mol->atoms[j].charge =charges[j]; }
		else
		{
			rewind(file);
			charges = get_array_real_from_fchk_gaussian_file(file, "Mulliken Charges  ", &n);
			if(n==mol->nAtoms && charges) { for(j=0;j<mol->nAtoms;j++) mol->atoms[j].charge =charges[j]; }
		}
	}
	if(charges) free(charges);
	rewind(file);
	zn = get_array_real_from_fchk_gaussian_file(file, "Nuclear charges ", &n);
	if(zn && n== j)
	{
		for(j=0;j<mol->nAtoms;j++) mol->atoms[j].charge0 =zn[j];
	}
	if(zn)free(zn);
	if(n!=j)
	{
  		fprintf(stderr,"Sorry\nI can not read the atomic numbers from this file");
  		return FALSE;
	}
	dipole = get_array_real_from_fchk_gaussian_file(file, "Dipole Moment  ", &n);
	for(i=0;i<3;i++) mol->dipole[i] = 0.0;
	if(n==3)
	{
		for(i=0;i<3;i++) mol->dipole[i] = dipole[i];
	}
	if(dipole) free(dipole);
 	fclose(file);
	return mol;
}
/********************************************************************************/
Molecule* readMoleculeFromGaussianFChkFile(char *fileName)
{
	Molecule* mol = getMoleculeFromGaussianFChkFile(fileName);
	readVibrationFromGaussianFChkFile(mol,fileName);
	return mol;
}
/********************************************************************************/
static boolean readVibrationFromGabeditFile(Molecule* mol, char* fileName)
{
	FILE* file = NULL;
	static char *t = NULL; 
	int i,j,k;
	int nModes = 0;
	int nProps = 0;

	if(!mol) return FALSE;
	if(!fileName) return FALSE; 
	if(t==NULL) t = malloc(BSIZE*sizeof(char));
	file = fopen(fileName,"rb");
	if(!file) { return FALSE; }
	/* compute number of modes */
	nModes = 0;
	if(goToStr(file, "[FREQ]"))
	{ 
		double d;
		while(!feof(file))
		{
    			if(!fgets(t,BSIZE, file)) break;
			deleteFirstSpaces(t);
			if(t[0]=='#') continue;
			if(sscanf(t,"%lf",&d)!=1) break;
			nModes++;
		}
	}
	if(nModes<1) { fclose(file); return FALSE;}
	/* compute number of props */
	nProps = 0;
	if(goToStr(file, "[INT]"))
	{ 
		double d1,d2,d3,d4;
		while(!feof(file))
		{
    			if(!fgets(t,BSIZE, file)) break;
			deleteFirstSpaces(t);
			if(t[0]=='#') continue;
			nProps = sscanf(t,"%lf %lf %lf %lf",&d1,&d2,&d3,&d4);
			break;
		}
	}
	initVibrations(mol, nModes, nProps);
	/* read harmonic frequencies */
	if(goToStr(file, "[FREQ]"))
	{ 
		double d;
		for(i=0;i<mol->vibration.nModes;i++)
		{
    			if(!fgets(t,BSIZE, file)) break;
			deleteFirstSpaces(t);
			if(t[0]=='#') { i--; continue;}
			if(sscanf(t,"%lf",&d)!=1) break;
			mol->vibration.modes[i].frequency = d;
		}
	}
	/* read masses */
	if(goToStr(file, "[MASS]"))
	{ 
		double d;
		for(i=0;i<mol->vibration.nModes;i++)
		{
    			if(!fgets(t,BSIZE, file)) break;
			deleteFirstSpaces(t);
			if(t[0]=='#') { i--; continue;}
			if(sscanf(t,"%lf",&d)!=1) break;
			mol->vibration.modes[i].mass = d;
		}
		if(i!=mol->vibration.nModes)
		{
			fprintf(stderr,"----------------- Error -------------------------------------------------\n");
			fprintf(stderr,"Error , I cannot read masses from %s file, check your file\n",fileName);
			fprintf(stderr,"-------------------------------------------------------------------------\n");
			exit(1);
		}
	}
	else
	{
		fprintf(stderr,"----------------- Warning -------------------------------------------------\n");
		fprintf(stderr,"Error , I cannot read masses from %s file, check your file\n",fileName);
		fprintf(stderr,"-------------------------------------------------------------------------\n");
	}
	/* read intensities */
	if(goToStr(file, "[INT]"))
	{ 
		double d[4];
		int n;
		for(i=0;i<mol->vibration.nModes;i++)
		{
    			if(!fgets(t,BSIZE, file)) break;
			deleteFirstSpaces(t);
			if(t[0]=='#') { i--; continue;}
			n = sscanf(t,"%lf %lf %lf %lf",&d[0],&d[1],&d[2],&d[3]);
			if(n<<mol->vibration.nProperties)
			for(j=0;j<n;j++) mol->vibration.modes[i].properties[j] = d[j];
		}
	}
	/* read vectors */
	if(goToStr(file, "[FR-NORM-COORD]"))
	{ 
		double d[3];
		for(i=0;i<mol->vibration.nModes;i++)
		{
    			if(!fgets(t,BSIZE, file)) break; /* vibration i title */
			deleteFirstSpaces(t);
			if(t[0]=='#') { i--; continue;}
			for(j=0;j<mol->nAtoms;j++) 
			{
    				if(!fgets(t,BSIZE, file)) break;
				deleteFirstSpaces(t);
				if(t[0]=='#') {j--; continue;}
				sscanf(t,"%lf %lf %lf",&d[0],&d[1],&d[2]);
				for(k=0;k<3;k++) mol->vibration.modes[i].vectors[k][j] = d[k];
			}
		}
	}
	fclose(file);
	sortFrequencies(mol);
	removeTransRotModes(mol);
	return TRUE;
}
/********************************************************************************/
static boolean readGradientFromGabeditFile(Molecule* mol, char* namefile)
{
	FILE* file = NULL;
	static char *t = NULL; 
	boolean Ok = FALSE;
	int n,ic,is;
#define SZ 50
	int i;
	char* pos;
	int nGeoms = 0;
	int nLabels = 0;

	if(!namefile) 
	{
		printf("Sorry I cannot read geometry namefile = NULL\n");
		exit(1);
	}
	if(t==NULL) t = malloc(BSIZE*sizeof(char));
	file = fopen(namefile,"rb");
	if(!file)
	{
		printf("Sorry I cannot open %s file\n",namefile);
		exit(1);
	}
	rewind(file);
	while(!feof(file))
  	{
    		if(!fgets(t,BSIZE, file)) break;
		deleteFirstSpaces(t);
		if(t[0]=='#') continue;
		uppercase(t);
		pos = strstr(t,"[GEOMS]");
		if(pos)
		{ 
			while(!feof(file))
			{
    				if(!fgets(t,BSIZE, file)) break;
				deleteFirstSpaces(t);
				if(t[0]=='#') continue;
				break;
			}
			if(2==sscanf(t,"%d%d",&nGeoms,&nLabels))
			{

				for(i=0; i<nLabels; i++)
    					if(!fgets(t,BSIZE, file)) break;
				if(i==nLabels) Ok = TRUE;
			}
			break;
		}
	}
	/* I read the gradient of the first geometry */
	for(i=0; i<nLabels; i++)
    		if(!fgets(t,BSIZE, file)) break;
	if(i!=nLabels) Ok = FALSE;
    	if(!fgets(t,BSIZE, file)) {Ok = FALSE; }
	if(Ok && 3==sscanf(t,"%d%d%d",&n,&ic,&is) && n>0 && is>0)
	{
		/*printf("Mult = %d\n",is);*/
		if(mol->nAtoms!= n || mol->spinMultiplicity != is || mol->totalCharge != ic)
		{
			printf("Sorry geometry in file %s and that in memory have not the same number of atoms or spin or charge\n",namefile);
			exit(1);
			Ok = FALSE;
		}
	}
	else Ok = FALSE;
	if(!Ok)
	{
		printf("Sorry I cannot read geometry from %s file\n",namefile);
		exit(1);
	}
	for(i=0; i<mol->nAtoms; i++)
	{
			if(!fgets(t,BSIZE,file))
			{
				printf("Sorry I cannot read geometry from %s file.\n",namefile);
				exit(1);
			}
			deleteFirstSpaces(t);
			if(t[0]=='#') { i--;continue;}
			get_gradients_one_atom(t, mol->atoms[i].gradient);
	}
	
	fclose(file);
	return Ok;
}
/********************************************************************************/
static Molecule** readMoleculesFromCChemIFile(char *fileName, boolean connections)
{
	char tmp[BSIZE];
	FILE *file;
	int nAtoms;
	Molecule** mols = NULL;
	int tcharge=0;
	int mult=1;
#define SZ 50
	char symbol[SZ];
	char mmType[SZ];
	char pdbType[SZ];
	char residueName[SZ];
	double X,Y,Z;
	double charge;
	int layer;
	double energy;

	file = fopen(fileName, "r");
	if(file ==NULL)
	{
	 	sprintf(tmp,"Sorry\nI can not open %s  file ",fileName);
	 	fprintf(stderr,"%s\n",tmp);
	 	exit(1);
		return NULL;
	}
	while(!feof(file))
  	{
    		if(!fgets(tmp,BSIZE, file)) break;
		deleteFirstSpaces(tmp);
		if(tmp[0]=='#') continue;
		uppercase(tmp);
		if(strstr(tmp,"GEOMETRY"))break;
	}
	if(!strstr(tmp,"GEOMETRY"))
	{
	 	sprintf(tmp,"Sorry\nI cannot read geometries from %s file ",fileName);
	 	fprintf(stderr,"%s\n",tmp);
	 	exit(1);
		return NULL;
	}
	int nGeoms = 0;
	while(!feof(file))
	{
		if(!fgets(tmp,BSIZE,file))break;
		energy=0;
		if(3>sscanf(tmp,"%d%d%d %lf",&nAtoms,&tcharge,&mult,&energy)) break;
		Molecule* mol = newMolecule();
		mol->totalCharge = tcharge;
		mol->spinMultiplicity = mult;
		mol->potentialEnergy = energy;
		mol->atoms = NULL;
		mol->nAtoms = nAtoms;
		mol->atoms = malloc(nAtoms*sizeof(Atom));
		int i;
		for(i=0; i<mol->nAtoms; i++) mol->atoms[i].typeConnections = malloc(mol->nAtoms*sizeof(int));
		int l;
		for(i=0; i<mol->nAtoms; i++)
		{
			int variable = 0;
			int ibeg = 12;
			if(!fgets(tmp,BSIZE,file))
			{
				fprintf(stderr,"Sorry I cannot read geometry from %s file.\n",fileName);
				exit(1);
			}
			deleteFirstSpaces(tmp);
			if(tmp[0]=='#') { i--;continue;}
    			sscanf(tmp,"%s %s %s %s %d %lf %d %d %lf %lf %lf",
					symbol,mmType,pdbType,residueName, 
					&mol->atoms[i].residueNumber,
					&charge,&layer,&variable,&X,&Y,&Z);
			symbol[0]=toupper(symbol[0]);
			l=strlen(symbol);
			if (l==2) symbol[1]=tolower(symbol[1]);

			mol->atoms[i].prop = propAtomGet(symbol);
			mol->atoms[i].mmType=strdup(mmType);
			mol->atoms[i].pdbType=strdup(pdbType);
			mol->atoms[i].residueName=strdup(residueName);
			mol->atoms[i].N=i+1;
			mol->atoms[i].layer=layer;
			mol->atoms[i].variable=variable;
			mol->atoms[i].show=TRUE;
			mol->atoms[i].coordinates[0] = X;
			mol->atoms[i].coordinates[1] = Y;
			mol->atoms[i].coordinates[2] = Z;
			mol->atoms[i].velocity[0] = 0;
			mol->atoms[i].velocity[1] = 0;
			mol->atoms[i].velocity[2] = 0;
			mol->atoms[i].gradient[0]=0;
			mol->atoms[i].gradient[1]=0;
			mol->atoms[i].gradient[2]=0;
			mol->atoms[i].charge = charge;
			mol->atoms[i].charge0 = charge;
			mol->atoms[i].electronegativity = 0;
			mol->atoms[i].hardness = 0;
			mol->atoms[i].width = mol->atoms[i].prop.covalentRadii;
			mol->atoms[i].mass = mol->atoms[i].prop.mass;
			mol->atoms[i].rho = 0.0;
	   		mol->atoms[i].U = 0.0;
			if(!get_connections_one_atom(tmp, mol->nAtoms, ibeg, mol->atoms[i].typeConnections))
			{
				/*
				fprintf(stderr,"Sorry I cannot read the connection for atom # %d from the %s file.\n",i+1,fileName);
				exit(1);
				*/
				fprintf(stderr,"Warning : I cannot read the connection for atom # %d from the %s file.\n",i+1,fileName);
			}
		}
		if(i==nAtoms) 
		{
			if(connections) setConnections(mol);
			else {
				createBondedMatrix(mol);
				fprintf(stderr,"Establishing connectivity : non bonded ...\n");
				setNonBondedConnections(mol);
				freeBondedMatrix(mol);
			}
			/* if all freezed, set all to variable */
			{
				int jv = 0;
				int i;
				for(i=0;i<mol->nAtoms;i++)
					if(!mol->atoms[i].variable) jv++;
				if(jv==mol->nAtoms) for(i=0;i<mol->nAtoms;i++) mol->atoms[i].variable = TRUE;
			}
			nGeoms++;
			mols = realloc(mols, nGeoms*sizeof(Molecule*));
			mols[nGeoms-1] = mol;
		}
	}
	if(nGeoms>0){
			mols = realloc(mols, (nGeoms+1)*sizeof(Molecule*));
			mols[nGeoms] = NULL;
		}
	fclose(file);
	return mols;
}
/********************************************************************************/
static Molecule** readMoleculesFromXYZFile(char *fileName, boolean connections)
{
	char tmp[BSIZE];
	FILE *file;
	char symbol[BSIZE];
	double x,y,z;
	int nAtoms;
	Molecule** mols = NULL;
	int charge=0;
	int mult=1;
	double energy;

	file = fopen(fileName, "r");
	if(file ==NULL)
	{
	 	sprintf(tmp,"Sorry\nI can not open %s  file ",fileName);
	 	printf("%s\n",tmp);
	 	exit(1);
		return NULL;
	}
	int nGeoms = 0;
	while(!feof(file))
	{
		if(!fgets(tmp,BSIZE,file))break;
		if(1!=sscanf(tmp,"%d",&nAtoms)) break;
		if(!fgets(tmp,BSIZE,file))break;
		int a,b;
		if(2>sscanf(tmp,"%d %d %lf",&a,&b,&energy)) { charge=a; mult=b;} 
		Molecule* mol = newMolecule();
		mol->totalCharge = charge;
		mol->spinMultiplicity = mult;
		mol->potentialEnergy = energy;
		mol->atoms = NULL;
		mol->nAtoms = nAtoms;
		mol->atoms = malloc(nAtoms*sizeof(Atom));
		int j;
		for(j=0;j<nAtoms;j++)
		{
			if(!fgets(tmp,BSIZE,file))break;
			if(4==sscanf(tmp,"%s %lf %lf %lf",symbol, &x, &y, &z)) { 
			mol->atoms[j].prop = propAtomGet(symbol);
			mol->atoms[j].mmType=strdup(symbol);
			mol->atoms[j].pdbType=strdup(symbol);
			mol->atoms[j].residueName=strdup(symbol);
			mol->atoms[j].N=j+1;
			mol->atoms[j].layer=HIGH_LAYER;
			mol->atoms[j].variable=TRUE;
			mol->atoms[j].show=TRUE;
			mol->atoms[j].residueNumber=0;

			mol->atoms[j].coordinates[0]=x;
			mol->atoms[j].coordinates[1]=y;
			mol->atoms[j].coordinates[2]=y;

			mol->atoms[j].gradient[0]=0;
			mol->atoms[j].gradient[1]=0;
			mol->atoms[j].gradient[2]=0;

			mol->atoms[j].velocity[0]=0.0;
			mol->atoms[j].velocity[1]=0.0;
			mol->atoms[j].velocity[2]=0.0;

	   		mol->atoms[j].charge=0.0;
	   		mol->atoms[j].charge0=0.0;
	   		mol->atoms[j].electronegativity=0.0;
	   		mol->atoms[j].hardness=0.0;
	   		mol->atoms[j].width= mol->atoms[j].prop.covalentRadii;
	   		mol->atoms[j].mass= mol->atoms[j].prop.mass;
	   		mol->atoms[j].rho=0.0;
	   		mol->atoms[j].U = 0.0;
			} 
		}
		if(j==nAtoms) 
		{
			int i;
        		for(i=0; i<mol->nAtoms; i++) mol->atoms[i].typeConnections = NULL;
			resetTypeConnections(mol);
			if(connections) setConnections(mol);
			else {
				createBondedMatrix(mol);
				printf(("Establishing connectivity : non bonded ...\n"));
				setNonBondedConnections(mol);
				freeBondedMatrix(mol);
			}
			/* if all freezed, set all to variable */
			{
				int jv = 0;
				int i;
				for(i=0;i<mol->nAtoms;i++)
					if(!mol->atoms[i].variable) jv++;
				if(jv==mol->nAtoms) for(i=0;i<mol->nAtoms;i++) mol->atoms[i].variable = TRUE;
			}
			nGeoms++;
			mols = realloc(mols, nGeoms*sizeof(Molecule*));
			mols[nGeoms-1] = mol;
		}
	}
	if(nGeoms>0) {
			mols = realloc(mols, (nGeoms+1)*sizeof(Molecule*));
			mols[nGeoms] = NULL;
		}
	fclose(file);
	return mols;
}
/********************************************************************************/
static Molecule** readMoleculesFromGabeditFile(char* namefile, boolean connections)
{
	FILE* file = NULL;
	static char *t = NULL; 
	boolean Ok = FALSE;
	int n,ic,is;
	Molecule** mols = NULL;
#define SZ 50
	char symbol[SZ];
	char mmType[SZ];
	char pdbType[SZ];
	char residueName[SZ];
	double X,Y,Z;
	double charge;
	int layer;
	char* pos;
	int nGeoms = 0;
	int nLabels = 0;
	int iEnergy;
	double convEnergy = 1.0;

	if(!namefile) 
	{
		fprintf(stderr,"Sorry I cannot read geometry namefile = NULL\n");
		exit(1);
	}
	if(t==NULL) t = malloc(BSIZE*sizeof(char));
	file = fopen(namefile,"rb");
	if(!file)
	{
		fprintf(stderr,"Sorry I cannot open %s file\n",namefile);
		exit(1);
	}
	rewind(file);
	iEnergy=-1;
	while(!feof(file))
  	{
    		if(!fgets(t,BSIZE, file)) break;
		deleteFirstSpaces(t);
		if(t[0]=='#') continue;
		uppercase(t);
		pos = strstr(t,"[GEOMS]");
		if(pos)
		{ 
			while(!feof(file))
			{
    				if(!fgets(t,BSIZE, file)) break;
				deleteFirstSpaces(t);
				if(t[0]=='#') continue;
				break;
			}
			if(2==sscanf(t,"%d%d",&nGeoms,&nLabels))
			{

				int i;
				for(i=0; i<nLabels; i++)
				{
    					if(!fgets(t,BSIZE, file)) break;
					if(strstr(t,"nergy")) 
					{
						iEnergy=i;
						if(strstr(t,"tree")) convEnergy = AUTOKCAL;
						if(strstr(t,"eV")) convEnergy = AUTOKCAL/(AUTOEV);
					}
				}
				if(i==nLabels) Ok = TRUE;
			}
			break;
		}
	}
	mols = (Molecule**) malloc((nGeoms+1)*sizeof(Molecule*));
	int ig;
	for(ig=0; ig<nGeoms; ig++) mols[ig] = newMolecule();
	mols[nGeoms] = NULL;

	for(ig=0; ig<nGeoms; ig++) {
		int ie = 0;
		double energy=0;
		double e;
		Molecule* mol = mols[ig];
		/* I read a geometry */
		int i;
		for(i=0; i<nLabels; i++) 
		{
			if(!fgets(t,BSIZE, file)) break;
			if(iEnergy == i) energy = atof(t)*convEnergy;
		}
		if(i!=nLabels) Ok = FALSE;
    		if(!fgets(t,BSIZE, file)) {Ok = FALSE; }
		if(Ok) ie = sscanf(t,"%d%d%d%lf",&n,&ic,&is,&e);
		if(Ok && 3<=ie && n>0 && is>0)
		{
			/*printf("Mult = %d\n",is);*/
			mol->nAtoms = n;
			mol->spinMultiplicity = is;
			mol->totalCharge = ic;

			if(ie==4) mol->potentialEnergy = e;
			else mol->potentialEnergy = energy;

			mol->atoms = malloc(mol->nAtoms*sizeof(Atom));
			for(i=0; i<mol->nAtoms; i++) mol->atoms[i].typeConnections = malloc(mol->nAtoms*sizeof(int));
		}
		else Ok = FALSE;
		if(!Ok)
		{
			fprintf(stderr,"Sorry I cannot read geometry from %s file\n",namefile);
			exit(1);
		}
		for(i=0; i<mol->nAtoms; i++)
		{
			int variable = 0;
			int ibeg = 12;
			int l;
			if(!fgets(t,BSIZE,file))
			{
				fprintf(stderr,"Sorry I cannot read geometry from %s file.\n",namefile);
				exit(1);
			}
			deleteFirstSpaces(t);
			if(t[0]=='#') { i--;continue;}
    			sscanf(t,"%s %s %s %s %d %lf %d %d %lf %lf %lf", symbol,mmType,pdbType,residueName, 
					&mol->atoms[i].residueNumber,
					&charge,&layer,&variable,&X,&Y,&Z);
			symbol[0]=toupper(symbol[0]);
			l=strlen(symbol);
			if (l==2) symbol[1]=tolower(symbol[1]);

			mol->atoms[i].prop = propAtomGet(symbol);
			mol->atoms[i].mmType=strdup(mmType);
			mol->atoms[i].pdbType=strdup(pdbType);
			mol->atoms[i].residueName=strdup(residueName);
			mol->atoms[i].N=i+1;
			mol->atoms[i].layer=layer;
			mol->atoms[i].variable=variable;
			mol->atoms[i].show=TRUE;
			mol->atoms[i].coordinates[0] = X;
			mol->atoms[i].coordinates[1] = Y;
			mol->atoms[i].coordinates[2] = Z;
			mol->atoms[i].gradient[0]=0;
			mol->atoms[i].gradient[1]=0;
			mol->atoms[i].gradient[2]=0;
			mol->atoms[i].velocity[0]=0.0;
			mol->atoms[i].velocity[1]=0.0;
			mol->atoms[i].velocity[2]=0.0;
			mol->atoms[i].charge = charge;
			mol->atoms[i].charge0 = charge;
			mol->atoms[i].electronegativity = 0;
			mol->atoms[i].hardness = 0;
			mol->atoms[i].width = mol->atoms[i].prop.covalentRadii;
			mol->atoms[i].mass = mol->atoms[i].prop.mass;
			mol->atoms[i].rho = 0.0;
	   		mol->atoms[i].U = 0.0;
			if(!get_connections_one_atom(t, mol->nAtoms, ibeg, mol->atoms[i].typeConnections))
			{
				/*
				printf("Sorry I cannot read the connection for atom # %d from the %s file.\n",i+1,namefile);
				exit(1);
				*/
				fprintf(stderr,"Warning : I cannot read the connection for atom # %d from the %s file.\n",i+1,namefile);
			}
			get_gradients_one_atom(t, mol->atoms[i].gradient);
		}
		if(connections) setConnections(mol);
		else {
			createBondedMatrix(mol);
			printf(("Establishing connectivity : non bonded ...\n"));
			setNonBondedConnections(mol);
			freeBondedMatrix(mol);
		}
		/* if all freezed, set all to variable */
		{
			int jv = 0;
			for(i=0;i<mol->nAtoms;i++) if(!mol->atoms[i].variable) jv++;
			if(jv==mol->nAtoms) for(i=0;i<mol->nAtoms;i++) mol->atoms[i].variable = TRUE;
		}
	}
	fclose(file);
	/* for(ig=0; ig<nGeoms; ig++) fprintf(stdout,"Mol %d Energy = %f\n",ig,mols[ig]->potentialEnergy);*/
	return mols;
}
/*************************************************************************************/
Molecule** readMolecules(char* fileName, boolean connections)
{
	char* fileNameToRead = malloc(BSIZE*sizeof(char));
	CCHEMIFileType type;
	Molecule** mols = NULL;

	if(!fileName) return mols;

	fileNameToRead = malloc(BSIZE*sizeof(char));
	getFileNameToRead(fileName, fileNameToRead);

	type = getTypeFile(fileNameToRead);

	//printf("type %d\n",type);

	if(type == FILETYPE_CCHEMI) mols = readMoleculesFromCChemIFile(fileNameToRead, connections);
	else if(type == FILETYPE_GABEDIT) mols = readMoleculesFromGabeditFile(fileNameToRead, connections);
	else if(type == FILETYPE_XYZ) mols = readMoleculesFromXYZFile(fileNameToRead, connections);
	else
	{
		fprintf(stderr,"Sorry, I cannot read severla geometries from the molecule from %s file\n",fileNameToRead);
		fprintf(stderr,"       Only Gabedit & XYZ are supported for reading severla geometries\n");
		exit(1);
	}

	free(fileNameToRead);
	return mols;
}
/********************************************************************************/
static Molecule* getMoleculeFromGabeditFile_GEOMS(char* namefile)
{
	FILE* file = NULL;
	static char *t = NULL; 
	boolean Ok = FALSE;
	int n,ic,is;
	Molecule* mol = newMolecule();
#define SZ 50
	char symbol[SZ];
	char mmType[SZ];
	char pdbType[SZ];
	char residueName[SZ];
	double X,Y,Z;
	double charge;
	int layer;
	int i,l;
	char* pos;
	int nGeoms = 0;
	int nLabels = 0;
	int iEnergy = -1;
	double convEnergy = 1.0;
	int iDipole = -1;
	double convDipole = 1.0;

	if(!namefile) 
	{
		printf("Sorry I cannot read geometry namefile = NULL\n");
		exit(1);
	}
	if(t==NULL) t = malloc(BSIZE*sizeof(char));
	file = fopen(namefile,"rb");
	if(!file)
	{
		printf("Sorry I cannot open %s file\n",namefile);
		exit(1);
	}
	rewind(file);
	while(!feof(file))
  	{
    		if(!fgets(t,BSIZE, file)) break;
		deleteFirstSpaces(t);
		if(t[0]=='#') continue;
		uppercase(t);
		pos = strstr(t,"[GEOMS]");
		if(pos)
		{ 
			while(!feof(file))
			{
    				if(!fgets(t,BSIZE, file)) break;
				deleteFirstSpaces(t);
				if(t[0]=='#') continue;
				break;
			}
			if(2==sscanf(t,"%d%d",&nGeoms,&nLabels))
			{

				for(i=0; i<nLabels; i++)
				{
    					if(!fgets(t,BSIZE, file)) break;
					if(strstr(t,"nergy")) 
					{
						iEnergy=i;
						if(strstr(t,"tree")) convEnergy = AUTOKCAL;
						if(strstr(t,"eV")) convEnergy = AUTOKCAL/(AUTOEV);
					}
					if(strstr(t,"pole")) 
					{
						iDipole=i;
						if(strstr(t,"au") || strstr(t,"atomi")) convDipole = AUTODEB;
					}
				}
				if(i==nLabels) Ok = TRUE;
			}
			break;
		}
	}
	/* I read the first geometry */
	for(i=0; i<nLabels; i++) 
	{
		if(!fgets(t,BSIZE, file)) break;
		if(iEnergy == i) mol->potentialEnergy = atof(t)*convEnergy;
		if(iDipole == i) 
		{
			int k;
			sscanf(t,"%lf %lf %lf",&mol->dipole[0], &mol->dipole[1], &mol->dipole[2]);
			for(k=0;k<3;k++) mol->dipole[k] *= convDipole;
		}
	}
	if(i!=nLabels) Ok = FALSE;
    	if(!fgets(t,BSIZE, file)) {Ok = FALSE; }
	if(Ok && 3==sscanf(t,"%d%d%d",&n,&ic,&is) && n>0 && is>0)
	{
		/*printf("Mult = %d\n",is);*/
		mol->nAtoms = n;
		mol->spinMultiplicity = is;
		mol->totalCharge = ic;
		mol->atoms = malloc(mol->nAtoms*sizeof(Atom));
		for(i=0; i<mol->nAtoms; i++) mol->atoms[i].typeConnections = malloc(mol->nAtoms*sizeof(int));
	}
	else Ok = FALSE;
	if(!Ok)
	{
		printf("Sorry I cannot read geometry from %s file\n",namefile);
		exit(1);
	}
	for(i=0; i<mol->nAtoms; i++)
	{
			int variable = 0;
			int ibeg = 12;
			if(!fgets(t,BSIZE,file))
			{
				printf("Sorry I cannot read geometry from %s file.\n",namefile);
				exit(1);
			}
			deleteFirstSpaces(t);
			if(t[0]=='#') { i--;continue;}
    			sscanf(t,"%s %s %s %s %d %lf %d %d %lf %lf %lf", symbol,mmType,pdbType,residueName, 
					&mol->atoms[i].residueNumber,
					&charge,&layer,&variable,&X,&Y,&Z);
			symbol[0]=toupper(symbol[0]);
			l=strlen(symbol);
			if (l==2) symbol[1]=tolower(symbol[1]);

			mol->atoms[i].prop = propAtomGet(symbol);
			mol->atoms[i].mmType=strdup(mmType);
			mol->atoms[i].pdbType=strdup(pdbType);
			mol->atoms[i].residueName=strdup(residueName);
			mol->atoms[i].N=i+1;
			mol->atoms[i].layer=layer;
			mol->atoms[i].variable=variable;
			mol->atoms[i].show=TRUE;
			mol->atoms[i].coordinates[0] = X;
			mol->atoms[i].coordinates[1] = Y;
			mol->atoms[i].coordinates[2] = Z;
			mol->atoms[i].gradient[0]=0;
			mol->atoms[i].gradient[1]=0;
			mol->atoms[i].gradient[2]=0;
			mol->atoms[i].velocity[0]=0.0;
			mol->atoms[i].velocity[1]=0.0;
			mol->atoms[i].velocity[2]=0.0;
			mol->atoms[i].charge = charge;
			mol->atoms[i].charge0 = charge;
			mol->atoms[i].electronegativity = 0;
			mol->atoms[i].hardness = 0;
			mol->atoms[i].width = mol->atoms[i].prop.covalentRadii;
			mol->atoms[i].mass = mol->atoms[i].prop.mass;
			mol->atoms[i].rho = 0.0;
	   		mol->atoms[i].U = 0.0;
			if(!get_connections_one_atom(t, mol->nAtoms, ibeg, mol->atoms[i].typeConnections))
			{
				/*
				printf("Sorry I cannot read the connection for atom # %d from the %s file.\n",i+1,namefile);
				exit(1);
				*/
				fprintf(stderr,"Warning : I cannot read the connection for atom # %d from the %s file.\n",i+1,namefile);
			}
			get_gradients_one_atom(t, mol->atoms[i].gradient);
	}
	
	fclose(file);
	/* if all freezed, set all to variable */
	{
		int j = 0;
		for(i=0;i<mol->nAtoms;i++)
			if(!mol->atoms[i].variable) j++;
		if(j==mol->nAtoms)
		for(i=0;i<mol->nAtoms;i++)
			mol->atoms[i].variable = TRUE;
	}
	return mol;
}
/********************************************************************************/
static Molecule* getMoleculeFromGabeditFile_FRCOORD(char* namefile)
{
	FILE* file = NULL;
	static char *t = NULL; 
	Molecule* mol = newMolecule();
#define SZ 50
	char symbol[SZ];
	double X,Y,Z;
	int i,l;

	if(!namefile) 
	{
		printf("Sorry I cannot read geometry namefile = NULL\n");
		exit(1);
	}
	if(t==NULL) t = malloc(BSIZE*sizeof(char));
	file = fopen(namefile,"rb");
	if(!file)
	{
		printf("Sorry I cannot open %s file\n",namefile);
		exit(1);
	}
	rewind(file);
	while(!feof(file))
  	{
    		if(!fgets(t,BSIZE, file)) break;
		deleteFirstSpaces(t);
		if(t[0]=='#') continue;
		uppercase(t);
		if(strstr(t,"[FR-COORD]")) break;
	}
	mol->nAtoms = 0;
	mol->spinMultiplicity = 1;
	mol->totalCharge = 0;
	i = -1;
	while(!feof(file))
  	{
    		if(!fgets(t,BSIZE, file)) break;
		deleteFirstSpaces(t);
		if(t[0]=='#') continue;
		uppercase(t);
		if(strstr(t,"[")) break;
    		if(4!=sscanf(t,"%s %lf %lf %lf", symbol, &X,&Y,&Z)) break;
		symbol[0]=toupper(symbol[0]);
		l=strlen(symbol);
		if (l>=2) symbol[1]=tolower(symbol[1]);
		if (l>=3) symbol[2]=tolower(symbol[2]);
		i++;
		mol->atoms = realloc(mol->atoms,(i+1)*sizeof(Atom));
		mol->atoms[i].prop = propAtomGet(symbol);
		mol->atoms[i].mmType=strdup(symbol);
		mol->atoms[i].pdbType=strdup(symbol);
		mol->atoms[i].residueName=strdup("UNK");
		mol->atoms[i].N=i+1;
		mol->atoms[i].layer=HIGH_LAYER;
		mol->atoms[i].variable=TRUE;
		mol->atoms[i].show=TRUE;
		mol->atoms[i].residueNumber=0;
		mol->atoms[i].coordinates[0] = X*BOHRTOANG;
		mol->atoms[i].coordinates[1] = Y*BOHRTOANG;
		mol->atoms[i].coordinates[2] = Z*BOHRTOANG;
		mol->atoms[i].gradient[0]=0.0;
		mol->atoms[i].gradient[1]=0.0;
		mol->atoms[i].gradient[2]=0.0;
		mol->atoms[i].velocity[0]=0.0;
		mol->atoms[i].velocity[1]=0.0;
		mol->atoms[i].velocity[2]=0.0;
		mol->atoms[i].charge = 0.0;
		mol->atoms[i].charge0 = 0.0;
		mol->atoms[i].electronegativity = 0;
		mol->atoms[i].hardness = 0;
		mol->atoms[i].width = mol->atoms[i].prop.covalentRadii;
		mol->atoms[i].mass = mol->atoms[i].prop.mass;
		mol->atoms[i].rho = 0.0;
	   	mol->atoms[i].U = 0.0;
	}
	mol->nAtoms = i+1;
	printf("nAtoms = %d\n",mol->nAtoms);
	fclose(file);
	/* if all freezed, set all to variable */
	{
		int j = 0;
		for(i=0;i<mol->nAtoms;i++) if(!mol->atoms[i].variable) j++;
		if(j==mol->nAtoms) for(i=0;i<mol->nAtoms;i++) mol->atoms[i].variable = TRUE;
	}
        for(i=0; i<mol->nAtoms; i++) mol->atoms[i].typeConnections = NULL;
	resetTypeConnections(mol);
	setConnections(mol);
	return mol;
}
/********************************************************************************/
static Molecule* getMoleculeFromGabeditFile(char* namefile)
{
	FILE* file = NULL;
	Molecule* mol = NULL;
	boolean geoms = FALSE;

	if(!namefile) 
	{
		printf("Sorry I cannot read geometry namefile = NULL\n");
		exit(1);
	}
	file = fopen(namefile,"rb");
	if(!file)
	{
		printf("Sorry I cannot open %s file\n",namefile);
		exit(1);
	}
        if(goToStr(file, "[GEOMS]")) geoms = TRUE;
	fclose(file);
	if(geoms) return getMoleculeFromGabeditFile_GEOMS(namefile);
	else return getMoleculeFromGabeditFile_FRCOORD(namefile);

	return mol;
}
/********************************************************************************/
Molecule* readMoleculeFromGabeditFile(char *fileName)
{
	Molecule* mol = getMoleculeFromGabeditFile(fileName);
	readVibrationFromGabeditFile(mol, fileName);
	readDipolesDerivativesFromGabeditFile(mol, fileName);
	return mol;
	
}
/********************************************************************************/
static void readBoxes(Molecule* mol, FILE* file, char* t)
{
	char* pos;
	rewind(file);
	mol->boxes.lengths[0] = mol->boxes.lengths[1] = mol->boxes.lengths[2]= 0.0;
	while(!feof(file))
  	{
    		if(!fgets(t,BSIZE, file)) break;
		deleteFirstSpaces(t);
		if(t[0]=='#') continue;
		uppercase(t);
		pos = strstr(t,"BOXES");
		if(pos && strstr(t,"="))
		{ 
			double x,y,z;
			double a,b,g;
			int n =0;
			pos = strstr(t,"=")+1;
			n =sscanf(pos,"%lf %lf %lf %lf %lf %lf",&x,&y,&z,&a,&b,&g);
			if(n==1) initBoxes(&mol->boxes, x,x,x,90.0,90.0,90.0); 
			if(n==2) initBoxes(&mol->boxes, x,y,y,90.0,90.0,90.0);
			if(n==3) initBoxes(&mol->boxes, x,y,z,90.0,90.0,90.0);
			if(n==6) initBoxes(&mol->boxes, x,y,z,a,b,g);
			break;
		}
	}
}
/********************************************************************************/
static void readWall(Molecule* mol, FILE* file, char* t)
{
	char* pos;
	rewind(file);
	initWallParameters(&mol->wall, 0, 1.0,2);
	while(!feof(file))
  	{
    		if(!fgets(t,BSIZE, file)) break;
		deleteFirstSpaces(t);
		if(t[0]=='#') continue;
		uppercase(t);
		pos = strstr(t,"WALL");
		if(pos && strstr(t,"="))
		{ 
			double E0;
			double rho;
			int nc;
			int n =0;
			pos = strstr(t,"=") + 1;
			n =sscanf(pos,"%lf %lf %d",&E0,&rho,&nc);
			//printf("t=%s\n",t);
			//printf("pos=%s\n",pos);
			if(n==3 && nc%2==0) initWallParameters(&mol->wall, E0, rho,nc);
			else  { 
				fprintf(stderr,"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n");
				fprintf(stderr,"Error during the reading of Wall parameters\n");
				fprintf(stderr,"You must give E0(au), rho (cutoff radius in angstrom) and nc(even integer)\n");
				fprintf(stderr,"Example : Wall=1000.0 10.0 6\n");
				fprintf(stderr,"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n");
				exit(1);
			}
			break;
		}
	}
}
/********************************************************************************/
static Molecule* readGeom(char* namefile)
{
	FILE* file = NULL;
	static char *t = NULL; 
	boolean Ok = FALSE;
	int n,ic,is;
	Molecule* mol = newMolecule();
#define SZ 50
	char symbol[SZ];
	char mmType[SZ];
	char pdbType[SZ];
	char residueName[SZ];
	double X,Y,Z;
	double charge;
	int layer;
	int i,l;
	char* pos;

	if(!namefile) 
	{
		printf("Sorry I cannot read geometry namefile = NULL\n");
		exit(1);
	}
	if(t==NULL) t = malloc(BSIZE*sizeof(char));
	file = fopen(namefile,"rb");
	if(!file)
	{
		printf("Sorry I cannot open %s file\n",namefile);
		exit(1);
	}
	rewind(file);
	while(!feof(file))
  	{
    		if(!fgets(t,BSIZE, file)) break;
		deleteFirstSpaces(t);
		if(t[0]=='#') continue;
		uppercase(t);
		pos = strstr(t,"GEOMETRY");
		if(pos)
		{ 
			while(!feof(file))
			{
    				if(!fgets(t,BSIZE, file)) break;
				deleteFirstSpaces(t);
				if(t[0]=='#') continue;
				break;
			}
			if(3==sscanf(t,"%d%d%d",&n,&ic,&is) && n>0 && is>0)
			{
				/*
				printf("Mult = %d\n",is);
				printf("nAtoms = %d\n",n);
				*/
				mol->nAtoms = n;
				mol->spinMultiplicity = is;
				mol->totalCharge = ic;
				mol->atoms = malloc(mol->nAtoms*sizeof(Atom));
				for(i=0; i<mol->nAtoms; i++) mol->atoms[i].typeConnections = malloc(mol->nAtoms*sizeof(int));
				Ok = TRUE;
			}
			break;
		}
	}
	if(!Ok)
	{
		printf("Sorry I cannot read geometry from %s file\n",namefile);
		exit(1);
	}
	for(i=0; i<mol->nAtoms; i++)
	{
			int variable = 0;
			int ibeg = 12;
			if(!fgets(t,BSIZE,file))
			{
				printf("Sorry I cannot read geometry from %s file.\n",namefile);
				exit(1);
			}
			deleteFirstSpaces(t);
			if(t[0]=='#') { i--;continue;}
    			sscanf(t,"%s %s %s %s %d %lf %d %d %lf %lf %lf",
					symbol,mmType,pdbType,residueName, 
					&mol->atoms[i].residueNumber,
					&charge,&layer,&variable,&X,&Y,&Z);
			symbol[0]=toupper(symbol[0]);
			l=strlen(symbol);
			if (l==2) symbol[1]=tolower(symbol[1]);

			mol->atoms[i].prop = propAtomGet(symbol);
			mol->atoms[i].mmType=strdup(mmType);
			mol->atoms[i].pdbType=strdup(pdbType);
			mol->atoms[i].residueName=strdup(residueName);
			mol->atoms[i].N=i+1;
			mol->atoms[i].layer=layer;
			mol->atoms[i].variable=variable;
			mol->atoms[i].show=TRUE;
			mol->atoms[i].coordinates[0] = X;
			mol->atoms[i].coordinates[1] = Y;
			mol->atoms[i].coordinates[2] = Z;
			mol->atoms[i].velocity[0] = 0;
			mol->atoms[i].velocity[1] = 0;
			mol->atoms[i].velocity[2] = 0;
			mol->atoms[i].charge = charge;
			mol->atoms[i].charge0 = charge;
			mol->atoms[i].electronegativity = 0;
			mol->atoms[i].hardness = 0;
			mol->atoms[i].width = mol->atoms[i].prop.covalentRadii;
			mol->atoms[i].mass = mol->atoms[i].prop.mass;
			mol->atoms[i].rho = 0.0;
	   		mol->atoms[i].U = 0.0;
			if(!get_connections_one_atom(t, mol->nAtoms, ibeg, mol->atoms[i].typeConnections))
			{
				/*
				printf("Sorry I cannot read the connection for atom # %d from the %s file.\n",i+1,namefile);
				exit(1);
				*/
				fprintf(stderr,"Warning : I cannot read the connection for atom # %d from the %s file.\n",i+1,namefile);
			}
	}
	readBoxes(mol, file, t);
	readWall(mol, file, t);
	
	fclose(file);
	return mol;
}
/********************************************************************************/
boolean readVelocities(Molecule* mol, char* namefile)
{
	FILE* file = NULL;
	static char *t = NULL; 
	boolean Ok = FALSE;
	char* pos = NULL;

	if(!namefile) 
	{
		printf("Sorry I cannot read velocities namefile = NULL !\n");
		exit(1);
	}
	if(!mol)
	{
		printf("Sorry I cannot read velocities from the %s file !\n",namefile);
		exit(1);
	}
	if(mol->nAtoms<1)
	{
		printf("Sorry I cannot read velocities, number of atoms = %d !\n",mol->nAtoms);
		exit(1);
	}
	if(t==NULL) t = malloc(BSIZE*sizeof(char));
	file = fopen(namefile,"rb");
	if(!file)
	{
		printf("Sorry I cannot open %s file\n",namefile);
		exit(1);
	}
	rewind(file);
	while(!feof(file) && !pos)
  	{
    		if(!fgets(t,BSIZE, file)) break;
		deleteFirstSpaces(t);
		if(t[0]=='#') continue;
		uppercase(t);
		pos = strstr(t,"VELOCITIES");
	}
	if(pos)
	{ 
		int ia = 0;
		double v[3] = {0};
		int ic = 0;
		Ok = TRUE;
		while(!feof(file))
		{
			while(!feof(file))
			{
    				if(!fgets(t,BSIZE, file)) { Ok = FALSE; break;}
				deleteFirstSpaces(t);
				if(t[0]=='#') continue;
				break;
			}
			if(3==sscanf(t,"%lf %lf %lf",&v[0], &v[1],&v[2]))
			{
				for(ic=0;ic<3;ic++) mol->atoms[ia].velocity[ic] = v[ic];
				ia++;
				if(ia==mol->nAtoms) break;
			}
			else {Ok = FALSE; break;}
		}
	}
	if(!Ok && pos)
	{
		printf("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n");
		printf("Sorry I cannot read velocities !\n");
		printf("Check your data in %s file\n",namefile);
		printf("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n");
		exit(1);
	}
	fclose(file);
	return Ok;
}
/********************************************************************************/
boolean readMasses(Molecule* mol, char* namefile)
{
	FILE* file = NULL;
	static char *t = NULL; 
	boolean Ok = FALSE;
	char* pos = NULL;

	if(!namefile) 
	{
		printf("Sorry I cannot read masses, namefile = NULL\n");
		exit(1);
	}
	if(!mol)
	{
		printf("Sorry I cannot read masses from the %s file\n",namefile);
		exit(1);
	}
	if(mol->nAtoms<1)
	{
		printf("Sorry I cannot read masses, number of atoms = %d !\n",mol->nAtoms);
		exit(1);
	}
	if(t==NULL) t = malloc(BSIZE*sizeof(char));
	file = fopen(namefile,"rb");
	if(!file)
	{
		printf("Sorry I cannot open %s file\n",namefile);
		exit(1);
	}
	rewind(file);
	while(!feof(file) && !pos)
  	{
    		if(!fgets(t,BSIZE, file)) break;
		deleteFirstSpaces(t);
		if(t[0]=='#') continue;
		uppercase(t);
		pos = strstr(t,"MASSES");
	}
	if(pos)
	{ 
		int ia = 0;
		double m = 0;
		Ok = TRUE;
		while(!feof(file))
		{
			while(!feof(file))
			{
    				if(!fgets(t,BSIZE, file)) { Ok = FALSE; break;}
				deleteFirstSpaces(t);
				if(t[0]=='#') continue;
				break;
			}
			if(1==sscanf(t,"%lf",&m))
			{
				if(m==0) mol->atoms[ia].mass = mol->atoms[ia].prop.mass;
				else if(m>0) mol->atoms[ia].mass = m;
				else
				{
					int im = (int)(-m);
					int j;
					for(j=0;j<mol->atoms[ia].prop.nIsotopes;j++)
					{
						if(im==mol->atoms[ia].prop.iMass[j]) {mol->atoms[ia].mass = mol->atoms[ia].prop.rMass[j]; break;}
					}
					if(j>=mol->atoms[ia].prop.nIsotopes) { 
					printf("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n");
					printf("I cannot find the mass of the isotope %d for atom %d(%s)\n",im,ia,mol->atoms[ia].prop.symbol);
						Ok = FALSE; break;
					}
				}
				ia++;
				if(ia==mol->nAtoms) break;
			}
			else {Ok = FALSE; break;}
		}
	}
	if(!Ok && pos)
	{
		printf("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n");
		printf("Sorry I cannot read masses !\n");
		printf("Check your data in %s file\n",namefile);
		printf("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n");
		exit(1);
	}
	fclose(file);
	return Ok;
}
/*************************************************************************************/
Molecule* getMolecule(char* fileName, boolean connections)
{
	int i;
	Molecule* molecule = readGeom(fileName);
	readVelocities(molecule,fileName);
	readMasses(molecule,fileName);

	/* if all freezed, set all to variable */
	{
		int j = 0;
		for(i=0;i<molecule->nAtoms;i++)
			if(!molecule->atoms[i].variable) j++;
		if(j==molecule->nAtoms)
		for(i=0;i<molecule->nAtoms;i++)
			molecule->atoms[i].variable = TRUE;
	}
	return molecule;
}
/*************************************************************************************/
CCHEMIFileType getTypeFile(char* fileName)
{
	char* t = NULL;
	CCHEMIFileType type = FILETYPE_UNKNOWN;
	FILE* file = fopen(fileName, "rb");

	t = malloc(BSIZE*sizeof(char));
	rewind(file);
    	fgets(t,BSIZE,file);
	uppercase(t);
        if(strstr(t, "ENTERING" )) type = FILETYPE_GAUSSIANOUT;
	else if(strstr( t, "[MOLDEN FORMAT]" )) type = FILETYPE_MOLDEN;
	else if(strstr( t, "[GABEDIT FORMAT]" )) type = FILETYPE_GABEDIT;
	else if(strstr( t, "GAMESS" )) type = FILETYPE_GAMESSOUT;
	else if(atoi(t)>0 && !strstr(t,"**********")) type = FILETYPE_XYZ;
	rewind(file);
	if(goToStr(file, "GEOMETRY"))
	{
		int n, ic, is;
		if(3==fscanf(file,"%d%d%d",&n,&ic,&is) && n>0 && is>0) type = FILETYPE_CCHEMI;
	}
	rewind(file);
	if(type == FILETYPE_UNKNOWN)
	{
		while(!feof(file))
		{
    			if(!fgets(t,BSIZE,file)) break;
			uppercase(t);
			if(strstr(t,"PROGRAM SYSTEM MOLPRO")) { type = FILETYPE_MOLPROOUT; break; }
			if(strstr(t,"GAMESS VERSION") || strstr(t,"PC GAMESS")) { type = FILETYPE_GAMESSOUT; break; }//FireFLy7.1
			if(strstr(t,"FIREFLY VERSION")) { type = FILETYPE_GAMESSOUT; break; } // FireFly8.1
			if(strstr(t,"Welcome to Q-Chem")) { type = FILETYPE_QCHEMOUT; break; }
			if(strstr(t,"Northwest Computational Chemistry Package")) { type = FILETYPE_NWCHEMOUT; break; }
			if(strstr(t,"TURBOMOLE GmbH")) { type = FILETYPE_TURBOMOLEOUT; break; }
        		if(strstr(t, "ENTERING GAUSSIAN" )) { type = FILETYPE_GAUSSIANOUT; break; }
		}
	}
	rewind(file);
	if(type == FILETYPE_UNKNOWN)
	{
		while(!feof(file))
		{
    			if(!fgets(t,BSIZE,file)) break;
			if(strstr(t,"$orca_hessian_file")) { type = FILETYPE_ORCAHESSIAN; break; }
		}
	}
	rewind(file);
	if(type == FILETYPE_UNKNOWN)
	{
		while(!feof(file))
		{
    			if(!fgets(t,BSIZE,file)) break;
			if(strstr(t,"* O   R   C   A *")) { type = FILETYPE_ORCAOUT; break; }
		}
	}
	rewind(file);
	if(type == FILETYPE_UNKNOWN)
	{
		while(!feof(file))
		{
    			if(!fgets(t,BSIZE,file)) break;
			if(strstr(t,"IOpCl")) { type = FILETYPE_GAUSSIANFCHK; break; }
		}
	}
	rewind(file);
	if(type == FILETYPE_UNKNOWN)
	{
		while(!feof(file))
		{
    			if(!fgets(t,BSIZE,file)) break;
			if(strstr(t,"VASP")) { type = FILETYPE_VASPOUT; break; }
		}
	}
	rewind(file);
	if(type == FILETYPE_UNKNOWN)
	{
		while(!feof(file))
		{
    			if(!fgets(t,BSIZE,file)) break;
			if(strstr(t,"GAMESS"))
			{
    				if(!fgets(t,BSIZE,file)) break;
				if(strstr(t,"FROM IOWA STATE UNIVERSITY"))
				type = FILETYPE_GAMESSOUT;
				break;
			}
		}
	}
	rewind(file);
	if(type == FILETYPE_UNKNOWN)
	{
    		if(fgets(t,BSIZE,file) && strstr(t,"START OF MOPAC FILE")) { type = FILETYPE_MOPACAUX;}
	}
	rewind(file);
	if(type == FILETYPE_UNKNOWN)
	{
  		while(!feof(file) )
  		{
 			if(!fgets(t, BSIZE, file))break;
			if( strstr(t,"MOPAC DONE")) { type = FILETYPE_MOPACOUT; break; }
			if( strstr(t,"MOPAC2012")) { type = FILETYPE_MOPACOUT; break; }
		}
	}
	rewind(file);
	if(type == FILETYPE_UNKNOWN)
	{
  		while(!feof(file) )
  		{
 			if(!fgets(t, BSIZE, file))break;
			if( strstr(t,"VARIABLE") && strstr(t,"FUNCTION")) { type = FILETYPE_MOPACSCAN; break; }
		}
	}
	rewind(file);
	if(type == FILETYPE_UNKNOWN)
	{
  		while(!feof(file) )    
  		{
 			if(!fgets(t, BSIZE, file))break;
                 	if(strstr(t,"INTRINSIC REACTION COORDINATE") ) { type = FILETYPE_MOPACIRC; break; }
		}
	}
 	free(t);
	fclose(file);
	return type;

	
}
/*************************************************************************************/
static void getFileNameToRead(char* fileName, char* fileNameToRead)
{
	char* t = malloc(BSIZE*sizeof(char));
	FILE* file = fopen(fileName, "rb");
	sprintf(fileNameToRead,"%s",fileName);
	if(goToStr(file, "GEOMETRY"))
	{
		int n, ic, is;
    		if(fgets(t,BSIZE,file) && 3!=sscanf(t,"%d%d%d",&n,&ic,&is)) sscanf(t,"%s",fileNameToRead);
	}
	free(t);
}
/*************************************************************************************/
Molecule* readMoleculeFromCChemIFile(char* fileName, boolean connections)
{
	Molecule* mol = getMolecule(fileName, connections);
	readVibrationFromGabeditFile(mol, fileName);
	return mol;
}
/*************************************************************************************/
static void getFreqMinFreqMax(char* fileName, double* freqMin, double* freqMax)
{
	FILE* file = fopen(fileName,"rb");
	*freqMin = -1.0;
	*freqMax = -1.0;
        readOneReal(file,"freqMin",freqMin);
        readOneReal(file,"freqMax",freqMax);
	fclose(file);
}
/********************************************************************************/
static void getChargesFromMolproOutputFile(Molecule* mol, FILE* file)
{
  	char t[BSIZE];
  	char dump[100];
  	char d[100];
  	char sign[10];
	int i;

/*
 Population analysis by basis function type

 Unique atom        s        p        d        f        g    Total    Charge
   1  O       3.64454  4.63505  0.01261  0.00000  0.00000  8.29219  - 0.29219
   2  H       0.76538  0.08852  0.00000  0.00000  0.00000  0.85390  + 0.14610
   3  H       0.76538  0.08852  0.00000  0.00000  0.00000  0.85390  + 0.14610
*/
  	while(!feof(file) )
	{
    		if(!fgets(t,BSIZE,file)) break;
		if(strstr(t,"Population analysis by basis function type"))
		{
    			if(!fgets(t,BSIZE,file)) break;
    			if(!fgets(t,BSIZE,file)) break;
    			if(!fgets(t,BSIZE,file)) break;
			for(i=0;i<mol->nAtoms;i++)
			{
    				if(!fgets(t,BSIZE,file)) break;
				if(sscanf(t,"%s %s %s %s %s %s %s %s %s %s",dump, dump ,dump, dump, dump, dump, dump, dump, sign, d)==6)
				{
					mol->atoms[i].charge = atof(d);
					if(strstr(sign,"-")) mol->atoms[i].charge = -mol->atoms[i].charge;
				}
				else break;
			}
			//break; // To read last charges
		}
	}
}

/********************************************************************************/
static int compute_molpro_modes_masses(Molecule* mol, FILE* file)
{
	int i,j,c;
	rewind(file);
 	char t[BSIZE];
 	char sdum1[100];
 	char sdum2[100];
 	char sdum3[100];
	double mass;

 	while(!feof(file))
	{
    		fgets(t,BSIZE,file);
 		// Atom  1: H         Mass   1.00794
	 	if (strstr(t,"Atom") && strstr( t,"Mass") )
		{
			for(i=0;i<strlen(t);i++) if(t[i]==':') t[i]=' ';
			int nf = sscanf(t,"%s %d %s %s %lf", sdum1,&j, sdum2,sdum3,&mass);
			if(nf==5 && j>=1 && j<=mol->nAtoms) mol->atoms[j-1].mass = mass;
		}
	}
        /* compute the effective masses */
        /* modes already mass-weighted in Molpro output file*/
        for(i=0;i<mol->vibration.nModes;i++)
        {
                double m = 0;
                for(j=0;j<mol->nAtoms;j++)
                {
                        double r2 = 0;
                        for(c=0;c<3;c++) r2 += mol->vibration.modes[i].vectors[c][j]*mol->vibration.modes[i].vectors[c][j];
                        //m += r2/(mol->atoms[j].mass);// don't divide by mass[j]
                        m += r2;
                }
                if(m<=0) m = 1;
                m = 1.0/m;
                for(j=0;j<mol->nAtoms;j++)
                {
                        //double r =sqrt(m)/sqrt(mol->atoms[j].mass); // don't divide by mass[j]
                        double r =sqrt(m);
                        for(c=0;c<3;c++) mol->vibration.modes[i].vectors[c][j] *= r;
                }
		mol->vibration.modes[i].mass = m;
        }
	return 0;
}
/********************************************************************************/
static int read_molpro_modes_str(Molecule* mol, FILE* file, char* str, int jmin)
{
 	char t[BSIZE];
 	char sdum1[BSIZE];
 	char sdum2[BSIZE];
 	boolean OK;
 	int taille=BSIZE;
	int i;
	int j;
	int k;
	int c;
	int ne;
	int nf;
	double freq[5];
	double IRIntensity[5];

	if(mol->nAtoms<1) return 2;

 	OK=FALSE;
 	while(!feof(file))
	{
    		fgets(t,taille,file);
	 	if ( strstr( t,str) )
	  	{
			OK = TRUE;
			break;
	  	}
	}

	if(!OK) return 1;

	j = jmin;
  	while(!feof(file))
  	{
		if(!fgets(t,taille,file)) { return 2; }
		if(atof(t)==0) if(!fgets(t,taille,file)) break;
    		if(!fgets(t,taille,file)) break;
		if(!strstr(t,"Wavenumbers")) break;
		nf = sscanf(t,"%s %s %lf %lf %lf %lf %lf", sdum1,sdum2, &freq[0],&freq[1],&freq[2],&freq[3],&freq[4]);
		nf -= 2;
		if(strstr(str,"imaginary")) for(k=0;k<nf;k++) freq[k] = -freq[k];
    		if(!fgets(t,taille,file)) break;
		sscanf(t,"%s %s %lf %lf %lf %lf %lf",
				sdum1,sdum2,
				&IRIntensity[0],&IRIntensity[1],&IRIntensity[2],&IRIntensity[3],&IRIntensity[4]);

    		if(!fgets(t,taille,file)) break;

		for(k=0;k<nf;k++)
		{
			mol->vibration.modes[j+k].frequency = freq[k];
			mol->vibration.modes[j+k].properties[0] = IRIntensity[k];
			mol->vibration.modes[j+k].mass = 1.0;
		}

		for(i=0;i<mol->nAtoms;i++)
		{
			for(c=0;c<3;c++)
			{
    				if(!fgets(t,taille,file)) break;
				ne = sscanf(t,"%s %lf %lf %lf %lf %lf",sdum1,
					&mol->vibration.modes[j  ].vectors[c][i],
					&mol->vibration.modes[j+1].vectors[c][i],
					&mol->vibration.modes[j+2].vectors[c][i],
					&mol->vibration.modes[j+3].vectors[c][i],
					&mol->vibration.modes[j+4].vectors[c][i]);
				if(ne <2) { return 2; }
			}

		}
		j+=nf;

	}
	mol->vibration.nModes = j;
	return 0;
}
/********************************************************************************/
static boolean readVibrationFromMolproOutputFile(Molecule* mol, char* fileName)
{
	int normalModes = 0;
	int normalModesImag = 0;
 	FILE* file = fopen(fileName, "rb");
	int nProps = 1;
	if(file ==NULL)
	{
		fprintf(stderr,"Sorry\nI can not open %s Molpro output file\n",fileName);
		return FALSE;
	}
	if(mol->nAtoms<1) 
	{
		if(file) fclose(file);
		return FALSE;
	}

	initVibrations(mol, 3*mol->nAtoms, nProps);
	normalModes = read_molpro_modes_str(mol, file, "Normal Modes",0);
	if(normalModes<=1) normalModesImag = read_molpro_modes_str(mol, file, "Normal Modes of imaginary frequencies",mol->vibration.nModes);
	if(normalModes>1 || normalModesImag>1)
	{
			char buffer[BSIZE];
			freeVibrations(mol);
  			fprintf(stderr,"%s",buffer);
			if(file) fclose(file);
			return FALSE;
	}
/*
	mol->vibration.nModes = j;
        mol->vibration.modes = realloc(mol->vibration.modes, mol->vibration.nModes*sizeof(VibMode));
	for(j=0;j< mol->vibration.nModes;j++)
	{
			for(i=0;i< mol->nAtoms;i++)
			for(c=0;c<3;c++)
			mol->vibration.modes[j].vectors[c][i] *= sqrt(mol->vibration.modes[j].mass);
	}
*/
	compute_molpro_modes_masses(mol, file);
	sortFrequencies(mol);
	removeTransRotModes(mol);
	if(file) fclose(file);
	return 1;
}
/********************************************************************************/
Molecule* getMoleculeFromMolproOutputFile(char *fileName, int numgeometry)
{
	char *t;
	boolean OK;
	char *AtomCoord[5];
	FILE *file;
	int i;
	int j=0;
	int l;
	int numgeom;
	char dum[100];
	Molecule* mol = NULL;
	int kk;
	char AtomCharge[100];
	int idummy;



	for(i=0;i<5;i++) AtomCoord[i]=malloc(BSIZE*sizeof(char));
  
	t=malloc(BSIZE*sizeof(char));
 	file = fopen(fileName, "rb");
	if(file ==NULL)
	{
		free(t);
		printf(("Sorry\nI can not open %s Molpro output file\n"),fileName);
		exit(1);
		return NULL;
	}
	numgeom = 0;
	do 
	{
		OK=FALSE;
		while(!feof(file)){
			fgets(t,BSIZE,file);
			uppercase(t);
			if ( strstr(t," ATOMIC COORDINATE"))
			{
	  			fgets(t,BSIZE,file);
	  			fgets(t,BSIZE,file);
	  			fgets(t,BSIZE,file);
 				numgeom++;
				if((int)numgeom == numgeometry ) { OK = TRUE; break; }
				if(numgeometry<0 ) { OK = TRUE; break; }
	  		}
		}
		if(!OK && (numgeom == 0) ){
			free(t);
			printf(("Sorry\nI can not open read geometry from %s Molpro output file\n"),fileName);
			exit(1);
			return NULL;
		}
		if(!OK)break;

		//printf("Begin freMole\n");
		if(mol) freeMolecule(mol);
		mol = newMolecule();
		mol->atoms = NULL;
		mol->nAtoms = 0;
		//printf("numgeom=%d numgeometry=%d\n",numgeom,numgeometry);

		j=-1;
		while(!feof(file) )
		{
			fgets(t,BSIZE,file);
			strDeleten(t);
			if (isABackspace(t)) break;
			if ( !strcmp(t,"\n")) break;
			if ( !strcmp(t,"\r\n")) break;
			j++;
			//printf("j=%d t = %s\n",j,t);
			//sscanf(t,"%s %s %s %s %s",AtomCoord[0],dum, AtomCoord[1], AtomCoord[2],AtomCoord[3]);
			kk = sscanf(t,"%d %s %s %s %s %s %s",&idummy, AtomCoord[0],AtomCharge,AtomCoord[1], AtomCoord[2],AtomCoord[3], dum);
                        if(kk==7) sscanf(t,"%d %s %s %s %s %s %s",&idummy, AtomCoord[0],AtomCharge,dum, AtomCoord[1], AtomCoord[2],AtomCoord[3]);
			{
				int k;
				for(k=0;k<(int)strlen(AtomCoord[0]);k++) if(isdigit(AtomCoord[0][k])) AtomCoord[0][k] = ' ';
				deleteAllSpaces(AtomCoord[0]);
			}

			AtomCoord[0][0]=toupper(AtomCoord[0][0]);
			l=strlen(AtomCoord[0]);
			if (l==2) AtomCoord[0][1]=tolower(AtomCoord[0][1]);
			if (atoi(AtomCharge) == 0) sprintf(AtomCoord[0],"X");
			mol->atoms = realloc(mol->atoms,(j+1)*sizeof(Atom));

			mol->atoms[j].prop = propAtomGet(AtomCoord[0]);
			mol->atoms[j].mmType=strdup(AtomCoord[0]);
			mol->atoms[j].pdbType=strdup(AtomCoord[0]);
			mol->atoms[j].residueName=strdup(AtomCoord[0]);
			mol->atoms[j].N=j+1;
			mol->atoms[j].layer=HIGH_LAYER;
			mol->atoms[j].variable=TRUE;
			mol->atoms[j].show=TRUE;
			mol->atoms[j].residueNumber=0;
	   		mol->atoms[j].charge=0.0;
	   		mol->atoms[j].charge0=0.0;
	   		mol->atoms[j].electronegativity=0.0;
	   		mol->atoms[j].hardness=0.0;
	   		mol->atoms[j].width=mol->atoms[j].prop.covalentRadii;
	   		mol->atoms[j].mass= mol->atoms[j].prop.mass;
	   		mol->atoms[j].rho=0.0;
	   		mol->atoms[j].U = 0.0;

			mol->atoms[j].coordinates[0]=atof(AtomCoord[1])*BOHRTOANG;
			mol->atoms[j].coordinates[1]=atof(AtomCoord[2])*BOHRTOANG;
			mol->atoms[j].coordinates[2]=atof(AtomCoord[3])*BOHRTOANG;

			mol->atoms[j].gradient[0]=0.0;
			mol->atoms[j].gradient[1]=0.0;
			mol->atoms[j].gradient[2]=0.0;

			mol->atoms[j].velocity[0]=0.0;
			mol->atoms[j].velocity[1]=0.0;
			mol->atoms[j].velocity[2]=0.0;
        		mol->atoms[j].typeConnections = NULL;
			if(j==mol->nAtoms-1) break;
		}
		if(OK && numgeometry>=0) break;
		mol->nAtoms = j+1;
	}while(!feof(file));
	//printf("end read molecule\n");

	rewind(file);
	getChargesFromMolproOutputFile(mol, file);
	fclose(file);
	free(t);
	for(i=0;i<5;i++) free(AtomCoord[i]);
	//printf("begin resetTypeConnections\n");
        for(i=0; i<mol->nAtoms; i++) mol->atoms[i].typeConnections = NULL;
	resetTypeConnections(mol);
	//printf("end resetTypeConnections\n");
	int connections = 1;
	if(connections) setConnections(mol);
	else {
		createBondedMatrix(mol);
		printf(("Establishing connectivity : non bonded ...\n"));
		setNonBondedConnections(mol);
		freeBondedMatrix(mol);
	}

	/* if all freezed, set all to variable */
	{
		int j = 0;
		for(i=0;i<mol->nAtoms;i++)
			if(!mol->atoms[i].variable) j++;
		if(j==mol->nAtoms)
		for(i=0;i<mol->nAtoms;i++)
			mol->atoms[i].variable = TRUE;
	}
	return mol;
}
/********************************************************************************/
Molecule* readMoleculeFromMolproOutputFile(char *fileName, int numgeometry)
{
	//printf("Begin redMolecule from Molpro output file\n");
	Molecule* mol = getMoleculeFromMolproOutputFile(fileName, numgeometry);
	//printf("End redMolecule from Molpro output file\n");
	readVibrationFromMolproOutputFile(mol, fileName);
	return mol;
	
}
/*************************************************************************************/
Molecule* readMolecule(char* fileName, boolean connections)
{
	char* fileNameToRead = malloc(BSIZE*sizeof(char));
	CCHEMIFileType type;
	Molecule* mol = NULL;
	double freqMin;
	double freqMax;

	if(!fileName) return mol;

	fileNameToRead = malloc(BSIZE*sizeof(char));
	getFileNameToRead(fileName, fileNameToRead);

	type = getTypeFile(fileNameToRead);

	//printf("type %d\n",type);

	if(type == FILETYPE_CCHEMI) mol = readMoleculeFromCChemIFile(fileNameToRead, connections);
	else if(type == FILETYPE_GABEDIT) mol = readMoleculeFromGabeditFile(fileNameToRead);
	else if(type == FILETYPE_GAUSSIANOUT) mol = readMoleculeFromGaussianOutputFile(fileNameToRead,-1);
	else if(type == FILETYPE_GAUSSIANFCHK) mol = readMoleculeFromGaussianFChkFile(fileNameToRead);
	else if(type == FILETYPE_GAMESSOUT) mol = readMoleculeFromGamessOutputFile(fileNameToRead,-1);
	else if(type == FILETYPE_MOPACOUT) mol = readMoleculeFromMopacOutputFile(fileNameToRead,-1);
	else if(type == FILETYPE_MOPACAUX) mol = readMoleculeFromMopacAuxFile(fileNameToRead,-1);
	else if(type == FILETYPE_ORCAOUT) mol = readMoleculeFromOrcaOutputFile(fileNameToRead,-1);
	else if(type == FILETYPE_ORCAHESSIAN) mol = readMoleculeFromOrcaHessianFile(fileNameToRead);
	else if(type == FILETYPE_MOLPROOUT) mol = readMoleculeFromMolproOutputFile(fileNameToRead,-1);
	else
	{
		fprintf(stderr,"Sorry, I cannot read the molecule from %s file\n",fileNameToRead);
		exit(1);
	}

	if(connections) setConnections(mol);
	else {
		createBondedMatrix(mol);
		printf(("Establishing connectivity : non bonded ...\n"));
		setNonBondedConnections(mol);
		freeBondedMatrix(mol);
	}
	getFreqMinFreqMax(fileName, &freqMin, &freqMax);
	if(freqMin>-1 || freqMax>-1) mol->klass->removeFrequencies(mol, freqMin,freqMax);
	//if(type == FILETYPE_GAUSSIANOUT ) already read with ir intenisties
	//if(type == FILETYPE_ORCAOUT ) .... 
	//if(type == FILETYPE_ORCAHESSIAN )  already read with ir intenisties
	//if(type == FILETYPE_GABEDIT) already read with ir intenisties

	free(fileNameToRead);
	return mol;
}
/********************************************************************************/
void exchangeGeometries(Molecule* mol1, Molecule* mol2)
{
	int i,j;
	if(mol1->nAtoms!=mol2->nAtoms)
	{
		printf("Sorry I cannot exahnge the geometries... not same number of atoms\n");
		exit(1);
	}
	for(i=0;i<mol1->nAtoms;i++)
	for(j=0;j<3;j++)
	{
		double t = mol1->atoms[i].coordinates[j];
		mol1->atoms[i].coordinates[j]= mol2->atoms[i].coordinates[j];
		mol2->atoms[i].coordinates[j]= t;
	}
}
/*****************************************************************************/
static boolean addGeometryToGabedit(Molecule* molecule,FILE* file)
{
	int j,k;
	int nc;

	if(!molecule) return FALSE;
	fprintf(file,"%d %d %d %lf\n",molecule->nAtoms, molecule->totalCharge, molecule->spinMultiplicity,molecule->potentialEnergy);
	for(j=0;j<molecule->nAtoms;j++)
	{
		nc = 0;
		for(k=0;k<molecule->nAtoms;k++) if(molecule->atoms[j].typeConnections[k]>0) nc++;
		fprintf(file," %s %s %s %s %d %0.12lf %d %d %0.12lf %0.12lf %0.12lf %d ", 
				molecule->atoms[j].prop.symbol,
				molecule->atoms[j].mmType,
				molecule->atoms[j].pdbType,
				molecule->atoms[j].residueName,
				molecule->atoms[j].residueNumber,
				molecule->atoms[j].charge,
				molecule->atoms[j].layer,
				molecule->atoms[j].variable,
				molecule->atoms[j].coordinates[0],
				molecule->atoms[j].coordinates[1],
				molecule->atoms[j].coordinates[2],
				nc
				);
		for(k=0;k<molecule->nAtoms;k++) 
		{
	 		int nk = molecule->atoms[k].N-1;
			if(molecule->atoms[j].typeConnections[nk]>0) 
			fprintf(file," %d %d", nk+1,molecule->atoms[j].typeConnections[nk]);
		}
		save_gradients_one_atom(file, molecule->atoms[j].gradient);
		fprintf(file,"\n");
	}
	return TRUE;

}
/*****************************************************************************/
static boolean addFirstDerivativeToFile(Molecule* molecule, FILE* file)
{
	if(molecule->vibration.nModes>0 && molecule->vibration.nProperties>=5)
	{
		int i,j;
		int nOld = molecule->vibration.nProperties-3;
		fprintf(file,"\n#xyz    i       Values[au cm^1/2]\n");
		fprintf(file,"First derivatives\n");
		
		for(i=0;i<molecule->vibration.nModes;i++) 
		{
			char axes[] = {'X','Y','Z'};
			for(j=0;j<3;j++)
			if(fabs(molecule->vibration.modes[i].properties[j+nOld])>1e-10)
				fprintf(file,"%c\t%d\t%0.10lf\n", axes[j], i+1, molecule->vibration.modes[i].properties[j+nOld]);
		}
		fprintf(file,"END\n");
		return TRUE;
	}

	return FALSE;
}
/*****************************************************************************/
static boolean addVibrationToFile(Molecule* molecule, FILE* file)
{
	int i,j,k;

	if(!molecule) return FALSE;
	if(molecule->vibration.nModes<1) return FALSE;

	fprintf(file,"[FREQ]\n");
	for(i=0;i<molecule->vibration.nModes;i++) 
		fprintf(file,"%f\n", molecule->vibration.modes[i].frequency);

	if(molecule->vibration.nProperties>0)
	{
		fprintf(file,"[INT]\n");
		for(i=0;i<molecule->vibration.nModes;i++) 
		{
			for(j=0;j<molecule->vibration.nProperties;j++) 
				fprintf(file,"%0.14f ", molecule->vibration.modes[i].properties[j]);
			fprintf(file,"\n");
		}
	}

	fprintf(file,"[MASS]\n");
	for(i=0;i<molecule->vibration.nModes;i++) 
		fprintf(file,"%0.14f\n", molecule->vibration.modes[i].mass);

	fprintf(file,"[FR-COORD]\n");
	for(j=0;j<molecule->nAtoms;j++)
	{
		fprintf(file," %s %0.14f %0.14f %0.14f\n", 
			molecule->atoms[j].prop.symbol,
			molecule->atoms[j].coordinates[0]*ANGTOBOHR,
			molecule->atoms[j].coordinates[1]*ANGTOBOHR,
			molecule->atoms[j].coordinates[2]*ANGTOBOHR
			);
	}
	fprintf(file,"[FR-NORM-COORD]\n");
	for(i=0;i<molecule->vibration.nModes;i++) 
	{
		fprintf(file,"vibration %d\n",i+1);
		for(j=0;j<molecule->nAtoms;j++)
		{
			for(k=0;k<3;k++) fprintf(file,"%0.14f ",  molecule->vibration.modes[i].vectors[k][j]);
			fprintf(file,"\n");
		}
		
	}
	molecule->klass->addFirstDerivativeToFile(molecule, file);
	return TRUE;

}
/*****************************************************************************/
static boolean saveMoleculeTypeSave(char* fileName, char* typeSave, Molecule* molecule)
{
	FILE* file = NULL;
	int j;
	int form = 1;

	printf("Save molecule in %s\n",fileName);
	if(!molecule) return FALSE;

 	file = fopen(fileName, typeSave);

	if(!file) return FALSE;

	fprintf(file,"[Gabedit Format]\n");
	fprintf(file,"[GEOCONV]\n");
	fprintf(file,"energy\n");
	fprintf(file,"%f\n",molecule->potentialEnergy);
	fprintf(file,"max-force\n");
	fprintf(file,"%f\n",0.0);
	fprintf(file,"rms-force\n");
	fprintf(file,"%f\n",0.0);

	fprintf(file,"\n");
	fprintf(file,"[GEOMETRIES]\n");
	{
		fprintf(file,"%d\n",molecule->nAtoms);
		fprintf(file,"\n");
		for(j=0;j<molecule->nAtoms;j++)
		fprintf(file," %s %0.14f %0.14f %0.14f\n", 
				molecule->atoms[j].prop.symbol,
				molecule->atoms[j].coordinates[0],
				molecule->atoms[j].coordinates[1],
				molecule->atoms[j].coordinates[2]
				);
	}
	fprintf(file,"\n");
	fprintf(file,"[GEOMS] %d\n",form);
	fprintf(file,"%d 3\n",1);
	fprintf(file,"energy kcal/mol 1\n");
	fprintf(file,"deltaE K 1\n");
	fprintf(file,"Dipole Debye 3\n");
	//molecule->klass->computeDipole(molecule);
	{
		fprintf(file,"%0.14f\n",molecule->potentialEnergy);
		fprintf(file,"0\n");
		fprintf(file,"%0.14f %0.14f %0.14f\n",molecule->dipole[0],molecule->dipole[1],molecule->dipole[2]);
		molecule->klass->addGeometryToGabedit(molecule,file);
	}
	addVibrationToFile(molecule, file);
	fclose(file);
	return TRUE;

}
/*****************************************************************************/
static boolean saveMolecule(Molecule* molecule, char* fileName)
{
	return saveMoleculeTypeSave(fileName, "w", molecule);
}
/******************************************************************************/
static void save_atom_hin_file(FILE* file, char* name, int atomNumber, char* atomPDBType, char* atomMMType,
		double x, double y, double z, char* symbol, double charge,
		int N, int* connection, int* connectionType)
{
	int i;
        fprintf(file,"%s %d ",name,atomNumber);
        fprintf(file,"%s ",atomPDBType);
        fprintf(file,"%s ",symbol); 
        fprintf(file,"%s - ",atomMMType); 
        fprintf(file,"%0.14f ",charge); 
        fprintf(file,"%0.14f ",x); 
        fprintf(file,"%0.14f ",y); 
        fprintf(file,"%0.14f ",z); 
	if(N>0)
	{
        	fprintf(file,"%d ",N); 
		for(i=0;i<N;i++)
		{
			if(connectionType[i]==3) fprintf(file,"%d t ",connection[i]); 
			else if(connectionType[i]==2) fprintf(file,"%d d ",connection[i]); 
			else fprintf(file,"%d s ",connection[i]); 
		}
	}
        fprintf(file,"\n"); 
}
/******************************************************************************/
static boolean saveHIN(Molecule* mol, char* fileName)
{
	int i,n,k;
 	FILE* file = fopen(fileName, "w");
	int nAtoms;
	Atom* atoms;
	int* connection;
	int* connectionType;
	if(!file) return FALSE;
	nAtoms = mol->nAtoms;
	atoms = mol->atoms;

	fprintf(file,"forcefield Amber99\n");
	fprintf(file,"sys 0 0 1\n");
	fprintf(file,"view 40 0.1272 55 15 0.247224 0.3713666 0.8949677 -0.8641704 0.5022867 0.0302929 -0.4382806 -0.7808937 0.4451014 6.191 0.64575 -54.754\n");
	fprintf(file,"seed -1108\n");
	fprintf(file,"mol 1\n");

	connection = malloc(nAtoms*sizeof(int));
	connectionType = malloc(nAtoms*sizeof(int));

	for(i=0;i<nAtoms;i++)
	{
		n = 0;
		if(atoms[i].typeConnections)
		for(k=0;k<nAtoms;k++)
		{
			if(i==k) continue;
			if(atoms[i].typeConnections[k]>0)
			{
				connection[n] = k+1;
				connectionType[n] = atoms[i].typeConnections[k];
				n++;
			}
		}
		save_atom_hin_file(file,"ATOM",i+1,atoms[i].pdbType, atoms[i].mmType,
		atoms[i].coordinates[0], atoms[i].coordinates[1], atoms[i].coordinates[2],
		atoms[i].prop.symbol, atoms[i].charge,n,connection, connectionType);
	}
	fprintf(file,"endmol 1\n");
	fclose(file);
	free(connection);
	free(connectionType);
	return TRUE;
}

/*****************************************************************************/
static boolean saveMol2(Molecule* mol, char* fileName)
{
	int i,n,j;
 	FILE* file = fopen(fileName, "w");
	if(!file) return FALSE;
	n = 0;
	for(i=0;i<mol->nAtoms;i++)
        if(mol->atoms[i].typeConnections)
        for(j=i+1;j<mol->nAtoms;j++)
                if(mol->atoms[i].typeConnections[j]) n++;

	fprintf(file,"@<TRIPOS>MOLECULE\n");
	fprintf(file,"MOL2  : Made in CChemI. mol2 file\n");
	fprintf(file," %10d %10d %10d\n",mol->nAtoms,n,1);
	fprintf(file,"SMALL\n");
	fprintf(file,"GASTEIGER\n");
	fprintf(file,"\n");
	fprintf(file,"@<TRIPOS>ATOM\n");
      	for (i=0;i<mol->nAtoms;i++)
	{
       		fprintf(file,"%7d%1s%-6s%12.4f%10.4f%10.4f%1s%-5s%4d%1s%-8s%10.4f\n",
               		i+1,"",mol->atoms[i].prop.symbol,
				mol->atoms[i].coordinates[0],
				mol->atoms[i].coordinates[1],
				mol->atoms[i].coordinates[2],
				"",mol->atoms[i].prop.symbol,1," ","LIG111",mol->atoms[i].charge);
	}
	fprintf(file,"@<TRIPOS>BOND\n");
	n = 0;
	for(i=0;i<mol->nAtoms;i++)
        if(mol->atoms[i].typeConnections)
        for(j=i+1;j<mol->nAtoms;j++)
                if(mol->atoms[i].typeConnections[j]) 
		{
			n++;
			fprintf(file,"%6d%6d%6d%3s%2d\n",n+1, i+1, j+1, "",mol->atoms[i].typeConnections[j]);
		}


	fclose(file);
	return TRUE;
}
/*****************************************************************************/
static boolean addGeometry(Molecule* molecule,FILE* file)
{
	int j,k;
	int nc;

	if(!molecule) return FALSE;

	fprintf(file,"# Geometry, nAtoms, charge, spin multiplicity.\n");
	fprintf(file,"# symbol, MMType, pdbType, residueName, numResidue, charge, layer, variable, x(Ang),y,z, nconn, num1, typ1, num2, typ2,...\n");
	fprintf(file,"Geometry\n");
	fprintf(file,"%d %d %d\n",molecule->nAtoms, molecule->totalCharge, molecule->spinMultiplicity);
	for(j=0;j<molecule->nAtoms;j++)
	{
		nc = 0;
		for(k=0;k<molecule->nAtoms;k++) if(molecule->atoms[j].typeConnections[k]>0) nc++;
		fprintf(file," %s %s %s %s %d %0.12lf %d %d %0.12lf %0.12lf %0.12lf %d ", 
				molecule->atoms[j].prop.symbol,
				molecule->atoms[j].mmType,
				molecule->atoms[j].pdbType,
				molecule->atoms[j].residueName,
				molecule->atoms[j].residueNumber,
				molecule->atoms[j].charge,
				molecule->atoms[j].layer,
				molecule->atoms[j].variable,
				molecule->atoms[j].coordinates[0],
				molecule->atoms[j].coordinates[1],
				molecule->atoms[j].coordinates[2],
				nc
				);
		for(k=0;k<molecule->nAtoms;k++) 
		{
	 		int nk = molecule->atoms[k].N-1;
			if(molecule->atoms[j].typeConnections[nk]>0) 
			fprintf(file," %d %d", nk+1,molecule->atoms[j].typeConnections[nk]);
		}
		fprintf(file,"\n");
	}
	return TRUE;

}
/*****************************************************************************/
static boolean addMolecule(Molecule* molecule,FILE* file)
{
	int j,k;
	int nc;

	if(!molecule) return FALSE;

	fprintf(file,"Geometry\n");
	fprintf(file,"%d %d %d\n",molecule->nAtoms, molecule->totalCharge, molecule->spinMultiplicity);
	for(j=0;j<molecule->nAtoms;j++)
	{
		nc = 0;
		for(k=0;k<molecule->nAtoms;k++) if(molecule->atoms[j].typeConnections[k]>0) nc++;
		fprintf(file," %s %s %s %s %d %0.12lf %d %d %0.12lf %0.12lf %0.12lf %d ", 
				molecule->atoms[j].prop.symbol,
				molecule->atoms[j].mmType,
				molecule->atoms[j].pdbType,
				molecule->atoms[j].residueName,
				molecule->atoms[j].residueNumber,
				molecule->atoms[j].charge,
				molecule->atoms[j].layer,
				molecule->atoms[j].variable,
				molecule->atoms[j].coordinates[0],
				molecule->atoms[j].coordinates[1],
				molecule->atoms[j].coordinates[2],
				nc
				);
		for(k=0;k<molecule->nAtoms;k++) 
		{
	 		int nk = molecule->atoms[k].N-1;
			if(molecule->atoms[j].typeConnections[nk]>0) 
			fprintf(file," %d %d", nk+1,molecule->atoms[j].typeConnections[nk]);
		}
		fprintf(file,"\n");
	}
	return TRUE;

}
/*****************************************************************************/
static boolean addVelocities(Molecule* molecule,FILE* file)
{
	int j;

	if(!molecule) return FALSE;
	fprintf(file,"# Velocities, vx(Ang/AKMA-time) vy(Ang/AKMA-time) vz(Ang/AKMA-time): 1fs = %f AKMA\n",fsInAKMA);
	fprintf(file,"Velocities\n");
	for(j=0;j<molecule->nAtoms;j++)
	{
		fprintf(file,"%0.12le %0.12le %0.12le\n", 
				molecule->atoms[j].velocity[0],
				molecule->atoms[j].velocity[1],
				molecule->atoms[j].velocity[2]
				);
	}
	return TRUE;
}
/*****************************************************************************/
static boolean saveFrequencies(Molecule* molecule, char* fileName, int nModes, double* frequencies, double** modes, double* reducedMasses, double* IRIntensities)
{
	FILE* file = NULL;
	int i,j;

	printf("Save molecule in %s\n",fileName);
	if(!molecule) return FALSE;

 	file = fopen(fileName, "w");

	if(!file) return FALSE;

	fprintf(file,"[Gabedit Format]\n");
	fprintf(file,"[Atoms] Angs\n");
	for(j=0;j<molecule->nAtoms;j++)
	{
		fprintf(file," %s %d %d %0.8f %0.8f %0.8f\n", 
			molecule->atoms[j].prop.symbol,
			j+1,
			molecule->atoms[j].prop.atomicNumber,
			molecule->atoms[j].coordinates[0],
			molecule->atoms[j].coordinates[1],
			molecule->atoms[j].coordinates[2]
			);
	}
	fprintf(file,"[FREQ]\n");
	for(i=0;i<nModes;i++) fprintf(file,"%f\n", frequencies[i]);
	fprintf(file,"[INT]\n");
	for(i=0;i<nModes;i++) fprintf(file,"%f\n", IRIntensities[i]);
	fprintf(file,"[MASS]\n");
	for(i=0;i<nModes;i++) fprintf(file,"%f\n", reducedMasses[i]);
	fprintf(file,"[FR-COORD]\n");
	for(j=0;j<molecule->nAtoms;j++)
	{
		fprintf(file," %s %f %f %f\n", 
			molecule->atoms[j].prop.symbol,
			molecule->atoms[j].coordinates[0]*ANGTOBOHR,
			molecule->atoms[j].coordinates[1]*ANGTOBOHR,
			molecule->atoms[j].coordinates[2]*ANGTOBOHR
			);
	}
	fprintf(file,"[FR-NORM-COORD]\n");
	for(i=0;i<nModes;i++)
	{
		fprintf(file,"vibration %d\n",i+1);
		for(j=0;j<molecule->nAtoms;j++)
			fprintf(file,"%f %f %f\n", modes[3*j+0][i], modes[3*j+1][i], modes[3*j+2][i]);
		
	}
	fclose(file);
	saveMoleculeTypeSave(fileName, "a", molecule);
	return TRUE;

}
/****************************************************************************************************************************************************/
static void removeTranslation(Molecule* molecule)
{
	double vtot[3] = {0,0,0};
	int i;
	int j;
	double mass = 1.0;
	double totMass = 0.0;
	Atom* atoms = molecule->atoms;
	int nAtoms = molecule->nAtoms;
	for ( i = 0; i < nAtoms; i++)
	{
		mass = atoms[i].mass;
		totMass += mass;
		for ( j = 0; j < 3; j++)
		{
			vtot[j] += mass*atoms[i].velocity[j];
		}
	}

	for ( j = 0; j < 3; j++)
		vtot[j] /= totMass;

	for ( i = 0; i < nAtoms; i++)
		for ( j = 0; j < 3; j++)
			atoms[i].velocity[j] -= vtot[j];
	/* check */
	/*
	for ( j = 0; j < 3; j++)
		vtot[j] = 0;
	for ( i = 0; i < nAtoms; i++)
	{
		mass = molecule.atoms[i].mass;
		for ( j = 0; j < 3; j++)
		{
			vtot[j] += mass*molecule.atoms[i].velocity[j];
		}
	}
	printf("Trans velocity = %f %f %f\n",vtot[0], vtot[1], vtot[2]);
	*/
}
/*********************************************************************************/
static void removeRotation(Molecule* molecule)
{
	double vtot[3] = {0,0,0};
	double cm[3] = {0,0,0};
	double L[3] = {0,0,0};
	int i;
	int j;
	int k;
	double mass = 1.0;
	double totMass = 0.0;
	double cdel[3];
	double vAng[3]={0,0,0};
	double tensor[3][3];
	double invTensor[3][3];
        double xx, xy,xz,yy,yz,zz;
	/* find the center of mass coordinates  and total velocity*/
	Atom* atoms = molecule->atoms;
	int nAtoms = molecule->nAtoms;


	for ( i = 0; i < nAtoms; i++)
	{
		mass = atoms[i].mass;
		totMass += mass;
		for ( j = 0; j < 3; j++)
			cm[j] += mass*atoms[i].coordinates[j];
		for ( j = 0; j < 3; j++)
			vtot[j] += mass*atoms[i].velocity[j];
	}


	for ( j = 0; j < 3; j++)
		cm[j] /= totMass;
	for ( j = 0; j < 3; j++)
		vtot[j] /= totMass;

	/*   compute the angular momentum  */
	for ( i = 0; i < nAtoms; i++)
	{
		mass = atoms[i].mass;
		for ( j = 0; j < 3; j++)
			L[j] += (
				atoms[i].coordinates[(j+1)%3]*atoms[i].velocity[(j+2)%3]
			      - atoms[i].coordinates[(j+2)%3]*atoms[i].velocity[(j+1)%3]
			      )*mass;
	}
	for ( j = 0; j < 3; j++)
		L[j] -= (
			cm[(j+1)%3]*vtot[(j+2)%3]
		      - cm[(j+2)%3]*vtot[(j+1)%3]
			      )*totMass;

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
	for ( i = 0; i < nAtoms; i++)
	{
		mass = atoms[i].mass;
		for ( j = 0; j < 3; j++)
			cdel[j] = atoms[i].coordinates[j]-cm[j];
		xx +=  cdel[0]*cdel[0]*mass;
		xy +=  cdel[0]*cdel[1]*mass;
		xz +=  cdel[0]*cdel[2]*mass;
		yy +=  cdel[1]*cdel[1]*mass;
		yz +=  cdel[1]*cdel[2]*mass;
		zz +=  cdel[2]*cdel[2]*mass;
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
	if(nAtoms>1)
	{
		double U0[3];
		double U1[3];
		for ( j = 0; j < 3; j++)U0[j] = atoms[0].coordinates[j];
		for ( j = 0; j < 3; j++)U1[j] = atoms[1].coordinates[j];
		//printf("!!!!!!!!!!!I cannot invert the rotational Tensor : linear molecule!\n");
		computeAngularVelocitiesForALinearMolecule(U0, U1, tensor, L, vAng);
		//printf("Angular velocityi before rotation = %f %f %f\n",vAng[0], vAng[1], vAng[2]);
	}
	/*  eliminate any rotation about the system center of mass */
	for ( i = 0; i < nAtoms; i++)
	{
		for ( j = 0; j < 3; j++)
			cdel[j] = atoms[i].coordinates[j]-cm[j];
		for ( j = 0; j < 3; j++)
			atoms[i].velocity[j] += 
				cdel[(j+1)%3]*vAng[(j+2)%3]-
				cdel[(j+2)%3]*vAng[(j+1)%3];
	}

/* Check */
/*
	for ( j = 0; j < 3; j++) L[j] = 0;
	for ( i = 0; i < nAtoms; i++)
	{
		mass = atoms[i].mass;
		for ( j = 0; j < 3; j++)
			L[j] += (
				atoms[i].coordinates[(j+1)%3]*atoms[i].velocity[(j+2)%3]
			      - atoms[i].coordinates[(j+2)%3]*atoms[i].velocity[(j+1)%3]
			      )*mass;
	}
	for ( j = 0; j < 3; j++)
		L[j] -= (
			cm[(j+1)%3]*vtot[(j+2)%3]
		      - cm[(j+2)%3]*vtot[(j+1)%3]
			      )*totMass;

	if(linear) computeAngularVelocitiesForALinearMolecule(molecule, tensor, L, vAng);
	else
	{
		for ( j = 0; j < 3; j++)
		{
			vAng[j] = 0;
			for ( k = 0; k < 3; k++)
				vAng[j] += invTensor[j][k]*L[k];
		}
	}

	printf("Angular velocity = %f %f %f\n",vAng[0], vAng[1], vAng[2]);
*/

}
/*********************************************************************************/
static void removeTranslationAndRotation(Molecule* molecule)
{
	molecule->klass->removeTranslation(molecule);
	molecule->klass->removeRotation(molecule);
}
/*********************************************************************************/
static void removeTranslationCluster(Molecule** molecules, int nMols)
{
	double vtot[3] = {0,0,0};
	int im;
	int i;
	int j;
	double mass = 1.0;
	double totMass = 0.0;

	for ( im = 0; im < nMols; im++)
	{
		Atom* atoms = molecules[im]->atoms;
		int nAtoms = molecules[im]->nAtoms;
		for ( i = 0; i < nAtoms; i++)
		{
			mass = atoms[i].mass;
			totMass += mass;
			for ( j = 0; j < 3; j++)
			{
				vtot[j] += mass*atoms[i].velocity[j];
			}
		}
	}

	for ( j = 0; j < 3; j++)
		vtot[j] /= totMass;

	for ( im = 0; im < nMols; im++)
	{
		Atom* atoms = molecules[im]->atoms;
		int nAtoms = molecules[im]->nAtoms;
		for ( i = 0; i < nAtoms; i++)
			for ( j = 0; j < 3; j++)
				atoms[i].velocity[j] -= vtot[j];
	}
}
/*********************************************************************************/
static void removeRotationCluster(Molecule** molecules, int nMols)
{
	double vtot[3] = {0,0,0};
	double cm[3] = {0,0,0};
	double L[3] = {0,0,0};
	int im;
	int i;
	int j;
	int k;
	double mass = 1.0;
	double totMass = 0.0;
	double cdel[3];
	double vAng[3]={0,0,0};
	double tensor[3][3];
	double invTensor[3][3];
        double xx, xy,xz,yy,yz,zz;
	/* find the center of mass coordinates  and total velocity*/

	for ( im = 0; im < nMols; im++)
	{
		Atom* atoms = molecules[im]->atoms;
		int nAtoms = molecules[im]->nAtoms;
		for ( i = 0; i < nAtoms; i++)
		{
			mass = atoms[i].mass;
			totMass += mass;
			for ( j = 0; j < 3; j++)
				cm[j] += mass*atoms[i].coordinates[j];
			for ( j = 0; j < 3; j++)
				vtot[j] += mass*atoms[i].velocity[j];
		}
	}


	for ( j = 0; j < 3; j++) cm[j] /= totMass;
	for ( j = 0; j < 3; j++) vtot[j] /= totMass;

	/*   compute the angular momentum  */
	for ( im = 0; im < nMols; im++)
	{
		Atom* atoms = molecules[im]->atoms;
		int nAtoms = molecules[im]->nAtoms;
		for ( i = 0; i < nAtoms; i++)
		{
			mass = atoms[i].mass;
			for ( j = 0; j < 3; j++)
			L[j] += (
				atoms[i].coordinates[(j+1)%3]*atoms[i].velocity[(j+2)%3]
			      - atoms[i].coordinates[(j+2)%3]*atoms[i].velocity[(j+1)%3]
			      )*mass;
		}
	}
	for ( j = 0; j < 3; j++)
		L[j] -= (
			cm[(j+1)%3]*vtot[(j+2)%3]
		      - cm[(j+2)%3]*vtot[(j+1)%3]
			      )*totMass;

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
	for ( im = 0; im < nMols; im++)
	{
		Atom* atoms = molecules[im]->atoms;
		int nAtoms = molecules[im]->nAtoms;
		for ( i = 0; i < nAtoms; i++)
		{
			mass = atoms[i].mass;
			for ( j = 0; j < 3; j++)
				cdel[j] = atoms[i].coordinates[j]-cm[j];
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
	if(molecules[0]->nAtoms>1)
	{
		double U0[3];
		double U1[3];
		Atom* atoms = molecules[0]->atoms;
		for ( j = 0; j < 3; j++)U0[j] = atoms[0].coordinates[j];
		for ( j = 0; j < 3; j++)U1[j] = atoms[1].coordinates[j];
		//printf("!!!!!!!!!!!I cannot invert the rotational Tensor : linear molecule!\n");
		computeAngularVelocitiesForALinearMolecule(U0, U1, tensor, L, vAng);
		//printf("Angular velocityi before rotation = %f %f %f\n",vAng[0], vAng[1], vAng[2]);
	}
	/*  eliminate any rotation about the system center of mass */
	for ( im = 0; im < nMols; im++)
	{
		Atom* atoms = molecules[im]->atoms;
		int nAtoms = molecules[im]->nAtoms;
		for ( i = 0; i < nAtoms; i++)
		{
			for ( j = 0; j < 3; j++)
				cdel[j] = atoms[i].coordinates[j]-cm[j];
			for ( j = 0; j < 3; j++)
				atoms[i].velocity[j] += 
				cdel[(j+1)%3]*vAng[(j+2)%3]-
				cdel[(j+2)%3]*vAng[(j+1)%3];
		}
	}
}
/*********************************************************************************/
static void removeTranslationAndRotationCluster(Molecule** molecules, int nMols)
{
	if(nMols<1) return;
	molecules[0]->klass->removeTranslationCluster(molecules,nMols);
	molecules[0]->klass->removeRotationCluster(molecules,nMols);
}
/*********************************************************************************/
static void removeTranslationForceCluster(Molecule** molecules, int nMols, double** vectors)
{
	// Vectors = FORCE or Moment
	// atot = acc. or veloc. 
	double atot[3] = {0,0,0};
	int im;
	int i;
	int j;
	double tot = 0.0;

	for ( im = 0; im < nMols; im++)
	{
		Atom* atoms = molecules[im]->atoms;
		int nAtoms = molecules[im]->nAtoms;
		for ( i = 0; i < nAtoms; i++)
		{
			tot += atoms[i].mass;
			for ( j = 0; j < 3; j++)
			{
				atot[j] += vectors[im][3*i+j];
			}
		}
	}

	for ( j = 0; j < 3; j++)
		atot[j] /= tot;

	for ( im = 0; im < nMols; im++)
	{
		int nAtoms = molecules[im]->nAtoms;
		for ( i = 0; i < nAtoms; i++)
			for ( j = 0; j < 3; j++)
				vectors[im][3*i+j] -=   molecules[im]->atoms[i].mass*atot[j];
	}
}
/*********************************************************************************/
static void removeRotationForceCluster(Molecule** molecules, int nMols, double** vectors)
{
	double atot[3] = {0,0,0};
	double cm[3] = {0,0,0};
	double L[3] = {0,0,0};
	int im;
	int i;
	int j;
	int k;
	double mass = 1.0;
	double totMass = 0.0;
	double tot = 0.0;
	double cdel[3];
	double vAng[3]={0,0,0};
	double tensor[3][3];
	double invTensor[3][3];
        double xx, xy,xz,yy,yz,zz;
	/* find the center of mass coordinates  and total vectors*/

	for ( im = 0; im < nMols; im++)
	{
		Atom* atoms = molecules[im]->atoms;
		int nAtoms = molecules[im]->nAtoms;
		for ( i = 0; i < nAtoms; i++)
		{
			mass = atoms[i].mass;
			totMass += mass;
			tot += 1;
			for ( j = 0; j < 3; j++)
				cm[j] += mass*atoms[i].coordinates[j];
			for ( j = 0; j < 3; j++)
				atot[j] += vectors[im][3*i+j];
		}
	}


	for ( j = 0; j < 3; j++) cm[j] /= totMass;
	//for ( j = 0; j < 3; j++) atot[j] /= totMass;

	/*   compute the angular momentum  */
	for ( im = 0; im < nMols; im++)
	{
		Atom* atoms = molecules[im]->atoms;
		int nAtoms = molecules[im]->nAtoms;
		for ( i = 0; i < nAtoms; i++)
		{
			mass = atoms[i].mass;
			for ( j = 0; j < 3; j++)
			L[j] += (
				atoms[i].coordinates[(j+1)%3]*vectors[im][3*i+(j+2)%3]
			      - atoms[i].coordinates[(j+2)%3]*vectors[im][3*i+(j+1)%3]
			      );
		}
	}
	for ( j = 0; j < 3; j++)
		L[j] -= (
			cm[(j+1)%3]*atot[(j+2)%3]
		      - cm[(j+2)%3]*atot[(j+1)%3]
			      );
			      //)*totMass;

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
	for ( im = 0; im < nMols; im++)
	{
		Atom* atoms = molecules[im]->atoms;
		int nAtoms = molecules[im]->nAtoms;
		for ( i = 0; i < nAtoms; i++)
		{
			mass = atoms[i].mass;
			for ( j = 0; j < 3; j++)
				cdel[j] = atoms[i].coordinates[j]-cm[j];
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
	if(molecules[0]->nAtoms>1)
	{
		double U0[3];
		double U1[3];
		Atom* atoms = molecules[0]->atoms;
		for ( j = 0; j < 3; j++)U0[j] = atoms[0].coordinates[j];
		for ( j = 0; j < 3; j++)U1[j] = atoms[1].coordinates[j];
		//printf("!!!!!!!!!!!I cannot invert the rotational Tensor : linear molecule!\n");
		computeAngularVelocitiesForALinearMolecule(U0,U1, tensor, L, vAng);
		//printf("Angular vectors before rotation = %f %f %f\n",vAng[0], vAng[1], vAng[2]);
	}
	/*  eliminate any rotation about the system center of mass */
	for ( im = 0; im < nMols; im++)
	{
		Atom* atoms = molecules[im]->atoms;
		int nAtoms = molecules[im]->nAtoms;
		for ( i = 0; i < nAtoms; i++)
		{
			for ( j = 0; j < 3; j++)
				cdel[j] = atoms[i].coordinates[j]-cm[j];
			for ( j = 0; j < 3; j++)
				vectors[im][3*i+j] += 
				(cdel[(j+1)%3]*vAng[(j+2)%3]-
				cdel[(j+2)%3]*vAng[(j+1)%3])
				*atoms[i].mass
				;
		}
	}
}
/*********************************************************************************/
static void removeTranslationAndRotationForceCluster(Molecule** molecules, int nMols, double** vectors)
{
	if(nMols<1) return;
	molecules[0]->klass->removeTranslationForceCluster(molecules, nMols, vectors);
	molecules[0]->klass->removeRotationForceCluster(molecules, nMols, vectors);
}
/*********************************************************************************/
static void removeTranslationForce(Molecule* molecule, double* f)
{
	// Vectors = FORCE or Moment
	// atot = acc. or veloc. 
	double atot[3] = {0,0,0};
	int i;
	int j;
	double tot = 0.0;

	{
		Atom* atoms = molecule->atoms;
		int nAtoms = molecule->nAtoms;
		for ( i = 0; i < nAtoms; i++)
		{
			tot += atoms[i].mass;
			for ( j = 0; j < 3; j++)
			{
				atot[j] += f[3*i+j];
			}
		}
	}

	for ( j = 0; j < 3; j++)
		atot[j] /= tot;

	{
		int nAtoms = molecule->nAtoms;
		for ( i = 0; i < nAtoms; i++)
			for ( j = 0; j < 3; j++)
				f[3*i+j] -=   molecule->atoms[i].mass*atot[j];
	}
}
/*********************************************************************************/
static void removeRotationForce(Molecule* molecule, double* f)
{
	double atot[3] = {0,0,0};
	double cm[3] = {0,0,0};
	double L[3] = {0,0,0};
	int i;
	int j;
	int k;
	double mass = 1.0;
	double totMass = 0.0;
	double tot = 0.0;
	double cdel[3];
	double vAng[3]={0,0,0};
	double tensor[3][3];
	double invTensor[3][3];
        double xx, xy,xz,yy,yz,zz;
	/* find the center of mass coordinates  and total f*/

	{
		Atom* atoms = molecule->atoms;
		int nAtoms = molecule->nAtoms;
		for ( i = 0; i < nAtoms; i++)
		{
			mass = atoms[i].mass;
			totMass += mass;
			tot += 1;
			for ( j = 0; j < 3; j++)
				cm[j] += mass*atoms[i].coordinates[j];
			for ( j = 0; j < 3; j++)
				atot[j] += f[3*i+j];
		}
	}


	for ( j = 0; j < 3; j++) cm[j] /= totMass;

	/*   compute the angular momentum  */
	{
		Atom* atoms = molecule->atoms;
		int nAtoms = molecule->nAtoms;
		for ( i = 0; i < nAtoms; i++)
		{
			mass = atoms[i].mass;
			for ( j = 0; j < 3; j++)
			L[j] += (
				atoms[i].coordinates[(j+1)%3]*f[3*i+(j+2)%3]
			      - atoms[i].coordinates[(j+2)%3]*f[3*i+(j+1)%3]
			      );
		}
	}
	for ( j = 0; j < 3; j++)
		L[j] -= (
			cm[(j+1)%3]*atot[(j+2)%3]
		      - cm[(j+2)%3]*atot[(j+1)%3]
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
	{
		Atom* atoms = molecule->atoms;
		int nAtoms = molecule->nAtoms;
		for ( i = 0; i < nAtoms; i++)
		{
			mass = atoms[i].mass;
			for ( j = 0; j < 3; j++)
				cdel[j] = atoms[i].coordinates[j]-cm[j];
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
	if(molecule->nAtoms>1)
	{
		double U0[3];
		double U1[3];
		Atom* atoms = molecule->atoms;
		for ( j = 0; j < 3; j++)U0[j] = atoms[0].coordinates[j];
		for ( j = 0; j < 3; j++)U1[j] = atoms[1].coordinates[j];
		//printf("!!!!!!!!!!!I cannot invert the rotational Tensor : linear molecule!\n");
		computeAngularVelocitiesForALinearMolecule(U0,U1, tensor, L, vAng);
		//printf("Angular f before rotation = %f %f %f\n",vAng[0], vAng[1], vAng[2]);
	}
	/*  eliminate any rotation about the system center of mass */
	{
		Atom* atoms = molecule->atoms;
		int nAtoms = molecule->nAtoms;
		for ( i = 0; i < nAtoms; i++)
		{
			for ( j = 0; j < 3; j++)
				cdel[j] = atoms[i].coordinates[j]-cm[j];
			for ( j = 0; j < 3; j++)
				f[3*i+j] += 
				(cdel[(j+1)%3]*vAng[(j+2)%3]-
				cdel[(j+2)%3]*vAng[(j+1)%3])
				*atoms[i].mass
				;
		}
	}
}
/*********************************************************************************/
static void removeTranslationAndRotationForce(Molecule* molecule, double*f)
{
	molecule->klass->removeTranslationForce(molecule, f);
	molecule->klass->removeRotationForce(molecule, f);
}
/*********************************************************************************/
static void removeTranslationMoments(Molecule** molecules, int nMols, double*** P)
{
	double vtot[3] = {0,0,0};
	int im;
	int i;
	int j;
	double totMass = 0.0;

	for ( im = 0; im < nMols; im++)
	{
		int nAtoms = molecules[im]->nAtoms;
		for ( i = 0; i < nAtoms; i++)
		{
			totMass += molecules[im]->atoms[i].mass;
			for ( j = 0; j < 3; j++)
			{
				vtot[j] += P[j][im][i];
			}
		}
	}
	for ( j = 0; j < 3; j++) vtot[j] /= totMass;


	for ( im = 0; im < nMols; im++)
	{
		int nAtoms = molecules[im]->nAtoms;
		for ( i = 0; i < nAtoms; i++)
			for ( j = 0; j < 3; j++)
				P[j][im][i] -= molecules[im]->atoms[i].mass*vtot[j];
	}
}
/*********************************************************************************/
static void removeRotationMoments(Molecule** molecules, int nMols, double*** P)
{
	double ptot[3] = {0,0,0};
	double cm[3] = {0,0,0};
	double L[3] = {0,0,0};
	int im;
	int i;
	int j;
	int k;
	double mass = 1.0;
	double totMass = 0.0;
	double cdel[3];
	double pAng[3]={0,0,0};
	double tensor[3][3];
	double invTensor[3][3];
        double xx, xy,xz,yy,yz,zz;
	boolean linear = FALSE;
	double n = 0;
	/* find the center of mass coordinates  and total vectors*/

	for ( im = 0; im < nMols; im++)
	{
		Atom* atoms = molecules[im]->atoms;
		int nAtoms = molecules[im]->nAtoms;
		for ( i = 0; i < nAtoms; i++)
		{
			mass = atoms[i].mass;
			totMass += mass;
			n++;
			for ( j = 0; j < 3; j++)
				cm[j] += mass*atoms[i].coordinates[j];
			for ( j = 0; j < 3; j++)
				ptot[j] += P[j][im][i];
		}
	}


	for ( j = 0; j < 3; j++) cm[j] /= totMass;
	//for ( j = 0; j < 3; j++) ptot[j] /= n;

	/*   compute the angular momentum  */
	for ( im = 0; im < nMols; im++)
	{
		Atom* atoms = molecules[im]->atoms;
		int nAtoms = molecules[im]->nAtoms;
		for ( i = 0; i < nAtoms; i++)
		{
			//mass = atoms[i].mass;
			for ( j = 0; j < 3; j++)
			L[j] += (
				atoms[i].coordinates[(j+1)%3]*P[(j+2)%3][im][i]
			      - atoms[i].coordinates[(j+2)%3]*P[(j+1)%3][im][i]
			      //)*mass;
			      );
		}
	}
	for ( j = 0; j < 3; j++)
		L[j] -= (
			cm[(j+1)%3]*ptot[(j+2)%3]
		      - cm[(j+2)%3]*ptot[(j+1)%3]
			      );
			      //)*totMass;

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
	for ( im = 0; im < nMols; im++)
	{
		Atom* atoms = molecules[im]->atoms;
		int nAtoms = molecules[im]->nAtoms;
		for ( i = 0; i < nAtoms; i++)
		{
			mass = atoms[i].mass;
			for ( j = 0; j < 3; j++)
				cdel[j] = atoms[i].coordinates[j]-cm[j];

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
			pAng[j] = 0;
			for ( k = 0; k < 3; k++)
				pAng[j] += invTensor[j][k]*L[k];
		}
		printf("Angular vectors before rotation = %f %f %f\n",pAng[0], pAng[1], pAng[2]);
	}
	else
	if(molecules[0]->nAtoms>1)
	{
		double U0[3];
		double U1[3];
		Atom* atoms = molecules[0]->atoms;
		for ( j = 0; j < 3; j++)U0[j] = atoms[0].coordinates[j];
		for ( j = 0; j < 3; j++)U1[j] = atoms[1].coordinates[j];
		//printf("!!!!!!!!!!!I cannot invert the rotational Tensor : linear molecule!\n");
		computeAngularVelocitiesForALinearMolecule(U0, U1, tensor, L, pAng);
		printf("Angular vectors before rotation = %f %f %f\n",pAng[0], pAng[1], pAng[2]);
		linear = TRUE;
	}
	/*  eliminate any rotation about the system center of mass */
	for ( im = 0; im < nMols; im++)
	{
		Atom* atoms = molecules[im]->atoms;
		int nAtoms = molecules[im]->nAtoms;
		for ( i = 0; i < nAtoms; i++)
		{
			for ( j = 0; j < 3; j++)
				cdel[j] = atoms[i].coordinates[j]-cm[j];
			for ( j = 0; j < 3; j++)
				P[j][im][i] += 
				(cdel[(j+1)%3]*pAng[(j+2)%3]-
				cdel[(j+2)%3]*pAng[(j+1)%3])*atoms[i].mass;
		}
	}
/* Check */
	for ( j = 0; j < 3; j++) L[j] = 0;
	for ( im = 0; im < nMols; im++)
	{
		Atom* atoms = molecules[im]->atoms;
		int nAtoms = molecules[im]->nAtoms;
		for ( i = 0; i < nAtoms; i++)
		{
			//mass = atoms[i].mass;
			for ( j = 0; j < 3; j++)
			L[j] += (
				atoms[i].coordinates[(j+1)%3]*P[(j+2)%3][im][i]
			      - atoms[i].coordinates[(j+2)%3]*P[(j+1)%3][im][i]
			      //)*mass;
			      );
		}
	}
	for ( j = 0; j < 3; j++)
		L[j] -= (
			cm[(j+1)%3]*ptot[(j+2)%3]
		      - cm[(j+2)%3]*ptot[(j+1)%3]
			      );
			      //)*totMass;

	if(linear) 
	{
		double U0[3];
		double U1[3];
		Atom* atoms = molecules[0]->atoms;
		for ( j = 0; j < 3; j++) U0[j] = atoms[0].coordinates[j];
		for ( j = 0; j < 3; j++) U1[j] = atoms[1].coordinates[j];
		computeAngularVelocitiesForALinearMolecule(U0, U1, tensor, L, pAng);
		printf("Linear molecule Angular velocity = %f %f %f\n",pAng[0], pAng[1], pAng[2]);
	}
	else
	{
		for ( j = 0; j < 3; j++)
		{
			pAng[j] = 0;
			for ( k = 0; k < 3; k++)
				pAng[j] += invTensor[j][k]*L[k];
		}
		printf("Angular velocity = %f %f %f\n",pAng[0], pAng[1], pAng[2]);
	}

}
/*********************************************************************************/
static void removeTranslationAndRotationMoments(Molecule** molecules, int nMols, double*** P)
{
	if(nMols<1) return;
	molecules[0]->klass->removeTranslationMoments(molecules, nMols, P);
	molecules[0]->klass->removeRotationMoments(molecules, nMols, P);
}
/*********************************************************************************/
/* vectors = one vectors of 3*nAtoms elements */
static void removeRotationAcceleration(Molecule* molecule, double* a)
{
	int nAtoms = molecule->nAtoms;
	Atom* atoms = molecule->atoms;
	double cm[3] = {0,0,0};
	double L[3] = {0,0,0};
	int i;
	int j;
	int k;
	double mass = 1.0;
	double totMass = 0.0;
	double cdel[3];
	double fAng[3]={0,0,0};
	double tensor[3][3];
	double invTensor[3][3];
        double xx, xy,xz,yy,yz,zz;
	boolean linear = FALSE;
	double vectot[3] = {0,0,0};

	for ( i = 0; i < nAtoms; i++)
	{
		mass = atoms[i].mass;
		totMass += mass;
		for ( j = 0; j < 3; j++) vectot[j] += mass*a[3*i+j];
		for ( j = 0; j < 3; j++) cm[j] += mass*atoms[i].coordinates[j];
	}

	for ( j = 0; j < 3; j++) cm[j] /= totMass;
	for ( j = 0; j < 3; j++) vectot[j] /= totMass;

	/*   compute the angular momentum  */
	for ( i = 0; i < nAtoms; i++)
	{
		mass = atoms[i].mass;
		for ( j = 0; j < 3; j++)
			L[j] += (
				atoms[i].coordinates[(j+1)%3]*a[3*i+(j+2)%3]
			      - atoms[i].coordinates[(j+2)%3]*a[3*i+(j+1)%3]
			      )*mass;
	}
	for ( j = 0; j < 3; j++)
		L[j] -= (
			cm[(j+1)%3]*vectot[(j+2)%3]
		      - cm[(j+2)%3]*vectot[(j+1)%3]
			      )*totMass;

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
	for ( i = 0; i < nAtoms; i++)
	{
		mass = atoms[i].mass;
		for ( j = 0; j < 3; j++) cdel[j] = atoms[i].coordinates[j]-cm[j];
		xx +=  cdel[0]*cdel[0]*mass;
		xy +=  cdel[0]*cdel[1]*mass;
		xz +=  cdel[0]*cdel[2]*mass;
		yy +=  cdel[1]*cdel[1]*mass;
		yz +=  cdel[1]*cdel[2]*mass;
		zz +=  cdel[2]*cdel[2]*mass;
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
			fAng[j] = 0;
			for ( k = 0; k < 3; k++)
				fAng[j] += invTensor[j][k]*L[k];
		}
	}
	else
	if(molecule->nAtoms>1)
	{
		double U0[3];
		double U1[3];
		Atom* atoms = molecule->atoms;
		for ( j = 0; j < 3; j++)U0[j] = atoms[0].coordinates[j];
		for ( j = 0; j < 3; j++)U1[j] = atoms[1].coordinates[j];
		//printf("!!!!!!!!!!!I cannot invert the rotational Tensor : linear molecule!\n");
		computeAngularVelocitiesForALinearMolecule(U0, U1, tensor, L, fAng);
		//printf("Angular forces before rotation = %f %f %f\n",fAng[0], fAng[1], fAng[2]);
		linear = TRUE;
	}
	/*  eliminate any rotation about the system center of mass */
	for ( i = 0; i < nAtoms; i++)
	{
		for ( j = 0; j < 3; j++) cdel[j] = atoms[i].coordinates[j]-cm[j];
		for ( j = 0; j < 3; j++)
			a[3*i+j] += 
				(cdel[(j+1)%3]*fAng[(j+2)%3]-
				cdel[(j+2)%3]*fAng[(j+1)%3])
				;
	}

/* Check */
	for ( j = 0; j < 3; j++) L[j] = 0;
	for ( i = 0; i < nAtoms; i++)
	{
		mass = atoms[i].mass;
		for ( j = 0; j < 3; j++)
			L[j] += (
				atoms[i].coordinates[(j+1)%3]*a[3*i+(j+2)%3]
			      - atoms[i].coordinates[(j+2)%3]*a[3*i+(j+1)%3]
			      )*mass;
	}
	for ( j = 0; j < 3; j++)
		L[j] -= (
			cm[(j+1)%3]*vectot[(j+2)%3]
		      - cm[(j+2)%3]*vectot[(j+1)%3]
			      )*totMass;
	if(linear)
	{
		double U0[3];
		double U1[3];
		Atom* atoms = molecule->atoms;
		for ( j = 0; j < 3; j++)U0[j] = atoms[0].coordinates[j];
		for ( j = 0; j < 3; j++)U1[j] = atoms[1].coordinates[j];
		computeAngularVelocitiesForALinearMolecule(U0,U1, tensor, L, fAng);
	}
	else
	{
		for ( j = 0; j < 3; j++)
		{
			fAng[j] = 0;
			for ( k = 0; k < 3; k++)
				fAng[j] += invTensor[j][k]*L[k];
		}
	}

	printf("Angular forces = %f %f %f\n",fAng[0], fAng[1], fAng[2]);
}
/*********************************************************************************/
/* a = one vectors of 3*nAtoms elements */
static void removeTranslationAcceleration(Molecule* molecule, double* a)
{
	int nAtoms = molecule->nAtoms;
	double vectot[3] = {0,0,0};
	int i,j;
	Atom* atoms = molecule->atoms;
	double totMass = 0;
	for ( i = 0; i < nAtoms; i++)
	{
		for ( j = 0; j < 3; j++) vectot[j] += atoms[i].mass*a[3*i+j];
		totMass += atoms[i].mass;
	}

	for ( j = 0; j < 3; j++) vectot[j] /= totMass;
	for ( i = 0; i < nAtoms; i++)
		for ( j = 0; j < 3; j++) a[3*i+j] -= vectot[j];
}
/*********************************************************************************/
/* vectors = one vectors of 3*nAtoms elements */
static void removeTranslationAndRotationAcceleration(Molecule* molecule, double* a)
{
	molecule->klass->removeTranslationAcceleration(molecule,a);
	molecule->klass->removeRotationAcceleration(molecule,a);
}
/*****************************************************************************/
static boolean printMolecule(Molecule* molecule, FILE* file)
{
	int j;
	if(!molecule) return FALSE;
	if(!file) return FALSE;

	{
		int k,nc;
		fprintf(file,"%f\n",molecule->potentialEnergy);
		fprintf(file,"0\n");
		fprintf(file,"%d %d %d\n",molecule->nAtoms, molecule->totalCharge, molecule->spinMultiplicity);
		for(j=0;j<molecule->nAtoms;j++)
		{
			nc = 0;
			for(k=0;k<molecule->nAtoms;k++) if(molecule->atoms[j].typeConnections[k]>0) nc++;
			fprintf(file," %s %s %s %s %d %f %d %d %f %f %f %d ", 
				molecule->atoms[j].prop.symbol,
				molecule->atoms[j].mmType,
				molecule->atoms[j].pdbType,
				molecule->atoms[j].residueName,
				molecule->atoms[j].residueNumber,
				molecule->atoms[j].charge,
				molecule->atoms[j].layer,
				molecule->atoms[j].variable,
				molecule->atoms[j].coordinates[0],
				molecule->atoms[j].coordinates[1],
				molecule->atoms[j].coordinates[2],
				nc
				);
			for(k=0;k<molecule->nAtoms;k++) 
			{
		 		int nk = molecule->atoms[k].N-1;
				if(molecule->atoms[j].typeConnections[nk]>0) 
				fprintf(file," %d %d", nk+1,molecule->atoms[j].typeConnections[nk]);
			}
			fprintf(file,"\n");
		}
	}
	return TRUE;

}
/********************************************************************************/
static double getKelvin(Molecule* molecule)
{
	int nFree = molecule->nFree;
	if(nFree<1) return 0;
	double kin = molecule->klass->getKineticEnergy(molecule);
	return 2*kin / ( nFree * Kb);
}
/*****************************************************************************/
static void scaleVelocities(Molecule* molecule, double temperature)
{
	double kelvin = molecule->klass->getKelvin(molecule);
	double scale = 1.0;
	int i,j;
	if(temperature<=0) return;
	if(kelvin<=0) return;

	scale = sqrt(temperature/kelvin);
#ifdef DEBUG
	printf("temp = %f kelvin = %f scale = %f\n",temperature, kelvin, scale);
#endif
	for(i = 0;i<molecule->nAtoms; i++) 
		if(molecule->atoms[i].variable)
		for ( j = 0; j < 3; j++)
			molecule->atoms[i].velocity[j] *= scale;
}
/*****************************************************************************/
static void setMaxwellVelocities(Molecule* molecule, double temperature)
{
	int i,j;
	for ( i = 0; i < molecule->nAtoms; i++)
	{
        	if(!molecule->atoms[i].variable) 
		for ( j = 0; j < 3; j++) molecule->atoms[i].velocity[j] = 0.0;
		else
		{
			double speed = maxwel(molecule->atoms[i].mass,temperature);
			getRandVect(speed, molecule->atoms[i].velocity);
		}
	}
	molecule->klass->scaleVelocities(molecule, temperature);
}
/*****************************************************************************/
static boolean setMaxwellVelocitiesIfNull(Molecule* molecule, double temperature)
{
	
	double ekin = molecule->klass->getKineticEnergy(molecule);
	if(fabs(ekin)>1e-14) return FALSE;
	setMaxwellVelocities(molecule, temperature);
	return TRUE;
}
/*****************************************************************************/
static boolean resetConstraints(Molecule* molecule, Constraints constraints)
{
	int i;
	int nvariables = 0;
        int nAtoms = molecule->nAtoms;

        molecule->constraints = constraints;
	setRattleConstraintsParameters(molecule);

        for ( i = 0; i < nAtoms; i++)
                	if(molecule->atoms[i].variable) nvariables +=1;
        if(nvariables==0) 
	{
		nvariables = nAtoms;
        	for ( i = 0; i < nAtoms; i++)
                	molecule->atoms[i].variable = TRUE;
	}
        molecule->nFree = 3* nvariables;
        molecule->nFree -= molecule->numberOfRattleConstraintsTerms;

        if(nvariables==nAtoms) molecule->nFree -=6;
        if(nvariables==2 && 2==nAtoms) molecule->nFree +=1;
        if(nvariables==nAtoms-1) molecule->nFree -=3;
        if(nvariables==nAtoms-2) molecule->nFree -=1;
	if( molecule->nFree<1)  {
		printf("nFree =%d < 1\n",molecule->nFree);
		exit(1);
	}
	printf("#Number of free coordinates =%d \n",molecule->nFree);
	return TRUE;
}
/**********************************************************************/
static void setRattleConstraintsParameters(Molecule* m)
{
	int i;
	int j;
	int k;
	int a1,a2,a3;
	double r2;
	double d;
	int numberOfRattleConstraintsTerms = 0;
	double* rattleConstraintsTerms[RATTLEDIM];

	m->numberOfRattleConstraintsTerms = 0;
	for( i=0; i<RATTLEDIM;i++) m->rattleConstraintsTerms[i] = NULL;

	if(m->nAtoms<1) return;

	if(m->constraints==NOCONSTRAINTS) return;
	numberOfRattleConstraintsTerms = m->numberOf2Connections;
	if(m->constraints==BONDSANGLESCONSTRAINTS) 
		numberOfRattleConstraintsTerms += m->numberOf3Connections;

	if(numberOfRattleConstraintsTerms<1) return;
	for( i=0; i<RATTLEDIM;i++)
       		rattleConstraintsTerms[i] = malloc(numberOfRattleConstraintsTerms*sizeof(double));


	/* 1=a1, 2=a2, 3=r2a1a2 */
	/* RATTLEDIM 	3 */
	j = 0;
	for ( i = 0; i < m->numberOf2Connections; i++)
	{
		a1 = m->connected2[0][i];
		a2 = m->connected2[1][i];
		if(!m->atoms[a1].variable &&!m->atoms[a2].variable) continue;
		r2 = 0;
		for (k=0;k<3;k++)
		{
			d = m->atoms[a1].coordinates[k]-m->atoms[a2].coordinates[k];
			r2 +=d*d;
		}
		rattleConstraintsTerms[0][j] = a1;
		rattleConstraintsTerms[1][j] = a2;
		rattleConstraintsTerms[2][j] = r2;
		j++;
	}
	if(m->constraints==BONDSANGLESCONSTRAINTS)
	{
		int a1p, a2p;
		int* nConnections = NULL;
		int* nAngles = NULL;
       		nConnections = malloc(m->nAtoms*sizeof(int));
       		nAngles = malloc(m->nAtoms*sizeof(int));
		for ( i = 0; i < m->nAtoms; i++)
		{
			nConnections[i] = 0;
			nAngles[i] = 0;
		}
		for ( i = 0; i < m->nAtoms; i++)
		if(m->atoms[i].typeConnections)
		{
			for ( k = 0; k < m->nAtoms; k++)
				if(i!=k && m->atoms[i].typeConnections[m->atoms[k].N-1]>0) nConnections[i]++;
			/* printf("%d %s nCon=%d\n",i,m->atoms[i].mmType,nConnections[i]);*/
		}
		for ( i = 0; i < m->numberOf3Connections; i++)
		{
			a1 = m->connected3[0][i];
			a2 = m->connected3[1][i];
			a3 = m->connected3[2][i];
			if(!m->atoms[a1].variable &&!m->atoms[a3].variable) continue;
			if(nAngles[a2]>=2*nConnections[a2]-3) continue;
			for (k=0;k<j;k++)
			{
				a1p = (int)rattleConstraintsTerms[0][k];
				a2p = (int)rattleConstraintsTerms[1][k];
				if(a1p==a1 && a2p==a3) break;
				if(a1p==a3 && a2p==a1) break;
			}
			if(k!=j) continue;

			nAngles[a2]++;
			r2 = 0;
			for (k=0;k<3;k++)
			{
				d = m->atoms[a1].coordinates[k]-m->atoms[a3].coordinates[k];
				r2 +=d*d;
			}
			rattleConstraintsTerms[0][j] = a1;
			rattleConstraintsTerms[1][j] = a3;
			rattleConstraintsTerms[2][j] = r2;
			j++;
		}
		/*
		for ( i = 0; i < m->nAtoms; i++)
		{
			printf("%d %s nAngle = %d 2*nCon-3=%d\n",i,m->atoms[i].mmType,nAngles[i],2*nConnections[i]-3);
		}
		*/
       		if(nConnections) free(nConnections);
       		if(nAngles) free(nAngles);
	}

	if(j<1)
	{
		numberOfRattleConstraintsTerms=0;
		for( i=0; i<RATTLEDIM;i++)
		{
       			free(rattleConstraintsTerms[i]);
       			rattleConstraintsTerms[i] = NULL;
		}
	}
	else if(numberOfRattleConstraintsTerms!=j)
	{
		numberOfRattleConstraintsTerms=j;
		for( i=0; i<RATTLEDIM;i++)
		{
       			rattleConstraintsTerms[i] = 
				realloc(rattleConstraintsTerms[i],numberOfRattleConstraintsTerms*sizeof(double));
		}
	}
	m->numberOfRattleConstraintsTerms = numberOfRattleConstraintsTerms;
	for( i=0; i<RATTLEDIM;i++)
       		m->rattleConstraintsTerms[i] = rattleConstraintsTerms[i]; 

/*
	for ( i = 0; i < m->numberOfRattleConstraintsTerms; i++)
	{
			a1 = (int)rattleConstraintsTerms[0][i];
			a2 = (int)rattleConstraintsTerms[1][i];
			r2 = rattleConstraintsTerms[2][i];
			printf("%d  %d %s %s r2= %f\n",
				a1,a2,
				m->atoms[a1].mmType,
				m->atoms[a2].mmType,
				r2);
	}
*/
}
/*****************************************************************************/
static int getRanInt(int* index, int n, int M)
{
	boolean accepted=TRUE;
	int ir = 0;
	int itmax = 1000; 
	int it = 0;
	do
	{
		int i;
		it++;
		ir = rand()%M;
		accepted = TRUE;
		for(i=0;i<n;i++) if(ir==index[i]) {accepted = FALSE; break;}
	}while (!accepted && it<itmax);
	if(it>=itmax)
	{
		fprintf(stderr,"Error : number of call for rand() > %d\n",itmax);
		fprintf(stderr,"Program stopped in getRanInt of Molecule file\n");
		exit(1);
	}
	return ir;
}
/*****************************************************************************/
static boolean setRandomFragments(Molecule* molecule)
{
	int n = 0;
	double bd = 0;
	int itmax = 1000; 
	int nRes=0;
	int** listByResidue = NULL;
	int* nByResidue = NULL;
	int nAtoms = molecule->nAtoms;
	Molecule* mol0 = NULL;
	int i;

	if(nAtoms<=1) return FALSE;

	nRes=0;
	for (i = 0; i<nAtoms;i++) if(molecule->atoms[i].residueNumber>nRes) nRes++;
	nRes++;
	if(nRes<2) return FALSE;

	mol0 = molecule->klass->copy(molecule);

	nByResidue = malloc(nRes*sizeof(int));
	listByResidue = malloc(nRes*sizeof(int*));
	for (n = 0; n<nRes;n++) listByResidue[n] = malloc(nAtoms*sizeof(int));

	for (n = 0; n<nRes;n++) nByResidue[n]=0;

	for (n = 0; n<nRes;n++) 
	for (i = 0; i<nAtoms;i++) 
	{
		if(molecule->atoms[i].residueNumber==n) 
		{
			int k = nByResidue[n];
			listByResidue[n][k] =  i;
			nByResidue[n]++;
		}
	}
	// don't change n=0
	for (n = 1; n<nRes;n++)
	{
		boolean accepted;
		int it=0;
		do{
			double direction[3];
			double v[3];
			double O[3];
			double a1a2[3];
			double distance = 0;
			double MRot[3][3];
			int i,j,k;
			Atom* a1 = NULL;
			Atom* a2 = NULL;
			int k1=rand()%nByResidue[0];
			int i1=listByResidue[0][k1];
			int k2=rand()%nByResidue[n];
			int i2=listByResidue[n][k2];
			a1 = &mol0->atoms[i1];
			a2 = &mol0->atoms[i2];
			distance = 0;
			for (k=0;k<3;k++) { double dij = a1->coordinates[k]-a2->coordinates[k]; distance +=dij*dij;}
			for (k=0;k<3;k++) a1a2[k] = a2->coordinates[k]-a1->coordinates[k]; 
			distance = sqrt(distance);
			getRandDirection(direction);
			//printf("distance=%f\n",distance);
			for (k=0;k<3;k++) direction[k]*= distance;

			//double mod=0; for (k=0;k<3;k++) mod+= direction[k]*direction[k]; mod=sqrt(mod); printf("mod=%f\n",mod);

			for (k=0;k<3;k++) v[k] = direction[k]-a1a2[k];

			for(j=0;j<nByResidue[n];j++)
			{
				int jj=listByResidue[n][j];
				for (k=0;k<3;k++) molecule->atoms[jj].coordinates[k] = mol0->atoms[jj].coordinates[k] + v[k];
			}
			get3DRandMatrix(MRot);
			
			for (k=0;k<3;k++) O[k] = molecule->atoms[i2].coordinates[k];
			for(j=0;j<nByResidue[n];j++)
			{
				int l;
				int jj=listByResidue[n][j];
				for (k=0;k<3;k++) v[k] = molecule->atoms[jj].coordinates[k]-O[k];
				for (l=0;l<3;l++) 
				{
					molecule->atoms[jj].coordinates[l] = O[l];
					for (k=0;k<3;k++) molecule->atoms[jj].coordinates[l] += v[k]*MRot[l][k];
				}
			}

			accepted = TRUE;
			for(j=0;j<nByResidue[n] && accepted;j++)
			for(i=0;i<nAtoms && accepted;i++)
			{
				if(molecule->atoms[i].residueNumber>=n) continue;
				int ii=i;
				int jj=listByResidue[n][j];
				a1 = &molecule->atoms[ii];
				a2 = &molecule->atoms[jj];

				double distance = 0;
				for (k=0;k<3;k++) { double dij = a1->coordinates[k]-a2->coordinates[k]; distance +=dij*dij;}
				distance = sqrt(distance)/BOHRTOANG;
				bd = a1->prop.covalentRadii + a2->prop.covalentRadii;
				if(distance<0.6*bd)
				{
					accepted = FALSE;
					break;
				}
			}
			it++;
		}while(!accepted && it<itmax);
	}
	if(nByResidue) free(nByResidue);
	for (n = 0; n<nRes;n++) if(listByResidue[n]) free(listByResidue[n]);
	if(listByResidue) free(listByResidue);

	mol0->klass->free(mol0);
	//setConnections(molecule);
	return TRUE;
}
/*****************************************************************************/
static boolean setRandomPositions(Molecule* molecule)
{
	int* index = NULL;
	int n = 0;
	int ir;
	double bd = 0;
	int ic = 0;
	int it = 0;
	int itmax = 1000; 
	int i;
	int j;

	if(molecule->nAtoms<=1) return FALSE;
	index = malloc(molecule->nAtoms*sizeof(int));
	for(i=0;i<molecule->nAtoms;i++) index[i] = -1;

	ir = rand()%molecule->nAtoms;
	index[n] = ir;
	for(j = 0;j<3;j++) molecule->atoms[index[n]].coordinates[j] = 0;
	for (n = 1; n<molecule->nAtoms;n++)
	{
		boolean accepted;
		Atom* a1 = NULL;
		Atom* a2 = NULL;
		it = 0;
		ir = getRanInt(index, n, molecule->nAtoms);
		index[n] = ir;
		do{
			int j;
			int i;
			double phi;
			double theta;
			ic = rand()%n;
			a2 = &molecule->atoms[index[ic]];
			a1 = &molecule->atoms[index[n]];
			bd = a1->prop.covalentRadii + a2->prop.covalentRadii;
			bd *= BOHRTOANG;
			bd *=0.9;
			//for(j = 0;j<3;j++) a1->coordinates[j] =  a2->coordinates[j]+bd*(0.5+rand()/(double)(RAND_MAX))*(rand()%2==0?1:-1);
			phi = rand()/(double)RAND_MAX*2*M_PI;
			theta = rand()/(double)RAND_MAX*M_PI;
			j = 0;
			a1->coordinates[j] =  a2->coordinates[j]+bd*sin(theta)*cos(phi);
			j = 1;
			a1->coordinates[j] =  a2->coordinates[j]+bd*sin(theta)*sin(phi);
			j = 2;
			a1->coordinates[j] =  a2->coordinates[j]+bd*cos(theta);
			accepted = TRUE;
			for(i=0;i<n;i++)
			{
				int k;
				a2 = &molecule->atoms[index[i]];
				double distance = 0;
				for (k=0;k<3;k++) { double dij = a1->coordinates[k]-a2->coordinates[k]; distance +=dij*dij;}
				distance = sqrt(distance)/BOHRTOANG;
				bd = a1->prop.covalentRadii + a2->prop.covalentRadii;
				//if(distance<0.6*bd)
				if(distance<0.9*bd)
				{
					accepted = FALSE;
					break;
				}
			}
			it++;
		}while(!accepted && it<itmax);
		if(it>=itmax)
		{
			fprintf(stderr,"Error : number of call for rand() > %d\n",itmax);
			fprintf(stderr,"Program stopped in setRandomMolecule of Molecule file\n");
			exit(1);
		}
	}
	//setConnections(molecule);
	return TRUE;

}
/*****************************************************************************/
static boolean setRandomPositionsChain(Molecule* molecule)
{
	int n = 0;
	int ir;
	double bd = 0;
	int ic = 0;
	int it = 0;
	int itmax = 1000*molecule->nAtoms; 
	int j;

	if(molecule->nAtoms<=1) return FALSE;

	for(j = 0;j<3;j++) molecule->atoms[0].coordinates[j] = 0;
	for (n = 1; n<molecule->nAtoms;n++)
	{
		boolean accepted;
		Atom* a1 = NULL;
		Atom* a2 = NULL;
		it = 0;
		ir = n;
		do{
			int j;
			int i;
			double phi;
			double theta;
			ic = n-1;
			a2 = &molecule->atoms[ic];
			a1 = &molecule->atoms[ir];
			bd = a1->prop.covalentRadii + a2->prop.covalentRadii;
			bd *= BOHRTOANG;
			//bd *=0.75;
			bd *=0.9;
			phi = rand()/(double)RAND_MAX*2*M_PI;
			theta = rand()/(double)RAND_MAX*M_PI;
			j = 0;
			a1->coordinates[j] =  a2->coordinates[j]+bd*sin(theta)*cos(phi);
			j = 1;
			a1->coordinates[j] =  a2->coordinates[j]+bd*sin(theta)*sin(phi);
			j = 2;
			a1->coordinates[j] =  a2->coordinates[j]+bd*cos(theta);
			accepted = TRUE;
			for(i=0;i<n;i++)
			{
				int k;
				a2 = &molecule->atoms[i];
				double distance = 0;
				for (k=0;k<3;k++) { double dij = a1->coordinates[k]-a2->coordinates[k]; distance +=dij*dij;}
				distance = sqrt(distance)/BOHRTOANG;
				bd = a1->prop.covalentRadii + a2->prop.covalentRadii;
				//if(distance<0.65*bd)
				if(distance<0.9*bd)
				{
					accepted = FALSE;
					break;
				}
			}
			it++;
		}while(!accepted && it<itmax);
		if(it>=itmax)
		{
			fprintf(stderr,"Error : number of call for rand() > %d\n",itmax);
			fprintf(stderr,"Program stopped in setRandomMolecule of Molecule file\n");
			exit(1);
		}
	}
	return TRUE;

}
/****************************************************************************************************************************/
static void swap2Double(double* a, double *b)
{
	double c = *a;
	*a = *b;
	*b = c;
}
/****************************************************************************************************************************/
static void swap2Int(int* a, int *b)
{
	int c = *a;
	*a = *b;
	*b = c;
}
/****************************************************************************************************************************/
static  int getDistancesFromCM(Molecule* mol, double** pdists, int** pindex, double Cm[], char* symbol)
{
	double* dists = NULL;
	int* index = NULL;
	double mt =0;
	int nAtoms = mol->nAtoms;
	Atom* atoms = mol->atoms;
	*pdists = dists;
	*pindex = index;
	if(mol->nAtoms<1) return 0;
	dists  = malloc(nAtoms*sizeof(double));
	index  = malloc(nAtoms*sizeof(int));
	int k=0;
	int i;
	for(i=0;i<nAtoms;i++) if(!strcmp(symbol,atoms[i].prop.symbol)) { dists[k] = 0;k++;}
	k=0;
	for(i=0;i<nAtoms;i++) if(!strcmp(symbol,atoms[i].prop.symbol)) { index[k] = i;k++;}
	int nA=k;

	for(k=0;k<3;k++) Cm[k]=0;
	for(i=0;i<nAtoms;i++) for(k=0;k<3;k++) Cm[k] += atoms[i].mass*atoms[i].coordinates[k];
	for(i=0;i<nAtoms;i++) mt += atoms[i].mass;
	if(mt>0) for(k=0;k<3;k++) Cm[k] /= mt;

	k=0;
	for(i=0;i<nAtoms;i++) 
	{
		if(strcmp(symbol,atoms[i].prop.symbol)) continue;
		dists[k] = 0;
		int c;
		for(c=0;c<3;c++) dists[k] += (atoms[i].coordinates[c]-Cm[c])*(atoms[i].coordinates[c]-Cm[c]);
		k++;
	}

	for(i=0;i<nA-1;i++) 
	{
		int k=i;
		int j;
		for(j=i+1;j<nA;j++) 
			if(dists[j]<dists[k]) k =j;
		if(k!=i) swap2Double(&dists[i],&dists[k]);
		if(k!=i) swap2Int(&index[i],&index[k]);
	}
	*pdists = dists;
	*pindex = index;
	return nA;
}
/****************************************************************************************************************************/
/*
 * Citation: The Journal of Chemical Physics 138, 214303 (2013);
 * https://doi.org/10.1063/1.4807091
 */
static  Molecule* makeSphereCutSpliceCrossover(Molecule* mol1, Molecule* mol2, int *err)
{
	double* dists1;
	double* dists2;
	int* index1;
	int* index2;
	double Cm1[3];
	double Cm2[3];
	int nTypes=0;

	if(mol1->nAtoms != mol2->nAtoms) return NULL;

	char** types = getTypesOfAtoms(mol1, &nTypes);
	int nTypes2=0;
	char** types2 = getTypesOfAtoms(mol2, &nTypes2);
	if(nTypes2 != nTypes)
	{
		freeTypesOfAtoms(types,nTypes);
		freeTypesOfAtoms(types2,nTypes2);
		return NULL;
	}
	freeTypesOfAtoms(types2,nTypes2);

	*err = -1;
	Molecule* child = mol1->klass->copy(mol1);
	child->potentialEnergy =1e14;
	int it;
	for(it=0;it<nTypes;it++)
	{
		int n1 = getDistancesFromCM(mol1, &dists1, &index1,Cm1,types[it]);
		if(!dists1) {child->klass->free(child); return NULL; }
		int n2 = getDistancesFromCM(mol2, &dists2, &index2,Cm2,types[it]);
		if(!dists2) {child->klass->free(child); return NULL;}
		if(n1!=n2) { child->klass->free(child); return NULL; }
		*err = 1;
		int i;
		for(i=0;i<n1-1;i++)
		{
			int k1=index1[i];
			int k2=index2[i];
			fprintf(stderr,"check %d  %d \n",k1,k2);
			fprintf(stderr,"symbol %s\n",mol1->atoms[k1].prop.symbol);
			fprintf(stderr,"dis1 %f %f dist2 %f %f\n", dists1[i],dists1[i+1],dists2[i],dists2[i+1]);
			if(dists1[i]<dists2[i+1] && dists2[i]<dists1[i+1])
			{
				//if(rand()/(double)(RAND_MAX)>0.5)
				{
					child->atoms[k1] = mol2->atoms[k2];
					int c;
					for(c=0;c<3;c++) child->atoms[k1].coordinates[c] += -Cm2[c]+Cm1[c];
					fprintf(stderr,"take atom %d from second parent place in %d \n",k2,k1);
					*err = 0;
				}
			}
		}
		if(dists1) free(dists1); dists1=NULL;
		if(dists2) free(dists2); dists2=NULL;
		if(index1) free(index1); index1=NULL;
		if(index2) free(index2); index2=NULL;
	}
	freeTypesOfAtoms(types,nTypes);

	/*
	fprintf(stderr,"-----------------------------------------------------------\n");
	fprintf(stderr,"Parent 1\n"); printMolecule(mol1, stderr);
	fprintf(stderr,"Parent 2\n"); printMolecule(mol2, stderr);
	fprintf(stderr,"Molecule after crossing\n"); printMolecule(child, stderr);
	fprintf(stderr,"-----------------------------------------------------------\n");
	fflush(stderr);
	*/
	return child;
}
/****************************************************************************************************************************/
static  void freeTypesOfAtoms(char** types, int nTypes)
{

	int i;
        for (i=0;i<nTypes;i++) if(types[i]) free(types[i]);
        if(types) free(types);
}
/****************************************************************************************************************************/
static  char** getTypesOfAtoms(Molecule* mol, int* pnTypes)
{
	char** types = malloc(mol->nAtoms*sizeof(char*));
	int i;
        for (i=0;i<mol->nAtoms;i++) types[i]=NULL;

	int nTypes = 1;
	types[0]=strdup(mol->atoms[0].prop.symbol);
        for (i=1;i<mol->nAtoms;i++)
	{
		boolean newType=TRUE; 
		int it;
        	for (it=0;it<nTypes;it++)
		if(!strcmp(mol->atoms[i].prop.symbol,types[it]))
		{
			newType=FALSE; 
			break;
		}
		if(newType)
		{
			types[nTypes]=strdup(mol->atoms[i].prop.symbol);
			nTypes++;
		}
	}
	types = realloc(types,nTypes*sizeof(char*));
	*pnTypes = nTypes;
	return types;
}
/****************************************************************************************************************************/
static  void getListProjection(Molecule* mol, double axis[], double** pprojections, int** pindex, double Cm[])
{
	double* projections = NULL;
	int* index = NULL;
	double mt =0;
	int nAtoms = mol->nAtoms;
	Atom* atoms = mol->atoms;
	*pprojections = projections;
	*pindex = index;
	if(mol->nAtoms<1) return;
	projections  = malloc(nAtoms*sizeof(double));
	index  = malloc(nAtoms*sizeof(int));
	int i,k;
	for(i=0;i<nAtoms;i++) projections[i] = 0;
	for(i=0;i<nAtoms;i++) index[i] = i;
	for(k=0;k<3;k++) Cm[k]=0;
	for(i=0;i<nAtoms;i++) for(k=0;k<3;k++) Cm[k] += atoms[i].coordinates[k];
	for(i=0;i<nAtoms;i++) mt += atoms[i].mass;
	if(mt>0) for(k=0;k<3;k++) Cm[k] /= mt;

	for(i=0;i<nAtoms;i++) 
		for(k=0;k<3;k++) 
			projections[i] += (atoms[i].coordinates[k]-Cm[k])*axis[k];

	for(i=0;i<nAtoms;i++) 
	{
		int k=i;
		int j;
		for(j=i+1;j<nAtoms;j++) 
			if(projections[j]>projections[k]) k =j;
		if(k!=i) swap2Double(&projections[i],&projections[k]);
		if(k!=i) swap2Int(&index[i],&index[k]);
	}
	*pprojections = projections;
	*pindex = index;
}
/****************************************************************************************************************************/
/*
 * Citation: The Journal of Chemical Physics 138, 214303 (2013);
 * https://doi.org/10.1063/1.4807091
 */
static  Molecule* makePlaneCutSpliceCrossoverWithoutTestBonds(Molecule* mol1, Molecule* mol2, int *err)
{
	double* projections1;
	double* projections2;
	int* index1;
	int* index2;
	double Cm1[3];
	double Cm2[3];
	int nAtoms = mol1->nAtoms;
	double axis[3];
	double rho = 1.0;
	double theta = rand()/(double)(RAND_MAX)*M_PI;
	double phi = 2.0*rand()/(double)(RAND_MAX)*M_PI;

	axis[0]=rho*cos(theta)*cos(phi);
	axis[1]=rho*cos(theta)*sin(phi);
	axis[2]=rho*sin(theta);

	//fprintf(stderr,"axis= %f %f %f\n",axis[0],axis[1],axis[2]); fflush(stderr);

	*err = -1;
	if(mol1->nAtoms != mol2->nAtoms) return NULL;
	getListProjection(mol1, axis, &projections1, &index1,  Cm1);
	if(!projections1) return NULL;
	getListProjection(mol2, axis, &projections2, &index2,  Cm2);
	if(!projections2) return NULL;

	int nTypes=0;
	char** types = getTypesOfAtoms(mol1, &nTypes);
	int nTypes2=0;
	char** types2 = getTypesOfAtoms(mol2, &nTypes2);
	if(nTypes2 != nTypes)
	{
		freeTypesOfAtoms(types,nTypes);
		freeTypesOfAtoms(types2,nTypes2);
		if(projections1) free(projections1);
		if(projections2) free(projections2);
		if(index1) free(index1);
		if(index2) free(index2);
		return NULL;
	}
	freeTypesOfAtoms(types2,nTypes2);

	Molecule* child = mol1->klass->copy(mol1);
	child->potentialEnergy =1e14;
	// TO REMOVE return child;


	int it;
	for(it=0;it<nTypes;it++)
	{
		int na1 = 0; 
		int i;
		for(i=0;i<nAtoms;i++) if(!strcmp(mol1->atoms[i].prop.symbol,types[it])) na1++;
		int na2 = 0; 
		for(i=0;i<nAtoms;i++) if(!strcmp(mol1->atoms[i].prop.symbol,types[it])) na2++;
		int n1 = na1/2+na1%2;
		int n2 = na2-n1;
		if(n2<0) n2=0;
		int it2=0;
		for(i=nAtoms-1;i>=0;i--) 
		{
			int k2=index2[i];
			if(!strcmp(mol2->atoms[k2].prop.symbol,types[it])) 
			{
				int it1=0;
				it2++;
				int j;
				for(j=nAtoms-1;j>=0;j--) 
				{
					int k1=index1[j];
					if(!strcmp(mol1->atoms[k1].prop.symbol,types[it])) 
					{
						it1++;
						if(it1==it2)
						{
							int *con = child->atoms[k1].typeConnections;
							child->atoms[k1]=getCopyAtom(&mol2->atoms[k2]);
							child->atoms[k1].typeConnections=con;
							int c;
							for(c=0;c<3;c++) child->atoms[k1].coordinates[c] += Cm1[c]-Cm2[c];
							break;
						}
					}
				}
				if(it2==n2) break;
			}
		}
	}
	freeTypesOfAtoms(types,nTypes);
	if(projections1) free(projections1);
	if(projections2) free(projections2);
	if(index1) free(index1);
	if(index2) free(index2);

	/*
	fprintf(stderr,"-----------------------------------------------------------\n");
	fprintf(stderr,"Parent 1\n"); printMolecule(mol1, stderr);
	fprintf(stderr,"Parent 2\n"); printMolecule(mol2, stderr);
	fprintf(stderr,"Molecule after crossing\n"); printMolecule(child, stderr);
	fprintf(stderr,"-----------------------------------------------------------\n");
	fflush(stderr);
	*/
	return child;
}
static  Molecule* makePlaneCutSpliceCrossover(Molecule* mol1, Molecule* mol2, int *err)
{
	int itmax = 1000; 
	Molecule* child = NULL;
	boolean accepted=FALSE;
	int it;
	for(it=0;it<itmax && !accepted ;it++)
	{
		//fprintf(stderr,"it=%d\n",it); fflush(stderr);
		if(child) child->klass->free(child);
		child = makePlaneCutSpliceCrossoverWithoutTestBonds(mol1,mol2, err);
		if(!child) continue;
		if(!child->klass->smallDistance(child)) { accepted = TRUE; break;}
	}
	return child;
}
/********************************************************************************/
static boolean makeCenterOfMassSphericalMutation(Molecule* molecule)
{
	int itmax = 1000; 
	int it = 0;
	Atom* atoms = molecule->atoms;
	double C[3];
	double Cm[3];
	int k;

	if(molecule->nAtoms<=1) return FALSE;
	for (k=0;k<3;k++) Cm[k] = 0;
	double mt = 0;
	int i;
	for (i = 0; i<molecule->nAtoms;i++) 
	{
		int k;
		for (k=0;k<3;k++) 
			Cm[k] += atoms[i].mass*atoms[i].coordinates[k];
		mt += atoms[i].mass; 
	}
	if(mt>0) for (k=0;k<3;k++) Cm[k] /= mt;

	boolean accepted=FALSE;
	do{
		int ir = rand()%molecule->nAtoms;
		double rho = 0;
		int k;
		for (k=0;k<3;k++) 
		{ 
			double dij = atoms[ir].coordinates[k]-Cm[k];
			rho +=dij*dij;
		}
		rho = sqrt(rho);
		double theta = rand()/(double)(RAND_MAX)*M_PI;
		double phi = rand()/(double)(RAND_MAX)*2.0*M_PI;
		C[0] = rho*cos(phi)*cos(theta);
		C[1] = rho*sin(phi)*cos(theta);
		C[2] = rho*sin(theta);

		int imin=-1;
		double dmin=-1;
		int i;
		for (i = 0; i<molecule->nAtoms;i++)
		{
			if(i==ir) continue;
			double distance = 0;
			int k;
			for (k=0;k<3;k++) 
			{ 
				double dij = atoms[i].coordinates[k]-C[k]; 
				distance +=dij*dij;
			}
			if(dmin<0) { dmin = distance; imin=i;}
			if(distance<dmin) { dmin = distance; imin=i;}
		}
		if(imin>-1) 
		{
			double bd = atoms[ir].prop.covalentRadii + atoms[imin].prop.covalentRadii;
			bd *= BOHRTOANG;
			bd *=0.6;
			if(sqrt(dmin)>bd) accepted=TRUE;
		}
		if(accepted) for (k=0;k<3;k++) atoms[ir].coordinates[k] = C[k];
		it++;
	}while(!accepted && it<itmax);
	//fprintf(stderr,"accepte = %d, it/itmax =%d/%d\n",accepted,it,itmax); fflush(stderr);
	//setConnections(molecule);
	return accepted;

}
/********************************************************************************/
static boolean setRandomOrientation(Molecule* molecule)
{
	double axis1[3];
	double axis2[3]; 
	double axis3[3];

	if(molecule->nAtoms<1) return FALSE;

	double rho = 1.0;
	double theta = rand()/(double)(RAND_MAX)*M_PI;
	double phi = 2.0*rand()/(double)(RAND_MAX)*M_PI;

	axis1[0]=rho*cos(theta)*cos(phi);
	axis1[1]=rho*cos(theta)*sin(phi);
	axis1[2]=rho*sin(theta);

	theta = rand()/(double)(RAND_MAX)*M_PI;
	phi = 2.0*rand()/(double)(RAND_MAX)*M_PI;

	axis2[0]=rho*cos(theta)*cos(phi);
	axis2[1]=rho*cos(theta)*sin(phi);
	axis2[2]=rho*sin(theta);

	int is=0;
	int k;
	for(k=0;k<3;k++) if(axis1[k]!=0) is=k;

	double ps=0;
	for(k=0;k<3;k++) if(k!=is) ps += axis1[k]*axis2[k];
	axis2[is]=-ps/axis1[is];

	for(k=0;k<3;k++) axis3[k] = axis1[(k+1)%3]*axis2[(k+2)%3]-axis1[(k+2)%3]*axis2[(k+1)%3];
	return setGeometryToAxes(molecule,axis1, axis2, axis3);
}
/********************************************************************************/
static boolean setGeometryToAxes(Molecule* molecule, double axis1[], double axis2[], double axis3[])
{
	Atom* atoms = molecule->atoms;
	int nAtoms = molecule->nAtoms;
	double tensor[3][3];
	double invTensor[3][3];
	double** C;
	int i,j,k;
	boolean ok = FALSE;

	if(nAtoms<1) return ok;

	C = malloc(nAtoms*sizeof(double*));
	for(i=0;i<nAtoms;i++) C[i] = malloc(3*sizeof(double));

	tensor[0][0] = axis1[0];
	tensor[0][1] = axis1[1];
	tensor[0][2] = axis1[2];

	tensor[1][0] = axis2[0];
	tensor[1][1] = axis2[1];
	tensor[1][2] = axis2[2];

	tensor[2][0] = axis3[0];
	tensor[2][1] = axis3[1];
	tensor[2][2] = axis3[2];

	if(InverseTensor(tensor,invTensor))
	{
		double Cm[3];
		double mt = 0;
		for(k=0;k<3;k++) Cm[k]=0;
		for(i = 0;i<nAtoms;i++) for(k=0;k<3;k++) Cm[k] += atoms[i].mass*atoms[i].coordinates[k];
		for(i = 0;i<nAtoms;i++) mt += atoms[i].mass;
		if(mt>0) for(k=0;k<3;k++) Cm[k] /= mt;

		for(i = 0;i<nAtoms;i++)
		{
			for(j=0;j<3;j++)
			{
				C[i][j] = 0.0;
				for(k=0;k<3;k++) C[i][j] += invTensor[k][j]*(atoms[i].coordinates[k]-Cm[k]);
			}
		}
		for(i = 0;i<nAtoms;i++) for(k=0;k<3;k++) atoms[i].coordinates[k]=C[i][k]+Cm[k];
		ok = TRUE;
	}
	for(i=0;i<nAtoms;i++) if(C[i]) free(C[i]);
	if(C) free(C);
	return ok;
}
/********************************************************************************/
static boolean makeLocalSphericalMutation(Molecule* molecule, double rate)
{
	int itmax = 1000; 
	int it = 0;
	Atom* atoms = molecule->atoms;
	double C[3];

	if(molecule->nAtoms<=1) return FALSE;
	boolean accepted=FALSE;
	do{
		int ir = rand()%molecule->nAtoms;
		double dmin=-1;
		int i;
		for (i = 0; i<molecule->nAtoms;i++)
		{
			if(i==ir) continue;
			double distance = 0;
			int k;
			for (k=0;k<3;k++) 
			{ 
				double dij = atoms[i].coordinates[k]-atoms[ir].coordinates[k]; 
				distance +=dij*dij;
			}
			if(dmin<0) dmin = distance;
			if(distance<dmin) dmin = distance;
		}
		double rho = (-1+2*rand()/(double)(RAND_MAX))*rate*sqrt(dmin);
		double theta = rand()/(double)(RAND_MAX)*M_PI;
		double phi = rand()/(double)(RAND_MAX)*2.0*M_PI;
		C[0] = rho*cos(phi)*cos(theta);
		C[1] = rho*sin(phi)*cos(theta);
		C[2] = rho*sin(theta);

		int imin=-1;
		dmin=-1;
		for (i = 0; i<molecule->nAtoms;i++)
		{
			if(i==ir) continue;
			double distance = 0;
			int k;
			for (k=0;k<3;k++) 
			{ 
				double dij = atoms[i].coordinates[k]-C[k]; 
				distance +=dij*dij;
			}
			if(dmin<0) { dmin = distance; imin=i;}
			if(distance<dmin) { dmin = distance; imin=i;}
		}
		if(imin>-1) 
		{
			double bd = atoms[ir].prop.covalentRadii + atoms[imin].prop.covalentRadii;
			bd *= BOHRTOANG;
			bd *=0.6;
			if(sqrt(dmin)>bd) accepted=TRUE;
		}
		int k;
		if(accepted) for (k=0;k<3;k++) atoms[ir].coordinates[k] = C[k];
		it++;
	}while(!accepted && it<itmax);
	//setConnections(molecule);
	return accepted;

}
/********************************************************************************/
static boolean readGeometry(Molecule* molecule, char* namefile)
{
	FILE* file = NULL;
	static char *t = NULL; 
	boolean Ok = FALSE;
	int n,ic,is;
	Molecule* mol = newMolecule();
#define SZ 50
	char symbol[SZ];
	char mmType[SZ];
	char pdbType[SZ];
	char residueName[SZ];
	double X,Y,Z;
	double charge;
	int layer;
	int i,l;
	char* pos;

	if(!namefile) 
	{
		printf("Sorry I cannot read geometry namefile = NULL\n");
		exit(1);
	}
	if(t==NULL) t = malloc(BSIZE*sizeof(char));
	file = fopen(namefile,"rb");
	if(!file)
	{
		printf("Sorry I cannot open %s file\n",namefile);
		exit(1);
	}
	rewind(file);
	while(!feof(file))
  	{
    		if(!fgets(t,BSIZE, file)) break;
		deleteFirstSpaces(t);
		if(t[0]=='#') continue;
		uppercase(t);
		pos = strstr(t,"GEOMETRY");
		if(pos)
		{ 
			while(!feof(file))
			{
    				if(!fgets(t,BSIZE, file)) break;
				deleteFirstSpaces(t);
				if(t[0]=='#') continue;
				break;
			}
			if(3==sscanf(t,"%d%d%d",&n,&ic,&is) && n>0 && is>0)
			{
				/*
				printf("Mult = %d\n",is);
				printf("nAtoms = %d\n",n);
				*/
				mol->nAtoms = n;
				mol->spinMultiplicity = is;
				mol->totalCharge = ic;
				mol->atoms = malloc(mol->nAtoms*sizeof(Atom));
				for(i=0; i<mol->nAtoms; i++) mol->atoms[i].typeConnections = malloc(mol->nAtoms*sizeof(int));
				Ok = TRUE;
			}
			break;
		}
	}
	if(!Ok)
	{
		printf("Sorry I cannot read geometry from %s file\n",namefile);
	 	return FALSE;
	//	exit(1);
	}
	for(i=0; i<mol->nAtoms; i++)
	{
			int variable = 0;
			int ibeg = 12;
			if(!fgets(t,BSIZE,file))
			{
				printf("Sorry I cannot read geometry from %s file.\n",namefile);
				return FALSE;
				//exit(1);
			}
			deleteFirstSpaces(t);
			if(t[0]=='#') { i--;continue;}
    			sscanf(t,"%s %s %s %s %d %lf %d %d %lf %lf %lf",
					symbol,mmType,pdbType,residueName, 
					&mol->atoms[i].residueNumber,
					&charge,&layer,&variable,&X,&Y,&Z);
			symbol[0]=toupper(symbol[0]);
			l=strlen(symbol);
			if (l==2) symbol[1]=tolower(symbol[1]);

			mol->atoms[i].prop = propAtomGet(symbol);
			mol->atoms[i].mmType=strdup(mmType);
			mol->atoms[i].pdbType=strdup(pdbType);
			mol->atoms[i].residueName=strdup(residueName);
			mol->atoms[i].N=i+1;
			mol->atoms[i].layer=layer;
			mol->atoms[i].variable=variable;
			mol->atoms[i].show=TRUE;
			mol->atoms[i].coordinates[0] = X;
			mol->atoms[i].coordinates[1] = Y;
			mol->atoms[i].coordinates[2] = Z;
			mol->atoms[i].velocity[0] = 0;
			mol->atoms[i].velocity[1] = 0;
			mol->atoms[i].velocity[2] = 0;
			mol->atoms[i].charge = charge;
			mol->atoms[i].charge0 = charge;
			mol->atoms[i].electronegativity = 0;
			mol->atoms[i].hardness = 0;
			mol->atoms[i].width = mol->atoms[i].prop.covalentRadii;
			mol->atoms[i].mass = mol->atoms[i].prop.mass;
			mol->atoms[i].rho = 0.0;
	   		mol->atoms[i].U = 0.0;
			if(!get_connections_one_atom(t, mol->nAtoms, ibeg, mol->atoms[i].typeConnections))
			{
				/*
				printf("Sorry I cannot read the connection for atom # %d from the %s file.\n",i+1,namefile);
				exit(1);
				*/
				fprintf(stderr,"Warning : I cannot read the connection for atom # %d from the %s file.\n",i+1,namefile);
			}
	}
	readBoxes(mol, file, t);
	readWall(mol, file, t);
	
	fclose(file);
	mol->potentialEnergy = molecule->potentialEnergy;
	*molecule = *mol;
	return TRUE;
}
/*******************************************************************************************************************/
static void setBondHardness(Molecule* molecule, int nBonds, char** atomTypes1, char** atomTypes2, double* hardness)
{
	int i;
	int j;
	int k;
	int n;
	double* bondHardness = NULL;
	int numberOfAtoms = molecule->nAtoms;
	Atom* atoms = molecule->atoms;
	
	if(molecule->bondHardness) free(molecule->bondHardness);

	if(numberOfAtoms<2) return;
	bondHardness = malloc(numberOfAtoms*(numberOfAtoms+1)/2*sizeof(double));

	k=0;
	//printf("begin setBondHardness\n");
	for (  i = 0; i < numberOfAtoms; i++ )
		for (  j = i + 1; j < numberOfAtoms; j++ )
		{
			bondHardness[k] = 0;
			for(n=0;n<nBonds;n++)
			{
				if(
				(!strcmp(atomTypes1[n],atoms[i].mmType) && !strcmp(atomTypes2[n],atoms[j].mmType))
				||
				(!strcmp(atomTypes1[n],atoms[j].mmType) && !strcmp(atomTypes2[n],atoms[i].mmType))
				)
				{
					bondHardness[k] = hardness[n]/AUTOEV;
				}
			}
			k++;
		}
	molecule->bondHardness = bondHardness;
	//printf("end setBondHardness\n");
	/* printing for test*/
}
/*******************************************************************************************************************/
static void setHardness(Molecule* molecule, int nTypes, char** atomTypes, double* hardness)
{
	int i;
	int n;
	int numberOfAtoms = molecule->nAtoms;
	Atom* atoms = molecule->atoms;
	
	for (  i = 0; i < numberOfAtoms; i++ )
	{
		atoms[i].hardness = 0;
		for(n=0;n<nTypes;n++)
		{
			if(!strcmp(atomTypes[n],atoms[i].mmType)) atoms[i].hardness = hardness[n]/AUTOEV;//convert in Hartree
			//printf("atomTypes = %s mmType = %s\n",atomTypes[n],atoms[i].mmType);
		}
	}
	/* printing for test*/
}
/*******************************************************************************************************************/
static void setElectronegativity(Molecule* molecule, int nTypes, char** atomTypes, double* electronegativity)
{
	int i;
	int n;
	int numberOfAtoms = molecule->nAtoms;
	Atom* atoms = molecule->atoms;
	
	for(i = 0; i < numberOfAtoms; i++)
	{
		atoms[i].electronegativity = 0;
		for(n=0;n<nTypes;n++)
			if(!strcmp(atomTypes[n],atoms[i].mmType))
				atoms[i].electronegativity = electronegativity[n]/AUTOEV; // convert in Hartree
	}
	/* printing for test*/

}
/*******************************************************************************************************************/
static void setWidth(Molecule* molecule, int nTypes, char** atomTypes, double* width)
{
	int i;
	int n;
	int numberOfAtoms = molecule->nAtoms;
	Atom* atoms = molecule->atoms;
	
	for(i = 0; i < numberOfAtoms; i++)
	{
		atoms[i].width = 0;
		for(n=0;n<nTypes;n++)
			if(!strcmp(atomTypes[n],atoms[i].mmType))
				atoms[i].width = width[n];
	}
	/* printing for test*/
}
/*******************************************************************************************************************/
static void setCharge0(Molecule* molecule, int nTypes, char** atomTypes, double* charge0)
{
	int i;
	int n;
	int numberOfAtoms = molecule->nAtoms;
	Atom* atoms = molecule->atoms;
	
	for (  i = 0; i < numberOfAtoms; i++ )
	{
		atoms[i].charge0 = 0;
		for(n=0;n<nTypes;n++)
		{
			if(!strcmp(atomTypes[n],atoms[i].mmType)) atoms[i].charge0 = charge0[n];
			//printf("atomTypes = %s mmType = %s\n",atomTypes[n],atoms[i].mmType);
		}
	}
	/* printing for test*/
}
/*******************************************************************************************************************/
static void setup_EEM(Molecule* molecule, double** A, double* B)
{
	int i;
	int j;
	int numberOfAtoms = molecule->nAtoms;
	double rijx, rijy, rijz;
	double rij;
	double r0;
	
	for (  i = 0; i < numberOfAtoms; i++ )
	{
		A[i][i] = molecule->atoms[i].hardness;
		B[i] =-molecule->atoms[i].electronegativity;
		for (  j = i + 1; j < numberOfAtoms; j++ )
		{
			rijx = molecule->atoms[i].coordinates[0] - molecule->atoms[j].coordinates[0];
			rijy = molecule->atoms[i].coordinates[1] - molecule->atoms[j].coordinates[1];
			rijz = molecule->atoms[i].coordinates[2] - molecule->atoms[j].coordinates[2];
			rij = sqrt( rijx * rijx + rijy * rijy + rijz * rijz );
		
			A[i][j] = A[j][i] = 0;
			r0 = sqrt(
				 2*molecule->atoms[i].width*molecule->atoms[i].width
				+2*molecule->atoms[j].width*molecule->atoms[j].width
				);
			//printf("r0=%f\n",r0);
			if(fabs(rij)>1e-10 && r0>1e-10) A[i][j] = A[j][i] = erf(rij/r0)/(rij*ANGTOBOHR);// in Hartree
			else if(fabs(rij)>1e-10) A[i][j] = A[j][i] = 1/(rij*ANGTOBOHR);// in Hartree
		}
	}
	//  Fill in the constraints
	for (  i = 0; i < numberOfAtoms; i++ )
	{
		j = numberOfAtoms;
		A[i][j] = A[j][i] = 1;
	}
	i = j = numberOfAtoms;
	A[i][j] = 0;

	for (  i = 0; i < numberOfAtoms; i++ )
		B[i] += molecule->atoms[i].hardness*molecule->atoms[i].charge0;
	i = numberOfAtoms;
	B[i] = molecule->totalCharge;
}
/***********************************************************************************************/
static double getEnergyEEM(Molecule* molecule)
{
	int ia,ja;
	int numberOfAtoms = molecule->nAtoms;
	double energy = 0;

	for (  ia = 0; ia < numberOfAtoms; ia++ )
	{
		double Delta = molecule->atoms[ia].charge0 - molecule->atoms[ia].charge;
		double mu    = -molecule->atoms[ia].electronegativity;
		energy += mu*Delta;
	}
	for (  ia = 0; ia < numberOfAtoms; ia++ )
	{
		double Deltai =  molecule->atoms[ia].charge0 - molecule->atoms[ia].charge;
		double eta   = molecule->atoms[ia].hardness;
		energy += 0.5*Deltai*Deltai*eta;
		for (  ja = ia + 1; ja < numberOfAtoms; ja++ ) // ja>ia, for 0.5
		{
			double rijx, rijy, rijz, rij;
			double r0;
			double etaij;
			double Deltaj =  molecule->atoms[ja].charge0 - molecule->atoms[ja].charge;
			
			rijx = molecule->atoms[ia].coordinates[0] - molecule->atoms[ja].coordinates[0];
			rijy = molecule->atoms[ia].coordinates[1] - molecule->atoms[ja].coordinates[1];
			rijz = molecule->atoms[ia].coordinates[2] - molecule->atoms[ja].coordinates[2];
			rij = sqrt( rijx * rijx + rijy * rijy + rijz * rijz );
		
			r0 = sqrt(2*molecule->atoms[ia].width*molecule->atoms[ia].width+2*molecule->atoms[ja].width*molecule->atoms[ja].width);

			//printf("r0=%f\n",r0);
			etaij = 0;
			if(fabs(rij)>1e-10 && r0>1e-10) etaij = erf(rij/r0)/(rij*ANGTOBOHR);// in Hartree
			else if(fabs(rij)>1e-10) etaij = 1/(rij*ANGTOBOHR);// in Hartree

			etaij *=  Deltai*Deltaj;
			energy += etaij;
		}
	}
	return energy*AUTOKCAL;
}
/***********************************************************************************************/
static double getEnergyACKS2(Molecule* molecule)
{
	int ia,ja;
	int numberOfAtoms = molecule->nAtoms;
	double energy = 0;

	energy = getEnergyEEM(molecule)/AUTOKCAL;

	for (  ia = 0; ia < numberOfAtoms; ia++ ) 
	{
		double Delta = molecule->atoms[ia].charge0 - molecule->atoms[ia].charge;
		double U = molecule->atoms[ia].U;
		energy += -U*Delta;
	}
	if(molecule->bondHardness)
	{
		int k =0;
		double** A = malloc(numberOfAtoms*sizeof(double*));
		for(ia=0;ia<numberOfAtoms;ia++) A[ia] = malloc(numberOfAtoms*sizeof(double));
		for(ia=0;ia<numberOfAtoms;ia++) for(ja=0;ja<numberOfAtoms;ja++)  A[ia][ja] = 0; 

		k = 0;
		for(ia=0;ia<numberOfAtoms;ia++) 
		for(ja=ia+1;ja<numberOfAtoms;ja++) 
		{
			//printf("ia = %d ja = %d k = %d kappa = %f\n",ia, ja, k, molecule->bondHardness[k]);
			if(fabs(molecule->bondHardness[k]) >1e-10)
			{
				double bsoft = 1/molecule->bondHardness[k];
				A[ia][ja] += bsoft;
				A[ja][ia] += bsoft;
				A[ia][ia] -= bsoft;
				A[ja][ja] -= bsoft;
			}
			k++;
		}
		for (  ia = 0; ia < numberOfAtoms; ia++ )
		{
			double Ui =  molecule->atoms[ia].U;
			for (  ja = 0; ja < numberOfAtoms; ja++ ) // ja>ia, for 0.5
			{
				double Xij = A[ia][ja];
				double Uj =  molecule->atoms[ja].U;
				energy += 0.5*Ui*Uj*Xij;
			}
		}
		for(ia=0;ia<numberOfAtoms;ia++) free(A[ia]);
		free(A);
	}
	return energy*AUTOKCAL;
}
/***********************************************************************************************/
static void setChargesEEM(Molecule* molecule)
{
	double* B;
	double** A;
	double* values;
	int ia;
	int i,j;
	int numberOfAtoms = molecule->nAtoms;
	int n = molecule->nAtoms+1;
	B = malloc(n*sizeof(double));
	A = malloc(n*sizeof(double*));
	for(i=0;i<n;i++) A[i] = malloc(n*sizeof(double));

	for(i=0;i<n;i++) for(j=0;j<n;j++)  A[i][j] = 0; 
	for(i=0;i<n;i++) B[i] = 0;
	setup_EEM(molecule, A, B);

	/* test solveSymEqQL */
	/*
	for(i=0;i<n;i++) B[i] = 1;
	for(i=0;i<n;i++) for(j=0;j<n;j++) A[i][j] = 1.0;
	A[1][1] = -1.0;
	B[1] = 2;
	*/
	
	// to check
	/*
	printf("EEM Matrix = \n");
	printMatrixDouble(A, n, n);
	*/

	values =  malloc(n*sizeof(double));

	// TO CHANGE
	// test solveSymEqQL n = 2;
        solveSymEqQL(n, A, B, values);
	//for(i=0;i<n;i++) printf("Val %d = %f\n",i,values[i]);
	//double S = 0; for(i=0;i<n-1;i++) S+= values[i]; printf("Sum charges=%f\n",S);

	for (  ia = 0; ia < numberOfAtoms; ia++ ) molecule->atoms[ia].charge = values[ia];
	for (  ia = 0; ia < numberOfAtoms; ia++ ) molecule->atoms[ia].U = 0.0;
	free(values);
	free(B);
	for(i=0;i<n;i++) free(A[i]);
	free(A);
}
/***********************************************************************************************/
static void setChargesACKS2(Molecule* molecule)
{
	double* B;
	double** A;
	double* values;
	int i,j,k;
	int ia,ja;
	int numberOfAtoms = molecule->nAtoms;
	int n = 2*(molecule->nAtoms+1);
	B = malloc(n*sizeof(double));
	A = malloc(n*sizeof(double*));
	for(i=0;i<n;i++) A[i] = malloc(n*sizeof(double));
	for(i=0;i<n;i++) B[i] = 0;
	for(i=0;i<n;i++) for(j=0;j<n;j++)  A[i][j] = 0; 

	setup_EEM(molecule, A, B);
	//  Fill in the constraints
	for(ia=0;ia<numberOfAtoms;ia++) B[ia+numberOfAtoms+1] = molecule->atoms[ia].charge0; 
	for(ia=0;ia<numberOfAtoms;ia++) A[2*numberOfAtoms+1][ia+numberOfAtoms+1] = 1;
	for(ia=0;ia<numberOfAtoms;ia++) A[ia+numberOfAtoms+1][2*numberOfAtoms+1] = 1;

	// Fill in off diagonal identity matrix bloks
	for(ia=0;ia<numberOfAtoms;ia++) A[ia][ia+numberOfAtoms+1] = 1; 
	for(ia=0;ia<numberOfAtoms;ia++) A[ia+numberOfAtoms+1][ia] = 1; 
	// Add the bond hardness terms as off-diagonal bond-softness parameters
	//printf("Begin bond hardness\n");
	k = 0;
	if(molecule->bondHardness)
	for(ia=0;ia<numberOfAtoms;ia++) 
	for(ja=ia+1;ja<numberOfAtoms;ja++) 
	{
		//printf("ia = %d ja = %d k = %d kappa = %f\n",ia, ja, k, molecule->bondHardness[k]);
		if(fabs(molecule->bondHardness[k]) >1e-10)
		{
			i = numberOfAtoms+1+ia;
			j = numberOfAtoms+1+ja;
			double bsoft = 1/molecule->bondHardness[k];
			A[i][j] += bsoft;
			A[j][i] += bsoft;
			A[i][i] -= bsoft;
			A[j][j] -= bsoft;
		}
		k++;
	}
	// to check
	/*
	printf("ACKS2 = \n");
	printMatrixDouble(A, n, n);
	*/

	values =  malloc(n*sizeof(double));
        solveSymEqQL(n, A, B, values);
	for (  ia = 0; ia < numberOfAtoms; ia++ ) molecule->atoms[ia].charge = values[ia];
	for (  ia = 0; ia < numberOfAtoms; ia++ ) molecule->atoms[ia].U = values[numberOfAtoms+1+ia];
	free(values);
	free(B);
	for(i=0;i<n;i++) free(A[i]);
	free(A);
}
/********************************************************************************/
static boolean saveGeometry(Molecule* molecule, char* fileName)
{
 	FILE* file = fopen(fileName, "w");
	boolean ok=TRUE;
	if(!file) return FALSE;
	ok = addGeometry(molecule,file);
	fclose(file);
	return ok;
}
/********************************************************************************/
static boolean saveGeometryAndVelocities(Molecule* molecule, char* fileName)
{
 	FILE* file = fopen(fileName, "w");
	boolean ok=TRUE;
	if(!file) return FALSE;
	ok = addGeometry(molecule,file);
	if(ok)
	{
		ok = addVelocities(molecule,file);
	}
	fclose(file);
	return ok;
}
/*****************************************************************************/
static void copyGradients(Molecule* mol, double* g[])
{
        int i,k;
        if(!mol) return;
        for(i=0;i<mol->nAtoms;i++)
                for(k=0;k<3;k++)
                        g[k][i] = mol->atoms[i].gradient[k];
}
/*****************************************************************************/
static void removeTransRotModes(Molecule* mol)
{
	int i;
	int k;
	int nToRemove;
	if(!mol ||  mol->vibration.nModes<1) return;
	nToRemove = 6;
	if(isLinear(mol)) nToRemove = 5;
	if(mol->vibration.nModes != 3*mol->nAtoms) return;

	sortFrequencies(mol);
	for(i=0;i< mol->vibration.nModes-nToRemove;i++)
	{
		k = i+nToRemove;
		mol->vibration.modes[i] =  mol->vibration.modes[k];
	}
	mol->vibration.nModes -= nToRemove;
	mol->vibration.modes =  realloc(mol->vibration.modes,mol->vibration.nModes*sizeof(VibMode));
}
/*****************************************************************************/
static void sortFrequencies(Molecule* mol)
{
	int i;
	int j;
	int k;
	VibMode dum;
	if(!mol ||  mol->vibration.nModes<1) return;
	for(i=0;i< mol->vibration.nModes;i++)
	{
		k = i;
		for(j=i+1;j< mol->vibration.nModes;j++)
			if( mol->vibration.modes[j].frequency< mol->vibration.modes[k].frequency) k = j;
		if(k==i) continue;
		/* swap i and k modes */
		dum =  mol->vibration.modes[i];
		mol->vibration.modes[i] =  mol->vibration.modes[k];
		mol->vibration.modes[k] =  dum;
	}
}
/*****************************************************************************/
static void removeFrequencies(Molecule* mol, double freqMin, double freqMax)
{
	int i;
	int j;
	int k;
	int iMin = 0;
	int iMax= 0;
	int nM = 0;
	if(!mol ||  mol->vibration.nModes<1) return;
	if(freqMin<0 && freqMax<0) return;
	if(freqMin>0 && freqMax>0 && freqMin>freqMax) 
	{
		double f = freqMin;
		freqMin = freqMax;
		freqMax = f;
	}
	sortFrequencies(mol);
	iMin = 0;
	if(freqMin>0)
	for(i=0;i< mol->vibration.nModes;i++)
		if( mol->vibration.modes[i].frequency < freqMin) iMin = i+1;
	iMax = mol->vibration.nModes-1;
	if(freqMax>0)
	for(i=mol->vibration.nModes-1;i>0;i--) 
		if( mol->vibration.modes[i].frequency > freqMax) iMax = i-1;

	//fprintf(stderr,"iMin = %d\n",iMin);
	//fprintf(stderr,"iMax = %d\n",iMax);
	for(i=0;i< iMin;i++)
        {
                mol->vibration.modes[i].frequency = 0;
                mol->vibration.modes[i].mass = 1;
                for(j=0;j<3;j++) if(mol->vibration.modes[i].vectors[j]) free(mol->vibration.modes[i].vectors[j]);
                if(mol->vibration.modes[i].properties) free(mol->vibration.modes[i].properties);
        }
	for(i=iMax+1;i<mol->vibration.nModes;i++)
        {
                mol->vibration.modes[i].frequency = 0;
                mol->vibration.modes[i].mass = 1;
                for(j=0;j<3;j++) if(mol->vibration.modes[i].vectors[j]) free(mol->vibration.modes[i].vectors[j]);
                if(mol->vibration.modes[i].properties) free(mol->vibration.modes[i].properties);
        }
	nM = iMax-iMin+1;
	k = 0;
	for(i=iMin;i<=iMax;i++)
	{
		mol->vibration.modes[k] =  mol->vibration.modes[i];
		k++;
	}

	mol->vibration.nModes = nM;
	if(nM>0) mol->vibration.modes = realloc(mol->vibration.modes, mol->vibration.nModes*sizeof(VibMode));
	else if(mol->vibration.modes) free(mol->vibration.modes);
}
/*****************************************************************************/
static void removeNotSelectedFrequencies(Molecule* mol, boolean* selected)
{
	int i;
	int j;
	int n = 0;
	if(!mol ||  mol->vibration.nModes<1) return;
	for(i=0;i< mol->vibration.nModes;i++) if(selected[i]) n++;

	if(n==0 || n == mol->vibration.nModes) return;
	for(i=0;i<mol->vibration.nModes;i++)
        {
		if(selected[i]) continue;

                for(j=0;j<3;j++) if(mol->vibration.modes[i].vectors[j]) free(mol->vibration.modes[i].vectors[j]);
                if(mol->vibration.modes[i].properties) free(mol->vibration.modes[i].properties);

		for(j=i+1;j<mol->vibration.nModes;j++)
			if(selected[j]){ 
				int k;
				selected[j] = FALSE;
				selected[i] = TRUE;
				mol->vibration.modes[i] =  mol->vibration.modes[j];
                		for(k=0;k<3;k++) mol->vibration.modes[j].vectors[k] = NULL;
				mol->vibration.modes[j].properties = NULL;
				break;
			}
        }
	mol->vibration.nModes = n;
	if(n>0) mol->vibration.modes = realloc(mol->vibration.modes, mol->vibration.nModes*sizeof(VibMode));
	else if(mol->vibration.modes) free(mol->vibration.modes);
}
/****************************************************************************************************************************************************/
static int computeIR(Molecule* mol, double *F, double* dmuX[3])
{
	int nAtoms = mol->nAtoms;
	int i;
	int j;
	int k;
	int c;
	double* frequencies = NULL;
	double** modes = NULL;


	freeVibrations(mol);
	initVibrations(mol, 3*nAtoms, 1);
	frequencies = malloc(3*nAtoms*sizeof(double));
	modes = malloc(3*nAtoms*sizeof(double*));
	for(i=0;i<3*nAtoms;i++) modes[i] = malloc(3*nAtoms*sizeof(double));

	//printf("begin diag\n");
	eigenQL(3*nAtoms, F, frequencies, modes);

	//printf("end eigneQL\n");
	for(i=0;i<mol->vibration.nModes;i++) mol->vibration.modes[i].properties[0] = 0.0;
	
	/* convert in atomic unit  from kcal/Ang^2/amu */
	for(i=0;i<3*nAtoms;i++) frequencies[i] *= 1.59360150e-03*0.529177*0.529177*5.48579911e-04; 
	/* convert frequencies in cm-1 */
	for(i=0;i<3*nAtoms;i++) 
		if(frequencies[i]>0) frequencies[i] = sqrt(frequencies[i])*219474.63633664;
		else frequencies[i] = -sqrt(-frequencies[i])*219474.63633664;

	/* compute the IR intensities */
	for(i=0;i<nAtoms;i++)
	for(k=0;k<3;k++)
	{
		int id=3*i+k;
		double IRI = 0;
		double D[3] = {0,0,0};
		int kp;
		for(c = 0;c<3;c++)
		for(j=0;j<nAtoms;j++)
		for(kp = 0;kp<3;kp++) 
		{
			int jd = 3*j+kp;
			double Lji = modes[jd][id];
			double a=dmuX[c][jd]*Lji/sqrt(mol->atoms[j].mass);
			D[c]+=a;
		}
		IRI = 0;
		for(c = 0;c<3;c++)  IRI+= D[c]*D[c];
		mol->vibration.modes[id].properties[0] = IRI;
	}
	/* Intensities in 1 (D/Ang)^2 amu^-1 = 42.255 km/mol=171.65 cm^-2 atm^-1 at 0 C and 1 atm */
	/* Refs : D. Porezag and M. R. Pederson, Phys. Rev. B 54, 7830 (1996). and Y. Yamaguchi el al., J. Chem. Phys. 84,2262(1986)*/
	/* conversion in km/mol*/
	for(i=0;i<mol->vibration.nModes;i++) mol->vibration.modes[i].properties[0] *= 42.255;

	/* compute the reduced mass */
	for(i=0;i<3*nAtoms;i++) 
	{
		double m = 0;
		for(j=0;j<mol->nAtoms;j++)
		{
			double r2 = 0;
			for(c=0;c<3;c++) r2+= modes[3*j+c][i]*modes[3*j+c][i];
			m+= r2/(mol->atoms[j].mass); 
		}
		if(m<=0) m = 1;
		m = 1/m;
		for(j=0;j<mol->nAtoms;j++)
		{
			double r =sqrt(m)/sqrt(mol->atoms[j].mass);
			for(c=0;c<3;c++) modes[3*j+c][i]*=r;
		}

		//printf("%f %f\n",(*frequencies)[i],m);
		mol->vibration.modes[i].mass = m;
	}
	for(i=0;i<mol->vibration.nModes;i++) 
	{
		for(c=0;c<3;c++) mol->vibration.modes[i].frequency = frequencies[i];
		for(j=0;j<mol->nAtoms;j++)
			for(c=0;c<3;c++) mol->vibration.modes[i].vectors[c][j] = modes[3*j+c][i];
	}
	if(frequencies) free(frequencies);
	if(modes) for(i=0;i<3*nAtoms;i++) if(modes[i]) free(modes[i]);
	if(modes) free(modes);
	sortFrequencies(mol);
	removeTransRotModes(mol);
	return 3*nAtoms;
}
/*****************************************************************************************************/
static int readEnergyDipoleGradFromGabeditFile(Molecule* mol, char* inputFileName, int index, boolean readGrad)
{
	char* fileName = NULL;
	char* prefixName = NULL;
	int i;

	if(!mol) return 1;
	if(mol->nAtoms<1) return 1;
	if(!inputFileName) return 1;

	prefixName = strdup_printf("%sFreq",getSuffixNameFile(inputFileName));

	fileName = strdup_printf("%s_%d.gab",prefixName,index);
	i = readEnergyAndDipoleFromGabeditFile(fileName, &mol->potentialEnergy, mol->dipole);
	printf("FileName = %s E = %f D = %f %f %f\n", fileName, mol->potentialEnergy, mol->dipole[0],  mol->dipole[1],  mol->dipole[2]);
	if(readGrad) mol->klass->readGradientFromGabeditFile(mol, fileName);

	if(fileName) free(fileName);
	return i;
}
/*****************************************************************************/
static int computeGradientsFromFiles(Molecule* mol, char* inputFileName, int indexBegin, double dx)
{
        int i;
        int k;
        int nAtoms;
	int index = indexBegin;
	double Ep, Em;

	if(!inputFileName) return indexBegin;
	nAtoms = mol->nAtoms;

        for(i=0;i<nAtoms;i++)
        for(k=0;k<3;k++)
        {
                mol->atoms[i].coordinates[k] += dx;
		readEnergyDipoleGradFromGabeditFile(mol, inputFileName, index, FALSE);
		Ep = mol->potentialEnergy;
		index++;

                mol->atoms[i].coordinates[k] -= 2*dx;
		readEnergyDipoleGradFromGabeditFile(mol, inputFileName, index, FALSE);
		Em = mol->potentialEnergy;
		index++;
		mol->atoms[i].gradient[k] = (Ep-Em)/dx/2;

                mol->atoms[i].coordinates[k] += dx;
        }
	return index;
}
/*****************************************************************************/
static int computeFrequenciesFromFiles(Molecule* mol, char* inputFileName, double dx)
{
	int i;
	int j;
	int k;
	int c;
	int id,jd,index,idx;
	double* F;
	double* gp[3];
	double* gm[3];
	double* dmuX[3];
	double Dp[3];
	double Dm[3];
	int nAtoms;
	int ret;

	if(!mol || mol->nAtoms<1) return 0;

	nAtoms = mol->nAtoms;
	for(k=0;k<3;k++) gp[k] = malloc(nAtoms*sizeof(double));
	for(k=0;k<3;k++) gm[k] = malloc(nAtoms*sizeof(double));
	for(k=0;k<3;k++) dmuX[k] = malloc(3*nAtoms*sizeof(double));

	F = malloc(3*nAtoms*(3*nAtoms+1)/2*sizeof(double));

	index = 0;
	for(i=0;i<nAtoms;i++)
	for(k=0;k<3;k++)
	{
		id=3*i+k;

		mol->atoms[i].coordinates[k] += dx;
		index = computeGradientsFromFiles(mol, inputFileName, index, dx);
		copyGradients(mol, gp);
		readEnergyDipoleGradFromGabeditFile(mol, inputFileName, index,FALSE);
		for(c = 0;c<3;c++)  Dp[c] = mol->dipole[c];
		index++;
		
		mol->atoms[i].coordinates[k] -= 2*dx;
		index = computeGradientsFromFiles(mol, inputFileName, index, dx);
		copyGradients(mol, gm);
		readEnergyDipoleGradFromGabeditFile(mol, inputFileName, index,FALSE);
		for(c = 0;c<3;c++)  Dm[c] = mol->dipole[c];
		index++;

		for(c = 0;c<3;c++) dmuX[c][id] = (Dp[c]-Dm[c])/dx/2;
		mol->atoms[i].coordinates[k] += dx;
		
		for(j=0;j<=i;j++)
		{
			double invm = 1.0/sqrt( mol->atoms[i].mass* mol->atoms[j].mass);
			for(c = 0;c<3;c++) 
			{
				jd = 3*j+c;
				//printf("id = %d jd = %d\n",id,jd);
				if(jd>id) continue;
				idx = jd + id*(id+1)/2;
				//printf("index = %d i = %d k = %d j = %d c = %d\n",index,i,k,j,c);
				F[idx] = (gp[c][j]-gm[c][j])/dx/2; 
				F[idx] *= invm;
			}
		}
	}
	//printf("F\n");
	/*
	for(id=0;id<3*nAtoms;id++)
	{
		for(jd=0;jd<=id;jd++) 
		{
			index = jd + id*(id+1)/2;
			printf("%14.8f",F[index]);
		}
		printf("\n");
	}
	*/

	for(k=0;k<3;k++) free(gp[k]);
	for(k=0;k<3;k++) free(gm[k]);
	ret = computeIR(mol, F, dmuX);
	free(F);
	for(k=0;k<3;k++) free(dmuX[k]);
	return ret;
}
/*******************************************************************************************************************/
static int computeGradientsOneStepFromFiles(Molecule* mol, char* inputFileName, int indexBegin, double D[], double dx)
{
        int i;
        int k;
        int nAtoms;
	int index = indexBegin;
	double Ep, E0;
	int c;

	if(!mol || mol->nAtoms<1) return indexBegin;
	if(!inputFileName) return indexBegin;

	nAtoms = mol->nAtoms;

	readEnergyDipoleGradFromGabeditFile(mol, inputFileName, index,FALSE);
	E0 = mol->potentialEnergy;
	for(c = 0;c<3;c++)  D[c] = mol->dipole[c];
	index++;
        for(i=0;i<nAtoms;i++)
        for(k=0;k<3;k++)
        {
                mol->atoms[i].coordinates[k] += dx;
		readEnergyDipoleGradFromGabeditFile(mol, inputFileName, index,FALSE);
		Ep = mol->potentialEnergy;
		index++;

		mol->atoms[i].gradient[k] = (Ep-E0)/dx;

                mol->atoms[i].coordinates[k] -= dx;
        }
	return index;
}
/*******************************************************************************************************************/
static int computeFrequenciesOneStepFromFiles(Molecule* mol, char* inputFileName, double dx)
{
	int i;
	int j;
	int k;
	int c;
	int id,jd,index,idx;
	double* F;
	double* gp[3];
	double* g0[3];
	double* dmuX[3];
	double Dp[3];
	double D0[3];
	int nAtoms;
	int ret;

	if(!mol || mol->nAtoms<1) return 0;
	nAtoms = mol->nAtoms;
	for(k=0;k<3;k++) gp[k] = malloc(nAtoms*sizeof(double));
	for(k=0;k<3;k++) g0[k] = malloc(nAtoms*sizeof(double));
	for(k=0;k<3;k++) dmuX[k] = malloc(3*nAtoms*sizeof(double));

	F = malloc(3*nAtoms*(3*nAtoms+1)/2*sizeof(double));

	index = 0;
	index = computeGradientsOneStepFromFiles(mol, inputFileName, index, D0, dx);
	copyGradients(mol, g0);

	for(i=0;i<nAtoms;i++)
	for(k=0;k<3;k++)
	{
		id=3*i+k;
		mol->atoms[i].coordinates[k] += dx;
		index = computeGradientsOneStepFromFiles(mol, inputFileName, index, Dp, dx);
		copyGradients(mol, gp);
		
		for(c = 0;c<3;c++) dmuX[c][id] = (Dp[c]-D0[c])/dx;
		mol->atoms[i].coordinates[k] -= dx;
		
		for(j=0;j<=i;j++)
		{
			double invm = 1.0/sqrt( mol->atoms[i].mass* mol->atoms[j].mass);
			for(c = 0;c<3;c++) 
			{
				jd = 3*j+c;
				//printf("id = %d jd = %d\n",id,jd);
				if(jd>id) continue;
				idx = jd + id*(id+1)/2;
				//printf("index = %d i = %d k = %d j = %d c = %d\n",index,i,k,j,c);
				F[idx] = (gp[c][j]-g0[c][j])/dx; 
				F[idx] *= invm;
			}
		}
	}
	//printf("F\n");
	/*
	for(id=0;id<3*nAtoms;id++)
	{
		for(jd=0;jd<=id;jd++) 
		{
			index = jd + id*(id+1)/2;
			printf("%14.8f",F[index]);
		}
		printf("\n");
	}
	*/

	for(k=0;k<3;k++) free(gp[k]);
	for(k=0;k<3;k++) free(g0[k]);
	ret = computeIR(mol, F, dmuX);
	free(F);
	for(k=0;k<3;k++) free(dmuX[k]);
	return ret;
}
/*****************************************************************************/
static int computeFrequenciesFromGradFiles(Molecule* mol, char* inputFileName, double dx)
{
	int i;
	int j;
	int k;
	int c;
	int id,jd,index,idx;
	double* F;
	double* gp[3];
	double* gm[3];
	double* dmuX[3];
	double Dp[3];
	double Dm[3];
	int nAtoms;
	int ret;

	if(!mol || mol->nAtoms<1) return 0;
	nAtoms = mol->nAtoms;
	for(k=0;k<3;k++) gp[k] = malloc(nAtoms*sizeof(double));
	for(k=0;k<3;k++) gm[k] = malloc(nAtoms*sizeof(double));
	for(k=0;k<3;k++) dmuX[k] = malloc(3*nAtoms*sizeof(double));

	F = malloc(3*nAtoms*(3*nAtoms+1)/2*sizeof(double));

	index = 0;
	for(i=0;i<nAtoms;i++)
	for(k=0;k<3;k++)
	{
		id=3*i+k;
		mol->atoms[i].coordinates[k] += dx;
		readEnergyDipoleGradFromGabeditFile(mol, inputFileName, index,TRUE);
		copyGradients(mol, gp);
		for(c = 0;c<3;c++)  Dp[c] = mol->dipole[c];
		index++;
		
		mol->atoms[i].coordinates[k] -= 2*dx;
		readEnergyDipoleGradFromGabeditFile(mol, inputFileName, index,TRUE);
		copyGradients(mol, gm);
		for(c = 0;c<3;c++)  Dm[c] = mol->dipole[c];
		index++;
		for(c = 0;c<3;c++) dmuX[c][id] = (Dp[c]-Dm[c])/dx/2;

		mol->atoms[i].coordinates[k] += dx;
		
		for(j=0;j<=i;j++)
		{
			double invm = 1.0/sqrt( mol->atoms[i].mass* mol->atoms[j].mass);
			for(c = 0;c<3;c++) 
			{
				jd = 3*j+c;
				//printf("id = %d jd = %d\n",id,jd);
				if(jd>id) continue;
				idx = jd + id*(id+1)/2;
				//printf("index = %d i = %d k = %d j = %d c = %d\n",index,i,k,j,c);
				F[idx] = (gp[c][j]-gm[c][j])/dx/2; 
				F[idx] *= invm;
			}
		}
	}
	//printf("F\n");
	/*
	for(id=0;id<3*nAtoms;id++)
	{
		for(jd=0;jd<=id;jd++) 
		{
			index = jd + id*(id+1)/2;
			printf("%14.8f",F[index]);
		}
		printf("\n");
	}
	*/

	for(k=0;k<3;k++) free(gp[k]);
	for(k=0;k<3;k++) free(gm[k]);
	ret = computeIR(mol, F, dmuX);
	free(F);
	for(k=0;k<3;k++) free(dmuX[k]);
	return ret;
}
/*******************************************************************************************************************/
static int computeFrequenciesOneStepFromGradFiles(Molecule* mol, char* inputFileName, double dx)
{
	int i;
	int j;
	int k;
	int c;
	int id,jd,index,idx;
	double* F;
	double* gp[3];
	double* g0[3];
	double* dmuX[3];
	double Dp[3];
	double D0[3];
	int nAtoms;
	int ret;

	if(!mol || mol->nAtoms<1) return 0;
	nAtoms = mol->nAtoms;
	for(k=0;k<3;k++) gp[k] = malloc(nAtoms*sizeof(double));
	for(k=0;k<3;k++) g0[k] = malloc(nAtoms*sizeof(double));
	for(k=0;k<3;k++) dmuX[k] = malloc(3*nAtoms*sizeof(double));

	F = malloc(3*nAtoms*(3*nAtoms+1)/2*sizeof(double));

	index = 0;
	readEnergyDipoleGradFromGabeditFile(mol, inputFileName, index,TRUE);
	copyGradients(mol, g0);
	for(c = 0;c<3;c++)  D0[c] = mol->dipole[c];
	index++;

	for(i=0;i<nAtoms;i++)
	for(k=0;k<3;k++)
	{
		id=3*i+k;
		mol->atoms[i].coordinates[k] += dx;
		readEnergyDipoleGradFromGabeditFile(mol, inputFileName, index,TRUE);
		copyGradients(mol, gp);
		for(c = 0;c<3;c++)  Dp[c] = mol->dipole[c];
		index++;
		
		for(c = 0;c<3;c++) dmuX[c][id] = (Dp[c]-D0[c])/dx;
		mol->atoms[i].coordinates[k] -= dx;
		
		for(j=0;j<=i;j++)
		{
			double invm = 1.0/sqrt( mol->atoms[i].mass* mol->atoms[j].mass);
			for(c = 0;c<3;c++) 
			{
				jd = 3*j+c;
				//printf("id = %d jd = %d\n",id,jd);
				if(jd>id) continue;
				idx = jd + id*(id+1)/2;
				//printf("index = %d i = %d k = %d j = %d c = %d\n",index,i,k,j,c);
				F[idx] = (gp[c][j]-g0[c][j])/dx; 
				F[idx] *= invm;
			}
		}
	}
	//printf("F\n");
	/*
	for(id=0;id<3*nAtoms;id++)
	{
		for(jd=0;jd<=id;jd++) 
		{
			index = jd + id*(id+1)/2;
			printf("%14.8f",F[index]);
		}
		printf("\n");
	}
	*/

	for(k=0;k<3;k++) free(gp[k]);
	for(k=0;k<3;k++) free(g0[k]);
	ret = computeIR(mol, F, dmuX);
	free(F);
	for(k=0;k<3;k++) free(dmuX[k]);
	return ret;
}
/*****************************************************************************/
static boolean addInputToCChemIFile(FILE* fileOut, char* inputFileName, char* runType)
{
	boolean Ok = TRUE;
	char buffer[1024];
	int i;
	FILE* fileIn = fopen(inputFileName,"r");
	if(!fileIn) return FALSE;
        while(!feof(fileIn))
        {
               	if(!fgets(buffer,BSIZE,fileIn)){ Ok = FALSE; break;}
		for(i=0;i<strlen(buffer);i++) if(buffer[i] != ' ') { break;}
		if(buffer[i]=='#') 
		{
			char* str = strdup(buffer);
			uppercase(str);
			if(!strstr(str,"GEOMETRY")) fprintf(fileOut,"%s",buffer);
			free(str);
		}
		else
		{
			char* str = strdup(buffer);
			sscanf(buffer,"%s",str);
			uppercase(str);
			if(!strcmp(str,"GEOMETRY"))
			{
				int nA = 0;
               			if(!fgets(buffer,BSIZE,fileIn)){ Ok = FALSE; break;}
				sscanf(buffer,"%d",&nA);
				for(i=0;i<nA;i++) if(!fgets(buffer,BSIZE,fileIn)) { Ok = FALSE; break;}
				if(!Ok) break;
				
			}
			else if(strstr(str,"RUNTYPE")) fprintf(fileOut,"%s%s\n","RunType=",runType);
			else fprintf(fileOut,"%s",buffer);
			free(str);
		}
	}
	fclose(fileIn);
	return Ok;
}
/*****************************************************************************/
static int createCChemIFile(Molecule* mol, char* inputFileName, int index, char* runType)
{
	char* fileName = NULL;
	FILE* file = NULL;
	char* prefixName = NULL;

	if(!mol || mol->nAtoms<1) return 1;
	if(!inputFileName) return 1;

	prefixName = strdup_printf("%sFreq",getSuffixNameFile(inputFileName));
	fileName = strdup_printf("%s_%d.ici",prefixName,index);
	file = fopen(fileName,"w");
	addInputToCChemIFile(file, inputFileName, runType);
	mol->klass->addGeometry(mol,file);
	fclose(file);

	if(fileName) free(fileName);
	return 0;
}
/******************************************************************************************************/
static int generateCChemIFilesGradients(Molecule* mol, char* inputFileName, int indexBegin, double dx)
{
        int i;
        int k;
        int nAtoms;
	int index = indexBegin;

	if(!mol || mol->nAtoms<1) return indexBegin;
	if(!inputFileName) return indexBegin;

	nAtoms = mol->nAtoms;

        for(i=0;i<nAtoms;i++)
        for(k=0;k<3;k++)
        {
                mol->atoms[i].coordinates[k] += dx;
		createCChemIFile(mol, inputFileName, index,"Energy");
		index++;

                mol->atoms[i].coordinates[k] -= 2*dx;
		createCChemIFile(mol, inputFileName, index,"Energy");
		index++;

                mol->atoms[i].coordinates[k] += dx;
        }
	return index;
}
/******************************************************************************************************/
static int generateCChemIFilesForFrequencies(Molecule* mol, char* inputFileName, double dx)
{
	int i;
	int k;
	int index;
	int nAtoms;

	if(!mol || mol->nAtoms<1) return 1;
	if(!inputFileName) return 1;

	nAtoms = mol->nAtoms;
	printf("nAtoms = %d\n",nAtoms);

	index = 0;
	for(i=0;i<nAtoms;i++)
	for(k=0;k<3;k++)
	{
		mol->atoms[i].coordinates[k] += dx;
		index = generateCChemIFilesGradients(mol, inputFileName, index, dx);
		createCChemIFile(mol, inputFileName, index,"Energy");
		index++;
		
		mol->atoms[i].coordinates[k] -= 2*dx;
		index = generateCChemIFilesGradients(mol, inputFileName, index, dx);
		createCChemIFile(mol, inputFileName, index,"Energy");
		index++;

		mol->atoms[i].coordinates[k] += dx;
	}
	printf("see %sFreq_*.ici generated files\n",inputFileName);
	return 0;
}
/******************************************************************************************************/
static int generateCChemIFilesOneStepGradients(Molecule* mol, char* inputFileName, int indexBegin, double dx)
{
        int i;
        int k;
        int nAtoms;
	int index = indexBegin;

	if(!mol || mol->nAtoms<1) return indexBegin;
	if(!inputFileName) return indexBegin;

	nAtoms = mol->nAtoms;

	createCChemIFile(mol, inputFileName, index,"Energy");
	index++;
        for(i=0;i<nAtoms;i++)
        for(k=0;k<3;k++)
        {
                mol->atoms[i].coordinates[k] += dx;
		createCChemIFile(mol, inputFileName, index,"Energy");
		index++;

                mol->atoms[i].coordinates[k] -= dx;
        }
	return index;
}
/******************************************************************************************************/
static int generateCChemIFilesOneStepForFrequencies(Molecule* mol, char* inputFileName, double dx)
{
	int i;
	int k;
	int index;
	int nAtoms;

	if(!mol || mol->nAtoms<1) return 1;
	if(!inputFileName) return 1;

	nAtoms = mol->nAtoms;
	printf("nAtoms = %d\n",nAtoms);

	index = 0;
	index = generateCChemIFilesOneStepGradients(mol, inputFileName, index, dx);

	for(i=0;i<nAtoms;i++)
	for(k=0;k<3;k++)
	{
		mol->atoms[i].coordinates[k] += dx;
		index = generateCChemIFilesOneStepGradients(mol, inputFileName, index, dx);
		
		mol->atoms[i].coordinates[k] -= dx;
	}
	printf("see %sFreq_*.ici generated files\n",inputFileName);
	return 0;
}
/******************************************************************************************************/
static int generateCChemIGradFilesForFrequencies(Molecule* mol, char* inputFileName, double dx)
{
	int i;
	int k;
	int index;
	int nAtoms;

	if(!mol || mol->nAtoms<1) return 1;
	if(!inputFileName) return 1;

	nAtoms = mol->nAtoms;
	printf("nAtoms = %d\n",nAtoms);

	index = 0;
	for(i=0;i<nAtoms;i++)
	for(k=0;k<3;k++)
	{
		mol->atoms[i].coordinates[k] += dx;
		createCChemIFile(mol, inputFileName, index,"Gradient");
		index++;
		
		mol->atoms[i].coordinates[k] -= 2*dx;
		createCChemIFile(mol, inputFileName, index,"Gradient");
		index++;

		mol->atoms[i].coordinates[k] += dx;
	}
	printf("see %sFreq_*.ici generated files\n",inputFileName);
	return 0;
}
/******************************************************************************************************/
static int generateCChemIGradFilesOneStepForFrequencies(Molecule* mol, char* inputFileName, double dx)
{
	int i;
	int k;
	int index;
	int nAtoms;

	if(!mol || mol->nAtoms<1) return 1;
	if(!inputFileName) return 1;

	nAtoms = mol->nAtoms;
	printf("nAtoms = %d\n",nAtoms);

	index = 0;
	createCChemIFile(mol, inputFileName, index,"Gradient");
	index++;

	for(i=0;i<nAtoms;i++)
	for(k=0;k<3;k++)
	{
		mol->atoms[i].coordinates[k] += dx;
		createCChemIFile(mol, inputFileName, index,"Gradient");
		index++;
		
		mol->atoms[i].coordinates[k] -= dx;
	}
	printf("see %sFreq_*.ici generated files\n",inputFileName);
	return 0;
}
/******************************************************************************************************/
static void computeFixedNormalModesSampling(Molecule* mol, double* quantumNumbers, boolean changeGeom)
{
	int i,k,j;
	int nAtoms;
	double* Q = NULL;
	double* dQ = NULL;
	int ntr = 6;
	int nModes = 0;

	if(!mol || mol->nAtoms<1 || mol->vibration.nModes<1 ) return;
	nAtoms = mol->nAtoms;
	nModes = mol->vibration.nModes;
	if(mol->klass->isLinear(mol)) ntr = 5;
	if(nModes != 3*nAtoms) ntr = 0;
	

	Q = malloc(nModes*sizeof(double));
	dQ = malloc(nModes*sizeof(double));
	for(i=0;i<nModes;i++) Q[i] = 0.0;
	for(i=0;i<nModes;i++) dQ[i] = 0.0;

	// Q and dQ in AU
	for(i=ntr;i<nModes;i++) 
	{
		double R = (normal()+1)/2.0;
		double fi = fabs(mol->vibration.modes[i].frequency)/AUTOCM1;
		double Ei = (quantumNumbers[i]+0.5)*mol->vibration.modes[i].frequency/AUTOCM1; 
		double Ai = 0;
		if(Ei<0) Ei = 0;
		if(fi>1e-10) Ai = sqrt(2.0*fabs(Ei))/fi;
		Q[i] = Ai*cos(2*PI*R); 
		dQ[i] = -fi*Ai*sin(2*PI*R); 
	}
	// Q and dQ in AU
	for(i=ntr;i<nModes;i++)  Q[i]  /= sqrt(mol->vibration.modes[i].mass*AMUTOAU);
	for(i=ntr;i<nModes;i++)  dQ[i] /= sqrt(mol->vibration.modes[i].mass*AMUTOAU);

	// Q in Ang
	for(i=ntr;i<nModes;i++)  Q[i] *= BOHRTOANG;
	// dQ in Ang/AKMA-time 
	for(i=ntr;i<nModes;i++)  dQ[i] *= BOHRTOANG/(AUTOfs*fsInAKMA);
	
	for(j=0;j<nAtoms;j++) 
	for(k=0;k<3;k++) mol->atoms[j].velocity[k] = 0.0;

	for(i=0;i<nModes;i++) 
	{
		for(j=0;j<nAtoms;j++) 
		for(k=0;k<3;k++) 
		{
	 		 mol->atoms[j].velocity[k] += mol->vibration.modes[i].vectors[k][j]*dQ[i];
			if(changeGeom) mol->atoms[j].coordinates[k] +=  mol->vibration.modes[i].vectors[k][j]*Q[i];
		}
	}
	free(Q);
	free(dQ);
}
/*****************************************************************************************************************************************************/
static void computeHarmonicVelocitiesCoordinates(Molecule* mol, double T, int numMode, boolean changeGeom)
{
	int i;
	double* quantumNumbers = NULL;
	int nModes = 0;

	if(!mol || mol->nAtoms<1 || mol->vibration.nModes<1) return;
	nModes = mol->vibration.nModes;
	if(T<=0) T = 300;
	if(numMode>nModes) 
	{
		fprintf(stderr,"Error in HarmonicVelocityModes keyword, the mode number must be <=i %d\n",nModes);
		return;
	}
	quantumNumbers = malloc(nModes*sizeof(double));
	for(i=0;i<nModes;i++) quantumNumbers[i] = 0.0;
	if(numMode>=0) quantumNumbers[numMode-1] = 1.0;

	computeFixedNormalModesSampling(mol, quantumNumbers, changeGeom);

	mol->klass->removeTranslationAndRotation(mol);
	mol->klass->resetConstraints(mol, mol->constraints);
	mol->klass->scaleVelocities(mol, T);

	free(quantumNumbers);
}
/*****************************************************************************************************************************************************/
static double* getDeltaTable(Molecule* mol, double delta, boolean reducedCoordinates)
{
	double* deltas = malloc(mol->vibration.nModes*sizeof(double));
	int j;
	if(!reducedCoordinates)
		for(j=0;j<mol->vibration.nModes;j++) deltas[j] = delta/BOHRTOANG;
	else
	{
		double conv = delta*sqrt(AUTOCM1/AMUTOAU);
		for(j=0;j<mol->vibration.nModes;j++) 
			deltas[j] = conv/sqrt(mol->vibration.modes[j].frequency*mol->vibration.modes[j].mass);
	}
	return deltas;
}
/********************************************************************************/
/*  See 
	Calculation of NMR and EPR parameters: theory and applications
	By Martin Kaupp, Michael Buhl, Vladimir G. Malkin
	Published by Wiley-VCH, 2004,    ISBN 3527307796, 9783527307791
	Page 163, Equation 10.39
   See also : Toyama et al, J. Mol. Spec. 13,193 (1964), Eq. 7
 */
/********************************************************************************/
static double* get_centrifuge_parameters(Molecule* mol)
{
	int i;
	int j;
	int mode;
	double I[3] = {0,0,0};
	double alpha = 0;
	double beta = 0;
	double* akOverI = NULL;
	double a = 0;
	Atom* atoms = NULL;
	int nAtoms;
	int nModes;

	if(!mol || mol->nAtoms<1 || mol->vibration.nModes<1) return NULL;
	atoms = mol->atoms;
	nAtoms = mol->nAtoms;
	nModes = mol->vibration.nModes;

	akOverI = malloc(nModes*sizeof(double));

	for(i=0;i<nAtoms;i++)
	{
		for(j=0;j<3;j++)
		{
			alpha = atoms[i].coordinates[(j+1)%3];
			beta  = atoms[i].coordinates[(j+2)%3];
			I[j] += (alpha*alpha+beta*beta)*atoms[i].mass;
		}
	}
	/* printf("Ix = %lf Iyy = %lf Izz = %lf\n",I[0],I[1],I[2]);*/

	for(mode = 0;mode<nModes;mode++)
	{
		akOverI[mode] = 0;
		for(j=0;j<3;j++)
		{
			if(fabs(I[j])<1e-8) continue;
			a = 0;
			for(i=0;i<nAtoms;i++)
			{
				a += 
				atoms[i].coordinates[(j+1)%3]*mol->vibration.modes[mode].vectors[(j+1)%3][i]
			       +atoms[i].coordinates[(j+2)%3]*mol->vibration.modes[mode].vectors[(j+2)%3][i];
				/*
				printf ("%lf %lf %lf %lf \n",
				atoms[i].coordinates[(j+1)%3],mol->vibration.modes[mode].vectors[(j+1)%3][i],
				atoms[i].coordinates[(j+2)%3],mol->vibration.modes[mode].vectors[(j+2)%3][i]
						);
						*/
			}
			a *=2;
			a *=sqrt(mol->vibration.modes[mode].mass);
			/* printf("k = %d ak = %lf\n",mode+1,a);*/
			akOverI[mode] += a/I[j];
		}
	}
	return akOverI;
}
/*****************************************************************************************************************************************************/
static int createQFFOneCChemIFile(Molecule* mol, char* inputFileName, int index, double*** geoms, int mode, double delta)
{
	char* fileName = NULL;
	FILE* file = NULL;
	char* prefixName = NULL;
	char ad = '+';
	int i,j;
	Molecule* cmol = NULL;

	if(!mol || mol->nAtoms<1) return 1;
	if(!inputFileName) return 1;
	cmol = mol->klass->copy(mol);
	if(!cmol || cmol->nAtoms<1) return 1;

	if(delta<0) ad = '-';

	prefixName = strdup_printf("%sQFF",getSuffixNameFile(inputFileName));
	fileName = strdup_printf("%s_%d.ici",prefixName,index);
	file = fopen(fileName,"w");
	addInputToCChemIFile(file, inputFileName, "Energy");

	if(mode>-1) fprintf(file,"\n#Mode: Freq= %0.12lf Mass= %0.12lf Q= Qeq %c %0.12lf\n\n",
				mol->vibration.modes[mode].frequency,
				mol->vibration.modes[mode].mass,
				ad,
				fabs(delta)
				);
	for(i=0;i<mol->nAtoms;i++)
	for(j=0;j<3;j++)
		cmol->atoms[i].coordinates[j] = geoms[index][i][j];

	cmol->klass->addGeometry(cmol,file);
	fclose(file);

	if(fileName) free(fileName);
	cmol->klass->free(cmol);
	return 0;
}
/*******************************************************************************************************************************/
static int generateQFFCChemIFiles(Molecule* mol, char* inputFileName, double delta, boolean reducedCoordinates, int ordre)
{
	int index =0;
	int j;
	double* deltas = NULL;
	int nGeoms = 0;
	double*** geoms = getGeomsQFF(mol, inputFileName, delta,  reducedCoordinates, ordre, &deltas, &nGeoms);
	index = 0;
	createQFFOneCChemIFile(mol, inputFileName, index, geoms, -1, 0);
	for(j=0;j<mol->vibration.nModes;j++)
	{
		index++; createQFFOneCChemIFile(mol, inputFileName, index, geoms, j,  3*deltas[j]);
		index++; createQFFOneCChemIFile(mol, inputFileName, index, geoms, j,  2*deltas[j]);
		index++; createQFFOneCChemIFile(mol, inputFileName, index, geoms, j,   deltas[j]);
		index++; createQFFOneCChemIFile(mol, inputFileName, index, geoms, j,  -deltas[j]);
		index++; createQFFOneCChemIFile(mol, inputFileName, index, geoms, j, -2*deltas[j]);
		index++; createQFFOneCChemIFile(mol, inputFileName, index, geoms, j, -3*deltas[j]);
	}
	index++;
	for(    ;index<nGeoms;index++) createQFFOneCChemIFile(mol, inputFileName, index, geoms, -1, 0);
	printf("see %sQFF_*.ici generated files\n",inputFileName);
	for(index=0;index<nGeoms;index++) 
	{
		int i;
		for(i=0;i<mol->nAtoms;i++) free(geoms[index][i]);
		free(geoms[index]);
	}
	{
        	char* fileNameOut = strdup_printf("%sModes.gab",getSuffixNameFile(inputFileName));
		mol->klass->save(mol, fileNameOut);
		if(fileNameOut) free(fileNameOut);
	}
	return 0;
}
/*****************************************************************************************************************************************************/
static int getQFFOneGeom(Molecule* mol, int mode1, int mode2, double delta1, double delta2, double akOverI1, double** geom)
{
	int i,j;

	if(!mol || mol->nAtoms<1) return 1;

	delta1 *= BOHRTOANG;
	delta2 *= BOHRTOANG;
	for(i=0;i<mol->nAtoms;i++)
	for(j=0;j<3;j++)
		geom[i][j] = mol->atoms[i].coordinates[j];

	for(i=0;i<mol->nAtoms;i++)
	for(j=0;j<3;j++)
	{
		if(mode1>-1) geom[i][j] += mol->vibration.modes[mode1].vectors[j][i]*delta1;
		if(mode2>-1) geom[i][j] += mol->vibration.modes[mode2].vectors[j][i]*delta2;
	}
	return 0;
}
/*****************************************************************************************************************************************************/
static int getQFFOneGeom3(Molecule* mol, int mode1, int mode2, int mode3, double delta1, double delta2, double delta3, double** geom)
{
	int i,j;

	if(!mol || mol->nAtoms<1) return 1;

	delta1 *= BOHRTOANG;
	delta2 *= BOHRTOANG;
	delta3 *= BOHRTOANG;
	for(i=0;i<mol->nAtoms;i++)
	for(j=0;j<3;j++)
		geom[i][j] = mol->atoms[i].coordinates[j];

	for(i=0;i<mol->nAtoms;i++)
	for(j=0;j<3;j++)
	{
		if(mode1>-1) geom[i][j] += mol->vibration.modes[mode1].vectors[j][i]*delta1;
		if(mode2>-1) geom[i][j] += mol->vibration.modes[mode2].vectors[j][i]*delta2;
		if(mode3>-1) geom[i][j] += mol->vibration.modes[mode3].vectors[j][i]*delta3;
	}
	return 0;
}
/*******************************************************************************************************************************/
static int getQFFOneGeom4(Molecule* mol, int mode1, int mode2, int mode3, int mode4, double delta1, double delta2, double delta3, double delta4, double** geom)
{
	int i,j;

	if(!mol || mol->nAtoms<1) return 1;

	delta1 *= BOHRTOANG;
	delta2 *= BOHRTOANG;
	delta3 *= BOHRTOANG;
	delta4 *= BOHRTOANG;
	for(i=0;i<mol->nAtoms;i++)
	for(j=0;j<3;j++)
		geom[i][j] = mol->atoms[i].coordinates[j];

	for(i=0;i<mol->nAtoms;i++)
	for(j=0;j<3;j++)
	{
		if(mode1>-1) geom[i][j] += mol->vibration.modes[mode1].vectors[j][i]*delta1;
		if(mode2>-1) geom[i][j] += mol->vibration.modes[mode2].vectors[j][i]*delta2;
		if(mode3>-1) geom[i][j] += mol->vibration.modes[mode3].vectors[j][i]*delta3;
		if(mode4>-1) geom[i][j] += mol->vibration.modes[mode4].vectors[j][i]*delta4;
	}
	return 0;
}
/*****************************************************************************************************************************************************/
static double*** getGeomsQFF(Molecule* mol, char* inputFileName, double delta, boolean reducedCoordinates, int ordre, double** pDeltas, int* nGeoms)
{
	double*** geoms;
	int index =0;
	double* akOverI = NULL;
	int nAll = 0;
	double* deltas = NULL;
	int i,j,k,l;
	akOverI = get_centrifuge_parameters(mol);
	if(!akOverI) return NULL;
	deltas = getDeltaTable(mol, delta, reducedCoordinates);
	if(!deltas) return NULL;
	nAll = 1 +6* mol->vibration.nModes;
	if(ordre>=2) nAll += 6*mol->vibration.nModes*(mol->vibration.nModes-1); 
	if(ordre>=3) nAll += 8*mol->vibration.nModes*(mol->vibration.nModes-1)*(mol->vibration.nModes-2)/6;
	if(ordre>=4) nAll += 16*mol->vibration.nModes*(mol->vibration.nModes-1)*(mol->vibration.nModes-2)*(mol->vibration.nModes-3)/24;

	geoms = malloc(nAll*sizeof(double**));
	for(k=0;k<nAll;k++)
	{
		geoms[k] =  malloc(mol->nAtoms*sizeof(double*));
	       for(i=0;i<mol->nAtoms;i++) geoms[k][i] =  malloc(3*sizeof(double));
	}

	index = 0;
	getQFFOneGeom(mol, -1, -1, 0.0, 0.0, 0.0,geoms[index]);
	for(j=0;j<mol->vibration.nModes;j++)
	{
		index++; getQFFOneGeom(mol,  j, -1,   3*deltas[j], 0,akOverI[j],geoms[index]);
		index++; getQFFOneGeom(mol,  j, -1,   2*deltas[j], 0,akOverI[j],geoms[index]);
		index++; getQFFOneGeom(mol,  j, -1,     deltas[j], 0,akOverI[j],geoms[index]);
		index++; getQFFOneGeom(mol,  j, -1,    -deltas[j], 0,akOverI[j],geoms[index]);
		index++; getQFFOneGeom(mol,  j, -1,  -2*deltas[j], 0,akOverI[j],geoms[index]);
		index++; getQFFOneGeom(mol,  j, -1,  -3*deltas[j], 0,akOverI[j],geoms[index]);
	}
	if(akOverI) free(akOverI);
	if(ordre>=2)
	for(j=0;j<mol->vibration.nModes;j++)
	{
		for(i=0;i<j;i++)
		{
			index++; getQFFOneGeom(mol, j, i,  deltas[j],   deltas[i],0, geoms[index]);
			index++; getQFFOneGeom(mol, j, i,  deltas[j],  -deltas[i],0, geoms[index]);
			index++; getQFFOneGeom(mol, j, i, -deltas[j],   deltas[i],0, geoms[index]);
			index++; getQFFOneGeom(mol, j, i, -deltas[j],  -deltas[i],0, geoms[index]);
		}
		for(i=0;i<mol->vibration.nModes;i++)
		{
			if(i==j) continue;

			index++; getQFFOneGeom(mol,  j, i,  deltas[j],  3*deltas[i],0, geoms[index]);
			index++; getQFFOneGeom(mol,  j, i,  deltas[j], -3*deltas[i],0, geoms[index]);
			index++; getQFFOneGeom(mol,  j, i, -deltas[j],  3*deltas[i],0, geoms[index]);
			index++; getQFFOneGeom(mol,  j, i, -deltas[j], -3*deltas[i],0, geoms[index]);
		}
	}
	if(ordre>=3)
	for(j=0;j<mol->vibration.nModes;j++)
	{
		for(i=0;i<j;i++)
		for(k=0;k<i;k++)
		{
			/* 0  deltas[j],   deltas[i],  deltas[k] */
			index++; getQFFOneGeom3(mol, j, i, k,  deltas[j], deltas[i], deltas[k], geoms[index]);
			/* 1  deltas[j],   deltas[i], -deltas[k] */
			index++; getQFFOneGeom3(mol, j, i, k,  deltas[j], deltas[i], -deltas[k], geoms[index]);
			/* 2  deltas[j],  -deltas[i],  deltas[k] */
			index++; getQFFOneGeom3(mol, j, i, k,  deltas[j], -deltas[i], deltas[k], geoms[index]);
			/* 3  deltas[j],  -deltas[i], -deltas[k] */
			index++; getQFFOneGeom3(mol, j, i, k,  deltas[j], -deltas[i], -deltas[k], geoms[index]);
			/* 4 -deltas[j],  deltas[i],  deltas[k] */
			index++; getQFFOneGeom3(mol, j, i, k,  -deltas[j], deltas[i], deltas[k], geoms[index]);
			/* 5 -deltas[j],  deltas[i], -deltas[k] */
			index++; getQFFOneGeom3(mol, j, i, k,  -deltas[j], deltas[i], -deltas[k], geoms[index]);
			/* 6 -deltas[j], -deltas[i],  deltas[k] */
			index++; getQFFOneGeom3(mol, j, i, k,  -deltas[j], -deltas[i], deltas[k], geoms[index]);
			/* 7 -deltas[j], -deltas[i], -deltas[k]- */
			index++; getQFFOneGeom3(mol, j, i, k,  -deltas[j], -deltas[i], -deltas[k], geoms[index]);
		}
	}
	if(ordre>=4)
	for(j=0;j<mol->vibration.nModes;j++)
	{
		for(i=0;i<j;i++)
		for(k=0;k<i;k++)
		for(l=0;l<k;l++)
		{
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

			index++; getQFFOneGeom4(mol, j, i, k, l,  deltas[j], deltas[i], deltas[k], deltas[l], geoms[index]);
			index++; getQFFOneGeom4(mol, j, i, k, l,  deltas[j], deltas[i], deltas[k], -deltas[l], geoms[index]);
			index++; getQFFOneGeom4(mol, j, i, k, l,  deltas[j], deltas[i], -deltas[k], deltas[l], geoms[index]);
			index++; getQFFOneGeom4(mol, j, i, k, l,  deltas[j], deltas[i], -deltas[k], -deltas[l], geoms[index]);
			index++; getQFFOneGeom4(mol, j, i, k, l,  deltas[j], -deltas[i], deltas[k], deltas[l], geoms[index]);
			index++; getQFFOneGeom4(mol, j, i, k, l,  deltas[j], -deltas[i], deltas[k], -deltas[l], geoms[index]);
			index++; getQFFOneGeom4(mol, j, i, k, l,  deltas[j], -deltas[i], -deltas[k], deltas[l], geoms[index]);
			index++; getQFFOneGeom4(mol, j, i, k, l,  deltas[j], -deltas[i], -deltas[k], -deltas[l], geoms[index]);

			index++; getQFFOneGeom4(mol, j, i, k, l, -deltas[j], deltas[i], deltas[k], deltas[l], geoms[index]);
			index++; getQFFOneGeom4(mol, j, i, k, l, -deltas[j], deltas[i], deltas[k], -deltas[l], geoms[index]);
			index++; getQFFOneGeom4(mol, j, i, k, l, -deltas[j], deltas[i], -deltas[k], deltas[l], geoms[index]);
			index++; getQFFOneGeom4(mol, j, i, k, l, -deltas[j], deltas[i], -deltas[k], -deltas[l], geoms[index]);
			index++; getQFFOneGeom4(mol, j, i, k, l, -deltas[j], -deltas[i], deltas[k], deltas[l], geoms[index]);
			index++; getQFFOneGeom4(mol, j, i, k, l, -deltas[j], -deltas[i], deltas[k], -deltas[l], geoms[index]);
			index++; getQFFOneGeom4(mol, j, i, k, l, -deltas[j], -deltas[i], -deltas[k], deltas[l], geoms[index]);
			index++; getQFFOneGeom4(mol, j, i, k, l, -deltas[j], -deltas[i], -deltas[k], -deltas[l], geoms[index]);
		}
	}
	*pDeltas = deltas;
	*nGeoms = nAll;
	//printf("index = %d\n",index);
	//printf("nAll = %d\n",nAll);
	return geoms;
}
#define EPSILON 1.0E-12
#define SQU(x,y,z) ((x)*(x) + (y)*(y) + (z)*(z))

/************************************************************************************************************/
/* Jacobi diagonalisation of 3x3 symmetric matrix */
/* matrix mat stored like   0 3 5    
                              1 4
                                2   */
static void jacobi(double *mat, double evec[3][3])
{
  
	double t,s,u;
	double a;
	evec[0][1] = evec[0][2] = evec[1][0] = 0.0;  /* unity matrix */
	evec[1][2] = evec[2][0] = evec[2][1] = 0.0;  /* unity matrix */
	evec[0][0] = evec[1][1] = evec[2][2] = 1.0;

	/* do jacobi sweep */
	while(SQU(mat[3],mat[4],mat[5]) > EPSILON)
	{
		/* set mat[3] to zero */
		if (mat[3]*mat[3] > EPSILON)
		{
			t = (mat[1]-mat[0])/(2.*mat[3]);
	 		t = (t>0) ? 1./(t+sqrt(t*t+1.)) : -1./(-t+sqrt(t*t+1.));
	 		s = t/(sqrt(t*t+1.));
	  		u = s*t/(s+t);

	  		mat[0] -= t*mat[3];
	  		mat[1] += t*mat[3];
	  		a = mat[5];
	  		mat[5] -= s*(mat[4]+u*mat[5]);
	  		mat[4] += s*(  a   -u*mat[4]);
	  		mat[3] = 0.;
	  
	 		a = evec[0][0];
	  		evec[0][0] -= s*(evec[0][1]+u*evec[0][0]);
	  		evec[0][1] += s*(   a   -u*evec[0][1]);

	  		a = evec[1][0];
	  		evec[1][0] -= s*(evec[1][1]+u*evec[1][0]);
	  		evec[1][1] += s*(   a   -u*evec[1][1]);

	 		a = evec[2][0];
	  		evec[2][0] -= s*(evec[2][1]+u*evec[2][0]);
	  		evec[2][1] += s*(   a   -u*evec[2][1]);
		}

		/* set mat[5] to zero */
		if (mat[5]*mat[5] > EPSILON)
		{
			t = (mat[2]-mat[0])/(2.*mat[5]);
	  		t = (t>0) ? 1./(t+sqrt(t*t+1.)) : -1./(-t+sqrt(t*t+1.));
	  		s = t/(sqrt(t*t+1.));
	  		u = s*t/(s+t);

	  		mat[0] -= t*mat[5];
	  		mat[2] += t*mat[5];
	  		a = mat[3];
	  		mat[3] -= s*(mat[4]+u*mat[3]);
	  		mat[4] += s*(  a   -u*mat[4]);
	  		mat[5] = 0.;
	  
	  		a = evec[0][0];
	  		evec[0][0] -= s*(evec[0][2]+u*evec[0][0]);
	  		evec[0][2] += s*(   a   -u*evec[0][2]);

	  		a = evec[1][0];
	  		evec[1][0] -= s*(evec[1][2]+u*evec[1][0]);
	  		evec[1][2] += s*(   a   -u*evec[1][2]);

	  		a = evec[2][0];
	  		evec[2][0] -= s*(evec[2][2]+u*evec[2][0]);
	  		evec[2][2] += s*(   a   -u*evec[2][2]);
		}

		/* set mat[4] to zero */
		if (mat[4]*mat[4] > EPSILON)
		{
	  		t = (mat[2]-mat[1])/(2.*mat[4]);
	  		t = (t>0) ? 1./(t+sqrt(t*t+1.)) : -1./(-t+sqrt(t*t+1.));
	  		s = t/(sqrt(t*t+1.));
	  		u = s*t/(s+t);

	  		mat[1] -= t*mat[4];
	  		mat[2] += t*mat[4];
	  		a = mat[3];
	  		mat[3] -= s*(mat[5]+u*mat[3]);
	  		mat[5] += s*(  a   -u*mat[5]);
	  		mat[4] = 0.;

	  		a = evec[0][1];
	  		evec[0][1] -= s*(evec[0][2]+u*evec[0][1]);
	  		evec[0][2] += s*(   a   -u*evec[0][2]);

	  		a = evec[1][1];
	  		evec[1][1] -= s*(evec[1][2]+u*evec[1][1]);
	  		evec[1][2] += s*(   a   -u*evec[1][2]);

	  		a = evec[2][1];
	  		evec[2][1] -= s*(evec[2][2]+u*evec[2][1]);
	  		evec[2][2] += s*(   a   -u*evec[2][2]);
		}
	}
}
/************************************************************************************************************/
/************************************************************************************************************/
static void swap(int i,int j,double* mat, double vecs[3][3])
{
	int k;
	double t;
	t = mat[i];
	mat[i] = mat[j];
	mat[j] = t;
	for(k=0;k<3;k++)
	{
		t = vecs[k][i];
		vecs[k][i] = vecs[k][j];
		vecs[k][j] = t;
	}
}
/************************************************************************************************************/
/* inertial moment of a molecule */
/* matrix mat stored like   0 3 5    
                              1 4
                                2   */
static void buildStandardOrientation(Molecule* mol,double* centerOfGravity, int* numberOfEquivalentAxes, double* inertialMoment, double axes[3][3])
{
	int i,j;
	double mat[6];
	double m;
	double mtot = 0.0;
	double x,y,z;
	double principalAxisTolerance = 1e-6;
	x = y = z =0.0;
	for (i=0;i<mol->nAtoms;i++)	  /* center of gravity and total mass */
	{
		/* m = sqrt(prime[atomList->type]);*/
		m = fabs(mol->atoms[i].mass);
		x += m* mol->atoms[i].coordinates[0];
		y += m* mol->atoms[i].coordinates[1];
		z += m* mol->atoms[i].coordinates[2];
		mtot += m;
	  }
	centerOfGravity[0] = x/mtot;
	centerOfGravity[1] = y/mtot;
	centerOfGravity[2] = z/mtot;
	
	for(i=0;i<6;i++) mat[i]=0.0;

	for (i=0;i<mol->nAtoms;i++)	  /* build up inertial tensor */
	{
		x = (mol->atoms[i].coordinates[0] -= centerOfGravity[0]);
		y = (mol->atoms[i].coordinates[1] -= centerOfGravity[1]);
		z = (mol->atoms[i].coordinates[2] -= centerOfGravity[2]);
		/* m = sqrt(prime[atomList->type]);*/
		m = fabs(mol->atoms[i].mass);
		mat[0] += m*(y*y+z*z);
		mat[1] += m*(x*x+z*z);
		mat[2] += m*(x*x+y*y);
		mat[3] -= m*(x*y);
		mat[4] -= m*(y*z);
		mat[5] -= m*(x*z);
	}

	jacobi(mat,axes);	  /* diagonalize tensor */

	/* sort eigenvalues */
	if (mat[0]<mat[1]) swap(0,1,mat,axes);
	if (mat[1]<mat[2]) swap(1,2,mat,axes);
	if (mat[0]<mat[1]) swap(0,1,mat,axes);

	inertialMoment[0] = mat[0];
	inertialMoment[1] = mat[1];
	inertialMoment[2] = mat[2];

	/* normalize moments if not pointlike */
	if (mat[0] > 1.E-8) 
	{
		mat[1] /= mat[0];
		mat[2] /= mat[0];
		mat[0] = 1.0;
	}

	if ((mat[1]-mat[2])*(mat[1]-mat[2])< principalAxisTolerance*principalAxisTolerance) swap(0,2,mat,axes); 

	/* determin number of equivalent axes */
	*numberOfEquivalentAxes = 1;  
	if ((mat[0]-mat[1])*(mat[0]-mat[1]) < principalAxisTolerance*principalAxisTolerance) (*numberOfEquivalentAxes)++; /* 2 axes equiv. */

	if ((mat[0]-mat[2])*(mat[0]-mat[2]) < principalAxisTolerance*principalAxisTolerance) (*numberOfEquivalentAxes)++; /* 3 axes equiv. */
	if (mat[2] < principalAxisTolerance) *numberOfEquivalentAxes = -*numberOfEquivalentAxes; /* linear or point */

	/* multiply atom vectors v_i by eigenvector matrix v_i' = v_i * V
	   to rotate molecule - principal axes will be equivalent to coordinate
	   axes
	 */
	for (i=0;i<mol->nAtoms;i++)	  /* perform rotation */
	{
		x = mol->atoms[i].coordinates[0];
		y = mol->atoms[i].coordinates[1];
		for (j=0;j<3;j++)
			mol->atoms[i].coordinates[j] = x*axes[0][j] + y*axes[1][j] + mol->atoms[i].coordinates[2]*axes[2][j];
	}
}
/*************************************************************************************************************************/
static char* saveFirstDerivatives(char* inputFileName, Molecule* mol)
{
	if(mol && mol->vibration.nModes>0 && mol->vibration.nProperties>=5)
	{
        	char* fileNameOut = strdup_printf("%sFirstDerivatives.txt",getSuffixNameFile(inputFileName));
        	FILE* file = fopen(fileNameOut,"w");
        	fprintf(stdout,"First derivatives saved in %s file\n", fileNameOut);
		mol->klass->addFirstDerivativeToFile(mol, file);
        	fclose(file);
		return fileNameOut;
	}
	return NULL;

}
/*************************************************************************************************************************/
static void computePseudoInertia(Molecule* mol, double*pI, double* pI4)
{
	int nAtoms = mol->nAtoms;
	if(nAtoms<1)
	{
		*pI = 0;
		*pI4 = 0;
		return;
	}
	int nTypes;
	char** types = getTypesOfAtoms(mol, &nTypes);
	int* m = malloc(nAtoms*sizeof(int));
	int* m4 = malloc(nAtoms*sizeof(int));
	double C[3];
	int i,j;

        for (j=0;j<3;j++) C[j] = 0;
        for (i=0;i<mol->nAtoms;i++)
	{
		m[i]=1;
		int it;
        	for (it=0;it<nTypes;it++)
		if(!strcmp(mol->atoms[i].prop.symbol,types[it]))
		{
			m[i]=it+1;
			m4[i]=nTypes-it;
			break;
		}
	}
	freeTypesOfAtoms(types,nTypes);

        for (i=0;i<mol->nAtoms;i++) for (j=0;j<3;j++) C[j] += mol->atoms[i].mass*mol->atoms[i].coordinates[j];
	double mt =0;
        for (i=0;i<mol->nAtoms;i++) mt += mol->atoms[i].mass;
	if(mt>0) for (j=0;j<3;j++) C[j] /= mt;

	double mt2 =0;
        for (i=0;i<mol->nAtoms;i++) mt2 += m[i];
	double mt4 =0;
        for (i=0;i<mol->nAtoms;i++) mt4 += m4[i];

	double I=0;
	double I4=0;
        for (i=0;i<mol->nAtoms;i++)
        {
		double d =0;
                for (j=0;j<3;j++)
                {
                	double dif= mol->atoms[i].coordinates[j]-C[j];
                	d += dif*dif;
                }
		I+= m[i]*d;
		I4+= m4[i]*d;
        }
	if(mt2>0) I /= mt2;
	if(mt4>0) I4 /= mt4;
	*pI = I;
	*pI4 = I4;
        if(m) free(m);
        if(m4) free(m4);
}
/****************************************************************************************************************************/
/* Jijun Zhao, Ruili Shi, Linwei Sai, Xiaoming Huang & Yan Su (2016)
 * Comprehensive genetic algorithm for ab initio global optimisation of clusters, Molecular Simulation,
 * 42:10, 809-819, DOI: 10.1080/08927022.2015.1121386
*/
static  boolean similarInertia(Molecule* mol1, Molecule* mol2, double precision)
{

	double I1,I2;
	double I41,I42;
	computePseudoInertia(mol1, &I1, &I41);
	computePseudoInertia(mol2, &I2, &I42);
	if(fabs(I2-I1)<precision && fabs(I42-I41)<precision) return TRUE;
	return FALSE;
}
