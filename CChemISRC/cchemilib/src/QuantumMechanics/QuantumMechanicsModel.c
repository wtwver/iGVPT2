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

/* QuantumMechanicsModel.c */
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>

#include "../Utils/Utils.h"
#include "../Utils/AtomsProp.h"
#include "../Utils/Constants.h"
#include "../Utils/Types.h"
#include "../Utils/QL.h"
#include "../QuantumMechanics/QuantumMechanicsModel.h"

static void freeQuantumMechanicsModel(QuantumMechanicsModel* qmModel);
static int computeQMFrequenciesNumeric(QuantumMechanicsModel* qmModel, double** frequencies, double*** modes, double** reducedMasses, double** IRIntensities);
static int computeQMFrequenciesAnalytic(QuantumMechanicsModel* qmModel, double** frequencies, double*** modes, double** reducedMasses, double** IRIntensities);
static int computeQMFrequencies(QuantumMechanicsModel* qmModel, double** frequencies, double*** modes, double** reducedMasses, double** IRIntensities);
static QuantumMechanicsModel copyQuantumMechanicsModel(QuantumMechanicsModel* f);
static void setRattleConstraintsParameters(QuantumMechanicsModel* quantumMechanicsModel);
static void removeFragmentedMolecules(QuantumMechanicsModel** qmModels, int* pnumberOfGeometries, double* energies, FILE* logfile);
static void removeSimilarInertiaGeometries(QuantumMechanicsModel** qmModels, int* pnumberOfGeometries, double* energies, FILE* logfile, double tol);
static void removeSimilarBondsGeometries(QuantumMechanicsModel** qmModels, int* pnumberOfGeometries, double* energies, FILE* logfile, double sTol, double distTol);
static void removeSmallDistanceMolecules(QuantumMechanicsModel** qmModels, int* pnumberOfGeometries, double* energies, FILE* logfile);
static void sortByInertia(QuantumMechanicsModel** qmModels, int* pnumberOfGeometries, double* energies);
static void cutByInertia(QuantumMechanicsModel** qmModels, int* pnumberOfGeometries, double* energies, int newNGeoms, FILE* logfile);
static QuantumMechanicsModel** getQuantumMechanicsRDConfo(QuantumMechanicsModel* qmModel, int numberOfGeometries, boolean chain, boolean saveFirstGeom);
static QuantumMechanicsModel** getQuantumMechanicsRDFConfo(QuantumMechanicsModel* qmModel, int numberOfGeometries, boolean saveFirstGeom);
static int computeIR(QuantumMechanicsModel* qmModel, double *F, double** dmu, double** frequencies, double*** modes, double** reducedMasses, double** IRIntensities);
/***********************************************************************************************************************/
static QuantumMechanicsModel** getQuantumMechanicsRDConfo(QuantumMechanicsModel* qmModel, int numberOfGeometries, boolean chain, boolean saveFirstGeom)
{
	int i;
	char* str = NULL;
	QuantumMechanicsModel** geometries = NULL;

	if(qmModel->molecule.nAtoms<1) return NULL;
	if(numberOfGeometries<2) return NULL;
	geometries = malloc(numberOfGeometries*sizeof(QuantumMechanicsModel*));
	fflush(stdout);
	for (i = 0; i < numberOfGeometries; i++ )
	{
		geometries[i] = NULL;
		if(i>0 || !saveFirstGeom) 
		{
			if(chain) qmModel->molecule.klass->setRandomPositionsChain(&qmModel->molecule);
			else qmModel->molecule.klass->setRandomPositions(&qmModel->molecule);
		}
		qmModel->klass->calculateEnergy(qmModel);
		if(str) free(str);
		str = strdup_printf(("Geometry # %d Potential energy \t=  %0.4f"), i+1, qmModel->molecule.potentialEnergy);
		printf("%s ",str);
		geometries[i] = malloc(sizeof(QuantumMechanicsModel));
		*geometries[i] = qmModel->klass->copy(qmModel);

		double I2,I4;
		geometries[i]->molecule.klass->computePseudoInertia(&geometries[i]->molecule, &I2, &I4);
		str = strdup_printf((" \tI2=%0.4f \tI4=%0.4f"), I2,I4);
		printf("%s\n",str);
		fflush(stdout);
	}
	if(str) free(str);
	return geometries;
}
/***********************************************************************************************************************/
static QuantumMechanicsModel** getQuantumMechanicsRDFConfo(QuantumMechanicsModel* qmModel, int numberOfGeometries, boolean saveFirstGeom)
{
	int i;
	char* str = NULL;
	QuantumMechanicsModel** geometries = NULL;

	if(qmModel->molecule.nAtoms<1) return NULL;
	if(numberOfGeometries<2) return NULL;
	geometries = malloc(numberOfGeometries*sizeof(QuantumMechanicsModel*));
	fflush(stdout);
	for (i = 0; i < numberOfGeometries; i++ )
	{
		geometries[i] = NULL;
		geometries[i] = malloc(sizeof(QuantumMechanicsModel));
		*geometries[i] = qmModel->klass->copy(qmModel);
		if(i>0 || !saveFirstGeom) 
		{
			qmModel->molecule.klass->setRandomFragments(&geometries[i]->molecule);
		}
		qmModel->klass->calculateEnergy(geometries[i]);
		if(str) free(str);
		str = strdup_printf(("Geometry # %d Potential energy \t=  %0.4f"), i+1, geometries[i]->molecule.potentialEnergy);
		printf("%s ",str);

		double I2,I4;
		geometries[i]->molecule.klass->computePseudoInertia(&geometries[i]->molecule, &I2, &I4);
		str = strdup_printf((" \tI2=%0.4f \tI4=%0.4f"), I2,I4);
		printf("%s\n",str);
		fflush(stdout);
	}
	if(str) free(str);
	return geometries;
}
/********************************************************************************************************************************************************/
static void sortByInertia(QuantumMechanicsModel** qmModels, int* pnumberOfGeometries, double* energies)
{
	int i,j;
	int numberOfGeometries = *pnumberOfGeometries;
	double* s = NULL;
	Molecule* mol;
	double I2,I4;
	double t;
	QuantumMechanicsModel* qm;
	double smax=-1;
	int nG=0; 

	if(numberOfGeometries<2) return;
	s = malloc(numberOfGeometries*sizeof(double));

	smax=-1;
	for (i = 0; i < numberOfGeometries; i++ ) 
	{
			s[i]= -1;
			if(!qmModels[i]) continue;
			mol = &qmModels[i]->molecule;
			mol->klass->computePseudoInertia(mol, &I2, &I4);
			s[i]= I4+I2;
			if(smax<0 || smax<s[i]) smax=s[i];
	}
	//fprintf(stderr,"End calcul s, smax=%f\n",smax); fflush(stderr);
	
	for (i = 0; i < numberOfGeometries; i++ ) if(!qmModels[i]) s[i]= smax+1e10;

	//if(qmModels && qmModels[0] && qmModels[0]->molecule.nAtoms<1) return;
	for (i = 0; i < numberOfGeometries; i++ )
	{
		int k=i;
		for (j = i+1; j < numberOfGeometries; j++ ) if(s[j]<s[k]) k=j;
		if(k!=i)
		{
			t = s[i];
			s[i]=s[k];
			s[k]=t;

			if(energies)
			{
				t = energies[i];
				energies[i]=energies[k];
				energies[k]=t;
			}

			qm=qmModels[i];
			qmModels[i]=qmModels[k];
			qmModels[k]=qm;
		}
	}

	nG=0; 
	for (i = 0; i < numberOfGeometries; i++ ) if(qmModels[i])  nG++;

	free(s);
	*pnumberOfGeometries = nG;
}
/********************************************************************************************************************************************************/
static void cutByInertia(QuantumMechanicsModel** qmModels, int* pnumberOfGeometries, double* energies, int newNGeoms, FILE* logfile)
{
	int i,j;
	int numberOfGeometries = *pnumberOfGeometries;
	int nG=0; 
	int* toRemove = NULL;
	int iBegin;

	if(numberOfGeometries<newNGeoms) return;
	//fprintf(stderr,"begin sorting\n");
	sortByInertia(qmModels, &numberOfGeometries, energies);
	//fprintf(stderr,"end sorting, nG=%d\n",numberOfGeometries); fflush(stderr);

	*pnumberOfGeometries = numberOfGeometries;
	if(numberOfGeometries<=newNGeoms) return;
	int nStep=numberOfGeometries/newNGeoms+1;

	//fprintf(stderr,"nStep=%d\n",nStep); fflush(stderr);
	//fprintf(stderr,"newNGeoms=%d\n",newNGeoms); fflush(stderr);

	toRemove = malloc(numberOfGeometries*sizeof(int));
	for (i = 0; i < numberOfGeometries; i++ ) toRemove[i]=1;

	nG=0;
	for (i = 0; i < numberOfGeometries; i+=nStep)
	{
		toRemove[i]=0;
		nG++;
		if(nG==newNGeoms) break;
	}
	//fprintf(stderr,"nG avant while=%d\n",nG); fflush(stderr);

	iBegin=0;
	while(nG<newNGeoms)
	{
		iBegin++;
		for (i = iBegin; i < numberOfGeometries; i+=nStep)
		{
			if(toRemove[i]==0) continue;
			toRemove[i]=0;
			nG++;
			if(nG==newNGeoms) break;
		}
	}
	//fprintf(stderr,"nG=%d\n",nG); fflush(stderr);

	//if(qmModels && qmModels[0] && qmModels[0]->molecule.nAtoms<1) return;
	for (i = 0; i < numberOfGeometries; i++ )
	{
		if(toRemove[i]==1 && qmModels[i]) 
		{
			qmModels[i]->klass->free(qmModels[i]);
			qmModels[i] = NULL;
		}
	}
	nG=numberOfGeometries;
	for (i = 0; i < nG; i++ ) 
	if(!qmModels[i]) 
	{
		for (j = i; j < nG-1; j++ ) 
		{
			qmModels[j]=qmModels[j+1];
			toRemove[j]=toRemove[j+1];
			if(energies) energies[j]=energies[j+1];
		}
		i--;
		nG--;
	}
	free(toRemove);
	fprintf(stderr,"After cutting by inertia, the new number of geometries = %d\n",nG);
	fflush(stderr);
	fprintf(logfile,"After cutting by inertia, the new number of geometries = %d\n",nG);
	fflush(logfile);
	*pnumberOfGeometries = nG;
}
/***********************************************************************************************************************/
static void removeSmallDistanceMolecules(QuantumMechanicsModel** qmModels, int* pnumberOfGeometries, double* energies, FILE* logfile)
{
	int i,j;
	int numberOfGeometries = *pnumberOfGeometries;
	int* toRemove = NULL;
	Molecule* mol;

	if(numberOfGeometries<2) return;
	toRemove = malloc(numberOfGeometries*sizeof(int));
	for (i = 0; i < numberOfGeometries; i++ ) toRemove[i]=0;

	//if(qmModels && qmModels[0] && qmModels[0]->molecule.nAtoms<1) return;
	for (i = 0; i < numberOfGeometries; i++ )
	{
		if(!qmModels[i]) { toRemove[i]=1; continue;}
		if(toRemove[i]==1) continue;

		mol = &qmModels[i]->molecule;
		if(mol->klass->smallDistance(mol)) 
		{
			if(energies) fprintf(logfile,"Mol num %4d : small distance, E = %0.8f\n",i+1,energies[i]);
			else fprintf(logfile,"Mol num %4d : small distance fragment\n",i+1);
			fflush(logfile);
			toRemove[i]=1;
		}
	}
	for (i = 0; i < numberOfGeometries; i++ )
	{
		if(toRemove[i]==1 && qmModels[i]) 
		{
			qmModels[i]->klass->free(qmModels[i]);
			qmModels[i] = NULL;
		}
	}
	int nG0=0; 
	for (i = 0; i < numberOfGeometries; i++ ) if(!qmModels[i])  nG0++;

	int nG=numberOfGeometries;
	for (i = 0; i < nG; i++ ) 
	if(!qmModels[i]) 
	{
		for (j = i; j < nG-1; j++ ) 
		{
			qmModels[j]=qmModels[j+1];
			toRemove[j]=toRemove[j+1];
			if(energies) energies[j]=energies[j+1];
		}
		i--;
		nG--;
	}
	free(toRemove);
	fprintf(stderr,"Number of removed geometries with small distance %d ; ",numberOfGeometries-nG);
	fprintf(stderr,"new number of geometries %d\n",nG);
	fflush(stderr);
	fprintf(logfile,"Number of removed geometries with small distance %d ; ",numberOfGeometries-nG);
	fprintf(logfile,"new number of geometries %d\n",nG);
	fflush(logfile);
	*pnumberOfGeometries = nG;
}
/***********************************************************************************************************************/
static void removeSimilarBondsGeometries(QuantumMechanicsModel** qmModels, int* pnumberOfGeometries, double* energies, FILE* logfile, double sTol, double distTol)
{
	int i,j;
	int numberOfGeometries = *pnumberOfGeometries;
	int* toRemove = NULL;
	Molecule* mol1;
	Molecule* mol2;

	if(numberOfGeometries<2) return;
	toRemove = malloc(numberOfGeometries*sizeof(int));
	for (i = 0; i < numberOfGeometries; i++ ) toRemove[i]=0;

	//if(qmModels && qmModels[0] && qmModels[0]->molecule.nAtoms<1) return;
	for (i = 0; i < numberOfGeometries; i++ )
	{
		if(!qmModels[i]) { toRemove[i]=1; continue;}
		if(toRemove[i]==1) continue;

		mol1 = &qmModels[i]->molecule;
		for (j = i+1; j < numberOfGeometries; j++ )
		{
			mol2 = &qmModels[j]->molecule;
			if(mol1->klass->similarBonds(mol1, mol2,sTol,distTol)) 
			{
				double s,maxDiffDistance;
				double I22,I42;
				s= mol1->klass->getSimilatityByBonds(mol1, mol2,&maxDiffDistance);
				mol2->klass->computePseudoInertia(mol2, &I22, &I42);
				fprintf(logfile,"Two (%4d, %4d) similar bonds : s=%+0.4f \tmax difference between distances=%+0.4f\n",
						i+1,j+1,s,maxDiffDistance);
				fflush(logfile);
				toRemove[j]=1;
			}
		}
	}
	for (i = 0; i < numberOfGeometries; i++ )
	{
		if(toRemove[i]==1 && qmModels[i]) 
		{
			qmModels[i]->klass->free(qmModels[i]);
			qmModels[i] = NULL;
		}
	}
	int nG0=0; 
	for (i = 0; i < numberOfGeometries; i++ ) if(!qmModels[i])  nG0++;

	int nG=numberOfGeometries;
	for (i = 0; i < nG; i++ ) 
	if(!qmModels[i]) 
	{
		for (j = i; j < nG-1; j++ ) 
		{
			qmModels[j]=qmModels[j+1];
			toRemove[j]=toRemove[j+1];
			if(energies) energies[j]=energies[j+1];
		}
		i--;
		nG--;
	}
	free(toRemove);
	fprintf(stderr,"Number of removed geometries due to similar bonds %d ; ",numberOfGeometries-nG);
	fprintf(stderr,"new number of geometries %d\n",nG);
	fflush(stderr);
	fprintf(logfile,"Number of removed geometries due to similar bonds %d ; ",numberOfGeometries-nG);
	fprintf(logfile,"new number of geometries %d\n",nG);
	fflush(logfile);
	*pnumberOfGeometries = nG;
}
/***********************************************************************************************************************/
static void removeFragmentedMolecules(QuantumMechanicsModel** qmModels, int* pnumberOfGeometries, double* energies, FILE* logfile)
{
	int i,j;
	int numberOfGeometries = *pnumberOfGeometries;
	int* toRemove = NULL;
	Molecule* mol;

	if(numberOfGeometries<2) return;
	toRemove = malloc(numberOfGeometries*sizeof(int));
	for (i = 0; i < numberOfGeometries; i++ ) toRemove[i]=0;

	//if(qmModels && qmModels[0] && qmModels[0]->molecule.nAtoms<1) return;
	for (i = 0; i < numberOfGeometries; i++ )
	{
		if(!qmModels[i]) { toRemove[i]=1; continue;}
		if(toRemove[i]==1) continue;

		mol = &qmModels[i]->molecule;
		if(!mol->klass->oneFragment(mol)) 
		{
			if(energies) fprintf(logfile,"Mol num %4d : more one fragment, E = %0.8f\n",i+1,energies[i]);
			else fprintf(logfile,"Mol num %4d : more one fragment\n",i+1);
			fflush(logfile);
			toRemove[i]=1;
		}
	}
	for (i = 0; i < numberOfGeometries; i++ )
	{
		if(toRemove[i]==1 && qmModels[i]) 
		{
			qmModels[i]->klass->free(qmModels[i]);
			qmModels[i] = NULL;
		}
	}
	int nG0=0; 
	for (i = 0; i < numberOfGeometries; i++ ) if(!qmModels[i])  nG0++;

	int nG=numberOfGeometries;
	for (i = 0; i < nG; i++ ) 
	if(!qmModels[i]) 
	{
		for (j = i; j < nG-1; j++ ) 
		{
			qmModels[j]=qmModels[j+1];
			toRemove[j]=toRemove[j+1];
			if(energies) energies[j]=energies[j+1];
		}
		i--;
		nG--;
	}
	free(toRemove);
	fprintf(stderr,"Number of removed fragmented geometries %d ; ",numberOfGeometries-nG);
	fprintf(stderr,"new number of geometries %d\n",nG);
	fflush(stderr);
	fprintf(logfile,"Number of removed fragmented geometries %d ; ",numberOfGeometries-nG);
	fprintf(logfile,"new number of geometries %d\n",nG);
	fflush(logfile);
	*pnumberOfGeometries = nG;
}
/***********************************************************************************************************************/
static void removeSimilarInertiaGeometries(QuantumMechanicsModel** qmModels, int* pnumberOfGeometries, double* energies, FILE* logfile, double tol)
{
	int i,j;
	int numberOfGeometries = *pnumberOfGeometries;
	int* toRemove = NULL;
	Molecule* mol1;
	Molecule* mol2;

	if(numberOfGeometries<2) return;
	toRemove = malloc(numberOfGeometries*sizeof(int));
	for (i = 0; i < numberOfGeometries; i++ ) toRemove[i]=0;

	//if(qmModels && qmModels[0] && qmModels[0]->molecule.nAtoms<1) return;
	for (i = 0; i < numberOfGeometries; i++ )
	{
		if(!qmModels[i]) { toRemove[i]=1; continue;}
		if(toRemove[i]==1) continue;

		mol1 = &qmModels[i]->molecule;
		for (j = i+1; j < numberOfGeometries; j++ )
		{
			mol2 = &qmModels[j]->molecule;
			if(mol1->klass->similarInertia(mol1, mol2,tol)) 
			{
				double I21,I41;
				double I22,I42;
				mol1->klass->computePseudoInertia(mol1, &I21, &I41);
				mol2->klass->computePseudoInertia(mol2, &I22, &I42);
				fprintf(logfile,"Two (%4d, %4d) similar geometries : I2=%+0.4f %+0.4f %+0.4f \tI4=%+0.4f %+0.4f %+0.4f\n",
						i+1,j+1,I21,I22,I22-I21,I41,I42,I42-I41);
				fflush(logfile);
				toRemove[j]=1;
			}
		}
	}
	for (i = 0; i < numberOfGeometries; i++ )
	{
		if(toRemove[i]==1 && qmModels[i]) 
		{
			qmModels[i]->klass->free(qmModels[i]);
			qmModels[i] = NULL;
		}
	}
	int nG0=0; 
	for (i = 0; i < numberOfGeometries; i++ ) if(!qmModels[i])  nG0++;

	int nG=numberOfGeometries;
	for (i = 0; i < nG; i++ ) 
	if(!qmModels[i]) 
	{
		for (j = i; j < nG-1; j++ ) 
		{
			qmModels[j]=qmModels[j+1];
			toRemove[j]=toRemove[j+1];
			if(energies) energies[j]=energies[j+1];
		}
		i--;
		nG--;
	}
	free(toRemove);
	fprintf(stderr,"Number of removed geometries due to similar inertia %d ; ",numberOfGeometries-nG);
	fprintf(stderr,"new number of geometries %d\n",nG);
	fflush(stderr);
	fprintf(logfile,"Number of removed geometries due to similar inertia %d ; ",numberOfGeometries-nG);
	fprintf(logfile,"new number of geometries %d\n",nG);
	fflush(logfile);
	*pnumberOfGeometries = nG;
}

/**********************************************************************/
static void setRattleConstraintsParameters(QuantumMechanicsModel* quantumMechanicsModel)
{
	Molecule* m = &quantumMechanicsModel->molecule;
	m->klass->resetConstraints(m, m->constraints);
}
/**********************************************************************/
QuantumMechanicsModel newQuantumMechanicsModel(char* method, char* dirName, char* nameCommand, char* N2P2Dir, char* tmModule, Constraints constraints, FILE* logfile)
{
	QuantumMechanicsModel qmModel;

	qmModel.molecule = *(newMolecule());

	qmModel.molecule.constraints = constraints;

	qmModel.klass = malloc(sizeof(QuantumMechanicsModelClass));
	qmModel.klass->calculateHessian = NULL;
	qmModel.klass->calculateGradient = NULL;
	qmModel.klass->calculateEnergy = NULL;
	qmModel.klass->free = freeQuantumMechanicsModel;
	qmModel.klass->computeQMFrequenciesNumeric = computeQMFrequenciesNumeric;
	qmModel.klass->computeQMFrequenciesAnalytic = computeQMFrequenciesAnalytic;
	qmModel.klass->computeQMFrequencies = computeQMFrequencies;
	qmModel.klass-> computeIR = computeIR;
	qmModel.klass->copy = copyQuantumMechanicsModel;
	qmModel.klass->setRattleConstraintsParameters = setRattleConstraintsParameters;
	qmModel.klass->removeFragmentedMolecules = removeFragmentedMolecules;
	qmModel.klass->removeSmallDistanceMolecules = removeSmallDistanceMolecules;
	qmModel.klass->removeSimilarInertiaGeometries = removeSimilarInertiaGeometries;
	qmModel.klass-> removeSimilarBondsGeometries = removeSimilarBondsGeometries;
	qmModel.klass-> sortByInertia = sortByInertia;
	qmModel.klass-> cutByInertia = cutByInertia;
	qmModel.klass-> getQuantumMechanicsRDConfo = getQuantumMechanicsRDConfo;
	qmModel.klass-> getQuantumMechanicsRDFConfo = getQuantumMechanicsRDFConfo;

	qmModel.firstRun = TRUE;
	qmModel.addD3Correction = FALSE;
	qmModel.addWallCorrection = FALSE;
	qmModel.logfile = logfile;
	qmModel.H4Parameters = NULL;
	qmModel.SRBParameters = NULL;
	qmModel.dx=1e-3;

	qmModel.N2P2Dir=NULL;
	if(N2P2Dir) qmModel.N2P2Dir=strdup(N2P2Dir);
	qmModel.interfaceLibN2P2 = NULL;
	qmModel.interfaceLibN2P2ES = NULL;

#ifdef ENABLE_PYTHON
	qmModel.tmModule=NULL;
	if(tmModule) qmModel.tmModule=strdup(tmModule);
	qmModel.interfaceTM = NULL;
#endif

	qmModel.method = NULL;
	if(method) qmModel.method = strdup(method);
	if(dirName) qmModel.workDir = strdup(dirName);
	else qmModel.workDir = strdup_printf("%s%stmp",cchemiDirectory(),DIR_SEPARATOR_S);
	if(nameCommand) qmModel.nameCommand = strdup(nameCommand);
	else qmModel.nameCommand = strdup("/opt/mopac/MOPAC2012");
	return qmModel;

}
/**********************************************************************/
static void freeQuantumMechanicsModel(QuantumMechanicsModel* qmModel)
{

	qmModel->molecule.klass->free(&qmModel->molecule);

	if(qmModel->klass != NULL)
	{
		free(qmModel->klass);
		qmModel->klass = NULL;
	}
	if(qmModel->method != NULL)
	{
		free(qmModel->method);
		qmModel->method = NULL;
	}
	if(qmModel->workDir != NULL)
	{
		free(qmModel->workDir);
		qmModel->workDir = NULL;
	}
	if(qmModel->nameCommand != NULL)
	{
		free(qmModel->nameCommand);
		qmModel->nameCommand = NULL;
	}
	if(qmModel->N2P2Dir != NULL)
	{
		free(qmModel->N2P2Dir);
		qmModel->N2P2Dir = NULL;
	}
	/* NON
	if(qmModel->interfaceLibN2P2 != NULL)
	{
		interfaceCChemIDestroy(qmModel->interfaceLibN2P2);
		qmModel->interfaceLibN2P2 = NULL;
	}
	if(qmModel->interfaceLibN2P2ES != NULL)
	{
		interfaceCChemIESDestroy(qmModel->interfaceLibN2P2ES);
		qmModel->interfaceLibN2P2ES = NULL;
	}
        */
#ifdef ENABLE_PYTHON
	if(qmModel->tmModule != NULL)
	{
		free(qmModel->tmModule);
		qmModel->tmModule = NULL;
	}
	if(qmModel->interfaceTM != NULL)
	{
		interfaceTMDestroy(qmModel->interfaceTM);
		qmModel->interfaceTM = NULL;
	}
#endif
}
/*****************************************************************************/
static QuantumMechanicsModel copyQuantumMechanicsModel(QuantumMechanicsModel* f)
{
	QuantumMechanicsModel qmModel = newQuantumMechanicsModel(NULL,NULL,NULL,NULL, NULL, NOCONSTRAINTS, stdout);

	qmModel.molecule = *(f->molecule.klass->copy(&f->molecule));
	qmModel.method = NULL;
	if(f->method) qmModel.method = strdup(f->method);
	qmModel.workDir = NULL;
	if(f->workDir) qmModel.workDir = strdup(f->workDir);

	qmModel.N2P2Dir = NULL;
	if(f->N2P2Dir) qmModel.N2P2Dir = strdup(f->N2P2Dir);
	qmModel.interfaceLibN2P2 = NULL;
	qmModel.interfaceLibN2P2ES = NULL;
	if(f->interfaceLibN2P2) qmModel.interfaceLibN2P2 = f->interfaceLibN2P2;// to do  make a copy of interfaceLibN2P2
	if(f->interfaceLibN2P2ES) qmModel.interfaceLibN2P2ES = f->interfaceLibN2P2ES;// to do make a copy of interfaceLibN2P2ES

#ifdef ENABLE_PYTHON
	qmModel.tmModule = NULL;
	if(f->tmModule) qmModel.tmModule = strdup(f->tmModule);
	qmModel.interfaceTM = NULL;
	if(f->interfaceTM) qmModel.interfaceTM = f->interfaceTM;// make a copy of interfaceTM
#endif

	qmModel.nameCommand = NULL;
	if(f->nameCommand) qmModel.nameCommand = strdup(f->nameCommand);

	qmModel.klass->calculateHessian = f->klass->calculateHessian;
	qmModel.klass->calculateGradient = f->klass->calculateGradient;
	qmModel.klass->calculateEnergy = f->klass->calculateEnergy;
	qmModel.logfile = f->logfile;
	qmModel.firstRun = f->firstRun;
	qmModel.addD3Correction = f->addD3Correction;
	qmModel.addWallCorrection = f->addWallCorrection;
	qmModel.H4Parameters = NULL;
	qmModel.SRBParameters = NULL;
	if(f->SRBParameters) *(qmModel.SRBParameters) = *(f->SRBParameters);

	return qmModel;
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
static void sortFrequencies(int nModes, double* frequencies, double** modes, double* reducedMasses, double* IRIntensities)
{
	int i;
	int j;
	int k;
	double dum;
	if(nModes<1 || !frequencies || !modes || !reducedMasses || !IRIntensities) return;
	for(i=0;i<nModes;i++)
	{
		k = i;
		for(j=i+1;j<nModes;j++)
			if(frequencies[j]<frequencies[k]) k = j;
		if(k==i) continue;
		/* swap i and k modes */
		dum = frequencies[i];
		frequencies[i] = frequencies[k];
		frequencies[k] = dum;
		dum = reducedMasses[i];
		reducedMasses[i] = reducedMasses[k];
		reducedMasses[k] = dum;
		dum = IRIntensities[i];
		IRIntensities[i] = IRIntensities[k];
		IRIntensities[k] = dum;
		for(j=0;j<nModes;j++)
		{
			dum =  modes[j][i];
			modes[j][i] = modes[j][k];
			modes[j][k] = dum;
		}
	}
}
/****************************************************************************************************************************************************/
static int computeIR(QuantumMechanicsModel* qmModel, double *F, double** dmu, double** frequencies, double*** modes, double** reducedMasses, double** IRIntensities)
{
	Molecule* mol = &qmModel->molecule;
	int nAtoms = mol->nAtoms;
	int i;
	int j;
	int k;
	int c;
	*frequencies = malloc(3*nAtoms*sizeof(double));
	*reducedMasses = malloc(3*nAtoms*sizeof(double));
	*IRIntensities = malloc(3*nAtoms*sizeof(double));
	*modes = malloc(3*nAtoms*sizeof(double*));
	for(i=0;i<3*nAtoms;i++) (*modes)[i] = malloc(3*nAtoms*sizeof(double));

	//printf("begin diag\n");
	eigenQL(3*nAtoms, F, *frequencies, *modes);

	//printf("end eigneQL\n");
	for(i=0;i<3*nAtoms;i++) (*IRIntensities)[i] = 0.0;
	
	/* convert in atomic unit  from kcal/Ang^2/amu */
	for(i=0;i<3*nAtoms;i++) (*frequencies)[i] *= 1.59360150e-03*0.529177*0.529177*5.48579911e-04; 
	/* convert frequencies in cm-1 */
	for(i=0;i<3*nAtoms;i++) 
		if( (*frequencies)[i]>0) (*frequencies)[i] = sqrt((*frequencies)[i])*219474.63633664;
		else (*frequencies)[i] = -sqrt(-(*frequencies)[i])*219474.63633664;

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
			double Lji = (*modes)[jd][id];
			double a=dmu[c][jd]*Lji/sqrt(mol->atoms[j].mass);
			D[c]+=a;
		}
		IRI = 0;
		for(c = 0;c<3;c++)  IRI+= D[c]*D[c];
		(*IRIntensities)[id] = IRI;
	}
	/* Intensities in 1 (D/Ang)^2 amu^-1 = 42.255 km/mol=171.65 cm^-2 atm^-1 at 0 C and 1 atm */
	/* Refs : D. Porezag and M. R. Pederson, Phys. Rev. B 54, 7830 (1996). and Y. Yamaguchi el al., J. Chem. Phys. 84,2262(1986)*/
	/* conversion in km/mol*/
	for(i=0;i<3*nAtoms;i++) (*IRIntensities)[i] *= 42.255;

	/* compute the reduced mass */
	for(i=0;i<3*nAtoms;i++) 
	{
		double m = 0;
		for(j=0;j<mol->nAtoms;j++)
		{
			double r2 = 0;
			for(c=0;c<3;c++) r2+= (*modes)[3*j+c][i]*(*modes)[3*j+c][i];
			m+= r2/(mol->atoms[j].mass); 
		}
		if(m<=0) m = 1;
		m = 1/m;
		for(j=0;j<mol->nAtoms;j++)
		{
			double r =sqrt(m)/sqrt(mol->atoms[j].mass);
			for(c=0;c<3;c++) (*modes)[3*j+c][i]*=r;
		}

		//printf("%f %f\n",(*frequencies)[i],m);
		(*reducedMasses)[i] = m;
	}
	sortFrequencies(3*nAtoms, *frequencies, *modes, *reducedMasses, *IRIntensities);
	return 3*nAtoms;
}
/*****************************************************************************/
static int computeQMFrequenciesNumeric(QuantumMechanicsModel* qmModel, double** frequencies, double*** modes, double** reducedMasses, double** IRIntensities)
{
	int i;
	int j;
	int k;
	int c;
	int id,jd,index;
	double* F;
	double* gp[3];
	double* gm[3];
	double** dmu;
	Molecule* mol;
	double Dp[3];
	double Dm[3];
	int nAtoms;
	boolean show = FALSE;
	double dx;
	int ret;

	if(!qmModel || qmModel->molecule.nAtoms<1) return 0;
	dx = qmModel->dx;

	mol = &qmModel->molecule;
	nAtoms = mol->nAtoms;
	if(nAtoms>4) show = TRUE;
	printf("nAtoms = %d\n",nAtoms);
	for(k=0;k<3;k++) gp[k] = malloc(nAtoms*sizeof(double));
	for(k=0;k<3;k++) gm[k] = malloc(nAtoms*sizeof(double));
	dmu = malloc(3*sizeof(double*));
	for(k=0;k<3;k++) dmu[k] = malloc(3*nAtoms*sizeof(double));

	F = malloc(3*nAtoms*(3*nAtoms+1)/2*sizeof(double));

	index = 0;
	for(i=0;i<nAtoms;i++)
	for(k=0;k<3;k++)
	{
		id=3*i+k;
		if(show && k ==0) { printf("Computing of derivatives for atom # %d/%d\n",i+1,nAtoms);};
		mol->atoms[i].coordinates[k] += dx;
		qmModel->klass->calculateGradient(qmModel);
		copyGradients(mol, gp);
		for(c = 0;c<3;c++)  Dp[c] = mol->dipole[c];
		
		mol->atoms[i].coordinates[k] -= 2*dx;
		qmModel->klass->calculateGradient(qmModel);
		copyGradients(mol, gm);
		for(c = 0;c<3;c++)  Dm[c] = mol->dipole[c];
		for(c = 0;c<3;c++) dmu[c][id] = (Dp[c]-Dm[c])/dx/2;
		mol->atoms[i].coordinates[k] += dx;

		for(j=0;j<=i;j++)
		{
			double invm = 1.0/sqrt( mol->atoms[i].mass* mol->atoms[j].mass);
			for(c = 0;c<3;c++) 
			{
				jd = 3*j+c;
				//printf("id = %d jd = %d\n",id,jd);
				if(jd>id) continue;
				index = jd + id*(id+1)/2;
				//printf("index = %d i = %d k = %d j = %d c = %d\n",index,i,k,j,c);
				F[index] = (gp[c][j]-gm[c][j])/dx/2; 
				F[index] *= invm;
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
	ret = computeIR(qmModel, F, dmu, frequencies, modes, reducedMasses, IRIntensities);
	free(F);
	for(k=0;k<3;k++) free(dmu[k]);
	if(dmu) free(dmu);
	return ret;

}
/*****************************************************************************/
static int computeQMFrequenciesAnalytic(QuantumMechanicsModel* qmModel, double** frequencies, double*** modes, double** reducedMasses, double** IRIntensities)
{
	double* F;
	double** dmu;
	int ret;
	int k;

	if(!qmModel || qmModel->molecule.nAtoms<1) return 0;
	if(!qmModel->interfaceLibN2P2 || !qmModel->interfaceLibN2P2ES) return 0;
	qmModel->klass->calculateHessian(qmModel,&F, &dmu);
	ret = computeIR(qmModel, F, dmu, frequencies, modes, reducedMasses, IRIntensities);
	free(F);
	for(k=0;k<3;k++) free(dmu[k]);
	if(dmu) free(dmu);
	return ret;

}
static int computeQMFrequencies(QuantumMechanicsModel* qmModel, double** frequencies, double*** modes, double** reducedMasses, double** IRIntensities)
{
	int ret = 0;
	if(qmModel->SRBParameters || qmModel->addD3Correction || qmModel->addWallCorrection || qmModel->H4Parameters || !qmModel->interfaceLibN2P2 || !qmModel->interfaceLibN2P2ES)
		ret = computeQMFrequenciesNumeric(qmModel, frequencies, modes, reducedMasses, IRIntensities);
	else 
		ret = computeQMFrequenciesAnalytic(qmModel, frequencies, modes, reducedMasses, IRIntensities);

	return ret;

}
