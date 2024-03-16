/********************************************************************************
 cchemi is an interface to ab initio computational chemistry programs 
 designed for add them many functionalities non available in these packages.
 Copyright (C) 2020 Abdulrahman Allouche (University Lyon 1)

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

#ifndef __CCHEMILIB_MOLECULE_H__
#define __CCHEMILIB_MOLECULE_H__

#include <stdio.h>
#include "../Utils/Types.h"
#include "../Molecule/Atom.h"

#define RATTLEDIM	        3 /* a1 a2 r12 */

typedef struct _VibMode VibMode;
typedef struct _Vibration Vibration;
typedef struct _Molecule  Molecule;
typedef struct _MoleculeClass  MoleculeClass;

struct _VibMode
{
	double frequency;
	double mass;
	double* properties;
	double* vectors[3];
};
struct _Vibration
{
	VibMode* modes;
	int nModes;
	int nProperties;
};

struct _MoleculeClass
{
	boolean (*saveGeometry)(Molecule* molecule, char* fileName);
	boolean (*saveGeometryAndVelocities)(Molecule* molecule, char* fileName);
	void (*scaleVelocities)(Molecule* molecule, double temperature);
	void (*setMaxwellVelocities)(Molecule* molecule, double temperature);
	boolean (*setMaxwellVelocitiesIfNull)(Molecule* molecule, double temperature);
	boolean (*print)(Molecule* molecule, FILE* file);
	double (*getGradientNorm)(Molecule* molecule);
	double (*getKelvin)(Molecule* molecule);
	double (*getKineticEnergy)(Molecule* molecule);
	void (*computeDipole)(Molecule* molecule);
	double (*getKineticEnergyMass)(Molecule* molecule, double* masses);
	void (*removeTranslation)(Molecule* molecule);
	void (*removeRotation)(Molecule* molecule);
	void (*removeTranslationAndRotation)(Molecule* molecule);

	void (*removeTranslationCluster)(Molecule** molecules, int nMols);
	void (*removeRotationCluster)(Molecule** molecules, int nMols);
	void (*removeTranslationAndRotationCluster)(Molecule** molecules, int nMols);

	void (*removeTranslationForceCluster)(Molecule** molecules, int nMols, double** vectors);
	void (*removeRotationForceCluster)(Molecule** molecules, int nMols, double** vectors);
	void (*removeTranslationAndRotationForceCluster)(Molecule** molecules, int nMols, double** vectors);

	void (*removeTranslationForce)(Molecule* molecule, double* f);
	void (*removeRotationForce)(Molecule* molecule, double* f);
	void (*removeTranslationAndRotationForce)(Molecule* molecule, double* f);

	void (*removeTranslationMoments)(Molecule** molecules, int nMols, double*** P);
	void (*removeRotationMoments)(Molecule** molecules, int nMols, double*** P);
	void (*removeTranslationAndRotationMoments)(Molecule** molecules, int nMols, double*** P);

	void (*removeTranslationAcceleration)(Molecule* molecule,double* a);
	void (*removeRotationAcceleration)(Molecule* molecule, double* a);
	void (*removeTranslationAndRotationAcceleration)(Molecule* molecule, double* a);
	boolean (*resetConstraints)(Molecule* molecule, Constraints constraints);
	void (*setRattleConstraintsParameters)(Molecule* molecule);
	boolean (*setRandomPositions)(Molecule* molecule);
	boolean (*setRandomFragments)(Molecule* molecule);
	boolean (*setRandomPositionsChain)(Molecule* molecule);
	boolean (*addGeometry)(Molecule* molecule,FILE* file);
	boolean (*addMolecule)(Molecule* molecule,FILE* file);
	boolean (*addVelocities)(Molecule* molecule,FILE* file);
	boolean (*addFirstDerivativeToFile)(Molecule* molecule, FILE* file);
	boolean (*readGeometry)(Molecule* molecule,char* fileName);
	void (*copyChargeInCharge0)(Molecule* molecule);
	void (*setBondHardness)(Molecule* molecule, int nBonds, char** atomTypes1, char** atomTypes2, double* hardness);
	void (*setHardness)(Molecule* molecule, int nTypes, char** atomTypes, double* hardness);
	void (*setElectronegativity)(Molecule* molecule, int nTypes, char** atomTypes, double* electronegativity);
	void (*setWidth)(Molecule* molecule, int nTypes, char** atomTypes, double* width);
	void (*setCharge0)(Molecule* molecule, int nTypes, char** atomTypes, double* charge0);
	void (*setChargesACKS2)(Molecule* molecule);
	void (*setChargesEEM)(Molecule* molecule);
	double (*getEnergyEEM)(Molecule* molecule);
	double (*getEnergyACKS2)(Molecule* molecule);
	boolean (*isLinear)(Molecule* molecule);
	boolean (*save)(Molecule* molecule, char* fileName);
	boolean (*addGeometryToGabedit)(Molecule* molecule,FILE* file);
	boolean (*saveFrequencies)(Molecule* molecule, char* fileName, int nModes, double* frequencies, double** modes, double* reducedMasses, double* IRIntensities);
	int (*computeIR)(Molecule* mol, double *F, double* dmuX[3]);
	void (*removeFrequencies)(Molecule* molecule, double freqMin, double freqMax);
	void (*removeNotSelectedFrequencies)(Molecule* molecule, boolean* selected);
	boolean (*readGradientFromGabeditFile)(Molecule* molecule, char* namefile);
	int (*computeFrequenciesFromFiles)(Molecule* mol, char* inputFileName, double dx);
	int (*computeFrequenciesOneStepFromFiles)(Molecule* mol, char* inputFileName, double dx);
	int (*computeFrequenciesFromGradFiles)(Molecule* mol, char* inputFileName, double dx);
	int (*computeFrequenciesOneStepFromGradFiles)(Molecule* mol, char* inputFileName, double dx);
	int (*generateCChemIFilesForFrequencies)(Molecule* mol, char* inputFileName, double dx);
	int (*generateCChemIFilesOneStepForFrequencies)(Molecule* mol, char* inputFileName, double dx);
	int (*generateCChemIGradFilesForFrequencies)(Molecule* mol, char* inputFileName, double dx);
	int (*generateCChemIGradFilesOneStepForFrequencies)(Molecule* mol, char* inputFileName, double dx);


	void (*free)(Molecule* molecule);
	Molecule* (*copy)(Molecule* molecule);

	void (*computeHarmonicVelocitiesCoordinates)(Molecule* mol, double T, int numMode, boolean changeGeom);
	int (*generateQFFCChemIFiles)(Molecule* mol, char* inputFileName, double delta, boolean reducedCoordinates, int ordre);
	double*** (*getGeomsQFF)(Molecule* mol, char* inputFileName, double delta, boolean reducedCoordinates, int ordre, double** pDeltas, int* nGeoms);
	boolean (*read)(Molecule* mol, char *fileName);
	void (*readGeomFromMopacOutputFile)(Molecule* mol, char *fileName, int numgeometry);
	void (*readGeomFromGamessOutputFile)(Molecule* mol, char *fileName, int numgeometry);
	void (*readGeomFromOrcaOutputFile)(Molecule* mol, char* fileName, int numgeometry);
	void (*readGeomFromOpenBabelOutputFile)(Molecule* mol, char* fileName, int numgeometry);
	void (*readGeomFromGaussianOutputFile)(Molecule* mol, char *fileName, int numgeometry);
	void (*readGeomFromMopacAuxFile)(Molecule* mol, char* fileName, int numgeometry);
	boolean (*saveMol2)(Molecule* mol, char* fileName);
	boolean (*saveHIN)(Molecule* mol, char* fileName);
	void (*removeTransRotModes)(Molecule* mol);
	void (*setConnections)(Molecule* mol);
	boolean (*resetMMTypes)(Molecule* mol, char* type);
	boolean (*buildMMTypes)(Molecule* mol, FILE* file);
        boolean (*fit)(Molecule* molToFit, Molecule* molRef, double u[3][3]);
        boolean (*fitRMSD)(Molecule* molToFit, Molecule* molRef, boolean centerRef);
        void (*center)(Molecule* mol, double C[]);
        void (*buildStandardOrientation)(Molecule* mol,double* centerOfGravity, int* numberOfEquivalentAxes, double* inertialMoment, double axes[3][3]);
	boolean (*compare)(Molecule* mol1, Molecule* mol2, boolean warning);
	char* (*saveFirstDerivatives)(char* inputFileName, Molecule* mol);

	double* (*getDeltaTable)(Molecule* mol, double delta, boolean reducedCoordinates);

	int (*getQFFOneGeom)(Molecule* mol, int mode1, int mode2, double delta1, double delta2, double akOverI1, double** geom);
	int (*getQFFOneGeom3)(Molecule* mol, int mode1, int mode2, int mode3, double delta1, double delta2, double delta3, double** geom);
	int (*getQFFOneGeom4)(Molecule* mol, int mode1, int mode2, int mode3, int mode4, double delta1, double delta2, double delta3, double delta4, double** geom);
	void (*computePseudoInertia)(Molecule* mol, double*pI, double* pI4);
	boolean (*similarInertia)(Molecule* mol1, Molecule* mol2, double precision);
	boolean (*similarBonds)(Molecule* mol1, Molecule* mol2, double sTol, double distMaxTol);
	boolean (*makeLocalSphericalMutation)(Molecule* molecule, double rate);
	boolean (*makeCenterOfMassSphericalMutation)(Molecule* molecule);
	boolean (*setGeometryToAxes)(Molecule* molecule, double axis1[], double axis2[], double axis3[]);
	boolean (*setRandomOrientation)(Molecule* molecule);
	Molecule* (*makeSphereCutSpliceCrossover)(Molecule* mol1, Molecule* mol2, int *err);
	Molecule* (*makePlaneCutSpliceCrossover)(Molecule* mol1, Molecule* mol2, int *err);
	boolean (*smallDistance)(Molecule* mol);
	boolean (*oneFragment)(Molecule* mol);
	double* (*getDistancesBetweenAtoms)(Molecule* mol);
	double (*getSimilatityByBonds)(Molecule* mol1, Molecule* mol2, double* pmaxDifference);

};

typedef struct 
{
	double lengths[3];
	double alpha;
	double beta;
	double gamma;
	double lVectors[3][3];
	double reciprocVectors[3][3];
}Boxes;

/* Spherical wall E0 *(1-exp(-r^2/rho^2)**n */
typedef struct
{
        double E0;
        double rho;
        int n;
}WallParameters;

struct _Molecule
{
	MoleculeClass* klass;
	int nAtoms;
	Atom* atoms;
	double potentialEnergy;
	double dipole[3];
	int numberOf2Connections;
	int* connected2[2];
	int numberOf3Connections;
	int* connected3[3];
	int numberOf4Connections;
	int* connected4[4];
	int numberOfNonBonded;
	int* nonBonded[2];

	double* bondHardness;// nAtoms*(nAtoms+1)/2 elements

	int spinMultiplicity;
	int totalCharge;
	Constraints constraints;
	int nFree;
	int numberOfRattleConstraintsTerms;
	double* rattleConstraintsTerms[RATTLEDIM];
	Boxes boxes;
	WallParameters wall;
	Vibration vibration;
	
};

Molecule* newMolecule();

Molecule* readMolecule(char* fileName, boolean connections);
Molecule* readMoleculeFromCChemIFile(char* fileName, boolean connections);
Molecule* readMoleculeFromMopacOutputFile(char *fileName, int numgeometry);
Molecule* readMoleculeFromMopacAuxFile(char *fileName, int numgeometry);
Molecule* readMoleculeFromGamessOutputFile(char *fileName, int numgeometry);
Molecule* readMoleculeFromOrcaOutputFile(char* fileName, int numgeometry);
Molecule* readMoleculeFromOrcaHessianFile(char* fileName);
Molecule* readMoleculeFromGaussianOutputFile(char *fileName, int numgeometry);
Molecule* readMoleculeFromGabeditFile(char* namefile);
Molecule* readMoleculeFromGaussianFChkFile(char *fileName);
Molecule* readMoleculeFromOrcaHessianFile(char *fileName);
Molecule** readMolecules(char* fileName, boolean connections);

void exchangeGeometries(Molecule* mol1, Molecule* mol2);

#endif /* __CCHEMILIB_MOLECULE_H__ */

