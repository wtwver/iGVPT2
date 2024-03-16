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

#ifndef __CCHEMILIB_QUANTUMMECHANICSMODEL_H__
#define __CCHEMILIB_QUANTUMMECHANICSMODEL_H__

#include "../Molecule/Molecule.h"
#include "../EmpriricalCorrections/HydrogenBondCorrection.h"
#include "../EmpriricalCorrections/ShortRangeBasisSetCorrection.h"
#include "../InterfaceLibN2P2/InterfaceCChemIC.h"
#include "../InterfaceTensorMol/InterfaceTM.h"


typedef struct _QuantumMechanicsModel  QuantumMechanicsModel;
typedef struct _QuantumMechanicsModelClass  QuantumMechanicsModelClass;
typedef struct _QuantumMechanicsModelOptions  QuantumMechanicsModelOptions;
#define RATTLEDIM	        3 /* a1 a2 r12 */

struct _QuantumMechanicsModel
{
	Molecule molecule;
	QuantumMechanicsModelClass* klass;
	char* method;
	char* workDir;
	char* nameCommand;
	char* N2P2Dir;
	FILE* logfile;
	boolean firstRun;
	boolean addD3Correction;
	boolean addWallCorrection;
	HyhrogenBondCorrectionParameters* H4Parameters;
	ShortRangeBasisSetCorrectionParameters* SRBParameters;
	void* interfaceLibN2P2;
	void* interfaceLibN2P2ES;
#ifdef ENABLE_PYTHON
	InterfaceTM* interfaceTM;
#endif
	double dx;
	char* tmModule;

};
struct _QuantumMechanicsModelClass
{
	void (*calculateHessian)(QuantumMechanicsModel* qmModel, double **F, double*** dmu);
	void (*calculateGradient)(QuantumMechanicsModel* qmModel);
	void (*calculateEnergy)(QuantumMechanicsModel* qmModel);
	void (*free)(QuantumMechanicsModel* qmModel);
	int (*computeQMFrequenciesNumeric)(QuantumMechanicsModel* qmModel, double** frequencies, double*** modes, double** reducedMasses, double** IRIntensities);
	int (*computeQMFrequenciesAnalytic)(QuantumMechanicsModel* qmModel, double** frequencies, double*** modes, double** reducedMasses, double** IRIntensities);
	int (*computeQMFrequencies)(QuantumMechanicsModel* qmModel, double** frequencies, double*** modes, double** reducedMasses, double** IRIntensities);
	QuantumMechanicsModel (*copy)(QuantumMechanicsModel* qmModel);
	void (*setRattleConstraintsParameters)(QuantumMechanicsModel* qmModel);
	void (*removeFragmentedMolecules)(QuantumMechanicsModel** qmModels, int* pnumberOfGeometries, double* energies, FILE* logfile);
	void (*removeSmallDistanceMolecules)(QuantumMechanicsModel** qmModels, int* pnumberOfGeometries, double* energies, FILE* logfile);
	void (*removeSimilarInertiaGeometries)(QuantumMechanicsModel** qmModels, int* pnumberOfGeometries, double* energies, FILE* logfile, double tol);
	 void (*removeSimilarBondsGeometries)(QuantumMechanicsModel** qmModels, int* pnumberOfGeometries, double* energies, FILE* logfile, double sTol, double distTol);
	void (*sortByInertia)(QuantumMechanicsModel** qmModels, int* pnumberOfGeometries, double* energies);
	void (*cutByInertia)(QuantumMechanicsModel** qmModels, int* pnumberOfGeometries, double* energies, int newNGeoms, FILE* logfile);
	QuantumMechanicsModel** (*getQuantumMechanicsRDConfo)(QuantumMechanicsModel* qmModel, int numberOfGeometries, boolean chain, boolean saveFirstGeom);
	QuantumMechanicsModel** (*getQuantumMechanicsRDFConfo)(QuantumMechanicsModel* qmModel, int numberOfGeometries, boolean saveFirstGeom);
	int (*computeIR)(QuantumMechanicsModel* qmModel, double *F, double** dmu, double** frequencies, double*** modes, double** reducedMasses, double** IRIntensities);
};


QuantumMechanicsModel newQuantumMechanicsModel(char* method, char* dirName, char* nameCommand, char* N2P2Dir,  char* tmModule, Constraints constraints, FILE* logfile);

#endif /* __CCHEMILIB_QUANTUMMECHANICSMODEL_H__ */

