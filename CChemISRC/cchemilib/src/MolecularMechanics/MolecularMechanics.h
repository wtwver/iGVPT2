
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

#ifndef __CCHEMILIB_MOLECULARMECHANICS_H__
#define __CCHEMILIB_MOLECULARMECHANICS_H__

#include "../Molecule/Molecule.h"
#include "../MolecularMechanics/ForceField.h"

typedef struct _AmberParameters  AmberParameters;
typedef struct _AmberAtomTypes  AmberAtomTypes;
typedef struct _AmberBondStretchTerms  AmberBondStretchTerms;
typedef struct _AmberBondHardnessTerms  AmberBondHardnessTerms;
typedef struct _AmberAngleBendTerms  AmberAngleBendTerms;
typedef struct _AmberStrBendTerms  AmberStrBendTerms;
typedef struct _AmberDihedralAngleTerms  AmberDihedralAngleTerms;
typedef struct _AmberImproperTorsionTerms  AmberImproperTorsionTerms;
typedef struct _AmberOutOfPlaneTerms  AmberOutOfPlaneTerms;
typedef struct _AmberVdw612Terms  AmberVdw612Terms;
typedef struct _AmberVdw714Terms  AmberVdw714Terms;
typedef struct _AmberHydrogenBonded1012Terms  AmberHydrogenBonded1012Terms;
typedef struct _AmberHydrogenBondedMorseTerms  AmberHydrogenBondedMorseTerms;
typedef struct _AmberSuttonChenTerms  AmberSuttonChenTerms;
typedef struct _AmberPairWiseTerms AmberPairWiseTerms;

/************************************/
struct _AmberAtomTypes
{
	char* name;
	char* symbol;
	int number;
	double mass;
	double polarisability;
	double charge;
	double electronegativity;
	double hardness;
	double width;
};
/************************************/
struct _AmberBondStretchTerms
{
	int numbers[2];
	double equilibriumDistance;
	double forceConstant;
	double h3;
	double h4;
	double h5;
	double h6;
	int type;
};
/************************************/
struct _AmberAngleBendTerms
{
	int numbers[3];
	double equilibriumAngle;
	double forceConstant;
	double h3;
	double h4;
	double h5;
	double h6;
};
/************************************/
struct _AmberStrBendTerms
{
	int numbers[3];
	double forceConstant12;
	double forceConstant23;
	double Re12;
	double Re23;
	double Theta123;
};
/************************************/
struct _AmberDihedralAngleTerms
{
	int numbers[4];
	int nSomme;
	double* divisor;
	double* barrier;
	double* phase;
	double* n;
};
/************************************/
struct _AmberImproperTorsionTerms
{
	int numbers[4];
	double barrier;
	double phase;
	double n;
};
/************************************/
struct _AmberOutOfPlaneTerms
{
	int numbers[4];
	int type;
	double force;
	double h3;
	double h4;
	double h5;
	double h6;
};
/************************************/
struct _AmberVdw612Terms
{
	int number;
	double r;
	double epsilon;
};
/************************************/
struct _AmberVdw714Terms
{
	int number;
	double r;
	double epsilon;
	double gamma;
	double delta;
};
/************************************/
struct _AmberHydrogenBonded1012Terms
{
	int numbers[2];
	double c;
	double d;
};
/************************************/
struct _AmberHydrogenBondedMorseTerms
{
	int numbers[2];
	double force;
	double Re;
	double De; 
};
/************************************/
struct _AmberSuttonChenTerms
{
	int numbers[2];
	double epsilon;
	double a;
	double C;
	double n;
	double m;
};
/************************************/
struct _AmberPairWiseTerms
{
	int numbers[2];
	double a;
	double beta;
	double c4;
	double c6;
	double c8;
	double c10;
	double b;
};
/************************************/
struct _AmberBondHardnessTerms
{
	int numbers[2];
	double kappa;
};
/************************************/
struct _AmberParameters
{
	int numberOfTypes;
	AmberAtomTypes* atomTypes;

	int numberOfStretchTerms;
	AmberBondStretchTerms* bondStretchTerms;

	int numberOfBendTerms;
	AmberAngleBendTerms* angleBendTerms;

	int numberOfStrBendTerms;
	AmberStrBendTerms* strBendTerms;

	int numberOfDihedralTerms;
	AmberDihedralAngleTerms* dihedralAngleTerms;

	int numberOfImproperTorsionTerms;
	AmberImproperTorsionTerms* improperTorsionTerms;

	int numberOfOutOfPlaneTerms;
	AmberOutOfPlaneTerms* outOfPlaneTerms;

	int numberOfVdw612;
	AmberVdw612Terms* vdw612Terms;

	int numberOfVdw714;
	AmberVdw714Terms* vdw714Terms;

	int numberOfHydrogenBonded1012;
	AmberHydrogenBonded1012Terms* hydrogenBonded1012Terms;

	int numberOfHydrogenBondedMorse;
	AmberHydrogenBondedMorseTerms* hydrogenBondedMorseTerms;

	int numberOfSuttonChen;
	AmberSuttonChenTerms* suttonChenTerms;

	int numberOfPairWise;
	AmberPairWiseTerms* pairWiseTerms;

	int numberOfHardnessTerms;
	AmberBondHardnessTerms* bondHardnessTerms;
};
/************************************/
ForceField createAmberModel(Molecule* mol,ForceFieldOptions forceFieldOptions, FILE* logfile);
ForceField createPairWiseModel(Molecule* mol,ForceFieldOptions forceFieldOptions, FILE* logfile);
void loadAmberParameters();
void saveAmberParameters();
AmberParameters* getPointerAmberParameters();
void setPointerAmberParameters(AmberParameters* ptr);
AmberParameters newAmberParameters();
char** getListMMTypes(int* nlist);
void initCLForceField (ForceField* forceField);
void setH4CorrectionMM(FILE* file, ForceField* forceField);

#endif /* __CCHEMILIB_MOLECULARMECHANICS_H__ */

