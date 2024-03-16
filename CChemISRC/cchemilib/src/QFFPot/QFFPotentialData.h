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

/* QFFPotentialData.h */

#ifndef __QFF_PotentialData_H__
#define __QFF_PotentialData_H__
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "../Molecule/Molecule.h"

typedef struct _QFFPotentialData        QFFPotentialData;
typedef struct _QFFPotentialDataClass   QFFPotentialDataClass;
typedef struct _QFFPotnMR   QFFPotnMR;
typedef struct _QFFPotParameters   QFFPotParameters;

/************************************/
struct _QFFPotnMR
{
	int numbers[4];
	double energy;
	double grad;
};
/************************************/
struct _QFFPotParameters
{
	int numberOf1MR;
	QFFPotnMR* qff1MR;

	int numberOf2MR;
	QFFPotnMR* qff2MR;

	int numberOf3MR;
	QFFPotnMR* qff3MR;

	int numberOf4MR;
	QFFPotnMR* qff4MR;
};
struct _QFFPotentialData
{
	Molecule molecule;
	int nFrequencies;
	double ***modes;
	double *effectiveMasses;
	double *gradients;
	double **hessian;
	double ***cubic;
	double ****quartic;
	double ***coriolis;
	double *Be;
	double ***MWModes;
        QFFPotentialDataClass* klass;
	QFFPotParameters* qffPotParameters;
};

struct _QFFPotentialDataClass
{
        void (*readData)(QFFPotentialData* qffPotentialData, char* fileName);
        void (*convertToAU)(QFFPotentialData* qffPotentialData);
        void (*convertToAU2)(QFFPotentialData* qffPotentialData);
        void (*print)(QFFPotentialData* qffPotentialData);
        void (*free)(QFFPotentialData* qffPotentialData);
        void (*rotModes)(QFFPotentialData* qffPotentialData, double u[3][3]);
        void (*computeQFFParameters)(QFFPotentialData* qffPotentialData);
};

QFFPotentialData newQFFPotentialData(int n, int nAtoms);

#endif /* __QFF_PotentialData_H__ */
