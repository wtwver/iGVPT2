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

/* QFFModel.h */

#ifndef __QFF_Model_H__
#define __QFF_Model_H__
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "../QFFPot/QFFPotentialData.h"
#include "../QFFPot/QFFPropertiesData.h"
#include "../Molecule/Molecule.h"

typedef struct _QFFModel        QFFModel;
typedef struct _QFFModelClass   QFFModelClass;

struct _QFFModel
{
	Molecule molecule;
	QFFPotentialData vData;
	QFFPropertiesData pData;
        QFFModelClass* klass;
	double* Q;
	double* gradQ;
	double* velocity;
	boolean* variable;
	double diffStep;
        double* frequencies;
        double** modes;
        double* reducedMasses;
        double* IRIntensities;
	int typeCalcul; // 0 = energy, 1 = gradients, 2 = harmonic ferquencies 
	boolean showWarning; // 0 = not show, 1 = show 
	boolean numGradients;
	boolean xyz;

};

struct _QFFModelClass
{
        void (*readData)(QFFModel* qffModel, char* geometryFileName, char* qffFileName);
        void (*convertToAU)(QFFModel* qffModel);
        void (*convertToAU2)(QFFModel* qffModel);
        void (*computeEnergy)(QFFModel* qffModel);
        void (*computeGradients)(QFFModel* qffModel);
        int (*computeFrequencies)(QFFModel* qffModel);
        void (*compareFrequencies)(QFFModel* qffModel,FILE*);
        void (*computeDipole)(QFFModel* qffModel);
        void (*compute)(QFFModel* qffModel);
        void (*calculateGradient)(QFFModel* qffModel);
        void (*runAlea)(QFFModel* qffModel);
        void (*free)(QFFModel* qffModel);
	double (*getKineticEnergy)(QFFModel* qffModel);
	double (*getKelvin)(QFFModel* qffModel);
	void (*scaleVelocities)(QFFModel* qffModel, double temperature);
	void (*setMaxwellVelocities)(QFFModel* qffModel, double temperature);
	boolean (*setMaxwellVelocitiesIfNull)(QFFModel* qffModel, double temperature);
	void (*printModesAndVelocities)(QFFModel* qffModel,FILE* file);
	void (*computeQFFParameters)(QFFModel* qffModel);
};

QFFModel newQFFModel();

#endif /* __QFF_Model_H__ */
