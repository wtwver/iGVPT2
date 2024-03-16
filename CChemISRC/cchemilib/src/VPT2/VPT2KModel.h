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

/* VPT2KModel.h */

#ifndef __VPT2_K_Model_H__
#define __VPT2_K_Model_H__
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "../VPT2/VPT2PotentialData.h"
#include "../VPT2/VPT2PropertiesData.h"
#include "../Utils/Utils.h"
#include "../Utils/QL.h"
#include "../Utils/Constants.h"


typedef struct 
{
	int nFrequencies;
	int nStates;
	double *harmonicEnergies;
	double **Xi;
	double*** D;
	int** v;
	double C;
	double **H;
	double *eigenValues;
	double **eigenVectors;
}VKAnharmonic;

typedef struct 
{
	int nFrequencies;
	int nStates;
	double *harmonic;
	double *anHarmonic;
	double maxFrequencyDifference11Resonance;
	double thresholds11Numerators[2];
	double parameters11Resonance[3];

}PropertiesKAnharmonic;


typedef struct _VPT2KModel        VPT2KModel;
typedef struct _VPT2KModelClass   VPT2KModelClass;

struct _VPT2KModel
{
	VKAnharmonic vAnharmonic;
	PropertiesKAnharmonic pAnharmonic;
	VPT2PotentialData vData;
	VPT2PropertiesData pData;
        VPT2KModelClass* klass;
};

struct _VPT2KModelClass
{
        void (*readData)(VPT2KModel* vpt2KModel, char* fileName);
        void (*computeAnharmonic)(VPT2KModel* vpt2KModel);
};

VPT2KModel newVPT2KModel();

#endif /* __VPT2_K_Model_H__ */
