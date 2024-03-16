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

/* VPT2PotentialData.h */

#ifndef __VPT2_PotentialData_H__
#define __VPT2_PotentialData_H__
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "../Utils/Utils.h"
#include "../Utils/QL.h"
#include "../Utils/Constants.h"

typedef struct _VPT2PotentialData        VPT2PotentialData;
typedef struct _VPT2PotentialDataClass   VPT2PotentialDataClass;

typedef enum
{
  MODEL_VPT2 = 0,/* pure VPT2 */
  MODEL_DCPT2, /* DCPT2 method JCTC Bloino, (2012) 8, 1015  */
  MODEL_HDCPT2, /* HDCPT2 hybrid DCPT2 method JCTC Bloino, (2012) 8, 1015  */
  MODEL_GVPT2, /* VPT2/GVPT2 of Gaussian : remove resonnace terms in intensities calculation  */
  MODEL_VPT2K
} VPT2ModelType;

typedef struct 
{
	VPT2ModelType type;
	double alphaHDCPT2;
	double betaHDCPT2;
}VPT2VModel;

struct _VPT2PotentialData
{
	int nFrequencies;
	double *gradients;
	double **hessian;
	double ***cubic;
	double ****quartic;
	double ***coriolis;
	double *Be;
	VPT2VModel model;
        double maxFrequencyDifferenceFermi;
        double parametersResonance[3];/* 0 and 1 for I and II Martin, 2 for Z Krasnoshchekov */
        VPT2PotentialDataClass* klass;
};

struct _VPT2PotentialDataClass
{
        void (*readData)(VPT2PotentialData* vpt2PotentialData, char* fileName);
};

VPT2PotentialData newVPT2PotentialData(int n);

#endif /* __VPT2_PotentialData_H__ */
