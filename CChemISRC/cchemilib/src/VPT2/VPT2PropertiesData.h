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

/* VPT2PropertiesData.h */

#ifndef __VPT2_PropertiesData_H__
#define __VPT2_PropertiesData_H__

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "../Utils/Utils.h"
#include "../Utils/QL.h"
#include "../Utils/Constants.h"


typedef enum
{
  MODEL_PROP_VPT2 = 0,/* pure VPT2 */
  MODEL_PROP_DCPT2, /* DCPT2 method JCTC Bloino, (2012) 8, 1015  */
  MODEL_PROP_HDCPT2, /* HDCPT2 hybrid DCPT2 method JCTC Bloino, (2012) 8, 1015  */
  MODEL_PROP_GVPT2, /* VPT2/GVPT2 of Gaussian : remove resonnace terms in intensities calculation  */
  MODEL_PROP_GVPT2S, /* VPT2/GVPT2 of Gaussian : remove resonnace terms in intensities calculation  + sum on all basis set */
  MODEL_PROP_VPT2K
} VPT2PropModelType;

typedef struct 
{
	VPT2PropModelType type;
	double alphaHDCPT2;
	double betaHDCPT2;
}VPT2PropModel;

typedef struct 
{
	int nFrequencies;
	int nDim; /* 3 = vector, 6 = tensor */
	double **first;
	double ***second;
	double ****cubic;
	VPT2PropModel model;
}PropertiesData;

typedef struct _VPT2PropertiesData        VPT2PropertiesData;
typedef struct _VPT2PropertiesDataClass   VPT2PropertiesDataClass;

struct _VPT2PropertiesData
{
	int nFrequencies;
	int nDim; /* 3 = vector, 6 = tensor */
	double **first;
	double ***second;
	double ****cubic;
	VPT2PropModel model;
        VPT2PropertiesDataClass* klass;

};

struct _VPT2PropertiesDataClass
{
        void (*readData)(VPT2PropertiesData* vpt2PropertiesData, char* fileName);
};

VPT2PropertiesData newVPT2PropertiesData(int n, int nDim);

#endif /* __VPT2_PropertiesData_H__ */
