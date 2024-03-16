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

/* VPT2Model.h */

#ifndef __VPT2_Model_H__
#define __VPT2_Model_H__
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "VPT2PotentialData.h"
#include "VPT2PropertiesData.h"

typedef struct 
{
	int nFrequencies;
	double *harmonicFrequencies;
	double *fundamentals;
	double *overtones;
	double **combinationBands;
	double **Xi;

}VAnharmonic;

typedef struct 
{
	int nFrequencies;
	double *harmonic;
	double *fundamentals;
	double *overtones;
	double **combinationBands;
	double maxFrequencyDifference11Resonance;
	double thresholds11Numerators[2];
	double parameters11Resonance[3];

}PropertiesAnharmonic;



typedef struct _VPT2Model        VPT2Model;
typedef struct _VPT2ModelClass   VPT2ModelClass;

struct _VPT2Model
{
	VAnharmonic vAnharmonic;
	PropertiesAnharmonic pAnharmonic;
	VPT2PotentialData vData;
	VPT2PropertiesData pData;
        VPT2ModelClass* klass;
};

struct _VPT2ModelClass
{
        void (*readData)(VPT2Model* vpt2Model, char* fileName);
        void (*computeAnharmonic)(VPT2Model* vpt2Model);
};

VPT2Model newVPT2Model();

#endif /* __VPT2_Model_H__ */
