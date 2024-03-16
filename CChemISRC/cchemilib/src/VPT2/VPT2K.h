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

/* VPT2Potential.h */

#ifndef __VPT2_K_H__
#define __VPT2_K_H__
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

typedef struct 
{
	int nFrequencies;
	double *gradients;
	double **hessian;
	double ***cubic;
	double ****quartic;
	double ***coriolis;
	double *Be;
	VPT2Model model;
}VData;

typedef struct 
{
	int nFrequencies;
	double *harmonicFrequencies;
	double *fundamentals;
	double *overtones;
	double **combinationBands;
	double **Xi;

}VAnharmonic;



typedef struct _VPT2Potential        VPT2Potential;
typedef struct _VPT2PotentialClass   VPT2PotentialClass;

struct _VPT2Potential
{
	VData data;
	VAnharmonic anharmonic;
        VPT2PotentialClass* klass;

};

struct _VPT2PotentialClass
{
        void (*readData)(VPT2Potential* vpt2Potential, char* fileName);
        void (*computeAnharmonic)(VPT2Potential* vpt2Potential);
};

VPT2Potential newVPT2Potential();

#endif
