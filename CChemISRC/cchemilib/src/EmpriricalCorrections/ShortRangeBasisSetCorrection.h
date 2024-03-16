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

#ifndef __CCHEMILIB_SRBCORRECTIONS_H__
#define __CCHEMILIB_SRBCORRECTIONS_H__

/* Reference: S. Grimme The Journal of Chemical Physics 148, 064104 (2018); https://doi.org/10.1063/1.5012601 */

#include "../Utils/Types.h"
#include "../Molecule/Molecule.h"

typedef struct _ShortRangeBasisSetCorrectionParameters  ShortRangeBasisSetCorrectionParameters;
//==============================================================================
// Constants
//==============================================================================
// A*exp(beta*rij)*erfc(w x)+C*exp(gamma*rij)*erf(w x)
// erfc(x) = 1 - erf(x)
// derive de erf(x)= 2/sqrt(M_PI)*exp(-x*x)
typedef struct _SRBBond
{
	char symbol1[10];
	char symbol2[10];
	double A;
	double beta;
	double C;
	double gamma;
        double omega;
} SRBBond;

struct _ShortRangeBasisSetCorrectionParameters
{
	char method[40];
	int nBonds;
	SRBBond* sRBBonds;
};

double getSRBCorrection(Molecule* molecule, ShortRangeBasisSetCorrectionParameters* parameters, boolean addGradient);
int readShortRangeBasisSetCorrectionParameters(ShortRangeBasisSetCorrectionParameters* parameters, char* fileName);
int setShortRangeBasisSetCorrectionParameters(ShortRangeBasisSetCorrectionParameters* parameters, char* fileName, char* method);

#endif /* __CCHEMILIB_SRBCORRECTIONS_H__ */

