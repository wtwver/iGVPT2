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

#ifndef __CCHEMILIB_DISPERSIONCORRECTION_H__
#define __CCHEMILIB_DISPERSIONCORRECTION_H__

/* Grimme at al, JCP, 132, 154104(2010) */

#include "../Utils/Types.h"
#include "../Molecule/Molecule.h"

typedef struct _DispersionParameters  DispersionParameters;
//==============================================================================
// Constants
//==============================================================================

struct _DispersionParameters
{
	char method[40];
	double k1;
	double k2;
	double k3;
	double rs6;
	double alp6;
	double rs8;
	double alp8;
	double s6;
	double s8;
	double rthr;
};
double getD3Correction(Molecule* molecule, char* method, boolean addGradient);

#endif /* __CCHEMILIB_DISPERSIONCORRECTION_H__ */

