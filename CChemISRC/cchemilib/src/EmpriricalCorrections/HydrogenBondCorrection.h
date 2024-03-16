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

#ifndef __CCHEMILIB_HYDROGENBONDCORRECTIONS_H__
#define __CCHEMILIB_HYDROGENBONDCORRECTIONS_H__

/* Reference: J. Rezac, P. Hobza J. Chem. Theory Comput. 8, 141-151 (2012)*/

#include "../Utils/Types.h"
#include "../Molecule/Molecule.h"

typedef struct _HyhrogenBondCorrectionParameters  HyhrogenBondCorrectionParameters;
//==============================================================================
// Constants
//==============================================================================
typedef enum
{
  HYDROGEN = 1,
  CARBON = 6,
  NITROGEN = 7,
  OXYGEN = 8
} HAtomsTypes;

struct _HyhrogenBondCorrectionParameters
{
	char method[40];
/* H4 correction*/
	double para_OH_O;
	double para_OH_N;
	double para_NH_O;
	double para_NH_N;
	double multiplier_WH_O;
	double multiplier_NH4;
	double multiplier_COO;
/* HH repulsion*/
	double HH_REPULSION_k;
	double HH_REPULSION_e;
	double HH_REPULSION_r0;
/*Cutoff for the correction, more distant donor-acceptor pairs do not contribute*/
	double HB_R_CUTOFF;
/* Short-range cutoff, closer donor-acceptor pairs do not contribute*/
	double HB_R_0;
/* max. X-H covalent bond distance*/
	double MAX_XH_BOND;
};
/*H4+HH*/
double getH4Correction(Molecule* molecule, HyhrogenBondCorrectionParameters* parameters, boolean addGradient);
int readHydrogenBondCorrectionParameters(HyhrogenBondCorrectionParameters* parameters, char* fileName);
int setHydrogenBondCorrectionParameters(HyhrogenBondCorrectionParameters* parameters, char* fileName, char* method);

#endif /* __CCHEMILIB_HYDROGENBONDCORRECTIONS_H__ */

