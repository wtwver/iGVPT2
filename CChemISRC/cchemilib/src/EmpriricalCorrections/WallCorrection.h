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

#ifndef __CCHEMILIB_WALLCORRECTION_H__
#define __CCHEMILIB_WALLCORRECTION_H__

#include "../Utils/Types.h"
#include "../Molecule/Molecule.h"

/* Spherical wall E0 *(1-exp(-r^2/rho^2)**n */
double getWallCorrection(Molecule* mol, boolean addGradient);

#endif /* __CCHEMILIB_WALLCORRECTION_H__ */

