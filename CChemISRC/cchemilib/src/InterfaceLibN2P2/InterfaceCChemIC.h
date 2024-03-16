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

#ifndef INTERFACECCHEMIC_H
#define INTERFACECCHEMIC_H

#include "../Molecule/Molecule.h"

#include "DerivativesIC.h"

/* Energy and forces */
void *newInterfaceCChemIDef();
void *newInterfaceCChemI(char* directory, double  cflength, double  cfenergy, int myRank);
void interfaceCChemIInitialize(void* interface, char* directory, double cflength, double cfenergy, int myRank);
int interfaceCChemIComputeForces(void* interface, int nAtoms,  char** symbols, double** coordinates, double** forces, double* energy);
int interfaceCChemIGetEnergy(void* interface, int nAtoms, char** symbols, double** coordinates, double* energy);
int interfaceCChemIComputeGradients(void* interface, Molecule* mol);
int interfaceCChemIComputeEnergy(void* interface, Molecule* mol);
int interfaceCChemIComputeHessian(void* interface, Molecule* mol, double** hessian);
#ifdef HIGH_DERIVATIVES
DerivativesIC* interfaceCChemIComputeHighDerivatives(void* interface, Molecule* mol, int order, int method);
#endif
void interfaceCChemIDestroy(void *interface);
int runN2P2(Molecule* mol, char* directory, double cflength,double cfenergy, int myRank, int computeGradient);

/* charges and dipoles */

void *newInterfaceCChemIESDef();
void *newInterfaceCChemIES(char* directory, double       cflength, double       cfdipole, int          myRank);
void interfaceCChemIESInitialize(void* interface, char* directory, double       cflength, double       cfdipole, int          myRank);
int interfaceCChemIESGetChargeAndDipole(void* interface, int nAtoms, char** symbols, double** coordinates, double* charge, double* dipole);
int interfaceCChemIESComputeChargeAndDipole(void* interface, Molecule* mol, double* charge);
int interfaceCChemIESComputedDipole(void* interface, Molecule* mol, double*** dmu, double* charge);
#ifdef HIGH_DERIVATIVES
int interfaceCChemIESComputed2Dipole(void* interface, Molecule* mol, double*** dmu, double* charge);
DerivativesIC** interfaceCChemIComputeDipoleHighDerivatives(void* interface, Molecule* mol, int order, int method);
#endif
void interfaceCChemIESDestroy(void *interface);
int runN2P2ES(Molecule* mol, char* directory, double cflength,double cfdipole, double* charge, int myRank);

#endif
