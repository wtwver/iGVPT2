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

#ifndef INTERFACETMCCHEMI_H
#define INTERFACTMECCHEMI_H


#include "../Molecule/Molecule.h"

#ifdef ENABLE_PYTHON
#include <Python.h>
typedef struct _InterfaceTM InterfaceTM;

struct _InterfaceTM
{
	PyObject *pName;
	PyObject *pModule;
	PyObject *pGetManager;
	PyObject *pGetEnergy;
	PyObject *pGetEnergyAndForces;
	PyObject *pOptGeom;
	PyObject *pSetData;
	PyObject *pSetCoordinates;
	PyObject *pCloseSession;
	PyObject *pManager;
	PyObject *pa;
	int initialized;
};
InterfaceTM *newInterfaceTM(char* tmModule);
int interfaceTMComputeGradients(InterfaceTM* interfaceTM, Molecule* mol);
int interfaceTMComputeEnergy(InterfaceTM* interfaceTM, Molecule* mol);
int interfaceTMOpt(InterfaceTM* interfaceTM, Molecule* mol);
void interfaceTMDestroy(InterfaceTM* interfaceTM);

#endif
int runTM(Molecule* mol, char* moduleName, int computeGradients);
int runTensorMol(Molecule* mol, char* moduleName, int computeGradients);
int runOptTensorMol(Molecule* mol, char* moduleName);


#endif
