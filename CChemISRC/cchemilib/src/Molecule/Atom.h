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

#ifndef __CCHEMILIB_ATOM_H__
#define __CCHEMILIB_ATOM_H__

#include "../Utils/Types.h"

typedef struct _Atom
{
	double coordinates[3];
	double gradient[3];
	double velocity[3];
	double charge; // partial charge
	double charge0; // oxidation
	double electronegativity; // chi
	double hardness; // eta
	double width; // eta
	double mass;
	SAtomsProp prop;
	char* mmType;
	char* pdbType;
	char* residueName;
	int residueNumber;
	int N;

	boolean show;
	boolean variable;
	CCHEMILayerType layer;
	int* typeConnections;
	double rho;
	double U;// for ACKS2
}Atom;

double getDistance(Atom *a1,Atom* a2);
double getAngle(Atom *a1,Atom* a2,Atom* a3);
double getTorsion(Atom* a1,Atom* a2,Atom* a3,Atom* a4);
Atom getCopyAtom(Atom* atom);

#endif /* __CCHEMILIB_ATOM_H__ */

