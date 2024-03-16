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

#ifndef __CCHEMI_PDBTEMPLATE_H__
#define __CCHEMI_PDBTEMPLATE_H__

#include "../Molecule/Molecule.h"

typedef struct _PDBTemplate  PDBTemplate;
typedef struct _PDBResidueTemplate  PDBResidueTemplate;
typedef struct _PDBTypeTemplate  PDBTypeTemplate;

/************************************/
struct _PDBTypeTemplate
{
	char* pdbType;
	char* mmType;
	double charge;
};
/************************************/
struct _PDBResidueTemplate
{
	char* residueName;
	int numberOfTypes;
	PDBTypeTemplate* typeTemplates;
};
/************************************/
struct _PDBTemplate
{
	int numberOfResidues;
	PDBResidueTemplate* residueTemplates;
};
PDBTemplate* LoadPersonalPDBTpl();
PDBTemplate* LoadDefaultPDBTpl();
void LoadPDBTpl();
char* getMMTypeFromPDBTpl(char* residueName,char* pdbType,double* charge);
int getHydrogensFromPDBTpl(char* residueName,char* pdbType, char** hAtoms);
PDBTemplate* getPointerPDBTemplate();
void setPointerPDBTemplate(PDBTemplate* ptr);
char** getListPDBTypes(char* residueName, int* nlist);
void buildMMTypesFromPDB(Molecule* mol, boolean setCharge);

#endif /* __CCHEMI_PDBTEMPLATE_H__ */

