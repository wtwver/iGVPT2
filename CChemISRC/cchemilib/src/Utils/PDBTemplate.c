/********************************************************************************
 cchemi is an interface to ab initio computational chemistry programs 
 designed for add them many functionalities non available in these packages.
 Copyright (C) 2010 Abdulrahman Allouche (University Lyon 1)

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

/* PDBTemplate.c */
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>

#include "../Utils/Types.h"
#include "../Utils/Utils.h"
#include "../Utils/PDBTemplate.h"
#include "../Utils/LoadPDBTemplate.h"
#include "../Utils/CreateDefaultPDBTpl.h"
#include "../Utils/SavePDBTemplate.h"
#include "../Molecule/Molecule.h"

static	PDBTemplate* staticPDBTemplate = NULL;
/************************************************************/
PDBTemplate* freePDBTpl(PDBTemplate* pdbTemplate)
{
	int i;
	int j;
	if(!pdbTemplate)
		return NULL;
	for(i=0;i<pdbTemplate->numberOfResidues;i++)
	{
		if(pdbTemplate->residueTemplates[i].residueName)
			free(pdbTemplate->residueTemplates[i].residueName);

		for(j=0;j<pdbTemplate->residueTemplates[i].numberOfTypes;j++)
		{
			if(pdbTemplate->residueTemplates[i].typeTemplates[j].pdbType)
				free(pdbTemplate->residueTemplates[i].typeTemplates[j].pdbType);
			if(pdbTemplate->residueTemplates[i].typeTemplates[j].mmType)
				free(pdbTemplate->residueTemplates[i].typeTemplates[j].mmType);
		}
		if(pdbTemplate->residueTemplates[i].typeTemplates)
			free(pdbTemplate->residueTemplates[i].typeTemplates);
	}
	if(pdbTemplate->residueTemplates)
		free(pdbTemplate->residueTemplates);
	free(pdbTemplate);
	return NULL;
}
/************************************************************/
PDBTemplate* LoadPersonalPDBTpl()
{
	char* filename = strdup_printf("%s%sPersonalPDBTemplate.tpl",
			cchemiDirectory(), DIR_SEPARATOR_S);
	PDBTemplate* pdbTemplate = NULL;

	pdbTemplate = malloc(sizeof(PDBTemplate));
	if(!readPDBTemplate(pdbTemplate,filename))
	{
		free(pdbTemplate);
		pdbTemplate = NULL;
	}

	free(filename);
	return pdbTemplate;
}
/************************************************************/
PDBTemplate* LoadDefaultPDBTpl()
{
	PDBTemplate* pdbTemplate = NULL;
	char* filename = strdup_printf("%s%sDefaultPDBTemplate.tpl",
			cchemiDirectory(), DIR_SEPARATOR_S);

	pdbTemplate = malloc(sizeof(PDBTemplate));
	if(!readPDBTemplate(pdbTemplate,filename))
	{
		free(pdbTemplate);
		pdbTemplate = NULL;
	}

	free(filename);
	return pdbTemplate;
}
/************************************************************/
void LoadPDBTpl()
{
	if(staticPDBTemplate)
		staticPDBTemplate = freePDBTpl(staticPDBTemplate);
	staticPDBTemplate = LoadPersonalPDBTpl();
	if(!staticPDBTemplate)
	{
		staticPDBTemplate = LoadDefaultPDBTpl();
		if(!staticPDBTemplate)
		{

			if(CreateDefaultPDBTpl())
				staticPDBTemplate = LoadDefaultPDBTpl();
		}
	}
}
/************************************************************/
static int getResiduePDBTplNumber(char* residueName)
{
	int i;
	if(!staticPDBTemplate)
		return -1;

	for(i=0;i<staticPDBTemplate->numberOfResidues;i++)
	{
		if(!strcmp(staticPDBTemplate->residueTemplates[i].residueName,residueName))
			return i;
	}
	return -1;
}
/************************************************************/
static char* getmmType(int residueNumber, char* pdbType,double* charge)
{
	int j;
	PDBTypeTemplate* typeTemplates = 
		staticPDBTemplate->residueTemplates[residueNumber].typeTemplates;
	int numberOfTypes = staticPDBTemplate->residueTemplates[residueNumber].numberOfTypes;
	char* mmType = strdup("UNK");
	for(j=0;j<numberOfTypes;j++)
	{
		if(!strcmp(pdbType,typeTemplates[j].pdbType))
		{
			free(mmType);
			mmType = strdup(typeTemplates[j].mmType);
			*charge = typeTemplates[j].charge;
			return mmType;
		}
	}
	return mmType;
}
/************************************************************/
char* getMMTypeFromPDBTpl(char* residueName,char* pdbType,double* charge)
{
	char* mmType = strdup("UNK");
	int residueNumber = -1;
	*charge = 0;

	if(!staticPDBTemplate) return mmType;
	residueNumber = getResiduePDBTplNumber(residueName);
	if(residueNumber==-1)
	{
		residueNumber = getResiduePDBTplNumber("ALLRESIDUE");
		if(residueNumber==-1)
			return mmType;
		else
			return getmmType(residueNumber,pdbType,charge);
	}
	else
	{
		mmType = getmmType(residueNumber,pdbType,charge);
		if(!strcmp(mmType,"UNK"))
		{
			residueNumber = getResiduePDBTplNumber("ALLRESIDUE");
			if(residueNumber==-1)
				return mmType;
			else
				return getmmType(residueNumber,pdbType,charge);
		}
		else
			return mmType;
	}
}
/************************************************************/
static int getHydrogens(int residueNumber, char* pdbType, char** hAtoms)
{
	int j;
	int k;
	PDBTypeTemplate* typeTemplates = 
		staticPDBTemplate->residueTemplates[residueNumber].typeTemplates;
	int numberOfTypes = staticPDBTemplate->residueTemplates[residueNumber].numberOfTypes;
	int nH = 0;
	for(j=0;j<numberOfTypes;j++)
	{
		if(!strcmp(pdbType,typeTemplates[j].pdbType))
		{
			for(k=j+1;k<numberOfTypes;k++)
			{
				if(!typeTemplates[k].pdbType) break;
				if(typeTemplates[k].pdbType[0]!='H') break;
				sprintf(hAtoms[nH],"%s",typeTemplates[k].pdbType);
				nH++;
				if(nH>10) break;
			}
			return nH;
		}
	}
	return nH;
}
/************************************************************/
int getHydrogensFromPDBTpl(char* residueName,char* pdbType, char** hAtoms)
{
	int nH = 0;
	int residueNumber = -1;

	if(!pdbType) return nH;
	if(pdbType && (pdbType[0]=='H' || pdbType[0]=='h')) return nH;
	if(!staticPDBTemplate) return nH;
	residueNumber = getResiduePDBTplNumber(residueName);
	if(residueNumber==-1)
	{
		residueNumber = getResiduePDBTplNumber("ALLRESIDUE");
		if(residueNumber==-1)
			return nH;
		else
			return getHydrogens(residueNumber,pdbType, hAtoms);
	}
	else
	{
		nH = getHydrogens(residueNumber,pdbType, hAtoms);
		if(nH == 0)
		{
			residueNumber = getResiduePDBTplNumber("ALLRESIDUE");
			if(residueNumber==-1)
				return nH;
			else
				return getHydrogens(residueNumber,pdbType, hAtoms);
		}
		else
			return nH;
	}
	return nH;
}
/************************************************************/
PDBTemplate* getPointerPDBTemplate()
{
	return staticPDBTemplate;
}
/************************************************************/
void setPointerPDBTemplate(PDBTemplate* ptr)
{
	staticPDBTemplate = ptr;
}
/************************************************************/
char** getListPDBTypes(char* residueName, int* nlist)
{
	char** t = NULL;
	int j;
	int residueNumber = getResiduePDBTplNumber(residueName);
	PDBTypeTemplate* typeTemplates = NULL;
	int numberOfTypes = 0;

	*nlist = 0;

	if(residueNumber==-1) residueNumber = getResiduePDBTplNumber("ALLRESIDUE");
	if(residueNumber==-1) return NULL;

	typeTemplates = staticPDBTemplate->residueTemplates[residueNumber].typeTemplates;
	numberOfTypes = staticPDBTemplate->residueTemplates[residueNumber].numberOfTypes;
	t = malloc(numberOfTypes*sizeof(char*));
	for(j=0;j<numberOfTypes;j++)
		t[j] = strdup(typeTemplates[j].pdbType);
	*nlist = numberOfTypes;
	return t;
}
/**********************************************************************************************************/
void buildMMTypesFromPDB(Molecule* mol, boolean setCharge)
{
	int i;
	double charge;
	LoadPDBTpl();
	for(i=0;i<mol->nAtoms;i++)
	{
		char* mmType = getMMTypeFromPDBTpl(mol->atoms[i].residueName,mol->atoms[i].pdbType,&charge);
		
		if(mmType)
		{
			if(mol->atoms[i].mmType) free(mol->atoms[i].mmType);
			mol->atoms[i].mmType = mmType;
			if(setCharge) mol->atoms[i].charge = charge;
		}
	}
}
