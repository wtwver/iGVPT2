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

/* SavePDBTemplate.c */
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "../Utils/Types.h"
#include "../Utils/PDBTemplate.h"
#include "../Utils/Utils.h"
/************************************************************/
static void saveTitleResidueTpl(FILE* file)
{
	fprintf(file,"Begin Title\n");
	fprintf(file,"	Residue        : PDB type atom  Amber type atom  charge of atom\n");
	fprintf(file,"End\n");
}
/**********************************************************************/
static void saveResiduesList(PDBTemplate* pdbTemplate, FILE* file)
{
	int i;

	fprintf(file,"Begin Residue List\n");

	for(i=0;i<pdbTemplate->numberOfResidues;i++)
		fprintf(file,"%s\n",pdbTemplate->residueTemplates[i].residueName);

	fprintf(file,"End\n");
}
/**********************************************************************/
static void saveOneResidue(PDBTemplate* pdbTemplate, FILE* file,int residueNumber)
{
	int i;
	int numberOfTypes = pdbTemplate->residueTemplates[residueNumber].numberOfTypes;
	PDBTypeTemplate *typeTemplates = pdbTemplate->residueTemplates[residueNumber].typeTemplates;

	fprintf(file,"Begin %s Residue\n",pdbTemplate->residueTemplates[residueNumber].residueName);
	for(i=0;i<numberOfTypes;i++)
		fprintf(file,"%s %s %f\n",
				typeTemplates[i].pdbType,
				typeTemplates[i].mmType,
				typeTemplates[i].charge
				);
	fprintf(file,"End\n");
}
/**********************************************************************/
static void saveAllResidues(PDBTemplate* pdbTemplate, FILE* file)
{
	int i;
	int n = pdbTemplate->numberOfResidues;
	for(i=0;i<n;i++)
		saveOneResidue(pdbTemplate,file,i);
}
/**********************************************************************/
boolean savePDBTemplate(PDBTemplate* pdbTemplate,char* filename)
{
	FILE* file;
	file = fopen(filename,"w");

	if(file == NULL)
		return FALSE;
	else
	{
		saveTitleResidueTpl(file);
		saveResiduesList(pdbTemplate,file);
		saveAllResidues(pdbTemplate,file);
		fclose(file);
	}
	return TRUE;
}
/**********************************************************************/
