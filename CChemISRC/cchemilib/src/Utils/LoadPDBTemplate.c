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

/* LoadPDBTemplate.c */
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>

#include "../Utils/Types.h"
#include "../Utils/PDBTemplate.h"
#include "../Utils/Utils.h"
#include "../Utils/Constants.h"
/**********************************************************************/
char** getResiduesList(FILE* file,int* nResidue)
{
	char** t = (char**)malloc(sizeof(char*));
	
	char dump[BSIZE];
	int len = BSIZE;
	int n;
	boolean Ok = FALSE;

	*nResidue = 0;
	fseek(file, 0L, SEEK_SET);
    	if(fgets(dump,len,file))
	while(!feof(file))
	{
		if(fgets(dump,len,file))
		{
			if(strstr(dump,"Begin Residue List"))
			{
				Ok = TRUE;
				break;
			}
		}
	}
	if(!Ok)
		return NULL;
	n = 0;
	while(!feof(file))
	{
		if(fgets(dump,len,file))
		{
			if(strstr(dump,"End"))
				break;
		}
		t = realloc(t,(n+1)*sizeof(char*));
		t[n] = (char*)malloc(BSIZE*sizeof(char));
		sscanf(dump,"%s",t[n]);
		n++;
	}
	if(n==0)
	{
		free(t);
		return NULL;
	}
	*nResidue = n;

	return t;
}
/**********************************************************************/
void setResiduesList(PDBTemplate* pdbTemplate, FILE* file)
{
	int nResidue;
	char** t = getResiduesList(file,&nResidue);

	pdbTemplate->numberOfResidues = nResidue;
	pdbTemplate->residueTemplates = NULL;
	if(nResidue>0 && t)
	{
		int i;
		pdbTemplate->residueTemplates = malloc(nResidue*sizeof(PDBResidueTemplate));
		for(i=0;i<nResidue;i++)
		{
			pdbTemplate->residueTemplates[i].residueName = strdup(t[i]);
			pdbTemplate->residueTemplates[i].numberOfTypes = 0;
			pdbTemplate->residueTemplates[i].typeTemplates = NULL;
			free(t[i]);
		}
		free(t);
	}
}
/**********************************************************************/
void setOneResidue(PDBTemplate* pdbTemplate, FILE* file,int residueNumber)
{
	char title[BSIZE];
	char pdbType[BSIZE];
	char mmType[BSIZE];
	char charge[BSIZE];
	char dump[BSIZE];
	int len = BSIZE;
	int n = 0;
	PDBTypeTemplate *typeTemplates = NULL;
	boolean Ok = FALSE;

	sprintf(title,"Begin %s Residue",pdbTemplate->residueTemplates[residueNumber].residueName);
	fseek(file, 0L, SEEK_SET);
	while(!feof(file))
	{
		if(fgets(dump,len,file))
		{
			if(strstr(dump,title))
			{
				Ok = TRUE;
				break;
			}
		}
	}
	if(!Ok)
		return;
	n = 0;
	typeTemplates = malloc(sizeof(PDBTypeTemplate));
	while(!feof(file))
	{
		if(fgets(dump,len,file))
		{
			if(strstr(dump,"End"))
				break;
		}
		sscanf(dump,"%s %s %s",pdbType, mmType, charge);
		/*printf("pdbType = %s mmType = %s charge = %s\n",pdbType, mmType, charge);*/
		typeTemplates[n].pdbType = strdup(pdbType);
		typeTemplates[n].mmType = strdup(mmType);
		typeTemplates[n].charge = atof(charge);

		n++;
		typeTemplates = realloc(typeTemplates,(n+1)*sizeof(PDBTypeTemplate));
	}
	if(n==0)
	{
		free(typeTemplates);
		pdbTemplate->residueTemplates[residueNumber].numberOfTypes = 0;
		pdbTemplate->residueTemplates[residueNumber].typeTemplates = NULL;
		return;
	}
	typeTemplates = realloc(typeTemplates,n*sizeof(PDBTypeTemplate));

	pdbTemplate->residueTemplates[residueNumber].numberOfTypes = n;
	pdbTemplate->residueTemplates[residueNumber].typeTemplates = typeTemplates;
}
/**********************************************************************/
void setAllResidues(PDBTemplate* pdbTemplate, FILE* file)
{
	int i;
	int n = pdbTemplate->numberOfResidues;
	/* printf("numberOfResidues = %d\n",pdbTemplate->numberOfResidues);*/
	for(i=0;i<n;i++)
	{
		/* printf("i = %d\n",i);*/
		setOneResidue(pdbTemplate,file,i);
	}
}
/**********************************************************************/
boolean readPDBTemplate(PDBTemplate* pdbTemplate,char* filename)
{
	FILE* file;
	file = fopen(filename,"rb");

	/* printf("Read Default TPL file %s \n",filename);*/
	if(file == NULL)
		return FALSE;
	else
	{
		/* printf("Read List Of residue\n");*/
		setResiduesList(pdbTemplate,file);
		/* printf("End Read List Of residue\n");*/
		setAllResidues(pdbTemplate,file);
		/* printf("End Set All Residue\n");*/
		fclose(file);
	}
	return TRUE;
}
/**********************************************************************/
