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

/*HydrogenBond.c*/
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include "../Utils/Types.h"
#include "../Utils/Utils.h"

static int nAtomsCanConnect = 6;
static char** atomsCanConnect = NULL;
static double minDistance = 1.50; /* in Agnstrom */
static double maxDistance = 3.15; /* in Agnstrom */ 
static double minAngle = 145.0;
static double maxAngle = 215.0;

/************************************************************************/
void initHBondsProperties()
{
	nAtomsCanConnect = 6;
	atomsCanConnect = malloc(nAtomsCanConnect*sizeof(char*));
	atomsCanConnect[0] = strdup("N");
	atomsCanConnect[1] = strdup("O");
	atomsCanConnect[2] = strdup("F");
	atomsCanConnect[3] = strdup("Cl");
	atomsCanConnect[4] = strdup("Br");
	atomsCanConnect[5] = strdup("I");
	minDistance = 1.50; /* in Agnstrom */
	maxDistance = 3.15; /* in Agnstrom */ 
	minAngle = 145.0;
	maxAngle = 215.0;
}
/******************************************************************/
void saveHBondsProperties()
{
	char *hbondsfile;
	FILE *file;
	int i;

	hbondsfile = strdup_printf("%s%shbonds.txt",cchemiDirectory(),DIR_SEPARATOR_S);

	file = fopen(hbondsfile, "w");
	if(!file) return;

 	fprintf(file,"%f\n",minDistance);
 	fprintf(file,"%f\n",maxDistance);
 	fprintf(file,"%f\n",minAngle);
 	fprintf(file,"%f\n",maxAngle);
 	fprintf(file,"%d\n",nAtomsCanConnect);
	for(i=0;i<nAtomsCanConnect;i++) fprintf(file,"%s\n",atomsCanConnect[i]);
	fclose(file);

	free(hbondsfile);
}
/******************************************************************/
void readHBondsProperties()
{
	char *hbondsfile;
	FILE *file;
	int n;
	int i;

	initHBondsProperties();
	hbondsfile = strdup_printf("%s%shbonds.txt",cchemiDirectory(),DIR_SEPARATOR_S);

	file = fopen(hbondsfile, "r");
	if(!file) return;

 	n = fscanf(file,"%lf\n",&minDistance);
	if(n != 1) { initHBondsProperties(); return ; fclose(file); free(hbondsfile);}

 	n = fscanf(file,"%lf\n",&maxDistance);
	if(n != 1) { initHBondsProperties(); return ; fclose(file); free(hbondsfile);}

 	n = fscanf(file,"%lf\n",&minAngle);
	if(n != 1) { initHBondsProperties(); return ; fclose(file); free(hbondsfile);}
 	n = fscanf(file,"%lf\n",&maxAngle);
	if(n != 1) { initHBondsProperties(); return ; fclose(file); free(hbondsfile);}

 	n = fscanf(file,"%d\n",&nAtomsCanConnect);
	if(n != 1 || nAtomsCanConnect<0 ) { initHBondsProperties(); return ; fclose(file); free(hbondsfile);}

	for(i=0;i<nAtomsCanConnect;i++)
	{
 		n = fscanf(file,"%s\n",atomsCanConnect[i]);
		if(n != 1) { initHBondsProperties(); return ; fclose(file); free(hbondsfile);}
		deleteLastSpaces(atomsCanConnect[i]);
		deleteFirstSpaces(atomsCanConnect[i]);
		strDeleten(atomsCanConnect[i]);
	}

	fclose(file);

	free(hbondsfile);
}
/************************************************************************/
double getMinDistanceHBonds()
{
	return minDistance;
}
/************************************************************************/
double getMaxDistanceHBonds()
{
	return maxDistance;
}
/************************************************************************/
double getMinAngleHBonds()
{
	return minAngle;
}
/************************************************************************/
double getMaxAngleHBonds()
{
	return maxAngle;
}
/************************************************************************/
boolean atomCanDoHydrogenBond(char* symbol)
{
	int k;
	for(k=0;k<nAtomsCanConnect;k++) if(strcmp(symbol,atomsCanConnect[k])==0) return TRUE;
	return FALSE;
}
/************************************************************************/
