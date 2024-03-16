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

/* TreeMolecule.c */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <ctype.h>

#include "../Utils/Utils.h"
#include "../Utils/AtomsProp.h"
#include "../Utils/Constants.h"
#include "../Utils/Types.h"
#include "../Utils/QL.h"
#include "../Molecule/Molecule.h"
#include "../Utils/CalculTypesAmber.h"
#include "../Utils/PDBTemplate.h"
#include "../Molecule/TreeMolecule.h"


/************************************************************************/
static void freeStack(TreeMolecule* treeMolecule)
{
	if(treeMolecule->inStack) free(treeMolecule->inStack);
	treeMolecule->inStack = NULL;
}
/************************************************************************/
static void initStack(TreeMolecule* treeMolecule)
{
	int i;
	/* if(treeMolecule->inStack) freeStack(treeMolecule);*/
	treeMolecule->inStack = malloc(treeMolecule->nAtoms*sizeof(boolean));
	for(i=0;i<treeMolecule->nAtoms;i++) treeMolecule->inStack[i] = FALSE;
}
/************************************************************************/
static void freeConnections(TreeMolecule* treeMolecule)
{
	int i;
	if(!treeMolecule->connected) return;
	for(i=0;i<treeMolecule->nAtoms;i++)
	{
		if(treeMolecule->connected[i]) free(treeMolecule->connected[i]);
	}
	free(treeMolecule->connected);
	treeMolecule->connected = NULL;
	treeMolecule->nAtoms = 0;
}
/************************************************************************/
static void disconnect(TreeMolecule* treeMolecule, int n1, int n2)
{
	int i;
	int k;
	if(treeMolecule->nAtoms<1) return;
	for(k=1;k<=treeMolecule->connected[n1][0];k++)
	{
		if(treeMolecule->connected[n1][k]==n2)
		{
			for(i=k;i<treeMolecule->connected[n1][0];i++)
				treeMolecule->connected[n1][i]=treeMolecule->connected[n1][i+1];
			treeMolecule->connected[n1][0]--;
			break;
		}
	}
	for(k=1;k<=treeMolecule->connected[n2][0];k++)
	{
		if(treeMolecule->connected[n2][k]==n1)
		{
			for(i=k;i<treeMolecule->connected[n2][0];i++)
				treeMolecule->connected[n2][i]=treeMolecule->connected[n2][i+1];
			treeMolecule->connected[n2][0]--;
			break;
		}
	}
}
/************************************************************************/
static void initConnections(TreeMolecule* treeMolecule, Molecule*  mol)
{
	int i;
	int j;
	int k;
	int nj;
	if(treeMolecule->nAtoms<1) return;
	treeMolecule->connected = malloc(treeMolecule->nAtoms*sizeof(int*));
	for(i=0;i<treeMolecule->nAtoms;i++)
	{
		treeMolecule->connected[i] = malloc((treeMolecule->nAtoms+1)*sizeof(int));
		treeMolecule->connected[i][0] = 0;
	}
	for(i=0;i<treeMolecule->nAtoms;i++)
	if(mol->atoms[i].typeConnections)
	for(j=0;j<treeMolecule->nAtoms;j++)
	{
		if(i==j) continue;
		nj = mol->atoms[j].N-1;
		if(mol->atoms[i].typeConnections[nj]>0)
		{
			treeMolecule->connected[i][0]++;
			k = treeMolecule->connected[i][0];
			treeMolecule->connected[i][k]=j;
		}
	}
}
/************************************************************************/
void initTreeMolecule(TreeMolecule* treeMolecule, Molecule*  mol, int ringSize)
{
	treeMolecule->done = FALSE;
	treeMolecule->bonds = 0;
	treeMolecule->ringSize = ringSize;
	treeMolecule->nAtoms = mol->nAtoms;
	if(mol->nAtoms<1) return;
	initStack(treeMolecule);
	initConnections(treeMolecule, mol);
}
/************************************************************************/
void freeTreeMolecule(TreeMolecule* treeMolecule)
{
	freeConnections(treeMolecule);
	freeStack(treeMolecule);
}
/************************************************************************/
boolean inRingTreeMolecule(TreeMolecule* treeMolecule,int currentAtom, int rootAtom)
{
	int i;
	treeMolecule->inStack[currentAtom] = TRUE;
	if (treeMolecule->done) return TRUE;
	else if ( ( currentAtom == rootAtom ) && ( treeMolecule->bonds == treeMolecule->ringSize ) ) return TRUE;
	else if ( ( currentAtom == rootAtom ) && ( treeMolecule->bonds > 2 ) && ( treeMolecule->ringSize < 3 ) ) return TRUE;
	if ( treeMolecule->bonds < treeMolecule->ringSize )
	{
		int numberOfConnections = treeMolecule->connected[ currentAtom ][ 0 ];
		for (i = 1; i <= numberOfConnections; i++ )
		{
			int newAtom = treeMolecule->connected[currentAtom][i];
			if ( ! ( treeMolecule->inStack[newAtom] ) )
			{
				treeMolecule->bonds++;
				treeMolecule->done = inRingTreeMolecule(treeMolecule, newAtom, rootAtom);
			}
			if (treeMolecule->done) return TRUE;
		}
	}
	treeMolecule->inStack[currentAtom] = FALSE;
	treeMolecule->bonds--;
	return FALSE;
}
/************************************************************************/
static boolean isConnected(TreeMolecule* treeMolecule, int i, int j)
{
	int k;
	for(k=0;k<treeMolecule->connected[i][0];k++)
		if(treeMolecule->connected[i][k+1]==j) return TRUE;
	return FALSE;
}
/************************************************************************/
int* getRingTreeMolecule(TreeMolecule* treeMolecule)
{
	int i;
	int n= 0;
 	int* ringAtoms = NULL;
	int k;
	int j;
	if(treeMolecule->inStack && treeMolecule->ringSize>1)
	{
 		ringAtoms = malloc(treeMolecule->ringSize*sizeof(int));
		for(i=0;i<treeMolecule->nAtoms;i++)
		{
			if(treeMolecule->inStack[i])
			{
				ringAtoms[n] = i;
				n++;
				if(n>=treeMolecule->ringSize) break;
			}
		}
	}
	if(ringAtoms)
	{
		for(i=1;i<n;i++)
			if(ringAtoms[i]<ringAtoms[0])
			{
				int t = ringAtoms[i];
				ringAtoms[i] = ringAtoms[0];
				ringAtoms[0] = t;
			}

		for(i=0;i<n-2;i++)
		{
			k = i+1;
			for(j=i+1;j<n;j++)
				if(isConnected(treeMolecule,ringAtoms[i],ringAtoms[j]))
				{
					k = j;
					break;
				}
			if(k!=(i+1))
			{
				int t = ringAtoms[i+1];
				ringAtoms[i+1] = ringAtoms[k];
				ringAtoms[k] = t;
			}
		}
	}
	return ringAtoms;
}
/********************************************************************************/
void getCentreRingTreeMolecule(TreeMolecule* treeMolecule, Molecule* mol,int i, int j, double C[])
{

	int k;
	int c;
	int* num  = NULL;
	int n;

	for(c=0;c<3;c++) C[c] = 0;
	if(mol->nAtoms != treeMolecule->nAtoms) return;
	n = 4;
	initTreeMolecule(treeMolecule, mol, n-1);
	if(inRingTreeMolecule(treeMolecule,j, i) ) num = getRingTreeMolecule(treeMolecule);
	if(!num)
	{
		n++;
		initTreeMolecule(treeMolecule, mol, n-1);
		if(inRingTreeMolecule(treeMolecule,j, i)) num = getRingTreeMolecule(treeMolecule);
	}
	if(!num)
	{
		n++;
		initTreeMolecule(treeMolecule, mol, n-1);
		if(inRingTreeMolecule(treeMolecule,j, i)) num = getRingTreeMolecule(treeMolecule);
	}
	if(num)
	{
		for(c=0;c<3;c++) C[c] += mol->atoms[j].coordinates[c];
		for(k=0;k<n-1;k++) 
		{
			for(c=0;c<3;c++) C[c] += mol->atoms[num[k]].coordinates[c];
		}
		for(c=0;c<3;c++) C[c] /= n;
		if(num)free(num);
	}
}
/********************************************************************************/
boolean inGroupTreeMolecule(TreeMolecule* treeMolecule,int currentAtom, int nEx1, int nEx2, int nEx3)
{
	int i;
	//int end = FALSE;
	treeMolecule->inStack[currentAtom] = TRUE;
	if ( currentAtom == nEx1  || currentAtom == nEx2  ||  currentAtom == nEx2 )
	{
		treeMolecule->done = TRUE;
		return TRUE;
	}
	{
		int numberOfConnections = treeMolecule->connected[ currentAtom ][ 0 ];
		for (i = 1; i <= numberOfConnections; i++ )
		{
			int newAtom = treeMolecule->connected[currentAtom][i];
			if ( ! ( treeMolecule->inStack[newAtom] ) )
			{
				treeMolecule->bonds++;
				//end = inGroupTreeMolecule(treeMolecule, newAtom, nEx1, nEx2, nEx3);
				inGroupTreeMolecule(treeMolecule, newAtom, nEx1, nEx2, nEx3);
			}
		}
	}
	return FALSE;
}
/************************************************************************/
int* getListGroupe(int* nGroupAtoms, Molecule*  mol, int i1, int i2, int i3, int i4)
{
	int i;
	int nG = 0;
	int nEx1 = 0;
	int nEx2 = 0;
	int nEx3 = 0;
	TreeMolecule treeMolecule;
	int* listGroupAtoms = NULL;
	int n = 0;

	*nGroupAtoms = 0;
	if(mol->nAtoms<2) return NULL;
	if(i1<0 || i2<0) return NULL;
	initTreeMolecule(&treeMolecule, mol, 6);
	if(i3>=0 && i4>=0) 
	{
		if(i4>=mol->nAtoms)
		{
			nG = i3;
			nEx1 = i2;
			nEx2 = i1;
			nEx3 = i4-mol->nAtoms;
		}
		else
		{
			nG = i4;
			nEx1 = i3;
			nEx2 = i2;
			nEx3 = i1;
		}
	}
	else if(i3>=0) 
	{
		nG = i3;
		nEx1 = i2;
		nEx2 = i1;
		nEx3 = i4;
	}
	else 
	{
		nG = i2;
		nEx1 = i1;
		nEx2 = i3;
		nEx3 = i4;
	}
	disconnect(&treeMolecule, nG, nEx1);
	inGroupTreeMolecule(&treeMolecule,nG, nEx1, nEx2, nEx3);
	if(treeMolecule.done) return NULL;
	
	/*
	printf("end = %d\n",treeMolecule.done);
	printf("nex = %d n = %d\n",mol->atoms[nEx1].N, mol->atoms[nG].N);
	for(i=0;i<treeMolecule.nAtoms;i++) printf("%d %d\n",mol->atoms[i].N,treeMolecule.inStack[i]);
	*/

	for(i=0;i<treeMolecule.nAtoms;i++) 
		if(treeMolecule.inStack[i] && i!=nG) n++;
	if(n==0) return NULL;
	listGroupAtoms = malloc(n*sizeof(int));
	for(i=0;i<n;i++) listGroupAtoms[i]=-1;
	n = 0;
	for(i=0;i<treeMolecule.nAtoms;i++) 
		if(treeMolecule.inStack[i] && i!=nG) listGroupAtoms[n++]=i;

	freeTreeMolecule(&treeMolecule);

	*nGroupAtoms = n;
	return listGroupAtoms;
}
