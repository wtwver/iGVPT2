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

/* CalculTypesAmber.c */
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "../Utils/Types.h"
#include "../Utils/Constants.h"
#include "../Utils/Timer.h"
#include "../Utils/Utils.h"
#include "../Utils/AtomsProp.h"
#include "../Utils/SList.h"
#include "../Utils/HydrogenBond.h"
#include "../Molecule/Molecule.h"

/**********************************************************************************************************/
typedef struct _MoleculeData MoleculeData;
struct _MoleculeData
{
	Molecule* mol;
	SList* bondsList;
	int* numberOfDoubleBonds;
	int* numberOfTripleBonds;
	int **connected;
	boolean *inStack;
	boolean doneRing;
	int nBondsRing;
};
/************************************************************************/
static boolean bonded(Molecule* mol, int i,int j)
{
	double distance;
	int k;
	double dif[3];
	for(k=0;k<3;k++) dif[k] = mol->atoms[i].coordinates[k] - mol->atoms[j].coordinates[k];
	distance = 0;
	for(k=0;k<3;k++) distance += dif[k]*dif[k];
	distance = sqrt(distance); 
	distance *= ANGTOBOHR;
	if(distance<(mol->atoms[i].prop.covalentRadii + mol->atoms[j].prop.covalentRadii)) return TRUE;
	else return FALSE;
}
/************************************************************************/
static boolean hbonded(Molecule* mol, int nAtoms, int i,int j)
{
	double minDistanceH;
	double maxDistanceH;
	double minDistanceH2;
	double maxDistanceH2;
	double minAngleH;
	double maxAngleH;
	double distance2;
	double angle;
	double A[3];
	double B[3];
	double dif[3];

	int k;
	int kH;
	int kO;
	int c;

	if(!strcmp(mol->atoms[i].prop.symbol,"H"))
	{
		kH = i;
		kO = j;
		if(!atomCanDoHydrogenBond(mol->atoms[j].prop.symbol)) return FALSE;
	}
	else
	{
		if(!strcmp(mol->atoms[j].prop.symbol,"H"))
		{
			kH = j;
			kO = i;
			if(!atomCanDoHydrogenBond(mol->atoms[i].prop.symbol)) return FALSE;
		}
		else return FALSE;
	}
	minDistanceH = getMinDistanceHBonds();
	minDistanceH2 = minDistanceH*minDistanceH;

	maxDistanceH = getMaxDistanceHBonds();
	maxDistanceH2 = maxDistanceH*maxDistanceH;

	minAngleH = getMinAngleHBonds();
	maxAngleH = getMaxAngleHBonds();

	distance2 = 0;
	for(k=0;k<3;k++) distance2 += dif[k]*dif[k];

	if(distance2<minDistanceH2 || distance2>maxDistanceH2) return FALSE;

	for(k=0;k<nAtoms;k++)
	{
		if(k==kH) continue;
		if(k==kO) continue;
		/* angle kO, kH, connection to kH */
		if(!bonded(mol, kH,k)) continue;

		for(c=0;c<3;c++) A[c] = mol->atoms[kO].coordinates[c]-mol->atoms[kH].coordinates[c];
		for(c=0;c<3;c++) B[c] = mol->atoms[k].coordinates[c]-mol->atoms[kH].coordinates[c];

		angle = getAngleVectors(A,B);
		if(angle>=minAngleH &&angle<=maxAngleH) return TRUE;
	}
	return FALSE;
}
/************************************************************************/
static void freeBonds(SList* bondsList)
{
	if(!bondsList) return;
	bondsList->klass->destroy(bondsList);
}
/************************************************************************/
static void setMultipleBonds(MoleculeData* m, int* nBonds)
{
	SList* bondsList = m->bondsList;
	Molecule* mol = m->mol;
        LNode *target;

	if(!bondsList) return;

        target = bondsList->head;
	while(target)
	{
		BondType* data=(BondType*)target->data;
		int i = data->n1;
		int j = data->n2;
		if(data->bondType == CCHEMI_BONDTYPE_HYDROGEN) continue;
		if(
		 nBonds[i] < mol->atoms[i].prop.maximumBondValence -1 &&
		 nBonds[j] < mol->atoms[j].prop.maximumBondValence -1 
		)
		{
			data->bondType = CCHEMI_BONDTYPE_TRIPLE;
			nBonds[i] += 2;
			nBonds[j] += 2;
		}
		else
		if(
		 nBonds[i] < mol->atoms[i].prop.maximumBondValence &&
		 nBonds[j] < mol->atoms[j].prop.maximumBondValence 
		)
		{
			data->bondType = CCHEMI_BONDTYPE_DOUBLE;
			nBonds[i] += 1;
			nBonds[j] += 1;
		}

		target = target->next;
	}
}
/************************************************************************/
static boolean buildBonds(MoleculeData* m, boolean buildHbonds)
{
	int i;
	int j;
	Molecule* mol = m->mol;
	int nAtoms = mol->nAtoms;
	int* nBonds = NULL;
	m->bondsList = newSList();
	if(nAtoms<1) return FALSE;

	nBonds =malloc(nAtoms*sizeof(int));
	for(i = 0;i<nAtoms;i++) nBonds[i] = 0;
	for(i = 0;i<nAtoms;i++)
	{
		for(j=i+1;j<nAtoms;j++)
			if(bonded(mol, i,j))
			{
				BondType* A=malloc(sizeof(BondType));
				A->n1 = i;
				A->n2 = j;
				nBonds[i]++;
				nBonds[j]++;
				A->bondType = CCHEMI_BONDTYPE_SINGLE;
				m->bondsList->klass->append(m->bondsList,A);
			}
		        else
			if(buildHbonds && hbonded(mol, nAtoms, i,j))
			{
				BondType* A=malloc(sizeof(BondType));
				A->n1 = i;
				A->n2 = j;
				A->bondType = CCHEMI_BONDTYPE_HYDROGEN;
				m->bondsList->klass->append(m->bondsList,A);
			}
	  }
	setMultipleBonds(m, nBonds);
	free(nBonds);
	nBonds = NULL;
	return TRUE;
}
/************************************************************************/
static void freeStack(MoleculeData* m)
{
	if(!m) return;
	if(m->inStack) free(m->inStack);
	m->inStack = NULL;
}
/************************************************************************/
static void initStack(MoleculeData* m)
{
	int i;
	if(!m) return;
	freeStack(m);
	m->inStack = malloc(m->mol->nAtoms*sizeof(boolean));
	for(i=0;i<m->mol->nAtoms;i++) m->inStack[i] = FALSE;
}
/************************************************************************/
static void freeConnections(MoleculeData* m)
{
	int i;
	if(!m->connected) return;
	for(i=0;i< m->mol->nAtoms;i++)
	{
		if(m->connected[i]) free(m->connected[i]);
	}
	free(m->connected);
	m->connected = NULL;
}
/************************************************************************/
/*
static void printConnections(MoleculeData* m)
{
	int nAtoms = m->nAtoms;
	int i = 0;
	int j = 0;
	int** connected = m->connected;
	for(i=0;i<nAtoms;i++)
	{
		printf("Nc = %d : ", i+1);
		for(j=0;j<connected[i][0];j++)
			printf(" %d ",connected[i][j+1]+1);
		printf("\n");
	}
}
*/
/************************************************************************/
static void buildConnections(MoleculeData* m)
{
	int nAtoms = m->mol->nAtoms;
	SList* bondsList = m->bondsList;
	int i = 0;
	int k = 0;
	int** connected = NULL;
        LNode *target;
	if(!bondsList) return;
	connected = malloc(nAtoms*sizeof(int*));
	for(i=0;i<nAtoms;i++)
	{
		connected[i] = malloc((nAtoms+1)*sizeof(int));
		connected[i][0] = 0;
	}
        target = bondsList->head;
	while(target)
	{
		BondType* data=(BondType*)target->data;
		int i = data->n1;
		int j = data->n2;

		connected[i][0]++;
		connected[j][0]++;

		k = connected[i][0];
		connected[i][k]=j;

		k = connected[j][0];
		connected[j][k]=i;
		target = target->next;
	}
	m->connected = connected;
}
/************************************************************************/
static boolean inRingRecursive(MoleculeData* m, int currentAtom, int rootAtom, int ringSize, boolean begin)
{
	int i;
	int** connected = m->connected;

	if(!begin) m->inStack[currentAtom] = TRUE;

	if (m->doneRing) return TRUE;
	else if ( ( currentAtom == rootAtom ) && ( m->nBondsRing == ringSize ) ) return TRUE;
	else if ( ( currentAtom == rootAtom ) && ( m->nBondsRing > 2 ) && ( ringSize < 3 ) ) return TRUE;
	if ( m->nBondsRing < ringSize )
	{
		int numberOfConnections = connected[ currentAtom ][ 0 ];
		for (i = 1; i <= numberOfConnections; i++ )
		{
			int newAtom = connected[currentAtom][i];
			if ( ! ( m->inStack[newAtom] ) )
			{
				m->nBondsRing++;
				m->doneRing = inRingRecursive( m, newAtom, rootAtom, ringSize, FALSE);
			}
			if (m->doneRing) return TRUE;
		}
	}
	m->inStack[currentAtom] = FALSE;
	m->nBondsRing--;
	return FALSE;
}
/************************************************************************/
static boolean atomInRing(MoleculeData* m, int numAtom, int ringSize)
{
	m->doneRing = FALSE;
	m->nBondsRing = 0;
	initStack(m);
	if(!m->connected) return FALSE;
	return inRingRecursive( m, numAtom,  numAtom, ringSize, TRUE);
}
/************************************************************************/
static boolean allCarbon(MoleculeData* m)
{
	Molecule* mol = m->mol;
	int i;
	for(i=0;i<mol->nAtoms;i++)
		if(m->inStack[i] && strcmp(mol->atoms[i].prop.symbol,"C" )) return FALSE;
	return TRUE;
}
/************************************************************************/
static char* subString(char* str, int begin, int end)
{
	char* res = NULL;
	int i;
	int l = strlen(str);
	int l2 = 0;
	if(l<begin) return NULL;
	if(l<end) end = l;
	if(end<0) end = l;
	l2 = end - begin + 1;
	res= malloc((l2+1)*sizeof(char));
	for(i=0;i<l2;i++) res[i] = str[i+begin];
	res[l2] = '\0';
	return res;
}
/************************************************************************/
static char getCharFromString(char* str, int index)
{
	int l = 0;
	if(!str) return '\0';
	l = strlen(str);
	if(l<=index) return '\0';
	return str[index];
}
/************************************************************************/
static boolean bondedType(MoleculeData* m, int atom1, int atom2, CChemIBondType type)
{
	SList* bondsList= m->bondsList;
        LNode *target;
	if(!bondsList) return FALSE;
        target = bondsList->head;
	while(target)
	{
		BondType* data=(BondType*)target->data;
		int i = data->n1;
		int j = data->n2;
		if( (i== atom1 && j == atom2 ) || (j== atom1 && i == atom2 ) )
		{
			if(data->bondType == type ) return TRUE;
		}
		target = target->next;
	}
	return FALSE;
}
/************************************************************************/
static int numberOfConnectedTypes( MoleculeData* m, int atomA, char* name )
{
	Molecule* mol = m->mol;
	int numberOfConnections = m->connected[ atomA ][ 0 ];
	int numberMatching = 0;
	int i;
	if(!name) return 0;
	if (!strcmp(name,"*")) return numberOfConnections;
	for ( i = 1; i <= numberOfConnections; i++ )
	{
		int atomB = m->connected[ atomA ][ i ];
		if ( !strcmp(mol->atoms[atomB].prop.symbol,name ))
			numberMatching++;
	}
	return( numberMatching );
}
/************************************************************************/
int sp(MoleculeData*m, int atomNumber)
{
	int i;
	int spHybridation = 0;
	int totalNumberOfBonds = 0;


	if(!m->connected) return 0;
	totalNumberOfBonds = m->connected[ atomNumber ][ 0 ];
	if ( m->numberOfTripleBonds[atomNumber] >= 1 ) return 1;
	else if ( m->numberOfDoubleBonds[atomNumber] >= 2 ) return 1;
	else if ( m->numberOfDoubleBonds[atomNumber] == 1 ) spHybridation += 2;
	for (  i = 1; i < totalNumberOfBonds; i++ )
	{
		int bondedToAtomNumber = m->connected[ atomNumber ][ i ];
		if (m->numberOfTripleBonds[bondedToAtomNumber] >= 1 ) return 1;
		else if ( bondedType( m, bondedToAtomNumber, atomNumber, CCHEMI_BONDTYPE_DOUBLE )) spHybridation += 2;
	}
	if ( ( spHybridation == 2 ) && ( totalNumberOfBonds == 3 ) ) return 2;
	else if ( ( spHybridation == 4 ) && ( totalNumberOfBonds == 2 ) ) return 1;
	return ( totalNumberOfBonds - 1 );
}
/************************************************************************/
static boolean isConnectedTo( MoleculeData* m, int atomA, char* expression, boolean initialize )
{
	static boolean* inStack = NULL;
	int i;
	int numberOfConnections;
	char* rootAtom;
	char* restOfExpression;
	char firstChar;
	char lastChar;
	int openClose = 0;
	int begin = 0;
	int multiplicity = 1;
	int bondOrder = 0;
	char* index;
	SList* stack = NULL;
	Molecule* mol = m->mol;

	if ( initialize )
	{
		inStack = realloc(inStack, mol->nAtoms*sizeof(boolean));	
		for( i = 0; i < mol->nAtoms; i++ )
				inStack[ i ] = FALSE;
	}
	inStack[ atomA ] = TRUE;
	if ( expression == NULL ) return TRUE;
	if(!m->connected) return FALSE;
	numberOfConnections = m->connected[ atomA ][ 0 ];
	/* printf("numberOfConnections = %d\n",numberOfConnections);*/

	rootAtom = strdup(expression);
	restOfExpression = NULL;
	index = strstr(expression, "(" );
	if ( index != expression )
	{  /* find the new root atom*/
		if ( index != NULL) 
		{
			int n0 = index - expression;
			int n1 = strlen(expression)-1;
			rootAtom = subString(expression, 0, n0-1); /* "N", e.g. */
			restOfExpression = subString(expression, n0, n1); /*  "(C)(C)", e.g. */
		}
		else
			restOfExpression = NULL;
		/*
		if(!strcmp(mol->atoms[atomA].prop.symbol,"N"))
		{
			if(expression) printf("expression = %s\n", expression);
			printf("rootAtom = %s\n", rootAtom);
			if(restOfExpression) printf("restOfExpression = %s\n", restOfExpression);
			else printf("restOfExpression = NULL\n");
		}
		*/

		firstChar = getCharFromString(rootAtom, 0);
		if(rootAtom) lastChar = getCharFromString(rootAtom, strlen(rootAtom) - 1);
		else lastChar = '\0';
		if ( firstChar != '\0' )
		{
			char* rA = rootAtom;
			if ( firstChar == '-' )
			{
				bondOrder = 1;
				rootAtom = subString(rootAtom, 1, -1);
				if(rA) free(rA);
			}
			else if ( firstChar == '=' )
			{
				bondOrder = 2;
				rootAtom = subString(rootAtom, 1, -1);
				if(rA) free(rA);
			}
			else if ( firstChar == '~')
			{
				bondOrder = 2;
				rootAtom = subString(rootAtom, 1, -1);
				if(rA) free(rA);
			}
			else if ( firstChar == '#' )
			{
				bondOrder = 3;
				rootAtom = subString(rootAtom, 1, -1);
				if(rA) free(rA);
			}
		}
		if ( lastChar  != '\0')
		{
			char* rA = rootAtom;
			if ( lastChar == '2' )
			{
				multiplicity = 2;
				rootAtom = subString(rootAtom, 0, strlen(rootAtom) - 2);
				if(rA) free(rA);
			}
			else if ( lastChar == '3' )
			{
				multiplicity = 3;
				rootAtom = subString(rootAtom, 0, strlen(rootAtom) - 2);
				if(rA) free(rA);
			}
			else if ( lastChar == '4'  )
			{
				multiplicity = 4;
				rootAtom = subString(rootAtom, 0, strlen(rootAtom) - 2);
				if(rA) free(rA);
			}
		}
		if ( multiplicity > numberOfConnectedTypes( m, atomA, rootAtom ) ) return FALSE;
		for (i = 1; i <= numberOfConnections; i++ )
		{
			int atomB = m->connected[ atomA ][ i ];
			if ( inStack[ atomB ] ) continue;
			if ( !strcmp(mol->atoms[atomB ].prop.symbol,rootAtom) || ( !strcmp(rootAtom,"*" ) ) )
			{
				if ( isConnectedTo( m, atomB, restOfExpression, FALSE ) )
				{
					if ( bondOrder == 0 ) return TRUE;
					else if ( bondOrder == 3 )
					{
						if ( ( bondedType( m, atomA, atomB, CCHEMI_BONDTYPE_TRIPLE ) ) || 
						      ( sp(m,  atomA ) == 1 ) || ( sp(m, atomB ) == 1 ) )
								return TRUE;
					}
					else if ( bondOrder == 2 )
					{
						if ( ( bondedType( m, atomA, atomB, CCHEMI_BONDTYPE_DOUBLE ) ) || 
						     ( sp(m, atomA ) == 2 ) || ( sp(m, atomB ) == 2 ) )
								return TRUE;
					}
					else if ( bondOrder == 1 )
					{
						if ( sp(m, atomA ) + sp(m, atomB ) >= 4 ) return TRUE;
					}
				}
			}
		}
		return FALSE;
	} 
	/* below, push (C)(C) onto stack */
	restOfExpression = strdup(expression);
	stack = newSList();
	openClose = 0;
	for (i = 0; i < strlen(restOfExpression); i++ )
	{
		if ( restOfExpression[i] == '(' )
		{
			if ( openClose++ == 0 ) begin = i + 1;
		}
		else if ( restOfExpression[i] == ')' )
		{
			if ( --openClose == 0 ) stack->klass->append(stack,subString(expression, begin, i-1) );
		}
	}
	while ( stack  && stack->size > 0 )
	{
		char* newExpression = NULL;
		LNode* last =  stack->klass->last(stack);
		if(!last || !last->data)
		{
			fprintf(stderr,"Error problem with last stack, stakc size = %d\n",stack->size);
			exit(1);
		}
		newExpression = (char*)(last->data);
		stack->klass->remove(stack, newExpression);
		if ( ! ( isConnectedTo( m, atomA, newExpression, FALSE ) ) ) 
		{
			stack->klass->destroy(stack);
			return FALSE;
		}
		// ? free(newExpression);
	}
	return TRUE;
}
/**********************************************************************************************************/
static char* getAmberTypeOfAtom(MoleculeData* m, int atomNumber)
{
	Molecule* mol = m->mol;
	/* printf("Atom number = %d symbol = %s\n",atomNumber, mol->atoms[atomNumber].prop.symbol);*/
	if ( !strcmp(mol->atoms[atomNumber].prop.symbol,"H" ))
	{
		if (  isConnectedTo( m, atomNumber, "N(*4)", TRUE ) ) return "H";
		else if (  isConnectedTo( m, atomNumber, "N(H)(C(N(H2)))", TRUE ) ) return "H";
		else if ( (  isConnectedTo( m, atomNumber, "N(C(N2))", TRUE ) ) &&
			 isConnectedTo( m, atomNumber, "N(C(N(H2))(N(H2)))", TRUE ) ) return "H";
		else if (  isConnectedTo( m, atomNumber, "C", TRUE ) ) 
		{
			int numberOfConnections = m->connected[ atomNumber ][ 0 ];
			int i;
			for (i = 1; i <= numberOfConnections; i++ )
			{
				int newAtom = m->connected[atomNumber][i];
				if(!strcmp(mol->atoms[newAtom].prop.symbol,"C" ))
				{
					if(atomInRing( m, newAtom, 6 ) && allCarbon(m)) return "HA";
					else if (  isConnectedTo( m, newAtom, "N2", TRUE ) ) return "H5";
					else if (  isConnectedTo( m, newAtom, "N(C)", TRUE ) ) return "H4";
					else if (  isConnectedTo( m, newAtom, "N(H3)", TRUE ) ) return "HP";
					else if (  isConnectedTo( m, newAtom, "N(*)", TRUE ) ) return "H1";
					else if (  isConnectedTo( m, newAtom, "S", TRUE ) ) return "H1";
					else if (  isConnectedTo( m, newAtom, "O", TRUE ) ) return "H1";
				}
			}
			return "HC";
		}
		else if (  isConnectedTo( m, atomNumber, "O(H)", TRUE ) ) return "HW";
		else if (  isConnectedTo( m, atomNumber, "O", TRUE ) ) return "HO";
		else if (  isConnectedTo( m, atomNumber, "S", TRUE ) ) return "HS";
		/*
		else if (  isConnectedTo( m, atomNumber, "N(H)(C(~*))", TRUE ) ) return "H2";
		else if (  isConnectedTo( m, atomNumber, "N(H)(C(=N)(-C))", TRUE ) ) return "H2";
		else if (  isConnectedTo( m, atomNumber, "N(H)(C(=N)(-N))", TRUE ) ) return "H2";
		else if (  isConnectedTo( m, atomNumber, "N(H2)", TRUE ) ) return "H2";
		*/
		else return "H";
	}
	else if ( !strcmp(mol->atoms[atomNumber].prop.symbol,"C" ))
	{
		/* printf("sp = %d\n",sp(m, atomNumber ));*/
		if ( sp(m, atomNumber ) == 2 )
		{
			if ( ( atomInRing( m, atomNumber, 5 ) ) && ( atomInRing( m, atomNumber, 6 ) ) )
			{
				if (  isConnectedTo( m, atomNumber, "C(N)", TRUE ) ) return "CB";
				else if (  isConnectedTo( m, atomNumber, "N(Fe)", TRUE ) ) return "CC";
				else if (  isConnectedTo( m, atomNumber, "N", TRUE ) ) return "CN";
				else return "CT";
			}
			else if ( atomInRing( m, atomNumber, 5 ) )
			{
				if (  isConnectedTo( m, atomNumber, "(C)(C)(C(N(Fe)))", TRUE ) ) return "CB";
				else if ( (  isConnectedTo( m, atomNumber, "(H)(N2)", TRUE ) )
					&& (  isConnectedTo( m, atomNumber, "(N)(H)(N(C(N(C(N(C))))))", TRUE ) ) )
							return "CK";
				else if (  isConnectedTo( m, atomNumber,"(N2)(H)", TRUE ) ) return "CR";
				else if (  isConnectedTo( m, atomNumber,"(C(H2))(N)", TRUE ) ) return "CC";
				else if (  isConnectedTo( m, atomNumber,"C(H2)", TRUE ) ) return "C*";
				else if (  isConnectedTo( m, atomNumber,"N(~*)(*)", TRUE ) ) return "CW";
				else if (  isConnectedTo( m, atomNumber,"N", TRUE ) ) return "CV";
				else return "CT";
			}
			else if (atomInRing( m, atomNumber, 6 ) )
			{
				if (  isConnectedTo( m, atomNumber, "O(H)", TRUE ) ) return "C";
				else if (  isConnectedTo( m, atomNumber, "=O", TRUE ) ) return "C";
				else if (  isConnectedTo( m, atomNumber, "(N2)(H)", TRUE ) ) return "CQ";
				else if (  isConnectedTo( m, atomNumber, "N(H2)", TRUE ) ) return "CA";
				else if (  isConnectedTo( m, atomNumber, "C(N(C(=O)))", TRUE ) ) return "CM";
				else if (  isConnectedTo( m, atomNumber, "N(C(=O))", TRUE ) ) return "CM";
				else return "CA";
			}
			else if (  isConnectedTo( m, atomNumber, "(~O2)", TRUE ) ) return "C";
			else if (  isConnectedTo( m, atomNumber, "=O", TRUE ) ) return "C";
			else if (  isConnectedTo( m, atomNumber, "(=C)(-S)", TRUE ) ) return "CY";
			/*else if (  isConnectedTo( m, atomNumber, "(=C)", TRUE ) ) return "CX";*/
			else if (  isConnectedTo( m, atomNumber, "(N3)", TRUE ) ) return "CA";
			else return "CT";
		}
		else if ( sp(m, atomNumber ) == 3 )
		{
			if (  isConnectedTo( m, atomNumber, "(N3)", TRUE ) ) return "CA";
			else if (  isConnectedTo( m, atomNumber, "=O", TRUE ) ) return "C";
			else if (  isConnectedTo( m, atomNumber, "(=C)(-S)", TRUE ) ) return "CY";
			/* else if (  isConnectedTo( m, atomNumber, "(=C)", TRUE ) ) return "CT";*/
			else return "CT";
		}
		else { 
			if (  isConnectedTo( m, atomNumber, "(N3)", TRUE ) ) return "CA";
			else if (  isConnectedTo( m, atomNumber, "=O", TRUE ) ) return "C";
			else if (  isConnectedTo( m, atomNumber, "(=C)(-S)", TRUE ) ) return "CY";
			/* else if (  isConnectedTo( m, atomNumber, "(=C)", TRUE ) ) return "CX";*/
			else if (  isConnectedTo( m, atomNumber, "#C", TRUE ) ) return "CZ";
			else return "CT";
		}
	}
	else if ( !strcmp(mol->atoms[atomNumber].prop.symbol,"N" ) )
	{
		if ( sp(m, atomNumber ) < 3 )
		{
			if ( atomInRing( m, atomNumber, 5 ) )
			{
				if ( (  isConnectedTo( m, atomNumber, "(C3)", TRUE ) ) && 
				     (  isConnectedTo( m, atomNumber, "(C)(C(=*))(C(=*))", TRUE ) ) ) return "N*";
				else if (  isConnectedTo( m, atomNumber, "(C3)", TRUE ) ) return "N*";
				else if (  isConnectedTo( m, atomNumber, "Fe", TRUE ) ) return "NA";
				else if (  isConnectedTo( m, atomNumber, "H", TRUE ) ) return "NA";
				else return "NB";
			}
			else if (atomInRing( m, atomNumber, 6 ) )
			{
				if (  isConnectedTo( m, atomNumber, "H", TRUE ) ) return "NA";
				else if ( (  isConnectedTo( m, atomNumber, "(H)(C2)", TRUE ) ) 
					&& (  isConnectedTo( m, atomNumber, "(H)(C(=O))(C(=O))", TRUE ) ) ) return "NA";
				else if (  isConnectedTo( m, atomNumber, "(H)(C(=O))(C(=N))", TRUE ) ) return "NA";
				else if (  isConnectedTo( m, atomNumber, "(C3)", TRUE ) ) return "N*";
				else return "NC";
			}
			else if (  isConnectedTo( m, atomNumber, "-C(=O)", TRUE ) ) 
			{
				return "N";
			}
			else if (  isConnectedTo( m, atomNumber, "-C(-C(=O))", TRUE ) )
			{
				return "N";
			}
			else if (  isConnectedTo( m, atomNumber, "C(N2)", TRUE ) ) return "N2";
			else if (  isConnectedTo( m, atomNumber, "(H2)(C~*)", TRUE ) ) return "N2";
			else if (  isConnectedTo( m, atomNumber, "(H2)(C=*)", TRUE ) ) return "N2";
			else return "NT";
		}
		else if ( sp(m, atomNumber ) == 3 ) return "N3";
		else return "NT";
	}
	else if ( !strcmp(mol->atoms[atomNumber].prop.symbol,"O" ))
	{
		if (  isConnectedTo( m, atomNumber, "~C(~O2)", TRUE ) ) 
		{
			if (  isConnectedTo( m, atomNumber, "H", TRUE ) ) return "OH";
			else return "O2";
		}
		else if (  isConnectedTo( m, atomNumber, "H2", TRUE ) ) return "OW";
		else if (  isConnectedTo( m, atomNumber, "H", TRUE ) ) return "OH";
		else if (  isConnectedTo( m, atomNumber, "=C", TRUE ) ) return "O";
		else if (  isConnectedTo( m, atomNumber, "(H)(C(=O))", TRUE ) ) return "OH";
		else if (  isConnectedTo( m, atomNumber, "-C2", TRUE ) ) return "OS";
		else if (  isConnectedTo( m, atomNumber, "(C(=O))(C)", TRUE ) ) return "OS";
		else if (  isConnectedTo( m, atomNumber, "(C)(P)", TRUE ) ) return "OS";
		else if (  isConnectedTo( m, atomNumber, "P", TRUE ) ) return "O2";
		else return "OS";
	}
	else if ( !strcmp(mol->atoms[atomNumber].prop.symbol,"F" ) ) return "F";
	else if ( !strcmp(mol->atoms[atomNumber].prop.symbol,"Na" ) ) return "IP";
	else if ( !strcmp(mol->atoms[atomNumber].prop.symbol,"Mg" )) return "MG";
	else if ( !strcmp(mol->atoms[atomNumber].prop.symbol,"P" ) ) return "P";
	else if ( !strcmp(mol->atoms[atomNumber].prop.symbol,"S" ) )
	{
		if (  isConnectedTo( m, atomNumber, "H", TRUE ) ) return "SH";
		else return "S";
	}
	else if ( !strcmp(mol->atoms[atomNumber].prop.symbol,"Cl" ))
	{
		if ( !m->connected ) return "CL";
		if ( m->connected[ atomNumber ][ 0 ] > 0 ) return "CL";
		else return "IM";
	}
	else if ( !strcmp(mol->atoms[atomNumber].prop.symbol,"Ca" ) ) return "C0";
	else if ( !strcmp(mol->atoms[atomNumber].prop.symbol,"Fe" ) ) return "FE";
	else if ( !strcmp(mol->atoms[atomNumber].prop.symbol,"Cu" ) ) return "CU";
	else if ( !strcmp(mol->atoms[atomNumber].prop.symbol,"Br" ) ) return "BR";
	else if ( !strcmp(mol->atoms[atomNumber].prop.symbol,"I" ) ) return "I";

	return mol->atoms[atomNumber].prop.symbol;
}
/************************************************************************/
static int* getNumberOfBonds(MoleculeData* m, CChemIBondType type)
{
	int i;
	int j;
	SList* bondsList = m->bondsList;
	int nAtoms = m->mol->nAtoms;
	int* numberOfBonds = NULL;
        LNode *target;

	if(nAtoms<1) return NULL;
	numberOfBonds =malloc(nAtoms*sizeof(int));
	for(i = 0;i<nAtoms;i++) numberOfBonds[i] = 0;

        target = bondsList->head;
	while(target)
	{
		BondType* data=(BondType*)target->data;
		i = data->n1;
		j = data->n2;
		if(data->bondType == type) 
		{
			numberOfBonds[i]++;
			numberOfBonds[j]++;
		}
		target = target->next;
	}
	return numberOfBonds;
}
/**********************************************************************************************************/
static MoleculeData getMyMoleculeData(Molecule* mol)
{
	MoleculeData m;
	m.mol = mol;

	m.bondsList = NULL;
	m.connected = NULL;
	m.numberOfDoubleBonds = NULL;
	m.numberOfTripleBonds = NULL;
	m.inStack = NULL;
	m.doneRing = FALSE;
	m.nBondsRing = 0;
	buildBonds(&m, FALSE);
	m.numberOfDoubleBonds = getNumberOfBonds(&m, CCHEMI_BONDTYPE_DOUBLE);
	m.numberOfTripleBonds = getNumberOfBonds(&m, CCHEMI_BONDTYPE_TRIPLE);
	buildConnections(&m);
	/* printConnections(&m);*/
	return m;
}
/**********************************************************************************************************/
static void freeMoleculeData(MoleculeData* m)
{
	free(m->numberOfDoubleBonds);
	free(m->numberOfTripleBonds);
	freeConnections(m);
	freeBonds(m->bondsList);
}
/**********************************************************************************************************/
void calculAmberTypes(Molecule* mol)
{
	int i;
	MoleculeData m;
	if(mol->nAtoms<1) return;
	m = getMyMoleculeData(mol);

	printf("End getMyMoleculeData\n");
	for(i=0;i<mol->nAtoms;i++)
	{
 
		/*
		printf("Atom n = %d symb = %s\n", i+1, mol->atoms[i].prop.symbol);
		if( atomInRing( &m, i, 5 )) printf("atom number %d is in pentagon\n",i+1);
		if( atomInRing( &m, i, 6 )) printf("atom number %d is in hexagon\n",i+1);
		if (  isConnectedTo( &m, i, "C", TRUE ) ) printf("atom number %d  have %s type \n",i+1,"HC");
		printf("mm %s\n", mol->atoms[i].mmType);
		*/

		/* printf("i = %d\n",i);*/
		if(mol->atoms[i].mmType) free(mol->atoms[i].mmType);
		mol->atoms[i].mmType = strdup(getAmberTypeOfAtom(&m, i));
	}
	freeMoleculeData(&m);
}
