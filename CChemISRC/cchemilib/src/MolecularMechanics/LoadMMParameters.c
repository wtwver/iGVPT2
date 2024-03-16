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

/* LoadMMParameters.c */
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>

#include "../Utils/Utils.h"
#include "../Utils/AtomsProp.h"
#include "../Utils/Constants.h"
#include "../Utils/Types.h"
#include "../Molecule/Molecule.h"
#include "../MolecularMechanics/MolecularMechanics.h"

static char atomTypesTitle[]       = "Begin  INPUT FOR ATOM TYPES, MASSE AND POLARISABILITIES";
static char bondStretchTitle[]     = "Begin INPUT FOR BOND LENGTH PARAMETERS";
static char angleBendTitle[]       = "Begin INPUT FOR BOND ANGLE PARAMETERS";
static char strBendTitle[]       = "Begin INPUT FOR STRETCH-BEND PARAMETERS";
static char hydrogenBonded1012Title[]  = "Begin INPUT FOR H-BOND 10-12 POTENTIAL PARAMETERS";
static char hydrogenBondedMorseTitle[]  = "Begin INPUT FOR H-BOND MORSE POTENTIAL PARAMETERS";
static char suttonChenTitle[]      = "Begin INPUT FOR SUTTON-CHEN POTENTIAL PARAMETERS";
static char improperTorsionTitle[] ="Begin INPUT FOR IMPROPER DIHEDRAL PARAMETERS";
static char outOfPlaneTitle[]      ="Begin INPUT FOR OUT OF PLANE PARAMETERS";
static char vdw612Title[]       ="Begin INPUT FOR THE NON-BONDED 6-12 POTENTIAL PARAMETERS";
static char vdw714Title[]       ="Begin INPUT FOR THE NON-BONDED 7-14 POTENTIAL PARAMETERS";
static char dihedralAngleTitle[]   = "Begin INPUT FOR DIHEDRAL PARAMETERS";
static char pairWiseTitle[]        = "Begin INPUT FOR PAIR WISE PARAMETERS";
static char bondHardnessTitle[]    = "Begin INPUT FOR BOND HARDNESS PARAMETERS";

/**********************************************************************/
static boolean readAmberTypes(AmberParameters* amberParameters, FILE* file)
{
	char t[BSIZE];
	char dumpName[BSIZE];
	char dumpSymb[BSIZE];
	int len = BSIZE;
	boolean Ok = FALSE;
	int n = 0;
	int nt;
	AmberAtomTypes* types = NULL;
	/* Search Begin INPUT FOR  ATOM TYPES */ 

	rewind(file);
	while(!feof(file))
	{
		if(fgets(t,len,file))
		{
			if(strstr(t,atomTypesTitle))
			{
				Ok = TRUE;
				break;
			}
		}
	}
	if(!Ok)
		return FALSE;

	types = malloc(sizeof(AmberAtomTypes));
	n = 0;
	Ok = FALSE;
	while(!feof(file))
	{
		if(fgets(t,len,file))
		{
			if(strstr(t,"End") || strstr(t,"END"))
			{
				Ok = TRUE;
				break;
			}
		}
		else 
		{
			Ok = FALSE;
			break;
		}



		nt = sscanf(t,"%s %s %d %lf %lf %lf %lf %lf %lf",
			dumpName,
			dumpSymb,
			&types[n].number,
			&types[n].mass,
			&types[n].polarisability, &types[n].charge, &types[n].electronegativity, &types[n].hardness, &types[n].width);
		if(nt<9) types[n].width = 0.0;
		if(nt<8) types[n].hardness = 0.0;
		if(nt<7) types[n].electronegativity = 0.0;
		if(nt<6) types[n].charge = 0.0;
		if(nt<5) types[n].polarisability = 0.0;
		if(nt<4) types[n].mass = 1.0;

	      	types[n].name = strdup(dumpName);	
	      	types[n].symbol = strdup(dumpSymb);
		types[n].number--;

		n++;
		types = realloc(types,(n+1)*sizeof(AmberAtomTypes));
	}
	if(n==0 || !Ok )
		free(types);
	else
	{
		amberParameters->numberOfTypes = n;
		amberParameters->atomTypes = types;
	}
	/* printing for test*/
	/*
	printf("umber of types = %d \n",amberParameters->numberOfTypes);
	for(n=0;n<amberParameters->numberOfTypes;n++)
	{
		printf("%s\t %d\t",
				amberParameters->atomTypes[n].name,
				amberParameters->atomTypes[n].number
				);
	}
	printf("\n");
	*/

	return TRUE;
			

}
/**********************************************************************/
static boolean readAmberBondStretchTerms(AmberParameters* amberParameters,FILE* file)
{
	char t[BSIZE];
	int len = BSIZE;
	boolean Ok = FALSE;
	int n = 0;
	AmberBondStretchTerms* terms = NULL;
	int nt = 0;

	/* Search Begin INPUT FOR  ATOM TYPES */ 

	rewind(file);
	while(!feof(file))
	{
		if(fgets(t,len,file))
		{
			if(strstr(t,bondStretchTitle))
			{
				Ok = TRUE;
				break;
			}
		}
	}
	if(!Ok)
		return FALSE;

	terms = malloc(sizeof(AmberBondStretchTerms));
	n = 0;
	Ok = FALSE;
	while(!feof(file))
	{
		if(fgets(t,len,file))
		{
			if(strstr(t,"End") || strstr(t,"END"))
			{
				Ok = TRUE;
				break;
			}
		}
		else 
		{
			Ok = FALSE;
			break;
		}

		// add anharmonic terms
		//sscanf(t,"%d %d %lf %lf", &terms[n].numbers[0], &terms[n].numbers[1], &terms[n].forceConstant, &terms[n].equilibriumDistance);
		terms[n].h3 = 0;
		terms[n].h4 = 0;
		terms[n].h5 = 0;
		terms[n].h6 = 0;
		nt = sscanf(t,"%d %d %lf %lf %lf %lf %lf %lf", &terms[n].numbers[0], &terms[n].numbers[1], 
				&terms[n].forceConstant, &terms[n].equilibriumDistance,
				&terms[n].h3, &terms[n].h4,
				&terms[n].h5, &terms[n].h6
				);
		if(nt<8) terms[n].h6 = 0;
		if(nt<7) terms[n].h5 = 0;
		if(nt<6) terms[n].h4 = 0;
		if(nt<5) terms[n].h3 = 0;
		terms[n].type = 0;// harmonic
		uppercase(t);
		if(strstr(t,"MORSE")) terms[n].type = 1;//Morse, h3 = De in kcal/mol
		if( terms[n].type ==1 && fabs(terms[n].h3)<1e-13)
		{
			fprintf(stderr,"===================================================================\n");
			fprintf(stderr,"Error in Molecular mechanics parameters file : De = 0 for Morse type \n");
			fprintf(stderr,"Check Strech terms  for type numbers : %d %d \n",terms[n].numbers[0],terms[n].numbers[1]);
			fprintf(stderr,"===================================================================\n");
			exit(1);
		}

	      	terms[n].numbers[0]--;
	      	terms[n].numbers[1]--;
		if(terms[n].numbers[0]>terms[n].numbers[1])
		{
			int t = terms[n].numbers[0];
			terms[n].numbers[0] = terms[n].numbers[1];
			terms[n].numbers[1] = t;
		}

		n++;
		terms = realloc(terms,(n+1)*sizeof(AmberBondStretchTerms));
	}
	if(n==0 || !Ok )
		free(terms);
	else
	{
		amberParameters->numberOfStretchTerms = n;
		amberParameters->bondStretchTerms = terms;
	}
	/* printing for test*/
	/*
	printf("number of bonds = %d \n",amberParameters->numberOfStretchTerms);
	for(n=0;n<amberParameters->numberOfStretchTerms;n++)
	{
		printf("%d %d %f %f\n",
				amberParameters->bondStretchTerms[n].numbers[0],
				amberParameters->bondStretchTerms[n].numbers[1],
				amberParameters->bondStretchTerms[n].forceConstant,
				amberParameters->bondStretchTerms[n].equilibriumDistance
				);
	}
	printf("\n");
	*/

	return TRUE;
			

}
/**********************************************************************/
static boolean readAmberAngleBendTerms(AmberParameters* amberParameters,FILE* file)
{
	char t[BSIZE];
	int len = BSIZE;
	boolean Ok = FALSE;
	int n = 0;
	AmberAngleBendTerms* terms = NULL;
	int nt = 0;

	/* Search Begin INPUT FOR  ATOM TYPES */ 

	rewind(file);
	while(!feof(file))
	{
		if(fgets(t,len,file))
		{
			if(strstr(t,angleBendTitle))
			{
				Ok = TRUE;
				break;
			}
		}
	}
	if(!Ok)
		return FALSE;

	terms = malloc(sizeof(AmberAngleBendTerms));
	n = 0;
	Ok = FALSE;
	while(!feof(file))
	{
		if(fgets(t,len,file))
		{
			if(strstr(t,"End") || strstr(t,"END"))
			{
				Ok = TRUE;
				break;
			}
		}
		else 
		{
			Ok = FALSE;
			break;
		}

		// ad anharmonic
		//sscanf(t,"%d %d %d %lf %lf", &terms[n].numbers[0], &terms[n].numbers[1], &terms[n].numbers[2], &terms[n].forceConstant, &terms[n].equilibriumAngle);
                terms[n].h3 = 0;
                terms[n].h4 = 0;
                terms[n].h5 = 0;
                terms[n].h6 = 0;
                nt = sscanf(t,"%d %d %d %lf %lf %lf %lf %lf %lf", &terms[n].numbers[0], &terms[n].numbers[1], &terms[n].numbers[2],
                                &terms[n].forceConstant, &terms[n].equilibriumAngle,
                                &terms[n].h3, &terms[n].h4,
                                &terms[n].h5, &terms[n].h6
                                );
		if(nt<9) terms[n].h6 = 0;
		if(nt<8) terms[n].h5 = 0;
		if(nt<7) terms[n].h4 = 0;
		if(nt<6) terms[n].h3 = 0;

		terms[n].numbers[0]--;
		terms[n].numbers[1]--;
		terms[n].numbers[2]--;
		if(terms[n].numbers[0]>terms[n].numbers[2])
		{
			int t = terms[n].numbers[0];
			terms[n].numbers[0] = terms[n].numbers[2];
			terms[n].numbers[2] = t;
		}
		

		n++;
		terms = realloc(terms,(n+1)*sizeof(AmberAngleBendTerms));
	}
	if(n==0 || !Ok )
		free(terms);
	else
	{
		amberParameters->numberOfBendTerms = n;
		amberParameters->angleBendTerms = terms;
	}
	/* printing for test*/
	/*
	printf("number of bonds = %d \n",amberParameters->numberOfBendTerms);
	for(n=0;n<amberParameters->numberOfBendTerms;n++)
	{
		printf("%d %d %d %f %f\n",
				amberParameters->angleBendTerms[n].numbers[0],
				amberParameters->angleBendTerms[n].numbers[1],
				amberParameters->angleBendTerms[n].numbers[2],
				amberParameters->angleBendTerms[n].forceConstant,
				amberParameters->angleBendTerms[n].equilibriumAngle
				);
	}
	printf("\n");
	*/

	return TRUE;
			

}
/**********************************************************************/
static boolean readAmberStrBendTerms(AmberParameters* amberParameters,FILE* file)
{
	char t[BSIZE];
	int len = BSIZE;
	boolean Ok = FALSE;
	int n = 0;
	AmberStrBendTerms* terms = NULL;
	int nt = 0;

	/* Begin INPUT FOR STRETCH-BEND PARAMETER */ 

	rewind(file);
	while(!feof(file))
	{
		if(fgets(t,len,file))
		{
			if(strstr(t,strBendTitle))
			{
				Ok = TRUE;
				break;
			}
		}
	}
	if(!Ok)
		return FALSE;

	terms = malloc(sizeof(AmberStrBendTerms));
	n = 0;
	Ok = FALSE;
	while(!feof(file))
	{
		if(fgets(t,len,file))
		{
			if(strstr(t,"End") || strstr(t,"END"))
			{
				Ok = TRUE;
				break;
			}
		}
		else 
		{
			Ok = FALSE;
			break;
		}

                nt = sscanf(t,"%d %d %d %lf %lf", &terms[n].numbers[0], &terms[n].numbers[1], &terms[n].numbers[2],
                                &terms[n].forceConstant12, &terms[n].forceConstant23
                                );
		if(nt<5) terms[n].forceConstant23 = terms[n].forceConstant12;

		terms[n].numbers[0]--;
		terms[n].numbers[1]--;
		terms[n].numbers[2]--;
		if(terms[n].numbers[0]>terms[n].numbers[2])
		{
			int t = terms[n].numbers[0];
			terms[n].numbers[0] = terms[n].numbers[2];
			terms[n].numbers[2] = t;
		}
		

		n++;
		terms = realloc(terms,(n+1)*sizeof(AmberStrBendTerms));
	}
	if(n==0 || !Ok )
		free(terms);
	else
	{
		amberParameters->numberOfStrBendTerms = n;
		amberParameters->strBendTerms = terms;
	}
	/* printing for test*/
	
	/*
	printf("number of bonds = %d \n",amberParameters->numberOfStrBendTerms);
	for(n=0;n<amberParameters->numberOfStrBendTerms;n++)
	{
		printf("%d %d %d %f %f\n",
				amberParameters->strBendTerms[n].numbers[0],
				amberParameters->strBendTerms[n].numbers[1],
				amberParameters->strBendTerms[n].numbers[2],
				amberParameters->strBendTerms[n].forceConstant12,
				amberParameters->strBendTerms[n].forceConstant23
				);
	}
	printf("\n");
	*/
	
	return TRUE;
}
/**********************************************************************/
static boolean readAmberDihedralAngleTerms(AmberParameters* amberParameters,FILE* file)
{
	char t[BSIZE];
	int len = BSIZE;
	boolean Ok = FALSE;
	int n = 0;
	AmberDihedralAngleTerms* terms = NULL;
	double divisor = 1;
	double barrier = 0;
	double phase = 0;
	double nN = 0;
	int d;

	/* Search Begin INPUT FOR  DIHEDRAL PARAMETERS */

	rewind(file);
	while(!feof(file))
	{
		if(fgets(t,len,file))
		{
			if(strstr(t,dihedralAngleTitle))
			{
				Ok = TRUE;
				break;
			}
		}
	}
	if(!Ok)
		return FALSE;

	terms = malloc(sizeof(AmberDihedralAngleTerms));

	n = 0;
	Ok = FALSE;
	while(!feof(file))
	{
		if(fgets(t,len,file))
		{
			if(strstr(t,"End") || strstr(t,"END"))
			{
				Ok = TRUE;
				break;
			}
		}
		else 
		{
			Ok = FALSE;
			break;
		}

		terms[n].nSomme = 1;
		terms[n].divisor = malloc(sizeof(double));
		terms[n].barrier = malloc(sizeof(double));
		terms[n].phase   = malloc(sizeof(double));
		terms[n].n       = malloc(sizeof(double));

		sscanf(t,"%d %d %d %d %lf %lf %lf %lf",
			&terms[n].numbers[0],
			&terms[n].numbers[1],
			&terms[n].numbers[2],
			&terms[n].numbers[3],
			&divisor,
			&barrier,
			&phase,
			&nN);

		terms[n].divisor[0] = divisor;
		terms[n].barrier[0] = barrier;
		terms[n].phase[0]   = phase;
		terms[n].n[0]       = fabs(nN);

	    terms[n].numbers[0]--;
	    terms[n].numbers[1]--;
	    terms[n].numbers[2]--;
	    terms[n].numbers[3]--;

	    if(terms[n].numbers[0]>terms[n].numbers[3])
	    {
		int t = terms[n].numbers[0];
		terms[n].numbers[0] = terms[n].numbers[3];
		terms[n].numbers[3] = t;

		t =  terms[n].numbers[1];
		terms[n].numbers[1] = terms[n].numbers[2];
		terms[n].numbers[2] = t;
	     }

		Ok = TRUE;
		while(!feof(file) && nN<0)
		{
			if(!fgets(t,len,file))
			{
				Ok = FALSE;
				break;
			}

			terms[n].nSomme++;
			terms[n].divisor = realloc(terms[n].divisor,terms[n].nSomme*sizeof(double));
			terms[n].barrier = realloc(terms[n].barrier,terms[n].nSomme*sizeof(double));
			terms[n].phase   = realloc(terms[n].phase,terms[n].nSomme*sizeof(double));
			terms[n].n       = realloc(terms[n].n,terms[n].nSomme*sizeof(double));

			sscanf(t,"%d %d %d %d %lf %lf %lf %lf",
					  &d,&d,&d,&d,
					  &divisor,&barrier,&phase,&nN);

			terms[n].divisor[terms[n].nSomme-1] = divisor;
			terms[n].barrier[terms[n].nSomme-1] = barrier;
			terms[n].phase[terms[n].nSomme-1]   = phase;
			terms[n].n[terms[n].nSomme-1]       = fabs(nN);
		}
		if(!Ok)
			break;
		n++;
		terms = realloc(terms,(n+1)*sizeof(AmberDihedralAngleTerms));
	}
	if(n==0 || !Ok )
		free(terms);
	else
	{
		amberParameters->numberOfDihedralTerms = n;
		amberParameters->dihedralAngleTerms = terms;
	}
	/* printing for test*/
	/*	
	printf("number of dihedral torsion terms = %d \n",
			amberParameters->numberOfDihedralTerms);

	for(n=0;n<amberParameters->numberOfDihedralTerms;n++)
	{
		int j;
		printf("%d %d %d %d \t",
				amberParameters->dihedralAngleTerms[n].numbers[0],
				amberParameters->dihedralAngleTerms[n].numbers[1],
				amberParameters->dihedralAngleTerms[n].numbers[2],
				amberParameters->dihedralAngleTerms[n].numbers[3]
			);
		for(j=0;j<amberParameters->dihedralAngleTerms[n].nSomme;j++)
		{
			printf("%f %f %f %f\t",
				amberParameters->dihedralAngleTerms[n].divisor[j],
				amberParameters->dihedralAngleTerms[n].barrier[j],
				amberParameters->dihedralAngleTerms[n].phase[j],
				amberParameters->dihedralAngleTerms[n].n[j]
				);
		}
		printf("\n");
	}
	printf("\n");
	*/	

	return TRUE;
			

}
/**********************************************************************/
static boolean readAmberImproperTorsionTerms(AmberParameters* amberParameters,FILE* file)
{
	char t[BSIZE];
	int len = BSIZE;
	boolean Ok = FALSE;
	int n = 0;
	AmberImproperTorsionTerms* terms = NULL;

	/* Search Begin INPUT FOR  ATOM TYPES */ 

	rewind(file);
	while(!feof(file))
	{
		if(fgets(t,len,file))
		{
			if(strstr(t,improperTorsionTitle))
			{
				Ok = TRUE;
				break;
			}
		}
	}
	if(!Ok)
		return FALSE;

	terms = malloc(sizeof(AmberImproperTorsionTerms));
	n = 0;
	Ok = FALSE;
	while(!feof(file))
	{
		if(fgets(t,len,file))
		{
			if(strstr(t,"End") || strstr(t,"END"))
			{
				Ok = TRUE;
				break;
			}
		}
		else 
		{
			Ok = FALSE;
			break;
		}



		sscanf(t,"%d %d %d %d %lf %lf %lf",
			&terms[n].numbers[0],
			&terms[n].numbers[1],
			&terms[n].numbers[2],
			&terms[n].numbers[3],
			&terms[n].barrier,
			&terms[n].phase,
			&terms[n].n);

	    terms[n].numbers[0]--;
	    terms[n].numbers[1]--;
	    terms[n].numbers[2]--;
	    terms[n].numbers[3]--;
	    if(terms[n].numbers[0]>terms[n].numbers[3])
	    {
		int t = terms[n].numbers[0];
		terms[n].numbers[0] = terms[n].numbers[3];
		terms[n].numbers[3] = t;

		t =  terms[n].numbers[1];
		terms[n].numbers[1] = terms[n].numbers[2];
		terms[n].numbers[2] = t;
	    }

		n++;
		terms = realloc(terms,(n+1)*sizeof(AmberImproperTorsionTerms));
	}
	if(n==0 || !Ok )
		free(terms);
	else
	{
		amberParameters->numberOfImproperTorsionTerms = n;
		amberParameters->improperTorsionTerms = terms;
	}
	/* printing for test*/
	/*
	printf("number of improper torsion terms = %d \n",
			amberParameters->numberOfImproperTorsionTerms);

	for(n=0;n<amberParameters->numberOfImproperTorsionTerms;n++)
	{
		printf("%d %d %d %d %f %f %f\n",
				amberParameters->improperTorsionTerms[n].numbers[0],
				amberParameters->improperTorsionTerms[n].numbers[1],
				amberParameters->improperTorsionTerms[n].numbers[2],
				amberParameters->improperTorsionTerms[n].numbers[3],
				amberParameters->improperTorsionTerms[n].barrier,
				amberParameters->improperTorsionTerms[n].phase,
				amberParameters->improperTorsionTerms[n].n
				);
	}
	printf("\n");
	*/

	return TRUE;
			

}
/**********************************************************************/
static boolean readAmberOutOfPlaneTerms(AmberParameters* amberParameters,FILE* file)
{
	char t[BSIZE];
	int len = BSIZE;
	boolean Ok = FALSE;
	int n = 0;
	AmberOutOfPlaneTerms* terms = NULL;
	int nt = 0;

	/* Search Begin INPUT FOR OUT OF PLANE PARAMETERS */ 

	rewind(file);
	while(!feof(file))
	{
		if(fgets(t,len,file))
		{
			if(strstr(t,outOfPlaneTitle))
			{
				Ok = TRUE;
				break;
			}
		}
	}
	if(!Ok)
		return FALSE;

	terms = malloc(sizeof(AmberOutOfPlaneTerms));
	n = 0;
	Ok = FALSE;
	while(!feof(file))
	{
		if(fgets(t,len,file))
		{
			if(strstr(t,"End") || strstr(t,"END"))
			{
				Ok = TRUE;
				break;
			}
		}
		else 
		{
			Ok = FALSE;
			break;
		}



		terms[n].h3 = 0.0;
		terms[n].h4 = 0.0;
		terms[n].h5 = 0.0;
		terms[n].h6 = 0.0;
		nt = sscanf(t,"%d %d %d %d %lf %lf %lf %lf %lf",
			&terms[n].numbers[0],
			&terms[n].numbers[1],
			&terms[n].numbers[2],
			&terms[n].numbers[3],
			&terms[n].force,
			&terms[n].h3,
			&terms[n].h4,
			&terms[n].h5,
			&terms[n].h6);
		if(nt<9) terms[n].h6 = 0.0;
		if(nt<8) terms[n].h5 = 0.0;
		if(nt<7) terms[n].h4 = 0.0;
		if(nt<6) terms[n].h3 = 0.0;
		terms[n].type = 0;
		if(strstr(t,"W-D-C") || strstr(t,"Cross")) terms[n].type = 1;

	    terms[n].numbers[0]--;
	    terms[n].numbers[1]--;
	    terms[n].numbers[2]--;
	    terms[n].numbers[3]--;
	/*
	    if(terms[n].numbers[0]>terms[n].numbers[3])
	    {
		int t = terms[n].numbers[0];
		terms[n].numbers[0] = terms[n].numbers[3];
		terms[n].numbers[3] = t;

		t =  terms[n].numbers[1];
		terms[n].numbers[1] = terms[n].numbers[2];
		terms[n].numbers[2] = t;
	    }
	*/

		n++;
		terms = realloc(terms,(n+1)*sizeof(AmberOutOfPlaneTerms));
	}
	if(n==0 || !Ok ) free(terms);
	else
	{
		amberParameters->numberOfOutOfPlaneTerms = n;
		amberParameters->outOfPlaneTerms = terms;
	}
	/* printing for test*/
	/*
	printf("number of improper torsion terms = %d \n",
			amberParameters->numberOfOutOfPlaneTerms);

	for(n=0;n<amberParameters->numberOfOutOfPlaneTerms;n++)
	{
		printf("%d %d %d %d %f %f %f\n",
				amberParameters->outOfPlaneTerms[n].numbers[0],
				amberParameters->outOfPlaneTerms[n].numbers[1],
				amberParameters->outOfPlaneTerms[n].numbers[2],
				amberParameters->outOfPlaneTerms[n].numbers[3],
				amberParameters->outOfPlaneTerms[n].barrier,
				amberParameters->outOfPlaneTerms[n].phase,
				amberParameters->outOfPlaneTerms[n].n
				);
	}
	printf("\n");
	*/

	return TRUE;
			

}
/**********************************************************************/
static boolean readAmberHydrogenBonded1012Terms(AmberParameters* amberParameters,FILE* file)
{
	char t[BSIZE];
	int len = BSIZE;
	boolean Ok = FALSE;
	int n = 0;
	AmberHydrogenBonded1012Terms* terms = NULL;

	/* Search Begin INPUT FOR  ATOM TYPES */ 

	rewind(file);
	while(!feof(file))
	{
		if(fgets(t,len,file))
		{
			if(strstr(t,hydrogenBonded1012Title))
			{
				Ok = TRUE;
				break;
			}
		}
	}
	if(!Ok)
		return FALSE;

	terms = malloc(sizeof(AmberHydrogenBonded1012Terms));
	n = 0;
	Ok = FALSE;
	while(!feof(file))
	{
		if(fgets(t,len,file))
		{
			if(strstr(t,"End") || strstr(t,"END"))
			{
				Ok = TRUE;
				break;
			}
		}
		else 
		{
			Ok = FALSE;
			break;
		}

		sscanf(t,"%d %d %lf %lf",
				&terms[n].numbers[0],
				&terms[n].numbers[1],
				&terms[n].c,
				&terms[n].d);

		terms[n].numbers[0]--;
		terms[n].numbers[1]--;
		if(terms[n].numbers[0]>terms[n].numbers[1])
		{
			int t = terms[n].numbers[0];
			terms[n].numbers[0] = terms[n].numbers[1];
			terms[n].numbers[1] = t;
		}

		n++;
		terms = realloc(terms,(n+1)*sizeof(AmberHydrogenBonded1012Terms));
	}
	if(n==0 || !Ok )
		free(terms);
	else
	{
		amberParameters->numberOfHydrogenBonded1012 = n;
		amberParameters->hydrogenBonded1012Terms = terms;
	}
	/* printing for test*/
	/*
	printf("number of hydrogen bonds terms = %d \n",amberParameters->numberOfHydrogenBonded1012);
	for(n=0;n<amberParameters->numberOfHydrogenBonded1012;n++)
	{
		printf("%d %d %f %f\n",
				amberParameters->hydrogenBonded1012Terms[n].numbers[0],
				amberParameters->hydrogenBonded1012Terms[n].numbers[1],
				amberParameters->hydrogenBonded1012Terms[n].c,
				amberParameters->hydrogenBonded1012Terms[n].d
				);
	}
	printf("\n");
	*/

	return TRUE;
			

}
/**********************************************************************/
static boolean readAmberHydrogenBondedMorseTerms(AmberParameters* amberParameters,FILE* file)
{
	char t[BSIZE];
	int len = BSIZE;
	boolean Ok = FALSE;
	int n = 0;
	AmberHydrogenBondedMorseTerms* terms = NULL;

	/* Search Begin INPUT FOR  ATOM TYPES */ 

	rewind(file);
	while(!feof(file))
	{
		if(fgets(t,len,file))
		{
			if(strstr(t,hydrogenBondedMorseTitle))
			{
				Ok = TRUE;
				break;
			}
		}
	}
	if(!Ok)
		return FALSE;

	terms = malloc(sizeof(AmberHydrogenBondedMorseTerms));
	n = 0;
	Ok = FALSE;
	while(!feof(file))
	{
		if(fgets(t,len,file))
		{
			if(strstr(t,"End") || strstr(t,"END"))
			{
				Ok = TRUE;
				break;
			}
		}
		else 
		{
			Ok = FALSE;
			break;
		}

		sscanf(t,"%d %d %lf %lf %lf",
				&terms[n].numbers[0],
				&terms[n].numbers[1],
				&terms[n].force,
				&terms[n].Re,
				&terms[n].De
				);

		terms[n].numbers[0]--;
		terms[n].numbers[1]--;
		if(terms[n].numbers[0]>terms[n].numbers[1])
		{
			int t = terms[n].numbers[0];
			terms[n].numbers[0] = terms[n].numbers[1];
			terms[n].numbers[1] = t;
		}

		n++;
		terms = realloc(terms,(n+1)*sizeof(AmberHydrogenBondedMorseTerms));
	}
	if(n==0 || !Ok )
		free(terms);
	else
	{
		amberParameters->numberOfHydrogenBondedMorse = n;
		amberParameters->hydrogenBondedMorseTerms = terms;
	}
	/* printing for test*/
	/*
	printf("number of hydrogen bonds terms = %d \n",amberParameters->numberOfHydrogenBondedMorse);
	for(n=0;n<amberParameters->numberOfHydrogenBondedMorse;n++)
	{
		printf("%d %d %f %f %f\n",
				amberParameters->hydrogenBondedMorseTerms[n].numbers[0],
				amberParameters->hydrogenBondedMorseTerms[n].numbers[1],
				amberParameters->hydrogenBondedMorseTerms[n].force,
				amberParameters->hydrogenBondedMorseTerms[n].Re,
				amberParameters->hydrogenBondedMorseTerms[n].De,
				);
	}
	printf("\n");
	*/

	return TRUE;
			

}
/**********************************************************************/
static boolean readAmberSuttonChenTerms(AmberParameters* amberParameters,FILE* file)
{
	char t[BSIZE];
	int len = BSIZE;
	boolean Ok = FALSE;
	int n = 0;
	AmberSuttonChenTerms* terms = NULL;

	rewind(file);
	while(!feof(file))
	{
		if(fgets(t,len,file))
		{
			if(strstr(t,suttonChenTitle))
			{
				Ok = TRUE;
				break;
			}
		}
	}
	if(!Ok)
		return FALSE;

	terms = malloc(sizeof(AmberSuttonChenTerms));
	n = 0;
	Ok = FALSE;
	while(!feof(file))
	{
		if(fgets(t,len,file))
		{
			if(strstr(t,"End") || strstr(t,"END"))
			{
				Ok = TRUE;
				break;
			}
		}
		else 
		{
			Ok = FALSE;
			break;
		}

		sscanf(t,"%d %d %lf %lf %lf %lf %lf",
				&terms[n].numbers[0],
				&terms[n].numbers[1],
				&terms[n].epsilon,
				&terms[n].a,
				&terms[n].C,
				&terms[n].n,
				&terms[n].m
				);

		terms[n].numbers[0]--;
		terms[n].numbers[1]--;
		if(terms[n].numbers[0]>terms[n].numbers[1])
		{
			int t = terms[n].numbers[0];
			terms[n].numbers[0] = terms[n].numbers[1];
			terms[n].numbers[1] = t;
		}

		n++;
		terms = realloc(terms,(n+1)*sizeof(AmberSuttonChenTerms));
	}
	if(n==0 || !Ok )
		free(terms);
	else
	{
		amberParameters->numberOfSuttonChen = n;
		amberParameters->suttonChenTerms = terms;
	}
	/* printing for test*/
	/*
	printf("number of Sutton Chen terms = %d \n",amberParameters->numberOfSuttonChen);
	for(n=0;n<amberParameters->numberOfSuttonChen;n++)
	{
		printf("%d %d %f %f\n",
				amberParameters->suttonChenTerms[n].numbers[0],
				amberParameters->suttonChenTerms[n].numbers[1],
				amberParameters->suttonChenTerms[n].c,
				amberParameters->suttonChenTerms[n].d
				);
	}
	printf("\n");
	*/

	return TRUE;
			

}
/**********************************************************************/
static boolean readAmberVdw612Terms(AmberParameters* amberParameters,FILE* file)
{
	char t[BSIZE];
	int len = BSIZE;
	boolean Ok = FALSE;
	int n = 0;
	AmberVdw612Terms* terms = NULL;

	/* Search Begin INPUT FOR  NON-BONDED  */ 
	Ok = FALSE;
	rewind(file);
	while(!feof(file))
	{
		if(fgets(t,len,file))
		{
			if(strstr(t,vdw612Title))
			{
				Ok = TRUE;
				break;
			}
		}
	}
	if(!Ok)
		return FALSE;

	terms = malloc(sizeof(AmberVdw612Terms));
	n = 0;
	Ok = FALSE;
	while(!feof(file))
	{
		if(fgets(t,len,file))
		{
			if(strstr(t,"End") || strstr(t,"END"))
			{
				Ok = TRUE;
				break;
			}
		}
		else 
		{
			Ok = FALSE;
			break;
		}

		sscanf(t,"%d %lf %lf",
			&terms[n].number,
			&terms[n].r,
			&terms[n].epsilon);
		//printf("t=%s\n",t);
		
		terms[n].number--;
		n++;
		terms = realloc(terms,(n+1)*sizeof(AmberVdw612Terms));
	}

	if(n==0 || !Ok )
		free(terms);
	else
	{
		amberParameters->numberOfVdw612 = n;
		amberParameters->vdw612Terms = terms;
	}
	/* printing for test*/
	/*
	printf("number of non bended terms = %d \n",amberParameters->numberOfVdw612);
	for(n=0;n<amberParameters->numberOfVdw612;n++)
	{
		printf("%d %f %f\n",
				amberParameters->vdw612Terms[n].number,
				amberParameters->vdw612Terms[n].r,
				amberParameters->vdw612Terms[n].epsilon
				);
	}
	printf("\n");
	*/

	return TRUE;
			

}
/**********************************************************************/
static boolean readAmberVdw714Terms(AmberParameters* amberParameters,FILE* file)
{
	char t[BSIZE];
	int len = BSIZE;
	boolean Ok = FALSE;
	int n = 0;
	AmberVdw714Terms* terms = NULL;
	int nt = 0;

	/* Search Begin INPUT FOR  NON-BONDED 7-14  */ 
	Ok = FALSE;
	rewind(file);
	while(!feof(file))
	{
		if(fgets(t,len,file))
		{
			if(strstr(t,vdw714Title))
			{
				Ok = TRUE;
				break;
			}
		}
	}
	if(!Ok)
		return FALSE;

	terms = malloc(sizeof(AmberVdw714Terms));
	n = 0;
	Ok = FALSE;
	while(!feof(file))
	{
		if(fgets(t,len,file))
		{
			if(strstr(t,"End") || strstr(t,"END"))
			{
				Ok = TRUE;
				break;
			}
		}
		else 
		{
			Ok = FALSE;
			break;
		}

		nt = sscanf(t,"%d %lf %lf %lf %lf", &terms[n].number, &terms[n].r, &terms[n].epsilon, &terms[n].gamma, &terms[n].delta);
		if(nt<5) terms[n].delta = 0.07/2.0;
		if(nt<4) terms[n].gamma = 0.12/2.0;
		//printf("t=%s\n",t);
		
		terms[n].number--;
		n++;
		terms = realloc(terms,(n+1)*sizeof(AmberVdw714Terms));
	}

	if(n==0 || !Ok )
		free(terms);
	else
	{
		amberParameters->numberOfVdw714 = n;
		amberParameters->vdw714Terms = terms;
	}
	/* printing for test*/
	/*
	printf("number of non bended terms = %d \n",amberParameters->numberOfVdw714);
	for(n=0;n<amberParameters->numberOfVdw714;n++)
	{
		printf("%d %f %f\n",
				amberParameters->vdw714Terms[n].number,
				amberParameters->vdw714Terms[n].r,
				amberParameters->vdw714Terms[n].epsilon,
				amberParameters->vdw714Terms[n].gamma,
				amberParameters->vdw714Terms[n].delta
				);
	}
	printf("\n");
	*/

	return TRUE;
			

}
/**********************************************************************/
static boolean readAmberPairWiseTerms(AmberParameters* amberParameters,FILE* file)
{
	char t[BSIZE];
	int len = BSIZE;
	boolean Ok = FALSE;
	int n = 0;
	AmberPairWiseTerms* terms = NULL;

	/* Search Begin INPUT FOR PAIR WIZE  */ 
	Ok = FALSE;
	rewind(file);
	while(!feof(file))
	{
		if(fgets(t,len,file))
		{
			if(strstr(t,pairWiseTitle))
			{
				Ok = TRUE;
				break;
			}
		}
	}
	if(!Ok)
		return FALSE;

	terms = malloc(sizeof(AmberPairWiseTerms));
	n = 0;
	Ok = FALSE;
	while(!feof(file))
	{
		if(fgets(t,len,file))
		{
			if(strstr(t,"End") || strstr(t,"END"))
			{
				Ok = TRUE;
				break;
			}
		}
		else 
		{
			Ok = FALSE;
			break;
		}



		sscanf(t,"%d %d %lf %lf %lf %lf %lf %lf %lf",
			&terms[n].numbers[0],
			&terms[n].numbers[1],
			&terms[n].a,
			&terms[n].beta,
			&terms[n].c4,
			&terms[n].c6,
			&terms[n].c8,
			&terms[n].c10,
			&terms[n].b
			);
		
		terms[n].numbers[0]--;
		terms[n].numbers[1]--;
		n++;
		terms = realloc(terms,(n+1)*sizeof(AmberPairWiseTerms));
	}

	if(n==0 || !Ok )
		free(terms);
	else
	{
		amberParameters->numberOfPairWise = n;
		amberParameters->pairWiseTerms = terms;
	}
	/* printing for test*/
	/*
	printf("number of pair wise terms = %d \n",amberParameters->numberOfPairWise);
	for(n=0;n<amberParameters->numberOfPairWise;n++)
	{
		printf("%d %d %f %f %f %f %d\n",
				amberParameters->pairWiseTerms[n].numbers[0],
				amberParameters->pairWiseTerms[n].numbers[1],
				amberParameters->pairWiseTerms[n].a,
				amberParameters->pairWiseTerms[n].beta,
				amberParameters->pairWiseTerms[n].c6,
				amberParameters->pairWiseTerms[n].c8,
				amberParameters->pairWiseTerms[n].c10,
				amberParameters->pairWiseTerms[n].b,
				);
	}
	printf("\n");
	*/

	return TRUE;
			

}
/**********************************************************************/
static boolean readAmberBondHardnessTerms(AmberParameters* amberParameters,FILE* file)
{
	char t[BSIZE];
	int len = BSIZE;
	boolean Ok = FALSE;
	int n = 0;
	AmberBondHardnessTerms* terms = NULL;
	//int nt = 0;

	//printf("readAmberBondHardnessTerms\n");
	rewind(file);
	while(!feof(file))
	{
		if(fgets(t,len,file))
		{
			if(strstr(t,bondHardnessTitle))
			{
				Ok = TRUE;
				break;
			}
		}
	}
	if(!Ok)
		return FALSE;

	terms = malloc(sizeof(AmberBondHardnessTerms));
	n = 0;
	Ok = FALSE;
	while(!feof(file))
	{
		if(fgets(t,len,file))
		{
			if(strstr(t,"End") || strstr(t,"END"))
			{
				Ok = TRUE;
				break;
			}
		}
		else 
		{
			Ok = FALSE;
			break;
		}

		//nt = sscanf(t,"%d %d %lf", &terms[n].numbers[0], &terms[n].numbers[1], &terms[n].kappa);
		sscanf(t,"%d %d %lf", &terms[n].numbers[0], &terms[n].numbers[1], &terms[n].kappa);

	      	terms[n].numbers[0]--;
	      	terms[n].numbers[1]--;
		if(terms[n].numbers[0]>terms[n].numbers[1])
		{
			int t = terms[n].numbers[0];
			terms[n].numbers[0] = terms[n].numbers[1];
			terms[n].numbers[1] = t;
		}
		n++;
		terms = realloc(terms,(n+1)*sizeof(AmberBondHardnessTerms));
	}
	if(n==0 || !Ok )
		free(terms);
	else
	{
		amberParameters->numberOfHardnessTerms = n;
		amberParameters->bondHardnessTerms = terms;
	}
	/* printing for test*/
	
	/*
	printf("number of bonds = %d \n",amberParameters->numberOfHardnessTerms);
	for(n=0;n<amberParameters->numberOfHardnessTerms;n++)
	{
		printf("%d %d %f\n",
				amberParameters->bondHardnessTerms[n].numbers[0],
				amberParameters->bondHardnessTerms[n].numbers[1],
				amberParameters->bondHardnessTerms[n].kappa
				);
	}
	printf("\n");
	*/


	return TRUE;
			

}
/**********************************************************************/
boolean readAmberParameters(AmberParameters* amberParameters,char* filename)
{
	FILE* file;
	file = fopen(filename,"r");

	if(file == NULL)
		return FALSE;
	else
	{
		readAmberTypes(amberParameters,file);
		readAmberBondStretchTerms(amberParameters,file);
		readAmberAngleBendTerms(amberParameters,file);
		readAmberStrBendTerms(amberParameters,file);
		readAmberDihedralAngleTerms(amberParameters,file);
		readAmberImproperTorsionTerms(amberParameters,file);
		readAmberOutOfPlaneTerms(amberParameters,file);
		readAmberHydrogenBonded1012Terms(amberParameters,file);
		readAmberHydrogenBondedMorseTerms(amberParameters,file);
		readAmberSuttonChenTerms(amberParameters,file);
		readAmberVdw612Terms(amberParameters,file);
		readAmberVdw714Terms(amberParameters,file);
		readAmberPairWiseTerms(amberParameters,file);
		readAmberBondHardnessTerms(amberParameters,file);
		fclose(file);
	}
	return TRUE;
}
/**********************************************************************/
