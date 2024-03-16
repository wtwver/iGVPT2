
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

/* MolecularMechanics.c */
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>

#include "../Utils/Utils.h"
#include "../Utils/Timer.h"
#ifdef ENABLE_CL
#include "../Utils/CLProp.h"
#endif
#include "../Utils/AtomsProp.h"
#include "../Utils/Constants.h"
#include "../Utils/Types.h"
#include "../Molecule/Molecule.h"
#include "../MolecularMechanics/MolecularMechanics.h"
#include "../MolecularMechanics/LoadMMParameters.h"
#include "../MolecularMechanics/CreateMolecularMechanicsFile.h"
#include "../EmpriricalCorrections/WallCorrection.h"
#ifdef ENABLE_CL
#include "../MolecularMechanics/MolecularMechanicsCL.h"
#endif


#define POSEPS 1e-10
static AmberParameters* staticAmberParameters = NULL;


/* static void calculateGradientNumericAmber(ForceField* forceField);*/
static void calculateGradientAmber(ForceField* forceField);
static void calculateEnergyAmber(ForceField* forceField);
static double calculateEnergyTmpAmber(ForceField* forceField,Molecule* m);
static void printEnergies(ForceField* forceField);
double calculateEnergyCoulombAmber(ForceField* forceField,Molecule* molecule);

void dessine();

/*****************************************************************************/
void setH4CorrectionMM(FILE* file, ForceField* forceField) 
{
	char* fileName = NULL;

	if(readOneString(file,"H4Correction",&fileName) && fileName) 
	{
		char tmp[BSIZE];
		sprintf(tmp,"%s",fileName);
		uppercase(tmp);
		if(!strstr(tmp,"NONE"))
		{
			HyhrogenBondCorrectionParameters parameters;
			if(!strstr(tmp,"DEFAULT")) setHydrogenBondCorrectionParameters(&parameters, fileName, NULL);
			else setHydrogenBondCorrectionParameters(&parameters, NULL, "MM");
			forceField->H4Parameters = malloc(sizeof(HyhrogenBondCorrectionParameters));
			*(forceField->H4Parameters) = parameters;
			return;
		}
	}
	forceField->H4Parameters= NULL;
}
/**********************************************************************/
static void addH4Correction(ForceField* forceField,boolean addGradient)
{
	if(forceField->H4Parameters)
		forceField->molecule.potentialEnergy += getH4Correction(&forceField->molecule, forceField->H4Parameters, addGradient);
}
/**********************************************************************/
static double getH4Energy(ForceField* forceField, Molecule* molecule)
{
	if(forceField->H4Parameters) return getH4Correction(molecule, forceField->H4Parameters, FALSE);
	else return 0;
}
/**********************************************************************/
static void addWallCorrection(ForceField* forceField, boolean addGradient)
{
	if(forceField->options.addWallCorrection)
	{
		//double eOld = forceField->molecule.potentialEnergy;
		forceField->molecule.potentialEnergy += getWallCorrection(&forceField->molecule, addGradient);
		//fprintf(stdout," Wall Correction = %f\n", forceField->molecule.potentialEnergy-eOld);
	}
}
/**********************************************************************/
AmberParameters newAmberParameters()
{
	AmberParameters amberParameters;

	amberParameters.numberOfTypes = 0;
	amberParameters.atomTypes = NULL;

	amberParameters.numberOfStretchTerms = 0;
	amberParameters.bondStretchTerms = NULL;

	amberParameters.numberOfBendTerms = 0;
	amberParameters.angleBendTerms = NULL;

	amberParameters.numberOfStrBendTerms = 0;
	amberParameters.strBendTerms = NULL;

	amberParameters.numberOfDihedralTerms = 0;
	amberParameters.dihedralAngleTerms = NULL;

	amberParameters.numberOfImproperTorsionTerms = 0;
	amberParameters.improperTorsionTerms = NULL;

	amberParameters.numberOfOutOfPlaneTerms = 0;
	amberParameters.outOfPlaneTerms = NULL;

	amberParameters.numberOfVdw612 = 0;
	amberParameters.vdw612Terms = NULL;

	amberParameters.numberOfVdw714 = 0;
	amberParameters.vdw714Terms = NULL;

	amberParameters.numberOfHydrogenBonded1012 = 0;
	amberParameters.hydrogenBonded1012Terms = NULL;

	amberParameters.numberOfHydrogenBondedMorse = 0;
	amberParameters.hydrogenBondedMorseTerms = NULL;

	amberParameters.numberOfSuttonChen = 0;
	amberParameters.suttonChenTerms = NULL;

	amberParameters.numberOfPairWise = 0;
	amberParameters.pairWiseTerms = NULL;

	amberParameters.numberOfHardnessTerms = 0;
	amberParameters.bondHardnessTerms = NULL;


	return amberParameters;
	
}
/**********************************************************************/
static void freeAmberParameters(AmberParameters* amberParameters)
{
	int i;

	for(i=0;i<amberParameters->numberOfTypes;i++)
		if(amberParameters->atomTypes[i].name)
			free(amberParameters->atomTypes[i].name);

	amberParameters->numberOfTypes = 0;
	if(amberParameters->atomTypes )
		free(amberParameters->atomTypes );

	amberParameters->atomTypes = NULL;

	amberParameters->numberOfStretchTerms = 0;
	if(amberParameters->bondStretchTerms)
		free(amberParameters->bondStretchTerms);
	amberParameters->bondStretchTerms = NULL;

	amberParameters->numberOfBendTerms = 0;
	if(amberParameters->angleBendTerms) free(amberParameters->angleBendTerms);
	amberParameters->angleBendTerms = NULL;

	amberParameters->numberOfStrBendTerms = 0;
	if(amberParameters->strBendTerms) free(amberParameters->strBendTerms);
	amberParameters->strBendTerms = NULL;

	for(i=0;i<amberParameters->numberOfDihedralTerms;i++)
	{
		if(amberParameters->dihedralAngleTerms[i].divisor)
			free(amberParameters->dihedralAngleTerms[i].divisor);
		if(amberParameters->dihedralAngleTerms[i].barrier)
			free(amberParameters->dihedralAngleTerms[i].barrier);
		if(amberParameters->dihedralAngleTerms[i].phase)
			free(amberParameters->dihedralAngleTerms[i].phase);
		if(amberParameters->dihedralAngleTerms[i].n)
			free(amberParameters->dihedralAngleTerms[i].n);

	}

	amberParameters->numberOfDihedralTerms = 0;
	if(amberParameters->dihedralAngleTerms)
		free(amberParameters->dihedralAngleTerms);
	amberParameters->dihedralAngleTerms = NULL;

	amberParameters->numberOfImproperTorsionTerms = 0;
	if(amberParameters->improperTorsionTerms) free(amberParameters->improperTorsionTerms);
	amberParameters->improperTorsionTerms = NULL;

	amberParameters->numberOfOutOfPlaneTerms = 0;
	if(amberParameters->outOfPlaneTerms) free(amberParameters->outOfPlaneTerms);
	amberParameters->outOfPlaneTerms = NULL;

	amberParameters->numberOfVdw612 = 0;
	if(amberParameters->vdw612Terms) free(amberParameters->vdw612Terms);
	amberParameters->vdw612Terms = NULL;

	amberParameters->numberOfVdw714 = 0;
	if(amberParameters->vdw714Terms) free(amberParameters->vdw714Terms);
	amberParameters->vdw714Terms = NULL;

	amberParameters->numberOfHydrogenBonded1012 = 0;
	if(amberParameters->hydrogenBonded1012Terms) free(amberParameters->hydrogenBonded1012Terms);
	amberParameters->hydrogenBonded1012Terms = NULL;

	amberParameters->numberOfHydrogenBondedMorse = 0;
	if(amberParameters->hydrogenBondedMorseTerms) free(amberParameters->hydrogenBondedMorseTerms);
	amberParameters->hydrogenBondedMorseTerms = NULL;

	amberParameters->numberOfSuttonChen = 0;
	if(amberParameters->suttonChenTerms)
		free(amberParameters->suttonChenTerms);
	amberParameters->suttonChenTerms = NULL;

}
/**********************************************************************/
static double getChargeFromMM(AmberParameters* amberParameters, char* type)
{
        int i;
        int nTypes = amberParameters->numberOfTypes;
        AmberAtomTypes* types = amberParameters->atomTypes;
        int len = strlen(type);

        if(strcmp(type,"X")==0) return 0.0;
        for(i=0;i<nTypes;i++)
        {
                if(len == (int)strlen(types[i].name) && strstr(types[i].name,type))
                        return types[i].charge;

        }
        return 0.0;
}
/**********************************************************************/
static void scalChargesToSetSumEqTotalCharge(ForceField* forceField)
{
	int i;
	double sp = 0;
	double sm = 0;
	int np = 0;
	int nm = 0;
	Molecule* m = &forceField->molecule;
	double totalCharge = m->totalCharge;
	double scalp = 1.0;
	double scalm = 1.0;

	for (  i = 0; i < m->nAtoms; i++ )
	{
		if(m->atoms[i].charge>0) {sp += m->atoms[i].charge; np++;}
		if(m->atoms[i].charge<0) {sm += m->atoms[i].charge; nm++;}
	}
	if(fabs(totalCharge)>0 && fabs(sp+sm)>0) scalp = scalm = totalCharge/(sp+sm);
	if(fabs(totalCharge)==0)
	{
		if(fabs(sm)>fabs(sp)) scalm = (totalCharge-sp)/sm;
		else if(fabs(sp)!=0) scalp  = (totalCharge-sm)/sp;
	}
	if(fabs(scalm-1)>1e-10 || fabs(scalp-1)>1e-10) 
	{
		printf("Warning : charges stacletd to set sum equal to total charge\n");
		for (  i = 0; i < m->nAtoms; i++ )
		{
			if(m->atoms[i].charge>0) m->atoms[i].charge *= scalp;
			if(m->atoms[i].charge<0) m->atoms[i].charge *= scalm;
		}
	}
	sp = 0.0;
	sm = 0.0;
	for (  i = 0; i < m->nAtoms; i++ )
	{
		if(m->atoms[i].charge>0) {sp += m->atoms[i].charge;}
		if(m->atoms[i].charge<0) {sm += m->atoms[i].charge;}
	}
	printf("\nSum of partial charges = %14.8f , total charge = %14.8f\n\n",sp+sm,totalCharge);

}
/**********************************************************************/
static void setChargesFromMMParameters(AmberParameters* amberParameters, ForceField* forceField)
{
	int i;
	Molecule* m = &forceField->molecule;

	if(strcmp(forceField->options.chargesType,"MM")) return;

	printf("forceField->options.chargesType = %s\n",forceField->options.chargesType);
	printf("Warning : I use the MM charge \n");

	
	for (  i = 0; i < m->nAtoms; i++ )
	{
		m->atoms[i].charge = getChargeFromMM(amberParameters, m->atoms[i].mmType);
	}
	scalChargesToSetSumEqTotalCharge(forceField);
}
/**********************************************************************/
static void setCharges(ForceField* forceField)
{
	if(strstr(forceField->options.chargesType,"SCALED")) scalChargesToSetSumEqTotalCharge(forceField);
	if(strstr(forceField->options.chargesType,"EEM")) forceField->molecule.klass->setChargesEEM(&forceField->molecule);
	if(strstr(forceField->options.chargesType,"ACKS2")) forceField->molecule.klass->setChargesACKS2(&forceField->molecule);
}
/**********************************************************************/
static int getNumberType(AmberParameters* amberParameters, char* type)
{
	int i;
	int nTypes = amberParameters->numberOfTypes;
	AmberAtomTypes* types = amberParameters->atomTypes;
	int len = strlen(type);

	if(strcmp(type,"X")==0)
		return -1;
	for(i=0;i<nTypes;i++)
	{
		if(len == (int)strlen(types[i].name) && strstr(types[i].name,type))
			return types[i].number;

	}
	return -2;
}
/**********************************************************************/
static ForceField newAmberModel()
{
	ForceField forceField = newForceField();

	forceField.klass->calculateGradient = calculateGradientAmber;
	/*forceField.klass->calculateGradient = calculateGradientNumericAmber;*/
	forceField.klass->calculateEnergy = calculateEnergyAmber;
	forceField.klass->calculateEnergyTmp = calculateEnergyTmpAmber;
	forceField.klass->printEnergies = printEnergies;

	forceField.options.type = AMBER;
	forceField.options.coulomb = TRUE;
	forceField.options.hydrogenBonded1012 = TRUE;
	forceField.options.hydrogenBonded612 = FALSE;
	forceField.options.hydrogenBondedMorse = TRUE;
	forceField.options.suttonChen = TRUE;
	forceField.options.improperTorsion = TRUE;
	forceField.options.outOfPlane = TRUE;

	return forceField;

}
/**********************************************************************/
static ForceField newPairWiseModel()
{
	ForceField forceField = newForceField();

	forceField.klass->calculateGradient = calculateGradientAmber;
	forceField.klass->calculateEnergy = calculateEnergyAmber;
	forceField.klass->calculateEnergyTmp = calculateEnergyTmpAmber;
	forceField.klass->printEnergies = printEnergies;

	forceField.options.type = PAIRWISE;

	forceField.options.coulomb = TRUE;
	forceField.options.vdw612 = TRUE;
	forceField.options.vdw714 = FALSE;

	forceField.options.bondStretch = FALSE;
	forceField.options.angleBend = FALSE;
	forceField.options.strBend = FALSE;
	forceField.options.dihedralAngle = FALSE;
	forceField.options.improperTorsion = FALSE;
	forceField.options.outOfPlane = FALSE;
	forceField.options.hydrogenBonded612 = FALSE;
	forceField.options.hydrogenBonded1012 = FALSE;
	forceField.options.hydrogenBondedMorse = FALSE;
	forceField.options.suttonChen = FALSE;


	return forceField;

}
/**********************************************************************/
/*
static boolean isIonic(char* mmType)
{
	if(!strcmp(mmType,"Li")) return TRUE;
	if(!strcmp(mmType,"Na")) return TRUE;
	if(!strcmp(mmType,"K")) return TRUE;
	if(!strcmp(mmType,"Rb")) return TRUE;
	if(!strcmp(mmType,"Cs")) return TRUE;
	if(!strcmp(mmType,"Ca")) return TRUE;
	if(!strcmp(mmType,"Sr")) return TRUE;
	if(!strcmp(mmType,"Ba")) return TRUE;
	if(!strcmp(mmType,"Zn")) return TRUE;
	if(!strcmp(mmType,"IB")) return TRUE;
	if(!strcmp(mmType,"Cl")) return TRUE;
	return FALSE;
}
*/
/**********************************************************************/
static boolean getStretchParameters(AmberParameters* amberParameters,
				int a1Type, int a2Type, 
				double* forceOrAlpha, double* equilibriumDistance,
				double* h3OrDe, double* h4, double* h5, double* h6,
				int* type
				)
{
	int i;
	forceOrAlpha[0] = 0.0;
	equilibriumDistance[0] = 0.0;
	h3OrDe[0] = 0.0;
	h4[0] = 0.0;
	h5[0] = 0.0;
	h6[0] = 0.0;
	*type = 0;

	if(a1Type>a2Type)
	{
		int t;
		t = a1Type;
		a1Type = a2Type;
		a2Type = t;
	}

	for(i=0;i<amberParameters->numberOfStretchTerms;i++)
	{
		if(
			a1Type == amberParameters->bondStretchTerms[i].numbers[0]  &&
			a2Type == amberParameters->bondStretchTerms[i].numbers[1] 
		)
		{
			*type = amberParameters->bondStretchTerms[i].type;
			if(*type==0)
			{
				forceOrAlpha[0]       = amberParameters->bondStretchTerms[i].forceConstant;// force for harmonic
			}
			else
			{
				double De =  amberParameters->bondStretchTerms[i].h3;
				if(fabs(De)>1e-13) forceOrAlpha[0] = sqrt(amberParameters->bondStretchTerms[i].forceConstant/2.0/De);// alpha
				else forceOrAlpha[0]       = 0.0;
			}
			equilibriumDistance[0] = amberParameters->bondStretchTerms[i].equilibriumDistance;
			h3OrDe[0] = amberParameters->bondStretchTerms[i].h3; // h3 for harmonic, De for Morse
			h4[0] = amberParameters->bondStretchTerms[i].h4;
			h5[0] = amberParameters->bondStretchTerms[i].h5;
			h6[0] = amberParameters->bondStretchTerms[i].h6;
			return TRUE;
		}
	}
	return FALSE;
}
/**********************************************************************/
static boolean getBendParameters(AmberParameters* amberParameters,int a1Type, int a2Type, int a3Type,
	       	double* forceConstant, double* equilibriumAngle,
		double* h3, double* h4, double* h5, double* h6
		)
{
	int i;
	forceConstant[0] = 0.0;
	equilibriumAngle[0] = 0.0;
	h3[0] = 0.0;
	h4[0] = 0.0;
	h5[0] = 0.0;
	h6[0] = 0.0;

	if(a1Type>a3Type)
	{
		int t;
		t = a1Type;
		a1Type = a3Type;
		a3Type = t;
	}

	for(i=0;i<amberParameters->numberOfBendTerms;i++)
	{
		if(
			a1Type == amberParameters->angleBendTerms[i].numbers[0]  &&
			a2Type == amberParameters->angleBendTerms[i].numbers[1]  &&
			a3Type == amberParameters->angleBendTerms[i].numbers[2] 
		)
		{
			forceConstant[0]       = amberParameters->angleBendTerms[i].forceConstant;
			equilibriumAngle[0]    = amberParameters->angleBendTerms[i].equilibriumAngle;
			h3[0]    = amberParameters->angleBendTerms[i].h3;
			h4[0]    = amberParameters->angleBendTerms[i].h4;
			h5[0]    = amberParameters->angleBendTerms[i].h5;
			h6[0]    = amberParameters->angleBendTerms[i].h6;
			return TRUE;
		}
	}
	return FALSE;
}
/**********************************************************************/
static boolean getStrBendParameters(AmberParameters* amberParameters,int a1Type, int a2Type, int a3Type,
	       	double* forceConstant12, double* forceConstant23
		)
{
	int i;
	forceConstant12[0] = 0.0;
	forceConstant23[0] = 0.0;

	if(a1Type>a3Type)
	{
		int t;
		t = a1Type;
		a1Type = a3Type;
		a3Type = t;
	}

	for(i=0;i<amberParameters->numberOfStrBendTerms;i++)
	{
		if(
			a1Type == amberParameters->strBendTerms[i].numbers[0]  &&
			a2Type == amberParameters->strBendTerms[i].numbers[1]  &&
			a3Type == amberParameters->strBendTerms[i].numbers[2] 
		)
		{
			forceConstant12[0]       = amberParameters->strBendTerms[i].forceConstant12;
			forceConstant23[0]       = amberParameters->strBendTerms[i].forceConstant23;
			return TRUE;
		}
	}
	return FALSE;
}
/**********************************************************************/
static boolean getHydrogenBonded1012Parameters(AmberParameters* amberParameters, int a1Type, int a2Type, double c[], double d[])
{
	int i;

	c[0] = 0.0;
	d[0] = 0.0;

	for(i=0;i<amberParameters->numberOfHydrogenBonded1012;i++)
	{
		if(
			(a1Type == amberParameters->hydrogenBonded1012Terms[i].numbers[0] &&
			a2Type == amberParameters->hydrogenBonded1012Terms[i].numbers[1])
			||
			(a2Type == amberParameters->hydrogenBonded1012Terms[i].numbers[0] &&
			a1Type == amberParameters->hydrogenBonded1012Terms[i].numbers[1])
		  )
		{
			c[0]    = amberParameters->hydrogenBonded1012Terms[i].c;
			d[0]    = amberParameters->hydrogenBonded1012Terms[i].d;
			return TRUE;
		}
	}
	return FALSE;
}
/**********************************************************************/
static boolean getHydrogenBondedMorseParameters(AmberParameters* amberParameters, int a1Type, int a2Type, double alpha[], double Re[], double De[])
{
	int i;

	alpha[0] = 0.0;
	Re[0] = 0.0;
	De[0] = 0.0;

	for(i=0;i<amberParameters->numberOfHydrogenBondedMorse;i++)
	{
		if(
			(a1Type == amberParameters->hydrogenBondedMorseTerms[i].numbers[0] &&
			a2Type == amberParameters->hydrogenBondedMorseTerms[i].numbers[1])
			||
			(a2Type == amberParameters->hydrogenBondedMorseTerms[i].numbers[0] &&
			a1Type == amberParameters->hydrogenBondedMorseTerms[i].numbers[1])
		  )
		{
			De[0]    = amberParameters->hydrogenBondedMorseTerms[i].De;
			if(fabs(De[0])>1e-13) alpha[0] = sqrt(fabs(amberParameters->hydrogenBondedMorseTerms[i].force/2.0/De[0]));
			Re[0]    = amberParameters->hydrogenBondedMorseTerms[i].Re;
			return TRUE;
		}
	}
	return FALSE;
}
/**********************************************************************/
static boolean getSuttonChenParameters(AmberParameters* amberParameters, int a1Type, int a2Type, double epsilon[], double a[], double C[], double n[], double m[] )
{
	int i;
	epsilon[0] = 0.0;
	a[0] = 0.0;
	C[0] = 0.0;
	n[0] = 0.0;
	m[0] = 0.0;

	for(i=0;i<amberParameters->numberOfSuttonChen;i++)
	{
		if(
			a1Type == amberParameters->suttonChenTerms[i].numbers[0] &&
			a2Type == amberParameters->suttonChenTerms[i].numbers[1]
		  )
		{
			epsilon[0]    = amberParameters->suttonChenTerms[i].epsilon;
			a[0]    = amberParameters->suttonChenTerms[i].a;
			C[0]    = amberParameters->suttonChenTerms[i].C;
			n[0]    = amberParameters->suttonChenTerms[i].n;
			m[0]    = amberParameters->suttonChenTerms[i].m;
			return TRUE;
		}
	}
	return FALSE;
}
/**********************************************************************/
static boolean getVdw612Parameters(AmberParameters* amberParameters, int atomType, double* r, double* epsilon )
{

	int i;
	r[0] = 1.0;
	epsilon[0] = 0.0;
	
	for(i=0;i<amberParameters->numberOfVdw612;i++)
	{
	//printf("r = %f eps = %f\n",amberParameters->vdw612Terms[i].r, amberParameters->vdw612Terms[i].epsilon);
		if(
			atomType == amberParameters->vdw612Terms[i].number
		  )
		{
			r[0]       = amberParameters->vdw612Terms[i].r;
			epsilon[0]    = amberParameters->vdw612Terms[i].epsilon;
			/*printf("r = %f eps = %f\n",r[0],epsilon[0]);*/
			return TRUE;
		}
	}

	return FALSE;
}
/********************************************************************************************************************************************************/
static boolean getHydrogenBonded612Parameters(AmberParameters* amberParameters, int a1Type, int a2Type, double c[], double d[] , FILE* logFile, char* type1, char* type2)
{
	double equilibriumDistance, epsilon, epsilonProduct, ri,rj;
	double Bij , Aij;

	c[0] = 0.0;
	d[0] = 0.0;
	if ( ! ( getVdw612Parameters(amberParameters, a1Type, &equilibriumDistance, &epsilon ) ) )
	{
		fprintf(logFile, ("**** couldn't find non bonded parameters for %s \n"),type1);
		fflush(logFile);
		return FALSE;
	}
		
	epsilonProduct = sqrt(fabs(epsilon));
	ri = equilibriumDistance;
	if(!getVdw612Parameters(amberParameters, a2Type, &equilibriumDistance, &epsilon )) 
	{
		fprintf(logFile, ("**** couldn't find non bonded parameters for %s \n"), type2);
		fflush(logFile);
		return FALSE;
	}
	epsilonProduct *= sqrt(fabs(epsilon));
	rj = equilibriumDistance;
	Bij = ( ri + rj ) * ( ri + rj );
	Bij = Bij * Bij * Bij;
	Aij = Bij * Bij * epsilonProduct;
	Bij *= epsilonProduct * 2.0;
	c[0] = Aij;
	d[0] = Bij;
	return TRUE;
}
/**********************************************************************/
static boolean getVdw714Parameters(AmberParameters* amberParameters, int atomType, double* r, double* epsilon , double* gamma, double* delta)
{

	int i;
	r[0] = 1.0;
	epsilon[0] = 0.0;
	gamma[0] = 0.07/2.0;
	delta[0] = 0.12/2.0;
	
	for(i=0;i<amberParameters->numberOfVdw714;i++)
	{
	//printf("r = %f eps = %f\n",amberParameters->vdw714Terms[i].r, amberParameters->vdw714Terms[i].epsilon);
		if( atomType == amberParameters->vdw714Terms[i].number)
		{
			r[0]       = amberParameters->vdw714Terms[i].r;
			epsilon[0]    = amberParameters->vdw714Terms[i].epsilon;
			gamma[0]    = amberParameters->vdw714Terms[i].gamma;
			delta[0]    = amberParameters->vdw714Terms[i].delta;
			/*printf("r = %f eps = %f\n",r[0],epsilon[0]);*/
			return TRUE;
		}
	}

	return FALSE;
}
/**********************************************************************/
static boolean getPairWiseParameters(AmberParameters* amberParameters,
	       	int a1Type, int a2Type,
		double* a, double* beta,
	       	 double* c4, double* c6, double* c8, double* c10, double* b)
{

	int i;

	a[0]    = 0.0;
	beta[0] = 1.0;
	c4[0]   = 0.0;
	c6[0]   = 0.0;
	c8[0]   = 0.0;
	c10[0]   = 0.0;
	b[0]    = 1.0;
	for(i=0;i<amberParameters->numberOfPairWise;i++)
	{
		if(
			(
			a1Type == amberParameters->pairWiseTerms[i].numbers[0] &&
			a2Type == amberParameters->pairWiseTerms[i].numbers[1] 
			) ||
			(
			a1Type == amberParameters->pairWiseTerms[i].numbers[1] &&
			a2Type == amberParameters->pairWiseTerms[i].numbers[0]
			)

		  )
		{
			a[0]    = amberParameters->pairWiseTerms[i].a;
			beta[0]    = amberParameters->pairWiseTerms[i].beta;
			c4[0]    = amberParameters->pairWiseTerms[i].c4;
			c6[0]    = amberParameters->pairWiseTerms[i].c6;
			c8[0]    = amberParameters->pairWiseTerms[i].c8;
			c10[0]    = amberParameters->pairWiseTerms[i].c10;
			b[0]    = amberParameters->pairWiseTerms[i].b;
			return TRUE;
		}
	}

	return FALSE;
}
/*********************************************************************************************************************************/
static int getNumberOutOfPlaneParameters( AmberParameters* amberParameters, int a1Type, int a2Type, int a3Type, int a4Type)
{
	int i;
	int a1Typet;
	int a2Typet;
	int a3Typet;
	int a4Typet;

	int a1Typet2;
	int a2Typet2;
	int a3Typet2;
	int a4Typet2;
	boolean btype;
	boolean Ok;
	int types[4];
	int k;

	a1Typet = a4Type;
	a2Typet = a2Type;
	a3Typet = a3Type;
	a4Typet = a1Type;

	a1Typet2 = a3Type;
	a2Typet2 = a2Type;
	a3Typet2 = a1Type;
	a4Typet2 = a4Type;

	/*
	for(i=0;i<amberParameters->numberOfOutOfPlaneTerms;i++)
	{
		for(k=0;k<4;k++) printf("%d ", amberParameters->outOfPlaneTerms[i].numbers[k]);
		printf("\n");
	}
	*/
	/* Je cherche d'abord sans les -1 */
	for(i=0;i<amberParameters->numberOfOutOfPlaneTerms;i++)
	{

		types[0] = a1Type;
		types[1] = a2Type;
		types[2] = a3Type;
		types[3] = a4Type;

		Ok = TRUE;
		for(k=0;k<4;k++)
		{
			btype = (types[k] == amberParameters->outOfPlaneTerms[i].numbers[k]);
			if(!btype)
			{
				Ok = FALSE;
				break;
			}
		}
		if(!Ok)
		{
			types[0] = a1Typet;
			types[1] = a2Typet;
			types[2] = a3Typet;
			types[3] = a4Typet;
			Ok = TRUE;
			for(k=0;k<4;k++)
			{
				btype = (types[k] == amberParameters->outOfPlaneTerms[i].numbers[k]);
				if(!btype)
				{
					Ok = FALSE;
					break;
				}
			}
		}
		if(!Ok)
		{
			types[0] = a1Typet2;
			types[1] = a2Typet2;
			types[2] = a3Typet2;
			types[3] = a4Typet2;
			Ok = TRUE;
			for(k=0;k<4;k++)
			{
				btype = (types[k] == amberParameters->outOfPlaneTerms[i].numbers[k]);
				if(!btype)
				{
					Ok = FALSE;
					break;
				}
			}
		}

			 
		if(Ok)
		{
			return i;
		}
	}
	/* Je cherche d'abord avec les -1 */
	for(i=0;i<amberParameters->numberOfOutOfPlaneTerms;i++)
	{

		types[0] = a1Type;
		types[1] = a2Type;
		types[2] = a3Type;
		types[3] = a4Type;

		Ok = TRUE;
		for(k=0;k<4;k++)
		{
			btype = 
			(amberParameters->outOfPlaneTerms[i].numbers[k] == -1) || 
			(types[k] == amberParameters->outOfPlaneTerms[i].numbers[k]);
			if(!btype)
			{
				Ok = FALSE;
				break;
			}
		}
		if(!Ok)
		{
			types[0] = a1Typet;
			types[1] = a2Typet;
			types[2] = a3Typet;
			types[3] = a4Typet;
			Ok = TRUE;
			for(k=0;k<4;k++)
			{
				btype = 
				(amberParameters->outOfPlaneTerms[i].numbers[k] == -1) || 
				(types[k] == amberParameters->outOfPlaneTerms[i].numbers[k]);
				if(!btype)
				{
					Ok = FALSE;
					break;
				}
			}
		}
		if(!Ok)
		{
			types[0] = a1Typet2;
			types[1] = a2Typet2;
			types[2] = a3Typet2;
			types[3] = a4Typet2;
			Ok = TRUE;
			for(k=0;k<4;k++)
			{
				btype = 
				(amberParameters->outOfPlaneTerms[i].numbers[k] == -1) || 
				(types[k] == amberParameters->outOfPlaneTerms[i].numbers[k]);
				if(!btype)
				{
					Ok = FALSE;
					break;
				}
			}
		}


			 
		if(Ok)
		{
			return i;
		}
	}

	return -1;
}
/**********************************************************************/
static boolean addOutOfPlaneParameters(AmberParameters* amberParameters, int a1, int a2, int a3, int a4, int* atomTypes, double** outOfPlaneTerms, int *pNumberOfOutOfPlaneTerms)
{
	int numberOfOutOfPlaneTerms = *pNumberOfOutOfPlaneTerms;
	
	int a1Type, a2Type, a3Type,a4Type;
	int k;
	a1Type = atomTypes[a1];
	a2Type = atomTypes[a2];
	a3Type = atomTypes[a3];
	a4Type = atomTypes[a4];
	//printf("type a1 a2 a3 a4 = %d %d %d %d\n",a1Type, a2Type, a3Type, a4Type);
	k = getNumberOutOfPlaneParameters(amberParameters, a1Type, a2Type, a3Type, a4Type);
	//printf("k =%d\n",k);
	if(k>-1)
	{
		outOfPlaneTerms[0][numberOfOutOfPlaneTerms] = a1;
		outOfPlaneTerms[1][numberOfOutOfPlaneTerms] = a2;
		outOfPlaneTerms[2][numberOfOutOfPlaneTerms] = a3;
		outOfPlaneTerms[3][numberOfOutOfPlaneTerms] = a4;
		outOfPlaneTerms[4][numberOfOutOfPlaneTerms] = amberParameters->outOfPlaneTerms[k].type;// 0 = ALLINGER, 1 = WDC
		outOfPlaneTerms[5][numberOfOutOfPlaneTerms] = amberParameters->outOfPlaneTerms[k].force;
		outOfPlaneTerms[6][numberOfOutOfPlaneTerms] = amberParameters->outOfPlaneTerms[k].h3;
		outOfPlaneTerms[7][numberOfOutOfPlaneTerms] = amberParameters->outOfPlaneTerms[k].h4;
		outOfPlaneTerms[8][numberOfOutOfPlaneTerms] = amberParameters->outOfPlaneTerms[k].h5;
		outOfPlaneTerms[9][numberOfOutOfPlaneTerms] = amberParameters->outOfPlaneTerms[k].h6;
		numberOfOutOfPlaneTerms++;
		*pNumberOfOutOfPlaneTerms = numberOfOutOfPlaneTerms;
		return TRUE;
	}
	return FALSE;
}
/**********************************************************************/
static void setOutOfPlaneParameters(AmberParameters* amberParameters, ForceField* forceField, int* atomTypes)
{
	int i;
	int j;
	int a1,a2,a3,a4;
	Molecule* m = &forceField->molecule;
	double* outOfPlaneTerms[OUTOFPLANEDIM];
	int numberOfOutOfPlaneTerms = 0;
	int a4j;
	int nc = 0;

	forceField->numberOfOutOfPlaneTerms = 0;
	for( i=0; i<OUTOFPLANEDIM;i++) forceField->outOfPlaneTerms[i] = NULL;
	if(!forceField->options.outOfPlane) return;

	/*  6 terms 1=a1, 2=a2, 3=a3, 4=a4, 5=force  6=h3 7=h4 8=h5 9=h6*/
	/*  OUTOFPLANEDIM 9 */

	for( i=0; i<OUTOFPLANEDIM;i++) outOfPlaneTerms[i] =  malloc(3*m->numberOf3Connections*sizeof(double)); 

	numberOfOutOfPlaneTerms = 0;

	for (  i = 0; i < m->numberOf3Connections; i++ )
	{
		a1 = m->connected3[0][i];
		a2 = m->connected3[1][i];
		a3 = m->connected3[2][i];
		a4 = -1;
		nc = 0;
		for (  j = 0; j < m->numberOf2Connections; j++ ) if(m->connected2[0][j] == a2 || m->connected2[1][j] == a2) nc++;
		//printf("a2 = %d nc = %d\n",a2,nc);
		if(nc!=3) continue;

		for (  j = 0; j < m->numberOf2Connections; j++ )
		{
			if(m->connected2[0][j] != a2 && m->connected2[1][j] != a2) continue;
			a4j = m->connected2[0][j]; 
			if(m->connected2[0][j] == a2) a4j = m->connected2[1][j];
			if(a4j == a1) continue;
			if(a4j == a3) continue;
			a4 = a4j;
			break;
		}
		//printf("a1 a2 a3 a4 = %d %d %d %d\n",a1, a2, a3, a4);
		if(a4<0) continue;

		addOutOfPlaneParameters(amberParameters, a1, a2, a3, a4, atomTypes, outOfPlaneTerms, &numberOfOutOfPlaneTerms);
		addOutOfPlaneParameters(amberParameters, a3, a2, a1, a4, atomTypes, outOfPlaneTerms, &numberOfOutOfPlaneTerms);
		addOutOfPlaneParameters(amberParameters, a4, a2, a3, a1, atomTypes, outOfPlaneTerms, &numberOfOutOfPlaneTerms);
	}

	forceField->numberOfOutOfPlaneTerms = numberOfOutOfPlaneTerms;
	for( i=0; i<OUTOFPLANEDIM;i++) forceField->outOfPlaneTerms[i] = outOfPlaneTerms[i]; 
	//printf("numberOfOutOfPlaneTerms=%d\n",numberOfOutOfPlaneTerms);
	//printf("numberOfOutOfPlaneParams=%d\n",amberParameters->numberOfOutOfPlaneTerms);
}
/**********************************************************************/
static int getNumberDihedralParameters( AmberParameters* amberParameters,
		int a1Type, int a2Type, int a3Type, int a4Type,
		int *n)
{
	int i;
	int a1Typet;
	int a2Typet;
	int a3Typet;
	int a4Typet;
	boolean btype;
	boolean Ok;
	int types[4];
	int k;

	*n = 0;

	a1Typet = a4Type;
	a2Typet = a3Type;
	a3Typet = a2Type;
	a4Typet = a1Type;

	/* Je cherche d'abord sans les -1 */
	for(i=0;i<amberParameters->numberOfDihedralTerms;i++)
	{

		types[0] = a1Type;
		types[1] = a2Type;
		types[2] = a3Type;
		types[3] = a4Type;

		Ok = TRUE;
		for(k=0;k<4;k++)
		{
			btype = (types[k] == amberParameters->dihedralAngleTerms[i].numbers[k]);
			if(!btype)
			{
				Ok = FALSE;
				break;
			}
		}
		if(!Ok)
		{
			types[0] = a1Typet;
			types[1] = a2Typet;
			types[2] = a3Typet;
			types[3] = a4Typet;
			Ok = TRUE;
			for(k=0;k<4;k++)
			{
				btype = (types[k] == amberParameters->dihedralAngleTerms[i].numbers[k]);
				if(!btype)
				{
					Ok = FALSE;
					break;
				}
			}
		}

			 
		if(Ok)
		{
			*n =i;
			return amberParameters->dihedralAngleTerms[i].nSomme;
		}
	}
	/* Je cherche d'abord avec les -1 */
	for(i=0;i<amberParameters->numberOfDihedralTerms;i++)
	{

		types[0] = a1Type;
		types[1] = a2Type;
		types[2] = a3Type;
		types[3] = a4Type;

		Ok = TRUE;
		for(k=0;k<4;k++)
		{
			btype = 
			(amberParameters->dihedralAngleTerms[i].numbers[k] == -1) || 
			(types[k] == amberParameters->dihedralAngleTerms[i].numbers[k]);
			if(!btype)
			{
				Ok = FALSE;
				break;
			}
		}
		if(!Ok)
		{
			types[0] = a1Typet;
			types[1] = a2Typet;
			types[2] = a3Typet;
			types[3] = a4Typet;
			Ok = TRUE;
			for(k=0;k<4;k++)
			{
				btype = 
				(amberParameters->dihedralAngleTerms[i].numbers[k] == -1) || 
				(types[k] == amberParameters->dihedralAngleTerms[i].numbers[k]);
				if(!btype)
				{
					Ok = FALSE;
					break;
				}
			}
		}

			 
		if(Ok)
		{
			*n =i;
			return amberParameters->dihedralAngleTerms[i].nSomme;
		}
	}

	return 0;
}

/**********************************************************************/
static boolean canHydrogenBond(AmberParameters* amberParameters, int a1Type, int a2Type )
{
	AmberAtomTypes* types = amberParameters->atomTypes;
	int nTypes = amberParameters->numberOfTypes;
	int a1=-1,a2=-1;
	int i;

	if( a1Type<0 || a2Type<0) return FALSE;
	for(i=0;i<nTypes;i++)
	{
		if(types[i].number==a1Type) a1 = i;
		if(types[i].number==a2Type) a2 = i;
	}
	if( a1<0 || a2<0) return FALSE;
	if( types[a1].name[0] == 'H' || types[a2].name[0] == 'H') return TRUE;

	return FALSE;
}
/**********************************************************************/
static boolean canSuttonChen(ForceField* forceField, int a)
{

	int i;
	for (  i = 0; i < forceField->numberOfSuttonChen; i++ )
	{
		int b1,b2;
		b1 = (int)forceField->suttonChenTerms[0][i];
		b2 = (int)forceField->suttonChenTerms[1][i];
		if(a==b1 || a==b2) return TRUE;
	}

	return FALSE;
}
/**********************************************************************/
static void setRattleConstraintsParameters(ForceField* forceField)
{
	Molecule* m = &forceField->molecule;
	m->klass->resetConstraints(m,forceField->options.rattleConstraints);
}
/**********************************************************************/
static void setStretchParameters(AmberParameters* amberParameters,ForceField* forceField,int* atomTypes)
{
	int i;
	int a1,a2;
	int type;
	int a1Type, a2Type;
	double forceOrAlpha, equilibriumDistance;
	double h3OrDe,h4,h5,h6;
	Molecule* m = &forceField->molecule;
	int numberOfStretchTerms = 0;
	double* bondStretchTerms[STRETCHDIM];

	if(!forceField->options.bondStretch) return;

	numberOfStretchTerms = m->numberOf2Connections;
	for( i=0; i<STRETCHDIM;i++)
       		bondStretchTerms[i] = malloc(numberOfStretchTerms*sizeof(double));

	/* 1=a1, 2=a2, 3=Force, 4=Re  h3 h4 h5 h6*/
	/* STRETCHDIM 	8 */
	for ( i = 0; i < numberOfStretchTerms; i++ )
	{
		a1 = m->connected2[0][i];
		a2 = m->connected2[1][i];
		a1Type = atomTypes[a1];
		a2Type = atomTypes[a2];
		
		if ( ! (getStretchParameters(amberParameters, a1Type, a2Type,&forceOrAlpha,&equilibriumDistance,&h3OrDe,&h4,&h5,&h6,&type) ) )
		{
			/*
			char l1 = m->atoms[a1].mmType[0];
			char l2 = m->atoms[a2].mmType[0];
			*/
			fprintf(forceField->logfile,  ("**** couldn't find stretch parameters for %s-%s(%d-%d) "), 
				m->atoms[a1].mmType,m->atoms[a2].mmType,a1Type, a2Type);
			fprintf(forceField->logfile, "\n");
			fflush(forceField->logfile);
			bondStretchTerms[0][i] = type;
			bondStretchTerms[1][i] = a1;
			bondStretchTerms[2][i] = a2;
			bondStretchTerms[3][i] = forceOrAlpha;
			bondStretchTerms[4][i] = equilibriumDistance;
			bondStretchTerms[5][i] = h3OrDe;
			bondStretchTerms[6][i] = h4;
			bondStretchTerms[7][i] = h5;
			bondStretchTerms[8][i] = h6;
		}
		else
		{
			bondStretchTerms[0][i] = type;
			bondStretchTerms[1][i] = a1;
			bondStretchTerms[2][i] = a2;
			bondStretchTerms[3][i] = forceOrAlpha;
			bondStretchTerms[4][i] = equilibriumDistance;
			bondStretchTerms[5][i] = h3OrDe;
			bondStretchTerms[6][i] = h4;
			bondStretchTerms[7][i] = h5;
			bondStretchTerms[8][i] = h6;
		}
	}

	forceField->numberOfStretchTerms = numberOfStretchTerms;
	for( i=0; i<STRETCHDIM;i++)
       		forceField->bondStretchTerms[i] = bondStretchTerms[i]; 
}
/**********************************************************************/
static void setBendParameters(AmberParameters* amberParameters,ForceField* forceField,int* atomTypes)
{
	int i;
	int a1,a2,a3;
	int a1Type, a2Type, a3Type;
	Molecule* m = &forceField->molecule;
	int numberOfBendTerms = 0;
	double* angleBendTerms[BENDDIM];
	double forceConstant, equilibriumAngle;
	double h3,h4,h5,h6;

	forceField-> numberOfBendTerms = 0;
	for( i=0; i<BENDDIM;i++) forceField->angleBendTerms[i] = NULL;
	if(!forceField->options.angleBend) return;

	numberOfBendTerms =  m->numberOf3Connections;
	for( i=0; i<BENDDIM;i++)
		angleBendTerms[i] =  malloc(numberOfBendTerms*sizeof(double)); 

	/* 5 terms 1=a1, 2=a2, 3=a3, 4=Force, 5=angle  h3 h4 h5 h6*/
	/* BENDDIM 9 */
	for ( i = 0; i < numberOfBendTerms; i++ )
	{
		a1 = m->connected3[0][i];
		a2 = m->connected3[1][i];
		a3 = m->connected3[2][i];
		a1Type = atomTypes[a1];
		a2Type = atomTypes[a2];
		a3Type = atomTypes[a3];

		if ( ! ( getBendParameters(amberParameters, a1Type, a2Type, a3Type,&forceConstant,&equilibriumAngle, &h3, &h4, &h5, &h6) ) )
		{
			/*
			char l1 = m->atoms[a1].mmType[0];
			char l2 = m->atoms[a2].mmType[0];
			char l3 = m->atoms[a3].mmType[0];
			*/
			fprintf(forceField->logfile, ("**** couldn't find bend parameters for %s-%s-%s "),
			m->atoms[a1].mmType,m->atoms[a2].mmType,m->atoms[a3].mmType);
			fprintf(forceField->logfile, "\n");
			fflush(forceField->logfile);
			/*
			forceConstant = 60.0;
			equilibriumAngle = 115.0;
			if(!strcmp(m->atoms[a2].mmType,"CT"))
			{
				forceConstant = 50.0;
				equilibriumAngle = 109.0;
			}
			else
			if(l1=='H' || l2=='H' || l3=='H')
			{
				forceConstant = 50.0;
				equilibriumAngle = 120.0;
			}
			if(isIonic( m->atoms[a1].mmType) || isIonic( m->atoms[a2].mmType) ||  isIonic( m->atoms[a3].mmType))
			{
				forceConstant = 0;
			}
			fprintf(forceField->logfile, ("-> I set force to %f and equilibrium angle to %f\n"), forceConstant, equilibriumAngle);
			fflush(forceField->logfile);
			*/
			angleBendTerms[0][i] = a1;
			angleBendTerms[1][i] = a2;
			angleBendTerms[2][i] = a3;
			angleBendTerms[3][i] = forceConstant;
			angleBendTerms[4][i] = equilibriumAngle;
			angleBendTerms[5][i] = h3;
			angleBendTerms[6][i] = h4;
			angleBendTerms[7][i] = h5;
			angleBendTerms[8][i] = h6;
		}
		else
		{
			angleBendTerms[0][i] = a1;
			angleBendTerms[1][i] = a2;
			angleBendTerms[2][i] = a3;
			angleBendTerms[3][i] = forceConstant;
			angleBendTerms[4][i] = equilibriumAngle;
			angleBendTerms[5][i] = h3;
			angleBendTerms[6][i] = h4;
			angleBendTerms[7][i] = h5;
			angleBendTerms[8][i] = h6;
		}
	}

	forceField-> numberOfBendTerms = numberOfBendTerms;
	for( i=0; i<BENDDIM;i++)
       		forceField->angleBendTerms[i] = angleBendTerms[i]; 

}
/**********************************************************************/
static void setStrBendParameters(AmberParameters* amberParameters,ForceField* forceField,int* atomTypes)
{
	int i;
	int a1,a2,a3;
	int a1Type, a2Type, a3Type;
	int type1, type2;
	Molecule* m = &forceField->molecule;
	int numberOfBendTerms = 0;
	int numberOfStrBendTerms = 0;
	double* strBendTerms[STRBENDDIM];
	double forceConstant, equilibriumAngle;
	double forceConstant12, forceConstant23;
	double alpha12, alpha23;
	double Re12, Re23;
	int nTerms = 0;
	double h3,h4,h5,h6;
	double h3OrDe1, h3OrDe2;

	if(!forceField->options.strBend) return;

	numberOfBendTerms =  m->numberOf3Connections;
	numberOfStrBendTerms =  m->numberOf3Connections;
	for( i=0; i<STRBENDDIM;i++) strBendTerms[i] =  malloc(numberOfStrBendTerms*sizeof(double)); 

	/* 8 terms 1=a1, 2=a2, 3=a3, 4=Force12, 5=Force23, 6=Re12 7=Re23 8=Angle0*/
	/* STRBENDDIM 8 */
	for ( i = 0; i < numberOfBendTerms; i++ )
	{
		a1 = m->connected3[0][i];
		a2 = m->connected3[1][i];
		a3 = m->connected3[2][i];
		a1Type = atomTypes[a1];
		a2Type = atomTypes[a2];
		a3Type = atomTypes[a3];

		if ( ( getBendParameters(amberParameters, a1Type, a2Type, a3Type,&forceConstant,&equilibriumAngle, &h3, &h4, &h5, &h6) ) )
		{
			if ( getStrBendParameters(amberParameters, a1Type, a2Type, a3Type,&forceConstant12,&forceConstant23)
			     &&
			     getStretchParameters(amberParameters, a1Type, a2Type,&alpha12,&Re12,&h3OrDe1,&h4,&h5,&h6,&type1)
			     && 
			     getStretchParameters(amberParameters, a2Type, a3Type,&alpha23,&Re23,&h3OrDe2,&h4,&h5,&h6,&type2) 
			)
			{
				strBendTerms[0][nTerms] = type1;
				strBendTerms[1][nTerms] = type2;
				strBendTerms[2][nTerms] = a1;
				strBendTerms[3][nTerms] = a2;
				strBendTerms[4][nTerms] = a3;
				strBendTerms[5][nTerms] = forceConstant12;
				strBendTerms[6][nTerms] = forceConstant23;
				strBendTerms[10][nTerms] = 0.0;
				strBendTerms[11][nTerms] = 0.0;
				/*
				if(type1==0) strBendTerms[5][nTerms] = forceConstant12;
				else{
					double De=h3OrDe1;
					strBendTerms[5][nTerms] = 0;
					if(fabs(De)>1e-13) strBendTerms[5][nTerms] = sqrt(forceConstant12/2.0/De);
					strBendTerms[10][nTerms] = De;
				}
				if(type2==0) strBendTerms[6][nTerms] = forceConstant23;
				else {
					double De=h3OrDe2;
					strBendTerms[6][nTerms] = 0;
					if(fabs(De)>1e-13) strBendTerms[6][nTerms] = sqrt(forceConstant23/2.0/De);
					strBendTerms[11][nTerms] = De;
				}
				*/
				strBendTerms[7][nTerms] = Re12;
				strBendTerms[8][nTerms] = Re23;
				strBendTerms[9][nTerms] = equilibriumAngle;
				nTerms++;
			}
		}
	}

	forceField->numberOfStrBendTerms = nTerms;
	for( i=0; i<STRBENDDIM;i++)
       		forceField->strBendTerms[i] = strBendTerms[i]; 
	//printf("numberOfStrBendTerms =%d\n",numberOfStrBendTerms);
}
/**********************************************************************/
static void setDihedralParameters(AmberParameters* amberParameters,ForceField* forceField,int* atomTypes)
{
	int i;
	int j;
	int k;
	int l;
	int a1,a2,a3,a4;
	int a1Type, a2Type, a3Type,a4Type;
	Molecule* m = &forceField->molecule;
	double* dihedralAngleTerms[DIHEDRALDIM];
	int numberOfDihedralTerms = 0;
	int dim;

	if(!forceField->options.dihedralAngle) return;
	/*  8 terms 1=a1, 2=a2, 3=a3, 4=a4, 5=Idiv, 6=Pk, 7=Phase, 8=Pn */
	/*  DIHEDRALDIM	8 */

	for( i=0; i<DIHEDRALDIM;i++)
		dihedralAngleTerms[i] =  malloc(4*m->numberOf4Connections*sizeof(double)); 

	numberOfDihedralTerms = 0;

	for (  i = 0; i < m->numberOf4Connections; i++ )
	{
		a1 = m->connected4[0][i];
		a2 = m->connected4[1][i];
		a3 = m->connected4[2][i];
		a4 = m->connected4[3][i];

		a1Type = atomTypes[a1];
		a2Type = atomTypes[a2];
		a3Type = atomTypes[a3];
		a4Type = atomTypes[a4];

		dim = getNumberDihedralParameters(amberParameters, a1Type, a2Type, a3Type, a4Type,&k);
		if(dim>0)
		{
			for(j=0;j<dim;j++)
			{
				dihedralAngleTerms[0][numberOfDihedralTerms] = a1;
				dihedralAngleTerms[1][numberOfDihedralTerms] = a2;
				dihedralAngleTerms[2][numberOfDihedralTerms] = a3;
				dihedralAngleTerms[3][numberOfDihedralTerms] = a4;
				dihedralAngleTerms[4][numberOfDihedralTerms] = 
					amberParameters->dihedralAngleTerms[k].divisor[j];
				dihedralAngleTerms[5][numberOfDihedralTerms] = 
					amberParameters->dihedralAngleTerms[k].barrier[j];
				dihedralAngleTerms[6][numberOfDihedralTerms] = 
					amberParameters->dihedralAngleTerms[k].phase[j];
				dihedralAngleTerms[7][numberOfDihedralTerms] = 
					amberParameters->dihedralAngleTerms[k].n[j];

				numberOfDihedralTerms++;
				if(numberOfDihedralTerms>4*m->numberOf4Connections)
				{
					for( l=0; l<DIHEDRALDIM;l++)
					{
						dihedralAngleTerms[l] =  
						realloc(dihedralAngleTerms[l],numberOfDihedralTerms*sizeof(double)); 
					}

				}
			}
		}		
	}

	forceField-> numberOfDihedralTerms = numberOfDihedralTerms;
	for( i=0; i<DIHEDRALDIM;i++)
       		forceField->dihedralAngleTerms[i] = dihedralAngleTerms[i]; 
}

/**********************************************************************/
static int getNumberImproperTorsionParameters( AmberParameters* amberParameters, int a1Type, int a2Type, int a3Type, int a4Type)
{
	int i;
	int a1Typet;
	int a2Typet;
	int a3Typet;
	int a4Typet;
	boolean btype;
	boolean Ok;
	int types[4];
	int k;

	a1Typet = a4Type;
	a2Typet = a3Type;
	a3Typet = a2Type;
	a4Typet = a1Type;

	/* Je cherche d'abord sans les -1 */
	for(i=0;i<amberParameters->numberOfImproperTorsionTerms;i++)
	{

		types[0] = a1Type;
		types[1] = a2Type;
		types[2] = a3Type;
		types[3] = a4Type;

		Ok = TRUE;
		for(k=0;k<4;k++)
		{
			btype = (types[k] == amberParameters->improperTorsionTerms[i].numbers[k]);
			if(!btype)
			{
				Ok = FALSE;
				break;
			}
		}
		if(!Ok)
		{
			types[0] = a1Typet;
			types[1] = a2Typet;
			types[2] = a3Typet;
			types[3] = a4Typet;
			Ok = TRUE;
			for(k=0;k<4;k++)
			{
				btype = (types[k] == amberParameters->improperTorsionTerms[i].numbers[k]);
				if(!btype)
				{
					Ok = FALSE;
					break;
				}
			}
		}

			 
		if(Ok)
		{
			return i;
		}
	}
	/* Je cherche d'abord avec les -1 */
	for(i=0;i<amberParameters->numberOfImproperTorsionTerms;i++)
	{

		types[0] = a1Type;
		types[1] = a2Type;
		types[2] = a3Type;
		types[3] = a4Type;

		Ok = TRUE;
		for(k=0;k<4;k++)
		{
			btype = 
			(amberParameters->improperTorsionTerms[i].numbers[k] == -1) || 
			(types[k] == amberParameters->improperTorsionTerms[i].numbers[k]);
			if(!btype)
			{
				Ok = FALSE;
				break;
			}
		}
		if(!Ok)
		{
			types[0] = a1Typet;
			types[1] = a2Typet;
			types[2] = a3Typet;
			types[3] = a4Typet;
			Ok = TRUE;
			for(k=0;k<4;k++)
			{
				btype = 
				(amberParameters->improperTorsionTerms[i].numbers[k] == -1) || 
				(types[k] == amberParameters->improperTorsionTerms[i].numbers[k]);
				if(!btype)
				{
					Ok = FALSE;
					break;
				}
			}
		}

			 
		if(Ok)
		{
			return i;
		}
	}

	return -1;
}
/**********************************************************************/
static boolean addImproperTorsionParameters(AmberParameters* amberParameters, int a1, int a2, int a3, int a4, int* atomTypes, double** improperTorsionTerms, int *pNumberOfImproperTorsionTerms)
{
	int i = *pNumberOfImproperTorsionTerms;
	
	int a1Type, a2Type, a3Type,a4Type;
	int n;
	a1Type = atomTypes[a1];
	a2Type = atomTypes[a2];
	a3Type = atomTypes[a3];
	a4Type = atomTypes[a4];
	printf("type a1 a2 a3 a4 = %d %d %d %d\n",a1Type, a2Type, a3Type, a4Type);
	n = getNumberImproperTorsionParameters(amberParameters, a1Type, a2Type, a3Type, a4Type);
	printf("n =%d\n",n);
	if(n>-1)
	{
		improperTorsionTerms[0][i] = a1;
		improperTorsionTerms[1][i] = a2;
		improperTorsionTerms[2][i] = a3;
		improperTorsionTerms[3][i] = a4;
		improperTorsionTerms[4][i] = amberParameters->improperTorsionTerms[n].barrier;
		improperTorsionTerms[5][i] = amberParameters->improperTorsionTerms[n].phase;
		improperTorsionTerms[6][i] = amberParameters->improperTorsionTerms[n].n;
		i++;
		*pNumberOfImproperTorsionTerms = i;
		return TRUE;
	}
	return FALSE;
}
/**********************************************************************/
static void setImproperTorionParameters(AmberParameters* amberParameters, ForceField* forceField,int* atomTypes)
{
	int i;
	int a1,a2,a3,a4;
	Molecule* m = &forceField->molecule;
	int numberOfImproperTorsionTerms = 0;
	double* improperTorsionTerms[IMPROPERDIHEDRALDIM];

	if(!forceField->options.improperTorsion) return;
	/*  8 terms 1=a1, 2=a2, 3=a3, 4=a4, 5=Idiv, 6=Pk, 7=Phase, 8=Pn */
	/*  IMPROPERDIHEDRALDIM	8 */

	for( i=0; i<IMPROPERDIHEDRALDIM;i++) improperTorsionTerms[i] =  malloc(6*m->numberOf3Connections*sizeof(double)); 

	numberOfImproperTorsionTerms = 0;

	for (  i = 0; i < m->numberOf3Connections; i++ )
	{
		int nc = 0;
		int j;
		int a4j;
		a1 = m->connected3[0][i];
		a2 = m->connected3[1][i];
		a3 = m->connected3[2][i];
		a4 = -1;
		
		for (  j = 0; j < m->numberOf2Connections; j++ ) if(m->connected2[0][j] == a2 || m->connected2[1][j] == a2) nc++;
		printf("a2 = %d nc = %d\n",a2,nc);
		if(nc!=3) continue;

		for (  j = 0; j < m->numberOf2Connections; j++ )
		{
			if(m->connected2[0][j] != a2 && m->connected2[1][j] != a2) continue;
			a4j = m->connected2[0][j]; 
			if(m->connected2[0][j] == a2) a4j = m->connected2[1][j];
			if(a4j == a1) continue;
			if(a4j == a3) continue;
			a4 = a4j;
			break;
		}
		printf("a1 a2 a3 a4 = %d %d %d %d\n",a1, a2, a3, a4);
		if(a4<0) continue;
		addImproperTorsionParameters(amberParameters, a1, a3, a2, a4, atomTypes, improperTorsionTerms, &numberOfImproperTorsionTerms);
		addImproperTorsionParameters(amberParameters, a3, a1, a2, a4, atomTypes, improperTorsionTerms, &numberOfImproperTorsionTerms);
		addImproperTorsionParameters(amberParameters, a1, a4, a2, a3, atomTypes, improperTorsionTerms, &numberOfImproperTorsionTerms);
		addImproperTorsionParameters(amberParameters, a4, a1, a2, a3, atomTypes, improperTorsionTerms, &numberOfImproperTorsionTerms);
		addImproperTorsionParameters(amberParameters, a3, a4, a2, a1, atomTypes, improperTorsionTerms, &numberOfImproperTorsionTerms);
		addImproperTorsionParameters(amberParameters, a4, a3, a2, a1, atomTypes, improperTorsionTerms, &numberOfImproperTorsionTerms);

	}

	forceField-> numberOfImproperTorsionTerms = numberOfImproperTorsionTerms;
	for( i=0; i<IMPROPERDIHEDRALDIM;i++) forceField->improperTorsionTerms[i] = improperTorsionTerms[i]; 

}
/**********************************************************************/
static void setHydrogenBondedParameters(AmberParameters* amberParameters,ForceField* forceField,int* atomTypes)
{
	int numberOfHydrogenBonded = 0;
	int i,j;
	int a1,a2,a4,ax;
	int a1Type,a2Type,a4Type;
	Molecule* m = &forceField->molecule;
	double C, D; // C = alpha, D = Re for Morse, C = A12 and D = A6 for 6-12 , C = A12 and D = A10 for 10-12
	double De;
	double* hydrogenBondedTerms[HYDROGENBONDEDDIM];
	if(!forceField->options.hydrogenBonded612 && !forceField->options.hydrogenBonded1012 && !forceField->options.hydrogenBondedMorse) return;


	for( i=0; i<HYDROGENBONDEDDIM;i++)
		hydrogenBondedTerms[i] =  malloc((m->numberOfNonBonded+m->numberOf4Connections)*sizeof(double));

	for ( i = 0; i < m->numberOfNonBonded; i++ )
	{
		boolean suttonChena1 = FALSE;
		boolean suttonChena2 = FALSE;
		a1 = m->nonBonded[0][i];
		a2 = m->nonBonded[1][i];
		ax = -1;

		a1Type = atomTypes[a1];
		a2Type = atomTypes[a2];
		suttonChena1 = canSuttonChen(forceField, a1);
		suttonChena2 = canSuttonChen(forceField, a2);
		if ( suttonChena1 || suttonChena2 ) continue;

		if ( canHydrogenBond( amberParameters, a1Type, a2Type ) )
		{ 
			De = 0.0;
			if(forceField->options.hydrogenBonded612) getHydrogenBonded612Parameters(amberParameters, a1Type, a2Type, &C, &D, forceField->logfile, m->atoms[a1].mmType, m->atoms[a2].mmType);
			else if(forceField->options.hydrogenBonded1012)  getHydrogenBonded1012Parameters(amberParameters, a1Type, a2Type, &C, &D );
			else getHydrogenBondedMorseParameters(amberParameters, a1Type, a2Type, &C, &D, &De );

			if(fabs(C)<1e-13 && fabs(D)<1e-13) continue;
			if(forceField->options.hbDirectional) 
			{
				int aH = a1;
				ax = -1;
				if(m->atoms[a2].mmType[0]=='H') aH = a2;
				for (  j = 0; j < m->numberOf2Connections; j++ ) 
				{
					if(m->connected2[0][j] == aH) { ax = m->connected2[1][j]; break;};
					if(m->connected2[1][j] == aH) { ax = m->connected2[0][j]; break;};
				}
			}
			hydrogenBondedTerms[0][numberOfHydrogenBonded] = a1;
			hydrogenBondedTerms[1][numberOfHydrogenBonded] = a2;
			hydrogenBondedTerms[2][numberOfHydrogenBonded] = ax;
			hydrogenBondedTerms[3][numberOfHydrogenBonded] = C;
			hydrogenBondedTerms[4][numberOfHydrogenBonded] = D;
			hydrogenBondedTerms[5][numberOfHydrogenBonded] = De;
			numberOfHydrogenBonded++;
		}
	}
	/* now 1/2 non bonded */
	for (  i = 0; i < m->numberOf4Connections; i++ )
	{
		boolean suttonChena1 = FALSE;
		boolean suttonChena4 = FALSE;

		a1 = m->connected4[0][i];
		a4 = m->connected4[3][i];
		ax = -1;
		if(forceField->options.hbDirectional) 
		{
			if(m->atoms[a1].mmType[0] =='H') ax =  m->connected4[1][i];
			else ax =  m->connected4[2][i];
		}

		a1Type = atomTypes[a1];
		a4Type = atomTypes[a4];
		suttonChena1 = canSuttonChen(forceField, a1);
		suttonChena4 = canSuttonChen(forceField, a4);
		if ( suttonChena1 || suttonChena4 ) continue;
		if ( canHydrogenBond( amberParameters, a1Type, a4Type ) )
		{ 
			De = 0;
			if(forceField->options.hydrogenBonded612) 
			{
				getHydrogenBonded612Parameters(amberParameters, a1Type, a4Type, &C, &D, forceField->logfile, m->atoms[a1].mmType, m->atoms[a4].mmType);
				C /= 2.0;
				D /= 2.0;
			}
			else if(forceField->options.hydrogenBonded1012)  
			{
				getHydrogenBonded1012Parameters(amberParameters, a1Type, a4Type, &C, &D );
				C /= 2.0;
				D /= 2.0;
			}
			else 
			{
				getHydrogenBondedMorseParameters(amberParameters, a1Type, a4Type, &C, &D, &De );
				De /= 2;
			}
			if(fabs(C)<1e-13 && fabs(D)<1e-13) continue;
			hydrogenBondedTerms[0][numberOfHydrogenBonded] = a1;
			hydrogenBondedTerms[1][numberOfHydrogenBonded] = a4;
			hydrogenBondedTerms[2][numberOfHydrogenBonded] = ax;
			hydrogenBondedTerms[3][numberOfHydrogenBonded] = C;
			hydrogenBondedTerms[4][numberOfHydrogenBonded] = D;
			hydrogenBondedTerms[5][numberOfHydrogenBonded] = De;
			numberOfHydrogenBonded++;
		}
	}

	if(numberOfHydrogenBonded==0)
		for( i=0; i<HYDROGENBONDEDDIM;i++)
		{
			free(hydrogenBondedTerms[i]);
			hydrogenBondedTerms[i] = NULL;
		}
	else
		for( i=0; i<HYDROGENBONDEDDIM;i++)
		{
			hydrogenBondedTerms[i] = 
				realloc(hydrogenBondedTerms[i],numberOfHydrogenBonded*sizeof(double));
		}

	forceField-> numberOfHydrogenBonded = numberOfHydrogenBonded;
	for( i=0; i<HYDROGENBONDEDDIM;i++)
       		forceField->hydrogenBondedTerms[i] = hydrogenBondedTerms[i]; 
	printf("numberOfHydrogenBonded=%d\n",numberOfHydrogenBonded);

}
/**********************************************************************/
static void setVdw612Parameters(AmberParameters* amberParameters, ForceField* forceField,int* atomTypes)
{
	int numberOfVdw612 = 0;
	int i;
	int a1,a2,a4;
	int a1Type,a2Type,a4Type;
	Molecule* m = &forceField->molecule;
	double equilibriumDistance, epsilon;
	double epsilonProduct;
	double ri, rj;
	double Aij, Bij;
	double* vdw612Terms[VDW612DIM];
	boolean useHydrogenBonded = forceField->options.hydrogenBonded612 || forceField->options.hydrogenBonded1012 || forceField->options.hydrogenBondedMorse;

	if(!forceField->options.vdw612) return;

	/* 5 terms 1=a1, 2=a2, 3=Aij, 4=Bij*/
	/* VDW612DIM 4 */
	for( i=0; i<VDW612DIM;i++) vdw612Terms[i] =  malloc((m->numberOfNonBonded+m->numberOf4Connections)*sizeof(double)); 

	for ( i = 0; i < m->numberOfNonBonded; i++ )
	{
		boolean suttonChena1 = FALSE;
		boolean suttonChena2 = FALSE;

		a1 = m->nonBonded[0][i];
		a2 = m->nonBonded[1][i];
		a1Type = atomTypes[a1];
		a2Type = atomTypes[a2];

		suttonChena1 = canSuttonChen(forceField, a1);
		suttonChena2 = canSuttonChen(forceField, a2);
		if ( suttonChena1 || suttonChena2 ) continue;

		if ( !useHydrogenBonded || !canHydrogenBond(amberParameters, a1Type, a2Type ) )
		{ 
			if ( ! ( getVdw612Parameters(amberParameters, a1Type, &equilibriumDistance, &epsilon ) ) )
			{
				fprintf(forceField->logfile, ("**** couldn't find non bonded parameters for %s \n"),m->atoms[a1].mmType);
				fflush(forceField->logfile);
			}
		
			epsilonProduct = sqrt(fabs(epsilon));
			ri = equilibriumDistance;
			/*printf("r1 = %f eps1 = %f\n",equilibriumDistance,epsilon);*/

			getVdw612Parameters(amberParameters, a2Type, &equilibriumDistance, &epsilon );
			/*printf("r2 = %f eps2 = %f\n",equilibriumDistance,epsilon);*/
			epsilonProduct *= sqrt(fabs(epsilon));
			rj = equilibriumDistance;
			Bij = ( ri + rj ) * ( ri + rj );
			Bij = Bij * Bij * Bij;
			Aij = Bij * Bij * epsilonProduct;
			Bij *= epsilonProduct * 2.0;
			if(fabs(Aij)<1e-13 && fabs(Bij)<1e-13) continue;
			/*
			if(a1Type==64 || a2Type==64)
			{
				printf("Aij = %f\n",Aij);
				printf("Bij = %f\n",Bij);
			}
			*/

			vdw612Terms[0][numberOfVdw612] = a1;
			vdw612Terms[1][numberOfVdw612] = a2;
			vdw612Terms[2][numberOfVdw612] = Aij;
			vdw612Terms[3][numberOfVdw612] = Bij;
			numberOfVdw612++;
		}
	}

	/* now 1/2 non bonded */
	for (  i = 0; i < m->numberOf4Connections; i++ )
	{
		boolean suttonChena1 = FALSE;
		boolean suttonChena4 = FALSE;

		a1 = m->connected4[0][i];
		a4 = m->connected4[3][i];
		a1Type = atomTypes[a1];
		a4Type = atomTypes[a4];
		suttonChena1 = canSuttonChen(forceField, a1);
		suttonChena4 = canSuttonChen(forceField, a4);
		if ( suttonChena1 || suttonChena4 ) continue;
		epsilonProduct = 0;
		ri = 0;
		rj = 0;
	        if ( getVdw612Parameters(amberParameters, a1Type, &equilibriumDistance, &epsilon ) )
		{
	        	epsilonProduct = sqrt(fabs(epsilon));
	        	ri = equilibriumDistance;
			/*printf("r1 = %f eps1 = %f\n",equilibriumDistance,epsilon);*/
		}
		else
		{
			epsilonProduct = 0;
		}

	        if ( getVdw612Parameters( amberParameters, a4Type, &equilibriumDistance, &epsilon ) )
		{
	        	epsilonProduct *= sqrt(fabs(epsilon));
	        	rj = equilibriumDistance;
			/*printf("r2 = %f eps2 = %f\n",equilibriumDistance,epsilon);*/
		}
		else
		{
			epsilonProduct = 0;
		}

	       	Bij = ( ri + rj ) * ( ri + rj );
	       	Bij = Bij * Bij * Bij;
	       	Aij = Bij * Bij * epsilonProduct / 2.0;
	       	Bij *= epsilonProduct;
		if(fabs(Aij)<1e-13 && fabs(Bij)<1e-13) continue;

		/*
			Aij = 0;
			Bij = 0;
		*/

		vdw612Terms[0][numberOfVdw612] = a1;
		vdw612Terms[1][numberOfVdw612] = a4;
		vdw612Terms[2][numberOfVdw612] = Aij;
		vdw612Terms[3][numberOfVdw612] = Bij;
		numberOfVdw612++;
	}
	if(numberOfVdw612==0)
		for( i=0; i<VDW612DIM;i++)
		{
			free(vdw612Terms[i]);
			vdw612Terms[i] = NULL;
		}
	else
		for( i=0; i<VDW612DIM;i++)
			vdw612Terms[i] = 
				realloc(vdw612Terms[i],numberOfVdw612*sizeof(double));

	forceField-> numberOfVdw612 = numberOfVdw612;
	for( i=0; i<VDW612DIM;i++) forceField->vdw612Terms[i] = vdw612Terms[i]; 
	//printf("numberOfVdw612=%d\n",numberOfVdw612);
}
/**********************************************************************/
static void setVdw714Parameters(AmberParameters* amberParameters, ForceField* forceField,int* atomTypes)
{
	int numberOfVdw714 = 0;
	int i;
	int a1,a2,a4;
	int a1Type,a2Type,a4Type;
	Molecule* m = &forceField->molecule;
	double R01, epsilon1;
	double R02, epsilon2;
	double gamma1, delta1;
	double gamma2, delta2;
	double gamma, delta;
	double epsilonij;
	double Rij0;
	double* vdw714Terms[VDW714DIM];
	boolean useHydrogenBonded = forceField->options.hydrogenBonded612 || forceField->options.hydrogenBonded1012 || forceField->options.hydrogenBondedMorse;
	double term;

	if(!forceField->options.vdw714) return;
	/* 5 terms 1=a1, 2=a2, 3=Aij, 4=Bij*/
	/* VDW714DIM 4 */
	for( i=0; i<VDW714DIM;i++) vdw714Terms[i] =  malloc((m->numberOfNonBonded+m->numberOf4Connections)*sizeof(double)); 

	for ( i = 0; i < m->numberOfNonBonded; i++ )
	{
		boolean suttonChena1 = FALSE;
		boolean suttonChena2 = FALSE;

		a1 = m->nonBonded[0][i];
		a2 = m->nonBonded[1][i];
		a1Type = atomTypes[a1];
		a2Type = atomTypes[a2];

		suttonChena1 = canSuttonChen(forceField, a1);
		suttonChena2 = canSuttonChen(forceField, a2);
		if ( suttonChena1 || suttonChena2 ) continue;

		if ( !useHydrogenBonded || !canHydrogenBond(amberParameters, a1Type, a2Type ) )
		{ 
			if ( ! ( getVdw714Parameters(amberParameters, a1Type, &R01, &epsilon1, &gamma1, &delta1 ) ) ) 
			{
				fprintf(forceField->logfile, ("**** couldn't find non bonded parameters for %s \n"),m->atoms[a1].mmType);
				fflush(forceField->logfile);
				continue;
			}
			if(!getVdw714Parameters(amberParameters, a2Type, &R02, &epsilon2, &gamma2, &delta2 ))
			{
				fprintf(forceField->logfile, ("**** couldn't find non bonded parameters for %s \n"),m->atoms[a2].mmType);
				fflush(forceField->logfile);
				continue;
			}

			term = sqrt(fabs(epsilon1))+sqrt(fabs(epsilon2));
			epsilonij = 4*( epsilon1*epsilon2)/(term*term);
			if(fabs(epsilonij)<1e-13) continue;
			Rij0 = (R01*R01*R01+R02*R02*R02)/(R01*R01+R02*R02);
			if(fabs(Rij0)<1e-13) continue;
			gamma = gamma1 + gamma2;
			delta = delta1 + delta2;

			vdw714Terms[0][numberOfVdw714] = a1;
			vdw714Terms[1][numberOfVdw714] = a2;
			vdw714Terms[2][numberOfVdw714] = epsilonij;
			vdw714Terms[3][numberOfVdw714] = Rij0;
			vdw714Terms[4][numberOfVdw714] = gamma;
			vdw714Terms[5][numberOfVdw714] = delta;
			numberOfVdw714++;
		}
	}

	/* now 1/2 non bonded */
	for (  i = 0; i < m->numberOf4Connections; i++ )
	{
		boolean suttonChena1 = FALSE;
		boolean suttonChena4 = FALSE;

		a1 = m->connected4[0][i];
		a4 = m->connected4[3][i];
		a1Type = atomTypes[a1];
		a4Type = atomTypes[a4];
		suttonChena1 = canSuttonChen(forceField, a1);
		suttonChena4 = canSuttonChen(forceField, a4);
		if ( suttonChena1 || suttonChena4 ) continue;
	        if (! getVdw714Parameters(amberParameters, a1Type, &R01, &epsilon1, &gamma1, &delta1 ) )
		{
			fprintf(forceField->logfile, ("**** couldn't find non bonded parameters for %s \n"),m->atoms[a1].mmType);
			fflush(forceField->logfile);
			continue;
		}
	        if (! getVdw714Parameters(amberParameters, a4Type, &R02, &epsilon2, &gamma2, &delta2 ) )
		{
			fprintf(forceField->logfile, ("**** couldn't find non bonded parameters for %s \n"),m->atoms[a4].mmType);
			fflush(forceField->logfile);
			continue;
		}

		term = sqrt(fabs(epsilon1))+sqrt(fabs(epsilon2));
		epsilonij = 4*( epsilon1*epsilon2)/(term*term);
		if(fabs(epsilonij)<1e-13) continue;
		Rij0 = (R01*R01*R01+R02*R02*R02)/(R01*R01+R02*R02);
		if(fabs(Rij0)<1e-13) continue;
		gamma = gamma1 + gamma2;
		delta = delta1 + delta2;

		vdw714Terms[0][numberOfVdw714] = a1;
		vdw714Terms[1][numberOfVdw714] = a4;
		vdw714Terms[2][numberOfVdw714] = epsilonij;
		vdw714Terms[3][numberOfVdw714] = Rij0;
		vdw714Terms[4][numberOfVdw714] = gamma;
		vdw714Terms[5][numberOfVdw714] = delta;
		numberOfVdw714++;
	}
	if(numberOfVdw714==0)
	for( i=0; i<VDW714DIM;i++)
	{
		free(vdw714Terms[i]);
		vdw714Terms[i] = NULL;
	}
	else
	for( i=0; i<VDW714DIM;i++) vdw714Terms[i] = realloc(vdw714Terms[i],numberOfVdw714*sizeof(double));

	forceField-> numberOfVdw714 = numberOfVdw714;
	for( i=0; i<VDW714DIM;i++) forceField->vdw714Terms[i] = vdw714Terms[i]; 
	printf("numberOfVdw714=%d\n",numberOfVdw714);
}
/**********************************************************************/
static void setCoulombParameters(AmberParameters* amberParameters, ForceField* forceField,int* atomTypes)
{
	int numberOfCoulomb = 0;
	int i;
	int a1,a2,a4;
	Molecule* m = &forceField->molecule;
	double* coulombTerms[COULOMBDIM];
	double permittivityScale = 1, permittivity = 1;
	double coulombFactor = 332.05382 / ( permittivity * permittivityScale );
	boolean useCoulomb = forceField->options.coulomb;

	if(!forceField->options.coulomb) return;

	/* 5 terms 1=a1, 2=a2, 3=Aij, 4=Bij*/
	/* COULOMBDIM 4 */
	for( i=0; i<COULOMBDIM;i++) coulombTerms[i] =  malloc((m->numberOfNonBonded+m->numberOf4Connections)*sizeof(double)); 

	if(useCoulomb)
	for ( i = 0; i < m->numberOfNonBonded; i++ )
	{
		boolean suttonChena1 = FALSE;
		boolean suttonChena2 = FALSE;

		a1 = m->nonBonded[0][i];
		a2 = m->nonBonded[1][i];

		suttonChena1 = canSuttonChen(forceField, a1);
		suttonChena2 = canSuttonChen(forceField, a2);
		if ( suttonChena1 || suttonChena2 ) continue;
		if(fabs(m->atoms[a1].charge*m->atoms[a2].charge)<1e-13) continue;

		coulombTerms[0][numberOfCoulomb] = a1;
		coulombTerms[1][numberOfCoulomb] = a2;
		//coulombTerms[2][numberOfCoulomb] = 1.0*m->atoms[a1].charge*m->atoms[a2].charge*coulombFactor;
		coulombTerms[2][numberOfCoulomb] = 1.0*coulombFactor;
		//coulombTerms[2][numberOfCoulomb] = 0.57*coulombFactor;
		numberOfCoulomb++;
	}

	/* now 1/2 non bonded */
	if(useCoulomb)
	for (  i = 0; i < m->numberOf4Connections; i++ )
	{
		boolean suttonChena1 = FALSE;
		boolean suttonChena4 = FALSE;

		a1 = m->connected4[0][i];
		a4 = m->connected4[3][i];
		suttonChena1 = canSuttonChen(forceField, a1);
		suttonChena4 = canSuttonChen(forceField, a4);
		if ( suttonChena1 || suttonChena4 ) continue;
		if(fabs(m->atoms[a1].charge*m->atoms[a4].charge)<1e-13) continue;
		coulombTerms[0][numberOfCoulomb] = a1;
		coulombTerms[1][numberOfCoulomb] = a4;
		//coulombTerms[2][numberOfCoulomb] = 1.0/(double)1.2*m->atoms[a1].charge*m->atoms[a4].charge*coulombFactor;
		//coulombTerms[2][numberOfCoulomb] = 1.0/(double)1.2*coulombFactor;
		coulombTerms[2][numberOfCoulomb] = 0.57*coulombFactor;
		numberOfCoulomb++;
	}
	if(numberOfCoulomb==0)
		for( i=0; i<COULOMBDIM;i++)
		{
			free(coulombTerms[i]);
			coulombTerms[i] = NULL;
		}
	else
		for( i=0; i<COULOMBDIM;i++) coulombTerms[i] = realloc(coulombTerms[i],numberOfCoulomb*sizeof(double));

	forceField-> numberOfCoulomb = numberOfCoulomb;
	for( i=0; i<COULOMBDIM;i++) forceField->coulombTerms[i] = coulombTerms[i]; 
	//printf("numberOfCoulomb=%d\n",numberOfCoulomb);
}
/**********************************************************************/
static void setSuttonChenParameters(AmberParameters* amberParameters, ForceField* forceField,int* atomTypes)
{
	int numberOfSuttonChen = 0;
	int i,j;
	int a1,a2;
	int a1Type,a2Type;
	Molecule* mol = &forceField->molecule;
	double epsilon, a, C, n, m;
	double* suttonChenTerms[SUTTONCHENDIM];
	int nA;
	if(!forceField->options.suttonChen) return;

	/* 7 terms 1=a1, 2=a2, 3=epsilon, 4=a, 5=C , 6=n, 7=m */
	/* SUTTONCHENDIM 7 */
	nA = mol->nAtoms*( mol->nAtoms-1)/2;
	for( i=0; i<SUTTONCHENDIM;i++) suttonChenTerms[i] =  malloc(nA*sizeof(double)); 

	for ( i = 0; i < mol->nAtoms; i++ )
	{
		a1 = i;
		for ( j = i+1; j < mol->nAtoms; j++ )
		{
			a2 = j;

			a1Type = atomTypes[a1];
			a2Type = atomTypes[a2];

			if (  ( getSuttonChenParameters(amberParameters, a1Type, a2Type, &epsilon, &a,&C,&n,&m ) ) )
			{
				suttonChenTerms[0][numberOfSuttonChen] = a1;
				suttonChenTerms[1][numberOfSuttonChen] = a2;
				suttonChenTerms[2][numberOfSuttonChen] = epsilon;
				suttonChenTerms[3][numberOfSuttonChen] = a;
				suttonChenTerms[4][numberOfSuttonChen] = C;
				suttonChenTerms[5][numberOfSuttonChen] = n;
				suttonChenTerms[6][numberOfSuttonChen] = m;
				numberOfSuttonChen++;
			}
                }
	}

	if(numberOfSuttonChen==0)
		for( i=0; i<SUTTONCHENDIM;i++)
		{
			free(suttonChenTerms[i]);
			suttonChenTerms[i] = NULL;
		}
	else
		for( i=0; i<SUTTONCHENDIM;i++)
			suttonChenTerms[i] = 
				realloc(suttonChenTerms[i],numberOfSuttonChen*sizeof(double));

	forceField-> numberOfSuttonChen = numberOfSuttonChen;
	for( i=0; i<SUTTONCHENDIM;i++)
       		forceField->suttonChenTerms[i] = suttonChenTerms[i]; 
}
/**********************************************************************/
static void setPairWiseParameters(AmberParameters* amberParameters, ForceField* forceField,int* atomTypes)
{
	int numberOfPairWise = 0;
	int i;
	int j;
	int a1,a2;
	int a1Type,a2Type;
	Molecule* m = &forceField->molecule;
	double a, beta, c4, c6, c8, c10, b;
	double* pairWiseTerms[PAIRWISEDIM];

	numberOfPairWise = m->nAtoms*(m->nAtoms-1)/2;

	/* PAIRWISEDIM 8 */
	for( i=0; i<PAIRWISEDIM;i++)
		pairWiseTerms[i] =  
			malloc((numberOfPairWise)*sizeof(double)); 

	numberOfPairWise = 0;
	for ( i = 0; i < m->nAtoms; i++ )
	for ( j = i+1; j < m->nAtoms; j++ )
	{
		a1 = i;
		a2 = j;

		a1Type = atomTypes[a1];
		a2Type = atomTypes[a2];

		if ( ! ( getPairWiseParameters(amberParameters, a1Type,a2Type,&a, &beta,&c4, &c6,&c8, &c10,&b) ) )
		{
			fprintf(forceField->logfile,  ("**** couldn't find pair wise parameters for %s-%s\n"),
					m->atoms[a1].mmType, m->atoms[a2].mmType);
			fflush(forceField->logfile);
		}
		
			pairWiseTerms[0][numberOfPairWise] = a1;
			pairWiseTerms[1][numberOfPairWise] = a2;
			pairWiseTerms[2][numberOfPairWise] = a;
			pairWiseTerms[3][numberOfPairWise] = beta;
			pairWiseTerms[4][numberOfPairWise] = c4;
			pairWiseTerms[5][numberOfPairWise] = c6;
			pairWiseTerms[6][numberOfPairWise] = c8;
			pairWiseTerms[7][numberOfPairWise] = c10;
			pairWiseTerms[8][numberOfPairWise] = b;
			numberOfPairWise++;
	}

	if(numberOfPairWise==0)
		for( i=0; i<PAIRWISEDIM;i++)
		{
			free(pairWiseTerms[i]);
			pairWiseTerms[i] = NULL;
		}
	else
		for( i=0; i<PAIRWISEDIM;i++)
			pairWiseTerms[i] = 
				realloc(pairWiseTerms[i],numberOfPairWise*sizeof(double));

	forceField-> numberOfPairWise = numberOfPairWise;
	for( i=0; i<PAIRWISEDIM;i++)
       		forceField->pairWiseTerms[i] = pairWiseTerms[i]; 
}
/**********************************************************************/
static void setACKS2Parameters(AmberParameters* amberParameters, ForceField* forceField)
{

	/*
	 if(
	strcmp(forceField->options.chargesType,"ACKS2")
	&& strcmp(forceField->options.chargesType,"EEM")
	) return;
	*/
        int nTypes = amberParameters->numberOfTypes;
        AmberAtomTypes* types = amberParameters->atomTypes;
	char** mmTypes1 = NULL;
	char** mmTypes2 = NULL;
	double* dum = NULL;
	Molecule* molecule = &forceField->molecule;
	int nBonds = 0;
	int i,j;

	if(nTypes<1) return;

	mmTypes1 = malloc(nTypes*sizeof(char*));
	dum = malloc(nTypes*sizeof(double));

        for(i=0;i<nTypes;i++) mmTypes1[i] = strdup(types[i].name);
        for(i=0;i<nTypes;i++) dum[i] = types[i].electronegativity;
        molecule->klass->setElectronegativity(molecule, nTypes, mmTypes1, dum);
        for(i=0;i<nTypes;i++) dum[i] = types[i].hardness;
        molecule->klass->setHardness(molecule, nTypes, mmTypes1, dum);
        for(i=0;i<nTypes;i++) dum[i] = types[i].width;
        molecule->klass->setWidth(molecule, nTypes, mmTypes1, dum);
        for(i=0;i<nTypes;i++) dum[i] = types[i].charge;
        molecule->klass->setCharge0(molecule, nTypes, mmTypes1, dum);

        for(i=0;i<nTypes;i++) if( mmTypes1 && mmTypes1[i]) free(mmTypes1[i]);
        if( mmTypes1) free(mmTypes1);
	if(dum) free(dum);
	
	nBonds = amberParameters->numberOfHardnessTerms;
	if(nBonds<1) return;
	dum = malloc(nBonds*sizeof(double));
	mmTypes1 = malloc(nBonds*sizeof(char*));
	mmTypes2 = malloc(nBonds*sizeof(char*));
        for(i=0;i<nBonds;i++) mmTypes1[i] = NULL;
        for(i=0;i<nBonds;i++) mmTypes2[i] = NULL;
        for(i=0;i<nBonds;i++) 
	{
		dum[i] = amberParameters->bondHardnessTerms[i].kappa;
		for(j=0;j<nTypes;j++)
        	{
                        if(amberParameters->bondHardnessTerms[i].numbers[0]==types[j].number)  mmTypes1[i] = strdup(types[j].name);
			if(amberParameters->bondHardnessTerms[i].numbers[1]==types[j].number)  mmTypes2[i] = strdup(types[j].name);
        	}
		if(mmTypes1[i] == NULL) mmTypes1[i] = strdup("X");
		if(mmTypes2[i] == NULL) mmTypes2[i] = strdup("X");
	}
        molecule->klass->setBondHardness(molecule, nBonds, mmTypes1, mmTypes2, dum);
        for(i=0;i<nBonds;i++) if( mmTypes1 && mmTypes1[i]) free(mmTypes1[i]);
        if( mmTypes1) free(mmTypes1);
        for(i=0;i<nBonds;i++) if( mmTypes2 && mmTypes2[i]) free(mmTypes2[i]);
        if( mmTypes2) free(mmTypes2);
	if(dum) free(dum);
}
/**********************************************************************/
/**********************************************************************/
static void setAtomTypes(AmberParameters* amberParameters,ForceField* forceField,int* atomTypes)
{
	Molecule* m = &forceField->molecule;
	int nAtoms = m->nAtoms;
	int i;
	for(i=0;i<nAtoms;i++)
	{ 
		/* printf("Atom %s=",m->atoms[i].mmType); */
		atomTypes[i] = getNumberType(amberParameters, m->atoms[i].mmType);
		/*
		{
			int j;
			int nTypes = amberParameters->numberOfTypes;
			AmberAtomTypes* types = amberParameters->atomTypes;
			char* type = m->atoms[i].mmType;
			int len = strlen(type);

			if(strcmp(type,"X")==0)
				printf("-1\n");
			for(j=0;j<nTypes;j++)
			{
					if(len == (int)strlen(types[j].name) && 
						strstr(types[j].name,type))
						printf(" %d \n",types[j].number);
			}
		}
		*/
	}
	
}
/**********************************************************************/
static char* getFileNameParameters()
{
	FILE* file;
	/* check if I can read parameters from a file in Default folder */
	char* fileName = strdup_printf("MolecularMechanics.prm");
	file = fopen(fileName,"r");
	if(file)
	{
		fclose(file);
		return fileName;
	}
	free(fileName);
	return strdup_printf("%s%sMolecularMechanics.prm",cchemiDirectory(), DIR_SEPARATOR_S);
}
/**********************************************************************/
static void setAmberParameters(ForceField* forceField)
{
	Molecule* m = &forceField->molecule;
	int* atomTypes = malloc(m->nAtoms*sizeof(int));
	AmberParameters amberParameters;

	if(staticAmberParameters && staticAmberParameters->numberOfTypes >0 )
		amberParameters = *staticAmberParameters;
	else
	{
		char* defaultFileName = getFileNameParameters();

		amberParameters =  newAmberParameters();
		if(!readAmberParameters(&amberParameters,defaultFileName))
		{
			free(defaultFileName);
			return;
		}

		staticAmberParameters = malloc(sizeof(AmberParameters));
		*staticAmberParameters = amberParameters;

		free(defaultFileName);

	}

	setAtomTypes(&amberParameters,forceField,atomTypes);
	setACKS2Parameters(&amberParameters, forceField);
	setChargesFromMMParameters(&amberParameters,forceField);
	setCharges(forceField);
	setStretchParameters(&amberParameters,forceField,atomTypes);
	setBendParameters(&amberParameters,forceField,atomTypes);
	setStrBendParameters(&amberParameters,forceField,atomTypes);
	setDihedralParameters(&amberParameters, forceField,atomTypes);
	setImproperTorionParameters(&amberParameters,forceField,atomTypes);
	//printf("options.outOfPlane=%d\n",forceField->options.outOfPlane);
	setOutOfPlaneParameters(&amberParameters,forceField,atomTypes);
	//printf("forceField->options.hydrogenBonded=%d\n",forceField->options.hydrogenBonded);
	setHydrogenBondedParameters(&amberParameters,forceField,atomTypes);
	setSuttonChenParameters(&amberParameters,forceField,atomTypes);
	setVdw612Parameters(&amberParameters,forceField,atomTypes);
	setVdw714Parameters(&amberParameters,forceField,atomTypes);
	//printf("forceField->options.coulomb=%d\n",forceField->options.coulomb);
	setCoulombParameters(&amberParameters,forceField,atomTypes);
	setRattleConstraintsParameters(forceField);
	/*
	freeAmberParameters(&amberParameters);
	*/
}
/**********************************************************************/
static void setAllPairWiseParameters(ForceField* forceField)
{
	Molecule* m = &forceField->molecule;
	int* atomTypes = malloc(m->nAtoms*sizeof(int));
	AmberParameters amberParameters;



	if(staticAmberParameters && staticAmberParameters->numberOfTypes >0 )
		amberParameters = *staticAmberParameters;
	else
	{
		char* defaultFileName = getFileNameParameters();

		amberParameters =  newAmberParameters();
		if(!readAmberParameters(&amberParameters,defaultFileName))
		{
			free(defaultFileName);
			return;
		}

		staticAmberParameters = malloc(sizeof(AmberParameters));
		*staticAmberParameters = amberParameters;

		free(defaultFileName);

	}
	

	setAtomTypes(&amberParameters,forceField,atomTypes);
	
	setPairWiseParameters(&amberParameters,forceField,atomTypes);

	setRattleConstraintsParameters(forceField);

	/*
	freeAmberParameters(&amberParameters);
	*/
}
#ifndef ENABLE_CL
/**********************************************************************/
static void calculateGradientBondAmber(ForceField* forceField)
{
	int i;
	Molecule* m = &forceField->molecule;
	double* bondStretchTerms[STRETCHDIM];
	int numberOfStretchTerms = forceField->numberOfStretchTerms;
	double energy = 0;

	for( i=0; i<STRETCHDIM;i++)
       		bondStretchTerms[i] = forceField->bondStretchTerms[i];

#ifdef ENABLE_OMP
#pragma omp parallel for private(i) 
#endif
	for ( i = 0; i < numberOfStretchTerms; i++ )
	{
		int ai, aj;
		double rijx, rijy, rijz, forceOrAlpha, equilibriumDistance, term;
		double h3OrDe,h4,h5,h6;
		double gradix, gradiy, gradiz;
		double bondLength;
		double diff = 0;
		double diff2 = 0;
		double diff3 = 0;
		double diff4 = 0;
		double diff5 = 0;
		double diff6 = 0;
		double der = 0;
		int type;

		type = bondStretchTerms[0][i];
		ai = (int)bondStretchTerms[1][i];
		aj = (int)bondStretchTerms[2][i];
		forceOrAlpha = bondStretchTerms[3][i];
		equilibriumDistance = bondStretchTerms[4][i];
		h3OrDe = bondStretchTerms[5][i];
		h4 = bondStretchTerms[6][i];
		h5 = bondStretchTerms[7][i];
		h6 = bondStretchTerms[8][i];

		rijx = m->atoms[ai].coordinates[0] - m->atoms[aj].coordinates[0];
		rijy = m->atoms[ai].coordinates[1] - m->atoms[aj].coordinates[1];
		rijz = m->atoms[ai].coordinates[2] - m->atoms[aj].coordinates[2];

		bondLength = sqrt( rijx * rijx + rijy * rijy + rijz * rijz );

		if ( bondLength < POSEPS ) bondLength = POSEPS;

		if(type==0)
		{
			diff = bondLength - equilibriumDistance;
			diff2 = diff*diff;
			diff3 = diff2*diff;
			diff4 = diff3*diff;
			diff5 = diff4*diff;
			diff6 = diff5*diff;
			der = forceOrAlpha/ bondLength;

			term = 2*der*diff;
			term += 3*h3OrDe*der*diff2;
			term += 4*h4*der*diff3;
			term += 5*h5*der*diff4;
			term += 6*h6*der*diff5;
		}
		else
		{
			double alpha = forceOrAlpha;
			double De = h3OrDe;
			double X = exp(-alpha*(bondLength-equilibriumDistance));
			double t = 1-X;
			term = 2*De*t*X*alpha/bondLength;
		}
		gradix = term * rijx;
		gradiy = term * rijy;
		gradiz = term * rijz;
#ifdef ENABLE_OMP
#pragma omp critical
#endif
		{
			m->atoms[ai].gradient[0] += gradix;
			m->atoms[ai].gradient[1] += gradiy;
			m->atoms[ai].gradient[2] += gradiz;
		
			m->atoms[aj].gradient[0] -= gradix;
			m->atoms[aj].gradient[1] -= gradiy;
			m->atoms[aj].gradient[2] -= gradiz;
		}
		if(type==0)
		{
			energy += forceOrAlpha * diff2;
			energy += h3OrDe * forceOrAlpha * diff3;
			energy += h4 * forceOrAlpha * diff4;
			energy += h5 * forceOrAlpha * diff5;
			energy += h6 * forceOrAlpha * diff6;
		}
		else
		{
			double alpha = forceOrAlpha;
			double De = h3OrDe;
			double X = exp(-alpha*(bondLength-equilibriumDistance));
			double t = 1-X;
			energy += De*t*t;
		}
	} 
	m->potentialEnergy += energy;
}
/**********************************************************************/
static void calculateGradientBendAmber(ForceField* forceField)
{
	int i;

	Molecule* m = &forceField->molecule;
	double* angleBendTerms[BENDDIM];
	static double D2R = 1/(RADTODEG);
	int numberOfBendTerms = forceField->numberOfBendTerms;
	double energy = 0;

	for( i=0; i<BENDDIM;i++)
		angleBendTerms[i] = forceField->angleBendTerms[i]; 

#ifdef ENABLE_OMP
#pragma omp parallel for private(i)
#endif
	for ( i = 0; i < numberOfBendTerms; i++ )
	{
		int ai, aj, ak;
		double term, term1, term2;
		double thetaRad;
		double delta = 1e-10;
		double cosine;
		double diff = 0;
		double diff2 = 0;
		double diff3 = 0;
		double diff4 = 0;

		double rijx, rijy, rijz;
		double rkjx, rkjy, rkjz;
		double rij2, rkj2;

		double gradix, gradiy, gradiz;
		double gradjx, gradjy, gradjz;
		double gradkx, gradky, gradkz;

		double rijDotrkj;
		double xp,yp,zp,rp;

		ai = (int)angleBendTerms[0][i];
		aj = (int)angleBendTerms[1][i];
		ak = (int)angleBendTerms[2][i];

		rijx = m->atoms[ai].coordinates[0] - m->atoms[aj].coordinates[0];
		rijy = m->atoms[ai].coordinates[1] - m->atoms[aj].coordinates[1];
		rijz = m->atoms[ai].coordinates[2] - m->atoms[aj].coordinates[2];
		rij2 = rijx * rijx + rijy * rijy + rijz * rijz;

		rkjx = m->atoms[ak].coordinates[0] - m->atoms[aj].coordinates[0];
		rkjy = m->atoms[ak].coordinates[1] - m->atoms[aj].coordinates[1];
		rkjz = m->atoms[ak].coordinates[2] - m->atoms[aj].coordinates[2];
		rkj2 = rkjx * rkjx + rkjy * rkjy + rkjz * rkjz;

		if (rij2==0 || rkj2==0) continue;

		rijDotrkj = rijx * rkjx + rijy * rkjy + rijz * rkjz;

		// corss rij^rkj
		xp = rkjy*rijz - rkjz*rijy;
		yp = rkjz*rijx - rkjx*rijz;
		zp = rkjx*rijy - rkjy*rijx;
		rp = sqrt(xp*xp + yp*yp + zp*zp);
		if ( rp < delta )
		{
			fprintf(forceField->logfile, "cut rp\n");
			fflush(forceField->logfile);
			rp = delta;
		printf("rp==0\n");
		}
	
	        cosine = rijDotrkj / sqrt(rij2*rkj2);
		if(cosine>1) cosine = 1;
		if(cosine<-1) cosine = -1;
                thetaRad = acos(cosine);
		diff =  thetaRad - D2R*angleBendTerms[4][i];
		diff2 = diff*diff;
		diff3 = diff2*diff;
		diff4 = diff3*diff;
		term =  angleBendTerms[3][i]*diff*
		(2.0+3.0*angleBendTerms[5][i]*diff+4.0*angleBendTerms[6][i]*diff2+5.0*angleBendTerms[7][i]*diff3+6.0*angleBendTerms[8][i]*diff4);

		term1 = -term/(rij2*rp);
		term2 = +term/(rkj2*rp);

		gradix = term1 * (rijy*zp-rijz*yp);
		gradiy = term1 * (rijz*xp-rijx*zp);
		gradiz = term1 * (rijx*yp-rijy*xp);

		gradkx = term2 * (rkjy*zp-rkjz*yp);
		gradky = term2 * (rkjz*xp-rkjx*zp);
		gradkz = term2 * (rkjx*yp-rkjy*xp);
			
		gradjx = -gradix -gradkx;
		gradjy = -gradiy -gradky;
		gradjz = -gradiz -gradkz;
			
#ifdef ENABLE_OMP
#pragma omp critical
#endif
		{
			m->atoms[ai].gradient[0] += gradix;
			m->atoms[ai].gradient[1] += gradiy;
			m->atoms[ai].gradient[2] += gradiz;
			
			m->atoms[aj].gradient[0] += gradjx;
			m->atoms[aj].gradient[1] += gradjy;
			m->atoms[aj].gradient[2] += gradjz;
			
			m->atoms[ak].gradient[0] += gradkx;
			m->atoms[ak].gradient[1] += gradky;
			m->atoms[ak].gradient[2] += gradkz;

/* energy */
			energy += angleBendTerms[3][i]*diff2*
			(1.0+angleBendTerms[5][i]*diff+angleBendTerms[6][i]*diff2+angleBendTerms[7][i]*diff3+angleBendTerms[8][i]*diff4);
		}
		//printf("Theta Equ diff = %f %f %f\n", thetaRad, D2R*angleBendTerms[4][i], diff);
		//printf("h = %f %f %f %f\n",angleBendTerms[5][i],angleBendTerms[6][i],angleBendTerms[7][i],angleBendTerms[8][i]);

	} 
	//printf("energy bend = %f\n",energy);
	m->potentialEnergy += energy;
}
/**********************************************************************/
static void calculateGradientStrBendAmber(ForceField* forceField)
{
	int i;

	Molecule* m = &forceField->molecule;
	double* strBendTerms[STRBENDDIM];
	static double D2R = 1/(RADTODEG);
	int numberOfStrBendTerms = forceField->numberOfStrBendTerms;
	double energy = 0;

	for( i=0; i<STRBENDDIM;i++) strBendTerms[i] = forceField->strBendTerms[i]; 

	//printf("numberOfStrBendTerms=%d\n",numberOfStrBendTerms);
#ifdef ENABLE_OMP
#pragma omp parallel for private(i)
#endif
	for ( i = 0; i < numberOfStrBendTerms; i++ )
	{
		int ai, aj, ak;
		//double thetaDeg, thetaRad, cosTheta;
		double thetaDeg;
		double absTheta;
		double delta = 1e-10;

		double rijx, rijy, rijz;
		double rkjx, rkjy, rkjz;
		double rij, rkj;
		double rij2, rkj2;
		double xp,yp,zp,rp;

		double gradix, gradiy, gradiz;
		double gradjx, gradjy, gradjz;
		double gradkx, gradky, gradkz;

		double diffAngle, diffR12, diffR23;
		double term1, term2;
		double termR, term1Angle, term2Angle;
		double termAngleix, termAngleiy, termAngleiz;
		double termAnglekx, termAngleky, termAnglekz;
		double termBondix, termBondiy,termBondiz;
		double termBondkx, termBondky,termBondkz;

		ai = (int)strBendTerms[2][i];
		aj = (int)strBendTerms[3][i];
		ak = (int)strBendTerms[4][i];

		thetaDeg = getAngle(&m->atoms[ai], &m->atoms[aj], &m->atoms[ak]);
		//thetaRad = thetaDeg * D2R;
		absTheta = fabs( thetaDeg );
		//cosTheta = cos( thetaRad );

		if ( ( absTheta > delta ) && ( absTheta < 180.0 - delta ) )
		{

			rijx = m->atoms[ai].coordinates[0] - m->atoms[aj].coordinates[0];
			rijy = m->atoms[ai].coordinates[1] - m->atoms[aj].coordinates[1];
			rijz = m->atoms[ai].coordinates[2] - m->atoms[aj].coordinates[2];

			rkjx = m->atoms[ak].coordinates[0] - m->atoms[aj].coordinates[0];
			rkjy = m->atoms[ak].coordinates[1] - m->atoms[aj].coordinates[1];
			rkjz = m->atoms[ak].coordinates[2] - m->atoms[aj].coordinates[2];

			rij2 = rijx * rijx + rijy * rijy + rijz * rijz;
			rij = sqrt(rij2);

			rkj2 = rkjx * rkjx + rkjy * rkjy + rkjz * rkjz;
			rkj = sqrt(rkj2);

			xp = rkjy*rijz - rkjz*rijy;
			yp = rkjz*rijx - rkjx*rijz;
			zp = rkjx*rijy - rkjy*rijx;
			rp = sqrt(xp*xp + yp*yp + zp*zp);
			if ( rp < delta )
			{
				fprintf(forceField->logfile, "cut rp\n");
				fflush(forceField->logfile);
				rp = delta;
			}
		//printf("Const = %f %f %f %f %f\n",strBendTerms[5][i],strBendTerms[6][i],strBendTerms[7][i],strBendTerms[8][i],strBendTerms[9][i]);
			// Terms for the bond angle derivative
          		diffAngle = D2R*(thetaDeg - strBendTerms[9][i]);
               		term1 = -1.0 / (rij2*rp);
               		term2 =  1.0 / (rkj2*rp);
               		termAngleix = term1 * (rijy*zp-rijz*yp);
               		termAngleiy = term1 * (rijz*xp-rijx*zp);
               		termAngleiz = term1 * (rijx*yp-rijy*xp);

               		termAnglekx = term2 * (rkjy*zp-rkjz*yp);
               		termAngleky = term2 * (rkjz*xp-rkjx*zp);
               		termAnglekz = term2 * (rkjx*yp-rkjy*xp);
			//  Terms for the bond length derivatives
          		diffR12 = rij - strBendTerms[7][i];
          		diffR23 = rkj - strBendTerms[8][i];
               		term1 = 1.0 / rij;
               		term2 = 1.0 / rkj;
               		termBondix = term1 * rijx;
               		termBondiy = term1 * rijy;
               		termBondiz = term1 * rijz;

               		termBondkx = term2 * rkjx;
               		termBondky = term2 * rkjy;
               		termBondkz = term2 * rkjz;

               		termR = strBendTerms[5][i]*diffR12 + strBendTerms[6][i]*diffR23;
               		term1Angle = strBendTerms[5][i] * diffAngle;
               		term2Angle = strBendTerms[6][i] * diffAngle;

			gradix = term1Angle * termBondix + termR*termAngleix;
			gradiy = term1Angle * termBondiy + termR*termAngleiy;
			gradiz = term1Angle * termBondiz + termR*termAngleiz;
			
			gradkx = term2Angle * termBondkx + termR*termAnglekx;
			gradky = term2Angle * termBondky + termR*termAngleky;
			gradkz = term2Angle * termBondkz + termR*termAnglekz;

			gradjx = -gradix-gradkx;
			gradjy = -gradiy-gradky;
			gradjz = -gradiz-gradkz;
			
#ifdef ENABLE_OMP
#pragma omp critical
#endif
			{
			m->atoms[ai].gradient[0] += gradix;
			m->atoms[ai].gradient[1] += gradiy;
			m->atoms[ai].gradient[2] += gradiz;
			
			m->atoms[aj].gradient[0] += gradjx;
			m->atoms[aj].gradient[1] += gradjy;
			m->atoms[aj].gradient[2] += gradjz;
			
			m->atoms[ak].gradient[0] += gradkx;
			m->atoms[ak].gradient[1] += gradky;
			m->atoms[ak].gradient[2] += gradkz;
			}
/* energy */
			//printf("ener str-bend  in grad= %f\n",termR*diffAngle);
			energy += termR*diffAngle;
		}
	} 
	m->potentialEnergy += energy;
}
/**********************************************************************/
void calculateGradientDihedralAmber(ForceField* forceField)
{

	int i;

	Molecule* m = &forceField->molecule;
	double* dihedralAngleTerms[DIHEDRALDIM];
	static double D2R = 1/(RADTODEG);
	int numberOfDihedralTerms = forceField->numberOfDihedralTerms;
	double phiDeg;
	double energy = 0;

	for(i=0;i<DIHEDRALDIM;i++)
		dihedralAngleTerms[i] = forceField->dihedralAngleTerms[i];

#ifdef ENABLE_OMP
#pragma omp parallel for private(i) 
#endif
	for (  i = 0; i < numberOfDihedralTerms; i++ )
	{
		int ai, aj, ak, al;
		int j;

		double rjix, rjiy, rjiz;
		double rkjx, rkjy, rkjz;
		double rkix, rkiy, rkiz;
		double rljx, rljy, rljz;
		double rlkx, rlky, rlkz;

		double gradix, gradiy, gradiz;
		double gradjx, gradjy, gradjz;
		double gradkx, gradky, gradkz;
		double gradlx, gradly, gradlz;

		double rkj;
		double xt, yt, zt;
		double xu, yu, zu;
		double xtu, ytu, ztu;
		double rt2, ru2, rtru;
		double cosine1, sine1, cosineN, sineN, cosold, sinold;
		double cosPhase, sinPhase;
		double dedxt, dedyt, dedzt;
		double dedxu, dedyu, dedzu;
		double dedphi;
		int n;
		double vn;

		ai = (int)dihedralAngleTerms[0][i];
		aj = (int)dihedralAngleTerms[1][i];
		ak = (int)dihedralAngleTerms[2][i];
		al = (int)dihedralAngleTerms[3][i];

		rjix = m->atoms[aj].coordinates[0] - m->atoms[ai].coordinates[0];
		rjiy = m->atoms[aj].coordinates[1] - m->atoms[ai].coordinates[1];
		rjiz = m->atoms[aj].coordinates[2] - m->atoms[ai].coordinates[2];

		rkjx = m->atoms[ak].coordinates[0] - m->atoms[aj].coordinates[0];
		rkjy = m->atoms[ak].coordinates[1] - m->atoms[aj].coordinates[1];
		rkjz = m->atoms[ak].coordinates[2] - m->atoms[aj].coordinates[2];

		rlkx = m->atoms[al].coordinates[0] - m->atoms[ak].coordinates[0];
		rlky = m->atoms[al].coordinates[1] - m->atoms[ak].coordinates[1];
		rlkz = m->atoms[al].coordinates[2] - m->atoms[ak].coordinates[2];


		xt = rjiy*rkjz - rkjy*rjiz;
		yt = rjiz*rkjx - rkjz*rjix;
		zt = rjix*rkjy - rkjx*rjiy;

		xu = rkjy*rlkz - rlky*rkjz;
		yu = rkjz*rlkx - rlkz*rkjx;
		zu = rkjx*rlky - rlkx*rkjy;

		xtu = yt*zu - yu*zt;
		ytu = zt*xu - zu*xt;
		ztu = xt*yu - xu*yt;

		rt2 = xt*xt + yt*yt + zt*zt;
		ru2 = xu*xu + yu*yu + zu*zu;

		rtru = sqrt(rt2 * ru2);

		rkj = sqrt(rkjx*rkjx + rkjy*rkjy + rkjz*rkjz);
		cosine1 = 1.0;
		sine1   = 0.0;

		if (rtru <1e-10) rtru = 1e-10;
		if (rt2 <1e-10) rt2 = 1e-10;
		if (ru2 <1e-10) ru2 = 1e-10;

		cosine1 = (xt*xu + yt*yu + zt*zu) / rtru;
		sine1 = (rkjx*xtu + rkjy*ytu + rkjz*ztu) / (rkj*rtru);

		n = (int)dihedralAngleTerms[7][i];
		cosPhase = cos(D2R*dihedralAngleTerms[6][i]);
		sinPhase = sin(D2R*dihedralAngleTerms[6][i]);
		vn = dihedralAngleTerms[5][i]/dihedralAngleTerms[4][i];

/*
     compute the multiple angle trigonometry and the phase terms
*/
		
		cosineN = cosine1;
		sineN   = sine1;

		for(j=2;j<=n;j++)
		{
		   cosold = cosineN;
		   sinold = sineN;
		   cosineN = cosine1*cosold - sine1*sinold;
		   sineN   = cosine1*sinold + sine1*cosold;
		}

		dedphi = vn*n*(cosineN*sinPhase-sineN*cosPhase);

/*
     chain rule terms for first derivative components
*/

		rkix = m->atoms[ak].coordinates[0] - m->atoms[ai].coordinates[0];
		rkiy = m->atoms[ak].coordinates[1] - m->atoms[ai].coordinates[1];
		rkiz = m->atoms[ak].coordinates[2] - m->atoms[ai].coordinates[2];

		rljx = m->atoms[al].coordinates[0] - m->atoms[aj].coordinates[0];
		rljy = m->atoms[al].coordinates[1] - m->atoms[aj].coordinates[1];
		rljz = m->atoms[al].coordinates[2] - m->atoms[aj].coordinates[2];

		dedxt = dedphi * (yt*rkjz - rkjy*zt) / (rt2*rkj);
		dedyt = dedphi * (zt*rkjx - rkjz*xt) / (rt2*rkj);
		dedzt = dedphi * (xt*rkjy - rkjx*yt) / (rt2*rkj);

		dedxu = -dedphi * (yu*rkjz - rkjy*zu) / (ru2*rkj);
		dedyu = -dedphi * (zu*rkjx - rkjz*xu) / (ru2*rkj);
		dedzu = -dedphi * (xu*rkjy - rkjx*yu) / (ru2*rkj);
/*

     compute first derivative components for this angle
*/

		gradix = rkjz*dedyt - rkjy*dedzt;
		gradiy = rkjx*dedzt - rkjz*dedxt;
		gradiz = rkjy*dedxt - rkjx*dedyt;

		gradjx = rkiy*dedzt - rkiz*dedyt + rlkz*dedyu - rlky*dedzu;
		gradjy = rkiz*dedxt - rkix*dedzt + rlkx*dedzu - rlkz*dedxu;
		gradjz = rkix*dedyt - rkiy*dedxt + rlky*dedxu - rlkx*dedyu;

		gradkx = rjiz*dedyt - rjiy*dedzt + rljy*dedzu - rljz*dedyu;
		gradky = rjix*dedzt - rjiz*dedxt + rljz*dedxu - rljx*dedzu;
		gradkz = rjiy*dedxt - rjix*dedyt + rljx*dedyu - rljy*dedxu;

		gradlx = rkjz*dedyu - rkjy*dedzu;
		gradly = rkjx*dedzu - rkjz*dedxu;
		gradlz = rkjy*dedxu - rkjx*dedyu;

#ifdef ENABLE_OMP
#pragma omp critical
#endif
		{
		m->atoms[ai].gradient[0] += gradix;
		m->atoms[ai].gradient[1] += gradiy;
		m->atoms[ai].gradient[2] += gradiz;

		m->atoms[aj].gradient[0] += gradjx;
		m->atoms[aj].gradient[1] += gradjy;
		m->atoms[aj].gradient[2] += gradjz;

		m->atoms[ak].gradient[0] += gradkx;
		m->atoms[ak].gradient[1] += gradky;
		m->atoms[ak].gradient[2] += gradkz;

		m->atoms[al].gradient[0] += gradlx;
		m->atoms[al].gradient[1] += gradly;
		m->atoms[al].gradient[2] += gradlz;
		}
		phiDeg = getTorsion(  &m->atoms[ai] ,&m->atoms[aj], &m->atoms[ak], &m->atoms[al]);
		//phiDeg = acos(cosine1)*RADTODEG;
		energy += dihedralAngleTerms[5][i]/dihedralAngleTerms[4][i] * 
		( 1.0 + cos( D2R*(dihedralAngleTerms[7][i] * phiDeg - dihedralAngleTerms[6][i] )) );
	}
	m->potentialEnergy += energy;
}
/**********************************************************************/
static void calculateGradientImproperTorsion(ForceField* forceField)
{
	int i;
	Molecule* m = &forceField->molecule;
	double forceConstant, equilibriumAngle;
	double* improperTorsionTerms[IMPROPERDIHEDRALDIM];
	static double D2R = 1/(RADTODEG);
	double energy = 0;
	double phiDeg;
	int numberOfImproperTorsionTerms = forceField->numberOfImproperTorsionTerms;

	for(i=0;i<IMPROPERDIHEDRALDIM;i++) improperTorsionTerms[i] = forceField->improperTorsionTerms[i];

#ifdef ENABLE_OMP
#pragma omp parallel for private(i) 
#endif
	for (  i = 0; i < numberOfImproperTorsionTerms; i++ )
	{
		int ai, aj, ak, al;
		int j;

		double rjix, rjiy, rjiz;
		double rkjx, rkjy, rkjz;
		double rkix, rkiy, rkiz;
		double rljx, rljy, rljz;
		double rlkx, rlky, rlkz;

		double gradix, gradiy, gradiz;
		double gradjx, gradjy, gradjz;
		double gradkx, gradky, gradkz;
		double gradlx, gradly, gradlz;

		double rkj;
		double xt, yt, zt;
		double xu, yu, zu;
		double xtu, ytu, ztu;
		double rt2, ru2, rtru;
		double cosine1, sine1, cosineN, sineN, cosold, sinold;
		double cosPhase, sinPhase;
		double dedxt, dedyt, dedzt;
		double dedxu, dedyu, dedzu;
		double dedphi;
		int n;

		ai = (int)improperTorsionTerms[0][i];
		aj = (int)improperTorsionTerms[1][i];
		ak = (int)improperTorsionTerms[2][i];
		al = (int)improperTorsionTerms[3][i];

		forceConstant = improperTorsionTerms[4][i];
		n = (int)improperTorsionTerms[6][i];
                equilibriumAngle = improperTorsionTerms[5][i];

		rjix = m->atoms[aj].coordinates[0] - m->atoms[ai].coordinates[0];
		rjiy = m->atoms[aj].coordinates[1] - m->atoms[ai].coordinates[1];
		rjiz = m->atoms[aj].coordinates[2] - m->atoms[ai].coordinates[2];

		rkjx = m->atoms[ak].coordinates[0] - m->atoms[aj].coordinates[0];
		rkjy = m->atoms[ak].coordinates[1] - m->atoms[aj].coordinates[1];
		rkjz = m->atoms[ak].coordinates[2] - m->atoms[aj].coordinates[2];

		rlkx = m->atoms[al].coordinates[0] - m->atoms[ak].coordinates[0];
		rlky = m->atoms[al].coordinates[1] - m->atoms[ak].coordinates[1];
		rlkz = m->atoms[al].coordinates[2] - m->atoms[ak].coordinates[2];


		xt = rjiy*rkjz - rkjy*rjiz;
		yt = rjiz*rkjx - rkjz*rjix;
		zt = rjix*rkjy - rkjx*rjiy;

		xu = rkjy*rlkz - rlky*rkjz;
		yu = rkjz*rlkx - rlkz*rkjx;
		zu = rkjx*rlky - rlkx*rkjy;

		xtu = yt*zu - yu*zt;
		ytu = zt*xu - zu*xt;
		ztu = xt*yu - xu*yt;

		rt2 = xt*xt + yt*yt + zt*zt;
		ru2 = xu*xu + yu*yu + zu*zu;

		rtru = sqrt(rt2 * ru2);

		rkj = sqrt(rkjx*rkjx + rkjy*rkjy + rkjz*rkjz);
		cosine1 = 1.0;
		sine1   = 0.0;

		if (rtru <1e-10) rtru = 1e-10;
		if (rt2 <1e-10) rt2 = 1e-10;
		if (ru2 <1e-10) ru2 = 1e-10;

		cosine1 = (xt*xu + yt*yu + zt*zu) / rtru;
		sine1 = (rkjx*xtu + rkjy*ytu + rkjz*ztu) / (rkj*rtru);

		cosPhase = cos(D2R*equilibriumAngle);
		sinPhase = sin(D2R*equilibriumAngle);

/*
     compute the multiple angle trigonometry and the phase terms
*/
		
		cosineN = cosine1;
		sineN   = sine1;

		for(j=2;j<=n;j++)
		{
		   cosold = cosineN;
		   sinold = sineN;
		   cosineN = cosine1*cosold - sine1*sinold;
		   sineN   = cosine1*sinold + sine1*cosold;
		}

		dedphi = forceConstant*n*(cosineN*sinPhase-sineN*cosPhase);

/*
     chain rule terms for first derivative components
*/

		rkix = m->atoms[ak].coordinates[0] - m->atoms[ai].coordinates[0];
		rkiy = m->atoms[ak].coordinates[1] - m->atoms[ai].coordinates[1];
		rkiz = m->atoms[ak].coordinates[2] - m->atoms[ai].coordinates[2];

		rljx = m->atoms[al].coordinates[0] - m->atoms[aj].coordinates[0];
		rljy = m->atoms[al].coordinates[1] - m->atoms[aj].coordinates[1];
		rljz = m->atoms[al].coordinates[2] - m->atoms[aj].coordinates[2];

		dedxt = dedphi * (yt*rkjz - rkjy*zt) / (rt2*rkj);
		dedyt = dedphi * (zt*rkjx - rkjz*xt) / (rt2*rkj);
		dedzt = dedphi * (xt*rkjy - rkjx*yt) / (rt2*rkj);

		dedxu = -dedphi * (yu*rkjz - rkjy*zu) / (ru2*rkj);
		dedyu = -dedphi * (zu*rkjx - rkjz*xu) / (ru2*rkj);
		dedzu = -dedphi * (xu*rkjy - rkjx*yu) / (ru2*rkj);
/*

     compute first derivative components for this angle
*/

		gradix = rkjz*dedyt - rkjy*dedzt;
		gradiy = rkjx*dedzt - rkjz*dedxt;
		gradiz = rkjy*dedxt - rkjx*dedyt;

		gradjx = rkiy*dedzt - rkiz*dedyt + rlkz*dedyu - rlky*dedzu;
		gradjy = rkiz*dedxt - rkix*dedzt + rlkx*dedzu - rlkz*dedxu;
		gradjz = rkix*dedyt - rkiy*dedxt + rlky*dedxu - rlkx*dedyu;

		gradkx = rjiz*dedyt - rjiy*dedzt + rljy*dedzu - rljz*dedyu;
		gradky = rjix*dedzt - rjiz*dedxt + rljz*dedxu - rljx*dedzu;
		gradkz = rjiy*dedxt - rjix*dedyt + rljx*dedyu - rljy*dedxu;

		gradlx = rkjz*dedyu - rkjy*dedzu;
		gradly = rkjx*dedzu - rkjz*dedxu;
		gradlz = rkjy*dedxu - rkjx*dedyu;

#ifdef ENABLE_OMP
#pragma omp critical
#endif
		{
		m->atoms[ai].gradient[0] += gradix;
		m->atoms[ai].gradient[1] += gradiy;
		m->atoms[ai].gradient[2] += gradiz;

		m->atoms[aj].gradient[0] += gradjx;
		m->atoms[aj].gradient[1] += gradjy;
		m->atoms[aj].gradient[2] += gradjz;

		m->atoms[ak].gradient[0] += gradkx;
		m->atoms[ak].gradient[1] += gradky;
		m->atoms[ak].gradient[2] += gradkz;

		m->atoms[al].gradient[0] += gradlx;
		m->atoms[al].gradient[1] += gradly;
		m->atoms[al].gradient[2] += gradlz;
		}
		phiDeg = getTorsion(  &m->atoms[ai] ,&m->atoms[aj], &m->atoms[ak], &m->atoms[al]);
		//phiDeg = acos(cosine1)*RADTODEG;
		energy += forceConstant * ( 1.0 + cos( D2R*(n * phiDeg - equilibriumAngle )) );
	}
	//printf("ImproperTorsion energy =%f\n",energy);
	m->potentialEnergy += energy;
}
/*****************************************************************************************/
void calculateGradientOutOfPlaneAmber(ForceField* forceField)
{

	int i;

	Molecule* m = &forceField->molecule;
	double* outOfPlaneTerms[OUTOFPLANEDIM];
	int numberOfOutOfPlaneTerms = forceField->numberOfOutOfPlaneTerms;
	double energy = 0;

	for(i=0;i<OUTOFPLANEDIM;i++) outOfPlaneTerms[i] = forceField->outOfPlaneTerms[i];

#ifdef ENABLE_OMP
#pragma omp parallel for private(i) 
#endif
	for (  i = 0; i < numberOfOutOfPlaneTerms; i++ )
	{
		int ai, aj, ak, al;

		double rijx, rijy, rijz;
		double rkjx, rkjy, rkjz;
		double rljx, rljy, rljz;

		double gradix, gradiy, gradiz;
		double gradjx, gradjy, gradjz;
		double gradkx, gradky, gradkz;
		double gradlx, gradly, gradlz;

		//double angle, e, e2, c2;
		double angle,  e2, c2;
		double dot;
		double rilx, rily, rilz;
		double rklx, rkly, rklz;
		double rij2 =0.0, rkj2 = 0.0, ril2 = 0.0, rkl2 = 0.0, rlj2 = 0.0;
		double bkk2, cosine;
		double dt, dt2, dt3, dt4;
		double deddt, dEdCos;
		double term;
		double dc2dix, dc2diy, dc2diz;
		double dc2dkx, dc2dky, dc2dkz;
		double dc2dlx, dc2dly, dc2dlz;
		double de2dix, de2diy, de2diz;
		double de2dkx, de2dky, de2dkz;
		double de2dlx, de2dly, de2dlz;

		boolean wdc = (outOfPlaneTerms[4][i]==1);
		double force = outOfPlaneTerms[5][i];
		double h3 = outOfPlaneTerms[6][i];
		double h4 = outOfPlaneTerms[7][i];
		double h5 = outOfPlaneTerms[8][i];
		double h6 = outOfPlaneTerms[9][i];

		ai = (int)outOfPlaneTerms[0][i];
		aj = (int)outOfPlaneTerms[1][i];
		ak = (int)outOfPlaneTerms[2][i];
		al = (int)outOfPlaneTerms[3][i];

		rijx = m->atoms[ai].coordinates[0] - m->atoms[aj].coordinates[0];
		rijy = m->atoms[ai].coordinates[1] - m->atoms[aj].coordinates[1];
		rijz = m->atoms[ai].coordinates[2] - m->atoms[aj].coordinates[2];

		rkjx = m->atoms[ak].coordinates[0] - m->atoms[aj].coordinates[0];
		rkjy = m->atoms[ak].coordinates[1] - m->atoms[aj].coordinates[1];
		rkjz = m->atoms[ak].coordinates[2] - m->atoms[aj].coordinates[2];

		rljx = m->atoms[al].coordinates[0] - m->atoms[aj].coordinates[0];
		rljy = m->atoms[al].coordinates[1] - m->atoms[aj].coordinates[1];
		rljz = m->atoms[al].coordinates[2] - m->atoms[aj].coordinates[2];

		rilx = m->atoms[ai].coordinates[0] - m->atoms[al].coordinates[0];
		rily = m->atoms[ai].coordinates[1] - m->atoms[al].coordinates[1];
		rilz = m->atoms[ai].coordinates[2] - m->atoms[al].coordinates[2];

		rklx = m->atoms[ak].coordinates[0] - m->atoms[al].coordinates[0];
		rkly = m->atoms[ak].coordinates[1] - m->atoms[al].coordinates[1];
		rklz = m->atoms[ak].coordinates[2] - m->atoms[al].coordinates[2];

		if(wdc) //  'W-D-C'
		{
		     	rij2 = rijx*rijx + rijy*rijy + rijz*rijz;
               		rkj2 = rkjx*rkjx + rkjy*rkjy + rkjz*rkjz;
               		dot = rijx*rkjx+rijy*rkjy+rijz*rkjz;
               		c2 = rij2*rkj2 - dot*dot;
		}
		else // ALLINGER
		{
			ril2 = rilx*rilx + rily*rily + rilz*rilz;
               		rkl2 = rklx*rklx + rkly*rkly + rklz*rklz;
               		dot = rilx*rklx + rily*rkly + rilz*rklz;
               		c2 = ril2*rkl2 - dot*dot;
		}
		// energy
		e2 = rljx*(rijy*rkjz-rijz*rkjy) + rljy*(rijz*rkjx-rijx*rkjz)+ rljz*(rijx*rkjy-rijy*rkjx);
            	rlj2 = rljx*rljx + rljy*rljy + rljz*rljz;
            	if (rlj2==0.0 || c2 == 0.0) continue;
               	bkk2 = rlj2 - e2*e2/c2;
               	cosine = sqrt(bkk2/rlj2);
		if(cosine>1) cosine = 1;
		if(cosine<-1) cosine = -1;
               	angle = acos(cosine);
               	dt = angle;
               	dt2 = dt * dt;
               	dt3 = dt2 * dt;
               	dt4 = dt2 * dt2;
               	//e = force * dt2*(1.0+h3*dt+h4*dt2+h5*dt3+h6*dt4);
               	deddt = force * dt * (2.0 + 3.0*h3*dt + 4.0*h4*dt2 +5.0*h5*dt3 + 6.0*h6*dt4);
               	dEdCos = -deddt / sqrt(c2*bkk2);
		if(e2<0) dEdCos = deddt / sqrt(c2*bkk2);
		//  first derivative components
		if(wdc)
		{
			term = e2 / c2;
                  	dc2dix = (rijx*rkj2-rkjx*dot) * term;
                  	dc2diy = (rijy*rkj2-rkjy*dot) * term;
                  	dc2diz = (rijz*rkj2-rkjz*dot) * term;
                  	dc2dkx = (rkjx*rij2-rijx*dot) * term;
                  	dc2dky = (rkjy*rij2-rijy*dot) * term;
                  	dc2dkz = (rkjz*rij2-rijz*dot) * term;
                  	dc2dlx = 0.0;
                  	dc2dly = 0.0;
                  	dc2dlz = 0.0;

		}
		else
		{
			term = e2 / c2;
                  	dc2dix = (rilx*rkl2-rklx*dot) * term;
                  	dc2diy = (rily*rkl2-rkly*dot) * term;
                  	dc2diz = (rilz*rkl2-rklz*dot) * term;
                  	dc2dkx = (rklx*ril2-rilx*dot) * term;
                  	dc2dky = (rkly*ril2-rily*dot) * term;
                  	dc2dkz = (rklz*ril2-rilz*dot) * term;
                  	dc2dlx = -dc2dix - dc2dkx;
                  	dc2dly = -dc2diy - dc2dky;
                  	dc2dlz = -dc2diz - dc2dkz;
		}
		term = e2 / rlj2;
               	de2dix = rljy*rkjz - rljz*rkjy;
               	de2diy = rljz*rkjx - rljx*rkjz;
               	de2diz = rljx*rkjy - rljy*rkjx;
               	de2dkx = rijy*rljz - rijz*rljy;
               	de2dky = rijz*rljx - rijx*rljz;
               	de2dkz = rijx*rljy - rijy*rljx;
               	de2dlx = rkjy*rijz - rkjz*rijy + rljx*term;
               	de2dly = rkjz*rijx - rkjx*rijz + rljy*term;
               	de2dlz = rkjx*rijy - rkjy*rijx + rljz*term;


	       	gradix = dEdCos * (dc2dix+de2dix);
               	gradiy = dEdCos * (dc2diy+de2diy);
               	gradiz = dEdCos * (dc2diz+de2diz);
               	gradkx = dEdCos * (dc2dkx+de2dkx);
               	gradky = dEdCos * (dc2dky+de2dky);
               	gradkz = dEdCos * (dc2dkz+de2dkz);
               	gradlx = dEdCos * (dc2dlx+de2dlx);
               	gradly = dEdCos * (dc2dly+de2dly);
               	gradlz = dEdCos * (dc2dlz+de2dlz);
               	gradjx = -gradix - gradkx - gradlx;
               	gradjy = -gradiy - gradky - gradly;
               	gradjz = -gradiz - gradkz - gradlz;
#ifdef ENABLE_OMP
#pragma omp critical
#endif
		{
			m->atoms[ai].gradient[0] += gradix;
			m->atoms[ai].gradient[1] += gradiy;
			m->atoms[ai].gradient[2] += gradiz;

			m->atoms[aj].gradient[0] += gradjx;
			m->atoms[aj].gradient[1] += gradjy;
			m->atoms[aj].gradient[2] += gradjz;

			m->atoms[ak].gradient[0] += gradkx;
			m->atoms[ak].gradient[1] += gradky;
			m->atoms[ak].gradient[2] += gradkz;

			m->atoms[al].gradient[0] += gradlx;
			m->atoms[al].gradient[1] += gradly;
			m->atoms[al].gradient[2] += gradlz;
		}
		energy += e2;
	}
	m->potentialEnergy += energy;
}
/**********************************************************************/
void calculateGradientVdw612Amber(ForceField* forceField)
{
	int i;
	double energy = 0;

	Molecule* m = &forceField->molecule;
	double* vdw612Terms[VDW612DIM];
	int numberOfVdw612 = forceField->numberOfVdw612;


	//printf("calculateGradientVdw612Amber : useCoulomb=%d\n",useCoulomb);
	for(i=0;i<VDW612DIM;i++)
		vdw612Terms[i] = forceField->vdw612Terms[i];

	/* non-bonded part */
#ifdef ENABLE_OMP
#pragma omp parallel for private(i) 
#endif
	for (  i = 0; i < numberOfVdw612; i++ )
	{
		int ai, aj;
		double rijx, rijy, rijz;
		double gradix, gradiy, gradiz;
		double rij2, rij;
		double rijm3;
		double rijm6, rijm8;
		double  term3;
		double A,B;

		ai       = (int)vdw612Terms[0][i];
		aj       = (int)vdw612Terms[1][i];
		A = vdw612Terms[2][i];
		B = vdw612Terms[3][i];

		rijx =  m->atoms[ai].coordinates[0] - m->atoms[aj].coordinates[0];
		rijy =  m->atoms[ai].coordinates[1] - m->atoms[aj].coordinates[1];
		rijz =  m->atoms[ai].coordinates[2] - m->atoms[aj].coordinates[2];

		rij2 = rijx * rijx + rijy * rijy + rijz * rijz;
		if ( rij2 < POSEPS*POSEPS ) rij2 = POSEPS*POSEPS;	

		rij = sqrt( rij2 );
		rijm3 = 1.0/(rij2 * rij);
		rijm6 = rijm3 * rijm3;
		rijm8 = rijm6 / rij2;
		
		term3 = -6*(2 * A * rijm6 - B)*rijm8;
		gradix = term3 * rijx;
		gradiy = term3 * rijy;
		gradiz = term3 * rijz;
#ifdef ENABLE_OMP
#pragma omp critical
#endif
		{
			m->atoms[ai].gradient[0] += gradix;
			m->atoms[ai].gradient[1] += gradiy;
			m->atoms[ai].gradient[2] += gradiz;

			m->atoms[aj].gradient[0] -= gradix;
			m->atoms[aj].gradient[1] -= gradiy;
			m->atoms[aj].gradient[2] -= gradiz;
		}
		energy += (A*rijm6 - B) * rijm6;
		if(rij<1) printf("Warning : rij2 = %f atom# = %d atom# = %d\n",rij,ai+1,aj+1);
	//	printf("rij2 = %f ai = %d aj = %d\n",rij2,ai+1,aj+1);
	}  
	m->potentialEnergy += energy;
}
/**********************************************************************/
void calculateGradientVdw714Amber(ForceField* forceField)
{
	int i;
	double energy = 0;

	Molecule* m = &forceField->molecule;
	double* vdw714Terms[VDW714DIM];
	int numberOfVdw714 = forceField->numberOfVdw714;


	//printf("calculateGradientVdw714Amber : useCoulomb=%d\n",useCoulomb);
	for(i=0;i<VDW714DIM;i++) vdw714Terms[i] = forceField->vdw714Terms[i];

	/* non-bonded part */
#ifdef ENABLE_OMP
#pragma omp parallel for private(i) 
#endif
	for (  i = 0; i < numberOfVdw714; i++ )
	{
		double rijx, rijy, rijz;
		double gradix, gradiy, gradiz;
		double gradjx, gradjy, gradjz;
		double rij2, rij;
		int ai, aj;
		double epsilon, gamma, delta, R0;
		double rho;
		double e,de;
		//double R02,R04,R07;
		double R02,R07;
		double rij6,rij7;
		double tau, dtau, tau7, gtau;

		ai       = (int)vdw714Terms[0][i];
		aj       = (int)vdw714Terms[1][i];

		epsilon  = vdw714Terms[2][i];
		R0       = vdw714Terms[3][i];
		gamma    = vdw714Terms[4][i];
		delta    = vdw714Terms[5][i];

		rijx =  m->atoms[ai].coordinates[0] - m->atoms[aj].coordinates[0];
		rijy =  m->atoms[ai].coordinates[1] - m->atoms[aj].coordinates[1];
		rijz =  m->atoms[ai].coordinates[2] - m->atoms[aj].coordinates[2];

		rij2 = rijx * rijx + rijy * rijy + rijz * rijz;
		if ( rij2 < POSEPS*POSEPS ) rij2 = POSEPS*POSEPS;	

		rij = sqrt( rij2 );
	printf("rij = %f R0 = %f epsilon = %f gamma = %f delta = %f\n",rij,R0,epsilon,gamma,delta);
		R02 = R0*R0;
		//R04 = R02*R02;
		R07 = R02*R02*R02*R0;
		rij6 = rij2*rij2*rij2;
		rij7 = rij6*rij;
		rho = rij7+gamma*R07;

                tau = (delta+1.0) / (rij + delta*R0);
                tau7 = pow(tau,7);
                dtau = tau / (delta+1.0);
                gtau = epsilon*tau7*rij6*(gamma+1.0)*(R07/rho)*(R07/rho);
                e = epsilon*tau7*R07*((gamma+1.0)*R07/rho-2.0);
                de = 7.0 * (dtau*e+gtau);
		de = de/rij;

		gradix = de * rijx;
		gradiy = de * rijy;
		gradiz = de * rijz;

		gradjx = -gradix;
		gradjy = -gradiy;
		gradjz = -gradiz;
#ifdef ENABLE_OMP
#pragma omp critical
#endif
		{
		m->atoms[ai].gradient[0] += gradix;
		m->atoms[ai].gradient[1] += gradiy;
		m->atoms[ai].gradient[2] += gradiz;

		m->atoms[aj].gradient[0] += gradjx;
		m->atoms[aj].gradient[1] += gradjy;
		m->atoms[aj].gradient[2] += gradjz;

		energy += e;
		}
		if(rij<1) printf("Warning : rij2 = %f atom# = %d atom# = %d\n",rij,ai+1,aj+1);
	//	printf("rij2 = %f ai = %d aj = %d\n",rij2,ai+1,aj+1);
	}  
	//printf("E in grad 7-14 %f\n",energy);
	m->potentialEnergy += energy;
}
/**********************************************************************/
void calculateGradientVdw714AmberOld(ForceField* forceField)
{
	int i;
	double energy = 0;

	Molecule* m = &forceField->molecule;
	double* vdw714Terms[VDW714DIM];
	int numberOfVdw714 = forceField->numberOfVdw714;


	//printf("calculateGradientVdw714Amber : useCoulomb=%d\n",useCoulomb);
	for(i=0;i<VDW714DIM;i++) vdw714Terms[i] = forceField->vdw714Terms[i];

	/* non-bonded part */
#ifdef ENABLE_OMP
#pragma omp parallel for private(i) 
#endif
	for (  i = 0; i < numberOfVdw714; i++ )
	{
		double rijx, rijy, rijz;
		double gradix, gradiy, gradiz;
		double gradjx, gradjy, gradjz;
		double rij2, rij;
		int ai, aj;
		double epsilon, gamma, delta, R0;
		double R02;
		double d,d2,d4,d7;
		double sr,sr2,sr4,sr7,sr8;
		double rho,rho2,rho4,rho5, rho7;
		double srijg7;
		double e7,e14;
		double term, term7, term14, term14_1, term14_2;

		ai       = (int)vdw714Terms[0][i];
		aj       = (int)vdw714Terms[1][i];

		epsilon  = vdw714Terms[2][i];
		R0       = vdw714Terms[3][i];
		gamma    = vdw714Terms[4][i];
		delta    = vdw714Terms[5][i];

		rijx =  m->atoms[ai].coordinates[0] - m->atoms[aj].coordinates[0];
		rijy =  m->atoms[ai].coordinates[1] - m->atoms[aj].coordinates[1];
		rijz =  m->atoms[ai].coordinates[2] - m->atoms[aj].coordinates[2];

		rij2 = rijx * rijx + rijy * rijy + rijz * rijz;
		if ( rij2 < POSEPS*POSEPS ) rij2 = POSEPS*POSEPS;	

		rij = sqrt( rij2 );
		R02 = R0*R0;

		rho = rij/R0;
		rho2 = rho*rho;
		rho4 = rho2*rho2;
		rho5 = rho4*rho;
		rho7 = rho5*rho2;

		d = (1+delta);
		d2 = d*d;
		d4 = d2*d2;
		d7 = d4*d2*d;

		sr = 1.0/(rho+delta);
		sr2 = sr*sr;
		sr4 = sr2*sr2;
		sr7 = sr4*sr2*sr;
		sr8 = sr4*sr4;
		srijg7=1.0/(rho7+gamma);
		
		e7  = epsilon*d7*sr7;
		e14 = e7*(1+gamma)*srijg7;
		e7  = 2*e7;

		term7 = -14*epsilon*d7*sr8/(rij*R0);
		term14  = 7*epsilon*d7*(1+gamma)*sr7*srijg7;
		term14_1 = term14*sr/(rij*R0);
		term14_2 = term14*srijg7*rho5/R02;
		term = term7 + term14_1 + term14_2;
		term = - term;

		gradix = term * rijx;
		gradiy = term * rijy;
		gradiz = term * rijz;

		gradjx = -gradix;
		gradjy = -gradiy;
		gradjz = -gradiz;
#ifdef ENABLE_OMP
#pragma omp critical
#endif
		{
		m->atoms[ai].gradient[0] += gradix;
		m->atoms[ai].gradient[1] += gradiy;
		m->atoms[ai].gradient[2] += gradiz;

		m->atoms[aj].gradient[0] += gradjx;
		m->atoms[aj].gradient[1] += gradjy;
		m->atoms[aj].gradient[2] += gradjz;

		energy += e14-e7;
		}
		if(rij<1) printf("Warning : rij2 = %f atom# = %d atom# = %d\n",rij,ai+1,aj+1);
	//	printf("rij2 = %f ai = %d aj = %d\n",rij2,ai+1,aj+1);
	}  
	//printf("E in grad 7-14 %f\n",energy);
	m->potentialEnergy += energy;
}
/**********************************************************************/
void calculateGradientCoulombAmberEEMACKS2(ForceField* forceField)
{
	double dx = forceField->options.dx;
	Molecule* m = &forceField->molecule;
	int nAtoms = m->nAtoms;
	int i,k;
	double Ep=0, Em=0;

	for(i=0;i<nAtoms;i++)
	for(k=0;k<3;k++)
	{
		m->atoms[i].coordinates[k] += dx;
		setCharges(forceField);
		Ep = calculateEnergyCoulombAmber(forceField,m);

		m->atoms[i].coordinates[k] -= 2*dx;
		setCharges(forceField);
		Em = calculateEnergyCoulombAmber(forceField,m);

		m->atoms[i].gradient[k] += (Ep-Em)/dx/2;
		m->atoms[i].coordinates[k] += dx;
	}
	setCharges(forceField);
	m->potentialEnergy += calculateEnergyCoulombAmber(forceField,m);
}
/**********************************************************************/
void calculateGradientCoulombAmber(ForceField* forceField)
{
	int i;
	double energy = 0;

	Molecule* m = &forceField->molecule;
	double* coulombTerms[COULOMBDIM];
	int numberOfCoulomb = forceField->numberOfCoulomb;

	if(!strcmp(forceField->options.chargesType,"EEM") || 
	!strcmp(forceField->options.chargesType,"ACKS2"))
	{
		calculateGradientCoulombAmberEEMACKS2(forceField);
		return;
	}

	for(i=0;i<COULOMBDIM;i++) coulombTerms[i] = forceField->coulombTerms[i];

	/* non-bonded part */
#ifdef ENABLE_OMP
#pragma omp parallel for private(i,ai,aj) 
#endif
	for (  i = 0; i < numberOfCoulomb; i++ )
	{
		double rijx, rijy, rijz;
		double gradix, gradiy, gradiz;
		double rij2, rij;
		double rij3;
		double  term3;
		int ai, aj;
		double C;

		ai       = (int)coulombTerms[0][i];
		aj       = (int)coulombTerms[1][i];
		C = coulombTerms[2][i]*m->atoms[ai].charge*m->atoms[aj].charge;

		rijx =  m->atoms[ai].coordinates[0] - m->atoms[aj].coordinates[0];
		rijy =  m->atoms[ai].coordinates[1] - m->atoms[aj].coordinates[1];
		rijz =  m->atoms[ai].coordinates[2] - m->atoms[aj].coordinates[2];

		rij2 = rijx * rijx + rijy * rijy + rijz * rijz;
		if ( rij2 < POSEPS*POSEPS ) rij2 = POSEPS*POSEPS;	

		rij = sqrt( rij2 );
		rij3 = 1.0/(rij2 * rij);
		
		term3 = -C * rij3;
		gradix = term3 * rijx;
		gradiy = term3 * rijy;
		gradiz = term3 * rijz;
#ifdef ENABLE_OMP
#pragma omp critical
#endif
		{
			m->atoms[ai].gradient[0] += gradix;
			m->atoms[ai].gradient[1] += gradiy;
			m->atoms[ai].gradient[2] += gradiz;

			m->atoms[aj].gradient[0] -= gradix;
			m->atoms[aj].gradient[1] -= gradiy;
			m->atoms[aj].gradient[2] -= gradiz;
		}
		energy += C / rij;
		if(rij<1) printf("Warning : rij2 = %f atom# = %d atom# = %d\n",rij,ai+1,aj+1);
	//	printf("rij2 = %f ai = %d aj = %d\n",rij2,ai+1,aj+1);
	}  
	m->potentialEnergy += energy;
}
/**********************************************************************/
static double calculateRhoSuttonChenAmber(ForceField* forceField,Molecule* molecule)
{
	int i;
	double rijx, rijy, rijz;
	double rij2, rijm,rijn;
	int ai, aj;
	double term = 0;
	double energyAtt = 0;
	double energyRep = 0;

	Molecule* m = molecule;
	double* suttonChenTerms[SUTTONCHENDIM];
	int numberOfSuttonChen = forceField->numberOfSuttonChen;

	for(i=0;i<SUTTONCHENDIM;i++)
		suttonChenTerms[i] = forceField->suttonChenTerms[i];

	for (  i = 0; i < m->nAtoms; i++ ) m->atoms[i].rho = 0;

#ifdef ENABLE_OMP
#pragma omp parallel for private(i,ai,aj) 
#endif
	for (  i = 0; i < numberOfSuttonChen; i++ )
	{
		ai       = (int)suttonChenTerms[0][i];
		aj       = (int)suttonChenTerms[1][i];

		rijx =  m->atoms[ai].coordinates[0] - m->atoms[aj].coordinates[0];
		rijy =  m->atoms[ai].coordinates[1] - m->atoms[aj].coordinates[1];
		rijz =  m->atoms[ai].coordinates[2] - m->atoms[aj].coordinates[2];

		rij2 = rijx * rijx + rijy * rijy + rijz * rijz;
		if ( rij2 < POSEPS*POSEPS ) rij2 = POSEPS*POSEPS;	

		rijn = pow( suttonChenTerms[3][i]*suttonChenTerms[3][i]/rij2, suttonChenTerms[5][i]/2);
		rijm = pow( suttonChenTerms[3][i]*suttonChenTerms[3][i]/rij2, suttonChenTerms[6][i]/2);
		//printf("ai=%d aj=%d rijm = %f\n",ai,aj,rijm);
		
		term =  suttonChenTerms[4][i]* suttonChenTerms[4][i]*suttonChenTerms[2][i]*suttonChenTerms[2][i]*rijm;
		m->atoms[ai].rho += term;
		m->atoms[aj].rho += term;
		energyRep += suttonChenTerms[2][i]* rijn;
	}  
	for (  i = 0; i < m->nAtoms; i++ ) 
	{
		if(m->atoms[i].rho>0)   
		{
			m->atoms[i].rho = sqrt(m->atoms[i].rho);
			energyAtt +=  m->atoms[i].rho; 
			m->atoms[i].rho = 1/m->atoms[i].rho;
		}
		else  
		{
			m->atoms[i].rho  = 0;
		}
		//printf("i = %d rho = %f\n",i,m->atoms[i].rho);
	}
	/*
	printf("eneryRep = %f\n",energyRep);
	printf("eneryAtt = %f\n",-energyAtt);
	printf("eneryAll = %f\n",energyRep-energyAtt);
	*/
	//printf("enery ev par atom = %f\n",(energyRep-energyAtt)/ m->nAtoms*0.04336410);
	return energyRep - energyAtt;
}
/**********************************************************************/
void calculateGradientSuttonChenAmber(ForceField* forceField)
{
	int i;
	double rijx, rijy, rijz;
	double gradix, gradiy, gradiz;
	double rij2, rijn, rijm;
	double  term;
	int ai, aj;
	double energy = 0;

	Molecule* m = &forceField->molecule;
	double* suttonChenTerms[SUTTONCHENDIM];
	int numberOfSuttonChen = forceField->numberOfSuttonChen;

	for(i=0;i<SUTTONCHENDIM;i++)
		suttonChenTerms[i] = forceField->suttonChenTerms[i];

	energy = calculateRhoSuttonChenAmber(forceField,m);
#ifdef ENABLE_OMP
#pragma omp parallel for private(i,ai,aj) 
#endif
	for (  i = 0; i < numberOfSuttonChen; i++ )
	{
		ai       = (int)suttonChenTerms[0][i];
		aj       = (int)suttonChenTerms[1][i];

		rijx =  m->atoms[ai].coordinates[0] - m->atoms[aj].coordinates[0];
		rijy =  m->atoms[ai].coordinates[1] - m->atoms[aj].coordinates[1];
		rijz =  m->atoms[ai].coordinates[2] - m->atoms[aj].coordinates[2];

		rij2 = rijx * rijx + rijy * rijy + rijz * rijz;
		if ( rij2 < POSEPS*POSEPS ) rij2 = POSEPS*POSEPS;	

		rijn = pow( suttonChenTerms[3][i]*suttonChenTerms[3][i]/rij2, suttonChenTerms[5][i]/2);
		rijm = pow( suttonChenTerms[3][i]*suttonChenTerms[3][i]/rij2, suttonChenTerms[6][i]/2);
		
		//printf("rijn = %f rijm = %f\n", rijn,rijm);
		term =  suttonChenTerms[2][i]/rij2*
			(suttonChenTerms[5][i]*rijn-suttonChenTerms[4][i]*suttonChenTerms[4][i]*suttonChenTerms[2][i]*suttonChenTerms[6][i]/2*rijm*(m->atoms[ai].rho+m->atoms[aj].rho));
		//	(suttonChenTerms[5][i]*rijn);
		term = - term;
		gradix = term * rijx;
		gradiy = term * rijy;
		gradiz = term * rijz;
		//printf("ai = %d force=%f %f %f\n", ai, gradix, gradiy, gradiz);
#ifdef ENABLE_OMP
#pragma omp critical
#endif
		{
			m->atoms[ai].gradient[0] += gradix;
			m->atoms[ai].gradient[1] += gradiy;
			m->atoms[ai].gradient[2] += gradiz;

			m->atoms[aj].gradient[0] -= gradix;
			m->atoms[aj].gradient[1] -= gradiy;
			m->atoms[aj].gradient[2] -= gradiz;
		}
	}  
	//printf("Energy Dans grad = %f\n",energy);
	m->potentialEnergy += energy;
}
/**********************************************************************/
static void addDerivCosBeta(ForceField* forceField, int ai, int aj, int ak, double factor)
{
	Molecule* m = &forceField->molecule;

	double cosBeta, sinBeta;
	double delta = 1e-13;

	double rijx, rijy, rijz;
	double rkjx, rkjy, rkjz;
	//double rij, rkj;
	//double rij;
	double rij2, rkj2;
	double xp,yp,zp,rp;

	double term1, term2;
	double termAngleix, termAngleiy, termAngleiz;
	double termAnglekx, termAngleky, termAnglekz;
	double dot;


	rijx = m->atoms[ai].coordinates[0] - m->atoms[aj].coordinates[0];
	rijy = m->atoms[ai].coordinates[1] - m->atoms[aj].coordinates[1];
	rijz = m->atoms[ai].coordinates[2] - m->atoms[aj].coordinates[2];

	rkjx = m->atoms[ak].coordinates[0] - m->atoms[aj].coordinates[0];
	rkjy = m->atoms[ak].coordinates[1] - m->atoms[aj].coordinates[1];
	rkjz = m->atoms[ak].coordinates[2] - m->atoms[aj].coordinates[2];

	rij2 = rijx * rijx + rijy * rijy + rijz * rijz;
	//rij = sqrt(rij2);

	rkj2 = rkjx * rkjx + rkjy * rkjy + rkjz * rkjz;
	//rkj = sqrt(rkj2);

	dot = rijx*rkjx + rijy*rkjy + rijz*rkjz;
	cosBeta = dot / sqrt(rij2*rkj2);
	sinBeta = sqrt(fabs(1.0-cosBeta*cosBeta));


	xp = rkjy*rijz - rkjz*rijy;
	yp = rkjz*rijx - rkjx*rijz;
	zp = rkjx*rijy - rkjy*rijx;
	rp = sqrt(xp*xp + yp*yp + zp*zp);
	if ( rp < delta )
	{
		fprintf(forceField->logfile, "cut rp\n");
		fflush(forceField->logfile);
		rp = delta;
	}
       	term1 = sinBeta/rij2/rp*factor;
       	term2 =-sinBeta/rkj2/rp*factor;

        termAngleix = term1 * (rijy*zp-rijz*yp);
        termAngleiy = term1 * (rijz*xp-rijx*zp);
        termAngleiz = term1 * (rijx*yp-rijy*xp);

      	termAnglekx = term2 * (rkjy*zp-rkjz*yp);
      	termAngleky = term2 * (rkjz*xp-rkjx*zp);
      	termAnglekz = term2 * (rkjx*yp-rkjy*xp);

	m->atoms[ai].gradient[0] += termAngleix;
	m->atoms[ai].gradient[1] += termAngleiy;
	m->atoms[ai].gradient[2] += termAngleiz;

	m->atoms[ak].gradient[0] += termAnglekx;
	m->atoms[ak].gradient[1] += termAngleky;
	m->atoms[ak].gradient[2] += termAnglekz;

	m->atoms[aj].gradient[0] -= termAngleix + termAnglekx;
	m->atoms[aj].gradient[1] -= termAngleiy + termAngleky;
	m->atoms[aj].gradient[2] -= termAngleiz + termAnglekz;
}
/**********************************************************************/
static double getCosBeta(Molecule* m, int ai, int ax, int aj)
{
	double cosBeta = 1.0;
	if(ai>=0 && ax>=0 && aj>=0)
	{
		double rix =  m->atoms[ai].coordinates[0] - m->atoms[ax].coordinates[0];
		double riy =  m->atoms[ai].coordinates[1] - m->atoms[ax].coordinates[1];
		double riz =  m->atoms[ai].coordinates[2] - m->atoms[ax].coordinates[2];
		double ri2 = rix*rix+riy*riy+riz*riz;
		double rjx =  m->atoms[aj].coordinates[0] - m->atoms[ax].coordinates[0];
		double rjy =  m->atoms[aj].coordinates[1] - m->atoms[ax].coordinates[1];
		double rjz =  m->atoms[aj].coordinates[2] - m->atoms[ax].coordinates[2];
		double rj2 = rjx*rjx+rjy*rjy+rjz*rjz;
		cosBeta = (rix*rjx + riy*rjy + riz*rjz) / sqrt(ri2*rj2);
		if(cosBeta<-1) cosBeta = -1.0;
		if(cosBeta>1) cosBeta = 1.0;
	}
	return cosBeta;
}
/*********************************************************************/
static void calculateGradientHydrogenBonded612Amber(ForceField* forceField)
{
	int i;
	double energy = 0;

	Molecule* m = &forceField->molecule;
	double* hydrogenBondedTerms[HYDROGENBONDEDDIM];
	int numberOfHydrogenBonded =  forceField->numberOfHydrogenBonded;

	for(i=0;i<HYDROGENBONDEDDIM;i++)
		hydrogenBondedTerms[i] = forceField->hydrogenBondedTerms[i];

	/* Hydrogen-bonded part */
#ifdef ENABLE_OMP
#pragma omp parallel for private(i) 
#endif
	for (  i = 0; i < numberOfHydrogenBonded; i++ )
	{
		int ai, aj,ax;
		double e12, e6;

		double rijx, rijy, rijz;

		double gradix, gradiy, gradiz;
		double gradjx, gradjy, gradjz;
		double cosBeta = 1.0;

		double Cij, Dij, rij2,  rij4, rij6, rij8, rij12, rij14;
		double  term1, term2, term3;

		ai = (int)hydrogenBondedTerms[0][i];
		aj = (int)hydrogenBondedTerms[1][i];
		ax = (int)hydrogenBondedTerms[2][i];
		Cij = hydrogenBondedTerms[3][i];
		Dij = hydrogenBondedTerms[4][i];

		cosBeta = getCosBeta(m, ai, ax, aj);

		rijx =  m->atoms[ai].coordinates[0] - m->atoms[aj].coordinates[0];
		rijy =  m->atoms[ai].coordinates[1] - m->atoms[aj].coordinates[1];
		rijz =  m->atoms[ai].coordinates[2] - m->atoms[aj].coordinates[2];

		rij2 = rijx * rijx + rijy * rijy + rijz * rijz;

		if ( rij2 < POSEPS*POSEPS ) rij2 = POSEPS*POSEPS;	

		rij4 = rij2 * rij2;
		rij6 = rij4 * rij2;
		rij8 = rij4 * rij4;
		rij12 = rij8 * rij4;
		rij14 = rij12 * rij2;
		term1 = -12.0*Cij / rij14;
		term2 = +6*Dij / rij8*cosBeta;

		term3 = term1 + term2;
		gradix = term3 * rijx;
		gradiy = term3 * rijy;
		gradiz = term3 * rijz;
		gradjx = -gradix;
		gradjy = -gradiy;
		gradjz = -gradiz;
#ifdef ENABLE_OMP
#pragma omp critical
#endif
		{
			m->atoms[ai].gradient[0] += gradix;
			m->atoms[ai].gradient[1] += gradiy;
			m->atoms[ai].gradient[2] += gradiz;

			m->atoms[aj].gradient[0] += gradjx;
			m->atoms[aj].gradient[1] += gradjy;
			m->atoms[aj].gradient[2] += gradjz;
		}
		e12 =  Cij / rij12;
		e6 =  -Dij / (rij6);
		energy += e12 + cosBeta*e6;
		if(ax>-1) addDerivCosBeta(forceField, ai, ax, aj, e6);
	}
	m->potentialEnergy += energy;
}
/*********************************************************************/
static void calculateGradientHydrogenBonded1012Amber(ForceField* forceField)
{
	int i;
	double energy = 0;

	Molecule* m = &forceField->molecule;
	double* hydrogenBondedTerms[HYDROGENBONDEDDIM];
	int numberOfHydrogenBonded =  forceField->numberOfHydrogenBonded;

	for(i=0;i<HYDROGENBONDEDDIM;i++)
		hydrogenBondedTerms[i] = forceField->hydrogenBondedTerms[i];

	/* Hydrogen-bonded part */
#ifdef ENABLE_OMP
#pragma omp parallel for private(i) 
#endif
	for (  i = 0; i < numberOfHydrogenBonded; i++ )
	{
		int ai, aj,ax;
		double e12, e10;

		double rijx, rijy, rijz;

		double gradix, gradiy, gradiz;
		double gradjx, gradjy, gradjz;
		double cosBeta = 1.0;

		//double Cij, Dij, rij2,  rij4, rij6, rij8, rij12, rij14;
		double Cij, Dij, rij2,  rij4, rij8, rij12, rij14;
		double  term1, term2, term3;

		ai = (int)hydrogenBondedTerms[0][i];
		aj = (int)hydrogenBondedTerms[1][i];
		ax = (int)hydrogenBondedTerms[2][i];
		Cij = hydrogenBondedTerms[3][i];
		Dij = hydrogenBondedTerms[4][i];
		cosBeta = getCosBeta(m, ai, ax, aj);

		rijx =  m->atoms[ai].coordinates[0] - m->atoms[aj].coordinates[0];
		rijy =  m->atoms[ai].coordinates[1] - m->atoms[aj].coordinates[1];
		rijz =  m->atoms[ai].coordinates[2] - m->atoms[aj].coordinates[2];

		rij2 = rijx * rijx + rijy * rijy + rijz * rijz;

		if ( rij2 < POSEPS*POSEPS ) rij2 = POSEPS*POSEPS;	

		rij4 = rij2 * rij2;
		//rij6 = rij4 * rij2;
		rij8 = rij4 * rij4;
		rij12 = rij8 * rij4;
		rij14 = rij12 * rij2;
		term1 = -12.0*Cij / rij14;
		term2 = +10.0*Dij / rij12*cosBeta;

		term3 = term1 + term2;
		gradix = term3 * rijx;
		gradiy = term3 * rijy;
		gradiz = term3 * rijz;
		gradjx = -gradix;
		gradjy = -gradiy;
		gradjz = -gradiz;
#ifdef ENABLE_OMP
#pragma omp critical
#endif
		{
			m->atoms[ai].gradient[0] += gradix;
			m->atoms[ai].gradient[1] += gradiy;
			m->atoms[ai].gradient[2] += gradiz;

			m->atoms[aj].gradient[0] += gradjx;
			m->atoms[aj].gradient[1] += gradjy;
			m->atoms[aj].gradient[2] += gradjz;
		}
		e12 =  Cij / rij12;
		e10 =  -Dij / (rij8*rij2);
		energy += e12 + cosBeta*e10;
		if(ax>-1) addDerivCosBeta(forceField, ai, ax, aj, e10);
	}
	m->potentialEnergy += energy;
}
/*********************************************************************/
static void calculateGradientHydrogenBondedMorseAmber(ForceField* forceField)
{
	int i;
	double energy = 0;

	Molecule* m = &forceField->molecule;
	double* hydrogenBondedTerms[HYDROGENBONDEDDIM];
	int numberOfHydrogenBonded =  forceField->numberOfHydrogenBonded;

	for(i=0;i<HYDROGENBONDEDDIM;i++) hydrogenBondedTerms[i] = forceField->hydrogenBondedTerms[i];

	/* Hydrogen-bonded part */
#ifdef ENABLE_OMP
#pragma omp parallel for private(i) 
#endif
	for (  i = 0; i < numberOfHydrogenBonded; i++ )
	{
		int ai, aj,ax;
		double e;

		double rijx, rijy, rijz;

		double gradix, gradiy, gradiz;
		double gradjx, gradjy, gradjz;
		double cosBeta = 1.0;

		double alpha, Re, De, rij,  rij2;
		double  X, term2, term3;

		ai = (int)hydrogenBondedTerms[0][i];
		aj = (int)hydrogenBondedTerms[1][i];
		ax = (int)hydrogenBondedTerms[2][i];
		alpha = hydrogenBondedTerms[3][i];
		Re = hydrogenBondedTerms[4][i];
		De = hydrogenBondedTerms[5][i];
		cosBeta = getCosBeta(m, ai, ax, aj);

		rijx =  m->atoms[ai].coordinates[0] - m->atoms[aj].coordinates[0];
		rijy =  m->atoms[ai].coordinates[1] - m->atoms[aj].coordinates[1];
		rijz =  m->atoms[ai].coordinates[2] - m->atoms[aj].coordinates[2];

		rij2 = rijx * rijx + rijy * rijy + rijz * rijz;

		if ( rij2 < POSEPS*POSEPS ) rij2 = POSEPS*POSEPS;	

		rij = sqrt(rij2);
		X = exp(-alpha*(rij-Re));
		term2 = 1-X;
		term3 = 2*De*term2*X*alpha/rij*cosBeta;

		gradix = term3 * rijx;
		gradiy = term3 * rijy;
		gradiz = term3 * rijz;
		gradjx = -gradix;
		gradjy = -gradiy;
		gradjz = -gradiz;
#ifdef ENABLE_OMP
#pragma omp critical
#endif
		{
			m->atoms[ai].gradient[0] += gradix;
			m->atoms[ai].gradient[1] += gradiy;
			m->atoms[ai].gradient[2] += gradiz;

			m->atoms[aj].gradient[0] += gradjx;
			m->atoms[aj].gradient[1] += gradjy;
			m->atoms[aj].gradient[2] += gradjz;
		}
		e = De*(X*X-2*X);
		energy += cosBeta*e;
		if(ax>-1) addDerivCosBeta(forceField, ai, ax, aj, e);
	}
	m->potentialEnergy += energy;
}
/*********************************************************************/
static void calculateGradientHydrogenBondedAmber(ForceField* forceField)
{
	if(forceField->options.hydrogenBonded612) calculateGradientHydrogenBonded612Amber(forceField);
	else if(forceField->options.hydrogenBonded1012) calculateGradientHydrogenBonded1012Amber(forceField);
	else if(forceField->options.hydrogenBondedMorse) calculateGradientHydrogenBondedMorseAmber(forceField);
}
/**********************************************************************/
static void calculateGradientPairWise(ForceField* forceField)
{
	int i;
	int ai, aj;
	double energy = 0;

	double rijx, rijy, rijz;

	double gradix, gradiy, gradiz;
	double gradjx, gradjy, gradjz;

	double permittivityScale = 1, permittivity = 1;
	double coulombFactor;
	double rij2, rij;
	double rij3,rij5;
	double coulombTerm;
	double rij6, rij7, rij8, rij9, rij10, rij11, rij12;
	double  term1, term4, term6, term8, term10, termAll;
	double A, Beta, C4, C6, C8, C10,b;
	double s, sp, fact, br, brk, ebr;
	int n, k;
	double  B6, B8, B10;

	boolean useCoulomb = forceField->options.coulomb;
	boolean useVanderWals = forceField->options.vdw612;
	Molecule* m = &forceField->molecule;
	double* pairWiseTerms[PAIRWISEDIM];
	int numberOfPairWise = forceField->numberOfPairWise;

	for(i=0;i<PAIRWISEDIM;i++)
		pairWiseTerms[i] = forceField->pairWiseTerms[i];

	/* non-bonded part */
	coulombFactor = 332.05382 / ( permittivity * permittivityScale );
	B6 = 0;
	B8 = 0;
	B10 = 0;
	for (  i = 0; i < numberOfPairWise; i++ )
	{
		ai       = (int)pairWiseTerms[0][i];
		aj       = (int)pairWiseTerms[1][i];
		A        = pairWiseTerms[2][i];
		Beta     = pairWiseTerms[3][i];
		C4       = pairWiseTerms[4][i];
		C6       = pairWiseTerms[5][i];
		C8       = pairWiseTerms[6][i];
		C10      = pairWiseTerms[7][i];
		b        = pairWiseTerms[8][i];

		rijx =  m->atoms[ai].coordinates[0] - m->atoms[aj].coordinates[0];
		rijy =  m->atoms[ai].coordinates[1] - m->atoms[aj].coordinates[1];
		rijz =  m->atoms[ai].coordinates[2] - m->atoms[aj].coordinates[2];

		rij2 = rijx * rijx + rijy * rijy + rijz * rijz;
		if ( rij2 < POSEPS*POSEPS ) rij2 = POSEPS*POSEPS;	

		rij = sqrt( rij2 );
		rij3 = rij2 * rij;
		rij5 = rij3 * rij2;
		rij6 = rij3 * rij3;
		rij7 = rij6 * rij;
		rij8 = rij7 * rij;
		rij9 = rij8 * rij;
		rij10 = rij9 * rij;
		rij11 = rij10 * rij;
		rij12 = rij11 * rij;
		if(useCoulomb)
			coulombTerm = ( m->atoms[ai].charge*m->atoms[aj].charge * coulombFactor ) / rij3;
		else
			coulombTerm = 0.0;
		
		/*term1 = -A*Beta/rij*exp(-Beta*rij);*/
		term1 = A*Beta/rij*exp(-Beta*rij);

		br = b*rij;
		ebr = exp(-b*rij);

		term4 =   0.0;
		if(useVanderWals && fabs(C4)>1e-12)
		{
			fact = 1.0;
			s = 1.0;
			n = 2;
			brk = 1.0;
			for(k=1;k<2*n;k++)
			{
				fact *= k;
				brk *= br;
				s += brk/fact;
			}
			sp = s*b;
			fact *=2*n;
			brk *= br;
			s += brk/fact;
			term4 =   b*C4*ebr*s/rij5
				-(2*n)*C4*(1-ebr*s)/rij6
				-C4*ebr/rij5*sp;
		}

		term6 =   0.0;
		if(useVanderWals && fabs(C6)>1e-12)
		{
			fact = 1.0;
			s = 1.0;
			n = 3;
			brk = 1.0;
			for(k=1;k<2*n;k++)
			{
				fact *= k;
				brk *= br;
				s += brk/fact;
			}
			sp = s*b;
			fact *=2*n;
			brk *= br;
			s += brk/fact;
			term6 =   b*C6*ebr*s/rij7
				-(2*n)*C6*(1-ebr*s)/rij8
				-C6*ebr/rij7*sp;
		}
		term8 =   0.0;
		if(useVanderWals && fabs(C8)>1e-12)
		{
			fact = 1.0;
			s = 1.0;
			n = 4;
			brk = 1.0;
			for(k=1;k<2*n;k++)
			{
				fact *= k;
				brk *= br;
				s += brk/fact;
			}
			sp = s*b;
			fact *=2*n;
			brk *= br;
			s += brk/fact;
			term8 =   b*C8*ebr*s/rij9
				-(2*n)*C8*(1-ebr*s)/rij10
				-C8*ebr/rij9*sp;
		}

		term10 =   0.0;
		if(useVanderWals && fabs(C10)>1e-12)
		{
			fact = 1.0;
			s = 1.0;
			n = 5;
			brk = 1.0;
			for(k=1;k<2*n;k++)
			{
				fact *= k;
				brk *= br;
				s += brk/fact;
			}
			sp = s*b;

			fact *=2*n;
			brk *= br;
			s += brk/fact;
			term10 =   b*C10*ebr*s/rij11
				-(2*n)*C10*(1-ebr*s)/rij12
				-C10*ebr/rij11*sp;
		}

		//termAll = term1 - term6 - term8 - term10 + coulombTerm;
		termAll = term1 + term4 + term6 + term8 + term10 + coulombTerm;
		termAll = - termAll;


		gradix = termAll * rijx;
		gradiy = termAll * rijy;
		gradiz = termAll * rijz;
		gradjx = - gradix;
		gradjy = - gradiy;
		gradjz = - gradiz;
#ifdef ENABLE_OMP
#pragma omp critical
#endif
		{
			m->atoms[ai].gradient[0] += gradix;
			m->atoms[ai].gradient[1] += gradiy;
			m->atoms[ai].gradient[2] += gradiz;

			m->atoms[aj].gradient[0] += gradjx;
			m->atoms[aj].gradient[1] += gradjy;
			m->atoms[aj].gradient[2] += gradjz;
		}
		if(useCoulomb) coulombTerm = ( m->atoms[ai].charge*m->atoms[aj].charge * coulombFactor ) / rij;
		else coulombTerm = 0.0;

		B6  = 0;
		B8  = 0;
		B10 = 0;
		/* printf("A = %f Beta = %f qi = %f qj = %f rij = %f\n",A,Beta,chargei,chargej,rij);*/
		if(useVanderWals)
		{
			double fact = 1.0;
			double s = 1.0;
			double br = b*rij;
			double brk = 1.0;
			int k;

			if(fabs(C6)>1e-12)
			{
				for(k=1;k<=2*3;k++)
				{
					fact *= k;
					brk *= br;
					s += brk/fact;
				}
				B6 = C6*(1-exp(-br)*s);
			}

			if(fabs(C8)>1e-12)
			{
				fact = 1.0;
				s = 1.0;
				br = b*rij;
				brk = 1.0;
				for(k=1;k<=2*4;k++)
				{
					fact *= k;
					brk *= br;
					s += brk/fact;
				}
				B8 = C8*(1-exp(-br)*s);
			}

			if(fabs(C10)>1e-12)
			{
				fact = 1.0;
				s = 1.0;
				br = b*rij;
				brk = 1.0;
				for(k=1;k<=2*5;k++)
				{
					fact *= k;
					brk *= br;
					s += brk/fact;
				}
				B10 = C10*(1-exp(-br)*s);
			}
		}

		energy += A*exp(-Beta*rij)
			- B6 / rij6 
			- B8 / rij8 
			- B10 / rij10 
			+ coulombTerm;
	}  
	m->potentialEnergy += energy;
}
#endif /* ENABLE_CL */
/**********************************************************************/
static void calculateGradientAmberAnalytic(ForceField* forceField)
{
#ifdef DEBUG
        TimerType timer;
        TimerType timer2;
#endif
#ifdef ENABLE_CL
	CLProp clProp = getCLProp();
	//size_t local = 64;
	size_t global;
	int i;
	int j;
	cl_int err;

#ifdef DEBUG
        timer_init(timer);
       	timer_start( timer );
#endif
	forceField->molecule.potentialEnergy = 0;
#ifdef DEBUG
	fprintf(forceField->logfile, "Begin calculateGradientAmber\n");
#endif
/*
	global = forceField->nMaxTerms;
	clEnqueueNDRangeKernel(clProp.command_queue, forceField->initVelocitiesKernel, 1, NULL, &global, NULL, 0, NULL, NULL);
	clFinish(clProp.command_queue);
*/
//	if(forceField->numberOfVdw612<forceField->nMaxTerms)
	{
#ifdef DEBUG
        timer_init(timer2);
       	timer_start( timer2 );
#endif
	global = forceField->nMaxTerms;
	err = clEnqueueNDRangeKernel(clProp.command_queue, forceField->initEnergyKernel, 1, NULL, &global, NULL, 0, NULL, NULL);
	if(err != CL_SUCCESS)
	{
		printErrorCLRun(err);
        	fprintf(forceField->logfile, "I cannot execute initEnergyKernel\n");
		fflush(forceField->logfile);
		exit(1);
	}
	clFinish(clProp.command_queue);
#ifdef DEBUG
       	timer_stop(timer2);
        fprintf(forceField->logfile, "time initEnergy (s) = %f\n", timer_get(timer2)*1e-6);
	fflush(forceField->logfile);

        timer_init(timer2);
       	timer_start( timer2 );
#endif
	//global = forceField->nBlockGradientBuffer*forceField->molecule.nAtoms;
	global = forceField->molecule.nAtoms;
	global = forceField->nBlockGradientBuffer;
	//global = forceField->nMaxTerms;
	global = forceField->nBlockGradientBuffer*forceField->molecule.nAtoms;
	//global = 512;
	err = clEnqueueNDRangeKernel(clProp.command_queue, forceField->initGradientsKernel, 1, NULL, &global, NULL, 0, NULL, NULL);
	if(err != CL_SUCCESS)
	{
		printErrorCLRun(err);
        	fprintf(forceField->logfile, "I cannot execute initGradientsKernel\n");
		fflush(forceField->logfile);
		exit(1);
	}
	clFinish(clProp.command_queue);
#ifdef DEBUG
       	timer_stop(timer2);
        fprintf(forceField->logfile, "time initGrad (s) = %f\n", timer_get(timer2)*1e-6);
	fflush(forceField->logfile);
#endif
	}
	if(forceField->numberOfVdw612>0)
	{
#ifdef DEBUG
        timer_init(timer2);
       	timer_start( timer2 );
#endif
	//size_t local = 128;
	//global = ((forceField->numberOfVdw612+local-1)/local)*local;
	//clEnqueueNDRangeKernel(clProp.command_queue, forceField->addGradientVdw612AmberKernel, 1, NULL, &global, &local, 0, NULL, NULL);
	global = forceField->numberOfVdw612;
	global = forceField->molecule.nAtoms;
	global = 64;
	global = forceField->numberOfVdw612;
	clEnqueueNDRangeKernel(clProp.command_queue, forceField->addGradientVdw612AmberKernel, 1, NULL, &global, NULL, 0, NULL, NULL);
	clFinish(clProp.command_queue);
#ifdef DEBUG
       	timer_stop(timer2);
        fprintf(forceField->logfile, "time Vdw612 (s) = %f\n", timer_get(timer2)*1e-6);
	fflush(forceField->logfile);
#endif
	}
	if(forceField->numberOfVdw714>0)
	{
#ifdef DEBUG
        timer_init(timer2);
       	timer_start( timer2 );
#endif
	//size_t local = 128;
	//global = ((forceField->numberOfVdw714+local-1)/local)*local;
	//clEnqueueNDRangeKernel(clProp.command_queue, forceField->addGradientVdw714AmberKernel, 1, NULL, &global, &local, 0, NULL, NULL);
	global = forceField->numberOfVdw714;
	global = forceField->molecule.nAtoms;
	global = 64;
	global = forceField->numberOfVdw714;
	clEnqueueNDRangeKernel(clProp.command_queue, forceField->addGradientVdw714AmberKernel, 1, NULL, &global, NULL, 0, NULL, NULL);
	clFinish(clProp.command_queue);
#ifdef DEBUG
       	timer_stop(timer2);
        fprintf(forceField->logfile, "time Vdw714 (s) = %f\n", timer_get(timer2)*1e-6);
	fflush(forceField->logfile);
#endif
	}
	if(forceField->numberOfSuttonChen>0)
	{
	global = forceField->numberOfSuttonChen;
	clEnqueueNDRangeKernel(clProp.command_queue, forceField->initRhoKernel, 1, NULL, &global, NULL, 0, NULL, NULL);
	clEnqueueNDRangeKernel(clProp.command_queue, forceField->computeRhoKernel, 1, NULL, &global, NULL, 0, NULL, NULL);
	global = forceField->molecule.nAtoms;
	clEnqueueNDRangeKernel(clProp.command_queue, forceField->reduceRhoKernel, 1, NULL, &global, NULL, 0, NULL, NULL);
	global = forceField->numberOfSuttonChen;
#ifdef DEBUG
        timer_init(timer2);
       	timer_start( timer2 );
#endif
	//size_t local = 128;
	//global = ((forceField->numberOfVdw612+local-1)/local)*local;
	//clEnqueueNDRangeKernel(clProp.command_queue, forceField->addGradientSuttonChenAmberKernel, 1, NULL, &global, &local, 0, NULL, NULL);
	global = forceField->numberOfSuttonChen;
	global = forceField->molecule.nAtoms;
	global = 64;
	global = forceField->numberOfSuttonChen;
	clEnqueueNDRangeKernel(clProp.command_queue, forceField->addGradientSuttonChenKernel, 1, NULL, &global, NULL, 0, NULL, NULL);
	clFinish(clProp.command_queue);
#ifdef DEBUG
       	timer_stop(timer2);
        fprintf(forceField->logfile, "time Vdw612 (s) = %f\n", timer_get(timer2)*1e-6);
	fflush(forceField->logfile);
#endif
	}


/*
	global = forceField->molecule.nAtoms;
	//global = forceField->nMaxTerms;
	err = clEnqueueNDRangeKernel(clProp.command_queue, forceField->initVelocitiesKernel, 1, NULL, &global, NULL, 0, NULL, NULL);
	if(err != CL_SUCCESS)
	{
		printErrorCLRun(err);
        	fprintf(forceField->logfile, "I cannot execute initVelocitiesKernel\n");
		fflush(forceField->logfile);
		exit(1);
	}
	clFinish(clProp.command_queue);
*/

	if(forceField->numberOfStretchTerms>0)
	{
#ifdef DEBUG
        timer_init(timer2);
       	timer_start( timer2 );
#endif
	//size_t local = 16;
	//global = ((forceField->numberOfStretchTerms+local-1)/local)*local;
	global = forceField->numberOfStretchTerms;
	//global = forceField->nMaxTerms;
	//fprintf(forceField->logfile, "Call addGradientBondAmberKernel \n");
	//fflush(forceField->logfile);
	err = clEnqueueNDRangeKernel(clProp.command_queue, forceField->addGradientBondAmberKernel, 1, NULL, &global, NULL, 0, NULL, NULL);
	if(err != CL_SUCCESS)
	{
		printErrorCLRun(err);
        	fprintf(forceField->logfile, "I cannot execute addGradientBondAmberKernel\n");
		fflush(forceField->logfile);
		exit(1);
	}
	clFinish(clProp.command_queue);
#ifdef DEBUG
       	timer_stop(timer2);
        fprintf(forceField->logfile, "time Bond (s) = %f\n", timer_get(timer2)*1e-6);
	fflush(forceField->logfile);
#endif
	}

	if(forceField->numberOfBendTerms>0)
	{
#ifdef DEBUG
        timer_init(timer2);
       	timer_start( timer2 );
#endif
	global = forceField->numberOfBendTerms;
	//global = forceField->nMaxTerms;
	clEnqueueNDRangeKernel(clProp.command_queue, forceField->addGradientBendAmberKernel, 1, NULL, &global, NULL, 0, NULL, NULL);
	clFinish(clProp.command_queue);
#ifdef DEBUG
       	timer_stop(timer2);
        fprintf(forceField->logfile, "time Bend (s) = %f\n", timer_get(timer2)*1e-6);
	fflush(forceField->logfile);
#endif
	}

	if(forceField->numberOfDihedralTerms>0)
	{
#ifdef DEBUG
        timer_init(timer2);
       	timer_start( timer2 );
#endif
	global = forceField->numberOfDihedralTerms;
	//global = forceField->nMaxTerms;
	clEnqueueNDRangeKernel(clProp.command_queue, forceField->addGradientDihedralAmberKernel, 1, NULL, &global, NULL, 0, NULL, NULL);
	clFinish(clProp.command_queue);
#ifdef DEBUG
       	timer_stop(timer2);
        fprintf(forceField->logfile, "time Dihedral (s) = %f\n", timer_get(timer2)*1e-6);
	fflush(forceField->logfile);
#endif
	}


	if(forceField->numberOfHydrogenBonded>0)
	{
#ifdef DEBUG
        timer_init(timer2);
       	timer_start( timer2 );
#endif
	global = forceField->numberOfHydrogenBonded;
	//global = forceField->nMaxTerms;
	clEnqueueNDRangeKernel(clProp.command_queue, forceField->addGradientHydrogenBondedAmberKernel, 1, NULL, &global, NULL, 0, NULL, NULL);
	clFinish(clProp.command_queue);
#ifdef DEBUG
       	timer_stop(timer2);
        fprintf(forceField->logfile, "time HydrogenBonded (s) = %f\n", timer_get(timer2)*1e-6);
	fflush(forceField->logfile);
#endif
	}

	if(forceField->numberOfPairWise>0)
	{
	global = forceField->numberOfPairWise;
	//global = forceField->nMaxTerms;
	clEnqueueNDRangeKernel(clProp.command_queue, forceField->addGradientPairWiseKernel, 1, NULL, &global, NULL, 0, NULL, NULL);
	clFinish(clProp.command_queue);
	}

#ifdef DEBUG
        timer_init(timer2);
       	timer_start( timer2 );
#endif
	global = forceField->molecule.nAtoms;
	global = 512;
	clEnqueueNDRangeKernel(clProp.command_queue, forceField->reduceGradientsKernel, 1, NULL, &global, NULL, 0, NULL, NULL);
	clFinish(clProp.command_queue);
#ifdef DEBUG
       	timer_stop(timer2);
        fprintf(forceField->logfile, "time reduceGrad (s) = %f\n", timer_get(timer2)*1e-6);
	fflush(forceField->logfile);

        timer_init(timer2);
       	timer_start( timer2 );
#endif
	clEnqueueReadBuffer(clProp.command_queue, forceField->gradientBufferCL, CL_TRUE, 0, sizeof(cl_float4)*forceField->molecule.nAtoms, forceField->gradientBufferCPU, 0, NULL, NULL);
	for( i=0; i<forceField->molecule.nAtoms;i++)
	{
		if(!forceField->molecule.atoms[i].variable)
			for(j=0;j<3;j++)
				forceField->molecule.atoms[i].gradient[j] = 0.0;
		else
			for(j=0;j<3;j++)
				forceField->molecule.atoms[i].gradient[j] = forceField->gradientBufferCPU[i].s[j];
	}
	clEnqueueReadBuffer(clProp.command_queue, forceField->energyBufferCL, CL_TRUE, 0, sizeof(cl_float)*forceField->nMaxTerms, forceField->energyBufferCPU, 0, NULL, NULL);
	forceField->molecule.potentialEnergy = 0;
	for(i=0;i<forceField->nMaxTerms;i++) forceField->molecule.potentialEnergy+=forceField->energyBufferCPU[i];
#ifdef DEBUG
       	timer_stop(timer2);
        fprintf(forceField->logfile, "time read buffer (s) = %f\n", timer_get(timer2)*1e-6);
	fflush(forceField->logfile);
       	timer_stop(timer);
        fprintf(forceField->logfile, "time (s) = %f\n", timer_get(timer)*1e-6);
	fflush(forceField->logfile);
#endif
#else // ENABLE_CL
	int i;
	int j;
	Molecule* m = &forceField->molecule;

	forceField->molecule.potentialEnergy = 0;

#ifdef DEBUG
        timer_init(timer);
       	timer_start( timer );
#endif

	for(j=0;j<3;j++)
		for( i=0; i<m->nAtoms;i++)
			m->atoms[i].gradient[j] = 0.0;

	if(!strstr(forceField->options.chargesType,"BEGIN")) setCharges(forceField);
	calculateGradientBondAmber(forceField);
	calculateGradientBendAmber(forceField);
	calculateGradientStrBendAmber(forceField);
#ifdef DEBUG
        timer_init(timer2);
       	timer_start( timer2 );
#endif
	calculateGradientDihedralAmber(forceField);
#ifdef DEBUG
       	timer_stop(timer2);
        fprintf(forceField->logfile, "time Dihedral (s) = %f\n", timer_get(timer2)*1e-6);
	fflush(forceField->logfile);
#endif
	calculateGradientImproperTorsion(forceField);
#ifdef DEBUG
        timer_init(timer2);
       	timer_start( timer2 );
#endif
	calculateGradientVdw612Amber(forceField);
#ifdef DEBUG
       	timer_stop(timer2);
        fprintf(forceField->logfile, "time Vdw612 (s) = %f\n", timer_get(timer2)*1e-6);
	fflush(forceField->logfile);
#endif

#ifdef DEBUG
        timer_init(timer2);
       	timer_start( timer2 );
#endif
	calculateGradientVdw714Amber(forceField);
#ifdef DEBUG
       	timer_stop(timer2);
        fprintf(forceField->logfile, "time Vdw714 (s) = %f\n", timer_get(timer2)*1e-6);
	fflush(forceField->logfile);
#endif

#ifdef DEBUG
        timer_init(timer2);
       	timer_start( timer2 );
#endif
	calculateGradientSuttonChenAmber(forceField);
#ifdef DEBUG
       	timer_stop(timer2);
        fprintf(forceField->logfile, "time SuttonChen (s) = %f\n", timer_get(timer2)*1e-6);
	fflush(forceField->logfile);
#endif

	calculateGradientCoulombAmber(forceField);
	calculateGradientHydrogenBondedAmber(forceField);
	calculateGradientPairWise(forceField);
	addH4Correction(forceField, TRUE);
	addWallCorrection(forceField,TRUE);

	for( i=0; i<m->nAtoms;i++)
	{
		if(!m->atoms[i].variable)
			for(j=0;j<3;j++)
				m->atoms[i].gradient[j] = 0.0;
	}
#ifdef DEBUG
       	timer_stop(timer);
        fprintf(forceField->logfile, "time (s) = %f\n", timer_get(timer)*1e-6);
	fflush(forceField->logfile);
#endif
#endif // ENABLE_CL
}
/**********************************************************************/
static void calculateGradientAmber(ForceField* forceField)
{
	if(!forceField->options.numeric) calculateGradientAmberAnalytic(forceField);
	else forceField->klass->calculateGradientNumeric(forceField);
}
/**********************************************************************/
/*
static void calculateGradientNumericAmber(ForceField* forceField)
{
	int i;
	int j;
	Molecule* m = &forceField->molecule;
	double h=0.0001;
	double E1;
	double E2;

	for(j=0;j<3;j++)
		for( i=0; i<m->nAtoms;i++)
		{
			m->atoms[i].coordinates[j] += h;
			E1 = calculateEnergyTmpAmber(forceField,&m);
			m->atoms[i].coordinates[j] -= h+h;
			E2 = calculateEnergyTmpAmber(forceField,&m);
			m->atoms[i].coordinates[j] += h;
			m->atoms[i].gradient[j] = (E1-E2)/2/h;
		}


}
*/
#ifndef ENABLE_CL
/**********************************************************************/
static double calculateEnergyBondAmber(ForceField* forceField,Molecule* molecule)
{
	int i;

	Molecule* m = molecule;
	double* bondStretchTerms[STRETCHDIM];
	int numberOfStretchTerms = forceField->numberOfStretchTerms;
	double energy = 0.0;

	for( i=0; i<STRETCHDIM;i++) bondStretchTerms[i] = forceField->bondStretchTerms[i];


#ifdef ENABLE_OMP
#pragma omp parallel for private(i) reduction(+:energy)
#endif
	for (  i = 0; i < numberOfStretchTerms; i++ )
	{
		int ai, aj;
		double rijx, rijy, rijz, forceOrAlpha, equilibriumDistance;
		double h3OrDe, h4, h5, h6;
		double bondLength;
		double diff, diff2, diff3, diff4, diff5, diff6;
		int type = 0;

		type = (int)bondStretchTerms[0][i];
		ai = (int)bondStretchTerms[1][i];
		aj = (int)bondStretchTerms[2][i];
		forceOrAlpha = bondStretchTerms[3][i];
		equilibriumDistance = bondStretchTerms[4][i];
		h3OrDe = bondStretchTerms[5][i];
		h4 = bondStretchTerms[6][i];
		h5 = bondStretchTerms[7][i];
		h6 = bondStretchTerms[8][i];
		
		rijx =  m->atoms[ai].coordinates[0] - m->atoms[aj].coordinates[0];
		rijy =  m->atoms[ai].coordinates[1] - m->atoms[aj].coordinates[1];
		rijz =  m->atoms[ai].coordinates[2] - m->atoms[aj].coordinates[2];

		bondLength = sqrt( rijx * rijx + rijy * rijy + rijz * rijz );
		diff = bondLength - equilibriumDistance;

		if(type==0)
		{
			diff2 = diff*diff;
			diff3 = diff2*diff;
			diff4 = diff3*diff;
			diff5 = diff4*diff;
			diff6 = diff5*diff;
			energy += forceOrAlpha * diff2;
			energy += h3OrDe * forceOrAlpha * diff3;
			energy += h4 * forceOrAlpha * diff4;
			energy += h5 * forceOrAlpha * diff5;
			energy += h6 * forceOrAlpha * diff6;
		}
		else
		{
			double alpha = forceOrAlpha;
			double De = h3OrDe;
			double X = exp(-alpha*(bondLength-equilibriumDistance));
			double t = 1-X;
			energy += De*t*t;
		}
	} 
	return energy;
}
/**********************************************************************/
static double calculateEnergyBendAmber(ForceField* forceField,Molecule* molecule)
{
	int i;
	double energy = 0.0;
	static double D2R = 1/(RADTODEG);

	Molecule* m = molecule;
	double* angleBendTerms[BENDDIM];
	int numberOfBendTerms = forceField->numberOfBendTerms;

	for( i=0; i<BENDDIM;i++)
		angleBendTerms[i] = forceField->angleBendTerms[i]; 

#ifdef ENABLE_OMP
#pragma omp parallel for private(i) reduction(+:energy)
#endif
	for (  i = 0; i < numberOfBendTerms; i++ )
	{
//HERE
		int ai, aj, ak;
		double thetaRad;
		double cosine;
		double diff = 0;
		double diff2 = 0;
		double diff3 = 0;
		double diff4 = 0;

		double rijx, rijy, rijz;
		double rkjx, rkjy, rkjz;
		double rij2, rkj2;

		double rijDotrkj;

		ai = (int)angleBendTerms[0][i];
		aj = (int)angleBendTerms[1][i];
		ak = (int)angleBendTerms[2][i];

		rijx = m->atoms[ai].coordinates[0] - m->atoms[aj].coordinates[0];
		rijy = m->atoms[ai].coordinates[1] - m->atoms[aj].coordinates[1];
		rijz = m->atoms[ai].coordinates[2] - m->atoms[aj].coordinates[2];
		rij2 = rijx * rijx + rijy * rijy + rijz * rijz;

		rkjx = m->atoms[ak].coordinates[0] - m->atoms[aj].coordinates[0];
		rkjy = m->atoms[ak].coordinates[1] - m->atoms[aj].coordinates[1];
		rkjz = m->atoms[ak].coordinates[2] - m->atoms[aj].coordinates[2];
		rkj2 = rkjx * rkjx + rkjy * rkjy + rkjz * rkjz;

		if (rij2==0 || rkj2==0) continue;

		rijDotrkj = rijx * rkjx + rijy * rkjy + rijz * rkjz;

	        cosine = rijDotrkj / sqrt(rij2*rkj2);
		if(cosine>1) cosine = 1;
		if(cosine<-1) cosine = -1;
                thetaRad = acos(cosine);
		diff =  thetaRad - D2R*angleBendTerms[4][i];
		diff2 = diff*diff;
		diff3 = diff2*diff;
		diff4 = diff3*diff;

		energy += angleBendTerms[3][i]*diff2*
		(1.0+angleBendTerms[5][i]*diff+angleBendTerms[6][i]*diff2+angleBendTerms[7][i]*diff3+angleBendTerms[8][i]*diff4);
	} 
	//printf("Energy bend = %f\n",energy);
	return energy;
}
/**********************************************************************/
static double calculateEnergyStrBendAmber(ForceField* forceField,Molecule* molecule)
{
	int i;
	int ai, aj, ak;
	double thetaDeg;
	double term;
	double energy = 0.0;
	static double D2R = 1/(RADTODEG);
	double forceConstant12;
	double forceConstant23;
	double Re12;
	double Re23;
	double angle123;
	double R12;
	double R23;
	double rx, ry, rz;


	Molecule* m = molecule;
	double* strBendTerms[STRBENDDIM];
	int numberOfStrBendTerms = forceField->numberOfStrBendTerms;

	for( i=0; i<STRBENDDIM;i++)
		strBendTerms[i] = forceField->strBendTerms[i]; 

#ifdef ENABLE_OMP
#pragma omp parallel for private(i,ai,aj,ak,thetaDeg,term) reduction(+:energy)
#endif
	for (  i = 0; i < numberOfStrBendTerms; i++ )
	{
		ai = (int)strBendTerms[2][i];
		aj = (int)strBendTerms[3][i];
		ak = (int)strBendTerms[4][i];
		forceConstant12 = strBendTerms[5][i];
		forceConstant23 = strBendTerms[6][i];
		Re12 = strBendTerms[7][i];
		Re23 = strBendTerms[8][i];
		angle123 = strBendTerms[9][i];
		thetaDeg = getAngle(  &m->atoms[ai] ,&m->atoms[aj], &m->atoms[ak]);
	        rx = m->atoms[ai].coordinates[0] - m->atoms[aj].coordinates[0];
                ry = m->atoms[ai].coordinates[1] - m->atoms[aj].coordinates[1];
                rz = m->atoms[ai].coordinates[2] - m->atoms[aj].coordinates[2];
                R12 = sqrt( rx * rx + ry * ry + rz * rz );

	        rx = m->atoms[ak].coordinates[0] - m->atoms[aj].coordinates[0];
                ry = m->atoms[ak].coordinates[1] - m->atoms[aj].coordinates[1];
                rz = m->atoms[ak].coordinates[2] - m->atoms[aj].coordinates[2];
                R23 = sqrt( rx * rx + ry * ry + rz * rz );
		term = D2R*(thetaDeg - angle123)*(forceConstant12*(R12 - Re12)+forceConstant23*(R23 - Re23));
		energy += term;
	
		/*
		fprintf(forceField->logfile, "f =%f t0 = %f  t= %f e= %f\n",
			strBendTerms[5][i],
			strBendTerms[6][i],
			thetaDeg,
			energy);
		fflush(forceField->logfile);
		*/
	

	} 
	return energy;
}
/**********************************************************************/
double calculateEnergyDihedralAmber(ForceField* forceField,Molecule* molecule)
{
	int i;
	int ai, aj, ak, al;
	double phiDeg;
	Molecule* m = molecule;
	double* dihedralAngleTerms[DIHEDRALDIM];
	int numberOfDihedralTerms = forceField->numberOfDihedralTerms;
	double energy = 0.0;
	static double D2R = 1/(RADTODEG);

	for(i=0;i<DIHEDRALDIM;i++)
		dihedralAngleTerms[i] = forceField->dihedralAngleTerms[i];

#ifdef ENABLE_OMP
#pragma omp parallel for private(i,ai,aj,ak,al,phiDeg) reduction(+:energy)
#endif
	for (  i = 0; i < numberOfDihedralTerms; i++ )
	{
		ai = (int)dihedralAngleTerms[0][i];
		aj = (int)dihedralAngleTerms[1][i];
		ak = (int)dihedralAngleTerms[2][i];
		al = (int)dihedralAngleTerms[3][i];

		phiDeg = getTorsion(  &m->atoms[ai] ,&m->atoms[aj], &m->atoms[ak], &m->atoms[al]);

		energy += dihedralAngleTerms[5][i]/dihedralAngleTerms[4][i] * 
		( 1.0 + cos( D2R*(dihedralAngleTerms[7][i] * phiDeg - dihedralAngleTerms[6][i] )) );
	}
	return energy;
}
/**********************************************************************/
static double calculateEnergyImproperTorsionAmber(ForceField* forceField,Molecule* molecule)
{
//HERE
	int i;
	Molecule* m = &forceField->molecule;
	double forceConstant, equilibriumAngle;
	int numberOfImproperTorsionTerms = forceField->numberOfImproperTorsionTerms;
	double* improperTorsionTerms[IMPROPERDIHEDRALDIM];
	static double D2R = 1/(RADTODEG);
	double energy = 0;
	double phiDeg;

	for(i=0;i<IMPROPERDIHEDRALDIM;i++) improperTorsionTerms[i] = forceField->improperTorsionTerms[i];

#ifdef ENABLE_OMP
#pragma omp parallel for private(i) 
#endif
	for (  i = 0; i < numberOfImproperTorsionTerms; i++ )
	{
		int ai, aj, ak, al;
		int n;

		ai = (int)improperTorsionTerms[0][i];
		aj = (int)improperTorsionTerms[1][i];
		ak = (int)improperTorsionTerms[2][i];
		al = (int)improperTorsionTerms[3][i];

		forceConstant = improperTorsionTerms[4][i];
		n = (int)improperTorsionTerms[6][i];
                equilibriumAngle = improperTorsionTerms[5][i];

		phiDeg = getTorsion(  &m->atoms[ai] ,&m->atoms[aj], &m->atoms[ak], &m->atoms[al]);
		//printf("PhiDeg force equi =%f %f %f\n",phiDeg,forceConstant,equilibriumAngle);
		energy += forceConstant * ( 1.0 + cos( D2R*(n * phiDeg - equilibriumAngle )) );
	}
	//printf("ImproperTorsion energy =%f\n",energy);
	//printf("numberOfImproperTorsionTerms  =%d\n",numberOfImproperTorsionTerms);
	m->potentialEnergy += energy;
	return energy;
}
/**********************************************************************/
double calculateEnergyVdw612Amber(ForceField* forceField,Molecule* molecule)
{
	int i;
	int ai, aj;
	double rij2, rij6;
	double rijx, rijy, rijz;
	Molecule* m = molecule;
	double* vdw612Terms[VDW612DIM];
	int numberOfVdw612 = forceField->numberOfVdw612;
	double energy = 0.0;

	//printf("calculateEnergyVdw612Amber : useCoulomb=%d\n",useCoulomb);
	for(i=0;i<VDW612DIM;i++)
		vdw612Terms[i] = forceField->vdw612Terms[i];

	/* now for non-bonded term */
	/*fprintf(forceField->logfile, "number of Non Bonded terms = %d\n",numberOfVdw612);*/
	//fflush(forceField->logfile);
#ifdef ENABLE_OMP
#pragma omp parallel for private(i,ai,aj,rijx,rijy,rijz,rij2,rij6,coulombTerm) reduction(+:energy)
#endif
	for (  i = 0; i < numberOfVdw612; i++ )
	{
		ai     = (int)vdw612Terms[0][i];
		aj     = (int)vdw612Terms[1][i];

		rijx =  m->atoms[ai].coordinates[0] -  m->atoms[aj].coordinates[0];
		rijy =  m->atoms[ai].coordinates[1] -  m->atoms[aj].coordinates[1];
		rijz =  m->atoms[ai].coordinates[2] -  m->atoms[aj].coordinates[2];

		rij2 = (rijx * rijx + rijy * rijy + rijz * rijz);
		if ( rij2 < POSEPS*POSEPS ) rij2 = POSEPS*POSEPS;	
		rij2 = 1.0/rij2;
		rij6 = rij2 * rij2 * rij2;

		energy += (vdw612Terms[2][i]* rij6 - vdw612Terms[3][i]) * rij6;

		//if(1/rij2<2) printf("rij2 = %f ai = %d aj = %d\n",1/rij2,ai+1,aj+1);
		//printf("rij2 = %f ai = %d aj = %d\n",1/rij2,ai+1,aj+1);
		//printf("Coulomterm = %f\n",coulombTerm);
		/*
		fprintf(forceField->logfile, "A =%f B = %f  r= %f e= %f\n",
			Aij,Bij ,rij,energy);
		fflush(forceField->logfile);
		*/
	}  
	/* fprintf(forceField->logfile, "Non Bonded energy = %f\n",energy);*/
	//fflush(forceField->logfile);
	return energy;
}
/**********************************************************************/
double calculateEnergyVdw714Amber(ForceField* forceField,Molecule* molecule)
{
	int i;
	double energy = 0;

	Molecule* m = &forceField->molecule;
	double* vdw714Terms[VDW714DIM];
	int numberOfVdw714 = forceField->numberOfVdw714;


	for(i=0;i<VDW714DIM;i++) vdw714Terms[i] = forceField->vdw714Terms[i];

#ifdef ENABLE_OMP
#pragma omp parallel for private(i,ai,aj) 
#endif
	for (  i = 0; i < numberOfVdw714; i++ )
	{
		double rijx, rijy, rijz;
		double rij2, rij;
		int ai, aj;
		double epsilon, gamma, delta, R0;
		//double R02;
		double d,d2,d4,d7;
		//double sr,sr2,sr4,sr7,sr8;
		double sr,sr2,sr4,sr7;
		//double rho,rho2,rho4,rho5, rho7;
		double rho,rho2,rho4, rho7;
		double srijg7;
		double e7,e14;

		ai       = (int)vdw714Terms[0][i];
		aj       = (int)vdw714Terms[1][i];

		epsilon  = vdw714Terms[2][i];
		R0       = vdw714Terms[3][i];
		gamma    = vdw714Terms[4][i];
		delta    = vdw714Terms[5][i];

		rijx =  m->atoms[ai].coordinates[0] - m->atoms[aj].coordinates[0];
		rijy =  m->atoms[ai].coordinates[1] - m->atoms[aj].coordinates[1];
		rijz =  m->atoms[ai].coordinates[2] - m->atoms[aj].coordinates[2];

		rij2 = rijx * rijx + rijy * rijy + rijz * rijz;
		if ( rij2 < POSEPS*POSEPS ) rij2 = POSEPS*POSEPS;	

		rij = sqrt( rij2 );
		//R02 = R0*R0;

		rho = rij/R0;
		rho2 = rho*rho;
		rho4 = rho2*rho2;
		//rho5 = rho4*rho;
		rho7 = rho4*rho2*rho;

		d = (1+delta);
		d2 = d*d;
		d4 = d2*d2;
		d7 = d4*d2*d;

		sr = 1.0/(rho+delta);
		sr2 = sr*sr;
		sr4 = sr2*sr2;
		sr7 = sr4*sr2*sr;
		//sr8 = sr4*sr4;
		srijg7=1.0/(rho7+gamma);
		
		e7 = epsilon*d7*sr7;
		e14=e7*(1+gamma)*srijg7;
		e7 = 2*e7;

#ifdef ENABLE_OMP
#pragma omp critical
#endif
		{
		energy += e14-e7;
		}
		if(rij<1) printf("Warning : rij2 = %f atom# = %d atom# = %d\n",rij,ai+1,aj+1);
	//	printf("rij2 = %f ai = %d aj = %d\n",rij2,ai+1,aj+1);
	}  
	//printf("E in ener 7-14 %f\n",energy);

	return energy;
}
/**********************************************************************/
double calculateEnergyCoulombAmber(ForceField* forceField,Molecule* molecule)
{
	int i;
	Molecule* m = molecule;
	double* coulombTerms[COULOMBDIM];
	int numberOfCoulomb = forceField->numberOfCoulomb;
	double energy = 0.0;

/*
	if(!strcmp(forceField->options.chargesType,"EEM")) return forceField->molecule.klass->getEnergyEEM(&forceField->molecule);
	if(!strcmp(forceField->options.chargesType,"ACKS2")) return forceField->molecule.klass->getEnergyACKS2(&forceField->molecule);
*/

	for(i=0;i<COULOMBDIM;i++) coulombTerms[i] = forceField->coulombTerms[i];

#ifdef ENABLE_OMP
#pragma omp parallel for private(i,ai,aj,rijx,rijy,rijz,rij2,rij6,coulombTerm) reduction(+:energy)
#endif
	for (  i = 0; i < numberOfCoulomb; i++ )
	{
		int ai, aj;
		double rij2;
		double rijx, rijy, rijz;
		double C;
	
		ai     = (int)coulombTerms[0][i];
		aj     = (int)coulombTerms[1][i];
		C = coulombTerms[2][i]*m->atoms[ai].charge*m->atoms[aj].charge;

		rijx =  m->atoms[ai].coordinates[0] -  m->atoms[aj].coordinates[0];
		rijy =  m->atoms[ai].coordinates[1] -  m->atoms[aj].coordinates[1];
		rijz =  m->atoms[ai].coordinates[2] -  m->atoms[aj].coordinates[2];

		rij2 = (rijx * rijx + rijy * rijy + rijz * rijz);
		if ( rij2 < POSEPS*POSEPS ) rij2 = POSEPS*POSEPS;	
		rij2 = 1.0/rij2;

		//printf("ai aj c=%d %d %f\n",ai,aj,coulombTerms[2][i]);
		energy += C  * sqrt( rij2 );
	}  
	//printf("numberOfCoulomb=%d\n",numberOfCoulomb);
	//printf("EC=%f\n",energy);
	return energy;
}
/**********************************************************************/
double calculateEnergySuttonChenAmber(ForceField* forceField,Molecule* molecule)
{
	double energy = 0;
	Molecule* m = molecule;

	energy = calculateRhoSuttonChenAmber(forceField,m);
	return energy;
}
/**********************************************************************/
static double calculateEnergyHydrogenBonded612Amber(ForceField* forceField,Molecule* molecule)
{
	int i;
	Molecule* m = molecule;
	double* hydrogenBondedTerms[HYDROGENBONDEDDIM];
	int numberOfHydrogenBonded =  forceField->numberOfHydrogenBonded;
	double energy = 0.0;

	//printf("Begin energ Hydrogen-bonded\n");
	for(i=0;i<HYDROGENBONDEDDIM;i++) hydrogenBondedTerms[i] = forceField->hydrogenBondedTerms[i];

	/* Hydrogen-bonded term */
#ifdef ENABLE_OMP
#pragma omp parallel for private(i) reduction(+:energy)
#endif
	for (  i = 0; i < numberOfHydrogenBonded; i++ )
	{
		int ai, aj,ax;
		double rij2, rij4, rij8, rij12;
		double rijx, rijy, rijz;
		double Cij, Dij;
		double cosBeta = 1.0;
		ai = (int)hydrogenBondedTerms[0][i];
		aj = (int)hydrogenBondedTerms[1][i];
		ax = (int)hydrogenBondedTerms[2][i];
		Cij = hydrogenBondedTerms[3][i];
		Dij = hydrogenBondedTerms[4][i];

		cosBeta = getCosBeta(m, ai, ax, aj);

		rijx =  m->atoms[ai].coordinates[0] - m->atoms[aj].coordinates[0];
		rijy =  m->atoms[ai].coordinates[1] - m->atoms[aj].coordinates[1];
		rijz =  m->atoms[ai].coordinates[2] - m->atoms[aj].coordinates[2];

		rij2 = rijx * rijx + rijy * rijy + rijz * rijz;
		if ( rij2 < POSEPS*POSEPS )
		{
			fprintf(forceField->logfile, "i = %d j = %d\n",ai,aj);
			fflush(forceField->logfile);
			rij2 = POSEPS*POSEPS;	
		}
		rij4 = rij2 * rij2;
                rij8 = rij4 * rij4;
                rij12 = rij8 * rij4;

		energy += Cij / rij12 - cosBeta*Dij / (rij4*rij2);

		/*
		fprintf(forceField->logfile, "C =%f D = %f  r= %f e= %f\n", Cij,Dij ,sqrt(rij2),energy);
		fflush(forceField->logfile);
		*/

	}  
	//printf("End energ Hydrogen-bonded\n");

	return energy;		
}
/**********************************************************************/
static double calculateEnergyHydrogenBonded1012Amber(ForceField* forceField,Molecule* molecule)
{
	int i;
	Molecule* m = molecule;
	double* hydrogenBondedTerms[HYDROGENBONDEDDIM];
	int numberOfHydrogenBonded =  forceField->numberOfHydrogenBonded;
	double energy = 0.0;

	//printf("Begin energ Hydrogen-bonded\n");
	for(i=0;i<HYDROGENBONDEDDIM;i++) hydrogenBondedTerms[i] = forceField->hydrogenBondedTerms[i];

	/* Hydrogen-bonded term */
#ifdef ENABLE_OMP
#pragma omp parallel for private(i,ai,aj,Cij,Dij,rijx,rijy,rijz,rij2,rij4,rij6,rij10,rij12) reduction(+:energy)
#endif
	for (  i = 0; i < numberOfHydrogenBonded; i++ )
	{
		int ai, aj,ax;
		double rij2, rij4, rij8, rij12;
		double rijx, rijy, rijz;
		double Cij, Dij;
		double cosBeta = 1.0;
		ai = (int)hydrogenBondedTerms[0][i];
		aj = (int)hydrogenBondedTerms[1][i];
		ax = (int)hydrogenBondedTerms[2][i];
		Cij = hydrogenBondedTerms[3][i];
		Dij = hydrogenBondedTerms[4][i];

		cosBeta = getCosBeta(m, ai, ax, aj);

		rijx =  m->atoms[ai].coordinates[0] - m->atoms[aj].coordinates[0];
		rijy =  m->atoms[ai].coordinates[1] - m->atoms[aj].coordinates[1];
		rijz =  m->atoms[ai].coordinates[2] - m->atoms[aj].coordinates[2];

		rij2 = rijx * rijx + rijy * rijy + rijz * rijz;
		if ( rij2 < POSEPS*POSEPS )
		{
			fprintf(forceField->logfile, "i = %d j = %d\n",ai,aj);
			fflush(forceField->logfile);
			rij2 = POSEPS*POSEPS;	
		}
		rij4 = rij2 * rij2;
                rij8 = rij4 * rij4;
                rij12 = rij8 * rij4;

		energy += Cij / rij12 - cosBeta*Dij / (rij8*rij2);

		/*
		fprintf(forceField->logfile, "C =%f D = %f  r= %f e= %f\n", Cij,Dij ,sqrt(rij2),energy);
		fflush(forceField->logfile);
		*/

	}  
	//printf("End energ Hydrogen-bonded\n");

	return energy;		
}
/**********************************************************************/
static double calculateEnergyHydrogenBondedMorseAmber(ForceField* forceField,Molecule* molecule)
{
	// Ref : J. Phys. Chem. B 1997, 101, 4851-4859
	int i;
	Molecule* m = molecule;
	double* hydrogenBondedTerms[HYDROGENBONDEDDIM];
	int numberOfHydrogenBonded =  forceField->numberOfHydrogenBonded;
	double energy = 0.0;

	//printf("Begin energ Hydrogen-bonded\n");
	for(i=0;i<HYDROGENBONDEDDIM;i++) hydrogenBondedTerms[i] = forceField->hydrogenBondedTerms[i];

	/* Hydrogen-bonded term */
#ifdef ENABLE_OMP
#pragma omp parallel for private(i,ai,aj,Cij,Dij,rijx,rijy,rijz,rij2,rij4,rij6,rij10,rij12) reduction(+:energy)
#endif
	for (  i = 0; i < numberOfHydrogenBonded; i++ )
	{
		int ai, aj,ax;
		double rij, rij2;
		double rijx, rijy, rijz;
		double alpha, Re,De;
		double cosBeta = 1.0;
		double X = 0.0;
		ai = (int)hydrogenBondedTerms[0][i];
		aj = (int)hydrogenBondedTerms[1][i];
		ax = (int)hydrogenBondedTerms[2][i];
		alpha = hydrogenBondedTerms[3][i];
		Re = hydrogenBondedTerms[4][i];
		De = hydrogenBondedTerms[5][i];

		cosBeta = getCosBeta(m, ai, ax, aj);

		rijx =  m->atoms[ai].coordinates[0] - m->atoms[aj].coordinates[0];
		rijy =  m->atoms[ai].coordinates[1] - m->atoms[aj].coordinates[1];
		rijz =  m->atoms[ai].coordinates[2] - m->atoms[aj].coordinates[2];

		rij2 = rijx * rijx + rijy * rijy + rijz * rijz;
		if ( rij2 < POSEPS*POSEPS )
		{
			fprintf(forceField->logfile, "i = %d j = %d\n",ai,aj);
			fflush(forceField->logfile);
			rij2 = POSEPS*POSEPS;	
		}
		rij = sqrt(rij2);
		X = exp(-alpha*(rij-Re));
		energy += De*(X*X-2*X)*cosBeta;
		/*
		fprintf(forceField->logfile, "alpha =%f Re = %f  De= %f r = %f  cosBeta = %f energy= %f\n", alpha,Re,De, rij,cosBeta, energy);
		fflush(forceField->logfile);
		*/

	}  
	//printf("End energ Hydrogen-bonded\n");
	return energy;		
}
/**********************************************************************/
static double calculateEnergyHydrogenBondedAmber(ForceField* forceField,Molecule* molecule)
{
	if(forceField->options.hydrogenBonded612) return calculateEnergyHydrogenBonded612Amber(forceField,molecule);
	else if(forceField->options.hydrogenBonded1012) return calculateEnergyHydrogenBonded1012Amber(forceField,molecule);
	else if(forceField->options.hydrogenBondedMorse) return calculateEnergyHydrogenBondedMorseAmber(forceField,molecule);
	return 0;
}
/**********************************************************************/
double calculateEnergyOutOfPlaneAmber(ForceField* forceField,Molecule* molecule)
{
	int i;
	Molecule* m = molecule;
	double energy = 0.0;
	double* outOfPlaneTerms[OUTOFPLANEDIM];
	int numberOfOutOfPlaneTerms = forceField->numberOfOutOfPlaneTerms;

	for(i=0;i<OUTOFPLANEDIM;i++) outOfPlaneTerms[i] = forceField->outOfPlaneTerms[i];

#ifdef ENABLE_OMP
#pragma omp parallel for private(i) 
#endif
	for (  i = 0; i < numberOfOutOfPlaneTerms; i++ )
	{
		int ai, aj, ak, al;

		double rijx, rijy, rijz;
		double rkjx, rkjy, rkjz;
		double rljx, rljy, rljz;

		double angle, e2, c2;
		double dot;
		double rilx, rily, rilz;
		double rklx, rkly, rklz;
		double rij2, rkj2, ril2, rkl2, rlj2;
		double bkk2, cosine;
		double dt, dt2, dt3, dt4;

		boolean wdc = (outOfPlaneTerms[4][i]==1);
		double force = outOfPlaneTerms[5][i];
		double h3 = outOfPlaneTerms[6][i];
		double h4 = outOfPlaneTerms[7][i];
		double h5 = outOfPlaneTerms[8][i];
		double h6 = outOfPlaneTerms[9][i];

		ai = (int)outOfPlaneTerms[0][i];
		aj = (int)outOfPlaneTerms[1][i];
		ak = (int)outOfPlaneTerms[2][i];
		al = (int)outOfPlaneTerms[3][i];

		rijx = m->atoms[ai].coordinates[0] - m->atoms[aj].coordinates[0];
		rijy = m->atoms[ai].coordinates[1] - m->atoms[aj].coordinates[1];
		rijz = m->atoms[ai].coordinates[2] - m->atoms[aj].coordinates[2];

		rkjx = m->atoms[ak].coordinates[0] - m->atoms[aj].coordinates[0];
		rkjy = m->atoms[ak].coordinates[1] - m->atoms[aj].coordinates[1];
		rkjz = m->atoms[ak].coordinates[2] - m->atoms[aj].coordinates[2];

		rljx = m->atoms[al].coordinates[0] - m->atoms[aj].coordinates[0];
		rljy = m->atoms[al].coordinates[1] - m->atoms[aj].coordinates[1];
		rljz = m->atoms[al].coordinates[2] - m->atoms[aj].coordinates[2];

		rilx = m->atoms[ai].coordinates[0] - m->atoms[al].coordinates[0];
		rily = m->atoms[ai].coordinates[1] - m->atoms[al].coordinates[1];
		rilz = m->atoms[ai].coordinates[2] - m->atoms[al].coordinates[2];

		rklx = m->atoms[ak].coordinates[0] - m->atoms[al].coordinates[0];
		rkly = m->atoms[ak].coordinates[1] - m->atoms[al].coordinates[1];
		rklz = m->atoms[ak].coordinates[2] - m->atoms[al].coordinates[2];

		if(wdc) //  'W-D-C'
		{
		     	rij2 = rijx*rijx + rijy*rijy + rijz*rijz;
               		rkj2 = rkjx*rkjx + rkjy*rkjy + rkjz*rkjz;
               		dot = rijx*rkjx+rijy*rkjy+rijz*rkjz;
               		c2 = rij2*rkj2 - dot*dot;
		}
		else // ALLINGER
		{
			ril2 = rilx*rilx + rily*rily + rilz*rilz;
               		rkl2 = rklx*rklx + rkly*rkly + rklz*rklz;
               		dot = rilx*rklx + rily*rkly + rilz*rklz;
               		c2 = ril2*rkl2 - dot*dot;
		}
		// energy
		e2 = rljx*(rijy*rkjz-rijz*rkjy) + rljy*(rijz*rkjx-rijx*rkjz)+ rljz*(rijx*rkjy-rijy*rkjx);
            	rlj2 = rljx*rljx + rljy*rljy + rljz*rljz;
            	if (rlj2==0.0 || c2 == 0.0) continue;
               	bkk2 = rlj2 - e2*e2/c2;
               	cosine = sqrt(bkk2/rlj2);
		if(cosine>1) cosine = 1;
		if(cosine<-1) cosine = -1;
               	angle = acos(cosine);
               	dt = angle;
               	dt2 = dt * dt;
               	dt3 = dt2 * dt;
               	dt4 = dt2 * dt2;
               	energy += force * dt2*(1.0+h3*dt+h4*dt2+h5*dt3+h6*dt4);
	}
	return energy;
}
/**********************************************************************/
static double calculateEnergyPairWise(ForceField* forceField,Molecule* molecule)
{
	int i;
	int ai, aj;
	double rij2, rij4, rij6, rij8, rij10;
	double coulombTerm;
	double rijx, rijy, rijz;
	double rij;
	double permittivityScale = 1, permittivity = 1;
	double coulombFactor;
	Molecule* m = molecule;
	double* pairWiseTerms[PAIRWISEDIM];
	int numberOfPairWise = forceField->numberOfPairWise;
	boolean useCoulomb = forceField->options.coulomb;
	boolean useVanderWals = forceField->options.vdw612;
	double energy = 0.0;
	double A, Beta;
	double  B4, B6, B8, B10;
	double c4, c6, c8, c10, b;

	for(i=0;i<PAIRWISEDIM;i++)
		pairWiseTerms[i] = forceField->pairWiseTerms[i];

	/* now for non-bonded term */
	coulombFactor = 332.05382/ ( permittivity * permittivityScale );
	/* fprintf(forceField->logfile, "number of Non Bonded terms = %d\n",numberOfPairWise);*/
	// fflush(forceField->logfile);
#ifdef ENABLE_OMP
#pragma omp parallel for private(i,A,Beta,c4,c6,c8,c10,b,rijx,rijy,rijz,rij2,rij,rij6,rij8,rij10,coulombTerm,B6,B8,B10) reduction(+:energy)
#endif
	for (  i = 0; i < numberOfPairWise; i++ )
	{
		ai     = (int)pairWiseTerms[0][i];
		aj     = (int)pairWiseTerms[1][i];
		A      = pairWiseTerms[2][i];
		Beta   = pairWiseTerms[3][i];
		c4     = pairWiseTerms[4][i];
		c6     = pairWiseTerms[5][i];
		c8     = pairWiseTerms[6][i];
		c10    = pairWiseTerms[7][i];
		b      = pairWiseTerms[8][i];

		rijx =  m->atoms[ai].coordinates[0] - m->atoms[aj].coordinates[0];
		rijy =  m->atoms[ai].coordinates[1] - m->atoms[aj].coordinates[1];
		rijz =  m->atoms[ai].coordinates[2] - m->atoms[aj].coordinates[2];

		rij2 = rijx * rijx + rijy * rijy + rijz * rijz;

		rij = sqrt( rij2 );
		rij4 = rij2 * rij2;
		rij6 = rij2 * rij2 * rij2;
		rij8 = rij6* rij2;
		rij10 = rij8 * rij2;

		if(useCoulomb)
			coulombTerm = ( m->atoms[ai].charge*m->atoms[aj].charge * coulombFactor ) / rij;
		else
			coulombTerm = 0.0;

		B4  = 0;
		B6  = 0;
		B8  = 0;
		B10 = 0;
		/* printf("A = %f Beta = %f qi = %f qj = %f rij = %f\n",A,Beta,chargei,chargej,rij);*/
		if(useVanderWals)
		{
			double fact = 1.0;
			double s = 1.0;
			double br = b*rij;
			double brk = 1.0;
			int k;

			if(fabs(c6)>1e-12)
			{
				for(k=1;k<=2*3;k++)
				{
					fact *= k;
					brk *= br;
					s += brk/fact;
				}
				B6 = c6*(1-exp(-br)*s);
			}
			if(fabs(c4)>1e-12)
			{
				fact = 1.0;
				s = 1.0;
				br = b*rij;
				brk = 1.0;
				for(k=1;k<=2*2;k++)
				{
					fact *= k;
					brk *= br;
					s += brk/fact;
				}
				B4 = c4*(1-exp(-br)*s);
			}

			if(fabs(c8)>1e-12)
			{
				fact = 1.0;
				s = 1.0;
				br = b*rij;
				brk = 1.0;
				for(k=1;k<=2*4;k++)
				{
					fact *= k;
					brk *= br;
					s += brk/fact;
				}
				B8 = c8*(1-exp(-br)*s);
			}

			if(fabs(c10)>1e-12)
			{
				fact = 1.0;
				s = 1.0;
				br = b*rij;
				brk = 1.0;
				for(k=1;k<=2*5;k++)
				{
					fact *= k;
					brk *= br;
					s += brk/fact;
				}
				B10 = c10*(1-exp(-br)*s);
			}
		}
					


		energy += A*exp(-Beta*rij)
			- B4 / rij4 
			- B6 / rij6 
			- B8 / rij8 
			- B10 / rij10 
			+ coulombTerm;
	}  
	return energy;
}
#endif /* ENABLE_CL*/

/**********************************************************************/
static void calculateEnergyAmber(ForceField* forceField)
{
	Molecule* m = &forceField->molecule;

	forceField->molecule.potentialEnergy =  calculateEnergyTmpAmber(forceField,m);
}
/**********************************************************************/
static double calculateEnergyTmpAmber(ForceField* forceField,Molecule* molecule)
{
	double energy = 0.0;
#ifdef DEBUG
        TimerType timer;
#endif

#ifdef ENABLE_CL
	CLProp clProp = getCLProp();
	size_t local = 64;
	size_t global;
	int i;

	fprintf(forceField->logfile, "Begin calculateEnergyTmpAmber\n");
	fflush(forceField->logfile);
	global = forceField->nMaxTerms;
	clEnqueueNDRangeKernel(clProp.command_queue, forceField->initEnergyKernel, 1, NULL, &global, NULL, 0, NULL, NULL);
	clFinish(clProp.command_queue);

#ifdef DEBUG
        timer_init(timer);
       	timer_start( timer );
#endif
	if(forceField->numberOfStretchTerms>0)
	{
	global = forceField->numberOfStretchTerms;
	clEnqueueNDRangeKernel(clProp.command_queue, forceField->addEnergyBondAmberKernel, 1, NULL, &global, NULL, 0, NULL, NULL);
	clFinish(clProp.command_queue);
	}

	if(forceField->numberOfBendTerms>0)
	{
	global = forceField->numberOfBendTerms;
	clEnqueueNDRangeKernel(clProp.command_queue, forceField->addEnergyBendAmberKernel, 1, NULL, &global, NULL, 0, NULL, NULL);
	clFinish(clProp.command_queue);
	}

	if(forceField->numberOfDihedralTerms>0)
	{
	global = forceField->numberOfDihedralTerms;
	clEnqueueNDRangeKernel(clProp.command_queue, forceField->addEnergyDihedralAmberKernel, 1, NULL, &global, NULL, 0, NULL, NULL);
	clFinish(clProp.command_queue);
	}

	if(forceField->numberOfVdw612>0)
	{
	global = forceField->numberOfVdw612;
	clEnqueueNDRangeKernel(clProp.command_queue, forceField->addEnergyVdw612AmberKernel, 1, NULL, &global, NULL, 0, NULL, NULL);
	clFinish(clProp.command_queue);
	}

	if(forceField->numberOfSuttonChen>0)
	{
	global = forceField->numberOfSuttonChen;
	clEnqueueNDRangeKernel(clProp.command_queue, forceField->initRhoKernel, 1, NULL, &global, NULL, 0, NULL, NULL);
	clEnqueueNDRangeKernel(clProp.command_queue, forceField->computeRhoKernel, 1, NULL, &global, NULL, 0, NULL, NULL);
	global = forceField->molecule.nAtoms;
	clEnqueueNDRangeKernel(clProp.command_queue, forceField->reduceRhoKernel, 1, NULL, &global, NULL, 0, NULL, NULL);
	global = forceField->numberOfSuttonChen;
	clEnqueueNDRangeKernel(clProp.command_queue, forceField->addEnergySuttoChenAmberKernel, 1, NULL, &global, NULL, 0, NULL, NULL);
	clFinish(clProp.command_queue);
	}

	if(forceField->numberOfHydrogenBonded>0)
	{
	global = forceField->numberOfHydrogenBonded;
	clEnqueueNDRangeKernel(clProp.command_queue, forceField->addEnergyHydrogenBondedAmberKernel, 1, NULL, &global, NULL, 0, NULL, NULL);
	clFinish(clProp.command_queue);
	}

	if(forceField->numberOfPairWise>0)
	{
	global = forceField->numberOfPairWise;
	clEnqueueNDRangeKernel(clProp.command_queue, forceField->addEnergyPairWiseKernel, 1, NULL, &global, NULL, 0, NULL, NULL);
	clFinish(clProp.command_queue);
	}


	clEnqueueReadBuffer(clProp.command_queue, forceField->energyBufferCL, CL_TRUE, 0, sizeof(cl_float)*forceField->nMaxTerms, forceField->energyBufferCPU, 0, NULL, NULL);
	energy = 0;
	for(i=0;i<forceField->nMaxTerms;i++) energy+=forceField->energyBufferCPU[i];
#ifdef DEBUG
       	timer_stop(timer);
        fprintf(forceField->logfile, "time (s) = %f\n", timer_get(timer)*1e-6);
	fflush(forceField->logfile);
#endif
#else // ENABLE_CL
#ifdef DEBUG
        timer_init(timer);
       	timer_start( timer );
#endif
	if(!strstr(forceField->options.chargesType,"BEGIN")) setCharges(forceField);
	energy +=calculateEnergyBondAmber(forceField,molecule);
	energy +=calculateEnergyBendAmber(forceField,molecule);
	energy +=calculateEnergyStrBendAmber(forceField,molecule);
	energy +=calculateEnergyDihedralAmber(forceField,molecule);
	energy +=calculateEnergyImproperTorsionAmber(forceField,molecule);
	energy +=calculateEnergyOutOfPlaneAmber(forceField,molecule);
	energy +=calculateEnergyVdw612Amber(forceField,molecule);
	energy +=calculateEnergySuttonChenAmber(forceField,molecule);
	energy +=calculateEnergyCoulombAmber(forceField,molecule);
	energy +=calculateEnergyHydrogenBondedAmber(forceField,molecule);
	energy +=calculateEnergyPairWise(forceField,molecule);
	energy +=getH4Energy(forceField, molecule);
#ifdef DEBUG
       	timer_stop(timer);
        fprintf(forceField->logfile, "time (s) = %f\n", timer_get(timer)*1e-6);
	fflush(forceField->logfile);
#endif
#endif // ENABLE_CL

	return energy;		
}
/**********************************************************************/
#ifdef ENABLE_CL
void initCLForceField (ForceField* forceField)
{
	CLProp clProp = getCLProp();
	cl_int err;
	int i;
	//cl_int2* clint2 = NULL;
	cl_int4* clint4 = NULL;
	cl_int8* clint8 = NULL;
	cl_float* clfloat = NULL;
	cl_float2* clfloat2 = NULL;
	cl_float4* clfloat4 = NULL;
	cl_float8* clfloat8 = NULL;
	int numberOfStretchTerms = forceField->numberOfStretchTerms;
	double* bondStretchTerms[STRETCHDIM];
	int numberOfBendTerms = forceField->numberOfBendTerms;
	double* angleBendTerms[BENDDIM];
	int numberOfDihedralTerms = forceField->numberOfDihedralTerms;
	double* dihedralAngleTerms[DIHEDRALDIM];
	double* vdw612Terms[VDW612DIM];
	int numberOfVdw612 = forceField->numberOfVdw612;
	double* hydrogenBondedTerms[HYDROGENBONDEDDIM];
	int numberOfHydrogenBonded =  forceField->numberOfHydrogenBonded;
	double* pairWiseTerms[PAIRWISEDIM];
	int numberOfPairWise = forceField->numberOfPairWise;
	size_t local = 64;
	size_t global;
	int nMaxTerms = 0;
	int useCoulomb = forceField->options.coulomb;
	int* gradCounters = malloc(forceField->molecule.nAtoms*sizeof(int));
	int* suttonChenCounters = malloc(forceField->molecule.nAtoms*sizeof(int));
	int ea,eb,ec;
	int maxCounters=0;
	int maxSuttonChenCounters=0;
	int maxRattleCounters=0;
	int numberOfRattleConstraintsTerms = 0;
	double* rattleConstraintsTerms[RATTLEDIM];

	numberOfRattleConstraintsTerms = forceField->numberOfRattleConstraintsTerms;
	for( i=0; i<RATTLEDIM;i++) rattleConstraintsTerms[i] = forceField->rattleConstraintsTerms[i];

        for( i=0; i<forceField->molecule.nAtoms;i++) gradCounters[i] = 0;

        for( i=0; i<STRETCHDIM;i++) bondStretchTerms[i] = forceField->bondStretchTerms[i];
	for( i=0; i<BENDDIM;i++) angleBendTerms[i] = forceField->angleBendTerms[i]; 
	for(i=0;i<DIHEDRALDIM;i++) dihedralAngleTerms[i] = forceField->dihedralAngleTerms[i];
	for(i=0;i<VDW612DIM;i++) vdw612Terms[i] = forceField->vdw612Terms[i];
	for(i=0;i<HYDROGENBONDEDDIM;i++) hydrogenBondedTerms[i] = forceField->hydrogenBondedTerms[i];
	for(i=0;i<PAIRWISEDIM;i++) pairWiseTerms[i] = forceField->pairWiseTerms[i];
	if(numberOfStretchTerms>nMaxTerms) nMaxTerms = numberOfStretchTerms;
	if(numberOfBendTerms>nMaxTerms) nMaxTerms = numberOfBendTerms ;
	if(numberOfDihedralTerms>nMaxTerms) nMaxTerms = numberOfDihedralTerms;
	if(numberOfVdw612>nMaxTerms) nMaxTerms = numberOfVdw612;
	if(numberOfHydrogenBonded>nMaxTerms) nMaxTerms =numberOfHydrogenBonded ;
	if(numberOfPairWise>nMaxTerms) nMaxTerms =numberOfPairWise ;
	forceField->nMaxTerms = nMaxTerms;

	// create a program from the kernel source code
	forceField->programMM = clCreateProgramWithSource(clProp.context,1,(const char **) &mmCLSource, NULL, &err);
#ifdef DEBUG
	fprintf(forceField->logfile, "err = %d\n",err);
	fflush(forceField->logfile);
#endif
	// compile the program
	if (clBuildProgram(forceField->programMM, 0, NULL, NULL, NULL, NULL) != CL_SUCCESS)
	{
		char build[2048];
		fprintf(forceField->logfile, "Error building MM CL program\n");
		fflush(forceField->logfile);
		clGetProgramBuildInfo(forceField->programMM, clProp.device_id, CL_PROGRAM_BUILD_LOG, 2048, build, NULL);
		fprintf(forceField->logfile, "Build Log:\n%s\n",build);
		fflush(forceField->logfile);
		exit(1);
	}

	forceField->atomsCPU = malloc(forceField->molecule.nAtoms*sizeof(cl_float8));
#ifdef DEBUG
	fprintf(forceField->logfile, "End malloc\n");
	fflush(forceField->logfile);
#endif
	for(i=0;i<forceField->molecule.nAtoms;i++)
	{
		forceField->atomsCPU[i].s[0] = forceField->molecule.atoms[i].coordinates[0];
		forceField->atomsCPU[i].s[1] = forceField->molecule.atoms[i].coordinates[1];
		forceField->atomsCPU[i].s[2] = forceField->molecule.atoms[i].coordinates[2];
		forceField->atomsCPU[i].s[3] = forceField->molecule.atoms[i].charge;
		forceField->atomsCPU[i].s[4] = forceField->molecule.atoms[i].mass;
		// velocity
		forceField->atomsCPU[i].s[5] = 0;
		forceField->atomsCPU[i].s[6] = 0;
		forceField->atomsCPU[i].s[7] = 0;
	}
	printf("create buffers for the input and ouput\n");
	forceField->atomsCL = clCreateBuffer(clProp.context, CL_MEM_READ_WRITE, sizeof(cl_float8) * forceField->molecule.nAtoms, NULL, &ea);
	if(ea!=CL_SUCCESS)
	{
		fprintf(forceField->logfile, "I cannot create atoms Buffer\n");
		fflush(forceField->logfile);
		exit(1);
	}
	clEnqueueWriteBuffer(clProp.command_queue, forceField->atomsCL, CL_TRUE, 0, sizeof(cl_float8) * forceField->molecule.nAtoms, forceField->atomsCPU, 0, NULL, NULL);
#ifdef DEBUG
	fprintf(forceField->logfile, "End atomsCL\n");
	fflush(forceField->logfile);

	fprintf(forceField->logfile, "nMaxTerms=%d\n",nMaxTerms);
	fflush(forceField->logfile);
#endif

	forceField->energyBufferCPU = malloc(nMaxTerms*sizeof(cl_float));
	forceField->energyBufferCL = clCreateBuffer(clProp.context, CL_MEM_READ_WRITE, sizeof(cl_float)*nMaxTerms, NULL, &ea);
	if(ea!=CL_SUCCESS)
	{
		fprintf(forceField->logfile, "I cannot create energies Buffer\n");
		fflush(forceField->logfile);
		exit(1);
	}

	if(numberOfStretchTerms>0)
	{
	clint4 = malloc(numberOfStretchTerms*sizeof(cl_int4));
	clfloat2 = malloc(numberOfStretchTerms*sizeof(cl_float2));
        for( i=0; i<forceField->molecule.nAtoms;i++) gradCounters[i] = 0;
	for(i=0;i<numberOfStretchTerms;i++)
	{
		clint4[i].s[0] = (int) bondStretchTerms[0][i];
		clint4[i].s[1] = (int) bondStretchTerms[1][i];
		clint4[i].s[2] = gradCounters[clint4[i].s[0]]++;
		clint4[i].s[3] = gradCounters[clint4[i].s[1]]++;
		
		clfloat2[i].s[0] = bondStretchTerms[2][i];
		clfloat2[i].s[1] = bondStretchTerms[3][i];
	}
	forceField->bondIndexCL = clCreateBuffer(clProp.context, CL_MEM_READ_ONLY, sizeof(cl_int4)*numberOfStretchTerms, NULL, &ea);
	forceField->bondTermsCL = clCreateBuffer(clProp.context, CL_MEM_READ_ONLY, sizeof(cl_float2)*numberOfStretchTerms, NULL, &eb);
	if(ea!=CL_SUCCESS||eb!=CL_SUCCESS)
	{
		fprintf(forceField->logfile, "I cannot bonds Buffer\n");
		fflush(forceField->logfile);
		exit(1);
	}
	ea = clEnqueueWriteBuffer(clProp.command_queue, forceField->bondIndexCL, CL_TRUE, 0, sizeof(cl_int4) * numberOfStretchTerms, clint4, 0, NULL, NULL);
	eb = clEnqueueWriteBuffer(clProp.command_queue, forceField->bondTermsCL, CL_TRUE, 0, sizeof(cl_float2) * numberOfStretchTerms, clfloat2, 0, NULL, NULL);
	if(ea!=CL_SUCCESS||eb!=CL_SUCCESS)
	{
		fprintf(forceField->logfile, "I cannot write bonds data in GPU Buffer\n");
		fflush(forceField->logfile);
		exit(1);
	}
	free(clint4);
	free(clfloat2);
#ifdef DEBUG
	fprintf(forceField->logfile, "End numberOfStretchTerms\n");
	fflush(forceField->logfile);
#endif
        for( i=0; i<forceField->molecule.nAtoms;i++) if(maxCounters<gradCounters[i]) maxCounters = gradCounters[i];
	}

	if(numberOfBendTerms>0)
	{
	clint8 = malloc(numberOfBendTerms*sizeof(cl_int8));
	clfloat2 = malloc(numberOfBendTerms*sizeof(cl_float2));
        for( i=0; i<forceField->molecule.nAtoms;i++) gradCounters[i] = 0;
	for(i=0;i<numberOfBendTerms;i++)
	{
		clint8[i].s[0] = (int) angleBendTerms[0][i];
		clint8[i].s[1] = (int) angleBendTerms[1][i];
		clint8[i].s[2] = (int) angleBendTerms[2][i];
		clint8[i].s[3] =  gradCounters[clint8[i].s[0]]++;
		clint8[i].s[4] =  gradCounters[clint8[i].s[1]]++;
		clint8[i].s[5] =  gradCounters[clint8[i].s[2]]++;
		clint8[i].s[6] =  -1;
		clint8[i].s[7] =  -1;

		clfloat2[i].s[0] = angleBendTerms[3][i];
		clfloat2[i].s[1] = angleBendTerms[4][i];
	}
	forceField->bendIndexCL = clCreateBuffer(clProp.context, CL_MEM_READ_ONLY, sizeof(cl_int8)*numberOfBendTerms, NULL, NULL);
	forceField->bendTermsCL = clCreateBuffer(clProp.context, CL_MEM_READ_ONLY, sizeof(cl_float2)*numberOfBendTerms, NULL, NULL);
	clEnqueueWriteBuffer(clProp.command_queue, forceField->bendIndexCL, CL_TRUE, 0, sizeof(cl_int8) * numberOfBendTerms, clint8, 0, NULL, NULL);
	clEnqueueWriteBuffer(clProp.command_queue, forceField->bendTermsCL, CL_TRUE, 0, sizeof(cl_float2) * numberOfBendTerms, clfloat2, 0, NULL, NULL);
	free(clint8);
	free(clfloat2);
        for( i=0; i<forceField->molecule.nAtoms;i++) if(maxCounters<gradCounters[i]) maxCounters = gradCounters[i];
#ifdef DEBUG
	fprintf(forceField->logfile, "End Bend\n");
	fflush(forceField->logfile);
#endif
	}

	if(numberOfDihedralTerms>0)
	{
	clint8 = malloc(numberOfDihedralTerms*sizeof(cl_int8));
	clfloat4 = malloc(numberOfDihedralTerms*sizeof(cl_float4));
        for( i=0; i<forceField->molecule.nAtoms;i++) gradCounters[i] = 0;
	for(i=0;i<numberOfDihedralTerms;i++)
	{
		clint8[i].s[0] = (int) dihedralAngleTerms[0][i];
		clint8[i].s[1] = (int) dihedralAngleTerms[1][i];
		clint8[i].s[2] = (int) dihedralAngleTerms[2][i];
		clint8[i].s[3] = (int) dihedralAngleTerms[3][i];
		clint8[i].s[4] =  gradCounters[clint8[i].s[0]]++;
		clint8[i].s[5] =  gradCounters[clint8[i].s[1]]++;
		clint8[i].s[6] =  gradCounters[clint8[i].s[2]]++;
		clint8[i].s[7] =  gradCounters[clint8[i].s[3]]++;

		clfloat4[i].s[0] = dihedralAngleTerms[4][i];
		clfloat4[i].s[1] = dihedralAngleTerms[5][i];
		clfloat4[i].s[2] = dihedralAngleTerms[6][i];
		clfloat4[i].s[3] = dihedralAngleTerms[7][i];
	}
	forceField->dihedralIndexCL = clCreateBuffer(clProp.context, CL_MEM_READ_ONLY, sizeof(cl_int8)*numberOfDihedralTerms, NULL, NULL);
	forceField->dihedralTermsCL = clCreateBuffer(clProp.context, CL_MEM_READ_ONLY, sizeof(cl_float4)*numberOfDihedralTerms, NULL, NULL);
	clEnqueueWriteBuffer(clProp.command_queue, forceField->dihedralIndexCL, CL_TRUE, 0, sizeof(cl_int8) * numberOfDihedralTerms, clint8, 0, NULL, NULL);
	clEnqueueWriteBuffer(clProp.command_queue, forceField->dihedralTermsCL, CL_TRUE, 0, sizeof(cl_float4) * numberOfDihedralTerms, clfloat4, 0, NULL, NULL);
	free(clint8);
	free(clfloat4);
        for( i=0; i<forceField->molecule.nAtoms;i++) if(maxCounters<gradCounters[i]) maxCounters = gradCounters[i];
#ifdef DEBUG
	fprintf(forceField->logfile, "End dihed\n");
	fflush(forceField->logfile);
#endif
	}


	forceField->improperTorsionIndexCL = 0;


	if(numberOfVdw612>0)
	{
	clint4 = malloc(numberOfVdw612*sizeof(cl_int4));
	clfloat4 = malloc(numberOfVdw612*sizeof(cl_float4));
        for( i=0; i<forceField->molecule.nAtoms;i++) gradCounters[i] = 0;
	for(i=0;i<numberOfVdw612;i++)
	{
		clint4[i].s[0] = (int) vdw612Terms[0][i];
		clint4[i].s[1] = (int) vdw612Terms[1][i];
		clint4[i].s[2] =  gradCounters[clint4[i].s[0]]++;
		clint4[i].s[3] =  gradCounters[clint4[i].s[1]]++;

		clfloat4[i].s[0] = vdw612Terms[2][i];
		clfloat4[i].s[1] = vdw612Terms[3][i];
		clfloat4[i].s[2] = vdw612Terms[4][i];
		clfloat4[i].s[3] = 0;
	}
	forceField->vdw612IndexCL = clCreateBuffer(clProp.context, CL_MEM_READ_ONLY, sizeof(cl_int4)*numberOfVdw612, NULL, &ea);
	forceField->vdw612TermsCL = clCreateBuffer(clProp.context, CL_MEM_READ_ONLY, sizeof(cl_float4)*numberOfVdw612, NULL, &eb);
	if(ea!=CL_SUCCESS ||eb!=CL_SUCCESS)
	{
		fprintf(forceField->logfile, "I cannot create Buffers for non-bonded terms\n");
		fflush(forceField->logfile);
		exit(1);
	}
	clEnqueueWriteBuffer(clProp.command_queue, forceField->vdw612IndexCL, CL_TRUE, 0, sizeof(cl_int4) * numberOfVdw612, clint4, 0, NULL, NULL);
	clEnqueueWriteBuffer(clProp.command_queue, forceField->vdw612TermsCL, CL_TRUE, 0, sizeof(cl_float4) * numberOfVdw612, clfloat4, 0, NULL, NULL);
	free(clint4);
	free(clfloat4);
        for( i=0; i<forceField->molecule.nAtoms;i++) if(maxCounters<gradCounters[i]) maxCounters = gradCounters[i];
#ifdef DEBUG
	fprintf(forceField->logfile, "End NonBond\n");
	fflush(forceField->logfile);
#endif
	}
	if(numberOfSuttonChen>0)
	{
	clint4 = malloc(numberOfSuttonChen*sizeof(cl_int4));
	clfloat8 = malloc(numberOfSuttonChen*sizeof(cl_float8));
        for( i=0; i<forceField->molecule.nAtoms;i++) suttonChenCounters[i] = 0;
	for(i=0;i<numberOfSuttonChen;i++)
	{
		clint4[i].s[0] = (int) suttonChenTerms[0][i];
		clint4[i].s[1] = (int) suttonChenTerms[1][i];
		clint4[i].s[2] =  suttonChenCounters[clint4[i].s[0]]++;
		clint4[i].s[3] =  suttonChenCounters[clint4[i].s[1]]++;

		clfloat8[i].s[0] = suttonChenTerms[2][i];
		clfloat8[i].s[1] = suttonChenTerms[3][i];
		clfloat8[i].s[2] = suttonChenTerms[4][i];
		clfloat8[i].s[3] = suttonChenTerms[5][i];
		clfloat8[i].s[4] = suttonChenTerms[6][i];
		clfloat8[i].s[5] = 0;
		clfloat8[i].s[6] = 0;
		clfloat8[i].s[7] = 0;
	}
	forceField->suttonChenIndexCL = clCreateBuffer(clProp.context, CL_MEM_READ_ONLY, sizeof(cl_int4)*numberOfSuttonChen, NULL, &ea);
	forceField->suttonChenTermsCL = clCreateBuffer(clProp.context, CL_MEM_READ_ONLY, sizeof(cl_float8)*numberOfSuttonChen, NULL, &eb);
	if(ea!=CL_SUCCESS ||eb!=CL_SUCCESS)
	{
		fprintf(forceField->logfile, "I cannot create Buffers for non-bonded terms\n");
		fflush(forceField->logfile);
		exit(1);
	}
	clEnqueueWriteBuffer(clProp.command_queue, forceField->suttonChenIndexCL, CL_TRUE, 0, sizeof(cl_int4) * numberOfSuttonChen, clint4, 0, NULL, NULL);
	clEnqueueWriteBuffer(clProp.command_queue, forceField->suttonChenTermsCL, CL_TRUE, 0, sizeof(cl_float8) * numberOfSuttonChen, clfloat8, 0, NULL, NULL);
	free(clint4);
	free(clfloat8);
	
        for( i=0; i<forceField->molecule.nAtoms;i++) if(maxSuttonChenCounters<suttonChenCounters[i]) maxSuttonChenCounters = suttonChenCounters[i];
#ifdef DEBUG
	fprintf(forceField->logfile, "End NonBond\n");
	fflush(forceField->logfile);
#endif
	free(suttonChenCounters);
	forceField->nBlockRhoBuffer = maxSuttonChenCounters;
	forceField->rhoBufferCL = clCreateBuffer(clProp.context,  CL_MEM_READ_WRITE, sizeof(cl_float)*forceField->molecule.nAtoms*forceField->nBlockRhoBuffer, NULL, &ea);
	}

	if(numberOfHydrogenBonded>0)
	{
	clint4 = malloc(numberOfHydrogenBonded*sizeof(cl_int4));
	clfloat2 = malloc(numberOfHydrogenBonded*sizeof(cl_float2));
        for( i=0; i<forceField->molecule.nAtoms;i++) gradCounters[i] = 0;
	for(i=0;i<numberOfHydrogenBonded;i++)
	{
		clint4[i].s[0] = (int) hydrogenBondedTerms[0][i];
		clint4[i].s[1] = (int) hydrogenBondedTerms[1][i];
		clint4[i].s[2] =  gradCounters[clint4[i].s[0]]++;
		clint4[i].s[3] =  gradCounters[clint4[i].s[1]]++;

		clfloat2[i].s[0] = hydrogenBondedTerms[2][i];
		clfloat2[i].s[1] = hydrogenBondedTerms[3][i];
	}
	forceField->hydrogenBondedIndexCL = clCreateBuffer(clProp.context, CL_MEM_READ_ONLY, sizeof(cl_int4)*numberOfHydrogenBonded, NULL, NULL);
	forceField->hydrogenBondedTermsCL = clCreateBuffer(clProp.context, CL_MEM_READ_ONLY, sizeof(cl_float2)*numberOfHydrogenBonded, NULL, NULL);
	clEnqueueWriteBuffer(clProp.command_queue, forceField->hydrogenBondedIndexCL, CL_TRUE, 0, sizeof(cl_int4) * numberOfHydrogenBonded, clint4, 0, NULL, NULL);
	clEnqueueWriteBuffer(clProp.command_queue, forceField->hydrogenBondedTermsCL, CL_TRUE, 0, sizeof(cl_float2) * numberOfHydrogenBonded, clfloat2, 0, NULL, NULL);
	free(clint4);
	free(clfloat2);
        for( i=0; i<forceField->molecule.nAtoms;i++) if(maxCounters<gradCounters[i]) maxCounters = gradCounters[i];
#ifdef DEBUG
	fprintf(forceField->logfile, "End numberOfHydrogenBonded\n");
	fflush(forceField->logfile);
#endif
	}

	if(numberOfPairWise>0)
	{
	clint4 = malloc(numberOfPairWise*sizeof(cl_int4));
	clfloat8 = malloc(numberOfPairWise*sizeof(cl_float8));
        for( i=0; i<forceField->molecule.nAtoms;i++) gradCounters[i] = 0;
	for(i=0;i<numberOfPairWise;i++)
	{
		clint4[i].s[0] = (int) pairWiseTerms[0][i];
		clint4[i].s[1] = (int) pairWiseTerms[1][i];
		clint4[i].s[2] =  gradCounters[clint4[i].s[0]]++;
		clint4[i].s[3] =  gradCounters[clint4[i].s[1]]++;

		clfloat8[i].s[0] = pairWiseTerms[2][i];
		clfloat8[i].s[1] = pairWiseTerms[3][i];
		clfloat8[i].s[2] = pairWiseTerms[4][i];
		clfloat8[i].s[3] = pairWiseTerms[5][i];
		clfloat8[i].s[4] = pairWiseTerms[6][i];
		clfloat8[i].s[5] = pairWiseTerms[7][i];
		clfloat8[i].s[6] = 0;
		clfloat8[i].s[7] = 0;
	}
	forceField->pairWiseIndexCL = clCreateBuffer(clProp.context, CL_MEM_READ_ONLY, sizeof(cl_int4)*numberOfPairWise, NULL, NULL);
	forceField->pairWiseTermsCL = clCreateBuffer(clProp.context, CL_MEM_READ_ONLY, sizeof(cl_float8)*numberOfPairWise, NULL, NULL);
	clEnqueueWriteBuffer(clProp.command_queue, forceField->pairWiseIndexCL, CL_TRUE, 0, sizeof(cl_int4) * numberOfPairWise, clint4, 0, NULL, NULL);
	clEnqueueWriteBuffer(clProp.command_queue, forceField->pairWiseTermsCL, CL_TRUE, 0, sizeof(cl_float8) * numberOfPairWise, clfloat8, 0, NULL, NULL);
	free(clint4);
	free(clfloat8);
        for( i=0; i<forceField->molecule.nAtoms;i++) if(maxCounters<gradCounters[i]) maxCounters = gradCounters[i];
#ifdef DEBUG
	fprintf(forceField->logfile, "End numberOfPairWise\n");
	fflush(forceField->logfile);
#endif
	}
	if(numberOfRattleConstraintsTerms>0)
	{
	cl_int* rattleCounters = malloc(forceField->molecule.nAtoms*sizeof(cl_int));
	clint4 = malloc(numberOfRattleConstraintsTerms*sizeof(cl_int4));
	clfloat = malloc(numberOfRattleConstraintsTerms*sizeof(cl_float));
        for( i=0; i<forceField->molecule.nAtoms;i++) rattleCounters[i] = 0;
	for(i=0;i<numberOfRattleConstraintsTerms;i++)
	{
		clint4[i].s[0] = (int) rattleConstraintsTerms[0][i];
		clint4[i].s[1] = (int) rattleConstraintsTerms[1][i];
		clint4[i].s[2] = rattleCounters[clint4[i].s[0]]++;
		clint4[i].s[3] = rattleCounters[clint4[i].s[1]]++;
		
		clfloat[i] = rattleConstraintsTerms[2][i];
	}
	forceField->rattleConstraintsIndexCL = clCreateBuffer(clProp.context, CL_MEM_READ_ONLY, sizeof(cl_int4)*numberOfRattleConstraintsTerms, NULL, &ea);
	forceField->rattleConstraintsTermsCL = clCreateBuffer(clProp.context, CL_MEM_READ_ONLY, sizeof(cl_float)*numberOfRattleConstraintsTerms, NULL, &eb);
	if(ea!=CL_SUCCESS||eb!=CL_SUCCESS)
	{
		fprintf(forceField->logfile, "I cannot bonds Buffer\n");
		fflush(forceField->logfile);
		exit(1);
	}
	ea = clEnqueueWriteBuffer(clProp.command_queue, forceField->rattleConstraintsIndexCL, CL_TRUE, 0, sizeof(cl_int4) * numberOfRattleConstraintsTerms, clint4, 0, NULL, NULL);
	eb = clEnqueueWriteBuffer(clProp.command_queue, forceField->rattleConstraintsTermsCL, CL_TRUE, 0, sizeof(cl_float) * numberOfRattleConstraintsTerms, clfloat, 0, NULL, NULL);
	if(ea!=CL_SUCCESS||eb!=CL_SUCCESS)
	{
		fprintf(forceField->logfile, "I cannot write bonds data in GPU Buffer\n");
		fflush(forceField->logfile);
		exit(1);
	}
	free(clint4);
	free(clfloat);
#ifdef DEBUG
	fprintf(forceField->logfile, "End numberOfRattleConstraintsTerms0\n");
	fflush(forceField->logfile);
#endif
	maxRattleCounters  = 1;
        for( i=0; i<forceField->molecule.nAtoms;i++) if(maxRattleCounters<rattleCounters[i]) maxRattleCounters = rattleCounters[i];
	free(rattleCounters);
	forceField->nBlockRattleBuffer = maxRattleCounters;
#ifdef DEBUG
	fprintf(forceField->logfile, "nBlockRattleBuffer=%d\n",maxRattleCounters);
	fflush(forceField->logfile);
#endif
	forceField->rattledeltaPositionBufferCL = clCreateBuffer(clProp.context, CL_MEM_READ_WRITE, sizeof(cl_float4)*forceField->nBlockRattleBuffer*forceField->molecule.nAtoms, NULL, &ec);
	forceField->rattledeltaVelocityBufferCL = clCreateBuffer(clProp.context, CL_MEM_READ_WRITE, sizeof(cl_float4)*forceField->nBlockRattleBuffer*forceField->molecule.nAtoms, NULL, &ec);
	forceField->rattleMovedBufferCL = clCreateBuffer(clProp.context, CL_MEM_READ_WRITE, sizeof(cl_int)*forceField->molecule.nAtoms, NULL, &ec);
	forceField->rattleUpdateBufferCL = clCreateBuffer(clProp.context, CL_MEM_READ_WRITE, sizeof(cl_int)*forceField->nBlockRattleBuffer*forceField->molecule.nAtoms, NULL, &ec);
	forceField->rattleDoneBufferCL = clCreateBuffer(clProp.context, CL_MEM_READ_WRITE, sizeof(cl_int), NULL, &ec);
#ifdef DEBUG
	fprintf(forceField->logfile, "End numberOfRattleConstraintsTerms\n");
	fflush(forceField->logfile);
#endif
	}

	free(gradCounters);
	forceField->nBlockGradientBuffer = maxCounters;
	//forceField->nBlockGradientBuffer = forceField->molecule.nAtoms;
	forceField->gradientBufferCL = clCreateBuffer(clProp.context,  CL_MEM_READ_WRITE, sizeof(cl_float4)*forceField->molecule.nAtoms*forceField->nBlockGradientBuffer, NULL, &ea);
	//maxCounters *= forceField->molecule.nAtoms;
#ifdef DEBUG
	fprintf(forceField->logfile, "nBlockGradientBuffer=%d\n",forceField->nBlockGradientBuffer);
	fprintf(forceField->logfile, "nAtoms=%d\n",forceField->molecule.nAtoms);
	fflush(forceField->logfile);
#endif
	//maxCounters = forceField->molecule.nAtoms*forceField->molecule.nAtoms;
	//forceField->gradientBufferCL = clCreateBuffer(clProp.context,  CL_MEM_READ_WRITE, sizeof(cl_float4)*maxCounters, NULL, &ea);
	if(ea!=CL_SUCCESS)
	{
		fprintf(forceField->logfile, "I cannot create the gradient Buffer\n");
		fflush(forceField->logfile);
		exit(1);
	}
	forceField->gradientBufferCPU = malloc(forceField->molecule.nAtoms*sizeof(cl_float4));

	// create kernels
	forceField->initEnergyKernel =  clCreateKernel(forceField->programMM, "initEnergy", &err);
	clSetKernelArg(forceField->initEnergyKernel, 0, sizeof(cl_mem),  &forceField->energyBufferCL);
	clSetKernelArg(forceField->initEnergyKernel, 1, sizeof(cl_float),  &forceField->nMaxTerms);
#ifdef DEBUG
	fprintf(forceField->logfile, "err = %d\n",err);
	fflush(forceField->logfile);
#endif


	if(forceField->numberOfStretchTerms>0)
	{
	forceField->addEnergyBondAmberKernel = clCreateKernel(forceField->programMM, "addEnergyBondAmber", &err);
	clSetKernelArg(forceField->addEnergyBondAmberKernel, 0, sizeof(cl_mem),  &forceField->energyBufferCL);
	clSetKernelArg(forceField->addEnergyBondAmberKernel, 1, sizeof(cl_mem),  &forceField->atomsCL);
	clSetKernelArg(forceField->addEnergyBondAmberKernel, 2, sizeof(cl_mem),  &forceField->bondIndexCL);
	clSetKernelArg(forceField->addEnergyBondAmberKernel, 3, sizeof(cl_mem),  &forceField->bondTermsCL);
	clSetKernelArg(forceField->addEnergyBondAmberKernel, 4, sizeof(cl_int),  &forceField->molecule.nAtoms);
	clSetKernelArg(forceField->addEnergyBondAmberKernel, 5, sizeof(cl_int),  &forceField->numberOfStretchTerms);
#ifdef DEBUG
	fprintf(forceField->logfile, "addEnergyBondAmber err = %d\n",err);
	fflush(forceField->logfile);
#endif
	size_t global = forceField->numberOfStretchTerms;
#ifdef DEBUG
	fprintf(forceField->logfile, "Call addEnergyBondAmber =====>0\n");
	fflush(forceField->logfile);
#endif
	err = clEnqueueNDRangeKernel(clProp.command_queue, forceField->addEnergyBondAmberKernel, 1, NULL, &global, NULL, 0, NULL, NULL);
	if(err != CL_SUCCESS)
	{
		printErrorCLRun(err);
        	fprintf(forceField->logfile, "I cannot execute addEnergyBondAmberKernel\n");
		fflush(forceField->logfile);
		exit(1);
	}
	}

	if(forceField->numberOfBendTerms>0)
	{
	forceField->addEnergyBendAmberKernel = clCreateKernel(forceField->programMM, "addEnergyBendAmber", &err);
	clSetKernelArg(forceField->addEnergyBendAmberKernel, 0, sizeof(cl_mem),  &forceField->energyBufferCL);
	clSetKernelArg(forceField->addEnergyBendAmberKernel, 1, sizeof(cl_mem),  &forceField->atomsCL);
	clSetKernelArg(forceField->addEnergyBendAmberKernel, 2, sizeof(cl_mem),  &forceField->bendIndexCL);
	clSetKernelArg(forceField->addEnergyBendAmberKernel, 3, sizeof(cl_mem),  &forceField->bendTermsCL);
	clSetKernelArg(forceField->addEnergyBendAmberKernel, 4, sizeof(cl_int),  &forceField->molecule.nAtoms);
	clSetKernelArg(forceField->addEnergyBendAmberKernel, 5, sizeof(cl_int),  &forceField->numberOfBendTerms);
#ifdef DEBUG
	fprintf(forceField->logfile, "addEnergyBendAmber err = %d\n",err);
	fflush(forceField->logfile);
#endif
	}

	if(forceField->numberOfDihedralTerms>0)
	{
	forceField->addEnergyDihedralAmberKernel =  clCreateKernel(forceField->programMM, "addEnergyDihedralAmber", &err);
	clSetKernelArg(forceField->addEnergyDihedralAmberKernel, 0, sizeof(cl_mem),  &forceField->energyBufferCL);
	clSetKernelArg(forceField->addEnergyDihedralAmberKernel, 1, sizeof(cl_mem),  &forceField->atomsCL);
	clSetKernelArg(forceField->addEnergyDihedralAmberKernel, 2, sizeof(cl_mem),  &forceField->dihedralIndexCL);
	clSetKernelArg(forceField->addEnergyDihedralAmberKernel, 3, sizeof(cl_mem),  &forceField->dihedralTermsCL);
	clSetKernelArg(forceField->addEnergyDihedralAmberKernel, 4, sizeof(cl_int),  &forceField->molecule.nAtoms);
	clSetKernelArg(forceField->addEnergyDihedralAmberKernel, 5, sizeof(cl_int),  &forceField->numberOfDihedralTerms);
#ifdef DEBUG
	fprintf(forceField->logfile, "addEnergyDihedralAmber err = %d\n",err);
	fflush(forceField->logfile);
#endif
	}

	if(numberOfVdw612>0)
	{
#ifdef DEBUG
	fprintf(forceField->logfile, "Use Coulomb = %d\n",forceField->options.coulomb);
	fflush(forceField->logfile);
#endif
	forceField->addEnergyVdw612AmberKernel =  clCreateKernel(forceField->programMM, "addEnergyVdw612Amber", &err);
	clSetKernelArg(forceField->addEnergyVdw612AmberKernel, 0, sizeof(cl_mem),  &forceField->energyBufferCL);
	clSetKernelArg(forceField->addEnergyVdw612AmberKernel, 1, sizeof(cl_mem),  &forceField->atomsCL);
	clSetKernelArg(forceField->addEnergyVdw612AmberKernel, 2, sizeof(cl_mem),  &forceField->vdw612IndexCL);
	clSetKernelArg(forceField->addEnergyVdw612AmberKernel, 3, sizeof(cl_mem),  &forceField->vdw612TermsCL);
	clSetKernelArg(forceField->addEnergyVdw612AmberKernel, 4, sizeof(cl_int),  &forceField->molecule.nAtoms);
	clSetKernelArg(forceField->addEnergyVdw612AmberKernel, 5, sizeof(cl_int),  &forceField->numberOfVdw612);
	clSetKernelArg(forceField->addEnergyVdw612AmberKernel, 6, sizeof(cl_int),  &useCoulomb);
#ifdef DEBUG
	fprintf(forceField->logfile, "addEnergyVdw612Amber err = %d\n",err);
	fflush(forceField->logfile);
#endif
	}
	if(numberOfSuttonChen>0)
	{
	forceField->addEnergySuttonChenAmberKernel =  clCreateKernel(forceField->programMM, "addEnergySuttonChenAmber", &err);
	clSetKernelArg(forceField->addEnergySuttonChenAmberKernel, 0, sizeof(cl_mem),  &forceField->energyBufferCL);
	clSetKernelArg(forceField->addEnergySuttonChenAmberKernel, 1, sizeof(cl_mem),  &forceField->atomsCL);
	clSetKernelArg(forceField->addEnergySuttonChenAmberKernel, 2, sizeof(cl_mem),  &forceField->suttonChenIndexCL);
	clSetKernelArg(forceField->addEnergySuttonChenAmberKernel, 3, sizeof(cl_mem),  &forceField->suttonChenTermsCL);
	clSetKernelArg(forceField->addEnergySuttonChenAmberKernel, 4, sizeof(cl_mem),  &forceField->rhoSuttonChenCL);
	clSetKernelArg(forceField->addEnergySuttonChenAmberKernel, 5, sizeof(cl_int),  &forceField->molecule.nAtoms);
	clSetKernelArg(forceField->addEnergySuttonChenAmberKernel, 6, sizeof(cl_int),  &forceField->numberOfSuttonChen);
#ifdef DEBUG
	fprintf(forceField->logfile, "addEnergyVdw612Amber err = %d\n",err);
	fflush(forceField->logfile);
#endif
	}

	if(forceField->numberOfHydrogenBonded>0)
	{
	forceField->addEnergyHydrogenBondedAmberKernel =  clCreateKernel(forceField->programMM, "addEnergyHydrogenBondedAmber", &err);
	clSetKernelArg(forceField->addEnergyHydrogenBondedAmberKernel, 0, sizeof(cl_mem),  &forceField->energyBufferCL);
	clSetKernelArg(forceField->addEnergyHydrogenBondedAmberKernel, 1, sizeof(cl_mem),  &forceField->atomsCL);
	clSetKernelArg(forceField->addEnergyHydrogenBondedAmberKernel, 2, sizeof(cl_mem),  &forceField->hydrogenBondedIndexCL);
	clSetKernelArg(forceField->addEnergyHydrogenBondedAmberKernel, 3, sizeof(cl_mem),  &forceField->hydrogenBondedTermsCL);
	clSetKernelArg(forceField->addEnergyHydrogenBondedAmberKernel, 4, sizeof(cl_int),  &forceField->molecule.nAtoms);
	clSetKernelArg(forceField->addEnergyHydrogenBondedAmberKernel, 5, sizeof(cl_int),  &forceField->numberOfHydrogenBonded);
#ifdef DEBUG
	fprintf(forceField->logfile, "addEnergyHydrogenBondedAmber err = %d\n",err);
	fflush(forceField->logfile);
#endif
	}

	if(numberOfPairWise>0)
	{
	forceField->addEnergyPairWiseKernel =  clCreateKernel(forceField->programMM, "addEnergyPairWise", &err);
	clSetKernelArg(forceField->addEnergyPairWiseKernel, 0, sizeof(cl_mem),  &forceField->energyBufferCL);
	clSetKernelArg(forceField->addEnergyPairWiseKernel, 1, sizeof(cl_mem),  &forceField->atomsCL);
	clSetKernelArg(forceField->addEnergyPairWiseKernel, 2, sizeof(cl_mem),  &forceField->pairWiseIndexCL);
	clSetKernelArg(forceField->addEnergyPairWiseKernel, 3, sizeof(cl_mem),  &forceField->pairWiseTermsCL);
	clSetKernelArg(forceField->addEnergyPairWiseKernel, 4, sizeof(cl_int),  &forceField->molecule.nAtoms);
	clSetKernelArg(forceField->addEnergyPairWiseKernel, 5, sizeof(cl_int),  &forceField->numberOfPairWise);
	clSetKernelArg(forceField->addEnergyPairWiseKernel, 6, sizeof(cl_int),  &forceField->options.coulomb);
	clSetKernelArg(forceField->addEnergyPairWiseKernel, 7, sizeof(cl_int),  &forceField->options.vanderWals);
#ifdef DEBUG
	fprintf(forceField->logfile, "addEnergyPairWise err = %d\n",err);
	fflush(forceField->logfile);
#endif
	}
	if(numberOfSuttonChen>0)
	{
	forceField->reduceRhoKernel =  clCreateKernel(forceField->programMM, "reduceRho", &err);
	clSetKernelArg(forceField->reduceRhoKernel, 0, sizeof(cl_mem),  &forceField->rhoBufferCL);
	clSetKernelArg(forceField->reduceRhoKernel, 1, sizeof(cl_int),  &forceField->molecule.nAtoms);
	clSetKernelArg(forceField->reduceRhoKernel, 2, sizeof(cl_int),  &forceField->nBlockRhoBuffer);
#ifdef DEBUG
	fprintf(forceField->logfile, "reduceRhoKernel err = %d\n",err);
	fflush(forceField->logfile);
#endif
	forceField->initRhoKernel =  clCreateKernel(forceField->programMM, "initRho", &err);
	clSetKernelArg(forceField->initRhoKernel, 0, sizeof(cl_mem),  &forceField->rhoBufferCL);
	clSetKernelArg(forceField->initRhoKernel, 1, sizeof(cl_int),  &forceField->molecule.nAtoms);
	clSetKernelArg(forceField->initRhoKernel, 2, sizeof(cl_int),  &forceField->nBlockRhoBuffer);
#ifdef DEBUG
	fprintf(forceField->logfile, "initRhoKernel err = %d\n",err);
	fflush(forceField->logfile);
#endif
	forceField->computeRhoKernel =  clCreateKernel(forceField->programMM, "computeRhoKernel", &err);
	clSetKernelArg(forceField->computeRhoKernel, 0, sizeof(cl_mem),  &forceField->rhoBufferCL);
	clSetKernelArg(forceField->computeRhoKernel, 1, sizeof(cl_int),  &forceField->molecule.nAtoms);
	clSetKernelArg(forceField->computeRhoKernel, 2, sizeof(cl_int),  &forceField->nBlockRhoBuffer);
#ifdef DEBUG
	fprintf(forceField->logfile, "computeRhoKernel err = %d\n",err);
	fflush(forceField->logfile);
#endif
	}

	forceField->reduceGradientsKernel =  clCreateKernel(forceField->programMM, "reduceGradients", &err);
	clSetKernelArg(forceField->reduceGradientsKernel, 0, sizeof(cl_mem),  &forceField->gradientBufferCL);
	clSetKernelArg(forceField->reduceGradientsKernel, 1, sizeof(cl_int),  &forceField->molecule.nAtoms);
	clSetKernelArg(forceField->reduceGradientsKernel, 2, sizeof(cl_int),  &forceField->nBlockGradientBuffer);
#ifdef DEBUG
	fprintf(forceField->logfile, "reduceGradients err = %d\n",err);
	fflush(forceField->logfile);
#endif

	forceField->initGradientsKernel =  clCreateKernel(forceField->programMM, "initGradients", &err);
	clSetKernelArg(forceField->initGradientsKernel, 0, sizeof(cl_mem),  &forceField->gradientBufferCL);
	clSetKernelArg(forceField->initGradientsKernel, 1, sizeof(cl_int),  &forceField->molecule.nAtoms);
	clSetKernelArg(forceField->initGradientsKernel, 2, sizeof(cl_int),  &forceField->nBlockGradientBuffer);
#ifdef DEBUG
	fprintf(forceField->logfile, "initGradients err = %d\n",err);
	fflush(forceField->logfile);
#endif

	forceField->initVelocitiesKernel =  clCreateKernel(forceField->programMM, "initVelocities", &err);
	clSetKernelArg(forceField->initVelocitiesKernel, 0, sizeof(cl_mem),  &forceField->atomsCL);
	clSetKernelArg(forceField->initVelocitiesKernel, 1, sizeof(cl_int),  &forceField->molecule.nAtoms);
#ifdef DEBUG
	fprintf(forceField->logfile, "initVelocities err = %d\n",err);
	fflush(forceField->logfile);
#endif
	//global = forceField->molecule.nAtoms;
	//clEnqueueNDRangeKernel(clProp.command_queue, forceField->initVelocitiesKernel, 1, NULL, &global, &local, 0, NULL, NULL);
	//clFinish(clProp.command_queue);

	if(forceField->numberOfStretchTerms>0)
	{
	forceField->addGradientBondAmberKernel =  clCreateKernel(forceField->programMM, "addGradientBondAmber", &err);
	if(err!=CL_SUCCESS)
	{
		fprintf(forceField->logfile, "I cannot create the addGradientBondAmberKernel kernel\n");
		fflush(forceField->logfile);
		exit(1);
	}
	err |= clSetKernelArg(forceField->addGradientBondAmberKernel, 0, sizeof(cl_mem),  &forceField->energyBufferCL);
	err |= clSetKernelArg(forceField->addGradientBondAmberKernel, 1, sizeof(cl_mem),  &forceField->gradientBufferCL);
	err |= clSetKernelArg(forceField->addGradientBondAmberKernel, 2, sizeof(cl_mem),  &forceField->atomsCL);
	err |= clSetKernelArg(forceField->addGradientBondAmberKernel, 3, sizeof(cl_mem),  &forceField->bondIndexCL);
	err |= clSetKernelArg(forceField->addGradientBondAmberKernel, 4, sizeof(cl_mem),  &forceField->bondTermsCL);
	err |= clSetKernelArg(forceField->addGradientBondAmberKernel, 5, sizeof(cl_int),  &forceField->molecule.nAtoms);
	err |= clSetKernelArg(forceField->addGradientBondAmberKernel, 6, sizeof(cl_int),  &forceField->numberOfStretchTerms);
#ifdef DEBUG
	fprintf(forceField->logfile, "addGradientBondAmber err = %d\n",err);
	fprintf(forceField->logfile, "Call addGradientBondAmberKernel =====>0\n");
	fflush(forceField->logfile);
#endif
	size_t global = forceField->numberOfStretchTerms;
	err = clEnqueueNDRangeKernel(clProp.command_queue, forceField->addGradientBondAmberKernel, 1, NULL, &global, NULL, 0, NULL, NULL);
	if(err != CL_SUCCESS)
	{
		printErrorCLRun(err);
        	fprintf(forceField->logfile, "I cannot execute addGradientBondAmberKernel\n");
		fflush(forceField->logfile);
		exit(1);
	}
	}

	if(forceField->numberOfBendTerms>0)
	{
	forceField->addGradientBendAmberKernel =  clCreateKernel(forceField->programMM, "addGradientBendAmber", &err);
	clSetKernelArg(forceField->addGradientBendAmberKernel, 0, sizeof(cl_mem),  &forceField->energyBufferCL);
	clSetKernelArg(forceField->addGradientBendAmberKernel, 1, sizeof(cl_mem),  &forceField->gradientBufferCL);
	clSetKernelArg(forceField->addGradientBendAmberKernel, 2, sizeof(cl_mem),  &forceField->atomsCL);
	clSetKernelArg(forceField->addGradientBendAmberKernel, 3, sizeof(cl_mem),  &forceField->bendIndexCL);
	clSetKernelArg(forceField->addGradientBendAmberKernel, 4, sizeof(cl_mem),  &forceField->bendTermsCL);
	clSetKernelArg(forceField->addGradientBendAmberKernel, 5, sizeof(cl_int),  &forceField->molecule.nAtoms);
	clSetKernelArg(forceField->addGradientBendAmberKernel, 6, sizeof(cl_int),  &forceField->numberOfBendTerms);
#ifdef DEBUG
	fprintf(forceField->logfile, "addGradientBendAmber err = %d\n",err);
	fflush(forceField->logfile);
#endif
	}

	if(forceField->numberOfDihedralTerms>0)
	{
	forceField->addGradientDihedralAmberKernel =  clCreateKernel(forceField->programMM, "addGradientDihedralAmber", &err);
	clSetKernelArg(forceField->addGradientDihedralAmberKernel, 0, sizeof(cl_mem),  &forceField->energyBufferCL);
	clSetKernelArg(forceField->addGradientDihedralAmberKernel, 1, sizeof(cl_mem),  &forceField->gradientBufferCL);
	clSetKernelArg(forceField->addGradientDihedralAmberKernel, 2, sizeof(cl_mem),  &forceField->atomsCL);
	clSetKernelArg(forceField->addGradientDihedralAmberKernel, 3, sizeof(cl_mem),  &forceField->dihedralIndexCL);
	clSetKernelArg(forceField->addGradientDihedralAmberKernel, 4, sizeof(cl_mem),  &forceField->dihedralTermsCL);
	clSetKernelArg(forceField->addGradientDihedralAmberKernel, 5, sizeof(cl_int),  &forceField->molecule.nAtoms);
	clSetKernelArg(forceField->addGradientDihedralAmberKernel, 6, sizeof(cl_int),  &forceField->numberOfDihedralTerms);
#ifdef DEBUG
	fprintf(forceField->logfile, "addGradientDihedralAmber err = %d\n",err);
	fflush(forceField->logfile);
#endif
	}

	if(forceField->numberOfImproperTorsionTerms>0)
	{
	forceField->addGradientImproperTorsionKernel =  clCreateKernel(forceField->programMM, "addGradientImproperTorsion", &err);
#ifdef DEBUG
	fprintf(forceField->logfile, "addGradientImproperTorsion err = %d\n",err);
	fflush(forceField->logfile);
#endif
	}

	if(numberOfVdw612>0)
	{
	forceField->addGradientVdw612AmberKernel =  clCreateKernel(forceField->programMM, "addGradientVdw612Amber", &err);
	clSetKernelArg(forceField->addGradientVdw612AmberKernel, 0, sizeof(cl_mem),  &forceField->energyBufferCL);
	clSetKernelArg(forceField->addGradientVdw612AmberKernel, 1, sizeof(cl_mem),  &forceField->gradientBufferCL);
	clSetKernelArg(forceField->addGradientVdw612AmberKernel, 2, sizeof(cl_mem),  &forceField->atomsCL);
	clSetKernelArg(forceField->addGradientVdw612AmberKernel, 3, sizeof(cl_mem),  &forceField->vdw612IndexCL);
	clSetKernelArg(forceField->addGradientVdw612AmberKernel, 4, sizeof(cl_mem),  &forceField->vdw612TermsCL);
	clSetKernelArg(forceField->addGradientVdw612AmberKernel, 5, sizeof(cl_int),  &forceField->molecule.nAtoms);
	clSetKernelArg(forceField->addGradientVdw612AmberKernel, 6, sizeof(cl_int),  &forceField->numberOfVdw612);
	clSetKernelArg(forceField->addGradientVdw612AmberKernel, 7, sizeof(cl_int),  &forceField->options.coulomb);
#ifdef DEBUG
	fprintf(forceField->logfile, "addGradientVdw612Amber err = %d\n",err);
	fflush(forceField->logfile);
#endif
	}
	if(numberOfVdw714>0)
	{
	forceField->addGradientVdw714AmberKernel =  clCreateKernel(forceField->programMM, "addGradientVdw714Amber", &err);
	clSetKernelArg(forceField->addGradientVdw714AmberKernel, 0, sizeof(cl_mem),  &forceField->energyBufferCL);
	clSetKernelArg(forceField->addGradientVdw714AmberKernel, 1, sizeof(cl_mem),  &forceField->gradientBufferCL);
	clSetKernelArg(forceField->addGradientVdw714AmberKernel, 2, sizeof(cl_mem),  &forceField->atomsCL);
	clSetKernelArg(forceField->addGradientVdw714AmberKernel, 3, sizeof(cl_mem),  &forceField->vdw714IndexCL);
	clSetKernelArg(forceField->addGradientVdw714AmberKernel, 4, sizeof(cl_mem),  &forceField->vdw714TermsCL);
	clSetKernelArg(forceField->addGradientVdw714AmberKernel, 5, sizeof(cl_int),  &forceField->molecule.nAtoms);
	clSetKernelArg(forceField->addGradientVdw714AmberKernel, 6, sizeof(cl_int),  &forceField->numberOfVdw714);
	clSetKernelArg(forceField->addGradientVdw714AmberKernel, 7, sizeof(cl_int),  &forceField->options.coulomb);
#ifdef DEBUG
	fprintf(forceField->logfile, "addGradientVdw714Amber err = %d\n",err);
	fflush(forceField->logfile);
#endif
	}
	if(numberOfSuttonChen>0)
	{
	forceField->addGradientSuttonChenAmberKernel =  clCreateKernel(forceField->programMM, "addGradientSuttonChenAmber", &err);
	clSetKernelArg(forceField->addGradientSuttonChenAmberKernel, 0, sizeof(cl_mem),  &forceField->energyBufferCL);
	clSetKernelArg(forceField->addGradientSuttonChenAmberKernel, 1, sizeof(cl_mem),  &forceField->gradientBufferCL);
	clSetKernelArg(forceField->addGradientSuttonChenAmberKernel, 2, sizeof(cl_mem),  &forceField->atomsCL);
	clSetKernelArg(forceField->addGradientSuttonChenAmberKernel, 3, sizeof(cl_mem),  &forceField->suttonChenIndexCL);
	clSetKernelArg(forceField->addGradientSuttonChenAmberKernel, 4, sizeof(cl_mem),  &forceField->suttonChenTermsCL);
	clSetKernelArg(forceField->addGradientSuttonChenAmberKernel, 6, sizeof(cl_mem),  &forceField->suttonChenRhoCL);
	clSetKernelArg(forceField->addGradientSuttonChenAmberKernel, 6, sizeof(cl_int),  &forceField->molecule.nAtoms);
	clSetKernelArg(forceField->addGradientSuttonChenAmberKernel, 7, sizeof(cl_int),  &forceField->numberOfVdw612);
#ifdef DEBUG
	fprintf(forceField->logfile, "addGradientSuttonChenAmberKernel err = %d\n",err);
	fflush(forceField->logfile);
#endif
	}

	if(forceField->numberOfHydrogenBonded>0)
	{
	forceField->addGradientHydrogenBondedAmberKernel =  clCreateKernel(forceField->programMM, "addGradientHydrogenBondedAmber", &err);
	clSetKernelArg(forceField->addGradientHydrogenBondedAmberKernel, 0, sizeof(cl_mem),  &forceField->energyBufferCL);
	clSetKernelArg(forceField->addGradientHydrogenBondedAmberKernel, 1, sizeof(cl_mem),  &forceField->gradientBufferCL);
	clSetKernelArg(forceField->addGradientHydrogenBondedAmberKernel, 2, sizeof(cl_mem),  &forceField->atomsCL);
	clSetKernelArg(forceField->addGradientHydrogenBondedAmberKernel, 3, sizeof(cl_mem),  &forceField->hydrogenBondedIndexCL);
	clSetKernelArg(forceField->addGradientHydrogenBondedAmberKernel, 4, sizeof(cl_mem),  &forceField->hydrogenBondedTermsCL);
	clSetKernelArg(forceField->addGradientHydrogenBondedAmberKernel, 5, sizeof(cl_int),  &forceField->molecule.nAtoms);
	clSetKernelArg(forceField->addGradientHydrogenBondedAmberKernel, 6, sizeof(cl_int),  &forceField->numberOfHydrogenBonded);
#ifdef DEBUG
	fprintf(forceField->logfile, "addGradientHydrogenBondedAmber err = %d\n",err);
	fflush(forceField->logfile);
#endif
	}

	if(numberOfPairWise>0)
	{
	forceField->addGradientPairWiseKernel =  clCreateKernel(forceField->programMM, "addGradientPairWise", &err);
	clSetKernelArg(forceField->addGradientPairWiseKernel, 0, sizeof(cl_mem),  &forceField->energyBufferCL);
	clSetKernelArg(forceField->addGradientPairWiseKernel, 1, sizeof(cl_mem),  &forceField->gradientBufferCL);
	clSetKernelArg(forceField->addGradientPairWiseKernel, 2, sizeof(cl_mem),  &forceField->atomsCL);
	clSetKernelArg(forceField->addGradientPairWiseKernel, 3, sizeof(cl_mem),  &forceField->pairWiseIndexCL);
	clSetKernelArg(forceField->addGradientPairWiseKernel, 4, sizeof(cl_mem),  &forceField->pairWiseTermsCL);
	clSetKernelArg(forceField->addGradientPairWiseKernel, 5, sizeof(cl_int),  &forceField->molecule.nAtoms);
	clSetKernelArg(forceField->addGradientPairWiseKernel, 6, sizeof(cl_int),  &forceField->numberOfPairWise);
	clSetKernelArg(forceField->addGradientPairWiseKernel, 7, sizeof(cl_int),  &forceField->options.coulomb);
	clSetKernelArg(forceField->addGradientPairWiseKernel, 8, sizeof(cl_int),  &forceField->options.vanderWals);
#ifdef DEBUG
	fprintf(forceField->logfile, " addGradientPairWise err = %d\n",err);
	fflush(forceField->logfile);
#endif
	}



	// set the argument list for the kernel command
	// cl_float4 pos => pos.s[0], pos.s[1], ... 
};
#endif
/**********************************************************************/
ForceField createAmberModel (Molecule* mol, ForceFieldOptions forceFieldOptions,FILE* logfile)
{
	ForceField forceField = newAmberModel();

	forceField.molecule = *(mol->klass->copy(mol));
	
	forceField.options = forceFieldOptions;
	forceField.logfile = logfile;
	setRattleConstraintsParameters(&forceField);

	fprintf(forceField.logfile, ("Setting of Parameters ...\n"));
	fflush(forceField.logfile);

	setAmberParameters(&forceField);
#ifdef ENABLE_CL
	initCLForceField (&forceField);
#endif

	return forceField;
}
/**********************************************************************/
ForceField createPairWiseModel(Molecule* mol, ForceFieldOptions forceFieldOptions, FILE* logfile)
{
	ForceField forceField = newPairWiseModel();
	

	forceField.molecule = *(mol->klass->copy(mol));
	
	forceField.options = forceFieldOptions;
	forceField.logfile = logfile;
	setRattleConstraintsParameters(&forceField);

	forceField.options.bondStretch = FALSE;
	forceField.options.angleBend = FALSE;
	forceField.options.dihedralAngle = FALSE;
	forceField.options.improperTorsion = FALSE;
	forceField.options.vdw612 = FALSE;
	forceField.options.hydrogenBonded612 = FALSE;
	forceField.options.hydrogenBonded1012 = FALSE;
	forceField.options.hydrogenBondedMorse = FALSE;


	fprintf(forceField.logfile, ("Setting of Parameters ...\n"));
	fflush(forceField.logfile);

	setAllPairWiseParameters(&forceField);

	return forceField;
}
/**********************************************************************/
void saveAmberParameters()
{
	createMMFile();
}
/**********************************************************************/
void loadAmberParameters()
{
	AmberParameters amberParameters;

	char* defaultFileName = getFileNameParameters();

	if(staticAmberParameters != NULL) freeAmberParameters(staticAmberParameters);

	amberParameters =  newAmberParameters();

	if(!readAmberParameters(&amberParameters,defaultFileName))
	{
		createMMFile();
		if(!readAmberParameters(&amberParameters,defaultFileName))
		{
			free(defaultFileName);
			return;
		}
	}

	staticAmberParameters = malloc(sizeof(AmberParameters));
	*staticAmberParameters = amberParameters;

	free(defaultFileName);
}
/**********************************************************************/
AmberParameters* getPointerAmberParameters()
{
	return staticAmberParameters;
}
/**********************************************************************/
void setPointerAmberParameters(AmberParameters* ptr)
{
	staticAmberParameters = ptr;
}
/********************************************************************************/
char** getListMMTypes(int* nlist)
{

	char** t = NULL;
	
	int i;

	*nlist = 0;

	if(!staticAmberParameters || staticAmberParameters->numberOfTypes<=0)
		return NULL;
	
	t = malloc(staticAmberParameters->numberOfTypes*sizeof(char*));
	*nlist = staticAmberParameters->numberOfTypes;

	for(i=0;i<staticAmberParameters->numberOfTypes;i++)
		t[i] = strdup(staticAmberParameters->atomTypes[i].name);

	return t;
}
/**********************************************************************/
static void printEnergies(ForceField* forceField)
{
	double totalEnergy = 0.0;
	double energy = 0.0;
	Molecule* molecule = &forceField->molecule;
	int j;

#ifdef DEBUG
        TimerType timer;
#endif

#ifdef ENABLE_CL
	CLProp clProp = getCLProp();
	size_t local = 64;
	size_t global;
	int i;

	fprintf(forceField->logfile, "Begin calculateEnergyTmpAmber\n");
	fflush(forceField->logfile);
	global = forceField->nMaxTerms;
	clEnqueueNDRangeKernel(clProp.command_queue, forceField->initEnergyKernel, 1, NULL, &global, NULL, 0, NULL, NULL);
	clFinish(clProp.command_queue);

#ifdef DEBUG
        timer_init(timer);
       	timer_start( timer );
#endif
	if(forceField->numberOfStretchTerms>0)
	{
	global = forceField->numberOfStretchTerms;
	clEnqueueNDRangeKernel(clProp.command_queue, forceField->addEnergyBondAmberKernel, 1, NULL, &global, NULL, 0, NULL, NULL);
	clFinish(clProp.command_queue);
	}

	if(forceField->numberOfBendTerms>0)
	{
	global = forceField->numberOfBendTerms;
	clEnqueueNDRangeKernel(clProp.command_queue, forceField->addEnergyBendAmberKernel, 1, NULL, &global, NULL, 0, NULL, NULL);
	clFinish(clProp.command_queue);
	}

	if(forceField->numberOfDihedralTerms>0)
	{
	global = forceField->numberOfDihedralTerms;
	clEnqueueNDRangeKernel(clProp.command_queue, forceField->addEnergyDihedralAmberKernel, 1, NULL, &global, NULL, 0, NULL, NULL);
	clFinish(clProp.command_queue);
	}

	if(forceField->numberOfVdw612>0)
	{
	global = forceField->numberOfVdw612;
	clEnqueueNDRangeKernel(clProp.command_queue, forceField->addEnergyVdw612AmberKernel, 1, NULL, &global, NULL, 0, NULL, NULL);
	clFinish(clProp.command_queue);
	}

	if(forceField->numberOfSuttonChen>0)
	{
	global = forceField->numberOfSuttonChen;
	clEnqueueNDRangeKernel(clProp.command_queue, forceField->initRhoKernel, 1, NULL, &global, NULL, 0, NULL, NULL);
	clEnqueueNDRangeKernel(clProp.command_queue, forceField->computeRhoKernel, 1, NULL, &global, NULL, 0, NULL, NULL);
	global = forceField->molecule.nAtoms;
	clEnqueueNDRangeKernel(clProp.command_queue, forceField->reduceRhoKernel, 1, NULL, &global, NULL, 0, NULL, NULL);
	global = forceField->numberOfSuttonChen;
	clEnqueueNDRangeKernel(clProp.command_queue, forceField->addEnergySuttoChenAmberKernel, 1, NULL, &global, NULL, 0, NULL, NULL);
	clFinish(clProp.command_queue);
	}

	if(forceField->numberOfHydrogenBonded>0)
	{
	global = forceField->numberOfHydrogenBonded;
	clEnqueueNDRangeKernel(clProp.command_queue, forceField->addEnergyHydrogenBondedAmberKernel, 1, NULL, &global, NULL, 0, NULL, NULL);
	clFinish(clProp.command_queue);
	}

	if(forceField->numberOfPairWise>0)
	{
	global = forceField->numberOfPairWise;
	clEnqueueNDRangeKernel(clProp.command_queue, forceField->addEnergyPairWiseKernel, 1, NULL, &global, NULL, 0, NULL, NULL);
	clFinish(clProp.command_queue);
	}


	clEnqueueReadBuffer(clProp.command_queue, forceField->energyBufferCL, CL_TRUE, 0, sizeof(cl_float)*forceField->nMaxTerms, forceField->energyBufferCPU, 0, NULL, NULL);
	energy = 0;
	for(i=0;i<forceField->nMaxTerms;i++) energy+=forceField->energyBufferCPU[i];
#ifdef DEBUG
       	timer_stop(timer);
        fprintf(forceField->logfile, "time (s) = %f\n", timer_get(timer)*1e-6);
	fflush(forceField->logfile);
#endif
#else // ENABLE_CL
#ifdef DEBUG
        timer_init(timer);
       	timer_start( timer );
#endif
#endif

	for(j=0;j<80;j++) printf("-"); printf("\n");
	energy =calculateEnergyBondAmber(forceField,molecule);
	totalEnergy += energy;
	printf("%-30s %s %20.14f\n","Bond energy","=",energy);
	energy =calculateEnergyBendAmber(forceField,molecule);
	totalEnergy += energy;
	printf("%-30s %s %20.14f\n","Bend energy","=",energy);
	energy =calculateEnergyStrBendAmber(forceField,molecule);
	totalEnergy += energy;
	printf("%-30s %s %20.14f\n","Str-Bend energy","=",energy);
	energy =calculateEnergyDihedralAmber(forceField,molecule);
	totalEnergy += energy;
	printf("%-30s %s %20.14f\n","Dihedral energy","=",energy);
	energy =calculateEnergyImproperTorsionAmber(forceField,molecule);
	totalEnergy += energy;
	printf("%-30s %s %20.14f\n","Improper Torsion energy","=",energy);
	energy =calculateEnergyOutOfPlaneAmber(forceField,molecule);
	totalEnergy += energy;
	printf("%-30s %s %20.14f\n","Out Of Plane energy","=",energy);
	energy =calculateEnergyVdw612Amber(forceField,molecule);
	totalEnergy += energy;
	printf("%-30s %s %20.14f\n","VDW 6-12 energy","=",energy);
	energy =calculateEnergyVdw714Amber(forceField,molecule);
	totalEnergy += energy;
	printf("%-30s %s %20.14f\n","VDW 7-14 energy","=",energy);
	energy =calculateEnergyCoulombAmber(forceField,molecule);
	totalEnergy += energy;
	printf("%-30s %s %20.14f\n","Coulomb energy","=",energy);
	energy =calculateEnergySuttonChenAmber(forceField,molecule);
	totalEnergy += energy;
	printf("%-30s %s %20.14f\n","Sutton-Chen energy","=",energy);
	energy =calculateEnergyHydrogenBondedAmber(forceField,molecule);
	totalEnergy += energy;
	printf("%-30s %s %20.14f\n","Hydrogen Bonded energy","=",energy);
	energy =calculateEnergyPairWise(forceField,molecule);
	totalEnergy += energy;
	printf("%-30s %s %20.14f\n","Pair wise energy","=",energy);
	energy =getH4Energy(forceField, molecule);
	totalEnergy += energy;
	printf("%-30s %s %20.14f\n","Hydrogen Bonded-H4 energy","=",energy);

	printf("%-30s %s %20.14f\n","Total energy","=", totalEnergy);
	for(j=0;j<80;j++) printf("-"); printf("\n");
}
