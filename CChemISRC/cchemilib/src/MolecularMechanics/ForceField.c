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

/* ForceField.c */
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>

#include "../Utils/Utils.h"
#include "../Utils/CLProp.h"
#include "../Utils/AtomsProp.h"
#include "../Utils/Constants.h"
#include "../Utils/Types.h"
#include "../Molecule/Molecule.h"
#include "../MolecularMechanics/ForceField.h"
#include "../Utils/QL.h"

static void computeGradientNumeric(ForceField* forceField);
/**********************************************************************/
ForceField newForceField()
{
	int i;
	ForceField forceField;

	forceField.molecule = *(newMolecule());

	forceField.klass = malloc(sizeof(ForceFieldClass));
	forceField.klass->calculateGradient = NULL;
	forceField.klass->calculateGradientNumeric = computeGradientNumeric;
	forceField.klass->calculateEnergy = NULL;
	forceField.klass->calculateEnergyTmp = NULL;
	forceField.klass->printEnergies = NULL;


	forceField.numberOfRattleConstraintsTerms = 0;
	forceField.numberOfStretchTerms = 0;
	forceField.numberOfBendTerms = 0;
	forceField.numberOfStrBendTerms = 0;
	forceField.numberOfDihedralTerms = 0;
	forceField.numberOfImproperTorsionTerms = 0;
	forceField.numberOfOutOfPlaneTerms = 0;
	forceField.numberOfVdw612 = 0;
	forceField.numberOfVdw714 = 0;
	forceField.numberOfCoulomb = 0;
	forceField.numberOfHydrogenBonded = 0;
	forceField.numberOfSuttonChen = 0;

	for(i=0;i<STRETCHDIM;i++) forceField.bondStretchTerms[i] = NULL;

	for(i=0;i<BENDDIM;i++) forceField.angleBendTerms[i] = NULL;
	for(i=0;i<STRBENDDIM;i++) forceField.strBendTerms[i] = NULL;

	for(i=0;i<DIHEDRALDIM;i++) forceField.dihedralAngleTerms[i] = NULL;
	for(i=0;i<IMPROPERDIHEDRALDIM;i++) forceField.improperTorsionTerms[i] = NULL;
	for(i=0;i<OUTOFPLANEDIM;i++) forceField.outOfPlaneTerms[i] = NULL;

	for(i=0;i<VDW612DIM;i++) forceField.vdw612Terms[i] = NULL;
	for(i=0;i<VDW714DIM;i++) forceField.vdw714Terms[i] = NULL;
	for(i=0;i<COULOMBDIM;i++) forceField.coulombTerms[i] = NULL;
	for(i=0;i<HYDROGENBONDEDDIM;i++) forceField.hydrogenBondedTerms[i] = NULL;

	forceField.numberOfPairWise = 0;
	for(i=0;i<PAIRWISEDIM;i++)
		forceField.pairWiseTerms[i] = NULL;
	forceField.logfile = stdout;

	forceField.options.type = AMBER;
	forceField.options.coulomb = TRUE;
	forceField.options.numeric = FALSE;
	forceField.options.hydrogenBonded612 = FALSE;
	forceField.options.hydrogenBonded1012 = TRUE;
	forceField.options.hydrogenBondedMorse = FALSE;
	forceField.options.improperTorsion = TRUE;
	forceField.options.outOfPlane = TRUE;
	forceField.options.vdw612 = TRUE;
	forceField.options.vdw714 = TRUE;
	forceField.options.hbDirectional = FALSE;
	forceField.options.bondStretch = TRUE;
	forceField.options.angleBend = TRUE;
	forceField.options.rattleConstraints = NOCONSTRAINTS;
	forceField.options.dx = 1e-5;
	sprintf(forceField.options.chargesType,"DEFAULT");
	return forceField;

}
/**********************************************************************/
void freeForceField(ForceField* forceField)
{

	int i;
	forceField->molecule.klass->free(&forceField->molecule);

	if(forceField->klass != NULL)
	{
		free(forceField->klass);
		forceField->klass = NULL;
	}
	for(i=0;i<STRETCHDIM;i++)
		if(forceField->bondStretchTerms[i] !=NULL)
		{
			free(forceField->bondStretchTerms[i]);
			forceField->bondStretchTerms[i] = NULL;
		}
	for(i=0;i<BENDDIM;i++)
		if(forceField->angleBendTerms[i] != NULL)
		{
			free(forceField->angleBendTerms[i]);
			forceField->angleBendTerms[i] = NULL;
		}
	for(i=0;i<STRBENDDIM;i++)
		if(forceField->strBendTerms[i] != NULL)
		{
			free(forceField->strBendTerms[i]);
			forceField->strBendTerms[i] = NULL;
		}
	for(i=0;i<DIHEDRALDIM;i++)
		if(forceField->dihedralAngleTerms[i] != NULL)
		{
			free(forceField->dihedralAngleTerms[i]);
			forceField->dihedralAngleTerms[i] = NULL;
		}
	for(i=0;i<IMPROPERDIHEDRALDIM;i++)
	if(forceField->improperTorsionTerms[i] != NULL)
	{
		free(forceField->improperTorsionTerms[i]);
		forceField->improperTorsionTerms[i] = NULL;
	}

	for(i=0;i<OUTOFPLANEDIM;i++)
	if(forceField->outOfPlaneTerms[i] != NULL)
	{
		free(forceField->outOfPlaneTerms[i]);
		forceField->outOfPlaneTerms[i] = NULL;
	}

	for(i=0;i<VDW612DIM;i++)
		if(forceField->vdw612Terms[i] != NULL)
		{
			free(forceField->vdw612Terms[i]);
			forceField->vdw612Terms[i] = NULL;
		}
	for(i=0;i<VDW714DIM;i++)
		if(forceField->vdw714Terms[i] != NULL)
		{
			free(forceField->vdw714Terms[i]);
			forceField->vdw714Terms[i] = NULL;
		}


	for(i=0;i<COULOMBDIM;i++)
		if(forceField->coulombTerms[i] != NULL)
		{
			free(forceField->coulombTerms[i]);
			forceField->coulombTerms[i] = NULL;
		}
	for(i=0;i<HYDROGENBONDEDDIM;i++)
		if(forceField->hydrogenBondedTerms[i] != NULL)
		{
			free(forceField->hydrogenBondedTerms[i]);
			forceField->hydrogenBondedTerms[i] = NULL;
		}

	forceField->numberOfStretchTerms = 0;
	forceField->numberOfBendTerms = 0;
	forceField->numberOfStrBendTerms = 0;
	forceField->numberOfDihedralTerms = 0;
	forceField->numberOfImproperTorsionTerms = 0;
	forceField->numberOfVdw612 = 0;
	forceField->numberOfVdw714 = 0;
	forceField->numberOfCoulomb = 0;
	forceField->numberOfHydrogenBonded = 0;
	forceField->numberOfSuttonChen = 0;
	forceField->numberOfRattleConstraintsTerms = 0;

	for(i=0;i<PAIRWISEDIM;i++)
		if(forceField->pairWiseTerms[i] != NULL)
		{
			free(forceField->pairWiseTerms[i]);
			forceField->pairWiseTerms[i] = NULL;
		}
	forceField->numberOfPairWise = 0;
}
/*****************************************************************************/
ForceField copyForceField(ForceField* f)
{
	int i;
	int j;
	int k;
	ForceField forceField = newForceField();
#ifdef ENABLE_CL
	forceField = *f;
#endif

	forceField.molecule = *(f->molecule.klass->copy(&f->molecule));
	//printf("End copyMolecule\n");
	fflush(stdout);

	/* already in newForceField */
	/* forceField.klass = malloc(sizeof(ForceFieldClass));*/
	forceField.klass->calculateGradient = f->klass->calculateGradient;
	forceField.klass->calculateGradientNumeric = f->klass->calculateGradientNumeric;
	forceField.klass->calculateEnergy = f->klass->calculateEnergy;
	forceField.klass->calculateEnergyTmp = f->klass->calculateEnergyTmp;
	forceField.klass->printEnergies = f->klass->printEnergies;


	forceField.numberOfStretchTerms = f->numberOfStretchTerms;
	forceField.numberOfBendTerms = f->numberOfBendTerms;
	forceField.numberOfStrBendTerms = f->numberOfStrBendTerms;
	forceField.numberOfDihedralTerms = f->numberOfDihedralTerms;
	forceField.numberOfImproperTorsionTerms = f->numberOfImproperTorsionTerms;
	forceField.numberOfOutOfPlaneTerms = f->numberOfOutOfPlaneTerms;
	forceField.numberOfVdw612 = f->numberOfVdw612;
	forceField.numberOfVdw714 = f->numberOfVdw714;
	forceField.numberOfCoulomb = f->numberOfCoulomb;
	forceField.numberOfHydrogenBonded = f->numberOfHydrogenBonded;
	forceField.numberOfSuttonChen = f->numberOfSuttonChen;
	forceField.numberOfRattleConstraintsTerms = f->numberOfRattleConstraintsTerms;
	forceField.numberOfPairWise = f->numberOfPairWise;

	k = forceField.numberOfStretchTerms;
	if(k>0)
	for(i=0;i<STRETCHDIM;i++)
	{
		forceField.bondStretchTerms[i] = malloc(k*sizeof(double));
		for(j=0;j<k;j++) forceField.bondStretchTerms[i][j] = f->bondStretchTerms[i][j];
	}

	k = forceField.numberOfBendTerms;
	if(k>0)
	for(i=0;i<BENDDIM;i++)
	{
		forceField.angleBendTerms[i] = malloc(k*sizeof(double));
		for(j=0;j<k;j++) forceField.angleBendTerms[i][j] = f->angleBendTerms[i][j];
	}

	k = forceField.numberOfStrBendTerms;
	if(k>0)
	for(i=0;i<STRBENDDIM;i++)
	{
		forceField.strBendTerms[i] = malloc(k*sizeof(double));
		for(j=0;j<k;j++) forceField.strBendTerms[i][j] = f->strBendTerms[i][j];
	}


	k = forceField.numberOfDihedralTerms;
	if(k>0)
	for(i=0;i<DIHEDRALDIM;i++)
	{
		forceField.dihedralAngleTerms[i] = malloc(k*sizeof(double));
		for(j=0;j<k;j++) forceField.dihedralAngleTerms[i][j] = f->dihedralAngleTerms[i][j];
	}

	k = forceField.numberOfImproperTorsionTerms;
	if(k>0)
	for(i=0;i<IMPROPERDIHEDRALDIM;i++)
	{
		forceField.improperTorsionTerms[i] = malloc(k*sizeof(double));
		for(j=0;j<k;j++) forceField.improperTorsionTerms[i][j] = f->improperTorsionTerms[i][j];
	}

	k = forceField.numberOfOutOfPlaneTerms;
	if(k>0)
	for(i=0;i<OUTOFPLANEDIM;i++)
	{
		forceField.outOfPlaneTerms[i] = malloc(k*sizeof(double));
		for(j=0;j<k;j++) forceField.outOfPlaneTerms[i][j] = f->outOfPlaneTerms[i][j];
	}


	k = forceField.numberOfVdw612;
	if(k>0)
	for(i=0;i<VDW612DIM;i++)
	{
		forceField.vdw612Terms[i] = malloc(k*sizeof(double));
		for(j=0;j<k;j++) forceField.vdw612Terms[i][j] = f->vdw612Terms[i][j];
	}

	k = forceField.numberOfVdw714;
	if(k>0)
	for(i=0;i<VDW714DIM;i++)
	{
		forceField.vdw714Terms[i] = malloc(k*sizeof(double));
		for(j=0;j<k;j++) forceField.vdw714Terms[i][j] = f->vdw714Terms[i][j];
	}

	k = forceField.numberOfCoulomb;
	if(k>0)
	for(i=0;i<COULOMBDIM;i++)
	{
		forceField.coulombTerms[i] = malloc(k*sizeof(double));
		for(j=0;j<k;j++) forceField.coulombTerms[i][j] = f->coulombTerms[i][j];
	}


	k = forceField.numberOfHydrogenBonded;
	if(k>0)
	for(i=0;i<HYDROGENBONDEDDIM;i++)
	{
		forceField.hydrogenBondedTerms[i] = malloc(k*sizeof(double));
		for(j=0;j<k;j++) forceField.hydrogenBondedTerms[i][j] = f->hydrogenBondedTerms[i][j];
	}

	k = forceField.numberOfSuttonChen;
	if(k>0)
	for(i=0;i<SUTTONCHENDIM;i++)
	{
		forceField.suttonChenTerms[i] = malloc(k*sizeof(double));
		for(j=0;j<k;j++) forceField.suttonChenTerms[i][j] = f->suttonChenTerms[i][j];
	}

	k = forceField.numberOfPairWise = f->numberOfPairWise;
	if(k>0)
	for(i=0;i<PAIRWISEDIM;i++)
	{
		forceField.pairWiseTerms[i] = malloc(k*sizeof(double));
		for(j=0;j<k;j++) forceField.pairWiseTerms[i][j] = f->pairWiseTerms[i][j];
	}
	forceField.logfile = f->logfile;

	forceField.options.type = f->options.type;
	forceField.options.coulomb = f->options.coulomb;
	forceField.options.numeric = f->options.numeric;
	forceField.options.hydrogenBonded612 = f->options.hydrogenBonded612;
	forceField.options.hydrogenBonded1012 = f->options.hydrogenBonded1012;
	forceField.options.hydrogenBondedMorse = f->options.hydrogenBondedMorse;
	forceField.options.improperTorsion = f->options.improperTorsion;
	forceField.options.outOfPlane = f->options.outOfPlane;
	forceField.options.vdw612 = f->options.vdw612;
	forceField.options.vdw714 = f->options.vdw714;
	forceField.options.hbDirectional = f->options.hbDirectional;
	forceField.options.bondStretch = f->options.bondStretch;
	forceField.options.angleBend = f->options.angleBend;
	forceField.options.dihedralAngle = f->options.dihedralAngle;
	forceField.options.rattleConstraints = f->options.rattleConstraints;
	forceField.options.dx = f->options.dx;
	sprintf(forceField.options.chargesType,"%s",f->options.chargesType);

	return forceField;

}
/*****************************************************************************/
void setForceFieldOptions(FILE* file, ForceFieldOptions* forceFieldOptions)
{
	char* tmp = NULL;
	int itmp;
	//forceFieldOptions->type = PAIRWISE;
	forceFieldOptions->type = AMBER;
	forceFieldOptions->bondStretch = TRUE;
	forceFieldOptions->angleBend = TRUE;
	forceFieldOptions->dihedralAngle = TRUE;
	forceFieldOptions->improperTorsion = FALSE;
	forceFieldOptions->strBend = FALSE;
	forceFieldOptions->outOfPlane = FALSE;
	forceFieldOptions->vdw612 = TRUE;
	forceFieldOptions->vdw714 = FALSE;
	forceFieldOptions->hbDirectional = FALSE;
	forceFieldOptions->hydrogenBonded612 = FALSE;
	forceFieldOptions->hydrogenBonded1012 = FALSE;
	forceFieldOptions->hydrogenBondedMorse = FALSE;
	forceFieldOptions->coulomb = TRUE;
	forceFieldOptions->numeric = FALSE;
	forceFieldOptions->dx = 1e-5;
	sprintf(forceFieldOptions->chargesType,"DEFAULT");

	if(readOneInt(file,"ForceFieldType",&itmp)) forceFieldOptions->type = itmp;
	readOneBoolean(file,"ForceFieldUseBond",&forceFieldOptions->bondStretch);
	readOneBoolean(file,"ForceFieldUseBend",&forceFieldOptions->angleBend);
	readOneBoolean(file,"ForceFieldUseStrBend",&forceFieldOptions->strBend);
	readOneBoolean(file,"ForceFieldUseDihedral",&forceFieldOptions->dihedralAngle);
	readOneBoolean(file,"ForceFieldUseImproper",&forceFieldOptions->improperTorsion);
	readOneBoolean(file,"ForceFieldUseOutOfPlane",&forceFieldOptions->outOfPlane);
	readOneBoolean(file,"ForceFieldUseSuttonChen",&forceFieldOptions->suttonChen);
	readOneBoolean(file,"ForceFieldUseCoulomb",&forceFieldOptions->coulomb);
	readOneBoolean(file,"ForceFieldUseNumericalGradient",&forceFieldOptions->numeric);
	readOneString(file,"ForceFieldVanderWals",&tmp);
	if(tmp) uppercase(tmp);
	if(tmp && strstr(tmp,"NONE")) forceFieldOptions->vdw612=FALSE;
	if(tmp && strstr(tmp,"6-12")) { forceFieldOptions->vdw612=TRUE; forceFieldOptions->vdw714=FALSE;}
	if(tmp && strstr(tmp,"7-14")) { forceFieldOptions->vdw612=FALSE; forceFieldOptions->vdw714=TRUE;}
	if(tmp) free(tmp);

	tmp = NULL;
	readOneString(file,"Wall",&tmp);
	if(tmp)  forceFieldOptions->addWallCorrection = TRUE;
	if(tmp) free(tmp);
	tmp = NULL;

	readOneString(file,"ForceFieldHydrogenBonded",&tmp);
	if(tmp) uppercase(tmp);
	if(tmp && strstr(tmp,"NONE")) 
	{
		forceFieldOptions->hydrogenBonded612=FALSE; 
		forceFieldOptions->hydrogenBonded1012=FALSE; 
		forceFieldOptions->hydrogenBondedMorse=FALSE; 
		forceFieldOptions->hbDirectional = FALSE;
	}
	if(tmp && strstr(tmp,"10-12"))
	{
		forceFieldOptions->hydrogenBonded612=FALSE; 
		forceFieldOptions->hydrogenBonded1012=TRUE; 
		forceFieldOptions->hydrogenBondedMorse=FALSE; 
		forceFieldOptions->hbDirectional = FALSE;
	}
	if(tmp && strstr(tmp,"ANGLE-10-12"))
	{
		forceFieldOptions->hydrogenBonded612=FALSE;
		forceFieldOptions->hydrogenBonded1012=TRUE;
		forceFieldOptions->hydrogenBondedMorse=FALSE; 
		forceFieldOptions->hbDirectional = TRUE;
	}
	if(tmp && strstr(tmp,"6-12"))
	{
		forceFieldOptions->hydrogenBonded612=TRUE;
		forceFieldOptions->hydrogenBonded1012=FALSE;
		forceFieldOptions->hydrogenBondedMorse=FALSE; 
		forceFieldOptions->hbDirectional = FALSE;
	}
	if(tmp && strstr(tmp,"ANGLE-6-12"))
	{
		forceFieldOptions->hydrogenBonded612=TRUE;
		forceFieldOptions->hydrogenBonded1012=FALSE;
		forceFieldOptions->hydrogenBondedMorse=FALSE; 
		forceFieldOptions->hbDirectional = TRUE;
	}
	if(tmp && strstr(tmp,"MORSE"))
	{
		forceFieldOptions->hydrogenBonded612=FALSE;
		forceFieldOptions->hydrogenBonded1012=FALSE;
		forceFieldOptions->hydrogenBondedMorse=TRUE; 
		forceFieldOptions->hbDirectional = FALSE;
	}
	if(tmp && strstr(tmp,"ANGLE-MORSE"))
	{
		forceFieldOptions->hydrogenBonded612=FALSE;
		forceFieldOptions->hydrogenBonded1012=FALSE;
		forceFieldOptions->hydrogenBondedMorse=TRUE; 
		forceFieldOptions->hbDirectional = TRUE;
	}
	if(tmp) free(tmp);
	tmp = NULL;

	readOneString(file,"ForceFieldChargesType",&tmp);
	if(tmp) uppercase(tmp);
	if(tmp && strstr(tmp,"MM")) sprintf(forceFieldOptions->chargesType,"MM");
	if(tmp && strstr(tmp,"SCALED")) sprintf(forceFieldOptions->chargesType,"SCALED");
	if(tmp && strstr(tmp,"EEM")) sprintf(forceFieldOptions->chargesType,"EEM");
	if(tmp && strstr(tmp,"ACKH2")) sprintf(forceFieldOptions->chargesType,"ACKS2");
	if(tmp && strstr(tmp,"ACKS2")) sprintf(forceFieldOptions->chargesType,"ACKS2");
	if(tmp && strstr(tmp,"ACKS2-BEGIN")) sprintf(forceFieldOptions->chargesType,"ACKS2-BEGIN");
	if(tmp && strstr(tmp,"EEM-BEGIN")) sprintf(forceFieldOptions->chargesType,"EEM-BEGIN");
	if(tmp) free(tmp);
	tmp = NULL;
	

	forceFieldOptions->rattleConstraints = NOCONSTRAINTS;
	//forceFieldOptions->rattleConstraints = BONDSCONSTRAINTS;
	//forceFieldOptions->rattleConstraints = BONDSANGLESCONSTRAINTS;
	if(readOneInt(file,"Constraints",&itmp))forceFieldOptions->rattleConstraints = itmp;

	readOneReal(file,"dx",&forceFieldOptions->dx);

}
/**********************************************************************/
void updateGeometryVelocitiesCL(ForceField* forceField, Molecule* mol)
{
#ifdef ENABLE_CL
	int i;
	CLProp clProp = getCLProp();
	if(!mol) mol = &forceField->molecule;
	for(i=0;i<mol->nAtoms;i++)
	{
		forceField->atomsCPU[i].s[0] = mol->atoms[i].coordinates[0];
		forceField->atomsCPU[i].s[1] = mol->atoms[i].coordinates[1];
		forceField->atomsCPU[i].s[2] = mol->atoms[i].coordinates[2];
		forceField->atomsCPU[i].s[3] = mol->atoms[i].charge;
		forceField->atomsCPU[i].s[4] = mol->atoms[i].mass;
		//printf("Mass = %f\n",mol->atoms[i].mass);
		// velocity
		{
			forceField->atomsCPU[i].s[5] = mol->atoms[i].velocity[0];
			forceField->atomsCPU[i].s[6] = mol->atoms[i].velocity[1];
			forceField->atomsCPU[i].s[7] = mol->atoms[i].velocity[2];
		}
	}
	clEnqueueWriteBuffer(clProp.command_queue, forceField->atomsCL, CL_TRUE, 0, sizeof(cl_float8) * mol->nAtoms, forceField->atomsCPU, 0, NULL, NULL);
#endif
}
/**********************************************************************/
void updateGeometryCL(ForceField* forceField, Molecule* mol)
{
	updateGeometryVelocitiesCL(forceField, mol);
}
/**********************************************************************/
void getGeometryVelocitiesCL(ForceField* forceField, Molecule* mol)
{
#ifdef ENABLE_CL
	int i,j;
	CLProp clProp = getCLProp();
	if(!mol) mol = &forceField->molecule;
	clEnqueueReadBuffer(clProp.command_queue, forceField->atomsCL, CL_TRUE, 0, sizeof(cl_float8)* mol->nAtoms, forceField->atomsCPU, 0, NULL, NULL);
	for(i=0;i<mol->nAtoms;i++)
	{
		for(j=0;j<3;j++) mol->atoms[i].coordinates[j] = forceField->atomsCPU[i].s[j];
		// velocity
		for(j=0;j<3;j++) mol->atoms[i].velocity[j] = forceField->atomsCPU[i].s[j+5];
	}
#endif
}
/*****************************************************************************/
/* in Debye */
static void calculDipole(Molecule* mol, double D[])
{
	int i,k;
	for(k=0;k<3;k++) D[k] = 0;
	for(i=0;i<mol->nAtoms;i++)
	for(k=0;k<3;k++)
		D[k] += mol->atoms[i].charge*mol->atoms[i].coordinates[k];
	for(k=0;k<3;k++) D[k] *= ANGTOBOHR*AUTODEB;
}
/*****************************************************************************/
static void copyGradients(Molecule* mol, double* g[])
{
	int i,k;
	if(!mol) return;
	for(i=0;i<mol->nAtoms;i++)
		for(k=0;k<3;k++)
			g[k][i] = mol->atoms[i].gradient[k];
}
/*****************************************************************************/
static void sortFrequencies(int nModes, double* frequencies, double** modes, double* reducedMasses, double* IRIntensities)
{
	int i;
	int j;
	int k;
	double dum;
	if(nModes<1 || !frequencies || !modes || !reducedMasses || !IRIntensities) return;
	for(i=0;i<nModes;i++)
	{
		k = i;
		for(j=i+1;j<nModes;j++)
			if(frequencies[j]<frequencies[k]) k = j;
		if(k==i) continue;
		/* swap i and k modes */
		dum = frequencies[i];
		frequencies[i] = frequencies[k];
		frequencies[k] = dum;
		dum = reducedMasses[i];
		reducedMasses[i] = reducedMasses[k];
		reducedMasses[k] = dum;
		dum = IRIntensities[i];
		IRIntensities[i] = IRIntensities[k];
		IRIntensities[k] = dum;
		for(j=0;j<nModes;j++)
		{
			dum =  modes[j][i];
			modes[j][i] = modes[j][k];
			modes[j][k] = dum;
		}
	}
}

/*****************************************************************************/
int computeMMFrequencies(ForceField* forceField, double** frequencies, double*** modes, double** reducedMasses, double** IRIntensities)
{
	int i;
	int j;
	int k;
	int c;
	int id,jd,index;
	double* F;
	double* gp[3];
	double* gm[3];
	double* dmuX[3];
	Molecule* mol;
	int nAtoms;
	double Dp[3];
	double Dm[3];
	double dx;

	if(!forceField || forceField->molecule.nAtoms<1) return 0;
	dx = forceField->options.dx;

	//printf("Begin calcul F\n");

	mol = &forceField->molecule;
	nAtoms = mol->nAtoms;
	for(k=0;k<3;k++) gp[k] = malloc(nAtoms*sizeof(double));
	for(k=0;k<3;k++) gm[k] = malloc(nAtoms*sizeof(double));
	for(k=0;k<3;k++) dmuX[k] = malloc(3*nAtoms*sizeof(double));

	F = malloc(3*nAtoms*(3*nAtoms+1)/2*sizeof(double));

	//printf("End alloc calcul F\n");

	index = 0;
	for(i=0;i<nAtoms;i++)
	for(k=0;k<3;k++)
	{
		id=3*i+k;
		mol->atoms[i].coordinates[k] += dx;
		forceField->klass->calculateGradient(forceField);
		copyGradients(mol, gp);
		calculDipole(mol, Dp);

		mol->atoms[i].coordinates[k] -= 2*dx;
		forceField->klass->calculateGradient(forceField);
		copyGradients(mol, gm);
		calculDipole(mol, Dm);

		for(c = 0;c<3;c++) dmuX[c][id] = (Dp[c]-Dm[c])/dx/2;
		mol->atoms[i].coordinates[k] += dx;
		for(j=0;j<=i;j++)
		{
			double invm = 1.0/sqrt( mol->atoms[i].mass* mol->atoms[j].mass);
			for(c = 0;c<3;c++) 
			{
				jd = 3*j+c;
				if(jd>id) continue;
				index = jd + id*(id+1)/2;
				F[index] = (gp[c][j]-gm[c][j])/dx/2; 
				F[index] *= invm;
			}
		}
	}
	//printf("En calcul F\n");
	for(k=0;k<3;k++) free(gp[k]);
	for(k=0;k<3;k++) free(gm[k]);
	*frequencies = malloc(3*nAtoms*sizeof(double));
	*reducedMasses = malloc(3*nAtoms*sizeof(double));
	*IRIntensities = malloc(3*nAtoms*sizeof(double));
	*modes = malloc(3*nAtoms*sizeof(double*));
	for(i=0;i<3*nAtoms;i++) (*modes)[i] = malloc(3*nAtoms*sizeof(double));

	eigenQL(3*nAtoms, F, *frequencies, *modes);
	free(F);
	/* convert in atomic unit  from kcal/Ang^2/amu */
	for(i=0;i<3*nAtoms;i++) (*frequencies)[i] *= 1.59360150e-03*0.529177*0.529177*5.48579911e-04; 
	/* convert frequencies in cm-1 */
	for(i=0;i<3*nAtoms;i++) 
		if( (*frequencies)[i]>0) (*frequencies)[i] = sqrt((*frequencies)[i])*219474.63633664;
		else (*frequencies)[i] = -sqrt(-(*frequencies)[i])*219474.63633664;

	/* compute the IR intensities */
	for(i=0;i<nAtoms;i++)
	for(k=0;k<3;k++)
	{
		int id=3*i+k;
		double IRI = 0;
		double D[3] = {0,0,0};
		int kp;
		for(c = 0;c<3;c++)
		for(j=0;j<nAtoms;j++)
		for(kp = 0;kp<3;kp++) 
		{
			int jd = 3*j+kp;
			double Lji = (*modes)[jd][id];
			double a=dmuX[c][jd]*Lji/sqrt(mol->atoms[j].mass);
			D[c]+=a;
		}
		IRI = 0;
		for(c = 0;c<3;c++)  IRI+= D[c]*D[c];
		(*IRIntensities)[id] = IRI;
	}
	for(k=0;k<3;k++) free(dmuX[k]);
	/* Intensities in 1 (D/Ang)^2 amu^-1 = 42.255 km/mol=171.65 cm^-2 atm^-1 at 0 C and 1 atm */
	/* Refs : D. Porezag and M. R. Pederson, Phys. Rev. B 54, 7830 (1996). and Y. Yamaguchi el al., J. Chem. Phys. 84,2262(1986)*/
	/* conversion in km/mol*/
	for(i=0;i<3*nAtoms;i++) (*IRIntensities)[i] *= 42.255;

	/* compute the reduced mass */
	for(i=0;i<3*nAtoms;i++) 
	{
		double m = 0;
		for(j=0;j<mol->nAtoms;j++)
		{
			double r2 = 0;
			for(c=0;c<3;c++) r2+= (*modes)[3*j+c][i]*(*modes)[3*j+c][i];
			m+= r2/(mol->atoms[j].mass); 
		}
		if(m<=0) m = 1;
		m = 1/m;
		for(j=0;j<mol->nAtoms;j++)
		{
			double r =sqrt(m)/sqrt(mol->atoms[j].mass);
			for(c=0;c<3;c++) (*modes)[3*j+c][i]*=r;
		}

		(*reducedMasses)[i] = m;
	}
	sortFrequencies(3*nAtoms, *frequencies, *modes, *reducedMasses, *IRIntensities);

	return 3*nAtoms;

}
/*****************************************************************************/
static void computeGradientNumeric(ForceField* forceField)
{
	int i;
	int k;
	Molecule* mol;
	int nAtoms;
	double Ep, Em;
	double dx;

	if(!forceField || forceField->molecule.nAtoms<1) return;
	dx = forceField->options.dx;

	mol = &forceField->molecule;
	nAtoms = mol->nAtoms;
	for(i=0;i<nAtoms;i++)
	for(k=0;k<3;k++)
		mol->atoms[i].gradient[k] = 0.0;

	for(i=0;i<nAtoms;i++)
	for(k=0;k<3;k++)
	{
		mol->atoms[i].coordinates[k] += dx;
		forceField->klass->calculateEnergy(forceField);
		Ep = mol->potentialEnergy;
		

		mol->atoms[i].coordinates[k] -= 2*dx;
		forceField->klass->calculateEnergy(forceField);
		Em = mol->potentialEnergy;

		mol->atoms[i].gradient[k] = (Ep-Em)/dx/2;
		mol->atoms[i].coordinates[k] += dx;
	}
	forceField->klass->calculateEnergy(forceField);
	//printf("End computeGradientNumeric\n");
}
/*****************************************************************************/
