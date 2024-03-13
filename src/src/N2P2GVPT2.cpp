/*****************************************************************************************
 iGVPT2 is a program for computing anharmonic corrections to vibration frequencies, 
 based on force field expansion of the potential energy surface in normal mode coordinates.
 iGVPT2 supports several computation chemistry packages(CCP).

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
*****************************************************************************************/

/**********************************************************************
N2P2GVPT2.cpp - calculate GVPT2
***********************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <time.h>
#include <locale.h>

#ifdef __cplusplus
extern "C" {
#include <Utils/Types.h>
#include <Utils/Utils.h>
#include <JobControl/Job.h>
#include <QuarticForceField/QFFnMR.h>
#include <MolecularMechanics/MolecularMechanics.h>

#include <QuantumMechanics/QuantumMechanics.h>
#ifdef HIGH_DERIVATIVES
#include <QuarticForceField/QFFN2P2.h>
#endif

#include <VPT2/VPT2.h>
#include <Utils/Constants.h>
}
#endif


#include <iomanip>
#include <iostream>
#ifndef _MSC_VER
  #include <unistd.h>
#endif
#include <cstring>

using namespace std;

static bool setACKS2Parameters(char* inputFileName, Molecule* mol);
static QuantumMechanicsModel n2p2(char*  filename);
static bool resetCharges(char* inputFileName, QuantumMechanicsModel& n2p2Model);
static void computeDipole(Molecule* mol, double** geom, double* dipole);

/**********************************************************************************************************************/
static int initMol(char* inputFileName, int* pOrdre, double** pDeltas, QuantumMechanicsModel& n2p2Model, bool* pACKS2Charges)
{
	Molecule* cchemiMol = &n2p2Model.molecule;
	double delta = 0.5;
	boolean reducedCoordinates = TRUE;
	int ordre = 2;
	int nAll;
	FILE* file = fopen(inputFileName,"rb");
	double* deltas = NULL;
	bool ACKS2Charges = FALSE;

	readOneBoolean(file,(char*)"QFFReducedCoordinates",&reducedCoordinates);
	readOneInt(file,(char*)"QFFnModes",&ordre);
	if(!reducedCoordinates) delta = 1e-2;
	readOneReal(file,(char*)"QFFDelta",&delta);
	fclose(file);
	printf("delta = %f ",delta);
	if(reducedCoordinates) printf("reducedCoordiantes\n");
	else printf(" Angshtrom\n");
	printf("nAtoms = %d\n",cchemiMol->nAtoms);
	printf("nModes = %d\n",cchemiMol->vibration.nModes);
	deltas = cchemiMol->klass->getDeltaTable(cchemiMol, delta, reducedCoordinates);


	if(cchemiMol) ACKS2Charges = setACKS2Parameters(inputFileName, cchemiMol);

	if(resetCharges(inputFileName, n2p2Model)) ACKS2Charges = FALSE;

	*pOrdre = ordre;
	*pDeltas = deltas;
	*pACKS2Charges = ACKS2Charges;
	//mol->klass->free(mol);
	nAll = 1 +6* cchemiMol->vibration.nModes;
	if(ordre>=2) nAll += 6*cchemiMol->vibration.nModes*(cchemiMol->vibration.nModes-1); 
	if(ordre>=3) nAll += 8*cchemiMol->vibration.nModes*(cchemiMol->vibration.nModes-1)*(cchemiMol->vibration.nModes-2)/6;
	if(ordre>=4) nAll += 16*cchemiMol->vibration.nModes*(cchemiMol->vibration.nModes-1)*(cchemiMol->vibration.nModes-2)*(cchemiMol->vibration.nModes-3)/24;
	return nAll;
}
/**********************************************************************************************************************/
static void computeEnergyAndDipole(QuantumMechanicsModel& n2p2Model, bool ACKS2Charges, double** geom, double* energies, double* dipoles[3], int k)
{
	Molecule* cchemiMol = &n2p2Model.molecule;
	double dipole[3];
	{
		for(int i=0;i<cchemiMol->nAtoms;i++)
		{
			for(int xyz=0;xyz<3;xyz++)
				cchemiMol->atoms[i].coordinates[xyz] = geom[i][xyz];
		}
		n2p2Model.klass->calculateEnergy(&n2p2Model);
		energies[k] = n2p2Model.molecule.potentialEnergy;
	//	energies[k] /= KCALTOKJ;
		if(ACKS2Charges)
		{
			//cerr<<" I reset the ACKS2 charges "<<endl;
			cchemiMol->klass->setChargesACKS2(cchemiMol);
		}
                //cerr<<"charge0"<<cchemiMol->atoms[0].charge<<endl;
		computeDipole(cchemiMol, geom, dipole);
		energies[k] /= AUTOKCAL;
		for(int xyz=0;xyz<3;xyz++) dipoles[xyz][k] = dipole[xyz];
		for(int xyz=0;xyz<3;xyz++)  dipoles[xyz][k] /=  AUTODEB;

		/*
		cout<<"Energy = "<<energies[k]<<endl;
		FOR_ATOMS_OF_MOL(atom, n2p2Model) cerr << " x: " << atom->x() << " y: " << atom->y() << " z: " << atom->z() << endl;
		*/
	}
}
/**********************************************************************************************************************/
static void getEnergiesAndDipoles(double* deltas, QuantumMechanicsModel& n2p2Model, bool ACKS2Charges, int ordre, double* energies, double* dipoles[3])
{
	double** geom;
	Molecule* mol = NULL;
	Molecule molecule = *(n2p2Model.molecule.klass->copy(&n2p2Model.molecule));
	int i,j,k,l;
	int index;
	mol = &molecule;

	geom =  new double* [mol->nAtoms];
        for(i=0;i<mol->nAtoms;i++) geom[i] =  new double[3];

	index = 0;
	mol->klass->getQFFOneGeom(mol, -1, -1, 0.0, 0.0, 0.0,geom);
	computeEnergyAndDipole(n2p2Model, ACKS2Charges, geom, energies, dipoles, index);
	for(j=0;j<mol->vibration.nModes;j++)
	{
		index++; mol->klass->getQFFOneGeom(mol,  j, -1,   3*deltas[j], 0,0,geom);
		computeEnergyAndDipole( n2p2Model, ACKS2Charges, geom, energies, dipoles, index);

		index++; mol->klass->getQFFOneGeom(mol,  j, -1,   2*deltas[j], 0,0,geom);
		computeEnergyAndDipole( n2p2Model, ACKS2Charges, geom, energies, dipoles, index);

		index++; mol->klass->getQFFOneGeom(mol,  j, -1,     deltas[j], 0,0,geom);
		computeEnergyAndDipole( n2p2Model, ACKS2Charges, geom, energies, dipoles, index);

		index++; mol->klass->getQFFOneGeom(mol,  j, -1,    -deltas[j], 0,0,geom);
		computeEnergyAndDipole( n2p2Model, ACKS2Charges, geom, energies, dipoles, index);

		index++; mol->klass->getQFFOneGeom(mol,  j, -1,  -2*deltas[j], 0,0,geom);
		computeEnergyAndDipole( n2p2Model, ACKS2Charges, geom, energies, dipoles, index);

		index++; mol->klass->getQFFOneGeom(mol,  j, -1,  -3*deltas[j], 0,0,geom);
		computeEnergyAndDipole( n2p2Model, ACKS2Charges, geom, energies, dipoles, index);
	}
	if(ordre>=2)
	for(j=0;j<mol->vibration.nModes;j++)
	{
		for(i=0;i<j;i++)
		{
			index++; mol->klass->getQFFOneGeom(mol, j, i,  deltas[j],   deltas[i],0, geom);
			computeEnergyAndDipole( n2p2Model, ACKS2Charges, geom, energies, dipoles, index);

			index++; mol->klass->getQFFOneGeom(mol, j, i,  deltas[j],  -deltas[i],0, geom);
			computeEnergyAndDipole( n2p2Model, ACKS2Charges, geom, energies, dipoles, index);

			index++; mol->klass->getQFFOneGeom(mol, j, i, -deltas[j],   deltas[i],0, geom);
			computeEnergyAndDipole( n2p2Model, ACKS2Charges, geom, energies, dipoles, index);

			index++; mol->klass->getQFFOneGeom(mol, j, i, -deltas[j],  -deltas[i],0, geom);
			computeEnergyAndDipole( n2p2Model, ACKS2Charges, geom, energies, dipoles, index);
		}
		for(i=0;i<mol->vibration.nModes;i++)
		{
			if(i==j) continue;

			index++; mol->klass->getQFFOneGeom(mol,  j, i,  deltas[j],  3*deltas[i],0, geom);
			computeEnergyAndDipole( n2p2Model, ACKS2Charges, geom, energies, dipoles, index);

			index++; mol->klass->getQFFOneGeom(mol,  j, i,  deltas[j], -3*deltas[i],0, geom);
			computeEnergyAndDipole( n2p2Model, ACKS2Charges, geom, energies, dipoles, index);

			index++; mol->klass->getQFFOneGeom(mol,  j, i, -deltas[j],  3*deltas[i],0, geom);
			computeEnergyAndDipole( n2p2Model, ACKS2Charges, geom, energies, dipoles, index);

			index++; mol->klass->getQFFOneGeom(mol,  j, i, -deltas[j], -3*deltas[i],0, geom);
			computeEnergyAndDipole( n2p2Model, ACKS2Charges, geom, energies, dipoles, index);
		}
	}
	if(ordre>=3)
	for(j=0;j<mol->vibration.nModes;j++)
	{
		for(i=0;i<j;i++)
		for(k=0;k<i;k++)
		{
			/* 0  deltas[j],   deltas[i],  deltas[k] */
			index++; mol->klass->getQFFOneGeom3(mol, j, i, k,  deltas[j], deltas[i], deltas[k], geom);
			computeEnergyAndDipole( n2p2Model, ACKS2Charges, geom, energies, dipoles, index);
			/* 1  deltas[j],   deltas[i], -deltas[k] */
			index++; mol->klass->getQFFOneGeom3(mol, j, i, k,  deltas[j], deltas[i], -deltas[k], geom);
			computeEnergyAndDipole( n2p2Model, ACKS2Charges, geom, energies, dipoles, index);
			/* 2  deltas[j],  -deltas[i],  deltas[k] */
			index++; mol->klass->getQFFOneGeom3(mol, j, i, k,  deltas[j], -deltas[i], deltas[k], geom);
			computeEnergyAndDipole( n2p2Model, ACKS2Charges, geom, energies, dipoles, index);
			/* 3  deltas[j],  -deltas[i], -deltas[k] */
			index++; mol->klass->getQFFOneGeom3(mol, j, i, k,  deltas[j], -deltas[i], -deltas[k], geom);
			computeEnergyAndDipole( n2p2Model, ACKS2Charges, geom, energies, dipoles, index);
			/* 4 -deltas[j],  deltas[i],  deltas[k] */
			index++; mol->klass->getQFFOneGeom3(mol, j, i, k,  -deltas[j], deltas[i], deltas[k], geom);
			computeEnergyAndDipole( n2p2Model, ACKS2Charges, geom, energies, dipoles, index);
			/* 5 -deltas[j],  deltas[i], -deltas[k] */
			index++; mol->klass->getQFFOneGeom3(mol, j, i, k,  -deltas[j], deltas[i], -deltas[k], geom);
			computeEnergyAndDipole( n2p2Model, ACKS2Charges, geom, energies, dipoles, index);
			/* 6 -deltas[j], -deltas[i],  deltas[k] */
			index++; mol->klass->getQFFOneGeom3(mol, j, i, k,  -deltas[j], -deltas[i], deltas[k], geom);
			computeEnergyAndDipole( n2p2Model, ACKS2Charges, geom, energies, dipoles, index);
			/* 7 -deltas[j], -deltas[i], -deltas[k]- */
			index++; mol->klass->getQFFOneGeom3(mol, j, i, k,  -deltas[j], -deltas[i], -deltas[k], geom);
			computeEnergyAndDipole( n2p2Model, ACKS2Charges, geom, energies, dipoles, index);
		}
	}
	if(ordre>=4)
	for(j=0;j<mol->vibration.nModes;j++)
	{
		for(i=0;i<j;i++)
		for(k=0;k<i;k++)
		for(l=0;l<k;l++)
		{
				/* VIJKL */
				/* 0    deltas[j],   deltas[i],   deltas[k],  deltas[l] */ //+
				/* 1    deltas[j],   deltas[i],   deltas[k], -deltas[l] */ //-
				/* 2    deltas[j],   deltas[i],  -deltas[k],  deltas[l] */ //-
				/* 3    deltas[j],   deltas[i],  -deltas[k], -deltas[l] */ //+
				/* 4    deltas[j],  -deltas[i],   deltas[k],  deltas[l] */ //-
				/* 5    deltas[j],  -deltas[i],   deltas[k], -deltas[l] */ //+
				/* 6    deltas[j],  -deltas[i],  -deltas[k],  deltas[l] */ //+
				/* 7    deltas[j],  -deltas[i],  -deltas[k], -deltas[l] */ //-
				/* 8   -deltas[j],   deltas[i],   deltas[k],  deltas[l] */ //-
				/* 9   -deltas[j],   deltas[i],   deltas[k], -deltas[l] */ //+
				/* 10  -deltas[j],   deltas[i],  -deltas[k],  deltas[l] */ //+
				/* 11  -deltas[j],   deltas[i],  -deltas[k], -deltas[l] */ //-
				/* 12  -deltas[j],  -deltas[i],   deltas[k],  deltas[l] */ //+
				/* 13  -deltas[j],  -deltas[i],   deltas[k], -deltas[l] */ //-
				/* 14  -deltas[j],  -deltas[i],  -deltas[k],  deltas[l] */ //-
				/* 15  -deltas[j],  -deltas[i],  -deltas[k], -deltas[l] */ //+

			index++; mol->klass->getQFFOneGeom4(mol, j, i, k, l,  deltas[j], deltas[i], deltas[k], deltas[l], geom);
			computeEnergyAndDipole( n2p2Model, ACKS2Charges, geom, energies, dipoles, index);
			index++; mol->klass->getQFFOneGeom4(mol, j, i, k, l,  deltas[j], deltas[i], deltas[k], -deltas[l], geom);
			computeEnergyAndDipole( n2p2Model, ACKS2Charges, geom, energies, dipoles, index);
			index++; mol->klass->getQFFOneGeom4(mol, j, i, k, l,  deltas[j], deltas[i], -deltas[k], deltas[l], geom);
			computeEnergyAndDipole( n2p2Model, ACKS2Charges, geom, energies, dipoles, index);
			index++; mol->klass->getQFFOneGeom4(mol, j, i, k, l,  deltas[j], deltas[i], -deltas[k], -deltas[l], geom);
			computeEnergyAndDipole( n2p2Model, ACKS2Charges, geom, energies, dipoles, index);
			index++; mol->klass->getQFFOneGeom4(mol, j, i, k, l,  deltas[j], -deltas[i], deltas[k], deltas[l], geom);
			computeEnergyAndDipole( n2p2Model, ACKS2Charges, geom, energies, dipoles, index);
			index++; mol->klass->getQFFOneGeom4(mol, j, i, k, l,  deltas[j], -deltas[i], deltas[k], -deltas[l], geom);
			computeEnergyAndDipole( n2p2Model, ACKS2Charges, geom, energies, dipoles, index);
			index++; mol->klass->getQFFOneGeom4(mol, j, i, k, l,  deltas[j], -deltas[i], -deltas[k], deltas[l], geom);
			computeEnergyAndDipole( n2p2Model, ACKS2Charges, geom, energies, dipoles, index);
			index++; mol->klass->getQFFOneGeom4(mol, j, i, k, l,  deltas[j], -deltas[i], -deltas[k], -deltas[l], geom);
			computeEnergyAndDipole( n2p2Model, ACKS2Charges, geom, energies, dipoles, index);

			index++; mol->klass->getQFFOneGeom4(mol, j, i, k, l, -deltas[j], deltas[i], deltas[k], deltas[l], geom);
			computeEnergyAndDipole( n2p2Model, ACKS2Charges, geom, energies, dipoles, index);
			index++; mol->klass->getQFFOneGeom4(mol, j, i, k, l, -deltas[j], deltas[i], deltas[k], -deltas[l], geom);
			computeEnergyAndDipole( n2p2Model, ACKS2Charges, geom, energies, dipoles, index);
			index++; mol->klass->getQFFOneGeom4(mol, j, i, k, l, -deltas[j], deltas[i], -deltas[k], deltas[l], geom);
			computeEnergyAndDipole( n2p2Model, ACKS2Charges, geom, energies, dipoles, index);
			index++; mol->klass->getQFFOneGeom4(mol, j, i, k, l, -deltas[j], deltas[i], -deltas[k], -deltas[l], geom);
			computeEnergyAndDipole( n2p2Model, ACKS2Charges, geom, energies, dipoles, index);
			index++; mol->klass->getQFFOneGeom4(mol, j, i, k, l, -deltas[j], -deltas[i], deltas[k], deltas[l], geom);
			computeEnergyAndDipole( n2p2Model, ACKS2Charges, geom, energies, dipoles, index);
			index++; mol->klass->getQFFOneGeom4(mol, j, i, k, l, -deltas[j], -deltas[i], deltas[k], -deltas[l], geom);
			computeEnergyAndDipole( n2p2Model, ACKS2Charges, geom, energies, dipoles, index);
			index++; mol->klass->getQFFOneGeom4(mol, j, i, k, l, -deltas[j], -deltas[i], -deltas[k], deltas[l], geom);
			computeEnergyAndDipole( n2p2Model, ACKS2Charges, geom, energies, dipoles, index);
			index++; mol->klass->getQFFOneGeom4(mol, j, i, k, l, -deltas[j], -deltas[i], -deltas[k], -deltas[l], geom);
			computeEnergyAndDipole( n2p2Model, ACKS2Charges, geom, energies, dipoles, index);
		}
	}
	molecule.klass->free(&molecule);
}
/*****************************************************************************/
static bool setACKS2Parameters(char* inputFileName, Molecule* mol)
{
        ForceField forceField;
        ForceFieldOptions forceFieldOptions;
        FILE* file = fopen(inputFileName,"rb");

        setForceFieldOptions(file, &forceFieldOptions);
        mol->klass->buildMMTypes(mol, file);
/* Molecule to read */
        forceField = createAmberModel(mol,forceFieldOptions, stdout);
        fclose(file);
	//cerr<<"forceField.options.chargesType="<<forceField.options.chargesType<<endl;
	if(strstr(forceField.options.chargesType,"ACKS2")) 
	{
		*mol  = *(forceField.molecule.klass->copy(&forceField.molecule));
		mol->klass->setChargesACKS2(mol);
		cerr<<"Warning : I computed the ACKS2 charges"<<endl;
		if(strstr(forceField.options.chargesType,"BEGIN")) return FALSE;
		cerr<<"Warning : I computed these charges at each step"<<endl;
		return TRUE;
	}

	return false;
}
/*****************************************************************************/
static QuantumMechanicsModel n2p2(char*  filename)
{
	QuantumMechanicsModel qmModel;
	Constraints constraints = NOCONSTRAINTS;
	Molecule* mol = readMolecule(filename,TRUE);
	char* n2p2Dir=strdup(getenv("PWD"));
	FILE* file = fopen(filename,"rb");
	readOneString(file,(char*)"HDNNDir",&n2p2Dir);
	readOneString(file,(char*)"n2p2Dir",&n2p2Dir);
	fclose(file);
	qmModel = createN2P2Model(mol, n2p2Dir, constraints, stderr);
	qmModel.klass->calculateEnergy(&qmModel);

	cout << "\nFINAL ENERGY:" << setw(20) << fixed << setprecision(14) << qmModel.molecule.potentialEnergy << endl;
	return qmModel;
}
/****************************************************************************************/
static bool resetCharges(char* inputFileName, QuantumMechanicsModel& n2p2Model)
{
	Molecule* mol = &n2p2Model.molecule;
	int i =0;
	for (i = 0; i<mol->nAtoms;i++) mol->atoms[i].charge = 0;
	return TRUE;
}
/****************************************************************************************/
static void computeDipole(Molecule* mol, double** geom, double* dipole)
{
       int i,k;
        for(k=0;k<3;k++) dipole[k] = 0;
        for(i=0;i<mol->nAtoms;i++)
        for(k=0;k<3;k++)
               dipole[k] += mol->atoms[i].charge*geom[i][k];
        for(k=0;k<3;k++) dipole[k] *= ANGTOBOHR*AUTODEB;
}
/****************************************************************************************/
static int n2p2GVPT2Numeric(char* inputFileName)
{
	Molecule* cchemiMol = NULL;
	char* fileNameDerivatives;
	QuantumMechanicsModel n2p2Model;
	int ordre;
	double* deltas = NULL;
	int nAll = 0;
	double* energies = NULL;
	double* dipoles[3] =  {NULL, NULL, NULL};

        double elaps;
        double elapsAll = 0;
	double time0 = get_cpu_time();
	double time1 = 0;
	bool ACKS2Charges = FALSE;

	n2p2Model = n2p2(inputFileName);
	cchemiMol = &n2p2Model.molecule; 
	nAll = initMol(inputFileName, &ordre, &deltas, n2p2Model,&ACKS2Charges);
	energies = new double [nAll];
	for(int xyz=0;xyz<3;xyz++)   dipoles[xyz] = new double [nAll];

	time1 = get_cpu_time();
	elaps = time1-time0;
	elapsAll += elaps;
	cerr<<"CPU time for initialisation "<<elaps<<" seconds"<<endl;

	time0 = time1;
	getEnergiesAndDipoles(deltas, n2p2Model, ACKS2Charges, ordre, energies, dipoles);
	time1 = get_cpu_time();
	elaps = time1-time0;
	elapsAll += elaps;
	cerr<<"CPU time for computing of energies and dipoles "<<elaps<<" seconds"<<endl;

	time0 = time1;

	n2p2Model = n2p2(inputFileName);
	cchemiMol = &n2p2Model.molecule;

	fileNameDerivatives = computeQFFDerivatives(inputFileName, cchemiMol, ordre, deltas, nAll, energies, dipoles);
	time1 = get_cpu_time();
	elaps = time1-time0;
	elapsAll += elaps;
	cerr<<"CPU time to compute derivatives "<<elaps<<" seconds"<<endl;
	time0 = time1;

	vpt2(fileNameDerivatives);
	time1 = get_cpu_time();
	elaps = time1-time0;
	elapsAll += elaps;
	cerr<<"CPU time to compute frequencies and intensities "<<elaps<<" seconds"<<endl;
	cerr<<endl<<"All time "<<elapsAll<<" seconds"<<endl;

	/* 
	finalize();
	*/

	return 0;
}
/****************************************************************************************/
#ifdef HIGH_DERIVATIVES
static int n2p2GVPT2Analytic(char* inputFileName, int order, int method)
{
	Molecule* cchemiMol = NULL;
	char* fileNameDerivatives;
	QuantumMechanicsModel n2p2Model;

        double elaps;
        double elapsAll = 0;
	double time0 = get_cpu_time();
	double time1 = 0;
	int orderDipole = order;

	n2p2Model = n2p2(inputFileName);
	cchemiMol = &n2p2Model.molecule; 

	time1 = get_cpu_time();
	elaps = time1-time0;
	elapsAll += elaps;
	cerr<<"CPU time for initialisation "<<elaps<<" seconds"<<endl;

	time0 = time1;


	if(order>=4) orderDipole--;
	fileNameDerivatives = computeN2P2QFFAnalyticDerivatives(inputFileName,  cchemiMol, order,  orderDipole, method, n2p2Model.N2P2Dir);
	time1 = get_cpu_time();
	elaps = time1-time0;
	elapsAll += elaps;
	cerr<<"CPU time to compute derivatives "<<elaps<<" seconds"<<endl;
	time0 = time1;

	vpt2(fileNameDerivatives);
	time1 = get_cpu_time();
	elaps = time1-time0;
	elapsAll += elaps;
	cerr<<"CPU time to compute frequencies and intensities "<<elaps<<" seconds"<<endl;
	cerr<<endl<<"All time "<<elapsAll<<" seconds"<<endl;

	/* 
	finalize();
	*/

	return 0;
}
/**********************************************************************************************************************/
static bool getN2P2HighDerivativesMethod(char* inputFileName, int* pOrdre, int*pMethod)
{
	int ordre = 2;
	int method = -1;
	FILE* file = fopen(inputFileName,"rb");

	if(!file) return (method<0);
	readOneInt(file,(char*)"QFFnModes",&ordre);
	readOneInt(file,(char*)"HighDerivativesMethod",&method);
	fclose(file);
	*pOrdre = ordre; 
	*pMethod = method;
	if(method>=0) *pOrdre = ordre+1; // order = nMR order : 3=> 3 index => QFF (derivatives of energy to order 4)
	return (method<0);
}
#endif
/**********************************************************************************************************************/
int n2p2GVPT2(char* inputFileName)
{
	int method;
	int order = 3;
	getN2P2HighDerivativesMethod(inputFileName, &order, &method);
	if(method<0) return n2p2GVPT2Numeric(inputFileName);
#ifdef HIGH_DERIVATIVES
	else return n2p2GVPT2Analytic(inputFileName, order, method);
#else
	else return 1;
#endif
	
}
