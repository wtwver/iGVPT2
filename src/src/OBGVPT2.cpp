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
obgvpt2.cpp - calculate GVPT2
***********************************************************************/

#ifdef WIN32
#define USING_OBDLL
#endif

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
#include <VPT2/VPT2.h>
#include <Utils/Constants.h>
}
#endif


#include <openbabel/babelconfig.h>
#include <openbabel/base.h>
#include <openbabel/mol.h>
#include <openbabel/obconversion.h>
#include <openbabel/forcefield.h>
#include <openbabel/chargemodel.h>
#include <iomanip>
#ifndef _MSC_VER
  #include <unistd.h>
#endif
#include <cstring>

using namespace std;
using namespace OpenBabel;

static bool setACKS2Parameters(char* inputFileName, Molecule* mol);
static OBMol ob(char*  filename, OBForceField*& pFF);
static bool resetCharges(char* inputFileName, OBMol& obMol,  Molecule* mol);
static void computeDipole(Molecule* mol, double** geom, double* dipole);

/**********************************************************************************************************************/
static int initMol(char* inputFileName, Molecule** molecule, int* pOrdre, double** pDeltas, OBForceField** ppFF, OBMol* pobMol, bool* pACKS2Charges)
{
	Molecule* cchemiMol = readMolecule(inputFileName,TRUE);
	double delta = 0.5;
	boolean reducedCoordinates = TRUE;
	int ordre = 2;
	int nAll;
	FILE* file = fopen(inputFileName,"rb");
	double* deltas = NULL;
	char* fileNamePrefix;
	char* fileNameHin;
	OBMol obMol;
	OBForceField* pFF = NULL;
	bool ACKS2Charges = FALSE;

	readOneBoolean(file,"QFFReducedCoordinates",&reducedCoordinates);
	
	readOneInt(file,"QFFnModes",&ordre);
	if(!reducedCoordinates) delta = 1e-2;
	readOneReal(file,"QFFDelta",&delta);
	fclose(file);
	printf("delta = %f ",delta);
	if(reducedCoordinates) printf("reducedCoordiantes\n");
	else printf(" Angshtrom\n");
	printf("nAtoms = %d\n",cchemiMol->nAtoms);
	printf("nModes = %d\n",cchemiMol->vibration.nModes);
	deltas = cchemiMol->klass->getDeltaTable(cchemiMol, delta, reducedCoordinates);


	if(cchemiMol) ACKS2Charges = setACKS2Parameters(inputFileName, cchemiMol);

	fileNamePrefix = getSuffixNameFile(inputFileName);
	fileNameHin = strdup_printf("%sOne.hin",fileNamePrefix);
	cchemiMol->klass->saveHIN(cchemiMol,fileNameHin);
	obMol =  ob(fileNameHin,pFF);
	if(resetCharges(inputFileName, obMol,  cchemiMol)) ACKS2Charges = FALSE;

	*pOrdre = ordre;
	*pDeltas = deltas;
	*molecule = cchemiMol;
	*ppFF = pFF;
	*pobMol = obMol;
	*pACKS2Charges = ACKS2Charges;
	//mol->klass->free(mol);
	nAll = 1 +6* cchemiMol->vibration.nModes;
	if(ordre>=2) nAll += 6*cchemiMol->vibration.nModes*(cchemiMol->vibration.nModes-1); 
	if(ordre>=3) nAll += 8*cchemiMol->vibration.nModes*(cchemiMol->vibration.nModes-1)*(cchemiMol->vibration.nModes-2)/6;
	if(ordre>=4) nAll += 16*cchemiMol->vibration.nModes*(cchemiMol->vibration.nModes-1)*(cchemiMol->vibration.nModes-2)*(cchemiMol->vibration.nModes-3)/24;
	return nAll;
}
/**********************************************************************************************************************/
static void computeEnergyAndDipole(Molecule* cchemiMol, OBForceField* pFF, OBMol& obMol, bool ACKS2Charges, double** geom, double* energies, double* dipoles[3], int k)
{
	double dipole[3];
	{
		int i=0;
		FOR_ATOMS_OF_MOL(atom, obMol)
		{
			atom->SetVector (geom[i][0], geom[i][1], geom[i][2]);
               		if(ACKS2Charges) 
			{
				cchemiMol->atoms[i].coordinates[0] = geom[i][0];
				cchemiMol->atoms[i].coordinates[1] = geom[i][1];
				cchemiMol->atoms[i].coordinates[2] = geom[i][2];
			}
			i++;
		}
		pFF->SetCoordinates(obMol);
		energies[k] = pFF->Energy(false);
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
		FOR_ATOMS_OF_MOL(atom, obMol) cerr << " x: " << atom->x() << " y: " << atom->y() << " z: " << atom->z() << endl;
		*/
	}
}
/**********************************************************************************************************************/
static void getEnergiesAndDipoles(Molecule* cchemiMol, double* deltas, OBForceField* pFF, OBMol& obMol, bool ACKS2Charges, int ordre, double* energies, double* dipoles[3])
{
	double** geom;
	Molecule* mol = cchemiMol;
	int i,j,k,l;
	int index;

	geom =  new double* [mol->nAtoms];
        for(i=0;i<mol->nAtoms;i++) geom[i] =  new double[3];

	index = 0;
	mol->klass->getQFFOneGeom(mol, -1, -1, 0.0, 0.0, 0.0,geom);
	computeEnergyAndDipole(mol, pFF, obMol, ACKS2Charges, geom, energies, dipoles, index);
	for(j=0;j<mol->vibration.nModes;j++)
	{
		index++; mol->klass->getQFFOneGeom(mol,  j, -1,   3*deltas[j], 0,0,geom);
		computeEnergyAndDipole(mol, pFF, obMol, ACKS2Charges, geom, energies, dipoles, index);

		index++; mol->klass->getQFFOneGeom(mol,  j, -1,   2*deltas[j], 0,0,geom);
		computeEnergyAndDipole(mol, pFF, obMol, ACKS2Charges, geom, energies, dipoles, index);

		index++; mol->klass->getQFFOneGeom(mol,  j, -1,     deltas[j], 0,0,geom);
		computeEnergyAndDipole(mol, pFF, obMol, ACKS2Charges, geom, energies, dipoles, index);

		index++; mol->klass->getQFFOneGeom(mol,  j, -1,    -deltas[j], 0,0,geom);
		computeEnergyAndDipole(mol, pFF, obMol, ACKS2Charges, geom, energies, dipoles, index);

		index++; mol->klass->getQFFOneGeom(mol,  j, -1,  -2*deltas[j], 0,0,geom);
		computeEnergyAndDipole(mol, pFF, obMol, ACKS2Charges, geom, energies, dipoles, index);

		index++; mol->klass->getQFFOneGeom(mol,  j, -1,  -3*deltas[j], 0,0,geom);
		computeEnergyAndDipole(mol, pFF, obMol, ACKS2Charges, geom, energies, dipoles, index);
	}
	if(ordre>=2)
	for(j=0;j<mol->vibration.nModes;j++)
	{
		for(i=0;i<j;i++)
		{
			index++; mol->klass->getQFFOneGeom(mol, j, i,  deltas[j],   deltas[i],0, geom);
			computeEnergyAndDipole(mol, pFF, obMol, ACKS2Charges, geom, energies, dipoles, index);

			index++; mol->klass->getQFFOneGeom(mol, j, i,  deltas[j],  -deltas[i],0, geom);
			computeEnergyAndDipole(mol, pFF, obMol, ACKS2Charges, geom, energies, dipoles, index);

			index++; mol->klass->getQFFOneGeom(mol, j, i, -deltas[j],   deltas[i],0, geom);
			computeEnergyAndDipole(mol, pFF, obMol, ACKS2Charges, geom, energies, dipoles, index);

			index++; mol->klass->getQFFOneGeom(mol, j, i, -deltas[j],  -deltas[i],0, geom);
			computeEnergyAndDipole(mol, pFF, obMol, ACKS2Charges, geom, energies, dipoles, index);
		}
		for(i=0;i<mol->vibration.nModes;i++)
		{
			if(i==j) continue;

			index++; mol->klass->getQFFOneGeom(mol,  j, i,  deltas[j],  3*deltas[i],0, geom);
			computeEnergyAndDipole(mol, pFF, obMol, ACKS2Charges, geom, energies, dipoles, index);

			index++; mol->klass->getQFFOneGeom(mol,  j, i,  deltas[j], -3*deltas[i],0, geom);
			computeEnergyAndDipole(mol, pFF, obMol, ACKS2Charges, geom, energies, dipoles, index);

			index++; mol->klass->getQFFOneGeom(mol,  j, i, -deltas[j],  3*deltas[i],0, geom);
			computeEnergyAndDipole(mol, pFF, obMol, ACKS2Charges, geom, energies, dipoles, index);

			index++; mol->klass->getQFFOneGeom(mol,  j, i, -deltas[j], -3*deltas[i],0, geom);
			computeEnergyAndDipole(mol, pFF, obMol, ACKS2Charges, geom, energies, dipoles, index);
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
			computeEnergyAndDipole(mol, pFF, obMol, ACKS2Charges, geom, energies, dipoles, index);
			/* 1  deltas[j],   deltas[i], -deltas[k] */
			index++; mol->klass->getQFFOneGeom3(mol, j, i, k,  deltas[j], deltas[i], -deltas[k], geom);
			computeEnergyAndDipole(mol, pFF, obMol, ACKS2Charges, geom, energies, dipoles, index);
			/* 2  deltas[j],  -deltas[i],  deltas[k] */
			index++; mol->klass->getQFFOneGeom3(mol, j, i, k,  deltas[j], -deltas[i], deltas[k], geom);
			computeEnergyAndDipole(mol, pFF, obMol, ACKS2Charges, geom, energies, dipoles, index);
			/* 3  deltas[j],  -deltas[i], -deltas[k] */
			index++; mol->klass->getQFFOneGeom3(mol, j, i, k,  deltas[j], -deltas[i], -deltas[k], geom);
			computeEnergyAndDipole(mol, pFF, obMol, ACKS2Charges, geom, energies, dipoles, index);
			/* 4 -deltas[j],  deltas[i],  deltas[k] */
			index++; mol->klass->getQFFOneGeom3(mol, j, i, k,  -deltas[j], deltas[i], deltas[k], geom);
			computeEnergyAndDipole(mol, pFF, obMol, ACKS2Charges, geom, energies, dipoles, index);
			/* 5 -deltas[j],  deltas[i], -deltas[k] */
			index++; mol->klass->getQFFOneGeom3(mol, j, i, k,  -deltas[j], deltas[i], -deltas[k], geom);
			computeEnergyAndDipole(mol, pFF, obMol, ACKS2Charges, geom, energies, dipoles, index);
			/* 6 -deltas[j], -deltas[i],  deltas[k] */
			index++; mol->klass->getQFFOneGeom3(mol, j, i, k,  -deltas[j], -deltas[i], deltas[k], geom);
			computeEnergyAndDipole(mol, pFF, obMol, ACKS2Charges, geom, energies, dipoles, index);
			/* 7 -deltas[j], -deltas[i], -deltas[k]- */
			index++; mol->klass->getQFFOneGeom3(mol, j, i, k,  -deltas[j], -deltas[i], -deltas[k], geom);
			computeEnergyAndDipole(mol, pFF, obMol, ACKS2Charges, geom, energies, dipoles, index);
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
			computeEnergyAndDipole(mol, pFF, obMol, ACKS2Charges, geom, energies, dipoles, index);
			index++; mol->klass->getQFFOneGeom4(mol, j, i, k, l,  deltas[j], deltas[i], deltas[k], -deltas[l], geom);
			computeEnergyAndDipole(mol, pFF, obMol, ACKS2Charges, geom, energies, dipoles, index);
			index++; mol->klass->getQFFOneGeom4(mol, j, i, k, l,  deltas[j], deltas[i], -deltas[k], deltas[l], geom);
			computeEnergyAndDipole(mol, pFF, obMol, ACKS2Charges, geom, energies, dipoles, index);
			index++; mol->klass->getQFFOneGeom4(mol, j, i, k, l,  deltas[j], deltas[i], -deltas[k], -deltas[l], geom);
			computeEnergyAndDipole(mol, pFF, obMol, ACKS2Charges, geom, energies, dipoles, index);
			index++; mol->klass->getQFFOneGeom4(mol, j, i, k, l,  deltas[j], -deltas[i], deltas[k], deltas[l], geom);
			computeEnergyAndDipole(mol, pFF, obMol, ACKS2Charges, geom, energies, dipoles, index);
			index++; mol->klass->getQFFOneGeom4(mol, j, i, k, l,  deltas[j], -deltas[i], deltas[k], -deltas[l], geom);
			computeEnergyAndDipole(mol, pFF, obMol, ACKS2Charges, geom, energies, dipoles, index);
			index++; mol->klass->getQFFOneGeom4(mol, j, i, k, l,  deltas[j], -deltas[i], -deltas[k], deltas[l], geom);
			computeEnergyAndDipole(mol, pFF, obMol, ACKS2Charges, geom, energies, dipoles, index);
			index++; mol->klass->getQFFOneGeom4(mol, j, i, k, l,  deltas[j], -deltas[i], -deltas[k], -deltas[l], geom);
			computeEnergyAndDipole(mol, pFF, obMol, ACKS2Charges, geom, energies, dipoles, index);

			index++; mol->klass->getQFFOneGeom4(mol, j, i, k, l, -deltas[j], deltas[i], deltas[k], deltas[l], geom);
			computeEnergyAndDipole(mol, pFF, obMol, ACKS2Charges, geom, energies, dipoles, index);
			index++; mol->klass->getQFFOneGeom4(mol, j, i, k, l, -deltas[j], deltas[i], deltas[k], -deltas[l], geom);
			computeEnergyAndDipole(mol, pFF, obMol, ACKS2Charges, geom, energies, dipoles, index);
			index++; mol->klass->getQFFOneGeom4(mol, j, i, k, l, -deltas[j], deltas[i], -deltas[k], deltas[l], geom);
			computeEnergyAndDipole(mol, pFF, obMol, ACKS2Charges, geom, energies, dipoles, index);
			index++; mol->klass->getQFFOneGeom4(mol, j, i, k, l, -deltas[j], deltas[i], -deltas[k], -deltas[l], geom);
			computeEnergyAndDipole(mol, pFF, obMol, ACKS2Charges, geom, energies, dipoles, index);
			index++; mol->klass->getQFFOneGeom4(mol, j, i, k, l, -deltas[j], -deltas[i], deltas[k], deltas[l], geom);
			computeEnergyAndDipole(mol, pFF, obMol, ACKS2Charges, geom, energies, dipoles, index);
			index++; mol->klass->getQFFOneGeom4(mol, j, i, k, l, -deltas[j], -deltas[i], deltas[k], -deltas[l], geom);
			computeEnergyAndDipole(mol, pFF, obMol, ACKS2Charges, geom, energies, dipoles, index);
			index++; mol->klass->getQFFOneGeom4(mol, j, i, k, l, -deltas[j], -deltas[i], -deltas[k], deltas[l], geom);
			computeEnergyAndDipole(mol, pFF, obMol, ACKS2Charges, geom, energies, dipoles, index);
			index++; mol->klass->getQFFOneGeom4(mol, j, i, k, l, -deltas[j], -deltas[i], -deltas[k], -deltas[l], geom);
			computeEnergyAndDipole(mol, pFF, obMol, ACKS2Charges, geom, energies, dipoles, index);
		}
	}
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
static OBMol ob(char*  filename, OBForceField*& pFF)
{
	int verbose = 1;
	string ff = "MMFF94";
	OBConversion conv;

	OBFormat *format_in = conv.FormatFromExt(filename);
    
	if (!format_in || !conv.SetInFormat(format_in))
	{
    		cerr << ": cannot read input format!" << endl;
    		exit (-1);
	}

	ifstream ifs;
	ofstream ofs;

	// Read the file
	ifs.open(filename);
	if (!ifs)
	{
		cerr << ": cannot read input file!" << endl;
		exit (-1);
	}

	pFF = OBForceField::FindForceField(ff);
	if (!pFF)
	{
		cerr << ": could not find forcefield '" << ff << "'." <<endl;
		exit (-1);
	}
	pFF->SetLogFile(&cout);
	if (verbose) pFF->SetLogLevel(OBFF_LOGLVL_HIGH);
	else pFF->SetLogLevel(OBFF_LOGLVL_NONE);

	OBMol mol;
	double energy;

	mol.Clear();
	if (!conv.Read(&mol, &ifs) || mol.Empty()) 
	{
		cerr << " I cannot read the molecule." << endl;
		exit (-1);
	}

	if (!pFF->Setup(mol))
	{
		cerr << ": could not setup force field." << endl;
		exit (-1);
	}
   
	//energy = pFF->Energy(true);// with gradient
	energy = pFF->Energy(false);
	/*
	if (!isfinite(energy))
	{
		cerr << " Title: " << mol.GetTitle() << endl;
		FOR_ATOMS_OF_MOL(atom, mol)
		{
			cerr << " x: " << atom->x() << " y: " << atom->y() << " z: " << atom->z() << endl;
		}
	}
	*/
	cout << "\nFINAL ENERGY:" << setw(20) << fixed << setprecision(14) << energy << " " <<pFF->GetUnit()<<endl;
	/*
	cout << "Gradients:" <<" "<<mol.NumAtoms()<< " atoms "<<endl;
	FOR_ATOMS_OF_MOL(atom, mol)
	{
		vector3 g = pFF->GetGradient(&*atom);
		cout 
		<< " " << setw(20) << fixed << setprecision(14) << -g[0]
		<< " " << setw(20) << fixed << setprecision(14) << -g[1]
		<< " " << setw(20) << fixed << setprecision(14) << -g[2]
		<< " " << setw(20) << atom->GetType()
		<<endl;
	}
	cout<<endl;
	*/
	pFF->SetLogLevel(OBFF_LOGLVL_NONE);
	return mol;
}
/****************************************************************************************/
static bool resetCharges(char* inputFileName, OBMol& obMol,  Molecule* mol)
{
	boolean mmff94Charges = FALSE; 
	FILE* file = fopen(inputFileName,"rb");
	string s = "mmff94Charges";
	char* tmp = strdup(s.c_str());
	readOneBoolean(file,tmp,&mmff94Charges);
        fclose(file);
	if(mmff94Charges) cerr<<"Warning : I use mmff94 charges"<<endl;
	else return FALSE;
	int i =0;
	OBChargeModel *mmffCharges = OBChargeModel::FindType("mmff94");
	std::vector<double> partialCharges;
	if (mmffCharges && mmffCharges->ComputeCharges(obMol)) {
		partialCharges = mmffCharges->GetPartialCharges();
		FOR_ATOMS_OF_MOL(atom, obMol) 
		{
		mol->atoms[i].charge = partialCharges[i];
		cerr <<  mol->atoms[i].prop.symbol<<" " 
		     <<" " << mol->atoms[i].coordinates[0]
		     <<" " << mol->atoms[i].coordinates[1]
		     <<" " << mol->atoms[i].coordinates[2]
		     <<" Charge = " << mol->atoms[i].charge<<endl;
		i++;
		}
	}
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
int obgvpt2(char* inputFileName)
{
	Molecule* cchemiMol = NULL;
	double*** geoms = NULL;
	char* fileNamePrefix;
	char* fileNameHin;
	char* fileNameDerivatives;
	OBMol obMol;
	OBForceField* pFF = NULL;
	int ordre;
	double* deltas = NULL;
	double dipole[3];
	int nAll = 0;
	double* energies = NULL;
	double* dipoles[3] =  {NULL, NULL, NULL};

        double elaps;
        double elapsAll = 0;
	double time0 = get_cpu_time();
	double time1 = 0;
	bool ACKS2Charges = FALSE;

	nAll = initMol(inputFileName, &cchemiMol, &ordre, &deltas,&pFF, &obMol,&ACKS2Charges);
	energies = new double [nAll];
	for(int xyz=0;xyz<3;xyz++)   dipoles[xyz] = new double [nAll];

	time1 = get_cpu_time();
	elaps = time1-time0;
	elapsAll += elaps;
	cerr<<"CPU time for initialisation "<<elaps<<" seconds"<<endl;

	time0 = time1;
	getEnergiesAndDipoles(cchemiMol, deltas, pFF, obMol, ACKS2Charges, ordre, energies, dipoles);
	time1 = get_cpu_time();
	elaps = time1-time0;
	elapsAll += elaps;
	cerr<<"CPU time for computing of energies and dipoles "<<elaps<<" seconds"<<endl;

	time0 = time1;
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

/**********************************************************************************************************************/
static double getQuartic(Molecule* cchemiMol, int i, int j, double* deltas, double V0, double** VI, double* VIJ)
{
	static double f4cm1 = AUTOCM1*AUTOCM1*AUTOCM1;
	double quarticEnergiesIIJJ = 0;
	{
		double massi = cchemiMol->vibration.modes[i].mass;
		double massj = cchemiMol->vibration.modes[j].mass;
		double freqi = cchemiMol->vibration.modes[i].frequency;
		double freqj = cchemiMol->vibration.modes[j].frequency;
		double mdi = massi*AMUTOAU*deltas[i]*deltas[i];
		double mdj = massj*AMUTOAU*deltas[j]*deltas[j];
		double f = 1.0/(mdi*mdj)*f4cm1;
		f = f/(freqi*freqj);
		{
			quarticEnergiesIIJJ = f*(
				   (VIJ[0]+VIJ[2]+VIJ[1]+VIJ[3])
				-2*(VI[i][2]+VI[i][3]+VI[j][2]+VI[j][3])
				+4*V0
			);
		}
	}
	return quarticEnergiesIIJJ*0.25;
}
/**********************************************************************************************************************/
static double getCubicIJJ(Molecule* cchemiMol, int i, int j, double* deltas, double V0, double** VI, double* VIJ)
{
	static double f3cm1 = AUTOCM1*sqrt(AUTOCM1)*AUTOCM1;
	double cubicTerm = 0;
	static double precision = 1e-12;
	// tijj
	if(i!=j)
        {
		double massi = cchemiMol->vibration.modes[i].mass;
		double massj = cchemiMol->vibration.modes[j].mass;
		double freqi = cchemiMol->vibration.modes[i].frequency;
		double freqj = cchemiMol->vibration.modes[j].frequency;
                double fj = 1.0/(2.0*deltas[j]*deltas[j]*massj*AMUTOAU);
                fj = fj/freqj;
                {
                        double fi,f;
                        fi = 1.0/(deltas[i]*sqrt(massi*AMUTOAU));
                        fi = fi/sqrt(freqi);
                        f = fi*fj*f3cm1;
                        cubicTerm = f*(
                                 (VIJ[0]+VIJ[1]-2*VI[i][2])
                                -(VIJ[2]+VIJ[3]-2*VI[i][3])
                        );
         
                }
		cubicTerm = cubicTerm*cubicTerm*0.25/fabs(freqi-2*freqj+precision);
        }
	return cubicTerm;
}
/**********************************************************************************************************************/
static double getCubicIJK(Molecule* cchemiMol, int i, int j, int k, double* deltas, double* VIJK)
{
	static double f3cm1 = AUTOCM1*sqrt(AUTOCM1)*AUTOCM1;
	double cubicTerm = 0;
	static double precision = 1e-12;
	// tijk
	if(i!=j && i!=k && j!=k)
        {
		double massi = cchemiMol->vibration.modes[i].mass;
		double massj = cchemiMol->vibration.modes[j].mass;
		double massk = cchemiMol->vibration.modes[k].mass;
		double freqi = cchemiMol->vibration.modes[i].frequency;
		double freqj = cchemiMol->vibration.modes[j].frequency;
		double freqk = cchemiMol->vibration.modes[k].frequency;
                double fi = 1.0/(deltas[i]*sqrt(massi*AMUTOAU));
                double fj = 1.0/(deltas[j]*sqrt(massj*AMUTOAU));
                double fk = 1.0/(deltas[k]*sqrt(massk*AMUTOAU));
                double f;
                fi /= sqrt(freqi);
                fj /= sqrt(freqj);
                fk /= sqrt(freqk);
                f = fi*fj*fk/8*f3cm1;

		cubicTerm = f*(
                                +VIJK[0]
                                -VIJK[1]
                                -VIJK[2]
                                +VIJK[3]
                                -VIJK[4]
                                +VIJK[5]
                                +VIJK[6]
                                -VIJK[7]
			);

         
		cubicTerm = cubicTerm*cubicTerm*0.25/fabs(freqi-freqj-freqk+precision);
        }
	return cubicTerm;
}
/***************************************************************************************************************/
static void setSelectedBy2CM(Molecule* cchemiMol, double* deltas, OBForceField* pFF, OBMol& obMol, bool ACKS2Charges, int ordre, boolean* targeted, boolean* activated,double Ecut,
double maxFrequencyDifferenceFermi)
{
	double** geom;
	Molecule* mol = cchemiMol;
	int i,j,k;
	int index;
	double* dipoles[3];
	double** VI;
	double VIJ[4];
	int nS = 0;
	double V0 = 0;

	if(Ecut<1e-10) return;

	for(j=0;j<mol->vibration.nModes;j++) if(activated[j]) nS++;
	if(nS==mol->vibration.nModes) return;

        for(i=0;i<3;i++) dipoles[i] =  new double[8];

	geom =  new double* [mol->nAtoms];
        for(i=0;i<mol->nAtoms;i++) geom[i] =  new double[3];

	VI =  new double* [mol->vibration.nModes];
	for(i=0;i<mol->vibration.nModes;i++) VI[i] =  new double[7];

	index = 0;
	mol->klass->getQFFOneGeom(mol, -1, -1, 0.0, 0.0, 0.0,geom);
	computeEnergyAndDipole(mol, pFF, obMol, ACKS2Charges, geom, &V0, dipoles, index);
	//cerr<<"End calcul V0"<<endl;

	for(j=0;j<mol->vibration.nModes;j++)
	{
		index = -1;
		index++; mol->klass->getQFFOneGeom(mol,  j, -1,   3*deltas[j], 0,0,geom);
		computeEnergyAndDipole(mol, pFF, obMol, ACKS2Charges, geom, VI[j], dipoles, index);

		index++; mol->klass->getQFFOneGeom(mol,  j, -1,   2*deltas[j], 0,0,geom);
		computeEnergyAndDipole(mol, pFF, obMol, ACKS2Charges, geom, VI[j], dipoles, index);

		index++; mol->klass->getQFFOneGeom(mol,  j, -1,     deltas[j], 0,0,geom);
		computeEnergyAndDipole(mol, pFF, obMol, ACKS2Charges, geom, VI[j], dipoles, index);

		index++; mol->klass->getQFFOneGeom(mol,  j, -1,    -deltas[j], 0,0,geom);
		computeEnergyAndDipole(mol, pFF, obMol, ACKS2Charges, geom, VI[j], dipoles, index);

		index++; mol->klass->getQFFOneGeom(mol,  j, -1,  -2*deltas[j], 0,0,geom);
		computeEnergyAndDipole(mol, pFF, obMol, ACKS2Charges, geom, VI[j], dipoles, index);

		index++; mol->klass->getQFFOneGeom(mol,  j, -1,  -3*deltas[j], 0,0,geom);
		computeEnergyAndDipole(mol, pFF, obMol, ACKS2Charges, geom, VI[j], dipoles, index);
	}
	//cerr<<"End calcul VI"<<endl;

	for(j=0;j<mol->vibration.nModes;j++)
	{
		if(!targeted[j]) continue;
		for(i=0;i<mol->vibration.nModes;i++)
		{
			if(activated[i]) continue;
			double freqi = mol->vibration.modes[i].frequency;
			double freqj = mol->vibration.modes[j].frequency;
			index = -1;
			index++; mol->klass->getQFFOneGeom(mol, j, i,  deltas[j],   deltas[i],0, geom);
			computeEnergyAndDipole(mol, pFF, obMol, ACKS2Charges, geom, VIJ, dipoles, index);

			index++; mol->klass->getQFFOneGeom(mol, j, i,  deltas[j],  -deltas[i],0, geom);
			computeEnergyAndDipole(mol, pFF, obMol, ACKS2Charges, geom, VIJ, dipoles, index);

			index++; mol->klass->getQFFOneGeom(mol, j, i, -deltas[j],   deltas[i],0, geom);
			computeEnergyAndDipole(mol, pFF, obMol, ACKS2Charges, geom, VIJ, dipoles, index);

			index++; mol->klass->getQFFOneGeom(mol, j, i, -deltas[j],  -deltas[i],0, geom);
			computeEnergyAndDipole(mol, pFF, obMol, ACKS2Charges, geom, VIJ, dipoles, index);
			double E = getQuartic(cchemiMol, j, i, deltas, V0, VI, VIJ);
			//if(fabs(E)>Ecut) cerr<<" i="<<j<<" j="<<i<<" E="<<E<<endl;
			if(fabs(E)>Ecut) { activated[i] = true; continue;}
			if(fabs(freqj-2*freqi)>maxFrequencyDifferenceFermi)  continue;
			E =  getCubicIJJ(cchemiMol, j, i, deltas, V0, VI, VIJ);
			//cerr<<" Cubic term"<< " E = "<<E<<endl;
			if(fabs(E)>Ecut)
			{ 
				cerr<<" mode j activated by tijj term"<< " E = "<<E<<endl;
				activated[i] = true; 
				continue;
			}
		}
	}

	//cerr<<"End calcul VIJ"<<endl;
        for(i=0;i<3;i++) free(dipoles[i]);
        for(i=0;i<mol->nAtoms;i++) free(geom[i]);
        free(geom);
	//cerr<<"End free geom"<<endl;
	for(i=0;i<mol->vibration.nModes;i++) free(VI[i]);
	free(VI);
}
/***************************************************************************************************************/
static void setSelectedBy3CM(Molecule* cchemiMol, double* deltas, OBForceField* pFF, OBMol& obMol, bool ACKS2Charges, int ordre, boolean* targeted, boolean* activated,double Ecut, double maxFrequencyDifferenceFermi)
{
	double** geom;
	Molecule* mol = cchemiMol;
	int i,j,k;
	int index;
	double* dipoles[3];
	double VIJK[8];
	int nS = 0;

	if(Ecut<1e-10) return;

	cerr<<" Selection by 3MC, please wait....."<<endl;
	for(j=0;j<mol->vibration.nModes;j++) if(activated[j]) nS++;
	if(nS==mol->vibration.nModes) return;

        for(i=0;i<3;i++) dipoles[i] =  new double[8];

	geom =  new double* [mol->nAtoms];
        for(i=0;i<mol->nAtoms;i++) geom[i] =  new double[3];

	for(j=0;j<mol->vibration.nModes;j++)
	{
		if(!targeted[j]) continue;
		for(i=0;i<mol->vibration.nModes;i++)
		{
			if(i==j) continue;
          		for(k=0;k<mol->vibration.nModes;k++)
                	{
				if(i==k) continue;
				if(j==k) continue;
				if(activated[i] && activated[k]) continue;
				double freqi = mol->vibration.modes[i].frequency;
				double freqj = mol->vibration.modes[j].frequency;
				double freqk = mol->vibration.modes[k].frequency;
				if(fabs(freqj-freqi-freqk)>maxFrequencyDifferenceFermi)  continue;

				index = -1;
                        	/* 0  deltas[j],   deltas[i],  deltas[k] */
                        	index++; mol->klass->getQFFOneGeom3(mol, j, i, k,  deltas[j], deltas[i], deltas[k], geom);
                        	computeEnergyAndDipole(mol, pFF, obMol, ACKS2Charges, geom, VIJK, dipoles, index);

                        	/* 1  deltas[j],   deltas[i], -deltas[k] */
                        	index++; mol->klass->getQFFOneGeom3(mol, j, i, k,  deltas[j], deltas[i], -deltas[k], geom);
                        	computeEnergyAndDipole(mol, pFF, obMol, ACKS2Charges, geom, VIJK, dipoles, index);

                        	/* 2  deltas[j],  -deltas[i],  deltas[k] */
                        	index++; mol->klass->getQFFOneGeom3(mol, j, i, k,  deltas[j], -deltas[i], deltas[k], geom);
                        	computeEnergyAndDipole(mol, pFF, obMol, ACKS2Charges, geom, VIJK, dipoles, index);

                        	/* 3  deltas[j],  -deltas[i], -deltas[k] */
                        	index++; mol->klass->getQFFOneGeom3(mol, j, i, k,  deltas[j], -deltas[i], -deltas[k], geom);
                        	computeEnergyAndDipole(mol, pFF, obMol, ACKS2Charges, geom, VIJK, dipoles, index);

                        	/* 4 -deltas[j],  deltas[i],  deltas[k] */
                        	index++; mol->klass->getQFFOneGeom3(mol, j, i, k,  -deltas[j], deltas[i], deltas[k], geom);
                        	computeEnergyAndDipole(mol, pFF, obMol, ACKS2Charges, geom, VIJK, dipoles, index);

                        	/* 5 -deltas[j],  deltas[i], -deltas[k] */
                        	index++; mol->klass->getQFFOneGeom3(mol, j, i, k,  -deltas[j], deltas[i], -deltas[k], geom);
                        	computeEnergyAndDipole(mol, pFF, obMol, ACKS2Charges, geom, VIJK, dipoles, index);

                        	/* 6 -deltas[j], -deltas[i],  deltas[k] */
                        	index++; mol->klass->getQFFOneGeom3(mol, j, i, k,  -deltas[j], -deltas[i], deltas[k], geom);
                        	computeEnergyAndDipole(mol, pFF, obMol, ACKS2Charges, geom, VIJK, dipoles, index);

                        	/* 7 -deltas[j], -deltas[i], -deltas[k]- */
                        	index++; mol->klass->getQFFOneGeom3(mol, j, i, k,  -deltas[j], -deltas[i], -deltas[k], geom);
                        	computeEnergyAndDipole(mol, pFF, obMol, ACKS2Charges, geom, VIJK, dipoles, index);

				double E =  getCubicIJK(cchemiMol, j, i, k, deltas, VIJK);
				//cerr<<" Cubic term ijk"<< " E = "<<E<< "\t\twi-wj-wk="<<freqj-freqi-freqk<<endl;
				if(fabs(E)>Ecut)
				{ 
					cerr<<" mode j,k activated by cubic ijk term"<< " E = "<<E<<endl;
					activated[i] = true; 
					activated[k] = true; 
				}
                	}
		}
	}

	//cerr<<"End calcul VIJ"<<endl;
        for(i=0;i<3;i++) free(dipoles[i]);
        for(i=0;i<mol->nAtoms;i++) free(geom[i]);
        free(geom);
	//cerr<<"End free geom"<<endl;
}
/**********************************************************************************************************************/
/*
static void setSelectedByFermi(Molecule* cchemiMol, boolean* targeted, boolean* activated,double Ecut)
{
	int i,j,k;

	for(i=0;i<cchemiMol->vibration.nModes;i++)
	{
		if(!targeted[i]) continue;
		for(j=0;j<cchemiMol->vibration.nModes;j++)
		for(k=0;k<cchemiMol->vibration.nModes;k++)
		{
			if(i==j) continue;
			if(i==k) continue;
			if(j==k) continue;
			double freqi = cchemiMol->vibration.modes[i].frequency;
			double freqj = cchemiMol->vibration.modes[j].frequency;
			double freqk = cchemiMol->vibration.modes[k].frequency;
			if(fabs(freqi-freqj-freqk)<=Ecut) { activated[j] = true; activated[k] = true; }
		}
	}
	for(i=0;i<cchemiMol->vibration.nModes;i++)
	{
		if(!targeted[i]) continue;
		for(j=0;j<cchemiMol->vibration.nModes;j++)
		{
			if(i==j) continue;
			double freqi = cchemiMol->vibration.modes[i].frequency;
			double freqj = cchemiMol->vibration.modes[j].frequency;
			if(fabs(freqi-2*freqj)<=Ecut) { activated[j] = true; }
		}
	}
}
*/
/**********************************************************************************************************************/
static void setSelectedByDarlingDennison(Molecule* cchemiMol, boolean* targeted, boolean* activated,double Ecut)
{
	int i,j;

	for(i=0;i<cchemiMol->vibration.nModes;i++)
	{
		//if(!targeted[i]) continue;
		if(!activated[i]) continue;
		for(j=0;j<cchemiMol->vibration.nModes;j++)
		{
			if(i==j) continue;
			double freqi = cchemiMol->vibration.modes[i].frequency;
			double freqj = cchemiMol->vibration.modes[j].frequency;
			if(fabs(freqi-freqj)<=Ecut) { activated[j] = true; }
		}
	}
}
/**********************************************************************************************************************/
static void initActivatedSelectedModes(char* inputFileName, Molecule* cchemiMol, boolean* targeted, boolean* activated,
	double& ECut2CM,
	double& ECut3CM,
	double& ECutDarlingDennison,
	double& maxFrequencyDifferenceFermi)
{
	double ETargetMin = -1;
	double ETargetMax = -1;
	int ITargetMin = -1;
	int ITargetMax = -1;
	FILE* file = fopen(inputFileName,"rb");

	readOneReal(file,"ETargetMin",&ETargetMin);
	readOneReal(file,"ETargetMax",&ETargetMax);
	readOneInt(file,"ITargetMin",&ITargetMin);
	readOneInt(file,"ITargetMax",&ITargetMax);

	ECut2CM = 10.0;
	ECutDarlingDennison = 1;
	readOneReal(file,"ECut2CM",&ECut2CM);
	ECut3CM = ECut2CM;
	readOneReal(file,"ECut3CM",&ECut3CM);
	readOneReal(file,"ECutDarlingDennison",&ECutDarlingDennison);
	readOneReal(file,"maxFrequencyDifferenceFermi",&maxFrequencyDifferenceFermi);


	fclose(file);
	int nModes = cchemiMol->vibration.nModes;

	for(int i=0;i<nModes;i++) targeted[i]=false;

	if(ITargetMin<0 && ITargetMax<0 && ETargetMin<0 && ETargetMax<0) 
	{
		for(int i=0;i<nModes;i++) targeted[i]=true;
	}
	else if(ITargetMin>0 || ITargetMax> 0) 
	{
		ITargetMin--;
		ITargetMax--;
		if(ITargetMin<0) ITargetMin = 0;
		if(ITargetMax<0) ITargetMax = nModes-1;
		for(int i=ITargetMin;i<=ITargetMax;i++) targeted[i]=true;
	}
	else if(ETargetMin>=0 || ETargetMax>=0) 
	{
		if(ETargetMin<0) 
		{
			ETargetMin = cchemiMol->vibration.modes[0].frequency;
			for(int i=1;i<nModes;i++) 
				if(ETargetMin>cchemiMol->vibration.modes[i].frequency)
					ETargetMin = cchemiMol->vibration.modes[i].frequency;
		}
		if(ETargetMax<0) 
		{
			ETargetMax = cchemiMol->vibration.modes[0].frequency;
			for(int i=1;i<nModes;i++) 
				if(ETargetMax<cchemiMol->vibration.modes[i].frequency)
					ETargetMax = cchemiMol->vibration.modes[i].frequency;
		}
		for(int i= 0;i<nModes;i++) 
			if(
			cchemiMol->vibration.modes[i].frequency>=ETargetMin 
			&&
			cchemiMol->vibration.modes[i].frequency<=ETargetMax
			)
			targeted[i]=true;
	}
	for(int i=0;i<nModes;i++) activated[i]=targeted[i];
	return;
}
/**********************************************************************************************************************/
static int getNumberOfTrue(boolean* states, int nStates)
{
	int n = 0;

	for(int i=0;i<nStates;i++)
		if(states[i])n++;

	return n;
}
/**********************************************************************************************************************/
int selectModesForTargetModes(char* inputFileName)
{
	Molecule* cchemiMol = NULL;
	OBMol obMol;
	OBForceField* pFF = NULL;
	int ordre;
	double* deltas = NULL;
	int nAll = 0;
	boolean* targeted = NULL;
	boolean* activated = NULL;
	int nModes = 0;

        double elaps;
        double elapsAll = 0;
	double time0 = get_cpu_time();
	double time1 = 0;
	bool ACKS2Charges = FALSE;
	double ECut2CM = 10.0;
	double ECut3CM = 10.0;
	double ECutDarlingDennison = 1;
	double maxFrequencyDifferenceFermi = 200.0;// cm-1

	nAll = initMol(inputFileName, &cchemiMol, &ordre, &deltas,&pFF, &obMol,&ACKS2Charges);

	nModes = cchemiMol->vibration.nModes;
	targeted = new boolean [nModes];
	activated = new boolean [nModes];
	initActivatedSelectedModes(inputFileName, cchemiMol, targeted, activated, ECut2CM, ECut3CM, ECutDarlingDennison, maxFrequencyDifferenceFermi);
	cout<<"Number of targeted modes "<<getNumberOfTrue(targeted, nModes)
	<<" from "<<nModes<< " modes "<<endl;

	time1 = get_cpu_time();
	elaps = time1-time0;
	elapsAll += elaps;
	cerr<<"CPU time for initialisation "<<elaps<<" seconds"<<endl;

	time0 = time1;
	setSelectedBy2CM(cchemiMol, deltas, pFF, obMol, ACKS2Charges, ordre, targeted, activated,ECut2CM,maxFrequencyDifferenceFermi);
	cout<<"Number of activated modes by 2CM "<<getNumberOfTrue(activated, nModes)
	<<" from "<<nModes<< " modes "<<endl;
	time1 = get_cpu_time();
	elaps = time1-time0;
	elapsAll += elaps;
	cerr<<"CPU time for computing of activated modes by 2CM "<<elaps<<" seconds"<<endl;

	time0 = time1;
	//setSelectedByFermi(cchemiMol, targeted, activated,ECutFermi);
	setSelectedBy3CM(cchemiMol, deltas, pFF, obMol, ACKS2Charges, ordre,targeted, activated,ECut3CM,maxFrequencyDifferenceFermi);
	cout<<"Number of activated modes by 2CM + 3CM "<<getNumberOfTrue(activated, nModes)
	<<" from "<<nModes<< " modes "<<endl;
	time1 = get_cpu_time();
	elaps = time1-time0;
	elapsAll += elaps;
	cerr<<"CPU time for computing of activated modes by 3CM resonances "<<elaps<<" seconds"<<endl;

	time0 = time1;
	setSelectedByDarlingDennison(cchemiMol, targeted, activated,ECutDarlingDennison);
	cout<<"Number of activated modes by 2CM + 3CM + Darling-Dennison resonances "<<getNumberOfTrue(activated, nModes)
	<<" from "<<nModes<< " modes "<<endl;
	time1 = get_cpu_time();
	elaps = time1-time0;
	elapsAll += elaps;
	cerr<<"CPU time for computing of activated modes by Darling-Dennison resonances "<<elaps<<" seconds"<<endl;

	cerr<<endl<<"All time "<<elapsAll<<" seconds"<<endl;
	cchemiMol->klass->removeNotSelectedFrequencies(cchemiMol, activated);
	char* fileNameOut = strdup_printf("%sSelected.gab",getSuffixNameFile(inputFileName));
	cchemiMol->klass->save(cchemiMol, fileNameOut);

	/* 
	finalize();
	*/

	return 0;
}
