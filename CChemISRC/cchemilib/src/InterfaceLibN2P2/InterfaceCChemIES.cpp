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

// writen by A.R. Allouche based on InterfaceCChemIES of n2p2 software

//#ifndef NOMPI
//#include <mpi.h>
//#include "mpi-extra.h"
//#endif
#include "InterfaceCChemIES.h"
#include "Atom.h"
#include "Element.h"
#include "utility.h"
#include <limits>    // std::numeric_limits
#include <cmath>
#include <string>
#include <cstring>
#include <cstdlib>
#include <iostream>
#include "../Molecule/Molecule.h"
#include "../Utils/Constants.h"

#ifdef HIGH_DERIVATIVES
#include "Derivatives.h"
#include "DerivativesIC.h"
#endif

#define TOLCUTOFF 1.0E-2

using namespace std;
using namespace nnp;

InterfaceCChemIES::InterfaceCChemIES() : 
                                     initialized (false),
                                     cflength    (1.0  ),
                                     cfdipole    (1.0  ),
				     myRank      (0    )
{
	directory = string(getenv("PWD"));
}
InterfaceCChemIES::InterfaceCChemIES(string directory,
                                 double       cflength,
                                 double       cfdipole,
                                 int          myRank)
{
	initialize(directory, cflength, cfdipole, myRank);
}

void InterfaceCChemIES::initialize(string directory,
                                 double       cflength,
                                 double       cfdipole,
                                 int          myRank)
{
    this->directory = directory;
    this->cflength = cflength;
    this->cfdipole = cfdipole;
    this->myRank = myRank;
    log.writeToStdout = false;
    string dir(directory);
    Mode::initialize();
    loadSettingsFile(dir + "/inputES.nn");
    setupGeneric();
    setupSymmetryFunctionScalingNone();
    setupSymmetryFunctionStatistics(false, false, false, false);
    setupNeuralNetworkWeights(dir + "/weightsES.%03d.data");

    structure.setElementMap(elementMap);

    initialized = true;
}
void InterfaceCChemIES::readStructureFromFile(char* const& fileName)
{
    ifstream file;
    file.open(fileName);
    structure.reset();
    structure.setElementMap(elementMap);
    structure.readFromFile(string(fileName));
    file.close();

    return;
}


void InterfaceCChemIES::process()
{
#ifdef NOSFGROUPS
    calculateSymmetryFunctions(structure, false);
#else
    calculateSymmetryFunctionGroups(structure, false);
#endif
    calculateAtomicChargesNeuralNetworks(structure,false);
    calculateCharge(structure);
    calculateDipole(structure);

    return;
}
bool InterfaceCChemIES::computeChargeAndDipole(int nAtoms, char** symbols, double** coordinates, double* charge, double* dipole)
{
    double cd = 1.0/(AUTODEB);
    double cl = ANGTOBOHR;

    if(!initialized) initialize(directory,cl, cd, myRank);
    if(structure.atoms.size() != nAtoms) setStructure(nAtoms, symbols, coordinates);
    else setCoordinates(nAtoms, symbols, coordinates);

    structure.calculateNeighborList(maxCutoffRadius);
#ifdef NOSFGROUPS
    calculateSymmetryFunctions(structure, false);
#else
    calculateSymmetryFunctionGroups(structure, false);
#endif
    calculateAtomicChargesNeuralNetworks(structure,false);
    if(charge) 
    {
	calculateCharge(structure);
        *charge = structure.charge;
    }
    if(dipole) 
    {
        calculateDipole(structure);
    	for(size_t i=0;i<3;i++) dipole[i] = structure.dipole[i] / cfdipole;
    }
    return true;
}
#ifdef HIGH_DERIVATIVES
bool InterfaceCChemIES::computeDipoleHighDerivatives(Molecule* mol, Derivatives*& deriv, int order, int method)
{
//	cerr<<"Begin Grad"<<endl;
    double cd = 1.0/(AUTODEB);
    double cl = ANGTOBOHR;
    // Conversion in atomic unit
    if(!initialized) initialize(directory,cl, cd, myRank);


//	cerr<<"Begin setSC"<<endl;

    if(structure.atoms.size() != mol->nAtoms) setStructureFromMolecule(mol);
    else setCoordinatesFromMolecule(mol);
// DEBUG setStructureFromMolecule(mol);

   if (normalize) structure.toNormalizedUnits(meanEnergy, convEnergy, convLength,convDipole);
	//cerr<<"Begin calculateNeighborList"<<endl;

    structure.calculateNeighborList(maxCutoffRadius);
	//cerr<<"End calculateNeighborList"<<endl;
    if(!structure.hasSymmetryFunctions)
    {
#ifdef NOSFGROUPS
	//cerr<<"Begin calculateSymmetryFunctions("<<endl;
    calculateSymmetryFunctions(structure, true);
#else
	//cerr<<"Begin calculateSymmetryFunctionsGroup("<<endl;
    calculateSymmetryFunctionGroups(structure, true);
#endif
    }
	//cerr<<"Begin nnp"<<endl;
    calculateAtomicChargesNeuralNetworks(structure, true);
    calculateCharge(structure);
    calculateDipole(structure);
    deriv = Mode::computeHighDerivativesDipole(structure, order, method);


    if (normalize && deriv)
    {
    	for(int c=0;c<3;c++)
		deriv[c].toPhysicalUnits(convDipole, convLength);

        structure.toPhysicalUnits(meanEnergy, convEnergy, convLength, convDipole);
    }
    // conversion in cchemilib unit : Debye, Ang
    for(int c=0;c<3;c++)
    	deriv[c].toPhysicalUnits(cfdipole, cflength);
    return true;
}
#endif
double InterfaceCChemIES::computeChargeAndDipole(Molecule* mol)
{
    double cd = 1.0/(AUTODEB);
    double cl = ANGTOBOHR;

    if(!initialized) initialize(directory,cl, cd, myRank);

    if(structure.atoms.size() != mol->nAtoms) setStructureFromMolecule(mol);
    else setCoordinatesFromMolecule(mol);

	/*
    for (size_t i =  0; i < structure.atoms.size(); ++i)
    {
	for(size_t c=0;c<3;c++) cout<<structure.atoms[i].r[c]<<" ";
	cout<<endl;
	for(size_t c=0;c<3;c++) cout<<mol->atoms[i].coordinates[c]*cflength<<" ";
	cout<<endl;
	for(size_t c=0;c<3;c++) cout<<mol->atoms[i].coordinates[c]<<" ";
	cout<<endl;
    }
	*/

    structure.calculateNeighborList(maxCutoffRadius);
#ifdef NOSFGROUPS
    calculateSymmetryFunctions(structure, false);
#else
    calculateSymmetryFunctionGroups(structure, false);
#endif
    calculateAtomicChargesNeuralNetworks(structure,false);
    calculateCharge(structure);
    calculateDipole(structure);
    for(size_t i=0;i<3;i++) mol->dipole[i] = structure.dipole[i] / cfdipole;
    size_t j =0;
    for (size_t i =  0; i < structure.atoms.size(); ++i,++j)
    {
	if(!strcmp(mol->atoms[j].prop.symbol,"TV") || !strcmp(mol->atoms[j].prop.symbol,"Tv"))
	{
		i--;
		continue;
	}
	mol->atoms[j].charge = structure.atoms[i].charge;
    }
    return structure.charge;
}
double InterfaceCChemIES::computedDipole(Molecule* mol, double*** dmu)
{
    double cd = 1.0/(AUTODEB);
    double cl = ANGTOBOHR;

    if(!initialized) initialize(directory,cl, cd, myRank);

    if(structure.atoms.size() != mol->nAtoms) setStructureFromMolecule(mol);
    else setCoordinatesFromMolecule(mol);

    structure.calculateNeighborList(maxCutoffRadius);
#ifdef NOSFGROUPS
    calculateSymmetryFunctions(structure, true);// true for derivative of dG/dxi
#else
    calculateSymmetryFunctionGroups(structure, true);
#endif
    calculateAtomicChargesNeuralNetworks(structure,true);// true for derivative dQ/dG
    calculateCharge(structure);
    calculateDipole(structure);
    *dmu = Mode::computedDipole(structure);
    for(size_t i=0;i<3;i++) mol->dipole[i] = structure.dipole[i] / cfdipole;
    size_t j =0;
    for (size_t i =  0; i < structure.atoms.size(); ++i,++j)
    {
	if(!strcmp(mol->atoms[j].prop.symbol,"TV") || !strcmp(mol->atoms[j].prop.symbol,"Tv"))
	{
		i--;
		continue;
	}
	mol->atoms[j].charge = structure.atoms[i].charge;
    }
    for (size_t i =  0; i < structure.atoms.size(); ++i) 
    	for (size_t l =  0; l < 3; ++l) 
    	for (size_t k =  0; k < 3; ++k) 
		(*dmu)[k][3*i+l] *= cflength/cfdipole;

    return structure.charge;
}
double InterfaceCChemIES::run(Molecule* mol, string directory, double cflength,double cfdipole, int myRank)
{
    initialize(directory, cflength, cfdipole,myRank);
    for (size_t i =  0; i < structure.atoms.size(); ++i)
    {
	for(size_t c=0;c<3;c++) structure.atoms[i].r[c] = mol->atoms[i].coordinates[c]*cflength;
    }
    structure.calculateNeighborList(maxCutoffRadius);
#ifdef NOSFGROUPS
    calculateSymmetryFunctions(structure, false);
#else
    calculateSymmetryFunctionGroups(structure, false);
#endif
    calculateAtomicChargesNeuralNetworks(structure,false);
    calculateCharge(structure);
    calculateDipole(structure);
    for(size_t i=0;i<3;i++) mol->dipole[i] = structure.dipole[i] / cfdipole;
    size_t j =0;
    for (size_t i =  0; i < structure.atoms.size(); ++i,++j)
    {
	if(!strcmp(mol->atoms[j].prop.symbol,"TV") || !strcmp(mol->atoms[j].prop.symbol,"Tv"))
	{
		i--;
		continue;
	}
	mol->atoms[j].charge = structure.atoms[i].charge;
    }
    return structure.charge;
    return true;
}
void InterfaceCChemIES::setStructure(int nAtoms, char** symbols, double** coordinates)
{
    size_t         iBoxVector = 0;
    size_t         i = 0;
    size_t         c = 0;
    double totalCharge = 0;
    structure.reset();
    structure.setElementMap(elementMap);
    for(size_t i=0;i<nAtoms;i++)
    {
	if(!strcmp(symbols[i],"TV") || !strcmp(symbols[i],"Tv"))
	{
            	if (iBoxVector > 2)
            	{
                	throw runtime_error("ERROR: Too many box vectors.\n");
            	}
    		for(size_t c=0;c<3;c++)  structure.box[iBoxVector][c] = coordinates[c][i]*cflength;
            	iBoxVector++;
            	if (iBoxVector == 3)
            	{
                	 structure.isPeriodic = true;
                	if (structure.box[0][1] > numeric_limits<double>::min() ||
                    	 structure.box[0][2] > numeric_limits<double>::min() ||
                    	 structure.box[1][0] > numeric_limits<double>::min() ||
                    	 structure.box[1][2] > numeric_limits<double>::min() ||
                    	 structure.box[2][0] > numeric_limits<double>::min() ||
                    	 structure.box[2][1] > numeric_limits<double>::min())
                	{
                    		 structure.isTriclinic = true;
                	}
                	 structure.calculateInverseBox();
                	 structure.calculateVolume();
            	}
	}
	else
	{
            	 structure.atoms.push_back(Atom());
            	 structure.atoms.back().index =  structure.numAtoms;
            	 structure.atoms.back().indexStructure =  structure.index;
            	 structure.atoms.back().tag            =  structure.numAtoms;
    		for(size_t c=0;c<3;c++)  structure.atoms.back().r[c]    =  coordinates[c][i]*cflength;
            	 structure.atoms.back().element        =  structure.elementMap[symbols[i]];// a changer
            	 structure.atoms.back().charge         = 0.0;
    		for(size_t c=0;c<3;c++)  structure.atoms.back().fRef[c] = 0.0;
            	 structure.atoms.back().numNeighborsPerElement.resize(structure.numElements, 0);
            	structure.numAtoms++;
            	 structure.numAtomsPerElement[elementMap[symbols[i]]]++;// A changer
    		totalCharge += structure.atoms.back().charge;
        }
    }
    structure.energyRef = 0.0;
    structure.chargeRef = totalCharge;
    structure.charge = 0;
    for(size_t i=0;i<3;i++) structure.dipoleRef[i] = 0;
    for(size_t i=0;i<3;i++) structure.dipole[i] = 0;

    for (size_t i = 0; i < structure.numElements; i++)
    {
        if (structure.numAtomsPerElement[i] > 0)
        {
            structure.numElementsPresent++;
        }
    }

    if (structure.isPeriodic)
    {
        for (vector<Atom>::iterator it = structure.atoms.begin(); it != structure.atoms.end(); ++it)
        {
            structure.remap((*it));
        }
    }
}
void InterfaceCChemIES::setCoordinates(int nAtoms, char** symbols, double** coordinates)
{
    int iBoxVector = 0;
    Atom* ai = NULL;
    int ia=0;
    structure.clearNeighborList();
    for(size_t i=0;i<nAtoms;i++)
    {
	if(!strcmp(symbols[i],"TV") || !strcmp(symbols[i],"Tv"))
	{
            	if (iBoxVector > 2)
            	{
                	throw runtime_error("ERROR: Too many box vectors.\n");
            	}
    		for(size_t c=0;c<3;c++)  structure.box[iBoxVector][c] = coordinates[c][i]*cflength;
            	iBoxVector++;
            	if (iBoxVector == 3)
            	{
                	 structure.isPeriodic = true;
                	if (structure.box[0][1] > numeric_limits<double>::min() ||
                    	 structure.box[0][2] > numeric_limits<double>::min() ||
                    	 structure.box[1][0] > numeric_limits<double>::min() ||
                    	 structure.box[1][2] > numeric_limits<double>::min() ||
                    	 structure.box[2][0] > numeric_limits<double>::min() ||
                    	 structure.box[2][1] > numeric_limits<double>::min())
                	{
                    		 structure.isTriclinic = true;
                	}
                	 structure.calculateInverseBox();
                	 structure.calculateVolume();
            	}
	}
	else
	{
        	// Set pointer to atom.
        	ai = &(structure.atoms.at(ia));
    		for(size_t c=0;c<3;c++)  ai->r[c] = coordinates[c][i]*cflength;
            	 ai->numNeighborsPerElement.resize(structure.numElements, 0);
    		for(size_t c=0;c<3;c++)  ai->fRef[c] = 0.0;
		ia++;
        }
    }
    structure.energyRef = 0;
    for(size_t i=0;i<3;i++) structure.dipoleRef[i] = 0;
    for(size_t i=0;i<3;i++) structure.dipole[i] = 0;
}
void InterfaceCChemIES::setCoordinatesFromMolecule(Molecule* mol)
{
    int iBoxVector = 0;
    Atom* ai = NULL;
    int ia=0;
	structure.clearNeighborList();
    for(size_t i=0;i<mol->nAtoms;i++)
    {
	if(!strcmp(mol->atoms[i].prop.symbol,"TV") || !strcmp(mol->atoms[i].prop.symbol,"Tv"))
	{
            	if (iBoxVector > 2)
            	{
                	throw runtime_error("ERROR: Too many box vectors.\n");
            	}
    		for(size_t c=0;c<3;c++)  structure.box[iBoxVector][c] = mol->atoms[i].coordinates[c]*cflength;
            	iBoxVector++;
            	if (iBoxVector == 3)
            	{
                	 structure.isPeriodic = true;
                	if (structure.box[0][1] > numeric_limits<double>::min() ||
                    	 structure.box[0][2] > numeric_limits<double>::min() ||
                    	 structure.box[1][0] > numeric_limits<double>::min() ||
                    	 structure.box[1][2] > numeric_limits<double>::min() ||
                    	 structure.box[2][0] > numeric_limits<double>::min() ||
                    	 structure.box[2][1] > numeric_limits<double>::min())
                	{
                    		 structure.isTriclinic = true;
                	}
                	 structure.calculateInverseBox();
                	 structure.calculateVolume();
            	}
	}
	else
	{
        	// Set pointer to atom.
        	ai = &(structure.atoms.at(ia));
    		for(size_t c=0;c<3;c++)  ai->r[c] = mol->atoms[i].coordinates[c]*cflength;
    		ai->charge = mol->atoms[i].charge;
            	 ai->numNeighborsPerElement.resize(structure.numElements, 0);
    		for(size_t c=0;c<3;c++)  ai->fRef[c] = 0.0;
		ia++;
        }
    }
    for(size_t i=0;i<3;i++) structure.dipoleRef[i] = mol->dipole[i];
    for(size_t i=0;i<3;i++) structure.dipole[i] = 0;
}
void InterfaceCChemIES::setStructureFromMolecule(Molecule* mol)
{
    size_t         iBoxVector = 0;
    size_t         i = 0;
    size_t         c = 0;
    double totalCharge = 0;
    structure.reset();
    structure.setElementMap(elementMap);
    for(size_t i=0;i<mol->nAtoms;i++)
    {
	if(!strcmp(mol->atoms[i].prop.symbol,"TV") || !strcmp(mol->atoms[i].prop.symbol,"Tv"))
	{
            	if (iBoxVector > 2)
            	{
                	throw runtime_error("ERROR: Too many box vectors.\n");
            	}
    		for(size_t c=0;c<3;c++)  structure.box[iBoxVector][c] = mol->atoms[i].coordinates[c]*cflength;
            	iBoxVector++;
            	if (iBoxVector == 3)
            	{
                	 structure.isPeriodic = true;
                	if (structure.box[0][1] > numeric_limits<double>::min() ||
                    	 structure.box[0][2] > numeric_limits<double>::min() ||
                    	 structure.box[1][0] > numeric_limits<double>::min() ||
                    	 structure.box[1][2] > numeric_limits<double>::min() ||
                    	 structure.box[2][0] > numeric_limits<double>::min() ||
                    	 structure.box[2][1] > numeric_limits<double>::min())
                	{
                    		 structure.isTriclinic = true;
                	}
                	 structure.calculateInverseBox();
                	 structure.calculateVolume();
            	}
	}
	else
	{
            	 structure.atoms.push_back(Atom());
            	 structure.atoms.back().index =  structure.numAtoms;
            	 structure.atoms.back().indexStructure =  structure.index;
            	 structure.atoms.back().tag            =  structure.numAtoms;
    		for(size_t c=0;c<3;c++)  structure.atoms.back().r[c]    =  mol->atoms[i].coordinates[c]*cflength;
            	 structure.atoms.back().element        =  structure.elementMap[mol->atoms[i].prop.symbol];// a changer
            	 structure.atoms.back().charge         = mol->atoms[i].charge;
    		for(size_t c=0;c<3;c++)  structure.atoms.back().fRef[c] = -mol->atoms[i].gradient[c]*cfdipole/cflength;
            	 structure.atoms.back().numNeighborsPerElement.resize(structure.numElements, 0);
            	 structure.numAtoms++;
            	 structure.numAtomsPerElement[elementMap[mol->atoms[i].prop.symbol]]++;// A changer
    		totalCharge += mol->atoms[i].charge;
        }
    }
    structure.energyRef = 0; 
    structure.charge = totalCharge;
    structure.chargeRef = totalCharge;
    for(size_t i=0;i<3;i++) structure.dipoleRef[i] = mol->dipole[i];
    for(size_t i=0;i<3;i++) structure.dipole[i] = mol->dipole[i];

    for (size_t i = 0; i < structure.numElements; i++)
    {
        if (structure.numAtomsPerElement[i] > 0)
        {
            structure.numElementsPresent++;
        }
    }

    if (structure.isPeriodic)
    {
        for (vector<Atom>::iterator it = structure.atoms.begin(); it != structure.atoms.end(); ++it)
        {
            structure.remap((*it));
        }
    }
}
#ifdef HIGH_DERIVATIVES
double InterfaceCChemIES::computed2Dipole(Molecule* mol, double*** dmu)
{
    double cd = 1.0/(AUTODEB);
    double cl = ANGTOBOHR;

    if(!initialized) initialize(directory,cl, cd, myRank);

    if(structure.atoms.size() != mol->nAtoms) setStructureFromMolecule(mol);
    else setCoordinatesFromMolecule(mol);

    structure.calculateNeighborList(maxCutoffRadius);
#ifdef NOSFGROUPS
    calculateSymmetryFunctions(structure, true);// true for derivative of dG/dxi
#else
    calculateSymmetryFunctionGroups(structure, true);
#endif
    calculateAtomicChargesNeuralNetworks(structure,true);// true for derivative dQ/dG
    calculateCharge(structure);
    calculateDipole(structure);
    size_t derivSize = 0;
    *dmu = Mode::computed2Dipole(structure, derivSize);
    for(size_t i=0;i<3;i++) mol->dipole[i] = structure.dipole[i] / cfdipole;
    size_t j =0;
    for (size_t i =  0; i < structure.atoms.size(); ++i,++j)
    {
	if(!strcmp(mol->atoms[j].prop.symbol,"TV") || !strcmp(mol->atoms[j].prop.symbol,"Tv"))
	{
		i--;
		continue;
	}
	mol->atoms[j].charge = structure.atoms[i].charge;
    }
    	double cv = convLength*convLength/convDipole;
    	for (size_t i =  0; i < derivSize; ++i) 
    		for (size_t k =  0; k < 3; ++k) 
			(*dmu)[k][i] *= cv;

    return structure.charge;
}
#ifdef HIGH_DERIVATIVES
static DerivativesIC* copyDerivatives(Derivatives& deriv, int size)
{
	DerivativesIC* derivS = new DerivativesIC;
	derivS->value = deriv();
	int order = deriv.getMaxOrder();
	if(order>=1)
	{
    		derivS->df = new double[size];
    		for (int i = 0; i < size; i++) 
    				derivS->df[i]  = deriv(i);
	}
	if(order>=2)
	{
		derivS->d2f = new double*[size];
		for (int i = 0; i < size; i++) 
			derivS->d2f[i]  = new double[i+1];

		for (int i = 0; i < size; i++) 
			for (int j = 0; j <=i; j++) 
				derivS->d2f[i][j]  = deriv(i,j);
	}
	if(order>=3)
	{
		derivS->d3f = new double**[size];
		for (int i = 0; i < size; i++) 
			derivS->d3f[i]  = new double*[i+1];

		for (int i = 0; i < size; i++) 
			for (int j = 0; j <=i; j++) 
				derivS->d3f[i][j]  = new double[j+1];

		for (int i = 0; i < size; i++) 
			for (int j = 0; j <=i; j++) 
				for (int k = 0; k <=j; k++) 
					derivS->d3f[i][j][k]  = deriv(i,j,k);
	}
	if(order>=4)
	{
		derivS->d4f = new double***[size];
		for (int i = 0; i < size; i++) 
			derivS->d4f[i]  = new double**[i+1];

		for (int i = 0; i < size; i++) 
			for (int j = 0; j <=i; j++) 
				derivS->d4f[i][j]  = new double*[j+1];

		for (int i = 0; i < size; i++) 
			for (int j = 0; j <=i; j++) 
				for (int k = 0; k <=j; k++) 
					derivS->d4f[i][j][k]  = new double[k+1];

		for (int i = 0; i < size; i++) 
			for (int j = 0; j <=i; j++) 
				for (int k = 0; k <=j; k++) 
					for (int l = 0; l <=k; l++) 
					derivS->d4f[i][j][k][l]  = deriv(i,j,k,l);
	}
	return derivS;
}
#endif
extern "C" {
	void *newInterfaceCChemIESDef()
	{
		return (void *)(new InterfaceCChemIES);
	}
	void *newInterfaceCChemIES(char* directory,
                                 double       cflength,
                                 double       cfdipole,
                                 int          myRank)
	{
		size_t l=strlen(directory)+strlen("/inputES.nn")+2;
		char* fn = (char*) malloc(l*sizeof(char));
		FILE* file;
		sprintf(fn,"%s%s",directory,"/inputES.nn");
    		file = fopen(fn,"r");
		free(fn);
		if(!file)
			return NULL;
		else
		{
			fclose(file);
			return (void *)(new InterfaceCChemIES(string(directory), cflength, cfdipole, myRank));
		}
	}
	void interfaceCChemIESInitialize(void* interface, char* directory,
                                 double       cflength,
                                 double       cfdipole,
                                 int          myRank)
	{
		if (interface != NULL)
		{
			InterfaceCChemIES *classPtr = static_cast<InterfaceCChemIES *>(interface);
			classPtr->initialize (string(directory), cflength, cfdipole, myRank);
		}
	}
	DerivativesIC** interfaceCChemIComputeDipoleHighDerivatives(void* interface, Molecule* mol, int order, int method)
	{
		bool err=false;
		DerivativesIC** derivS = NULL;
		if (interface != NULL)
		{
			InterfaceCChemIES *classPtr = static_cast<InterfaceCChemIES *>(interface);
			Derivatives* deriv = NULL;
			err = classPtr->computeDipoleHighDerivatives(mol, deriv, order, method);
			if(err!=0 && deriv)
			{
				derivS = new DerivativesIC* [3];
				for(int i=0;i<3;i++)
					derivS[i] = copyDerivatives(deriv[i], 3*mol->nAtoms);

				delete [] deriv; 
			}
		}
		return derivS;
	}
	int interfaceCChemIESComputed2Dipole(void* interface, Molecule* mol, double*** dmu, double* charge)
	{
		bool err=false;
		if (interface != NULL)
		{
			InterfaceCChemIES *classPtr = static_cast<InterfaceCChemIES *>(interface);
			double c= classPtr->computed2Dipole (mol,dmu);
			if(charge) *charge = c;
			err = true;
		}
		return (err?0:1);
	}
#endif
	int interfaceCChemIESGetChargeAndDipole(void* interface, int nAtoms, char** symbols, double** coordinates, double* charge, double* dipole)
	{
		bool err=false;
		if (interface != NULL)
		{
			InterfaceCChemIES *classPtr = static_cast<InterfaceCChemIES *>(interface);
			err = classPtr->computeChargeAndDipole (nAtoms, symbols, coordinates, charge, dipole);
		}
		return (err?0:1);
	}
	int interfaceCChemIESComputeChargeAndDipole(void* interface, Molecule* mol, double* charge)
	{
		bool err=false;
		if (interface != NULL)
		{
			InterfaceCChemIES *classPtr = static_cast<InterfaceCChemIES *>(interface);
			double c= classPtr->computeChargeAndDipole (mol);
			if(charge) *charge = c;
			err = true;
		}
		return (err?0:1);
	}
	int interfaceCChemIESComputedDipole(void* interface, Molecule* mol, double*** dmu, double* charge)
	{
		bool err=false;
		if (interface != NULL)
		{
			InterfaceCChemIES *classPtr = static_cast<InterfaceCChemIES *>(interface);
			double c= classPtr->computedDipole (mol,dmu);
			if(charge) *charge = c;
			err = true;
		}
		return (err?0:1);
	}
	void interfaceCChemIESDestroy(void *interface)
	{
		if (interface != NULL)
		{
			delete (static_cast<InterfaceCChemIES *>(interface));
		}
	}
	int runN2P2ES(Molecule* mol, char* directory, double cflength,double cfdipole, double* charge, int myRank)
	{
		void *interface = newInterfaceCChemIES(directory, cflength, cfdipole, myRank);
		int err = 0;
		err = interfaceCChemIESComputeChargeAndDipole(interface, mol, charge);
		interfaceCChemIESDestroy(interface);
		return err;
	}
}
