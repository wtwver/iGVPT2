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

// writen by A.R. Allouche based on InterfaceCChemI of n2p2 software

//#ifndef NOMPI
//#include <mpi.h>
//#include "mpi-extra.h"
//#endif
#include "InterfaceCChemI.h"
#include "Atom.h"
#include "Element.h"
#include "utility.h"
#include "utilityDerivatives.h"
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

InterfaceCChemI::InterfaceCChemI() : myRank      (0    ),
                                     initialized (false),
                                     cflength    (1.0  ),
                                     cfenergy    (1.0  )
{
	directory = string(getenv("PWD"));
}
InterfaceCChemI::InterfaceCChemI(string directory,
                                 double       cflength,
                                 double       cfenergy,
                                 int          myRank)
{
	initialize(directory,cflength, cfenergy, myRank);
}

void InterfaceCChemI::initialize(string directory,
                                 double       cflength,
                                 double       cfenergy,
                                 int          myRank)
{
    this->directory = directory;
    this->cflength = cflength;
    this->cfenergy = cfenergy;
    this->myRank = myRank;
    log.writeToStdout = false;
    string dir(directory);
    Mode::initialize();
    loadSettingsFile(dir + "/input.nn");
    setupGeneric();
    setupSymmetryFunctionScaling(dir + "/scaling.data");
    setupSymmetryFunctionStatistics(false, false, false, false);
    setupNeuralNetworkWeights(dir + "/weights.%03d.data");

    structure.setElementMap(elementMap);

    initialized = true;
}
void InterfaceCChemI::readStructureFromFile(char* const& fileName)
{
    ifstream file;
    file.open(fileName);
    structure.reset();
    structure.setElementMap(elementMap);
    structure.readFromFile(string(fileName));
    removeEnergyOffset(structure);
    if (normalize)
    {
        double convDipole;
        structure.toNormalizedUnits(meanEnergy, convEnergy, convLength, convDipole);
    }
    file.close();

    return;
}


void InterfaceCChemI::process()
{
#ifdef NOSFGROUPS
    calculateSymmetryFunctions(structure, true);
#else
    calculateSymmetryFunctionGroups(structure, true);
#endif
    calculateAtomicNeuralNetworks(structure, true);
    calculateEnergy(structure);
    if (normalize)
    {
        structure.energy = physicalEnergy(structure, false);
    }
    addEnergyOffset(structure, false);

    return;
}

bool InterfaceCChemI::computeForces(int nAtoms, char** symbols, double** coordinates, double** forces, double* energy)
{
    double ce = 1.0/(AUTOKCAL);
    double cl = ANGTOBOHR;
    if(!initialized) initialize(directory,cl, ce, myRank);

    if(structure.atoms.size() != nAtoms) setStructure(nAtoms, symbols, coordinates);
    else setCoordinates(nAtoms, symbols, coordinates);
    if (normalize) structure.toNormalizedUnits(meanEnergy, convEnergy, convLength,1.0);
    double const cfforce = cflength / cfenergy;

    structure.calculateNeighborList(maxCutoffRadius);
#ifdef NOSFGROUPS
    calculateSymmetryFunctions(structure, true);
#else
    calculateSymmetryFunctionGroups(structure, true);
#endif
    calculateAtomicNeuralNetworks(structure, true);
    calculateEnergy(structure);
    calculateForces(structure);
    if (normalize)
    {
        structure.toPhysicalUnits(meanEnergy, convEnergy, convLength,1.0);
    }
    addEnergyOffset(structure, false);
    //addEnergyOffset(structure, true);
    Atom const* a = NULL;
    for (size_t i =  0; i < structure.atoms.size(); ++i)
    {
        a = &(structure.atoms.at(i));
	for(size_t c=0; c<3;c++) forces[c][i] = a->f[c]*cfforce;
    }
    *energy = structure.energy / cfenergy;
    return true;
}
bool InterfaceCChemI::computeHessian(Molecule* mol, double** hessian)
{
//	cerr<<"Begin Grad"<<endl;
    double ce = 1.0/(AUTOKCAL);
    double cl = ANGTOBOHR;
    if(!initialized) initialize(directory,cl, ce, myRank);

//	cerr<<"Begin setSC"<<endl;

    if(structure.atoms.size() != mol->nAtoms) setStructureFromMolecule(mol);
    else setCoordinatesFromMolecule(mol);
// DEBUG setStructureFromMolecule(mol);

   if (normalize) structure.toNormalizedUnits(meanEnergy, convEnergy, convLength,1.0);
   double const cfforce = cflength / cfenergy;
   double const cfhessian = cflength*cflength/ cfenergy;
	//cerr<<"Begin calculateNeighborList"<<endl;

    structure.calculateNeighborList(maxCutoffRadius);
	//cerr<<"End calculateNeighborList"<<endl;
#ifdef NOSFGROUPS
	//cerr<<"Begin calculateSymmetryFunctions("<<endl;
    calculateSymmetryFunctions(structure, true);
#else
	//cerr<<"Begin calculateSymmetryFunctionsGroup("<<endl;
    calculateSymmetryFunctionGroups(structure, true);
#endif
	//cerr<<"Begin nnp"<<endl;
    calculateAtomicNeuralNetworks(structure, true);
    calculateEnergy(structure);
    calculateForces(structure);
    size_t derivSize;
    double* deriv = Mode::computeHessian(structure, derivSize);
    if (normalize)
    {
        toPhysicalUnits(deriv, derivSize, convEnergy, convLength);
        structure.toPhysicalUnits(meanEnergy, convEnergy, convLength, 1.0);
    }

    addEnergyOffset(structure, false);
    //addEnergyOffset(structure, true);
    Atom const* a = NULL;
    int ia=0;
	//cerr<<"Begin copyeF"<<endl;
    for(size_t i=0;i<mol->nAtoms;i++)
    {
	if(!strcmp(mol->atoms[i].prop.symbol,"TV") || !strcmp(mol->atoms[i].prop.symbol,"Tv")) continue;
	{
        	// Set pointer to atom.
        	a = &(structure.atoms.at(ia));
    		for(size_t c=0;c<3;c++) mol->atoms[i].gradient[c] = -a->f[c]*cfforce;
		ia++;
        }
    }
    for(size_t i=0;i<derivSize;i++) deriv[i] *= cfhessian;
    mol->potentialEnergy = structure.energy / cfenergy;
    *hessian = deriv;
    return true;
}
#ifdef HIGH_DERIVATIVES
bool InterfaceCChemI::computeHighDerivatives(Molecule* mol, Derivatives& deriv, int order, int method)
{
//	cerr<<"Begin Grad"<<endl;
    double ce = 1.0/(AUTOKCAL);
    double cl = ANGTOBOHR;
    // Conversion in atomic unit
    if(!initialized) initialize(directory,cl, ce, myRank);


//	cerr<<"Begin setSC"<<endl;

    if(structure.atoms.size() != mol->nAtoms) setStructureFromMolecule(mol);
    else setCoordinatesFromMolecule(mol);
// DEBUG setStructureFromMolecule(mol);

   if (normalize) structure.toNormalizedUnits(meanEnergy, convEnergy, convLength,1.0);
	//cerr<<"Begin calculateNeighborList"<<endl;
    //cerr<<"After Norm convEnergy="<<convEnergy<<endl;
    //cerr<<"convLength="<<convLength<<endl;

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
    calculateAtomicNeuralNetworks(structure, true);
    calculateEnergy(structure);
    calculateForces(structure);
    deriv = Mode::computeHighDerivatives(structure, order, method);
    if (normalize)
    {
	deriv.toPhysicalUnits(convEnergy, convLength);
        structure.toPhysicalUnits(meanEnergy, convEnergy, convLength, 1.0);
    }
    addEnergyOffset(structure, false);
    addEnergyOffset(structure, true);
    //cerr<<"After Phys convEnergy="<<convEnergy<<endl;
    //cerr<<"convLength="<<convLength<<endl;
    //cerr<<"cflength="<<cflength<<endl;
    //cerr<<"cfenergy="<<cfenergy<<endl;
    // conversion in cchemilib unit : kcal, Ang
    deriv.toPhysicalUnits(cfenergy, cflength);
    return true;
}
#endif
bool InterfaceCChemI::computeEnergy(int nAtoms, char** symbols, double** coordinates, double* energy)
{
    double ce = 1.0/(AUTOKCAL);
    double cl = ANGTOBOHR;
    if(!initialized) initialize(directory,cl, ce, myRank);
    if(structure.atoms.size() != nAtoms) setStructure(nAtoms, symbols, coordinates);
    else setCoordinates(nAtoms, symbols, coordinates);
    if (normalize) structure.toNormalizedUnits(meanEnergy, convEnergy, convLength,1.0);

    structure.calculateNeighborList(maxCutoffRadius);
#ifdef NOSFGROUPS
    calculateSymmetryFunctions(structure, false);
#else
    calculateSymmetryFunctionGroups(structure, false);
#endif
    calculateAtomicNeuralNetworks(structure, false);
    calculateEnergy(structure);
    if (normalize)
    {
        structure.toPhysicalUnits(meanEnergy, convEnergy, convLength,1.0);
    }
    addEnergyOffset(structure, false);
    //addEnergyOffset(structure, true);
    *energy = structure.energy / cfenergy;
    return true;
}
bool InterfaceCChemI::computeGradients(Molecule* mol)
{
//	cerr<<"Begin Grad"<<endl;
    double ce = 1.0/(AUTOKCAL);
    double cl = ANGTOBOHR;
    if(!initialized) initialize(directory,cl, ce, myRank);

//	cerr<<"Begin setSC"<<endl;

    if(structure.atoms.size() != mol->nAtoms) setStructureFromMolecule(mol);
    else setCoordinatesFromMolecule(mol);
// DEBUG setStructureFromMolecule(mol);

   if (normalize) structure.toNormalizedUnits(meanEnergy, convEnergy, convLength,1.0);
   double const cfforce = cflength / cfenergy;
	//cerr<<"Begin calculateNeighborList"<<endl;

    structure.calculateNeighborList(maxCutoffRadius);
	//cerr<<"End calculateNeighborList"<<endl;
#ifdef NOSFGROUPS
	//cerr<<"Begin calculateSymmetryFunctions("<<endl;
    calculateSymmetryFunctions(structure, true);
#else
	//cerr<<"Begin calculateSymmetryFunctionsGroup("<<endl;
    calculateSymmetryFunctionGroups(structure, true);
#endif
	//cerr<<"Begin nnp"<<endl;
    calculateAtomicNeuralNetworks(structure, true);
    calculateEnergy(structure);
    calculateForces(structure);
    if (normalize) structure.toPhysicalUnits(meanEnergy, convEnergy, convLength,1.0);
    addEnergyOffset(structure, false);
    //addEnergyOffset(structure, true);
    Atom const* a = NULL;
    int ia=0;
	//cerr<<"Begin copyeF"<<endl;
    for(size_t i=0;i<mol->nAtoms;i++)
    {
	if(!strcmp(mol->atoms[i].prop.symbol,"TV") || !strcmp(mol->atoms[i].prop.symbol,"Tv")) continue;
	{
        	// Set pointer to atom.
        	a = &(structure.atoms.at(ia));
    		for(size_t c=0;c<3;c++) mol->atoms[i].gradient[c] = -a->f[c]*cfforce;
		ia++;
        }
    }
    mol->potentialEnergy = structure.energy / cfenergy;
    return true;
}
bool InterfaceCChemI::computeEnergy(Molecule* mol)
{
    double ce = 1.0/(AUTOKCAL);
    double cl = ANGTOBOHR;
    if(!initialized) initialize(directory,cl, ce, myRank);

    if(structure.atoms.size() != mol->nAtoms) setStructureFromMolecule(mol);
    else setCoordinatesFromMolecule(mol);
   if (normalize) structure.toNormalizedUnits(meanEnergy, convEnergy, convLength,1.0);

    structure.calculateNeighborList(maxCutoffRadius);
#ifdef NOSFGROUPS
    calculateSymmetryFunctions(structure, false);
#else
    calculateSymmetryFunctionGroups(structure, false);
#endif
    calculateAtomicNeuralNetworks(structure, false);
    calculateEnergy(structure);
    if (normalize) structure.toPhysicalUnits(meanEnergy, convEnergy, convLength,1.0);
    addEnergyOffset(structure, false);
    //addEnergyOffset(structure, true);
    mol->potentialEnergy = structure.energy / cfenergy;
    return true;
}
bool InterfaceCChemI::run(Molecule* mol, string directory, double cflength,double cfenergy, int myRank)
{
    initialize(directory, cflength, cfenergy,myRank);
    for (size_t i =  0; i < structure.atoms.size(); ++i)
    {
	for(size_t c=0;c<3;c++) structure.atoms[i].r[c] = mol->atoms[i].coordinates[c]*cflength;
    }
    if (normalize) structure.toNormalizedUnits(meanEnergy, convEnergy, convLength,1.0);
    structure.calculateNeighborList(maxCutoffRadius);
#ifdef NOSFGROUPS
    calculateSymmetryFunctions(structure, true);
#else
    calculateSymmetryFunctionGroups(structure, true);
#endif
    calculateAtomicNeuralNetworks(structure, true);
    calculateEnergy(structure);
    if (normalize)
    {
        structure.toPhysicalUnits(meanEnergy, convEnergy, convLength,1.0);
    }
    addEnergyOffset(structure, false);
    addEnergyOffset(structure, true);
    mol->potentialEnergy = structure.energy / cfenergy;
    return true;
}
void InterfaceCChemI::setStructure(int nAtoms, char** symbols, double** coordinates)
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
    removeEnergyOffset(structure);
}
void InterfaceCChemI::setCoordinates(int nAtoms, char** symbols, double** coordinates)
{
    int iBoxVector = 0;
    double totalCharge = 0;
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
    		ai->charge = 0;
            	 ai->numNeighborsPerElement.resize(structure.numElements, 0);
    		for(size_t c=0;c<3;c++)  ai->fRef[c] = 0.0;
		ia++;
    		totalCharge += ai->charge;
        }
    }
    structure.energyRef = 0;
    structure.chargeRef = totalCharge;
    removeEnergyOffset(structure);
}
void InterfaceCChemI::setCoordinatesFromMolecule(Molecule* mol)
{
    int iBoxVector = 0;
    double totalCharge = 0;
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
    		totalCharge += mol->atoms[i].charge;
        }
    }
    structure.energyRef = mol->potentialEnergy*cfenergy;
    structure.chargeRef = totalCharge;
    removeEnergyOffset(structure);
}
void InterfaceCChemI::setStructureFromMolecule(Molecule* mol)
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
    		for(size_t c=0;c<3;c++)  structure.atoms.back().fRef[c] = -mol->atoms[i].gradient[c]*cfenergy/cflength;
            	 structure.atoms.back().numNeighborsPerElement.resize(structure.numElements, 0);
            	 structure.numAtoms++;
            	 structure.numAtomsPerElement[elementMap[mol->atoms[i].prop.symbol]]++;// A changer
    		totalCharge += mol->atoms[i].charge;
        }
    }
    structure.energyRef = mol->potentialEnergy*cfenergy; 
    structure.chargeRef = totalCharge;

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
    removeEnergyOffset(structure);
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
	void *newInterfaceCChemIDef()
	{
		return (void *)(new InterfaceCChemI);
	}
	void *newInterfaceCChemI(char* directory,
                                 double       cflength,
                                 double       cfenergy,
                                 int          myRank)
	{
		return (void *)(new InterfaceCChemI(string(directory), cflength, cfenergy, myRank));
	}
	void interfaceCChemIInitialize(void* interface, char* directory,
                                 double       cflength,
                                 double       cfenergy,
                                 int          myRank)
	{
		if (interface != NULL)
		{
			InterfaceCChemI *classPtr = static_cast<InterfaceCChemI *>(interface);
			classPtr->initialize (string(directory), cflength, cfenergy, myRank);
		}
	}
	int interfaceCChemIComputeForces(void* interface, int nAtoms, char** symbols, double** coordinates, double** forces, double* energy)
	{
		bool err=false;
		if (interface != NULL)
		{
			InterfaceCChemI *classPtr = static_cast<InterfaceCChemI *>(interface);
			err = classPtr->computeForces (nAtoms, symbols, coordinates, forces, energy);
		}
		return (err?0:1);
	}
	int interfaceCChemIGetEnergy(void* interface, int nAtoms, char** symbols, double** coordinates, double* energy)
	{
		bool err=false;
		if (interface != NULL)
		{
			InterfaceCChemI *classPtr = static_cast<InterfaceCChemI *>(interface);
			err = classPtr->computeEnergy (nAtoms, symbols, coordinates, energy);
		}
		return (err?0:1);
	}
	int interfaceCChemIComputeHessian(void* interface, Molecule* mol, double** hessian)
	{
		bool err=false;
		if (interface != NULL)
		{
			InterfaceCChemI *classPtr = static_cast<InterfaceCChemI *>(interface);
			err = classPtr->computeHessian(mol, hessian);
		}
		return (err?0:1);
	}
#ifdef HIGH_DERIVATIVES
	DerivativesIC* interfaceCChemIComputeHighDerivatives(void* interface, Molecule* mol, int order, int method)
	{
		bool err=false;
		DerivativesIC* derivS = NULL;
		if (interface != NULL)
		{
			InterfaceCChemI *classPtr = static_cast<InterfaceCChemI *>(interface);
			Derivatives deriv;
			err = classPtr->computeHighDerivatives(mol, deriv, order, method);
			if(err!=0)
				derivS = copyDerivatives(deriv, 3*mol->nAtoms);
		}
		return derivS;
	}
#endif
	int interfaceCChemIComputeGradients(void* interface, Molecule* mol)
	{
		bool err=false;
		if (interface != NULL)
		{
			InterfaceCChemI *classPtr = static_cast<InterfaceCChemI *>(interface);
			err = classPtr->computeGradients(mol);
		}
		return (err?0:1);
	}
	int interfaceCChemIComputeEnergy(void* interface, Molecule* mol)
	{
		bool err=false;
		if (interface != NULL)
		{
			InterfaceCChemI *classPtr = static_cast<InterfaceCChemI *>(interface);
			err = classPtr->computeEnergy (mol);
		}
		return (err?0:1);
	}
	void interfaceCChemIDestroy(void *interface)
	{
		if (interface != NULL)
		{
			delete (static_cast<InterfaceCChemI *>(interface));
		}
	}
	int runN2P2(Molecule* mol, char* directory, double cflength,double cfenergy, int myRank, int computeGradients)
	{
		void *interface = newInterfaceCChemI(directory, cflength, cfenergy, myRank);
		int err = 0;
		if(computeGradients) interfaceCChemIComputeGradients(interface, mol);
		else err = interfaceCChemIComputeEnergy(interface, mol);
		interfaceCChemIDestroy(interface);
		return err;
	}
}
