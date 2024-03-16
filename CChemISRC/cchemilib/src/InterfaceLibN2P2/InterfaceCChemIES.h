/********************************************************************************
 cchemi is an interface to ab initio computational chemistry programs 
 designed for add them many functionalities non available in these packages.
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
********************************************************************************/

// writen by A.R. Allouche based on InterfaceLammps of n2p2 software

#ifndef INTERFACECCHEMIES_H
#define INTERFACECCHEMIES_H

#include "Mode.h"
#include "Structure.h"
#include <cstddef> // std::size_t
#include <string> 
#include "../Molecule/Molecule.h"

#ifdef HIGH_DERIVATIVES
#include "Derivatives.h"
#endif

namespace nnp
{

class InterfaceCChemIES : public Mode
{
public:
    InterfaceCChemIES();
    InterfaceCChemIES(std::string directory, double cflength, double cfdipole, int myRank);
    /** Initialize the CCHEMI interface.
     *
     * @param[in] directory Directory containing NNP data files (weights,
     *                      scaling, settings).
     * @param[in] cflength Length unit conversion factor.// CCheimI use kcal/mol for energy & Angstorom for distance
     * @param[in] cfdiople Dipole unit conversion factor.
     * @param[in] myRank MPI process rank (passed on to structure index).// not used presently
     */
    void   initialize(std::string directory,
                      double       cflength,
                      double       cfdipole,
                      int          myRank);
    /** Calculate symmetry functions, atomic neural networks and sum of local
     * energy contributions.
     */
    void   process();
    /** Return sum of local energy contributions.
     *
     * @return Sum of local energy contributions.
     */
    void readStructureFromFile(char* const& fileName);
#ifdef HIGH_DERIVATIVES
    bool computeDipoleHighDerivatives(Molecule* mol, Derivatives*& deriv, int order, int method);
#endif
    bool computeChargeAndDipole(int nAtoms, char** symbols, double** coordinates, double* charge, double* dipole);
    double computeChargeAndDipole(Molecule* mol);
    double computedDipole(Molecule* mol, double*** dmu);
    double computed2Dipole(Molecule* mol, double*** dmu);
    double run(Molecule* mol, std::string directory, double cflength,double cfdipole, int myRank);
    bool isInitialized() const;
    void setStructure(int nAtoms, char** symbols, double** coordinates);
    void setCoordinates(int nAtoms, char** symbols, double** coordinates);
    void setStructureFromMolecule(Molecule* mol);
    void setCoordinatesFromMolecule(Molecule* mol);

protected:
    std::string directory;
    /// Process rank.
    int       myRank;
    /// Initialization state.
    bool      initialized;
    /// Corresponds to CCHEMI `cflength` keyword.
    double    cflength;
    /// Corresponds to CCHEMI `cfenergy` keyword.
    double    cfdipole;
    /// Structure containing local atoms.
    Structure structure;
};

//////////////////////////////////
// Inlined function definitions //
//////////////////////////////////

inline bool InterfaceCChemIES::isInitialized() const
{
    return initialized;
}

}

#endif
