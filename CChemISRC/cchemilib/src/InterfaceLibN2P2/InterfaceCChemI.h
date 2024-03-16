// n2p2 - A neural network potential package
// Copyright (C) 2018 Andreas Singraber (University of Vienna)
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <https://www.gnu.org/licenses/>.

// writen by A.R. Allouche based on InterfaceLammps

#ifndef INTERFACECCHEMI_H
#define INTERFACECCHEMI_H

#include "Mode.h"
#include "Structure.h"
#include <cstddef> // std::size_t
#include <string> 
#include "../Molecule/Molecule.h"

namespace nnp
{

class InterfaceCChemI : public Mode
{
public:
    InterfaceCChemI();
    InterfaceCChemI(std::string directory, double cflength, double cfenergy, int myRank);
    /** Initialize the CCHEMI interface.
     *
     * @param[in] directory Directory containing NNP data files (weights,
     *                      scaling, settings).
     * @param[in] cflength Length unit conversion factor.// CCheimI use kcal/mol for energy & Angstorom for distance
     * @param[in] cfenergy Energy unit conversion factor.
     * @param[in] myRank MPI process rank (passed on to structure index).// not used presently
     */
    void   initialize(std::string directory,
                      double       cflength,
                      double       cfenergy,
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
    bool computeForces(int nAtoms, char** symbols, double** coordinates, double** forces, double* energy);
    bool computeEnergy(int nAtoms, char** symbols, double** coordinates, double* energy);
    bool computeHessian(Molecule* mol, double** hessian);
    bool computeHighDerivatives(Molecule* mol, Derivatives& deriv, int order, int method);
    bool computeGradients(Molecule* mol);
    bool computeEnergy(Molecule* mol);
    bool run(Molecule* mol, std::string directory, double cflength,double cfenergy, int myRank);
    bool   isInitialized() const;
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
    double    cfenergy;
    /// Structure containing local atoms.
    Structure structure;
};

//////////////////////////////////
// Inlined function definitions //
//////////////////////////////////

inline bool InterfaceCChemI::isInitialized() const
{
    return initialized;
}

}

#endif
