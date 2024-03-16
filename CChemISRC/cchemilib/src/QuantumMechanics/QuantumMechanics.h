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

#ifndef __CCHEMILIB_QUANTUMMECHANICS_H__
#define __CCHEMILIB_QUANTUMMECHANICS_H__

#include "../Molecule/Molecule.h"
#include "../QuantumMechanics/QuantumMechanicsModel.h"

QuantumMechanicsModel createMopacModel(Molecule* mol, char* keywords,char* dirName,  char* nameCommand, Constraints constraints,FILE* logfile);
QuantumMechanicsModel createFireFlyModel(Molecule* mol, char* keywords, char* dirName,  char* nameCommand, Constraints constraints,FILE* logfile);
QuantumMechanicsModel createOrcaModel(Molecule* mol, char* keywords, char* dirName,  char* nameCommand, Constraints constraints,FILE* logfile);
QuantumMechanicsModel createGaussianModel(Molecule* mol, char* keywords, char* dirName,  char* nameCommand, Constraints constraints,FILE* logfile);
QuantumMechanicsModel createGenericModel(Molecule* mol, char* keywords, char* dirName,  char* nameCommand, Constraints constraints,FILE* logfile);
QuantumMechanicsModel createOpenBabelModel(Molecule* mol, char* keywords, char* dirName,  char* nameCommand, Constraints constraints,FILE* logfile);
QuantumMechanicsModel createN2P2Model(Molecule* mol, char* N2P2Dir, Constraints constraints,FILE* logfile);
QuantumMechanicsModel createTMModel(Molecule* mol, char* tmModule, Constraints constraints,FILE* logfile);
void setH4Correction(FILE* file,QuantumMechanicsModel* qmModel);
void setSRBCorrection(FILE* file,QuantumMechanicsModel* qmModel);

#endif /* __CCHEMILIB_QUANTUMMECHANICS_H__ */

