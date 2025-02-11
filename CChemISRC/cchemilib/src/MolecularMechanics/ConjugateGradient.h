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

#ifndef __CCHEMILIB_CONJUGATEGRADIENT_H__
#define __CCHEMILIB_CONJUGATEGRADIENT_H__

#include "../MolecularMechanics/ForceField.h"
typedef struct _ConjugateGradient  ConjugateGradient;
typedef struct 
{
	double gradientNorm;
	int maxIterations;
	int updateFrequency;
	int maxLines;
	double initialStep;
	int method;
}ConjugateGradientOptions;

struct _ConjugateGradient
{
	ForceField* forceField;
	int numberOfAtoms;
	int updateFrequency;
	int maxIterations;
	int updateNumber;
	int maxLine;
	double epsilon;
	double initialStep;
	double rmsDeplacment;
	double maxDeplacment;
	double gradientNorm;
	double lastInitialBracket;
	double initialBracket;
	double term;
	FILE* logfile;

	double* lastGradient[3];
	double* direction[3];
	Molecule* temporaryMolecule;
};
void	freeConjugateGradient(ConjugateGradient* conjugateGradient);
void	runConjugateGradient(ConjugateGradient* conjugateGradient, ForceField* forceField, ConjugateGradientOptions conjugateGradientOptions);
void setCGOptions(FILE* file, ConjugateGradientOptions* conjugateGradientOptions);

#endif /* __CCHEMILIB_CONJUGATEGRADIENT_H__ */

