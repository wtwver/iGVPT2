
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


#ifndef __CCHEMILIB_STEEPESTDESCENT_H__
#define __CCHEMILIB_STEEPESTDESCENT_H__

#include "../MolecularMechanics/ForceField.h"


typedef struct _SteepestDescent  SteepestDescent;

struct _SteepestDescent
{
	ForceField* forceField;
	int numberOfAtoms;
	int updateFrequency;
	int maxIterations;
	int updateNumber;
	double epsilon;
	double initialStep;
	double rmsDeplacment;
	double maxDeplacment;
	double gradientNorm;
	int maxLines;
	FILE* logfile;

};
void	freeSteepestDescent(SteepestDescent* steepestDescent);
void	runSteepestDescent(SteepestDescent* steepestDescent, ForceField* forceField, int updateFrequency, int maxIterations, double epsilon,int maxLines);

#endif /* __CCHEMILIB_STEEPESTDESCENT_H__ */

