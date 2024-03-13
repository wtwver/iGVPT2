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
igvpt2.cpp - calculate GVPT2
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
#include <VPT2/VPT2.h>
#include <Utils/Constants.h>
}
#endif

#include "JobiGVPT2.h"

int main(int argc, char *argv[])
{
	int r = 0;
	char* inputFileName = NULL;
	JobiGVPT2 job;

	if(argc<2)
	{
		printf("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n");
		printf("Usage : iGVPT2 inputFileName.ici\n");
		printf("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n");
		exit(1);
	}
	inputFileName = argv[1];
	setlocale(LC_ALL,"C");
	setlocale(LC_NUMERIC,"C");

	initAll(argc, argv);
	readRessources();
	job = newJobiGVPT2(inputFileName);
	runJobiGVPT2(&job);
	finalize();
	return r;
}
