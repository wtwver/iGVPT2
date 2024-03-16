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

/* cchemi.c */

#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <time.h>
#include <locale.h>

#include <Utils/Types.h>
#include <Utils/Utils.h>
#include <JobControl/Job.h>

int main(int argc, char *argv[])
{
	char* inputFileName = NULL;
	Job job;
	if(argc<2)
	{
		printf("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n");
		printf("Usage : cchemi inputFileName.ici\n");
		printf("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n");
		exit(1);
	}
	inputFileName = argv[1];

/*
	g_setenv("LANG","en_US",TRUE);
	g_setenv("GDM_LANG","en_US",TRUE);
*/
	setlocale(LC_ALL,"C");
	setlocale(LC_NUMERIC,"C");

	initAll(argc, argv);
	readRessources();
	job = newJob(inputFileName);
	runJob(&job);
	finalize();

	return 0;
}
