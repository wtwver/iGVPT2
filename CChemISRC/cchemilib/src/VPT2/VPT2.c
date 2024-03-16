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

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "../VPT2/VPT2Model.h"
#include "../VPT2/VPT2KModel.h"
#include "../Utils/Utils.h"
#include "../Utils/Constants.h"
boolean vpt2k(char* fileName)
{
	char* tmp = NULL;
	FILE* file = fopen(fileName, "rb");
	readOneString(file,"VPT2Model",&tmp);
	fclose(file);
	if(tmp && strstr(tmp,"VPT2+K")) return TRUE;
	return FALSE;
}
boolean gvpt2(char* fileName)
{
	char* tmp = NULL;
	FILE* file = fopen(fileName, "rb");
	readOneString(file,"VPT2Model",&tmp);
	fclose(file);
	if(tmp && strstr(tmp,"GVPT2")) return TRUE;
	return FALSE;
}
int vpt2(char* fileName)
{
	if(!fileName)
	{
		printf("Vous devez fournir le nom du fichier de donnees\n");
		return 1;
	}
	if(!vpt2k(fileName) && !gvpt2(fileName))
	{
		VPT2Model  vpt2Model = newVPT2Model();
		vpt2Model.klass->readData(&vpt2Model, fileName);
		vpt2Model.klass->computeAnharmonic(&vpt2Model);
	}
	else
	{
		VPT2KModel  vpt2KModel = newVPT2KModel();
		vpt2KModel.klass->readData(&vpt2KModel, fileName);
		vpt2KModel.klass->computeAnharmonic(&vpt2KModel);
	}

	return 0;
}
