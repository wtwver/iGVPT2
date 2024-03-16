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

/* InterfacTM.c */


#ifndef OS_WIN32
#include <unistd.h>
#endif

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#ifdef ENABLE_MPI
#include <mpi.h>
#endif

#include "../Utils/Utils.h"
#include "../Utils/AtomsProp.h"
#include "../Utils/Constants.h"
#include "../Utils/Types.h"
#include "../Molecule/Molecule.h"
#include "InterfaceTM.h"

#ifdef ENABLE_PYTHON

PyObject* vectorToList_Str(char** data, int size)
{
	int i;
	PyObject* listObj = PyList_New(size);
	if (!listObj) fprintf(stderr," Unable to allocate memory for Python list");
	for (i = 0; i < size; i++)
	{
		PyObject *num = PyUnicode_FromString(data[i]);
		if (!num) {
			Py_DECREF(listObj);
			fprintf(stderr," Unable to allocate memory for Python list");
		}
		PyList_SET_ITEM(listObj, i, num);
	}
	return listObj;
}

PyObject* vectorToList_Float(double* data, int size)
{
	int i;
	PyObject* listObj = PyList_New(size);
	if (!listObj) fprintf(stderr," Unable to allocate memory for Python list");
	for (i = 0; i < size; i++)
	{
		PyObject *num = PyFloat_FromDouble(data[i]);
		if (!num) {
			Py_DECREF(listObj);
			fprintf(stderr," Unable to allocate memory for Python list");
		}
		PyList_SET_ITEM(listObj, i, num);
	}
	return listObj;
}
// PyObject -> Vector
double* listToVector_Float(PyObject* incoming)
{
	double* data = NULL;
	if (PyList_Check(incoming))
	{
		int i;
		Py_ssize_t j = 0;
		int size = PyList_Size(incoming);
		data = malloc(size*sizeof(double));
		for(j = 0, i=0; j < PyList_Size(incoming); j++,i++)
				data[i] = PyFloat_AsDouble(PyList_GetItem(incoming, j));
	}
	else fprintf(stderr," Passed PyObject pointer was not a list or tuple!");

	return data;
}
static int initManager(InterfaceTM* interfaceTM, Molecule* mol)
{
	int i;

	/* fprintf(stderr,"initManager ================== \n");*/
	if (interfaceTM->pSetData && PyCallable_Check(interfaceTM->pSetData))
	{
		PyObject *pname = PyUnicode_FromString("Test");
		PyObject *pnatoms = PyLong_FromLongLong((long long) mol->nAtoms);
		PyObject *pmult = PyLong_FromLongLong((long long) mol->spinMultiplicity);
		PyObject *pcharge = PyLong_FromLongLong((long long) mol->totalCharge);
		double *X = malloc(mol->nAtoms*sizeof(double));
		double *Y = malloc(mol->nAtoms*sizeof(double));
		double *Z = malloc(mol->nAtoms*sizeof(double));
		char** symbols = malloc(mol->nAtoms*sizeof(char*));
		PyObject *pArgs = PyTuple_New(8);
		PyObject* pX = NULL;
		PyObject* pY = NULL;
		PyObject* pZ = NULL;
		PyObject* pSymbols = NULL;

		for(i=0;i<mol->nAtoms;i++) symbols[i] = strdup(mol->atoms[i].prop.symbol);
		pSymbols = vectorToList_Str(symbols, mol->nAtoms);

		for(i=0;i<mol->nAtoms;i++) X[i] = mol->atoms[i].coordinates[0];
		pX = vectorToList_Float(X, mol->nAtoms);
		for(i=0;i<mol->nAtoms;i++) Y[i] = mol->atoms[i].coordinates[1];
		pY = vectorToList_Float(Y, mol->nAtoms);
		for(i=0;i<mol->nAtoms;i++) Z[i] = mol->atoms[i].coordinates[2];
		pZ = vectorToList_Float(Z, mol->nAtoms);


                PyTuple_SetItem(pArgs, 0, pname);
                PyTuple_SetItem(pArgs, 1, pnatoms);
                PyTuple_SetItem(pArgs, 2, pcharge);
                PyTuple_SetItem(pArgs, 3, pmult);
                PyTuple_SetItem(pArgs, 4, pSymbols);
                PyTuple_SetItem(pArgs, 5, pX);
                PyTuple_SetItem(pArgs, 6, pY);
                PyTuple_SetItem(pArgs, 7, pZ);
		interfaceTM->pa = PyObject_CallObject(interfaceTM->pSetData, pArgs);

		if(interfaceTM->pa && interfaceTM->pGetManager && PyCallable_Check(interfaceTM->pGetManager))
		{
			PyObject *pArgs = PyTuple_New(1);
                	PyTuple_SetItem(pArgs, 0, interfaceTM->pa);
			/* fprintf(stderr,"Call getManager\n");*/
			interfaceTM->pManager = PyObject_CallObject(interfaceTM->pGetManager, pArgs);
			/* if(interfaceTM->pManager) fprintf(stderr,"pManager != NULL\n");*/
			if(!interfaceTM->pManager) 
			{
				fprintf(stderr,"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n");
				fprintf(stderr,"I cannot get Manager from tmModule\n");
				fprintf(stderr,"Program stopped\n");
				fprintf(stderr,"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n");
			}
			//if(pArgs) Py_DECREF(pArgs);
		}

		if(pSymbols) Py_DECREF(pSymbols);
		if(pX) Py_DECREF(pX);
		if(pY) Py_DECREF(pY);
		if(pZ) Py_DECREF(pZ);
		for(i=0;i<mol->nAtoms;i++) if(symbols[i]) free(symbols[i]);
		free(symbols);
		free(X);
		free(Y);
		free(Z);
        }

	interfaceTM->initialized = 1;
	/* fprintf(stderr,"END initManager ================== \n");*/
	return 0;
}
static int setCoordinates(InterfaceTM* interfaceTM, Molecule* mol)
{

	if(!interfaceTM->pa) initManager(interfaceTM, mol);
	if(interfaceTM->pa && interfaceTM->pManager) 
	{
               	/* fprintf(stderr,"ok pa & pManage in setCoordinates\n");*/
		if (interfaceTM->pSetCoordinates && PyCallable_Check(interfaceTM->pSetCoordinates))
		{
			int i;
			double *X = malloc(mol->nAtoms*sizeof(double));
			double *Y = malloc(mol->nAtoms*sizeof(double));
			double *Z = malloc(mol->nAtoms*sizeof(double));
			PyObject* pX = NULL;
			PyObject* pY = NULL;
			PyObject* pZ = NULL;
			PyObject* a= NULL;

			for(i=0;i<mol->nAtoms;i++) X[i] = mol->atoms[i].coordinates[0];
			pX = vectorToList_Float(X, mol->nAtoms);
			for(i=0;i<mol->nAtoms;i++) Y[i] = mol->atoms[i].coordinates[1];
			pY = vectorToList_Float(Y, mol->nAtoms);
			for(i=0;i<mol->nAtoms;i++) Z[i] = mol->atoms[i].coordinates[2];
			pZ = vectorToList_Float(Z, mol->nAtoms);

                	/* fprintf(stderr,"ok pGetEnergy in interfaceTMComputeEnergy\n");*/
			a= PyObject_CallFunctionObjArgs(interfaceTM->pSetCoordinates, interfaceTM->pa,  pX, pY, pZ,  NULL);
			if(pX) Py_DECREF(pX);
			if(pY) Py_DECREF(pY);
			if(pZ) Py_DECREF(pZ);
            		if (a != NULL)
			{
                		/* fprintf(stderr,"Ok setCoordinates\n");*/
				return 0;
			}
			else
			{
                		fprintf(stderr,"a = NULL in setCoordinates\n");
				return 1;
			}
		}
	}
	return 0;
}


/*
static void printDebug(int nAtoms, char** symbols, double** coordinates, double** forces, double* energy)
{
	int i;
	fprintf(stderr," printDebug\n");
	fprintf(stderr," nAtoms = %d\n", nAtoms);
	fprintf(stderr," Energy = %f\n", *energy);
	for(i=0;i<nAtoms;i++)
	{
		fprintf(stderr," %s %f %f %f %f %f %f\n",symbols[i],
		coordinates[i][0], coordinates[i][1], coordinates[i][2], 
		forces[i][0], forces[i][1], forces[i][2]); 
	}
}
static void printDebug2(Molecule* mol)
{
	int i;
	fprintf(stderr," printDebug\n");
	fprintf(stderr," nAtoms = %d\n", mol->nAtoms);
	fprintf(stderr," Energy = %f\n", mol->potentialEnergy);
	for(i=0;i<mol->nAtoms;i++)
	{
		fprintf(stderr," %s %f %f %f %f %f %f\n",mol->atoms[i].prop.symbol,
		mol->atoms[i].coordinates[0],
		mol->atoms[i].coordinates[1],
		mol->atoms[i].coordinates[2],
		mol->atoms[i].gradient[0],
		mol->atoms[i].gradient[1],
		mol->atoms[i].gradient[2]);
	}
}
*/
static void interfaceTMInitialize(InterfaceTM* interfaceTM, char* tmModule)
{
	if(!Py_IsInitialized()) Py_Initialize();
	/*
	fprintf(stderr," End Py_Initialize\n");
	fprintf(stderr,"tmModule =%s\n",tmModule);
	*/
	interfaceTM->pName = PyUnicode_DecodeFSDefault(tmModule);
	if(interfaceTM->pName) 
	{
		/* fprintf(stderr," ok pName\n");*/
		interfaceTM->pModule = PyImport_Import(interfaceTM->pName);
		Py_DECREF(interfaceTM->pName);
		interfaceTM->pName = NULL;
	}
	if (interfaceTM->pModule)
	{
		/* fprintf(stderr," ok pModule\n");*/
        	interfaceTM->pGetManager = PyObject_GetAttrString(interfaceTM->pModule, "getManager");
        	/* pGet is a new reference */
		/* if(interfaceTM->pGetManager) fprintf(stderr," ok getManager\n");*/
		if(!interfaceTM->pGetManager) fprintf(stderr," problem getManager\n");
	}
	if (interfaceTM->pGetManager)
	{
        	interfaceTM->pGetEnergy = PyObject_GetAttrString(interfaceTM->pModule, "getEnergy");
		/* if(interfaceTM->pGetEnergy) fprintf(stderr,"ok getEnergy\n");*/
		if(!interfaceTM->pGetEnergy) fprintf(stderr," problem getEnergy\n");
	}
	if (interfaceTM->pGetEnergy)
	{
        	interfaceTM->pGetEnergyAndForces = PyObject_GetAttrString(interfaceTM->pModule, "getEnergyAndForces");
		/* if(interfaceTM->pGetEnergyAndForces) fprintf(stderr,"ok getEnergyAndForces\n");*/
		if(!interfaceTM->pGetEnergyAndForces) fprintf(stderr," problem with getEnergyAndForces\n");
	}
	if (interfaceTM->pGetEnergyAndForces)
	{
        	interfaceTM->pSetData = PyObject_GetAttrString(interfaceTM->pModule, "setData");
		/* if(interfaceTM->pSetData) fprintf(stderr,"ok setData\n");*/
		if(!interfaceTM->pSetData) fprintf(stderr," problem with setData\n");
	}
	if (interfaceTM->pSetData)
	{
        	interfaceTM->pSetCoordinates = PyObject_GetAttrString(interfaceTM->pModule, "setCoordinates");
		/* if(interfaceTM->pSetCoordinates) fprintf(stderr,"ok pSetCoordinates\n");*/
		if(!interfaceTM->pSetCoordinates) fprintf(stderr," problem with pSetCoordinates\n");
	}
	if (interfaceTM->pSetCoordinates)
	{
        	interfaceTM->pCloseSession = PyObject_GetAttrString(interfaceTM->pModule, "closeSession");
		/* if(interfaceTM->pCloseSession) fprintf(stderr,"ok pCloseSession\n");*/
		if(!interfaceTM->pCloseSession) fprintf(stderr," problem with pCloseSession\n");
	}
	if (interfaceTM->pCloseSession)
	{
        	interfaceTM->pOptGeom = PyObject_GetAttrString(interfaceTM->pModule, "optGeom");
		/* if(interfaceTM->pOptGeom) fprintf(stderr,"ok pOptGeom\n");*/
		if(!interfaceTM->pOptGeom) fprintf(stderr," problem with pOptGeom\n");
	}
	interfaceTM->initialized = 1;
}
InterfaceTM *newInterfaceTM(char* tmModule)
{
	InterfaceTM *interfaceTM= malloc(sizeof(InterfaceTM));
	interfaceTM->pName = NULL;
	interfaceTM->pModule = NULL;
	interfaceTM->pManager = NULL;
	interfaceTM->pGetManager = NULL;
	interfaceTM->pGetEnergy = NULL;
	interfaceTM->pGetEnergyAndForces = NULL;
	interfaceTM->pSetData = NULL;
	interfaceTM->pSetCoordinates = NULL;
	interfaceTM->pOptGeom = NULL;
	interfaceTM->pa = NULL;
	interfaceTM->initialized = 0;
	interfaceTMInitialize(interfaceTM,tmModule);
	return interfaceTM;
}
int interfaceTMComputeGradients(InterfaceTM* interfaceTM, Molecule* mol)
{
	if(!interfaceTM->pa) initManager(interfaceTM, mol);
	if(interfaceTM->pa && interfaceTM->pManager) 
	{
               	/* fprintf(stderr,"ok pa & pManage in interfaceTMComputeGradients\n");*/
		if (interfaceTM->pGetEnergyAndForces && PyCallable_Check(interfaceTM->pGetEnergyAndForces))
		{
			PyObject* listRes = NULL;
                	/* fprintf(stderr,"ok pGetEnergyAndForces in interfaceTMComputeGradients\n");*/
			/*
			PyObject *pArgs = PyTuple_New(2);
                	PyTuple_SetItem(pArgs, 0, interfaceTM->pa );
                	PyTuple_SetItem(pArgs, 1, interfaceTM->pManager );
			PyObject* pEner = PyObject_CallObject(interfaceTM->pGetEnergy, pArgs);
			Py_DECREF(pArgs);
			*/
			setCoordinates(interfaceTM, mol);
               		/* fprintf(stderr,"End setCoordinates in interfaceTMComputeGradients\n");*/
			listRes = PyObject_CallFunctionObjArgs(interfaceTM->pGetEnergyAndForces, interfaceTM->pa,  interfaceTM->pManager,  NULL);
               		/* fprintf(stderr,"End PyObject_CallFunctionObjArgs in interfaceTMComputeGradients\n");*/
            		if (listRes != NULL)
			{
				int i,j,k;
				double energy = PyFloat_AsDouble(PyList_GetItem(listRes, 0));
				double dipole[3];
				for(k=0;k<3;k++) dipole[k] = PyFloat_AsDouble(PyList_GetItem(listRes, k+1));
                		/* fprintf(stderr,"Result of call: Ener = %0.14lf dipole= %0.14lf %0.14lf %0.14lf\n", energy,dipole[0], dipole[1],dipole[2]);*/
                		//Py_DECREF(pEner);
				mol->potentialEnergy = energy*AUTOKCAL;
				for(k=0;k<3;k++) mol->dipole[k] = dipole[k]*AUTODEB;
				j=4;
				for(i=0;i<mol->nAtoms;i++)
					for(k=0;k<3;k++) 
					{
						mol->atoms[i].gradient[k] = AUTOKCAL/BOHRTOANG*PyFloat_AsDouble(PyList_GetItem(listRes, j));
						j++;
					}
				for(i=0;i<mol->nAtoms;i++)
				{
					mol->atoms[i].charge =PyFloat_AsDouble(PyList_GetItem(listRes, j));
					j++;
				}

        		}
			else
                		fprintf(stderr,"listRed of pGetEnergyAndForces = NULL in interfaceTMComputeGradients\n");
		}
	}
	/* fprintf(stderr," interfaceTMComputeGradients\n");*/
	return 0;
}
int interfaceTMComputeEnergy(InterfaceTM* interfaceTM, Molecule* mol)
{
	if(!interfaceTM->pa) initManager(interfaceTM, mol);
	if(interfaceTM->pa && interfaceTM->pManager) 
	{
               	/* fprintf(stderr,"ok pa & pManage in interfaceTMComputeEnergy\n");*/
		if (interfaceTM->pGetEnergy && PyCallable_Check(interfaceTM->pGetEnergy))
		{
			PyObject* listRes = NULL;
                	/* fprintf(stderr,"ok pGetEnergy in interfaceTMComputeEnergy\n");*/
			/*
			PyObject *pArgs = PyTuple_New(2);
                	PyTuple_SetItem(pArgs, 0, interfaceTM->pa );
                	PyTuple_SetItem(pArgs, 1, interfaceTM->pManager );
			PyObject* pEner = PyObject_CallObject(interfaceTM->pGetEnergy, pArgs);
			Py_DECREF(pArgs);
			*/
			setCoordinates(interfaceTM, mol);
               		/* fprintf(stderr,"End setCoordinates in interfaceTMComputeEnergy\n");*/
			listRes = PyObject_CallFunctionObjArgs(interfaceTM->pGetEnergy, interfaceTM->pa,  interfaceTM->pManager,  NULL);
            		if (listRes != NULL)
			{
				int i,j,k;
				double energy = PyFloat_AsDouble(PyList_GetItem(listRes, 0));
				double dipole[3];
				for(k=0;k<3;k++) dipole[k] = PyFloat_AsDouble(PyList_GetItem(listRes, k+1));
                		/* fprintf(stderr,"Result of call: Ener = %0.14lf dipole= %0.14lf %0.14lf %0.14lf\n", energy,dipole[0], dipole[1],dipole[2]);*/
				mol->potentialEnergy = energy*AUTOKCAL;
				for(k=0;k<3;k++) mol->dipole[k] = dipole[k]*AUTODEB;
				for(i=0;i<mol->nAtoms;i++) for(k=0;k<3;k++) mol->atoms[i].gradient[k] = 0;
				j=4;
				for(i=0;i<mol->nAtoms;i++)
				{
					mol->atoms[i].charge =PyFloat_AsDouble(PyList_GetItem(listRes, j));
					j++;
				}

        		}
			else
			{
                		fprintf(stderr,"pEnergy = NULL in interfaceTMComputeEnergy\n");
				return 1;
			}
		}
	}
	/* fprintf(stderr," interfaceTMComputeEnergy\n");*/
	return 0;
}
int interfaceTMOpt(InterfaceTM* interfaceTM, Molecule* mol)
{
	if(!interfaceTM->pa) initManager(interfaceTM, mol);
	if(interfaceTM->pa && interfaceTM->pManager) 
	{
               	/* fprintf(stderr,"ok pa & pManage in interfaceTMComputeGradients\n");*/
		if (interfaceTM->pOptGeom && PyCallable_Check(interfaceTM->pOptGeom))
		{
			PyObject* listRes = NULL;
                	/* fprintf(stderr,"ok pOptGeom in interfaceTMOpt\n");*/
			setCoordinates(interfaceTM, mol);
			listRes = PyObject_CallFunctionObjArgs(interfaceTM->pOptGeom, interfaceTM->pa,  interfaceTM->pManager,  NULL);
            		if (listRes != NULL)
			{
				int i,j,k;
				double energy = PyFloat_AsDouble(PyList_GetItem(listRes, 0));
				double dipole[3];
				for(k=0;k<3;k++) dipole[k] = PyFloat_AsDouble(PyList_GetItem(listRes, k+1));
                		/* fprintf(stderr,"Result of call: Ener = %0.14lf dipole= %0.14lf %0.14lf %0.14lf\n", energy,dipole[0], dipole[1],dipole[2]);*/
                		//Py_DECREF(pEner);
				mol->potentialEnergy = energy*AUTOKCAL;
				for(k=0;k<3;k++) mol->dipole[k] = dipole[k]*AUTODEB;
				j=4;
				for(i=0;i<mol->nAtoms;i++)
					for(k=0;k<3;k++) 
					{
						mol->atoms[i].coordinates[k] = PyFloat_AsDouble(PyList_GetItem(listRes, j));
						j++;
					}
				for(i=0;i<mol->nAtoms;i++)
					for(k=0;k<3;k++) 
					{
						mol->atoms[i].gradient[k] = AUTOKCAL/BOHRTOANG*PyFloat_AsDouble(PyList_GetItem(listRes, j));
						j++;
					}
				for(i=0;i<mol->nAtoms;i++)
				{
					mol->atoms[i].charge =PyFloat_AsDouble(PyList_GetItem(listRes, j));
					j++;
				}

        		}
			else
                		fprintf(stderr,"listRed of pOptGeom = NULL in iinterfaceTMOpt\n");
		}
	}
	/* fprintf(stderr," interfaceTMOpt\n");*/
	return 0;
}
static int interfaceCloseSession(InterfaceTM* interfaceTM)
{
	if(interfaceTM->pa && interfaceTM->pManager) 
	{
		if (interfaceTM->pCloseSession && PyCallable_Check(interfaceTM->pCloseSession))
			PyObject_CallFunctionObjArgs(interfaceTM->pCloseSession, NULL);
	}
	/* fprintf(stderr," interfaceCloseSession\n");*/
	return 0;
}
void interfaceTMDestroy(InterfaceTM* interfaceTM)
{
	interfaceCloseSession(interfaceTM);
	/* Py_FinalizeEx();*/
}
int runTensorMol(Molecule* mol, char* moduleName, int computeGradients)
{
	InterfaceTM* interfaceTM = newInterfaceTM(moduleName);
	int err = 0;
	if(computeGradients) interfaceTMComputeGradients(interfaceTM, mol);
	else err = interfaceTMComputeEnergy(interfaceTM, mol);
	interfaceTMDestroy(interfaceTM);
	/* fprintf(stderr," runTensorMol\n");*/
	return err;
}
int runTM(Molecule* mol, char* moduleName, int computeGradients)
{
	return runTensorMol(mol,moduleName,computeGradients);
}
int runOptTensorMol(Molecule* mol, char* moduleName)
{
	InterfaceTM* interfaceTM = newInterfaceTM(moduleName);
	int err = interfaceTMOpt(interfaceTM, mol);
	interfaceTMDestroy(interfaceTM);
	/* fprintf(stderr," runOptTensorMol\n");*/
	return err;
}
#else
int runTM(Molecule* mol, char* moduleName, int computeGradients)
{
	exit(1);/* compilatation with python is required for TensorMol interface */
	return 1;
}
int runOptTensorMol(Molecule* mol, char* moduleName)
{
	exit(1);/* compilatation with python is required for TensorMol interface */
	return 1;
}
#endif  /* ENABLE_PYTHON*/
