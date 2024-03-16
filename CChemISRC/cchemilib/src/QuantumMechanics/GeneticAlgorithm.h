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

#ifndef __CCHEMILIB_GNETICALGORITHM_H__
#define __CCHEMILIB_GNETICALGORITHM_H__

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>

#include "../Utils/Utils.h"
#include "../Utils/AtomsProp.h"
#include "../Utils/Constants.h"
#include "../Utils/Types.h"
#include "../Molecule/Molecule.h"
#include "../QuantumMechanics/QuantumMechanicsModel.h"

typedef enum { GA_CROSS_PLANE=0, GA_CROSS_SPHERICAL } CCHEMIGACrossType;
typedef enum { GA_MUTATION_SPHERICAL_LOCAL=0, GA_MUTATION_SPHERICAL_GLOBAL } CCHEMIGAMutationType;

typedef struct _GeneticAlgorithm  GeneticAlgorithm;
typedef struct _GeneticAlgorithmClass  GeneticAlgorithmClass;

struct _GeneticAlgorithm
{
/*
 *		maxGens: Maximum number of generations to evaluate.
 *		popSize: Size of the population
 *		pCross: Probability of cross between individuals.
 *		pMutat: Probability of mutation in an idividual.
 *		low: Vector that containts the lower bounds of the solution space.
 *		up: Vector that containts the upper bounds of the solution space.
 */
	GeneticAlgorithmClass* klass;
	QuantumMechanicsModel** pop;
	QuantumMechanicsModel** childs;
	int* parent1;
	int* parent2;
	int maxGens;
	int popSize;
	int nChilds;
	double pCross;
	double pMutation;
	FILE* logfile;
	double* fun;
	double* FitnessFun;
	double emean;
	double sigmax;
	double fmin;
	int numGen;
	int nGeoms;

	char* inputFileName;
	char* suffixFileName;
	double* pseudoInertia;
	double* childsPseudoInertia;
	CCHEMIGACrossType crossType;
	CCHEMIGAMutationType mutationType;
	FILE* fileAllEnergies;
	FILE* fileSelectedEnergies;
	boolean removeSimilarInertia;
	double inertiaTol;
	boolean removeFragmented;
	boolean removeSmallDistance;
	boolean removeSimilarBonds;
	double sTol;
	double distMaxTol ;
	boolean chain ;
	boolean saveFirstGeom ;
	int nTimesGeoms;
	boolean constraints ;
	char* command;
	char* QMKeys;
	char* dirName;
};
struct _GeneticAlgorithmClass
{
	void(*computeChildsEnergies)(GeneticAlgorithm* ga);
	void(*computeEnergies)(GeneticAlgorithm* ga);
	void(*saveGeometries)(GeneticAlgorithm* ga);
	void(*evaluatePopulation)(GeneticAlgorithm* ga);
	void(*elitist)(GeneticAlgorithm* ga);
	void(*applyMutation)(GeneticAlgorithm* ga);
	void(*makeSelection)(GeneticAlgorithm* ga);
	void(*makeCrossover)(GeneticAlgorithm* ga);
	void(*free)(GeneticAlgorithm* ga);
	void(*run)(GeneticAlgorithm* ga);
};
GeneticAlgorithm* newGeneticAlgorithm(char* inputFileName);
#endif /* __CCHEMILIB_GNETICALGORITHM_H__ */
