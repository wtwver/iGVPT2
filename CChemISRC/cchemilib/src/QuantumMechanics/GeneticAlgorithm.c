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

/* GeneticAlgorithm.c */

#include "QuantumMechanics.h"
#include "GeneticAlgorithm.h"

static void freeGeneticAlgorithm(GeneticAlgorithm* ga);
static void evaluatePopulation(GeneticAlgorithm* ga);
static void saveGeometries(GeneticAlgorithm* ga);
static void elitist(GeneticAlgorithm* ga);
static void applyMutation(GeneticAlgorithm* ga);
static void makeSelection(GeneticAlgorithm* ga);
static void makeCrossover(GeneticAlgorithm* ga);
static void run(GeneticAlgorithm* ga);
static boolean runOneOptGeneric(QuantumMechanicsModel* qm, char* fileNamePrefix, FILE* logfile);
static boolean runGenericFiles(int numberOfGeometries, QuantumMechanicsModel** qm, char* fileNamePrefix, FILE* logfile);
static void computeEnergies(GeneticAlgorithm* ga);
static void computeChildsEnergies(GeneticAlgorithm* ga);
static void computePseudoInertia(GeneticAlgorithm* ga);
static void computePopPseudoInertia(GeneticAlgorithm* ga);
static void computeChildsPseudoInertia(GeneticAlgorithm* ga);
static void printPopEnergies(GeneticAlgorithm* ga);
static void printChildsEnergies(GeneticAlgorithm* ga);
static void printEnergies(GeneticAlgorithm* ga);
static void savePopGeometries(GeneticAlgorithm* ga);
static void saveChildsGeometries(GeneticAlgorithm* ga);
static void savePopEnergies(GeneticAlgorithm* ga, FILE* file);
static void saveChildsEnergies(GeneticAlgorithm* ga, FILE* file);
/*********************************************************************************************************************/
static void freeGeneticAlgorithm(GeneticAlgorithm* ga)
{
	if(ga->fun) free(ga->fun);
	if(ga->FitnessFun) free(ga->FitnessFun);
	if(ga->pseudoInertia) free(ga->pseudoInertia);
	if(ga->childsPseudoInertia) free(ga->childsPseudoInertia);
	if(ga->suffixFileName) free(ga->suffixFileName);
	if(ga->inputFileName) free(ga->inputFileName);
}
/*********************************************************************************************************************/
static void computeEnergies(GeneticAlgorithm* ga)
{
	runGenericFiles(ga->popSize, ga->pop, ga->suffixFileName,ga->logfile);
	// TEST, TO REMOVE
	/*
	computePopPseudoInertia(ga); 
	for (int i=0; i<ga->popSize; i++) 
		ga->pop[i]->molecule.potentialEnergy=ga->pseudoInertia[i];
		*/
}
/*********************************************************************************************************************/
static void computeChildsEnergies(GeneticAlgorithm* ga)
{
	runGenericFiles(ga->nChilds, ga->childs, ga->suffixFileName, ga->logfile);
	// TEST, TO REMOVE
	/*
	computeChildsPseudoInertia(ga); 
	for (int i=0; i<ga->nChilds; i++) 
		if(ga->childs[i])
			ga->childs[i]->molecule.potentialEnergy=ga->childsPseudoInertia[i];
			*/
}
/*********************************************************************************************************************/
static void computePopPseudoInertia(GeneticAlgorithm* ga)
{
	int i;
	for (i=0; i<ga->popSize; i++) 
	{
		Molecule* mol =  &ga->pop[i]->molecule;
		double I2;
		double I4;
		mol->klass->computePseudoInertia(mol, &I2, &I4);
		ga->pseudoInertia[i] = I4+I2;
	}
}
/*********************************************************************************************************************/
static void computeChildsPseudoInertia(GeneticAlgorithm* ga)
{
	int i;
	for (i=0; i<ga->nChilds; i++) 
	{
		if(!ga->childs[i]) continue;
		Molecule* mol =  &ga->childs[i]->molecule;
		double I2;
		double I4;
		mol->klass->computePseudoInertia(mol, &I2, &I4);
		ga->childsPseudoInertia[i] = I4+I2;
	}
}
/*********************************************************************************************************************/
static void computePseudoInertia(GeneticAlgorithm* ga)
{
	computePopPseudoInertia(ga);
	computeChildsPseudoInertia(ga);
}
/*********************************************************************************************************************/
static void printPopEnergies(GeneticAlgorithm* ga)
{
	fprintf(ga->logfile,"Energies of generation # %-12d : ",ga->numGen);
	int i;
	for (i=0; i<ga->popSize; i++) 
	{
		Molecule* mol =  &ga->pop[i]->molecule;
		fprintf(ga->logfile,"%-14.10lf ",mol->potentialEnergy);
		if((i+1)%5==0) fprintf(ga->logfile,"\n %39s"," ");
	}
	fprintf(ga->logfile,"\n");
	fflush(ga->logfile);
}
/*********************************************************************************************************************/
static void printChildsEnergies(GeneticAlgorithm* ga)
{
	fprintf(ga->logfile,"Childs energies of generation # %-5d : ",ga->numGen);
	int i;
	for (i=0; i<ga->nChilds; i++) 
	{
		if(!ga->childs[i]) continue;
		Molecule* mol =  &ga->childs[i]->molecule;
		fprintf(ga->logfile,"%-14.10lf ",mol->potentialEnergy);
		if((i+1)%5==0) fprintf(ga->logfile,"\n %39s"," ");
	}
	fprintf(ga->logfile,"\n");
	fflush(ga->logfile);
}
/*********************************************************************************************************************/
static void printEnergies(GeneticAlgorithm* ga)
{
	printPopEnergies(ga);
	printChildsEnergies(ga);
}
/*********************************************************************************************************************/
static FILE* openFile(GeneticAlgorithm* ga, char* name, char* message)
{
	char* suff = ga->suffixFileName;
	FILE* file=NULL;
	char tmp[BSIZE];
	sprintf(tmp,"%s_%s_GA.txt",suff,name);
	file = fopen(tmp,"w");
	if(message) fprintf(ga->logfile,"%s %s\n",message,tmp);
	return file;
}
/*********************************************************************************************************************/
static void savePopEnergies(GeneticAlgorithm* ga,FILE* file)
{
	int i;
	for (i=0; i<ga->popSize; i++) 
	{
		Molecule* mol =  &ga->pop[i]->molecule;
		fprintf(file,"%-20.12lf\n",mol->potentialEnergy);
	}
	fflush(file);
}
/*********************************************************************************************************************/
static void saveChildsEnergies(GeneticAlgorithm* ga,FILE* file)
{
	int i;
	for (i=0; i<ga->nChilds; i++) 
	{
		if(!ga->childs[i]) continue;
		Molecule* mol =  &ga->childs[i]->molecule;
		fprintf(file,"%-20.12lf\n",mol->potentialEnergy);
	}
	fflush(file);
}
/*********************************************************************************************************************/
static void evaluatePopulation(GeneticAlgorithm* ga)
{
	int i;
	for (i=0; i<ga->popSize; i++) ga->fun[i] = ga->pop[i]->molecule.potentialEnergy;
	// Evaluate the population
	double emax =  ga->pop[0]->molecule.potentialEnergy;
	double emin =  ga->pop[0]->molecule.potentialEnergy;
	for (i=0; i<ga->popSize; i++) 
		if(emin>ga->pop[i]->molecule.potentialEnergy) 
			emin=ga->pop[i]->molecule.potentialEnergy;
	for (i=0; i<ga->popSize; i++) 
		if(emax<ga->pop[i]->molecule.potentialEnergy) 
			emax=ga->pop[i]->molecule.potentialEnergy;

	double de=emax-emin;
	if(de<=0) de =1;
	// Population fitness
	for (i=0; i<ga->popSize; i++) 
	{
		// Roy L. Johnston Dalton Trans. 2003, 4193-4207 
		double rho = (ga->pop[i]->molecule.potentialEnergy-emin)/de;
		ga->FitnessFun[i] = exp(-3.0*rho);
		/*
		ga->FitnessFun[i] = 1-0.7*rho;
		ga->FitnessFun[i] = 0.5*(1-tanh(2*rho-1));
		*/
	}
	//cout<<"End Calcul "<<endl;

	ga->fmin = ga->fun[0];
	for (i=1; i<ga->popSize; i++) 
		if(ga->fmin > ga->fun[i]) 
			ga->fmin = ga->fun[i]; 


}
/*********************************************************************************************************************/
/*
static void saveGeometries(GeneticAlgorithm* ga)
{
	char* suff = ga->suffixFileName;
	char tmp[BSIZE];
	for (int i=0; i<ga->popSize; i++) 
	{
		sprintf(tmp,"%s%d_GA.gab",suff,i+1);
		Molecule* mol= &ga->pop[i]->molecule;
		mol->klass->save(mol,tmp);
	}
}
*/

/*********************************************************************************************************************/
static void saveGeometriesInGabeditFile(GeneticAlgorithm* ga,char* fileNameGab)
{
	FILE* file = NULL;
	int form = 1;
	int nG=0;
	int* index =NULL;
	int i;

	if(ga->popSize<1) return;
	if(!ga->pop) return;
	for(i=0;i<ga->popSize;i++) if(ga->pop[i]) nG++;
	if(nG<1) return;
	//fprintf(stderr,"nG=%d\n",nG); fflush(stderr);
 	file = fopen(fileNameGab, "w");
	if(!file) return;
	// index by energy values
	index = malloc(ga->popSize*sizeof(int));
	for(i=0;i<ga->popSize;i++) index[i] = i;

	for(i=0;i<ga->popSize;i++) 
	{
		int ii=index[i];
		if(!ga->pop[ii]) continue;
		int k=i;
		int j;
		for(j=i+1;j<ga->popSize;j++) 
		{
			int jj=index[j];
			if(!ga->pop[jj]) continue;
			if(ga->pop[jj]->molecule.potentialEnergy<ga->pop[index[k]]->molecule.potentialEnergy) k=j;

		}
		if(k!=i)
		{
			int t=index[i];
			index[i]= index[k];
			index[k]= t;
		}
	}
	//fprintf(stderr,"end sorting\n"); fflush(stderr);


	fprintf(file,"[Gabedit Format]\n");
	fprintf(file,"[GEOCONV]\n");
	fprintf(file,"energy\n");
	for(i=0;i<ga->popSize;i++) if(ga->pop[index[i]]) fprintf(file,"%f\n",ga->pop[index[i]]->molecule.potentialEnergy);
	fprintf(file,"max-force\n");
	for(i=0;i<nG;i++) fprintf(file,"0.0\n");
	fprintf(file,"rms-force\n");
	for(i=0;i<nG;i++) fprintf(file,"0.0\n");
	fprintf(file,"\n");
	fprintf(file,"[GEOMETRIES]\n");
	for (i=0; i<ga->popSize; i++) 
	{
		int id=index[i];
		Molecule* mol;
		if(!ga->pop[id]) continue;
		mol =  &ga->pop[id]->molecule;

		fprintf(file,"%d\n",mol->nAtoms);
		fprintf(file,"%d %d\n",mol->totalCharge, mol->spinMultiplicity);
		int j;
		for(j=0;j<mol->nAtoms;j++)
		fprintf(file," %s %0.8f %0.8f %0.8f\n", 
				mol->atoms[j].prop.symbol,
				mol->atoms[j].coordinates[0],
				mol->atoms[j].coordinates[1],
				mol->atoms[j].coordinates[2]
				);
	}
	//fprintf(stderr,"end save first part\n"); fflush(stderr);
	fprintf(file,"\n");
	fprintf(file,"[GEOMS] %d\n",form);
	fprintf(file,"%d 2\n",nG);
	fprintf(file,"energy kcal/mol 1\n");
	fprintf(file,"deltaE eV 1\n");
	int k=-1;
	double e0=0;
	for (i=0; i<ga->popSize; i++) 
	{
		int id=index[i];
		Molecule* mol;
		if(!ga->pop[id]) continue;
		mol =  &ga->pop[id]->molecule;
		if(k<0){ 
			k=i; 
			e0=ga->pop[index[k]]->molecule.potentialEnergy;
		}
		fprintf(file,"%f\n",mol->potentialEnergy);

		//if(k>=0) fprintf(file,"%f\n",(mol->potentialEnergy-e0)*503.21892494);// in K
		if(k>=0) fprintf(file,"%f\n",(mol->potentialEnergy-e0)*0.04336410);// in eV
		else fprintf(file,"0\n");
		mol->klass->addGeometryToGabedit(mol, file);
	}
	//fprintf(stderr,"end save second part\n"); fflush(stderr);
	fclose(file);
	k = -1;
	free(index);
	fprintf(ga->logfile,"Population of the GA are saved in %s\n",fileNameGab);
}
/*********************************************************************************************************************/
static void saveGeometries(GeneticAlgorithm* ga)
{
	char* suff = ga->suffixFileName;
	char tmp[BSIZE];
	sprintf(tmp,"%s_Gen%d_GA.gab",suff,ga->numGen);
	saveGeometriesInGabeditFile(ga,tmp); 
}
/*********************************************************************************************************************/
static void saveFinalGeometries(GeneticAlgorithm* ga)
{
	char tmp[BSIZE];
	int nG=0;
	int i;

	if(ga->popSize<1) return;
	if(!ga->pop) return;
	for(i=0;i<ga->popSize;i++) if(ga->pop[i]) nG++;
	if(nG<1) return;
	
	QuantumMechanicsModel* qm=NULL;
	int nGeoms=ga->popSize;
	for(i=0;i<ga->popSize;i++) { if(ga->pop[i]) qm=ga->pop[i]; break;}
	nGeoms =ga->popSize;
	if(qm)
	{
		if(ga->removeSimilarInertia) 
			qm->klass->removeSimilarInertiaGeometries(ga->pop, &nGeoms, ga->fun,ga->logfile,ga->inertiaTol); 
		if(ga->removeFragmented) 
			qm->klass->removeFragmentedMolecules(ga->pop, &nGeoms, ga->fun, ga->logfile);
		if(ga->removeSmallDistance) 
			qm->klass->removeSmallDistanceMolecules(ga->pop, &nGeoms, ga->fun, ga->logfile);
		if(ga->removeSimilarBonds) 
			qm->klass->removeSimilarBondsGeometries(ga->pop, &nGeoms, ga->fun,ga->logfile,ga->sTol, ga->distMaxTol);
	}
	ga->popSize=nGeoms;

	//fprintf(stderr,"nG=%d\n",nG); fflush(stderr);

	sprintf(tmp,"%s_FinalGen_GA.gab",ga->suffixFileName);
	saveGeometriesInGabeditFile(ga,tmp); 
}
/*********************************************************************************************************************/
static void saveGeometry(Molecule* mol, int nGeoms, char* suff)
{
	char fname[BSIZE];
	sprintf(fname,"%s_SelectedGeom%d_GA.gab",suff, nGeoms);
	mol->klass->save(mol,fname);
}
/*********************************************************************************************************************/
static void savePopGeometries(GeneticAlgorithm* ga)
{
	char* suff = ga->suffixFileName;
	char fname[BSIZE];
	int ibegin= 1;
	int i;
	for (i=0; i<ga->popSize; i++) 
	{
		Molecule* mol = &ga->pop[i]->molecule;
		sprintf(fname,"%s_GeomNumber%d_GA.gab",suff,ibegin+i);
		mol->klass->save(mol,fname);
	}
}
/*********************************************************************************************************************/
static void saveChildsGeometries(GeneticAlgorithm* ga)
{
	char* suff = ga->suffixFileName;
	char fname[BSIZE];
	int ibegin= ga->popSize+(ga->numGen-1)*ga->nChilds+1;
	int i;
	for (i=0; i<ga->nChilds; i++) 
	{
		if(!ga->childs[i]) continue;
		Molecule* mol = &ga->childs[i]->molecule;
		sprintf(fname,"%s_GeomNumber%d_GA.gab",suff,ibegin+i);
		mol->klass->save(mol,fname);
	}
}
/*********************************************************************************************************************/
static void elitist(GeneticAlgorithm* ga)
{
	// Elitist
	int ic;
	for(ic=0;ic<ga->nChilds;ic++)
	{
		if(!ga->childs[ic]) continue;
		Molecule* mol = &ga->childs[ic]->molecule;
		double childEnergy = ga->childs[ic]->molecule.potentialEnergy;
		int ir=-1;
		int i;
		for (i=0; i<ga->popSize; i++) 
			if(mol->klass->similarInertia(mol, &ga->pop[i]->molecule, 0.04)) { 
				ir = i; 
				break;
			}

		if(ir==-1)
		{
			int imax = 0;
			int i;
			for (i=1; i<ga->popSize; i++) 
				if(ga->pop[i]->molecule.potentialEnergy > ga->pop[imax]->molecule.potentialEnergy) imax = i;
			if(ga->pop[imax]->molecule.potentialEnergy > childEnergy) ir = imax;
		}
		//fprintf(ga->logfile,"ir=%d\n",ir); fflush(ga->logfile);

		if(ir>-1 && childEnergy<ga->pop[ir]->molecule.potentialEnergy)
		{
			ga->nGeoms++;
			saveGeometry(mol, ga->nGeoms,ga->suffixFileName);
			fflush(ga->logfile);
			fprintf(ga->fileSelectedEnergies,"%-20.12lf\n",mol->potentialEnergy);

			ga->pop[ir]->klass->free(ga->pop[ir]);
			ga->pop[ir] = ga->childs[ic];
			ga->childs[ic]=NULL;
		}
		else
		{
			ga->childs[ic]->klass->free(ga->childs[ic]);
			ga->childs[ic]=NULL;
		}
	}
	ga->klass->evaluatePopulation(ga);

}
/*********************************************************************************************************************/
static void applyMutation(GeneticAlgorithm* ga)
{
	//int itr = ga->numGen-1;
	int i;
	for (i=0; i<ga->nChilds; i++)
	{
		double r = rand()/(double)(RAND_MAX);
		if(r < ga->pMutation) 
		{
			if(ga->mutationType ==GA_MUTATION_SPHERICAL_LOCAL)
				ga->childs[i]->molecule.klass->makeLocalSphericalMutation(&ga->childs[i]->molecule,0.05);
			else if(ga->mutationType ==GA_MUTATION_SPHERICAL_GLOBAL)
				ga->childs[i]->molecule.klass->makeCenterOfMassSphericalMutation(&ga->childs[i]->molecule);
		}
	}
}
/*********************************************************************************************************************/
static void makeSelection(GeneticAlgorithm* ga)
{
	double p[ga->popSize];
	double q[ga->popSize];
	double fit = 0;
	double cProb = 0;
	double r;
	//int crossPos[ga->popSize];

	fit = 0;
	int i;
	for (i=0; i<ga->popSize; i++) fit += ga->FitnessFun[i];
	q[0] = 0;
	cProb = 0;
	// Cummulative probability of individuals
	for (i=0; i<ga->popSize; i++)
	{
		p[i] = ga->FitnessFun[i]/fit;
		cProb += p[i];
		q[i] = cProb;
	}
	// Selection process
	for (i=0; i<ga->nChilds; i++)
	{
		int j;
		r = rand()/(double)(RAND_MAX);
		if (r < q[0]) ga->parent1[i]=i;
		else for (j=1; j<ga->popSize; j++) if ( (r > q[j-1]) && (r < q[j]) ) ga->parent1[i] = j;
		for (j=0; j<ga->popSize; j++) 
		{
			if(ga->parent1[i] == j) continue;
			ga->parent2[i]=j;
			boolean alreadyUsed=FALSE;
			int k;
			for (k=0; k<i; k++) 
			{
				if(ga->parent1[i] == ga->parent1[k] && ga->parent2[i] == ga->parent2[k] ) { alreadyUsed=TRUE; break;}
				if(ga->parent1[i] == ga->parent2[k] && ga->parent2[i] == ga->parent1[k] ) { alreadyUsed=TRUE; break;}
			}
			if(!alreadyUsed) break;
		}

	}
}
/*********************************************************************************************************************/
static void makeCrossover(GeneticAlgorithm* ga)
{
	// Arithmetic Crossover
	Molecule* molChild = NULL;
	QuantumMechanicsModel* p1=NULL;
	QuantumMechanicsModel* p2=NULL;
	int err = 0;
	//for (int i=0; i<ga->nChilds; i++) { fprintf(stderr,"Child %d P1=%d P2=%d\n",i+1,ga->parent1[i], ga->parent2[i]); fflush(stderr); }

	int i;
	for (i=0; i<ga->nChilds; i++)
	{
		p1=ga->pop[ga->parent1[i]];
		p2=ga->pop[ga->parent2[i]];
		if(ga->childs[i]) ga->childs[i]->klass->free(ga->childs[i]);
		ga->childs[i] = malloc(sizeof(QuantumMechanicsModel));
		Molecule* mol1 = &p1->molecule;
		Molecule* mol2 = &p2->molecule;
		err = 0;
		if(ga->crossType == GA_CROSS_PLANE) molChild = mol1->klass->makePlaneCutSpliceCrossover(mol1, mol2, &err);
		else if(ga->crossType == GA_CROSS_SPHERICAL ) molChild = mol1->klass->makeSphereCutSpliceCrossover(mol1, mol2, &err);

		//if(err!=0) continue;
		*ga->childs[i] = createGenericModel(molChild, p1->method, p1->workDir, p1->nameCommand, p1->molecule.constraints, p1->logfile);
	}
}
/*********************************************************************************************************************/
static void printLine(FILE* logfile, char* s)
{
	int i;
        for(i=0;i<120;i++) 
		fprintf(logfile,"%s",s); 
	fprintf(logfile,"\n");
}
/*********************************************************************************************************************/
static void run(GeneticAlgorithm* ga)
{
	ga->fileAllEnergies = openFile(ga, "AllEnergies", "All energies saved in");
	ga->fileSelectedEnergies = openFile(ga, "SelectedEnergies", "Selected energies saven in");
	// TO DO 
	ga->klass->computeEnergies(ga);
	computePseudoInertia(ga);

	//fprintf(stderr,"evaluatePopulation\n"); fflush(stderr);
	ga->klass->evaluatePopulation(ga);
	//fprintf(stderr,"end evaluatePopulation\n"); fflush(stderr);
	savePopGeometries(ga);
	savePopEnergies(ga,ga->fileAllEnergies);
	savePopEnergies(ga,ga->fileSelectedEnergies);
	
	// Evaluate Generations
	//fprintf(stderr,"initial value %f\n",ga->fmin); fflush(stderr);
	ga->klass->saveGeometries(ga);
	//fprintf(stderr,"end save geometries\n"); fflush(stderr);
	printLine(ga->logfile, "=");
        fprintf(ga->logfile,"Step # %d/%d Fmin= %f\n",0,ga->maxGens,ga->fmin);
	printLine(ga->logfile, "=");
	int itr;
	for (itr=0; itr<ga->maxGens; itr++)
	{
		/*
                fprintf(stderr,"Step # %d/%d Fmin= %f\n",itr,ga->maxGens,ga->fmin);
                fprintf(stderr,"Begin makeSel\n");
		fflush(stderr);
		*/
		ga->numGen = itr+1;
		ga->klass->makeSelection(ga);

		ga->klass->makeCrossover(ga);
		ga->klass->applyMutation(ga);
                fprintf(ga->logfile,"After crossing & mutation\n");
                fprintf(ga->logfile,"=========================\n");
		ga->klass->computeChildsEnergies(ga);
		computeChildsPseudoInertia(ga);
		printEnergies(ga);
		/*
                fprintf(stderr,"Begin makeCross\n"); fflush(stderr);
		ga->klass->makeCrossover(ga);
                fprintf(ga->logfile,"After crossing\n");
                fprintf(ga->logfile,"==============\n");
		ga->klass->computeChildsEnergies(ga);
		computeChildsPseudoInertia(ga);
		printEnergies(ga);

                fprintf(stderr,"Begin applyMutation\n"); fflush(stderr);
		ga->klass->applyMutation(ga);
                fprintf(ga->logfile,"After mutation\n");
                fprintf(ga->logfile,"==============\n");
		ga->klass->computeChildsEnergies(ga);
		printEnergies(ga);
		*/

		saveChildsEnergies(ga,ga->fileAllEnergies);

		// before elitist to take all geoms inluding if not selected 
		saveChildsGeometries(ga);

                //fprintf(stderr,"Begin elitist\n"); fflush(stderr);
		ga->klass->elitist(ga);
		ga->klass->saveGeometries(ga);
                fprintf(ga->logfile,"Step # %d/%d Fmin= %f\n",itr+1,ga->maxGens,ga->fmin);
		printLine(ga->logfile, "=");
                fflush(ga->logfile);
	}
	saveFinalGeometries(ga);
	if(ga->fileAllEnergies) fclose(ga->fileAllEnergies);
	if(ga->fileSelectedEnergies) fclose(ga->fileSelectedEnergies);
	ga->fileAllEnergies = NULL;
	ga->fileSelectedEnergies = NULL;
}
/*********************************************************************************************************************/
static boolean getEnergyGeneric(char* fileNameOut, double* energy)
{
	FILE* file = NULL;
	char buffer[1024];
	int i;
 	file = fopen(fileNameOut, "r");
	if(!file) return FALSE;
	if(!fgets(buffer,BSIZE,file)) { fclose(file); return FALSE;}/* first line for energy in Hartree*/

	for(i=0;i<strlen(buffer);i++) if(buffer[i]=='D' || buffer[i]=='d') buffer[i] ='E';
	if(sscanf(buffer,"%lf",energy)==1)
	{
		fclose(file);
		*energy *=AUTOKCAL;
		return TRUE;
	}
	fclose(file);
	return FALSE;
}
/***********************************************************************************************************************/
static boolean runOneOptGeneric(QuantumMechanicsModel* qm, char* fileNamePrefix, FILE* logfile)
{
	char* keyWords = "OPT";
	char* genericCommand = qm->nameCommand;
	FILE* file = NULL;
	FILE* fileSH = NULL;
	char* fileNameIn = NULL;
	char* fileNameOut = NULL;
	char* fileNameLog = NULL;
	char* fileNameSH = NULL;
	char buffer[1024];
	int ir;
	int type = 0;
	Molecule* mol = &qm->molecule;
	mol->potentialEnergy=1e10;
	fflush(logfile);
#ifdef OS_WIN32
	char c='%';
#endif
	if(!qm) { fprintf(logfile,"Error qm=NULL\n"); fflush(logfile);}
	if(!qm) return FALSE;
	if(qm->molecule.nAtoms<1) { fprintf(logfile,"Error nAtoms<1\n"); fflush(logfile);}
	if(qm->molecule.nAtoms<1) return FALSE;
#ifndef OS_WIN32
	if(fileNamePrefix[0]=='/') fileNameSH = strdup_printf("%sGeneOne.sh",fileNamePrefix);
	else fileNameSH = strdup_printf("./%sGeneOne.sh",fileNamePrefix);
#else
	fileNameSH = strdup_printf("%sGeneOne.bat",fileNamePrefix);
#endif
 	fileSH = fopen(fileNameSH, "w");
	if(!fileSH) 
	{
		fprintf(logfile,"I cannot create %s\n",fileNameSH); fflush(logfile);
		return FALSE;
	}

	fileNameIn = strdup_printf("%sOne.inp",fileNamePrefix);
	fileNameOut = strdup_printf("%sOne.out",fileNamePrefix);

	fileNameLog = strdup_printf("%sOne.log",fileNamePrefix);

 	file = fopen(fileNameIn, "w");
	if(!file) 
	{
		fprintf(logfile,"I cannot create %s\n",fileNameIn); fflush(logfile);
 		if(fileNameIn) free(fileNameIn);
 		if(fileNameOut) free(fileNameOut);
 		if(fileNameSH) free(fileNameSH);
 		if(fileNameLog) free(fileNameLog);
		return FALSE;
	}

	if(strstr(keyWords,"Opt")) type = 2;
	if(strstr(keyWords,"OPT")) type = 2;
	if(strstr(keyWords,"ENGRAD")) type = 1;
	fprintf(file,"%d\n",type);
	mol->klass->addMolecule(mol,file);
	fclose(file);
	/*
	{
		char* str = NULL;
		if(strstr(keyWords,"OPT")) str = strdup_printf("Minimization by Generic/%s ... Please wait",genericCommand);
		else str = strdup_printf("Computing of energy by Generic/%s .... Please wait ",genericCommand);
		fprintf(logfile,"%s\n",str);
		fflush(logfile);
		if(str) free(str);
	}
	*/
	//fprintf(stdout,"1=====>%s %s %s\n",genericCommand,fileNameIn,fileNameOut);
#ifndef OS_WIN32
	//fprintf(stdout,"2=====>%s %s %s\n",genericCommand,fileNameIn,fileNameOut);
	//fprintf(fileSH,"%s %s %s",genericCommand,fileNameIn,fileNameOut);
	fprintf(fileSH,"#!/bin/bash\n");
	fprintf(fileSH,"%s %s %s >& %s",genericCommand,fileNameIn,fileNameOut, fileNameLog);
	fclose(fileSH);

	//sprintf(buffer,"cat %s",fileNameSH); fprintf(stdout,"%s\n",buffer);system(buffer); fflush(stdout);
	//sprintf(buffer,"cat %s",fileNameIn); ir=system(buffer); fprintf(logfile,"ir=%d\n",ir); fflush(logfile);

	sprintf(buffer,"chmod u+x %s",fileNameSH);
	ir=system(buffer);
	if(ir!=0) { fprintf(logfile,"ir=%d\n",ir); fflush(logfile); }
	ir=system(fileNameSH);
	if(ir!=0) { fprintf(logfile,"ir=%d\n",ir); fflush(logfile); }

	//sprintf(buffer,"cat %s",fileNameOut); ir=system(buffer); fprintf(logfile,"ir=%d\n",ir); fflush(logfile);
#else
	//fprintf(fileSH,"\"%s\" \"%s\" \"%s\"",genericCommand,fileNameIn,fileNameOut);
	fprintf(fileSH,"\"%s\" \"%s\" \"%s\" >& \"%s\"",genericCommand,fileNameIn,fileNameOut, fileNameLog);
	fclose(fileSH);
	sprintf(buffer,"\"%s\"",fileNameSH);
	system(buffer);
#endif
	boolean ok=TRUE;
	if(getEnergyGeneric(fileNameOut,&mol->potentialEnergy))
	{
		fprintf(logfile,"Energy by Generic = %f\n", mol->potentialEnergy);
		fflush(logfile);
		ok=mol->klass->readGeometry(mol,fileNameOut);
	}
	else ok=FALSE;

 	remove(fileNameIn);
 	remove(fileNameOut);
 	remove(fileNameSH);
 	remove(fileNameLog);
 	if(fileNameIn) free(fileNameIn);
 	if(fileNameOut) free(fileNameOut);
 	if(fileNameSH) free(fileNameSH);
 	if(fileNameLog) free(fileNameLog);
	fflush(logfile);
	return ok;
}
/*************************************************************************************************************************************************************/
static boolean runGenericFiles(int numberOfGeometries, QuantumMechanicsModel** qm, char* fileNamePrefix, FILE* logfile)
{
	int i;
	int nG = 0;
	int nM = 0;
	char buffer[BSIZE];
	for(i=0;i<numberOfGeometries;i++)
	{
		if(!qm[i]) continue;
		nG++;
		fprintf(logfile,"Minimization by Generic of geometry n = %d... Please wait\n", i+1);
		fflush(logfile);
		sprintf(buffer,"%s%d",fileNamePrefix,i);
		if(runOneOptGeneric(qm[i], buffer, logfile)) 
		{
			nM++;
		}
		else
		{
			qm[i]->klass->free(qm[i]);
			qm[i] =NULL;
		}
		fflush(logfile);

	}
	/*
	if(nM==nG) return TRUE;
	return FALSE;
	*/
	fprintf(logfile,"Number of generic runs with errors = %d\n", nG-nM); fflush(logfile);
	fprintf(logfile,"-------------------------------------------\n"); fflush(logfile);
	return (nM>0);

}
/*********************************************************************************************************************/
static boolean checkSizes(GeneticAlgorithm* ga)
{
	if(ga->nChilds<1) ga->nChilds = 1;
	if(ga->popSize<2) return FALSE;
	// can be extended to popSize*(popSize-1)/2;
	if(ga->nChilds>ga->popSize-1) ga->nChilds = ga->popSize-1;
	return TRUE;
}

/*********************************************************************************************************************/
static void initTables(GeneticAlgorithm* ga)
{
	ga->fun =(double*)malloc(ga->popSize*sizeof(double));
	ga->FitnessFun =(double*)malloc(ga->popSize*sizeof(double));
	ga->pseudoInertia =(double*)malloc(ga->popSize*sizeof(double));
	ga->childsPseudoInertia =(double*)malloc(ga->popSize*sizeof(double));
	if(ga->nChilds>0) ga->parent1 =(int*)malloc(ga->nChilds*sizeof(int));
	if(ga->nChilds>0) ga->parent2 =(int*)malloc(ga->nChilds*sizeof(int));
	if(ga->nChilds>0) ga->childs =(QuantumMechanicsModel**)malloc(ga->nChilds*sizeof(QuantumMechanicsModel*));
	int i;
	for(i=0;i<ga->nChilds;i++) ga->childs[i] = NULL;
}
/*********************************************************************************************************************/
static void setMethods(GeneticAlgorithm* ga)
{
	ga->klass = (GeneticAlgorithmClass*)malloc(sizeof(GeneticAlgorithmClass));

	ga->klass->free = freeGeneticAlgorithm;
	ga->klass->run = run;
	ga->klass->evaluatePopulation = evaluatePopulation;
	ga->klass->saveGeometries = saveGeometries ;
	ga->klass->elitist = elitist ;
	ga->klass->applyMutation = applyMutation ;
	ga->klass->makeSelection = makeSelection ;
	ga->klass->makeCrossover = makeCrossover;
	ga->klass->computeEnergies = computeEnergies ;
	ga->klass->computeChildsEnergies = computeChildsEnergies ;
}
/*********************************************************************************************************************/
static void readParameters(GeneticAlgorithm* ga)
{
	char* dirName = NULL;
	char* QMKeys = NULL;
	int maxGen=100;
	int nChilds;
	double pCross=0.5;
	double pMut=0.08;
	FILE* logfile = stdout;
	FILE* file = fopen(ga->inputFileName,"rb");
	char* genericCommand=strdup("runGeneric");
	CCHEMIGACrossType crossType =  GA_CROSS_PLANE;
	CCHEMIGAMutationType mutationType  = GA_MUTATION_SPHERICAL_LOCAL;
	char* tmp = NULL;
	char* model = NULL;
	int popSize=-1;
	boolean chain=FALSE;
	boolean saveFirstGeom=FALSE;
	int nTimesGeoms=1;
	boolean removeSimilarInertia = FALSE;
	boolean removeFragmented = FALSE;
	boolean removeSmallDistance = FALSE;
	Constraints constraints = NOCONSTRAINTS;
	double inertiaTol = 0.04; // recomended byjun Zhao et al (2016)
       	//Comprehensive genetic algorithm for ab initio global optimisation of clusters, Molecular Simulation,
	// 42:10, 809-819, DOI: 10.1080/08927022.2015.1121386	     

	boolean removeSimilarBonds = FALSE;
	double sTol=0.02;
	double distMaxTol=0.7;
	// sTol = 0.02 , distMaxTol = 0.7 Ang, recommanded in Jorgensen et al JCTC, 2017
//Mathias S. Jørgensen , Michael N. Groves, and Bjørk Hammer
//J. Chem. Theory Comput., 2017, 13 (3), pp 1486–1493
//DOI: 10.1021/acs.jctc.6b01119


	model = strdup("GENERIC");
	if(readOneString(file,"Model",&model)) 
	{
		uppercase(model);
		if(!strstr(model,"GENERIC")) { 
			fprintf(stderr,"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n");
			fprintf(stderr," Sorry, The genetic algorithim works only with Generc model\n"); 
			fprintf(stderr,"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n");
			fflush(stderr); 
			exit(1);
		}
	}

	if(readOneString(file,"CrossType",&tmp) || readOneString(file,"GACrossType",&tmp)) 
	{
		uppercase(tmp);
		if(strstr(tmp,"SPHER")) crossType =  GA_CROSS_SPHERICAL;
	}
	if(readOneString(file,"MutationType",&tmp) || readOneString(file,"GAMutationType",&tmp)) 
	{
		uppercase(tmp);
		if(strstr(tmp,"GLOBAL")) mutationType =  GA_MUTATION_SPHERICAL_GLOBAL;
		if(strstr(tmp,"CM")) mutationType =  GA_MUTATION_SPHERICAL_GLOBAL;
	}

	readOneString(file,"genericCommand",&genericCommand);
	if(!readOneString(file,"QMKeys",&QMKeys)) QMKeys = strdup("NONE");
	readOneInt(file,"popSize",&popSize);
	readOneInt(file,"maxGen",&maxGen);
	readOneInt(file,"maxGeneration",&maxGen);
	readOneInt(file,"nChilds",&nChilds);
	readOneInt(file,"numberOfChilds",&nChilds);
	readOneReal(file,"pCross",&pCross);
	readOneReal(file,"pMut",&pMut);
	readOneReal(file,"mutationRate",&pMut);
	readOneReal(file,"pMutation",&pMut);
	if(readOneBoolean(file,"removeSimilarInertia",&removeSimilarInertia) && removeSimilarInertia) 
		readOneReal(file,"InertiaTol",&inertiaTol);

	readOneBoolean(file,"removeFragmented",&removeFragmented);
	readOneBoolean(file,"removeDissociated",&removeFragmented);
	readOneBoolean(file,"removeSmallDistance",&removeSmallDistance);

	if(readOneBoolean(file,"removeSimilarBonds",&removeSimilarBonds) && removeSimilarBonds)
	{
		readOneReal(file,"sTol",&sTol);
		readOneReal(file,"distMaxTol",&distMaxTol);
	}

	readOneBoolean(file,"RDChain",&chain);
	readOneBoolean(file,"RDSaveFirstGeom",&saveFirstGeom);
	readOneInt(file,"nTimesGeoms",&nTimesGeoms);

	ga->chain = chain ;
	ga->saveFirstGeom =saveFirstGeom ;
	ga-> nTimesGeoms=nTimesGeoms;

	dirName = strdup(getenv("PWD"));

	ga->removeSimilarInertia = removeSimilarInertia;
	ga->inertiaTol = inertiaTol;
	ga->removeFragmented = removeFragmented;
	ga->removeSmallDistance = removeSmallDistance;
	ga->removeSimilarBonds = removeSimilarBonds;
	ga->sTol = sTol;
	ga->distMaxTol = distMaxTol ;

	ga->maxGens = maxGen;
	ga->popSize = popSize;
	// pCross not used, The number of childs is fixed by nChilds (pCross=nChilds/popSize)
	ga->pCross = pCross;
	ga->pMutation = pMut;
	ga->crossType = crossType;
	ga->mutationType = mutationType;
	ga->pop = NULL;
	ga->nChilds = nChilds;

	ga->numGen = 0;
	ga->nGeoms = popSize;
	ga->fmin=0;

	ga->suffixFileName = getSuffixNameFile(ga->inputFileName);

	ga->logfile=logfile;
	ga->fileAllEnergies=NULL;
	ga->fileSelectedEnergies=NULL;

	ga->constraints = constraints;
	ga->command = genericCommand;
	ga->QMKeys = QMKeys;
	ga->dirName = dirName;
}

/*********************************************************************************************************************/
static boolean popFromInputFiles(GeneticAlgorithm* ga)
{

	//fprintf(stderr,"Begin readMolecules\n"); fflush(stderr);
	Molecule** mols = mols = readMolecules(ga->inputFileName,FALSE);
	int popSize=0;
	while(mols && mols[popSize] != NULL) popSize++;
	//fprintf(stderr,"End readMolecules popSize=%d\n",popSize); fflush(stderr);

	if(popSize<2) return FALSE;
	ga->popSize=popSize;

	int i;
	ga->pop = malloc(ga->popSize*sizeof(QuantumMechanicsModel*));
	for (i = 0; i < ga->popSize; i++ )
	{
		//fprintf(stderr,"create mol n %d\n",i+1); fflush(stderr);
		ga->pop[i] = malloc(sizeof(QuantumMechanicsModel));
		*ga->pop[i] = createGenericModel(mols[i], ga->QMKeys, ga->dirName, ga->command, ga->constraints, ga->logfile);
	}
	return TRUE;

}
/*********************************************************************************************************************/
static boolean popRandom(GeneticAlgorithm* ga)
{
	Molecule* mol = readMolecule(ga->inputFileName,TRUE);
	double* energies = NULL;

	QuantumMechanicsModel qmModel = createGenericModel(mol, ga->QMKeys, ga->dirName, ga->command, ga->constraints, ga->logfile);

/*
	setH4Correction(file,&qmModel);
	readOneBoolean(file,"addD3Correction",&qmModel.addD3Correction);
	checkWallCorrection(file, &qmModel);
*/

	int nOld = ga->popSize*ga->nTimesGeoms;
	if(ga->nTimesGeoms>1) ga->popSize = nOld;
	ga->pop = qmModel.klass->getQuantumMechanicsRDConfo(&qmModel, ga->popSize, ga->chain, ga->saveFirstGeom);

	if(ga->removeSimilarInertia) 
		qmModel.klass->removeSimilarInertiaGeometries(ga->pop, &ga->popSize, energies,ga->logfile,ga->inertiaTol);
	if(ga->removeFragmented) 
		qmModel.klass->removeFragmentedMolecules(ga->pop, &ga->popSize, energies, ga->logfile);
	if(ga->removeSmallDistance) 
		qmModel.klass->removeSmallDistanceMolecules(ga->pop, &ga->popSize, energies, ga->logfile);
	if(ga->removeSimilarBonds) 
		qmModel.klass->removeSimilarBondsGeometries(ga->pop, &ga->popSize, energies,ga->logfile,ga->sTol, ga->distMaxTol);

	if(ga->nTimesGeoms>1) qmModel.klass->cutByInertia(ga->pop, &ga->popSize, energies,nOld/ga->nTimesGeoms,ga->logfile);

	return TRUE;

}
/*********************************************************************************************************************/
GeneticAlgorithm* newGeneticAlgorithm(char* inputFileName)
{
	GeneticAlgorithm* ga= malloc(sizeof(GeneticAlgorithm));

	if(inputFileName) ga->inputFileName=strdup(inputFileName);
	else ga->inputFileName=strdup("Pop.ici");

	readParameters(ga);
	if(ga->popSize>1)
	{
		if(!popRandom(ga)) return NULL;
	}
	else
	{
		if(!popFromInputFiles(ga)) return NULL;
	}
	if(!checkSizes(ga)) return NULL;
	initTables(ga);
	setMethods(ga);
	return ga;

}
