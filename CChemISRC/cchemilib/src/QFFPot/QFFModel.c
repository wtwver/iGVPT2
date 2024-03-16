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

/* QFFModel.c */
#include <math.h>
#include <time.h>
#include "../QFFPot/QFFModel.h"
#include "../Utils/Utils.h"
#include "../Utils/AtomsProp.h"
#include "../Utils/Constants.h"
#include "../Utils/QL.h"

//static boolean printMax = FALSE;
//static boolean printMax = TRUE;

static void readData(QFFModel* qffModel, char* geometryFileName, char* qffFileName);
static void computeEnergy(QFFModel* qffModel);
static void computeDipole(QFFModel* qffModel);
static void computeQ(QFFModel* qffModel);
static void computeNumericGradients(QFFModel* qffModel);
static void computeAnalyticGradients(QFFModel* qffModel);
static int computeFrequencies(QFFModel* qffModel);
static void compareFrequencies(QFFModel* qffModel,FILE* file);
static void compute(QFFModel* qffModel);
static void runAlea(QFFModel* qffModel);
static void convertToAU(QFFModel* qffModel);
static void convertToAU2(QFFModel* qffModel);
static void computeQFFParameters(QFFModel* qffModel);
/*****************************************************************************/
static void printModesAndVelocities(QFFModel* qffModel,FILE* file)
{
	int i;
	fprintf(file,"QValues\n");
	for(i=0;i<qffModel->vData.nFrequencies;i++)
		fprintf(file,"%d %0.10e\n", i+1, qffModel->Q[i]);
	fprintf(file,"END\n");
	fprintf(file,"QVelocities\n");
	for(i=0;i<qffModel->vData.nFrequencies;i++)
		fprintf(file,"%d %0.10e\n", i+1, qffModel->velocity[i]);
	fprintf(file,"END\n");
}
/*****************************************************************************/
static double getKineticEnergy(QFFModel* qffModel)
{
	double ekin = 0;
	int i;
	for(i=0;i<qffModel->vData.nFrequencies;i++)
	{
			double mass = 1.0; // Mass-weigted normal coordinates
			//mass = qffModel->vData.effectiveMasses[i]*AMUTOAU;// to be use with convertToAU2
			//mass = qffModel->vData.effectiveMasses[i];
	
			ekin += qffModel->velocity[i]*
			        qffModel->velocity[i]*
			        mass;
	}
	ekin /=2;
	return ekin;
}
/********************************************************************************/
static double getKelvin(QFFModel* qffModel)
{
	int nFree = qffModel->vData.nFrequencies;
	int i;
	for(i=0;i<qffModel->vData.nFrequencies;i++) if(!qffModel->variable[i]) nFree--;
	if(nFree<1) return 0;
	double kin = qffModel->klass->getKineticEnergy(qffModel);
	return 2*kin / ( nFree * KbInAU);
}
/*****************************************************************************/
static void scaleVelocities(QFFModel* qffModel, double temperature)
{
	double kelvin = qffModel->klass->getKelvin(qffModel);
	double scale = 1.0;
	int i;
	if(temperature<=0) return;
	if(kelvin<=0) return;

	scale = sqrt(temperature/kelvin);
#ifdef DEBUG
	printf("temp = %f kelvin = %f scale = %f\n",temperature, kelvin, scale);
#endif
	for(i=0;i<qffModel->vData.nFrequencies;i++)
		if(qffModel->variable[i]) qffModel->velocity[i] *= scale;
}
/*****************************************************************************/
static void setMaxwellVelocities(QFFModel* qffModel, double temperature)
{
	int i;
	for(i=0;i<qffModel->vData.nFrequencies;i++)
	{
        	if(!qffModel->variable[i]) qffModel->velocity[i] = 0.0;
		else
		{
			double mass = 1.0; // Mass-weigted normal coordinates
			double speed;
			//mass = qffModel->vData.effectiveMasses[i]*AMUTOAU;
			//mass = qffModel->vData.effectiveMasses[i];
			//fprintf(stderr,"mass = %f\n",mass);
			speed = maxwel(mass,temperature);
			qffModel->velocity[i]=speed*drandom();
		}
	}
	qffModel->klass->scaleVelocities(qffModel, temperature);
}
/*****************************************************************************/
static boolean setMaxwellVelocitiesIfNull(QFFModel* qffModel, double temperature)
{
	
	double ekin = qffModel->klass->getKineticEnergy(qffModel);

	fprintf(stderr,"Ekin0=%f\n",ekin);
	if(fabs(ekin)>1e-14) return FALSE;
	setMaxwellVelocities(qffModel, temperature);
	return TRUE;
}
/*****************************************************************************/
static void freeQFFModel(QFFModel* qffModel)
{
	int nAtoms = qffModel->molecule.nAtoms;

	qffModel->vData.klass->free(&qffModel->vData);
	qffModel->pData.klass->free(&qffModel->pData);

	freeVectorDouble(&qffModel->frequencies);
	freeVectorDouble(&qffModel->reducedMasses);
	freeVectorDouble(&qffModel->IRIntensities);
	freeVectorDouble(&qffModel->Q);
	freeMatrixDouble(&qffModel->modes, 3*nAtoms);
}
/*****************************************************************************/
static void convertToAU(QFFModel* qffModel)
{
	qffModel->vData.klass->convertToAU(&qffModel->vData);
	qffModel->pData.klass->convertToAU(&qffModel->pData);
}
/*****************************************************************************/
static void convertToAU2(QFFModel* qffModel)
{
	qffModel->vData.klass->convertToAU2(&qffModel->vData);
	qffModel->pData.klass->convertToAU2(&qffModel->pData);
}
/*****************************************************************************/
static void computeQFFParameters(QFFModel* qffModel)
{
	qffModel->vData.klass->computeQFFParameters(&qffModel->vData);
	//qffModel->pData.klass->computeQFFParameters(&qffModel->pData);// not yet implemented
}
/*****************************************************************************/
static void readData(QFFModel* qffModel, char* geometryFileName, char* qffFileName)
{
	FILE* inputFile;
	Molecule* molFit;
	Molecule* molRef;
	//double C[3];
	double u[3][3];

        inputFile = fopen(qffFileName,"rb");
	if(!inputFile)
	{
		fprintf(stderr, "==========================================================\n");
		fprintf(stderr, "Sorry, I cannot open the %s file\n", qffFileName);
		fprintf(stderr, "==========================================================\n");
		exit(1);
	}
	fclose(inputFile);
	if(geometryFileName != NULL)
	{
        inputFile = fopen(geometryFileName,"rb");
	if(!inputFile)
	{
		fprintf(stderr, "==========================================================\n");
		fprintf(stderr, "Sorry, I cannot open the %s file\n", geometryFileName);
		fprintf(stderr, "==========================================================\n");
		exit(1);
	}
	fclose(inputFile);
	}

	qffModel->vData.klass->readData(&qffModel->vData, qffFileName);
	qffModel->pData.klass->readData(&qffModel->pData, qffFileName);
	qffModel->molecule.klass->read(&qffModel->molecule, geometryFileName);
	if(qffModel->Q) free(qffModel->Q);
	qffModel->Q=NULL;
	qffModel->gradQ=NULL;
	qffModel->velocity=NULL;
	qffModel->variable=NULL;
	if(qffModel->vData.nFrequencies>0) qffModel->Q = malloc(qffModel->vData.nFrequencies*sizeof(double));
	if(qffModel->Q) computeQ(qffModel);
	qffModel->diffStep = 0.01;
	qffModel->typeCalcul = 0;
	qffModel->showWarning = FALSE;
	qffModel->numGradients = FALSE;
	qffModel->xyz = FALSE;
	readOneRealFromAFile(geometryFileName,"diffStep",&qffModel->diffStep);
	readOneIntFromAFile(geometryFileName,"typeCalcul",&qffModel->typeCalcul);
	readOneIntFromAFile(geometryFileName,"showWarning",&qffModel->showWarning);
	readOneIntFromAFile(geometryFileName,"numGradients",&qffModel->numGradients);
	readOneIntFromAFile(geometryFileName,"XYZ",&qffModel->xyz);
        qffModel->frequencies = NULL;
        qffModel->modes = NULL;
        qffModel->reducedMasses = NULL;
        qffModel->IRIntensities = NULL;
	if(qffModel->numGradients) qffModel->klass->computeGradients = computeNumericGradients;
	if(qffModel->vData.nFrequencies>0) 
	{
		int i;
		qffModel->gradQ = newVectorDouble(qffModel->vData.nFrequencies);
		initVectorDouble(qffModel->gradQ, qffModel->vData.nFrequencies, 0.0);
		qffModel->velocity = newVectorDouble(qffModel->vData.nFrequencies);
		qffModel->variable = malloc(qffModel->vData.nFrequencies*sizeof(boolean));
		for(i=0;i<qffModel->vData.nFrequencies;i++) qffModel->variable[i] = TRUE;
	}
	if(qffModel->Q) 
	{
		FILE* file = fopen(geometryFileName,"rb");
		readVectorReal(file,"QValues",qffModel->vData.nFrequencies, qffModel->Q);
		if(qffModel->velocity) readVectorReal(file,"QVelocities",qffModel->vData.nFrequencies, qffModel->velocity);
		fclose(file);
	}

	qffModel->klass->calculateGradient = qffModel->klass->computeGradients;

	//fprintf(stderr,"diffStep=%f\n",qffModel->diffStep);


	molRef = &qffModel->molecule;
	molFit = &qffModel->vData.molecule;

	molFit->klass->fit(molFit,molRef,u);
	//qffModel->vData.klass->rotModes(&qffModel->vData,u);

//	printf("End read readData\n");
}
/**********************************************************************/
static void computeQ(QFFModel* qffModel)
{
	int i,k,c;
	for(i=0;i<qffModel->vData.nFrequencies;i++) 
	{
		qffModel->Q[i] = 0.0;
		for(k=0;k<qffModel->molecule.nAtoms;k++)
		{
			for(c=0;c<3;c++)
			{
				double dx = qffModel->molecule.atoms[k].coordinates[c]-qffModel->vData.molecule.atoms[k].coordinates[c];
				//printf("iatom = %d dx = %f mi = %f mEffe = %f \n",k,dx,qffModel->vData.molecule.atoms[k].mass,qffModel->vData.effectiveMasses[i]);
				//printf("mode %d  atom = %d axis %d\n",i,k,c);
				//if(qffModel->vData.modes) printf("v=%f m=%f em=%f\n",qffModel->vData.modes[i][k][c],qffModel->vData.molecule.atoms[k].mass,qffModel->vData.effectiveMasses[i]);
				//else printf("mode=NULL\n");
				
				qffModel->Q[i] += qffModel->vData.MWModes[i][k][c]*dx;
			}
		}
	}
	// coordinates are in Ang in molecule
	for(i=0;i<qffModel->vData.nFrequencies;i++) 
				qffModel->Q[i] *= ANGTOBOHR;

}
/**********************************************************************/
static void computeEnergyXYZ(QFFModel* qffModel)
{
	int i,j,k,l;
	double E = 0;
	double E2 = 0;
	double E3 = 0;
	double E4 = 0;

	if(qffModel->Q) computeQ(qffModel);

	E2 = 0.0;
	for(i=0;i<qffModel->vData.nFrequencies;i++)
		E2 += qffModel->vData.hessian[i][i]*qffModel->Q[i]*qffModel->Q[i];
	E2 /= 2.0;

	E3 = 0.0;
	for(i=0;i<qffModel->vData.nFrequencies;i++)
	for(j=0;j<qffModel->vData.nFrequencies;j++)
	for(k=0;k<qffModel->vData.nFrequencies;k++)
		E3 += qffModel->vData.cubic[i][j][k]*qffModel->Q[i]*qffModel->Q[j]*qffModel->Q[k];
	E3 /= 6.0;

	E4 = 0.0;
	for(i=0;i<qffModel->vData.nFrequencies;i++)
	for(j=0;j<qffModel->vData.nFrequencies;j++)
	for(k=0;k<qffModel->vData.nFrequencies;k++)
	for(l=0;l<qffModel->vData.nFrequencies;l++)
		E4 += qffModel->vData.quartic[i][j][k][l]*qffModel->Q[i]*qffModel->Q[j]*qffModel->Q[k]*qffModel->Q[l];
	E4 /= 24.0;

	/*
	E3*=-1;
	E4*=-1;
	*/
	// TO REMOVE
	//E3 = 0;
        //E4 = 0;

	// in reduced form : use k instead Phi
	/*
	for(i=0;i<qffModel->vData.nFrequencies;i++)
	for(j=0;j<=i;j++)
	for(k=0;k<=j;k++)
		E += qffModel->vData.cubic[i][j][k]*qffModel->Q[i]*qffModel->Q[j]*qffModel->Q[k];

	for(i=0;i<qffModel->vData.nFrequencies;i++)
	for(j=0;j<=i;j++)
	for(k=0;k<=j;k++)
	for(l=0;l<=k;l++)
		E += qffModel->vData.quartic[i][j][k][l]*qffModel->Q[i]*qffModel->Q[j]*qffModel->Q[k]*qffModel->Q[l];
	*/

	//printf("E2=%0.12e E32=%0.12e E4=%0.12e\n",E2,E3,E4);

	E = E2 + E3 + E4;
	
	if(E<0) 
	{
		char symbol[100];
		fprintf(stderr,"%d\n",qffModel->molecule.nAtoms);
		fprintf(stderr,"Error E<0\n");
		for(k=0;k<qffModel->molecule.nAtoms;k++)
		{
			fprintf(stderr,"%s %f %f %f\n",
				symbol,
				qffModel->molecule.atoms[k].coordinates[0],
				qffModel->molecule.atoms[k].coordinates[1],
				qffModel->molecule.atoms[k].coordinates[2]);
			qffModel->molecule.atoms[k].prop = propAtomGet(symbol);
		}
		exit(1);
	}
	qffModel->molecule.potentialEnergy = E;
	//printf("E(cm^-1) = %0.12f\n",E);
}
/**********************************************************************/
/*
static boolean addWallEnergyQFF(QFFModel* qffModel)
{
	if(!qffModel) return FALSE;
	if(!(qffModel->molecule.wall.E0>0)) return FALSE;
	else{
		return TRUE;// TO REMOVE
		double E0 = qffModel->molecule.wall.E0;
		//double srho2 =1/(qffModel->molecule.wall.rho*qffModel->molecule.wall.rho);
		int n =qffModel->molecule.wall.n;
		int i;
		double energy = 0;
		for(i=0;i<qffModel->vData.nFrequencies;i++) 
		{
			double Qmax = 2.0*1/sqrt(sqrt(qffModel->vData.hessian[i][i]));
			double srho2 = 1/Qmax/Qmax;
			//energy += E0*pow(1-exp(-srho2*qffModel->Q[i]*qffModel->Q[i]),n);
			energy += E0*pow(1-exp(-srho2*qffModel->Q[i]*qffModel->Q[i]),n);
			fprintf(stderr,"mode %d Qmax = %f Q = %f\n",i+1, 1/sqrt(sqrt(qffModel->vData.hessian[i][i])),qffModel->Q[i]);
		}

		fprintf(stderr,"EWall=%0.10fEOld = %0.10f\n",energy,qffModel->molecule.potentialEnergy);
		qffModel->molecule.potentialEnergy += energy;
		return TRUE;
	}
}
*/
/**********************************************************************/
/*
static boolean addLEnergyQFF(QFFModel* qffModel)
{
	if(!qffModel) return FALSE;
	else{
		int j;
		for(j=0;j<qffModel->vData.nFrequencies;j++) 
		{
			double Qmax = 1/sqrt(sqrt(qffModel->vData.hessian[j][j]));
			if(qffModel->Q[j]<-Qmax || qffModel->Q[j]>Qmax)
			{
				int i;
				double E1MR = 0.0;
				double Qi;
				int c;
				int n;
				QFFPotParameters* qffPotParameters = qffModel->vData.qffPotParameters;
				int m;
				for(m=0;m<qffPotParameters->numberOf1MR;m++) 
				{
					i = qffPotParameters->qff1MR[m].numbers[0];
					n = qffPotParameters->qff1MR[m].numbers[1];
					Qi=1; for(c=0;c<n;c++) Qi*= qffModel->Q[i];
					E1MR +=  qffPotParameters->qff1MR[m].energy*Qi;
				}
				qffModel->molecule.potentialEnergy = E1MR;
				printf("i= %d Epot = %f Qi = %f Qmax  = %f\n",j,E1MR,qffModel->Q[j],Qmax);
				if(E1MR<0) exit(1);
				return TRUE;
			}
		}
		/ *
		int i;
		int k = 0;
		double energy = 0;
		for(i=0;i<qffModel->vData.nFrequencies;i++) 
		{
			double Qmax = 1/sqrt(sqrt(qffModel->vData.hessian[i][i]));
			QFFPotParameters* qffPotParameters = qffModel->vData.qffPotParameters;
			if(qffModel->Q[i]<-Qmax || qffModel->Q[i]>Qmax)
			{
				double Emax = 0;
				int m;
				int sign = 1;
				double xmin,xmax,x;
				double y, ymax;
				int l;
				k++;
				if(qffModel->Q[i]<-Qmax) sign = -1;
				for(m=0;m<qffPotParameters->numberOf1MR;m++) 
				{
					int j = qffPotParameters->qff1MR[m].numbers[0];
					int n = qffPotParameters->qff1MR[m].numbers[1];
					int c;
					double Qi=1; 
					if(i!=j) continue;
					for(c=0;c<n;c++) Qi*= qffModel->Q[j];
					Emax +=  qffPotParameters->qff1MR[m].energy*Qi;
				}
				printf("Emax = %f\n",Emax);
				ymax = Emax;
				y = 0;
				for(l=0;l<qffModel->vData.nFrequencies;l++) 
				{
					sign = 1;
					if(qffModel->Q[l]<0) sign = -1;
					xmax = Qmax*sign;
					xmin = 0.9*xmax;
					y += ymax*(qffModel->Q[l]-xmin)/(xmax-xmin);
				}
				energy += y;
				printf("i= %d Epot = %f Qi = %f Qmax  = %f\n",i,y,qffModel->Q[i],Qmax);
				if(Emax<0) exit(1);
				break;
			}
		}
		if(k>0) 
		{
			qffModel->molecule.potentialEnergy = energy;
			printf("Epot = %f\n",qffModel->molecule.potentialEnergy);
			return TRUE;
		}
		* /
		return FALSE;
	}
}
*/
/**********************************************************************/
static boolean computeEnergyQFF(QFFModel* qffModel)
{
	int i,j,k,l;
	int m;
	int n;
	int c;
	double Qi;
	double Qj;
	double Qk;
	double E = 0;
	double E1MR = 0;
	double E2MR = 0;
	double E3MR = 0;
	double E4MR = 0;
	//double f = qffModel->molecule.wall.rho;
	QFFPotParameters* qffPotParameters = qffModel->vData.qffPotParameters;
	double* QOld;

	if(!qffPotParameters) return FALSE;
//	if(addLEnergyQFF(qffModel)) return TRUE;

//	if((qffModel->molecule.wall.E0>0)) fprintf(stderr,"f=%0.12f\n",f);

	QOld = malloc(qffModel->vData.nFrequencies*sizeof(double));
	for(i=0;i<qffModel->vData.nFrequencies;i++) 
	{
		double Qmax = 1/sqrt(sqrt(qffModel->vData.hessian[i][i]));
		double srho2 =  qffModel->molecule.wall.rho*1/Qmax/Qmax;
		double f = exp(-srho2*qffModel->Q[i]*qffModel->Q[i]);
		if(srho2*qffModel->Q[i]*qffModel->Q[i]>100) f = 0;
		//printf("i=%d f = %f ex=%f srho2= %f Qmax = %f Qi = %f\n",i,f,srho2*qffModel->Q[i]*qffModel->Q[i],srho2,Qmax,qffModel->Q[i]);
		QOld[i] = qffModel->Q[i];
		qffModel->Q[i] *= f;
	}


	E1MR = 0.0;
	for(m=0;m<qffPotParameters->numberOf1MR;m++) 
	{
		i = qffPotParameters->qff1MR[m].numbers[0];
		n = qffPotParameters->qff1MR[m].numbers[1];
		Qi=1; for(c=0;c<n;c++) Qi*= qffModel->Q[i];
		E1MR +=  qffPotParameters->qff1MR[m].energy*Qi;
	}
	E2MR = 0.0;
	for(m=0;m<qffPotParameters->numberOf2MR;m++) 
	{
		i = qffPotParameters->qff2MR[m].numbers[0];
		j = qffPotParameters->qff2MR[m].numbers[1];
		n = qffPotParameters->qff2MR[m].numbers[2];
		if(n==-2)
		{
			// ii jj
			n = 2;
			Qi=1; for(c=0;c<n;c++) Qi *= qffModel->Q[i];
			Qj = qffModel->Q[j]*qffModel->Q[j];
			E2MR +=  qffPotParameters->qff2MR[m].energy*Qi*Qj;
		}
		else
		{
			Qi=1; for(c=0;c<n;c++) Qi *= qffModel->Q[i];
			Qj = qffModel->Q[j];
			E2MR +=  qffPotParameters->qff2MR[m].energy*Qi*Qj;
		}
	}
	E3MR = 0.0;
	for(m=0;m<qffPotParameters->numberOf3MR;m++) 
	{
		i = qffPotParameters->qff3MR[m].numbers[0];
		j = qffPotParameters->qff3MR[m].numbers[1];
		k = qffPotParameters->qff3MR[m].numbers[2];
		n = qffPotParameters->qff3MR[m].numbers[3];
		Qi=1; for(c=0;c<n;c++) Qi*= qffModel->Q[i];
		Qj = qffModel->Q[j];
		Qk = qffModel->Q[k];
		E3MR +=  qffPotParameters->qff3MR[m].energy*Qi*Qj*Qk;
	}
	E4MR = 0.0;
	for(m=0;m<qffPotParameters->numberOf4MR;m++) 
	{
		i = qffPotParameters->qff4MR[m].numbers[0];
		j = qffPotParameters->qff4MR[m].numbers[1];
		k = qffPotParameters->qff4MR[m].numbers[2];
		l = qffPotParameters->qff4MR[m].numbers[3];
		E4MR +=  qffPotParameters->qff4MR[m].energy*qffModel->Q[i]*qffModel->Q[j]*qffModel->Q[k]*qffModel->Q[l];
	}

	E =  E1MR + E2MR + E3MR + E4MR ;

	for(i=0;i<qffModel->vData.nFrequencies;i++) qffModel->Q[i] = QOld[i];
	free(QOld);

	qffModel->molecule.potentialEnergy = E;
	//addWallEnergyQFF(qffModel);
	//fprintf(stderr,"QFF: E(cm-1)=%f\n",E*AUTOCM1);
	//fprintf(stderr,"Q1=%f\n",qffModel->Q[0]);
	if(qffModel->molecule.potentialEnergy<0) 
	{
		fprintf(stderr,"Error E<0 : E = %f\n",qffModel->molecule.potentialEnergy);
		for(k=0;k<qffModel->vData.nFrequencies;k++)
		{
			fprintf(stderr,"Mode #%d Q=%0.10f\n",
				k+1,
				qffModel->Q[k]); 
		}
		exit(1);
	}
	return TRUE;
}
/**********************************************************************/
static void computeEnergyQ(QFFModel* qffModel)
{
	int i,j,k,l;
	double E = 0;
	double E2 = 0;
	double E3 = 0;
	double E4 = 0;
	if(computeEnergyQFF(qffModel)) return;

	E2 = 0.0;
	for(i=0;i<qffModel->vData.nFrequencies;i++)
		E2 += qffModel->vData.hessian[i][i]*qffModel->Q[i]*qffModel->Q[i];
	E2 /= 2.0;

	E3 = 0.0;
	for(i=0;i<qffModel->vData.nFrequencies;i++)
	for(j=0;j<qffModel->vData.nFrequencies;j++)
	for(k=0;k<qffModel->vData.nFrequencies;k++)
		E3 += qffModel->vData.cubic[i][j][k]*qffModel->Q[i]*qffModel->Q[j]*qffModel->Q[k];
	E3 /= 6.0;

	E4 = 0.0;
	for(i=0;i<qffModel->vData.nFrequencies;i++)
	for(j=0;j<qffModel->vData.nFrequencies;j++)
	for(k=0;k<qffModel->vData.nFrequencies;k++)
	for(l=0;l<qffModel->vData.nFrequencies;l++)
		E4 += qffModel->vData.quartic[i][j][k][l]*qffModel->Q[i]*qffModel->Q[j]*qffModel->Q[k]*qffModel->Q[l];
	E4 /= 24.0;

	E = E2 + E3 + E4;

	//fprintf(stderr,"E(cm-1)=%f\n",E*AUTOCM1);
	//fprintf(stderr,"Q1=%f\n",qffModel->Q[0]);
	if(E<0) 
	{
		fprintf(stderr,"Error E<0\n");
		for(k=0;k<qffModel->vData.nFrequencies;k++)
		{
			fprintf(stderr,"Mode #%d Q=%0.10f\n",
				k+1,
				qffModel->Q[k]); 
		}
		exit(1);
	}
	qffModel->molecule.potentialEnergy = E;
}
/**********************************************************************/
static void computeEnergy(QFFModel* qffModel)
{
	if(qffModel->xyz) computeEnergyXYZ(qffModel);
	else computeEnergyQ(qffModel);
}
/**********************************************************************/
/*
static void invGradientsXYZ(QFFModel* qffModel)
{
	int a,c;
	Molecule* mol = &qffModel->molecule;
        for(a=0;a<mol->nAtoms;a++)
        for(c=0;c<3;c++)
                mol->atoms[a].gradient[c] = - mol->atoms[a].gradient[c];
}
static void invGradientsQ(QFFModel* qffModel)
{
	int i ;
	for(i=0;i<qffModel->vData.nFrequencies;i++)
                qffModel->gradQ[i] = -qffModel->gradQ[i];
}
*/
/**********************************************************************/
static void computeAnalyticGradientsXYZ(QFFModel* qffModel)
{
        int i,j,k,l;
        int a;
	int c;
        Molecule* mol;
        int nAtoms;
	double dE2, dE3, dE4;

        if(!qffModel || qffModel->molecule.nAtoms<1) return;

//	printf("dx = %f\n",dx);
        mol = &qffModel->molecule;
        nAtoms = mol->nAtoms;

	//printf("Analytic\n");
        qffModel->klass->computeEnergy(qffModel);
        qffModel->klass->computeDipole(qffModel);
//	if(qffModel->molecule.potentialEnergy<0)return;

        for(a=0;a<nAtoms;a++)
        for(c=0;c<3;c++)
                mol->atoms[a].gradient[c] = 0.0;

        for(a=0;a<nAtoms;a++)
        for(c=0;c<3;c++)
	{
		dE2 = 0.0;
		for(i=0;i<qffModel->vData.nFrequencies;i++)
			dE2 += qffModel->vData.hessian[i][i]*qffModel->Q[i]*qffModel->vData.MWModes[i][a][c];
                mol->atoms[a].gradient[c] += dE2;
	}

        for(a=0;a<nAtoms;a++)
        for(c=0;c<3;c++)
	{
		dE3 = 0.0;
		for(i=0;i<qffModel->vData.nFrequencies;i++)
		for(j=0;j<qffModel->vData.nFrequencies;j++)
		for(k=0;k<qffModel->vData.nFrequencies;k++)
		{
			dE3 += qffModel->vData.cubic[i][j][k]*(
			qffModel->vData.MWModes[i][a][c]*qffModel->Q[j]*qffModel->Q[k]+
			qffModel->vData.MWModes[j][a][c]*qffModel->Q[k]*qffModel->Q[i]+
			qffModel->vData.MWModes[k][a][c]*qffModel->Q[i]*qffModel->Q[j]
			);
		}
		dE3 /= 6.0;
                mol->atoms[a].gradient[c] += dE3;
	}

        for(a=0;a<nAtoms;a++)
        for(c=0;c<3;c++)
	{
		dE4 = 0.0;
		for(i=0;i<qffModel->vData.nFrequencies;i++)
		for(j=0;j<qffModel->vData.nFrequencies;j++)
		for(k=0;k<qffModel->vData.nFrequencies;k++)
		for(l=0;l<qffModel->vData.nFrequencies;l++)
		{
			dE4 += qffModel->vData.quartic[i][j][k][l]*(
			qffModel->vData.MWModes[i][a][c]*qffModel->Q[j]*qffModel->Q[k]*qffModel->Q[l]+
			qffModel->vData.MWModes[j][a][c]*qffModel->Q[k]*qffModel->Q[l]*qffModel->Q[i]+
			qffModel->vData.MWModes[k][a][c]*qffModel->Q[l]*qffModel->Q[i]*qffModel->Q[j]+
			qffModel->vData.MWModes[l][a][c]*qffModel->Q[i]*qffModel->Q[j]*qffModel->Q[k]
			);
		}
		dE4 /= 24.0;
                mol->atoms[a].gradient[c] += dE4;
	}

        //printf("End computeGradientNumeric\n");
}
/*
static boolean addWallGradientsQFF(QFFModel* qffModel)
{

	if(!qffModel) return FALSE;
	if(!(qffModel->molecule.wall.E0>0)) return FALSE;
	else
	{
		double E0 = qffModel->molecule.wall.E0;
		double srho2 =1/(qffModel->molecule.wall.rho*qffModel->molecule.wall.rho);
		int n =qffModel->molecule.wall.n;
		//double energy = 0;
		int i;
		for(i=0;i<qffModel->vData.nFrequencies;i++) 
		{
			double r2 = qffModel->Q[i]*qffModel->Q[i];
			double ex = 0; 
			double exn = 0;
			ex = exp(-r2*srho2);
			exn = pow(1-ex,n-1);
			//energy += E0*exn*(1-ex);
			exn  = 2*n*srho2*exn*ex*E0;
	       	 	qffModel->gradQ[i] += qffModel->Q[i]*exn;
		}
		//printf("EWall=%f\n",energy);
		//return energy;
		return TRUE;
	}

}
*/
static boolean computeAnalyticGradientsQFF(QFFModel* qffModel)
{
	int i,j,k,l;
	int m;
	int n;
	int c;
	double Qi;
	QFFPotParameters* qffPotParameters = qffModel->vData.qffPotParameters;

	if(!qffPotParameters) return FALSE;
        qffModel->klass->computeEnergy(qffModel);
        qffModel->klass->computeDipole(qffModel);
        for(i=0;i<qffModel->vData.nFrequencies;i++) qffModel->gradQ[i]= 0.0;

	// 1MR
	for(m=0;m<qffPotParameters->numberOf1MR;m++) 
	{
		i = qffPotParameters->qff1MR[m].numbers[0];
		n = qffPotParameters->qff1MR[m].numbers[1];
		Qi=1; for(c=0;c<n-1;c++) Qi*= qffModel->Q[i];
		qffModel->gradQ[i] += qffPotParameters->qff1MR[m].grad*Qi;
	}

	//2MR
	for(m=0;m<qffPotParameters->numberOf2MR;m++) 
	{
		i = qffPotParameters->qff2MR[m].numbers[0];
		j = qffPotParameters->qff2MR[m].numbers[1];
		n = qffPotParameters->qff2MR[m].numbers[2];
		if(n==-2)
		{
			qffModel->gradQ[i] += qffPotParameters->qff2MR[m].grad*qffModel->Q[i]*qffModel->Q[j]*qffModel->Q[j];
			qffModel->gradQ[j] += qffPotParameters->qff2MR[m].grad*qffModel->Q[j]*qffModel->Q[i]*qffModel->Q[i];
		}
		else
		{
			Qi=1; for(c=0;c<n-1;c++) Qi*= qffModel->Q[i];
			qffModel->gradQ[i] += qffPotParameters->qff2MR[m].grad*Qi*qffModel->Q[j];
			Qi=1; for(c=0;c<n;c++) Qi*= qffModel->Q[i];
			qffModel->gradQ[j] += qffPotParameters->qff2MR[m].grad/n*Qi;
		}
	}
	//3MR
	for(m=0;m<qffPotParameters->numberOf3MR;m++) 
	{
		i = qffPotParameters->qff3MR[m].numbers[0];
		j = qffPotParameters->qff3MR[m].numbers[1];
		k = qffPotParameters->qff3MR[m].numbers[2];
		n = qffPotParameters->qff3MR[m].numbers[3];
		Qi=1; for(c=0;c<n-1;c++) Qi*= qffModel->Q[i];
		qffModel->gradQ[i] += qffPotParameters->qff3MR[m].grad*Qi*qffModel->Q[j]*qffModel->Q[k];
	}
	//4MR
	for(m=0;m<qffPotParameters->numberOf4MR;m++) 
	{
		i = qffPotParameters->qff4MR[m].numbers[0];
		j = qffPotParameters->qff4MR[m].numbers[1];
		k = qffPotParameters->qff4MR[m].numbers[2];
		l = qffPotParameters->qff4MR[m].numbers[3];
		qffModel->gradQ[i] += qffPotParameters->qff4MR[m].grad*qffModel->Q[j]*qffModel->Q[k]*qffModel->Q[l];
	}

	//addWallGradientsQFF(qffModel);

	/*
       fprintf(stderr,"End computeAnalyticGradientsQFF\n");
	for(i=0;i<qffModel->vData.nFrequencies;i++)
			fprintf(stderr,"Mode #%d Q=%0.10f grad=%0.10f\n", i+1, qffModel->Q[i],qffModel->gradQ[i]); 
	*/
	return TRUE;
}
/**********************************************************************/
static void computeAnalyticGradientsQ(QFFModel* qffModel)
{
        int i,j,k,l;
	double dE3, dE4;

	if(computeAnalyticGradientsQFF(qffModel)) return;
        if(!qffModel || qffModel->vData.nFrequencies<1) return;

	//printf("dx = %f\n",dx);

	//printf("Analytic\n");
        qffModel->klass->computeEnergy(qffModel);
        qffModel->klass->computeDipole(qffModel);
//	if(qffModel->molecule.potentialEnergy<0) return;

        for(i=0;i<qffModel->vData.nFrequencies;i++) qffModel->gradQ[i]= 0.0;

        for(i=0;i<qffModel->vData.nFrequencies;i++) 
		qffModel->gradQ[i] = qffModel->vData.hessian[i][i]*qffModel->Q[i];

        for(i=0;i<qffModel->vData.nFrequencies;i++) 
	{
		dE3 = 3*qffModel->vData.cubic[i][i][i]*qffModel->Q[i]*qffModel->Q[i];
		for(j=0;j<qffModel->vData.nFrequencies;j++)
			if(i!=j) dE3 += 2*(qffModel->vData.cubic[i][j][i]+qffModel->vData.cubic[i][i][j]+qffModel->vData.cubic[j][i][i])*qffModel->Q[j]*qffModel->Q[i];

		for(j=0;j<qffModel->vData.nFrequencies;j++)
		for(k=0;k<qffModel->vData.nFrequencies;k++)
		{
			if(i==j||i==k) continue;
			dE3 += (
				qffModel->vData.cubic[i][j][k]
				+qffModel->vData.cubic[j][i][k]
				+qffModel->vData.cubic[j][k][i]
				) *qffModel->Q[j]*qffModel->Q[k];
		}
		dE3 /= 6.0;
		qffModel->gradQ[i] += dE3;
	}

	for(i=0;i<qffModel->vData.nFrequencies;i++)
	{
		dE4 = 4*qffModel->vData.quartic[i][i][i][i]*qffModel->Q[i]*qffModel->Q[i]*qffModel->Q[i];
		for(j=0;j<qffModel->vData.nFrequencies;j++)
		if(i!=j) dE4 += 3*(	qffModel->vData.quartic[j][i][i][i]
					+qffModel->vData.quartic[i][j][i][i]
					+qffModel->vData.quartic[i][i][j][i]
					+qffModel->vData.quartic[i][i][i][j]
					)*qffModel->Q[i]*qffModel->Q[i]*qffModel->Q[j];

		for(j=0;j<qffModel->vData.nFrequencies;j++)
		for(k=0;k<qffModel->vData.nFrequencies;k++)
		{
			if(i==j||i==k) continue;
			dE4 += 2*(
					qffModel->vData.quartic[i][i][j][k]+
					qffModel->vData.quartic[i][j][i][k]+
					qffModel->vData.quartic[j][i][i][k]+
					qffModel->vData.quartic[i][j][k][i]+
					qffModel->vData.quartic[j][k][i][i]+
					qffModel->vData.quartic[j][i][k][i]
				   )*qffModel->Q[i]*qffModel->Q[j]*qffModel->Q[k];
		}

		for(j=0;j<qffModel->vData.nFrequencies;j++)
		for(k=0;k<qffModel->vData.nFrequencies;k++)
		for(l=0;l<qffModel->vData.nFrequencies;l++)
		{
			if(i==j||i==k|| i==l) continue;
			dE4 += 
				(
				qffModel->vData.quartic[i][j][k][l]+
				qffModel->vData.quartic[j][i][k][l]+
				qffModel->vData.quartic[j][k][i][l]+
				qffModel->vData.quartic[j][k][l][i]
				)*qffModel->Q[j]*qffModel->Q[k]*qffModel->Q[l];
		}
		dE4 /= 24.0;
		qffModel->gradQ[i] += dE4;
	}

	/*
       fprintf(stderr,"End computeAnalyticGradientsQ\n");
	for(i=0;i<qffModel->vData.nFrequencies;i++)
			fprintf(stderr,"Mode #%d Q=%0.10f grad=%0.10f\n", i+1, qffModel->Q[i],qffModel->gradQ[i]); 
	*/
}
static void computeAnalyticGradients(QFFModel* qffModel)
{

	if(qffModel->xyz) computeAnalyticGradientsXYZ(qffModel);
	else computeAnalyticGradientsQ(qffModel);
	/*
	if(qffModel->molecule.potentialEnergy<0)
	{
		if(qffModel->xyz)invGradientsXYZ(qffModel);
		else invGradientsQ(qffModel);
		qffModel->molecule.potentialEnergy *= -1;
	}
	*/
}
/**********************************************************************/
static void computeNumericGradientsXYZ(QFFModel* qffModel)
{
        int i;
        int k;
        Molecule* mol;
        int nAtoms;
        double Ep, Em;
        double dx;
	double conv = 1/ANGTOBOHR;

        if(!qffModel || qffModel->molecule.nAtoms<1) return;
        qffModel->klass->computeEnergy(qffModel);
        qffModel->klass->computeDipole(qffModel);
//	if(qffModel->molecule.potentialEnergy<0) return;

        dx = qffModel->diffStep;

	//printf("dx = %f\n",dx);
        mol = &qffModel->molecule;
        nAtoms = mol->nAtoms;
	//fprintf(stderr,"computeNumericGradients\n");

        for(i=0;i<nAtoms;i++)
        for(k=0;k<3;k++)
        {
		//printf("i=%d k = %d\n",i,k);
                mol->atoms[i].coordinates[k] += dx;
                qffModel->klass->computeEnergy(qffModel);
                Ep = mol->potentialEnergy;

                mol->atoms[i].coordinates[k] -= 2*dx;
                qffModel->klass->computeEnergy(qffModel);
                Em = mol->potentialEnergy;

                mol->atoms[i].gradient[k] = (Ep-Em)/dx/2*conv;
                mol->atoms[i].coordinates[k] += dx;
        }
        //printf("End computeGradientNumeric\n");
}
/**********************************************************************/
static void computeNumericGradientsQ(QFFModel* qffModel)
{
        int i;
        double Ep, Em;
        double dx;

        if(!qffModel || qffModel->vData.nFrequencies<1) return;
        qffModel->klass->computeEnergy(qffModel);
        qffModel->klass->computeDipole(qffModel);
//	if(qffModel->molecule.potentialEnergy<0) return;

        dx = qffModel->diffStep;

	printf("dx = %f\n",dx);
	//fprintf(stderr,"computeNumericGradientsQ\n");

        for(i=0;i<qffModel->vData.nFrequencies;i++)
        {
		//printf("i=%d k = %d\n",i,k);
		qffModel->Q[i] += dx;
                qffModel->klass->computeEnergy(qffModel);
                Ep = qffModel->molecule.potentialEnergy;

		qffModel->Q[i] -= 2*dx;
                qffModel->klass->computeEnergy(qffModel);
                Em =qffModel->molecule.potentialEnergy;
		qffModel->Q[i] += dx;

		qffModel->gradQ[i] =(Ep-Em)/dx/2;
        }
	/*
       fprintf(stderr,"End computeNumericGradientsQ\n");
	for(i=0;i<qffModel->vData.nFrequencies;i++)
			fprintf(stderr,"Mode #%d Q=%0.10f grad=%0.10f\n", i+1, qffModel->Q[i],qffModel->gradQ[i]); 
	*/
}
static void computeNumericGradients(QFFModel* qffModel)
{

	if(qffModel->xyz)computeNumericGradientsXYZ(qffModel);
	else computeNumericGradientsQ(qffModel);
	/*
	if(qffModel->molecule.potentialEnergy<0)
	{
		if(qffModel->xyz)invGradientsXYZ(qffModel);
		else invGradientsQ(qffModel);
		qffModel->molecule.potentialEnergy *= -1;
	}
	*/
}
/**********************************************************************/
static void computeDipole(QFFModel* qffModel)
{
	int i,j,k;
	int c;
	double dipole[3];
	for(c=0;c<3;c++) dipole[c] = qffModel->pData.zero[c];

	for(i=0;i<qffModel->vData.nFrequencies;i++)
	for(c=0;c<3;c++) dipole[c] += qffModel->pData.first[c][i]*qffModel->Q[i];

	for(i=0;i<qffModel->vData.nFrequencies;i++)
	for(j=0;j<qffModel->vData.nFrequencies;j++)
	for(c=0;c<3;c++) dipole[c] += 0.5*qffModel->pData.second[c][i][j]*qffModel->Q[i]*qffModel->Q[j];

	for(i=0;i<qffModel->vData.nFrequencies;i++)
	for(j=0;j<qffModel->vData.nFrequencies;j++)
	for(k=0;k<qffModel->vData.nFrequencies;k++)
	for(c=0;c<3;c++) dipole[c] += qffModel->pData.cubic[c][i][j][k]*qffModel->Q[i]*qffModel->Q[j]*qffModel->Q[k]/6;

	for(c=0;c<3;c++) qffModel->molecule.dipole[c] = dipole[c];
}
/*****************************************************************************/
static void copyGradients(Molecule* mol, double* g[])
{
	int i,k;
	if(!mol) return;
	for(i=0;i<mol->nAtoms;i++)
		for(k=0;k<3;k++)
			g[k][i] = mol->atoms[i].gradient[k];
}
/*****************************************************************************/
static void copyDipole(Molecule* mol, double d[])
{
	int k;
	if(!mol) return;
	for(k=0;k<3;k++)
		d[k] = mol->dipole[k];
}
/*****************************************************************************/
static void sortFrequencies(int nModes, double* frequencies, double** modes, double* reducedMasses, double* IRIntensities)
{
	int i;
	int j;
	int k;
	double dum;
	if(nModes<1 || !frequencies || !modes || !reducedMasses || !IRIntensities) return;
	for(i=0;i<nModes;i++)
	{
		k = i;
		for(j=i+1;j<nModes;j++)
			if(frequencies[j]<frequencies[k]) k = j;
		if(k==i) continue;
		/* swap i and k modes */
		dum = frequencies[i];
		frequencies[i] = frequencies[k];
		frequencies[k] = dum;
		dum = reducedMasses[i];
		reducedMasses[i] = reducedMasses[k];
		reducedMasses[k] = dum;
		dum = IRIntensities[i];
		IRIntensities[i] = IRIntensities[k];
		IRIntensities[k] = dum;
		for(j=0;j<nModes;j++)
		{
			dum =  modes[j][i];
			modes[j][i] = modes[j][k];
			modes[j][k] = dum;
		}
	}
}
/*****************************************************************************/
static int computeFrequencies(QFFModel* qffModel)
{
	int i;
	int j;
	int k;
	int c;
	int id,jd,index;
	double* F;
	double* gp[3];
	double* gm[3];
	double* dmuX[3];
	Molecule* mol;
	int nAtoms;
	double Dp[3];
	double Dm[3];
	double dx;
	double* frequencies; 
	double** modes; 
	double* reducedMasses; 
	double* IRIntensities;
	/* Intensities in 1 (D/Ang)^2 amu^-1 = 42.255 km/mol=171.65 cm^-2 atm^-1 at 0 C and 1 atm */
	/* Refs : D. Porezag and M. R. Pederson, Phys. Rev. B 54, 7830 (1996). and Y. Yamaguchi el al., J. Chem. Phys. 84,2262(1986)*/
	double AUTokmmolM1 = AUTODEB*AUTODEB*ANGTOBOHR*ANGTOBOHR*AMUTOAU*42.255;

	if(!qffModel || qffModel->molecule.nAtoms<1) return 0;
	dx = qffModel->diffStep;

	//fprintf(stderr,"diffStep of freq =%f\n",qffModel->diffStep);

	//printf("Begin calcul F\n");

	mol = &qffModel->molecule;
	nAtoms = mol->nAtoms;
	for(k=0;k<3;k++) gp[k] = malloc(nAtoms*sizeof(double));
	for(k=0;k<3;k++) gm[k] = malloc(nAtoms*sizeof(double));
	for(k=0;k<3;k++) dmuX[k] = malloc(3*nAtoms*sizeof(double));

	F = malloc(3*nAtoms*(3*nAtoms+1)/2*sizeof(double));

	//printf("End alloc calcul F\n");

	index = 0;
	for(i=0;i<nAtoms;i++)
	for(k=0;k<3;k++)
	{
		id=3*i+k;
		mol->atoms[i].coordinates[k] += dx;
		qffModel->klass->computeGradients(qffModel);
		copyGradients(mol, gp);
		qffModel->klass->computeDipole(qffModel);
		copyDipole(mol, Dp);

		mol->atoms[i].coordinates[k] -= 2*dx;
		qffModel->klass->computeGradients(qffModel);
		copyGradients(mol, gm);
		qffModel->klass->computeDipole(qffModel);
		copyDipole(mol, Dm);
		mol->atoms[i].coordinates[k] += dx;

		for(c = 0;c<3;c++) dmuX[c][id] = (Dp[c]-Dm[c])/dx/2/ANGTOBOHR;
		for(j=0;j<=i;j++)
		{
			double invm = 1.0/sqrt( mol->atoms[i].mass* mol->atoms[j].mass);
			invm/=AMUTOAU;
			for(c = 0;c<3;c++) 
			{
				jd = 3*j+c;
				if(jd>id) continue;
				index = jd + id*(id+1)/2;
				F[index] = (gp[c][j]-gm[c][j])/dx/2/ANGTOBOHR; 
				F[index] *= invm;
			}
		}
	}
	qffModel->klass->computeGradients(qffModel);
	qffModel->klass->computeDipole(qffModel);

	//printf("En calcul F\n");
	for(k=0;k<3;k++) free(gp[k]);
	for(k=0;k<3;k++) free(gm[k]);
	frequencies = malloc(3*nAtoms*sizeof(double));
	reducedMasses = malloc(3*nAtoms*sizeof(double));
	IRIntensities = malloc(3*nAtoms*sizeof(double));
	modes = malloc(3*nAtoms*sizeof(double*));
	for(i=0;i<3*nAtoms;i++) modes[i] = malloc(3*nAtoms*sizeof(double));

	eigenQL(3*nAtoms, F, frequencies, modes);
	free(F);
	/* convert in atomic unit  from kcal/Ang^2/amu */
	//for(i=0;i<3*nAtoms;i++) frequencies[i] *= 1.59360150e-03*0.529177*0.529177*5.48579911e-04; 
	/* convert frequencies in cm-1 */
	for(i=0;i<3*nAtoms;i++) 
		if( frequencies[i]>0) frequencies[i] = sqrt(frequencies[i])*AUTOCM1;
		else frequencies[i] = -sqrt(-frequencies[i])*AUTOCM1;

	/* compute the IR intensities */
	for(i=0;i<nAtoms;i++)
	for(k=0;k<3;k++)
	{
		int id=3*i+k;
		double IRI = 0;
		double D[3] = {0,0,0};
		int kp;
		for(c = 0;c<3;c++)
		for(j=0;j<nAtoms;j++)
		for(kp = 0;kp<3;kp++) 
		{
			int jd = 3*j+kp;
			double Lji = modes[jd][id];
			double a=dmuX[c][jd]*Lji/sqrt(mol->atoms[j].mass*AMUTOAU);
			D[c]+=a;
		}
		IRI = 0;
		for(c = 0;c<3;c++)  IRI+= D[c]*D[c];
		IRIntensities[id] = IRI;
	}
	for(k=0;k<3;k++) free(dmuX[k]);
	/* conversion in km/mol*/

	for(i=0;i<3*nAtoms;i++) IRIntensities[i] *= AUTokmmolM1;

	/* compute the reduced mass */
	for(i=0;i<3*nAtoms;i++) 
	{
		double m = 0;
		for(j=0;j<mol->nAtoms;j++)
		{
			double r2 = 0;
			for(c=0;c<3;c++) r2+= modes[3*j+c][i]*modes[3*j+c][i];
			m+= r2/(mol->atoms[j].mass); 
		}
		if(m<=0) m = 1;
		m = 1/m;
		for(j=0;j<mol->nAtoms;j++)
		{
			double r =sqrt(m)/sqrt(mol->atoms[j].mass);
			for(c=0;c<3;c++) modes[3*j+c][i]*=r;
		}

		reducedMasses[i] = m;
	}
	sortFrequencies(3*nAtoms, frequencies, modes, reducedMasses, IRIntensities);
	qffModel->frequencies = frequencies;
	qffModel->modes = modes;
	qffModel->reducedMasses = reducedMasses;
	qffModel->IRIntensities = IRIntensities;

	return 3*nAtoms;

}
/*****************************************************************************/
static void compareFrequencies(QFFModel* qffModel,FILE* file)
{
	int i,j,k,c;
	double diffFreq;
	double diffMass;
	double overlapModes;
	int nAtoms = qffModel->molecule.nAtoms;
	int ntr = 6;
	double sumDiffFreq=0;
	fprintf(file,"%22s %22s %22s %22s %22s %22s\n",
	"Calc. Freq. [cm^-1]","IR Int(km/mol)",  "Orig. Freq. [cm^-1]", "diffFreq[cm^-1]", "diffMass[amu]", "overlapModes");
	
	if(nAtoms<3) ntr=5;
	for(i=0;i<ntr;i++) fprintf(file,"%10.3lf Rot/Trans\n",qffModel->frequencies[i]);
	for(i=ntr;i<3*nAtoms;i++) 
	{
		j = i-ntr;
		diffFreq = qffModel->frequencies[i]-sqrt(qffModel->vData.hessian[j][j])*AUTOCM1;
		diffMass = qffModel->reducedMasses[i]-qffModel->vData.effectiveMasses[j];
		overlapModes = 0;
		for(k=0;k<nAtoms;k++) 
		for(c=0;c<3;c++) 
			overlapModes += qffModel->modes[3*k+c][i]*qffModel->vData.modes[j][k][c];
		fprintf(file,"%22.3lf ",qffModel->frequencies[i]);
		fprintf(file,"%22.3lf ",qffModel->IRIntensities[i]);
		fprintf(file,"%22.3lf ",sqrt(qffModel->vData.hessian[j][j])*AUTOCM1);
		fprintf(file,"%22.3lf ",diffFreq);
		fprintf(file,"%22.3lf ",diffMass);
		fprintf(file,"%22.3lf ",overlapModes);
		fprintf(file,"\n");
		sumDiffFreq += fabs(diffFreq);
	}
	if(3*nAtoms-ntr>0) fprintf(file,"\naveAbsDiff=%10.3lf\n",sumDiffFreq/(3*nAtoms-ntr));
	fprintf(file,"\n");

	
}
/*****************************************************************************/
static void compute(QFFModel* qffModel)
{
	if(qffModel->typeCalcul==0)
	{
		qffModel->klass->computeEnergy(qffModel);
		qffModel->klass->computeDipole(qffModel);
	}
	else if(qffModel->typeCalcul==1)
		qffModel->klass->computeGradients(qffModel);
	else 
	{
		qffModel->klass->computeFrequencies(qffModel);
		qffModel->klass->compareFrequencies(qffModel,stdout);
	}
}
static void runAlea(QFFModel* qffModel)
{
	int i,j;
	double max=50;
	srand((unsigned int)(time(NULL)));
	for(j=0;j<qffModel->vData.nFrequencies*1000;j++) 
	{
		for(i=0;i<qffModel->vData.nFrequencies;i++) 
			qffModel->Q[i] = 2*max*rand()/(double)RAND_MAX-max;
		qffModel->klass->compute(qffModel);
	}

}
/**********************************************************************/
QFFModel newQFFModel()
{
	QFFModel qffModel;
	qffModel.klass = malloc(sizeof(QFFModelClass));
	qffModel.klass->readData = readData;
	qffModel.klass->convertToAU = convertToAU;
	qffModel.klass->convertToAU2 = convertToAU2;
	qffModel.klass->computeQFFParameters = computeQFFParameters;
	qffModel.klass->free = freeQFFModel;
	qffModel.klass->computeEnergy = computeEnergy;
	qffModel.klass->computeGradients = computeAnalyticGradients;
	qffModel.klass->computeFrequencies = computeFrequencies;
	qffModel.klass->compareFrequencies = compareFrequencies;
	qffModel.klass->compute = compute;
	qffModel.klass->runAlea = runAlea;
	qffModel.klass->computeDipole = computeDipole;
	qffModel.klass->getKineticEnergy = getKineticEnergy;
	qffModel.klass->getKelvin = getKelvin;
	qffModel.klass->scaleVelocities = scaleVelocities;
	qffModel.klass->setMaxwellVelocities = setMaxwellVelocities;
	qffModel.klass->setMaxwellVelocitiesIfNull = setMaxwellVelocitiesIfNull;
	qffModel.klass->printModesAndVelocities = printModesAndVelocities;

	qffModel.vData = newQFFPotentialData(0,0);
	qffModel.pData = newQFFPropertiesData(0,0);
	qffModel.molecule = *(newMolecule());
	qffModel.Q=NULL;
//	printf("End newQFFModel\n");
	return qffModel;
}
/**********************************************************************/
