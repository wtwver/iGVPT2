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

#ifndef __CCHEMILIB_FORCEFIELD_H__
#define __CCHEMILIB_FORCEFIELD_H__

#include <stdio.h>
#include "../Molecule/Molecule.h"
#include "../EmpriricalCorrections/HydrogenBondCorrection.h"

#define STRETCHDIM 		9 /* type(=0) a1 a2 Force Re  h3 h4 h5 h6 where hk=anharomonic, or type(=1) a1 a2 alpha Re De if Morse*/
#define BENDDIM 		9 /* a1 a2 a3 Force Theta0 h3 h4 h5 h6 where hk=anharomonic*/
#define DIHEDRALDIM 		8 /*a1 a2 a3 a4 Idiv Pk Phase Pn */
#define IMPROPERDIHEDRALDIM	7 /*a1 a2 a3 a4 Pk Phase Pn */
#define OUTOFPLANEDIM	        10 /*a1 a2 a3 a4 type K h3 h4 h5 h6 where type = 0 (ALLINGER) or 1(Wilson-Decius-Cross) hk = cubic, quartic, .... */
#define VDW612DIM 		4 /* a1 a2 Aij Bij */
#define VDW714DIM 		6 /* a1 a2 epsilon R0 gamma delta */
#define COULOMBDIM 		3 /* a1 a2 CoulombFactor */
#define HYDROGENBONDEDDIM 	6 /* a1 a2 ax Cij Dij if 6-12 or 10-12 or a1 a2 ax alpha Re De if Morse */
#define SUTTONCHENDIM	        7 /* a1 a2 epsilon a C n m*/

#define STRBENDDIM 		12 /* type1 type2 a1 a2 a3 Force12 Force23 Re12 Re23 Theta123 (De12, De23 if Morse), type = 0(harmonic) or 1(Morse) */

#define PAIRWISEDIM 	9 /* a1 a2 A Beta C4 C6 C8 C10 b  : 
			     potential = A*exp(-Beta*r)-Somme C2n*f2n/r**(2*n) + Zi*Zj/r 
			     f2n = 1 - exp(-b*r)*somme((b*r)**k/k!,k=0..2*n)
			   */


typedef struct _ForceField  ForceField;
typedef struct _ForceFieldClass  ForceFieldClass;

typedef enum
{
  AMBER,
  UFF,
  PAIRWISE
} ForceFieldTypes;

typedef struct _ForceFieldOptions
{
	ForceFieldTypes type;
	
	boolean bondStretch;/* For Amber */
	boolean angleBend;/* For Amber */
	boolean strBend;/* For Amber */ // not implemented on GPU
	boolean dihedralAngle;/* For Amber */
	boolean improperTorsion;/* For Amber */
	boolean outOfPlane;/* For Amber */ // not implemented on GPU
	boolean hydrogenBonded612;/* For Amber */
	boolean hydrogenBonded1012;/* For Amber */
	boolean hydrogenBondedMorse;/* For Amber */
	boolean coulomb; /* For Amber and Pair-Wise */
	boolean vdw612;
	boolean vdw714;
	boolean hbDirectional;
	boolean suttonChen;/* For Amber */
	char chargesType[50];
	Constraints rattleConstraints;
	boolean numeric;
	boolean addWallCorrection;
	double dx;
}ForceFieldOptions;

struct _ForceField
{
	Molecule molecule;
	ForceFieldClass* klass;
	int numberOfStretchTerms;
	int numberOfBendTerms;
	int numberOfStrBendTerms;
	int numberOfDihedralTerms;
	int numberOfImproperTorsionTerms;
	int numberOfOutOfPlaneTerms;
	int numberOfVdw612;
	int numberOfVdw714;
	int numberOfCoulomb;
	int numberOfHydrogenBonded;
	int numberOfSuttonChen;
	int numberOfRattleConstraintsTerms;

	double* bondStretchTerms[STRETCHDIM];
	double* angleBendTerms[BENDDIM];
	double* strBendTerms[STRBENDDIM];
	double* dihedralAngleTerms[DIHEDRALDIM];
	double* improperTorsionTerms[IMPROPERDIHEDRALDIM];
	double* outOfPlaneTerms[OUTOFPLANEDIM];
	double* vdw612Terms[VDW612DIM];
	double* vdw714Terms[VDW714DIM];
	double* coulombTerms[COULOMBDIM];
	double* hydrogenBondedTerms[HYDROGENBONDEDDIM];
	double* suttonChenTerms[SUTTONCHENDIM];


	int numberOfPairWise;
	double* pairWiseTerms[PAIRWISEDIM];

	HyhrogenBondCorrectionParameters* H4Parameters;

	FILE* logfile;

	ForceFieldOptions options;
#ifdef ENABLE_CL
	cl_program programMM;

	cl_kernel initEnergyKernel;
	cl_kernel addEnergyBondAmberKernel;
	cl_kernel addEnergyBendAmberKernel;
	cl_kernel addEnergyDihedralAmberKernel;
	cl_kernel addEnergyImproperTorsionAmberKernel;
	cl_kernel addEnergyVdw612AmberKernel;
	cl_kernel addEnergyVdw714AmberKernel;
	cl_kernel addEnergyHydrogenBondedAmberKernel;
	cl_kernel addEnergySuttonChenAmberKernel;
	cl_kernel addEnergyPairWiseKernel;
	cl_kernel initGradientsKernel;
	cl_kernel initVelocitiesKernel;
	cl_kernel addGradientBondAmberKernel;
	cl_kernel addGradientBendAmberKernel;
	cl_kernel addGradientDihedralAmberKernel;
	cl_kernel addGradientImproperTorsionKernel;
	cl_kernel addGradientVdw612AmberKernel;
	cl_kernel addGradientVdw714AmberKernel;
	cl_kernel addGradientHydrogenBondedAmberKernel;
	cl_kernel addGradientSuttonChenAmberKernel;
	cl_kernel addGradientPairWiseKernel;
	cl_kernel reduceGradientsKernel;
	cl_kernel reduceRhoKernel;
	cl_kernel initRhoKernel;
	cl_kernel computeRhoKernel;

	cl_mem energyBufferCL;
	cl_mem gradientBufferCL;
	cl_mem rhoBufferCL;
	cl_mem rattledeltaPositionBufferCL;
	cl_mem rattledeltaVelocityBufferCL;
	cl_mem rattleMovedBufferCL;
	cl_mem rattleUpdateBufferCL;
	cl_mem rattleDoneBufferCL;
	cl_mem atomsCL; /* GPU mem */

	cl_mem bondIndexCL;
	cl_mem bendIndexCL;
	cl_mem dihedralIndexCL;
	cl_mem improperTorsionIndexCL;
	cl_mem nonBondedIndexCL;
	cl_mem  hydrogenBondedIndexCL;
	cl_mem  pairWiseIndexCL;
	cl_mem  rattleConstraintsIndexCL;

	cl_mem bondTermsCL;
	cl_mem bendTermsCL;
	cl_mem dihedralTermsCL;
	cl_mem improperTorsionTermsCL;
	cl_mem vdw612TermsCL;
	cl_mem vdw714TermsCL;
	cl_mem  hydrogenBondedTermsCL;
	cl_mem  pairWiseTermsCL;
	cl_mem rattleConstraintsTermsCL;

	
	int nMaxTerms;
	cl_float* energyBufferCPU; /* CPU */
	cl_float8* atomsCPU; /* CPU mem */

	int nBlockGradientBuffer;
	int nBlockRhoBuffer;
	cl_float4* gradientBufferCPU; /* CPU mem */

	int nBlockRattleBuffer;

	/* cl_float4 pos => pos.s[0], pos.s[1], ... */
	
	
#endif
};
struct _ForceFieldClass
{
	void (*calculateGradient)(ForceField* forceField);
	void (*calculateGradientNumeric)(ForceField* forceField);
	void (*calculateEnergy)(ForceField* forceField);
	double (*calculateEnergyTmp)(ForceField* forceField,Molecule* m);
	void (*printEnergies)(ForceField* forceField);
};


ForceField newForceField();
void freeForceField(ForceField* forceField);
ForceField copyForceField(ForceField* forceField);
void setForceFieldOptions(FILE* file, ForceFieldOptions* forceFieldOptions);
void updateGeometryCL(ForceField* forceField, Molecule* mol);
void updateGeometryVelocitiesCL(ForceField* forceField, Molecule* mol);
void getGeometryVelocitiesCL(ForceField* forceField, Molecule* mol);
int computeMMFrequencies(ForceField* forceField, double** frequencies, double*** modes, double** reducedMasses, double** IRIntensities);
#endif /* __CCHEMILIB_FORCEFIELD_H__ */

