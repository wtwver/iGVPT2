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

#ifndef __CCHEMILIB_UTILS_H__
#define __CCHEMILIB_UTILS_H__

#include "../Utils/Types.h"
#include "../Molecule/Molecule.h"
#include <stdarg.h>

char* strdup_vprintf(const char* format, va_list ap);
char* strdup_printf(const char* format, ...);

void timing(double* cpu,double *sys);
#ifdef OS_WIN32
void addUnitDisk(FILE* file, CONST char* name);
#endif /* OS_WIN32 */
char* getTimeStr();
boolean isABackspace(char *st);
void waiting(double tsecond);
void debug(char *fmt,...);
char* getLineChars(char c,int n);
char* catFile(char* namefile,boolean tabulation);
CONST char *cchemiDirectory();
char *getSuffixNameFile(const char* allname);
char* getHomeDir();
void uppercase(char *);
void lowercase(char *);
void strfreev (char **str);
char** split(char *str);
void addToPath(char*);
void deleteLastSpaces(char* str);
void deleteFirstSpaces(char* str);
void deleteAllSpaces(char* str);
void strDeleten(char* str);
char* getToStr(char* str,char* end);
boolean isInteger(CONST char *t);
boolean isFloat(CONST char *t);
char * mystrcasestr(CONST char *haystack, CONST char *needle);
boolean readOneReal(FILE* file, char* tag, double*value);
boolean readOneRealFromAFile(char* namefile, char* tag, double* value);
boolean readOneInt(FILE* file, char* tag, int*value);
boolean readOneIntFromAFile(char* namefile, char* tag, int* value);
boolean readOneBoolean(FILE* file, char* tag, boolean*value);
boolean readOneBooleanFromAFile(char* namefile, char* tag, boolean* value);
boolean readOneStringFromAFile(char* namefile, char* tag, int* value);
boolean readOneString(FILE* file, char* tag, char**value);

void userInstallVerify();
void readRessources();
void initAll(int argc, char * argv[]);
void finalize();
double drandom();
double normal();
Molecule* getFixedNormalModeSampling(char* inputFileName, int nModes, double* frequencies, double** modes, double* reducedMasses, double* quantumNumbers);
void printHarmonicVelocities(char* inputFileName, int nModes, double* frequencies, double** modes,double* reducedMasses);
void addHarmonicVelocities(char* inputFileName, int nModes, double* frequencies, double** modes, double* reducedMasses, double* IRIntensities);
double* newVectorDouble(int n);
void freeVectorDouble(double** v);
double** newMatrixDouble(int nrows, int ncolumns);
void freeMatrixDouble(double*** M, int nrows);
void printMatrixDouble(double** M, int nrows, int ncolumns);
void printCubeDouble(double*** C, int nrows, int ncolumns, int nslices);
double*** newCubeDouble(int nrows, int ncolumns, int nslices);
void freeCubeDouble(double**** C, int nrows, int ncolumns);
void initVectorDouble(double* v, int n, double val);
void initMatrixDouble(double** M, int nrows, int ncolumns, double val);
void initCubeDouble(double*** C, int nrows, int ncolumns, int nslices, double val);
int* newVectorInt(int n);
void freeVectorInt(int** v);
int** newMatrixInt(int nrows, int ncolumns);
void freeMatrixInt(int*** M, int nrows);
int*** newCubeInt(int nrows, int ncolumns, int nslices);
void freeCubeInt(int**** C, int nrows, int ncolumns);
void initMatrixInt(int** M, int nrows, int ncolumns, int val);
void initCubeInt(int*** C, int nrows, int ncolumns, int nslices, int val);
int**** newQuarticInt(int nrows, int ncolumns, int nslices, int nl);
void initQuarticInt(int**** C, int nrows, int ncolumns, int nslices, int nl, int val);
void freeQuarticInt(int***** C, int nrows, int ncolumns, int nslices);
double***** newQuinticDouble(int nrows, int ncolumns, int nslices, int nl, int n5);
void initQuinticDouble(double***** C, int nrows, int ncolumns, int nslices, int nl, int n5, double val);
void freeQuinticDouble(double****** C, int nrows, int ncolumns, int nslices, int nl);

void printVectorDoubleCutOff(double* C, int n, double cutoff);

void printMatrixDoubleCutOff(double** M, int nrows, int ncolumns, double cutoff);
void printCubeDouble(double*** C, int nrows, int ncolumns, int nslices);
void printCubeDoubleCutOff(double*** C, int nrows, int ncolumns, int nslices, double cutoff);
void printQuarticDouble(double**** C, int nrows, int ncolumns, int nslices, int nl);
void printQuarticDoubleCutOff(double**** C, int nrows, int ncolumns, int nslices, int nl, double cutoff);
double erfinv( double y );
void getRandVect(double len, double V[]);
double	maxwel(double mass, double temperature);
boolean InverseTensor(double mat[3][3],double invmat[3][3]);
void setMDOptions(FILE* file, int* updateFrequency, 
double* heatTime, double*equiTime, double* runTime, double* coolTime, 
double* heatTemp,  double*equiTemp, double*runTemp, double*coolTemp, 
double* stepSize, MDIntegratorType* integrator, MDThermostatType* thermostat, double* friction, double* omegaMax, int* Nf, double* collide, double* qNH);
boolean readVectorReal(FILE* file, char* tag, int n, double*values);
void computeAngularVelocitiesForALinearMolecule(double* R1, double* R2, double  inert[3][3], double*L, double* vAng);
int readEnergyAndDipoleFromGabeditFile(char* fileName, double energy[], double dipole[]);
double**** newQuarticDouble(int nrows, int ncolumns, int nslices, int nl);
void printQuarticDouble(double**** C, int nrows, int ncolumns, int nslices, int nl);
void printQuarticDoubleCutOff(double**** C, int nrows, int ncolumns, int nslices, int nl, double cutoff);
void initQuarticDouble(double**** C, int nrows, int ncolumns, int nslices, int nl, double val);
void freeQuarticDouble(double***** C, int nrows, int ncolumns, int nslices);
void symmetrizeMatrixDouble(double** M, int nrows, int ncolumns, double cutOff);
void symmetrizeCubeDouble(double*** C, int nrows, int ncolumns,  int nslices, double cutOff);
void symmetrizeQuarticDouble(double**** Q, int nrows, int ncolumns,  int nslices, int nq, double cutOff);
boolean goToStr(FILE* file, char* tag);
void computeFrequenciesFromFilesDlg(char* inputFileName, boolean oneStep);
void computeFrequenciesFromGradFilesDlg(char* inputFileName, boolean oneStep);
void generateCChemIFilesForFrequenciesDlg(char* inputFileName, boolean oneStep);
void generateCChemIGradFilesForFrequenciesDlg(char* inputFileName, boolean oneStep);
int get_num_orbitals_from_aux_mopac_file(FILE* file, char* blockName,  int* begin, int* end);
char** get_one_block_from_aux_mopac_file(FILE* file, char* blockName,  int* n);
char** free_one_string_table(char** table, int n);
void get_dipole_from_mopac_aux_file(FILE* file, double dipole[]);
void get_dipole_from_gamess_output_file(FILE* file, double dipole[]);
void get_dipole_from_turbomole_output_file(FILE* file, double dipole[]);
void get_dipole_from_gaussian_output_file(FILE* file, double dipole[]);
void get_dipole_from_molpro_output_file(FILE* file, double dipole[]);
void get_dipole_from_dalton_output_file(FILE* file, double dipole[]);
void get_dipole_from_orca_output_file(FILE* file, double dipole[]);
void get_dipole_from_vasp_output_file(FILE* file, double dipole[]);
void get_dipole_from_nwchem_output_file(FILE* file, double dipole[]);
void get_dipole_from_psicode_output_file(FILE* file, double dipole[]);
void get_dipole_from_qchem_output_file(FILE* file, double dipole[]);
void get_dipole_from_mopac_output_file(FILE* file, double dipole[]);
void changeDInE(char *st);
int get_one_int_from_fchk_gaussian_file(FILE* file, char* blockName);
double get_one_real_from_fchk_gaussian_file(FILE* file, char* blockName);
int* get_array_int_from_fchk_gaussian_file(FILE* file, char* blockName, int* nElements);
double* get_array_real_from_fchk_gaussian_file(FILE* file, char* blockName, int* nElements);
char** get_array_string_from_fchk_gaussian_file(FILE* file, char* blockName, int* nElements);
double getMaxQuarticIJKL(double**** C, int i, int j, int k, int l);
double getMaxCubicIJK(double*** C, int i, int j, int k);
double getMaxMatrixIJ(double** C, int i, int j);
double getMaxCubicIJ(double*** C, int a, int i, int j);
double getMaxQuarticIJK(double**** C, int a, int i, int j, int k);
boolean readMatrixReal(FILE* file, char* tag, int nrows, int ncolumns, double**values);
boolean readCubeReal(FILE* file, char* tag, int nrows, int ncolumns, int nslices, double***values);
boolean readQuarticReal(FILE* file, char* tag, int nrows, int ncolumns, int nslices, int nq, double****values);
void getCrossProduct(double * V1,double * V2, double *PV);
double getDistancePoints(double* P1,double* P2);
double getNorm(double* V);
double getDotProduct(double* V1,double* V2);
double getAngleVectors(double* V1,double* V2);
double get_wall_time();
double get_cpu_time();
void get3DRandMatrix(double M[3][3]);
void getRandDirection(double D[]);

#endif /* __CCHEMILIB_UTILS_H__ */

