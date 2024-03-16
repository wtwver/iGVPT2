/* HydrogenBondCorrection.c */
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

/* Reference: J. Rezac, P. Hobza J. Chem. Theory Comput. 8, 141-151 (2012)*/

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

#include "HydrogenBondCorrection.h"
#include "../Utils/Constants.h"
#include "../Utils/Utils.h"

static int setDefParameters(HyhrogenBondCorrectionParameters* parameters, char* method);
static void printParameters(HyhrogenBondCorrectionParameters* parameters, FILE* file);
/**********************************************************************/
static void printDefParameters(char* method, FILE* file)
{
	HyhrogenBondCorrectionParameters pars;
	setDefParameters(&pars, method);
	printParameters(&pars, file);
}
/**********************************************************************/
static int setDefParameters(HyhrogenBondCorrectionParameters* parameters, char* method)
{
	sprintf(parameters->method,"%s",method);
/* cutoffs */
	parameters->HB_R_CUTOFF = 5.5;
	parameters->HB_R_0 = 1.5;
	parameters->MAX_XH_BOND=1.15;

	if(strstr(method,"PM6"))
	{
/* H4 repulsion*/
		parameters->para_OH_O = 2.32;
		parameters->para_OH_N = 3.10;
		parameters->para_NH_O = 1.07;
		parameters->para_NH_N = 2.01;
		parameters->multiplier_WH_O = 0.42;
		parameters->multiplier_COO = 1.41;
		parameters->multiplier_NH4 = 3.61;
/* HH repulsion*/
		parameters->HH_REPULSION_k = 0.4;
		parameters->HH_REPULSION_e = 12.7;
		parameters->HH_REPULSION_r0 = 2.3;
	}
	else if(strstr(method,"SCC-DFTB"))
	{
/* H4 repulsion*/
		parameters->para_OH_O = 1.11;
		parameters->para_OH_N = 2.58;
		parameters->para_NH_O = 0.80;
		parameters->para_NH_N = 2.01;
		parameters->multiplier_WH_O = 1.32;
		parameters->multiplier_COO = 1.22;
		parameters->multiplier_NH4 = 2.33;
/* HH repulsion*/
		parameters->HH_REPULSION_k = 0.4;
		parameters->HH_REPULSION_e = 12.7;
		parameters->HH_REPULSION_r0 = 2.3;
	}
	else if(strstr(method,"RM1"))
	{
/* H4 repulsion*/
		parameters->para_OH_O = 3.76;
		parameters->para_OH_N = 3.90;
		parameters->para_NH_O = 3.14;
		parameters->para_NH_N = 2.95;
		parameters->multiplier_WH_O = 0.94;
		parameters->multiplier_COO = 1.10;
		parameters->multiplier_NH4 = 1.21;
/* HH repulsion*/
		parameters->HH_REPULSION_k = 0.4;
		parameters->HH_REPULSION_e = 12.7;
		parameters->HH_REPULSION_r0 = 2.3;
	}
	else if(strstr(method,"OM3"))
	{
/* H4 repulsion*/
		parameters->para_OH_O = 1.95;
		parameters->para_OH_N = 1.64;
		parameters->para_NH_O = 0.93;
		parameters->para_NH_N = 1.35;
		parameters->multiplier_WH_O = 0.50;
		parameters->multiplier_COO = 1.63;
		parameters->multiplier_NH4 = 0.9;
/* HH repulsion*/
		parameters->HH_REPULSION_k = 0.4;
		parameters->HH_REPULSION_e = 12.7;
		parameters->HH_REPULSION_r0 = 2.3;
	}
	else if(strstr(method,"AM1"))
	{
/* H4 repulsion*/
		parameters->para_OH_O = 4.89;
		parameters->para_OH_N = 6.23;
		parameters->para_NH_O = 2.54;
		parameters->para_NH_N = 4.56;
		parameters->multiplier_WH_O = 0.49;
		parameters->multiplier_COO = 1.08;
		parameters->multiplier_NH4 = 2.78;
/* HH repulsion*/
		parameters->HH_REPULSION_k = 0.4;
		parameters->HH_REPULSION_e = 12.7;
		parameters->HH_REPULSION_r0 = 2.3;
	}
	else if(strstr(method,"PM3"))
	{
/* H4 repulsion*/
		parameters->para_OH_O = 2.71;
		parameters->para_OH_N = 4.37;
		parameters->para_NH_O = 2.29;
		parameters->para_NH_N = 3.86;
		parameters->multiplier_WH_O = 0.91;
		parameters->multiplier_COO = 0.89;
		parameters->multiplier_NH4 = 2.54;
/* HH repulsion*/
		parameters->HH_REPULSION_k = 0.4;
		parameters->HH_REPULSION_e = 12.7;
		parameters->HH_REPULSION_r0 = 2.3;
	}
	else if(strstr(method,"MM"))
	{
	// copie from PM3
/* H4 repulsion*/
		parameters->para_OH_O = 2.71;
		parameters->para_OH_N = 4.37;
		parameters->para_NH_O = 2.29;
		parameters->para_NH_N = 3.86;
		parameters->multiplier_WH_O = 0.91;
		parameters->multiplier_COO = 0.89;
		parameters->multiplier_NH4 = 2.54;
/* HH repulsion*/
		parameters->HH_REPULSION_k = 0.4;
		parameters->HH_REPULSION_e = 12.7;
		parameters->HH_REPULSION_r0 = 2.3;
	}
	else
	{
		fprintf(stderr,"I cannot set HB parameters, method %s is unknown\n"
	 	"The known methods are : SCC-DFTB, OM3, PM6, AM1, RM1, PM3\n"
		,method
		);
		fprintf(stderr,"You can also give your parameters using an input file\n"
	 	"Required parameters : OH_O, OH_N, NH_O, NH_N, WH_O, COO, NH4\n"
		"Format of file \n"
		"OH_O= value1\n"
		"OH_N= value2\n"
		"NH_O= value3\n"
		"NH_N= value4\n"
		"WH_O= value5\n"
		"COO= value5\n"
		"NH4= value7\n"
		);
		fprintf(stderr,"Here are the default parameters for several method\n");
		printDefParameters("SCC-DFTB", stderr);
		printDefParameters("OM3", stderr);
		printDefParameters("PM6", stderr);
		printDefParameters("AM1", stderr);
		printDefParameters("RM1", stderr);
		printDefParameters("PM3", stderr);
		exit(1);
		return 1;
	}
	printParameters(parameters, stdout);
	return 0;
}
/****************************************************************************************************/
/* H-H repulsion calculation*/
static double getHHRep(HyhrogenBondCorrectionParameters* parameters, Molecule* molecule, boolean addGradient)
{
	double e_corr_sum = 0; 
	int i, j, k;
	double r;
	double d_rad;
	double g[3];

	// Iterate over H atoms twice
	for(i=0;i<molecule->nAtoms;i++)
	{
		if( molecule->atoms[i].prop.atomicNumber == HYDROGEN)
		for (j = 0; j < i; j++)
		{ 
			if ( molecule->atoms[j].prop.atomicNumber == HYDROGEN)
			{
				r = getDistance(&molecule->atoms[i],&molecule->atoms[j]);
				e_corr_sum += parameters->HH_REPULSION_k * (1.0 - 1.0/(1.0 + exp(-parameters->HH_REPULSION_e *(r/parameters->HH_REPULSION_r0 - 1.0))));

			if (addGradient)
				{
				d_rad = (1.0 / pow(1.0 + exp(- parameters->HH_REPULSION_e*(r/ parameters->HH_REPULSION_r0-1.0)), 2) *  parameters->HH_REPULSION_e/ parameters->HH_REPULSION_r0 * exp(- parameters->HH_REPULSION_e*(r/ parameters->HH_REPULSION_r0-1.0))) *  parameters->HH_REPULSION_k;

				for(k=0;k<3;k++) g[k] =  (molecule->atoms[i].coordinates[k]-molecule->atoms[j].coordinates[k])/ r* d_rad;

				for(k=0;k<3;k++) molecule->atoms[i].gradient[k] -= g[k];
				for(k=0;k<3;k++) molecule->atoms[j].gradient[k] += g[k];
				}
			}
		}
	}
	return e_corr_sum;
}
/*********************************************************************************************************************/
static double getCovalentRadii(int z)
{
// Covalent radii (indexed by proton number), zero -> not available
static const double covalent_radii[119] = {0.0, 0.37, 0.32, 1.34, 0.9, 0.82, 0.77, 
	0.75, 0.73, 0.71, 0.69, 1.54, 1.3, 1.18, 1.11, 1.06, 1.02, 0.99, 0.97,
	1.96, 1.74, 1.44, 1.36, 1.25, 1.27, 1.39, 1.25, 1.26, 1.21, 1.38, 1.31,
	1.26, 1.22, 1.19, 1.16, 1.14, 1.1, 2.11, 1.92, 1.62, 1.48, 1.37, 1.45,
	1.56, 1.26, 1.35, 1.31, 1.53, 1.48, 1.44, 1.41, 1.38, 1.35, 1.33, 1.3,
	2.25, 1.98, 1.69, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
	0.0, 0.0, 1.6, 1.5, 1.38, 1.46, 1.59, 1.28, 1.37, 1.28, 1.44, 1.49, 0.0,
	0.0, 1.46, 0.0, 0.0, 1.45, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
	0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
	0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
	return covalent_radii[z];

}
/*********************************************************************************************************************/
/* Continuous valence contribution of a pair of atoms*/
static double cvalence_contribution(Atom* a, Atom* b)
{
	double r;
	double ri, rj;
	double r0, r1;
	double x;
	int ia = a->prop.atomicNumber;
	int ib = b->prop.atomicNumber;
	ri = getCovalentRadii(ia);
	rj = getCovalentRadii(ib);
	r0 = ri + rj;
	r1 = r0 * 1.6;
	r = getDistance(a,b);
	if (r == 0.0) return 0.0;
	if (r >= r1) return 0.0;
	if (r <= r0) return 1.0;
	x = (r - r0) / (r1 - r0);
	return 1.0 - (-20.0*pow(x,7) + 70.0*pow(x,6) -84.0*pow(x,5) + 35.0*pow(x,4));
}

/*********************************************************************************************************************/
/* Continuous valence contribution of a pair of atoms - derivative in the internal COOrdinate */
static double cvalence_contribution_d(Atom* a, Atom* b)
{
	double r;
	double ri, rj;
	double r0, r1;
	double x;
	int ia = a->prop.atomicNumber;
	int ib = b->prop.atomicNumber;
	ri = getCovalentRadii(ia);
	rj = getCovalentRadii(ib);
	r0 = ri + rj;
	r1 = r0 * 1.6;
	r = getDistance(a,b);
	if (r == 0.0) return 0.0;
	if (r >= r1) return 0.0;
	if (r <= r0) return 0.0;
	x = (r - r0) / (r1 - r0);
	return -(-140.0*pow(x,6) + 420.0*pow(x,5) - 420.0*pow(x,4) + 140.0*pow(x,3)) / (r1 - r0);
}
/*********************************************************************************************************************/
/* addVectorToGrad */
void addVectorToGrad(Molecule* molecule, int i, double add[], double factor)
{
	int c;
	for(c=0;c<3;c++) molecule->atoms[i].gradient[c] += add[c]*factor;
}
/**********************************************************************************************************/
/* H4 correction calculation*/
static double getH4(HyhrogenBondCorrectionParameters* parameters, Molecule* molecule, boolean addGradient)
{
	double e_corr_sum = 0;

	/* H-bond description*/
	int d_i, a_i, h_i; /* donor, acceptor, hydrogen indices*/
	double rda; /* donor-acceptor distance*/
	double rdh, rah;
	double angle;


	/* Energy terms:*/
	double e_para = 0;
	double e_bond_switch;
	double e_radial;
	double e_angular;
	double e_scale_w;
	double e_scale_chd;
	double e_scale_cha;
	double e_corr; 

	/* Derivatives*/
	double d_radial;
	double d_radial_d[3];
	double d_radial_a[3];

	double d_angular;
	double d_angular_d[3];
	double d_angular_h[3];
	double d_angular_a[3];

	double d_bs;
	double d_bs_d[3];
	double d_bs_a[3];
	double d_bs_h[3];

	double g[3];

	/* Scaling derivatives*/
	double sign_wat = 1;

	int o1;
	int o2;
	int cc;

	double cv_O1 = 0;
	double cv_O2 = 0;
	double cv_cc = 0;

	double f_O1;
	double f_O2;
	double f_cc;

	int i, j, k, c; 
	double rih, rjh;
	double x, xd, xd2, a, d;
	double slope, v, fv, fv2;
	double rdhs, ravgs;
	double factor;

	for(i=0;i<molecule->nAtoms;i++)
	{
		if( molecule->atoms[i].prop.atomicNumber == NITROGEN || molecule->atoms[i].prop.atomicNumber == OXYGEN)
		{
		for (j = 0; j < i; j++)
		{ 
		if( molecule->atoms[j].prop.atomicNumber == NITROGEN || molecule->atoms[j].prop.atomicNumber == OXYGEN)
		{
			// Calculate donor-acceptor distance
			rda = getDistance(&molecule->atoms[i],&molecule->atoms[j]);
			// Continue only WHen in range WHere correction acts
			if (rda > parameters->HB_R_0 && rda < parameters->HB_R_CUTOFF)
			{
				// Iterate over hydrogens
				for(h_i=0;h_i<molecule->nAtoms;h_i++)
				{ 
					if ( molecule->atoms[h_i].prop.atomicNumber == HYDROGEN)
					{
					// Distances to hydrogen
					rih = getDistance(&molecule->atoms[i],&molecule->atoms[h_i]);

					rjh = getDistance(&molecule->atoms[j],&molecule->atoms[h_i]);
					angle = M_PI - getAngle(&molecule->atoms[i], &molecule->atoms[h_i], &molecule->atoms[j])*(1/RADTODEG);
					if (angle < M_PI/2) {
						// filterd out everything but corrected H-bonds
						// Determine donor and acceptor - donor is the closer one
						if (rih <= rjh)
						{
							d_i = i;
							a_i = j;
							rdh = rih;
							rah = rjh;
						}
						else
						{
							d_i = j;
							a_i = i;
							rdh = rjh;
							rah = rih;
						}

						// Radial term
						e_radial = -0.00303407407407313510 * pow(rda,7) +
						            0.07357629629627092382 * pow(rda,6) +
						           -0.70087111111082800452 * pow(rda,5) +
						            3.25309629629461749545 * pow(rda,4) +
						           -7.20687407406838786983 * pow(rda,3) +
						            5.31754666665572184314 * pow(rda,2) +
						            3.40736000001102778967 * rda +
						           -4.68512000000450434811;

						// Radial gradient
						if (addGradient)
						{
							// In rDA COOrdinate
							d_radial = -0.02123851851851194655 * pow(rda,6) +
							            0.44145777777762551519 * pow(rda,5) +
							           -3.50435555555413991158 * pow(rda,4) +
							           13.01238518517846998179 * pow(rda,3) +
							          -21.62062222220516360949 * pow(rda,2) +
							           10.63509333331144368628 * rda +
							            3.40736000001102778967;

							// Cartesian gradients on D and A atoms
							for(k=0;k<3;k++)
								d_radial_d[k] = (molecule->atoms[d_i].coordinates[k] - molecule->atoms[a_i].coordinates[k])/rda * d_radial;

							for(k=0;k<3;k++) d_radial_a[k] = -d_radial_d[k];
						}

						// Angular term
						a = angle/(M_PI/2.0);
						x = -20.0*pow(a,7) + 70.0*pow(a,6) - 84.0*pow(a,5) + 35.0*pow(a,4);
						e_angular = 1.0 - x*x;

						// Angular gradient
						if (addGradient)
						{
							xd = (-140.0*pow(a,6) + 420.0*pow(a,5) - 420.0*pow(a,4) + 140.0*pow(a,3)) / (M_PI/2.0);
							d_angular = -xd * 2.0 * x;

							// Dot product of bond vectors
							d = (molecule->atoms[d_i].coordinates[0] - molecule->atoms[h_i].coordinates[0])
							    *(molecule->atoms[a_i].coordinates[0] - molecule->atoms[h_i].coordinates[0]) + 
							    (molecule->atoms[d_i].coordinates[1] - molecule->atoms[h_i].coordinates[1])
							    *(molecule->atoms[a_i].coordinates[1] - molecule->atoms[h_i].coordinates[1]) + 
							    (molecule->atoms[d_i].coordinates[2] - molecule->atoms[h_i].coordinates[2])
							    *(molecule->atoms[a_i].coordinates[2] - molecule->atoms[h_i].coordinates[2]);

							x = -d_angular / sqrt(1.0 - (d*d) / (rdh*rdh) / (rah*rah));
							
							// Donor atom
							for(k=0;k<3;k++)
							d_angular_d[k] = x * -((molecule->atoms[a_i].coordinates[k] - molecule->atoms[h_i].coordinates[k])/rdh/rah 
 									- (molecule->atoms[d_i].coordinates[k] - molecule->atoms[h_i].coordinates[k])*d/pow(rdh,3)/rah);
							// Acceptor atom
							for(k=0;k<3;k++)
							d_angular_a[k] = x * -((molecule->atoms[d_i].coordinates[k] - molecule->atoms[h_i].coordinates[k])/rdh/rah 
 									- (molecule->atoms[a_i].coordinates[k] - molecule->atoms[h_i].coordinates[k])*d/pow(rah,3)/rdh);

							// Hydrogen
							for(k=0;k<3;k++)
							d_angular_h[k] = -d_angular_d[k] - d_angular_a[k];
						}

						// Energy coefficient
						if ( molecule->atoms[d_i].prop.atomicNumber == OXYGEN &&  molecule->atoms[a_i].prop.atomicNumber == OXYGEN)     
														e_para = parameters->para_OH_O;
						if ( molecule->atoms[d_i].prop.atomicNumber == OXYGEN &&  molecule->atoms[a_i].prop.atomicNumber == NITROGEN)     
														e_para = parameters->para_OH_N;
						if ( molecule->atoms[d_i].prop.atomicNumber == NITROGEN &&  molecule->atoms[a_i].prop.atomicNumber == OXYGEN)     
														e_para = parameters->para_NH_O;
						if ( molecule->atoms[d_i].prop.atomicNumber == NITROGEN &&  molecule->atoms[a_i].prop.atomicNumber == NITROGEN)     
														e_para = parameters->para_NH_N;
						// Bond switching
						if (rdh > 1.15)
						{
							rdhs = rdh - 1.15;
							ravgs = 0.5*rdh + 0.5*rah - 1.15;
							x = rdhs/ravgs;
							e_bond_switch = 1.0-(-20.0*pow(x,7) + 70.0*pow(x,6) - 84.0*pow(x,5) + 35.0*pow(x,4));

							// Gradient
							if (addGradient)
							{
								d_bs = -(-140.0*pow(x,6) + 420.0*pow(x,5) - 420.0*pow(x,4) + 140.0*pow(x,3));

								xd = d_bs / ravgs;
								xd2 = 0.5 * d_bs * -x / ravgs;

								for(k=0;k<3;k++)
								d_bs_d[k] = (molecule->atoms[d_i].coordinates[k] - molecule->atoms[h_i].coordinates[k])/rdh * xd 
									+ (molecule->atoms[d_i].coordinates[k] - molecule->atoms[h_i].coordinates[k])/rdh * xd2;

								for(k=0;k<3;k++)
									d_bs_a[k] = (molecule->atoms[a_i].coordinates[k] - molecule->atoms[h_i].coordinates[k])/rah * xd2;

								for(k=0;k<3;k++) d_bs_h[k] = -d_bs_d[k] + -d_bs_a[k];
							}
						} 
						else
						 {
							// No switching, no gradient
							e_bond_switch = 1.0;
							if (addGradient)
							{
								for(k=0;k<3;k++) d_bs_d[k] = 0;
								for(k=0;k<3;k++) d_bs_a[k] = 0;
								for(k=0;k<3;k++) d_bs_h[k] = 0;
							}
						}
						// Water scaling
						e_scale_w = 1.0;
						if ( molecule->atoms[d_i].prop.atomicNumber == OXYGEN &&  molecule->atoms[a_i].prop.atomicNumber == OXYGEN)
						{
							// Count hydrogens and other atoms in vicinity
							double hydrogens = 0.0;
							double others = 0.0;
							for (k = 0; k < molecule->nAtoms; k++) 
							{
								if ( molecule->atoms[k].prop.atomicNumber  == HYDROGEN)
								{
									hydrogens += cvalence_contribution(& molecule->atoms[d_i],&molecule->atoms[k]);
								} else 
								{
									others += cvalence_contribution(&molecule->atoms[d_i],&molecule->atoms[k]);
								}
							}

							// If it is water
							if (hydrogens >= 1.0 )
							{
								sign_wat = 1.0;
								slope = parameters->multiplier_WH_O - 1.0;
								v = hydrogens;
								fv = 0.0;
								if (v > 1.0 && v <= 2.0)
								{
									fv = v - 1.0;
									sign_wat = 1.0;
								}
								if (v > 2.0 && v < 3.0)
								{
									fv = 3.0 - v;
									sign_wat = -1.0;
								}
								fv2 = 1.0 - others;
								if (fv2 < 0.0) fv2 = 0.0;
								e_scale_w = 1.0 + slope * fv * fv2;
							}
						}

						// Charged groups
						e_scale_chd = 1.0;
						e_scale_cha = 1.0;

						// Scaled groups: NR4+
						if (1 && molecule->atoms[d_i].prop.atomicNumber == NITROGEN)
						{
							slope = parameters->multiplier_NH4 - 1.0;
							v = 0.0;
							for (k = 0; k <  molecule->nAtoms; k++) v += cvalence_contribution(&molecule->atoms[d_i],&molecule->atoms[k]);
							if (v > 3.0) v = v - 3.0; else v = 0.0;
							e_scale_chd = 1.0 + slope * v;
						}

						// Scaled groups: COO-
						f_O1 = 0.0;
						f_O2 = 0.0;
						f_cc = 0.0;

						o1 = a_i;
						o2 = -1;
						cc = -1;
						if (molecule->atoms[a_i].prop.atomicNumber == OXYGEN)
						{
							slope = parameters->multiplier_COO - 1.0;

							// Search for closest C atom
							double cdist = 9.9e9;
							cv_O1 = 0.0;
							for (k = 0; k <  molecule->nAtoms; k++)
							{
								v = cvalence_contribution(& molecule->atoms[o1],&molecule->atoms[k]);
								cv_O1 += v; // Sum O1 valence
								if (v > 0.0 && molecule->atoms[k].prop.atomicNumber == CARBON && getDistance(&molecule->atoms[o1], &molecule->atoms[k]) < cdist) {
									cdist =  getDistance(&molecule->atoms[o1], &molecule->atoms[k]);
									cc = k;
								}
							}

							// If C found, look for the second O
							if (cc != -1) {
								double odist = 9.9e9;
								cv_cc = 0.0;
								for (k = 0; k <  molecule->nAtoms; k++)
								{
									v = cvalence_contribution(&molecule->atoms[cc],&molecule->atoms[k]);
									cv_cc += v;
									if (v > 0.0 && k != o1 &&molecule->atoms[k].prop.atomicNumber == OXYGEN && getDistance(&molecule->atoms[cc], &molecule->atoms[k]) < odist) {
										odist = getDistance(&molecule->atoms[cc], &molecule->atoms[k]);
										o2 = k;
									}
								}
							}


							// O1-C-O2 triad:
							if (o2 != -1) {
								// Get O2 valence
								cv_O2 = 0.0;
								for (k = 0; k <  molecule->nAtoms; k++)
										cv_O2 += cvalence_contribution(&molecule->atoms[o2], &molecule->atoms[k]);

								f_O1 = 1.0 - fabs(1.0 - cv_O1);
								if (f_O1 < 0.0) f_O1 = 0.0;

								f_O2 = 1.0 - fabs(1.0 - cv_O2);
								if (f_O2 < 0.0) f_O2 = 0.0;

								f_cc = 1.0 - fabs(3.0 - cv_cc);
								if (f_cc < 0.0) f_cc = 0.0;

								e_scale_cha = 1.0 + slope * f_O1 * f_O2 * f_cc;
							}

						}

						// Final energy
						e_corr = e_para * e_radial * e_angular * e_bond_switch * e_scale_w * e_scale_chd * e_scale_cha;
						e_corr_sum += e_corr;
				
						// Total gradient
						// radial
						factor = e_para * e_angular * e_bond_switch * e_scale_w * e_scale_chd * e_scale_cha;
						addVectorToGrad(molecule,d_i, d_radial_d, factor);
						addVectorToGrad(molecule,a_i, d_radial_a, factor);

						// angular
						factor =  e_para * e_radial * e_bond_switch * e_scale_w * e_scale_chd * e_scale_cha;
						addVectorToGrad(molecule,d_i, d_angular_d, factor);
						addVectorToGrad(molecule,a_i, d_angular_a, factor);
						addVectorToGrad(molecule,h_i, d_angular_h, factor);
						// bond_switch
						factor =   e_para * e_radial * e_angular * e_scale_w * e_scale_chd * e_scale_cha;
						addVectorToGrad(molecule,d_i, d_bs_d, factor);
						addVectorToGrad(molecule,a_i, d_bs_a, factor);
						addVectorToGrad(molecule,h_i, d_bs_h, factor);
						// water scaling
						if (addGradient && e_scale_w != 1.0) 
						{
							slope = parameters->multiplier_WH_O - 1.0;
							for (k = 0; k <  molecule->nAtoms; k++)
							{ 
								if (k != d_i) 
								{
									x = getDistance(&molecule->atoms[d_i], &molecule->atoms[k]);
									if (molecule->atoms[k].prop.atomicNumber == HYDROGEN) 
									{
									xd = cvalence_contribution_d(&molecule->atoms[d_i], &molecule->atoms[k]) * sign_wat;
									for(c=0;c<3;c++) g[c] =(molecule->atoms[d_i].coordinates[c]-molecule->atoms[k].coordinates[c]) * -xd/x * slope;

									addVectorToGrad(molecule, d_i, g, -e_para * e_radial * e_angular * e_bond_switch * e_scale_chd * e_scale_cha);
									addVectorToGrad(molecule, k, g, e_para * e_radial * e_angular * e_bond_switch * e_scale_chd * e_scale_cha);
									}
									else 
									{
									xd = cvalence_contribution_d(&molecule->atoms[d_i], &molecule->atoms[k]);
									for(c=0;c<3;c++) g[c] =(molecule->atoms[d_i].coordinates[c]-molecule->atoms[k].coordinates[c]) * xd/x * slope;
									addVectorToGrad(molecule, d_i, g, -e_para * e_radial * e_angular * e_bond_switch * e_scale_chd * e_scale_cha);
									addVectorToGrad(molecule, k, g, e_para * e_radial * e_angular * e_bond_switch * e_scale_chd * e_scale_cha);
									}
								}
							}
						}
						// scaled groups: NR4+
						if (addGradient && e_scale_chd != 1.0) 
						{
							slope = parameters->multiplier_NH4 - 1.0;
							for (k = 0; k <  molecule->nAtoms; k++)
							{ 
								if (k != d_i)
								{
									x = getDistance(&molecule->atoms[d_i], &molecule->atoms[k]);
									xd = cvalence_contribution_d(&molecule->atoms[d_i], &molecule->atoms[k]);
									for(c=0;c<3;c++) g[c] =(molecule->atoms[d_i].coordinates[c]-molecule->atoms[k].coordinates[c]) * -xd/x * slope;
									addVectorToGrad(molecule, d_i, g, -e_para * e_radial * e_angular * e_bond_switch * e_scale_cha * e_scale_w);
									addVectorToGrad(molecule, k, g, e_para * e_radial * e_angular * e_bond_switch * e_scale_cha * e_scale_w);
								}
							}
						}
						// scaled groups: COO-
						if (addGradient && f_O1 * f_O2 * f_cc != 0.0)
						{
							slope = parameters->multiplier_COO - 1.0;
							// Atoms around O1
							for (k = 0; k <  molecule->nAtoms; k++)
							{
								if (k != o1) 
								{
									xd = cvalence_contribution_d(&molecule->atoms[o1], &molecule->atoms[k]);
									if (xd != 0.0)
									{
										x = getDistance(&molecule->atoms[o1], &molecule->atoms[k]);
										if (cv_O1 > 1.0) xd *= -1.0; 
										xd *= f_O2 * f_cc;
										for(c=0;c<3;c++) g[c] =(molecule->atoms[o1].coordinates[c]-molecule->atoms[k].coordinates[c]) * -xd/x * slope;
										 addVectorToGrad(molecule, o1, g, -e_para * e_radial * e_angular * e_bond_switch * e_scale_chd * e_scale_w);
										 addVectorToGrad(molecule, k, g, e_para * e_radial * e_angular * e_bond_switch * e_scale_chd * e_scale_w);
									}
								}
							}
							slope = parameters->multiplier_COO - 1.0;
							// Atoms around O2
							for (k = 0; k <  molecule->nAtoms; k++)
							{ 
								if (k != o2) 
								{
									xd = cvalence_contribution_d(&molecule->atoms[o2], &molecule->atoms[k]);
									if (xd != 0.0) 
									{
										x =  getDistance(&molecule->atoms[o2], &molecule->atoms[k]);
										if (cv_O2 > 1.0) xd *= -1.0; 
										xd *= f_O1 * f_cc;
										for(c=0;c<3;c++) g[c] =(molecule->atoms[o2].coordinates[c]-molecule->atoms[k].coordinates[c]) * -xd/x * slope;
										addVectorToGrad(molecule, o2, g, -e_para * e_radial * e_angular * e_bond_switch * e_scale_chd * e_scale_w);
										addVectorToGrad(molecule, k, g, e_para * e_radial * e_angular * e_bond_switch * e_scale_chd * e_scale_w);
									}
								}
							}
							slope = parameters->multiplier_COO - 1.0;
							for (k = 0; k <  molecule->nAtoms; k++)
							{ 
								if (k != cc) 
								{
									xd = cvalence_contribution_d(&molecule->atoms[cc], &molecule->atoms[k]);
									if (xd != 0.0) 
									{
										x =  getDistance(&molecule->atoms[cc], &molecule->atoms[k]);
										if (cv_cc > 3.0) xd *= -1.0; 
										xd *= f_O1 * f_O2;
										for(c=0;c<3;c++) g[c] =(molecule->atoms[cc].coordinates[c]-molecule->atoms[k].coordinates[c]) * -xd/x * slope;
										addVectorToGrad(molecule, cc, g, -e_para * e_radial * e_angular * e_bond_switch * e_scale_chd * e_scale_w);
										addVectorToGrad(molecule, k, g, e_para * e_radial * e_angular * e_bond_switch * e_scale_chd * e_scale_w);
									}
								}
							}
						}
					}
				}
			}
		}
		
	}
	}
	}}

	return e_corr_sum;
}
/**********************************************************************/
static boolean getOneDouble(char* buffer, char* name, double* pval)
{
	char* st = strstr(buffer,name);
	if(st) 	
	{
		char* beg = strstr(buffer,"=");
		int k = 0;
		if(beg) k=sscanf(beg+1,"%lf",pval);
		else k=sscanf(st+strlen(st),"%lf",pval);
		if(k==1) return TRUE;
	}
	return FALSE;
}
/**********************************************************************/
static boolean getOneString(char* buffer, char* name, char*pval)
{
	char* st = strstr(buffer,name);
	if(st) 	
	{
		char* beg = strstr(buffer,"=");
		int k = 0;
		if(beg) k=sscanf(beg+1,"%s",pval);
		else k=sscanf(st+strlen(st),"%s",pval);
		if(k==1) return TRUE;
	}
	return FALSE;
}
/**********************************************************************/
static void printParameters(HyhrogenBondCorrectionParameters* parameters, FILE* file)
{
	fprintf(file,"---------------------------------------------------\n");
	fprintf(file,"Hyhrogen-Bond Correction Parameters:\n");
	fprintf(file,"Method=%s\n",parameters->method);
	fprintf(file,"OH_O=%f\n",parameters->para_OH_O);
	fprintf(file,"OH_N=%f\n",parameters->para_OH_N);
	fprintf(file,"NH_O=%f\n",parameters->para_NH_O);
	fprintf(file,"NH_N=%f\n",parameters->para_NH_N);
	fprintf(file,"WH_O=%f\n",parameters->multiplier_WH_O);
	fprintf(file,"COO=%f\n",parameters->multiplier_COO);
	fprintf(file,"NH4=%f\n",parameters->multiplier_COO);
	fprintf(file,"HB_R_CUTOFF=%f\n", parameters->HB_R_CUTOFF);
	fprintf(file,"HB_R_0=%f\n", parameters->HB_R_0);
	fprintf(file,"MAX_XH_BOND=%f\n", parameters->MAX_XH_BOND);
	fprintf(file,"HH_REPULSION_K=%f\n", parameters->HH_REPULSION_k);
	fprintf(file,"HH_REPULSION_E=%f\n", parameters->HH_REPULSION_e);
	fprintf(file,"HH_REPULSION_R0=%f\n", parameters->HH_REPULSION_r0);
	fprintf(file,"---------------------------------------------------\n");
}
/**********************************************************************/
int readHydrogenBondCorrectionParameters(HyhrogenBondCorrectionParameters* parameters, char* fileName)
{
	char buffer[BSIZE];
	FILE* file = fopen(fileName,"r");
	sprintf(parameters->method,"%s","Generic");
/* default cutoffs */
	parameters->HB_R_CUTOFF = 5.5;
	parameters->HB_R_0 = 1.5;
	parameters->MAX_XH_BOND=1.15;
/* default HH repulsion*/
	parameters->HH_REPULSION_k = 0.4;
	parameters->HH_REPULSION_e = 12.7;
	parameters->HH_REPULSION_r0 = 2.3;

/* init H4 repulsion tp -1, to check after read*/
	parameters->para_OH_O = -1;
	parameters->para_OH_N = -1;
	parameters->para_NH_O = -1;
	parameters->para_NH_N = -1;
	parameters->multiplier_WH_O = -1;
	parameters->multiplier_COO = -1;
	parameters->multiplier_NH4 = -1;
	while(file && !feof(file))
	{
		if(!fgets(buffer,BSIZE,file))break;
		uppercase(buffer);
		if(!getOneDouble(buffer, "OH_O", &parameters->para_OH_O))
		if(!getOneDouble(buffer, "OH_N", &parameters->para_OH_N))
		if(!getOneDouble(buffer, "NH_O", &parameters->para_NH_O))
		if(!getOneDouble(buffer, "NH_N", &parameters->para_NH_N))
		if(!getOneDouble(buffer, "WH_O", &parameters->multiplier_WH_O))
		if(!getOneDouble(buffer, "COO", &parameters->multiplier_COO))
		if(!getOneDouble(buffer, "NH4", &parameters->multiplier_NH4))
		if(!getOneDouble(buffer, "HB_R_CUTOFF", &parameters->HB_R_CUTOFF))
		if(!getOneDouble(buffer, "HB_R_0", &parameters->HB_R_0))
		if(!getOneDouble(buffer, "MAX_XH_BOND", &parameters->MAX_XH_BOND))
		if(!getOneDouble(buffer, "HH_REPULSION_K", &parameters->HH_REPULSION_k))
		if(!getOneDouble(buffer, "HH_REPULSION_E", &parameters->HH_REPULSION_e))
		if(!getOneDouble(buffer, "HH_REPULSION_R0", &parameters->HH_REPULSION_r0))
		getOneString(buffer, "METHOD", parameters->method);
	}
	if(!file)
	{
		fprintf(stderr,"Sorry I caanot opent %s file\n",fileName);
		exit(1);
	}
	if(
		parameters->para_OH_O<0 || 
		parameters->para_OH_N<0 || 
		parameters->para_NH_O<0 || 
		parameters->para_NH_N<0 || 
		parameters->multiplier_WH_O<0  ||
		parameters->multiplier_COO<0 ||
		parameters->multiplier_NH4<0
	)
	{
		fprintf(stderr,"Parameters =\n");
		printParameters(parameters, stderr);

		fprintf(stderr,"I cannot read all required parameters\n"
	 	"Required parameters : OH_O, OH_N, NH_O, NH_N, WH_O, COO, NH4\n"
		"Format of file \n"
		"OH_O= value1\n"
		"OH_N= value2\n"
		"NH_O= value3\n"
		"NH_N= value4\n"
		"WH_O= value5\n"
		"COO= value5\n"
		"NH4= value7\n"
		);
		fprintf(stderr,"Here are the default parameters for several method\n");
		printDefParameters("SCC-DFTB", stderr);
		printDefParameters("OM3", stderr);
		printDefParameters("PM6", stderr);
		printDefParameters("AM1", stderr);
		printDefParameters("RM1", stderr);
		printDefParameters("PM3", stderr);
		exit(1);
		return 1;
	}
	printParameters(parameters, stdout);
	return 0;
}
/**********************************************************************/
int setHydrogenBondCorrectionParameters(HyhrogenBondCorrectionParameters* parameters, char* fileName, char* method)
{
	if(fileName) return readHydrogenBondCorrectionParameters(parameters, fileName);
	return setDefParameters(parameters, method);
	return 1;
}
/****************************************************************************************************/
double getH4Correction(Molecule* molecule, HyhrogenBondCorrectionParameters* parameters, boolean addGradient)
{
	double e;
	if(!parameters) return 0; 
	e = getHHRep(parameters, molecule, addGradient);
	/* printf("HH=%f\n",e);*/
	e +=getH4(parameters, molecule, addGradient);
	/* printf("HH+H4=%f\n",e);*/
	return e;
}
