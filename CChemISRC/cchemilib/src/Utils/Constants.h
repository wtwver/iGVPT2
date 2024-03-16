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


#ifndef __CCHEMILIB_CONSTANTS_H__
#define __CCHEMILIB_CONSTANTS_H__

/*
static gdouble hbar = 6.62606891e-34/2/PI;// PRL 1998 NIST
static gdouble e = 1.602176462e-19;
static gdouble a0 = 0.5291772083e-10;
static gdouble me = 9.10938188e-31;
static gdouble c =299792458.0;
static gdouble OneUnderfourPIEpsilon0 = 1e-7*299792458.0*299792458.0;
static gdouble NA =6.0221415e23; //(6.0221415(10)x10^23; from NIST
static gdouble kb =1.3806505e-23; //1.380 6505(24) x 10^-23 J K-1  from NIST
*/

#define BSIZE 1024
#define BOHRTOANG   0.5291772083
#define ANGTOBOHR  (1.0/BOHRTOANG)
#define RADTODEG   57.29577951308232090712
#define AUTODEB  2.54158059
#define AUTOEV  27.21138469
#define PI   3.14159265358979323846
#define AMUTOAU 1822.88848121
#define AUTOCM1 219474.63633664 
#define AUTOKCAL 627.509544796
#define KCALTOKJ 4.18400000
/* double AUTOfs = a0^2 me hbar^-1*1e15 */
#define AUTOfs 0.0241888427177214
#define fsInAU (1/AUTOfs)
#define KbInAU 3.1668114e-6
/* double fsInAKMA = 1.0/sqrt(1e-10*1e-10*1.6605655e-27*6.022045e23/4184.0)/1e15;*/
#define fsInAKMA 0.020454828110640
#define kcalcmM1 349.75511054
/* Kb = 6.022045e23*1.380662e-23/4184.0 */
/* Kb in AKMA */
#define Kb  1.9871914e-3
/* hbar in AKMA */
/* 6.62606957e-34*1.43932636e+20*0.020454828110640*1e15/(2pi) */
#define hbar   (6.62606957e-34*1.43932636e20*0.020454828110640*1e15/2.0/PI)



#define DEGTORAD   0.017453293 
/* c en m/s */
#define slight (299792458.0) 
#define NAvogadro (6.02214129e23)
#define hPlank (6.6260693e-34)
#define epsilon0 (8.854187817e-12)


#endif /* __CCHEMILIB_CONSTANTS_H__ */
