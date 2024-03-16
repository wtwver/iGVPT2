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


#ifndef __CCHEMILIB_TYPES_H__
#define __CCHEMILIB_TYPES_H__

#define MAJOR_VERSION    0
#define MINOR_VERSION    0
#define MICRO_VERSION    2

#define TRUE 1
#define FALSE 0
typedef int boolean;
#ifdef DISABLE_CONST
#define CONST
#else
#define CONST const
#endif

#ifdef OS_WIN32

#define DIR_SEPARATOR '\\'
#define DIR_SEPARATOR_S "\\"

#else  /* !OS_WIN32 */
/* Unix */

#define DIR_SEPARATOR '/'
#define DIR_SEPARATOR_S "/"

#endif /* !OS_WIN32 */

typedef enum
{
 LOW_LAYER=0, MEDIUM_LAYER, HIGH_LAYER
}CCHEMILayerType;

typedef enum
{
	FILETYPE_UNKNOWN = 0,
	FILETYPE_CCHEMI,
	FILETYPE_GAUSSIANOUT,
	FILETYPE_GAUSSIANFCHK,
	FILETYPE_MOLDEN,
	FILETYPE_GABEDIT,
	FILETYPE_GAMESSOUT,
	FILETYPE_XYZ,
	FILETYPE_MOLPROOUT,
	FILETYPE_QCHEMOUT,
	FILETYPE_NWCHEMOUT,
	FILETYPE_ORCAOUT,
	FILETYPE_ORCAHESSIAN,
	FILETYPE_VASPOUT,
	FILETYPE_MOPACOUT,
	FILETYPE_MOPACAUX,
	FILETYPE_MOPACSCAN,
	FILETYPE_MOPACIRC,
	FILETYPE_TURBOMOLEOUT,
}CCHEMIFileType;

#define MAXISOTOP 10
typedef struct _SAtomsProp
{
	char *name;
	char *symbol;
	int atomicNumber;
	double covalentRadii;
	double bondOrderRadii;
	double vanDerWaalsRadii;
	double radii;
	int maximumBondValence;
	double mass;
	double electronegativity;
	int color[3];
	int nIsotopes;
	int iMass[MAXISOTOP];
	double rMass[MAXISOTOP];
	double abundances[MAXISOTOP];
}SAtomsProp;

typedef enum
{
  VERLET = 0,
  BEEMAN = 1,
  STOCHASTIC = 2,
  LANGEVIN = 3,
  QTB = 4,
  MARTYNATUCKERMAN = 5
} MDIntegratorType;
typedef enum
{
  NONE = 0,
  ANDERSEN = 1,
  BERENDSEN = 2,
  BUSSI = 3,
  NOSEHOOVER = 4
} MDThermostatType;
typedef enum
{
  PIMDTHERMOSTATNONE = 0,
  PIMDTHERMOSTATLANGEVIN = 1,
  PIMDTHERMOSTATNOSEHOOVER = 2,
  PIMDTHERMOSTATQTB = 3
} PIMDThermostatType;
typedef enum
{
  PIMDTRANSFORMATIONNONE = 0,
  PIMDTRANSFORMATIONSTAGING = 1,
  PIMDTRANSFORMATIONNORMALMODE = 2
} PIMDTransformationType;

typedef enum
{
  NOCONSTRAINTS = 0,
  BONDSCONSTRAINTS = 1,
  BONDSANGLESCONSTRAINTS = 2
} Constraints;
typedef enum
{
  CCHEMI_BONDTYPE_SINGLE = 0,
  CCHEMI_BONDTYPE_DOUBLE,
  CCHEMI_BONDTYPE_TRIPLE,
  CCHEMI_BONDTYPE_HYDROGEN,
} CChemIBondType;
typedef struct _BondType BondType;
struct _BondType
{
        int n1;
        int n2;
        CChemIBondType bondType;
};

#ifdef ENABLE_CL
#include<CL/cl.h>
typedef struct _CLProp
{
       cl_context context;
       cl_command_queue command_queue;
	cl_device_id device_id;
}CLProp;
#endif

#endif /* __CCHEMILIB_TYPES_H__ */

