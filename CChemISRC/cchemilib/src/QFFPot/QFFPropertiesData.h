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

/* QFFPropertiesData.h */

#ifndef __QFF_PropertiesData_H__
#define __QFF_PropertiesData_H__

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

typedef struct _QFFPropertiesData        QFFPropertiesData;
typedef struct _QFFPropertiesDataClass   QFFPropertiesDataClass;

struct _QFFPropertiesData
{
	int nFrequencies;
	int nDim; /* 3 = vector, 6 = tensor */
	double* zero;
	double **first;
	double ***second;
	double ****cubic;
        QFFPropertiesDataClass* klass;

};

struct _QFFPropertiesDataClass
{
        void (*readData)(QFFPropertiesData* qffPropertiesData, char* fileName);
        void (*convertToAU)(QFFPropertiesData* qffPropertiesData);
        void (*convertToAU2)(QFFPropertiesData* qffPropertiesData);
        void (*free)(QFFPropertiesData* qffPropertiesData);
        void (*print)(QFFPropertiesData* qffPropertiesData);
};

QFFPropertiesData newQFFPropertiesData(int n, int nDim);

#endif /* __QFF_PropertiesData_H__ */
