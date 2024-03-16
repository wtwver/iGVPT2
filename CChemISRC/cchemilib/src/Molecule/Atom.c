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

/* Atom.c */
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "../Utils/Constants.h"
#include "../Utils/AtomsProp.h"
#include "Atom.h"

/*****************************************************************/
Atom getCopyAtom(Atom *atom)
{
	Atom a;
	a=*atom;
	a.prop = propAtomGet(atom->prop.symbol);
	a.mmType = strdup(atom->mmType);
	a.pdbType = strdup(atom->pdbType);
	a.residueName = strdup(atom->residueName);
	a.typeConnections = NULL;
	return a;
}
/*****************************************************************/
double getDistance(Atom *a1,Atom* a2)
{
	double* C1 = a1->coordinates;
	double* C2 = a2->coordinates;

	double x, y, z;
	
        x = C1[ 0 ] - C2[ 0 ];
       	y = C1[ 1 ] - C2[ 1 ];
       	z = C1[ 2 ] - C2[ 2 ];

	return sqrt( x * x + y * y + z * z );
}
/*****************************************************************/
double getAngle(Atom *a1,Atom* a2,Atom* a3)
{
	double* C1 = a1->coordinates;
	double* C2 = a2->coordinates;
	double* C3 = a3->coordinates;

	double x12, x32, y12, y32, z12, z32, l12, l32, dp;
	
        x12 = C1[ 0 ] - C2[ 0 ];
       	y12 = C1[ 1 ] - C2[ 1 ];
       	z12 = C1[ 2 ] - C2[ 2 ];
       	x32 = C3[ 0 ] - C2[ 0 ];
       	y32 = C3[ 1 ] - C2[ 1 ];
       	z32 = C3[ 2 ] - C2[ 2 ];

       	l12 = sqrt( x12 * x12 + y12 * y12 + z12 * z12 );
       	l32 = sqrt( x32 * x32 + y32 * y32 + z32 * z32 );
        if( l12 == 0.0 )
	{
               	return 0.0;
        }
        if( l32 == 0.0 )
	{
               	return 0.0;
       	}
        dp = ( x12 * x32 + y12 * y32 + z12 * z32 ) / (l12 * l32 );
	if ( dp < -1.0 )
		dp = -1.0;
	else if ( dp > 1.0 )
		dp = 1.0;
    return RADTODEG * acos(dp);
}
/*****************************************************************/
double getTorsion(Atom* a1,Atom* a2,Atom* a3,Atom* a4)
{
	double* C1 = a1->coordinates;
	double* C2 = a2->coordinates;
	double* C3 = a3->coordinates;
	double* C4 = a4->coordinates;

	double   xij, yij, zij;
       	double   xkj, ykj, zkj;
	double   xkl, ykl, zkl;
      	double   dx, dy, dz;
        double   gx, gy, gz;
        double   bi, bk;
        double   ct, d, ap, app, bibk;

        xij = C1[ 0 ] - C2[ 0 ];
        yij = C1[ 1 ] - C2[ 1 ];
        zij = C1[ 2 ] - C2[ 2 ];
        xkj = C3[ 0 ] - C2[ 0 ];
        ykj = C3[ 1 ] - C2[ 1 ];
        zkj = C3[ 2 ] - C2[ 2 ];
        xkl = C3[ 0 ] - C4[ 0 ];
        ykl = C3[ 1 ] - C4[ 1 ];
        zkl = C3[ 2 ] - C4[ 2 ];

        dx = yij * zkj - zij * ykj;
        dy = zij * xkj - xij * zkj;
        dz = xij * ykj - yij * xkj;
        gx = zkj * ykl - ykj * zkl;
        gy = xkj * zkl - zkj * xkl;
        gz = ykj * xkl - xkj * ykl;

        bi = dx * dx + dy * dy + dz * dz;
        bk = gx * gx + gy * gy + gz * gz;
        ct = dx * gx + dy * gy + dz * gz;
		bibk = bi * bk;
		
		if ( bibk < 1.0e-6 )	
			return 0;
        
		ct = ct / sqrt( bibk );
        
		if( ct < -1.0 )
                ct = -1.0;
        else if( ct > 1.0 )
                ct = 1.0;

        ap = acos( ct );
        
		d  = xkj*(dz*gy-dy*gz) + ykj*(dx*gz-dz*gx) + zkj*(dy*gx-dx*gy);
        
		if( d < 0.0 )
                ap = -ap;
        
		ap = PI - ap;
       	app = 180.0 * ap / PI;
       	if( app > 180.0 )
               	app = app - 360.0;
        return( app );
}
