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

const char* mmCLSource = 

__kernel void initEnergy(__global float* energy, int nTerms)
{
	size_t i = get_global_id(0);
	if(i>=nTerms) return;
	energy[i] = 0;
}
__kernel void addEnergyBondAmber(__global float* energy, __global float8* atoms, __global int2* bondIndex, __global float2* bondTerms, int nAtoms, int nTerms)
{
	size_t i = get_global_id(0);
	int ai,aj;
	float term;
	float4 rij;

	if(i>=nTerms) return;

	ai = bondIndex[i].s0;
	aj = bondIndex[i].s1;

	rij.s0 = atoms[ai].s0  - atoms[aj].s0;
	rij.s1 = atoms[ai].s1  - atoms[aj].s1;
	rij.s2 = atoms[ai].s2  - atoms[aj].s2;

	term = sqrt( rij.s0*rij.s0 + rij.s1*rij.s1 + rij.s2*rij.s2) - bondTerms[i].s1;
	energy[i] += bondTerms[i].s0 * term * term;
}
__kernel void addEnergyBendAmber(__global float* energy, __global float8* atoms, __global int4* bondIndex, __global float2* bondTerms, int nAtoms, int nTerms)
{
	size_t i = get_global_id(0);
	int ai,aj,ak;
	float x12, x32, y12, y32, z12, z32, l12, l32, dp;
	float D2RxD2R = 1/(57.29577951308232090712 * 57.29577951308232090712);
	float term;
	float thetaDeg;

	if(i>=nTerms) return;

	ai = bondIndex[i].s0;
	aj = bondIndex[i].s1;
	ak = bondIndex[i].s2;

        x12 = atoms[ai].s0 - atoms[aj].s0;
        y12 = atoms[ai].s1 - atoms[aj].s1;
        z12 = atoms[ai].s2 - atoms[aj].s2;

        x32 = atoms[ak].s0 - atoms[aj].s0;
        y32 = atoms[ak].s1 - atoms[aj].s1;
        z32 = atoms[ak].s2 - atoms[aj].s2;

       	l12 = sqrt( x12 * x12 + y12 * y12 + z12 * z12 );
       	l32 = sqrt( x32 * x32 + y32 * y32 + z32 * z32 );
        if( l12 != 0.0  &&  l32 != 0.0)
	{
		dp = ( x12 * x32 + y12 * y32 + z12 * z32 ) / (l12 * l32 );
		if ( dp < -1.0 ) dp = -1.0;
		else if ( dp > 1.0 ) dp = 1.0;
    		thetaDeg = 57.29577951308232090712 * acos(dp);
	}

	term = thetaDeg - bondTerms[i].s1;
	term *= term *  bondTerms[i].s0;
	term *= D2RxD2R;

	energy[i] += term;
}
__kernel void addEnergyDihedralAmber(__global float* energy, __global float8* atoms, __global int4* bondIndex, __global float4* bondTerms, int nAtoms, int nTerms)
{
	size_t i = get_global_id(0);
	int ai, aj, ak, al;
	float phiDeg= 0;
	float	D2R = 1.0/57.29577951308232090712;
	float	PI = 3.14159265358979323846;
	float   xij, yij, zij;
       	float   xkj, ykj, zkj;
	float   xkl, ykl, zkl;
      	float   dx, dy, dz;
        float   gx, gy, gz;
        float   bi, bk;
        float   ct, d, ap, app, bibk;

	if(i>=nTerms) return;

	ai = bondIndex[i].s0;
	aj = bondIndex[i].s1;
	ak = bondIndex[i].s2;
	al = bondIndex[i].s3;

        xij = atoms[ai].s0 -  atoms[aj].s0;
        yij = atoms[ai].s1 -  atoms[aj].s1;
        zij = atoms[ai].s2 -  atoms[aj].s2;

        xkj = atoms[ak].s0 -  atoms[aj].s0;
        ykj = atoms[ak].s1 -  atoms[aj].s1;
        zkj = atoms[ak].s2 -  atoms[aj].s2;

        xkl = atoms[ak].s0 -  atoms[al].s0;
        ykl = atoms[ak].s1 -  atoms[al].s1;
        zkl = atoms[ak].s2 -  atoms[al].s2;

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
		
	phiDeg = 0.0;
	if ( bibk >= 1.0e-6 )
	{
		ct = ct / sqrt( bibk );
        
		if( ct < -1.0 ) ct = -1.0;
        	else if( ct > 1.0 ) ct = 1.0;
        	ap = acos( ct );
        
		d  = xkj*(dz*gy-dy*gz) + ykj*(dx*gz-dz*gx) + zkj*(dy*gx-dx*gy);
        
		if( d < 0.0 ) ap = -ap;
        
		ap = PI - ap;
    		app = 180.0 * ap / PI;
       		if( app > 180.0 ) app = app - 360.0;
		phiDeg = app;
	}

	energy[i] += bondTerms[i].s1/bondTerms[i].s0 * ( 1.0 + cos( D2R*(bondTerms[i].s3 * phiDeg - bondTerms[i].s2 )) );
}
__kernel void addEnergyImproperTorsionAmber(__global float* energy, __global float8* atoms, __global int4* bondIndex, __global float4* bondTerms, int nAtoms, int numberOfStretchTerms)
{
	size_t i = get_global_id(0);
}
__kernel void addEnergyfNonBondedAmber(__global float* energy, __global float8* atoms, __global int2* bondIndex, __global float4* bondTerms, int nAtoms, int nTerms, int useCoulomb)
{
	size_t i = get_global_id(0);
	int ai, aj;
	float rij2, rij6, rij12, coulombTerm, factorNonBonded;
	float rijx, rijy, rijz;
	float chargei, chargej, Aij, Bij, rij;
	float permittivityScale = 1, permittivity = 1;
	float coulombFactor;

	if(i>=nTerms) return;

	coulombFactor = 332.05382 / ( permittivity * permittivityScale );

	ai = bondIndex[i].s0;
	aj = bondIndex[i].s1;

	Aij    = bondTerms [i].s0;
	Bij    = bondTerms [i].s1;
	factorNonBonded = bondTerms[i].s2;

	chargei = atoms[ai].s3;
	chargej =  atoms[aj].s3;


	rijx =  atoms[ai].s0-atoms[aj].s0;
	rijy =  atoms[ai].s1-atoms[aj].s1;
	rijz =  atoms[ai].s2-atoms[aj].s2;

	rij2 = rijx * rijx + rijy * rijy + rijz * rijz;
	rij = sqrt( rij2 );
	rij6 = rij2 * rij2 * rij2;
	rij12 = rij6 * rij6;

	if(useCoulomb!=0) coulombTerm = ( chargei * chargej * coulombFactor*factorNonBonded ) / rij;
	else coulombTerm = 0.0;

	energy[i] += Aij / rij12 - Bij / rij6 + coulombTerm;
}
__kernel void addEnergyHydrogenBondedAmber(__global float* energy, __global float8* atoms, __global int2* bondIndex, __global float2* bondTerms, int nAtoms, int nTerms)
{
	size_t i = get_global_id(0);
	int ai, aj;
	float rij2, rij6, rij12;
	float rijx, rijy, rijz;
	float rij4, rij10;
	float Cij, Dij;

	if(i>=nTerms) return;

	ai = bondIndex[i].s0;
	aj = bondIndex[i].s1;

	Cij    = bondTerms [i].s0;
	Dij    = bondTerms [i].s1;

	rijx =  atoms[ai].s0-atoms[aj].s0;
	rijy =  atoms[ai].s1-atoms[aj].s1;
	rijz =  atoms[ai].s2-atoms[aj].s2;

	rij2 = rijx * rijx + rijy * rijy + rijz * rijz;
	if ( rij2 < 1.0e-2 ) rij2 = 1.0e-2;	
	rij4 = rij2 * rij2;
	rij6 = rij4 * rij2;
	rij10 = rij6 * rij4;
	rij12 = rij10 * rij2;

	energy[i] += Cij / rij12 - Dij / rij10;
}
__kernel void addEnergyPairWise(__global float* energy, __global float8* atoms, __global int2* bondIndex, __global float8* bondTerms, int nAtoms, int nTerms, int useCoulomb, int useVanderWals)
{
	size_t i = get_global_id(0);
	int ai, aj;
	float rij2, rij6, rij8, rij10;
	float coulombTerm;
	float rijx, rijy, rijz;
	float chargei, chargej, rij;
	float permittivityScale = 1, permittivity = 1;
	float coulombFactor;
	float A, Beta;
	float  B6, B8, B10;
	float c6, c8, c10, b;

	if(i>=nTerms) return;

	coulombFactor = 332.05382/ ( permittivity * permittivityScale );

	ai = bondIndex[i].s0;
	aj = bondIndex[i].s1;

	A    = bondTerms [i].s0;
	Beta = bondTerms [i].s1;
	c6   = bondTerms [i].s2;
	c8   = bondTerms [i].s3;
	c10  = bondTerms [i].s4;
	b    = bondTerms [i].s5;

	chargei = atoms[ai].s3;
	chargej = atoms[aj].s3;

	rijx =  atoms[ai].s0-atoms[aj].s0;
	rijy =  atoms[ai].s1-atoms[aj].s1;
	rijz =  atoms[ai].s2-atoms[aj].s2;

	rij2 = rijx * rijx + rijy * rijy + rijz * rijz;
	if(rij2<1e-2) rij = 1e-2;

	rij = sqrt( rij2 );
	rij6 = rij2 * rij2 * rij2;
	rij8 = rij6* rij2;
	rij10 = rij8 * rij2;

	if(useCoulomb!=0) coulombTerm = ( chargei * chargej * coulombFactor ) / rij;
	else coulombTerm = 0.0;

	B6  = 0;
	B8  = 0;
	B10 = 0;
	if(useVanderWals!=0)
	{
		float fact = 1.0;
		float s = 1.0;
		float br = b*rij;
		float brk = 1.0;
		int k;

		if(fabs(c6)>1e-12)
		{
			for(k=1;k<=2*3;k++)
			{
				fact *= k;
				brk *= br;
				s += brk/fact;
			}
			B6 = c6*(1-exp(-br)*s);
		}

		if(fabs(c8)>1e-12)
		{
			fact = 1.0;
			s = 1.0;
			br = b*rij;
			for(k=1;k<=2*4;k++)
			{
				fact *= k;
				brk *= br;
				s += brk/fact;
			}
			B8 = c8*(1-exp(-br)*s);
		}

		if(fabs(c10)>1e-12)
		{
			fact = 1.0;
			s = 1.0;
			br = b*rij;
			for(k=1;k<=2*5;k++)
			{
				fact *= k;
				brk *= br;
				s += brk/fact;
			}
			B10 = c10*(1-exp(-br)*s);
		}
	}
				
	energy[i] += A*exp(-Beta*rij)
		- B6 / rij6 
		- B8 / rij8 
		- B10 / rij10 
		+ coulombTerm;
}
__kernel void initVelocities( __global float8* atoms, int nAtoms)
{
	size_t i = get_global_id(0);
	if(i>=nAtoms) return;
	atoms[i].s5 = 0;
	atoms[i].s6 = 0;
	atoms[i].s7 = 0;
}
__kernel void addGradientBondAmber(__global float* energy, __global float8* atoms, __global int2* bondIndex, __global float2* bondTerms, int nAtoms, int nTerms)
{
	size_t i = get_global_id(0);
	int ai, aj;
	float rijx, rijy, rijz, forceConstant, equilibriumDistance, term;
	float forceix, forceiy, forceiz;
	float bondLength;

	if(i>=nTerms) return;

	ai = bondIndex[i].s0;
	aj = bondIndex[i].s1;

	forceConstant = bondTerms [i].s0;
	equilibriumDistance = bondTerms [i].s1;

	rijx =  atoms[ai].s0-atoms[aj].s0;
	rijy =  atoms[ai].s1-atoms[aj].s1;
	rijz =  atoms[ai].s2-atoms[aj].s2;

	bondLength = sqrt( rijx * rijx + rijy * rijy + rijz * rijz );

	if ( bondLength < 1.0e-10 ) bondLength = 1.0e-10;

	term = - 2*forceConstant * ( bondLength - equilibriumDistance ) / bondLength;
	forceix = term * rijx;
	forceiy = term * rijy;
	forceiz = term * rijz;

	atoms[ai].s5 -= forceix;
	atoms[ai].s6 -= forceiy;
	atoms[ai].s7 -= forceiz;
		
	atoms[aj].s5 += forceix;
	atoms[aj].s6 += forceiy;
	atoms[aj].s7 += forceiz;
}
__kernel void addGradientBendAmber(__global float* energy, __global float8* atoms, __global int4* bondIndex, __global float2* bondTerms, int nAtoms, int nTerms)
{
	size_t i = get_global_id(0);
	float	PI = 3.14159265358979323846;
	float	D2R = 1.0/57.29577951308232090712;
	int ai, aj, ak;
	float x12, y12, z12, x32, y32, z32, l12, l32, dp;

	float term;
	float thetaRad, cosTheta;
	float denominator, absTheta;
	float delta = 1e-10;

	float rijx, rijy, rijz;
	float rkjx, rkjy, rkjz;
	float rij2, rij, rkj2, rkj,rij3, rkj3;
	float denominatori, denominatork;

	float forceix, forceiy, forceiz;
	float forcejx, forcejy, forcejz;
	float forcekx, forceky, forcekz;

	float rijDotrkj;
	float term2ix, term2iy, term2iz;
	float term2jx, term2jy, term2jz;
	float term2kx, term2ky, term2kz;

	if(i>=nTerms) return;

	ai = bondIndex[i].s0;
	aj = bondIndex[i].s1;
	ak = bondIndex[i].s2;

        x12 = atoms[ai].s0 - atoms[aj].s0;
        y12 = atoms[ai].s1 - atoms[aj].s1;
        z12 = atoms[ai].s2 - atoms[aj].s2;

        x32 = atoms[ak].s0 - atoms[aj].s0;
        y32 = atoms[ak].s1 - atoms[aj].s1;
        z32 = atoms[ak].s2 - atoms[aj].s2;

       	l12 = sqrt( x12 * x12 + y12 * y12 + z12 * z12 );
       	l32 = sqrt( x32 * x32 + y32 * y32 + z32 * z32 );
        if( l12 != 0.0  &&  l32 != 0.0)
	{
		dp = ( x12 * x32 + y12 * y32 + z12 * z32 ) / (l12 * l32 );
		if ( dp < -1.0 ) dp = -1.0;
		else if ( dp > 1.0 ) dp = 1.0;
    		thetaRad = acos(dp);
	}

	absTheta = fabs( thetaRad );
	cosTheta = cos( thetaRad );

	if ( ( absTheta > delta ) && ( absTheta < PI - delta ) )
	{
		denominator = sin(thetaRad);
			
		if ( denominator < 1.0e-10 ) denominator = 1.0e-10;

		term = 2*bondIndex[i].s0 * (thetaRad/PI*180.0 - bondIndex[i].s1) / denominator;
		term *= D2R;

        	rijx = atoms[ai].s0 - atoms[aj].s0;
        	rijy = atoms[ai].s1 - atoms[aj].s1;
        	rijz = atoms[ai].s2 - atoms[aj].s2;

        	rkjx = atoms[ak].s0 - atoms[aj].s0;
        	rkjy = atoms[ak].s1 - atoms[aj].s1;
        	rkjz = atoms[ak].s2 - atoms[aj].s2;

		rij2 = rijx * rijx + rijy * rijy + rijz * rijz;
		rij = sqrt( rij2 );

		rkj2 = rkjx * rkjx + rkjy * rkjy + rkjz * rkjz;
		rkj = sqrt( rkj2 );

		rijDotrkj = rijx * rkjx + rijy * rkjy + rijz * rkjz;

		rij3 = rij2 * rij;
		rkj3 = rkj2 * rkj;

		denominatori = rij3 * rkj;
		if ( denominatori < 1.0e-10 ) denominatori = 1.0e-10;

		denominatork = rij * rkj3;
		if ( denominatork < 1.0e-10 ) denominatork = 1.0e-10;
			
		term2ix = ( rij2 * rkjx - rijDotrkj * rijx ) / denominatori;
		term2iy = ( rij2 * rkjy - rijDotrkj * rijy ) / denominatori;
		term2iz = ( rij2 * rkjz - rijDotrkj * rijz ) / denominatori;
			
		term2kx = ( rkj2 * rijx - rijDotrkj * rkjx ) / denominatork;
		term2ky = ( rkj2 * rijy - rijDotrkj * rkjy ) / denominatork;
		term2kz = ( rkj2 * rijz - rijDotrkj * rkjz ) / denominatork;
			
		term2jx = - term2ix - term2kx;
		term2jy = - term2iy - term2ky;
		term2jz = - term2iz - term2kz;
			
		forceix = term * term2ix;
		forceiy = term * term2iy;
		forceiz = term * term2iz;
			
		forcejx = term * term2jx;
		forcejy = term * term2jy;
		forcejz = term * term2jz;
			
		forcekx = term * term2kx;
		forceky = term * term2ky;
		forcekz = term * term2kz;

		atoms[ai].s5 -= forceix;
		atoms[ai].s6 -= forceiy;
		atoms[ai].s7 -= forceiz;
		
		atoms[aj].s5 -= forcejx;
		atoms[aj].s6 -= forcejy;
		atoms[aj].s7 -= forcejz;

		atoms[aj].s5 -= forcekx;
		atoms[aj].s6 -= forceky;
		atoms[aj].s7 -= forcekz;
	}
}
__kernel void addGradientDihedralAmber(__global float* energy, __global float8* atoms, __global int4* bondIndex, __global float4* bondTerms, int nAtoms, int nTerms)
{
	size_t i = get_global_id(0);
	int ai,aj,ak,al;
	int j;
	float	D2R = 1.0/57.29577951308232090712;

	float rjix, rjiy, rjiz;
	float rkjx, rkjy, rkjz;
	float rkix, rkiy, rkiz;
	float rljx, rljy, rljz;
	float rlkx, rlky, rlkz;

	float forceix, forceiy, forceiz;
	float forcejx, forcejy, forcejz;
	float forcekx, forceky, forcekz;
	float forcelx, forcely, forcelz;

	float rkj;
	float xt, yt, zt;
	float xu, yu, zu;
	float xtu, ytu, ztu;
	float rt2, ru2, rtru;
	float cosine1, sine1, cosineN, sineN, cosold, sinold;
	float cosPhase, sinPhase;
	float dedxt, dedyt, dedzt;
	float dedxu, dedyu, dedzu;
	float dedphi;
	int n;
	float vn;

	if(i>=nTerms) return;

	ai = bondIndex[i].s0;
	aj = bondIndex[i].s1;
	ak = bondIndex[i].s2;
	al = bondIndex[i].s3;

        rjix = atoms[aj].s0 -  atoms[ai].s0;
        rjiy = atoms[aj].s1 -  atoms[ai].s1;
        rjiz = atoms[aj].s2 -  atoms[ai].s2;

        rkjx = atoms[ak].s0 -  atoms[aj].s0;
        rkjy = atoms[ak].s1 -  atoms[aj].s1;
        rkjz = atoms[ak].s2 -  atoms[aj].s2;

        rlkx = atoms[al].s0 -  atoms[ak].s0;
        rlky = atoms[al].s1 -  atoms[ak].s1;
        rlkz = atoms[al].s2 -  atoms[ak].s2;

	xt = rjiy*rkjz - rkjy*rjiz;
	yt = rjiz*rkjx - rkjz*rjix;
	zt = rjix*rkjy - rkjx*rjiy;

	xu = rkjy*rlkz - rlky*rkjz;
	yu = rkjz*rlkx - rlkz*rkjx;
	zu = rkjx*rlky - rlkx*rkjy;

	xtu = yt*zu - yu*zt;
	ytu = zt*xu - zu*xt;
	ztu = xt*yu - xu*yt;

	rt2 = xt*xt + yt*yt + zt*zt;
	ru2 = xu*xu + yu*yu + zu*zu;

	rtru = sqrt(rt2 * ru2);

	rkj = sqrt(rkjx*rkjx + rkjy*rkjy + rkjz*rkjz);
	cosine1 = 1.0;
	sine1   = 0.0;

	if (rtru <1e-10) rtru = 1e-10;
	if (rt2 <1e-10) rt2 = 1e-10;
	if (ru2 <1e-10) ru2 = 1e-10;

	cosine1 = (xt*xu + yt*yu + zt*zu) / rtru;
	sine1 = (rkjx*xtu + rkjy*ytu + rkjz*ztu) / (rkj*rtru);

	n = (int)bondTerms[i].s3;
	cosPhase = cos(D2R*bondTerms[i].s2);
	sinPhase = cos(D2R*bondTerms[i].s2);
	vn = bondTerms[i].s1/bondTerms[i].s0;

	cosineN = cosine1;
	sineN   = sine1;

	for(j=2;j<=n;j++)
	{
	   cosold = cosineN;
	   sinold = sineN;
	   cosineN = cosine1*cosold - sine1*sinold;
	   sineN   = cosine1*sinold + sine1*cosold;
	}

	dedphi = vn*n*(cosineN*sinPhase-sineN*cosPhase);

        rkix = atoms[ak].s0 -  atoms[ai].s0;
        rkiy = atoms[ak].s1 -  atoms[ai].s1;
        rkiz = atoms[ak].s2 -  atoms[ai].s2;

        rljx = atoms[al].s0 -  atoms[aj].s0;
        rljy = atoms[al].s1 -  atoms[aj].s1;
        rljz = atoms[al].s2 -  atoms[aj].s2;

	dedxt = dedphi * (yt*rkjz - rkjy*zt) / (rt2*rkj);
	dedyt = dedphi * (zt*rkjx - rkjz*xt) / (rt2*rkj);
	dedzt = dedphi * (xt*rkjy - rkjx*yt) / (rt2*rkj);

	dedxu = -dedphi * (yu*rkjz - rkjy*zu) / (ru2*rkj);
	dedyu = -dedphi * (zu*rkjx - rkjz*xu) / (ru2*rkj);
	dedzu = -dedphi * (xu*rkjy - rkjx*yu) / (ru2*rkj);

	forceix = rkjz*dedyt - rkjy*dedzt;
	forceiy = rkjx*dedzt - rkjz*dedxt;
	forceiz = rkjy*dedxt - rkjx*dedyt;

	forcejx = rkiy*dedzt - rkiz*dedyt + rlkz*dedyu - rlky*dedzu;
	forcejy = rkiz*dedxt - rkix*dedzt + rlkx*dedzu - rlkz*dedxu;
	forcejz = rkix*dedyt - rkiy*dedxt + rlky*dedxu - rlkx*dedyu;

	forcekx = rjiz*dedyt - rjiy*dedzt + rljy*dedzu - rljz*dedyu;
	forceky = rjix*dedzt - rjiz*dedxt + rljz*dedxu - rljx*dedzu;
	forcekz = rjiy*dedxt - rjix*dedyt + rljx*dedyu - rljy*dedxu;

	forcelx = rkjz*dedyu - rkjy*dedzu;
	forcely = rkjx*dedzu - rkjz*dedxu;
	forcelz = rkjy*dedxu - rkjx*dedyu;

	atoms[ai].s5 += forceix;
	atoms[ai].s6 += forceiy;
	atoms[ai].s7 += forceiz;

	atoms[aj].s5 += forcejx;
	atoms[aj].s6 += forcejy;
	atoms[aj].s7 += forcejz;

	atoms[ak].s5 += forcekx;
	atoms[ak].s6 += forceky;
	atoms[ak].s7 += forcekz;

	atoms[al].s5 += forcelx;
	atoms[al].s6 += forcely;
	atoms[al].s7 += forcelz;
}
__kernel void addGradientImproperTorsion(__global float* energy, __global float8* atoms, __global int4* bondIndex, __global float4* bondTerms, int nAtoms, int nTerms)
{
	size_t i = get_global_id(0);
}
__kernel void addGradientNonBondedAmber(__global float* energy, __global float8* atoms, __global int2* bondIndex, __global float4* bondTerms, int nAtoms, int nTerms, int useCoulomb)
{
	size_t i = get_global_id(0);
	int ai, aj;

	float rijx, rijy, rijz;

	float forceix, forceiy, forceiz;
	float forcejx, forcejy, forcejz;

	float permittivityScale = 1, permittivity = 1;
	float coulombFactor, factorNonBonded;
	float rij2, rij;
	float rij3;
	float chargei, chargej, coulombTerm;
	float Aij, Bij, rij6, rij7, rij14, rij8;
	float  term1, term2, term3;

	if(i>=nTerms) return;

	coulombFactor = 332.05382 / ( permittivity * permittivityScale );

	ai = bondIndex[i].s0;
	aj = bondIndex[i].s1;

	Aij    = bondTerms [i].s0;
	Bij    = bondTerms [i].s1;
	factorNonBonded = bondTerms[i].s2;

	chargei = atoms[ai].s3;
	chargej =  atoms[aj].s3;

	rijx =  atoms[ai].s0-atoms[aj].s0;
	rijy =  atoms[ai].s1-atoms[aj].s1;
	rijz =  atoms[ai].s2-atoms[aj].s2;

	rij2 = rijx * rijx + rijy * rijy + rijz * rijz;
	if ( rij2 < 1.0e-2 ) rij2 = 1.0e-2;	

	rij = sqrt( rij2 );
	rij3 = rij2 * rij;
	rij6 = rij3 * rij3;
	rij7 = rij6 * rij;
	rij8 = rij7 * rij;
	rij14 = rij7 * rij7;
	if(useCoulomb!=0) coulombTerm = ( chargei * chargej * coulombFactor*factorNonBonded ) / rij3;
	else coulombTerm = 0.0;

	term1 = 12 * Aij / rij14;
	term2 = 6 * Bij / rij8;
	term3 = term1 - term2 + coulombTerm;
	forceix = term3 * rijx;
	forceiy = term3 * rijy;
	forceiz = term3 * rijz;
	forcejx = - forceix;
	forcejy = - forceiy;
	forcejz = - forceiz;

	atoms[ai].s5 -= forceix;
	atoms[ai].s6 -= forceiy;
	atoms[ai].s7 -= forceiz;

	atoms[aj].s5 -= forcejx;
	atoms[aj].s6 -= forcejy;
	atoms[aj].s7 -= forcejz;
}
__kernel void addGradientHydrogenBondedAmber(__global float* energy, __global float8* atoms, __global int2* bondIndex, __global float2* bondTerms, int nAtoms, int nTerms)
{
	size_t i = get_global_id(0);

	int ai, aj;

	float rijx, rijy, rijz;

	float forceix, forceiy, forceiz;
	float forcejx, forcejy, forcejz;

	float Cij, Dij, rij2,  rij4, rij8, rij12, rij14;
	float  term1, term2, term3;

	if(i>=nTerms) return;

	ai = bondIndex[i].s0;
	aj = bondIndex[i].s1;

	Cij    = bondTerms [i].s0;
	Dij    = bondTerms [i].s1;

	rijx =  atoms[ai].s0-atoms[aj].s0;
	rijy =  atoms[ai].s1-atoms[aj].s1;
	rijz =  atoms[ai].s2-atoms[aj].s2;

	rij2 = rijx * rijx + rijy * rijy + rijz * rijz;

	if ( rij2 < 1.0e-2 ) rij2 = 1.0e-2;	

	rij4 = rij2 * rij2;
	rij8 = rij4 * rij4;
	rij12 = rij8 * rij4;
	rij14 = rij12 * rij2;
	term1 = Cij / rij14;
	term2 = Dij / rij12;
	term3 = term1 - term2;
	forceix = term3 * rijx;
	forceiy = term3 * rijy;
	forceiz = term3 * rijz;
	forcejx = - forceix;
	forcejy = - forceiy;
	forcejz = - forceiz;

	atoms[ai].s5 -= forceix;
	atoms[ai].s6 -= forceiy;
	atoms[ai].s7 -= forceiz;

	atoms[aj].s5 -= forcejx;
	atoms[aj].s6 -= forcejy;
	atoms[aj].s7 -= forcejz;
}
__kernel void addGradientPairWise(__global float* energy, __global float8* atoms, __global int2* bondIndex, __global float8* bondTerms, int nAtoms, int nTerms, int useCoulomb, int useVanderWals)
{
	size_t i = get_global_id(0);
	int ai, aj;

	float rijx, rijy, rijz;

	float forceix, forceiy, forceiz;
	float forcejx, forcejy, forcejz;

	float permittivityScale = 1, permittivity = 1;
	float coulombFactor;
	float rij2, rij;
	float rij3;
	float chargei, chargej, coulombTerm;
	float rij6, rij7, rij8, rij9, rij10, rij11, rij12;
	float  term1, term6, term8, term10, termAll;
	float A, Beta, C6, C8, C10,b;
	float s, sp, fact, br, brk, ebr;
	int n, k;

	if(i>=nTerms) return;

	coulombFactor = 332.05382 / ( permittivity * permittivityScale );
	ai = bondIndex[i].s0;
	aj = bondIndex[i].s1;

	A    = bondTerms [i].s0;
	Beta = bondTerms [i].s1;
	C6   = bondTerms [i].s2;
	C8   = bondTerms [i].s3;
	C10  = bondTerms [i].s4;
	b    = bondTerms [i].s5;

	chargei = atoms[ai].s3;
	chargej = atoms[aj].s3;

	rijx =  atoms[ai].s0-atoms[aj].s0;
	rijy =  atoms[ai].s1-atoms[aj].s1;
	rijz =  atoms[ai].s2-atoms[aj].s2;

	rij2 = rijx * rijx + rijy * rijy + rijz * rijz;
	if ( rij2 < 1.0e-2 ) rij2 = 1.0e-2;	

	rij = sqrt( rij2 );
	rij3 = rij2 * rij;
	rij6 = rij3 * rij3;
	rij7 = rij6 * rij;
	rij8 = rij7 * rij;
	rij9 = rij8 * rij;
	rij10 = rij9 * rij;
	rij11 = rij10 * rij;
	rij12 = rij11 * rij;
	if(useCoulomb!=0) coulombTerm = ( chargei * chargej * coulombFactor ) / rij3;
	else coulombTerm = 0.0;
		
		
	term1 = -A*Beta/rij*exp(-Beta*rij);

	br = b*rij;
	ebr = exp(-b*rij);

	term6 =   0.0;
	if(useVanderWals!=0 && fabs(C6)>1e-12)
	{
		fact = 1.0;
		s = 1.0;
		n = 3;
		brk = 1.0;
		for(k=1;k<2*n;k++)
		{
			fact *= k;
			brk *= br;
			s += brk/fact;
		}
		sp = s*b;
		fact *=2*n;
		brk *= br;
		s += brk/fact;
		term6 =   b*C6*ebr*s/rij7
			-(2*n)*C6*(1-ebr*s)/rij8
			-C6*ebr/rij7*sp;
	}
	term8 =   0.0;
	if(useVanderWals!=1 && fabs(C8)>1e-12)
	{
		fact = 1.0;
		s = 1.0;
		n = 4;
		brk = 1.0;
		for(k=1;k<2*n;k++)
		{
			fact *= k;
			brk *= br;
			s += brk/fact;
		}
		sp = s*b;
		fact *=2*n;
		brk *= br;
		s += brk/fact;
		term8 =   b*C8*ebr*s/rij9
			-(2*n)*C8*(1-ebr*s)/rij10
			-C8*ebr/rij9*sp;
	}

	term10 =   0.0;
	if(useVanderWals!=0 && fabs(C10)>1e-12)
	{
		fact = 1.0;
		s = 1.0;
		n = 5;
		brk = 1.0;
		for(k=1;k<2*n;k++)
		{
			fact *= k;
			brk *= br;
			s += brk/fact;
		}
		sp = s*b;

		fact *=2*n;
		brk *= br;
		s += brk/fact;
		term10 =   b*C10*ebr*s/rij11
			-(2*n)*C10*(1-ebr*s)/rij12
			-C10*ebr/rij11*sp;
	}

	termAll = term1 - term6 - term8 - term10 + coulombTerm;


	forceix = termAll * rijx;
	forceiy = termAll * rijy;
	forceiz = termAll * rijz;
	forcejx = - forceix;
	forcejy = - forceiy;
	forcejz = - forceiz;

	atoms[ai].s5 -= forceix;
	atoms[ai].s6 -= forceiy;
	atoms[ai].s7 -= forceiz;

	atoms[aj].s5 -= forcejx;
	atoms[aj].s6 -= forcejy;
	atoms[aj].s7 -= forcejz;
}
