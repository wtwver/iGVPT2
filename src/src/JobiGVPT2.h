/*****************************************************************************************
 iGVPT2 is a program for computing anharmonic corrections to vibration frequencies, 
 based on force field expansion of the potential energy surface in normal mode coordinates.
 iGVPT2 supports several computation chemistry packages(CCP).

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
*****************************************************************************************/

#ifndef __IGVPT2_JOBIGVPT2_H__
#define __IGVPT2_JOBIGVPT2_H__

struct JobiGVPT2
{
	char* inputFileName;
};
JobiGVPT2 newJobiGVPT2(char* inputFileName);
void freeJobiGVPT2(JobiGVPT2* job);
void runJobiGVPT2(JobiGVPT2* job);
#endif /*  __IGVPT2_JOBIGVPT2_H__ */

