/*
Copyright (c) 2020 Richard King

The StressRefine library is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

The StressRefine library is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

The terms of the GNU General Public License are explained in the file COPYING,
also available at <https://www.gnu.org/licenses/>

*/

//////////////////////////////////////////////////////////////////////
//
// SRmachDep.h: interface for the SRmachDep class.
// POSSIBLY MACHINE DEPENDENT FUNCTIONS
//
//////////////////////////////////////////////////////////////////////

#if !defined(SRMACHDEP_INCLUDED)
#define SRMACHDEP_INCLUDED


#define SRBUFSIZE 8192

enum dlgfilemode{newmode,openmode,saveasmode};

#ifdef _WINDOWS
#endif


class SRstring;
class SRmachDep  
{
public:
	static void GetTime(SRstring &t, bool timeOnly = false);
	static void ProgReturn();
	static bool CreateDir(char* name);
	static double availMemCheck();
};

#endif // !defined(SRMACHDEP_INCLUDED)
