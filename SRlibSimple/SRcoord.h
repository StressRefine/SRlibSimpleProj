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
// SRcoord.h: interface for the SRcoord class: Coordinate system definitions
//

#if !defined(SRCOORD_INCLUDED)
#define SRCOORD_INCLUDED

enum SRcoordType {cartesian,spherical,cylindrical};

#include "SRmath.h"
#include "SRstring.h"

class SRcoord
{
public:
	SRcoord();
	void Create(double x0, double y0, double z0, SRvec3 p13, SRvec3 p3);
	void Create(double x0, double y0, double z0);
	void CalculateBasisVectors(SRvec3& p, SRvec3 &e1l, SRvec3 &e2l, SRvec3 &e3l);
	void CalculateLocalDirection(SRvec3& p, int dof, SRvec3 &el);
	void GetRotationMatrix(bool toLocal, SRvec3& p, SRmat33& R);
	void Copy(SRcoord& c2);
	SRcoordType GetType(){ return type; };
	void GetPos(double &x, double &y, double &z, SRvec3& pos);
	void operator =(SRcoord& c2);
	void setName(const char *nameIn);
	const char* GetName();
	void VecTransform(SRvec3 p, SRvec3 &v);
	int checkParallelToGcs(int dof);
	void setType(SRcoordType typeIn){type = typeIn;};

	//data:
private:
	SRstring name;
	SRcoordType type;
	SRvec3 origin;
	SRvec3 e1, e2, e3;
	bool gcsAligned;
	int uid;
	int otherCoordid;
	bool gcsaligned;
};

#endif // !defined(SRCOORD_INCLUDED)
