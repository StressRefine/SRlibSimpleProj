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
// SRnode.h: interface for the SRnode class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(SRNODE_INCLUDED)
#define SRNODE_INCLUDED

class SRconstraint;

class SRnode  
{
public:
	SRnode();
	void Create(int useridt, double xt, double yt, double zt){ userId = useridt; pos.Assign(xt, yt, zt); forceId = -1; };
	int GetId(){ return id; };
	int GetUserid(){ return userId; };
	SRvec3& Position(){ return pos; };
	void SetPosition(SRvec3 newPos){ pos.Copy(newPos); };
	double GetXyz(int i){ return pos.d[i]; };
	int GetGlobalFunctionNumber(){ return globalFunctionNumber; };
	void PutGlobalFunctionNumber(int g){globalFunctionNumber = g; };
	bool isMidSide() { return (midSideEdgeOwner != -1); };
	int GetMidSideEdgeOwner() { return midSideEdgeOwner; };
	bool isOrphan() { return (firstElementOwner == -1); };
	void SetAsMidside(int edgeId) { midSideEdgeOwner = edgeId; };
	void SetFirstElementOwner(int e) { firstElementOwner = e; };
	int GetFirstElementOwner() { return firstElementOwner; };
	SRconstraint* GetConstraint();
	void SetConstraintId(int i){ constraintId = i; };
	void SetUserId(int u){ userId = u; };
	bool checkSacr(){ return sacrificial; };
	int getConstraintId(){ return constraintId; };
	void setId(bool i){ id = i; };
	SRvec3& getPos(){ return pos; };
	void setPos(SRvec3 &that){ pos.Copy(that); };
	bool GetDisp(SRvec3& disp);

private:
	int id;
	int userId; //user original nodes numbers in case non-contiguous
	SRvec3 pos;
	int globalFunctionNumber;
	int midSideEdgeOwner;
	int firstElementOwner;
	int constraintId;
	bool sacrificial;
	int forceId;
	int dispCoordid;
};

#endif //!defined(SRNODE_INCLUDED)
