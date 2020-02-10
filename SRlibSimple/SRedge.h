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
// SRedge.h: interface for the SRedge class.
//
//
//////////////////////////////////////////////////////////////////////

	

#if !defined(SREDGE_INCLUDED)
#define SREDGE_INCLUDED

class SRedge;
class SRface;
class SRforce;

class SRedgeUtil
{
public:
	//static functions used during creation of global edges
	static bool EdgeMatch(int n1, int n2, int gn1, int gn2, int& direction);
	static int GlobalEdgeMatchN2(int n1, int n2, int& direction);
	static int GlobalEdgeMatch(int n1, int n2, int& direction);
};


class SRlocalEdge
{
public:
	SRlocalEdge::SRlocalEdge();
	int GetPOrder();
	SRedge* GetEdge();
	int GetGlobalEdgeId(){ return globalEdgeId; };
	int GetMidNodeId();
	int GetDirection(){ return direction; };
	void Assign(int id, int directiont){ globalEdgeId = id; direction = directiont; };
private:
	int globalEdgeId;
	int direction;
};

class SRboundaryFaceData
{
public:
	int faceId;
	int localEdgeId;
};

class SRedge
{
	friend class SRedgeUtil;
public:
	SRedge();
	SRconstraint* GetConstraint();
	void GetDisplacement(double r, SRvec3 &disp);
	void Clear();
	void GetForceValue(double r, SRforce* force, double forceVal[]);
	void Create(int n1, int n2, int& midnode, int idt);
	SRnode* GetNode(int localnodenum);
	int GetNodeId(int i);
	int GetMidNodeId(){ return midnodeId; };
	int GetNodeOrMidNodeId(int n)
	{
		if (n < 2)
			return GetNodeId(n);
		else
			return midnodeId;
	};
	int GetMidNodeUserId();
	int GetPorder(){ return pOrder; };
	int GetPrevPorder(){ return prevpOrder; };
	int GetNumGlobalFunctions(){ return globalFunctionNumbers.GetNum(); };
	int GetGlobalFunctionNumber(int i){ return globalFunctionNumbers.Get(i); };
	void FillMappingNodes();
	void AllocateGlobalFunctionNumbers(int n){ globalFunctionNumbers.Allocate(n); };
	void putPorder(int p);
	void AssignGlobalFunctionNumbers(int& fun, int pmin, int pmax);
	bool isStraight(){ return straight; };
	void Position(double r, SRvec3& p);
	bool PChanged(){ return pOrder != prevpOrder; };
	int GetConstraintId(){ return constraintId; };
	void SetConstraintId(int c){ constraintId = c; };
	double GetSize(){ return size; };
	void SetThin(bool thinin){ thin = thinin; };
	void Straighten(double fraction = 1.0);
	void ScaleStraightenVsEdgeLength();
	void basisFunctions(double r, double* basis);
	bool isSacrificial(){ return sacrificial; };
	void setSacrificial(bool tf = true){ sacrificial = tf; };
	void getDisp(double r, SRvec3& disp);
	int GetId(){ return id; };
	double getXnode(int i){ return xnode[i]; };
	double getYnode(int i){ return ynode[i]; };
	double getZnode(int i){ return znode[i]; };
	bool isOrphan();
	void ProcessForce(SRforce* force, SRvec3& ResF);
	void PutBoundaryFaceId(int f, int lej);
	bool checkKinkOK(SRface* face0, int lej0, SRface* face1, int lej1);
	double getStraightenFraction(){ return straightenFraction; };
	int getNumBoundaryfaceData(){ return boundaryfaceData.GetNum(); };
	SRboundaryFaceData& getBoundaryFaceData(int i){ return boundaryfaceData.Get(i); };

private:
	int id;
	int nodeIds[2];
	int pOrder;
	int prevpOrder;
	SRintVector globalFunctionNumbers;
	int midnodeId;
	int constraintId;
	SRvector <SRboundaryFaceData> boundaryfaceData;
	double xnode[3], ynode[3], znode[3];
	double size;
	bool straight;
	bool straightened;
	bool thin;
	bool sacrificial;
	int forceId;
	double initialBow;
	double straightenFraction;
};

#endif //!defined(SREDGE_INCLUDED)
