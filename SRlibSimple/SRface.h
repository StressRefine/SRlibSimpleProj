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
// SRface.h: interface for the SRface class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(SRFACE_INCLUDED)
#define SRFACE_INCLUDED

#include "SRutil.h"

class SRlocalEdge;
class SRconstraint;
class SRelement;
class SRforce;

class SRfaceUtil
{
public:
	//static functions used during creation of global faces
	static int GetLocalEdge(int n1, int n2);
	static int GlobalFaceMatch(int& gn1, int& gn2, int& gn3, int& gn4, int n1, int n2, int n3, int n4);
	static bool QuadFaceMatch(int n1, int n2, int n3, int n4, int g1, int g2, int g3, int g4, int& gn1, int& gn2, int& gn3, int& gn4);
	static bool FaceMatch(int nv[], int gv[], int gnv[]);
	static bool FaceMatch(int n1, int n2, int n3, int n4, int g1, int g2, int g3, int g4, int& gn1, int& gn2, int& gn3, int& gn4);
	static void GetLocalEdgeNodes(int nej, int lej, int& n1, int& n2);
	static void QuadLocalEdgeNodes(int lej, int& n1, int& n2);
};

class SRface  
{
	friend class SRfaceUtil;

public:
	SRface();
	SRconstraint* GetConstraint();
	void GetLocalEdgeNodes(int lej, int &n1, int &n2);
	void QuadFaceNaturalCoordinatesFromEdge(int lej, double rej, double& rf, double& sf);
	void naturalCoordinatesFromEdge(int lej, double rej, double& rf, double& sf, bool useDirection = false);
	void Clear();
	void ProcessForce(SRforce *force, SRvec3& ResF);
	bool GetForceValue(SRforce *force, double rf, double sf, double forceVal[3]);
	void GetSummedForceValue(double rf, double sf, double forceVal[3]);
	SRnode* GetNode(int localnodenum);

	SRedge* GetEdge(int localedgenum);
	void Create(int idt, int n1, int n2, int n3, int n4 = -1);
	int GetNumNodes(){ return (nodeIds[3] == -1 ? 3 : 4); };
	int GetNodeId(int i);
	int GetNumLocalEdges(){ return localEdges.GetNum(); };
	int GetLocalEdgePOrder(int lej){ return localEdges.Get(lej).GetPOrder(); }
	int GetLocalEdgeGlobalId(int lej){ return localEdges.Get(lej).GetGlobalEdgeId(); };
	int GetLocalEdgeDirection(int lej){	return localEdges.Get(lej).GetDirection(); };
	int GetNumGlobalFunctions(){ return globalFunctionNumbers.GetNum();	};
	int GetGlobalFunctionNumber(int i){ return globalFunctionNumbers.Get(i); };
	int GetConstraintId(){ return constraintId; };
	int GetElementOwner(int i){ return elementOwners[i]; };
	void setElementOwner(int i, int e){ elementOwners[i] = e; };
	int GetElementLocalFace(int i){ return elementLocalFace[i]; };
	int GetMidNodeId(int lej){ return localEdges.GetPointer(lej)->GetMidNodeId(); };
	bool hasForce(){ return !forceIds.isEmpty(); };
	bool hasConstraint(){ return constraintId != -1; };
	double GetTractionJump(){ return tractionJump; };
	void SetTractionJump(double jump){ tractionJump = jump; };
	void NodeNaturalCoordinates(int lnode, double &rf, double &sf);
	void NaturalCoordinatesNearNode(int lnode, double &rf, double &sf);
	void LocalNodePosition(int lnode, SRvec3& pos);
	bool GetFlipNormal(){ return flipNormal; };
	void SetFlipNormal();
	void FillMappingNodes();
	double UnitTriad(double rf, double sf, SRvec3 &p, SRvec3 &e1, SRvec3 &e2, SRvec3 &e3, bool detjonly = false);
	void QuadShapeDerivs(double rf, double sf);
	void PutElementOwner(int n, int el){ elementOwners[n] = el; };
	void PutElementLocalFace(int n, int lf){ elementLocalFace[n] = lf; }
	void AssignLocalEdge(int lej, int gej, int direction){ localEdges.GetPointer(lej)->Assign(gej, direction); };
	bool IsBoundaryFace(){ return (elementOwners[1] == -1); };
	void AllocateGlobalFunctionNumbers(int n){ globalFunctionNumbers.Allocate(n); };
	void PutGlobalFunctionNumber(int lfun, int gfun){ globalFunctionNumbers.Put(lfun, gfun); }
	bool isFlat();
	void NaturalCoordinatesNearMidedge(int lej, double& r, double& s);
	void Position(double rf, double sf, SRvec3& p);
	void OutwardNormal(double rf, double sf, SRvec3& norm, bool checkOutwardLocally = true);
	double Jacobian(double rf, double sf);
	double GetSize(){ return size; };
	int GetId(){ return id; };
	bool natCoordsAtCorner(int cornerId, double& r, double& s);
	int GetNumNodesTotal(){ return 2 * localEdges.GetNum(); };
	int midNodeMatch(int mid);
	void RotateVec(int coordid, double rf, double sf, double forceVal[]);
	int GetNodeOrMidNodeId(int i);
	SRnode* GetNodeOrMidnode(int localnodenum);
	void OutwardNormalLinear(double rf, double sf, SRvec3& norm);
	double UnitTriadLinear(double rf, double sf, SRvec3 &p, SRvec3 &e1, SRvec3 &e2, SRvec3 &e3);
	void QuadLinearShapeDerivs(double rf, double sf);

	void pushBackForceids(int id){ forceIds.PushBack(id); };
	void setConstraintId(int id){ constraintId = id; };
	int* getNodeIds(){ return nodeIds; };
	SRlocalEdge* getLocalEdge(int e){ return localEdges.GetPointer(e); };
	int getMultifaceForceGroupId(){ return multifaceForceGroupId; };
	void setMultifaceForceGroupId(int id){ multifaceForceGroupId = id; };
	int getForceId(int i){ return forceIds.Get(i); };
	int getNumForceIds(){ return forceIds.GetNum(); };
	void freeForceIds(){ forceIds.Free(); };

private:
	int nodeIds[4];
	SRvector <SRlocalEdge> localEdges;
	SRintVector globalFunctionNumbers;
	int constraintId;
	int elementOwners[2];
	int elementLocalFace[2];
	SRintVector forceIds;
	double tractionJump;
	bool flipNormal;
	double xnode[8], ynode[8], znode[8];
	double dNdrf[8], dNdsf[8];
	double size;
	int flat;
	int id;
	int multifaceForceGroupId;
#ifdef _DEBUG
	int nodeUids[8]; //(for debugging)
#endif
};

#endif // !defined(SRFACE_INCLUDED)
