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
// SRmap.h: interface for the SRmap class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(SRMAP_INCLUDED)
#define SRMAP_INCLUDED

class SRelement;
class SRnode;
class SRface;
class SRedge;
class SRcoord;

class SRmap
{
public:
	SRmap();
	void Setup();
	int EdgeShapeFunctions(double re, double N[]);
	void FaceShapeFunctions(SRface* face, double rf, double sf, double N[]);
	void FaceShapeFunctionsLinear(SRface* face, double rf, double sf, double N[]);
	void QuadShapeFunctions(double rf, double sf, double N[]);
	double BrickTetMap(double rb, double sb, double tb, double& r, double& s, double& t);
	double QuadTriMap(double rq, double sq, double& r, double& s);
	double TetMidEdgeMappingDeriv(int n0, int n1, double dldra[], double L[]);
	void AssignBrickMidNodes(int midnode, int node0, int node1);
	double EdgeTangent(SRedge* edge, double r, SRvec3& dxdr, bool normalize = false);
	void ElementCentroid(SRelement* elem, SRvec3& pos);
	void ElementCentroid(SRelement* elem, double& r, double& s, double& t);
	void BrickWedgeQuadInit();
	void FaceCentroid(SRface* face, double& rc, double& sc);
	void FaceCentroid(SRface* face, SRvec3& pos);
	void ElementNaturalCoordsAtMidedge(SRelement* elem, int localEdge, double& r, double& s, double& t);
	void ElementNaturalCoordsFromEdge(SRelement* elem, int localEdge, double re, double& r, double& s, double& t);
	void ElementNaturalCoordsFromEdge(SRelement* elem, SRvec3& edgePos, double& r, double& s, double& t);
	void ElementNaturalCoordsFromFace(SRelement* elem, SRvec3& facePos, double& r, double& s, double& t);
	void ElementNaturalCoordsFromFace(SRelement* elem, int lface, double rf, double sf, double& r, double& s, double& t);
	double EdgeArcLength(SRedge* edge, double r);
	void TetVolumeCoords(double r, double s, double t, double lv[]);
	void TriangleAreaCoords(double r, double s, double lv[]);
	void TetVolumeCoords(double r, double s, double t, double& l1, double& l2, double& l3, double& l4);
	void TriangleAreaCoords(double r, double s, double& l1, double& l2, double& l3);
	SRvec3& GetBrickNodeCoord(int i){ return brickNodes[i]; };
	double GetTriNoder(int i){ return trinoder[i]; };
	double GetTriNodes(int i){ return trinodes[i]; };
	double GetQuadNoder(int i){ return quadnoder[i]; };
	double GetQuadNodes(int i){ return quadnodes[i]; };
	int GetQuadEdgeLocalNode(int ledge, int lnode){	return quadEdgeLocalNodes[ledge][lnode];};
	int GetTriEdgeLocalNode(int ledge, int lnode){ return triEdgeLocalNodes[ledge][lnode]; };
	int GetTetFaceLocalNode(int lface, int lnode){ return tetFaceLocalNodes[lface][lnode]; };
	int GetTetEdgeLocalNode(int ledge, int lnode){ return tetEdgeLocalNodes[ledge][lnode]; };
	int GetWedgeEdgeLocalNode(int ledge, int lnode){ return wedgeEdgeLocalNodes[ledge][lnode]; };
	int GetBrickEdgeLocalNode(int ledge, int lnode){ return brickEdgeLocalNodes[ledge][lnode]; };
	double GetDLtdr(int i){ return dLtdr[i]; };
	double GetDLtds(int i){ return dLtds[i]; };
	double GetDLtdt(int i){ return dLtdt[i]; };
	double* GetDLtdrVec(){ return dLtdr; };
	double* GetDLtdsVec(){ return dLtds; };
	double* GetDLtdtVec(){ return dLtdt; };
	double ApproxFaceArea(SRface *face);
	double FaceArea(SRface *face);
	bool ElementInverseMap(SRelement *elem, SRvec3 &x, double &r, double &s, double &t, double tolin);
	int GetFaceEdgeLocalNode(int Nedge, int lej, int node);
	void QuadLinearShapeFunctions(double rf, double sf, double N[]);
	void BrickLinearShapeFunctions(double r, double s, double t, double N[]);
	void WedgeLinearShapeFunctions(double r, double s, double t, double N[]);

	void EdgeShapeDerives(double r, double* dNdr);
	void ElementFaceNaturalNormal(SRelement* elem, int lface, SRvec3& norm);
	void FaceEdgeNaturalNormal(SRface* face, int lej, SRvec3& norm);

	SRvec3& getBrickNode(int i){ return brickNodes[i]; };
	SRvec3& getWedgeCorners(int i){ return wedgeCorners[i]; };
	SRvec3& getTetCorners(int i){ return tetCorners[i]; };
	int getQuadFaceEdgeNodes(int i, int j){ return quadFaceEdgeNodes[i][j]; };
	int getTriFaceEdgeNodes(int i, int j){ return triFaceEdgeNodes[i][j]; };
	int getQuadEdgeLocalNodes(int i, int j){ return quadEdgeLocalNodes[i][j]; };
	SRvec3& getBrickNatNormal(int i){ return brickNatNormals[i]; };
	SRvec3& getWedgeNatNormals(int i){ return wedgeNatNormals[i]; };
	SRvec3& getTetNatNormals(int i){ return tetNatNormals[i]; };
	SRvec3& getQuadNatNormals(int i){ return quadNatNormals[i]; };
	SRvec3& getTriNatNormals(int i){ return triNatNormals[i]; };
	bool wasSetupCalled(){ return setupWasCalled; };
	void setSetupWasCalled(bool tf)
	{
		setupWasCalled = tf;
	}

private:
	double dLtdr[4], dLtds[4], dLtdt[4];
	SRvec3 brickNodes[20];
	SRvec3 wedgeCorners[6];
	SRvec3 tetCorners[4];
	double trinoder[6];
	double trinodes[6];
	double quadnoder[8];
	double quadnodes[8];
	int quadEdgeLocalNodes[4][2];
	int triEdgeLocalNodes[3][2];
	int tetFaceLocalNodes[4][3];
	int tetEdgeLocalNodes[6][2];
	int wedgeEdgeLocalNodes[9][2];
	int brickEdgeLocalNodes[12][2];
	int quadFaceEdgeNodes[4][3];
	int triFaceEdgeNodes[3][3];
	SRvec3 brickNatNormals[6];
	SRvec3 wedgeNatNormals[5];
	SRvec3 tetNatNormals[6];
	SRvec3 quadNatNormals[4];
	SRvec3 triNatNormals[3];
	bool setupWasCalled;
};

#endif // !defined(SRMAP_INCLUDED)
