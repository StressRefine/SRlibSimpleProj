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
// SRforce.h: interface for the SRforce class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(SRFORCE_INCLUDED)
#define SRFORCE_INCLUDED

class SRforce;
class SRvolumeForce;

enum SRforceType { nodalForce, edgeForce, faceForce, inactiveForce };
enum SRvolumeForceType { gravity, centrifugal };

class SRforce
{
public:
	SRforce();
	SRforceType GetType(){ return type; };
	void setType(SRforceType t)
	{
		type = t;
	};
	void setElemId(int id)
	{
		elemId = id;
	};
	bool isPressure(){ return pressure; };
	void setPressure(bool tf)
	{
		pressure = tf;
	};
	void AllocateForceVals(int n, int m){ forceVals.Allocate(n, m); };
	void zeroForceVals(){ forceVals.Zero(); };
	bool isGcs(){ return coordId == -1; };
	int GetCoordId(){ return coordId; };
	int GetEntityId(){ return entityId; };
	void SetEntityId(int i){ entityId = i; };
	double GetForceVal(int i, int j){ return forceVals.Get(i, j); };
	int GetNumForces(){ return forceVals.getNumCols(); };
	void Copy(SRforce& that, bool copyForceVals = true);
	void Clear();
	void AddNodalForce(SRforce& that, bool SummingForceVals = false);
	void setForceVal(int i, int j, double val){ forceVals.Put(i, j, val); };
	void catForceVal(int i, int j, double val){ forceVals.PlusAssign(i, j, val); };
	void setCoordId(int id){ coordId = id; };
	void setNv(int i, int n){ nv[i] = n; };
	int getNv(int i){ return nv[i]; };
	int getElemid(){ return elemId; };
	int* getNvVec(){ return nv; };
	void setFaceFromNodal(bool tf){ faceFromNodal = tf; };
	bool getfaceFromNodal(){ return faceFromNodal; };
	void allocateStiffMat(int n){ StiffMat.Allocate(n); };
	void allocateStiffDiag(int n){ stiffDiag.Allocate(n); };
	void catStiffMat(int i, double v){ StiffMat.PlusAssign(i, v); };
	void putStiffDiag(int i, int v){ stiffDiag.Put(i, v); };
	void freeStiffMat(){ StiffMat.Free(); };
	void freeStiffDiag(){ stiffDiag.Free(); };
	bool isFaceFromNodal(){ return faceFromNodal; };
	double* getStiffMatd(){ return StiffMat.d; };
	int getStiffDiag(int i){ return stiffDiag.Get(i); };

private:
	SRforceType type;
	int coordId;
	int entityId; //node, edge, or face.
	int elemId; //element id for pressure on faces
	int nv[4]; //nodes at corners of face for pressure on faces
	SRdoubleMatrix forceVals;
	bool pressure;
	bool faceFromNodal;
	SRdoubleVector StiffMat;
	SRintVector stiffDiag;
};


class SRthermalForce
{
public:
	double GetTemp(SRelement* elem, double r, double s, double t);
	bool CeTMult(SRmaterial* mat, double eTx, double eTy, double eTz, double ceT[]);
	void Process();
	bool isConstantTemp(){ return constantTemp; };
	void setConstantTemp(bool tf){constantTemp = tf; };
	double getTemp(){ return temp; };
	void setTemp(double t){ temp = t; };
	double getNodalTemp(int i){ return nodalTemp.Get(i); };
	void setNodalTemp(int i, double t){ nodalTemp.Put(i, t); };
	void allocateNodalTemp(int n){ nodalTemp.Allocate(n); };
	int getNumNodalTemp(){ return nodalTemp.GetNum(); };
	void putNodalTemp(int i, double t){ nodalTemp.Put(i, t); };

private:
	bool constantTemp;
	double temp;
	SRdoubleVector nodalTemp;
};

class SRvolumeForce
{
public:
	SRvolumeForce(){ g1 = g2 = g3 = 0.0; omega2 = 0.0; }
	void GetForceValue(SRelement* elem, SRvec3& p, double val[]);
	SRvolumeForceType GetType(){ return type; };
	void setType(SRvolumeForceType typeIn){ type = typeIn; };
	void setG(SRvec3& g)
	{
		g1 = g.d[0];
		g2 = g.d[1];
		g3 = g.d[2];
	};
	double getG1(){ return g1; };
	double getG2(){ return g2; };
	double getG3(){ return g3; };
	double getAlpha(){ return alpha; };
	void setOmega2(double o2){ omega2 = o2; };
	double getOmega2(){ return omega2; };
	void setAlpha(double a){ alpha = a; };
	void setAxis(SRvec3& a){ axis = a; };
	void setOrigin(SRvec3& o){ origin = o; };
	SRvec3& getAxis(){ return axis; };
	SRvec3& getOrigin(){ return origin; };
	void allocateElList(int n){ elList.Allocate(n); };
	void putElList(int i, int e){ elList.Put(i, e); };
	int getElList(int i){ return elList.Get(i); };
	int getElListLen(){ return elList.GetNum(); };

private:
	SRintVector elList;
	SRvolumeForceType type;
	//for gravity loads:
	double g1, g2, g3;
	//for centrifugal loads:
	double omega2;
	double alpha;
	SRvec3 axis;
	SRvec3 origin;
};

#endif // !defined(SRFORCE_INCLUDED)
