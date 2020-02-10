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
// SRconstraint.h: interface for the SRconstraint class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(SRCONSTRAINT_INCLUDED)
#define SRCONSTRAINT_INCLUDED

enum SRconstraintType{ nodalCon, faceCon, inactiveCon };

class SRconstraint
{
public:
	SRconstraint();
	void ProcessFaceConstraint();
	void FillFaceEnforcedDispCoeffs(int dof, double *dispvec);
	void Copy(SRconstraint& that);
	void AddNodalConstraint(SRconstraint& that, bool SummingEnfd = false);
	void operator = (SRconstraint& c2){ Copy(c2); };
	void allocateEnforcedDisplacementData(int n){ enforcedDisplacementData.Allocate(n, 3); };
	void PutEnforcedDisplacementData(int n, int dof, double val);
	double GetEdgeEnforcedDisp(double re, int dof);
	double GetFaceEnforcedDisp(double rf, double sf, int dof);
	double GetEnforcedDisp(int nodeNum, int dof);
	int GetNumEnforcedDisp(){ return enforcedDisplacementData.getNumCols(); };
	bool allDofsConstrained();
	void Clear();
	void SetConstrainedDof(int i, bool tf = true){ constrainedDof[i] = tf; };
	void SetType(SRconstraintType typeIn){ type = typeIn; };
	void SetEntityId(int i){ entityId = i; };
	bool IsConstrainedDof(int i){ return constrainedDof[i]; };
	int GetEntityId(){ return entityId; };
	SRconstraintType GetType(){ return type; };
	SRcoord* GetCoord();
	bool hasEnforcedDisp(){ return !enforcedDisplacementData.isEmpty(); };
	int GetCoordId(){ return coordId; };
	void SetCoordId(int id){ coordId = id; };
	bool isGcs(){ return (coordId == -1); };
	void getDisp(int n, SRvec3& enfd);
	double getDisp(int n, int d);
	int getId(){ return id; };
	void setId(int idIn){ id = idIn; };

private:
	bool constrainedDof[3];
	int id;
	int entityId; //node, edge, or face
	SRconstraintType type;
	int coordId;
	SRdoubleMatrix enforcedDisplacementData;
};


#endif // !defined(SRCONSTRAINT_INCLUDED)
