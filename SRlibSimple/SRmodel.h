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
// SRmodel.h: interface for the SRmodel class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(SRMODEL_INCLUDED)
#define SRMODEL_INCLUDED

#include "SRpoundDefines.h"
#include "SRmachDep.h"
#include "SRstring.h"
#include "SRutil.h"
#include "SRmath.h"
#include "SRmap.h"
#include "SRcoord.h"
#include "SRmaterial.h"
#include "SRbasis.h"
#include "SRnode.h"
#include "SRedge.h"
#include "SRface.h"
#include "SRelement.h"
#include "SRconstraint.h"
#include "SRforce.h"
#include "SRfile.h"
#include "SRerrorCheck.h"

enum SRentityType{ nodeType, edgeType, faceType, elementType };

#define MAXNODEFACEOWNERS 101

class SRElementData
{
	friend class SRmodel;

public:
	SRdoubleMatrix& Getdbdx()
	{ 
		return dbdx; 
	}
	SRdoubleMatrix& Getdbdy()
	{
		return dbdy; 
	}
	SRdoubleMatrix& Getdbdz()
	{
		return dbdz; 
	}
	double* GetIntwt()
	{ 
		return intWt.d; 
	};
	double* GetBasisVec()
	{ 
		return basisVec.d;
	};
	double* Getdbasisdr()
	{
		return dbasisdr.d; 
	};
	double* Getdbasisds()
	{ 
		return dbasisds.d; 
	};
	double* Getdbasisdt()
	{
		return dbasisdt.d;
	};
	double* GetEnfdForceVec()
	{ 
		return enfdForceVec.d;
	};
	SRdoubleVector& getElementStiffness()
	{ 
		return elementStiffness;
	};

private:
	SRdoubleMatrix dbdx, dbdy, dbdz;
	SRdoubleVector basisVec;
	SRdoubleVector dbasisdr;
	SRdoubleVector dbasisds;
	SRdoubleVector dbasisdt;
	SRdoubleVector intWt;
	SRdoubleVector enfdForceVec;
	SRdoubleVector elementStiffness;
};

class SRFaceForceGroup
{
public:
	SRFaceForceGroup(){ nfaceFunMax = 0; fromNodalForces = false; };
	void addFace(SRface* face)
	{
		faceIds.PushBack(face->GetId());
		if (face->GetNumNodesTotal() > nfaceFunMax)
			nfaceFunMax = face->GetNumNodesTotal();
	};

	SRintVector faceIds;
	SRintVector nodeIds;
	int nfaceFunMax;
	SRdoubleMatrix ForceDof;
	SRdoubleVector smoothStiff;
	SRintMatrix faceFunLoc;
	bool fromNodalForces;
};


class SRmodel
{
public:

	SRmodel();
	void FillGlobalFaces();

	int GetNumNodes(){ return nodes.GetNum(); };
	SRnode* GetNode(int i){ return nodes.GetPointer(i); };

	int GetNumEdges(){ return edges.GetNum(); };
	SRedge* GetEdge(int i){ return edges.GetPointer(i); };

	int GetNumFaces(){ return faces.GetNum(); };
	SRface* GetFace(int i){ return faces.GetPointer(i); };

	int GetNumElements() { return elements.GetNum(); };
	SRelement* GetElement(int i){ return elements.GetPointer(i); };

	int GetNumMaterials() { return materials.GetNum(); };
	SRmaterial* GetMaterial(int i){ return materials.GetPointer(i); };

	int GetNumConstraints(){ return constraints.GetNum(); };
	SRconstraint* GetConstraint(int i){ return constraints.GetPointer(i); };

	int GetNumForces(){ return forces.GetNum(); };
	int GetNumAllocatedForces(){ return forces.GetNumAllocated(); };
	SRforce* GetForce(int i){ return forces.GetPointer(i); };
	void freeForce(int i){ forces.Free(i); };
	void packForces(){ forces.packNulls(); };

	int GetNumCoords(){ return Coords.GetNum(); };
	SRcoord* GetCoord(int i){ return Coords.GetPointer(i); };

	SRthermalForce* GetThermalForce() { return thermalForce; };
	void allocateThermalForce() { thermalForce = ALLOCATEMEMORY SRthermalForce; };


	void CreateElemEdges(SRelement *elem, int nnodes, int inputNodes[]);
	void CreateElem(int id, int userid, int nnodes, int nodes[], SRmaterial* mat);

	void PutNumNodeFaces(int i, int n){ numNodeFaces.Put(i, n); };
	int GetNumNodeFaces(int i){ return numNodeFaces[i]; };
	void PutNodeFace(int i, int j, int f){ nodeFaces.Put(i, j, f); };
	int GetNodeFace(int i, int f){ return nodeFaces.Get(i, f); };
	void PutNumNodeEdges(int i, int n){ numNodeEdges.Put(i, n); };
	int GetNumNodeEdges(int i){ return numNodeEdges[i]; };
	int GetNodeEdge(int i, int ej){ return nodeEdges.Get(i, ej); };
	void PutNodeEdge(int i, int j, int ej){ nodeEdges.Put(i, j, ej); };

	bool getAllFacesFlat(){ return allFacesFlat; };
	void setAllFacesFlat(bool tf){ allFacesFlat = tf; };

	double GetSize(){ return size; };
	void SetSize(double s){ size = s; };
	bool checkOrphanNode(int id){ return GetNode(id)->isOrphan(); };

	void CleanUp(bool partial = false);
	void allocateElementData();
	void FreeElementData();
	void allocateSmallElementData(bool anyLcsEnfd = false);
	void FreeSmallElementData();
	SRElementData* GetElementData();
	void setSolutionVector(double *s);
	double *getSolutionVector();
	double GetDisplacementCoeff(int gfun, int dof);
	void PutFunctionEquation(int fun, int dof, int id);
	int GetFunctionEquation(int fun, int dof);
	double GetEnforcedDisplacement(int fun, int dof);
	void PutEnforcedDisplacement(int fun, int dof, double e);
	int getMaxNodeUid(){ return maxNodeUid; };
	void setMaxNodeUid(int i){ maxNodeUid = i; };
	bool getPartialFlatten(){ return partialFlatten; };
	void setPartialFlatten(bool tf){ partialFlatten = tf; };
	int getMaxFlattenedElUid(){ return maxFlattenedElUId; };
	void setMaxFlattenedElUid(int i){ maxFlattenedElUId = i; };
	double getMaxFlattened(){ return maxFlattened; };
	void setMaxFlattened(double i){ maxFlattened = i; };
	int GetmaxNumElementFunctions(){ return maxNumElementFunctions; };
	void setmaxNumElementFunctions(int i){ maxNumElementFunctions = i; };
	double getVolume(){ return volume; };
	void setVolume(double v){ volume = v; };
	void updateVolume(double v){ volume += v; };
	double* getEnforcedVec(){ return elData.enfdForceVec.d; };
	int getNumVolumeForces(){ return volumeForces.GetNum(); };
	SRvolumeForce* getVolumeForce(int v){ return volumeForces.GetPointer(v); };
	void resizeConstraints(int n){ constraints.Allocate(n); };
	int GetSmoothFunctionEquation(int fun);
	void PutSmoothFunctionEquation(int fun, int v);
	void checkElementMapping();
	int elemFaceFind(SRelement* elem, int nv[], int gno[]);

	bool getAllMatsHaveAllowable(){ return allMatsHaveAllowable; };
	double getMaxAllowableAnyActiveMat();
	void AllocateDofVectors();
	void AllocateSmoothFunctionEquations(int n);
	void allocateElements(int nel);
	void allocateFaces(int n);
	void allocateEdges(int n);
	void allocateNodes(int n);
	void allocateConstraints(int n);
	void allocateForces(int n);
	void allocateVolumeForces(int n);
	void allocateFaceForceGroups(int n);
	int getNumFaceForceGroups();
	SRFaceForceGroup* addFaceForceGroup();
	SRFaceForceGroup* getFaceForceGroup(int i);
	void freeConstraints(int n);
	void packConstraints();
	SRcoord* addCoord();
	SRmaterial* addMat();
	SRnode* addNode();
	SRedge* addEdge();
	SRface* addFace();
	SRelement* addElement();
	SRconstraint* addConstraint();
	SRforce* addForce();
	SRvolumeForce* addVolumeForce() { return volumeForces.Add(); };
	void setRepFileName(char* s);
	void allocateNodeEdges(int n);
	void allocateNodeFaces(int n);
	void freeNodeEdges();
	void freeNodeFaces();
	void SetsimpleElements();
	bool UseSimpleElements();
	void NumberGlobalFunctions();
	int GetNumFunctions();
	void initializeErrorCheck();
	bool checkForSmallMaxStress();
	void setupErrorCheck(bool finalAdapt, double stressMaxIn, double ErrorToleranceIn, int maxPorderIn, int maxPJumpIn, int maxPorderLowStressIn);
	int CheckAutoSacrificialElements();
	void CleanUpErrorCheck();
	bool FindNextP(SRelement* elem, int& p);
	double FindElementError(SRelement* elem);
	void SetLowStressTolerance(double tol);
	void SetLowStressToleranceFinalAdapt(double tol);
	bool PreProcessPenaltyConstraints();
	void ProcessConstraints();
	int NumberEquations();
	double* GetElementStiffnessVector();

	//1 instance of each utility class:
	SRbasis basis;
	SRmath math;
	SRmap map;
	SRerrorCheck errorChecker;

	SRfile repFile;
	SRfile logFile;
	SRfile outFile;

private:
	double size;
	double volume;

	double* solutionVector;

	SRintVector numNodeEdges;
	SRintMatrix nodeEdges;
	SRintVector numNodeFaces;
	SRintMatrix nodeFaces;
	SRpointerVector <SRnode> nodes;
	SRpointerVector <SRedge> edges;
	SRpointerVector <SRface> faces;
	SRpointerVector <SRconstraint> constraints;
	SRpointerVector <SRcoord> Coords;
	SRpointerVector <SRmaterial> materials;
	SRpointerVector <SRforce> forces;
	SRpointerVector <SRvolumeForce> volumeForces;
	SRpointerVector <SRelement> elements;
	SRthermalForce* thermalForce;
	SRpointerVector <SRFaceForceGroup> faceForceGroups;
	SRpointerVector <SRElProperty> elProps;
	SRintMatrix functionEquations;
	SRdoubleMatrix enforcedDisplacements;
	SRintVector smoothFunctionEquations;

	//scratch space needed by elements:
	SRElementData elData;
	int maxNumElementFunctions;
	int numFunctions;
	int numEquations;
	
	bool allFacesFlat;
	bool partialFlatten;

	double maxFlattened;

	int maxNodeUid;
	int maxFlattenedElUId;
	int maxPinModel;
	double errorMax;
	double errorSmoothRawAtMax;
	double errorFaceJumpAtMax;
	bool allMatsHaveAllowable;
	double maxAllowableAnyActiveMat;
	bool anyBricks;
	bool anyWedges;
	bool simpleElements;
};

#endif //!defined(SRMODEL_INCLUDED)
