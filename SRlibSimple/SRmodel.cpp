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
// SRmodel.cpp: implementation of the SRmodel class.
//
//////////////////////////////////////////////////////////////////////

#include <stdlib.h>
#include "SRmodel.h"

//One global instance of "model""
SRmodel model;

//////////////////////////////////////////////////////////////////////

SRmodel::SRmodel()
{
	allFacesFlat = false;
	partialFlatten = false;

	maxFlattenedElUId = 0;
	maxPinModel = 2;
	
	size = 0.0;
	maxFlattened = 0.0;

	thermalForce = NULL;
	math.Setup();

	materials.Allocate(MAXNMAT);
	Coords.Allocate(MAXNCOORD);

	anyBricks = false;
	anyWedges = false;

	simpleElements = false;
};

void SRmodel::FillGlobalFaces()
{
	//fill global faces in the model
	//note:
		//fills class variable faces
		//faces are created with node numbers, elementOwners, elementLocalFace numbers,
		//and localEdges

	//make sure faces were allocated:
	if (faces.GetNum() == 0)
	{
		//worst case allocation:
		//worst case for number of faces:
		int nface;
		int nel = GetNumElements();
		if (anyBricks)
			nface = 6 * nel;
		else if (anyWedges)
			nface = 5 * nel;
		else
			nface = 4 * nel;
		faces.Allocate(nface);
	}

	int nnode = GetNumNodes();
	numNodeFaces.Allocate(nnode);
	nodeFaces.Allocate(nnode, MAXNODEFACEOWNERS);
	
	int lface, n1, n2, n3, n4, gn1, gn2, gn3, gn4, gface;
	SRface* face;
	SRlocalFace* ellocFace;
	int direction;
	SRelement* elem;
	for (int el = 0; el < elements.GetNum(); el++)
	{
		elem = GetElement(el);
		int nlocface = elem->GetNumLocalFaces();
		for (lface = 0; lface < nlocface; lface++)
		{
			elem->GetFaceNodes(lface, n1, n2, n3, n4);
			gface = SRfaceUtil::GlobalFaceMatch(gn1, gn2, gn3, gn4, n1, n2, n3, n4);
			ellocFace = GetElement(el)->GetLocalFace(lface);
			if (gface == -1)
			{
				//face not found, add new global face:
				gface = faces.GetNum();
				face = faces.Add();
				face->Create(gface,n1, n2, n3, n4);
				face->PutElementOwner(0, el);
				face->PutElementLocalFace(0, lface);
				//assign local edges to the new face:
				int gej = SRedgeUtil::GlobalEdgeMatch(n1, n2, direction);
				face->AssignLocalEdge(0, gej, direction);
				gej = SRedgeUtil::GlobalEdgeMatch(n2, n3, direction);
				face->AssignLocalEdge(1, gej, direction);
				if (n4 == -1)
				{
					gej = SRedgeUtil::GlobalEdgeMatch(n1, n3, direction);
					face->AssignLocalEdge(2, gej, direction);
				}
				else
				{
					gej = SRedgeUtil::GlobalEdgeMatch(n4, n3, direction);
					face->AssignLocalEdge(2, gej, direction);

					gej = SRedgeUtil::GlobalEdgeMatch(n1, n4, direction);
					face->AssignLocalEdge(3, gej, direction);
				}
				gn1 = 0;
				gn2 = 1;
				gn3 = 2;
				if (n4 == -1)
					gn4 = -1;
				else
					gn4 = 3;
			}
			else
			{
				face = GetFace(gface);
				face->PutElementOwner(1, el);
				face->PutElementLocalFace(1, lface);
			}
			ellocFace->PutGlobalNodeOrder(gn1, 0);
			ellocFace->PutGlobalNodeOrder(gn2, 1);
			ellocFace->PutGlobalNodeOrder(gn3, 2);
			if (gn4 != -1)
				ellocFace->PutGlobalNodeOrder(gn4, 3);
			ellocFace->PutGlobalFaceId(gface);
		}
	}

	//assign boundaryFaceId to all edges of boundary faces;
	//overwrite the numnodefaces and nodefaces arrays with boundary faces only
	int nbdryface = 0;
	numNodeFaces.Zero();
	for (int f = 0; f < faces.GetNum(); f++)
	{
		SRface* face = GetFace(f);
		if (!face->IsBoundaryFace())
			continue;
		nbdryface++;
		for (int e = 0; e < face->GetNumLocalEdges(); e++)
		{
			SRedge* edge = face->GetEdge(e);
			edge->PutBoundaryFaceId(f,e);
		}
		for (int i = 0; i < face->GetNumNodes(); i++)
		{
			int nid = face->GetNodeId(i);
			int n = numNodeFaces.Get(nid);
			PutNodeFace(nid, n, f);
			numNodeFaces.PlusAssign(nid, 1);
		}
	}
}

void SRmodel::CreateElemEdges(SRelement *elem, int nnodes, int inputNodes[])
{
	//Fill contribution of this element to global edges
	//input:
		//elem = pointer to eleme
		//nnodes = number of input nodes for this element
		//inputNodes = list of input nodes for this element
	//note: 
		//this routine also assigns the global edge ids to the element's local edges
		//fills class variable edges
		//new edges created with node numbers (corner and midside), and localEdges

	int ncorner = 4;
	int nej = 6;
	if (elem->GetType() == wedge)
	{
		ncorner = 6;
		nej = 9;
	}
	else if (elem->GetType() == brick)
	{
		ncorner = 8;
		nej = 12;
	}

	for (int lej = 0; lej < nej; lej++)
	{
		int n1, n2;
		elem->GetEdgeNodeIds(lej, n1, n2);
		int mid = inputNodes[lej + ncorner];
		int direction;
		SRedge* edge;
		int gej = SRedgeUtil::GlobalEdgeMatch(n1, n2, direction);
		if (gej == -1)
		{
			//edge not found:
			gej = GetNumEdges();
			edge = edges.Add();
			edge->Create(n1, n2, mid, gej);
			direction = 1;
			SRnode* node = GetNode(mid);
			node->SetAsMidside(gej);
			node->SetFirstElementOwner(elem->GetId());
		}
		elem->AssignLocalEdge(lej, gej, direction);
	}
}

void SRmodel::CleanUp(bool partial)
{
	//miscelaneous memory and disk cleanup at end of run

	FreeElementData();

	if (partial)
		return;

	nodes.Free();
	edges.Free();
	faces.Free();
	constraints.Free();
	Coords.Free();
	materials.Free();
	forces.Free();
	for (int e = 0; e < GetNumElements(); e++)
	{
		SRelement *elem = GetElement(e);
		elem->Cleanup();
	}
	elements.Free();
	volumeForces.Free();
	if (thermalForce != NULL)
	{
		DELETEMEMORY thermalForce;
		thermalForce = NULL;
	}
	functionEquations.Free();
	enforcedDisplacements.Free();
	smoothFunctionEquations.Free();
}

void SRmodel::allocateElementData()
{
	//allocate memory intensive data related to basis functions
	//for use by element stiffness routines

	//Note:
		//allocate and these data before element stiffness loop

	//worst case size scratch space for elements
	int nintMax = 0;
	for (int i = 0; i < elements.GetNum(); i++)
	{
		SRelement* elem = GetElement(i);
		int nint = math.CountGaussPoints(elem);
		if (nint > nintMax)
			nintMax = nint;
	}
	int nfunMax = maxNumElementFunctions;
	int neq = 3 * nfunMax;
	int stiffnessLength = (neq)*(neq + 1) / 2;


	SRElementData* eld = &elData;
	eld->dbdx.Allocate(nintMax, nfunMax);
	eld->dbdy.Allocate(nintMax, nfunMax);
	eld->dbdz.Allocate(nintMax, nfunMax);
	eld->elementStiffness.Allocate(stiffnessLength);
}

void SRmodel::FreeElementData()
{
	//free memory intensive data related to basis functions
	//for use by element stiffness routines

	//Note:
		//free these data after element stiffness loop
	SRElementData* eld = &elData;
	eld->dbdx.Free();
	eld->dbdy.Free();
	eld->dbdz.Free();
}

void SRmodel::allocateSmallElementData(bool anyLcsEnfd)
{
	//allocate less memory intensive data 
	//for use by element stiffness and stress routines

	//note:
	//these data must persist through element stiffness, error checking, and postprocessing,
	//so allocate at top of adaptivity loop.
	//but this must be allocated after call to NumberGlobalFunctions so nfunmax is known
	//data must not be freed until after final postprocessing outside adapt loop

	//worst case size scratch space for elements
	int nintMax = 0;
	for (int i = 0; i < elements.GetNum(); i++)
	{
		SRelement* elem = GetElement(i);
		int nint = math.CountGaussPoints(elem);
		if (nint > nintMax)
			nintMax = nint;
	}
	int nfunMax = maxNumElementFunctions;

	SRElementData* eld = &elData;
	eld->basisVec.Allocate(nfunMax);
	eld->dbasisdr.Allocate(nfunMax);
	eld->dbasisds.Allocate(nfunMax);
	eld->dbasisdt.Allocate(nfunMax);
	eld->intWt.Allocate(nintMax);

	if (anyLcsEnfd)
		eld->enfdForceVec.Allocate(numEquations);

}

void SRmodel::FreeSmallElementData()
{
	//note:
	//these data must persist through element stiffness, error checking, and postprocessing,
	//so free after final postprocessing outside adapt loop

	SRElementData* eld = &elData;
	eld->basisVec.Free();
	eld->dbasisdr.Free();
	eld->dbasisds.Free();
	eld->dbasisdt.Free();
	eld->intWt.Free();
	eld->enfdForceVec.Free();
}

SRElementData* SRmodel::GetElementData()
{
	return &elData;
};

double SRmodel::GetDisplacementCoeff(int gfun, int dof)
{
	//look up displacement coefficient for a global function and dof
	//input:
	//gfun = global function number
	//dof = degree of freedom
	//return:
	//displacement coefficient

	int eq = GetFunctionEquation(gfun, dof);
	if (eq >= 0)
		return solutionVector[eq];
	else
		return GetEnforcedDisplacement(gfun, dof);
}

SRcoord* SRmodel::addCoord()
{
	return Coords.Add();
}

SRmaterial* SRmodel::addMat()
{
	return materials.Add();
}

SRnode* SRmodel::addNode()
{
	return nodes.Add();
}

SRedge* SRmodel::addEdge()
{
	return edges.Add();
}

SRface* SRmodel::addFace()
{
	return faces.Add();
}

SRelement* SRmodel::addElement()
{
	return elements.Add();
}

SRconstraint* SRmodel::addConstraint()
{
	return constraints.Add();
}

SRforce* SRmodel::addForce()
{
	return forces.Add();
}

void SRmodel::AllocateDofVectors()
{
	functionEquations.Allocate(numFunctions, 3);
	enforcedDisplacements.Allocate(numFunctions, 3);
}

void SRmodel::AllocateSmoothFunctionEquations(int n)
{
	smoothFunctionEquations.Allocate(n);
}

void SRmodel::allocateElements(int nel)
{
	elements.Allocate(nel);
}

void SRmodel::allocateEdges(int nel)
{
	//do nothing if edges were already allocated, e.g. in allocateElements
	if (edges.num < nel)
		edges.Allocate(nel);
}
void SRmodel::allocateFaces(int n)
{
	faces.Allocate(n);
}
void SRmodel::allocateNodes(int n)
{
	nodes.Allocate(n);
}
void SRmodel::allocateConstraints(int n)
{
	constraints.Allocate(n);
}
void SRmodel::allocateForces(int n)
{
	forces.Allocate(n);
}

void SRmodel::allocateVolumeForces(int n)
{
	volumeForces.Allocate(n);
}

void SRmodel::freeConstraints(int n)
{
	constraints.Free(n);
}
void SRmodel::packConstraints()
{
	constraints.packNulls();
}

void SRmodel::allocateFaceForceGroups(int n)
{
	faceForceGroups.Allocate(n);
}

int SRmodel::getNumFaceForceGroups()
{
	return faceForceGroups.GetNum();
}

SRFaceForceGroup* SRmodel::addFaceForceGroup()
{
	return faceForceGroups.Add();
}

SRFaceForceGroup* SRmodel::getFaceForceGroup(int i)
{
	return faceForceGroups.GetPointer(i);
}

void SRmodel::setRepFileName(char* s)
{
	repFile.setFileName(s);
}

void SRmodel::allocateNodeEdges(int n)
{
	numNodeEdges.Allocate(n);
	nodeEdges.Allocate(n, 20);
}

void SRmodel::allocateNodeFaces(int n)
{
	numNodeFaces.Allocate(n);
	nodeFaces.Allocate(n, 20);
}

void SRmodel::freeNodeEdges()
{
	numNodeEdges.Free();
	nodeEdges.Free();
}

void SRmodel::freeNodeFaces()
{
	numNodeFaces.Free();
	nodeFaces.Free();
}
void SRmodel::setSolutionVector(double *s)
{
	solutionVector = s; 
};
double* SRmodel::getSolutionVector()
{
	return solutionVector;
};

void SRmodel::PutFunctionEquation(int fun, int dof, int id)
{
	functionEquations.Put(fun, dof, id);
};

int SRmodel::GetFunctionEquation(int fun, int dof)
{
	return functionEquations.Get(fun, dof);
};
double SRmodel::GetEnforcedDisplacement(int fun, int dof)
{
	return enforcedDisplacements.Get(fun, dof);
};

void SRmodel::PutEnforcedDisplacement(int fun, int dof, double e)
{
	enforcedDisplacements.Put(fun, dof, e);
};

int SRmodel::GetSmoothFunctionEquation(int fun)
{
	return smoothFunctionEquations.Get(fun);
};

void SRmodel::PutSmoothFunctionEquation(int fun,int v)
{
	smoothFunctionEquations.Put(fun,v);
};

void SRmodel::CreateElem(int id, int userid, int nnodes, int nodes[], SRmaterial* mat)
{
	SRelement* elem = GetElement(id);
	elem->Create(userid, nnodes, nodes, mat);
}

void SRmodel::SetsimpleElements()
{
	simpleElements = true;
};

bool SRmodel::UseSimpleElements()
{
	return simpleElements;
};

void SRmodel::NumberGlobalFunctions()
{
	//assign global function numbers to all functions in model

	SRnode* node;
	SRedge* edge;
	SRface* face;
	SRelement* elem;
	int i, n, nint;
	int fun;
	int j, k, pej;
	int lfun, gfun;

	//allocate space for function numbers for edges, faces, and elements:
	for (i = 0; i < GetNumEdges(); i++)
	{
		edge = GetEdge(i);
		n = 1 + edge->GetPorder();
		edge->AllocateGlobalFunctionNumbers(n);
	}
	for (i = 0; i < GetNumFaces(); i++)
	{
		face = GetFace(i);
		n = basis.CountFaceTotalFunctions(face, j);
		face->AllocateGlobalFunctionNumbers(n);
	}
	for (i = 0; i < GetNumElements(); i++)
	{
		elem = GetElement(i);
		n = basis.CountElementFunctions(elem, nint);
		elem->AllocateGlobalFunctionNumbers(n);
	}

	//for nodes, assign a function number unless it is a midside node
	//or orphan:
	fun = 0;
	for (i = 0; i < GetNumNodes(); i++)
	{
		node = GetNode(i);
		if (!node->isMidSide() && !node->isOrphan())
		{
			node->PutGlobalFunctionNumber(fun);
			fun++;
		}
	}

	//edges:
	//assign p2 funs 1st for easier mapping of p2 soln to adapted solns (e.g. itsolv):
	for (i = 0; i < GetNumEdges(); i++)
	{
		edge = GetEdge(i);
		edge->AssignGlobalFunctionNumbers(fun, 2, 2);
		//assign the same function number to the corresponding global node.
		//but fun was incremented in edge->AssignGlobalFunctionNumbers so subtract 1:
		SRnode* node = GetNode(edge->GetMidNodeId());
		node->PutGlobalFunctionNumber(fun - 1);
	}
	for (i = 0; i < GetNumEdges(); i++)
	{
		edge = GetEdge(i);
		edge->AssignGlobalFunctionNumbers(fun, 3, 8);
	}

	//faces:
	int nf = GetNumFaces();
	int nn, ne;
	for (i = 0; i < nf; i++)
	{
		lfun = 0;
		face = GetFace(i);
		//corner functions:
		nn = face->GetNumNodes();
		for (j = 0; j < nn; j++)
		{
			gfun = face->GetNode(j)->GetGlobalFunctionNumber();
			face->PutGlobalFunctionNumber(lfun, gfun);
			lfun++;
		}

		//edge p2 functions:
		ne = face->GetNumLocalEdges();
		for (j = 0; j < ne; j++)
		{
			edge = face->GetEdge(j);
			gfun = edge->GetGlobalFunctionNumber(2);
			face->PutGlobalFunctionNumber(lfun, gfun);
			lfun++;
		}

		//edge higher functions:
		ne = face->GetNumLocalEdges();
		for (j = 0; j < ne; j++)
		{
			edge = face->GetEdge(j);
			pej = edge->GetPorder();
			for (k = 3; k <= pej; k++)
			{
				gfun = edge->GetGlobalFunctionNumber(k);
				face->PutGlobalFunctionNumber(lfun, gfun);
				lfun++;
			}
		}

		n = face->GetNumGlobalFunctions();
		for (j = lfun; j < n; j++)
		{
			face->PutGlobalFunctionNumber(lfun, fun);
			lfun++;
			fun++;
		}
	}
	//elements:
	maxNumElementFunctions = 0;
	for (i = 0; i < GetNumElements(); i++)
	{
		lfun = 0;
		elem = GetElement(i);
		//nodes:
		n = elem->GetNumNodes();
		for (j = 0; j < n; j++)
		{
			gfun = elem->GetNode(j)->GetGlobalFunctionNumber();
			elem->PutGlobalFunctionNumbers(lfun, gfun);
			lfun++;
		}

		//edge p2 functions:
		for (j = 0; j < elem->GetNumLocalEdges(); j++)
		{
			edge = elem->GetEdge(j);
			gfun = edge->GetGlobalFunctionNumber(2);
			elem->PutGlobalFunctionNumbers(lfun, gfun);
			lfun++;
		}

		//edge higher functions:
		for (j = 0; j < elem->GetNumLocalEdges(); j++)
		{
			edge = elem->GetEdge(j);
			pej = edge->GetPorder();
			for (k = 3; k <= pej; k++)
			{
				gfun = edge->GetGlobalFunctionNumber(k);
				elem->PutGlobalFunctionNumbers(lfun, gfun);
				lfun++;
			}
		}

		//face functions:
		for (j = 0; j < elem->GetNumLocalFaces(); j++)
		{
			face = elem->GetFace(j);
			n = basis.CountFaceTotalFunctions(face, nint);
			//skip number of edge and node functions on face:
			for (k = (n - nint); k < n; k++)
			{
				gfun = face->GetGlobalFunctionNumber(k);
				elem->PutGlobalFunctionNumbers(lfun, gfun);
				lfun++;
			}
		}

		//internal functions:
		n = basis.CountElementFunctions(elem, nint);
		if (n > maxNumElementFunctions)
			maxNumElementFunctions = n;
		for (j = 0; j < nint; j++)
		{
			elem->PutGlobalFunctionNumbers(lfun, fun);
			lfun++;
			fun++;
		}
	}

	numFunctions = fun;
	setmaxNumElementFunctions(maxNumElementFunctions);
}

int SRmodel::GetNumFunctions()
{
	return numFunctions;
}

double SRmodel::getMaxAllowableAnyActiveMat()
{
	return maxAllowableAnyActiveMat;
};

void SRmodel::initializeErrorCheck()
{
	errorChecker.Initialize();
}

bool SRmodel::checkForSmallMaxStress()
{
	return errorChecker.getSmallMaxStressDetected();
}

void SRmodel::setupErrorCheck(bool finalAdapt, double stressMaxIn, double ErrorToleranceIn, int maxPorderIn, int maxPJumpIn, int maxPorderLowStressIn)
{
	errorChecker.SetUp(finalAdapt, stressMaxIn, ErrorToleranceIn, maxPorderIn, maxPJumpIn, maxPorderLowStressIn);
}

int SRmodel::CheckAutoSacrificialElements()
{
	return errorChecker.AutoSacrificialElements();
}

void SRmodel::CleanUpErrorCheck()
{
	errorChecker.CleanUp();
}

bool SRmodel::FindNextP(SRelement* elem, int& p)
{
	return errorChecker.FindRequiredPOrder(elem, p);
}
double SRmodel::FindElementError(SRelement* elem)
{
	return errorChecker.FindError(elem);
}

void SRmodel::SetLowStressTolerance(double tol)
{
	errorChecker.SetLowStressTol(tol);
}

void SRmodel::SetLowStressToleranceFinalAdapt(double tol)
{
	errorChecker.setLowStressTolFinalAdapt(tol);
}

bool SRmodel::PreProcessPenaltyConstraints()
{
	//see if any penalty constraints are needed; if so, calibrate the penalty constant

	bool anyLcsEnfd = false;

	SRconstraint* con;
	//fill faceLcsCoonstraints for elements with non-gcs constraints on boundary faces:
	for (int i = 0; i < model.GetNumConstraints(); i++)
	{
		con = model.GetConstraint(i);
		if (con->GetType() != faceCon)
			continue;
		if (!con->isGcs())
		{
			int f = con->GetEntityId();
			SRface* face = model.GetFace(f);
			if (!face->IsBoundaryFace())
				ERROREXIT; //not boundary face, non-gcs not supported
			int eId = face->GetElementOwner(0);
			SRelement* elem = model.GetElement(eId);
			int lface = face->GetElementLocalFace(0);
			elem->CalibratePenalty();
			elem->AddfaceLCSConstraint(lface);
			if (con->hasEnforcedDisp())
				anyLcsEnfd = true;

		}
	}

	//fill nodeLcsCoonstraints for elements with non-gcs constraints on nodes:
	for (int i = 0; i < model.GetNumConstraints(); i++)
	{
		con = model.GetConstraint(i);
		if (con->GetType() != nodalCon)
			continue;
		if (!con->isGcs())
		{
			int n = con->GetEntityId();
			bool found = false;
			for (int f = 0; f < model.GetNumFaces(); f++)
			{
				SRface* face = model.GetFace(f);
				if (!face->IsBoundaryFace())
					continue;
				int eid = face->GetElementOwner(0);
				SRelement* elem = model.GetElement(eid);
				//find local node corresponding to node:
				for (int l = 0; l < elem->GetNumNodesTotal(); l++)
				{
					if (elem->getNodeOrMidNodeId(l) == n)
					{
						elem->AddNodeLcSConstraints(l);
						elem->CalibratePenalty();
						found = true;
						break;
					}
				}
				if (found)
					break;
			}
			if (!found)
				ERROREXIT;//lcs constraint not on boundary- not supported
			if (con->hasEnforcedDisp())
				anyLcsEnfd = true;
		}
	}

	return anyLcsEnfd;
}

void SRmodel::ProcessConstraints()
{
	//process constraints
	//notes:
		//number functionEquations(fun,dof) as negative for a gcs constrained global dof.
		//the number is -1 if there is no enforced displacement, else
		//it is the number in the enforcedDisplacements matrix where the enforced displacement is stored
		//class variables functionEquations and enforcedDisplacements are updated
		//lcs constraints are skipped because they are handled via penalty

	int i, dof, gfun;
	SRconstraint* con;
	SRnode* node;

	AllocateDofVectors();

	for (i = 0; i < GetNumConstraints(); i++)
	{
		con = GetConstraint(i);

		//skip lcs constraint, they are handled via penalty:
		if (!con->isGcs())
			continue;
		int eid = con->GetEntityId();
		if (con->GetType() == nodalCon)
		{
			node = model.GetNode(eid);
			gfun = node->GetGlobalFunctionNumber();
			for (dof = 0; dof < 3; dof++)
			{
				if (con->IsConstrainedDof(dof))
				{
					model.PutFunctionEquation(gfun, dof, -1);
					if (con->hasEnforcedDisp())
					{
						double enfd = con->getDisp(0, dof);
						model.PutEnforcedDisplacement(gfun, dof, enfd);
					}
				}
			}
		}
		else if (con->GetType() == faceCon)
		{
			con->ProcessFaceConstraint();
		}
	}

	for (i = 0; i < model.GetNumElements(); i++)
	{
		SRelement* elem = model.GetElement(i);
		for (int f = 0; f < elem->GetNumFunctions(); f++)
		{
			int gfun = elem->GetFunctionNumber(f);
			for (dof = 0; dof < 3; dof++)
			{
				if (model.GetFunctionEquation(gfun, dof) < 0)
				{
					elem->SetConstrained();
					break;
				}
			}
			if (elem->isConstrained())
				break;
		}
	}
}

int SRmodel::NumberEquations()
{
	//global equation numbering
	//note:
		//fills class variable functionEquations
		//ProcessConstraints has to be called 1st. 
		//it assigns a negative number to global dofs that are constrained
		//this routine will only assign an equation number to unconstrained
		//global dofs

	int eq = 0;
	for (int gfun = 0; gfun < model.GetNumFunctions(); gfun++)
	{
		for (int dof = 0; dof < 3; dof++)
		{
			//unconstrained dofs have temporarily been assigned to 0.
			//constrained dofs are assigned a negative number
			int eqt = model.GetFunctionEquation(gfun, dof);
			if (eqt == 0)
			{
				model.PutFunctionEquation(gfun, dof, eq);
				eq++;
			}
		}
	}
	numEquations = eq;
	return numEquations;
}

void SRmodel::checkElementMapping()
{

	//check mapping of all elements at all gauss points needed this p-pass.
	//straighten edges and mark associated elements as sacrificial if mapping is bad.

	int nel = model.GetNumElements();
	bool anyFail = false;
	int e;
	for (e = 0; e < nel; e++)
	{
		SRelement* elem = model.GetElement(e);
		elem->FillMappingNodes();
		if (!elem->testMapping())
			anyFail = true;
	}

	if (!anyFail)
		return;

	for (e = 0; e < nel; e++)
	{
		SRelement* elem = model.GetElement(e);
		elem->checkForStraightenedEdges();
	}

	for (e = 0; e < nel; e++)
	{
		SRelement* elem = model.GetElement(e);
		if (!elem->isFlattened())
			continue;
		elem->FillMappingNodes();
		if (!elem->testMapping())
		{
			//LOGPRINT("checkElementMapping elem failed: %d", elem->userId);
			anyFail = true;
		}
	}
	//flattening of an element may have invalidated an adjacent element
	if (!anyFail)
		return;

	//flattening of an element may have invalidated an adjacent element
	for (e = 0; e < nel; e++)
	{
		SRelement* elem = model.GetElement(e);
		if (!elem->isFlattened())
			continue;
		elem->FillMappingNodes();
		if (!elem->testMapping())
			anyFail = true;
	}
	if (anyFail)
	{
		model.setPartialFlatten(false);
		for (e = 0; e < nel; e++)
		{
			SRelement* elem = model.GetElement(e);
			if (!elem->isFlattened())
				continue;
			elem->FillMappingNodes();
			if (!elem->testMapping())
				anyFail = true;
		}
		for (e = 0; e < nel; e++)
		{
			SRelement* elem = model.GetElement(e);
			elem->checkForStraightenedEdges();
		}
	}
}

int SRmodel::elemFaceFind(SRelement* elem, int nv[], int gno[])
{
	//find the face on an element that matches nodes in nv array
	//input:
		//elem = pointer to element
		//nv = corner nodes of faces. nv[3] = -1 for tri face
	//output:
		//gno = global node orders of the nodes on the face corresponding to nv
	//return
		//id of the global faces in the model with corner nodes nv; -1 if no match
	for (int l = 0; l < elem->GetNumLocalFaces(); l++)
	{
		SRface* face = elem->GetFace(l);
		if (SRfaceUtil::FaceMatch(nv, face->getNodeIds(), gno))
			return elem->GetLocalFaceGlobalId(l);
	}
	return -1;
}

double* SRmodel::GetElementStiffnessVector()
{
	//look up the element stiffness Vector for a processor number
	//input:
		//processorNum = processor number, 0 if not multi-threaded

	SRElementData* eld = GetElementData();
	return eld->elementStiffness.GetVector();
}









