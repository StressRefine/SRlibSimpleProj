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

#include "SRmodel.h"

extern SRmodel model;

static int nodeBuf[20];

static double* LTStiffbuf;

//////////////////////////////////////////////////////////////////////
//
// globalWrappers.c: function wrappers for calling stressRefine library from c and fortran
//
//////////////////////////////////////////////////////////////////////

//C WRAPPERS:

void SetNumNodes(int numNodes)
{
	//allocate space for nodes
	//input:
		//numNodes = total number of nodes in the mesh
	model.allocateNodes(numNodes);
}

int numberGlobalFunctions()
{
	//calculate number of global basis functions in model
	//return:
		//number of global basis functions
	model.NumberGlobalFunctions();
	return model.GetNumFunctions();
}


int createNode(int uid, double x, double y, double z)
{
	//create a node
	//input:
		//uid = node userId
		//x,y,z = node coordinates in global C.S.

	int id = model.GetNumNodes();
	SRnode* node = model.addNode();
	node->setId(id);
	node->SetUserId(uid);
	SRvec3 pos;
	pos.Assign(x, y, z);
	node->SetPosition(pos);
	return id;
}


void SetNumElements(int numBricks, int numWedges, int numTets)
{
	//allocate space for elements, edges, and faces
	//input:
		//numBricks = number of brick (hexahedral) elements
		//numWedges = number of wedge (pentahedral) elements
		//numTets = number of tetrahedral elements

	//worst case allocation, does not account for sharing. this
	//is not wasteful because only pointers are allocated
	int nedge = 12 * numBricks + 9 * numWedges + 6 * numTets;
	int nface = 6 * numBricks + 5 * numWedges + 4 * numTets;
	int nelem = numBricks + numWedges + numTets;

	model.allocateEdges(nedge);
	model.allocateFaces(nface);
	model.allocateElements(nelem);
	LTStiffbuf = NULL;
}

int createElementIso(int uid, int numnodes, int* nodeids, double E, double nu)
{
	//create an isotropic p-element from the definition
	//of a conventional quadratic element
	//input:
		//uid = user Id of the element
		//numnodes = total number (corner and midside) of nodes in the
			//quadratic element. supported values are 10, 15, 20
			//(tets, wedges, and bricks)
		//nodeIds = list of global ids of the corner and midside nodes
		//E = young's modulus
		//nu = Poisson's ratio
	//return:
		//0 for successful completion
		//-1: unsupported number of nodes (linear element)
		//-2: unsupported number of nodes (general)
		//-3: too many unique materials

	if (numnodes == 4 || numnodes == 6 || numnodes == 8)
		return -1;
	if (numnodes != 10 && numnodes != 15 && numnodes != 20)
		return -2;

	SRmaterial* thisMat = NULL;

	bool newMat = true;
	if (model.GetNumMaterials() != 0)
	{
		//see if this element's material matches the first in the model:
		SRmaterial* mat0 = model.GetMaterial(0);
		double E0, nu0;
		E0 = mat0->getE();
		nu0 = mat0->getNu();
		if (fabs(E0 - E) < TINY && fabs(nu0 - nu) < TINY)
		{
			thisMat = mat0;
			newMat = false;
		}
	}
	if(newMat)
	{
		int nmat = model.GetNumMaterials();
		if (nmat > MAXNMAT - 1)
			return -3;
		thisMat = model.addMat();
		thisMat->assignIso(E, nu);
	}
	thisMat->updateNumElements();

	int nelem = model.GetNumElements();
	SRelement* elem = model.addElement();
	elem->setId(nelem);

	for (int i = 0; i < numnodes; i++)
	{
		SRnode* node = model.GetNode(nodeids[i]);
		if (node->GetFirstElementOwner() == -1)
			node->SetFirstElementOwner(elem->GetId());
	}

	//void Create(int userid, int nnodes, int nodest[], SRmaterial * mat);
	elem->Create(uid, numnodes, nodeids, thisMat);

	model.CreateElemEdges(elem, numnodes, nodeids);

	return 0;
}

void createGlobalFaces()
{
	//create global faces for model
	int nnode = model.GetNumNodes();
	model.allocateNodeFaces(nnode);

	model.FillGlobalFaces();
}

double* CalculateStiffnessMatrix(int id, bool upperTriangle, double* globalEnfdVec)
{
	//calculate global stiffness matrix for an element
	//input:
		//id = element id
		//upperTriangle = true if output stiffness matrix (symmetric) is to be stored as upper trianlge
						//false for lower triangle
		//input:
		//globalEnfdVec = storage for global enforced displacement vector
		//(0 to use internal storage or if enforced displacement not needed)
	//return: element stiffness matrix in format specified by upperTriangle flag
	//note: if globalEnfdVec was provided, it is updated with this Element's contribution

	SRelement* elem = model.GetElement(id);
	int numEq;
	double* kelLower = elem->CalculateStiffnessMatrix(numEq, globalEnfdVec);
	if (upperTriangle)
		return kelLower;
	else
	{
		if (LTStiffbuf == NULL)
		{
			int maxelemEq = 3 * MAXELEMBASISFUNCTIONS;
			int maxElemSize = maxelemEq * (maxelemEq + 1) / 2;
			LTStiffbuf = new double[maxelemEq];
		}
		//copy LT to UT
		int utloc = 0;
		for (int row = 0; row < numEq; row++)
		{
			for (int col = 0; col <= row; col++)
			{
				int ltloc = elem->GetStiffnessLocationAboveDiag(col, row);
				LTStiffbuf[ltloc] = kelLower[utloc];
				utloc++;
			}
		}
		return LTStiffbuf;
	}
}

int CheckAutoSacrificialElements()
{
	//check for sacrificial elements in model (due to point loads, reentrant corners, etc)
	//return:
		//number of sacrificial elements
	return model.CheckAutoSacrificialElements();
}

int GetNumFunctions()
{
	//get number of global functions in model
	//return:
		//number of functions
	return model.GetNumFunctions();
}

void initializeErrorCheck()
{
	//initilialize error checking
	model.initializeErrorCheck();
}

bool checkForSmallMaxStress()
{
	//see if max stress in model is small, e.g. free thermal expansion case
	//return:
		//true if max stress is small
	return model.checkForSmallMaxStress();
}

void setupErrorCheck(bool finalAdapt, double stressMaxIn, double ErrorToleranceIn, int maxPorderIn, int maxPJumpIn, int maxPorderLowStressIn)
{
	//set up for error checking
	//input:
		//finalAdapt = true if this is the final adaptivity pass
		//stressMaxIn = max von Mises stress in model
		//ErrorToleranceIn = relative error tolerance, e.g 0.05 for 5 percent
		//maxPorderIn = max polynomial order allowed in model
		//maxPJumpIn = max change in polynomial order allowed this pass
		//maxPorderLowStressIn = max polynomial order allowed in model for elements with low stress
	model.setupErrorCheck(finalAdapt, stressMaxIn, ErrorToleranceIn, maxPorderIn, maxPJumpIn, maxPorderLowStressIn);
}

void CleanUpErrorCheck()
{
	//perform memory cleanup tasks for error checker
	model.CleanUpErrorCheck();
}
bool FindNextP(int elId, int* pNext)
{
	//find the polynomial order for an element for the next solution pass
	//input:
		//elId = element Id
	//output:
		//pNext = polynomial order
	//return
		//true if polynomial order for the element was increased

	int p;
	SRelement* elem = model.GetElement(elId);
	bool pup = model.FindNextP(elem, p);
	*pNext = p;
	return pup;
}
double FindElementError(int elId)
{
	//find the relative stress error for an element
	//input:
		//elId = element Id
	//return
		//relative stress error

	SRelement* elem = model.GetElement(elId);
	return model.FindElementError(elem);
}
void SetLowStressTolerance(double tol)
{
	//set the relative tolerance for which elements polynomial orders are kept below maxPorderLowStressIn
	//(see setupErrorCheck)
	//input:
		//tol = relativie tolerance, e.g. to keep p below maxPorderLowStressIn for elements with stress less than
			//50% of max in model, tol = 0.5
	model.SetLowStressTolerance(tol);
}

void SetLowStressToleranceFinalAdapt(double tol)
{
	//same as SetLowStressTolerance but for final adaptive solution pass
	model.SetLowStressToleranceFinalAdapt(tol);
}

bool PreProcessPenaltyConstraints()
{
	//pre process constraints for which penalty method will be used
	//(all constraints not in gcs)
	//return:
		//true if any non-gcs constraints were found
	return model.PreProcessPenaltyConstraints();
}

void ProcessConstraints()
{
	//process all constraints for model
	model.ProcessConstraints();
}

int numberEquations()
{
	//number all equations in model.
	//note: 
		//fills internal variable functionEquations in model class
		//ProcessConstraints has to be called 1st. 
		//this assigns a negative number to global dofs that are constrained
		//this routine will only assign an equation number to unconstrained
		//global dofs

	return model.NumberEquations();
}

int GetFunctionEquation(int fun, int dof)
{
	//get the euqation number for a global function and dof
	//return
		//equation number, or -1 for a constrained dof
	return model.GetFunctionEquation(fun, dof);
}

void checkElementMapping()
{
	//check mapping of all elements in model,
	//and perform partial element flattening
	//if any elements with invalid mapping are encountered
	model.checkElementMapping();
}

void allocateConstraints(int ncon)
{
	//allocate space for storage of constraints
	//input:
		//ncon = number of constraints
	model.allocateConstraints(ncon);
}

void inputNodalConstraint(int nid, bool constraineddof[3], double enforcedDisp[3])
{
	//input a nodal constraint
	//input:
		//nid = node Id
		//constraineddof[dof] = true if dof is constrained
		//enforcedDisp[dof] = enforced displacement value if dof is constrained (0 if not enforced)
	SRconstraint* con = model.addConstraint();
	con->SetType(nodalCon);
	con->SetEntityId(nid);
	bool anyNZenfd = false;
	for (int i = 0; i < 3; i++)
	{
		con->SetConstrainedDof(i, constraineddof[i]);
		if (fabs(enforcedDisp[i]) > TINY)
			anyNZenfd = true;
	}
	if (anyNZenfd)
	{
		con->allocateEnforcedDisplacementData(1);
		for (int i = 0; i < 3; i++)
			con->PutEnforcedDisplacementData(0, i, enforcedDisp[i]);
	}
}

int addFaceConstraint(int elemId, int* nidv, bool constraineddof[3], int coordid)
{
	//add a face constraint to the model
	//input:
		//elemId = element Id for element owning the face
		//constraineddof[dof] = true if dof is constrained
		//nidv = vector of face nodes. For a trianlge, nidv[3] should be set to -1
		//coordid = coordinate system Id for the constraint, 0 for gcs
	//return:
		//constraintId for the new constraint
	//note:
		//enforced displacements for the nodes on the face will be input separately in 
		//inputFaceNodeEnfd
	int gno[4];
	SRconstraint* con = model.addConstraint();
	con->SetType(faceCon);
	SRelement* elem = model.GetElement(elemId);
	int fid = model.elemFaceFind(elem, nidv, gno);
	if (fid == -1)
		return -1;
	con->SetEntityId(fid);
	for (int i = 0; i < 3; i++)
		con->SetConstrainedDof(i, constraineddof[i]);
	con->SetCoordId(coordid);
	return model.GetNumConstraints() - 1;
}

void inputFaceNodeEnfd(int constraintId, int localNodeNum, double enforcedDisp[3])
{
	//input enforced displacemnet for a node on a face constraint
	//input:
		//constraintId = constraint Id for the constraint, returned by addFaceConstraint
		//localNodeNum = local node number on the face
		//enforcedDisp[dof] = enforced displacement value if dof is constrained (0 if not enforced)

	SRconstraint* con = model.GetConstraint(constraintId);
	int faceId = con->GetEntityId();
	SRface* face = model.GetFace(faceId);
	if (!con->hasEnforcedDisp())
		con->allocateEnforcedDisplacementData(face->GetNumNodesTotal());
	for (int i = 0; i < 3; i++)
		con->PutEnforcedDisplacementData(localNodeNum, i, enforcedDisp[i]);
}

void calculateRawStress(int elemId, double r, double s, double t, double* stress)
{
	//directly calculate stress in an element from displacements and basis functions
	//input:
		//elemId = element Id
		//r,s,t = natural coordinates in element
	//output:
		//stress = 6 stress components
		//convention: stress[6] = { "xx", "yy", "zz", "xy", "xz", "yz" };

	SRelement* elem = model.GetElement(elemId);
	double strain[6], etx, ety, etz;
	elem->CalculateRawStrain(r, s, t, strain, etx, ety, etz);
	elem->StraintoStress(r, s, t, strain, stress);
}

void calculateSmoothedStress(int elemId, double r, double s, double t, double* stress)
{
	//calculate stress in an element for which smoothed strains had previously been calculated
	//input:
		//elemId = element Id
		//r,s,t = natural coordinates in element
	//output:
		//stress = 6 stress components
		//convention: stress[6] = { "xx", "yy", "zz", "xy", "xz", "yz" };

	SRelement* elem = model.GetElement(elemId);
	double strain[6];
	elem->GetSmoothedStrain(r, s, t, strain);
	elem->StraintoStress(r, s, t, strain, stress);
}

void mapSetup()
{
	//set up quadratic mapping for all elements in model
	model.map.Setup();
}

int ElementBasisDerivs(int elemId, double r, double s, double t, double* dbasisdr, double* dbasisds, double* dbasisdt)
{
	//calculate derivatives of element basis functions with respect to natural coordinates
	//input:
		//elemId = element Id
		//r,s,t = natural coordinates in element
	//output:
		//dbasisdr = derivatives of element basis functions with respect to r
		//dbasisds = derivatives of element basis functions with respect to s
		//dbasisdt = derivatives of element basis functions with respect to t
	//return:
		//number of basis functions

	double* basis = NULL;
	SRelement* elem = model.GetElement(elemId);
	return model.basis.ElementBasisFuncs(r, s, t, elem, basis, dbasisdr, dbasisds, dbasisdt, derivonly);
}

int ElementBasisFuncs(int elemId, double r, double s, double t, double* basis)
{
	//calculate element basis functions
	//input:
		//elemId = element Id
		//r,s,t = natural coordinates in element
	//output:
		//basis = element basis functions
	//return:
		//number of basis functions
	double* dbasisdr = NULL;
	double* dbasisds = NULL;
	double* dbasisdt = NULL;
	SRelement* elem = model.GetElement(elemId);
	return model.basis.ElementBasisFuncs(r, s, t, elem, basis, dbasisdr, dbasisds, dbasisdt, basisonly);
}

//FORTRAN WRAPPERS:

void fsetnumnodes_(int* numNodes)
{
	//allocate space for nodes
	//input:
		//numNodes = total number of nodes in the mesh
	SetNumNodes(*numNodes);
}

int fnumberGlobalFunctions()
{
	//calculate number of global basis functions in model
	//return:
		//number of global basis functions

	return numberGlobalFunctions();
}

int fcreatenode_(int* uid, double* x, double* y, double* z)
{
	//create a node
	//input:
		//uid = node userId
		//x,y,z = node coordinates in global C.S.
	return createNode(*uid, *x, *y, *z);
}

void fsetnumelements_(int* numBricks, int* numWedges, int* numTets)
{
	//allocate space for elements, edges, and faces
	//input:
		//numBricks = number of brick (hexahedral) elements
		//numWedges = number of wedge (pentahedral) elements
		//numTets = number of tetrahedral elements

	SetNumElements(*numBricks, *numWedges, *numTets);
}

extern int fcreateelementiso_(int* uid, int* numnodes, int* nodeids, double* E, double* nu)
{
	//create an isotropic p-element from the definition
	//of a conventional quadratic element
	//input:
		//uid = user Id of the element
		//numnodes = total number (corner and midside) of nodes in the
			//quadratic element. supported values are 10, 15, 20
			//(tets, wedges, and bricks)
		//nodeIds = list of global ids of the corner and midside nodes
		//E = young's modulus
		//nu = Poisson's ratio
	//return:
		//0 for successful completion
		//-1: unsupported number of nodes (linear element)
		//-2: unsupported number of nodes (general)
		//-3: too many unique materials
	//Note: it is assumed that nodeIds are 1-based (1 is the id of the first node in the model)

	for (int i = 0; i < *numnodes; i++)
		nodeBuf[i] = nodeids[i] - 1;

	return createElementIso(*uid, *numnodes, nodeBuf, *E, *nu);
}

void fcreateGlobalFaces_()
{
	//create global faces for model
	createGlobalFaces();
}

int fnumberGlobalFunctions_()
{
	//calculate number of global basis functions in model
	//return:
		//number of global basis functions

	return 0;
}

void fSetupErrorCheck_(int* finalAdapt, double* stressMaxIn, double* ErrorToleranceIn, int* maxPorderIn, int* maxPJumpIn, int* maxPorderLowStressIn)
{
	//set up for error checking
	//input:
		//finalAdapt = true if this is the final adaptivity pass
		//stressMaxIn = max von Mises stress in model
		//ErrorToleranceIn = relative error tolerance, e.g 0.05 for 5 percent
		//maxPorderIn = max polynomial order allowed in model
		//maxPJumpIn = max change in polynomial order allowed this pass
		//maxPorderLowStressIn = max polynomial order allowed in model for elements with low stress

	model.setupErrorCheck(*finalAdapt, *stressMaxIn, *ErrorToleranceIn, *maxPorderIn, *maxPJumpIn, *maxPorderLowStressIn);
}


double* fCalculateStiffnessMatrix_(int* id, int* upperTriangle, double* globalEnfdVec)
{
	//input:
		//id = element id
		//upperTriangle = 1 if output stiffness matrix (symmetric) is to be stored as upper trianlge
						//0 for lower triangle
	//globalEnfdVec = storage for global enforced displacement vector
	//(0 to use internal storage or if enforced displacement not needed)
	//return: element stiffness matrix in format specified by upperTriangle flag
	//note: if globalEnfdVec was provided, it is updated with this Element's contribution
	bool ut = (*upperTriangle == 1);
	int idc = *id - 1;
	return CalculateStiffnessMatrix(idc, ut, globalEnfdVec);
}

int fCheckAutoSacrificialElements_()
{
	//check for sacrificial elements in model (due to point loads, reentrant corners, etc)
	//return:
		//number of sacrificial elements

	if (model.CheckAutoSacrificialElements())
		return 1;
	else
		return 0;
}

int fGetNumFunctions_()
{
	//get number of global functions in model
	//return:
		//number of functions

	return model.GetNumFunctions();
}

void fInitializeErrorCheck_()
{
	//initilialize error checking

	model.initializeErrorCheck();
}

int fCheckForSmallMaxStress_()
{
	//see if max stress in model is small, e.g. free thermal expansion case
	//return:
		//true if max stress is small

	if (model.checkForSmallMaxStress())
		return 1;
	else
		return 0;
}

void fCleanUpErrorCheck_()
{
	//perform memory cleanup tasks for error checker

	model.CleanUpErrorCheck();
}

int fFindNextP_(int* elId, int* pNext)
{
	//find the polynomial order for an element for the next solution pass
	//input:
		//elId = element Id
	//output:
		//pNext = polynomial order
	//return
		//1 if polynomial order for the element was increased else 0

	SRelement* elem = model.GetElement(*elId);
	int p;
	bool ret = model.FindNextP(elem, p);
	*pNext = p;
	if (ret)
		return 1;
	else
		return 0;

}

double fFindElementError_(int* elId)
{
	//find the relative stress error for an element
	//input:
		//elId = element Id
	//return
		//relative stress error


	SRelement* elem = model.GetElement(*elId);
	return model.FindElementError(elem);
}

void fSetLowStressTolerance_(double* tol)
{
	//set the relative tolerance for which elements polynomial orders are kept below maxPorderLowStressIn
	//(see setupErrorCheck)
	//input:
		//tol = relativie tolerance, e.g. to keep p below maxPorderLowStressIn for elements with stress less than
			//50% of max in model, tol = 0.5

	model.SetLowStressTolerance(*tol);
}

void fSetLowStressToleranceFinalAdapt_(double* tol)
{
	//same as fSetLowStressTolerance_ but for final adaptive solution pass

	model.SetLowStressToleranceFinalAdapt(*tol);
}

int fPreProcessPenaltyConstraints_()
{
	//pre process constraints for which penalty method will be used
	//(all constraints not in gcs)
	//return:
		//true if any non-gcs constraints were found

	if (model.PreProcessPenaltyConstraints())
		return 1;
	else
		return 0;
}

void fProcessConstraints_()
{
	//process all constraints for model

	model.ProcessConstraints();
}

int fnumberEquations_()
{
	//number all equations in model.
	//note: 
		//fills internal variable functionEquations in model class
		//fProcessConstraints_ has to be called 1st. 
		//this assigns a negative number to global dofs that are constrained
		//this routine will only assign an equation number to unconstrained
		//global dofs


	return model.NumberEquations();
}

int fGetFunctionEquation_(int* fun, int* dof)
{
	//get the euqation number for a global function and dof
	//return
		//equation number, or -1 for a constrained dof

	return model.GetFunctionEquation(*fun, *dof);
}

void fcheckElementMapping_()
{
	//check mapping of all elements in model,
	//and perform partial element flattening
	//if any elements with invalid mapping are encountered

	model.checkElementMapping();
}

void fallocateConstraints_(int* ncon)
{
	//allocate space for storage of constraints
	//input:
		//ncon = number of constraints

	model.allocateConstraints(*ncon);
}

void finputNodalConstraint_(int* nidf, int constraineddoff[3], double enforcedDisp[3])
{
	//input a nodal constraint
	//input:
		//nid = node Id
		//constraineddof[dof] = true if dof is constrained
		//enforcedDisp[dof] = enforced displacement value if dof is constrained (0 if not enforced)

	bool constraineddof[3];
	for (int i = 0; i < 3; i++)
		constraineddof[i] = (constraineddoff[i] == 1);

	int nid = *nidf - 1;
	inputNodalConstraint(nid, constraineddof, enforcedDisp);
}

int faddFaceConstraint_(int *elemidf, int* nidv, int* constraineddoff, int* coordidf)
{
	//add a face constraint to the model
	//input:
		//elemId = element Id for element owning the face
		//constraineddof[dof] = true if dof is constrained
		//nidv = vector of face nodes. For a trianlge, nidv[3] should be set to -1
		//coordid = coordinate system Id for the constraint, 0 for gcs
	//return:
		//constraintId for the new constraint
	//note:
		//enforced displacements for the nodes on the face will be input separately in 
		//inputFaceNodeEnfd

	int coordid = (*coordidf) - 1;
	int elemId = (*elemidf) - 1;
	bool constraineddof[3];
	for (int i = 0; i < 3; i++)
		constraineddof[i] = (constraineddoff[i] == 1);
	return addFaceConstraint(elemId, nidv, constraineddof, coordid);
}

void finputFaceNodeEnfd_(int* constraintId, int* localNodeNumf, double* enforcedDisp)
{
	//input enforced displacemnet for a node on a face constraint
	//input:
		//constraintId = constraint Id for the constraint, returned by addFaceConstraint
		//localNodeNum = local node number on the face
		//enforcedDisp[dof] = enforced displacement value if dof is constrained (0 if not enforced)

	int localnodenum = (*localNodeNumf) = 1;
	inputFaceNodeEnfd(*constraintId, localnodenum, enforcedDisp);
}

void fcalculateRawStress_(int* elemId, double* r, double* s, double* t, double* stress)
{
	//directly calculate stress in an element from displacements and basis functions
	//input:
		//elemId = element Id
		//r,s,t = natural coordinates in element
	//output:
		//stress = 6 stress components
		//convention: stress(6) = { "xx", "yy", "zz", "xy", "xz", "yz" };


	calculateRawStress(*elemId, *r, *s, *t, stress);
}

void fcalculateSmoothedStress_(int* elemId, double* r, double* s, double* t, double* stress)
{
	//calculate stress in an element for which smoothed strains had previously been calculated
	//input:
		//elemId = element Id
		//r,s,t = natural coordinates in element
	//output:
		//stress = 6 stress components
		//convention: stress(6) = { "xx", "yy", "zz", "xy", "xz", "yz" };


	calculateSmoothedStress(*elemId, *r, *s, *t, stress);
}

void fmapSetup_()
{
	//set up quadratic mapping for all elements in model
	model.map.Setup();
}

int fElementBasisDerivs_(int* elemIdf, double* r, double* s, double* t, double* dbasisdr, double* dbasisds, double* dbasisdt)
{
	//calculate derivatives of element basis functions with respect to natural coordinates
	//input:
		//elemId = element Id
		//r,s,t = natural coordinates in element
	//output:
		//dbasisdr = derivatives of element basis functions with respect to r
		//dbasisds = derivatives of element basis functions with respect to s
		//dbasisdt = derivatives of element basis functions with respect to t
	//return:
		//number of basis functions

	int elemId = *elemIdf - 1;
	return ElementBasisDerivs(elemId, *r, *s, *t, dbasisdr, dbasisdr, dbasisdt);
}

int fElementBasisFuncs_(int* elemIdf, double* r, double* s, double* t, double* basis)
{
	//calculate element basis functions
	//input:
		//elemId = element Id
		//r,s,t = natural coordinates in element
	//output:
		//basis = element basis functions
	//return:
		//number of basis functions

	int elemId = *elemIdf - 1;
	return ElementBasisFuncs(elemId, *r, *s, *t, basis);
}
















