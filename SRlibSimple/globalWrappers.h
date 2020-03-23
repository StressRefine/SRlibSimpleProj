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
// globalWrappers.h: declaration of function wrappers for calling stressRefine library
// from c and fortran
//
//////////////////////////////////////////////////////////////////////

#if !defined(GLOBALWRAPPERS_INCLUDED)
#define GLOBALWRAPPERS_INCLUDED

//C-calleable
void SetNumNodes(int numNodes);
int createNode(int uid, double x, double y, double z);
void SetNumElements(int numBricks, int numWedges, int numTets);
int createElementIso(int uid, int numnodes, int* nodeids, double E, double nu);
void createGlobalFaces();
int numberGlobalFunctions();
void mapSetup();
double* CalculateStiffnessMatrix(int id, bool upperTriangle, double* globalEnfdVec);
int CheckAutoSacrificialElements();
int GetNumFunctions();
void initializeErrorCheck();
bool checkForSmallMaxStress();
void setupErrorCheck(bool finalAdapt, double stressMaxIn, double ErrorToleranceIn, int maxPorderIn, int maxPJumpIn, int maxPorderLowStressIn);
void CleanUpErrorCheck();
bool FindNextP(int elId, int* pNext);
double FindElementError(int elId);
void SetLowStressTolerance(double tol);
void SetLowStressToleranceFinalAdapt(double tol);
bool PreProcessPenaltyConstraints();
void ProcessConstraints();
int numberEquations();
int GetFunctionEquation(int fun, int dof);
void checkElementMapping();
void allocateConstraints(int ncon);
void inputNodalConstraint(int nid, bool constraineddof[3], double enforcedDisp[3]);
int addFaceConstraint(int elemId, int* nidv, bool constraineddof[3], int coordid);
void inputFaceNodeEnfd(int constraintId, int localNodeNum, double enforcedDisp[3]);
void calculateRawStress(int elemId, double r, double s, double t, double* stress);
void calculateSmoothedStress(int elemId, double r, double s, double t, double* stress);
int ElementBasisDerivs(int elemId, double r, double s, double t, double* dbasisdr, double* dbasisds, double* dbasisdt);
int ElementBasisFuncs(int elemId, double r, double s, double t, double* basis);

//Fortran-calleable
void fsetnumnodes_(int* numNodes);
int fcreatenode_(int* uid, double* x, double* y, double* z);
void fsetnumelements_(int* numBricks, int* numWedges, int* numTets);
int fcreateelementiso_(int* uid, int* numnodes, int* nodeids, double* E, double* nu);
void fcreateGlobalFaces_();
int fnumberGlobalFunctions_();
void fmapSetup_();
double* fCalculateStiffnessMatrix_(int* id, int* upperTriangle, double* globalEnfdVec);
int fGetNumFunctions_();
void fInitializeErrorCheck_();
int fCheckForSmallMaxStress_();
void fSetupErrorCheck_(int* finalAdapt, double* stressMaxIn, double* ErrorToleranceIn, int* maxPorderIn, int* maxPJumpIn, int* maxPorderLowStressIn);
void fCleanUpErrorCheck_();
int fFindNextP_(int* elId, int* pNext);
double fFindElementError_(int* elId);
void fSetLowStressTolerance_(double* tol);
void fSetLowStressToleranceFinalAdapt_(double* tol);
int fPreProcessPenaltyConstraints_();
void fProcessConstraints_();
int fnumberEquations_();
int fGetFunctionEquation_(int* fun, int* dof);
void fcheckElementMapping_();
void fallocateConstraints_(int* ncon);
void finputNodalConstraint_(int* nid, int constraineddof[3], double enforcedDisp[3]);
void faddFaceConstraint_(int* elemIdf, int* nidv, int *constraineddoff, int* coordidf);
void finputFaceNodeEnfd_(int* constraintId, int* localNodeNum, double* enforcedDisp);
void fcalculateRawStress_(int* elemId, double* r, double* s, double* t, double* stress);
int fElementBasisDerivs_(int* elemIdf, double* r, double* s, double* t, double* dbasisdr, double* dbasisds, double* dbasisdt);
int fElementBasisFuncs_(int* elemIdf, double* r, double* s, double* t, double* basis);

#endif // !defined(GLOBALWRAPPERS_INCLUDED)
