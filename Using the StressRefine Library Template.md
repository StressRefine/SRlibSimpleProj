Using the StressRefine Library to Add p-adaptivity to a Finite Element Code 
 
Note: whereever functions calls are mentioned below, the functions are defined in SrlibSimple in file globalWrappers.cpp 
There are two steps:
1. Setup so element and stress recovery routines work for arbitrary p-order
2. Wrap an adaptive loop around your solution routine

I will illustrate how to carry out these steps by showing modifications to the pseudo-code of a conventional finite element program in italics.

Notes on programming languages:
The StressRefine library was designed to be called directly from C++
*wrappers have been provided so all the needed functions can be also called from c
*wrappers have been provided so all the needed functions can be also called from fortran using
linux conventions, for example, for a routine cfun(), a wrapper cfun_() is provided, and the function is call from fortran using “call cfun()”. In the descriptions below, I always give the name of the c-function. The corresponding Fortran-calleable functions has a lower case “f” in front. For example, to set the number of nodes in the model, the c function is SetNumNodes, while from fortran it would be “call fSetNumNodes”.

1. Setup for p-adaptive elements
Here is the pseudo-code for finite element setup and linear solution
* Read mesh information. Nodal Coordinates, Element definition (node ids and materials), material properties, loads and constraints
* element bookkeeping- relate local node number to global node id for each element (this is often called the ID array in textbooks)
* process loads
* process constraints- determine which global degrees of freedom are constrained, and determine global equation numbers for the unconstrained dofs, and assemble enforced displacement vector
* calculate element stiffnesses
* Assemble global stiffness and decomp
* backsolve load vector to get global displacement vector
* Evaluate stress: loop over elements
..* download global displacement to local displacements in element
..* call element stress routine
* Postprocess and output in desired format

Here are the modifications needed to make the elements work for arbitrary p-orders, shown in italics

__These steps are only done once to set up the problem:__
* Read mesh information. Nodal Coordinates, Element definition (node ids and materials), material properties, loads and constraints

* _set number of nodes in model_: SetNumNodes(n)
* _input nodal coordinates: for each node_: createNode(int uid, double x, double y, double z)
* _Enter Number of elements of each type (bricks, wedges, tets)_: SetNumElements(int numBricks, int numWedges, int numTets)
..* Note: this must be called even if your are modifying your own elements. stressRefine stores informaton needed for the basis function routines in the elements
* _Use element definitions to create edges and faces for the model_
..* _Note: this must be called even if your are modifying your own elements. stressRefine stores informaton needed for the basis function routines in the elements
..* _For each element_, createElementIso(int uid, int numnodes, int* nodeids, double E, double nu);
..* _after all elements have been created, create global faces_: createGlobalFaces()
* _input nodal and face constraints_: allocateConstraints(int ncon)._ Then for each nodal constraint_: inputNodalConstraint(int nid, bool constraineddof[3], double enforcedDisp[3]);
_for each face constraint_: addFaceConstraint(int elemId, int* nidv, bool constraineddof[3], bool anyenfd);
_This returns the constraintid that was added. If there are enforced displacements,
for each corner node of the face_, inputFaceNodeEnfd(int constraintId, int localNodeNum, double enforcedDisp[3]);
_where_ constraintId _was returned by_ addFaceConstraint

* _If there are any constraints in coordinates systems other than the global coordinate system (gcs), they need to be preprocessed so they can be handled via a penalty method:_ PreProcessPenaltyConstraints

__These steps must be repeated for every solution pass:__
* element bookkeeping- relate local node number to global node id for each element (this is often called the ID array in textbooks)
..* _use the p-orders for the element’s edges to determine number of functions and relate element local function numbers to global functions numbers. routines
* process loads_
..* _No changes are required in processing loads as long as they are nodel, or constant, linear or quadratic distributed surface loads or body loads. They may be processed the same way as for conventional quadratic elements_
* process constraints- _determine which global degrees of freedom are constrained, and determine global equation numbers for the unconstrained dofs, and assemble enforced displacement vector. This is only for constraints in gcs: ProcessConstraints
(non-gcs constraints are handled by the above call to PreProcessPenaltyConstraints and during element stiffness calculation)_
* number equations: _fill up array functionsequations for each function and degree of freedom in the model:_ numberEquations.
_This returns total number of equations. Equation numbers corresponding to each function, dof can then be accessed using_ GetFunctionEquation(fun, dof)_. This returns an equation number if the function, dof is unconstrained, else -1._
* calculate element stiffnesses
..* _either modify your code’s element routines to use p-adaptive displacement functions instead of conventional shape functions (discussed below) OR_
..* _calculate element stiffnesses using stressRefine elements. Since these use stressRefine’s internal quadratic mapping, first call_ mapSetup_, then for each element:_
 CalculateStiffnessMatrix (id, upperTriangle)
..*_NOTES:_
1. _The StressRefine Element Stiffness routine will automatically handle the penalty constraints at the element level. If you are using your own stiffness routines they will have to be modified to handle this._
2. _If using the The StressRefine Element Stiffness routine, it is recommended to call checkElementMapping before entering the loop to calculate elements stiffnesses. This will perform element “partial flattening” link!! if any elements have invalid mapping_
* Assemble global stiffness and decomp
..*_use the modifications above to do the element bookkeeping_
..* _No effect on direct solvers_
..* _For iterative solvers, preconditioners that work with conventional elements may not perform as well with higher order elements._
* backsolve load vector to get global displacement vector
* Evaluate stress: loop over elements
..* download global displacement to local displacements in element
..* call element stress routine
1. _either modify your code’s element stress routines to use p-adaptive displacement functions instead of conventional shape functions (discussed below) OR_
2. _evaluate stresses with stressRefine routine. routine_ void calculateRawStress(int elemId, double r, double s, double t, double* stress)
_Note: hierarchical basis functions cannot be evaluated directly at corners, they are singular. Either evaluate them away from the corner and project OR_
_evaluate them at optimal sampling points and smooth to corners and midnodes (recommended, more accurate): see _ post.globalstrainsmooth  _in SRwithMklProj_
_This calculates the smoothed strain vector for each function of each element. You will need to modify _globalstrainsmooth_ if you are using a different solver than Intel Mkl Pardiso._
_Then to calculate the stress from the smoothed strains at any point in an element:_
void calculateSmoothedStress(int elemId, double r, double s, double t, double* stress);
Postprocess and output in desired format
2. **Adaptive Solution loop**
* Initialization: perform the steps from “read mesh information…” to “PreProcessPenaltyConstraints”
**LOOP:**
* Perform the solution steps, starting with “element bookkeeping”, down to “Evaluate stress: loop over elements”
* Calculate errors from solution
..* setup: setupErrorCheck(bool finalAdapt, double stressMax, double ErrorTolerance, int maxPorder, int maxPJump, int maxPorderLowStress); This routine calculates the tractions at sample points for each local face of each element
..* for each element: FindElementError(int elId); This returns the element error relative to the max stress in model
* Calcuate new require p-orders: _for each element:_ FindNextP(int elId, int* pNext) _returns true if the element needed to be adapted else false._ Note: if FindNextP returns false, the element stiffness matrix for the element with elId does not need to be recomputed in the next pass. The stiffness matrix from the previous solution pass can be reused.
* _if no element was adapted, break loop_
**UNTIL CONVERGED**

**Modifying Element Stiffness Routines for P-adaptivity**
(a similar modification is done for stress Routines)
_Wherever your element routine calculates derivatives of the basis functions for displacement, this needs to be replaced with the routines from the stressRefine library:_ int ElementBasisDerivs(int elemId, double r, double s, double t, double* dbasisdr, double* dbasisds, double* dbasisdt)

