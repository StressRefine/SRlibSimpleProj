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
// SRelement.cpp: implementation of the SRelement class.
//
//////////////////////////////////////////////////////////////////////

#include "SRmodel.h"

extern SRmodel model;

SRBTC::SRBTC()
{
	BTC11 = BTC12 = BTC13 = BTC14 = BTC15 = BTC16 = 0.0;
	BTC21 = BTC22 = BTC23 = BTC24 = BTC25 =  BTC26 = 0.0;
	BTC21 = BTC32 = BTC33 = BTC34 = BTC35 = BTC36 = 0.0;
}


SRface * SRlocalFace::GetFace()
{
	//get pointer to global face corresponding to this local face of an element
	//return:
		//pointer to global face

	return model.GetFace(globalFaceId);
}


SRelement::SRelement()
{
	type = tet;
	sacrificial = false;
	svmMax = 0.0;
	flattened = false;
	flattenedHighStress = false;
	stiffLength = 0;
	constrained = false;
	elVol = 0.0;
	flattenFraction = 0.0;
	approxVol = 0.0;
};

void SRelement::Create(int useridt, int nnodes, int nodest[], SRmaterial* mat)
{
	if (nnodes == 10)
	{
		Create(tet, useridt, nnodes, nodest, mat);
	}
	else if (nnodes == 15)
	{
		Create(wedge, useridt, nnodes, nodest, mat);
	}
	else if (nnodes == 20)
	{
		Create(brick, useridt, nnodes, nodest, mat);
	}

}
void SRelement::Create(SRelementType typet, int useridt, int nnodes, int nodest[], SRmaterial* mat)
{
	//Create element. store nodes, material, edges, faces in element arrays
	//input:
		//typet = element type
		//useridt = user ID
		//nodest = vector of node ids
		//mat = material property

	type = typet;
	userId = useridt;
	elMat = mat;
	int ncorner, nej;
	penaltyStiffness = 0.0;
	if (type == tet)
	{
		ncorner = 4;
		nej = 6;
		localEdges.Allocate(nej);
		localFaces.Allocate(4);

	}
	else if (type == wedge)
	{
		ncorner = 6;
		nej = 9;
		localEdges.Allocate(nej);
		localFaces.Allocate(5);
	}
	else
	{
		ncorner = 8;
		nej = 12;
		localEdges.Allocate(nej);
		localFaces.Allocate(6);
	}

	nodeIds.Allocate(ncorner);
	for (int i = 0; i < ncorner; i++)
		nodeIds.Put(i, nodest[i]);
	pChanged = true;
	error = 0.0;
}

void SRelement::GetStress(double r, double s, double t, double stress[])
{
	//look up smoothed stress in an element at a point
	//input:
		//r,s,t = natural coordinates
	//output:
		//stress[6] = vector of strains at this location

	double strain[6];
	GetSmoothedStrain(r, s, t, strain);
	StraintoStress(r, s, t, strain, stress);
}

void SRelement::PutRawStrain(int i, double *strain)
{
	//store raw strains at a gauss point of this elements
	//input:
		//i = gauss point number
		//strain = strain tensor
	//note:
		//stores strain in element rawStrains matrix

	for (int e = 0; e < 6; e++)
		rawStrains.Put(i, e, strain[e]);
}


void SRelement::GetEdgeLocalNodeIds(int lej, int& ln0, int& ln1)
{
	//find the element local node ids of the two ends of a local edge
	//input:
		//lej = local edge number
	//output:
		//ln0, ln1 = element local node ids of the ends of the local edge

	if (type == tet)
	{
		ln0 = model.map.GetTetEdgeLocalNode(lej, 0);
		ln1 = model.map.GetTetEdgeLocalNode(lej, 1);
	}
	else if (type == wedge)
	{
		ln0 = model.map.GetWedgeEdgeLocalNode(lej, 0);
		ln1 = model.map.GetWedgeEdgeLocalNode(lej, 1);
	}
	else
	{
		ln0 = model.map.GetBrickEdgeLocalNode(lej, 0);
		ln1 = model.map.GetBrickEdgeLocalNode(lej, 1);
	}
}


void SRelement::GetEdgeNodeIds(int lej, int &n1, int &n2)
{
	//get the node ids at the ends of a local edge of an element.
	//this is for use before global edge assignment, so this cannot be looked
	//up with localEdges.
	//input:
		//lej = local edge number
	//output:
		//n1,n2 = nodes at ends of edge

	int ln0, ln1;
	GetEdgeLocalNodeIds(lej, ln0, ln1);
	n1 = nodeIds.Get(ln0);
	n2 = nodeIds.Get(ln1);
}

int SRelement::GetLocalEdgeNodeId(int lej, int lnode)
{ 
	//look up the global node id corresponding to a local node of a local edge

	return localEdges.Get(lej).GetEdge()->GetNodeId(lnode);
}

void SRelement::GetFaceNodes(int lface, int &n1, int &n2, int &n3, int &n4)
{
	//get the node numbers at the corners of a local face of an element.
	//this is for use before global face assignment, so this cannot be looked
	//up with localFaces.
	//input:
		//lface = local face number
	//output:
		//n1,n2,n3,n4 = nodes at corners of face; n4 = -1 for tri face

	if (type == tet)
	{
		if (lface == 0)
		{
			//local face 1 = 2-3-4 (face for which L1=0)
			n1 = nodeIds.Get(1);
			n2 = nodeIds.Get(2);
			n3 = nodeIds.Get(3);
		}
		else if (lface == 1)
		{
			//local face 2 = 1-3-4 (face for which L2=0)
			n1 = nodeIds.Get(0);
			n2 = nodeIds.Get(2);
			n3 = nodeIds.Get(3);
		}
		else if (lface == 2)
		{
			//local face 3 = 1-2-4 (face for which L3=0)
			n1 = nodeIds.Get(0);
			n2 = nodeIds.Get(1);
			n3 = nodeIds.Get(3);
		}
		else if (lface == 3)
		{
			//local face 4 = 1-2-3 (face for which L4=0)
			n1 = nodeIds.Get(0);
			n2 = nodeIds.Get(1);
			n3 = nodeIds.Get(2);
		}
		n4 = -1;
	}
	else if (type == wedge)
	{
		GetWedgeFaceNodes(lface, n1, n2, n3, n4);
	}
	else
	{
		GetBrickFaceNodes(lface, n1, n2, n3, n4);
	}
}

SRedge* SRelement::GetEdge(int localedgenum)
{
	//look up global edge corresponding to element local edge
	//input:
		//localedgenum = element local edge number
	//return
		//pointer to global edge

	int gej = localEdges.GetPointer(localedgenum)->GetGlobalEdgeId();
	return model.GetEdge(gej);
}

SRface* SRelement::GetFace(int localfacenum)
{
	//look up global face corresponding to element local face
	//input:
		//localfacenum = element local face number
	//return
		//pointer to global face

	return localFaces.GetPointer(localfacenum)->GetFace();
}


SRnode* SRelement::GetNode(int localnodenum)
{
	//get the pointer to the global node corresponding to local node number
	//input:
		//localnodenum = local number in element
	//return:
		//pointer to the global node

	int nodeid = nodeIds.Get(localnodenum);
	return model.GetNode(nodeid);
}

void SRelement::colFunLoopIsoOrtho(double* dbdxv, double* dbdyv, double* dbdzv, double w, int rowfun, SRBTC& btc, double* stiff)
{
	//implement code inside "colfun" loop for calculating stiffness matrix
	//isotropic/orthotropic material version
	//input:
		// dbdxv, dbdyv, dbdzv: vectors containing derivatives of all basis functions w.r.t. x,y,z at current gauss point
		// w = integration weight at current gauss point
		// rowfun = current row function number
	//modified:
		// stiff = element stiffness matrix
	//scratch:
		//btc = "BT*C" = transpose of strain-displacement times elasticity matrix for current row function
	//note:
		//this uses "vanilla" colfun loops.
		//I tried using cblas_daxpy instead but it is slower- have to use inc 3 for elstiff updates 
		//separate tuning tests showed daxpy for inc 3 becomes faster for vector length > 300
		//for p8 brick, nfun = 192 = max vector length, avg length is 96. so daxpy is slower

	double dbdxi = dbdxv[rowfun];
	double dbdyi = dbdyv[rowfun];
	double dbdzi = dbdzv[rowfun];
	FillBTC(w, dbdxi, dbdyi, dbdzi, btc);

	int nfun = GetNumFunctions();
	int colfun;
	int rowcolloc;
	double a, b, c;

	//rowdof 0:
	int roweq = rowfun * 3;
	int rowcolloc0 = stiffDiag.Get(roweq);
	//coldof 0
	rowcolloc = rowcolloc0;
	a = btc.BTC11;
	b = btc.BTC14;
	c = btc.BTC15;
	for (colfun = rowfun; colfun < nfun; colfun++)
	{
		stiff[rowcolloc] += (a * dbdxv[colfun] + b * dbdyv[colfun] + c * dbdzv[colfun]);
		rowcolloc += 3;
	}
	//coldof 1
	rowcolloc = rowcolloc0 + 1;
	a = btc.BTC12;
	b = btc.BTC14;
	for (colfun = rowfun; colfun < nfun; colfun++)
	{
		stiff[rowcolloc] += (a * dbdyv[colfun] + b * dbdxv[colfun]);
		rowcolloc += 3;
	}
	//coldof 2
	rowcolloc = rowcolloc0 + 2;
	a = btc.BTC13;
	b = btc.BTC15;
	for (colfun = rowfun; colfun < nfun; colfun++)
	{
		stiff[rowcolloc] += (a * dbdzv[colfun] + b * dbdxv[colfun]);
		rowcolloc += 3;
	}

	//rowdof 1:
	//note: 1st coldof is skipped for colfun = rowfun, rowdof 1 because it is below diagonal
	roweq++;
	rowcolloc0 = stiffDiag.Get(roweq);
	//coldof 0
	//stiffdiag points to coldof 1 for this fun, so skip 2 to get to coldof 0 for rowfun + 1:
	rowcolloc = rowcolloc0 + 2;
	a = btc.BTC21;
	b = btc.BTC24;
	for (colfun = rowfun + 1; colfun < nfun; colfun++)
	{
		stiff[rowcolloc] += (a * dbdxv[colfun] + b * dbdyv[colfun]);
		rowcolloc += 3;
	}
	//coldof 1
	rowcolloc = rowcolloc0;
	a = btc.BTC22;
	b = btc.BTC24;
	c = btc.BTC26;
	for (colfun = rowfun; colfun < nfun; colfun++)
	{
		stiff[rowcolloc] += (a * dbdyv[colfun] + b * dbdxv[colfun] + c * dbdzv[colfun]);
		rowcolloc += 3;
	}
	//coldof 2
	rowcolloc = rowcolloc0 + 1;
	a = btc.BTC23;
	b = btc.BTC26;
	for (colfun = rowfun; colfun < nfun; colfun++)
	{
		stiff[rowcolloc] += (a * dbdzv[colfun] + b * dbdyv[colfun]);
		rowcolloc += 3;
	}

	//rowdof 2:
	//note: 1st and 2nd coldofs are skipped for colfun = rowfun, rowdof 2 because they are below diagonal
	roweq++;
	rowcolloc0 = stiffDiag.Get(roweq);
	//coldof 0
	//stiffdiag points to coldof 2 for this fun, so skip 1 to get to coldof 0 for rowfun + 1:
	rowcolloc = rowcolloc0 + 1;
	a = btc.BTC31;
	b = btc.BTC35;
	for (colfun = rowfun + 1; colfun < nfun; colfun++)
	{
		stiff[rowcolloc] += (a * dbdxv[colfun] + b * dbdzv[colfun]);
		rowcolloc += 3;
	}
	//coldof 1
	//stiffdiag points to coldof 2 for this fun, so skip 2 to get to coldof 1 for rowfun + 1:
	rowcolloc = rowcolloc0 + 2;
	a = btc.BTC32;
	b = btc.BTC36;
	for (colfun = rowfun + 1; colfun < nfun; colfun++)
	{
		stiff[rowcolloc] += (a * dbdyv[colfun] + b * dbdzv[colfun]);
		rowcolloc += 3;
	}
	//coldof 2
	rowcolloc = rowcolloc0;
	a = btc.BTC33;
	b = btc.BTC35;
	c = btc.BTC36;
	for (colfun = rowfun; colfun < nfun; colfun++)
	{
		stiff[rowcolloc] += (a * dbdzv[colfun] + b * dbdxv[colfun] + c * dbdyv[colfun]);
		rowcolloc += 3;
	}
}

void SRelement::colFunLoopGenAniso(double* dbdxv, double* dbdyv, double* dbdzv, double w, int rowfun, SRBTC& btc, double* stiff)
{
	//implement code inside "colfun" loop for calculating stiffness matrix
	//isotropic material version
	//input:
		// dbdxv, dbdyv, dbdzv: vectors containing derivatives of all basis functions w.r.t. x,y,z at current gauss point
		// w = integration weight at current gauss point
		// rowfun = current row function number
	//modified:
		// stiff = element stiffness matrix
	//scratch:
		//btc = "BT*C" = transpose of strain-displacement times elasticity matrix for current row function
	double kel33[3][3];

	double dbdxi = dbdxv[rowfun];
	double dbdyi = dbdyv[rowfun];
	double dbdzi = dbdzv[rowfun];
	FillBTC(w, dbdxi, dbdyi, dbdzi, btc);
	int nfun = GetNumFunctions();
	for (int colfun = rowfun; colfun < nfun; colfun++)
	{
		double dbdxj = dbdxv[colfun];
		double dbdyj = dbdyv[colfun];
		double dbdzj = dbdzv[colfun];
		FillKel33GenAnisoRowWise(btc, rowfun, colfun, dbdxj, dbdyj, dbdzj, kel33);

		int roweq = rowfun * 3;
		for (int rowdof = 0; rowdof < 3; rowdof++)
		{
			int coldof0 = 0;
			//colfun = rowfun, special case because LT of Kel33 not stored:
			if (colfun == rowfun)
				coldof0 = rowdof;
			int rowcolloc = stiffDiag.Get(roweq);
			for (int coldof = coldof0; coldof < 3; coldof++)
			{
				int coleq = colfun * 3 + coldof;
				rowcolloc = GetStiffnessLocation(roweq, coleq);
				stiff[rowcolloc] += kel33[rowdof][coldof];
			}
			roweq++;
		}
	}
}

double* SRelement::CalculateStiffnessMatrix(int& numEq, double* globalEnfdVec)
{
	//Calculate Element Stiffness Matrix for a three dimensional element
	//input:
		//globalEnfdVec = storage for global enforced displacement vector
		//(NULL to use internal storage or if enforced displacement not needed)
	//output:
		//numEq = number of equations in element stiffmatrix
	//return:
		//Element Stiffness Matrix stored symmetrically (upper triangular)
	double* stiff = NULL;

	int i, neq;
	int nfun = globalFunctionNumbers.GetNum();
	int nint, gp;
	SRBTC btc;
	double w, dbdx, dbdy, dbdz;
	nint = model.math.FillGaussPoints(this);

	numEq = neq = 3 * nfun;
	FillStiffDiag(neq);
	int len = neq * (neq + 1) / 2;

	if (len > stiffnessMatrix.GetNum())
		stiffnessMatrix.Allocate(len);
	//don't need to recalculate stiffness if p not changed:
	stiff = stiffnessMatrix.GetVector();
	if (!pChanged)
		return stiff;

	for (i = 0; i < len; i++)
		stiff[i] = 0.0;

	SRElementData* eldata = model.GetElementData();

	SRdoubleMatrix& dbdxm = eldata->Getdbdx();
	SRdoubleMatrix& dbdym = eldata->Getdbdy();
	SRdoubleMatrix& dbdzm = eldata->Getdbdz();
	double* wv = eldata->GetIntwt();
	basisVec = eldata->GetBasisVec();
	dbasisdr = eldata->Getdbasisdr();
	dbasisds = eldata->Getdbasisds();
	dbasisdt = eldata->Getdbasisdt();

	elVol = 0.0;
	for (gp = 0; gp < nint; gp++)
	{
		double r, s, t;
		model.math.GetGP3d(gp, r, s, t, w);
		double detj = FillMapping(r, s, t);
		if (detj < approxVol * SMALL)
		{

			REPPRINT(" bad Jacobian. Element: %d", GetUserid());
			ERROREXIT;
		}
		w *= detj;
		elVol += w;
		wv[gp] = w;
		FillBasisFuncs(r, s, t, derivonly);

		for (int fun = 0; fun < nfun; fun++)
		{
			double dbdr = dbasisdr[fun];
			double dbds = dbasisds[fun];
			double dbdt = dbasisdt[fun];
			XyzDerivatives(dbdr, dbds, dbdt, dbdx, dbdy, dbdz);
			dbdxm.Put(gp, fun, dbdx);
			dbdym.Put(gp, fun, dbdy);
			dbdzm.Put(gp, fun, dbdz);
		}
	}

	bool isGenAniso = (elMat->GetType() == genAniso) || model.UseSimpleElements();

	for (gp = 0; gp < nint; gp++)
	{
		double* dbdxv = dbdxm.GetRow(gp);
		double* dbdyv = dbdym.GetRow(gp);
		double* dbdzv = dbdzm.GetRow(gp);
		w = wv[gp];
		if (!isGenAniso)
		{
			for (int rowfun = 0; rowfun < nfun; rowfun++)
				colFunLoopIsoOrtho(dbdxv, dbdyv, dbdzv, w, rowfun, btc, stiff);
		}
		else
		{
			for (int rowfun = 0; rowfun < nfun; rowfun++)
				colFunLoopGenAniso(dbdxv, dbdyv, dbdzv, w, rowfun, btc, stiff);
		}
	}

	if (hasLcsConstraint())
		AddPenaltyToStiffnessandEnfDisp(stiff, globalEnfdVec);


	stiffLength = len;

	return stiff;

}

void SRelement::AddPenaltyToStiffnessandEnfDisp(double *stiff, double* globalEnfdVec)
{
	//handle constraints in non-gcs coordinates with penalty method
	//input:
		//globalEnfdVec = storage for global enforced displacement vector
		//(0 to use internal storage or if enforced displacement not needed)
	//output:
		//stiff updated with contribution of penalty
	//notes:
		//if there are any enforced displacements, global force vector is modified
		//this routine if called if there is a non-gcs edge or face constraint on one of the boundary faces of this element

	int faceToElFuns[MAXFACEBASISFUNCTIONS];
	int condof, rowfun, colfun, rowdof, coldof, roweq, coleq, elrowfun, elcolfun;
	SRface *face;
	SRconstraint *con;
	double rf, sf, w, bi, bj, enfdisp = 0.0, detj, bibj, bibjelB;
	int i, nint, nfun, symloc;
	double *basisv = basisVec;
	SRvec3 elocal;
	double elA, elB;
	SRvec3 pos;


	double* forceVec = NULL;
	if (globalEnfdVec != 0)
		forceVec = globalEnfdVec;
	else
		model.GetElementData()->GetEnfdForceVec();

	for (int facon = 0; facon < faceLCSConstraints.GetNum(); facon++)
	{
		int lface = faceLCSConstraints.Get(facon);
		face = localFaces.GetPointer(lface)->GetFace();
		con = face->GetConstraint();
		if (con == NULL)
			ERROREXIT;
		SRcoord* coord = con->GetCoord();
		//map face basis function numbers to element:
		nfun = MapFaceToElFuns(lface, faceToElFuns);
		nint = model.math.FillGaussPoints(face, -1);// -1 means p not being provided, use face pmax
		int gfun, geq;
		for (condof = 0; condof < 3; condof++)
		{
			if (!con->IsConstrainedDof(condof))
				continue;
			for (i = 0; i < nint; i++)
			{
				model.math.GetGP2d(i, rf, sf, w);
				model.basis.FaceBasisFuncs(rf, sf, face, basisv);
				detj = face->Jacobian(rf, sf);
				w *= detj;
				double pw = penaltyStiffness *w;
				face->Position(rf, sf, pos);
				coord->CalculateLocalDirection(pos, condof, elocal);
				if (con->hasEnforcedDisp())
					enfdisp = con->GetFaceEnforcedDisp(rf, sf, condof);
				for (rowfun = 0; rowfun< nfun; rowfun++)
				{
					gfun = face->GetGlobalFunctionNumber(rowfun);
					bj = basisv[rowfun];
					bj *= pw;
					elrowfun = faceToElFuns[rowfun];
					//contribution to enforced displacement:
					if (con->hasEnforcedDisp())
					{
						for (rowdof = 0; rowdof < 3; rowdof++)
						{
							elB = elocal.d[rowdof];
							geq = model.GetFunctionEquation(gfun, rowdof);
							if (geq >= 0)
								forceVec[geq] += (enfdisp*bj*elB);
						}
					}
					for (colfun = 0; colfun < nfun; colfun++)
					{
						bi = basisv[colfun];
						elcolfun = faceToElFuns[colfun];
						bibj = bj*bi;
						for (rowdof = 0; rowdof < 3; rowdof++)
						{
							roweq = elrowfun * 3 + rowdof;
							elB = elocal.d[rowdof];
							bibjelB = bibj*elB;
							for (coldof = 0; coldof<3; coldof++)
							{
								elA = elocal.d[coldof];
								coleq = elcolfun * 3 + coldof;
								if (roweq > coleq)
									continue;
								symloc = GetStiffnessLocationAboveDiag(roweq, coleq);
								stiff[symloc] += (bibjelB*elA);
							}
						}
					}
				}
			}
		}
	}
	for (int ncon = 0; ncon < nodeLcSConstraints.GetNum(); ncon++)
	{
		int lnode = nodeLcSConstraints.Get(ncon);
		SRnode* node = getNodeOrMidNode(lnode);
		con = node->GetConstraint();
		SRcoord* coord = con->GetCoord();
		int elfun = lnode;
		if (node->isMidSide())
		{
			int eid = node->GetMidSideEdgeOwner();
			int lej = localEdgeMatch(eid);
			elfun = MapEdgeP2FnToElFun(lej);
		}
		for (condof = 0; condof < 3; condof++)
		{
			if (!con->IsConstrainedDof(condof))
				continue;
			coord->CalculateLocalDirection(node->Position(), condof, elocal);
			if (con->hasEnforcedDisp())
				enfdisp = con->getDisp(0, condof);
			int gfun;
			if (!node->isMidSide())
				gfun = node->GetGlobalFunctionNumber();
			else
			{
				int eid = node->GetMidSideEdgeOwner();
				SRedge* edge = model.GetEdge(eid);
				gfun = edge->GetGlobalFunctionNumber(2);
			}
			for (rowdof = 0; rowdof < 3; rowdof++)
			{
				elB = elocal.d[rowdof] * penaltyStiffness;
				//contribution to enforced displacement:
				if (con->hasEnforcedDisp())
				{
					int geq = model.GetFunctionEquation(gfun, rowdof);
					if (geq >= 0)
						forceVec[geq] += (enfdisp*elB);
				}

				roweq = elfun * 3 + rowdof;
				for (coldof = 0; coldof < 3; coldof++)
				{
					elA = elocal.d[coldof];
					coleq = elfun * 3 + coldof;
					if (roweq > coleq)
						continue;
					symloc = GetStiffnessLocationAboveDiag(roweq, coleq);
					stiff[symloc] += (elA*elB);
				}
			}
		}
	}
}

void SRelement::CalibratePenalty()
{
	//calibrate penalty stiffness for constraints in curvilinear coordinates
	//note:
		//stores class variable penaltyStiffness = stiffness of a thin layer of same material as element
		//duplicate calls ignored

	if (penaltyStiffness > TINY)
		return;
	penaltyStiffness = 0.0;
	double E = elMat->MatScale();

	//penalty multiplier. this makes penalty stiffness equivalent to
	//that of a layer, same material as element,
	//thickness = (element size)/penaltyMult
	double penaltyMult = 1.0e4;

	penaltyStiffness = (penaltyMult*E / size);
}

void SRelement::FillKel33GenAnisoRowWise(SRBTC& btc, int rowfun, int colfun, double dbdx, double dbdy, double dbdz, double kel33[3][3])
{
	//Fill 3x3 portion of elemental stiffness matrix for a general Anisotropic material,
	//compatible with rowwise storage of elements
	//input:
		//btc = strain-displacement matrixe times elasticity matrix
		//rowfun = row function number
		//colfun = column function number
		//dbdx, dbdy, dbdz = derivatives of basis function
	//output:
		//fills kel33[3][3]

	kel33[0][0] = btc.BTC11*dbdx + btc.BTC14*dbdy + btc.BTC15*dbdz;
	kel33[0][1] = btc.BTC12*dbdy + btc.BTC14*dbdx + btc.BTC16*dbdz;
	kel33[0][2] = btc.BTC13*dbdz + btc.BTC15*dbdx + btc.BTC16*dbdy;
	kel33[1][1] = btc.BTC22*dbdy + btc.BTC24*dbdx + btc.BTC26*dbdz;
	kel33[1][2] = btc.BTC23*dbdz + btc.BTC25*dbdx + btc.BTC26*dbdy;
	kel33[2][2] = btc.BTC33*dbdz + btc.BTC35*dbdx + btc.BTC36*dbdy;
	if (rowfun != colfun)
	{
		kel33[1][0] = btc.BTC21*dbdx + btc.BTC24*dbdy + btc.BTC25*dbdz;
		kel33[2][0] = btc.BTC31*dbdx + btc.BTC34*dbdy + btc.BTC35*dbdz;
		kel33[2][1] = btc.BTC32*dbdy + btc.BTC34*dbdx + btc.BTC36*dbdz;
	}
}

void SRelement::StraintoStress(double r, double s, double t, double strain[], double stress[])
{
	//calculate stress given strain at a point of this elements
	//input:
		//r,s,t = natural coordinate at this point
		//strain = strain tensor
	//output:
		//stress = stress tensor

	SRthermalForce *tf = model.GetThermalForce();
	double mechStrain[6];
	int i;
	for(i = 0; i < 6; i++)
		mechStrain[i] = strain[i];
	if(tf != NULL)
	{
		double deltaTemp, etx, ety, etz;
		deltaTemp = tf->GetTemp(this, r, s, t) - elMat->getTref();
		etx = deltaTemp*elMat->getAlphax();
		ety = deltaTemp*elMat->getAlphay();
		etz = deltaTemp*elMat->getAlphaz();
		mechStrain[0] -= etx;
		mechStrain[1] -= ety;
		mechStrain[2] -= etz;
	}
	if (elMat->GetType() == iso)
	{
		double c11 = elMat->getC11();
		double lambda = elMat->getLambda();
		double G = elMat->getG();
		stress[0] = c11*mechStrain[0] + lambda*(mechStrain[1] + mechStrain[2]);
		stress[1] = c11*mechStrain[1] + lambda*(mechStrain[0] + mechStrain[2]);
		stress[2] = c11*mechStrain[2] + lambda*(mechStrain[0] + mechStrain[1]);
		stress[3] = G*mechStrain[3];
		stress[4] = G*mechStrain[4];
		stress[5] = G*mechStrain[5];
	}
	else if (elMat->GetType() == ortho)
	{
		SRcij& cij = elMat->getOrthoCij();
		stress[0] = cij.c11*mechStrain[0] + cij.c12*mechStrain[1] + cij.c13*mechStrain[2];
		stress[1] = cij.c12*mechStrain[0] + cij.c22*mechStrain[1] + cij.c23*mechStrain[2];
		stress[2] = cij.c13*mechStrain[0] + cij.c23*mechStrain[1] + cij.c33*mechStrain[2];
		stress[3] = cij.c44*mechStrain[3];
		stress[4] = cij.c55*mechStrain[4];
		stress[5] = cij.c66*mechStrain[5];
	}
	else if (elMat->GetType() == genAniso)
	{
		SRgenAnisoCij& gcij = elMat->getGenAnisoCij();
		for(int i = 0; i < 6; i++)
		{
			stress[i] = 0.0;
			for(int j = 0; j < 6; j++)
				stress[i] += (gcij.getC(i,j) * mechStrain[j]);
		}
	}

	for (int i = 0; i < 6; i++)
	{
		double ae = fabs(stress[i]);
		if (ae > maxStress)
			maxStress = ae;
	}

}

void SRelement::NodeNaturalCoords(int lnode, double &r, double &s, double &t)
{
	//look up natural coordinate position of a local node of an element (corner or midside node)
	//input:
		//lnode = local node number
	//output:
		//r,s,t = natural coordinates in element

	if (lnode >= nodeIds.GetNum())
	{
		//midside node
		int lej = lnode - nodeIds.GetNum();
		model.map.ElementNaturalCoordsAtMidedge(this, lej, r, s, t);
		return;
	}

	if (type == tet)
	{
		if (lnode == 0)
		{
			r = -1.0;
			s = 0.0;
			t = 0.0;
		}
		else if (lnode == 1)
		{
			r = 1.0;
			s = 0.0;
			t = 0.0;
		}
		else if (lnode == 2)
		{
			r = 0.0;
			s = SQRT3;
			t = 0.0;
		}
		else if (lnode == 3)
		{
			r = 0.0;
			s = SQRT3OVER3;
			t = TWOSQRTTWOTHIRDS;
		}
		else
			ERROREXIT;
	}
	else if (type == wedge)
		GetWedgeNodePos(lnode, r, s, t);
	else
		GetBrickNodePos(lnode, r, s, t);
}

void SRelement::Clear()
{
	//free memory of an element that is no longer needed
	localEdges.Free();
	localFaces.Free();
	globalFunctionNumbers.Free();
}

void SRelement::DownloadDisplacement()
{
	//download from global solution to local displacements of this element
	//note:
		//fills up class variable dispEl

	int elfun;
	int dof, gfun, leq;
	int nfun = globalFunctionNumbers.GetNum();
	dispEl.Allocate(3 * nfun);
	double *dispElVec = dispEl.GetVector();
	for (elfun = 0; elfun < nfun; elfun++)
	{
		gfun = globalFunctionNumbers.Get(elfun);
		for (dof = 0; dof < 3; dof++)
		{
			leq = 3 * elfun + dof;
			double disp = model.GetDisplacementCoeff(gfun, dof);
			dispElVec[leq] = disp;
		}
	}
}

int SRelement::GetMaxPorder()
{
	//get max p order of this element

	int i, p;
	maxPorder = 0;
	for(i = 0; i < localEdges.GetNum(); i++)
	{
		p = localEdges.Get(i).GetPOrder();
		if(p > maxPorder)
			maxPorder = p;
	}
	return maxPorder;
}

double SRelement::CalculateRawStrain(double r, double s, double t, double strain[], double& etx, double& ety, double& etz)
{
	//calculate raw strain in an element at a point
	//input:
		//r,s,t = natural coordinates in element
		//dbasisdr, dbasisds, dbasisdt = element basis derivative vectors if already calculated, else NULL
	//output:
		//strain[6] = vector of strains at this location
		//etx, ety, etz = thermal strain at this location
	//return:
		//maximum thermal strain
	//Note: if thermal loading is present, strain is total strain

	int fun, nfun, eq;
	double dbdx, dbdy, dbdz, dbdr, dbds, dbdt, uel, vel, wel, ex, ey, ez, gamxy, gamxz, gamyz;

	nfun = globalFunctionNumbers.GetNum();

	ex = 0.0;
	ey = 0.0;
	ez = 0.0;
	gamxy = 0.0;
	gamxz = 0.0;
	gamyz = 0.0;
	double *dispElVec = dispEl.GetVector();
	for (fun = 0; fun < nfun; fun++)
	{
		dbdr = dbasisdr[fun];
		dbds = dbasisds[fun];
		dbdt = dbasisdt[fun];
		XyzDerivatives(dbdr, dbds, dbdt, dbdx, dbdy, dbdz);
		eq = 3 * fun;
		uel = dispElVec[eq];
		vel = dispElVec[eq + 1];
		wel = dispElVec[eq + 2];
		ex += dbdx*uel;
		ey += dbdy*vel;
		ez += dbdz*wel;
		gamxy += (dbdy*uel + dbdx*vel);
		gamxz += (dbdz*uel + dbdx*wel);
		gamyz += (dbdz*vel + dbdy*wel);
	}

	strain[0] = ex;
	strain[1] = ey;
	strain[2] = ez;
	strain[3] = gamxy;
	strain[4] = gamxz;
	strain[5] = gamyz;

	//max thermal strain component at this point:
	double eT = 0.0;
	etx = ety = etz = 0.0;
	SRthermalForce* tf = model.GetThermalForce();
	if(tf != NULL)
	{
		double deltaTemp, etx, ety, etz;
		deltaTemp = tf->GetTemp(this, r, s, t) - elMat->getTref();
		etx = deltaTemp*elMat->getAlphax();
		ety = deltaTemp*elMat->getAlphay();
		etz = deltaTemp*elMat->getAlphaz();
		eT = fabs(etx);
		ety = fabs(ety);
		etz = fabs(etz);
		if(ety > eT)
			eT = ety;
		if(etz > eT)
			eT = etz;
	}

	return eT;
}

int SRelement::MapFaceToElFuns(int lface, int FaceToElFuns[])
{
	//get basis function numbers on a face to corresponding local element basis function numbers
	//for this element that owns the face
	//input:
		//lface = local face number
	//output:
		//FaceToElFuns[face function number] = element function number
	//return:
		//number of functions on the face

	int i, lej, lfaceej, elfun = 0, facefun = 0, gf, ge, n, elf, elfun0;
	SRface* face = GetFace(lface);
	for (i = 0; i < face->GetNumNodes(); i++)
	{
		gf = face->GetNodeId(i);
		for (elfun = 0; elfun < nodeIds.GetNum(); elfun++)
		{
			ge = nodeIds.Get(elfun);
			if (ge == gf)
			{
				FaceToElFuns[facefun] = elfun;
				facefun++;
				break;
			}
		}
	}

	elfun0 = nodeIds.GetNum();

	//p2 edge funs of the face:
	for (lfaceej = 0; lfaceej < face->GetNumLocalEdges(); lfaceej++)
	{
		gf = face->GetLocalEdgeGlobalId(lfaceej);
		elfun = elfun0;
		for (lej = 0; lej < localEdges.GetNum(); lej++)
		{
			ge = localEdges.GetPointer(lej)->GetGlobalEdgeId();
			if (ge == gf)
			{
				FaceToElFuns[facefun] = elfun;
				facefun++;
			}
			elfun++;
		}
	}

	elfun0 = elfun;

	//higher edge funs of the face:
	for (lfaceej = 0; lfaceej < face->GetNumLocalEdges(); lfaceej++)
	{
		gf = face->GetLocalEdgeGlobalId(lfaceej);
		int nf = face->GetLocalEdgePOrder(lfaceej);
		elfun = elfun0;
		for (lej = 0; lej < localEdges.GetNum(); lej++)
		{
			ge = localEdges.GetPointer(lej)->GetGlobalEdgeId();
			if (ge == gf)
			{
				for (i = 3; i < nf; i++)
				{
					FaceToElFuns[facefun] = elfun;
					facefun++;
					elfun++;
				}
			}
			else
			{
				int nfl = localEdges.GetPointer(lej)->GetEdge()->GetNumGlobalFunctions() - 3;
				if (nfl > 0)
					elfun += nfl;
			}
		}
	}

	for (elf = 0; elf < lface; elf++)
	{
		SRface* elface = GetFace(elf);
		n = model.basis.CountFaceInternalFunctions(elface);
		elfun += n;
	}

	n = model.basis.CountFaceInternalFunctions(face);
	for (i = 0; i < n; i++)
	{
		FaceToElFuns[facefun] = elfun;
		facefun++;
		elfun++;
	}

	n = model.basis.CountFaceTotalFunctions(face, i);

	return facefun;
}

int SRelement::MapEdgeToElFuns(int ledge, int EdgeToElFuns[])
{
	//get basis function numbers on an edge to corresponding local element basis function numbers
	//for this element that owns the edge
	//input:
		//ledge = local edge number
	//output:
		//EdgeToElFuns[edge function number] = element function number

	SRedge* edge = GetEdge(ledge);
	int edgefun = 0;
	for (int i = 0; i < 2; i++)
	{
		int gf = edge->GetNodeId(i);
		for (int elfun = 0; elfun < nodeIds.GetNum(); elfun++)
		{
			int ge = nodeIds.Get(elfun);
			if (ge == gf)
			{
				EdgeToElFuns[edgefun] = elfun;
				edgefun++;
				break;
			}
		}
	}

	int elfun = nodeIds.GetNum();
	//midedge fn:
	for (int e = 0; e < ledge; e++)
		elfun++;
	EdgeToElFuns[edgefun] = elfun;
	elfun++;
	edgefun++;

	elfun = nodeIds.GetNum() + localEdges.GetNum();

	for (int e = 0; e < ledge; e++)
	{
		int nf = localEdges.Get(e).GetEdge()->GetNumGlobalFunctions();
		if (nf > 2)
			elfun += (nf - 2);
	}
	for (int f = 2; f < edge->GetNumGlobalFunctions(); f++)
	{
		EdgeToElFuns[edgefun] = elfun;
		elfun++;
		edgefun++;
	}

	return edge->GetNumGlobalFunctions();
}

int SRelement::MapEdgeP2FnToElFun(int ledge)
{
	//get element's local edge p2 function numbers on an local edge
	//input:
	//ledge = local edge number
	//return:
	//local element function number for the edge's p2 function

	int elfun = nodeIds.GetNum() + ledge;
	return elfun;
}

void SRelement::FillMappingNodes()
{
	//fill in mapping node vectors xnode, ynode, znode
	//from nodes and midnodes of element

	//note:
		//fills class variables xnode, ynode, znode

	int numcorner, nummidside;
	numcorner = GetNumNodes();
	nummidside = GetNumLocalEdges();
	numNodesTotal = numcorner + nummidside;

	xnodev.Allocate(numNodesTotal);
	ynodev.Allocate(numNodesTotal);
	znodev.Allocate(numNodesTotal);
	xnode = xnodev.d;
	ynode = ynodev.d;
	znode = znodev.d;
	dNdrv.Allocate(numNodesTotal);
	dNdsv.Allocate(numNodesTotal);
	dNdtv.Allocate(numNodesTotal);
	dNdr = dNdrv.d;
	dNds = dNdsv.d;
	dNdt = dNdtv.d;

	int i;
	for (i = 0; i < numcorner; i++)
		NodeLocalCopy(i, GetNodeId(i));
	for (i = 0; i < nummidside; i++)
		NodeLocalCopy(i + numcorner, GetLocalEdgeMidNodeId(i));

	//Element size; approximate as average length of all element
	//edges:
	size = 0.0;
	SRnode* n0;
	SRnode* n1;
	for (i = 0; i < nummidside; i++)
	{
		n0 = GetLocalEdgeNode(i, 0);
		n1 = GetLocalEdgeNode(i, 1);
		size += model.math.Distance(n0->Position(), n1->Position());
	}
	size /= nummidside;

	//approx volume is volume of BB. use as scale in Jacobian checks:
	double x, y, z, xmin, xmax, ymin, ymax, zmin, zmax;
	xmin = BIG;
	xmax = -BIG;
	ymin = BIG;
	ymax = -BIG;
	zmin = BIG;
	zmax = -BIG;
	for (i = 0; i < numcorner; i++)
	{
		x = xnode[i];
		y = ynode[i];
		z = znode[i];
		if (x < xmin)
			xmin = x;
		if (x > xmax)
			xmax = x;
		if (y < ymin)
			ymin = y;
		if (y > ymax)
			ymax = y;
		if (z < zmin)
			zmin = z;
		if (z > zmax)
			zmax = z;
	}
	approxVol = (xmax - xmin)*(ymax - ymin)*(zmax - zmin);
	if (approxVol < TINY)
		ERROREXIT;
}

void SRelement::approxCentroid(SRvec3& elcen)
{
	//approx centroid is average corner position
	//use this when mapping has not been setup
	
	int nn = GetNumCorners();
	for (int i = 0; i < nn; i++)
	{
		SRnode* node = GetNode(i);
		elcen.PlusAssign(node->getPos());
	}
	double scale = 1.0 / nn;
	elcen.Scale(scale);
}


void SRelement::NodeLocalCopy(int lnode, int gnode)
{
	//copy nodal coordinates into local xnode, ynode, znode vectors
	//input:
		//lnode = element local node number
		//gnode = global node id
	//note:
		//fills class variables xnode, ynode, znode

	SRnode* node;
	node = model.GetNode(gnode);
	xnode[lnode] = node->GetXyz(0);
	ynode[lnode] = node->GetXyz(1);
	znode[lnode] = node->GetXyz(2);
}

double SRelement::FillMapping(double r, double s, double t, bool detJonly)
{
	//Calculate mapping at r,s,t for an element
	//input:
		//r,s,t = natural coordinates
		//detJonly = true if only need detJ else false
	//return:
		//determinant of Jacobian of mapping
	//note:
		//fill class variable elJac.

	FillJacobian(r, s, t);

	if (detJonly)
		return elJac.Determinant();
	//invert mapping:
	double detj = elJac.Invert();
	return detj;
}

void SRelement::FillJacobian(double r, double s, double t, bool fillElJac)
{
	//fill jacobian for an element; support routine for FillMapping
	//input:
		//r,s,t = natural coordinates in mapElem
		//fillEljac = true to fill jacobian (eljac), false to just fill mapping derivatives
	//note:
		//fills class variable elJac = jacobian

	if (GetType() == tet)
		FillTetJacobian(r, s, t, fillElJac);
	else if (GetType() == wedge)
		FillWedgeJacobian(r, s, t, fillElJac);
	else if (GetType() == brick)
		FillBrickJacobian(r, s, t, fillElJac);
}

void SRelement::FillTetJacobian(double r, double s, double t, bool fillElJac)
{
	//calculate the jacobian of a tet
	//input:
		//r,s,t = natural coordinate
		//fillEljac = true to fill jacobian (eljac), false to just fill mapping derivatives
	//note:
		//fill class variables dNdr, dNds, dNdt, and elJac.

	double L[4];
	int i;
	model.map.TetVolumeCoords(r, s, t, L);

	//corner mapping functions:
	for (i = 0; i < 4; i++)
	{
		dNdr[i] = model.map.GetDLtdr(i) * (4.0*L[i] - 1.0);
		dNds[i] = model.map.GetDLtds(i) * (4.0*L[i] - 1.0);
		dNdt[i] = model.map.GetDLtdt(i) * (4.0*L[i] - 1.0);
	}


	//mid-edge mapping functions:
	double *dLtdr = model.map.GetDLtdrVec();
	double *dLtds = model.map.GetDLtdsVec();
	double *dLtdt = model.map.GetDLtdtVec();

	dNdr[4] = model.map.TetMidEdgeMappingDeriv(0, 1, dLtdr, L);
	dNds[4] = model.map.TetMidEdgeMappingDeriv(0, 1, dLtds, L);
	dNdt[4] = model.map.TetMidEdgeMappingDeriv(0, 1, dLtdt, L);
	dNdr[5] = model.map.TetMidEdgeMappingDeriv(1, 2, dLtdr, L);
	dNds[5] = model.map.TetMidEdgeMappingDeriv(1, 2, dLtds, L);
	dNdt[5] = model.map.TetMidEdgeMappingDeriv(1, 2, dLtdt, L);
	dNdr[6] = model.map.TetMidEdgeMappingDeriv(0, 2, dLtdr, L);
	dNds[6] = model.map.TetMidEdgeMappingDeriv(0, 2, dLtds, L);
	dNdt[6] = model.map.TetMidEdgeMappingDeriv(0, 2, dLtdt, L);
	dNdr[7] = model.map.TetMidEdgeMappingDeriv(0, 3, dLtdr, L);
	dNds[7] = model.map.TetMidEdgeMappingDeriv(0, 3, dLtds, L);
	dNdt[7] = model.map.TetMidEdgeMappingDeriv(0, 3, dLtdt, L);
	dNdr[8] = model.map.TetMidEdgeMappingDeriv(1, 3, dLtdr, L);
	dNds[8] = model.map.TetMidEdgeMappingDeriv(1, 3, dLtds, L);
	dNdt[8] = model.map.TetMidEdgeMappingDeriv(1, 3, dLtdt, L);
	dNdr[9] = model.map.TetMidEdgeMappingDeriv(2, 3, dLtdr, L);
	dNds[9] = model.map.TetMidEdgeMappingDeriv(2, 3, dLtds, L);
	dNdt[9] = model.map.TetMidEdgeMappingDeriv(2, 3, dLtdt, L);

	if (fillElJac)
		FillJac();
}

void SRelement::FillJac()
{
	//fill the jacobian at a natural coordinate point in an element
	//note:
		//class variables xnode, ynode, znode and dNdr, dNds, dNdr must be filled first
		//fills the class variable elJac

	elJac.Zero();
	int i;
	for (i = 0; i< numNodesTotal; i++)
	{
		elJac.rows[0].d[0] += xnode[i] * dNdr[i]; //dxdr
		elJac.rows[1].d[0] += xnode[i] * dNds[i]; //dxds
		elJac.rows[2].d[0] += xnode[i] * dNdt[i]; //dxdt
		elJac.rows[0].d[1] += ynode[i] * dNdr[i]; //dydr
		elJac.rows[1].d[1] += ynode[i] * dNds[i]; //dyds
		elJac.rows[2].d[1] += ynode[i] * dNdt[i]; //dydt
		elJac.rows[0].d[2] += znode[i] * dNdr[i]; //dzdr
		elJac.rows[1].d[2] += znode[i] * dNds[i]; //dzds
		elJac.rows[2].d[2] += znode[i] * dNdt[i]; //dzdt
	}
}

int SRelement::ShapeFunctionsLinear(double r, double s, double t, double N[])
{
	if (type == tet)
	{
		model.map.TetVolumeCoords(r, s, t, N);
		return 4;
	}
	else if (type == wedge)
	{
		model.map.WedgeLinearShapeFunctions(r, s, t, N);
		return 6;
	}
	else
	{
		model.map.BrickLinearShapeFunctions(r, s, t, N);
		return 8;
	}
}

void SRelement::PositionLinear(double r, double s, double t, SRvec3 &p)
{
	//determine element position at r,s,t using linear mapping
	//input:
		//r,s,t = natural coordinates in element
	//output:
		//p = position vector

	int i;
	double N[8];
	int nn = ShapeFunctionsLinear(r, s, t, N);
	p.Zero();
	for (i = 0; i < nn; i++)
	{
		p.d[0] += N[i] * xnode[i];
		p.d[1] += N[i] * ynode[i];
		p.d[2] += N[i] * znode[i];
	}
}

void SRelement::Position(double r, double s, double t, SRvec3 &p)
{
	//determine element position at r,s,t using quadratic mapping
	//input:
		//r,s,t = natural coordinates in element
	//output:
		//p = position vector

	int i;
	double N[20];
	ShapeFunctions(r, s, t, N);
	p.Zero();
	for (i = 0; i < numNodesTotal; i++)
	{
		p.d[0] += N[i] * xnode[i];
		p.d[1] += N[i] * ynode[i];
		p.d[2] += N[i] * znode[i];
	}
}

int SRelement::ShapeFunctions(double r, double s, double t, double N[])
{
	//calculate the shape functions for an element at a natural coordinate point
	//input:
		//r, s, t = natural coordinates
	//output:
	//N = shape functions
	//return:
		//number of functions

	double L[4], l04, l14;
	int i;
	if (GetType() == tet)
	{
		model.map.TetVolumeCoords(r, s, t, L);
		for (i = 0; i < 4; i++)
			N[i] = L[i] * (2.0*L[i] - 1.0);
		l04 = 4.0*L[0];
		l14 = 4.0*L[1];
		N[4] = l04*L[1];
		N[5] = l14*L[2];
		N[6] = l04*L[2];
		N[7] = l04*L[3];
		N[8] = l14*L[3];
		N[9] = 4.0*L[2] * L[3];
		return 10;
	}
	else if (GetType() == wedge)
		return WedgeShapeFns(r, s, t, N);
	else
		return BrickShapeFns(r, s, t, N);
}

void SRelement::XyzDerivatives(double dfdr, double dfds, double dfdt, double &dfdx, double &dfdy, double &dfdz)
{
	//Calculate x,y,z derivatives of a function given r,s,t derivatives and mapping
	//input:
		//dfdr,dfds,dfdt = derivatives of function in natural coordinates
	//output:
		//dfdx,dfdy,dfdz = derivatives of function in global coordinates
	//note:
		//must first have called FillMapping to
		//store  jacobian inverse j11, etc, in "elJac"

	dfdx = (elJac.rows[0].d[0] * dfdr) + (elJac.rows[0].d[1] * dfds) + (elJac.rows[0].d[2] * dfdt);
	dfdy = (elJac.rows[1].d[0] * dfdr) + (elJac.rows[1].d[1] * dfds) + (elJac.rows[1].d[2] * dfdt);
	dfdz = (elJac.rows[2].d[0] * dfdr) + (elJac.rows[2].d[1] * dfds) + (elJac.rows[2].d[2] * dfdt);
}

double* SRelement::FillBasisFuncs(double r, double s, double t, SRBasisCallType calltype)
{
	//calculate basis functions and derivatives
	//input:
		//r,s,t = natural coordinates
	//return:
		//basv = pointer to basisVec vector
	//note: fills class variables basisVec, etc
	model.basis.ElementBasisFuncs(r, s, t, this, basisVec, dbasisdr, dbasisds, dbasisdt, calltype);
	return basisVec;
}

void SRelement::FreeRawStrains()
{
	rawStrains.Free();
	dispEl.Free();
};


void SRelement::Cleanup()
{
	nodeIds.Free();
	localEdges.Free();
	localFaces.Free();
	globalFunctionNumbers.Free();
	rawStrains.Free();
	nodeLcSConstraints.Free();
	faceLCSConstraints.Free();
	stiffnessMatrix.Free();
	stiffDiag.Free();
	dispEl.Free();
	for (int i = 0; i < 6; i++)
		smoothedStrains[i].Free();
}

void SRelement::FillStiffDiag(int neq)
{
	//fill the element storage ids of the diagonal of this element
	//note:
		//fills class variable stiffDiag

	if (stiffDiag.num == neq)
		return;
	stiffDiag.Allocate(neq);
	int nz = 0;
	for (int r = 0; r < neq; r++)
	{
		stiffDiag.Put(r, nz);
		nz += (neq - r);
	}
}

int SRelement::GetStiffnessLocationAboveDiag(int row, int col)
{
	//given the row, column of a stiffness matrix, return the storage location
	//for a symmetric element stored as upper triangle
	//note:
		//only call this if row, col are on or above the diagonal of the element, 
		//else call the more general function GetStiffnessLocation

	return stiffDiag.Get(row) + col - row;
}

int SRelement::GetStiffnessLocation(int row, int col)
{
	//given the row, column of a stiffness matrix, return the storage location
	//for a symmetric element stored as upper triangle
	//note:
		//to avoid if statement in an inner loop, call GetStiffnessLocationAboveDiag
		//if row, col are on or above the diagonal of the element

	if (col >= row)
		return stiffDiag.Get(row) + col - row;
	else
		return stiffDiag.Get(col) + row - col;
}


void SRelement::FillBTC(double intwt, double dbdx, double dbdy, double dbdz, SRBTC& btc)
{
	//fill matrix "BTC" = strain-displacement matrixe times elasticity matrix
	//input:
		//intwt = integration weight at this gauss point
		//dbdx,dbdy,dbdz = derivatives of basis functions at this element function and gauss point
	//output:
		//fills BTC = B-T*C

	if (elMat->GetType() == iso)
	{
		double c11 = intwt * elMat->getC11();
		double lambda = intwt * elMat->getLambda();
		double G = intwt * elMat->getG();
		btc.BTC11 = c11*dbdx;
		btc.BTC12 = btc.BTC13 = lambda*dbdx;
		btc.BTC14 = G*dbdy;
		btc.BTC15 = G*dbdz;
		btc.BTC21 = btc.BTC23 = lambda*dbdy;
		btc.BTC22 = c11*dbdy;
		btc.BTC24 = G*dbdx;
		btc.BTC26 = G*dbdz;
		btc.BTC31 = btc.BTC32 = lambda*dbdz;
		btc.BTC33 = c11*dbdz;
		btc.BTC35 = G*dbdx;
		btc.BTC36 = G*dbdy;
	}
	else if (elMat->GetType() == ortho)
	{
		SRcij& cij = elMat->getOrthoCij();
		btc.BTC11 = intwt * cij.c11*dbdx;
		btc.BTC12 = intwt * cij.c12*dbdx;
		btc.BTC13 = intwt *  cij.c13*dbdx;
		btc.BTC14 = intwt *  cij.c44*dbdy;
		btc.BTC15 = intwt *  cij.c55*dbdz;
		btc.BTC21 = intwt *  cij.c12*dbdy;
		btc.BTC22 = intwt *  cij.c22*dbdy;
		btc.BTC23 = intwt *  cij.c23*dbdy;
		btc.BTC24 = intwt *  cij.c44*dbdx;
		btc.BTC26 = intwt *  cij.c66*dbdz;
		btc.BTC31 = intwt *  cij.c13*dbdz;
		btc.BTC32 = intwt *  cij.c23*dbdz;
		btc.BTC33 = intwt *  cij.c33*dbdz;
		btc.BTC35 = intwt *  cij.c55*dbdx;
		btc.BTC36 = intwt *  cij.c66*dbdy;
	}
	else
	{
		SRgenAnisoCij& gcij = elMat->getGenAnisoCij();
		btc.BTC11 = intwt * (gcij.getC(0,0) * dbdx + gcij.getC(3,0) * dbdy + gcij.getC(4,0) * dbdz);
		btc.BTC12 = intwt * (gcij.getC(0,1) * dbdx + gcij.getC(3,1) * dbdy + gcij.getC(4,1) * dbdz);
		btc.BTC13 = intwt * (gcij.getC(0,2) * dbdx + gcij.getC(3,2) * dbdy + gcij.getC(4,2) * dbdz);
		btc.BTC14 = intwt * (gcij.getC(0,3) * dbdx + gcij.getC(3,3) * dbdy + gcij.getC(4,3) * dbdz);
		btc.BTC15 = intwt * (gcij.getC(0,4) * dbdx + gcij.getC(3,4) * dbdy + gcij.getC(4,4) * dbdz);
		btc.BTC16 = intwt * (gcij.getC(0,5) * dbdx + gcij.getC(3,5) * dbdy + gcij.getC(4,5) * dbdz);
		btc.BTC21 = intwt * (gcij.getC(1,0) * dbdy + gcij.getC(3,0) * dbdx + gcij.getC(5,0) * dbdz);
		btc.BTC22 = intwt * (gcij.getC(1,1) * dbdy + gcij.getC(3,1) * dbdx + gcij.getC(5,1) * dbdz);
		btc.BTC23 = intwt * (gcij.getC(1,2) * dbdy + gcij.getC(3,2) * dbdx + gcij.getC(5,2) * dbdz);
		btc.BTC24 = intwt * (gcij.getC(1,3) * dbdy + gcij.getC(3,3) * dbdx + gcij.getC(5,3) * dbdx);
		btc.BTC25 = intwt * (gcij.getC(1,4) * dbdy + gcij.getC(3,4) * dbdx + gcij.getC(5,4) * dbdx);
		btc.BTC26 = intwt * (gcij.getC(1,5) * dbdy + gcij.getC(3,5) * dbdz + gcij.getC(5,5) * dbdz);
		btc.BTC31 = intwt * (gcij.getC(2, 0) * dbdz + gcij.getC(4, 0) * dbdx + gcij.getC(5, 0) * dbdy);
		btc.BTC32 = intwt * (gcij.getC(2,1) * dbdz + gcij.getC(4,1) * dbdx + gcij.getC(5,1) * dbdy);
		btc.BTC33 = intwt * (gcij.getC(2,2) * dbdz + gcij.getC(4,2) * dbdx + gcij.getC(5,2) * dbdy);
		btc.BTC34 = intwt * (gcij.getC(2,3) * dbdz + gcij.getC(4,3) * dbdx + gcij.getC(5,3) * dbdy);
		btc.BTC35 = intwt * (gcij.getC(2,4) * dbdz + gcij.getC(4,4) * dbdx + gcij.getC(5,4) * dbdy);
		btc.BTC36 = intwt * (gcij.getC(2,5) * dbdz + gcij.getC(4,5) * dbdx + gcij.getC(5,5) * dbdy);
	}
}

void SRelement::rescaleFlattenFraction()
{
	//rescale flattenfreaction after edges have been straightened

	if (!flattened)
		return;
	flattenFraction = 0.0;
	for (int l = 0; l < localEdges.GetNum(); l++)
	{
		SRedge* edge = GetEdge(l);
		if (edge->getStraightenFraction() > flattenFraction)
			flattenFraction = edge->getStraightenFraction();
	}
}

bool SRelement::checkForStraightenedEdges()
{
	//see if this element has any edges that were straightened to fix bad mapping
	//return:
		//true if any edges were straightened else false
	//note:
		//if any edges were straightened, the element is treated as sacrificial
		//since the geometry is not being followed properly, solution is inaccurate in this element regardless of p order

	for (int l = 0; l < localEdges.GetNum(); l++)
	{
		SRedge* edge = GetEdge(l);
		if (edge->getStraightenFraction() > approxVol*SMALL)
		{
			flattened = true;
			if (edge->getStraightenFraction() > flattenFraction)
				flattenFraction = edge->getStraightenFraction();
		}
		if (edge->getStraightenFraction() > 0.2)
		{
			FillMappingNodes();
			sacrificial = true;
			return true;
		}
	}

	if (flattenFraction > model.getMaxFlattened())
	{
		model.setMaxFlattened(flattenFraction);
		model.setMaxFlattenedElUid(userId);
	}
	//also refill the mapping nodes of the element's faces:
	for (int i = 0; i < localFaces.GetNum(); i++)
	{
		SRface* face = GetFace(i);
		face->FillMappingNodes();
	}

	return false;
}

bool SRelement::testMapping(bool okToFail)
{
	//test the mapping of this element at integration points 
	//corresponding to local p-order
	//if any unacceptable Jacobian is encountered, straighten the element edges
	//the smallest amount possible to make the mapping valid

	//input:
		//oktoFail = true if it is ok to do the straightening for recovery
			//else bad Jacobian is fatal error
	//return:
		//true if mapping is ok at all integration points else false
	//note:
		//if any straightening was needed on the element edges, class
		//variable "flattenFraction" is set to a nonzero number
		
	minJac = BIG;
	double fractions[7] = { 0.05, 0.1, 0.2, 0.4, 0.6, 0.8, 1.0 };
	int nint = model.math.FillGaussPoints(this);
	double maxStraightenFraction = 0.0;
	bool flattenedThisPass = false;
	int maxTrys = 7;
	bool mappingOk = true;
	if (!model.getPartialFlatten())
	{
		maxTrys = 1;
		fractions[0] = 1.0;
	}
	for (int gp = 0; gp < nint; gp++)
	{
		double r, s, t, w;
		model.math.GetGP3d(gp, r, s, t, w);
		double detj = FillMapping(r, s, t, true);//detj only is true
		if (detj < minJac)
			minJac = detj;
		if (detj < approxVol*SMALL)
		{
			mappingOk = false;
			flattened = true;
			flattenedThisPass = true;
			if (!okToFail)
			{
				REPPRINT(" bad Jacobian. Element: %d", userId);
				ERROREXIT;
			}

			double straightenFraction;
			for (int trynum = 0; trynum < maxTrys; trynum++)
			{
				straightenFraction = fractions[trynum];
				if (straightenFraction > maxStraightenFraction)
					maxStraightenFraction = straightenFraction;
				for (int l = 0; l < localEdges.GetNum(); l++)
				{
					SRedge* edge = GetEdge(l);
					edge->Straighten(straightenFraction);
				}
				FillMappingNodes();
				detj = FillMapping(r, s, t, true);//detj only is true
				if (detj > approxVol*SMALL)
					break;
			}
		}
	}
	if (flattenedThisPass)
	{
		flattenFraction = maxStraightenFraction;
	}
	return mappingOk;
}

void SRelement::SetBasisData()
{
	//set pointers to store basis function data for this element
	SRElementData *eldata = model.GetElementData();
	basisVec = eldata->GetBasisVec();
	dbasisdr = eldata->Getdbasisdr();
	dbasisds = eldata->Getdbasisds();
	dbasisdt = eldata->Getdbasisdt();
}

SRvec3 SRelement::GetDisp(double r, double s, double t)
{
	//calculate the displacement at a natural coordinate in an element
	//input:
		//r,s,t = natural coordinates
	//return
		//disp = 3 dof vector of displacement at r,s,t
	//note:
		//DownloadDisplacement must be called before entering loop that calls this routine
		
	SRElementData *eldata = model.GetElementData();
	basisVec = eldata->GetBasisVec();
	int nfun = globalFunctionNumbers.GetNum();
	FillBasisFuncs(r, s, t, basisonly);
	SRvec3 disp;
	for (int fun = 0; fun < nfun; fun++)
	{
		double bf = basisVec[fun];
		int gfun = globalFunctionNumbers.Get(fun);
		for (int dof = 0; dof < 3; dof++)
		{
			double u = model.GetDisplacementCoeff(gfun, dof);
			disp.d[dof] += u*bf;
		}
	}
	return disp;
}

bool SRelement::nodeDistCheck(SRvec3& pos, double radius)
{
	//check the distance between all nodes of the element and a position.
	//input:
		//pos = position
		//radius = test radius
	//return
		//true if any nodes are within radius of pos, else false
	bool anYNodeInside = false;
	for (int n = 0; n < numNodesTotal; n++)
	{
		double d = pos.Distance(xnode[n], ynode[n], znode[n]);
		if (d <= radius)
		{
			anYNodeInside = true;
			break;
		}
	}
	return anYNodeInside;
}

void SRelement::checkMaterialStressMax(double s)
{
	if (isSacrificial())
		return;
	if (s > elMat->GetMaxStress())
		elMat->PutMaxStress(s);
}

void SRelement::checkMaterialStrainMax(double e)
{
	if (sacrificial)
		return;
	if (e > elMat->GetMaxStrain())
		elMat->PutMaxStrain(e);
}

int SRelement::GetNumCorners()
{
	if (type == brick)
		return 8;
	else if (type == tet)
		return 4;
	else
		return 6;
}

int SRelement::GetLocalCornerIdFromGlobalFaceNodeId(int faceCornerId)
{
	for (int i = 0; i < GetNumCorners(); i++)
	{
		if (GetNodeId(i) == faceCornerId)
			return i;
	}
	return -1;
}

double SRelement::minCornerDistance(SRvec3& pos)
{
	double minDist = BIG;
	for (int i = 0; i < GetNumCorners(); i++)
	{
		double dist = pos.Distance(GetNode(i)->getPos());
		if (dist < minDist)
			minDist = dist;
	}
	return minDist;
}

double SRelement::centroidDistance(SRvec3& pos)
{
	SRvec3 elcen;
	approxCentroid(elcen);
	return elcen.Distance(pos);
}

int SRelement::getNodeOrMidNodeId(int n)
{
	int ncorner = nodeIds.GetNum();
	if (n < ncorner)
		return nodeIds.Get(n);
	else
	{
		int lej = n - ncorner;
		return GetLocalEdgeMidNodeId(lej);
	}
}

SRnode* SRelement::getNodeOrMidNode(int n)
{
	int nid = getNodeOrMidNodeId(n);
	return model.GetNode(nid);
}

int SRelement::localEdgeMatch(int gId)
{
	for (int i = 0; i < localEdges.GetNum(); i++)
	{
		if (GetLocalEdgeGlobalId(i) == gId)
			return i;
	}
	ERROREXIT;
	return -1;
}

bool SRelement::checkSlopeKink()
{
	//check if there is a slope kink between any adjacent boundary faces belonging to this element
	//slope kink is if outward normals to the faces near midnode of share edge are not parallel within
	//tolerance of 10 degrees
	//return:
		//true if any kink found

	double dotmin = 1.0;
	for (int f = 0; f < localFaces.GetNum(); f++)
	{
		SRface* face = GetFace(f);
		if (!face->IsBoundaryFace())
			continue;
		for (int e = 0; e < face->GetNumLocalEdges(); e++)
		{

			SRedge* edge = face->GetEdge(e);
			SRvec3 norm, pos1, v21;
			double rf, sf;
			face->NaturalCoordinatesNearMidedge(e, rf, sf);
			face->Position(rf, sf, pos1);
			face->OutwardNormal(rf, sf, norm);
			for (int bf = 0; bf < edge->getNumBoundaryfaceData(); bf++)
			{
				int fid = edge->getBoundaryFaceData(bf).faceId;
				if (fid == face->GetId())
					continue;
				int lej = edge->getBoundaryFaceData(bf).localEdgeId;
				SRface* face2 = model.GetFace(fid);
				SRvec3 norm2, pos2;
				face2->NaturalCoordinatesNearMidedge(lej, rf, sf);
				face2->Position(rf, sf, pos2);
				face2->OutwardNormal(rf, sf, norm2);
				pos1.Subtract(pos2, v21);
				v21.Normalize();
				double normdot = norm.Dot(norm2);
				if (normdot < 0.0)
					continue;
				if ((norm2.Dot(v21) > 0.0) && (norm.Dot(v21) < 0.0))
				{
					//reentrant if vector from point on face2 to point on face 1
					//points in same direction as normal to face2 and in opposite direction to normal to face1
					if (normdot < dotmin)
					{
						dotmin = normdot;
					}
				}
			}
		}
	}

	if (dotmin < 0.984)//tol ~10 degrees
		return true;
	else
		return false;
}

void SRelement::DownloadSmoothedStrains(double *strainvecAllComps[6])
{
	//download smoothed strains from global solution to local strains of this element
	//note:
	//fills up class variable smoothedStrains

	int n = globalFunctionNumbers.GetNum();
	int c, i, f, eq;
	for (i = 0; i < 6; i++)
		smoothedStrains[i].Allocate(n);
	for (c = 0; c < 6; c++)
	{
		double *strainVec = strainvecAllComps[c];
		for (i = 0; i < n; i++)
		{
			f = globalFunctionNumbers.Get(i);
			eq = model.GetSmoothFunctionEquation(f);
			if (eq >= 0)
				smoothedStrains[c].Put(i, strainVec[eq]);
		}
	}
}

void SRelement::GetSmoothedStrain(double r, double s, double t, double strain[])
{
	//look up smoothed strain in an element at a point
	//input:
	//r,s,t = natural coordinates
	//output:
	//strain[6] = vector of strains at this location

	int i, n, c;
	double* basis = basisVec;

	for (c = 0; c < 6; c++)
	{
		strain[c] = 0.0;
		n = model.basis.ElementBasisFuncs(r, s, t, this, basis);
		for (i = 0; i < n; i++)
			strain[c] += (smoothedStrains[c].Get(i)*basis[i]);
	}
}

bool SRelement::hasLcsConstraint()
{
	return (nodeLcSConstraints.GetNum() != 0 || faceLCSConstraints.GetNum() != 0);
};

void SRelement::AddfaceLCSConstraint(int lface)
{
	faceLCSConstraints.pushBack(lface);
};

int SRelement::GetNumNodes()
{
	return nodeIds.GetNum();
};

int SRelement::GetNumNodesTotal()
{
	return nodeIds.GetNum() + localEdges.GetNum();

};

int SRelement::GetNodeId(int i)
{
	return nodeIds.Get(i);
};

int SRelement::GetNumLocalEdges()
{
	return localEdges.GetNum();
};

int SRelement::GetLocalEdgeGlobalId(int i)
{
	return localEdges.Get(i).GetGlobalEdgeId();
};

int SRelement::GetLocalEdgePOrder(int i)
{
	return localEdges.Get(i).GetPOrder();
};

int SRelement::GetLocalEdgeMidNodeId(int i)
{
	return localEdges.Get(i).GetMidNodeId();
};

int SRelement::GetMaterialId()
{
	return elMat->getId();
};

int SRelement::GetLocalEdgeDirection(int lej)
{
	return localEdges.Get(lej).GetDirection();
};

int SRelement::GetNumLocalFaces()
{
	return localFaces.GetNum(); 
};

int SRelement::GetLocalFaceGlobalId(int i)
{
	return localFaces.Get(i).GetGlobalFaceId(); 
};

int SRelement::GetLocalFaceGlobalNodeOrder(int lface, int lnode)
{
	return localFaces.Get(lface).globalNodeOrder[lnode]; 
};

int SRelement::GetNumFunctions()
{
	return globalFunctionNumbers.GetNum(); 
};

int SRelement::GetFunctionNumber(int i)
{
	return globalFunctionNumbers.Get(i); 
};

int SRelement::getStiffDiag(int eq)
{
	return stiffDiag.Get(eq); 
};

SRlocalFace* SRelement::GetLocalFace(int lface)
{
	return localFaces.GetPointer(lface); 
};

SRmaterial* SRelement::GetMaterial()
{
	return elMat; 
};

void SRelement::AllocateRawStrains(int n)
{
	rawStrains.Allocate(n, 6); 
};

void SRelement::AllocateGlobalFunctionNumbers(int n)
{
	globalFunctionNumbers.Allocate(n); 
};

void SRelement::PutGlobalFunctionNumbers(int lfun, int gfun)
{
	globalFunctionNumbers.Put(lfun, gfun); 
};

void SRelement::AssignLocalEdge(int lej, int gej, int direction)
{
	localEdges.GetPointer(lej)->Assign(gej, direction); 
};

void SRelement::AddNodeLcSConstraints(int l)
{
	nodeLcSConstraints.pushBack(l); 
};

void SRelement::SetMaterial(SRmaterial* mat)
{
	elMat = mat; 
};

void SRelement::deleteStiffnessMaxtrix()
{
	stiffnessMatrix.Free(); 
};

double SRelement::GetRawStrain(int i, int e)
{
	return rawStrains.Get(i, e); 
};

double* SRelement::GetStiffnessMatrix()
{
	return stiffnessMatrix.GetVector();
};

double SRelement::getXnode(int n)
{
	return xnode[n];
};

double SRelement::getYnode(int n)
{
	return ynode[n];
};

double SRelement::getZnode(int n)
{
	return znode[n];
};

double SRelement::getdNdr(int i)
{
	return dNdr[i];
};

double SRelement::getdNds(int i)
{
	return dNds[i];
};

double SRelement::getdNdt(int i)
{
	return dNdt[i];
};

SRnode* SRelement::GetLocalEdgeNode(int lej, int localnodenum)
{
	return localEdges.Get(lej).GetEdge()->GetNode(localnodenum);
};

