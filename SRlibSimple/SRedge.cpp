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
// SRedge.cpp: implementation of the SRedge class.
//
//////////////////////////////////////////////////////////////////////

#include "SRmodel.h"

extern SRmodel model;

//////////////////////////////////////////////////////////////////////////
// SRedgeUtil static functions are used during creation of global edges //
//////////////////////////////////////////////////////////////////////////


bool SRedgeUtil::EdgeMatch(int n1, int n2, int g1, int g2, int &direction)
{
	//see if edge with nodes n1,n2 matches edges with nodes g1,g2
	//input:
		//n1, n2 = nodes at end of edge being tested
		//g1, g2 = nodes at end of edge being compared against
	//output:
		//direction = 1 if edges have same orientation else -1
	//return:
		//true if edges match else false
	if (n1 == g1 && n2 == g2)
	{
		direction = 1;
		return true;
	}
	if (n1 == g2 && n2 == g1)
	{
		direction = -1;
		return true;
	}
	return false;
}

int SRedgeUtil::GlobalEdgeMatchN2(int n1, int n2, int &direction)
{
	//find which global edges has nodes n1,n2.
	//input:
		//n1,n2 = nodes at ends of an edges
	//output:
		//direction = 1 if edges direction matches global else -1
	//return:
		//edge id; -1 if not found
	int ej, g1, g2;
	for (ej = 0; ej < model.GetNumEdges(); ej++)
	{
		g1 = model.GetEdge(ej)->nodeIds[0];
		g2 = model.GetEdge(ej)->nodeIds[1];
		if (SRedgeUtil::EdgeMatch(n1, n2, g1, g2, direction))
			return ej;
	}

	return -1;
}

int SRedgeUtil::GlobalEdgeMatch(int n1, int n2, int &direction)
{
	//find which global edge has nodes n1,n2.
	//input:
		//n1,n2 = nodes at ends of an edge
	//output:
		//direction = 1 if edges direction matches global else -1
	//return:
		//edge id; -1 if not found
	int ej, g1, g2;
	int nej = model.GetNumNodeEdges(n1);
	if (nej == -1)
	{
		//for global node n1, it either has no nodeEdges assigned yet,
		//or it overflowed its nodeEdges array, need to search all model edges:
		for (ej = 0; ej < model.GetNumEdges(); ej++)
		{
			SRedge* edge = model.GetEdge(ej);
			g1 = edge->nodeIds[0];
			g2 = edge->nodeIds[1];
			if (SRedgeUtil::EdgeMatch(n1, n2, g1, g2, direction))
				return ej;
		}
	}
	else
	{
		for (ej = 0; ej < nej; ej++)
		{
			int nodeEj = model.GetNodeEdge(n1, ej);
			SRedge* edge = model.GetEdge(nodeEj);
			g1 = edge->nodeIds[0];
			g2 = edge->nodeIds[1];
			if (SRedgeUtil::EdgeMatch(n1, n2, g1, g2, direction))
				return nodeEj;
		}
	}

	return -1;
}

///////////////////////////////////////////////////////////////
//             SRlocalEdge member functions                  //
///////////////////////////////////////////////////////////////

SRlocalEdge::SRlocalEdge()
{
	globalEdgeId = -1;
	direction = -1;
}

SRedge* SRlocalEdge::GetEdge()
{
	return model.GetEdge(globalEdgeId);
}

int SRlocalEdge::GetPOrder()
{
	int p = GetEdge()->GetPorder();
	return p;
}

int SRlocalEdge::GetMidNodeId()
{
	return model.GetEdge(globalEdgeId)->GetMidNodeId();
}


///////////////////////////////////////////////////////////////
//               SRedge member functions                     //
///////////////////////////////////////////////////////////////

SRedge::SRedge()
{
	pOrder = 0;
	straight = false;
	thin = false;
	straightenFraction = 0.0;
	sacrificial = false;
	initialBow = -1.0;
}

void SRedge::Create(int n1, int n2, int& midnode, int idt)
{
	//input:
		//n1, n2 are corner node ids
		//midnode = mid node id
		//idt = id of this edge

	id = idt;
	nodeIds[0] = n1;
	nodeIds[1] = n2;
	SRvec3 p1 = model.GetNode(n1)->Position();
	SRvec3 p2 = model.GetNode(n2)->Position();
	constraintId = -1;
	SRvec3 pmid, pave;
	size = model.math.Distance(p1, p2);
	pave = p1;
	pave.PlusAssign(p2);
	pave.Scale(0.5);
	if (midnode == -1)
	{
		midnode = model.GetNumNodes();
		SRnode *node = model.addNode();
		int uid = model.getMaxNodeUid() + 1;
		model.setAllFacesFlat(uid);
		node->Create(uid, pave.d[0], pave.d[1], pave.d[2]);
		node->setId(midnode);
		pmid = pave;
	}
	else
		pmid = model.GetNode(midnode)->Position();
	midnodeId = midnode;
	double bow = pave.Distance(pmid);
	straight = true;
	if (bow > RELSMALL*size)
	{
		straight = false;
		//better estimate of curved size:
		size = model.math.Distance(p1, pmid);
		size += model.math.Distance(p2, pmid);
	}
	prevpOrder = 2;//min p order
	int n = model.GetNumNodeEdges(n1);
	if (n != -1)
	{
		if (n == 20)
			model.PutNumNodeEdges(n1, -1);
		else
		{
			model.PutNodeEdge(n1, n, id);
			n++;
			model.PutNumNodeEdges(n1, n);
		}
	}
	n = model.GetNumNodeEdges(n2);
	if (n != -1)
	{
		if (n == 20)
			model.PutNumNodeEdges(n2, -1);
		else
		{
			model.PutNodeEdge(n2, n, id);
			n++;
			model.PutNumNodeEdges(n2, n);
		}
	}
}

int SRedge::GetNodeId(int i)
{
	if (i == 2)
		return midnodeId;
	else
		return nodeIds[i]; };


void SRedge::putPorder(int p)
{
	//set this edges p order to p
	//input:
		//p = new p order to set
	//note:
		//for sacrifical edges, do not exceed 2

	prevpOrder = pOrder;
	pOrder = p;
	if (sacrificial)
	{
		if (pOrder > 2)
			pOrder = 2;
	}
};


void SRedge::AssignGlobalFunctionNumbers(int& fun, int pmin, int pmax)
{
	//assign global function numbers to an edge
	//input
		//fun = first function to assign
		//pmin = min p order to assign function number to
		//pmax = max p order to assign function number to
	// (mod) fun = next function to be assigned

	int pej = MATHMIN(pOrder, pmax);
	int gfun = GetNode(0)->GetGlobalFunctionNumber();
	globalFunctionNumbers.Put(0, gfun);
	gfun = GetNode(1)->GetGlobalFunctionNumber();
	globalFunctionNumbers.Put(1, gfun);
	for (int j = pmin; j <= pej; j++)
	{
		globalFunctionNumbers.Put(j, fun);
		fun++;
	}
}


SRconstraint * SRedge::GetConstraint()
{
	//get constraint associated with this edge
	//return:
		//pointer to constraint, NULL if none.
	if (constraintId == -1)
		return NULL;
	else
		return model.GetConstraint(constraintId);
}



SRnode* SRedge::GetNode(int localnodenum)
{
	//get global node corresponding to local node of this edge
	//input:
		//localnodenum = local node number
	//return:
		//pointer to global node

	int nodeid = nodeIds[localnodenum];
	return model.GetNode(nodeid);
};


void SRedge::GetForceValue(double r, SRforce *force, double forceVal[])
{
	//fill up 3 dof force vector at natural coordinate r for an edge force
	//input:
		//r = natural coordinate on edge
		//force = pointer to force
	//output:
		//forceVal = 3 dof force vector

	double N[3];
	N[2] = 1.0 - r*r;
	N[0] = 0.5*(1.0 - r) - 0.5*N[2];
	N[1] = 0.5*(1.0 + r) - 0.5*N[2];
	for(int dof = 0; dof < 3; dof++)
	{
		forceVal[dof] = 0.0;
		for (int j = 0; j < 3; j++)
		{
			double ft = force->GetForceVal(j, dof);
			forceVal[dof] += N[j] * ft;
		}
	}
}


void SRedge::Clear()
{
	//free memory for an edge that is no longer needed
	globalFunctionNumbers.Free();
}

void SRedge::GetDisplacement(double r, SRvec3 &disp)
{
	//get displacement vector at a natural coordinate on an edge.
	//input:
		//r = natural coordinate
	//output:
		//disp = displacement vector

	int i, dof, fun;
	double d;
	double basis1d[10];
	disp.Zero();
	basisFunctions(r, basis1d);
	for (i = 0; i <globalFunctionNumbers.GetNum(); i++)
	{
		fun = globalFunctionNumbers.Get(i);
		for(dof = 0; dof < 3; dof++)
		{
			d = model.GetDisplacementCoeff(fun, dof);
			disp.d[dof] += (d*basis1d[i]);
		}
	}
}

void SRedge::FillMappingNodes()
{
	//fill in mapping nodes for an edge
	//note:
		//fills class variables xnode, ynode, znode, and size

	SRnode* node;
	int i;
	for (i = 0; i < 2; i++)
	{
		node = GetNode(i);
		xnode[i] = node->GetXyz(0);
		ynode[i] = node->GetXyz(1);
		znode[i] = node->GetXyz(2);
	}
	i = GetMidNodeId();
	node = model.GetNode(i);
	xnode[2] = node->GetXyz(0);
	ynode[2] = node->GetXyz(1);
	znode[2] = node->GetXyz(2);
}

void SRedge::Position(double r, SRvec3& p)
{
	//calculate the position on an edge at a natural coodinate
	//input:
		//r = natural coodinate
	//output:
		//p =  position of centroid

	double N[3];
	int i, n;
	p.Zero();
	n = model.map.EdgeShapeFunctions(r, N);
	for (i = 0; i < n; i++)
	{
		p.d[0] += xnode[i] * N[i];
		p.d[1] += ynode[i] * N[i];
		p.d[2] += znode[i] * N[i];
	}
}

void SRedge::basisFunctions(double r, double* basis)
{
	//calculate edge basis functions at a natural coordinate
	//input:
		//r = natural coordinate on edge
	//output:
		//basis = vector of basis functions for 
	model.basis.Basisfuns1d(pOrder, r, basis);
}

void SRedge::Straighten(double fraction)
{
	//straighten edge to fix element mapping issues, e.g. for slivery elements.

	if (straight)
		return;

	if (fraction <= straightenFraction*(1.0 + SMALL))
		return;

	SRvec3 p1 = GetNode(0)->getPos();
	SRvec3 p2 = GetNode(1)->getPos();
	SRnode* midnode = model.GetNode(midnodeId);
	SRvec3 pmid;
	pmid.Copy(midnode->getPos());
	SRvec3 pAve;
	pAve.Copy(p1);
	pAve.PlusAssign(p2);
	pAve.Scale(0.5);
	double bow = pmid.Distance(pAve);
	if (initialBow < 0.0)
		initialBow = bow;
	SRvec3 toAve;
	pAve.Subtract(pmid, toAve);
	toAve.Scale(fraction);
	pmid.PlusAssign(toAve);

	midnode->setPos(pmid);
	straightenFraction = fraction;
	if (straightenFraction > 0.2)
		sacrificial = true;
}

void SRedge::ScaleStraightenVsEdgeLength()
{
	//scale straightenFraction so it improves with mesh refinement:
	//note:
		//updates class variable straightenFraction
		//(straightenFraction is 0 for not straightened, 1 for completely straightened)

	if (straightenFraction < SMALL)
		return;
	SRvec3 p1 = model.GetNode(nodeIds[0])->getPos();
	SRvec3 p2 = model.GetNode(nodeIds[1])->getPos();
	double chordLength = model.math.Distance(p1, p2);
	//scale straightenFraction so it improves with mesh refinement:
	straightenFraction *= (initialBow / chordLength);
}

void SRedge::getDisp(double r, SRvec3& disp)
{
	//calculate the displacement at a natural coordinate on an edge
	//input:
		//r = natural coordinat
	//output
		//disp = 3 dof displacement vector at r
	disp.Zero();
	double basis[9];
	basisFunctions(r, basis);
	int nfun = globalFunctionNumbers.GetNum();
	for (int fun = 0; fun < nfun; fun++)
	{
		double bf = basis[fun];
		int gfun = globalFunctionNumbers.Get(fun);
		for (int dof = 0; dof < 3; dof++)
		{
			double u = model.GetDisplacementCoeff(gfun, dof);
			disp.d[dof] += u*bf;
		}
	}
}

int SRedge::GetMidNodeUserId()
{
	if (midnodeId != -1)
		return model.GetNode(midnodeId)->GetUserid();
	else
		return -1;
}

bool SRedge::isOrphan()
{
	if (GetNode(0)->isOrphan() && GetNode(1)->isOrphan())
		return true;
	else return false;
}

void SRedge::ProcessForce(SRforce *force, SRvec3& resF)
{
	//fill up contribution of an edge force to global force vector
	//input:
		//force=pointer to force
		//ResF = resultant force vector on model
	//output:
		//ResF = updated

	int i,dof, nint, gp, eq, gfun;
	SRgaussPoint* g1d;
	double* globalForce = model.getSolutionVector();
	double r, w, basis[10], ds, bw, forceVal[3];
	nint = pOrder + 1;
	g1d = model.math.GetGaussPoints1d(nint);
	FillMappingNodes();
	SRvec3 tmpRes;
	for (gp = 0; gp < nint; gp++)
	{
		r = g1d[gp].x;
		w = g1d[gp].w;
		model.basis.Basisfuns1d(pOrder, r, basis);
		ds = model.map.EdgeArcLength(this, r);
		GetForceValue(r, force, forceVal);
		w *= ds;
		tmpRes.Assign(forceVal);
		tmpRes.Scale(w);
		resF.PlusAssign(tmpRes);
		for (i = 0; i < globalFunctionNumbers.GetNum(); i++)
		{
			gfun = globalFunctionNumbers.Get(i);
			bw = basis[i] * w;
			for (dof = 0; dof < 3; dof++)
			{
				eq = model.GetFunctionEquation(gfun, dof);
				if (eq >= 0)
					globalForce[eq] += (bw*forceVal[dof]);
			}
		}
	}
}

void SRedge::PutBoundaryFaceId(int f, int lej)
{
	//add a face id to boundary face list for this edge
	//input:
		//f = global face id
		//lej = local edge number
	//note:
		//appends f to boundaryfaceData list for this edge
		//if not already on it

	bool found = false;
	for (int i = 0; i < boundaryfaceData.GetNum(); i++)
	{
		if (boundaryfaceData.Get(i).faceId == f)
		{
			found = true;
			break;
		}
	}
	if(!found)
	{
		SRboundaryFaceData tmp;
		tmp.faceId = f;
		tmp.localEdgeId = lej;
		boundaryfaceData.pushBack(tmp);
	}
}

bool SRedge::checkKinkOK(SRface* face0, int lej0, SRface* face1, int lej1)
{
	//check if kink in normals to two faces that share this edge is less than tolerance
	//input:
		//face0 = 1st face number
		//lej0 = local edge number on face0 corresponding to this edge
		//face1 = 2nd face number
		//lej1 = local edge number on face1 corresponding to this edge
	//return:
		//true if kink in normals is within tolerance, else false

	double r, s;
	SRvec3 norm0, norm1;
	face0->NaturalCoordinatesNearMidedge(lej0, r, s);
	face0->OutwardNormal(r, s, norm0);
	face1->NaturalCoordinatesNearMidedge(lej1, r, s);
	face1->OutwardNormal(r, s, norm1);
	double dot = norm0.Dot(norm1);
	if (dot > 0.9848)
		return true;//normal kink is less than tolerance 
	return false;
}
