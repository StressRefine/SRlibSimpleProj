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
// SRforce.cpp: implementation of the SRforce class.
//
//////////////////////////////////////////////////////////////////////

#include "SRmodel.h"

#ifdef _DEBUG
#undef THIS_FILE
static char THIS_FILE[]=__FILE__;
#endif

extern SRmodel model;

SRforce::SRforce()
{
	pressure = false;
	coordId = -1;
	elemId = -1;
	for (int i = 0; i < 4; i++)
		nv[i] = -1;
	entityId = -1;
	faceFromNodal = false;
}

void SRforce::Copy(SRforce& that, bool copyForceVals)
{
	//copy force "that" to this one
	//input:
		//that = force to copy
		//copyForceVals = true to copy data as well as attributes

	type = that.type;
	pressure = that.isPressure();
	coordId = that.coordId;
	entityId = that.entityId;
	if (copyForceVals)
		forceVals.Copy(that.forceVals);
}

void SRforce::Clear()
{
	type = nodalForce;
	coordId = -1;
	entityId = -1;
	elemId = -1;
	pressure = false;
	forceVals.Free();
}

void SRforce::AddNodalForce(SRforce& that, bool SummingForceVals)
{
	//append force values from force "that" to this one if they
	//differ or "SummingForceVals" is true

	int ndof = 3;
	if (pressure)
		ndof = 1;
	for (int dof = 0; dof < ndof; dof++)
	{
		double f = GetForceVal(0, dof);
		double thatf = that.GetForceVal(0, dof);
		if ((fabs(f - thatf) > TINY) || SummingForceVals)
		{
			f += thatf;
			forceVals.Put(0, dof, f);
		}
	}
}

void SRvolumeForce::GetForceValue(SRelement *elem,SRvec3 &p,double val[])
{
    //return the value of a volume force in an element at position p
    //input:
        //elem = pointer to element
        //p = position
    //output
        //val = 3 dof force values at p
    //note:
		//p is only used for centrifugal, gravity is not dependent on position

	double rho, rhoOmega2;
	int mid = elem->GetMaterialId();
	SRmaterial *mat = model.GetMaterial(mid);
	rho = mat->GetRho();
	SRvec3 R, omR;
	if (type == gravity)
	{
		val[0] = rho*g1;
		val[1] = rho*g2;
		val[2] = rho*g3;
	}
	else if (type == centrifugal)
	{
		rhoOmega2 = rho*omega2;
		p.Subtract(origin, R);
		axis.Cross(R, omR);
		axis.Cross(omR, R);
		val[0] = -rhoOmega2*R.d[0];
		val[1] = -rhoOmega2*R.d[1];
		val[2] = -rhoOmega2*R.d[2];
		val[0] -= alpha*omR.d[0];
		val[1] -= alpha*omR.d[1];
		val[2] -= alpha*omR.d[2];
	}
}

void SRthermalForce::Process()
{
	//fill up contribution of thermal force to the global force vector

	int el, i, fun, gfun, geq, gp, dof;
	int nint;
	double* globalForce;
	double alphax, alphay, alphaz, r, s, t, w, detj, dbdx, dbdy, dbdz, dhdr, dhds, dhdt;
	double temprst, ceT[6];
	SRdoubleVector bv, dbdr, dbds, dbdt;
	bv.Allocate(model.GetmaxNumElementFunctions());
	dbdr.Allocate(model.GetmaxNumElementFunctions());
	dbds.Allocate(model.GetmaxNumElementFunctions());
	dbdt.Allocate(model.GetmaxNumElementFunctions());
	double* basisVec = bv.GetVector();
	double* dbasisdr = dbdr.GetVector();
	double* dbasisds = dbds.GetVector();
	double* dbasisdt = dbdt.GetVector();
	bool fullceT;
	double eTx,eTy,eTz;
	SRvec3 bTceT;
	SRelement* elem;
	SRmaterial* mat;
	globalForce = model.getSolutionVector();
	for(el = 0; el < model.GetNumElements(); el++)
	{
		elem = model.GetElement(el);
		mat = elem->GetMaterial();
		if(constantTemp)
		{
			double deltaTemp = temp - mat->getTref();
			eTx = mat->getAlphax()*deltaTemp;
			eTy = mat->getAlphay()*deltaTemp;
			eTz = mat->getAlphaz()*deltaTemp;
			//compute C*eT:
			fullceT = CeTMult(mat, eTx, eTy, eTz, ceT);
		}
		nint = model.math.FillGaussPoints(elem);
		for(gp = 0; gp < nint; gp++)
		{
			model.math.GetGP3d(gp, r, s, t, w);
			detj = elem->FillMapping(r, s, t);
			w *= detj;
			if(!constantTemp)
			{
				temprst = GetTemp(elem, r, s, t);
				double deltaTemp = temprst - mat->getTref();
				eTx = mat->getAlphax()*deltaTemp;
				eTy = mat->getAlphay()*deltaTemp;
				eTz = mat->getAlphaz()*deltaTemp;
				//compute C*eT:
				fullceT = CeTMult(mat, eTx, eTy, eTz, ceT);
			}
			model.basis.ElementBasisFuncs(r, s, t, elem, basisVec, dbasisdr, dbasisds, dbasisdt, derivonly);
			for (fun = 0; fun < elem->GetNumFunctions(); fun++)
			{
				dhdr = dbasisdr[fun];
				dhds = dbasisds[fun];
				dhdt = dbasisdt[fun];
				elem->XyzDerivatives(dhdr, dhds, dhdt, dbdx, dbdy, dbdz);
				//compute B-transpose*C*eT:
				bTceT.Zero();
				bTceT.d[0] = dbdx*ceT[0];
				bTceT.d[1] = dbdy*ceT[1];
				bTceT.d[2] = dbdz*ceT[2];
				if(fullceT)
				{
					bTceT.d[0] += (dbdy*ceT[3] + dbdz*ceT[4]);
					bTceT.d[1] += (dbdx*ceT[3] + dbdz*ceT[5]);
					bTceT.d[2] += (dbdx*ceT[4] + dbdy*ceT[5]);
				}
				gfun = elem->GetFunctionNumber(fun);
				for(dof = 0; dof < 3; dof++)
				{
					geq = model.GetFunctionEquation(gfun, dof);
					if (geq >= 0)
						globalForce[geq] += (bTceT.d[dof] * w);
				}
			}
		}
	}
}

bool SRthermalForce::CeTMult(SRmaterial *mat, double eTx, double eTy, double eTz, double ceT[])
{
	//multiply C-matrix times eT for this thermal force
	//input::
		//mat = pointer to material
		//eTx = alphax-T=thermal normal strain
		//eTy = alphay-T=thermal normal strain
		//eTz = alphaz-T=thermal normal strain
	//output:
		//ceT[6] = C*eT
	//return:
		//fullCeT indicates whether Ce[3,4,5] are 0

	double C[3][3],eT[3];
	int i,j;
	eT[0] = eTx;
	eT[1] = eTy;
	eT[2] = eTz;
	for(i = 0;i < 6; i++)
		ceT[i] = 0.0;
	if (mat->GetType() == iso)
	{
		double c11 = mat->getC11();
		double lambda = mat->getLambda();
		double G = mat->getG();
		C[0][0] = c11;
		C[0][1] = lambda;
		C[0][2] = lambda;
		C[1][0] = lambda;
		C[1][1] = c11;
		C[1][2] = lambda;
		C[2][0] = lambda;
		C[2][1] = lambda;
		C[2][2] = c11;
	}
	else if (mat->GetType() == ortho)
	{
		SRcij& ocij = mat->getOrthoCij();
		C[0][0] = ocij.c11;
		C[0][1] = ocij.c12;
		C[0][2] = ocij.c13;
		C[1][0] = ocij.c12;
		C[1][1] = ocij.c22;
		C[1][2] = ocij.c23;
		C[2][0] = ocij.c13;
		C[2][1] = ocij.c23;
		C[2][2] = ocij.c33;
	}
	else if (mat->GetType() == genAniso)
	{
		for(i = 0; i < 6; i++)
		{
			ceT[i] = 0.0;
			for(j = 0; j < 3; j++)
				ceT[i] += (mat->getGenAnisoCij().getC(i,j) * eT[j]);
		}
		return true;
	}

	for(i = 0; i < 3; i++)
	{
		for(j = 0; j < 3; j++)
			ceT[i] += (C[i][j] * eT[j]);
	}
	return false;
}

double SRthermalForce::GetTemp(SRelement *elem, double r, double s, double t)
{
	//get temperature in an element due to this thermal force at natural coordinate point
	//input:
		//elem = pointer to element
		//r,s,t = natural coordinates
	//return:
		//temperature

	int i, nn, ne, id;
	double temprst, tempi;
	SRedge* edge;
	if(constantTemp)
		return temp;
	temprst = 0.0;
	double N[20];
	elem->ShapeFunctions(r, s, t, N);
	nn = elem->GetNumNodes();
	ne = elem->GetNumLocalEdges();
	for(i = 0; i < nn; i++)
	{
		id = elem->GetNodeId(i);
		tempi = nodalTemp.Get(id);
		temprst += (tempi*N[i]);
	}
	for(i = 0; i < ne; i++)
	{
		edge = elem->GetEdge(i);
		id = edge->GetMidNodeId();
		tempi = nodalTemp.Get(id);
		temprst += (tempi*N[i+nn]);
	}
	return temprst;
}

