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
// SRnode.cpp: implementation of the SRnode class.
//
//////////////////////////////////////////////////////////////////////

#include "SRmodel.h"

extern SRmodel model;

SRnode::SRnode()
{
	globalFunctionNumber = -1;
	userId = -1;
	firstElementOwner = -1;
	midSideEdgeOwner = -1;
	constraintId = -1;
	sacrificial = false;
	dispCoordid = -1;
}

SRconstraint* SRnode::GetConstraint()
{
	//return pointer to constraint acting on this node, NULL if not constrained
	if (constraintId == -1)
		return NULL;
	else
		return model.GetConstraint(constraintId);
}

bool SRnode::GetDisp(SRvec3& disp)
{
	//look up the displacement result for this node
	//output:
		//disp = displacement vector for this node
	//return:
		//true if node has displacement

	if (isOrphan())
	{
		disp.Zero();
		return false;
	}
	if (isMidSide())
	{
		SRedge* edge = model.GetEdge(midSideEdgeOwner);
		edge->getDisp(0.0, disp);
		return true;
	}
	int gfun = globalFunctionNumber;
	for (int dof = 0; dof < 3; dof++)
	{
		double d = model.GetDisplacementCoeff(gfun, dof);
		disp.d[dof] = d;
	}
	return true;
}

