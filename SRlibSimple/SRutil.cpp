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
// SRutil.cpp: implementation of the SRutil class.
//
//////////////////////////////////////////////////////////////////////

#include "SRmodel.h"

extern SRmodel model;

static bool errorExitCalled = false;

void SRutil::ErrorExit(const char *file, int line)
{
	//print error messages and shut down
	//when a fatal error occurs
	//input:
	//file = filename where error occurred
	//line = line number where error occurred

	if (errorExitCalled)
		exit(0); //prevent recursion
	errorExitCalled = true;
	REPPRINT("StressRefine library abnormal termination");

	exit(0);
}

void SRintVector::PushBack(int v)
{
	//append integer "v" to the end of this vector

	SRintVector tmp;
	tmp.Copy(*this);
	num++;
	Allocate(num);
	for (int i = 0; i < tmp.GetNum(); i++)
		d[i] = tmp.d[i];
	d[num - 1] = v;
};

void SRdoubleVector::PushBack(double v)
{
	//append double "v" to the end of this vector

	SRdoubleVector tmp;
	tmp.Copy(*this);
	num++;
	Allocate(num);
	for (int i = 0; i < tmp.GetNum(); i++)
		d[i] = tmp.d[i];
	d[num - 1] = v;
};

void SRdoubleMatrix::PlusAssign(SRdoubleMatrix& that)
{
	//add contains of matrix "that" to this matrix
	//if this matrix has not yet been allocated, 
	//zero it and assign the contents of "that" to it

	if (n==0)
		Allocate(that.n, that.m);
	else
	{
		if (n != that.n || m != that.m)
			ERROREXIT;
	}
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < m; j++)
			d[i][j] += that.d[i][j];
	}
}

void SRdoubleMatrix::Copy(SRdoubleMatrix& that)
{
	//copy contains of matrix "that" into this matrix

	Allocate(that.n, that.m);
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < m; j++)
			d[i][j] = that.d[i][j];
	}
}

int SRIntCompareFunc(const void *v1, const void *v2)
{
	//compare function for sorting SRintVector

	return *((int *)v1) - *((int *)v2);
}

int SRintVector::Find(int intIn)
{
	//find an integer in this vector using binary search
	//input:
		//intIn = integer value
	//return:
		//location where intIn resides in this vector, -1 if not found
	//note:
		//Sort must be called before using this routine

	int* intOutPtr = (int *)bsearch(&intIn, d, num, sizeof(int), SRIntCompareFunc);
	if (intOutPtr == NULL)
		return -1;
	else
		return intOutPtr - d;
}

void SRintVector::Sort()
{
	//binary sort this vector

	qsort((void *)d, num, sizeof(int), SRIntCompareFunc);
}
