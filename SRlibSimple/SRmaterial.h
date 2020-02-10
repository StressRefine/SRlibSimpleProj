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
// SRmaterial.h: interface for the SRmaterial class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(SRMATERIAL_INCLUDED)
#define SRMATERIAL_INCLUDED

#include "SRstring.h"


enum SRmaterialType {iso, ortho, genAniso};

struct SRcij
{
	double c11;
	double c12;
	double c13;
	double c22;
	double c23;
	double c33;
	double c44;
	double c55;
	double c66;
};

class SRgenAnisoCij
{
	friend class SRmaterial;
public:
    bool symCheck();
	inline double getC(int i, int j){ return c[i][j]; };
	inline void setC(int i, int j, double v){ c[i][j] = v; };
private:
	double c[6][6];
};


class SRmaterial
{
public:
	SRmaterial();
	~SRmaterial(){};
	void setId(int i){ id = i; };
	void setName(SRstring& nameIn);
	void assignTempRho(double trefIn, double ax, double ay, double az, double rhoIn, double allowable);
	void assignIso(double et, double nut);
	void assignOrtho(SRcij& cij);
	void assignGenAniso(SRgenAnisoCij& gcij);
	char* GetName(){ return name.str; };
	SRmaterialType GetType(){ return type; };
	double GetRho(){ return rho; };
	double MatScale();
	double GetAllowableStress(){ return allowableStress; };
	bool isAllowableAssigned(){ return allowableAssigned; };
	void GetIso(double& et, double& nut){ et = E; nut = nu; };
	double GetMaxSvm(){ return maxSvm; };
	void PutMaxSvm(double s);
	double GetMaxStress(){ return maxStress; };
	void PutMaxStress(double s){ maxStress = s; };
	double GetMaxStrain(){ return maxStrain; };
	void PutMaxStrain(double e){ maxStrain = e; };
	bool GetHighStressWarned() { return highStressWarned; };
	void SetHighStressWarned(bool tf = true) { highStressWarned = tf; };
	double GetVolPercentYielded();
	void SetVolPercentYielded(double p) { volPercentYielded = p; };
	void AddToVolPercentYielded(double p) { volPercentYielded += p; };
	bool isActive(){ return (numElements > 0); };
	bool diffElast(SRmaterial* that);
	int getId(){ return id; };
	SRgenAnisoCij& getGenAnisoCij() { return genAnisoCij; };
	SRcij& getOrthoCij() { return orthoCij; };
	double getAlphax(){ return alphax; };
	double getAlphay(){ return alphay; };
	double getAlphaz(){ return alphaz; };
	double getRho(){ return rho; };
	double getE(){ return E; };
	double getC11(){ return c11; };
	double getLambda(){ return lambda; };
	double getNu(){ return nu; };
	double getG(){ return G; };
	double getTref(){ return tref; };
	double getAllowableStress(){ return allowableStress; };
	void updateNumElements() {numElements++; };
	int getNumElements() {return numElements; };

private:
	int id;
	int uid;
	SRstring name;
	SRmaterialType type;
	//density:
	double rho;
	//coefficient of thermal expansion:
	double alphax;
	double alphay;
	double alphaz;
	//elastic properties:
	double E;
	double nu;
	double lambda;
	double G;
	double c11; //=lambda+2G for iso
	double tref;
	double allowableStress;
	SRcij orthoCij;
	SRgenAnisoCij genAnisoCij;
	bool allowableAssigned;
	bool highStressWarned;
	double maxSvm;
	double maxStress; //max vm or in any component
	double maxStrain;
	double svmavg;
	double volPercentYielded;
	int numElements;
	int numelyielded;
};

struct SRElProperty
{
	int uid;
	int matid;
	int matuid;
};


#endif //!defined(SRMATERIAL_INCLUDED)
