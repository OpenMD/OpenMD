/* This is a fun little patch to do molecular restraints for
   thermodynamic integration of solids.  Use only 
   if you really know what you are doing.
*/

#ifndef _RESTRAINTS_H_
#define _RESTRAINTS_H_

#include <stdlib.h>
#include <iostream>
#include <vector>

#include "primitives/Atom.hpp"
#include "brains/SimState.hpp"

//#include "brains/SimInfo.hpp"
//#include "io/ReadWrite.hpp"

class Restraints{

 public:
  Restraints(double lambdaVal, double lambdaExp);
  ~Restraints();

  void Calc_rVal(double position[3], double refPosition[3]);
  void Calc_body_thetaVal(double matrix[3][3], double refUnit[3]);
  void Calc_body_omegaVal(double matrix[3][3], double zAngle);
  double Calc_Restraint_Forces(vector<StuntDouble*> vecParticles);
  //  void Store_Init_Info(vector<StuntDouble*> vecParticles);
  void Write_zAngle_File(vector<StuntDouble*> vecParticles, 
			 int currTime,
			 int nIntObj);
  double getVharm();

 private:
  char moleculeName[15];

  int i, j;

  double scaleLam;
  double delRx, delRy, delRz;
  double theta, omega;
  double vProj0[3];
  double vProjDist;
  double uTx, uTy, uTz, vTx, vTy, vTz;
  double ub0x, ub0y, ub0z, vb0x, vb0y, vb0z;
  double kTheta, kOmega, kDist;
  double restraintFrc[3];
  double restraintTrq[3];
  double normalize;
  double dVdrx, dVdry, dVdrz;
  double dVdux, dVduy, dVduz;
  double dVdvx, dVdvy, dVdvz;
  double harmPotent;
  double lambdaValue;
  double lambdaK;

  char *token;
  char fileName[200];
  char angleName[200];
  char inLine[1000];
  char inValue[200];
  char springName[200];
};

#endif
