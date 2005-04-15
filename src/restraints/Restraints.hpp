/*
 * Copyright (c) 2005 The University of Notre Dame. All Rights Reserved.
 *
 * The University of Notre Dame grants you ("Licensee") a
 * non-exclusive, royalty free, license to use, modify and
 * redistribute this software in source and binary code form, provided
 * that the following conditions are met:
 *
 * 1. Acknowledgement of the program authors must be made in any
 *    publication of scientific results based in part on use of the
 *    program.  An acceptable form of acknowledgement is citation of
 *    the article in which the program was described (Matthew
 *    A. Meineke, Charles F. Vardeman II, Teng Lin, Christopher
 *    J. Fennell and J. Daniel Gezelter, "OOPSE: An Object-Oriented
 *    Parallel Simulation Engine for Molecular Dynamics,"
 *    J. Comput. Chem. 26, pp. 252-271 (2005))
 *
 * 2. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 *
 * 3. Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in the
 *    documentation and/or other materials provided with the
 *    distribution.
 *
 * This software is provided "AS IS," without a warranty of any
 * kind. All express or implied conditions, representations and
 * warranties, including any implied warranty of merchantability,
 * fitness for a particular purpose or non-infringement, are hereby
 * excluded.  The University of Notre Dame and its licensors shall not
 * be liable for any damages suffered by licensee as a result of
 * using, modifying or distributing the software or its
 * derivatives. In no event will the University of Notre Dame or its
 * licensors be liable for any lost revenue, profit or data, or for
 * direct, indirect, special, consequential, incidental or punitive
 * damages, however caused and regardless of the theory of liability,
 * arising out of the use of or inability to use software, even if the
 * University of Notre Dame has been advised of the possibility of
 * such damages.
 */

#ifndef _RESTRAINTS_H_
#define _RESTRAINTS_H_

#include <stdlib.h>
#include <vector>
#include "primitives/Atom.hpp"
#include "brains/SimInfo.hpp"
#include "io/RestReader.hpp"
#include "math/SquareMatrix3.hpp"
#include "math/Vector3.hpp"

namespace oopse {
  
  class Restraints{
    
  public:
    Restraints(SimInfo * info, double lambdaVal, double lambdaExp);
    ~Restraints();
    
    void Calc_rVal(Vector3d &position, double refPosition[3]);
    void Calc_body_thetaVal(RotMat3x3d &matrix, double refUnit[3]);
    void Calc_body_omegaVal(double zAngle);
    double Calc_Restraint_Forces();
    double getVharm() { return harmPotent; }
    
  private:
    SimInfo * info_;
    RestReader* restRead_;
    
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

} // end namespace oopse

#endif
