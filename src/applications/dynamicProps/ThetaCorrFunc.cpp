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

#include "applications/dynamicProps/ThetaCorrFunc.hpp"


namespace oopse {
  ThetaCorrFunc::ThetaCorrFunc(SimInfo* info, const std::string& filename, const std::string& sele1, const std::string& sele2)
    : ParticleTimeCorrFunc(info, filename, sele1, sele2, DataStorage::dslPosition){

      setCorrFuncType("ThetaCorrFunc");
      setOutputName(getPrefix(dumpFilename_) + ".tcorr");

// Turn on COM calculation in block snapshot
      bool ncp = true;
      bsMan_->needCOMprops(ncp);
    }

  RealType ThetaCorrFunc::calcCorrVal(int frame1, int frame2, StuntDouble* sd1, StuntDouble* sd2) {
    RealType r0;
    RealType rt;
    RealType csth;
    Vector3d coord_t0;
    Vector3d coord_t;

    Vector3d r1 = sd1->getPos(frame1);
    Vector3d r2 = sd2->getPos(frame2);
    Vector3d com1 = sd1->getCOM(frame1);
    Vector3d com2 = sd2->getCOM(frame2);

    coord_t0 = r1 - com1;
    coord_t  = r2 - com2;
    
    r0 = sqrt(dot(coord_t0, coord_t0));
    rt = sqrt(dot(coord_t,coord_t));

    csth = dot(coord_t0,coord_t)/(r0*rt);
    return csth;
  }

}
