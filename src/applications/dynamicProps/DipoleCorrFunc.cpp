/*
 * Copyright (c) 2005 The University of Notre Dame. All Rights Reserved.
 *
 * The University of Notre Dame grants you ("Licensee") a
 * non-exclusive, royalty free, license to use, modify and
 * redistribute this software in source and binary code form, provided
 * that the following conditions are met:
 *
 * 1. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright
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
 *
 * SUPPORT OPEN SCIENCE!  If you use OpenMD or its source code in your
 * research, please cite the appropriate papers when you publish your
 * work.  Good starting points are:
 *                                                                      
 * [1]  Meineke, et al., J. Comp. Chem. 26, 252-271 (2005).             
 * [2]  Fennell & Gezelter, J. Chem. Phys. 124, 234104 (2006).          
 * [3]  Sun, Lin & Gezelter, J. Chem. Phys. 128, 234107 (2008).          
 * [4]  Kuang & Gezelter,  J. Chem. Phys. 133, 164101 (2010).
 * [5]  Vardeman, Stocker & Gezelter, J. Chem. Theory Comput. 7, 834 (2011).
 */

#include "applications/dynamicProps/DipoleCorrFunc.hpp"
#include "primitives/Atom.hpp"
#include "types/MultipoleAdapter.hpp"
#include "utils/simError.h"

namespace OpenMD {
  DipoleCorrFunc::DipoleCorrFunc(SimInfo* info, const std::string& filename, const std::string& sele1, const std::string& sele2, long long int memSize)
    : ParticleTimeCorrFunc(info, filename, sele1, sele2, DataStorage::dslAmat | DataStorage::dslDipole, memSize){

      setCorrFuncType("Dipole Correlation Function");
      setOutputName(getPrefix(dumpFilename_) + ".dcorr");

    }

  RealType DipoleCorrFunc::calcCorrVal(int frame1, int frame2, StuntDouble* sd1,  StuntDouble* sd2) {

    Vector3d v1 = sd1->getDipole(frame1);
    Vector3d v2 = sd2->getDipole(frame2);

    return dot(v1, v2)/(v1.length()*v2.length());
  }


  void DipoleCorrFunc::validateSelection(const SelectionManager& seleMan) {
    StuntDouble* sd;
    int i;    
    for (sd = seleMan1_.beginSelected(i); sd != NULL; 
         sd = seleMan1_.nextSelected(i)) {

      AtomType* at = static_cast<Atom*>(sd)->getAtomType();
      MultipoleAdapter ma = MultipoleAdapter(at);

      if (!ma.isDipole()) {
	sprintf(painCave.errMsg,
                "DipoleCorrFunc::validateSelection Error: selected atoms do\n"
                "\tnot have a dipole\n");
	painCave.isFatal = 1;
	simError();        
      }
    }    
  }
}


