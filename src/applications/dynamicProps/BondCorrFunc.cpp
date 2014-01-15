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

#include "applications/dynamicProps/BondCorrFunc.hpp"

namespace OpenMD {
  BondCorrFunc::BondCorrFunc(SimInfo* info, const std::string& filename, 
                             const std::string& sele1, long long int memSize)
    : InteractionTimeCorrFunc(info, filename, sele1, 
                              DataStorage::dslPosition, memSize){

      setCorrFuncType("Bond Correlation Function");
      setOutputName(getPrefix(dumpFilename_) + ".bondcorr");
      nSelectedBonds_ =   seleMan1_.getBondSelectionCount();  
  }
  

  void BondCorrFunc::correlateFrames(int frame1, int frame2) {
    snapshot1 = bsMan_->getSnapshot(frame1);
    snapshot2 = bsMan_->getSnapshot(frame2);
    assert(snapshot1 && snapshot2);

    RealType time1 = snapshot1->getTime();
    RealType time2 = snapshot2->getTime();

    int timeBin = int ((time2 - time1) /deltaTime_ + 0.5);
    count_[timeBin] += nSelectedBonds_;    

    int i;

    for (bond = seleMan1_.beginSelectedBond(i);
         bond != NULL ;
         bond = seleMan1_.nextSelectedBond(i)) {

      RealType corrVal = calcCorrVal(frame1, frame2, bond);
      histogram_[timeBin] += corrVal;    
    }    
  }

  RealType BondCorrFunc::calcCorrVal(int frame1, int frame2, ShortRangeInteraction* sri) {
    bond = dynamic_cast<Bond*>(sri);
    re = bond->getBondType()->getEquilibriumBondLength();
    val1 = bond->getValue(frame1) - re;
    val2 = bond->getValue(frame2) - re;
    return val1 * val2;
  }
}

