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
 
#include "applications/dynamicProps/CrossTimeCorrFunc.hpp"

namespace oopse {

  CrossTimeCorrFunc::CrossTimeCorrFunc(SimInfo * info, const  std::string & filename,
				       const std :: string & sele1, const std :: string & sele2, int storageLayout) 
    : TimeCorrFunc(info, filename, sele1, sele2, storageLayout) {
    
      nSelected1_ = seleMan1_.getSelectionCount();  
      nSelected2_ = seleMan2_.getSelectionCount();  
      nSelectedPairs_ = nSelected1_ * nSelected2_;
    }

  void CrossTimeCorrFunc::correlateFrames(int frame1, int frame2) {
    Snapshot* snapshot1 = bsMan_->getSnapshot(frame1);
    Snapshot* snapshot2 = bsMan_->getSnapshot(frame2);
    assert(snapshot1 && snapshot2);

    double time1 = snapshot1->getTime();
    double time2 = snapshot2->getTime();

    int timeBin = int ((time2 - time1) /deltaTime_ + 0.5);
    count_[timeBin] += nSelectedPairs_ ;    

    int i;
    int j;
    StuntDouble* sd1;
    StuntDouble* sd2;
    for (sd1 = seleMan1_.beginSelected(i); sd1 != NULL; sd1 = seleMan1_.nextSelected(i)) {

      for (sd2 = seleMan2_.beginSelected(j); sd2 != NULL; sd2 = seleMan2_.nextSelected(j)) {
	double corrVal = calcCorrVal(frame1, frame2, sd1, sd2);
	histogram_[timeBin] += corrVal;    
      }            
    }
  }

}
