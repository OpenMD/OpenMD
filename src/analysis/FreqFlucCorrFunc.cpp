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

#include "analysis/FreqFlucCorrFunc.hpp"
#include "primitives/Atom.hpp"
#include "types/MultipoleAdapter.hpp"
#include "utils/simError.h"

namespace OpenMD {
  FreqFlucCorrFunc::FreqFlucCorrFunc(SimInfo* info, const std::string& filename,
                                     const std::string& sele1, 
                                     const std::string& sele2, 
                                     long long int memSize)
    : ParticleTimeCorrFunc(info, filename, sele1, sele2,  
                           DataStorage::dslElectricField | 
                           DataStorage::dslAmat | DataStorage::dslDipole, 
                           memSize){
    
    setCorrFuncType("FreqFluc Correlation Function");
    setOutputName(getPrefix(dumpFilename_) + ".ffcorr");
    
  }

  RealType FreqFlucCorrFunc::calcCorrVal(int frame1, int frame2, 
                                         StuntDouble* sd1,  StuntDouble* sd2) {

    Vector3d e1;    
    Vector3d u1;
    Vector3d e2;    
    Vector3d u2;

    e1 = sd1->getElectricField(frame1);
    
    if (sd1->isRigidBody()) {
      u1 = sd1->getA(frame1).getRow(2);
    } else {
      AtomType* at = static_cast<Atom*>(sd1)->getAtomType();
      MultipoleAdapter ma = MultipoleAdapter(at);
      
      if (ma.isDipole()) {
        u1 = sd1->getDipole(frame1);
      }
    }

    e2 = sd2->getElectricField(frame2);
    
    if (sd2->isRigidBody()) {
      u2 = sd2->getA(frame2).getRow(2);
    } else {
      AtomType* at = static_cast<Atom*>(sd2)->getAtomType();
      MultipoleAdapter ma = MultipoleAdapter(at);
      
      if (ma.isDipole()) {
        u2 = sd2->getDipole(frame2);
      }
    }
    
    return (dot(u1, e1) - mean_) * (dot(u2, e2) - mean_);
  }


  void FreqFlucCorrFunc::validateSelection(const SelectionManager& seleMan) {
    StuntDouble* sd;
    int i;    
    for (sd = seleMan1_.beginSelected(i); sd != NULL; 
         sd = seleMan1_.nextSelected(i)) {

      if (!sd->isRigidBody()) {

        AtomType* at = static_cast<Atom*>(sd)->getAtomType();
        MultipoleAdapter ma = MultipoleAdapter(at);
        
        if (!ma.isDipole()) {
          sprintf(painCave.errMsg,
                  "FreqFlucCorrFunc::validateSelection Error: selection is not a RigidBody or does\n"
                  "\tnot have a dipole\n");
          painCave.isFatal = 1;
          simError();        
        }
      }
    }    
  }

  void FreqFlucCorrFunc::preCorrelate() {
    mean_ = 0.0;
    RealType sum(0.0);
    int count(0);
    StuntDouble* sd1;
    Vector3d e1;
    Vector3d u1;
    int ii;
    std::cerr << "preCorrelating to compute mean values\n";
    // dump files can be enormous, so read them in block-by-block:
    int nblocks = bsMan_->getNBlocks();

    for (int i = 0; i < nblocks; ++i) {
      bsMan_->loadBlock(i);
      assert(bsMan_->isBlockActive(i));      
      SnapshotBlock block1 = bsMan_->getSnapshotBlock(i);
      for (int j = block1.first; j < block1.second; ++j) {

        for (sd1 = seleMan1_.beginSelected(ii); sd1 != NULL; 
	     sd1 = seleMan1_.nextSelected(ii)) {

          e1 = sd1->getElectricField(j);
          
          if (sd1->isRigidBody()) {
            u1 = sd1->getA(j).getRow(2);
          } else {
            AtomType* at = static_cast<Atom*>(sd1)->getAtomType();
            MultipoleAdapter ma = MultipoleAdapter(at);
            
            if (ma.isDipole()) {
              u1 = sd1->getDipole(j);
            }
          }

          RealType value = dot(u1, e1);

          sum += value;                         
          count++;
        }
      }
      bsMan_->unloadBlock(i);
    }

    mean_ = sum / RealType(count);
    std::cerr << "done with preCorrelation\n";
  }
}


