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

#include "applications/dynamicProps/ActionCorrFunc.hpp"
#include "utils/Constants.hpp"
#include "utils/Revision.hpp"

namespace OpenMD {

  // We need all of the positions, velocities, etc. so that we can
  // recalculate pressures and actions on the fly:
  ActionCorrFunc::ActionCorrFunc(SimInfo* info, const std::string& filename, 
				 const std::string& sele1, 
				 const std::string& sele2)
    : MPFrameTimeCorrFunc<Mat3x3d>(info, filename, sele1, sele2, 
                                   DataStorage::dslPosition | 
                                   DataStorage::dslVelocity |
                                   DataStorage::dslForce ){
    
      setCorrFuncType("ActionCorrFunc");      
      setOutputName(getPrefix(dumpFilename_) + ".action");

      // We'll need the force manager to compute forces for the average pressure
      forceMan_ = new ForceManager(info);
      
      // We'll need thermo to compute the pressures from the virial
      thermo_ =  new Thermo(info);

      pSum_ = 0.0;
      vSum_ = 0.0;
      nsamp_ = 0;
    }


  void ActionCorrFunc::computeProperty(int frame1) {

    StuntDouble* sd1;
    int i;
    
    // do the forces:
    forceMan_->calcForces();
    // call thermo to get the pressure and volume.
    pSum_ += thermo_->getPressure();
    RealType vol1 =  thermo_->getVolume();
    vSum_ += vol1;
    nsamp_++;

    Mat3x3d actionTensor1(0.0);
    for (sd1 = seleMan1_.beginSelected(i); sd1 != NULL;
         sd1 = seleMan1_.nextSelected(i)) {

      Vector3d r1 = sd1->getPos(frame1);
      Vector3d v1 = sd1->getVel(frame1);
      RealType m = sd1->getMass();

      actionTensor1 += m*outProduct(r1, v1);
    }
    actionTensor_[frame1] = actionTensor1 / vol1;
    
  }
  
  void ActionCorrFunc::correlateFrames(int frame1, int frame2, int timeBin) {

    Mat3x3d corrTensor(0.0);
    RealType thisTerm;

    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 3; j++) {

        thisTerm = actionTensor_[frame2](i, j) - actionTensor_[frame1](i, j);
        
        if (i == j)
          thisTerm -= avePress_ * (times_[frame2] - times_[frame1]);
 
        corrTensor(i, j) += thisTerm * thisTerm;
      }
    }

    histogram_[timeBin] += corrTensor;    
    count_[timeBin]++;    
  }

  void ActionCorrFunc::preCorrelate() {

    this->preCorrelate();

    avePress_ = pSum_ / ( Constants::pressureConvert * (RealType)nsamp_);
    aveVol_ = vSum_ / (RealType)nsamp_;
  } 
}
