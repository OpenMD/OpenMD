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
 * [3]  Sun, Lin & Gezelter, J. Chem. Phys. 128, 24107 (2008).          
 * [4]  Vardeman & Gezelter, in progress (2009).                        
 */

#include "applications/dynamicProps/StressCorrFunc.hpp"
#include "utils/PhysicalConstants.hpp"
#include "brains/ForceManager.hpp"
#include "brains/Thermo.hpp"

namespace OpenMD {

  // We need all of the positions, velocities, etc. so that we can
  // recalculate pressures and actions on the fly:
  StressCorrFunc::StressCorrFunc(SimInfo* info, const std::string& filename, 
				 const std::string& sele1, 
				 const std::string& sele2)
    : FrameTimeCorrFunc(info, filename, sele1, sele2, 
			DataStorage::dslPosition | 
			DataStorage::dslVelocity |
			DataStorage::dslForce ){

      setCorrFuncType("StressCorrFunc");
      setOutputName(getPrefix(dumpFilename_) + ".action");
      histogram_.resize(nTimeBins_); 
      count_.resize(nTimeBins_);
    }

  void StressCorrFunc::correlateFrames(int frame1, int frame2) {
    Snapshot* snapshot1 = bsMan_->getSnapshot(frame1);
    Snapshot* snapshot2 = bsMan_->getSnapshot(frame2);
    assert(snapshot1 && snapshot2);

    RealType time1 = snapshot1->getTime();
    RealType time2 = snapshot2->getTime();
    RealType vol1 = snapshot1->getVolume();
    RealType vol2 = snapshot2->getVolume();
       
    int timeBin = int ((time2 - time1) /deltaTime_ + 0.5);

    //std::cerr << "times = " << time1 << " " << time2 << "\n";
    //std::cerr << "vols = " << vol1 << " " << vol2 << "\n";

    int i;
    int j;

    StuntDouble* sd1;

    Mat3x3d actionTensor1(0.0);
    //std::cerr << "at1 = " << actionTensor1 << "\n";
    Mat3x3d actionTensor2(0.0);
    //std::cerr << "at2 = " << actionTensor2 << "\n";

    for (sd1 = seleMan1_.beginSelected(i); sd1 != NULL;
         sd1 = seleMan1_.nextSelected(i)) {
      //std::cerr << "found a SD\n";
      Vector3d r1 = sd1->getPos(frame1);
      //std::cerr << "r1 = " << r1 << "\n";
      Vector3d v1 = sd1->getVel(frame1);
      //std::cerr << "v1 = " << v1 << "\n";
      Vector3d r2 = sd1->getPos(frame2);
      //std::cerr << "r2 = " << r2 << "\n";
      Vector3d v2 = sd1->getVel(frame2);
      //std::cerr << "v2 = " << v2 << "\n";

      RealType m = sd1->getMass();

      //std::cerr << "m = " << m << "\n";

      actionTensor1 += m*outProduct(r1, v1);
      actionTensor2 += m*outProduct(r2, v2);
    }

    actionTensor1 /= vol1;
    //std::cerr << "at1 = " << actionTensor1 << "\n";
    actionTensor2 /= vol2;
    //std::cerr << "at2 = " << actionTensor2 << "\n";

    Mat3x3d corrTensor(0.0);
    //std::cerr << "ct = " << corrTensor << "\n";
    RealType thisTerm;

    for (i = 0; i < 3; i++) {
      for (j = 0; j < 3; j++) {
       
        //std::cerr << "i, j = " << i << " " << j << "\n"; 
        if (i == j) {
          thisTerm = (actionTensor2(i, j) - actionTensor1(i, j) - avePress_ *(time2-time1));
          std::cerr << "at1, at2 = " << actionTensor1(i,j) << " " << actionTensor2(i,j) << " p = " << avePress_ << "\n";
        } else {
          thisTerm = (actionTensor2(i, j) - actionTensor1(i, j));
        }
        
        //std::cerr << "thisTerm = " << thisTerm << "\n"; 
        corrTensor(i, j) += thisTerm * thisTerm;
      }
    }

    //std::cerr << "ct = " << corrTensor << "\n";
    //std::cerr << "hist = " << histogram_[timeBin] << "\n";
    histogram_[timeBin] += corrTensor;    
    count_[timeBin]++;    
    
  }

  void StressCorrFunc::postCorrelate() {
    for (int i =0 ; i < nTimeBins_; ++i) {
      if (count_[i] > 0) {
        histogram_[i] /= count_[i];
      }
    }
  }


  void StressCorrFunc::preCorrelate() {
    // Fill the histogram with empty 3x3 matrices:
    std::fill(histogram_.begin(), histogram_.end(), Mat3x3d(0.0));
    // count array set to zero
    std::fill(count_.begin(), count_.end(), 0);

    SimInfo::MoleculeIterator mi;
    Molecule* mol;
    Molecule::AtomIterator ai;
    Atom* atom;

    // We'll need the force manager to compute forces for the average pressure
    ForceManager* forceMan = new ForceManager(info_);
    forceMan->init();

    // We'll need thermo to compute the pressures from the virial
    Thermo* thermo =  new Thermo(info_);

    // prepare the averages
    RealType pSum = 0.0;
    RealType vSum = 0.0;
    int nsamp = 0;

    // dump files can be enormous, so read them in block-by-block:
    int nblocks = bsMan_->getNBlocks();
    for (int i = 0; i < nblocks; ++i) {
      std::cerr << "block = " << i << "\n";
      bsMan_->loadBlock(i);
      assert(bsMan_->isBlockActive(i));      
      SnapshotBlock block1 = bsMan_->getSnapshotBlock(i);
      for (int j = block1.first; j < block1.second; ++j) {

	// go snapshot-by-snapshot through this block:
        Snapshot* snap = bsMan_->getSnapshot(j);
       
        // update the positions and velocities of the atoms belonging
        // to rigid bodies:

        updateFrame(j);

	// do the forces:
	//forceMan->calcForces(true, true);
	// call thermo to get the pressure and volume.
        pSum += thermo->getPressure();
        vSum += thermo->getVolume();
        nsamp++;
        
      }
      bsMan_->unloadBlock(i);
    }

    avePress_ = pSum / ( PhysicalConstants::pressureConvert * (RealType)nsamp);
    aveVol_ = vSum / (RealType)nsamp;
    std::cout << "pAve = " << avePress_ << " vAve = " << aveVol_ << "\n";
  }   


  void StressCorrFunc::writeCorrelate() {
    std::ofstream ofs(getOutputFileName().c_str());

    if (ofs.is_open()) {

      ofs << "#" << getCorrFuncType() << "\n";
      ofs << "#time\tcorrTensor\txx\txy\txz\tyx\tyy\tyz\tzx\tzy\tzz\n";

      for (int i = 0; i < nTimeBins_; ++i) {
        ofs << time_[i] << "\t" << 
	  histogram_[i](0,0) << "\t" <<
	  histogram_[i](0,1) << "\t" <<
	  histogram_[i](0,2) << "\t" <<
	  histogram_[i](1,0) << "\t" <<
	  histogram_[i](1,1) << "\t" <<
	  histogram_[i](1,2) << "\t" <<
	  histogram_[i](2,0) << "\t" <<
	  histogram_[i](2,1) << "\t" <<
	  histogram_[i](2,2) << "\t" << "\n";
      }
            
    } else {
      sprintf(painCave.errMsg,
              "StressCorrFunc::writeCorrelate Error: fail to open %s\n", getOutputFileName().c_str());
      painCave.isFatal = 1;
      simError();        
    }

    ofs.close();    
  }

}
