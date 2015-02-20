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

#include "applications/dynamicProps/MultipassCorrFunc.hpp"
#include "utils/simError.h"
#include "primitives/Molecule.hpp"
using namespace std;
namespace OpenMD {

  MultipassCorrFunc::MultipassCorrFunc(SimInfo* info, const string& filename, 
                                       const string& sele1, const string& sele2,
                                       int storageLayout)
    : info_(info), storageLayout_(storageLayout),
      dumpFilename_(filename), selectionScript1_(sele1), 
      selectionScript2_(sele2), evaluator1_(info), evaluator2_(info), 
      seleMan1_(info), seleMan2_(info) {
    
    int nAtoms = info->getNGlobalAtoms();
    int nRigidBodies = info->getNGlobalRigidBodies();

    // Request maximum needed storage for the simulation (including of
    // whatever was passed down by the individual correlation
    // function).

    storageLayout_ = info->getStorageLayout() | storageLayout;

    DumpReader reader(info_, dumpFilename_);    
    
    evaluator1_.loadScriptString(selectionScript1_);
    evaluator2_.loadScriptString(selectionScript2_);
    
    //if selection is static, we only need to evaluate it once
    if (!evaluator1_.isDynamic()) {
      seleMan1_.setSelectionSet(evaluator1_.evaluate());
      validateSelection(seleMan1_);
    }
    
    if (!evaluator2_.isDynamic()) {
      seleMan2_.setSelectionSet(evaluator2_.evaluate());
      validateSelection(seleMan2_);
    }
    
    /**@todo Fix Me */
    Globals* simParams = info_->getSimParams();
    if (simParams->haveSampleTime()){
      deltaTime_ = simParams->getSampleTime();
    } else {
      sprintf(painCave.errMsg,
              "MultipassCorrFunc::writeCorrelate Error: can not figure out deltaTime\n");
      painCave.isFatal = 1;
      simError();  
    }

  }

  void MultipassCorrFunc::preCorrelate() {
    Molecule* mol;
    RigidBody* rb;
    SimInfo::MoleculeIterator mi;
    Molecule::RigidBodyIterator rbIter;
    StuntDouble* sd;

    int index, isd1, isd2;
    
    fill(histogram_.begin(), histogram_.end(), 0.0);
    fill(count_.begin(), count_.end(), 0);

    DumpReader reader(info_, dumpFilename_);    

    nFrames_ = reader.getNFrames();
    times_.resize(nFrames_);
    
    for (int istep = 0; istep < nFrames_; istep++) {
      reader.readFrame(istep);
      currentSnapshot_ = info_->getSnapshotManager()->getCurrentSnapshot();
      times_[istep] = currentSnapshot_->getTime();

      // update the positions of atoms which belong to the rigidbodies
      for (mol = info_->beginMolecule(mi); mol != NULL;
           mol = info_->nextMolecule(mi)) {
        for (rb = mol->beginRigidBody(rbIter); rb != NULL;
             rb = mol->nextRigidBody(rbIter)) {
          rb->updateAtoms();
        }
      }
      
      if (evaluator1_.isDynamic()) {
        seleMan1_.setSelectionSet(evaluator1_.evaluate());
      }
      
      if (evaluator2_.isDynamic()) {
        seleMan2_.setSelectionSet(evaluator2_.evaluate());
      }      

      for (sd = seleMan1_.beginSelected(isd1); sd != NULL;
           sd = seleMan1_.nextSelected(isd1)) {

        sele1ToIndex_[istep].push_back(sd->getGlobalIndex());
        index = sele1ToIndex_[istep].size();
        computeProperty1(istep, sd, index);
      }
      
      for (sd = seleMan2_.beginSelected(isd2); sd != NULL;
           sd = seleMan2_.nextSelected(isd2)) {

        sele2ToIndex_[istep].push_back(sd->getGlobalIndex());
        index = sele2ToIndex_[istep].size();
                
        computeProperty2(istep, sd, index);
      }
    }     
  }

  
  void MultipassCorrFunc::doCorrelate() {
    preCorrelate();
    correlation();   
    postCorrelate();
    writeCorrelate();
  }

  void MultipassCorrFunc::correlation() {
    std::vector<int> s1;
    std::vector<int> s2;

    std::vector<int>::iterator i1;
    std::vector<int>::iterator i2;

    
    for (int i = 0; i < nFrames_; ++i) {

      RealType time1 = times_[i];
      s1 = sele1ToIndex_[i];
      
      for(int j  = i; j < nFrames_; ++j) {

        // Perform a sanity check on the actual configuration times to
        // make sure the configurations are spaced the same amount the
        // sample time said they were spaced:
               
        RealType time2 = times_[j];

        if ( fabs( (time2 - time1) - (j-i)*deltaTime_ ) > 1.0e-4 ) {
          sprintf(painCave.errMsg,
                  "MultipassCorrFunc::correlateBlocks Error: sampleTime (%f)\n"
                  "\tin %s does not match actual time-spacing between\n"
                  "\tconfigurations %d (t = %f) and %d (t = %f).\n",
                  deltaTime_, dumpFilename_.c_str(), i, time1, j, time2);
          painCave.isFatal = 1;
          simError();  
        }
        
        int timeBin = int ((time2 - time1) / deltaTime_ + 0.5);        
        s2 = sele2ToIndex_[j];

        for (i1 = s1.begin(), i2 = s2.begin();
             i1 != s1.end() && i2 != s2.end(); ++i1, ++i2){
          
          // If the selections are dynamic, they might not have the
          // same objects in both frames, so we need to roll either of
          // the selections until we have the same object to
          // correlate.
          
          while ( *i1 < *i2 && i1 != s1.end()) {
            ++i1;
          }
          
          while ( *i2 < *i1 && i2 != s2.end() ) {
            ++i2;
          }
          
          if ( i1 == s1.end() || i2 == s2.end() ) break;

          RealType corrVal = calcCorrVal(i, j, *i1, *i2);
          histogram_[timeBin] += corrVal;
          count_[timeBin]++;
          
        }
      }
    }
  }


  
  void MultipassCorrFunc::postCorrelate() {
    for (int i =0 ; i < nTimeBins_; ++i) {
      if (count_[i] > 0) {
	histogram_[i] /= count_[i];
      } else {
        histogram_[i] = 0;
      }
    }
  }


  void MultipassCorrFunc::updateFrame(int frame){
    Molecule* mol;
    RigidBody* rb;
    SimInfo::MoleculeIterator mi;
    Molecule::RigidBodyIterator rbIter;
    /** @todo need improvement */    
    if (storageLayout_ & DataStorage::dslPosition) {
      for (mol = info_->beginMolecule(mi); mol != NULL; 
           mol = info_->nextMolecule(mi)) {

        //change the positions of atoms which belong to the rigidbodies
        for (rb = mol->beginRigidBody(rbIter); rb != NULL; 
             rb = mol->nextRigidBody(rbIter)) {
          rb->updateAtoms(frame);
        }
      }        
    }

    if (storageLayout_ & DataStorage::dslVelocity) {
      for (mol = info_->beginMolecule(mi); mol != NULL; 
           mol = info_->nextMolecule(mi)) {
        
        //change the positions of atoms which belong to the rigidbodies
        for (rb = mol->beginRigidBody(rbIter); rb != NULL; 
             rb = mol->nextRigidBody(rbIter)) {
          rb->updateAtomVel(frame);
        }
      }      
    }    
  }



  void MultipassCorrFunc::writeCorrelate() {
    ofstream ofs(outputFilename_.c_str());

    if (ofs.is_open()) {

      ofs << "#" << getCorrFuncType() << "\n";
      ofs << "#selection script1: \"" << selectionScript1_ ;
      ofs << "\"\tselection script2: \"" << selectionScript2_ << "\"\n";
      ofs << "#extra information: " << extra_ << "\n";
      ofs << "#time\tcorrVal\n";

      for (int i = 0; i < nTimeBins_; ++i) {
	ofs << times_[i]-times_[0] << "\t" << histogram_[i] << "\n";
      }
            
    } else {
      sprintf(painCave.errMsg,
	      "MultipassCorrFunc::writeCorrelate Error: fail to open %s\n", 
              outputFilename_.c_str());
      painCave.isFatal = 1;
      simError();        
    }

    ofs.close();    
  }

}
