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

#include "applications/dynamicProps/HBondJump.hpp"
#include <algorithm>

namespace OpenMD {
  HBondJump::HBondJump(SimInfo* info, const std::string& filename,
                       const std::string& sele1, const std::string& sele2,
                       double rCut, double thetaCut)
    : MultipassCorrFunc(info, filename, sele1, sele2,
                        DataStorage::dslPosition | DataStorage::dslAmat ),
      rCut_(rCut), thetaCut_(thetaCut) {
    
    setCorrFuncType("HBondJump");
    setOutputName(getPrefix(dumpFilename_) + ".jump");

    //nFrames_ is initialized in MultipassCorrFunc:
    bondList_.resize(nFrames_);
  }
  
  void HBondJump::preCorrelate() {
    Molecule* mol1;
    Molecule* mol2;
    RigidBody* rb1;
    Molecule::HBondDonor* hbd1;
    Molecule::HBondDonor* hbd2;
    std::vector<Molecule::HBondDonor*>::iterator hbdi;
    std::vector<Molecule::HBondDonor*>::iterator hbdj;
    std::vector<Atom*>::iterator hbai;
    std::vector<Atom*>::iterator hbaj;
    Atom* hba1;
    Atom* hba2;
    SimInfo::MoleculeIterator mi;
    Molecule::RigidBodyIterator rbIter;
    Vector3d dPos;
    Vector3d aPos;
    Vector3d hPos;
    Vector3d DH;
    Vector3d DA;
    RealType DAdist, DHdist, theta, ctheta;
    int ii, jj;
    int hInd, aInd, index;
    std::vector<int>::iterator ind;
   
    for (int istep = 0; istep < nFrames_; istep++) {
      reader_->readFrame(istep);
      currentSnapshot_ = info_->getSnapshotManager()->getCurrentSnapshot();
      times_[istep] = currentSnapshot_->getTime();

      // update the positions of atoms which belong to the rigidbodies
      for (mol1 = info_->beginMolecule(mi); mol1 != NULL;
           mol1 = info_->nextMolecule(mi)) {
        for (rb1 = mol1->beginRigidBody(rbIter); rb1 != NULL;
             rb1 = mol1->nextRigidBody(rbIter)) {
          rb1->updateAtoms();
        }
      }
      
      if (evaluator1_.isDynamic()) {
        seleMan1_.setSelectionSet(evaluator1_.evaluate());
      }
      
      if (evaluator2_.isDynamic()) {
        seleMan2_.setSelectionSet(evaluator2_.evaluate());
      }      

      bondList_[istep].resize(seleMan2_.getSelectionCount());
      
      for (mol1 = seleMan1_.beginSelectedMolecule(ii);
           mol1 != NULL; mol1 = seleMan1_.nextSelectedMolecule(ii)) {
        
        for (mol2 = seleMan2_.beginSelectedMolecule(jj);
             mol2 != NULL; mol2 = seleMan2_.nextSelectedMolecule(jj)) {
          
          // loop over the possible donors in molecule 1:
          for (hbd1 = mol1->beginHBondDonor(hbdi); hbd1 != NULL;
               hbd1 = mol1->nextHBondDonor(hbdi)) {
            dPos = hbd1->donorAtom->getPos(); 
            hPos = hbd1->donatedHydrogen->getPos();
            DH = hPos - dPos; 
            currentSnapshot_->wrapVector(DH);
            DHdist = DH.length();

            // loop over the possible acceptors in molecule 2:
            for (hba2 = mol2->beginHBondAcceptor(hbaj); hba2 != NULL;
                 hba2 = mol2->nextHBondAcceptor(hbaj)) {
              aPos = hba2->getPos();
              DA = aPos - dPos;              
              currentSnapshot_->wrapVector(DA);
              DAdist = DA.length();

              // Distance criteria: are the donor and acceptor atoms
              // close enough?
              if (DAdist < rCut_) {

                ctheta = dot(DH, DA) / (DHdist * DAdist);
                theta = acos(ctheta) * 180.0 / M_PI;

                // Angle criteria: are the D-H and D-A and vectors close?
                if (theta < thetaCut_) {
                  // molecule 1 is a Hbond donor:
                  hInd = hbd1->donatedHydrogen->getGlobalIndex();

                  ind = std::find(indices_.begin(), indices_.end(), hInd);
                  if (ind == indices_.end()) {
                    index = indices_.size();
                    indices_.push_back(hInd);
                  } else {
                    index = std::distance(indices_.begin(), ind);
                  }

                  aInd = hba2->getGlobalIndex();
                  bondList_[istep][index].insert(aInd);                  
                }
              }            
            }            
          }

          // now loop over the possible acceptors in molecule 1:
          for (hba1 = mol1->beginHBondAcceptor(hbai); hba1 != NULL;
               hba1 = mol1->nextHBondAcceptor(hbai)) {
            aPos = hba1->getPos();
            
            // loop over the possible donors in molecule 2:
            for (hbd2 = mol2->beginHBondDonor(hbdj); hbd2 != NULL;
               hbd2 = mol2->nextHBondDonor(hbdj)) {
              dPos = hbd2->donorAtom->getPos();

              DA = aPos - dPos;
              currentSnapshot_->wrapVector(DA);
              DAdist = DA.length();
              
              // Distance criteria: are the donor and acceptor atoms
              // close enough?
              if (DAdist < rCut_) {
                hPos = hbd2->donatedHydrogen->getPos();
                DH = hPos - dPos; 
                currentSnapshot_->wrapVector(DH);
                DHdist = DH.length();
                ctheta = dot(DH, DA) / (DHdist * DAdist);
                theta = acos(ctheta) * 180.0 / M_PI;
                // Angle criteria: are the D-H and D-A and vectors close?
                if (theta < thetaCut_) {
                  // molecule 1 is a Hbond acceptor:
                  hInd = hbd2->donatedHydrogen->getGlobalIndex();
                  aInd = hba1->getGlobalIndex();
                  ind = std::find(indices_.begin(), indices_.end(), hInd);
                  if (ind == indices_.end()) {
                    index = indices_.size();
                    indices_.push_back(hInd);
                  } else {
                    index = std::distance(indices_.begin(), ind);
                  }
                  bondList_[istep][index].insert(aInd);
                }                
              }
            }
          }
        }
      }
    }   
  }
  
  
  void HBondJump::correlation() {
    std::vector<int>::iterator i1;
    std::set<int>::iterator ki;
    std::set<int>::iterator kj;
          
    RealType corrVal;
    
    for (int i = 0; i < nFrames_; ++i) {

      RealType time1 = times_[i];           
      
      for(int j  = i; j < nFrames_; ++j) {

        // Perform a sanity check on the actual configuration times to
        // make sure the configurations are spaced the same amount the
        // sample time said they were spaced:
               
        RealType time2 = times_[j];

        if ( fabs( (time2 - time1) - (j-i)*deltaTime_ ) > 1.0e-4 ) {
          sprintf(painCave.errMsg,
                  "HBondJump::correlation Error: sampleTime (%f)\n"
                  "\tin %s does not match actual time-spacing between\n"
                  "\tconfigurations %d (t = %f) and %d (t = %f).\n",
                  deltaTime_, dumpFilename_.c_str(), i, time1, j, time2);
          painCave.isFatal = 1;
          simError();  
        }
        
        int timeBin = int ((time2 - time1) / deltaTime_ + 0.5);        

        for (i1 = indices_.begin(); i1 != indices_.end(); ++i1) {

          corrVal = 0.0;
          
          for (ki = bondList_[i][*i1].begin();
               ki != bondList_[i][*i1].end(); ++ki) {
            
            for (kj = bondList_[j][*i1].begin();
                 kj != bondList_[j][*i1].end(); ++kj) {
              
              if ( *ki == *kj ) corrVal += 1;
            }
          }
          
          histogram_[timeBin] += corrVal;
          count_[timeBin] += bondList_[i][*i1].size();
        }
      }
    }
  }
}
