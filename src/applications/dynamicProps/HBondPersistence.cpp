/*
 * Copyright (c) 2004-2020 The University of Notre Dame. All Rights Reserved.
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
 * [1] Meineke, et al., J. Comp. Chem. 26, 252-271 (2005).
 * [2] Fennell & Gezelter, J. Chem. Phys. 124, 234104 (2006).
 * [3] Sun, Lin & Gezelter, J. Chem. Phys. 128, 234107 (2008).
 * [4] Vardeman, Stocker & Gezelter, J. Chem. Theory Comput. 7, 834 (2011).
 * [5] Kuang & Gezelter, Mol. Phys., 110, 691-701 (2012).
 * [6] Lamichhane, Gezelter & Newman, J. Chem. Phys. 141, 134109 (2014).
 * [7] Lamichhane, Newman & Gezelter, J. Chem. Phys. 141, 134110 (2014).
 * [8] Bhattarai, Newman & Gezelter, Phys. Rev. B 99, 094106 (2019).
 */

#include "applications/dynamicProps/HBondPersistence.hpp"
#include "utils/Constants.hpp"
#include <algorithm>

namespace OpenMD {
  HBondPersistence::HBondPersistence(SimInfo* info, const std::string& filename,
                                     const std::string& sele1,
                                     const std::string& sele2,
                                     double OOcut, double thetaCut,
                                     double OHcut)
    : TimeCorrFunc<RealType>(info, filename, sele1, sele2,
                             DataStorage::dslPosition |
                             DataStorage::dslAmat ),
    OOCut_(OOcut), thetaCut_(thetaCut), OHCut_(OHcut) {
    
    setCorrFuncType("HBondPersistence");
    setOutputName(getPrefix(dumpFilename_) + ".HBpersistence");
    
    std::stringstream params;
    params << " OOcut = " << OOCut_
           << ", thetacut = " << thetaCut_
           << ", OHcut = " << OHCut_;    
    const std::string paramString = params.str();
    setParameterString( paramString );
    
    //nFrames_ is initialized in MultipassCorrFunc:
    GIDtoDonor_.resize(nFrames_);
    DonorToGID_.resize(nFrames_);
    acceptor_.resize(nFrames_);
  }
  
  void HBondPersistence::computeFrame(int istep) {
    Molecule* mol1;
    Molecule* mol2;
    std::vector<Molecule::HBondDonor*>::iterator hbdi;
    Molecule::HBondDonor* hbd;
    std::vector<Atom*>::iterator hbai;
    Atom* hba;
    Vector3d dPos;
    Vector3d aPos;
    Vector3d hPos;
    Vector3d DH;
    Vector3d DA;
    Vector3d HA;
    Vector3d uDA;
    RealType DAdist, DHdist, HAdist, theta, ctheta;
    int ii, jj;
    int hInd, aInd, index;

    // Map of atomic global IDs to donor atoms:
    GIDtoDonor_[istep].resize( info_->getNGlobalAtoms(), -1);

    if (!uniqueSelections_) {
      seleMan2_ = seleMan1_;
    }
    
    if (evaluator1_.isDynamic()) {
      seleMan1_.setSelectionSet(evaluator1_.evaluate());
    }
    
    if (uniqueSelections_ && evaluator2_.isDynamic()) {
      seleMan2_.setSelectionSet(evaluator2_.evaluate());
    }      

    for (mol1 = seleMan1_.beginSelectedMolecule(ii);
         mol1 != NULL; mol1 = seleMan1_.nextSelectedMolecule(ii)) {

      for (mol2 = seleMan2_.beginSelectedMolecule(jj);
           mol2 != NULL; mol2 = seleMan2_.nextSelectedMolecule(jj)) {

        // loop over the possible donors in molecule 1:
        for (hbd = mol1->beginHBondDonor(hbdi); hbd != NULL;
             hbd = mol1->nextHBondDonor(hbdi)) {
          dPos = hbd->donorAtom->getPos(); 
          hPos = hbd->donatedHydrogen->getPos();
          DH = hPos - dPos; 
          currentSnapshot_->wrapVector(DH);
          DHdist = DH.length();
          
          hInd = hbd->donatedHydrogen->getGlobalIndex();
          aInd = -1;
          
          // loop over the possible acceptors in molecule 2:
          for (hba = mol2->beginHBondAcceptor(hbai); hba != NULL;
               hba = mol2->nextHBondAcceptor(hbai)) {
            aPos = hba->getPos();
            DA = aPos - dPos;
            currentSnapshot_->wrapVector(DA);
            DAdist = DA.length();
            
            // Distance criteria: are the donor and acceptor atoms
            // close enough?
            if (DAdist < OOCut_) {
              HA = aPos - hPos;
              currentSnapshot_->wrapVector(HA);
              HAdist = HA.length();
              
              ctheta = dot(DH, DA) / (DHdist * DAdist);
              theta = acos(ctheta) * 180.0 / Constants::PI;
              
              // Angle criteria: are the D-H and D-A and vectors close?
              if (theta < thetaCut_ && HAdist < OHCut_) {
                // molecule 2 is a Hbond acceptor:
                aInd = hba->getGlobalIndex();
                
                index = acceptor_[istep].size();
                GIDtoDonor_[istep][hInd] = index;
                
                acceptor_[istep].push_back(aInd);
                DonorToGID_[istep].push_back(hInd);
              }
            }            
          }                    
        }

        // loop over the possible donors in molecule 2:
        for (hbd = mol2->beginHBondDonor(hbdi); hbd != NULL;
             hbd = mol2->nextHBondDonor(hbdi)) {
          dPos = hbd->donorAtom->getPos(); 
          hPos = hbd->donatedHydrogen->getPos();
          DH = hPos - dPos; 
          currentSnapshot_->wrapVector(DH);
          DHdist = DH.length();
          
          hInd = hbd->donatedHydrogen->getGlobalIndex();
          aInd = -1;
          
          // loop over the possible acceptors in molecule 1:
          for (hba = mol1->beginHBondAcceptor(hbai); hba != NULL;
               hba = mol1->nextHBondAcceptor(hbai)) {
            aPos = hba->getPos();
            DA = aPos - dPos;
            currentSnapshot_->wrapVector(DA);
            DAdist = DA.length();
            
            // Distance criteria: are the donor and acceptor atoms
            // close enough?
            if (DAdist < OOCut_) {
              HA = aPos - hPos;
              currentSnapshot_->wrapVector(HA);
              HAdist = HA.length();

              ctheta = dot(DH, DA) / (DHdist * DAdist);
              theta = acos(ctheta) * 180.0 / Constants::PI;
              
              // Angle criteria: are the D-H and D-A and vectors close?
              if (theta < thetaCut_ && HAdist < OHCut_) {
                // molecule 1 is a Hbond acceptor:
                aInd = hba->getGlobalIndex();
                index = acceptor_[istep].size();
                GIDtoDonor_[istep][hInd] = index;
                
                acceptor_[istep].push_back(aInd);
                DonorToGID_[istep].push_back(hInd);                
              }
            }            
          }          
        }        
      }
    }
  }
  
  
  void HBondPersistence::correlation() {
    std::vector<int> s1;
    std::vector<int>::iterator i1;

    RealType corrVal;
    int index1, index2, count, gid;

    for (int i = 0; i < nFrames_; ++i) {

      RealType time1 = times_[i];           
      s1 = DonorToGID_[i];
            
      for(int j  = i; j < nFrames_; ++j) {

        // Perform a sanity check on the actual configuration times to
        // make sure the configurations are spaced the same amount the
        // sample time said they were spaced:
               
        RealType time2 = times_[j];

        if ( fabs( (time2 - time1) - (j-i)*deltaTime_ ) > 1.0e-4 ) {
          sprintf(painCave.errMsg,
                  "HBondPersistence::correlation Error: sampleTime (%f)\n"
                  "\tin %s does not match actual time-spacing between\n"
                  "\tconfigurations %d (t = %f) and %d (t = %f).\n",
                  deltaTime_, dumpFilename_.c_str(), i, time1, j, time2);
          painCave.isFatal = 1;
          simError();  
        }
        
        int timeBin = int ((time2 - time1) / deltaTime_ + 0.5);        

        corrVal = 0.0;
        count = s1.size();

        // loop over the H-bond donors found in frame i:
        
        for (i1 = s1.begin(); i1 != s1.end(); ++i1) {

          // gid is the global ID of H-bond donor index1 in frame i:
          gid = *i1;
          
          index1 = GIDtoDonor_[i][gid];
          
          // find matching donor in frame j:
          index2 = GIDtoDonor_[j][gid];

          // sometimes the donor doesn't have a Hydrogen bond in a
          // given frame, so the index will default to -1:
          
          if ( index2 < 0 ) {
            corrVal += 0;
          } else {
            
            if (acceptor_[i][index1] == acceptor_[j][index2] )
              corrVal += 1;
            else
              corrVal += 0;
          }
        }          
        histogram_[timeBin] += corrVal;
        count_[timeBin] += count;
      }
    }
  }

  void HBondPersistence::postCorrelate() {
    for (unsigned int i =0 ; i < nTimeBins_; ++i) {
      if (count_[i] > 0) {
	histogram_[i] /= count_[i];
      } else {
        histogram_[i] = 0;
      }
    }
  }
  
}
