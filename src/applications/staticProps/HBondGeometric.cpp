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
 
#include "applications/staticProps/HBondGeometric.hpp"
#include "utils/simError.h"
#include "io/DumpReader.hpp"
#include "primitives/Molecule.hpp"
#include "utils/Constants.hpp"

#include <vector>

namespace OpenMD {

  HBondGeometric::HBondGeometric(SimInfo* info, 
                                 const std::string& filename, 
                                 const std::string& sele1,
                                 const std::string& sele2,
                                 double rCut, double thetaCut, int nbins) :
    StaticAnalyser(info, filename, nbins),
    selectionScript1_(sele1), seleMan1_(info), evaluator1_(info),
    selectionScript2_(sele2), seleMan2_(info), evaluator2_(info) {
    
    setOutputName(getPrefix(filename) + ".hbg");

    ff_ = info_->getForceField();

    evaluator1_.loadScriptString(sele1);
    if (!evaluator1_.isDynamic()) {
      seleMan1_.setSelectionSet(evaluator1_.evaluate());
    }
    evaluator2_.loadScriptString(sele2);
    if (!evaluator2_.isDynamic()) {
      seleMan2_.setSelectionSet(evaluator2_.evaluate());
    }


    // Set up cutoff values:
    nBins_ = nbins;
    rCut_ = rCut;
    thetaCut_ = thetaCut;
    
    nHBonds_.resize(nBins_);
    nDonor_.resize(nBins_);
    nAcceptor_.resize(nBins_);

    initializeHistogram();
  }

  HBondGeometric::~HBondGeometric() {
    nHBonds_.clear();
    nDonor_.clear();
    nAcceptor_.clear(); 
  }
  
  void HBondGeometric::initializeHistogram() {
    std::fill(nHBonds_.begin(),   nHBonds_.end(),   0);
    std::fill(nDonor_.begin(),    nDonor_.end(),    0);
    std::fill(nAcceptor_.begin(), nAcceptor_.end(), 0);
    nSelected_ = 0;
  }
  
  void HBondGeometric::process() {
    Molecule* mol1;
    Molecule* mol2;
    Molecule::HBondDonor* hbd1;
    Molecule::HBondDonor* hbd2;
    std::vector<Molecule::HBondDonor*>::iterator hbdi;
    std::vector<Molecule::HBondDonor*>::iterator hbdj;
    std::vector<Atom*>::iterator hbai;
    std::vector<Atom*>::iterator hbaj;
    Atom* hba1;
    Atom* hba2;
    Vector3d dPos;
    Vector3d aPos;
    Vector3d hPos;
    Vector3d DH;
    Vector3d DA;
    RealType DAdist, DHdist, theta, ctheta;
    int ii, jj;
    int nHB, nA, nD;

    DumpReader reader(info_, dumpFilename_);    
    int nFrames = reader.getNFrames();
    frameCounter_ = 0;

    for (int istep = 0; istep < nFrames; istep += step_) {
      reader.readFrame(istep);
      frameCounter_++;
      currentSnapshot_ = info_->getSnapshotManager()->getCurrentSnapshot();
     
      if  (evaluator1_.isDynamic()) {
        seleMan1_.setSelectionSet(evaluator1_.evaluate());
      }
      if  (evaluator2_.isDynamic()) {
        seleMan2_.setSelectionSet(evaluator2_.evaluate());
      }
      
      for (mol1 = seleMan1_.beginSelectedMolecule(ii);
           mol1 != NULL; mol1 = seleMan1_.nextSelectedMolecule(ii)) {

        // We're collecting statistics on the molecules in selection 1:
        nHB = 0;
        nA = 0;
        nD = 0;
        
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
                theta = acos(ctheta) * 180.0 / Constants::PI;

                // Angle criteria: are the D-H and D-A and vectors close?
                if (theta < thetaCut_) {
                  // molecule 1 is a Hbond donor:
                  nHB++;
                  nD++;
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
                theta = acos(ctheta) * 180.0 / Constants::PI;
                // Angle criteria: are the D-H and D-A and vectors close?
                if (theta < thetaCut_) {
                  // molecule 1 is a Hbond acceptor:
                  nHB++;
                  nA++;
                }                
              }
            }
          }
        }                 
        collectHistogram(nHB, nA, nD);
      }
    }
    writeHistogram();
  }
 
        
  void HBondGeometric::collectHistogram(int nHB, int nA, int nD) {
    nHBonds_[nHB] += 1;
    nAcceptor_[nA] += 1;
    nDonor_[nD] += 1;
    nSelected_++;
  }


  void HBondGeometric::writeHistogram() {
        
    std::ofstream osq(getOutputFileName().c_str());

    if (osq.is_open()) {
      
      osq << "# HydrogenBonding Statistics\n";
      osq << "# selection1: (" << selectionScript1_ << ")"
          << "\tselection2: (" << selectionScript2_ <<  ")\n";
      osq << "# molecules in selection1: " << nSelected_ << "\n";
      osq << "# nHBonds\tnAcceptor\tnDonor\tp(nHBonds)\tp(nAcceptor)\tp(nDonor)\n";
      // Normalize by number of frames and write it out:
      for (int i = 0; i < nBins_; ++i) {
        osq << i;
        osq << "\t" << nHBonds_[i];
        osq << "\t" << nAcceptor_[i];
        osq << "\t" << nDonor_[i];
	osq << "\t" << (RealType) (nHBonds_[i]) / nSelected_;
        osq << "\t" << (RealType) (nAcceptor_[i]) / nSelected_;
        osq << "\t" << (RealType) (nDonor_[i]) / nSelected_;
        osq << "\n";
      }
      osq.close();
      
    } else {
      sprintf(painCave.errMsg, "HBondGeometric: unable to open %s\n", 
              (getOutputFileName() + "q").c_str());
      painCave.isFatal = 1;
      simError();  
    }
  }
}

      


