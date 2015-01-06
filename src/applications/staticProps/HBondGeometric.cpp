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
 * [4] Kuang & Gezelter,  J. Chem. Phys. 133, 164101 (2010).
 * [4] , Stocker & Gezelter, J. Chem. Theory Comput. 7, 834 (2011). 
 */
 
#include "applications/staticProps/HBondGeometric.hpp"
#include "utils/simError.h"
#include "io/DumpReader.hpp"
#include "primitives/Molecule.hpp"
#include "utils/NumericConstant.hpp"

#include <vector>

namespace OpenMD {

  HBondGeometric::HBondGeometric(SimInfo* info, 
                                 const std::string& filename, 
                                 const std::string& sele1,
                                 const std::string& sele2,
                                 double rCut, double thetaCut, int nbins) :
    StaticAnalyser(info, filename),
    selectionScript1_(sele1), evaluator1_(info), seleMan1_(info),
    selectionScript2_(sele2), evaluator2_(info), seleMan2_(info){
    
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

    rCut_ = rCut;
    thetaCut_ = thetaCut;
    nBins_ = nbins;

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
    Molecule* mol;
    StuntDouble* sd1;
    StuntDouble* sd2;
    RigidBody* rb1;
    RigidBody* rb2;
    SimInfo::MoleculeIterator mi;
    Molecule::RigidBodyIterator rbIter;
    Molecule::IntegrableObjectIterator ioi;
    int ii, jj;
    std::string rbName;
    std::vector<Atom *> atoms1;
    std::vector<Atom *> atoms2;
    std::vector<Atom *>::iterator ai1;
    std::vector<Atom *>::iterator ai2;
    Vector3d O1pos, O2pos;
    Vector3d H1apos, H1bpos, H2apos, H2bpos;
    int nHB, nA, nD;

    DumpReader reader(info_, dumpFilename_);    
    int nFrames = reader.getNFrames();
    frameCounter_ = 0;

    for (int istep = 0; istep < nFrames; istep += step_) {
      reader.readFrame(istep);
      frameCounter_++;
      currentSnapshot_ = info_->getSnapshotManager()->getCurrentSnapshot();
     
      // update the positions of atoms which belong to the rigidbodies
      
      for (mol = info_->beginMolecule(mi); mol != NULL; 
           mol = info_->nextMolecule(mi)) {
        for (rb1 = mol->beginRigidBody(rbIter); rb1 != NULL; 
             rb1 = mol->nextRigidBody(rbIter)) {
          rb1->updateAtoms();
        }        
      }           
      
      if  (evaluator1_.isDynamic()) {
        seleMan1_.setSelectionSet(evaluator1_.evaluate());
      }
      if  (evaluator2_.isDynamic()) {
        seleMan2_.setSelectionSet(evaluator2_.evaluate());
      }
      
      for (sd1 = seleMan1_.beginSelected(ii); sd1 != NULL; sd1 = seleMan1_.nextSelected(ii)) {
        if (sd1->isRigidBody()) {
          rb1 = dynamic_cast<RigidBody*>(sd1);
          atoms1 = rb1->getAtoms();
          
          int nH = 0;
          int nO = 0;
          
          for (ai1 = atoms1.begin(); ai1 != atoms1.end(); ++ai1) {
            std::string atName =  (*ai1)->getType();
            // query the force field for the AtomType associated with this
            // atomTypeName:
            AtomType* at = ff_->getAtomType(atName);
            // get the chain of base types for this atom type:
            std::vector<AtomType*> ayb = at->allYourBase();
            // use the last type in the chain of base types for the name:
            std::string bn = ayb[ayb.size()-1]->getName();
            
            bool isH = bn.compare("H") == 0 ? true : false;
            bool isO = bn.compare("O") == 0 ? true : false;
            
            if (isO && nO == 0) {
              O1pos = (*ai1)->getPos();
              nO++;
            }
            if (isH) {
              if (nH == 0) {
                H1apos =  (*ai1)->getPos();
              }
              if (nH == 1) {
                H1bpos =  (*ai1)->getPos();
              }
              nH++;
            }
          }
        }


        nHB = 0;
        nA = 0;
        nD = 0;
        
        for (sd2 = seleMan2_.beginSelected(jj); sd2 != NULL; sd2 = seleMan2_.nextSelected(jj)) {

          if (sd1 == sd2) continue;
           
          if (sd2->isRigidBody()) {
            rb2 = dynamic_cast<RigidBody*>(sd2);
            atoms2 = rb2->getAtoms();
            
            int nH = 0;
            int nO = 0;
            
            for (ai2 = atoms2.begin(); ai2 != atoms2.end(); ++ai2) {
              std::string atName =  (*ai2)->getType();
              // query the force field for the AtomType associated with this
              // atomTypeName:
              AtomType* at = ff_->getAtomType(atName);
              // get the chain of base types for this atom type:
              std::vector<AtomType*> ayb = at->allYourBase();
              // use the last type in the chain of base types for the name:
              std::string bn = ayb[ayb.size()-1]->getName();
              
              bool isH = bn.compare("H") == 0 ? true : false;
              bool isO = bn.compare("O") == 0 ? true : false;

              if (isO && nO == 0) {
                O2pos = (*ai2)->getPos();
                  nO++;
              }
              if (isH) {
                if (nH == 0) {
                  H2apos =  (*ai2)->getPos();
                }
                if (nH == 1) {
                  H2bpos =  (*ai2)->getPos();
                }
                nH++;
              }
            }
            
            // Do our testing:
            Vector3d Odiff = O2pos - O1pos;
            currentSnapshot_->wrapVector(Odiff);
            RealType Odist = Odiff.length();
            if (Odist < rCut_) {
              // OH vectors:
              Vector3d HO1a = H1apos - O1pos;
              Vector3d HO1b = H1bpos - O1pos;
              Vector3d HO2a = H2apos - O2pos;
              Vector3d HO2b = H2bpos - O2pos;
              // wrapped in case a molecule is split across boundaries:
              currentSnapshot_->wrapVector(HO1a);
              currentSnapshot_->wrapVector(HO1b);
              currentSnapshot_->wrapVector(HO2a);
              currentSnapshot_->wrapVector(HO2a);
              // cos thetas:
              RealType ctheta1a = dot(HO1a, Odiff) / (Odist * HO1a.length());
              RealType ctheta1b = dot(HO1b, Odiff) / (Odist * HO1b.length());
              RealType ctheta2a = dot(HO2a, -Odiff) / (Odist * HO2a.length());
              RealType ctheta2b = dot(HO2b, -Odiff) / (Odist * HO2b.length());

              RealType theta1a = acos(ctheta1a) * 180.0 / M_PI;
              RealType theta1b = acos(ctheta1b) * 180.0 / M_PI;
              RealType theta2a = acos(ctheta2a) * 180.0 / M_PI;
              RealType theta2b = acos(ctheta2b) * 180.0 / M_PI;

              if (theta1a < thetaCut_) {
                // molecule 1 is a Hbond donor:
                nHB++;
                nD++;
              }
              if (theta1b < thetaCut_) {
                // molecule 1 is a Hbond donor:
                nHB++;
                nD++;
              }
              if (theta2a < thetaCut_) {
                // molecule 1 is a Hbond acceptor:
                nHB++;
                nA++;
              }
              if (theta2b < thetaCut_) {
                // molecule 1 is a Hbond acceptor:
                nHB++;
                nA++;
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
    cerr << "nSelected = " << nSelected_ << "\n";

    if (osq.is_open()) {
      
      osq << "# HydrogenBonding Statistics\n";
      osq << "# selection1: (" << selectionScript1_ << ")"
          << "\tselection2: (" << selectionScript2_ <<  ")\n";
      osq << "# p(nHBonds)\tp(nAcceptor)\tp(nDonor)\n";
      // Normalize by number of frames and write it out:
      for (int i = 0; i < nBins_; ++i) {
        osq << i;
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

      


