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
 * [4] Kuang & Gezelter,  J. Chem. Phys. 133, 164101 (2010).
 * [4] , Stocker & Gezelter, J. Chem. Theory Comput. 7, 834 (2011). *
 *  Created by J. Daniel Gezelter on 09/26/06.
 *  @author  J. Daniel Gezelter
 *  @version $Id$
 *
 */
 
#include "applications/staticProps/BondOrderParameter.hpp"
#include "utils/simError.h"
#include "io/DumpReader.hpp"
#include "primitives/Molecule.hpp"
#include "utils/NumericConstant.hpp"
#include "math/Wigner3jm.hpp"

using namespace MATPACK;
namespace OpenMD {

  BondOrderParameter::BondOrderParameter(SimInfo* info, 
                                         const std::string& filename, 
                                         const std::string& sele,
                                         double rCut, int nbins) : StaticAnalyser(info, filename), selectionScript_(sele), evaluator_(info), seleMan_(info){
    
    setOutputName(getPrefix(filename) + ".bo");

    evaluator_.loadScriptString(sele);
    if (!evaluator_.isDynamic()) {
      seleMan_.setSelectionSet(evaluator_.evaluate());
    }

    // Set up cutoff radius and order of the Legendre Polynomial:

    rCut_ = rCut;
    nBins_ = nbins;
    Qcount_.resize(lMax_+1);
    Wcount_.resize(lMax_+1);

    // Q can take values from 0 to 1

    MinQ_ = 0.0;
    MaxQ_ = 1.1;
    deltaQ_ = (MaxQ_ - MinQ_) / nbins;

    // W_6 for icosahedral clusters is 11 / sqrt(4199) = 0.169754, so we'll
    // use values for MinW_ and MaxW_ that are slightly larger than this:

    MinW_ = -1.1;
    MaxW_ = 1.1;
    deltaW_ = (MaxW_ - MinW_) / nbins;

    // Make arrays for Wigner3jm
    RealType* THRCOF = new RealType[2*lMax_+1];
    // Variables for Wigner routine
    RealType lPass, m1Pass, m2m, m2M;
    int error, mSize;
    mSize = 2*lMax_+1;

    for (int l = 0; l <= lMax_; l++) {
      lPass = (RealType)l;
      for (int m1 = -l; m1 <= l; m1++) {
        m1Pass = (RealType)m1;

        std::pair<int,int> lm = std::make_pair(l, m1);
        
        // Zero work array
        for (int ii = 0; ii < 2*l + 1; ii++){
          THRCOF[ii] = 0.0;
        }

        // Get Wigner coefficients
        Wigner3jm(lPass, lPass, lPass, 
                  m1Pass, m2m, m2M, 
                  THRCOF, mSize, error);
        
        m2Min[lm] = (int)floor(m2m);
        m2Max[lm] = (int)floor(m2M);
        
        for (int mmm = 0; mmm <= (int)(m2M - m2m); mmm++) {
          w3j[lm].push_back(THRCOF[mmm]);
        }
      }
    }
    delete [] THRCOF;
    THRCOF = NULL;
  }
  
  BondOrderParameter::~BondOrderParameter() {
    Q_histogram_.clear();
    W_histogram_.clear();
    for (int l = 0; l <= lMax_; l++) {
      for (int m = -l; m <= l; m++) {
        w3j[std::make_pair(l,m)].clear();
      }
    }
    w3j.clear();
    m2Min.clear();
    m2Max.clear();
  }
  
  void BondOrderParameter::initializeHistogram() {
    for (int bin = 0; bin < nBins_; bin++) {
      for (int l = 0; l <= lMax_; l++) {
        Q_histogram_[std::make_pair(bin,l)] = 0;
        W_histogram_[std::make_pair(bin,l)] = 0;
      }
    }
  }

  void BondOrderParameter::process() {
    Molecule* mol;
    Atom* atom;
    RigidBody* rb;
    int myIndex;
    SimInfo::MoleculeIterator mi;
    Molecule::RigidBodyIterator rbIter;
    Molecule::AtomIterator ai;
    StuntDouble* sd;
    Vector3d vec;
    RealType costheta;
    RealType phi;
    RealType r;
    std::map<std::pair<int,int>,ComplexType> q;
    std::vector<RealType> q_l;
    std::vector<RealType> q2;
    std::vector<ComplexType> w;
    std::vector<ComplexType> w_hat;
    std::map<std::pair<int,int>,ComplexType> QBar;
    std::vector<RealType> Q2;
    std::vector<RealType> Q;
    std::vector<ComplexType> W;
    std::vector<ComplexType> W_hat;
    int nBonds, Nbonds;
    SphericalHarmonic sphericalHarmonic;
    int i;

    DumpReader reader(info_, dumpFilename_);    
    int nFrames = reader.getNFrames();
    frameCounter_ = 0;

    q_l.resize(lMax_+1);
    q2.resize(lMax_+1);
    w.resize(lMax_+1);
    w_hat.resize(lMax_+1);

    Q2.resize(lMax_+1);
    Q.resize(lMax_+1);
    W.resize(lMax_+1);
    W_hat.resize(lMax_+1);
    Nbonds = 0;

    for (int istep = 0; istep < nFrames; istep += step_) {
      reader.readFrame(istep);
      frameCounter_++;
      currentSnapshot_ = info_->getSnapshotManager()->getCurrentSnapshot();
      
      if (evaluator_.isDynamic()) {
        seleMan_.setSelectionSet(evaluator_.evaluate());
      }

      // update the positions of atoms which belong to the rigidbodies

      for (mol = info_->beginMolecule(mi); mol != NULL; 
           mol = info_->nextMolecule(mi)) {
        for (rb = mol->beginRigidBody(rbIter); rb != NULL; 
             rb = mol->nextRigidBody(rbIter)) {
          rb->updateAtoms();
        }        
      }           
            
      // outer loop is over the selected StuntDoubles:

      for (sd = seleMan_.beginSelected(i); sd != NULL; 
           sd = seleMan_.nextSelected(i)) {

        myIndex = sd->getGlobalIndex();
        nBonds = 0;
        
        for (int l = 0; l <= lMax_; l++) {
          for (int m = -l; m <= l; m++) {
            q[std::make_pair(l,m)] = 0.0;
          }
        }
        
        // inner loop is over all other atoms in the system:
        
        for (mol = info_->beginMolecule(mi); mol != NULL; 
             mol = info_->nextMolecule(mi)) {
          for (atom = mol->beginAtom(ai); atom != NULL; 
               atom = mol->nextAtom(ai)) {

            if (atom->getGlobalIndex() != myIndex) {

              vec = sd->getPos() - atom->getPos();       

              if (usePeriodicBoundaryConditions_) 
                currentSnapshot_->wrapVector(vec);
              
              // Calculate "bonds" and build Q_lm(r) where 
              //      Q_lm = Y_lm(theta(r),phi(r))                
              // The spherical harmonics are wrt any arbitrary coordinate
              // system, we choose standard spherical coordinates 
              
              r = vec.length();
              
              // Check to see if neighbor is in bond cutoff 
              
              if (r < rCut_) { 
                costheta = vec.z() / r; 
                phi = atan2(vec.y(), vec.x());

                for (int l = 0; l <= lMax_; l++) {
                  sphericalHarmonic.setL(l);
                  for(int m = -l; m <= l; m++){
                    sphericalHarmonic.setM(m);
                    q[std::make_pair(l,m)] += sphericalHarmonic.getValueAt(costheta, phi);
                  }
                }
                nBonds++;
              }  
            }
          }
        }
        
        
        for (int l = 0; l <= lMax_; l++) {
          q2[l] = 0.0;
          for (int m = -l; m <= l; m++){
            q[std::make_pair(l,m)] /= (RealType)nBonds; 
            q2[l] += norm(q[std::make_pair(l,m)]);
          }
          q_l[l] = sqrt(q2[l] * 4.0 * NumericConstant::PI / (RealType)(2*l + 1));
        }
        
        // Find Third Order Invariant W_l
    
        for (int l = 0; l <= lMax_; l++) {
          w[l] = 0.0;
          for (int m1 = -l; m1 <= l; m1++) {
            std::pair<int,int> lm = std::make_pair(l, m1);
            for (int mmm = 0; mmm <= (m2Max[lm] - m2Min[lm]); mmm++) {
              int m2 = m2Min[lm] + mmm;
              int m3 = -m1-m2;
              w[l] += w3j[lm][mmm] * q[lm] * 
                q[std::make_pair(l,m2)] *  q[std::make_pair(l,m3)];
            }
          }
          
          w_hat[l] = w[l] / pow(q2[l], RealType(1.5));
        }

        collectHistogram(q_l, w_hat);
        
        Nbonds += nBonds;
        for (int l = 0; l <= lMax_;  l++) {
          for (int m = -l; m <= l; m++) {
            QBar[std::make_pair(l,m)] += (RealType)nBonds*q[std::make_pair(l,m)];
          }
        }
      }
    }
      
    // Normalize Qbar2
    for (int l = 0; l <= lMax_; l++) {
      for (int m = -l; m <= l; m++){
        QBar[std::make_pair(l,m)] /= Nbonds;
      }
    }
    
    // Find second order invariant Q_l
    
    for (int l = 0; l <= lMax_; l++) {
      Q2[l] = 0.0;
      for (int m = -l; m <= l; m++){
        Q2[l] += norm(QBar[std::make_pair(l,m)]);
      }
      Q[l] = sqrt(Q2[l] * 4.0 * NumericConstant::PI / (RealType)(2*l + 1));
    }
    
    // Find Third Order Invariant W_l
    
    for (int l = 0; l <= lMax_; l++) {
      W[l] = 0.0;
      for (int m1 = -l; m1 <= l; m1++) {
        std::pair<int,int> lm = std::make_pair(l, m1);
        for (int mmm = 0; mmm <= (m2Max[lm] - m2Min[lm]); mmm++) {
          int m2 = m2Min[lm] + mmm;
          int m3 = -m1-m2;
          W[l] += w3j[lm][mmm] * QBar[lm] * 
            QBar[std::make_pair(l,m2)] * QBar[std::make_pair(l,m3)];
        }
      }
      
      W_hat[l] = W[l] / pow(Q2[l], RealType(1.5));
    }
    
    writeOrderParameter(Q, W_hat);    
  }

  void BondOrderParameter::collectHistogram(std::vector<RealType> q, 
                                            std::vector<ComplexType> what) {

    for (int l = 0; l <= lMax_; l++) {
      if (q[l] >= MinQ_ && q[l] < MaxQ_) {
        int qbin = int((q[l] - MinQ_) / deltaQ_);
        Q_histogram_[std::make_pair(qbin,l)] += 1;
        Qcount_[l]++;      
      } else {
        sprintf( painCave.errMsg,
                 "q_l value outside reasonable range\n");
        painCave.severity = OPENMD_ERROR;
        painCave.isFatal = 1;
        simError();  
      }
    }

    for (int l = 0; l <= lMax_; l++) {
      if (real(what[l]) >= MinW_ && real(what[l]) < MaxW_) {
        int wbin = int((real(what[l]) - MinW_) / deltaW_);
        W_histogram_[std::make_pair(wbin,l)] += 1;
        Wcount_[l]++;      
      } else {
        sprintf( painCave.errMsg,
                 "Re[w_hat] value (%lf) outside reasonable range\n", real(what[l]));
        painCave.severity = OPENMD_ERROR;
        painCave.isFatal = 1;
        simError();  
      }
    }

  }  


  void BondOrderParameter::writeOrderParameter(std::vector<RealType> Q, 
                                               std::vector<ComplexType> What) {
    
    std::ofstream osq((getOutputFileName() + "q").c_str());

    if (osq.is_open()) {
      
      osq << "# Bond Order Parameters\n";
      osq << "# selection: (" << selectionScript_ << ")\n";
      osq << "# \n";
      for (int l = 0; l <= lMax_; l++) {
        osq << "# <Q_" << l << ">: " << Q[l] << "\n";
      }
      // Normalize by number of frames and write it out:
      for (int i = 0; i < nBins_; ++i) {
        RealType Qval = MinQ_ + (i + 0.5) * deltaQ_;               
        osq << Qval;
        for (int l = 0; l <= lMax_; l++) {

          osq << "\t" << (RealType)Q_histogram_[std::make_pair(i,l)]/(RealType)Qcount_[l]/deltaQ_;
        }
        osq << "\n";
      }

      osq.close();

    } else {
      sprintf(painCave.errMsg, "BondOrderParameter: unable to open %s\n", 
              (getOutputFileName() + "q").c_str());
      painCave.isFatal = 1;
      simError();  
    }

    std::ofstream osw((getOutputFileName() + "w").c_str());

    if (osw.is_open()) {
      osw << "# Bond Order Parameters\n";
      osw << "# selection: (" << selectionScript_ << ")\n";
      osw << "# \n";
      for (int l = 0; l <= lMax_; l++) {
        osw << "# <W_" << l << ">: " << real(What[l]) << "\t" << imag(What[l]) << "\n";
      }
      // Normalize by number of frames and write it out:
      for (int i = 0; i < nBins_; ++i) {
        RealType Wval = MinW_ + (i + 0.5) * deltaW_;               
        osw << Wval;
        for (int l = 0; l <= lMax_; l++) {

          osw << "\t" << (RealType)W_histogram_[std::make_pair(i,l)]/(RealType)Wcount_[l]/deltaW_;
        }
        osw << "\n";
      }

      osw.close();
    } else {
      sprintf(painCave.errMsg, "BondOrderParameter: unable to open %s\n", 
              (getOutputFileName() + "w").c_str());
      painCave.isFatal = 1;
      simError();  
    }
       
  }
}
