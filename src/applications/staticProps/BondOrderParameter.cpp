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


/* Creates orientational bond order parameters as outlined by
 *     Bond-orientaional order in liquids and glasses, Steinhart,Nelson,Ronchetti
 *     Phys Rev B, 28,784,1983
 * 
 */
 
#include "applications/staticProps/BondOrderParameter.hpp"
#include "utils/simError.h"
#include "io/DumpReader.hpp"
#include "primitives/Molecule.hpp"
#include "utils/NumericConstant.hpp"

namespace oopse {

  BondOrderParameter::BondOrderParameter(SimInfo* info, 
                                         const std::string& filename, 
                                         const std::string& sele,
                                         double rCut, int lMax, int nbins) : StaticAnalyser(info, filename), selectionScript_(sele), evaluator_(info), seleMan_(info){
    
    setOutputName(getPrefix(filename) + ".bo");

    evaluator_.loadScriptString(sele);
    if (!evaluator_.isDynamic()) {
      seleMan_.setSelectionSet(evaluator_.evaluate());
    }

    // Set up cutoff radius and order of the Legendre Polynomial:

    lMax_ = lMax;
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

    MinW_ = -0.18;
    MaxW_ = 0.18;
    deltaW_ = (MaxW_ - MinW_) / nbins;
  }

  BondOrderParameter::~BondOrderParameter() {
    Q_histogram_.clear();
    W_histogram_.clear();
  }

  void BondOrderParameter::initalizeHistogram() {
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
    RealType dist;
    std::map<std::pair<int,int>,ComplexType> q;
    std::vector<RealType> q_l;
    std::map<std::pair<int,int>,ComplexType> QBar;
    std::vector<RealType> Q2;
    std::vector<RealType> Q;
    std::vector<ComplexType> W;
    std::vector<ComplexType> W_hat;
    int nBonds, Nbonds;
    SphericalHarmonic sphericalHarmonic;
    int i, j;
    // Make arrays for Wigner3jm
    double* THRCOF = new double[2*lMax_+1];
    // Variables for Wigner routine
    double lPass, m1Pass, m2Min, m2Max;
    int error, m1, m2, m3, mSize;
    mSize = 2*lMax_+1;

    DumpReader reader(info_, dumpFilename_);    
    int nFrames = reader.getNFrames();
    frameCounter_ = 0;

    q_l.resize(lMax_+1);
    Q2.resize(lMax_+1);
    Q.resize(lMax_+1);
    W.resize(lMax_+1);
    W_hat.resize(lMax_+1);

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
          q_l[l] = 0.0;
          for(int m = -l; m <= l; m++) { 
            q_l[l] += norm(q[std::make_pair(l,m)]);
          }     
          q_l[l] *= 4.0*NumericConstant::PI/(RealType)(2*l + 1);
          q_l[l] = sqrt(q_l[l])/(RealType)nBonds;
        }
        collectHistogram(q_l);
        
        Nbonds += nBonds;
        for (int l = 0; l <= lMax_;  l++) {
          for (int m = -l; m <= l; m++) {
            QBar[std::make_pair(l,m)] += q[std::make_pair(l,m)];
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
      lPass = (double)l;
      for (int m1 = -l; m1 <= l; m1++) {
        // Zero work array
        for (int ii = 0; ii < 2*l + 1; ii++){
          THRCOF[ii] = 0.0;
        }
        // Get Wigner coefficients
        m1Pass = (double)m1;
        
        Wigner3jm(&lPass, &lPass, &lPass, 
                  &m1Pass, &m2Min, &m2Max, 
                  THRCOF, &mSize, &error);
        
        for (int mmm = 0; mmm < (int)(m2Max - m2Min); mmm++) {
          m2 = (int)floor(m2Min) + mmm;
          m3 = -m1-m2;
          W[l] += THRCOF[mmm] * 
            QBar[std::make_pair(l,m1)] * 
            QBar[std::make_pair(l,m2)] * 
            QBar[std::make_pair(l,m3)];
        }
      }
      
      W_hat[l] = W[l] / pow(Q2[l], 1.5);
    }
    
    writeOrderParameter(Q, W_hat);    
  }

  void BondOrderParameter::collectHistogram(std::vector<RealType> q) {

    for (int l = 0; l <= lMax_; l++) {
      if (q[l] >= MinQ_ && q[l] < MaxQ_) {
        int qbin = (q[l] - MinQ_) / deltaQ_;
        Q_histogram_[std::make_pair(qbin,l)] += 1;
        Qcount_[l]++;      
      } else {
        sprintf( painCave.errMsg,
                 "q_l value outside reasonable range\n");
        painCave.severity = OOPSE_ERROR;
        painCave.isFatal = 1;
        simError();  
      }
    }

  }  


  void BondOrderParameter::writeOrderParameter(std::vector<RealType> Q, std::vector<ComplexType> What) {
    
    std::ofstream os(getOutputFileName().c_str());
    
    if (os.is_open()) {
      
      os << "# Bond Order Parameters\n";
      os << "# selection: (" << selectionScript_ << ")\n";
      for (int l = 0; l <= lMax_; l++) {
        os << "# \n";
        os << "# <Q_" << l << ">: " << Q[l] << "\n";
        os << "# <W_" << l << ">: " << real(What[l]) << "\n";
      }
      // Normalize by number of frames and write it out:
      for (int i = 0; i < nBins_; ++i) {
        RealType Qval = MinQ_ + (i + 0.5) * deltaQ_;               
        os << Qval;
        for (int l = 0; l <= lMax_; l++) {
          os << "\t" << (RealType)Q_histogram_[std::make_pair(i,l)] / (RealType)Qcount_[l];
        }
        os << "\n";
      }

      os.close();

    } else {
      sprintf(painCave.errMsg, "BondOrderParameter: unable to open %s\n", 
              getOutputFileName().c_str());
      painCave.isFatal = 1;
      simError();  
    }
  }
}
