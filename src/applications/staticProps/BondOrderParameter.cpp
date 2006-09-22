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
#include "math/SphericalHarmonic.hpp"

namespace oopse {

  BondOrderParameter::BondOrderParameter(SimInfo* info, 
                                         const std::string& filename, 
                                         const std::string& sele,
                                         double rCut, int lNumber, int nbins) : StaticAnalyser(info, filename), selectionScript_(sele), evaluator_(info), seleMan_(info){
    
    setOutputName(getPrefix(filename) + ".bo");

    evaluator_.loadScriptString(sele);
    if (!evaluator_.isDynamic()) {
      seleMan_.setSelectionSet(evaluator_.evaluate());
    }

    // Set up cutoff radius and order of the Legendre Polynomial:

    lNumber_ = lNumber;
    rCut_ = rCut;
    mSize_ = 2*lNumber_+1;   

    // Q can take values from 0 to 1

    MinQ_ = 0.0;
    MaxQ_ = 3.0;
    deltaQ_ = (MaxQ_ - MinQ_) / nbins;
    Q_histogram_.resize(nbins);

    // W_6 for icosahedral clusters is 11 / sqrt(4199) = 0.169754, so we'll
    // use values for MinW_ and MaxW_ that are slightly larger than this:

    MinW_ = -0.18;
    MaxW_ = 0.18;
    deltaW_ = (MaxW_ - MinW_) / nbins;
    W_histogram_.resize(nbins);

  }

  BondOrderParameter::~BondOrderParameter() {
    Q_histogram_.clear();
    W_histogram_.clear();
  }

  void BondOrderParameter::initalizeHistogram() {
    std::fill(Q_histogram_.begin(), Q_histogram_.end(), 0);
    std::fill(W_histogram_.begin(), W_histogram_.end(), 0);
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
    std::map<int,ComplexType> q_lm;
    std::map<int,ComplexType> QBar_lm;
    RealType QSq_l;
    RealType Q_l;
    ComplexType W_l;
    ComplexType W_l_hat;
    int nBonds, Nbonds;
    SphericalHarmonic sphericalHarmonic;
    int i, j;
    // Make arrays for Wigner3jm
    double* THRCOF = new double[mSize_];
    // Variables for Wigner routine
    double l_ = (double)lNumber_;
    double m1Pass, m2Min, m2Max;
    int error, m1, m2, m3;

    // Set the l for the spherical harmonic, it doesn't change
    sphericalHarmonic.setL(lNumber_);

    DumpReader reader(info_, dumpFilename_);    
    int nFrames = reader.getNFrames();
    frameCounter_ = 0;

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
        for (int m = -lNumber_; m <= lNumber_; m++) {
          q_lm[m] = 0.0;
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
                
                for(int m = -lNumber_; m <= lNumber_; m++){
                  sphericalHarmonic.setM(m);
                  q_lm[m] += sphericalHarmonic.getValueAt(costheta, phi);
                }
                nBonds++;
              }  
            }
          }
        }
        RealType ql = 0.0;
        for(int m=-lNumber_; m<=lNumber_; m++) {          
          ql += norm(QBar_lm[m]);
        }        
        ql *= 4.0*NumericConstant::PI/(RealType)(2*lNumber_+1);
        collectHistogram(sqrt(ql)/(RealType)nBonds);

        Nbonds += nBonds;
        for (int m=-lNumber_; m<=lNumber_; m++) {
          QBar_lm[m] += q_lm[m];
        }
      }
    }

    // Normalize Qbar2
    for (int m = -lNumber_;m <= lNumber_; m++){
      QBar_lm[m] /= Nbonds;
    }
    
    // Find second order invariant Q_l
    
    QSq_l = 0.0;
    for (int m = -lNumber_; m <= lNumber_; m++){
      QSq_l += norm(QBar_lm[m]);
    }
    
    std::cout << "qsl = " << QSq_l << "\n";
    Q_l = sqrt(QSq_l * 4.0 * NumericConstant::PI / (RealType)(2*lNumber_ + 1));
    
    // Find Third Order Invariant W_l
    
    W_l = 0.0;
    for (int m1 = -lNumber_; m1 <= lNumber_; m1++) {
      // Zero work array
      for (int ii = 0; ii < mSize_; ii++){
        THRCOF[ii] = 0.0;
      }
      // Get Wigner coefficients
      m1Pass = (double)m1;
      
      Wigner3jm(&l_, &l_, &l_, 
                &m1Pass, &m2Min, &m2Max, 
                THRCOF, &mSize_, &error);
      
      for (int mmm = 0; mmm < (int)(m2Max - m2Min); mmm++) {
        m2 = (int)floor(m2Min) + mmm;
        m3 = -m1-m2;
        W_l += THRCOF[mmm] * QBar_lm[m1] * QBar_lm[m2] * QBar_lm[m3];
      }
    }
    
    W_l_hat = W_l / pow(QSq_l, 1.5);               
    
    writeOrderParameter(Q_l, real(W_l_hat));    
  }

  void BondOrderParameter::collectHistogram(RealType q_l) {

    if (q_l >= MinQ_ && q_l < MaxQ_) {
      int qbin = (q_l - MinQ_) / deltaQ_;
      Q_histogram_[qbin] += 1;
      Qcount_++;
      sumQ_ += q_l;
      sumQ2_ += q_l * q_l;
    } else {
      sprintf( painCave.errMsg,
               "q_l value outside reasonable range\n");
      painCave.severity = OOPSE_ERROR;
      painCave.isFatal = 1;
      simError();  
    }

  }  


  void BondOrderParameter::writeOrderParameter(RealType ql, RealType Wlhat) {

    std::ofstream os(getOutputFileName().c_str());

    if (os.is_open()) {
      
      os << "# Bond Order Parameters\n";
      os << "# selection: (" << selectionScript_ << ")\n";
      os << "# \n";
      os << "# <Q_" << lNumber_ << ">: " << ql << "\n";
      os << "# <W_" << lNumber_ << ">: " << Wlhat << "\n";
      // Normalize by number of frames and write it out:
      for (int i = 0; i < Q_histogram_.size(); ++i) {
        RealType Qval = MinQ_ + (i + 0.5) * deltaQ_;
        osq << Qval << "\t" << (RealType)Q_histogram_[i] / (RealType)Qcount_ << "\n";
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
