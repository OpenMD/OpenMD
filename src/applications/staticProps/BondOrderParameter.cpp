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
#include "math/RealSphericalHarmonic.hpp"
namespace oopse {


  BondOrderParameter::BondOrderParameter(SimInfo* info, 
                                         const std::string& filename, 
                                         const std::string& sele,
                                         double rCut, int lNumber, int nbins)
    : StaticAnalyser(info, filename), selectionScript_(sele), 
      evaluator_(info), seleMan_(info){
    
    setOutputName(getPrefix(filename) + ".obo");

    evaluator_.loadScriptString(sele);
    if (!evaluator_.isDynamic()) {
      seleMan_.setSelectionSet(evaluator_.evaluate());
    }

    // Set up cutoff radius and order of the Legendre Polynomial:

    lNumber_ = lNumber;
    rCut_ = rCut;
    mSize_ = 2*lNumber_+1;    

    // Set the l for the spherical harmonic, it doesn't change

    sphericalHarmonic.setL(lNumber_);

    delta_Q = 1.0 / nbins;
    delta_W = 2.0 / nbins;

    Q_histogram_.resize(nbins);
    W_histogram_.resize(nbins);

  }

  void BondOrderParameter::initalizeHistogram() {
    std::fill(Q_histogram_.begin(), Q_histogram_.end(), 0);
    std::fill(W_histogram_.begin(), W_histogram_.end(), 0);
  }

  void BondOrderParameter::process() {
    Molecule* mol;
    Atom* atom;
    RigidBody* rb;
    SimInfo::MoleculeIterator mi;
    Molecule::RigidBodyIterator rbIter;
    Molecule::AtomIterator ai;
    StuntDouble* sd;
    RealType theta;
    RealType phi;
    RealType r;
    RealType dist;
    std::map<int, RealType> QBar_lm;
    RealType QSq_l;
    RealType Q_l;
    int nBonds;
    RealSphericalHarmonic sphericalHarmonic;
    int i, j;
  
  
    DumpReader reader(info_, dumpFilename_);    
    int nFrames = reader.getNFrames();


    for (int istep = 0; istep < nFrames; istep += step_) {
      reader.readFrame(istep);
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

        // For this central atom, zero out nBonds and QBar_lm

        nBonds = 0;
       
        for (int m = -lNumber_; m <= lNumber_; m++) {
          QBar_lm[m] = 0.0;
        }
        
        // inner loop is over all other atoms in the system:
        
        for (mol = info_->beginMolecule(mi); mol != NULL; 
             mol = info_->nextMolecule(mi)) {
          for (atom = mol->beginAtom(ai); atom != NULL; 
               atom = mol->nextAtom(ai)) {


            Vector3d vec = sd->getPos() - atom->getPos();       
            currentSnapshot_->wrapVector(vec);
            
            // Calculate "bonds" and build Q_lm(r) where 
            //      Q_lm = Y_lm(theta(r),phi(r))                
            // The spherical harmonics are wrt any arbitrary coordinate
            // system, we choose standard spherical coordinates 
            
            r = sqrt(pow(vec.x(),2)+pow(vec.y(),2)+pow(vec.z(),2));
            
            // Check to see if neighbor is in bond cutoff 
            
            if (r < rCut_) {  		
              theta = atan2(vec.y(), vec.x());
              phi = acos(vec.z()/r);
              for(int m = -lNumber_; m <= lNumber_; m++){
                sphericalHarmonic.setM(m);
                QBar_lm[m] += sphericalHarmonic.getValueAt(theta,phi);
              }
              nBonds++;
            }  
          }
        }
        
        // Normalize Qbar
        for (int m = -lNumber_;m <= lNumber_; m++){
          QBar_lm[m] /= nBonds;	
        }

        // Find second order invariant Q_l

        QSq_l = 0.0;
        for (int m = -lNumber_; m <= lNumber_; m++){
          QSq_l += pow(QBar_lm[m], 2);	
        }
        Q_l = sqrt(QSq_l*(4.0 * NumericConstant::PI / (2.0*(RealType)lNumber_ + 1)));
     
        // Find Third Order Invariant W_l

        // Make arrays for Wigner3jm
        double* THRCOF = new double[mSize_];
        // Variables for Wigner routine
        double l_ = (double)lNumber_;
        double m2Min, m2Max;
        int error, m1, m2, m3;
        
        W_l_ = 0.0;
        for (int m1 = -lNumber_; m1 <= lNumber_; m1++) {
          // Zero work array
          for (int ii = 0; ii < mSize_; ii+){
            THRCOF[i] = 0.0;
          }
          // Get Wigner coefficients
          Wigner3jm(&l_, &l_, &l_, &(double)m1, &m2Min, &m2Max, THRCOF, &mSize_, &error);
          for (int m_index = 1; i < (int)(m2Max - m2Min-1.0); m_index++) {
            m2 = floor(m2Min) + m_index - 1;
            m3 = -m1-m2;
            W_l_ += THRCOF[m_index]*QBar_lm[m1+lNumber_]*QBar_lm[m2+lNumber_]*QBar_lm[m3+lNumber_];
          }
        }

        W_l_hat = W_l_ / pow(QSq_l, 1.5);

        // accumulate histogram data for Q_l and W_l_hat:

        collectHistogram(Q_l, W_l_hat);
                
      }
    }
    
    // Normalize by number of frames
    for (int m = -lNumber_; m <= lNumber_; m++){
      QBar_lm[m] /=  nFrames;	
    }
    
     
     
     
     
    writeOrderParameter();
    
  }

  
  void BondOrderParameter::processHistogram() {
    
    int nPairs = getNPairs();
    RealType volume = info_->getSnapshotManager()->getCurrentSnapshot()->getVolume();
    RealType pairDensity = nPairs /volume * 2.0;
    RealType pairConstant = ( 4.0 * NumericConstant::PI * pairDensity ) / 3.0;

    for(int i = 0 ; i < histogram_.size(); ++i){

      RealType rLower = i * deltaR_;
      RealType rUpper = rLower + deltaR_;
      RealType volSlice = ( rUpper * rUpper * rUpper ) - ( rLower * rLower * rLower );
      RealType nIdeal = volSlice * pairConstant;

      avgGofr_[i] += histogram_[i] / nIdeal;    
    }

  }

  void BondOrderParameter::collectHistogram(RealType Q_l, RealType W_l_hat) {

    if (Q_l < Max_Q) {
      int whichBin = Q_l / deltaQ_;
      Q_histogram_[whichBin] += 1;
    }
    if (W_l_hat < Max_W) {
      int whichBin = W_l_hat / deltaW_;
      W_histogram_[whichBin] += 1;
    }
  }
  

  void BondOrderParameter::writeOrderParameter() {

    std::ofstream os(getOutputFileName().c_str());
    os << "#Bond Order Parameter\n";
    os << "#selection: (" << selectionScript_ << ")\n";

    for (std::size_t i = 0; i < orderParams_.size(); ++i) {
      os <<  orderParams_[i].p2 << "\t"
         <<  orderParams_[i].director[0] << "\t"
         <<  orderParams_[i].director[1] << "\t"
         <<  orderParams_[i].director[2] << "\t"
         <<  orderParams_[i].angle << "\n";

    }
  }



}
