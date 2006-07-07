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


BondOrderParameter::BondOrderParameter(SimInfo* info, const std::string& filename, const std::string& sele1,
 const std::string& sele2, double rCut, int lNumber)
  : StaticAnalyser(info, filename),
    selectionScript1_(sele1), evaluator1_(info), 
    seleMan1_(info){

    setOutputName(getPrefix(filename) + ".obo");
        
    evaluator1_.loadScriptString(sele1);
    evaluator2_.loadScriptString(sele2);

    if (!evaluator1_.isDynamic()) {
      seleMan1_.setSelectionSet(evaluator1_.evaluate());
    }else {
        sprintf( painCave.errMsg,
                 "--sele1 must be static selection\n");
        painCave.severity = OOPSE_ERROR;
        painCave.isFatal = 1;
        simError();  
    }

/* Set up cutoff radius and type of order parameter we are calcuating*/
	lNumber_ = lNumber;
	rCut_ = rCut;
	mSize_ = 2*lNumber_+1;

  int i;
  int j;
  StuntDouble* sd1;
  StuntDouble* sd2;
  for (sd1 = seleMan1_.beginSelected(i), sd2 = seleMan1_.beginSelected(j);
     sd1 != NULL && sd2 != NULL;
     sd1 = seleMan1_.nextSelected(i), sd2 = seleMan2_.nextSelected(j)) {

     sdPairs_.push_back(std::make_pair(sd1, sd2));
  }

    
  }

void BondOrderParameter::process() {
  Molecule* mol;
  RigidBody* rb;
  SimInfo::MoleculeIterator mi;
  Molecule::RigidBodyIterator rbIter;
  RealType theta;
  RealType phi;
  RealType r;
  RealType dist;
  RealType* QBar_lm;
  int nBonds;
  RealSphericalHarmonic sphericalHarmonic;
  
  
  DumpReader reader(info_, dumpFilename_);    
  int nFrames = reader.getNFrames();

   /*Set the l for the spherical harmonic, it doesn't change*/
	sphericalHarmonic.setL(lNumber_);

  for (int i = 0; i < nFrames; i += step_) {
    reader.readFrame(i);
    currentSnapshot_ = info_->getSnapshotManager()->getCurrentSnapshot();
	nBonds = 0;
    
    for (mol = info_->beginMolecule(mi); mol != NULL; mol = info_->nextMolecule(mi)) {
        //change the positions of atoms which belong to the rigidbodies
        for (rb = mol->beginRigidBody(rbIter); rb != NULL; rb = mol->nextRigidBody(rbIter)) {
            rb->updateAtoms();
        }
        
    }      

      /* Calculate "bonds" and build Q_lm(r) where Q_lm = Y_lm(theta(r),phi(r)) */
      for (std::vector<std::pair<StuntDouble*, StuntDouble*> >::iterator j = sdPairs_.begin(); j != sdPairs_.end(); ++j) {
          Vector3d vec = j->first->getPos() - j->second->getPos();
          currentSnapshot_->wrapVector(vec);
 /* The spherical harmonics are wrt any arbitray coordiate sysetm, 
  * we choose standard spherical coordinates */
  		r = sqrt(pow(vec.x(),2)+pow(vec.y(),2)+pow(vec.z(),2));
  		
  		/* Check to see if neighbor is in bond cuttoff*/
  		if (r<rCut_){  		
  	    theta = atan(vec.y()/vec.x());
  	    phi = acos(vec.z()/r);
  	    for(int m = -lNumber_; m <= lNumber_; m++){
  	    	sphericalHarmonic.setM(m);
  	    	QBar_lm(m) += sphericalHarmonic.getValueAt(theta,phi);
  	    }
  	    nBonds++;
  		}  
      }
      
     /*Normalize Qbar*/
     for (int m = -lNumber_;m <= lNumber_; m++){
      QBar_lm(m) = QBar_lm(m)/nBonds;	
     }
      
    
  }
  
  /*Normalize by number of frames*/
    for (int m = -lNumber_;m <= lNumber_; m++){
      QBar_lm(m) = QBar_lm(m)/nFrames;	
     }
     
     
     
     /* Find second order invariant Q_l*/
     
       for (int m = -lNumber_;m <= lNumber_; m++){
     	 QSq_l += pow(QBar_lm(m),2);	
     	}
     Q_l = sqrt((4*NumericConstant::PI/lNumber_+1)*QSq_l);
     
     /* Find Third Order Invariant W_l*/
     for (int m = -lNumber_;m<= lNumber_;m++){
     	
     	
     }
     
     
  writeOrderParameter();
  
}

void BondOrderParameter::initalizeHistogram() {
    std::fill(histogram_.begin(), histogram_.end(), 0);
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

  void BondOrderParameter::collectHistogram(StuntDouble* sd1, StuntDouble* sd2) {

    if (sd1 == sd2) {
      return;
    }
    
    Vector3d pos1 = sd1->getPos();
    Vector3d pos2 = sd2->getPos();
    Vector3d r12 = pos2 - pos1;
    currentSnapshot_->wrapVector(r12);

    RealType distance = r12.length();

    if (distance < len_) {
      int whichBin = distance / deltaR_;
      histogram_[whichBin] += 2;
    }
  }






void BondOrderParameter::writeOrderParameter() {

    std::ofstream os(getOutputFileName().c_str());
    os << "#radial distribution function\n";
    os<< "#selection1: (" << selectionScript1_ << ")\t";
    os << "selection2: (" << selectionScript2_ << ")\n";
    os << "#p2\tdirector_x\tdirector_y\tdiretor_z\tangle(degree)\n";    

    for (std::size_t i = 0; i < orderParams_.size(); ++i) {
        os <<  orderParams_[i].p2 << "\t"
            <<  orderParams_[i].director[0] << "\t"
            <<  orderParams_[i].director[1] << "\t"
            <<  orderParams_[i].director[2] << "\t"
            <<  orderParams_[i].angle << "\n";

    }

}

}

