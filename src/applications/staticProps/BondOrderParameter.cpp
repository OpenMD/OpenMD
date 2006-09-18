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
      for (sd2 = seleMan1_.beginSelected(j),sd2
	     sdPairs_.push_back(std::make_pair(sd1, sd2));
	   }


    }

    void BondOrderParameter::process
      () {
      Molecule* mol;
      RigidBody* rb;
      SimInfo::MoleculeIterator mi;
      Molecule::RigidBodyIterator rbIter;
      RealType theta;
      RealType phi;
      RealType r;
      RealType dist;
      RealType* QBar_lm;
      RealType QSq_l;
      int nBonds;
      int m, m_index;
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


	/* Setup QBar */
	QBar_lm = new double[mSize_];

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
	    for(int m_index = 0; m_index < mSize_; m_index++){
	      sphericalHarmonic.setM(m_index-lNumber_);
	      QBar_lm(m_index) += sphericalHarmonic.getValueAt(theta,phi);
	    }
	    nBonds++;
	  }
	}

	/*Normalize Qbar by number of Bonds*/
	for ( int m_index = 0;m_index < mSize_; m_index++){
	  QBar_lm(m_index) = QBar_lm(m_index)/nBonds;
	}


      }

      /*Normalize by number of frames*/
      for ( int m_index = 0;m_index < mSize_; m_index++){
	QBar_lm(m_index) = QBar_lm(m_index)/nFrames;
      }



      /* Find second order invariant Q_l*/

      for (int m_index = 0 ;m_index <= sizeM_; m++){
	QSq_l += pow(QBar_lm(m),2);
      }
      Q_l_ = sqrt((4*NumericConstant::PI/lNumber_+1)*QSq_l);

      /* Find Third Order Invariant W_l*/

      /* Make arrays for Wigner3jm */
      double* THRCOF = new double[mSize_];
      /* Variables for Wigner routine */
      double l_ = (double)lNumber_;
      double m2Min;
      double m2Max;
      int error;
      int m1;
      int m2;
      int m3;

      for (int m1 = -lNumber_;m <= lNumber_;m1++){
	/* Zero work array */
	for (i=0; i<mSize_;i++){
	  THRCOF[i] = 0.0;      
	}
	/* Get wigner coefficients */
	Wigner3jm(&l_,&l_,&l_,&(double)m1,&m2Min,&m2Max&,THRCOF,&mSize_,&error);
	for (m_index=1; i<(m2Max-M2Min-1.0);m_index++){
	  m2 = floor(m2Min) + m_index - 1;
	  m3 = -m1-m2;
	  W_l_ += THRCOF(m_index)*QBar_lm(m1+lNumber_)*QBar_lm(m2+lNumber_)*QBar_lm(m3+lNumber_);
	}
      }


      writeOrderParameter();

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

