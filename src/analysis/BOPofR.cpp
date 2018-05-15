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
 * [4] , Stocker & Gezelter, J. Chem. Theory Comput. 7, 834 (2011). *
 *  Created by J. Daniel Gezelter on 09/26/06.
 *  @author  J. Daniel Gezelter
 *  @version $Id$
 *
 */
 
#include "analysis/BOPofR.hpp"
#include "utils/simError.h"
#include "utils/Revision.hpp"
#include "io/DumpReader.hpp"
#include "primitives/Molecule.hpp"
#include "utils/Constants.hpp"
#include "math/Wigner3jm.hpp"
#include "brains/Thermo.hpp"

using namespace MATPACK;
namespace OpenMD {

  BOPofR::BOPofR(SimInfo* info) : NonSpatialStatistics(info, sele, nbins),
				  SingletType() {
    
    string prefixFileName = info_->getPrefixFileName();
    setOutputName(prefixFileName + ".bo");
    setAnalysisType("Bond Order Parameter(r)");

    // need to deal with evaluator_
    evaluator_.loadScriptString(sele);
    if (!evaluator_.isDynamic()) {
      seleMan_.setSelectionSet(evaluator_.evaluate());
    }
    
    std::stringstream params;
    params << " rcut = " << rCut_
           << ", len = " << len_
           << ", nbins = " << nBins_;
    const std::string paramString = params.str();
    setParameterString( paramString );


    deltaR_ = len_/nBins_;

    WofR = new OutputData;
    WofR->units =  "Unitless";
    WofR->title =  "OrderParam";
    WofR->dataType = odtReal;
    WofR->dataHandling = odhAverage;
    WofR->accumulator.reserve(nBins_);
    for (unsigned int i = 0; i < nBins_; i++) 
      WofR->accumulator.push_back( new Accumulator() );
    data_.push_back(WofR);

    QofR = new OutputData;
    QofR->units =  "Unitless";
    QofR->title =  "OrderParam";
    QofR->dataType = odtReal;
    QofR->dataHandling = odhAverage;
    QofR->accumulator.reserve(nBins_);
    for (unsigned int i = 0; i < nBins_; i++) 
      QofR->accumulator.push_back( new Accumulator() );
    data_.push_back(QofR);

        
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
    
    q_l_.resize(lMax_+1);
    q2_.resize(lMax_+1);
    w_.resize(lMax_+1);
    w_hat_.resize(lMax_+1);

    Q2_.resize(lMax_+1);
    Q_.resize(lMax_+1);
    W_.resize(lMax_+1);
    W_hat_.resize(lMax_+1);

    thermo_ = new Thermo(info_);
    frameCounter_ = 0;
    
  }
  
  BOPofR::~BOPofR() {
    delete thermo_;
  /*
    std::cerr << "Freeing stuff" << std::endl;
    for (int l = 0; l <= lMax_; l++) {
      for (int m = -l; m <= l; m++) {
        w3j[std::make_pair(l,m)].clear();
      }
    }
	std::cerr << "w3j made free...." << std::endl;
   for (int bin = 0; bin < nBins_; bin++) {
		QofR_[bin].clear();
		WofR_[bin].clear();
		RCount_[bin].clear();
    }
	std::cout << "R arrays made free...." << std::endl;
   w3j.clear();
    m2Min.clear();
    m2Max.clear();
    RCount_.clear();
    WofR_.clear();
    QofR_.clear();
 */
  }

  
  
  void BOPofR::processFrame(int istep) {
    Molecule* mol;
    Atom* atom;
    int myIndex;
    SimInfo::MoleculeIterator mi;
    Molecule::AtomIterator ai;
    StuntDouble* sd;
    Vector3d vec;
    RealType costheta;
    RealType phi;
    RealType r;
    Vector3d rCOM;
    RealType distCOM;
    Vector3d pos;
    Vector3d CenterOfMass;
    int nBonds;
    SphericalHarmonic sphericalHarmonic;
    int i;
    bool usePeriodicBoundaryConditions_ = info_->getSimParams()->getUsePeriodicBoundaryConditions();

    
    
    CenterOfMass = thermo_->getCom();
    if (evaluator_.isDynamic()) {
      seleMan_.setSelectionSet(evaluator_.evaluate());
    }

    // outer loop is over the selected StuntDoubles:
    
    for (sd = seleMan_.beginSelected(i); sd != NULL; 
	 sd = seleMan_.nextSelected(i)) {
      
      myIndex = sd->getGlobalIndex();
      
      nBonds = 0;
        
      for (int l = 0; l <= lMax_; l++) {
	for (int m = -l; m <= l; m++) {
	  q_[std::make_pair(l,m)] = 0.0;
	}
      }
      pos = sd->getPos();
      rCOM = CenterOfMass - pos;
      if (usePeriodicBoundaryConditions_) 
	currentSnapshot_->wrapVector(rCOM);
      distCOM = rCOM.length();
      
      // inner loop is over all other atoms in the system:
      
      for (mol = info_->beginMolecule(mi); mol != NULL; 
	   mol = info_->nextMolecule(mi)) {
	for (atom = mol->beginAtom(ai); atom != NULL; 
	     atom = mol->nextAtom(ai)) {
	  
	  if (atom->getGlobalIndex() != myIndex) {
	    vec = pos - atom->getPos();       
	    
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
		  q_[std::make_pair(l,m)] += sphericalHarmonic.getValueAt(costheta, phi);
		}
	      }
	      nBonds++;
	    }  
	  }
	}
      }
      
      
      for (int l = 0; l <= lMax_; l++) {
	q2_[l] = 0.0;
	for (int m = -l; m <= l; m++){
	  q_[std::make_pair(l,m)] /= (RealType)nBonds;            
	  q2_[l] += norm(q_[std::make_pair(l,m)]);
	}
	q_l_[l] = sqrt(q2_[l] * 4.0 * Constants::PI / (RealType)(2*l + 1));
      }
      
      // Find Third Order Invariant W_l
      
      for (int l = 0; l <= lMax_; l++) {
	w_[l] = 0.0;
	for (int m1 = -l; m1 <= l; m1++) {
	  std::pair<int,int> lm = std::make_pair(l, m1);
	  for (int mmm = 0; mmm <= (m2Max[lm] - m2Min[lm]); mmm++) {
	    int m2 = m2Min[lm] + mmm;
	    int m3 = -m1-m2;
	    w_[l] += w3j[lm][mmm] * q_[lm] * 
	      q_[std::make_pair(l,m2)] *  q_[std::make_pair(l,m3)];
	  }
	}
	
	w_hat_[l] = w_[l] / pow(q2_[l], RealType(1.5));
      }
      
      collectHistogram(q_l_, w_hat_, distCOM);
      
      //  printf( "%s  %18.10g %18.10g %18.10g %18.10g \n", sd->getType().c_str(),pos[0],pos[1],pos[2],real(w_hat[6]));
      
    }
  }
  

  void BOPofR::processStuntDouble(StuntDouble* sd, int bin){
    // fill in later
  }

  IcosahedralOfR::IcosahedralOfR(SimInfo* info, 
                                 const std::string& sele, double rCut, 
                                 unsigned int nbins, RealType len) :
    BOPofR(info, sele, rCut, nbins, len) {
    
    setAnalysisType("Icosahedral Bond Order Parameter(r)");
  }

  void IcosahedralOfR::collectHistogram(std::vector<RealType> q, 
                                        std::vector<ComplexType> what, 
                                        RealType distCOM) {
    if ( distCOM < len_){
      // Figure out where this distance goes...
      int whichBin = int(distCOM / deltaR_);
      dynamic_cast<Accumulator* >(counts_->accumulator[whichBin])->add(1);

      if(real(what[6]) < -0.15){
	dynamic_cast<Accumulator* >(WofR->accumulator[whichBin])->add(1);
      }
      if(q[6] > 0.5){
	dynamic_cast<Accumulator* >(QofR->accumulator[whichBin])->add(1);
      }
    }      
  }

  FCCOfR::FCCOfR(SimInfo* info, 
                 const std::string& sele, double rCut, 
                 unsigned int nbins, RealType len) :
    BOPofR(info) {
    setAnalysisType("FCC Bond Order Parameter(r)");
  }
  
  
  void FCCOfR::collectHistogram(std::vector<RealType> q, 
                                std::vector<ComplexType> what, 
                                RealType distCOM) {
    
    if ( distCOM < len_){
      // Figure out where this distance goes...
      int whichBin = int(distCOM / deltaR_);
      dynamic_cast<Accumulator* >(counts_->accumulator[whichBin])->add(1);

      if(real(what[4]) < -0.12){
	dynamic_cast<Accumulator* >(WofR->accumulator[whichBin])->add(1);
      }
    }      
  }

  void FCCOfR::processStuntDouble(StuntDouble* sd, int bin) {
    // Fill in later
  }

  void IcosahedralOfR::processStuntDouble(StuntDouble* sd, int bin) {
    // Fill in later
  }
  
  void IcosahedralOfR::writeOutput() {
    
    Revision rev; 
    std::ofstream osq((getOutputFileName() + "qr").c_str());

    if (osq.is_open()) {
      osq << "# " << getAnalysisType() << "\n";
      osq << "# OpenMD " << rev.getFullRevision() << "\n";
      osq << "# " << rev.getBuildDate() << "\n";
      osq << "# selection script: \"" << selectionScript_  << "\"\n";
      if (!paramString_.empty())
        osq << "# parameters: " << paramString_ << "\n";
      osq << "#r\tQ(r)\tW(r)\n";
      
      // Normalize by number of frames and write it out:

      int nWR;
      int nQR;
      int nSele;
      RealType aveWR;
      RealType aveQR;
      
      for (unsigned int j = 0; j < nBins_; j++) {        
       	  
	nWR = WofR->accumulator[j]->count();
	nQR = QofR->accumulator[j]->count();

	nSele = counts_->accumulator[j]->count();

	aveWR = (RealType)nWR / nSele;
	aveQR = (RealType)nQR / nSele;

	RealType Rval = (j + 0.5) * deltaR_;

	osq << Rval;
        if (nSele == 0){
          osq << "\t" << 0;
          osq << "\n";
        }else{
          osq << "\t" << aveQR << "\t" << aveWR;        
          osq << "\n";
	}
      }
      
      osq.close();
      
    } else {
      sprintf(painCave.errMsg, "IcosahedralOfR: unable to open %s\n", 
              (getOutputFileName() + "q").c_str());
      painCave.isFatal = 1;
      simError();  
    }    
    
  }
  
  void FCCOfR::writeOutput() {
    
    Revision rev; 
    std::ofstream osq((getOutputFileName() + "qr").c_str());

    if (osq.is_open()) {
      osq << "# " << getAnalysisType() << "\n";
      osq << "# OpenMD " << rev.getFullRevision() << "\n";
      osq << "# " << rev.getBuildDate() << "\n";
      osq << "# selection script: \"" << selectionScript_  << "\"\n";
      if (!paramString_.empty())
        osq << "# parameters: " << paramString_ << "\n";
      osq << "#r\tW(r)\n";
      
      // Normalize by number of frames and write it out:

      int nWR;
      int nSele;
      RealType aveWR;
      
      for (unsigned int j = 0; j < nBins_; j++) {        
       	  
	nWR = WofR->accumulator[j]->count();

	nSele = counts_->accumulator[j]->count();

	aveWR = (RealType)nWR / nSele;

	RealType Rval = (j + 0.5) * deltaR_;

	osq << Rval;
        if (nSele == 0){
          osq << "\t" << 0;
          osq << "\n";
        }else{
          osq << "\t" << "\t" << aveWR;        
          osq << "\n";
	}
      }
      
      osq.close();
      
    } else {
      sprintf(painCave.errMsg, "IcosahedralOfR: unable to open %s\n", 
              (getOutputFileName() + "q").c_str());
      painCave.isFatal = 1;
      simError();  
    }
	
  }
}
