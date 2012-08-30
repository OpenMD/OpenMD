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
 *  @version $Id: BondOrderParameter.cpp 1442 2010-05-10 17:28:26Z gezelter $
 *
 */
 
#include "applications/staticProps/TetrahedralityParam.hpp"
#include "utils/simError.h"
#include "io/DumpReader.hpp"
#include "primitives/Molecule.hpp"
#include "utils/NumericConstant.hpp"
#include <vector>

namespace OpenMD {

  TetrahedralityParam::TetrahedralityParam(SimInfo* info, 
                                           const std::string& filename, 
                                           const std::string& sele,
                                           double rCut, int nbins) : StaticAnalyser(info, filename), selectionScript_(sele), evaluator_(info), seleMan_(info){
    
    setOutputName(getPrefix(filename) + ".q");

    evaluator_.loadScriptString(sele);
    if (!evaluator_.isDynamic()) {
      seleMan_.setSelectionSet(evaluator_.evaluate());
    }

    // Set up cutoff radius:

    rCut_ = rCut;
    nBins_ = nbins;

    Q_histogram_.resize(nBins_);

    // Q can take values from 0 to 1

    MinQ_ = 0.0;
    MaxQ_ = 1.1;
    deltaQ_ = (MaxQ_ - MinQ_) / nbins;

  }

  TetrahedralityParam::~TetrahedralityParam() {
    Q_histogram_.clear(); 
  }
  
  void TetrahedralityParam::initializeHistogram() {
    std::fill(Q_histogram_.begin(), Q_histogram_.end(), 0);
  }
  
  void TetrahedralityParam::process() {
    Molecule* mol;
    StuntDouble* sd;
    StuntDouble* sd2;
    StuntDouble* sdi;
    StuntDouble* sdj;
    RigidBody* rb;
    int myIndex;
    SimInfo::MoleculeIterator mi;
    Molecule::RigidBodyIterator rbIter;
    Molecule::IntegrableObjectIterator ioi;
    Vector3d vec;
    Vector3d ri, rj, rk, rik, rkj, dposition, tposition;
    RealType r;
    RealType cospsi;
    RealType Qk;
    std::vector<std::pair<RealType,StuntDouble*> > myNeighbors;
    int isd;

    DumpReader reader(info_, dumpFilename_);    
    int nFrames = reader.getNFrames();
    frameCounter_ = 0;

    Distorted_.clear();
    Tetrahedral_.clear();

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

      for (sd = seleMan_.beginSelected(isd); sd != NULL; 
           sd = seleMan_.nextSelected(isd)) {
	
        myIndex = sd->getGlobalIndex();
	Qk = 1.0;

	myNeighbors.clear();
                
        // inner loop is over all StuntDoubles in the system:
        
        for (mol = info_->beginMolecule(mi); mol != NULL; 
             mol = info_->nextMolecule(mi)) {

          for (sd2 = mol->beginIntegrableObject(ioi); sd2 != NULL; 
               sd2 = mol->nextIntegrableObject(ioi)) {
	    
            if (sd2->getGlobalIndex() != myIndex) {
	      
              vec = sd->getPos() - sd2->getPos();       
	      
              if (usePeriodicBoundaryConditions_) 
                currentSnapshot_->wrapVector(vec);
	      
              r = vec.length();             

              // Check to see if neighbor is in bond cutoff 
              
              if (r < rCut_) { 
		
		myNeighbors.push_back(std::make_pair(r,sd2));
	      }
	    }
	  }
	}

	// Sort the vector using predicate and std::sort
	std::sort(myNeighbors.begin(), myNeighbors.end());

	//std::cerr << myNeighbors.size() <<  " neighbors within " << rCut_  << " A" << " \n";
	
	// Use only the 4 closest neighbors to do the rest of the work:
	
	int nbors =  myNeighbors.size()> 4 ? 4 : myNeighbors.size();
	int nang = int (0.5 * (nbors * (nbors - 1)));

	rk = sd->getPos();
	//std::cerr<<nbors<<endl;
	for (int i = 0; i < nbors-1; i++) {	  

	  sdi = myNeighbors[i].second;
	  ri = sdi->getPos();
	  rik = rk - ri;
	  if (usePeriodicBoundaryConditions_) 
	    currentSnapshot_->wrapVector(rik);
	  
	  rik.normalize();

	  for (int j = i+1; j < nbors; j++) {	    

	    sdj = myNeighbors[j].second;
	    rj = sdj->getPos();
	    rkj = rk - rj;
	    if (usePeriodicBoundaryConditions_) 
	      currentSnapshot_->wrapVector(rkj);
	    rkj.normalize();
	    
	    cospsi = dot(rik,rkj);

	    //std::cerr << "cos(psi) = " << cospsi << " \n";

	    // Calculates scaled Qk for each molecule using calculated angles from 4 or fewer nearest neighbors.
	    Qk = Qk - (pow(cospsi + 1.0 / 3.0, 2) * 2.25 / nang);
	    //std::cerr<<Qk<<"\t"<<nang<<endl;
	  }
	}
	//std::cerr<<nang<<endl;
	if (nang > 0) {
	  collectHistogram(Qk);

	// Saves positions of StuntDoubles & neighbors with distorted coordination (low Qk value)
	if ((Qk < 0.55) && (Qk > 0.45)) {
	  //std::cerr<<Distorted_.size()<<endl;
	  Distorted_.push_back(sd);
	  //std::cerr<<Distorted_.size()<<endl;
	  dposition = sd->getPos();
	  //std::cerr << "distorted position \t" << dposition << "\n";
	}

	// Saves positions of StuntDoubles & neighbors with tetrahedral coordination (high Qk value)
	if (Qk > 0.05) { 

	  Tetrahedral_.push_back(sd);

	  tposition = sd->getPos();
	  //std::cerr << "tetrahedral position \t" << tposition << "\n";
	}

	//std::cerr<<Tetrahedral_.size()<<endl;


	}

      }
    }
    
    writeOrderParameter();
    std::cerr << "number of distorted StuntDoubles = " << Distorted_.size() << "\n";
    std::cerr << "number of tetrahedral StuntDoubles = " << Tetrahedral_.size() << "\n";
  }
        
  void TetrahedralityParam::collectHistogram(RealType Qk) {

    if (Qk > MinQ_ && Qk < MaxQ_) {
     
      int whichBin = int((Qk - MinQ_) / deltaQ_);
      Q_histogram_[whichBin] += 1;
    }
  }    


  void TetrahedralityParam::writeOrderParameter() {
    
    int nSelected = 0;

    for (int i = 0; i < nBins_; ++i) {
      nSelected = nSelected + int(Q_histogram_[i]*deltaQ_);
    }
    
    std::ofstream osq((getOutputFileName() + "Q").c_str());

    if (osq.is_open()) {
      
      osq << "# Tetrahedrality Parameters\n";
      osq << "# selection: (" << selectionScript_ << ")\n";
      osq << "# \n";
      // Normalize by number of frames and write it out:
      for (int i = 0; i < nBins_; ++i) {
        RealType Qval = MinQ_ + (i + 0.5) * deltaQ_;               
        osq << Qval;
	osq << "\t" << (RealType) (Q_histogram_[i]/deltaQ_)/nSelected;        
        osq << "\n";
      }

      osq.close();
      
    }else {
      sprintf(painCave.errMsg, "TetrahedralityParam: unable to open %s\n", 
              (getOutputFileName() + "q").c_str());
      painCave.isFatal = 1;
      simError();  
    }

    DumpReader reader(info_, dumpFilename_);    
    int nFrames = reader.getNFrames();

    if (nFrames == 1) {

    std::vector<StuntDouble*>::iterator iter;
    std::ofstream osd((getOutputFileName() + "dxyz").c_str());

    if (osd.is_open()) {

      osd << Distorted_.size() << "\n";

      osd << "1000000.00000000;    34.52893134     0.00000000     0.00000000;     0.00000000    34.52893134     0.00000000;     0.00000000     0.00000000    34.52893134" << "\n";
      
      for (iter = Distorted_.begin(); iter != Distorted_.end(); ++iter) {
	Vector3d position;
	position = (*iter)->getPos();
	osd << "O  " << "\t";
	  for (unsigned int z = 0; z < position.size(); z++) {
	    osd << position[z] << "  " << "\t";
	  }
	  osd << "\n";
      }
      osd.close();
    }


    std::ofstream ost((getOutputFileName() + "txyz").c_str());
    
    if (ost.is_open()) {

      ost << Tetrahedral_.size() << "\n";

      ost << "1000000.00000000;    34.52893134     0.00000000     0.00000000;     0.00000000    34.52893134     0.00000000;     0.00000000     0.00000000    34.52893134" << "\n";
      
      for (iter = Tetrahedral_.begin(); iter != Tetrahedral_.end(); ++iter) {

	Vector3d position;

	position = (*iter)->getPos();

	ost << "O  " << "\t";

	  for (unsigned int z = 0; z < position.size(); z++) {

	    ost << position[z] << "  " << "\t";
	  }

	  ost << "\n";

      }

      ost.close();
    }
    
  }
}
}
      


