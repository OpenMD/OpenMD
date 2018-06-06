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
 * [4]  Kuang & Gezelter,  J. Chem. Phys. 133, 164101 (2010).
 * [5]  Vardeman, Stocker & Gezelter, J. Chem. Theory Comput. 7, 834 (2011).
 * [6]  Kuang & Gezelter, Mol. Phys., 110, 691-701 (2012).
 */
 
#include "analysis/TetrahedralityParamZ.hpp"
#include "utils/simError.h"
#include "io/DumpReader.hpp"
#include "primitives/Molecule.hpp"
#include <vector>
#include <algorithm>
#include <fstream>

using namespace std;
namespace OpenMD {
  TetrahedralityParamZ::TetrahedralityParamZ(SimInfo* info) 
    : ObjectAnalyzer(info) {

    string prefixFileName = info->getPrefixFileName();
    setOutputName(prefixFileName + ".Qz");
    
    evaluator1_.loadScriptString(sele1);
    if (!evaluator1_.isDynamic()) {
      seleMan1_.setSelectionSet(evaluator1_.evaluate());
    }
    evaluator2_.loadScriptString(sele2);
    if (!evaluator2_.isDynamic()) {
      seleMan2_.setSelectionSet(evaluator2_.evaluate());
    }
    

    tetrahedrality = new OutputData;
    tetrahedrality->units =  "Unitless";
    tetrahedrality->title =  "Tetrahedrality";
    tetrahedrality->dataType = odtReal;
    tetrahedrality->dataHandling = odhAverage;
    tetrahedrality->accumulator.reserve(nBins_);
    for (unsigned int i = 0; i < nBins_; i++) 
      tetrahedrality->accumulator.push_back( new Accumulator() );
    data_.push_back(tetrahedrality);

  }
  
  TetrahedralityParamZ::~TetrahedralityParamZ() {
  }

    
  void TetrahedralityParamZ::processFrame(int istep) {
    StuntDouble* sd;
    StuntDouble* sd2;
    StuntDouble* sdi;
    StuntDouble* sdj;
    int myIndex;
    Vector3d vec;
    Vector3d ri, rj, rk, rik, rkj;
    RealType r;
    RealType cospsi;
    RealType Qk;
    std::vector<std::pair<RealType,StuntDouble*> > myNeighbors;
    int isd1;
    int isd2;
    bool usePeriodicBoundaryConditions_ = info_->getSimParams()->getUsePeriodicBoundaryConditions();

    RealType z;
    hmat_ = currentSnapshot_->getHmat();
    for (unsigned int i = 0; i < nBins_; i++) {
      z = (((RealType)i + 0.5) / (RealType)nBins_) * hmat_(axis_,axis_);
      dynamic_cast<Accumulator*>(z_->accumulator[i])->add(z);
    }
    
    vector<RealType> zBox(nBins_, 0.0);
    vector<RealType> sliceQ(nBins_, 0.0);
    vector<int> sliceCount(nBins_, 0);

      
    zBox.push_back(hmat_(axis_,axis_));
    
    RealType halfBoxZ = hmat_(axis_,axis_) / 2.0;      
    
    if (evaluator1_.isDynamic()) {
      seleMan1_.setSelectionSet(evaluator1_.evaluate());
    }
    
    if (evaluator2_.isDynamic()) {
      seleMan2_.setSelectionSet(evaluator2_.evaluate());
    }
    // outer loop is over the selected StuntDoubles:
    for (sd = seleMan1_.beginSelected(isd1); sd != NULL;
	 sd = seleMan1_.nextSelected(isd1)) {
      
      myIndex = sd->getGlobalIndex();
      
      Qk = 1.0;	  
      myNeighbors.clear();       
      
      for (sd2 = seleMan2_.beginSelected(isd2); sd2 != NULL;
	   sd2 = seleMan2_.nextSelected(isd2)) {
	
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
      
      // Sort the vector using predicate and std::sort
      std::sort(myNeighbors.begin(), myNeighbors.end());
      
      // Use only the 4 closest neighbors to do the rest of the work:
      
      //int nbors =  myNeighbors.size()> 4 ? 4 : myNeighbors.size();
      int nbors = myNeighbors.size();
      int nang = int (0.5 * (nbors * (nbors - 1)));
      
      rk = sd->getPos();
      
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
          
	  // Calculates scaled Qk for each molecule using calculated
	  // angles from 4 or fewer nearest neighbors.
	  Qk = Qk - (pow(cospsi + 1.0 / 3.0, 2) * 9.0 / (4.0 * nang) );         
	}
      }
      
      if (nang > 0) {
	if (usePeriodicBoundaryConditions_)
	  currentSnapshot_->wrapVector(rk);
	
	int binNo = int(nBins_ * (halfBoxZ + rk[axis_]) / hmat_(axis_,axis_));
	sliceQ[binNo] += Qk;
	sliceCount[binNo] += 1;
      }  
    }

    
    for (unsigned int i = 0; i < nBins_; i++) {
      if (sliceCount[i] > 0) {
        RealType tetra = sliceQ[i] / sliceCount[i];
        dynamic_cast<Accumulator *>(tetrahedrality->accumulator[i])->add(tetra);
	dynamic_cast<Accumulator *>(counts_->accumulator[i])->add(1);
      }
    }
  }//processFrame

  void TetrahedralityParamZ::processStuntDouble(StuntDouble* sd, int bin) {
    // Fill in later
  }
  

}



