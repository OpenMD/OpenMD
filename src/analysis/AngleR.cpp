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
 */

/* Calculates Angle(R) for DirectionalAtoms*/

#include <algorithm>
#include <fstream>
#include <sstream>
#include <cmath>

#include "analysis/AngleR.hpp"
#include "utils/simError.h"
#include "utils/Revision.hpp"
#include "io/DumpReader.hpp"
#include "primitives/Molecule.hpp"
#include "brains/Thermo.hpp"
#include "utils/Revision.hpp"

namespace OpenMD {

  AngleR::AngleR(SimInfo* info,
                 const std::string& sele, RealType len, int nrbins)
    : NonSpatialStatistics(info, sele, nrbins) : SingletType() {
    
    evaluator_.loadScriptString(sele);
    if (!evaluator_.isDynamic()) {
      seleMan_.setSelectionSet(evaluator_.evaluate());
    }
        
    deltaR_ = len_ /nRBins_;
    
    setAnalysisType("radial density function Angle(r)");
    string prefixFileName = info_->getPrefixFileName();
    setOutputName(prefixFileName + ".AngleR");
    std::stringstream params;
    params << " len = " << len_
           << ", nrbins = " << nRBins_;

    angleR = new OutputData;
    angleR->units = "Unitless";
    angleR->title = "Angle R";
    angleR->dataType = odtReal;
    angleR->dataHandling = odhAverage;
    angleR->accumulator.reserve(nRBins_);
    for (unsigned int i = 0; i < nRBins_; i++) 
      angleR->accumulator.push_back( new Accumulator() );
    data_.push_back(angleR);


  }

  AngleR::~AngleR(){
  }
  
  
  void AngleR::processFrame(int istep) {
    
    int i;    
    StuntDouble* sd;
    Thermo thermo(info_);
    Vector3d CenterOfMass = thermo.getCom();      
    
    if (evaluator_.isDynamic()) {
      seleMan_.setSelectionSet(evaluator_.evaluate());
    }
    
    //determine which atom belongs to which slice
    for (sd = seleMan_.beginSelected(i); sd != NULL;
	 sd = seleMan_.nextSelected(i)) {
      Vector3d pos = sd->getPos();
      Vector3d r1 = CenterOfMass - pos;
      // only do this if the stunt double actually has a vector associated
      // with it
      if (sd->isDirectional()) {
	Vector3d uz = sd->getA().transpose() * V3Z;
	// std::cerr << "pos = " << pos << " uz = " << uz << "\n";
	RealType distance = r1.length();
	
	uz.normalize();
	r1.normalize();
	RealType cosangle = dot(r1, uz);
	
	if (distance < len_) {
	  int whichBin = int(distance / deltaR_);
	  // update accumulators here
	  dynamic_cast<Accumulator *>(angleR->accumulator[whichBin])->add(cosangle);
	  
	}
      } 
    }
  }
  
  void AngleR::processStuntDouble(StuntDouble* sd, int bin) {
    // Fill in later
  }

  void AngleR::writeOutput() {
    std::ofstream ofs(outputFilename_.c_str());
    if (ofs.is_open()) {
      
      Revision rev;
      
      ofs << "# " << getAnalysisType() << "\n";
      ofs << "# OpenMD " << rev.getFullRevision() << "\n";
      ofs << "# " << rev.getBuildDate() << "\n";
      ofs << "# selection script: \"" << selectionScript_  << "\"\n";
      if (!paramString_.empty())
        ofs << "# parameters: " << paramString_ << "\n";

      RealType aveAngleR = 0.0;	
      
      for (unsigned int j = 0; j < nRBins_; j++) {        
	      	
	dynamic_cast<Accumulator*>(angleR->accumulator[j])->getAverage(aveAngleR);	  
	  
	RealType r = deltaR_ * (j + 0.5);
	ofs << r << "\t" << aveAngleR << "\n";
      }
        
    } else {
	
      sprintf(painCave.errMsg, "AngleR: unable to open %s\n",
	      outputFilename_.c_str());
      painCave.isFatal = 1;
      simError();  
    }
    
    ofs.close();
  }
  
}

