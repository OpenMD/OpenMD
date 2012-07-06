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
 * [4]  Kuang & Gezelter,  J. Chem. Phys. 133, 164101 (2010).
 * [5]  Vardeman, Stocker & Gezelter, J. Chem. Theory Comput. 7, 834 (2011).
 */

/* Calculates Angle(R) for DirectionalAtoms*/
#include <algorithm>
#include <fstream>
#include "applications/staticProps/AngleR.hpp"
#include "utils/simError.h"
#include "utils/NumericConstant.hpp"
#include "io/DumpReader.hpp"
#include "primitives/Molecule.hpp"
#include "brains/Thermo.hpp"
#include <math.h>

namespace OpenMD {

     
  AngleR::AngleR(SimInfo* info, const std::string& filename, const std::string& sele, RealType len, int nrbins)
    : StaticAnalyser(info, filename), selectionScript_(sele),  evaluator_(info), seleMan_(info), len_(len), nRBins_(nrbins){
    
    
    evaluator_.loadScriptString(sele);
    if (!evaluator_.isDynamic()) {
      seleMan_.setSelectionSet(evaluator_.evaluate());
    }
    
    
    deltaR_ = len_ /nRBins_;
    
    histogram_.resize(nRBins_);
    count_.resize(nRBins_);
    avgAngleR_.resize(nRBins_);
    setOutputName(getPrefix(filename) + ".AngleR");
  }
  

  void AngleR::process() {
    
    DumpReader reader(info_, dumpFilename_);    
    int nFrames = reader.getNFrames();
    nProcessed_ = nFrames/step_;
    
    std::fill(avgAngleR_.begin(), avgAngleR_.end(), 0.0);     
    std::fill(histogram_.begin(), histogram_.end(), 0.0);
    std::fill(count_.begin(), count_.end(), 0);
    
    for (int istep = 0; istep < nFrames; istep += step_) {
      
      int i;    
      StuntDouble* sd;
      reader.readFrame(istep);
      currentSnapshot_ = info_->getSnapshotManager()->getCurrentSnapshot();
      Thermo thermo(info_);
      Vector3d CenterOfMass = thermo.getCom();      
      
      
      if (evaluator_.isDynamic()) {
	seleMan_.setSelectionSet(evaluator_.evaluate());
      }
      
      //determine which atom belongs to which slice
      for (sd = seleMan_.beginSelected(i); sd != NULL; sd = seleMan_.nextSelected(i)) {
	Vector3d pos = sd->getPos();
	Vector3d r1 = CenterOfMass - pos;
	// only do this if the stunt double actually has a vector associated
	// with it
	if (sd->isDirectional()) {
	  Vector3d dipole = sd->getA().getColumn(2);
	  // std::cerr << "pos = " << pos << " dipole = " << dipole << "\n";
	  RealType distance = r1.length();
	  
	  dipole.normalize();
	  r1.normalize();
	  RealType cosangle = dot(r1, dipole);
	  
	  if (distance < len_) {
	    int whichBin = distance / deltaR_;
	    histogram_[whichBin] += cosangle;
	    count_[whichBin] += 1;
	  }
	}
	
      }

    }

    processHistogram(); 
    writeAngleR();
  }



  void AngleR::processHistogram() {

    for(unsigned int i = 0 ; i < histogram_.size(); ++i){

      if (count_[i] > 0)
	avgAngleR_[i] += histogram_[i] / count_[i];    
      else 
	avgAngleR_[i] = 0.0;

    std::cerr << " count = " << count_[i] << " avgAngle = " << avgAngleR_[i] << "\n";  
    }

  }

 

  void AngleR::writeAngleR() {
    std::ofstream rdfStream(outputFilename_.c_str());
    if (rdfStream.is_open()) {
      rdfStream << "#radial density function Angle(r)\n";
      rdfStream << "#r\tcorrValue\n";
      for (unsigned int i = 0; i < avgAngleR_.size(); ++i) {
	RealType r = deltaR_ * (i + 0.5);
	rdfStream << r << "\t" << avgAngleR_[i] << "\n";
      }
        
    } else {

      sprintf(painCave.errMsg, "AngleR: unable to open %s\n", outputFilename_.c_str());
      painCave.isFatal = 1;
      simError();  
    }

    rdfStream.close();
  }

}

