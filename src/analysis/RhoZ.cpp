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
 *  Created by Charles F. Vardeman II on 11/26/05.
 *  @author  Charles F. Vardeman II 
 *  @version $Id$
 *
 */

/* Calculates Rho(Z) for density profile of liquid slab. */

#include <algorithm>
#include <fstream>
#include "analysis/RhoZ.hpp"
#include "utils/simError.h"
#include "io/DumpReader.hpp"
#include "primitives/Molecule.hpp"
namespace OpenMD {
  
  RhoZ::RhoZ(SimInfo* info, 
	     const std::string& sele, int nzbins, int axis)
    : SlabStatistics(info, sele, nzbins, axis), selectionScript_(sele), 
      evaluator_(info), seleMan_(info), axis_(axis) {

    string prefixFileName = info->getPrefixFileName();
    setOutputName(prefixFileName + ".RhoZ");
    
    evaluator_.loadScriptString(sele);
    if (!evaluator_.isDynamic()) {
      seleMan_.setSelectionSet(evaluator_.evaluate());
    }       

    density = new OutputData;
    density->units =  "g cm^-3";
    density->title =  "Density";
    density->dataType = odtReal;
    density->dataHandling = odhAverage;
    density->accumulator.reserve(nBins_);
    for (unsigned int i = 0; i < nBins_; i++) 
      density->accumulator.push_back( new Accumulator() );
    data_.push_back(density);
 
    usePeriodicBoundaryConditions_ = 
      info_->getSimParams()->getUsePeriodicBoundaryConditions();
  }

  RhoZ::~RhoZ() {
  }
  
  void RhoZ::processFrame(int istep) {
    RealType z;

    
    hmat_ = currentSnapshot_->getHmat();
    for (unsigned int i = 0; i < nBins_; i++) {
      z = (((RealType)i + 0.5) / (RealType)nBins_) * hmat_(axis_,axis_);
      dynamic_cast<Accumulator*>(z_->accumulator[i])->add(z);
    }
    volume_ = currentSnapshot_->getVolume();
    
    StuntDouble* sd;
    int i;

    vector<RealType> binMass(nBins_, 0.0);
    vector<int> binDof(nBins_, 0);

    if (evaluator_.isDynamic()) {
      seleMan_.setSelectionSet(evaluator_.evaluate());
    } 

    for (sd = seleMan_.beginSelected(i); sd != NULL; 
         sd = seleMan_.nextSelected(i)) {
      
      Vector3d pos = sd->getPos(); 
      RealType m = sd->getMass();
      
      int bin = getBin(pos);      
      
      binMass[bin] += m;
      binDof[bin] += 3;
      
      if (sd->isDirectional()) {
	if (sd->isLinear()) {
	  binDof[bin] += 2;
	} else {
	  binDof[bin] += 3;
	}
      }
    }
    for (unsigned int i = 0; i < nBins_; i++) {
      
      if (binDof[i] > 0) {
        RealType den = binMass[i] * nBins_ * Constants::densityConvert 
          / volume_;
        
        dynamic_cast<Accumulator *>(density->accumulator[i])->add(den);
        dynamic_cast<Accumulator *>(counts_->accumulator[i])->add(1);
      }
    }    
  }

  void RhoZ::processStuntDouble(StuntDouble* sd, int bin) {
    // Fill in later
  }
        
}

