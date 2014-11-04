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

#include <algorithm>
#include <functional>
#include "applications/sequentialProps/ContactAngle2.hpp"
#include "utils/simError.h"
#include "io/DumpReader.hpp"
#include "primitives/Molecule.hpp"
#include "utils/NumericConstant.hpp"
#include "utils/PhysicalConstants.hpp"
#include "math/Polynomial.hpp"

namespace OpenMD {

  ContactAngle2::ContactAngle2(SimInfo* info, const std::string& filename, 
                               const std::string& sele, RealType solidZ,
                               RealType threshDens, int nrbins, int nzbins)
    : SequentialAnalyzer(info, filename), selectionScript_(sele), 
      evaluator_(info), seleMan_(info), solidZ_(solidZ),
      threshDens_(threshDens), nRBins_(nrbins), nZBins_(nzbins) {
    
    setOutputName(getPrefix(filename) + ".ca2");
    
    evaluator_.loadScriptString(sele);
    
    if (!evaluator_.isDynamic()) {
      seleMan_.setSelectionSet(evaluator_.evaluate());
    }            
  }

  void ContactAngle2::doFrame() {
    StuntDouble* sd;
    int i;

    // set up the bins for density analysis

    Mat3x3d hmat = info_->getSnapshotManager()->getCurrentSnapshot()->getHmat();
    RealType len = std::min(hmat(0, 0), hmat(1, 1));
    RealType zLen = hmat(2,2);
    RealType dr = len / (RealType) nRBins_;
    RealType dz = zLen / (RealType) nZBins_;

    std::vector<std::vector<RealType> > histo;
    histo.resize(nRBins_);
    for (int i = 0 ; i < nRBins_; ++i) {
      histo[i].resize(nZBins_);
    }
    for (unsigned int i = 0; i < histo.size(); ++i){
      std::fill(histo[i].begin(), histo[i].end(), 0.0);
    }      
        
    if (evaluator_.isDynamic()) {
      seleMan_.setSelectionSet(evaluator_.evaluate());
    }
    

    RealType mtot = 0.0;
    Vector3d com(V3Zero);
    RealType mass;
    
    for (sd = seleMan_.beginSelected(i); sd != NULL;
         sd = seleMan_.nextSelected(i)) {      
      mass = sd->getMass();
      mtot += mass;
      com += sd->getPos() * mass;
    }

    com /= mtot;

    // now that we have the centroid, we can make cylindrical density maps
    Vector3d pos;
    RealType r;
    RealType z;
    
    for (sd = seleMan_.beginSelected(i); sd != NULL;
         sd = seleMan_.nextSelected(i)) {      
      pos = sd->getPos() - com;
      
      r = sqrt(pow(pos.x(), 2) + pow(pos.y(), 2));
      z = pos.z() - solidZ_;
      
      int whichRBin = int(r / dr);
      int whichZBin = int(z/ dz);
      
      if ((r <= len) && (z <= zLen)) 
        histo[whichRBin][whichZBin] += sd->getMass();
      
    }
    
    for(unsigned int i = 0 ; i < histo.size(); ++i){

      RealType rL = i * dr;
      RealType rU = rL + dr;
      RealType volSlice = NumericConstant::PI * dz * (( rU*rU ) - ( rL*rL ));

      for (unsigned int j = 0; j < histo[i].size(); ++j){
        histo[i][j] *= PhysicalConstants::densityConvert / volSlice;
      }
    }

    for (unsigned int i = 0; i < histo.size(); ++i) {
      RealType ther = dr * (i + 0.5);
      for(unsigned int j = 0; j < histo[i].size(); ++j) {
        if (histo[i][j] <= threshDens_) {
          RealType thez = dz * (j + 0.5);
          cerr << ther << "\t" << thez << "\n";
          break;
        }
      }
    }

    // values_.push_back( acos(maxct)*(180.0/M_PI) );
    
  }   
}


