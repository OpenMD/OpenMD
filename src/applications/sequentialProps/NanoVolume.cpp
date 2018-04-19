/* Copyright (c) 2006, 2009, 2010 The University of Notre Dame. All Rights Reserved.
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
 *
 *  NanoVolume.cpp
 *
 *  Created by Charles F. Vardeman II on 14 Dec 2006.
 *  @author  Charles F. Vardeman II
 *  @version $Id$
 *
 */

#include "applications/sequentialProps/NanoVolume.hpp"
#if defined(HAVE_QHULL)
#include "math/ConvexHull.hpp"
#include "math/AlphaHull.hpp"
#endif
#include "utils/simError.h"
#include "io/DumpReader.hpp"
#include "primitives/Molecule.hpp"

namespace OpenMD {

  NanoVolume::NanoVolume(SimInfo* info,
			 const std::string& sele)
    : NonSpatialStatistics(info, sele, 1), selectionScript_(sele),
      seleMan_(info),evaluator_(info) {
    
    string prefixFileName = info->getPrefixFileName();
    setOutputName(prefixFileName + ".avol");
    
    osq_.open(getOutputFileName().c_str());
    
    evaluator_.loadScriptString(sele);
    if (!evaluator_.isDynamic()) {
      seleMan_.setSelectionSet(evaluator_.evaluate());
    }

    theAtoms_.reserve(info_->getNGlobalAtoms());
    time_.clear();
    volume_.clear();
  
  }

  NanoVolume::~NanoVolume() {
  }
    
  
  void NanoVolume::processFrame(int istep) {
#if defined(HAVE_QHULL)
    StuntDouble* sd;
    Vector3d vec;
    int i;
    
    // Do convex hull for now - alpha has issues with perfect structures
    //AlphaHull* thishull = new AlphaHull(2.0);
    ConvexHull* thishull = new ConvexHull();
    
    // Clear pos vector between each frame.
    theAtoms_.clear();
    
    if (evaluator_.isDynamic()) {
      seleMan_.setSelectionSet(evaluator_.evaluate());
    }
    
    // outer loop is over the selected StuntDoubles:
    
    for (sd = seleMan_.beginSelected(i); sd != NULL;
	 sd = seleMan_.nextSelected(i)) {      
      theAtoms_.push_back(sd);      
    }
    
    /* variant below for single atoms, not StuntDoubles:
       for (mol = info_->beginMolecule(mi); mol != NULL; 
       mol = info_->nextMolecule(mi)) {
       for (atom = mol->beginAtom(ai); atom != NULL; 
       atom = mol->nextAtom(ai)) {
       theAtoms_.push_back(atom);
       }
       }
    */
    
    // Generate convex hull for this frame.
    thishull->computeHull(theAtoms_);

    RealType currentVolume = thishull->getVolume();
    volume_.push_back(currentVolume);
    
    RealType currentSurfaceArea = thishull->getArea();
    surfaceArea_.push_back(currentSurfaceArea);
    
    RealType currentTime = currentSnapshot_->getTime();
    time_.push_back(currentTime);

#else
    sprintf(painCave.errMsg, "NanoVolume: qhull support was not compiled in!\n");
    painCave.isFatal = 1;
    simError();  
#endif
    
  }

  void NanoVolume::processStuntDouble(StuntDouble* sd, int bin) {
    // Fill in later
  }

  void NanoVolume::writeOutput() {
    osq_.precision(7);
    if (osq_.is_open()){
      for (std::vector<int>::size_type i = 0; i != volume_.size(); i++) {
	osq_ << time_[i] << "\t" << volume_[i]
	     << "\t"  << surfaceArea_[i] << std::endl;  
      }     
    }
  }

}
  
