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

#include <algorithm>
#include <fstream>
#include "applications/staticProps/GofXyz.hpp"
#include "utils/simError.h"
#include "primitives/Molecule.hpp"
#include "types/MultipoleAdapter.hpp"

namespace OpenMD {

  GofXyz::GofXyz(SimInfo* info, const std::string& filename, const std::string& sele1, const std::string& sele2, const std::string& sele3, RealType len, int nrbins)
    : RadialDistrFunc(info, filename, sele1, sele2), evaluator3_(info), seleMan3_(info), len_(len), halfLen_(len/2), nRBins_(nrbins) {
      setOutputName(getPrefix(filename) + ".gxyz");

      evaluator3_.loadScriptString(sele3);
      if (!evaluator3_.isDynamic()) {
        seleMan3_.setSelectionSet(evaluator3_.evaluate());
      }    

      deltaR_ =  len_ / nRBins_;
    
      histogram_.resize(nRBins_);
      for (int i = 0 ; i < nRBins_; ++i) {
        histogram_[i].resize(nRBins_);
        for(int j = 0; j < nRBins_; ++j) {
	  histogram_[i][j].resize(nRBins_);
        }
      }   
   
    }


  void GofXyz::preProcess() {
    for (int i = 0 ; i < nRBins_; ++i) {
      histogram_[i].resize(nRBins_);
      for(int j = 0; j < nRBins_; ++j) {
	std::fill(histogram_[i][j].begin(), histogram_[i][j].end(), 0);
      }
    }   
  }


  void GofXyz::initalizeHistogram() {
    //calculate the center of mass of the molecule of selected stuntdouble in selection1

    if (!evaluator3_.isDynamic()) {
      seleMan3_.setSelectionSet(evaluator3_.evaluate());
    }    

    assert(seleMan1_.getSelectionCount() == seleMan3_.getSelectionCount());
    
    //dipole direction of selection3 and position of selection3 will be used to determine the y-z plane
    //v1 = s3 -s1, 
    //z = origin.dipole
    //x = v1 X z
    //y = z X x 
    rotMats_.clear();

    int i;
    int j;
    StuntDouble* sd1;
    StuntDouble* sd3;
    
    for (sd1 = seleMan1_.beginSelected(i), sd3 = seleMan3_.beginSelected(j); 
	 sd1 != NULL || sd3 != NULL;
	 sd1 = seleMan1_.nextSelected(i), sd3 = seleMan3_.nextSelected(j)) {

      Vector3d r3 = sd3->getPos();
      Vector3d r1 = sd1->getPos();
      Vector3d v1 =  r3 - r1;
      if (usePeriodicBoundaryConditions_)
        info_->getSnapshotManager()->getCurrentSnapshot()->wrapVector(v1);

      AtomType* atype1 = static_cast<Atom*>(sd1)->getAtomType();
      MultipoleAdapter ma1 = MultipoleAdapter(atype1);

      Vector3d zaxis;
      if (ma1.isDipole()) 
        zaxis = sd1->getDipole();
      else
        zaxis = sd1->getA().transpose() * V3Z;

      Vector3d xaxis = cross(v1, zaxis);
      Vector3d yaxis = cross(zaxis, xaxis);

      xaxis.normalize();
      yaxis.normalize();
      zaxis.normalize();

      RotMat3x3d rotMat;
      rotMat.setRow(0, xaxis);
      rotMat.setRow(1, yaxis);
      rotMat.setRow(2, zaxis);
        
      rotMats_.insert(std::map<int, RotMat3x3d>::value_type(sd1->getGlobalIndex(), rotMat));
    }

  }

  void GofXyz::collectHistogram(StuntDouble* sd1, StuntDouble* sd2) {

    Vector3d pos1 = sd1->getPos();
    Vector3d pos2 = sd2->getPos();
    Vector3d r12 = pos2 - pos1;
    if (usePeriodicBoundaryConditions_)
      currentSnapshot_->wrapVector(r12);

    std::map<int, RotMat3x3d>::iterator i = rotMats_.find(sd1->getGlobalIndex());
    assert(i != rotMats_.end());
    
    Vector3d newR12 = i->second * r12;
    // x, y and z's possible values range -halfLen_ to halfLen_
    int xbin = (newR12.x() + halfLen_) / deltaR_;
    int ybin = (newR12.y() + halfLen_) / deltaR_;
    int zbin = (newR12.z() + halfLen_) / deltaR_;

    if (xbin < nRBins_ && xbin >=0 &&
        ybin < nRBins_ && ybin >= 0 &&
        zbin < nRBins_ && zbin >=0 ) {
      ++histogram_[xbin][ybin][zbin];
    }
    
  }

  void GofXyz::writeRdf() {
    std::ofstream rdfStream(outputFilename_.c_str(), std::ios::binary);
    if (rdfStream.is_open()) {
      //rdfStream << "#g(x, y, z)\n";
      //rdfStream << "#selection1: (" << selectionScript1_ << ")\t";
      //rdfStream << "selection2: (" << selectionScript2_ << ")\n";
      //rdfStream << "#nRBins = " << nRBins_ << "\t maxLen = " << len_ << "deltaR = " << deltaR_ <<"\n";
      for (unsigned int i = 0; i < histogram_.size(); ++i) { 
	for(unsigned int j = 0; j < histogram_[i].size(); ++j) { 
	  for(unsigned int k = 0;k < histogram_[i][j].size(); ++k) {
	    rdfStream.write(reinterpret_cast<char *>(&histogram_[i][j][k] ),
                            sizeof(histogram_[i][j][k] ));
	  }
	}
      }
        
    } else {

      sprintf(painCave.errMsg, "GofXyz: unable to open %s\n", outputFilename_.c_str());
      painCave.isFatal = 1;
      simError();  
    }

    rdfStream.close();
  }

}
