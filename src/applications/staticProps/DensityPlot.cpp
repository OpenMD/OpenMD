 /*
 * Copyright (c) 2004-2020 The University of Notre Dame. All Rights Reserved.
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
 * [1] Meineke, et al., J. Comp. Chem. 26, 252-271 (2005).
 * [2] Fennell & Gezelter, J. Chem. Phys. 124, 234104 (2006).
 * [3] Sun, Lin & Gezelter, J. Chem. Phys. 128, 234107 (2008).
 * [4] Vardeman, Stocker & Gezelter, J. Chem. Theory Comput. 7, 834 (2011).
 * [5] Kuang & Gezelter, Mol. Phys., 110, 691-701 (2012).
 * [6] Lamichhane, Gezelter & Newman, J. Chem. Phys. 141, 134109 (2014).
 * [7] Lamichhane, Newman & Gezelter, J. Chem. Phys. 141, 134110 (2014).
 * [8] Bhattarai, Newman & Gezelter, Phys. Rev. B 99, 094106 (2019).
 */

#include <algorithm>
#include <functional>
#include <memory>

#include "applications/staticProps/DensityPlot.hpp"
#include "utils/simError.h"
#include "io/DumpReader.hpp"
#include "primitives/Molecule.hpp"
#include "utils/Constants.hpp"
#include "types/LennardJonesAdapter.hpp"

using namespace std;

namespace OpenMD {

  
  DensityPlot::DensityPlot(SimInfo* info, const std::string& filename, 
			   const std::string& sele, const std::string& cmSele,
			   RealType len, int nrbins) 
    : StaticAnalyser(info, filename, nrbins), 
      len_(len), halfLen_(len/2), nRBins_(nrbins),
      selectionScript_(sele), seleMan_(info), evaluator_(info), 
      cmSelectionScript_(cmSele), cmSeleMan_(info), cmEvaluator_(info) {

    setOutputName(getPrefix(filename) + ".density");
    
    deltaR_ = len_ /nRBins_;  
    histogram_.resize(nRBins_);
    density_.resize(nRBins_);
    
    std::fill(histogram_.begin(), histogram_.end(), 0);  
    
    evaluator_.loadScriptString(sele);

    if (!evaluator_.isDynamic()) {
      seleMan_.setSelectionSet(evaluator_.evaluate());
    }

    cmEvaluator_.loadScriptString(cmSele);
    if (!cmEvaluator_.isDynamic()) {
      cmSeleMan_.setSelectionSet(cmEvaluator_.evaluate());
    }    
  }

  void DensityPlot::process() {

    bool usePeriodicBoundaryConditions_ = info_->getSimParams()->getUsePeriodicBoundaryConditions();

    DumpReader reader(info_, dumpFilename_);    
    int nFrames = reader.getNFrames();
    for (int i = 0; i < nFrames; i += step_) {
      reader.readFrame(i);
      currentSnapshot_ = info_->getSnapshotManager()->getCurrentSnapshot();
      
      if (evaluator_.isDynamic()) {
	seleMan_.setSelectionSet(evaluator_.evaluate());
      }

      if (cmEvaluator_.isDynamic()) {
	cmSeleMan_.setSelectionSet(cmEvaluator_.evaluate());
      }

      Vector3d origin = calcNewOrigin();

      Mat3x3d hmat = currentSnapshot_->getHmat();
      RealType slabVolume = deltaR_ * hmat(0, 0) * hmat(1, 1);
      int k; 
      for (StuntDouble* sd = seleMan_.beginSelected(k); sd != NULL; 
	   sd = seleMan_.nextSelected(k)) {


        if (!sd->isAtom()) {
          sprintf( painCave.errMsg, 
		   "Can not calculate electron density if it is not atom\n");
          painCave.severity = OPENMD_ERROR;
          painCave.isFatal = 1;
          simError(); 
        }
            
        Atom* atom = static_cast<Atom*>(sd);
	std::shared_ptr<GenericData> data = atom->getAtomType()->getPropertyByName("nelectron");
        if (data == nullptr) {
          sprintf( painCave.errMsg, "Can not find Parameters for nelectron\n");
          painCave.severity = OPENMD_ERROR;
          painCave.isFatal = 1;
          simError(); 
        }
            
	std::shared_ptr<DoubleGenericData> doubleData = std::dynamic_pointer_cast<DoubleGenericData>(data);
        if (doubleData == nullptr) {
          sprintf( painCave.errMsg,
                   "Can not cast GenericData to DoubleGenericData\n");
          painCave.severity = OPENMD_ERROR;
          painCave.isFatal = 1;
          simError();   
        }
            
        RealType nelectron = doubleData->getData();
        LennardJonesAdapter lja = LennardJonesAdapter(atom->getAtomType());
        RealType sigma = lja.getSigma() * 0.5;
        RealType sigma2 = sigma * sigma;
            
        Vector3d pos = sd->getPos() - origin;
        for (int j =0; j < nRBins_; ++j) {
          Vector3d tmp(pos);
          RealType zdist =j * deltaR_ - halfLen_;
          tmp[2] += zdist;
          if (usePeriodicBoundaryConditions_) 
            currentSnapshot_->wrapVector(tmp);
              
          RealType wrappedZdist = tmp.z() + halfLen_;
          if (wrappedZdist < 0.0 || wrappedZdist > len_) {
            continue;
          }
              
          int which = int(wrappedZdist / deltaR_);
          density_[which] += nelectron * exp(-zdist*zdist/(sigma2*2.0)) /(slabVolume* sqrt(2*Constants::PI*sigma*sigma));
              
        }            
      }        
    }
  
    int nProcessed = nFrames /step_;
    std::transform(density_.begin(), density_.end(), density_.begin(), 
		   std::bind(std::divides<RealType>(), placeholders::_1, nProcessed)); 
    writeDensity();
  }

  Vector3d DensityPlot::calcNewOrigin() {

    int i;
    Vector3d newOrigin(0.0);
    RealType totalMass = 0.0;
    for (StuntDouble* sd = seleMan_.beginSelected(i); sd != NULL; 
	 sd = seleMan_.nextSelected(i)) {
      RealType mass = sd->getMass();
      totalMass += mass;
      newOrigin += sd->getPos() * mass;        
    }
    newOrigin /= totalMass;
    return newOrigin;
  }

  void DensityPlot::writeDensity() {
    std::ofstream ofs(outputFilename_.c_str(), std::ios::binary);
    if (ofs.is_open()) {
      ofs << "#g(x, y, z)\n";
      ofs << "#selection: (" << selectionScript_ << ")\n";
      ofs << "#cmSelection:(" << cmSelectionScript_ << ")\n";
      ofs << "#nRBins = " << nRBins_ << "\t maxLen = " 
	  << len_ << "\tdeltaR = " << deltaR_ <<"\n";
      for (unsigned int i = 0; i < histogram_.size(); ++i) {
        ofs << i*deltaR_ - halfLen_ <<"\t" << density_[i]<< std::endl;
      }        
    } else {

      sprintf(painCave.errMsg, "DensityPlot: unable to open %s\n", 
	      outputFilename_.c_str());
      painCave.isFatal = 1;
      simError();  
    }

    ofs.close();


  }

}


