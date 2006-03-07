/*
 * Copyright (c) 2005 The University of Notre Dame. All Rights Reserved.
 *
 * The University of Notre Dame grants you ("Licensee") a
 * non-exclusive, royalty free, license to use, modify and
 * redistribute this software in source and binary code form, provided
 * that the following conditions are met:
 *
 * 1. Acknowledgement of the program authors must be made in any
 *    publication of scientific results based in part on use of the
 *    program.  An acceptable form of acknowledgement is citation of
 *    the article in which the program was described (Matthew
 *    A. Meineke, Charles F. Vardeman II, Teng Lin, Christopher
 *    J. Fennell and J. Daniel Gezelter, "OOPSE: An Object-Oriented
 *    Parallel Simulation Engine for Molecular Dynamics,"
 *    J. Comput. Chem. 26, pp. 252-271 (2005))
 *
 * 2. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 *
 * 3. Redistributions in binary form must reproduce the above copyright
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
 *
 *  RhoZ.cpp
 *  OOPSE-2.0
 *
 *  Created by Charles F. Vardeman II on 11/26/05.
 *  @author  Charles F. Vardeman II 
 *  @version $Id: RhoZ.cpp,v 1.4 2006-03-07 16:43:52 gezelter Exp $
 *
 */

/* Calculates Rho(Z) for density profile of liquid slab. */

#include <algorithm>
#include <fstream>
#include "applications/staticProps/RhoZ.hpp"
#include "utils/simError.h"
#include "io/DumpReader.hpp"
#include "primitives/Molecule.hpp"
namespace oopse {
  
  RhoZ::RhoZ(SimInfo* info, const std::string& filename, const std::string& sele, double len, int nrbins)
    : StaticAnalyser(info, filename), selectionScript_(sele),  evaluator_(info), seleMan_(info), len_(len), nRBins_(nrbins){

    evaluator_.loadScriptString(sele);
    if (!evaluator_.isDynamic()) {
      seleMan_.setSelectionSet(evaluator_.evaluate());
    }
      
    
    deltaR_ = len_ /(nRBins_);
    
    sliceSDLists_.resize(nRBins_);
    density_.resize(nRBins_);
    
    setOutputName(getPrefix(filename) + ".RhoZ");
  }

  void RhoZ::process() {
    DumpReader reader(info_, dumpFilename_);    
    int nFrames = reader.getNFrames();
    nProcessed_ = nFrames/step_;

    for (int istep = 0; istep < nFrames; istep += step_) {

      int i;    
      for (i=0; i < nRBins_; i++) {
        sliceSDLists_[i].clear();
      }

      StuntDouble* sd;
      reader.readFrame(istep);
      currentSnapshot_ = info_->getSnapshotManager()->getCurrentSnapshot();
      
      double sliceVolume = currentSnapshot_->getVolume() /nRBins_;
      //assume simulation box will never change
      //Mat3x3d hmat = currentSnapshot_->getHmat();
      double halfBoxZ_ = len_ / 2.0;      
        
        if (evaluator_.isDynamic()) {
          seleMan_.setSelectionSet(evaluator_.evaluate());
        }

        //wrap the stuntdoubles into a cell      
        for (sd = seleMan_.beginSelected(i); sd != NULL; sd = seleMan_.nextSelected(i)) {
            Vector3d pos = sd->getPos();
            currentSnapshot_->wrapVector(pos);
            sd->setPos(pos);
        }

        //determine which atom belongs to which slice
        for (sd = seleMan_.beginSelected(i); sd != NULL; sd = seleMan_.nextSelected(i)) {
           Vector3d pos = sd->getPos();
           //int binNo = (pos.z() /deltaR_) - 1;
           int binNo = (pos.z() + halfBoxZ_) /deltaR_   ;
           //std::cout << "pos.z = " << pos.z() << " halfBoxZ_ = " << halfBoxZ_ << " deltaR_ = "  << deltaR_ << " binNo = " << binNo << "\n";
           sliceSDLists_[binNo].push_back(sd);
        }

        //loop over the slices to calculate the densities
        for (i = 0; i < nRBins_; i++) {
            double totalMass = 0;
            for (int k = 0; k < sliceSDLists_[i].size(); ++k) {
                totalMass += sliceSDLists_[i][k]->getMass();
            }
            density_[i] += totalMass/sliceVolume;
        }
    }

    writeDensity();

  }

 
    
  void RhoZ::writeDensity() {
    std::ofstream rdfStream(outputFilename_.c_str());
    if (rdfStream.is_open()) {
      rdfStream << "#RhoZ\n";
      rdfStream << "#nFrames:\t" << nProcessed_ << "\n";
      rdfStream << "#selection: (" << selectionScript_ << ")\n";
      rdfStream << "#z\tdensity\n";
      for (int i = 0; i < density_.size(); ++i) {
        double r = deltaR_ * (i + 0.5);
        rdfStream << r << "\t" << 1.660535*density_[i]/nProcessed_ << "\n";
      }
      
    } else {
      
      sprintf(painCave.errMsg, "RhoZ: unable to open %s\n", outputFilename_.c_str());
      painCave.isFatal = 1;
      simError();  
    }
    
    rdfStream.close();
  }
  
}

