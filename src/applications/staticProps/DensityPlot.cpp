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
 */

#include <algorithm>
#include "applications/staticProps/DensityPlot.hpp"
#include "utils/simError.h"
#include "io/DumpReader.hpp"
#include "primitives/Molecule.hpp"
#include "utils/NumericConstant.hpp"
namespace oopse {


DensityPlot::DensityPlot(SimInfo* info, const std::string& filename, const std::string& sele, double len, int nrbins)
  : StaticAnalyser(info, filename), selectionScript_(sele), evaluator_(info), seleMan_(info), 
    len_(len), nRBins_(nrbins), halfLen_(len/2)     {

    setOutputName(getPrefix(filename) + ".density");

    deltaR_ = len_ /nRBins_;  
    histogram_.resize(nRBins_);
    density_.resize(nRBins_);
    
    std::fill(histogram_.begin(), histogram_.end(), 0);  
    
    evaluator_.loadScriptString(sele);

    if (!evaluator_.isDynamic()) {
      seleMan_.setSelectionSet(evaluator_.evaluate());
    }
    
  }

void DensityPlot::process() {
  Molecule* mol;
  RigidBody* rb;
  SimInfo::MoleculeIterator mi;
  Molecule::RigidBodyIterator rbIter;
  
  DumpReader reader(info_, dumpFilename_);    
  int nFrames = reader.getNFrames();

  for (int i = 0; i < nFrames; i += step_) {
    reader.readFrame(i);
    currentSnapshot_ = info_->getSnapshotManager()->getCurrentSnapshot();

    for (mol = info_->beginMolecule(mi); mol != NULL; mol = info_->nextMolecule(mi)) {
        //change the positions of atoms which belong to the rigidbodies
        for (rb = mol->beginRigidBody(rbIter); rb != NULL; rb = mol->nextRigidBody(rbIter)) {
            rb->updateAtoms();
        }
        
    }
    
    Vector3d cm = info_->getCom();

    if (evaluator_.isDynamic()) {
	seleMan_.setSelectionSet(evaluator_.evaluate());
    }

    int i;        
    for (StuntDouble* sd = seleMan_.beginSelected(i); sd != NULL; sd = seleMan_.nextSelected(i)) {
        Vector3d pos = sd->getPos() - cm;
        currentSnapshot_->wrapVector(pos);

        int which = (pos.z() + halfLen_) / deltaR_;        
        if (which < nRBins_ && which >=0 ) {
          ++histogram_[which];
        }        
    }
        
  }
  
  //std::vector<int>::iterator it = std::max_element(histogram_.begin(), histogram_.end());
  //int maxnum = *it;

  int nProcessed = nFrames / step_;
  int ntot = info_->getNGlobalAtoms();
  int maxnum = ntot * nProcessed;
  std::transform(histogram_.begin(), histogram_.end(), density_.begin(), std::bind2nd(std::divides<double>(), maxnum));
  writeDensity();
  
}

void DensityPlot::writeDensity() {
    std::ofstream ofs(outputFilename_.c_str(), std::ios::binary);
    if (ofs.is_open()) {
      ofs << "#g(x, y, z)\n";
      ofs << "#selection: (" << selectionScript_ << ")\t";
      ofs << "#nRBins = " << nRBins_ << "\t maxLen = " << len_ << "deltaR = " << deltaR_ <<"\n";
      for (int i = 0; i < histogram_.size(); ++i) {
          ofs << i*deltaR_ - halfLen_ <<"\t" << density_[i]<< std::endl;
      }        
    } else {

      sprintf(painCave.errMsg, "DensityPlot: unable to open %s\n", outputFilename_.c_str());
      painCave.isFatal = 1;
      simError();  
    }

    ofs.close();


}

}


