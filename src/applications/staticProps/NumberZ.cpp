/*
 * Copyright (c) 2004-2021 The University of Notre Dame. All Rights Reserved.
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

/*
 * Computes the number density distribution along preferred axis for the
 * selected molecules
 */

#include "applications/staticProps/NumberZ.hpp"

#include <algorithm>
#include <fstream>

#include "brains/Thermo.hpp"
#include "io/DumpReader.hpp"
#include "primitives/Molecule.hpp"
#include "utils/simError.h"

namespace OpenMD {

  NumberZ::NumberZ(SimInfo* info, const std::string& filename,
                   const std::string& sele, int nzbins, int axis) :
      StaticAnalyser(info, filename, nzbins),
      selectionScript_(sele), evaluator_(info), seleMan_(info), thermo_(info),
      axis_(axis) {
    
    evaluator_.loadScriptString(sele);
    if (!evaluator_.isDynamic()) {
      seleMan_.setSelectionSet(evaluator_.evaluate());
    }

    numberZ_.resize(nBins_);

    switch (axis_) {
    case 0:
      axisLabel_ = "x";
      break;
    case 1:
      axisLabel_ = "y";
      break;
    case 2:
    default:
      axisLabel_ = "z";
      break;
    }

    setOutputName(getPrefix(filename) + ".NumberZ");
  }

  void NumberZ::process() {
    Molecule* mol;
    int ii;

    bool usePeriodicBoundaryConditions_ =
        info_->getSimParams()->getUsePeriodicBoundaryConditions();

    DumpReader reader(info_, dumpFilename_);
    int nFrames = reader.getNFrames();
    nProcessed_ = nFrames / step_;

    for (int istep = 0; istep < nFrames; istep += step_) {
      reader.readFrame(istep);
      currentSnapshot_ = info_->getSnapshotManager()->getCurrentSnapshot();

      Mat3x3d hmat = currentSnapshot_->getHmat();
      zBox_.push_back(hmat(axis_, axis_));

      RealType halfBoxZ_ = hmat(axis_, axis_) / 2.0;
      RealType area = 0.0;
      switch (axis_) {
      case 0:
        area = currentSnapshot_->getYZarea();
        break;
      case 1:
        area = currentSnapshot_->getXZarea();
        break;
      case 2:
      default:
        area = currentSnapshot_->getXYarea();
        break;
      }

      areas_.push_back(area);

      if (evaluator_.isDynamic()) {
        seleMan_.setSelectionSet(evaluator_.evaluate());
      }

      for (mol = seleMan_.beginSelectedMolecule(ii); mol != NULL;
           mol = seleMan_.nextSelectedMolecule(ii)) {
        Vector3d pos = mol->getCom();
        if (usePeriodicBoundaryConditions_) currentSnapshot_->wrapVector(pos);
        // shift molecules by half a box to have bins start at 0
        int binNo = int(nBins_ * (halfBoxZ_ + pos[axis_]) / hmat(axis_, axis_));
        numberZ_[binNo]++;
      }
    }
    writeNumberZ();
  }

  void NumberZ::writeNumberZ() {
    // compute average box length:
    std::vector<RealType>::iterator j;
    RealType zSum = 0.0;
    for (j = zBox_.begin(); j != zBox_.end(); ++j) {
      zSum += *j;
    }
    RealType zAve = zSum / zBox_.size();
    
    RealType areaSum = 0.0;
    for (j = areas_.begin(); j != areas_.end(); ++j) {
      areaSum += *j;
    }
    RealType areaAve = areaSum / areas_.size();

    std::ofstream rdfStream(outputFilename_.c_str());
    if (rdfStream.is_open()) {
      rdfStream << "#NumberZ "
                << "\n";
      rdfStream << "#selection: (" << selectionScript_ << ")\n";
      rdfStream << "#" << axisLabel_ << "\tnumber\n";
      RealType binNumber;
      for (unsigned int i = 0; i < numberZ_.size(); ++i) {
        RealType z = zAve * (i + 0.5) / nBins_;

        RealType volSlice = areaAve * zAve / nBins_; 

        binNumber = numberZ_[i] / (volSlice * nProcessed_);

        rdfStream << z << "\t" << binNumber << "\n";
      }

    } else {
      sprintf(painCave.errMsg, "NumberZ: unable to open %s\n",
              outputFilename_.c_str());
      painCave.isFatal = 1;
      simError();
    }

    rdfStream.close();
  }
}  // namespace OpenMD
