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
#include "applications/staticProps/ChargeDensityZ.hpp"
#include "types/FixedChargeAdapter.hpp"
#include "types/FluctuatingChargeAdapter.hpp"
#include "utils/simError.h"
#include "io/DumpReader.hpp"
#include "primitives/Molecule.hpp"
#include "brains/Thermo.hpp"
#include "primitives/Molecule.hpp"
#include "utils/Constants.hpp"

namespace OpenMD {


  ChargeDensityZ::ChargeDensityZ(SimInfo* info, const std::string& filename,
                                 const std::string& sele, int nzbins, int axis)
    : StaticAnalyser(info, filename, nzbins), selectionScript_(sele),
      evaluator_(info), seleMan_(info), thermo_(info), axis_(axis) {

    evaluator_.loadScriptString(sele);
    if (!evaluator_.isDynamic()) {
      seleMan_.setSelectionSet(evaluator_.evaluate());
    }

    // fixed number of bins



    densityZ_.resize(nBins_);
    absDensityZ_.resize(nBins_);
    std::fill(densityZ_.begin(), densityZ_.end(), 0);
    std::fill(absDensityZ_.begin(), absDensityZ_.end(), 0);

    switch(axis_) {
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

    setOutputName(getPrefix(filename) + ".ChargeDensityZ");
  }
  void ChargeDensityZ::process() {
    StuntDouble* sd;
    int ii;
    averageCharge_ = 0;
    int countSD(0);

    bool usePeriodicBoundaryConditions_ = info_->getSimParams()->getUsePeriodicBoundaryConditions();

    DumpReader reader(info_, dumpFilename_);
    int nFrames = reader.getNFrames();
    nProcessed_ = nFrames/step_;

    bool findAverage = true;
    for (int i = 1; i < nFrames; i += step_) {

      reader.readFrame(i);
      currentSnapshot_ = info_->getSnapshotManager()->getCurrentSnapshot();

      if (evaluator_.isDynamic()) {
	       seleMan_.setSelectionSet(evaluator_.evaluate());
      }

      for (sd = seleMan_.beginSelected(ii); sd != NULL; sd = seleMan_.nextSelected(ii)) {
        Vector3d pos = sd->getPos();
        if (usePeriodicBoundaryConditions_)
          currentSnapshot_->wrapVector(pos);
        sd->setPos(pos);
      }

      Mat3x3d hmat = currentSnapshot_->getHmat();
      zBox_.push_back(hmat(axis_,axis_));

      RealType halfBoxZ_ = hmat(axis_,axis_) / 2.0;

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
        RealType q = 0.0;
        Atom* atom = static_cast<Atom*>(sd);

        AtomType* atomType = atom->getAtomType();

        if (sd->isAtom()) {
          FixedChargeAdapter fca = FixedChargeAdapter(atomType);
          if ( fca.isFixedCharge() ) {
            q += fca.getCharge();
          }

          FluctuatingChargeAdapter fqa = FluctuatingChargeAdapter(atomType);
          if ( fqa.isFluctuatingCharge() ) {
            q += atom->getFlucQPos();
          }
        }

        if(findAverage == true){
          averageCharge_ += q;
          countSD++;
        }
        else{
        int obanum(0);
        RealType sigma(0);
        std::vector<AtomType*> atChain = atomType->allYourBase();
        std::vector<AtomType*>::iterator i;
        for (i = atChain.begin(); i != atChain.end(); ++i) {
          obanum = etab.GetAtomicNum((*i)->getName().c_str());
          if (obanum != 0) {
            sigma = etab.GetVdwRad(obanum);
            break;
          }
        }
        if (obanum == 0) {
          sprintf( painCave.errMsg,
                   "Could not find atom type in default element.txt\n");
          painCave.severity = OPENMD_ERROR;
          painCave.isFatal = 1;
          simError();
        }



        RealType sigma2 = sigma * sigma;
        Vector3d pos = sd->getPos();

        // shift molecules by half a box to have bins start at 0
        int binNo = int(nBins_ * (halfBoxZ_ + pos[axis_]) / hmat(axis_,axis_));
        RealType zdist = pos[axis_] - binNo * (hmat(axis_,axis_)/nBins_) ;
        densityZ_[binNo] += (q - (averageCharge_/countSD)) * exp(-zdist*zdist/(sigma2*2.0)) / sqrt(2*Constants::PI*sigma*sigma);
        absDensityZ_[binNo] += abs(q - (averageCharge_/countSD)) * exp(-zdist*zdist/(sigma2*2.0)) / sqrt(2*Constants::PI*sigma*sigma);
        }
      }

        if(findAverage == true){
          i--;
          findAverage = false;

      }
    }

    writeDensity();



  }



  void ChargeDensityZ::writeDensity() {
    // compute average box length:
    std::vector<RealType>::iterator j;
    RealType zSum = 0.0;
    for (j = zBox_.begin(); j != zBox_.end(); ++j) {
      zSum += *j;
    }
    RealType zAve = zSum / zBox_.size();

    std::ofstream rdfStream(outputFilename_.c_str());
    if (rdfStream.is_open()) {
      rdfStream << "#ChargeDensityZ "<<"\n";
      rdfStream << "#selection: (" << selectionScript_ << ")\n";
      rdfStream << "#" << axisLabel_ << "\tchargeDensity\tAbsolute chargeDensity\n";
      for (unsigned int i = 0; i < densityZ_.size(); ++i) {
        RealType z = zAve * (i+0.5)/densityZ_.size();
        rdfStream << z << "\t"
                  << densityZ_[i] / (nProcessed_)<<"\t"
                  << absDensityZ_[i] / (nProcessed_)<<"\n";
      }

    } else {

      sprintf(painCave.errMsg, "ChargeDensityZ: unable to open %s\n",
	      outputFilename_.c_str());
      painCave.isFatal = 1;
      simError();
    }

    rdfStream.close();
  }

}
