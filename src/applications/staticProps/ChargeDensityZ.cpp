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
                                 const std::string& sele, int nzbins, RealType vRadius, int axis)
    : StaticAnalyser(info, filename, nzbins), selectionScript_(sele),
      evaluator_(info), seleMan_(info), thermo_(info), axis_(axis), vRadius_(vRadius) {

    evaluator_.loadScriptString(sele);
    if (!evaluator_.isDynamic()) {
      seleMan_.setSelectionSet(evaluator_.evaluate());
    }

    // fixed number of bins


    densityZAverageAllFrame_.resize(nBins_);
    densityFlucZAverageAllFrame_.resize(nBins_);
    absDensityFlucZAverageAllFrame_.resize(nBins_);
    std::fill(densityFlucZAverageAllFrame_.begin(), densityFlucZAverageAllFrame_.end(), 0);
    std::fill(densityZAverageAllFrame_.begin(), densityZAverageAllFrame_.end(), 0);
    std::fill(absDensityFlucZAverageAllFrame_.begin(), absDensityFlucZAverageAllFrame_.end(), 0);


    densityFlucZAverageFirstFrame_.resize(nBins_);
    absDensityFlucZAverageFirstFrame_.resize(nBins_);
    std::fill(densityFlucZAverageFirstFrame_.begin(), densityFlucZAverageFirstFrame_.end(), 0);
    std::fill(absDensityFlucZAverageFirstFrame_.begin(), absDensityFlucZAverageFirstFrame_.end(), 0);

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
    bool usePeriodicBoundaryConditions_ = info_->getSimParams()->getUsePeriodicBoundaryConditions();


    DumpReader reader(info_, dumpFilename_);
    int nFrames = reader.getNFrames();
    nProcessed_ = nFrames/step_;




//test using global index starts hereby



reader.readFrame(1);
currentSnapshot_ = info_->getSnapshotManager()->getCurrentSnapshot();

Mat3x3d hmat = currentSnapshot_->getHmat();
zBox_.push_back(hmat(axis_,axis_));

std::vector<RealType>::iterator j;
RealType zSum = 0.0;
for (j = zBox_.begin(); j != zBox_.end(); ++j) {
  zSum += *j;
}

RealType zAve = zSum / zBox_.size();
set<int> axisValues {0,1,2};
axisValues.erase(axis_);
set<int>:: iterator axis_it;
RealType area = 1.0;
for( axis_it = axisValues.begin(); axis_it!=axisValues.end(); ++axis_it){
  area *= hmat(*axis_it,*axis_it);
}

RealType sliceVolume = zAve/densityZAverageAllFrame_.size() * area;

for (int i = 1; i < nFrames; i += step_) {
  reader.readFrame(i);
  currentSnapshot_ = info_->getSnapshotManager()->getCurrentSnapshot();

  if (evaluator_.isDynamic()) {
     seleMan_.setSelectionSet(evaluator_.evaluate());
  }

  int kk;
  for (StuntDouble* sd = seleMan_.beginSelected(kk); sd != NULL;
  sd = seleMan_.nextSelected(kk)) {

    if (!sd->isAtom()) {
      sprintf( painCave.errMsg,
      "Can not calculate charge density distribution if it is not atom\n");
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
    int obanum(0);
    RealType sigma(0);
    std::string atomName;
    std::vector<AtomType*> atChain = atomType->allYourBase();
    std::vector<AtomType*>::iterator i;
    for (i = atChain.begin(); i != atChain.end(); ++i) {
      atomName = (*i)->getName().c_str();
      obanum = etab.GetAtomicNum((*i)->getName().c_str());
      if (obanum != 0) {
        sigma = etab.GetVdwRad(obanum);
        break;
      }
    }
    if (obanum == 0) {
      sigma = vRadius_;
    }



    Vector3d pos = sd->getPos();
    if (usePeriodicBoundaryConditions_)
      currentSnapshot_->wrapVector(pos);
    sd->setPos(pos);

    pos = sd->getPos();

    int globalIndex = sd->getGlobalIndex();



    atomNameGlobalIndex_[globalIndex] = atomName;
    averageChargeUsingGlobalIndex_[globalIndex] = (countUsingGlobalIndex_[globalIndex]*
                                                  averageChargeUsingGlobalIndex_[globalIndex] + q) / (countUsingGlobalIndex_[globalIndex] + 1);
    totalChargeUsingGlobalIndex_[globalIndex].push_back(q);
    zPosUsingGlobalIndex_[globalIndex].push_back(pos[axis_]);
    vanderRUsingGlobalIndex_[globalIndex] = sigma;
    countUsingGlobalIndex_[globalIndex]++;


  }


  for (auto const element: atomNameGlobalIndex_ ){
    int g_Index = element.first;
    std::string a_name = element.second;
    averageChargeForEachType_[a_name] = ( SDCount_[a_name] * averageChargeForEachType_[a_name] +
                                        totalChargeUsingGlobalIndex_[g_Index].front() ) / ( SDCount_[a_name] + 1 );
   SDCount_[a_name]++;


  }


  }
  RealType epsilon = 1e-10;
  for (auto const element: atomNameGlobalIndex_ ){
    int key = element.first;
    std::string a_name = element.second;

    RealType averageCharge = averageChargeUsingGlobalIndex_[key];
    RealType averageChargeUsingFirstFrame = averageChargeForEachType_[a_name];
    RealType v_radius = vanderRUsingGlobalIndex_[key];
    RealType v_radius2 = v_radius * v_radius;
    for(size_t index = 0; index < zPosUsingGlobalIndex_[key].size(); ++index){
      RealType z_pos = zPosUsingGlobalIndex_[key][index];
      RealType charge = totalChargeUsingGlobalIndex_[key][index];

      for (size_t i = 0; i < densityFlucZAverageAllFrame_.size(); ++i) {
        RealType z = -zAve/2 + i * zAve/densityZAverageAllFrame_.size();
        RealType zdist = z_pos - z;
        RealType chargeDiff = fabs(charge - averageCharge) > epsilon ? charge - averageCharge : 0;
        RealType chargeDiffUsingFirstFrameAverage = fabs(charge - averageChargeUsingFirstFrame) > epsilon ? charge - averageChargeUsingFirstFrame : 0;

        RealType gaussian_scale = exp(-zdist*zdist/(v_radius2*2.0)) / (sliceVolume *(sqrt(2*Constants::PI*v_radius*v_radius)));
        densityZAverageAllFrame_[i] += charge * gaussian_scale ;

        densityFlucZAverageAllFrame_[i] += chargeDiff * gaussian_scale;
        absDensityFlucZAverageAllFrame_[i] += abs(chargeDiff) * gaussian_scale;

        densityFlucZAverageFirstFrame_[i] += chargeDiffUsingFirstFrameAverage * gaussian_scale;
        absDensityFlucZAverageFirstFrame_[i] += abs(chargeDiffUsingFirstFrameAverage) * gaussian_scale;
    }

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
      rdfStream << "#" << axisLabel_ << "\tchargeDensity\tchargeDensityFluctuations_Average_on_all_frame\t"
                                    <<"Absolute_chargeDensityFluctuations_Average_on_all_frame\t"
                                    <<"chargeDensityFluctuations_Average_On_first_frame\t"
                                    <<"Absolute_chargeDensityFluctuations_Average_on_first_frame\n";

      for (unsigned int i = 0; i < densityFlucZAverageAllFrame_.size(); ++i) {
        RealType z = zAve * (i+0.5)/densityFlucZAverageAllFrame_.size();
        rdfStream << z << "\t"
                  <<densityZAverageAllFrame_[i]/nProcessed_<<"\t"
                  << densityFlucZAverageAllFrame_[i] / (nProcessed_)<<"\t"
                  << absDensityFlucZAverageAllFrame_[i] / (nProcessed_)<<"\t"
                  << densityFlucZAverageFirstFrame_[i] / (nProcessed_)<<"\t"
                  << absDensityFlucZAverageFirstFrame_[i] / (nProcessed_)<<"\n";
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
