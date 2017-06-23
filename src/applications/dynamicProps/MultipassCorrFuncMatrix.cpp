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

#include "applications/dynamicProps/MultipassCorrFuncMatrix.hpp"
#include "utils/simError.h"
#include "utils/Revision.hpp"
#include "primitives/Molecule.hpp"

using namespace std;
namespace OpenMD {

  MultipassCorrFuncMatrix::MultipassCorrFuncMatrix(SimInfo* info, const string& filename,
                                       const string& sele1, const string& sele2,
                                       int storageLayout)
    : storageLayout_(storageLayout), info_(info), dumpFilename_(filename),
      seleMan1_(info_), seleMan2_(info_),
      selectionScript1_(sele1), selectionScript2_(sele2),
      evaluator1_(info_), evaluator2_(info_), autoCorrFunc_(false) {

    // Request maximum needed storage for the simulation (including of
    // whatever was passed down by the individual correlation
    // function).

    storageLayout_ = info->getStorageLayout() | storageLayout;

    reader_ = new DumpReader(info_, dumpFilename_);

    uniqueSelections_ = (sele1.compare(sele2) != 0) ? true : false;

    evaluator1_.loadScriptString(selectionScript1_);
    //if selection is static, we only need to evaluate it once
    if (!evaluator1_.isDynamic()) {
      seleMan1_.setSelectionSet(evaluator1_.evaluate());
      validateSelection(seleMan1_);
    }

    if (uniqueSelections_) {
      evaluator2_.loadScriptString(selectionScript2_);
      if (!evaluator2_.isDynamic()) {
        seleMan2_.setSelectionSet(evaluator2_.evaluate());
        validateSelection(seleMan2_);
      }
    }

    Globals* simParams = info_->getSimParams();
    if (simParams->haveSampleTime()){
      deltaTime_ = simParams->getSampleTime();
    } else {
      sprintf(painCave.errMsg,
              "MultipassCorrFuncMatrix Error: can not figure out deltaTime\n");
      painCave.isFatal = 1;
      simError();
    }

    sprintf(painCave.errMsg, "Scanning for frames.");
    painCave.isFatal = 0;
    painCave.severity=OPENMD_INFO;
    simError();

    nFrames_ = reader_->getNFrames();
    nTimeBins_ = nFrames_;
    histogram_.resize(nTimeBins_, 0.0);
    count_.resize(nTimeBins_, 0);

    times_.resize(nFrames_);
    sele1ToIndex_.resize(nFrames_);
    if (uniqueSelections_) {
      sele2ToIndex_.resize(nFrames_);
    }
  }

  void MultipassCorrFuncMatrix::preCorrelate() {

    for (int istep = 0; istep < nFrames_; istep++) {
      reader_->readFrame(istep);
      currentSnapshot_ = info_->getSnapshotManager()->getCurrentSnapshot();
      times_[istep] = currentSnapshot_->getTime();
      computeFrame(istep);
    }

  }

  void MultipassCorrFuncMatrix::computeFrame(int istep) {
    StuntDouble* sd;

    int isd1, isd2;
    unsigned int index;

    if (evaluator1_.isDynamic()) {
      seleMan1_.setSelectionSet(evaluator1_.evaluate());
    }

    if (uniqueSelections_ && evaluator2_.isDynamic()) {
      seleMan2_.setSelectionSet(evaluator2_.evaluate());
    }

    for (sd = seleMan1_.beginSelected(isd1); sd != NULL;
         sd = seleMan1_.nextSelected(isd1)) {

      index = computeProperty1(istep, sd);

      if (index == sele1ToIndex_[istep].size()) {
        sele1ToIndex_[istep].push_back(sd->getGlobalIndex());
      } else {
        sele1ToIndex_[istep].resize(index+1);
        sele1ToIndex_[istep][index] = sd->getGlobalIndex();
      }

      if (!uniqueSelections_) {
        index = computeProperty2(istep, sd);
      }
    }

    if (uniqueSelections_) {
      for (sd = seleMan2_.beginSelected(isd2); sd != NULL;
           sd = seleMan2_.nextSelected(isd2)) {

        if (autoCorrFunc_) {
          index = computeProperty1(istep, sd);
        } else {
          index = computeProperty2(istep, sd);
        }
        if (index == sele2ToIndex_[istep].size()) {
          sele2ToIndex_[istep].push_back(sd->getGlobalIndex());
        } else {
          sele2ToIndex_[istep].resize(index+1);
          sele2ToIndex_[istep][index] = sd->getGlobalIndex();
        }
      }
    }
  }


  void MultipassCorrFuncMatrix::doCorrelate() {

    painCave.isFatal = 0;
    painCave.severity=OPENMD_INFO;
    sprintf(painCave.errMsg, "Starting pre-correlate scan.");
    simError();
    preCorrelate();

    sprintf(painCave.errMsg, "Calculating correlation function.");
    simError();
    correlation();

    sprintf(painCave.errMsg, "Doing post-correlation calculations.");
    simError();
    postCorrelate();

    sprintf(painCave.errMsg, "Writing output.");
    simError();
    writeCorrelate();
  }

  void MultipassCorrFuncMatrix::correlation() {

    for (int i =0 ; i < nTimeBins_; ++i) {
      Mat3x3d Mat3Zero(0.0);//check this. To overwrite histogram with zeros
      histogram_[i] = Mat3Zero;
      //histogram_[i] = 0.0;
      count_[i] = 0;
    }

    for (int i = 0; i < nFrames_; ++i) {

      RealType time1 = times_[i];

      for(int j  = i; j < nFrames_; ++j) {

        // Perform a sanity check on the actual configuration times to
        // make sure the configurations are spaced the same amount the
        // sample time said they were spaced:

        RealType time2 = times_[j];

        if ( fabs( (time2 - time1) - (j-i)*deltaTime_ ) > 1.0e-4 ) {
          sprintf(painCave.errMsg,
                  "MultipassCorrFuncMatrix::correlateBlocks Error: sampleTime (%f)\n"
                  "\tin %s does not match actual time-spacing between\n"
                  "\tconfigurations %d (t = %f) and %d (t = %f).\n",
                  deltaTime_, dumpFilename_.c_str(), i, time1, j, time2);
          painCave.isFatal = 1;
          simError();
        }

        int timeBin = int ((time2 - time1) / deltaTime_ + 0.5);
        correlateFrames(i,j, timeBin);
      }
    }
  }

  void MultipassCorrFuncMatrix::correlateFrames(int frame1, int frame2, int timeBin) {
    std::vector<int> s1;
    std::vector<int> s2;

    std::vector<int>::iterator i1;
    std::vector<int>::iterator i2;

    Mat3x3d corrValMatrix(0.0);//check this

    s1 = sele1ToIndex_[frame1];

    if (uniqueSelections_)
       s2 = sele2ToIndex_[frame2];
    else
       s2 = sele1ToIndex_[frame2];

    for (i1 = s1.begin(), i2 = s2.begin();
         i1 != s1.end() && i2 != s2.end(); ++i1, ++i2){

      // If the selections are dynamic, they might not have the
      // same objects in both frames, so we need to roll either of
      // the selections until we have the same object to
      // correlate.

      while ( i1 != s1.end() && *i1 < *i2 ) {
        ++i1;
      }

      while ( i2 != s2.end() && *i2 < *i1 ) {
        ++i2;
      }

      if ( i1 == s1.end() || i2 == s2.end() ) break;

      corrValMatrix = calcCorrVal(frame1, frame2, i1 - s1.begin(), i2 - s2.begin());
      histogram_[timeBin].add(corrValMatrix);//check this. changes values of this to this plus corrValMatrix
      count_[timeBin]++;

    }
  }

  void MultipassCorrFuncMatrix::postCorrelate() {
    sumForces_.div(nFrames_); //gets the average of the forces_
    sumTorques_.div(nFrames_); //gets the average of the torques_
    Mat3x3d correlationOfAverages_ = outProduct(sumForces_, sumTorques_);
    for (int i =0 ; i < nTimeBins_; ++i) {
      if (count_[i] > 0) {
	       histogram_[i].div(count_[i]);//divides matrix by count_[i]
         histogram_[i].sub(correlationOfAverages_);//the outerProduct correlation of the averages is subtracted from the correlation value
      } else {
        Mat3x3d Mat3Zero(0.0);//check this. To overwrite histogram with zeros
        histogram_[i] = Mat3Zero;
      }
    }



  }

  void MultipassCorrFuncMatrix::writeCorrelate() {
    ofstream ofs(outputFilename_.c_str());

    if (ofs.is_open()) {

      Revision r;

      ofs << "# " << getCorrFuncType() << "\n";
      ofs << "# OpenMD " << r.getFullRevision() << "\n";
      ofs << "# " << r.getBuildDate() << "\n";
      ofs << "# selection script1: \"" << selectionScript1_ ;
      ofs << "\"\tselection script2: \"" << selectionScript2_ << "\"\n";
      if (!paramString_.empty())
        ofs << "# parameters: " << paramString_ << "\n";
      ofs << "#time\tcorrVal\n";

      for (int i = 0; i < nTimeBins_; ++i) {
	       ofs << times_[i]-times_[0] << "\t";
         for (int j = 0; j < 3; j++) {
           for (int k = 0; k < 3; k++) {
             ofs << histogram_[i](j,k) << '\t';
           }
         }
         ofs << '\n';
      }

    } else {
      sprintf(painCave.errMsg,
	      "MultipassCorrFuncMatrix::writeCorrelate Error: fail to open %s\n",
              outputFilename_.c_str());
      painCave.isFatal = 1;
      simError();
    }

    ofs.close();
  }

}
