/*
 * Copyright (c) 2004-present, The University of Notre Dame. All rights
 * reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 * 1. Redistributions of source code must retain the above copyright notice,
 *    this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright notice,
 *    this list of conditions and the following disclaimer in the documentation
 *    and/or other materials provided with the distribution.
 *
 * 3. Neither the name of the copyright holder nor the names of its
 *    contributors may be used to endorse or promote products derived from
 *    this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 *
 * SUPPORT OPEN SCIENCE!  If you use OpenMD or its source code in your
 * research, please cite the following paper when you publish your work:
 *
 * [1] Drisko et al., J. Open Source Softw. 9, 7004 (2024).
 *
 * Good starting points for code and simulation methodology are:
 *
 * [2] Meineke, et al., J. Comp. Chem. 26, 252-271 (2005).
 * [3] Fennell & Gezelter, J. Chem. Phys. 124, 234104 (2006).
 * [4] Sun, Lin & Gezelter, J. Chem. Phys. 128, 234107 (2008).
 * [5] Vardeman, Stocker & Gezelter, J. Chem. Theory Comput. 7, 834 (2011).
 * [6] Kuang & Gezelter, Mol. Phys., 110, 691-701 (2012).
 * [7] Lamichhane, Gezelter & Newman, J. Chem. Phys. 141, 134109 (2014).
 * [8] Bhattarai, Newman & Gezelter, Phys. Rev. B 99, 094106 (2019).
 * [9] Drisko & Gezelter, J. Chem. Theory Comput. 20, 4986-4997 (2024).
 */

#include "applications/dynamicProps/TimeCorrFunc.hpp"

#include <memory>

#include "math/DynamicVector.hpp"
#include "primitives/Molecule.hpp"
#include "utils/Revision.hpp"
#include "utils/simError.h"

using namespace std;
namespace OpenMD {

  template<typename T>
  TimeCorrFunc<T>::TimeCorrFunc(SimInfo* info, const string& filename,
                                const string& sele1, const string& sele2) :
      info_(info),
      dumpFilename_(filename), seleMan1_(info_), seleMan2_(info_),
      selectionScript1_(sele1), selectionScript2_(sele2), evaluator1_(info_),
      evaluator2_(info_), autoCorrFunc_(false), doSystemProperties_(false),
      doMolecularProperties_(false), doObjectProperties_(false),
      doBondProperties_(false), allowTimeFuzz_(false) {
    reader_ = new DumpReader(info_, dumpFilename_);

    uniqueSelections_ = (sele1.compare(sele2) != 0) ? true : false;

    Globals* simParams = info_->getSimParams();
    if (simParams->haveSampleTime()) {
      deltaTime_ = simParams->getSampleTime();
    } else {
      snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
               "TimeCorrFunc Error: can not figure out deltaTime\n");
      painCave.isFatal = 1;
      simError();
    }

    snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH, "Scanning for frames.");
    painCave.isFatal  = 0;
    painCave.severity = OPENMD_INFO;
    simError();

    nFrames_   = reader_->getNFrames();
    nTimeBins_ = nFrames_;

    T zeroType(0.0);
    histogram_.resize(nTimeBins_, zeroType);
    count_.resize(nTimeBins_, 0);

    times_.resize(nFrames_);
    sele1ToIndex_.resize(nFrames_);
    GIDtoSele1_.resize(nFrames_);
    selection1StartFrame_.resize(nFrames_);

    if (uniqueSelections_) { 
      sele2ToIndex_.resize(nFrames_); 
      GIDtoSele2_.resize(nFrames_);
      selection2StartFrame_.resize(nFrames_);
    }

    progressBar_ = std::make_unique<ProgressBar>();
  }

  template<typename T>
  void TimeCorrFunc<T>::preCorrelate() {
    evaluator1_.loadScriptString(selectionScript1_);
    // if selection is static, we only need to evaluate it once
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

    progressBar_->clear();

    for (int istep = 0; istep < nFrames_; istep++) {
      reader_->readFrame(istep);
      currentSnapshot_ = info_->getSnapshotManager()->getCurrentSnapshot();
      times_[istep]    = currentSnapshot_->getTime();

      progressBar_->setStatus(istep + 1, nFrames_);
      progressBar_->update();

      computeFrame(istep);
    }

    RealType dt2Avg(0.0);
    for (int istep = 1; istep < nFrames_; istep++) {
      RealType dt = times_[istep] - times_[istep - 1];
      dtMean_ += dt;
      dt2Avg += dt * dt;
    }
    dtMean_ /= RealType(nFrames_ - 1);
    dt2Avg /= RealType(nFrames_ - 1);
    dtSigma_ = std::sqrt(dt2Avg - dtMean_ * dtMean_);

    if (dtSigma_ > 1.0e-6) {
      snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
               "TimeCorrFunc::preCorrelate: sampleTime (%f) does not match\n"
               "\tthe mean spacing between configurations (%f), with\n"
               "\tsigma (%f). Proceeding with the mean value.\n",
               deltaTime_, dtMean_, dtSigma_);
      allowTimeFuzz_    = true;
      painCave.isFatal  = 0;
      painCave.severity = OPENMD_INFO;
      simError();
    }
  }

  template<typename T>
  void TimeCorrFunc<T>::computeFrame(int istep) {
    Molecule* mol;
    StuntDouble* sd;
    Bond* bond;
    int imol1, imol2, isd1, isd2, ibond1, ibond2;
    unsigned int index;

    if (evaluator1_.isDynamic()) {
      seleMan1_.setSelectionSet(evaluator1_.evaluate());
      validateSelection(seleMan1_);
    }

    if (uniqueSelections_ && evaluator2_.isDynamic()) {
      seleMan2_.setSelectionSet(evaluator2_.evaluate());
      validateSelection(seleMan2_);
    }

    if (doSystemProperties_) {
      computeProperty1(istep);
      if (!autoCorrFunc_) computeProperty2(istep);
    }

    if (doMolecularProperties_) {
      // Map of molecular global IDs to selections
      if (selectionModeRestart_) 
	GIDtoSele1_[istep].resize(info_->getNGlobalMolecules(), -1);

      for (mol = seleMan1_.beginSelectedMolecule(imol1); mol != NULL;
           mol = seleMan1_.nextSelectedMolecule(imol1)) {
        index = computeProperty1(istep, mol);

        if (index == sele1ToIndex_[istep].size()) {
          sele1ToIndex_[istep].push_back(mol->getGlobalIndex());
        } else {
          sele1ToIndex_[istep].resize(index + 1);
          sele1ToIndex_[istep][index] = mol->getGlobalIndex();
        }
	if (selectionModeRestart_) {
	  GIDtoSele1_[istep][mol->getGlobalIndex()] = index;
	  
	  if (istep == 0) {
	    selection1StartFrame_[istep].push_back(istep);
	  } else {
	    int prevIndex = GIDtoSele1_[istep-1][mol->getGlobalIndex()];
	    if (prevIndex == -1)
	      selection1StartFrame_[istep].push_back(istep);
	    else
	      selection1StartFrame_[istep].push_back(selection1StartFrame_[istep-1][prevIndex]);
	  }
	}

        if (!uniqueSelections_) { index = computeProperty2(istep, mol); }
      }

      if (uniqueSelections_) {
	// Map of molecular global IDs to selections
	if (selectionModeRestart_) 	  
	  GIDtoSele2_[istep].resize(info_->getNGlobalMolecules(), -1);

        for (mol = seleMan2_.beginSelectedMolecule(imol2); mol != NULL;
             mol = seleMan2_.nextSelectedMolecule(imol2)) {
          if (autoCorrFunc_) {
            index = computeProperty1(istep, mol);
          } else {
            index = computeProperty2(istep, mol);
          }
          if (index == sele2ToIndex_[istep].size()) {
            sele2ToIndex_[istep].push_back(mol->getGlobalIndex());
          } else {
            sele2ToIndex_[istep].resize(index + 1);
            sele2ToIndex_[istep][index] = mol->getGlobalIndex();
          }
	  if (selectionModeRestart_) {
	    
	    GIDtoSele2_[istep][mol->getGlobalIndex()] = index;
	    if (istep == 0) {
	      selection2StartFrame_[istep].push_back(istep);
	    } else {
	      int prevIndex = GIDtoSele2_[istep-1][mol->getGlobalIndex()];
	      if (prevIndex == -1)
		selection2StartFrame_[istep].push_back(istep);
	      else 	      
		selection2StartFrame_[istep].push_back(selection2StartFrame_[istep-1][prevIndex]);
	    }
	  }	  
        }
      }
    }

    if (doObjectProperties_) {
      // Map of StuntDouble global IDs to selections
      if (selectionModeRestart_)	
	GIDtoSele1_[istep].resize(info_->getNGlobalIntegrableObjects(), -1);
      
      for (sd = seleMan1_.beginSelected(isd1); sd != NULL;
           sd = seleMan1_.nextSelected(isd1)) {
        index = computeProperty1(istep, sd);

        if (index == sele1ToIndex_[istep].size()) {
          sele1ToIndex_[istep].push_back(sd->getGlobalIndex());
        } else {
          sele1ToIndex_[istep].resize(index + 1);
          sele1ToIndex_[istep][index] = sd->getGlobalIndex();
        }

	if (selectionModeRestart_) {

	  GIDtoSele1_[istep][sd->getGlobalIndex()] = index;
	  if (istep == 0) {
	    selection1StartFrame_[istep].push_back(istep);
	  } else {
	    int prevIndex = GIDtoSele1_[istep-1][sd->getGlobalIndex()];
	    if (prevIndex == -1)
	      selection1StartFrame_[istep].push_back(istep);
	    else
	      selection1StartFrame_[istep].push_back(selection1StartFrame_[istep-1][prevIndex]);
	  }
	}

        if (!uniqueSelections_) { index = computeProperty2(istep, sd); }
      }

      if (uniqueSelections_) {
	// Map of StuntDouble global IDs to selections
	if (selectionModeRestart_)
	  GIDtoSele2_[istep].resize(info_->getNGlobalIntegrableObjects(), -1);
	
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
            sele2ToIndex_[istep].resize(index + 1);
            sele2ToIndex_[istep][index] = sd->getGlobalIndex();
          }
	  if (selectionModeRestart_) {

	    GIDtoSele2_[istep][sd->getGlobalIndex()] = index;
	    if (istep == 0) {
	      selection2StartFrame_[istep].push_back(istep);
	    } else {
	      int prevIndex = GIDtoSele2_[istep-1][sd->getGlobalIndex()];
	      if (prevIndex == -1)
		selection2StartFrame_[istep].push_back(istep);
	      else
		selection2StartFrame_[istep].push_back(selection2StartFrame_[istep-1][prevIndex]);
	    }
	  }
        }
      }
    }

    if (doBondProperties_) {
      // Map of bond global IDs to selections
      if (selectionModeRestart_) 
	GIDtoSele1_[istep].resize(info_->getNGlobalBonds(), -1);

      for (bond = seleMan1_.beginSelectedBond(ibond1); bond != NULL;
           bond = seleMan1_.nextSelectedBond(ibond1)) {
        index = computeProperty1(istep, bond);

        if (index == sele1ToIndex_[istep].size()) {
          sele1ToIndex_[istep].push_back(bond->getGlobalIndex());
        } else {
          sele1ToIndex_[istep].resize(index + 1);
          sele1ToIndex_[istep][index] = bond->getGlobalIndex();
        }

	if (selectionModeRestart_) {
	  
	  GIDtoSele1_[istep][bond->getGlobalIndex()] = index;
	  if (istep == 0) {
	    selection1StartFrame_[istep].push_back(istep);
	  } else {
	    int prevIndex = GIDtoSele1_[istep-1][bond->getGlobalIndex()];
	    if (prevIndex == -1)
	      selection1StartFrame_[istep].push_back(istep);
	    else
	      selection1StartFrame_[istep].push_back(selection1StartFrame_[istep-1][prevIndex]);
	  }
	}
	
        if (!uniqueSelections_) { index = computeProperty2(istep, bond); }
      }

      if (uniqueSelections_) {
	// Map of bond global IDs to selections
	if (selectionModeRestart_) 	 
	  GIDtoSele2_[istep].resize(info_->getNGlobalBonds(), -1);

        for (bond = seleMan2_.beginSelectedBond(ibond2); bond != NULL;
             bond = seleMan2_.nextSelectedBond(ibond2)) {
          if (autoCorrFunc_) {
            index = computeProperty1(istep, bond);
          } else {
            index = computeProperty2(istep, bond);
          }
          if (index == sele2ToIndex_[istep].size()) {
            sele2ToIndex_[istep].push_back(bond->getGlobalIndex());
          } else {
            sele2ToIndex_[istep].resize(index + 1);
            sele2ToIndex_[istep][index] = bond->getGlobalIndex();
          }
	  if (selectionModeRestart_) {

	    GIDtoSele2_[istep][bond->getGlobalIndex()] = index;
	    if (istep == 0) {
	      selection2StartFrame_[istep].push_back(istep);
	    } else {
	      int prevIndex = GIDtoSele2_[istep-1][bond->getGlobalIndex()];
	      if (prevIndex == -1)
		selection2StartFrame_[istep].push_back(istep);
	      else	      
		selection2StartFrame_[istep].push_back(selection2StartFrame_[istep-1][prevIndex]);
	    }
	  }
        }
      }
    }
  }

  template<typename T>
  void TimeCorrFunc<T>::doCorrelate() {
    painCave.isFatal  = 0;
    painCave.severity = OPENMD_INFO;
    snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
             "Starting pre-correlate scan.");
    simError();
    preCorrelate();

    snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
             "Calculating correlation function.");
    simError();
    correlation();

    snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
             "Doing post-correlation calculations.");
    simError();
    postCorrelate();

    snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH, "Writing output.");
    simError();
    writeCorrelate();
  }

  template<typename T>
  void TimeCorrFunc<T>::correlation() {
    T zeroType(0.0);
    for (unsigned int i = 0; i < nTimeBins_; ++i) {
      histogram_[i] = zeroType;
      count_[i]     = 0;
    }

    progressBar_->clear();
    RealType samples = 0.5 * (nFrames_ + 1) * nFrames_;
    int visited      = 0;

    for (int i = 0; i < nFrames_; ++i) {
      RealType time1 = times_[i];

      for (int j = i; j < nFrames_; ++j) {
        visited++;
        progressBar_->setStatus(visited, samples);
        progressBar_->update();

        // Perform a sanity check on the actual configuration times to
        // make sure the configurations are spaced the same amount the
        // sample time said they were spaced:

        RealType time2 = times_[j];

        if (std::fabs((time2 - time1) - (j - i) * dtMean_) >
            6 * dtSigma_ * (j - i)) {
          snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
                   "TimeCorrFunc::correlateBlocks: mean sampleTime (%f)\n"
                   "\tin %s does not match actual time-spacing between\n"
                   "\tconfigurations %d (t = %f) and %d (t = %f).\n",
                   dtMean_, dumpFilename_.c_str(), i, time1, j, time2);
          if (allowTimeFuzz_) {
            painCave.isFatal  = 0;
            painCave.severity = OPENMD_INFO;
          } else {
            painCave.isFatal  = 1;
            painCave.severity = OPENMD_ERROR;
          }
          simError();
        }

        int timeBin = int((time2 - time1) / dtMean_ + 0.5);
        correlateFrames(i, j, timeBin);
      }
    }
  }

  /*
template<typename T>
void TimeCorrFunc<T>::validateSelection(SelectionManager& seleMan) {
}
*/

  template<typename T>
  void TimeCorrFunc<T>::correlateFrames(int frame1, int frame2, int timeBin) {
    std::vector<int> s1;
    std::vector<int> s2;

    std::vector<int>::iterator i1;
    std::vector<int>::iterator i2;

    T corrVal(0.0);
    T zeroType(0.0);

    if (doSystemProperties_) {
      corrVal = calcCorrVal(frame1, frame2);
      histogram_[timeBin] += corrVal;
      count_[timeBin]++;

    } else {
      s1 = sele1ToIndex_[frame1];

      if (uniqueSelections_)
        s2 = sele2ToIndex_[frame2];
      else
        s2 = sele1ToIndex_[frame2];

      for (i1 = s1.begin(), i2 = s2.begin(); i1 != s1.end() && i2 != s2.end();
           ++i1, ++i2) {
        // If the selections are dynamic, they might not have the
        // same objects in both frames, so we need to roll either of
        // the selections until we have the same object to
        // correlate.

        while (i1 != s1.end() && *i1 < *i2) {
          ++i1;
        }

        while (i2 != s2.end() && *i2 < *i1) {
          ++i2;
        }

        if (i1 == s1.end() || i2 == s2.end()) break;

	if (selectionModeRestart_) {
	  
	  int ssf1 = selection1StartFrame_[frame1][i1 - s1.begin()];
	  int ssf2 = uniqueSelections_ ?
	    selection2StartFrame_[frame2][i2 - s2.begin()] :
	    selection1StartFrame_[frame2][i2 - s2.begin()];


	  if (ssf1 != ssf2) break;
	  corrVal = calcCorrVal(frame1, frame2,
				i1 - s1.begin(), i2 - s2.begin());	  
	} else {
	  
	  corrVal = calcCorrVal(frame1, frame2,
				i1 - s1.begin(), i2 - s2.begin());        
	}
        histogram_[timeBin] += corrVal;
        count_[timeBin]++;
      }
    }
  }

  template<typename T>
  void TimeCorrFunc<T>::postCorrelate() {
    T zeroType(0.0);
    for (unsigned int i = 0; i < nTimeBins_; ++i) {
      if (count_[i] > 0) {
        histogram_[i] /= count_[i];
      } else {
        histogram_[i] = zeroType;
      }
    }
  }

  template<typename T>
  void TimeCorrFunc<T>::validateSelection(SelectionManager&) {}

  template<typename T>
  void TimeCorrFunc<T>::writeCorrelate() {
    ofstream ofs(outputFilename_.c_str());

    if (ofs.is_open()) {
      Revision r;

      ofs << "# " << getCorrFuncType() << "\n";
      ofs << "# OpenMD " << r.getFullRevision() << "\n";
      ofs << "# " << r.getBuildDate() << "\n";
      ofs << "# selection script1: \"" << selectionScript1_;
      ofs << "\"\tselection script2: \"" << selectionScript2_ << "\"\n";
      if (!paramString_.empty())
        ofs << "# parameters: " << paramString_ << "\n";
      if (!labelString_.empty())
        ofs << "#time\t" << labelString_ << "\n";
      else
        ofs << "#time\tcorrVal\n";
      
      for (unsigned int i = 0; i < nTimeBins_; ++i) {
	if (count_[i] > 0) 
	  ofs << times_[i] - times_[0] << "\t" << histogram_[i] << "\n";
      }
      
    } else {
      snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
               "TimeCorrFunc::writeCorrelate Error: fail to open %s\n",
               outputFilename_.c_str());
      painCave.isFatal = 1;
      simError();
    }

    ofs.close();
  }

  // Template specialization of writeCorrelate for Vector3d
  template<>
  void TimeCorrFunc<Vector3d>::writeCorrelate() {
    ofstream ofs(outputFilename_.c_str());

    if (ofs.is_open()) {
      Revision r;

      ofs << "# " << getCorrFuncType() << "\n";
      ofs << "# OpenMD " << r.getFullRevision() << "\n";
      ofs << "# " << r.getBuildDate() << "\n";
      ofs << "# selection script1: \"" << selectionScript1_;
      ofs << "\"\tselection script2: \"" << selectionScript2_ << "\"\n";
      if (!paramString_.empty())
        ofs << "# parameters: " << paramString_ << "\n";
      if (!labelString_.empty())
        ofs << "#time\t" << labelString_ << "\n";
      else
        ofs << "#time\tcorrVal\n";

      for (unsigned int i = 0; i < nTimeBins_; ++i) {
	if (count_[i] > 0) {
	
	  ofs << times_[i] - times_[0] << "\t";
	  for (int j = 0; j < 3; j++) {
	    ofs << histogram_[i](j) << '\t';
	  }
	  ofs << '\n';
	  
	}
      }

    } else {
      snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
               "TimeCorrFunc::writeCorrelate Error: fail to open %s\n",
               outputFilename_.c_str());
      painCave.isFatal = 1;
      simError();
    }

    ofs.close();
  }

  // Template specialization of writeCorrelate for Mat3x3d
  template<>
  void TimeCorrFunc<Mat3x3d>::writeCorrelate() {
    ofstream ofs(outputFilename_.c_str());

    if (ofs.is_open()) {
      Revision r;

      ofs << "# " << getCorrFuncType() << "\n";
      ofs << "# OpenMD " << r.getFullRevision() << "\n";
      ofs << "# " << r.getBuildDate() << "\n";
      ofs << "# selection script1: \"" << selectionScript1_;
      ofs << "\"\tselection script2: \"" << selectionScript2_ << "\"\n";
      if (!paramString_.empty())
        ofs << "# parameters: " << paramString_ << "\n";
      if (!labelString_.empty())
        ofs << "#time\t" << labelString_ << "\n";
      else
        ofs << "#time\tcorrVal\n";

      for (unsigned int i = 0; i < nTimeBins_; ++i) {

	if (count_[i] > 0) {
	  
	  ofs << times_[i] - times_[0] << "\t";
	  for (int j = 0; j < 3; j++) {
	    for (int k = 0; k < 3; k++) {
	      ofs << histogram_[i](j, k) << '\t';
	    }
	  }
	  ofs << '\n';
	}
      }

    } else {
      snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
               "TimeCorrFunc::writeCorrelate Error: fail to open %s\n",
               outputFilename_.c_str());
      painCave.isFatal = 1;
      simError();
    }

    ofs.close();
  }

  // it is necessary to keep the constructor definitions here or the code wont
  // be generated and linking issues will occur. Blame templating
  template<typename T>
  CrossCorrFunc<T>::CrossCorrFunc(SimInfo* info, const std::string& filename,
                                  const std::string& sele1,
                                  const std::string& sele2) :
      TimeCorrFunc<T>(info, filename, sele1, sele2) {
    this->autoCorrFunc_ = false;
  }

  template<typename T>
  AutoCorrFunc<T>::AutoCorrFunc(SimInfo* info, const std::string& filename,
                                const std::string& sele1,
                                const std::string& sele2) :
      TimeCorrFunc<T>(info, filename, sele1, sele2) {
    this->autoCorrFunc_ = true;
  }

  template<typename T>
  SystemACF<T>::SystemACF(SimInfo* info, const std::string& filename,
                          const std::string& sele1, const std::string& sele2) :
      AutoCorrFunc<T>(info, filename, sele1, sele2) {
    this->autoCorrFunc_          = true;
    this->doSystemProperties_    = true;
    this->doMolecularProperties_ = false;
    this->doObjectProperties_    = false;
    this->doBondProperties_      = false;
  }

  template<typename T>
  SystemCCF<T>::SystemCCF(SimInfo* info, const std::string& filename,
                          const std::string& sele1, const std::string& sele2) :
      CrossCorrFunc<T>(info, filename, sele1, sele2) {
    this->autoCorrFunc_          = false;
    this->doSystemProperties_    = true;
    this->doMolecularProperties_ = false;
    this->doObjectProperties_    = false;
    this->doBondProperties_      = false;
  }

  template<typename T>
  ObjectACF<T>::ObjectACF(SimInfo* info, const std::string& filename,
                          const std::string& sele1, const std::string& sele2) :
      AutoCorrFunc<T>(info, filename, sele1, sele2) {
    this->autoCorrFunc_          = true;
    this->doSystemProperties_    = false;
    this->doMolecularProperties_ = false;
    this->doObjectProperties_    = true;
    this->doBondProperties_      = false;
  }

  template<typename T>
  ObjectCCF<T>::ObjectCCF(SimInfo* info, const std::string& filename,
                          const std::string& sele1, const std::string& sele2) :
      CrossCorrFunc<T>(info, filename, sele1, sele2) {
    this->autoCorrFunc_          = false;
    this->doSystemProperties_    = false;
    this->doMolecularProperties_ = false;
    this->doObjectProperties_    = true;
    this->doBondProperties_      = false;
  }

  template<typename T>
  MoleculeACF<T>::MoleculeACF(SimInfo* info, const std::string& filename,
                              const std::string& sele1,
                              const std::string& sele2) :
      AutoCorrFunc<T>(info, filename, sele1, sele2) {
    this->autoCorrFunc_          = true;
    this->doSystemProperties_    = false;
    this->doMolecularProperties_ = true;
    this->doObjectProperties_    = false;
    this->doBondProperties_      = false;
  }

  template<typename T>
  MoleculeCCF<T>::MoleculeCCF(SimInfo* info, const std::string& filename,
                              const std::string& sele1,
                              const std::string& sele2) :
      CrossCorrFunc<T>(info, filename, sele1, sele2) {
    this->autoCorrFunc_          = false;
    this->doSystemProperties_    = false;
    this->doMolecularProperties_ = true;
    this->doObjectProperties_    = false;
    this->doBondProperties_      = false;
  }

  template class AutoCorrFunc<RealType>;
  template class TimeCorrFunc<RealType>;
  template class CrossCorrFunc<RealType>;

  template class AutoCorrFunc<Vector3d>;
  template class TimeCorrFunc<Vector3d>;
  template class CrossCorrFunc<Vector3d>;

  template class AutoCorrFunc<Mat3x3d>;
  template class TimeCorrFunc<Mat3x3d>;
  template class CrossCorrFunc<Mat3x3d>;

  template class TimeCorrFunc<DynamicVector<RealType>>;

  template class AutoCorrFunc<Vector<RealType, 4>>;
  template class TimeCorrFunc<Vector<RealType, 4>>;
  template class CrossCorrFunc<Vector<RealType, 4>>;

  template class SystemACF<RealType>;
  template class SystemACF<Vector3d>;
  template class SystemACF<Mat3x3d>;

  template class SystemCCF<RealType>;
  template class SystemCCF<Vector3d>;
  template class SystemCCF<Mat3x3d>;

  template class ObjectACF<RealType>;
  template class ObjectACF<Vector3d>;
  template class ObjectACF<Mat3x3d>;

  template class ObjectCCF<RealType>;
  template class ObjectCCF<Vector3d>;
  template class ObjectCCF<Mat3x3d>;
  template class ObjectCCF<int>;

  template class MoleculeACF<RealType>;
  template class MoleculeACF<Vector3d>;
  template class MoleculeACF<Mat3x3d>;
  template class MoleculeACF<Vector<RealType, 4>>;

  template class MoleculeCCF<RealType>;
  template class MoleculeCCF<Vector3d>;
  template class MoleculeCCF<Mat3x3d>;
}  // namespace OpenMD
