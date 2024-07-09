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

#ifndef APPLICATIONS_DYNAMICPROPS_TIMECORRFUNC_HPP
#define APPLICATIONS_DYNAMICPROPS_TIMECORRFUNC_HPP

#include <string>
#include <vector>

#include "applications/dynamicProps/DynamicProperty.hpp"
#include "brains/SimInfo.hpp"
#include "io/DumpReader.hpp"
#include "primitives/StuntDouble.hpp"
#include "selection/SelectionEvaluator.hpp"
#include "selection/SelectionManager.hpp"
#include "utils/ProgressBar.hpp"

namespace OpenMD {

  //! Computes a correlation function by scanning a trajectory once to
  //! precompute quantities to be correlated

  template<typename T>
  class TimeCorrFunc : public DynamicProperty {
  public:
    TimeCorrFunc(SimInfo* info, const std::string& filename,
                 const std::string& sele1, const std::string& sele2);

    virtual ~TimeCorrFunc() { delete reader_; }
    virtual void doCorrelate();

    const std::string& getCorrFuncType() const { return corrFuncType_; }

    void setCorrFuncType(const std::string& type) { corrFuncType_ = type; }

    void setParameterString(const std::string& params) {
      paramString_ = params;
    }

    void setLabelString(const std::string& label) { labelString_ = label; }

  protected:
    virtual void preCorrelate();
    virtual void correlation();
    virtual void postCorrelate();
    virtual void computeFrame(int frame);
    virtual void validateSelection(SelectionManager& seleMan);
    virtual void correlateFrames(int frame1, int frame2, int timeBin);
    virtual void writeCorrelate();

    // The pure virtual functions that must be implemented.

    // For System Properties:
    virtual void computeProperty1(int frame)      = 0;
    virtual void computeProperty2(int frame)      = 0;
    virtual T calcCorrVal(int frame1, int frame2) = 0;
    // For Molecular Properties
    virtual int computeProperty1(int frame, Molecule* mol) = 0;
    virtual int computeProperty2(int frame, Molecule* mol) = 0;
    // For Object
    virtual int computeProperty1(int frame, StuntDouble* sd) = 0;
    virtual int computeProperty2(int frame, StuntDouble* sd) = 0;
    // For Bond properties
    virtual int computeProperty1(int frame, Bond* b) = 0;
    virtual int computeProperty2(int frame, Bond* b) = 0;
    // For everything except System properties:
    virtual T calcCorrVal(int frame1, int frame2, int id1, int id2) = 0;

    RealType deltaTime_;
    unsigned int nTimeBins_;
    int nFrames_;
    std::vector<T> histogram_;
    std::vector<int> count_;
    std::vector<RealType> times_;
    RealType dtMean_{};
    RealType dtSigma_{};

    SimInfo* info_ {nullptr};
    DumpReader* reader_;
    std::string dumpFilename_;
    SelectionManager seleMan1_;
    SelectionManager seleMan2_;

    Snapshot* currentSnapshot_;

    std::string selectionScript1_;
    std::string selectionScript2_;

    SelectionEvaluator evaluator1_;
    SelectionEvaluator evaluator2_;

    bool uniqueSelections_;
    bool autoCorrFunc_;
    bool doSystemProperties_;
    bool doMolecularProperties_;
    bool doObjectProperties_;
    bool doAtomicProperties_;
    bool doBondProperties_;
    bool allowTimeFuzz_;

    std::string corrFuncType_;
    std::string paramString_;
    std::string labelString_;

    std::vector<std::vector<int>> sele1ToIndex_;
    std::vector<std::vector<int>> sele2ToIndex_;

    ProgressBarPtr progressBar_;
  };

  template<typename T>
  class AutoCorrFunc : public TimeCorrFunc<T> {
  public:
    AutoCorrFunc(SimInfo* info, const std::string& filename,
                 const std::string& sele1, const std::string& sele2);

  protected:
    virtual void computeProperty1(int frame)                 = 0;
    virtual int computeProperty1(int frame, Molecule* mol)   = 0;
    virtual int computeProperty1(int frame, StuntDouble* sd) = 0;
    virtual int computeProperty1(int frame, Bond* bond)      = 0;

    virtual void computeProperty2(int) {}
    virtual int computeProperty2(int, Molecule*) { return -1; }
    virtual int computeProperty2(int, StuntDouble*) { return -1; }
    virtual int computeProperty2(int, Bond*) { return -1; }
  };

  template<typename T>
  class CrossCorrFunc : public TimeCorrFunc<T> {
  public:
    CrossCorrFunc(SimInfo* info, const std::string& filename,
                  const std::string& sele1, const std::string& sele2);

  protected:
    virtual void computeProperty1(int frame)                 = 0;
    virtual int computeProperty1(int frame, Molecule* mol)   = 0;
    virtual int computeProperty1(int frame, StuntDouble* sd) = 0;
    virtual int computeProperty1(int frame, Bond* bond)      = 0;
    virtual void computeProperty2(int frame)                 = 0;
    virtual int computeProperty2(int frame, Molecule* mol)   = 0;
    virtual int computeProperty2(int frame, StuntDouble* sd) = 0;
    virtual int computeProperty2(int frame, Bond* bond)      = 0;
  };

  template<typename T>
  class SystemACF : public AutoCorrFunc<T> {
  public:
    SystemACF(SimInfo* info, const std::string& filename,
              const std::string& sele1, const std::string& sele2);

  protected:
    virtual void computeProperty1(int frame)      = 0;
    virtual T calcCorrVal(int frame1, int frame2) = 0;

    virtual int computeProperty1(int, Molecule*) { return -1; }
    virtual int computeProperty1(int, StuntDouble*) { return -1; }
    virtual int computeProperty1(int, Bond*) { return -1; }

    virtual void computeProperty2(int) {}
    virtual int computeProperty2(int, Molecule*) { return -1; }
    virtual int computeProperty2(int, StuntDouble*) { return -1; }
    virtual int computeProperty2(int, Bond*) { return -1; }

    T calcCorrVal(int, int, int, int) { return T(0.0); }
  };

  template<typename T>
  class SystemCCF : public CrossCorrFunc<T> {
  public:
    SystemCCF(SimInfo* info, const std::string& filename,
              const std::string& sele1, const std::string& sele2);

  protected:
    virtual void computeProperty1(int frame)      = 0;
    virtual void computeProperty2(int frame)      = 0;
    virtual T calcCorrVal(int frame1, int frame2) = 0;

    virtual int computeProperty1(int, Molecule*) { return -1; }
    virtual int computeProperty1(int, StuntDouble*) { return -1; }
    virtual int computeProperty1(int, Bond*) { return -1; }

    virtual int computeProperty2(int, Molecule*) { return -1; }
    virtual int computeProperty2(int, StuntDouble*) { return -1; }
    virtual int computeProperty2(int, Bond*) { return -1; }

    T calcCorrVal(int, int, int, int) { return T(0.0); }
  };

  template<typename T>
  class ObjectACF : public AutoCorrFunc<T> {
  public:
    ObjectACF(SimInfo* info, const std::string& filename,
              const std::string& sele1, const std::string& sele2);

  protected:
    virtual int computeProperty1(int frame, StuntDouble* sd)        = 0;
    virtual T calcCorrVal(int frame1, int frame2, int id1, int id2) = 0;

    virtual void computeProperty1(int) { return; }
    virtual int computeProperty1(int, Molecule*) { return -1; }
    virtual int computeProperty1(int, Bond*) { return -1; }

    virtual void computeProperty2(int) { return; }
    virtual int computeProperty2(int, Molecule*) { return -1; }
    virtual int computeProperty2(int, StuntDouble*) { return -1; }
    virtual int computeProperty2(int, Bond*) { return -1; }

    virtual T calcCorrVal(int, int) { return T(0.0); }
  };

  template<typename T>
  class ObjectCCF : public CrossCorrFunc<T> {
  public:
    ObjectCCF(SimInfo* info, const std::string& filename,
              const std::string& sele1, const std::string& sele2);

  protected:
    virtual int computeProperty1(int frame, StuntDouble* sd)        = 0;
    virtual int computeProperty2(int frame, StuntDouble* sd)        = 0;
    virtual T calcCorrVal(int frame1, int frame2, int id1, int id2) = 0;

    virtual void computeProperty1(int) { return; }
    virtual int computeProperty1(int, Molecule*) { return -1; }
    virtual int computeProperty1(int, Bond*) { return -1; }

    virtual void computeProperty2(int) {}
    virtual int computeProperty2(int, Molecule*) { return -1; }
    virtual int computeProperty2(int, Bond*) { return -1; }

    virtual T calcCorrVal(int, int) { return T(0.0); }
  };

  template<typename T>
  class MoleculeACF : public AutoCorrFunc<T> {
  public:
    MoleculeACF(SimInfo* info, const std::string& filename,
                const std::string& sele1, const std::string& sele2);

  protected:
    virtual int computeProperty1(int frame, Molecule* mol)          = 0;
    virtual T calcCorrVal(int frame1, int frame2, int id1, int id2) = 0;

    virtual void computeProperty1(int) { return; }
    virtual int computeProperty1(int, StuntDouble*) { return -1; }
    virtual int computeProperty1(int, Bond*) { return -1; }

    virtual void computeProperty2(int) { return; }
    virtual int computeProperty2(int, Molecule*) { return -1; }
    virtual int computeProperty2(int, StuntDouble*) { return -1; }
    virtual int computeProperty2(int, Bond*) { return -1; }

    virtual T calcCorrVal(int, int) { return T(0.0); }
  };

  template<typename T>
  class MoleculeCCF : public CrossCorrFunc<T> {
  public:
    MoleculeCCF(SimInfo* info, const std::string& filename,
                const std::string& sele1, const std::string& sele2);

  protected:
    virtual int computeProperty1(int frame, Molecule* mol)          = 0;
    virtual int computeProperty2(int frame, Molecule* mol)          = 0;
    virtual T calcCorrVal(int frame1, int frame2, int id1, int id2) = 0;

    virtual void computeProperty1(int) { return; }
    virtual int computeProperty1(int, StuntDouble*) { return -1; }
    virtual int computeProperty1(int, Bond*) { return -1; }

    virtual void computeProperty2(int) {}
    virtual int computeProperty2(int, StuntDouble*) { return -1; }
    virtual int computeProperty2(int, Bond*) { return -1; }

    virtual T calcCorrVal(int, int) { return T(0.0); }
  };
}  // namespace OpenMD

#endif
