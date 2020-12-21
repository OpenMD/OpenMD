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


  //! Computes a correlation function by scanning a trajectory once to precompute quantities to be correlated

  template<typename T>
  class TimeCorrFunc : public DynamicProperty {
  public:
    TimeCorrFunc(SimInfo* info, const std::string& filename,
                 const std::string& sele1, const std::string& sele2,
                 int storageLayout);
    
    virtual ~TimeCorrFunc(){ }
    virtual void doCorrelate();
    
    const std::string& getCorrFuncType() const {
      return corrFuncType_;
    }
    
    void setCorrFuncType(const std::string& type) {
      corrFuncType_ = type;
    }
    
    void setParameterString(const std::string& params) {
      paramString_ = params;
    }

    void setLabelString(const std::string& label) {
      labelString_ = label;
    }

    
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
    virtual void computeProperty1(int frame) = 0;
    virtual void computeProperty2(int frame) = 0;
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

    int storageLayout_;

    RealType deltaTime_;
    unsigned int nTimeBins_;
    int nFrames_;
    std::vector<T> histogram_;
    std::vector<int> count_;
    std::vector<RealType> times_;
   
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

    std::string corrFuncType_;
    std::string paramString_;
    std::string labelString_;

    std::vector<std::vector<int> > sele1ToIndex_;
    std::vector<std::vector<int> > sele2ToIndex_;

    ProgressBarPtr progressBar_;

  };

  template<typename T>
  class AutoCorrFunc : public TimeCorrFunc<T> {
  public:
    AutoCorrFunc(SimInfo* info, const std::string& filename,
                 const std::string& sele1, const std::string& sele2,
                 int storageLayout);
    
  protected:
    virtual void computeProperty1(int frame) = 0;
    virtual int computeProperty1(int frame, Molecule* mol) = 0;
    virtual int computeProperty1(int frame, StuntDouble* sd) = 0;
    virtual int computeProperty1(int frame, Bond* bond) = 0;
    
    virtual void computeProperty2(int frame) {}
    virtual int computeProperty2(int frame, Molecule* mol) { return -1; }
    virtual int computeProperty2(int frame, StuntDouble* sd) { return -1; }
    virtual int computeProperty2(int frame, Bond* bond) { return -1; }

  };

  template<typename T>
  class CrossCorrFunc : public TimeCorrFunc<T> {
  public:
     CrossCorrFunc(SimInfo* info, const std::string& filename,
                   const std::string& sele1, const std::string& sele2,
                   int storageLayout);

  protected:
    virtual void computeProperty1(int frame) = 0;
    virtual int computeProperty1(int frame, Molecule* mol) = 0;
    virtual int computeProperty1(int frame, StuntDouble* sd) = 0;
    virtual int computeProperty1(int frame, Bond* bond) = 0;
    virtual void computeProperty2(int frame) = 0;
    virtual int computeProperty2(int frame, Molecule* mol) = 0;
    virtual int computeProperty2(int frame, StuntDouble* sd) = 0;
    virtual int computeProperty2(int frame, Bond* bond) = 0;
  };

  template<typename T>
  class SystemACF : public AutoCorrFunc<T> {
  public:
    SystemACF(SimInfo* info, const std::string& filename,
              const std::string& sele1, const std::string& sele2,
              int storageLayout);
  protected:
    virtual void computeProperty1(int frame) = 0;
    virtual T calcCorrVal(int frame1, int frame2) = 0;

    virtual int computeProperty1(int frame, Molecule* mol) { return -1; }
    virtual int computeProperty1(int frame, StuntDouble* sd) { return -1; }
    virtual int computeProperty1(int frame, Bond* bond) { return -1; }

    virtual void computeProperty2(int frame) {}
    virtual int computeProperty2(int frame, Molecule* mol) { return -1; }
    virtual int computeProperty2(int frame, StuntDouble* sd) { return -1; }
    virtual int computeProperty2(int frame, Bond* bond) { return -1; }
    
    T calcCorrVal(int frame1, int frame2, int id1, int id2) { return T(0.0); }

  };

  template<typename T>
  class SystemCCF : public CrossCorrFunc<T> {
  public:
    SystemCCF(SimInfo* info, const std::string& filename,
              const std::string& sele1, const std::string& sele2,
              int storageLayout);
  protected:
    virtual void computeProperty1(int frame) = 0;
    virtual void computeProperty2(int frame) = 0;
    virtual T calcCorrVal(int frame1, int frame2) = 0;

    virtual int computeProperty1(int frame, Molecule* mol) { return -1; }
    virtual int computeProperty1(int frame, StuntDouble* sd) { return -1; }
    virtual int computeProperty1(int frame, Bond* bond) { return -1; }

    virtual int computeProperty2(int frame, Molecule* mol) { return -1; }
    virtual int computeProperty2(int frame, StuntDouble* sd) { return -1; }
    virtual int computeProperty2(int frame, Bond* bond) { return -1; }
    
    T calcCorrVal(int frame1, int frame2, int id1, int id2) { return T(0.0); }

  };
  
  template<typename T>
  class ObjectACF : public AutoCorrFunc<T> {
  public:
    ObjectACF(SimInfo* info, const std::string& filename,
              const std::string& sele1, const std::string& sele2,
              int storageLayout);
  protected:
    virtual int computeProperty1(int frame, StuntDouble* sd) = 0;
    virtual T calcCorrVal(int frame1, int frame2, int id1, int id2) = 0;
    
    virtual void computeProperty1(int frame) { return; }
    virtual int computeProperty1(int frame, Molecule* mol) { return -1; }
    virtual int computeProperty1(int frame, Bond* bond) { return -1; }

    virtual void computeProperty2(int frame) { return; }
    virtual int computeProperty2(int frame, Molecule* mol) { return -1; }
    virtual int computeProperty2(int frame, StuntDouble* sd) { return -1; }
    virtual int computeProperty2(int frame, Bond* bond) { return -1; }
    
    virtual T calcCorrVal(int frame1, int frame2) {return T(0.0);}
  };

  template<typename T>
  class ObjectCCF : public CrossCorrFunc<T> {
  public:
    ObjectCCF(SimInfo* info, const std::string& filename,
              const std::string& sele1, const std::string& sele2,
              int storageLayout);
  protected:
    virtual int computeProperty1(int frame, StuntDouble* sd) = 0;
    virtual int computeProperty2(int frame, StuntDouble* sd) = 0;
    virtual T calcCorrVal(int frame1, int frame2, int id1, int id2) = 0;
    
    virtual void computeProperty1(int frame) {return;}
    virtual int computeProperty1(int frame, Molecule* mol) { return -1; }
    virtual int computeProperty1(int frame, Bond* bond) { return -1; }

    virtual void computeProperty2(int frame) {}
    virtual int computeProperty2(int frame, Molecule* mol) { return -1; }
    virtual int computeProperty2(int frame, Bond* bond) { return -1; }
    
    virtual T calcCorrVal(int frame1, int frame2) {return T(0.0);}
    
  };
  
  template<typename T>
  class MoleculeACF : public AutoCorrFunc<T> {
  public:
    MoleculeACF(SimInfo* info, const std::string& filename,
                const std::string& sele1, const std::string& sele2,
                int storageLayout);
  protected:
    virtual int computeProperty1(int frame, Molecule* mol) = 0;
    virtual T calcCorrVal(int frame1, int frame2, int id1, int id2) = 0;
    
    virtual void computeProperty1(int frame) { return; }
    virtual int computeProperty1(int frame, StuntDouble* sd) { return -1; }
    virtual int computeProperty1(int frame, Bond* bond) { return -1; }
    
    virtual void computeProperty2(int frame) { return; }
    virtual int computeProperty2(int frame, Molecule* mol) { return -1; }
    virtual int computeProperty2(int frame, StuntDouble* sd) { return -1; }
    virtual int computeProperty2(int frame, Bond* bond) { return -1; }
    
    virtual T calcCorrVal(int frame1, int frame2) {return T(0.0);}
  };
  
  template<typename T>
  class MoleculeCCF : public CrossCorrFunc<T> {
  public:
    MoleculeCCF(SimInfo* info, const std::string& filename,
                const std::string& sele1, const std::string& sele2,
                int storageLayout);
  protected:
    virtual int computeProperty1(int frame, Molecule* mol) = 0;
    virtual int computeProperty2(int frame, Molecule* mol) = 0;
    virtual T calcCorrVal(int frame1, int frame2, int id1, int id2) = 0;
    
    virtual void computeProperty1(int frame) {return;}
    virtual int computeProperty1(int frame, StuntDouble* sd) { return -1; }
    virtual int computeProperty1(int frame, Bond* bond) { return -1; }
    
    virtual void computeProperty2(int frame) {}
    virtual int computeProperty2(int frame, StuntDouble* sd) { return -1; }
    virtual int computeProperty2(int frame, Bond* bond) { return -1; }
    
    virtual T calcCorrVal(int frame1, int frame2) {return T(0.0);}
    
  };
  
  
}
#endif
