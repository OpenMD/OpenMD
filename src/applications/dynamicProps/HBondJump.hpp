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

#ifndef APPLICATIONS_DYNAMICPROPS_HBONDJUMP_HPP
#define APPLICATIONS_DYNAMICPROPS_HBONDJUMP_HPP

#include "applications/dynamicProps/TimeCorrFunc.hpp"

#define HONKING_LARGE_VALUE 1.0e10

namespace OpenMD {

  class HBondJump : public TimeCorrFunc<RealType> {
  public:
    HBondJump(SimInfo* info, const std::string& filename,
              const std::string& sele1, const std::string& sele2, double OOCut,
              double thetaCut, double OHCut);

  protected:
    virtual void correlation();
    virtual void computeFrame(int istep);

    virtual void computeProperty1(int) { return; }
    virtual void computeProperty2(int) { return; }
    virtual int computeProperty1(int, Molecule*) { return -1; }
    virtual int computeProperty1(int, StuntDouble*) { return -1; }
    virtual int computeProperty1(int, Bond*) { return -1; }
    virtual int computeProperty2(int, Molecule*) { return -1; }
    virtual int computeProperty2(int, StuntDouble*) { return -1; }
    virtual int computeProperty2(int, Bond*) { return -1; }

    virtual RealType calcCorrVal(int, int, int, int) { return 0.0; }
    virtual RealType calcCorrVal(int, int) { return 0.0; }

    virtual void postCorrelate();

    virtual int registerHydrogen(int frame, int hIndex);
    virtual void findHBonds(int frame);
    bool isHBond(Vector3d donorPos, Vector3d hPos, Vector3d acceptorPos);
    void registerHydrogenBond(int frame, int index, int hIndex, int aIndex);
    void processNonOverlapping(int frame, SelectionManager& sman1,
                               SelectionManager& sman2);
    void processOverlapping(int frame, SelectionManager& sman);

    std::vector<std::vector<int>> GIDtoH_;
    std::vector<std::vector<int>> hydrogen_;
    std::vector<std::vector<int>> acceptor_;
    std::vector<std::vector<int>> lastAcceptor_;
    std::vector<std::vector<bool>> selected_;
    std::vector<std::vector<int>> acceptorStartFrame_;

    RealType OOCut_;
    RealType thetaCut_;
    RealType OHCut_;

    SelectionManager sele1_minus_common_;
    SelectionManager sele2_minus_common_;
    SelectionManager common_;
  };

  class HBondJumpZ : public HBondJump {
  public:
    HBondJumpZ(SimInfo* info, const std::string& filename,
               const std::string& sele1, const std::string& sele2, double OOCut,
               double thetaCut, double OHCut, int nZbins, int axis = 2);
    virtual int registerHydrogen(int frame, int hIndex);
    virtual void findHBonds(int frame);
    virtual void correlation();
    virtual void postCorrelate();
    virtual void writeCorrelate();

  private:
    std::vector<std::vector<RealType>> histogram_;
    std::vector<std::vector<int>> counts_;
    std::vector<std::vector<int>> zbin_;
    unsigned int nZBins_;
    int axis_;
    std::string axisLabel_;
  };

  class HBondJumpR : public HBondJump {
  public:
    HBondJumpR(SimInfo* info, const std::string& filename,
               const std::string& sele1, const std::string& sele2,
               const std::string& sele3, double OOCut, RealType thetaCut,
               RealType OHCut, RealType len, int nRbins);
    virtual int registerHydrogen(int frame, int hIndex);
    virtual void findHBonds(int frame);
    virtual void correlation();
    virtual void postCorrelate();
    virtual void writeCorrelate();

  private:
    RealType len_;
    unsigned int nRBins_;
    RealType deltaR_;
    SelectionManager seleMan3_;
    std::string selectionScript3_;
    SelectionEvaluator evaluator3_;

    std::vector<std::vector<RealType>> histogram_;
    std::vector<std::vector<int>> counts_;
    std::vector<std::vector<int>> rbin_;
  };
}  // namespace OpenMD

#endif
