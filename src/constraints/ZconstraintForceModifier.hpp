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

#ifndef CONSTRAINTS_ZCONSTRAINTFORCEMODIFIER_HPP
#define CONSTRAINTS_ZCONSTRAINTFORCEMODIFIER_HPP

#include <list>
#include <string>
#include <vector>

#include "brains/ForceModifier.hpp"
#include "constraints/ZconsStruct.hpp"
#include "io/ZConsWriter.hpp"

namespace OpenMD {

  class ZConstraintForceModifier : public ForceModifier {
  public:
    ZConstraintForceModifier(SimInfo* info);
    ~ZConstraintForceModifier();

    void modifyForces() override;

    RealType getZConsTime() { return zconsTime_; }
    std::string getZConsOutput() { return zconsOutput_; }

  private:
    RealType getZTargetPos(int index);
    void update();
    void calcTotalMassMovingZMols();
    bool isZMol(Molecule* mol);
    void updateZPos();
    bool checkZConsState();
    void zeroVelocity();
    bool haveFixedZMols();
    void doZconstraintForce();
    RealType getZFOfFixedZMols(Molecule* mol, StuntDouble* sd,
                               RealType totalForce);
    RealType getZFOfMovingMols(Molecule* mol, RealType totalForce);
    bool haveMovingZMols();
    void doHarmonic();
    RealType getHFOfFixedZMols(Molecule* mol, StuntDouble* sd,
                               RealType totalForce);
    RealType getHFOfUnconsMols(Molecule* mol, RealType totalForce);

    // void updateCantPos();

    std::list<ZconstraintMol> movingZMols_; /**< moving zconstraint molecules */
    std::list<ZconstraintMol> fixedZMols_;  /**< fixed zconstraint molecules */
    std::vector<Molecule*> unzconsMols_;    /**< free molecules */

    RealType zconsTime_;
    std::string zconsOutput_;
    RealType zconsTol_;
    bool usingSMD_;
    RealType zconsFixingTime_;
    RealType zconsGap_;
    bool usingZconsGap_;
    RealType dt_;

    const static int whichDirection = 2;

    std::map<int, ZconstraintParam> allZMolIndices_;

    Snapshot* currSnapshot_;
    RealType currZconsTime_;

    RealType totMassMovingZMols_;
    RealType totMassUnconsMols_; /**< mass of unconstrained molecules
                                     in the system (never changes) */

    ZConsWriter* fzOut;
    const RealType infiniteTime;
  };
}  // namespace OpenMD

#endif  // CONSTRAINTS_ZCONSTRAINTFORCEMODIFIER_HPP
