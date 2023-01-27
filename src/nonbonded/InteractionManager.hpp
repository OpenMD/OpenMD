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

#ifndef NONBONDED_INTERACTIONMANAGER_HPP
#define NONBONDED_INTERACTIONMANAGER_HPP

#include <memory>

#include "brains/SimInfo.hpp"
#include "nonbonded/EAM.hpp"
#include "nonbonded/Electrostatic.hpp"
#include "nonbonded/GB.hpp"
#include "nonbonded/InversePowerSeries.hpp"
#include "nonbonded/LJ.hpp"
#include "nonbonded/MAW.hpp"
#include "nonbonded/Mie.hpp"
#include "nonbonded/Morse.hpp"
#include "nonbonded/RepulsivePower.hpp"
#include "nonbonded/SC.hpp"
#include "nonbonded/Sticky.hpp"
#include "nonbonded/SwitchingFunction.hpp"
#include "types/AtomType.hpp"
// #include "flucq/FluctuatingChargeForces.hpp"

using namespace std;

namespace OpenMD {

  /**
   * @class InteractionManager
   * InteractionManager is responsible for
   * keeping track of the non-bonded interactions (C++)
   */
  class InteractionManager {
  public:
    InteractionManager();
    virtual ~InteractionManager() = default;
    void setSimInfo(SimInfo* info) { info_ = info; }
    void initialize();

    // Fortran support routines

    void doPrePair(InteractionData& idat);
    void doPreForce(SelfData& sdat);
    void doPair(InteractionData& idat);
    void doSkipCorrection(InteractionData& idat);
    void doSelfCorrection(SelfData& sdat);
    void doSurfaceTerm(bool slabGeometry, int axis, RealType& surfacePot);
    void doReciprocalSpaceSum(RealType& recipPot);
    void setCutoffRadius(RealType rCut);
    RealType getSuggestedCutoffRadius(int* atid1);
    RealType getSuggestedCutoffRadius(AtomType* atype);

  private:
    bool initialized_ {false};

    void setupElectrostatics();

    SimInfo* info_ {nullptr};

    // std::shared_ptr<FluctuatingChargeForces> flucq_;

    std::shared_ptr<LJ> lj_ {nullptr};
    std::shared_ptr<GB> gb_ {nullptr};
    std::shared_ptr<Sticky> sticky_ {nullptr};
    std::shared_ptr<EAM> eam_ {nullptr};
    std::shared_ptr<SC> sc_ {nullptr};
    std::shared_ptr<Morse> morse_ {nullptr};
    std::shared_ptr<Electrostatic> electrostatic_ {nullptr};
    std::shared_ptr<RepulsivePower> repulsivePower_ {nullptr};
    std::shared_ptr<Mie> mie_ {nullptr};
    std::shared_ptr<MAW> maw_ {nullptr};
    std::shared_ptr<InversePowerSeries> inversePowerSeries_ {nullptr};

    map<int, AtomType*> typeMap_;
    /**
     * Each pair of atom types can have multiple interactions, so the
     * natural data structures are a map between the pair, and a set
     * of non-bonded interactions:
     *
     *  map<pair<AtomType*, AtomType*>, set<NonBondedInteraction*> >
     * interactions_;
     *
     * Pair creation turns out to be inefficient, and map searching
     * isn't necessary.  Instead of AtomType* sort keys, we now use
     * the AtomType idents (atids) which are ints to access the vector
     * locations.  iHash_ contains largely the same information as
     * interactions_, but in a way that doesn't require a set iterator
     * inside the main pair loop.
     */
    vector<vector<set<NonBondedInteractionPtr>>> interactions_;
    vector<vector<int>> iHash_;

    /* sHash_ contains the self-interaction version of iHash_ */
    vector<int> sHash_;
  };
}  // namespace OpenMD

#endif
