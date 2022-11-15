/*
 * Copyright (c) 2004-2022, The University of Notre Dame. All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 * 1. Redistributions of source code must retain the above copyright notice,
 * this list of conditions and the following disclaimer.
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

#ifndef NONBONDED_EAM_HPP
#define NONBONDED_EAM_HPP

#include <memory>

#include "brains/ForceField.hpp"
#include "math/CubicSpline.hpp"
#include "math/Vector3.hpp"
#include "nonbonded/Electrostatic.hpp"
#include "nonbonded/NonBondedInteraction.hpp"
#include "types/EAMAdapter.hpp"

namespace OpenMD {

  struct EAMAtomData {
    CubicSplinePtr rho;
    CubicSplinePtr F;
    CubicSplinePtr Z;
    CubicSplinePtr phiCC;
    CubicSplinePtr phiCV;
    RealType rcut;
    RealType nValence;
    RealType nMobile;
    bool isFluctuatingCharge;
  };

  struct EAMInteractionData {
    CubicSplinePtr phi;
    CubicSplinePtr phiCC;
    RealType Ci = 0.0;  // to zero out the CV interaction in the eam-explict
                        // interaction for EAMTable and EAMZhou
    RealType Cj = 0.0;  // to zero out the CV interaction in the eam-explict
                        // interaction for EAMTable and EAMZhou
    RealType rcut;
    bool explicitlySet;
  };

  enum EAMMixingMethod {
    eamJohnson,
    eamDaw,
    eamDream1,
    eamDream2,
    eamUnknownMix
  };

  class EAM : public MetallicInteraction {
  public:
    EAM();
    void setForceField(ForceField* ff) { forceField_ = ff; };
    void setElectrostatic(Electrostatic* el) { electrostatic_ = el; };
    void setSimulatedAtomTypes(AtomTypeSet& simtypes) {
      simTypes_ = simtypes;
      initialize();
    };
    void addType(AtomType* atomType);
    void addExplicitInteraction(AtomType* atype1, AtomType* atype2, RealType dr,
                                int nr, std::vector<RealType> phiAB);

    void addExplicitInteraction(AtomType* atype1, AtomType* atype2, RealType re,
                                RealType alpha, RealType beta, RealType A,
                                RealType B, RealType kappa, RealType lambda);

    void addExplicitInteraction(AtomType* atype1, AtomType* atype2, RealType re,
                                RealType alpha, RealType A, RealType Ci,
                                RealType Cj);

    RealType fastPower(RealType x, int y);

    // RealType ZhouPhiCoreCore(RealType r, RealType re, RealType A,
    // RealType alpha, RealType kappa);
    // RealType ZhouPhiCoreValence(RealType r, RealType re, RealType B,
    // RealType beta, RealType lambda);
    //
    // RealType ZhouPhi(RealType r, RealType re, RealType A, RealType B,
    // RealType alpha, RealType beta, RealType kappa,
    // RealType lambda);
    //
    // RealType ZhouRho(RealType r, RealType re, RealType fe, RealType beta,
    // RealType lambda);

    RealType PhiCoreCore(RealType r, RealType re, RealType A, RealType alpha,
                         RealType kappa);
    RealType PhiCoreValence(RealType r, RealType re, RealType B, RealType beta,
                            RealType lambda);

    RealType Phi(RealType r, RealType re, RealType A, RealType B,
                 RealType alpha, RealType beta, RealType kappa,
                 RealType lambda);

    RealType Rho(RealType r, RealType re, RealType fe, RealType beta,
                 RealType lambda);
    RealType gFunc(RealType q, RealType nV, RealType nM);
    RealType gPrime(RealType q, RealType nV, RealType nM);
    RealType Zhou2001Functional(RealType rho, RealType rhoe,
                                std::vector<RealType> Fn,
                                std::vector<RealType> F, RealType Fe,
                                RealType eta);
    RealType Zhou2004Functional(RealType rho, RealType rhoe, RealType rhos,
                                std::vector<RealType> Fn,
                                std::vector<RealType> F, RealType Fe,
                                RealType eta, RealType rhol, RealType rhoh);
    RealType Zhou2005Functional(RealType rho, RealType rhoe, RealType rhos,
                                std::vector<RealType> Fn,
                                std::vector<RealType> F, RealType F3plus,
                                RealType F3minus, RealType Fe, RealType eta);
    RealType Zhou2005OxygenFunctional(RealType rho,
                                      std::vector<RealType> OrhoLimits,
                                      std::vector<RealType> OrhoE,
                                      std::vector<std::vector<RealType>> OF);
    RealType RoseFunctional(RealType rho, RealType rhoe, RealType F0);

    void calcDensity(InteractionData& idat);
    void calcFunctional(SelfData& sdat);
    void calcForce(InteractionData& idat);

    virtual string getName() { return name_; }
    virtual int getHash() { return EAM_INTERACTION; }
    virtual RealType getSuggestedCutoffRadius(
        pair<AtomType*, AtomType*> atypes);
    void setCutoffRadius(RealType rCut);

  private:
    void initialize();
    CubicSplinePtr getPhi(AtomType* atomType1, AtomType* atomType2);

    bool initialized_;
    bool haveCutoffRadius_;
    set<int> EAMtypes; /**< The set of AtomType idents that are EAM types */
    vector<int>
        EAMtids; /**< The mapping from AtomType ident -> EAM type ident */
    vector<EAMAtomData>
        EAMdata; /**< The EAM atomic data indexed by EAM type ident */
    vector<vector<EAMInteractionData>>
        MixingMap; /**< The mixing parameters between two EAM types */
    int nEAM_;

    ForceField* forceField_;
    Electrostatic* electrostatic_;
    AtomTypeSet simTypes_;
    RealType pre11_;
    RealType eamRcut_;
    // RealType oss_;
    Vector3d rhat;

    EAMMixingMethod mixMeth_;
    string name_;
  };
}  // namespace OpenMD

#endif
