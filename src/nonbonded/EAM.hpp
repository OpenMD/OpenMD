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

#ifndef NONBONDED_EAM_HPP
#define NONBONDED_EAM_HPP

#include <memory>

#include "nonbonded/NonBondedInteraction.hpp"
#include "nonbonded/Electrostatic.hpp"
#include "types/EAMAdapter.hpp"
#include "brains/ForceField.hpp"
#include "math/Vector3.hpp"
#include "math/CubicSpline.hpp"

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
    RealType rcut;
    bool explicitlySet;
  };

  enum EAMMixingMethod{
    eamJohnson,
    eamDaw,
    eamDream1,
    eamDream2,
    eamUnknownMix
  };
  
  class EAM : public MetallicInteraction {
  public:
    EAM();
    void setForceField(ForceField *ff) {forceField_ = ff;};
    void setElectrostatic(Electrostatic *el) { electrostatic_ = el;};    
    void setSimulatedAtomTypes(set<AtomType*> &simtypes) {simTypes_ = simtypes; initialize();};
    void addType(AtomType* atomType);
    void addExplicitInteraction(AtomType* atype1, AtomType* atype2,
                                RealType dr, int nr,
                                std::vector<RealType> phiAB);
    
    void addExplicitInteraction(AtomType* atype1, AtomType* atype2,
                                RealType re, RealType alpha, RealType beta,
                                RealType A, RealType B, RealType kappa,
                                RealType lambda);

    RealType fastPower(RealType x, int y);
    RealType ZhouPhiCoreCore(RealType r, RealType re, RealType A,
                             RealType alpha, RealType kappa);
    RealType ZhouPhiCoreValence(RealType r, RealType re, RealType B,
                                RealType beta, RealType lambda);
    
    RealType ZhouPhi(RealType r, RealType re, RealType A, RealType B,
                     RealType alpha, RealType beta, RealType kappa,
                     RealType lambda);
    
    RealType ZhouRho(RealType r, RealType re, RealType fe,
                     RealType beta, RealType lambda);
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
                                std::vector<RealType> F,
                                RealType F3plus, RealType F3minus,
                                RealType Fe, RealType eta);
    RealType Zhou2005OxygenFunctional(RealType rho,
                                      std::vector<RealType> OrhoLimits,
                                      std::vector<RealType> OrhoE,
                                      std::vector<std::vector<RealType> > OF);
    RealType RoseFunctional(RealType rho, RealType rhoe, RealType F0);
    
                                      
    void calcDensity(InteractionData &idat);
    void calcFunctional(SelfData &sdat);
    void calcForce(InteractionData &idat);
    
    virtual string getName() { return name_; }
    virtual int getHash() { return EAM_INTERACTION; }
    virtual RealType getSuggestedCutoffRadius(pair<AtomType*,AtomType*> atypes);
    void setCutoffRadius( RealType rCut );

  private:
    void initialize();
    CubicSplinePtr getPhi(AtomType* atomType1, AtomType* atomType2);
    
    bool initialized_;
    bool haveCutoffRadius_;
    set<int> EAMtypes;         /**< The set of AtomType idents that are EAM types */
    vector<int> EAMtids;       /**< The mapping from AtomType ident -> EAM type ident */
    vector<EAMAtomData> EAMdata; /**< The EAM atomic data indexed by EAM type ident */
    vector<vector<EAMInteractionData> > MixingMap;  /**< The mixing parameters between two EAM types */
    int nEAM_;

    ForceField* forceField_;
    Electrostatic* electrostatic_;
    set<AtomType*> simTypes_;
    RealType pre11_;
    RealType eamRcut_;
    RealType oss_;
    Vector3d rhat;

    EAMMixingMethod mixMeth_;
    string name_;

  };
}


#endif
