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

#ifndef NONBONDED_ELECTROSTATIC_HPP
#define NONBONDED_ELECTROSTATIC_HPP

#include "brains/ForceField.hpp"
#include "brains/SimInfo.hpp"
#include "flucq/FluctuatingChargeForces.hpp"
#include "math/CubicSpline.hpp"
#include "math/SquareMatrix3.hpp"
#include "nonbonded/NonBondedInteraction.hpp"
#include "primitives/Atom.hpp"
#include "types/AtomType.hpp"

namespace OpenMD {

  struct ElectrostaticAtomData {
    bool is_Charge;
    bool is_Dipole;
    bool is_Quadrupole;
    bool is_Fluctuating;
    bool uses_SlaterIntramolecular;
    RealType fixedCharge;
    RealType hardness;
    RealType electronegativity;
    int slaterN;
    RealType slaterZeta;
    Vector3d dipole;
    Mat3x3d quadrupole;
  };

  enum ElectrostaticSummationMethod {
    esm_HARD,
    esm_SWITCHING_FUNCTION,
    esm_SHIFTED_POTENTIAL,
    esm_SHIFTED_FORCE,
    esm_TAYLOR_SHIFTED,
    esm_REACTION_FIELD,
    esm_EWALD_FULL,
    esm_EWALD_PME, /**< PME  Ewald methods aren't supported yet */
    esm_EWALD_SPME /**< SPME Ewald methods aren't supported yet */
  };

  enum ElectrostaticScreeningMethod { UNDAMPED, DAMPED };

  class Electrostatic : public ElectrostaticInteraction {
  public:
    Electrostatic();
    ~Electrostatic();
    void setForceField(ForceField* ff);
    void setSimulatedAtomTypes(AtomTypeSet& simtypes);
    void setSimInfo(SimInfo* info) { info_ = info; };
    void addType(AtomType* atomType);
    virtual void calcForce(InteractionData& idat);
    virtual void calcSelfCorrection(SelfData& sdat);
    virtual string getName() { return name_; }
    virtual RealType getSuggestedCutoffRadius(
        pair<AtomType*, AtomType*> atypes);
    void setCutoffRadius(RealType rCut);
    void setElectrostaticSummationMethod(ElectrostaticSummationMethod esm);
    void setElectrostaticScreeningMethod(ElectrostaticScreeningMethod sm);
    void setDampingAlpha(RealType alpha);
    void setReactionFieldDielectric(RealType dielectric);
    void calcSurfaceTerm(bool slabGeometry, int axis, RealType& pot);
    void ReciprocalSpaceSum(RealType& pot);

    // Used by EAM to compute local fields:
    RealType getFieldFunction(RealType r);

    // Utility routine
    void getSitePotentials(Atom* a1, Atom* a2, bool excluded, RealType& spot1,
                           RealType& spot2);

  private:
    void initialize();
    string name_;
    bool initialized_;
    bool haveCutoffRadius_;
    bool haveDampingAlpha_;
    bool haveDielectric_;
    bool haveElectroSplines_;

    int nElectro_;
    int nFlucq_;

    set<int>
        Etypes; /**< The set of AtomType idents that are Electrostatic types */
    vector<int>
        Etids; /**< The mapping from AtomType ident -> Electrostatic ident */
    set<int>
        FQtypes; /**< The set of AtomType idents that are fluctuating types */
    vector<int>
        FQtids; /**< The mapping from AtomType ident -> fluctuating ident */
    vector<ElectrostaticAtomData>
        ElectrostaticMap; /**< data about Electrostatic types */
    vector<vector<CubicSplinePtr>>
        Jij; /**< Coulomb integral for two fq types */

    SimInfo* info_ {nullptr};
    ForceField* forceField_;
    FluctuatingChargeForces* flucQ_;
    AtomTypeSet simTypes_;
    RealType cutoffRadius_;
    RealType pre11_;
    RealType pre12_;
    RealType pre22_;
    RealType pre14_;
    RealType pre24_;
    RealType pre44_;
    RealType v01, v11, v21, v22, v31, v32, v41, v42, v43;
    RealType dv01, dv11, dv21, dv22, dv31, dv32, dv41, dv42, dv43;
    RealType v11or, v22or, v31or, v32or, v42or, v43or;
    RealType chargeToC_;
    RealType angstromToM_;
    RealType debyeToCm_;
    int np_;
    ElectrostaticSummationMethod summationMethod_;
    ElectrostaticScreeningMethod screeningMethod_;
    map<string, ElectrostaticSummationMethod> summationMap_;
    map<string, ElectrostaticScreeningMethod> screeningMap_;
    RealType dampingAlpha_;
    RealType dielectric_;
    RealType preRF_;
    RealType selfMult1_;
    RealType selfMult2_;
    RealType selfMult4_;

    CubicSplinePtr v01s;
    CubicSplinePtr v11s;
    CubicSplinePtr v21s;
    CubicSplinePtr v22s;
    CubicSplinePtr v31s;
    CubicSplinePtr v32s;
    CubicSplinePtr v41s;
    CubicSplinePtr v42s;
    CubicSplinePtr v43s;

    ElectrostaticAtomData data1;
    ElectrostaticAtomData data2;
    RealType C_a, C_b;  // Charges
    Vector3d D_a, D_b;  // Dipoles (space-fixed)
    Mat3x3d Q_a, Q_b;   // Quadrupoles (space-fixed)

    RealType ri;                                // Distance utility scalar
    RealType rdDa, rdDb;                        // Dipole utility scalars
    Vector3d rxDa, rxDb;                        // Dipole utility vectors
    RealType rdQar, rdQbr, trQa, trQb;          // Quadrupole utility scalars
    Vector3d Qar, Qbr, rQa, rQb, rxQar, rxQbr;  // Quadrupole utility vectors
    RealType pref;

    RealType DadDb, trQaQb, DadQbr, DbdQar;  // Cross-interaction scalars
    RealType rQaQbr;
    Vector3d DaxDb, DadQb, DbdQa, DaxQbr, DbxQar;  // Cross-interaction vectors
    Vector3d rQaQb, QaQbr, QaxQb, rQaxQbr;
    Mat3x3d QaQb;  // Cross-interaction matrices

    RealType U;      // Potential
    Vector3d F;      // Force
    Vector3d Ta;     // Torque on site a
    Vector3d Tb;     // Torque on site b
    Vector3d Ea;     // Electric field at site a
    Vector3d Eb;     // Electric field at site b
    RealType Pa;     // Site potential at site a
    RealType Pb;     // Site potential at site b
    RealType dUdCa;  // fluctuating charge force at site a
    RealType dUdCb;  // fluctuating charge force at site a

    // Indirect interactions mediated by the reaction field.
    RealType indirect_Pot;  // Potential
    Vector3d indirect_F;    // Force
    Vector3d indirect_Ta;   // Torque on site a
    Vector3d indirect_Tb;   // Torque on site b

    // Excluded potential that is still computed for fluctuating charges
    RealType excluded_Pot;

    RealType rfContrib, coulInt;

    // spline for coulomb integral
    CubicSplinePtr J;
    Vector3d rhat;

    // logicals

    bool a_is_Charge;
    bool a_is_Dipole;
    bool a_is_Quadrupole;
    bool a_is_Fluctuating;
    bool a_uses_SlaterIntra;

    bool b_is_Charge;
    bool b_is_Dipole;
    bool b_is_Quadrupole;
    bool b_is_Fluctuating;
    bool b_uses_SlaterIntra;
  };
}  // namespace OpenMD

#endif
