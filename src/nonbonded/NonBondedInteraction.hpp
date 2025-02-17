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

#ifndef NONBONDED_NONBONDEDINTERACTION_HPP
#define NONBONDED_NONBONDEDINTERACTION_HPP

#include <memory>

#include "math/SquareMatrix3.hpp"
#include "types/AtomType.hpp"

using namespace std;
namespace OpenMD {
  /**
   * The InteractionFamily enum.
   *
   * This is used to sort different types of non-bonded interaction
   * and to prevent multiple interactions in the same family from
   * being applied to any given pair of atom types.
   */
  enum InteractionFamily {
    NO_FAMILY = 0, /**< No family defined */
    VANDERWAALS_FAMILY =
        1, /**< Long-range dispersion and short-range pauli repulsion */
    ELECTROSTATIC_FAMILY = 2, /**< Coulombic and point-multipole interactions */
    METALLIC_EMBEDDING_FAMILY =
        3, /**< Transition metal interactions involving electron density */
    METALLIC_PAIR_FAMILY =
        4, /**< Transition metal interactions involving pair potentials */
    HYDROGENBONDING_FAMILY = 5, /**< Short-range directional interactions */
    BONDED_FAMILY = 6, /**< directly bonded 1-2, 1-3, or 1-4 interactions */
    N_INTERACTION_FAMILIES = 7
  };

  /**
   * Boolean flags for the iHash_ and sHash_ data structures. These
   * are used to greatly increase the speed of looking up the
   * low-level interaction for any given pair:
   */
  const static int ELECTROSTATIC_INTERACTION      = (1 << 0);
  const static int LJ_INTERACTION                 = (1 << 1);
  const static int EAM_INTERACTION                = (1 << 2);
  const static int SC_INTERACTION                 = (1 << 3);
  const static int STICKY_INTERACTION             = (1 << 4);
  const static int GB_INTERACTION                 = (1 << 5);
  const static int MORSE_INTERACTION              = (1 << 6);
  const static int REPULSIVEPOWER_INTERACTION     = (1 << 7);
  const static int MAW_INTERACTION                = (1 << 8);
  const static int MIE_INTERACTION                = (1 << 9);
  const static int BUCKINGHAM_INTERACTION         = (1 << 10);
  const static int INVERSEPOWERSERIES_INTERACTION = (1 << 11);

  using potVec = Vector<RealType, N_INTERACTION_FAMILIES>;

  /**
   * The InteractionData struct.
   *
   * This is used to pass data and references to data to specific non-bonded
   * interactions for force calculations. Not all of the struct
   * members are utilized by any given interaction.
   */
  struct InteractionData {
    int atid1 = -1;   /**< atomType ident for atom 1 */
    int atid2 = -1;   /**< atomType ident for atom 2 */
    Vector3d d {};    /**< interatomic vector (already wrapped into box) */
    RealType rij {};  /**< interatomic separation */
    RealType r2 {};   /**< square of rij */
    RealType rcut {}; /**< cutoff radius for this interaction */
    bool shiftedPot {false}; /**< shift the potential up inside the cutoff? */
    bool shiftedForce {
        false};      /**< shifted forces smoothly inside the cutoff? */
    RealType sw {};  /**< switching function value at rij */
    int topoDist {}; /**< topological distance between atoms */
    bool excluded {
        false}; /**< is this excluded from direct pairwise interactions? (some
                   indirect interactions may still apply) */
    bool sameRegion {
        false}; /**< are these atoms specified to be in the same region? */
    RealType vdwMult {};     /**< multiplier for van der Waals interactions */
    RealType electroMult {}; /**< multiplier for electrostatic interactions */
    potVec pot {};           /**< total potential */
    potVec excludedPot {};   /**< potential energy excluded from the overall
                                calculation */
    RealType vpair {};       /**< pair potential */
    bool doParticlePot {false}; /**< should we bother with the particle pot? */
    bool doElectricField {
        false}; /**< should we bother with the electric field? */
    bool doSitePotential {
        false}; /**< should we bother with electrostatic site potential */
    bool isSelected {false};    /**< one of the particles has been selected for
                                   selection potential */
    potVec selePot {};          /**< potential energies of the selected sites */
    RealType particlePot1 {};   /**< particle potential for atom1 */
    RealType particlePot2 {};   /**< particle potential for atom2 */
    Vector3d f1 {};             /**< force between the two atoms */
    RotMat3x3d A1 {};           /**< rotation matrix of first atom */
    RotMat3x3d A2 {};           /**< rotation matrix of second atom */
    Vector3d D_1 {};            /**< dipole vector of first atom */
    Vector3d D_2 {};            /**< dipole vector of first atom */
    Mat3x3d Q_1 {};             /**< quadrupole tensor of first atom */
    Mat3x3d Q_2 {};             /**< quadrupole tensor of first atom */
    Vector3d t1 {};             /**< torque on first atom */
    Vector3d t2 {};             /**< torque on second atom */
    RealType rho1 {};           /**< total electron density at first atom */
    RealType rho2 {};           /**< total electron density at second atom */
    RealType frho1 {};          /**< density functional at first atom */
    RealType frho2 {};          /**< density functional at second atom */
    RealType dfrho1 {};         /**< derivative of functional for atom 1 */
    RealType dfrho2 {};         /**< derivative of functional for atom 2 */
    RealType flucQ1 {};         /**< fluctuating charge on atom1 */
    RealType flucQ2 {};         /**< fluctuating charge on atom2 */
    RealType dVdFQ1 {};         /**< fluctuating charge force on atom1 */
    RealType dVdFQ2 {};         /**< fluctuating charge force on atom2 */
    Vector3d eField1 {};        /**< electric field on first atom */
    Vector3d eField2 {};        /**< electric field on second atom */
    RealType skippedCharge1 {}; /**< charge skipped on atom1 in pairwise
                                   interaction loop with atom2 */
    RealType skippedCharge2 {}; /**< charge skipped on atom2 in pairwise
                                   interaction loop with atom1 */
    RealType sPot1 {};          /**< site potential on first atom */
    RealType sPot2 {};          /**< site potential on second atom */
  };

  /**
   * The SelfData struct.
   *
   * This is used to pass data for the self-interaction or
   * derived information on a single atom after a pass through all
   * other interactions. This is used by electrostatic methods that
   * have long-range corrections involving interactions with a medium
   * or a boundary and also by specific metal interactions for
   * electron density functional calculations. Not all of the struct
   * members are utilized by any given self interaction.
   */
  struct SelfData {
    int atid {}; /**< atomType ident for the atom */
    RealType
        skippedCharge {}; /**< charge skipped in pairwise interaction loop */
    potVec selfPot {};    /**< total potential (including embedding energy) */
    bool doParticlePot {false}; /**< should we bother with the particle pot? */
    RealType
        particlePot {};    /**< contribution to potential from this particle */
    RealType rho {};       /**< electron density */
    RealType frho {};      /**< value of density functional for atom */
    RealType dfrhodrho {}; /**< derivative of density functional for atom */
    RealType flucQ {};     /**< current value of atom's fluctuating charge */
    RealType flucQfrc {};  /**< fluctuating charge derivative */
    bool isSelected {
        false}; /**< this site has been selected for selection potential */
    potVec selePot {}; /**< potential energy of the selected site */
  };

  /**
   * The basic interface for non-bonded interactions.
   */
  class NonBondedInteraction {
  public:
    NonBondedInteraction() {}
    virtual ~NonBondedInteraction() {}
    virtual void calcForce(InteractionData& idat) = 0;
    virtual InteractionFamily getFamily()         = 0;
    virtual int getHash()                         = 0;
    virtual RealType getSuggestedCutoffRadius(
        pair<AtomType*, AtomType*> atypes) = 0;
    virtual string getName()               = 0;
  };

  using NonBondedInteractionPtr = std::shared_ptr<NonBondedInteraction>;

  /**
   * The basic interface for van der Waals interactions.
   */
  class VanDerWaalsInteraction : public NonBondedInteraction {
  public:
    VanDerWaalsInteraction() : NonBondedInteraction() {}
    virtual ~VanDerWaalsInteraction() {}
    virtual InteractionFamily getFamily() { return VANDERWAALS_FAMILY; }
  };

  /**
   * The basic interface for electrostatic interactions.
   */
  class ElectrostaticInteraction : public NonBondedInteraction {
  public:
    ElectrostaticInteraction() : NonBondedInteraction() {}
    virtual ~ElectrostaticInteraction() {}
    virtual void calcSelfCorrection(SelfData& sdat) = 0;
    virtual InteractionFamily getFamily() { return ELECTROSTATIC_FAMILY; }
    virtual int getHash() { return ELECTROSTATIC_INTERACTION; }
  };

  /**
   * The basic interface for metallic interactions.
   */
  class MetallicInteraction : public NonBondedInteraction {
  public:
    MetallicInteraction() : NonBondedInteraction() {}
    virtual ~MetallicInteraction() {}
    virtual void calcDensity(InteractionData& idat) = 0;
    virtual void calcFunctional(SelfData& sdat)     = 0;
    virtual InteractionFamily getFamily() { return METALLIC_EMBEDDING_FAMILY; }
  };

  /**
   * The basic interface for hydrogen bonding interactions.
   */
  class HydrogenBondingInteraction : public NonBondedInteraction {
  public:
    HydrogenBondingInteraction() : NonBondedInteraction() {}
    virtual ~HydrogenBondingInteraction() {}
    virtual InteractionFamily getFamily() { return HYDROGENBONDING_FAMILY; }
  };
}  // namespace OpenMD

#endif
