/*
 * Copyright (c) 2005 The University of Notre Dame. All Rights Reserved.
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
 * [1]  Meineke, et al., J. Comp. Chem. 26, 252-271 (2005).             
 * [2]  Fennell & Gezelter, J. Chem. Phys. 124, 234104 (2006).          
 * [3]  Sun, Lin & Gezelter, J. Chem. Phys. 128, 24107 (2008).          
 * [4]  Kuang & Gezelter,  J. Chem. Phys. 133, 164101 (2010).
 * [5]  Vardeman, Stocker & Gezelter, J. Chem. Theory Comput. 7, 834 (2011).
 */
 
#ifndef NONBONDED_NONBONDEDINTERACTION_HPP
#define NONBONDED_NONBONDEDINTERACTION_HPP

#include "types/AtomType.hpp"
#include "math/SquareMatrix3.hpp"

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
    NO_FAMILY = 0,             /**< No family defined */
    VANDERWAALS_FAMILY = 1,    /**< Long-range dispersion and short-range pauli repulsion */
    ELECTROSTATIC_FAMILY = 2,  /**< Coulombic and point-multipole interactions */
    METALLIC_FAMILY = 3,       /**< Transition metal interactions involving electron density */
    HYDROGENBONDING_FAMILY = 4,/**< Short-range directional interactions */
    N_INTERACTION_FAMILIES = 5
  };

  typedef Vector<RealType, N_INTERACTION_FAMILIES> potVec;

  /**
   * The InteractionData struct.  
   *
   * This is used to pass pointers to data to specific non-bonded
   * interactions for force calculations.  Not all of the struct
   * members are utilized by any given interaction.
   */
  struct InteractionData {
    pair<AtomType*, AtomType*> atypes; /**< pair of atom types interacting */
    Vector3d* d;              /**< interatomic vector (already wrapped into box) */
    RealType* rij;            /**< interatomic separation */
    RealType* r2;             /**< square of rij */
    RealType* rcut;           /**< cutoff radius for this interaction */
    bool shiftedPot;          /**< shift the potential up inside the cutoff? */
    bool shiftedForce;        /**< shifted forces smoothly inside the cutoff? */
    RealType* sw;             /**< switching function value at rij */
    int* topoDist;            /**< topological distance between atoms */
    bool excluded;            /**< is this excluded from *direct* pairwise interactions? (some indirect interactions may still apply) */
    RealType* vdwMult;        /**< multiplier for van der Waals interactions */
    RealType* electroMult;    /**< multiplier for electrostatic interactions */
    potVec* pot;              /**< total potential */
    potVec* excludedPot;      /**< potential energy excluded from the overall calculation */
    RealType* vpair;          /**< pair potential */
    bool doParticlePot;       /**< should we bother with the particle pot? */
    RealType* particlePot1;   /**< pointer to particle potential for atom1 */
    RealType* particlePot2;   /**< pointer to particle potential for atom2 */
    Vector3d* f1;             /**< force between the two atoms */
    Mat3x3d* eFrame1;         /**< pointer to electrostatic frame for atom 1 */
    Mat3x3d* eFrame2;         /**< pointer to electrostatic frame for atom 2 */
    RotMat3x3d* A1;           /**< pointer to rotation matrix of first atom */
    RotMat3x3d* A2;           /**< pointer to rotation matrix of second atom */
    Vector3d* t1;             /**< pointer to torque on first atom */
    Vector3d* t2;             /**< pointer to torque on second atom */
    RealType* rho1;           /**< total electron density at first atom */
    RealType* rho2;           /**< total electron density at second atom */
    RealType* frho1;          /**< density functional at first atom */
    RealType* frho2;          /**< density functional at second atom */
    RealType* dfrho1;         /**< derivative of functional for atom 1 */
    RealType* dfrho2;         /**< derivative of functional for atom 2 */
    RealType* flucQ1;         /**< fluctuating charge on atom1 */
    RealType* flucQ2;         /**< fluctuating charge on atom2 */
    RealType* dVdFQ1;         /**< fluctuating charge force on atom1 */
    RealType* dVdFQ2;         /**< fluctuating charge force on atom2 */
    Vector3d* eField1;        /**< pointer to electric field on first atom */
    Vector3d* eField2;        /**< pointer to electric field on second atom */
    RealType* skippedCharge1; /**< charge skipped on atom1 in pairwise interaction loop with atom2 */
    RealType* skippedCharge2; /**< charge skipped on atom2 in pairwise interaction loop with atom1 */
  };
  
  /** 
   * The SelfData struct.
   * 
   * This is used to pass pointers to data for the self-interaction or
   * derived information on a single atom after a pass through all
   * other interactions.  This is used by electrostatic methods that
   * have long-range corrections involving interactions with a medium
   * or a boundary and also by specific metal interactions for
   * electron density functional calculations.  Not all of the struct
   * members are utilized by any given self interaction.
   */
  struct SelfData {
    AtomType* atype;        /**< pointer to AtomType of the atom */
    Mat3x3d* eFrame;        /**< pointer to electrostatic frame for atom */
    RealType* skippedCharge;/**< charge skipped in pairwise interaction loop */
    potVec* pot;            /**< total potential */
    bool doParticlePot;     /**< should we bother with the particle pot? */
    RealType* particlePot;  /**< contribution to potential from this particle */
    Vector3d* t;            /**< pointer to resultant torque on atom */
    RealType* rho;          /**< electron density */
    RealType* frho;         /**< value of density functional for atom */
    RealType* dfrhodrho;    /**< derivative of density functional for atom */
    RealType* flucQ;	    /**< current value of atom's fluctuating charge */
    RealType* dVdFQ;	    /**< fluctuating charge derivative */
  };
  
    
  /**
   * The basic interface for non-bonded interactions.  
   */
  class NonBondedInteraction {
  public:
    NonBondedInteraction() {}
    virtual ~NonBondedInteraction() {}
    virtual void calcForce(InteractionData &idat) = 0;
    virtual InteractionFamily getFamily() = 0;
    virtual RealType getSuggestedCutoffRadius(pair<AtomType*, AtomType*> atypes) = 0;
    virtual string getName() =  0;
  };    

  /**
   * The basic interface for van der Waals interactions.
   */
  class VanDerWaalsInteraction : public NonBondedInteraction {
  public:
    VanDerWaalsInteraction() : NonBondedInteraction() { }
    virtual ~VanDerWaalsInteraction() {}
    virtual InteractionFamily getFamily() {return VANDERWAALS_FAMILY;}
  };    

  /**
   * The basic interface for electrostatic interactions.
   */
  class ElectrostaticInteraction : public NonBondedInteraction {
  public:
    ElectrostaticInteraction() : NonBondedInteraction() { }
    virtual ~ElectrostaticInteraction() {}
    virtual void calcSelfCorrection(SelfData &sdat) = 0;
    virtual InteractionFamily getFamily() {return ELECTROSTATIC_FAMILY;}    
  };    

  /**
   * The basic interface for metallic interactions.
   */
  class MetallicInteraction : public NonBondedInteraction {
  public:
    MetallicInteraction() : NonBondedInteraction() { }
    virtual ~MetallicInteraction() {}
    virtual void calcDensity(InteractionData &idat) = 0;
    virtual void calcFunctional(SelfData &sdat) = 0;
    virtual InteractionFamily getFamily() {return METALLIC_FAMILY;}
  };
          
  /**
   * The basic interface for hydrogen bonding interactions.
   */
  class HydrogenBondingInteraction : public NonBondedInteraction {
  public:
    HydrogenBondingInteraction() : NonBondedInteraction() { }
    virtual ~HydrogenBondingInteraction() {}
    virtual InteractionFamily getFamily() {return HYDROGENBONDING_FAMILY;}
  };    

} //end namespace OpenMD
#endif
