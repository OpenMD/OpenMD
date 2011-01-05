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
 * [4]  Vardeman & Gezelter, in progress (2009).                        
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
    NO_FAMILY,             /**< No family defined */
    VANDERWAALS_FAMILY,    /**< Long-range dispersion and short-range pauli repulsion */
    ELECTROSTATIC_FAMILY,  /**< Coulombic and point-multipole interactions */
    METALLIC_FAMILY,       /**< Transition metal interactions involving electron density */
    HYDROGENBONDING_FAMILY /**< Short-range directional interactions */
  };

  /**
   * The InteractionData struct.  
   *
   * This is used to pass data to specific non-bonded interactions for
   * force calculations.  Not all of the struct members are utilized
   * by any given interaction.
   */
  struct InteractionData {
    AtomType* atype1;     /**< pointer to AtomType of first atom */
    AtomType* atype2;     /**< pointer to AtomType of second atom */
    Vector3d d;           /**< interatomic vector (already wrapped into box) */
    RealType rij;         /**< interatomic separation (precomputed) */
    RealType r2;          /**< square of rij (precomputed) */
    RealType rcut;        /**< cutoff radius for this interaction */
    RealType sw;          /**< switching function value at rij (precomputed) */
    RealType vdwMult;     /**< multiplier for van der Waals interactions */
    RealType electroMult; /**< multiplier for electrostatic interactions */
    RealType pot[4];      /**< total potential */
    RealType vpair[4];    /**< pair potential */
    Vector3d f1;          /**< force between the two atoms */
    Mat3x3d eFrame1;      /**< pointer to electrostatic frame for first atom */
    Mat3x3d eFrame2;      /**< pointer to electrostatic frame for second atom*/
    RotMat3x3d A1;        /**< pointer to rotation matrix of first atom */
    RotMat3x3d A2;        /**< pointer to rotation matrix of second atom */
    Vector3d t1;          /**< pointer to torque on first atom */
    Vector3d t2;          /**< pointer to torque on second atom */
    RealType rho1;        /**< electron density at first atom */
    RealType rho2;        /**< electron density at second atom */
    RealType dfrho1;      /**< derivative of density functional for atom 1 */
    RealType dfrho2;      /**< derivative of density functional for atom 2 */
    RealType fshift1;     /**< indirect potential contribution from atom 1 */
    RealType fshift2;     /**< indirect potential contribution from atom 2 */
  };
  
  /** 
   * The SkipCorrectionData struct.
   *
   * This is used to pass data for corrections due to "skipped" pairs
   * of atoms.  These are atoms that are topologically bonded to each
   * other, but which have indirect effects by way of long range
   * interactions.  In the normal pair interaction loop, these atoms
   * are skipped entirely, so a second pass must be made to compute
   * their indirect interactions on each other.
   */
  struct SkipCorrectionData {
    AtomType* atype1;      /**< pointer to AtomType of first atom */
    AtomType* atype2;      /**< pointer to AtomType of second atom */
    Vector3d d;            /**< interatomic vector (already wrapped into box) */
    RealType rij;          /**< interatomic separation (precomputed) */
    RealType skippedCharge1; /**< charge skipped in normal pairwise interaction loop */
    RealType skippedCharge2; /**< charge skipped in normal pairwise interaction loop */
    RealType sw;           /**< switching function value at rij (precomputed) */
    RealType electroMult;  /**< multiplier for electrostatic interactions */
    RealType pot[4];       /**< total potential */
    RealType vpair[4];     /**< pair potential */
    Vector3d f1;           /**< force correction */
    Mat3x3d eFrame1;       /**< pointer to electrostatic frame for first atom */
    Mat3x3d eFrame2;       /**< pointer to electrostatic frame for second atom*/
    Vector3d t1;           /**< pointer to torque on first atom */
    Vector3d t2;           /**< pointer to torque on second atom */
  };

  /** 
   * The SelfCorrectionData struct.
   * 
   * This is used to pass data for the self-interaction (used by
   * electrostatic methods) that have long-range corrections involving
   * interactions with a medium or a boundary.
   */
  struct SelfCorrectionData {
    AtomType* atype;        /**< pointer to AtomType of the atom */
    Mat3x3d eFrame;        /**< pointer to electrostatic frame for first atom */
    RealType skippedCharge; /**< charge skipped in normal pairwise interaction loop */
    RealType pot[4];       /**< total potential contribution */
    Vector3d t;            /**< pointer to resultant torque on atom */
  };



  /**
   * The DensityData struct.  
   *
   * This is used to pass data to specific metallic interactions for
   * electron density calculations.  
   */
  
  struct DensityData {
    AtomType* atype1;     /**< pointer to AtomType of first atom */
    AtomType* atype2;     /**< pointer to AtomType of second atom */
    RealType rij;         /**< interatomic separation (precomputed) */
    RealType rho_i_at_j;  /**< electron density at second atom due to first */
    RealType rho_j_at_i;  /**< electron density at first atom due to second */
  };
  
  /**
   * The FunctionalData struct.  
   *
   * This is used to pass data to specific metallic interactions for
   * electron density functional calculations.  
   */
  
  struct FunctionalData {
    AtomType* atype;     /**< pointer to AtomType of the atom */
    RealType rho;        /**< electron density (precomputed) */
    RealType frho;       /**< value of density functional for the atom */
    RealType dfrhodrho;  /**< derivative of density functional for the atom */
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
    virtual RealType getSuggestedCutoffRadius(AtomType* at1, AtomType* at2) = 0;
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
    virtual void calcSkipCorrection(SkipCorrectionData &skdat) = 0;
    virtual void calcSelfCorrection(SelfCorrectionData &scdat) = 0;
    virtual InteractionFamily getFamily() {return ELECTROSTATIC_FAMILY;}    
  };    

  /**
   * The basic interface for metallic interactions.
   */
  class MetallicInteraction : public NonBondedInteraction {
  public:
    MetallicInteraction() : NonBondedInteraction() { }
    virtual ~MetallicInteraction() {}
    virtual void calcDensity(DensityData &ddat) = 0;
    virtual void calcFunctional(FunctionalData &fdat) = 0;
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
