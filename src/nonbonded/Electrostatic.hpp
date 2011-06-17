/*
 * Copyright (c) 2009 The University of Notre Dame. All Rights Reserved.
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
 
#ifndef NONBONDED_ELECTROSTATIC_HPP
#define NONBONDED_ELECTROSTATIC_HPP

#include "nonbonded/NonBondedInteraction.hpp"
#include "types/AtomType.hpp"
#include "UseTheForce/ForceField.hpp"
#include "math/SquareMatrix3.hpp"
#include "math/CubicSpline.hpp"
#include "brains/SimInfo.hpp"

namespace OpenMD {

  struct ElectrostaticAtomData {
    bool is_Charge;
    bool is_Dipole;
    bool is_SplitDipole;
    bool is_Quadrupole;
    RealType charge;
    RealType dipole_moment;
    RealType split_dipole_distance;
    Vector3d quadrupole_moments;
  };
  
  enum ElectrostaticSummationMethod{
    esm_HARD,
    esm_SWITCHING_FUNCTION,
    esm_SHIFTED_POTENTIAL,
    esm_SHIFTED_FORCE,
    esm_REACTION_FIELD,
    esm_EWALD_FULL,  /**< Ewald methods aren't supported yet */
    esm_EWALD_PME,   /**< Ewald methods aren't supported yet */
    esm_EWALD_SPME   /**< Ewald methods aren't supported yet */
  };

  enum ElectrostaticScreeningMethod{
    UNDAMPED,
    DAMPED
  };
  
  class Electrostatic : public ElectrostaticInteraction {
    
  public:    
    Electrostatic();
    void setForceField(ForceField *ff) {forceField_ = ff;};
    void setSimInfo(SimInfo* info) {info_ = info;};
    void addType(AtomType* atomType);
    virtual void calcForce(InteractionData &idat);
    virtual void calcSkipCorrection(InteractionData &idat);
    virtual void calcSelfCorrection(SelfData &sdat);
    virtual string getName() {return name_;}
    virtual RealType getSuggestedCutoffRadius(pair<AtomType*, AtomType*> atypes);
    void setCutoffRadius( RealType rCut );
    void setSwitchingRadius( RealType rSwitch );
    void setElectrostaticSummationMethod( ElectrostaticSummationMethod esm );
    void setElectrostaticScreeningMethod( ElectrostaticScreeningMethod sm );
    void setDampingAlpha( RealType alpha );
    void setReactionFieldDielectric( RealType dielectric );

  private:
    void initialize();
    string name_;
    bool initialized_;
    bool haveCutoffRadius_;
    bool haveDampingAlpha_;
    bool haveDielectric_;
    bool haveElectroSpline_;
    std::map<int, AtomType*> ElectrostaticList;
    std::map<AtomType*, ElectrostaticAtomData> ElectrostaticMap;
    SimInfo* info_;
    ForceField* forceField_;
    RealType cutoffRadius_;
    RealType cutoffRadius2_;
    RealType pre11_;
    RealType pre12_;
    RealType pre22_;
    RealType pre14_;
    RealType chargeToC_;
    RealType angstromToM_;
    RealType debyeToCm_;
    int np_;
    ElectrostaticSummationMethod summationMethod_;    
    ElectrostaticScreeningMethod screeningMethod_;
    map<string, ElectrostaticSummationMethod> summationMap_;
    map<string, ElectrostaticScreeningMethod> screeningMap_;
    RealType dampingAlpha_;
    RealType alpha2_;
    RealType alpha4_;
    RealType alpha6_;
    RealType alpha8_;
    RealType dielectric_;
    RealType constEXP_;
    RealType rcuti_;
    RealType rcuti2_;
    RealType rcuti3_;
    RealType rcuti4_;
    RealType alphaPi_;
    RealType invRootPi_;
    RealType rrf_;
    RealType rt_;
    RealType rrfsq_;
    RealType preRF_;
    RealType preRF2_;
    RealType erfcVal_;
    RealType derfcVal_;
    CubicSpline* erfcSpline_;
    RealType c1_;
    RealType c2_;
    RealType c3_;
    RealType c4_;
    RealType c5_;
    RealType c6_;
    RealType c1c_;
    RealType c2c_;
    RealType c3c_;
    RealType c4c_;
    RealType c5c_;
    RealType c6c_;
    RealType one_third_;    
  };
}

                               
#endif
