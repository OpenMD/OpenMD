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
 * [3]  Sun, Lin & Gezelter, J. Chem. Phys. 128, 234107 (2008).          
 * [4]  Kuang & Gezelter,  J. Chem. Phys. 133, 164101 (2010).
 * [5]  Vardeman, Stocker & Gezelter, J. Chem. Theory Comput. 7, 834 (2011).
 */
 
#ifndef NONBONDED_ELECTROSTATIC_HPP
#define NONBONDED_ELECTROSTATIC_HPP

#include "nonbonded/NonBondedInteraction.hpp"
#include "types/AtomType.hpp"
#include "brains/ForceField.hpp"
#include "math/SquareMatrix3.hpp"
#include "math/CubicSpline.hpp"
#include "brains/SimInfo.hpp"

namespace OpenMD {

  struct ElectrostaticAtomData {
    bool is_Charge;
    bool is_Dipole;
    bool is_Quadrupole;
    bool is_Fluctuating;
    RealType fixedCharge;
    RealType hardness;
    RealType electronegativity;
    int slaterN;
    RealType slaterZeta;
    Vector3d dipole;
    Mat3x3d  quadrupole;
  };
      
  enum ElectrostaticSummationMethod{
    esm_HARD,
    esm_SWITCHING_FUNCTION,
    esm_SHIFTED_POTENTIAL,
    esm_SHIFTED_FORCE,
    esm_TAYLOR_SHIFTED,
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
    virtual void calcSelfCorrection(SelfData &sdat);
    virtual string getName() {return name_;}
    virtual RealType getSuggestedCutoffRadius(pair<AtomType*, AtomType*> atypes);
    void setCutoffRadius( RealType rCut );
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
    bool haveElectroSplines_;
    std::map<int, AtomType*> ElectrostaticList;
    std::map<AtomType*, ElectrostaticAtomData> ElectrostaticMap;
    map<pair<AtomType*, AtomType*>, CubicSpline*> Jij;  /** coulomb integral */
    SimInfo* info_;
    ForceField* forceField_;
    RealType cutoffRadius_;
    RealType pre11_;
    RealType pre12_;
    RealType pre22_;
    RealType pre14_;
    RealType pre24_;
    RealType pre44_;
    RealType v01, v11, v21, v22, v31, v32, v41, v42, v43;
    RealType dv01, dv11, dv21, dv22, dv31, dv32, dv41, dv42, dv43;
    RealType v01or, v11or, v21or, v22or, v31or, v32or, v41or, v42or, v43or;
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
    RealType selfMult_;

    CubicSpline* v01s;
    CubicSpline* v11s;
    CubicSpline* v21s;
    CubicSpline* v22s;
    CubicSpline* v31s;
    CubicSpline* v32s;
    CubicSpline* v41s;
    CubicSpline* v42s;
    CubicSpline* v43s;

    /*
    CubicSpline* dv01s;
    CubicSpline* dv11s;
    CubicSpline* dv21s;
    CubicSpline* dv22s;
    CubicSpline* dv31s;
    CubicSpline* dv32s;
    CubicSpline* dv41s;
    CubicSpline* dv42s;
    CubicSpline* dv43s;
    */
  };
}

                               
#endif
