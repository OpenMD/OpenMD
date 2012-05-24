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
 * [2]  Fennell & Gezelter, J. Chem. Phys. 124 234104 (2006).          
 * [3]  Sun, Lin & Gezelter, J. Chem. Phys. 128, 24107 (2008).          
 * [4]  Kuang & Gezelter,  J. Chem. Phys. 133, 164101 (2010).
 * [5]  Vardeman, Stocker & Gezelter, J. Chem. Theory Comput. 7, 834 (2011).
 */

#include <stdio.h>
#include <string.h>

#include <cmath>
#include "nonbonded/Electrostatic.hpp"
#include "utils/simError.h"
#include "types/NonBondedInteractionType.hpp"
#include "types/FixedChargeAdapter.hpp"
#include "types/FluctuatingChargeAdapter.hpp"
#include "types/MultipoleAdapter.hpp"
#include "io/Globals.hpp"
#include "nonbonded/SlaterIntegrals.hpp"
#include "utils/PhysicalConstants.hpp"


namespace OpenMD {
  
  Electrostatic::Electrostatic(): name_("Electrostatic"), initialized_(false),
                                  forceField_(NULL), info_(NULL), 
                                  haveCutoffRadius_(false),
                                  haveDampingAlpha_(false), 
                                  haveDielectric_(false),
                                  haveElectroSpline_(false)
  {}
  
  void Electrostatic::initialize() { 
    
    Globals* simParams_ = info_->getSimParams();

    summationMap_["HARD"]               = esm_HARD;
    summationMap_["NONE"]               = esm_HARD;
    summationMap_["SWITCHING_FUNCTION"] = esm_SWITCHING_FUNCTION;
    summationMap_["SHIFTED_POTENTIAL"]  = esm_SHIFTED_POTENTIAL; 
    summationMap_["SHIFTED_FORCE"]      = esm_SHIFTED_FORCE;     
    summationMap_["REACTION_FIELD"]     = esm_REACTION_FIELD;    
    summationMap_["EWALD_FULL"]         = esm_EWALD_FULL;        
    summationMap_["EWALD_PME"]          = esm_EWALD_PME;         
    summationMap_["EWALD_SPME"]         = esm_EWALD_SPME;        
    screeningMap_["DAMPED"]             = DAMPED;
    screeningMap_["UNDAMPED"]           = UNDAMPED;

    // these prefactors convert the multipole interactions into kcal / mol
    // all were computed assuming distances are measured in angstroms
    // Charge-Charge, assuming charges are measured in electrons
    pre11_ = 332.0637778;
    // Charge-Dipole, assuming charges are measured in electrons, and
    // dipoles are measured in debyes
    pre12_ = 69.13373;
    // Dipole-Dipole, assuming dipoles are measured in debyes
    pre22_ = 14.39325;
    // Charge-Quadrupole, assuming charges are measured in electrons, and
    // quadrupoles are measured in 10^-26 esu cm^2
    // This unit is also known affectionately as an esu centi-barn.
    pre14_ = 69.13373;
    
    // conversions for the simulation box dipole moment
    chargeToC_ = 1.60217733e-19;
    angstromToM_ = 1.0e-10;
    debyeToCm_ = 3.33564095198e-30;
    
    // number of points for electrostatic splines
    np_ = 100;
    
    // variables to handle different summation methods for long-range 
    // electrostatics:
    summationMethod_ = esm_HARD;    
    screeningMethod_ = UNDAMPED;
    dielectric_ = 1.0;
    one_third_ = 1.0 / 3.0;
  
    // check the summation method:
    if (simParams_->haveElectrostaticSummationMethod()) {
      string myMethod = simParams_->getElectrostaticSummationMethod();
      toUpper(myMethod);
      map<string, ElectrostaticSummationMethod>::iterator i;
      i = summationMap_.find(myMethod);
      if ( i != summationMap_.end() ) {
        summationMethod_ = (*i).second;
      } else {
        // throw error
        sprintf( painCave.errMsg,
                 "Electrostatic::initialize: Unknown electrostaticSummationMethod.\n"
                 "\t(Input file specified %s .)\n"
                 "\telectrostaticSummationMethod must be one of: \"hard\",\n"
                 "\t\"shifted_potential\", \"shifted_force\", or \n"
                 "\t\"reaction_field\".\n", myMethod.c_str() );
        painCave.isFatal = 1;
        simError();
      }
    } else {
      // set ElectrostaticSummationMethod to the cutoffMethod:
      if (simParams_->haveCutoffMethod()){
        string myMethod = simParams_->getCutoffMethod();
        toUpper(myMethod);
        map<string, ElectrostaticSummationMethod>::iterator i;
        i = summationMap_.find(myMethod);
        if ( i != summationMap_.end() ) {
          summationMethod_ = (*i).second;
        }
      }
    }
    
    if (summationMethod_ == esm_REACTION_FIELD) {        
      if (!simParams_->haveDielectric()) {
        // throw warning
        sprintf( painCave.errMsg,
                 "SimInfo warning: dielectric was not specified in the input file\n\tfor the reaction field correction method.\n"
                 "\tA default value of %f will be used for the dielectric.\n", dielectric_);
        painCave.isFatal = 0;
        painCave.severity = OPENMD_INFO;
        simError();
      } else {
        dielectric_ = simParams_->getDielectric();       
      }
      haveDielectric_ = true;
    }
    
    if (simParams_->haveElectrostaticScreeningMethod()) {
      string myScreen = simParams_->getElectrostaticScreeningMethod();
      toUpper(myScreen);
      map<string, ElectrostaticScreeningMethod>::iterator i;
      i = screeningMap_.find(myScreen);
      if ( i != screeningMap_.end()) {
        screeningMethod_ = (*i).second;
      } else {
        sprintf( painCave.errMsg,
                 "SimInfo error: Unknown electrostaticScreeningMethod.\n"
                 "\t(Input file specified %s .)\n"
                 "\telectrostaticScreeningMethod must be one of: \"undamped\"\n"
                 "or \"damped\".\n", myScreen.c_str() );
        painCave.isFatal = 1;
        simError();
      }
    }

    // check to make sure a cutoff value has been set:
    if (!haveCutoffRadius_) {
      sprintf( painCave.errMsg, "Electrostatic::initialize has no Default "
               "Cutoff value!\n");
      painCave.severity = OPENMD_ERROR;
      painCave.isFatal = 1;
      simError();
    }
           
    if (screeningMethod_ == DAMPED) {      
      if (!simParams_->haveDampingAlpha()) {
        // first set a cutoff dependent alpha value
        // we assume alpha depends linearly with rcut from 0 to 20.5 ang
        dampingAlpha_ = 0.425 - cutoffRadius_* 0.02;
        if (dampingAlpha_ < 0.0) dampingAlpha_ = 0.0;
        
        // throw warning
        sprintf( painCave.errMsg,
                 "Electrostatic::initialize: dampingAlpha was not specified in the input file.\n"
                 "\tA default value of %f (1/ang) will be used for the cutoff of\n\t%f (ang).\n", 
                 dampingAlpha_, cutoffRadius_);
        painCave.severity = OPENMD_INFO;
        painCave.isFatal = 0;
        simError();
      } else {
        dampingAlpha_ = simParams_->getDampingAlpha();
      }
      haveDampingAlpha_ = true;
    }

    // find all of the Electrostatic atom Types:
    ForceField::AtomTypeContainer* atomTypes = forceField_->getAtomTypes();
    ForceField::AtomTypeContainer::MapTypeIterator i;
    AtomType* at;
    
    for (at = atomTypes->beginType(i); at != NULL; 
         at = atomTypes->nextType(i)) {
      
      if (at->isElectrostatic())
        addType(at);
    }
    

    cutoffRadius2_ = cutoffRadius_ * cutoffRadius_;
    rcuti_ = 1.0 / cutoffRadius_;
    rcuti2_ = rcuti_ * rcuti_;
    rcuti3_ = rcuti2_ * rcuti_;
    rcuti4_ = rcuti2_ * rcuti2_;

    if (screeningMethod_ == DAMPED) {
      
      alpha2_ = dampingAlpha_ * dampingAlpha_;
      alpha4_ = alpha2_ * alpha2_;
      alpha6_ = alpha4_ * alpha2_;
      alpha8_ = alpha4_ * alpha4_;
      
      constEXP_ = exp(-alpha2_ * cutoffRadius2_);
      invRootPi_ = 0.56418958354775628695;
      alphaPi_ = 2.0 * dampingAlpha_ * invRootPi_;

      c1c_ = erfc(dampingAlpha_ * cutoffRadius_) * rcuti_;
      c2c_ = alphaPi_ * constEXP_ * rcuti_ + c1c_ * rcuti_;
      c3c_ = 2.0 * alphaPi_ * alpha2_ + 3.0 * c2c_ * rcuti_;
      c4c_ = 4.0 * alphaPi_ * alpha4_ + 5.0 * c3c_ * rcuti2_;
      c5c_ = 8.0 * alphaPi_ * alpha6_ + 7.0 * c4c_ * rcuti2_;
      c6c_ = 16.0 * alphaPi_ * alpha8_ + 9.0 * c5c_ * rcuti2_;
    } else {
      c1c_ = rcuti_;
      c2c_ = c1c_ * rcuti_;
      c3c_ = 3.0 * c2c_ * rcuti_;
      c4c_ = 5.0 * c3c_ * rcuti2_;
      c5c_ = 7.0 * c4c_ * rcuti2_;
      c6c_ = 9.0 * c5c_ * rcuti2_;
    }
   
    if (summationMethod_ == esm_REACTION_FIELD) {
      preRF_ = (dielectric_ - 1.0) / 
        ((2.0 * dielectric_ + 1.0) * cutoffRadius2_ * cutoffRadius_);
      preRF2_ = 2.0 * preRF_;
    }
    
    // Add a 2 angstrom safety window to deal with cutoffGroups that
    // have charged atoms longer than the cutoffRadius away from each
    // other.  Splining may not be the best choice here.  Direct calls
    // to erfc might be preferrable.

    RealType dx = (cutoffRadius_ + 2.0) / RealType(np_ - 1);
    RealType rval;
    vector<RealType> rvals;
    vector<RealType> yvals;
    for (int i = 0; i < np_; i++) {
      rval = RealType(i) * dx;
      rvals.push_back(rval);
      yvals.push_back(erfc(dampingAlpha_ * rval));
    }
    erfcSpline_ = new CubicSpline();
    erfcSpline_->addPoints(rvals, yvals);
    haveElectroSpline_ = true;

    initialized_ = true;
  }
      
  void Electrostatic::addType(AtomType* atomType){

    ElectrostaticAtomData electrostaticAtomData;
    electrostaticAtomData.is_Charge = false;
    electrostaticAtomData.is_Dipole = false;
    electrostaticAtomData.is_SplitDipole = false;
    electrostaticAtomData.is_Quadrupole = false;

    FixedChargeAdapter fca = FixedChargeAdapter(atomType);

    if (fca.isFixedCharge()) {
      electrostaticAtomData.is_Charge = true;
      electrostaticAtomData.fixedCharge = fca.getCharge();
    }

    MultipoleAdapter ma = MultipoleAdapter(atomType);
    if (ma.isMultipole()) {
      if (ma.isDipole()) {
        electrostaticAtomData.is_Dipole = true;
        electrostaticAtomData.dipole_moment = ma.getDipoleMoment();
      }
      if (ma.isSplitDipole()) {
        electrostaticAtomData.is_SplitDipole = true;
        electrostaticAtomData.split_dipole_distance = ma.getSplitDipoleDistance();
      }
      if (ma.isQuadrupole()) {
        // Quadrupoles in OpenMD are set as the diagonal elements
        // of the diagonalized traceless quadrupole moment tensor.
        // The column vectors of the unitary matrix that diagonalizes 
        // the quadrupole moment tensor become the eFrame (or the
        // electrostatic version of the body-fixed frame.
        electrostaticAtomData.is_Quadrupole = true;
        electrostaticAtomData.quadrupole_moments = ma.getQuadrupoleMoments();
      }
    }
    
    FluctuatingChargeAdapter fqa = FluctuatingChargeAdapter(atomType);

    if (fqa.isFluctuatingCharge()) {
      electrostaticAtomData.is_Fluctuating = true;
      electrostaticAtomData.electronegativity = fqa.getElectronegativity();
      electrostaticAtomData.hardness = fqa.getHardness();
      electrostaticAtomData.slaterN = fqa.getSlaterN();
      electrostaticAtomData.slaterZeta = fqa.getSlaterZeta();
    }

    pair<map<int,AtomType*>::iterator,bool> ret;    
    ret = ElectrostaticList.insert( pair<int,AtomType*>(atomType->getIdent(),
                                                        atomType) );
    if (ret.second == false) {
      sprintf( painCave.errMsg,
               "Electrostatic already had a previous entry with ident %d\n",
               atomType->getIdent() );
      painCave.severity = OPENMD_INFO;
      painCave.isFatal = 0;
      simError();         
    }
    
    ElectrostaticMap[atomType] = electrostaticAtomData;   

    // Now, iterate over all known types and add to the mixing map:
    
    map<AtomType*, ElectrostaticAtomData>::iterator it;
    for( it = ElectrostaticMap.begin(); it != ElectrostaticMap.end(); ++it) {
      AtomType* atype2 = (*it).first;
      ElectrostaticAtomData eaData2 = (*it).second;
      if (eaData2.is_Fluctuating && electrostaticAtomData.is_Fluctuating) {
        
        RealType a = electrostaticAtomData.slaterZeta;
        RealType b = eaData2.slaterZeta;
        int m = electrostaticAtomData.slaterN;
        int n = eaData2.slaterN;

        // Create the spline of the coulombic integral for s-type
        // Slater orbitals.  Add a 2 angstrom safety window to deal
        // with cutoffGroups that have charged atoms longer than the
        // cutoffRadius away from each other.

        RealType rval;
        RealType dr = (cutoffRadius_ + 2.0) / RealType(np_ - 1);
        vector<RealType> rvals;
        vector<RealType> J1vals;
        vector<RealType> J2vals;
        for (int i = 0; i < np_; i++) {
          rval = RealType(i) * dr;
          rvals.push_back(rval);
          J1vals.push_back( sSTOCoulInt( a, b, m, n, rval * PhysicalConstants::angstromsToBohr ) );
          J2vals.push_back( sSTOCoulInt( b, a, n, m, rval * PhysicalConstants::angstromsToBohr ) );
        }

        CubicSpline* J1 = new CubicSpline();
        J1->addPoints(rvals, J1vals);
        CubicSpline* J2 = new CubicSpline();
        J2->addPoints(rvals, J2vals);
        
        pair<AtomType*, AtomType*> key1, key2;
        key1 = make_pair(atomType, atype2);
        key2 = make_pair(atype2, atomType);
        
        Jij[key1] = J1;
        Jij[key2] = J2;
      }
    }
 
    return;
  }
  
  void Electrostatic::setCutoffRadius( RealType rCut ) {
    cutoffRadius_ = rCut;
    rrf_ = cutoffRadius_;
    haveCutoffRadius_ = true;
  }

  void Electrostatic::setSwitchingRadius( RealType rSwitch ) {
    rt_ = rSwitch;
  }
  void Electrostatic::setElectrostaticSummationMethod( ElectrostaticSummationMethod esm ) {
    summationMethod_ = esm;
  }
  void Electrostatic::setElectrostaticScreeningMethod( ElectrostaticScreeningMethod sm ) {
    screeningMethod_ = sm;
  }
  void Electrostatic::setDampingAlpha( RealType alpha ) {
    dampingAlpha_ = alpha;
    haveDampingAlpha_ = true;
  }
  void Electrostatic::setReactionFieldDielectric( RealType dielectric ){
    dielectric_ = dielectric;
    haveDielectric_ = true;
  }

  void Electrostatic::calcForce(InteractionData &idat) {

    // utility variables.  Should clean these up and use the Vector3d and 
    // Mat3x3d to replace as many as we can in future versions:

    RealType q_i, q_j, mu_i, mu_j, d_i, d_j;
    RealType qxx_i, qyy_i, qzz_i;
    RealType qxx_j, qyy_j, qzz_j;
    RealType cx_i, cy_i, cz_i;
    RealType cx_j, cy_j, cz_j;
    RealType cx2, cy2, cz2;
    RealType ct_i, ct_j, ct_ij, a1;
    RealType riji, ri, ri2, ri3, ri4;
    RealType pref, vterm, epot, dudr;
    RealType vpair(0.0);
    RealType scale, sc2;
    RealType pot_term, preVal, rfVal;
    RealType c2ri, c3ri, c4rij, cti3, ctj3, ctidotj;
    RealType preSw, preSwSc;
    RealType c1, c2, c3, c4;
    RealType erfcVal(1.0), derfcVal(0.0);
    RealType BigR;
    RealType two(2.0), three(3.0);

    Vector3d Q_i, Q_j;
    Vector3d ux_i, uy_i, uz_i;
    Vector3d ux_j, uy_j, uz_j;
    Vector3d dudux_i, duduy_i, duduz_i;
    Vector3d dudux_j, duduy_j, duduz_j;
    Vector3d rhatdot2, rhatc4;
    Vector3d dVdr;

    // variables for indirect (reaction field) interactions for excluded pairs:
    RealType indirect_Pot(0.0);
    RealType indirect_vpair(0.0);
    Vector3d indirect_dVdr(V3Zero);
    Vector3d indirect_duduz_i(V3Zero), indirect_duduz_j(V3Zero);

    pair<RealType, RealType> res;
    
    if (!initialized_) initialize();
    
    ElectrostaticAtomData data1 = ElectrostaticMap[idat.atypes.first];
    ElectrostaticAtomData data2 = ElectrostaticMap[idat.atypes.second];
    
    // some variables we'll need independent of electrostatic type:

    riji = 1.0 /  *(idat.rij) ;
    Vector3d rhat =  *(idat.d)   * riji;

    // logicals

    bool i_is_Charge = data1.is_Charge;
    bool i_is_Dipole = data1.is_Dipole;
    bool i_is_SplitDipole = data1.is_SplitDipole;
    bool i_is_Quadrupole = data1.is_Quadrupole;

    bool j_is_Charge = data2.is_Charge;
    bool j_is_Dipole = data2.is_Dipole;
    bool j_is_SplitDipole = data2.is_SplitDipole;
    bool j_is_Quadrupole = data2.is_Quadrupole;
    
    if (i_is_Charge) {
      q_i = data1.fixedCharge;
      if (idat.excluded) {
        *(idat.skippedCharge2) += q_i;
      }
    }

    if (i_is_Dipole) {
      mu_i = data1.dipole_moment;
      uz_i = idat.eFrame1->getColumn(2);
      
      ct_i = dot(uz_i, rhat);

      if (i_is_SplitDipole) 
        d_i = data1.split_dipole_distance;
      
      duduz_i = V3Zero;
    }
    
    if (i_is_Quadrupole) {
      Q_i = data1.quadrupole_moments;
      qxx_i = Q_i.x();
      qyy_i = Q_i.y();
      qzz_i = Q_i.z();
      
      ux_i = idat.eFrame1->getColumn(0);
      uy_i = idat.eFrame1->getColumn(1);
      uz_i = idat.eFrame1->getColumn(2);

      cx_i = dot(ux_i, rhat);
      cy_i = dot(uy_i, rhat);
      cz_i = dot(uz_i, rhat); 

      dudux_i = V3Zero;
      duduy_i = V3Zero;
      duduz_i = V3Zero;
    }

    if (j_is_Charge) {
      q_j = data2.fixedCharge;
      if (idat.excluded) {
        *(idat.skippedCharge1) += q_j;
      }
    }


    if (j_is_Dipole) {
      mu_j = data2.dipole_moment;
      uz_j = idat.eFrame2->getColumn(2);
      
      ct_j = dot(uz_j, rhat);

      if (j_is_SplitDipole) 
        d_j = data2.split_dipole_distance;
      
      duduz_j = V3Zero;
    }
    
    if (j_is_Quadrupole) {
      Q_j = data2.quadrupole_moments;
      qxx_j = Q_j.x();
      qyy_j = Q_j.y();
      qzz_j = Q_j.z();
      
      ux_j = idat.eFrame2->getColumn(0);
      uy_j = idat.eFrame2->getColumn(1);
      uz_j = idat.eFrame2->getColumn(2);

      cx_j = dot(ux_j, rhat);
      cy_j = dot(uy_j, rhat);
      cz_j = dot(uz_j, rhat);

      dudux_j = V3Zero;
      duduy_j = V3Zero;
      duduz_j = V3Zero;
    }
    
    epot = 0.0;
    dVdr = V3Zero;
    
    if (i_is_Charge) {
      
      if (j_is_Charge) {
        if (screeningMethod_ == DAMPED) {
          // assemble the damping variables
          //res = erfcSpline_->getValueAndDerivativeAt( *(idat.rij) );
          //erfcVal = res.first;
          //derfcVal = res.second;

          erfcVal = erfc(dampingAlpha_ * *(idat.rij));
          derfcVal = - alphaPi_ * exp(-alpha2_ * *(idat.r2));

          c1 = erfcVal * riji;
          c2 = (-derfcVal + c1) * riji;
        } else {
          c1 = riji;
          c2 = c1 * riji;
        }

        preVal =  *(idat.electroMult) * pre11_ * q_i * q_j;
        
        if (summationMethod_ == esm_SHIFTED_POTENTIAL) {
          vterm = preVal * (c1 - c1c_);
          dudr  = - *(idat.sw)  * preVal * c2;

        } else if (summationMethod_ == esm_SHIFTED_FORCE)  {
          vterm = preVal * ( c1 - c1c_ + c2c_*( *(idat.rij)  - cutoffRadius_) );
          dudr  =  *(idat.sw)  * preVal * (c2c_ - c2);

        } else if (summationMethod_ == esm_REACTION_FIELD) {
          rfVal = preRF_ *  *(idat.rij)  *  *(idat.rij);

          vterm = preVal * ( riji + rfVal );             
          dudr  =  *(idat.sw)  * preVal * ( 2.0 * rfVal - riji ) * riji;
          
          // if this is an excluded pair, there are still indirect
          // interactions via the reaction field we must worry about:

          if (idat.excluded) {
            indirect_vpair += preVal * rfVal;
            indirect_Pot += *(idat.sw) * preVal * rfVal;
            indirect_dVdr += *(idat.sw)  * preVal * two * rfVal  * riji * rhat;
          }
          
        } else {

          vterm = preVal * riji * erfcVal;           
          dudr  = -  *(idat.sw)  * preVal * c2;

        }

        vpair += vterm;
        epot +=  *(idat.sw)  * vterm;
        dVdr += dudr * rhat;                 
      }

      if (j_is_Dipole) {
        // pref is used by all the possible methods
        pref =  *(idat.electroMult) * pre12_ * q_i * mu_j;
        preSw =  *(idat.sw)  * pref;

        if (summationMethod_ == esm_REACTION_FIELD) {
          ri2 = riji * riji;
          ri3 = ri2 * riji;
    
          vterm = - pref * ct_j * ( ri2 - preRF2_ *  *(idat.rij)  );
          vpair += vterm;
          epot +=  *(idat.sw)  * vterm;

          dVdr +=  -preSw * (ri3 * (uz_j - three * ct_j * rhat) - preRF2_*uz_j);
          duduz_j += -preSw * rhat * (ri2 - preRF2_ *  *(idat.rij) );   

          // Even if we excluded this pair from direct interactions,
          // we still have the reaction-field-mediated charge-dipole
          // interaction:

          if (idat.excluded) {
            indirect_vpair += pref * ct_j * preRF2_ * *(idat.rij);
            indirect_Pot += preSw * ct_j * preRF2_ * *(idat.rij);
            indirect_dVdr += preSw * preRF2_ * uz_j;
            indirect_duduz_j += preSw * rhat * preRF2_ *  *(idat.rij); 
          }
                       
        } else {
          // determine the inverse r used if we have split dipoles
          if (j_is_SplitDipole) {
            BigR = sqrt( *(idat.r2) + 0.25 * d_j * d_j);
            ri = 1.0 / BigR;
            scale =  *(idat.rij)  * ri;
          } else {
            ri = riji;
            scale = 1.0;
          }
          
          sc2 = scale * scale;

          if (screeningMethod_ == DAMPED) {
            // assemble the damping variables
            //res = erfcSpline_->getValueAndDerivativeAt( *(idat.rij) );
            //erfcVal = res.first;
            //derfcVal = res.second;
            erfcVal = erfc(dampingAlpha_ * *(idat.rij));
            derfcVal = - alphaPi_ * exp(-alpha2_ * *(idat.r2));
            c1 = erfcVal * ri;
            c2 = (-derfcVal + c1) * ri;
            c3 = -2.0 * derfcVal * alpha2_ + 3.0 * c2 * ri;
          } else {
            c1 = ri;
            c2 = c1 * ri;
            c3 = 3.0 * c2 * ri;
          }
             
          c2ri = c2 * ri;

          // calculate the potential
          pot_term =  scale * c2;
          vterm = -pref * ct_j * pot_term;
          vpair += vterm;
          epot +=  *(idat.sw)  * vterm;
             
          // calculate derivatives for forces and torques

          dVdr += -preSw * (uz_j * c2ri - ct_j * rhat * sc2 * c3);
          duduz_j += -preSw * pot_term * rhat;

        }
      }

      if (j_is_Quadrupole) {
        // first precalculate some necessary variables
        cx2 = cx_j * cx_j;
        cy2 = cy_j * cy_j;
        cz2 = cz_j * cz_j;
        pref =   *(idat.electroMult) * pre14_ * q_i * one_third_;
          
        if (screeningMethod_ == DAMPED) {
          // assemble the damping variables
          //res = erfcSpline_->getValueAndDerivativeAt( *(idat.rij) );
          //erfcVal = res.first;
          //derfcVal = res.second;
          erfcVal = erfc(dampingAlpha_ * *(idat.rij));
          derfcVal = - alphaPi_ * exp(-alpha2_ * *(idat.r2));
          c1 = erfcVal * riji;
          c2 = (-derfcVal + c1) * riji;
          c3 = -2.0 * derfcVal * alpha2_ + 3.0 * c2 * riji;
          c4 = -4.0 * derfcVal * alpha4_ + 5.0 * c3 * riji * riji;
        } else {
          c1 = riji;
          c2 = c1 * riji;
          c3 = 3.0 * c2 * riji;
          c4 = 5.0 * c3 * riji * riji;
        }

        // precompute variables for convenience
        preSw =  *(idat.sw)  * pref;
        c2ri = c2 * riji;
        c3ri = c3 * riji;
        c4rij = c4 *  *(idat.rij) ;
        rhatdot2 = two * rhat * c3;
        rhatc4 = rhat * c4rij;

        // calculate the potential
        pot_term = ( qxx_j * (cx2*c3 - c2ri) + 
                     qyy_j * (cy2*c3 - c2ri) + 
                     qzz_j * (cz2*c3 - c2ri) );
        vterm = pref * pot_term;
        vpair += vterm;
        epot +=  *(idat.sw)  * vterm;
                
        // calculate derivatives for the forces and torques

        dVdr += -preSw * ( qxx_j* (cx2*rhatc4 - (two*cx_j*ux_j + rhat)*c3ri) +
                           qyy_j* (cy2*rhatc4 - (two*cy_j*uy_j + rhat)*c3ri) +
                           qzz_j* (cz2*rhatc4 - (two*cz_j*uz_j + rhat)*c3ri));
                           
        dudux_j += preSw * qxx_j * cx_j * rhatdot2;
        duduy_j += preSw * qyy_j * cy_j * rhatdot2;
        duduz_j += preSw * qzz_j * cz_j * rhatdot2;
      }
    }
    
    if (i_is_Dipole) {

      if (j_is_Charge) {
        // variables used by all the methods
        pref =  *(idat.electroMult) * pre12_ * q_j * mu_i;
        preSw =  *(idat.sw)  * pref;

        if (summationMethod_ == esm_REACTION_FIELD) {

          ri2 = riji * riji;
          ri3 = ri2 * riji;

          vterm = pref * ct_i * ( ri2 - preRF2_ *  *(idat.rij)  );
          vpair += vterm;
          epot +=  *(idat.sw)  * vterm;
          
          dVdr += preSw * (ri3 * (uz_i - three * ct_i * rhat) - preRF2_ * uz_i);
          
          duduz_i += preSw * rhat * (ri2 - preRF2_ *  *(idat.rij) );

          // Even if we excluded this pair from direct interactions,
          // we still have the reaction-field-mediated charge-dipole
          // interaction:

          if (idat.excluded) {
            indirect_vpair += -pref * ct_i * preRF2_ * *(idat.rij);
            indirect_Pot += -preSw * ct_i * preRF2_ * *(idat.rij);
            indirect_dVdr += -preSw * preRF2_ * uz_i;
            indirect_duduz_i += -preSw * rhat * preRF2_ *  *(idat.rij); 
          }
             
        } else {
          
          // determine inverse r if we are using split dipoles
          if (i_is_SplitDipole) {
            BigR = sqrt( *(idat.r2) + 0.25 * d_i * d_i);
            ri = 1.0 / BigR;
            scale =  *(idat.rij)  * ri;
          } else {
            ri = riji;
            scale = 1.0;
          }
          
          sc2 = scale * scale;
            
          if (screeningMethod_ == DAMPED) {
            // assemble the damping variables
            //res = erfcSpline_->getValueAndDerivativeAt( *(idat.rij) );
            //erfcVal = res.first;
            //derfcVal = res.second;
            erfcVal = erfc(dampingAlpha_ * *(idat.rij));
            derfcVal = - alphaPi_ * exp(-alpha2_ * *(idat.r2));
            c1 = erfcVal * ri;
            c2 = (-derfcVal + c1) * ri;
            c3 = -2.0 * derfcVal * alpha2_ + 3.0 * c2 * ri;
          } else {
            c1 = ri;
            c2 = c1 * ri;
            c3 = 3.0 * c2 * ri;
          }
          
          c2ri = c2 * ri;
               
          // calculate the potential
          pot_term = c2 * scale;
          vterm = pref * ct_i * pot_term;
          vpair += vterm;
          epot +=  *(idat.sw)  * vterm;

          // calculate derivatives for the forces and torques
          dVdr += preSw * (uz_i * c2ri - ct_i * rhat * sc2 * c3);
          duduz_i += preSw * pot_term * rhat;
        }
      }

      if (j_is_Dipole) {
        // variables used by all methods
        ct_ij = dot(uz_i, uz_j);

        pref =  *(idat.electroMult) * pre22_ * mu_i * mu_j;
        preSw =  *(idat.sw)  * pref;

        if (summationMethod_ == esm_REACTION_FIELD) {
          ri2 = riji * riji;
          ri3 = ri2 * riji;
          ri4 = ri2 * ri2;

          vterm = pref * ( ri3 * (ct_ij - 3.0 * ct_i * ct_j) -
                           preRF2_ * ct_ij );
          vpair += vterm;
          epot +=  *(idat.sw)  * vterm;
             
          a1 = 5.0 * ct_i * ct_j - ct_ij;
             
          dVdr += preSw * three * ri4 * (a1 * rhat - ct_i * uz_j - ct_j * uz_i);

          duduz_i += preSw * (ri3 * (uz_j - three * ct_j * rhat) - preRF2_*uz_j);
          duduz_j += preSw * (ri3 * (uz_i - three * ct_i * rhat) - preRF2_*uz_i);

          if (idat.excluded) {
            indirect_vpair +=  - pref * preRF2_ * ct_ij;
            indirect_Pot +=    - preSw * preRF2_ * ct_ij;
            indirect_duduz_i += -preSw * preRF2_ * uz_j;
            indirect_duduz_j += -preSw * preRF2_ * uz_i;
          }

        } else {
          
          if (i_is_SplitDipole) {
            if (j_is_SplitDipole) {
              BigR = sqrt( *(idat.r2) + 0.25 * d_i * d_i + 0.25 * d_j * d_j);
            } else {
              BigR = sqrt( *(idat.r2) + 0.25 * d_i * d_i);
            }
            ri = 1.0 / BigR;
            scale =  *(idat.rij)  * ri;
          } else {
            if (j_is_SplitDipole) {
              BigR = sqrt( *(idat.r2) + 0.25 * d_j * d_j);
              ri = 1.0 / BigR;
              scale =  *(idat.rij)  * ri;
            } else {
              ri = riji;
              scale = 1.0;
            }
          }
          if (screeningMethod_ == DAMPED) {
            // assemble damping variables
            //res = erfcSpline_->getValueAndDerivativeAt( *(idat.rij) );
            //erfcVal = res.first;
            //derfcVal = res.second;
            erfcVal = erfc(dampingAlpha_ * *(idat.rij));
            derfcVal = - alphaPi_ * exp(-alpha2_ * *(idat.r2));
            c1 = erfcVal * ri;
            c2 = (-derfcVal + c1) * ri;
            c3 = -2.0 * derfcVal * alpha2_ + 3.0 * c2 * ri;
            c4 = -4.0 * derfcVal * alpha4_ + 5.0 * c3 * ri * ri;
          } else {
            c1 = ri;
            c2 = c1 * ri;
            c3 = 3.0 * c2 * ri;
            c4 = 5.0 * c3 * ri * ri;
          }

          // precompute variables for convenience
          sc2 = scale * scale;
          cti3 = ct_i * sc2 * c3;
          ctj3 = ct_j * sc2 * c3;
          ctidotj = ct_i * ct_j * sc2;
          preSwSc = preSw * scale;
          c2ri = c2 * ri;
          c3ri = c3 * ri;
          c4rij = c4 *  *(idat.rij) ;

          // calculate the potential 
          pot_term = (ct_ij * c2ri - ctidotj * c3);
          vterm = pref * pot_term;
          vpair += vterm;
          epot +=  *(idat.sw)  * vterm;

          // calculate derivatives for the forces and torques
          dVdr += preSwSc * ( ctidotj * rhat * c4rij  - 
                              (ct_i*uz_j + ct_j*uz_i + ct_ij*rhat) * c3ri);
          
          duduz_i += preSw * (uz_j * c2ri - ctj3 * rhat);
          duduz_j += preSw * (uz_i * c2ri - cti3 * rhat);
        }
      }
    }

    if (i_is_Quadrupole) {
      if (j_is_Charge) {
        // precompute some necessary variables
        cx2 = cx_i * cx_i;
        cy2 = cy_i * cy_i;
        cz2 = cz_i * cz_i;

        pref =  *(idat.electroMult) * pre14_ * q_j * one_third_;

        if (screeningMethod_ == DAMPED) {
          // assemble the damping variables
          //res = erfcSpline_->getValueAndDerivativeAt( *(idat.rij) );
          //erfcVal = res.first;
          //derfcVal = res.second;
          erfcVal = erfc(dampingAlpha_ * *(idat.rij));
          derfcVal = - alphaPi_ * exp(-alpha2_ * *(idat.r2));
          c1 = erfcVal * riji;
          c2 = (-derfcVal + c1) * riji;
          c3 = -2.0 * derfcVal * alpha2_ + 3.0 * c2 * riji;
          c4 = -4.0 * derfcVal * alpha4_ + 5.0 * c3 * riji * riji;
        } else {
          c1 = riji;
          c2 = c1 * riji;
          c3 = 3.0 * c2 * riji;
          c4 = 5.0 * c3 * riji * riji;
        }
          
        // precompute some variables for convenience
        preSw =  *(idat.sw)  * pref;
        c2ri = c2 * riji;
        c3ri = c3 * riji;
        c4rij = c4 *  *(idat.rij) ;
        rhatdot2 = two * rhat * c3;
        rhatc4 = rhat * c4rij;

        // calculate the potential
        pot_term = ( qxx_i * (cx2 * c3 - c2ri) + 
                     qyy_i * (cy2 * c3 - c2ri) + 
                     qzz_i * (cz2 * c3 - c2ri) );
        
        vterm = pref * pot_term;
        vpair += vterm;
        epot +=  *(idat.sw)  * vterm;

        // calculate the derivatives for the forces and torques

        dVdr += -preSw * (qxx_i* (cx2*rhatc4 - (two*cx_i*ux_i + rhat)*c3ri) +
                          qyy_i* (cy2*rhatc4 - (two*cy_i*uy_i + rhat)*c3ri) +
                          qzz_i* (cz2*rhatc4 - (two*cz_i*uz_i + rhat)*c3ri));

        dudux_i += preSw * qxx_i * cx_i *  rhatdot2;
        duduy_i += preSw * qyy_i * cy_i *  rhatdot2;
        duduz_i += preSw * qzz_i * cz_i *  rhatdot2;
      }
    }


    if (!idat.excluded) {
      *(idat.vpair) += vpair;
      (*(idat.pot))[ELECTROSTATIC_FAMILY] += epot;
      *(idat.f1) += dVdr;
      
      if (i_is_Dipole || i_is_Quadrupole) 
        *(idat.t1) -= cross(uz_i, duduz_i);
      if (i_is_Quadrupole) {
        *(idat.t1) -= cross(ux_i, dudux_i);
        *(idat.t1) -= cross(uy_i, duduy_i);
      }
      
      if (j_is_Dipole || j_is_Quadrupole) 
        *(idat.t2) -= cross(uz_j, duduz_j);
      if (j_is_Quadrupole) {
        *(idat.t2) -= cross(uz_j, dudux_j);
        *(idat.t2) -= cross(uz_j, duduy_j);
      }

    } else {

      // only accumulate the forces and torques resulting from the
      // indirect reaction field terms.

      *(idat.vpair) += indirect_vpair;
      (*(idat.pot))[ELECTROSTATIC_FAMILY] += indirect_Pot;
      *(idat.f1) += indirect_dVdr;
      
      if (i_is_Dipole) 
        *(idat.t1) -= cross(uz_i, indirect_duduz_i);
      if (j_is_Dipole) 
        *(idat.t2) -= cross(uz_j, indirect_duduz_j);
    }


    return;
  }  
    
  void Electrostatic::calcSelfCorrection(SelfData &sdat) {
    RealType mu1, preVal, chg1, self;
    
    if (!initialized_) initialize();

    ElectrostaticAtomData data = ElectrostaticMap[sdat.atype];
   
    // logicals
    bool i_is_Charge = data.is_Charge;
    bool i_is_Dipole = data.is_Dipole;

    if (summationMethod_ == esm_REACTION_FIELD) {
      if (i_is_Dipole) {
        mu1 = data.dipole_moment;          
        preVal = pre22_ * preRF2_ * mu1 * mu1;
        (*(sdat.pot))[ELECTROSTATIC_FAMILY] -= 0.5 * preVal;
        
        // The self-correction term adds into the reaction field vector
        Vector3d uz_i = sdat.eFrame->getColumn(2);
        Vector3d ei = preVal * uz_i;

        // This looks very wrong.  A vector crossed with itself is zero.
        *(sdat.t) -= cross(uz_i, ei);
      }
    } else if (summationMethod_ == esm_SHIFTED_FORCE || summationMethod_ == esm_SHIFTED_POTENTIAL) {
      if (i_is_Charge) {        
        chg1 = data.fixedCharge;
        if (screeningMethod_ == DAMPED) {
          self = - 0.5 * (c1c_ + alphaPi_) * chg1 * (chg1 + *(sdat.skippedCharge)) * pre11_;
        } else {         
          self = - 0.5 * rcuti_ * chg1 * (chg1 +  *(sdat.skippedCharge)) * pre11_;
        }
        (*(sdat.pot))[ELECTROSTATIC_FAMILY] += self;
      }
    }
  }

  RealType Electrostatic::getSuggestedCutoffRadius(pair<AtomType*, AtomType*> atypes) {
    // This seems to work moderately well as a default.  There's no
    // inherent scale for 1/r interactions that we can standardize.
    // 12 angstroms seems to be a reasonably good guess for most
    // cases.
    return 12.0;
  }
}
