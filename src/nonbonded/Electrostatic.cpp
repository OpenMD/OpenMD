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
#include "math/erfc.hpp"
#include "math/SquareMatrix.hpp"

namespace OpenMD {
  
  Electrostatic::Electrostatic(): name_("Electrostatic"), initialized_(false),
                                  forceField_(NULL), info_(NULL), 
                                  haveCutoffRadius_(false),
                                  haveDampingAlpha_(false), 
                                  haveDielectric_(false),
                                  haveElectroSplines_(false)
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
    // Dipole-Dipole, assuming dipoles are measured in Debye
    pre22_ = 14.39325;
    // Charge-Quadrupole, assuming charges are measured in electrons, and
    // quadrupoles are measured in 10^-26 esu cm^2
    // This unit is also known affectionately as an esu centibarn.
    pre14_ = 69.13373;
    // Dipole-Quadrupole, assuming dipoles are measured in debyes and
    // quadrupoles in esu centibarns:
    pre24_ = 14.39325;
    // Quadrupole-Quadrupole, assuming esu centibarns:
    pre44_ = 14.39325;

    // conversions for the simulation box dipole moment
    chargeToC_ = 1.60217733e-19;
    angstromToM_ = 1.0e-10;
    debyeToCm_ = 3.33564095198e-30;
    
    // Default number of points for electrostatic splines
    np_ = 140;
    
    // variables to handle different summation methods for long-range 
    // electrostatics:
    summationMethod_ = esm_HARD;    
    screeningMethod_ = UNDAMPED;
    dielectric_ = 1.0;
  
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
                 "Electrostatic::initialize: dampingAlpha was not specified in the\n"
                 "\tinput file.  A default value of %f (1/ang) will be used for the\n"
                 "\tcutoff of %f (ang).\n", 
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
    
    if (summationMethod_ == esm_REACTION_FIELD) {
      preRF_ = (dielectric_ - 1.0) / 
        ((2.0 * dielectric_ + 1.0) * pow(cutoffRadius_,3) );
    }
    
    RealType b0c, b1c, b2c, b3c, b4c, b5c;
    RealType db0c_1, db0c_2, db0c_3, db0c_4, db0c_5;
    RealType a2, expTerm, invArootPi;
    
    RealType r = cutoffRadius_;
    RealType r2 = r * r;

    if (screeningMethod_ == DAMPED) {      
      a2 = dampingAlpha_ * dampingAlpha_;
      invArootPi = 1.0 / (dampingAlpha_ * sqrt(M_PI));    
      expTerm = exp(-a2 * r2);
      // values of Smith's B_l functions at the cutoff radius:
      b0c = erfc(dampingAlpha_ * r) / r;
      b1c = (      b0c     + 2.0*a2     * expTerm * invArootPi) / r2;
      b2c = (3.0 * b1c + pow(2.0*a2, 2) * expTerm * invArootPi) / r2;
      b3c = (5.0 * b2c + pow(2.0*a2, 3) * expTerm * invArootPi) / r2;
      b4c = (7.0 * b3c + pow(2.0*a2, 4) * expTerm * invArootPi) / r2;
      b5c = (9.0 * b4c + pow(2.0*a2, 5) * expTerm * invArootPi) / r2;
      selfMult_ = b0c + a2 * invArootPi;
    } else {
      a2 = 0.0;
      b0c = 1.0 / r;
      b1c = (      b0c) / r2;
      b2c = (3.0 * b1c) / r2;
      b3c = (5.0 * b2c) / r2;
      b4c = (7.0 * b3c) / r2;
      b5c = (9.0 * b4c) / r2;
      selfMult_ = b0c;
    }

    // higher derivatives of B_0 at the cutoff radius:
    db0c_1 = -r * b1c;
    db0c_2 =     -b1c + r2 * b2c;
    db0c_3 =          3.0*r*b2c  - r2*r*b3c;
    db0c_4 =          3.0*b2c  - 6.0*r2*b3c     + r2*r2*b4c;
    db0c_5 =                    -15.0*r*b3c + 10.0*r2*r*b4c - r2*r2*r*b5c;
    
    // working variables for the splines:
    RealType ri, ri2;
    RealType b0, b1, b2, b3, b4, b5;
    RealType db0_1, db0_2, db0_3, db0_4, db0_5;
    RealType f0;
    RealType g0, g1, g2, g3, g4;
    RealType h1, h2, h3, h4;
    RealType s2, s3, s4;
    RealType t3, t4;
    RealType u4;

    // working variables for Taylor expansion:
    RealType rmRc, rmRc2, rmRc3, rmRc4;

    // Approximate using splines using a maximum of 0.1 Angstroms
    // between points.
    int nptest = int((cutoffRadius_ + 2.0) / 0.1);
    np_ = (np_ > nptest) ? np_ : nptest;
  
    // Add a 2 angstrom safety window to deal with cutoffGroups that
    // have charged atoms longer than the cutoffRadius away from each
    // other.  Splining is almost certainly the best choice here.
    // Direct calls to erfc would be preferrable if it is a very fast
    // implementation.

    RealType dx = (cutoffRadius_ + 2.0) / RealType(np_);

    // Storage vectors for the computed functions    
    vector<RealType> rv;
    vector<RealType> v01v, v02v;
    vector<RealType> v11v, v12v, v13v;
    vector<RealType> v21v, v22v, v23v, v24v;
    vector<RealType> v31v, v32v, v33v, v34v, v35v;
    vector<RealType> v41v, v42v, v43v, v44v, v45v, v46v;

    for (int i = 1; i < np_ + 1; i++) {
      r = RealType(i) * dx;
      rv.push_back(r);

      ri = 1.0 / r;
      ri2 = ri * ri;

      r2 = r * r;
      expTerm = exp(-a2 * r2);

      // Taylor expansion factors (no need for factorials this way):
      rmRc = r - cutoffRadius_;
      rmRc2 = rmRc  * rmRc / 2.0;
      rmRc3 = rmRc2 * rmRc / 3.0;
      rmRc4 = rmRc3 * rmRc / 4.0;

      // values of Smith's B_l functions at r:
      if (screeningMethod_ == DAMPED) {             
        b0 = erfc(dampingAlpha_ * r) * ri;
        b1 = (      b0 +     2.0*a2     * expTerm * invArootPi) * ri2;
        b2 = (3.0 * b1 + pow(2.0*a2, 2) * expTerm * invArootPi) * ri2;
        b3 = (5.0 * b2 + pow(2.0*a2, 3) * expTerm * invArootPi) * ri2;
        b4 = (7.0 * b3 + pow(2.0*a2, 4) * expTerm * invArootPi) * ri2;
        b5 = (9.0 * b4 + pow(2.0*a2, 5) * expTerm * invArootPi) * ri2;
      } else {
        b0 = ri;
        b1 = (      b0) * ri2;
        b2 = (3.0 * b1) * ri2;
        b3 = (5.0 * b2) * ri2;
        b4 = (7.0 * b3) * ri2;
        b5 = (9.0 * b4) * ri2;
      }
                
      // higher derivatives of B_0 at r:
      db0_1 = -r * b1;
      db0_2 =     -b1 + r2 * b2;
      db0_3 =          3.0*r*b2   - r2*r*b3;
      db0_4 =          3.0*b2   - 6.0*r2*b3     + r2*r2*b4;
      db0_5 =                    -15.0*r*b3 + 10.0*r2*r*b4 - r2*r2*r*b5;


      switch (summationMethod_) {
      case esm_SHIFTED_FORCE:
        f0 = b0 - b0c - rmRc*db0c_1;
        
        g0 = db0_1 - db0c_1;
        g1 = g0 - rmRc *db0c_2;
        g2 = g1 - rmRc2*db0c_3;
        g3 = g2 - rmRc3*db0c_4;
        g4 = g3 - rmRc4*db0c_5;
        
        h1 = db0_2 - db0c_2;
        h2 = h1 - rmRc *db0c_3;
        h3 = h2 - rmRc2*db0c_4;
        h4 = h3 - rmRc3*db0c_5;
        
        s2 = db0_3 - db0c_3;
        s3 = s2 - rmRc *db0c_4;
        s4 = s3 - rmRc2*db0c_5;
        
        t3 = db0_4 - db0c_4;
        t4 = t3 - rmRc *db0c_5;
        
        u4 = db0_5 - db0c_5;
        break;

      case esm_SHIFTED_POTENTIAL:
        f0 = b0 - b0c;
        
        g0 = db0_1;
        g1 = db0_1 - db0c_1;
        g2 = g1 - rmRc *db0c_2;
        g3 = g2 - rmRc2*db0c_3;
        g4 = g3 - rmRc3*db0c_4;

        h1 = db0_2;
        h2 = db0_2 - db0c_2;
        h3 = h2 - rmRc *db0c_3;
        h4 = h3 - rmRc2*db0c_4;
        
        s2 = db0_3;
        s3 = db0_3 - db0c_3;
        s4 = s3 - rmRc *db0c_4;

        t3 = db0_4;
        t4 = db0_4 - db0c_4;
        
        u4 = db0_5;
        break;

      case esm_SWITCHING_FUNCTION:
      case esm_HARD:
        f0 = b0;
        
        g0 = db0_1;
        g1 = g0;
        g2 = g1;
        g3 = g2;
        g4 = g3;
        
        h1 = db0_2;
        h2 = h1;
        h3 = h2;
        h4 = h3;
        
        s2 = db0_3;
        s3 = s2;
        s4 = s3;
        
        t3 = db0_4;
        t4 = t3;
        
        u4 = db0_5;
        break;

      case esm_REACTION_FIELD:

        // following DL_POLY's lead for shifting the image charge potential:
        f0 = b0  + preRF_ * r2
          - (b0c + preRF_ * cutoffRadius_ * cutoffRadius_);

        g0 = db0_1 + preRF_ * 2.0 * r;
        g1 = g0;
        g2 = g1;
        g3 = g2;
        g4 = g3;
 
        h1 = db0_2 + preRF_ * 2.0;
        h2 = h1;
        h3 = h2;
        h4 = h3;

        s2 = db0_3;
        s3 = s2;
        s4 = s3;
        
        t3 = db0_4;
        t4 = t3;
        
        u4 = db0_5;        
        break;
                
      case esm_EWALD_FULL:
      case esm_EWALD_PME:
      case esm_EWALD_SPME:
      default :
        map<string, ElectrostaticSummationMethod>::iterator i;
        std::string meth;
        for (i = summationMap_.begin(); i != summationMap_.end(); ++i) {
          if ((*i).second == summationMethod_) meth = (*i).first;
        }
        sprintf( painCave.errMsg,
                 "Electrostatic::initialize: electrostaticSummationMethod %s \n"
                 "\thas not been implemented yet. Please select one of:\n"
                 "\t\"hard\", \"shifted_potential\", or \"shifted_force\"\n", 
                 meth.c_str() );
        painCave.isFatal = 1;
        simError();
        break;       
      }

      v01 = f0;
      v02 = g0;

      v11 = g1;
      v12 = g1 * ri;
      v13 = h1 - v12;

      v21 = g2 * ri;
      v22 = h2 - v21;
      v23 = v22 * ri;
      v24 = s2 - 3.0*v23;        

      v31 = (h3 - g3 * ri) * ri;
      v32 = s3 - 3.0*v31;
      v33 = v31 * ri;
      v34 = v32 * ri;
      v35 = t3 - 6.0*v34 - 3.0*v33;

      v41 = (h4 - g4 * ri) * ri2;
      v42 = s4 * ri - 3.0*v41;
      v43 = t4 - 6.0*v42 - 3.0*v41;
      v44 = v42 * ri;
      v45 = v43 * ri;
      v46 = u4 - 10.0*v45 - 15.0*v44;

      // Add these computed values to the storage vectors for spline creation:
      v01v.push_back(v01);
      v02v.push_back(v02);

      v11v.push_back(v11);
      v12v.push_back(v12);
      v13v.push_back(v13);

      v21v.push_back(v21);
      v22v.push_back(v22);
      v23v.push_back(v23);
      v24v.push_back(v24);

      v31v.push_back(v31);
      v32v.push_back(v32);
      v33v.push_back(v33);
      v34v.push_back(v34);
      v35v.push_back(v35);
      
      v41v.push_back(v41);
      v42v.push_back(v42);
      v43v.push_back(v43);
      v44v.push_back(v44);
      v45v.push_back(v45);
      v46v.push_back(v46);
    }

    // construct the spline structures and fill them with the values we've
    // computed:

    v01s = new CubicSpline();
    v01s->addPoints(rv, v01v);
    v02s = new CubicSpline();
    v02s->addPoints(rv, v02v);

    v11s = new CubicSpline();
    v11s->addPoints(rv, v11v);
    v12s = new CubicSpline();
    v12s->addPoints(rv, v12v);
    v13s = new CubicSpline();
    v13s->addPoints(rv, v13v);

    v21s = new CubicSpline();
    v21s->addPoints(rv, v21v);
    v22s = new CubicSpline();
    v22s->addPoints(rv, v22v);
    v23s = new CubicSpline();
    v23s->addPoints(rv, v23v);
    v24s = new CubicSpline();
    v24s->addPoints(rv, v24v);

    v31s = new CubicSpline();
    v31s->addPoints(rv, v31v);
    v32s = new CubicSpline();
    v32s->addPoints(rv, v32v);
    v33s = new CubicSpline();
    v33s->addPoints(rv, v33v);
    v34s = new CubicSpline();
    v34s->addPoints(rv, v34v);
    v35s = new CubicSpline();
    v35s->addPoints(rv, v35v);

    v41s = new CubicSpline();
    v41s->addPoints(rv, v41v);
    v42s = new CubicSpline();
    v42s->addPoints(rv, v42v);
    v43s = new CubicSpline();
    v43s->addPoints(rv, v43v);
    v44s = new CubicSpline();
    v44s->addPoints(rv, v44v);
    v45s = new CubicSpline();
    v45s->addPoints(rv, v45v);
    v46s = new CubicSpline();
    v46s->addPoints(rv, v46v);

    haveElectroSplines_ = true;

    initialized_ = true;
  }
      
  void Electrostatic::addType(AtomType* atomType){

    ElectrostaticAtomData electrostaticAtomData;
    electrostaticAtomData.is_Charge = false;
    electrostaticAtomData.is_Dipole = false;
    electrostaticAtomData.is_Quadrupole = false;
    electrostaticAtomData.is_Fluctuating = false;

    FixedChargeAdapter fca = FixedChargeAdapter(atomType);

    if (fca.isFixedCharge()) {
      electrostaticAtomData.is_Charge = true;
      electrostaticAtomData.fixedCharge = fca.getCharge();
    }

    MultipoleAdapter ma = MultipoleAdapter(atomType);
    if (ma.isMultipole()) {
      if (ma.isDipole()) {
        electrostaticAtomData.is_Dipole = true;
        electrostaticAtomData.dipole = ma.getDipole();
      }
      if (ma.isQuadrupole()) {
        electrostaticAtomData.is_Quadrupole = true;
        electrostaticAtomData.quadrupole = ma.getQuadrupole();
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
        vector<RealType> Jvals;
        // don't start at i = 0, as rval = 0 is undefined for the
        // slater overlap integrals.
        for (int i = 1; i < np_+1; i++) {
          rval = RealType(i) * dr;
          rvals.push_back(rval);
          Jvals.push_back(sSTOCoulInt( a, b, m, n, rval * 
                                       PhysicalConstants::angstromToBohr ) * 
                          PhysicalConstants::hartreeToKcal );
        }
        
        CubicSpline* J = new CubicSpline();
        J->addPoints(rvals, Jvals);
        
        pair<AtomType*, AtomType*> key1, key2;
        key1 = make_pair(atomType, atype2);
        key2 = make_pair(atype2, atomType);
        
        Jij[key1] = J;
        Jij[key2] = J;
      }
    }
 
    return;
  }
  
  void Electrostatic::setCutoffRadius( RealType rCut ) {
    cutoffRadius_ = rCut;
    haveCutoffRadius_ = true;
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

    RealType C_a, C_b;  // Charges 
    Vector3d D_a, D_b;  // Dipoles (space-fixed)
    Mat3x3d  Q_a, Q_b;  // Quadrupoles (space-fixed)

    RealType ri;                                 // Distance utility scalar
    RealType rdDa, rdDb;                         // Dipole utility scalars
    Vector3d rxDa, rxDb;                         // Dipole utility vectors
    RealType rdQar, rdQbr, trQa, trQb;           // Quadrupole utility scalars
    Vector3d Qar, Qbr, rQa, rQb, rxQar, rxQbr;   // Quadrupole utility vectors
    RealType pref;

    RealType DadDb, trQaQb, DadQbr, DbdQar;       // Cross-interaction scalars
    RealType rQaQbr;
    Vector3d DaxDb, DadQb, DbdQa, DaxQbr, DbxQar; // Cross-interaction vectors
    Vector3d rQaQb, QaQbr, QaxQb, rQaxQbr;
    Mat3x3d  QaQb;                                // Cross-interaction matrices

    RealType U(0.0);  // Potential
    Vector3d F(0.0);  // Force
    Vector3d Ta(0.0); // Torque on site a
    Vector3d Tb(0.0); // Torque on site b
    Vector3d Ea(0.0); // Electric field at site a
    Vector3d Eb(0.0); // Electric field at site b
    RealType dUdCa(0.0); // fluctuating charge force at site a
    RealType dUdCb(0.0); // fluctuating charge force at site a
    
    // Indirect interactions mediated by the reaction field.
    RealType indirect_Pot(0.0);  // Potential
    Vector3d indirect_F(0.0);    // Force
    Vector3d indirect_Ta(0.0);   // Torque on site a
    Vector3d indirect_Tb(0.0);   // Torque on site b

    // Excluded potential that is still computed for fluctuating charges
    RealType excluded_Pot(0.0);

    RealType rfContrib, coulInt;
    
    // spline for coulomb integral
    CubicSpline* J;

    if (!initialized_) initialize();
    
    ElectrostaticAtomData data1 = ElectrostaticMap[idat.atypes.first];
    ElectrostaticAtomData data2 = ElectrostaticMap[idat.atypes.second];
    
    // some variables we'll need independent of electrostatic type:

    ri = 1.0 /  *(idat.rij);
    Vector3d rhat =  *(idat.d)  * ri;
      
    // logicals

    bool a_is_Charge = data1.is_Charge;
    bool a_is_Dipole = data1.is_Dipole;
    bool a_is_Quadrupole = data1.is_Quadrupole;
    bool a_is_Fluctuating = data1.is_Fluctuating;

    bool b_is_Charge = data2.is_Charge;
    bool b_is_Dipole = data2.is_Dipole;
    bool b_is_Quadrupole = data2.is_Quadrupole;
    bool b_is_Fluctuating = data2.is_Fluctuating;

    // Obtain all of the required radial function values from the
    // spline structures:
    
    // needed for fields (and forces):
    if (a_is_Charge || b_is_Charge) {
      v02 = v02s->getValueAt( *(idat.rij) );
    }
    if (a_is_Dipole || b_is_Dipole) {
      v12 = v12s->getValueAt( *(idat.rij) );
      v13 = v13s->getValueAt( *(idat.rij) );
    }
    if (a_is_Quadrupole || b_is_Quadrupole) {
      v23 = v23s->getValueAt( *(idat.rij) );
      v24 = v24s->getValueAt( *(idat.rij) );
    }

    // needed for potentials (and torques):
    if (a_is_Charge && b_is_Charge) {
      v01 = v01s->getValueAt( *(idat.rij) );
    }
    if ((a_is_Charge && b_is_Dipole) || (b_is_Charge && a_is_Dipole)) {
      v11 = v11s->getValueAt( *(idat.rij) );
    }
    if ((a_is_Charge && b_is_Quadrupole) || (b_is_Charge && a_is_Quadrupole)) {
      v21 = v21s->getValueAt( *(idat.rij) );
      v22 = v22s->getValueAt( *(idat.rij) );
    } else if (a_is_Dipole && b_is_Dipole) {
      v21 = v21s->getValueAt( *(idat.rij) );
      v22 = v22s->getValueAt( *(idat.rij) );
      v23 = v23s->getValueAt( *(idat.rij) );
      v24 = v24s->getValueAt( *(idat.rij) );
    }
    if ((a_is_Dipole && b_is_Quadrupole) || 
        (b_is_Dipole && a_is_Quadrupole)) {
      v31 = v31s->getValueAt( *(idat.rij) );
      v32 = v32s->getValueAt( *(idat.rij) );
      v33 = v33s->getValueAt( *(idat.rij) );
      v34 = v34s->getValueAt( *(idat.rij) );
      v35 = v35s->getValueAt( *(idat.rij) );
    }
    if (a_is_Quadrupole && b_is_Quadrupole) {
      v41 = v41s->getValueAt( *(idat.rij) );
      v42 = v42s->getValueAt( *(idat.rij) );
      v43 = v43s->getValueAt( *(idat.rij) );
      v44 = v44s->getValueAt( *(idat.rij) );
      v45 = v45s->getValueAt( *(idat.rij) );
      v46 = v46s->getValueAt( *(idat.rij) );
    }

    // calculate the single-site contributions (fields, etc).
    
    if (a_is_Charge) {
      C_a = data1.fixedCharge;
      
      if (a_is_Fluctuating) {
        C_a += *(idat.flucQ1);
      }
      
      if (idat.excluded) {
        *(idat.skippedCharge2) += C_a;
      } else {
        // only do the field if we're not excluded:
        Eb -= C_a *  pre11_ * v02 * rhat;
      }
    }
    
    if (a_is_Dipole) {
      D_a = *(idat.dipole1);
      rdDa = dot(rhat, D_a);
      rxDa = cross(rhat, D_a);
      if (!idat.excluded) 
        Eb -=  pre12_ * (v13 * rdDa * rhat + v12 * D_a);
    }
    
    if (a_is_Quadrupole) {
      Q_a = *(idat.quadrupole1);
      trQa =  Q_a.trace();
      Qar =   Q_a * rhat;
      rQa = rhat * Q_a;
      rdQar = dot(rhat, Qar);
      rxQar = cross(rhat, Qar);
      if (!idat.excluded) 
        Eb -= pre14_ * ((trQa * rhat + 2.0 * Qar) * v23 + rdQar * rhat * v24);
    }
    
    if (b_is_Charge) {
      C_b = data2.fixedCharge;
      
      if (b_is_Fluctuating) 
        C_b += *(idat.flucQ2);
      
      if (idat.excluded) {
        *(idat.skippedCharge1) += C_b;
      } else {
        // only do the field if we're not excluded:
        Ea += C_b *  pre11_ * v02 * rhat;
      }
    }
    
    if (b_is_Dipole) {
      D_b = *(idat.dipole2);
      rdDb = dot(rhat, D_b);
      rxDb = cross(rhat, D_b);
      if (!idat.excluded) 
        Ea += pre12_ * (v13 * rdDb * rhat + v12 * D_b);
    }
    
    if (b_is_Quadrupole) {
      Q_b = *(idat.quadrupole2);
      trQb =  Q_b.trace();
      Qbr =   Q_b * rhat;
      rQb = rhat * Q_b;
      rdQbr = dot(rhat, Qbr);
      rxQbr = cross(rhat, Qbr);
      if (!idat.excluded) 
        Ea += pre14_ * ((trQb * rhat + 2.0 * Qbr) * v23 + rdQbr * rhat * v24);
    }
    
    if ((a_is_Fluctuating || b_is_Fluctuating) && idat.excluded) {
      J = Jij[idat.atypes];
    }    
    
    if (a_is_Charge) {     
      
      if (b_is_Charge) {
        pref =  pre11_ * *(idat.electroMult);          
        U  += C_a * C_b * pref * v01;
        F  += C_a * C_b * pref * v02 * rhat;
        
        // If this is an excluded pair, there are still indirect
        // interactions via the reaction field we must worry about:

        if (summationMethod_ == esm_REACTION_FIELD && idat.excluded) {
          rfContrib = preRF_ * pref * C_a * C_b * *(idat.r2);
          indirect_Pot += rfContrib;
          indirect_F   += rfContrib * 2.0 * ri * rhat;
        }
        
        // Fluctuating charge forces are handled via Coulomb integrals
        // for excluded pairs (i.e. those connected via bonds) and
        // with the standard charge-charge interaction otherwise.

        if (idat.excluded) {          
          if (a_is_Fluctuating || b_is_Fluctuating) {
            coulInt = J->getValueAt( *(idat.rij) );
            if (a_is_Fluctuating)  dUdCa += coulInt * C_b;
            if (b_is_Fluctuating)  dUdCb += coulInt * C_a;
            excluded_Pot += C_a * C_b * coulInt;
          }          
        } else {
          if (a_is_Fluctuating) dUdCa += C_b * pref * v01;
          if (a_is_Fluctuating) dUdCb += C_a * pref * v01;
        }
      }

      if (b_is_Dipole) {
        pref =  pre12_ * *(idat.electroMult);        
        U  += C_a * pref * v11 * rdDb;
        F  += C_a * pref * (v13 * rdDb * rhat + v12 * D_b);
        Tb += C_a * pref * v11 * rxDb;

        if (a_is_Fluctuating) dUdCa += pref * v11 * rdDb;

        // Even if we excluded this pair from direct interactions, we
        // still have the reaction-field-mediated charge-dipole
        // interaction:

        if (summationMethod_ == esm_REACTION_FIELD && idat.excluded) {
          rfContrib = C_a * pref * preRF_ * 2.0 * *(idat.rij);
          indirect_Pot += rfContrib * rdDb;
          indirect_F   += rfContrib * D_b / (*idat.rij);
          indirect_Tb  += C_a * pref * preRF_ * rxDb;
        }
      }

      if (b_is_Quadrupole) {
        pref = pre14_ * *(idat.electroMult);
        U  +=  C_a * pref * (v21 * trQb + v22 * rdQbr);
        F  +=  C_a * pref * (trQb * rhat + 2.0 * Qbr) * v23;
        F  +=  C_a * pref * rdQbr * rhat * v24;
        Tb +=  C_a * pref * 2.0 * rxQbr * v22;

        if (a_is_Fluctuating) dUdCa += pref * (v21 * trQb + v22 * rdQbr);
      }
    }

    if (a_is_Dipole) {

      if (b_is_Charge) {
        pref = pre12_ * *(idat.electroMult);

        U  -= C_b * pref * v11 * rdDa;
        F  -= C_b * pref * (v13 * rdDa * rhat + v12 * D_a);
        Ta -= C_b * pref * v11 * rxDa;

        if (b_is_Fluctuating) dUdCb -= pref * v11 * rdDa;

        // Even if we excluded this pair from direct interactions,
        // we still have the reaction-field-mediated charge-dipole
        // interaction:
        if (summationMethod_ == esm_REACTION_FIELD && idat.excluded) {
          rfContrib = C_b * pref * preRF_ * 2.0 * *(idat.rij);
          indirect_Pot -= rfContrib * rdDa;
          indirect_F   -= rfContrib * D_a / (*idat.rij);
          indirect_Ta  -= C_b * pref * preRF_ * rxDa;
        }
      }

      if (b_is_Dipole) {
        pref = pre22_ * *(idat.electroMult);
        DadDb = dot(D_a, D_b);
        DaxDb = cross(D_a, D_b);

        U  -= pref * (DadDb * v21 + rdDa * rdDb * v22);
        F  -= pref * (DadDb * rhat + rdDb * D_a + rdDa * D_b)*v23;
        F  -= pref * (rdDa * rdDb) * v24 * rhat;
        Ta += pref * ( v21 * DaxDb - v22 * rdDb * rxDa);
        Tb += pref * (-v21 * DaxDb - v22 * rdDa * rxDb);

        // Even if we excluded this pair from direct interactions, we
        // still have the reaction-field-mediated dipole-dipole
        // interaction:
        if (summationMethod_ == esm_REACTION_FIELD && idat.excluded) {
          rfContrib = -pref * preRF_ * 2.0;
          indirect_Pot += rfContrib * DadDb;
          indirect_Ta  += rfContrib * DaxDb;
          indirect_Tb  -= rfContrib * DaxDb;
        }

      }

      if (b_is_Quadrupole) {
        pref = pre24_ * *(idat.electroMult);
        DadQb = D_a * Q_b;
        DadQbr = dot(D_a, Qbr);
        DaxQbr = cross(D_a, Qbr);

        U  -= pref * ((trQb*rdDa + 2.0*DadQbr)*v31 + rdDa*rdQbr*v32);
        F  -= pref * (trQb*D_a + 2.0*DadQb) * v33;
        F  -= pref * (trQb*rdDa*rhat + 2.0*DadQbr*rhat + D_a*rdQbr 
                      + 2.0*rdDa*rQb)*v34;
        F  -= pref * (rdDa * rdQbr * rhat * v35);
        Ta += pref * ((-trQb*rxDa + 2.0 * DaxQbr)*v31 - rxDa*rdQbr*v32);
        Tb += pref * ((2.0*cross(DadQb, rhat) - 2.0*DaxQbr)*v31 
                      - 2.0*rdDa*rxQbr*v32);
      }
    }

    if (a_is_Quadrupole) {
      if (b_is_Charge) {
        pref = pre14_ * *(idat.electroMult);
        U  += C_b * pref * (v21 * trQa + v22 * rdQar);
        F  += C_b * pref * (trQa * rhat + 2.0 * Qar) * v23;
        F  += C_b * pref * rdQar * rhat * v24;
        Ta += C_b * pref * 2.0 * rxQar * v22;

        if (b_is_Fluctuating) dUdCb += pref * (v21 * trQa + v22 * rdQar);
      }
      if (b_is_Dipole) {
        pref = pre24_ * *(idat.electroMult);
        DbdQa = D_b * Q_a;
        DbdQar = dot(D_b, Qar);
        DbxQar = cross(D_b, Qar);

        U  += pref * ((trQa*rdDb + 2.0*DbdQar)*v31 + rdDb*rdQar*v32);
        F  += pref * (trQa*D_b + 2.0*DbdQa) * v33;
        F  += pref * (trQa*rdDb*rhat + 2.0*DbdQar*rhat + D_b*rdQar 
                      + 2.0*rdDb*rQa)*v34;
        F  += pref * (rdDb * rdQar * rhat * v35);
        Ta += pref * ((-2.0*cross(DbdQa, rhat) + 2.0*DbxQar)*v31 
                      + 2.0*rdDb*rxQar*v32);
        Tb += pref * ((trQa*rxDb - 2.0 * DbxQar)*v31 + rxDb*rdQar*v32);
      }
      if (b_is_Quadrupole) {
        pref = pre44_ * *(idat.electroMult);  // yes
        QaQb = Q_a * Q_b;
        trQaQb = QaQb.trace();
        rQaQb = rhat * QaQb;
        QaQbr = QaQb * rhat;
        QaxQb = cross(Q_a, Q_b);
        rQaQbr = dot(rQa, Qbr);
        rQaxQbr = cross(rQa, Qbr);
        
        U  += pref * (trQa * trQb + 2.0 * trQaQb) * v41;
        U  += pref * (trQa * rdQbr + trQb * rdQar  + 4.0 * rQaQbr) * v42;
        U  += pref * (rdQar * rdQbr) * v43;

        F  += pref * rhat * (trQa * trQb + 2.0 * trQaQb)*v44;
        F  += pref * rhat * (trQa * rdQbr + trQb * rdQar + 4.0 * rQaQbr)*v45;
        F  += pref * rhat * (rdQar * rdQbr)*v46;

        F  += pref * 2.0 * (trQb*rQa + trQa*rQb)*v44;
        F  += pref * 4.0 * (rQaQb + QaQbr)*v44;
        F  += pref * 2.0 * (rQa*rdQbr + rdQar*rQb)*v45;

        Ta += pref * (- 4.0 * QaxQb * v41);
        Ta += pref * (- 2.0 * trQb * cross(rQa, rhat) 
                      + 4.0 * cross(rhat, QaQbr) 
                      - 4.0 * rQaxQbr) * v42;
        Ta += pref * 2.0 * cross(rhat,Qar) * rdQbr * v43;


        Tb += pref * (+ 4.0 * QaxQb * v41);
        Tb += pref * (- 2.0 * trQa * cross(rQb, rhat) 
                      - 4.0 * cross(rQaQb, rhat) 
                      + 4.0 * rQaxQbr) * v42;
        // Possible replacement for line 2 above:
        //             + 4.0 * cross(rhat, QbQar) 

        Tb += pref * 2.0 * cross(rhat,Qbr) * rdQar * v43;

        //  cerr << " tsum = " << Ta + Tb - cross(  *(idat.d) , F ) << "\n";
      }
    }

    if (idat.doElectricField) {
      *(idat.eField1) += Ea * *(idat.electroMult);
      *(idat.eField2) += Eb * *(idat.electroMult);
    }

    if (a_is_Fluctuating) *(idat.dVdFQ1) += dUdCa * *(idat.sw);
    if (b_is_Fluctuating) *(idat.dVdFQ2) += dUdCb * *(idat.sw);

    if (!idat.excluded) {
      
      *(idat.vpair) += U;
      (*(idat.pot))[ELECTROSTATIC_FAMILY] += U * *(idat.sw);
      *(idat.f1) += F * *(idat.sw);
      
      if (a_is_Dipole || a_is_Quadrupole) 
        *(idat.t1) += Ta * *(idat.sw);

      if (b_is_Dipole || b_is_Quadrupole) 
        *(idat.t2) += Tb * *(idat.sw);
      
    } else {

      // only accumulate the forces and torques resulting from the
      // indirect reaction field terms.

      *(idat.vpair) += indirect_Pot;      
      (*(idat.excludedPot))[ELECTROSTATIC_FAMILY] +=  excluded_Pot;
      (*(idat.pot))[ELECTROSTATIC_FAMILY] += *(idat.sw) * indirect_Pot;
      *(idat.f1) += *(idat.sw) * indirect_F;
      
      if (a_is_Dipole || a_is_Quadrupole) 
        *(idat.t1) += *(idat.sw) * indirect_Ta;
            
      if (b_is_Dipole || b_is_Quadrupole) 
        *(idat.t2) += *(idat.sw) * indirect_Tb;
    }
    return;
  }
    
  void Electrostatic::calcSelfCorrection(SelfData &sdat) {

    if (!initialized_) initialize();

    ElectrostaticAtomData data = ElectrostaticMap[sdat.atype];
    
    // logicals
    bool i_is_Charge = data.is_Charge;
    bool i_is_Dipole = data.is_Dipole;
    bool i_is_Fluctuating = data.is_Fluctuating;
    RealType C_a = data.fixedCharge;   
    RealType self, preVal, DadDa;
    
    if (i_is_Fluctuating) {
      C_a += *(sdat.flucQ);
      // dVdFQ is really a force, so this is negative the derivative
      *(sdat.dVdFQ) -=  *(sdat.flucQ) * data.hardness + data.electronegativity;
      (*(sdat.excludedPot))[ELECTROSTATIC_FAMILY] += (*sdat.flucQ) * 
        (*(sdat.flucQ) * data.hardness * 0.5 + data.electronegativity);
    }

    switch (summationMethod_) {
    case esm_REACTION_FIELD:
      
      if (i_is_Charge) {
        // Self potential [see Wang and Hermans, "Reaction Field
        // Molecular Dynamics Simulation with Friedman’s Image Charge
        // Method," J. Phys. Chem. 99, 12001-12007 (1995).]
        preVal = pre11_ * preRF_ * C_a * C_a;
        (*(sdat.pot))[ELECTROSTATIC_FAMILY] -= 0.5 * preVal / cutoffRadius_;
      }

      if (i_is_Dipole) {
        DadDa = data.dipole.lengthSquare();
        (*(sdat.pot))[ELECTROSTATIC_FAMILY] -= pre22_ * preRF_ * DadDa;
      }
      
      break;
      
    case esm_SHIFTED_FORCE:
    case esm_SHIFTED_POTENTIAL:
      if (i_is_Charge) {
        self = - selfMult_ * C_a * (C_a + *(sdat.skippedCharge)) * pre11_;
        (*(sdat.pot))[ELECTROSTATIC_FAMILY] += self;
      }
      break;
    default:
      break;
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
