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
 * [3]  Sun, Lin & Gezelter, J. Chem. Phys. 128, 234107 (2008).          
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
    summationMap_["TAYLOR_SHIFTED"]     = esm_TAYLOR_SHIFTED;     
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
    np_ = 100;
    
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
                 "\t\"shifted_potential\", \"shifted_force\",\n"
                 "\t\"taylor_shifted\", or \"reaction_field\".\n", 
                 myMethod.c_str() );
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

    Etypes.clear();
    Etids.clear();
    FQtypes.clear();
    FQtids.clear();
    ElectrostaticMap.clear();
    Jij.clear();
    nElectro_ = 0;
    nFlucq_ = 0;

    Etids.resize( forceField_->getNAtomType(), -1);
    FQtids.resize( forceField_->getNAtomType(), -1);

    set<AtomType*>::iterator at;
    for (at = simTypes_.begin(); at != simTypes_.end(); ++at) {     
      if ((*at)->isElectrostatic()) nElectro_++;
      if ((*at)->isFluctuatingCharge()) nFlucq_++;
    }
    
    Jij.resize(nFlucq_);

    for (at = simTypes_.begin(); at != simTypes_.end(); ++at) {
      if ((*at)->isElectrostatic()) addType(*at);
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
    RealType ric = 1.0 / r;
    RealType ric2 = ric * ric;

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
      //selfMult1_ = - 2.0 * a2 * invArootPi;
      //selfMult2_ = - 4.0 * a2 * a2 * invArootPi / 3.0;
      //selfMult4_ = - 8.0 * a2 * a2 * a2 * invArootPi / 5.0;
      // Half the Smith self piece:
      selfMult1_ = - a2 * invArootPi;
      selfMult2_ = - 2.0 * a2 * a2 * invArootPi / 3.0;
      selfMult4_ = - 4.0 * a2 * a2 * a2 * invArootPi / 5.0;
    } else {
      a2 = 0.0;
      b0c = 1.0 / r;
      b1c = (      b0c) / r2;
      b2c = (3.0 * b1c) / r2;
      b3c = (5.0 * b2c) / r2;
      b4c = (7.0 * b3c) / r2;
      b5c = (9.0 * b4c) / r2;
      selfMult1_ = 0.0;
      selfMult2_ = 0.0;
      selfMult4_ = 0.0;
    }

    // higher derivatives of B_0 at the cutoff radius:
    db0c_1 = -r * b1c;
    db0c_2 =     -b1c + r2 * b2c;
    db0c_3 =          3.0*r*b2c  - r2*r*b3c;
    db0c_4 =          3.0*b2c  - 6.0*r2*b3c     + r2*r2*b4c;
    db0c_5 =                    -15.0*r*b3c + 10.0*r2*r*b4c - r2*r2*r*b5c;   
    
    selfMult1_ -= b0c;
    selfMult2_ += (db0c_2 + 2.0*db0c_1*ric) /  3.0;
    selfMult4_ -= (db0c_4 + 4.0*db0c_3*ric) / 15.0;

    // working variables for the splines:
    RealType ri, ri2;
    RealType b0, b1, b2, b3, b4, b5;
    RealType db0_1, db0_2, db0_3, db0_4, db0_5;
    RealType f, fc, f0;
    RealType g, gc, g0, g1, g2, g3, g4;
    RealType h, hc, h1, h2, h3, h4;
    RealType s, sc, s2, s3, s4;
    RealType t, tc, t3, t4;
    RealType u, uc, u4;

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
    vector<RealType> v01v;
    vector<RealType> v11v;
    vector<RealType> v21v, v22v;
    vector<RealType> v31v, v32v;
    vector<RealType> v41v, v42v, v43v;

    /*
    vector<RealType> dv01v;
    vector<RealType> dv11v;
    vector<RealType> dv21v, dv22v;
    vector<RealType> dv31v, dv32v;
    vector<RealType> dv41v, dv42v, dv43v;
    */

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

      f = b0;
      fc = b0c;
      f0 = f - fc - rmRc*db0c_1;

      g = db0_1;        
      gc = db0c_1;
      g0 = g - gc;
      g1 = g0 - rmRc *db0c_2;
      g2 = g1 - rmRc2*db0c_3;
      g3 = g2 - rmRc3*db0c_4;
      g4 = g3 - rmRc4*db0c_5;

      h = db0_2;      
      hc = db0c_2;
      h1 = h - hc;
      h2 = h1 - rmRc *db0c_3;
      h3 = h2 - rmRc2*db0c_4;
      h4 = h3 - rmRc3*db0c_5;

      s = db0_3;      
      sc = db0c_3;
      s2 = s - sc;
      s3 = s2 - rmRc *db0c_4;
      s4 = s3 - rmRc2*db0c_5;

      t = db0_4;      
      tc = db0c_4;
      t3 = t - tc;
      t4 = t3 - rmRc *db0c_5;
      
      u = db0_5;         
      uc = db0c_5;
      u4 = u - uc;

      // in what follows below, the various v functions are used for
      // potentials and torques, while the w functions show up in the
      // forces.

      switch (summationMethod_) {
      case esm_SHIFTED_FORCE:
                
        v01 = f - fc - rmRc*gc;
        v11 = g - gc - rmRc*hc;
        v21 = g*ri - gc*ric - rmRc*(hc - gc*ric)*ric;
        v22 = h - g*ri - (hc - gc*ric) - rmRc*(sc - (hc - gc*ric)*ric);
        v31 = (h-g*ri)*ri - (hc-gc*ric)*ric - rmRc*(sc-2.0*(hc-gc*ric)*ric)*ric;
        v32 = (s - 3.0*(h-g*ri)*ri) - (sc - 3.0*(hc-gc*ric)*ric) 
          - rmRc*(tc - 3.0*(sc-2.0*(hc-gc*ric)*ric)*ric);
        v41 = (h - g*ri)*ri2 - (hc - gc*ric)*ric2 
          - rmRc*(sc - 3.0*(hc-gc*ric)*ric)*ric2;
        v42 = (s-3.0*(h-g*ri)*ri)*ri - (sc-3.0*(hc-gc*ric)*ric)*ric 
          - rmRc*(tc - (4.0*sc - 9.0*(hc - gc*ric)*ric)*ric)*ric;
        
        v43 = (t - 3.0*(2.0*s - 5.0*(h - g*ri)*ri)*ri) 
          - (tc - 3.0*(2.0*sc - 5.0*(hc - gc*ric)*ric)*ric)
          - rmRc*(uc-3.0*(2.0*tc - (7.0*sc - 15.0*(hc - gc*ric)*ric)*ric)*ric);

        dv01 = g - gc;
        dv11 = h - hc;
        dv21 = (h - g*ri)*ri - (hc - gc*ric)*ric;
        dv22 = (s - (h - g*ri)*ri) - (sc - (hc - gc*ric)*ric);        
        dv31 = (s - 2.0*(h-g*ri)*ri)*ri - (sc - 2.0*(hc-gc*ric)*ric)*ric;
        dv32 = (t - 3.0*(s-2.0*(h-g*ri)*ri)*ri) 
          - (tc - 3.0*(sc-2.0*(hc-gc*ric)*ric)*ric);
        dv41 = (s - 3.0*(h - g*ri)*ri)*ri2 - (sc - 3.0*(hc - gc*ric)*ric)*ric2;
        dv42 = (t - (4.0*s - 9.0*(h-g*ri)*ri)*ri)*ri 
          - (tc - (4.0*sc - 9.0*(hc-gc*ric)*ric)*ric)*ric; 
        dv43 = (u - 3.0*(2.0*t - (7.0*s - 15.0*(h - g*ri)*ri)*ri)*ri)
          - (uc - 3.0*(2.0*tc - (7.0*sc - 15.0*(hc - gc*ric)*ric)*ric)*ric);
        
        break;

      case esm_TAYLOR_SHIFTED:
        
        v01 = f0;
        v11 = g1;
        v21 = g2 * ri;
        v22 = h2 - v21;
        v31 = (h3 - g3 * ri) * ri;
        v32 = s3 - 3.0*v31;
        v41 = (h4 - g4 * ri) * ri2;
        v42 = s4 * ri - 3.0*v41;
        v43 = t4 - 6.0*v42 - 3.0*v41;

        dv01 = g0;
        dv11 = h1;
        dv21 = (h2 - g2*ri)*ri;
        dv22 = (s2 - (h2 - g2*ri)*ri);
        dv31 = (s3 - 2.0*(h3-g3*ri)*ri)*ri;
        dv32 = (t3 - 3.0*(s3-2.0*(h3-g3*ri)*ri)*ri); 
        dv41 = (s4 - 3.0*(h4 - g4*ri)*ri)*ri2;
        dv42 = (t4 - (4.0*s4 - 9.0*(h4-g4*ri)*ri)*ri)*ri; 
        dv43 = (u4 - 3.0*(2.0*t4 - (7.0*s4 - 15.0*(h4 - g4*ri)*ri)*ri)*ri);

        break;

      case esm_SHIFTED_POTENTIAL:

        v01 = f - fc;
        v11 = g - gc;
        v21 = g*ri - gc*ric;
        v22 = h - g*ri - (hc - gc*ric);
        v31 = (h-g*ri)*ri - (hc-gc*ric)*ric;
        v32 = (s - 3.0*(h-g*ri)*ri) - (sc - 3.0*(hc-gc*ric)*ric);
        v41 = (h - g*ri)*ri2 - (hc - gc*ric)*ric2; 
        v42 = (s-3.0*(h-g*ri)*ri)*ri - (sc-3.0*(hc-gc*ric)*ric)*ric;        
        v43 = (t - 3.0*(2.0*s - 5.0*(h - g*ri)*ri)*ri) 
          - (tc - 3.0*(2.0*sc - 5.0*(hc - gc*ric)*ric)*ric);

        dv01 = g;
        dv11 = h;
        dv21 = (h - g*ri)*ri;
        dv22 = (s - (h - g*ri)*ri);
        dv31 = (s - 2.0*(h-g*ri)*ri)*ri;
        dv32 = (t - 3.0*(s-2.0*(h-g*ri)*ri)*ri); 
        dv41 = (s - 3.0*(h - g*ri)*ri)*ri2;
        dv42 = (t - (4.0*s - 9.0*(h-g*ri)*ri)*ri)*ri; 
        dv43 = (u - 3.0*(2.0*t - (7.0*s - 15.0*(h - g*ri)*ri)*ri)*ri);

        break;

      case esm_SWITCHING_FUNCTION:
      case esm_HARD:

        v01 = f;
        v11 = g;
        v21 = g*ri;
        v22 = h - g*ri;
        v31 = (h-g*ri)*ri;
        v32 = (s - 3.0*(h-g*ri)*ri);
        v41 = (h - g*ri)*ri2; 
        v42 = (s-3.0*(h-g*ri)*ri)*ri;        
        v43 = (t - 3.0*(2.0*s - 5.0*(h - g*ri)*ri)*ri);

        dv01 = g;
        dv11 = h;
        dv21 = (h - g*ri)*ri;
        dv22 = (s - (h - g*ri)*ri);
        dv31 = (s - 2.0*(h-g*ri)*ri)*ri;
        dv32 = (t - 3.0*(s-2.0*(h-g*ri)*ri)*ri); 
        dv41 = (s - 3.0*(h - g*ri)*ri)*ri2;
        dv42 = (t - (4.0*s - 9.0*(h-g*ri)*ri)*ri)*ri; 
        dv43 = (u - 3.0*(2.0*t - (7.0*s - 15.0*(h - g*ri)*ri)*ri)*ri);

        break;

      case esm_REACTION_FIELD:
        
        // following DL_POLY's lead for shifting the image charge potential:
        f = b0 + preRF_ * r2;
        fc = b0c + preRF_ * cutoffRadius_ * cutoffRadius_;

        g = db0_1 + preRF_ * 2.0 * r;        
        gc = db0c_1 + preRF_ * 2.0 * cutoffRadius_;

        h = db0_2 + preRF_ * 2.0;
        hc = db0c_2 + preRF_ * 2.0;

        v01 = f - fc;
        v11 = g - gc;
        v21 = g*ri - gc*ric;
        v22 = h - g*ri - (hc - gc*ric);
        v31 = (h-g*ri)*ri - (hc-gc*ric)*ric;
        v32 = (s - 3.0*(h-g*ri)*ri) - (sc - 3.0*(hc-gc*ric)*ric);
        v41 = (h - g*ri)*ri2 - (hc - gc*ric)*ric2; 
        v42 = (s-3.0*(h-g*ri)*ri)*ri - (sc-3.0*(hc-gc*ric)*ric)*ric;        
        v43 = (t - 3.0*(2.0*s - 5.0*(h - g*ri)*ri)*ri) 
          - (tc - 3.0*(2.0*sc - 5.0*(hc - gc*ric)*ric)*ric);

        dv01 = g;
        dv11 = h;
        dv21 = (h - g*ri)*ri;
        dv22 = (s - (h - g*ri)*ri);
        dv31 = (s - 2.0*(h-g*ri)*ri)*ri;
        dv32 = (t - 3.0*(s-2.0*(h-g*ri)*ri)*ri); 
        dv41 = (s - 3.0*(h - g*ri)*ri)*ri2;
        dv42 = (t - (4.0*s - 9.0*(h-g*ri)*ri)*ri)*ri; 
        dv43 = (u - 3.0*(2.0*t - (7.0*s - 15.0*(h - g*ri)*ri)*ri)*ri);

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

      // Add these computed values to the storage vectors for spline creation:
      v01v.push_back(v01);
      v11v.push_back(v11);
      v21v.push_back(v21);
      v22v.push_back(v22);
      v31v.push_back(v31);
      v32v.push_back(v32);      
      v41v.push_back(v41);
      v42v.push_back(v42);
      v43v.push_back(v43);
      /*
      dv01v.push_back(dv01);
      dv11v.push_back(dv11);
      dv21v.push_back(dv21);
      dv22v.push_back(dv22);
      dv31v.push_back(dv31);
      dv32v.push_back(dv32);      
      dv41v.push_back(dv41);
      dv42v.push_back(dv42);
      dv43v.push_back(dv43);
      */
    }

    // construct the spline structures and fill them with the values we've
    // computed:

    v01s = new CubicSpline();
    v01s->addPoints(rv, v01v);
    v11s = new CubicSpline();
    v11s->addPoints(rv, v11v);
    v21s = new CubicSpline();
    v21s->addPoints(rv, v21v);
    v22s = new CubicSpline();
    v22s->addPoints(rv, v22v);
    v31s = new CubicSpline();
    v31s->addPoints(rv, v31v);
    v32s = new CubicSpline();
    v32s->addPoints(rv, v32v);
    v41s = new CubicSpline();
    v41s->addPoints(rv, v41v);
    v42s = new CubicSpline();
    v42s->addPoints(rv, v42v);
    v43s = new CubicSpline();
    v43s->addPoints(rv, v43v);

    /*
    dv01s = new CubicSpline();
    dv01s->addPoints(rv, dv01v);
    dv11s = new CubicSpline();
    dv11s->addPoints(rv, dv11v);
    dv21s = new CubicSpline();
    dv21s->addPoints(rv, dv21v);
    dv22s = new CubicSpline();
    dv22s->addPoints(rv, dv22v);
    dv31s = new CubicSpline();
    dv31s->addPoints(rv, dv31v);
    dv32s = new CubicSpline();
    dv32s->addPoints(rv, dv32v);
    dv41s = new CubicSpline();
    dv41s->addPoints(rv, dv41v);
    dv42s = new CubicSpline();
    dv42s->addPoints(rv, dv42v);
    dv43s = new CubicSpline();
    dv43s->addPoints(rv, dv43v);
    */

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

    int atid = atomType->getIdent();
    int etid = Etypes.size();
    int fqtid = FQtypes.size();

    pair<set<int>::iterator,bool> ret;    
    ret = Etypes.insert( atid );
    if (ret.second == false) {
      sprintf( painCave.errMsg,
               "Electrostatic already had a previous entry with ident %d\n",
               atid);
      painCave.severity = OPENMD_INFO;
      painCave.isFatal = 0;
      simError();         
    }
    
    Etids[ atid ] = etid;
    ElectrostaticMap.push_back(electrostaticAtomData);

    if (electrostaticAtomData.is_Fluctuating) {
      ret = FQtypes.insert( atid );
      if (ret.second == false) {
        sprintf( painCave.errMsg,
                 "Electrostatic already had a previous fluctuating charge entry with ident %d\n",
                 atid );
        painCave.severity = OPENMD_INFO;
        painCave.isFatal = 0;
        simError();         
      }
      FQtids[atid] = fqtid;
      Jij[fqtid].resize(nFlucq_);

      // Now, iterate over all known fluctuating and add to the coulomb integral map:
      
      std::set<int>::iterator it;
      for( it = FQtypes.begin(); it != FQtypes.end(); ++it) {     
        int etid2 = Etids[ (*it) ];
        int fqtid2 = FQtids[ (*it) ];
        ElectrostaticAtomData eaData2 = ElectrostaticMap[ etid2 ];
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
        Jij[fqtid][fqtid2] = J;
        Jij[fqtid2].resize( nFlucq_ );
        Jij[fqtid2][fqtid] = J;
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

    if (!initialized_) initialize();
    
    data1 = ElectrostaticMap[Etids[idat.atid1]];
    data2 = ElectrostaticMap[Etids[idat.atid2]];

    U = 0.0;  // Potential
    F.zero();  // Force
    Ta.zero(); // Torque on site a
    Tb.zero(); // Torque on site b
    Ea.zero(); // Electric field at site a
    Eb.zero(); // Electric field at site b
    dUdCa = 0.0; // fluctuating charge force at site a
    dUdCb = 0.0; // fluctuating charge force at site a
    
    // Indirect interactions mediated by the reaction field.
    indirect_Pot = 0.0;   // Potential
    indirect_F.zero();    // Force
    indirect_Ta.zero();   // Torque on site a
    indirect_Tb.zero();   // Torque on site b

    // Excluded potential that is still computed for fluctuating charges
    excluded_Pot= 0.0;


    // some variables we'll need independent of electrostatic type:

    ri = 1.0 /  *(idat.rij);
    rhat =  *(idat.d)  * ri;
      
    // logicals

    a_is_Charge = data1.is_Charge;
    a_is_Dipole = data1.is_Dipole;
    a_is_Quadrupole = data1.is_Quadrupole;
    a_is_Fluctuating = data1.is_Fluctuating;

    b_is_Charge = data2.is_Charge;
    b_is_Dipole = data2.is_Dipole;
    b_is_Quadrupole = data2.is_Quadrupole;
    b_is_Fluctuating = data2.is_Fluctuating;

    // Obtain all of the required radial function values from the
    // spline structures:
    
    // needed for fields (and forces):
    if (a_is_Charge || b_is_Charge) {
      v01s->getValueAndDerivativeAt( *(idat.rij), v01, dv01);
    }
    if (a_is_Dipole || b_is_Dipole) {
      v11s->getValueAndDerivativeAt( *(idat.rij), v11, dv11);
      v11or = ri * v11;
    }
    if (a_is_Quadrupole || b_is_Quadrupole ||  (a_is_Dipole && b_is_Dipole)) {
      v21s->getValueAndDerivativeAt( *(idat.rij), v21, dv21);
      v22s->getValueAndDerivativeAt( *(idat.rij), v22, dv22);
      v22or = ri * v22;
    }      

    // needed for potentials (and forces and torques):
    if ((a_is_Dipole && b_is_Quadrupole) || 
        (b_is_Dipole && a_is_Quadrupole)) {
      v31s->getValueAndDerivativeAt( *(idat.rij), v31, dv31);
      v32s->getValueAndDerivativeAt( *(idat.rij), v32, dv32);
      v31or = v31 * ri;
      v32or = v32 * ri;
    }
    if (a_is_Quadrupole && b_is_Quadrupole) {
      v41s->getValueAndDerivativeAt( *(idat.rij), v41, dv41);
      v42s->getValueAndDerivativeAt( *(idat.rij), v42, dv42);
      v43s->getValueAndDerivativeAt( *(idat.rij), v43, dv43);
      v42or = v42 * ri;
      v43or = v43 * ri;
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
        Eb -= C_a *  pre11_ * dv01 * rhat;
      }
    }
    
    if (a_is_Dipole) {
      D_a = *(idat.dipole1);
      rdDa = dot(rhat, D_a);
      rxDa = cross(rhat, D_a);
      if (!idat.excluded) 
        Eb -=  pre12_ * ((dv11-v11or) * rdDa * rhat + v11or * D_a);
    }
    
    if (a_is_Quadrupole) {
      Q_a = *(idat.quadrupole1);
      trQa =  Q_a.trace();
      Qar =   Q_a * rhat;
      rQa = rhat * Q_a;
      rdQar = dot(rhat, Qar);
      rxQar = cross(rhat, Qar);
      if (!idat.excluded) 
        Eb -= pre14_ * (trQa * rhat * dv21 + 2.0 * Qar * v22or 
                        + rdQar * rhat * (dv22 - 2.0*v22or));
    }
    
    if (b_is_Charge) {
      C_b = data2.fixedCharge;
      
      if (b_is_Fluctuating) 
        C_b += *(idat.flucQ2);
      
      if (idat.excluded) {
        *(idat.skippedCharge1) += C_b;
      } else {
        // only do the field if we're not excluded:
        Ea += C_b *  pre11_ * dv01 * rhat;
      }
    }
    
    if (b_is_Dipole) {
      D_b = *(idat.dipole2);
      rdDb = dot(rhat, D_b);
      rxDb = cross(rhat, D_b);
      if (!idat.excluded) 
        Ea += pre12_ * ((dv11-v11or) * rdDb * rhat + v11or * D_b);
    }
    
    if (b_is_Quadrupole) {
      Q_b = *(idat.quadrupole2);
      trQb =  Q_b.trace();
      Qbr =   Q_b * rhat;
      rQb = rhat * Q_b;
      rdQbr = dot(rhat, Qbr);
      rxQbr = cross(rhat, Qbr);
      if (!idat.excluded) 
        Ea += pre14_ * (trQb * rhat * dv21 + 2.0 * Qbr * v22or 
                        + rdQbr * rhat * (dv22 - 2.0*v22or));
    }
    
    if ((a_is_Fluctuating || b_is_Fluctuating) && idat.excluded) {
      J = Jij[FQtids[idat.atid1]][FQtids[idat.atid2]];
    }    
    
    if (a_is_Charge) {     
      
      if (b_is_Charge) {
        pref =  pre11_ * *(idat.electroMult);      
        U  += C_a * C_b * pref * v01;
        F  += C_a * C_b * pref * dv01 * rhat;
        
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
        F  += C_a * pref * ((dv11 - v11or) * rdDb * rhat + v11or * D_b);
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
        F  +=  C_a * pref * (trQb * dv21 * rhat + 2.0 * Qbr * v22or);
        F  +=  C_a * pref * rdQbr * rhat * (dv22 - 2.0*v22or);
        Tb +=  C_a * pref * 2.0 * rxQbr * v22;

        if (a_is_Fluctuating) dUdCa += pref * (v21 * trQb + v22 * rdQbr);
      }
    }

    if (a_is_Dipole) {

      if (b_is_Charge) {
        pref = pre12_ * *(idat.electroMult);

        U  -= C_b * pref * v11 * rdDa;
        F  -= C_b * pref * ((dv11-v11or) * rdDa * rhat + v11or * D_a);
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
        F  -= pref * (dv21 * DadDb * rhat + v22or * (rdDb * D_a + rdDa * D_b));
        F  -= pref * (rdDa * rdDb) * (dv22 - 2.0*v22or) * rhat;
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
        F  -= pref * (trQb*D_a + 2.0*DadQb) * v31or;
        F  -= pref * (trQb*rdDa + 2.0*DadQbr) * (dv31-v31or) * rhat;
        F  -= pref * (D_a*rdQbr + 2.0*rdDa*rQb) * v32or;
        F  -= pref * (rdDa * rdQbr * rhat * (dv32-3.0*v32or));
        Ta += pref * ((-trQb*rxDa + 2.0 * DaxQbr)*v31 - rxDa*rdQbr*v32);
        Tb += pref * ((2.0*cross(DadQb, rhat) - 2.0*DaxQbr)*v31 
                      - 2.0*rdDa*rxQbr*v32);
      }
    }

    if (a_is_Quadrupole) {
      if (b_is_Charge) {
        pref = pre14_ * *(idat.electroMult);
        U  += C_b * pref * (v21 * trQa + v22 * rdQar);
        F  += C_b * pref * (trQa * rhat * dv21 + 2.0 * Qar * v22or);
        F  += C_b * pref * rdQar * rhat * (dv22 - 2.0*v22or);
        Ta += C_b * pref * 2.0 * rxQar * v22;

        if (b_is_Fluctuating) dUdCb += pref * (v21 * trQa + v22 * rdQar);
      }
      if (b_is_Dipole) {
        pref = pre24_ * *(idat.electroMult);
        DbdQa = D_b * Q_a;
        DbdQar = dot(D_b, Qar);
        DbxQar = cross(D_b, Qar);

        U  += pref * ((trQa*rdDb + 2.0*DbdQar)*v31 + rdDb*rdQar*v32);
        F  += pref * (trQa*D_b + 2.0*DbdQa) * v31or;
        F  += pref * (trQa*rdDb + 2.0*DbdQar) * (dv31-v31or) * rhat;
        F  += pref * (D_b*rdQar + 2.0*rdDb*rQa) * v32or;
        F  += pref * (rdDb * rdQar * rhat * (dv32-3.0*v32or));
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

        F  += pref * rhat * (trQa * trQb + 2.0 * trQaQb)*dv41;
        F  += pref*rhat*(trQa*rdQbr + trQb*rdQar + 4.0*rQaQbr)*(dv42-2.0*v42or);
        F  += pref * rhat * (rdQar * rdQbr)*(dv43 - 4.0*v43or);

        F  += pref * 2.0 * (trQb*rQa + trQa*rQb) * v42or;
        F  += pref * 4.0 * (rQaQb + QaQbr) * v42or;
        F  += pref * 2.0 * (rQa*rdQbr + rdQar*rQb) * v43or;

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

    ElectrostaticAtomData data = ElectrostaticMap[Etids[sdat.atid]];
    
    // logicals
    bool i_is_Charge = data.is_Charge;
    bool i_is_Dipole = data.is_Dipole;
    bool i_is_Quadrupole = data.is_Quadrupole;
    bool i_is_Fluctuating = data.is_Fluctuating;
    RealType C_a = data.fixedCharge;   
    RealType self(0.0), preVal, DdD, trQ, trQQ;

    if (i_is_Dipole) {
      DdD = data.dipole.lengthSquare();
    }
        
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
        (*(sdat.pot))[ELECTROSTATIC_FAMILY] -= pre22_ * preRF_ * DdD;
      }
      
      break;
      
    case esm_SHIFTED_FORCE:
    case esm_SHIFTED_POTENTIAL:
    case esm_TAYLOR_SHIFTED:
      if (i_is_Charge) 
        self += selfMult1_ * pre11_ * C_a * (C_a + *(sdat.skippedCharge));      
      if (i_is_Dipole) 
        self += selfMult2_ * pre22_ * DdD;      
      if (i_is_Quadrupole) {
        trQ = data.quadrupole.trace();
        trQQ = (data.quadrupole * data.quadrupole).trace();
        self += selfMult4_ * pre44_ * (2.0*trQQ + trQ*trQ);
        if (i_is_Charge)
          self -= selfMult2_ * pre14_ * 2.0 * C_a * trQ;
      }
      (*(sdat.pot))[ELECTROSTATIC_FAMILY] += self;      
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
