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
 * [3]  Sun, Lin & Gezelter, J. Chem. Phys. 128, 234107 (2008).          
 * [4]  Kuang & Gezelter,  J. Chem. Phys. 133, 164101 (2010).
 * [5]  Vardeman, Stocker & Gezelter, J. Chem. Theory Comput. 7, 834 (2011).
 */

#include <stdio.h>
#include <string.h>

#include <cmath>
#include "nonbonded/GB.hpp"
#include "utils/simError.h"
#include "types/LennardJonesAdapter.hpp"
#include "types/GayBerneAdapter.hpp"

using namespace std;
namespace OpenMD {

  /* GB is the Gay-Berne interaction for ellipsoidal particles.  The original 
   * paper (for identical uniaxial particles) is:
   *    J. G. Gay and B. J. Berne, J. Chem. Phys., 74, 3316-3319 (1981).
   * A more-general GB potential for dissimilar uniaxial particles:
   *    D. J. Cleaver, C. M. Care, M. P. Allen and M. P. Neal, Phys. Rev. E, 
   *    54, 559-567 (1996).
   * Further parameterizations can be found in:
   *    A. P. J. Emerson, G. R. Luckhurst and S. G. Whatling, Mol. Phys., 
   *    82, 113-124 (1994).
   * And a nice force expression:
   *    G. R. Luckhurst and R. A. Stephens, Liq. Cryst. 8, 451-464 (1990).
   * Even clearer force and torque expressions:
   *    P. A. Golubkov and P. Y. Ren, J. Chem. Phys., 125, 64103 (2006).
   * New expressions for cross interactions of strength parameters:
   *    J. Wu, X. Zhen, H. Shen, G. Li, and P. Ren, J. Chem. Phys., 
   *    135, 155104 (2011).
   *
   * In this version of the GB interaction, each uniaxial ellipsoidal type 
   * is described using a set of 6 parameters:
   *  d:  range parameter for side-by-side (S) and cross (X) configurations
   *  l:  range parameter for end-to-end (E) configuration
   *  epsilon_X:  well-depth parameter for cross (X) configuration
   *  epsilon_S:  well-depth parameter for side-by-side (S) configuration
   *  epsilon_E:  well depth parameter for end-to-end (E) configuration
   *  dw: "softness" of the potential
   * 
   * Additionally, there are two "universal" paramters to govern the overall 
   * importance of the purely orientational (nu) and the mixed
   * orientational / translational (mu) parts of strength of the interactions. 
   * These parameters have default or "canonical" values, but may be changed
   * as a force field option:
   * nu_: purely orientational part : defaults to 1
   * mu_: mixed orientational / translational part : defaults to 2
   */


  GB::GB() : name_("GB"), initialized_(false), mu_(2.0), nu_(1.0), forceField_(NULL) {}
    
  void GB::initialize() {    
    
    ForceFieldOptions& fopts = forceField_->getForceFieldOptions();
    mu_ = fopts.getGayBerneMu();
    nu_ = fopts.getGayBerneNu();
    ForceField::AtomTypeContainer* atomTypes = forceField_->getAtomTypes();
    ForceField::AtomTypeContainer::MapTypeIterator i;
    AtomType* at;

    // GB handles all of the GB-GB interactions as well as GB-LJ cross
    // interactions:

    for (at = atomTypes->beginType(i); at != NULL;
         at = atomTypes->nextType(i)) {
      
      LennardJonesAdapter lja = LennardJonesAdapter(at);
      GayBerneAdapter gba = GayBerneAdapter(at);

      if (gba.isGayBerne() || lja.isLennardJones())
        addType(at);
    }
   
    initialized_ = true;
  }
      
  void GB::addType(AtomType* atomType){
    // add it to the map:

    pair<map<int,AtomType*>::iterator,bool> ret;    
    ret = GBMap.insert( pair<int, AtomType*>(atomType->getIdent(), atomType) );
    if (ret.second == false) {
      sprintf( painCave.errMsg,
               "GB already had a previous entry with ident %d\n",
               atomType->getIdent() );
      painCave.severity = OPENMD_INFO;
      painCave.isFatal = 0;
      simError();         
    }
    
    RealType d1, l1, eX1, eS1, eE1, dw1;
        
    LennardJonesAdapter lja1 = LennardJonesAdapter(atomType);
    GayBerneAdapter gba1 = GayBerneAdapter(atomType);
    if (gba1.isGayBerne()) {
      d1 = gba1.getD();
      l1 = gba1.getL();
      eX1 = gba1.getEpsX();
      eS1 = gba1.getEpsS();
      eE1 = gba1.getEpsE();
      dw1 = gba1.getDw();
    } else if (lja1.isLennardJones()) {
      d1 = lja1.getSigma() / sqrt(2.0);
      l1 = d1;
      eX1 = lja1.getEpsilon();
      eS1 = eX1;
      eE1 = eX1;
      dw1 = 1.0;      
    } else {
      sprintf( painCave.errMsg,
               "GB::addType was passed an atomType (%s) that does not\n"
               "\tappear to be a Gay-Berne or Lennard-Jones atom.\n",
               atomType->getName().c_str());
      painCave.severity = OPENMD_ERROR;
      painCave.isFatal = 1;
      simError();
    }
      

    // Now, iterate over all known types and add to the mixing map:
    
    map<int, AtomType*>::iterator it;
    for( it = GBMap.begin(); it != GBMap.end(); ++it) {
      
      AtomType* atype2 = (*it).second;
      LennardJonesAdapter lja2 = LennardJonesAdapter(atype2);
      GayBerneAdapter gba2 = GayBerneAdapter(atype2);
      RealType d2, l2, eX2, eS2, eE2, dw2;
      
      if (gba2.isGayBerne()) {
        d2 = gba2.getD();
        l2 = gba2.getL();
        eX2 = gba2.getEpsX();
        eS2 = gba2.getEpsS();
        eE2 = gba2.getEpsE();
        dw2 = gba2.getDw();
      } else if (lja2.isLennardJones()) {
        d2 = lja2.getSigma() / sqrt(2.0);
        l2 = d2;
        eX2 = lja2.getEpsilon();
        eS2 = eX2;
        eE2 = eX2;
        dw2 = 1.0;
      } else {
        sprintf( painCave.errMsg,
                 "GB::addType found an atomType (%s) that does not\n"
                 "\tappear to be a Gay-Berne or Lennard-Jones atom.\n",
                 atype2->getName().c_str());
        painCave.severity = OPENMD_ERROR;
        painCave.isFatal = 1;
        simError();
      }

                       
      GBInteractionData mixer1, mixer2;     
      
      //  Cleaver paper uses sqrt of squares to get sigma0 for
      //  mixed interactions.
            
      mixer1.sigma0 = sqrt(d1*d1 + d2*d2);
      mixer1.xa2 = (l1*l1 - d1*d1)/(l1*l1 + d2*d2);
      mixer1.xai2 = (l2*l2 - d2*d2)/(l2*l2 + d1*d1);
      mixer1.x2 = (l1*l1 - d1*d1) * (l2*l2 - d2*d2) /
        ((l2*l2 + d1*d1) * (l1*l1 + d2*d2));

      mixer2.sigma0 = mixer1.sigma0;
      // xa2 and xai2 for j-i pairs are reversed from the same i-j pairing.
      // Swapping the particles reverses the anisotropy parameters:
      mixer2.xa2 = mixer1.xai2;
      mixer2.xai2 = mixer1.xa2;
      mixer2.x2 = mixer1.x2;
 
      // assumed LB mixing rules for now:

      mixer1.dw = 0.5 * (dw1 + dw2);
      mixer1.eps0 = sqrt(eX1 * eX2);

      mixer2.dw = mixer1.dw;
      mixer2.eps0 = mixer1.eps0;

      RealType mi = RealType(1.0)/mu_;
      
      mixer1.xpap2  = (pow(eS1, mi) - pow(eE1, mi)) / (pow(eS1, mi) + pow(eE2, mi));
      mixer1.xpapi2 = (pow(eS2, mi) - pow(eE2, mi)) / (pow(eS2, mi) + pow(eE1, mi));
      mixer1.xp2    = (pow(eS1, mi) - pow(eE1, mi)) * (pow(eS2, mi) - pow(eE2, mi))  / 
        (pow(eS2, mi) + pow(eE1, mi)) / (pow(eS1, mi) + pow(eE2, mi)) ;

      // xpap2 and xpapi2 for j-i pairs are reversed from the same i-j pairing.
      // Swapping the particles reverses the anisotropy parameters:
      mixer2.xpap2 = mixer1.xpapi2;
      mixer2.xpapi2 = mixer1.xpap2;
      mixer2.xp2 = mixer1.xp2;

      // only add this pairing if at least one of the atoms is a Gay-Berne atom

      if (gba1.isGayBerne() || gba2.isGayBerne()) {

        pair<AtomType*, AtomType*> key1, key2;
        key1 = make_pair(atomType, atype2);
        key2 = make_pair(atype2, atomType);
        
        MixingMap[key1] = mixer1;
        if (key2 != key1) {
          MixingMap[key2] = mixer2;
        }
      }
    }      
  }
   
  void GB::calcForce(InteractionData &idat) {

    if (!initialized_) initialize();
    
    GBInteractionData mixer = MixingMap[idat.atypes];

    RealType sigma0 = mixer.sigma0;
    RealType dw     = mixer.dw;
    RealType eps0   = mixer.eps0;  
    RealType x2     = mixer.x2;    
    RealType xa2    = mixer.xa2;   
    RealType xai2   = mixer.xai2;  
    RealType xp2    = mixer.xp2;   
    RealType xpap2  = mixer.xpap2; 
    RealType xpapi2 = mixer.xpapi2;

    // cerr << "atypes = " << idat.atypes.first->getName() << " " << idat.atypes.second->getName() << "\n";
    // cerr << "sigma0 = " <<mixer.sigma0 <<"\n";
    // cerr << "dw     = " <<mixer.dw <<"\n";
    // cerr << "eps0   = " <<mixer.eps0 <<"\n";  
    // cerr << "x2     = " <<mixer.x2 <<"\n";    
    // cerr << "xa2    = " <<mixer.xa2 <<"\n";   
    // cerr << "xai2   = " <<mixer.xai2 <<"\n";  
    // cerr << "xp2    = " <<mixer.xp2 <<"\n";   
    // cerr << "xpap2  = " <<mixer.xpap2 <<"\n"; 
    // cerr << "xpapi2 = " <<mixer.xpapi2 <<"\n";

    Vector3d ul1 = idat.A1->getRow(2);
    Vector3d ul2 = idat.A2->getRow(2);

    // cerr << "ul1 = " <<ul1<<"\n";
    // cerr << "ul2 = " <<ul2<<"\n";

    RealType a, b, g;

    // This is not good.  We should store this in the mixing map, and not
    // query atom types in calc force.
    bool i_is_LJ = idat.atypes.first->isLennardJones();
    bool j_is_LJ = idat.atypes.second->isLennardJones();
    
    if (i_is_LJ) {
      a = 0.0;
      ul1 = V3Zero;
    } else {
      a = dot(*(idat.d), ul1);
    }

    if (j_is_LJ) {
      b = 0.0;
      ul2 = V3Zero;
    } else {
      b = dot(*(idat.d), ul2);
    }

    if (i_is_LJ || j_is_LJ) 
      g = 0.0;
    else
      g = dot(ul1, ul2);

    RealType au = a / *(idat.rij);
    RealType bu = b / *(idat.rij);
    
    RealType au2 = au * au;
    RealType bu2 = bu * bu;
    RealType g2 = g * g;

    RealType H  = (xa2 * au2 + xai2 * bu2 - 2.0*x2*au*bu*g)  / (1.0 - x2*g2);
    RealType Hp = (xpap2*au2 + xpapi2*bu2 - 2.0*xp2*au*bu*g) / (1.0 - xp2*g2);

    // cerr << "au2 = " << au2 << "\n";
    // cerr << "bu2 = " << bu2 << "\n";
    // cerr << "g2 = " << g2 << "\n";
    // cerr << "H = " << H << "\n";
    // cerr << "Hp = " << Hp << "\n";

    RealType sigma = sigma0 / sqrt(1.0 - H);
    RealType e1 = 1.0 / sqrt(1.0 - x2*g2);
    RealType e2 = 1.0 - Hp;
    RealType eps = eps0 * pow(e1,nu_) * pow(e2,mu_);
    RealType BigR = dw*sigma0 / (*(idat.rij) - sigma + dw*sigma0);
    
    RealType R3 = BigR*BigR*BigR;
    RealType R6 = R3*R3;
    RealType R7 = R6 * BigR;
    RealType R12 = R6*R6;
    RealType R13 = R6*R7;

    RealType U = *(idat.vdwMult) * 4.0 * eps * (R12 - R6);

    RealType s3 = sigma*sigma*sigma;
    RealType s03 = sigma0*sigma0*sigma0;

    // cerr << "vdwMult = " << *(idat.vdwMult) << "\n";
    // cerr << "eps = " << eps <<"\n";
    // cerr << "mu = " << mu_ << "\n";
    // cerr << "R12 = " << R12 << "\n";
    // cerr << "R6 = " << R6 << "\n";
    // cerr << "R13 = " << R13 << "\n";
    // cerr << "R7 = " << R7 << "\n";
    // cerr << "e2 = " << e2 << "\n";
    // cerr << "rij = " << *(idat.rij) << "\n";
    // cerr << "s3 = " << s3 << "\n";
    // cerr << "s03 = " << s03 << "\n";
    // cerr << "dw = " << dw << "\n";

    RealType pref1 = - *(idat.vdwMult) * 8.0 * eps * mu_ * (R12 - R6) / 
      (e2 * *(idat.rij));

    RealType pref2 = *(idat.vdwMult) * 8.0 * eps * s3 * (6.0*R13 - 3.0*R7) /
      (dw*  *(idat.rij) * s03);

    RealType dUdr = - (pref1 * Hp + pref2 * (sigma0 * sigma0 *  
                                             *(idat.rij) / s3 + H));
    
    RealType dUda = pref1 * (xpap2*au - xp2*bu*g) / (1.0 - xp2 * g2) 
      + pref2 * (xa2 * au - x2 *bu*g) / (1.0 - x2 * g2);
    
    RealType dUdb = pref1 * (xpapi2*bu - xp2*au*g) / (1.0 - xp2 * g2) 
      + pref2 * (xai2 * bu - x2 *au*g) / (1.0 - x2 * g2);

    RealType dUdg = 4.0 * eps * nu_ * (R12 - R6) * x2 * g / (1.0 - x2*g2)
      + 8.0 * eps * mu_ * (R12 - R6) * (xp2*au*bu - Hp*xp2*g) / 
      (1.0 - xp2 * g2) / e2 + 8.0 * eps * s3 * (3.0 * R7 - 6.0 * R13) * 
      (x2 * au * bu - H * x2 * g) / (1.0 - x2 * g2) / (dw * s03);

    // cerr << "pref = " << pref1 << " " << pref2 << "\n";
    // cerr << "dU = " << dUdr << " " << dUda <<" " << dUdb << " " << dUdg << "\n";

    Vector3d rhat = *(idat.d) / *(idat.rij);   
    Vector3d rxu1 = cross(*(idat.d), ul1);
    Vector3d rxu2 = cross(*(idat.d), ul2);
    Vector3d uxu = cross(ul1, ul2);

    (*(idat.pot))[VANDERWAALS_FAMILY] += U *  *(idat.sw);
    *(idat.f1) += (dUdr * rhat + dUda * ul1 + dUdb * ul2) * *(idat.sw);
    *(idat.t1) += (dUda * rxu1 - dUdg * uxu) * *(idat.sw);
    *(idat.t2) += (dUdb * rxu2 + dUdg * uxu) * *(idat.sw);
    *(idat.vpair) += U;

    // cerr << "f1 term = " <<  (dUdr * rhat + dUda * ul1 + dUdb * ul2) * *(idat.sw) << "\n";
    // cerr << "t1 term = " << (dUda * rxu1 - dUdg * uxu) * *(idat.sw) << "\n";
    // cerr << "t2 term = " << (dUdb * rxu2 + dUdg * uxu) * *(idat.sw) << "\n";
    // cerr << "vp term = " << U << "\n";

    return;

  }

  RealType GB::getSuggestedCutoffRadius(pair<AtomType*, AtomType*> atypes) {
    if (!initialized_) initialize();   

    RealType cut = 0.0;

    LennardJonesAdapter lja1 = LennardJonesAdapter(atypes.first);
    GayBerneAdapter gba1 = GayBerneAdapter(atypes.first);
    LennardJonesAdapter lja2 = LennardJonesAdapter(atypes.second);
    GayBerneAdapter gba2 = GayBerneAdapter(atypes.second);

    if (gba1.isGayBerne()) {
      RealType d1 = gba1.getD();
      RealType l1 = gba1.getL();
      // sigma is actually sqrt(2)*l  for prolate ellipsoids 
      cut = max(cut, RealType(2.5) * sqrt(RealType(2.0)) * max(d1, l1));
    } else if (lja1.isLennardJones()) {
      cut = max(cut, RealType(2.5) * lja1.getSigma());
    }

    if (gba2.isGayBerne()) {
      RealType d2 = gba2.getD();
      RealType l2 = gba2.getL();
      cut = max(cut, RealType(2.5) * sqrt(RealType(2.0)) * max(d2, l2));
    } else if (lja2.isLennardJones()) {
      cut = max(cut, RealType(2.5) * lja2.getSigma());
    }
   
    return cut;
  }
}

