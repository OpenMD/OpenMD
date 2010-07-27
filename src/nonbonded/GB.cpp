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

#include <stdio.h>
#include <string.h>

#include <cmath>
#include "nonbonded/GB.hpp"
#include "nonbonded/LJ.hpp"
#include "utils/simError.h"

using namespace std;
namespace OpenMD {

  bool GB::initialized_ = false;
  RealType GB::mu_ = 2.0;
  RealType GB::nu_ = 1.0;
  ForceField* GB::forceField_ = NULL;
  map<int, AtomType*> GB::GBMap;
  map<pair<AtomType*, AtomType*>, GBInteractionData> GB::MixingMap;
  
  GB* GB::_instance = NULL;

  GB* GB::Instance() {
    if (!_instance) {
      _instance = new GB();
    }
    return _instance;
  }

  GayBerneParam GB::getGayBerneParam(AtomType* atomType) {
    
    // Do sanity checking on the AtomType we were passed before
    // building any data structures:
    if (!atomType->isGayBerne()) {
      sprintf( painCave.errMsg,
               "GB::getGayBerneParam was passed an atomType (%s) that does\n"
               "\tnot appear to be a Gay-Berne atom.\n",
               atomType->getName().c_str());
      painCave.severity = OPENMD_ERROR;
      painCave.isFatal = 1;
      simError();
    }
    
    DirectionalAtomType* daType = dynamic_cast<DirectionalAtomType*>(atomType);
    GenericData* data = daType->getPropertyByName("GayBerne");
    if (data == NULL) {
      sprintf( painCave.errMsg, "GB::getGayBerneParam could not find\n"
               "\tGay-Berne parameters for atomType %s.\n", 
               daType->getName().c_str());
      painCave.severity = OPENMD_ERROR;
      painCave.isFatal = 1;
      simError(); 
    }
    
    GayBerneParamGenericData* gbData = dynamic_cast<GayBerneParamGenericData*>(data);
    if (gbData == NULL) {
      sprintf( painCave.errMsg,
               "GB::getGayBerneParam could not convert GenericData to\n"
               "\tGayBerneParamGenericData for atom type %s\n", 
               daType->getName().c_str());
      painCave.severity = OPENMD_ERROR;
      painCave.isFatal = 1;
      simError();          
    }
    
    return gbData->getData();
  }

  void GB::initialize() {    
    ForceField::AtomTypeContainer* atomTypes = forceField_->getAtomTypes();
    ForceField::AtomTypeContainer::MapTypeIterator i;
    AtomType* at;

    // GB handles all of the GB-GB interactions as well as GB-LJ cross
    // interactions:

    for (at = atomTypes->beginType(i); at != NULL;
         at = atomTypes->nextType(i)) {
      
      if (at->isGayBerne() || at->isLennardJones())
        addType(at);
    }
   
    initialized_ = true;
  }
      
  void GB::addType(AtomType* atomType){
    // add it to the map:
    AtomTypeProperties atp = atomType->getATP();    

    pair<map<int,AtomType*>::iterator,bool> ret;    
    ret = GBMap.insert( pair<int, AtomType*>(atp.ident, atomType) );
    if (ret.second == false) {
      sprintf( painCave.errMsg,
               "GB already had a previous entry with ident %d\n",
               atp.ident);
      painCave.severity = OPENMD_INFO;
      painCave.isFatal = 0;
      simError();         
    }
    
    RealType d1, l1, e1, er1, dw1;
    
    if (atomType->isGayBerne()) {
      GayBerneParam gb1 = getGayBerneParam(atomType);
      d1 = gb1.GB_d;
      l1 = gb1.GB_l;
      e1 = gb1.GB_eps;
      er1 = gb1.GB_eps_ratio;
      dw1 = gb1.GB_dw;
    } else if (atomType->isLennardJones()) {
      d1 = LJ::Instance()->getSigma(atomType) / sqrt(2.0);
      e1 = LJ::Instance()->getEpsilon(atomType);
      l1 = d1;
      er1 = 1.0;
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
      
      RealType d2, l2, e2, er2, dw2;
      
      if (atype2->isGayBerne()) {
        GayBerneParam gb2 = getGayBerneParam(atype2);
        d2 = gb2.GB_d;
        l2 = gb2.GB_l;
        e2 = gb2.GB_eps;
        er2 = gb2.GB_eps_ratio;
        dw2 = gb2.GB_dw;
      } else if (atype2->isLennardJones()) {
        d2 = LJ::Instance()->getSigma(atype2) / sqrt(2.0);
        e2 = LJ::Instance()->getEpsilon(atype2);
        l2 = d2;
        er2 = 1.0;
        dw2 = 1.0;
      } 
                       
      GBInteractionData mixer;        
      
      //  Cleaver paper uses sqrt of squares to get sigma0 for
      //  mixed interactions.
            
      mixer.sigma0 = sqrt(d1*d1 + d2*d2);
      mixer.xa2 = (l1*l1 - d1*d1)/(l1*l1 + d2*d2);
      mixer.xai2 = (l2*l2 - d2*d2)/(l2*l2 + d1*d1);
      mixer.x2 = (l1*l1 - d1*d1) * (l2*l2 - d2*d2) /
        ((l2*l2 + d1*d1) * (l1*l1 + d2*d2));
 
      // assumed LB mixing rules for now:

      mixer.dw = 0.5 * (dw1 + dw2);
      mixer.eps0 = sqrt(e1 * e2);
      
      RealType er = sqrt(er1 * er2);
      RealType ermu = pow(er,(1.0 / mu_));
      RealType xp = (1.0 - ermu) / (1.0 + ermu);
      RealType ap2 = 1.0 / (1.0 + ermu);
      
      mixer.xp2 = xp * xp;
      mixer.xpap2 = xp * ap2;
      mixer.xpapi2 = xp / ap2;

      // only add this pairing if at least one of the atoms is a Gay-Berne atom

      if (atomType->isGayBerne() || atype2->isGayBerne()) {

        pair<AtomType*, AtomType*> key1, key2;
        key1 = make_pair(atomType, atype2);
        key2 = make_pair(atype2, atomType);
        
        MixingMap[key1] = mixer;
        if (key2 != key1) {
          MixingMap[key2] = mixer;
        }
      }
    }      
  }
  

  RealType GB::getGayBerneCut(int atid) { 
    if (!initialized_) initialize();
    std::map<int, AtomType*> :: const_iterator it;
    it = GBMap.find(atid);
    if (it == GBMap.end()) {
      sprintf( painCave.errMsg,
               "GB::getGayBerneCut could not find atid %d in GBMap\n",
               (atid));
      painCave.severity = OPENMD_ERROR;
      painCave.isFatal = 1;
      simError();          
    }

    AtomType* atype = it->second;

    RealType gbCut;
    
    if (atype->isGayBerne()) {
      GayBerneParam gb = getGayBerneParam(atype);

      // sigma is actually sqrt(2) * l for prolate ellipsoids
      gbCut = 2.5 * sqrt(2.0) * max(gb.GB_l, gb.GB_d);

    } else if (atype->isLennardJones()) {
      gbCut = 2.5 * LJ::Instance()->getSigma(atype);
    }
    
    return gbCut;
  }

 
  void GB::calcForce(AtomType* at1, AtomType* at2, Vector3d d, 
                     RealType r, RealType r2, RealType sw, 
                     RealType vdwMult, RealType &vpair, RealType &pot, 
                     RotMat3x3d A1, RotMat3x3d A2, Vector3d &f1,
                     Vector3d &t1, Vector3d &t2) {

    if (!initialized_) initialize();
    
    pair<AtomType*, AtomType*> key = make_pair(at1, at2);
    GBInteractionData mixer = MixingMap[key];

    RealType sigma0 = mixer.sigma0;
    RealType dw     = mixer.dw;
    RealType eps0   = mixer.eps0;  
    RealType x2     = mixer.x2;    
    RealType xa2    = mixer.xa2;   
    RealType xai2   = mixer.xai2;  
    RealType xp2    = mixer.xp2;   
    RealType xpap2  = mixer.xpap2; 
    RealType xpapi2 = mixer.xpapi2;

    Vector3d ul1 = A1.getColumn(2);
    Vector3d ul2 = A2.getColumn(2);

    RealType a, b, g;

    bool i_is_LJ = at1->isLennardJones();
    bool j_is_LJ = at2->isLennardJones();

    if (i_is_LJ) {
      a = 0.0;
      ul1 = V3Zero;
    } else {
      a = dot(d, ul1);
    }

    if (j_is_LJ) {
      b = 0.0;
      ul2 = V3Zero;
    } else {
      b = dot(d, ul2);
    }

    if (i_is_LJ || j_is_LJ) 
      g = 0.0;
    else
      g = dot(ul1, ul2);

    RealType au = a / r;
    RealType bu = b / r;
    
    RealType au2 = au * au;
    RealType bu2 = bu * bu;
    RealType g2 = g * g;
    
    RealType H  = (xa2 * au2 + xai2 * bu2 - 2.0*x2*au*bu*g)  / (1.0 - x2*g2);
    RealType Hp = (xpap2*au2 + xpapi2*bu2 - 2.0*xp2*au*bu*g) / (1.0 - xp2*g2);

    RealType sigma = sigma0 / sqrt(1.0 - H);
    RealType e1 = 1.0 / sqrt(1.0 - x2*g2);
    RealType e2 = 1.0 - Hp;
    RealType eps = eps0 * pow(e1,nu_) * pow(e2,mu_);
    RealType BigR = dw*sigma0 / (r - sigma + dw*sigma0);
    
    RealType R3 = BigR*BigR*BigR;
    RealType R6 = R3*R3;
    RealType R7 = R6 * BigR;
    RealType R12 = R6*R6;
    RealType R13 = R6*R7;

    RealType U = vdwMult * 4.0 * eps * (R12 - R6);

    RealType s3 = sigma*sigma*sigma;
    RealType s03 = sigma0*sigma0*sigma0;

    RealType pref1 = - vdwMult * 8.0 * eps * mu_ * (R12 - R6) / (e2 * r);

    RealType pref2 = vdwMult * 8.0 * eps * s3 * (6.0*R13 - 3.0*R7) /(dw*r*s03);

    RealType dUdr = - (pref1 * Hp + pref2 * (sigma0*sigma0*r/s3 + H));
    
    RealType dUda = pref1 * (xpap2*au - xp2*bu*g) / (1.0 - xp2 * g2) 
      + pref2 * (xa2 * au - x2 *bu*g) / (1.0 - x2 * g2);
    
    RealType dUdb = pref1 * (xpapi2*bu - xp2*au*g) / (1.0 - xp2 * g2) 
      + pref2 * (xai2 * bu - x2 *au*g) / (1.0 - x2 * g2);

    RealType dUdg = 4.0 * eps * nu_ * (R12 - R6) * x2 * g / (1.0 - x2*g2)
      + 8.0 * eps * mu_ * (R12 - R6) * (xp2*au*bu - Hp*xp2*g) / 
      (1.0 - xp2 * g2) / e2 + 8.0 * eps * s3 * (3.0 * R7 - 6.0 * R13) * 
      (x2 * au * bu - H * x2 * g) / (1.0 - x2 * g2) / (dw * s03);
    

    Vector3d rhat = d / r;   
    Vector3d rxu1 = cross(d, ul1);
    Vector3d rxu2 = cross(d, ul2);
    Vector3d uxu = cross(ul1, ul2);
    
    pot += U*sw;
    f1 += dUdr * rhat + dUda * ul1 + dUdb * ul2;    
    t1 += dUda * rxu1 - dUdg * uxu;
    t2 += dUdb * rxu2 - dUdg * uxu;
    vpair += U*sw;

    return;

  }

  void GB::do_gb_pair(int *atid1, int *atid2, RealType *d, RealType *r, 
                      RealType *r2, RealType *sw, RealType *vdwMult,
                      RealType *vpair, RealType *pot, RealType *A1,
                      RealType *A2, RealType *f1, RealType *t1, RealType *t2) {
    
    if (!initialized_) initialize();
    
    AtomType* atype1 = GBMap[*atid1];
    AtomType* atype2 = GBMap[*atid2];
    
    Vector3d disp(d);
    Vector3d frc(f1);
    Vector3d trq1(t1);
    Vector3d trq2(t2);
    RotMat3x3d Ai(A1);
    RotMat3x3d Aj(A2);
   
    // Fortran has the opposite matrix ordering from c++, so we'll use 
    // transpose here.  When we finish the conversion to C++, this wrapper 
    // will disappear, as will the transpose below:

    calcForce(atype1, atype2, disp, *r, *r2, *sw, *vdwMult, *vpair, *pot, 
              Ai.transpose(), Aj.transpose(), frc, trq1, trq1);
      
    f1[0] = frc.x();
    f1[1] = frc.y();
    f1[2] = frc.z();

    t1[0] = trq1.x();
    t1[1] = trq1.y();
    t1[2] = trq1.z();

    t2[0] = trq2.x();
    t2[1] = trq2.y();
    t2[2] = trq2.z();

    return;    
  }
}

extern "C" {
  
#define fortranGetGayBerneCut FC_FUNC(getgaybernecut, GETGAYBERNECUT)
#define fortranDoGBPair FC_FUNC(do_gb_pair, DO_GB_PAIR)
  
  RealType fortranGetGayBerneCut(int* atid) {
    return OpenMD::GB::Instance()->getGayBerneCut(*atid);
  }

  void fortranDoGBPair(int *atid1, int *atid2, RealType *d, RealType *r, 
                       RealType *r2, RealType *sw, RealType *vdwMult, 
                       RealType *vpair, RealType *pot, RealType *A1, 
                       RealType *A2, RealType *f1, RealType *t1, RealType *t2){
    
    return OpenMD::GB::Instance()->do_gb_pair(atid1, atid2, d, r, r2, sw,
                                              vdwMult, vpair, pot, A1, A2, f1,
                                              t1, t2);
  }
}
