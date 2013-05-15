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
#include "nonbonded/SHAPES.hpp"
#include "nonbonded/LJ.hpp"
#include "utils/simError.h"

using namespace std;
namespace OpenMD {
 
  SHAPES::SHAPES() {
    initialized_ = false;
    lMax_ = 64;
    mMax_ = 64;
    forceField_ = NULL;
  }
  
  void SHAPES::initialize() {    
    
    ForceFieldOptions& fopts = forceField_->getForceFieldOptions();
    ForceField::AtomTypeContainer* atomTypes = forceField_->getAtomTypes();
    ForceField::AtomTypeContainer::MapTypeIterator i;
    AtomType* at;

    // SHAPES handles all of the SHAPES-SHAPES interactions as well as
    // SHAPES-LJ cross interactions:
    
    for (at = atomTypes->beginType(i); at != NULL;
         at = atomTypes->nextType(i)) {
      
      if (at->isShape()) 
        addShape(dynamic_cast<ShapeAtomType*>(at));

      if (at->isLennardJones())
        addLJ(at);
      
    }
    
    initialized_ = true;
  }
  
  void SHAPES::addShape(ShapeAtomType* atomType){
    // add it to the map:
    AtomTypeProperties atp = atomType->getATP();    

    if (atomType->isShape() ) {
      pair<map<int,ShapeAtomType*>::iterator, bool> ret;
      ret = ShapesMap.insert( pair<int, ShapeAtomType*>(atp.ident, atomType));
      if (ret.second == false) {
        sprintf( painCave.errMsg,
                 "SHAPES already had a previous entry with ident %d\n",
                 atp.ident);
        painCave.severity = OPENMD_INFO;
        painCave.isFatal = 0;
        simError();         
      }   
      
      ShapesMap.insert( pair<int, ShapeAtomType*>(atp.ident, static_cast<ShapeAtomType*>(atomType)) );
      
    } else if (atomType->isLennardJones()) {
      RealType d1 = getLJSigma(atomType) / sqrt(2.0);
      RealType e1 = getLJEpsilon(atomType);
    } else {
      sprintf( painCave.errMsg,
               "SHAPES::addType was passed an atomType (%s) that does not\n"
               "\tappear to be a SHAPES or Lennard-Jones atom.\n",
               atomType->getName().c_str());
      painCave.severity = OPENMD_ERROR;
      painCave.isFatal = 1;
      simError();
    }   
  }
  

  LJParam SHAPES::getLJParam(AtomType* atomType) {
    
    // Do sanity checking on the AtomType we were passed before
    // building any data structures:
    if (!atomType->isLennardJones()) {
      sprintf( painCave.errMsg,
               "SHAPES::getLJParam was passed an atomType (%s) that does not\n"
               "\tappear to be a Lennard-Jones atom.\n",
               atomType->getName().c_str());
      painCave.severity = OPENMD_ERROR;
      painCave.isFatal = 1;
      simError();
    }
    
    GenericData* data = atomType->getPropertyByName("LennardJones");
    if (data == NULL) {
      sprintf( painCave.errMsg, "SHAPES::getLJParam could not find Lennard-Jones\n"
               "\tparameters for atomType %s.\n", atomType->getName().c_str());
      painCave.severity = OPENMD_ERROR;
      painCave.isFatal = 1;
      simError(); 
    }
    
    LJParamGenericData* ljData = dynamic_cast<LJParamGenericData*>(data);
    if (ljData == NULL) {
      sprintf( painCave.errMsg,
               "SHAPES::getLJParam could not convert GenericData to LJParam for\n"
               "\tatom type %s\n", atomType->getName().c_str());
      painCave.severity = OPENMD_ERROR;
      painCave.isFatal = 1;
      simError();          
    }
    
    return ljData->getData();
  }
  
  RealType SHAPES::getLJEpsilon(AtomType* atomType) {    
    LJParam ljParam = getLJParam(atomType);
    return ljParam.epsilon;
  }
  RealType SHAPES::getLJSigma(AtomType* atomType) {    
    LJParam ljParam = getLJParam(atomType);
    return ljParam.sigma;
  }

  RealType SHAPES::getGayBerneCut(int atid) { 
    if (!initialized_) initialize();
    std::map<int, AtomType*> :: const_iterator it;
    it = SHAPESMap.find(atid);
    if (it == SHAPESMap.end()) {
      sprintf( painCave.errMsg,
               "SHAPES::getGayBerneCut could not find atid %d in SHAPESMap\n",
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
      gbCut = 2.5 * sqrt(2.0) * max(gb.SHAPES_l, gb.SHAPES_d);

    } else if (atype->isLennardJones()) {
      gbCut = 2.5 * LJ::Instance()->getSigma(atype);
    }
    
    return gbCut;
  }

 
  void SHAPES::calcForce(AtomType* at1, AtomType* at2, Vector3d d, 
                         RealType r, RealType r2, RealType sw, 
                         RealType &vpair, RealType &pot, 
                         RotMat3x3d A1, RotMat3x3d A2, Vector3d &f1,
                         Vector3d &t1, Vector3d &t2) {

    if (!initialized_) initialize();
    
    pair<AtomType*, AtomType*> key = make_pair(at1, at2);
    SHAPESInteractionData mixer = MixingMap[key];

    RealType r3 = r2 * r;
    RealType r5 = r3 * r2;
    
    Vector3d drdi  = -d / r;
    Vector3d drdui = V3Zero;
    Vector3d drdj = d / r;
    Vector3d drduj = V3Zero;
    
    bool i_is_LJ = at1->isLennardJones();
    bool j_is_LJ = at2->isLennardJones();

    RealType sigma_i;
    RealType s_i;
    RealType eps_i;
    Vector3d dsigmaidr;
    Vector3d disgmaidu;
    Vector3d dsidr;
    Vector3d dsidu;
    Vector3d depsidr;
    Vector3d depsidu;

    if (i_is_LJ) {
      sigma_i = LJ::Instance()->getSigma(at1);
      s_i = sigma_i;
      epsilon_i = LJ::Instance()->getEpsilon(at1);
      dsigmaidr = V3Zero;
      dsigmaidu = V3Zero;
      dsidr = V3Zero;
      dsidu = V3Zero;
      depsidr = V3Zero;
      depsidu = V3Zero;
    } else {
      
      // rotate the inter-particle separation into the two different
      // body-fixed coordinate systems:

      Vector3d ri = A1 * d;
      
      RealType xi = ri.x() / r;
      RealType yi = ri.y() / r;
      RealType zi = ri.z() / r;
      RealType xi2 = xi * xi;
      RealType yi2 = yi * yi;
      RealType zi2 = zi * zi;
      RealType cti = zi / r;

      if (cti > 1.0) cti = 1.0;
      if (cti < -1.0_dp) cti = -1.0;

      Vector3d dctidr(-zi * xi / r3, 
                      -zi * yi / r3, 
                      1.0 / r - zi2 / r3);

      Vector3d dctidu(yi / r, 
                      -zi / r,
                      0.0);

      // this is an attempt to try to truncate the singularity when
      // sin(theta) is near 0.0:

      RealType sti2 = 1.0 - cti*cti;
      if (fabs(sti2) < 1.0e-12) {
        RealType proji = sqrt(r * 1.0e-12);
        Vector3d dcpidx(1.0 / proji,
                        0.0,

                        // pickup the ball here!
                        
        dcpidx = 1.0_dp / proji
          dcpidy = 0.0_dp
          dcpidux = xi / proji
          dcpiduy = 0.0_dp
          dspidx = 0.0_dp
          dspidy = 1.0_dp / proji
          dspidux = 0.0_dp
          dspiduy = yi / proji
       else
          proji = sqrt(xi2 + yi2)
          proji3 = proji*proji*proji
          dcpidx = 1.0_dp / proji - xi2 / proji3
          dcpidy = - xi * yi / proji3
          dcpidux = xi / proji - (xi2 * xi) / proji3
          dcpiduy = - (xi * yi2) / proji3
          dspidx = - xi * yi / proji3
          dspidy = 1.0_dp / proji - yi2 / proji3
          dspidux = - (yi * xi2) / proji3
          dspiduy = yi / proji - (yi2 * yi) / proji3
       endif

       cpi = xi / proji
       dcpidz = 0.0_dp
       dcpiduz = 0.0_dp

       spi = yi / proji
       dspidz = 0.0_dp
       dspiduz = 0.0_dp




    RealType sigma0 = mixer.sigma0;
    RealType dw     = mixer.dw;
    RealType eps0   = mixer.eps0;  
    RealType x2     = mixer.x2;    
    RealType xa2    = mixer.xa2;   
    RealType xai2   = mixer.xai2;  
    RealType xp2    = mixer.xp2;   
    RealType xpap2  = mixer.xpap2; 
    RealType xpapi2 = mixer.xpapi2;

    Vector3d ul1 = A1.getRow(2);
    Vector3d ul2 = A2.getRow(2);

    RealType a, b, g;


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

  void SHAPES::do_gb_pair(int *atid1, int *atid2, RealType *d, RealType *r, 
                      RealType *r2, RealType *sw, RealType *vdwMult,
                      RealType *vpair, RealType *pot, RealType *A1,
                      RealType *A2, RealType *f1, RealType *t1, RealType *t2) {
    
    if (!initialized_) initialize();
    
    AtomType* atype1 = SHAPESMap[*atid1];
    AtomType* atype2 = SHAPESMap[*atid2];
    
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
              Ai, Aj, frc, trq1, trq1);
      
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
#define fortranDoSHAPESPair FC_FUNC(do_gb_pair, DO_SHAPES_PAIR)
  
  RealType fortranGetGayBerneCut(int* atid) {
    return OpenMD::SHAPES::Instance()->getGayBerneCut(*atid);
  }

  void fortranDoSHAPESPair(int *atid1, int *atid2, RealType *d, RealType *r, 
                       RealType *r2, RealType *sw, RealType *vdwMult, 
                       RealType *vpair, RealType *pot, RealType *A1, 
                       RealType *A2, RealType *f1, RealType *t1, RealType *t2){
    
    return OpenMD::SHAPES::Instance()->do_gb_pair(atid1, atid2, d, r, r2, sw,
                                              vdwMult, vpair, pot, A1, A2, f1,
                                              t1, t2);
  }
}
