/*
 * Copyright (c) 2004-2020 The University of Notre Dame. All Rights Reserved.
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
 * [1] Meineke, et al., J. Comp. Chem. 26, 252-271 (2005).
 * [2] Fennell & Gezelter, J. Chem. Phys. 124, 234104 (2006).
 * [3] Sun, Lin & Gezelter, J. Chem. Phys. 128, 234107 (2008).
 * [4] Vardeman, Stocker & Gezelter, J. Chem. Theory Comput. 7, 834 (2011).
 * [5] Kuang & Gezelter, Mol. Phys., 110, 691-701 (2012).
 * [6] Lamichhane, Gezelter & Newman, J. Chem. Phys. 141, 134109 (2014).
 * [7] Lamichhane, Newman & Gezelter, J. Chem. Phys. 141, 134110 (2014).
 * [8] Bhattarai, Newman & Gezelter, Phys. Rev. B 99, 094106 (2019).
 */

#include <cstdio>
#include <cstring>
#include <cmath>

#include "nonbonded/MAW.hpp"
#include "utils/simError.h"

using namespace std;

namespace OpenMD {

  MAW::MAW() : initialized_(false), forceField_(NULL), name_("MAW") {}
  
  void MAW::initialize() {    

    MAWtypes.clear();
    MAWtids.clear();
    MixingMap.clear();
    MAWtids.resize( forceField_->getNAtomType(), -1);

    ForceField::NonBondedInteractionTypeContainer* nbiTypes = forceField_->getNonBondedInteractionTypes();
    ForceField::NonBondedInteractionTypeContainer::MapTypeIterator j;
    ForceField::NonBondedInteractionTypeContainer::KeyType keys;
    NonBondedInteractionType* nbt;
    int mtid1, mtid2;

    for (nbt = nbiTypes->beginType(j); nbt != NULL; 
         nbt = nbiTypes->nextType(j)) {
      
      if (nbt->isMAW()) {
        keys = nbiTypes->getKeys(j);
        AtomType* at1 = forceField_->getAtomType(keys[0]);
        if (at1 == NULL) {
          sprintf( painCave.errMsg,
                   "MAW::initialize could not find AtomType %s\n"
                   "\tto for for %s - %s interaction.\n",
                   keys[0].c_str(), keys[0].c_str(), keys[1].c_str());
          painCave.severity = OPENMD_ERROR;
          painCave.isFatal = 1;
          simError();          
        }

        AtomType* at2 = forceField_->getAtomType(keys[1]);
        if (at2 == NULL) {
          sprintf( painCave.errMsg,
                   "MAW::initialize could not find AtomType %s\n"
                   "\tfor %s - %s nonbonded interaction.\n",
                   keys[1].c_str(), keys[0].c_str(), keys[1].c_str());
          painCave.severity = OPENMD_ERROR;
          painCave.isFatal = 1;
          simError();          
        }

        int atid1 = at1->getIdent();
        if (MAWtids[atid1] == -1) {
          mtid1 = MAWtypes.size();
          MAWtypes.insert(atid1);
          MAWtids[atid1] = mtid1;         
        }
        int atid2 = at2->getIdent();
        if (MAWtids[atid2] == -1) {
          mtid2 = MAWtypes.size();
          MAWtypes.insert(atid2);
          MAWtids[atid2] = mtid2;
        }

        MAWInteractionType* mit = dynamic_cast<MAWInteractionType*>(nbt);

        if (mit == NULL) {
          sprintf( painCave.errMsg,
                   "MAW::initialize could not convert NonBondedInteractionType\n"
                   "\tto MAWInteractionType for %s - %s interaction.\n", 
                   at1->getName().c_str(),
                   at2->getName().c_str());
          painCave.severity = OPENMD_ERROR;
          painCave.isFatal = 1;
          simError();          
        }
        
        RealType De = mit->getD();
        RealType beta = mit->getBeta();
        RealType Re = mit->getR();
        RealType ca1 = mit->getCA1();
        RealType cb1 = mit->getCB1();
        
        addExplicitInteraction(at1, at2, 
                               De, beta, Re, ca1, cb1);
      }
    }  
    initialized_ = true;
  }
  
  void MAW::addExplicitInteraction(AtomType* atype1, AtomType* atype2, 
                                   RealType De, RealType beta, RealType Re, 
                                   RealType ca1, RealType cb1) {

    MAWInteractionData mixer;
    mixer.De = De;
    mixer.beta = beta;
    mixer.Re = Re;
    mixer.ca1 = ca1;
    mixer.cb1 = cb1;
    mixer.j_is_Metal = atype2->isMetal();

    int mtid1 = MAWtids[atype1->getIdent()];
    int mtid2 = MAWtids[atype2->getIdent()];
    int nM = MAWtypes.size();

    MixingMap.resize(nM);
    MixingMap[mtid1].resize(nM);
    
    MixingMap[mtid1][mtid2] = mixer;
    if (mtid2 != mtid1) {
      MixingMap[mtid2].resize(nM);
      mixer.j_is_Metal = atype1->isMetal();
      MixingMap[mtid2][mtid1] = mixer;
    }    
  }
  
  void MAW::calcForce(InteractionData &idat) {

    if (!initialized_) initialize();

    MAWInteractionData &mixer = MixingMap[MAWtids[idat.atid1]][MAWtids[idat.atid2]];

    RealType myPot = 0.0;
    RealType myPotC = 0.0;
    RealType myDeriv = 0.0;
    RealType myDerivC = 0.0;
    
    RealType D_e = mixer.De;
    RealType R_e = mixer.Re;
    RealType beta = mixer.beta;
    RealType ca1 = mixer.ca1;
    RealType cb1 = mixer.cb1;   
    
    Vector3d r;
    RotMat3x3d Atrans;
    if (mixer.j_is_Metal) {
      // rotate the inter-particle separation into the two different 
      // body-fixed coordinate systems:
      r = idat.A1 * idat.d;
      Atrans = idat.A1.transpose();
    } else {
      // negative sign because this is the vector from j to i:       
      r = - idat.A2 * idat.d;
      Atrans = idat.A2.transpose();
    }
    
    // V(r) = D_e exp(-a(r-re)(exp(-a(r-re))-2)
    
    RealType expt     = -beta*( idat.rij - R_e);
    RealType expfnc   = exp(expt);
    RealType expfnc2  = expfnc*expfnc;
    
    RealType exptC = 0.0;
    RealType expfncC = 0.0;
    RealType expfnc2C = 0.0;
    
    myPot  = D_e * (expfnc2  - 2.0 * expfnc);
    myDeriv   = 2.0 * D_e * beta * (expfnc - expfnc2);
    
    if (idat.shiftedPot || idat.shiftedForce) {
      exptC     = -beta*( idat.rcut  - R_e);
      expfncC   = exp(exptC);
      expfnc2C  = expfncC*expfncC;
    }
    
    if (idat.shiftedPot) {
      myPotC = D_e * (expfnc2C - 2.0 * expfncC);
      myDerivC = 0.0;
    } else if (idat.shiftedForce) {
      myPotC = D_e * (expfnc2C - 2.0 * expfncC);
      myDerivC  = 2.0 * D_e * beta * (expfncC - expfnc2C);
      myPotC += myDerivC * ( idat.rij  -  idat.rcut );
    } else {
      myPotC = 0.0;
      myDerivC = 0.0;
    }
    
    RealType x = r.x();
    RealType y = r.y();
    RealType z = r.z();
    RealType x2 = x * x;
    RealType z2 = z * z;
    
    RealType r3 = idat.r2 *  idat.rij;
    RealType r4 = idat.r2 *  idat.r2;
    
    // angular modulation of morse part of potential to approximate 
    // the squares of the two HOMO lone pair orbitals in water:
    //********************** old form*************************
    // s = 1 / (4 pi)
    // ta1 = (s - pz)^2 = (1 - sqrt(3)*cos(theta))^2 / (4 pi)
    // b1 = px^2 = 3 * (sin(theta)*cos(phi))^2 / (4 pi)   
    //********************** old form*************************
    // we'll leave out the 4 pi for now
    
    // new functional form just using the p orbitals.
    // Vmorse(r)*[a*p_x + b p_z + (1-a-b)]
    // which is 
    // Vmorse(r)*[a sin^2(theta) cos^2(phi) + b cos(theta) + (1-a-b)]
    // Vmorse(r)*[a*x2/r2 + b*z/r + (1-a-b)]
    
    RealType Vmorse = (myPot - myPotC);
    RealType Vang = ca1 * x2 / idat.r2 + 
      cb1 * z /  idat.rij  + (0.8 - ca1 / 3.0);
    
    RealType pot_temp = idat.vdwMult * Vmorse * Vang;
    idat.vpair += pot_temp;
    idat.pot[VANDERWAALS_FAMILY] += idat.sw * pot_temp;
    if (idat.isSelected)
      idat.selePot[VANDERWAALS_FAMILY] += idat.sw * pot_temp;
    Vector3d dVmorsedr = (myDeriv - myDerivC) * Vector3d(x, y, z) /  idat.rij ;
    
    Vector3d dVangdr = Vector3d(-2.0 * ca1 * x2 * x / r4 + 2.0 * ca1 * x / idat.r2 - cb1 * x * z / r3,
                                -2.0 * ca1 * x2 * y / r4 - cb1 * z * y / r3,
                                -2.0 * ca1 * x2 * z / r4 + cb1 /  idat.rij           - cb1 * z2  / r3);
    
    // chain rule to put these back on x, y, z
    
    Vector3d dvdr = Vang * dVmorsedr + Vmorse * dVangdr;
    
    // Torques for Vang using method of Price:
    // S. L. Price, A. J. Stone, and M. Alderton, Mol. Phys. 52, 987 (1984).
    
    Vector3d dVangdu = Vector3d(cb1 * y /  idat.rij,
                                2.0 * ca1 * x * z / idat.r2 - cb1 * x /  idat.rij,
                                -2.0 * ca1 * y * x / idat.r2);
    
    // do the torques first since they are easy:
    // remember that these are still in the body fixed axes    
    
    Vector3d trq = idat.vdwMult * Vmorse * dVangdu * idat.sw;
    
    // go back to lab frame using transpose of rotation matrix:
    
    if (mixer.j_is_Metal) {
      idat.t1 += Atrans * trq;
    } else {
      idat.t2 += Atrans * trq;
    }
    
    // Now, on to the forces (still in body frame of water)
    
    Vector3d ftmp = idat.vdwMult * idat.sw * dvdr;
    
    // rotate the terms back into the lab frame:
    Vector3d flab;
    if (mixer.j_is_Metal) {
      flab = Atrans * ftmp;
    } else {
      flab = - Atrans * ftmp;
    }
    
    idat.f1 += flab;
    
    return;    
  }

  RealType MAW::getSuggestedCutoffRadius(pair<AtomType*, AtomType*> atypes) {
    if (!initialized_) initialize();   
    int atid1 = atypes.first->getIdent();
    int atid2 = atypes.second->getIdent();
    int mtid1 = MAWtids[atid1];
    int mtid2 = MAWtids[atid2];
    
    if ( mtid1 == -1 || mtid2 == -1) return 0.0;
    else {      
      MAWInteractionData mixer = MixingMap[mtid1][mtid2];
      RealType R_e = mixer.Re;
      RealType beta = mixer.beta;
      // This value of the r corresponds to an energy about 1.48% of 
      // the energy at the bottom of the Morse well.  For comparison, the
      // Lennard-Jones function is about 1.63% of it's minimum value at
      // a distance of 2.5 sigma.
      return (4.9 + beta * R_e) / beta;
    }
  }
}

