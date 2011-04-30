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

#include "nonbonded/MAW.hpp"
#include "utils/simError.h"

using namespace std;

namespace OpenMD {

  MAW::MAW() : name_("MAW"), initialized_(false), forceField_(NULL), 
                   shiftedPot_(false), shiftedFrc_(false) {}
  
  void MAW::initialize() {    

    ForceField::NonBondedInteractionTypeContainer* nbiTypes = forceField_->getNonBondedInteractionTypes();
    ForceField::NonBondedInteractionTypeContainer::MapTypeIterator j;
    NonBondedInteractionType* nbt;

    for (nbt = nbiTypes->beginType(j); nbt != NULL; 
         nbt = nbiTypes->nextType(j)) {
      
      if (nbt->isMAW()) {
        pair<AtomType*, AtomType*> atypes = nbt->getAtomTypes();
        
        GenericData* data = nbt->getPropertyByName("MAW");
        if (data == NULL) {
          sprintf( painCave.errMsg, "MAW::initialize could not find\n"
                   "\tMAW parameters for %s - %s interaction.\n", 
                   atypes.first->getName().c_str(),
                   atypes.second->getName().c_str());
          painCave.severity = OPENMD_ERROR;
          painCave.isFatal = 1;
          simError(); 
        }
        
        MAWData* mawData = dynamic_cast<MAWData*>(data);
        if (mawData == NULL) {
          sprintf( painCave.errMsg,
                   "MAW::initialize could not convert GenericData to\n"
                   "\tMAWData for %s - %s interaction.\n", 
                   atypes.first->getName().c_str(),
                   atypes.second->getName().c_str());
          painCave.severity = OPENMD_ERROR;
          painCave.isFatal = 1;
          simError();          
        }
        
        MAWParam mawParam = mawData->getData();

        RealType De = mawParam.De;
        RealType beta = mawParam.beta;
        RealType Re = mawParam.Re;
        RealType ca1 = mawParam.ca1;
        RealType cb1 = mawParam.cb1;
        
        addExplicitInteraction(atypes.first, atypes.second, 
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

    pair<AtomType*, AtomType*> key1, key2;
    key1 = make_pair(atype1, atype2);
    key2 = make_pair(atype2, atype1);
    
    MixingMap[key1] = mixer;
    if (key2 != key1) {
      MixingMap[key2] = mixer;
    }    
  }
  
  void MAW::calcForce(InteractionData &idat) {

    if (!initialized_) initialize();
    
    map<pair<AtomType*, AtomType*>, MAWInteractionData>::iterator it;
    it = MixingMap.find( *(idat.atypes) );
    if (it != MixingMap.end()) {
      MAWInteractionData mixer = (*it).second;
      
      RealType myPot = 0.0;
      RealType myPotC = 0.0;
      RealType myDeriv = 0.0;
      RealType myDerivC = 0.0;
      
      RealType D_e = mixer.De;
      RealType R_e = mixer.Re;
      RealType beta = mixer.beta;
      RealType ca1 = mixer.ca1;
      RealType cb1 = mixer.cb1;

      bool j_is_Metal = idat.atypes->second->isMetal();

      Vector3d r;
      RotMat3x3d Atrans;
      if (j_is_Metal) {
        // rotate the inter-particle separation into the two different 
        // body-fixed coordinate systems:
        r = *(idat.A1) * *(idat.d);
        Atrans = idat.A1->transpose();
      } else {
        // negative sign because this is the vector from j to i:       
        r = -*(idat.A2) * *(idat.d);
        Atrans = idat.A2->transpose();
      }
      
      // V(r) = D_e exp(-a(r-re)(exp(-a(r-re))-2)
      
      RealType expt     = -beta*( *(idat.rij) - R_e);
      RealType expfnc   = exp(expt);
      RealType expfnc2  = expfnc*expfnc;
      
      RealType exptC = 0.0;
      RealType expfncC = 0.0;
      RealType expfnc2C = 0.0;

      myPot  = D_e * (expfnc2  - 2.0 * expfnc);
      myDeriv   = 2.0 * D_e * beta * (expfnc - expfnc2);
      
      if (MAW::shiftedPot_ || MAW::shiftedFrc_) {
        exptC     = -beta*( *(idat.rcut)  - R_e);
        expfncC   = exp(exptC);
        expfnc2C  = expfncC*expfncC;
      }
      
      if (MAW::shiftedPot_) {
        myPotC = D_e * (expfnc2C - 2.0 * expfncC);
        myDerivC = 0.0;
      } else if (MAW::shiftedFrc_) {
        myPotC = D_e * (expfnc2C - 2.0 * expfncC);
        myDerivC  = 2.0 * D_e * beta * (expfnc2C - expfnc2C);
        myPotC += myDerivC * ( *(idat.rij)  -  *(idat.rcut) );
      } else {
        myPotC = 0.0;
        myDerivC = 0.0;
      }

      RealType x = r.x();
      RealType y = r.y();
      RealType z = r.z();
      RealType x2 = x * x;
      RealType y2 = y * y;
      RealType z2 = z * z;

      RealType r3 = *(idat.r2) *  *(idat.rij) ;
      RealType r4 = *(idat.r2) *  *(idat.r2);

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
      RealType Vang = ca1 * x2 / *(idat.r2) + 
        cb1 * z /  *(idat.rij)  + (0.8 - ca1 / 3.0);
      
      RealType pot_temp = *(idat.vdwMult) * Vmorse * Vang;
      *(idat.vpair) += pot_temp;
      idat.pot[VANDERWAALS_FAMILY] += *(idat.sw) * pot_temp;
           
      Vector3d dVmorsedr = (myDeriv - myDerivC) * Vector3d(x, y, z) /  *(idat.rij) ;
      
      Vector3d dVangdr = Vector3d(-2.0 * ca1 * x2 * x / r4 + 2.0 * ca1 * x / *(idat.r2) - cb1 * x * z / r3,
                                  -2.0 * ca1 * x2 * y / r4 - cb1 * z * y / r3,
                                  -2.0 * ca1 * x2 * z / r4 + cb1 /  *(idat.rij)           - cb1 * z2  / r3);
      
      // chain rule to put these back on x, y, z

      Vector3d dvdr = Vang * dVmorsedr + Vmorse * dVangdr;

      // Torques for Vang using method of Price:
      // S. L. Price, A. J. Stone, and M. Alderton, Mol. Phys. 52, 987 (1984).
      
      Vector3d dVangdu = Vector3d(cb1 * y /  *(idat.rij) ,
                                  2.0 * ca1 * x * z / *(idat.r2) - cb1 * x /  *(idat.rij),
                                  -2.0 * ca1 * y * x / *(idat.r2));

      // do the torques first since they are easy:
      // remember that these are still in the body fixed axes    

      Vector3d trq = *(idat.vdwMult) * Vmorse * dVangdu * *(idat.sw);

      // go back to lab frame using transpose of rotation matrix:
      
      if (j_is_Metal) {
        *(idat.t1) += Atrans * trq;
      } else {
        *(idat.t2) += Atrans * trq;
      }

      // Now, on to the forces (still in body frame of water)

      Vector3d ftmp = *(idat.vdwMult) * *(idat.sw) * dvdr;

      // rotate the terms back into the lab frame:
      Vector3d flab;
      if (j_is_Metal) {
        flab = Atrans * ftmp;
      } else {
        flab = - Atrans * ftmp;
      }
      
      *(idat.f1) += flab;
    }
    return;
    
  }
    
  RealType MAW::getSuggestedCutoffRadius(pair<AtomType*, AtomType*> atypes) {
    if (!initialized_) initialize();   
    map<pair<AtomType*, AtomType*>, MAWInteractionData>::iterator it;
    it = MixingMap.find(atypes);
    if (it == MixingMap.end()) 
      return 0.0;
    else  {
      MAWInteractionData mixer = (*it).second;

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

