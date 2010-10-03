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
#include "nonbonded/EAM.hpp"
#include "utils/simError.h"
#include "types/NonBondedInteractionType.hpp"


namespace OpenMD {

  EAM::EAM() : name_("EAM"), initialized_(false), forceField_(NULL), 
               mixMeth_(eamJohnson), eamRcut_(0.0) {}
  
  EAMParam EAM::getEAMParam(AtomType* atomType) {
    
    // Do sanity checking on the AtomType we were passed before
    // building any data structures:
    if (!atomType->isEAM()) {
      sprintf( painCave.errMsg,
               "EAM::getEAMParam was passed an atomType (%s) that does not\n"
               "\tappear to be an embedded atom method (EAM) atom.\n",
               atomType->getName().c_str());
      painCave.severity = OPENMD_ERROR;
      painCave.isFatal = 1;
      simError();
    }
    
    GenericData* data = atomType->getPropertyByName("EAM");
    if (data == NULL) {
      sprintf( painCave.errMsg, "EAM::getEAMParam could not find EAM\n"
               "\tparameters for atomType %s.\n", 
               atomType->getName().c_str());
      painCave.severity = OPENMD_ERROR;
      painCave.isFatal = 1;
      simError(); 
    }
    
    EAMParamGenericData* eamData = dynamic_cast<EAMParamGenericData*>(data);
    if (eamData == NULL) {
      sprintf( painCave.errMsg,
               "EAM::getEAMParam could not convert GenericData to EAMParam for\n"
               "\tatom type %s\n", atomType->getName().c_str());
      painCave.severity = OPENMD_ERROR;
      painCave.isFatal = 1;
      simError();          
    }
    
    return eamData->getData();
  }

  CubicSpline* EAM::getZ(AtomType* atomType) {    
    EAMParam eamParam = getEAMParam(atomType);
    int nr = eamParam.nr;
    RealType dr = eamParam.dr;
    vector<RealType> rvals;
    
    for (int i = 0; i < nr; i++) rvals.push_back(RealType(i) * dr);
      
    CubicSpline* cs = new CubicSpline();
    cs->addPoints(rvals, eamParam.Z);
    return cs;
  }

  RealType EAM::getRcut(AtomType* atomType) {    
    EAMParam eamParam = getEAMParam(atomType);
    return eamParam.rcut;
  }

  CubicSpline* EAM::getRho(AtomType* atomType) {    
    EAMParam eamParam = getEAMParam(atomType);
    int nr = eamParam.nr;
    RealType dr = eamParam.dr;
    vector<RealType> rvals;
    
    for (int i = 0; i < nr; i++) rvals.push_back(RealType(i) * dr);
      
    CubicSpline* cs = new CubicSpline();
    cs->addPoints(rvals, eamParam.rho);
    return cs;
  }

  CubicSpline* EAM::getF(AtomType* atomType) {    
    EAMParam eamParam = getEAMParam(atomType);
    int nrho = eamParam.nrho;
    RealType drho = eamParam.drho;
    vector<RealType> rhovals;
    vector<RealType> scaledF;
    
    for (int i = 0; i < nrho; i++) {
      rhovals.push_back(RealType(i) * drho);
      scaledF.push_back( eamParam.F[i] * 23.06054 );
    }
      
    CubicSpline* cs = new CubicSpline();
    cs->addPoints(rhovals, scaledF);
    return cs;
  }
  
  CubicSpline* EAM::getPhi(AtomType* atomType1, AtomType* atomType2) {    
    EAMParam eamParam1 = getEAMParam(atomType1);
    EAMParam eamParam2 = getEAMParam(atomType2);
    CubicSpline* z1 = getZ(atomType1);
    CubicSpline* z2 = getZ(atomType2);

    // make the r grid:


    // we need phi out to the largest value we'll encounter in the radial space;
    
    RealType rmax = 0.0;
    rmax = max(rmax, eamParam1.rcut);
    rmax = max(rmax, eamParam1.nr * eamParam1.dr);

    rmax = max(rmax, eamParam2.rcut);
    rmax = max(rmax, eamParam2.nr * eamParam2.dr);

    // use the smallest dr (finest grid) to build our grid:

    RealType dr = min(eamParam1.dr, eamParam2.dr); 

    int nr = int(rmax/dr + 0.5);

    vector<RealType> rvals;
    for (int i = 0; i < nr; i++) rvals.push_back(RealType(i*dr));

    // construct the pair potential:

    vector<RealType> phivals;
    RealType phi;
    RealType r;
    RealType zi, zj;

    phivals.push_back(0.0);

    for (int i = 1; i < rvals.size(); i++ ) {
      r = rvals[i];

      // only use z(r) if we're inside this atom's cutoff radius,
      // otherwise, we'll use zero for the charge.  This effectively
      // means that our phi grid goes out beyond the cutoff of the
      // pair potential

      zi = r <= eamParam1.rcut ? z1->getValueAt(r) : 0.0;
      zj = r <= eamParam2.rcut ? z2->getValueAt(r) : 0.0;

      phi = 331.999296 * (zi * zj) / r;

      phivals.push_back(phi);
    }
      
    CubicSpline* cs = new CubicSpline();
    cs->addPoints(rvals, phivals);
    return cs;
  }

  void EAM::initialize() { 

    // set up the mixing method:
    ForceFieldOptions& fopts = forceField_->getForceFieldOptions();
    string EAMMixMeth = fopts.getEAMMixingMethod();
    toUpper(EAMMixMeth);
   
    if (EAMMixMeth == "JOHNSON") 
      mixMeth_ = eamJohnson;    
    else if (EAMMixMeth == "DAW")
      mixMeth_ = eamDaw;
    else
      mixMeth_ = eamUnknown;
      
    // find all of the EAM atom Types:
    ForceField::AtomTypeContainer* atomTypes = forceField_->getAtomTypes();
    ForceField::AtomTypeContainer::MapTypeIterator i;
    AtomType* at;

    for (at = atomTypes->beginType(i); at != NULL; 
         at = atomTypes->nextType(i)) {
      
      if (at->isEAM())
        addType(at);
    }
    
    // find all of the explicit EAM interactions (setfl):
    ForceField::NonBondedInteractionTypeContainer* nbiTypes = forceField_->getNonBondedInteractionTypes();
    ForceField::NonBondedInteractionTypeContainer::MapTypeIterator j;
    NonBondedInteractionType* nbt;

    for (nbt = nbiTypes->beginType(j); nbt != NULL; 
         nbt = nbiTypes->nextType(j)) {
      
      if (nbt->isEAM()) {
        
        pair<AtomType*, AtomType*> atypes = nbt->getAtomTypes();
        
        GenericData* data = nbt->getPropertyByName("EAM");
        if (data == NULL) {
          sprintf( painCave.errMsg, "EAM::rebuildMixingMap could not find\n"
                   "\tEAM parameters for %s - %s interaction.\n", 
                   atypes.first->getName().c_str(),
                   atypes.second->getName().c_str());
          painCave.severity = OPENMD_ERROR;
          painCave.isFatal = 1;
          simError(); 
        }
        
        EAMMixingData* eamData = dynamic_cast<EAMMixingData*>(data);
        if (eamData == NULL) {
          sprintf( painCave.errMsg,
                   "EAM::rebuildMixingMap could not convert GenericData to\n"
                   "\tEAMMixingData for %s - %s interaction.\n", 
                   atypes.first->getName().c_str(),
                   atypes.second->getName().c_str());
          painCave.severity = OPENMD_ERROR;
          painCave.isFatal = 1;
          simError();          
        }
        
        EAMMixingParam eamParam = eamData->getData();

        vector<RealType> phiAB = eamParam.phi;
        RealType dr = eamParam.dr;
        int nr = eamParam.nr;

        addExplicitInteraction(atypes.first, atypes.second, dr, nr, phiAB);
      }
    }  
    initialized_ = true;
  }
      


  void EAM::addType(AtomType* atomType){

    EAMAtomData eamAtomData;
    
    eamAtomData.rho = getRho(atomType);
    eamAtomData.F = getF(atomType);
    eamAtomData.Z = getZ(atomType);
    eamAtomData.rcut = getRcut(atomType);

    // add it to the map:
    AtomTypeProperties atp = atomType->getATP();    

    pair<map<int,AtomType*>::iterator,bool> ret;    
    ret = EAMlist.insert( pair<int, AtomType*>(atp.ident, atomType) );
    if (ret.second == false) {
      sprintf( painCave.errMsg,
               "EAM already had a previous entry with ident %d\n",
               atp.ident);
      painCave.severity = OPENMD_INFO;
      painCave.isFatal = 0;
      simError();         
    }

    EAMMap[atomType] = eamAtomData;
    
    // Now, iterate over all known types and add to the mixing map:
    
    map<AtomType*, EAMAtomData>::iterator it;
    for( it = EAMMap.begin(); it != EAMMap.end(); ++it) {
      
      AtomType* atype2 = (*it).first;

      EAMInteractionData mixer;
      mixer.phi = getPhi(atomType, atype2);
      mixer.explicitlySet = false;

      pair<AtomType*, AtomType*> key1, key2;
      key1 = make_pair(atomType, atype2);
      key2 = make_pair(atype2, atomType);
      
      MixingMap[key1] = mixer;
      if (key2 != key1) {
        MixingMap[key2] = mixer;
      }
    }      
    return;
  }
  
  void EAM::addExplicitInteraction(AtomType* atype1, AtomType* atype2, 
                                   RealType dr, int nr,
                                   vector<RealType> phiVals) {
    
    // in case these weren't already in the map
    addType(atype1);
    addType(atype2);

    EAMInteractionData mixer;
    CubicSpline* cs = new CubicSpline();
    vector<RealType> rVals;

    for (int i = 0; i < nr; i++) rVals.push_back(i * dr);

    cs->addPoints(rVals, phiVals);
    mixer.phi = cs;
    mixer.explicitlySet = true;

    pair<AtomType*, AtomType*> key1, key2;
    key1 = make_pair(atype1, atype2);
    key2 = make_pair(atype2, atype1);
    
    MixingMap[key1] = mixer;
    if (key2 != key1) {
      MixingMap[key2] = mixer;
    }    
    return;
  }

  void EAM::calcDensity(DensityData ddat) {
    
    if (!initialized_) initialize();
    
    EAMAtomData data1 = EAMMap[ddat.atype1];
    EAMAtomData data2 = EAMMap[ddat.atype2];

    if (ddat.rij < data1.rcut) 
      ddat.rho_i_at_j = data1.rho->getValueAt(ddat.rij);

    if (ddat.rij < data2.rcut) 
      ddat.rho_j_at_i = data2.rho->getValueAt(ddat.rij);

    return;
  }

  void EAM::calcFunctional(FunctionalData fdat) {

    if (!initialized_) initialize();

    EAMAtomData data1 = EAMMap[fdat.atype];
        
    pair<RealType, RealType> result = data1.F->getValueAndDerivativeAt(fdat.rho);

    fdat.frho = result.first;
    fdat.dfrhodrho = result.second;
    return;
  }

 
  void EAM::calcForce(InteractionData idat) {

    if (!initialized_) initialize();

    pair<RealType, RealType> res;
    
    if (idat.rij < eamRcut_) {

      EAMAtomData data1 = EAMMap[idat.atype1];
      EAMAtomData data2 = EAMMap[idat.atype2];

      // get type-specific cutoff radii

      RealType rci = data1.rcut;
      RealType rcj = data2.rcut;
      
      RealType rha, drha, rhb, drhb;
      RealType pha, dpha, phb, dphb;
      RealType phab, dvpdr;
      RealType drhoidr, drhojdr, dudr;
      
      if (idat.rij < rci) {
        res = data1.rho->getValueAndDerivativeAt(idat.rij);
        rha = res.first;
        drha = res.second;

        res = MixingMap[make_pair(idat.atype1, idat.atype1)].phi->getValueAndDerivativeAt(idat.rij);
        pha = res.first;
        dpha = res.second;
      }

      if (idat.rij < rcj) {
        res = data2.rho->getValueAndDerivativeAt(idat.rij);
        rhb = res.first;
        drhb = res.second;

        res = MixingMap[make_pair(idat.atype2, idat.atype2)].phi->getValueAndDerivativeAt(idat.rij);
        phb = res.first;
        dphb = res.second;
      }

      phab = 0.0;
      dvpdr = 0.0;

      switch(mixMeth_) {
      case eamJohnson:
       
        if (idat.rij < rci) {
          phab = phab + 0.5 * (rhb / rha) * pha;
          dvpdr = dvpdr + 0.5*((rhb/rha)*dpha + 
                               pha*((drhb/rha) - (rhb*drha/rha/rha)));
        }

        if (idat.rij < rcj) {
          phab = phab + 0.5 * (rha / rhb) * phb;
          dvpdr = dvpdr + 0.5 * ((rha/rhb)*dphb + 
                                 phb*((drha/rhb) - (rha*drhb/rhb/rhb)));
        }

        break;

      case eamDaw:
        res = MixingMap[make_pair(idat.atype1,idat.atype2)].phi->getValueAndDerivativeAt(idat.rij);
        phab = res.first;
        dvpdr = res.second;

        break;
      case eamUnknown:
      default:

        sprintf(painCave.errMsg,
                "EAM::calcForce hit a mixing method it doesn't know about!\n"
                );
        painCave.severity = OPENMD_ERROR;
        painCave.isFatal = 1;
        simError();        
          
      }
      
      drhoidr = drha;
      drhojdr = drhb;

      dudr = drhojdr*idat.dfrho1 + drhoidr*idat.dfrho2 + dvpdr; 

      idat.f1 = idat.d * dudr / idat.rij;
        
      // particle_pot is the difference between the full potential 
      // and the full potential without the presence of a particular
      // particle (atom1).
      //
      // This reduces the density at other particle locations, so
      // we need to recompute the density at atom2 assuming atom1
      // didn't contribute.  This then requires recomputing the
      // density functional for atom2 as well.
      //
      // Most of the particle_pot heavy lifting comes from the
      // pair interaction, and will be handled by vpair.
     
      idat.fshift1 = data1.F->getValueAt( idat.rho1 - rhb );
      idat.fshift2 = data1.F->getValueAt( idat.rho2 - rha );

      idat.pot += phab;

      idat.vpair += phab;
    }

    return;
    
  }

  RealType EAM::getSuggestedCutoffRadius(AtomType* at1, AtomType* at2) {
    if (!initialized_) initialize();   

    RealType cut = 0.0;

    map<AtomType*, EAMAtomData>::iterator it;

    it = EAMMap.find(at1);
    if (it != EAMMap.end()) {
      EAMAtomData data1 = (*it).second;
      cut = data1.rcut;
    }

    it = EAMMap.find(at2);
    if (it != EAMMap.end()) {
      EAMAtomData data2 = (*it).second;
      if (data2.rcut > cut)
        cut = data2.rcut;
    }

    return cut;
  }
}

