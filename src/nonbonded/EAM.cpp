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
#include "nonbonded/EAM.hpp"
#include "utils/simError.h"
#include "types/NonBondedInteractionType.hpp"


namespace OpenMD {

  EAM::EAM() : name_("EAM"), initialized_(false), forceField_(NULL), 
               mixMeth_(eamJohnson), eamRcut_(0.0), haveCutoffRadius_(false) {}
  
  CubicSpline* EAM::getPhi(AtomType* atomType1, AtomType* atomType2) {   
    EAMAdapter ea1 = EAMAdapter(atomType1);
    EAMAdapter ea2 = EAMAdapter(atomType2);
    CubicSpline* z1 = ea1.getZ();
    CubicSpline* z2 = ea2.getZ();

    // Thise prefactors convert the charge-charge interactions into
    // kcal / mol all were computed assuming distances are measured in
    // angstroms Charge-Charge, assuming charges are measured in
    // electrons.  Matches value in Electrostatics.cpp
    pre11_ = 332.0637778;

    // make the r grid:

    // we need phi out to the largest value we'll encounter in the radial space;
    
    RealType rmax = 0.0;
    rmax = max(rmax, ea1.getRcut());
    rmax = max(rmax, ea1.getNr() * ea1.getDr());

    rmax = max(rmax, ea2.getRcut());
    rmax = max(rmax, ea2.getNr() * ea2.getDr());

    // use the smallest dr (finest grid) to build our grid:

    RealType dr = min(ea1.getDr(), ea2.getDr()); 

    int nr = int(rmax/dr + 0.5);

    vector<RealType> rvals;
    for (int i = 0; i < nr; i++) rvals.push_back(RealType(i*dr));

    // construct the pair potential:

    vector<RealType> phivals;
    RealType phi;
    RealType r;
    RealType zi, zj;

    phivals.push_back(0.0);

    for (unsigned int i = 1; i < rvals.size(); i++ ) {
      r = rvals[i];

      // only use z(r) if we're inside this atom's cutoff radius,
      // otherwise, we'll use zero for the charge.  This effectively
      // means that our phi grid goes out beyond the cutoff of the
      // pair potential

      zi = r <= ea1.getRcut() ? z1->getValueAt(r) : 0.0;
      zj = r <= ea2.getRcut() ? z2->getValueAt(r) : 0.0;

      phi = pre11_ * (zi * zj) / r;

      phivals.push_back(phi);
    }
      
    CubicSpline* cs = new CubicSpline();
    cs->addPoints(rvals, phivals);
    return cs;
  }

  void EAM::setCutoffRadius( RealType rCut ) {
    eamRcut_ = rCut;
    haveCutoffRadius_ = true;
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
    EAMtypes.clear();
    EAMtids.clear();
    EAMdata.clear();
    MixingMap.clear();
    nEAM_ = 0;
    
    EAMtids.resize( forceField_->getNAtomType(), -1);

    set<AtomType*>::iterator at;
    for (at = simTypes_.begin(); at != simTypes_.end(); ++at) {
      if ((*at)->isEAM()) nEAM_++;
    }
    EAMdata.resize(nEAM_);
    MixingMap.resize(nEAM_);

    for (at = simTypes_.begin(); at != simTypes_.end(); ++at) {
      if ((*at)->isEAM()) addType(*at);
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

    EAMAdapter ea = EAMAdapter(atomType);
    EAMAtomData eamAtomData;

    eamAtomData.rho = ea.getRho();
    eamAtomData.F = ea.getF();
    eamAtomData.Z = ea.getZ();
    eamAtomData.rcut = ea.getRcut();
    eamAtomData.isFluctuating = atomType->isFluctuatingCharge();
      
    // add it to the map:
    int atid = atomType->getIdent();
    int eamtid = EAMtypes.size();

    pair<set<int>::iterator,bool> ret;    
    ret = EAMtypes.insert( atid );
    if (ret.second == false) {
      sprintf( painCave.errMsg,
               "EAM already had a previous entry with ident %d\n",
               atid);
      painCave.severity = OPENMD_INFO;
      painCave.isFatal = 0;
      simError();         
    }

    if (eamAtomData.isFluctuating) {
      // compute charge to rho scaling:
      RealType z0 = eamAtomData.Z->getValueAt(0.0);
      RealType dr = ea.getDr();
      RealType rmax = max(eamAtomData.rcut, ea.getNr() * dr);
      int nr = int(rmax/dr + 0.5);
      RealType r;
      RealType sum(0.0);

      for (int i = 0; i < nr; i++) {
        r = RealType(i*dr);
        sum += r * r * eamAtomData.rho->getValueAt(r) * dr;       
      } 
      sum *= 4.0 * M_PI;
      eamAtomData.qToRhoScaling = sum / z0;
    }


    EAMtids[atid] = eamtid;
    EAMdata[eamtid] = eamAtomData;
    MixingMap[eamtid].resize(nEAM_);
    
    // Now, iterate over all known types and add to the mixing map:
    
    std::set<int>::iterator it;
    for( it = EAMtypes.begin(); it != EAMtypes.end(); ++it) {
      
      int eamtid2 = EAMtids[ (*it) ];
      AtomType* atype2 = forceField_->getAtomType( (*it) );

      EAMInteractionData mixer;
      mixer.phi = getPhi(atomType, atype2);
      mixer.explicitlySet = false;

      MixingMap[eamtid2].resize( nEAM_ );
      
      MixingMap[eamtid][eamtid2] = mixer;
      if (eamtid2 != eamtid) {
        MixingMap[eamtid2][eamtid] = mixer;
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

    int eamtid1 = EAMtids[ atype1->getIdent() ];
    int eamtid2 = EAMtids[ atype2->getIdent() ];
    
    MixingMap[eamtid1][eamtid2] = mixer;
    if (eamtid2 != eamtid1) {
      MixingMap[eamtid2][eamtid1] = mixer;
    }    
    return;
  }

  void EAM::calcDensity(InteractionData &idat) {
    
    if (!initialized_) initialize();
    
    EAMAtomData &data1 = EAMdata[EAMtids[idat.atid1]];
    EAMAtomData &data2 = EAMdata[EAMtids[idat.atid2]];

    if (haveCutoffRadius_) 
      if ( *(idat.rij) > eamRcut_) return;
    
    if ( *(idat.rij) < data1.rcut) {
      if (data1.isFluctuating) {
        *(idat.rho2) += (1.0 -  *(idat.flucQ1) * data1.qToRhoScaling ) * 
          data1.rho->getValueAt( *(idat.rij) );
      } else {
        *(idat.rho2) += data1.rho->getValueAt( *(idat.rij));
      }
    }
      
    if ( *(idat.rij) < data2.rcut) {
      if (data2.isFluctuating) {
        *(idat.rho1) += (1.0 -  *(idat.flucQ2) * data2.qToRhoScaling ) *
          data2.rho->getValueAt( *(idat.rij) );
      } else {
        *(idat.rho1) += data2.rho->getValueAt( *(idat.rij));
      }
    }
    
    return;  
  }
  
  void EAM::calcFunctional(SelfData &sdat) {
    
    if (!initialized_) initialize();

    EAMAtomData &data1 = EAMdata[ EAMtids[sdat.atid] ];
            
    data1.F->getValueAndDerivativeAt( *(sdat.rho), *(sdat.frho), *(sdat.dfrhodrho) );

    (*(sdat.pot))[METALLIC_FAMILY] += *(sdat.frho);
    if (sdat.doParticlePot) {
      *(sdat.particlePot) += *(sdat.frho);
    }

    return;
  }

 
  void EAM::calcForce(InteractionData &idat) {

    if (!initialized_) initialize();

    if (haveCutoffRadius_) 
      if ( *(idat.rij) > eamRcut_) return;
   

    int eamtid1 = EAMtids[idat.atid1];
    int eamtid2 = EAMtids[idat.atid2];
    
    EAMAtomData &data1 = EAMdata[eamtid1];
    EAMAtomData &data2 = EAMdata[eamtid2];
    
    // get type-specific cutoff radii
    
    RealType rci = data1.rcut;
    RealType rcj = data2.rcut;
    
    RealType rha(0.0), drha(0.0), rhb(0.0), drhb(0.0);
    RealType pha(0.0), dpha(0.0), phb(0.0), dphb(0.0);
    RealType phab(0.0), dvpdr(0.0);
    RealType drhoidr, drhojdr, dudr;
    
    if ( *(idat.rij) < rci) {
      data1.rho->getValueAndDerivativeAt( *(idat.rij), rha, drha);
      CubicSpline* phi = MixingMap[eamtid1][eamtid1].phi;
      phi->getValueAndDerivativeAt( *(idat.rij), pha, dpha);
      if (data1.isFluctuating) {
        *(idat.dVdFQ1) -= *(idat.dfrho2) * rha * data1.qToRhoScaling;
      }
    }
    
    if ( *(idat.rij) < rcj) {
      data2.rho->getValueAndDerivativeAt( *(idat.rij), rhb, drhb );
      CubicSpline* phi = MixingMap[eamtid2][eamtid2].phi;
      phi->getValueAndDerivativeAt( *(idat.rij), phb, dphb);
      if (data2.isFluctuating) {
        *(idat.dVdFQ2) -= *(idat.dfrho1) * rhb * data2.qToRhoScaling;
      }
    }

    switch(mixMeth_) {
    case eamJohnson:
      
      if ( *(idat.rij) < rci) {
        phab = phab + 0.5 * (rhb / rha) * pha;
        dvpdr = dvpdr + 0.5*((rhb/rha)*dpha + 
                             pha*((drhb/rha) - (rhb*drha/rha/rha)));
      }
      
      
      
      if ( *(idat.rij) < rcj) {
        phab = phab + 0.5 * (rha / rhb) * phb;
        dvpdr = dvpdr + 0.5 * ((rha/rhb)*dphb + 
                               phb*((drha/rhb) - (rha*drhb/rhb/rhb)));
      }
      
      break;
      
    case eamDaw:
      MixingMap[eamtid1][eamtid2].phi->getValueAndDerivativeAt( *(idat.rij), phab, dvpdr);
      
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
    
    dudr = drhojdr* *(idat.dfrho1) + drhoidr* *(idat.dfrho2) + dvpdr; 
    
    *(idat.f1) += *(idat.d) * dudr / *(idat.rij);

        
    if (idat.doParticlePot) {
      // particlePot is the difference between the full potential and
      // the full potential without the presence of a particular
      // particle (atom1).
      //
      // This reduces the density at other particle locations, so we
      // need to recompute the density at atom2 assuming atom1 didn't
      // contribute.  This then requires recomputing the density
      // functional for atom2 as well.
      
      *(idat.particlePot1) += data2.F->getValueAt( *(idat.rho2) - rha ) 
        - *(idat.frho2);
      
      *(idat.particlePot2) += data1.F->getValueAt( *(idat.rho1) - rhb) 
        - *(idat.frho1);
    }
    
    (*(idat.pot))[METALLIC_FAMILY] += phab;
    
    *(idat.vpair) += phab;
  
    return;
    
  }

  RealType EAM::getSuggestedCutoffRadius(pair<AtomType*, AtomType*> atypes) {
    if (!initialized_) initialize();   

    RealType cut = 0.0;

    int atid1 = atypes.first->getIdent();
    int atid2 = atypes.second->getIdent();
    int eamtid1 = EAMtids[atid1];
    int eamtid2 = EAMtids[atid2];
    
    if (eamtid1 != -1) {
      EAMAtomData data1 = EAMdata[eamtid1];
      cut = data1.rcut;
    }

    if (eamtid2 != -1) {
      EAMAtomData data2 = EAMdata[eamtid2];
      if (data2.rcut > cut)
        cut = data2.rcut;
    }
    
    return cut;
  }
}

