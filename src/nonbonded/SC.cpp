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
#include "nonbonded/SC.hpp"
#include "utils/simError.h"
#include "types/NonBondedInteractionType.hpp"

namespace OpenMD {


  SC::SC() : name_("SC"), initialized_(false), forceField_(NULL), 
             scRcut_(0.0), np_(3000) {}
  
  SCParam SC::getSCParam(AtomType* atomType) {
    
    // Do sanity checking on the AtomType we were passed before
    // building any data structures:
    if (!atomType->isSC()) {
      sprintf( painCave.errMsg,
               "SC::getSCParam was passed an atomType (%s) that does not\n"
               "\tappear to be a Sutton-Chen (SC) atom.\n",
               atomType->getName().c_str());
      painCave.severity = OPENMD_ERROR;
      painCave.isFatal = 1;
      simError();
    }
    
    GenericData* data = atomType->getPropertyByName("SC");
    if (data == NULL) {
      sprintf( painCave.errMsg, "SC::getSCParam could not find SC\n"
               "\tparameters for atomType %s.\n", 
               atomType->getName().c_str());
      painCave.severity = OPENMD_ERROR;
      painCave.isFatal = 1;
      simError(); 
    }
    
    SCParamGenericData* scData = dynamic_cast<SCParamGenericData*>(data);
    if (scData == NULL) {
      sprintf( painCave.errMsg,
               "SC::getSCParam could not convert GenericData to SCParamGenericData\n"
               "\tfor atom type %s\n", atomType->getName().c_str());
      painCave.severity = OPENMD_ERROR;
      painCave.isFatal = 1;
      simError();          
    }
    
    return scData->getData();
  }

  RealType SC::getC(AtomType* atomType) {    
    SCParam scParam = getSCParam(atomType);
    return scParam.c;
  }

  RealType SC::getM(AtomType* atomType) {    
    SCParam scParam = getSCParam(atomType);
    return scParam.m;
  }

  RealType SC::getM(AtomType* atomType1, AtomType* atomType2) {    
    RealType m1 = getM(atomType1);
    RealType m2 = getM(atomType2);
    return 0.5 * (m1 + m2);
  }

  RealType SC::getN(AtomType* atomType) {    
    SCParam scParam = getSCParam(atomType);
    return scParam.n;
  }

  RealType SC::getN(AtomType* atomType1, AtomType* atomType2) {    
    RealType n1 = getN(atomType1);
    RealType n2 = getN(atomType2);
    return 0.5 * (n1 + n2);
  }

  RealType SC::getAlpha(AtomType* atomType) {    
    SCParam scParam = getSCParam(atomType);
    return scParam.alpha;
  }

  RealType SC::getAlpha(AtomType* atomType1, AtomType* atomType2) {    
    RealType alpha1 = getAlpha(atomType1);
    RealType alpha2 = getAlpha(atomType2);

    ForceFieldOptions& fopts = forceField_->getForceFieldOptions();
    std::string DistanceMix = fopts.getDistanceMixingRule();
    toUpper(DistanceMix);

    if (DistanceMix == "GEOMETRIC") 
      return sqrt(alpha1 * alpha2);
    else 
      return 0.5 * (alpha1 + alpha2);
  }

  RealType SC::getEpsilon(AtomType* atomType) {    
    SCParam scParam = getSCParam(atomType);
    return scParam.epsilon;
  }

  RealType SC::getEpsilon(AtomType* atomType1, AtomType* atomType2) {    
    RealType epsilon1 = getEpsilon(atomType1);
    RealType epsilon2 = getEpsilon(atomType2);
    return sqrt(epsilon1 * epsilon2);
  }

  void SC::initialize() { 
    // find all of the SC atom Types:
    ForceField::AtomTypeContainer* atomTypes = forceField_->getAtomTypes();
    ForceField::AtomTypeContainer::MapTypeIterator i;
    AtomType* at;

    for (at = atomTypes->beginType(i); at != NULL; 
         at = atomTypes->nextType(i)) {
      if (at->isSC())
        addType(at);
    }    
    initialized_ = true;
  }
  
  

  void SC::addType(AtomType* atomType){

    SCAtomData scAtomData;
    
    scAtomData.c = getC(atomType);
    scAtomData.m = getM(atomType);
    scAtomData.n = getN(atomType);
    scAtomData.alpha = getAlpha(atomType);
    scAtomData.epsilon = getEpsilon(atomType);
    scAtomData.rCut = 2.0 * scAtomData.alpha;

    // add it to the map:
    AtomTypeProperties atp = atomType->getATP();    

    pair<map<int,AtomType*>::iterator,bool> ret;    
    ret = SClist.insert( pair<int, AtomType*>(atp.ident, atomType) );
    if (ret.second == false) {
      sprintf( painCave.errMsg,
               "SC already had a previous entry with ident %d\n",
               atp.ident);
      painCave.severity = OPENMD_INFO;
      painCave.isFatal = 0;
      simError();         
    }

    SCMap[atomType] = scAtomData;
    
    // Now, iterate over all known types and add to the mixing map:
    
    map<AtomType*, SCAtomData>::iterator it;
    for( it = SCMap.begin(); it != SCMap.end(); ++it) {
      
      AtomType* atype2 = (*it).first;
      
      SCInteractionData mixer;

      mixer.alpha = getAlpha(atomType, atype2);
      mixer.rCut = 2.0 * mixer.alpha;
      mixer.epsilon = getEpsilon(atomType, atype2);
      mixer.m = getM(atomType, atype2);
      mixer.n = getN(atomType, atype2);

      RealType dr = mixer.rCut / (np_ - 1);
      vector<RealType> rvals;
      vector<RealType> vvals;
      vector<RealType> phivals;
    
      rvals.push_back(0.0);
      vvals.push_back(0.0);
      phivals.push_back(0.0);

      for (int k = 1; k < np_; k++) {
        RealType r = dr * k;
        rvals.push_back(r);
        vvals.push_back( mixer.epsilon * pow(mixer.alpha/r, mixer.n) );
        phivals.push_back( pow(mixer.alpha/r, mixer.m) );
      }

      mixer.vCut = mixer.epsilon * pow(mixer.alpha/mixer.rCut, mixer.n);
    
      CubicSpline* V = new CubicSpline();
      V->addPoints(rvals, vvals);
      
      CubicSpline* phi = new CubicSpline();
      phi->addPoints(rvals, phivals);
      
      mixer.V = V;
      mixer.phi = phi;

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
  
  void SC::addExplicitInteraction(AtomType* atype1, AtomType* atype2, 
                                  RealType epsilon, RealType m, RealType n,
                                  RealType alpha) {
    
    // in case these weren't already in the map
    addType(atype1);
    addType(atype2);

    SCInteractionData mixer;

    mixer.epsilon = epsilon;
    mixer.m = m;
    mixer.n = n;
    mixer.alpha = alpha;
    mixer.rCut = 2.0 * mixer.alpha;
    
    RealType dr = mixer.rCut / (np_ - 1);
    vector<RealType> rvals;
    vector<RealType> vvals;
    vector<RealType> phivals;
    
    rvals.push_back(0.0);
    vvals.push_back(0.0);
    phivals.push_back(0.0);
    
    for (int k = 1; k < np_; k++) {
      RealType r = dr * k;
      rvals.push_back(r);
      vvals.push_back( mixer.epsilon * pow(mixer.alpha/r, mixer.n) );
      phivals.push_back( pow(mixer.alpha/r, mixer.m) );
    }
    
    mixer.vCut = mixer.epsilon * pow(mixer.alpha/mixer.rCut, mixer.n);
    
    CubicSpline* V = new CubicSpline();
    V->addPoints(rvals, vvals);
    
    CubicSpline* phi = new CubicSpline();
    phi->addPoints(rvals, phivals);
    
    mixer.V = V;
    mixer.phi = phi;
    
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

  void SC::calcDensity(DensityData ddat) {
    
    if (!initialized_) initialize();
    
    SCInteractionData mixer = MixingMap[make_pair(ddat.atype1, ddat.atype2)];

    RealType rcij = mixer.rCut;

    if (ddat.rij < rcij) {
      ddat.rho_i_at_j = mixer.phi->getValueAt(ddat.rij);
      ddat.rho_j_at_i = ddat.rho_i_at_j;
    } else {
      ddat.rho_i_at_j = 0.0;
      ddat.rho_j_at_i = 0.0;
    }

    return;
  }

  void SC::calcFunctional(FunctionalData fdat) {

    if (!initialized_) initialize();

    SCAtomData data1 = SCMap[fdat.atype];
    
    fdat.frho = - data1.c * data1.epsilon * sqrt(fdat.rho);
    fdat.dfrhodrho = 0.5 * fdat.frho / fdat.rho;
    
    return;
  }
  
 
  void SC::calcForce(InteractionData idat) {
    
    if (!initialized_) initialize();
    
    SCAtomData data1 = SCMap[idat.atype1];
    SCAtomData data2 = SCMap[idat.atype2];

    SCInteractionData mixer = MixingMap[make_pair(idat.atype1, idat.atype2)];

    RealType rcij = mixer.rCut;

    if (idat.rij < rcij) {
      RealType vcij = mixer.vCut;
      
      pair<RealType, RealType> res;
      
      res = mixer.phi->getValueAndDerivativeAt(idat.rij);
      RealType rhtmp = res.first;
      RealType drhodr = res.second;
      
      res = mixer.V->getValueAndDerivativeAt(idat.rij);
      RealType vptmp = res.first;
      RealType dvpdr = res.second;
      
      RealType pot_temp = vptmp - vcij;
      idat.vpair += pot_temp;
      
      RealType dudr = drhodr * (idat.dfrho1 + idat.dfrho2) + dvpdr;
      
      idat.f1 += idat.d * dudr / idat.rij;
        
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
      
      idat.fshift1 = - data1.c * data1.epsilon * sqrt(idat.rho1 - rhtmp);
      idat.fshift2 = - data2.c * data2.epsilon * sqrt(idat.rho2 - rhtmp);
      
      idat.pot += pot_temp;
    }
      
    return;    
  }

  RealType SC::getSuggestedCutoffRadius(AtomType* at1, AtomType* at2) {
    if (!initialized_) initialize();   
    pair<AtomType*, AtomType*> key = make_pair(at1, at2); 
    map<pair<AtomType*, AtomType*>, SCInteractionData>::iterator it;
    it = MixingMap.find(key);
    if (it == MixingMap.end()) 
      return 0.0;
    else  {
      SCInteractionData mixer = (*it).second;
      return mixer.rCut;
    }
  }
}
