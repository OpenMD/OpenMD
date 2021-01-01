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

#include "nonbonded/SC.hpp"
#include "utils/simError.h"
#include "types/NonBondedInteractionType.hpp"

namespace OpenMD {


  SC::SC() : name_("SC"), initialized_(false), forceField_(NULL),  np_(3000) {}
  
  SC::~SC() {
    initialized_ = false;

    MixingMap.clear();
    SCtypes.clear();
    SCdata.clear();
    SCtids.clear();
  }
        
  RealType SC::getM(AtomType* atomType1, AtomType* atomType2) {    
    SuttonChenAdapter sca1 = SuttonChenAdapter(atomType1);
    SuttonChenAdapter sca2 = SuttonChenAdapter(atomType2);
    RealType m1 = sca1.getM();
    RealType m2 = sca2.getM();
    return 0.5 * (m1 + m2);
  }

  RealType SC::getN(AtomType* atomType1, AtomType* atomType2) {    
    SuttonChenAdapter sca1 = SuttonChenAdapter(atomType1);
    SuttonChenAdapter sca2 = SuttonChenAdapter(atomType2);
    RealType n1 = sca1.getN();
    RealType n2 = sca2.getN();
    return 0.5 * (n1 + n2);
  }

  RealType SC::getAlpha(AtomType* atomType1, AtomType* atomType2) {    
    SuttonChenAdapter sca1 = SuttonChenAdapter(atomType1);
    SuttonChenAdapter sca2 = SuttonChenAdapter(atomType2);
    RealType alpha1 = sca1.getAlpha();
    RealType alpha2 = sca2.getAlpha();

    ForceFieldOptions& fopts = forceField_->getForceFieldOptions();
    std::string DistanceMix = fopts.getDistanceMixingRule();
    toUpper(DistanceMix);

    if (DistanceMix == "GEOMETRIC") 
      return sqrt(alpha1 * alpha2);
    else 
      return 0.5 * (alpha1 + alpha2);
  }

  RealType SC::getEpsilon(AtomType* atomType1, AtomType* atomType2) {   
    SuttonChenAdapter sca1 = SuttonChenAdapter(atomType1);
    SuttonChenAdapter sca2 = SuttonChenAdapter(atomType2);
    RealType epsilon1 = sca1.getEpsilon();
    RealType epsilon2 = sca2.getEpsilon();
    return sqrt(epsilon1 * epsilon2);
  }

  void SC::initialize() {      
    // find all of the SC atom Types:
    SCtypes.clear();
    SCtids.clear();
    SCdata.clear();
    MixingMap.clear();
    nSC_ = 0;

    SCtids.resize( forceField_->getNAtomType(), -1);

    set<AtomType*>::iterator at;
    for (at = simTypes_.begin(); at != simTypes_.end(); ++at) {
      if ((*at)->isSC()) nSC_++;
    }
    SCdata.resize(nSC_);
    MixingMap.resize(nSC_);
    for (at = simTypes_.begin(); at != simTypes_.end(); ++at) {
      if ((*at)->isSC()) addType((*at));
    }
    initialized_ = true;
  }
  
  

  void SC::addType(AtomType* atomType){

    SuttonChenAdapter sca = SuttonChenAdapter(atomType);
    SCAtomData scAtomData;
    
    scAtomData.c = sca.getC();
    scAtomData.m = sca.getM();
    scAtomData.n = sca.getN();
    scAtomData.alpha = sca.getAlpha();
    scAtomData.epsilon = sca.getEpsilon();
    scAtomData.rCut = 2.0 * scAtomData.alpha;
 
    // add it to the map:
    int atid = atomType->getIdent();
    int sctid = SCtypes.size();

    pair<set<int>::iterator,bool> ret;    
    ret = SCtypes.insert( atid );
    if (ret.second == false) {
      sprintf( painCave.errMsg,
               "SC already had a previous entry with ident %d\n",
               atid );
      painCave.severity = OPENMD_INFO;
      painCave.isFatal = 0;
      simError();         
    }
    
    SCtids[atid] = sctid;
    SCdata[sctid] = scAtomData;
    MixingMap[sctid].resize(nSC_);
    
    // Now, iterate over all known types and add to the mixing map:
    
    std::set<int>::iterator it;
    for( it = SCtypes.begin(); it != SCtypes.end(); ++it) {
      
      int sctid2 = SCtids[ (*it) ];
      AtomType* atype2 = forceField_->getAtomType( (*it) );
      
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

      MixingMap[sctid2].resize( nSC_ );
      
      MixingMap[sctid][sctid2] = mixer;
      if (sctid2 != sctid) {
        MixingMap[sctid2][sctid] = mixer;
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

    int sctid1 = SCtids[ atype1->getIdent() ];
    int sctid2 = SCtids[ atype2->getIdent() ];

    MixingMap[sctid1][sctid2] = mixer;
    if (sctid2 != sctid1) {
      MixingMap[sctid2][sctid1] = mixer;
    }    
    return;
  }

  void SC::calcDensity(InteractionData &idat) {
    
    if (!initialized_) initialize();
    int sctid1 = SCtids[idat.atid1];
    int sctid2 = SCtids[idat.atid2];
    
    SCInteractionData &mixer = MixingMap[sctid1][sctid2];

    RealType rcij = mixer.rCut;

    if ( idat.rij  < rcij) {
      RealType rho = mixer.phi->getValueAt( idat.rij );
      idat.rho1 += rho;
      idat.rho2 += rho;
    } 
    
    return;
  }

  void SC::calcFunctional(SelfData &sdat) {

    if (!initialized_) initialize();

    SCAtomData &data1 = SCdata[SCtids[sdat.atid]];
   
    RealType u = - data1.c * data1.epsilon * sqrt( sdat.rho );
    sdat.frho = u;
    sdat.dfrhodrho = 0.5 * sdat.frho / sdat.rho;

    sdat.selfPot[METALLIC_EMBEDDING_FAMILY] += u;
    
    if (sdat.isSelected) 
      sdat.selePot[METALLIC_EMBEDDING_FAMILY] += u;
   
    if (sdat.doParticlePot) {
      sdat.particlePot += u;
    }

    return;
  }
  
 
  void SC::calcForce(InteractionData &idat) {
    
    if (!initialized_) initialize();
    
    int &sctid1 = SCtids[idat.atid1];
    int &sctid2 = SCtids[idat.atid2];

    SCAtomData &data1 = SCdata[sctid1];
    SCAtomData &data2 = SCdata[sctid2];

    SCInteractionData &mixer = MixingMap[sctid1][sctid2];

    RealType rcij = mixer.rCut;

    if ( idat.rij  < rcij) {
      RealType vcij = mixer.vCut; 
      RealType rhtmp, drhodr, vptmp, dvpdr;
      
      mixer.phi->getValueAndDerivativeAt( idat.rij, rhtmp, drhodr );      
      mixer.V->getValueAndDerivativeAt( idat.rij, vptmp, dvpdr);
      
      RealType pot_temp = vptmp - vcij;
      idat.vpair += pot_temp;
      
      RealType dudr = drhodr * ( idat.dfrho1 + idat.dfrho2 ) + dvpdr;
      
      idat.f1 += idat.d * dudr / idat.rij ;
        
      if (idat.doParticlePot) {
        // particlePot is the difference between the full potential and
        // the full potential without the presence of a particular
        // particle (atom1).
        //
        // This reduces the density at other particle locations, so we
        // need to recompute the density at atom2 assuming atom1 didn't
        // contribute.  This then requires recomputing the density
        // functional for atom2 as well.
        
        idat.particlePot1 -= data2.c * data2.epsilon * 
          sqrt( idat.rho2 - rhtmp) + idat.frho2;

        idat.particlePot2 -= data1.c * data1.epsilon * 
          sqrt( idat.rho1 - rhtmp) + idat.frho1;
      }
      
      idat.pot[METALLIC_PAIR_FAMILY] += pot_temp;
      
      if (idat.isSelected) 
        idat.selePot[METALLIC_PAIR_FAMILY] += pot_temp;
      
    }
      
    return;    
  }

  RealType SC::getSuggestedCutoffRadius(pair<AtomType*, AtomType*> atypes) {
    if (!initialized_) initialize();   

    int atid1 = atypes.first->getIdent();
    int atid2 = atypes.second->getIdent();
    int &sctid1 = SCtids[atid1];
    int &sctid2 = SCtids[atid2];
    
    if (sctid1 == -1 || sctid2 == -1) {
      return 0.0;
    } else {
      return MixingMap[sctid1][sctid2].rCut;
    }
  }
}
