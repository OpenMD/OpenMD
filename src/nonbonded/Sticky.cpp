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
#include "nonbonded/Sticky.hpp"
#include "nonbonded/LJ.hpp"
#include "utils/simError.h"

using namespace std;
namespace OpenMD {
  
  Sticky::Sticky() : name_("Sticky"), initialized_(false), forceField_(NULL) {}
  
  StickyParam Sticky::getStickyParam(AtomType* atomType) {
    
    // Do sanity checking on the AtomType we were passed before
    // building any data structures:
    if (!atomType->isSticky() && !atomType->isStickyPower()) {
      sprintf( painCave.errMsg,
               "Sticky::getStickyParam was passed an atomType (%s) that does\n"
               "\tnot appear to be a Sticky atom.\n",
               atomType->getName().c_str());
      painCave.severity = OPENMD_ERROR;
      painCave.isFatal = 1;
      simError();
    }
    
    DirectionalAtomType* daType = dynamic_cast<DirectionalAtomType*>(atomType);
    GenericData* data = daType->getPropertyByName("Sticky");
    if (data == NULL) {
      sprintf( painCave.errMsg, "Sticky::getStickyParam could not find\n"
               "\tSticky parameters for atomType %s.\n", 
               daType->getName().c_str());
      painCave.severity = OPENMD_ERROR;
      painCave.isFatal = 1;
      simError(); 
    }
    
    StickyParamGenericData* stickyData = dynamic_cast<StickyParamGenericData*>(data);
    if (stickyData == NULL) {
      sprintf( painCave.errMsg,
               "Sticky::getStickyParam could not convert GenericData to\n"
               "\tStickyParamGenericData for atom type %s\n", 
               daType->getName().c_str());
      painCave.severity = OPENMD_ERROR;
      painCave.isFatal = 1;
      simError();          
    }
    
    return stickyData->getData();
  }
  
  void Sticky::initialize() {    
    
    ForceFieldOptions& fopts = forceField_->getForceFieldOptions();
    ForceField::AtomTypeContainer* atomTypes = forceField_->getAtomTypes();
    ForceField::AtomTypeContainer::MapTypeIterator i;
    AtomType* at;

    // Sticky handles all of the Sticky-Sticky interactions

    for (at = atomTypes->beginType(i); at != NULL;
         at = atomTypes->nextType(i)) {
      
      if (at->isSticky() || at->isStickyPower())
        addType(at);
    }
    
    initialized_ = true;
  }
      
  void Sticky::addType(AtomType* atomType){
    // add it to the map:
    AtomTypeProperties atp = atomType->getATP();    
    
    pair<map<int,AtomType*>::iterator,bool> ret;    
    ret = StickyMap.insert( pair<int, AtomType*>(atp.ident, atomType) );
    if (ret.second == false) {
      sprintf( painCave.errMsg,
               "Sticky already had a previous entry with ident %d\n",
               atp.ident);
      painCave.severity = OPENMD_INFO;
      painCave.isFatal = 0;
      simError();         
    }
    
    RealType w0i, v0i, v0pi, rli, rui, rlpi, rupi;
    
    StickyParam sticky1 = getStickyParam(atomType);

    // Now, iterate over all known types and add to the mixing map:
    
    map<int, AtomType*>::iterator it;
    for( it = StickyMap.begin(); it != StickyMap.end(); ++it) {
      
      AtomType* atype2 = (*it).second;
    
      StickyParam sticky2 = getStickyParam(atype2);

      StickyInteractionData mixer;        
           
      // Mixing two different sticky types is silly, but if you want 2
      // sticky types in your simulation, we'll let you do it with the
      // Lorentz- Berthelot mixing rules (which happen to do the right thing 
      // when atomType and atype2 happen to be the same.
      
      mixer.rl   = 0.5 * ( sticky1.rl + sticky2.rl );
      mixer.ru   = 0.5 * ( sticky1.ru + sticky2.ru );
      mixer.rlp  = 0.5 * ( sticky1.rlp + sticky2.rlp );
      mixer.rup  = 0.5 * ( sticky1.rup + sticky2.rup );
      mixer.rbig = max(mixer.ru, mixer.rup);
      mixer.w0  = sqrt( sticky1.w0   * sticky2.w0  );
      mixer.v0  = sqrt( sticky1.v0   * sticky2.v0  );
      mixer.v0p = sqrt( sticky1.v0p  * sticky2.v0p );
      mixer.isPower = atomType->isStickyPower() && atype2->isStickyPower();

      CubicSpline* s = new CubicSpline();
      s->addPoint(mixer.rl, 1.0);
      s->addPoint(mixer.ru, 0.0);
      mixer.s = s;

      CubicSpline* sp = new CubicSpline();
      sp->addPoint(mixer.rlp, 1.0);
      sp->addPoint(mixer.rup, 0.0);
      mixer.sp = sp;

      
      pair<AtomType*, AtomType*> key1, key2;
      key1 = make_pair(atomType, atype2);
      key2 = make_pair(atype2, atomType);
      
      MixingMap[key1] = mixer;
      if (key2 != key1) {
        MixingMap[key2] = mixer;
      }
    }
  }

  /**
   * This function does the sticky portion of the SSD potential
   * [Chandra and Ichiye, Journal of Chemical Physics 111, 2701
   * (1999)].  The Lennard-Jones and dipolar interaction must be
   * handled separately.  We assume that the rotation matrices have
   * already been calculated and placed in the A1 & A2 entries in the
   * idat structure.
   */
  
  void Sticky::calcForce(InteractionData &idat) {
   
    if (!initialized_) initialize();
    
    pair<AtomType*, AtomType*> key = make_pair(idat.atype1, idat.atype2);
    map<pair<AtomType*, AtomType*>, StickyInteractionData>::iterator it;
    it = MixingMap.find(key);
    if (it != MixingMap.end()) {

      StickyInteractionData mixer = (*it).second;
      
      RealType w0  = mixer.w0;
      RealType v0  = mixer.v0; 
      RealType v0p = mixer.v0p;
      RealType rl  = mixer.rl; 
      RealType ru  = mixer.ru; 
      RealType rlp = mixer.rlp;
      RealType rup = mixer.rup;
      RealType rbig = mixer.rbig;
      bool isPower = mixer.isPower;
      
      if (idat.rij <= rbig) {
        
        RealType r3 = idat.r2 * idat.rij;
        RealType r5 = r3 * idat.r2;
        
        RotMat3x3d A1trans = idat.A1.transpose();
        RotMat3x3d A2trans = idat.A2.transpose();
        
        // rotate the inter-particle separation into the two different
        // body-fixed coordinate systems:
        
        Vector3d ri = idat.A1 * idat.d;
        
        // negative sign because this is the vector from j to i:
        
        Vector3d rj = - idat.A2 * idat.d;
        
        RealType xi = ri.x();
        RealType yi = ri.y();
        RealType zi = ri.z();
        
        RealType xj = rj.x();
        RealType yj = rj.y();
        RealType zj = rj.z();
        
        RealType xi2 = xi * xi;
        RealType yi2 = yi * yi;
        RealType zi2 = zi * zi;
        
        RealType xj2 = xj * xj;
        RealType yj2 = yj * yj;
        RealType zj2 = zj * zj;     
        
        // calculate the switching info. from the splines
        
        RealType s = 0.0;
        RealType dsdr = 0.0;
        RealType sp = 0.0;
        RealType dspdr = 0.0;
        
        if (idat.rij < ru) {
          if (idat.rij < rl) {
            s = 1.0;
            dsdr = 0.0;
          } else {          
            // we are in the switching region 
            
            pair<RealType, RealType> res = mixer.s->getValueAndDerivativeAt(idat.rij);
            s = res.first;
            dsdr = res.second;
          }
        }
        
        if (idat.rij < rup) {
          if (idat.rij < rlp) {
            sp = 1.0;
            dspdr = 0.0;
          } else {
            // we are in the switching region 
            
            pair<RealType, RealType> res =mixer.sp->getValueAndDerivativeAt(idat.rij);
            sp = res.first;
            dspdr = res.second;
          }
        }
        
        RealType wi = 2.0*(xi2-yi2)*zi / r3;
        RealType wj = 2.0*(xj2-yj2)*zj / r3;
        RealType w = wi+wj;
        
        
        RealType zif = zi/idat.rij - 0.6;
        RealType zis = zi/idat.rij + 0.8;
        
        RealType zjf = zj/idat.rij - 0.6;
        RealType zjs = zj/idat.rij + 0.8;
        
        RealType wip = zif*zif*zis*zis - w0;
        RealType wjp = zjf*zjf*zjs*zjs - w0;
        RealType wp = wip + wjp;
        
        Vector3d dwi(4.0*xi*zi/r3  - 6.0*xi*zi*(xi2-yi2)/r5, 
                     - 4.0*yi*zi/r3  - 6.0*yi*zi*(xi2-yi2)/r5,
                     2.0*(xi2-yi2)/r3  - 6.0*zi2*(xi2-yi2)/r5);
        
        Vector3d dwj(4.0*xj*zj/r3  - 6.0*xj*zj*(xj2-yj2)/r5,
                     - 4.0*yj*zj/r3  - 6.0*yj*zj*(xj2-yj2)/r5,
                     2.0*(xj2-yj2)/r3  - 6.0*zj2*(xj2-yj2)/r5);
        
        RealType uglyi = zif*zif*zis + zif*zis*zis;
        RealType uglyj = zjf*zjf*zjs + zjf*zjs*zjs;

        Vector3d dwip(-2.0*xi*zi*uglyi/r3, 
                      -2.0*yi*zi*uglyi/r3,
                      2.0*(1.0/idat.rij - zi2/r3)*uglyi);
        
        Vector3d dwjp(-2.0*xj*zj*uglyj/r3,
                      -2.0*yj*zj*uglyj/r3,
                      2.0*(1.0/idat.rij - zj2/r3)*uglyj);
        
        Vector3d dwidu(4.0*(yi*zi2 + 0.5*yi*(xi2-yi2))/r3,
                       4.0*(xi*zi2 - 0.5*xi*(xi2-yi2))/r3,
                       - 8.0*xi*yi*zi/r3);
        
        Vector3d dwjdu(4.0*(yj*zj2 + 0.5*yj*(xj2-yj2))/r3,
                       4.0*(xj*zj2 - 0.5*xj*(xj2-yj2))/r3,
                       - 8.0*xj*yj*zj/r3);
        
        Vector3d dwipdu(2.0*yi*uglyi/idat.rij,
                        -2.0*xi*uglyi/idat.rij,
                        0.0);
        
        Vector3d dwjpdu(2.0*yj*uglyj/idat.rij,
                        -2.0*xj*uglyj/idat.rij,
                        0.0);
        
        if (isPower) {
          RealType frac1 = 0.25;
          RealType frac2 = 0.75;      
          RealType wi2 = wi*wi;
          RealType wj2 = wj*wj;
          // sticky power has no w' function:
          w = frac1 * wi * wi2 + frac2*wi + frac1*wj*wj2 + frac2*wj + v0p; 
          wp = 0.0;
          dwi = frac1*3.0*wi2*dwi + frac2*dwi;
          dwj = frac1*3.0*wj2*dwi + frac2*dwi;
          dwip = V3Zero;
          dwjp = V3Zero;
          dwidu = frac1*3.0*wi2*dwidu + frac2*dwidu;
          dwidu = frac1*3.0*wj2*dwjdu + frac2*dwjdu;
          dwipdu = V3Zero;
          dwjpdu = V3Zero;
          sp = 0.0;
          dspdr = 0.0;
        }
        
        idat.vpair[2] += 0.5*(v0*s*w + v0p*sp*wp);
        idat.pot[2] += 0.5*(v0*s*w + v0p*sp*wp)*idat.sw;
        
        // do the torques first since they are easy:
        // remember that these are still in the body-fixed axes
        
        Vector3d ti = 0.5*idat.sw*(v0*s*dwidu + v0p*sp*dwipdu);
        Vector3d tj = 0.5*idat.sw*(v0*s*dwjdu + v0p*sp*dwjpdu);
        
        // go back to lab frame using transpose of rotation matrix:
        
        idat.t1 += A1trans * ti;
        idat.t2 += A2trans * tj;
        
        // Now, on to the forces:
        
        // first rotate the i terms back into the lab frame:
        
        Vector3d radcomi = (v0 * s * dwi + v0p * sp * dwip) * idat.sw;
        Vector3d radcomj = (v0 * s * dwj + v0p * sp * dwjp) * idat.sw;
        
        Vector3d fii = A1trans * radcomi;
        Vector3d fjj = A2trans * radcomj;
        
        // now assemble these with the radial-only terms:
        
        idat.f1 += 0.5 * ((v0*dsdr*w + v0p*dspdr*wp) * idat.d / 
                          idat.rij + fii - fjj);
        
      }
    }
    
    return;      
  }

  RealType Sticky::getSuggestedCutoffRadius(AtomType* at1, AtomType* at2) {
    if (!initialized_) initialize();   
    pair<AtomType*, AtomType*> key = make_pair(at1, at2); 
    map<pair<AtomType*, AtomType*>, StickyInteractionData>::iterator it;
    it = MixingMap.find(key);
    if (it == MixingMap.end()) 
      return 0.0;
    else  {
      StickyInteractionData mixer = (*it).second;
      return mixer.rbig;
    }
  }
}
