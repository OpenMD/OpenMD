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

#include "nonbonded/Sticky.hpp"
#include "nonbonded/LJ.hpp"
#include "types/StickyAdapter.hpp"
#include "utils/simError.h"

using namespace std;
namespace OpenMD {
  
  Sticky::Sticky() : initialized_(false), forceField_(NULL), name_("Sticky") {}
    
  void Sticky::initialize() {    
    Stypes.clear();
    Stids.clear();
    MixingMap.clear();
    nSticky_=0;

    Stids.resize( forceField_->getNAtomType(), -1);

    // Sticky handles all of the Sticky-Sticky interactions

    set<AtomType*>::iterator at;
    for (at = simTypes_.begin(); at != simTypes_.end(); ++at) {
      if ((*at)->isSticky()) nSticky_++;
    }

    MixingMap.resize(nSticky_);

    for (at = simTypes_.begin(); at != simTypes_.end(); ++at) {
      if ((*at)->isSticky()) addType( *at );
    }
    
    initialized_ = true;
  }
      
  void Sticky::addType(AtomType* atomType){
    StickyAdapter sticky1 = StickyAdapter(atomType);

    // add it to the map:
    
    int atid = atomType->getIdent();
    int stid = Stypes.size();
  
    pair<set<int>::iterator,bool> ret;    
    ret = Stypes.insert( atid );
    if (ret.second == false) {
      sprintf( painCave.errMsg,
               "Sticky already had a previous entry with ident %d\n",
               atid) ;
      painCave.severity = OPENMD_INFO;
      painCave.isFatal = 0;
      simError();         
    }

    Stids[atid] = stid;
    MixingMap[stid].resize( nSticky_ );    


    // Now, iterate over all known types and add to the mixing map:
    
    std::set<int>::iterator it; 
    for( it = Stypes.begin(); it != Stypes.end(); ++it) {
      

      int stid2 = Stids[ (*it) ];
      AtomType* atype2 = forceField_->getAtomType( (*it) );
      StickyAdapter sticky2 = StickyAdapter(atype2);

      StickyInteractionData mixer;        
           
      // Mixing two different sticky types is silly, but if you want 2
      // sticky types in your simulation, we'll let you do it with the
      // Lorentz- Berthelot mixing rules (which happen to do the right thing 
      // when atomType and atype2 happen to be the same.
      
      mixer.rl   = 0.5 * ( sticky1.getRl() + sticky2.getRl() );
      mixer.ru   = 0.5 * ( sticky1.getRu() + sticky2.getRu() );
      mixer.rlp  = 0.5 * ( sticky1.getRlp() + sticky2.getRlp() );
      mixer.rup  = 0.5 * ( sticky1.getRup() + sticky2.getRup() );
      mixer.rbig = max(mixer.ru, mixer.rup);
      mixer.w0  = sqrt( sticky1.getW0()   * sticky2.getW0()  );
      mixer.v0  = sqrt( sticky1.getV0()   * sticky2.getV0()  );
      mixer.v0p = sqrt( sticky1.getV0p()  * sticky2.getV0p() );
      mixer.isPower = sticky1.isStickyPower() && sticky2.isStickyPower();

      CubicSpline* s = new CubicSpline();
      s->addPoint(mixer.rl, 1.0);
      s->addPoint(mixer.ru, 0.0);
      mixer.s = s;

      CubicSpline* sp = new CubicSpline();
      sp->addPoint(mixer.rlp, 1.0);
      sp->addPoint(mixer.rup, 0.0);
      mixer.sp = sp;

      MixingMap[stid2].resize( nSticky_ );

      MixingMap[stid][stid2] = mixer;
      if (stid2 != stid) {
        MixingMap[stid2][stid] = mixer;
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
    
    StickyInteractionData &mixer = MixingMap[Stids[idat.atid1]][Stids[idat.atid2]];
    
    RealType w0  = mixer.w0;
    RealType v0  = mixer.v0; 
    RealType v0p = mixer.v0p;
    RealType rl  = mixer.rl; 
    RealType ru  = mixer.ru; 
    RealType rlp = mixer.rlp;
    RealType rup = mixer.rup;
    RealType rbig = mixer.rbig;
    bool isPower = mixer.isPower;
    
    if ( idat.rij <= rbig) {
      
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
      
      if ( idat.rij < ru) {
        if ( idat.rij < rl) {
          s = 1.0;
          dsdr = 0.0;
        } else {          
          // we are in the switching region           
          mixer.s->getValueAndDerivativeAt(idat.rij, s, dsdr);
        }
      }
      
      if (idat.rij < rup) {
        if ( idat.rij < rlp) {
          sp = 1.0;
          dspdr = 0.0;
        } else {
          // we are in the switching region           
          mixer.sp->getValueAndDerivativeAt( idat.rij, sp, dspdr);
        }
      }
      
      
      RealType wi = 2.0*(xi2-yi2)*zi / r3;
      RealType wj = 2.0*(xj2-yj2)*zj / r3;
      RealType w = wi+wj;
      
      
      RealType zif = zi/ idat.rij  - 0.6;
      RealType zis = zi/ idat.rij  + 0.8;
      
      RealType zjf = zj/ idat.rij  - 0.6;
      RealType zjs = zj/ idat.rij  + 0.8;
      
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
                    2.0*(1.0/ idat.rij  - zi2/r3)*uglyi);
      
      Vector3d dwjp(-2.0*xj*zj*uglyj/r3,
                    -2.0*yj*zj*uglyj/r3,
                    2.0*(1.0/ idat.rij  - zj2/r3)*uglyj);
      
      Vector3d dwidu(4.0*(yi*zi2 + 0.5*yi*(xi2-yi2))/r3,
                     4.0*(xi*zi2 - 0.5*xi*(xi2-yi2))/r3,
                     - 8.0*xi*yi*zi/r3);
      
      Vector3d dwjdu(4.0*(yj*zj2 + 0.5*yj*(xj2-yj2))/r3,
                     4.0*(xj*zj2 - 0.5*xj*(xj2-yj2))/r3,
                     - 8.0*xj*yj*zj/r3);
      
      Vector3d dwipdu(2.0*yi*uglyi/ idat.rij ,
                      -2.0*xi*uglyi/ idat.rij ,
                      0.0);
      
      Vector3d dwjpdu(2.0*yj*uglyj/ idat.rij ,
                      -2.0*xj*uglyj/ idat.rij ,
                      0.0);
      
      if (isPower) {
        cerr << "This is probably an error!\n";
        RealType frac1 = 0.25;
        RealType frac2 = 0.75;      
        RealType wi2 = wi*wi;
        RealType wj2 = wj*wj;
        // sticky power has no w' function:
        w = frac1 * wi * wi2 + frac2*wi + frac1*wj*wj2 + frac2*wj + v0p; 
        wp = 0.0;
        dwi = frac1*RealType(3.0)*wi2*dwi + frac2*dwi;
        dwj = frac1*RealType(3.0)*wj2*dwi + frac2*dwi;
        dwip = V3Zero;
        dwjp = V3Zero;
        dwidu = frac1*RealType(3.0)*wi2*dwidu + frac2*dwidu;
        dwidu = frac1*RealType(3.0)*wj2*dwjdu + frac2*dwjdu;
        dwipdu = V3Zero;
        dwjpdu = V3Zero;
        sp = 0.0;
        dspdr = 0.0;
      }
      
      idat.vpair += 0.5 * (v0*s*w + v0p*sp*wp);
      idat.pot[HYDROGENBONDING_FAMILY] += 0.5 * (v0*s*w + v0p*sp*wp)* idat.sw;
      if (idat.isSelected)
        idat.selePot[HYDROGENBONDING_FAMILY] += 0.5 * (v0*s*w +
						       v0p*sp*wp)* idat.sw;
      
      // do the torques first since they are easy:
      // remember that these are still in the body-fixed axes
      
      Vector3d ti = 0.5 * idat.sw *(v0*s*dwidu + v0p*sp*dwipdu);
      Vector3d tj = 0.5 * idat.sw *(v0*s*dwjdu + v0p*sp*dwjpdu);
      
      // go back to lab frame using transpose of rotation matrix:
      
      idat.t1 += A1trans * ti;
      idat.t2 += A2trans * tj;
      
      // Now, on to the forces:
      
      // first rotate the i terms back into the lab frame:
      
      Vector3d radcomi = (v0 * s * dwi + v0p * sp * dwip) *  idat.sw;
      Vector3d radcomj = (v0 * s * dwj + v0p * sp * dwjp) *  idat.sw;
      
      Vector3d fii = A1trans * radcomi;
      Vector3d fjj = A2trans * radcomj;
      
      // now assemble these with the radial-only terms:
      
      idat.f1 += 0.5 * ((v0*dsdr*w + v0p*dspdr*wp) * idat.d /
			idat.rij  + fii - fjj);
      
    }
    
    return;      
  }
  
  RealType Sticky::getSuggestedCutoffRadius(pair<AtomType*, AtomType*> atypes) {
    if (!initialized_) initialize();   
    int atid1 = atypes.first->getIdent();
    int atid2 = atypes.second->getIdent();
    int stid1 = Stids[atid1];
    int stid2 = Stids[atid2];
    
    if (stid1 == -1 || stid2 == -1) return 0.0;
    else {      
      return MixingMap[stid1][stid2].rbig;
    }
  }
}
