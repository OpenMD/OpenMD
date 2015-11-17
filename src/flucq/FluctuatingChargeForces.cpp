/*
 * Copyright (c) 2012 The University of Notre Dame. All Rights Reserved.
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
 
#ifdef IS_MPI
#include <mpi.h>
#endif

#include "flucq/FluctuatingChargeForces.hpp"
#include "types/FluctuatingChargeAdapter.hpp"
#include "math/RealSymmetricTridiagonal.hpp"

namespace OpenMD {

  FluctuatingChargeForces::FluctuatingChargeForces(SimInfo* info) : 
    info_(info), initialized_(false) {
  }

  void FluctuatingChargeForces::initialize(){
    FQtypes.clear();
    FQtids.clear();
    FQtids.resize( forceField_->getNAtomType(), -1);

    set<AtomType*>::iterator at;
    for (at = simTypes_.begin(); at != simTypes_.end(); ++at) 
      if ((*at)->isFluctuatingCharge()) addType(*at);
    
    initialized_ = true;
  }

  void FluctuatingChargeForces::getSelfInteraction(int atid, RealType charge, 
                                                   RealType &potential, 
                                                   RealType &force) {
    if (!initialized_) initialize();

    data = FQMap[FQtids[atid]];

    if (data.hasMultipleMinima) {
      int nDiabats = data.diabaticStates.size();
      RealType k = data.curvature;
      RealType c = data.coupling;

      if (nDiabats == 2){
        RealType q1 = data.diabaticStates[0].first;
        RealType e1 = data.diabaticStates[0].second;
        RealType q2 = data.diabaticStates[1].first;
        RealType e2 = data.diabaticStates[1].second;

        RealType v1 = e1 + 0.5 * k * (charge-q1)*(charge-q1);
        RealType v1p = k*(charge-q1);
        RealType v2 = e2 + 0.5 * k * (charge-q2)*(charge-q2);
        RealType v2p = k*(charge-q2);
        potential += (v1 + v2 - sqrt(4*c*c+ v1*v1 - 2.0*v1*v2 + v2*v2))/2.0;
        force -= (v1p + v2p - ((v1 - v2) * (v1p - v2p)) /
                  sqrt(4*c*c + v1*v1 - 2*v1*v2 + v2*v2))/2.0;
      } else {
        DynamicVector<RealType> diagonals(nDiabats);
        DynamicVector<RealType> subdiagonals(nDiabats, c);
        DynamicVector<RealType> vp(nDiabats);
        DynamicVector<RealType> eigenvalues(nDiabats);
        DynamicRectMatrix<RealType> eigenvectors(nDiabats, nDiabats);
        RealType q;
        RealType e; 

        for (int i = 0; i < nDiabats; i++) {
          q = data.diabaticStates[i].first;
          e = data.diabaticStates[i].second;
          diagonals(i) = e + 0.5 * k * (charge-q)*(charge-q);
          vp(i) = k*(charge-q);
        }

        RealSymmetricTridiagonal<RealType> eig(diagonals, subdiagonals);
        eig.getEigenvalues( eigenvalues );
        eig.getEigenvectors( eigenvectors );

        potential += eigenvalues(0);
        for (int i = 0; i < nDiabats; i++) {
          // all of the couplings are constant, so the condon
          // approximation and Hellmann-Feynman let us obtain the
          // forces easily:
          force -= pow(eigenvectors(0, i), 2) * vp(i);
        }              
      }
    } else {
      RealType Jii = data.hardness;
      RealType chi = data.electronegativity;
      force -=  charge * Jii + chi;
      potential += charge * (charge * Jii * 0.5 + chi);
    }
    return;
  }

  void FluctuatingChargeForces::addType(AtomType* atomType) {
    FluctuatingChargeAtomData data;
    FluctuatingChargeAdapter fqa = FluctuatingChargeAdapter(atomType);
    if (fqa.isFluctuatingCharge()) {
      if (fqa.hasMultipleMinima()) {
        data.hasMultipleMinima = true;      
        data.curvature = fqa.getCurvature();
        data.coupling = fqa.getCoupling();
        data.diabaticStates = fqa.getDiabaticStates();
      } else {
        data.hasMultipleMinima = false;      
        data.electronegativity = fqa.getElectronegativity();
        data.hardness = fqa.getHardness();
        data.slaterN = fqa.getSlaterN();
        data.slaterZeta = fqa.getSlaterZeta();
      }
      int atid = atomType->getIdent();
      int fqtid = FQtypes.size();

      pair<set<int>::iterator,bool> ret;    
      ret = FQtypes.insert( atid );
      if (ret.second == false) {
        sprintf( painCave.errMsg,
                 "FluctuatingChargeForces already had a previous fluctuating "
                 "charge entry with ident %d\n",
                 atid );
        painCave.severity = OPENMD_INFO;
        painCave.isFatal = 0;
        simError();         
      }
      FQtids[atid] = fqtid;
      FQMap.push_back(data);
    }
  }
}
