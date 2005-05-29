/*
 * Copyright (c) 2005 The University of Notre Dame. All Rights Reserved.
 *
 * The University of Notre Dame grants you ("Licensee") a
 * non-exclusive, royalty free, license to use, modify and
 * redistribute this software in source and binary code form, provided
 * that the following conditions are met:
 *
 * 1. Acknowledgement of the program authors must be made in any
 *    publication of scientific results based in part on use of the
 *    program.  An acceptable form of acknowledgement is citation of
 *    the article in which the program was described (Matthew
 *    A. Meineke, Charles F. Vardeman II, Teng Lin, Christopher
 *    J. Fennell and J. Daniel Gezelter, "OOPSE: An Object-Oriented
 *    Parallel Simulation Engine for Molecular Dynamics,"
 *    J. Comput. Chem. 26, pp. 252-271 (2005))
 *
 * 2. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 *
 * 3. Redistributions in binary form must reproduce the above copyright
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
 */
 
#include "NPTi.hpp"
#include "brains/SimInfo.hpp"
#include "brains/Thermo.hpp"
#include "integrators/NPT.hpp"
#include "primitives/Molecule.hpp"
#include "utils/OOPSEConstant.hpp"
#include "utils/simError.h"

namespace oopse {

  // Basic isotropic thermostating and barostating via the Melchionna
  // modification of the Hoover algorithm:
  //
  //    Melchionna, S., Ciccotti, G., and Holian, B. L., 1993,
  //       Molec. Phys., 78, 533.
  //
  //           and
  //
  //    Hoover, W. G., 1986, Phys. Rev. A, 34, 2499.

  NPTi::NPTi ( SimInfo *info) : NPT(info){

  }

  void NPTi::evolveEtaA() {
    eta += dt2 * ( instaVol * (instaPress - targetPressure) /
		   (OOPSEConstant::pressureConvert*NkBT*tb2));
    oldEta = eta;
  }

  void NPTi::evolveEtaB() {

    prevEta = eta;
    eta = oldEta + dt2 * ( instaVol * (instaPress - targetPressure) /
			   (OOPSEConstant::pressureConvert*NkBT*tb2));
  }

  void NPTi::calcVelScale() {
    vScale = chi + eta;
  }

  void NPTi::getVelScaleA(Vector3d& sc, const Vector3d& vel) {
    sc = vel * vScale;
  }

  void NPTi::getVelScaleB(Vector3d& sc, int index ){
    sc = oldVel[index] * vScale;    
  }


  void NPTi::getPosScale(const Vector3d& pos, const Vector3d& COM,
			 int index, Vector3d& sc){
    /**@todo*/
    sc  = (oldPos[index] + pos)/2.0 -COM;
    sc *= eta;
  }

  void NPTi::scaleSimBox(){

    double scaleFactor;

    scaleFactor = exp(dt*eta);

    if ((scaleFactor > 1.1) || (scaleFactor < 0.9)) {
      sprintf( painCave.errMsg,
	       "NPTi error: Attempting a Box scaling of more than 10 percent"
	       " check your tauBarostat, as it is probably too small!\n"
	       " eta = %lf, scaleFactor = %lf\n", eta, scaleFactor
	       );
      painCave.isFatal = 1;
      simError();
    } else {
      Mat3x3d hmat = currentSnapshot_->getHmat();
      hmat *= scaleFactor;
      currentSnapshot_->setHmat(hmat);
    }

  }

  bool NPTi::etaConverged() {

    return ( fabs(prevEta - eta) <= etaTolerance );
  }

  double NPTi::calcConservedQuantity(){

    chi= currentSnapshot_->getChi();
    integralOfChidt = currentSnapshot_->getIntegralOfChiDt();
    loadEta();
    // We need NkBT a lot, so just set it here: This is the RAW number
    // of integrableObjects, so no subtraction or addition of constraints or
    // orientational degrees of freedom:
    NkBT = info_->getNGlobalIntegrableObjects()*OOPSEConstant::kB *targetTemp;

    // fkBT is used because the thermostat operates on more degrees of freedom
    // than the barostat (when there are particles with orientational degrees
    // of freedom).  
    fkBT = info_->getNdf()*OOPSEConstant::kB *targetTemp;    
    
    double conservedQuantity;
    double Energy;
    double thermostat_kinetic;
    double thermostat_potential;
    double barostat_kinetic;
    double barostat_potential;

    Energy =thermo.getTotalE();

    thermostat_kinetic = fkBT* tt2 * chi * chi / (2.0 * OOPSEConstant::energyConvert);

    thermostat_potential = fkBT* integralOfChidt / OOPSEConstant::energyConvert;


    barostat_kinetic = 3.0 * NkBT * tb2 * eta * eta /(2.0 * OOPSEConstant::energyConvert);

    barostat_potential = (targetPressure * thermo.getVolume() / OOPSEConstant::pressureConvert) /
      OOPSEConstant::energyConvert;

    conservedQuantity = Energy + thermostat_kinetic + thermostat_potential +
      barostat_kinetic + barostat_potential;
    
    return conservedQuantity;
  }

  void NPTi::loadEta() {
    Mat3x3d etaMat = currentSnapshot_->getEta();
    eta = etaMat(0,0);
    //if (fabs(etaMat(1,1) - eta) >= oopse::epsilon || fabs(etaMat(1,1) - eta) >= oopse::epsilon || !etaMat.isDiagonal()) {
    //    sprintf( painCave.errMsg,
    //             "NPTi error: the diagonal elements of  eta matrix are not the same or etaMat is not a diagonal matrix");
    //    painCave.isFatal = 1;
    //    simError();
    //}
  }

  void NPTi::saveEta() {
    Mat3x3d etaMat(0.0);
    etaMat(0, 0) = eta;
    etaMat(1, 1) = eta;
    etaMat(2, 2) = eta;
    currentSnapshot_->setEta(etaMat);
  }
}
