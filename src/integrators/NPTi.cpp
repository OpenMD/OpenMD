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
 * [4]  Kuang & Gezelter,  J. Chem. Phys. 133, 164101 (2010).
 * [5]  Vardeman, Stocker & Gezelter, J. Chem. Theory Comput. 7, 834 (2011).
 */
 
#include "NPTi.hpp"
#include "brains/SimInfo.hpp"
#include "brains/Thermo.hpp"
#include "integrators/NPT.hpp"
#include "primitives/Molecule.hpp"
#include "utils/PhysicalConstants.hpp"
#include "utils/simError.h"

namespace OpenMD {

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
		   (PhysicalConstants::pressureConvert*NkBT*tb2));
    oldEta = eta;
  }

  void NPTi::evolveEtaB() {

    prevEta = eta;
    eta = oldEta + dt2 * ( instaVol * (instaPress - targetPressure) /
			   (PhysicalConstants::pressureConvert*NkBT*tb2));
  }

  void NPTi::calcVelScale() {
    vScale = thermostat.first + eta;
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
    sc  = (oldPos[index] + pos)/(RealType)2.0 -COM;
    sc *= eta;
  }

  void NPTi::scaleSimBox(){

    RealType scaleFactor;

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
      Mat3x3d hmat = snap->getHmat();
      hmat *= scaleFactor;
      snap->setHmat(hmat);
    }

  }

  bool NPTi::etaConverged() {

    return ( fabs(prevEta - eta) <= etaTolerance );
  }

  RealType NPTi::calcConservedQuantity(){

    thermostat = snap->getThermostat();
    loadEta();
    // We need NkBT a lot, so just set it here: This is the RAW number
    // of integrableObjects, so no subtraction or addition of constraints or
    // orientational degrees of freedom:
    NkBT = info_->getNGlobalIntegrableObjects()*PhysicalConstants::kB *targetTemp;

    // fkBT is used because the thermostat operates on more degrees of freedom
    // than the barostat (when there are particles with orientational degrees
    // of freedom).  
    fkBT = info_->getNdf()*PhysicalConstants::kB *targetTemp;    
    
    RealType conservedQuantity;
    RealType Energy;
    RealType thermostat_kinetic;
    RealType thermostat_potential;
    RealType barostat_kinetic;
    RealType barostat_potential;

    Energy =thermo.getTotalEnergy();

    thermostat_kinetic = fkBT* tt2 * thermostat.first * 
      thermostat.first / (2.0 * PhysicalConstants::energyConvert);

    thermostat_potential = fkBT* thermostat.second / PhysicalConstants::energyConvert;


    barostat_kinetic = 3.0 * NkBT * tb2 * eta * eta /(2.0 * PhysicalConstants::energyConvert);

    barostat_potential = (targetPressure * thermo.getVolume() / PhysicalConstants::pressureConvert) /
      PhysicalConstants::energyConvert;

    conservedQuantity = Energy + thermostat_kinetic + thermostat_potential +
      barostat_kinetic + barostat_potential;
    
    return conservedQuantity;
  }

  void NPTi::loadEta() {
    Mat3x3d etaMat = snap->getBarostat();
    eta = etaMat(0,0);
    //if (fabs(etaMat(1,1) - eta) >= OpenMD::epsilon || fabs(etaMat(1,1) - eta) >= OpenMD::epsilon || !etaMat.isDiagonal()) {
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
    snap->setBarostat(etaMat);
  }
}
