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
 
#include "brains/SimInfo.hpp"
#include "brains/Thermo.hpp"
#include "integrators/IntegratorCreator.hpp"
#include "integrators/NPrT.hpp"
#include "primitives/Molecule.hpp"
#include "utils/Constants.hpp"
#include "utils/simError.h"

namespace OpenMD {
  NPrT::NPrT(SimInfo* info) : NPT(info) {
    Globals* simParams = info_->getSimParams();
    if (!simParams->haveSurfaceTension()) {
      sprintf(painCave.errMsg,
              "If you use the NPT integrator, you must set tauBarostat.\n");
      painCave.severity = OPENMD_ERROR;
      painCave.isFatal = 1;
      simError();
    } else {
      surfaceTension_= simParams->getSurfaceTension()* Constants::surfaceTensionConvert * Constants::energyConvert;

      // Default value of privilegedAxis is "z"
      if (simParams->getPrivilegedAxis() == "x")
	axis_ = 0;
      else if (simParams->getPrivilegedAxis() == "y")
	axis_ = 1;
      else if (simParams->getPrivilegedAxis() == "z")
	axis_ = 2;
	     
      // Compute complementary axes to the privileged axis
      axis1_ = (axis_ + 1) % 3;
      axis2_ = (axis_ + 2) % 3;
    }

  }
  void NPrT::evolveEtaA() {
    Mat3x3d hmat = snap->getHmat();
    RealType hz = hmat(axis_, axis_);
    RealType Axy = hmat(axis1_,axis1_) * hmat(axis2_, axis2_);
    RealType sx = -hz * (press(axis1_, axis1_) - targetPressure/Constants::pressureConvert);
    RealType sy = -hz * (press(axis2_, axis2_) - targetPressure/Constants::pressureConvert);
    eta(axis1_,axis1_) -= dt2* Axy * (sx - surfaceTension_) / (NkBT*tb2);
    eta(axis2_,axis2_) -= dt2* Axy * (sy - surfaceTension_) / (NkBT*tb2);
    eta(axis_,axis_) += dt2 *  instaVol * (press(axis_, axis_) - targetPressure/Constants::pressureConvert) / (NkBT*tb2);
    oldEta_ = eta;  
  }

  void NPrT::evolveEtaB() {
    Mat3x3d hmat = snap->getHmat();
    RealType hz = hmat(axis_, axis_);
    RealType Axy = hmat(axis1_,axis1_) * hmat(axis2_, axis2_);
    prevEta_ = eta;
    RealType sx = -hz * (press(axis1_, axis1_) - targetPressure/Constants::pressureConvert);
    RealType sy = -hz * (press(axis2_, axis2_) - targetPressure/Constants::pressureConvert);
    eta(axis1_,axis1_) = oldEta_(axis1_, axis1_) - dt2 * Axy * (sx -surfaceTension_) / (NkBT*tb2);
    eta(axis2_,axis2_) = oldEta_(axis2_, axis2_) - dt2 * Axy * (sy -surfaceTension_) / (NkBT*tb2);
    eta(axis_,axis_) = oldEta_(axis_, axis_) + dt2 *  instaVol *
	    (press(axis_, axis_) - targetPressure/Constants::pressureConvert) / (NkBT*tb2);
  }

  void NPrT::calcVelScale(){

    for (int i = 0; i < 3; i++ ) {
      for (int j = 0; j < 3; j++ ) {
	vScale_(i, j) = eta(i, j);

	if (i == j) {
	  vScale_(i, j) += thermostat.first;
	}
      }
    }
  }

  void NPrT::getVelScaleA(Vector3d& sc, const Vector3d& vel){
    sc = vScale_ * vel;
  }

  void NPrT::getVelScaleB(Vector3d& sc, int index ) {
    sc = vScale_ * oldVel[index];
  }

  void NPrT::getPosScale(const Vector3d& pos, const Vector3d& COM, int index, Vector3d& sc) {

    /**@todo */
    Vector3d rj = (oldPos[index] + pos)/(RealType)2.0 -COM;
    sc = eta * rj;
  }

  void NPrT::scaleSimBox(){
    Mat3x3d scaleMat;

    scaleMat(axis1_, axis1_) = exp(dt*eta(axis1_, axis1_));
    scaleMat(axis2_, axis2_) = exp(dt*eta(axis2_, axis2_));    
    scaleMat(axis_, axis_) = exp(dt*eta(axis_, axis_));
    Mat3x3d hmat = snap->getHmat();
    hmat = hmat *scaleMat;
    snap->setHmat(hmat);

  }

  bool NPrT::etaConverged() {
    int i;
    RealType diffEta, sumEta;

    sumEta = 0;
    for(i = 0; i < 3; i++) {
      sumEta += pow(prevEta_(i, i) - eta(i, i), 2);
    }
    
    diffEta = sqrt( sumEta / 3.0 );

    return ( diffEta <= etaTolerance );
  }

  RealType NPrT::calcConservedQuantity(){
    thermostat = snap->getThermostat();
    loadEta();
    
    // We need NkBT a lot, so just set it here: This is the RAW number
    // of integrableObjects, so no subtraction or addition of constraints or
    // orientational degrees of freedom:
    NkBT = info_->getNGlobalIntegrableObjects()*Constants::kB *targetTemp;

    // fkBT is used because the thermostat operates on more degrees of freedom
    // than the barostat (when there are particles with orientational degrees
    // of freedom).  
    fkBT = info_->getNdf()*Constants::kB *targetTemp;    
    

    RealType totalEnergy = thermo.getTotalEnergy();

    RealType thermostat_kinetic = fkBT * tt2 * thermostat.first * thermostat.first /(2.0 * Constants::energyConvert);

    RealType thermostat_potential = fkBT* thermostat.second / Constants::energyConvert;

    SquareMatrix<RealType, 3> tmp = eta.transpose() * eta;
    RealType trEta = tmp.trace();
    
    RealType barostat_kinetic = NkBT * tb2 * trEta /(2.0 * Constants::energyConvert);

    RealType barostat_potential = (targetPressure * thermo.getVolume() / Constants::pressureConvert) /Constants::energyConvert;

    Mat3x3d hmat = snap->getHmat();
    RealType area = hmat(axis1_,axis1_) * hmat(axis2_, axis2_);

    RealType conservedQuantity = totalEnergy + thermostat_kinetic
      + thermostat_potential + barostat_kinetic + barostat_potential
      - surfaceTension_ * area/ Constants::energyConvert;

    return conservedQuantity;

  }

  void NPrT::loadEta() {
    eta= snap->getBarostat();

    //if (!eta.isDiagonal()) {
    //    sprintf( painCave.errMsg,
    //             "NPrT error: the diagonal elements of eta matrix are not the same or etaMat is not a diagonal matrix");
    //    painCave.isFatal = 1;
    //    simError();
    //}
  }

  void NPrT::saveEta() {
    snap->setBarostat(eta);
  }

}


