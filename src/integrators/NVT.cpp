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
 
#include "integrators/NVT.hpp"
#include "primitives/Molecule.hpp"
#include "utils/simError.h"
#include "utils/OOPSEConstant.hpp"

namespace oopse {

NVT::NVT(SimInfo* info) : VelocityVerletIntegrator(info), chiTolerance_ (1e-6) {

    Globals* simParams = info_->getSimParams();

    if (simParams->getUseInitXSstate()) {
        Snapshot* currSnapshot = info_->getSnapshotManager()->getCurrentSnapshot();
        currSnapshot->setChi(0.0);
        currSnapshot->setIntegralOfChiDt(0.0);
    }
    
    if (!simParams->haveTargetTemp()) {
        sprintf(painCave.errMsg, "You can't use the NVT integrator without a targetTemp_!\n");
        painCave.isFatal = 1;
        painCave.severity = OOPSE_ERROR;
        simError();
    } else {
        targetTemp_ = simParams->getTargetTemp();
    }

    // We must set tauThermostat_.

    if (!simParams->haveTauThermostat()) {
        sprintf(painCave.errMsg, "If you use the constant temperature\n"
                                     "\tintegrator, you must set tauThermostat_.\n");

        painCave.severity = OOPSE_ERROR;
        painCave.isFatal = 1;
        simError();
    } else {
        tauThermostat_ = simParams->getTauThermostat();
    }

    update();
}

void NVT::doUpdate() {
    oldVel_.resize(info_->getNIntegrableObjects());
    oldJi_.resize(info_->getNIntegrableObjects());    
}
void NVT::moveA() {
    SimInfo::MoleculeIterator i;
    Molecule::IntegrableObjectIterator  j;
    Molecule* mol;
    StuntDouble* integrableObject;
    Vector3d Tb;
    Vector3d ji;
    double mass;
    Vector3d vel;
    Vector3d pos;
    Vector3d frc;

    double chi = currentSnapshot_->getChi();
    double integralOfChidt = currentSnapshot_->getIntegralOfChiDt();
    
    // We need the temperature at time = t for the chi update below:

    double instTemp = thermo.getTemperature();

    for (mol = info_->beginMolecule(i); mol != NULL; mol = info_->nextMolecule(i)) {
        for (integrableObject = mol->beginIntegrableObject(j); integrableObject != NULL;
               integrableObject = mol->nextIntegrableObject(j)) {

        vel = integrableObject->getVel();
        pos = integrableObject->getPos();
        frc = integrableObject->getFrc();

        mass = integrableObject->getMass();

        // velocity half step  (use chi from previous step here):
        //vel[j] += dt2 * ((frc[j] / mass ) * OOPSEConstant::energyConvert - vel[j]*chi);
        vel += dt2 *OOPSEConstant::energyConvert/mass*frc - dt2*chi*vel;
        
        // position whole step
        //pos[j] += dt * vel[j];
        pos += dt * vel;

        integrableObject->setVel(vel);
        integrableObject->setPos(pos);

        if (integrableObject->isDirectional()) {

            //convert the torque to body frame
            Tb = integrableObject->lab2Body(integrableObject->getTrq());

            // get the angular momentum, and propagate a half step

            ji = integrableObject->getJ();

            //ji[j] += dt2 * (Tb[j] * OOPSEConstant::energyConvert - ji[j]*chi);
            ji += dt2*OOPSEConstant::energyConvert*Tb - dt2*chi *ji;
            rotAlgo->rotate(integrableObject, ji, dt);

            integrableObject->setJ(ji);
        }
    }

    }
    
    rattle->constraintA();

    // Finally, evolve chi a half step (just like a velocity) using
    // temperature at time t, not time t+dt/2

    
    chi += dt2 * (instTemp / targetTemp_ - 1.0) / (tauThermostat_ * tauThermostat_);
    integralOfChidt += chi * dt2;

    currentSnapshot_->setChi(chi);
    currentSnapshot_->setIntegralOfChiDt(integralOfChidt);
}

void NVT::moveB() {
    SimInfo::MoleculeIterator i;
    Molecule::IntegrableObjectIterator  j;
    Molecule* mol;
    StuntDouble* integrableObject;
    
    Vector3d Tb;
    Vector3d ji;    
    Vector3d vel;
    Vector3d frc;
    double mass;
    double instTemp;
    int index;
    // Set things up for the iteration:

    double chi = currentSnapshot_->getChi();
    double oldChi = chi;
    double  prevChi;
    double integralOfChidt = currentSnapshot_->getIntegralOfChiDt();

    index = 0;
    for (mol = info_->beginMolecule(i); mol != NULL; mol = info_->nextMolecule(i)) {
        for (integrableObject = mol->beginIntegrableObject(j); integrableObject != NULL;
               integrableObject = mol->nextIntegrableObject(j)) {
                oldVel_[index] = integrableObject->getVel();
                oldJi_[index] = integrableObject->getJ();                

                ++index;    
        }
           
    }

    // do the iteration:

    for(int k = 0; k < maxIterNum_; k++) {
        index = 0;
        instTemp = thermo.getTemperature();

        // evolve chi another half step using the temperature at t + dt/2

        prevChi = chi;
        chi = oldChi + dt2 * (instTemp / targetTemp_ - 1.0) / (tauThermostat_ * tauThermostat_);

        for (mol = info_->beginMolecule(i); mol != NULL; mol = info_->nextMolecule(i)) {
            for (integrableObject = mol->beginIntegrableObject(j); integrableObject != NULL;
                   integrableObject = mol->nextIntegrableObject(j)) {

                frc = integrableObject->getFrc();
                vel = integrableObject->getVel();

                mass = integrableObject->getMass();

                // velocity half step
                //for(j = 0; j < 3; j++)
                //    vel[j] = oldVel_[3*i+j] + dt2 * ((frc[j] / mass ) * OOPSEConstant::energyConvert - oldVel_[3*i + j]*chi);
                vel = oldVel_[index] + dt2/mass*OOPSEConstant::energyConvert * frc - dt2*chi*oldVel_[index];
            
                integrableObject->setVel(vel);

                if (integrableObject->isDirectional()) {

                    // get and convert the torque to body frame

                    Tb =  integrableObject->lab2Body(integrableObject->getTrq());

                    //for(j = 0; j < 3; j++)
                    //    ji[j] = oldJi_[3*i + j] + dt2 * (Tb[j] * OOPSEConstant::energyConvert - oldJi_[3*i+j]*chi);
                    ji = oldJi_[index] + dt2*OOPSEConstant::energyConvert*Tb - dt2*chi *oldJi_[index];

                    integrableObject->setJ(ji);
                }


                ++index;
            }
        }
    

        rattle->constraintB();

        if (fabs(prevChi - chi) <= chiTolerance_)
            break;

    }

    integralOfChidt += dt2 * chi;

    currentSnapshot_->setChi(chi);
    currentSnapshot_->setIntegralOfChiDt(integralOfChidt);
}


double NVT::calcConservedQuantity() {

    double chi = currentSnapshot_->getChi();
    double integralOfChidt = currentSnapshot_->getIntegralOfChiDt();
    double conservedQuantity;
    double fkBT;
    double Energy;
    double thermostat_kinetic;
    double thermostat_potential;
    
    fkBT = info_->getNdf() *OOPSEConstant::kB *targetTemp_;

    Energy = thermo.getTotalE();

    thermostat_kinetic = fkBT * tauThermostat_ * tauThermostat_ * chi * chi / (2.0 * OOPSEConstant::energyConvert);

    thermostat_potential = fkBT * integralOfChidt / OOPSEConstant::energyConvert;

    conservedQuantity = Energy + thermostat_kinetic + thermostat_potential;

    return conservedQuantity;
}


}//end namespace oopse
