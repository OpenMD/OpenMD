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
 
#include <math.h>
#include <iostream>

#ifdef IS_MPI
#include <mpi.h>
#endif //is_mpi

#include "brains/Thermo.hpp"
#include "primitives/Molecule.hpp"
#include "utils/simError.h"
#include "utils/PhysicalConstants.hpp"

namespace OpenMD {

  RealType Thermo::getKinetic() {
    SimInfo::MoleculeIterator miter;
    std::vector<StuntDouble*>::iterator iiter;
    Molecule* mol;
    StuntDouble* integrableObject;    
    Vector3d vel;
    Vector3d angMom;
    Mat3x3d I;
    int i;
    int j;
    int k;
    RealType mass;
    RealType kinetic = 0.0;
    RealType kinetic_global = 0.0;
    
    for (mol = info_->beginMolecule(miter); mol != NULL; mol = info_->nextMolecule(miter)) {
      for (integrableObject = mol->beginIntegrableObject(iiter); integrableObject != NULL; 
	   integrableObject = mol->nextIntegrableObject(iiter)) {
        
	mass = integrableObject->getMass();
	vel = integrableObject->getVel();
        
	kinetic += mass * (vel[0]*vel[0] + vel[1]*vel[1] + vel[2]*vel[2]);
        
	if (integrableObject->isDirectional()) {
	  angMom = integrableObject->getJ();
	  I = integrableObject->getI();

	  if (integrableObject->isLinear()) {
	    i = integrableObject->linearAxis();
	    j = (i + 1) % 3;
	    k = (i + 2) % 3;
	    kinetic += angMom[j] * angMom[j] / I(j, j) + angMom[k] * angMom[k] / I(k, k);
	  } else {                        
	    kinetic += angMom[0]*angMom[0]/I(0, 0) + angMom[1]*angMom[1]/I(1, 1) 
	      + angMom[2]*angMom[2]/I(2, 2);
	  }
	}
            
      }
    }
    
#ifdef IS_MPI

    MPI_Allreduce(&kinetic, &kinetic_global, 1, MPI_REALTYPE, MPI_SUM,
                  MPI_COMM_WORLD);
    kinetic = kinetic_global;

#endif //is_mpi

    kinetic = kinetic * 0.5 / PhysicalConstants::energyConvert;

    return kinetic;
  }

  RealType Thermo::getPotential() {
    RealType potential = 0.0;
    Snapshot* curSnapshot = info_->getSnapshotManager()->getCurrentSnapshot();
    RealType shortRangePot_local =  curSnapshot->statData[Stats::SHORT_RANGE_POTENTIAL] ;

    // Get total potential for entire system from MPI.

#ifdef IS_MPI

    MPI_Allreduce(&shortRangePot_local, &potential, 1, MPI_REALTYPE, MPI_SUM,
                  MPI_COMM_WORLD);
    potential += curSnapshot->statData[Stats::LONG_RANGE_POTENTIAL];

#else

    potential = shortRangePot_local + curSnapshot->statData[Stats::LONG_RANGE_POTENTIAL];

#endif // is_mpi

    return potential;
  }

  RealType Thermo::getTotalE() {
    RealType total;

    total = this->getKinetic() + this->getPotential();
    return total;
  }

  RealType Thermo::getTemperature() {
    
    RealType temperature = ( 2.0 * this->getKinetic() ) / (info_->getNdf()* PhysicalConstants::kb );
    return temperature;
  }

  RealType Thermo::getVolume() { 
    Snapshot* curSnapshot = info_->getSnapshotManager()->getCurrentSnapshot();
    return curSnapshot->getVolume();
  }

  RealType Thermo::getPressure() {

    // Relies on the calculation of the full molecular pressure tensor


    Mat3x3d tensor;
    RealType pressure;

    tensor = getPressureTensor();

    pressure = PhysicalConstants::pressureConvert * (tensor(0, 0) + tensor(1, 1) + tensor(2, 2)) / 3.0;

    return pressure;
  }

  RealType Thermo::getPressure(int direction) {

    // Relies on the calculation of the full molecular pressure tensor

	  
    Mat3x3d tensor;
    RealType pressure;

    tensor = getPressureTensor();

    pressure = PhysicalConstants::pressureConvert * tensor(direction, direction);

    return pressure;
  }

  Mat3x3d Thermo::getPressureTensor() {
    // returns pressure tensor in units amu*fs^-2*Ang^-1
    // routine derived via viral theorem description in:
    // Paci, E. and Marchi, M. J.Phys.Chem. 1996, 100, 4314-4322
    Mat3x3d pressureTensor;
    Mat3x3d p_local(0.0);
    Mat3x3d p_global(0.0);

    SimInfo::MoleculeIterator i;
    std::vector<StuntDouble*>::iterator j;
    Molecule* mol;
    StuntDouble* integrableObject;    
    for (mol = info_->beginMolecule(i); mol != NULL; mol = info_->nextMolecule(i)) {
      for (integrableObject = mol->beginIntegrableObject(j); integrableObject != NULL; 
	   integrableObject = mol->nextIntegrableObject(j)) {

	RealType mass = integrableObject->getMass();
	Vector3d vcom = integrableObject->getVel();
	p_local += mass * outProduct(vcom, vcom);         
      }
    }
    
#ifdef IS_MPI
    MPI_Allreduce(p_local.getArrayPointer(), p_global.getArrayPointer(), 9, MPI_REALTYPE, MPI_SUM, MPI_COMM_WORLD);
#else
    p_global = p_local;
#endif // is_mpi

    RealType volume = this->getVolume();
    Snapshot* curSnapshot = info_->getSnapshotManager()->getCurrentSnapshot();
    Mat3x3d tau = curSnapshot->statData.getTau();

    pressureTensor =  (p_global + PhysicalConstants::energyConvert* tau)/volume;
    
    return pressureTensor;
  }


  void Thermo::saveStat(){
    Snapshot* currSnapshot = info_->getSnapshotManager()->getCurrentSnapshot();
    Stats& stat = currSnapshot->statData;
    
    stat[Stats::KINETIC_ENERGY] = getKinetic();
    stat[Stats::POTENTIAL_ENERGY] = getPotential();
    stat[Stats::TOTAL_ENERGY] = stat[Stats::KINETIC_ENERGY]  + stat[Stats::POTENTIAL_ENERGY] ;
    stat[Stats::TEMPERATURE] = getTemperature();
    stat[Stats::PRESSURE] = getPressure();
    stat[Stats::VOLUME] = getVolume();      

    Mat3x3d tensor =getPressureTensor();
    stat[Stats::PRESSURE_TENSOR_XX] = tensor(0, 0);      
    stat[Stats::PRESSURE_TENSOR_XY] = tensor(0, 1);      
    stat[Stats::PRESSURE_TENSOR_XZ] = tensor(0, 2);      
    stat[Stats::PRESSURE_TENSOR_YX] = tensor(1, 0);      
    stat[Stats::PRESSURE_TENSOR_YY] = tensor(1, 1);      
    stat[Stats::PRESSURE_TENSOR_YZ] = tensor(1, 2);      
    stat[Stats::PRESSURE_TENSOR_ZX] = tensor(2, 0);      
    stat[Stats::PRESSURE_TENSOR_ZY] = tensor(2, 1);      
    stat[Stats::PRESSURE_TENSOR_ZZ] = tensor(2, 2);      


    Globals* simParams = info_->getSimParams();

    if (simParams->haveTaggedAtomPair() && 
        simParams->havePrintTaggedPairDistance()) {
      if ( simParams->getPrintTaggedPairDistance()) {
        
        std::pair<int, int> tap = simParams->getTaggedAtomPair();
        Vector3d pos1, pos2, rab;

#ifdef IS_MPI        
        std::cerr << "tap = " << tap.first << "  " << tap.second << std::endl;

	int mol1 = info_->getGlobalMolMembership(tap.first);
	int mol2 = info_->getGlobalMolMembership(tap.second);
        std::cerr << "mols = " << mol1 << " " << mol2 << std::endl;

        int proc1 = info_->getMolToProc(mol1);
        int proc2 = info_->getMolToProc(mol2);

        std::cerr << " procs = " << proc1 << " " <<proc2 <<std::endl;

	RealType data[3];
        if (proc1 == worldRank) {
          StuntDouble* sd1 = info_->getIOIndexToIntegrableObject(tap.first);
          std::cerr << " on proc " << proc1 << ", sd1 has global index= " << sd1->getGlobalIndex() << std::endl;
          pos1 = sd1->getPos();
          data[0] = pos1.x();
          data[1] = pos1.y();
          data[2] = pos1.z();          
          MPI_Bcast(data, 3, MPI_REALTYPE, proc1, MPI_COMM_WORLD);
        } else {
          MPI_Bcast(data, 3, MPI_REALTYPE, proc1, MPI_COMM_WORLD);
          pos1 = Vector3d(data);
        }


        if (proc2 == worldRank) {
          StuntDouble* sd2 = info_->getIOIndexToIntegrableObject(tap.second);
          std::cerr << " on proc " << proc2 << ", sd2 has global index= " << sd2->getGlobalIndex() << std::endl;
          pos2 = sd2->getPos();
          data[0] = pos2.x();
          data[1] = pos2.y();
          data[2] = pos2.z();          
          MPI_Bcast(data, 3, MPI_REALTYPE, proc2, MPI_COMM_WORLD);
        } else {
          MPI_Bcast(data, 3, MPI_REALTYPE, proc2, MPI_COMM_WORLD);
          pos2 = Vector3d(data);
        }
#else
        StuntDouble* at1 = info_->getIOIndexToIntegrableObject(tap.first);
        StuntDouble* at2 = info_->getIOIndexToIntegrableObject(tap.second);
        pos1 = at1->getPos();
        pos2 = at2->getPos();
#endif        
        rab = pos2 - pos1;
        currSnapshot->wrapVector(rab);
        stat[Stats::TAGGED_PAIR_DISTANCE] =  rab.length();
      }
    }
      
    /**@todo need refactorying*/
    //Conserved Quantity is set by integrator and time is set by setTime
    
  }



 Vector3d Thermo::getBoxDipole() {
    Snapshot* currSnapshot = info_->getSnapshotManager()->getCurrentSnapshot();
    SimInfo::MoleculeIterator miter;
    std::vector<Atom*>::iterator aiter;
    Molecule* mol;
    Atom* atom;
    RealType charge;
    RealType moment(0.0);
    Vector3d ri(0.0);
    Vector3d dipoleVector(0.0);
    Vector3d nPos(0.0);
    Vector3d pPos(0.0);
    RealType nChg(0.0);
    RealType pChg(0.0);
    int nCount = 0;
    int pCount = 0;

    RealType chargeToC = 1.60217733e-19;
    RealType angstromToM = 1.0e-10;    RealType debyeToCm = 3.33564095198e-30;

    for (mol = info_->beginMolecule(miter); mol != NULL;
         mol = info_->nextMolecule(miter)) {

      for (atom = mol->beginAtom(aiter); atom != NULL;
           atom = mol->nextAtom(aiter)) {

        if (atom->isCharge() ) {
          charge = 0.0;
          GenericData* data = atom->getAtomType()->getPropertyByName("Charge");
          if (data != NULL) {

            charge = (dynamic_cast<DoubleGenericData*>(data))->getData();
            charge *= chargeToC;

            ri = atom->getPos();
            currSnapshot->wrapVector(ri);
            ri *= angstromToM;

            if (charge < 0.0) {
              nPos += ri;
              nChg -= charge;
              nCount++;
            } else if (charge > 0.0) {
              pPos += ri;
              pChg += charge;
              pCount++;
            }
          }
        }

        if (atom->isDipole() ) {
          Vector3d u_i = atom->getElectroFrame().getColumn(2);
          GenericData* data = dynamic_cast<DirectionalAtomType*>(atom->getAtomType())->getPropertyByName("Dipole");
          if (data != NULL) {
            moment = (dynamic_cast<DoubleGenericData*>(data))->getData();

            moment *= debyeToCm;
            dipoleVector += u_i * moment;
          }
        }
      }
    }


#ifdef IS_MPI
    RealType pChg_global, nChg_global;
    int pCount_global, nCount_global;
    Vector3d pPos_global, nPos_global, dipVec_global;

    MPI_Allreduce(&pChg, &pChg_global, 1, MPI_REALTYPE, MPI_SUM,
                  MPI_COMM_WORLD);
    pChg = pChg_global;
    MPI_Allreduce(&nChg, &nChg_global, 1, MPI_REALTYPE, MPI_SUM,
                  MPI_COMM_WORLD);
    nChg = nChg_global;
    MPI_Allreduce(&pCount, &pCount_global, 1, MPI_INTEGER, MPI_SUM,
                  MPI_COMM_WORLD);
    pCount = pCount_global;
    MPI_Allreduce(&nCount, &nCount_global, 1, MPI_INTEGER, MPI_SUM,
                  MPI_COMM_WORLD);
    nCount = nCount_global;
    MPI_Allreduce(pPos.getArrayPointer(), pPos_global.getArrayPointer(), 3,
                  MPI_REALTYPE, MPI_SUM, MPI_COMM_WORLD);
    pPos = pPos_global;
    MPI_Allreduce(nPos.getArrayPointer(), nPos_global.getArrayPointer(), 3,
                  MPI_REALTYPE, MPI_SUM, MPI_COMM_WORLD);
    nPos = nPos_global;
    MPI_Allreduce(dipoleVector.getArrayPointer(),
                  dipVec_global.getArrayPointer(), 3,
                  MPI_REALTYPE, MPI_SUM, MPI_COMM_WORLD);
    dipoleVector = dipVec_global;
#endif //is_mpi

    // first load the accumulated dipole moment (if dipoles were present)
    Vector3d boxDipole = dipoleVector;
    // now include the dipole moment due to charges
    // use the lesser of the positive and negative charge totals
    RealType chg_value = nChg <= pChg ? nChg : pChg;

    // find the average positions
    if (pCount > 0 && nCount > 0 ) {
      pPos /= pCount;
      nPos /= nCount;
    }

    // dipole is from the negative to the positive (physics notation)
    boxDipole += (pPos - nPos) * chg_value;

    return boxDipole;
  }


} //end namespace OpenMD
