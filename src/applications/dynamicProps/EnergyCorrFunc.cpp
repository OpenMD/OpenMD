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

 /* Uses the Helfand-moment method for calculating thermal
  * conductivity using the relation kappa = (N,V)lim(t)->inf 1/(2*k_B*T^2*V*t) <[G_K(t)-G_K(0)]^2>
  * where G_K is the Helfand moment for thermal conductivity definded as
  * G_K(t) = sum_{a=1}{^N} x_a(E_a-<E_a>) and E_a is defined to be
  *	E_a = p_2^2/(2*m)+1/2 sum_{b.ne.a} u(r_ab) where p is momentum and u is pot energy for the 
  * particle pair a-b. This routine calculates E_a, <E_a> and does the correlation
  * <[G_K(t)-G_K(0)]^2>.
  * See Viscardy et al. JCP 126, 184513 (2007)
  */



#include "applications/dynamicProps/EnergyCorrFunc.hpp"
#include "utils/PhysicalConstants.hpp"
#include "brains/ForceManager.hpp"
#include "brains/Thermo.hpp"

namespace OpenMD {

  // We need all of the positions, velocities, etc. so that we can
  // recalculate pressures and actions on the fly:
  EnergyCorrFunc::EnergyCorrFunc(SimInfo* info, const std::string& filename, 
				 const std::string& sele1, 
				 const std::string& sele2,
                                 long long int memSize)
    : FrameTimeCorrFunc(info, filename, sele1, sele2, 
			DataStorage::dslPosition | 
			DataStorage::dslVelocity |
			DataStorage::dslForce |
			DataStorage::dslTorque |			
			DataStorage::dslParticlePot,
			memSize){

      setCorrFuncType("EnergyCorrFunc");
      setOutputName(getPrefix(dumpFilename_) + ".moment");
      histogram_.resize(nTimeBins_); 
      count_.resize(nTimeBins_);
    }

  void EnergyCorrFunc::correlateFrames(int frame1, int frame2) {
    SimInfo::MoleculeIterator mi1;
    SimInfo::MoleculeIterator mi2;
    Molecule::IntegrableObjectIterator mj1;
    Molecule::IntegrableObjectIterator mj2;
    Molecule* mol1;
    Molecule* mol2;
    Molecule::AtomIterator ai1;
    Molecule::AtomIterator ai2;
    Atom* atom1;
    Atom* atom2;
    std::vector<RealType> particleEnergies1;
    std::vector<RealType> particleEnergies2;
    std::vector<Vector3d> atomPositions1;
    std::vector<Vector3d> atomPositions2;
    int thisAtom1, thisAtom2;

    Snapshot* snapshot1 = bsMan_->getSnapshot(frame1);
    Snapshot* snapshot2 = bsMan_->getSnapshot(frame2);
    assert(snapshot1 && snapshot2);

    RealType time1 = snapshot1->getTime();
    RealType time2 = snapshot2->getTime();
       
    int timeBin = int ((time2 - time1) /deltaTime_ + 0.5);

    // now do the correlation

    particleEnergies1 = E_a_[frame1];
    particleEnergies2 = E_a_[frame2];

    updateFrame(frame1);       
    atomPositions1.clear();
    for (mol1 = info_->beginMolecule(mi1); mol1 != NULL; 
         mol1 = info_->nextMolecule(mi1)) {
      for(atom1 = mol1->beginAtom(ai1); atom1 != NULL; 
          atom1 = mol1->nextAtom(ai1)) {        
        atomPositions1.push_back(atom1->getPos(frame1));
      }
    }
    updateFrame(frame2);       
    atomPositions2.clear();
    for (mol2 = info_->beginMolecule(mi2); mol2 != NULL; 
         mol2 = info_->nextMolecule(mi2)) {
      for(atom2 = mol2->beginAtom(ai2); atom2 != NULL; 
          atom2 = mol2->nextAtom(ai2)) {        
        atomPositions2.push_back(atom2->getPos(frame2));
      }
    }

    thisAtom1 = 0;

    for (mol1 = info_->beginMolecule(mi1); mol1 != NULL; 
         mol1 = info_->nextMolecule(mi1)) {
      for(atom1 = mol1->beginAtom(ai1); atom1 != NULL; 
          atom1 = mol1->nextAtom(ai1)) {
        
        Vector3d r1 = atomPositions1[thisAtom1];
        RealType energy1 = particleEnergies1[thisAtom1] - AvgE_a_[thisAtom1];

        thisAtom2 = 0;

        for (mol2 = info_->beginMolecule(mi2); mol2 != NULL; 
             mol2 = info_->nextMolecule(mi2)) {
          for(atom2 = mol2->beginAtom(ai2); atom2 != NULL; 
              atom2 = mol2->nextAtom(ai2)) {
            
            Vector3d r2 = atomPositions2[thisAtom2];
            RealType energy2 = particleEnergies2[thisAtom2] - AvgE_a_[thisAtom2];

            Vector3d deltaPos = (r2-r1);            
            RealType Eprod = energy2*energy1;
            
            histogram_[timeBin][0] += deltaPos.x()*deltaPos.x() * Eprod;
            histogram_[timeBin][1] += deltaPos.y()*deltaPos.y() * Eprod;
            histogram_[timeBin][2] += deltaPos.z()*deltaPos.z() * Eprod;
                                           
            thisAtom2++;                    
          }
        }
        
        thisAtom1++;
      } 
    }
    
    count_[timeBin]++;
    
  }

  void EnergyCorrFunc::postCorrelate() {
    for (int i =0 ; i < nTimeBins_; ++i) {
      if (count_[i] > 0) {
        histogram_[i] /= count_[i];
      }
    }
  }

  void EnergyCorrFunc::preCorrelate() {
    // Fill the histogram with empty 3x3 matrices:
    std::fill(histogram_.begin(), histogram_.end(), 0.0);
    // count array set to zero
    std::fill(count_.begin(), count_.end(), 0);

    SimInfo::MoleculeIterator mi;
    Molecule::IntegrableObjectIterator mj;
    Molecule* mol;
    Molecule::AtomIterator ai;
    Atom* atom;
    std::vector<RealType > particleEnergies;

    // We'll need the force manager to compute forces for the average pressure
    ForceManager* forceMan = new ForceManager(info_);

    forceMan->init();

    // We'll need thermo to compute the pressures from the virial
    Thermo* thermo =  new Thermo(info_);

    // prepare the averages
    RealType pSum = 0.0;
    RealType vSum = 0.0;
    int nsamp = 0;

    // dump files can be enormous, so read them in block-by-block:
    int nblocks = bsMan_->getNBlocks();
    bool firsttime = true;
    int junkframe = 0;
    for (int i = 0; i < nblocks; ++i) {
      bsMan_->loadBlock(i);
      assert(bsMan_->isBlockActive(i));      
      SnapshotBlock block1 = bsMan_->getSnapshotBlock(i);
      for (int j = block1.first; j < block1.second; ++j) {

	// go snapshot-by-snapshot through this block:
        Snapshot* snap = bsMan_->getSnapshot(j);
        
        // update the positions and velocities of the atoms belonging
        // to rigid bodies:
        
        updateFrame(j);        
        
	// do the forces:
        
        forceMan->calcForces(); 
        
        int index = 0;
        
        for (mol = info_->beginMolecule(mi); mol != NULL; 
             mol = info_->nextMolecule(mi)) {
          for(atom = mol->beginAtom(ai); atom != NULL; atom = mol->nextAtom(ai)) {
            RealType mass = atom->getMass();
            Vector3d vel = atom->getVel(j);
            RealType kinetic = mass * (vel[0]*vel[0] + vel[1]*vel[1] + 
                                       vel[2]*vel[2]) / PhysicalConstants::energyConvert;
            RealType potential =  atom->getParticlePot(j);
            RealType eatom = (kinetic + potential)/2.0;
            particleEnergies.push_back(eatom);
            if(firsttime)
              {
                AvgE_a_.push_back(eatom);
              } else {
              /* We assume the the number of atoms does not change.*/
              AvgE_a_[index] += eatom;
            }
            index++;
          }	 
        }
        firsttime = false;
        E_a_.push_back(particleEnergies);
      }
      
      bsMan_->unloadBlock(i);
    }
    
    int nframes =  bsMan_->getNFrames();
    for (int i = 0; i < AvgE_a_.size(); i++){
      AvgE_a_[i] /= nframes;
    }
    
  }
  


  void EnergyCorrFunc::writeCorrelate() {
    std::ofstream ofs(getOutputFileName().c_str());

    if (ofs.is_open()) {

      ofs << "#" << getCorrFuncType() << "\n";
      ofs << "#time\tK_x\tK_y\tK_z\n";

      for (int i = 0; i < nTimeBins_; ++i) {
        ofs << time_[i] << "\t" << 
	      histogram_[i].x() << "\t" <<
	      histogram_[i].y() << "\t" <<
	      histogram_[i].z() << "\t" << "\n";
      }
            
    } else {
      sprintf(painCave.errMsg,
              "EnergyCorrFunc::writeCorrelate Error: fail to open %s\n", getOutputFileName().c_str());
      painCave.isFatal = 1;
      simError();        
    }

    ofs.close();    
    
  }

}
