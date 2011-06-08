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
				 const std::string& sele2)
    : FrameTimeCorrFunc(info, filename, sele1, sele2, 
			DataStorage::dslPosition | 
			DataStorage::dslVelocity |
			DataStorage::dslForce |
			DataStorage::dslTorque |			
			DataStorage::dslParticlePot			){

      setCorrFuncType("EnergyCorrFunc");
      setOutputName(getPrefix(dumpFilename_) + ".moment");
      histogram_.resize(nTimeBins_); 
      count_.resize(nTimeBins_);
    }

  void EnergyCorrFunc::correlateFrames(int frame1, int frame2) {
    Snapshot* snapshot1 = bsMan_->getSnapshot(frame1);
    Snapshot* snapshot2 = bsMan_->getSnapshot(frame2);
    assert(snapshot1 && snapshot2);

    RealType time1 = snapshot1->getTime();
    RealType time2 = snapshot2->getTime();
       
    int timeBin = int ((time2 - time1) /deltaTime_ + 0.5);

    Vector3d G_t_frame1 = G_t_[frame1];
    Vector3d G_t_frame2 = G_t_[frame2];
    
    
    RealType diff_x = G_t_frame1.x()-G_t_frame2.x();
    RealType diff_x_sq = diff_x * diff_x;
    
    RealType diff_y = G_t_frame1.y()-G_t_frame2.y();
    RealType diff_y_sq = diff_y * diff_y;
    
    RealType diff_z = G_t_frame1.z()-G_t_frame2.z();
    RealType diff_z_sq = diff_z*diff_z;    
    
    histogram_[timeBin][0] += diff_x_sq;
    histogram_[timeBin][1] += diff_y_sq;
    histogram_[timeBin][2] += diff_z_sq;
    
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

    forceMan->initialize();

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
            Vector3d vel = atom->getVel();
            RealType kinetic = mass * (vel[0]*vel[0] + vel[1]*vel[1] + 
                                       vel[2]*vel[2]) / PhysicalConstants::energyConvert;
            RealType potential =  atom->getParticlePot();
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
    
    int frame = 0;
    
    // Do it again to compute G^(kappa)(t) for x,y,z
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

        // this needs to be updated to the frame value:
        particleEnergies = E_a_[j];
        
        int thisAtom = 0;
        Vector3d G_t;
        
        for (mol = info_->beginMolecule(mi); mol != NULL; 
             mol = info_->nextMolecule(mi)) {
          for(atom = mol->beginAtom(ai); atom != NULL; atom = mol->nextAtom(ai)) {
                        
            Vector3d pos = atom->getPos();

            G_t[0] += pos.x()*(particleEnergies[thisAtom]-AvgE_a_[thisAtom]);
            G_t[1] += pos.y()*(particleEnergies[thisAtom]-AvgE_a_[thisAtom]);
            G_t[2] += pos.z()*(particleEnergies[thisAtom]-AvgE_a_[thisAtom]);
            
            thisAtom++;                    
          }
        }

        G_t_.push_back(G_t);
        //std::cerr <<"Frame: " << j <<"\t" << G_t << std::endl;
      } 
      bsMan_->unloadBlock(i);
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
