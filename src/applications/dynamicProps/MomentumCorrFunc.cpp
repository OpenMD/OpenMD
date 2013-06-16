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
 * [3]  Sun, Lin & Gezelter, J. Chem. Phys. 128, 234107 (2008).          
 * [4]  Kuang & Gezelter,  J. Chem. Phys. 133, 164101 (2010).
 * [5]  Vardeman, Stocker & Gezelter, J. Chem. Theory Comput. 7, 834 (2011).
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

#include "applications/dynamicProps/MomentumCorrFunc.hpp"
#include "utils/PhysicalConstants.hpp"
#include "primitives/Molecule.hpp"

namespace OpenMD {

  // We need all of the positions, velocities, etc. so that we can
  // recalculate pressures and actions on the fly:
  MomentumCorrFunc::MomentumCorrFunc(SimInfo* info, const std::string& filename, 
				 const std::string& sele1, 
				 const std::string& sele2,
                                 long long int memSize)
    : FrameTimeCorrFunc(info, filename, sele1, sele2, 
			DataStorage::dslPosition | 
			DataStorage::dslVelocity,
			memSize){

      setCorrFuncType("MomentumCorrFunc");
      setOutputName(getPrefix(dumpFilename_) + ".momcorr");
      histogram_.resize(nTimeBins_); 
      count_.resize(nTimeBins_);
    }

  void MomentumCorrFunc::correlateFrames(int frame1, int frame2) {
    SimInfo::MoleculeIterator mi1;
    SimInfo::MoleculeIterator mi2;
    Molecule* mol1;
    Molecule* mol2;
    Molecule::AtomIterator ai1;
    Molecule::AtomIterator ai2;
    Atom* atom1;
    Atom* atom2;

    std::vector<Vector3d> atomPositions1;
    std::vector<Vector3d> atomPositions2;
    std::vector<Vector3d> atomVelocity1;
    std::vector<Vector3d> atomVelocity2;
    int thisAtom1, thisAtom2;

    Snapshot* snapshot1 = bsMan_->getSnapshot(frame1);
    Snapshot* snapshot2 = bsMan_->getSnapshot(frame2);

    assert(snapshot1 && snapshot2);

    RealType time1 = snapshot1->getTime();
    RealType time2 = snapshot2->getTime();
       
    int timeBin = int ((time2 - time1) /deltaTime_ + 0.5);

    updateFrame(frame1);       
    atomPositions1.clear();
    for (mol1 = info_->beginMolecule(mi1); mol1 != NULL; 
         mol1 = info_->nextMolecule(mi1)) {
      for(atom1 = mol1->beginAtom(ai1); atom1 != NULL; 
          atom1 = mol1->nextAtom(ai1)) {        
        atomPositions1.push_back(atom1->getPos(frame1));
        atomVelocity1.push_back(atom1->getVel(frame1));
      }
    }
    updateFrame(frame2);       
    atomPositions2.clear();
    for (mol2 = info_->beginMolecule(mi2); mol2 != NULL; 
         mol2 = info_->nextMolecule(mi2)) {
      for(atom2 = mol2->beginAtom(ai2); atom2 != NULL; 
          atom2 = mol2->nextAtom(ai2)) {        
        atomPositions2.push_back(atom2->getPos(frame2));
        atomVelocity2.push_back(atom2->getVel(frame2));
      }
    }

    thisAtom1 = 0;

    for (mol1 = info_->beginMolecule(mi1); mol1 != NULL; 
         mol1 = info_->nextMolecule(mi1)) {
      for(atom1 = mol1->beginAtom(ai1); atom1 != NULL; 
          atom1 = mol1->nextAtom(ai1)) {
        
        Vector3d r1 = atomPositions1[thisAtom1];
        Vector3d p1 = atom1->getMass() * atomVelocity1[thisAtom1];

     	thisAtom2 = 0;

        for (mol2 = info_->beginMolecule(mi2); mol2 != NULL; 
             mol2 = info_->nextMolecule(mi2)) {
          for(atom2 = mol2->beginAtom(ai2); atom2 != NULL; 
              atom2 = mol2->nextAtom(ai2)) {
            
            Vector3d r2 = atomPositions2[thisAtom2];
            Vector3d p2 = atom2->getMass()  * atomVelocity1[thisAtom1];

            Vector3d deltaPos = (r2-r1);
            Vector3d dp2( deltaPos.x() * deltaPos.x(),
                          deltaPos.y() * deltaPos.y(),
                          deltaPos.z() * deltaPos.z());
            Vector3d pprod( p1.x() * p2.x(),
                            p1.y() * p2.y(),
                            p1.z() * p2.z());
            
            histogram_[timeBin] += outProduct(dp2, pprod);
                                           
            thisAtom2++;                    
          }
        }
        
        thisAtom1++;
      } 
    }
    
    count_[timeBin]++;
    
  }

  void MomentumCorrFunc::postCorrelate() {
    for (int i =0 ; i < nTimeBins_; ++i) {
      if (count_[i] > 0) {
        histogram_[i] /= count_[i];
      }
    }
  }

  void MomentumCorrFunc::preCorrelate() {
    // Fill the histogram with empty 3x3 matrices:
    std::fill(histogram_.begin(), histogram_.end(), Mat3x3d(0.0));
    // count array set to zero
    std::fill(count_.begin(), count_.end(), 0);
  }
  


  void MomentumCorrFunc::writeCorrelate() {
    std::ofstream ofs(getOutputFileName().c_str());

    if (ofs.is_open()) {

      ofs << "#" << getCorrFuncType() << "\n";
      ofs << "#time\tcorrTensor\txx\txy\txz\tyx\tyy\tyz\tzx\tzy\tzz\n";

      for (int i = 0; i < nTimeBins_; ++i) {
        ofs << time_[i] << "\t" << 
          histogram_[i](0,0) << "\t" <<
          histogram_[i](0,1) << "\t" <<
          histogram_[i](0,2) << "\t" <<
          histogram_[i](1,0) << "\t" <<
          histogram_[i](1,1) << "\t" <<
          histogram_[i](1,2) << "\t" <<
          histogram_[i](2,0) << "\t" <<
          histogram_[i](2,1) << "\t" <<
          histogram_[i](2,2) << "\t" << "\n";
      }
            
    } else {
      sprintf(painCave.errMsg,
              "MomentumCorrFunc::writeCorrelate Error: fail to open %s\n", getOutputFileName().c_str());
      painCave.isFatal = 1;
      simError();        
    }

    ofs.close();    
    
  }

}

