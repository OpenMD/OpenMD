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

#include <fstream>
#include "applications/staticProps/MultipoleSum.hpp"
#include "primitives/Atom.hpp"
#include "primitives/Molecule.hpp"
#include "types/MultipoleAdapter.hpp"
#include "utils/simError.h"
#include "io/DumpReader.hpp"

namespace OpenMD {

  MultipoleSum::MultipoleSum(SimInfo* info, const std::string& filename, 
                             const std::string& sele1, RealType rmax, 
			     int nrbins)
    : StaticAnalyser(info, filename, nrbins), nRBins_(nrbins), rMax_(rmax),
      selectionScript1_(sele1), seleMan1_(info), evaluator1_(info) {
    
    setOutputName(getPrefix(filename) + ".multipoleSum");
    
    evaluator1_.loadScriptString(sele1);
    if (!evaluator1_.isDynamic()) {
      seleMan1_.setSelectionSet(evaluator1_.evaluate());
    }
    
    aveDlength_.clear();
    aveDlength_.resize(nRBins_, 0.0);
    aveQlength_.clear();
    aveQlength_.resize(nRBins_, 0.0 );
    aveDcount_.clear();
    aveDcount_.resize(nRBins_, 0.0);
    aveQcount_.clear();
    aveQcount_.resize(nRBins_, 0.0 );
    aveDproj_.clear();
    aveDproj_.resize(nRBins_, 0.0);
    deltaR_ = rMax_ / nRBins_;
  }

  void MultipoleSum::process() {
    Molecule* mol;
    SimInfo::MoleculeIterator miter;
    vector<Atom*>::iterator aiter;
    Atom* atom;
    StuntDouble* sd1;
    int i1;
    Vector3d pos1;
    Vector3d ri;
    std::vector<RealType> dipoleHist(nRBins_, 0.0); 
    std::vector<RealType> qpoleHist(nRBins_, 0.0); 
    std::vector<int> lengthCount(nRBins_, 0);
    std::vector<Vector3d> totalDipole; 
    std::vector<Mat3x3d> totalQpole; 
    std::vector<int> dipoleCount; 
    std::vector<int> qpoleCount; 
    std::vector<RealType> dipoleProjection;
    Vector3d dipole;
    Mat3x3d qpole;
    bool usePeriodicBoundaryConditions_ = info_->getSimParams()->getUsePeriodicBoundaryConditions();

    DumpReader reader(info_, dumpFilename_);    
    int nFrames = reader.getNFrames();

    for (int i = 0; i < nFrames; i += step_) {
      reader.readFrame(i);
      currentSnapshot_ = info_->getSnapshotManager()->getCurrentSnapshot();
      
      if  (evaluator1_.isDynamic()) {
        seleMan1_.setSelectionSet(evaluator1_.evaluate());
      }

      for (sd1 = seleMan1_.beginSelected(i1); sd1 != NULL; 
           sd1 = seleMan1_.nextSelected(i1)) {
      
        pos1 = sd1->getPos();

        totalDipole.clear();
        totalDipole.resize(nRBins_, V3Zero); 
        dipoleCount.clear();
        dipoleCount.resize(nRBins_, 0); 
        totalQpole.clear();
        totalQpole.resize(nRBins_, M3Zero); 
        qpoleCount.clear();
        qpoleCount.resize(nRBins_, 0); 
	dipoleProjection.clear();
	dipoleProjection.resize(nRBins_, 0.0);

        for (mol = info_->beginMolecule(miter); mol != NULL; 
             mol = info_->nextMolecule(miter)) {
          
          for (atom = mol->beginAtom(aiter); atom != NULL;
               atom = mol->nextAtom(aiter)) {

            // ri is vector difference between central site and this atom:
            ri = atom->getPos() - pos1;

            if (usePeriodicBoundaryConditions_)
              currentSnapshot_->wrapVector(ri);
            
            dipole = V3Zero;
            qpole = M3Zero;
            AtomType* atype2 = atom->getAtomType();
            MultipoleAdapter ma2 = MultipoleAdapter(atype2);

            if (ma2.isDipole()) 
              dipole = atom->getDipole();
            if (ma2.isQuadrupole()) 
              qpole = atom->getQuadrupole();

            RealType distance = ri.length();
	    std::size_t bin = int(distance / deltaR_);
            // this multipole is contained within the cutoff spheres that are 
            // larger than the bin:
            if (bin < nRBins_) {
              for (std::size_t j = bin; j < nRBins_; j++) {              
                totalDipole[j] += dipole;
                dipoleCount[j]++;
                totalQpole[j] += qpole;
                qpoleCount[j]++;
              }           
            }
          }
        }
	Vector3d myDipole = sd1->getDipole();
	  
        for (std::size_t j = 0; j < nRBins_; j++) {              
	  RealType myProjection = dot(myDipole, totalDipole[j]) / myDipole.length();

          RealType dipoleLength = totalDipole[j].length();
          RealType Qtrace = totalQpole[j].trace();
          RealType Qddot = doubleDot(totalQpole[j], totalQpole[j]);
          RealType qpoleLength =  2.0*(3.0*Qddot - Qtrace*Qtrace);
          dipoleHist[j] += dipoleLength;
          qpoleHist[j] += qpoleLength;
          aveDcount_[j] += dipoleCount[j];
          aveQcount_[j] += qpoleCount[j];
          lengthCount[j] += 1;
	  dipoleProjection[j] += myProjection;
        }
      }
    }
    
    int nSelected = seleMan1_.getSelectionCount();
    for (std::size_t j = 0; j < nRBins_; j++) {
      if (lengthCount[j] > 0) {
        aveDlength_[j] = dipoleHist[j] / RealType(lengthCount[j]);
        aveQlength_[j] = qpoleHist[j] / RealType(lengthCount[j]);
        aveDcount_[j] /= RealType(nSelected) ;
        aveQcount_[j] /= RealType(nSelected) ;
	aveDproj_[j] = dipoleProjection[j] / RealType(lengthCount[j]);
      } else {
        aveDlength_[j] = 0.0;
        aveQlength_[j] = 0.0;
        aveDcount_[j] = 0.0;
        aveQcount_[j] = 0.0;
	aveDproj_[j] = 0.0;
     }
    }
    writeOut();
  }

  void MultipoleSum::writeOut() {

    ofstream os(getOutputFileName().c_str());
    os << "#multipole sum\n";
    os << "#selection1: (" << selectionScript1_ << ")\t";
    os << "#r\taveDlength\taveDdensity\taveDproj\taveQlength\taveQdensity\n";
    
    for (std::size_t i = 0; i < nRBins_; ++i) {
      RealType r = deltaR_ * i;
      os << r << "\t" << aveDlength_[i] << "\t" 
         << aveDlength_[i] / aveDcount_[i] << "\t"
	 << aveDproj_[i] << "\t" 
         << aveQlength_[i] << "\t"
         << aveQlength_[i] / aveQcount_[i] << "\n";
    }
    os.close();
  }
}


