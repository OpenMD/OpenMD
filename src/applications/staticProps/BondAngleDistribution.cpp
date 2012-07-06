/*
 * Copyright (c) 2007 The University of Notre Dame. All Rights Reserved.
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
 * [4] Kuang & Gezelter,  J. Chem. Phys. 133, 164101 (2010).
 * [4] , Stocker & Gezelter, J. Chem. Theory Comput. 7, 834 (2011). *
 *  Created by J. Daniel Gezelter on 07/27/07.
 *  @author  J. Daniel Gezelter
 *  @version $Id$
 *
 */
 
#include "applications/staticProps/BondAngleDistribution.hpp"
#include "utils/simError.h"
#include "io/DumpReader.hpp"
#include "primitives/Molecule.hpp"
#include "utils/NumericConstant.hpp"

namespace OpenMD {

  BondAngleDistribution::BondAngleDistribution(SimInfo* info, 
                                              const std::string& filename, 
                                              const std::string& sele,
                                              double rCut, int nbins) : StaticAnalyser(info, filename), 
                                              selectionScript_(sele), 
			                                        evaluator_(info), seleMan_(info){
    
    setOutputName(getPrefix(filename) + ".bad");
    
    evaluator_.loadScriptString(sele);
    if (!evaluator_.isDynamic()) {
      seleMan_.setSelectionSet(evaluator_.evaluate());
    }
    
    // Set up cutoff radius and order of the Legendre Polynomial:

    rCut_ = rCut;
    nBins_ = nbins;


    // Theta can take values from 0 to 180    
    deltaTheta_ = (180.0) / nBins_;
    histogram_.resize(nBins_);
  }
  
  BondAngleDistribution::~BondAngleDistribution() {
    histogram_.clear();
  }
  
  void BondAngleDistribution::initalizeHistogram() {
    for (int bin = 0; bin < nBins_; bin++) {      
      histogram_[bin] = 0;
    }
  }
  
  void BondAngleDistribution::process() {
    Molecule* mol;
    Atom* atom;
    RigidBody* rb;
    int myIndex;
    SimInfo::MoleculeIterator mi;
    Molecule::RigidBodyIterator rbIter;
    Molecule::AtomIterator ai;
    StuntDouble* sd;
    Vector3d vec;
    std::vector<Vector3d> bondvec;
    std::vector<RealType> bonddist;
    RealType r;    
    int nBonds;    
    int i;
    
    DumpReader reader(info_, dumpFilename_);    
    int nFrames = reader.getNFrames();
    frameCounter_ = 0;
    
    nTotBonds_ = 0;
    
    for (int istep = 0; istep < nFrames; istep += step_) {
      reader.readFrame(istep);
      frameCounter_++;
      currentSnapshot_ = info_->getSnapshotManager()->getCurrentSnapshot();
      
      if (evaluator_.isDynamic()) {
        seleMan_.setSelectionSet(evaluator_.evaluate());
      }

      // update the positions of atoms which belong to the rigidbodies
      
      for (mol = info_->beginMolecule(mi); mol != NULL; 
           mol = info_->nextMolecule(mi)) {
        for (rb = mol->beginRigidBody(rbIter); rb != NULL; 
             rb = mol->nextRigidBody(rbIter)) {
          rb->updateAtoms();
        }        
      }           
            
      // outer loop is over the selected StuntDoubles:

      for (sd = seleMan_.beginSelected(i); sd != NULL; 
           sd = seleMan_.nextSelected(i)) {

        myIndex = sd->getGlobalIndex();
        nBonds = 0;
        bondvec.clear();
        
        // inner loop is over all other atoms in the system:
        
        for (mol = info_->beginMolecule(mi); mol != NULL; 
             mol = info_->nextMolecule(mi)) {
          for (atom = mol->beginAtom(ai); atom != NULL; 
               atom = mol->nextAtom(ai)) {

            if (atom->getGlobalIndex() != myIndex) {

              vec = sd->getPos() - atom->getPos();       

              if (usePeriodicBoundaryConditions_) 
                currentSnapshot_->wrapVector(vec);
              
              // Calculate "bonds" and make a pair list 
              
              r = vec.length();
              
              // Check to see if neighbor is in bond cutoff 
              
              if (r < rCut_) { 
                // Add neighbor to bond list's
                bondvec.push_back(vec);
                nBonds++;
                nTotBonds_++;
              }  
            }
          }
          
          
          for (int i = 0; i < nBonds-1; i++ ){
            Vector3d vec1 = bondvec[i];
            vec1.normalize();
            for(int j = i+1; j < nBonds; j++){
              Vector3d vec2 = bondvec[j];
              
              vec2.normalize();
             
              RealType theta = acos(dot(vec1,vec2))*180.0/NumericConstant::PI;
               
              
              if (theta > 180.0){
                theta = 360.0 - theta;
              }
              int whichBin = theta/deltaTheta_;
              
              histogram_[whichBin] += 2;
            }
          }           
	      }
      }
    }
    
    
    writeBondAngleDistribution();    
  }
  

  void BondAngleDistribution::writeBondAngleDistribution() {
    
    std::ofstream osbad(getOutputFileName().c_str());


    RealType norm = (RealType)nTotBonds_*((RealType)nTotBonds_-1.0)/2.0;
    if (osbad.is_open()) {
      
      // Normalize by number of frames and write it out:
      for (int i = 0; i < nBins_; ++i) {
        RealType Thetaval = i * deltaTheta_;               
        osbad << Thetaval;        
        osbad << "\t" << (RealType)histogram_[i]/norm/frameCounter_;
        
        osbad << "\n";
      }

      osbad.close();

    } else {
      sprintf(painCave.errMsg, "BondAngleDistribution: unable to open %s\n", 
              (getOutputFileName() + "q").c_str());
      painCave.isFatal = 1;
      simError();  
    }
       
  }
    }
