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
 * [4] Kuang & Gezelter,  J. Chem. Phys. 133, 164101 (2010).
 */


#include <algorithm>
#include <fstream>
#include "applications/staticProps/RNEMDStats.hpp"
#include "primitives/Molecule.hpp"
#include "utils/PhysicalConstants.hpp"

namespace OpenMD {
  
  RNEMDZ::RNEMDZ(SimInfo* info, const std::string& filename, 
                 const std::string& sele, int nzbins)
    : SlabStatistics(info, filename, sele, nzbins) {
        
    setOutputName(getPrefix(filename) + ".rnemdZ");
    
    temperature = new OutputData;
    temperature->units =  "K";
    temperature->title =  "Temperature";
    temperature->dataType = odtReal;
    temperature->dataHandling = odhAverage;
    temperature->accumulator.reserve(nBins_);
    for (int i = 0; i < nBins_; i++) 
      temperature->accumulator.push_back( new Accumulator() );
    data_.push_back(temperature);
    
    velocity = new OutputData;
    velocity->units = "angstroms/fs";
    velocity->title =  "Velocity";  
    velocity->dataType = odtVector3;
    velocity->dataHandling = odhAverage;
    velocity->accumulator.reserve(nBins_);
    for (int i = 0; i < nBins_; i++) 
      velocity->accumulator.push_back( new VectorAccumulator() );
    data_.push_back(velocity);
    
    density = new OutputData;
    density->units =  "g cm^-3";
    density->title =  "Density";
    density->dataType = odtReal;
    density->dataHandling = odhAverage;
    density->accumulator.reserve(nBins_);
    for (int i = 0; i < nBins_; i++) 
      density->accumulator.push_back( new Accumulator() );
    data_.push_back(density);
  }

  void RNEMDZ::processFrame(int istep) {
    Molecule* mol;
    RigidBody* rb;
    StuntDouble* sd;
    SimInfo::MoleculeIterator mi;
    Molecule::RigidBodyIterator rbIter;
    int i;

    vector<RealType> binMass(nBins_, 0.0);
    vector<Vector3d> binVel(nBins_, V3Zero);
    vector<RealType> binKE(nBins_, 0.0);
    vector<int> binDof(nBins_, 0);
    vector<int> binCount(nBins_, 0);
    
    
    for (mol = info_->beginMolecule(mi); mol != NULL; 
         mol = info_->nextMolecule(mi)) {
      
      // change the positions of atoms which belong to the rigidbodies
      
      for (rb = mol->beginRigidBody(rbIter); rb != NULL; 
           rb = mol->nextRigidBody(rbIter)) {
        rb->updateAtoms();
      }
    }
    
    if (evaluator_.isDynamic()) {
      seleMan_.setSelectionSet(evaluator_.evaluate());
    }
    
    // loop over the selected atoms:
    
    for (sd = seleMan_.beginSelected(i); sd != NULL; 
         sd = seleMan_.nextSelected(i)) {
      
      // figure out where that object is:
      Vector3d pos = sd->getPos(); 
      currentSnapshot_->wrapVector(pos);

      int bin = getBin(pos);
      binCount[bin]++;

      RealType m = sd->getMass();
      binMass[bin] += m;
      Vector3d vel = sd->getVel();
      binVel[bin] += vel;
      binKE[bin] += 0.5 * (m * vel.lengthSquare());
      binDof[bin] += 3;
      
      if (sd->isDirectional()) {
        Vector3d angMom = sd->getJ();
        Mat3x3d I = sd->getI();
        if (sd->isLinear()) {
          int i = sd->linearAxis();
          int j = (i + 1) % 3;
          int k = (i + 2) % 3;
          binKE[bin] += 0.5 * (angMom[j] * angMom[j] / I(j, j) + 
                               angMom[k] * angMom[k] / I(k, k));
          binDof[bin] += 2;
        } else {
          binKE[bin] += 0.5 * (angMom[0] * angMom[0] / I(0, 0) +
                               angMom[1] * angMom[1] / I(1, 1) +
                               angMom[2] * angMom[2] / I(2, 2));
          binDof[bin] += 3;
        }
      }
    }
    
    for (int i = 0; i < nBins_; i++) {
      if (binDof[i] > 0) {
        RealType temp = 2.0 * binKE[i] / (binDof[i] * PhysicalConstants::kb *
                                          PhysicalConstants::energyConvert);
        RealType den = binMass[i] * nBins_ * PhysicalConstants::densityConvert 
          / volume_;
        Vector3d vel = binVel[i] / RealType(binCount[i]);
        dynamic_cast<Accumulator *>(temperature->accumulator[i])->add(temp);
        dynamic_cast<VectorAccumulator *>(velocity->accumulator[i])->add(vel);
        dynamic_cast<Accumulator *>(density->accumulator[i])->add(den);
        dynamic_cast<Accumulator *>(counts_->accumulator[i])->add(1);
      }
    }
  }
  
  void RNEMDZ::processStuntDouble(StuntDouble* sd, int bin) {
  }

  RNEMDR::RNEMDR(SimInfo* info, const std::string& filename, 
                 const std::string& sele, int nrbins)
    : ShellStatistics(info, filename, sele, nrbins) {
        

    setOutputName(getPrefix(filename) + ".rnemdR");
    
    temperature = new OutputData;
    temperature->units =  "K";
    temperature->title =  "Temperature";
    temperature->dataType = odtReal;
    temperature->dataHandling = odhAverage;
    temperature->accumulator.reserve(nBins_);
    for (int i = 0; i < nBins_; i++) 
      temperature->accumulator.push_back( new Accumulator() );
    data_.push_back(temperature);
    
    angularVelocity = new OutputData;
    angularVelocity->units = "angstroms^2/fs";
    angularVelocity->title =  "Velocity";  
    angularVelocity->dataType = odtVector3;
    angularVelocity->dataHandling = odhAverage;
    angularVelocity->accumulator.reserve(nBins_);
    for (int i = 0; i < nBins_; i++) 
      angularVelocity->accumulator.push_back( new VectorAccumulator() );
    data_.push_back(angularVelocity);
    
    density = new OutputData;
    density->units =  "g cm^-3";
    density->title =  "Density";
    density->dataType = odtReal;
    density->dataHandling = odhAverage;
    density->accumulator.reserve(nBins_);
    for (int i = 0; i < nBins_; i++) 
      density->accumulator.push_back( new Accumulator() );
    data_.push_back(density);
  }

  void RNEMDR::processStuntDouble(StuntDouble* sd, int bin) {
    RealType mass = sd->getMass();
    Vector3d vel = sd->getVel();
    Vector3d rPos = sd->getPos() - coordinateOrigin_;
    Vector3d aVel = cross(rPos, vel);

    RealType KE = 0.5 * (mass * vel.lengthSquare());
    int dof = 3;

    if (sd->isDirectional()) {
      Vector3d angMom = sd->getJ();
      Mat3x3d I = sd->getI();
      if (sd->isLinear()) {
        int i = sd->linearAxis();
        int j = (i + 1) % 3;
        int k = (i + 2) % 3;
        KE += 0.5 * (angMom[j] * angMom[j] / I(j, j) + 
                     angMom[k] * angMom[k] / I(k, k));
        dof += 2;
      } else {
        KE += 0.5 * (angMom[0] * angMom[0] / I(0, 0) +
                     angMom[1] * angMom[1] / I(1, 1) +
                     angMom[2] * angMom[2] / I(2, 2));
        dof += 3;
      }
    }
    
    RealType temp = 2.0 * KE / (dof * PhysicalConstants::kb *
                                PhysicalConstants::energyConvert);

    RealType rinner = (RealType)bin * binWidth_;
    RealType router = (RealType)(bin+1) * binWidth_;
    RealType den = mass * 3.0 * PhysicalConstants::densityConvert
      / (4.0 * M_PI * (pow(router,3) - pow(rinner,3)));  
    
    dynamic_cast<Accumulator *>(temperature->accumulator[bin])->add(temp);
    dynamic_cast<VectorAccumulator *>(angularVelocity->accumulator[bin])->add(aVel);
    dynamic_cast<Accumulator *>(density->accumulator[bin])->add(den);

  }
}

