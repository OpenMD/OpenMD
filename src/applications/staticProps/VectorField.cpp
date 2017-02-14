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
 * [6]  Kuang & Gezelter, Mol. Phys., 110, 691-701 (2012).
 */
 
#include "applications/staticProps/VectorField.hpp"
#include "utils/simError.h"
#include "io/DumpReader.hpp"
#include "primitives/Molecule.hpp"
#include "utils/NumericConstant.hpp"
#include <vector>
#include <algorithm>
#include <fstream>

using namespace std;
namespace OpenMD {
  VectorField::VectorField(SimInfo* info,  
                           const std::string& filename, 
			   const std::string& sele1,
			   RealType voxelSize) 
    : StaticAnalyser(info, filename, 1), 
      selectionScript1_(sele1),  
      seleMan1_(info), evaluator1_(info), voxelSize_(voxelSize){
    
    evaluator1_.loadScriptString(sele1);
    if (!evaluator1_.isDynamic()) {
      seleMan1_.setSelectionSet(evaluator1_.evaluate());
    }

    Mat3x3d hmat = info->getSnapshotManager()->getCurrentSnapshot()->getHmat();

    nBins_(0) = int(hmat(0,0) / voxelSize);
    nBins_(1) = int(hmat(1,1) / voxelSize);
    nBins_(2) = int(hmat(2,2) / voxelSize);


    Vector3d V3Zero(0.0 , 0.0, 0.0);
    //Build the vector field histogram and count histogram, 
    // fill the vector field with the zero vector and
    // fill the count histogram with int(0).
    hist_.resize(nBins_(0));
    count_.resize(nBins_(0));
    for (int i = 0 ; i < nBins_(0); ++i) {
      hist_[i].resize(nBins_(1));
      count_[i].resize(nBins_(1));
      for(int j = 0; j < nBins_(1); ++j) {
        hist_[i][j].resize(nBins_(2));
        count_[i][j].resize(nBins_(2));
	for(int k = 0; k < nBins_(2); k++){
	  hist_[i][j][k] = V3Zero;
	  count_[i][j][k] = 0.0;
	}
      }
    }

    setOutputName(getPrefix(filename) + ".vectorField");
  }
  
  VectorField::~VectorField() {
 
  }
    
  void VectorField::process() {
    Molecule* mol;
    StuntDouble* sd;
    RigidBody* rb;
    SimInfo::MoleculeIterator mi;
    Molecule::RigidBodyIterator rbIter;
    int isd1;
    bool usePeriodicBoundaryConditions_ = info_->getSimParams()->getUsePeriodicBoundaryConditions();

    
    DumpReader reader(info_, dumpFilename_);    
    int nFrames = reader.getNFrames();
    
    // loop over frames of the .dump file
    for (int istep = 0; istep < nFrames; istep += step_) {
      reader.readFrame(istep);
      currentSnapshot_ = info_->getSnapshotManager()->getCurrentSnapshot();
      
      Mat3x3d hmat = currentSnapshot_->getHmat();
      Vector3d halfBox = Vector3d(hmat(0,0), hmat(1,1), hmat(2,2)) / 2.0;
      Mat3x3d invBox = currentSnapshot_->getInvHmat();
     
      
      if (evaluator1_.isDynamic()) {
        seleMan1_.setSelectionSet(evaluator1_.evaluate());
      }
            
      // update the positions of atoms which belong to the rigidbodies
      for (mol = info_->beginMolecule(mi); mol != NULL;
           mol = info_->nextMolecule(mi)) {
        for (rb = mol->beginRigidBody(rbIter); rb != NULL;
             rb = mol->nextRigidBody(rbIter)) {
          rb->updateAtoms();
	}
      }
	  
      
      //Loop over the selected StuntDoubles:
      for (sd = seleMan1_.beginSelected(isd1); sd != NULL;
           sd = seleMan1_.nextSelected(isd1)) {

	//Get the velocity of the sd
	Vector3d vel = sd->getVel();
	
	
	//Get the position of the sd
	Vector3d pos = sd->getPos();
       
	
	//Wrap the sd back into the box, positions now range from (-boxl/2, boxl/2)
	if (usePeriodicBoundaryConditions_){ 
	  currentSnapshot_->wrapVector(pos); 
	  sd->setPos(pos);
	}
	
	//Convert to a scaled position vector, range (-1/2, 1/2)
	// want range to be (0,1), so add 1/2
	Vector3d sPos = invBox * pos;
	for (int i = 0; i < 3; i++) {
	  sPos(i) = sPos(i) + 0.5;
	  //Treat the case where the sd is found on the edge of the box 
	  if (sPos(i) >= 1.0) sPos(i) = sPos(i) - 1.0;
	}
	
	
	//To determine which bins the sd belongs to,  
	// multiply scaled position by nBins
	int xbin = int( sPos.x() * nBins_(0) );
	int ybin = int( sPos.y() * nBins_(1) ); 
	int zbin = int( sPos.z() * nBins_(2) );
       
	
	//add the velocity of the sd to the current total
	if ( (xbin < int(nBins_(0))) && (xbin >= 0) ){ 
	  if ( (ybin < int(nBins_(1))) && (ybin >= 0) ){
	    if ( (zbin < int(nBins_(2))) && (zbin >= 0) ){
	      hist_[xbin][ybin][zbin] = hist_[xbin][ybin][zbin] + vel;
	      ++count_[xbin][ybin][zbin]; 
	    }
	  }
	}
	    
      }// seleMan_1        	    
    }// dumpFile frames
    writeVectorField();
  }//Process 
  

  void VectorField::writeVectorField() {
    // Need to write the output file as (x \t y \t z \t u \t v \t w \n) format
    // where (x,y,z) is the location of the center of the voxel, and (u,v,w) is the velocity vector
    
    Mat3x3d hmat = info_->getSnapshotManager()->getCurrentSnapshot()->getHmat();
    
    // normalize by total number of elements in each voxel:                                                                    
    for(unsigned int i = 0; i < hist_.size(); ++i) {
      for(unsigned int j = 0; j < hist_[i].size(); ++j) {
        for(unsigned int k = 0; k < hist_[i][j].size(); ++k) {
	  if (count_[i][j][k] > 0.0) {
	      hist_[i][j][k] = hist_[i][j][k] / count_[i][j][k];
	  }
	}
      }
    }
    
    std::ofstream VectorFieldstream(outputFilename_.c_str());
    if (VectorFieldstream.is_open()) {
      VectorFieldstream <<  "# Vector Field output file format (x,y,z) (Vx,Vy,Vz)\n";
      VectorFieldstream <<  "# where (x,y,z) is the location of the center of the voxel and (Vx,Vy,Vz) is the \n";
      VectorFieldstream <<  "# average velocity vector for that voxel. \n";
      

      for (std::size_t k = 0; k < hist_[0][0].size(); ++k) {
	for(std::size_t j = 0; j < hist_[0].size(); ++j) {
	  for(std::size_t i = 0; i < hist_.size(); ++i) {
	    VectorFieldstream << float(hmat(0,0)* (float(i)/nBins_(0) - 0.5)) <<  "\t" << float(hmat(1,1) * (float(j)/nBins_(1) - 0.5)) << "\t" <<  float(hmat(2,2) * (float(k)/nBins_(2) - 0.5)) << "\t" << float(hist_[i][j][k](0)) << "\t" << float(hist_[i][j][k](1)) << "\t" << float(hist_[i][j][k](2)) << "\n";
      
	  }
	}
      }
      
    } //if:VectorFieldstream    

  }//writeVectorField() 

} //openmd
