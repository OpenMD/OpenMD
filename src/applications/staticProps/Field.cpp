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
 
#include "applications/staticProps/Field.hpp"
#include "utils/simError.h"
#include "io/DumpReader.hpp"
#include "primitives/Molecule.hpp"
#include "utils/NumericConstant.hpp"
#include <vector>
#include <algorithm>
#include <fstream>
#include <math.h>

using namespace std;
namespace OpenMD {
  
  Field::Field(SimInfo* info,const std::string& filename, 
	       const std::string& sele, RealType voxelSize) 
    : StaticAnalyser(info, filename, 1), 
      selectionScript_(sele),  
      seleMan_(info), evaluator_(info), voxelSize_(voxelSize){

    RealType nObjects_ = (info_->getNGlobalAtoms()+info_->getNGlobalRigidBodies());
    std::cout << nObjects_ << "\n";
    
    evaluator_.loadScriptString(sele);
    if (!evaluator_.isDynamic()) {
      seleMan_.setSelectionSet(evaluator_.evaluate());
    }

    Mat3x3d hmat = info->getSnapshotManager()->getCurrentSnapshot()->getHmat();

    nBins_(0) = int(hmat(0,0) / voxelSize);
    nBins_(1) = int(hmat(1,1) / voxelSize);
    nBins_(2) = int(hmat(2,2) / voxelSize);

    //Build the field histogram and dens histogram, 
    // fill the field and the dens histogram with int(0).
    field_.resize(nBins_(0));
    dens_.resize(nBins_(0));
    for (int i = 0 ; i < nBins_(0); ++i) {
      field_[i].resize(nBins_(1));
      dens_[i].resize(nBins_(1));
      for(int j = 0; j < nBins_(1); ++j) {
        field_[i][j].resize(nBins_(2));
        dens_[i][j].resize(nBins_(2));
	for(int k = 0; k < nBins_(2); k++){
	  field_[i][j][k] = 0.0;
	  dens_[i][j][k] = 0.0;
	}
      }
    }

    setOutputName(getPrefix(filename) + ".field");
  }
  
  Field::~Field() {
 
  }
    
  void Field::process() {
    Molecule* mol;
    StuntDouble* sd;
    RigidBody* rb;
    SimInfo::MoleculeIterator mi;
    Molecule::RigidBodyIterator rbIter;
    int isd;
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
      RealType volume = currentSnapshot_->getVolume();
      RealType lenX_ = hmat(0,0);
      RealType lenY_ = hmat(1,1);
      RealType lenZ_ = hmat(2,2);

      RealType x, y, z, dx, dy, dz;
      RealType reff, rcut;
      int di, dj, dk, ibin, jbin, kbin;
      int igrid, jgrid, kgrid;
      Vector3d scaled;

      // d(x,y,z) = the width of each bin
      dx = lenX_ / nBins_(0);
      dy = lenY_ / nBins_(1);
      dz = lenZ_ / nBins_(2);
      
      // reff will be used in the gaussian weighting of the
      // should be based on r_{effective}
      RealType reffective = pow( (volume / nObjects_), 1./3.);
      rcut = 3.0 * reffective;

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
	  
      
      //Loop over the selected StuntDoubles:
      for (sd = seleMan_.beginSelected(isd); sd != NULL;
           sd = seleMan_.nextSelected(isd)) {

	/* Here is where abstraction might need to happen, or templating:
	   the goal is to be able to querry some scalar property here, to be 
	   accumulated into the histogram of voxels.
	   
	   RealType scalarVal = sd->getScalarVal();
	*/
	RealType scalarVal = 0.0;
	
	//Get the position of the sd
	Vector3d pos = sd->getPos();
       
	
	//Wrap the sd back into the box, positions now range from
	// (-boxl/2, boxl/2)
	if (usePeriodicBoundaryConditions_){ 
	  currentSnapshot_->wrapVector(pos); 
	  sd->setPos(pos);
	}
	
	//Convert to a scaled position vector, range (-1/2, 1/2)
	// want range to be (0,1), so add 1/2
	Vector3d scaled = invBox * pos;

	// wrap the vector back into the unit box by subtracting                
	// integer box numbers                                                  
	for (int j = 0; j < 3; j++) {
	  scaled[j] -= roundMe(scaled[j]);
	  scaled[j] += 0.5;
	  // Handle the special case when an object is exactly on               
	  // the boundary (a scaled coordinate of 1.0 is the same as            
	  // scaled coordinate of 0.0)                                          
	  if (scaled[j] >= 1.0) scaled[j] -= 1.0;
	}

	//find ijk-indices of voxel that atom is in,
	// multiply scaled positions by nBins
	ibin = nBins_(0) * scaled.x();
	jbin = nBins_(1) * scaled.y();
	kbin = nBins_(2) * scaled.z();
	
	// di = the magnitude of distance (in x-dimension) that we 
	// should loop through to add the velocity density to.
	di = (int) (rcut / dx);
	dj = (int) (rcut / dy);
	dk = (int) (rcut / dz);
	
	
	for (int i = -di; i <= di; i++) {
	  igrid = ibin + i;
	  while (igrid >= int(nBins_(0))) { igrid -= int(nBins_(0)); }
	  while (igrid < 0) { igrid += int(nBins_(0)); }
	  
	  x = lenX_ * (RealType(i) / RealType(nBins_(0)) );
	  
	  for (int j = -dj; j <= dj; j++) {
	    jgrid = jbin + j;
	    while (jgrid >= int(nBins_(1))) {jgrid -= int(nBins_(1));}
	    while (jgrid < 0) {jgrid += int(nBins_(1));}
	    
	    y = lenY_ * (RealType(j) / RealType(nBins_(1)));
	    
	    for (int k = -dk; k <= dk; k++) {
	      kgrid = kbin + k;
	      while (kgrid >= int(nBins_(2))) {kgrid -= int(nBins_(2));}
	      while (kgrid < 0) {kgrid += int(nBins_(2));}
	      
	      z = lenZ_ * (RealType(k) / RealType(nBins_(2)));
	      
	      RealType dist = sqrt(x*x + y*y + z*z);
	      
	      dens_[igrid][jgrid][kgrid] += getDensity(dist, reff, rcut);
	      field_[igrid][jgrid][kgrid] += dens_[igrid][jgrid][kgrid] * scalarVal;
	    }//k loop
	  }//j loop
	} //i loop
    
    
      }// seleMan_        	    
    }// dumpFile frames
    writeField();
    
  }// void Field::process()

  void Field::writeField() {
    // Need to write the output file as (x \t y \t z \t <scalar>) format
    // where (x,y,z) is the location of the center of the voxel, and <scalar>
    // is the average scalar value of the voxel
    
    Mat3x3d hmat = info_->getSnapshotManager()->getCurrentSnapshot()->getHmat();
    
    // normalize by total number of elements in each voxel:                   
    for(unsigned int i = 0; i < field_.size(); ++i) {
      for(unsigned int j = 0; j < field_[i].size(); ++j) {
        for(unsigned int k = 0; k < field_[i][j].size(); ++k) {
	  if (dens_[i][j][k] > 0.0) {
	      field_[i][j][k] = field_[i][j][k] / dens_[i][j][k];
	  }
	}
      }
    }

    //Write out the raw data for the scalar field
    std::ofstream Fieldstream(outputFilename_.c_str());
    if (Fieldstream.is_open()) {
      Fieldstream <<  "# Vector Field output file format (x,y,z) <scalar> \n";
      Fieldstream <<  "# where (x,y,z) is the location of the center of the voxel and <scalar> is the \n";
      Fieldstream <<  "# average scalar for that voxel. \n";
      
      for (std::size_t k = 0; k < field_[0][0].size(); ++k) {
	for(std::size_t j = 0; j < field_[0].size(); ++j) {
	  for(std::size_t i = 0; i < field_.size(); ++i) {
	    Fieldstream << float(hmat(0,0)* (float(i)/nBins_(0) - 0.5)) <<  "\t" << float(hmat(1,1) * (float(j)/nBins_(1) - 0.5)) << "\t" <<  float(hmat(2,2) * (float(k)/nBins_(2) - 0.5)) << "\t" << float(field_[i][j][k]) << "\n";
	  }
	}
      }
      
    } //if:Fieldstream    

    //Write out the python script which will be interpreted by Mayavi2
    string tab = "     ";
    string pythonFilename = outputFilename_ + ".py";
    std::ofstream pythonScriptStream(pythonFilename.c_str());
    if (pythonScriptStream.is_open()) {
      pythonScriptStream << "#!/opt/local/bin/python\n\n";
      pythonScriptStream << "__author__ = \"Patrick Louden (plouden@nd.edu)\" \n";
      pythonScriptStream << "__copyright__ = \"Copyright (c) 2016 by the University of Notre Dame\" \n";
      pythonScriptStream << "__license__ = \"OpenMD\"\n\n";
      pythonScriptStream << "import numpy as np\n";
      pythonScriptStream << "from mayavi.mlab import * \n\n";
      pythonScriptStream << "def plotField(inputFileName): \n";
      pythonScriptStream << tab + "inputFile = open(inputFileName, 'r') \n";
      pythonScriptStream << tab + "x = np.array([]) \n";
      pythonScriptStream << tab + "y = np.array([]) \n";
      pythonScriptStream << tab + "z = np.array([]) \n";
      pythonScriptStream << tab + "scalar = np.array([]) \n\n";
      pythonScriptStream << tab + "for line in inputFile:\n";
      pythonScriptStream << tab + tab + "if line.split()[0] != \"#\": \n";
      pythonScriptStream << tab + tab + tab + "x = np.append(x, float(line.strip().split()[0])) \n";
      pythonScriptStream << tab + tab + tab + "y = np.append(y, float(line.strip().split()[1])) \n";
      pythonScriptStream << tab + tab + tab + "z = np.append(z, float(line.strip().split()[2])) \n";
      pythonScriptStream << tab + tab + tab + "scalar = np.append(scalar, float(line.strip().split()[3])) \n\n";
      pythonScriptStream << tab + "obj = quiver3d(x, y, z, scalar, line_width=2, scale_factor=3) \n";
      pythonScriptStream << tab + "return obj \n\n";
      pythonScriptStream << "plotField(\"";
      pythonScriptStream << outputFilename_.c_str();
      pythonScriptStream << "\")";
    }// pythonScriptStream

  }//writeField() 

  RealType Field::getDensity(RealType r, RealType sigma, RealType rcut) {
    RealType sigma2 = sigma*sigma;
    RealType dens = exp(-r*r/(sigma2*2.0)) /
      (pow(2.0*NumericConstant::PI*sigma2, 3));
    RealType dcut = exp(-rcut*rcut/(sigma2*2.0)) /
      (pow(2.0*NumericConstant::PI*sigma2, 3));
    if (r < rcut)
      return dens - dcut;
    else
      return 0.0;
  }










  //Now for the VectorField portion
  VectorField::VectorField(SimInfo* info,const std::string& filename, 
			   const std::string& sele, RealType voxelSize) 
    : Field(info, filename, sele, voxelSize), selectionScript_(sele),  
      seleMan_(info), evaluator_(info), voxelSize_(voxelSize){
    
    RealType nObjects_ = (info_->getNGlobalAtoms()+info_->getNGlobalRigidBodies());
    std::cout << nObjects_ << "\n";
    
    evaluator_.loadScriptString(sele);
    if (!evaluator_.isDynamic()) {
      seleMan_.setSelectionSet(evaluator_.evaluate());
    }

    Mat3x3d hmat = info->getSnapshotManager()->getCurrentSnapshot()->getHmat();

    nBins_(0) = int(hmat(0,0) / voxelSize);
    nBins_(1) = int(hmat(1,1) / voxelSize);
    nBins_(2) = int(hmat(2,2) / voxelSize);


    Vector3d V3Zero(0.0 , 0.0, 0.0);
    //Build the vector field histogram and dens histogram, 
    // fill the vector field with the zero vector and
    // fill the dens histogram with int(0).
    vectorField_.resize(nBins_(0));
    dens_.resize(nBins_(0));
    for (int i = 0 ; i < nBins_(0); ++i) {
      vectorField_[i].resize(nBins_(1));
      dens_[i].resize(nBins_(1));
      for(int j = 0; j < nBins_(1); ++j) {
        vectorField_[i][j].resize(nBins_(2));
        dens_[i][j].resize(nBins_(2));
	for(int k = 0; k < nBins_(2); k++){
	  vectorField_[i][j][k] = V3Zero;
	  dens_[i][j][k] = 0.0;
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
    int isd;
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
      RealType volume = currentSnapshot_->getVolume();
      RealType lenX_ = hmat(0,0);
      RealType lenY_ = hmat(1,1);
      RealType lenZ_ = hmat(2,2);

      RealType x, y, z, dx, dy, dz;
      RealType reff, rcut;
      int di, dj, dk, ibin, jbin, kbin;
      int igrid, jgrid, kgrid;
      Vector3d scaled;

      // d(x,y,z) = the width of each bin
      dx = lenX_ / nBins_(0);
      dy = lenY_ / nBins_(1);
      dz = lenZ_ / nBins_(2);

      // reff will be used in the gaussian weighting of the
      // should be based on r_{effective}
      RealType reffective = pow( (volume / nObjects_), 1./3.);
      rcut = 3.0 * reffective;


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
	  
      
      //Loop over the selected StuntDoubles:
      for (sd = seleMan_.beginSelected(isd); sd != NULL;
           sd = seleMan_.nextSelected(isd)) {

	/* The following should be abstracted or templated into any vector
	   quantity which could be querried, velocity, dipole moment, etc.
	   
	   Vector3d vectorVal = sd->getvectorVal();
	*/
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
	Vector3d scaled = invBox * pos;

	// wrap the vector back into the unit box by subtracting                                                                                                                                                       
	// integer box numbers                                                                                                                                                                                         
	for (int j = 0; j < 3; j++) {
	  scaled[j] -= roundMe(scaled[j]);
	  scaled[j] += 0.5;
	  // Handle the special case when an object is exactly on                                                                                                                                                      
	  // the boundary (a scaled coordinate of 1.0 is the same as                                                                                                                                                   
	  // scaled coordinate of 0.0)                                                                                                                                                                                 
	  if (scaled[j] >= 1.0) scaled[j] -= 1.0;
	}

	//find ijk-indices of voxel that atom is in,
	// multiply scaled positions by nBins
	ibin = nBins_(0) * scaled.x();
	jbin = nBins_(1) * scaled.y();
	kbin = nBins_(2) * scaled.z();
	
	// di = the magnitude of distance (in x-dimension) that we 
	// should loop through to add the velocity density to.
	di = (int) (rcut / dx);
	dj = (int) (rcut / dy);
	dk = (int) (rcut / dz);
	
	
	for (int i = -di; i <= di; i++) {
	  igrid = ibin + i;
	  while (igrid >= int(nBins_(0))) { igrid -= int(nBins_(0)); }
	  while (igrid < 0) { igrid += int(nBins_(0)); }
	  
	  x = lenX_ * (RealType(i) / RealType(nBins_(0)) );
	  
	  for (int j = -dj; j <= dj; j++) {
	    jgrid = jbin + j;
	    while (jgrid >= int(nBins_(1))) {jgrid -= int(nBins_(1));}
	    while (jgrid < 0) {jgrid += int(nBins_(1));}
	    
	    y = lenY_ * (RealType(j) / RealType(nBins_(1)));
	    
	    for (int k = -dk; k <= dk; k++) {
	      kgrid = kbin + k;
	      while (kgrid >= int(nBins_(2))) {kgrid -= int(nBins_(2));}
	      while (kgrid < 0) {kgrid += int(nBins_(2));}
	      
	      z = lenZ_ * (RealType(k) / RealType(nBins_(2)));
	      
	      RealType dist = sqrt(x*x + y*y + z*z);
	      
	      dens_[igrid][jgrid][kgrid] += getDensity(dist, reff, rcut);
	      vectorField_[igrid][jgrid][kgrid] += dens_[igrid][jgrid][kgrid] * vel;
	    }//k loop
	  }//j loop
	} //i loop
    
    
      }// seleMan_1        	    
    }// dumpFile frames
    writeVectorField();
    
  }// void VectorField::process()

  void VectorField::writeVectorField() {
    // Need to write the output file as (x \t y \t z \t u \t v \t w \n) format
    // where (x,y,z) is the location of the center of the voxel, and (u,v,w) is the velocity vector
    
    Mat3x3d hmat = info_->getSnapshotManager()->getCurrentSnapshot()->getHmat();
    
    // normalize by total number of elements in each voxel:                                                                    
    for(unsigned int i = 0; i < vectorField_.size(); ++i) {
      for(unsigned int j = 0; j < vectorField_[i].size(); ++j) {
        for(unsigned int k = 0; k < vectorField_[i][j].size(); ++k) {
	  if (dens_[i][j][k] > 0.0) {
	      vectorField_[i][j][k] = vectorField_[i][j][k] / dens_[i][j][k];
	  }
	}
      }
    }
    
    std::ofstream VectorFieldstream(outputFilename_.c_str());
    if (VectorFieldstream.is_open()) {
      VectorFieldstream <<  "# Vector Field output file format (x,y,z) (Vx,Vy,Vz)\n";
      VectorFieldstream <<  "# where (x,y,z) is the location of the center of the voxel and (Vx,Vy,Vz) is the \n";
      VectorFieldstream <<  "# average velocity vector for that voxel. \n";
      
      for (std::size_t k = 0; k < vectorField_[0][0].size(); ++k) {
	for(std::size_t j = 0; j < vectorField_[0].size(); ++j) {
	  for(std::size_t i = 0; i < vectorField_.size(); ++i) {
	    VectorFieldstream << float(hmat(0,0)* (float(i)/nBins_(0) - 0.5)) <<  "\t" << float(hmat(1,1) * (float(j)/nBins_(1) - 0.5)) << "\t" <<  float(hmat(2,2) * (float(k)/nBins_(2) - 0.5)) << "\t" << float(vectorField_[i][j][k](0)) << "\t" << float(vectorField_[i][j][k](1)) << "\t" << float(vectorField_[i][j][k](2)) << "\n";
      
	  }
	}
      }
      
    } //if:VectorFieldstream    

    string tab = "     ";
    string pythonFilename = outputFilename_ + ".py";
    std::ofstream pythonScriptStream(pythonFilename.c_str());
    if (pythonScriptStream.is_open()) {
      pythonScriptStream << "#!/opt/local/bin/python\n\n";
      pythonScriptStream << "__author__ = \"Patrick Louden (plouden@nd.edu)\" \n";
      pythonScriptStream << "__copyright__ = \"Copyright (c) 2016 by the University of Notre Dame\" \n";
      pythonScriptStream << "__license__ = \"OpenMD\"\n\n";
      pythonScriptStream << "import numpy as np\n";
      pythonScriptStream << "from mayavi.mlab import * \n\n";
      pythonScriptStream << "def plotVectorField(inputFileName): \n";
      pythonScriptStream << tab + "inputFile = open(inputFileName, 'r') \n";
      pythonScriptStream << tab + "x = np.array([]) \n";
      pythonScriptStream << tab + "y = np.array([]) \n";
      pythonScriptStream << tab + "z = np.array([]) \n";
      pythonScriptStream << tab + "vx = np.array([]) \n";
      pythonScriptStream << tab + "vy = np.array([]) \n";
      pythonScriptStream << tab + "vz = np.array([]) \n\n";
      pythonScriptStream << tab + "for line in inputFile:\n";
      pythonScriptStream << tab + tab + "if line.split()[0] != \"#\": \n";
      pythonScriptStream << tab + tab + tab + "x = np.append(x, float(line.strip().split()[0])) \n";
      pythonScriptStream << tab + tab + tab + "y = np.append(y, float(line.strip().split()[1])) \n";
      pythonScriptStream << tab + tab + tab + "z = np.append(z, float(line.strip().split()[2])) \n";
      pythonScriptStream << tab + tab + tab + "vx = np.append(vx, float(line.strip().split()[3])) \n";
      pythonScriptStream << tab + tab + tab + "vy = np.append(vy, float(line.strip().split()[4])) \n";
      pythonScriptStream << tab + tab + tab + "vz = np.append(vz, float(line.strip().split()[5])) \n\n";
      pythonScriptStream << tab + "obj = quiver3d(x, y, z, vx, vy, vz, line_width=2, scale_factor=3) \n";
      pythonScriptStream << tab + "return obj \n\n";
      pythonScriptStream << "plotVectorField(\"";
      pythonScriptStream << outputFilename_.c_str();
      pythonScriptStream << "\")";
    }// pythonScriptStream


  }//writeVectorField() 

  RealType VectorField::getVectorDensity(RealType r, RealType sigma, RealType rcut) {
    RealType sigma2 = sigma*sigma;
    RealType dens = exp(-r*r/(sigma2*2.0)) /
      (pow(2.0*NumericConstant::PI*sigma2, 3));
    RealType dcut = exp(-rcut*rcut/(sigma2*2.0)) /
      (pow(2.0*NumericConstant::PI*sigma2, 3));
    if (r < rcut)
      return dens - dcut;
    else
      return 0.0;
  }
  


} //openmd
