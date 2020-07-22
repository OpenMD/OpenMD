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

#include <vector>
#include <algorithm>
#include <fstream>
#include <cmath>

#include "applications/staticProps/Field.hpp"
#include "utils/simError.h"
#include "utils/Revision.hpp"
#include "io/DumpReader.hpp"
#include "primitives/Molecule.hpp"
#include "utils/Constants.hpp"
#include "types/FixedChargeAdapter.hpp"
#include "types/FluctuatingChargeAdapter.hpp"
#include "types/MultipoleAdapter.hpp"
#include "primitives/RigidBody.hpp"


using namespace std;
namespace OpenMD {

  template<class T>
  Field<T>::Field(SimInfo* info,const std::string& filename, 
	       const std::string& sele, RealType voxelSize) 
    : StaticAnalyser(info, filename, 1), 
      selectionScript_(sele),  
      seleMan_(info), evaluator_(info), voxelSize_(voxelSize){
    
    evaluator_.loadScriptString(sele);
    if (!evaluator_.isDynamic()) {
      seleMan_.setSelectionSet(evaluator_.evaluate());
    }

    int nSelected = seleMan_.getSelectionCount();

    Mat3x3d box;
    Mat3x3d invBox;

    usePeriodicBoundaryConditions_ = info_->getSimParams()->getUsePeriodicBoundaryConditions();

    if (!usePeriodicBoundaryConditions_) {
      box = snap_->getBoundingBox();
      invBox = snap_->getInvBoundingBox();
    } else {
      box = snap_->getHmat();
      invBox = snap_->getInvHmat();
    }

    // reff will be used in the gaussian weighting of the
    // should be based on r_{effective}
    snap_ = info_->getSnapshotManager()->getCurrentSnapshot();
    if (nSelected == 0) {
      reffective_ = 2.0 * voxelSize;
    } else {
      reffective_ = pow( (snap_->getVolume() / nSelected), 1./3.);
    }
    rcut_ = 3.0 * reffective_;    

    Vector3d A = box.getColumn(0);
    Vector3d B = box.getColumn(1);
    Vector3d C = box.getColumn(2);

    // Required for triclinic cells
    Vector3d AxB = cross(A, B);
    Vector3d BxC = cross(B, C);
    Vector3d CxA = cross(C, A);

    // unit vectors perpendicular to the faces of the triclinic cell:
    AxB.normalize();
    BxC.normalize();
    CxA.normalize();

    // A set of perpendicular lengths in triclinic cells:
    RealType Wa = abs(dot(A, BxC));
    RealType Wb = abs(dot(B, CxA));
    RealType Wc = abs(dot(C, AxB));
    
    nBins_.x() = int( Wa / voxelSize_ );
    nBins_.y() = int( Wb / voxelSize_ );
    nBins_.z() = int( Wc / voxelSize_ );
    
    //Build the field histogram and dens histogram and fill with zeros.
    field_.resize(nBins_.x());
    dens_.resize(nBins_.x());   
    for (int i = 0 ; i < nBins_.x(); ++i) {
      field_[i].resize(nBins_.y());
      dens_[i].resize(nBins_.y());
      for(int j = 0; j < nBins_.y(); ++j) {
        field_[i][j].resize(nBins_.z());
        dens_[i][j].resize(nBins_.z());

        std::fill(dens_[i][j].begin(), dens_[i][j].end(), 0.0);
        std::fill(field_[i][j].begin(), field_[i][j].end(), T(0));

      }
    }

    setOutputName(getPrefix(filename) + ".field");
  }

  template<class T>
  Field<T>::~Field() {
  }

  template<class T>
  void Field<T>::process() {

    DumpReader reader(info_, dumpFilename_);    
    int nFrames = reader.getNFrames();
    nProcessed_ = nFrames/step_;

    for (int istep = 0; istep < nFrames; istep += step_) {
      reader.readFrame(istep);
      snap_ = info_->getSnapshotManager()->getCurrentSnapshot();
      processFrame(istep);
    }
    postProcess();    
    writeField();
    writeVisualizationScript();
  }

  template<class T>
  void Field<T>::processFrame(int istep) {

    Mat3x3d box;
    Mat3x3d invBox;
    int di, dj, dk;
    int ibin, jbin, kbin;
    int igrid, jgrid, kgrid;
    Vector3d scaled;
    RealType x, y, z;
    RealType weight, Wa, Wb, Wc;

    if (!usePeriodicBoundaryConditions_) {
      box = snap_->getBoundingBox();
      invBox = snap_->getInvBoundingBox();
    } else {
      box = snap_->getHmat();
      invBox = snap_->getInvHmat();
    }    
    Vector3d A = box.getColumn(0);
    Vector3d B = box.getColumn(1);
    Vector3d C = box.getColumn(2);

    // Required for triclinic cells
    Vector3d AxB = cross(A, B);
    Vector3d BxC = cross(B, C);
    Vector3d CxA = cross(C, A);

    // unit vectors perpendicular to the faces of the triclinic cell:
    AxB.normalize();
    BxC.normalize();
    CxA.normalize();

    // A set of perpendicular lengths in triclinic cells:
    Wa = abs(dot(A, BxC));
    Wb = abs(dot(B, CxA));
    Wc = abs(dot(C, AxB));
        
    if (evaluator_.isDynamic()) {
      seleMan_.setSelectionSet(evaluator_.evaluate());
    }
    
    // di = the magnitude of distance (in x-dimension) that we 
    // should loop through to add the density to.
    di = (int) (rcut_ * nBins_.x() / Wa);
    dj = (int) (rcut_ * nBins_.y() / Wb);
    dk = (int) (rcut_ * nBins_.z() / Wc);

    int isd;
    StuntDouble* sd;

    //Loop over the selected StuntDoubles:
    for (sd = seleMan_.beginSelected(isd); sd != NULL;
	 sd = seleMan_.nextSelected(isd)) {
            
      //Get the position of the sd
      Vector3d pos = sd->getPos();

      //Wrap the sd back into the box, positions now range from
      // (-boxl/2, boxl/2)
      if (usePeriodicBoundaryConditions_){ 
	snap_->wrapVector(pos); 
      }

      // scaled positions relative to the box vectors
      scaled = invBox * pos;
      
      // wrap the vector back into the unit box by subtracting integer box 
      // numbers
      for (int j = 0; j < 3; j++) {
        // Convert to a scaled position vector, range (-1/2, 1/2)
        // want range to be (0,1), so add 1/2
        scaled[j] -= roundMe(scaled[j]);
        scaled[j] += 0.5;
        // Handle the special case when an object is exactly on the
        // boundary (a scaled coordinate of 1.0 is the same as
        // scaled coordinate of 0.0)
        if (scaled[j] >= 1.0) scaled[j] -= 1.0;
      }
      
      // find xyz-indices of cell that stuntDouble is in.
      ibin = (int) (nBins_.x() * scaled.x());
      jbin = (int) (nBins_.y() * scaled.y());
      kbin = (int) (nBins_.z() * scaled.z());
      
      for (int i = -di; i <= di; i++) {
	igrid = ibin + i;
	while (igrid >= int(nBins_.x())) { igrid -= int(nBins_.x()); }
	while (igrid < 0) { igrid += int(nBins_.x()); }
	
	x = Wa * (RealType(i) / RealType(nBins_.x()) );
	
	for (int j = -dj; j <= dj; j++) {
	  jgrid = jbin + j;
	  while (jgrid >= int(nBins_.y())) {jgrid -= int(nBins_.y());}
	  while (jgrid < 0) {jgrid += int(nBins_.y());}
	  
	  y = Wb * (RealType(j) / RealType(nBins_.y()));
	  
	  for (int k = -dk; k <= dk; k++) {
	    kgrid = kbin + k;
	    while (kgrid >= int(nBins_.z())) {kgrid -= int(nBins_.z());}
	    while (kgrid < 0) {kgrid += int(nBins_.z());}
	    
	    z = Wc * (RealType(k) / RealType(nBins_.z()));
	    
	    RealType dist = sqrt(x*x + y*y + z*z);

            weight = getDensity(dist, reffective_, rcut_);
	    dens_[igrid][jgrid][kgrid] += weight;
            field_[igrid][jgrid][kgrid] += weight*getValue(sd);
            
	  }
	}
      }
    }
  }  

  template<class T>
  RealType Field<T>::getDensity(RealType r, RealType sigma, RealType rcut) {
    RealType sigma2 = sigma*sigma;
    RealType dens = exp(-r*r/(sigma2*2.0)) /
      (pow(2.0*Constants::PI*sigma2, 3));
    RealType dcut = exp(-rcut*rcut/(sigma2*2.0)) /
      (pow(2.0*Constants::PI*sigma2, 3));
    if (r < rcut)
      return dens - dcut;
    else
      return 0.0;
  }

  template<class T>
  void Field<T>::postProcess() {
    for(unsigned int i = 0; i < field_.size(); ++i) {
      for(unsigned int j = 0; j < field_[i].size(); ++j) {
        for(unsigned int k = 0; k < field_[i][j].size(); ++k) {
	  if (dens_[i][j][k] > 0.0) {
            field_[i][j][k] /= dens_[i][j][k];
	  }
	}
      }
    }
  }

  template<class T>
  void Field<T>::writeField() {    
    Mat3x3d hmat = info_->getSnapshotManager()->getCurrentSnapshot()->getHmat();
    Vector3d scaled(0.0);

    
    std::ofstream fs(outputFilename_.c_str());
    if (fs.is_open()) {
      
      Revision r;
      fs << "# " << getAnalysisType() << "\n";
      fs << "# OpenMD " << r.getFullRevision() << "\n";
      fs << "# " << r.getBuildDate() << "\n";
      fs << "# selection script: \"" << selectionScript_ << "\"\n";
      fs << "#\n";
      fs << "# Vector Field output file format has coordinates of the center\n";
      fs << "# of the voxel followed by the value of field (either vector or\n";
      fs << "# scalar).\n";

      for (std::size_t k = 0; k < field_[0][0].size(); ++k) {
        scaled.z() = RealType(k)/nBins_.z() - 0.5;
	for(std::size_t j = 0; j < field_[0].size(); ++j) {
          scaled.y() = RealType(j)/nBins_.y() - 0.5;
	  for(std::size_t i = 0; i < field_.size(); ++i) {
            scaled.y() = RealType(j)/nBins_.y() - 0.5;

            Vector3d voxelLoc = hmat * scaled;

            fs << voxelLoc.x() << "\t";
            fs << voxelLoc.y() << "\t";
            fs << voxelLoc.z() << "\t";
            
            fs << writeValue( field_[i][j][k] );

            fs << std::endl;            
	  }
	}
      }      
    }
  }

  template<class T>
  void Field<T>::writeVisualizationScript() {
    string t = "     ";
    string pythonFilename = outputFilename_ + ".py";
    std::ofstream pss(pythonFilename.c_str());
    if (pss.is_open()) {
      pss << "#!/opt/local/bin/python\n\n";
      pss << "__author__ = \"Patrick Louden (plouden@nd.edu)\" \n";
      pss << "__copyright__ = \"Copyright (c) 2004-2020 The University of Notre Dame. All Rights Reserved.\" \n";
      pss << "__license__ = \"OpenMD\"\n\n";
      pss << "import numpy as np\n";
      pss << "from mayavi.mlab import * \n\n";
      pss << "def plotField(inputFileName): \n";
      pss << t + "inputFile = open(inputFileName, 'r') \n";
      pss << t + "x = np.array([]) \n";
      pss << t + "y = np.array([]) \n";
      pss << t + "z = np.array([]) \n";
      if (typeid(T) == typeid(int)) {
	pss << t + "scalarVal = np.array([]) \n\n";
      }
      if (typeid(T) == typeid(Vector3d)) {
	pss << t + "vectorVal.x = np.array([]) \n";
	pss << t + "vectorVal.y = np.array([]) \n";
	pss << t + "vectorVal.z = np.array([]) \n\n";
      }
      pss << t + "for line in inputFile:\n";
      pss << t + t + "if line.split()[0] != \"#\": \n";
      pss << t + t + t + "x = np.append(x, float(line.strip().split()[0])) \n";
      pss << t + t + t + "y = np.append(y, float(line.strip().split()[1])) \n";
      pss << t + t + t + "z = np.append(z, float(line.strip().split()[2])) \n";
       if (typeid(T) == typeid(int)) {
	pss << t + t + t + "scalarVal = np.append(scalarVal, float(line.strip().split()[3])) \n\n";
	pss << t + "obj = quiver3d(x, y, z, scalarVal, line_width=2, scale_factor=3) \n"; 
      }
      if (typeid(T) == typeid(Vector3d)) {
	pss << t + t + t + "vectorVal.x = np.append(vectorVal.x, float(line.strip().split()[3])) \n";
	pss << t + t + t + "vectorVal.y = np.append(vectorVal.y, float(line.strip().split()[4])) \n";
	pss << t + t + t + "vextorVal.z = np.append(vectorVal.z, float(line.strip().split()[5])) \n\n";
	pss << t + "obj = quiver3d(x, y, z, vectorVal.x, vectorVal.y, vectorVal.z, line_width=2, scale_factor=3) \n";
      }
      pss << t + "return obj \n\n";
      pss << "plotField(\"";
      pss << outputFilename_.c_str();
      pss << "\")";
    }
  }

  template <>
  std::string Field<RealType>::writeValue(RealType v) {
    std::stringstream str;
    if (std::isinf(v) || std::isnan(v)) {      
      sprintf( painCave.errMsg,
               "Field detected a numerical error.\n");
        painCave.isFatal = 1;
        simError();
    }
    str << v;
    return str.str();      
  }

  template <>
  std::string Field<Vector3d>::writeValue(Vector3d v) {
    std::stringstream str;
    if (std::isinf(v[0]) || std::isnan(v[0]) || 
        std::isinf(v[1]) || std::isnan(v[1]) || 
        std::isinf(v[2]) || std::isnan(v[2]) ) {      
      sprintf( painCave.errMsg,
               "Field detected a numerical error.\n");
        painCave.isFatal = 1;
        simError();
    }
    str << v[0] << "\t" << v[1] << "\t" << v[2];
    return str.str();      
  }
  
  DensityField::DensityField(SimInfo* info, const std::string& filename,
                             const std::string& sele1, RealType voxelSize) :
    Field<RealType>(info, filename, sele1, voxelSize) {
    setOutputName(getPrefix(filename) + ".densityField");    
  }
  
  DensityField::~DensityField() {
  }

  RealType DensityField::getValue(StuntDouble* sd) {
    return sd->getMass();
  }
  
  ChargeField:: ChargeField(SimInfo* info, const std::string& filename,
                            const std::string& sele1, RealType voxelSize) :
    Field<RealType>(info, filename, sele1, voxelSize) {
    setOutputName(getPrefix(filename) + ".chargeField");    
  }
  
  RealType ChargeField::getValue(StuntDouble* sd) {
    RealType charge = 0.0;
    Atom* atom = static_cast<Atom*>(sd);

    AtomType* atomType = atom->getAtomType();

    FixedChargeAdapter fca = FixedChargeAdapter(atomType);
    if ( fca.isFixedCharge() ) {
      charge = fca.getCharge();
    }
    
    FluctuatingChargeAdapter fqa = FluctuatingChargeAdapter(atomType);
    if ( fqa.isFluctuatingCharge() ) {
      charge += atom->getFlucQPos();
    }
    return charge;
  }

  
  VelocityField::VelocityField(SimInfo* info, const std::string& filename,
                               const std::string& sele1, RealType voxelSize) :
    Field<Vector3d>(info, filename, sele1, voxelSize) {
    setOutputName(getPrefix(filename) + ".velocityField");    
  }
  
  Vector3d VelocityField::getValue(StuntDouble* sd) {
    return sd->getVel();
  }

  DipoleField::DipoleField(SimInfo* info, const std::string& filename,
                           const std::string& sele1, RealType voxelSize) :
    Field<Vector3d>(info, filename, sele1, voxelSize) {
    setOutputName(getPrefix(filename) + ".dipoleField");    
  }
  
  Vector3d DipoleField::getValue(StuntDouble* sd) {
    const RealType eAtoDebye = 4.8032045444;
    Vector3d dipoleVector(0.0);
    
    if (sd->isDirectionalAtom()) {
      dipoleVector += static_cast<DirectionalAtom*>(sd)->getDipole();
    }

    if (sd->isRigidBody()) {
      RigidBody* rb = static_cast<RigidBody*>(sd);
      Atom* atom;
      RigidBody::AtomIterator ai;
      for (atom = rb->beginAtom(ai); atom != NULL; atom = rb->nextAtom(ai)) {
        RealType charge(0.0);
        AtomType* atomType = atom->getAtomType();

        FixedChargeAdapter fca = FixedChargeAdapter(atomType);
        if ( fca.isFixedCharge() ) {
          charge = fca.getCharge();
        }
        
        FluctuatingChargeAdapter fqa = FluctuatingChargeAdapter(atomType);
        if ( fqa.isFluctuatingCharge() ) {
          charge += atom->getFlucQPos();
        }
        
        Vector3d pos = atom->getPos();        
        dipoleVector += pos * charge * eAtoDebye;

        if (atom->isDipole()) {
          dipoleVector += atom->getDipole();
        }
      }
    }
    return dipoleVector;
  }
} //openmd
