/* Copyright (c) 2006, 2009, 2010 The University of Notre Dame. All Rights Reserved.
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
 * [4]  Kuang & Gezelter,  J. Chem. Phys. 133, 164101 (2010).
 * [5]  Vardeman, Stocker & Gezelter, J. Chem. Theory Comput. 7, 834 (2011).
 */

#include "applications/staticProps/NanoLength.hpp"
#include "utils/simError.h"
#include "io/DumpReader.hpp"
#include "primitives/Molecule.hpp"
#include "utils/NumericConstant.hpp"

using namespace OpenMD;

bool pairComparator( const evIndex& l, const evIndex& r) { 
  return l.first < r.first; 
}

NanoLength::NanoLength(SimInfo* info,
                       const std::string& filename,
                       const std::string& sele)
  : StaticAnalyser(info, filename), selectionScript_(sele), evaluator_(info), seleMan_(info) {
  setOutputName(getPrefix(filename) + ".length");
  
  osq.open(getOutputFileName().c_str());

  evaluator_.loadScriptString(sele);
  if (!evaluator_.isDynamic()) {
    seleMan_.setSelectionSet(evaluator_.evaluate());
  }
  frameCounter_ = 0;
}

void NanoLength::process() {
  Molecule* mol;
  RigidBody* rb;
  SimInfo::MoleculeIterator mi;
  Molecule::RigidBodyIterator rbIter;
  Molecule::AtomIterator ai;
  StuntDouble* sd;
  Vector3d vec;
  int i;
  
  DumpReader reader(info_, dumpFilename_);
  int nFrames = reader.getNFrames();
  frameCounter_ = 0;

  theAtoms_.reserve(info_->getNGlobalAtoms());

  for (int istep = 0; istep < nFrames; istep += step_) {
    reader.readFrame(istep);
    frameCounter_++;
    currentSnapshot_ = info_->getSnapshotManager()->getCurrentSnapshot();
    RealType time = currentSnapshot_->getTime();
    
    // Clear pos vector between each frame.
    theAtoms_.clear();
    
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
      theAtoms_.push_back(sd);      
    }
    
    RealType rodLength = getLength(theAtoms_);
    
    osq.precision(7);
    if (osq.is_open()){
      osq << time << "\t" << rodLength << std::endl;      
    }
  }
}
    
RealType NanoLength::getLength(std::vector<StuntDouble*> atoms) {
  Vector3d COM(0.0);
  RealType mass = 0.0;
  RealType mtmp;
  for (std::vector<StuntDouble*>::iterator i = atoms.begin(); 
       i != atoms.end(); ++i) {
    mtmp = (*i)->getMass();
    mass += mtmp;
    COM += (*i)->getPos() * mtmp;
  }
  COM /= mass;
  
  // Moment of Inertia calculation
  Mat3x3d Itmp(0.0);    
  for (std::vector<StuntDouble*>::iterator i = atoms.begin(); 
       i != atoms.end(); ++i) {
    
    Mat3x3d IAtom(0.0);  
    mtmp = (*i)->getMass();
    Vector3d delta = (*i)->getPos() - COM;
    IAtom -= outProduct(delta, delta) * mtmp;
    RealType r2 = delta.lengthSquare();
    IAtom(0, 0) += mtmp * r2;
    IAtom(1, 1) += mtmp * r2;
    IAtom(2, 2) += mtmp * r2;
    Itmp += IAtom;
  }
  
  //diagonalize 
  Vector3d evals;
  Mat3x3d evects;
  Mat3x3d::diagonalize(Itmp, evals, evects);
  
  // we need to re-order the axes so that the smallest moment of
  // inertia (which corresponds to the long axis of the rod) is
  // along the z-axis. We'll just reverse the order of the three
  // axes.  Python has an argsort function, but we had to invent our
  // own:
  
  std::vector<evIndex> evals_prime;
  for (int i = 0; i < 3; i++) 
    evals_prime.push_back(std::make_pair(evals[i], i));
  std::sort(evals_prime.begin(), evals_prime.end(), pairComparator);
  
  RotMat3x3d A;
  Mat3x3d I;
  
  for (int i = 0; i < 3; i++) {
    int index = evals_prime[2-i].second;
    A.setColumn(i, evects.getColumn(index));
    I(i,i) = evals[index];
  }
  
  // now project the delta from the center of mass onto the long
  // axis of the object
  
  Vector3d longAxis = A.getColumn(2);
  RealType axisLength = longAxis.length();
  RealType projmin = 0.0;
  RealType projmax = 0.0;
  
  for (std::vector<StuntDouble*>::iterator i = atoms.begin(); 
       i != atoms.end(); ++i) {
    Vector3d delta = (*i)->getPos() - COM;
    RealType projection = dot(delta, longAxis) / axisLength;
    if (projection > projmax) projmax = projection;
    if (projection < projmin) projmin = projection;      
  }
  
  return projmax - projmin;
}


