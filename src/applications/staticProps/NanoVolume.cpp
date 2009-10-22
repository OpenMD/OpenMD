/* Copyright (c) 2006 The University of Notre Dame. All Rights Reserved.
 *
 * The University of Notre Dame grants you ("Licensee") a
 * non-exclusive, royalty free, license to use, modify and
 * redistribute this software in source and binary code form, provided
 * that the following conditions are met:
 *
 * 1. Acknowledgement of the program authors must be made in any
 *    publication of scientific results based in part on use of the
 *    program.  An acceptable form of acknowledgement is citation of
 *    the article in which the program was described (Matthew
 *    A. Meineke, Charles F. Vardeman II, Teng Lin, Christopher
 *    J. Fennell and J. Daniel Gezelter, "OOPSE: An Object-Oriented
 *    Parallel Simulation Engine for Molecular Dynamics,"
 *    J. Comput. Chem. 26, pp. 252-271 (2005))
 *
 * 2. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 *
 * 3. Redistributions in binary form must reproduce the above copyright
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
 *
 *  NanoVolume.cpp
 *
 *  Created by Charles F. Vardeman II on 14 Dec 2006.
 *  @author  Charles F. Vardeman II
 *  @version $Id: NanoVolume.cpp,v 1.9 2009-10-22 14:19:26 gezelter Exp $
 *
 */

#include "applications/staticProps/NanoVolume.hpp"
#include "math/ConvexHull.hpp"
#include "utils/simError.h"
#include "io/DumpReader.hpp"
#include "primitives/Molecule.hpp"
#include "utils/NumericConstant.hpp"

using namespace oopse;

NanoVolume::NanoVolume(SimInfo* info,
                       const std::string& filename,
                       const std::string& sele)
  : StaticAnalyser(info, filename), selectionScript_(sele), evaluator_(info), seleMan_(info) {
  setOutputName(getPrefix(filename) + ".avol");
  
  evaluator_.loadScriptString(sele);
  if (!evaluator_.isDynamic()) {
    seleMan_.setSelectionSet(evaluator_.evaluate());
  }
  frameCounter_ = 0;
  totalVolume_ = 0.0;
}

void NanoVolume::process() {
#if defined(HAVE_QHULL)
  Molecule* mol;
  Atom* atom;
  RigidBody* rb;
  int myIndex;
  SimInfo::MoleculeIterator mi;
  Molecule::RigidBodyIterator rbIter;
  Molecule::AtomIterator ai;
  StuntDouble* sd;
  Vector3d vec;
  int i,j;

#ifdef HAVE_QHULL
  ConvexHull* thishull = new ConvexHull();
#endif

  DumpReader reader(info_, dumpFilename_);
  int nFrames = reader.getNFrames();
  frameCounter_ = 0;

  theAtoms_.reserve(info_->getNGlobalAtoms());

  for (int istep = 0; istep < nFrames; istep += step_) {
    reader.readFrame(istep);
    frameCounter_++;
    currentSnapshot_ = info_->getSnapshotManager()->getCurrentSnapshot();
    
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
      myIndex = sd->getGlobalIndex();
      
    }
    
    /*
    for (mol = info_->beginMolecule(mi); mol != NULL; 
                 mol = info_->nextMolecule(mi)) {
              for (atom = mol->beginAtom(ai); atom != NULL; 
                   atom = mol->nextAtom(ai)) {
                     theAtoms_.push_back(atom);
              }
            }
    */
    // Generate convex hull for this frame.
    thishull->computeHull(theAtoms_);
  //  totalVolume_ += hull->getVolume();		
  }
  //RealType avgVolume = totalVolume_/(RealType) frameCounter_;
  //std::cout.precision(7);
  //std::cout  << avgVolume << std::endl;
/*
  std::ofstream osq(getOutputFileName().c_str());
  osq.precision(7);
  if (osq.is_open()){
      osq << avgVolume << std::endl;

  }
  osq.close();
*/
#else
  sprintf(painCave.errMsg, "NanoVolume: qhull support was not compiled in!\n");
  painCave.isFatal = 1;
  simError();  

#endif

}
