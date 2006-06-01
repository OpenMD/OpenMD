/*
 * Copyright (c) 2005 The University of Notre Dame. All Rights Reserved.
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
 */

#include "applications/staticProps/RippleOP.hpp"
#include "utils/simError.h"
#include "io/DumpReader.hpp"
#include "primitives/Molecule.hpp"
#include "utils/NumericConstant.hpp"
namespace oopse {


  RippleOP::RippleOP(SimInfo* info, const std::string& filename, const std::string& sele1, const std::string& sele2)
    : StaticAnalyser(info, filename),
      selectionScript1_(sele1), selectionScript2_(sele2), evaluator1_(info), evaluator2_(info), 
      seleMan1_(info), seleMan2_(info){
    
    setOutputName(getPrefix(filename) + ".rp2");
    
    evaluator1_.loadScriptString(sele1);
    evaluator2_.loadScriptString(sele2);
    
    if (!evaluator1_.isDynamic()) {
      seleMan1_.setSelectionSet(evaluator1_.evaluate());
    }else {
      sprintf( painCave.errMsg,
	       "--sele1 must be static selection\n");
      painCave.severity = OOPSE_ERROR;
      painCave.isFatal = 1;
      simError();  
    }
    
    if (!evaluator2_.isDynamic()) {
      seleMan2_.setSelectionSet(evaluator2_.evaluate());
    }else {
      sprintf( painCave.errMsg,
	       "--sele2 must be static selection\n");
      painCave.severity = OOPSE_ERROR;
      painCave.isFatal = 1;
      simError();  
    }
    
    if (seleMan1_.getSelectionCount() != seleMan2_.getSelectionCount() ) {
      sprintf( painCave.errMsg,
	       "The number of selected Stuntdoubles are not the same in --sele1 and sele2\n");
      painCave.severity = OOPSE_ERROR;
      painCave.isFatal = 1;
      simError();  

    }
    
    int i;
    int j;
    StuntDouble* sd1;
    StuntDouble* sd2;
    for (sd1 = seleMan1_.beginSelected(i), sd2 = seleMan2_.beginSelected(j);
     sd1 != NULL && sd2 != NULL;
	 sd1 = seleMan1_.nextSelected(i), sd2 = seleMan2_.nextSelected(j)) {
      
      sdPairs_.push_back(std::make_pair(sd1, sd2));
    }
    
  }
  
  void RippleOP::process() {
    Molecule* mol;
    RigidBody* rb;
    SimInfo::MoleculeIterator mi;
    Molecule::RigidBodyIterator rbIter;
  
    DumpReader reader(info_, dumpFilename_);    
    int nFrames = reader.getNFrames();
    
    for (int i = 0; i < nFrames; i += step_) {
      reader.readFrame(i);
      currentSnapshot_ = info_->getSnapshotManager()->getCurrentSnapshot();
      
      
      for (mol = info_->beginMolecule(mi); mol != NULL; mol = info_->nextMolecule(mi)) {
	//change the positions of atoms which belong to the rigidbodies
	for (rb = mol->beginRigidBody(rbIter); rb != NULL; rb = mol->nextRigidBody(rbIter)) {
	  rb->updateAtoms();
	}
      }      
      
      Mat3x3d orderTensorHead(0.0), orderTensorTail(0.0);
      for (std::vector<std::pair<StuntDouble*, StuntDouble*> >::iterator j = sdPairs_.begin(); j != sdPairs_.end(); ++j) {
	Vector3d vecHead = j->first->getElectroFrame().getColumn(2);
	Vector3d vecTail = j->second->getA().getColumn(2);
	orderTensorHead +=outProduct(vecHead, vecHead);
	orderTensorTail +=outProduct(vecTail, vecTail);
      }
    
      orderTensorHead /= sdPairs_.size();
      orderTensorHead -= (RealType)(1.0/3.0) * Mat3x3d::identity();  
      
      orderTensorTail /= sdPairs_.size();
      orderTensorTail -= (RealType)(1.0/3.0) * Mat3x3d::identity();  
      
      Vector3d eigenvaluesHead, eigenvaluesTail;
      Mat3x3d eigenvectorsHead, eigenvectorsTail;    
      Mat3x3d::diagonalize(orderTensorHead, eigenvaluesHead, eigenvectorsHead);
      Mat3x3d::diagonalize(orderTensorTail, eigenvaluesTail, eigenvectorsTail);
      
      int which;
      RealType maxEval = 0.0;
      for(int k = 0; k< 3; k++){
	if(fabs(eigenvaluesHead[k]) > maxEval){
	  which = k;
	  maxEval = fabs(eigenvaluesHead[k]);
	}
      }
      RealType p2Head = 1.5 * maxEval;
      maxEval = 0.0;
      for(int k = 0; k< 3; k++){
	if(fabs(eigenvaluesTail[k]) > maxEval){
	  which = k;
	  maxEval = fabs(eigenvaluesTail[k]);
	}
      }
      RealType p2Tail = 1.5 * maxEval;
      
      //the eigen vector is already normalized in SquareMatrix3::diagonalize
      Vector3d directorHead = eigenvectorsHead.getColumn(which);
      if (directorHead[0] < 0) {
	directorHead.negate();
      }
      Vector3d directorTail = eigenvectorsTail.getColumn(which);
      if (directorTail[0] < 0) {
	directorTail.negate();
      }   

      OrderParam paramHead, paramTail;
      paramHead.p2 = p2Head;
      paramHead.director = directorHead;
      paramTail.p2 = p2Tail;
      paramTail.director = directorTail;
      
      orderParamsHead_.push_back(paramHead);
      orderParamsTail_.push_back(paramTail);       
      
    }

    OrderParam sumOPHead, errsumOPHead;
    OrderParam sumOPTail, errsumOPTail;

    sumOPHead.p2 = 0.0;
    errsumOPHead.p2 = 0.0;
    for (std::size_t i = 0; i < orderParamsHead_.size(); ++i) {
      sumOPHead.p2 += orderParamsHead_[i].p2;
      sumOPHead.director[0] += orderParamsHead_[i].director[0];
      sumOPHead.director[1] += orderParamsHead_[i].director[1];
      sumOPHead.director[2] += orderParamsHead_[i].director[2];
    }

    avgOPHead.p2 = sumOPHead.p2 / (RealType)orderParamsHead_.size();
    avgOPHead.director[0] = sumOPHead.director[0] / (RealType)orderParamsHead_.size();
    avgOPHead.director[1] = sumOPHead.director[1] / (RealType)orderParamsHead_.size();
    avgOPHead.director[2] = sumOPHead.director[2] / (RealType)orderParamsHead_.size();
    for (std::size_t i = 0; i < orderParamsHead_.size(); ++i) {
      errsumOPHead.p2 += pow((orderParamsHead_[i].p2 - avgOPHead.p2), 2);
    }
    errOPHead.p2 = sqrt(errsumOPHead.p2 / (RealType)orderParamsHead_.size());

    sumOPTail.p2 = 0.0;
    errsumOPTail.p2 = 0.0;
    for (std::size_t i = 0; i < orderParamsTail_.size(); ++i) {
      sumOPTail.p2 += orderParamsTail_[i].p2;
      sumOPTail.director[0] += orderParamsTail_[i].director[0];
      sumOPTail.director[1] += orderParamsTail_[i].director[1];
      sumOPTail.director[2] += orderParamsTail_[i].director[2];
    }

    avgOPTail.p2 = sumOPTail.p2 / (RealType)orderParamsTail_.size();
    avgOPTail.director[0] = sumOPTail.director[0] / (RealType)orderParamsTail_.size();
    avgOPTail.director[1] = sumOPTail.director[1] / (RealType)orderParamsTail_.size();
    avgOPTail.director[2] = sumOPTail.director[2] / (RealType)orderParamsTail_.size();
    for (std::size_t i = 0; i < orderParamsTail_.size(); ++i) {
      errsumOPTail.p2 += pow((orderParamsTail_[i].p2 - avgOPTail.p2), 2);
    }
    errOPTail.p2 = sqrt(errsumOPTail.p2 / (RealType)orderParamsTail_.size());
    
    writeP2();
    
  }
  
  void RippleOP::writeP2() {
    
    std::ofstream os(getOutputFileName().c_str());
    os<< "#selection1: (" << selectionScript1_ << ")\n";
    os << "#p2\terror\tdirector_x\tdirector_y\tdiretor_z\n";    
    
    os <<  avgOPHead.p2 << "\t"
       <<  errOPHead.p2 << "\t"
       <<  avgOPHead.director[0] << "\t"
       <<  avgOPHead.director[1] << "\t"
       <<  avgOPHead.director[2] << "\n";
    
    os << "selection2: (" << selectionScript2_ << ")\n";
    os << "#p2\terror\tdirector_x\tdirector_y\tdiretor_z\n";
    
    os <<  avgOPTail.p2 << "\t"
       <<  errOPTail.p2 << "\t"
       <<  avgOPTail.director[0] << "\t"
       <<  avgOPTail.director[1] << "\t"
       <<  avgOPTail.director[2] << "\n";
  }
  
}

