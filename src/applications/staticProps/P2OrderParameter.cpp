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

#include "applications/staticProps/P2OrderParameter.hpp"
#include "utils/simError.h"
#include "io/DumpReader.hpp"
#include "primitives/Molecule.hpp"
#include "utils/NumericConstant.hpp"

using namespace std;
namespace OpenMD {

  P2OrderParameter::P2OrderParameter(SimInfo* info, const string& filename, 
                                     const string& sele1)
  : StaticAnalyser(info, filename), doVect_(true), doOffset_(false),
    selectionScript1_(sele1), seleMan1_(info), seleMan2_(info),
    evaluator1_(info), evaluator2_(info) {
    
    setOutputName(getPrefix(filename) + ".p2");
    
    evaluator1_.loadScriptString(sele1);
  }

  P2OrderParameter::P2OrderParameter(SimInfo* info, const string& filename, 
                                     const string& sele1, const string& sele2)
  : StaticAnalyser(info, filename), doVect_(false), doOffset_(false),
    selectionScript1_(sele1), selectionScript2_(sele2), seleMan1_(info), 
    seleMan2_(info), evaluator1_(info), evaluator2_(info) {
    
    setOutputName(getPrefix(filename) + ".p2");
    
    evaluator1_.loadScriptString(sele1);
    evaluator2_.loadScriptString(sele2);    
  }

  P2OrderParameter::P2OrderParameter(SimInfo* info, const string& filename, 
                                     const string& sele1, int seleOffset)
  : StaticAnalyser(info, filename), doVect_(false), doOffset_(true), 
    selectionScript1_(sele1), seleMan1_(info), seleMan2_(info), 
    evaluator1_(info), evaluator2_(info), seleOffset_(seleOffset) {
    
    setOutputName(getPrefix(filename) + ".p2");
    
    evaluator1_.loadScriptString(sele1);
  }

  void P2OrderParameter::process() {
    Molecule* mol;
    RigidBody* rb;
    SimInfo::MoleculeIterator mi;
    Molecule::RigidBodyIterator rbIter;
    StuntDouble* sd1;
    StuntDouble* sd2;
    int ii; 
    int jj;
    int vecCount;
 
    DumpReader reader(info_, dumpFilename_);    
    int nFrames = reader.getNFrames();

    for (int i = 0; i < nFrames; i += step_) {
      reader.readFrame(i);
      currentSnapshot_ = info_->getSnapshotManager()->getCurrentSnapshot();
      
      for (mol = info_->beginMolecule(mi); mol != NULL; 
           mol = info_->nextMolecule(mi)) {
        //change the positions of atoms which belong to the rigidbodies
        for (rb = mol->beginRigidBody(rbIter); rb != NULL; 
             rb = mol->nextRigidBody(rbIter)) {
          rb->updateAtoms();
        }        
      }      

      Mat3x3d orderTensor(0.0);
      vecCount = 0;

      seleMan1_.setSelectionSet(evaluator1_.evaluate());
      
      if (doVect_) {
        
        for (sd1 = seleMan1_.beginSelected(ii); sd1 != NULL; 
             sd1 = seleMan1_.nextSelected(ii)) {
          if (sd1->isDirectional()) {
            Vector3d vec = sd1->getA().transpose()*V3Z;
            
            vec.normalize();
            orderTensor += outProduct(vec, vec);
            vecCount++;
          }
        }
  
        orderTensor /= vecCount;

      } else {

        if (doOffset_) {

          for (sd1 = seleMan1_.beginSelected(ii); sd1 != NULL;
               sd1 = seleMan1_.nextSelected(ii)) {

            // This will require careful rewriting if StaticProps is
            // ever parallelized.  For an example, see
            // Thermo::getTaggedAtomPairDistance

            int sd2Index = sd1->getGlobalIndex() + seleOffset_;
            sd2 = info_->getIOIndexToIntegrableObject(sd2Index);
            
            Vector3d vec = sd1->getPos() - sd2->getPos();
            
            if (usePeriodicBoundaryConditions_)
              currentSnapshot_->wrapVector(vec);
            
            vec.normalize();
            
            orderTensor +=outProduct(vec, vec);
            vecCount++;
          }
          
          orderTensor /= vecCount;
        } else {
          
          seleMan2_.setSelectionSet(evaluator2_.evaluate());
          
          if (seleMan1_.getSelectionCount() != seleMan2_.getSelectionCount() ) {
            sprintf( painCave.errMsg,
                     "In frame %d, the number of selected StuntDoubles are\n"
                     "\tnot the same in --sele1 and sele2\n", i);
            painCave.severity = OPENMD_INFO;
            painCave.isFatal = 0;
            simError();            
          }
          
          for (sd1 = seleMan1_.beginSelected(ii), 
                 sd2 = seleMan2_.beginSelected(jj);
               sd1 != NULL && sd2 != NULL;
               sd1 = seleMan1_.nextSelected(ii), 
                 sd2 = seleMan2_.nextSelected(jj)) {
            
            Vector3d vec = sd1->getPos() - sd2->getPos();
            
            if (usePeriodicBoundaryConditions_)
              currentSnapshot_->wrapVector(vec);

            vec.normalize();
            
            orderTensor +=outProduct(vec, vec);
            vecCount++;
          }
          
          orderTensor /= vecCount;
        }
      }
      
      if (vecCount == 0) {
          sprintf( painCave.errMsg,
                   "In frame %d, the number of selected vectors was zero.\n"
                   "\tThis will not give a meaningful order parameter.", i);
          painCave.severity = OPENMD_ERROR;
          painCave.isFatal = 1;
          simError();        
      }

      orderTensor -= (RealType)(1.0/3.0) * Mat3x3d::identity();  
      
      Vector3d eigenvalues;
      Mat3x3d eigenvectors;    

      Mat3x3d::diagonalize(orderTensor, eigenvalues, eigenvectors);
      
      int which(-1);
      RealType maxEval = 0.0;
      for(int k = 0; k< 3; k++){
        if(fabs(eigenvalues[k]) > maxEval){
          which = k;
          maxEval = fabs(eigenvalues[k]);
        }
      }
      RealType p2 = 1.5 * maxEval;
      
      //the eigen vector is already normalized in SquareMatrix3::diagonalize
      Vector3d director = eigenvectors.getColumn(which);
      if (director[0] < 0) {
        director.negate();
      }   

      RealType angle = 0.0;
      vecCount = 0;
      
      if (doVect_) {
        for (sd1 = seleMan1_.beginSelected(ii); sd1 != NULL; 
             sd1 = seleMan1_.nextSelected(ii)) {
          if (sd1->isDirectional()) {
            Vector3d vec = sd1->getA().transpose()*V3Z;
            vec.normalize();
            angle += acos(dot(vec, director));
            vecCount++;
          }
        }
        angle = angle/(vecCount*NumericConstant::PI)*180.0;
        
      } else {
        if (doOffset_) {

          for (sd1 = seleMan1_.beginSelected(ii); sd1 != NULL;
               sd1 = seleMan1_.nextSelected(ii)) {
            
            // This will require careful rewriting if StaticProps is
            // ever parallelized.  For an example, see
            // Thermo::getTaggedAtomPairDistance
            
            int sd2Index = sd1->getGlobalIndex() + seleOffset_;
            sd2 = info_->getIOIndexToIntegrableObject(sd2Index);
            
            Vector3d vec = sd1->getPos() - sd2->getPos();
            if (usePeriodicBoundaryConditions_)
              currentSnapshot_->wrapVector(vec);
            vec.normalize();          
            angle += acos(dot(vec, director)) ;
            vecCount++;
          }
          angle = angle / (vecCount * NumericConstant::PI) * 180.0;

        } else {

          for (sd1 = seleMan1_.beginSelected(ii), 
                 sd2 = seleMan2_.beginSelected(jj);
               sd1 != NULL && sd2 != NULL;
               sd1 = seleMan1_.nextSelected(ii), 
                 sd2 = seleMan2_.nextSelected(jj)) {
            
            Vector3d vec = sd1->getPos() - sd2->getPos();
            if (usePeriodicBoundaryConditions_)
              currentSnapshot_->wrapVector(vec);
            vec.normalize();          
            angle += acos(dot(vec, director)) ;
            vecCount++;
          }
          angle = angle / (vecCount * NumericConstant::PI) * 180.0;
        }
      }

      OrderParam param;
      param.p2 = p2;
      param.director = director;
      param.angle = angle;

      orderParams_.push_back(param);       
    
    }
    
    writeP2();
    
  }

  void P2OrderParameter::writeP2() {

    ofstream os(getOutputFileName().c_str());
    os << "#radial distribution function\n";
    os<< "#selection1: (" << selectionScript1_ << ")\t";
    if (!doVect_) {
      os << "selection2: (" << selectionScript2_ << ")\n";
    }
    os << "#p2\tdirector_x\tdirector_y\tdiretor_z\tangle(degree)\n";    

    for (size_t i = 0; i < orderParams_.size(); ++i) {
      os <<  orderParams_[i].p2 << "\t"
         <<  orderParams_[i].director[0] << "\t"
         <<  orderParams_[i].director[1] << "\t"
         <<  orderParams_[i].director[2] << "\t"
         <<  orderParams_[i].angle << "\n";

    }

  }

}

