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
 
#include <cstring>
#include "visitors/ReplacementVisitor.hpp"
#include "primitives/DirectionalAtom.hpp"
#include "primitives/RigidBody.hpp"

namespace OpenMD { 

  void ReplacementVisitor::addReplacedAtomName(const std::string& repName) {
    myTypes_.insert(repName);
  }

  ReplacementVisitor::~ReplacementVisitor(){
    delete sites_;
    myTypes_.clear();
  }
   
  bool ReplacementVisitor::isReplacedAtom(const std::string &atomType) {
    std::set<std::string>::iterator strIter;
    strIter = myTypes_.find(atomType);
    return strIter != myTypes_.end() ? true : false;
  }

  void ReplacementVisitor::addSite(const std::string &name, 
                                   const Vector3d &refPos) {
    AtomInfo* atomInfo = new AtomInfo();
    atomInfo->atomTypeName = name;
    atomInfo->pos = refPos;
    sites_->addAtomInfo(atomInfo);
  }
  void ReplacementVisitor::addSite(const std::string &name, 
                                   const Vector3d &refPos, 
                                   const Vector3d &refVec) {
    AtomInfo* atomInfo = new AtomInfo();
    atomInfo->atomTypeName = name;
    atomInfo->pos = refPos;
    atomInfo->vec = refVec;
    atomInfo->hasVector = true;
    sites_->addAtomInfo(atomInfo);
  }
  
  void ReplacementVisitor::visit(DirectionalAtom *datom) {
    
    RotMat3x3d   A;
    RotMat3x3d   Atrans;
    Mat3x3d      I;
    Vector3d     pos;
    Vector3d     vel;
    Vector3d     frc;
    Vector3d     trq;
    Vector3d     j;
    Mat3x3d      skewMat;

    Vector3d     newVec;
    AtomInfo *   atomInfo;
    AtomData *   atomData;
    GenericData *data;
    bool         haveAtomData;
    
    //if atom is not one of our recognized atom types, just skip it
    if (!isReplacedAtom(datom->getType())) 
      return;
    
    data = datom->getPropertyByName("ATOMDATA");
    
    if (data != NULL) {
      atomData = dynamic_cast<AtomData *>(data);
      
      if (atomData == NULL) {
        std::cerr << "can not get Atom Data from " << datom->getType() << std::endl;
        atomData = new AtomData;
        haveAtomData = false;
      } else
        haveAtomData = true;
    } else {
      atomData = new AtomData;
      haveAtomData = false;
    }
        
    pos = datom->getPos();
    vel = datom->getVel();

    j   = datom->getJ();
    I   = datom->getI();
    A   = datom->getA();

    skewMat(0, 0) =  0;
    skewMat(0, 1) =  j[2] / I(2, 2);
    skewMat(0, 2) = -j[1] / I(1, 1);    
    skewMat(1, 0) = -j[2] / I(2, 2);
    skewMat(1, 1) =  0;
    skewMat(1, 2) =  j[0] / I(0, 0);    
    skewMat(2, 0) =  j[1] / I(1, 1);
    skewMat(2, 1) = -j[0] / I(0, 0);
    skewMat(2, 2) =  0;
    Mat3x3d mat = (A * skewMat).transpose();
    
    // We need A^T to convert from body-fixed to space-fixed:
    Atrans = A.transpose();
    
    AtomInfo* siteInfo;
    std::vector<AtomInfo*>::iterator iter;
    
    for( siteInfo = sites_->beginAtomInfo(iter); siteInfo; 
         siteInfo = sites_->nextAtomInfo(iter) ) {

      newVec = Atrans * siteInfo->pos;    
      
      atomInfo = new AtomInfo;
      atomInfo->atomTypeName = siteInfo->atomTypeName;
      atomInfo->pos = pos + newVec;
      
      if (siteInfo->hasVector) {
        newVec = Atrans * siteInfo->vec;
        atomInfo->vec = newVec;
      } else {
        atomInfo->vec = V3Zero;
      }

      atomInfo->vel = vel + mat * siteInfo->pos;
      atomInfo->hasVelocity = true;
            
      atomData->addAtomInfo(atomInfo);
    }
    if (!haveAtomData) {
      atomData->setID("ATOMDATA");
      datom->addProperty(atomData);
    }
    
    setVisited(datom);
  }
  
  const std::string ReplacementVisitor::toString() {
    char   buffer[65535];
    std::string result;
    
    sprintf(buffer,
            "------------------------------------------------------------------\n");
    result += buffer;
    
    sprintf(buffer, "Visitor name: %s\n", visitorName.c_str());
    result += buffer;
    
    sprintf(buffer,
            "Visitor Description: replace atom with other sites\n");
    result += buffer;
  
    sprintf(buffer,
            "------------------------------------------------------------------\n");
    result += buffer;
    
    return result;
  }
}
