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
 * [3]  Sun, Lin & Gezelter, J. Chem. Phys. 128, 24107 (2008).          
 * [4]  Kuang & Gezelter,  J. Chem. Phys. 133, 164101 (2010).
 * [5]  Vardeman, Stocker & Gezelter, J. Chem. Theory Comput. 7, 834 (2011).
 */
 
#include "visitors/RigidBodyVisitor.hpp"
#include "primitives/RigidBody.hpp"


namespace OpenMD {

  void LipidHeadVisitor::visit(RigidBody* rb){
    Vector3d pos;
    Vector3d u(0, 0, 1);
    Vector3d newVec;
    GenericData* data;
    AtomData* atomData;
    AtomInfo* atomInfo;
    bool haveAtomData;
    RotMat3x3d rotMatrix;

    if(!canVisit(rb->getType()))
      return;

    pos = rb->getPos();
    rotMatrix = rb->getA();
    //matVecMul3(rotMatrix, u, newVec);
    newVec = rotMatrix * u;

    data = rb->getPropertyByName("ATOMDATA");

    if(data != NULL){
      atomData = dynamic_cast<AtomData*>(data);  
      
      if(atomData == NULL){
	std::cerr << "can not get Atom Data from " << rb->getType() << std::endl;
	
	atomData = new AtomData; 
	haveAtomData = false;      
	
      } else
	haveAtomData = true;
      
    } else {
      atomData = new AtomData;
      haveAtomData = false;
      
    }

    atomInfo = new AtomInfo;
    atomInfo->atomTypeName = "X";
    atomInfo->pos[0] = pos[0];
    atomInfo->pos[1] = pos[1];
    atomInfo->pos[2] = pos[2];
    atomInfo->vec[0] = newVec[0];
    atomInfo->vec[1] = newVec[1];
    atomInfo->vec[2] = newVec[2];

    atomData->addAtomInfo(atomInfo);

    if(!haveAtomData){
      atomData->setID("ATOMDATA");
      rb->addProperty(atomData);
    }

  }


  void LipidHeadVisitor::addLipidHeadName(const std::string& name){
    lipidHeadName.insert(name);

  }

  bool LipidHeadVisitor::canVisit(const std::string& name){
    return lipidHeadName.find(name) != lipidHeadName.end() ? true : false;

  }


  const  std::string LipidHeadVisitor::toString(){
    char buffer[65535];
    std::string result;
    std::set<std::string>::iterator i;

    sprintf(buffer ,"------------------------------------------------------------------\n");
    result += buffer;

    sprintf(buffer ,"Visitor name: %s\n", visitorName.c_str());
    result += buffer;

    //print the ignore type list
    sprintf(buffer , "lipidHeadName list contains below types:\n");
    result += buffer;

    for(i = lipidHeadName.begin(); i != lipidHeadName.end(); ++i){
      sprintf(buffer ,"%s\t", i->c_str());
      result += buffer;
    }

    sprintf(buffer ,"\n");
    result += buffer;

    sprintf(buffer ,"------------------------------------------------------------------\n");
    result += buffer;

    return result;

  }


  void RBCOMVisitor::visit(RigidBody* rb){
    AtomData* atomData;
    AtomInfo* atomInfo;
    Vector3d pos;

    pos = rb->getPos();
    atomInfo = new AtomInfo;
    atomInfo->atomTypeName = "X";
    atomInfo->pos[0] = pos[0];
    atomInfo->pos[1] = pos[1];
    atomInfo->pos[2] = pos[2];
    atomInfo->vec[0] = 0;
    atomInfo->vec[1] = 0;
    atomInfo->vec[2] = 0;

    atomData = new AtomData; 
    atomData->setID("ATOMDATA");
    atomData->addAtomInfo(atomInfo);

    rb->addProperty(atomData);

  }


  const  std::string RBCOMVisitor::toString(){
    char buffer[65535];
    std::string result;

    sprintf(buffer ,"------------------------------------------------------------------\n");
    result += buffer;

    sprintf(buffer ,"Visitor name: %s\n", visitorName.c_str());
    result += buffer;

    //print the ignore type list
    sprintf(buffer , "Visitor Description: add a pseudo atom at the center of the mass of the rigidbody\n");
    result += buffer;

    sprintf(buffer ,"------------------------------------------------------------------\n");
    result += buffer;

    return result;

  }

}//namespace OpenMD

