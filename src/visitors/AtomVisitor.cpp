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
 
#include <cstring>
#include "visitors/AtomVisitor.hpp"
#include "primitives/DirectionalAtom.hpp"
#include "primitives/RigidBody.hpp"

namespace oopse {
void BaseAtomVisitor::visit(RigidBody *rb) {
//vector<Atom*> myAtoms;
//vector<Atom*>::iterator atomIter;

//myAtoms = rb->getAtoms();

//for(atomIter = myAtoms.begin(); atomIter != myAtoms.end(); ++atomIter)
//  (*atomIter)->accept(this);
    }

void BaseAtomVisitor::setVisited(Atom *atom) {
    GenericData *data;
    data = atom->getPropertyByName("VISITED");

    //if visited property is not existed, add it as new property
    if (data == NULL) {
        data = new GenericData();
        data->setID("VISITED");
        atom->addProperty(data);
    }
}

bool BaseAtomVisitor::isVisited(Atom *atom) {
    GenericData *data;
    data = atom->getPropertyByName("VISITED");
    return data == NULL ? false : true;
}

bool SSDAtomVisitor::isSSDAtom(const std::string&atomType) {
    std::set<std::string>::iterator strIter;
    strIter = ssdAtomType.find(atomType);
    return strIter != ssdAtomType.end() ? true : false;
}

void SSDAtomVisitor::visit(DirectionalAtom *datom) {
    std::vector<AtomInfo*>atoms;

    //we need to convert SSD into 4 differnet atoms
    //one oxygen atom, two hydrogen atoms and one pseudo atom which is the center of the mass
    //of the water with a dipole moment
    Vector3d h1(0.0, -0.75695, 0.5206);
    Vector3d h2(0.0, 0.75695, 0.5206);
    Vector3d ox(0.0, 0.0, -0.0654);
    Vector3d u(0, 0, 1);
    RotMat3x3d   rotMatrix;
    RotMat3x3d   rotTrans;
    AtomInfo *   atomInfo;
    Vector3d     pos;
    Vector3d     newVec;
    Quat4d       q;
    AtomData *   atomData;
    GenericData *data;
    bool         haveAtomData;

    //if atom is not SSD atom, just skip it
    if (!isSSDAtom(datom->getType()))
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
    q = datom->getQ();
    rotMatrix = datom->getA();

    // We need A^T to convert from body-fixed to space-fixed:
    //transposeMat3(rotMatrix, rotTrans);
    rotTrans = rotMatrix.transpose();

    //center of mass of the water molecule
    //matVecMul3(rotTrans, u, newVec);
    newVec = rotTrans * u;

    atomInfo = new AtomInfo;
    atomInfo->AtomType = "X";
    atomInfo->pos[0] = pos[0];
    atomInfo->pos[1] = pos[1];
    atomInfo->pos[2] = pos[2];
    atomInfo->dipole[0] = newVec[0];
    atomInfo->dipole[1] = newVec[1];
    atomInfo->dipole[2] = newVec[2];

    atomData->addAtomInfo(atomInfo);

    //oxygen
    //matVecMul3(rotTrans, ox, newVec);
    newVec = rotTrans * ox;

    atomInfo = new AtomInfo;
    atomInfo->AtomType = "O";
    atomInfo->pos[0] = pos[0] + newVec[0];
    atomInfo->pos[1] = pos[1] + newVec[1];
    atomInfo->pos[2] = pos[2] + newVec[2];
    atomInfo->dipole[0] = 0.0;
    atomInfo->dipole[1] = 0.0;
    atomInfo->dipole[2] = 0.0;
    atomData->addAtomInfo(atomInfo);

    //hydrogen1
    //matVecMul3(rotTrans, h1, newVec);
    newVec = rotTrans * h1;
    atomInfo = new AtomInfo;
    atomInfo->AtomType = "H";
    atomInfo->pos[0] = pos[0] + newVec[0];
    atomInfo->pos[1] = pos[1] + newVec[1];
    atomInfo->pos[2] = pos[2] + newVec[2];
    atomInfo->dipole[0] = 0.0;
    atomInfo->dipole[1] = 0.0;
    atomInfo->dipole[2] = 0.0;
    atomData->addAtomInfo(atomInfo);

    //hydrogen2
    //matVecMul3(rotTrans, h2, newVec);
    newVec = rotTrans * h2;
    atomInfo = new AtomInfo;
    atomInfo->AtomType = "H";
    atomInfo->pos[0] = pos[0] + newVec[0];
    atomInfo->pos[1] = pos[1] + newVec[1];
    atomInfo->pos[2] = pos[2] + newVec[2];
    atomInfo->dipole[0] = 0.0;
    atomInfo->dipole[1] = 0.0;
    atomInfo->dipole[2] = 0.0;
    atomData->addAtomInfo(atomInfo);

    //add atom data into atom's property

    if (!haveAtomData) {
        atomData->setID("ATOMDATA");
        datom->addProperty(atomData);
    }

    setVisited(datom);
}

const std::string SSDAtomVisitor::toString() {
    char   buffer[65535];
    std::string result;

    sprintf(buffer,
            "------------------------------------------------------------------\n");
    result += buffer;

    sprintf(buffer, "Visitor name: %s\n", visitorName.c_str());
    result += buffer;

    sprintf(buffer,
            "Visitor Description: Convert SSD into 4 different atoms\n");
    result += buffer;

    sprintf(buffer,
            "------------------------------------------------------------------\n");
    result += buffer;

    return result;
}

bool LinearAtomVisitor::isLinearAtom(const std::string& atomType){
    std::set<std::string>::iterator strIter;
    strIter = linearAtomType.find(atomType);

    return strIter != linearAtomType.end() ? true : false;
}

void LinearAtomVisitor::visit(DirectionalAtom* datom){
    std::vector<AtomInfo*> atoms;
    //we need to convert linear into 4 different atoms
    Vector3d c1(0.0, 0.0, -1.8);
    Vector3d c2(0.0, 0.0, -0.6);
    Vector3d c3(0.0, 0.0,  0.6);
    Vector3d c4(0.0, 0.0,  1.8);
    RotMat3x3d rotMatrix;
    RotMat3x3d rotTrans;
    AtomInfo* atomInfo;
    Vector3d pos;
    Vector3d newVec;
    Quat4d q;
    AtomData* atomData;
    GenericData* data;
    bool haveAtomData;

    //if atom is not SSD atom, just skip it
    if(!isLinearAtom(datom->getType()))
        return;

    data = datom->getPropertyByName("ATOMDATA");
    if(data != NULL){
        atomData = dynamic_cast<AtomData*>(data);  
        if(atomData == NULL){
            std::cerr << "can not get Atom Data from " << datom->getType() << std::endl;
            atomData = new AtomData; 
            haveAtomData = false;      
        } else {
            haveAtomData = true;
        }
    } else {
        atomData = new AtomData;
        haveAtomData = false;
    }
   
  
    pos = datom->getPos();
    q = datom->getQ();
    rotMatrix = datom->getA();

    // We need A^T to convert from body-fixed to space-fixed:  
    rotTrans = rotMatrix.transpose();

    newVec = rotTrans * c1;
    atomInfo = new AtomInfo;
    atomInfo->AtomType = "C";
    atomInfo->pos[0] = pos[0] + newVec[0];
    atomInfo->pos[1] = pos[1] + newVec[1];
    atomInfo->pos[2] = pos[2] + newVec[2];
    atomInfo->dipole[0] = 0.0;
    atomInfo->dipole[1] = 0.0;
    atomInfo->dipole[2] = 0.0;
    atomData->addAtomInfo(atomInfo);

    newVec = rotTrans * c2;
    atomInfo = new AtomInfo;
    atomInfo->AtomType = "C";
    atomInfo->pos[0] = pos[0] + newVec[0];
    atomInfo->pos[1] = pos[1] + newVec[1];
    atomInfo->pos[2] = pos[2] + newVec[2];
    atomInfo->dipole[0] = 0.0;
    atomInfo->dipole[1] = 0.0;
    atomInfo->dipole[2] = 0.0;
    atomData->addAtomInfo(atomInfo);

    newVec = rotTrans * c3;
    atomInfo = new AtomInfo;
    atomInfo->AtomType = "C";
    atomInfo->pos[0] = pos[0] + newVec[0];
    atomInfo->pos[1] = pos[1] + newVec[1];
    atomInfo->pos[2] = pos[2] + newVec[2];
    atomInfo->dipole[0] = 0.0;
    atomInfo->dipole[1] = 0.0;
    atomInfo->dipole[2] = 0.0;
    atomData->addAtomInfo(atomInfo);

    newVec = rotTrans * c4;
    atomInfo = new AtomInfo;
    atomInfo->AtomType = "C";
    atomInfo->pos[0] = pos[0] + newVec[0];
    atomInfo->pos[1] = pos[1] + newVec[1];
    atomInfo->pos[2] = pos[2] + newVec[2];
    atomInfo->dipole[0] = 0.0;
    atomInfo->dipole[1] = 0.0;
    atomInfo->dipole[2] = 0.0;
    atomData->addAtomInfo(atomInfo);

    //add atom data into atom's property

    if(!haveAtomData){
        atomData->setID("ATOMDATA");
        datom->addProperty(atomData);
    }

    setVisited(datom);

}

const std::string LinearAtomVisitor::toString(){
  char buffer[65535];
  std::string result;
  
  sprintf(buffer ,"------------------------------------------------------------------\n");
  result += buffer;

  sprintf(buffer ,"Visitor name: %s\n", visitorName.c_str());
  result += buffer;

  sprintf(buffer , "Visitor Description: Convert linear into 4 different atoms\n");
  result += buffer;

  sprintf(buffer ,"------------------------------------------------------------------\n");
  result += buffer;

  return result;
}

//----------------------------------------------------------------------------//

void DefaultAtomVisitor::visit(Atom *atom) {
    AtomData *atomData;
    AtomInfo *atomInfo;
    Vector3d  pos;

    if (isVisited(atom))
        return;

    atomInfo = new AtomInfo;

    atomData = new AtomData;
    atomData->setID("ATOMDATA");

    pos = atom->getPos();
    atomInfo->AtomType = atom->getType();
    atomInfo->pos[0] = pos[0];
    atomInfo->pos[1] = pos[1];
    atomInfo->pos[2] = pos[2];
    atomInfo->dipole[0] = 0.0;
    atomInfo->dipole[1] = 0.0;
    atomInfo->dipole[2] = 0.0;

    atomData->addAtomInfo(atomInfo);

    atom->addProperty(atomData);

    setVisited(atom);
}

void DefaultAtomVisitor::visit(DirectionalAtom *datom) {
    AtomData *atomData;
    AtomInfo *atomInfo;
    Vector3d  pos;
    Vector3d  u;

    if (isVisited(datom))
        return;

    pos = datom->getPos();
    u = datom->getElectroFrame().getColumn(2);

    atomData = new AtomData;
    atomData->setID("ATOMDATA");
    atomInfo = new AtomInfo;

    atomInfo->AtomType = datom->getType();
    atomInfo->pos[0] = pos[0];
    atomInfo->pos[1] = pos[1];
    atomInfo->pos[2] = pos[2];
    atomInfo->dipole[0] = u[0];
    atomInfo->dipole[1] = u[1];
    atomInfo->dipole[2] = u[2];

    atomData->addAtomInfo(atomInfo);

    datom->addProperty(atomData);

    setVisited(datom);
}

const std::string DefaultAtomVisitor::toString() {
    char   buffer[65535];
    std::string result;

    sprintf(buffer,
            "------------------------------------------------------------------\n");
    result += buffer;

    sprintf(buffer, "Visitor name: %s\n", visitorName.c_str());
    result += buffer;

    sprintf(buffer,
            "Visitor Description: copy atom infomation into atom data\n");
    result += buffer;

    sprintf(buffer,
            "------------------------------------------------------------------\n");
    result += buffer;

    return result;
}
} //namespace oopse
