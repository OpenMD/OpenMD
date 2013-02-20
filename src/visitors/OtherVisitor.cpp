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
#include "selection/SelectionManager.hpp" 
#include "visitors/OtherVisitor.hpp"
#include "primitives/DirectionalAtom.hpp"
#include "primitives/RigidBody.hpp"
#include "primitives/Molecule.hpp"
#include "brains/SimInfo.hpp"
#include "brains/Thermo.hpp"

namespace OpenMD {

  void WrappingVisitor::visit(Atom *atom) {
    internalVisit(atom);
  }

  void WrappingVisitor::visit(DirectionalAtom *datom) {
    internalVisit(datom);
  }

  void WrappingVisitor::visit(RigidBody *rb) {
    internalVisit(rb);
  }

  void WrappingVisitor::internalVisit(StuntDouble *sd) {
    GenericData *                     data;
    AtomData *                        atomData;
    AtomInfo *                        atomInfo;
    std::vector<AtomInfo *>::iterator i;

    data = sd->getPropertyByName("ATOMDATA");

    if (data != NULL) {
      atomData = dynamic_cast<AtomData *>(data);

      if (atomData == NULL)
	return;
    } else
      return;

    Snapshot* currSnapshot = info->getSnapshotManager()->getCurrentSnapshot();
    
    for( atomInfo = atomData->beginAtomInfo(i); atomInfo; 
         atomInfo = atomData->nextAtomInfo(i) ) {

      Vector3d newPos = atomInfo->pos - origin_;
      currSnapshot->wrapVector(newPos);
      atomInfo->pos = newPos;

    }
  }

  void WrappingVisitor::update() {
    if (useCom_){
      Thermo thermo(info);
      origin_ = thermo.getCom();
    }
  }
  
  const std::string WrappingVisitor::toString() {
    char        buffer[65535];
    std::string result;

    sprintf(buffer,
            "------------------------------------------------------------------\n");
    result += buffer;

    sprintf(buffer, "Visitor name: %s\n", visitorName.c_str());
    result += buffer;

    sprintf(buffer,
            "Visitor Description: wrapping atoms back to periodic box\n");
    result += buffer;

    sprintf(buffer,
            "------------------------------------------------------------------\n");
    result += buffer;

    return result;
  }

  //----------------------------------------------------------------------------//

  ReplicateVisitor::ReplicateVisitor(SimInfo *info, Vector3i opt) :
    BaseVisitor() {
      this->info = info;
      visitorName = "ReplicateVisitor";
      this->replicateOpt = opt;

      //generate the replicate directions
      for( int i = 0; i <= replicateOpt[0]; i++ ) {
        for( int j = 0; j <= replicateOpt[1]; j++ ) {
	  for( int k = 0; k <= replicateOpt[2]; k++ ) {
	    //skip original frame
	    if (i == 0 && j == 0 && k == 0) {
	      continue; 
	    } else {
	      dir.push_back(Vector3i(i, j, k));
	    }
	  }
        }
      }
    
    }

  void ReplicateVisitor::visit(Atom *atom) {
    internalVisit(atom);
  }

  void ReplicateVisitor::visit(DirectionalAtom *datom) {
    internalVisit(datom);
  }

  void ReplicateVisitor::visit(RigidBody *rb) {
    internalVisit(rb);
  }

  void ReplicateVisitor::internalVisit(StuntDouble *sd) {
    GenericData *          data;
    AtomData *             atomData;

    //if there is not atom data, just skip it
    data = sd->getPropertyByName("ATOMDATA");

    if (data != NULL) {
      atomData = dynamic_cast<AtomData *>(data);

      if (atomData == NULL) {
	return;
      }
    } else {
      return;
    }

    Snapshot* currSnapshot = info->getSnapshotManager()->getCurrentSnapshot();
    Mat3x3d box = currSnapshot->getHmat();

    std::vector<AtomInfo *> atomInfoList = atomData->getData();

    replicate(atomInfoList, atomData, box);
  }

  void ReplicateVisitor::replicate(std::vector<AtomInfo *>&infoList, AtomData *data, const Mat3x3d& box) {
    AtomInfo* newAtomInfo;
    std::vector<Vector3i>::iterator dirIter;
    std::vector<AtomInfo *>::iterator i;

    for( dirIter = dir.begin(); dirIter != dir.end(); ++dirIter ) {
      for( i = infoList.begin(); i != infoList.end(); i++ ) {
	newAtomInfo = new AtomInfo();
	*newAtomInfo = *(*i);

	for( int j = 0; j < 3; j++ )
	  newAtomInfo->pos[j] += (*dirIter)[0]*box(j, 0) + (*dirIter)[1]*box(j, 1) + (*dirIter)[2]*box(j, 2);

	data->addAtomInfo(newAtomInfo);
      }
    } // end for(dirIter)  
  }

  const std::string ReplicateVisitor::toString() {
    char                            buffer[65535];
    std::string                     result;
    std::set<std::string>::iterator i;

    sprintf(buffer,
            "------------------------------------------------------------------\n");
    result += buffer;

    sprintf(buffer, "Visitor name: %s\n", visitorName.c_str());
    result += buffer;

    sprintf(buffer,
            "Visitor Description: replicate the atoms in different direction\n");
    result += buffer;

    //print the replicate direction
    sprintf(buffer, "repeatX = %d:\n", replicateOpt[0]);
    result += buffer;

    sprintf(buffer, "repeatY = %d:\n", replicateOpt[1]);
    result += buffer;

    sprintf(buffer, "repeatZ = %d:\n", replicateOpt[2]);
    result += buffer;

    sprintf(buffer,
            "------------------------------------------------------------------\n");
    result += buffer;

    return result;
  }

  //------------------------------------------------------------------------//

  XYZVisitor::XYZVisitor(SimInfo *info) : BaseVisitor(), seleMan(info), 
                                          evaluator(info), doPositions_(true),
                                          doVelocities_(false), 
                                          doForces_(false), doVectors_(false),
                                          doCharges_(false) {
    this->info = info;
    visitorName = "XYZVisitor";
    
    evaluator.loadScriptString("select all");
    
    if (!evaluator.isDynamic()) {
      seleMan.setSelectionSet(evaluator.evaluate());
    }
  }
  
  XYZVisitor::XYZVisitor(SimInfo *info, const std::string& script) :
    BaseVisitor(), seleMan(info), evaluator(info), doPositions_(true),
    doVelocities_(false), doForces_(false), doVectors_(false),
    doCharges_(false) {
    
    this->info = info;
    visitorName = "XYZVisitor";
    
    evaluator.loadScriptString(script);
    
    if (!evaluator.isDynamic()) {
      seleMan.setSelectionSet(evaluator.evaluate());
    }
  }
    
  void XYZVisitor::visit(Atom *atom) {
    if (isSelected(atom))
      internalVisit(atom);
  }

  void XYZVisitor::visit(DirectionalAtom *datom) {
    if (isSelected(datom))
      internalVisit(datom);
  }
  
  void XYZVisitor::visit(RigidBody *rb) {
    if (isSelected(rb))
      internalVisit(rb);
  }
  
  void XYZVisitor::update() {
    //if dynamic, we need to re-evaluate the selection
    if (evaluator.isDynamic()) {
      seleMan.setSelectionSet(evaluator.evaluate());
    }
  }
  
  void XYZVisitor::internalVisit(StuntDouble *sd) {
    GenericData *                     data;
    AtomData *                        atomData;
    AtomInfo *                        atomInfo;
    std::vector<AtomInfo *>::iterator i;
    char                              buffer[1024];
    
    //if there is not atom data, just skip it
    data = sd->getPropertyByName("ATOMDATA");
    
    if (data != NULL) {
      atomData = dynamic_cast<AtomData *>(data);
      
      if (atomData == NULL)
	return;
    } else
      return;

    for( atomInfo = atomData->beginAtomInfo(i); atomInfo;
         atomInfo = atomData->nextAtomInfo(i) ) {
     
      std::string line;
      sprintf(buffer, "%s", atomInfo->atomTypeName.c_str());
      line += buffer;
      
      if (doPositions_){
        sprintf(buffer, "%15.8f%15.8f%15.8f", atomInfo->pos[0], 
                atomInfo->pos[1], atomInfo->pos[2]);
        line += buffer;
      }      
      if (doCharges_ && atomInfo->hasCharge) {
        sprintf(buffer, "%15.8f", atomInfo->charge);
        line += buffer;
      }
      if (doVectors_ && atomInfo->hasVector) {
        sprintf(buffer, "%15.8f%15.8f%15.8f", atomInfo->vec[0], 
                atomInfo->vec[1], atomInfo->vec[2]);
        line += buffer;
      }
      if (doVelocities_ && atomInfo->hasVelocity) {
        sprintf(buffer, "%15.8f%15.8f%15.8f", atomInfo->vel[0], 
                atomInfo->vel[1], atomInfo->vel[2]);
        line += buffer;
      }
      if (doForces_ && atomInfo->hasForce) {
        sprintf(buffer, "%15.8f%15.8f%15.8f", atomInfo->frc[0], 
                atomInfo->frc[1], atomInfo->frc[2]);
        line += buffer;
      }      
      frame.push_back(line);
    }    
  }

  bool XYZVisitor::isSelected(StuntDouble *sd) {
    return seleMan.isSelected(sd);
  }

  void XYZVisitor::writeFrame(std::ostream &outStream) {
    std::vector<std::string>::iterator i;
    char buffer[1024];
    
    if (frame.size() == 0)
      std::cerr << "Current Frame does not contain any atoms" << std::endl;
    
    //total number of atoms  
    outStream << frame.size() << std::endl;
    
    //write comment line
    Snapshot* currSnapshot = info->getSnapshotManager()->getCurrentSnapshot();
    Mat3x3d box = currSnapshot->getHmat();
    
    sprintf(buffer,
            "%15.8f;%15.8f%15.8f%15.8f;%15.8f%15.8f%15.8f;%15.8f%15.8f%15.8f",
            currSnapshot->getTime(),
            box(0, 0), box(0, 1), box(0, 2),
            box(1, 0), box(1, 1), box(1, 2),
            box(2, 0), box(2, 1), box(2, 2));
    
    outStream << buffer << std::endl;
    
    for( i = frame.begin(); i != frame.end(); ++i )
      outStream << *i << std::endl;
  }
  
  std::string XYZVisitor::trimmedName(const std::string&atomTypeName) {    
    return atomTypeName.substr(0, atomTypeName.find('-'));
  }
  
  const std::string XYZVisitor::toString() {
    char        buffer[65535];
    std::string result;
    
    sprintf(buffer,
            "------------------------------------------------------------------\n");
    result += buffer;
    
    sprintf(buffer, "Visitor name: %s\n", visitorName.c_str());
    result += buffer;
    
    sprintf(buffer,
            "Visitor Description: assemble the atom data and output xyz file\n");
    result += buffer;

    sprintf(buffer,
            "------------------------------------------------------------------\n");
    result += buffer;
    
    return result;
  }
  
  //----------------------------------------------------------------------------//
  
  void PrepareVisitor::internalVisit(Atom *atom) {
    GenericData *data;
    AtomData *   atomData;
    
    //if visited property is  existed, remove it
    data = atom->getPropertyByName("VISITED");
    
    if (data != NULL) {
      atom->removeProperty("VISITED");
    }
    
    //remove atomdata
    data = atom->getPropertyByName("ATOMDATA");

    if (data != NULL) {
      atomData = dynamic_cast<AtomData *>(data);

      if (atomData != NULL)
	atom->removeProperty("ATOMDATA");
    }
  }

  void PrepareVisitor::internalVisit(RigidBody *rb) {
    GenericData* data;
    AtomData* atomData;
    std::vector<Atom *> myAtoms;
    std::vector<Atom *>::iterator atomIter;

    //if visited property is  existed, remove it
    data = rb->getPropertyByName("VISITED");

    if (data != NULL) {
      rb->removeProperty("VISITED");
    }

    //remove atomdata
    data = rb->getPropertyByName("ATOMDATA");

    if (data != NULL) {
      atomData = dynamic_cast<AtomData *>(data);

      if (atomData != NULL)
	rb->removeProperty("ATOMDATA");
    }

    myAtoms = rb->getAtoms();

    for( atomIter = myAtoms.begin(); atomIter != myAtoms.end(); ++atomIter )
      internalVisit(*atomIter);
  }

  const std::string PrepareVisitor::toString() {
    char buffer[65535];
    std::string result;
 
    sprintf(buffer,
	    "------------------------------------------------------------------\n");
    result += buffer;
 
    sprintf(buffer, "Visitor name: %s", visitorName.c_str());
    result += buffer;
 
    sprintf(buffer,
	    "Visitor Description: prepare for operation of other vistors\n");
    result += buffer;
 
    sprintf(buffer,
	    "------------------------------------------------------------------\n");
    result += buffer;
 
    return result;
  }

  //----------------------------------------------------------------------------//

  WaterTypeVisitor::WaterTypeVisitor() {
    visitorName = "WaterTypeVisitor";
    waterTypeList.insert("TIP3P_RB_0");
    waterTypeList.insert("TIP4P_RB_0");
    waterTypeList.insert("TIP4P-Ew_RB_0");
    waterTypeList.insert("TIP5P_RB_0");
    waterTypeList.insert("TIP5P-E_RB_0");
    waterTypeList.insert("SPCE_RB_0");
    waterTypeList.insert("SPC_RB_0");
  }

  void WaterTypeVisitor::visit(RigidBody *rb) {
    std::string rbName;
    std::vector<Atom *> myAtoms;
    std::vector<Atom *>::iterator atomIter;
    GenericData* data;
    AtomData* atomData;
    AtomInfo* atomInfo;
    std::vector<AtomInfo *>::iterator i;

    rbName = rb->getType();
    
    if (waterTypeList.find(rbName) != waterTypeList.end()) {
      myAtoms = rb->getAtoms();

      for( atomIter = myAtoms.begin(); atomIter != myAtoms.end();
	   ++atomIter ) {
	data = (*atomIter)->getPropertyByName("ATOMDATA");
        
	if (data != NULL) {
	  atomData = dynamic_cast<AtomData *>(data);
          
	  if (atomData == NULL)
	    continue;
	} else
	  continue;
        
	for( atomInfo = atomData->beginAtomInfo(i); atomInfo;
	     atomInfo = atomData->nextAtomInfo(i) ) {
	  atomInfo->atomTypeName = trimmedName(atomInfo->atomTypeName);
	} //end for(atomInfo)
      }     //end for(atomIter)
    }         //end if (waterTypeList.find(rbName) != waterTypeList.end())
  }

  std::string WaterTypeVisitor::trimmedName(const std::string&atomTypeName) {
    return atomTypeName.substr(0, atomTypeName.find('_'));
  }

  const std::string WaterTypeVisitor::toString() {
    char buffer[65535];
    std::string result;

    sprintf(buffer,
            "------------------------------------------------------------------\n");
    result += buffer;

    sprintf(buffer, "Visitor name: %s\n", visitorName.c_str());
    result += buffer;

    sprintf(buffer,
            "Visitor Description: Replace the atom type in water model\n");
    result += buffer;

    sprintf(buffer,
            "------------------------------------------------------------------\n");
    result += buffer;

    return result;
  }

} //namespace OpenMD
