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
 
#include <cmath>
#include "visitors/ZconsVisitor.hpp"
#include "primitives/Molecule.hpp"
#include "utils/StringUtils.hpp"
#include "types/ZconsStamp.hpp"
namespace OpenMD {

  ZConsVisitor::ZConsVisitor(SimInfo* info) : BaseVisitor(), info_(info), zconsReader_(NULL){

    visitorName = "ZConsVisitor";
    currSnapshot_ = info_->getSnapshotManager()->getCurrentSnapshot();
    Globals* simParam = info_->getSimParams();

    if (simParam->haveZconsTime()){
      zconsTime_ = simParam->getZconsTime();
    }
    else{
      sprintf(painCave.errMsg,
	      "ZConstraint error: If you use a ZConstraint,\n"
	      "\tyou must set zconsTime.\n");
      painCave.isFatal = 1;
      simError();
    }

    if (simParam->haveZconsTol()){
      zconsTol_ = simParam->getZconsTol();
    }
    else{
      zconsTol_ = 0.01;
      sprintf(painCave.errMsg,
	      "ZConstraint Warning: Tolerance for z-constraint method is not specified.\n"
	      "\tOpenMD will use a default value of %f.\n"
	      "\tTo set the tolerance, use the zconsTol variable.\n",
	      zconsTol_);
      painCave.isFatal = 0;
      simError();      
    }    
         
    int nZconstraints = simParam->getNZconsStamps();
    std::vector<ZConsStamp*> stamp = simParam->getZconsStamps();
    for (int i = 0; i < nZconstraints; i++){
      int zmolIndex = stamp[i]->getMolIndex();
      zmolStates_.insert(std::make_pair(zmolIndex, zsMoving));
    }


    //fill zatomToZmol_ array
    /** @todo only works for single version now*/
    std::map<int, ZConsState>::iterator j;
    for (j = zmolStates_.begin(); j != zmolStates_.end(); ++j) {
      Molecule* mol = info_->getMoleculeByGlobalIndex(j->first);
      assert(mol != NULL);
      Molecule::AtomIterator ai;
      Atom* at;
      for (at = mol->beginAtom(ai); at != NULL; at = mol->nextAtom(ai)) {
	zatomToZmol_.insert(std::make_pair(at->getGlobalIndex(), mol->getGlobalIndex()));
      }
    }

    zconsFilename_ = getPrefix(info_->getFinalConfigFileName()) + ".fz"; 
    
    zconsReader_ = new ZConsReader(info);

    if (zconsReader_->hasNextFrame())
      zconsReader_->readNextFrame();
  
  }

  ZConsVisitor::~ZConsVisitor(){
    if(!zconsReader_) 
      delete zconsReader_;
  
  }

  void ZConsVisitor::visit(Atom* atom){
    std::string prefix;
    if(isZconstraint(atom->getGlobalIndex(), prefix))
      internalVisit(atom, prefix);
  }

  void ZConsVisitor::visit(DirectionalAtom* datom){
    std::string prefix;

    if(isZconstraint(datom->getGlobalIndex(), prefix))
      internalVisit(datom, prefix);
  }

  void ZConsVisitor::visit(RigidBody* rb){
    std::string prefix;
    std::vector<Atom*> atoms;

    atoms = rb->getAtoms();

    if(isZconstraint(atoms[0]->getGlobalIndex(), prefix))
      internalVisit(rb, prefix);
  }

  void ZConsVisitor::update(){
    Vector3d com;
    std::map<int, ZConsState>::iterator i;
    for ( i = zmolStates_.begin(); i != zmolStates_.end(); ++i) {
      i->second = zsMoving;
    }
     
    readZconsFile(currSnapshot_->getTime());

    const std::vector<ZconsData>& fixedZmolData = zconsReader_->getFixedZMolData();
    std::vector<ZconsData>::const_iterator j;
    for (j = fixedZmolData.begin(); j != fixedZmolData.end(); ++j) {
      std::map<int, ZConsState>::iterator k = zmolStates_.find(j->zmolIndex);
      assert(k != zmolStates_.end());
      k->second = zsFixed;
    }
    
  }

  void ZConsVisitor::readZconsFile(RealType time) {
    RealType tempTime;
    while(zconsReader_->hasNextFrame()){
      tempTime = zconsReader_->getCurTime();
      if(tempTime >= time) {
	return;
      }
        
      zconsReader_->readNextFrame();
    } 
  }

  void ZConsVisitor::internalVisit(StuntDouble* sd, const std::string& prefix){
    GenericData* data;
    AtomData* atomData;
    AtomInfo* atomInfo;
    std::vector<AtomInfo*>::iterator iter;

    //if there is not atom data, just skip it
    data = sd->getPropertyByName("ATOMDATA");
    if(data != NULL){
      atomData = dynamic_cast<AtomData*>(data);  
      if(atomData == NULL)
	return;
    }
    else
      return;

    for(atomInfo  = atomData->beginAtomInfo(iter); atomInfo; atomInfo = atomData->nextAtomInfo(iter))
      (atomInfo->atomTypeName).insert(0, prefix);
  }


  bool ZConsVisitor::isZconstraint(int atomIndex, std::string& prefix){
    std::string prefixString[] = {"ZF", "ZM"};
    std::map<int, int>::iterator i = zatomToZmol_.find(atomIndex);
    if (i ==  zatomToZmol_.end() ){
      prefix = "";
      return false;
    } else {

      std::map<int, ZConsState>::iterator j = zmolStates_.find(i->second);
      assert(j !=zmolStates_.end());
      prefix = prefixString[j->second];
      return true;
    }
  }

  const std::string ZConsVisitor::toString(){
    char buffer[65535];
    std::string result;

    sprintf(buffer ,"------------------------------------------------------------------\n");
    result += buffer;

    sprintf(buffer ,"Visitor name: %s\n", visitorName.c_str());
    result += buffer;

    sprintf(buffer , "number of zconstraint molecule: %d\n", (int) zmolStates_.size());
    result += buffer;

    sprintf(buffer , "zconstraint tolerance = %lf\n", zconsTol_);
    result += buffer;

    sprintf(buffer , "zconstraint sample time = %lf\n", zconsTime_);
    result += buffer;

    sprintf(buffer , "zconstraint output filename = %s\n", zconsFilename_.c_str());
    result += buffer;

    std::map<int, ZConsState>::iterator i;
    int j = 0;
    for ( i = zmolStates_.begin(); i != zmolStates_.end(); ++i) {
      sprintf(buffer ,"zconstraint molecule[%d] = %d\n", j++, i->first);
      result += buffer;
    }
  
    sprintf(buffer ,"------------------------------------------------------------------\n");
    result += buffer;

    return result;
  }


}//namespace OpenMD
