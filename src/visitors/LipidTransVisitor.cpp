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

#include "visitors/LipidTransVisitor.hpp"
#include "types/MultipoleAdapter.hpp"
#include "utils/simError.h"

namespace OpenMD {
  LipidTransVisitor::LipidTransVisitor(SimInfo* info, const std::string& originSeleScript, const std::string& refSeleScript) 
    : BaseVisitor(), info_(info), originEvaluator_(info), originSeleMan_(info), refEvaluator_(info), refSeleMan_(info), refSd_(NULL) {

 
      visitorName = "LipidTransVisitor";
    
      originEvaluator_.loadScriptString(originSeleScript);            
      if (!originEvaluator_.isDynamic()) {  
        originSeleMan_.setSelectionSet(originEvaluator_.evaluate());
        if (originSeleMan_.getSelectionCount() == 1) {
	  int i;
	  originDatom_ = dynamic_cast<DirectionalAtom*>(originSeleMan_.beginSelected(i));
	  if (originDatom_ ==  NULL) {
	    sprintf(painCave.errMsg, "LipidTransVisitor: origin selection must select an directional atom");
	    painCave.isFatal = 1;
	    simError();                  
	  }
        } else {
	  sprintf(painCave.errMsg, "LipidTransVisitor: origin selection must select an directional atom");
	  painCave.isFatal = 1;
	  simError();                  
            
        }
      }

      refEvaluator_.loadScriptString(refSeleScript);
      if (!refEvaluator_.isDynamic()) {  
        refSeleMan_.setSelectionSet(refEvaluator_.evaluate());
        if (refSeleMan_.getSelectionCount() == 1) {
	  int i;
	  refSd_ = refSeleMan_.beginSelected(i);
            
        } else {
	  //error
            
        }
      }
    
    }

  void LipidTransVisitor::update() {

    Vector3d ref = refSd_->getPos();
    origin_ = originDatom_->getPos();
    Vector3d v1 =  ref - origin_;
    info_->getSnapshotManager()->getCurrentSnapshot()->wrapVector(v1);

    MultipoleAdapter ma = MultipoleAdapter(originDatom_->getAtomType());
    Vector3d zaxis;
    if (ma.isDipole() ) {
      zaxis = originDatom_->getDipole();
    } else {
      zaxis = originDatom_->getA().transpose()*V3Z;
    }

    Vector3d xaxis = cross(v1, zaxis);
    Vector3d yaxis = cross(zaxis, xaxis);

    xaxis.normalize();
    yaxis.normalize();
    zaxis.normalize();
        
    rotMat_.setRow(0, xaxis);
    rotMat_.setRow(1, yaxis);
    rotMat_.setRow(2, zaxis);


  }

  void LipidTransVisitor::internalVisit(StuntDouble *sd) {
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

    Snapshot* currSnapshot = info_->getSnapshotManager()->getCurrentSnapshot();
    
    for( atomInfo = atomData->beginAtomInfo(i); atomInfo; atomInfo = atomData->nextAtomInfo(i) ) {

      Vector3d tmp= atomInfo->pos - origin_;
      currSnapshot->wrapVector(tmp);
      atomInfo->pos = rotMat_ * tmp;;
      atomInfo->vec = rotMat_ * atomInfo->vec;
    }
  }

  const std::string LipidTransVisitor::toString() {
    char        buffer[65535];
    std::string result;

    sprintf(buffer,
            "------------------------------------------------------------------\n");
    result += buffer;

    sprintf(buffer, "Visitor name: %s\n", visitorName.c_str());
    result += buffer;

    sprintf(buffer,
            "Visitor Description: rotate the whole system\n");
    result += buffer;

    sprintf(buffer,
            "------------------------------------------------------------------\n");
    result += buffer;

    return result;
  }

}
