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

#include <algorithm>
#include <functional>
#include "applications/sequentialProps/CenterOfMass.hpp"
#include "utils/simError.h"
#include "io/DumpReader.hpp"
#include "primitives/Molecule.hpp"
#include "utils/NumericConstant.hpp"
namespace OpenMD {

  CenterOfMass::CenterOfMass(SimInfo* info, const std::string& filename, 
                                   const std::string& sele)
    : SequentialAnalyzer(info, filename), selectionScript_(sele), 
      seleMan_(info), evaluator_(info) {
    
    setOutputName(getPrefix(filename) + ".dens");
    
    evaluator_.loadScriptString(sele);
    
    if (!evaluator_.isDynamic()) {
      seleMan_.setSelectionSet(evaluator_.evaluate());
    }            
  }

  void CenterOfMass::doFrame() {
    StuntDouble* sd;
    int i;
    
    if (evaluator_.isDynamic()) {
	seleMan_.setSelectionSet(evaluator_.evaluate());
    }


    RealType mtot = 0.0;
    Vector3d com(V3Zero);
    RealType mass;
    
    for (sd = seleMan_.beginSelected(i); sd != NULL;
         sd = seleMan_.nextSelected(i)) {      
      mass = sd->getMass();
      mtot += mass;
      com += sd->getPos() * mass;
    }

    com /= mtot;
        
    values_.push_back( com );
  }
  
  void CenterOfMass::writeSequence() {
    std::ofstream ofs(outputFilename_.c_str(), std::ios::binary);
    
    if (ofs.is_open()) {

      ofs << "#" << getSequenceType() << "\n";
      ofs << "#extra information: " << extra_ << "\n";
      ofs << "#time\tvalue\n";

      for (unsigned int i = 0; i < times_.size(); ++i) {
        ofs << times_[i]
            << "\t" << values_[i].x()
            << "\t" << values_[i].y()
            << "\t" << values_[i].z()
            << "\n";
      }
      
    } else {
      sprintf(painCave.errMsg,
              "CenterOfMass::writeSequence Error: fail to open %s\n", 
              outputFilename_.c_str());
      painCave.isFatal = 1;
      simError();        
    }
    
    ofs.close();    
  }

}


