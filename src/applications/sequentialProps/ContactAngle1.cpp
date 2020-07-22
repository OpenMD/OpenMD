/*
 * Copyright (c) 2004-2020 The University of Notre Dame. All Rights Reserved.
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
 * [1] Meineke, et al., J. Comp. Chem. 26, 252-271 (2005).
 * [2] Fennell & Gezelter, J. Chem. Phys. 124, 234104 (2006).
 * [3] Sun, Lin & Gezelter, J. Chem. Phys. 128, 234107 (2008).
 * [4] Vardeman, Stocker & Gezelter, J. Chem. Theory Comput. 7, 834 (2011).
 * [5] Kuang & Gezelter, Mol. Phys., 110, 691-701 (2012).
 * [6] Lamichhane, Gezelter & Newman, J. Chem. Phys. 141, 134109 (2014).
 * [7] Lamichhane, Newman & Gezelter, J. Chem. Phys. 141, 134110 (2014).
 * [8] Bhattarai, Newman & Gezelter, Phys. Rev. B 99, 094106 (2019).
 */

#include <algorithm>
#include <functional>
#include <sstream>
#include "utils/Constants.hpp"
#include "applications/sequentialProps/ContactAngle1.hpp"
#include "utils/simError.h"
#include "io/DumpReader.hpp"
#include "primitives/Molecule.hpp"
#include "math/Polynomial.hpp"

namespace OpenMD {

  ContactAngle1::ContactAngle1(SimInfo* info, const std::string& filename, 
                               const std::string& sele1,
                               const std::string& sele2, RealType solidZ,
                               RealType dropletRadius)
    : SequentialAnalyzer(info, filename, sele1, sele2), solidZ_(solidZ),
      dropletRadius_(dropletRadius) {
    
    setOutputName(getPrefix(filename) + ".ca1");

    std::stringstream params;
    params << " solid Z = " << solidZ_
           << ", droplet radius = " << dropletRadius_;
    
    const std::string paramString = params.str();
    setParameterString( paramString );
  }

  void ContactAngle1::doFrame(int frame) {
    StuntDouble* sd;
    int i;
    
    if (evaluator1_.isDynamic()) {
      seleMan1_.setSelectionSet(evaluator1_.evaluate());
    }
    
    
    RealType mtot = 0.0;
    Vector3d com(V3Zero);
    RealType mass;
    
    for (sd = seleMan1_.beginSelected(i); sd != NULL;
         sd = seleMan1_.nextSelected(i)) {      
      mass = sd->getMass();
      mtot += mass;
      com += sd->getPos() * mass;
    }
    
    com /= mtot;

    RealType dz = com.z() - solidZ_;

    if (dz < 0.0) {
      sprintf(painCave.errMsg, 
              "ContactAngle1: Z-center of mass of selection, %lf, was\n"
              "\tlocated below the solid reference plane, %lf\n",
              com.z(), solidZ_);
      painCave.isFatal = 1;
      painCave.severity = OPENMD_ERROR;
      simError();
    }

    if (dz > dropletRadius_) {
      values_.push_back(180.0);
    } else {
    
      RealType k = pow(2.0, -4.0/3.0) * dropletRadius_;
      
      RealType z2 = dz*dz;
      RealType z3 = z2 * dz;
      RealType k2 = k*k;
      RealType k3 = k2*k;
      
      Polynomial<RealType> poly;
      poly.setCoefficient(4,      z3 +      k3);
      poly.setCoefficient(3,  8.0*z3 +  8.0*k3);
      poly.setCoefficient(2, 24.0*z3 + 18.0*k3);
      poly.setCoefficient(1, 32.0*z3          );
      poly.setCoefficient(0, 16.0*z3 - 27.0*k3);
      vector<RealType> realRoots = poly.FindRealRoots();

      RealType ct;
      
      vector<RealType>::iterator ri;


      RealType maxct = -1.0;
      for (ri = realRoots.begin(); ri !=realRoots.end(); ++ri) {
        ct = *ri;
        if (ct > 1.0)  ct = 1.0;
        if (ct < -1.0) ct = -1.0;

        // use the largest magnitude of ct that it finds:
        if (ct > maxct) {
          maxct = ct;
        }                  
      }
      
      values_.push_back( acos(maxct)*(180.0/Constants::PI) );
    }
  }    
}


