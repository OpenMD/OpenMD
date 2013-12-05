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
 
#include "config.h"
#include <cmath>
#include "primitives/GhostBend.hpp"
#include "primitives/DirectionalAtom.hpp"
namespace OpenMD {

  /**@todo still a lot left to improve*/
  void GhostBend::calcForce(RealType& angle, bool doParticlePot) {
    DirectionalAtom* ghostAtom = static_cast<DirectionalAtom*>(atoms_[1]);
    
    Vector3d pos1 = atoms_[0]->getPos();
    Vector3d pos2 = ghostAtom->getPos();

    Vector3d r21 = pos1 - pos2;   
    RealType d21 = r21.length();
    
    RealType d21inv = 1.0 / d21;
   
    // we need the transpose of A to get the lab fixed vector:
    Vector3d r23 = ghostAtom->getA().transpose().getColumn(2);
    RealType d23 = r23.length();
    
    RealType d23inv = 1.0 / d23;
    
    RealType cosTheta = dot(r21, r23) / (d21 * d23);

    //check roundoff     
    if (cosTheta > 1.0) {
      cosTheta = 1.0;
    } else if (cosTheta < -1.0) {
      cosTheta = -1.0;
    }
    
    RealType theta = acos(cosTheta);

    RealType dVdTheta;
    
    bendType_->calcForce(theta, potential_, dVdTheta);
    
    RealType sinTheta = sqrt(1.0 - cosTheta * cosTheta);
    
    if (fabs(sinTheta) < 1.0E-6) {
      sinTheta = 1.0E-6;
    }
    
    RealType commonFactor1 = dVdTheta / sinTheta * d21inv;
    RealType commonFactor2 = dVdTheta / sinTheta * d23inv;
    
    Vector3d force1 = commonFactor1 * (r23 * d23inv - r21*d21inv*cosTheta);
    Vector3d force3 = commonFactor2 * (r21 * d21inv - r23*d23inv*cosTheta);

    // Total force in current bend is zero

    atoms_[0]->addFrc(force1);
    ghostAtom->addFrc(-force1);

    ghostAtom->addTrq( cross(r23, force3) );    
    if(doParticlePot) {
      atoms_[0]->addParticlePot(potential_);
      ghostAtom->addParticlePot(potential_);
    }

    angle = theta /M_PI * 180.0;
   
  }  
} //end namespace OpenMD

