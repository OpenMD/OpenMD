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
 
#include "primitives/GhostBend.hpp"
#include "primitives/DirectionalAtom.hpp"
namespace oopse {

  /**@todo still a lot left to improve*/
  void GhostBend::calcForce(RealType& angle) {
    DirectionalAtom* ghostAtom = static_cast<DirectionalAtom*>(atom2_);
    
    Vector3d pos1 = atom1_->getPos();
    Vector3d pos2 = ghostAtom->getPos();
    
    Vector3d r12 = pos1 - pos2;
    RealType d12 = r12.length();
    
    RealType d12inv = 1.0 / d12;
    
    Vector3d r32 = ghostAtom->getElectroFrame().getColumn(2);
    RealType d32 = r32.length();
    
    RealType d32inv = 1.0 / d32;
    
    RealType cosTheta = dot(r12, r32) / (d12 * d32);
    
    //check roundoff     
    if (cosTheta > 1.0) {
      cosTheta = 1.0;
    } else if (cosTheta < -1.0) {
      cosTheta = -1.0;
    }
    
    RealType theta = acos(cosTheta);
    
    RealType firstDerivative;
    
    bendType_->calcForce(theta, potential_, firstDerivative);
    
    RealType sinTheta = sqrt(1.0 - cosTheta * cosTheta);
    
    if (fabs(sinTheta) < 1.0E-12) {
      sinTheta = 1.0E-12;
    }
    
    RealType commonFactor1 = -firstDerivative / sinTheta * d12inv;
    RealType commonFactor2 = -firstDerivative / sinTheta * d32inv;
    
    Vector3d force1 = commonFactor1*(r12*(d12inv*cosTheta) - r32*d32inv);
    Vector3d force3 = commonFactor2*(r32*(d32inv*cosTheta) - r12*d12inv);
    atom1_->addFrc(force1);
    ghostAtom->addFrc(-force1);
    /**@todo test correctness */
    ghostAtom->addTrq(cross(r32, force3) );    

    atom1_->addParticlePot(potential_);
    ghostAtom->addParticlePot(potential_);

    angle = theta /M_PI * 180.0;
    
  }  
} //end namespace oopse

