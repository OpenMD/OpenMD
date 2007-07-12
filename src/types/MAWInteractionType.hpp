/*
 * Copyright (c) 2007 The University of Notre Dame. All Rights Reserved.
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
 
#ifndef TYPES_MAWINTERACTIONTYPE_HPP
#define TYPES_MAWINTERACTIONTYPE_HPP

#include "types/NonBondedInteractionType.hpp"

namespace oopse {
  /**
   * @class MAWInteractionType 
   *
   
   * MAWInteractionType (Metal-Angular-Water) is one of the basic
   * Metal-to-NonMetal interaction types.
   *
   * Formula is V = D0*exp(-beta0*(r-r0))*(exp(-beta0(r-r0)-2) +
   *            gamma*D0*(exp(-betaH*(r-r0))*(1+alpha*cos(theta))*sin(phi)^2
   *
   * The spherical coordinates are defined in the body-fixed frame
   * of a rigid-body water molecule (HO bonds are on the Y-Z plane)
   * and the dipole vector of the water molecule points along the
   * Z-axis.  A metal atom's position is uniquely defined by a set
   * of spherical polar coordinates (r, theta, phi) in the
   * body-fixed frame of each water molecule.
   */
  class MAWInteractionType : public NonBondedInteractionType {
    
  public:
    
    MAWInteractionType(RealType myD0, RealType myBeta0, RealType myR0,
		       RealType myBetaH, RealType myGamma, RealType myAlpha){
      D0 = myD0;
      beta0 = myBeta0;
      r0 = myR0;
      betaH = myBetaH;
      gamma = myGamma;
      alpha = myAlpha;
    }
    
    virtual void tellFortran(int atid1, int atid2) {
      mnmit.MNMInteractionType = MNM_MAW;
      mnmit.metal_atid = atid1;
      mnmit.nonmetal_atid = atid2;
      mnmit.R0 = r0;
      mnmit.D0 = D0;
      mnmit.beta0 = beta0;
      mnmit.betaH = betaH;
      mnmit.alpha = alpha;
      mnmit.gamma = gamma;
      
      addMNMInteraction(&mnmit);
    }
    
  private:    
    RealType D0;
    RealType beta0;
    RealType r0;
    RealType betaH;
    RealType gamma;
    RealType alpha;    
  };
}
#endif
