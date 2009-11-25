/*
 * Copyright (c) 2007 The University of Notre Dame. All Rights Reserved.
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
 * [4]  Vardeman & Gezelter, in progress (2009).                        
 */

/**
 * @file NonBondedInteractionType.hpp
 * @author    gezelter
 * @date  07/04/2007
 * @version 1.0
 */

#ifndef TYPES_NONBONDEDINTERACTIONTYPE_HPP
#define TYPES_NONBONDEDINTERACTIONTYPE_HPP

#define __OPENMD_C
#include "UseTheForce/DarkSide/fMnMInteractions.h"
#include "UseTheForce/DarkSide/MetalNonMetal_interface.h"

namespace OpenMD {
  
  /**
   * @class NonBondedInteractionType NonBondedInteractionType.hpp "types/NonBondedInteractionType.hpp" 
   *
   * NonBondedInteractionType class is responsible for keeping track
   * of parameters for some special non-bonded interactions.  No
   * calculations are done by NonBondedInteractionTypes (at least not
   * yet).
   */
  class NonBondedInteractionType {
  public:
    NonBondedInteractionType(){
      /* set all of the values in case we pass down to fortran
       * without setting some unimportant ones... 
       */
      mnmit.MNMInteractionType = 0;
      mnmit.metal_atid = -1;
      mnmit.nonmetal_atid = -1;
      mnmit.R0 = 0.0;
      mnmit.D0 = 0.0;
      mnmit.beta0 = 0.0;
      mnmit.betaH = 0.0;
      mnmit.ca1 = 0.0;
      mnmit.cb1 = 0.0;
      mnmit.sigma = 0.0;
      mnmit.epsilon = 0.0;
    }
    virtual ~NonBondedInteractionType(){}
    
    /**
     * in metal-nonmetal interactions atid1 is always the metallic atom type.
     */
    virtual void tellFortran(int atid1, int atid2){}
    
  protected:
    MNMtype mnmit;
    
  };  
  
} //end namespace OpenMD
#endif //TYPES_NONBONDEDINTERACTIONTYPE_HPP

