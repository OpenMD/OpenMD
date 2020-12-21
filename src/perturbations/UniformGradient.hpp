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
 

/*! \file perturbations/UniformGradient.hpp
    \brief Uniform Electric Field Gradient perturbation
*/

#ifndef PERTURBATIONS_UNIFORMGRADIENT_HPP
#define PERTURBATIONS_UNIFORMGRADIENT_HPP

#include "perturbations/Perturbation.hpp"
#include "brains/SimInfo.hpp"

namespace OpenMD {
   
  //! Applies a uniform electric field gradient to the system
  /*! The gradient is applied as an external perturbation.  The user specifies

  \code{.unparsed}
    uniformGradientStrength = c;
    uniformGradientDirection1 = (a1, a2, a3)
    uniformGradientDirection2 = (b1, b2, b3);
  \endcode 
    
    in the .omd file where the two direction vectors, \f$ \mathbf{a} \f$ 
    and \f$ \mathbf{b} \f$ are unit vectors, and the value of \f$ g \f$ 
    is in units of \f$ V / \AA^2 \f$

    The electrostatic potential corresponding to this uniform gradient is

    \f$ \phi(\mathbf{r})  = - \frac{g}{2} \left[ 
    \left(a_1 b_1 - \frac{\cos\psi}{3}\right) x^2
    + (a_1 b_2 + a_2 b_1) x y + (a_1 b_3 + a_3 b_1) x z +
    + (a_2 b_1 + a_1 b_2) y x 
    + \left(a_2 b_2 - \frac{\cos\psi}{3}\right) y^2
    + (a_2 b_3 + a_3 b_2) y z + (a_3 b_1 + a_1 b_3) z x 
    + (a_3 b_2 + a_2 b_3) z y
    + \left(a_3 b_3 - \frac{\cos\psi}{3}\right) z^2 \right] \f$

    where \f$ \cos \psi = \mathbf{a} \cdot \mathbf{b} \f$.  Note that
    this potential grows unbounded and is not periodic.  For these reasons,
    care should be taken in using a Uniform Gradient with point charges.

    The corresponding field is:

    \f$ \mathbf{E} = \frac{g}{2} \left( \begin{array}{c} 
    2\left(a_1 b_1 - \frac{\cos\psi}{3}\right) x +  (a_1 b_2 + a_2 b_1) y 
    + (a_1 b_3 + a_3 b_1) z  \\
    (a_2 b_1 + a_1 b_2)  x + 2 \left(a_2 b_2 - \frac{\cos\psi}{3}\right) y 
    + (a_2 b_3 + a_3 b_2) z \\
    (a_3 b_1 + a_1 b_3) x +  (a_3 b_2 + a_2 b_3) y 
    + 2 \left(a_3 b_3 - \frac{\cos\psi}{3}\right) z \end{array} \right) \f$

    The field also grows unbounded and is not periodic.  For these reasons,
    care should be taken in using a Uniform Gradient with point dipoles.

    The corresponding field gradient is:

    \f$ \nabla \mathbf{E} = \frac{g}{2} \left( \begin{array}{ccc}  
    2\left(a_1 b_1 - \frac{\cos\psi}{3}\right)  &  
    (a_1 b_2 + a_2 b_1) & (a_1 b_3 + a_3 b_1) \\
    (a_2 b_1 + a_1 b_2)  & 2 \left(a_2 b_2 - \frac{\cos\psi}{3}\right) & 
    (a_2 b_3 + a_3 b_2) \\
    (a_3 b_1 + a_1 b_3) & (a_3 b_2 + a_2 b_3) & 
    2 \left(a_3 b_3 - \frac{\cos\psi}{3}\right) \end{array} \right) \f$

    which is uniform everywhere.

    The uniform field gradient applies a force on charged atoms, 
    \f$ \mathbf{F} = C \mathbf{E}(\mathbf{r}) \f$.  
    For dipolar atoms, the gradient applies both a potential, 
    \f$ U = -\mathbf{D} \cdot \mathbf{E}(\mathbf{r}) \f$, a force, 
    \f$ \mathbf{F} = \mathbf{D} \cdot \nabla \mathbf{E} \f$, and a torque,
    \f$ \mathbf{\tau} = \mathbf{D} \times \mathbf{E}(\mathbf{r}) \f$.
    
    For quadrupolar atoms, the uniform field gradient exerts a potential, 
    \f$ U = - \mathsf{Q}:\nabla \mathbf{E} \f$, and a torque 
    \f$ \mathbf{F} = 2 \mathsf{Q} \times \nabla \mathbf{E} \f$
    
  */ 
  class UniformGradient : public Perturbation {
    

  public:
    UniformGradient(SimInfo* info);
    
  protected:
    virtual void initialize();
    virtual void applyPerturbation();
    
  private:
    bool initialized;
    bool doUniformGradient;
    bool doParticlePot;
    Globals* simParams;
    SimInfo* info_ {nullptr};
    Mat3x3d Grad_;
    Vector3d a_, b_;
    RealType g_, cpsi_;    
  };


} //end namespace OpenMD
#endif

