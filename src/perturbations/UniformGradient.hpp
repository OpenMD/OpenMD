/*
 * Copyright (c) 2014 The University of Notre Dame. All Rights Reserved.
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
 

/*! \file perturbations/UniformGradient.hpp
    \brief Uniform Electric Field Gradient perturbation
*/

#ifndef PERTURBATIONS_UNIFORMGRADIENT_HPP
#define PERTURBATIONS_UNIFORMGRADIENT_HPP

#include "perturbations/Perturbation.hpp"
#include "brains/SimInfo.hpp"

namespace OpenMD {

  struct UniGradPars {
    RealType a;
    RealType b;
    RealType c;
    RealType alpha;
    RealType beta;
  };
   
  //! Applies a uniform electric field gradient to the system
  /*! The gradient is applied as an external perturbation.  The user specifies

  \code{.unparsed}
    uniformGradient = (a, b, c, alpha, beta);
  \endcode 
    
    in the .md file where the values of a, b, c, alpha, beta are in units of
    \f$ V / \AA^2 \f$

    The electrostatic potential corresponding to this uniform gradient is

    \f$ \phi(\mathbf{r})  = - a x y - b x z - c y z - \alpha x^2 / 2 - \beta y^2 / 2 + (\alpha + \beta) z^2 / 2 \f$

    which grows unbounded and is not periodic.  For these reasons,
    care should be taken in using a Uniform Gradient with point charges.

    The corresponding field is:

    \f$ \mathbf{E} = \left( \array{c} \alpha x + a y + b z  \\ a x + \beta y + c z \\ b x + c y - (\alpha + \beta) z \end{array} \right) \f$

    The field also grows unbounded and is not periodic.  For these reasons,
    care should be taken in using a Uniform Gradient with point dipoles.

    The corresponding field gradient is:

    \f$ \nabla \mathbf{E} = \left( \array{ccc} \alpha & a & b \\ a & \beta & c \\ b & c & -(\alpha + \beta) \end{array} \right) \f$

    which is uniform everywhere.

    The uniform field gradient applies a force on charged atoms, 
    \f$ \mathbf{F} = C \mathbf{E}(\mathbf{r}) \f$.  
    For dipolar atoms, the gradient applies both a potential, 
    \f$ U = -\mathbf{D} \cdot \mathbf{E}(\mathbf{r}) \f$, a force, 
    \f$ \mathbf{F} = \mathbf{D} \cdot \nabla \mathbf{E} \f$, and a torque,
    \f$ \mathbf{\tau} = \mathbf{D} \times \mathbf{E}(\mathbf{r}) \f$.
    
    For quadrupolar atoms, the uniform field gradient exerts a potential, 
    \f$ U = - \mathsf{Q}:\nabla \mathbf{E} $\f, and a torque 
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
    SimInfo* info_;
    UniGradPars pars_;
    Mat3x3d Grad_;
    Vector3d v1_, v2_, v3_;
  };


} //end namespace OpenMD
#endif

