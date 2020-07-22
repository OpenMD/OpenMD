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

#ifndef APPLICATIONS_DYNAMICPROPS_COLLECTIVEDIPOLEDISPLACEMENT_HPP
#define APPLICATIONS_DYNAMICPROPS_COLLECTIVEDIPOLEDISPLACEMENT_HPP

#include "applications/dynamicProps/TimeCorrFunc.hpp"
#include "brains/Thermo.hpp"

namespace OpenMD {
  //! Calculates the collective dipole displacement function
  /*! This time correlation function is the Helfand moment conjugate
      to the current density. Helfand moments are used to calculate
      the Einstein-Helfand relations for transport that are formally
      equivalent to Green-Kubo expressions using a related flux.  In
      this case, the flux,

      \f[ \mathbf{J}(t) = \sum_{i=1}^{N} q_i \mathbf{v}_{\mathrm{cm},i}(t) \f] 

      is normally used to calculate an ionic conductivity, 

      \f[ \sigma = \frac{1}{3V k_b T} \int_0^\infty \left< \mathbf{J}(0) \cdot \mathbf{J}(t) \right> dt \f]

      The cm subscript denotes center of mass locations for all molecules.

      This class computes the collective translational dipole moment,

      \f[ \mathbf{M}_\mathrm{trans}(t) = \sum_{i=1}^{N} q_i \mathbf{r}_{\mathrm{cm},i}(t) \f]

      as well as total contributions to the system's net dipole moment

      \f[ \mathbf{M}_\mathrm{tot}(t) =  \sum_{i=1}^{N} \sum_{a} q_{ia} \mathbf{r}_{ia}(t) = \sum_{i=1}^{N} q_i \mathbf{r}_{\mathrm{cq},i}(t) \f]

      where cq denotes the molecular center of charge.  It also
      calculates the rotational contribution,

      \f[ \mathbf{M}_\mathrm{rot}(t) = \sum_{i=1}^{N} q_i \left[ \mathbf{r}_{\mathrm{cq},i}(t) - \mathbf{r}_{\mathrm{cm},i}(t) \right] \f]

      The correlation functions are the displacements of these terms
      from their values at an earlier time,

      \f[ \left< \left| \mathbf{M}_\mathrm{trans}(t) - \mathbf{M}_\mathrm{trans}(0) \right|^2 \right> \f]

      and identical quantities for the total and rotational contributions.
  */ 
  class CollectiveDipoleDisplacement : public SystemACF<Vector3d> {
  public:
    CollectiveDipoleDisplacement(SimInfo* info, const std::string& filename,
                                 const std::string& sele1,
                                 const std::string& sele2);
    
  private:
    virtual void computeProperty1(int frame);
    virtual Vector3d calcCorrVal(int frame1, int frame2);
    
    Thermo* thermo_;
    
    std::vector<Vector3d> CRcm_;
    std::vector<Vector3d> CRtot_;
    std::vector<Vector3d> CRrot_;
  };
}
#endif
