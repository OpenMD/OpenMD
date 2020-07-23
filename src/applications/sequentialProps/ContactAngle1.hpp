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

#ifndef APPLICATIONS_SEQUENTIALPROPS_CONTACTANGLE1_HPP
#define APPLICATIONS_SEQUENTIALPROPS_CONTACTANGLE1_HPP
#include "applications/sequentialProps/SequentialAnalyzer.hpp"

using namespace std;
namespace OpenMD {

  /// Calculates the contact angle of a droplet with a surface
  /// using a spherical cap approximation for the droplet.

  /**  The position of the spherical cap relative to the surface plane
       is determined by the center-of-mass position of the selection,
       and this method assumes a uniform density in the droplet.  The
       angle of intersection between the surface of the spherical cap
       and the plane defines the contact angle, which is related to
       the center of mass height by:

       \f$ z_\mathrm{cm} = (2)^{-4/3} R_0 \left(
       \frac{1-\cos\theta}{2+\cos\theta}\right)^{1/3}
       \frac{3+\cos\theta}{2+\cos\theta} \f$

       where \f$z_\mathrm{cm}\f$ is measured relative to the planar
       surface, and \f$R_0\f$ is the radius of the free spherical
       droplet.

       This method was first proposed in:

       J. Hautman and M.L. Klein, Phys. Rev. Lett. 67(13), 1763 (1991).
       DOI: 10.1103/PhysRevLett.67.1763

       This Analyzer requires statement of the reference height of the
       solid surface, solidZ, and \f$R_0\f$, the dropletRadius.
 
  */ 
  class ContactAngle1 : public SequentialAnalyzer{
  public:
    ContactAngle1(SimInfo* info, const std::string& filename, 
                  const std::string& sele1, const std::string& sele2,
                  RealType solidZ, RealType dropletRadius);

    virtual void doFrame(int frame);
    
  private:
    RealType solidZ_;
    RealType dropletRadius_;       
  };
}

#endif



