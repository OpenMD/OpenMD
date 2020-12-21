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

/*! \file FreqFlucCorrFunc.hpp
  \brief Frequency Fluctuation Correlation Function
*/

#ifndef APPLICATIONS_DYNAMICPROPS_FREQFLUCCORRFUNC_HPP
#define APPLICATIONS_DYNAMICPROPS_FREQFLUCCORRFUNC_HPP

#include <string>
#include <vector>

#include "applications/dynamicProps/TimeCorrFunc.hpp"
#include "brains/SimInfo.hpp"
#include "primitives/StuntDouble.hpp"
#include "selection/SelectionManager.hpp"
#include "utils/StaticAccumulator.hpp"

namespace OpenMD {

  //! Frequency Fluctuation Correlation Function
  /*! See <http://dx.doi.org/10.1021/jp010798o>

    Williams, Loring, and Fayer, J. Phys. Chem. B 2001, 105, 4068-4071 for details

    Makes an assumption that the frequency depends on the
    orientation of the dipole relative to the local Electric
    Field.  Requires a preCorrelate to compute the mean of
    \f$ = \langle \mathbf{E} \cdot \mathbf{u} \rangle \f$
  */

  class FreqFlucCorrFunc : public ObjectACF<RealType> {
  public:
    FreqFlucCorrFunc(SimInfo* info, const std::string& filename,
                     const std::string& sele1, const std::string& sele2);

  private:
    virtual int computeProperty1(int frame, StuntDouble* sd);
    virtual RealType calcCorrVal(int frame1, int frame2, int id1, int id2);
    virtual void validateSelection(SelectionManager& seleMan);

    std::vector< std::vector<RealType> > ue_;
    RealAccumulator ueStats_ {};
  };
}

#endif
