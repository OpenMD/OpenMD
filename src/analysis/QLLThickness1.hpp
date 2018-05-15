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
 * [6]  Kuang & Gezelter, Mol. Phys., 110, 691-701 (2012).
 */

#ifndef ANALYSIS_QLLTHICKNESS1_HPP
#define ANALYSIS_QLLTHICKNESS1_HPP
#include "selection/SelectionEvaluator.hpp"
#include "selection/SelectionManager.hpp"
#include "analysis/SequentialAnalyzer.hpp"
#include "math/Vector3.hpp"

namespace OpenMD {

  /**
   * @class QLLThickness1
   * @brief QLLThickness1
   *
   * Computes the thickness of the QLL as determined by the local tetrahedral 
   * order parameter Q as introduced in:
   *
   * "The Thickness of a Liquid Layer on the Free Surface of Ice as Obtained
   * From Computer Simulation" by M. Conde, C. Vega, and A. Patrykiejew
   * J. Chem. Phys, 2008, 129 (014702).
   *
   * The threshold q value (q_t) has the default of that for the 
   * TIP4P/Ice model, q_t = 0.9076 as determined by Conde2008.
   *
   */
  class QLLThickness1 : public SequentialAnalyzer{
  public:
    QLLThickness1(SimInfo* info, const std::string& filename, 
                         const std::string& sele1, const std::string& sele2, 
                         double rCut, RealType qt=0.9076);
    
    virtual ~QLLThickness1();
    virtual void doFrame(int frame);
    
  private:
    RealType rCut_;
    RealType qt_;
    unsigned int nLiquid_;
  };
}
#endif

