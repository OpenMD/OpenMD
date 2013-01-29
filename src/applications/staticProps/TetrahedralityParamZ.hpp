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
 * [3]  Sun, Lin & Gezelter, J. Chem. Phys. 128, 24107 (2008).          
 * [4]  Kuang & Gezelter,  J. Chem. Phys. 133, 164101 (2010).
 * [5]  Vardeman, Stocker & Gezelter, J. Chem. Theory Comput. 7, 834 (2011).
 * [6]  Kuang & Gezelter, Mol. Phys., 110, 691-701 (2012).
 */

#ifndef APPLICATIONS_STATICPROPS_TETRAHEDRALITYPARAMZ_HPP
#define APPLICATIONS_STATICPROPS_TETRAHEDRALITYPARAMZ_HPP
#include "selection/SelectionEvaluator.hpp"
#include "selection/SelectionManager.hpp"
#include "applications/staticProps/StaticAnalyser.hpp"
#include "math/Vector3.hpp"

namespace OpenMD {

  /**
   * @class TetrahedralityParamZ
   * @brief Tetrahedrality ParameterZ
   *
   * Computes local tetrahedral order parameter Q as introduced in:
   *
   *  "A new order parameter for tetrahedral configurations," by P.-L. Chau and 
   *    A.J. Hardwick, Mol. Phys. 93, pp. 511-518 (1998).
   *
   *
   * Note that we use a rescaled version of the tetrahedral order
   * parameter 'Q' such that a perfectly tetrahedral configuration has a Q value
   * of 1 and an ideal gas configuration has a Q value of 0. This rescaled
   * version of the tetrahedrality parameter was first introduced in:
   *
   *  "Relationship between structural order and the anomalies of liquid water,"
   *    by J.R. Errington and P.G. Debenedetti, Nature 409, pp. 318-321 (2001).
   * 
   *
   * Characterization of the spatial correlations of the the local
   * order parameter Q are done according to the procedure outlined
   * in:
   *
   *   "Space-time correlations in the orientational order parameter and the
   *    orientational entropy of water," by P. Kumar, S.V. Buldyrev, and 
   *    H.E. Stanley, arXiv:0807.4699v1 [cond-mat.soft] 29 Jul 2008.
   *
   */
  class TetrahedralityParamZ : public StaticAnalyser{
  public:
    TetrahedralityParamZ(SimInfo* info, const std::string& filename, 
                         const std::string& sele, double rCut, int nzbins);
    int getNZBins(){
      return nZBins_;
    }

    virtual ~TetrahedralityParamZ();
    virtual void process();
    
  private:
    void writeQz();

    Snapshot* currentSnapshot_;
    std::string selectionScript_;
    SelectionManager seleMan_;
    SelectionEvaluator evaluator_;
    RealType rCut_;
    int nZBins_;
    std::vector<RealType> zBox_;
    std::vector<RealType> sliceQ_;
    std::vector<int> sliceCount_;
  };
}
#endif

