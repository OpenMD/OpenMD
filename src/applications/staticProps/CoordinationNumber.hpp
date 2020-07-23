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

#ifndef APPLICATIONS_STATICPROPS_COORDINATIONNUMBER_HPP
#define APPLICATIONS_STATICPROPS_COORDINATIONNUMBER_HPP

#include <string>
#include <vector>
#include "selection/SelectionEvaluator.hpp"
#include "selection/SelectionManager.hpp"
#include "applications/staticProps/StaticAnalyser.hpp"

using namespace std;
namespace OpenMD {

  /**
   * @class CoordinationNumber
   * @brief CoordinationNumber
   *
   * Computes a distribution of coordination numbers defined as the number of
   * atoms in \a sele2 that are within \a rCut of the atom in \a sele1
   *
   * Note that extra parameters must be declared:
   *
   *   \param rCut cutoff radius for finding lists of nearest neighbors
   *   \param sele1 selection of StuntDoubles used for the distribution
   *   \param sele2 selection of StuntDoubles used for nearest neighbor computation
   */
  class CoordinationNumber : public StaticAnalyser {
    
  public:
    CoordinationNumber(SimInfo* info, const std::string& filename,
                       const std::string& sele1, const std::string& sele2,
                       RealType rCut, int bins);

    virtual ~CoordinationNumber();
    virtual void process();
    virtual void writeOutput();

  protected:
    virtual RealType computeCoordination(int a, vector<vector<int> > neighbors);
    RealType rCut_;
    int bins_;
    
    std::string sele1_;
    SelectionManager seleMan1_;
    SelectionEvaluator evaluator1_;

    std::string sele2_;
    SelectionManager seleMan2_;
    SelectionEvaluator evaluator2_;

    int selectionCount1_;
    int selectionCount2_;
    int nnMax_;
    RealType delta_;
    int count_;
    std::vector<RealType>  histogram_;
  };

  /**
   * @class SCN
   * @brief Secondary Coordinate Number
   *
   * Computes a distribution of secondary coordination numbers where
   * each atom is assigned the mean coordination number of the
   * neighboring atoms.
   */
  class SCN : public CoordinationNumber {
    
  public:
    SCN(SimInfo* info, const std::string& filename, const std::string& sele1,
        const std::string& sele2, RealType rCut, int bins);

    virtual ~SCN();
    virtual RealType computeCoordination(int a, vector<vector<int> > neighbors);
  };
    
  /**
   * @class GCN
   * @brief Generalized Coordinate Number
   *
   * Computes a distribution of generalized coordinate numbers as
   * described in:
   *
   *   "Finding optimal surface sites on heterogeneous catalysts by
   *    counting nearest neighbors," by F. Calle-Vallejo et al.,
   *    Science 350(6257) pp. 185-189 (2015).
   *    http://dx.doi.org/10.1126/science.aab3501
   */
  class GCN : public CoordinationNumber {
    
  public:
    GCN(SimInfo* info, const std::string& filename, const std::string& sele1,
        const std::string& sele2, RealType rCut, int bins);

    virtual ~GCN();
    virtual RealType computeCoordination(int a, vector<vector<int> > neighbors);
  };

}
#endif
