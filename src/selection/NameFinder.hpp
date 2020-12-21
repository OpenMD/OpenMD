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

#ifndef SELECTION_NAMEFINDER_HPP
#define SELECTION_NAMEFINDER_HPP

#include <map>
#include <memory>
#include <set>
#include <string>

#include "brains/SimInfo.hpp"
#include "selection/SelectionSet.hpp"

namespace OpenMD {

  class TreeNode{
  public:
    std::string name;
    SelectionSet bs;
    std::map<std::string, std::shared_ptr<TreeNode>> children;
  };

  using TreeNodePtr = std::shared_ptr<TreeNode>;

  class NameFinder {
  public:
    NameFinder(SimInfo* info);
    SelectionSet match(const std::string &name);

  private:
    void loadNames();
    void matchMolecule(const std::string &molName, SelectionSet &bs);
    void matchStuntDouble(const std::string &molName, const std::string &sdName, SelectionSet &bs);
    void matchRigidAtoms(const std::string &molName, const std::string &rbName, const std::string &rbAtomName, SelectionSet &bs);
    void matchBond(const std::string &molName, const std::string &bondName,  SelectionSet &bs);
    void matchBend(const std::string &molName, const std::string &bendName, SelectionSet &bs);
    void matchTorsion(const std::string &molName, const std::string &torsionName, SelectionSet &bs);
    void matchInversion(const std::string &molName, const std::string &inversionName, SelectionSet &bs);

    void matchInternalIndex(const std::string &name, int internalIndex, SelectionSet &bs);

    TreeNodePtr createNode(TreeNodePtr parent, const std::string &name);
    std::vector<TreeNodePtr> getAllChildren(TreeNodePtr node);
    std::vector<TreeNodePtr> getMatchedChildren(TreeNodePtr node, const std::string &name);
    bool isMatched(const std::string &str, const std::string &wildcard);

    bool isInteger(const std::string &str);

    SimInfo* info_ {nullptr};
    vector<int> nObjects_;
    TreeNodePtr root_ {nullptr};
  };
}
#endif
