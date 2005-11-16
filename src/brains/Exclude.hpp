/*
 * Copyright (c) 2005 The University of Notre Dame. All Rights Reserved.
 *
 * The University of Notre Dame grants you ("Licensee") a
 * non-exclusive, royalty free, license to use, modify and
 * redistribute this software in source and binary code form, provided
 * that the following conditions are met:
 *
 * 1. Acknowledgement of the program authors must be made in any
 *    publication of scientific results based in part on use of the
 *    program.  An acceptable form of acknowledgement is citation of
 *    the article in which the program was described (Matthew
 *    A. Meineke, Charles F. Vardeman II, Teng Lin, Christopher
 *    J. Fennell and J. Daniel Gezelter, "OOPSE: An Object-Oriented
 *    Parallel Simulation Engine for Molecular Dynamics,"
 *    J. Comput. Chem. 26, pp. 252-271 (2005))
 *
 * 2. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 *
 * 3. Redistributions in binary form must reproduce the above copyright
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
 */
 
#ifndef BRAINS_EXCLUDE_HPP

#define BRAINS_EXCLUDE_HPP

#include <iostream>
#include <set>
#include <utility>
#include <vector>

namespace oopse {

  /**
   * @class Exclude Exclude.hpp "brains/Exclude.hpp"
   * @brief Exclude class maintains the exclusion list of the simulation by the global indices of the 
   * atoms. The exclusion list will be passed to fortran in a plain array.
   */
  class Exclude {
  public:

    Exclude() : modified_(false) {}


    /** Adds a pair into this Exclude class */
    void addPair(int i, int j);

    void addPairs(std::set<int>& set1, std::set<int>& set2);
    template<typename IterType1, typename IterType2>
    void addPairs(IterType1 iter1_first, IterType1 iter1_last, IterType2 iter2_first, IterType2 iter2_last);

    /** Remove a pair from Exclude class */
    void removePair(int i, int j);

    void removePairs(std::set<int>& set1, std::set<int>& set2);
    template<typename IterType1, typename IterType2>
    void removePairs(IterType1 iter1_first, IterType1 iter1_last, IterType2 iter2_first, IterType2 iter2_last);

    /** Checks whether pair (i, j) is in this Exclude class */
    bool hasPair(int i, int j);

    /** Returns the number of exclusion pair */
    int getSize();

    /** Returns the size of exclusion list */
    int getSizeOfExcludeList();

    /** Returns the exclusion pairs in a plain array*/
    int *getExcludeList();

    /** write out the exclusion list to an ostream */
    friend std::ostream& operator <<(std::ostream& o, Exclude& e);

  private:

    std::set < std::pair<int, int> > excludeSet_;
    std::vector <int> excludeList_;
    bool modified_;
  };

}      //end namespace oopse

#endif // __EXCLUDE_H__
