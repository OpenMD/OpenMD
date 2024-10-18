/*
 * Copyright (c) 2004-present, The University of Notre Dame. All rights
 * reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 * 1. Redistributions of source code must retain the above copyright notice,
 *    this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright notice,
 *    this list of conditions and the following disclaimer in the documentation
 *    and/or other materials provided with the distribution.
 *
 * 3. Neither the name of the copyright holder nor the names of its
 *    contributors may be used to endorse or promote products derived from
 *    this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 *
 * SUPPORT OPEN SCIENCE!  If you use OpenMD or its source code in your
 * research, please cite the following paper when you publish your work:
 *
 * [1] Drisko et al., J. Open Source Softw. 9, 7004 (2024).
 *
 * Good starting points for code and simulation methodology are:
 *
 * [2] Meineke, et al., J. Comp. Chem. 26, 252-271 (2005).
 * [3] Fennell & Gezelter, J. Chem. Phys. 124, 234104 (2006).
 * [4] Sun, Lin & Gezelter, J. Chem. Phys. 128, 234107 (2008).
 * [5] Vardeman, Stocker & Gezelter, J. Chem. Theory Comput. 7, 834 (2011).
 * [6] Kuang & Gezelter, Mol. Phys., 110, 691-701 (2012).
 * [7] Lamichhane, Gezelter & Newman, J. Chem. Phys. 141, 134109 (2014).
 * [8] Bhattarai, Newman & Gezelter, Phys. Rev. B 99, 094106 (2019).
 * [9] Drisko & Gezelter, J. Chem. Theory Comput. 20, 4986-4997 (2024).
 */

/**
 * @file next_combination.hpp
 * @author    tlin
 * @date  10/27/2004
 * @version 1.0
 */

#ifndef UTILS_NEXT_COMBINATION_HPP
#define UTILS_NEXT_COMBINATION_HPP

#include <iostream>
#include <iterator>
#include <vector>

namespace OpenMD {

  /**
   * @brief STL next_permuation-like combination sequence generator.
   * Given the first and last iterator of a sequence, next_combination
   * iteratively generates all possible combinations.
   * @return if more combination is availiable, otherwise return false
   * @param iterContainer iterator container
   * @param first the first iterator
   * @param last the last iterator
   * @note first and last must be random access iterators and iterContainer must
   * be the container of random access iterators . And all of the iteratos in
   * iterContainer must be within the range [first, last)
   *
   * @code
   * std::vector<int> iv;
   * iv.push_back(1);
   * iv.push_back(8);
   * std::vector<std::vector<int>::iterator> ic;
   * while(next_combination(ic, iv.begin(), iv.end())) {
   *     for (i =  ic.begin(); i < ic.end(); ++i) {
   *         std::cout << **i << "\t";
   *     }
   *     std::cout << std::endl;
   * }
   * //output
   * //1
   * //8
   * //1  8
   * @endcode
   */
  template<class RandomAccessIterator>
  bool next_combination(std::vector<RandomAccessIterator>& iterContainer,
                        RandomAccessIterator first, RandomAccessIterator last) {
    if (first == last) { return false; }

    RandomAccessIterator endIter = --last;
    typename std::vector<RandomAccessIterator>::iterator i =
        iterContainer.end();

    if (iterContainer.empty()) {
      // if sequence is empty, we insert the first iterator
      iterContainer.insert(iterContainer.end(), first);
      return true;
    } else if (*(--i) != endIter) {
      // if the last iterator in iterContainer does not reaches the end, just
      // increase its iterator by 1
      ++(*i);
      return true;
    } else {  // the last iterator in iterContainer does not reaches the end

      // starts at the end of the sequence and works its way towards the front,
      // looking for two consecutive members of the sequence where the
      // difference between them is greater than one. For example , if the
      // sequence contains 1, 5, 8, 9 (total number is 10, first is 0 and the
      // last is 10 (due to STL's half open range)). At the end of while loop, j
      // will point to 5, and i will point to 8, next combination should be 1,
      // 6, 7, 8. If j is less than zero, it means it already reaches the last
      // combination of current size. For instance, sequence may contain 6, 7,
      // 8, 9 at this time, we need to increase the size of combination to 5
      typename std::vector<RandomAccessIterator>::iterator j = i;
      --j;
      while (j >= iterContainer.begin() && *i == *j + 1) {
        --i;
        --j;
      };

      RandomAccessIterator raIter;
      if (j - iterContainer.begin() <
          0) {  // reaches the last combination of current size
        // half open range
        if (last - first + 1 == int(iterContainer.size())) {
          // if the current size equals to total number, done
          return false;
        } else {
          // push the first currentSize+1 element into sequence

          for (i = iterContainer.begin(), raIter = first;
               i != iterContainer.end(); ++i, ++raIter) {
            *i = raIter;
          }
          iterContainer.insert(iterContainer.end(), raIter);

          return true;
        }

      } else {
        ++(*j);
        raIter = *j;
        for (; i != iterContainer.end(); ++i) {
          ++raIter;
          *i = raIter;
        }
        return true;
      }
    }
  }  // end next_combination

}  // namespace OpenMD

#endif  // UTILS_NEXT_COMBINATION_HPP
