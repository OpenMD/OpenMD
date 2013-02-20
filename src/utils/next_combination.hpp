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
 */
 
/**
 * @file next_combination.hpp
 * @author    tlin
 * @date  10/27/2004
 * @version 1.0
 */ 

#ifndef UTILS_NEXT_COMBINATION_HPP
#define UTILS_NEXT_COMBINATION_HPP

#include <vector>
#include <iterator>
#include <iostream>
namespace OpenMD {

  /**
   * @brief STL next_permuationtation like combination sequence generator.
   * Given the first and last iterator of a sequence, next_combination iteratively generates all 
   * possible combinations.
   * @return if more combination is availiable, otherwise return false
   * @param iterContainer iterator container
   * @param first the first iterator
   * @param last the last iterator
   * @note first and last must be random access iterators and iterContainer must be the container of  
   * random access iterators . And all of the iteratos in iterContainer must be within the range [first, last)
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
  bool next_combination(std::vector<RandomAccessIterator>& iterContainer, RandomAccessIterator first, RandomAccessIterator last) {
    if (first == last) {
      return false;
    }
    
    RandomAccessIterator endIter = --last;
    typename std::vector<RandomAccessIterator>::iterator i = iterContainer.end();
    
    if (iterContainer.empty()) {
      //if sequence is empty, we insert the first iterator
      iterContainer.insert(iterContainer.end(), first);
      return true;
    } else if (*(--i) != endIter){
      //if the last iterator in iterContainer does not reaches the end, just increase its iterator by 1 
      ++(*i);
      return true;
    } else {// the last iterator in iterContainer does not reaches the end

      //starts at the end of the sequence and works its way towards the front, looking for two 
      //consecutive members of the sequence where the difference between them is greater 
      //than one. For example , if the sequence contains 1, 5, 8, 9 (total number is 10, first is 0
      //and the last is 10 (due to STL's half open range)). At the end of while
      //loop, j will point to 5, and i will point to 8, next combination should be 1, 6, 7, 8.
      //If j is less than zero, it means it already reaches the last combination of current size.
      //For instance, sequence may contain 6, 7, 8, 9 at this time, we need to increase the size
      // of combination to 5
      typename std::vector<RandomAccessIterator>::iterator j = i;
      j--;
      while( j >= iterContainer.begin() && *i == *j + 1){
	i--;
	j--;
      };

      RandomAccessIterator raIter;
      if (j - iterContainer.begin() < 0) { //reaches the last combination of current size
	//half open range
	if (last  - first + 1  == iterContainer.size()) {
	  //if the current size equals to total number, done
	  return false;
	} else {

	  //push the first currentSize+1 element into sequence 
                
	  for(i = iterContainer.begin(), raIter= first; i  != iterContainer.end(); ++i, ++raIter) {
	    *i = raIter;
	  }    
	  iterContainer.insert(iterContainer.end(), raIter);
                
	  return true; 
	}
            
      } else {
	++(*j);
	raIter = *j;
	for(; i != iterContainer.end(); ++i) {
	  ++raIter;
	  *i = raIter;
	}
	return true;
      }
    }
  } //end next_combination

} //end namespace OpenMD
#endif //UTILS_NEXT_COMBINATION_HPP

