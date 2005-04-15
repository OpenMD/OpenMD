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
 
/**
 * @file LocalIndexManager.hpp
 * @author    tlin
 * @date  11/02/2004
 * @version 1.0
 */ 

#ifndef UTILS_LOCALINDEXMANAGER_HPP
#define UTILS_LOCALINDEXMANAGER_HPP
#include <algorithm>
#include <iostream>
#include <cassert>
#include <list>
#include <utility>

namespace oopse {

  /**
   * @class IndexListContainer
   * @brief 
   * @todo documentation
   */
  class IndexListContainer{
  public:
    static const int MAX_INTEGER = 2147483647;
    typedef std::list<std::pair<int, int> >::iterator  IndexListContainerIterator;

        
    IndexListContainer(int minIndex = 0, int maxIndex = MAX_INTEGER) 
      :minIndex_(minIndex), maxIndex_(maxIndex) {
            
	indexContainer_.push_back(std::make_pair(minIndex, maxIndex));
      }
        
    int pop() {
            
      if (indexContainer_.empty()) {
	std::cerr << "" << std::endl;
      }

      IndexListContainerIterator i = indexContainer_.begin();
      int result;
            
      result = indexContainer_.front().first;
            
      if (indexContainer_.front().first  == indexContainer_.front().second) {
	indexContainer_.pop_front();
      } else if (indexContainer_.front().first < indexContainer_.front().second) {
	++indexContainer_.front().first;
      } else {
	std::cerr << "" << std::endl;
      }

      return result;
    }


    /**
     *
     */
    void insert(int index) {
      IndexListContainerIterator insertPos = internalInsert(index, index);            
      merge(insertPos);
    }

    /**
     * Reclaims an index range
     * @param beginIndex
     * @param endIndex
     */
    void insert(int beginIndex, int endIndex) {
      IndexListContainerIterator insertPos = internalInsert(beginIndex, endIndex); 
      merge(insertPos);            
    }

    /**
     * Reclaims an index array.
     * @param indices
     */
    void insert(std::vector<int>& indices){

      if (indices.empty()) {
	return;
      }
            
      std::sort(indices.begin(), indices.end());
      std::unique(indices.begin(), indices.end());

      std::vector<int>::iterator i;
      IndexListContainerIterator insertPos;
      int beginIndex;

      beginIndex = indices[0];
        
      for ( i = indices.begin() + 1 ; i != indices.end(); ++i) {
	if (*i != *(i -1) + 1) {
	  insert(beginIndex, *(i-1));
	  beginIndex = *i;
	}
      }
            
            
    }
        
    std::vector<int> getIndicesBefore(int index) {
      std::vector<int> result;
      IndexListContainerIterator i;

      for(i = indexContainer_.begin(); i != indexContainer_.end(); ++i) {
	if ((*i).first > index) {
	  //we locate the node whose minimum index is greater that index
	  indexContainer_.erase(indexContainer_.begin(), i);
	  break;
	} else if ((*i).second < index) {
	  //The biggest index current node hold is less than the index we want
	  for (int j = (*i).first; j <= (*i).second; ++j) {
	    result.push_back(j);
	  }
	  continue;
	} else if ((*i).first == (*i).second) {
	  //the index happen to equal to a node which only contains one index
	  result.push_back((*i).first);
	  indexContainer_.erase(indexContainer_.begin(), i);
	  break;
	} else {

	  for (int j = (*i).first; j < index; ++j) {
	    result.push_back(j);
	  }
                        
	  (*i).first =  index;                   
	  indexContainer_.erase(indexContainer_.begin(), i);
	  break;
	}
                
      }

      return result;

    }
        
    int getMaxIndex() {
      return maxIndex_;
    }
                
  private:
        
    IndexListContainerIterator internalInsert(int beginIndex, int endIndex) {

      if (beginIndex > endIndex) {
	std::swap(beginIndex, endIndex);
	std::cerr << "" << std::endl;
      }

      if (endIndex > maxIndex_) {
	std::cerr << "" << std::endl;
      }

                      
      IndexListContainerIterator j;

      IndexListContainerIterator i = indexContainer_.begin();  
      for (; i != indexContainer_.end(); ++i) {
	if ((*i).first > endIndex) {
	  indexContainer_.insert(i, std::make_pair(beginIndex, endIndex));
	  return --i;
	} else if ((*i).second < beginIndex) {
	  continue;
	} else {
	  std::cerr << "" << std::endl;
	}
      }

      indexContainer_.push_back(std::make_pair(beginIndex, endIndex));
      return --indexContainer_.end();

    }

    void merge(IndexListContainerIterator i) {
      IndexListContainerIterator j;

      //check whether current node can be merged with its previous node            
      if ( i != indexContainer_.begin()) {
	j = i;
	--j;
	if (j != indexContainer_.begin() && (*j).second + 1 == (*i).first) {
	  (*i).first = (*j).first;
	  indexContainer_.erase(j);
	}
      }

      //check whether current node can be merged with its next node            
      if ( i != indexContainer_.end()) {
	j = i;
	++j;

	if (j != indexContainer_.end() && (*i).second + 1 == (*j).first) {
	  (*i).second = (*j).second;
	  indexContainer_.erase(j);
	}
      }

    }
    int minIndex_;
    int maxIndex_;        
    std::list<std::pair<int, int> > indexContainer_;
  };


  /**
   * @class LocalIndexManager LocalIndexManager.hpp "utils/LocalIndexManager.hpp"
   * @brief
   */
  class LocalIndexManager {
  public:
        
    int getNextAtomIndex() {
      return atomIndexContainer_.pop();
    }

    std::vector<int> getAtomIndicesBefore(int index) {
      return atomIndexContainer_.getIndicesBefore(index);
    }
        
    void releaseAtomIndex(int index) {
      atomIndexContainer_.insert(index);
    }

    void releaseAtomIndex(int beginIndex, int endIndex) {
      atomIndexContainer_.insert(beginIndex, endIndex);
    }

    void releaseAtomIndex(std::vector<int> indices) {
      atomIndexContainer_.insert(indices);
    }

    int getNextRigidBodyIndex() {
      return rigidBodyIndexContainer_.pop();
    }

    std::vector<int> getRigidBodyIndicesBefore(int index) {
      return rigidBodyIndexContainer_.getIndicesBefore(index);
    }
        
    void releaseRigidBodyIndex(int index) {
      rigidBodyIndexContainer_.insert(index);
    }

    void releaseRigidBodyIndex(int beginIndex, int endIndex) {
      rigidBodyIndexContainer_.insert(beginIndex, endIndex);
    }

    void releaseRigidBodyIndex(std::vector<int> indices) {
      rigidBodyIndexContainer_.insert(indices);
    }
 
  private:

    IndexListContainer atomIndexContainer_;
    IndexListContainer rigidBodyIndexContainer_;
  };

} //end namespace oopse
#endif //UTILS_LOCALINDEXMANAGER_HPP
