/*
 * Copyright (C) 2000-2004  Object Oriented Parallel Simulation Engine (OOPSE) project
 * 
 * Contact: oopse@oopse.org
 * 
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public License
 * as published by the Free Software Foundation; either version 2.1
 * of the License, or (at your option) any later version.
 * All we ask is that proper credit is given for our work, which includes
 * - but is not limited to - adding the above copyright notice to the beginning
 * of your source code files, and to any copyright notice that you may distribute
 * with programs based on this work.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 * 
 * You should have received a copy of the GNU Lesser General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
 *
 */

/**
 * @file GenerateCombination.hpp
 * @author    tlin
 * @date  10/27/2004
 * @version 1.0
 */ 

#ifndef UTILS_NEXT_COMBINATION_HPP
#define UTILS_NEXT_COMBINATION_HPP

#include <vector>
#include <iterator>
#include <iostream>
namespace oopse {

/**
 * @fn bool next_combination(IteratorContainer<RandomAccessIterator>& iterContainer, RandomAccessIterator first, RandomAccessIterator last) 
 * @brief STL next_permuationtation like combination sequence generator.
 * Given the first and last iterator of a sequence, next_combination iteratively generates all possible combination.
 * @param iterContainer iterator container
 * @param first the first iterator
 * @param last the last iterator
 * @note first and last must be random access iterators and iterContainer must be the container which
 * element is iterator. And all of the iteratos in iterContainer must be within the range [first, last)
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
 */
template<class RandomAccessIterator, template<typename ELEM, typename = std::allocator<ELEM> > class IteratorContainer>
bool next_combination(IteratorContainer<RandomAccessIterator>& iterContainer, RandomAccessIterator first, RandomAccessIterator last) {
    if (first == last) {
        return false;
    }
    
    RandomAccessIterator endIter = --last;
    typename IteratorContainer<RandomAccessIterator>::iterator i = iterContainer.end();
    
    if (iterContainer.empty()) {
        //if sequence is empty, we insert the first iterator
        iterContainer.insert(iterContainer.end(), first);
        return true;
    } else if (*(--i) != endIter){
        //if the last iterator in iterContainer does not reaches the end, just increment it 
        ++(*i);
        return true;
    } else {// the last iterator in iterContainer does not reaches the end

        //starts at the end of the sequence and works its way towards the front, looking for two 
        //consecutive members of the sequence where the difference between them is greater 
        //than one. For example , if the sequence contains 1, 5, 8, 9 (total number is 10, begin
        //index is 0, therefore 9 is the end index, and the current size is 4). At the end of while
        //loop, j will point to 5, and i will point to 8, next combination should be 1, 6, 7, 8.
        //If j is less than zero, it means it already reaches the last combination of current size.
        //For instance, sequence may contain 6, 7, 8, 9 at this time, we need to increase the size
        // of combination to 5
        typename IteratorContainer<RandomAccessIterator>::iterator j = i;
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

bool replaceWildCard(std::vector<std::vector<std::string>::iterator>& cont,
                                             std::vector<std::string>& sequence, std::vector<std::string>& result,
                                             const std::string& wildCard = "X") {
    if (cont.size() > sequence.size()) {
        std::cerr << "the size of iterator container is greater than the size of sequence";
    }

    bool hasMoreCombination = next_combination(cont, sequence.begin(), sequence.end());
    if (hasMoreCombination) {
        result.clear();
        result.insert(result.begin(), sequence.size(), wildCard);
        std::vector<std::vector<std::string>::iterator>::iterator i;
        for ( i = cont.begin(); i != cont.end(); i++){
            result[*i - sequence.begin()] = **i;
        }
    }

     return hasMoreCombination;
    
}

} //end namespace oopse
#endif //UTILS_NEXT_COMBINATION_HPP

