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

#ifndef UTILS_GENERATECOMBINATION_HPP
#define UTILS_GENERATECOMBINATION_HPP

#include <vector>

namespace oopse {
//this algorithm can be modified to be compatible with STL container
//Basically, changes totNum to last iteratir and beginIndex to first iterator, sequence from
//a vector<int> to a vector<iterator>, STL compatible version will be committed next time
bool next_combination(std::vector<int> sequence, int totNum, int beginIndex = 0) {
    int endIndex = beginIndex + totNum - 1;
    int currentSize;
    currentSize = sequence.size();
    
    if (currentSize == 0) {
        //if sequence is empty, we push the begin index into it
        sequence.push_back(beginIndex);
        return true;
    } else if (sequence[currentSize - 1] !=  endIndex) {
        //if the last element of the sequence does not reach the end index, just increment it 
        sequence[currentSize - 1] ++;
        return true;
    } else if (totNum == 1) {
        //if total number is one, and the last element already reaches the end index
        return false;
    } else {// At this point, the sequence contains at least two elements, and the last one already reaches
                // end index
        int i;
        int j;

        //starts at the end of the sequence and works its way towards the front, looking for two 
        //consecutive members of the sequence where the difference between them is greater 
        //than one. For example , if the sequence contains 1, 5, 8, 9 (total number is 10, begin
        //index is 0, therefore 9 is the end index, and the current size is 4). At the end of while
        //loop, j will point to 5, and i will point to 8, next combination should be 1, 6, 7, 8.
        //If j is less than zero, it means it already reaches the last combination of current size.
        //For instance, sequence may contain 6, 7, 8, 9 at this time, we need to increase the size
        // of combination to 5
        i = currentSize;
        j = i - 1;
        while( j >= 0 && sequence[i] == sequence[j] + 1){
            i--;
            j = i -1;
        };

        if (j < 0) { //reaches the last combination of current size

            if (currentSize == totNum) {
                //if the current size equals to total number, done
                return false;
            } else {

                //push the first currentSize+1 element into sequence 
                for(int i = 0; i < currentSize; i++)
                    sequence[i] = beginIndex + i;
                sequence.push_back(beginIndex + currentSize);
                
		return true; 
	    }            
        } else {
            ++sequence[j];

            for(int k = j + 1; k < currentSize; k++)
                sequence[k] = sequence[k-1] + 1; 

            return true;
        }
        
        
    }
    
}

} //end namespace oopse
#endif //UTILS_GENERATECOMBINATION_HPP
