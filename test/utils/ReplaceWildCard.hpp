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
 * @file ReplaceWildCard.hpp
 * @author    tlin
 * @date  10/27/2004
 * @version 1.0
 */ 

#ifndef UTILS_REPLACEWILDCARD_HPP
#define UTILS_REPLACEWILDCARD_HPP

#include <iostream>
#include <vector>

namespace oopse{

//use -1 to represent the wild card, it is easy and cheap to handle integer
const int WildCard = -1;

std::vector<std::vector<int> > ReplaceWildCard(int beginIndex, int endIndex, int nWildCard);
std::vector<std::vector<int> > join( const std::vector<int>& firstPart, const std::vector<std::vector<int> >& secondPart);

/**
 * Driver function for replacing
 * @code
 * std::vector<std::vector<int> > result;
 * result = ReplaceAll(3);
 * /contents in result
 * //0, 1, 2
 * //-1, 1, 2
 * //0, -1, 2
 * //0, 1, -1
 * //-1, -1, 2
 * //-1, 1, -1 
 * //0, -1, -1
 * //-1, -1, -1
 * @endcode
 */
std::vector<std::vector<int> > ReplaceAll(int num) {
    std::vector<std::vector<int> > results;
    std::vector<std::vector<int> > v;
    
    for (int i = 0; i <= num; i++) {
        v = ReplaceWildCard(0, num -1, i);
        results.insert(results.end(), v.begin(), v.end());
    }
    return results;
}

/** Replace a sequence with n wildcards, returns all of the possible replacement*/
std::vector<std::vector<int> > ReplaceWildCard(int beginIndex, int endIndex, int nWildCard) {
    std::vector<std::vector<int> > results;

    int len = endIndex + 1 - beginIndex;
    int nRecursive = len - nWildCard;

    if (nWildCard == 0) {
        //if the number of the wild card is zero, just return the whole sequence
        std::vector<int> singleResult;

        for(int i = beginIndex; i <= endIndex; i++)
            singleResult.push_back(i);

        results.push_back(singleResult);
        return results;
    }
        
    if (len < nWildCard) {
        //give warning message
        std::cerr << "Error" << std::endl;
        exit(1);
    } else if (len == nWildCard) {
        //if the lengths are the same, we only have one choice
        //replace all of the indices with WildCard
        std::vector<int> singleResult(nWildCard, WildCard);
        results.push_back(singleResult);
        return results;
    } else {
        //we need to recursively calling ReplaceWildCard
        std::vector<int> firstPart;
        std::vector<std::vector<int> > secondPart;
        std::vector<std::vector<int> > sequences;

        for (int i = 0; i <=nRecursive; i ++) {
            
            firstPart.clear();
            for(int j = 0; j < i; ++j) {
                firstPart.push_back(beginIndex + j);
            }       
            firstPart.push_back(WildCard);
            
            secondPart = ReplaceWildCard(beginIndex + i + 1, endIndex, nWildCard - 1); 
            sequences = join(firstPart, secondPart);
            results.insert(results.end(), sequences.begin(), sequences.end());              
        }

        return results;
    }
}

std::vector<std::vector<int> > join( const std::vector<int>& firstPart, const std::vector<std::vector<int> >& secondPart){
    std::vector<std::vector<int> > results(secondPart.size());

    for (int i = 0; i < secondPart.size(); i++) {
        results[i].insert(results[i].end(), firstPart.begin(), firstPart.end());
        results[i].insert(results[i].end(), secondPart[i].begin(), secondPart[i].end());
    }

    return results;
}


}//end namespace std

#endif //UTILS_REPLACEWILDCARD_HPP
