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

#include <algorithm>
#include <cassert>
#include <string>

#include "utils/BitSet.hpp"
#include "utils/Algorithm.hpp"

namespace oopse {
int BitSet::countBits() {
#ifdef __RWSTD    
    //For the compiler(Sun, MSVC6.0) binding with RougeWave STL Library, we need to use old-style
    // std::count which is error-prone.
    int count = 0;
    std::count(bitset_.begin(), bitset_.end(), true, count);
    return count;
#else
    return std::count(bitset_.begin(), bitset_.end(), true);
#endif
}

void BitSet::flip(int fromIndex, int toIndex) {
    assert(fromIndex <= toIndex);
    assert(fromIndex >=0);
    assert(toIndex <= size());
    std::vector<char>::iterator first = bitset_.begin() + fromIndex;
    std::vector<char>::iterator last = bitset_.begin() + toIndex;

    std::transform(first, last, first, std::logical_not<bool>());
        
}

BitSet BitSet::get(int fromIndex, int toIndex) {
    assert(fromIndex <= toIndex);
    assert(fromIndex >=0);
    assert(toIndex <= size());
    std::vector<char>::iterator first = bitset_.begin() + fromIndex;
    std::vector<char>::iterator last = bitset_.begin() + toIndex;

    BitSet result;
    std::copy(first, last, std::back_inserter(result.bitset_));
    return result;
}

bool BitSet::none() {
    std::vector<char>::iterator i = std::find(bitset_.begin(), bitset_.end(), true);
    return i == bitset_.end() ? true : false;
}
    
int BitSet::nextOffBit(int fromIndex) const {
    if (fromIndex <= -1) {
        //in case -1 or other negative number is passed to this function
        return -1;
    }
    
    ++fromIndex;
    while (fromIndex < size()) {
        if (!bitset_[fromIndex]) {
            return fromIndex;
        }
        ++fromIndex;
    }

    return -1;
}

int BitSet::nextOnBit(int fromIndex) const {
    if (fromIndex <= -1) {
        //in case -1 or other negative number is passed to this function
        return -1;
    }

    ++fromIndex;
    while (fromIndex < size()) {
        if (bitset_[fromIndex]) {
            return fromIndex;
        }
        ++fromIndex;
    }

    return -1;
}

void BitSet::andOperator (const BitSet& bs) {
    assert(size() == bs.size());

    std::transform(bs.bitset_.begin(), bs.bitset_.end(), bitset_.begin(), bitset_.begin(), std::logical_and<bool>());
}

void BitSet::orOperator (const BitSet& bs) {
    assert(size() == bs.size());
    std::transform(bs.bitset_.begin(), bs.bitset_.end(), bitset_.begin(), bitset_.begin(), std::logical_or<bool>());    
}

void BitSet::xorOperator (const BitSet& bs) {
    assert(size() == bs.size());
    std::transform(bs.bitset_.begin(), bs.bitset_.end(), bitset_.begin(), bitset_.begin(), oopse::logical_xor<bool>());        
}
   
void BitSet::setBits(int fromIndex, int toIndex, bool value) {
    assert(fromIndex <= toIndex);
    assert(fromIndex >=0);
    assert(toIndex <= size());
    std::vector<char>::iterator first = bitset_.begin() + fromIndex;
    std::vector<char>::iterator last = bitset_.begin() + toIndex;
    std::fill(first, last, value);
}

void BitSet::resize(int nbits) {
    int oldSize = size();
    bitset_.resize(nbits);
    if (nbits > oldSize) {
        std::fill(bitset_.begin()+oldSize, bitset_.begin()+nbits+1, false);
    }
}

BitSet operator| (const BitSet& bs1, const BitSet& bs2) {
    assert(bs1.size() == bs2.size());

    BitSet result(bs1);
    result |= bs2;
    return result;
}

BitSet operator& (const BitSet& bs1, const BitSet& bs2) {
    assert(bs1.size() == bs2.size());

    BitSet result(bs1);
    result &= bs2;
    return result;
}

BitSet operator^ (const BitSet& bs1, const BitSet& bs2) {
    assert(bs1.size() == bs2.size());

    BitSet result(bs1);
    result ^= bs2;
    return result;
}

BitSet operator- (const BitSet& bs1, const BitSet& bs2) {
    assert(bs1.size() == bs2.size());

    BitSet result(bs1);
    result -= bs2;
    return result;
}

bool operator== (const BitSet & bs1, const BitSet &bs2) {
    assert(bs1.size() == bs2.size());
    return std::equal(bs1.bitset_.begin(), bs1.bitset_.end(), bs2.bitset_.begin());
}

std::istream& operator>> ( std::istream& is, const BitSet& bs) {

    return is;
}

std::ostream& operator<< ( std::ostream& os, const BitSet& bs) {
    for (int i = 0; i < bs.bitset_.size(); ++i) {
        std::string val = bs[i] ? "true" : "false";
        os << "BitSet[" << i <<"] = " << val << std::endl; 
    }
    
    return os;
}

}
