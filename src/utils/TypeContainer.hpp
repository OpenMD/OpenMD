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
  * @file TypeContainer.hpp
  * @author tlin
  * @date 11/04/2004
  * @time 13:51am
  * @version 1.0
  */

#ifndef UTILS_TYPECONTAINER_HPP
#define UTILS_TYPECONTAINER_HPP
#include <algorithm>
#include <map>
#include <vector>

#include "utils/Utility.hpp"

namespace oopse {

    /**
     * @class TypeContainer TypeContainer.hpp "utils/TypeContainer.hpp"
     */
    template<class ElemType, int SIZE>
    class TypeContainer {
        public:
            
            typedef ElemType* ElemPtr;
            typedef std::vector<std::string> KeyType;
            typedef typename KeyType::iterator KeyTypeIterator;
            typedef std::pair<int, ElemPtr> ValueType;
            typedef std::map<KeyType, ValueType> MapType;
            typedef typename std::map<KeyType, ValueType>::iterator MapTypeIterator;

            TypeContainer() : index_(0) {}
            
            ~TypeContainer() {
                MapTypeIterator i;
                for (i = data_.begin(); i != data_.end(); ++i) {
                    delete (i->second).second;
                }
                data_.clear();
            }
            
            bool add(KeyType& keys, ElemPtr elem) {
                assert(keys.size() == SIZE);
                assert(elem);
                return data_.insert(MapType::value_type(keys, std::make_pair(index_++,elem))).second;
            }

            /** Exact Match */
            ElemPtr find(KeyType& keys) {
                assert(keys.size() == SIZE);
                MapTypeIterator i;
                ValueType foundType;
                
                i = data_.find(keys);
                if (i != data_.end()) {
                    return (i->second).second;
                }

                KeyType reversedKeys = keys;
                std::reverse(reversedKeys.begin(), reversedKeys.end());

                i = data_.find(reversedKeys);
                if (i != data_.end()) {
                    return (i->second).second;
                } else {
                    return NULL;
                }
                
            }

            /**
             * @todo 
             */
            ElemPtr find(KeyType& keys, const std::string& wildCard) {
                assert(keys.size() == SIZE);
                
                std::vector<KeyTypeIterator> iterCont;
                KeyType replacedKey;
                MapTypeIterator i;
                std::vector<ValueType> foundTypes;
                
                while (replaceWithWildCard(iterCont, keys, replacedKey, wildCard)) {                    
                    i = data_.find(replacedKey);
                    if (i != data_.end()) {
                        foundTypes.push_back(i->second);
                    }
                }

                //reverse the order of keys
                KeyType reversedKeys = keys;
                std::reverse(reversedKeys.begin(), reversedKeys.end());

                //if the reversedKeys is the same as keys, just skip it
                if (reversedKeys != keys) {

                    
                    //empty the iterator container
                    iterCont.clear();

                    while (replaceWithWildCard(iterCont, reversedKeys, replacedKey, wildCard)) {
                        i = data_.find(replacedKey);
                        if (i != data_.end()) {
                            foundTypes.push_back(i->second);
                        }
                    }

                    //replaceWithWildCard can not generate this particular sequence, we have to 
                    //do it manually                
                    KeyType allWildCards(SIZE, wildCard);
                    i = data_.find(replacedKey);
                    if (i != data_.end()) {
                        foundTypes.push_back(i->second);
                    }

                }

                typename std::vector<ValueType>::iterator j;                
                j = std::min_element(foundTypes.begin(), foundTypes.end());

                return j == foundTypes.end() ? NULL : j->second;
            }

            unsigned int size() {
                return data_.size();
            }

            ElemPtr beginType(MapTypeIterator& i) {
                i = data_.begin();
                return i  != data_.end() ? (i->second).second : NULL;
            }

            ElemPtr nextType(MapTypeIterator& i) {
                ++i;
                return i  != data_.end() ? (i->second).second : NULL;
            }
            
        private:
            int index_;
            MapType data_;        
            
    };


}//end namespace oopse

#endif //UTILS_TYPECONTAINER_HPP
