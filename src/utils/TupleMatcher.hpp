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
 * @file TupleMatcher.hpp
 * @author    tlin
 * @date  10/27/2004
 * @version 1.0
 */ 

#ifndef UTILS_TUPLEMATCHER_HPP
#define UTILS_TUPLEMATCHER_HPP

#include <map>
namespace oopse {

    template<class TupleType, class ReturnType>
    class TupleMatcher {
        static ReturnType match(const std::map<TupleType, ReturnType>& container, const TupleType& t) {
            typename std::map<TupleType, ReturnType>::iterator i;

            i = container.find(t);
            if (i != container.end()) {
                return i->second;
            } else {
                return NULL;
            }
        }
    };
    
}
#endif //UTILS_TUPLEMATCHER_HPP