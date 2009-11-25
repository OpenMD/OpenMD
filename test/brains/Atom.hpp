/*
 * Copyright (C) 2000-2009  The Open Molecular Dynamics Engine (OpenMD) project
 * 
 * Contact: gezelter@openscience.org
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
 * @file Atom.hpp
 * @author    tlin
 * @date  10/22/2004
 * @version 1.0
 */ 

#ifndef PRIMITIVES_ATOM_HPP
#define PRIMITIVES_ATOM_HPP

#include "primitives/StuntDouble.hpp"
#include "types/AtomType.hpp"

namespace OpenMD{
    class Atom : public StuntDouble {
        public:
            Atom(AtomType* at);
            /**
             * Returns the inertia tensor of this stuntdouble
             * @return the inertia tensor of this stuntdouble
             */ 
            virtual Mat3x3d getI();

            /**
             * Returns the gradient of this stuntdouble
             * @return the inertia tensor of this stuntdouble
             */ 
            virtual std::vector<double> getGrad();

            virtual void accept(BaseVisitor* v);

            /** 
             * Returns the AtomType of this Atom.
             * @return the atom type of this atom
             */
            AtomType* getAtomType() {
                return atomType_;
            }
            
            //forward  functions of AtomType class
            bool    isCharge()  {
                return atomType_->isCharge(); 
            }
            
            bool    isDirectional() {
                return atomType_->isDirectional(); 
            }

            bool    isDipole()  { 
                return atomType_->isDipole(); 
            }
            
            bool    isGayBerne()  {
                return atomType_->isGayBerne(); 
            }
            
            bool    isSticky()  { 
                return atomType_->isSticky(); 
            }

            bool    isShape()  { 
                return atomType_->isShape(); 
            }            

        private:
            AtomType* atomType_;
       
    };

}//namepace OpenMD

#endif //PRIMITIVES_ATOM_HPP
