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
 * @file DirectionalAtom.hpp
 * @author    tlin
 * @date  10/23/2004
 * @version 1.0
 */ 

#ifndef PRIMITIVES_DIRECTIONALATOM_HPP
#define PRIMITIVES_DIRECTIONALATOM_HPP

#include "primitives/Atom.hpp"

namespace oopse{
    class DirectionalAtom : public Atom {
        public:
            DirectionalAtom();
            /**
             * Returns the inertia tensor of this stuntdouble
             * @return the inertia tensor of this stuntdouble
             * @see #setI
             */ 
            virtual Mat3x3d getI();

            /**
             * Sets the inertia tensor of this stuntdouble
             * @param trq new inertia tensor
             * @see #getI
             */      
            virtual void setI(Mat3x3d& I);

            /**
             * Returns the gradient of this stuntdouble
             * @return the inertia tensor of this stuntdouble
             * @see #setI
             */ 
            virtual std::vector<double> getGrad();

            virtual void accept(BaseVisitor* v);

        protected:
            Mat3x3d inertiaTensor_;       
    };

}//namepace oopse

#endif //PRIMITIVES_DIRECTIONALATOM_HPP
