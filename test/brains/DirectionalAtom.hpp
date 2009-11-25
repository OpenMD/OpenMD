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
 * @file DirectionalAtom.hpp
 * @author    tlin
 * @date  10/23/2004
 * @version 1.0
 */ 

#ifndef PRIMITIVES_DIRECTIONALATOM_HPP
#define PRIMITIVES_DIRECTIONALATOM_HPP

#include "primitives/Atom.hpp"
#include "types/DirectionalAtomType.hpp"
namespace OpenMD{
    class DirectionalAtom : public Atom {
        public:
            DirectionalAtom(DirectionalAtomType* dAtomType);
            /**
             * Returns the inertia tensor of this stuntdouble
             * @return the inertia tensor of this stuntdouble
             */ 
            virtual Mat3x3d getI();


           /**
             * Sets  the previous rotation matrix of this stuntdouble
             * @param a  new rotation matrix 
             */         
           virtual void setPrevA(const RotMat3x3d& a);
           
           /**
             * Sets  the current rotation matrix of this stuntdouble
             * @param a  new rotation matrix 
             */         
            virtual void setA(const RotMat3x3d& a);

           /**
             * Sets  the rotation matrix of this stuntdouble in specified snapshot
             * @param a rotation matrix to be set 
             * @param snapshotNo 
             * @see #getA
             */         
            virtual void setA(const RotMat3x3d& a, int snapshotNo);

            /** 
             * Left multiple rotation matrix by another rotation matrix 
             * @param m a rotation matrix
             */
            void rotateBy(const RotMat3x3d& m);
            
            /** Sets the internal unit frame of this stuntdouble by three euler angles */
            void setUnitFrameFromEuler(double phi, double theta, double psi);

            /**
             * Returns the gradient of this stuntdouble
             * @return the gradient of this stuntdouble
             */ 
            virtual std::vector<double> getGrad();

            virtual void accept(BaseVisitor* v);

        protected:
            Mat3x3d inertiaTensor_;   /**< inertial tensor */    
            RotMat3x3d sU_;               /**< body fixed standard unit vector */
    };

}//namespace OpenMD

#endif //PRIMITIVES_DIRECTIONALATOM_HPP
