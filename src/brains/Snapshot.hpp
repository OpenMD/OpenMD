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
  * @file Snapshot.hpp
  * @author tlin
  * @date 10/20/2004
  * @time 23:56am
  * @version 1.0
  */
#ifndef BRAINS_SNAPSHOT_HPP
#define BRAINS_SNAPSHOT_HPP

#include <vector>

#include "math/Vector3.hpp"
#include "math/SquareMatrix3.hpp"

using namespace std;

namespace oopse{

    /**
     * @class Snapshot Snapshot.hpp "brains/Snapshot.hpp"
     * @brief Snapshot class is a repository class for storing dynamic data during 
     *  Simulation
     * @see SimData
     */
    class Snapshot {
        public:

            Snapshot(int i) {

            }

            Snapshot(const Snapshot& s);

            Snapshot& operator =(const Snapshot& s);
            
            /** Returns the id of this Snapshot */
            int getID() {
                return id_;
            }

            /** Sets the id of this Snapshot */
            void setID(int id) {
                id_ = id;
            }

            /** */
            Snapshot* clone() {
                return new Snapshot(*this);
            }


            //template<typename T>
            //static typename T::ElemPointerType getArrayPointer(vector<T>& v) {
            //    return v[0]->getArrayPointer();
            //}

            static double* getArrayPointer(vector<Vector3d>& v) {
                return v[0].getArrayPointer();
            }
            
            static double* getArrayPointer(vector<RotMat3x3d>& v) {
                return v[0].getArrayPointer();
            }
            
            static double* getArrayPointer(vector<double>& v) {
                assert(v.size() > 0);
                return &(v[0]);
            }
            
            vector<Vector3d> pos;
            vector<Vector3d> vel;
            vector<Vector3d> frc;
            vector<Vector3d> trq;
            vector<RotMat3x3d> Amat;
            vector<Vector3d> mu;
            vector<Vector3d> ul;
            vector<double> zAngle;
            
        private:

            int id_; /**< identification number of the snapshot */
    };

}
#endif //BRAINS_SNAPSHOT_HPP
