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
  * @file SimModel.hpp
  * @author tlin
  * @date 10/24/2004
  * @time 22:17am
  * @version 1.0
  */

#ifndef BRAINS_SIMMODEL_HPP
#define BRAINS_SIMMODEL_HPP

namespace oopse {

class SimModel {
    public:

    private:
        PropertyMap properties_;
        SnapshotManager* snapshotMan_;
        vector<Molecule*> molecules_;
};


} //namespace oopse
#endif //BRAINS_SIMMODEL_HPP
  
