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
  * @file SimSnapshotManager.hpp
  * @author tlin
  * @date 10/20/2004
  * @time 23:56am
  * @version 1.0
  */
#ifndef BRAINS_SNAPSHOTMANAGER_HPP
#define BRAINS_SNAPSHOTMANAGER_HPP

#include "brains/SnapshotManager.hpp"

namespace oopse{

    /**
     * @class SimSnapshotManager SimSnapshotManager.hpp "brains/SimSnapshotManager.hpp"
     * @brief SimSnapshotManager class is the concrete snapshot manager for actual simulation
     * SimSnapshotManager only maintains two snapshots. 
     * @see PropSimSnapshotManager
     */
    class SimSnapshotManager : public SnapshotManager {
        public:
            SimSnapshotManager(SimModel* model);
            virtual bool advance();

            virtual Snapshot* getSnapshot(int id);
            
            int getCapacity();

            virtual void setCapacity(int capacity);
            
        private:
    };

}
#endif //BRAINS_SNAPSHOTMANAGER_HPP
