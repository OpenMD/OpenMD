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
 * @file FileReader.hpp
 * @author Teng Lin
 * @date 10/14/2004
 * @version 1.0
 */

#ifndef IO_FILEREADER_HPP
#define IO_FILEREADER_HPP

#include <fstream>

namespace oopse{
    /**
     * @class FileReader FileReader.hpp "io/FileReader.hpp"
     * @brief an single/parallel File Reader class
     * @warning In parallel mode, master node will read the whole file
     * into memory, and then pass it to other slave nodes.
     * 
     */
    class FileReader {
        public:

            /**
             *
             */
            FileReader(const string& filename);

            /**
             *
             */
            virtual ~FileReader();

            /**
             *
             */
            int open();

            /*
             *
             */
            int read();

            /**
             *
             */
            int close();
            
        protected:
            ifstream input;
    };

}
#endif //IO_FILEREADER_HPP