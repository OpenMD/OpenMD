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
 
#define _LARGEFILE_SOURCE64
#define _FILE_OFFSET_BITS 64

#include "io/StatWriter.hpp"
#include "utils/simError.h"

namespace oopse {
StatWriter::StatWriter( const std::string& filename, const StatsBitSet& mask) : mask_(mask){

#ifdef IS_MPI
  if(worldRank == 0 ){
#endif // is_mpi

    statfile_.open(filename.c_str(), std::ios::out | std::ios::trunc );
    
    if( !statfile_ ){
      
      sprintf( painCave.errMsg,
           "Could not open \"%s\" for stat output.\n",
           filename.c_str());
      painCave.isFatal = 1;
      simError();
    }

    writeTitle();
    
#ifdef IS_MPI
  }

  sprintf( checkPointMsg,
               "Sucessfully opened output file for stating.\n");
  MPIcheckPoint();
#endif // is_mpi

}

StatWriter::~StatWriter( ){

#ifdef IS_MPI
  if(worldRank == 0 ){
#endif // is_mpi

    statfile_.close();

#ifdef IS_MPI
  }
#endif // is_mpi
}

void StatWriter::writeTitle() {


#ifdef IS_MPI
    if(worldRank == 0 ){
#endif // is_mpi

        //write title
        statfile_ << "#";
        for (int i =0; i < mask_.size(); ++i) {
            if (mask_[i]) {
            statfile_ << "\t" << Stats::getTitle(i) << "(" << Stats::getUnits(i) << ")";
            }
        }
        statfile_ << std::endl;

#ifdef IS_MPI
  }
#endif // is_mpi    
}

void StatWriter::writeStat(const Stats& s){

#ifdef IS_MPI
    if(worldRank == 0 ){
#endif // is_mpi

        statfile_.precision(8);

        for (int i =0; i < mask_.size(); ++i) {
            if (mask_[i]) {
                statfile_ << "\t" << s[i];
            }
        }
        statfile_ << std::endl;

        statfile_.flush();

#ifdef IS_MPI
    }
#endif // is_mpi
}

}
