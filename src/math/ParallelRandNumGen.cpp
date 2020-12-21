/*
 * Copyright (c) 2004-2020 The University of Notre Dame. All Rights Reserved.
 *
 * The University of Notre Dame grants you ("Licensee") a
 * non-exclusive, royalty free, license to use, modify and
 * redistribute this software in source and binary code form, provided
 * that the following conditions are met:
 *
 * 1. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright
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
 *
 * SUPPORT OPEN SCIENCE!  If you use OpenMD or its source code in your
 * research, please cite the appropriate papers when you publish your
 * work.  Good starting points are:
 *
 * [1] Meineke, et al., J. Comp. Chem. 26, 252-271 (2005).
 * [2] Fennell & Gezelter, J. Chem. Phys. 124, 234104 (2006).
 * [3] Sun, Lin & Gezelter, J. Chem. Phys. 128, 234107 (2008).
 * [4] Vardeman, Stocker & Gezelter, J. Chem. Theory Comput. 7, 834 (2011).
 * [5] Kuang & Gezelter, Mol. Phys., 110, 691-701 (2012).
 * [6] Lamichhane, Gezelter & Newman, J. Chem. Phys. 141, 134109 (2014).
 * [7] Lamichhane, Newman & Gezelter, J. Chem. Phys. 141, 134110 (2014).
 * [8] Bhattarai, Newman & Gezelter, Phys. Rev. B 99, 094106 (2019).
 */

#ifdef IS_MPI
#include <mpi.h>
#endif

#include "math/ParallelRandNumGen.hpp"
#include "utils/MemoryUtils.hpp"

namespace OpenMD {

  int ParallelRandNumGen::nCreatedRNG_ = 0;

  ParallelRandNumGen::ParallelRandNumGen(const uint32& oneSeed) {

    unsigned long seed = oneSeed;

#ifdef IS_MPI
    const int primaryNode = 0;
    MPI_Bcast(&seed, 1, MPI_UNSIGNED_LONG, primaryNode, MPI_COMM_WORLD); 
#endif

    if (seed != oneSeed) {
      sprintf(painCave.errMsg,
	      "Using different seed to initialize ParallelRandNumGen.\n");
      painCave.isFatal = 1;;
      simError();
    }

    int nProcessors;
#ifdef IS_MPI
    MPI_Comm_size( MPI_COMM_WORLD, &nProcessors);
    MPI_Comm_rank( MPI_COMM_WORLD, &myRank_);

#else
    nProcessors = 1;
    myRank_ = 0;
#endif
    //In order to generate independent random number stream, the
    //actual seed used by random number generator is the seed passed
    //to the constructor plus the number of random number generators
    //which are already created.
    unsigned long newSeed = oneSeed + nCreatedRNG_;
    mtRand_ = MemoryUtils::make_unique<MTRand>(newSeed, nProcessors, myRank_);
    
    ++nCreatedRNG_;
  }

  ParallelRandNumGen::ParallelRandNumGen() {

    int nProcessors;
#ifdef IS_MPI
    MPI_Comm_size( MPI_COMM_WORLD, &nProcessors);
    MPI_Comm_rank( MPI_COMM_WORLD, &myRank_);
#else
    nProcessors = 1;
    myRank_ = 0;
#endif
    mtRand_ = MemoryUtils::make_unique<MTRand>(nProcessors, myRank_);

    seed();       /** @todo calling virtual function in constructor is
                      not a good design */
  }


  void ParallelRandNumGen::seed( const uint32 oneSeed ) {

    unsigned long seed = oneSeed;
#ifdef IS_MPI
    const int primaryNode = 0;
    MPI_Bcast(&seed, 1, MPI_UNSIGNED_LONG, primaryNode, MPI_COMM_WORLD); 
#endif
    if (seed != oneSeed) {
      sprintf(painCave.errMsg,
	      "Using different seed to initialize ParallelRandNumGen.\n");
      painCave.isFatal = 1;;
      simError();
    }
    
    unsigned long newSeed = oneSeed +nCreatedRNG_;
    mtRand_->seed(newSeed);
    
    ++nCreatedRNG_;
  }
        
  void ParallelRandNumGen::seed() {

    std::vector<uint32> bigSeed;

#ifdef IS_MPI
    int size;
    const int primaryNode = 0;
    if (worldRank == primaryNode) {
#endif

      bigSeed = mtRand_->generateSeeds();

#ifdef IS_MPI
      size = bigSeed.size();
      MPI_Bcast(&size, 1, MPI_INT, primaryNode, MPI_COMM_WORLD);
      MPI_Bcast(&bigSeed[0], size, MPI_UNSIGNED_LONG, primaryNode, 
                MPI_COMM_WORLD);
    }else {
      MPI_Bcast(&size, 1, MPI_INT, primaryNode, MPI_COMM_WORLD);
      bigSeed.resize(size);
      MPI_Bcast(&bigSeed[0], size, MPI_UNSIGNED_LONG, primaryNode, 
                MPI_COMM_WORLD);
    }
#endif
    
    if (bigSeed.size() == 1) {
      mtRand_->seed(bigSeed[0]);
    } else {
      mtRand_->seed(&bigSeed[0], bigSeed.size());
    }

    ++nCreatedRNG_;
  }        
}
