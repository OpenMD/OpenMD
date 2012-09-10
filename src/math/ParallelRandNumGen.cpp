/*
 * Copyright (c) 2005 The University of Notre Dame. All Rights Reserved.
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
 * [1]  Meineke, et al., J. Comp. Chem. 26, 252-271 (2005).             
 * [2]  Fennell & Gezelter, J. Chem. Phys. 124, 234104 (2006).          
 * [3]  Sun, Lin & Gezelter, J. Chem. Phys. 128, 24107 (2008).          
 * [4]  Kuang & Gezelter,  J. Chem. Phys. 133, 164101 (2010).
 * [5]  Vardeman, Stocker & Gezelter, J. Chem. Theory Comput. 7, 834 (2011).
 */

#include "math/ParallelRandNumGen.hpp"
#ifdef IS_MPI
#include <mpi.h>
#endif

namespace OpenMD {

  int ParallelRandNumGen::nCreatedRNG_ = 0;

  ParallelRandNumGen::ParallelRandNumGen(const uint32& oneSeed) {

    unsigned long seed = oneSeed;

#ifdef IS_MPI
    const int masterNode = 0;
    MPI::COMM_WORLD.Bcast(&seed, 1, MPI::UNSIGNED_LONG, masterNode); 
#endif

    if (seed != oneSeed) {
      sprintf(painCave.errMsg,
	      "Using different seed to initialize ParallelRandNumGen.\n");
      painCave.isFatal = 1;;
      simError();
    }

    int nProcessors;
#ifdef IS_MPI
    nProcessors = MPI::COMM_WORLD.Get_size();
    myRank_ = MPI::COMM_WORLD.Get_rank();
#else
    nProcessors = 1;
    myRank_ = 0;
#endif
    //In order to generate independent random number stream, the
    //actual seed used by random number generator is the seed passed
    //to the constructor plus the number of random number generators
    //which are already created.
    unsigned long newSeed = oneSeed + nCreatedRNG_;
    mtRand_ = new MTRand(newSeed, nProcessors, myRank_);
    
    ++nCreatedRNG_;
  }

  ParallelRandNumGen::ParallelRandNumGen() {

    std::vector<uint32> bigSeed;
    int nProcessors;
#ifdef IS_MPI
    nProcessors = MPI::COMM_WORLD.Get_size();
    myRank_ = MPI::COMM_WORLD.Get_rank();
#else
    nProcessors = 1;
    myRank_ = 0;
#endif
    mtRand_ = new MTRand(nProcessors, myRank_);

    seed();       /** @todo calling virtual function in constructor is
                      not a good design */
  }


  void ParallelRandNumGen::seed( const uint32 oneSeed ) {

    unsigned long seed = oneSeed;
#ifdef IS_MPI
    const int masterNode = 0;
    MPI::COMM_WORLD.Bcast(&seed, 1, MPI::UNSIGNED_LONG, masterNode); 
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
    const int masterNode = 0;
    if (worldRank == masterNode) {
#endif

      bigSeed = mtRand_->generateSeeds();

#ifdef IS_MPI
      size = bigSeed.size();
      MPI::COMM_WORLD.Bcast(&size, 1, MPI::INT, masterNode);
      MPI::COMM_WORLD.Bcast(&bigSeed[0], size, MPI::UNSIGNED_LONG, masterNode);
    }else {
      MPI::COMM_WORLD.Bcast(&size, 1, MPI::INT, masterNode);
      bigSeed.resize(size);
      MPI::COMM_WORLD.Bcast(&bigSeed[0], size, MPI::UNSIGNED_LONG, masterNode);
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
