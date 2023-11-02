/*
 * Copyright (c) 2004-present, The University of Notre Dame. All rights
 * reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 * 1. Redistributions of source code must retain the above copyright notice,
 *    this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright notice,
 *    this list of conditions and the following disclaimer in the documentation
 *    and/or other materials provided with the distribution.
 *
 * 3. Neither the name of the copyright holder nor the names of its
 *    contributors may be used to endorse or promote products derived from
 *    this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
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

#include "utils/RandNumGen.hpp"

#ifdef IS_MPI
#include <mpi.h>
#endif

#include <numeric>
#include <random>
#include <vector>

namespace OpenMD {
  namespace Utils {

    RandNumGen::RandNumGen(result_type seed) {
      result_type result;
      int nProcessors {1};

#ifdef IS_MPI
      int worldRank {};
      MPI_Status status;

      MPI_Comm_rank(MPI_COMM_WORLD, &worldRank);
      MPI_Comm_size(MPI_COMM_WORLD, &nProcessors);

      if (worldRank == 0) {
#endif
        std::vector<result_type> initalSequence(nProcessors);
        std::iota(initalSequence.begin(), initalSequence.end(), seed);

        // Generate the seed_seq to initialize the RNG on each processor
        std::seed_seq seq(initalSequence.begin(), initalSequence.end());
        std::vector<result_type> seeds(nProcessors);
        seq.generate(seeds.begin(), seeds.end());

        result = seeds[0];

#ifdef IS_MPI
        for (int index {1}; index < nProcessors; ++index)
          MPI_Send(&seeds[index], 1, MPI_UINT32_T, index, 10, MPI_COMM_WORLD);
      } else {
        MPI_Recv(&result, 1, MPI_UINT32_T, 0, 10, MPI_COMM_WORLD, &status);
      }
#endif

      engine = std::mt19937(result);
    }
  }  // namespace Utils
}  // namespace OpenMD
