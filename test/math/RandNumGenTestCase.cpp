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
 * research, please cite the following paper when you publish your work:
 *
 * [1] Drisko et al., J. Open Source Softw. 9, 7004 (2024).
 *
 * Good starting points for code and simulation methodology are:
 *
 * [2] Meineke, et al., J. Comp. Chem. 26, 252-271 (2005).
 * [3] Fennell & Gezelter, J. Chem. Phys. 124, 234104 (2006).
 * [4] Sun, Lin & Gezelter, J. Chem. Phys. 128, 234107 (2008).
 * [5] Vardeman, Stocker & Gezelter, J. Chem. Theory Comput. 7, 834 (2011).
 * [6] Kuang & Gezelter, Mol. Phys., 110, 691-701 (2012).
 * [7] Lamichhane, Gezelter & Newman, J. Chem. Phys. 141, 134109 (2014).
 * [8] Bhattarai, Newman & Gezelter, Phys. Rev. B 99, 094106 (2019).
 * [9] Drisko & Gezelter, J. Chem. Theory Comput. 20, 4986-4997 (2024).
 */

#include "math/RandNumGenTestCase.hpp"

// Registers the fixture into the 'registry'
CPPUNIT_TEST_SUITE_REGISTRATION(RandNumGenTestCase);

using namespace OpenMD;

void RandNumGenTestCase::testUniform() {
  MTRand randNumGen(823645754);

  const int N = 16;
  std::vector<unsigned long int> histogram(N, 0);
  const int num = 1000000;
  for (int i = 0; i < num; ++i) {
    ++histogram[randNumGen.randInt(
        N - 1)];  // rantInt returns an integer in [0, N-1]
  }

  int avg          = num / N;
  double tolerance = 0.01;
  for (int i = 0; i < num; ++i) {
    if ((histogram[i] - avg) / avg > tolerance) {}
  }
}

void RandNumGenTestCase::testGaussian() {
  MTRand randNumGen(823645754);
  double mean     = 3.0;
  double variance = 1.0;
  const int num   = 1000000;
  double interval = 0.001;
  unsigned long int histogram[1000];
  for (int i = 0; i < num; ++i) {
    int index =
        static_cast<int>(randNumGen.randNorm(mean, variance) / interval);
    ++histogram[index];
  }

  // fitting
}

void RandNumGenTestCase::testMPIRNG() {
#ifdef IS_MPI
  const int seed   = 324271632;
  const int nloops = 1000000;
  MPI_Status istatus;
  ParallelRandNumGen mpiRandNumGen(seed);
  const int primaryNode = 0;
  if (worldRank = primaryNode) {
    MTRand singleRandNumGen(seed);

    int nProcessors;
    MPI_Comm_size(MPI_COMM_WORLD, &nProcessors);
    std::vector<unsigned long int> mpiRandNums(nProcessors);
    std::vector<unsigned long int> singleRandNums(nProcessors);

    for (int i = 0; i < nloops; ++i) {
      mpiRandNums[primaryNode] = mpiRandNumGen.randInt();

      for (int j = 0; j < nProcessors; ++j) {
        if (j != primaryNode) {
          MPI_Recv(&mpiRandNums[j], 1, MPI_UNSIGNED_LONG, j, i, MPI_COMM_WORLD,
                   &istatus);
        }

        singleRandNums[j] = mpiRandNumGen.randInt();

        assert(mpiRandNums[j], singleRandNums[j]);
      }
    }

    mpiRandNumGen.randInt();

  } else {
    unsigned long int randNum;
    for (int i = 0; i < nloops; ++i) {
      randNum = mpiRandNumGen.randInt();
      MPI_Send(&randNum, 1, MPI_INT, primaryNode, i, MPI_COMM_WORLD);
    }
  }

#endif
}
