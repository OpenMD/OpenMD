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

#include <functional>
#include <cassert>
#include <fstream>
#include <algorithm>
#include <vector>
#include <iostream>

#include "math/SeqRandNumGen.hpp"
#ifdef IS_MPI
#include <mpi.h>
#include "math/ParallelRandNumGen.hpp"
#endif
using namespace OpenMD;
using namespace std;

void testUniform(){
    SeqRandNumGen randNumGen(823645754);
    const int N = 100;
    std::vector<unsigned long int> histogram(N, 0);
    const unsigned long int num = 10000000;
    for (int i = 0; i <num; ++i) {
	int index = randNumGen.randInt(N -1 );   
        ++histogram[index]; // rantInt returns an integer in [0, N-1]        
    }
    std::ofstream uniform("uniform.dat");
    double avg = num / N;
    double tolerance = 0.01*avg;
    for (int i = 0; i < N; ++i) {
        //assert((histogram[i] - avg) /avg <= tolerance);
        uniform << i << "\t" << histogram[i] << std::endl;
    }
}

void testGaussian(){      
    SeqRandNumGen randNumGen(823645754);
    double mean = 100.0;
    double variance = 1.0;
    const unsigned long int num = 1000000;
    double interval = 0.1;
    const int size = 2000;
    vector<unsigned long int> histogram(size , 0);
    vector<double> normalizedHistogram(size);
    for (unsigned long int i = 0; i < num; ++i) {
        int index = static_cast<int>(randNumGen.randNorm(mean, variance) / interval);
        ++histogram[index];        
    }

    std::transform(histogram.begin(), histogram.end(), normalizedHistogram.begin(), std::bind(std::divides<double>(), placeholders::_1, num));
    std::ofstream gaussian("gaussian.dat");
    for (int i = 0; i < size; ++i) {
        gaussian << i*interval << "\t" << normalizedHistogram[i] << std::endl;
    }    
}
#ifdef IS_MPI
void testParallelRandNumGen(){      
    const unsigned long int seed = 324271632;
    const unsigned long int nloops = 1000000;
    MPI_Status istatus;
    ParallelRandNumGen mpiRandNumGen(seed);
    const int primaryNode = 0;
    int myRank;
    MPI_Comm_rank( MPI_COMM_WORLD, &myRank );
    if (myRank == primaryNode) {

        SeqRandNumGen singleRandNumGen(seed);
        std::ofstream singleOs("single.dat");
	std::ofstream parallelOs("parallel.dat");
        int nProcessors;
        MPI_Comm_size(MPI_COMM_WORLD, &nProcessors);
        std::vector<unsigned long int> mpiRandNums(nProcessors);
        std::vector<unsigned long int> singleRandNums(nProcessors);

        for (unsigned long int i = 0; i < nloops; ++i) {
            mpiRandNums[primaryNode] = mpiRandNumGen.randInt();
        
            for (int j = 0; j < nProcessors; ++j) {
                if (j != primaryNode) {
                    MPI_Recv(&mpiRandNums[j], 1, MPI_UNSIGNED_LONG, j, i, MPI_COMM_WORLD, &istatus);
                }

                singleRandNums[j] = mpiRandNumGen.randInt();
            }

	    for (int j = 0; j < nProcessors; ++j) {
                singleOs << singleRandNums[j] << "\n";
	        parallelOs << singleRandNums[j] << "\n";
	    }
        }


       
    } else {

        unsigned long int randNum;
        for (unsigned long int i = 0; i < nloops; ++i) {
            randNum = mpiRandNumGen.randInt();
            MPI_Send(&randNum, 1, MPI_INT, primaryNode, i, MPI_COMM_WORLD);
        }

    }

}
#endif
int main(int argc, char* argv[]) {
#ifdef IS_MPI
    MPI_Init(&argc, &argv);
    std::cout << "begin test" << std::endl;
    if (worldRank == 0 ) {
        testUniform();
        testGaussian();
    }

    testParallelRandNumGen();
    
    MPI_Finalize();
#else
   std::cout << "begin test" <<std::endl;	
   testUniform();
   testGaussian();
#endif  
   std::cout << "test done" << std::endl;  
   return 0;
}
