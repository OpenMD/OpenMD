#include <cassert>
#include <fstream>
#include <algorithm>
#include "utils/simError.h"
#include "math/ParallelRandNumGen.hpp"

void RandNumGenTestCase::testUniform(){
    MTRand randNumGen(823645754);

    const int N = 16;
    std::vector<unsigned long int> histogram(N, 0);
    const int num = 1000000;
    for (int i = 0; i <num; ++i) {
        ++histogram[randNumGen.randInt(N -1 )]; // rantInt returns an integer in [0, N-1]        
    }

    std::ofstream uniform("uniform.dat")
    int avg = num / N;
    double tolerance = 0.01*avg;
    for (int i = 0; i < num; ++i) {
        assert((histogram[i] - avg) /avg <= tolerance);
        uniform << i << "\t" << histogram[i] << std::endl;
    }
}

void RandNumGenTestCase::testGaussian(){      
    MTRand randNumGen(823645754);
    double mean = 100.0;
    double variance = 1.0;
    const int num = 1000000;
    double interval = 0.1;
    const int size = 2000;
    vector<unsigned long int> histogram(size , 0);
    vector<double> normalizedHistogram(size);
    for (int i = 0; i < num; ++i) {
        int index = static_cast<int>(randNumGen.randNorm(mean, variance) / interval);
        ++histogram[index];        
    }

    std::transform(histogram.begin(), histogram.end(), normalizedHistogram.begin(), std::bind2nd(std::divides<double>(), num));
    std::ofstream gaussian("gaussian.dat");
    for (int i = 0; i < num; ++i) {
        gaussian << i << "\t" << normalizedHistogram[i] << std::endl;
    }    
}

void RandNumGenTestCase::testParallelRandNumGen(){      
    const int seed = 324271632;
    const int nloops = 1000000;
    MPI_Status istatus;
    ParallelRandNumGen mpiRandNumGen(seed);
    const int masterNode = 0;
    if (worldRank = masterNode) {

        MTRand singleRandNumGen(seed);

        int nProcessors;
        MPI_Comm_size(MPI_COMM_WORLD, &nProcessors);
        std::vector<unsigned long int> mpiRandNums(nProcessors);
        std::vector<unsigned long int> singleRandNums(nProcessors);

        for (int i = 0; i < nloops; ++i) {
            mpiRandNums[masterNode] = mpiRandNumGen.randInt();
        
            for (int j = 0; j < nProcessors; ++j) {
                if (j != masterNode) {
                    MPI_Recv(&mpiRandNums[j], 1, MPI_UNSIGNED_LONG, j, i, MPI_COMM_WORLD, &istatus);
                }

                singleRandNums[j] = mpiRandNumGen.randInt();
            }

            assert(mpiRandNums, singleRandNums);
        }


       
    } else {

        unsigned long int randNum;
        for (int i = 0; i < nloops; ++i) {
            randNum = mpiRandNumGen.randInt();
            MPI_Send(&randNum, 1, MPI_INT, masterNode, i, MPI_COMM_WORLD);
        }

    }

}


int main(int argc, char* argv[]) {

    MPI_Init(argc, argv);

    if (worldRank == 0 ) {
        testUniform();
        testGaussian();
    }

    testParallelRandNumGen();
    
    MPI_Finalize();
}
