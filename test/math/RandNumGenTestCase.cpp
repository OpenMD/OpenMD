#include "math/RandNumGenTestCase.hpp"

// Registers the fixture into the 'registry'
CPPUNIT_TEST_SUITE_REGISTRATION( RandNumGenTestCase );

using namespace OpenMD;

void RandNumGenTestCase::testUniform(){
    MTRand randNumGen(823645754);

    const int N = 16;
    std::vector<unsigned long int> histogram(N, 0);
    const int num = 1000000;
    for (int i = 0; i <num; ++i) {
        ++histogram[randNumGen.randInt(N -1 )]; // rantInt returns an integer in [0, N-1]        
    }

    int avg = num / N;
    double tolerance = 0.01;
    for (int i = 0; i < num; ++i) {
        if ((histogram[i] - avg) /avg > tolerance) {

        }
    }
}

void RandNumGenTestCase::testGaussian(){      
    MTRand randNumGen(823645754);
    double mean = 3.0;
    double variance = 1.0;
    const int num = 1000000;
    double interval = 0.001;
    unsigned long int histogram[1000];
    for (int i = 0; i < num; ++i) {
        int index = static_cast<int>(randNumGen.randNorm(mean, variance) / interval);
        ++histogram[index];        
    }

    //fitting
}

void RandNumGenTestCase::testMPIRNG(){      
#ifdef IS_MPI
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

                assert(mpiRandNums[j], singleRandNums[j]);
            }

            
        }

        mpiRandNumGen.randInt();

        
    } else {

        unsigned long int randNum;
        for (int i = 0; i < nloops; ++i) {
            randNum = mpiRandNumGen.randInt();
            MPI_Send(&randNum, 1, MPI_INT, masterNode, i, MPI_COMM_WORLD);
        }

    }

#endif
}

