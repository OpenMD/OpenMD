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
using namespace oopse;
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

    std::transform(histogram.begin(), histogram.end(), normalizedHistogram.begin(), std::bind2nd(std::divides<double>(), num));
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
    const int masterNode = 0;
    int myRank;
    MPI_Comm_rank( MPI_COMM_WORLD, &myRank );
    if (myRank == masterNode) {

        SeqRandNumGen singleRandNumGen(seed);
        std::ofstream singleOs("single.dat");
	std::ofstream parallelOs("parallel.dat");
        int nProcessors;
        MPI_Comm_size(MPI_COMM_WORLD, &nProcessors);
        std::vector<unsigned long int> mpiRandNums(nProcessors);
        std::vector<unsigned long int> singleRandNums(nProcessors);

        for (unsigned long int i = 0; i < nloops; ++i) {
            mpiRandNums[masterNode] = mpiRandNumGen.randInt();
        
            for (int j = 0; j < nProcessors; ++j) {
                if (j != masterNode) {
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
            MPI_Send(&randNum, 1, MPI_INT, masterNode, i, MPI_COMM_WORLD);
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
