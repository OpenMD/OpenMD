
#ifndef PARALLEL_PARALLEL_HPP
#define PARALLEL_PARALLEL_HPP
#include <vector>
#include <string>

//#include <protomol/type/Real.h>
//#include <protomol/parallel/ParallelType.h>

namespace OpenMD {
  class ScalarStructure;
  class Vector3DBlock;
  class GenericTopology;

  //____ Parallel
  class Parallel  {

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Constructors, destructors, assignment
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  private:
    Parallel();
    ~Parallel();
    Parallel(const Parallel &) {}

    Parallel &operator=(const Parallel &) {return *this;}

  public:
    static void init(int &argc, char ** &argv); ///< MPI_Init and sets id and
                                                ///< number of process.
    static void finalize();                     ///< MPI_Finialize.

    static bool initialized()          {return myInitialized;}

    static bool finalized()            {return myFinalized;}

    static bool ok()                   {return myInitialized && !myFinalized;}

    static int getId()                {return myId;}

    static int getMasterId()          {return myMasterId;}

    static int getNum()               {return myNum;}

    static int getAvailableId()       {return myAvailableId;}

    static int getAvailableNum()      {return myAvailableNum;}

    static bool isParallel()           {return myIsParallel;}

    static bool iAmMaster()            {return myIAmMaster;}

    static bool iAmSlave()             {return myIAmSlave;}


  public:

    /**
     * Broadcast a Vector3DBlock from the master to the rest.
     * Handles also the case if the receiving nodes have an emty
     * (other size) Vector3DBlock.
     */
    static void bcast(Vector3DBlock *coords);
    static void bcast(int &n);

    /**
     * Does the preprocessing, clearing the energies and forces for
     * all nodes except the master.
     */
    static void distribute(ScalarStructure *energy, Vector3DBlock *coords);

    /**
     * Does a reduction and a broadcast, such that all
     * nodes have the correct summed energies and forces locally.
     */
    static void reduce(ScalarStructure *energy, Vector3DBlock *coords);

    /// Takes care of sending a Vector3DBlock to a single node.
    static void send(Vector3DBlock *vect, int address);

    /// Recieves a Vector3DBlock from a node.
    static void recv(Vector3DBlock *vect, int address);

    /// Sends a Vector3DBlock, to one node, and overwrites it with one recieved
    /// from another node
    static void sendrecv_replace(Vector3DBlock *Vect, int sendaddr,
                                 int recvaddr);

    /// Sends a Real (or array of Reals) from one node to another
    static void send(RealType *data, int num, int address);

    /// Recieves a Real (or array of Reals) from another node
    static void recv(RealType *data, int num, int address);

    /// Sends one real and receives another real from one or two nodes (sends
    ///  to one, recieves to the other)
    static void sendrecv(Real *senddata, int sendnum, int sendaddr,
                         Real *recvdata, int recvnum,
                         int recvaddr);

    /// Sends a Real (or array of Reals) to one node and overwrites it with
    ///  data from another node.
    static void sendrecv_replace(Real *data, int num, int sendaddr,
                                 int recvaddr);

    /// Gathers a vector of Reals from the compute nodes and stores it in an
    /// array.
    static void gather(Real *data, int num, Real *data_array, int address);

    /// Gathers a vector of Reals from the compute nodes and stores it in an
    /// array, distributed to all the nodes.
    static void allgather(Real *data, int num, Real *data_array);

    // Partitioning
  public:
    /// Returns the number of packages over all forces.
    static unsigned int getNumberOfPackages(unsigned int n);

  public:
    /// true if the actual work package should computed, false skip it
    static bool next();

    static void resetNext() {myNext = 0; myNextRange[0] = 0; myNextRange[1] = 0;
                             myWorkState = getWorkState();}

    static void resetNext(const std::vector<int> &blocks);

  public:
    /// Disables MPI on the current node
    static void isolateNode();

    /// Enables MPI on the current node
    static void integrateNode();

  private:
    static void kill();
    static Parallel &instance();
    static bool ok(const std::string &err);
    static WorkState getWorkState();
    static void nextMaster();

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // My data members
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    static const bool isMPI;
  private:
    static bool myInitialized;   ///< If Parallel is initialized
    static bool myFinalized;     // /<If Parallel is finalized
    static int myId;            ///< Actual id of the node, [0,1,...,N-1]
    static int myMasterId;      ///< Id of the master
    static int myNum;           ///< Number of nodes, N
    static int myAvailableId;   ///< Actual master-slave id, -1 for the master
    static int myAvailableNum;  ///< Available number of nodes for computation
    static bool myIsParallel;    ///< If environment has more that 1 node
    static bool myIAmMaster;     ///< If this node is master
    static bool myIAmSlave;      ///< If this node is slave/worker
    static ParallelType myMode;          ///< Parallelization scheme
    static WorkState myWorkState;     ///< Actual work state

    static int myPipeSize;      ///< Number of add. work packages to push to
                                ///< slaves
    static bool myUseBarrier;    ///< Flag to signal usage of MPI_Barrier
    static int myMaxPackages;   ///< Number of max. packages per node per force

    static int *myBuffer;        ///< Bsend buffer
    static int myNext;          ///< Counter of next() calls
    static int myNextRange[2];  ///< Actual work package to work on [from,to]
    static std::vector<int> myDone;          ///< Master, keeps track of slave
                                             ///< have got already their work
    static std::vector<int> myBlockList;     ///< List of of the force
                                             ///< partitioning

    static int myRecv;          ///< Number of outstanding receives
    static int myI;             ///< Index of the actual work package (only
                                ///< master)
    static int myP;             ///< Node number

    static Parallel *obj;             ///< Instance

    static bool myIsolated;      ///< Saves isolation state
    static int myOldId;         ///< Saves processor id when switching MPI off
    static int myOldNum;        ///< Saves number of processors when switching
                                ///< MPI off
    static ParallelType myOldMode;  ///< Saves parallelization scheme when
                                    ///< switching MPI off
  };
  //____ INLINES
}
#endif /* PARAMETER_H */
