/**
 * @file Parallel.cpp
 * @author Charles Vardeman <cvardema.at.nd.edu>
 * @date 08/18/2010
 * @time 11:56am
 * @version 1.0
 *
 * @section LICENSE
 * Copyright (C) 2010 The University of Notre Dame. All Rights Reserved.
 *
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
 * [4]  Vardeman & Gezelter, in progress (2009).                  
 */


#include <stdlib.h>
#include "config.h"

#ifdef IS_MPI
#include <mpi.h>
#endif

#include <iostream>
#include <vector>
#include <algorithm>

#include "parallel/Parallel.hpp"


using namespace std;
using namespace OpenMD;

//#define DEBUG_PARALLEL


#ifdef SINGLE_PRECISION
#define MY_MPI_REAL MPI_FLOAT
#else
#define MY_MPI_REAL MPI_DOUBLE
#endif




//____ mpiAbort
static void mpiAbort();

void mpiAbort() {
  if (Parallel::ok()) {
#ifdef IS_MPI
    MPI::Comm_world.Abort(1);
#endif
  }
  exit(EXIT_FAILURE);
}

//____ mpiExit
static void mpiExit();

void mpiExit() {
  if (Parallel::ok())
    Parallel::finalize();
  exit(EXIT_SUCCESS);
}


//____ MPITypeTraits
template<typename T>
struct MPITypeTraits;

#ifdef IS_MPI
template<>
struct MPITypeTraits<RealType> {
  static const MPI::Datatype datatype;
};
const MPI_Datatype MPITypeTraits<RealType>::datatype = MY_MPI_REAL;

template<>
struct MPITypeTraits<int> {
  static const MPI::Datatype datatype;
};
const MPI::Datatype MPITypeTraits<int>::datatype = MPI_INT;


//____ allReduceScalar
template<typename T>
void allReduceScalar(T &begin) {
  T tmp = begin;
  /*
  MPI_Allreduce(&tmp, &begin, 1, MPITypeTraits<T>::datatype, MPI_SUM,
                (exludeMaster ? slaveComm : MPI_COMM_WORLD));
  */
  MPI::Comm_world.Allreduce(&tmp,&begin,1, MPITypeTraits<T>::datatype, MPI::SUM);
}

//____ allReduce
template<typename T>
void allReduce(T *begin, T *end) {
  vector<T> tmp(end - begin);
  copy(begin, end, tmp.begin());
  /*  MPI_Allreduce(&(tmp[0]), begin, (end - begin), MPITypeTraits<T>::datatype,
                MPI_SUM,
                (exludeMaster ? slaveComm : MPI_COMM_WORLD));
  */
   MPI::Comm_world.Allreduce(&tmp[0],&begin,1, MPITypeTraits<T>::datatype, MPI::SUM);
}

template<bool exludeMaster, bool dobarrier>
void allReduce(Vector3DBlock *coords) {
  allReduce<exludeMaster, dobarrier>(&(coords->begin()->x),
                                     &(coords->end()->x));
}

template<bool exludeMaster, bool dobarrier>
void allReduce(ScalarStructure *energies) {
  allReduce<exludeMaster, dobarrier>(
    &((*energies)[ScalarStructure::FIRST]),
    &((*energies)[ScalarStructure::LASTREDUCE]));
}

//____ broadcastScalar
template<bool exludeMaster, bool dobarrier, typename T>
void broadcastScalar(T &begin) {
  if (dobarrier)
    doBarrier<exludeMaster>();
  MPI_Bcast(&begin, 1, MPITypeTraits<T>::datatype, master,
            (exludeMaster ? slaveComm : MPI_COMM_WORLD));
}

//____ broadcast

template<bool exludeMaster, bool dobarrier, typename T>
void broadcast(T *begin, T *end) {
  if (dobarrier)
    doBarrier<exludeMaster>();
  MPI_Bcast(begin, (end - begin), MPITypeTraits<T>::datatype,
            (exludeMaster ? 0 : master),
            (exludeMaster ? slaveComm : MPI_COMM_WORLD));
}

template<bool exludeMaster, bool dobarrier>
void broadcast(Vector3DBlock *coords) {
  broadcast<exludeMaster, dobarrier>(&(coords->begin()->x),
                                     &(coords->end()->x));
}

template<bool exludeMaster, bool dobarrier>
void broadcast(ScalarStructure *energies) {
  broadcast<exludeMaster, dobarrier>(
    &((*energies)[ScalarStructure::FIRST]),
    &((*energies)[ScalarStructure::LASTREDUCE]));
}

#endif

//____ Parallel

#ifdef IS_MPI
const bool Parallel::isMPI = true;
#else
const bool Parallel::isMPI = false;
#endif

bool Parallel::myInitialized = false;
bool Parallel::myFinalized = false;
int Parallel::myId = 0;
int Parallel::myMasterId = master;
int Parallel::myNum = 1;
int Parallel::myAvailableId = 0;
int Parallel::myAvailableNum = 1;
bool Parallel::myIsParallel = Parallel::isMPI;
bool Parallel::myIAmMaster = true;
bool Parallel::myIAmSlave = true;
ParallelType Parallel::myMode = ParallelType::STATIC;
Parallel::WorkState Parallel::myWorkState = Parallel::SEQUENTIAL;

int Parallel::myPipeSize = 1;
bool Parallel::myUseBarrier = false;
int Parallel::myMaxPackages = -1;

int *Parallel::myBuffer = NULL;
int Parallel::myNext = 0;
int Parallel::myNextRange[2] = {
  0, 0
};
vector<int>  Parallel::myDone;
vector<int>  Parallel::myBlockList;

int Parallel::myRecv;
int Parallel::myI;
int Parallel::myP;

Parallel * Parallel::obj = NULL;

int Parallel::myOldId = 0;
int Parallel::myOldNum = 1;
ParallelType Parallel::myOldMode = ParallelType::STATIC;
bool Parallel::myIsolated = false;

//____ Parallel

Parallel::~Parallel() {
}

Parallel::Parallel() {
  myInitialized = false;
  myFinalized = false;
  myId = 0;
  myMasterId = master;
  myNum = 1;
  myAvailableId = 0;
  myAvailableNum = 1;
  myIsParallel = Parallel::isMPI;
  myIAmMaster = true;
  myIAmSlave = true;
  myMode = ParallelType::DYNAMIC;

  myPipeSize = 1;
  myUseBarrier = false;
  myMaxPackages = -1;

  myBuffer = NULL;
  myNext = 0;
  myNextRange[0] = 0;
  myNextRange[1] = 0;
  myDone.clear();
  myBlockList.clear();

  myOldId = 0;
  myOldNum = 1;
  myOldMode = ParallelType::STATIC;
  myIsolated = false;
}

Parallel &Parallel::instance() {
  // We have to do it ourself ... M$ problem ...
  if (obj == NULL) {
    obj = new Parallel();
    atexit(kill);
  }
  return *obj;
}

void Parallel::kill() {
  Parallel *p = obj;
  obj = NULL;
  p->~Parallel();
}

void Parallel::init(int &argc, char ** &argv) {
  instance();
  if (!myInitialized && !myFinalized) {
#ifdef IS_MPI
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &myNum);
    MPI_Comm_rank(MPI_COMM_WORLD, &myId);
    MPI_Barrier(MPI_COMM_WORLD);
    setOpenMDStartSerial(mpiStartSerial);
    setOpenMDEndSerial(mpiEndSerial);
#endif
    setOpenMDAbort(mpiAbort);
    setOpenMDExit(mpiExit);
    myInitialized = true;
    myIsParallel = (myNum > 1);
    myIAmMaster = (myId == myMasterId);
    report.setIAmMaster(iAmMaster());
    setMode(isParallel() ? ParallelType::DYNAMIC : ParallelType::STATIC);
    setPipeSize(myPipeSize);
    TimerStatistic::setParallel(myIsParallel);
    if (iAmMaster() && isParallel())
      myDone.resize(myNum);

    report << plain << (isMPI ? "Using MPI." : "No MPI compilation.") << endr;
  } else
  if (iAmMaster())
    report << recoverable << "MPI is" <<
    ((myInitialized) ? " " : " not ") << "initialized and is" <<
    ((myFinalized) ? " " : " not ") <<
    "finalized. Called [Parallel::init].\n" << endr;
}

void Parallel::finalize() {
  if (ok("Called [Parallel::finalize].")) {
    setOpenMDAbort(NULL);
    setOpenMDExit(NULL);
    setOpendMDStartSerial(NULL);
    setOpenMDEndSerial(NULL);
#ifdef IS_MPI
    if (iAmSlave())
      MPI_Comm_free(&slaveComm);

    if (myBuffer != NULL) {
      int size;
      MPI_Buffer_detach(myBuffer, &size);
      delete[] myBuffer;
      myBuffer = NULL;
    }
    FFTComplex::FFTComplexMPIFinalize();
    MPI_Finalize();
#endif
    myFinalized = true;
  }
}

bool Parallel::ok(const string &err) {
  if (ok())
    return true;

  if (iAmMaster())
    report << recoverable << "MPI is" << ((myInitialized) ? " " : " not ") <<
    "initialized and is" << ((myFinalized) ? " " : " not ") <<
    "finalized. " << err << "\n" << endr;
  return false;
}

void Parallel::setMode(ParallelType mode) {
  if (!ok("Called [Parallel::setMasterSlave]"))
    return;
#ifndef IS_MPI
  mode = ParallelType::STATIC;
#endif
  if (!isParallel() || mode == ParallelType::UNDEFINED)
    mode = ParallelType::STATIC;

  myMode = mode;
  myAvailableNum = (mode == ParallelType::MASTERSLAVE ? myNum - 1  : myNum);
  myIAmSlave = (mode == ParallelType::MASTERSLAVE ? (!myIAmMaster) :
                true);
  myAvailableId =
    (myNum == myAvailableNum ? myId :
     (myId == myMasterId ? -1 : (myId < myMasterId ? myId : myId - 1)));
#ifdef IS_MPI
  if (slaveComm != MPI_COMM_NULL)
    MPI_Comm_free(&slaveComm);

  // Create intracommunicator only with slaves
  MPI_Group worldGroup = MPI_GROUP_NULL;
  MPI_Group slaveGroup = MPI_GROUP_NULL;
  int excl[] = {
    myMasterId
  };

  MPI_Comm_group(MPI_COMM_WORLD, &worldGroup);
  MPI_Group_excl(worldGroup, myNum - myAvailableNum, excl, &slaveGroup);
  MPI_Comm_create(MPI_COMM_WORLD, slaveGroup, &slaveComm);
  MPI_Group_free(&worldGroup);
  MPI_Group_free(&slaveGroup);
  if (slaveComm != MPI_COMM_NULL)
    MPI_Comm_rank(slaveComm, &myAvailableId);
  FFTComplex::FFTComplexMPIInit(
    mode == ParallelType::MASTERSLAVE ? myMasterId : -1);
#endif
#ifdef DEBUG_PARALLEL
  report << allnodesserial << "ParallelMode (" << getId << ") : "
         << getMode().getString() << endr;
#endif
}

void Parallel::setPipeSize(int n) {
  if (!ok("Called [Parallel::setPipeSize]"))
    return;
  myPipeSize = max(n, 0);

#ifdef IS_MPI
  if (myBuffer != NULL) {
    int size;
    MPI_Buffer_detach(myBuffer, &size);
    delete[] myBuffer;
    myBuffer = NULL;
  }

  // Allocate buffer for Bsend()
  const int size = 3 * 2 * myNum * (myPipeSize + 1) + MPI_BSEND_OVERHEAD + 10;
  myBuffer = new int[size];
  for (int i = 0; i < size; i++)
    myBuffer[i] = 0;

  MPI_Buffer_attach(myBuffer, static_cast<int>(size * sizeof(int)));
#endif
}

void Parallel::setMaxPackages(int n) {
  if (!ok("Called [Parallel::setMaxPackages]"))
    return;
  myMaxPackages = n < 0 ? (isDynamic() ? 3 : 0) :
                  n;
}

Parallel::WorkState Parallel::getWorkState() {
  if (!isParallel())
    return SEQUENTIAL;
  else if (myMode == ParallelType::STATIC)
    return STATIC;
  else if (iAmMaster() && myMode == ParallelType::DYNAMIC)
    return MASTER;
  else
    return SLAVE;
}

void Parallel::sync() {
#ifdef IS_MPI
  MPI_Barrier(MPI_COMM_WORLD);
#endif
}

void Parallel::syncSlave() {
#ifdef IS_MPI
  if (!iAmSlave() || !isParallel())
    return;
  MPI_Barrier(slaveComm);
#endif
}

// Senda a vector 3DBlock over MPI as an array.
#ifdef IS_MPI
void Parallel::send(Vector3DBlock *vect, int address) {
  /* Create a C-style array large enough to hold all the real values
   * (3 / Vector3D)
   */
  int size = vect->size();
  Real *vectArray = new Real[3 * size];
  if (vectArray == 0) {
    cout << "Can't create Parallel::send() array! Quitting!" << endl;
    MPI_Abort(MPI_COMM_WORLD, 1);
  }

  /* Start dumping the Vector3D into the array.  This keeps us from having to
   * create an MPI struct, pack it,     * and in general, mess with all that.
   * It's inefficient, but it's at least a start. It's designed such that
   * the array looks like {x1 y1 z1 x2 y2 z2 ... xN yN zN} for N Vector3D's */
  for (int i = 0; i < size; i++) {
    vectArray[3 * i] = (*vect)[i][0];
    vectArray[3 * i + 1] = (*vect)[i][1];
    vectArray[3 * i + 2] = (*vect)[i][2];
  }

  /* Since it's an array of Reals, we can use the plain old Parallel::send
   * routine since MPI can handle
   * both single values and arrays with the same function call */
  send(vectArray, 3 * size, address);
  delete[] vectArray;
}

#else
void Parallel::send(Vector3DBlock *, int) {}

#endif

// Same philosophy as Parallel::send, except we're receiving an array
#ifdef IS_MPI
void Parallel::recv(Vector3DBlock *vect, int address) {
  /* Create an array of proper size. */
  int size = vect->size();
  Real *vectArray = new Real[3 * size];
  if (vectArray == 0) {
    cout << "Can't create Parallel::recv() array!  Quitting!" << endl;
    MPI_Abort(MPI_COMM_WORLD, 1);
  }
  /* Since it's an array, the Parallel::recv() call is sufficient (don't you
   * just love function overloading?) */
  recv(vectArray, 3 * size, address);
  /* Map the vector back onto an actual Vector3DBlock.  The array looks like
   * this: {x1 y1 z1 x2 y2 z2 ... xN yN zN} */
  for (int i = 0; i < size; i++) {
    (*vect)[i][0] = vectArray[3 * i];
    (*vect)[i][1] = vectArray[3 * i + 1];
    (*vect)[i][2] = vectArray[3 * i + 2];
  }

  delete[] vectArray;
}

#else
void Parallel::recv(Vector3DBlock *, int) {}

#endif

// Overwrites the Vector3D with a new one after sending it to another node
#ifdef IS_MPI
void Parallel::sendrecv_replace(Vector3DBlock *vect, int sendaddr,
                                int recvaddr) {
  int size = vect->size();
  Real *vectArray = new Real[3 * size];
  if (vectArray == 0) {
    cout << "Can't create Parallel::send() array! Quitting!" << endl;
    MPI_Abort(MPI_COMM_WORLD, 1);
  }
  // Familiar mapping strategy...
  for (int i = 0; i < size; i++) {
    vectArray[3 * i] = (*vect)[i][0];
    vectArray[3 * i + 1] = (*vect)[i][1];
    vectArray[3 * i + 2] = (*vect)[i][2];
  }

  // Reuse of function calls to actually handle the MPI calls
  sendrecv_replace(vectArray, 3 * size, sendaddr, recvaddr);
  // Map it back in to the Vector3D and nobody will ever know we mucked with
  // it... 8-)
  for (int i = 0; i < size; i++) {
    (*vect)[i][0] = vectArray[3 * i];
    (*vect)[i][1] = vectArray[3 * i + 1];
    (*vect)[i][2] = vectArray[3 * i + 2];
  }

  delete[] vectArray;
}

#else
void Parallel::sendrecv_replace(Vector3DBlock *, int, int) {}

#endif

#ifdef IS_MPI
void Parallel::send(Real *data, int num, int address) {
  // Just a nice wrapper that automatically selects the MPI datatype for you
  // and handles all the annoying things
  MPI_Send(data, num, MPITypeTraits<Real>::datatype, address, 0,
           MPI_COMM_WORLD);
}

#else
void Parallel::send(Real *, int, int) {}

#endif

#ifdef IS_MPI
void Parallel::recv(Real *data, int num, int address) {
  // Another MPI wrapper function....
  MPI_Status status;
  MPI_Recv(data, num, MPITypeTraits<Real>::datatype, address, 0,
           MPI_COMM_WORLD,
           &status);
}

#else
void Parallel::recv(Real *, int, int) {}

#endif

#ifdef IS_MPI
void Parallel::sendrecv(Real *senddata, int sendnum, int sendaddr,
                        Real *recvdata, int recvnum,
                        int recvaddr) {
  MPI_Status status;
  MPI_Sendrecv(senddata, sendnum, MPITypeTraits<Real>::datatype, sendaddr, 0,
               recvdata, recvnum, MPITypeTraits<Real>::datatype, recvaddr, 0,
               MPI_COMM_WORLD,
               &status);
}

#else
void Parallel::sendrecv(Real *, int, int, Real *, int, int) {}

#endif

#ifdef IS_MPI
void Parallel::sendrecv_replace(Real *data, int num, int sendaddr,
                                int recvaddr) {
  MPI_Status status;
  MPI_Sendrecv_replace(data, num, MPITypeTraits<Real>::datatype, sendaddr, 0,
                       recvaddr, 0, MPI_COMM_WORLD,
                       &status);
}

#else
void Parallel::sendrecv_replace(Real *, int, int, int) {}

#endif

#ifdef IS_MPI
void Parallel::gather(Real *data, int num, Real *data_array, int address) {
  MPI_Gather(data, num, MPITypeTraits<Real>::datatype, data_array, num,
             MPITypeTraits<Real>::datatype, address,
             MPI_COMM_WORLD);
}

#else
void Parallel::gather(Real *, int, Real *, int) {}

#endif

#ifdef IS_MPI
void Parallel::allgather(Real *data, int num, Real *data_array) {
  MPI_Allgather(data, num, MPITypeTraits<Real>::datatype, data_array, num,
                MPITypeTraits<Real>::datatype,
                MPI_COMM_WORLD);
}

#else
void Parallel::allgather(Real *, int, Real *) {}

#endif

#ifdef IS_MPI
void Parallel::reduceSlaves(Real *begin, Real *end) {
  if (!iAmSlave() || !isParallel())
    return;
  TimerStatistic::timer[TimerStatistic::COMMUNICATION].start();

  allReduce<true, true>(begin, end);

  TimerStatistic::timer[TimerStatistic::COMMUNICATION].stop();
}

#else
void Parallel::reduceSlaves(Real *, Real *) {}

#endif

#ifdef IS_MPI
void Parallel::distribute(ScalarStructure *energies, Vector3DBlock *coords) {
  if (!iAmMaster() && !energies->distributed()) {
    energies->clear();
    coords->zero();
  }-
  energies->distribute();
  coords->distribute();
}

#else
void Parallel::distribute(ScalarStructure *, Vector3DBlock *) {}

#endif

#ifdef IS_MPI
void Parallel::reduce(ScalarStructure *energies, Vector3DBlock *coords) {
  energies->reduce();
  coords->reduce();
  if (!isParallel() || energies->distributed())
    return;
  TimerStatistic::timer[TimerStatistic::COMMUNICATION].start();

  allReduce<false, true>(coords);
  allReduce<false, false>(energies);

  TimerStatistic::timer[TimerStatistic::COMMUNICATION].stop();
}

#else
void Parallel::reduce(ScalarStructure *, Vector3DBlock *) {}

#endif

#ifdef IS_MPI
void Parallel::bcast(Vector3DBlock *coords) {
  if (!isParallel())
    return;
  broadcast<false, true>(coords);
}

#else
void Parallel::bcast(Vector3DBlock *) {}

#endif

#ifdef IS_MPI
void Parallel::bcast(int &n) {
  if (!isParallel())
    return;
  broadcastScalar<false, true>(n);
}

#else
void Parallel::bcast(int &) {}

#endif

#ifdef IS_MPI
void Parallel::bcastSlaves(Real *begin, Real *end) {
  if (!iAmSlave() || !isParallel())
    return;
  TimerStatistic::timer[TimerStatistic::COMMUNICATION].start();

  broadcast<true, true>(begin, end);

  TimerStatistic::timer[TimerStatistic::COMMUNICATION].stop();
}

#else
void Parallel::bcastSlaves(Real *, Real *) {}

#endif

unsigned int Parallel::getNumberOfPackages(unsigned int n) {
  if (getMaxPackages() < 1 ||
      static_cast<unsigned int>(getAvailableNum()) >= n)
    return n;
  return min(n / getAvailableNum(),
             static_cast<unsigned int>(getMaxPackages()))
         * static_cast<unsigned int>(getAvailableNum());
}

void Parallel::resetNext(const vector<int> &blocks) {
  resetNext();
  if (!iAmMaster())
    return;
#ifdef IS_MPI
  if (blocks.empty())
    return;

  // Vector of ranges (n0,n1,n3, ... ,nM, -1, -1, ..., -1)
  // n0,n1,n3, ... ,nM are the numbers a slave will call next()
  // for one given force
  // A slave will call next() and check if the
  // actual number of next()-calls is inside the
  // range received by the master and compute if true.
  int n = 0;
  for (unsigned int i = 0; i < blocks.size(); i++)
    n += blocks[i];

  if (n < 1)
    return;

  myBlockList.resize(n + 1 + myAvailableNum + 1);
  for (unsigned int i = 0; static_cast<int>(i) <= n; i++)
    myBlockList[i] = i;

  // Adding ranges (-1,-1) to indicate that there is no more
  // to compute.
  for (unsigned int i = n + 1; i < myBlockList.size(); i++)
    myBlockList[i] = -1;

#ifdef DEBUG_PARALLEL
  report << allnodes << plain << "Blocks :";
  for (unsigned int i = 0; i < myBlockList.size(); i++)
    report << myBlockList[i] << " ";

  report << endr;
#endif
  // Fill up the pipe ...
  myI = 0;          // Actual range index
  myRecv = 1;     // Number of pending recieves
  // One more since one slave will send a request but not call
  // recv for the last range
  int stop = 0;     //
  myP = 0;
  myDone[myId] = 0;
  for (int j = 0; j < (myPipeSize + 1) * myNum; j++) {
    if (stop == myNum - 1)
      break;

    int myP = j % myNum;
    if (myP == myMasterId)
      continue;

    MPI_Bsend(&(myBlockList[myI]), 2, MPI_INT, myP, SEND_RANGE,
              MPI_COMM_WORLD);
#ifdef DEBUG_PARALLEL
    report << allnodes << plain << "Block[" << myI << "] = [" <<
    myBlockList[myI] << "," <<
    myBlockList[myI + 1] << "] to " << myP << "." << endr;
#endif
    myDone[myP] = (myBlockList[myI + 2] < 0 ? 1 : 0);
    if (myDone[myP] != 0)
      stop++;
    if (myDone[myP] == 0)
      myRecv++;
    myI++;
  }

  if (myMode == ParallelType::MASTERSLAVE)
    nextMaster();
#endif
}

void Parallel::nextMaster() {
#ifdef IS_MPI
#ifdef DEBUG_PARALLEL
  report << allnodes << plain << "Parallel::nextMaster Recv " << myRecv <<
  "." << endr;
#endif
  while (myRecv > 0) {
    int test = 0;
    MPI_Status status;
    char tmp[1];

    // Receiving the request for a new range from a slave.
    // We wait until we got a msg. Note that waiting for ANY_SOURCE
    // could lead to a starving of some nodes.
    while (!test) {
      myP = (1 + myP) % myNum;
      if (myP != myMasterId)
        MPI_Iprobe(myP, NEED_RANGE, MPI_COMM_WORLD, &test, &status);
    }

    MPI_Recv(tmp, 0, MPI_CHAR, myP, NEED_RANGE, MPI_COMM_WORLD, &status);
    myRecv--;
#ifdef DEBUG_PARALLEL
    report << allnodes << plain << "Recieve from " << myP << "." << endr;
    report << allnodes << plain << "Recv " << myRecv << "." << endr;
#endif

    // We skip to send a further range if the node got the last
    // range since the that node will terminate without
    // receiving any messages.
    if (myDone[myP] == 0) {
      MPI_Bsend(&(myBlockList[myI]), 2, MPI_INT, myP, SEND_RANGE,
                MPI_COMM_WORLD);
#ifdef DEBUG_PARALLEL
      report << allnodes << plain << "Block[" << myI << "] = [" <<
      myBlockList[myI] << "," <<
      myBlockList[myI + 1] << "] to " << myP << "." << endr;
#endif
      myDone[myP] = (myBlockList[myI + 2] < 0 ? 1 : 0);
      if (myDone[myP] == 0)
        myRecv++;
      myI++;
    }
  }
#endif

}

bool Parallel::next() {
  bool doNext = true;
  
#ifdef DEBUG_PARALLEL
  int oldNext = myNext;
#endif
  switch (myWorkState) {
#ifdef IS_MPI
    case MASTER:
      
    do
    {
      
    } while (condition);Next = false;
    do {
      for (int k = 0; k < myNum; ++k) {
        myP = (1 + myP) % myNum;
        if (myP == myMasterId)
          continue;
        
        int test = 0;
        MPI_Status status;
        char tmp[1];
        
        MPI_Iprobe(myP, NEED_RANGE, MPI_COMM_WORLD, &test, &status);
        
        if (test) {
          MPI_Recv(tmp, 0, MPI_CHAR, myP, NEED_RANGE, MPI_COMM_WORLD, &status);
          myRecv--;
#ifdef DEBUG_PARALLEL
          report << allnodes << plain << "Recieve from " << myP << "." <<
              endr;
#endif
          
          // We skip to send a further range if the node got the last
          // range since the that node will terminate without
          // receiving any messages.
          if (myDone[myP] == 0) {
            MPI_Bsend(&(myBlockList[myI]), 2, MPI_INT, myP, SEND_RANGE,
                      MPI_COMM_WORLD);
#if defined (DEBUG_PARALLEL)
            report << allnodes << plain << "Block[" << myI << "] = [" <<
                myBlockList[myI] << "," <<
                myBlockList[myI + 1] << "] to " << myP << "." << endr;
#endif
            myDone[myP] = (myBlockList[myI + 2] < 0 ? 1 : 0);
            if (myDone[myP] == 0)
              myRecv++;
            myI++;
          }
        }
      }
    } while (myBlockList[myI + 2] < 0 && myBlockList[myI + 0] >= 0 &&
             myBlockList[myI + 1] >= 0);

    if (myNextRange[1] == myNext) {
      myNextRange[0] = myBlockList[myI];
      myNextRange[1] = myBlockList[myI + 1];
      myDone[myId] = (myBlockList[myI + 2] < 0 ? 1 : 0);
#if defined (DEBUG_PARALLEL)
      report << allnodes << plain << "Block[" << myI << "] = [" <<
      myBlockList[myI] << "," <<
      myBlockList[myI + 1] << "] to " << getId() << "." << endr;
#endif
      myI++;
    }
    doNext = !(myNext < myNextRange[0] || myNextRange[1] < 0);
    myNext++;

#ifdef DEBUG_PARALLEL
    report << allnodes << plain << "Recv " << myRecv << "." << endr;
#endif
    if (myDone[myId] != 0)
      nextMaster();

    break;

  case Parallel::SLAVE:
    if (myNextRange[1] == myNext) {
      // We need a new range from the master
      char tmp[1];
      MPI_Status status;
      TimerStatistic::timer[TimerStatistic::IDLE].start();
#ifdef DEBUG_PARALLEL
      report << allnodes << plain << "Recv new block " << getId() << "." <<
      endr;
#endif
      MPI_Recv(myNextRange, 2, MPI_INT, myMasterId, SEND_RANGE,
               MPI_COMM_WORLD,
               &status);
#ifdef DEBUG_PARALLEL
      report << allnodes << plain << "Recv new block " << getId() << " [" <<
      myNextRange[0] << "," << myNextRange[1] << "]." << endr;
#endif
      // Asking for the next range such we have it when we need it next time.
      TimerStatistic::timer[TimerStatistic::IDLE].stop();
      if (myNextRange[1] >= 0) {
        MPI_Bsend(tmp, 0, MPI_CHAR, myMasterId, NEED_RANGE, MPI_COMM_WORLD);
#ifdef DEBUG_PARALLEL
        report << allnodes << plain << "Ask new block " << getId() << "." <<
        endr;
#endif
      }
      if (myNextRange[1] < 0)
        myWorkState = DONE;
    }

    doNext = !(myNext < myNextRange[0] || myNextRange[1] < 0);
    myNext++;

    break;
#endif

  case DONE:
    doNext = false;
    break;

  case STATIC:
    doNext = (myNext == myId);
    myNext = (myNext + 1) % myNum;
    break;

  case SEQUENTIAL:
  default:
    break;
  }

#ifdef DEBUG_PARALLEL
  report << allnodes << plain << "Next (" << getId() << ") : "
         << (bool)doNext << ", " << myWorkState << ", " << oldNext << ", ["
         << myNextRange[0] << "," << myNextRange[1] << "]" << endr;
#endif
  return doNext;
}

void Parallel::isolateNode() {
  if (myIsolated == false) {
    myOldId = myId;
    myOldNum = myNum;
    myOldMode = myMode;
    myIsolated = true;

    myId = 0;
    myNum = 1;
    myIsParallel = false;
    myIAmMaster = (myId == myMasterId);
  }
}

void Parallel::integrateNode() {
  if (myIsolated == true) {
    myId = myOldId;
    myNum = myOldNum;
    myIsolated = false;

    myIsParallel = (myNum > 1);
    myIAmMaster = (myId == myMasterId);
  }
}










