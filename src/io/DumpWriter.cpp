/*
 * Copyright (c) 2005 The University of Notre Dame. All Rights Reserved.
 *
 * The University of Notre Dame grants you ("Licensee") a
 * non-exclusive, royalty free, license to use, modify and
 * redistribute this software in source and binary code form, provided
 * that the following conditions are met:
 *
 * 1. Acknowledgement of the program authors must be made in any
 *    publication of scientific results based in part on use of the
 *    program.  An acceptable form of acknowledgement is citation of
 *    the article in which the program was described (Matthew
 *    A. Meineke, Charles F. Vardeman II, Teng Lin, Christopher
 *    J. Fennell and J. Daniel Gezelter, "OOPSE: An Object-Oriented
 *    Parallel Simulation Engine for Molecular Dynamics,"
 *    J. Comput. Chem. 26, pp. 252-271 (2005))
 *
 * 2. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 *
 * 3. Redistributions in binary form must reproduce the above copyright
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
 */
 
#include "io/DumpWriter.hpp"
#include "primitives/Molecule.hpp"
#include "utils/simError.h"
#include "io/basic_teebuf.hpp"
#include "io/gzstream.hpp"
#include "io/Globals.hpp"

#ifdef IS_MPI
#include <mpi.h>
#endif //is_mpi

namespace oopse {

  DumpWriter::DumpWriter(SimInfo* info) 
    : info_(info), filename_(info->getDumpFileName()), eorFilename_(info->getFinalConfigFileName()){

    Globals* simParams = info->getSimParams();
    needCompression_ = simParams->getCompressDumpFile();
    needForceVector_ = simParams->getDumpForceVector();

#ifdef HAVE_LIBZ
    if (needCompression_) {
        filename_ += ".gz";
        eorFilename_ += ".gz";
    }
#endif
    
#ifdef IS_MPI

      if (worldRank == 0) {
#endif // is_mpi


        dumpFile_ = createOStream(filename_);

        if (!dumpFile_) {
	  sprintf(painCave.errMsg, "Could not open \"%s\" for dump output.\n",
		  filename_.c_str());
	  painCave.isFatal = 1;
	  simError();
        }

#ifdef IS_MPI

      }

      sprintf(checkPointMsg, "Sucessfully opened output file for dumping.\n");
      MPIcheckPoint();

#endif // is_mpi

    }


  DumpWriter::DumpWriter(SimInfo* info, const std::string& filename) 
    : info_(info), filename_(filename){

    Globals* simParams = info->getSimParams();
    eorFilename_ = filename_.substr(0, filename_.rfind(".")) + ".eor";    

    needCompression_ = simParams->getCompressDumpFile();
    needForceVector_ = simParams->getDumpForceVector();

#ifdef HAVE_LIBZ
    if (needCompression_) {
        filename_ += ".gz";
        eorFilename_ += ".gz";
    }
#endif
    
#ifdef IS_MPI

      if (worldRank == 0) {
#endif // is_mpi


        dumpFile_ = createOStream(filename_);

        if (!dumpFile_) {
	  sprintf(painCave.errMsg, "Could not open \"%s\" for dump output.\n",
		  filename_.c_str());
	  painCave.isFatal = 1;
	  simError();
        }

#ifdef IS_MPI

      }

      sprintf(checkPointMsg, "Sucessfully opened output file for dumping.\n");
      MPIcheckPoint();

#endif // is_mpi

    }

  DumpWriter::~DumpWriter() {

#ifdef IS_MPI

    if (worldRank == 0) {
#endif // is_mpi

      delete dumpFile_;

#ifdef IS_MPI

    }

#endif // is_mpi

  }

  void DumpWriter::writeCommentLine(std::ostream& os, Snapshot* s) {

    double currentTime;
    Mat3x3d hmat;
    double chi;
    double integralOfChiDt;
    Mat3x3d eta;
    
    currentTime = s->getTime();
    hmat = s->getHmat();
    chi = s->getChi();
    integralOfChiDt = s->getIntegralOfChiDt();
    eta = s->getEta();
    
    os << currentTime << ";\t" 
       << hmat(0, 0) << "\t" << hmat(1, 0) << "\t" << hmat(2, 0) << ";\t" 
       << hmat(0, 1) << "\t" << hmat(1, 1) << "\t" << hmat(2, 1) << ";\t"
       << hmat(0, 2) << "\t" << hmat(1, 2) << "\t" << hmat(2, 2) << ";\t";

    //write out additional parameters, such as chi and eta

    os << chi << "\t" << integralOfChiDt << "\t;";

    os << eta(0, 0) << "\t" << eta(1, 0) << "\t" << eta(2, 0) << ";\t" 
       << eta(0, 1) << "\t" << eta(1, 1) << "\t" << eta(2, 1) << ";\t"
       << eta(0, 2) << "\t" << eta(1, 2) << "\t" << eta(2, 2) << ";";
        
    os << "\n";
  }

  void DumpWriter::writeFrame(std::ostream& os) {
    const int BUFFERSIZE = 2000;
    const int MINIBUFFERSIZE = 100;

    char tempBuffer[BUFFERSIZE];
    char writeLine[BUFFERSIZE];

    Quat4d q;
    Vector3d ji;
    Vector3d pos;
    Vector3d vel;
    Vector3d frc;
    Vector3d trq;

    Molecule* mol;
    StuntDouble* integrableObject;
    SimInfo::MoleculeIterator mi;
    Molecule::IntegrableObjectIterator ii;
  
    int nTotObjects;    
    nTotObjects = info_->getNGlobalIntegrableObjects();

#ifndef IS_MPI


    os << nTotObjects << "\n";
        
    writeCommentLine(os, info_->getSnapshotManager()->getCurrentSnapshot());

    for (mol = info_->beginMolecule(mi); mol != NULL; mol = info_->nextMolecule(mi)) {

      for (integrableObject = mol->beginIntegrableObject(ii); integrableObject != NULL; 
	   integrableObject = mol->nextIntegrableObject(ii)) { 
                

	pos = integrableObject->getPos();
	vel = integrableObject->getVel();

	sprintf(tempBuffer, "%s\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t",
		integrableObject->getType().c_str(), 
		pos[0], pos[1], pos[2],
		vel[0], vel[1], vel[2]);

	strcpy(writeLine, tempBuffer);

	if (integrableObject->isDirectional()) {
	  q = integrableObject->getQ();
	  ji = integrableObject->getJ();

	  sprintf(tempBuffer, "%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf", 
		  q[0], q[1], q[2], q[3],
		  ji[0], ji[1], ji[2]);
	  strcat(writeLine, tempBuffer);
	} else {
	  strcat(writeLine, "0.0\t0.0\t0.0\t0.0\t0.0\t0.0\t0.0");
	}

	if (needForceVector_) {
	  frc = integrableObject->getFrc();
	  trq = integrableObject->getTrq();
	  
	  sprintf(tempBuffer, "\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf", 
		  frc[0], frc[1], frc[2],
		  trq[0], trq[1], trq[2]);
	  strcat(writeLine, tempBuffer);
	}
	
	strcat(writeLine, "\n");
	os << writeLine;

      }
    }

    os.flush();
#else // is_mpi
    /*********************************************************************
     * Documentation?  You want DOCUMENTATION?
     * 
     * Why all the potatoes below?  
     *
     * To make a long story short, the original version of DumpWriter
     * worked in the most inefficient way possible.  Node 0 would 
     * poke each of the node for an individual atom's formatted data 
     * as node 0 worked its way down the global index. This was particularly 
     * inefficient since the method blocked all processors at every atom 
     * (and did it twice!).
     *
     * An intermediate version of DumpWriter could be described from Node
     * zero's perspective as follows:
     * 
     *  1) Have 100 of your friends stand in a circle.
     *  2) When you say go, have all of them start tossing potatoes at
     *     you (one at a time).
     *  3) Catch the potatoes.
     *
     * It was an improvement, but MPI has buffers and caches that could 
     * best be described in this analogy as "potato nets", so there's no 
     * need to block the processors atom-by-atom.
     * 
     * This new and improved DumpWriter works in an even more efficient 
     * way:
     * 
     *  1) Have 100 of your friend stand in a circle.
     *  2) When you say go, have them start tossing 5-pound bags of 
     *     potatoes at you.
     *  3) Once you've caught a friend's bag of potatoes,
     *     toss them a spud to let them know they can toss another bag.
     *
     * How's THAT for documentation?
     *
     *********************************************************************/
    const int masterNode = 0;

    int * potatoes;
    int myPotato;
    int nProc;
    int which_node;
    double atomData[19];
    int isDirectional;
    char MPIatomTypeString[MINIBUFFERSIZE];
    int msgLen; // the length of message actually recieved at master nodes
    int haveError;
    MPI_Status istatus;
    int nCurObj;
    
    // code to find maximum tag value
    int * tagub;
    int flag;
    int MAXTAG;
    MPI_Attr_get(MPI_COMM_WORLD, MPI_TAG_UB, &tagub, &flag);

    if (flag) {
      MAXTAG = *tagub;
    } else {
      MAXTAG = 32767;
    }

    if (worldRank == masterNode) { //master node (node 0) is responsible for writing the dump file

      // Node 0 needs a list of the magic potatoes for each processor;

      MPI_Comm_size(MPI_COMM_WORLD, &nProc);
      potatoes = new int[nProc];

      //write out the comment lines
      for(int i = 0; i < nProc; i++) {
	potatoes[i] = 0;
      }


      os << nTotObjects << "\n";
      writeCommentLine(os, info_->getSnapshotManager()->getCurrentSnapshot());

      for(int i = 0; i < info_->getNGlobalMolecules(); i++) {

	// Get the Node number which has this atom;

	which_node = info_->getMolToProc(i);

	if (which_node != masterNode) { //current molecule is in slave node
	  if (potatoes[which_node] + 1 >= MAXTAG) {
	    // The potato was going to exceed the maximum value, 
	    // so wrap this processor potato back to 0:         

	    potatoes[which_node] = 0;
	    MPI_Send(&potatoes[which_node], 1, MPI_INT, which_node, 0,
		     MPI_COMM_WORLD);
	  }

	  myPotato = potatoes[which_node];

	  //recieve the number of integrableObject in current molecule
	  MPI_Recv(&nCurObj, 1, MPI_INT, which_node, myPotato,
		   MPI_COMM_WORLD, &istatus);
	  myPotato++;

	  for(int l = 0; l < nCurObj; l++) {
	    if (potatoes[which_node] + 2 >= MAXTAG) {
	      // The potato was going to exceed the maximum value, 
	      // so wrap this processor potato back to 0:         

	      potatoes[which_node] = 0;
	      MPI_Send(&potatoes[which_node], 1, MPI_INT, which_node,
		       0, MPI_COMM_WORLD);
	    }

	    MPI_Recv(MPIatomTypeString, MINIBUFFERSIZE, MPI_CHAR,
		     which_node, myPotato, MPI_COMM_WORLD,
		     &istatus);

	    myPotato++;

	    MPI_Recv(atomData, 19, MPI_DOUBLE, which_node, myPotato,
		     MPI_COMM_WORLD, &istatus);
	    myPotato++;

	    MPI_Get_count(&istatus, MPI_DOUBLE, &msgLen);

	    if (msgLen == 13 || msgLen == 19)
	      isDirectional = 1;
	    else
	      isDirectional = 0;

	    // If we've survived to here, format the line:

	    if (!isDirectional) {
	      sprintf(writeLine, "%s\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t",
		      MPIatomTypeString, atomData[0],
		      atomData[1], atomData[2],
		      atomData[3], atomData[4],
		      atomData[5]);

	      strcat(writeLine,
		     "0.0\t0.0\t0.0\t0.0\t0.0\t0.0\t0.0");
	    } else {
	      sprintf(writeLine,
		      "%s\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf",
		      MPIatomTypeString,
		      atomData[0],
		      atomData[1],
		      atomData[2],
		      atomData[3],
		      atomData[4],
		      atomData[5],
		      atomData[6],
		      atomData[7],
		      atomData[8],
		      atomData[9],
		      atomData[10],
		      atomData[11],
		      atomData[12]);
	    }
	    
	    if (needForceVector_) {
	      if (!isDirectional) {
		sprintf(writeLine, "\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf",
			atomData[6], 
			atomData[7], 
			atomData[8],
			atomData[9],
			atomData[10],
			atomData[11]);
	      } else {
		sprintf(writeLine, "\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf",
			atomData[13], 
			atomData[14], 
			atomData[15],
			atomData[16],
			atomData[17],
			atomData[18]);
	      }
	    }

	    sprintf(writeLine, "\n");
	    os << writeLine;

	  } // end for(int l =0)

	  potatoes[which_node] = myPotato;
	} else { //master node has current molecule

	  mol = info_->getMoleculeByGlobalIndex(i);

	  if (mol == NULL) {
	    sprintf(painCave.errMsg, "Molecule not found on node %d!", worldRank);
	    painCave.isFatal = 1;
	    simError();
	  }
                
	  for (integrableObject = mol->beginIntegrableObject(ii); integrableObject != NULL; 
	       integrableObject = mol->nextIntegrableObject(ii)) {      

	    pos = integrableObject->getPos();
	    vel = integrableObject->getVel();

	    atomData[0] = pos[0];
	    atomData[1] = pos[1];
	    atomData[2] = pos[2];

	    atomData[3] = vel[0];
	    atomData[4] = vel[1];
	    atomData[5] = vel[2];

	    isDirectional = 0;

	    if (integrableObject->isDirectional()) {
	      isDirectional = 1;

	      q = integrableObject->getQ();
	      ji = integrableObject->getJ();

	      for(int j = 0; j < 6; j++) {
		atomData[j] = atomData[j];
	      }

	      atomData[6] = q[0];
	      atomData[7] = q[1];
	      atomData[8] = q[2];
	      atomData[9] = q[3];

	      atomData[10] = ji[0];
	      atomData[11] = ji[1];
	      atomData[12] = ji[2];
	    }

	    if (needForceVector_) {
	      frc = integrableObject->getFrc();
	      trq = integrableObject->getTrq();

	      if (!isDirectional) {
		atomData[6] = frc[0];
		atomData[7] = frc[1];
		atomData[8] = frc[2];
		atomData[9] = trq[0];
		atomData[10] = trq[1];
		atomData[11] = trq[2];
	      } else {
		atomData[13] = frc[0];
		atomData[14] = frc[1];
		atomData[15] = frc[2];
		atomData[16] = trq[0];
		atomData[17] = trq[1];
		atomData[18] = trq[2];
	      }
	    }

	    // If we've survived to here, format the line:

	    if (!isDirectional) {
	      sprintf(writeLine, "%s\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t",
		      integrableObject->getType().c_str(), atomData[0],
		      atomData[1], atomData[2],
		      atomData[3], atomData[4],
		      atomData[5]);

	      strcat(writeLine,
		     "0.0\t0.0\t0.0\t0.0\t0.0\t0.0\t0.0");
	    } else {
	      sprintf(writeLine,
		      "%s\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf",
		      integrableObject->getType().c_str(),
		      atomData[0],
		      atomData[1],
		      atomData[2],
		      atomData[3],
		      atomData[4],
		      atomData[5],
		      atomData[6],
		      atomData[7],
		      atomData[8],
		      atomData[9],
		      atomData[10],
		      atomData[11],
		      atomData[12]);
	    }

	    if (needForceVector_) {
	      if (!isDirectional) {
	      sprintf(writeLine, "\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf",
		      atomData[6],
		      atomData[7], 
		      atomData[8],
		      atomData[9],
		      atomData[10],
		      atomData[11]);
	      } else {
		sprintf(writeLine, "\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf",
			atomData[13],
			atomData[14], 
			atomData[15],
			atomData[16],
			atomData[17],
			atomData[18]);
	      }
	    }

	    sprintf(writeLine, "\n");
	    os << writeLine;

	  } //end for(iter = integrableObject.begin())
	}
      } //end for(i = 0; i < mpiSim->getNmol())

      os.flush();
        
      sprintf(checkPointMsg, "Sucessfully took a dump.\n");
      MPIcheckPoint();

      delete [] potatoes;
    } else {

      // worldRank != 0, so I'm a remote node.  

      // Set my magic potato to 0:

      myPotato = 0;

      for(int i = 0; i < info_->getNGlobalMolecules(); i++) {

	// Am I the node which has this integrableObject?
	int whichNode = info_->getMolToProc(i);
	if (whichNode == worldRank) {
	  if (myPotato + 1 >= MAXTAG) {

	    // The potato was going to exceed the maximum value, 
	    // so wrap this processor potato back to 0 (and block until
	    // node 0 says we can go:

	    MPI_Recv(&myPotato, 1, MPI_INT, 0, 0, MPI_COMM_WORLD,
		     &istatus);
	  }

	  mol = info_->getMoleculeByGlobalIndex(i);

                
	  nCurObj = mol->getNIntegrableObjects();

	  MPI_Send(&nCurObj, 1, MPI_INT, 0, myPotato, MPI_COMM_WORLD);
	  myPotato++;

	  for (integrableObject = mol->beginIntegrableObject(ii); integrableObject != NULL; 
	       integrableObject = mol->nextIntegrableObject(ii)) {

	    if (myPotato + 2 >= MAXTAG) {

	      // The potato was going to exceed the maximum value, 
	      // so wrap this processor potato back to 0 (and block until
	      // node 0 says we can go:

	      MPI_Recv(&myPotato, 1, MPI_INT, 0, 0, MPI_COMM_WORLD,
		       &istatus);
	    }

	    pos = integrableObject->getPos();
	    vel = integrableObject->getVel();

	    atomData[0] = pos[0];
	    atomData[1] = pos[1];
	    atomData[2] = pos[2];

	    atomData[3] = vel[0];
	    atomData[4] = vel[1];
	    atomData[5] = vel[2];

	    isDirectional = 0;

	    if (integrableObject->isDirectional()) {
	      isDirectional = 1;

	      q = integrableObject->getQ();
	      ji = integrableObject->getJ();

	      atomData[6] = q[0];
	      atomData[7] = q[1];
	      atomData[8] = q[2];
	      atomData[9] = q[3];

	      atomData[10] = ji[0];
	      atomData[11] = ji[1];
	      atomData[12] = ji[2];
	    }

	    if (needForceVector_) {
	      frc = integrableObject->getFrc();
	      trq = integrableObject->getTrq();
	      
	      if (!isDirectional) {
		atomData[6] = frc[0];
		atomData[7] = frc[1];
		atomData[8] = frc[2];
		
		atomData[9] = trq[0];
		atomData[10] = trq[1];
		atomData[11] = trq[2];
	      } else {
		atomData[13] = frc[0];
		atomData[14] = frc[1];
		atomData[15] = frc[2];
		
		atomData[16] = trq[0];
		atomData[17] = trq[1];
		atomData[18] = trq[2];
	      }
	    }

	    strncpy(MPIatomTypeString, integrableObject->getType().c_str(), MINIBUFFERSIZE);

	    // null terminate the  std::string before sending (just in case):
	    MPIatomTypeString[MINIBUFFERSIZE - 1] = '\0';

	    MPI_Send(MPIatomTypeString, MINIBUFFERSIZE, MPI_CHAR, 0,
		     myPotato, MPI_COMM_WORLD);

	    myPotato++;

	    if (isDirectional && needForceVector_) {
	      MPI_Send(atomData, 19, MPI_DOUBLE, 0, myPotato,
		       MPI_COMM_WORLD);
	    } else if (isDirectional) {
	      MPI_Send(atomData, 13, MPI_DOUBLE, 0, myPotato,
		       MPI_COMM_WORLD);
	    } else if (needForceVector_) {
	      MPI_Send(atomData, 12, MPI_DOUBLE, 0, myPotato,
		       MPI_COMM_WORLD);
	    } else {
	      MPI_Send(atomData, 6, MPI_DOUBLE, 0, myPotato,
		       MPI_COMM_WORLD);
	    }

	    myPotato++;
	  }
                    
	}
            
      }
      sprintf(checkPointMsg, "Sucessfully took a dump.\n");
      MPIcheckPoint();
    }

#endif // is_mpi

  }

  void DumpWriter::writeDump() {
    writeFrame(*dumpFile_);
  }

  void DumpWriter::writeEor() {
    std::ostream* eorStream;
    
#ifdef IS_MPI
    if (worldRank == 0) {
#endif // is_mpi

      eorStream = createOStream(eorFilename_);

#ifdef IS_MPI
    }
#endif // is_mpi    

    writeFrame(*eorStream);

#ifdef IS_MPI
    if (worldRank == 0) {
#endif // is_mpi
    delete eorStream;

#ifdef IS_MPI
    }
#endif // is_mpi  

  }


  void DumpWriter::writeDumpAndEor() {
    std::vector<std::streambuf*> buffers;
    std::ostream* eorStream;
#ifdef IS_MPI
    if (worldRank == 0) {
#endif // is_mpi

      buffers.push_back(dumpFile_->rdbuf());

      eorStream = createOStream(eorFilename_);

      buffers.push_back(eorStream->rdbuf());
        
#ifdef IS_MPI
    }
#endif // is_mpi    

    TeeBuf tbuf(buffers.begin(), buffers.end());
    std::ostream os(&tbuf);

    writeFrame(os);

#ifdef IS_MPI
    if (worldRank == 0) {
#endif // is_mpi
    delete eorStream;

#ifdef IS_MPI
    }
#endif // is_mpi  
    
  }

std::ostream* DumpWriter::createOStream(const std::string& filename) {

    std::ostream* newOStream;
#ifdef HAVE_LIBZ 
    if (needCompression_) {
        newOStream = new ogzstream(filename.c_str());
    } else {
        newOStream = new std::ofstream(filename.c_str());
    }
#else
    newOStream = new std::ofstream(filename.c_str());
#endif
    return newOStream;
}

}//end namespace oopse
