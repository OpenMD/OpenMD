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

#include "io/RestReader.hpp"

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <memory>
#include <sstream>
#include <string>

#include <sys/stat.h>
#include <sys/types.h>

#ifdef IS_MPI
#include <mpi.h>
#endif

#include "primitives/Molecule.hpp"
#include "restraints/MolecularRestraint.hpp"
#include "restraints/ObjectRestraint.hpp"
#include "utils/StringTokenizer.hpp"
#include "utils/simError.h"

namespace OpenMD {

  void RestReader::scanFile() {
    std::streampos prevPos;
    std::streampos currPos;

#ifdef IS_MPI

    if (worldRank == 0) {
#endif  // is_mpi

      inFile_->clear();
      currPos = inFile_->tellg();
      prevPos = currPos;

      bool foundOpenSnapshotTag = false;
      int lineNo                = 0;
      while (!foundOpenSnapshotTag && inFile_->getline(buffer, bufferSize)) {
        ++lineNo;

        std::string line = buffer;
        currPos          = inFile_->tellg();
        if (line.find("<Snapshot>") != std::string::npos) {
          foundOpenSnapshotTag = true;
          framePos_            = (long long)prevPos;
        }
        prevPos = currPos;
      }

#ifdef IS_MPI
    }
    MPI_Bcast(&framePos_, 1, MPI_LONG_LONG, 0, MPI_COMM_WORLD);
#endif  // is_mpi
  }

  void RestReader::readSet() {
    std::string line;

#ifndef IS_MPI

    inFile_->clear();
    inFile_->seekg(framePos_);

    std::istream& inputStream = *inFile_;
#else

    int primaryNode = 0;
    std::stringstream sstream;
    if (worldRank == primaryNode) {
      std::string sendBuffer;

      inFile_->clear();
      inFile_->seekg(framePos_);

      while (inFile_->getline(buffer, bufferSize)) {
        line = buffer;
        sendBuffer += line;
        sendBuffer += '\n';
        if (line.find("</Snapshot>") != std::string::npos) { break; }
      }

      int sendBufferSize = sendBuffer.size();
      MPI_Bcast(&sendBufferSize, 1, MPI_INT, primaryNode, MPI_COMM_WORLD);
      MPI_Bcast((void*)sendBuffer.c_str(), sendBufferSize, MPI_CHAR,
                primaryNode, MPI_COMM_WORLD);

      sstream.str(sendBuffer);
    } else {
      int sendBufferSize;
      MPI_Bcast(&sendBufferSize, 1, MPI_INT, primaryNode, MPI_COMM_WORLD);
      char* recvBuffer = new char[sendBufferSize + 1];
      assert(recvBuffer);
      recvBuffer[sendBufferSize] = '\0';
      MPI_Bcast(recvBuffer, sendBufferSize, MPI_CHAR, primaryNode,
                MPI_COMM_WORLD);
      sstream.str(recvBuffer);
      delete[] recvBuffer;
    }

    std::istream& inputStream = sstream;
#endif

    inputStream.getline(buffer, bufferSize);

    line = buffer;
    if (line.find("<Snapshot>") == std::string::npos) {
      snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
               "RestReader Error: can not find <Snapshot>\n");
      painCave.isFatal = 1;
      simError();
    }

    // read frameData
    readFrameProperties(inputStream);

    // read StuntDoubles
    readStuntDoubles(inputStream);

    inputStream.getline(buffer, bufferSize);
    line = buffer;
    if (line.find("</Snapshot>") == std::string::npos) {
      snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
               "RestReader Error: can not find </Snapshot>\n");
      painCave.isFatal = 1;
      simError();
    }
  }

  void RestReader::readReferenceStructure() {
    // We need temporary storage to keep track of all StuntDouble positions
    // in case some of the restraints are molecular (i.e. if they use
    // multiple SD positions to determine restrained orientations or positions:

    all_pos_.clear();
    all_pos_.resize(info_->getNGlobalIntegrableObjects());

    // Restraint files are just standard dump files, but with the reference
    // structure stored in the first frame (frame 0).
    // RestReader overloads readSet and explicitly handles all of the
    // ObjectRestraints in that method:

    scanFile();

    readSet();

    // all ObjectRestraints have been handled, now we have to worry about
    // molecular restraints:

    SimInfo::MoleculeIterator i;
    Molecule::IntegrableObjectIterator j;
    Molecule* mol;
    StuntDouble* sd;

    // no need to worry about parallel molecules, as molecules are not
    // split across processor boundaries.  Just loop over all molecules
    // we know about:

    for (mol = info_->beginMolecule(i); mol != NULL;
         mol = info_->nextMolecule(i)) {
      // is this molecule restrained?
      std::shared_ptr<GenericData> data = mol->getPropertyByName("Restraint");

      if (data != nullptr) {
        // make sure we can reinterpret the generic data as restraint data:

        std::shared_ptr<RestraintData> restData =
            std::dynamic_pointer_cast<RestraintData>(data);

        if (restData != nullptr) {
          // make sure we can reinterpet the restraint data as a
          // pointer to a MolecularRestraint:

          MolecularRestraint* mRest =
              dynamic_cast<MolecularRestraint*>(restData->getData());

          if (mRest == NULL) {
            snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
                     "Can not cast RestraintData to MolecularRestraint\n");
            painCave.severity = OPENMD_ERROR;
            painCave.isFatal  = 1;
            simError();
          } else {
            // now we need to pack the stunt doubles for the reference
            // structure:

            std::vector<Vector3d> ref;
            int count = 0;
            RealType mass, mTot;
            Vector3d COM(0.0);

            mTot = 0.0;
            // loop over the stunt doubles in this molecule in the order we
            // will be looping them in the restraint code:

            for (sd = mol->beginIntegrableObject(j); sd != NULL;
                 sd = mol->nextIntegrableObject(j)) {
              // push back the reference positions of the stunt
              // doubles from the *globally* sorted array of
              // positions:

              ref.push_back(all_pos_[sd->getGlobalIntegrableObjectIndex()]);
              mass  = sd->getMass();
              COM   = COM + mass * ref[count];
              mTot  = mTot + mass;
              count = count + 1;
            }
            COM /= mTot;
            mRest->setReferenceStructure(ref, COM);
          }
        }
      }
    }
  }

  void RestReader::parseDumpLine(const std::string& line) {
    StringTokenizer tokenizer(line);
    int nTokens;

    nTokens = tokenizer.countTokens();

    if (nTokens < 2) {
      snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
               "RestReader Error: Not enough Tokens.\n%s\n", line.c_str());
      painCave.isFatal = 1;
      simError();
    }

    int index = tokenizer.nextTokenAsInt();

    StuntDouble* sd = info_->getIOIndexToIntegrableObject(index);

    if (sd == NULL) { return; }

    std::string type = tokenizer.nextToken();
    int size         = type.size();

    Vector3d pos;
    Quat4d q;

    for (int i = 0; i < size; ++i) {
      switch (type[i]) {
      case 'p': {
        pos[0] = tokenizer.nextTokenAsDouble();
        pos[1] = tokenizer.nextTokenAsDouble();
        pos[2] = tokenizer.nextTokenAsDouble();
        break;
      }
      case 'v': {
        Vector3d vel;
        vel[0] = tokenizer.nextTokenAsDouble();
        vel[1] = tokenizer.nextTokenAsDouble();
        vel[2] = tokenizer.nextTokenAsDouble();
        break;
      }

      case 'q': {
        if (sd->isDirectional()) {
          q[0] = tokenizer.nextTokenAsDouble();
          q[1] = tokenizer.nextTokenAsDouble();
          q[2] = tokenizer.nextTokenAsDouble();
          q[3] = tokenizer.nextTokenAsDouble();

          RealType qlen = q.length();
          if (qlen < OpenMD::epsilon) {  // check quaternion is not equal to 0

            snprintf(
                painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
                "RestReader Error: initial quaternion error (q0^2 + q1^2 + "
                "q2^2 + "
                "q3^2) ~ 0\n");
            painCave.isFatal = 1;
            simError();
          }

          q.normalize();
        }
        break;
      }
      case 'j': {
        Vector3d ji;
        if (sd->isDirectional()) {
          ji[0] = tokenizer.nextTokenAsDouble();
          ji[1] = tokenizer.nextTokenAsDouble();
          ji[2] = tokenizer.nextTokenAsDouble();
        }
        break;
      }
      case 'f': {
        Vector3d force;
        force[0] = tokenizer.nextTokenAsDouble();
        force[1] = tokenizer.nextTokenAsDouble();
        force[2] = tokenizer.nextTokenAsDouble();
        break;
      }
      case 't': {
        Vector3d torque;
        torque[0] = tokenizer.nextTokenAsDouble();
        torque[1] = tokenizer.nextTokenAsDouble();
        torque[2] = tokenizer.nextTokenAsDouble();
        break;
      }
      default: {
        snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
                 "RestReader Error: %s is an unrecognized type\n",
                 type.c_str());
        painCave.isFatal = 1;
        simError();
        break;
      }
      }
      // keep the position in case we need it for a molecular restraint:

      all_pos_[index] = pos;

      // is this io restrained?
      std::shared_ptr<GenericData> data = sd->getPropertyByName("Restraint");

      if (data != nullptr) {
        // make sure we can reinterpret the generic data as restraint data:
        std::shared_ptr<RestraintData> restData =
            std::dynamic_pointer_cast<RestraintData>(data);
        if (restData != nullptr) {
          // make sure we can reinterpet the restraint data as a pointer to
          // an ObjectRestraint:
          ObjectRestraint* oRest =
              dynamic_cast<ObjectRestraint*>(restData->getData());
          if (oRest == NULL) {
            snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
                     "Can not cast RestraintData to ObjectRestraint\n");
            painCave.severity = OPENMD_ERROR;
            painCave.isFatal  = 1;
            simError();
          } else {
            if (sd->isDirectional()) {
              oRest->setReferenceStructure(pos, q.toRotationMatrix3());
            } else {
              oRest->setReferenceStructure(pos);
            }
          }
        }
      }
    }
  }

  void RestReader::readStuntDoubles(std::istream& inputStream) {
    inputStream.getline(buffer, bufferSize);
    std::string line(buffer);

    if (line.find("<StuntDoubles>") == std::string::npos) {
      snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
               "RestReader Error: Missing <StuntDoubles>\n");
      painCave.isFatal = 1;
      simError();
    }

    while (inputStream.getline(buffer, bufferSize)) {
      line = buffer;

      if (line.find("</StuntDoubles>") != std::string::npos) { break; }

      parseDumpLine(line);
    }
  }

  void RestReader::readFrameProperties(std::istream& inputStream) {
    inputStream.getline(buffer, bufferSize);
    std::string line(buffer);

    if (line.find("<FrameData>") == std::string::npos) {
      snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
               "RestReader Error: Missing <FrameData>\n");
      painCave.isFatal = 1;
      simError();
    }

    // restraints don't care about frame data (unless we need to wrap
    // coordinates, but we'll worry about that later), so
    // we'll just scan ahead until the end of the frame data:

    while (inputStream.getline(buffer, bufferSize)) {
      line = buffer;

      if (line.find("</FrameData>") != std::string::npos) { break; }
    }
  }

}  // namespace OpenMD
