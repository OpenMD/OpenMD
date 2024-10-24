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

#define _LARGEFILE_SOURCE64
#define _FILE_OFFSET_BITS 64

#include "io/DumpReader.hpp"

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <iostream>

#include <sys/stat.h>
#include <sys/types.h>

#ifdef IS_MPI
#include <mpi.h>
#endif

#include "brains/Thermo.hpp"
#include "primitives/Molecule.hpp"
#include "utils/StringTokenizer.hpp"
#include "utils/simError.h"

namespace OpenMD {

  DumpReader::DumpReader(SimInfo* info, const std::string& filename) :
      info_(info), filename_(filename), isScanned_(false), nframes_(0),
      needCOMprops_(false) {
#ifdef IS_MPI
    if (worldRank == 0) {
#endif

      inFile_ =
          std::ifstream(filename_.c_str(), ifstream::in | ifstream::binary);

      if (inFile_.fail()) {
        snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
                 "DumpReader: Cannot open file: %s\n", filename_.c_str());
        painCave.isFatal = 1;
        simError();
      }
#ifdef IS_MPI
    }
    strcpy(checkPointMsg, "Dump file opened for reading successfully.");
    errorCheckPoint();
#endif
  }

  DumpReader::~DumpReader() {
#ifdef IS_MPI
    strcpy(checkPointMsg, "Dump file closed successfully.");
    errorCheckPoint();
#endif
  }

  int DumpReader::getNFrames(void) {
    if (!isScanned_) scanFile();

    return nframes_;
  }

  void DumpReader::scanFile(void) {
    std::streampos prevPos;
    std::streampos currPos;

#ifdef IS_MPI
    if (worldRank == 0) {
#endif

      currPos                     = inFile_.tellg();
      prevPos                     = currPos;
      bool foundOpenSnapshotTag   = false;
      bool foundClosedSnapshotTag = false;

      int lineNo = 0;
      while (inFile_.getline(buffer, bufferSize)) {
        ++lineNo;

        std::string line = buffer;
        currPos          = inFile_.tellg();
        if (line.find("<Snapshot>") != std::string::npos) {
          if (foundOpenSnapshotTag) {
            snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
                     "DumpReader:<Snapshot> is multiply nested at line %d "
                     "in %s \n",
                     lineNo, filename_.c_str());
            painCave.isFatal = 1;
            simError();
          }
          foundOpenSnapshotTag   = true;
          foundClosedSnapshotTag = false;
          framePos_.push_back(prevPos);

        } else if (line.find("</Snapshot>") != std::string::npos) {
          if (!foundOpenSnapshotTag) {
            snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
                     "DumpReader:</Snapshot> appears before <Snapshot> at "
                     "line %d in %s \n",
                     lineNo, filename_.c_str());
            painCave.isFatal = 1;
            simError();
          }

          if (foundClosedSnapshotTag) {
            snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
                     "DumpReader:</Snapshot> appears multiply nested at "
                     "line %d in %s \n",
                     lineNo, filename_.c_str());
            painCave.isFatal = 1;
            simError();
          }
          foundClosedSnapshotTag = true;
          foundOpenSnapshotTag   = false;
        }
        prevPos = currPos;
      }

      // only found <Snapshot> for the last frame means the file is
      // corrupted, we should discard it and give a warning message
      if (foundOpenSnapshotTag) {
        snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
                 "DumpReader: last frame in %s is invalid\n",
                 filename_.c_str());
        painCave.isFatal = 0;
        simError();
        framePos_.pop_back();
      }

      nframes_ = framePos_.size();

      if (nframes_ == 0) {
        snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
                 "DumpReader: %s does not contain a valid frame\n",
                 filename_.c_str());
        painCave.isFatal = 1;
        simError();
      }

#ifdef IS_MPI
    }
    MPI_Bcast(&nframes_, 1, MPI_INT, 0, MPI_COMM_WORLD);
#endif  // is_mpi

    isScanned_ = true;
  }

  void DumpReader::readFrame(int whichFrame) {
    if (!isScanned_) scanFile();

    int asl  = info_->getSnapshotManager()->getAtomStorageLayout();
    int rbsl = info_->getSnapshotManager()->getRigidBodyStorageLayout();

    needPos_ =
        (asl & DataStorage::dslPosition || rbsl & DataStorage::dslPosition) ?
            true :
            false;
    needVel_ =
        (asl & DataStorage::dslVelocity || rbsl & DataStorage::dslVelocity) ?
            true :
            false;
    needQuaternion_ =
        (asl & DataStorage::dslAmat || asl & DataStorage::dslDipole ||
         asl & DataStorage::dslQuadrupole || rbsl & DataStorage::dslAmat ||
         rbsl & DataStorage::dslDipole || rbsl & DataStorage::dslQuadrupole) ?
            true :
            false;
    needAngMom_ = (asl & DataStorage::dslAngularMomentum ||
                   rbsl & DataStorage::dslAngularMomentum) ?
                      true :
                      false;

    // some dump files contain the efield, but we should only parse
    // and set the field if we have actually allocated memory for it
    readField_ = (asl & DataStorage::dslElectricField) ? true : false;

    readSet(whichFrame);

    if (needCOMprops_) {
      Thermo thermo(info_);
      Vector3d com;

      if (needPos_ && needVel_) {
        Vector3d comvel;
        Vector3d comw;
        thermo.getComAll(com, comvel);
        comw = thermo.getAngularMomentum();
      } else {
        com = thermo.getCom();
      }
    }
  }

  void DumpReader::readSet(int whichFrame) {
    std::string line;

#ifndef IS_MPI
    inFile_.clear();
    inFile_.seekg(framePos_[whichFrame]);

    std::istream& inputStream = inFile_;
#else

    int primaryNode = 0;
    std::stringstream sstream;
    if (worldRank == primaryNode) {
      std::string sendBuffer;

      inFile_.clear();
      inFile_.seekg(framePos_[whichFrame]);

      while (inFile_.getline(buffer, bufferSize)) {
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
               "DumpReader Error: can not find <Snapshot>\n");
      painCave.isFatal = 1;
      simError();
    }

    // read frameData
    readFrameProperties(inputStream);

    // read StuntDoubles
    int nSD = readStuntDoubles(inputStream);

    inputStream.getline(buffer, bufferSize);
    line = buffer;

    if (line.find("<SiteData>") != std::string::npos) {
      // read SiteData
      readSiteData(inputStream);
    } else {
      if (line.find("</Snapshot>") == std::string::npos) {
        snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
                 "DumpReader Error: can not find </Snapshot>\n");
        painCave.isFatal = 1;
        simError();
      }
    }

    if (nSD != info_->getNGlobalIntegrableObjects()) {
      snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
               "DumpReader Error: Number of parsed StuntDouble lines (%d)\n"
               "\tis not the same as the expected number of Objects (%d)\n",
               nSD, info_->getNGlobalIntegrableObjects());
      painCave.isFatal = 1;
      simError();
    }
  }

  void DumpReader::parseDumpLine(const std::string& line) {
    StringTokenizer tokenizer(line);
    int nTokens;

    nTokens = tokenizer.countTokens();

    if (nTokens < 2) {
      snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
               "DumpReader Error: Not enough Tokens.\n%s\n", line.c_str());
      painCave.isFatal = 1;
      simError();
    }

    int index = tokenizer.nextTokenAsInt();

    StuntDouble* sd = info_->getIOIndexToIntegrableObject(index);

    if (sd == NULL) { return; }
    std::string type = tokenizer.nextToken();
    int size         = type.size();

    size_t found;

    if (needPos_) {
      found = type.find("p");
      if (found == std::string::npos) {
        snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
                 "DumpReader Error: StuntDouble %d has no Position\n"
                 "\tField (\"p\") specified.\n%s\n",
                 index, line.c_str());
        painCave.isFatal = 1;
        simError();
      }
    }

    if (sd->isDirectional()) {
      if (needQuaternion_) {
        found = type.find("q");
        if (found == std::string::npos) {
          snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
                   "DumpReader Error: Directional StuntDouble %d has no\n"
                   "\tQuaternion Field (\"q\") specified.\n%s\n",
                   index, line.c_str());
          painCave.isFatal = 1;
          simError();
        }
      }
    }

    for (int i = 0; i < size; ++i) {
      switch (type[i]) {
      case 'p': {
        Vector3d pos;
        pos[0] = tokenizer.nextTokenAsDouble();
        pos[1] = tokenizer.nextTokenAsDouble();
        pos[2] = tokenizer.nextTokenAsDouble();
        if (needPos_) { sd->setPos(pos); }
        break;
      }
      case 'v': {
        Vector3d vel;
        vel[0] = tokenizer.nextTokenAsDouble();
        vel[1] = tokenizer.nextTokenAsDouble();
        vel[2] = tokenizer.nextTokenAsDouble();
        if (needVel_) { sd->setVel(vel); }
        break;
      }

      case 'q': {
        Quat4d q;
        if (sd->isDirectional()) {
          q[0] = tokenizer.nextTokenAsDouble();
          q[1] = tokenizer.nextTokenAsDouble();
          q[2] = tokenizer.nextTokenAsDouble();
          q[3] = tokenizer.nextTokenAsDouble();

          RealType qlen = q.length();
          if (qlen < OpenMD::epsilon) {  // check quaternion is not
            // equal to 0

            snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
                     "DumpReader Error: initial quaternion error "
                     "(q0^2 + q1^2 + q2^2 + q3^2) ~ 0\n");
            painCave.isFatal = 1;
            simError();
          }

          q.normalize();
          if (needQuaternion_) { sd->setQ(q); }
        }
        break;
      }
      case 'j': {
        Vector3d ji;
        if (sd->isDirectional()) {
          ji[0] = tokenizer.nextTokenAsDouble();
          ji[1] = tokenizer.nextTokenAsDouble();
          ji[2] = tokenizer.nextTokenAsDouble();
          if (needAngMom_) { sd->setJ(ji); }
        }
        break;
      }
      case 'f': {
        Vector3d force;
        force[0] = tokenizer.nextTokenAsDouble();
        force[1] = tokenizer.nextTokenAsDouble();
        force[2] = tokenizer.nextTokenAsDouble();
        sd->setFrc(force);
        break;
      }
      case 't': {
        Vector3d torque;
        torque[0] = tokenizer.nextTokenAsDouble();
        torque[1] = tokenizer.nextTokenAsDouble();
        torque[2] = tokenizer.nextTokenAsDouble();
        sd->setTrq(torque);
        break;
      }
      case 'u': {
        RealType particlePot;
        particlePot = tokenizer.nextTokenAsDouble();
        sd->setParticlePot(particlePot);
        break;
      }
      case 'c': {
        RealType flucQPos;
        flucQPos = tokenizer.nextTokenAsDouble();
        if (sd->isAtom())
          if (dynamic_cast<Atom*>(sd)->isFluctuatingCharge())
            sd->setFlucQPos(flucQPos);
        break;
      }
      case 'w': {
        RealType flucQVel;
        flucQVel = tokenizer.nextTokenAsDouble();
        if (sd->isAtom())
          if (dynamic_cast<Atom*>(sd)->isFluctuatingCharge())
            sd->setFlucQVel(flucQVel);
        break;
      }
      case 'g': {
        RealType flucQFrc;
        flucQFrc = tokenizer.nextTokenAsDouble();
        if (sd->isAtom())
          if (dynamic_cast<Atom*>(sd)->isFluctuatingCharge())
            sd->setFlucQFrc(flucQFrc);
        break;
      }
      case 'e': {
        Vector3d eField;
        eField[0] = tokenizer.nextTokenAsDouble();
        eField[1] = tokenizer.nextTokenAsDouble();
        eField[2] = tokenizer.nextTokenAsDouble();
        // only set the field if we have allocated memory for it
        if (readField_) sd->setElectricField(eField);
        break;
      }
      case 's': {
        RealType sPot;
        sPot = tokenizer.nextTokenAsDouble();
        sd->setSitePotential(sPot);
        break;
      }
      case 'd': {
        RealType density;
        density = tokenizer.nextTokenAsDouble();
        sd->setDensity(density);
        break;
      }

      default: {
        snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
                 "DumpReader Error: %s is an unrecognized type\n",
                 type.c_str());
        painCave.isFatal = 1;
        simError();
        break;
      }
      }
    }
    if (sd->isRigidBody()) {
      RigidBody* rb = static_cast<RigidBody*>(sd);
      if (needPos_) {
        // This should let us use various atom-based selections even
        // if we have only rigid bodies:
        rb->updateAtoms();
      }
      if (needVel_) { rb->updateAtomVel(); }
    }
  }

  void DumpReader::parseSiteLine(const std::string& line) {
    StringTokenizer tokenizer(line);
    int nTokens;

    nTokens = tokenizer.countTokens();

    if (nTokens < 1) {
      snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
               "DumpReader Error: Not enough Tokens.\n%s\n", line.c_str());
      painCave.isFatal = 1;
      simError();
    }

    /**
     * The first token is the global integrable object index.
     */

    int index       = tokenizer.nextTokenAsInt();
    StuntDouble* sd = info_->getIOIndexToIntegrableObject(index);
    if (sd == NULL) { return; }

    // StuntDoubles have a line even if there is nothing stored in the
    // site data:
    if (nTokens == 1) { return; }

    /**
     * Test to see if the next token is an integer or not.  If not,
     * we've got data on the integrable object itself.  If there is an
     * integer, we're parsing data for a site on a rigid body.
     */
    std::string indexTest = tokenizer.peekNextToken();
    std::istringstream i(indexTest);
    int siteIndex;
    if (i >> siteIndex) {
      // chew up this token and parse as an int:
      siteIndex = tokenizer.nextTokenAsInt();
      if (sd->isRigidBody()) {
        RigidBody* rb = static_cast<RigidBody*>(sd);

        // Sometimes site lines are inherited from other models, so
        // just ignore a site line that exceeds the number of atoms in
        // our RB:
        if (siteIndex >= static_cast<int>(rb->getNumAtoms())) { return; }

        sd = rb->getAtoms()[siteIndex];
      }
    }

    /**
     * The next token contains information on what follows.
     */
    std::string type = tokenizer.nextToken();
    int size         = type.size();

    for (int i = 0; i < size; ++i) {
      switch (type[i]) {
      case 'u': {
        RealType particlePot;
        particlePot = tokenizer.nextTokenAsDouble();
        sd->setParticlePot(particlePot);
        break;
      }
      case 'c': {
        RealType flucQPos;
        flucQPos = tokenizer.nextTokenAsDouble();
        if (sd->isAtom())
          if (dynamic_cast<Atom*>(sd)->isFluctuatingCharge())
            sd->setFlucQPos(flucQPos);
        break;
      }
      case 'w': {
        RealType flucQVel;
        flucQVel = tokenizer.nextTokenAsDouble();
        if (sd->isAtom())
          if (dynamic_cast<Atom*>(sd)->isFluctuatingCharge())
            sd->setFlucQVel(flucQVel);
        break;
      }
      case 'g': {
        RealType flucQFrc;
        flucQFrc = tokenizer.nextTokenAsDouble();
        if (sd->isAtom())
          if (dynamic_cast<Atom*>(sd)->isFluctuatingCharge())
            sd->setFlucQFrc(flucQFrc);
        break;
      }
      case 'e': {
        Vector3d eField;
        eField[0] = tokenizer.nextTokenAsDouble();
        eField[1] = tokenizer.nextTokenAsDouble();
        eField[2] = tokenizer.nextTokenAsDouble();
        if (readField_) sd->setElectricField(eField);
        break;
      }
      case 's': {
        RealType sPot;
        sPot = tokenizer.nextTokenAsDouble();
        sd->setSitePotential(sPot);
        break;
      }
      case 'd': {
        RealType dens;
        dens = tokenizer.nextTokenAsDouble();
        sd->setDensity(dens);
        break;
      }
      default: {
        snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
                 "DumpReader Error: %s is an unrecognized type\n",
                 type.c_str());
        painCave.isFatal = 1;
        simError();
        break;
      }
      }
    }
  }

  int DumpReader::readStuntDoubles(std::istream& inputStream) {
    inputStream.getline(buffer, bufferSize);
    std::string line(buffer);

    if (line.find("<StuntDoubles>") == std::string::npos) {
      snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
               "DumpReader Error: Missing <StuntDoubles>\n");
      painCave.isFatal = 1;
      simError();
    }

    int nSD = 0;

    while (inputStream.getline(buffer, bufferSize)) {
      line = buffer;

      if (line.find("</StuntDoubles>") != std::string::npos) { break; }

      parseDumpLine(line);
      nSD++;
    }

    return nSD;
  }

  void DumpReader::readSiteData(std::istream& inputStream) {
    std::string line(buffer);

    // We already found the starting <SiteData> tag or we wouldn't be
    // here, so just start parsing until we get to the ending
    // </SiteData> tag:

    while (inputStream.getline(buffer, bufferSize)) {
      line = buffer;

      if (line.find("</SiteData>") != std::string::npos) { break; }

      parseSiteLine(line);
    }
  }

  void DumpReader::readFrameProperties(std::istream& inputStream) {
    Snapshot* s = info_->getSnapshotManager()->getCurrentSnapshot();
    // We're about to overwrite all frame properties, so clear out any
    // derived properties from previous use:
    s->clearDerivedProperties();

    inputStream.getline(buffer, bufferSize);
    std::string line(buffer);

    if (line.find("<FrameData>") == std::string::npos) {
      snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
               "DumpReader Error: Missing <FrameData>\n");
      painCave.isFatal = 1;
      simError();
    }

    while (inputStream.getline(buffer, bufferSize)) {
      line = buffer;

      if (line.find("</FrameData>") != std::string::npos) { break; }

      StringTokenizer tokenizer(line, " ;\t\n\r{}:,");
      if (!tokenizer.hasMoreTokens()) {
        snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
                 "DumpReader Error: Not enough Tokens.\n%s\n", line.c_str());
        painCave.isFatal = 1;
        simError();
      }

      std::string propertyName = tokenizer.nextToken();
      if (propertyName == "Time") {
        RealType currTime = tokenizer.nextTokenAsDouble();
        s->setTime(currTime);
      } else if (propertyName == "Hmat") {
        Mat3x3d hmat;
        hmat(0, 0) = tokenizer.nextTokenAsDouble();
        hmat(0, 1) = tokenizer.nextTokenAsDouble();
        hmat(0, 2) = tokenizer.nextTokenAsDouble();
        hmat(1, 0) = tokenizer.nextTokenAsDouble();
        hmat(1, 1) = tokenizer.nextTokenAsDouble();
        hmat(1, 2) = tokenizer.nextTokenAsDouble();
        hmat(2, 0) = tokenizer.nextTokenAsDouble();
        hmat(2, 1) = tokenizer.nextTokenAsDouble();
        hmat(2, 2) = tokenizer.nextTokenAsDouble();
        s->setHmat(hmat);
      } else if (propertyName == "Thermostat") {
        pair<RealType, RealType> thermostat;
        thermostat.first  = tokenizer.nextTokenAsDouble();
        thermostat.second = tokenizer.nextTokenAsDouble();
        s->setThermostat(thermostat);
      } else if (propertyName == "Barostat") {
        Mat3x3d eta;
        eta(0, 0) = tokenizer.nextTokenAsDouble();
        eta(0, 1) = tokenizer.nextTokenAsDouble();
        eta(0, 2) = tokenizer.nextTokenAsDouble();
        eta(1, 0) = tokenizer.nextTokenAsDouble();
        eta(1, 1) = tokenizer.nextTokenAsDouble();
        eta(1, 2) = tokenizer.nextTokenAsDouble();
        eta(2, 0) = tokenizer.nextTokenAsDouble();
        eta(2, 1) = tokenizer.nextTokenAsDouble();
        eta(2, 2) = tokenizer.nextTokenAsDouble();
        s->setBarostat(eta);
      } else if (propertyName == "SPFData") {
        std::shared_ptr<SPFData> spfData = s->getSPFData();

        spfData->pos[0]   = tokenizer.nextTokenAsDouble();
        spfData->pos[1]   = tokenizer.nextTokenAsDouble();
        spfData->pos[2]   = tokenizer.nextTokenAsDouble();
        spfData->lambda   = tokenizer.nextTokenAsDouble();
        spfData->globalID = tokenizer.nextTokenAsInt();
      } else {
        snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
                 "DumpReader Error: %s is an invalid property in <FrameData>\n",
                 propertyName.c_str());
        painCave.isFatal = 0;
        simError();
      }
    }
  }
}  // namespace OpenMD
