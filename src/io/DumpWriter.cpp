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

#include "io/DumpWriter.hpp"

#ifdef IS_MPI
#include <mpi.h>
#endif

#include <config.h>

#include "io/Globals.hpp"
#include "io/basic_teebuf.hpp"
#include "primitives/Molecule.hpp"
#include "utils/simError.h"
#ifdef HAVE_ZLIB
#include "io/gzstream.hpp"
#endif

using namespace std;
namespace OpenMD {

  DumpWriter::DumpWriter(SimInfo* info) :
      info_(info), filename_(info->getDumpFileName()),
      eorFilename_(info->getFinalConfigFileName()) {
    Globals* simParams = info->getSimParams();
    needCompression_   = simParams->getCompressDumpFile();
    needForceVector_   = simParams->getOutputForceVector();
    needParticlePot_   = simParams->getOutputParticlePotential();
    needFlucQ_         = simParams->getOutputFluctuatingCharges();
    needElectricField_ = simParams->getOutputElectricField();
    needSitePotential_ = simParams->getOutputSitePotential();
    needDensity_       = simParams->getOutputDensity();

    if (needParticlePot_ || needFlucQ_ || needElectricField_ ||
        needSitePotential_ || needDensity_) {
      doSiteData_ = true;
    } else {
      doSiteData_ = false;
    }

    createDumpFile_ = true;
#ifdef HAVE_LIBZ
    if (needCompression_) {
      filename_ += ".gz";
      eorFilename_ += ".gz";
    }
#endif

#ifdef IS_MPI

    if (worldRank == 0) {
#endif  // is_mpi

      dumpFile_ = createOStream(filename_);

      if (!dumpFile_) {
        snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
                 "Could not open \"%s\" for dump output.\n", filename_.c_str());
        painCave.isFatal = 1;
        simError();
      }

#ifdef IS_MPI
    }

#endif  // is_mpi
  }

  DumpWriter::DumpWriter(SimInfo* info, const std::string& filename) :
      info_(info), filename_(filename) {
    Globals* simParams = info->getSimParams();
    eorFilename_       = filename_.substr(0, filename_.rfind(".")) + ".eor";

    needCompression_   = simParams->getCompressDumpFile();
    needForceVector_   = simParams->getOutputForceVector();
    needParticlePot_   = simParams->getOutputParticlePotential();
    needFlucQ_         = simParams->getOutputFluctuatingCharges();
    needElectricField_ = simParams->getOutputElectricField();
    needSitePotential_ = simParams->getOutputSitePotential();
    needDensity_       = simParams->getOutputDensity();

    if (needParticlePot_ || needFlucQ_ || needElectricField_ ||
        needSitePotential_ || needDensity_) {
      doSiteData_ = true;
    } else {
      doSiteData_ = false;
    }

    createDumpFile_ = true;
#ifdef HAVE_LIBZ
    if (needCompression_) {
      filename_ += ".gz";
      eorFilename_ += ".gz";
    }
#endif

#ifdef IS_MPI

    if (worldRank == 0) {
#endif  // is_mpi

      dumpFile_ = createOStream(filename_);

      if (!dumpFile_) {
        snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
                 "Could not open \"%s\" for dump output.\n", filename_.c_str());
        painCave.isFatal = 1;
        simError();
      }

#ifdef IS_MPI
    }

#endif  // is_mpi
  }

  DumpWriter::DumpWriter(SimInfo* info, const std::string& filename,
                         bool writeDumpFile) :
      info_(info),
      filename_(filename) {
    Globals* simParams = info->getSimParams();
    eorFilename_       = filename_.substr(0, filename_.rfind(".")) + ".eor";

    needCompression_   = simParams->getCompressDumpFile();
    needForceVector_   = simParams->getOutputForceVector();
    needParticlePot_   = simParams->getOutputParticlePotential();
    needFlucQ_         = simParams->getOutputFluctuatingCharges();
    needElectricField_ = simParams->getOutputElectricField();
    needSitePotential_ = simParams->getOutputSitePotential();
    needDensity_       = simParams->getOutputDensity();

    if (needParticlePot_ || needFlucQ_ || needElectricField_ ||
        needSitePotential_ || needDensity_) {
      doSiteData_ = true;
    } else {
      doSiteData_ = false;
    }

#ifdef HAVE_LIBZ
    if (needCompression_) {
      filename_ += ".gz";
      eorFilename_ += ".gz";
    }
#endif

#ifdef IS_MPI

    if (worldRank == 0) {
#endif  // is_mpi

      createDumpFile_ = writeDumpFile;
      if (createDumpFile_) {
        dumpFile_ = createOStream(filename_);

        if (!dumpFile_) {
          snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
                   "Could not open \"%s\" for dump output.\n",
                   filename_.c_str());
          painCave.isFatal = 1;
          simError();
        }
      }
#ifdef IS_MPI
    }

#endif  // is_mpi
  }

  DumpWriter::~DumpWriter() {
#ifdef IS_MPI

    if (worldRank == 0) {
#endif  // is_mpi
      if (createDumpFile_) {
        writeClosing(*dumpFile_);
        delete dumpFile_;
      }
#ifdef IS_MPI
    }

#endif  // is_mpi
  }

  void DumpWriter::writeFrameProperties(std::ostream& os, Snapshot* s) {
    char buffer[1024];

    os << "    <FrameData>\n";

    RealType currentTime = s->getTime();

    if (std::isinf(currentTime) || std::isnan(currentTime)) {
      snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
               "DumpWriter detected a numerical error writing the time");
      painCave.isFatal = 1;
      simError();
    }

    snprintf(buffer, 1024, "        Time: %.10g\n", currentTime);
    os << buffer;

    Mat3x3d hmat;
    hmat = s->getHmat();

    for (unsigned int i = 0; i < 3; i++) {
      for (unsigned int j = 0; j < 3; j++) {
        if (std::isinf(hmat(i, j)) || std::isnan(hmat(i, j))) {
          snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
                   "DumpWriter detected a numerical error writing the box");
          painCave.isFatal = 1;
          simError();
        }
      }
    }

    snprintf(
        buffer, 1024,
        "        Hmat: {{ %.10g, %.10g, %.10g }, { %.10g, %.10g, %.10g }, { "
        "%.10g, "
        "%.10g, %.10g }}\n",
        hmat(0, 0), hmat(1, 0), hmat(2, 0), hmat(0, 1), hmat(1, 1), hmat(2, 1),
        hmat(0, 2), hmat(1, 2), hmat(2, 2));
    os << buffer;

    pair<RealType, RealType> thermostat = s->getThermostat();

    if (std::isinf(thermostat.first) || std::isnan(thermostat.first) ||
        std::isinf(thermostat.second) || std::isnan(thermostat.second)) {
      snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
               "DumpWriter detected a numerical error writing the thermostat");
      painCave.isFatal = 1;
      simError();
    }
    snprintf(buffer, 1024, "  Thermostat: %.10g , %.10g\n", thermostat.first,
             thermostat.second);
    os << buffer;

    Mat3x3d eta;
    eta = s->getBarostat();

    for (unsigned int i = 0; i < 3; i++) {
      for (unsigned int j = 0; j < 3; j++) {
        if (std::isinf(eta(i, j)) || std::isnan(eta(i, j))) {
          snprintf(
              painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
              "DumpWriter detected a numerical error writing the barostat");
          painCave.isFatal = 1;
          simError();
        }
      }
    }

    snprintf(
        buffer, 1024,
        "    Barostat: {{ %.10g, %.10g, %.10g }, { %.10g, %.10g, %.10g }, { "
        "%.10g, "
        "%.10g, %.10g }}\n",
        eta(0, 0), eta(1, 0), eta(2, 0), eta(0, 1), eta(1, 1), eta(2, 1),
        eta(0, 2), eta(1, 2), eta(2, 2));
    os << buffer;

    // SPF Data
    std::shared_ptr<SPFData> spfData = s->getSPFData();

    if (std::isinf(spfData->pos[0]) || std::isnan(spfData->pos[0]) ||
        std::isinf(spfData->pos[1]) || std::isnan(spfData->pos[1]) ||
        std::isinf(spfData->pos[2]) || std::isnan(spfData->pos[2]) ||
        std::isinf(spfData->lambda) || std::isnan(spfData->lambda)) {
      snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
               "DumpWriter detected a numerical error writing the spf data "
               "structure");
      painCave.isFatal = 1;
      simError();
    }
    snprintf(buffer, 1024,
             "     SPFData: {{ %.10g, %.10g, %.10g }, %.10g, %d }\n",
             spfData->pos[0], spfData->pos[1], spfData->pos[2], spfData->lambda,
             spfData->globalID);
    os << buffer;

    os << "    </FrameData>\n";
  }

  void DumpWriter::writeFrame(std::ostream& os) {
#ifdef IS_MPI
    MPI_Status istatus;
#endif

    Molecule* mol;
    StuntDouble* sd;
    SimInfo::MoleculeIterator mi;
    Molecule::IntegrableObjectIterator ii;
    RigidBody::AtomIterator ai;

#ifndef IS_MPI
    os << "  <Snapshot>\n";

    writeFrameProperties(os, info_->getSnapshotManager()->getCurrentSnapshot());

    os << "    <StuntDoubles>\n";
    for (mol = info_->beginMolecule(mi); mol != NULL;
         mol = info_->nextMolecule(mi)) {
      for (sd = mol->beginIntegrableObject(ii); sd != NULL;
           sd = mol->nextIntegrableObject(ii)) {
        os << prepareDumpLine(sd);
      }
    }
    os << "    </StuntDoubles>\n";

    if (doSiteData_) {
      os << "    <SiteData>\n";
      for (mol = info_->beginMolecule(mi); mol != NULL;
           mol = info_->nextMolecule(mi)) {
        for (sd = mol->beginIntegrableObject(ii); sd != NULL;
             sd = mol->nextIntegrableObject(ii)) {
          int ioIndex = sd->getGlobalIntegrableObjectIndex();
          // do one for the IO itself
          os << prepareSiteLine(sd, ioIndex, 0);

          if (sd->isRigidBody()) {
            RigidBody* rb = static_cast<RigidBody*>(sd);
            int siteIndex = 0;
            for (Atom* atom = rb->beginAtom(ai); atom != NULL;
                 atom       = rb->nextAtom(ai)) {
              os << prepareSiteLine(atom, ioIndex, siteIndex);
              siteIndex++;
            }
          }
        }
      }
      os << "    </SiteData>\n";
    }
    os << "  </Snapshot>\n";

    os.flush();
    os.rdbuf()->pubsync();
#else

    const int primaryNode = 0;
    int worldRank;
    int nProc;

    MPI_Comm_size(MPI_COMM_WORLD, &nProc);
    MPI_Comm_rank(MPI_COMM_WORLD, &worldRank);

    if (worldRank == primaryNode) {
      os << "  <Snapshot>\n";
      writeFrameProperties(os,
                           info_->getSnapshotManager()->getCurrentSnapshot());
      os << "    <StuntDoubles>\n";
    }

    // every node prepares the dump lines for integrable objects belong to
    // itself
    std::string buffer;
    for (mol = info_->beginMolecule(mi); mol != NULL;
         mol = info_->nextMolecule(mi)) {
      for (sd = mol->beginIntegrableObject(ii); sd != NULL;
           sd = mol->nextIntegrableObject(ii)) {
        buffer += prepareDumpLine(sd);
      }
    }

    if (worldRank == primaryNode) {
      os << buffer;

      for (int i = 1; i < nProc; ++i) {
        // tell processor i to start sending us data:

        MPI_Bcast(&i, 1, MPI_INT, primaryNode, MPI_COMM_WORLD);

        // receive the length of the string buffer that was
        // prepared by processor i:
        int recvLength;
        MPI_Recv(&recvLength, 1, MPI_INT, i, MPI_ANY_TAG, MPI_COMM_WORLD,
                 &istatus);

        // create a buffer to receive the data
        char* recvBuffer = new char[recvLength];
        if (recvBuffer == NULL) {
        } else {
          // receive the data:
          MPI_Recv(recvBuffer, recvLength, MPI_CHAR, i, MPI_ANY_TAG,
                   MPI_COMM_WORLD, &istatus);
          // send it to the file:
          os << recvBuffer;
          // get rid of the receive buffer:
          delete[] recvBuffer;
        }
      }
    } else {
      int sendBufferLength = buffer.size() + 1;
      int myturn           = 0;
      for (int i = 1; i < nProc; ++i) {
        // wait for the primary node to call our number:
        MPI_Bcast(&myturn, 1, MPI_INT, primaryNode, MPI_COMM_WORLD);
        if (myturn == worldRank) {
          // send the length of our buffer:

          MPI_Send(&sendBufferLength, 1, MPI_INT, primaryNode, 0,
                   MPI_COMM_WORLD);

          // send our buffer:
          MPI_Send((void*)buffer.c_str(), sendBufferLength, MPI_CHAR,
                   primaryNode, 0, MPI_COMM_WORLD);
        }
      }
    }

    if (worldRank == primaryNode) { os << "    </StuntDoubles>\n"; }

    if (doSiteData_) {
      if (worldRank == primaryNode) { os << "    <SiteData>\n"; }
      buffer.clear();
      for (mol = info_->beginMolecule(mi); mol != NULL;
           mol = info_->nextMolecule(mi)) {
        for (sd = mol->beginIntegrableObject(ii); sd != NULL;
             sd = mol->nextIntegrableObject(ii)) {
          int ioIndex = sd->getGlobalIntegrableObjectIndex();
          // do one for the IO itself
          buffer += prepareSiteLine(sd, ioIndex, 0);

          if (sd->isRigidBody()) {
            RigidBody* rb = static_cast<RigidBody*>(sd);
            int siteIndex = 0;
            for (Atom* atom = rb->beginAtom(ai); atom != NULL;
                 atom       = rb->nextAtom(ai)) {
              buffer += prepareSiteLine(atom, ioIndex, siteIndex);
              siteIndex++;
            }
          }
        }
      }

      if (worldRank == primaryNode) {
        os << buffer;

        for (int i = 1; i < nProc; ++i) {
          // tell processor i to start sending us data:
          MPI_Bcast(&i, 1, MPI_INT, primaryNode, MPI_COMM_WORLD);

          // receive the length of the string buffer that was
          // prepared by processor i:
          int recvLength;
          MPI_Recv(&recvLength, 1, MPI_INT, i, MPI_ANY_TAG, MPI_COMM_WORLD,
                   &istatus);

          // create a buffer to receive the data
          char* recvBuffer = new char[recvLength];
          if (recvBuffer == NULL) {
          } else {
            // receive the data:
            MPI_Recv(recvBuffer, recvLength, MPI_CHAR, i, MPI_ANY_TAG,
                     MPI_COMM_WORLD, &istatus);
            // send it to the file:
            os << recvBuffer;
            // get rid of the receive buffer:
            delete[] recvBuffer;
          }
        }
      } else {
        int sendBufferLength = buffer.size() + 1;
        int myturn           = 0;
        for (int i = 1; i < nProc; ++i) {
          // wait for the primary node to call our number:
          MPI_Bcast(&myturn, 1, MPI_INT, primaryNode, MPI_COMM_WORLD);
          if (myturn == worldRank) {
            // send the length of our buffer:
            MPI_Send(&sendBufferLength, 1, MPI_INT, primaryNode, 0,
                     MPI_COMM_WORLD);
            // send our buffer:
            MPI_Send((void*)buffer.c_str(), sendBufferLength, MPI_CHAR,
                     primaryNode, 0, MPI_COMM_WORLD);
          }
        }
      }

      if (worldRank == primaryNode) { os << "    </SiteData>\n"; }
    }

    if (worldRank == primaryNode) {
      os << "  </Snapshot>\n";
      os.flush();
      os.rdbuf()->pubsync();
    }

#endif  // is_mpi
  }

  std::string DumpWriter::prepareDumpLine(StuntDouble* sd) {
    int index = sd->getGlobalIntegrableObjectIndex();
    std::string type("pv");
    std::string line;
    char tempBuffer[4096];

    Vector3d pos;
    Vector3d vel;
    pos = sd->getPos();

    if (std::isinf(pos[0]) || std::isnan(pos[0]) || std::isinf(pos[1]) ||
        std::isnan(pos[1]) || std::isinf(pos[2]) || std::isnan(pos[2])) {
      snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
               "DumpWriter detected a numerical error writing the position"
               " for object %d",
               index);
      painCave.isFatal = 1;
      simError();
    }

    vel = sd->getVel();

    if (std::isinf(vel[0]) || std::isnan(vel[0]) || std::isinf(vel[1]) ||
        std::isnan(vel[1]) || std::isinf(vel[2]) || std::isnan(vel[2])) {
      snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
               "DumpWriter detected a numerical error writing the velocity"
               " for object %d",
               index);
      painCave.isFatal = 1;
      simError();
    }

    snprintf(tempBuffer, 4096, "%18.10g %18.10g %18.10g %13e %13e %13e", pos[0],
             pos[1], pos[2], vel[0], vel[1], vel[2]);
    line += tempBuffer;

    if (sd->isDirectional()) {
      type += "qj";
      Quat4d q;
      Vector3d ji;
      q = sd->getQ();

      if (std::isinf(q[0]) || std::isnan(q[0]) || std::isinf(q[1]) ||
          std::isnan(q[1]) || std::isinf(q[2]) || std::isnan(q[2]) ||
          std::isinf(q[3]) || std::isnan(q[3])) {
        snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
                 "DumpWriter detected a numerical error writing the quaternion"
                 " for object %d",
                 index);
        painCave.isFatal = 1;
        simError();
      }

      ji = sd->getJ();

      if (std::isinf(ji[0]) || std::isnan(ji[0]) || std::isinf(ji[1]) ||
          std::isnan(ji[1]) || std::isinf(ji[2]) || std::isnan(ji[2])) {
        snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
                 "DumpWriter detected a numerical error writing the angular"
                 " momentum for object %d",
                 index);
        painCave.isFatal = 1;
        simError();
      }

      snprintf(tempBuffer, 4096, " %13e %13e %13e %13e %13e %13e %13e", q[0],
               q[1], q[2], q[3], ji[0], ji[1], ji[2]);
      line += tempBuffer;
    }

    if (needForceVector_) {
      type += "f";
      Vector3d frc = sd->getFrc();
      if (std::isinf(frc[0]) || std::isnan(frc[0]) || std::isinf(frc[1]) ||
          std::isnan(frc[1]) || std::isinf(frc[2]) || std::isnan(frc[2])) {
        snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
                 "DumpWriter detected a numerical error writing the force"
                 " for object %d",
                 index);
        painCave.isFatal = 1;
        simError();
      }
      snprintf(tempBuffer, 4096, " %13e %13e %13e", frc[0], frc[1], frc[2]);
      line += tempBuffer;

      if (sd->isDirectional()) {
        type += "t";
        Vector3d trq = sd->getTrq();
        if (std::isinf(trq[0]) || std::isnan(trq[0]) || std::isinf(trq[1]) ||
            std::isnan(trq[1]) || std::isinf(trq[2]) || std::isnan(trq[2])) {
          snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
                   "DumpWriter detected a numerical error writing the torque"
                   " for object %d",
                   index);
          painCave.isFatal = 1;
          simError();
        }
        snprintf(tempBuffer, 4096, " %13e %13e %13e", trq[0], trq[1], trq[2]);
        line += tempBuffer;
      }
    }

    snprintf(tempBuffer, 4096, "%10d %7s %s\n", index, type.c_str(),
             line.c_str());
    return std::string(tempBuffer);
  }

  std::string DumpWriter::prepareSiteLine(StuntDouble* sd, int ioIndex,
                                          int siteIndex) {
    int asl  = info_->getSnapshotManager()->getAtomStorageLayout();
    int rbsl = info_->getSnapshotManager()->getRigidBodyStorageLayout();
    int sl {};

    std::string id;
    std::string type;
    std::string line;
    char tempBuffer[4096];

    if (sd->isRigidBody()) {
      sl = rbsl;
      snprintf(tempBuffer, 4096, "%10d           ", ioIndex);
      id = std::string(tempBuffer);
    } else {
      sl = asl;
      snprintf(tempBuffer, 4096, "%10d %10d", ioIndex, siteIndex);
      id = std::string(tempBuffer);
    }

    if (needFlucQ_) {
      if (sl & DataStorage::dslFlucQPosition) {
        type += "c";
        RealType fqPos = sd->getFlucQPos();
        if (std::isinf(fqPos) || std::isnan(fqPos)) {
          snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
                   "DumpWriter detected a numerical error writing the"
                   " fluctuating charge for object %s",
                   id.c_str());
          painCave.isFatal = 1;
          simError();
        }
        snprintf(tempBuffer, 4096, " %13e ", fqPos);
        line += tempBuffer;
      }

      if (sl & DataStorage::dslFlucQVelocity) {
        type += "w";
        RealType fqVel = sd->getFlucQVel();
        if (std::isinf(fqVel) || std::isnan(fqVel)) {
          snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
                   "DumpWriter detected a numerical error writing the"
                   " fluctuating charge velocity for object %s",
                   id.c_str());
          painCave.isFatal = 1;
          simError();
        }
        snprintf(tempBuffer, 4096, " %13e ", fqVel);
        line += tempBuffer;
      }

      if (needForceVector_) {
        if (sl & DataStorage::dslFlucQForce) {
          type += "g";
          RealType fqFrc = sd->getFlucQFrc();
          if (std::isinf(fqFrc) || std::isnan(fqFrc)) {
            snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
                     "DumpWriter detected a numerical error writing the"
                     " fluctuating charge force for object %s",
                     id.c_str());
            painCave.isFatal = 1;
            simError();
          }
          snprintf(tempBuffer, 4096, " %13e ", fqFrc);
          line += tempBuffer;
        }
      }
    }

    if (needElectricField_) {
      if (sl & DataStorage::dslElectricField) {
        type += "e";
        Vector3d eField = sd->getElectricField();
        if (std::isinf(eField[0]) || std::isnan(eField[0]) ||
            std::isinf(eField[1]) || std::isnan(eField[1]) ||
            std::isinf(eField[2]) || std::isnan(eField[2])) {
          snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
                   "DumpWriter detected a numerical error writing the electric"
                   " field for object %s",
                   id.c_str());
          painCave.isFatal = 1;
          simError();
        }
        snprintf(tempBuffer, 4096, " %13e %13e %13e", eField[0], eField[1],
                 eField[2]);
        line += tempBuffer;
      }
    }

    if (needSitePotential_) {
      if (sl & DataStorage::dslSitePotential) {
        type += "s";
        RealType sPot = sd->getSitePotential();
        if (std::isinf(sPot) || std::isnan(sPot)) {
          snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
                   "DumpWriter detected a numerical error writing the"
                   " site potential for object %s",
                   id.c_str());
          painCave.isFatal = 1;
          simError();
        }
        snprintf(tempBuffer, 4096, " %13e ", sPot);
        line += tempBuffer;
      }
    }

    if (needParticlePot_) {
      if (sl & DataStorage::dslParticlePot) {
        type += "u";
        RealType particlePot = sd->getParticlePot();
        if (std::isinf(particlePot) || std::isnan(particlePot)) {
          snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
                   "DumpWriter detected a numerical error writing the particle "
                   " potential for object %s",
                   id.c_str());
          painCave.isFatal = 1;
          simError();
        }
        snprintf(tempBuffer, 4096, " %13e", particlePot);
        line += tempBuffer;
      }
    }

    if (needDensity_) {
      if (sl & DataStorage::dslDensity) {
        type += "d";
        RealType density = sd->getDensity();
        if (std::isinf(density) || std::isnan(density)) {
          snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
                   "DumpWriter detected a numerical error writing the density "
                   " for object %s",
                   id.c_str());
          painCave.isFatal = 1;
          simError();
        }
        snprintf(tempBuffer, 4096, " %13e", density);
        line += tempBuffer;
      }
    }

    snprintf(tempBuffer, 4096, "%s %7s %s\n", id.c_str(), type.c_str(),
             line.c_str());
    return std::string(tempBuffer);
  }

  void DumpWriter::writeDump() { writeFrame(*dumpFile_); }

  void DumpWriter::writeEor() {
    std::ostream* eorStream = NULL;

#ifdef IS_MPI
    if (worldRank == 0) {
#endif  // is_mpi

      eorStream = createOStream(eorFilename_);

#ifdef IS_MPI
    }
#endif

    writeFrame(*eorStream);

#ifdef IS_MPI
    if (worldRank == 0) {
#endif

      writeClosing(*eorStream);
      delete eorStream;

#ifdef IS_MPI
    }
#endif  // is_mpi
  }

  void DumpWriter::writeDumpAndEor() {
    std::vector<std::streambuf*> buffers;
    std::ostream* eorStream = NULL;
#ifdef IS_MPI
    if (worldRank == 0) {
#endif  // is_mpi
      buffers.push_back(dumpFile_->rdbuf());
      eorStream = createOStream(eorFilename_);
      buffers.push_back(eorStream->rdbuf());
#ifdef IS_MPI
    }
#endif  // is_mpi

    TeeBuf tbuf(buffers.begin(), buffers.end());
    std::ostream os(&tbuf);
    writeFrame(os);

#ifdef IS_MPI
    if (worldRank == 0) {
#endif  // is_mpi
      writeClosing(*eorStream);
      delete eorStream;
#ifdef IS_MPI
    }
#endif  // is_mpi
  }

  std::ostream* DumpWriter::createOStream(const std::string& filename) {
    std::ostream* newOStream;
#ifdef HAVE_ZLIB
    if (needCompression_) {
      newOStream = new ogzstream(filename.c_str());
    } else {
      newOStream = new std::ofstream(filename.c_str());
    }
#else
    newOStream = new std::ofstream(filename.c_str());
#endif
    // write out MetaData first
    (*newOStream) << "<OpenMD version=2>" << std::endl;
    (*newOStream) << "  <MetaData>" << std::endl;
    (*newOStream) << info_->getRawMetaData();
    (*newOStream) << "  </MetaData>" << std::endl;
    return newOStream;
  }

  void DumpWriter::writeClosing(std::ostream& os) {
    os << "</OpenMD>\n";
    os.flush();
  }

}  // namespace OpenMD
