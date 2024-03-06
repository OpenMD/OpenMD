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

/**
 * @file SimCreator.cpp
 * @author tlin
 * @date 11/03/2004
 * @version 1.0
 */

#include "brains/SimCreator.hpp"

#include <exception>
#include <iostream>
#include <sstream>
#include <string>

#ifdef IS_MPI
#include <mpi.h>
#endif

#include "antlr/ANTLRException.hpp"
#include "antlr/CharStreamException.hpp"
#include "antlr/MismatchedCharException.hpp"
#include "antlr/MismatchedTokenException.hpp"
#include "antlr/NoViableAltException.hpp"
#include "antlr/NoViableAltForCharException.hpp"
#include "antlr/RecognitionException.hpp"
#include "antlr/TokenStreamException.hpp"
#include "antlr/TokenStreamIOException.hpp"
#include "antlr/TokenStreamRecognitionException.hpp"
#include "brains/ForceField.hpp"
#include "brains/MoleculeCreator.hpp"
#include "brains/SimSnapshotManager.hpp"
#include "io/DumpReader.hpp"
#include "mdParser/MDLexer.hpp"
#include "mdParser/MDParser.hpp"
#include "mdParser/MDTreeParser.hpp"
#include "mdParser/SimplePreprocessor.hpp"
#include "types/DirectionalAdapter.hpp"
#include "types/EAMAdapter.hpp"
#include "types/FixedChargeAdapter.hpp"
#include "types/FluctuatingChargeAdapter.hpp"
#include "types/MultipoleAdapter.hpp"
#include "types/PolarizableAdapter.hpp"
#include "types/SuttonChenAdapter.hpp"
#include "utils/RandNumGen.hpp"
#include "utils/Revision.hpp"
#include "utils/Trim.hpp"
#include "utils/simError.h"

namespace OpenMD {

  Globals* SimCreator::parseFile(std::istream& rawMetaDataStream,
                                 const std::string& filename, int mdFileVersion,
                                 int startOfMetaDataBlock) {
    Globals* simParams = NULL;
    try {
      // Create a preprocessor that preprocesses md file into an ostringstream
      std::stringstream ppStream;
#ifdef IS_MPI
      int streamSize;
      const int primaryNode = 0;

      if (worldRank == primaryNode) {
        MPI_Bcast(&mdFileVersion, 1, MPI_INT, primaryNode, MPI_COMM_WORLD);
#endif
        SimplePreprocessor preprocessor;
        preprocessor.preprocess(rawMetaDataStream, filename,
                                startOfMetaDataBlock, ppStream);

#ifdef IS_MPI
        // broadcasting the stream size
        streamSize = ppStream.str().size() + 1;
        MPI_Bcast(&streamSize, 1, MPI_INT, primaryNode, MPI_COMM_WORLD);
        MPI_Bcast(static_cast<void*>(const_cast<char*>(ppStream.str().c_str())),
                  streamSize, MPI_CHAR, primaryNode, MPI_COMM_WORLD);
      } else {
        MPI_Bcast(&mdFileVersion, 1, MPI_INT, primaryNode, MPI_COMM_WORLD);

        // get stream size
        MPI_Bcast(&streamSize, 1, MPI_INT, primaryNode, MPI_COMM_WORLD);
        char* buf = new char[streamSize];
        assert(buf);

        // receive file content
        MPI_Bcast(buf, streamSize, MPI_CHAR, primaryNode, MPI_COMM_WORLD);

        ppStream.str(buf);
        delete[] buf;
      }
#endif
      // Create a scanner that reads from the input stream
      MDLexer lexer(ppStream);
      lexer.setFilename(filename);
      lexer.initDeferredLineCount();

      // Create a parser that reads from the scanner
      MDParser parser(lexer);
      parser.setFilename(filename);

      // Create an observer that synchorizes file name change
      FilenameObserver observer;
      observer.setLexer(&lexer);
      observer.setParser(&parser);
      lexer.setObserver(&observer);
      antlr::ASTFactory factory;
      parser.initializeASTFactory(factory);
      parser.setASTFactory(&factory);
      parser.mdfile();
      // Create a tree parser that reads information into Globals
      MDTreeParser treeParser;
      treeParser.initializeASTFactory(factory);
      treeParser.setASTFactory(&factory);
      simParams = treeParser.walkTree(parser.getAST());
    }

    catch (antlr::MismatchedCharException& e) {
      snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
               "Mismatched Character: %s in file %s at line %d, column %d\n",
               e.getMessage().c_str(), e.getFilename().c_str(), e.getLine(),
               e.getColumn());
      painCave.isFatal = 1;
      simError();
    } catch (antlr::MismatchedTokenException& e) {
      snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
               "Mismatched Token: %s in file %s at line %d, column %d\n",
               e.getMessage().c_str(), e.getFilename().c_str(), e.getLine(),
               e.getColumn());
      painCave.isFatal = 1;
      simError();
    } catch (antlr::NoViableAltForCharException& e) {
      snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
               "No Viable Alternative for Character: %s in file %s at line %d, "
               "column %d\n",
               e.getMessage().c_str(), e.getFilename().c_str(), e.getLine(),
               e.getColumn());
      painCave.isFatal = 1;
      simError();
    } catch (antlr::NoViableAltException& e) {
      snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
               "No Viable Alternative: %s in file %s at line %d, column %d\n",
               e.getMessage().c_str(), e.getFilename().c_str(), e.getLine(),
               e.getColumn());
      painCave.isFatal = 1;
      simError();
    }

    catch (antlr::TokenStreamRecognitionException& e) {
      snprintf(
          painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
          "Token Stream Recognition: %s in file %s at line %d, column %d\n",
          e.getMessage().c_str(), e.getFilename().c_str(), e.getLine(),
          e.getColumn());
      painCave.isFatal = 1;
      simError();
    }

    catch (antlr::TokenStreamIOException& e) {
      snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
               "Token Stream IO exception: %s\n", e.getMessage().c_str());
      painCave.isFatal = 1;
      simError();
    }

    catch (antlr::TokenStreamException& e) {
      snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
               "Token Stream exception: %s\n", e.getMessage().c_str());
      painCave.isFatal = 1;
      simError();
    } catch (antlr::RecognitionException& e) {
      snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
               "Recognition exception: %s in file %s at line %d, column %d\n",
               e.getMessage().c_str(), e.getFilename().c_str(), e.getLine(),
               e.getColumn());
      painCave.isFatal = 1;
      simError();
    } catch (antlr::CharStreamException& e) {
      snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
               "Character Stream exception: %s\n", e.getMessage().c_str());
      painCave.isFatal = 1;
      simError();
    } catch (OpenMDException& e) {
      snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH, "%s\n", e.what());
      painCave.isFatal = 1;
      simError();
    } catch (std::exception& e) {
      snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
               "parser exception: %s\n", e.what());
      painCave.isFatal = 1;
      simError();
    }

    simParams->setMDfileVersion(mdFileVersion);
    return simParams;
  }

  SimInfo* SimCreator::createSim(const std::string& mdFileName,
                                 bool loadInitCoords) {
    const int bufferSize = 65535;
    char buffer[bufferSize];
    int lineNo = 0;
    std::string mdRawData;
    int metaDataBlockStart = -1;
    int metaDataBlockEnd   = -1;
    streamoff mdOffset {};
    int mdFileVersion(2);

    // Create a string for embedding the version information in the MetaData
    std::string version;
    version.assign("## Last run using OpenMD version: ");
    version.append(OPENMD_VERSION_MAJOR);
    version.append(".");
    version.append(OPENMD_VERSION_MINOR);
    version.append(",");

    std::string rev(revision, strnlen(revision, 40));
    rev.append(40 - rev.length(), ' ');

    version.append(" revision: ");
    // If there's no GIT revision, just call this the RELEASE revision.
    if (!rev.empty()) {
      version.append(rev);
    } else {
      version.append("RELEASE");
    }

#ifdef IS_MPI
    const int primaryNode = 0;
    if (worldRank == primaryNode) {
#endif

      std::ifstream mdFile_;
      mdFile_.open(mdFileName.c_str(), ifstream::in | ifstream::binary);

      if (mdFile_.fail()) {
        snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
                 "SimCreator: Cannot open file: %s\n", mdFileName.c_str());
        painCave.isFatal = 1;
        simError();
      }

      mdFile_.getline(buffer, bufferSize);
      ++lineNo;
      std::string line = Utils::trimLeftCopy(buffer);
      std::size_t i    = CaseInsensitiveFind(line, "<OpenMD");
      if (i == string::npos) {
        // try the older file strings to see if that works:
        i = CaseInsensitiveFind(line, "<OOPSE");
      }

      if (i == string::npos) {
        // still no luck!
        snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
                 "SimCreator: File: %s is not a valid OpenMD file!\n",
                 mdFileName.c_str());
        painCave.isFatal = 1;
        simError();
      }

      // found the correct opening string, now try to get the file
      // format version number.

      StringTokenizer tokenizer(line, "=<> \t\n\r");
      std::string fileType = tokenizer.nextToken();
      toUpper(fileType);

      mdFileVersion = 0;

      if (fileType == "OPENMD") {
        while (tokenizer.hasMoreTokens()) {
          std::string token(tokenizer.nextToken());
          toUpper(token);
          if (token == "VERSION") {
            mdFileVersion = tokenizer.nextTokenAsInt();
            break;
          }
        }
      }

      bool startFound(false);
      bool endFound(false);
      // scan through the input stream and find MetaData tag
      while (mdFile_.getline(buffer, bufferSize)) {
        ++lineNo;

        std::string line = Utils::trimLeftCopy(buffer);
        if (metaDataBlockStart == -1) {
          std::size_t i = CaseInsensitiveFind(line, "<MetaData>");
          if (i != string::npos) {
            metaDataBlockStart = lineNo;
            startFound         = true;
            mdOffset           = mdFile_.tellg();
          }
        } else {
          std::size_t i = CaseInsensitiveFind(line, "</MetaData>");
          if (i != string::npos) {
            metaDataBlockEnd = lineNo;
            endFound         = true;
          }
        }
        if (startFound && endFound) break;
      }

      if (metaDataBlockStart == -1) {
        snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
                 "SimCreator: File: %s did not contain a <MetaData> tag!\n",
                 mdFileName.c_str());
        painCave.isFatal = 1;
        simError();
      }
      if (metaDataBlockEnd == -1) {
        snprintf(
            painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
            "SimCreator: File: %s did not contain a closed MetaData block!\n",
            mdFileName.c_str());
        painCave.isFatal = 1;
        simError();
      }

      mdFile_.clear();
      mdFile_.seekg(0);
      mdFile_.seekg(mdOffset);

      mdRawData.clear();

      bool foundVersion = false;

      for (int i = 0; i < metaDataBlockEnd - metaDataBlockStart - 1; ++i) {
        mdFile_.getline(buffer, bufferSize);
        std::string line = Utils::trimLeftCopy(buffer);
        std::size_t j =
            CaseInsensitiveFind(line, "## Last run using OpenMD Version");
        if (j != string::npos) {
          foundVersion = true;
          mdRawData += version;
        } else {
          mdRawData += buffer;
        }
        mdRawData += "\n";
      }

      if (!foundVersion) mdRawData += version + "\n";

      mdFile_.close();

#ifdef IS_MPI
    }
#endif

    std::stringstream rawMetaDataStream(mdRawData);

    // parse meta-data file
    Globals* simParams = parseFile(rawMetaDataStream, mdFileName, mdFileVersion,
                                   metaDataBlockStart + 1);

    // create the force field
    ForceField* ff = new ForceField(simParams->getForceField());

    if (ff == NULL) {
      snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
               "ForceField Factory can not create %s force field\n",
               simParams->getForceField().c_str());
      painCave.isFatal = 1;
      simError();
    }

    if (simParams->haveForceFieldFileName()) {
      ff->setForceFieldFileName(simParams->getForceFieldFileName());
    }

    std::string forcefieldFileName;
    forcefieldFileName = ff->getForceFieldFileName();

    if (simParams->haveForceFieldVariant()) {
      // If the force field has variant, the variant force field name will be
      // Base.variant.frc. For exampel EAM.u6.frc

      std::string variant = simParams->getForceFieldVariant();

      std::string::size_type pos = forcefieldFileName.rfind(".frc");
      variant                    = "." + variant;
      if (pos != std::string::npos) {
        forcefieldFileName.insert(pos, variant);
      } else {
        // If the default force field file name does not containt .frc suffix,
        // just append the .variant
        forcefieldFileName.append(variant);
      }
    }

    ff->parse(forcefieldFileName);
    // create SimInfo
    SimInfo* info = new SimInfo(ff, simParams);

    info->setRawMetaData(mdRawData);

    // gather parameters (SimCreator only retrieves part of the
    // parameters)
    gatherParameters(info, mdFileName);

    // divide the molecules and determine the global index of molecules
#ifdef IS_MPI
    divideMolecules(info);
#endif

    // create the molecules
    createMolecules(info);

    // find the storage layout

    computeStorageLayouts(info);

    int asl  = info->getAtomStorageLayout();
    int rbsl = info->getRigidBodyStorageLayout();
    int cgsl = info->getCutoffGroupStorageLayout();

    // allocate memory for DataStorage(circular reference, need to
    // break it)
    info->setSnapshotManager(new SimSnapshotManager(info, asl, rbsl, cgsl));

    // set the global index of atoms, rigidbodies and cutoffgroups
    //(only need to be set once, the global index will never change
    // again). Local indices of atoms and rigidbodies are already set
    // by MoleculeCreator class which actually delegates the
    // responsibility to LocalIndexManager.
    setGlobalIndex(info);

    // Although addInteractionPairs is called inside SimInfo's addMolecule
    // method, at that point atoms don't have the global index yet
    //(their global index are all initialized to -1).  Therefore we
    // have to call addInteractionPairs explicitly here. A way to work
    // around is that we can determine the beginning global indices of
    // atoms before they get created.
    SimInfo::MoleculeIterator mi;
    Molecule* mol;
    for (mol = info->beginMolecule(mi); mol != NULL;
         mol = info->nextMolecule(mi)) {
      info->addInteractionPairs(mol);
    }

    if (loadInitCoords) loadCoordinates(info, mdFileName);
    return info;
  }

  void SimCreator::gatherParameters(SimInfo* info, const std::string& mdfile) {
    // figure out the output file names
    std::string prefix;

#ifdef IS_MPI

    if (worldRank == 0) {
#endif  // is_mpi
      Globals* simParams = info->getSimParams();
      if (simParams->haveFinalConfig()) {
        prefix = getPrefix(simParams->getFinalConfig());
      } else {
        prefix = getPrefix(mdfile);
      }

      info->setFinalConfigFileName(prefix + ".eor");
      info->setDumpFileName(prefix + ".dump");
      info->setStatFileName(prefix + ".stat");
      info->setReportFileName(prefix + ".report");
      info->setRestFileName(prefix + ".zang");

#ifdef IS_MPI
    }
#endif
  }

#ifdef IS_MPI
  void SimCreator::divideMolecules(SimInfo* info) {
    RealType a;
    int nProcessors;
    std::vector<int> atomsPerProc;
    int nGlobalMols = info->getNGlobalMolecules();
    std::vector<int> molToProcMap(nGlobalMols, -1);  // default to an error
                                                     // condition:

    MPI_Comm_size(MPI_COMM_WORLD, &nProcessors);

    if (nProcessors > nGlobalMols) {
      snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
               "nProcessors (%d) > nMol (%d)\n"
               "\tThe number of processors is larger than\n"
               "\tthe number of molecules.  This will not result in a \n"
               "\tusable division of atoms for force decomposition.\n"
               "\tEither try a smaller number of processors, or run the\n"
               "\tsingle-processor version of OpenMD.\n",
               nProcessors, nGlobalMols);

      painCave.isFatal = 1;
      simError();
    }

    a = 3.0 * nGlobalMols / info->getNGlobalAtoms();

    // initialize atomsPerProc
    atomsPerProc.insert(atomsPerProc.end(), nProcessors, 0);

    if (worldRank == 0) {
      Utils::RandNumGenPtr myRandom = info->getRandomNumberGenerator();

      std::uniform_int_distribution<> processorDistribution {0,
                                                             nProcessors - 1};
      std::uniform_real_distribution<RealType> yDistribution {0, 1};

      RealType numerator   = info->getNGlobalAtoms();
      RealType denominator = nProcessors;
      RealType precast     = numerator / denominator;
      int nTarget          = (int)(precast + 0.5);
      int which_proc {0};

      for (int i = 0; i < nGlobalMols; i++) {
        // get the molecule stamp first
        int stampId                  = info->getMoleculeStampId(i);
        MoleculeStamp* moleculeStamp = info->getMoleculeStamp(stampId);
        int add_atoms                = moleculeStamp->getNAtoms();

        if (nProcessors > 1) {
          int done  = 0;
          int loops = 0;

          while (!done) {
            loops++;

            // Pick a processor at random
            which_proc = processorDistribution(*myRandom);

            // How many atoms does this processor have so far?
            int old_atoms = atomsPerProc[which_proc];
            int new_atoms = old_atoms + add_atoms;

            // If we've been through this loop too many times, we need
            // to just give up and assign the molecule to this processor
            // and be done with it.

            if (loops > 100) {
              snprintf(
                  painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
                  "There have been 100 attempts to assign molecule %d to an\n"
                  "\tunderworked processor, but there's no good place to\n"
                  "\tleave it. OpenMD is assigning it at random to processor "
                  "%d.\n",
                  i, which_proc);

              painCave.isFatal  = 0;
              painCave.severity = OPENMD_INFO;
              simError();

              molToProcMap[i] = which_proc;
              atomsPerProc[which_proc] += add_atoms;

              done = 1;
              continue;
            }

            // If we can add this molecule to this processor without sending
            // it above nTarget, then go ahead and do it:

            if (new_atoms <= nTarget) {
              molToProcMap[i] = which_proc;
              atomsPerProc[which_proc] += add_atoms;

              done = 1;
              continue;
            }

            // The only situation left is when new_atoms > nTarget.  We
            // want to accept this with some probability that dies off the
            // farther we are from nTarget

            // roughly:  x = new_atoms - nTarget
            //           Pacc(x) = exp(- a * x)
            // where a = penalty / (average atoms per molecule)

            RealType x = (RealType)(new_atoms - nTarget);
            RealType y = yDistribution(*myRandom);

            if (y < exp(-a * x)) {
              molToProcMap[i] = which_proc;
              atomsPerProc[which_proc] += add_atoms;

              done = 1;
              continue;
            } else {
              continue;
            }
          }
        } else {
          which_proc      = 0;
          molToProcMap[i] = which_proc;
          atomsPerProc[which_proc] += add_atoms;
        }
      }

      // Spray out this nonsense to all other processors:
      MPI_Bcast(&molToProcMap[0], nGlobalMols, MPI_INT, 0, MPI_COMM_WORLD);

    } else {
      // Listen to your marching orders from processor 0:
      MPI_Bcast(&molToProcMap[0], nGlobalMols, MPI_INT, 0, MPI_COMM_WORLD);
    }

    info->setMolToProcMap(molToProcMap);
    snprintf(checkPointMsg, MAX_SIM_ERROR_MSG_LENGTH,
             "Successfully divided the molecules among the processors.\n");
    errorCheckPoint();
  }

#endif

  void SimCreator::createMolecules(SimInfo* info) {
    MoleculeCreator molCreator;
    int stampId;

    for (int i = 0; i < info->getNGlobalMolecules(); i++) {
#ifdef IS_MPI
      if (info->getMolToProc(i) == worldRank) {
#endif

        stampId       = info->getMoleculeStampId(i);
        Molecule* mol = molCreator.createMolecule(
            info->getForceField(), info->getMoleculeStamp(stampId), stampId, i,
            info->getLocalIndexManager());

        info->addMolecule(mol);

#ifdef IS_MPI
      }
#endif
    }
  }

  void SimCreator::computeStorageLayouts(SimInfo* info) {
    Globals* simParams    = info->getSimParams();
    int nRigidBodies      = info->getNGlobalRigidBodies();
    AtomTypeSet atomTypes = info->getSimulatedAtomTypes();
    AtomTypeSet::iterator i;
    bool hasDirectionalAtoms     = false;
    bool hasFixedCharge          = false;
    bool hasDipoles              = false;
    bool hasQuadrupoles          = false;
    bool hasPolarizable          = false;
    bool hasFluctuatingCharge    = false;
    bool hasMetallic             = false;
    int atomStorageLayout        = 0;
    int rigidBodyStorageLayout   = 0;
    int cutoffGroupStorageLayout = 0;

    atomStorageLayout |= DataStorage::dslPosition;
    atomStorageLayout |= DataStorage::dslVelocity;
    atomStorageLayout |= DataStorage::dslForce;
    cutoffGroupStorageLayout |= DataStorage::dslPosition;

    for (i = atomTypes.begin(); i != atomTypes.end(); ++i) {
      DirectionalAdapter da        = DirectionalAdapter((*i));
      MultipoleAdapter ma          = MultipoleAdapter((*i));
      EAMAdapter ea                = EAMAdapter((*i));
      SuttonChenAdapter sca        = SuttonChenAdapter((*i));
      PolarizableAdapter pa        = PolarizableAdapter((*i));
      FixedChargeAdapter fca       = FixedChargeAdapter((*i));
      FluctuatingChargeAdapter fqa = FluctuatingChargeAdapter((*i));

      if (da.isDirectional()) { hasDirectionalAtoms = true; }
      if (ma.isDipole()) { hasDipoles = true; }
      if (ma.isQuadrupole()) { hasQuadrupoles = true; }
      if (ea.isEAM() || sca.isSuttonChen()) { hasMetallic = true; }
      if (fca.isFixedCharge()) { hasFixedCharge = true; }
      if (fqa.isFluctuatingCharge()) { hasFluctuatingCharge = true; }
      if (pa.isPolarizable()) { hasPolarizable = true; }
    }

    if (nRigidBodies > 0) {
      rigidBodyStorageLayout |= DataStorage::dslPosition;
      rigidBodyStorageLayout |= DataStorage::dslVelocity;
      rigidBodyStorageLayout |= DataStorage::dslForce;
      rigidBodyStorageLayout |= DataStorage::dslAmat;
      rigidBodyStorageLayout |= DataStorage::dslAngularMomentum;
      rigidBodyStorageLayout |= DataStorage::dslTorque;
    }
    if (hasDirectionalAtoms) {
      atomStorageLayout |= DataStorage::dslAmat;
      if (atomStorageLayout & DataStorage::dslVelocity) {
        atomStorageLayout |= DataStorage::dslAngularMomentum;
      }
      if (atomStorageLayout & DataStorage::dslForce) {
        atomStorageLayout |= DataStorage::dslTorque;
      }
    }
    if (hasDipoles) { atomStorageLayout |= DataStorage::dslDipole; }
    if (hasQuadrupoles) { atomStorageLayout |= DataStorage::dslQuadrupole; }
    if (hasFixedCharge || hasFluctuatingCharge) {
      atomStorageLayout |= DataStorage::dslSkippedCharge;
    }
    if (hasMetallic) {
      atomStorageLayout |= DataStorage::dslDensity;
      atomStorageLayout |= DataStorage::dslFunctional;
      atomStorageLayout |= DataStorage::dslFunctionalDerivative;
    }
    if (hasPolarizable) { atomStorageLayout |= DataStorage::dslElectricField; }
    if (hasFluctuatingCharge) {
      atomStorageLayout |= DataStorage::dslFlucQPosition;
      if (atomStorageLayout & DataStorage::dslVelocity) {
        atomStorageLayout |= DataStorage::dslFlucQVelocity;
      }
      if (atomStorageLayout & DataStorage::dslForce) {
        atomStorageLayout |= DataStorage::dslFlucQForce;
      }
    }

    // if the user has asked for them, make sure we've got the memory for the
    // objects defined.

    if (simParams->getOutputParticlePotential()) {
      atomStorageLayout |= DataStorage::dslParticlePot;
    }

    if (simParams->havePrintHeatFlux()) {
      if (simParams->getPrintHeatFlux()) {
        atomStorageLayout |= DataStorage::dslParticlePot;
      }
    }

    if (simParams->getOutputElectricField() | simParams->haveElectricField() |
        simParams->haveUniformField() |
        simParams->haveUniformGradientStrength() |
        simParams->haveUniformGradientDirection1() |
        simParams->haveUniformGradientDirection2() |
	simParams->getLightParameters()->getUseLight() ) {
      atomStorageLayout |= DataStorage::dslElectricField;
      rigidBodyStorageLayout |= DataStorage::dslElectricField;
    }

    if (simParams->getRNEMDParameters()->haveUseRNEMD()) {
      if (simParams->getRNEMDParameters()->getUseRNEMD()) {
        if (simParams->getRNEMDParameters()->requiresElectricField()) {
          atomStorageLayout |= DataStorage::dslElectricField;
          rigidBodyStorageLayout |= DataStorage::dslElectricField;
        }
      }
    }

    if (simParams->getOutputSitePotential()) {
      atomStorageLayout |= DataStorage::dslSitePotential;
      rigidBodyStorageLayout |= DataStorage::dslSitePotential;
    }

    if (simParams->getOutputFluctuatingCharges()) {
      atomStorageLayout |= DataStorage::dslFlucQPosition;
      atomStorageLayout |= DataStorage::dslFlucQVelocity;
      atomStorageLayout |= DataStorage::dslFlucQForce;
    }

    info->setAtomStorageLayout(atomStorageLayout);
    info->setRigidBodyStorageLayout(rigidBodyStorageLayout);
    info->setCutoffGroupStorageLayout(cutoffGroupStorageLayout);

    return;
  }

  void SimCreator::setGlobalIndex(SimInfo* info) {
    SimInfo::MoleculeIterator mi;
    Molecule::AtomIterator ai;
    Molecule::RigidBodyIterator ri;
    Molecule::CutoffGroupIterator ci;
    Molecule::BondIterator boi;
    Molecule::BendIterator bei;
    Molecule::TorsionIterator ti;
    Molecule::InversionIterator ii;
    Molecule::IntegrableObjectIterator ioi;
    Molecule* mol;
    Atom* atom;
    RigidBody* rb;
    CutoffGroup* cg;
    Bond* bond;
    Bend* bend;
    Torsion* torsion;
    Inversion* inversion;
    int beginAtomIndex;
    int beginRigidBodyIndex;
    int beginCutoffGroupIndex;
    int beginBondIndex;
    int beginBendIndex;
    int beginTorsionIndex;
    int beginInversionIndex;
#ifdef IS_MPI
    int nGlobalAtoms       = info->getNGlobalAtoms();
    int nGlobalRigidBodies = info->getNGlobalRigidBodies();
#endif

    beginAtomIndex = 0;
    // The rigid body indices begin immediately after the atom indices:
    beginRigidBodyIndex   = info->getNGlobalAtoms();
    beginCutoffGroupIndex = 0;
    beginBondIndex        = 0;
    beginBendIndex        = 0;
    beginTorsionIndex     = 0;
    beginInversionIndex   = 0;

    for (int i = 0; i < info->getNGlobalMolecules(); i++) {
#ifdef IS_MPI
      if (info->getMolToProc(i) == worldRank) {
#endif
        // stuff to do if I own this molecule
        mol = info->getMoleculeByGlobalIndex(i);

        // The local index(index in DataStorge) of the atom is important:
        for (atom = mol->beginAtom(ai); atom != NULL;
             atom = mol->nextAtom(ai)) {
          atom->setGlobalIndex(beginAtomIndex++);
        }

        for (rb = mol->beginRigidBody(ri); rb != NULL;
             rb = mol->nextRigidBody(ri)) {
          rb->setGlobalIndex(beginRigidBodyIndex++);
        }

        // The local index of other objects only depends on the order
        // of traversal:
        for (cg = mol->beginCutoffGroup(ci); cg != NULL;
             cg = mol->nextCutoffGroup(ci)) {
          cg->setGlobalIndex(beginCutoffGroupIndex++);
        }
        for (bond = mol->beginBond(boi); bond != NULL;
             bond = mol->nextBond(boi)) {
          bond->setGlobalIndex(beginBondIndex++);
        }
        for (bend = mol->beginBend(bei); bend != NULL;
             bend = mol->nextBend(bei)) {
          bend->setGlobalIndex(beginBendIndex++);
        }
        for (torsion = mol->beginTorsion(ti); torsion != NULL;
             torsion = mol->nextTorsion(ti)) {
          torsion->setGlobalIndex(beginTorsionIndex++);
        }
        for (inversion = mol->beginInversion(ii); inversion != NULL;
             inversion = mol->nextInversion(ii)) {
          inversion->setGlobalIndex(beginInversionIndex++);
        }

#ifdef IS_MPI
      } else {
        // stuff to do if I don't own this molecule

        int stampId          = info->getMoleculeStampId(i);
        MoleculeStamp* stamp = info->getMoleculeStamp(stampId);

        beginAtomIndex += stamp->getNAtoms();
        beginRigidBodyIndex += stamp->getNRigidBodies();
        beginCutoffGroupIndex +=
            stamp->getNCutoffGroups() + stamp->getNFreeAtoms();
        beginBondIndex += stamp->getNBonds();
        beginBendIndex += stamp->getNBends();
        beginTorsionIndex += stamp->getNTorsions();
        beginInversionIndex += stamp->getNInversions();
      }
#endif

    }  // end for(int i=0)

    // fill globalGroupMembership
    std::vector<int> globalGroupMembership(info->getNGlobalAtoms(), 0);
    for (mol = info->beginMolecule(mi); mol != NULL;
         mol = info->nextMolecule(mi)) {
      for (cg = mol->beginCutoffGroup(ci); cg != NULL;
           cg = mol->nextCutoffGroup(ci)) {
        for (atom = cg->beginAtom(ai); atom != NULL; atom = cg->nextAtom(ai)) {
          globalGroupMembership[atom->getGlobalIndex()] = cg->getGlobalIndex();
        }
      }
    }

#ifdef IS_MPI
    // Since the globalGroupMembership has been zero filled and we've only
    // poked values into the atoms we know, we can do an Allreduce
    // to get the full globalGroupMembership array (We think).
    // This would be prettier if we could use MPI_IN_PLACE like the MPI-2
    // docs said we could.
    std::vector<int> tmpGroupMembership(info->getNGlobalAtoms(), 0);
    MPI_Allreduce(&globalGroupMembership[0], &tmpGroupMembership[0],
                  nGlobalAtoms, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

    info->setGlobalGroupMembership(tmpGroupMembership);
#else
    info->setGlobalGroupMembership(globalGroupMembership);
#endif

    // fill molMembership
    std::vector<int> globalMolMembership(
        info->getNGlobalAtoms() + info->getNGlobalRigidBodies(), 0);

    for (mol = info->beginMolecule(mi); mol != NULL;
         mol = info->nextMolecule(mi)) {
      for (atom = mol->beginAtom(ai); atom != NULL; atom = mol->nextAtom(ai)) {
        globalMolMembership[atom->getGlobalIndex()] = mol->getGlobalIndex();
      }
      for (rb = mol->beginRigidBody(ri); rb != NULL;
           rb = mol->nextRigidBody(ri)) {
        globalMolMembership[rb->getGlobalIndex()] = mol->getGlobalIndex();
      }
    }

#ifdef IS_MPI
    std::vector<int> tmpMolMembership(
        info->getNGlobalAtoms() + info->getNGlobalRigidBodies(), 0);
    MPI_Allreduce(&globalMolMembership[0], &tmpMolMembership[0],
                  nGlobalAtoms + nGlobalRigidBodies, MPI_INT, MPI_SUM,
                  MPI_COMM_WORLD);

    info->setGlobalMolMembership(tmpMolMembership);
#else
    info->setGlobalMolMembership(globalMolMembership);
#endif

    // nIOPerMol holds the number of integrable objects per molecule
    // here the molecules are listed by their global indices.

    std::vector<int> nIOPerMol(info->getNGlobalMolecules(), 0);
    for (mol = info->beginMolecule(mi); mol != NULL;
         mol = info->nextMolecule(mi)) {
      nIOPerMol[mol->getGlobalIndex()] = mol->getNIntegrableObjects();
    }

#ifdef IS_MPI
    std::vector<int> numIntegrableObjectsPerMol(info->getNGlobalMolecules(), 0);
    MPI_Allreduce(&nIOPerMol[0], &numIntegrableObjectsPerMol[0],
                  info->getNGlobalMolecules(), MPI_INT, MPI_SUM,
                  MPI_COMM_WORLD);
#else
    std::vector<int> numIntegrableObjectsPerMol = nIOPerMol;
#endif

    std::vector<int> startingIOIndexForMol(info->getNGlobalMolecules());

    int startingIndex = 0;
    for (int i = 0; i < info->getNGlobalMolecules(); i++) {
      startingIOIndexForMol[i] = startingIndex;
      startingIndex += numIntegrableObjectsPerMol[i];
    }

    std::vector<StuntDouble*> IOIndexToIntegrableObject(
        info->getNGlobalIntegrableObjects(), (StuntDouble*)NULL);
    for (mol = info->beginMolecule(mi); mol != NULL;
         mol = info->nextMolecule(mi)) {
      int myGlobalIndex = mol->getGlobalIndex();
      int globalIO      = startingIOIndexForMol[myGlobalIndex];
      for (StuntDouble* sd = mol->beginIntegrableObject(ioi); sd != NULL;
           sd              = mol->nextIntegrableObject(ioi)) {
        sd->setGlobalIntegrableObjectIndex(globalIO);
        IOIndexToIntegrableObject[globalIO] = sd;
        globalIO++;
      }
    }

    info->setIOIndexToIntegrableObject(IOIndexToIntegrableObject);
  }

  void SimCreator::loadCoordinates(SimInfo* info,
                                   const std::string& mdFileName) {
    DumpReader reader(info, mdFileName);
    int nframes = reader.getNFrames();

    if (nframes > 0) {
      reader.readFrame(nframes - 1);
    } else {
      // invalid initial coordinate file
      snprintf(
          painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
          "Initial configuration file %s should at least contain one frame\n",
          mdFileName.c_str());
      painCave.isFatal = 1;
      simError();
    }
    // copy the current snapshot to previous snapshot
    info->getSnapshotManager()->advance();
  }

}  // namespace OpenMD
