/*
 * Copyright (c) 2005 The University of Notre Dame. All Rights Reserved.
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
 * [4]  Kuang & Gezelter,  J. Chem. Phys. 133, 164101 (2010).
 * [5]  Vardeman, Stocker & Gezelter, J. Chem. Theory Comput. 7, 834 (2011).
 */
 
/**
 * @file SimCreator.cpp
 * @author tlin
 * @date 11/03/2004
 * @version 1.0
 */
#include <exception>
#include <iostream>
#include <sstream>
#include <string>

#include "brains/MoleculeCreator.hpp"
#include "brains/SimCreator.hpp"
#include "brains/SimSnapshotManager.hpp"
#include "io/DumpReader.hpp"
#include "brains/ForceField.hpp"
#include "utils/simError.h"
#include "utils/StringUtils.hpp"
#include "math/SeqRandNumGen.hpp"
#include "mdParser/MDLexer.hpp"
#include "mdParser/MDParser.hpp"
#include "mdParser/MDTreeParser.hpp"
#include "mdParser/SimplePreprocessor.hpp"
#include "antlr/ANTLRException.hpp"
#include "antlr/TokenStreamRecognitionException.hpp"
#include "antlr/TokenStreamIOException.hpp"
#include "antlr/TokenStreamException.hpp"
#include "antlr/RecognitionException.hpp"
#include "antlr/CharStreamException.hpp"

#include "antlr/MismatchedCharException.hpp"
#include "antlr/MismatchedTokenException.hpp"
#include "antlr/NoViableAltForCharException.hpp"
#include "antlr/NoViableAltException.hpp"

#include "types/DirectionalAdapter.hpp"
#include "types/MultipoleAdapter.hpp"
#include "types/EAMAdapter.hpp"
#include "types/SuttonChenAdapter.hpp"
#include "types/PolarizableAdapter.hpp"
#include "types/FixedChargeAdapter.hpp"
#include "types/FluctuatingChargeAdapter.hpp"

#ifdef IS_MPI
#include "mpi.h"
#include "math/ParallelRandNumGen.hpp"
#endif

namespace OpenMD {
  
  Globals* SimCreator::parseFile(std::istream& rawMetaDataStream, const std::string& filename, int mdFileVersion, int startOfMetaDataBlock ){
    Globals* simParams = NULL;
    try {

      // Create a preprocessor that preprocesses md file into an ostringstream
      std::stringstream ppStream;
#ifdef IS_MPI            
      int streamSize;
      const int masterNode = 0;
      int commStatus;
      if (worldRank == masterNode) {
        commStatus = MPI_Bcast(&mdFileVersion, 1, MPI_INT, masterNode, MPI_COMM_WORLD);
#endif                 
        SimplePreprocessor preprocessor;
        preprocessor.preprocess(rawMetaDataStream, filename, startOfMetaDataBlock, ppStream);
                
#ifdef IS_MPI            
        //brocasting the stream size
        streamSize = ppStream.str().size() +1;
        commStatus = MPI_Bcast(&streamSize, 1, MPI_LONG, masterNode, MPI_COMM_WORLD);                   

        commStatus = MPI_Bcast(static_cast<void*>(const_cast<char*>(ppStream.str().c_str())), streamSize, MPI_CHAR, masterNode, MPI_COMM_WORLD); 
            
                
      } else {

        commStatus = MPI_Bcast(&mdFileVersion, 1, MPI_INT, masterNode, MPI_COMM_WORLD);

        //get stream size
        commStatus = MPI_Bcast(&streamSize, 1, MPI_LONG, masterNode, MPI_COMM_WORLD);   

        char* buf = new char[streamSize];
        assert(buf);
                
        //receive file content
        commStatus = MPI_Bcast(buf, streamSize, MPI_CHAR, masterNode, MPI_COMM_WORLD); 
                
        ppStream.str(buf);
        delete [] buf;

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

      
    catch(antlr::MismatchedCharException& e) {
      sprintf(painCave.errMsg, 
              "parser exception: %s %s:%d:%d\n",
              e.getMessage().c_str(),e.getFilename().c_str(), e.getLine(), e.getColumn());
      painCave.isFatal = 1;
      simError();           
    }
    catch(antlr::MismatchedTokenException &e) {
      sprintf(painCave.errMsg, 
              "parser exception: %s %s:%d:%d\n",
              e.getMessage().c_str(),e.getFilename().c_str(), e.getLine(), e.getColumn());
      painCave.isFatal = 1;
      simError();   
    }
    catch(antlr::NoViableAltForCharException &e) {
      sprintf(painCave.errMsg, 
              "parser exception: %s %s:%d:%d\n",
              e.getMessage().c_str(),e.getFilename().c_str(), e.getLine(), e.getColumn());
      painCave.isFatal = 1;
      simError();   
    }
    catch(antlr::NoViableAltException &e) {
      sprintf(painCave.errMsg, 
              "parser exception: %s %s:%d:%d\n",
              e.getMessage().c_str(),e.getFilename().c_str(), e.getLine(), e.getColumn());
      painCave.isFatal = 1;
      simError();   
    }
      
    catch(antlr::TokenStreamRecognitionException& e) {
      sprintf(painCave.errMsg, 
              "parser exception: %s %s:%d:%d\n",
              e.getMessage().c_str(),e.getFilename().c_str(), e.getLine(), e.getColumn());
      painCave.isFatal = 1;
      simError();   
    }
	
    catch(antlr::TokenStreamIOException& e) {
      sprintf(painCave.errMsg, 
              "parser exception: %s\n",
              e.getMessage().c_str());
      painCave.isFatal = 1;
      simError();
    }
	
    catch(antlr::TokenStreamException& e) {
      sprintf(painCave.errMsg, 
              "parser exception: %s\n",
              e.getMessage().c_str());
      painCave.isFatal = 1;
      simError();
    }        
    catch (antlr::RecognitionException& e) {
      sprintf(painCave.errMsg, 
              "parser exception: %s %s:%d:%d\n",
              e.getMessage().c_str(),e.getFilename().c_str(), e.getLine(), e.getColumn());
      painCave.isFatal = 1;
      simError();          
    }
    catch (antlr::CharStreamException& e) {
      sprintf(painCave.errMsg, 
              "parser exception: %s\n",
              e.getMessage().c_str());
      painCave.isFatal = 1;
      simError();        
    }
    catch (OpenMDException& e) {
      sprintf(painCave.errMsg, 
              "%s\n",
              e.getMessage().c_str());
      painCave.isFatal = 1;
      simError();
    }
    catch (std::exception& e) {
      sprintf(painCave.errMsg, 
              "parser exception: %s\n",
              e.what());
      painCave.isFatal = 1;
      simError();
    }

    simParams->setMDfileVersion(mdFileVersion);
    return simParams;
  }
  
  SimInfo*  SimCreator::createSim(const std::string & mdFileName, 
                                  bool loadInitCoords) {
    
    const int bufferSize = 65535;
    char buffer[bufferSize];
    int lineNo = 0;
    std::string mdRawData;
    int metaDataBlockStart = -1;
    int metaDataBlockEnd = -1;
    int i, j;
    streamoff mdOffset;
    int mdFileVersion;

    // Create a string for embedding the version information in the MetaData
    std::string version;
    version.assign("## Last run using OpenMD Version: ");
    version.append(OPENMD_VERSION_MAJOR);
    version.append(".");
    version.append(OPENMD_VERSION_MINOR);

    std::string svnrev;
    //convert a macro from compiler to a string in c++
    STR_DEFINE(svnrev, SVN_REV );
    version.append(" Revision: ");
    // If there's no SVN revision, just call this the RELEASE revision.
    if (!svnrev.empty()) {
      version.append(svnrev);
    } else {
      version.append("RELEASE");
    }
   
#ifdef IS_MPI            
    const int masterNode = 0;
    if (worldRank == masterNode) {
#endif 

      std::ifstream mdFile_;
      mdFile_.open(mdFileName.c_str(), ifstream::in | ifstream::binary);
      
      if (mdFile_.fail()) { 
        sprintf(painCave.errMsg, 
                "SimCreator: Cannot open file: %s\n", 
                mdFileName.c_str()); 
        painCave.isFatal = 1; 
        simError(); 
      } 

      mdFile_.getline(buffer, bufferSize);
      ++lineNo;
      std::string line = trimLeftCopy(buffer);
      i = CaseInsensitiveFind(line, "<OpenMD");
      if (static_cast<size_t>(i) == string::npos) {
        // try the older file strings to see if that works:
        i = CaseInsensitiveFind(line, "<OOPSE");
      }
      
      if (static_cast<size_t>(i) == string::npos) {
        // still no luck!
        sprintf(painCave.errMsg, 
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
            
      //scan through the input stream and find MetaData tag        
      while(mdFile_.getline(buffer, bufferSize)) {
        ++lineNo;
        
        std::string line = trimLeftCopy(buffer);
        if (metaDataBlockStart == -1) {
          i = CaseInsensitiveFind(line, "<MetaData>");
          if (i != string::npos) {
            metaDataBlockStart = lineNo;
            mdOffset = mdFile_.tellg();
          }
        } else {
          i = CaseInsensitiveFind(line, "</MetaData>");
          if (i != string::npos) {
            metaDataBlockEnd = lineNo;
          }
        }
      }

      if (metaDataBlockStart == -1) {
        sprintf(painCave.errMsg, 
                "SimCreator: File: %s did not contain a <MetaData> tag!\n", 
                mdFileName.c_str()); 
        painCave.isFatal = 1; 
        simError(); 
      }
      if (metaDataBlockEnd == -1) {
        sprintf(painCave.errMsg, 
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
        std::string line = trimLeftCopy(buffer);
        j = CaseInsensitiveFind(line, "## Last run using OpenMD Version");
        if (static_cast<size_t>(j) != string::npos) {
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

    //parse meta-data file
    Globals* simParams = parseFile(rawMetaDataStream, mdFileName, mdFileVersion,
                                   metaDataBlockStart + 1);
    
    //create the force field
    ForceField * ff = new ForceField(simParams->getForceField());

    if (ff == NULL) {
      sprintf(painCave.errMsg, 
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
      //If the force field has variant, the variant force field name will be
      //Base.variant.frc. For exampel EAM.u6.frc
      
      std::string variant = simParams->getForceFieldVariant();
      
      std::string::size_type pos = forcefieldFileName.rfind(".frc");
      variant = "." + variant;
      if (pos != std::string::npos) {
        forcefieldFileName.insert(pos, variant);
      } else {
        //If the default force field file name does not containt .frc suffix, just append the .variant
        forcefieldFileName.append(variant);
      }
    } 
    
    ff->parse(forcefieldFileName);
    //create SimInfo
    SimInfo * info = new SimInfo(ff, simParams);

    info->setRawMetaData(mdRawData);
     
    //gather parameters (SimCreator only retrieves part of the
    //parameters)
    gatherParameters(info, mdFileName);
    
    //divide the molecules and determine the global index of molecules
#ifdef IS_MPI
    divideMolecules(info);
#endif 
    
    //create the molecules
    createMolecules(info);
    
    //find the storage layout

    int storageLayout = computeStorageLayout(info);

    //allocate memory for DataStorage(circular reference, need to
    //break it)
    info->setSnapshotManager(new SimSnapshotManager(info, storageLayout));
    
    //set the global index of atoms, rigidbodies and cutoffgroups
    //(only need to be set once, the global index will never change
    //again). Local indices of atoms and rigidbodies are already set
    //by MoleculeCreator class which actually delegates the
    //responsibility to LocalIndexManager.
    setGlobalIndex(info);
    
    //Although addInteractionPairs is called inside SimInfo's addMolecule
    //method, at that point atoms don't have the global index yet
    //(their global index are all initialized to -1).  Therefore we
    //have to call addInteractionPairs explicitly here. A way to work
    //around is that we can determine the beginning global indices of
    //atoms before they get created.
    SimInfo::MoleculeIterator mi;
    Molecule* mol;
    for (mol= info->beginMolecule(mi); mol != NULL; mol = info->nextMolecule(mi)) {
      info->addInteractionPairs(mol);
    }
    
    if (loadInitCoords)
      loadCoordinates(info, mdFileName);    
    return info;
  }
  
  void SimCreator::gatherParameters(SimInfo *info, const std::string& mdfile) {
    
    //figure out the output file names
    std::string prefix;
    
#ifdef IS_MPI
    
    if (worldRank == 0) {
#endif // is_mpi
      Globals * simParams = info->getSimParams();
      if (simParams->haveFinalConfig()) {
        prefix = getPrefix(simParams->getFinalConfig());
      } else {
        prefix = getPrefix(mdfile);
      }
      
      info->setFinalConfigFileName(prefix + ".eor");
      info->setDumpFileName(prefix + ".dump");
      info->setStatFileName(prefix + ".stat");
      info->setRestFileName(prefix + ".zang");
      
#ifdef IS_MPI
      
    }
    
#endif
    
  }
  
#ifdef IS_MPI
  void SimCreator::divideMolecules(SimInfo *info) {
    RealType numerator;
    RealType denominator;
    RealType precast;
    RealType x;
    RealType y;
    RealType a;
    int old_atoms;
    int add_atoms;
    int new_atoms;
    int nTarget;
    int done;
    int i;
    int j;
    int loops;
    int which_proc;
    int nProcessors;
    std::vector<int> atomsPerProc;
    int nGlobalMols = info->getNGlobalMolecules();
    std::vector<int> molToProcMap(nGlobalMols, -1); // default to an error condition:
    
    nProcessors = MPI::COMM_WORLD.Get_size();
    
    if (nProcessors > nGlobalMols) {
      sprintf(painCave.errMsg,
              "nProcessors (%d) > nMol (%d)\n"
              "\tThe number of processors is larger than\n"
              "\tthe number of molecules.  This will not result in a \n"
              "\tusable division of atoms for force decomposition.\n"
              "\tEither try a smaller number of processors, or run the\n"
              "\tsingle-processor version of OpenMD.\n", nProcessors, nGlobalMols);
      
      painCave.isFatal = 1;
      simError();
    }
    
    int seedValue;
    Globals * simParams = info->getSimParams();
    SeqRandNumGen* myRandom; //divide labor does not need Parallel random number generator
    if (simParams->haveSeed()) {
      seedValue = simParams->getSeed();
      myRandom = new SeqRandNumGen(seedValue);
    }else {
      myRandom = new SeqRandNumGen();
    }   
    
    
    a = 3.0 * nGlobalMols / info->getNGlobalAtoms();
    
    //initialize atomsPerProc
    atomsPerProc.insert(atomsPerProc.end(), nProcessors, 0);
    
    if (worldRank == 0) {
      numerator = info->getNGlobalAtoms();
      denominator = nProcessors;
      precast = numerator / denominator;
      nTarget = (int)(precast + 0.5);
      
      for(i = 0; i < nGlobalMols; i++) {

        done = 0;
        loops = 0;
        
        while (!done) {
          loops++;
          
          // Pick a processor at random
          
          which_proc = (int) (myRandom->rand() * nProcessors);
          
          //get the molecule stamp first
          int stampId = info->getMoleculeStampId(i);
          MoleculeStamp * moleculeStamp = info->getMoleculeStamp(stampId);
          
          // How many atoms does this processor have so far?
          old_atoms = atomsPerProc[which_proc];
          add_atoms = moleculeStamp->getNAtoms();
          new_atoms = old_atoms + add_atoms;
          
          // If we've been through this loop too many times, we need
          // to just give up and assign the molecule to this processor
          // and be done with it. 
          
          if (loops > 100) {

            sprintf(painCave.errMsg,
                    "There have been 100 attempts to assign molecule %d to an\n"
                    "\tunderworked processor, but there's no good place to\n"
                    "\tleave it.  OpenMD is assigning it at random to processor %d.\n",
                    i, which_proc);
           
            painCave.isFatal = 0;
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
          
          x = (RealType)(new_atoms - nTarget);
          y = myRandom->rand();
          
          if (y < exp(- a * x)) {
            molToProcMap[i] = which_proc;
            atomsPerProc[which_proc] += add_atoms;
            
            done = 1;
            continue;
          } else {
            continue;
          }
        }
      }
      
      delete myRandom;

      // Spray out this nonsense to all other processors:
      MPI::COMM_WORLD.Bcast(&molToProcMap[0], nGlobalMols, MPI::INT, 0);
    } else {
      
      // Listen to your marching orders from processor 0:
      MPI::COMM_WORLD.Bcast(&molToProcMap[0], nGlobalMols, MPI::INT, 0);

    }
    
    info->setMolToProcMap(molToProcMap);
    sprintf(checkPointMsg,
            "Successfully divided the molecules among the processors.\n");
    errorCheckPoint();
  }
  
#endif
  
  void SimCreator::createMolecules(SimInfo *info) {
    MoleculeCreator molCreator;
    int stampId;
    
    for(int i = 0; i < info->getNGlobalMolecules(); i++) {
      
#ifdef IS_MPI
      
      if (info->getMolToProc(i) == worldRank) {
#endif
        
        stampId = info->getMoleculeStampId(i);
        Molecule * mol = molCreator.createMolecule(info->getForceField(), 
                                                   info->getMoleculeStamp(stampId),
                                                   stampId, i, 
                                                   info->getLocalIndexManager());
        
        info->addMolecule(mol);
        
#ifdef IS_MPI
        
      }
      
#endif
      
    } //end for(int i=0)   
  }
    
  int SimCreator::computeStorageLayout(SimInfo* info) {

    Globals* simParams = info->getSimParams();
    int nRigidBodies = info->getNGlobalRigidBodies();
    set<AtomType*> atomTypes = info->getSimulatedAtomTypes();
    set<AtomType*>::iterator i;
    bool hasDirectionalAtoms = false;
    bool hasFixedCharge = false;
    bool hasDipoles = false;    
    bool hasQuadrupoles = false;    
    bool hasPolarizable = false;    
    bool hasFluctuatingCharge = false;    
    bool hasMetallic = false;
    int storageLayout = 0;
    storageLayout |= DataStorage::dslPosition;
    storageLayout |= DataStorage::dslVelocity;
    storageLayout |= DataStorage::dslForce;

    for (i = atomTypes.begin(); i != atomTypes.end(); ++i) {

      DirectionalAdapter da = DirectionalAdapter( (*i) );
      MultipoleAdapter ma = MultipoleAdapter( (*i) );
      EAMAdapter ea = EAMAdapter( (*i) );
      SuttonChenAdapter sca = SuttonChenAdapter( (*i) );
      PolarizableAdapter pa = PolarizableAdapter( (*i) );
      FixedChargeAdapter fca = FixedChargeAdapter( (*i) );
      FluctuatingChargeAdapter fqa = FluctuatingChargeAdapter( (*i) );

      if (da.isDirectional()){
        hasDirectionalAtoms = true;
      }
      if (ma.isDipole()){
        hasDipoles = true;
      }
      if (ma.isQuadrupole()){
        hasQuadrupoles = true;
      }
      if (ea.isEAM() || sca.isSuttonChen()){
        hasMetallic = true;
      }
      if ( fca.isFixedCharge() ){
        hasFixedCharge = true;
      }
      if ( fqa.isFluctuatingCharge() ){
        hasFluctuatingCharge = true;
      }
      if ( pa.isPolarizable() ){
        hasPolarizable = true;
      }
    }
    
    if (nRigidBodies > 0 || hasDirectionalAtoms) {
      storageLayout |= DataStorage::dslAmat;
      if(storageLayout & DataStorage::dslVelocity) {
        storageLayout |= DataStorage::dslAngularMomentum;
      }
      if (storageLayout & DataStorage::dslForce) {
        storageLayout |= DataStorage::dslTorque;
      }
    }
    if (hasDipoles) {
      storageLayout |= DataStorage::dslDipole;
    }
    if (hasQuadrupoles) {
      storageLayout |= DataStorage::dslQuadrupole;
    }
    if (hasFixedCharge || hasFluctuatingCharge) {
      storageLayout |= DataStorage::dslSkippedCharge;
    }
    if (hasMetallic) {
      storageLayout |= DataStorage::dslDensity;
      storageLayout |= DataStorage::dslFunctional;
      storageLayout |= DataStorage::dslFunctionalDerivative;
    }
    if (hasPolarizable) {
      storageLayout |= DataStorage::dslElectricField;
    }
    if (hasFluctuatingCharge){
      storageLayout |= DataStorage::dslFlucQPosition;
      if(storageLayout & DataStorage::dslVelocity) {
        storageLayout |= DataStorage::dslFlucQVelocity;
      }
      if (storageLayout & DataStorage::dslForce) {
        storageLayout |= DataStorage::dslFlucQForce;
      }
    }
    
    // if the user has asked for them, make sure we've got the memory for the
    // objects defined.

    if (simParams->getOutputParticlePotential()) {
      storageLayout |= DataStorage::dslParticlePot;
    }

    if (simParams->havePrintHeatFlux()) {
      if (simParams->getPrintHeatFlux()) {
        storageLayout |= DataStorage::dslParticlePot;
      }
    }

    if (simParams->getOutputElectricField()) {
      storageLayout |= DataStorage::dslElectricField;
    }

    if (simParams->getOutputFluctuatingCharges()) {
      storageLayout |= DataStorage::dslFlucQPosition;
      storageLayout |= DataStorage::dslFlucQVelocity;
      storageLayout |= DataStorage::dslFlucQForce;
    }

    return storageLayout;
  }

  void SimCreator::setGlobalIndex(SimInfo *info) {
    SimInfo::MoleculeIterator mi;
    Molecule::AtomIterator ai;
    Molecule::RigidBodyIterator ri;
    Molecule::CutoffGroupIterator ci;
    Molecule::IntegrableObjectIterator  ioi;
    Molecule * mol;
    Atom * atom;
    RigidBody * rb;
    CutoffGroup * cg;
    int beginAtomIndex;
    int beginRigidBodyIndex;
    int beginCutoffGroupIndex;
    int nGlobalAtoms = info->getNGlobalAtoms();
    int nGlobalRigidBodies = info->getNGlobalRigidBodies();
    
    beginAtomIndex = 0;
    //rigidbody's index begins right after atom's
    beginRigidBodyIndex = info->getNGlobalAtoms();
    beginCutoffGroupIndex = 0;

    for(int i = 0; i < info->getNGlobalMolecules(); i++) {
      
#ifdef IS_MPI      
      if (info->getMolToProc(i) == worldRank) {
#endif        
        // stuff to do if I own this molecule
        mol = info->getMoleculeByGlobalIndex(i);

        //local index(index in DataStorge) of atom is important
        for(atom = mol->beginAtom(ai); atom != NULL; atom = mol->nextAtom(ai)) {
          atom->setGlobalIndex(beginAtomIndex++);
        }
        
        for(rb = mol->beginRigidBody(ri); rb != NULL;
            rb = mol->nextRigidBody(ri)) {
          rb->setGlobalIndex(beginRigidBodyIndex++);
        }
        
        //local index of cutoff group is trivial, it only depends on
        //the order of travesing
        for(cg = mol->beginCutoffGroup(ci); cg != NULL;
            cg = mol->nextCutoffGroup(ci)) {
          cg->setGlobalIndex(beginCutoffGroupIndex++);
        }        
        
#ifdef IS_MPI        
      }  else {

        // stuff to do if I don't own this molecule
        
        int stampId = info->getMoleculeStampId(i);
        MoleculeStamp* stamp = info->getMoleculeStamp(stampId);

        beginAtomIndex += stamp->getNAtoms();
        beginRigidBodyIndex += stamp->getNRigidBodies();
        beginCutoffGroupIndex += stamp->getNCutoffGroups() + stamp->getNFreeAtoms();
      }
#endif          

    } //end for(int i=0)  

    //fill globalGroupMembership
    std::vector<int> globalGroupMembership(info->getNGlobalAtoms(), 0);
    for(mol = info->beginMolecule(mi); mol != NULL; mol = info->nextMolecule(mi)) {        
      for (cg = mol->beginCutoffGroup(ci); cg != NULL; cg = mol->nextCutoffGroup(ci)) {
        
        for(atom = cg->beginAtom(ai); atom != NULL; atom = cg->nextAtom(ai)) {
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
    MPI::COMM_WORLD.Allreduce(&globalGroupMembership[0], 
                              &tmpGroupMembership[0], nGlobalAtoms,
                              MPI::INT, MPI::SUM);
    info->setGlobalGroupMembership(tmpGroupMembership);
#else
    info->setGlobalGroupMembership(globalGroupMembership);
#endif
    
    //fill molMembership
    std::vector<int> globalMolMembership(info->getNGlobalAtoms() + 
                                         info->getNGlobalRigidBodies(), 0);
    
    for(mol = info->beginMolecule(mi); mol != NULL; 
        mol = info->nextMolecule(mi)) {
      for(atom = mol->beginAtom(ai); atom != NULL; atom = mol->nextAtom(ai)) {
        globalMolMembership[atom->getGlobalIndex()] = mol->getGlobalIndex();
      }
      for (rb = mol->beginRigidBody(ri); rb != NULL; 
           rb = mol->nextRigidBody(ri)) {
        globalMolMembership[rb->getGlobalIndex()] = mol->getGlobalIndex();
      }
    }
    
#ifdef IS_MPI
    std::vector<int> tmpMolMembership(info->getNGlobalAtoms() + 
                                      info->getNGlobalRigidBodies(), 0);
    MPI::COMM_WORLD.Allreduce(&globalMolMembership[0], &tmpMolMembership[0], 
                              nGlobalAtoms + nGlobalRigidBodies,
                              MPI::INT, MPI::SUM);
    
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
    MPI::COMM_WORLD.Allreduce(&nIOPerMol[0], &numIntegrableObjectsPerMol[0], 
                              info->getNGlobalMolecules(), MPI::INT, MPI::SUM);
#else
    std::vector<int> numIntegrableObjectsPerMol = nIOPerMol;
#endif    

    std::vector<int> startingIOIndexForMol(info->getNGlobalMolecules());
    
    int startingIndex = 0;
    for (int i = 0; i < info->getNGlobalMolecules(); i++) {
      startingIOIndexForMol[i] = startingIndex;
      startingIndex += numIntegrableObjectsPerMol[i];
    }
    
    std::vector<StuntDouble*> IOIndexToIntegrableObject(info->getNGlobalIntegrableObjects(), (StuntDouble*)NULL);
    for (mol = info->beginMolecule(mi); mol != NULL; 
         mol = info->nextMolecule(mi)) {
      int myGlobalIndex = mol->getGlobalIndex();
      int globalIO = startingIOIndexForMol[myGlobalIndex];
      for (StuntDouble* sd = mol->beginIntegrableObject(ioi); sd != NULL;
           sd = mol->nextIntegrableObject(ioi)) {
        sd->setGlobalIntegrableObjectIndex(globalIO);
        IOIndexToIntegrableObject[globalIO] = sd;
        globalIO++;
      }
    }
      
    info->setIOIndexToIntegrableObject(IOIndexToIntegrableObject);
    
  }
  
  void SimCreator::loadCoordinates(SimInfo* info, const std::string& mdFileName) {
    Globals* simParams;

    simParams = info->getSimParams();
    
    DumpReader reader(info, mdFileName);
    int nframes = reader.getNFrames();

    if (nframes > 0) {
      reader.readFrame(nframes - 1);
    } else {
      //invalid initial coordinate file
      sprintf(painCave.errMsg, 
              "Initial configuration file %s should at least contain one frame\n",
              mdFileName.c_str());
      painCave.isFatal = 1;
      simError();
    }
    //copy the current snapshot to previous snapshot
    info->getSnapshotManager()->advance();
  }
  
} //end namespace OpenMD


