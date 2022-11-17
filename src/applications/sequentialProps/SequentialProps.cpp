/*
 * Copyright (c) 2004-2022, The University of Notre Dame. All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 * 1. Redistributions of source code must retain the above copyright notice,
 * this list of conditions and the following disclaimer.
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

#include <fstream>
#include <iostream>
#include <string>

#include "SequentialPropsCmd.hpp"
#include "applications/sequentialProps/COMVel.hpp"
#include "applications/sequentialProps/CenterOfMass.hpp"
#include "applications/sequentialProps/ContactAngle1.hpp"
#include "applications/sequentialProps/ContactAngle2.hpp"
#include "applications/sequentialProps/GCNSeq.hpp"
#include "applications/sequentialProps/SequentialAnalyzer.hpp"
#include "applications/sequentialProps/equipartitionTest.hpp"
#include "brains/SimCreator.hpp"
#include "brains/SimInfo.hpp"
#include "utils/StringUtils.hpp"
#include "utils/simError.h"

using namespace OpenMD;

int main(int argc, char* argv[]) {
  gengetopt_args_info args_info;

  // parse the command line option
  if (cmdline_parser(argc, argv, &args_info) != 0) { exit(1); }

  // get the dumpfile name and meta-data file name
  std::string dumpFileName = args_info.input_arg;

  std::string sele1;
  std::string sele2;

  // check the first selection argument, or set it to the environment
  // variable, or failing that, set it to "select all"

  if (args_info.sele1_given) {
    sele1 = args_info.sele1_arg;
  } else {
    char* sele1Env = getenv("SELECTION1");
    if (sele1Env) {
      sele1 = sele1Env;
    } else {
      sele1 = "select all";
    }
  }

  // check the second selection argument, or set it to the environment
  // variable, or failing that, set it to the first selection

  if (args_info.sele2_given) {
    sele2 = args_info.sele2_arg;
  } else {
    char* sele2Env = getenv("SELECTION2");
    if (sele2Env) {
      sele2 = sele2Env;
    } else {
      // If sele2 is not specified, then the default behavior
      // should be what is already intended for sele1
      sele2 = sele1;
    }
  }

  // parse md file and set up the system
  SimCreator creator;
  SimInfo* info = creator.createSim(dumpFileName, false);

  SequentialAnalyzer* analyzer = NULL;
  if (args_info.com_given) {
    analyzer = new CenterOfMass(info, dumpFileName, sele1, sele2);
  } else if (args_info.comvel_given) {
    analyzer = new COMVel(info, dumpFileName, sele1, sele2);
  } else if (args_info.testequi_given) {
    analyzer = new Equipartition(info, dumpFileName, sele1, sele2);
  } else if (args_info.gcn_given) {
    if (args_info.rcut_given) {
      analyzer = new GCNSeq(info, dumpFileName, sele1, sele2,
                            args_info.rcut_arg, args_info.nbins_arg);
    } else {
      snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
               "A cutoff radius (rcut) must be specified when calculating\n"
               "\tGeneralized Coordinate Number");
      painCave.severity = OPENMD_ERROR;
      painCave.isFatal  = 1;
      simError();
    }
  } else if (args_info.ca1_given) {
    RealType solidZ(0.0);
    if (args_info.referenceZ_given)
      solidZ = args_info.referenceZ_arg;
    else {
      snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
               "--referenceZ must be set if --ca1 is used\n");
      painCave.severity = OPENMD_ERROR;
      painCave.isFatal  = 1;
      simError();
    }
    RealType dropletR(0.0);
    if (args_info.dropletR_given)
      dropletR = args_info.dropletR_arg;
    else {
      snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
               "--dropletR must be set if --ca1 is used\n");
      painCave.severity = OPENMD_ERROR;
      painCave.isFatal  = 1;
      simError();
    }

    analyzer =
        new ContactAngle1(info, dumpFileName, sele1, sele2, solidZ, dropletR);
  } else if (args_info.ca2_given) {
    RealType solidZ(0.0);
    if (args_info.referenceZ_given)
      solidZ = args_info.referenceZ_arg;
    else {
      snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
               "--referenceZ must be set if --ca2 is used\n");
      painCave.severity = OPENMD_ERROR;
      painCave.isFatal  = 1;
      simError();
    }
    RealType centroidX(0.0);
    if (args_info.centroidX_given)
      centroidX = args_info.centroidX_arg;
    else {
      snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
               "--centroidX must be set if --ca2 is used\n");
      painCave.severity = OPENMD_ERROR;
      painCave.isFatal  = 1;
      simError();
    }
    RealType centroidY(0.0);
    if (args_info.centroidY_given)
      centroidY = args_info.centroidY_arg;
    else {
      snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
               "--centroidY must be set if --ca2 is used\n");
      painCave.severity = OPENMD_ERROR;
      painCave.isFatal  = 1;
      simError();
    }
    RealType threshDens(0.0);
    if (args_info.threshDens_given)
      threshDens = args_info.threshDens_arg;
    else {
      snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
               "--threshDens must be set if --ca2 is used\n");
      painCave.severity = OPENMD_ERROR;
      painCave.isFatal  = 1;
      simError();
    }
    RealType bufferLength(0.0);
    if (args_info.bufferLength_given)
      bufferLength = args_info.bufferLength_arg;
    else {
      snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
               "--bufferLength must be set if --ca2 is used\n");
      painCave.severity = OPENMD_ERROR;
      painCave.isFatal  = 1;
      simError();
    }

    analyzer = new ContactAngle2(info, dumpFileName, sele1, sele2, solidZ,
                                 centroidX, centroidY, threshDens, bufferLength,
                                 args_info.nbins_arg, args_info.nbins_z_arg);
  }

  if (args_info.output_given) { analyzer->setOutputName(args_info.output_arg); }

  analyzer->doSequence();

  delete analyzer;
  delete info;

  return 0;
}
