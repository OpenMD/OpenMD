/*
 * Copyright (c) 2004-2020 The University of Notre Dame. All Rights Reserved.
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
 * [1] Meineke, et al., J. Comp. Chem. 26, 252-271 (2005).
 * [2] Fennell & Gezelter, J. Chem. Phys. 124, 234104 (2006).
 * [3] Sun, Lin & Gezelter, J. Chem. Phys. 128, 234107 (2008).
 * [4] Vardeman, Stocker & Gezelter, J. Chem. Theory Comput. 7, 834 (2011).
 * [5] Kuang & Gezelter, Mol. Phys., 110, 691-701 (2012).
 * [6] Lamichhane, Gezelter & Newman, J. Chem. Phys. 141, 134109 (2014).
 * [7] Lamichhane, Newman & Gezelter, J. Chem. Phys. 141, 134110 (2014).
 * [8] Bhattarai, Newman & Gezelter, Phys. Rev. B 99, 094106 (2019).
 */

#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <cmath>
#include <iostream>
#include <string>
#include <map>
#include <fstream>

#include "thermalizerCmd.hpp"
#include "brains/Velocitizer.hpp"
#include "brains/Register.hpp"
#include "brains/SimInfo.hpp"
#include "brains/SimCreator.hpp"
#include "brains/Thermo.hpp"
#include "brains/ForceManager.hpp"
#include "io/DumpReader.hpp"
#include "io/DumpWriter.hpp"
#include "utils/StringUtils.hpp"
#include "utils/MemoryUtils.hpp"

using namespace OpenMD;

int main(int argc, char *argv []) {

  gengetopt_args_info args_info;
  std::string inputFileName;
  std::string outputFileName;

  // parse command line arguments
  if (cmdline_parser(argc, argv, &args_info) != 0) {
    cmdline_parser_print_help();
    exit(1);
  }

  // get input file name
  if (args_info.input_given) {
    inputFileName = args_info.input_arg;
  } else {
    if (args_info.inputs_num){
      inputFileName = args_info.inputs[0];
    }
    else{
      sprintf(painCave.errMsg,
              "No input file name was specified on the command line");
      painCave.severity = OPENMD_ERROR;
      painCave.isFatal = 1;
      simError();
    }
  }

  // get output file name
  outputFileName = args_info.output_arg;

  if (!outputFileName.compare(inputFileName)) {
    sprintf(painCave.errMsg,
            "Input and Output File names should be different!");
    painCave.severity = OPENMD_ERROR;
    painCave.isFatal = 1;
    simError();
  }

  // Parse the input file, set up the system, and read the last frame:
  SimCreator creator;
  SimInfo* info = creator.createSim(inputFileName, true);
  // communicate velocity information onto the atoms:
  info->update();

  // Important utility classes for computing system properties:
  Thermo thermo(info);

  // Remove in favor of std::MemoryUtils::make_unique<> when we switch to C++14 and above
  VelocitizerPtr veloSet {MemoryUtils::make_unique<Velocitizer>(info)};
  
  ForceManager* forceMan = new ForceManager(info);

  // Just in case we were passed a system that is on the move:
  veloSet->removeComDrift();
  forceMan->calcForces();

  RealType instPE = thermo.getPotential();
  RealType instKE = thermo.getKinetic();

  // Now that we have the information from the current frame, advance
  // the snapshot to make a modified frame:
  info->getSnapshotManager()->advance();

  // Create DumpWriter to hold the modified frame:
  DumpWriter* writer = new DumpWriter(info, outputFileName);
  if (writer == NULL) {
    sprintf(painCave.errMsg, "error in creating DumpWriter");
    painCave.severity = OPENMD_ERROR;
    painCave.isFatal = 1;
    simError();
  }

  // If resampling temperature, we call the randomizer method:
  if (args_info.temperature_given) {
    RealType temperature = args_info.temperature_arg;

    if (temperature < 0.0) {
      sprintf(painCave.errMsg, "Temperatures must be positive numbers.");
      painCave.severity = OPENMD_ERROR;
      painCave.isFatal = 1;
      simError();
    }

    veloSet->randomize(temperature);
  }

  // If resampling charge temperature, we call the randomizeChargeVelocity method
  if (args_info.chargetemperature_given) {
    RealType charge_temperature = args_info.chargetemperature_arg;

    if (charge_temperature < 0.0) {
      sprintf(painCave.errMsg, "Temperatures must be positive numbers.");
      painCave.severity = OPENMD_ERROR;
      painCave.isFatal = 1;
      simError();
    }

    veloSet->randomizeChargeVelocity(charge_temperature);
  }

  // If scaling total energy, scale only the kinetic:
  if (args_info.energy_given) {
    RealType energy = args_info.energy_arg;
    RealType epsilon = 1e-6;
    RealType lambda = 0.0;

    if (energy < instPE) {
      sprintf(painCave.errMsg,
              "Energy must be larger than current potential energy.");
      painCave.severity = OPENMD_ERROR;
      painCave.isFatal = 1;
      simError();
    }
    else {
      if (instKE >= epsilon) {
	lambda =  sqrt((energy - instPE) / instKE);
        veloSet->scale(lambda);
      }
      // If the current kinetic energy is close to zero, we will
      // sample velocities from a 10K distribution and then
      // subsequently scale from 10K to the desired energy.
      else {
        veloSet->randomize(10.0);
        instKE = thermo.getKinetic();
	lambda = sqrt( (energy - instPE) / instKE );
        veloSet->scale(lambda);
      }
    }
  }

  writer->writeDump();
  // deleting the writer will put the closing at the end of the dump file.
  delete writer;

  sprintf(painCave.errMsg,
          "A new OpenMD file called \"%s\" has been generated.\n",
          outputFileName.c_str());
  painCave.isFatal = 0;
  painCave.severity = OPENMD_INFO;
  simError();
  return 0;
}
