/*
 * Copyright (c) 2006 The University of Notre Dame. All Rights Reserved.
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
 * [3]  Sun, Lin & Gezelter, J. Chem. Phys. 128, 234107 (2008).          
 * [4]  Kuang & Gezelter,  J. Chem. Phys. 133, 164101 (2010).
 * [5]  Vardeman, Stocker & Gezelter, J. Chem. Theory Comput. 7, 834 (2011).
 */

/**
 * @file thermalizer.cpp
 * @author Dan Gezelter
 * @date 09/20/2006
 * @version 1.0
 */

#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <cmath>
#include <iostream>
#include <string>
#include <map>
#include <fstream>

#include "applications/thermalizer/thermalizerCmd.h"
#include "integrators/Velocitizer.hpp"
#include "brains/Register.hpp"
#include "brains/SimInfo.hpp"
#include "brains/SimCreator.hpp"
#include "io/DumpReader.hpp"
#include "io/DumpWriter.hpp"
#include "utils/StringUtils.hpp"

using namespace std;
using namespace OpenMD;

int main(int argc, char *argv []) {

  registerLattice();
    
  gengetopt_args_info args_info;
  std::string inputFileName;
  std::string outputFileName;
  RealType temperature;

  // parse command line arguments
  if (cmdline_parser(argc, argv, &args_info) != 0)
    exit(1);

  temperature = args_info.temperature_arg;
    
  if (temperature < 0.0) {
    sprintf(painCave.errMsg, "Temperatures must be positive numbers.");
    painCave.isFatal = 1;
    simError();
  }

  //get input file name
  if (args_info.inputs_num)
    inputFileName = args_info.inputs[0];
  else {
    sprintf(painCave.errMsg, "No input file name was specified "
            "on the command line");
    painCave.isFatal = 1;
    simError();
  }
  
  //parse md file and set up the system

  SimCreator creator;
  SimInfo* info = creator.createSim(inputFileName, false);
  Velocitizer* veloSet = new Velocitizer(info);
  DumpReader reader(info, inputFileName);
  // very important step:
  info->update();

  outputFileName = args_info.output_arg;

  if (!outputFileName.compare(inputFileName)) {
    sprintf(painCave.errMsg, "Input and Output File names should be different!");
    painCave.isFatal = 1;
    simError();
  }

  DumpWriter* writer = new DumpWriter(info, outputFileName);
    
  if (writer == NULL) {
    sprintf(painCave.errMsg, "error in creating DumpWriter");
    painCave.isFatal = 1;
    simError();
  }

  int nFrames = reader.getNFrames();

  for (int istep = 0; istep < nFrames; istep++) {
    reader.readFrame(istep);
    veloSet->velocitize(temperature);
    writer->writeDump();
  }

  // deleting the writer will put the closing at the end of the dump file.

  delete writer;

  sprintf(painCave.errMsg, "A new OpenMD file called \"%s\" has been "
          "generated.\n", outputFileName.c_str());
  painCave.isFatal = 0;
  simError();
  return 0;
}
