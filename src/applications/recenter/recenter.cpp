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

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iostream>
#include <map>
#include <string>

#include "applications/recenter/recenterCmd.hpp"
#include "brains/Register.hpp"
#include "brains/SimCreator.hpp"
#include "brains/SimInfo.hpp"
#include "brains/Thermo.hpp"
#include "io/DumpReader.hpp"
#include "io/DumpWriter.hpp"
#include "utils/StringUtils.hpp"

using namespace std;
using namespace OpenMD;

int main(int argc, char* argv[]) {
  registerLattice();

  gengetopt_args_info args_info;

  std::string inputFileName;
  std::string outputFileName;

  // parse command line arguments
  if (cmdline_parser(argc, argv, &args_info) != 0) exit(1);

  // get input file name
  if (args_info.inputs_num)
    inputFileName = args_info.inputs[0];
  else {
    snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
             "No input file name was specified "
             "on the command line");
    painCave.severity = OPENMD_ERROR;
    painCave.isFatal  = 1;
    simError();
  }

  // parse md file and set up the system

  SimCreator creator;
  SimInfo* info = creator.createSim(inputFileName, false);
  DumpReader reader(info, inputFileName);
  // very important step:
  info->update();

  outputFileName = args_info.output_arg;

  if (!outputFileName.compare(inputFileName)) {
    snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
             "Input and Output File names should be different!");
    painCave.severity = OPENMD_ERROR;
    painCave.isFatal  = 1;
    simError();
  }

  DumpWriter* writer = new DumpWriter(info, outputFileName);

  if (writer == NULL) {
    snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
             "error in creating DumpWriter");
    painCave.severity = OPENMD_ERROR;
    painCave.isFatal  = 1;
    simError();
  }

  int nFrames = reader.getNFrames();
  Vector3d COM;
  Vector3d pos;
  SimInfo::MoleculeIterator i;
  Molecule::IntegrableObjectIterator j;
  Molecule* mol;
  StuntDouble* sd;
  Thermo thermo(info);

  for (int istep = 0; istep < nFrames; istep++) {
    reader.readFrame(istep);
    COM = thermo.getCom();
    for (mol = info->beginMolecule(i); mol != NULL;
         mol = info->nextMolecule(i)) {
      for (sd = mol->beginIntegrableObject(j); sd != NULL;
           sd = mol->nextIntegrableObject(j)) {
        pos = sd->getPos();
        sd->setPos(pos - COM);
      }
    }
    writer->writeDump();
  }

  // deleting the writer will put the closing at the end of the dump file.

  delete writer;

  snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
           "A new OpenMD file called \"%s\" has been generated.\n",
           outputFileName.c_str());
  painCave.severity = OPENMD_INFO;
  painCave.isFatal  = 0;
  simError();

  return 0;
}
