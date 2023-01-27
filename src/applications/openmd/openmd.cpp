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

#include <fstream>
#include <iostream>
#include <locale>

#ifdef IS_MPI
#include <mpi.h>
#endif

#include "brains/Register.hpp"
#include "brains/SimCreator.hpp"
#include "brains/SimInfo.hpp"
#include "integrators/Integrator.hpp"
#include "integrators/IntegratorFactory.hpp"
#include "optimization/Constraint.hpp"
#include "optimization/Method.hpp"
#include "optimization/OptimizationFactory.hpp"
#include "optimization/PotentialEnergyObjectiveFunction.hpp"
#include "optimization/Problem.hpp"
#include "utils/CaseConversion.hpp"
#include "utils/Revision.hpp"
#include "utils/simError.h"

using namespace OpenMD;
using namespace QuantLib;

int main(int argc, char* argv[]) {
  // first things first, all of the initializations

#ifdef IS_MPI
  MPI_Init(&argc, &argv);  // the MPI communicators
#endif

  initSimError();  // the error handler

  Revision r;

#ifdef IS_MPI
  if (worldRank == 0) {
#endif
    std::cout
        << "  "
           "+------------------------------------------------------------------"
           "----"
           "----+\n"
        << "  |    ____                    __  ___ ____                        "
           "    "
           "      |\n"
        << "  |   / __ \\____  ___  ____   /  |/  // __ \\  The Open Molecular "
           "Dynamics    |\n"
        << "  |  / / / / __ \\/ _ \\/ __ \\ / /|_/ // / / /  Engine (formerly "
           "OOPSE).       |\n"
        << "  | / /_/ / /_/ /  __/ / / // /  / // /_/ /                        "
           "    "
           "      |\n"
        << "  | \\____/ .___/\\___/_/ /_//_/  /_//_____/    Copyright "
           "2004-2023 by "
           "the     |\n"
        << "  |     /_/                                   University of Notre "
           "Dame.      |\n"
        << "  |           http://openmd.org                                    "
           "    "
           "      |\n"
        << "  |                                                                "
           "    "
           "      |\n"
        << "  |   " << r.getFullRevision() << "       |\n"
        << "  |               " << r.getBuildDate()
        << "                       |\n"
        << "  |                                                                "
           "    "
           "      |\n"
        << "  | OpenMD is an OpenScience project.  All source code is "
           "available "
           "for any  |\n"
        << "  | use whatsoever under a BSD-style license.                      "
           "    "
           "      |\n"
        << "  |                                                                "
           "    "
           "      |\n"
        << "  | Support OpenScience!  If you use OpenMD or its source code in "
           "your "
           "      |\n"
        << "  | research, please cite the appropriate papers when you publish "
           "your "
           "work. |\n"
        << "  | Good starting points for code and simulation methodology are:  "
           "    "
           "      |\n"
        << "  |                                                                "
           "    "
           "      |\n"
        << "  | [1] Meineke, et al., J. Comp. Chem. 26, 252-271 (2005).        "
           "    "
           "      |\n"
        << "  | [2] Fennell & Gezelter, J. Chem. Phys. 124, 234104 (2006).     "
           "    "
           "      |\n"
        << "  | [3] Sun, Lin & Gezelter, J. Chem. Phys. 128, 234107 (2008).    "
           "    "
           "      |\n"
        << "  | [4] Vardeman, Stocker & Gezelter, J. Chem. Theory Comput. 7, "
           "834 "
           "(2011). |\n"
        << "  | [5] Kuang & Gezelter, Mol. Phys., 110, 691-701 (2012).         "
           "    "
           "      |\n"
        << "  | [6] Lamichhane, Gezelter & Newman, J. Chem. Phys. 141, 134109 "
           "(2014).    |\n"
        << "  | [7] Lamichhane, Newman & Gezelter, J. Chem. Phys. 141, 134110 "
           "(2014).    |\n"
        << "  | [8] Bhattarai, Newman & Gezelter, Phys. Rev. B 99, 094106 "
           "(2019).  "
           "      |\n"
        << "  "
           "+------------------------------------------------------------------"
           "----"
           "----+\n"
        << "\n";

    if (argc < 2) {
      strcpy(painCave.errMsg,
             "No meta-data file was specified on the command line.\n");
      painCave.isFatal = 1;
      simError();
    }
#ifdef IS_MPI
  }
#endif

  strcpy(checkPointMsg, "Successful number of arguments");
  errorCheckPoint();

  // register forcefields, integrators and minimizers
  registerAll();

  // create simulation model
  SimCreator creator;
  SimInfo* info = creator.createSim(argv[1]);

  Globals* simParams            = info->getSimParams();
  MinimizerParameters* miniPars = simParams->getMinimizerParameters();

  if (miniPars->getUseMinimizer() && simParams->haveEnsemble()) {
    snprintf(
        painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
        "Ensemble keyword can not co-exist with useMinimizer = \"true\" in the "
        "minimizer block\n");
    painCave.isFatal = 1;
    simError();
  }

  if (miniPars->getUseMinimizer()) {
    // create minimizer
    OptimizationMethod* myMinimizer =
        OptimizationFactory::getInstance().createOptimization(
            toUpperCopy(miniPars->getMethod()), info);

    if (myMinimizer == NULL) {
      snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
               "Optimization Factory can not create %s OptimizationMethod\n",
               miniPars->getMethod().c_str());
      painCave.isFatal = 1;
      simError();
    }

    ForceManager* fman = new ForceManager(info);
    fman->initialize();

    PotentialEnergyObjectiveFunction potObjf(info, fman);
    NoConstraint noConstraint {};
    DumpStatusFunction dsf(info);
    DynamicVector<RealType> initCoords = potObjf.setInitialCoords();
    Problem problem(potObjf, noConstraint, dsf, initCoords);

    int maxIter              = miniPars->getMaxIterations();
    int mssIter              = miniPars->getMaxStationaryStateIterations();
    RealType rEps            = miniPars->getRootEpsilon();
    RealType fEps            = miniPars->getFunctionEpsilon();
    RealType gnEps           = miniPars->getGradientNormEpsilon();
    RealType initialStepSize = miniPars->getInitialStepSize();

    EndCriteria endCriteria(maxIter, mssIter, rEps, fEps, gnEps);
    myMinimizer->minimize(problem, endCriteria, initialStepSize);

    delete myMinimizer;
  } else if (simParams->haveEnsemble()) {
    // create Integrator
    Integrator* myIntegrator =
        IntegratorFactory::getInstance().createIntegrator(
            toUpperCopy(simParams->getEnsemble()), info);

    if (myIntegrator == NULL) {
      snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
               "Integrator Factory can not create %s Integrator\n",
               simParams->getEnsemble().c_str());
      painCave.isFatal = 1;
      simError();
    }

    myIntegrator->integrate();
    delete myIntegrator;
  } else {
    snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
             "Integrator Factory can not create %s Integrator\n",
             simParams->getEnsemble().c_str());
    painCave.isFatal = 1;
    simError();
  }

  delete info;

  strcpy(checkPointMsg, "Great googly moogly!  It worked!");
  errorCheckPoint();

#ifdef IS_MPI
  MPI_Finalize();
#endif

  return 0;
}
