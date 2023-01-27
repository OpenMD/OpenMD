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

#include "utils/simError.h"

#include <config.h>

#include <cstdio>
#include <cstdlib>
#include <cstring>

#ifdef IS_MPI
#include <mpi.h>
#endif

int nChecks;

errorStruct painCave;

char checkPointMsg[MAX_SIM_ERROR_MSG_LENGTH];
int worldRank;

void initSimError(void) {
  painCave.errMsg[0]   = '\0';
  painCave.isFatal     = 0;
  painCave.severity    = OPENMD_ERROR;
  painCave.isEventLoop = 0;
  nChecks              = 0;
#ifdef IS_MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &worldRank);
#else
  worldRank = 0;
#endif
}

int simError(void) {
  char errorMsg[MAX_SIM_ERROR_MSG_LENGTH];

#ifdef IS_MPI
  int myError = 1;
  int isError;
  char nodeMsg[MAX_SIM_ERROR_MSG_LENGTH];
#endif

  strcpy(errorMsg, "OpenMD ");
  switch (painCave.severity) {
  case OPENMD_WARNING:
    strcat(errorMsg, "warning");
    break;
  case OPENMD_INFO:
    strcat(errorMsg, "info");
    break;
  default:
    if (painCave.isFatal) { strcat(errorMsg, "FATAL "); }
    strcat(errorMsg, "ERROR");
  }

#ifdef IS_MPI
  if (worldRank == 0) {
    if (painCave.isEventLoop) {
      snprintf(nodeMsg, MAX_SIM_ERROR_MSG_LENGTH, " (reported by MPI node %d)",
               worldRank);
      strncat(errorMsg, nodeMsg,
              MAX_SIM_ERROR_MSG_LENGTH - strlen(errorMsg) - 1);
      errorMsg[MAX_SIM_ERROR_MSG_LENGTH - 1] = '\0';
    }
#endif

    strcat(errorMsg, ":\n\t");
    strncat(errorMsg, painCave.errMsg,
            MAX_SIM_ERROR_MSG_LENGTH - strlen(errorMsg) - 1);
    errorMsg[MAX_SIM_ERROR_MSG_LENGTH - 1] = '\0';
    strcat(errorMsg, "\n");

    switch (painCave.severity) {
    case OPENMD_WARNING:
    case OPENMD_INFO:
      fprintf(stdout, "%s", errorMsg);
      break;
    default:
      fprintf(stderr, "%s", errorMsg);
    }

#ifdef IS_MPI
    if (painCave.isEventLoop) return 1;
  }
#endif

  if (painCave.isFatal) {
#ifdef IS_MPI
    MPI_Allreduce(&myError, &isError, 1, MPI_INT, MPI_LOR, MPI_COMM_WORLD);
    MPI_Finalize();
#endif
    exit(0);
  }
  return 1;
}

void errorCheckPoint(void) {
  int myError = 0;
  int isError = 0;

#ifdef IS_MPI
  MPI_Allreduce(&myError, &isError, 1, MPI_INT, MPI_LOR, MPI_COMM_WORLD);
#else
  isError   = myError;
#endif

  if (isError) {
#ifdef IS_MPI
    MPI_Finalize();
#endif
    exit(0);
  }

#ifdef CHECKPOINT_VERBOSE
  nChecks++;

#ifdef IS_MPI
  if (worldRank == 0) {
#endif

    fprintf(stderr, "Checkpoint #%d reached: %s\n", nChecks, checkPointMsg);
#ifdef IS_MPI
  }
#endif

#endif
}
