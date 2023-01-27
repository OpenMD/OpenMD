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

#include <cstdlib>
#include <iostream>

#ifdef _MSC_VER
#include <Windows.h>
#include <io.h>
#include <stdio.h>
#define isatty _isatty
#define fileno _fileno
#else
#include <unistd.h>

#include <cstdio>

#include <sys/ioctl.h>
#endif

#ifdef IS_MPI
#include <mpi.h>
#endif

#include "utils/ProgressBar.hpp"

using namespace std;

namespace OpenMD {

  const char* progressSpinner_ = "|/-\\";

  ProgressBar::ProgressBar() :
      value_(0.0), maximum_(-1.0), iteration_(0), start_(time(NULL)) {}

  void ProgressBar::clear() {
#ifdef IS_MPI
    int myRank;
    MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
    if (myRank == 0) {
#endif
      cout << endl;
      cout.flush();
#ifdef IS_MPI
    }
#endif
    iteration_ = 0;
    value_     = 0;
    maximum_   = -1;
    start_     = time(NULL);
  }

  void ProgressBar::update() {
#ifdef IS_MPI
    int myRank;
    MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
    if (myRank == 0) {
#endif

      // only do the progress bar if we are actually running in a tty:
      if (isatty(fileno(stdout)) && (getenv("SGE_TASK_ID") == NULL)) {
        // get the window width:

        int width = 0;
#ifdef _MSC_VER
        CONSOLE_SCREEN_BUFFER_INFO csbi;
        HANDLE hConsole = GetStdHandle(STD_OUTPUT_HANDLE);
        int ret         = GetConsoleScreenBufferInfo(hConsole, &csbi);
        if (ret) { width = csbi.dwSize.X - 1; }
#else
      struct winsize w;
      ioctl(fileno(stdout), TIOCGWINSZ, &w);
      width = w.ws_col;
#endif

        // handle the case when the width is returned as a nonsensical value.
        if (width <= 0) width = 80;

        // We'll use:
        // 31 characters for the completion estimate,
        //  6 for the % complete,
        //  2 characters for the open and closing brackets.

        int avail = width - 31 - 6 - 2;

        ++iteration_;

        if (maximum_ > 0.0) {
          // we know the maximum, so draw a progress bar

          RealType percent =
              std::min(std::max(value_ * 100.0 / maximum_, 1e-6), 100.0);

          int hashes = int(percent * avail / 100.0);

          // compute the best estimate of the ending time:
          time_t current_  = time(NULL);
          time_t end_      = static_cast<time_t>(start_ + (current_ - start_) *
                                                         (100.0 / percent));
          struct tm* ender = localtime(&end_);
          char buffer[22];
          strftime(buffer, 22, "%a %b %d @ %I:%M %p", ender);

#ifdef _MSC_VER
          csbi.dwCursorPosition.X = 0;
          SetConsoleCursorPosition(hConsole, csbi.dwCursorPosition);
#else
        cout << '\r';
#endif
          cout.width(3);
          cout << right << int(percent);
          cout.width(3);
          cout << "% [";
          cout.fill('#');
          if (hashes + 1 < avail) {
            cout.width(hashes + 1);
            cout << progressSpinner_[iteration_ & 3];
          } else {
            cout.width(avail);
            cout << '#';
          }
          cout.fill(' ');
          if (avail - hashes - 1 > 0) {
            cout.width(avail - hashes - 1);
            cout << ' ';
          }
          cout.width(11);
          cout << "] Estimate:";
          cout.width(22);
          cout << buffer;
        }
        cout.flush();
      }
#ifdef IS_MPI
    }
#endif
  }

  void ProgressBar::setStatus(RealType val, RealType max) {
    value_   = val;
    maximum_ = max;
  }
}  // namespace OpenMD
