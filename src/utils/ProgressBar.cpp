/*
 * Copyright (c) 2010 The University of Notre Dame. All Rights Reserved.
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
 * [4]  Vardeman, Stocker & Gezelter, in progress (2010).                        
 */

#include <string>
#include <sstream>
#include <sys/ioctl.h>
#define _POSIX_SOURCE
#include <unistd.h>
#ifdef IS_MPI
#include <mpi.h>
#endif //is_mpi
#include "utils/ProgressBar.hpp"

namespace OpenMD {

  const char * progressSpinner_ = "|/-\\";

  ProgressBar::ProgressBar() : value_(0.0), maximum_(-1.0), iteration_(0), start_(time(NULL)) {
  }
  
  void ProgressBar::clear() {
#ifdef IS_MPI
    if (MPI::COMM_WORLD.Get_rank() == 0) {
#endif
      printf("\n");
      fflush(stdout);
#ifdef IS_MPI
    }
#endif
    iteration_ = 0;
    value_ = 0;
    maximum_ = -1;
    start_ = time(NULL);
  }
  
  void ProgressBar::update() {
    int width = 80;
#ifdef IS_MPI
    if (MPI::COMM_WORLD.Get_rank() == 0) {
#endif

      // only do the progress bar if we are actually running in a tty:
      if (isatty(fileno(stdout))) {        
        // get the window width:
        struct winsize w;
        ioctl(fileno(stdout), TIOCGWINSZ, &w);
        width = w.ws_col;

        // We'll use:
        // 31 characters for the completion estimate, 
        //  6 for the % complete, 
        //  2 characters for the open and closing brackets.
        
        int avail = width - 31 - 6 - 2;
        
        ++iteration_;
        
        if (maximum_ > 0.0) {

          // we know the maximum, so draw a progress bar
          
          RealType percent = value_ * 100.0 / maximum_;
          int hashes = int(percent * avail / 100.0);
          std::string progressbar;
          progressbar.assign(hashes, '#');
          
          // add the spinner to the end of the progress bar:
          progressbar += progressSpinner_[iteration_ & 3];
          
          // compute the best estimate of the ending time:
          time_t current_ = time(NULL);
          time_t end_ = start_ + (current_ - start_) * (100.0/percent);
          struct tm * ender = localtime(&end_);
          char buffer[24];
          strftime(buffer, 24, "%a %b %d @ %I:%M %p", ender);
          
          std::stringstream fmt;
          fmt << "\r%3d%% [%-" << avail << "s] Estimate: %s";
          std::string st = fmt.str();
          
          printf(st.c_str(), int(percent),
                 progressbar.c_str(),
                 buffer);
          
        } else {
          // we don't know the maximum, so we can't draw a progress bar
          int center = (iteration_ % 48) + 1; // 50 spaces, minus 2
          std::string before;
          std::string after;
          before.assign(std::max(center - 2, 0), ' ');
          after.assign(std::min(center + 2, 50), ' ');
          
          printf("\r[%s###%s]            ",
                 before.c_str(), after.c_str());
        }
        fflush(stdout);
      }
#ifdef IS_MPI
    }
#endif
  }
  
  void ProgressBar::setStatus(RealType val, RealType max) {
    value_ = val;
    maximum_ = max;
  }  
}
