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
 * [4]  Kuang & Gezelter,  J. Chem. Phys. 133, 164101 (2010).
 * [5]  Vardeman, Stocker & Gezelter, J. Chem. Theory Comput. 7, 834 (2011).
 */

#include <iostream>
#include <cstdlib>

#ifdef _MSC_VER
#include <Windows.h>
#include <stdio.h>
#include <io.h>
#define isatty _isatty
#define fileno _fileno
#else
#include <cstdio>
#include <sys/ioctl.h>
#include <unistd.h>
#endif

#ifdef IS_MPI
#include <mpi.h>
#endif

#include "utils/ProgressBar.hpp"

using namespace std;

namespace OpenMD {

  const char * progressSpinner_ = "|/-\\";

  ProgressBar::ProgressBar() : value_(0.0), maximum_(-1.0), iteration_(0), start_(time(NULL)) {
  }
  
  void ProgressBar::clear() {
#ifdef IS_MPI
    if (MPI::COMM_WORLD.Get_rank() == 0) {
#endif
      cout << endl;
      cout.flush();
#ifdef IS_MPI
    }
#endif
    iteration_ = 0;
    value_ = 0;
    maximum_ = -1;
    start_ = time(NULL);
  }
  
  void ProgressBar::update() {

    int width;

#ifdef IS_MPI
    if (MPI::COMM_WORLD.Get_rank() == 0) {
#endif
      
      // only do the progress bar if we are actually running in a tty:
      if (isatty(fileno(stdout))  && (getenv("SGE_TASK_ID")==NULL)) {     
        // get the window width:

#ifdef _MSC_VER
        CONSOLE_SCREEN_BUFFER_INFO csbi;
        HANDLE hConsole = GetStdHandle( STD_OUTPUT_HANDLE );
        int ret = GetConsoleScreenBufferInfo(hConsole, &csbi);
        if(ret) {
          width = csbi.dwSize.X - 1;
        }
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
          
          RealType percent = value_ * 100.0 / maximum_;
          int hashes = int(percent * avail / 100.0);
          
          // compute the best estimate of the ending time:
          time_t current_ = time(NULL);
          time_t end_ = start_ + (current_ - start_) * (100.0/percent);
          struct tm * ender = localtime(&end_);
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
          if (hashes+1 < avail) {
            cout.width(hashes+1); 
            cout << progressSpinner_[iteration_ & 3];
          } else {
            cout.width(avail);
            cout << '#';            
          }
          cout.fill(' ');
          if (avail - hashes - 1  > 0) {
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
    value_ = val;
    maximum_ = max;
  }  
}
