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
 * research, please cite the following paper when you publish your work:
 *
 * [1] Drisko et al., J. Open Source Softw. 9, 7004 (2024).
 *
 * Good starting points for code and simulation methodology are:
 *
 * [2] Meineke, et al., J. Comp. Chem. 26, 252-271 (2005).
 * [3] Fennell & Gezelter, J. Chem. Phys. 124, 234104 (2006).
 * [4] Sun, Lin & Gezelter, J. Chem. Phys. 128, 234107 (2008).
 * [5] Vardeman, Stocker & Gezelter, J. Chem. Theory Comput. 7, 834 (2011).
 * [6] Kuang & Gezelter, Mol. Phys., 110, 691-701 (2012).
 * [7] Lamichhane, Gezelter & Newman, J. Chem. Phys. 141, 134109 (2014).
 * [8] Bhattarai, Newman & Gezelter, Phys. Rev. B 99, 094106 (2019).
 * [9] Drisko & Gezelter, J. Chem. Theory Comput. 20, 4986-4997 (2024).
 */

/**
 * @file DumpReader.hpp
 * @author tlin
 * @date 11/15/2004
 * @version 2.0
 */

#ifndef IO_DUMPREADER_HPP
#define IO_DUMPREADER_HPP

#include <cstdio>
#include <string>

#include "brains/SimInfo.hpp"
#include "primitives/StuntDouble.hpp"

namespace OpenMD {

  /**
   * @class DumpReader DumpReader.hpp "io/DumpReader.hpp"
   * @todo get rid of more junk code from DumpReader
   */
  class DumpReader final {
  public:
    DumpReader(SimInfo* info, const std::string& filename);
    ~DumpReader();

    /** Returns the number of frames in the dump file*/
    int getNFrames();

    void setNeedCOMprops(bool ncp) { needCOMprops_ = ncp; }

    void readFrame(int whichFrame);

  protected:
    void scanFile();
    void readSet(int whichFrame);
    void parseDumpLine(const std::string&);
    void parseSiteLine(const std::string&);
    void readFrameProperties(std::istream& inputStream);
    int readStuntDoubles(std::istream& inputStream);
    void readSiteData(std::istream& inputStream);

    SimInfo* info_ {nullptr};

    std::string filename_;
    bool isScanned_;

    int nframes_;

    std::ifstream inFile_;

    std::vector<std::streampos> framePos_;

    bool needPos_;
    bool needVel_;
    bool needQuaternion_;
    bool needAngMom_;
    bool needCOMprops_;
    bool readField_;

    const static int bufferSize = 4096;
    char buffer[bufferSize];
  };
}  // namespace OpenMD

#endif  // IO_DUMPREADER_HPP
