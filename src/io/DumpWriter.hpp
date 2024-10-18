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

#ifndef IO_DUMPWRITER_HPP
#define IO_DUMPWRITER_HPP

#define _LARGEFILE_SOURCE64
#ifndef _FILE_OFFSET_BITS
#define _FILE_OFFSET_BITS 64
#endif
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iostream>
#include <string>

#include <sys/stat.h>
#include <sys/types.h>

#include "brains/SimInfo.hpp"
#include "brains/Thermo.hpp"
#include "primitives/Atom.hpp"
#include "primitives/StuntDouble.hpp"

namespace OpenMD {

  /**
   * @class DumpWriter DumpWriter.hpp "io/DumpWriter.hpp"
   * @todo
   */
  class DumpWriter {
  public:
    DumpWriter(SimInfo* info);
    DumpWriter(SimInfo* info, const std::string& filename);
    DumpWriter(SimInfo* info, const std::string& filename, bool writeDumpFile);
    ~DumpWriter();

    void writeDumpAndEor();
    void writeDump();
    void writeEor();

  private:
    void writeFrame(std::ostream& os);
    void writeFrameProperties(std::ostream& os, Snapshot* s);
    std::string prepareDumpLine(StuntDouble* sd);
    std::string prepareSiteLine(StuntDouble* sd, int ioIndex, int siteIndex);
    std::ostream* createOStream(const std::string& filename);
    void writeClosing(std::ostream& os);

    SimInfo* info_ {nullptr};
    std::string filename_;
    std::ostream* dumpFile_ {};
    std::string eorFilename_;
    bool needCompression_;
    bool needForceVector_;
    bool needParticlePot_;
    bool needFlucQ_;
    bool needElectricField_;
    bool needSitePotential_;
    bool needDensity_;
    bool doSiteData_;
    bool createDumpFile_;
  };
}  // namespace OpenMD

#endif  // #define IO_DUMPWRITER_HPP
