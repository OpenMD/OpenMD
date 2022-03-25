/*
 * Copyright (c) 2004-2022, The University of Notre Dame. All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 * 1. Redistributions of source code must retain the above copyright notice,
 * this list of conditions and the following disclaimer.
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

/**
 * @file ifstrstream.hpp
 * @author Teng Lin
 * @date 10/14/2004
 * @version 1.0
 */

#ifndef IO_IFSTRSTREAM_HPP
#define IO_IFSTRSTREAM_HPP

#include <cassert>
#include <cstring>
#include <fstream>
#include <sstream>

namespace OpenMD {

  /**
   * @class ifstrstream ifstrstream.hpp "io/ifstrstream.hpp"
   * @brief ifstrstream class provides a stream interface to read data
   * from files.
   * <p>In single mode, it falls back to ifstream, as we don't need to
   * read the whole file into memory.  In parallel mode, the primary
   * node will read the whole file and broadcast it to other secondary
   * nodes.  After broadcasting, every node will fall back to
   * stringstream.</p>
   *
   * @code
   *       const int MAXLEN = 1024;
   *       char buffer[MAXLEN];
   *       ifstrstream in;
   *       in.open("Shapes.frc");
   *       if (in.is_open()) {
   *           in.getline(buffer, MAXLEN);
   *       }
   *       in.close();
   * @endcode
   */
  class ifstrstream : public std::basic_istream<char, std::char_traits<char>> {
  public:
    // traits
    using char_type   = char;
    using int_type    = std::char_traits<char>::int_type;
    using pos_type    = std::char_traits<char>::pos_type;
    using off_type    = std::char_traits<char>::off_type;
    using traits_type = std::char_traits<char>;

    using _Basic_ios = std::basic_ios<char, std::char_traits<char>>;
    using _Base      = std::basic_istream<char, std::char_traits<char>>;
    using _Buf       = std::basic_streambuf<char, std::char_traits<char>>;
    using _StringBuf = std::basic_stringbuf<char, std::char_traits<char>>;
    using _FileBuf   = std::basic_filebuf<char, std::char_traits<char>>;

    static const int FileNotExists = -1;
    static const int FileIOError   = -2;

  public:
    /**  Constructs an object of class ifstream.  */
    ifstrstream();

    /**
     * Explicit constructor
     * @param filename String containing the name of the file to be opened
     * @param mode Flags describing the requested i/o mode for the
     * file, default value is ios_base::in
     * @param checkFilename Flags indicating checking the file name in parallel
     */
    explicit ifstrstream(const char* filename,
                         std::ios_base::openmode mode = std::ios_base::in,
                         bool checkFilename           = false);

    /**
     * virtual destructor will close the file(in single mode) and
     * clear the stream buffer
     */
    ~ifstrstream();

    /**
     * Opens a file and associates a buffer with the specified file to
     * perform the i/o operations (single mode). The primary node reads
     * a file and broadcasts its content to the other secondary
     * nodes. After broadcasting, all nodes fall back to stringstream
     * (parallel mode).
     * @param filename String containing the name of the file to be opened
     * @param mode Flags describing the requested i/o mode for the file
     * @param checkFilename Flags indicating checking the file name in parallel
     */
    void open(const char* filename,
              std::ios_base::openmode mode = std::ios_base::in,
              bool checkFilename           = false);

    /**
     * Tests if the stream is currently associated with a valid buffer.
     * @return true if a file has successfully been opened (single
     * mode) or the whole file has been read and spread among all of
     * the processors (parallel mode), otherwise false is returned

     */
    bool is_open();

    /**
     * In single mode, closes a file. The stream's file buffer is
     * released from its association with the currently open file. In
     * parallel mode, clean up.
     */
    void close();

    /**
     * Gets the stream buffer object associated with the stream
     * @return A pointer to the stream buffer object (filebuf in
     * single mode, stringbuf in parallel mode) associated with the
     * stream.
     */
    _Buf* rdbuf();

  private:
    /**
     * Internal function used to open the file
     * @return true if succesfully opens a file (single mode) or gets the file
     * content (parallel mode) otherwise returns false
     * @param filename String containing the name of the file to be opened
     * @param mode Flags describing the requested i/o mode for the file
     * @param checkFilename Flags indicating checking the file name in parallel
     * @todo use try - catch syntax to make the program more readable
     */
    bool internalOpen(const char* filename, std::ios_base::openmode mode,
                      bool checkFilename);

    _StringBuf internalStringBuf_; /** internal stream buffer */
    _FileBuf internalFileBuf_;     /** internal stream buffer */
    bool isRead;                   /** file opened flag */
  };
}  // namespace OpenMD

#endif
