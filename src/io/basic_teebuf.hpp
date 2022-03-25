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

#ifndef UTILS_MULTISTREAMBUF_HPP
#define UTILS_MULTISTREAMBUF_HPP

#include <streambuf>
#include <vector>

namespace OpenMD {

  /**
   * @class basic_teebuf basic_teebuf.hpp "utils/basic_teebuf.hpp"
   * @brief As a subclass of basic_streambuf,  basic_teebuf can operate on
   * multiple stream simultaneously.
   * @code
   *    std::ofstream file1("file1");
   *    std::ofstream file2("file22");
   *    std::vector<std::streambuf*> buffers;
   *    buffers.push_back(file1.rdbuf());
   *    buffers.push_back(file2.rdbuf());
   *    teebuf tmp(buffers.begin(), buffers.end());
   *    std::ostream myOs(&tmp);
   *    myOs << "hello world";
   * @endcode
   */
  template<class CharT, class Traits = std::char_traits<CharT>>
  class basic_teebuf : public std::basic_streambuf<CharT, Traits> {
  public:
    using streambuf_type = std::basic_streambuf<CharT, Traits>;
    using traits_type    = Traits;
    using char_type      = CharT;
    using int_type       = typename traits_type::int_type;

    template<typename ForwardIterator>
    basic_teebuf(ForwardIterator begin, ForwardIterator end) :
        buffers_(begin, end) {}

  protected:
    int_type overflow(int_type c = traits_type::eof()) {
      // avoid writing eof to stream
      if (c == traits_type::eof()) { return traits_type::eof(); }

      typename std::vector<streambuf_type*>::iterator
          iter;  // typename is needed since it's a dependant name
      for (iter = buffers_.begin(); iter != buffers_.end(); ++iter) {
        if ((*iter)->sputc(c) == traits_type::eof()) {
          return traits_type::eof();
        }
      }

      return traits_type::not_eof(c);
    }

    int sync() {
      // flush buffer, checking return for eof.
      if (traits_type::eq_int_type(overflow(traits_type::eof()),
                                   traits_type::eof())) {
        return -1;
      }

      // flush streams
      typename std::vector<streambuf_type*>::iterator iter;
      for (iter = buffers_.begin(); iter != buffers_.end(); ++iter) {
        if ((*iter)->pubsync() == -1) { return -1; }
      }

      return 0;
    }

  private:
    std::vector<streambuf_type*> buffers_;
  };

  using TeeBuf = basic_teebuf<char>;
}  // namespace OpenMD

#endif
