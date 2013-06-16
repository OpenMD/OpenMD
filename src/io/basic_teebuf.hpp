/*
 * Copyright (c) 2005 The University of Notre Dame. All Rights Reserved.
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
 * [3]  Sun, Lin & Gezelter, J. Chem. Phys. 128, 234107 (2008).          
 * [4]  Kuang & Gezelter,  J. Chem. Phys. 133, 164101 (2010).
 * [5]  Vardeman, Stocker & Gezelter, J. Chem. Theory Comput. 7, 834 (2011).
 */

#ifndef UTILS_MULTISTREAMBUF_HPP
#define UTILS_MULTISTREAMBUF_HPP
#include <streambuf>
#include <vector>
namespace OpenMD {

  /**
   * @class basic_teebuf basic_teebuf.hpp "utils/basic_teebuf.hpp"
   * @brief As a subclass of basic_streambuf,  basic_teebuf can operate on multiple stream simultaneously.
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

  template <class CharT, class Traits = std::char_traits<CharT> > 
  class basic_teebuf: public std::basic_streambuf<CharT, Traits> { 
  public: 
    typedef std::basic_streambuf<CharT, Traits> streambuf_type; 
    typedef Traits traits_type; 
    typedef CharT char_type; 
    typedef typename traits_type::int_type int_type; 

    template <typename ForwardIterator>
    basic_teebuf(ForwardIterator begin, ForwardIterator end) : buffers_(begin, end){

    }

  protected:
    int_type overflow(int_type c = traits_type::eof()) {

      //avoid writing eof to stream
      if (c == traits_type::eof()) {
	return traits_type::eof();
      }
            
      typename std::vector<streambuf_type*>::iterator iter; //typename is needed since it's a dependant name
      for (iter = buffers_.begin(); iter != buffers_.end(); ++iter) {
	if ((*iter)->sputc(c) ==  traits_type::eof()) {
	  return  traits_type::eof();
	}
      }

      return traits_type::not_eof(c); 
    }

    int sync() { 

      //flush buffer, checking return for eof. 
      if (traits_type::eq_int_type(overflow(traits_type::eof()), traits_type::eof())) { 
	return -1; 
      } 

      //flush streams
      typename std::vector<streambuf_type*>::iterator iter;
      for (iter = buffers_.begin(); iter != buffers_.end(); ++iter) {
	if ((*iter)->pubsync() == -1) {
	  return -1;
	}
      }
            
      return 0; 
    } 
   
  private:
    std::vector<streambuf_type*> buffers_;
  };

  typedef basic_teebuf<char> TeeBuf; 

}
#endif
