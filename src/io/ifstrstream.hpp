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
 * [3]  Sun, Lin & Gezelter, J. Chem. Phys. 128, 24107 (2008).          
 * [4]  Vardeman & Gezelter, in progress (2009).                        
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
   * @brief ifstrstream class provides a stream interface to read data from files.
   * <p>In single mode, it falls back to ifstream. Don't need to read the whole file into memory.
   * In parallel mode, the master node will read the whole file and brocast it to other slave nodes.
   * After brocasting, every node will fall back to stringstream.</p>
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
class ifstrstream : public std::basic_istream<char, std::char_traits<char> > {
  public:
    //traits
    typedef char                     char_type;
    typedef std::char_traits<char>::int_type int_type;
    typedef std::char_traits<char>::pos_type pos_type;
    typedef std::char_traits<char>::off_type off_type;
    typedef std::char_traits<char>                    traits_type;

    typedef std::basic_ios<char, std::char_traits<char> >       _Basic_ios;
    typedef std::basic_istream<char, std::char_traits<char> >   _Base;   
    typedef std::basic_streambuf<char, std::char_traits<char> > _Buf;
    typedef std::basic_stringbuf<char, std::char_traits<char> > _StringBuf;
    typedef std::basic_filebuf<char, std::char_traits<char> >   _FileBuf;
  

    static const int FileNotExists = -1;
    static const int FileIOError = -2;
        
  public:
        
    /**  Constructs an object of class ifstream.  */
    ifstrstream();
        
    /**
     * Explicit constructor
     * @filename String containing the name of the file to be opened
     * @mode Flags describing the requested i/o mode for the file, default value is ios_base::in      
     * @checkFilename Flags indicating checking the file name in parallel
     */
    explicit ifstrstream(const char* filename, std::ios_base::openmode mode = std::ios_base::in, bool checkFilename = false);

    /**
     * virtual destructor will close the file(in single mode) and clear the stream buffer
     */
    ~ifstrstream();

    /**
     * Opens a file and associats a buffer with the specified file to perform the i/o operations 
     * (single mode). Master reads a file and brocasts its content to the other slave nodes. After
     * brocasting, every nodes fall back to stringstream (parallel mode).
     * @filename String containing the name of the file to be opened
     * @mode Flags describing the requested i/o mode for the file
     * @checkFilename Flags indicating checking the file name in parallel
     */
    void open(const char* filename, std::ios_base::openmode mode = std::ios_base::in, bool checkFilename = false);


    /**
     * Tests if the stream is currently associated with a valid  buffer.
     * @return true if a file has successfully been opened (single mode) or the whole file is read 
     * and spreaded among all of the processors (parallel mode),  otherwise false is returned
     */
    bool is_open ( );

    /**
     * In single mode, closes a file. The stream's file buffer is released from its association with
     * the currently open file. In parallel mode, clean the 
     */
    void close();
    
    /**
     * Gets the stream buffer object associated with the stream
     * @return   A pointer to the stream buffer object(filebuf in single mode, 
     * stringbuf in parallel mode) associated with the stream.
     */
  _Buf* rdbuf();

  private:

    /**
     * Internal function used to open the file
     * @return true if succesfully opens a file (single mode) or gets the file content (parallel mode)
     * otherwise return false
     * @filename String containing the name of the file to be opened
     * @mode Flags describing the requested i/o mode for the file
     * @todo use try - catch syntax to make the program more readable
     */
    bool internalOpen(const char* filename, std::ios_base::openmode mode, bool checkFilename);
  
  _StringBuf   internalStringBuf_; /** internal stream buffer */        
  _FileBuf     internalFileBuf_;    /** internal stream buffer */        
  bool isRead;        /** file opened flag */
};
  
}
#endif
