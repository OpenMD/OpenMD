/*
 * Copyright (c) 2005 The University of Notre Dame. All Rights Reserved.
 *
 * The University of Notre Dame grants you ("Licensee") a
 * non-exclusive, royalty free, license to use, modify and
 * redistribute this software in source and binary code form, provided
 * that the following conditions are met:
 *
 * 1. Acknowledgement of the program authors must be made in any
 *    publication of scientific results based in part on use of the
 *    program.  An acceptable form of acknowledgement is citation of
 *    the article in which the program was described (Matthew
 *    A. Meineke, Charles F. Vardeman II, Teng Lin, Christopher
 *    J. Fennell and J. Daniel Gezelter, "OOPSE: An Object-Oriented
 *    Parallel Simulation Engine for Molecular Dynamics,"
 *    J. Comput. Chem. 26, pp. 252-271 (2005))
 *
 * 2. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 *
 * 3. Redistributions in binary form must reproduce the above copyright
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
 */
 
/**
 * @file basic_ifstrstream.hpp
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

#ifdef IS_MPI
#include <mpi.h>
#endif

namespace oopse {

  /**
   * @class basic_ifstrstream basic_ifstrstream.hpp "io/basic_ifstrstream.hpp"
   * @brief basic_ifstrstream class provides a stream interface to read data from files.
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
  template <class _CharT, class _Traits,  class _Alloc>
  class basic_ifstrstream : public std::basic_istream<_CharT, _Traits> {
  public:
    //traits
    typedef _CharT                     char_type;
    typedef typename _Traits::int_type int_type;
    typedef typename _Traits::pos_type pos_type;
    typedef typename _Traits::off_type off_type;
    typedef _Traits                    traits_type;

    typedef std::basic_ios<_CharT, _Traits>                _Basic_ios;
    typedef std::basic_istream<_CharT, _Traits>            _Base;

#ifdef IS_MPI
    typedef std::basic_stringbuf<_CharT, _Traits, _Alloc>  _Buf;        
#else
    typedef std::basic_filebuf<_CharT, _Traits>            _Buf;
#endif

    static  const int FileNotExists = -1;
    static const int FileIOError = -2;
        
  public:
        
    /**  Constructs an object of class ifstream.  */
    basic_ifstrstream()
      :  std::basic_istream<_CharT, _Traits>(0),
	 internalBuf_(), isRead(false)  {
    
	this->init(&internalBuf_);
      }
        
    /**
     * Explicit constructor
     * @filename String containing the name of the file to be opened
     * @mode Flags describing the requested i/o mode for the file, default value is ios_base::in      
     * @checkFilename Flags indicating checking the file name in parallel
     */
    explicit basic_ifstrstream(const char* filename, std::ios_base::openmode mode = std::ios_base::in, bool checkFilename = false)
      :  std::basic_istream<_CharT, _Traits>(0),
	 internalBuf_(), isRead(false) {

	this->init(&internalBuf_);
	isRead =  internalOpen(filename,  mode | std::ios_base::in, checkFilename);
      }

    /**
     * virtual destructor will close the file(in single mode) and clear the stream buffer
     */
    ~basic_ifstrstream(){
      close();
    }

    /**
     * Opens a file and associats a buffer with the specified file to perform the i/o operations 
     * (single mode). Master reads a file and brocasts its content to the other slave nodes. After
     * brocasting, every nodes fall back to stringstream (parallel mode).
     * @filename String containing the name of the file to be opened
     * @mode Flags describing the requested i/o mode for the file
     * @checkFilename Flags indicating checking the file name in parallel
     */
    void open(const char* filename, std::ios_base::openmode mode = std::ios_base::in, bool checkFilename = false){

      if (!isRead ) {
	isRead = internalOpen(filename, mode, checkFilename);
      }
    }

    /**
     * Tests if the stream is currently associated with a valid  buffer.
     * @return true if a file has successfully been opened (single mode) or the whole file is read 
     * and spreaded among all of the processors (parallel mode),  otherwise false is returned
     */
    bool is_open ( ) {
#ifdef IS_MPI            
      return isRead; 
#else
      //single version fall back to ifstream
      return internalBuf_.is_open();
#endif
    }

    /**
     * In single mode, closes a file. The stream's file buffer is released from its association with
     * the currently open file. In parallel mode, clean the 
     */
    void close() {
#ifndef IS_MPI            
      //single version fall back to ifstream
      if (!internalBuf_.close())
	this->setstate(std::ios_base::failbit);
#endif             

      isRead = false;
    }
        
    /**
     * Gets the stream buffer object associated with the stream
     * @return   A pointer to the  stream buffe object(filebuf in single mode, stringbuf in 
     * parallel mode) associated with the stream.
     */
    _Buf* rdbuf() const{
      return static_cast<_Buf*>(&internalBuf_); 
    }

  private:

    /**
     * Internal function used to open the file
     * @return true if succesfully opens a file (single mode) or gets the file content (parallel mode)
     * otherwise return false
     * @filename String containing the name of the file to be opened
     * @mode Flags describing the requested i/o mode for the file
     * @todo use try - catch syntax to make the program more readable
     */
    bool internalOpen(const char* filename, std::ios_base::openmode mode, bool checkFilename){

#ifdef IS_MPI         
      int commStatus;
      long fileSize;
      char* fbuf;
      int filenameLen;
      int diffFilename;
      int error;
      int myRank;
      int masterNode;

      commStatus = MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
      masterNode = 0;
            
      if (myRank == masterNode) {

	if (checkFilename) {

	  //check the filename is the same
	  filenameLen = strlen(filename);
	  commStatus = MPI_Bcast(&filenameLen, 1, MPI_INT, masterNode, MPI_COMM_WORLD);     
	  commStatus = MPI_Bcast((void*)filename, filenameLen, MPI_CHAR, masterNode, MPI_COMM_WORLD);     

	  diffFilename = 0;
	  commStatus = MPI_Allreduce(&diffFilename, &error, 1,  MPI_INT, MPI_SUM,  MPI_COMM_WORLD);             

	  //if file names are different just return false
	  if (error > 0)
	    return false;   
	}
                
	std::ifstream fin(filename, mode);

	if (fin.is_open()) {
                    
	  fin.seekg(0, std::ios::end); 
	  fileSize = fin.tellg(); 
	  fin.seekg(0, std::ios::beg);
                    
	  // '\0' need one more char
	  fbuf = new char[fileSize+1];
                    
	  assert(fbuf != 0);

	  fin.read(fbuf, fileSize);

	  if (fin.fail())
	    fileSize = FileIOError;
                    
	  //brocasting the file size
	  commStatus = MPI_Bcast(&fileSize, 1, MPI_LONG, masterNode, MPI_COMM_WORLD);                   

	  if (fileSize < 0) {
	    fin.close();                    
	    delete fbuf;
                        
	    return false;
	  }
                    
	  // make a c-style  std::string and brocasting it
	  fbuf[fileSize] = '\0';
	  commStatus = MPI_Bcast(fbuf, fileSize + 1, MPI_CHAR, masterNode, MPI_COMM_WORLD); 

	  //close the file and delete the buffer
	  fin.close();      
	  internalBuf_.str(fbuf);
	  delete fbuf;
	}else{
	  fileSize = FileNotExists;
	  commStatus = MPI_Bcast(&fileSize, 1, MPI_LONG, masterNode, MPI_COMM_WORLD);   
	  return false;
	}
               
      } else{ //slave nodes

	//check file name
	if (checkFilename) {
	  commStatus = MPI_Bcast(&filenameLen, 1, MPI_INT, masterNode, MPI_COMM_WORLD);     

	  char * masterFilename = new char[filenameLen];
	  commStatus = MPI_Bcast(masterFilename, filenameLen, MPI_CHAR, masterNode, MPI_COMM_WORLD);     
        
	  if( strcmp(masterFilename, filename) == 0)
	    diffFilename = 0;
	  else
	    diffFilename = 1;

	  delete masterFilename;
                        
	  commStatus = MPI_Allreduce(&diffFilename, &error, 1,  MPI_INT,  MPI_SUM, MPI_COMM_WORLD);    

	  if (error > 0)
	    return false;                        
	}
	//get file size
	commStatus = MPI_Bcast(&fileSize, 1, MPI_LONG, masterNode, MPI_COMM_WORLD);   

	if (fileSize >= 0 ) {
	  fbuf = new char[fileSize+1];
	  assert(fbuf);

	  //receive file content
	  commStatus = MPI_Bcast(fbuf, fileSize + 1, MPI_CHAR, masterNode, MPI_COMM_WORLD); 

	  internalBuf_.str(fbuf);
	  delete fbuf;

	} else if (fileSize == FileNotExists ) {
	  return false;

	} else if (fileSize == FileIOError ) {
	  return false;
	}
      }

#else
      //in single version, fall back to ifstream
      if (!internalBuf_.open(filename, mode)) {
	this->setstate(std::ios_base::failbit);
	return false;
      }    

#endif
      this->clear();
      return true;
    }
        
    _Buf  internalBuf_; /** internal stream buffer */        
    bool isRead;                                                                    /** file opened flag */
  };

  typedef basic_ifstrstream<char, std::char_traits<char>, std::allocator<char> > ifstrstream;
}//namespace oopse
#endif //IO_IFSTRSTREAM_HPP
