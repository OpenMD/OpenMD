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
 
/**
 * @file ifstrstream.cpp
 * @author Teng Lin
 * @date 10/14/2004
 * @version 1.0
 */

#include "io/ifstrstream.hpp"

#ifdef IS_MPI
#include <mpi.h>
#endif

namespace OpenMD {
  
  /**  Constructs an object of class ifstream.  */
#ifdef IS_MPI
  ifstrstream::ifstrstream() 
    :  std::basic_istream<char, std::char_traits<char> >(0),
       internalStringBuf_(), isRead(false)  {    
    this->init(&internalStringBuf_);
  }
#else
  ifstrstream::ifstrstream() 
    :  std::basic_istream<char, std::char_traits<char> >(0),
       internalFileBuf_(), isRead(false)  {    
    this->init(&internalFileBuf_);
  }
#endif
  
  /**
   * Explicit constructor
   * @param filename String containing the name of the file to be opened
   * @param mode Flags describing the requested i/o mode for the file, default value is ios_base::in      
   * @param checkFilename Flags indicating checking the file name in parallel
   */
#ifdef IS_MPI
  ifstrstream::ifstrstream(const char* filename, std::ios_base::openmode mode, bool checkFilename)
    :  std::basic_istream<char, std::char_traits<char> >(0),
       internalStringBuf_(), isRead(false) {
    this->init(&internalStringBuf_);
    isRead =  internalOpen(filename,  mode | std::ios_base::in, checkFilename);
  }
#else
  ifstrstream::ifstrstream(const char* filename, std::ios_base::openmode mode, bool checkFilename)
    :  std::basic_istream<char, std::char_traits<char> >(0),
       internalFileBuf_(), isRead(false) {
    
    this->init(&internalFileBuf_);
    isRead =  internalOpen(filename,  mode | std::ios_base::in, checkFilename);
  }
#endif  
  /**
   * virtual destructor will close the file(in single mode) and clear the stream buffer
   */
  ifstrstream::~ifstrstream(){
    close();
  }
  
  /**
   * Opens a file and associats a buffer with the specified file to perform the i/o operations 
   * (single mode). Master reads a file and brocasts its content to the other slave nodes. After
   * brocasting, every nodes fall back to stringstream (parallel mode).
   * @param filename String containing the name of the file to be opened
   * @param mode Flags describing the requested i/o mode for the file
   * @param checkFilename Flags indicating checking the file name in parallel
   */
  void ifstrstream::open(const char* filename, std::ios_base::openmode mode, bool checkFilename){
    
    if (!isRead ) {
      isRead = internalOpen(filename, mode, checkFilename);
    }
  }

  /**
   * Tests if the stream is currently associated with a valid  buffer.
   * @return true if a file has successfully been opened (single mode) or the whole file is read 
   * and spreaded among all of the processors (parallel mode),  otherwise false is returned
   */
  bool ifstrstream::is_open ( ) {
#ifdef IS_MPI            
    return isRead; 
#else
    //single version fall back to ifstream
    return internalFileBuf_.is_open();
#endif
  }
  
  /**
   * In single mode, closes a file. The stream's file buffer is
   * released from its association with the currently open file. In
   * parallel mode, clean the
   */
  void ifstrstream::close() {
#ifndef IS_MPI            
    //single version fall back to ifstream
    if (!internalFileBuf_.close())
      this->setstate(std::ios_base::failbit);
#endif             
    
    isRead = false;
  }
  
  /**
   * Gets the stream buffer object associated with the stream
   * @return A pointer to the stream buffer object(filebuf in single
   * mode, stringbuf in parallel mode) associated with the stream.
   */
  std::basic_streambuf<char, std::char_traits<char> >*  ifstrstream::rdbuf() {
#ifdef IS_MPI
    return static_cast<_StringBuf*>(&internalStringBuf_); 
#else
    return static_cast<_FileBuf*>(&internalFileBuf_); 
#endif
  }  
  
  /**
     * Internal function used to open the file
     * @return true if succesfully opens a file (single mode) or gets the file content (parallel mode)
     * otherwise return false
     * @param filename String containing the name of the file to be opened
     * @param mode Flags describing the requested i/o mode for the file
     * @param checkFilename Flags indicating checking the file name in parallel
     * @todo use try - catch syntax to make the program more readable
     */

  bool ifstrstream::internalOpen(const char* filename, std::ios_base::openmode mode, bool checkFilename){
    
#ifdef IS_MPI         
    //int commStatus;
    long fileSize;
    char* fbuf;
    int filenameLen;
    int diffFilename;
    int error;
    int myRank;
    int masterNode;
    
    myRank =  MPI::COMM_WORLD.Get_rank();
    masterNode = 0;
    
    if (myRank == masterNode) {
      
      if (checkFilename) {
        
        //check the filename is the same
        filenameLen = strlen(filename);
        MPI::COMM_WORLD.Bcast(&filenameLen, 1, MPI::INT, masterNode);
        MPI::COMM_WORLD.Bcast((void*)filename, filenameLen, MPI::CHAR, 
                              masterNode);
        
        diffFilename = 0;
        MPI::COMM_WORLD.Allreduce(&diffFilename, &error, 1, MPI::INT, MPI::SUM);

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
        
        //broadcast the file size
        MPI::COMM_WORLD.Bcast(&fileSize, 1, MPI::LONG, masterNode);
        
        if (fileSize < 0) {
          fin.close();                    
          delete fbuf;
          
          return false;
        }
        
        // make a c-style  std::string and broadcast it
        fbuf[fileSize] = '\0';
        MPI::COMM_WORLD.Bcast(fbuf, fileSize + 1, MPI::CHAR, masterNode);
        
        //close the file and delete the buffer
        fin.close();      
        internalStringBuf_.str(fbuf);
        delete [] fbuf;
      }else{
        fileSize = FileNotExists;
        MPI::COMM_WORLD.Bcast(&fileSize, 1, MPI::LONG, masterNode);
        return false;
      }
      
    } else{ //slave nodes
      
      //check file name
      if (checkFilename) {
        MPI::COMM_WORLD.Bcast(&filenameLen, 1, MPI::INT, masterNode);
        
        char * masterFilename = new char[filenameLen];
        MPI::COMM_WORLD.Bcast(masterFilename, filenameLen, MPI::CHAR, 
                              masterNode);
        
        if( strcmp(masterFilename, filename) == 0)
          diffFilename = 0;
        else
          diffFilename = 1;
        
        delete masterFilename;
        
        MPI::COMM_WORLD.Allreduce(&diffFilename, &error, 1, MPI::INT, MPI::SUM);
        
        if (error > 0)
          return false;                        
      }
      //get file size
      MPI::COMM_WORLD.Bcast(&fileSize, 1, MPI::LONG, masterNode);
      
      if (fileSize >= 0 ) {
        fbuf = new char[fileSize+1];
        assert(fbuf);
        
        //receive file content
        MPI::COMM_WORLD.Bcast(fbuf, fileSize + 1, MPI::CHAR, masterNode);
        
        internalStringBuf_.str(fbuf);
        delete [] fbuf;
        
      } else if (fileSize == FileNotExists ) {
        return false;
        
      } else if (fileSize == FileIOError ) {
        return false;
      }
    }
    
#else
    //in single version, fall back to ifstream
    if (!internalFileBuf_.open(filename, mode)) {
      this->setstate(std::ios_base::failbit);
      return false;
    }    
    
#endif
    this->clear();
    return true;
  }
}
