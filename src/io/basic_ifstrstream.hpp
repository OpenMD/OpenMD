/*
 * Copyright (C) 2000-2004  Object Oriented Parallel Simulation Engine (OOPSE) project
 * 
 * Contact: oopse@oopse.org
 * 
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public License
 * as published by the Free Software Foundation; either version 2.1
 * of the License, or (at your option) any later version.
 * All we ask is that proper credit is given for our work, which includes
 * - but is not limited to - adding the above copyright notice to the beginning
 * of your source code files, and to any copyright notice that you may distribute
 * with programs based on this work.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 * 
 * You should have received a copy of the GNU Lesser General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
 *
 */

/**
 * @file basic_ifstrstream.hpp
 * @author Teng Lin
 * @date 10/14/2004
 * @version 1.0
 */

#ifndef IO_IFSTRSTREAM_HPP
#define IO_IFSTRSTREAM_HPP

#include <fstream>
#include <sstream>

namespace oopse {
using namespace std;
/**
 * @class basic_ifstrstream basic_ifstrstream.hpp "io/basic_ifstrstream.hpp"
 * @brief class provides a stream interface to read data from files.
 * <p>In single mode, it falls back to ifstream. Don't need to read the whole file into memory.
 * In parallel mode, the master node will read the whole file and brocast it to other slave nodes.
 * After brocasting, every node will fall back to stringstream.</p>
 *
 * @code
 *       const int MAXLEN = 1024;
 *       char buffer[MAXLEN];
 *       ifstrstream in;
 *       in.open("Shapes.frc");
 *       if (in.is_open) {
 *           in.getline(buffer, MAXLEN);
 *       }
 *       in.close();
 * @endcode
 */
template <class _CharT, class _Traits,  class _Alloc>
class basic_ifstrstream : public basic_istream<_CharT, _Traits> {
    public:
        //traits
        typedef _CharT                     char_type;
        typedef typename _Traits::int_type int_type;
        typedef typename _Traits::pos_type pos_type;
        typedef typename _Traits::off_type off_type;
        typedef _Traits                    traits_type;

        typedef basic_ios<_CharT, _Traits>                _Basic_ios;
        typedef basic_istream<_CharT, _Traits>            _Base;

#ifdef IS_MPI
         typedef basic_stringbuf<_CharT, _Traits, _Alloc>  _Buf;        
#else
        typedef basic_filebuf<_CharT, _Traits>            _Buf;
#endif

        static  const int FileNoExists = -1;
        static const int FileIOError = -2;
        
    public:
        
        /**  Constructs an object of class ifstream.  */
        basic_ifstrstream()
            : basic_ios<_CharT, _Traits>(),  basic_istream<_CharT, _Traits>(0),
              internalBuf_(NULL), isRead(false)  {

#ifdef IS_MPI         
            //in parallel mode, fall back to istringstream
            basic_stringbuf<_CharT, _Traits, _Alloc>* stringBuffer = new  basic_stringbuf<_CharT, _Traits, _Alloc>();
            internalBuf_ =  stringBuffer;
#else
            //in single version, fall back to ifstream
            basic_filebuf<_CharT, _Traits>* fileBuffer = new  basic_filebuf<_CharT, _Traits>();
            internalBuf_ =  fileBuffer;
#endif            

            this->init(internalBuf_);
            isRead = false;
        }
        
        /**
         * Explicit constructor
         * @filename String containing the name of the file to be opened
         * @mode Flags describing the requested i/o mode for the file, default value is ios_base::in         
         */
        explicit basic_ifstrstream(const char* filename, ios_base::openmode mode = ios_base::in)
            : basic_ios<_CharT, _Traits>(),  basic_istream<_CharT, _Traits>(0),
              internalBuf_(NULL), isRead(false) {

           isRead =  internalOpen(filename,  mode | ios_base::in);
         }

        /**
         *
         */
        ~basic_ifstrstream(){
            close();
            delete internalBuf_;
            internalBuf_ = NULL;
        }

        /**
         * Opens a file and associats a buffer with the specified file to perform the i/o operations 
         * (single mode). Master reads a file and brocasts its content to the other slave nodes. After
         * brocasting, every nodes fall back to stringstream (parallel mode).
         * @filename String containing the name of the file to be opened
         * @mode Flags describing the requested i/o mode for the file
         */
        void open(const char* filename, ios_base::openmode mode = ios_base::in){

            if (!isRead ) {
                isRead = internalOpen();
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
            return this->rdbuf()->is_open();
#endif
        }

        /**
         * In single mode, closes a file. The stream's file buffer is released from its association with
         * the currently open file. In parallel mode, clean the 
         */
        void close() {
#ifndef IS_MPI            
            //single version fall back to ifstream
            if (!this->rdbuf()->close())
                this->setstate(ios_base::failbit);
#endif             

            isRead = false;
        }
        
        /**
         * Gets the stream buffer object associated with the stream
         * @return   A pointer to the  stream buffe object(filebuf in single mode, stringbuf in 
         * parallel mode) associated with the stream.
         */
        _Buf* rdbuf() const{
            return const_cast<_Buf*>(internalBuf_); 
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
        bool internalOpen(const char* filename, ios_base::openmode mode){

#ifdef IS_MPI         
            int commStatus;
            long fileSize;
            char* fbuf;
            int filenameLen;
            int diffFilename;
            int error;

            if (worldRank == masterNode) {

                //check the filename is the same
                filenameLen = strlen(filename);
                commStatus = MPI_Bcast(filenameLen, 1, MPI_INT, masterNode, MPI_COMM_WORLD);     
                commStatus = MPI_Bcast(filename, filenameLen, MPI_CHAR, masterNode, MPI_COMM_WORLD);     

                diffFilename = 0;
                commStatus = MPI_Allreduce(sameFilename, error, 1,  MPI_INT,  MPI_COMM_WORLD);             

                //if file names are different just return false
                if (error > 0)
                        return false;   
                
                ifstream fin(filename, mod);
                basic_stringbuf<_CharT, _Traits, _Alloc>* sbuf;

                if (fin.is_open()) {
                    
                    fin.in.seekg(0, ios::end); 
                    fileSize = fin.tellg(); 

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
                    
                    // make a c-style string and brocasting it
                    fbuf[fileSize] = '\0';
                    commStatus = MPI_Bcast(fbuf, fileSize + 1, MPI_CHAR, masterNode, MPI_COMM_WORLD); 

                    //it is safe to delete null pointer
                    delete internalBuf_;

                    //initilaize istream
                    internalBuf_  = new basic_stringbuf<_CharT, _Traits, _Alloc>(fbuf, mode);
                    assert(internalBuf_);
                    this->init(internalBuf_);

                    //close the file and delete the buffer
                    fin.close();                    
                    delete fbuf;
                }else{
                    fileSize = FileNoExists;
                    commStatus = MPI_Bcast(&fileSize, 1, MPI_LONG, masterNode, MPI_COMM_WORLD);   
                    return false;
                }
               
             } else{
                    //check file name
                    commStatus = MPI_Bcast(filenameLen, 1, MPI_INT, masterNode, MPI_COMM_WORLD);     

                    char * masterFilename = new char[filenameLen];
                    commStatus = MPI_Bcast(masterFilename, filenameLen, MPI_CHAR, masterNode, MPI_COMM_WORLD);     
    
                    if( strcmp(masterFilename, filename) == 0)
                        diffFilename = 0;
                    else
                        diffFilename = 1;

                    delete masterFilename;
                    
                    commStatus = MPI_Allreduce(sameFilename, error, 1,  MPI_INT,  MPI_COMM_WORLD);    

                    if (error > 0)
                        return false;                        

                    //get file size
                    commStatus = MPI_Bcast(&fileSize, 1, MPI_LONG, masterNode, MPI_COMM_WORLD);   

                    if (fileSize >= 0 ) {
                        fbuf = new char[fileSize+1];
                        assert(fbuf);

                        //receive file content
                        commStatus = MPI_Bcast(fbuf, fileSize + 1, MPI_CHAR, masterNode, MPI_COMM_WORLD); 

                        //it is safe to delete null pointer
                        delete internalBuf_;

                        //initilaize istream                        
                        internalBuf_  = new basic_stringbuf<_CharT, _Traits, _Alloc>(fbuf, mode);
                        assert(internalBuf_);
                        this->init(internalBuf_);

                        delete fbuf;

                    } else if (fileSize == FileNoExists ) {
                        return false;

                    } else if (fileSize == FileIOError ) {
                        return false;
                    }
             }

#else
            //in single version, fall back to ifstream
            basic_filebuf<_CharT, _Traits>* fileBuffer = new  basic_filebuf<_CharT, _Traits>();

            this->init(fileBuffer);
            if (!fileBuffer->open(filename, mode))
                this->setstate(ios_base::failbit);

            //it is safe to delete null pointer
            delete internalBuf_;
 
           internalBuf_ =  fileBuffer;
#endif

            return true;
        }
        
        basic_streambuf<_CharT, _Traits>*  internalBuf_; /** internal stream buffer */        
        bool isRead;                                                                    /** file opened flag */
};

typedef basic_istringstream<char, char_traits<char>, allocator<char> > ifstringstream;
}//namespace oopse
#endif //IO_IFSTRSTREAM_HPP
