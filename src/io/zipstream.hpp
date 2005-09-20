/*
zipstream Library License:
--------------------------

The zlib/libpng License Copyright (c) 2003 Jonathan de Halleux.

This software is provided 'as-is', without any express or implied warranty. In no event will the authors be held liable for any damages arising from the use of this software.

Permission is granted to anyone to use this software for any purpose, including commercial applications, and to alter it and redistribute it freely, subject to the following restrictions:

1. The origin of this software must not be misrepresented; you must not claim that you wrote the original software. If you use this software in a product, an acknowledgment in the product documentation would be appreciated but is not required.

2. Altered source versions must be plainly marked as such, and must not be misrepresented as being the original software.

3. This notice may not be removed or altered from any source distribution

Author: Jonathan de Halleux, dehalleux@pelikhan.com, 2003
*/

#ifndef ZIPSTREAM_HPP
#define ZIPSTREAM_HPP

#include <vector>
#include <iostream>
#include <algorithm>
#include <sstream>
#include <zlib.h>
//#include <zlib/zlib.h>
//#include <zlib/zutil.h>

#ifdef HAVE_ZUTIL_H
#include <zutil.h>
#else
#ifdef MSDOS
#  define OS_CODE  0x00
#  if defined(__TURBOC__) || defined(__BORLANDC__)
#    if(__STDC__ == 1) && (defined(__LARGE__) || defined(__COMPACT__))
       /* Allow compilation with ANSI keywords only enabled */
       void _Cdecl farfree( void *block );
       void *_Cdecl farmalloc( unsigned long nbytes );
#    else
#     include <alloc.h>
#    endif
#  else /* MSC or DJGPP */
#    include <malloc.h>
#  endif
#endif

#ifdef OS2
#  define OS_CODE  0x06
#endif

#ifdef WIN32 /* Window 95 & Windows NT */
#  define OS_CODE  0x0b
#endif

#if defined(VAXC) || defined(VMS)
#  define OS_CODE  0x02
#  define F_OPEN(name, mode) \
     fopen((name), (mode), "mbc=60", "ctx=stm", "rfm=fix", "mrs=512")
#endif

#ifdef AMIGA
#  define OS_CODE  0x01
#endif

#if defined(ATARI) || defined(atarist)
#  define OS_CODE  0x05
#endif

#if defined(MACOS) || defined(TARGET_OS_MAC)
#  define OS_CODE  0x07
#  if defined(__MWERKS__) && __dest_os != __be_os && __dest_os != __win32_os
#    include <unix.h> /* for fdopen */
#  else
#    ifndef fdopen
#      define fdopen(fd,mode) NULL /* No fdopen() */
#    endif
#  endif
#endif

#ifdef __50SERIES /* Prime/PRIMOS */
#  define OS_CODE  0x0F
#endif

#ifdef TOPS20
#  define OS_CODE  0x0a
#endif

#if defined(_BEOS_) || defined(RISCOS)
#  define fdopen(fd,mode) NULL /* No fdopen() */
#endif

#if (defined(_MSC_VER) && (_MSC_VER > 600))
#  define fdopen(fd,type)  _fdopen(fd,type)
#endif


        /* Common defaults */

#ifndef OS_CODE
#  define OS_CODE  0x03  /* assume Unix */
#endif

#endif

namespace zlib_stream{

/// default gzip buffer size,
/// change this to suite your needs
const size_t default_buffer_size = 4096;

/// Compression strategy, see zlib doc.
enum EStrategy
{
	StrategyFiltered = 1,
	StrategyHuffmanOnly = 2,
	DefaultStrategy = 0
};

/** \brief A stream decorator that takes raw input and zips it to a ostream.

The class wraps up the inflate method of the zlib library 1.1.4 http://www.gzip.org/zlib/
*/
template<
	typename Elem, 
	typename Tr = std::char_traits<Elem>,
    typename ElemA = std::allocator<Elem>,
    typename ByteT = unsigned char,
    typename ByteAT = std::allocator<ByteT>
>	
class basic_zip_streambuf : public std::basic_streambuf<Elem, Tr> 
{
public:
	typedef std::basic_ostream<Elem, Tr>& ostream_reference;
    typedef ElemA char_allocator_type;
	typedef ByteT byte_type;
    typedef ByteAT byte_allocator_type;
	typedef byte_type* byte_buffer_type;
	typedef std::vector<byte_type, byte_allocator_type > byte_vector_type;
	typedef std::vector<char_type, char_allocator_type > char_vector_type;

    /** Construct a zip stream
     * More info on the following parameters can be found in the zlib documentation.
     */
    basic_zip_streambuf(
		ostream_reference ostream_,
		size_t level_,
		EStrategy strategy_,
		size_t window_size_,
		size_t memory_level_,
		size_t buffer_size_
		);
	
	~basic_zip_streambuf();

	int sync ();
    int_type overflow (int_type c);

	/** flushes the zip buffer and output buffer.

	This method should be called at the end of the compression. Calling flush multiple times, will lower the
	compression ratio.
	*/
	std::streamsize flush();
	/// returns a reference to the output stream
	ostream_reference get_ostream() const	{	return m_ostream;};
	/// returns the latest zlib error status
	int get_zerr() const					{	return m_err;};
	/// returns the crc of the input data compressed so far.
	long get_crc() const					{	return m_crc;};
	/// returns the size (bytes) of the input data compressed so far.
	long get_in_size() const				{	return m_zip_stream.total_in;};
	/// returns the size (bytes) of the compressed data so far.
	long get_out_size() const				{	return m_zip_stream.total_out;};
private:
	bool zip_to_stream( char_type*, std::streamsize);
	size_t fill_input_buffer();

	ostream_reference m_ostream;
	z_stream m_zip_stream;
    int m_err;
	byte_vector_type m_output_buffer;
	char_vector_type m_buffer; 
	long m_crc;
};

/** \brief A stream decorator that takes compressed input and unzips it to a istream.

The class wraps up the deflate method of the zlib library 1.1.4 http://www.gzip.org/zlib/
*/
template<
	typename Elem, 
	typename Tr = std::char_traits<Elem>,
    typename ElemA = std::allocator<Elem>,
    typename ByteT = unsigned char,
    typename ByteAT = std::allocator<ByteT>
>	
class basic_unzip_streambuf : 
	public std::basic_streambuf<Elem, Tr> 
{
public:
	typedef std::basic_istream<Elem, Tr>& istream_reference;
    typedef ElemA char_allocator_type;
	typedef ByteT byte_type;
    typedef ByteAT byte_allocator_type;
	typedef byte_type* byte_buffer_type;
	typedef std::vector<byte_type, byte_allocator_type > byte_vector_type;
	typedef std::vector<char_type, char_allocator_type > char_vector_type;

     /** Construct a unzip stream
     * More info on the following parameters can be found in the zlib documentation.
     */
	 basic_unzip_streambuf(
		istream_reference istream_,
		size_t window_size_,
		size_t read_buffer_size_,
		size_t input_buffer_size_
		);
	
	~basic_unzip_streambuf();

    int_type underflow();


	/// returns the compressed input istream
	istream_reference get_istream()	{	return m_istream;};
	/// returns the zlib stream structure
	z_stream& get_zip_stream()		{	return m_zip_stream;};
	/// returns the latest zlib error state
	int get_zerr() const					{	return m_err;};
	/// returns the crc of the uncompressed data so far 
	long get_crc() const					{	return m_crc;};
	/// returns the number of uncompressed bytes
	long get_out_size() const				{	return m_zip_stream.total_out;};
	/// returns the number of read compressed bytes
	long get_in_size() const				{	return m_zip_stream.total_in;};
private:
	void put_back_from_zip_stream();
	std::streamsize unzip_from_stream( char_type*, std::streamsize);

	size_t fill_input_buffer();

	istream_reference m_istream;
	z_stream m_zip_stream;
    int m_err;
	byte_vector_type m_input_buffer;
	char_vector_type m_buffer; 
	long m_crc;
};

/*! \brief Base class for zip ostreams

Contains a basic_zip_streambuf.
*/
template<
	typename Elem, 
	typename Tr = std::char_traits<Elem>,
    typename ElemA = std::allocator<Elem>,
    typename ByteT = unsigned char,
    typename ByteAT = std::allocator<ByteT>
>	
class basic_zip_ostreambase : virtual public std::basic_ios<Elem,Tr>
{
public:
	typedef std::basic_ostream<Elem, Tr>& ostream_reference;
	typedef basic_zip_streambuf<
        Elem,
        Tr,
        ElemA,
        ByteT,
        ByteAT
        > zip_streambuf_type;

    /** Construct a zip stream
     * More info on the following parameters can be found in the zlib documentation.
     */
	basic_zip_ostreambase( 
		ostream_reference ostream_,
		size_t level_,
		EStrategy strategy_,
		size_t window_size_,
		size_t memory_level_,
		size_t buffer_size_
		)
		: m_buf(ostream_,level_,strategy_,window_size_,memory_level_,buffer_size_)
	{
		init(&m_buf );
	};
	
	/// returns the underlying zip ostream object
	zip_streambuf_type* rdbuf() { return &m_buf; };

	/// returns the zlib error state
	int get_zerr() const					{	return m_buf.get_err();};
	/// returns the uncompressed data crc
	long get_crc() const					{	return m_buf.get_crc();};
	/// returns the compressed data size
	long get_out_size() const				{	return m_buf.get_out_size();};
	/// returns the uncompressed data size
	long get_in_size() const				{	return m_buf.get_in_size();};
private:
	zip_streambuf_type m_buf;
};

/*! \brief Base class for unzip istreams

Contains a basic_unzip_streambuf.
*/
template<
	typename Elem, 
	typename Tr = std::char_traits<Elem>,
    typename ElemA = std::allocator<Elem>,
    typename ByteT = unsigned char,
    typename ByteAT = std::allocator<ByteT>
>
class basic_zip_istreambase : virtual public std::basic_ios<Elem,Tr>
{
public:
	typedef std::basic_istream<Elem, Tr>& istream_reference;
	typedef basic_unzip_streambuf<
        Elem,
        Tr,
        ElemA,
        ByteT,
        ByteAT
        > unzip_streambuf_type;

	basic_zip_istreambase( 
		istream_reference ostream_,
		size_t window_size_,
		size_t read_buffer_size_,
		size_t input_buffer_size_
		)
		: m_buf(ostream_,window_size_, read_buffer_size_, input_buffer_size_)
	{
		init(&m_buf );
	};
	
	/// returns the underlying unzip istream object
	unzip_streambuf_type* rdbuf() { return &m_buf; };

	/// returns the zlib error state
	int get_zerr() const					{	return m_buf.get_zerr();};
	/// returns the uncompressed data crc
	long get_crc() const					{	return m_buf.get_crc();};
	/// returns the uncompressed data size
	long get_out_size() const				{	return m_buf.get_out_size();};
	/// returns the compressed data size
	long get_in_size() const				{	return m_buf.get_in_size();};
private:
	unzip_streambuf_type m_buf;
};

/*! \brief A zipper ostream

This class is a ostream decorator that behaves 'almost' like any other ostream.

At construction, it takes any ostream that shall be used to output of the compressed data.

When finished, you need to call the special method zflush or call the destructor 
to flush all the intermidiate streams.

Example:
\code
// creating the target zip string, could be a fstream
ostringstream ostringstream_;
// creating the zip layer
zip_ostream zipper(ostringstream_);

	
// writing data	
zipper<<f<<" "<<d<<" "<<ui<<" "<<ul<<" "<<us<<" "<<c<<" "<<dum;
// zip ostream needs special flushing...
zipper.zflush();
\endcode
*/
template<
	typename Elem, 
	typename Tr = std::char_traits<Elem>,
    typename ElemA = std::allocator<Elem>,
    typename ByteT = unsigned char,
    typename ByteAT = std::allocator<ByteT>
>	
class basic_zip_ostream : 
	public basic_zip_ostreambase<Elem,Tr,ElemA,ByteT,ByteAT>, 
	public std::basic_ostream<Elem,Tr>
{
public:
	typedef basic_zip_ostreambase<
        Elem,Tr,ElemA,ByteT,ByteAT> zip_ostreambase_type;
	typedef std::basic_ostream<Elem,Tr> ostream_type;

	/** Constructs a zipper ostream decorator
	 *
	 * \param ostream_ ostream where the compressed output is written
	 * \param is_gzip_ true if gzip header and footer have to be added
	 * \param level_ level of compression 0, bad and fast, 9, good and slower,
	 * \param strategy_ compression strategy
	 * \param window_size_ see zlib doc
	 * \param memory_level_ see zlib doc
	 * \param buffer_size_ the buffer size used to zip data

	 When is_gzip_ is true, a gzip header and footer is automatically added.
	 */
	basic_zip_ostream( 
		ostream_reference ostream_, 
        int open_mode = std::ios::out, 
		bool is_gzip_ = false,
		size_t level_ = Z_DEFAULT_COMPRESSION,
		EStrategy strategy_ = DefaultStrategy,
		size_t window_size_ = 15,
		size_t memory_level_ = 8,
		size_t buffer_size_ = default_buffer_size
		)
	: 
		zip_ostreambase_type(
            ostream_,
            level_,
            strategy_,
            window_size_,
            memory_level_,
            buffer_size_
            ), 
		m_is_gzip(is_gzip_),
		ostream_type(rdbuf())
	{
		if (m_is_gzip)
			add_header();
	};
	~basic_zip_ostream()
	{
		if (m_is_gzip)
			add_footer();
	}

	/// returns true if it is a gzip 
	bool is_gzip() const		{	return m_is_gzip;};
	/// flush inner buffer and zipper buffer
	basic_zip_ostream<Elem,Tr>& zflush()	
	{	
		flush(); rdbuf()->flush(); return *this; 
	};

private:
    static void put_long(ostream_reference out_, unsigned long x_);

	void add_header();
	void add_footer();
	bool m_is_gzip;
};

/*! \brief A zipper istream

This class is a istream decorator that behaves 'almost' like any other ostream.

At construction, it takes any istream that shall be used to input of the compressed data.

Simlpe example:
\code
// create a stream on zip string
istringstream istringstream_( ostringstream_.str());
// create unzipper istream
zip_istream unzipper( istringstream_);

// read and unzip
unzipper>>f_r>>d_r>>ui_r>>ul_r>>us_r>>c_r>>dum_r;
\endcode
*/
template<
	typename Elem, 
	typename Tr = std::char_traits<Elem>,
    typename ElemA = std::allocator<Elem>,
    typename ByteT = unsigned char,
    typename ByteAT = std::allocator<ByteT>
>
class basic_zip_istream : 
	public basic_zip_istreambase<Elem,Tr,ElemA,ByteT,ByteAT>, 
	public std::basic_istream<Elem,Tr>
{
public:
	typedef basic_zip_istreambase<
        Elem,Tr,ElemA,ByteT,ByteAT> zip_istreambase_type;
	typedef std::basic_istream<Elem,Tr> istream_type;
	typedef unsigned char byte_type;

	/** Construct a unzipper stream
	 *
	 * \param istream_ input buffer
	 * \param window_size_ 
	 * \param read_buffer_size_ 
	 * \param input_buffer_size_ 
	 */
	basic_zip_istream( 
		istream_reference istream_, 
		size_t window_size_ = 15,
		size_t read_buffer_size_ = default_buffer_size,
		size_t input_buffer_size_ = default_buffer_size
		)
	  : 
		zip_istreambase_type(istream_,window_size_, read_buffer_size_, input_buffer_size_), 
		istream_type(rdbuf()),
		m_is_gzip(false),
		m_gzip_crc(0),
		m_gzip_data_size(0)
	{
 	      if (rdbuf()->get_zerr()==Z_OK)
			  check_header();
	};

	/// returns true if it is a gzip file
	bool is_gzip() const				{	return m_is_gzip;};
	/// reads the gzip header
	void read_footer();
	/** return crc check result

	When you have finished reading the compressed data, call read_footer to read the uncompressed data crc.
	This method compares it to the crc of the uncompressed data.

	\return true if crc check is succesful 
	*/
	bool check_crc() const				{	return get_crc() == m_gzip_crc;};
	/// return data size check
	bool check_data_size() const		{	return get_out_size() == m_gzip_data_size;};

	/// return the crc value in the file
	long get_gzip_crc() const			{	return m_gzip_crc;};
	/// return the data size in the file 
	long get_gzip_data_size() const		{	return m_gzip_data_size;};
protected:
    static void read_long(istream_reference in_, unsigned long& x_);

	int check_header();
	bool m_is_gzip;
	unsigned long m_gzip_crc;
	unsigned long m_gzip_data_size;
};

/// A typedef for basic_zip_ostream<char>
typedef basic_zip_ostream<char> zip_ostream;
/// A typedef for basic_zip_ostream<wchar_t>
typedef basic_zip_ostream<wchar_t> zip_wostream;
/// A typedef for basic_zip_istream<char>
typedef basic_zip_istream<char> zip_istream;
/// A typedef for basic_zip_istream<wchart>
typedef basic_zip_istream<wchar_t> zip_wistream;

}; // zlib_sream



namespace zlib_stream{

namespace detail{
	const int gz_magic[2] = {0x1f, 0x8b}; /* gzip magic header */

	/* gzip flag byte */
	const int gz_ascii_flag =  0x01; /* bit 0 set: file probably ascii text */
	const int gz_head_crc    = 0x02; /* bit 1 set: header CRC present */
	const int gz_extra_field = 0x04; /* bit 2 set: extra field present */
	const int gz_orig_name  =  0x08; /* bit 3 set: original file name present */
	const int gz_comment    =  0x10; /* bit 4 set: file comment present */
	const int gz_reserved   =  0xE0; /* bits 5..7: reserved */	
}

	template<
		typename Elem, 
		typename Tr,
		typename ElemA,
		typename ByteT,
		typename ByteAT
	>
	basic_zip_streambuf<
		Elem,Tr,ElemA,ByteT,ByteAT
		>::basic_zip_streambuf(
		ostream_reference ostream_,
		size_t level_,
		EStrategy strategy_,
		size_t window_size_,
		size_t memory_level_,
		size_t buffer_size_
	)
	: 
		m_ostream(ostream_),
		m_output_buffer(buffer_size_,0),
		m_buffer(buffer_size_,0),
		m_crc(0)
	{
		m_zip_stream.zalloc=(alloc_func)0;
		m_zip_stream.zfree=(free_func)0;

		m_zip_stream.next_in=NULL;
		m_zip_stream.avail_in=0;
		m_zip_stream.avail_out=0;
		m_zip_stream.next_out=NULL;

		m_err=deflateInit2(
			&m_zip_stream, 
			std::min( 9, static_cast<int>(level_)),
			Z_DEFLATED,
			- static_cast<int>(window_size_), // <-- changed
			std::min( 9, static_cast<int>(memory_level_) ),
			static_cast<int>(strategy_)
			);
			
		setp( &(m_buffer[0]), &(m_buffer[m_buffer.size()-1]));
	};

	template<
		typename Elem, 
		typename Tr,
		typename ElemA,
		typename ByteT,
		typename ByteAT
	>
	basic_zip_streambuf<
		Elem,Tr,ElemA,ByteT,ByteAT
		>::~basic_zip_streambuf()
	{
		flush();
		m_ostream.flush();
		m_err=deflateEnd(&m_zip_stream);
	};

	template<
		typename Elem, 
		typename Tr,
		typename ElemA,
		typename ByteT,
		typename ByteAT
	>
	int basic_zip_streambuf<
		Elem,Tr,ElemA,ByteT,ByteAT
		>::sync ()
	{ 
		if ( pptr() && pptr() > pbase()) 
		{
			int c = overflow( EOF);

			if ( c == EOF)
				return -1;
        }

        return 0;
	}

	template<
		typename Elem, 
		typename Tr,
		typename ElemA,
		typename ByteT,
		typename ByteAT
	>
	typename basic_zip_streambuf<
		Elem,Tr,ElemA,ByteT,ByteAT
		>::int_type 
		basic_zip_streambuf<
			Elem,Tr,ElemA,ByteT,ByteAT
			>::overflow (
			typename basic_zip_streambuf<
				Elem,Tr,ElemA,ByteT,ByteAT
				>::int_type c
			)
	{ 
        int w = static_cast<int>(pptr() - pbase());
        if (c != EOF) {
             *pptr() = c;
             ++w;
         }
         if ( zip_to_stream( pbase(), w)) {
             setp( pbase(), epptr() - 1);
             return c;
         } else
             return EOF;
	}
	
	template<
		typename Elem, 
		typename Tr,
		typename ElemA,
		typename ByteT,
		typename ByteAT
	>
	bool basic_zip_streambuf<
		Elem,Tr,ElemA,ByteT,ByteAT
		>::zip_to_stream( 
		typename basic_zip_streambuf<
			Elem,Tr,ElemA,ByteT,ByteAT
			>::char_type* buffer_, 
			std::streamsize buffer_size_
		)
	{	
		std::streamsize written_byte_size=0, total_written_byte_size = 0;

		m_zip_stream.next_in=(byte_buffer_type)buffer_;
		m_zip_stream.avail_in=static_cast<uInt>(buffer_size_*sizeof(char_type));
		m_zip_stream.avail_out=static_cast<uInt>(m_output_buffer.size());
		m_zip_stream.next_out=&(m_output_buffer[0]);
		size_t remainder=0;

		// updating crc
		m_crc = crc32( 
			m_crc, 
			m_zip_stream.next_in,
			m_zip_stream.avail_in
			);		

		do
		{
			m_err = deflate(&m_zip_stream, 0);
	
			if (m_err == Z_OK  || m_err == Z_STREAM_END)
			{
				written_byte_size= 
					static_cast<std::streamsize>(m_output_buffer.size()) 
					- m_zip_stream.avail_out;
				total_written_byte_size+=written_byte_size;
				// ouput buffer is full, dumping to ostream
				m_ostream.write( 
					(const char_type*) &(m_output_buffer[0]), 
					static_cast<std::streamsize>( 
						written_byte_size/sizeof(char_type) 
						)
					);
												
				// checking if some bytes were not written.
				if ( (remainder = written_byte_size%sizeof(char_type))!=0)
				{
					// copy to the beginning of the stream
					memcpy(
						&(m_output_buffer[0]), 
						&(m_output_buffer[written_byte_size-remainder]),
						remainder);
					
				}
				
				m_zip_stream.avail_out=
					static_cast<uInt>(m_output_buffer.size()-remainder);
				m_zip_stream.next_out=&m_output_buffer[remainder];
			}
		} 
		while (m_zip_stream.avail_in != 0 && m_err == Z_OK);
	
		return m_err == Z_OK;
	};

	template<
		typename Elem, 
		typename Tr,
		typename ElemA,
		typename ByteT,
		typename ByteAT
	>
	std::streamsize basic_zip_streambuf<
		Elem,Tr,ElemA,ByteT,ByteAT
		>::flush()
	{
		std::streamsize written_byte_size=0, total_written_byte_size=0;

		size_t remainder=0;

		// updating crc
		m_crc = crc32( 
			m_crc, 
			m_zip_stream.next_in,
			m_zip_stream.avail_in
			);		

		do
		{
			m_err = deflate(&m_zip_stream, Z_FINISH);
			if (m_err == Z_OK || m_err == Z_STREAM_END)
			{
				written_byte_size=
					static_cast<std::streamsize>(m_output_buffer.size()) 
					- m_zip_stream.avail_out;
				total_written_byte_size+=written_byte_size;
				// ouput buffer is full, dumping to ostream
				m_ostream.write( 
					(const char_type*) &(m_output_buffer[0]), 
					static_cast<std::streamsize>( 
						written_byte_size/sizeof(char_type)*sizeof(byte_type) 
						)
					);
			
				// checking if some bytes were not written.
				if ( (remainder = written_byte_size%sizeof(char_type))!=0)
				{
					// copy to the beginning of the stream
					memcpy(
						&(m_output_buffer[0]), 
						&(m_output_buffer[written_byte_size-remainder]),
						remainder);
					
				}
				
				m_zip_stream.avail_out=static_cast<uInt>(m_output_buffer.size()-remainder);
				m_zip_stream.next_out=&m_output_buffer[remainder];
			}
		} while (m_err == Z_OK);

		m_ostream.flush();

		return total_written_byte_size;
	};


	template<
		typename Elem, 
		typename Tr,
		typename ElemA,
		typename ByteT,
		typename ByteAT
	>
    basic_unzip_streambuf<
		Elem,Tr,ElemA,ByteT,ByteAT
		>::basic_unzip_streambuf(
			istream_reference istream_,
			size_t window_size_,
			size_t read_buffer_size_,
			size_t input_buffer_size_
	)
	:  
		m_istream(istream_),
		m_input_buffer(input_buffer_size_),
		m_buffer(read_buffer_size_),
		m_crc(0)
	{
		// setting zalloc, zfree and opaque
		m_zip_stream.zalloc=(alloc_func)0;
		m_zip_stream.zfree=(free_func)0;

		m_zip_stream.next_in=NULL;
		m_zip_stream.avail_in=0;
		m_zip_stream.avail_out=0;
		m_zip_stream.next_out=NULL;
	
		m_err=inflateInit2(
			&m_zip_stream,
			-static_cast<int>(window_size_)
		);
		
		setg( &(m_buffer[0])+4,     // beginning of putback area
		    &(m_buffer[0])+4,     // read position
	        &(m_buffer[0])+4);    // end position    
	};

	template<
		typename Elem, 
		typename Tr,
		typename ElemA,
		typename ByteT,
		typename ByteAT
	>
    size_t basic_unzip_streambuf<
		Elem,Tr,ElemA,ByteT,ByteAT
		>::fill_input_buffer()
	{
		m_zip_stream.next_in=&(m_input_buffer[0]);
		m_istream.read( 
			(char_type*)(&(m_input_buffer[0])), 
			static_cast<std::streamsize>(m_input_buffer.size()/sizeof(char_type)) 
			); 
		return m_zip_stream.avail_in=m_istream.gcount()*sizeof(char_type);
	}


	template<
		typename Elem, 
		typename Tr,
		typename ElemA,
		typename ByteT,
		typename ByteAT
	>
    basic_unzip_streambuf<
		Elem,Tr,ElemA,ByteT,ByteAT
		>::~basic_unzip_streambuf()
	{
		inflateEnd(&m_zip_stream);
	};

	template<
		typename Elem, 
		typename Tr,
		typename ElemA,
		typename ByteT,
		typename ByteAT
	>
    typename basic_unzip_streambuf<
		Elem,Tr,ElemA,ByteT,ByteAT
		>::int_type 
		basic_unzip_streambuf<
			Elem,Tr,ElemA,ByteT,ByteAT
			>::underflow() 
	{ 
		if ( gptr() && ( gptr() < egptr()))
			return * reinterpret_cast<unsigned char *>( gptr());
     
       int n_putback = static_cast<int>(gptr() - eback());
       if ( n_putback > 4)
          n_putback = 4;
       memcpy( 
			&(m_buffer[0]) + (4 - n_putback), 
			gptr() - n_putback, 
			n_putback*sizeof(char_type)
			);
  
	   int num = unzip_from_stream( 
		   &(m_buffer[0])+4, 
		   static_cast<std::streamsize>((m_buffer.size()-4)*sizeof(char_type))
		   );
        if (num <= 0) // ERROR or EOF
           return EOF;
    
        // reset buffer pointers
        setg( &(m_buffer[0]) + (4 - n_putback),   // beginning of putback area
              &(m_buffer[0]) + 4,                 // read position
              &(m_buffer[0]) + 4 + num);          // end of buffer
    
         // return next character
         return* reinterpret_cast<unsigned char *>( gptr());    
     }

	template<
		typename Elem, 
		typename Tr,
		typename ElemA,
		typename ByteT,
		typename ByteAT
	>
    std::streamsize basic_unzip_streambuf<
		Elem,Tr,ElemA,ByteT,ByteAT
		>::unzip_from_stream( 
			char_type* buffer_, 
			std::streamsize buffer_size_
			)
	{
		m_zip_stream.next_out=(byte_buffer_type)buffer_;
		m_zip_stream.avail_out=static_cast<uInt>(buffer_size_*sizeof(char_type));
		size_t count =m_zip_stream.avail_in;

		do
		{
			if (m_zip_stream.avail_in==0)
				count=fill_input_buffer();

			if (m_zip_stream.avail_in)
			{
				m_err = inflate( &m_zip_stream,  Z_SYNC_FLUSH );
			}
		} while (m_err==Z_OK && m_zip_stream.avail_out != 0 && count != 0);

		// updating crc
		m_crc = crc32( 
			m_crc, 
			(byte_buffer_type)buffer_,
			buffer_size_ - m_zip_stream.avail_out/sizeof(char_type)
			);	
		std::streamsize n_read = buffer_size_ - m_zip_stream.avail_out/sizeof(char_type);
		
		// check if it is the end
		if (m_err==Z_STREAM_END)
			put_back_from_zip_stream();				
		
		return n_read;
	}
	
	template<
		typename Elem, 
		typename Tr,
		typename ElemA,
		typename ByteT,
		typename ByteAT
	>
	void basic_unzip_streambuf<
		Elem,Tr,ElemA,ByteT,ByteAT
		>::put_back_from_zip_stream()
	{
		if (m_zip_stream.avail_in==0)
			return;

		m_istream.clear( ios::goodbit );
		m_istream.seekg(
			-static_cast<int>(m_zip_stream.avail_in),
			ios_base::cur
			);

		m_zip_stream.avail_in=0;
	};

	template<
		typename Elem, 
		typename Tr,
		typename ElemA,
		typename ByteT,
		typename ByteAT
	>
	int basic_zip_istream<
		Elem,Tr,ElemA,ByteT,ByteAT
		>::check_header()
	{
	    int method; /* method byte */
	    int flags;  /* flags byte */
	    uInt len;
		int c;
		int err=0;
		z_stream& zip_stream = rdbuf()->get_zip_stream();

	    /* Check the gzip magic header */
		 for (len = 0; len < 2; len++) 
		 {
			c = (int)rdbuf()->get_istream().get();
			if (c != detail::gz_magic[len]) 
			{
			    if (len != 0) 
					rdbuf()->get_istream().unget();
			    if (c!= EOF) 
			    {
					rdbuf()->get_istream().unget();
			    }
		    
			    err = zip_stream.avail_in != 0 ? Z_OK : Z_STREAM_END;
			    m_is_gzip = false;
			    return err;
			}
		}
    
		m_is_gzip = true;
		method = (int)rdbuf()->get_istream().get();
		flags = (int)rdbuf()->get_istream().get();
		if (method != Z_DEFLATED || (flags & detail::gz_reserved) != 0) 
		{
			err = Z_DATA_ERROR;
			return err;
		}

	    /* Discard time, xflags and OS code: */
	    for (len = 0; len < 6; len++) 
			rdbuf()->get_istream().get();
	
	    if ((flags & detail::gz_extra_field) != 0) 
	    { 
			/* skip the extra field */
			len  =  (uInt)rdbuf()->get_istream().get();
			len += ((uInt)rdbuf()->get_istream().get())<<8;
			/* len is garbage if EOF but the loop below will quit anyway */
			while (len-- != 0 && rdbuf()->get_istream().get() != EOF) ;
	    }
	    if ((flags & detail::gz_orig_name) != 0) 
	    { 
			/* skip the original file name */
			while ((c = rdbuf()->get_istream().get()) != 0 && c != EOF) ;
		}
	    if ((flags & detail::gz_comment) != 0) 
	    {   
			/* skip the .gz file comment */
			while ((c = rdbuf()->get_istream().get()) != 0 && c != EOF) ;
		}
		if ((flags & detail::gz_head_crc) != 0) 
		{  /* skip the header crc */
			for (len = 0; len < 2; len++) 
				rdbuf()->get_istream().get();
		}
		err = rdbuf()->get_istream().eof() ? Z_DATA_ERROR : Z_OK;

		return err;
	}

	template<
		typename Elem, 
		typename Tr,
		typename ElemA,
		typename ByteT,
		typename ByteAT
	>
	void basic_zip_istream<
		Elem,Tr,ElemA,ByteT,ByteAT
		>::read_footer()
	{		
		if (m_is_gzip)
		{
			read_long( rdbuf()->get_istream(), m_gzip_crc );
			read_long( rdbuf()->get_istream(), m_gzip_data_size );
		}
	}

	template<
		typename Elem, 
		typename Tr,
		typename ElemA,
		typename ByteT,
		typename ByteAT
	>
    void basic_zip_ostream<
		Elem,Tr,ElemA,ByteT,ByteAT
		>::put_long( 
			typename basic_zip_ostream<
				Elem,Tr,ElemA,ByteT,ByteAT
				>::ostream_reference out_, 
			unsigned long x_
			)
    {
		static const int size_ul = sizeof(unsigned long);
		static const int size_c = sizeof(char_type);
        static const int n_end = size_ul/size_c;
		out_.write(reinterpret_cast<char_type const*>(&x_), n_end);
    }
   
	template<
		typename Elem, 
		typename Tr,
		typename ElemA,
		typename ByteT,
		typename ByteAT
	>
    void basic_zip_istream<
			Elem,Tr,ElemA,ByteT,ByteAT
			>::read_long(
				istream_reference in_, 
			unsigned long& x_
			)
	{
		static const int size_ul = sizeof(unsigned long);
		static const int size_c = sizeof(char_type);
        static const int n_end = size_ul/size_c;
		in_.read(reinterpret_cast<char*>(&x_),n_end);
	}
    
	template<
		typename Elem, 
		typename Tr,
		typename ElemA,
		typename ByteT,
		typename ByteAT
	>
	void basic_zip_ostream<
		Elem,Tr,ElemA,ByteT,ByteAT
		>::add_header()
	{
	    char_type zero=0;
	    
        rdbuf()->get_ostream()
			.put(static_cast<char_type>(detail::gz_magic[0]))
			.put(static_cast<char_type>(detail::gz_magic[1]))
			.put(static_cast<char_type>(Z_DEFLATED))
			.put(zero) //flags
			.put(zero).put(zero).put(zero).put(zero) // time
			.put(zero) //xflags
			.put(static_cast<char_type>(OS_CODE));
	}

	template<
		typename Elem, 
		typename Tr,
		typename ElemA,
		typename ByteT,
		typename ByteAT
	>
	void basic_zip_ostream<
		Elem,Tr,ElemA,ByteT,ByteAT
		>::add_footer()
	{
		put_long( rdbuf()->get_ostream(), rdbuf()->get_crc() );
		put_long( rdbuf()->get_ostream(), rdbuf()->get_in_size() ); 
	};

}; // zlib_sream


#endif

