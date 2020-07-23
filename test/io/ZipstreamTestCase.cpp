/*
 * Copyright (c) 2004-2020 The University of Notre Dame. All Rights Reserved.
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
 * [1] Meineke, et al., J. Comp. Chem. 26, 252-271 (2005).
 * [2] Fennell & Gezelter, J. Chem. Phys. 124, 234104 (2006).
 * [3] Sun, Lin & Gezelter, J. Chem. Phys. 128, 234107 (2008).
 * [4] Vardeman, Stocker & Gezelter, J. Chem. Theory Comput. 7, 834 (2011).
 * [5] Kuang & Gezelter, Mol. Phys., 110, 691-701 (2012).
 * [6] Lamichhane, Gezelter & Newman, J. Chem. Phys. 141, 134109 (2014).
 * [7] Lamichhane, Newman & Gezelter, J. Chem. Phys. 141, 134110 (2014).
 * [8] Bhattarai, Newman & Gezelter, Phys. Rev. B 99, 094106 (2019).
 */

#include <iostream>
#include <string>

#include "io/basic_ifstrstream.hpp"
#include "ZipstreamTestCase.hpp"

// Registers the fixture into the 'registry'

using namespace std;
using namespace OpenMD;
CPPUNIT_TEST_SUITE_REGISTRATION( ZipstreamTeseCase );

#include "zipstream.hpp"
#include <iostream>
#include <string>
#include <sstream>
#include <fstream>
#include <cmath>
#include <cstdlib>

using namespace std;
using namespace zlib_stream;

// a dummy class
class dummy
{
public:
	dummy(float f_ = 3.14f, int i_= 1) : f(f_), i(i_)
	{};
	void reset(){f=0.0;i=0;};

	friend ostream& operator << ( ostream& out, dummy const& d)
	{
		out<<" "<<d.f<<" "<<d.i;
		return out;
	}
	friend istream& operator >> ( istream& in, dummy& d)
	{
		in>>d.f>>d.i;
		return in;
	}

	friend wostream& operator << ( wostream& out, dummy const& d)
	{
		out<<L" "<<d.f<<L" "<<d.i;
		return out;
	}
	friend wistream& operator >> ( wistream& in, dummy& d)
	{
		in>>d.f>>d.i;
		return in;
	}
protected:
	float f;
	int i;
};

void ZipstreamTeseCase::test_buffer_to_buffer()
{
	const size_t n= 1024;
	char in_buffer[n]={'\0'};
	char zip_buffer[n]={'\0'};

	for (size_t i=0;i<n-1;++i)
		in_buffer[i]=static_cast<char>(48+i%48);


	ostringstream out;
	zip_ostream zipper(out);

	zipper.write(in_buffer, n);
	if (zipper.fail())
		cerr<<"failed to write stream"<<endl;

	zipper.zflush();


	istringstream in( out.str() );
	zip_istream unzipper(in);
	unzipper.read(zip_buffer, n);

	cerr<<"buffers equals: "<<endl
		<<"-----------------"<<endl
		 <<"\t crc source:"<<crc32(0L,(unsigned char*)in_buffer,n)<<endl
		 <<"\t crc result:"<<crc32(0L,(unsigned char*)zip_buffer,n)<<endl;

}

void ZipstreamTeseCase::test_wbuffer_to_wbuffer()
{
	const size_t n= 1024;
	wchar_t in_buffer[n]={'\0'};
	wchar_t zip_buffer[n]={'\0'};

	for (size_t i=0;i<n-1;++i)
		in_buffer[i]=static_cast<wchar_t>(48+i%48);


	wostringstream out;
	zip_wostream zipper(out);

	zipper.write(in_buffer, n);
	if (zipper.fail())
		cerr<<"failed to write stream"<<endl;

	zipper.zflush();


	wistringstream in( out.str() );
	zip_wistream unzipper(in);
	unzipper.read(zip_buffer, n);

	wcerr<<L"buffers equals: "<<endl
		 <<L"-----------------"<<endl
		 <<L"\t crc source:"<<crc32(0L,(unsigned char*)in_buffer,n*sizeof(wchar_t))<<endl
		 <<L"\t crc result:"<<crc32(0L,(unsigned char*)zip_buffer,n*sizeof(wchar_t))<<endl;
}

void ZipstreamTeseCase::test_string_string()
{
	// create some test values
	char c_r='-',c= 'a';
	string s_r="",s = "std::string";
	double d_r=0,d = asin( static_cast<double>(1.0) ) *2.0;
	float f_r=0, f = static_cast<float>(asin( 1.0 ) *2.0);
	unsigned int ui_r=0,ui = rand();
	unsigned long ul_r=0,ul = rand();
	unsigned short us_r=0, us = rand();
	dummy dum,dum_r(0,0);

	/*----------------------------------------------------------------------
	
	Testing string to string zipping/unzipping
	
	-------------------------------------------------------------------------*/
	// creating the target zip string, could be a fstream
	ostringstream ostringstream_(ios::out );
	// creating the zip layer
	zip_ostream zipper(ostringstream_);

	// writing data
	zipper<<f<<" "<<d<<" "<<ui<<" "<<ul<<" "<<us<<" "<<c<<" "<<dum;
	// zip ostream needs special flushing...
	zipper.zflush();

	// create a stream on zip string
	istringstream istringstream_( ostringstream_.str(), ios::in);
	// create unzipper istream
	zip_istream unzipper( istringstream_);

	// unzipping
	unzipper>>f_r>>d_r>>ui_r>>ul_r>>us_r>>c_r>>dum_r;

	// ouputing results
	cerr<<"tests string-string,  char:"<<endl
		<<"----------------------------"<<endl
		<<"double : "<<d<<" "<<d_r<<endl
		<<"char : "<<c<<" "<<c_r<<endl
		<<"float : "<<f<<" "<<f_r<<endl
		<<"unsigned int : "<<ui<<" "<<ui_r<<endl
		<<"unsigned long : "<<ul<<" "<<ul_r<<endl
		<<"unsigned short : "<<us<<" "<<us_r<<endl
		<<"dummy class: "<<dum<<" "<<dum_r<<endl
		<<endl
		;
}

void ZipstreamTeseCase::test_wstring_wstring()
{
	// create some test values
	char c_r='-',c= 'a';
	double d_r=0,d = asin( 1.0 ) *2.0;
	float f_r=0, f = static_cast<float>(asin( 1.0 ) *2.0);
	unsigned int ui_r=0,ui = rand();
	unsigned long ul_r=0,ul = rand();
	unsigned short us_r=0, us = rand();
	dummy dum,dum_r(0,0);

	/*

	Testing wide characters

	*/
	f_r=0.0f;
	d_r=0.0;
	ui_r=ul_r=us_r=0;
	dum_r.reset();
	// creating the target zip string, could be a fstream
	wostringstream wostringstream_;
	// creating the zip layer
	zip_wostream wzipper(wostringstream_);

	// writing data
	wzipper<<f<<L" "<<d<<L" "<<ui<<L" "<<ul<<L" "<<us<<L" "<<dum;
	// zip ostream needs special flushing...
	wzipper.zflush();


	// create a stream on zip string
	wistringstream wistringstream_( wostringstream_.str());
	// create unzipper istream
	zip_wistream wunzipper( wistringstream_ );

	// unzipping
	wunzipper>>f_r>>d_r>>ui_r>>ul_r>>us_r>>dum_r;

	// ouputing results
	cerr<<"tests string-string (wchar_t):"<<endl
		<<"------------------------------"<<endl
		<<"double : "<<d<<" "<<d_r<<endl
		<<"float : "<<f<<" "<<f_r<<endl
		<<"unsigned int : "<<ui<<" "<<ui_r<<endl
		<<"unsigned long : "<<ul<<" "<<ul_r<<endl
		<<"unsigned short : "<<us<<" "<<us_r<<endl
		<<"dummy class: "<<dum<<" "<<dum_r<<endl
		<<endl
		;
}

void ZipstreamTeseCase::test_file_file()
{
      bool add_gzip_header = true;
	// create some test values
	char c_r='-',c= 'a';
	const char* sz = "const char*";
	string s_r="",s = "std::string";
	double d_r=0,d = asin( 1.0 ) *2.0;
	float f_r=0, f = static_cast<float>(asin( 1.0 ) *2.0);
	unsigned int ui_r=0,ui = rand();
	unsigned long ul_r=0,ul = rand();
	unsigned short us_r=0, us = rand();
	dummy dum,dum_r(0,0);
	char sbuf[256]={'-'};
	long crc_z(0), crc_uz(0);
	long in_size_z(0), in_size_uz(0), out_size_z(0), out_size_uz(0);

	/*----------------------------------------------------------------------------
	
	Testing file to file unzipping

	------------------------------------------------------------------------------*/
	f_r=0.0f;
	d_r=0.0;
	ui_r=ul_r=us_r=0; dum_r.reset();

	{
		// creating the target zip string, could be a fstream
        ofstream ofstream_("test.zip",ios::out | ios::binary );
		// creating the zip layer
        zip_ostream fzipper(ofstream_, add_gzip_header);
		// writing data
		fzipper<<f<<" "<<d<<" "<<ui<<" "<<ul<<" "<<us<<" "<<c<<" "<<dum;
		// zip ostream needs special flushing...
		fzipper.zflush();
		crc_z = fzipper.rdbuf()->get_crc();
		in_size_z = fzipper.rdbuf()->get_in_size();
		out_size_z = fzipper.rdbuf()->get_out_size();
	}

	// create a stream on zip string
	ifstream ifstream_;
	ifstream_.open("test.zip", ios::in | ios::binary);
	if (!ifstream_.is_open())
	{
		cerr<<"Could not open file test.zip"<<endl;
	}
	// create unzipper istream
	zip_istream funzipper( ifstream_);

	// unzipping
	funzipper>>f_r>>d_r>>ui_r>>ul_r>>us_r>>c_r>>dum_r;
	funzipper.read_footer();

	crc_uz = funzipper.rdbuf()->get_crc();
	in_size_uz = funzipper.rdbuf()->get_in_size();
	out_size_uz = funzipper.rdbuf()->get_out_size();


	// ouputing results
	cerr<<"tests file-to-file (char, "<<(funzipper.is_gzip() ? "gzip" : "no gzip")<<"):"<<endl
		<<"------------------------------"<<endl
		<<"double : "<<d<<" "<<d_r<<endl
		<<"char : "<<c<<" "<<c_r<<endl
		<<"float : "<<f<<" "<<f_r<<endl
		<<"unsigned int : "<<ui<<" "<<ui_r<<endl
		<<"unsigned long : "<<ul<<" "<<ul_r<<endl
		<<"unsigned short : "<<us<<" "<<us_r<<endl
		<<"dummy class: "<<dum<<" "<<dum_r<<endl
		<<"crc (zipped, unzipped, from file): "<<crc_z<<" "<<crc_uz<<" 
"<<funzipper.get_gzip_crc()<<endl
		<<"uncompressed data size: "<<in_size_z<<" "<<out_size_uz<<" 
"<<funzipper.get_gzip_data_size()<<endl
		<<"compressed data size: "<<out_size_z<<" "<<in_size_uz<<endl
		<<"check_crc: "<<( funzipper.check_crc() ? "ok" : "failed")<<endl
		<<"check_data_size: "<<( funzipper.check_data_size() ? "ok" : "failed")<<endl
		;
}

void ZipstreamTeseCase::test_crc()
{
int i;

{
    ofstream outFile("test.gz", ios::out | ios::binary);
    zip_ostream zipOut(outFile,true);
    char buff[102400];
    for (i = 0; i < 102400; i++)
    {
        buff[i] = (char)48+i%48;
    }
    zipOut.write(buff, sizeof(buff));
}

ifstream inFile("test.gz", ios::in | ios::binary);

zip_istream unzipper(inFile);

char c;

for (i = 0; i < 102400; i++)
{
    unzipper >> c;
}

unzipper.read_footer();

std::cout << "DoMyTest() " 
<< ", outsize " << unzipper.get_out_size()
<< ", insize " << unzipper.get_in_size() 
<< ", z_err " << unzipper.get_zerr() 
<< ", crc " << unzipper.get_crc() 
<< ", gzip data size " << unzipper.get_gzip_data_size() << std::endl;
}
