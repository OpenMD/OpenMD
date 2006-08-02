/**********************************************************************
main.cpp - Main conversion program, command-line handling.

Copyright (C) 1998-2001 by OpenEye Scientific Software, Inc.
Some portions Copyright (C) 2001-2005 by Geoffrey R. Hutchison
Some portions Copyright (C) 2004-2005 by Chris Morley

This file is part of the Open Babel project.
For more information, see <http://openbabel.sourceforge.net/>

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation version 2 of the License.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.
***********************************************************************/

#include "config.h"
#if HAVE_IOSTREAM
	#include <iostream>
#elif HAVE_IOSTREAM_H
	#include <iostream.h>
#endif
#if HAVE_FSTREAM
	#include <fstream>
#elif HAVE_FSTREAM_H
	#include <fstream.h>
#endif
#if HAVE_SSTREAM
	#include <sstream>
#elif
	#include <sstream.h>
#endif

#include <string>
#include <map>
#if HAVE_CONIO_H
	#include <conio.h>
#endif

#if !HAVE_STRNCASECMP
extern "C" int strncasecmp(const char *s1, const char *s2, size_t n);
#endif

#include "openbabel/obconversion.hpp"
#include "brains/Register.hpp"

using namespace std;
using namespace OpenBabel;
using namespace oopse;
void DoOption(const char* p, OBConversion& Conv, OBConversion::Option_type typ, 
							int& arg, int argc, char *argv[]); 
void usage();
void help();

// There isn't a great way to do this -- we need to save argv[0] for usage()
static char *program_name;

int main(int argc,char *argv[])
{
  registerOBFormats();	
  OBConversion Conv(&cin, &cout); //default input and output are console 

//  string GenOptions;
  OBFormat* pInFormat = NULL;
  OBFormat* pOutFormat = NULL;
  vector<string> FileList, OutputFileList;
  string OutputFileName;
//  obMessageLevel filterLevel = obWarning; // 2 out of 5

  // Parse commandline
  bool gotInType = false, gotOutType = false;
  bool UseSavedOptions = false;

  char *oext;
  char *iext;
  string inputExt;
  string outputExt;

  //Save name of program without its path (and .exe)
  string pn(argv[0]);
  unsigned int pos;
#ifdef _WIN32
  pos = pn.find(".exe");
  if(pos!=string::npos)
    argv[0][pos]='\0';
#endif
  pos = pn.find_last_of("/\\");
  if(pos==string::npos)
    program_name=argv[0];
  else
    program_name=argv[0]+pos+1;

	const char* p;
	int arg;
	for (arg = 1; arg < argc; arg++)
	{
		if (argv[arg])
		{
			if (argv[arg][0] == '-')
			{
				switch (argv[arg][1])
				{

				case 'V':
				{
					cout << "Open Babel " << BABEL_VERSION << " -- " 
				       << __DATE__ << " -- " << __TIME__ << endl;
					exit(0);
				}

				case 'i':
					gotInType = true;
					iext = argv[arg] + 2;
					if(!*iext)
						iext = argv[++arg]; //space left after -i: use next argument

					if (strncasecmp(iext, "MIME", 4) == 0)
					{
						// get the MIME type from the next argument
						iext = argv[++arg];
						pInFormat = Conv.FormatFromMIME(iext);
					}
					else
					{
						//The ID provided by the OBFormat class is used as the identifying file extension
						pInFormat = Conv.FindFormat(iext);
					}
					if(pInFormat==NULL)
					{
						cerr << program_name << ": cannot read input format!" << endl;
						usage();
					}
                                inputExt = iext;
					break;
					
				case 'o':
					gotOutType = true;
					oext = argv[arg] + 2;
					if(!*oext)
						oext = argv[++arg]; //space left after -i: use next argument
					
					if (strncasecmp(oext, "MIME", 4) == 0)
					{
						// get the MIME type from the next argument
						oext = argv[++arg];
						pOutFormat = Conv.FormatFromMIME(oext);
					}
					else
						pOutFormat = Conv.FindFormat(oext);

					if(pOutFormat==NULL)
					{
						cerr << program_name << ": cannot write output format!" << endl;
						usage();
					}
                                outputExt = oext;
					break;				
				
				case '?':
				case 'H':
					if(isalnum(argv[arg][2]))
					{
						if(strncasecmp(argv[arg]+2,"all",3))
						{
							OBFormat* pFormat = Conv.FindFormat(argv[arg]+2);
							if(pFormat)
								{
									cout << argv[arg]+2 << "  " << pFormat->Description() << endl;
									if(strlen(pFormat->SpecificationURL()))
							cout << "Specification at: " << pFormat->SpecificationURL() << endl;
								}
							else
								cout << "Format type: " << argv[arg]+2 << " was not recognized" <<endl;
						}
						else
						{
							Formatpos pos;
							OBFormat* pFormat;
							const char* str=NULL;
							while(OBConversion::GetNextFormat(pos,str,pFormat))
								{
									if((pFormat->Flags() & NOTWRITABLE) && (pFormat->Flags() & NOTREADABLE))
										continue;
									cout << str << endl;
									const char* p = strchr(pFormat->Description(),'\n');
									cout << p+1; //second line of description
									if(strlen(pFormat->SpecificationURL()))
							cout << "Specification at: " << pFormat->SpecificationURL();
									cout << endl << endl;
								}
						}
					}
					else
						help();
						exit(0);					
					
					
				default: //single character general option
					p = argv[arg]+1;
					DoOption(p,Conv,OBConversion::GENOPTIONS,arg,argc,argv);
					break;
				}
			}
			else 
			{
				//filenames
				if(!gotOutType)
					FileList.push_back(argv[arg]);
				else
					OutputFileName = argv[arg];
			}
		}
	}
    
  //user didn't specify input and output format in commandline option
  //try to parse it from program name (pdb2mdin means input format is pdb, output format is mdin)
  string formatName(program_name);
  pos = formatName.find_first_of("2");
  if(pos!=string::npos) {
    if (!gotInType)
    {
        string tmpExt = formatName.substr(0, pos);
        pInFormat = Conv.FindFormat(tmpExt.c_str());
 	 if(pInFormat==NULL)
 	 {
 		cerr << program_name << ": cannot read input format!" << endl;
 		usage();
 	 } else 
        {
            gotInType = true;
            inputExt = tmpExt;
 	 }
    }

    if (!gotOutType)
    {
        string tmpExt = formatName.substr(pos+1, string::npos);
        pOutFormat = Conv.FindFormat(tmpExt.c_str());
   	 if(pOutFormat==NULL)
 	 {
 		cerr << program_name << ": cannot write output format!" << endl;
 		usage();
 	 }else {
            gotOutType = true;
            outputExt = tmpExt;
 	 }
    }
  }

    if(FileList.empty())
    {
    	cerr << "No input file or format spec!" <<endl;
    	usage();
    }

   if (OutputFileName.empty())
   {
        pos = FileList.front().rfind(".");
        if(pos==string::npos)
            OutputFileName = FileList.front()+ "." + outputExt;
        else
            OutputFileName = FileList.front().substr(0, pos) + "." + outputExt;            
    }
  
  Conv.SetInAndOutFormats(pInFormat,pOutFormat);


  // send info message to clog -- don't mess up cerr or cout for user programs
  int count = Conv.FullConvert(FileList, OutputFileName, OutputFileList);
  if ( count == 1 )
    clog << count << " molecule converted" << endl;
  else
    clog << count << " molecules converted" << endl;
  
  if(OutputFileList.size()>1)
  {
    clog << OutputFileList.size() << " files output. The first is " << OutputFileList[0] <<endl;
  }

#ifdef _DEBUG
  //CM keep window open
  cout << "Press any key to finish" <<endl;
  getch();
#endif
  
  return 0;
}

void DoOption(const char* p, OBConversion& Conv, OBConversion::Option_type typ, int& arg, int argc, char *argv[]) 
{
	while(p && *p) //can have multiple single char options
	{
		char ch[2]="?";
		*ch = *p++;
		const char* txt=NULL;				
		//Get the option text if needed
		int nParams = Conv.GetOptionParams(ch, typ);
		if(nParams)
		{
			if(*p)
			{
				txt = p; //use text immediately following the option letter
				p=NULL; //no more single char options
			}
			else if(arg<argc-1)
			{
				txt = argv[++arg]; //use text from next arg
				if(*txt=='-')
				{
					//...unless it is another option
					cerr << "Option -" << ch << " takes a parameter" << endl;
					exit(0);
				}
			}
		}
		Conv.AddOption(ch, typ, txt);
	}
}

void usage()
{
  cout << "Open Babel " << BABEL_VERSION << " -- " << __DATE__ << " -- "
       << __TIME__ << endl;
  cout << "Usage: " << program_name
       << " [-i<input-type>] <name> [-o<output-type>] <name>" << endl;
  cout << "Try  -H option for more information." << endl;
/*  
#ifdef _DEBUG
  //CM keep window open
  cout << "Press any key to finish" <<endl;
  getch();
#endif
*/  
  exit (0);
}

void help()
{
  cout << "Open Babel converts chemical structures from one file format to another"<< endl << endl;
  cout << "Usage: " << program_name << " <input spec> <output spec> [Options]" << endl << endl;
  cout << "Each spec can be a file whose extension decides the format." << endl;
  cout << "Optionally the format can be specified by preceding the file by" << endl;
  cout << "-i<format-type> e.g. -ipdb, for input and -o<format-type> for output" << endl << endl;
  cout << "See below for available format-types, which are the same as the " << endl;
  cout << "file extensions and are case independent." << endl; 
  cout << "If no input or output file is given stdin or stdout are used instead." << endl << endl; 
  cout << "More than one input file can be specified and their names can contain" <<endl;
  cout << "wildcard chars (* and ?).The molecules are aggregated in the output file.\n" << endl;
  cout << OBConversion::Description(); // Conversion options
  cout << "  -H Outputs this help text" << endl;
  cout << "  -Hxxx (xxx is file format ID e.g. -Hpdb) gives format info" <<endl; 
  cout << "  -Hall Outputs details of all formats" <<endl; 
  cout << "  -V Outputs version number" <<endl; 

  
  OBFormat* pDefault = OBConversion::GetDefaultFormat();
	if(pDefault)
		cout << pDefault->TargetClassDescription();// some more options probably for OBMol
  
	OBFormat* pAPI= OBConversion::FindFormat("obapi");
	if(pAPI)
		cout << pAPI->Description();

	cout << "The following file formats are recognized:" << endl;
  Formatpos pos;
  OBFormat* pFormat;
  const char* str=NULL;
  while(OBConversion::GetNextFormat(pos,str,pFormat))
  {
		if((pFormat->Flags() & NOTWRITABLE) && (pFormat->Flags() & NOTREADABLE))
			continue;
		cout << "  " << str << endl;
	}
  cout << "\nSee further specific info and options using -H<format-type>, e.g. -Hpdb" << endl;
}

