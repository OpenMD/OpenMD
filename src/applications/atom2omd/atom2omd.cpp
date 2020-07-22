/**********************************************************************
atom2md.cpp - OpenBabel-based conversion program to OpenMD file, 
              command-line handling.

Copyright (C) 1998-2001 by OpenEye Scientific Software, Inc.
Some portions Copyright (C) 2001-2006 by Geoffrey R. Hutchison
Some portions Copyright (C) 2004-2006 by Chris Morley
Some portions Copyright (C) 2004-2020 by J. Daniel Gezelter

This file is part of both the OpenMD and Open Babel projects.
For more information, see <http://openmd.net> and <http://openbabel.sourceforge.net/>

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation version 2 of the License.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.
***********************************************************************/

#include "config.h"

// used to set import/export for Cygwin DLLs
#ifdef WIN32
#define USING_OBDLL
#endif

#include <openbabel/babelconfig.h>

#include <iostream>
#include <fstream>
#include <sstream>

#include <string>
#include <map>
#if HAVE_CONIO_H
	#include <conio.h>
#endif
#include <cstdlib> // for exit() on Linux

#ifdef _MSC_VER
#define strncasecmp _strnicmp
#endif

#include <openbabel/obconversion.h>
#include <openbabel/plugin.h>

#include <cstring>

using namespace std;
using namespace OpenBabel;

void DoOption(const char* p, OBConversion& Conv, OBConversion::Option_type typ,
	      int& arg, int argc, char *argv[]); 
void usage();
void help();

// There isn't a great way to do this -- we need to save argv[0] for usage()
static char *program_name;

int main(int argc,char *argv[])
{
  OBConversion Conv(&cin, &cout); //default input and output are console 

  OBFormat* pInFormat = nullptr;
  OBFormat* pOutFormat = nullptr;
  bool outGzip = false;
  vector<string> FileList, OutputFileList;
  string OutputFileName;

  // Parse commandline
  bool gotInType = false, gotOutType = false;
  bool SplitOrBatch=false;

  char *oext = nullptr;
  char *iext = nullptr;

  // for use with command name to type conversion
  string inputExt;
  string outputExt;

  //Save name of program without its path (and .exe)
  string pn(argv[0]);
  string::size_type pos;
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
  for (arg = 1; arg < argc; ++arg)
    {
      if (argv[arg])
        {
          if (argv[arg][0] == '-')
            {
              char opchar[2]="?";
              opchar[0]=argv[arg][1];
              switch (opchar[0])
                {

                case 'V':
                  {
                    cout << program_name << ": part of OpenMD " << 
                      OPENMD_VERSION_MAJOR << "." << OPENMD_VERSION_MINOR << 
                      "." << OPENMD_VERSION_TINY << 
                      " and Open Babel " << BABEL_VERSION << " -- " 
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
                      exit(1);
                    }
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

                  if(pOutFormat == nullptr)
                    {
                      cerr << program_name << ": cannot write output format!" << endl;
                      usage();
                      exit(1);
                    }
                  break;
                  
                case 'O':
                  OutputFileName = argv[arg] + 2;
                  if(OutputFileName.empty())
                    OutputFileName = argv[++arg]; //space left after -O: use next argument
                  break;

                case 'L': //display a list of plugin type or classes
                  {
                    const char* param = nullptr;
                    if(argc>arg+1)
                      param = argv[arg+2];

                    // First assume first arg is a plugin type and
                    // param is a subtype, like babel -L ops gen3D
                    // or first arg is a plugin ID, like babel -L cml
                    OBPlugin* plugin;
                    if ((OBPlugin::GetPlugin("plugins", argv[arg+1]) &&
                         (plugin = OBPlugin::GetPlugin(argv[arg+1], param))) ||
                        (plugin = OBPlugin::GetPlugin(nullptr, argv[arg+1])))
                    {
                      //Output details of subtype
                      string txt;
                      plugin->Display(txt, "verbose", argv[arg+1]);
                      cout << "One of the " << plugin->TypeID() << '\n' << txt << endl;
                      return 0;
                    }
                    //...otherwise assume it is a plugin type, like babel -L forcefields
                    //Output list of subtypes
                    OBPlugin::List(argv[arg+1], param);
                    return 0;                  }
                case '?':
                case 'H':
                  if(isalnum(argv[arg][2]) || arg==argc-2)
                    {
                      if(strncasecmp(argv[arg]+2,"all",3))
                        {
                          OBFormat* pFormat
                            = (arg==argc-2) ? Conv.FindFormat(argv[arg+1]) : Conv.FindFormat(argv[arg]+2);
                          if(pFormat)
                            {
                              cout << argv[arg]+2 << "  " << pFormat->Description() << endl;
                              if(pFormat->Flags() & NOTWRITABLE)
                                cout << " This format is Read-only" << endl;
                              if(pFormat->Flags() & NOTREADABLE)
                                cout << " This format is Write-only" << endl;

                              if(strlen(pFormat->SpecificationURL()))
                                cout << "Specification at: " << pFormat->SpecificationURL() << endl;
                            }
                          else
                            cout << "Format type: " << argv[arg]+2 << " was not recognized" <<endl;
                        }
                      else
                        {
                          OBPlugin::List("formats","verbose");
                        }
                    }
                  else
                    help();
                  return 0;
                  
                case '-': //long option --name text
                  {
                    //Option's text is in the next and subsequent args, until one starts with -
                    char* nam = argv[arg]+2;
                    if(!strcasecmp(nam, "help")) //special case handled here
                    {
                      help();
                      return 0;
                    }
                    if(*nam != '\0') //Do nothing if name is empty
                      {
                        string txt;
                        while(arg<argc-1 && *argv[arg+1]!='-')
                          {
                            //use text from subsequent args
                            if(!txt.empty())txt += ' '; //..space separated if more than one
                            txt += argv[++arg]; 
                          }

                        // If a API directive, e.g.---errorlevel
                        // send to the pseudoformat "obapi" (without any leading -)
                        if(*nam=='-')
                          {
                            OBConversion apiConv;
                            OBFormat* pAPI= OBConversion::FindFormat("obapi");
                            if(pAPI)
                              {
                                apiConv.SetOutFormat(pAPI);
                                apiConv.AddOption(nam+1, OBConversion::GENOPTIONS, txt.c_str());
                                apiConv.Write(nullptr, &std::cout);
                              }
                          }
                        else
                          // Is a normal long option name, e.g --addtotitle
                          Conv.AddOption(nam,OBConversion::GENOPTIONS,txt.c_str());
                      }
                  }
                  break;
                  
                case 'm': //multiple output files
                  SplitOrBatch=true;
                  break;
					
                case 'a': //single character input option
                  p = argv[arg]+2;
                  DoOption(p,Conv,OBConversion::INOPTIONS,arg,argc,argv);
                  break;

                case 'x': //single character output option
                  p = argv[arg]+2;
                  DoOption(p,Conv,OBConversion::OUTOPTIONS,arg,argc,argv);
                  break;
					
                //Not essential, but allows these options to be before input filenames
                //since we know they take one parameter, and are the most likely options to be misplaced
                case 'f':
                case 'l':
                  p = argv[arg] + 2;
                  if(!*p)
                    p = argv[++arg]; //space left after -f: use next argument
                  Conv.AddOption(opchar, OBConversion::GENOPTIONS, p);
                  break;
                
                case ':':
                  //e.g. -:c1ccccc1. SMILES passed as a file name and handled in OBConversion
                  FileList.push_back(argv[arg]);
                  break;

                default: //single character general option
                  p = argv[arg]+1;
                  DoOption(p,Conv,OBConversion::GENOPTIONS,arg,argc,argv);
                  break;
                }
            }
          else //filenames
            FileList.push_back(argv[arg]);
        }
    }
  
  // user didn't specify input and output format in commandline option
  // try to parse it from program name (pdb2omd means input format is pdb, 
  // output format is omd)

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

  if(!gotOutType) //the last file is the output
    {
      if(FileList.empty())
        {
          cerr << "No output file or format spec!" << endl;
          usage();
        }
      OutputFileName = FileList.back();
      FileList.pop_back();
    }

#if defined(_WIN32) && defined(USING_DYNAMIC_LIBS)
  //Expand wildcards in input filenames and add to FileList
  vector<string> tempFileList(FileList);
  FileList.clear();
  vector<string>::iterator itr;
  for(itr=tempFileList.begin();itr!=tempFileList.end();++itr)
  {
    if((*itr)[0]=='-')
      FileList.push_back(*itr);
    else
      DLHandler::findFiles (FileList, *itr);
  }
#endif
  
  if (!gotInType)
    {
      if(FileList.empty())
        {
          cerr << "No input file or format spec or possibly a misplaced option.\n"
            "Most options must come after the input files. (-i -o -O -m can be anywhwere.)\n" <<endl;
          usage();
          exit(1);
        }
    }

  if (!gotOutType)
    {
      //check there is a valid output format, but the extension will be re-interpreted in OBConversion
      pOutFormat = Conv.FormatFromExt(OutputFileName.c_str(), outGzip);
      if (OutputFileName.empty() || pOutFormat==nullptr)
        {
          cerr << "Missing or unknown output file or format spec or possibly a misplaced option.\n"
            "Options, other than -i -o -O -m, must come after the input files.\n" <<endl;
          usage();
          exit(1);
        }
    }
  
    if(!Conv.SetInFormat(pInFormat)) //rely on autodetection for gzipped input
    {
      cerr << "Invalid input format" << endl;
      usage();
      exit(1);
    }
    if(!Conv.SetOutFormat(pOutFormat, outGzip))
    {
      cerr << "Invalid output format" << endl;
      usage();
      exit(1);
    }

  if(SplitOrBatch)
    {
      //Put * into output file name before extension (or ext.gz)
      if(OutputFileName.empty())
        {
          OutputFileName = "*.";
          OutputFileName += oext;
        }
      else
        {
          string::size_type pos = OutputFileName.rfind(".gz");
          if(pos==string::npos)
            pos = OutputFileName.rfind('.');
          else
            pos = OutputFileName.rfind('.',pos-1);
          if(pos==string::npos)
            OutputFileName += '*';
          else
            OutputFileName.insert(pos,"*");
        }
    }
  
  int count = Conv.FullConvert(FileList, OutputFileName, OutputFileList);
 
  Conv.ReportNumberConverted(count);

  if(OutputFileList.size()>1)
    {
      clog << OutputFileList.size() << " files output. The first is " << OutputFileList[0] <<endl;
    }

  std::string messageSummary = obErrorLog.GetMessageSummary();
  if (messageSummary.size())
    {
      clog << messageSummary << endl;
    }

#ifdef DEBUG
  //CM keep window open
  cout << "Press any key to finish" <<endl;
  getch();
#endif
  
  return 0;
}

void DoOption(const char* p, OBConversion& Conv,
	      OBConversion::Option_type typ, int& arg, int argc, char *argv[]) 
{
  //Unlike babel, cannot have multiple concatenated single char options
  //accepts: -sCCC -s CCC -s"CCC" -s CCC red -sCCC red
  char ch[2]="?";
  *ch = *p++;
  std::string txt;
  //Get the option text
  if(*p)
    txt = p; //use text immediately following the option letter, and keep looking

  while(arg<argc-1 && *argv[arg+1]!='-')
  {
    //use text from subsequent args
    if(!txt.empty())txt += ' '; //..space separated if more than one
    txt += argv[++arg]; 
  }
  Conv.AddOption(ch, typ, txt.c_str());
}

void usage()
{
  cout << program_name << ": part of OpenMD " << 
    OPENMD_VERSION_MAJOR << "." << OPENMD_VERSION_MINOR << "." << 
    OPENMD_VERSION_TINY << " and OpenBabel " << BABEL_VERSION << " -- " 
       << __DATE__ << " -- " << __TIME__ << endl;
  cout << "Usage: " << program_name
       << " [-i<input-type>] <name> [-o<output-type>] <name>" << endl;
  cout << "Try  -H option for more information." << endl;
 
#ifdef DEBUG
  //CM keep window open
  cout << "Press any key to finish" <<endl;
  getch();
#endif
  exit (0);
}

void help()
{
  cout << "Open Babel converts chemical structures from one file format to another"<< endl << endl;
  cout << "Usage: " << endl;
  cout << program_name << "[-i<input-type>] <infilename> [-o<output-type>] -O<outfilename> [Options]" << endl;
  cout << "The extension of a file decides the format, unless it is overridden" <<endl;
  cout << " by -i or -o options, e.g. -icml, or -o smi" << endl;
  cout << "See below for available format-types, which are the same as the " << endl;
  cout << "file extensions and are case independent." << endl; 
  cout << "If no input or output file is given stdin or stdout are used instead." << endl << endl; 
  cout << "More than one input file can be specified and their names can contain" <<endl;
  cout << "wildcard chars (* and ?). The format of each file can be different unless" <<endl;
  cout << "the -i option has been used, when they are all the same." <<endl;
  cout << "By default, the molecules are aggregated in the output file," << endl;
  cout << " but see -m option, Splitting, below.\n" << endl;

  cout << "Options, other than -i -o -O -m, must come after the input files.\n" <<endl;
  cout << OBConversion::Description(); // Conversion options
  cout << "-H Outputs this help text" << endl;
  cout << "-Hxxx (xxx is file format ID e.g. -Hcml) gives format info" <<endl; 
  cout << "-Hall Outputs details of all formats" <<endl; 
  cout << "-V Outputs version number" <<endl; 
  cout << "-L <category> Lists plugin classes of this category, e.g. <formats>" << endl;
  cout << "   Use just -L for a list of plugin categories." << endl; 
  cout << "   Use -L <ID> e.g. -L sdf for details of a format or other plugin." << endl; 
  cout << "-m Produces multiple output files, to allow:" <<endl;
  cout << "    Splitting: e.g.        " << program_name << " infile.mol -O new.smi -m" <<endl;
  cout << "      puts each molecule into new1.smi new2.smi etc" <<endl;
  cout << "    Batch conversion: e.g. " << program_name << " *.mol -osmi -m" <<endl;
  cout << "      converts each input file to a .smi file" << endl;
#ifdef _WIN32
  cout << "   In Windows these can also be done using the forms" <<endl;
  cout << "     " << program_name << " infile.mol -O new*.smi and " << program_name << " *.mol -O *.smi respectively.\n" <<endl;
#endif
  
  OBFormat* pDefault = OBConversion::GetDefaultFormat();
  if(pDefault)
    cout << pDefault->TargetClassDescription();// some more options probably for OBMol
  
  OBFormat* pAPI= OBConversion::FindFormat("obapi");
  if(pAPI)
    cout << pAPI->Description();
  
  cout << "To see a list of recognized file formats use\n  babel -L formats [read] [write]\n"
       << "To see details and specific options for a particular format, e.g CML, use\n  babel -L cml\n"
       << endl;
  //cout << "The following file formats are recognized:" << endl;
  //OBPlugin::List("formats");
  //cout << "\nSee further specific info and options using -H<format-type>, e.g. -Hcml" << endl;
}

