/*  -*- c++ -*-  */
#ifndef REPORT_H
#define REPORT_H

#include <iostream>
#include <string>
using std::cerr;

namespace ProtoMol {
  //________________________________________________________ Report
  namespace Report {
    class MyStreamer;

    //_________________________________________ debug
    /**
       Manipulator of MyStreamer, to print out debug
       information with different report level importance, where
       1, pretty important, 10 really unimportant.
     */
    class debug {
public:
      debug() : myLevel(1) {}
      debug(short l) : myLevel(l) {}
      MyStreamer &operator()(MyStreamer &stream) const;
private:
      short myLevel;
    };

    //_________________________________________ reportlevel
    /**
       reportlevel sets the maximal level of report, if less
       than 1 all debug information is suppressed.
     */
    class reportlevel {
public:
      reportlevel() : myReportlevel(0) {}
      reportlevel(short l) : myReportlevel(l) {}
      MyStreamer &operator()(MyStreamer &stream) const;
private:
      short myReportlevel;
    };

    //_________________________________________ MyStreamer
    /**
       MyStreamer wraps the output to a given stream and provides different
       manipulators a la std::cout and std::cerr. report is a global instance
       of MyStreamer and acts as std::cout and std::cerr. In parallel 
       environment, by default, only output from the master is passed to the
       actual stream. One can also pipe the putput to a file.
       @n
       It handles different report levels. Production code will by default
       suppress all debug output (print if level <= 0, reportlevel=0), where as
       debug compilation will by default print all information with level <= 1,
       reportlevel=1.
       @n
       NB! Do not put debug(<int>) inside forces, ok if it is inside
       initialization, otherwise use a type from oneAtomContraints.h to
       debug pair wise potentials @n @n

       report << debug(1) << "This is only debug inforamtion" << endr @n

       report << quit << plain << "Well, this storys gonna end soon ... " 
              << endr @n @n

       - plain              : -4
       - error              : -3
       - recoverable        : -2
       - warning            : -1
       - hint               :  0
       - debug(1)           :  1   important debug information
       - debug(2)           :  2   less important
       - debug(10)          : 10   really not important
       - aborting           : does an abort (MPI_Abort)
       - quit               : does an exit (MPI_Finalize)
       - endr               : end of report/MyStreamer output
       - allnodesserial     : output synchronized node by node including the
                              master
       - allslavesserial    : output synchronized node by node, only slaves
       - allnodes           : output unsynchronized from all nodes
       - donthint           : suppress all hint output
       - dohint             : enable hint output
     */
    class MyStreamer {
      //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      // Constructors, destructors, assignment
      //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
public:
      MyStreamer(std::ostream *);

      //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      // New methods of class MyStreamer
      //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      std::ostream *setStream(std::ostream *);
      void setAbort(bool);
      void setQuit(bool);
      void setAllNodes(bool);
      void setReportLevel(short level) {myReportLevel = level;}
      void setLevel(short level) {myLevel = level;}
      void setIAmMaster(bool master) {myIAmMaster = master;}

      //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      // New methods of class MyStreamer
      //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      void setf(std::ios::fmtflags flag);
      void setf(std::ios::fmtflags flag, std::ios::fmtflags mask);
      void precision(int prec);
      void width(int wide);
      void reset();

      //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      // New methods of class MyStreamer
      //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      MyStreamer &operator<<(bool);
      MyStreamer &operator<<(char);
      MyStreamer &operator<<(unsigned char);
      MyStreamer &operator<<(signed char);
      MyStreamer &operator<<(char *);
      MyStreamer &operator<<(const char *);
      MyStreamer &operator<<(unsigned char *);
      MyStreamer &operator<<(signed char *);
      MyStreamer &operator<<(short);
      MyStreamer &operator<<(unsigned short);
      MyStreamer &operator<<(int);
      MyStreamer &operator<<(unsigned int);
      MyStreamer &operator<<(long long);
      MyStreamer &operator<<(unsigned long long);
      MyStreamer &operator<<(long);
      MyStreamer &operator<<(unsigned long);
      MyStreamer &operator<<(float);
      MyStreamer &operator<<(double);
      MyStreamer &operator<<(const std::string &);
      MyStreamer &operator<<(const std::ostream &);
      MyStreamer &operator<<(std::ostream *);
      MyStreamer &operator<<(std::ostream &(*f)(std::ostream &));
      MyStreamer &operator<<(std::ios &(*f)(std::ios &));
      MyStreamer &operator<<(MyStreamer &(*f)(MyStreamer &));

      friend MyStreamer &allnodes(MyStreamer &stream);
      friend MyStreamer &plain(MyStreamer &stream);
      friend MyStreamer &hint(MyStreamer &stream);
      friend MyStreamer &quit(MyStreamer &stream);
      friend MyStreamer &aborting(MyStreamer &stream);
      friend MyStreamer &warning(MyStreamer &stream);
      friend MyStreamer &recoverable(MyStreamer &stream);
      friend MyStreamer &error(MyStreamer &stream);
      friend MyStreamer &allnodesserial(MyStreamer &stream);
      friend MyStreamer &allslavesserial(MyStreamer &stream);
      friend MyStreamer &endr(MyStreamer &stream);
      friend MyStreamer &dohint(MyStreamer &stream);
      friend MyStreamer &donthint(MyStreamer &stream);

private:
      bool print() const {
        return mySilentHint < 1 && myLevel <= myReportLevel &&
          (myAllNodes || myIAmMaster);
      }

      //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      // My data members
      //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
private:
      std::ostream *myStream;
      std::ios::fmtflags myResetFlags; // Resest/default flags
      bool myAbort;
      bool myQuit;
      bool myAllNodes;
      bool myAllNodesSerial;
      bool myAllSlavesSerial;
      bool myIAmMaster;
      bool myDoHint;
      short mySilentHint;
      short myReportLevel;
      short myLevel;
    };

    /// Plain output
    MyStreamer &plain(MyStreamer &stream);
    /// Error, abort
    MyStreamer &error(MyStreamer &stream);
    /// Recoverable error
    MyStreamer &recoverable(MyStreamer &stream);
    /// Warning
    MyStreamer &warning(MyStreamer &stream);
    /// Hint, can be suppressed by donthint/dohint
    MyStreamer &hint(MyStreamer &stream);

    MyStreamer &allnodes(MyStreamer &stream);
    MyStreamer &quit(MyStreamer &stream);
    MyStreamer &aborting(MyStreamer &stream);
    MyStreamer &allnodesserial(MyStreamer &stream);
    MyStreamer &allslavesserial(MyStreamer &stream);
    MyStreamer &endr(MyStreamer &stream);
    MyStreamer &dohint(MyStreamer &stream);
    MyStreamer &donthint(MyStreamer &stream);

    MyStreamer &operator<<(MyStreamer &stream, const reportlevel &rl);
    MyStreamer &operator<<(MyStreamer &stream, const debug &d);

    // Our global streamer
    extern MyStreamer report;
    //static MyStreamer report = &cerr;
  }
}
#endif /* REPORT_H */

