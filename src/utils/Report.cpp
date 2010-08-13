#include <protomol/base/Report.h>

#include <iomanip>

#include <protomol/base/SystemUtilities.h>

using namespace std;

namespace ProtoMol {
//____ Report
  namespace Report {
    MyStreamer &debug::operator()(MyStreamer &stream) const {
      stream.setLevel(myLevel);
      stream << "DEBUG:  ";
      return stream;
    }

    MyStreamer &reportlevel::operator()(MyStreamer &stream) const {
      stream.setReportLevel(myReportlevel);
      return stream;
    }
    //___________________________________________________________ 
    // Currently defaulting to std::cerr
    MyStreamer report(&(std::cerr));

    // MyStreamer
    MyStreamer::MyStreamer(ostream *a) {
      setStream(a);
    }

    // Allows the changing of the stream
    ostream *MyStreamer::setStream(ostream *a) {
      ostream *tmp = myStream;
      myStream = (a != NULL) ? a : &(cerr);
      myResetFlags = myStream->flags();
      myAbort = false;
      myQuit = false;
      myAllNodes = false;
      myAllNodesSerial = false;
      myAllSlavesSerial = false;
      myIAmMaster = true;
      myDoHint = true;
      mySilentHint = 0;
      myReportLevel = 0;
      myLevel = 0;
      return tmp;
    }

    void MyStreamer::setAbort(bool flag) {
      myAbort = flag;
    }

    void MyStreamer::setQuit(bool flag) {
      myQuit = flag;
    }

    void MyStreamer::setAllNodes(bool flag) {
      myAllNodes = flag;
    }

    void MyStreamer::setf(ios::fmtflags flag) {
      myStream->setf(flag);
    }

    void MyStreamer::setf(ios::fmtflags flag, ios::fmtflags mask) {
      myStream->setf(flag, mask);
    }

    void MyStreamer::precision(int prec) {
      myStream->precision(prec);
    }

    void MyStreamer::width(int wide) {
      myStream->width(wide);
    }

    void MyStreamer::reset() {
      myStream->flags(myResetFlags);
    }

    // MyStreamer
    // Standard Overloads

    MyStreamer &MyStreamer::operator<<(bool a) {
      if (print()){
        if (a) {
		  *myStream << "true";
		}else{
		  *myStream << "false";
		}
	  }
      return *this;
    }

    MyStreamer &MyStreamer::operator<<(char a) {
      if (print())
        *myStream << a;
      return *this;
    }

    MyStreamer &MyStreamer::operator<<(unsigned char a) {
      if (print())
        *myStream << a;
      return *this;
    }

    MyStreamer &MyStreamer::operator<<(signed char a) {
      if (print())
        *myStream << a;
      return *this;
    }

    MyStreamer &MyStreamer::operator<<(char *a) {
      if (print())
        *myStream << a;
      return *this;
    }

    MyStreamer &MyStreamer::operator<<(const char *a) {
      if (print())
        *myStream << a;
      return *this;
    }

    MyStreamer &MyStreamer::operator<<(unsigned char *a) {
      if (print())
        *myStream << a;
      return *this;
    }

    MyStreamer &MyStreamer::operator<<(signed char *a) {
      if (print())
        *myStream << a;
      return *this;
    }

    MyStreamer &MyStreamer::operator<<(short a) {
      if (print())
        *myStream << a;
      return *this;
    }

    MyStreamer &MyStreamer::operator<<(unsigned short a) {
      if (print())
        *myStream << a;
      return *this;
    }

    MyStreamer &MyStreamer::operator<<(int a) {
      if (print())
        *myStream << a;
      return *this;
    }

    MyStreamer &MyStreamer::operator<<(unsigned int a) {
      if (print())
        *myStream << a;
      return *this;
    }

    MyStreamer &MyStreamer::operator<<(long a) {
      if (print())
        *myStream << a;
      return *this;
    }

    MyStreamer &MyStreamer::operator<<(unsigned long a) {
      if (print())
        *myStream << a;
      return *this;
    }

    MyStreamer &MyStreamer::operator<<(long long a) {
      if (print())
#if (defined (_GLIBCPP_USE_LONG_LONG) || !defined (__GNUC__))
        *myStream << a;
#else
        *myStream << (long)a;
#endif
      return *this;
    }

    MyStreamer &MyStreamer::operator<<(unsigned long long a) {
      if (print())
#if (defined (_GLIBCPP_USE_LONG_LONG) || !defined (__GNUC__))
        *myStream << a;
#else
        *myStream << (unsigned long)a;
#endif
      return *this;
    }

    MyStreamer &MyStreamer::operator<<(float a) {
      if (print())
        *myStream << a;
      return *this;
    }

    MyStreamer &MyStreamer::operator<<(double a) {
      if (print())
        *myStream << a;
      return *this;
    }

    MyStreamer &MyStreamer::operator<<(const ostream &a) {
      if (print())
        *myStream << a;
      return *this;
    }

    MyStreamer &MyStreamer::operator<<(ostream *a) {
      if (print())
        *myStream << a;
      return *this;
    }

    MyStreamer &MyStreamer::operator<<(ios & (*f)(ios &)) {
      if (print())
        (*f)(*myStream);
      return *this;
    }

    MyStreamer &MyStreamer::operator<<(ostream & (*f)(ostream &)) {
      if (print())
        (*f)(*myStream);
      return *this;
    }

    MyStreamer &MyStreamer::operator<<(const string &a) {
      if (print())
        *myStream << a;
      return *this;
    }

    MyStreamer &MyStreamer::operator<<(MyStreamer & (*f)(MyStreamer &)) {
      (*f)(*this);
      return *this;
    }

    // Our Output Levels

    MyStreamer &plain(MyStreamer &stream) {
      stream.myLevel = -4;
      return stream;
    }

    MyStreamer &error(MyStreamer &stream) {
      stream.myAbort = true;
      stream.myLevel = -3;
      if (stream.print())
        stream << "Fatal Error:  ";
      return stream;
    }

    MyStreamer &recoverable(MyStreamer &stream) {
      stream.myLevel = -2;
      if (stream.print())
        stream << "Recoverable Error:  ";
      return stream;
    }

    MyStreamer &warning(MyStreamer &stream) {
      stream.myLevel = -1;
      if (stream.print())
        stream << "Warning:  ";
      return stream;
    }

    MyStreamer &hint(MyStreamer &stream) {
      stream.myLevel = 0;
      if (!stream.myDoHint)
        stream.mySilentHint++;
      if (stream.print())
        stream << "Hint:  ";
      return stream;
    }

    MyStreamer &operator<<(MyStreamer &stream, const reportlevel &rl) {
      return rl(stream);
    }

    MyStreamer &operator<<(MyStreamer &stream, const debug &d) {      
      return d(stream);
    }

    MyStreamer &allnodes(MyStreamer &stream) {
      stream.myAllNodes = true;
      return stream;
    }

    MyStreamer &quit(MyStreamer &stream) {
      stream.myQuit = true;
      return stream;
    }

    MyStreamer &aborting(MyStreamer &stream) {
      stream.myAbort = true;
      return stream;
    }

    MyStreamer &dohint(MyStreamer &stream) {
      stream.myDoHint = true;
      return stream;
    }

    MyStreamer &donthint(MyStreamer &stream) {
      stream.myDoHint = false;
      return stream;
    }

    // Our Output Controll
    MyStreamer &allnodesserial(MyStreamer &stream) {
      stream.myAllNodes = true;
      stream.myAllNodesSerial = true;
      stream << flush;
      protomolStartSerial(false);
      return stream;
    }

    MyStreamer &allslavesserial(MyStreamer &stream) {
      stream.myAllNodes = true;
      stream.myAllSlavesSerial = true;
      stream << flush;
      protomolStartSerial(true);
      return stream;
    }

    MyStreamer &endr(MyStreamer &stream) {
      stream << endl;
      if (!stream.myDoHint)
        stream.mySilentHint = max(stream.mySilentHint - 1, 0);
      if (stream.myQuit)
        if (!stream.myAllNodesSerial || !stream.myAllSlavesSerial)
          protomolExit();
        else
          protomolAbort();
      else if (stream.myAbort)
        protomolAbort();
      else if (stream.myAllNodesSerial)
        protomolEndSerial(false);
      else if (stream.myAllSlavesSerial)
        protomolEndSerial(true);
      stream.myAllSlavesSerial = false;
      stream.myAllNodesSerial = false;
      stream.myAllNodes = false;
      return stream;
    }
  }
}
