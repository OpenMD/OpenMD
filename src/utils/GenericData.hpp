#ifndef __GENERICDATA_H__
#define __GENERICDATA_H__

#include <algorithm>
#include <string>
#include <vector>

#define ZCONSTIME_ID            "ZCONSTIME"
#define ZCONSPARADATA_ID    "ZCONSPARA"
#define ZCONSFILENAME_ID     "ZCONSFILENAME"
#define ZCONSTOL_ID              "ZCONSTOL"
#define ZCONSFORCEPOLICY_ID "ZCONSFORCEPOLICY"
#define ZCONSGAP_ID "ZCONSGAP"
#define ZCONSFIXTIME_ID "ZCONSFIXTIME"
#define ZCONSUSINGSMD_ID "ZCONSUSINGSMD"

#define CHIVALUE_ID "CHIVALUE"
#define INTEGRALOFCHIDT_ID "INTEGRALOFCHIDT"
#define ETAVALUE_ID "ETAVALUE"

using namespace std;
////////////////////////////////////////////////////////////////////////////////
//Declaration of GenericData
////////////////////////////////////////////////////////////////////////////////
class GenericData
{
  public:
    GenericData();
    GenericData(const GenericData& rhs)  {  id = rhs.getID(); }
    GenericData& operator =(const GenericData& rhs);
    virtual ~GenericData() {}

    const string& getID() const {  return id; }
    void setID(const string& rhs)     {  id = rhs;  }

  protected:
    string id;
};

/** 
 * Something we can improve it here is to use template
 */
////////////////////////////////////////////////////////////////////////////////
//Declaration of IntData
////////////////////////////////////////////////////////////////////////////////
class IntData : public GenericData{

  public:

    double getData()         { return data; }
    void setData(int rhs) { data = rhs;  }

  protected:
    int data;
};

////////////////////////////////////////////////////////////////////////////////
//Declaration of DoubleData
////////////////////////////////////////////////////////////////////////////////
class DoubleData : public GenericData{

  public:

    double getData()         { return data; }
    void setData(double rhs) { data = rhs;  }

  protected:
    double data;
};

////////////////////////////////////////////////////////////////////////////////
//Declaration of StringData
////////////////////////////////////////////////////////////////////////////////
class StringData : public GenericData{

  public:
    const string& getData() const  {  return data; }
    void setData(const string& rhs) {  data = rhs;  }
  protected:
    string data;
};

////////////////////////////////////////////////////////////////////////////////
//Declaration of BoolData
////////////////////////////////////////////////////////////////////////////////
class BoolData : public GenericData{
  public:
    bool getData() const  {  return data; }
    void setData(const bool rhs) {  data = rhs;  }
  protected:
    bool data;

};

////////////////////////////////////////////////////////////////////////////////
//Declaration of ZConsParaData
////////////////////////////////////////////////////////////////////////////////
struct ZConsParaItem {
  int zconsIndex;
  bool havingZPos;
  double zPos;
  double kRatio;
  bool havingCantVel;
  double cantVel;
};

class ZConsParaData : public GenericData{

  public:
    ZConsParaData();
    void addItem(ZConsParaItem& item) {data.push_back(item);}
    vector<ZConsParaItem>* getData() {return &data;}
    void setData(vector<ZConsParaItem>& theData) {data = theData;}
    void sortByIndex();
    bool isIndexUnique();

  private:
    vector<ZConsParaItem> data;
  };

class ZConsParaSortCriterion{
  public:
    bool operator ()(const ZConsParaItem& item1, const ZConsParaItem& item2){
      return item1.zconsIndex < item2.zconsIndex;
    }

};

////////////////////////////////////////////////////////////////////////////////
//Declaration of IntData
////////////////////////////////////////////////////////////////////////////////
class DoubleArrayData : public GenericData{

 public:
   vector<double> getData() const  {  return data; }
   void setData(double* source, int num){
    data.clear();
    for(int i = 0; i < num; i++)
     data.push_back(source[i]);
   }
 protected:
   vector<double> data;
};

////////////////////////////////////////////////////////////////////////////////
//Declaration of AtomData
////////////////////////////////////////////////////////////////////////////////
struct AtomInfo : public GenericData{
  public:
    string AtomType;
    double pos[3];
    double dipole[3];  
};

class AtomData : public GenericData{
  public:
    ~AtomData();
    void addAtomInfo(AtomInfo* info) {data.push_back(info);}
    void clearAllAtomInfo();
    AtomInfo* beginAtomInfo(vector<AtomInfo*>::iterator& i){
      i = data.begin();
      return i != data.end()? *i : NULL;
    }
    AtomInfo* nextAtomInfo(vector<AtomInfo*>::iterator& i){
      ++i;
      return i != data.end()? *i: NULL;
    }
    vector<AtomInfo*> getData() {return data;}
    int getSize() {return data.size();}
  protected:
    vector<AtomInfo*> data;
};

#endif
