#ifndef CONSTRAINT_ZCONSDATA_HPP
#define CONSTRAINT_ZCONSDATA_HPP
#include <algorithm>

#include "utils/GenericData.hpp"

#define ZCONSTIME_ID            "ZCONSTIME"
#define ZCONSPARADATA_ID    "ZCONSPARA"
#define ZCONSFILENAME_ID     "ZCONSFILENAME"
#define ZCONSTOL_ID              "ZCONSTOL"
#define ZCONSFORCEPOLICY_ID "ZCONSFORCEPOLICY"
#define ZCONSGAP_ID "ZCONSGAP"
#define ZCONSFIXTIME_ID "ZCONSFIXTIME"
#define ZCONSUSINGSMD_ID "ZCONSUSINGSMD"



using namespace std;
using namespace oopse;

struct ZConsParaItem {
    int zconsIndex;
    bool havingZPos;
    double zPos;
    double kRatio;
    bool havingCantVel;
    double cantVel;
};


class ZConsParaSortCriterion{
    public:
        bool operator ()(const ZConsParaItem& item1, const ZConsParaItem& item2){
            return item1.zconsIndex < item2.zconsIndex;
        }

};

class ZConsParaData : public GenericData{

    public:
        ZConsParaData(const string& id = ZCONSPARADATA_ID) : GenericData(id) {}
            
        void addItem(ZConsParaItem& item) {data.push_back(item);}

        vector<ZConsParaItem>* getData() {return &data;}

        void setData(vector<ZConsParaItem>& theData) {data = theData;}

        void sortByIndex() {
            sort(data.begin(), data.end(), ZConsParaSortCriterion());
        }

        bool isIndexUnique() {

            for(int i = 0; i < (int)(data.size() - 1); i++) {
                for(int j = i + 1; j < (int)(data.size()); j++) {
                    if(data[i].zconsIndex == data[j].zconsIndex) {
                        return false;  
                    }
                }
            }
            
            return true;
        }

    private:
        vector<ZConsParaItem> data;
};

#endif //CONSTRAINT_ZCONSDATA_HPP
