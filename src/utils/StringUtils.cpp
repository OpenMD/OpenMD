#include "StringUtils.hpp"

string UpperCase(const string& S) {
  string uc = S;
  unsigned int n = uc.size();
  for (unsigned int j = 0; j < n; j++) {   
    char sj = uc[j];
    if (sj >= 'a' && sj <= 'z') uc[j] = (char)(sj - ('a' - 'A'));
  }
  return uc;
}

string LowerCase(const string& S) {
  string lc = S;
  unsigned int n = lc.size();
  for (unsigned int j = 0; j < n; j++) {
    char sj = lc[j];
    if (sj >= 'A' && sj <= 'Z') lc[j] = (char)(sj + ('a' - 'A'));
  }
  return lc;
}
