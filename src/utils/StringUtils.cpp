#include "utils/StringUtils.hpp"

using namespace std;
using namespace oopse;

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

int findBegin(istream theStream, char* startText ){
  const int MAXLEN = 1024;
  char readLine[MAXLEN];   
  int foundText = 0;
  int lineNum;
  char* the_token;
  char* eof_test;

  // rewind the stream
  theStream.seekg (0, ios::beg);
  lineNum = 0;
  
  if (!theStream.eof()) {
    theStream.getline(readLine, MAXLEN);
    lineNum++;
  } else {
    printf( "Error fast forwarding stream: stream is empty.\n");
    return -1;
  }
  
  while ( !foundText ) {
    
    if (theStream.eof()) {
      printf("Error fast forwarding stream at line %d: "
             "stream ended unexpectedly.\n", lineNum);
      return -1;
    }

    the_token = strtok( readLine, " ,;\t" );

    if (!strcasecmp("begin", the_token)) {
      the_token = TrimSpaces(strtok( NULL, " ,;\t" ));
      foundText = !strcasecmp( startText, the_token );
    }

    if (!foundText) {
      if (!theStream.eof()) {
        theStream.getline(readLine, MAXLEN);
        lineNum++;
      } else {
        printf( "Error fast forwarding stream at line %d: "
                "stream ended unexpectedly.\n", lineNum);
        return -1;
      }
    }
  }
  return lineNum;
}

int countTokens(char *line, char *delimiters) {
  /* PURPOSE: RETURN A COUNT OF THE NUMBER OF TOKENS ON THE LINE. */
  
  char *working_line;   /* WORKING COPY OF LINE. */
  int ntokens;          /* NUMBER OF TOKENS FOUND IN LINE. */
  char *strtok_ptr;     /* POINTER FOR STRTOK. */

  strtok_ptr= working_line= strdup(line);

  ntokens=0;
  while (strtok(strtok_ptr,delimiters)!=NULL)
    {
      ntokens++;
      strtok_ptr=NULL;
    }

  free(working_line);
  return(ntokens);
}

char* TrimSpaces(char *str) {
  size_t len;
  char *right, *left;
  
  if (strlen(str) == 0) return(str);
  
  /* Trim whitespace from left side */
  for (left = str; isspace(*left); left++);
  
  /* Trim whitespace from right side */
  if ((len = strlen(left)))
    {
      right = left + (len - 1);
      
      while (isspace(*right))
        {
          *right = '\0';
          right--;
        }
    }
  
  /* Only do the str copy if their was spaces to the left */
  if (left != str)
    strcpy(str, left);
  
  return (str);
}

int isEndLine(char *line) {
  char *working_line;
  char *foo;
  
  working_line = strdup(line);
  
  foo = strtok(working_line, " ,;\t");

 if (!strcasecmp(foo, "end")) return 1;
 
 return 0;
}
