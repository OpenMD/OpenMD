#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>
#include "simError.h"
#include "StringUtils.h"

int fastForwardToBeginBlock(FILE* theFile, char* blockName ){
  int foundText = 0;
  int lineNum;
  char* the_token;
  char* junk_token;
  char* eof_test;
  char readLine[500];
  
  rewind( theFile );
  lineNum = 0;
  
  eof_test = fgets( readLine, sizeof(readLine), theFile );
  lineNum++;
  if( eof_test == NULL ){
    sprintf( painCave.errMsg,
             "fastForwardToBeginBlock: File Is Empty!\n");
    painCave.isFatal = 1;
    simError();
    return -1;
  }

  while( !foundText ){

    if( eof_test == NULL ){
      sprintf( painCave.errMsg,
               "fastForwardToBeginBlock: File Ended Unexpectedly!\n"
               "\tAt Line: %d\n", lineNum);
      painCave.isFatal = 1;
      simError();
      return -1;
    }

    if (isBeginLine(readLine)) {
      junk_token = strtok( readLine, " ,;\t" );
      the_token  = TrimSpaces(strtok( NULL, " ,;\t" ));
      foundText = !strcasecmp( blockName, the_token );      
    }
    
    if( !foundText ){
      eof_test = fgets( readLine, sizeof(readLine), theFile );
      lineNum++;
      
      
      if (eof_test == NULL) {
        sprintf( painCave.errMsg,
                 "fastForwardToBeginBlock: File Ended Unexpectedly!\n"
                 "\tAt Line: %d\n", lineNum);
        painCave.isFatal = 1;
        simError();
        return -1;
      }
    }
  }

  return lineNum;
}

int isBeginLine(char *line) {
  /* checks to see if the first token on the line is "begin" */
  char *working_line;
  char *foo;
  
  working_line = strdup(line);
  
  foo = strtok(working_line, " ,;\t");
  
  if (!strcasecmp(foo, "begin")) return 1;
 
 return 0;
}

int isEndLine(char *line) {
  /* checks to see if the first token on the line is "end" */
  char *working_line;
  char *foo;
  
  working_line = strdup(line);
  
  foo = strtok(working_line, " ,;\t");
  
  if (!strcasecmp(foo, "end")) return 1;
 
 return 0;
}

/**
 * Removes left and right spaces from a string
 *
 * @param str  String to trim
 *
 * @return  char* to the trimed string
 */
char * TrimSpaces(char *str) {
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
