#ifndef __STRINGUTILS_H__
#define __STRINGUTILS_H__

#ifdef __cplusplus
extern "C" {
#endif
  
  extern int count_tokens(char *line, char *delimiters);
  int isEndLine(char *line);
  int isBeginLine(char *line);
  char* TrimSpaces(char *str);

#ifdef __cplusplus
}
#endif

#endif
