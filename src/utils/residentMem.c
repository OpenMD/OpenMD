#include <config.h> 
#include <string.h> 
#include <stdio.h>

double residentMem () {

  FILE* procresults;
  char buf[150];
  char* foo;
  long int myRSS, totRSS;
  char* pscommand;

  pscommand = strdup("PS");

#if PSTYPE == BSD
  strcat(pscommand, " ax -o rss");
#else
  strcat(pscommand, " -ef -o rss");
#endif

  procresults = popen(pscommand, "r");

  totRSS = 0;
  while (!feof(procresults)) {
    fgets( buf, 150, procresults);
    if (!strcmp(buf, " RSS")) continue;
        
    foo = strtok(buf, " ,;\t");
    myRSS = atoi(foo);
    totRSS += myRSS;
  } 
  pclose(procresults);

  return(totRSS);

}
