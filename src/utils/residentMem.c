#include <config.h> 
#include <string.h> 
#include <stdio.h>
#include <stdlib.h>
#ifdef __sgi
#include <unistd.h>
#endif

/* 
 * returns an estimate of the resident memory size in kB 
 */
RealType residentMem () {

  FILE* procresults;
  char buf[150];
  char* foo;
  long int myRSS, totRSS;
  char pscommand[150];
  char* psPath;

  STR_DEFINE(psPath, PSCOMMAND );
 
  
  strncpy(pscommand, psPath, strlen(psPath)+1);

#ifdef PSTYPE_IS_BSD
  strcat(pscommand, " ax -o rss");
#else
#ifdef PSTYPE_IS_POSIX
  strcat(pscommand, " -ef -o rss");
#else 
  printf("Unknown ps syntax!\n");
#endif
#endif

  /* printf("doing %s\n", pscommand); */

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

#ifdef __sgi
  
  
  totRSS *= getpagesize() / 1024;
#endif

  
  totRSS *= 1024;
  return(totRSS);

}
