#ifndef __MODEL_PARSE_H__
#define __MODEL_PARSE_H__


// the following are for the defines hash table

extern void insert_define( char* defined, char* definition );
extern char* get_definition( char* defined );
#define DEFINED_BUFFER_SIZE 500
extern int is_defined( char* text );

// checks to see if matched word is a reserved word.
extern int res_word( char* text ); 

// removes the reserved word list from memory
extern void kill_lists( void );

#endif
