#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "BASSyacc.h"
#include "io/BASS_parse.h"
#include "utils/simError.h"
#ifdef IS_MPI
#define __is_lex__
#include "io/mpiBASS.h"
#endif

#define HASH_SIZE 211 // the size of the hash table
#define SHIFT 4 // the bit shift for the hash index 


//*** Global functions, variables, and structures ************

unsigned short is_initialized = 0; // tells whether to init the linked list

// reserved word elements for the hash table
struct res_element{
  char word[100];        //the reserved word
  int token;             //the token to return;
  struct res_element* next;
};

struct res_element** reserved = NULL; //the table of reserved words

void initialize_res_list(void); // function to initialize the list

int hash ( char* text ); // generates the hash key

void add_res_element( char* text, int the_token ); //adds elements to the table

// the elements of the defined Hash table
struct defined_element{
  char defined[100];
  char definition[DEFINED_BUFFER_SIZE];
  struct defined_element* next;
};

struct defined_element** defined_list = NULL; // the defined hash table



//*** The Functions *******************************************


/*
  function to initialize the list of reserved words into memory.
*/

void initialize_res_list(){
  
  register unsigned short i; // loop counter

  reserved = ( struct res_element ** ) 
    calloc( HASH_SIZE, sizeof( struct hash__element* ) );

  for( i=0; i < HASH_SIZE; i++ ){
    
    reserved[i] = NULL;
  }

  // add the reserved words

  add_res_element( "molecule", MOLECULE );
  add_res_element( "atom", ATOM );
  add_res_element( "bond", BOND );
  add_res_element( "bend", BEND );
  add_res_element( "torsion", TORSION );
  add_res_element( "position", POSITION );
  add_res_element( "members", MEMBERS );
  add_res_element( "constraint", CONSTRAINT );
  add_res_element( "component", COMPONENT );
  add_res_element( "zConstraint", ZCONSTRAINT );
  add_res_element( "orientation", ORIENTATION );
  add_res_element( "rigidBody", RIGIDBODY );
  add_res_element( "cutoffGroup", CUTOFFGROUP );
  
  is_initialized = 1; // set the initialization boolean to true
}


/*
   checks for reserved words.
   If a reserved word is found, returns the token,
   else returns 0.
*/

int res_word( char* text ){
  
  int matched = 0;
  int key; // the hash key
  struct res_element* current_ptr; // points to the current hash element
  struct defined_element* def_ptr; // points to the current define element

  if( !is_initialized ){ // initialize the list if not already done

    initialize_res_list();
  }
  
  // get the hash key

  key = hash( text );
  current_ptr = reserved[key];

  // walk through possible mutiple entries

  while( current_ptr != NULL ){
    
    matched = !strcmp( text, current_ptr->word );
    if( matched ) return current_ptr->token; // search successful

    current_ptr = current_ptr->next;
  }

  // if no keywords turn up, check the word against the list of defines
  
  if( defined_list != NULL ){

    def_ptr = defined_list[key];
    
    while( def_ptr != NULL ){

      matched = !strcmp( text, def_ptr->defined );
      if( matched ) return DEFINED;
      
      def_ptr = def_ptr->next;
    }
  }

  return 0; //search failed
}

/*
 * This is used to test whether a given word is defined
 * returns 1 if true, 0 if False.
 */

int is_defined( char* text ){

  int matched = 0;
  int key;
  struct defined_element* def_ptr; // points to the current define element
  
  key = hash( text );

  if( defined_list != NULL ){
        def_ptr = defined_list[key];
    
    while( def_ptr != NULL ){
      
      matched = !strcmp( text, def_ptr->defined );
      if( matched ) return 1; // search succesful
      
      def_ptr = def_ptr->next;
    }
  }

  return 0; //search failed
}

/*
 * The next bit returns the text to substitute for any given define.
 */

char* get_definition( char* defined ){
  int key;
  int matched = 0;
  struct defined_element* def_ptr;
  
  if( defined_list != NULL ){

    key = hash( defined );
    def_ptr = defined_list[key];
    
    while( def_ptr != NULL ){
      
      matched = !strcmp( defined, def_ptr->defined );
      if( matched ) return def_ptr->definition;
      
      def_ptr = def_ptr->next;
    }
  }
  
  // search failed, therefore there is an error

  sprintf( painCave.errMsg, "%s was not found in the defined list\n", 
	   defined );
  painCave.isFatal = 1;
  simError();
  return NULL;
}
 
/*
 * Add a define statement to the hash table
 */

void insert_define( char* defined, char* definition ){
  
  int i, key;
  struct defined_element* element;

  if( defined_list == NULL ){
    
    defined_list = ( struct defined_element** )
      calloc( HASH_SIZE, sizeof( struct defined_element* ) );
    
    for( i=0; i<HASH_SIZE; i++ ){

      defined_list[i] = NULL;
    }
  }

  key = hash( defined );
  element = (struct defined_element* )
    malloc( sizeof( struct defined_element ) );

  // fill the element
  
  strcpy( element->defined, defined );
  strcpy( element->definition, definition );

  // add the element to the table

  element->next = defined_list[key];
  defined_list[key] = element;
}

/*
 * adds an element to the hash table
 */

void add_res_element( char* text, int the_token ){
  
  int key; // the calculated hash key;
  struct res_element* element; // the element being added

  key = hash( text );
  element = (struct res_element *) malloc( sizeof( struct res_element ) );
  
  // fill the element

  strcpy( element->word, text );
  element->token = the_token;

  // add the element to the table

  element->next = reserved[key];
  reserved[key] = element;
}

/*
  generartes the hash key from the given word.
*/

int hash ( char* text ){

  register unsigned short int i = 0; // loop counter
  int key = 0; // the hash key

  while( text[i] != '\0' ){

    key = ( ( key << SHIFT ) + text[i] ) % HASH_SIZE;
    
    i++;
  }
  
  if( key < 0 ){

    // if the key is less than zero, we've had an overflow error

    sprintf( painCave.errMsg,
	     "Meta-data parse error: There has been an overflow error in the hash key.");
    painCave.isFatal =1;
    simError();
  }
  
  return key;
}

/*
  function to clean up the memory after we are finished with the lists
*/

void kill_lists( ){
  
  struct res_element* current_res_ptr;
  struct res_element* temp_res_ptr;

  struct defined_element* current_def_ptr;
  struct defined_element* temp_def_ptr;

  register unsigned short int i;

  if( reserved != NULL ){

    for( i=0; i < HASH_SIZE; i++){
      
      current_res_ptr = reserved[i];
      
      while( current_res_ptr != NULL ){
	
	temp_res_ptr = current_res_ptr->next;
	free( current_res_ptr );
	current_res_ptr = temp_res_ptr;
      }
    }    

    free( reserved );
  }

  if( defined_list != NULL ){

    for( i=0; i < HASH_SIZE; i++){
      
      current_def_ptr = defined_list[i];
      
      while( current_def_ptr != NULL ){
	
	temp_def_ptr = current_def_ptr->next;
	free( current_def_ptr );
	current_def_ptr = temp_def_ptr;
      }
    }    

    free( defined_list );
  }

  reserved = NULL;
  defined_list = NULL;

  is_initialized = 0; // initialization is now false
}
