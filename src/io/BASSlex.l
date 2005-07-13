
/* define some search patterns */

digit			[0-9]
letter			[A-Za-z]
natural			{digit}+
signedNat		("+"|"-")?{natural}
number			{signedNat}("."{natural})?([Ee]{signedNat})?
identifier		("_"|{letter})({letter}|{digit}|"-"|"_")*

/* define start states */

%x PRE_P
%x INCL
%x DEF
%x IFDEF
%x IFNDEF
%x ESCAPE
%x CHECK_ESCAPE

/* what to put at the top of the lex code */

%{
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "io/BASSyacc.h"
#include "io/BASS_parse.h"
#include "utils/simError.h"
#ifdef IS_MPI
#define __is_lex__
#include "mpiBASS.h"
#endif

typedef unsigned short int r_short;
typedef	unsigned short	my_short;

extern void change_in_file( FILE* in_file );

/* the following is used by the include start state*/

#define MAX_BUFFER_DEPTH 10
YY_BUFFER_STATE buffer_stack[MAX_BUFFER_DEPTH]; /* a stack of the include buffers*/
int buffer_stack_ptr = 0;
struct filename_list{
  char my_name[300];
  struct filename_list* next;
};
struct filename_list* yyfile_name;
struct filename_list* temp_yyfile_name;
int yylineno_stack[MAX_BUFFER_DEPTH];


/* the following is a check against the define buffer length*/

void check_def_buff( char* defined, int index );

/*these are used by the ifdef and ifndef statements*/

int escape_stack_ptr = 0; /*keeps track of the escape stack*/


%}

%option yylineno

%%

{signedNat}			{
				  yylval.i_val = atoi( yytext );
				  return INTEGER;
				}

{number}			{
				  yylval.d_val = atof( yytext );
				  return DOUBLE;
				}

{identifier}			{
				  int token;
				  token = res_word( yytext );
				  
				  if( token == DEFINED ){
				  
				    if( buffer_stack_ptr >= MAX_BUFFER_DEPTH ){
				      fprintf( stderr,
				               "Maximum buffer depth exceeded for %s.\n",
					       yytext );
				      exit(1);
				    }
				    
				    buffer_stack[buffer_stack_ptr] = YY_CURRENT_BUFFER;
				    buffer_stack_ptr++;

				    yy_scan_string( get_definition( yytext ) );
				  }
				  else if( token ){
				    return token;
				  }
				  else{
				    yylval.s_ptr = strdup( yytext );
				    return IDENTIFIER;
				  }
				}

\".*\"				{
				  /* little routine to strip off the quotes */
				  
				  my_short i;
				  i = 0; /* index*/
				  while( yytext[i+1] != '\"' ){

				    yytext[i] = yytext[i+1];
				    i++;
				  }
				  yytext[i] = '\0';
				  yylval.s_ptr = strdup( yytext );
				  return QUOTED_STRING;
				}

\(.*\)				{
				  /* little routine to pass a list string */
				  
				  my_short i;
				  i = 0; /* index*/
				  while( yytext[i+1] != ')' ){

				    yytext[i] = yytext[i+1];
				    i++;
				  }
				  yytext[i] = '\0';
				  yylval.s_ptr = strdup( yytext );
				  return LIST_STRING;
				}

\[{natural}\]			{
				  int index;
				  sscanf(yytext, "[%d]", &index);
				  yylval.i_val = index;
				  return ARRAY_INDEX;
				}

[ \t\n]+			/* ignore whitespace */;

"/*"				{
				  /* ignore comments */
				  my_short done;
				  char c;

				  done = 0;
				  while( !done ){
				    
				    c = input();
				    while( c != '*' ){
				      
				      c = input();
				    }
				    while( c == '*' ){
				      
				      c = input();
				    }
				    if( c == '/' ) done = 1;
				  }
				}

"//".*\n			/* ignore comments */;

"#"				BEGIN(PRE_P);

.				{
				  /* pass everything else to yacc*/
				  return yytext[0];
				}

<PRE_P>"include"		BEGIN(INCL);
<PRE_P>"define"			BEGIN(DEF);
<PRE_P>"ifdef"			BEGIN(IFDEF);
<PRE_P>"ifndef"			BEGIN(IFNDEF);
<PRE_P>"endif"			/* do nothing */;
<PRE_P>\n			BEGIN(INITIAL);

<INCL>[ \t]*			/* eat white space */
<INCL>\".*\"			{
				  char foo_name[300];
				  
				  /* little routine to strip off the quotes */

				  my_short i;
				  i = 0; /* index*/
				  while( yytext[i+1] != '\"' ){

				    yytext[i] = yytext[i+1];
				    i++;
				  }
				  
				  yytext[i] = '\0';
				  strcpy( foo_name, yytext );

				  /* now we have the include file name*/
				  
				  if( buffer_stack_ptr >= MAX_BUFFER_DEPTH ){
				    
				    fprintf( stderr, "Includes nested too deeply\n" );
				    exit(1);
				  }

				  buffer_stack[buffer_stack_ptr] = YY_CURRENT_BUFFER;
				  yylineno_stack[buffer_stack_ptr] = yylineno;
				  buffer_stack_ptr++;

				  yyin = fopen( yytext, "r" );
				  if( yyin == NULL ){
				    fprintf( stderr, "Unable to include file %s\n", yytext );
				    exit(1);
				  }

				  yy_switch_to_buffer( yy_create_buffer( yyin, YY_BUF_SIZE ) );
				  yylineno = 0;
				  
				  temp_yyfile_name = (struct filename_list* )malloc( sizeof( struct filename_list ) );
				  temp_yyfile_name->next = yyfile_name;
				  yyfile_name = temp_yyfile_name;
				  strcpy( yyfile_name->my_name, foo_name );


				  BEGIN(INITIAL);
				}
<INCL>\n			BEGIN(INITIAL);
<INCL>.				/* ignore everything else */

<DEF>[ \t]*			/* eat white space */;
<DEF>{identifier}		{
				  char c;
				  char definition[ DEFINED_BUFFER_SIZE ];
				  short int done;
				  short int c_done; /* a done marker for the comments*/
				  int def_ptr;
				  
				  /* initialize the definition buffer*/

				  for( def_ptr = 0; def_ptr < DEFINED_BUFFER_SIZE; def_ptr++ ){
				    definition[def_ptr] = '\0';
				  }
				  
				  def_ptr =0;
				  done = 0;
				  c_done = 0;
				 
				  while( !done ){
				    
				    c = input();
				    if( c == '\"' ){
				    
				      /* shove the whole quoted string into the macro*/
				    
				      definition[def_ptr] = c;
				      /*fprintf( stderr, "%c", c );*/
				      def_ptr++;
				      check_def_buff( yytext, def_ptr );
				      
				      c = input();
				      while( c != '\"' ){
				      
				        definition[def_ptr] = c;
				        /*fprintf( stderr, "%c", c );*/
					def_ptr++;
				        check_def_buff( yytext, def_ptr );

					c = input();
				      }
				      definition[def_ptr] = c;
				      /*fprintf( stderr, "%c", c );*/
				      def_ptr++;
				      check_def_buff( yytext, def_ptr );
				      c = input();
				    }
				    
				    /* handle comments*/

				    if( c == '/' ){
				      c = input();
				      switch( c ){

				        case '/':
				          while( c != '\n' && c != '\\' ){
					    c = input();
					  }
					  break;
				      	
					case '*':
					  c_done = 0;
					  while( !c_done ){
					    c = input();
					    while( c != '*' ){
					      c = input();
					    }
					    while( c == '*' ){
					      c = input();
					    }
					    if( c == '/' ) c_done = 1;
					  }
					  c = input();
					  break;

					default:
					  /* the '/' char was a normal symbol*/
					  definition[def_ptr] = '/';
				          /*fprintf( stderr, "%c", c );*/
					  def_ptr++;
					  check_def_buff( yytext, def_ptr );
					  break;
				      }
				    }
					
				    
				    if( c == '\n' ){
				      done = 1;
				    }
				    
				    else{
				      
				      /* check for the line wrap character '\'*/
				      
				      if( c == '\\' ){
				        
					/* skip the rest of the line until the line return*/

					c = input();
					while( c != '\n' ){
					  c = input();
					}
				      }

				      else{
				        
					/* we now know the character is a good one */

					definition[def_ptr] = c;
				        /*fprintf( stderr, "%c", c );*/
					def_ptr++;
					check_def_buff( yytext, def_ptr );
				      }
				    }
				  }  
				  
				  insert_define( yytext, definition );
				  BEGIN(INITIAL);
				}
<DEF>\n				BEGIN(INITIAL);
<DEF>.				/* ignore everything else */;

<IFDEF>[ \t]*			/* eat white space */;
<IFDEF>{identifier}		{
				  if( !is_defined( yytext ) ){
				    escape_stack_ptr++;
				    BEGIN(ESCAPE);
				  }
				}
<IFDEF>\n			BEGIN(INITIAL);
<IFDEF>.			/* ignore everything else */;

<IFNDEF>[ \t]*			/* eat the white space */;
<IFNDEF>{identifier}		{
				  if( is_defined( yytext ) ){
				    escape_stack_ptr++;
				    BEGIN(ESCAPE);
				  }
				}
<IFNDEF>\n			BEGIN(INITIAL);
<IFNDEF>.			/* ignore everything else */;

<ESCAPE>\".*\"			/* do nothing */;
<ESCAPE>"//".*\n		/* ignore comments */;
<ESCAPE>"/*"			{
                                  /* ignore comments */
				  my_short done;
				  char c;

				  done = 0;
				  while( !done ){

				    c = input();
				    while( c != '*' ){
				    
				      c = input();
				    }
				    
				    while( c == '*' ){
				    
				      c = input();
				    }
				    
				    if( c == '/' ) done = 1;
				  }
				}
<ESCAPE>"#"			BEGIN(CHECK_ESCAPE);
<ESCAPE>.			/* ignore everything else */;

<CHECK_ESCAPE>[ \t]		/* ignore whitespace */;
<CHECK_ESCAPE>"ifdef"		{ escape_stack_ptr++; BEGIN(ESCAPE); }
<CHECK_ESCAPE>"ifndef"		{ escape_stack_ptr++; BEGIN(ESCAPE); }
<CHECK_ESCAPE>"endif"		{
				  escape_stack_ptr--;
				  if( escape_stack_ptr <= 0){
				    escape_stack_ptr = 0; /* just in case something flubbed */
				    BEGIN(INITIAL);
				  }
				}
<CHECK_ESCAPE>\n		BEGIN(ESCAPE);
<CHECK_ESCAPE>.			/* ignore everything else */;

<<EOF>> 			{
				  buffer_stack_ptr--;
				  if( buffer_stack_ptr < 0 ){
                                    /*reinitialize buffer_stack_ptr*/
				    buffer_stack_ptr = 0; 
				    yyterminate();
				  }
				  
				  else{
				  
				    yy_delete_buffer( YY_CURRENT_BUFFER );
				    yy_switch_to_buffer( buffer_stack[buffer_stack_ptr] );
				    yylineno = yylineno_stack[buffer_stack_ptr];

				    temp_yyfile_name = yyfile_name;
				    yyfile_name = temp_yyfile_name->next;
				    free( temp_yyfile_name );
				  }
				}
				
%%

void check_def_buff( char* defined, int index ){
  
  if( index >= DEFINED_BUFFER_SIZE ){
    
    sprintf( painCave.errMsg, "Define buffer size exceeded for %s\n", defined );
    painCave.isFatal = 1;
    simError();
  }
}

int yywrap(void){
  return 1;
}

void change_in_file( FILE* in_file ){

  yyin = in_file;
}
