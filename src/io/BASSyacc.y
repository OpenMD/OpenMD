
/* define some general tokens */

%token MOLECULE ATOM BOND BEND TORSION POSITION MEMBERS CONSTRAINT
%token COMPONENT START_INDEX DEFINED ORIENTATION ZCONSTRAINT RIGIDBODY
%token CUTOFFGROUP

/* more advanced tokens */

%union {
  int i_val;                 /* integer value */
  double d_val;              /* double value  */
  char * s_ptr;              /* string pointer */
  struct node_tag* node_ptr; /* pointer to the statement node tree */
}

%token <i_val> INTEGER
%token <i_val> ARRAY_INDEX

%token <d_val> DOUBLE

%token <s_ptr> IDENTIFIER
%token <s_ptr> QUOTED_STRING
%token <s_ptr> LIST_STRING

%type <node_ptr> stmt
%type <node_ptr> stmt_list
%type <node_ptr> assignment
%type <node_ptr> members
%type <node_ptr> constraint
%type <node_ptr> orientation
%type <node_ptr> position
%type <node_ptr> block
%type <node_ptr> molecule_block
%type <node_ptr> atom_block
%type <node_ptr> bond_block
%type <node_ptr> bend_block
%type <node_ptr> torsion_block
%type <node_ptr> component_block
%type <node_ptr> zconstraint_block
%type <node_ptr> rigidbody_block
%type <node_ptr> cutoffgroup_block
 

%{
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "io/node_list.h"
#include "io/make_nodes.h"
#include "io/parse_tree.h"
#include "utils/simError.h"
#ifdef IS_MPI
#define __is_lex__
#include "mpiBASS.h"
#endif

extern int yylineno;

struct filename_list{
  char my_name[300];
  struct filename_list* next;
};
extern struct filename_list* yyfile_name;

extern void change_in_file( FILE* in_file );
extern void yacc_model( char* file_name );
extern void kill_lists(void);

/* top of the node list */

struct node_tag* head_node;
struct node_tag* current_node;

%}

%%

program: 
	commands
	;

commands: /* NULL */
	| commands stmt		{ 
				  current_node->next_stmt = $2;
				  $2->prev_stmt = current_node;
				  current_node = $2;
				}
	;

stmt:
	  assignment		{ $$ = $1; }
	| members		{ $$ = $1; }
	| constraint		{ $$ = $1; }
	| orientation		{ $$ = $1; }
	| position		{ $$ = $1; }
	| block			{ $$ = $1; }
	;

assignment:
	  IDENTIFIER '=' INTEGER ';'
	  			{ $$ = assign_i( $1, $3 ); }
	| IDENTIFIER '=' DOUBLE ';'
				{ $$ = assign_d( $1, $3 ); }
	| IDENTIFIER '=' IDENTIFIER ';'
				{ $$ = assign_s( $1, $3 ); }
	| IDENTIFIER '=' QUOTED_STRING ';'
				{ $$ = assign_s( $1, $3 ); }
	;

members:
	  MEMBERS LIST_STRING ';'
				{ $$ = members( $2 ); }
	;

constraint:
	  CONSTRAINT LIST_STRING ';'
	  			{ $$ = constraint( $2 ); }
	;

orientation:
	  ORIENTATION LIST_STRING ';'
	  			{ $$ = orientation( $2 ); }
	;

position:
	  POSITION LIST_STRING ';'
	  			{ $$ = position( $2 ); }
	;

block:
	  molecule_block	{ $$ = $1; }
	| atom_block		{ $$ = $1; }
	| bond_block		{ $$ = $1; }
	| bend_block		{ $$ = $1; }
	| torsion_block		{ $$ = $1; }
	| zconstraint_block	{ $$ = $1; }
        | rigidbody_block       { $$ = $1; }
        | cutoffgroup_block     { $$ = $1; }
	| component_block	{ $$ = $1; }
	;

molecule_block:
	MOLECULE '{' stmt_list '}'
				{ $$ = molecule_blk( $3 ); }
	;

atom_block:
	ATOM ARRAY_INDEX '{' stmt_list '}'
				{ $$ = atom_blk( $2, $4 ); }
	;
	
bond_block:
	BOND ARRAY_INDEX '{' stmt_list '}'
				{ $$ = bond_blk( $2, $4 ); }
	;
	
bend_block:
	BEND ARRAY_INDEX '{' stmt_list '}'
				{ $$ = bend_blk( $2, $4 ); }
	;
	
torsion_block:
	TORSION ARRAY_INDEX '{' stmt_list '}'
				{ $$ = torsion_blk( $2, $4 ); }
	;

zconstraint_block:
	ZCONSTRAINT ARRAY_INDEX '{' stmt_list '}'
				{ $$ = zconstraint_blk( $2, $4 ); }
	;

rigidbody_block:
	RIGIDBODY ARRAY_INDEX '{' stmt_list '}'
				{ $$ = rigidbody_blk( $2, $4 ); }
	;
		
cutoffgroup_block:
	CUTOFFGROUP ARRAY_INDEX '{' stmt_list '}'
				{ $$ = cutoffgroup_blk( $2, $4 ); }
	;
		
component_block:
	COMPONENT '{' stmt_list '}'
				{ $$ = component_blk( $3 ); }
	;

stmt_list:
	  stmt			{ $$ = $1; }
	| stmt_list stmt	{ 
				  $1->next_stmt = $2;
				  $2->prev_stmt = $1;
				  $$ = $2;
				}
	;

%%

extern int yyerror( char *err_msg ){

  sprintf( painCave.errMsg, "OOPSE parse error in %s at line %d: %s\n", 
	   yyfile_name->my_name, yylineno + 1, err_msg );
  painCave.isFatal = 1;
  simError();
  return 0;
}

void yacc_BASS( char* file_name ){

  FILE* in_file;

  head_node = (struct node_tag* )malloc( sizeof( node ) );

  head_node->type = GLOBAL_HEAD;
  head_node->index =0;
  head_node->next_stmt = NULL;
  head_node->prev_stmt = NULL;
  head_node->stmt_list = NULL;

  current_node = head_node;

  in_file = fopen( file_name, "r" );
  if( in_file == NULL ){
    sprintf( painCave.errMsg, "yacc error: couldn't open file =>%s\n", 
	     file_name );
    painCave.isFatal = 1;
    simError();
  }

  yyfile_name = 
    (struct filename_list*)malloc( sizeof( struct filename_list ) );
  yyfile_name->next = NULL;
  strcpy( yyfile_name->my_name, file_name );
  change_in_file( in_file );
  
  yyparse();
  
#ifdef IS_MPI
  strcpy( checkPointMsg, "yyParse successful." );
  MPIcheckPoint();
  painCave.isEventLoop = 1;
#endif /* is_mpi*/
  
  fclose( in_file );
  kill_lists(); 
  
  pt_me( head_node );

#ifdef IS_MPI
  painCave.isEventLoop = 0;
#endif /* is_mpi*/
  
  kill_tree( head_node ); 
  head_node = NULL; 
}
