#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>

struct coords{
  double pos[3]; // cartesian coords
  double q[4];  // the quanternions
  char name[30];
};

struct linked_xyz{
  struct coords *r;
  double time;
  double Hmat[3][3];
  double HmatI[3][3];
  struct linked_xyz *next;
};

char *program_name; /*the name of the program */
int printWater = 0;
int printDipole = 0;

void usage(void);

void rotWrite( FILE* out_file, struct coords* theSSD );

void setRot( double A[3][3], double q[4] );

void rotMe( double A[3][3], double r[3] );
void wrapVector( double thePos[3], double Hmat[3][3], double HmatI[3][3]);
double matDet3(double a[3][3]);
void invertMat3(double a[3][3], double b[3][3]);
void matVecMul3(double m[3][3], double inVec[3], double outVec[3]);
double roundMe(double x);

int main(argc, argv)
     int argc;
     char *argv[];
{


  struct coords *out_coords;

  int i, j, k, l, m; /* loop counters */
  mode_t dir_mode = S_IRWXU;

  int lineCount = 0;  //the line number
  int newN; // the new number of atoms
  int nSSD=0; // the number of SSD per frame

  unsigned int nAtoms; /*the number of atoms in each time step */
  char read_buffer[1000]; /*the line buffer for reading */
  char *eof_test; /*ptr to see when we reach the end of the file */
  char *foo; /*the pointer to the current string token */
  FILE *in_file; /* the input file */
  FILE *out_file; /*the output file */
  char *out_prefix = NULL; /*the prefix of the output file */
  char out_name[500]; /*the output name */
  int have_outName =0;
  char in_name[500]; // the input file name
  char temp_name[500];
  char *root_name = NULL; /*the root name of the input file */
  unsigned int n_out = 0; /*keeps track of which output file is being written*/


  int skipFrames = 1;

  struct linked_xyz *current_frame;
  struct linked_xyz *temp_frame;

  int nframes =0;

  int wrap = 0;

  double cX = 0.0;
  double cY = 0.0;
  double cZ = 0.0;

  double mcX, mcY, mcZ;
  double newCX, newCY, newCZ;
  double dx, dy, dz;
  int mcount;

  int repeatX = 0;
  int repeatY = 0;
  int repeatZ = 0;

  int nRepeatX = 0;
  int nRepeatY = 0;
  int nRepeatZ = 0;

  struct coords* rCopy;
  
  int done;
  char current_flag;

  program_name = argv[0]; /*save the program name in case we need it*/
  
  for( i = 1; i < argc; i++){
    
    if(argv[i][0] =='-'){
      
      if(argv[i][1] == '-' ){
	
	// parse long word options
	
	if( !strcmp( argv[i], "--repeatX" ) ){
	  repeatX = 1;
	  i++;
	  nRepeatX = atoi( argv[i] );
	}

	else if( !strcmp( argv[i], "--repeatY" ) ){
	  repeatY = 1;
	  i++;
	  nRepeatY = atoi( argv[i] );
	}

	else if( !strcmp( argv[i], "--repeatZ" ) ){
	  repeatZ = 1;
	  i++;
	  nRepeatZ = atoi( argv[i] );
	}

	else{
	  fprintf( stderr, 
		   "Invalid option \"%s\"\n", argv[i] );
	  usage();
	}
      }
      
      else{

	// parse single character options
	
	done =0;
	j = 1;
	current_flag = argv[i][j];
	while( (current_flag != '\0') && (!done) ){
	  
	  switch(current_flag){
	    
	  case 'o':
	    // -o <prefix> => the output
	    
	    i++;
	    strcpy( out_name, argv[i] );
	    have_outName = 1;
	    done = 1;
	    break;
	    
	  case 'n':
	    // -n <#> => print every <#> frames
	    
	    i++;
	    skipFrames = atoi(argv[i]);
	    done = 1;
	    break;
	    
	  case 'h':
	    // -h => give the usage
	    
	    usage();
	    break;
	    
	  case 'd':
	    // print the dipoles
	    
	    printDipole = 1;
	    break;
	    
	  case 'w':
	    // print the waters
	    
	    printWater = 1;
	    break;
	    
	  case 'm':
	    // map to a periodic box
	    
	    wrap = 1;
	    break;
	    
	  default:
	    
	    fprintf(stderr, "Bad option \"-%c\"\n", current_flag);
	    usage();
	  }
	  j++;
	  current_flag = argv[i][j];
	}
      }
    }
    
    else{
      
      if( root_name != NULL ){
	fprintf( stderr, 
		 "Error at \"%s\", program does not currently support\n"
		 "more than one input file.\n"
		 "\n",
		 argv[i]);
	usage();
      }
      
      root_name = argv[i];
    }
  }
  
  if(root_name == NULL){
    usage();
  }

  sprintf( in_name, "%s", root_name );
  if( !have_outName ) sprintf( out_name, "%s.xyz", root_name );

  in_file = fopen(in_name, "r");
  if(in_file == NULL){
    printf("Cannot open file \"%s\" for reading.\n", in_name);
    exit(8);
  }
  
  out_file = fopen( out_name, "w" );
  if( out_file == NULL ){
    printf("Cannot open file \"%s\" for writing.\n", out_name);
    exit(8);
  }

  // start reading the first frame

  eof_test = fgets(read_buffer, sizeof(read_buffer), in_file);
  lineCount++;
  
  current_frame = (struct linked_xyz *)malloc(sizeof(struct linked_xyz));
  current_frame->next = NULL;
  

  while(eof_test != NULL){
    
    (void)sscanf(read_buffer, "%d", &nAtoms);
    current_frame->r = 
      (struct coords *)calloc(nAtoms, sizeof(struct coords));

    // read and the comment line and grab the time and box dimensions

    eof_test = fgets(read_buffer, sizeof(read_buffer), in_file);
    lineCount++;
    if(eof_test == NULL){
      printf("error in reading file at line: %d\n", lineCount);
      exit(8);
    }

    foo = strtok( read_buffer, " ,;\t" );
    sscanf( read_buffer, "%lf", &current_frame->time );
    
    foo = strtok(NULL, " ,;\t");
    if(foo == NULL){
      printf("error in reading file at line: %d\n", lineCount);
      exit(8);
    }
    (void)sscanf(foo, "%lf",&current_frame->Hmat[0][0]);

    foo = strtok(NULL, " ,;\t");
    if(foo == NULL){
      printf("error in reading file at line: %d\n", lineCount);
      exit(8);
    }
    (void)sscanf(foo, "%lf",&current_frame->Hmat[1][0]);

    foo = strtok(NULL, " ,;\t");
    if(foo == NULL){
      printf("error in reading file at line: %d\n", lineCount);
      exit(8);
    }
    (void)sscanf(foo, "%lf",&current_frame->Hmat[2][0]);

    foo = strtok(NULL, " ,;\t");
    if(foo == NULL){
      printf("error in reading file at line: %d\n", lineCount);
      exit(8);
    }
    (void)sscanf(foo, "%lf",&current_frame->Hmat[0][1]);

    foo = strtok(NULL, " ,;\t");
    if(foo == NULL){
      printf("error in reading file at line: %d\n", lineCount);
      exit(8);
    }
    (void)sscanf(foo, "%lf",&current_frame->Hmat[1][1]);

    foo = strtok(NULL, " ,;\t");
    if(foo == NULL){
      printf("error in reading file at line: %d\n", lineCount);
      exit(8);
    }
    (void)sscanf(foo, "%lf",&current_frame->Hmat[2][1]);

    foo = strtok(NULL, " ,;\t");
    if(foo == NULL){
      printf("error in reading file at line: %d\n", lineCount);
      exit(8);
    }
    (void)sscanf(foo, "%lf",&current_frame->Hmat[0][2]);

    foo = strtok(NULL, " ,;\t");
    if(foo == NULL){
      printf("error in reading file at line: %d\n", lineCount);
      exit(8);
    }
    (void)sscanf(foo, "%lf",&current_frame->Hmat[1][2]);

    foo = strtok(NULL, " ,;\t");
    if(foo == NULL){
      printf("error in reading file at line: %d\n", lineCount);
      exit(8);
    }
    (void)sscanf(foo, "%lf",&current_frame->Hmat[2][2]);


    // Find HmatI:

    invertMat3(current_frame->Hmat, current_frame->HmatI);


    for( i=0; i < nAtoms; i++){
      
      eof_test = fgets(read_buffer, sizeof(read_buffer), in_file);
      lineCount++;
      if(eof_test == NULL){
	printf("error in reading file at line: %d\n", lineCount);
	exit(8);
      }

      foo = strtok(read_buffer, " ,;\t");
      if ( !strcmp( "SSD", foo ) ) nSSD++;
      (void)strcpy(current_frame->r[i].name, foo); //copy the atom name 
			  
      // next we grab the positions 
      
      foo = strtok(NULL, " ,;\t");
      if(foo == NULL){
	printf("error in reading postition x from %s\n"
	       "natoms  = %d, line = %d\n",
	       in_name, nAtoms, lineCount );
	exit(8);
      }
      (void)sscanf( foo, "%lf", &current_frame->r[i].pos[0] );
      
      foo = strtok(NULL, " ,;\t");
      if(foo == NULL){
	printf("error in reading postition y from %s\n"
	     "natoms  = %d, line = %d\n",
	       in_name, nAtoms, lineCount );
	exit(8);
      }
      (void)sscanf( foo, "%lf", &current_frame->r[i].pos[1] );
      
      foo = strtok(NULL, " ,;\t");
      if(foo == NULL){
	printf("error in reading postition z from %s\n"
	       "natoms  = %d, line = %d\n",
	       in_name, nAtoms, lineCount );
	exit(8);
      }
      (void)sscanf( foo, "%lf", &current_frame->r[i].pos[2] );
    
      // get the velocities
      
      foo = strtok(NULL, " ,;\t");
      if(foo == NULL){
	printf("error in reading velocity x from %s\n"
	       "natoms  = %d, line = %d\n",
	       in_name, nAtoms, lineCount );
	exit(8);
      }
      
    
      foo = strtok(NULL, " ,;\t");
      if(foo == NULL){
	printf("error in reading velocity y from %s\n"
	       "natoms  = %d, line = %d\n",
	       in_name, nAtoms, lineCount );
	exit(8);
      }
      
      
      foo = strtok(NULL, " ,;\t");
      if(foo == NULL){
	printf("error in reading velocity z from %s\n"
	       "natoms  = %d, line = %d\n",
	       in_name, nAtoms, lineCount );
	exit(8);
      }
      
      // get the quaternions
      
      foo = strtok(NULL, " ,;\t");
      if(foo == NULL){
	printf("error in reading quaternion 0 from %s\n"
	       "natoms  = %d, line = %d\n",
	       in_name, nAtoms, lineCount );
	exit(8);
      }
      (void)sscanf( foo, "%lf", &current_frame->r[i].q[0] );
      
      foo = strtok(NULL, " ,;\t");
      if(foo == NULL){
	printf("error in reading quaternion 1 from %s\n"
	       "natoms  = %d, line = %d\n",
	       in_name, nAtoms, lineCount );
	exit(8);
      }
      (void)sscanf( foo, "%lf", &current_frame->r[i].q[1] );
      
      foo = strtok(NULL, " ,;\t");
      if(foo == NULL){
	printf("error in reading quaternion 2 from %s\n"
	       "natoms  = %d, line = %d\n",
	       in_name, nAtoms, lineCount );
	exit(8);
      }
      (void)sscanf( foo, "%lf", &current_frame->r[i].q[2] );

      foo = strtok(NULL, " ,;\t");
      if(foo == NULL){
	printf("error in reading quaternion 3 from %s\n"
	       "natoms  = %d, line = %d\n",
	       in_name, nAtoms, lineCount );
	exit(8);
      }
      (void)sscanf( foo, "%lf", &current_frame->r[i].q[3] );
      
      // get the angular velocities
      
      foo = strtok(NULL, " ,;\t");
      if(foo == NULL){
	printf("error in reading angular momentum jx from %s\n"
	       "natoms  = %d, line = %d\n",
	       in_name, nAtoms, lineCount );
	exit(8);
      }
      
      foo = strtok(NULL, " ,;\t");
      if(foo == NULL){
	printf("error in reading angular momentum jy from %s\n"
	       "natoms  = %d, line = %d\n",
	       in_name, nAtoms, lineCount );
	exit(8);
      }
      
      foo = strtok(NULL, " ,;\t");
      if(foo == NULL){
	printf("error in reading angular momentum jz from %s\n"
	       "natoms  = %d, line = %d\n",
	       in_name, nAtoms, lineCount );
	exit(8);
      }
    }

    // if necassary, wrap into the periodic image
    
    if(wrap){
      
      for( i=0; i<nAtoms; i++ ){
	
	/*
	if( !strcmp( current_frame->r[i].name, "HEAD" ) ){
	  
	  // find the center of the lipid
	  
	  mcX = 0.0;
	  mcY = 0.0;
	  mcZ = 0.0;
	  mcount = 0;

	  mcX += current_frame->r[i].x;
	  mcY += current_frame->r[i].y;
	  mcZ += current_frame->r[i].z;
	  mcount++;
	  i++;

	  while( strcmp( current_frame->r[i].name, "HEAD" ) &&
		 strcmp( current_frame->r[i].name, "SSD" ) ){
	    
	    mcX += current_frame->r[i].x;
	    mcY += current_frame->r[i].y;
	    mcZ += current_frame->r[i].z;
	    mcount++;
	    i++;
	  }
	  
	  mcX /= mcount;
	  mcY /= mcount;
	  mcZ /= mcount;

	  newCX = mcX;
	  newCY = mcY;
	  newCZ = mcZ;
	  
	  map( &newCX, &newCY, &newCZ,
	       cX, cY, cZ,
	       current_frame->boxX,
	       current_frame->boxY,
	       current_frame->boxZ );
	  
	  dx = mcX - newCX;
	  dy = mcY - newCY;
	  dz = mcZ - newCZ;

	  for(j=(i-mcount); j<mcount; j++ ){
	    
	    current_frame->r[j].x -= dx;
	    current_frame->r[j].y -= dy;
	    current_frame->r[j].z -= dz;
	  }
	  
	  i--; // so we don't skip the next atom
	}
	*/
	// else 
        
        wrapVector(current_frame->r[i].pos, 
                   current_frame->Hmat, 
                   current_frame->HmatI);
	
      }
    }   

    nframes++; 
    
    // write out here    
    
    if(printWater) newN = nAtoms + ( nSSD * 3 ); 
    else newN = nAtoms - nSSD;
    
    newN *= (nRepeatX+1);
    newN *= (nRepeatY+1);
    newN *= (nRepeatZ+1);

    if( !(nframes%skipFrames) ){
      fprintf( out_file,	
  	       "%d\n"
	       "%lf;\t%lf\t%lf\t%lf;\t%lf\t%lf\t%lf;\t%lf\t%lf\t%lf;\n",
	       newN, current_frame->time,
	       current_frame->Hmat[0][0], current_frame->Hmat[1][0],
	       current_frame->Hmat[2][0],
	       current_frame->Hmat[0][1], current_frame->Hmat[1][1],
	       current_frame->Hmat[2][1],
	       current_frame->Hmat[0][2], current_frame->Hmat[1][2],
	       current_frame->Hmat[2][2] );
    
      rCopy = (struct coords*)
	calloc( nAtoms, sizeof( struct coords ));

      for(i=0; i<nAtoms; i++){
	rCopy[i].q[0] = current_frame->r[i].q[0];
	rCopy[i].q[1] = current_frame->r[i].q[1];
	rCopy[i].q[2] = current_frame->r[i].q[2];
	rCopy[i].q[3] = current_frame->r[i].q[3];

	strcpy(rCopy[i].name, current_frame->r[i].name);
      }

      for(j=0; j<(nRepeatX+1); j++){
	for(k=0; k<(nRepeatY+1); k++){
	  for(l=0; l<(nRepeatZ+1); l++){

	    for(i=0; i<nAtoms; i++){

              for(m = 0; m < 3; m++) {
                rCopy[i].pos[m] = current_frame->r[i].pos[m] +
                  j * current_frame->Hmat[m][0] + 
                  k * current_frame->Hmat[m][1] + 
                  l * current_frame->Hmat[m][2];
              }
              
	      if( !strcmp( "SSD", rCopy[i].name ) ){
		
		rotWrite( out_file, &rCopy[i] );
	      }
	      
	      else if( !strcmp( "HEAD", rCopy[i].name ) ){
		
		rotWrite( out_file, &rCopy[i] );
	      }
	      else{
		//modify
		fprintf( out_file,
			 "%s\t%lf\t%lf\t%lf\t0.0\t0.0\t0.0\n",
			 rCopy[i].name,
			 rCopy[i].pos[0],
			 rCopy[i].pos[1],
			 rCopy[i].pos[2] );
	      }
	    }
	  }
	}
      }
      
      free( rCopy );
    }
    
    /*free up memory */
    
    temp_frame = current_frame->next;
    current_frame->next = NULL;
    
    if(temp_frame != NULL){
      
      free(temp_frame->r);
      free(temp_frame);
    }
    
    /* make a new frame */
    
    temp_frame = (struct linked_xyz *)malloc(sizeof(struct linked_xyz));
    temp_frame->next = current_frame;
    current_frame = temp_frame;
    
    nSSD = 0;

    eof_test = fgets(read_buffer, sizeof(read_buffer), in_file);
    lineCount++;
  }
  
  (void)fclose(in_file);
  
  return 0;
  
}



void rotWrite( FILE* out_file, struct coords* theDipole ){

  double u[3];
  double h1[3];
  double h2[3];
  double ox[3];
  
  double A[3][3];
  
 
  u[0] = 0.0; u[1] = 0.0; u[2] = 1.0;

  h1[0] = 0.0; h1[1] = -0.75695; h1[2] = 0.5206;

  h2[0] = 0.0; h2[1] = 0.75695;  h2[2] = 0.5206;

  ox[0] = 0.0; ox[1] = 0.0;      ox[2] = -0.0654;

  
  setRot( A, theDipole->q );
  rotMe( A, u );

  if( !strcmp( "SSD", theDipole->name )){
    
    rotMe( A, h1 );
    rotMe( A, h2 );
    rotMe( A, ox );

    if( printWater ){
      
      if(printDipole) {
	
	fprintf( out_file,
		 "X\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n",
		 theDipole->pos[0],
		 theDipole->pos[1],
		 theDipole->pos[2],
		 u[0], u[1], u[2] );

        fprintf( out_file,
                 "O\t%lf\t%lf\t%lf\t0.0\t0.0\t0.0\n",
                 theDipole->pos[0] + ox[0],
                 theDipole->pos[1] + ox[1],
                 theDipole->pos[2] + ox[2] );
        
        fprintf( out_file,
                 "H\t%lf\t%lf\t%lf\t0.0\t0.0\t0.0\n",
                 theDipole->pos[0] + h1[0],
                 theDipole->pos[1] + h1[1],
                 theDipole->pos[2] + h1[2] );
        
        fprintf( out_file,
                 "H\t%lf\t%lf\t%lf\t0.0\t0.0\t0.0\n",
                 theDipole->pos[0] + h2[0],
                 theDipole->pos[1] + h2[1],
                 theDipole->pos[2] + h2[2] );
	}
      else{
	
	fprintf( out_file,
		 "X\t%lf\t%lf\t%lf\n",
		 theDipole->pos[0],
		 theDipole->pos[1],
		 theDipole->pos[2] );

        fprintf( out_file,
                 "O\t%lf\t%lf\t%lf\n",
                 theDipole->pos[0] + ox[0],
                 theDipole->pos[1] + ox[1],
                 theDipole->pos[2] + ox[2] );
        
        fprintf( out_file,
                 "H\t%lf\t%lf\t%lf\n",
                 theDipole->pos[0] + h1[0],
                 theDipole->pos[1] + h1[1],
                 theDipole->pos[2] + h1[2] );
        
        fprintf( out_file,
                 "H\t%lf\t%lf\t%lf\n",
                 theDipole->pos[0] + h2[0],
                 theDipole->pos[1] + h2[1],
                 theDipole->pos[2] + h2[2] );
      }     
    }
  }
  
  else if( !strcmp( "HEAD", theDipole->name )) {
    
    if( printDipole ){
      
      fprintf( out_file,
	       "HEAD\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n",
	       theDipole->pos[0],
	       theDipole->pos[1],
	       theDipole->pos[2],
	       u[0], u[1], u[2] );
    }
    else{
     
      fprintf( out_file,
	       "HEAD\t%lf\t%lf\t%lf\n",
	       theDipole->pos[0],
	       theDipole->pos[1],
	       theDipole->pos[2] );
    }
  }
}


void setRot( double A[3][3], double the_q[4] ){

  double q0Sqr, q1Sqr, q2Sqr, q3Sqr;
  
  q0Sqr = the_q[0] * the_q[0];
  q1Sqr = the_q[1] * the_q[1];
  q2Sqr = the_q[2] * the_q[2];
  q3Sqr = the_q[3] * the_q[3];
  
  A[0][0] = q0Sqr + q1Sqr - q2Sqr - q3Sqr;
  A[0][1] = 2.0 * ( the_q[1] * the_q[2] + the_q[0] * the_q[3] );
  A[0][2] = 2.0 * ( the_q[1] * the_q[3] - the_q[0] * the_q[2] );
  
  A[1][0]  = 2.0 * ( the_q[1] * the_q[2] - the_q[0] * the_q[3] );
  A[1][1] = q0Sqr - q1Sqr + q2Sqr - q3Sqr;
  A[1][2] = 2.0 * ( the_q[2] * the_q[3] + the_q[0] * the_q[1] );
  
  A[2][0] = 2.0 * ( the_q[1] * the_q[3] + the_q[0] * the_q[2] );
  A[2][1] = 2.0 * ( the_q[2] * the_q[3] - the_q[0] * the_q[1] );
  A[2][2] = q0Sqr - q1Sqr -q2Sqr +q3Sqr;
}



void rotMe( double A[3][3], double r[3] ){
  
  int i,j;
  double rb[3]; // the body frame vector 
 
  for(i=0; i<3; i++ ){
    rb[i] = r[i];
  }
  
  for( i=0; i<3; i++ ){
    r[i] = 0.0;
    
    for( j=0; j<3; j++ ){
      r[i] += A[j][i] * rb[j];
    }
  }
}



/***************************************************************************
 * prints out the usage for the command line arguments, then exits.
 ***************************************************************************/

void usage(){
  (void)fprintf(stderr, 
		"The proper usage is: %s [options] <dumpfile>\n"
		"\n"
		"Options:\n"
		"\n"
		"   -h              Display this message\n"
		"   -o <out_name>   The output file (Defaults to <dumpfile>.xyz)\n"
		"   -d              print the dipole moments\n"
		"   -w              print the waters\n"
		"   -m              map to the periodic box\n"
		"   -n <#>          print every <#> frames\n" 
		"\n"
		"  --repeatX <#>    The number of images to repeat in the x direction\n"
		"  --repeatY <#>    The number of images to repeat in the y direction\n"
		"  --repeatZ <#>    The number of images to repeat in the z direction\n"
		
		"\n",
		
		program_name);
  exit(8);
}

void wrapVector( double thePos[3], double Hmat[3][3], double HmatInv[3][3]){

  int i;
  double scaled[3];

  // calc the scaled coordinates.
  
  matVecMul3(HmatInv, thePos, scaled);
  
  for(i=0; i<3; i++)
    scaled[i] -= roundMe(scaled[i]);
  
  // calc the wrapped real coordinates from the wrapped scaled coordinates
  
  matVecMul3(Hmat, scaled, thePos);
    
}

double matDet3(double a[3][3]) {
  int i, j, k;
  double determinant;

  determinant = 0.0;

  for(i = 0; i < 3; i++) {
    j = (i+1)%3;
    k = (i+2)%3;

    determinant += a[0][i] * (a[1][j]*a[2][k] - a[1][k]*a[2][j]);
  }
  return determinant;
}

void invertMat3(double a[3][3], double b[3][3]) {

  int  i, j, k, l, m, n;
  double determinant;

  determinant = matDet3( a );

  if (determinant == 0.0) {
    printf("Can't invert a matrix with a zero determinant!\n");
  }

  for (i=0; i < 3; i++) {
    j = (i+1)%3;
    k = (i+2)%3;
    for(l = 0; l < 3; l++) {
      m = (l+1)%3;
      n = (l+2)%3;

      b[l][i] = (a[j][m]*a[k][n] - a[j][n]*a[k][m]) / determinant;
    }
  }
}


void matVecMul3(double m[3][3], double inVec[3], double outVec[3]) {
  double a0, a1, a2;

  a0 = inVec[0];  a1 = inVec[1];  a2 = inVec[2];

  outVec[0] = m[0][0]*a0 + m[0][1]*a1 + m[0][2]*a2;
  outVec[1] = m[1][0]*a0 + m[1][1]*a1 + m[1][2]*a2;
  outVec[2] = m[2][0]*a0 + m[2][1]*a1 + m[2][2]*a2;
}

double roundMe( double x ){
  return ( x >= 0 ) ? floor( x + 0.5 ) : ceil( x - 0.5 );
}
