#ifndef __BASS_INTERFACE_H__
#define __BASS_INTERFACE_H__



typedef enum { MOLECULE, ATOM, BOND, BEND, TORSION, COMPONENT, 
	       POSITION, ASSIGNMENT, MEMBERS, CONSTRAINT, ORIENTATION,
	       ZCONSTRAINT, RIGIDBODY, CUTOFFGROUP, BLOCK_END } event_enum;


typedef struct{
  double x;
  double y;
  double z;
} position_event;

typedef struct{
  double phi;
  double theta;
  double psi;
} orientation_event;

typedef enum { STRING, INT, DOUBLE } interface_assign_type;

typedef struct{
  interface_assign_type asmt_type;
  char lhs[80];
  union{
    int ival;
    double dval;
    char sval[120];
  }rhs;
} assignment_event;

typedef struct{
  int nMembers;
  int *memberList;
} members_event;

typedef struct{
  event_enum event_type;
  char* err_msg;

  union{
    int               blk_index; // block index
    position_event    pos;
    orientation_event ornt; // use the same structure for orientation
    assignment_event  asmt;
    members_event     mbrs;
    double            cnstr; // the constraint value
  } evt;
} event;

#ifdef __cplusplus
extern "C" {
#endif

  int event_handler( event* the_event );

#ifdef __cplusplus
}
#endif


#endif // ifndef __BASS_INTERFACE_H__
