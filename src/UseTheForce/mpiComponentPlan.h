#ifdef __OOPSE_C
#ifndef USETHEFORCE_MPICOMPONENTPLAN_H
#define USETHEFORCE_MPICOMPONENTPLAN_H


/** This header provides dual access for mpiComponentPlan 
    structure in fortran and in c, C++. 
*/

typedef struct{
  int nMolGlobal;
  int nAtomsGlobal;
  int nGroupsGlobal;
  int nMolLocal;
  int nAtomsLocal;
  int nGroupsLocal;
  int myNode;
  int nProcessors;
  int rowComm;
  int columnComm;
  int nRows;
  int nColumns;
  int nAtomsInRow;
  int nAtomsInColumn;
  int nGroupsInRow;
  int nGroupsInColumn;
  int rowIndex;
  int columnIndex;
} mpiSimData;

#endif /* __MPICOMPONENTPLAN_H__*/

#endif /* __OOPSE_C*/


#ifdef __FORTRAN90
  type, public :: mpiComponentPlan
    sequence
    integer :: nMolGlobal = 0
    integer :: nAtomsGlobal = 0
    integer :: nGroupsGlobal = 0
    integer :: nMolLocal = 0
    integer :: nAtomsLocal = 0
    integer :: nGroupsLocal = 0
    integer :: myNode = 0
    integer :: nProcessors = 0
    integer :: rowComm = 0
    integer :: columnComm = 0
    integer :: nRows = 0
    integer :: nColumns = 0
    integer :: nAtomsInRow = 0
    integer :: nAtomsInColumn = 0
    integer :: nGroupsInRow = 0
    integer :: nGroupsInColumn = 0
    integer :: rowIndex = 0
    integer :: columnIndex = 0
  end type mpiComponentPlan

#endif
