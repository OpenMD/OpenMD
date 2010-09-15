!!
!! Copyright (c) 2005 The University of Notre Dame. All Rights Reserved.
!!
!! The University of Notre Dame grants you ("Licensee") a
!! non-exclusive, royalty free, license to use, modify and
!! redistribute this software in source and binary code form, provided
!! that the following conditions are met:
!!
!! 1. Redistributions of source code must retain the above copyright
!!    notice, this list of conditions and the following disclaimer.
!!
!! 2. Redistributions in binary form must reproduce the above copyright
!!    notice, this list of conditions and the following disclaimer in the
!!    documentation and/or other materials provided with the
!!    distribution.
!!
!! This software is provided "AS IS," without a warranty of any
!! kind. All express or implied conditions, representations and
!! warranties, including any implied warranty of merchantability,
!! fitness for a particular purpose or non-infringement, are hereby
!! excluded.  The University of Notre Dame and its licensors shall not
!! be liable for any damages suffered by licensee as a result of
!! using, modifying or distributing the software or its
!! derivatives. In no event will the University of Notre Dame or its
!! licensors be liable for any lost revenue, profit or data, or for
!! direct, indirect, special, consequential, incidental or punitive
!! damages, however caused and regardless of the theory of liability,
!! arising out of the use of or inability to use software, even if the
!! University of Notre Dame has been advised of the possibility of
!! such damages.
!!
!! SUPPORT OPEN SCIENCE!  If you use OpenMD or its source code in your
!! research, please cite the appropriate papers when you publish your
!! work.  Good starting points are:
!!                                                                      
!! [1]  Meineke, et al., J. Comp. Chem. 26, 252-271 (2005).             
!! [2]  Fennell & Gezelter, J. Chem. Phys. 124, 234104 (2006).          
!! [3]  Sun, Lin & Gezelter, J. Chem. Phys. 128, 24107 (2008).          
!! [4]  Vardeman & Gezelter, in progress (2009).
!!


!! MPI support for long range forces using force decomposition 
!! on a square grid of processors.
!! Corresponds to mpiSimulation.cpp for C++
!! mpi_module also contains a private interface for mpi f90 routines.
!!
!! @author Charles F. Vardeman II
!! @author Matthew Meineke
!! @version $Id$, $Date$, $Name: not supported by cvs2svn $, $Revision$

module mpiSimulation  
  use definitions
  use status
#ifdef IS_MPI
  use OpenMDMPI
  implicit none
  PRIVATE
#endif


  !! Include mpiComponentPlan type. mpiComponentPlan is a 
  !! dual header file for both c and fortran.
#define __FORTRAN90
#include "UseTheForce/mpiComponentPlan.h"
public :: setupSimParallel

#ifdef IS_MPI

  !! PUBLIC  Subroutines contained in this module
  !! gather and scatter are a generic interface
  !! to gather and scatter routines
  public :: gather, scatter
  public :: replanSimParallel
  public :: getNatomsInCol
  public :: getNatomsInRow
  public :: getNgroupsInCol
  public :: getNgroupsInRow
  public :: isMPISimSet
  public :: printComponentPlan
  public :: getMyNode

  !! PUBLIC  Subroutines contained in MPI module
  public :: mpi_bcast
  public :: mpi_allreduce
  !  public :: mpi_reduce
  public :: mpi_send
  public :: mpi_recv
  public :: mpi_get_processor_name
  public :: mpi_finalize

  !! PUBLIC mpi variables
  public :: mpi_comm_world
  public :: mpi_character
  public :: mpi_integer
  public :: mpi_lor
  public :: mpi_logical
  public :: mpi_real
  public :: mpi_double_precision
  public :: mpi_sum
  public :: mpi_max
  public :: mpi_status_size
  public :: mpi_any_source



  !! Safety logical to prevent access to ComponetPlan until
  !! set by C++.
  logical, save :: ComponentPlanSet = .false.


  !! generic mpi error declaration.
  integer, public :: mpi_err
  character(len = statusMsgSize) :: errMsg

#ifdef PROFILE
  public :: printCommTime
  public :: getCommTime
  real,save   :: commTime = 0.0
  real   :: commTimeInitial,commTimeFinal
#endif


  !! Tags used during force loop for parallel simulation
  integer, public, allocatable, dimension(:) :: AtomLocalToGlobal
  integer, public, allocatable, dimension(:) :: AtomRowToGlobal
  integer, public, allocatable, dimension(:) :: AtomColToGlobal
  integer, public, allocatable, dimension(:) :: AtomGlobalToLocal
  integer, public, allocatable, dimension(:) :: GroupLocalToGlobal
  integer, public, allocatable, dimension(:) :: GroupRowToGlobal
  integer, public, allocatable, dimension(:) :: GroupColToGlobal

  !! Logical set true if mpiSimulation has been initialized
  logical, save :: isSimSet = .false.


  type (mpiComponentPlan), save :: mpiSim

  !! gs_plan contains plans for gather and scatter routines
  type, public :: gs_plan
     private
     type (mpiComponentPlan), pointer :: gsComponentPlan => NULL()
     integer :: gsPlanSize !! size of this plan (nDim*nComponents)
     integer :: globalPlanSize !! size of all components in plan
     integer, dimension(:), pointer :: displs !! Displacements array for mpi indexed from 0.
     integer, dimension(:), pointer :: counts !! Counts array for mpi indexed from 0.
     integer :: myPlanComm  !! My communicator for this plan
     integer :: myPlanRank  !! My rank in this plan
     integer :: planNprocs  !! Number of processors in this plan
  end type gs_plan

  ! plans for different decompositions
  type (gs_plan), public, save :: plan_atom_row
  type (gs_plan), public, save :: plan_atom_row_3d
  type (gs_plan), public, save :: plan_atom_col
  type (gs_plan), public, save :: plan_atom_col_3d
  type (gs_plan),  public, save :: plan_atom_row_Rotation
  type (gs_plan),  public, save :: plan_atom_col_Rotation
  type (gs_plan),  public, save :: plan_group_row
  type (gs_plan),  public, save :: plan_group_col
  type (gs_plan),  public, save :: plan_group_row_3d
  type (gs_plan),  public, save :: plan_group_col_3d

  type (mpiComponentPlan), pointer :: simComponentPlan

  ! interface for gather scatter routines
  !! Generic interface for gather.
  !! Gathers an local array into row or column array
  !! Interface provided for integer, double and double
  !! rank 2 arrays.
  interface gather
     module procedure gather_integer
     module procedure gather_double
     module procedure gather_double_2d
  end interface

  !! Generic interface for scatter.
  !! Scatters a row or column array, adding componets
  !! and reducing them to a local nComponent array.
  !! Interface provided for double and double rank=2 arrays.

  interface scatter
     module procedure scatter_double
     module procedure scatter_double_2d
  end interface



contains

  !! Sets up mpiComponentPlan with structure passed from C++.
  subroutine setupSimParallel(thisComponentPlan, nAtomTags, atomTags, &
       nGroupTags, groupTags, status)
    !! Passed Arguments
    !! mpiComponentPlan struct from C
    type (mpiComponentPlan), intent(inout) :: thisComponentPlan
    !! Number of tags passed
    integer, intent(in) :: nAtomTags, nGroupTags
    !! Result status, 0 = normal, -1 = error
    integer, intent(out) :: status
    integer :: localStatus
    !! Global reference tag for local particles
    integer, dimension(nAtomTags), intent(inout) :: atomTags
    integer, dimension(nGroupTags), intent(inout) :: groupTags

    !write(*,*) 'mpiSim_mod thinks node', thisComponentPlan%myNode, &
    !     ' has atomTags(1) = ', atomTags(1)

    status = 0
    if (componentPlanSet) then
       return
    endif
    componentPlanSet = .true.

    !! copy c component plan to fortran   
    mpiSim = thisComponentPlan 

    call make_Force_Grid(mpiSim, localStatus)
    if (localStatus /= 0) then
       write(errMsg, *) 'An error in making the force grid has occurred'
       call handleError("setupSimParallel", errMsg)
       status = -1
       return
    endif

    call updateGridComponents(mpiSim, localStatus)
    if (localStatus /= 0) then
       write(errMsg,*) "Error updating grid components"
       call handleError("setupSimParallel", errMsg)
       status = -1
       return
    endif

    !! initialize gather and scatter plans used in this simulation
    call plan_gather_scatter(1, mpiSim%nAtomsLocal, &
         mpiSim, mpiSim%rowComm, plan_atom_row)
    call plan_gather_scatter(nDim, mpiSim%nAtomsLocal, &
         mpiSim, mpiSim%rowComm, plan_atom_row_3d)
    call plan_gather_scatter(9, mpiSim%nAtomsLocal, &
         mpiSim, mpiSim%rowComm, plan_atom_row_Rotation)
    call plan_gather_scatter(1, mpiSim%nGroupsLocal, &
         mpiSim, mpiSim%rowComm, plan_group_row)
    call plan_gather_scatter(nDim, mpiSim%nGroupsLocal, &
         mpiSim, mpiSim%rowComm, plan_group_row_3d)

    call plan_gather_scatter(1, mpiSim%nAtomsLocal, &
         mpiSim, mpiSim%columnComm, plan_atom_col)
    call plan_gather_scatter(nDim, mpiSim%nAtomsLocal, &
         mpiSim, mpiSim%columnComm, plan_atom_col_3d)
    call plan_gather_scatter(9, mpiSim%nAtomsLocal, &
         mpiSim, mpiSim%columnComm, plan_atom_col_Rotation)
    call plan_gather_scatter(1, mpiSim%nGroupsLocal, &
         mpiSim, mpiSim%columnComm, plan_group_col)
    call plan_gather_scatter(nDim, mpiSim%nGroupsLocal, &
         mpiSim, mpiSim%columnComm, plan_group_col_3d)

    !  Initialize tags    

    call setAtomTags(atomTags,localStatus)
    if (localStatus /= 0) then
       write(errMsg, *) 'An error in setting Atom Tags has occured'
       call handleError("setupSimParallel", errMsg)
       status = -1
       return
    endif


    call setGroupTags(groupTags,localStatus)
    if (localStatus /= 0) then
       write(errMsg, *) 'An error in setting Group Tags has occured'
       call handleError("setupSimParallel", errMsg)
       status = -1
       return
    endif

    isSimSet = .true.

    !    call printComponentPlan(mpiSim,0)
  end subroutine setupSimParallel

  subroutine replanSimParallel(thisComponentPlan,status)
    !  Passed Arguments
    !! mpiComponentPlan struct from C
    type (mpiComponentPlan), intent(inout) :: thisComponentPlan   
    integer, intent(out) :: status
    integer :: localStatus
    integer :: mpierror
    status = 0

    call updateGridComponents(thisComponentPlan,localStatus)
    if (localStatus /= 0) then
       status = -1
       return
    endif

    !! Unplan Gather Scatter plans
    call unplan_gather_scatter(plan_atom_row)
    call unplan_gather_scatter(plan_atom_row_3d)
    call unplan_gather_scatter(plan_atom_row_Rotation)
    call unplan_gather_scatter(plan_group_row)
    call unplan_gather_scatter(plan_group_row_3d)

    call unplan_gather_scatter(plan_atom_col)
    call unplan_gather_scatter(plan_atom_col_3d)
    call unplan_gather_scatter(plan_atom_col_Rotation)
    call unplan_gather_scatter(plan_group_col)
    call unplan_gather_scatter(plan_group_col_3d)

    !! initialize gather and scatter plans used in this simulation
    call plan_gather_scatter(1, mpiSim%nAtomsLocal, &
         mpiSim, mpiSim%rowComm, plan_atom_row)
    call plan_gather_scatter(nDim, mpiSim%nAtomsLocal, &
         mpiSim, mpiSim%rowComm, plan_atom_row_3d)
    call plan_gather_scatter(9, mpiSim%nAtomsLocal, &
         mpiSim, mpiSim%rowComm, plan_atom_row_Rotation)
    call plan_gather_scatter(1, mpiSim%nGroupsLocal, &
         mpiSim, mpiSim%rowComm, plan_group_row)
    call plan_gather_scatter(nDim, mpiSim%nGroupsLocal, &
         mpiSim, mpiSim%rowComm, plan_group_row_3d)

    call plan_gather_scatter(1, mpiSim%nAtomsLocal, &
         mpiSim, mpiSim%columnComm, plan_atom_col)
    call plan_gather_scatter(nDim, mpiSim%nAtomsLocal, &
         mpiSim, mpiSim%columnComm, plan_atom_col_3d)
    call plan_gather_scatter(9, mpiSim%nAtomsLocal, &
         mpiSim, mpiSim%columnComm, plan_atom_col_Rotation)
    call plan_gather_scatter(1, mpiSim%nGroupsLocal, &
         mpiSim, mpiSim%columnComm, plan_group_col)
    call plan_gather_scatter(nDim, mpiSim%nGroupsLocal, &
         mpiSim, mpiSim%columnComm, plan_group_col_3d)

  end subroutine replanSimParallel

  !! Updates number of row and column components for long range forces.
  subroutine updateGridComponents(thisComponentPlan, status)
    type (mpiComponentPlan) :: thisComponentPlan !! mpiComponentPlan

    !! Status return
    !! -  0 Success
    !! - -1 Failure
    integer, intent(out) :: status 
    integer :: nAtomsLocal
    integer :: nAtomsInRow = 0
    integer :: nAtomsInColumn = 0
    integer :: nGroupsLocal
    integer :: nGroupsInRow = 0
    integer :: nGroupsInColumn = 0
    integer :: mpiErrors

    status = 0
    if (.not. componentPlanSet) return

    if (thisComponentPlan%nAtomsLocal == 0) then
       status = -1
       return
    endif
    if (thisComponentPlan%nGroupsLocal == 0) then
       write(*,*) 'tcp%ngl = ', thisComponentPlan%nGroupsLocal
       status = -1
       return
    endif

    nAtomsLocal = thisComponentPlan%nAtomsLocal
    nGroupsLocal = thisComponentPlan%nGroupsLocal

    call mpi_allreduce(nAtomsLocal, nAtomsInRow, 1, mpi_integer, &
         mpi_sum, thisComponentPlan%rowComm, mpiErrors)
    if (mpiErrors /= 0) then
       status = -1
       return
    endif

    call mpi_allreduce(nAtomsLocal, nAtomsInColumn, 1, mpi_integer, &
         mpi_sum, thisComponentPlan%columnComm, mpiErrors)    
    if (mpiErrors /= 0) then
       status = -1
       return
    endif

    call mpi_allreduce(nGroupsLocal, nGroupsInRow, 1, mpi_integer, &
         mpi_sum, thisComponentPlan%rowComm, mpiErrors)
    if (mpiErrors /= 0) then
       status = -1
       return
    endif

    call mpi_allreduce(nGroupsLocal, nGroupsInColumn, 1, mpi_integer, &
         mpi_sum, thisComponentPlan%columnComm, mpiErrors)    
    if (mpiErrors /= 0) then
       status = -1
       return
    endif

    thisComponentPlan%nAtomsInRow = nAtomsInRow
    thisComponentPlan%nAtomsInColumn = nAtomsInColumn
    thisComponentPlan%nGroupsInRow = nGroupsInRow
    thisComponentPlan%nGroupsInColumn = nGroupsInColumn

  end subroutine updateGridComponents


  !! Creates a square force decomposition of processors into row and column
  !! communicators.
  subroutine make_Force_Grid(thisComponentPlan,status)
    type (mpiComponentPlan) :: thisComponentPlan
    integer, intent(out) :: status !! status returns -1 if error
    integer :: nColumnsMax !! Maximum number of columns
    integer :: nWorldProcessors !! Total number of processors in World comm.
    integer :: rowIndex !! Row for this processor.
    integer :: columnIndex !! Column for this processor.
    integer :: nRows !! Total number of rows.
    integer :: nColumns !! Total number of columns.
    integer :: mpiErrors !! MPI errors.
    integer :: rowCommunicator !! MPI row communicator.
    integer :: columnCommunicator
    integer :: myWorldRank
    integer :: i


    if (.not. ComponentPlanSet) return
    status = 0

    !! We make a dangerous assumption here that if numberProcessors is
    !! zero, then we need to get the information from MPI.
    if (thisComponentPlan%nProcessors == 0 ) then 
       call mpi_comm_size( MPI_COMM_WORLD, nWorldProcessors,mpiErrors)
       if ( mpiErrors /= 0 ) then
          status = -1 
          return
       endif
       call mpi_comm_rank( MPI_COMM_WORLD,myWorldRank,mpiErrors)
       if ( mpiErrors /= 0 ) then
          status = -1 
          return
       endif

    else
       nWorldProcessors = thisComponentPlan%nProcessors
       myWorldRank = thisComponentPlan%myNode
    endif

    nColumnsMax = nint(sqrt(real(nWorldProcessors,kind=dp)))

    do i = 1, nColumnsMax
       if (mod(nWorldProcessors,i) == 0) nColumns = i
    end do

    nRows = nWorldProcessors/nColumns

    rowIndex = myWorldRank/nColumns

    call mpi_comm_split(mpi_comm_world,rowIndex,0,rowCommunicator,mpiErrors)
    if ( mpiErrors /= 0 ) then
       write(errMsg, *) 'An error ',mpiErrors ,'occurred in splitting communicators'
       call handleError("makeForceGrid", errMsg)
       status = -1 
       return
    endif

    columnIndex = mod(myWorldRank,nColumns)
    call mpi_comm_split(mpi_comm_world,columnIndex,0,columnCommunicator,mpiErrors)
    if ( mpiErrors /= 0 ) then
       write(errMsg, *) "MPI comm split faild at columnCommunicator by error ",mpiErrors 
       call handleError("makeForceGrid", errMsg)
       status = -1 
       return
    endif

    ! Set appropriate components of thisComponentPlan
    thisComponentPlan%rowComm = rowCommunicator
    thisComponentPlan%columnComm = columnCommunicator
    thisComponentPlan%rowIndex = rowIndex
    thisComponentPlan%columnIndex = columnIndex
    thisComponentPlan%nRows = nRows
    thisComponentPlan%nColumns = nColumns

  end subroutine make_Force_Grid

  !! initalizes a gather scatter plan
  subroutine plan_gather_scatter( nDim, nObjects, thisComponentPlan, &
       thisComm, this_plan, status)  
    integer, intent(in) :: nDim !! Number of dimensions for gather scatter plan
    integer, intent(in) :: nObjects
    type (mpiComponentPlan), intent(in), target :: thisComponentPlan
    type (gs_plan), intent(out) :: this_plan !! MPI Component Plan
    integer, intent(in) :: thisComm !! MPI communicator for this plan

    integer :: arraySize !! size to allocate plan for
    integer, intent(out), optional :: status
    integer :: ierror
    integer :: i,junk

    if (present(status)) status = 0

    !! Set gsComponentPlan pointer 
    !! to the componet plan we want to use for this gather scatter plan.
    !! WARNING this could be dangerous since thisComponentPlan was origionally
    !! allocated in C++ and there is a significant difference between c and 
    !! f95 pointers....  
    this_plan%gsComponentPlan => thisComponentPlan

    ! Set this plan size for displs array.
    this_plan%gsPlanSize = nDim * nObjects

    ! Duplicate communicator for this plan
    call mpi_comm_dup(thisComm, this_plan%myPlanComm, mpi_err)
    if (mpi_err /= 0) then
       if (present(status)) status = -1
       return
    end if
    call mpi_comm_rank(this_plan%myPlanComm, this_plan%myPlanRank, mpi_err)
    if (mpi_err /= 0) then
       if (present(status)) status = -1
       return
    end if

    call mpi_comm_size(this_plan%myPlanComm, this_plan%planNprocs, mpi_err)

    if (mpi_err /= 0) then
       if (present(status)) status = -1
       return
    end if

    !! counts and displacements arrays are indexed from 0 to be compatable
    !! with MPI arrays.
    allocate (this_plan%counts(0:this_plan%planNprocs-1),STAT=ierror)
    if (ierror /= 0) then
       if (present(status)) status = -1
       return
    end if

    allocate (this_plan%displs(0:this_plan%planNprocs-1),STAT=ierror)
    if (ierror /= 0) then
       if (present(status)) status = -1
       return
    end if

    !! gather all the local sizes into a size # processors array.
    call mpi_allgather(this_plan%gsPlanSize,1,mpi_integer,this_plan%counts, &
         1,mpi_integer,thisComm,mpi_err)

    if (mpi_err /= 0) then
       if (present(status)) status = -1
       return
    end if

    !! figure out the total number of particles in this plan
    this_plan%globalPlanSize = sum(this_plan%counts)

    !! initialize plan displacements.
    this_plan%displs(0) = 0
    do i = 1, this_plan%planNprocs - 1,1
       this_plan%displs(i) = this_plan%displs(i-1) + this_plan%counts(i-1) 
    end do
  end subroutine plan_gather_scatter

  subroutine unplan_gather_scatter(this_plan)
    type (gs_plan), intent(inout) :: this_plan

    this_plan%gsComponentPlan => null()
    call mpi_comm_free(this_plan%myPlanComm,mpi_err)

    deallocate(this_plan%counts)
    deallocate(this_plan%displs)

  end subroutine unplan_gather_scatter

  subroutine gather_integer( sbuffer, rbuffer, this_plan, status)

    type (gs_plan), intent(inout) :: this_plan
    integer, dimension(:), intent(inout) :: sbuffer
    integer, dimension(:), intent(inout) :: rbuffer
    integer :: noffset
    integer, intent(out), optional :: status
    integer :: i

    if (present(status)) status = 0
    noffset = this_plan%displs(this_plan%myPlanRank)

    !    if (getmyNode() == 1) then
    !       write(*,*) "Node 0 printing allgatherv vars"
    !       write(*,*) "Noffset: ", noffset
    !       write(*,*) "PlanSize: ", this_plan%gsPlanSize
    !       write(*,*) "PlanComm: ", this_plan%myPlanComm
    !    end if

    call mpi_allgatherv(sbuffer,this_plan%gsPlanSize, mpi_integer, &
         rbuffer,this_plan%counts,this_plan%displs,mpi_integer, &
         this_plan%myPlanComm, mpi_err)

    if (mpi_err /= 0) then
       write(errMsg, *) "mpi_allgatherv failed by error message ",mpi_err 
       call handleError("gather_integer", errMsg)
       if (present(status)) status  = -1
    endif

  end subroutine gather_integer

  subroutine gather_double( sbuffer, rbuffer, this_plan, status)

    type (gs_plan), intent(in) :: this_plan
    real( kind = DP ), dimension(:), intent(inout) :: sbuffer
    real( kind = DP ), dimension(:), intent(inout) :: rbuffer
    integer :: noffset
    integer, intent(out), optional :: status


    if (present(status)) status = 0
    noffset = this_plan%displs(this_plan%myPlanRank)
#ifdef PROFILE
    call cpu_time(commTimeInitial)
#endif
#ifdef SINGLE_PRECISION
    call mpi_allgatherv(sbuffer,this_plan%gsPlanSize, mpi_real, &
         rbuffer,this_plan%counts,this_plan%displs,mpi_real, &
         this_plan%myPlanComm, mpi_err)
#else
    call mpi_allgatherv(sbuffer,this_plan%gsPlanSize, mpi_double_precision, &
         rbuffer,this_plan%counts,this_plan%displs,mpi_double_precision, &
         this_plan%myPlanComm, mpi_err)
#endif
#ifdef PROFILE
    call cpu_time(commTimeFinal)
    commTime = commTime + commTimeFinal - commTimeInitial
#endif

    if (mpi_err /= 0) then
       write(errMsg, *) "mpi_allgatherv failed by error message ",mpi_err 
       call handleError("gather_double", errMsg)
       if (present(status)) status  = -1
    endif

  end subroutine gather_double

  subroutine gather_double_2d( sbuffer, rbuffer, this_plan, status)

    type (gs_plan), intent(in) :: this_plan
    real( kind = DP ), dimension(:,:), intent(inout) :: sbuffer
    real( kind = DP ), dimension(:,:), intent(inout) :: rbuffer
    integer :: noffset,i,ierror
    integer, intent(out), optional :: status

    external  mpi_allgatherv

    if (present(status)) status = 0

    !    noffset = this_plan%displs(this_plan%me)
#ifdef PROFILE
    call cpu_time(commTimeInitial)
#endif

#ifdef SINGLE_PRECISION
    call mpi_allgatherv(sbuffer,this_plan%gsPlanSize, mpi_real, &
         rbuffer,this_plan%counts,this_plan%displs,mpi_real, &
         this_plan%myPlanComm, mpi_err)
#else
    call mpi_allgatherv(sbuffer,this_plan%gsPlanSize, mpi_double_precision, &
         rbuffer,this_plan%counts,this_plan%displs,mpi_double_precision, &
         this_plan%myPlanComm, mpi_err)
#endif

#ifdef PROFILE
    call cpu_time(commTimeFinal)
    commTime = commTime + commTimeFinal - commTimeInitial
#endif

    if (mpi_err /= 0) then
       write(errMsg, *) "mpi_allgatherv failed by error message ",mpi_err 
       call handleError("gather_double_2d", errMsg)
       if (present(status)) status = -1
    endif

  end subroutine gather_double_2d

  subroutine scatter_double( sbuffer, rbuffer, this_plan, status)

    type (gs_plan), intent(in) :: this_plan
    real( kind = DP ), dimension(:), intent(inout) :: sbuffer
    real( kind = DP ), dimension(:), intent(inout) :: rbuffer
    integer, intent(out), optional :: status
    external mpi_reduce_scatter

    if (present(status)) status = 0

#ifdef PROFILE
    call cpu_time(commTimeInitial)
#endif
#ifdef SINGLE_PRECISION
    call mpi_reduce_scatter(sbuffer,rbuffer, this_plan%counts, &
         mpi_real, MPI_SUM, this_plan%myPlanComm, mpi_err)
#else
    call mpi_reduce_scatter(sbuffer,rbuffer, this_plan%counts, &
         mpi_double_precision, MPI_SUM, this_plan%myPlanComm, mpi_err)
#endif
#ifdef PROFILE
    call cpu_time(commTimeFinal)
    commTime = commTime + commTimeFinal - commTimeInitial
#endif

    if (mpi_err /= 0) then
       write(errMsg, *) "mpi_reduce_scatter failed by error message ",mpi_err 
       call handleError("scatter_double", errMsg)
       if (present(status))  status = -1
    endif

  end subroutine scatter_double

  subroutine scatter_double_2d( sbuffer, rbuffer, this_plan, status)

    type (gs_plan), intent(in) :: this_plan
    real( kind = DP ), dimension(:,:), intent(inout) :: sbuffer
    real( kind = DP ), dimension(:,:), intent(inout) :: rbuffer
    integer, intent(out), optional :: status
    external mpi_reduce_scatter

    if (present(status)) status = 0
#ifdef PROFILE
    call cpu_time(commTimeInitial)
#endif
#ifdef SINGLE_PRECISION
    call mpi_reduce_scatter(sbuffer,rbuffer, this_plan%counts, &
         mpi_real, MPI_SUM, this_plan%myPlanComm, mpi_err)
#else
    call mpi_reduce_scatter(sbuffer,rbuffer, this_plan%counts, &
         mpi_double_precision, MPI_SUM, this_plan%myPlanComm, mpi_err)
#endif
#ifdef PROFILE
    call cpu_time(commTimeFinal)
    commTime = commTime + commTimeFinal - commTimeInitial
#endif

    if (mpi_err /= 0) then
       write(errMsg, *) "mpi_reduce_scatter failed by error message ",mpi_err 
       call handleError("scatter_double_2d", errMsg)
       if (present(status)) status = -1
    endif

  end subroutine scatter_double_2d

  subroutine setAtomTags(tags, status)
    integer, dimension(:) :: tags
    integer :: status

    integer :: alloc_stat

    integer :: nAtomsInCol
    integer :: nAtomsInRow
    integer :: i, glob, nAtomTags, maxGlobal
    

    status = 0
    ! allocate row arrays
    nAtomsInRow = getNatomsInRow(plan_atom_row)
    nAtomsInCol = getNatomsInCol(plan_atom_col)

    if(.not. allocated(AtomLocalToGlobal)) then
       allocate(AtomLocalToGlobal(size(tags)),STAT=alloc_stat)
       if (alloc_stat /= 0 ) then
          status = -1 
          return
       endif
    else
       deallocate(AtomLocalToGlobal)
       allocate(AtomLocalToGlobal(size(tags)),STAT=alloc_stat)
       if (alloc_stat /= 0 ) then
          status = -1 
          return
       endif

    endif

    AtomLocalToGlobal = tags

    if (.not. allocated(AtomRowToGlobal)) then
       allocate(AtomRowToGlobal(nAtomsInRow),STAT=alloc_stat)
       if (alloc_stat /= 0 ) then
          status = -1 
          return
       endif
    else
       deallocate(AtomRowToGlobal)
       allocate(AtomRowToGlobal(nAtomsInRow),STAT=alloc_stat)
       if (alloc_stat /= 0 ) then
          status = -1 
          return
       endif

    endif
    ! allocate column arrays
    if (.not. allocated(AtomColToGlobal)) then
       allocate(AtomColToGlobal(nAtomsInCol),STAT=alloc_stat)
       if (alloc_stat /= 0 ) then
          status = -1 
          return
       endif
    else
       deallocate(AtomColToGlobal)
       allocate(AtomColToGlobal(nAtomsInCol),STAT=alloc_stat)
       if (alloc_stat /= 0 ) then
          status = -1 
          return
       endif
    endif

    call gather(tags, AtomRowToGlobal, plan_atom_row)
    call gather(tags, AtomColToGlobal, plan_atom_col)

    nAtomTags = size(tags)
    maxGlobal = -1
    do i = 1, nAtomTags
       if (tags(i).gt.maxGlobal) maxGlobal = tags(i)
    enddo

    if(.not. allocated(AtomGlobalToLocal)) then
       allocate(AtomGlobalToLocal(maxGlobal),STAT=alloc_stat)
       if (alloc_stat /= 0 ) then
          status = -1 
          return
       endif
    else
       deallocate(AtomGlobalToLocal)
       allocate(AtomGlobalToLocal(maxGlobal),STAT=alloc_stat)
       if (alloc_stat /= 0 ) then
          status = -1 
          return
       endif

    endif

    do i = 1, nAtomTags
       glob = AtomLocalToGlobal(i)
       AtomGlobalToLocal(glob) = i
    enddo

  end subroutine setAtomTags

  subroutine setGroupTags(tags, status)
    integer, dimension(:) :: tags
    integer :: status

    integer :: alloc_stat

    integer :: nGroupsInCol
    integer :: nGroupsInRow

    status = 0

    nGroupsInRow = getNgroupsInRow(plan_group_row)
    nGroupsInCol = getNgroupsInCol(plan_group_col)

    if(allocated(GroupLocalToGlobal)) then
       deallocate(GroupLocalToGlobal)
    endif
    allocate(GroupLocalToGlobal(size(tags)),STAT=alloc_stat)
    if (alloc_stat /= 0 ) then
       status = -1 
       return
    endif

    GroupLocalToGlobal = tags

    if(allocated(GroupRowToGlobal)) then
       deallocate(GroupRowToGlobal)
    endif
    allocate(GroupRowToGlobal(nGroupsInRow),STAT=alloc_stat)
    if (alloc_stat /= 0 ) then
       status = -1 
       return
    endif

    if(allocated(GroupColToGlobal)) then
       deallocate(GroupColToGlobal)
    endif
    allocate(GroupColToGlobal(nGroupsInCol),STAT=alloc_stat)
    if (alloc_stat /= 0 ) then
       status = -1 
       return
    endif

    call gather(tags, GroupRowToGlobal, plan_group_row)
    call gather(tags, GroupColToGlobal, plan_group_col)

  end subroutine setGroupTags

  function getNatomsInCol(thisplan) result(nInCol)
    type (gs_plan), intent(in) :: thisplan
    integer :: nInCol
    nInCol = thisplan%gsComponentPlan%nAtomsInColumn
  end function getNatomsInCol

  function getNatomsInRow(thisplan) result(nInRow)
    type (gs_plan), intent(in) :: thisplan
    integer :: nInRow
    nInRow = thisplan%gsComponentPlan%nAtomsInRow
  end function getNatomsInRow

  function getNgroupsInCol(thisplan) result(nInCol)
    type (gs_plan), intent(in) :: thisplan
    integer :: nInCol
    nInCol = thisplan%gsComponentPlan%nGroupsInColumn
  end function getNgroupsInCol

  function getNgroupsInRow(thisplan) result(nInRow)
    type (gs_plan), intent(in) :: thisplan
    integer :: nInRow
    nInRow = thisplan%gsComponentPlan%nGroupsInRow
  end function getNgroupsInRow

  function isMPISimSet() result(isthisSimSet)
    logical :: isthisSimSet
    if (isSimSet) then 
       isthisSimSet = .true.
    else
       isthisSimSet = .false.
    endif
  end function isMPISimSet

  subroutine printComponentPlan(this_plan,printNode)

    type (mpiComponentPlan), intent(in) :: this_plan
    integer, optional :: printNode
    logical :: print_me = .false.

    if (present(printNode)) then
       if (printNode == mpiSim%myNode) print_me = .true.
    else
       print_me = .true.
    endif

    if (print_me) then
       write(default_error,*) "SetupSimParallel: writing component plan"

       write(default_error,*) "nMolGlobal: ", mpiSim%nMolGlobal
       write(default_error,*) "nAtomsGlobal: ", mpiSim%nAtomsGlobal
       write(default_error,*) "nAtomsLocal: ", mpiSim%nAtomsLocal
       write(default_error,*) "myNode: ", mpiSim%myNode
       write(default_error,*) "nProcessors: ", mpiSim%nProcessors
       write(default_error,*) "rowComm: ", mpiSim%rowComm
       write(default_error,*) "columnComm: ", mpiSim%columnComm
       write(default_error,*) "nRows: ", mpiSim%nRows
       write(default_error,*) "nColumns: ", mpiSim%nColumns
       write(default_error,*) "nAtomsInRow: ", mpiSim%nAtomsInRow
       write(default_error,*) "nAtomsInColumn: ", mpiSim%nAtomsInColumn
       write(default_error,*) "rowIndex: ", mpiSim%rowIndex
       write(default_error,*) "columnIndex: ", mpiSim%columnIndex
    endif
  end subroutine printComponentPlan

  function getMyNode() result(myNode)
    integer :: myNode
    myNode = mpiSim%myNode
  end function getMyNode

#ifdef PROFILE
  subroutine printCommTime()
    write(*,*) "MPI communication time is: ", commTime
  end subroutine printCommTime

  function getCommTime() result(comm_time)
    real :: comm_time
    comm_time = commTime
  end function getCommTime

#endif

#else
contains
  subroutine setupSimParallel(thisComponentPlan, nAtomTags, atomTags, &
       nGroupTags, groupTags, status)
  !! Passed Arguments
    !! mpiComponentPlan struct from C
    type (mpiComponentPlan), intent(inout) :: thisComponentPlan
    !! Number of tags passed
    integer, intent(in) :: nAtomTags, nGroupTags
    !! Result status, 0 = normal, -1 = error
    integer, intent(out) :: status
    integer :: localStatus
    !! Global reference tag for local particles
    integer, dimension(nAtomTags), intent(inout) :: atomTags
    integer, dimension(nGroupTags), intent(inout) :: groupTags
   
  end subroutine setupSimParallel

#endif


end module mpiSimulation

