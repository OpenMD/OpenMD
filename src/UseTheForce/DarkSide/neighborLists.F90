!!
!! Copyright (c) 2005 The University of Notre Dame. All Rights Reserved.
!!
!! The University of Notre Dame grants you ("Licensee") a
!! non-exclusive, royalty free, license to use, modify and
!! redistribute this software in source and binary code form, provided
!! that the following conditions are met:
!!
!! 1. Acknowledgement of the program authors must be made in any
!!    publication of scientific results based in part on use of the
!!    program.  An acceptable form of acknowledgement is citation of
!!    the article in which the program was described (Matthew
!!    A. Meineke, Charles F. Vardeman II, Teng Lin, Christopher
!!    J. Fennell and J. Daniel Gezelter, "OOPSE: An Object-Oriented
!!    Parallel Simulation Engine for Molecular Dynamics,"
!!    J. Comput. Chem. 26, pp. 252-271 (2005))
!!
!! 2. Redistributions of source code must retain the above copyright
!!    notice, this list of conditions and the following disclaimer.
!!
!! 3. Redistributions in binary form must reproduce the above copyright
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


!! Module neighborLists
!! Implements verlet neighbor lists for force modules.
!! Automagically expands neighbor list if size too small
!! up to maxAllocations times. If after maxAllocations we try to
!! expand the neighbor list, we get an error message and quit.
!! @author Charles F. Vardeman II
!! @author Matthew Meineke
!! @author J. Daniel Gezelter
!! @version $Id: neighborLists.F90,v 1.2 2005-01-12 22:40:45 gezelter Exp $, 

module neighborLists
  
  use definitions
#ifdef IS_MPI
  use mpiSimulation
#endif
  
  implicit none
  PRIVATE
  
  !--------------MODULE VARIABLES---------------------->
  !! Parameter for size > # of long range particles neighbor list
  !! should be.
  integer, parameter :: listMultiplier = 80
  !! Maximum number of times we should reallocate neighbor list.
  integer, parameter :: maxAllocations = 5
  !! Number of times we have allocated the neighbor list.
  integer, save       :: nAllocations = 0
  !! Pointer array to location in list for atom i.
  integer, dimension(:),public, pointer :: point => null()
  !! Neighbor list for atom i.
  integer, dimension(:),public, pointer :: list  => null()
  !! Position array of previous positions for check. Allocated first time
  !! into saveNeighborList.
  real( kind = dp ), dimension(:,:), allocatable, save  :: q0
  !! Current list size
  integer, save :: listSize
  !--------------MODULE ACCESS-------------------------->
  public :: expandNeighborList
  public :: checkNeighborList
  public :: saveNeighborList
  public :: getNeighborListSize
  
contains


  subroutine expandNeighborList(nGroups, error)
    integer, intent(out) :: error
    integer :: nGroups
    integer :: alloc_error
    integer :: oldSize = 0
    integer :: newSize = 0
    integer :: i
    integer, dimension(:), pointer :: newList => null()
    error = 0   

    !! First time through we should allocate point and list.
    !! If one is associated and one is not, something is wrong
    !! and return a error.

#ifdef IS_MPI !! // MPI
    if (.not. associated(point) .and. &
         .not. associated(list) ) then
       allocate(point(getNgroupsInRow(plan_group_row)+1),stat=alloc_error)
       if (alloc_error /= 0) then
          write(default_error,*) &
               "ExpandNeighborLists: Error in allocating MPI point"
           error = -1
          return
       end if
       allocate(list(listMultiplier * getNgroupsInCol(plan_group_col)),stat=alloc_error)
       if (alloc_error /= 0) then
          write(default_error,*) &
               "ExpandNeighborLists: Error in allocating MPI list"
          error = -1
          return
       end if
       listSize = size(list)
       nAllocations = nAllocations + 1
       return
    end if
#else !! // NONMPI
    if (.not. associated(point) .and. &
         .not. associated(list) ) then
       allocate(point(nGroups),stat=alloc_error)
       if (alloc_error /= 0) then
          write(default_error,*) &
               "ExpandNeighborLists: Error in allocating point"
          error = -1
          return
       end if
       allocate(list(listMultiplier * nGroups),stat=alloc_error)
       if (alloc_error /= 0) then
          write(default_error,*) &
               "ExpandNeighborLists: Error in allocating list"
          error = -1
          return
       end if
       listSize = size(list)
       nAllocations = nAllocations + 1
       return
    end if
#endif !! //MPI
    
    ! Expand the neighbor list
    
    ! Check to see if we have exceeded the maximum number of allocations.
    if (nAllocations > maxAllocations) then
       write(default_error,*) & 
            "ExpandNeighborList: exceeded maximum number of re-allocations"
       error = -1
       return
    else !! Expand the list.
       oldSize = size(list)


#ifdef IS_MPI !! MPI
       newSize = listMultiplier * getNgroupsInCol(plan_group_col) + oldSize
#else
       newSize = listMultiplier * nGroups + oldSize
#endif !! MPI

       

       allocate(newList(newSize), stat=alloc_error)
       if (alloc_error /= 0) then
          write(*,*) "Error allocating new neighborlist"
          error = -1
          return
       end if

       !! Copy old list to new list
       do i = 1, oldSize
          newList(i) = list(i)
       end do
       !! Free old list
       if(associated(list)) deallocate(list,stat=alloc_error)
       if (alloc_error /= 0) then
          error = -1
          return
       end if
       
       !! Point list at new list
       
       list => newList
    end if

    nAllocations = nAllocations + 1
    listSize = size(list)
  end subroutine expandNeighborList
  
  !! checks to see if any long range particle has moved
  !! through the neighbor list skin thickness.
  subroutine checkNeighborList(nGroups, q, listSkin, update_nlist)
    integer :: nGroups
    real(kind = dp), intent(in) :: listSkin
    real( kind = dp ), dimension(:,:)  :: q 
    integer :: i
    real( kind = DP ) :: dispmx 
    logical, intent(out) :: update_nlist
    real( kind = DP ) :: dispmx_tmp

    dispmx = 0.0E0_DP
    !! calculate the largest displacement of any atom in any direction
    
    !! If we have changed the particle idents, then we need to update    
    if (.not. allocated(q0)) then       
       update_nlist = .true.
       return
    end if

    if (size(q0,2) /= nGroups) then
       update_nlist = .true.
       return
    end if


#ifdef MPI 
    
    dispmx_tmp = 0.0E0_DP
    do i = 1, nGroups
       dispmx_tmp = max( abs ( q(1,i) - q0(1,i) ), dispmx_tmp )
       dispmx_tmp = max( abs ( q(2,i) - q0(2,i) ), dispmx_tmp )
       dispmx_tmp = max( abs ( q(3,i) - q0(3,i) ), dispmx_tmp )
    end do
    call mpi_allreduce(dispmx_tmp,dispmx,1,mpi_double_precision, & 
         mpi_max,mpi_comm_world,mpi_err)

#else
    
    dispmx = 0.0_DP
    do i = 1, nGroups
       dispmx = max( abs ( q(1,i) - q0(1,i) ), dispmx )
       dispmx = max( abs ( q(2,i) - q0(2,i) ), dispmx )
       dispmx = max( abs ( q(3,i) - q0(3,i) ), dispmx )
    end do

#endif   

    !! a conservative test of list skin crossings
    dispmx = 2.0E0_DP * sqrt (3.0E0_DP * dispmx * dispmx)   

    update_nlist = (dispmx.gt.listSkin)
 
   end subroutine checkNeighborList
  
  
  !! Saves neighbor list for comparison in check.
  !! Save_neighborList will work even if the number of
  !! local atoms has changed.
  subroutine saveNeighborList(nGroups, q)

    integer :: nGroups
    real(kind = dp ), dimension(3,nGroups), intent(in)  :: q
    integer :: list_size
    
    !! get size of list
    list_size = nGroups
    
    if (.not. allocated(q0)) then
       allocate(q0(3,list_size))
    else if( list_size > size(q0,2)) then
       deallocate(q0)
       allocate(q0(3,list_size))
    endif
    q0 = q
  end subroutine saveNeighborList
  
  
  function getNeighborListSize() result(returnListSize)
    integer :: returnListSize
    returnListSize = listSize
  end function getNeighborListSize
    
end module neighborLists
