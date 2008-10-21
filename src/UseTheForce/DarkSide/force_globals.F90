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


! Fortran interface to C entry plug.

module force_globals
  use definitions
#ifdef IS_MPI
  use mpiSimulation
#endif

  implicit none
  PRIVATE
#define __FORTRAN90
#include "UseTheForce/DarkSide/fInteractionMap.h"

  logical, save :: force_globals_initialized = .false.

#ifdef IS_MPI
  real( kind = dp ), allocatable, dimension(:,:), public :: q_Row
  real( kind = dp ), allocatable, dimension(:,:), public :: q_Col
  real( kind = dp ), allocatable, dimension(:,:), public :: q_group_Row
  real( kind = dp ), allocatable, dimension(:,:), public :: q_group_Col
  real( kind = dp ), allocatable, dimension(:,:), public :: eFrame_Row
  real( kind = dp ), allocatable, dimension(:,:), public :: eFrame_Col
  real( kind = dp ), allocatable, dimension(:,:), public :: A_Row
  real( kind = dp ), allocatable, dimension(:,:), public :: A_Col

  real( kind = dp ), allocatable, dimension(:,:), public :: pot_Row
  real( kind = dp ), allocatable, dimension(:,:), public :: pot_Col
  real( kind = dp ), allocatable, dimension(:,:), public :: pot_Temp
  real( kind = dp ), allocatable, dimension(:,:), public :: f_Row
  real( kind = dp ), allocatable, dimension(:,:), public :: f_Col
  real( kind = dp ), allocatable, dimension(:,:), public :: f_Temp
  real( kind = dp ), allocatable, dimension(:,:), public :: t_Row
  real( kind = dp ), allocatable, dimension(:,:), public :: t_Col
  real( kind = dp ), allocatable, dimension(:,:), public :: t_Temp
  real( kind = dp ), allocatable, dimension(:,:), public :: rf_Row
  real( kind = dp ), allocatable, dimension(:,:), public :: rf_Col
  real( kind = dp ), allocatable, dimension(:,:), public :: rf_Temp
  real( kind = dp ), allocatable, dimension(:),   public :: ppot_Row
  real( kind = dp ), allocatable, dimension(:),   public :: ppot_Col
  real( kind = dp ), allocatable, dimension(:),   public :: ppot_Temp

  integer, allocatable, dimension(:), public :: atid_Row
  integer, allocatable, dimension(:), public :: atid_Col
#endif

  integer, allocatable, dimension(:), public :: atid

  real( kind = dp ), allocatable, dimension(:,:), public :: rf
  real(kind = dp), dimension(9), public :: tau_Temp = 0.0_dp
  real(kind = dp), public :: virial_Temp = 0.0_dp

  public :: InitializeForceGlobals

contains

  subroutine InitializeForceGlobals(nlocal, thisStat)
    integer, intent(out) :: thisStat
    integer :: nAtomsInRow, nAtomsInCol
    integer :: nGroupsInRow, nGroupsInCol
    integer :: nlocal
    integer :: ndim = 3
    integer :: alloc_stat

    thisStat = 0

#ifdef IS_MPI
    nAtomsInRow = getNatomsInRow(plan_atom_row)
    nAtomsInCol = getNatomsInCol(plan_atom_col)
    nGroupsInRow = getNgroupsInRow(plan_group_row)
    nGroupsInCol = getNgroupsInCol(plan_group_col)

#endif

    call FreeForceGlobals()

#ifdef IS_MPI

    allocate(q_Row(ndim,nAtomsInRow),stat=alloc_stat)
    if (alloc_stat /= 0 ) then
       thisStat = -1
       return
    endif

    allocate(q_Col(ndim,nAtomsInCol),stat=alloc_stat)
    if (alloc_stat /= 0 ) then
       thisStat = -1
       return
    endif

    allocate(q_group_Row(ndim,nGroupsInRow),stat=alloc_stat)
    if (alloc_stat /= 0 ) then
       thisStat = -1
       return
    endif

    allocate(q_group_Col(ndim,nGroupsInCol),stat=alloc_stat)
    if (alloc_stat /= 0 ) then
       thisStat = -1
       return
    endif

    allocate(eFrame_Row(9,nAtomsInRow),stat=alloc_stat)
    if (alloc_stat /= 0 ) then
       thisStat = -1
       return
    endif

    allocate(eFrame_Col(9,nAtomsInCol),stat=alloc_stat)
    if (alloc_stat /= 0 ) then
       thisStat = -1
       return
    endif

    allocate(A_row(9,nAtomsInRow),stat=alloc_stat)
    if (alloc_stat /= 0 ) then
       thisStat = -1
       return
    endif

    allocate(A_Col(9,nAtomsInCol),stat=alloc_stat)
    if (alloc_stat /= 0 ) then
       thisStat = -1
       return
    endif

    allocate(pot_row(LR_POT_TYPES,nAtomsInRow),stat=alloc_stat)
    if (alloc_stat /= 0 ) then
       thisStat = -1
       return
    endif

    allocate(pot_Col(LR_POT_TYPES,nAtomsInCol),stat=alloc_stat)
    if (alloc_stat /= 0 ) then
       thisStat = -1
       return
    endif

    allocate(pot_Temp(LR_POT_TYPES,nlocal),stat=alloc_stat)
    if (alloc_stat /= 0 ) then
       thisStat = -1
       return
    endif

    allocate(f_Row(ndim,nAtomsInRow),stat=alloc_stat)
    if (alloc_stat /= 0 ) then
       thisStat = -1
       return
    endif

    allocate(f_Col(ndim,nAtomsInCol),stat=alloc_stat)
    if (alloc_stat /= 0 ) then
       thisStat = -1
       return
    endif

    allocate(f_Temp(ndim,nlocal),stat=alloc_stat)
    if (alloc_stat /= 0 ) then
       thisStat = -1
       return
    endif

    allocate(t_Row(ndim,nAtomsInRow),stat=alloc_stat)
    if (alloc_stat /= 0 ) then
       thisStat = -1
       return
    endif

    allocate(t_Col(ndim,nAtomsInCol),stat=alloc_stat)
    if (alloc_stat /= 0 ) then
       thisStat = -1
       return
    endif

    allocate(t_temp(ndim,nlocal),stat=alloc_stat)
    if (alloc_stat /= 0 ) then
       thisStat = -1
       return
    endif

    allocate(atid(nlocal),stat=alloc_stat)
    if (alloc_stat /= 0 ) then
       thisStat = -1
       return
    endif


    allocate(atid_Row(nAtomsInRow),stat=alloc_stat)
    if (alloc_stat /= 0 ) then
       thisStat = -1
       return
    endif

    allocate(atid_Col(nAtomsInCol),stat=alloc_stat)
    if (alloc_stat /= 0 ) then
       thisStat = -1
       return
    endif

    allocate(rf_Row(ndim,nAtomsInRow),stat=alloc_stat)
    if (alloc_stat /= 0 ) then
       thisStat = -1
       return
    endif

    allocate(rf_Col(ndim,nAtomsInCol),stat=alloc_stat)
    if (alloc_stat /= 0 ) then
       thisStat = -1
       return
    endif

    allocate(rf_Temp(ndim,nlocal),stat=alloc_stat)
    if (alloc_stat /= 0 ) then
       thisStat = -1
       return
    endif

    allocate(ppot_Row(ndim,nAtomsInRow),stat=alloc_stat)
    if (alloc_stat /= 0 ) then
       thisStat = -1
       return
    endif

    allocate(ppot_Col(ndim,nAtomsInCol),stat=alloc_stat)
    if (alloc_stat /= 0 ) then
       thisStat = -1
       return
    endif

    allocate(ppot_Temp(nlocal),stat=alloc_stat)
    if (alloc_stat /= 0 ) then
       thisStat = -1
       return
    endif

#else

    allocate(atid(nlocal),stat=alloc_stat)
    if (alloc_stat /= 0 ) then
       thisStat = -1
       return
    end if

#endif

    allocate(rf(ndim,nlocal),stat=alloc_stat)
    if (alloc_stat /= 0 ) then
       thisStat = -1
       return
    endif

    force_globals_initialized = .true.

  end subroutine InitializeForceGlobals

  subroutine FreeForceGlobals()

    !We free in the opposite order in which we allocate in.

    if (allocated(rf))         deallocate(rf)
#ifdef IS_MPI
    if (allocated(ppot_Temp))  deallocate(ppot_Temp)
    if (allocated(ppot_Col))   deallocate(ppot_Col)
    if (allocated(ppot_Row))   deallocate(ppot_Row)    
    if (allocated(rf_Temp))    deallocate(rf_Temp)
    if (allocated(rf_Col))     deallocate(rf_Col)
    if (allocated(rf_Row))     deallocate(rf_Row)    
    if (allocated(atid_Col))   deallocate(atid_Col)
    if (allocated(atid_Row))   deallocate(atid_Row)
    if (allocated(atid))       deallocate(atid)
    if (allocated(t_Temp))     deallocate(t_Temp)
    if (allocated(t_Col))      deallocate(t_Col)
    if (allocated(t_Row))      deallocate(t_Row)
    if (allocated(f_Temp))     deallocate(f_Temp)
    if (allocated(f_Col))      deallocate(f_Col)
    if (allocated(f_Row))      deallocate(f_Row)
    if (allocated(pot_Temp))   deallocate(pot_Temp)
    if (allocated(pot_Col))    deallocate(pot_Col)
    if (allocated(pot_Row))    deallocate(pot_Row)
    if (allocated(A_Col))      deallocate(A_Col)
    if (allocated(A_Row))      deallocate(A_Row)
    if (allocated(eFrame_Col))  deallocate(eFrame_Col)
    if (allocated(eFrame_Row))  deallocate(eFrame_Row)
    if (allocated(q_group_Col)) deallocate(q_group_Col)
    if (allocated(q_group_Row)) deallocate(q_group_Row)    
    if (allocated(q_Col))       deallocate(q_Col)
    if (allocated(q_Row))       deallocate(q_Row)    
#else    
    if (allocated(atid))       deallocate(atid)    
#endif

  end subroutine FreeForceGlobals

end module force_globals
