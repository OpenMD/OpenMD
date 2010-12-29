!!
!! Copyright (c) 2005, 2009 The University of Notre Dame. All Rights Reserved.
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


  !! Arrays for MPI storage with metallic potentials
  real( kind = dp),save, dimension(:), allocatable, public :: dfrhodrho_col
  real( kind = dp),save, dimension(:), allocatable, public :: dfrhodrho_row
  real( kind = dp),save, dimension(:), allocatable, public :: frho_row
  real( kind = dp),save, dimension(:), allocatable, public :: frho_col
  real( kind = dp),save, dimension(:), allocatable, public :: rho_row
  real( kind = dp),save, dimension(:), allocatable, public :: rho_col
  real( kind = dp),save, dimension(:), allocatable, public :: rho_tmp

  integer, allocatable, dimension(:), public :: c_idents_Row
  integer, allocatable, dimension(:), public :: c_idents_Col
#endif
  integer, allocatable, dimension(:), public :: c_idents_local

  real( kind = dp ), allocatable, dimension(:,:), public :: rf
  real(kind = dp), dimension(9), public :: tau_Temp = 0.0_dp
  real(kind = dp), public :: virial_Temp = 0.0_dp

!! Metal potentials
  real( kind = dp), dimension(:), allocatable, public :: frho
  real( kind = dp), dimension(:), allocatable, public :: rho
  real( kind = dp), dimension(:), allocatable, public :: dfrhodrho

  
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

    allocate(c_idents_local(nlocal),stat=alloc_stat)
    if (alloc_stat /= 0 ) then
       thisStat = -1
       return
    endif

    allocate(c_idents_Row(nAtomsInRow),stat=alloc_stat)
    if (alloc_stat /= 0 ) then
       thisStat = -1
       return
    endif

    allocate(c_idents_Col(nAtomsInCol),stat=alloc_stat)
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

    allocate(ppot_Row(nAtomsInRow),stat=alloc_stat)
    if (alloc_stat /= 0 ) then
       thisStat = -1
       return
    endif

    allocate(ppot_Col(nAtomsInCol),stat=alloc_stat)
    if (alloc_stat /= 0 ) then
       thisStat = -1
       return
    endif

    allocate(ppot_Temp(nlocal),stat=alloc_stat)
    if (alloc_stat /= 0 ) then
       thisStat = -1
       return
    endif

!! Array allocation for metallic potentials
    allocate(rho_tmp(nlocal),stat=alloc_stat)
    if (alloc_stat /= 0) then 
       thisStat = -1
       return
    end if

    allocate(frho_row(nAtomsInRow),stat=alloc_stat)
    if (alloc_stat /= 0) then 
       thisStat = -1
       return
    end if
    allocate(rho_row(nAtomsInRow),stat=alloc_stat)
    if (alloc_stat /= 0) then 
       thisStat = -1
       return
    end if
    allocate(dfrhodrho_row(nAtomsInRow),stat=alloc_stat)
    if (alloc_stat /= 0) then 
       thisStat = -1
       return
    end if

    ! Now do column arrays

    allocate(frho_col(nAtomsInCol),stat=alloc_stat)
    if (alloc_stat /= 0) then 
       thisStat = -1
       return
    end if
    allocate(rho_col(nAtomsInCol),stat=alloc_stat)
    if (alloc_stat /= 0) then 
       thisStat = -1
       return
    end if
    allocate(dfrhodrho_col(nAtomsInCol),stat=alloc_stat)
    if (alloc_stat /= 0) then 
       thisStat = -1
       return
    end if
    
#else

    allocate(c_idents_local(nlocal),stat=alloc_stat)
    if (alloc_stat /= 0 ) then
       thisStat = -1
       return
    end if
#endif

    allocate(frho(nlocal),stat=alloc_stat)
    if (alloc_stat /= 0 ) then
       thisStat = -1
       return
    end if
    allocate(rho(nlocal),stat=alloc_stat)
    if (alloc_stat /= 0 ) then
       thisStat = -1
       return
    end if
    allocate(dfrhodrho(nlocal),stat=alloc_stat)
    if (alloc_stat /= 0 ) then
       thisStat = -1
       return
    end if


    
    allocate(rf(ndim,nlocal),stat=alloc_stat)
    if (alloc_stat /= 0 ) then
       thisStat = -1
       return
    endif

    force_globals_initialized = .true.

  end subroutine InitializeForceGlobals

  subroutine FreeForceGlobals()

    !We free in the opposite order in which we allocate in.

    if (allocated(rf))           deallocate(rf)
#ifdef IS_MPI
    if (allocated(ppot_Temp))     deallocate(ppot_Temp)
    if (allocated(ppot_Col))      deallocate(ppot_Col)
    if (allocated(ppot_Row))      deallocate(ppot_Row)    
    if (allocated(rf_Temp))       deallocate(rf_Temp)
    if (allocated(rf_Col))        deallocate(rf_Col)
    if (allocated(rf_Row))        deallocate(rf_Row)    

    if (allocated(c_idents_Col))  deallocate(c_idents_Col)
    if (allocated(c_idents_Row))  deallocate(c_idents_Row)
    if (allocated(c_idents_local)) deallocate(c_idents_local)
    if (allocated(t_Temp))        deallocate(t_Temp)
    if (allocated(t_Col))         deallocate(t_Col)
    if (allocated(t_Row))         deallocate(t_Row)
    if (allocated(f_Temp))        deallocate(f_Temp)
    if (allocated(f_Col))         deallocate(f_Col)
    if (allocated(f_Row))         deallocate(f_Row)
    if (allocated(pot_Temp))      deallocate(pot_Temp)
    if (allocated(pot_Col))       deallocate(pot_Col)
    if (allocated(pot_Row))       deallocate(pot_Row)
    if (allocated(A_Col))         deallocate(A_Col)
    if (allocated(A_Row))         deallocate(A_Row)
    if (allocated(eFrame_Col))    deallocate(eFrame_Col)
    if (allocated(eFrame_Row))    deallocate(eFrame_Row)
    if (allocated(q_group_Col))   deallocate(q_group_Col)
    if (allocated(q_group_Row))   deallocate(q_group_Row)    
    if (allocated(q_Col))         deallocate(q_Col)
    if (allocated(q_Row))         deallocate(q_Row)    

    if (allocated(rho_tmp))       deallocate(rho_tmp)    
    if (allocated(frho_row))      deallocate(frho_row)    
    if (allocated(rho_row))       deallocate(rho_row)    
    if (allocated(dfrhodrho_row)) deallocate(dfrhodrho_row)    
    if (allocated(frho_col))      deallocate(frho_col)    
    if (allocated(rho_col))       deallocate(rho_col)    
    if (allocated(dfrhodrho_col)) deallocate(dfrhodrho_col)    

    
#else    
    if (allocated(c_idents_local))   deallocate(c_idents_local)    
    if (allocated(rho))        deallocate(rho)    
    if (allocated(frho))       deallocate(frho)    
    if (allocated(dfrhodrho))  deallocate(dfrhodrho)    

#endif

  end subroutine FreeForceGlobals

end module force_globals
