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

!! module defines atypes available to simulation

module interactions

  use vector_class
  use atype_module
  use lj

  implicit none
  private

#define __FORTRAN90
#include "UseTheForce/DarkSide/fInteractionMap.h"
  
  type, public :: Interaction
     integer :: InteractionHash
     real(kind=dp) :: rCut
  end type Interaction
  
  type(Interaction), public, dimension(:,:), allocatable :: InteractionMap
  
  !public :: addInteraction
  !public :: setInteractionHash
  !public :: getInteractionHash

  public :: createInteractionMap
  
  
contains
  
  subroutine createInteractionMap(status)
    integer :: nAtypes
    integer :: status
    integer :: i
    integer :: j

    if (.not. associated(atypes)) then
       call handleError("atype", "atypes was not present before call of createDefaultInteractionMap!")
       status = -1
       return
    endif
    
    nAtypes = getSize(atypes)
    
    if (nAtypes == 0) then
       status = -1
       return
    end if

    if (.not. allocated(InteractionMap)) then
       allocate(InteractionMap(nAtypes,nAtypes))
    endif
        
    do i = 1, nAtypes
       call getElementProperty(atypes, i, "is_LennardJones", i_is_LJ)
       call getElementProperty(atypes, i, "is_Electrostatic", i_is_Elect)
       call getElementProperty(atypes, i, "is_Sticky", i_is_Sticky)
       call getElementProperty(atypes, i, "is_StickyPower", i_is_StickyP)
       call getElementProperty(atypes, i, "is_GayBerne", i_is_GB)
       call getElementProperty(atypes, i, "is_EAM", i_is_EAM)
       call getElementProperty(atypes, i, "is_Shape", i_is_Shape)

       do j = i, nAtypes

          iHash = 0
          myRcut = 0.0_dp

          call getElementProperty(atypes, j, "is_LennardJones", j_is_LJ)
          call getElementProperty(atypes, j, "is_Electrostatic", j_is_Elect)
          call getElementProperty(atypes, j, "is_Sticky", j_is_Sticky)
          call getElementProperty(atypes, j, "is_StickyPower", j_is_StickyP)
          call getElementProperty(atypes, j, "is_GayBerne", j_is_GB)
          call getElementProperty(atypes, j, "is_EAM", j_is_EAM)
          call getElementProperty(atypes, j, "is_Shape", j_is_Shape)

          if (i_is_LJ .and. j_is_LJ) then
             iHash = ior(iHash, LJ_PAIR)
             


          endif



          if (i_is_Elect .and. j_is_Elect) iHash = ior(iHash, ELECTROSTATIC_PAIR)
          if (i_is_Sticky .and. j_is_Sticky) iHash = ior(iHash, STICKY_PAIR)
          if (i_is_StickyP .and. j_is_StickyP) iHash = ior(iHash, STICKYPOWER_PAIR)

          if (i_is_EAM .and. j_is_EAM) iHash = ior(iHash, EAM_PAIR)

          if (i_is_GB .and. j_is_GB) iHash = ior(iHash, GAYBERNE_PAIR)
          if (i_is_GB .and. j_is_LJ) iHash = ior(iHash, GAYBERNE_LJ)
          if (i_is_LJ .and. j_is_GB) iHash = ior(iHash, GAYBERNE_LJ)

          if (i_is_Shape .and. j_is_Shape) iHash = ior(iHash, SHAPE_PAIR)
          if (i_is_Shape .and. j_is_LJ) iHash = ior(iHash, SHAPE_LJ)
          if (i_is_LJ .and. j_is_Shape) iHash = ior(iHash, SHAPE_LJ)


          InteractionMap(i,j)%InteractionHash = iHash
          InteractionMap(j,i)%InteractionHash = iHash

       end do

    end do
  end subroutine createInteractionMap

end module interactions
