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

module dipole_dipole
  
  use force_globals
  use definitions
  use atype_module
  use vector_class
  use simulation
  use status
#ifdef IS_MPI
  use mpiSimulation
#endif
  implicit none

  PRIVATE

  real(kind=dp), parameter :: pre = 14.38362_dp

  public :: newDipoleType
  public :: do_dipole_pair
  public :: getDipoleMoment

  type :: MomentList
     integer :: c_ident
     real(kind=DP) :: dipole_moment = 0.0_DP
  end type MomentList

  type(MomentList), dimension(:),allocatable :: MomentMap

contains

  subroutine newDipoleType(c_ident, dipole_moment, status)
    integer,intent(in) :: c_ident
    real(kind=dp),intent(in) :: dipole_moment
    integer,intent(out) :: status
    integer :: nAtypes, myATID

    status = 0
    myATID = getFirstMatchingElement(atypes, "c_ident", c_ident)
    
    !! Be simple-minded and assume that we need a MomentMap that
    !! is the same size as the total number of atom types

    if (.not.allocated(MomentMap)) then
       
       nAtypes = getSize(atypes)
    
       if (nAtypes == 0) then
          status = -1
          return
       end if
       
       if (.not. allocated(MomentMap)) then
          allocate(MomentMap(nAtypes))
       endif
       
    end if

    if (myATID .gt. size(MomentMap)) then
       status = -1
       return
    endif
    
    ! set the values for MomentMap for this atom type:

    MomentMap(myATID)%c_ident = c_ident
    MomentMap(myATID)%dipole_moment = dipole_moment
    
  end subroutine newDipoleType

  function getDipoleMoment(atid) result (dm)
    integer, intent(in) :: atid
    integer :: localError
    real(kind=dp) :: dm
    
    if (.not.allocated(MomentMap)) then
       call handleError("dipole-dipole", "no MomentMap was present before first call of getDipoleMoment!")
       return
    end if
    
    dm = MomentMap(atid)%dipole_moment
  end function getDipoleMoment

  subroutine do_dipole_pair(atom1, atom2, d, rij, r2, sw, vpair, fpair, &
       pot, eFrame, f, t, do_pot)
    
    logical :: do_pot

    integer atom1, atom2, me1, me2, id1, id2
    integer :: localError
    real(kind=dp) :: rij, mu1, mu2
    real(kind=dp) :: dfact1, dfact2, dip2, r2, r3, r5
    real(kind=dp) :: dudx, dudy, dudz, dudu1x, dudu1y, dudu1z
    real(kind=dp) :: dudu2x, dudu2y, dudu2z, rdotu1, rdotu2, u1dotu2
    real(kind=dp) :: sw, vpair, vterm

    real( kind = dp ) :: pot
    real( kind = dp ), dimension(3) :: d, fpair
    real( kind = dp ), dimension(9,nLocal) :: eFrame
    real( kind = dp ), dimension(3,nLocal) :: f
    real( kind = dp ), dimension(3,nLocal) :: t
    
    real (kind = dp), dimension(3) :: ul1
    real (kind = dp), dimension(3) :: ul2


    if (.not.allocated(MomentMap)) then
       call handleError("dipole-dipole", "no MomentMap was present before first call of do_dipole_pair!")
       return
    end if


#ifdef IS_MPI
    me1 = atid_Row(atom1)
    ul1(1) = eFrame_Row(3,atom1)
    ul1(2) = eFrame_Row(6,atom1)
    ul1(3) = eFrame_Row(9,atom1)

    me2 = atid_Col(atom2)
    ul2(1) = eFrame_Col(3,atom2)
    ul2(2) = eFrame_Col(6,atom2)
    ul2(3) = eFrame_Col(9,atom2)
#else
    me1 = atid(atom1)
    ul1(1) = eFrame(3,atom1)
    ul1(2) = eFrame(6,atom1)
    ul1(3) = eFrame(9,atom1)

    me2 = atid(atom2)
    ul2(1) = eFrame(3,atom2)
    ul2(2) = eFrame(6,atom2)
    ul2(3) = eFrame(9,atom2)
#endif

    mu1 = MomentMap(me1)%dipole_moment
    mu2 = MomentMap(me2)%dipole_moment
    
    r3 = r2*rij
    r5 = r3*r2
    
    rdotu1 = d(1)*ul1(1) + d(2)*ul1(2) + d(3)*ul1(3)
    rdotu2 = d(1)*ul2(1) + d(2)*ul2(2) + d(3)*ul2(3)
    u1dotu2 = ul1(1)*ul2(1) + ul1(2)*ul2(2) + ul1(3)*ul2(3)
    
    dip2 = pre * mu1 * mu2
    dfact1 = 3.0d0*dip2 / r2
    dfact2 = 3.0d0*dip2 / r5
    
    vterm = dip2*((u1dotu2/r3) - 3.0d0*(rdotu1*rdotu2/r5))
    
    vpair = vpair + vterm
    
    if (do_pot) then
#ifdef IS_MPI 
       pot_row(atom1) = pot_row(atom1) + 0.5d0*vterm*sw
       pot_col(atom2) = pot_col(atom2) + 0.5d0*vterm*sw
#else
       pot = pot + vterm*sw
#endif
    endif
    
    dudx = (-dfact1 * d(1) * ((u1dotu2/r3) - &
         (5.0d0*(rdotu1*rdotu2)/r5)) -  &
            dfact2*(ul1(1)*rdotu2 + ul2(1)*rdotu1))*sw
    
    dudy = (-dfact1 * d(2) * ((u1dotu2/r3) - &
         (5.0d0*(rdotu1*rdotu2)/r5)) -  &
            dfact2*(ul1(2)*rdotu2 + ul2(2)*rdotu1))*sw
    
    dudz = (-dfact1 * d(3) * ((u1dotu2/r3) - &
         (5.0d0*(rdotu1*rdotu2)/r5)) -  &
         dfact2*(ul1(3)*rdotu2 + ul2(3)*rdotu1))*sw
    
    dudu1x = (dip2*((ul2(1)/r3) - (3.0d0*d(1)*rdotu2/r5)))*sw
    dudu1y = (dip2*((ul2(2)/r3) - (3.0d0*d(2)*rdotu2/r5)))*sw
    dudu1z = (dip2*((ul2(3)/r3) - (3.0d0*d(3)*rdotu2/r5)))*sw
    
    dudu2x = (dip2*((ul1(1)/r3) - (3.0d0*d(1)*rdotu1/r5)))*sw
    dudu2y = (dip2*((ul1(2)/r3) - (3.0d0*d(2)*rdotu1/r5)))*sw
    dudu2z = (dip2*((ul1(3)/r3) - (3.0d0*d(3)*rdotu1/r5)))*sw
    
    
#ifdef IS_MPI
    f_Row(1,atom1) = f_Row(1,atom1) + dudx
    f_Row(2,atom1) = f_Row(2,atom1) + dudy
    f_Row(3,atom1) = f_Row(3,atom1) + dudz
    
    f_Col(1,atom2) = f_Col(1,atom2) - dudx
    f_Col(2,atom2) = f_Col(2,atom2) - dudy
    f_Col(3,atom2) = f_Col(3,atom2) - dudz
    
    t_Row(1,atom1) = t_Row(1,atom1) - ul1(2)*dudu1z + ul1(3)*dudu1y
    t_Row(2,atom1) = t_Row(2,atom1) - ul1(3)*dudu1x + ul1(1)*dudu1z
    t_Row(3,atom1) = t_Row(3,atom1) - ul1(1)*dudu1y + ul1(2)*dudu1x
    
    t_Col(1,atom2) = t_Col(1,atom2) - ul2(2)*dudu2z + ul2(3)*dudu2y
    t_Col(2,atom2) = t_Col(2,atom2) - ul2(3)*dudu2x + ul2(1)*dudu2z
    t_Col(3,atom2) = t_Col(3,atom2) - ul2(1)*dudu2y + ul2(2)*dudu2x
#else
    f(1,atom1) = f(1,atom1) + dudx
    f(2,atom1) = f(2,atom1) + dudy
    f(3,atom1) = f(3,atom1) + dudz
    
    f(1,atom2) = f(1,atom2) - dudx
    f(2,atom2) = f(2,atom2) - dudy
    f(3,atom2) = f(3,atom2) - dudz
    
    t(1,atom1) = t(1,atom1) - ul1(2)*dudu1z + ul1(3)*dudu1y
    t(2,atom1) = t(2,atom1) - ul1(3)*dudu1x + ul1(1)*dudu1z
    t(3,atom1) = t(3,atom1) - ul1(1)*dudu1y + ul1(2)*dudu1x
    
    t(1,atom2) = t(1,atom2) - ul2(2)*dudu2z + ul2(3)*dudu2y
    t(2,atom2) = t(2,atom2) - ul2(3)*dudu2x + ul2(1)*dudu2z
    t(3,atom2) = t(3,atom2) - ul2(1)*dudu2y + ul2(2)*dudu2x
#endif
    
#ifdef IS_MPI
    id1 = AtomRowToGlobal(atom1)
    id2 = AtomColToGlobal(atom2)
#else
    id1 = atom1
    id2 = atom2
#endif

    if (molMembershipList(id1) .ne. molMembershipList(id2)) then
       
       fpair(1) = fpair(1) + dudx
       fpair(2) = fpair(2) + dudy
       fpair(3) = fpair(3) + dudz

    endif

    return
  end subroutine do_dipole_pair
  
end module dipole_dipole

