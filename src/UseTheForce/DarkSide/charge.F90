module charge_charge
  
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

  real(kind=dp), parameter :: pre = 332.06508_DP  
  logical, save :: haveChargeMap = .false.

  public::do_charge_pair

  type :: ChargeList
     real(kind=DP) :: charge = 0.0_DP
  end type ChargeList

  type(ChargeList), dimension(:), allocatable :: ChargeMap

contains
      
  subroutine createChargeMap(status)
    integer :: nAtypes
    integer :: status
    integer :: i
    real (kind=DP) :: thisCharge
    logical :: thisProperty
    
    status = 0
    
    nAtypes = getSize(atypes)
    
    if (nAtypes == 0) then
       status = -1
       return
    end if
    
    if (.not. allocated(ChargeMap)) then
       allocate(ChargeMap(nAtypes))
    endif

    do i = 1, nAtypes

       call getElementProperty(atypes, i, "is_Charge", thisProperty)

       if (thisProperty) then
          call getElementProperty(atypes, i, "charge", thisCharge)
          ChargeMap(i)%charge = thisCharge
       endif
       
    end do
    
    haveChargeMap = .true.
    
  end subroutine createChargeMap
    
  subroutine do_charge_pair(atom1, atom2, d, rij, r2, sw, vpair, fpair, &
       pot, f, do_pot)
    
    logical :: do_pot

    integer atom1, atom2, me1, me2, id1, id2
    integer :: localError
    real(kind=dp) :: rij, r2, q1, q2, sw, vpair
    real(kind=dp) :: drdx, drdy, drdz, dudr, fx, fy, fz
    real(kind=dp) :: vterm

    real( kind = dp ) :: pot
    real( kind = dp ), dimension(3) :: d, fpair
    real( kind = dp ), dimension(3,nLocal) :: f
    
    if (.not.haveChargeMap) then
       localError = 0
       call createChargeMap(localError)
       if ( localError .ne. 0 ) then
          call handleError("charge-charge", "ChargeMap creation failed!")
          return
       end if
    endif

#ifdef IS_MPI
    me1 = atid_Row(atom1)
    me2 = atid_Col(atom2)
#else
    me1 = atid(atom1)
    me2 = atid(atom2)
#endif

    q1 = ChargeMap(me1)%charge
    q2 = ChargeMap(me2)%charge

    vterm = pre * q1 * q2 / rij 

    dudr =  -sw * vterm  / rij

    drdx = d(1) / rij
    drdy = d(2) / rij
    drdz = d(3) / rij
    
    fx = dudr * drdx
    fy = dudr * drdy
    fz = dudr * drdz          

    vpair = vpair + vterm
    
#ifdef IS_MPI
    if (do_pot) then
       pot_Row(atom1) = pot_Row(atom1) + sw*vterm*0.5
       pot_Col(atom2) = pot_Col(atom2) + sw*vterm*0.5
    endif
    
    f_Row(1,atom1) = f_Row(1,atom1) + fx
    f_Row(2,atom1) = f_Row(2,atom1) + fy
    f_Row(3,atom1) = f_Row(3,atom1) + fz
    
    f_Col(1,atom2) = f_Col(1,atom2) - fx
    f_Col(2,atom2) = f_Col(2,atom2) - fy
    f_Col(3,atom2) = f_Col(3,atom2) - fz
    
#else

    if (do_pot) pot = pot + sw*vterm
    
    f(1,atom1) = f(1,atom1) + fx
    f(2,atom1) = f(2,atom1) + fy
    f(3,atom1) = f(3,atom1) + fz
    
    f(1,atom2) = f(1,atom2) - fx
    f(2,atom2) = f(2,atom2) - fy
    f(3,atom2) = f(3,atom2) - fz
        
#endif

#ifdef IS_MPI
    id1 = AtomRowToGlobal(atom1)
    id2 = AtomColToGlobal(atom2)
#else
    id1 = atom1
    id2 = atom2
#endif
    
    if (molMembershipList(id1) .ne. molMembershipList(id2)) then
       
       fpair(1) = fpair(1) + fx
       fpair(2) = fpair(2) + fy
       fpair(3) = fpair(3) + fz
       
    endif
    return
  end subroutine do_charge_pair
  
end module charge_charge
