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


module shapes

  use force_globals
  use definitions
  use atype_module
  use vector_class
  use simulation
  use status
  implicit none

  real(kind=dp), external :: getSigma
  real(kind=dp), external :: getEpsilon

  PRIVATE
#define __FORTRAN90
#include "UseTheForce/DarkSide/fInteractionMap.h"
  INTEGER, PARAMETER:: CHEBYSHEV_TN = 1
  INTEGER, PARAMETER:: CHEBYSHEV_UN = 2
  INTEGER, PARAMETER:: LAGUERRE     = 3
  INTEGER, PARAMETER:: HERMITE      = 4
  INTEGER, PARAMETER:: SH_COS       = 0
  INTEGER, PARAMETER:: SH_SIN       = 1

  logical, save :: haveShapeMap = .false.

  public :: do_shape_pair
  public :: newShapeType
  public :: complete_Shape_FF
  public :: destroyShapeTypes
  public :: getShapeCut

  type, private :: Shape
     integer :: atid
     integer :: nContactFuncs 
     integer :: nRangeFuncs 
     integer :: nStrengthFuncs 
     integer :: bigL
     integer :: bigM
     integer, pointer, dimension(:) :: ContactFuncLValue             => null()
     integer, pointer, dimension(:) :: ContactFuncMValue             => null()
     integer, pointer, dimension(:) :: ContactFunctionType           => null()
     real(kind=dp), pointer, dimension(:) :: ContactFuncCoefficient  => null()
     integer, pointer, dimension(:) :: RangeFuncLValue               => null()
     integer, pointer, dimension(:) :: RangeFuncMValue               => null()
     integer, pointer, dimension(:) :: RangeFunctionType             => null()
     real(kind=dp), pointer, dimension(:) :: RangeFuncCoefficient    => null()
     integer, pointer, dimension(:) :: StrengthFuncLValue            => null()
     integer, pointer, dimension(:) :: StrengthFuncMValue            => null()
     integer, pointer, dimension(:) :: StrengthFunctionType          => null()
     real(kind=dp), pointer, dimension(:) :: StrengthFuncCoefficient => null()
     logical :: isLJ 
     real ( kind = dp )  :: epsilon 
     real ( kind = dp )  :: sigma 
  end type Shape

  type, private :: ShapeList
     integer :: n_shapes = 0
     integer :: currentShape = 0
     type(Shape), pointer :: Shapes(:)      => null()
     integer, pointer     :: atidToShape(:) => null()
  end type ShapeList
  
  type(ShapeList), save :: ShapeMap
  
  integer :: lmax
  
contains  
  
  subroutine newShapeType(nContactFuncs, ContactFuncLValue, &
       ContactFuncMValue, ContactFunctionType, ContactFuncCoefficient, &
       nRangeFuncs, RangeFuncLValue, RangeFuncMValue, RangeFunctionType, &
       RangeFuncCoefficient, nStrengthFuncs, StrengthFuncLValue, &
       StrengthFuncMValue, StrengthFunctionType, StrengthFuncCoefficient, &
       c_ident, status)
    
    integer :: nContactFuncs 
    integer :: nRangeFuncs 
    integer :: nStrengthFuncs 
    integer :: shape_ident
    integer :: status
    integer :: c_ident
    integer :: myATID
    integer :: bigL
    integer :: bigM
    integer :: j, me, nShapeTypes, nLJTypes, ntypes, current, alloc_stat
    integer, pointer :: MatchList(:) => null()

    integer, dimension(nContactFuncs) :: ContactFuncLValue           
    integer, dimension(nContactFuncs) :: ContactFuncMValue           
    integer, dimension(nContactFuncs) :: ContactFunctionType         
    real(kind=dp), dimension(nContactFuncs) :: ContactFuncCoefficient
    integer, dimension(nRangeFuncs) :: RangeFuncLValue             
    integer, dimension(nRangeFuncs) :: RangeFuncMValue             
    integer, dimension(nRangeFuncs) :: RangeFunctionType           
    real(kind=dp), dimension(nRangeFuncs) :: RangeFuncCoefficient  
    integer, dimension(nStrengthFuncs) :: StrengthFuncLValue          
    integer, dimension(nStrengthFuncs) :: StrengthFuncMValue          
    integer, dimension(nStrengthFuncs) :: StrengthFunctionType        
    real(kind=dp), dimension(nStrengthFuncs) :: StrengthFuncCoefficient

    status = 0
    ! check to see if this is the first time into this routine...
    if (.not.associated(ShapeMap%Shapes)) then

       call getMatchingElementList(atypes, "is_Shape", .true., &
            nShapeTypes, MatchList)

       call getMatchingElementList(atypes, "is_LennardJones", .true., &
            nLJTypes, MatchList)

       ShapeMap%n_shapes = nShapeTypes + nLJTypes

       allocate(ShapeMap%Shapes(nShapeTypes + nLJTypes))

       ntypes = getSize(atypes)
 
       allocate(ShapeMap%atidToShape(ntypes))
    end if

    ShapeMap%currentShape = ShapeMap%currentShape + 1
    current = ShapeMap%currentShape

    call allocateShape(nContactFuncs, nRangeFuncs, nStrengthFuncs, &
         ShapeMap%Shapes(current), stat=alloc_stat)
    if (alloc_stat .ne. 0) then
       status = -1
       return
    endif

    myATID = getFirstMatchingElement(atypes, "c_ident", c_ident)

    ShapeMap%atidToShape(myATID)                     = current
    ShapeMap%Shapes(current)%atid                    = myATID
    ShapeMap%Shapes(current)%nContactFuncs           = nContactFuncs
    ShapeMap%Shapes(current)%nRangeFuncs             = nRangeFuncs
    ShapeMap%Shapes(current)%nStrengthFuncs          = nStrengthFuncs
    ShapeMap%Shapes(current)%ContactFuncLValue       = ContactFuncLValue
    ShapeMap%Shapes(current)%ContactFuncMValue       = ContactFuncMValue
    ShapeMap%Shapes(current)%ContactFunctionType     = ContactFunctionType
    ShapeMap%Shapes(current)%ContactFuncCoefficient  = ContactFuncCoefficient
    ShapeMap%Shapes(current)%RangeFuncLValue         = RangeFuncLValue
    ShapeMap%Shapes(current)%RangeFuncMValue         = RangeFuncMValue
    ShapeMap%Shapes(current)%RangeFunctionType       = RangeFunctionType
    ShapeMap%Shapes(current)%RangeFuncCoefficient    = RangeFuncCoefficient
    ShapeMap%Shapes(current)%StrengthFuncLValue      = StrengthFuncLValue
    ShapeMap%Shapes(current)%StrengthFuncMValue      = StrengthFuncMValue
    ShapeMap%Shapes(current)%StrengthFunctionType    = StrengthFunctionType
    ShapeMap%Shapes(current)%StrengthFuncCoefficient = StrengthFuncCoefficient

    bigL = -1
    bigM = -1

    do j = 1, ShapeMap%Shapes(current)%nContactFuncs
       if (ShapeMap%Shapes(current)%ContactFuncLValue(j) .gt. bigL) then
          bigL = ShapeMap%Shapes(current)%ContactFuncLValue(j)
       endif
       if (ShapeMap%Shapes(current)%ContactFuncMValue(j) .gt. bigM) then
          bigM = ShapeMap%Shapes(current)%ContactFuncMValue(j)
       endif
    enddo
    do j = 1, ShapeMap%Shapes(current)%nRangeFuncs
       if (ShapeMap%Shapes(current)%RangeFuncLValue(j) .gt. bigL) then
          bigL = ShapeMap%Shapes(current)%RangeFuncLValue(j)
       endif
       if (ShapeMap%Shapes(current)%RangeFuncMValue(j) .gt. bigM) then
          bigM = ShapeMap%Shapes(current)%RangeFuncMValue(j)
       endif
    enddo
    do j = 1, ShapeMap%Shapes(current)%nStrengthFuncs
       if (ShapeMap%Shapes(current)%StrengthFuncLValue(j) .gt. bigL) then
          bigL = ShapeMap%Shapes(current)%StrengthFuncLValue(j)
       endif
       if (ShapeMap%Shapes(current)%StrengthFuncMValue(j) .gt. bigM) then
          bigM = ShapeMap%Shapes(current)%StrengthFuncMValue(j)
       endif
    enddo

    ShapeMap%Shapes(current)%bigL                    = bigL
    ShapeMap%Shapes(current)%bigM                    = bigM

  end subroutine newShapeType

  subroutine allocateShape(nContactFuncs, nRangeFuncs, nStrengthFuncs, &
       myShape, stat)

    integer, intent(in) :: nContactFuncs, nRangeFuncs, nStrengthFuncs
    type(Shape), intent(inout) :: myShape
    integer, intent(out) :: stat
    integer :: alloc_stat

    stat = 0
    if (associated(myShape%contactFuncLValue)) then
       deallocate(myShape%contactFuncLValue)
    endif
    allocate(myShape%contactFuncLValue(nContactFuncs), stat = alloc_stat)
    if (alloc_stat .ne. 0) then
       stat = -1
       return
    endif
    if (associated(myShape%contactFuncMValue)) then
       deallocate(myShape%contactFuncMValue)
    endif
    allocate(myShape%contactFuncMValue(nContactFuncs), stat = alloc_stat)
    if (alloc_stat .ne. 0) then
       stat = -1
       return
    endif
    if (associated(myShape%contactFunctionType)) then
       deallocate(myShape%contactFunctionType)
    endif
    allocate(myShape%contactFunctionType(nContactFuncs), stat = alloc_stat)
    if (alloc_stat .ne. 0) then
       stat = -1
       return
    endif
    if (associated(myShape%contactFuncCoefficient)) then
       deallocate(myShape%contactFuncCoefficient)
    endif
    allocate(myShape%contactFuncCoefficient(nContactFuncs), stat = alloc_stat)
    if (alloc_stat .ne. 0) then
       stat = -1
       return
    endif

    if (associated(myShape%rangeFuncLValue)) then
       deallocate(myShape%rangeFuncLValue)
    endif
    allocate(myShape%rangeFuncLValue(nRangeFuncs), stat = alloc_stat)
    if (alloc_stat .ne. 0) then
       stat = -1
       return
    endif
    if (associated(myShape%rangeFuncMValue)) then
       deallocate(myShape%rangeFuncMValue)
    endif
    allocate(myShape%rangeFuncMValue(nRangeFuncs), stat = alloc_stat)
    if (alloc_stat .ne. 0) then
       stat = -1
       return
    endif
    if (associated(myShape%rangeFunctionType)) then
       deallocate(myShape%rangeFunctionType)
    endif
    allocate(myShape%rangeFunctionType(nRangeFuncs), stat = alloc_stat)
    if (alloc_stat .ne. 0) then
       stat = -1
       return
    endif
    if (associated(myShape%rangeFuncCoefficient)) then
       deallocate(myShape%rangeFuncCoefficient)
    endif
    allocate(myShape%rangeFuncCoefficient(nRangeFuncs), stat = alloc_stat)
    if (alloc_stat .ne. 0) then
       stat = -1
       return
    endif

    if (associated(myShape%strengthFuncLValue)) then
       deallocate(myShape%strengthFuncLValue)
    endif
    allocate(myShape%strengthFuncLValue(nStrengthFuncs), stat = alloc_stat)
    if (alloc_stat .ne. 0) then
       stat = -1
       return
    endif
    if (associated(myShape%strengthFuncMValue)) then
       deallocate(myShape%strengthFuncMValue)
    endif
    allocate(myShape%strengthFuncMValue(nStrengthFuncs), stat = alloc_stat)
    if (alloc_stat .ne. 0) then
       stat = -1
       return
    endif
    if (associated(myShape%strengthFunctionType)) then
       deallocate(myShape%strengthFunctionType)
    endif
    allocate(myShape%strengthFunctionType(nStrengthFuncs), stat = alloc_stat)
    if (alloc_stat .ne. 0) then
       stat = -1
       return
    endif
    if (associated(myShape%strengthFuncCoefficient)) then
       deallocate(myShape%strengthFuncCoefficient)
    endif
    allocate(myShape%strengthFuncCoefficient(nStrengthFuncs), stat=alloc_stat)
    if (alloc_stat .ne. 0) then
       stat = -1
       return
    endif

    return

  end subroutine allocateShape

  subroutine complete_Shape_FF(status)
    integer :: status
    integer :: i, j, l, m, lm, function_type
    real(kind=dp) :: thisDP, sigma
    integer :: alloc_stat, iTheta, iPhi, nSteps, nAtypes, myATID, current
    logical :: thisProperty

    status = 0
    if (ShapeMap%currentShape == 0) then
       call handleError("init_Shape_FF", "No members in ShapeMap")
       status = -1
       return
    end if

    nAtypes = getSize(atypes)

    if (nAtypes == 0) then
       status = -1
       return
    end if       

    ! atypes comes from c side
    do i = 1, nAtypes
    
       myATID = getFirstMatchingElement(atypes, 'c_ident', i)
       call getElementProperty(atypes, myATID, "is_LennardJones", thisProperty)
         
       if (thisProperty) then
          ShapeMap%currentShape = ShapeMap%currentShape + 1
          current = ShapeMap%currentShape

          ShapeMap%atidToShape(myATID) = current
          ShapeMap%Shapes(current)%atid = myATID

          ShapeMap%Shapes(current)%isLJ = .true.

          ShapeMap%Shapes(current)%epsilon = getEpsilon(myATID)
          ShapeMap%Shapes(current)%sigma = getSigma(myATID)

       endif

    end do

    haveShapeMap = .true.

!    do i = 1, ShapeMap%n_shapes
!       write(*,*) 'i = ', i, ' isLJ = ', ShapeMap%Shapes(i)%isLJ
!    end do

  end subroutine complete_Shape_FF

  function getShapeCut(atomID) result(cutValue)
    integer, intent(in) :: atomID
    real(kind=dp) :: cutValue, whoopdedoo

    !! this is just a placeholder for a cutoff value, hopefully we'll 
    !! develop a method to calculate a sensible value
    whoopdedoo = 9.0_dp

    cutValue = whoopdedoo

  end function getShapeCut

  subroutine do_shape_pair(atid1, atid2, d, rij, r2, sw, &
       vpair, pot, A1, A2, f1, t1, t2)

    INTEGER, PARAMETER:: LMAX         = 64
    INTEGER, PARAMETER:: MMAX         = 64

    integer, intent(in) :: atid1, atid2
    real (kind=dp), intent(inout) :: rij, r2
    real (kind=dp), dimension(3), intent(in) :: d
    real (kind=dp) :: pot, vpair, sw, dswdr
    real (kind=dp), dimension(9) :: A1, A2
    real (kind=dp), dimension(3) :: f1
    real (kind=dp), dimension(3) :: t1, t2

    real (kind=dp) :: r3, r5, rt2, rt3, rt5, rt6, rt11, rt12, rt126
    integer :: st1, st2
    integer :: l, m, lm, id1, id2, localError, function_type
    real (kind=dp) :: sigma_i, s_i, eps_i, sigma_j, s_j, eps_j
    real (kind=dp) :: coeff
    real (kind=dp) :: pot_temp

    real (kind=dp) :: dsigmaidx, dsigmaidy, dsigmaidz
    real (kind=dp) :: dsigmaidux, dsigmaiduy, dsigmaiduz
    real (kind=dp) :: dsigmajdx, dsigmajdy, dsigmajdz
    real (kind=dp) :: dsigmajdux, dsigmajduy, dsigmajduz

    real (kind=dp) :: dsidx, dsidy, dsidz
    real (kind=dp) :: dsidux, dsiduy, dsiduz
    real (kind=dp) :: dsjdx, dsjdy, dsjdz
    real (kind=dp) :: dsjdux, dsjduy, dsjduz

    real (kind=dp) :: depsidx, depsidy, depsidz
    real (kind=dp) :: depsidux, depsiduy, depsiduz
    real (kind=dp) :: depsjdx, depsjdy, depsjdz
    real (kind=dp) :: depsjdux, depsjduy, depsjduz

    real (kind=dp) :: xi, yi, zi, xj, yj, zj, xi2, yi2, zi2, xj2, yj2, zj2

    real (kind=dp) :: sti2, stj2

    real (kind=dp) :: proji, proji3, projj, projj3
    real (kind=dp) :: cti, ctj, cpi, cpj, spi, spj
    real (kind=dp) :: Phunc, sigma, s, eps, rtdenom, rt

    real (kind=dp) :: dctidx, dctidy, dctidz
    real (kind=dp) :: dctidux, dctiduy, dctiduz
    real (kind=dp) :: dctjdx, dctjdy, dctjdz
    real (kind=dp) :: dctjdux, dctjduy, dctjduz

    real (kind=dp) :: dcpidx, dcpidy, dcpidz
    real (kind=dp) :: dcpidux, dcpiduy, dcpiduz
    real (kind=dp) :: dcpjdx, dcpjdy, dcpjdz
    real (kind=dp) :: dcpjdux, dcpjduy, dcpjduz

    real (kind=dp) :: dspidx, dspidy, dspidz
    real (kind=dp) :: dspidux, dspiduy, dspiduz
    real (kind=dp) :: dspjdx, dspjdy, dspjdz
    real (kind=dp) :: dspjdux, dspjduy, dspjduz

    real (kind=dp) :: dPhuncdX, dPhuncdY, dPhuncdZ
    real (kind=dp) :: dPhuncdUx, dPhuncdUy, dPhuncdUz

    real (kind=dp) :: dsigmadxi, dsigmadyi, dsigmadzi
    real (kind=dp) :: dsigmaduxi, dsigmaduyi, dsigmaduzi
    real (kind=dp) :: dsigmadxj, dsigmadyj, dsigmadzj
    real (kind=dp) :: dsigmaduxj, dsigmaduyj, dsigmaduzj

    real (kind=dp) :: dsdxi, dsdyi, dsdzi
    real (kind=dp) :: dsduxi, dsduyi, dsduzi
    real (kind=dp) :: dsdxj, dsdyj, dsdzj
    real (kind=dp) :: dsduxj, dsduyj, dsduzj

    real (kind=dp) :: depsdxi, depsdyi, depsdzi
    real (kind=dp) :: depsduxi, depsduyi, depsduzi
    real (kind=dp) :: depsdxj, depsdyj, depsdzj
    real (kind=dp) :: depsduxj, depsduyj, depsduzj

    real (kind=dp) :: drtdxi, drtdyi, drtdzi
    real (kind=dp) :: drtduxi, drtduyi, drtduzi
    real (kind=dp) :: drtdxj, drtdyj, drtdzj
    real (kind=dp) :: drtduxj, drtduyj, drtduzj

    real (kind=dp) :: drdxi, drdyi, drdzi
    real (kind=dp) :: drduxi, drduyi, drduzi
    real (kind=dp) :: drdxj, drdyj, drdzj
    real (kind=dp) :: drduxj, drduyj, drduzj

    real (kind=dp) :: dvdxi, dvdyi, dvdzi
    real (kind=dp) :: dvduxi, dvduyi, dvduzi
    real (kind=dp) :: dvdxj, dvdyj, dvdzj
    real (kind=dp) :: dvduxj, dvduyj, dvduzj  

    real (kind=dp) :: fxi, fyi, fzi, fxj, fyj, fzj
    real (kind=dp) :: txi, tyi, tzi, txj, tyj, tzj
    real (kind=dp) :: fxii, fyii, fzii, fxij, fyij, fzij
    real (kind=dp) :: fxji, fyji, fzji, fxjj, fyjj, fzjj
    real (kind=dp) :: fxradial, fyradial, fzradial

    real (kind=dp) :: xihat, yihat, zihat, xjhat, yjhat, zjhat

    real (kind=dp) :: plm_i(0:LMAX,0:MMAX), dlm_i(0:LMAX,0:MMAX)
    real (kind=dp) :: plm_j(0:LMAX,0:MMAX), dlm_j(0:LMAX,0:MMAX)
    real (kind=dp) :: tm_i(0:MMAX), dtm_i(0:MMAX), um_i(0:MMAX), dum_i(0:MMAX)
    real (kind=dp) :: tm_j(0:MMAX), dtm_j(0:MMAX), um_j(0:MMAX), dum_j(0:MMAX)

    if (.not.haveShapeMap) then
       call handleError("calc_shape", "NO SHAPEMAP!!!!")
       return       
    endif

    !! We assume that the rotation matrices have already been calculated
    !! and placed in the A array.
    r3 = r2*rij
    r5 = r3*r2

    drdxi = -d(1) / rij
    drdyi = -d(2) / rij
    drdzi = -d(3) / rij
    drduxi = 0.0_dp
    drduyi = 0.0_dp
    drduzi = 0.0_dp

    drdxj = d(1) / rij
    drdyj = d(2) / rij
    drdzj = d(3) / rij
    drduxj = 0.0_dp
    drduyj = 0.0_dp
    drduzj = 0.0_dp

    ! use the atid to find the shape type (st) for each atom:
    st1 = ShapeMap%atidToShape(atid1)
    st2 = ShapeMap%atidToShape(atid2)   

    if (ShapeMap%Shapes(st1)%isLJ) then

       sigma_i = ShapeMap%Shapes(st1)%sigma
       s_i = ShapeMap%Shapes(st1)%sigma
       eps_i = ShapeMap%Shapes(st1)%epsilon
       dsigmaidx = 0.0_dp
       dsigmaidy = 0.0_dp
       dsigmaidz = 0.0_dp
       dsigmaidux = 0.0_dp
       dsigmaiduy = 0.0_dp
       dsigmaiduz = 0.0_dp
       dsidx = 0.0_dp
       dsidy = 0.0_dp
       dsidz = 0.0_dp
       dsidux = 0.0_dp
       dsiduy = 0.0_dp
       dsiduz = 0.0_dp
       depsidx = 0.0_dp
       depsidy = 0.0_dp
       depsidz = 0.0_dp
       depsidux = 0.0_dp
       depsiduy = 0.0_dp
       depsiduz = 0.0_dp
    else
       
       ! rotate the inter-particle separation into the two different
       ! body-fixed coordinate systems:

       xi = A1(1)*d(1) + A1(2)*d(2) + A1(3)*d(3)
       yi = A1(4)*d(1) + A1(5)*d(2) + A1(6)*d(3)
       zi = A1(7)*d(1) + A1(8)*d(2) + A1(9)*d(3)
       
       xihat = xi / rij
       yihat = yi / rij
       zihat = zi / rij
       xi2 = xi*xi
       yi2 = yi*yi
       zi2 = zi*zi             
       cti = zi / rij

       if (cti .gt. 1.0_dp) cti = 1.0_dp
       if (cti .lt. -1.0_dp) cti = -1.0_dp

       dctidx = - zi * xi / r3
       dctidy = - zi * yi / r3
       dctidz = 1.0_dp / rij - zi2 / r3
       dctidux = yi / rij ! - (zi * xi2) / r3
       dctiduy = -xi / rij !- (zi * yi2) / r3
       dctiduz = 0.0_dp !zi / rij - (zi2 * zi) / r3

       ! this is an attempt to try to truncate the singularity when
       ! sin(theta) is near 0.0:

       sti2 = 1.0_dp - cti*cti
       if (abs(sti2) .lt. 1.0e-12_dp) then
          proji = sqrt(rij * 1.0e-12_dp)
          dcpidx = 1.0_dp / proji
          dcpidy = 0.0_dp
          dcpidux = xi / proji
          dcpiduy = 0.0_dp
          dspidx = 0.0_dp
          dspidy = 1.0_dp / proji
          dspidux = 0.0_dp
          dspiduy = yi / proji
       else
          proji = sqrt(xi2 + yi2)
          proji3 = proji*proji*proji
          dcpidx = 1.0_dp / proji - xi2 / proji3
          dcpidy = - xi * yi / proji3
          dcpidux = xi / proji - (xi2 * xi) / proji3
          dcpiduy = - (xi * yi2) / proji3
          dspidx = - xi * yi / proji3
          dspidy = 1.0_dp / proji - yi2 / proji3
          dspidux = - (yi * xi2) / proji3
          dspiduy = yi / proji - (yi2 * yi) / proji3
       endif

       cpi = xi / proji
       dcpidz = 0.0_dp
       dcpiduz = 0.0_dp

       spi = yi / proji
       dspidz = 0.0_dp
       dspiduz = 0.0_dp

       call Associated_Legendre(cti, ShapeMap%Shapes(st1)%bigM, &
            ShapeMap%Shapes(st1)%bigL, LMAX, &
            plm_i, dlm_i)

       call Orthogonal_Polynomial(cpi, ShapeMap%Shapes(st1)%bigM, MMAX, &
            CHEBYSHEV_TN, tm_i, dtm_i)
       call Orthogonal_Polynomial(cpi, ShapeMap%Shapes(st1)%bigM, MMAX, &
            CHEBYSHEV_UN, um_i, dum_i)

       sigma_i = 0.0_dp
       s_i = 0.0_dp
       eps_i = 0.0_dp
       dsigmaidx = 0.0_dp
       dsigmaidy = 0.0_dp
       dsigmaidz = 0.0_dp
       dsigmaidux = 0.0_dp
       dsigmaiduy = 0.0_dp
       dsigmaiduz = 0.0_dp
       dsidx = 0.0_dp
       dsidy = 0.0_dp
       dsidz = 0.0_dp
       dsidux = 0.0_dp
       dsiduy = 0.0_dp
       dsiduz = 0.0_dp
       depsidx = 0.0_dp
       depsidy = 0.0_dp
       depsidz = 0.0_dp
       depsidux = 0.0_dp
       depsiduy = 0.0_dp
       depsiduz = 0.0_dp

       do lm = 1, ShapeMap%Shapes(st1)%nContactFuncs
          l = ShapeMap%Shapes(st1)%ContactFuncLValue(lm)
          m = ShapeMap%Shapes(st1)%ContactFuncMValue(lm)
          coeff = ShapeMap%Shapes(st1)%ContactFuncCoefficient(lm)
          function_type = ShapeMap%Shapes(st1)%ContactFunctionType(lm)

          if ((function_type .eq. SH_COS).or.(m.eq.0)) then
             Phunc = coeff * tm_i(m)
             dPhuncdX = coeff * dtm_i(m) * dcpidx
             dPhuncdY = coeff * dtm_i(m) * dcpidy
             dPhuncdZ = coeff * dtm_i(m) * dcpidz
             dPhuncdUx = coeff * dtm_i(m) * dcpidux
             dPhuncdUy = coeff * dtm_i(m) * dcpiduy
             dPhuncdUz = coeff * dtm_i(m) * dcpiduz
          else
             Phunc = coeff * spi * um_i(m-1)
             dPhuncdX = coeff * (spi * dum_i(m-1) * dcpidx + dspidx *um_i(m-1))
             dPhuncdY = coeff * (spi * dum_i(m-1) * dcpidy + dspidy *um_i(m-1))
             dPhuncdZ = coeff * (spi * dum_i(m-1) * dcpidz + dspidz *um_i(m-1))
             dPhuncdUx = coeff*(spi * dum_i(m-1)*dcpidux + dspidux *um_i(m-1))
             dPhuncdUy = coeff*(spi * dum_i(m-1)*dcpiduy + dspiduy *um_i(m-1))
             dPhuncdUz = coeff*(spi * dum_i(m-1)*dcpiduz + dspiduz *um_i(m-1))
          endif

          sigma_i = sigma_i + plm_i(m,l)*Phunc
!!$          write(*,*) 'dsigmaidux = ', dsigmaidux
!!$          write(*,*) 'Phunc = ', Phunc
          dsigmaidx = dsigmaidx + plm_i(m,l)*dPhuncdX + &
               Phunc * dlm_i(m,l) * dctidx
          dsigmaidy = dsigmaidy + plm_i(m,l)*dPhuncdY + &
               Phunc * dlm_i(m,l) * dctidy
          dsigmaidz = dsigmaidz + plm_i(m,l)*dPhuncdZ + &
               Phunc * dlm_i(m,l) * dctidz
          dsigmaidux = dsigmaidux + plm_i(m,l)* dPhuncdUx + &
               Phunc * dlm_i(m,l) * dctidux
          dsigmaiduy = dsigmaiduy + plm_i(m,l)* dPhuncdUy + &
               Phunc * dlm_i(m,l) * dctiduy
          dsigmaiduz = dsigmaiduz + plm_i(m,l)* dPhuncdUz + &
               Phunc * dlm_i(m,l) * dctiduz
!!$          write(*,*) 'dsigmaidux = ', dsigmaidux, '; dPhuncdUx = ', dPhuncdUx, &
!!$                     '; dctidux = ', dctidux, '; plm_i(m,l) = ', plm_i(m,l), &
!!$                     '; dlm_i(m,l) = ', dlm_i(m,l), '; m = ', m, '; l = ', l 
       end do

       do lm = 1, ShapeMap%Shapes(st1)%nRangeFuncs
          l = ShapeMap%Shapes(st1)%RangeFuncLValue(lm)
          m = ShapeMap%Shapes(st1)%RangeFuncMValue(lm)
          coeff = ShapeMap%Shapes(st1)%RangeFuncCoefficient(lm)
          function_type = ShapeMap%Shapes(st1)%RangeFunctionType(lm)

          if ((function_type .eq. SH_COS).or.(m.eq.0)) then
             Phunc = coeff * tm_i(m)
             dPhuncdX = coeff * dtm_i(m) * dcpidx
             dPhuncdY = coeff * dtm_i(m) * dcpidy
             dPhuncdZ = coeff * dtm_i(m) * dcpidz
             dPhuncdUx = coeff * dtm_i(m) * dcpidux
             dPhuncdUy = coeff * dtm_i(m) * dcpiduy
             dPhuncdUz = coeff * dtm_i(m) * dcpiduz
          else
             Phunc = coeff * spi * um_i(m-1)
             dPhuncdX = coeff * (spi * dum_i(m-1) * dcpidx + dspidx *um_i(m-1))
             dPhuncdY = coeff * (spi * dum_i(m-1) * dcpidy + dspidy *um_i(m-1))
             dPhuncdZ = coeff * (spi * dum_i(m-1) * dcpidz + dspidz *um_i(m-1))
             dPhuncdUx = coeff*(spi * dum_i(m-1)*dcpidux + dspidux *um_i(m-1))
             dPhuncdUy = coeff*(spi * dum_i(m-1)*dcpiduy + dspiduy *um_i(m-1))
             dPhuncdUz = coeff*(spi * dum_i(m-1)*dcpiduz + dspiduz *um_i(m-1))
          endif

          s_i = s_i + plm_i(m,l)*Phunc

          dsidx = dsidx + plm_i(m,l)*dPhuncdX + &
               Phunc * dlm_i(m,l) * dctidx
          dsidy = dsidy + plm_i(m,l)*dPhuncdY + &
               Phunc * dlm_i(m,l) * dctidy
          dsidz = dsidz + plm_i(m,l)*dPhuncdZ + &
               Phunc * dlm_i(m,l) * dctidz

          dsidux = dsidux + plm_i(m,l)* dPhuncdUx + &
               Phunc * dlm_i(m,l) * dctidux
          dsiduy = dsiduy + plm_i(m,l)* dPhuncdUy + &
               Phunc * dlm_i(m,l) * dctiduy
          dsiduz = dsiduz + plm_i(m,l)* dPhuncdUz + &
               Phunc * dlm_i(m,l) * dctiduz      

       end do

       do lm = 1, ShapeMap%Shapes(st1)%nStrengthFuncs
          l = ShapeMap%Shapes(st1)%StrengthFuncLValue(lm)
          m = ShapeMap%Shapes(st1)%StrengthFuncMValue(lm)
          coeff = ShapeMap%Shapes(st1)%StrengthFuncCoefficient(lm)
          function_type = ShapeMap%Shapes(st1)%StrengthFunctionType(lm)

          if ((function_type .eq. SH_COS).or.(m.eq.0)) then
             Phunc = coeff * tm_i(m)
             dPhuncdX = coeff * dtm_i(m) * dcpidx
             dPhuncdY = coeff * dtm_i(m) * dcpidy
             dPhuncdZ = coeff * dtm_i(m) * dcpidz
             dPhuncdUx = coeff * dtm_i(m) * dcpidux
             dPhuncdUy = coeff * dtm_i(m) * dcpiduy
             dPhuncdUz = coeff * dtm_i(m) * dcpiduz
          else
             Phunc = coeff * spi * um_i(m-1)
             dPhuncdX = coeff * (spi * dum_i(m-1) * dcpidx + dspidx *um_i(m-1))
             dPhuncdY = coeff * (spi * dum_i(m-1) * dcpidy + dspidy *um_i(m-1))
             dPhuncdZ = coeff * (spi * dum_i(m-1) * dcpidz + dspidz *um_i(m-1))
             dPhuncdUx = coeff*(spi * dum_i(m-1)*dcpidux + dspidux *um_i(m-1))
             dPhuncdUy = coeff*(spi * dum_i(m-1)*dcpiduy + dspiduy *um_i(m-1))
             dPhuncdUz = coeff*(spi * dum_i(m-1)*dcpiduz + dspiduz *um_i(m-1))
          endif

          eps_i = eps_i + plm_i(m,l)*Phunc

          depsidx = depsidx + plm_i(m,l)*dPhuncdX + &
               Phunc * dlm_i(m,l) * dctidx
          depsidy = depsidy + plm_i(m,l)*dPhuncdY + &
               Phunc * dlm_i(m,l) * dctidy
          depsidz = depsidz + plm_i(m,l)*dPhuncdZ + &
               Phunc * dlm_i(m,l) * dctidz

          depsidux = depsidux + plm_i(m,l)* dPhuncdUx + &
               Phunc * dlm_i(m,l) * dctidux
          depsiduy = depsiduy + plm_i(m,l)* dPhuncdUy + &
               Phunc * dlm_i(m,l) * dctiduy
          depsiduz = depsiduz + plm_i(m,l)* dPhuncdUz + &
               Phunc * dlm_i(m,l) * dctiduz      

       end do

    endif

    ! now do j:

    if (ShapeMap%Shapes(st2)%isLJ) then
       sigma_j = ShapeMap%Shapes(st2)%sigma
       s_j = ShapeMap%Shapes(st2)%sigma
       eps_j = ShapeMap%Shapes(st2)%epsilon
       dsigmajdx = 0.0_dp
       dsigmajdy = 0.0_dp
       dsigmajdz = 0.0_dp
       dsigmajdux = 0.0_dp
       dsigmajduy = 0.0_dp
       dsigmajduz = 0.0_dp
       dsjdx = 0.0_dp
       dsjdy = 0.0_dp
       dsjdz = 0.0_dp
       dsjdux = 0.0_dp
       dsjduy = 0.0_dp
       dsjduz = 0.0_dp
       depsjdx = 0.0_dp
       depsjdy = 0.0_dp
       depsjdz = 0.0_dp
       depsjdux = 0.0_dp
       depsjduy = 0.0_dp
       depsjduz = 0.0_dp
    else

       ! rotate the inter-particle separation into the two different
       ! body-fixed coordinate systems:
       ! negative sign because this is the vector from j to i:
       
       xj = -(A2(1)*d(1) + A2(2)*d(2) + A2(3)*d(3))
       yj = -(A2(4)*d(1) + A2(5)*d(2) + A2(6)*d(3))
       zj = -(A2(7)*d(1) + A2(8)*d(2) + A2(9)*d(3))

       xjhat = xj / rij
       yjhat = yj / rij
       zjhat = zj / rij
       xj2 = xj*xj
       yj2 = yj*yj
       zj2 = zj*zj
       ctj = zj / rij

       if (ctj .gt. 1.0_dp) ctj = 1.0_dp
       if (ctj .lt. -1.0_dp) ctj = -1.0_dp

       dctjdx = - zj * xj / r3
       dctjdy = - zj * yj / r3
       dctjdz = 1.0_dp / rij - zj2 / r3
       dctjdux = yj / rij !- (zi * xj2) / r3
       dctjduy = -xj / rij !- (zj * yj2) / r3
       dctjduz = 0.0_dp !zj / rij - (zj2 * zj) / r3

       ! this is an attempt to try to truncate the singularity when
       ! sin(theta) is near 0.0:

       stj2 = 1.0_dp - ctj*ctj
       if (abs(stj2) .lt. 1.0e-12_dp) then
          projj = sqrt(rij * 1.0e-12_dp)
          dcpjdx = 1.0_dp / projj 
          dcpjdy = 0.0_dp
          dcpjdux = xj / projj
          dcpjduy = 0.0_dp
          dspjdx = 0.0_dp
          dspjdy = 1.0_dp / projj
          dspjdux = 0.0_dp
          dspjduy = yj / projj
       else
          projj = sqrt(xj2 + yj2)
          projj3 = projj*projj*projj
          dcpjdx = 1.0_dp / projj - xj2 / projj3
          dcpjdy = - xj * yj / projj3
          dcpjdux = xj / projj - (xj2 * xj) / projj3
          dcpjduy = - (xj * yj2) / projj3
          dspjdx = - xj * yj / projj3
          dspjdy = 1.0_dp / projj - yj2 / projj3
          dspjdux = - (yj * xj2) / projj3
          dspjduy = yj / projj - (yj2 * yj) / projj3
       endif

       cpj = xj / projj
       dcpjdz = 0.0_dp
       dcpjduz = 0.0_dp

       spj = yj / projj
       dspjdz = 0.0_dp
       dspjduz = 0.0_dp


!       write(*,*) 'dcpdu = ' ,dcpidux, dcpiduy, dcpiduz
!       write(*,*) 'dcpdu = ' ,dcpjdux, dcpjduy, dcpjduz
       call Associated_Legendre(ctj, ShapeMap%Shapes(st2)%bigM, &
            ShapeMap%Shapes(st2)%bigL, LMAX, &
            plm_j, dlm_j)

       call Orthogonal_Polynomial(cpj, ShapeMap%Shapes(st2)%bigM, MMAX, &
            CHEBYSHEV_TN, tm_j, dtm_j)
       call Orthogonal_Polynomial(cpj, ShapeMap%Shapes(st2)%bigM, MMAX, &
            CHEBYSHEV_UN, um_j, dum_j)

       sigma_j = 0.0_dp
       s_j = 0.0_dp
       eps_j = 0.0_dp
       dsigmajdx = 0.0_dp
       dsigmajdy = 0.0_dp
       dsigmajdz = 0.0_dp
       dsigmajdux = 0.0_dp
       dsigmajduy = 0.0_dp
       dsigmajduz = 0.0_dp
       dsjdx = 0.0_dp
       dsjdy = 0.0_dp
       dsjdz = 0.0_dp
       dsjdux = 0.0_dp
       dsjduy = 0.0_dp
       dsjduz = 0.0_dp
       depsjdx = 0.0_dp
       depsjdy = 0.0_dp
       depsjdz = 0.0_dp
       depsjdux = 0.0_dp
       depsjduy = 0.0_dp
       depsjduz = 0.0_dp

       do lm = 1, ShapeMap%Shapes(st2)%nContactFuncs
          l = ShapeMap%Shapes(st2)%ContactFuncLValue(lm)
          m = ShapeMap%Shapes(st2)%ContactFuncMValue(lm)
          coeff = ShapeMap%Shapes(st2)%ContactFuncCoefficient(lm)
          function_type = ShapeMap%Shapes(st2)%ContactFunctionType(lm)

          if ((function_type .eq. SH_COS).or.(m.eq.0)) then
             Phunc = coeff * tm_j(m)
             dPhuncdX = coeff * dtm_j(m) * dcpjdx
             dPhuncdY = coeff * dtm_j(m) * dcpjdy
             dPhuncdZ = coeff * dtm_j(m) * dcpjdz
             dPhuncdUx = coeff * dtm_j(m) * dcpjdux
             dPhuncdUy = coeff * dtm_j(m) * dcpjduy
             dPhuncdUz = coeff * dtm_j(m) * dcpjduz
          else
             Phunc = coeff * spj * um_j(m-1)
             dPhuncdX = coeff * (spj * dum_j(m-1) * dcpjdx + dspjdx *um_j(m-1))
             dPhuncdY = coeff * (spj * dum_j(m-1) * dcpjdy + dspjdy *um_j(m-1))
             dPhuncdZ = coeff * (spj * dum_j(m-1) * dcpjdz + dspjdz *um_j(m-1))
             dPhuncdUx = coeff*(spj * dum_j(m-1)*dcpjdux + dspjdux *um_j(m-1))
             dPhuncdUy = coeff*(spj * dum_j(m-1)*dcpjduy + dspjduy *um_j(m-1))
             dPhuncdUz = coeff*(spj * dum_j(m-1)*dcpjduz + dspjduz *um_j(m-1))
          endif

          sigma_j = sigma_j + plm_j(m,l)*Phunc

          dsigmajdx = dsigmajdx + plm_j(m,l)*dPhuncdX + &
               Phunc * dlm_j(m,l) * dctjdx
          dsigmajdy = dsigmajdy + plm_j(m,l)*dPhuncdY + &
               Phunc * dlm_j(m,l) * dctjdy
          dsigmajdz = dsigmajdz + plm_j(m,l)*dPhuncdZ + &
               Phunc * dlm_j(m,l) * dctjdz

          dsigmajdux = dsigmajdux + plm_j(m,l)* dPhuncdUx + &
               Phunc * dlm_j(m,l) * dctjdux
          dsigmajduy = dsigmajduy + plm_j(m,l)* dPhuncdUy + &
               Phunc * dlm_j(m,l) * dctjduy
          dsigmajduz = dsigmajduz + plm_j(m,l)* dPhuncdUz + &
               Phunc * dlm_j(m,l) * dctjduz

       end do

       do lm = 1, ShapeMap%Shapes(st2)%nRangeFuncs
          l = ShapeMap%Shapes(st2)%RangeFuncLValue(lm)
          m = ShapeMap%Shapes(st2)%RangeFuncMValue(lm)
          coeff = ShapeMap%Shapes(st2)%RangeFuncCoefficient(lm)
          function_type = ShapeMap%Shapes(st2)%RangeFunctionType(lm)

          if ((function_type .eq. SH_COS).or.(m.eq.0)) then
             Phunc = coeff * tm_j(m)
             dPhuncdX = coeff * dtm_j(m) * dcpjdx
             dPhuncdY = coeff * dtm_j(m) * dcpjdy
             dPhuncdZ = coeff * dtm_j(m) * dcpjdz
             dPhuncdUx = coeff * dtm_j(m) * dcpjdux
             dPhuncdUy = coeff * dtm_j(m) * dcpjduy
             dPhuncdUz = coeff * dtm_j(m) * dcpjduz
          else
             Phunc = coeff * spj * um_j(m-1)
             dPhuncdX = coeff * (spj * dum_j(m-1) * dcpjdx + dspjdx *um_j(m-1))
             dPhuncdY = coeff * (spj * dum_j(m-1) * dcpjdy + dspjdy *um_j(m-1))
             dPhuncdZ = coeff * (spj * dum_j(m-1) * dcpjdz + dspjdz *um_j(m-1))
             dPhuncdUx = coeff*(spj * dum_j(m-1)*dcpjdux + dspjdux *um_j(m-1))
             dPhuncdUy = coeff*(spj * dum_j(m-1)*dcpjduy + dspjduy *um_j(m-1))
             dPhuncdUz = coeff*(spj * dum_j(m-1)*dcpjduz + dspjduz *um_j(m-1))
          endif

          s_j = s_j + plm_j(m,l)*Phunc

          dsjdx = dsjdx + plm_j(m,l)*dPhuncdX + &
               Phunc * dlm_j(m,l) * dctjdx
          dsjdy = dsjdy + plm_j(m,l)*dPhuncdY + &
               Phunc * dlm_j(m,l) * dctjdy
          dsjdz = dsjdz + plm_j(m,l)*dPhuncdZ + &
               Phunc * dlm_j(m,l) * dctjdz

          dsjdux = dsjdux + plm_j(m,l)* dPhuncdUx + &
               Phunc * dlm_j(m,l) * dctjdux
          dsjduy = dsjduy + plm_j(m,l)* dPhuncdUy + &
               Phunc * dlm_j(m,l) * dctjduy
          dsjduz = dsjduz + plm_j(m,l)* dPhuncdUz + &
               Phunc * dlm_j(m,l) * dctjduz

       end do

       do lm = 1, ShapeMap%Shapes(st2)%nStrengthFuncs
          l = ShapeMap%Shapes(st2)%StrengthFuncLValue(lm)
          m = ShapeMap%Shapes(st2)%StrengthFuncMValue(lm)
          coeff = ShapeMap%Shapes(st2)%StrengthFuncCoefficient(lm)
          function_type = ShapeMap%Shapes(st2)%StrengthFunctionType(lm)

          if ((function_type .eq. SH_COS).or.(m.eq.0)) then
             Phunc = coeff * tm_j(m)
             dPhuncdX = coeff * dtm_j(m) * dcpjdx
             dPhuncdY = coeff * dtm_j(m) * dcpjdy
             dPhuncdZ = coeff * dtm_j(m) * dcpjdz
             dPhuncdUz = coeff * dtm_j(m) * dcpjdux
             dPhuncdUy = coeff * dtm_j(m) * dcpjduy
             dPhuncdUz = coeff * dtm_j(m) * dcpjduz
          else
             Phunc = coeff * spj * um_j(m-1)
             dPhuncdX = coeff * (spj * dum_j(m-1) * dcpjdx + dspjdx *um_j(m-1))
             dPhuncdY = coeff * (spj * dum_j(m-1) * dcpjdy + dspjdy *um_j(m-1))
             dPhuncdZ = coeff * (spj * dum_j(m-1) * dcpjdz + dspjdz *um_j(m-1))
             dPhuncdUx = coeff*(spj * dum_j(m-1)*dcpjdux + dspjdux *um_j(m-1))
             dPhuncdUy = coeff*(spj * dum_j(m-1)*dcpjduy + dspjduy *um_j(m-1))
             dPhuncdUz = coeff*(spj * dum_j(m-1)*dcpjduz + dspjduz *um_j(m-1))
          endif

!          write(*,*) 'l,m = ', l, m, coeff, dPhuncdUx, dPhuncdUy, dPhuncdUz

          eps_j = eps_j + plm_j(m,l)*Phunc

          depsjdx = depsjdx + plm_j(m,l)*dPhuncdX + &
               Phunc * dlm_j(m,l) * dctjdx
          depsjdy = depsjdy + plm_j(m,l)*dPhuncdY + &
               Phunc * dlm_j(m,l) * dctjdy
          depsjdz = depsjdz + plm_j(m,l)*dPhuncdZ + &
               Phunc * dlm_j(m,l) * dctjdz

          depsjdux = depsjdux + plm_j(m,l)* dPhuncdUx + &
               Phunc * dlm_j(m,l) * dctjdux
          depsjduy = depsjduy + plm_j(m,l)* dPhuncdUy + &
               Phunc * dlm_j(m,l) * dctjduy
          depsjduz = depsjduz + plm_j(m,l)* dPhuncdUz + &
               Phunc * dlm_j(m,l) * dctjduz

       end do

    endif

    ! phew, now let's assemble the potential energy:

    sigma = 0.5*(sigma_i + sigma_j)
!    write(*,*) sigma_i, ' = sigma_i; ', sigma_j, ' = sigma_j'
    dsigmadxi = 0.5*dsigmaidx
    dsigmadyi = 0.5*dsigmaidy
    dsigmadzi = 0.5*dsigmaidz
    dsigmaduxi = 0.5*dsigmaidux
    dsigmaduyi = 0.5*dsigmaiduy
    dsigmaduzi = 0.5*dsigmaiduz

    dsigmadxj = 0.5*dsigmajdx
    dsigmadyj = 0.5*dsigmajdy
    dsigmadzj = 0.5*dsigmajdz
    dsigmaduxj = 0.5*dsigmajdux
    dsigmaduyj = 0.5*dsigmajduy
    dsigmaduzj = 0.5*dsigmajduz

    s = 0.5*(s_i + s_j)

    dsdxi = 0.5*dsidx
    dsdyi = 0.5*dsidy
    dsdzi = 0.5*dsidz
    dsduxi = 0.5*dsidux
    dsduyi = 0.5*dsiduy
    dsduzi = 0.5*dsiduz

    dsdxj = 0.5*dsjdx
    dsdyj = 0.5*dsjdy
    dsdzj = 0.5*dsjdz
    dsduxj = 0.5*dsjdux
    dsduyj = 0.5*dsjduy
    dsduzj = 0.5*dsjduz

    eps = sqrt(eps_i * eps_j)
!!$    write(*,*) 'dsidu = ', dsidux, dsiduy, dsiduz
!!$    write(*,*) 'dsigidu = ', dsigmaidux, dsigmaiduy, dsigmaiduz
!!$    write(*,*) sigma_j, ' is sigma j; ', s_j, ' is s j; ', eps_j, ' is eps j'
    depsdxi = eps_j * depsidx / (2.0_dp * eps)
    depsdyi = eps_j * depsidy / (2.0_dp * eps)
    depsdzi = eps_j * depsidz / (2.0_dp * eps)
    depsduxi = eps_j * depsidux / (2.0_dp * eps)
    depsduyi = eps_j * depsiduy / (2.0_dp * eps)
    depsduzi = eps_j * depsiduz / (2.0_dp * eps)

    depsdxj = eps_i * depsjdx / (2.0_dp * eps)
    depsdyj = eps_i * depsjdy / (2.0_dp * eps)
    depsdzj = eps_i * depsjdz / (2.0_dp * eps)
    depsduxj = eps_i * depsjdux / (2.0_dp * eps)
    depsduyj = eps_i * depsjduy / (2.0_dp * eps)
    depsduzj = eps_i * depsjduz / (2.0_dp * eps)

!!$    write(*,*) 'depsidu = ', depsidux, depsiduy, depsiduz

!!$    write(*,*) 'depsjdu = ', depsjdux, depsjduy, depsjduz
!!$    write(*,*) 'depsduj = ', depsduxj, depsduyj, depsduzj
!!$
!!$    write(*,*) 's, sig, eps = ', s, sigma, eps

    rtdenom = rij-sigma+s
    rt = s / rtdenom

    drtdxi = (dsdxi - rt * (drdxi - dsigmadxi + dsdxi)) / rtdenom
    drtdyi = (dsdyi - rt * (drdyi - dsigmadyi + dsdyi)) / rtdenom
    drtdzi = (dsdzi - rt * (drdzi - dsigmadzi + dsdzi)) / rtdenom
    drtduxi = (dsduxi - rt * (drduxi - dsigmaduxi + dsduxi)) / rtdenom
    drtduyi = (dsduyi - rt * (drduyi - dsigmaduyi + dsduyi)) / rtdenom
    drtduzi = (dsduzi - rt * (drduzi - dsigmaduzi + dsduzi)) / rtdenom
    drtdxj = (dsdxj - rt * (drdxj - dsigmadxj + dsdxj)) / rtdenom
    drtdyj = (dsdyj - rt * (drdyj - dsigmadyj + dsdyj)) / rtdenom
    drtdzj = (dsdzj - rt * (drdzj - dsigmadzj + dsdzj)) / rtdenom
    drtduxj = (dsduxj - rt * (drduxj - dsigmaduxj + dsduxj)) / rtdenom
    drtduyj = (dsduyj - rt * (drduyj - dsigmaduyj + dsduyj)) / rtdenom
    drtduzj = (dsduzj - rt * (drduzj - dsigmaduzj + dsduzj)) / rtdenom

!!$    write(*,*) 'drtd_i = ', drtdxi, drtdyi, drtdzi
!!$    write(*,*) 'drtdu_j = ', drtduxj, drtduyj, drtduzj

    rt2 = rt*rt
    rt3 = rt2*rt
    rt5 = rt2*rt3
    rt6 = rt3*rt3
    rt11 = rt5*rt6
    rt12 = rt6*rt6
    rt126 = rt12 - rt6

    pot_temp = 4.0_dp * eps * rt126

    vpair = vpair + pot_temp

    pot = pot + pot_temp*sw

!!$    write(*,*) 'drtdu, depsdu = ', drtduxi, depsduxi

    dvdxi = 24.0_dp*eps*(2.0_dp*rt11 - rt5)*drtdxi + 4.0_dp*depsdxi*rt126
    dvdyi = 24.0_dp*eps*(2.0_dp*rt11 - rt5)*drtdyi + 4.0_dp*depsdyi*rt126
    dvdzi = 24.0_dp*eps*(2.0_dp*rt11 - rt5)*drtdzi + 4.0_dp*depsdzi*rt126
    dvduxi = 24.0_dp*eps*(2.0_dp*rt11 - rt5)*drtduxi + 4.0_dp*depsduxi*rt126
    dvduyi = 24.0_dp*eps*(2.0_dp*rt11 - rt5)*drtduyi + 4.0_dp*depsduyi*rt126
    dvduzi = 24.0_dp*eps*(2.0_dp*rt11 - rt5)*drtduzi + 4.0_dp*depsduzi*rt126

    dvdxj = 24.0_dp*eps*(2.0_dp*rt11 - rt5)*drtdxj + 4.0_dp*depsdxj*rt126
    dvdyj = 24.0_dp*eps*(2.0_dp*rt11 - rt5)*drtdyj + 4.0_dp*depsdyj*rt126
    dvdzj = 24.0_dp*eps*(2.0_dp*rt11 - rt5)*drtdzj + 4.0_dp*depsdzj*rt126
    dvduxj = 24.0_dp*eps*(2.0_dp*rt11 - rt5)*drtduxj + 4.0_dp*depsduxj*rt126
    dvduyj = 24.0_dp*eps*(2.0_dp*rt11 - rt5)*drtduyj + 4.0_dp*depsduyj*rt126
    dvduzj = 24.0_dp*eps*(2.0_dp*rt11 - rt5)*drtduzj + 4.0_dp*depsduzj*rt126
!!$    write(*,*) 'drtduxi = ', drtduxi, ' depsduxi = ', depsduxi
    ! do the torques first since they are easy:
    ! remember that these are still in the body fixed axes

    txi = 0.0_dp
    tyi = 0.0_dp
    tzi = 0.0_dp

    txj = 0.0_dp
    tyj = 0.0_dp
    tzj = 0.0_dp

    txi = (dvduyi - dvduzi) * sw
    tyi = (dvduzi - dvduxi) * sw
    tzi = (dvduxi - dvduyi) * sw

    txj = (dvduyj - dvduzj) * sw
    tyj = (dvduzj - dvduxj) * sw
    tzj = (dvduxj - dvduyj) * sw

!!$    txi = dvduxi * sw
!!$    tyi = dvduyi * sw
!!$    tzi = dvduzi * sw
!!$
!!$    txj = dvduxj * sw
!!$    tyj = dvduyj * sw
!!$    tzj = dvduzj * sw

    write(*,*) 't1 = ', txi, tyi, tzi
    write(*,*) 't2 = ', txj, tyj, tzj

    ! go back to lab frame using transpose of rotation matrix:

    t1(1) = t1(1) + a1(1)*txi + a1(4)*tyi + a1(7)*tzi
    t1(2) = t1(2) + a1(2)*txi + a1(5)*tyi + a1(8)*tzi
    t1(3) = t1(3) + a1(3)*txi + a1(6)*tyi + a1(9)*tzi

    t2(1) = t2(1) + a2(1)*txj + a2(4)*tyj + a2(7)*tzj
    t2(2) = t2(2) + a2(2)*txj + a2(5)*tyj + a2(8)*tzj
    t2(3) = t2(3) + a2(3)*txj + a2(6)*tyj + a2(9)*tzj

    ! Now, on to the forces:

    ! first rotate the i terms back into the lab frame:

    fxi = -dvdxi * sw
    fyi = -dvdyi * sw
    fzi = -dvdzi * sw

    fxj = -dvdxj * sw
    fyj = -dvdyj * sw
    fzj = -dvdzj * sw

    fxii = a1(1)*fxi + a1(4)*fyi + a1(7)*fzi
    fyii = a1(2)*fxi + a1(5)*fyi + a1(8)*fzi
    fzii = a1(3)*fxi + a1(6)*fyi + a1(9)*fzi

    fxjj = a2(1)*fxj + a2(4)*fyj + a2(7)*fzj
    fyjj = a2(2)*fxj + a2(5)*fyj + a2(8)*fzj
    fzjj = a2(3)*fxj + a2(6)*fyj + a2(9)*fzj

    fxij = -fxii
    fyij = -fyii
    fzij = -fzii

    fxji = -fxjj
    fyji = -fyjj
    fzji = -fzjj

    fxradial = (fxii + fxji)
    fyradial = (fyii + fyji)
    fzradial = (fzii + fzji)
!!$    write(*,*) fxradial, ' is fxrad; ', fyradial, ' is fyrad; ', fzradial, 'is fzrad'

    f1(1) = f1(1) + fxradial
    f1(2) = f1(2) + fxradial
    f1(3) = f1(3) + fxradial

  end subroutine do_shape_pair

  SUBROUTINE Associated_Legendre(x, l, m, lmax, plm, dlm)        

    ! Purpose: Compute the associated Legendre functions 
    !          Plm(x) and their derivatives Plm'(x)
    ! Input :  x  --- Argument of Plm(x)
    !          l  --- Order of Plm(x),  l = 0,1,2,...,n
    !          m  --- Degree of Plm(x), m = 0,1,2,...,N
    !          lmax --- Physical dimension of PLM and DLM
    ! Output:  PLM(l,m) --- Plm(x)
    !          DLM(l,m) --- Plm'(x)
    !
    ! adapted from the routines in 
    ! COMPUTATION OF SPECIAL FUNCTIONS by Shanjie Zhang and Jianming Jin
    ! ISBN 0-471-11963-6
    !
    ! The original Fortran77 codes can be found here:
    ! http://iris-lee3.ece.uiuc.edu/~jjin/routines/routines.html

    real (kind=dp), intent(in) :: x
    integer, intent(in) :: l, m, lmax 
    real (kind=dp), dimension(0:lmax,0:m), intent(out) :: PLM, DLM
    integer :: i, j, ls
    real (kind=dp) :: xq, xs

    ! zero out both arrays:
    DO I = 0, m
       DO J = 0, l
          PLM(J,I) = 0.0_dp
          DLM(J,I) = 0.0_dp
       end DO
    end DO

    ! start with 0,0:
    PLM(0,0) = 1.0_DP

    ! x = +/- 1 functions are easy:
    IF (abs(X).EQ.1.0_DP) THEN
       DO I = 1, m
          PLM(0, I) = X**I
          DLM(0, I) = 0.5_DP*I*(I+1.0_DP)*X**(I+1)
       end DO
       DO J = 1, m
          DO I = 1, l
             IF (I.EQ.1) THEN
                DLM(I, J) = 1.0D+300
             ELSE IF (I.EQ.2) THEN
                DLM(I, J) = -0.25_DP*(J+2)*(J+1)*J*(J-1)*X**(J+1)
             ENDIF
          end DO
       end DO
       RETURN
    ENDIF

    LS = 1
    IF (abs(X).GT.1.0_DP) LS = -1
    XQ = sqrt(LS*(1.0_DP-X*X))
    XS = LS*(1.0_DP-X*X)

    DO I = 1, l
       PLM(I, I) = -LS*(2.0_DP*I-1.0_DP)*XQ*PLM(I-1, I-1)
    enddo

    DO I = 0, l
       PLM(I, I+1)=(2.0_DP*I+1.0_DP)*X*PLM(I, I)
    enddo

    DO I = 0, l
       DO J = I+2, m
          PLM(I, J)=((2.0_DP*J-1.0_DP)*X*PLM(I,J-1) - &
               (I+J-1.0_DP)*PLM(I,J-2))/(J-I)
       end DO
    end DO

    DLM(0, 0)=0.0_DP
    DO J = 1, m
       DLM(0, J)=LS*J*(PLM(0,J-1)-X*PLM(0,J))/XS
    end DO

    DO I = 1, l
       DO J = I, m
          DLM(I,J) = LS*I*X*PLM(I, J)/XS + (J+I)*(J-I+1.0_DP)/XQ*PLM(I-1, J)
       end DO
    end DO

    RETURN
  END SUBROUTINE Associated_Legendre


  subroutine Orthogonal_Polynomial(x, m, mmax, function_type, pl, dpl)

    ! Purpose: Compute orthogonal polynomials: Tn(x) or Un(x),
    !          or Ln(x) or Hn(x), and their derivatives
    ! Input :  function_type --- Function code
    !                 =1 for Chebyshev polynomial Tn(x)
    !                 =2 for Chebyshev polynomial Un(x)
    !                 =3 for Laguerre polynomial Ln(x)
    !                 =4 for Hermite polynomial Hn(x)
    !          n ---  Order of orthogonal polynomials
    !          x ---  Argument of orthogonal polynomials
    ! Output:  PL(n) --- Tn(x) or Un(x) or Ln(x) or Hn(x)
    !          DPL(n)--- Tn'(x) or Un'(x) or Ln'(x) or Hn'(x)
    !
    ! adapted from the routines in 
    ! COMPUTATION OF SPECIAL FUNCTIONS by Shanjie Zhang and Jianming Jin
    ! ISBN 0-471-11963-6
    !
    ! The original Fortran77 codes can be found here:
    ! http://iris-lee3.ece.uiuc.edu/~jjin/routines/routines.html

    real(kind=dp), intent(in) :: x
    integer, intent(in):: m, mmax
    integer, intent(in):: function_type
    real(kind=dp), dimension(0:mmax), intent(inout) :: pl, dpl

    real(kind=dp) :: a, b, c, y0, y1, dy0, dy1, yn, dyn
    integer :: k

    A = 2.0_DP
    B = 0.0_DP
    C = 1.0_DP
    Y0 = 1.0_DP
    Y1 = 2.0_DP*X
    DY0 = 0.0_DP
    DY1 = 2.0_DP
    PL(0) = 1.0_DP
    PL(1) = 2.0_DP*X
    DPL(0) = 0.0_DP
    DPL(1) = 2.0_DP
    IF (function_type.EQ.CHEBYSHEV_TN) THEN
       Y1 = X
       DY1 = 1.0_DP
       PL(1) = X
       DPL(1) = 1.0_DP
    ELSE IF (function_type.EQ.LAGUERRE) THEN
       Y1 = 1.0_DP-X
       DY1 = -1.0_DP
       PL(1) = 1.0_DP-X
       DPL(1) = -1.0_DP
    ENDIF
    DO K = 2, m
       IF (function_type.EQ.LAGUERRE) THEN
          A = -1.0_DP/K
          B = 2.0_DP+A
          C = 1.0_DP+A
       ELSE IF (function_type.EQ.HERMITE) THEN
          C = 2.0_DP*(K-1.0_DP)
       ENDIF
       YN = (A*X+B)*Y1-C*Y0
       DYN = A*Y1+(A*X+B)*DY1-C*DY0
       PL(K) = YN
       DPL(K) = DYN
       Y0 = Y1
       Y1 = YN
       DY0 = DY1
       DY1 = DYN
    end DO


    RETURN

  end subroutine Orthogonal_Polynomial

  subroutine deallocateShapes(this)
    type(Shape), pointer :: this

    if (associated( this%ContactFuncLValue)) then
       deallocate(this%ContactFuncLValue)
       this%ContactFuncLValue => null()
    end if

    if (associated( this%ContactFuncMValue)) then
       deallocate( this%ContactFuncMValue)
       this%ContactFuncMValue => null()
    end if
    if (associated( this%ContactFunctionType)) then
       deallocate(this%ContactFunctionType)
       this%ContactFunctionType => null()
    end if

    if (associated( this%ContactFuncCoefficient)) then 
       deallocate(this%ContactFuncCoefficient)
       this%ContactFuncCoefficient => null()
    end if

    if (associated( this%RangeFuncLValue)) then 
       deallocate(this%RangeFuncLValue)
       this%RangeFuncLValue => null()
    end if
    if (associated( this%RangeFuncMValue)) then 
       deallocate( this%RangeFuncMValue) 
       this%RangeFuncMValue => null()
    end if

    if (associated( this%RangeFunctionType)) then
       deallocate( this%RangeFunctionType) 
       this%RangeFunctionType => null()
    end if
    if (associated( this%RangeFuncCoefficient)) then
       deallocate(this%RangeFuncCoefficient) 
       this%RangeFuncCoefficient => null()
    end if

    if (associated( this%StrengthFuncLValue)) then
       deallocate(this%StrengthFuncLValue)
       this%StrengthFuncLValue => null()
    end if

    if (associated( this%StrengthFuncMValue )) then
       deallocate(this%StrengthFuncMValue)
       this%StrengthFuncMValue => null()
    end if

    if(associated( this%StrengthFunctionType)) then
       deallocate(this%StrengthFunctionType) 
       this%StrengthFunctionType => null()
    end if
    if (associated( this%StrengthFuncCoefficient )) then 
       deallocate(this%StrengthFuncCoefficient)
       this%StrengthFuncCoefficient => null()
    end if
  end subroutine deallocateShapes

  subroutine destroyShapeTypes
    integer :: i 
    type(Shape), pointer :: thisShape

    ! First walk through and kill the shape
    do i = 1,ShapeMap%n_shapes
       thisShape => ShapeMap%Shapes(i)
       call deallocateShapes(thisShape)
    end do

    ! set shape map to starting values
    ShapeMap%n_shapes = 0
    ShapeMap%currentShape = 0

    if (associated(ShapeMap%Shapes)) then
       deallocate(ShapeMap%Shapes)
       ShapeMap%Shapes => null()
    end if

    if (associated(ShapeMap%atidToShape)) then
       deallocate(ShapeMap%atidToShape)
       ShapeMap%atidToShape => null()
    end if


  end subroutine destroyShapeTypes


end module shapes
