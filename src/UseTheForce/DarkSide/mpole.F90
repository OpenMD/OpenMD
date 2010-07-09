subroutine doMultipolePair(atom1, atom2, d, rij, r2, rcut, sw, &
     vpair, fpair, pot, eFrame, f, t)

  !***************************************************************
  ! doMultipolePair evaluates the potential, forces, and torques 
  ! between point multipoles.  It is based on the ewald2 real
  ! space routine in the MDMULP code written by W. Smith and 
  ! available from the CCP5 website.
  !
  ! We're using the damped real space portion of the Ewald sum
  ! as the entire interaction.   Details on this portion of the
  ! sum can be found in:
  !
  ! "Point Multipoles in the Ewald Summation (Revisited)," by
  ! W. Smith, CCP5 Newsletter, 46, pp. 18-30 (1998).
  !
  !**************************************************************

  integer, intent(in) :: atom1, atom2
  real(kind=dp), intent(in) :: rij, r2, sw, rcut
  real(kind=dp), intent(in), dimension(3) :: d
  real(kind=dp), intent(inout) :: vpair
  real(kind=dp), intent(inout), dimension(3) :: fpair

  real( kind = dp ) :: pot
  real( kind = dp ), dimension(9,nLocal) :: eFrame
  real( kind = dp ), dimension(3,nLocal) :: f
  real( kind = dp ), dimension(3,nLocal) :: t
  logical :: i_is_Charge, i_is_Dipole, i_is_SplitDipole, i_is_Quadrupole
  logical :: j_is_Charge, j_is_Dipole, j_is_SplitDipole, j_is_Quadrupole
  integer :: me1, me2, id1, id2
  real (kind=dp) :: q_i, q_j, mu_i, mu_j, d_i, d_j
  real (kind=dp) :: qxx_i, qyy_i, qzz_i
  real (kind=dp) :: qxx_j, qyy_j, qzz_j

  real (kind=dp) :: epot, qii, alsq2, alsq2n, exp2a, ralpha 
  real (kind=dp), dimension(3) :: ftm, ttm_i, ttm_j
  real (kind=dp), dimension(3) :: qir, qjr, qidj, qjdi, qiqjr, qjqir
  real (kind=dp), dimension(3) :: dixdj, qixqj, dixr, djxr, rxqir, rxqjr
  real (kind=dp), dimension(3) :: dixqjr, djxqir, rxqidj, rxqjdi, rxqijr
  real (kind=dp), dimension(3) :: rxqjir, qjrxqir
  real (kind=dp), dimension(6) :: bn
  real (kind=dp), dimension(10) :: sc
  real (kind=dp), dimension(9) :: gl

  real(kind=dp) :: a1, a2, a3, a4, a5, p

  DATA A1,A2,A3/0.254829592,-0.284496736,1.421413741/
  DATA A4,A5,P/-1.453152027,1.061405429,0.3275911/

  if (.not.summationMethodChecked) then
     call checkSummationMethod()
  endif

#ifdef IS_MPI
  me1 = atid_Row(atom1)
  me2 = atid_Col(atom2)
#else
  me1 = atid(atom1)
  me2 = atid(atom2)
#endif

  ! logicals
  i_is_Charge = ElectrostaticMap(me1)%is_Charge
  i_is_Dipole = ElectrostaticMap(me1)%is_Dipole
  i_is_Quadrupole = ElectrostaticMap(me1)%is_Quadrupole

  j_is_Charge = ElectrostaticMap(me2)%is_Charge
  j_is_Dipole = ElectrostaticMap(me2)%is_Dipole
  j_is_Quadrupole = ElectrostaticMap(me2)%is_Quadrupole

  if (i_is_Charge) then
     q_i = ElectrostaticMap(me1)%charge
  else
     q_i = 0.0_dp
  endif

  if (j_is_Charge) then
     q_j = ElectrostaticMap(me2)%charge
  else
     q_j = 0.0_dp
  endif

  if (i_is_Dipole) then
     mu_i = ElectrostaticMap(me1)%dipole_moment
#ifdef IS_MPI
     d_i(1) = eFrame_Row(3,atom1)
     d_i(2) = eFrame_Row(6,atom1)
     d_i(3) = eFrame_Row(9,atom1)
#else
     d_i(1) = eFrame(3,atom1)
     d_i(2) = eFrame(6,atom1)
     d_i(3) = eFrame(9,atom1)
#endif
     d_i = d_i * mu_i
  else
     d_i = 0.0_dp
  endif

  if (j_is_Dipole) then
     mu_j = ElectrostaticMap(me2)%dipole_moment
#ifdef IS_MPI
     d_j(1) = eFrame_Col(3,atom2)
     d_j(2) = eFrame_Col(6,atom2)
     d_j(3) = eFrame_Col(9,atom2)
#else
     d_j(1) = eFrame(3,atom2)
     d_j(2) = eFrame(6,atom2)
     d_j(3) = eFrame(9,atom2)
#endif
     d_j = d_j * mu_j
  else
     d_j = 0.0_dp
  endif

  if (i_is_Quadrupole) then
     qxx_i = ElectrostaticMap(me1)%quadrupole_moments(1)
     qyy_i = ElectrostaticMap(me1)%quadrupole_moments(2)
     qzz_i = ElectrostaticMap(me1)%quadrupole_moments(3)
#ifdef IS_MPI
     qpole_i(1) = qxx_i * eFrame_Row(1,atom1)
     qpole_i(2) = qxx_i * eFrame_Row(4,atom1)
     qpole_i(3) = qxx_i * eFrame_Row(7,atom1)
     qpole_i(4) = qyy_i * eFrame_Row(2,atom1)
     qpole_i(5) = qyy_i * eFrame_Row(5,atom1)
     qpole_i(6) = qyy_i * eFrame_Row(8,atom1)
     qpole_i(7) = qzz_i * eFrame_Row(3,atom1)
     qpole_i(8) = qzz_i * eFrame_Row(6,atom1)
     qpole_i(9) = qzz_i * eFrame_Row(9,atom1)
#else    
     qpole_i(1) = qxx_i * eFrame(1,atom1)
     qpole_i(2) = qxx_i * eFrame(4,atom1)
     qpole_i(3) = qxx_i * eFrame(7,atom1)
     qpole_i(4) = qyy_i * eFrame(2,atom1)
     qpole_i(5) = qyy_i * eFrame(5,atom1)
     qpole_i(6) = qyy_i * eFrame(8,atom1)
     qpole_i(7) = qzz_i * eFrame(3,atom1)
     qpole_i(8) = qzz_i * eFrame(6,atom1)
     qpole_i(9) = qzz_i * eFrame(9,atom1)
#endif
  else
     qpole_i = 0.0_dp
  endif

  if (j_is_Quadrupole) then
     qxx_j = ElectrostaticMap(me2)%quadrupole_moments(1)
     qyy_j = ElectrostaticMap(me2)%quadrupole_moments(2)
     qzz_j = ElectrostaticMap(me2)%quadrupole_moments(3)
#ifdef IS_MPI
     qpole_j(1) = qxx_j * eFrame_Col(1,atom2)
     qpole_j(2) = qxx_j * eFrame_Col(4,atom2)
     qpole_j(3) = qxx_j * eFrame_Col(7,atom2)
     qpole_j(4) = qyy_j * eFrame_Col(2,atom2)
     qpole_j(5) = qyy_j * eFrame_Col(5,atom2)
     qpole_j(6) = qyy_j * eFrame_Col(8,atom2)
     qpole_j(7) = qzz_j * eFrame_Col(3,atom2)
     qpole_j(8) = qzz_j * eFrame_Col(6,atom2)
     qpole_j(9) = qzz_j * eFrame_Col(9,atom2)
#else    
     qpole_j(1) = qxx_j * eFrame(1,atom2)
     qpole_j(2) = qxx_j * eFrame(4,atom2)
     qpole_j(3) = qxx_j * eFrame(7,atom2)
     qpole_j(4) = qyy_j * eFrame(2,atom2)
     qpole_j(5) = qyy_j * eFrame(5,atom2)
     qpole_j(6) = qyy_j * eFrame(8,atom2)
     qpole_j(7) = qzz_j * eFrame(3,atom2)
     qpole_j(8) = qzz_j * eFrame(6,atom2)
     qpole_j(9) = qzz_j * eFrame(9,atom2)
#endif
  else
     qpole_j = 0.0_dp
  endif

  !
  !     INITIALISE TEMPORARY FORCE AND TORQUE ACCUMULATORS
  !

  ftm = 0.0
  ttm_i = 0.0
  ttm_j = 0.0

  ri = 1.0_dp / rij
  ri2 = ri * ri

  !
  !     CALCULATE BN FUNCTIONS
  !
  if (screeningMethod .eq. DAMPED) then
     ! assemble the damping variables
     call lookupUniformSpline1d(erfcSpline, rij, erfcVal, derfcVal)
  else
     dampingAlpha = 0.0_dp
     erfcVal = 1.0_dp
  endif

  rrtpi = 1.0 / sqrt(pi)
  alsq2 = 2.0 * dampingAlpha**2
  ralpi = 0.0
  if(dampingAlpha .gt. 0.0) ralpi = rrtpi / dampingAlpha

  alphar = dampingAlpha * rij
  t = 1.0 / (1.0 + p*alphar)
  exp2a = exp(-alphar**2)
  bn(1) = ((((a5*t+a4)*t+a3)*t+a2)*t+a1) * t * exp2a * ri
  alsq2n = ralpi

  do n = 1, 5
     bfac = real(n+n-1, kind=dp)
     alsq2n = alsq2*alsq2n     
     bn(n+1) = ri2 * (bfac*bn(n) + alsq2n*exp2a)
  enddo

  ralpha = dampingAlpha * rij
  bn(0) = erfcVal * ri
  alsq2 = 2.0d0 * dampingAlpha**2
  alsq2n = 0.0d0
  if (dampingAlpha .gt. 0.0_dp)  alsq2n = 1.0_dp / (sqrtpi*dampingAlpha)
  exp2a = exp(-ralpha**2)
  do m = 1, 5
     bfac = real(m+m-1, kind=dp)
     alsq2n = alsq2 * alsq2n
     bn(m) = (bfac*bn(m-1)+alsq2n*exp2a) / r2
  end do

  !
  !     CONSTRUCT AUXILIARY VECTORS
  !

  dixdj(1) = d_i(2)*d_j(3) - d_i(3)*d_j(2)
  dixdj(2) = d_i(3)*d_j(1) - d_i(1)*d_j(3)
  dixdj(3) = d_i(1)*d_j(2) - d_i(2)*d_j(1)

  dixr(1) = d_i(2)*d(3) - d_i(3)*d(2)
  dixr(2) = d_i(3)*d(1) - d_i(1)*d(3)
  dixr(3) = d_i(1)*d(2) - d_i(2)*d(1)

  djxr(1) = d_j(2)*d(3) - d_j(3)*d(2)
  djxr(2) = d_j(3)*d(1) - d_j(1)*d(3)
  djxr(3) = d_j(1)*d(2) - d_j(2)*d(1)

  qir(1) = qpole_i(1)*d(1) + qpole_i(4)*d(2) + qpole_i(7)*d(3)
  qir(2) = qpole_i(2)*d(1) + qpole_i(5)*d(2) + qpole_i(8)*d(3)
  qir(3) = qpole_i(3)*d(1) + qpole_i(6)*d(2) + qpole_i(9)*d(3)

  qjr(1) = qpole_j(1)*d(1) + qpole_j(4)*d(2) + qpole_j(7)*d(3)
  qjr(2) = qpole_j(2)*d(1) + qpole_j(5)*d(2) + qpole_j(8)*d(3)
  qjr(3) = qpole_j(3)*d(1) + qpole_j(6)*d(2) + qpole_j(9)*d(3)

  qiqjr(1) = qpole_i(1)*qjr(1) + qpole_i(4)*qjr(2) + qpole_i(7)*qjr(3)
  qiqjr(2) = qpole_i(2)*qjr(1) + qpole_i(5)*qjr(2) + qpole_i(8)*qjr(3)
  qiqjr(3) = qpole_i(3)*qjr(1) + qpole_i(6)*qjr(2) + qpole_i(9)*qjr(3)


  qjqir(1) = qpole_j(1)*QIR(1) + qpole_j(4)*QIR(2) + qpole_j(7)*QIR(3)
  qjqir(2) = qpole_j(2)*QIR(1) + qpole_j(5)*QIR(2) + qpole_j(8)*QIR(3)
  qjqir(3) = qpole_j(3)*QIR(1) + qpole_j(6)*QIR(2) + qpole_j(9)*QIR(3)

  qixqj(1) = qpole_i(2)*qpole_j(3) + qpole_i(5)*qpole_j(6) + &
       qpole_i(8)*qpole_j(9) - qpole_i(3)*qpole_j(2) - &
       qpole_i(6)*qpole_j(5) - qpole_i(9)*qpole_j(8)

  qixqj(2) = qpole_i(3)*qpole_j(1) + qpole_i(6)*qpole_j(4) + &
       qpole_i(9)*qpole_j(7) - qpole_i(1)*qpole_j(3) - &
       qpole_i(4)*qpole_j(6) - qpole_i(7)*qpole_j(9)

  qixqj(3) = qpole_i(1)*qpole_j(2) + qpole_i(4)*qpole_j(5) + &
       qpole_i(7)*qpole_j(8) - qpole_i(2)*qpole_j(1) - &
       qpole_i(5)*qpole_j(4) - qpole_i(8)*qpole_j(7)

  rxqir(1) = d(2)*qir(3) - d(3)*qir(2)
  rxqir(2) = d(3)*qir(1) - d(1)*qir(3)
  rxqir(3) = d(1)*qir(2) - d(2)*qir(1)

  rxqjr(1) = d(2)*qjr(3) - d(3)*qjr(2)
  rxqjr(2) = d(3)*qjr(1) - d(1)*qjr(3)
  rxqjr(3) = d(1)*qjr(2) - d(2)*qjr(1)

  rxqijr(1) = d(2)*qiqjr(3) - d(3)*qiqjr(2)
  rxqijr(2) = d(3)*qiqjr(1) - d(1)*qiqjr(3)
  rxqijr(3) = d(1)*qiqjr(2) - d(2)*qiqjr(1)

  rxqjir(1) = d(2)*qjqir(3) - d(3)*qjqir(2)
  rxqjir(2) = d(3)*qjqir(1) - d(1)*qjqir(3)
  rxqjir(3) = d(1)*qjqir(2) - d(2)*qjqir(1)

  qjrxqir(1) = qjr(2)*qir(3) - qjr(3)*qir(2)
  qjrxqir(2) = qjr(3)*qir(1) - qjr(1)*qir(3)
  qjrxqir(3) = qjr(1)*qir(2) - qjr(2)*qir(1)

  qidj(1) = qpole_i(1)*d_j(1) + qpole_i(4)*d_j(2) + qpole_i(7)*d_j(3)
  qidj(2) = qpole_i(2)*d_j(1) + qpole_i(5)*d_j(2) + qpole_i(8)*d_j(3)
  qidj(3) = qpole_i(3)*d_j(1) + qpole_i(6)*d_j(2) + qpole_i(9)*d_j(3)

  qjdi(1) = qpole_j(1)*d_i(1) + qpole_j(4)*d_i(2) + qpole_j(7)*d_i(3)
  qjdi(2) = qpole_j(2)*d_i(1) + qpole_j(5)*d_i(2) + qpole_j(8)*d_i(3)
  qjdi(3) = qpole_j(3)*d_i(1) + qpole_j(6)*d_i(2) + qpole_j(9)*d_i(3)

  dixqjr(1) = d_i(2)*qjr(3) - d_i(3)*qjr(2)
  dixqjr(2) = d_i(3)*qjr(1) - d_i(1)*qjr(3)
  dixqjr(3) = d_i(1)*qjr(2) - d_i(2)*qjr(1)

  djxqir(1) = d_j(2)*qir(3) - d_j(3)*qir(2)
  djxqir(2) = d_j(3)*qir(1) - d_j(1)*qir(3)
  djxqir(3) = d_j(1)*qir(2) - d_j(2)*qir(1)

  rxqidj(1) = d(2)*qidj(3) - d(3)*qidj(2)
  rxqidj(2) = d(3)*qidj(1) - d(1)*qidj(3)
  rxqidj(3) = d(1)*qidj(2) - d(2)*qidj(1)

  rxqjdi(1) = d(2)*qjdi(3) - d(3)*qjdi(2)
  rxqjdi(2) = d(3)*qjdi(1) - d(1)*qjdi(3)
  rxqjdi(3) = d(1)*qjdi(2) - d(2)*qjdi(1)

  !
  !     CALCULATE SCALAR PRODUCTS
  !

  qii = qpole_i(1) + qpole_i(5) + qpole_i(9)

  sc(1) = qpole_j(1) + qpole_j(5) + qpole_j(9)
  sc(2) = d_i(1)*d_j(1) + d_i(2)*d_j(2) + d_i(3)*d_j(3)
  sc(3) = d_i(1)*d(1) + d_i(2)*d(2) + d_i(3)*d(3)
  sc(4) = d_j(1)*d(1) + d_j(2)*d(2) + d_j(3)*d(3)
  sc(5) = qir(1)*d(1) + qir(2)*d(2) + qir(3)*d(3)
  sc(6) = qjr(1)*d(1) + qjr(2)*d(2) + qjr(3)*d(3)
  sc(7) = qir(1)*d_j(1) + qir(2)*d_j(2) + qir(3)*d_j(3)
  sc(8) = qjr(1)*d_i(1) + qjr(2)*d_i(2) + qjr(3)*d_i(3)
  sc(9) = qir(1)*qjr(1) + qir(2)*qjr(2) + qir(3)*qjr(3)
  sc(10) = qpole_i(1)*qpole_j(1) + qpole_i(2)*qpole_j(2) + &
       qpole_i(3)*qpole_j(3) + qpole_i(4)*qpole_j(4) + &
       qpole_i(5)*qpole_j(5) + qpole_i(6)*qpole_j(6) + &
       qpole_i(7)*qpole_j(7) + qpole_i(8)*qpole_j(8) + &
       qpole_i(9)*qpole_j(9)

  !
  !     CALCULATE GL FUNCTIONS FOR POTENTIAL
  !

  gl(1) =  q_i*q_j
  gl(2) =  q_j*sc(3) - q_i*sc(4)
  gl(3) =  q_i*sc(6) + q_j*sc(5) - sc(3)*sc(4)
  gl(4) =  sc(3)*sc(6) - sc(4)*sc(5)
  gl(5) =  sc(5)*sc(6)
  gl(7) =  sc(2) - q_j*qii - q_i*sc(1)
  gl(8) =  sc(4)*qii - sc(1)*sc(3) + 2.0*(sc(7) - sc(8))
  gl(9) =  2.0*sc(10) + sc(1)*qii
  gl(6) =  -(sc(1)*sc(5) + sc(6)*qii + 4.0*sc(9))

  !
  !     CALCULATE POTENTIAL AND VIRIAL
  !

  epot = bn(1)*gl(1) + bn(2)*(gl(2) + gl(7)) + bn(3)*(gl(3) &
       + gl(8) + gl(9)) + bn(4)*(gl(4) + gl(6)) + bn(5)*gl(5)

  vpair = vpair + epot
  
#ifdef IS_MPI
  pot_row(ELECTROSTATIC_POT,atom1) = pot_row(ELECTROSTATIC_POT,atom1) + &
       0.5_dp*epot*sw
  pot_col(ELECTROSTATIC_POT,atom2) = pot_col(ELECTROSTATIC_POT,atom2) + &
       0.5_dp*epot*sw
#else
  pot = pot + epot*sw
#endif

  !
  ! CALCULATE FORCE AND TORQUE COEFFICIENTS
  !

  gl(1) = bn(2)*gl(1) + bn(3)*(gl(2) + gl(7)) + bn(4)*(gl(3) &
       + gl(8) + gl(9)) + bn(5)*(gl(4) + gl(6)) + bn(6)*gl(5)
  gl(2) =  -q_j*bn(2) + (sc(4) + sc(1))*bn(3) - sc(6)*bn(4)
  gl(3) =  q_i*bn(2) + (sc(3) - qii)*bn(3) + sc(5)*bn(4)
  gl(4) =  2.0*bn(3)
  gl(5) = 2.0*( - q_j*bn(3) + (sc(1) + sc(4))*bn(4) - sc(6)* bn(5))
  gl(6) = 2.0*( - q_i*bn(3) + (qii - sc(3))*bn(4) - sc(5)*bn(5))
  gl(7) = 4.0*bn(4)

  !
  !     CALCULATE FORCES BETWEEN MULTIPOLES I AND J
  !

  ftm(1) = gl(1)*d(1) + gl(2)*d_i(1) + gl(3)*d_j(1) +  &
       gl(4)*(qjdi(1) - qidj(1)) + gl(5)*qir(1) +  &
       gl(6)*qjr(1) + gl(7)*(qiqjr(1) + qjqir(1))
  ftm(2) = gl(1)*d(2) + gl(2)*d_i(2) + gl(3)*d_j(2) +  &
       gl(4)*(qjdi(2) - qidj(2)) + gl(5)*qir(2) +  &
       gl(6)*qjr(2) + gl(7)*(qiqjr(2) + qjqir(2))
  ftm(3) = gl(1)*d(3) + gl(2)*d_i(3) + gl(3)*d_j(3) +  &
       gl(4)*(qjdi(3) - qidj(3)) + gl(5)*qir(3) +  &
       gl(6)*qjr(3) + gl(7)*(qiqjr(3) + qjqir(3))

  !
  !     CALCULATE TORQUES BETWEEN MULTIPOLES I AND J
  !

  ttm_i(1) =  - bn(2)*dixdj(1) + gl(2)*dixr(1) + gl(4)* &
       (dixqjr(1) + djxqir(1) + rxqidj(1) - 2.0*qixqj(1)) -  &
       gl(5)*rxqir(1) - gl(7)*(rxqijr(1) + qjrxqir(1))
  ttm_i(2) =  - bn(2)*dixdj(2) + gl(2)*dixr(2) + gl(4)* &
       (dixqjr(2) + djxqir(2) + rxqidj(2) - 2.0*qixqj(2)) -  &
       gl(5)*rxqir(2) - gl(7)*(rxqijr(2) + qjrxqir(2))
  ttm_i(3) =  - bn(2)*dixdj(3) + gl(2)*dixr(3) + gl(4)* &
       (dixqjr(3) + djxqir(3) + rxqidj(3) - 2.0*qixqj(3)) -  &
       gl(5)*rxqir(3) - gl(7)*(rxqijr(3) + qjrxqir(3))
  ttm_j(1) =  bn(2)*dixdj(1) + gl(3)*djxr(1) - gl(4)* &
       (dixqjr(1) + djxqir(1) + rxqjdi(1) - 2.0*qixqj(1)) -  &
       gl(6)*rxqjr(1) - gl(7)*(rxqjir(1) - qjrxqir(1))
  ttm_j(2) =  bn(2)*dixdj(2) + gl(3)*djxr(2) - gl(4)* &
       (dixqjr(2) + djxqir(2) + rxqjdi(2) - 2.0*qixqj(2)) -  &
       gl(6)*rxqjr(2) - gl(7)*(rxqjir(2) - qjrxqir(2))
  ttm_j(3) =  bn(2)*dixdj(3) + gl(3)*djxr(3) - gl(4)* &
       (dixqjr(3) + djxqir(3) + rxqjdi(3) - 2.0*qixqj(3)) -  &
       gl(6)*rxqjr(3) - gl(7)*(rxqjir(3) - qjrxqir(3))

  ftm = ftm*sw
  ttm_i = ttm_i*sw
  ttm_j = ttm_j*sw

#ifdef IS_MPI
  f_Row(1,atom1) = f_Row(1,atom1) + ftm(1)
  f_Row(2,atom1) = f_Row(2,atom1) + ftm(2)
  f_Row(3,atom1) = f_Row(3,atom1) + ftm(3)

  f_Col(1,atom2) = f_Col(1,atom2) - ftm(1)
  f_Col(2,atom2) = f_Col(2,atom2) - ftm(2)
  f_Col(3,atom2) = f_Col(3,atom2) - ftm(3)

  t_Row(1,atom1) = t_Row(1,atom1) + ttm_i(1)
  t_Row(2,atom1) = t_Row(2,atom1) + ttm_i(2)
  t_Row(3,atom1) = t_Row(3,atom1) + ttm_i(3)

  t_Col(1,atom2) = t_Col(1,atom2) + ttm_j(1)
  t_Col(2,atom2) = t_Col(2,atom2) + ttm_j(2)
  t_Col(3,atom2) = t_Col(3,atom2) + ttm_j(3)
#else
  f(1,atom1) = f(1,atom1) + ftm(1)
  f(2,atom1) = f(2,atom1) + ftm(2)
  f(3,atom1) = f(3,atom1) + ftm(3)

  f(1,atom2) = f(1,atom2) - ftm(1)
  f(2,atom2) = f(2,atom2) - ftm(2)
  f(3,atom2) = f(3,atom2) - ftm(3)

  t(1,atom1) = t(1,atom1) + ttm_i(1)
  t(2,atom1) = t(2,atom1) + ttm_i(2)
  t(3,atom1) = t(3,atom1) + ttm_i(3)

  t(1,atom2) = t(1,atom2) + ttm_j(1)
  t(2,atom2) = t(2,atom2) + ttm_j(2)
  t(3,atom2) = t(3,atom2) + ttm_j(3)
#endif

#ifdef IS_MPI
  id1 = AtomRowToGlobal(atom1)
  id2 = AtomColToGlobal(atom2)
#else
  id1 = atom1
  id2 = atom2
#endif

  if (molMembershipList(id1) .ne. molMembershipList(id2)) then

     fpair(1) = fpair(1) + ftm(1)
     fpair(2) = fpair(2) + ftm(2)
     fpair(3) = fpair(3) + ftm(3)

  endif

  return

end subroutine domultipolepair
