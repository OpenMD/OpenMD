subroutine do_electrostatic_new()
  
  ci = ElectrostaticMap(me1)%charge

  mui = ElectrostaticMap(me1)%dipole_moment
  Thetai = ElectrostaticMap(me1)%quadrupole_moments
  eFramei = eFrame(atom1)
  uz_i = eFramei(:,3)


  di = mui * uz_i
  qi = Thetai(:)*eFramei


  f = electric / dielec
  xji = x(jj) - x(ii)
  yji = y(jj) - y(ii)
  zji = z(jj) - z(ii)
  r2 = xji*xji + yji*yji + zji*zji
  r = sqrt(r2)
  rr1 = 1.0d0 / r
  rr3 = rr1 / r2
  rr5 = 3.0d0 * rr3 / r2
  rr7 = 5.0d0 * rr5 / r2
  rr9 = 7.0d0 * rr7 / r2
  rr11 = 9.0d0 * rr9 / r2
  
  !     construct auxiliary vectors
  
  qir(1) = qi(1)*xji + qi(4)*yji + qi(7)*zji
  qir(2) = qi(2)*xji + qi(5)*yji + qi(8)*zji
  qir(3) = qi(3)*xji + qi(6)*yji + qi(9)*zji
  qjr(1) = qj(1)*xji + qj(4)*yji + qj(7)*zji
  qjr(2) = qj(2)*xji + qj(5)*yji + qj(8)*zji
  qjr(3) = qj(3)*xji + qj(6)*yji + qj(9)*zji

  !     get intermediate variables for energy terms

  sc(2) = di(1)*dj(1) + di(2)*dj(2) + di(3)*dj(3)
  sc(3) = di(1)*xji + di(2)*yji + di(3)*zji
  sc(4) = dj(1)*xji + dj(2)*yji + dj(3)*zji
  sc(5) = qir(1)*xji + qir(2)*yji + qir(3)*zji
  sc(6) = qjr(1)*xji + qjr(2)*yji + qjr(3)*zji
  sc(7) = qir(1)*dj(1) + qir(2)*dj(2) + qir(3)*dj(3)
  sc(8) = qjr(1)*di(1) + qjr(2)*di(2) + qjr(3)*di(3)
  sc(9) = qir(1)*qjr(1) + qir(2)*qjr(2) + qir(3)*qjr(3)
  sc(10) = qi(1)*qj(1) + qi(2)*qj(2) + qi(3)*qj(3) &
       + qi(4)*qj(4) + qi(5)*qj(5) + qi(6)*qj(6) &
       + qi(7)*qj(7) + qi(8)*qj(8) + qi(9)*qj(9)

  !     calculate the gl functions for potential energy

  gl(0) = ci*cj
  gl(1) = cj*sc(3) - ci*sc(4)
  gl(2) = ci*sc(6) + cj*sc(5) - sc(3)*sc(4)
  gl(3) = sc(3)*sc(6) - sc(4)*sc(5)
  gl(4) = sc(5)*sc(6)
  gl(6) = sc(2)
  gl(7) = 2.0d0 * (sc(7)-sc(8))
  gl(8) = 2.0d0 * sc(10)
  gl(5) = -4.0d0 * sc(9)

  !    get the permanent multipole energy
 
  e = rr1*gl(0) + rr3*(gl(1)+gl(6)) &
       + rr5*(gl(2)+gl(7)+gl(8)) &
       + rr7*(gl(3)+gl(5)) + rr9*gl(4)

  !     intermediate variables for the multipole terms

  gf(1) = rr3*gl(0) + rr5*(gl(1)+gl(6)) &
       + rr7*(gl(2)+gl(7)+gl(8)) &
       + rr9*(gl(3)+gl(5)) + rr11*gl(4)
  gf(2) = -cj*rr3 + sc(4)*rr5 - sc(6)*rr7
  gf(3) = ci*rr3 + sc(3)*rr5 + sc(5)*rr7
  gf(4) = 2.0d0 * rr5
  gf(5) = 2.0d0 * (-cj*rr5+sc(4)*rr7-sc(6)*rr9)
  gf(6) = 2.0d0 * (-ci*rr5-sc(3)*rr7-sc(5)*rr9)
  gf(7) = 4.0d0 * rr7

  !     get the force

  ftm2(1) = gf(1)*xji + gf(2)*di(1) + gf(3)*dj(1)          &
       + gf(4)*(qjdi(1)-qidj(1)) + gf(5)*qir(1)     &
       + gf(6)*qjr(1) + gf(7)*(qiqjr(1)+qjqir(1))    
  ftm2(2) = gf(1)*yji + gf(2)*di(2) + gf(3)*dj(2)           &
       + gf(4)*(qjdi(2)-qidj(2)) + gf(5)*qir(2)     &
       + gf(6)*qjr(2) + gf(7)*(qiqjr(2)+qjqir(2))    
  ftm2(3) = gf(1)*zji + gf(2)*di(3) + gf(3)*dj(3)           &
       + gf(4)*(qjdi(3)-qidj(3)) + gf(5)*qir(3)     &
       + gf(6)*qjr(3) + gf(7)*(qiqjr(3)+qjqir(3))

  !     now perform the torque calculation

  !     get the interaction torques

  ttm2(1) = -rr3*dixdj(1) + gf(2)*dixr(1)-gf(5)*rxqir(1)  & 
       + gf(4)*(dixqjr(1)+djxqir(1)+rxqidj(1)-2.0d0*qixqj(1))&
       - gf(7)*(rxqijr(1)+qjrxqir(1))                         
  ttm2(2) = -rr3*dixdj(2) + gf(2)*dixr(2)-gf(5)*rxqir(2)  &
       + gf(4)*(dixqjr(2)+djxqir(2)+rxqidj(2)-2.0d0*qixqj(2))&
       - gf(7)*(rxqijr(2)+qjrxqir(2))                         
  ttm2(3) = -rr3*dixdj(3) + gf(2)*dixr(3)-gf(5)*rxqir(3)  &
       + gf(4)*(dixqjr(3)+djxqir(3)+rxqidj(3)-2.0d0*qixqj(3))&
       - gf(7)*(rxqijr(3)+qjrxqir(3))                         
  ttm3(1) = rr3*dixdj(1) + gf(3)*djxr(1) -gf(6)*rxqjr(1)  &
       - gf(4)*(dixqjr(1)+djxqir(1)+rxqjdi(1)-2.0d0*qixqj(1))&
       - gf(7)*(rxqjir(1)-qjrxqir(1))                         
  ttm3(2) = rr3*dixdj(2) + gf(3)*djxr(2) -gf(6)*rxqjr(2)  &
       - gf(4)*(dixqjr(2)+djxqir(2)+rxqjdi(2)-2.0d0*qixqj(2))&
       - gf(7)*(rxqjir(2)-qjrxqir(2))                         
  ttm3(3) = rr3*dixdj(3) + gf(3)*djxr(3) -gf(6)*rxqjr(3)  &
       - gf(4)*(dixqjr(3)+djxqir(3)+rxqjdi(3)-2.0d0*qixqj(3))&
       - gf(7)*(rxqjir(3)-qjrxqir(3))

  !     update force and torque on site j

  ftm1(1,j) = ftm1(1,j) + ftm2(1)
  ftm1(2,j) = ftm1(2,j) + ftm2(2)
  ftm1(3,j) = ftm1(3,j) + ftm2(3)
  ttm1(1,j) = ttm1(1,j) + ttm3(1)
  ttm1(2,j) = ttm1(2,j) + ttm3(2)
  ttm1(3,j) = ttm1(3,j) + ttm3(3)

  !     update force and torque on site i

  ftm1(1,i) = ftm1(1,i) - ftm2(1)
  ftm1(2,i) = ftm1(2,i) - ftm2(2)
  ftm1(3,i) = ftm1(3,i) - ftm2(3)
  ttm1(1,i) = ttm1(1,i) + ttm2(1)
  ttm1(2,i) = ttm1(2,i) + ttm2(2)
  ttm1(3,i) = ttm1(3,i) + ttm2(3)
