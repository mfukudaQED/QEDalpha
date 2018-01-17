!=============================================================================
!
!
!  2010.11.21
! -subroutines to calculate spinor bilinear integrations.
! - N
! - T, M, V, h
!
! 12.2
! -two electron integrals (calc_inttwoele_nmpq & inttwoele_mat)
! 
! 2011.1.4
! - E(R) (calc_intE_pq & calc_intE_mat)
!
! 1.6
! - integration subroutines for photon terms
!   Fourier transformation of j (calc_intF_pq & intF_mat)
!
! 1.7
! - calligraphic F (continuous momentum version, intcalF_mat)
! - utility routines for discretization (convert_label, convert_label_2)
!
! 1.8
! - calligraphic F (discretized momentum version, intcalFj_mat)
!
! 1.10
! -intN_ele -> added NEL=odd case
!
! 2.2
! - defined function intcalFj_dp_mat
!
! 3.24
! - modified calc_intE_mat (include a facotor *(Ze/4.d0/PI))
!
! 5.25
! - fixed bug in intcalFj_dp_mat (a factor in tilde Delta p was missing)
!
! 5.26
! - modified calc_intE_mat (changed definition of E)
!
! 6.3
! - definded function to calculate tilde(Delta p)_j (function DeltaPj(j))
! - subroutine calc_DeltaPj_and_P0j(j,DeltaPj,Pj)
!
! 11.4
! - function intTM_mat (T+M)
!
! 12.5 
! - removed function DeltaPj
!
! 2012.5.19
! - efficient routines for two-electron integrals!
!    calc_inttwoele_mat(inttwoele)
!      calc_inttwoele_pg_mat(NL,NS,pg,inttwoele_pg_mat)
!        inttwoele_pg(i,j,k,l,s_i,s_j,s_k,s_l,NL,NS,pg)
!
! 11.1
! Routines to calculate electric field integrals.  
!    calc_intE_mat(kk,posR,intE)
!      calc_intE_pg_mat(kk,posR,NL,NS,pg,intE_pg_mat)
!        intE_pg(kk,posR,i,j,s_i,s_j,NL,NS,pg)
!      
! 11.6
!  subroutine calc_mode_pj(j,P0j,vecPj,vec_pol)  (for calculation of dArad/dt) )
!
! 11.15
! routines to compute integration for constant outer electric current.
!    calc_intJcA_mat(posR,intJcA)
!      subroutine calc_intV2_pg_mat(posR,NL,NS,pg,intV2_pg_mat)
!        intV_pg(posR,i,j,s_i,s_j,NL,NS,pg) --> this is also used for nuclear integral
!
! To do
! -Taking spinor sum afterward is inefficient.
!  (it is possible to reduce calls of gauss_int_twoele etc.)
! -Calculation of Kramers pair integration (Use negative value of M.O. to represent KP)
!
! -impose symmetry ? h_nm = (h_mn)^* etc
! -Kramers symmetry
!
!============================================================
function intN_ele(NEL)
!============================================================
  use Precision
  implicit none
 
  real(kind=dp) :: intN_ele
  integer,intent(in) :: NEL
  
  complex(kind=dp) :: intN_mat
  complex(kind=dp) :: sum
  integer :: i
  
  sum = (0._dp,0._dp)
  if(mod(NEL,2).eq.0) then
     do i=1,NEL/2
        sum = sum +intN_mat(i,"+",i,"+")*2._dp
        !     write(*,*) i,intN_mat(i,"+",i,"+")
     end do
  else
     do i=1,(NEL-1)/2
        sum = sum +intN_mat(i,"+",i,"+")*2._dp
     end do
     i=(NEL+1)/2
     sum = sum +intN_mat(i,"+",i,"+")
  end if
     
  intN_ele = real(sum)
  return
end function intN_ele

!-----------------------------------------------------------------------------
! Routines for integration used in photon terms (calligraphic F)
!-----------------------------------------------------------------------------

!======================================================================
subroutine calc_DeltaPj_and_P0j(j,DeltaPj,P0j)
! calculate tilde(Delta p_j)   (the volume element for radiation momentum)
! and p0_j for given j.
! j denote discretized photon momentum and polarization collectively.
!
! use p0,th,phi mesh consistent with intcalFj_mat
!======================================================================
  use Precision
  use Constants  ! use P0MAX,NP0,NPTH,NPPHI
  implicit none
  
  integer,intent(in) :: j ! collective index which expresees photon momentum and polarization
  real(kind=dp),intent(out) :: DeltaPj,P0j

  real(kind=dp) :: p0,th,phi ! photon momentum vector in spherical coord.
  character(LEN=1) :: sig ! polarization (+ or -)  
  integer :: ip0,ith,iphi
  real(kind=dp) :: dp0,dth,dphi ! mesh width
  
  dp0  = P0MAX/NP0 ! P0MAX specified in GammaMatrix. We do not use p0=0
  dth  = PI/(NPTH-1) ! NPTH specified in GammaMatrix. th include both 0 and PI
  dphi = 2._dp*PI/NPPHI ! NPPHI specified in GammaMatrix. phi include 0 but not 2PI

  call convert_label_2(j,NP0,NPTH,NPPHI,ip0,ith,iphi,sig)
  ! (DeltaPj does not depend on sig and phi)

  p0 = dp0*ip0
  th = dth*(ith-1)  ! start from 0
  phi= dphi*(iphi-1) ! start from 0

  DeltaPj = sqrt(1._dp/(2._dp*PI**2*CCC)) *dp0*dth*dphi*p0**2*sin(th)/sqrt(2._dp*p0)  
  P0j = p0

  return
end subroutine calc_DeltaPj_and_P0j

!======================================================================
subroutine calc_mode_pj(j,P0j,vecPj,vecpolj)
! Calculate vec{p}, p0_j and e^k(polarization vector) for given j.
! (convert from spherical to cartesian)
! j denote discretized photon momentum and polarization collectively.
!
! 121106 
!======================================================================
  use Precision
  use Constants  ! use P0MAX,NP0,NPTH,NPPHI
  implicit none
  
  integer,intent(in) :: j ! collective index which expresees photon momentum and polarization
  real(kind=dp),intent(out) :: P0j,vecPj(3)
  complex(kind=dp),intent(out) :: vecpolj(3) ! (circular) polarization vector (depends on sig)
  
  real(kind=dp) :: p0,th,phi ! photon momentum vector in spherical coord.
  character(LEN=1) :: sig ! polarization (+ or -)  
  integer :: ip0,ith,iphi
  real(kind=dp) :: dp0,dth,dphi ! mesh width
  integer :: k
  real(kind=dp) :: sin_th,cos_th,sin_phi,cos_phi

  dp0  = P0MAX/NP0 ! P0MAX specified in GammaMatrix. We do not use p0=0
  dth  = PI/(NPTH-1) ! NPTH specified in GammaMatrix. th include both 0 and PI
  dphi = 2._dp*PI/NPPHI ! NPPHI specified in GammaMatrix. phi include 0 but not 2PI
  
  call convert_label_2(j,NP0,NPTH,NPPHI,ip0,ith,iphi,sig)
  
  p0 = dp0*ip0
  th = dth*(ith-1)  ! start from 0
  phi= dphi*(iphi-1) ! start from 0
  
  P0j = p0
  
  ! avoid numerical error due to the truncation of PI
  ! 0 <= th <= PI
  sin_th = sin(th)

  if(abs(th-PI/2._dp).lt.5.e-15_dp) then   
     cos_th = 0._dp
  else
     cos_th = cos(th)
  end if

  ! 0 <= phi < 2*PI
  if(abs(phi-PI).lt.5.e-15_dp) then
     sin_phi = 0._dp
  else
     sin_phi = sin(phi)
  end if

  if(abs(phi-PI/2._dp).lt.5.e-15_dp) then   
     cos_phi = 0._dp
  elseif(abs(phi-1.5_dp*PI).lt.5.e-15_dp) then
     cos_phi = 0._dp
  else
     cos_phi = cos(phi)
  end if
  
  vecPj(1) = p0*sin_th*cos_phi
  vecPj(2) = p0*sin_th*sin_phi
  vecPj(3) = p0*cos_th


!!$!  if(abs(th-1.57079632679490_dp).lt.5.e-15_dp) then
!!$  if(abs(th-PI/2._dp).lt.5.e-15_dp) then
!!$     vecPj(1) = p0*cos(phi)
!!$     vecPj(2) = p0*sin(phi)
!!$     vecPj(3) = 0._dp
!!$  else
!!$     vecPj(1) = p0*sin(th)*cos(phi)
!!$     vecPj(2) = p0*sin(th)*sin(phi)
!!$     vecPj(3) = p0*cos(th)
!!$  end if


  
!  write(*,*) p0,th,phi
!  write(*,*) vecPj

  if(sig.eq."+") then
!!$     vecpolj(1) = cos(phi)*cos(th) -IU*sin(phi)
!!$     vecpolj(2) = sin(phi)*cos(th) +IU*cos(phi)
!!$     vecpolj(3) = -sin(th)
     vecpolj(1) = cos_phi*cos_th -IU*sin_phi
     vecpolj(2) = sin_phi*cos_th +IU*cos_phi
     vecpolj(3) = -sin_th
  elseif(sig.eq."-") then
!!$     vecpolj(1) = cos(phi)*cos(th) +IU*sin(phi)
!!$     vecpolj(2) = sin(phi)*cos(th) -IU*cos(phi)
!!$     vecpolj(3) = -sin(th)
     vecpolj(1) = cos_phi*cos_th +IU*sin_phi
     vecpolj(2) = sin_phi*cos_th -IU*cos_phi
     vecpolj(3) = -sin_th
  else
     write(*,*) "sig should be + or -."
     stop
  end if
  do k=1,3
     vecpolj(k) = vecpolj(k)/sqrt(2._dp)
  end do

  return
end subroutine calc_mode_pj


!======================================================================
subroutine get_pj(j,P0j,thetaj,phij,dv)
! Calculate p0_j, theta_j, phi_j and dp0*dtheta*dphi
! j denote discretized photon momentum and polarization collectively.
!
! 130124
!======================================================================
  use Precision
  use Constants  ! use P0MAX,NP0,NPTH,NPPHI
  implicit none
  
  integer,intent(in) :: j ! collective index which expresees photon momentum and polarization
  real(kind=dp),intent(out) :: P0j,thetaj,phij,dv
  
  real(kind=dp) :: p0,th,phi ! photon momentum vector in spherical coord.
  character(LEN=1) :: sig ! polarization (+ or -)  
  integer :: ip0,ith,iphi
  real(kind=dp) :: dp0,dth,dphi ! mesh width

  dp0  = P0MAX/NP0 ! P0MAX specified in GammaMatrix. We do not use p0=0
  dth  = PI/(NPTH-1) ! NPTH specified in GammaMatrix. th include both 0 and PI
  dphi = 2._dp*PI/NPPHI ! NPPHI specified in GammaMatrix. phi include 0 but not 2PI

  dv = dp0*dth*dphi

  call convert_label_2(j,NP0,NPTH,NPPHI,ip0,ith,iphi,sig)
  
  p0 = dp0*ip0
  th = dth*(ith-1)  ! start from 0
  phi= dphi*(iphi-1) ! start from 0
  
  P0j    = p0
  thetaj = th
  phij   = phi
  
  return
end subroutine get_pj


!!$!======================================================================
!!$function DeltaPj(j)
!!$! tilde(Delta p_j)   (the volume element for radiation momentum)
!!$! j denote discretized photon momentum and polarization collectively.
!!$!
!!$! use p0,th,phi mesh consistent with intcalFj_mat
!!$!======================================================================
!!$  use Constants  ! use P0MAX,NP0,NPTH,NPPHI
!!$  implicit none
!!$
!!$  real(kind=dp) :: DeltaPj
!!$  integer,intent(in) :: j ! collective index which expresees photon momentum and polarization
!!$
!!$  real(kind=dp) :: p0,th,phi ! photon momentum vector in spherical coord.
!!$  character(LEN=1) :: sig ! polarization (+ or -)  
!!$  integer :: ip0,ith,iphi
!!$  real(kind=dp) :: dp0,dth,dphi ! mesh width
!!$
!!$  dp0  = P0MAX/NP0 ! P0MAX specified in GammaMatrix. We do not use p0=0
!!$  dth  = PI/(NPTH-1) ! NPTH specified in GammaMatrix. th include both 0 and PI
!!$  dphi = 2.d0*PI/NPPHI ! NPPHI specified in GammaMatrix. phi include 0 but not 2PI
!!$
!!$  call convert_label_2(j,NP0,NPTH,NPPHI,ip0,ith,iphi,sig)
!!$  ! (DeltaPj does not depend on sig and phi)
!!$
!!$  p0 = dp0*ip0
!!$  th = dth*(ith-1)  ! start from 0
!!$  phi= dphi*(iphi-1) ! start from 0
!!$
!!$  DeltaPj = sqrt(1.d0/(2.d0*PI**2*CCC)) *dp0*dth*dphi*p0**2*sin(th)/sqrt(2.d0*p0)  
!!$
!!$  return
!!$end function DeltaPj


!============================================================
function intcalFj_dp_mat(t,j,p,a,q,b)
! [Cal F]^j_paqb(t) x tilde(Delta p_j)   (multiplied by the volume element)
! j denote discretized photon momentum and polarization collectively.
! p,q : labels for molecular orbitals.
! a,b : labels for electron/positron ("+"/"-")
!============================================================
  use Precision
  use Constants  ! use P0MAX,NP0,NPTH,NPPHI
  implicit none

  complex(kind=dp) :: intcalFj_dp_mat
  real(kind=dp),intent(in) :: t  ! time
  integer,intent(in) :: j ! collective index which expresees photon momentum and polarization
  integer,intent(in) :: p,q
  character(LEN=1),intent(in) :: a,b

  complex(kind=dp) :: intcalF_mat
  real(kind=dp) :: p0,th,phi ! photon momentum vector in spherical coord.
  character(LEN=1) :: sig ! polarization (+ or -)  
  integer :: ip0,ith,iphi
  real(kind=dp) :: dp0,dth,dphi ! mesh width
  real(kind=dp) :: Delta_p

  dp0  = P0MAX/NP0 ! P0MAX specified in GammaMatrix. We do not use p0=0
  dth  = PI/(NPTH-1) ! NPTH specified in GammaMatrix. th include both 0 and PI
  dphi = 2._dp*PI/NPPHI ! NPPHI specified in GammaMatrix. phi include 0 but not 2PI

  call convert_label_2(j,NP0,NPTH,NPPHI,ip0,ith,iphi,sig)

  p0 = dp0*ip0
  th = dth*(ith-1)  ! start from 0
  phi= dphi*(iphi-1) ! start from 0

!  write(*,*) ip0,ith,iphi
!  write(*,*) p0,th/PI,phi/PI,sig

  Delta_p = sqrt(1._dp/(2._dp*PI**2*CCC)) *dp0*dth*dphi*p0**2*sin(th)/sqrt(2._dp*p0)  ! modified 110525

  intcalFj_dp_mat = intcalF_mat(t,p0,th,phi,sig,p,a,q,b)*Delta_p

  return
end function intcalFj_dp_mat

!============================================================
function intcalFj_mat(t,j,p,a,q,b)
! [Cal F]^j_paqb(t) 
! j denote discretized photon momentum and polarization collectively.
! p,q : labels for molecular orbitals.
! a,b : labels for electron/positron ("+"/"-")
!============================================================
  use Precision
  use Constants  ! use P0MAX,NP0,NPTH,NPPHI
  implicit none

  complex(kind=dp) :: intcalFj_mat
  real(kind=dp),intent(in) :: t  ! time
  integer,intent(in) :: j ! collective index which expresees photon momentum and polarization
  integer,intent(in) :: p,q
  character(LEN=1),intent(in) :: a,b

  complex(kind=dp) :: intcalF_mat
  real(kind=dp) :: p0,th,phi ! photon momentum vector in spherical coord.
  character(LEN=1) :: sig ! polarization (+ or -)  
  integer :: ip0,ith,iphi
  real(kind=dp) :: dp0,dth,dphi ! mesh width

  dp0  = P0MAX/NP0 ! P0MAX specified in GammaMatrix. We do not use p0=0
  dth  = PI/(NPTH-1) ! NPTH specified in GammaMatrix. th include both 0 and PI
  dphi = 2._dp*PI/NPPHI ! NPPHI specified in GammaMatrix. phi include 0 but not 2PI

  call convert_label_2(j,NP0,NPTH,NPPHI,ip0,ith,iphi,sig)

  p0 = dp0*ip0
  th = dth*(ith-1)  ! start from 0
  phi= dphi*(iphi-1) ! start from 0

!  write(*,*) ip0,ith,iphi
!  write(*,*) p0,th/PI,phi/PI,sig

  intcalFj_mat = intcalF_mat(t,p0,th,phi,sig,p,a,q,b)

  return
end function intcalFj_mat

!============================================================
function intcalF_mat(t,p0,th,phi,sig,p,a,q,b)
! [Cal F]_paqb(p,sigma,t) 
! p (photon momentum) specified in spherical coordinate (p0,th,phi)
! sig (+/-) specifies rotating direction of circular polarization (which depends on th and phi)
! p,q : labels for molecular orbitals.
! a,b : labels for electron/positron ("+"/"-")
!============================================================
  use Precision
  use DefineTypes
  use DiracOutput
  use Constants
  implicit none

  complex(kind=dp) :: intcalF_mat
  real(kind=dp),intent(in) :: t  ! time
  real(kind=dp),intent(in) :: p0,th,phi ! photon momentum vector in spherical coord.
  character(LEN=1),intent(in) :: sig ! polarization (+ or -)  
  integer,intent(in) :: p,q
  character(LEN=1),intent(in) :: a,b

  real(kind=dp) :: vecP(3) ! photon momentum vector in Cartesian coord.
  complex(kind=dp) :: vec_pol(3) ! (circular) polarization vector (depends on sig)
  integer :: k
  complex(kind=dp) :: intF_mat

  vecP(1) = p0*sin(th)*cos(phi)
  vecP(2) = p0*sin(th)*sin(phi)
  vecP(3) = p0*cos(th)

  if(sig.eq."+") then
     vec_pol(1) = cos(phi)*cos(th) -IU*sin(phi)
     vec_pol(2) = sin(phi)*cos(th) +IU*cos(phi)
     vec_pol(3) = -sin(th)
  elseif(sig.eq."-") then
     vec_pol(1) = cos(phi)*cos(th) +IU*sin(phi)
     vec_pol(2) = sin(phi)*cos(th) -IU*cos(phi)
     vec_pol(3) = -sin(th)
  else
     write(*,*) "sig should be + or - in intcalF_mat."
     stop
  end if
  do k=1,3
     vec_pol(k) = vec_pol(k)/sqrt(2._dp)
  end do

  intcalF_mat = (0._dp,0._dp)
  do k=1,3
     intcalF_mat = intcalF_mat +intF_mat(k,vecP,p,a,q,b)*vec_pol(k)
  end do
  intcalF_mat = intcalF_mat*exp(-IU*CCC*p0*t)  ! hbar=1

  return
end function intcalF_mat

!-----------------------------------------------------------------------------
! Routines for spinor version two-electron integration 
!-----------------------------------------------------------------------------
!============================================================
!-----------------------------------------------------------------------------
! This subroutine will not be used because it is too slow. 
!-----------------------------------------------------------------------------
function inttwoele_mat(n,a,m,b,p,c,q,d)
! Two-electron integration
! n,m,p,q : labels for molecular orbitals.
! a,b,c,d : labels for electron/positron ("+"/"-")
!============================================================
  use Precision
  use DefineTypes
  use DiracOutput
  use Constants
  implicit none
 
  complex(kind=dp) :: inttwoele_mat
  integer,intent(in) :: n,m,p,q
  character(LEN=1),intent(in) :: a,b,c,d

  type(primitive_gaussian) :: pg
  complex(kind=dp) :: c_n(4,NMAX_PG),c_m(4,NMAX_PG)
  complex(kind=dp) :: c_p(4,NMAX_PG),c_q(4,NMAX_PG)
  complex(kind=dp) :: int_nmpq
  integer :: i,j

  call copy_DiracOutput(n,a,m,b,pg,c_n,c_m)  ! set pg,c_n,c_m
  call copy_DiracOutput(p,c,q,d,pg,c_p,c_q)  ! set pg,c_p,c_q
  
  inttwoele_mat = (0._dp,0._dp)
  do i=1,4
     do j=1,4
        call calc_inttwoele_nmpq(i,i,j,j,NBS_L,NBS_S,pg,c_n,c_m,c_p,c_q,int_nmpq) 
        inttwoele_mat = inttwoele_mat +int_nmpq
     end do
  end do
  return
end function inttwoele_mat


!-----------------------------------------------------------------------------
! Routines for integration of spinor bilinears.
!-----------------------------------------------------------------------------
!============================================================
function intF_mat(k,vecP,p,a,q,b)
! F^k_paqb(p) : Fourier transform of current j (k-th component)
! p,q : labels for molecular orbitals.
! a,b : labels for electron/positron ("+"/"-")
!============================================================
  use Precision
  use DefineTypes
  use DiracOutput
  use Constants
  implicit none
 
  complex(kind=dp) :: intF_mat
  integer,intent(in) :: k 
  real(kind=dp),intent(in) :: vecP(3)
  integer,intent(in) :: p,q
  character(LEN=1),intent(in) :: a,b

  type(primitive_gaussian) :: pg
  complex(kind=dp) :: c_p(4,NMAX_PG),c_q(4,NMAX_PG)
  ! same gamma matrix (gamma^0 gamma^k) calculation as in intT_mat
  complex(kind=dp) :: intF_pq_14,intF_pq_23,intF_pq_32,intF_pq_41
  complex(kind=dp) :: intF_pq_13,intF_pq_24,intF_pq_31,intF_pq_42
  
  call copy_DiracOutput(p,a,q,b,pg,c_p,c_q)  ! set pg,c_p,c_q
  
  if(k.eq.1 .or. k.eq.2) then 
     call calc_intF_pq(1,4,vecP,NBS_L,NBS_S,pg,c_p,c_q,intF_pq_14) 
     call calc_intF_pq(2,3,vecP,NBS_L,NBS_S,pg,c_p,c_q,intF_pq_23) 
     call calc_intF_pq(3,2,vecP,NBS_L,NBS_S,pg,c_p,c_q,intF_pq_32) 
     call calc_intF_pq(4,1,vecP,NBS_L,NBS_S,pg,c_p,c_q,intF_pq_41) 
     if(k.eq.1) then
        intF_mat = intF_pq_14+intF_pq_23+intF_pq_32+intF_pq_41
     elseif(k.eq.2) then
        intF_mat = -IU*(intF_pq_14+intF_pq_32) +IU*(intF_pq_23+intF_pq_41)
     end if
  elseif(k.eq.3) then
     call calc_intF_pq(1,3,vecP,NBS_L,NBS_S,pg,c_p,c_q,intF_pq_13) 
     call calc_intF_pq(2,4,vecP,NBS_L,NBS_S,pg,c_p,c_q,intF_pq_24) 
     call calc_intF_pq(3,1,vecP,NBS_L,NBS_S,pg,c_p,c_q,intF_pq_31) 
     call calc_intF_pq(4,2,vecP,NBS_L,NBS_S,pg,c_p,c_q,intF_pq_42) 
     intF_mat = (intF_pq_13+intF_pq_31) -(intF_pq_24+intF_pq_42)
  else
     write(*,*) "k should be 1-3 in intF_mat"
     stop
  end if

  intF_mat = intF_mat*(Ze*CCC)
  return
end function intF_mat

!============================================================
subroutine calc_intE_mat(posR,p,a,q,b,intE_mat)
! Integration of propto (-grad V) at posR
! p,q : labels for molecular orbitals.
! a,b : labels for electron/positron ("+"/"-")
!
! used in function integrand_L in sub_intret.f90
! 110324 : modified
! 110526 : modified
!============================================================
  use Precision
  use DefineTypes
  use DiracOutput
  use Constants
  implicit none
 
  complex(kind=dp),intent(out) :: intE_mat(3)
  real(kind=dp),intent(in) :: posR(3)
  integer,intent(in) :: p,q
  character(LEN=1),intent(in) :: a,b

  type(primitive_gaussian) :: pg
  complex(kind=dp) :: c_p(4,NMAX_PG),c_q(4,NMAX_PG)
  complex(kind=dp) :: intE_pq(3)
  integer :: i,k

  call copy_DiracOutput(p,a,q,b,pg,c_p,c_q)  ! set pg,c_p,c_q
  
  intE_mat(1:3) = (0._dp,0._dp)
  
  do i=1,4
     call calc_intE_pq(i,i,posR,NBS_L,NBS_S,pg,c_p,c_q,intE_pq) 
     do k=1,3
        intE_mat(k) = intE_mat(k) +intE_pq(k)
     end do
  end do

  do k=1,3
!     intE_mat(k) = intE_mat(k) *(Ze/4._dp/PI)  ! 110324 redefined
     intE_mat(k) = intE_mat(k) *(-Ze/4._dp/PI)  ! 110526 redefined
  end do

  return
end subroutine calc_intE_mat

!============================================================
function intTM_mat(p,a,q,b)
! Integration of kinetic and mass terms T+M
! used for non-BO calculation
! p,q : labels for molecular orbitals.
! a,b : labels for electron/positron ("+"/"-")
!============================================================
  use Precision
  implicit none

  complex(kind=dp) :: intTM_mat
  integer,intent(in) :: p,q
  character(LEN=1),intent(in) :: a,b
  
  complex(kind=dp) :: intT_mat,intM_mat

  intTM_mat = intT_mat(p,a,q,b)+intM_mat(p,a,q,b)

  return
end function intTM_mat


!============================================================
function inth_mat(p,a,q,b)
! Integration of one-electron terms
! T+M+V
! p,q : labels for molecular orbitals.
! a,b : labels for electron/positron ("+"/"-")
!============================================================
  use Precision
  implicit none

  complex(kind=dp) :: inth_mat
  integer,intent(in) :: p,q
  character(LEN=1),intent(in) :: a,b

  complex(kind=dp) :: intT_mat,intM_mat,intV_mat

  inth_mat = intT_mat(p,a,q,b)+intM_mat(p,a,q,b)+intV_mat(p,a,q,b)

!  inth_mat = intT_mat(p,a,q,b)+intV_mat(p,a,q,b)
!  inth_mat = intM_mat(p,a,q,b)+intV_mat(p,a,q,b)
!  inth_mat = intT_mat(p,a,q,b)+intM_mat(p,a,q,b)
!  inth_mat = intM_mat(p,a,q,b)
!  inth_mat = intT_mat(p,a,q,b)
!  inth_mat = intV_mat(p,a,q,b)


  return
end function inth_mat


!============================================================
function intV_mat(p,a,q,b)
! Integration of nuclear attraction energy.
! Summed over nulcei.
! p,q : labels for molecular orbitals.
! a,b : labels for electron/positron ("+"/"-")
!============================================================
  use Precision
  use DefineTypes
  use DiracOutput
  use Constants
  implicit none
 
  complex(kind=dp) :: intV_mat
  integer,intent(in) :: p,q
  character(LEN=1),intent(in) :: a,b

  type(primitive_gaussian) :: pg
  complex(kind=dp) :: c_p(4,NMAX_PG),c_q(4,NMAX_PG)
  complex(kind=dp) :: intV_pq
  real(kind=dp) :: posRA(3)
  integer :: i,iN
  
  call copy_DiracOutput(p,a,q,b,pg,c_p,c_q)  ! set pg,c_p,c_q
  
  intV_mat = (0._dp,0._dp)
  do iN=1,NAT
     posRA(1) = xc(iN); posRA(2) = yc(iN); posRA(3) = zc(iN)
     do i=1,4
        call calc_intV_pq(i,i,posRA,NBS_L,NBS_S,pg,c_p,c_q,intV_pq) 
        intV_mat = intV_mat +intV_pq*cn(iN)  ! multipy by Za
     end do
  end do
  intV_mat = intV_mat*Ze
  return
end function intV_mat

!============================================================
function intM_mat(p,a,q,b)
! Integration of mass energy: me c^2 beta
! p,q : labels for molecular orbitals.
! a,b : labels for electron/positron ("+"/"-")
!============================================================
  use Precision
  use DefineTypes
  use DiracOutput
  use Constants
  implicit none
 
  complex(kind=dp) :: intM_mat
  integer,intent(in) :: p,q
  character(LEN=1),intent(in) :: a,b
  
  type(primitive_gaussian) :: pg
  complex(kind=dp) :: c_p(4,NMAX_PG),c_q(4,NMAX_PG)
  complex(kind=dp) :: intM_pq_11,intM_pq_22,intM_pq_33,intM_pq_44
  
  call copy_DiracOutput(p,a,q,b,pg,c_p,c_q)  ! set pg,c_p,c_q

  ! beta = gamma^0 = diag(1,1,-1,-1)

  call calc_intN_pq(1,1,NBS_L,NBS_S,pg,c_p,c_q,intM_pq_11) 
  call calc_intN_pq(2,2,NBS_L,NBS_S,pg,c_p,c_q,intM_pq_22) 
  call calc_intN_pq(3,3,NBS_L,NBS_S,pg,c_p,c_q,intM_pq_33) 
  call calc_intN_pq(4,4,NBS_L,NBS_S,pg,c_p,c_q,intM_pq_44) 

  intM_mat = (intM_pq_11 +intM_pq_22 -intM_pq_33 -intM_pq_44)*CCC**2

  ! if me c^2 is subtracted ( me c^2 (beta -I) )
  ! as is done in relativistic quantum chemistry calculation
!  intM_mat = 2._dp*(-intM_pq_33 -intM_pq_44)*CCC**2

  return
end function intM_mat

!============================================================
function intT_mat(p,a,q,b)
! Integration of kinetic energy: c.alpha.p = c.alpha.(-i hbar d/dx)
! p,q : labels for molecular orbitals.
! a,b : labels for electron/positron ("+"/"-")
!============================================================
  use Precision
  use DefineTypes
  use DiracOutput
  use Constants
  implicit none
 
  complex(kind=dp) :: intT_mat
  integer,intent(in) :: p,q
  character(LEN=1),intent(in) :: a,b
  
  type(primitive_gaussian) :: pg
  complex(kind=dp) :: c_p(4,NMAX_PG),c_q(4,NMAX_PG)
  complex(kind=dp) :: intT_pq_14(3),intT_pq_23(3),intT_pq_32(3),intT_pq_41(3)
  complex(kind=dp) :: intT_pq_13(3),intT_pq_24(3),intT_pq_31(3),intT_pq_42(3)
  complex(kind=dp) :: intT_mat1,intT_mat2,intT_mat3
  
  call copy_DiracOutput(p,a,q,b,pg,c_p,c_q)  ! set pg,c_p,c_q

  call calc_intT_pq(1,4,NBS_L,NBS_S,pg,c_p,c_q,intT_pq_14) 
  call calc_intT_pq(2,3,NBS_L,NBS_S,pg,c_p,c_q,intT_pq_23) 
  call calc_intT_pq(3,2,NBS_L,NBS_S,pg,c_p,c_q,intT_pq_32) 
  call calc_intT_pq(4,1,NBS_L,NBS_S,pg,c_p,c_q,intT_pq_41) 
  call calc_intT_pq(1,3,NBS_L,NBS_S,pg,c_p,c_q,intT_pq_13) 
  call calc_intT_pq(2,4,NBS_L,NBS_S,pg,c_p,c_q,intT_pq_24) 
  call calc_intT_pq(3,1,NBS_L,NBS_S,pg,c_p,c_q,intT_pq_31) 
  call calc_intT_pq(4,2,NBS_L,NBS_S,pg,c_p,c_q,intT_pq_42) 

  intT_mat1 = intT_pq_14(1)+intT_pq_23(1)+intT_pq_32(1)+intT_pq_41(1)
  intT_mat2 = -IU*(intT_pq_14(2)+intT_pq_32(2)) +IU*(intT_pq_23(2)+intT_pq_41(2))
  intT_mat3 = (intT_pq_13(3)+intT_pq_31(3)) -(intT_pq_24(3)+intT_pq_42(3))

  intT_mat = (intT_mat1+intT_mat2+intT_mat3)*(-IU*CCC)

  return
end function intT_mat

!============================================================
function intN_mat(p,a,q,b)
! p,q : labels for molecular orbitals.
! a,b : labels for electron/positron ("+"/"-")
!============================================================
  use Precision
  use DefineTypes
  use DiracOutput
  implicit none
 
  complex(kind=dp) :: intN_mat
  integer,intent(in) :: p,q
  character(LEN=1),intent(in) :: a,b

  type(primitive_gaussian) :: pg
  complex(kind=dp) :: c_p(4,NMAX_PG),c_q(4,NMAX_PG)
  complex(kind=dp) :: intN_pq
  integer :: i
  
  call copy_DiracOutput(p,a,q,b,pg,c_p,c_q)  ! set pg,c_p,c_q
  
  intN_mat = (0._dp,0._dp)
  do i=1,4
     call calc_intN_pq(i,i,NBS_L,NBS_S,pg,c_p,c_q,intN_pq) 
!     write(*,*) i,intN_pq
     intN_mat = intN_mat +intN_pq
  end do
  return
end function intN_mat

!====================================================================================================
!====================================================================================================
!====================================================================================================
!====================================================================================================
!====================================================================================================
!====================================================================================================

!-----------------------------------------------------------------------------
! Routines to calculate J(constant).A integrals.  121112
! (I^k_{JcA na mb}(vecR))
!
!    calc_intJcA_mat(posR,intJcA)
!      subroutine calc_intV2_pg_mat(posR,NL,NS,pg,intV2_pg_mat)
!        intV_pg(posR,i,j,s_i,s_j,NL,NS,pg) --> this is also used for nuclear integral
!
!-----------------------------------------------------------------------------

!=========================================================================
subroutine calc_intJcA_mat(posR,vecJ,intJcA)
! Ijca_{na,mb}
! = J^k int ds (psi(s)^+)_na gamma^0 gamma^k (psi(s))_mb /|s-R|   (sum over k=1~3)
! = Sum_{alpha,beta} J^k [gamma^0 gamma^k]_{alpha,beta}
!        int ds (psi(s)^+)_alpha_na (psi(s))_beta_mb /|s-R|  (sum over k=1~3)
! 2012.11.15
!   11.16 modified to include vecJ
!=========================================================================  
  use Precision
  use DefineTypes
  use DiracOutput
  use Constants
  implicit none
  
  real(kind=dp),intent(in) :: posR(3),vecJ(3)
  complex(kind=dp),intent(out) :: intJcA(-NBS:NBS,2,-NBS:NBS,2) ! IJcA_{na,mb}
  ! negative values account for Kramers partners. 0 is not used.

  type(primitive_gaussian) :: pg
  complex(kind=dp) :: c_pg(4,NMAX_PG)

  complex(kind=dp) :: intV2_pg_mat(NBS_S,NBS_S,4,4) ! V(ij)^{alpha,beta}

  integer :: kk ! vector index
  integer :: i,j ! p.g. indice
  integer :: n,m,a,b ! MO indice 
  integer :: alpha,beta ! spinor index

  complex(kind=dp) :: tmp1(NBS_S,-NBS:NBS,2,4,4) ! V(i m^b)^{alpha,beta}
  complex(kind=dp) :: tmp2(-NBS:NBS,2,-NBS:NBS,2,4,4) ! V(n^a m^b)^{alpha,beta}
  complex(kind=dp) :: sum
  complex(kind=dp) :: JkGam0k(4,4) ! sum_{k=1}^3 (J^k gamma^0 gamma^k )

  ! sum_{k=1}^3 (J^k gamma^0 gamma^k )
  do alpha=1,4
     do beta=1,4
        sum = (0._dp,0._dp)
        do kk=1,3
           sum = sum +vecJ(kk)*Gam0k(kk,alpha,beta)
        end do
        JkGam0k(alpha,beta) = sum
     end do
  end do

  call copy_DiracOutput_pg(pg) ! get pg

  ! V(i j)^{alpha,beta}
  call calc_intV2_pg_mat(posR,NBS_L,NBS_S,pg,intV2_pg_mat) 

  ! V(i m^b)^{alpha,beta}
  do alpha=1,4
     do beta=1,4
        do i=1,NBS_S
           do m = -NBS,NBS  !m=0 is not used.
              do b = 1,2
                 sum = (0._dp,0._dp)
                 call copy_DiracOutput_cp(m,b,c_pg)
                 do j=1,NBS_S
                    sum = sum +c_pg(beta,j)*intV2_pg_mat(i,j,alpha,beta)
                    ! c^beta_{m^b j} x V(i j)^{alpha,beta}
                 end do
                 tmp1(i,m,b,alpha,beta) = sum
              end do
           end do
        end do
     end do
  end do

  ! V(n^a m^b)^{alpha,beta}
  do alpha=1,4
     do beta=1,4
        do m = -NBS,NBS  !m=0 is not used.
           do b = 1,2
              do n = -NBS,NBS  !n=0 is not used.
                 do a = 1,2
                    sum = (0._dp,0._dp)
                    call copy_DiracOutput_cp(n,a,c_pg)
                    do i=1,NBS_S
                       sum = sum +conjg(c_pg(alpha,i))*tmp1(i,m,b,alpha,beta)
                       ! (c^alpha_{n^a i})^* x V(i m^b)^{alpha,beta}    
                    end do
                    tmp2(n,a,m,b,alpha,beta) = sum
                 end do
              end do
           end do
        end do
     end do
  end do
  
  ! IJcA_{na,mb}
  do m = -NBS,NBS  !m=0 is not used.
     do b = 1,2
        do n = -NBS,NBS  !n=0 is not used.
           do a = 1,2
                 
              sum = (0._dp,0._dp)
              do alpha=1,4
                 do beta=1,4
                    sum = sum + JkGam0k(alpha,beta)*tmp2(n,a,m,b,alpha,beta)
                 end do
              end do
              intJcA(n,a,m,b) = sum
              
           end do
        end do
     end do
  end do
  
  return
end subroutine calc_intJcA_mat

!!$!=========================================================================
!!$subroutine calc_intJcA_mat(posR,intJcA)
!!$! IJcA^k_{na,mb}
!!$! = int ds (psi(s)^+)_na gamma^0 gamma^k (psi(s))_mb /|s-R|  (k=1~3)
!!$! = Sum_{alpha,beta} [gamma^0 gamma^k]_{alpha,beta}
!!$!        int ds (psi(s)^+)_alpha_na (psi(s))_beta_mb /|s-R|  (k=1~3)
!!$! 2012.11.15
!!$!=========================================================================  
!!$  use Precision
!!$  use DefineTypes
!!$  use DiracOutput
!!$  use Constants
!!$  implicit none
!!$  
!!$  real(kind=dp),intent(in) :: posR(3)  
!!$  complex(kind=dp),intent(out) :: intJcA(3,-NBS:NBS,2,-NBS:NBS,2) ! IJcA^k_{na,mb}
!!$  ! negative values account for Kramers partners. 0 is not used.
!!$
!!$  type(primitive_gaussian) :: pg
!!$  complex(kind=dp) :: c_pg(4,NMAX_PG)
!!$
!!$  complex(kind=dp) :: intV2_pg_mat(NBS_S,NBS_S,4,4) ! V(ij)^{alpha,beta}
!!$
!!$  integer :: kk ! vector index
!!$  integer :: i,j ! p.g. indice
!!$  integer :: n,m,a,b ! MO indice 
!!$  integer :: alpha,beta ! spinor index
!!$
!!$  complex(kind=dp) :: tmp1(NBS_S,-NBS:NBS,2,4,4) ! V(i m^b)^{alpha,beta}
!!$  complex(kind=dp) :: tmp2(-NBS:NBS,2,-NBS:NBS,2,4,4) ! V(n^a m^b)^{alpha,beta}
!!$  complex(kind=dp) :: sum
!!$
!!$  call copy_DiracOutput_pg(pg) ! get pg
!!$
!!$  ! V(i j)^{alpha,beta}
!!$  call calc_intV2_pg_mat(posR,NBS_L,NBS_S,pg,intV2_pg_mat) 
!!$
!!$  ! V(i m^b)^{alpha,beta}
!!$  do alpha=1,4
!!$     do beta=1,4
!!$        do i=1,NBS_S
!!$           do m = -NBS,NBS  !m=0 is not used.
!!$              do b = 1,2
!!$                 sum = (0._dp,0._dp)
!!$                 call copy_DiracOutput_cp(m,b,c_pg)
!!$                 do j=1,NBS_S
!!$                    sum = sum +c_pg(beta,j)*intV2_pg_mat(i,j,alpha,beta)
!!$                    ! c^beta_{m^b j} x V(i j)^{alpha,beta}
!!$                 end do
!!$                 tmp1(i,m,b,alpha,beta) = sum
!!$              end do
!!$           end do
!!$        end do
!!$     end do
!!$  end do
!!$
!!$  ! V(n^a m^b)^{alpha,beta}
!!$  do alpha=1,4
!!$     do beta=1,4
!!$        do m = -NBS,NBS  !m=0 is not used.
!!$           do b = 1,2
!!$              do n = -NBS,NBS  !n=0 is not used.
!!$                 do a = 1,2
!!$                    sum = (0._dp,0._dp)
!!$                    call copy_DiracOutput_cp(n,a,c_pg)
!!$                    do i=1,NBS_S
!!$                       sum = sum +conjg(c_pg(alpha,i))*tmp1(i,m,b,alpha,beta)
!!$                       ! (c^alpha_{n^a i})^* x V(i m^b)^{alpha,beta}    
!!$                    end do
!!$                    tmp2(n,a,m,b,alpha,beta) = sum
!!$                 end do
!!$              end do
!!$           end do
!!$        end do
!!$     end do
!!$  end do
!!$  
!!$  ! IJcA^k_{na,mb}
!!$  ! gamma^0 gamma^k is a global variable Gam0k(k,alpha,beta)
!!$  do kk=1,3
!!$     do m = -NBS,NBS  !m=0 is not used.
!!$        do b = 1,2
!!$           do n = -NBS,NBS  !n=0 is not used.
!!$              do a = 1,2
!!$                 
!!$                 sum = (0._dp,0._dp)
!!$                 do alpha=1,4
!!$                    do beta=1,4
!!$                       sum = sum + Gam0k(kk,alpha,beta)*tmp2(n,a,m,b,alpha,beta)
!!$                    end do
!!$                 end do
!!$                 intJcA(kk,n,a,m,b) = sum
!!$
!!$              end do
!!$           end do
!!$        end do
!!$     end do
!!$  end do
!!$
!!$  return
!!$end subroutine calc_intJcA_mat

!=========================================================================
subroutine calc_intV2_pg_mat(posR,NL,NS,pg,intV2_pg_mat)
! Int ds (g^a(r)^+)_i (g^b(r))_j 1/|s-R| = V(ij)^{ab}
! posR : vecR of (V_{na mb}(vecR))
! Note that g depends on large or small components
! g^a = g^L if a=1,2
!     = g^S if a=3,4
! 2012.11.15
!=========================================================================  
  use Precision
  use DefineTypes
  implicit none
  
  real(kind=dp),intent(in) :: posR(3)  
  integer,intent(in) :: NL,NS
  type(primitive_gaussian),intent(in) :: pg
  complex(kind=dp),intent(out) :: intV2_pg_mat(NS,NS,4,4) ! V(ij)^{ab}

  integer :: i,j ! p.g. indice
  integer :: a,b ! spinor indice
  integer :: NLS
  complex(kind=dp) :: intV_pg

  intV2_pg_mat = (0._dp,0._dp)

  do a=1,4
     do b=1,4
        do i=1,NLS(a,NL,NS)
           do j=1,NLS(b,NL,NS)
              intV2_pg_mat(i,j,a,b) = intV_pg(posR,i,j,a,b,NL,NS,pg)
           end do
        end do
     end do
  end do

  return
end subroutine calc_intV2_pg_mat

!=========================================================================
function intV_pg(posR,i,j,s_i,s_j,NL,NS,pg)
! Int ds (g^si(r)^+)_i (g^sj(r))_j 1/|s-R|
! posR : vecR of (V_{na mb}(vecR))
! Note that g depends on large or small components
! g^si = g^L if si=1,2
!      = g^S if si=3,4
! 2012.11.15
!=========================================================================  
  use Precision
  use DefineTypes
  implicit none

  complex(kind=dp) :: intV_pg

  real(kind=dp),intent(in) :: posR(3)  
  integer,intent(in) :: i,j ! p.g. indice
  integer,intent(in) :: s_i,s_j ! spinor indice
  integer,intent(in) :: NL,NS
  type(primitive_gaussian),intent(in) :: pg
  
  integer :: NLS
  real(kind=dp) :: overlap,ef(3) 
  real(kind=dp) :: nucatt
  real(kind=dp) :: posA(3), posB(3) ! position of center of gaussian 
  real(kind=dp) :: alphaA, alphaB ! exponents
  integer :: nx,nbarx,ny,nbary,nz,nbarz
  real(kind=dp) :: norm_pg  ! function
  real(kind=dp) :: f_norm_pg

  call set_pg(s_i,i,NL,NS,pg,posA,alphaA,nx,ny,nz)
  call set_pg(s_j,j,NL,NS,pg,posB,alphaB,nbarx,nbary,nbarz)
  call gauss_int(posA,alphaA,nx,ny,nz,posB,alphaB,nbarx,nbary,nbarz,posR,overlap,nucatt,ef)
  f_norm_pg = norm_pg(alphaA,nx,ny,nz)*norm_pg(alphaB,nbarx,nbary,nbarz)

  intV_pg = nucatt*f_norm_pg

  return
end function intV_pg


!-----------------------------------------------------------------------------
! Routines to calculate electric field integrals.  121031
! (E^k_{na mb}(vecR))
!
!    calc_intE_mat(kk,posR,intE)
!      calc_intE_pg_mat(kk,posR,NL,NS,pg,intE_pg_mat)
!        intE_pg(kk,posR,i,j,s_i,s_j,NL,NS,pg)
!
!-----------------------------------------------------------------------------
!=========================================================================
Subroutine calc_intEk_mat(posR,intE)
! E^k_{na,mb}
! = -(Ze e)/(4pi) Sum_alpha
!        int ds (psi(s)^+)_alpha_na (psi(s))_alpha_mb (s-R)^k/|s-R|  (k=1~3)
! 2012.11.1
! named as calc_intEk_mat (with k after E) to distinguish from previously
! made routine calc_intE_mat (slow routines and should be replaced.)
!=========================================================================  
  use Precision
  use DefineTypes
  use DiracOutput
  use Constants
  implicit none
  
  real(kind=dp),intent(in) :: posR(3)  
  complex(kind=dp),intent(out) :: intE(3,-NBS:NBS,2,-NBS:NBS,2) ! E^k_{na,mb}
  ! negative values account for Kramers partners. 0 is not used.

  type(primitive_gaussian) :: pg
  complex(kind=dp) :: c_pg(4,NMAX_PG)

  complex(kind=dp) :: intE_pg_mat(3,NBS_S,NBS_S,4) ! E^k(ij)^{alpha}

  integer :: kk ! vector index
  integer :: i,j ! p.g. indice
  integer :: n,m,a,b ! MO indice 
  integer :: alpha ! spinor index

  complex(kind=dp) :: tmp1(3,NBS_S,-NBS:NBS,2,4) ! E^k(i m^b)^{alpha}
  complex(kind=dp) :: tmp2(3,-NBS:NBS,2,-NBS:NBS,2,4) ! E^k(n^a m^b)^{alpha}
  complex(kind=dp) :: sum

  call copy_DiracOutput_pg(pg) ! get pg

  ! E^k(i j)^{alpha}
  call calc_intE_pg_mat(posR,NBS_L,NBS_S,pg,intE_pg_mat) 

  ! E^k(i m^b)^{alpha}
  do kk=1,3
     do alpha=1,4
        do i=1,NBS_S
           do m = -NBS,NBS  !m=0 is not used.
              do b = 1,2
                 sum = (0._dp,0._dp)
                 call copy_DiracOutput_cp(m,b,c_pg)
                 do j=1,NBS_S
                    sum = sum +c_pg(alpha,j)*intE_pg_mat(kk,i,j,alpha)
                    ! c^alpha_{m^b j} x E^k(i j)^{alpha}
                 end do
                 tmp1(kk,i,m,b,alpha) = sum
              end do
           end do
        end do
     end do
  end do

  ! E^k(n^a m^b)^{alpha}
  do kk=1,3
     do alpha=1,4
        do m = -NBS,NBS  !m=0 is not used.
           do b = 1,2
              do n = -NBS,NBS  !n=0 is not used.
                 do a = 1,2
                    sum = (0._dp,0._dp)
                    call copy_DiracOutput_cp(n,a,c_pg)
                    do i=1,NBS_S
                       sum = sum +conjg(c_pg(alpha,i))*tmp1(kk,i,m,b,alpha)
                       ! (c^alpha_{n^a i})^* x E^k(i m^b)^{alpha}    
                    end do
                    tmp2(kk,n,a,m,b,alpha) = sum
                 end do
              end do
           end do
        end do
     end do
  end do
  
  ! E^k_{na,mb}
  do kk=1,3
     do m = -NBS,NBS  !m=0 is not used.
        do b = 1,2
           do n = -NBS,NBS  !n=0 is not used.
              do a = 1,2
                 
                 sum = (0._dp,0._dp)
                 do alpha=1,4
                    sum = sum +tmp2(kk,n,a,m,b,alpha)
                 end do
                 intE(kk,n,a,m,b) = sum *(-Ze/4._dp/PI)

              end do
           end do
        end do
     end do
  end do

  return
end subroutine calc_intEk_mat

!=========================================================================
subroutine calc_intE_pg_mat(posR,NL,NS,pg,intE_pg_mat)
! Int ds (g^a(r)^+)_i (g^a(r))_j (s-R)^k/|s-R|  (k=1~3)
! posR : vecR of (E^k_{na mb}(vecR))
! Note that g depends on large or small components
! g^a = g^L if a=1,2
!     = g^S if a=3,4
! compute only necessary spinor combinations
! 2012.11.1
!=========================================================================  
  use Precision
  use DefineTypes
  implicit none
  
  real(kind=dp),intent(in) :: posR(3)  
  integer,intent(in) :: NL,NS
  type(primitive_gaussian),intent(in) :: pg
  complex(kind=dp),intent(out) :: intE_pg_mat(3,NS,NS,4) ! E^k(ij)^{a}

  integer :: kk ! vector index
  integer :: i,j ! p.g. indice
  integer :: a ! spinor indice
  integer :: NLS
  complex(kind=dp) :: intE_pg

  intE_pg_mat = (0._dp,0._dp)

  do kk=1,3
     do a=1,4
        do i=1,NLS(a,NL,NS)
           do j=1,NLS(a,NL,NS)
              intE_pg_mat(kk,i,j,a) = intE_pg(kk,posR,i,j,a,a,NL,NS,pg)
!              write(*,"(4i6,2es20.10)") kk,i,j,a,intE_pg_mat(kk,i,j,a)
           end do
        end do
     end do
  end do
!  stop
  return
end subroutine calc_intE_pg_mat

!=========================================================================
function intE_pg(kk,posR,i,j,s_i,s_j,NL,NS,pg)
! Int ds (g^si(r)^+)_i (g^sj(r))_j (s-R)^kk/|s-R|
! kk : 1~3 (vector components), k of (E^k_{na mb}(vecR))
! posR : vecR of (E^k_{na mb}(vecR))
! Note that g depends on large or small components
! g^si = g^L if si=1,2
!      = g^S if si=3,4
! 2012.11.1
!=========================================================================  
  use Precision
  use DefineTypes
  implicit none

  complex(kind=dp) :: intE_pg

  integer,intent(in) :: kk ! vector index
  real(kind=dp),intent(in) :: posR(3)  
  integer,intent(in) :: i,j ! p.g. indice
  integer,intent(in) :: s_i,s_j ! spinor indice
  integer,intent(in) :: NL,NS
  type(primitive_gaussian),intent(in) :: pg
  
  integer :: NLS
  real(kind=dp) :: overlap,ef(3) 
  real(kind=dp) :: nucatt
  real(kind=dp) :: posA(3), posB(3) ! position of center of gaussian 
  real(kind=dp) :: alphaA, alphaB ! exponents
  integer :: nx,nbarx,ny,nbary,nz,nbarz
  real(kind=dp) :: norm_pg  ! function
  real(kind=dp) :: f_norm_pg

  call set_pg(s_i,i,NL,NS,pg,posA,alphaA,nx,ny,nz)
  call set_pg(s_j,j,NL,NS,pg,posB,alphaB,nbarx,nbary,nbarz)
  call gauss_int(posA,alphaA,nx,ny,nz,posB,alphaB,nbarx,nbary,nbarz,posR,overlap,nucatt,ef)
  f_norm_pg = norm_pg(alphaA,nx,ny,nz)*norm_pg(alphaB,nbarx,nbary,nbarz)

  intE_pg = ef(kk)*f_norm_pg

  return
end function intE_pg



!-----------------------------------------------------------------------------
! Routines to calculate two electron integrals. (n^a m^b | p^c q^d)
!
!    calc_inttwoele_mat(inttwoele)
!      calc_inttwoele_pg_mat(NL,NS,pg,inttwoele_pg_mat)
!        inttwoele_pg(i,j,k,l,s_i,s_j,s_k,s_l,NL,NS,pg)
!
!-----------------------------------------------------------------------------
!=========================================================================
subroutine calc_inttwoele_mat(inttwoele)
! Sum_alpha,beta
! int dr ds (psi^alpha(r)^+)_na (psi^alpha(r))_mb (1/|r-s|)(psi^beta(s)^+)_pc (psi^beta(s))_qd 
! 2012.5.19
!=========================================================================  
  use Precision
  use DefineTypes
  use DiracOutput
  use Constants
  implicit none
  
  complex(kind=dp),intent(out) :: inttwoele(-NBS:NBS,2,-NBS:NBS,2,-NBS:NBS,2,-NBS:NBS,2) 
  ! (n^a m^b|p^c q^d)
  ! negative values account for Kramers partners. 0 is not used.

  type(primitive_gaussian) :: pg
  complex(kind=dp) :: c_pg(4,NMAX_PG)

  complex(kind=dp) :: inttwoele_pg_mat(NBS_S,NBS_S,NBS_S,NBS_S,4,4) !(ij|kl)^ab

  integer :: i,j,k,l ! p.g. indice
  integer :: n,m,p,q,a,b,c,d ! MO indice 
  integer :: alpha,beta ! spinor indice

  complex(kind=dp) :: tmp1(NBS_S,NBS_S,NBS_S,-NBS:NBS,2,4,4) !(i j|k q^d,alpha,beta)
  complex(kind=dp) :: tmp2(NBS_S,NBS_S,-NBS:NBS,2,-NBS:NBS,2,4,4) !(i j|p^c q^d,alpha,beta)
  complex(kind=dp) :: tmp3(NBS_S,-NBS:NBS,2,-NBS:NBS,2,-NBS:NBS,2,4,4) !(i m^b|p^c q^d,alpha,beta)
  complex(kind=dp) :: tmp4(-NBS:NBS,2,-NBS:NBS,2,-NBS:NBS,2,-NBS:NBS,2,4,4) !(n^a m^b|p^c q^d,alpha,beta)
  complex(kind=dp) :: sum

  call copy_DiracOutput_pg(pg) ! get pg

  ! (i j|k l)^{alpha,beta}
  call calc_inttwoele_pg_mat(NBS_L,NBS_S,pg,inttwoele_pg_mat) 

  ! (i j|k q^d)^{alpha,beta}
  do alpha=1,4
     do beta=1,4
        do i=1,NBS_S
           do j=1,NBS_S
              do k=1,NBS_S
                 do q = -NBS,NBS  !q=0 is not used.
                    do d = 1,2
                       sum = (0._dp,0._dp)
                       call copy_DiracOutput_cp(q,d,c_pg)
                       do l=1,NBS_S
                          sum = sum +c_pg(beta,l)*inttwoele_pg_mat(i,j,k,l,alpha,beta)
                          ! c^beta_{q^d l} x (ij|kl)^{alpha beta}
                       end do
                       tmp1(i,j,k,q,d,alpha,beta) = sum
                    end do
                 end do
              end do
           end do
        end do
     end do
  end do

  ! (i j|p^c q^d)^{alpha,beta}
  do alpha=1,4
     do beta=1,4
        do i=1,NBS_S
           do j=1,NBS_S
              do q = -NBS,NBS  !q=0 is not used.
                 do d = 1,2
                    do p = -NBS,NBS  !p=0 is not used.
                       do c = 1,2
                          sum = (0._dp,0._dp)
                          call copy_DiracOutput_cp(p,c,c_pg)
                          do k=1,NBS_S
                             sum = sum +conjg(c_pg(beta,k))*tmp1(i,j,k,q,d,alpha,beta)
                             ! (c^beta_{p^c k})^* x (ij|k q^d)^{alpha beta}
                          end do
                          tmp2(i,j,p,c,q,d,alpha,beta) = sum
                       end do
                    end do
                 end do
              end do
           end do
        end do
     end do
  end do

  ! (i m^b|p^c q^d)^{alpha,beta}
  do alpha=1,4
     do beta=1,4
        do i=1,NBS_S
           do q = -NBS,NBS  !q=0 is not used.
              do d = 1,2
                 do p = -NBS,NBS  !p=0 is not used.
                    do c = 1,2
                       do m = -NBS,NBS  !m=0 is not used.
                          do b = 1,2
                             sum = (0._dp,0._dp)
                             call copy_DiracOutput_cp(m,b,c_pg)
                             do j=1,NBS_S
                                sum = sum +c_pg(alpha,j)*tmp2(i,j,p,c,q,d,alpha,beta)
                                ! c^alpha_{m^b j} x (i j|p^c q^d)^{alpha beta}
                             end do
                             tmp3(i,m,b,p,c,q,d,alpha,beta) = sum
                          end do
                       end do
                    end do
                 end do
              end do
           end do
        end do
     end do
  end do

  ! (n^a m^b|p^c q^d)^{alpha,beta}
  do alpha=1,4
     do beta=1,4
        do q = -NBS,NBS  !q=0 is not used.
           do d = 1,2
              do p = -NBS,NBS  !p=0 is not used.
                 do c = 1,2
                    do m = -NBS,NBS  !m=0 is not used.
                       do b = 1,2
                          do n = -NBS,NBS  !n=0 is not used.
                             do a = 1,2
                                sum = (0._dp,0._dp)
                                call copy_DiracOutput_cp(n,a,c_pg)
                                do i=1,NBS_S
                                   sum = sum +conjg(c_pg(alpha,i))*tmp3(i,m,b,p,c,q,d,alpha,beta)
                                   ! (c^alpha_{n^a i})^* x (i m^b|p^c q^d)^{alpha beta}
                                end do
                                tmp4(n,a,m,b,p,c,q,d,alpha,beta) = sum
                             end do
                          end do
                       end do
                    end do
                 end do
              end do
           end do
        end do
     end do
  end do

  ! (n^a m^b|p^c q^d)
  do q = -NBS,NBS  !q=0 is not used.
     do d = 1,2
        do p = -NBS,NBS  !p=0 is not used.
           do c = 1,2
              do m = -NBS,NBS  !m=0 is not used.
                 do b = 1,2
                    do n = -NBS,NBS  !n=0 is not used.
                       do a = 1,2
                          
                          sum = (0._dp,0._dp)
                          do alpha=1,4
                             do beta=1,4
                                sum = sum +tmp4(n,a,m,b,p,c,q,d,alpha,beta)
                             end do
                          end do
                          inttwoele(n,a,m,b,p,c,q,d) = sum

                       end do
                    end do
                 end do
              end do
           end do
        end do
     end do
  end do

  return
end subroutine calc_inttwoele_mat


!=========================================================================
subroutine calc_inttwoele_pg_mat(NL,NS,pg,inttwoele_pg_mat)
! Int dr ds (g^a(r)^+)_i (g^a(r))_j (1/|r-s|)(g^b(s)^+)_k (g^b(s))_l 
! Note that g depends on large or small components
! g^a = g^L if a=1,2
!     = g^S if a=3,4
! compute only necessary spinor combinations
! 4c integration for primitive gaussians.
! 2012.5.19
!=========================================================================  
  use Precision
  use DefineTypes
  implicit none
  
  integer,intent(in) :: NL,NS
  type(primitive_gaussian),intent(in) :: pg
  complex(kind=dp),intent(out) :: inttwoele_pg_mat(NS,NS,NS,NS,4,4) !(ij|kl)^ab
  
  integer :: i,j,k,l ! p.g. indice
  integer :: a,b ! spinor indice
  integer :: NLS
  complex(kind=dp) :: inttwoele_pg

  inttwoele_pg_mat = (0._dp,0._dp)

  do a=1,4
     do b=1,4
        do i=1,NLS(a,NL,NS)
           do j=1,NLS(a,NL,NS)
              do k=1,NLS(b,NL,NS)
                 do l=1,NLS(b,NL,NS)
                    inttwoele_pg_mat(i,j,k,l,a,b) = inttwoele_pg(i,j,k,l,a,a,b,b,NL,NS,pg)
                 end do
              end do
           end do
        end do
     end do
  end do

  return
end subroutine calc_inttwoele_pg_mat


!=========================================================================
function inttwoele_pg(i,j,k,l,s_i,s_j,s_k,s_l,NL,NS,pg)
! Int dr ds (g^si(r)^+)_i (g^sj(r))_j (1/|r-s|)(g^sk(s)^+)_k (g^sl(s))_l 
! 4c integration for primitive gaussians.
! Note that g depends on large or small components
! g^si = g^L if si=1,2
!      = g^S if si=3,4
! 2012.5.19
!=========================================================================  
  use Precision
  use DefineTypes
  implicit none

  complex(kind=dp) :: inttwoele_pg

  integer,intent(in) :: i,j,k,l ! p.g. indice
  integer,intent(in) :: s_i,s_j,s_k,s_l ! spinor indice
  integer,intent(in) :: NL,NS
  type(primitive_gaussian),intent(in) :: pg
  
  integer :: NLS
  real(kind=dp) :: twoele
  real(kind=dp) :: posA(3), posB(3), posC(3), posD(3) ! position of center of gaussian 
  real(kind=dp) :: alphaA, alphaB, alphaC, alphaD ! exponents
  integer :: n1(3),nbar1(3),n2(3),nbar2(3)  
  real(kind=dp) :: norm_pg  ! function
  real(kind=dp) :: f_norm_pg

  call set_pg(s_i,i,NL,NS,pg,posA,alphaA,n1(1),n1(2),n1(3))
  call set_pg(s_j,j,NL,NS,pg,posB,alphaB,nbar1(1),nbar1(2),nbar1(3))
  call set_pg(s_k,k,NL,NS,pg,posC,alphaC,n2(1),n2(2),n2(3))
  call set_pg(s_l,l,NL,NS,pg,posD,alphaD,nbar2(1),nbar2(2),nbar2(3))
  call gauss_int_twoele(posA,alphaA,n1,posB,alphaB,nbar1,posC,alphaC,n2,posD,alphaD,nbar2,twoele)  ! (AB|CD)

  f_norm_pg = norm_pg(alphaA,n1(1),n1(2),n1(3))*norm_pg(alphaB,nbar1(1),nbar1(2),nbar1(3)) &
       *norm_pg(alphaC,n2(1),n2(2),n2(3))*norm_pg(alphaD,nbar2(1),nbar2(2),nbar2(3))

  inttwoele_pg = twoele*f_norm_pg

  return
end function inttwoele_pg


!=========================================================================
!-----------------------------------------------------------------------------
! This subroutine will not be used because it is too slow. 
!-----------------------------------------------------------------------------
subroutine calc_inttwoele_nmpq(in,im,ip,iq,NL,NS,pg,c_n,c_m,c_p,c_q,int_nmpq)
! int dr ds (psi(r)^+)_in (psi(r))_im (1/|r-s|)(psi(s)^+)_ip (psi(s))_iq 
!-----------------------------------------------------------------------------
! Routines to calculate two electron integrals.
! Coefficients are summed but spinor indices are left unsummed.
!-----------------------------------------------------------------------------
!=========================================================================  
  use Precision
  use DefineTypes
  implicit none

  integer,intent(in) :: in,im,ip,iq ! spinor indice
  integer,intent(in) :: NL,NS
  type(primitive_gaussian),intent(in) :: pg
  complex(kind=dp),intent(in) :: c_n(4,NMAX_PG),c_m(4,NMAX_PG)
  complex(kind=dp),intent(in) :: c_p(4,NMAX_PG),c_q(4,NMAX_PG)
  
  complex(kind=dp),intent(out) :: int_nmpq
  
  integer :: NLS
  integer :: i,j,k,l
  real(kind=dp) :: twoele
  real(kind=dp) :: posA(3), posB(3), posC(3), posD(3) ! position of center of gaussian 
  real(kind=dp) :: alphaA, alphaB, alphaC, alphaD ! exponents
  integer :: n1(3),nbar1(3),n2(3),nbar2(3)  
  real(kind=dp) :: norm_pg  ! function
  real(kind=dp) :: f_norm_pg


  int_nmpq = (0._dp,0._dp)

  do i=1,NLS(in,NL,NS)
     do j=1,NLS(im,NL,NS)
        do k=1,NLS(ip,NL,NS)
           do l=1,NLS(iq,NL,NS)
              call set_pg(in,i,NL,NS,pg,posA,alphaA,n1(1),n1(2),n1(3))
              call set_pg(im,j,NL,NS,pg,posB,alphaB,nbar1(1),nbar1(2),nbar1(3))
              call set_pg(ip,k,NL,NS,pg,posC,alphaC,n2(1),n2(2),n2(3))
              call set_pg(iq,l,NL,NS,pg,posD,alphaD,nbar2(1),nbar2(2),nbar2(3))
              call gauss_int_twoele(posA,alphaA,n1,posB,alphaB,nbar1,posC,alphaC,n2,posD,alphaD,nbar2,twoele)  ! (AB|CD)
              f_norm_pg = norm_pg(alphaA,n1(1),n1(2),n1(3))*norm_pg(alphaB,nbar1(1),nbar1(2),nbar1(3)) &
                   *norm_pg(alphaC,n2(1),n2(2),n2(3))*norm_pg(alphaD,nbar2(1),nbar2(2),nbar2(3))
              int_nmpq = int_nmpq + conjg(c_n(in,i))*c_m(im,j)*conjg(c_p(ip,k))*c_q(iq,l)*twoele*f_norm_pg
           end do
        end do
     end do
  end do
  return
end subroutine calc_inttwoele_nmpq


!-----------------------------------------------------------------------------
! Routines to calculate gaussian integration for spinor bilinears.
! Coefficients are summed but spinor indices are left unsummed.
!
! These routines should be modified to faster routines.(121031)
!-----------------------------------------------------------------------------
!=========================================================================
subroutine calc_intF_pq(ip,iq,vecP,NL,NS,pg,c_p,c_q,intF_pq)
! int dr (psi^+)_ip (psi)_iq *exp(I*p.r/hbar)  (hbar=1)
!=========================================================================  
  use Precision
  use DefineTypes
  implicit none

  integer,intent(in) :: ip,iq ! spinor indice
  real(kind=dp),intent(in) :: vecP(3) ! wave number for fourier transformation  
  integer,intent(in) :: NL,NS
  type(primitive_gaussian),intent(in) :: pg
  complex(kind=dp),intent(in) :: c_p(4,NMAX_PG),c_q(4,NMAX_PG)
  
  complex(kind=dp),intent(out) :: intF_pq
  
  integer :: NLS
  integer :: i,j,k
  complex(kind=dp) :: fourier
  real(kind=dp) :: posA(3), posB(3) ! position of center of gaussian 
  real(kind=dp) :: alphaA, alphaB ! exponents
  integer :: nx,nbarx,ny,nbary,nz,nbarz
  real(kind=dp) :: norm_pg  ! function
  real(kind=dp) :: f_norm_pg


  intF_pq = (0._dp,0._dp)

  do i=1,NLS(ip,NL,NS)
     do j=1,NLS(iq,NL,NS)
        call set_pg(ip,i,NL,NS,pg,posA,alphaA,nx,ny,nz)
        call set_pg(iq,j,NL,NS,pg,posB,alphaB,nbarx,nbary,nbarz)
        call gauss_int_fourier(posA,alphaA,nx,ny,nz,posB,alphaB,nbarx,nbary,nbarz,vecP,fourier)
        f_norm_pg = norm_pg(alphaA,nx,ny,nz)*norm_pg(alphaB,nbarx,nbary,nbarz)
        intF_pq = intF_pq + conjg(c_p(ip,i))*c_q(iq,j)*fourier*f_norm_pg
     end do
  end do
  return
end subroutine calc_intF_pq


!=========================================================================
subroutine calc_intE_pq(ip,iq,posR,NL,NS,pg,c_p,c_q,intE_pq)
! Int dr (psi^+)_ip (psi)_iq *vec(r-R)/|r-R|^3
! calculate three components 
!=========================================================================  
  use Precision
  use DefineTypes
  implicit none

  integer,intent(in) :: ip,iq ! spinor indice
  real(kind=dp),intent(in) :: posR(3)  
  integer,intent(in) :: NL,NS
  type(primitive_gaussian),intent(in) :: pg
  complex(kind=dp),intent(in) :: c_p(4,NMAX_PG),c_q(4,NMAX_PG)
  
  complex(kind=dp),intent(out) :: intE_pq(3)
  
  integer :: NLS
  integer :: i,j,k
  real(kind=dp) :: overlap,ef(3) 
  real(kind=dp) :: nucatt
  real(kind=dp) :: posA(3), posB(3) ! position of center of gaussian 
  real(kind=dp) :: alphaA, alphaB ! exponents
  integer :: nx,nbarx,ny,nbary,nz,nbarz
  real(kind=dp) :: norm_pg  ! function
  real(kind=dp) :: f_norm_pg


  intE_pq(1:3) = (0._dp,0._dp)

  do i=1,NLS(ip,NL,NS)
     do j=1,NLS(iq,NL,NS)
        call set_pg(ip,i,NL,NS,pg,posA,alphaA,nx,ny,nz)
        call set_pg(iq,j,NL,NS,pg,posB,alphaB,nbarx,nbary,nbarz)
        call gauss_int(posA,alphaA,nx,ny,nz,posB,alphaB,nbarx,nbary,nbarz,posR,overlap,nucatt,ef)
        f_norm_pg = norm_pg(alphaA,nx,ny,nz)*norm_pg(alphaB,nbarx,nbary,nbarz)
        do k=1,3
           intE_pq(k) = intE_pq(k) + conjg(c_p(ip,i))*c_q(iq,j)*ef(k)*f_norm_pg
        end do
     end do
  end do
  return
end subroutine calc_intE_pq

!=========================================================================
subroutine calc_intV_pq(ip,iq,posRA,NL,NS,pg,c_p,c_q,intV_pq)
! Int (psi^+)_ip (psi)_iq /|r-R|
!=========================================================================  
  use Precision
  use DefineTypes
  implicit none

  integer,intent(in) :: ip,iq ! spinor indice
  real(kind=dp),intent(in) :: posRA(3)  ! position of nucleus
  integer,intent(in) :: NL,NS
  type(primitive_gaussian),intent(in) :: pg
  complex(kind=dp),intent(in) :: c_p(4,NMAX_PG),c_q(4,NMAX_PG)
  
  complex(kind=dp),intent(out) :: intV_pq
  
  integer :: NLS
  integer :: i,j
  real(kind=dp) :: overlap,ef(3) ! --> dummy variables here 
  real(kind=dp) :: nucatt
  real(kind=dp) :: posA(3), posB(3) ! position of center of gaussian 
  real(kind=dp) :: alphaA, alphaB ! exponents
  integer :: nx,nbarx,ny,nbary,nz,nbarz
  real(kind=dp) :: norm_pg  ! function
  real(kind=dp) :: f_norm_pg


  intV_pq = (0._dp,0._dp)

  do i=1,NLS(ip,NL,NS)
     do j=1,NLS(iq,NL,NS)
        call set_pg(ip,i,NL,NS,pg,posA,alphaA,nx,ny,nz)
        call set_pg(iq,j,NL,NS,pg,posB,alphaB,nbarx,nbary,nbarz)
        call gauss_int(posA,alphaA,nx,ny,nz,posB,alphaB,nbarx,nbary,nbarz,posRA,overlap,nucatt,ef)
        f_norm_pg = norm_pg(alphaA,nx,ny,nz)*norm_pg(alphaB,nbarx,nbary,nbarz)
        intV_pq = intV_pq + conjg(c_p(ip,i))*c_q(iq,j)*nucatt*f_norm_pg
     end do
  end do
  return
end subroutine calc_intV_pq

!=========================================================================
subroutine calc_intT_pq(ip,iq,NL,NS,pg,c_p,c_q,intT_pq)
! int (psi^+)_ip (d/dx_i psi)_iq
!=========================================================================  
  use Precision
  use DefineTypes
  implicit none

  integer,intent(in) :: ip,iq ! spinor indice
  integer,intent(in) :: NL,NS
  type(primitive_gaussian),intent(in) :: pg
  complex(kind=dp),intent(in) :: c_p(4,NMAX_PG),c_q(4,NMAX_PG)
  
  complex(kind=dp),intent(out) :: intT_pq(3)
  
  integer :: NLS
  integer :: i,j,k
  real(kind=dp) :: grad(3)
  real(kind=dp) :: posA(3), posB(3) ! position of center of gaussian 
  real(kind=dp) :: alphaA, alphaB ! exponents
  integer :: nx,nbarx,ny,nbary,nz,nbarz
  real(kind=dp) :: norm_pg  ! function
  real(kind=dp) :: f_norm_pg

  intT_pq(1:3) = (0._dp,0._dp)

  do i=1,NLS(ip,NL,NS)
     do j=1,NLS(iq,NL,NS)
        call set_pg(ip,i,NL,NS,pg,posA,alphaA,nx,ny,nz)
        call set_pg(iq,j,NL,NS,pg,posB,alphaB,nbarx,nbary,nbarz)
        call gauss_int_relKE(posA,alphaA,nx,ny,nz,posB,alphaB,nbarx,nbary,nbarz,grad)
        f_norm_pg = norm_pg(alphaA,nx,ny,nz)*norm_pg(alphaB,nbarx,nbary,nbarz)
        do k=1,3
           intT_pq(k) = intT_pq(k) + conjg(c_p(ip,i))*c_q(iq,j)*grad(k)*f_norm_pg
        end do
     end do
  end do
  return
end subroutine calc_intT_pq

!=========================================================================
subroutine calc_intN_pq(ip,iq,NL,NS,pg,c_p,c_q,intN_pq)
! int (psi^+)_ip (psi)_iq 
!=========================================================================  
  use Precision
  use DefineTypes
  implicit none

  integer,intent(in) :: ip,iq ! spinor indice
  integer,intent(in) :: NL,NS
  type(primitive_gaussian),intent(in) :: pg
  complex(kind=dp),intent(in) :: c_p(4,NMAX_PG),c_q(4,NMAX_PG)
  
  complex(kind=dp),intent(out) :: intN_pq
  
  integer :: NLS
  integer :: i,j
  real(kind=dp) :: overlap
  real(kind=dp) :: posA(3), posB(3) ! position of center of gaussian 
  real(kind=dp) :: alphaA, alphaB ! exponents
  integer :: nx,nbarx,ny,nbary,nz,nbarz
  real(kind=dp) :: norm_pg  ! function
  real(kind=dp) :: f_norm_pg

  intN_pq = (0._dp,0._dp)

  do i=1,NLS(ip,NL,NS)
     do j=1,NLS(iq,NL,NS)
        call set_pg(ip,i,NL,NS,pg,posA,alphaA,nx,ny,nz)
        call set_pg(iq,j,NL,NS,pg,posB,alphaB,nbarx,nbary,nbarz)
        call gauss_int_overlap(posA,alphaA,nx,ny,nz,posB,alphaB,nbarx,nbary,nbarz,overlap)
        f_norm_pg = norm_pg(alphaA,nx,ny,nz)*norm_pg(alphaB,nbarx,nbary,nbarz)
        intN_pq = intN_pq + conjg(c_p(ip,i))*c_q(iq,j)*overlap*f_norm_pg
     end do
  end do
  return
end subroutine calc_intN_pq

!-----------------------------------------------------------------------------
! utility routines
!-----------------------------------------------------------------------------
!=============================================================
subroutine set_pg(i,j,NL,NS,pg,pos,alpha,nx,ny,nz)
! i:1,2 -> large, i:3,4-> small
! jth L or S primitive gaussian 
!=============================================================
  use Precision
  use DefineTypes
  implicit none
  
  integer,intent(in) :: i,j,NL,NS
  type(primitive_gaussian),intent(in) :: pg

  real(kind=dp),intent(out) :: pos(3),alpha
  integer,intent(out) :: nx,ny,nz

  if(i.eq.1 .or. i.eq.2) then  ! Large component
     pos(1)=pg%xL(j); pos(2)=pg%yL(j); pos(3)=pg%zL(j); alpha=pg%aL(j); nx=pg%nxL(j); ny=pg%nyL(j); nz=pg%nzL(j)
  elseif(i.eq.3 .or. i.eq.4) then
     pos(1)=pg%xS(j); pos(2)=pg%yS(j); pos(3)=pg%zS(j); alpha=pg%aS(j); nx=pg%nxS(j); ny=pg%nyS(j); nz=pg%nzS(j)
  else
     write(*,*) "i should be 1,2,3,4 in subroutine set_pg."
     stop
  end if

  return
end subroutine set_pg
  

!=======================
function NLS(i,NL,NS)
!=======================
  use Precision
  implicit none
  integer :: NLS
  integer,intent(in) :: i,NL,NS

  if(i.eq.1 .or. i.eq.2) then
     NLS = NL
  elseif(i.eq.3 .or. i.eq.4) then
     NLS = NS
  else
     write(*,*) "i should be 1,2,3,4 in function NLS."
     stop
  end if

  return
end function NLS

!=============================================================
subroutine convert_label(j,NX,NY,NZ,ix,iy,iz)
!  Convert 1D label (j) to 3D label (ix,iy,iz)
! 3D label range is ix=(1:NX) etc.
! j = (iz-1)*NX*NY +(iy-1)*NX +ix
!=============================================================
  use Precision
  implicit none
  
  integer,intent(in)  :: j,NX,NY,NZ
  integer,intent(out) :: ix,iy,iz

  integer :: k

  if(j.lt.1 .or. j.gt.NX*NY*NZ) then
     write(*,*) "j should be 1~NX*NY*NZ."
     stop
  end if

  ix = mod(j,NX)
  if(ix.eq.0) then
     ix = NX
  end if

  k = (j-ix)/NX +1
  iy = mod(k,NY)
  if(iy.eq.0) then
     iy = NY
  end if

  iz = (j-ix-(iy-1)*NX)/NX/NY +1

  return
end subroutine convert_label

!=============================================================
subroutine convert_label_2(j,NX,NY,NZ,ix,iy,iz,sig)
!  Convert 1D label (j) to 3D label (ix,iy,iz) and sig (+/-)
! 3D label range is ix=(1:NX) etc.
! j = (iz-1)*NX*NY +(iy-1)*NX +ix for sig = +
! j = (iz-1)*NX*NY +(iy-1)*NX +ix +NX*NY*NZ for sig = -
! j=1~2*NX*NY*NZ
!=============================================================
  use Precision
  implicit none
  
  integer,intent(in)  :: j,NX,NY,NZ
  integer,intent(out) :: ix,iy,iz
  character(LEN=1),intent(out) :: sig

  integer :: k,NXYZ

  NXYZ = NX*NY*NZ
  if(j.ge.1 .and. j.le.NXYZ) then
     sig = "+"
     k = j
     call convert_label(k,NX,NY,NZ,ix,iy,iz)
  elseif(j.ge.(NXYZ+1) .and. j.le.(2*NXYZ)) then
     sig = "-"
     k = j-NXYZ
     call convert_label(k,NX,NY,NZ,ix,iy,iz)
  else
     write(*,*) "j should be 1<= j <=2*NX*NY*NZ."
     stop
  end if
  
  return
end subroutine convert_label_2
