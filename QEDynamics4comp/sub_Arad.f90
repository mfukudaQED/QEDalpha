! Last Change:26-Nov-2011.
!  function Arad_vec(t,i,x,y,z)     <coherent|Arad^i(t)|coherent>
!  function tau_Arad_mat(t,k,l,x,y,z,p,a,q,b)
!  function t_Arad_mat(t,i,x,y,z,p,a,q,b)
!  function tau_Arad_Qmat(t,k,l,x,y,z,nn,mm)
!  function t_Arad_Qmat(t,i,x,y,z,nn,mm)
!============================================================
function Arad_vec(t,i,x,y,z)
! written by Fukuda 17-Oct-2011
! <coherent|Arad^i(t)|coherent>
! t : time = timestep*DeltaT
! p (photon momentum) specified in spherical coordinate (p0,th,phi)
! sig (+/-) specifies rotating direction of circular polarization (which depends on th and phi)
!============================================================
  use Constants  ! use P0MAX,NP0,NPTH,NPPHI,Nph,ALPHA_COH
  implicit none

  complex(kind=8) :: Arad_vec
  complex(kind=8) :: Aradj_vec
  real(kind=8),intent(in) :: t  ! time
  integer,intent(in) :: i
  integer :: k
  real(kind=8),intent(in) :: x,y,z
  integer :: j ! collective index which expresees photon momentum and polarization

  real(kind=8) :: p0,th,phi ! photon momentum vector in spherical coord.
  character(LEN=1) :: sig ! polarization (+ or -)  
  integer :: ip0,ith,iphi
  real(kind=8) :: dp0,dth,dphi ! mesh width

  real(kind=8) :: vecP(3) ! photon momentum vector in Cartesian coord.
  complex(kind=8) :: vec_pol(3) ! (circular) polarization vector (depends on sig)
  complex(kind=8) :: tmpA
  real(kind=8) :: DeltaPj,P0j

  dp0  = P0MAX/NP0 ! P0MAX specified in GammaMatrix. We do not use p0=0
  dth  = PI/(NPTH-1) ! NPTH specified in GammaMatrix. th include both 0 and PI
  dphi = 2.d0*PI/NPPHI ! NPPHI specified in GammaMatrix. phi include 0 but not 2PI

  Aradj_vec = (0d0,0d0)

  do j=1,Nph

    call convert_label_2(j,NP0,NPTH,NPPHI,ip0,ith,iphi,sig)

    p0 = dp0*ip0
    th = dth*(ith-1)  ! start from 0
    phi= dphi*(iphi-1) ! start from 0

!    write(*,*) ip0,ith,iphi
!    write(*,*) p0,th/PI,phi/PI,sig

    vecP(1) = p0*sin(th)*cos(phi)
    vecP(2) = p0*sin(th)*sin(phi)
    vecP(3) = p0*cos(th)
!    write(*,*)'vecP:',vecP(1),vecP(2),vecP(3)

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
       vec_pol(k) = vec_pol(k)/sqrt(2.d0)
!       write(*,*)'vec_pol:',vec_pol(k)
    end do
    
    call calc_DeltaPj_and_P0j(j,DeltaPj,P0j)

    tmpA = vec_pol(i)*ALPHA_COH(j)*exp(IU*(vecP(1)*x +vecP(2)*y +vecP(3)*z -CCC*p0*t))
    Aradj_vec = Aradj_vec + DeltaPj*(tmpA + conjg(tmpA))
!    write(*,*)'ALPHA_COH(',j,')=',ALPHA_COH(j)
!    write(*,*)'tmpA:',tmpA
!    write(*,*)'DeltaPj(',j,')=',DeltaPj
!    write(*,*) j, Aradj_vec
!    Aradj_vec = Aradj_vec + DeltaPj * (exp(IU*(vecP(1)*x +vecP(2)*y +vecP(3)*z -CCC*p0*t)) *vec_pol(i) *ALPHA_COH(j) &
!                                        +exp(-IU*(vecP(1)*x +vecP(2)*y +vecP(3)*z -CCC*p0*t)) *conjg(vec_pol(i)) *conjg(ALPHA_COH(j)))
    
  end do

  Arad_vec = Aradj_vec

  return
end function Arad_vec

!============================================================
function tau_Arad_mat(t,k,l,x,y,z,p,a,q,b)
! written by Fukuda 17-Oct-2011
! t : time = timestep*DeltaT
! p,q : labels for molecular orbitals.
! a,b : labels for electron/positron ("+"/"-")
!============================================================

  use Constants  ! use CCC, Nph
  use DefineTypes
  use DiracOutput
  implicit none

  complex(kind=8) :: tau_Arad_mat
  real(kind=8),intent(in) :: t  ! time
  integer,intent(in) :: k,l
  integer,intent(in) :: p,q
  character(LEN=1),intent(in) :: a,b
  real(kind=8),intent(in) :: x,y,z

  complex(kind=8) :: j_mat,tau_mat,Arad_vec

   tau_Arad_mat = tau_mat(k,l,x,y,z,p,a,q,b) + j_mat(l,x,y,z,p,a,q,b)*Arad_vec(t,k,x,y,z)/CCC

  return
end function tau_Arad_mat

!============================================================
function t_Arad_mat(t,i,x,y,z,p,a,q,b)
! written by Fukuda 17-Oct-2011
! t : time = timestep*DeltaT
! i-th component of t(ab)_pq, i=1,3
! p,q : labels for molecular orbitals.
! a,b : labels for electron/positron ("+"/"-")
!============================================================
  use DefineTypes
  use DiracOutput
  use Constants  ! use Nph
  implicit none
 
  complex(kind=8) :: t_Arad_mat
  real(kind=8),intent(in) :: t  ! time
  integer,intent(in) :: i
  real(kind=8),intent(in) :: x,y,z
  integer,intent(in) :: p,q
  character(LEN=1),intent(in) :: a,b

  complex(kind=8) :: tau_Arad_mat

  if(i.lt.1 .or. i.gt.3) then
     write(*,*) "i should be 1-3 in t_mat."
     stop
  end if

  if(i.eq.1) then
     t_Arad_mat = -tau_Arad_mat(t,2,3,x,y,z,p,a,q,b) +tau_Arad_mat(t,3,2,x,y,z,p,a,q,b)
  elseif(i.eq.2) then
     t_Arad_mat =  tau_Arad_mat(t,1,3,x,y,z,p,a,q,b) -tau_Arad_mat(t,3,1,x,y,z,p,a,q,b)
  elseif(i.eq.3) then
     t_Arad_mat = -tau_Arad_mat(t,1,2,x,y,z,p,a,q,b) +tau_Arad_mat(t,2,1,x,y,z,p,a,q,b)
  else
     write(*,*) "i should be 1-3 in t_mat."
     stop
  end if
  
  return
end function t_Arad_mat

!============================================================
function tau_Arad_Qmat(t,k,l,x,y,z,nn,mm)
! written by Fukuda 4-Nov-2011
! t : time = timestep*DeltaT
! nn,mm : labels for molecular orbitals (KP) and +/-. 1~4*NBS
!============================================================
  use DiracOutput
  use Constants  ! use Nph
  implicit none

  complex(kind=8) :: tau_Arad_Qmat
  real(kind=8),intent(in) :: t  ! time = timestep * DeltaT
  integer,intent(in) :: k,l
  real(kind=8),intent(in) :: x,y,z
  integer,intent(in) :: nn,mm

  complex(kind=8) :: tau_Arad_mat
  integer :: n,m ! index for 1~2*NBS (+ or - specified)
  character(LEN=1) :: a,b

  call index_from_Qmat(nn,n,a)
  call index_from_Qmat(mm,m,b)

  tau_Arad_Qmat = tau_Arad_mat(t,k,l,x,y,z,n,a,m,b)

  return
end function tau_Arad_Qmat

!============================================================
function t_Arad_Qmat(t,i,x,y,z,nn,mm)
! written by Fukuda 4-Nov-2011
! t : time = timestep*DeltaT
! nn,mm : labels for molecular orbitals (KP) and +/-. 1~4*NBS
!============================================================
  use DiracOutput
  use Constants  ! use Nph
  implicit none

  complex(kind=8) :: t_Arad_Qmat
  real(kind=8),intent(in) :: t  ! time = timestep * DeltaT
  integer,intent(in) :: i
  real(kind=8),intent(in) :: x,y,z
  integer,intent(in) :: nn,mm

  complex(kind=8) :: t_Arad_mat
  integer :: n,m ! index for 1~2*NBS (+ or - specified)
  character(LEN=1) :: a,b
  
  call index_from_Qmat(nn,n,a)
  call index_from_Qmat(mm,m,b)

  t_Arad_Qmat = t_Arad_mat(t,i,x,y,z,n,a,m,b)

  return
end function t_Arad_Qmat

