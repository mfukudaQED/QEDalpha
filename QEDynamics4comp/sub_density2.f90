!=============================================================================
!
!
!  2010.11.11
! -subroutine for n_pq^(ab)  (function density_mat)
!
!  11.12
! -subroutine for j^k div j
!
!  11.13
! -subroutine for copying electron/positorn coefficient matrix to coefficient vector.
!  (used in  function ..._mat) copy_DiracOutput
! -stress tensor
! 
!  11.16  sub_density2.f90
! -define type primitive_gaussian (in module DefineTypes)
! -rewrite subroutines using this type.
! -spin density
!
!  11.18
! -spintorque and zeta force
!
! 2011.1.10
! -density_ele etc. -> added NEL=odd case
!
! 10.31
! -rho (charge density) = Ze x n_pq^(ab)
!
! 2012.5.19
! -copy_DiracOutput_pg, copy_DiracOutput_cp
!
!============================================================
function zeta_ele(i,x,y,z,NEL)
! i-th component of electron zeta force density
!============================================================
  implicit none

  complex(kind=8) :: zeta_ele
  integer,intent(in) :: i
  real(kind=8),intent(in) :: x,y,z
  integer,intent(in) :: NEL

  complex(kind=8) :: zeta_mat
  complex(kind=8) :: sum
  integer :: j
  
  sum = (0.d0,0.d0)
  if(mod(NEL,2).eq.0) then
     do j=1,NEL/2  ! assume Kramers-Restricted
        sum = sum +2.d0*zeta_mat(i,x,y,z,j,"+",j,"+")  ! Kramers pair
     end do
  else
     do j=1,(NEL-1)/2  
        sum = sum +2.d0*zeta_mat(i,x,y,z,j,"+",j,"+")  ! Kramers pair
     end do
     j=(NEL+1)/2
     sum = sum +zeta_mat(i,x,y,z,j,"+",j,"+")  
  end if

  zeta_ele = sum
  return
end function zeta_ele


!============================================================
function t_ele(i,x,y,z,NEL)
! i-th component of electron spintorque density
!============================================================
  implicit none

  complex(kind=8) :: t_ele
  integer,intent(in) :: i
  real(kind=8),intent(in) :: x,y,z
  integer,intent(in) :: NEL

  complex(kind=8) :: t_mat
  complex(kind=8) :: sum
  integer :: j

  sum = (0.d0,0.d0)
  if(mod(NEL,2).eq.0) then
     do j=1,NEL/2  ! assume Kramers-Restricted
        sum = sum +2.d0*t_mat(i,x,y,z,j,"+",j,"+")  ! Kramers pair
     end do
  else
     do j=1,(NEL-1)/2  
        sum = sum +2.d0*t_mat(i,x,y,z,j,"+",j,"+")  ! Kramers pair
     end do
     j=(NEL+1)/2
     sum = sum +t_mat(i,x,y,z,j,"+",j,"+")  
  end if
  
  t_ele = sum
  return
end function t_ele

!============================================================
function tau_ele(k,l,x,y,z,NEL)
! (k,l)-component of electronic stress tensor
!============================================================
  implicit none

  complex(kind=8) :: tau_ele
  integer,intent(in) :: k,l
  real(kind=8),intent(in) :: x,y,z
  integer,intent(in) :: NEL

  complex(kind=8) :: tau_mat
  complex(kind=8) :: sum
  integer :: j
  
  sum = (0.d0,0.d0)
  if(mod(NEL,2).eq.0) then
     do j=1,NEL/2  ! assume Kramers-Restricted
        sum = sum +2.d0*tau_mat(k,l,x,y,z,j,"+",j,"+")  ! Kramers pair
     end do
  else
     do j=1,(NEL-1)/2  
        sum = sum +2.d0*tau_mat(k,l,x,y,z,j,"+",j,"+")  ! Kramers pair
     end do
     j=(NEL+1)/2
     sum = sum +tau_mat(k,l,x,y,z,j,"+",j,"+")  
  end if
  tau_ele = sum
  return
end function tau_ele

!============================================================
function divj_ele(x,y,z,NEL)
! electron density
!============================================================
  implicit none

  real(kind=8) :: divj_ele
  real(kind=8),intent(in) :: x,y,z
  integer,intent(in) :: NEL

  complex(kind=8) :: divj_mat
  complex(kind=8) :: sum
  integer :: i
  
  sum = (0.d0,0.d0)
  do i=1,NEL/2  ! assume Kramers-Restricted
     sum = sum +2.d0*divj_mat(x,y,z,i,"+",i,"+")  ! Kramers pair
  end do

  divj_ele = real(sum)

  return
end function divj_ele

!============================================================
function j_ele(i,x,y,z,NEL)
! i-th component of electron current density
!============================================================
  implicit none

  real(kind=8) :: j_ele
  integer,intent(in) :: i
  real(kind=8),intent(in) :: x,y,z
  integer,intent(in) :: NEL

  complex(kind=8) :: j_mat
  complex(kind=8) :: sum
  integer :: j
  
  sum = (0.d0,0.d0)
  do j=1,NEL/2  ! assume Kramers-Restricted
     sum = sum +2.d0*j_mat(i,x,y,z,j,"+",j,"+")  ! Kramers pair
  end do

  j_ele = real(sum)

  return
end function j_ele

!============================================================
function s_ele(i,x,y,z,NEL)
! i-th component of electron current density
!============================================================
  implicit none

  real(kind=8) :: s_ele
  integer,intent(in) :: i
  real(kind=8),intent(in) :: x,y,z
  integer,intent(in) :: NEL

  complex(kind=8) :: s_mat
  complex(kind=8) :: sum
  integer :: j
  
  sum = (0.d0,0.d0)
  do j=1,NEL/2  ! assume Kramers-Restricted
     sum = sum +2.d0*s_mat(i,x,y,z,j,"+",j,"+")  ! Kramers pair
  end do

  s_ele = real(sum)

  return
end function s_ele

!============================================================
function rho_ele(x,y,z,NEL)
! electron charge density
!============================================================
  implicit none

  real(kind=8) :: rho_ele
  real(kind=8),intent(in) :: x,y,z
  integer,intent(in) :: NEL

  complex(kind=8) :: rho_mat
  complex(kind=8) :: sum
  integer :: i
  
  sum = (0.d0,0.d0)
  if(mod(NEL,2).eq.0) then
     do i=1,NEL/2  ! assume Kramers-Restricted
        sum = sum +2.d0*rho_mat(x,y,z,i,"+",i,"+")  ! Kramers pair
     end do
  else
     do i=1,(NEL-1)/2  
        sum = sum +2.d0*rho_mat(x,y,z,i,"+",i,"+")  ! Kramers pair
     end do
     i=(NEL+1)/2
     sum = sum +rho_mat(x,y,z,i,"+",i,"+")  
  end if

  rho_ele = real(sum)

  return
end function rho_ele

!============================================================
function density_ele(x,y,z,NEL)
! electron position probablility density
!============================================================
  implicit none

  real(kind=8) :: density_ele
  real(kind=8),intent(in) :: x,y,z
  integer,intent(in) :: NEL

  complex(kind=8) :: density_mat
  complex(kind=8) :: sum
  integer :: i
  
  sum = (0.d0,0.d0)
  if(mod(NEL,2).eq.0) then
     do i=1,NEL/2  ! assume Kramers-Restricted
        sum = sum +2.d0*density_mat(x,y,z,i,"+",i,"+")  ! Kramers pair
     end do
  else
     do i=1,(NEL-1)/2  
        sum = sum +2.d0*density_mat(x,y,z,i,"+",i,"+")  ! Kramers pair
     end do
     i=(NEL+1)/2
     sum = sum +density_mat(x,y,z,i,"+",i,"+")  
  end if

  density_ele = real(sum)

  return
end function density_ele

!----------------------------------------------------------------------------
!  zeta_mat    :zeta force density
!  t_mat       :spintorque density
!  tau_mat     :stress tensor density
!  divj_mat    :div of current density
!  j_mat       :current density
!  s_mat　　　  :spin density
!  rho_mat　　　:charge density
!  density_mat :position probability density
!----------------------------------------------------------------------------
!============================================================
function zeta_mat(i,x,y,z,p,a,q,b)
! i-th component of zeta(ab)_pq , i=1-3
! p,q : labels for molecular orbitals.
! a,b : labels for electron/positron ("+"/"-")
!============================================================
  use DefineTypes
  use DiracOutput
  implicit none
 
  complex(kind=8) :: zeta_mat
  integer,intent(in) :: i
  real(kind=8),intent(in) :: x,y,z
  integer,intent(in) :: p,q
  character(LEN=1),intent(in) :: a,b

  type(primitive_gaussian) :: pg
  complex(kind=8) :: c_p(4,NMAX_PG),c_q(4,NMAX_PG)

  if(i.lt.1 .or. i.gt.3) then
     write(*,*) "i should be 1-3 in zeta_mat."
     stop
  end if

  call copy_DiracOutput(p,a,q,b,pg,c_p,c_q)  ! set pg,c_p,c_q
  call calc_zeta_pq(i,x,y,z,NBS_L,NBS_S,pg,c_p,c_q,zeta_mat)
  
  return
end function zeta_mat

!============================================================
function t_mat(i,x,y,z,p,a,q,b)
! i-th component of t(ab)_pq, i=1,3
! p,q : labels for molecular orbitals.
! a,b : labels for electron/positron ("+"/"-")
!============================================================
  use DefineTypes
  use DiracOutput
  implicit none
 
  complex(kind=8) :: t_mat
  integer,intent(in) :: i
  real(kind=8),intent(in) :: x,y,z
  integer,intent(in) :: p,q
  character(LEN=1),intent(in) :: a,b

  type(primitive_gaussian) :: pg
  complex(kind=8) :: c_p(4,NMAX_PG),c_q(4,NMAX_PG)

  complex(kind=8) :: tau(3,3)
  integer :: k,l
  
  if(i.lt.1 .or. i.gt.3) then
     write(*,*) "i should be 1-3 in t_mat."
     stop
  end if

  call copy_DiracOutput(p,a,q,b,pg,c_p,c_q)  ! set pg,c_p,c_q
  do k=1,3
     do l=1,3
        call calc_tau_pq(k,l,x,y,z,NBS_L,NBS_S,pg,c_p,c_q,tau(k,l))
     end do
  end do
  if(i.eq.1) then
     t_mat = -tau(2,3)+tau(3,2)
  elseif(i.eq.2) then
     t_mat =  tau(1,3)-tau(3,1)
  elseif(i.eq.3) then
     t_mat = -tau(1,2)+tau(2,1)
  else
     write(*,*) "i should be 1-3 in t_mat."
     stop
  end if
  
  return
end function t_mat


!============================================================
function tau_mat(k,l,x,y,z,p,a,q,b)
! (k,l)-component of tau(ab)_pq , k,l=1-3
! p,q : labels for molecular orbitals.
! a,b : labels for electron/positron ("+"/"-")
!============================================================
  use DefineTypes
  use DiracOutput
  implicit none
 
  complex(kind=8) :: tau_mat
  integer,intent(in) :: k,l
  real(kind=8),intent(in) :: x,y,z
  integer,intent(in) :: p,q
  character(LEN=1),intent(in) :: a,b

  type(primitive_gaussian) :: pg
  complex(kind=8) :: c_p(4,NMAX_PG),c_q(4,NMAX_PG)

  if(k.lt.1 .or. k.gt.3 .or. l.lt.1 .or. l.gt.3) then
     write(*,*) "k and l should be 1-3 in tau_mat."
     stop
  end if

  call copy_DiracOutput(p,a,q,b,pg,c_p,c_q)  ! set pg,c_p,c_q
  call calc_tau_pq(k,l,x,y,z,NBS_L,NBS_S,pg,c_p,c_q,tau_mat)
  
  return
end function tau_mat

!============================================================
function divj_mat(x,y,z,p,a,q,b)
! p,q : labels for molecular orbitals.
! a,b : labels for electron/positron ("+"/"-")
!============================================================
  use DefineTypes
  use DiracOutput  
  implicit none
 
  complex(kind=8) :: divj_mat
  real(kind=8),intent(in) :: x,y,z
  integer,intent(in) :: p,q
  character(LEN=1),intent(in) :: a,b

  type(primitive_gaussian) :: pg
  complex(kind=8) :: c_p(4,NMAX_PG),c_q(4,NMAX_PG)

  call copy_DiracOutput(p,a,q,b,pg,c_p,c_q)  ! set pg,c_p,c_q
  call calc_divj_pq(x,y,z,NBS_L,NBS_S,pg,c_p,c_q,divj_mat)

  return
end function divj_mat


!============================================================
function j_mat(i,x,y,z,p,a,q,b)
! i-th component of j(ab)_pq , i=1-3
! p,q : labels for molecular orbitals.
! a,b : labels for electron/positron ("+"/"-")
!============================================================
  use DefineTypes
  use DiracOutput
  implicit none
 
  complex(kind=8) :: j_mat
  integer,intent(in) :: i
  real(kind=8),intent(in) :: x,y,z
  integer,intent(in) :: p,q
  character(LEN=1),intent(in) :: a,b

  type(primitive_gaussian) :: pg
  complex(kind=8) :: c_p(4,NMAX_PG),c_q(4,NMAX_PG)

  if(i.lt.1 .or. i.gt.3) then
     write(*,*) "i should be 1-3 in j_mat."
     stop
  end if

  call copy_DiracOutput(p,a,q,b,pg,c_p,c_q)  ! set pg,c_p,c_q
  call calc_j_pq(i,x,y,z,NBS_L,NBS_S,pg,c_p,c_q,j_mat)
  
  return
end function j_mat

!============================================================
function s_mat(i,x,y,z,p,a,q,b)
! i-th component of s(ab)_pq , i=1-3
! p,q : labels for molecular orbitals.
! a,b : labels for electron/positron ("+"/"-")
!============================================================
  use DefineTypes
  use DiracOutput
  implicit none
 
  complex(kind=8) :: s_mat
  integer,intent(in) :: i
  real(kind=8),intent(in) :: x,y,z
  integer,intent(in) :: p,q
  character(LEN=1),intent(in) :: a,b

  type(primitive_gaussian) :: pg
  complex(kind=8) :: c_p(4,NMAX_PG),c_q(4,NMAX_PG)

  if(i.lt.1 .or. i.gt.3) then
     write(*,*) "i should be 1-3 in s_mat."
     stop
  end if

  call copy_DiracOutput(p,a,q,b,pg,c_p,c_q)  ! set pg,c_p,c_q
  call calc_s_pq(i,x,y,z,NBS_L,NBS_S,pg,c_p,c_q,s_mat)
  
  return
end function s_mat

!============================================================
function rho_mat(x,y,z,p,a,q,b)
! Electron-positron charge density 
! p,q : labels for molecular orbitals.
! a,b : labels for electron/positron ("+"/"-")
!============================================================
  use DefineTypes
  use DiracOutput
  implicit none
 
  complex(kind=8) :: rho_mat
  real(kind=8),intent(in) :: x,y,z
  integer,intent(in) :: p,q
  character(LEN=1),intent(in) :: a,b

  type(primitive_gaussian) :: pg
  complex(kind=8) :: c_p(4,NMAX_PG),c_q(4,NMAX_PG)

  call copy_DiracOutput(p,a,q,b,pg,c_p,c_q)  ! set pg,c_p,c_q
  call calc_rho_pq(x,y,z,NBS_L,NBS_S,pg,c_p,c_q,rho_mat) 

  return
end function rho_mat


!============================================================
function density_mat(x,y,z,p,a,q,b)
! p,q : labels for molecular orbitals.
! a,b : labels for electron/positron ("+"/"-")
!============================================================
  use DefineTypes
  use DiracOutput
  implicit none
 
  complex(kind=8) :: density_mat
  real(kind=8),intent(in) :: x,y,z
  integer,intent(in) :: p,q
  character(LEN=1),intent(in) :: a,b

  type(primitive_gaussian) :: pg
  complex(kind=8) :: c_p(4,NMAX_PG),c_q(4,NMAX_PG)

  call copy_DiracOutput(p,a,q,b,pg,c_p,c_q)  ! set pg,c_p,c_q
  call calc_density_pq(x,y,z,NBS_L,NBS_S,pg,c_p,c_q,density_mat) 

  return
end function density_mat



!----------------------------------------------------------------------------
!  calc_zeta_pq
!  calc_tau_pq
!  calc_divj_pq
!  calc_j_pq
!  calc_s_pq
!  calc_rho_pq
!  calc_density_pq
!----------------------------------------------------------------------------
!============================================================
subroutine calc_zeta_pq(i,x,y,z,NL,NS,pg,c_p,c_q,zeta_pq_i)
! calculate i-th component of zeta_pq
!============================================================
  use DefineTypes
  use Constants
  implicit none
  
  integer,intent(in) :: i   !zeta^i
  real(kind=8),intent(in) :: x,y,z
  integer,intent(in) :: NL,NS  ! number of primitive gaussian for Large and Small components
  type(primitive_gaussian),intent(in) :: pg
  complex(kind=8),intent(in) :: c_p(4,NMAX_PG),c_q(4,NMAX_PG)

  complex(kind=8),intent(out) :: zeta_pq_i

  complex(kind=8) :: psi_p(4),psi_q(4)
  complex(kind=8) :: dpsi_p(3,4),dpsi_q(3,4)  ! 3:x,y,z
  integer :: m1,m2,m3,m4 ! spinor indice
  integer :: j,k,l

  call calc_psi(x,y,z,NL,NS,pg,c_p,psi_p)
  call calc_psi(x,y,z,NL,NS,pg,c_q,psi_q)
  call calc_dpsi(x,y,z,NL,NS,pg,c_p,dpsi_p)
  call calc_dpsi(x,y,z,NL,NS,pg,c_q,dpsi_q)

  zeta_pq_i = (0.d0,0.d0)
  do m1=1,4
     do m2=1,4
        do m3=1,4
           do m4=1,4
              zeta_pq_i = zeta_pq_i &
                   +conjg(dpsi_p(i,m1))*Gam0(m1,m2)*Gam(i,m2,m3)*Sigma(i,m3,m4)*psi_q(m4) & 
                   +conjg(psi_p(m1))*Gam0(m1,m2)*Gam(i,m2,m3)*Sigma(i,m3,m4)*dpsi_q(i,m4)
           end do
        end do
     end do
  end do
  zeta_pq_i = -(CCC/2.d0)*zeta_pq_i
  return
end subroutine calc_zeta_pq

!============================================================
subroutine calc_tau_pq(k,l,x,y,z,NL,NS,pg,c_p,c_q,tau_pq_kl)
! calculate (k,l) component of tau^kl (stress tensor w/o vector potential)
!============================================================
  use DefineTypes
  use Constants
  implicit none
  
  integer,intent(in) :: k,l  ! tau^kl
  real(kind=8),intent(in) :: x,y,z
  integer,intent(in) :: NL,NS  ! number of primitive gaussian for Large and Small components
  type(primitive_gaussian),intent(in) :: pg
  complex(kind=8),intent(in) :: c_p(4,NMAX_PG),c_q(4,NMAX_PG)

  complex(kind=8),intent(out) :: tau_pq_kl

  complex(kind=8) :: psi_p(4),psi_q(4)
  complex(kind=8) :: dpsi_p(3,4),dpsi_q(3,4)  ! 3:x,y,z
  integer :: m1,m2,m3 ! spinor indice

  call calc_psi(x,y,z,NL,NS,pg,c_p,psi_p)
  call calc_psi(x,y,z,NL,NS,pg,c_q,psi_q)
  call calc_dpsi(x,y,z,NL,NS,pg,c_p,dpsi_p)
  call calc_dpsi(x,y,z,NL,NS,pg,c_q,dpsi_q)

  tau_pq_kl = (0.d0,0.d0)
  do m1=1,4
     do m2=1,4
        do m3=1,4
           tau_pq_kl = tau_pq_kl &
                +conjg(psi_p(m1))*Gam0(m1,m2)*Gam(l,m2,m3)*dpsi_q(k,m3) &
                -conjg(dpsi_p(k,m1))*Gam0(m1,m2)*Gam(l,m2,m3)*psi_q(m3) 
        end do
     end do
  end do

!!$  do m1=1,4
!!$     write(*,'(2es16.6)',advance="no") psi_p(m1)
!!$  end do
!!$  do m1=1,4
!!$     write(*,'(2es16.6)',advance="no") dpsi_p(1,m1)
!!$  end do
!!$  do m1=1,4
!!$     write(*,'(2es16.6)',advance="no") dpsi_p(2,m1)
!!$  end do
!!$  do m1=1,4
!!$     write(*,'(2es16.6)',advance="no") dpsi_p(3,m1)
!!$  end do

  tau_pq_kl = (IU*CCC/2.d0)*tau_pq_kl
  return
end subroutine calc_tau_pq

!============================================================
subroutine calc_divj_pq(x,y,z,NL,NS,pg,c_p,c_q,divj_pq)
! calculate div j
!============================================================
  use DefineTypes
  use Constants
  implicit none
  
  real(kind=8),intent(in) :: x,y,z
  integer,intent(in) :: NL,NS  ! number of primitive gaussian for Large and Small components
  type(primitive_gaussian),intent(in) :: pg
  complex(kind=8),intent(in) :: c_p(4,NMAX_PG),c_q(4,NMAX_PG)

  complex(kind=8),intent(out) :: divj_pq

  complex(kind=8) :: psi_p(4),psi_q(4)
  complex(kind=8) :: dpsi_p(3,4),dpsi_q(3,4)  ! 3:x,y,z
  integer :: i,j,k,l

  call calc_psi(x,y,z,NL,NS,pg,c_p,psi_p)
  call calc_psi(x,y,z,NL,NS,pg,c_q,psi_q)
  call calc_dpsi(x,y,z,NL,NS,pg,c_p,dpsi_p)
  call calc_dpsi(x,y,z,NL,NS,pg,c_q,dpsi_q)

  divj_pq = (0.d0,0.d0)
  do i=1,3
     do j=1,4
        do k=1,4
           do l=1,4
              divj_pq = divj_pq &
                   +conjg(dpsi_p(i,j))*Gam0(j,k)*Gam(i,k,l)*psi_q(l) &
                   +conjg(psi_p(j))*Gam0(j,k)*Gam(i,k,l)*dpsi_q(i,l) 
           end do
        end do
     end do
  end do
  divj_pq = Ze*CCC*divj_pq
  return
end subroutine calc_divj_pq


!============================================================
subroutine calc_j_pq(i,x,y,z,NL,NS,pg,c_p,c_q,j_pq_i)
! calculate i-th component of j_pq
! output : j^i_pq (bilinear for charge current density) (i=1-3)
!============================================================
  use DefineTypes
  use Constants
  implicit none
  
  integer,intent(in) :: i   !j^i
  real(kind=8),intent(in) :: x,y,z
  integer,intent(in) :: NL,NS  ! number of primitive gaussian for Large and Small components
  type(primitive_gaussian),intent(in) :: pg
  complex(kind=8),intent(in) :: c_p(4,NMAX_PG),c_q(4,NMAX_PG)

  complex(kind=8),intent(out) :: j_pq_i

  complex(kind=8) :: psi_p(4),psi_q(4)
  integer :: j,k,l

  call calc_psi(x,y,z,NL,NS,pg,c_p,psi_p)
  call calc_psi(x,y,z,NL,NS,pg,c_q,psi_q)

  j_pq_i = (0.d0,0.d0)
  do j=1,4
     do k=1,4
        do l=1,4
           j_pq_i = j_pq_i +conjg(psi_p(j))*Gam0(j,k)*Gam(i,k,l)*psi_q(l) 
        end do
     end do
  end do
  j_pq_i = Ze*CCC*j_pq_i
  return
end subroutine calc_j_pq

!============================================================
subroutine calc_s_pq(i,x,y,z,NL,NS,pg,c_p,c_q,s_pq_i)
! calculate i-th component of s_pq
!============================================================
  use DefineTypes
  use Constants
  implicit none
  
  integer,intent(in) :: i   !s^i
  real(kind=8),intent(in) :: x,y,z
  integer,intent(in) :: NL,NS  ! number of primitive gaussian for Large and Small components
  type(primitive_gaussian),intent(in) :: pg
  complex(kind=8),intent(in) :: c_p(4,NMAX_PG),c_q(4,NMAX_PG)

  complex(kind=8),intent(out) :: s_pq_i

  complex(kind=8) :: psi_p(4),psi_q(4)
  integer :: j,k

  call calc_psi(x,y,z,NL,NS,pg,c_p,psi_p)
  call calc_psi(x,y,z,NL,NS,pg,c_q,psi_q)

  s_pq_i = (0.d0,0.d0)
  do j=1,4
     do k=1,4
        s_pq_i = s_pq_i +conjg(psi_p(j))*Sigma(i,j,k)*psi_q(k)
     end do
  end do
  s_pq_i = s_pq_i/2.d0
  return
end subroutine calc_s_pq

!============================================================
subroutine calc_rho_pq(x,y,z,NL,NS,pg,c_p,c_q,rho_pq)
! Output : rho_pq
!============================================================
  use DefineTypes
  use Constants
  implicit none
  
  real(kind=8),intent(in) :: x,y,z
  integer,intent(in) :: NL,NS  ! number of primitive gaussian for Large and Small components
  type(primitive_gaussian),intent(in) :: pg
  complex(kind=8),intent(in) :: c_p(4,NMAX_PG),c_q(4,NMAX_PG)

  complex(kind=8),intent(out) :: rho_pq

  complex(kind=8) :: psi_p(4),psi_q(4)
  integer :: j

  call calc_psi(x,y,z,NL,NS,pg,c_p,psi_p)
  call calc_psi(x,y,z,NL,NS,pg,c_q,psi_q)
  
  rho_pq = (0.d0,0.d0)
  do j=1,4
     rho_pq = rho_pq +conjg(psi_p(j))*psi_q(j)
  end do

  rho_pq = Ze * rho_pq

  return
end subroutine calc_rho_pq


!============================================================
subroutine calc_density_pq(x,y,z,NL,NS,pg,c_p,c_q,density_pq)
! output : n_pq
!============================================================
  use DefineTypes
  implicit none
  
  real(kind=8),intent(in) :: x,y,z
  integer,intent(in) :: NL,NS  ! number of primitive gaussian for Large and Small components
  type(primitive_gaussian),intent(in) :: pg
  complex(kind=8),intent(in) :: c_p(4,NMAX_PG),c_q(4,NMAX_PG)

  complex(kind=8),intent(out) :: density_pq

  complex(kind=8) :: psi_p(4),psi_q(4)
  integer :: j

  call calc_psi(x,y,z,NL,NS,pg,c_p,psi_p)
  call calc_psi(x,y,z,NL,NS,pg,c_q,psi_q)
  
  density_pq = (0.d0,0.d0)
  do j=1,4
     density_pq = density_pq +conjg(psi_p(j))*psi_q(j)
  end do
  return
end subroutine calc_density_pq


!----------------------------------------------------------------------------
!  calc_psi
!  calc_dpsi
!----------------------------------------------------------------------------
!============================================================
subroutine calc_psi(x,y,z,NL,NS,pg,c,psi)
! output : psi(4)
!============================================================
  use DefineTypes
  implicit none
  
  real(kind=8),intent(in) :: x,y,z
  integer,intent(in) :: NL,NS  ! number of primitive gaussian for Large and Small components
  type(primitive_gaussian),intent(in) :: pg
  complex(kind=8),intent(in) :: c(4,NMAX_PG)  ! N is defined in DefineTypes
  
  complex(kind=8),intent(out) :: psi(4)

  integer :: nx,ny,nz
  real(kind=8) :: posA(3),alpha
  real(kind=8) :: func_pg  ! primitive gaussian function (normalized)
  integer :: j

  psi(1:4) = (0.d0,0.d0)
  do j=1,NL
     posA(1)=pg%xL(j); posA(2)=pg%yL(j); posA(3)=pg%zL(j); alpha=pg%aL(j); nx=pg%nxL(j); ny=pg%nyL(j); nz=pg%nzL(j)
     psi(1) = psi(1) +c(1,j)*func_pg(x,y,z,posA,alpha,nx,ny,nz)
     psi(2) = psi(2) +c(2,j)*func_pg(x,y,z,posA,alpha,nx,ny,nz)
  end do
  do j=1,NS
     posA(1)=pg%xS(j); posA(2)=pg%yS(j); posA(3)=pg%zS(j); alpha=pg%aS(j); nx=pg%nxS(j); ny=pg%nyS(j); nz=pg%nzS(j)
     psi(3) = psi(3) +c(3,j)*func_pg(x,y,z,posA,alpha,nx,ny,nz)
     psi(4) = psi(4) +c(4,j)*func_pg(x,y,z,posA,alpha,nx,ny,nz)
  end do

end subroutine calc_psi


!============================================================
subroutine calc_dpsi(x,y,z,NL,NS,pg,c,dpsi)
! output : dpsidx,dpsidy,dpsidz
!============================================================
  use DefineTypes
  implicit none
  
  real(kind=8),intent(in) :: x,y,z
  integer,intent(in) :: NL,NS  ! number of primitive gaussian for Large and Small components
  type(primitive_gaussian),intent(in) :: pg
  complex(kind=8),intent(in) :: c(4,NMAX_PG)  ! NMAX_PG is defined in DefineTypes
  
  complex(kind=8),intent(out) :: dpsi(3,4) ! d_k psi

  integer :: nx,ny,nz
  real(kind=8) :: posA(3),alpha
  real(kind=8) :: func_dpg  ! derivative of primitive gaussian function (normalized)
  integer :: i,j

  dpsi(1:3,1:4) = (0.d0,0.d0)
  do i=1,3
     do j=1,NL
        posA(1)=pg%xL(j); posA(2)=pg%yL(j); posA(3)=pg%zL(j); alpha=pg%aL(j); nx=pg%nxL(j); ny=pg%nyL(j); nz=pg%nzL(j)
        dpsi(i,1) = dpsi(i,1) +c(1,j)*func_dpg(i,x,y,z,posA,alpha,nx,ny,nz)
        dpsi(i,2) = dpsi(i,2) +c(2,j)*func_dpg(i,x,y,z,posA,alpha,nx,ny,nz)
     end do
     do j=1,NS
        posA(1)=pg%xS(j); posA(2)=pg%yS(j); posA(3)=pg%zS(j); alpha=pg%aS(j); nx=pg%nxS(j); ny=pg%nyS(j); nz=pg%nzS(j)
        dpsi(i,3) = dpsi(i,3) +c(3,j)*func_dpg(i,x,y,z,posA,alpha,nx,ny,nz)
        dpsi(i,4) = dpsi(i,4) +c(4,j)*func_dpg(i,x,y,z,posA,alpha,nx,ny,nz)
     end do
  end do
  return
end subroutine calc_dpsi


!-------------------------------------------------------------------------
! copy_DiracOutput
!-------------------------------------------------------------------------
!=======================================================
subroutine copy_DiracOutput(p,a,q,b,pg,c_p,c_q)
! copy p.g. values and appropriate coefficients according to the value of p,q and a,b
! 2010.12.2 modified to include Kramers pair
!    If negative values are in p or q, they denote Kramers pair of -p or -q MO.
!    We do similar transformation for positron coefficients (just formally).
!=======================================================
  use DefineTypes
  use DiracOutput
  implicit none

  integer,intent(in) :: p,q
  character(LEN=1),intent(in) :: a,b
  type(primitive_gaussian),intent(out) :: pg
  complex(kind=8),intent(out) :: c_p(4,NMAX_PG),c_q(4,NMAX_PG)

  complex(kind=8) :: cLa_p(NBS_L),cLb_p(NBS_L),cSa_p(NBS_S),cSb_p(NBS_S)
  complex(kind=8) :: cLa_q(NBS_L),cLb_q(NBS_L),cSa_q(NBS_S),cSb_q(NBS_S)
  complex(kind=8) :: dLa_p(NBS_L),dLb_p(NBS_L),dSa_p(NBS_S),dSb_p(NBS_S)
  complex(kind=8) :: dLa_q(NBS_L),dLb_q(NBS_L),dSa_q(NBS_S),dSb_q(NBS_S)
  integer :: j

  do j=1,NBS_L
     pg%aL(j) = aa_L(j); pg%xL(j) = xx_L(j); pg%yL(j) = yy_L(j); pg%zL(j) = zz_L(j)
     pg%nxL(j) = nx_L(j); pg%nyL(j) = ny_L(j); pg%nzL(j) = nz_L(j)
  end do
  do j=1,NBS_S
     pg%aS(j) = aa_S(j); pg%xS(j) = xx_S(j); pg%yS(j) = yy_S(j); pg%zS(j) = zz_S(j)
     pg%nxS(j) = nx_S(j); pg%nyS(j) = ny_S(j); pg%nzS(j) = nz_S(j)
  end do

  ! c_La,c_Lb,c_Sa,c_Sb for electron, d_La,d_Lb,d_Sa,d_Sb for positron
  ! La->1, Lb->2, Sa->3, Sb->4, no matter components are really large or small

  ! Kramers pair (MO with bar) is obtained by transforming
  ! bar(c_La) = -(c_Lb)*
  ! bar(c_Lb) =  (c_La)*
  ! bar(c_Sa) = -(c_Sb)*
  ! bar(c_Sb) =  (c_Sa)*

  if(p.gt.0) then
     do j=1,NBS_L
        cLa_p(j) = c_La(j,p); cLb_p(j) = c_Lb(j,p)
        dLa_p(j) = d_La(j,p); dLb_p(j) = d_Lb(j,p)
     end do
     do j=1,NBS_S
        cSa_p(j) = c_Sa(j,p); cSb_p(j) = c_Sb(j,p)
        dSa_p(j) = d_Sa(j,p); dSb_p(j) = d_Sb(j,p)
     end do
  elseif(p.lt.0) then ! Kramers pair
     do j=1,NBS_L
        cLa_p(j) = -conjg(c_Lb(j,-p)); cLb_p(j) = conjg(c_La(j,-p))
        dLa_p(j) = -conjg(d_Lb(j,-p)); dLb_p(j) = conjg(d_La(j,-p))
     end do
     do j=1,NBS_S
        cSa_p(j) = -conjg(c_Sb(j,-p)); cSb_p(j) = conjg(c_Sa(j,-p))
        dSa_p(j) = -conjg(d_Sb(j,-p)); dSb_p(j) = conjg(d_Sa(j,-p))
     end do
  end if

  if(q.gt.0) then
     do j=1,NBS_L
        cLa_q(j) = c_La(j,q); cLb_q(j) = c_Lb(j,q)
        dLa_q(j) = d_La(j,q); dLb_q(j) = d_Lb(j,q)
     end do
     do j=1,NBS_S
        cSa_q(j) = c_Sa(j,q); cSb_q(j) = c_Sb(j,q)
        dSa_q(j) = d_Sa(j,q); dSb_q(j) = d_Sb(j,q)
     end do
  elseif(q.lt.0) then ! Kramers qair
     do j=1,NBS_L
        cLa_q(j) = -conjg(c_Lb(j,-q)); cLb_q(j) = conjg(c_La(j,-q))
        dLa_q(j) = -conjg(d_Lb(j,-q)); dLb_q(j) = conjg(d_La(j,-q))
     end do
     do j=1,NBS_S
        cSa_q(j) = -conjg(c_Sb(j,-q)); cSb_q(j) = conjg(c_Sa(j,-q))
        dSa_q(j) = -conjg(d_Sb(j,-q)); dSb_q(j) = conjg(d_Sa(j,-q))
     end do
  end if


  c_p=(0.d0,0.d0);  c_q=(0.d0,0.d0)

  if(a.eq."+") then
     do j=1,NBS_L
        c_p(1,j) = cLa_p(j); c_p(2,j) = cLb_p(j)
     end do
     do j=1,NBS_S
        c_p(3,j) = cSa_p(j); c_p(4,j) = cSb_p(j)
     end do
  elseif(a.eq."-") then
     do j=1,NBS_L
        c_p(1,j) = dLa_p(j); c_p(2,j) = dLb_p(j)
     end do
     do j=1,NBS_S
        c_p(3,j) = dSa_p(j); c_p(4,j) = dSb_p(j)
     end do
  else
     write(*,*) "Use + or - in function ..._mat."
     stop
  end if

  if(b.eq."+") then
     do j=1,NBS_L
        c_q(1,j) = cLa_q(j); c_q(2,j) = cLb_q(j)
     end do
     do j=1,NBS_S
        c_q(3,j) = cSa_q(j); c_q(4,j) = cSb_q(j)
     end do
  elseif(b.eq."-") then
     do j=1,NBS_L
        c_q(1,j) = dLa_q(j); c_q(2,j) = dLb_q(j)
     end do
     do j=1,NBS_S
        c_q(3,j) = dSa_q(j); c_q(4,j) = dSb_q(j)
     end do
  else
     write(*,*) "Use + or - in function ..._mat."
     stop
  end if

  return
end subroutine copy_DiracOutput


!=======================================================
subroutine copy_DiracOutput_cp(p,a,c_p)
! copy coefficients according to the value of p and a
! 2010.12.2 modified to include Kramers pair
!    If negative values are in p or q, they denote Kramers pair of -p or -q MO.
!    We do similar transformation for positron coefficients (just formally).
! 2012.5.19 modified from copy_DiracOutput
!    a is specified by 1(+,electron) or 2(-,positron)
!    to be used in modified twoele integration routines.
!=======================================================
  use DefineTypes
  use DiracOutput
  implicit none

  integer,intent(in) :: p,a
  complex(kind=8),intent(out) :: c_p(4,NMAX_PG)

  complex(kind=8) :: cLa_p(NBS_L),cLb_p(NBS_L),cSa_p(NBS_S),cSb_p(NBS_S)
  complex(kind=8) :: dLa_p(NBS_L),dLb_p(NBS_L),dSa_p(NBS_S),dSb_p(NBS_S)
  integer :: j

  ! c_La,c_Lb,c_Sa,c_Sb for electron, d_La,d_Lb,d_Sa,d_Sb for positron
  ! La->1, Lb->2, Sa->3, Sb->4, no matter components are really large or small

  ! Kramers pair (MO with bar) is obtained by transforming
  ! bar(c_La) = -(c_Lb)*
  ! bar(c_Lb) =  (c_La)*
  ! bar(c_Sa) = -(c_Sb)*
  ! bar(c_Sb) =  (c_Sa)*

  if(p.gt.0) then
     do j=1,NBS_L
        cLa_p(j) = c_La(j,p); cLb_p(j) = c_Lb(j,p)
        dLa_p(j) = d_La(j,p); dLb_p(j) = d_Lb(j,p)
     end do
     do j=1,NBS_S
        cSa_p(j) = c_Sa(j,p); cSb_p(j) = c_Sb(j,p)
        dSa_p(j) = d_Sa(j,p); dSb_p(j) = d_Sb(j,p)
     end do
  elseif(p.lt.0) then ! Kramers pair
     do j=1,NBS_L
        cLa_p(j) = -conjg(c_Lb(j,-p)); cLb_p(j) = conjg(c_La(j,-p))
        dLa_p(j) = -conjg(d_Lb(j,-p)); dLb_p(j) = conjg(d_La(j,-p))
     end do
     do j=1,NBS_S
        cSa_p(j) = -conjg(c_Sb(j,-p)); cSb_p(j) = conjg(c_Sa(j,-p))
        dSa_p(j) = -conjg(d_Sb(j,-p)); dSb_p(j) = conjg(d_Sa(j,-p))
     end do
  elseif(p.eq.0) then
     c_p = (0.d0,0.d0)
     !stop
     return
  end if

  c_p=(0.d0,0.d0)

  if(a.eq.1) then
     do j=1,NBS_L
        c_p(1,j) = cLa_p(j); c_p(2,j) = cLb_p(j)
     end do
     do j=1,NBS_S
        c_p(3,j) = cSa_p(j); c_p(4,j) = cSb_p(j)
     end do
  elseif(a.eq.2) then
     do j=1,NBS_L
        c_p(1,j) = dLa_p(j); c_p(2,j) = dLb_p(j)
     end do
     do j=1,NBS_S
        c_p(3,j) = dSa_p(j); c_p(4,j) = dSb_p(j)
     end do
  else
     write(*,*) "Use + or - in function ..._mat."
     stop
  end if

  return
end subroutine copy_DiracOutput_cp

!=======================================================
subroutine copy_DiracOutput_pg(pg)
! copy p.g. values
! 2012.5.19 modified from copy_DiracOutput
!=======================================================
  use DefineTypes
  use DiracOutput
  implicit none

  type(primitive_gaussian),intent(out) :: pg
  integer :: j

  do j=1,NBS_L
     pg%aL(j) = aa_L(j); pg%xL(j) = xx_L(j); pg%yL(j) = yy_L(j); pg%zL(j) = zz_L(j)
     pg%nxL(j) = nx_L(j); pg%nyL(j) = ny_L(j); pg%nzL(j) = nz_L(j)
  end do
  do j=1,NBS_S
     pg%aS(j) = aa_S(j); pg%xS(j) = xx_S(j); pg%yS(j) = yy_S(j); pg%zS(j) = zz_S(j)
     pg%nxS(j) = nx_S(j); pg%nyS(j) = ny_S(j); pg%nzS(j) = nz_S(j)
  end do

  return
end subroutine copy_DiracOutput_pg
