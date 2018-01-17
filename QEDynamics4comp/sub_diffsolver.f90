!================================================================================
! subroutines for solvers of differential equations of rigged QED
!
! 2011.11.29
!  -Euler method
!  -2nd order Runge-Kutta
!
!  -moved calc_dcalEdt and calc_dcalCdt from simQED.f90
!
! 11.30
!  -4th order Runge-Kutta
!
!================================================================================



!================================================================================
subroutine diffsol_euler(Time,calE,calC,TM_Qmat,Tnuc_mat, &
     & twoele_Qmat,nucele_Qmat,twonuc_mat,calFj0_Qmat,alpha, &
     & dcalE,dcalC)
!
! compute next step by Euler method.
! y_n+1 = y_n + h*f(t_n,y_n)
!
! y : calE and calC
! f : dcalEdt and dcalCdt
! h : DeltaT
!================================================================================
  Use DiracOutput
  use Constants
  use NucBasis

  implicit none

  real(kind=8),intent(in) :: time 
  complex(kind=8),intent(in) :: calE(4*NBS,4*NBS)  ! calE_PQ
  complex(kind=8),intent(in) :: calC(NBS_N,NBS_N)  ! calC_ij
  complex(kind=8),intent(in) :: TM_Qmat(4*NBS,4*NBS)  ! T_PQ + M_PQ
  complex(kind=8),intent(in) :: Tnuc_mat(NBS_N,NBS_N)  ! T_Nij
  complex(kind=8),intent(in) :: twoele_Qmat(4*NBS,4*NBS,4*NBS,4*NBS)  ! (PQ|RS)
  complex(kind=8),intent(in) :: nucele_Qmat(NBS_N,NBS_N,4*NBS,4*NBS)  ! (ij|PQ)
  complex(kind=8),intent(in) :: twonuc_mat(NBS_N,NBS_N,NBS_N,NBS_N)  ! (ij|kl)
  complex(kind=8),intent(in) :: calFj0_Qmat(4*NBS,4*NBS,Nph) ! calF_PQj(t=0)
  complex(kind=8),intent(in) :: alpha(Nph)   ! <coherent|hat(a)|coherent>
  
  complex(kind=8),intent(out) :: dcalE(4*NBS,4*NBS)  ! d(calE_PQ)
  complex(kind=8),intent(out) :: dcalC(NBS_N,NBS_N)  ! d(calC_ij)

  integer :: pp,qq,i,j
  complex(kind=8) :: dcalEdt,dcalCdt

  ! electron
  do pp=1,4*NBS
     do qq=1,4*NBS
        ! calculate derivative
        call calc_dcalEdt(pp,qq,time,calE,calC,TM_Qmat,twoele_Qmat,nucele_Qmat,calFj0_Qmat,alpha,dcalEdt)
        dcalE(pp,qq) = DeltaT*dcalEdt
     end do
  end do

  ! nucleus
  do i=1,NBS_N
     do j=1,NBS_N
        ! calculate derivative
        call calc_dcalCdt(i,j,time,calE,calC,Tnuc_mat,nucele_Qmat,twonuc_mat,alpha,dcalCdt)
        dcalC(i,j) = DeltaT*dcalCdt
     end do
  end do

  return
end subroutine diffsol_euler


!================================================================================
subroutine diffsol_rk2(Time,calE,calC,TM_Qmat,Tnuc_mat, &
     & twoele_Qmat,nucele_Qmat,twonuc_mat,calFj0_Qmat,alpha, &
     & dcalE,dcalC)
!
! compute next step by 2nd order Runge-Kutta.
! y_n+1 = y_n + k2
! k2 = h*f(t_n +h/2, y_n +k1/2)
! k1 = h*f(t_n,y_n)
!
! y : calE and calC
! f : dcalEdt and dcalCdt
! h : DeltaT
!================================================================================
  Use DiracOutput
  use Constants
  use NucBasis

  implicit none

  real(kind=8),intent(in) :: time 
  complex(kind=8),intent(in) :: calE(4*NBS,4*NBS)  ! calE_PQ
  complex(kind=8),intent(in) :: calC(NBS_N,NBS_N)  ! calC_ij
  complex(kind=8),intent(in) :: TM_Qmat(4*NBS,4*NBS)  ! T_PQ + M_PQ
  complex(kind=8),intent(in) :: Tnuc_mat(NBS_N,NBS_N)  ! T_Nij
  complex(kind=8),intent(in) :: twoele_Qmat(4*NBS,4*NBS,4*NBS,4*NBS)  ! (PQ|RS)
  complex(kind=8),intent(in) :: nucele_Qmat(NBS_N,NBS_N,4*NBS,4*NBS)  ! (ij|PQ)
  complex(kind=8),intent(in) :: twonuc_mat(NBS_N,NBS_N,NBS_N,NBS_N)  ! (ij|kl)
  complex(kind=8),intent(in) :: calFj0_Qmat(4*NBS,4*NBS,Nph) ! calF_PQj(t=0)
  complex(kind=8),intent(in) :: alpha(Nph)   ! <coherent|hat(a)|coherent>
  
  complex(kind=8),intent(out) :: dcalE(4*NBS,4*NBS)  ! d(calE_PQ)
  complex(kind=8),intent(out) :: dcalC(NBS_N,NBS_N)  ! d(calC_ij)

  integer :: pp,qq,i,j
  complex(kind=8) :: dcalEdt,dcalCdt

  real(kind=8) :: time_2 ! t_n +h/2
  complex(kind=8) :: calE_2(4*NBS,4*NBS)  ! y_n +k1/2 (calE)
  complex(kind=8) :: calC_2(NBS_N,NBS_N)  ! y_n +k1/2 (calC)

  !--------------------
  ! compute k1/2
  !--------------------
  do pp=1,4*NBS
     do qq=1,4*NBS
        call calc_dcalEdt(pp,qq,time,calE,calC,TM_Qmat,twoele_Qmat,nucele_Qmat,calFj0_Qmat,alpha,dcalEdt)
        dcalE(pp,qq) = DeltaT*dcalEdt/2.d0
     end do
  end do
  do i=1,NBS_N
     do j=1,NBS_N
        call calc_dcalCdt(i,j,time,calE,calC,Tnuc_mat,nucele_Qmat,twonuc_mat,alpha,dcalCdt)
        dcalC(i,j) = DeltaT*dcalCdt/2.d0
     end do
  end do

  !--------------------
  ! compute yn +k1/2
  !--------------------
  do pp=1,4*NBS
     do qq=1,4*NBS
        calE_2(pp,qq) = calE(pp,qq) + dcalE(pp,qq)
     end do
  end do
  do i=1,NBS_N
     do j=1,NBS_N
        calC_2(i,j) = calC(i,j) + dcalC(i,j)
     end do
  end do
  
  !--------------------
  ! compute k2
  !--------------------
  time_2 = time + DeltaT/2.d0

  do pp=1,4*NBS
     do qq=1,4*NBS
        call calc_dcalEdt(pp,qq,time_2,calE_2,calC_2,TM_Qmat,twoele_Qmat,nucele_Qmat,calFj0_Qmat,alpha,dcalEdt)
        dcalE(pp,qq) = DeltaT*dcalEdt
     end do
  end do
  do i=1,NBS_N
     do j=1,NBS_N
        call calc_dcalCdt(i,j,time_2,calE_2,calC_2,Tnuc_mat,nucele_Qmat,twonuc_mat,alpha,dcalCdt)
        dcalC(i,j) = DeltaT*dcalCdt
     end do
  end do
  
  return
end subroutine diffsol_rk2

!================================================================================
subroutine diffsol_rk4(Time,calE,calC,TM_Qmat,Tnuc_mat, &
     & twoele_Qmat,nucele_Qmat,twonuc_mat,calFj0_Qmat,alpha, &
     & dcalE,dcalC)
!
! compute next step by 4th order Runge-Kutta.
! y_n+1 = y_n + (1/6)*(k1 +2*k2 +2*k3 +k4)
! k1 = h*f(t_n     , y_n)
! k2 = h*f(t_n +h/2, y_n +k1/2)
! k3 = h*f(t_n +h/2, y_n +k2/2)
! k4 = h*f(t_n +h  , y_n +k3)
!
! y : calE and calC
! f : dcalEdt and dcalCdt
! h : DeltaT
!================================================================================
  Use DiracOutput
  use Constants
  use NucBasis

  implicit none

  real(kind=8),intent(in) :: time 
  complex(kind=8),intent(in) :: calE(4*NBS,4*NBS)  ! calE_PQ
  complex(kind=8),intent(in) :: calC(NBS_N,NBS_N)  ! calC_ij
  complex(kind=8),intent(in) :: TM_Qmat(4*NBS,4*NBS)  ! T_PQ + M_PQ
  complex(kind=8),intent(in) :: Tnuc_mat(NBS_N,NBS_N)  ! T_Nij
  complex(kind=8),intent(in) :: twoele_Qmat(4*NBS,4*NBS,4*NBS,4*NBS)  ! (PQ|RS)
  complex(kind=8),intent(in) :: nucele_Qmat(NBS_N,NBS_N,4*NBS,4*NBS)  ! (ij|PQ)
  complex(kind=8),intent(in) :: twonuc_mat(NBS_N,NBS_N,NBS_N,NBS_N)  ! (ij|kl)
  complex(kind=8),intent(in) :: calFj0_Qmat(4*NBS,4*NBS,Nph) ! calF_PQj(t=0)
  complex(kind=8),intent(in) :: alpha(Nph)   ! <coherent|hat(a)|coherent>
  
  complex(kind=8),intent(out) :: dcalE(4*NBS,4*NBS)  ! d(calE_PQ)
  complex(kind=8),intent(out) :: dcalC(NBS_N,NBS_N)  ! d(calC_ij)

  integer :: pp,qq,i,j
  complex(kind=8) :: dcalEdt,dcalCdt

  real(kind=8) :: time_2 ! t_n +h/2
  real(kind=8) :: time_3 ! t_n +h
  complex(kind=8) :: calE_tmp(4*NBS,4*NBS)
  complex(kind=8) :: calC_tmp(NBS_N,NBS_N)
  complex(kind=8) :: k1_calE(4*NBS,4*NBS),k2_calE(4*NBS,4*NBS),k3_calE(4*NBS,4*NBS),k4_calE(4*NBS,4*NBS)
  complex(kind=8) :: k1_calC(NBS_N,NBS_N),k2_calC(NBS_N,NBS_N),k3_calC(NBS_N,NBS_N),k4_calC(NBS_N,NBS_N)

  time_2 = time + DeltaT/2.d0 ! for k2,k3
  time_3 = time + DeltaT ! for k4

  !-----------------------------
  ! compute k1 = h*f(t_n, y_n)
  !-----------------------------
  do pp=1,4*NBS
     do qq=1,4*NBS
        call calc_dcalEdt(pp,qq,time,calE,calC,TM_Qmat,twoele_Qmat,nucele_Qmat,calFj0_Qmat,alpha,dcalEdt)
        k1_calE(pp,qq) = DeltaT*dcalEdt
     end do
  end do
  do i=1,NBS_N
     do j=1,NBS_N
        call calc_dcalCdt(i,j,time,calE,calC,Tnuc_mat,nucele_Qmat,twonuc_mat,alpha,dcalCdt)
        k1_calC(i,j) = DeltaT*dcalCdt
     end do
  end do

  !--------------------
  ! compute y_n +k1/2
  !--------------------
  do pp=1,4*NBS
     do qq=1,4*NBS
        calE_tmp(pp,qq) = calE(pp,qq) + k1_calE(pp,qq)/2.d0
     end do
  end do
  do i=1,NBS_N
     do j=1,NBS_N
        calC_tmp(i,j) = calC(i,j) + k1_calC(i,j)/2.d0
     end do
  end do
  
  !------------------------------------------------
  ! compute k2　= h*f(t_n +h/2, y_n +k1/2)
  !------------------------------------------------
  do pp=1,4*NBS
     do qq=1,4*NBS
        call calc_dcalEdt(pp,qq,time_2,calE_tmp,calC_tmp,TM_Qmat,twoele_Qmat,nucele_Qmat,calFj0_Qmat,alpha,dcalEdt)
        k2_calE(pp,qq) = DeltaT*dcalEdt
     end do
  end do
  do i=1,NBS_N
     do j=1,NBS_N
        call calc_dcalCdt(i,j,time_2,calE_tmp,calC_tmp,Tnuc_mat,nucele_Qmat,twonuc_mat,alpha,dcalCdt)
        k2_calC(i,j) = DeltaT*dcalCdt
     end do
  end do
  
  !--------------------
  ! compute y_n +k2/2
  !--------------------
  do pp=1,4*NBS
     do qq=1,4*NBS
        calE_tmp(pp,qq) = calE(pp,qq) + k2_calE(pp,qq)/2.d0
     end do
  end do
  do i=1,NBS_N
     do j=1,NBS_N
        calC_tmp(i,j) = calC(i,j) + k2_calC(i,j)/2.d0
     end do
  end do

  !------------------------------------------------
  ! compute k3　= h*f(t_n +h/2, y_n +k2/2)
  !------------------------------------------------
  do pp=1,4*NBS
     do qq=1,4*NBS
        call calc_dcalEdt(pp,qq,time_2,calE_tmp,calC_tmp,TM_Qmat,twoele_Qmat,nucele_Qmat,calFj0_Qmat,alpha,dcalEdt)
        k3_calE(pp,qq) = DeltaT*dcalEdt
     end do
  end do
  do i=1,NBS_N
     do j=1,NBS_N
        call calc_dcalCdt(i,j,time_2,calE_tmp,calC_tmp,Tnuc_mat,nucele_Qmat,twonuc_mat,alpha,dcalCdt)
        k3_calC(i,j) = DeltaT*dcalCdt
     end do
  end do

  !--------------------
  ! compute y_n +k3
  !--------------------
  do pp=1,4*NBS
     do qq=1,4*NBS
        calE_tmp(pp,qq) = calE(pp,qq) + k3_calE(pp,qq)
     end do
  end do
  do i=1,NBS_N
     do j=1,NBS_N
        calC_tmp(i,j) = calC(i,j) + k3_calC(i,j)
     end do
  end do

  !------------------------------------------------
  ! compute k4　= h*f(t_n +h, y_n +k3)
  !------------------------------------------------
  do pp=1,4*NBS
     do qq=1,4*NBS
        call calc_dcalEdt(pp,qq,time_3,calE_tmp,calC_tmp,TM_Qmat,twoele_Qmat,nucele_Qmat,calFj0_Qmat,alpha,dcalEdt)
        k4_calE(pp,qq) = DeltaT*dcalEdt
     end do
  end do
  do i=1,NBS_N
     do j=1,NBS_N
        call calc_dcalCdt(i,j,time_3,calE_tmp,calC_tmp,Tnuc_mat,nucele_Qmat,twonuc_mat,alpha,dcalCdt)
        k4_calC(i,j) = DeltaT*dcalCdt
     end do
  end do

  !------------------------------------------------
  ! compute  (1/6)*(k1 +2*k2 +2*k3 +k4)
  !------------------------------------------------
  do pp=1,4*NBS
     do qq=1,4*NBS
        dcalE(pp,qq) = (k1_calE(pp,qq) +2.d0*k2_calE(pp,qq) +2.d0*k3_calE(pp,qq) +k4_calE(pp,qq))/6.d0
     end do
  end do
  do i=1,NBS_N
     do j=1,NBS_N
        dcalC(i,j) = (k1_calC(i,j) +2.d0*k4_calC(i,j) +2.d0*k4_calC(i,j) +k4_calC(i,j))/6.d0
     end do
  end do

  return
end subroutine diffsol_rk4


!========================================================================================================
subroutine calc_dcalEdt(pp,qq,time,calE,calC,TM_Qmat,twoele_Qmat,nucele_Qmat,calFj0_Qmat,alpha,dcalEdt)
!========================================================================================================
  Use DiracOutput
  use Constants
  use NucBasis
  implicit none

  integer,intent(in) :: pp,qq ! compute (pp,qq) component of d(calE_PQ)/dt
!  integer,intent(in) :: it  ! it-th time step
  real(kind=8),intent(in) :: time 
  complex(kind=8),intent(in) :: calE(4*NBS,4*NBS)  ! calE_PQ
  complex(kind=8),intent(in) :: calC(NBS_N,NBS_N)  ! calC_ij
  complex(kind=8),intent(in) :: TM_Qmat(4*NBS,4*NBS)  ! T_PQ + M_PQ
  complex(kind=8),intent(in) :: twoele_Qmat(4*NBS,4*NBS,4*NBS,4*NBS)  ! (PQ|RS)
  complex(kind=8),intent(in) :: calFj0_Qmat(4*NBS,4*NBS,Nph) ! calF_PQj(t=0)
  complex(kind=8),intent(in) :: alpha(Nph)   ! <coherent|hat(a)|coherent>
  complex(kind=8),intent(in) :: nucele_Qmat(NBS_N,NBS_N,4*NBS,4*NBS)  ! (ij|PQ)

  complex(kind=8),intent(out) :: dcalEdt  ! d(calE_PQ)/dt

  real(kind=8) :: dpj ! tilde(Delta p)_j
  real(kind=8) :: p0j ! p_j
  real(kind=8) :: t ! time = it*DeltaT

  integer :: i,j,k,l
  integer :: rr,ss,nn,mm
  complex(kind=8) :: sum_te,sum_ne,sum,sum1,sum2
  complex(kind=8) :: dcalF

  complex(kind=8) :: I2_Qmat(4*NBS,4*NBS), I4_Qmat(4*NBS,4*NBS) 
  complex(kind=8) :: I_Qmat(4*NBS,4*NBS)  ! I_PQ = (I_1+I_2+I_3+I_4)_PQ = (TM +I2 +I4)_PQ


  I2_Qmat  = (0.d0,0.d0)
  I4_Qmat  = (0.d0,0.d0)

  ! it starts from 0 -> t starts from 0
!  t = DeltaT*it  ! DeltaT is defined in GammaMatrix and set in the main routine.
  t = time
  

  ! calculation of I4
  do nn=1,4*NBS
     do mm=1,4*NBS

        sum_te = (0.d0,0.d0)
        do rr=1,4*NBS
           do ss=1,4*NBS
              sum_te = sum_te + twoele_Qmat(nn,mm,rr,ss)*calE(rr,ss)
           end do
        end do
        sum_ne = (0.d0,0.d0)
        do i=1,NBS_N
           do j=1,NBS_N
              sum_ne = sum_ne + nucele_Qmat(i,j,nn,mm)*calC(i,j)  ! (PQ|ij)=(ij|PQ)
           end do
        end do
        I4_Qmat(nn,mm) = sum_te + sum_ne

     end do
  end do

  I_Qmat(:,:) = TM_Qmat(:,:) + I2_Qmat(:,:) + I4_Qmat(:,:)

  sum = (0.d0,0.d0)
  sum1 = (0.d0,0.d0)
  sum2 = (0.d0,0.d0)
  do rr=1,4*NBS
     sum = sum -I_Qmat(rr,pp)*calE(rr,qq) +I_Qmat(qq,rr)*calE(pp,rr)
!     sum1 = sum1 -I_Qmat(rr,pp)*calE(rr,qq)
!     sum2 = sum2 +I_Qmat(qq,rr)*calE(pp,rr)
  end do
!  write(10,"(6es16.6)") sum,sum1,sum2
  
  
!!$  ! Arad part
!!$  do j=1,Nph
!!$     call calc_DeltaPj_and_P0j(j,dpj,p0j)
!!$     sum_dcalF=(0.d0,0.d0)
!!$     do rr=1,4*NBS
!!$        dcalF = conjg(calFj0_Qmat(nn,rr,j))*exp(+IU*CCC*p0j*t)*calE(rr,mm)*conjg(alpha(j)) &
!!$             & +calFj0_Qmat(rr,nn,j)       *exp(-IU*CCC*p0j*t)*calE(rr,mm)*alpha(j) &
!!$             & -calFj0_Qmat(mm,rr,j)       *exp(-IU*CCC*p0j*t)*calE(nn,rr)*alpha(j) &
!!$             & -conjg(calFj0_Qmat(rr,mm,j))*exp(+IU*CCC*p0j*t)*calE(nn,rr)*conjg(alpha(j)) 
!!$        sum_dcalF = sum_dcalF +dcalF
!!$     end do
!!$!     write(*,"(3i6,2es16.6,4es16.6)") nn,mm,j,dpj,p0j,calFj0_Qmat(nn,rr,j),sum_calF
!!$!     write(*,"(3i6,2es16.6)") nn,mm,j,sum_calF
!!$     sum_calF = sum_calF -dpj*sum_dcalF
!!$!     write(*,"(3i6,4es16.6)") nn,mm,j,sum_dcalF,sum_calF
!!$  end do
!!$!  write(*,"(2i6,2es16.6)") nn,mm,sum_calF
  

  dcalEdt = -IU*sum
  
  return
end subroutine calc_dcalEdt


!========================================================================================================
subroutine calc_dcalCdt(i,j,time,calE,calC,Tnuc_mat,nucele_Qmat,twonuc_mat,alpha,dcalCdt)
!========================================================================================================
  Use DiracOutput
  use Constants
  use NucBasis
  implicit none

  integer,intent(in) :: i,j ! compute (i,j) component of d(calC_ij)/dt
!  integer,intent(in) :: it  ! it-th time step
  real(kind=8),intent(in) :: time 
  complex(kind=8),intent(in) :: calE(4*NBS,4*NBS)  ! calE_PQ
  complex(kind=8),intent(in) :: calC(NBS_N,NBS_N)  ! calC_ij
  complex(kind=8),intent(in) :: Tnuc_mat(NBS_N,NBS_N)  ! T_Nij
  complex(kind=8),intent(in) :: alpha(Nph)   ! <coherent|hat(a)|coherent>
  complex(kind=8),intent(in) :: nucele_Qmat(NBS_N,NBS_N,4*NBS,4*NBS)  ! (ij|PQ)
  complex(kind=8),intent(in) :: twonuc_mat(NBS_N,NBS_N,NBS_N,NBS_N)  ! (ij|kl)

  complex(kind=8),intent(out) :: dcalCdt  ! d(calC_ij)/dt

  real(kind=8) :: dpj ! tilde(Delta p)_j
  real(kind=8) :: p0j ! p_j
  real(kind=8) :: t ! time = it*DeltaT

  integer :: k,l
  integer :: ii,jj,rr,ss
  complex(kind=8) :: sum_tn,sum_ne,sum
  complex(kind=8) :: dcalF

  complex(kind=8) :: I4_mat(NBS_N,NBS_N) 
  complex(kind=8) :: I_mat(NBS_N,NBS_N)  ! I_ij = (I_1+I_2+I_3+I_4)_ij = (T +I4)_ij


  I4_mat  = (0.d0,0.d0)

  ! it starts from 0 -> t starts from 0
!  t = DeltaT*it  ! DeltaT is defined in GammaMatrix and set in the main routine.
  t = time
  
  ! calculation of I4
  do ii=1,NBS_N
     do jj=1,NBS_N

        sum_ne = (0.d0,0.d0)
        do rr=1,4*NBS
           do ss=1,4*NBS
              sum_ne = sum_ne + nucele_Qmat(ii,jj,rr,ss)*calE(rr,ss)
           end do
        end do
        sum_tn = (0.d0,0.d0)
        do k=1,NBS_N
           do l=1,NBS_N
              sum_tn = sum_tn + twonuc_mat(ii,jj,k,l)*calC(k,l)  
           end do
        end do
        I4_mat(ii,jj) = sum_tn + sum_ne

     end do
  end do

  I_mat(:,:) = Tnuc_mat(:,:) + I4_mat(:,:)

  sum = (0.d0,0.d0)
  do k=1,NBS_N
     sum = sum -I_mat(k,i)*calC(k,j) +I_mat(j,k)*calC(i,k)
  end do
  
  
  dcalCdt = -IU*sum
  
  return
end subroutine calc_dcalCdt

