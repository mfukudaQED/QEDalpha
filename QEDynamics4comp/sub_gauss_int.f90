!========================================================================
!
! functions:
!      func_dunpg(x,y,z,posA,alpha,nx,ny,nz,i)  (added 110110)
!      func_unpg(x,y,z,posA,alpha,nx,ny,nz)  (added 110110)
!      func_dpg(x,y,z,posA,alpha,nx,ny,nz,i)  (added 100209, modified 101112,110110)
!      func_pg(x,y,z,posA,alpha,nx,ny,nz)  (added 100208, modified 100209,110110)
!      norm_pg(alpha,nx,ny,nz)
!      norm_pg1(alpha,n)
!
! subroutines:
!      permutate_cccc(inttype,i,j,k,l,A,B,C,D,coefr,NBS0,NBS,factor)
!      find_int_type(A,B,C,D,inttype,num_term)
!
!      gauss_int_fourier(posA,alphaA,nx,ny,nz,posB,alphaB,nbarx,nbary,nbarz,posC,vecK,fourier) (added 101103)
!      gauss_int_twoele(posA,alphaA,n1,posB,alphaB,nbar1,posC,alphaC,n2,posD,alphaD,nbar2,twoele)
!      gauss_int(posA,alphaA,nx,ny,nz,posB,alphaB,nbarx,nbary,nbarz,posC,overlap,nucatt,ef)
!      gauss_int_moment(posA,alphaA,nx,ny,nz,posB,alphaB,nbarx,nbary,nbarz,posC,moment)
!      gauss_int_KE(posA,alphaA,nx,ny,nz,posB,alphaB,nbarx,nbary,nbarz,ene_kin)
!      gauss_int_relKE(posA,alphaA,nx,ny,nz,posB,alphaB,nbarx,nbary,nbarz,grad)
!      gauss_int_overlap(posA,alphaA,nx,ny,nz,posB,alphaB,nbarx,nbary,nbarz,overlap)
!
!      calc_fourier_integral(i,j,k,alphaP,posP,vecK,int_fourier)  (added 101103)
!      calc_moment_integral(i,j,k,alphaP,vecPC,int_mom)
!      calc_overlap_integral(i,j,k,alphaP,integral)
!      calc_R(vecPC,n_sum,list_R000j,nx_sum,ny_sum,nz_sum,table_R)
!      calc_PC_and_T(posC,posP,alphaP,vecPC,T)
!      calc_R000j(alphaP,T,n,list_R000j)
!      calc_FjT(T,n,list_FjT)
!      calc_d(alphaP,PAx,PBx,n,nbar,d)
!      calc_gaussP(posA,alphaA,posB,alphaB,c_PAB,posP,alphaP,vecPA,vecPB)
!========================================================================

!==================================================================
function func_dpg(i,x,y,z,posA,alpha,nx,ny,nz)
! i(=1,2,3=x,y,z) derivative of normalized gaussian function
! center at (Ax,Ay,Az)
!==================================================================
  implicit none

  real(kind=8) :: func_dpg
  integer,intent(in) :: i
  real(kind=8),intent(in) :: x,y,z,posA(3),alpha
  integer,intent(in) :: nx,ny,nz
  real(kind=8) :: norm_pg,func_dunpg

  func_dpg = norm_pg(alpha,nx,ny,nz) *func_dunpg(i,x,y,z,posA,alpha,nx,ny,nz)

  return
end function func_dpg

!==================================================================
function func_pg(x,y,z,posA,alpha,nx,ny,nz)
! normalized gaussian function
! center at (Ax,Ay,Az)
!==================================================================
  implicit none
  
  real(kind=8) :: func_pg
  real(kind=8),intent(in) :: x,y,z,posA(3),alpha
  integer,intent(in) :: nx,ny,nz
  real(kind=8) :: norm_pg,func_unpg
  
  func_pg = norm_pg(alpha,nx,ny,nz) *func_unpg(x,y,z,posA,alpha,nx,ny,nz)
  
  return
end function func_pg

!==================================================================
function func_dunpg(i,x,y,z,posA,alpha,nx,ny,nz)
! i(=1,2,3=x,y,z) derivative of UNnormalized gaussian function
! center at (Ax,Ay,Az)
!==================================================================
  implicit none

  real(kind=8) :: func_dunpg
  integer,intent(in) :: i
  real(kind=8),intent(in) :: x,y,z,posA(3),alpha
  integer,intent(in) :: nx,ny,nz
  real(kind=8) :: func_unpg

  if(i.eq.1) then  ! d/dx
     func_dunpg = nx*func_unpg(x,y,z,posA,alpha,nx-1,ny,nz) -2.d0*alpha*func_unpg(x,y,z,posA,alpha,nx+1,ny,nz)
  elseif(i.eq.2) then ! d/dy
     func_dunpg = ny*func_unpg(x,y,z,posA,alpha,nx,ny-1,nz) -2.d0*alpha*func_unpg(x,y,z,posA,alpha,nx,ny+1,nz)
  elseif(i.eq.3) then ! d/dz
     func_dunpg = nz*func_unpg(x,y,z,posA,alpha,nx,ny,nz-1) -2.d0*alpha*func_unpg(x,y,z,posA,alpha,nx,ny,nz+1)
  else
     stop "i should be 1 or 2 or 3 in fund_dunpg."
  end if

  return
end function func_dunpg

!==================================================================
function func_unpg(x,y,z,posA,alpha,nx,ny,nz)
! UNnormalized gaussian function
! center at (Ax,Ay,Az)
!==================================================================
  implicit none

  real(kind=8) :: func_unpg
  real(kind=8),intent(in) :: x,y,z,posA(3),alpha
  integer,intent(in) :: nx,ny,nz
  real(kind=8) :: rA2

  if((nx.lt.0).or.(ny.lt.0).or.(nz.lt.0)) then
     func_unpg = 0.d0
  else
     rA2 = (x-posA(1))**2+(y-posA(2))**2+(z-posA(3))**2
     func_unpg = (x-posA(1))**nx *(y-posA(2))**ny *(z-posA(3))**nz *exp(-alpha*rA2)
  end if

  return
end function func_unpg


!==================================================================
function norm_pg(alpha,nx,ny,nz)
!==================================================================
  implicit none
  real(kind=8) :: norm_pg
  real(kind=8),intent(in) :: alpha
  integer,intent(in) :: nx,ny,nz
  real(kind=8) :: norm_pg1

  norm_pg = norm_pg1(alpha,nx)*norm_pg1(alpha,ny)*norm_pg1(alpha,nz)
end function norm_pg

!==================================================================
function norm_pg1(alpha,n)
!==================================================================
  implicit none
  real(kind=8) :: norm_pg1
  real(kind=8),intent(in) :: alpha
  integer,intent(in) :: n

  real(kind=8) :: PI
  integer :: i

  PI = atan(1.d0)*4.d0

  norm_pg1 = (2.d0*alpha/PI)**(1.d0/4.d0)
  if(n.ge.1) then
     do i=1,n
        norm_pg1 = norm_pg1*(4.d0*alpha/(2.d0*i -1.d0))**(1.d0/2.d0)
     end do
  end if

end function norm_pg1

!========================================================================
subroutine permutate_cccc(inttype,i,j,k,l,A,B,C,D,coefr,NBS0,NBS,factor)
!========================================================================
  implicit none
  integer,intent(in) :: inttype,i,j,k,l,A,B,C,D,NBS,NBS0
  real(kind=8),intent(in) ::coefr(NBS0,NBS,2)
  real(kind=8),intent(out) :: factor

  if(inttype.eq.1) then  
     factor = coefr(A,i,1)*coefr(A,j,1)*coefr(A,k,1)*coefr(A,l,1)    ! (AA|AA)
  elseif(inttype.eq.2) then 
     factor = coefr(A,i,1)*coefr(A,j,1)*coefr(C,k,1)*coefr(C,l,1) &  ! (AA|CC)
            + coefr(C,i,1)*coefr(C,j,1)*coefr(A,k,1)*coefr(A,l,1)    ! (CC|AA)
  elseif(inttype.eq.3) then
     factor = coefr(A,i,1)*coefr(A,j,1)*coefr(C,k,1)*coefr(D,l,1) &  ! (AA|CD)
            + coefr(A,i,1)*coefr(A,j,1)*coefr(D,k,1)*coefr(C,l,1) &  ! (AA|DC)
            + coefr(C,i,1)*coefr(D,j,1)*coefr(A,k,1)*coefr(A,l,1) &  ! (CD|AA)
            + coefr(D,i,1)*coefr(C,j,1)*coefr(A,k,1)*coefr(A,l,1)    ! (DC|AA)
  elseif(inttype.eq.4) then
     factor = coefr(A,i,1)*coefr(B,j,1)*coefr(C,k,1)*coefr(C,l,1) &  ! (AB|CC)
            + coefr(B,i,1)*coefr(A,j,1)*coefr(C,k,1)*coefr(C,l,1) &  ! (BA|CC)
            + coefr(C,i,1)*coefr(C,j,1)*coefr(A,k,1)*coefr(B,l,1) &  ! (CC|AB)
            + coefr(C,i,1)*coefr(C,j,1)*coefr(B,k,1)*coefr(A,l,1)    ! (CC|BA)
  elseif(inttype.eq.5) then
     factor = coefr(A,i,1)*coefr(B,j,1)*coefr(A,k,1)*coefr(B,l,1) &  ! (AB|AB)
            + coefr(A,i,1)*coefr(B,j,1)*coefr(B,k,1)*coefr(A,l,1) &  ! (AB|BA)
            + coefr(B,i,1)*coefr(A,j,1)*coefr(A,k,1)*coefr(B,l,1) &  ! (BA|AB)
            + coefr(B,i,1)*coefr(A,j,1)*coefr(B,k,1)*coefr(A,l,1)    ! (BA|BA)
  elseif(inttype.eq.6) then
     factor = coefr(A,i,1)*coefr(B,j,1)*coefr(A,k,1)*coefr(D,l,1) &  ! (AB|AD)
            + coefr(A,i,1)*coefr(B,j,1)*coefr(D,k,1)*coefr(A,l,1) &  ! (AB|DA)
            + coefr(B,i,1)*coefr(A,j,1)*coefr(A,k,1)*coefr(D,l,1) &  ! (BA|AD)
            + coefr(B,i,1)*coefr(A,j,1)*coefr(D,k,1)*coefr(A,l,1) &  ! (BA|DA)
            + coefr(A,i,1)*coefr(D,j,1)*coefr(A,k,1)*coefr(B,l,1) &  ! (AD|AB)
            + coefr(A,i,1)*coefr(D,j,1)*coefr(B,k,1)*coefr(A,l,1) &  ! (AD|BA)
            + coefr(D,i,1)*coefr(A,j,1)*coefr(A,k,1)*coefr(B,l,1) &  ! (DA|AB)
            + coefr(D,i,1)*coefr(A,j,1)*coefr(B,k,1)*coefr(A,l,1)    ! (DA|BA)
  elseif(inttype.eq.7) then
     factor = coefr(A,i,1)*coefr(B,j,1)*coefr(C,k,1)*coefr(D,l,1) &  ! (AB|CD)
            + coefr(A,i,1)*coefr(B,j,1)*coefr(D,k,1)*coefr(C,l,1) &  ! (AB|DC)
            + coefr(B,i,1)*coefr(A,j,1)*coefr(C,k,1)*coefr(D,l,1) &  ! (BA|CD)
            + coefr(B,i,1)*coefr(A,j,1)*coefr(D,k,1)*coefr(C,l,1) &  ! (BA|DC)
            + coefr(C,i,1)*coefr(D,j,1)*coefr(A,k,1)*coefr(B,l,1) &  ! (CD|AB)
            + coefr(C,i,1)*coefr(D,j,1)*coefr(B,k,1)*coefr(A,l,1) &  ! (CD|BA)
            + coefr(D,i,1)*coefr(C,j,1)*coefr(A,k,1)*coefr(B,l,1) &  ! (DC|AB)
            + coefr(D,i,1)*coefr(C,j,1)*coefr(B,k,1)*coefr(A,l,1)    ! (DC|BA)
  else
     stop "inttype should be 1 to 7."
  end if
end subroutine permutate_cccc

!===========================================================
subroutine find_int_type(A,B,C,D,inttype,num_term)
! classify the type of  (AB|CD) integral 
! use symmetry A<->B, C<->D, AB<->CD
! assume A<=B, C<=D, A+B<=C+D, A<=C
!===========================================================
  implicit none
  integer,intent(in) :: A,B,C,D
  integer,intent(out) :: inttype, num_term
  ! num_term for check
  integer :: sum

  if(A.eq.B) then
     if(C.eq.D) then
        if(A.eq.C) then
           inttype = 1 ! (AA|AA)
           num_term = 1
        else
           inttype = 2 ! (AA|CC)
           num_term = 2
        end if
     else 
        inttype = 3 ! (AA|CD)
        num_term = 4
     end if
  else
     if(C.eq.D) then
        inttype = 4 ! (AB|CC)
        num_term = 4
     else
        if(A.eq.C) then
           if(B.eq.D) then
              inttype = 5  ! (AB|AB)
              num_term = 4
           else
              inttype = 6  ! (AB|AD)
              num_term = 8
           end if
        else
           inttype = 7  ! (AB|CD)
           num_term = 8
        end if
     end if
  end if

end subroutine find_int_type

!===============================================================================================
subroutine gauss_int_fourier(posA,alphaA,nx,ny,nz,posB,alphaB,nbarx,nbary,nbarz,vecK,fourier)
! [input]
! posA : position of primitive gaussian A (unnormalized) 
! alphaA : exponent of A
! nx,ny,nz : type of p.g. A 
! posB : position of primitive gaussian B (unnormalized) 
! alphaB : exponenet of B
! nbarx,nbary,nbarz : type of p.g. B
! vecK : wave number vector of Fourier transformation.
! 
! [output]
! fourier : fourier integral  (A|e^{ik.x}|B)  --> complex number
!===============================================================================================
  implicit none

  real(kind=8) :: PI
  
  real(kind=8),intent(in) :: posA(3), posB(3) ! position of center of gaussian 
  real(kind=8),intent(in) :: alphaA, alphaB ! exponents
  real(kind=8),intent(in) :: vecK(3) ! wave number
  integer,intent(in) :: nx,nbarx,ny,nbary,nz,nbarz

  complex(kind=8),intent(out) :: fourier

  real(kind=8) :: c_PAB ! prefactor when two gaussians are combined
  real(kind=8) :: alphaP, posP(3) ! exponents and position of combined gaussian
  real(kind=8) :: vecPA(3), vecPB(3) ! P-A, P-B

  integer :: i,j,k,l
  
  real(kind=8), allocatable, dimension(:,:,:) :: d,e,f
  
  integer :: n_sum
  integer :: nx_sum,ny_sum,nz_sum  ! upper bounds of N,L,M
  
  complex(kind=8) :: int_fourier
  real(kind=8) :: dijk
  complex(kind=8) :: sum
  real(kind=8) :: norm

  PI = atan(1.d0)*4.d0
  
  nx_sum = nx+nbarx
  ny_sum = ny+nbary
  nz_sum = nz+nbarz
  n_sum = nx_sum+ny_sum+nz_sum

  call calc_gaussP(posA,alphaA,posB,alphaB,c_PAB,posP,alphaP,vecPA,vecPB)
  
  norm = 1.d0 ! unnormalized

  allocate(d(0:nx,0:nbarx,0:nx_sum))
  allocate(e(0:ny,0:nbary,0:ny_sum))
  allocate(f(0:nz,0:nbarz,0:nz_sum))

  call calc_d(alphaP,vecPA(1),vecPB(1),nx,nbarx,d)
  call calc_d(alphaP,vecPA(2),vecPB(2),ny,nbary,e)
  call calc_d(alphaP,vecPA(3),vecPB(3),nz,nbarz,f)

  !-------------------------------------
  ! fourier integral
  !-------------------------------------
  sum=(0.d0,0.d0)
  do i=0,nx_sum
     do j=0,ny_sum
        do k=0,nz_sum
           call calc_fourier_integral(i,j,k,alphaP,posP,vecK,int_fourier)
           dijk = d(nx,nbarx,i)*e(ny,nbary,j)*f(nz,nbarz,k)
           sum = sum +dijk*int_fourier
!           write(*,*) sum
        end do
     end do
  end do
  
  sum = c_PAB*norm**2*sum
  fourier = sum

end subroutine gauss_int_fourier

!========================================================================
subroutine calc_fourier_integral(i,j,k,alphaP,posP,vecK,int_fourier)
! calculation of [NLM|e^{ik.x}]
!========================================================================
  implicit none
  integer,intent(in) :: i,j,k 
  real(kind=8),intent(in) :: alphaP,posP(3),vecK(3)
  complex(kind=8),intent(out) :: int_fourier
  real(kind=8) :: PI,k2,kp
  complex(kind=8) :: F000,Fijk
  complex(kind=8),parameter :: IU = (0.d0,1.d0)    ! imaginary unit

  PI = atan(1.d0)*4.d0
  
  k2 = vecK(1)**2+vecK(2)**2+vecK(3)**2  ! k.k
  kp = vecK(1)*posP(1)+vecK(2)*posP(2)+vecK(3)*posP(3) ! k.p
  F000 = (PI/alphaP)**(3.d0/2.d0)*exp(-k2/4.d0/alphaP)*exp(IU*kp)

  int_fourier = (IU*vecK(1))**i *(IU*vecK(2))**j *(IU*vecK(3))**k *F000

  return
end subroutine calc_fourier_integral


!=====================================================================================================
subroutine gauss_int_twoele(posA,alphaA,n1,posB,alphaB,nbar1,posC,alphaC,n2,posD,alphaD,nbar2,twoele)
! calculate two electron integral for primitive gaussians (AB|CD)
!=====================================================================================================
  implicit none
  
  real(kind=8),intent(in) :: posA(3),posB(3),posC(3),posD(3)
  real(kind=8),intent(in) :: alphaA,alphaB,alphaC,alphaD
  integer,intent(in) :: n1(3),nbar1(3),n2(3),nbar2(3)  
  real(kind=8),intent(out) :: twoele

  real(kind=8) :: PI
  integer :: i,j,k
  integer :: i1,j1,k1,i2,j2,k2
  
  real(kind=8), allocatable, dimension(:,:,:) :: d1,e1,f1,d2,e2,f2
  real(kind=8) :: d1ijk,d2ijk

  real(kind=8) :: posP(3),alphaP,c_PAB,vecPA(3),vecPB(3) ! mixing of A and B
  real(kind=8) :: posQ(3),alphaQ,c_QCD,vecQC(3),vecQD(3) ! mixing of C and D
  real(kind=8) :: lambda,vecPQ(3),alphaT,T  ! R_NLM for two electron integral
  real(kind=8), allocatable, dimension(:) :: list_R000j
  real(kind=8), allocatable, dimension(:,:,:) :: R
  integer :: n1_sum(3)  ! upper bounds of N,L,M
  integer :: n2_sum(3)  ! upper bounds of N',L',M'
  integer :: n_sum(3) ! upper bounds of N+N',L+L',M+M'
  integer :: j_max ! upper bound for j when construction R_NLMj
  real(kind=8) :: integral, twoele_tmp

  PI = atan(1.d0)*4.d0

  do i=1,3
     n1_sum(i) = n1(i)+nbar1(i)
     n2_sum(i) = n2(i)+nbar2(i)
     n_sum(i)  = n1_sum(i)+n2_sum(i)
  end do

  call calc_gaussP(posA,alphaA,posB,alphaB,c_PAB,posP,alphaP,vecPA,vecPB)  ! construct P 
  call calc_gaussP(posC,alphaC,posD,alphaD,c_QCD,posQ,alphaQ,vecQC,vecQD)  ! construct Q
  
!!$  write(*,*) "gaussian P"  
!!$  write(*,'(5es16.6)') posP(1),posP(2),posP(3),alphaP,c_PAB
!!$  write(*,*) "gaussian Q"  
!!$  write(*,'(5es16.6)') posQ(1),posQ(2),posQ(3),alphaQ,c_QCD

  ! construct d,e,f,d',e',f'
  allocate(d1(0:n1(1),0:nbar1(1),0:n1_sum(1)))
  allocate(e1(0:n1(2),0:nbar1(2),0:n1_sum(2)))
  allocate(f1(0:n1(3),0:nbar1(3),0:n1_sum(3)))
  allocate(d2(0:n2(1),0:nbar2(1),0:n2_sum(1)))
  allocate(e2(0:n2(2),0:nbar2(2),0:n2_sum(2)))
  allocate(f2(0:n2(3),0:nbar2(3),0:n2_sum(3)))

  call calc_d(alphaP,vecPA(1),vecPB(1),n1(1),nbar1(1),d1)
  call calc_d(alphaP,vecPA(2),vecPB(2),n1(2),nbar1(2),e1)
  call calc_d(alphaP,vecPA(3),vecPB(3),n1(3),nbar1(3),f1)
  call calc_d(alphaQ,vecQC(1),vecQD(1),n2(1),nbar2(1),d2)
  call calc_d(alphaQ,vecQC(2),vecQD(2),n2(2),nbar2(2),e2)
  call calc_d(alphaQ,vecQC(3),vecQD(3),n2(3),nbar2(3),f2)
  
!!$  ! output d's (for check)
!!$  write(*,*)
!!$  write(*,'(2a6,1a16)') "(i,","j)","d1(i,j,k)"
!!$  do i=0,n1(1)
!!$     do j=0,nbar1(1)
!!$        write(*,'(2i6)',advance="no") i,j
!!$        do k=0,n1_sum(1)
!!$           write(*,'(1es16.6)',advance="no") d1(i,j,k) 
!!$        end do
!!$        write(*,*)
!!$     end do
!!$     write(*,*)
!!$  end do
!!$  write(*,'(2a6,1a16)') "(i,","j)","e1(i,j,k)"
!!$  do i=0,n1(2)
!!$     do j=0,nbar1(2)
!!$        write(*,'(2i6)',advance="no") i,j
!!$        do k=0,n1_sum(2)
!!$           write(*,'(1es16.6)',advance="no") e1(i,j,k) 
!!$        end do
!!$        write(*,*)
!!$     end do
!!$     write(*,*)
!!$  end do
!!$  write(*,'(2a6,1a16)') "(i,","j)","f1(i,j,k)"
!!$  do i=0,n1(3)
!!$     do j=0,nbar1(3)
!!$        write(*,'(2i6)',advance="no") i,j
!!$        do k=0,n1_sum(3)
!!$           write(*,'(1es16.6)',advance="no") f1(i,j,k) 
!!$        end do
!!$        write(*,*)
!!$     end do
!!$     write(*,*)
!!$  end do
!!$
!!$  write(*,*)
!!$  write(*,'(2a6,1a16)') "(i,","j)","d2(i,j,k)"
!!$  do i=0,n2(1)
!!$     do j=0,nbar2(1)
!!$        write(*,'(2i6)',advance="no") i,j
!!$        do k=0,n2_sum(1)
!!$           write(*,'(1es16.6)',advance="no") d2(i,j,k) 
!!$        end do
!!$        write(*,*)
!!$     end do
!!$     write(*,*)
!!$  end do
!!$  write(*,'(2a6,1a16)') "(i,","j)","e2(i,j,k)"
!!$  do i=0,n2(2)
!!$     do j=0,nbar2(2)
!!$        write(*,'(2i6)',advance="no") i,j
!!$        do k=0,n2_sum(2)
!!$           write(*,'(1es16.6)',advance="no") e2(i,j,k) 
!!$        end do
!!$        write(*,*)
!!$     end do
!!$     write(*,*)
!!$  end do
!!$  write(*,'(2a6,1a16)') "(i,","j)","f2(i,j,k)"
!!$  do i=0,n2(3)
!!$     do j=0,nbar2(3)
!!$        write(*,'(2i6)',advance="no") i,j
!!$        do k=0,n2_sum(3)
!!$           write(*,'(1es16.6)',advance="no") f2(i,j,k) 
!!$        end do
!!$        write(*,*)
!!$     end do
!!$     write(*,*)
!!$  end do

  ! construct R_{N+N',L+L',M+M'}  (3.34)
  lambda = 2.d0*PI**(5.d0/2.d0)/(alphaP*alphaQ*sqrt(alphaP+alphaQ))  !M&D eq.(3.31)
  T = 0.d0
  do i=1,3
     vecPQ(i) = posP(i)-posQ(i)
     T = T +vecPQ(i)**2
  end do
  alphaT = alphaP*alphaQ/(alphaP+alphaQ)
  T = alphaT*T  ! M&D eq.(3.32)

!!$  write(*,'(6a16)')  'PQx','PQy','PQz','lambda','alphaT','T'
!!$  write(*,'(6es16.6)') vecPQ(1),vecPQ(2),vecPQ(3),lambda,alphaT,T

  allocate(R(0:n_sum(1),0:n_sum(2),0:n_sum(3))) ! for two-electron integral  (need R_{N+N',L,M} , R_{N,L+L',M} and R_{N,L,M+M'})
  j_max = n_sum(1)+n_sum(2)+n_sum(3) ! for two-electron integral
  allocate(list_R000j(0:j_max))
  call calc_R000j(alphaT,T,j_max,list_R000j)
  call calc_R(vecPQ,j_max,list_R000j,n_sum(1),n_sum(2),n_sum(3),R) 
    
!!$  write(*,*) "output R"
!!$  write(*,'(3a5,1a16)') "i","j","k","R(i,j,k)"
!!$  do i=0,n_sum(1) ! loop for all N
!!$     do j=0,n_sum(2) ! loop for all L
!!$        do k=0,n_sum(3) ! loop for all M
!!$           write(*,'(3i5,1es16.6)') i,j,k,R(i,j,k)
!!$        end do
!!$     end do
!!$  end do

  twoele = 0.d0
  do i1=0,n1_sum(1)  ! N
     do j1=0,n1_sum(2)  ! L
        do k1=0,n1_sum(3)  ! M
           d1ijk = d1(n1(1),nbar1(1),i1) * e1(n1(2),nbar1(2),j1) * f1(n1(3),nbar1(3),k1)
           twoele_tmp = 0.d0 
           do i2=0,n2_sum(1)  ! N'
              do j2=0,n2_sum(2) ! L'
                 do k2=0,n2_sum(3) ! M'
!                    integral = lambda*(-1.d0)**(i2+j2+k2)*R(i1+i2,j1+j2,k1+k2)  ! eq.(3.34)
                    integral = (-1.d0)**(i2+j2+k2)*R(i1+i2,j1+j2,k1+k2)  ! eq.(3.34)
                    d2ijk = d2(n2(1),nbar2(1),i2) * e2(n2(2),nbar2(2),j2) * f2(n2(3),nbar2(3),k2)
!                    twoele = twoele +d1ijk*d2ijk*integral
                    twoele_tmp = twoele_tmp +d2ijk*integral
                 end do
              end do
           end do
           twoele = twoele + twoele_tmp * d1ijk
        end do
     end do
  end do
!  twoele = c_PAB*c_QCD*twoele
  twoele = c_PAB*c_QCD*twoele * lambda

end subroutine gauss_int_twoele

!========================================================================
!program gauss_int
subroutine gauss_int(posA,alphaA,nx,ny,nz,posB,alphaB,nbarx,nbary,nbarz,posC,overlap,nucatt,ef)
! [input]
! posA : position of primitive gaussian A (unnormalized) 
! alphaA : exponent of A
! nx,ny,nz : type of p.g. A 
! posB : position of primitive gaussian B (unnormalized) 
! alphaB : exponenet of B
! nbarx,nbary,nbarz : type of p.g. B
! posC : position of nucleus (or where you want to evaluate electric field)
! 
! [output]
! overlap : overlap integral
! nucatt : nuclear attraction integral  (A|1/rC|B)  (w/o negative sign)
! ef(1~3) : x,y,z components of electric field
!========================================================================
  implicit none

  real(kind=8) :: PI
  !  integer :: nx,nbarx,ny,nbary,nz,nbarz ! angular momentum
  
  real(kind=8),intent(in) :: posA(3), posB(3) ! position of center of gaussian 
  real(kind=8),intent(in) :: alphaA, alphaB ! exponents
  real(kind=8),intent(in) :: posC(3) ! position of nucleus
  integer,intent(in) :: nx,nbarx,ny,nbary,nz,nbarz

  real(kind=8),intent(out) :: overlap,nucatt,ef(3)


  real(kind=8) :: c_PAB ! prefactor when two gaussians are combined
  real(kind=8) :: alphaP, posP(3) ! exponents and position of combined gaussian
  real(kind=8) :: vecPA(3), vecPB(3) ! P-A, P-B
  real(kind=8) :: vecPC(3) ! P-C
  
  integer :: i,j,k,l
  
  !  real(kind=8) :: d(0:nx,0:nbarx,0:nx+nbarx)
  ! real(kind=8) :: d(0:n_max,0:nbar_max,0:n_max+nbar_max)
  real(kind=8), allocatable, dimension(:,:,:) :: d,e,f
  real(kind=8), allocatable, dimension(:,:,:) :: R
!  real(kind=8), allocatable, dimension(:) :: list_FjT
  real(kind=8), allocatable, dimension(:) :: list_R000j
  
  real(kind=8) :: T
  
  integer :: n_sum
  integer :: nx_sum,ny_sum,nz_sum  ! upper bounds of N,L,M
  integer :: j_max ! upper bound for j when construction R_NLMj
!  integer :: n_max,nbar_max
  
  real(kind=8) :: sum, integral, dijk
  real(kind=8) :: norm
!  real(kind=8) :: ef(3)

  PI = atan(1.d0)*4.d0
  
!  nx = 2; nbarx = 2
!  ny = 0; nbary = 0
!  nz = 0; nbarz = 0
  nx_sum = nx+nbarx
  ny_sum = ny+nbary
  nz_sum = nz+nbarz
  n_sum = nx_sum+ny_sum+nz_sum
!  n_max = max(nx,ny,nz)     
!  nbar_max = max(nbarx,nbary,nbarz)
!!$  write(*,'(6a6)') 'nx',' nbarx','ny',' nbary','nz',' nbarz'
!!$  write(*,'(6i6)') nx,nbarx,ny,nbary,nz,nbarz
!  write(*,'(2a10)') 'n_max','nbar_max'
!  write(*,'(2i10)') n_max,nbar_max

!  alphaA = 0.5d0
!  alphaB = 1.0d0
!  alphaB = alphaA
!  alphaA = alphaB
!  posA(1) = 1.d0;  posA(2) = 2.d0;  posA(3) = 3.d0
!  posB(1) = 2.d0;  posB(2) = 3.d0;  posB(3) = 4.d0
!  posB(1:3) = posA(1:3)
!  posA(1:3) = posB(1:3)
  call calc_gaussP(posA,alphaA,posB,alphaB,c_PAB,posP,alphaP,vecPA,vecPB)
  
!  posC(1) = 3.d0; posC(2) = 4.d0; posC(3) = 5.d0
  
!!$  write(*,'(1a10,1f16.6)') "c_PAB",c_PAB
!!$  write(*,'(3a16)') "alphaA","alphaB","alphaP"
!!$  write(*,'(3f16.6)') alphaA, alphaB, alphaP
!!$  write(*,'(5a16)')"posA","posB","posP" ,"vecPA","vecPB"
!!$  do i=1,3
!!$     write(*,'(5f16.6)') posA(i), posB(i), posP(i), vecPA(i), vecPB(i)
!!$  end do
  
 
  norm = 1.d0 ! unnormalized
!  norm = (8.d0*alphaA**3/PI**3)**0.25d0 ! for s
!  norm = (128.d0*alphaA**5/PI**3)**0.25d0  ! for p
!  norm = (2048.d0*alphaA**7/PI**3)**0.25d0 ! for dxy, dxz, dyz
!  norm = (2048.d0*alphaA**7/PI**3/9.d0)**0.25d0 ! for dxx,dyy,dzz

  allocate(d(0:nx,0:nbarx,0:nx_sum))
  allocate(e(0:ny,0:nbary,0:ny_sum))
  allocate(f(0:nz,0:nbarz,0:nz_sum))

  call calc_d(alphaP,vecPA(1),vecPB(1),nx,nbarx,d)
  call calc_d(alphaP,vecPA(2),vecPB(2),ny,nbary,e)
  call calc_d(alphaP,vecPA(3),vecPB(3),nz,nbarz,f)
  
!!$  ! output d's (for check)
!!$  write(*,*)
!!$  write(*,'(2a6,1a16)') "(i,","j)","d(i,j,k)"
!!$  do i=0,nx
!!$     do j=0,nbarx
!!$        write(*,'(2i6)',advance="no") i,j
!!$        do k=0,nx_sum
!!$           write(*,'(1es16.6)',advance="no") d(i,j,k) 
!!$        end do
!!$        write(*,*)
!!$     end do
!!$     write(*,*)
!!$  end do
!!$  write(*,'(2a6,1a16)') "(i,","j)","e(i,j,k)"
!!$  do i=0,ny
!!$     do j=0,nbary
!!$        write(*,'(2i6)',advance="no") i,j
!!$        do k=0,ny_sum
!!$           write(*,'(1es16.6)',advance="no") e(i,j,k) 
!!$        end do
!!$        write(*,*)
!!$     end do
!!$     write(*,*)
!!$  end do
!!$  write(*,'(2a6,1a16)') "(i,","j)","f(i,j,k)"
!!$  do i=0,nz
!!$     do j=0,nbarz
!!$        write(*,'(2i6)',advance="no") i,j
!!$        do k=0,nz_sum
!!$           write(*,'(1es16.6)',advance="no") f(i,j,k) 
!!$        end do
!!$        write(*,*)
!!$     end do
!!$     write(*,*)
!!$  end do
  
  !-------------------------------------
  ! test overlapping integral
  !-------------------------------------
  sum = 0.d0
  do i=0,nx_sum
     do j=0,ny_sum
        do k=0,nz_sum
           call calc_overlap_integral(i,j,k,alphaP,integral)
           dijk = d(nx,nbarx,i)*e(ny,nbary,j)*f(nz,nbarz,k)
           sum = sum +dijk*integral
        end do
     end do
  end do

  sum = c_PAB*norm**2*sum
  overlap = sum

!!$  write(*,*) "overlapping integral"
!!$  write(*,'(1es16.6)') sum
!!$
!!$  write(*,*)
  !-------------------------------------
  ! test nuclear attraction integral
  !-------------------------------------

!  allocate(R(0:nx+nbarx,0:ny+nbary,0:nz+nbarz))
  allocate(R(0:nx_sum+1,0:ny_sum+1,0:nz_sum+1)) ! for Electric field (need R_{N+1,L,M} , R_{N,L+1,M} and R_{N,L,M+1})
!  allocate(R(0:nx_sum*2,0:ny_sum*2,0:nz_sum*2)) ! for two-electron integral  (need R_{N+N',L,M} , R_{N,L+L',M} and R_{N,L,M+M'})
  j_max =  n_sum+3 ! for Electric field
!  j_max = n_sum*2 ! for two-electron integral
  allocate(list_R000j(0:j_max))

  call calc_PC_and_T(posC,posP,alphaP,vecPC,T)
  call calc_R000j(alphaP,T,j_max,list_R000j)

!!$  write(*,'(1a10,3es16.6)') "vecPC",vecPC
!!$  write(*,'(1a10,1es16.6)') "T",T
!!$
!!$  write(*,'(1a5,1a16)') "j","R000j    "
!!$  do i=0,j_max
!!$!     write(*,*) i,list_FjT(i)
!!$     write(*,'(1i5,1es16.6)') i,list_R000j(i)
!!$  end do

! calc_R only depend on nx+nbarx, ny+nbary, nz+nbarz 
  call calc_R(vecPC,j_max,list_R000j,nx_sum+1,ny_sum+1,nz_sum+1,R) ! for Electric field
!  call calc_R(vecPC,j_max,list_R000j,nx_sum*2,ny_sum*2,nz_sum*2,R) ! for two-electron integral


!!$  write(*,*) "output R"
!!$  write(*,'(3a5,1a16)') "i","j","k","R(i,j,k)"
!!$  do i=0,nx+nbarx+1 ! loop for all N
!!$     do j=0,ny+nbary+1 ! loop for all L
!!$        do k=0,nz+nbarz+1 ! loop for all M
!!$           write(*,'(3i5,1es16.6)') i,j,k,R(i,j,k)
!!$        end do
!!$     end do
!!$  end do

  sum = 0.d0
  do i=0,nx_sum
     do j=0,ny_sum
        do k=0,nz_sum
           integral = (2.d0*PI/alphaP)*R(i,j,k)  !  [NLM | rc^-1]
           dijk = d(nx,nbarx,i)*e(ny,nbary,j)*f(nz,nbarz,k)
           sum = sum +dijk*integral
        end do
     end do
  end do

  sum = c_PAB*norm**2*sum
!!$  write(*,*) "nuclear attraction integral"
!!$  write(*,'(1es16.6)') sum
  nucatt = sum

  !-------------------------------------
  ! test electric field integral
  !-------------------------------------

  ef(:) = 0.d0
  do i=0,nx_sum
     do j=0,ny_sum
        do k=0,nz_sum
           dijk = d(nx,nbarx,i)*e(ny,nbary,j)*f(nz,nbarz,k)
           ef(1) = ef(1) +dijk*(-(2.d0*PI/alphaP)*R(i+1,j,k)) !  [NLM | xc rc^-3] (3.21)
           ef(2) = ef(2) +dijk*(-(2.d0*PI/alphaP)*R(i,j+1,k)) !  [NLM | yc rc^-3] (3.21)
           ef(3) = ef(3) +dijk*(-(2.d0*PI/alphaP)*R(i,j,k+1)) !  [NLM | zc rc^-3] (3.21)
        end do
     end do
  end do

  do i=1,3
     ef(i) = c_PAB*norm**2*ef(i)
  end do
!!$  write(*,*) "electric field integral"
!!$  do i=1,3
!!$     write(*,'(1es16.6)',advance="no")  ef(i)
!!$  end do
!!$  write(*,*)
!!$
!!$  write(*,*)
end subroutine gauss_int

!===============================================================================================
subroutine gauss_int_moment(posA,alphaA,nx,ny,nz,posB,alphaB,nbarx,nbary,nbarz,posC,moment)
! [input]
! posA : position of primitive gaussian A (unnormalized) 
! alphaA : exponent of A
! nx,ny,nz : type of p.g. A 
! posB : position of primitive gaussian B (unnormalized) 
! alphaB : exponenet of B
! nbarx,nbary,nbarz : type of p.g. B
! posC : position related to moment
! 
! [output]
! moment : moment integral  (A|xC|B)  (xC = x-Cx) , (A|yC|B), (A|zC|B)  p.223 of M&D
!===============================================================================================
  implicit none

  real(kind=8) :: PI
  
  real(kind=8),intent(in) :: posA(3), posB(3) ! position of center of gaussian 
  real(kind=8),intent(in) :: alphaA, alphaB ! exponents
  real(kind=8),intent(in) :: posC(3) 
  integer,intent(in) :: nx,nbarx,ny,nbary,nz,nbarz

  real(kind=8),intent(out) :: moment(3)

  real(kind=8) :: c_PAB ! prefactor when two gaussians are combined
  real(kind=8) :: alphaP, posP(3) ! exponents and position of combined gaussian
  real(kind=8) :: vecPA(3), vecPB(3) ! P-A, P-B
  real(kind=8) :: vecPC(3) ! P-C
  real(kind=8) :: T

  integer :: i,j,k,l
  
  real(kind=8), allocatable, dimension(:,:,:) :: d,e,f
  
  integer :: n_sum
  integer :: nx_sum,ny_sum,nz_sum  ! upper bounds of N,L,M
  
  real(kind=8) :: int_mom(3), dijk
  real(kind=8) :: sum(3)
  real(kind=8) :: norm

  PI = atan(1.d0)*4.d0
  
  nx_sum = nx+nbarx
  ny_sum = ny+nbary
  nz_sum = nz+nbarz
  n_sum = nx_sum+ny_sum+nz_sum

  call calc_gaussP(posA,alphaA,posB,alphaB,c_PAB,posP,alphaP,vecPA,vecPB)
  
  norm = 1.d0 ! unnormalized

  allocate(d(0:nx,0:nbarx,0:nx_sum))
  allocate(e(0:ny,0:nbary,0:ny_sum))
  allocate(f(0:nz,0:nbarz,0:nz_sum))

  call calc_d(alphaP,vecPA(1),vecPB(1),nx,nbarx,d)
  call calc_d(alphaP,vecPA(2),vecPB(2),ny,nbary,e)
  call calc_d(alphaP,vecPA(3),vecPB(3),nz,nbarz,f)

  call calc_PC_and_T(posC,posP,alphaP,vecPC,T)  ! compute vecPC. we do not use T
!  write(*,*) vecPC
!  stop

  !-------------------------------------
  ! moment integral
  !-------------------------------------
  do l=1,3
     sum(l)=0.d0
  end do
  do i=0,nx_sum
     do j=0,ny_sum
        do k=0,nz_sum
           call calc_moment_integral(i,j,k,alphaP,vecPC,int_mom)
           dijk = d(nx,nbarx,i)*e(ny,nbary,j)*f(nz,nbarz,k)
           do l=1,3
              sum(l) = sum(l) +dijk*int_mom(l)
           end do
        end do
     end do
  end do

  do l=1,3
     sum(l) = c_PAB*norm**2*sum(l)
!     write(*,*) l,sum(l)
     moment(l) = sum(l)
  end do

end subroutine gauss_int_moment

!========================================================================
subroutine calc_moment_integral(i,j,k,alphaP,vecPC,int_mom)
! calculation of 
! [NLM|xC], [NLM|yC], [NLM|zC] in McMurchie & Davidson Eq. (3.6) for moment integral
!========================================================================
  implicit none
  integer,intent(in) :: i,j,k 
  real(kind=8),intent(in) :: alphaP,vecPC(3)
  real(kind=8),intent(out) :: int_mom(3)
  real(kind=8) :: PI

  PI = atan(1.d0)*4.d0

  ! [NLM|xC]
  if((j.eq.0).and.(k.eq.0)) then
     if(i.eq.0) then
        int_mom(1) = vecPC(1)*(PI/alphaP)**1.5d0
     elseif(i.eq.1) then
        int_mom(1) = (PI/alphaP)**1.5d0
     else
        int_mom(1) = 0.d0
     end if
  else
     int_mom(1) = 0.d0
  end if

  ! [NLM|yC]
  if((i.eq.0).and.(k.eq.0)) then
     if(j.eq.0) then
        int_mom(2) = vecPC(2)*(PI/alphaP)**1.5d0
     elseif(j.eq.1) then
        int_mom(2) = (PI/alphaP)**1.5d0
     else
        int_mom(2) = 0.d0
     end if
  else
     int_mom(2) = 0.d0
  end if

  ! [NLM|zC]
  if((i.eq.0).and.(j.eq.0)) then
     if(k.eq.0) then
        int_mom(3) = vecPC(3)*(PI/alphaP)**1.5d0
     elseif(k.eq.1) then
        int_mom(3) = (PI/alphaP)**1.5d0
     else
        int_mom(3) = 0.d0
     end if
  else
     int_mom(3) = 0.d0
  end if

  return
end subroutine calc_moment_integral


!==================================================================================
subroutine gauss_int_KE(posA,alphaA,nx,ny,nz,posB,alphaB,nbarx,nbary,nbarz,ene_kin)
! calculate kinetic energy integral of p.g.
! use formula (2.36)(2.37) of M&D
! [input]
! posA : position of primitive gaussian A (unnormalized) 
! alphaA : exponent of A
! nx,ny,nz : type of p.g. A 
! posB : position of primitive gaussian B (unnormalized) 
! alphaB : exponenet of B
! nbarx,nbary,nbarz : type of p.g. B
! 
! [output]
! ene_kin : (-1/2)Lap   --> include the factor (1/2)
!==================================================================================
  implicit none

  real(kind=8) :: PI
  !  integer :: nx,nbarx,ny,nbary,nz,nbarz ! angular momentum
  
  real(kind=8),intent(in) :: posA(3), posB(3) ! position of center of gaussian 
  real(kind=8),intent(in) :: alphaA, alphaB ! exponents
  integer,intent(in) :: nx,nbarx,ny,nbary,nz,nbarz

  real(kind=8),intent(out) :: ene_kin


  real(kind=8) :: c_PAB ! prefactor when two gaussians are combined (E_IJ)
  real(kind=8) :: alphaP, posP(3) ! exponents and position of combined gaussian
  real(kind=8) :: vecPA(3), vecPB(3) ! P-A, P-B
  integer :: i,j,k,l
  
  real(kind=8), allocatable, dimension(:,:,:) :: d,e,f  ! d(n,nbar,N)
  
  integer :: n_sum
  integer :: nx_sum,ny_sum,nz_sum  ! upper bounds of N,L,M
  
  real(kind=8) :: sum, integral, dijk
  real(kind=8) :: t,txx,tyy,tzz ! Eq.(2.36)

  PI = atan(1.d0)*4.d0
  
  nx_sum = nx+nbarx
  ny_sum = ny+nbary
  nz_sum = nz+nbarz
  n_sum = nx_sum+ny_sum+nz_sum

  call calc_gaussP(posA,alphaA,posB,alphaB,c_PAB,posP,alphaP,vecPA,vecPB)
  
!  allocate(d(0:nx,0:nbarx,0:nx_sum))
!  allocate(e(0:ny,0:nbary,0:ny_sum))
!  allocate(f(0:nz,0:nbarz,0:nz_sum))
  allocate(d(0:nx+1,0:nbarx+1,0:nx_sum+2))
  allocate(e(0:ny+1,0:nbary+1,0:ny_sum+2))
  allocate(f(0:nz+1,0:nbarz+1,0:nz_sum+2))

  call calc_d(alphaP,vecPA(1),vecPB(1),nx+1,nbarx+1,d)
  call calc_d(alphaP,vecPA(2),vecPB(2),ny+1,nbary+1,e)
  call calc_d(alphaP,vecPA(3),vecPB(3),nz+1,nbarz+1,f)
  
  !-------------------------------------
  ! kinetic energy integral
  !-------------------------------------

!!$  do i=0,nx_sum
!!$     do j=0,ny_sum
!!$        do k=0,nz_sum
!!$           call calc_overlap_integral(i,j,k,alphaP,integral)
!!$           dijk = d(nx,nbarx,i)*e(ny,nbary,j)*f(nz,nbarz,k)
!!$           sum = sum +dijk*integral
!!$        end do
!!$     end do
!!$  end do
  ! We need only N=0,L=0,M=0 (p.224) (i=0,j=0,k=0 in my notation)
  ! use formula (2.37) note that d(-1,-1,0) etc are not defined.
  if(nx.eq.0 .and. nbarx.eq.0) then
     txx = 4.d0*alphaA*alphaB*d(nx+1,nbarx+1,0) 
  else if(nx.eq.0) then
     txx = 4.d0*alphaA*alphaB*d(nx+1,nbarx+1,0) -2.d0*nbarx*alphaA*d(nx+1,nbarx-1,0)
  else if(nbarx.eq.0) then
     txx = 4.d0*alphaA*alphaB*d(nx+1,nbarx+1,0) -2.d0*nx*alphaB*d(nx-1,nbarx+1,0)
  else
     txx = 4.d0*alphaA*alphaB*d(nx+1,nbarx+1,0) -2.d0*nbarx*alphaA*d(nx+1,nbarx-1,0) & 
          -2.d0*nx*alphaB*d(nx-1,nbarx+1,0) +nx*nbarx*d(nx-1,nbarx-1,0) 
  end if
  txx = txx*e(ny,nbary,0)*f(nz,nbarz,0)

  if(ny.eq.0 .and. nbary.eq.0) then
     tyy = 4.d0*alphaA*alphaB*e(ny+1,nbary+1,0) 
  else if(ny.eq.0) then
     tyy = 4.d0*alphaA*alphaB*e(ny+1,nbary+1,0) -2.d0*nbary*alphaA*e(ny+1,nbary-1,0)
  else if(nbary.eq.0) then
     tyy = 4.d0*alphaA*alphaB*e(ny+1,nbary+1,0) -2.d0*ny*alphaB*e(ny-1,nbary+1,0)
  else
     tyy = 4.d0*alphaA*alphaB*e(ny+1,nbary+1,0) -2.d0*nbary*alphaA*e(ny+1,nbary-1,0) & 
          -2.d0*ny*alphaB*e(ny-1,nbary+1,0) +ny*nbary*e(ny-1,nbary-1,0) 
  end if
  tyy = tyy*d(nx,nbarx,0)*f(nz,nbarz,0)

  if(nz.eq.0 .and. nbarz.eq.0) then
     tzz = 4.d0*alphaA*alphaB*f(nz+1,nbarz+1,0) 
  else if(nz.eq.0) then
     tzz = 4.d0*alphaA*alphaB*f(nz+1,nbarz+1,0) -2.d0*nbarz*alphaA*f(nz+1,nbarz-1,0)
  else if(nbarz.eq.0) then
     tzz = 4.d0*alphaA*alphaB*f(nz+1,nbarz+1,0) -2.d0*nz*alphaB*f(nz-1,nbarz+1,0)
  else
     tzz = 4.d0*alphaA*alphaB*f(nz+1,nbarz+1,0) -2.d0*nbarz*alphaA*f(nz+1,nbarz-1,0) & 
          -2.d0*nz*alphaB*f(nz-1,nbarz+1,0) +nz*nbarz*f(nz-1,nbarz-1,0) 
  end if
  tzz = tzz*d(nx,nbarx,0)*e(ny,nbary,0)
  
  t = (txx+tyy+tzz)*c_PAB

  ene_kin = t*(PI/alphaP)**1.5d0
  ene_kin = ene_kin/2.d0

end subroutine gauss_int_KE

!====================================================================================
subroutine gauss_int_relKE(posA,alphaA,nx,ny,nz,posB,alphaB,nbarx,nbary,nbarz,grad)
! calculate kinetic energy integral of p.g. for relativistic theory
! (or gradient)
! use formula (2.33)(2.34) of M&D
! [input]
! posA : position of primitive gaussian A (unnormalized) 
! alphaA : exponent of A
! nx,ny,nz : type of p.g. A 
! posB : position of primitive gaussian B (unnormalized) 
! alphaB : exponenet of B
! nbarx,nbary,nbarz : type of p.g. B
! 
! [output]
! grad(3)  : g*nabla g (see (2.29))
!
! 2010.1.14 programmed based on subroutine gauss_int_KE (non-rel version)
!====================================================================================
  implicit none

  real(kind=8) :: PI
  !  integer :: nx,nbarx,ny,nbary,nz,nbarz ! angular momentum
  
  real(kind=8),intent(in) :: posA(3), posB(3) ! position of center of gaussian 
  real(kind=8),intent(in) :: alphaA, alphaB ! exponents
  integer,intent(in) :: nx,nbarx,ny,nbary,nz,nbarz

  real(kind=8),intent(out) :: grad(3)


  real(kind=8) :: c_PAB ! prefactor when two gaussians are combined (E_IJ)
  real(kind=8) :: alphaP, posP(3) ! exponents and position of combined gaussian
  real(kind=8) :: vecPA(3), vecPB(3) ! P-A, P-B
  integer :: i,j,k,l
  
  real(kind=8), allocatable, dimension(:,:,:) :: d,e,f  ! d(n,nbar,N)
  
  integer :: n_sum
  integer :: nx_sum,ny_sum,nz_sum  ! upper bounds of N,L,M
  
  real(kind=8) :: sum, integral, dijk
!  real(kind=8) :: t,txx,tyy,tzz ! Eq.(2.36)
  real(kind=8) :: gx,gy,gz ! Eq.(2.34)

  PI = atan(1.d0)*4.d0
  
  nx_sum = nx+nbarx
  ny_sum = ny+nbary
  nz_sum = nz+nbarz
  n_sum = nx_sum+ny_sum+nz_sum

  call calc_gaussP(posA,alphaA,posB,alphaB,c_PAB,posP,alphaP,vecPA,vecPB)
  
!  allocate(d(0:nx,0:nbarx,0:nx_sum))
!  allocate(e(0:ny,0:nbary,0:ny_sum))
!  allocate(f(0:nz,0:nbarz,0:nz_sum))
  allocate(d(0:nx,0:nbarx+1,0:nx_sum+1))
  allocate(e(0:ny,0:nbary+1,0:ny_sum+1))
  allocate(f(0:nz,0:nbarz+1,0:nz_sum+1))

  call calc_d(alphaP,vecPA(1),vecPB(1),nx,nbarx+1,d)
  call calc_d(alphaP,vecPA(2),vecPB(2),ny,nbary+1,e)
  call calc_d(alphaP,vecPA(3),vecPB(3),nz,nbarz+1,f)

  !-------------------------------------
  ! kinetic energy integral
  !-------------------------------------

!!$  do i=0,nx_sum
!!$     do j=0,ny_sum
!!$        do k=0,nz_sum
!!$           call calc_overlap_integral(i,j,k,alphaP,integral)
!!$           dijk = d(nx,nbarx,i)*e(ny,nbary,j)*f(nz,nbarz,k)
!!$           sum = sum +dijk*integral
!!$        end do
!!$     end do
!!$  end do

  ! NOTES for non-rel version
  ! We need only N=0,L=0,M=0 (p.224) (i=0,j=0,k=0 in my notation)
  ! use formula (2.37) note that d(-1,-1,0) etc are not defined.
  ! NOTES for rel version
  ! We need only N=0,L=0,M=0 (p.224) (i=0,j=0,k=0 in my notation)
  ! use formula (2.34) note that d(0,-1,0) etc are not defined.

  if(nbarx.eq.0) then
     gx = -2.d0*alphaB*d(nx,nbarx+1,0)
  else
     gx = -2.d0*alphaB*d(nx,nbarx+1,0) +nbarx*d(nx,nbarx-1,0)
  end if
  gx = gx*e(ny,nbary,0)*f(nz,nbarz,0)*c_PAB

  if(nbary.eq.0) then
     gy = -2.d0*alphaB*e(ny,nbary+1,0)
  else
     gy = -2.d0*alphaB*e(ny,nbary+1,0) +nbary*e(ny,nbary-1,0)
  end if
  gy = gy*d(nx,nbarx,0)*f(nz,nbarz,0)*c_PAB

  if(nbarz.eq.0) then
     gz = -2.d0*alphaB*f(nz,nbarz+1,0)
  else
     gz = -2.d0*alphaB*f(nz,nbarz+1,0) +nbarz*f(nz,nbarz-1,0)
  end if
  gz = gz*d(nx,nbarx,0)*e(ny,nbary,0)*c_PAB

  grad(1) = gx*(PI/alphaP)**1.5d0 
  grad(2) = gy*(PI/alphaP)**1.5d0 
  grad(3) = gz*(PI/alphaP)**1.5d0 

end subroutine gauss_int_relKE

!========================================================================
subroutine gauss_int_overlap(posA,alphaA,nx,ny,nz,posB,alphaB,nbarx,nbary,nbarz,overlap)
! calculate overlap of p.g.
! [input]
! posA : position of primitive gaussian A (unnormalized) 
! alphaA : exponent of A
! nx,ny,nz : type of p.g. A 
! posB : position of primitive gaussian B (unnormalized) 
! alphaB : exponenet of B
! nbarx,nbary,nbarz : type of p.g. B
! 
! [output]
! overlap : overlap integral
!========================================================================
  implicit none

  real(kind=8) :: PI
  !  integer :: nx,nbarx,ny,nbary,nz,nbarz ! angular momentum
  
  real(kind=8),intent(in) :: posA(3), posB(3) ! position of center of gaussian 
  real(kind=8),intent(in) :: alphaA, alphaB ! exponents
  integer,intent(in) :: nx,nbarx,ny,nbary,nz,nbarz

  real(kind=8),intent(out) :: overlap


  real(kind=8) :: c_PAB ! prefactor when two gaussians are combined
  real(kind=8) :: alphaP, posP(3) ! exponents and position of combined gaussian
  real(kind=8) :: vecPA(3), vecPB(3) ! P-A, P-B
  integer :: i,j,k,l
  
  real(kind=8), allocatable, dimension(:,:,:) :: d,e,f
  
  integer :: n_sum
  integer :: nx_sum,ny_sum,nz_sum  ! upper bounds of N,L,M
  
  real(kind=8) :: sum, integral, dijk

  PI = atan(1.d0)*4.d0
  
  nx_sum = nx+nbarx
  ny_sum = ny+nbary
  nz_sum = nz+nbarz
  n_sum = nx_sum+ny_sum+nz_sum

  call calc_gaussP(posA,alphaA,posB,alphaB,c_PAB,posP,alphaP,vecPA,vecPB)
  
  allocate(d(0:nx,0:nbarx,0:nx_sum))
  allocate(e(0:ny,0:nbary,0:ny_sum))
  allocate(f(0:nz,0:nbarz,0:nz_sum))

  call calc_d(alphaP,vecPA(1),vecPB(1),nx,nbarx,d)
  call calc_d(alphaP,vecPA(2),vecPB(2),ny,nbary,e)
  call calc_d(alphaP,vecPA(3),vecPB(3),nz,nbarz,f)
  
  !-------------------------------------
  ! test overlapping integral
  !-------------------------------------
  sum = 0.d0
  do i=0,nx_sum
     do j=0,ny_sum
        do k=0,nz_sum
           call calc_overlap_integral(i,j,k,alphaP,integral)
           dijk = d(nx,nbarx,i)*e(ny,nbary,j)*f(nz,nbarz,k)
           sum = sum +dijk*integral
        end do
     end do
  end do

  sum = c_PAB*sum
  overlap = sum

end subroutine gauss_int_overlap

!!$!========================================================================
!!$subroutine calc_nucatt_integral(i,j,k,alphaP,integral)
!!$! calculation of [NLM | rc^-1] in McMurchie & Davidson Eq. (3.14)
!!$! for nuclear attraction integral
!!$!========================================================================
!!$  use global
!!$  implicit none
!!$  integer,intent(in) :: i,j,k 
!!$  real(kind=8),intent(in) :: alphaP
!!$  real(kind=8),intent(out) :: integral
!!$
!!$  if((i.eq.0).and.(j.eq.0).and.(k.eq.0)) then
!!$     integral = (PI/alphaP)**1.5d0
!!$  else
!!$     integral = 0.d0
!!$  end if
!!$  return
!!$end subroutine calc_nucatt_integral

!========================================================================
subroutine calc_overlap_integral(i,j,k,alphaP,integral)
! calculation of 
! [NLM | 1] in McMurchie & Davidson Eq. (3.4)  for overlapping integral
!========================================================================
  implicit none
  integer,intent(in) :: i,j,k 
  real(kind=8),intent(in) :: alphaP
  real(kind=8),intent(out) :: integral
  real(kind=8) :: PI

  PI = atan(1.d0)*4.d0

  if((i.eq.0).and.(j.eq.0).and.(k.eq.0)) then
     integral = (PI/alphaP)**1.5d0
  else
     integral = 0.d0
  end if
  return
end subroutine calc_overlap_integral




!====================================================================================
!subroutine calc_R(vecPC,n_sum,list_R000j,nx,nbarx,ny,nbary,nz,nbarz,table_R)
subroutine calc_R(vecPC,n_sum,list_R000j,nx_sum,ny_sum,nz_sum,table_R)
! calculate R_NMLj using Eq(4.4)~(4.8)
! only depend on nx+nbarx, ny+nbary, nz+nbarz 
!====================================================================================
  implicit none
  real(kind=8),intent(in) :: vecPC(3)
!  integer,intent(in) :: n_sum,nx,nbarx,ny,nbary,nz,nbarz
  integer,intent(in) :: n_sum,nx_sum,ny_sum,nz_sum
  real(kind=8),intent(in) :: list_R000j(0:n_sum)
!  real(kind=8),intent(out) :: table_R(0:nx+nbarx,0:ny+nbary,0:nz+nbarz)
  real(kind=8),intent(out) :: table_R(0:nx_sum,0:ny_sum,0:nz_sum)
  integer :: i,j,k,l
!  real(kind=8) :: table_Rj(0:nx+nbarx,0:ny+nbary,0:nz+nbarz,0:n_sum)
  real(kind=8) :: table_Rj(0:nx_sum,0:ny_sum,0:nz_sum,0:n_sum)
!  real(kind=8) :: R(0:nx+nbarx,0:ny+nbary,0:nz+nbarz)
  real(kind=8) :: R(0:nx_sum,0:ny_sum,0:nz_sum)
  real(kind=8) :: a,b,c

  a = vecPC(1); b = vecPC(2); c = vecPC(3)
!  write(*,*) a,b,c
  table_Rj = 0.d0

  ! N=0,L=0,M=0  (Eq.(4.4))
  do i=0,n_sum
     table_Rj(0,0,0,i) = list_R000j(i)
!     write(*,*) i,table_Rj(0,0,0,i)
  end do
  
  ! use loop variable i for N, j for L, k for M, l for j

  ! Generate nonzero M for N=0,L=0 using Eq.(4.6)
  if(nz_sum.ge.1) then  ! if nz=0,nbarz=0, we do not need to calculate nonzero M
     do k=0,nz_sum-1 ! loop for M
!        do l=0,n_sum ! loop for j
        do l=0,n_sum-(k+1) ! loop for j
           if(k.eq.0) then ! when M=0, there is no M=-1
              table_Rj(0,0,k+1,l) = c*table_Rj(0,0,k,l+1)
           else
              table_Rj(0,0,k+1,l) = c*table_Rj(0,0,k,l+1) +k*table_Rj(0,0,k-1,l+1)
           end if
!           write(*,'(4i5,1es16.6)') 0,0,k+1,l,table_Rj(0,0,k+1,l)
        end do
     end do
  end if

!  write(*,*)

  ! Generate nonzero L for N=0 and all M using Eq.(4.7)
  do k=0,nz_sum ! loop for all M
     if(ny_sum.ge.1) then ! if ny=0,nbary=0, we do not need to calculate nonzero L
        do j=0,ny_sum-1 ! loop for L
!           do l=0,n_sum ! loop for j
           do l=0,n_sum-(j+1+k) ! loop for j
              if(j.eq.0) then ! L=0
                 table_Rj(0,j+1,k,l) = b*table_Rj(0,j,k,l+1) 
              else
                 table_Rj(0,j+1,k,l) = b*table_Rj(0,j,k,l+1) +j*table_Rj(0,j-1,k,l+1)
              end if
!              write(*,'(4i5,1es16.6)') 0,j+1,k,l,table_Rj(0,j+1,k,l)
           end do
        end do
     end if
  end do
  
!  write(*,*)
!  stop

  ! Generate nonzero N for all L&M using Eq.(4.8)
  do k=0,nz_sum ! loop for all M
     do j=0,ny_sum ! loop for all L
        if(nx_sum.ge.1) then ! if nx=0,nbarx=0, we do not need to calculate nonzero N
           do i=0,nx_sum-1 ! loop for N
!              do l=0,n_sum ! loop for j
              do l=0,n_sum-(i+1+j+k) ! loop for j
                 if(i.eq.0) then
                    table_Rj(i+1,j,k,l) = a*table_Rj(i,j,k,l+1)
                 else
                    table_Rj(i+1,j,k,l) = a*table_Rj(i,j,k,l+1) +i*table_Rj(i-1,j,k,l+1)
                 end if
 !                write(*,'(4i5,1es16.6)') i+1,j,k,l,table_Rj(i+1,j,k,l)
              end do
           end do
        end if
     end do
  end do

  ! R_NLM = R_NLM0
  do i=0,nx_sum ! loop for all N
     do j=0,ny_sum ! loop for all L
        do k=0,nz_sum ! loop for all M
!           table_R(i,j,k) = table_Rj(i,j,k,0)
           R(i,j,k) = table_Rj(i,j,k,0)
!           write(*,'(3i5,2es16.6)') i,j,k,table_Rj(i,j,k,0),table_R(i,j,k)
!           write(*,'(3i5,2es16.6)') i,j,k,table_Rj(i,j,k,0),R(i,j,k)
        end do
     end do
  end do

  table_R(:,:,:) = R(:,:,:)
!  write(*,*)

!!$  do i=0,nx_sum ! loop for all N
!!$     do j=0,ny_sum ! loop for all L
!!$        do k=0,nz_sum ! loop for all M
!!$ !          write(*,'(3i5,2es16.6)') i,j,k,R(i,j,k),table_R(i,j,k)
!!$        end do
!!$     end do
!!$  end do

  return
end subroutine calc_R

!================================================
subroutine calc_PC_and_T(posC,posP,alphaP,vecPC,T)
! calculate P-C and T Eq.(4.2)
!================================================
  implicit none
  real(kind=8),intent(in) :: posC(3),posP(3),alphaP
  real(kind=8),intent(out) :: vecPC(3),T
  integer :: i
  real(kind=8) :: norm2

  norm2 = 0.d0
  do i=1,3
     vecPC(i) = posP(i)-posC(i)
     norm2 = norm2 +vecPC(i)**2
  end do
  T = alphaP*norm2 ! Eq.(4.2)
  return
end subroutine calc_PC_and_T

!====================================
subroutine calc_R000j(alphaP,T,n,list_R000j)
! calculate R000j for j=0 ~ nx+nbarx+ny+nbary+nz+nbarz
! Eq.(4.4)
! modified 120131 by Y. Ikeda
!====================================
  implicit none
  real(kind=8),intent(in) :: alphaP,T
  integer,intent(in) :: n
  real(kind=8),intent(out) :: list_R000j(0:n)
  
  real(kind=8) :: list_FjT(0:n)
  integer :: j
  
!  write(*,*) T
  call calc_FjT(T,n,list_FjT(0:n))

  do j=0,n
 !    list_FjT(j) = FjT(j,T)
     list_R000j(j) = (-2.d0*alphaP)**j *list_FjT(j)
  end do
  return
end subroutine calc_R000j

!====================================
subroutine calc_FjT(T,n,list_FjT)
! calculate FjT for j=0 ~ nx+nbarx+ny+nbary+nz+nbarz
! modified 120202 by Y. Ikeda
!   [1] I. Shavitt, in Methods in Computational Physics,
!       B. Alder, S. Fernbach, and M. Rotenberg, Eds.
!       (Academic, New York, 1963), Vol. 2, pp. 1-45.
!   [2] F. E. Harris, Int. J. Quant. Chem. 23 (1983) 1469
!   [3] B. A. Mamedov J. Math. Chem. 36 (2004) 301
!   are referred
!====================================
  implicit none
  real(kind=8),intent(in) :: T
  integer,intent(in) :: n
  real(kind=8),intent(out) :: list_FjT(0:n)
  
  real(kind=8) :: Factor, FnT_old
  integer :: j, k
  real(kind=8), parameter :: PI  = 3.14159265358979323846264338327950288419716939937510d0
  real(kind=8), parameter :: EPS = 1.d-18

! for small T, only the first term is used
  if (T <= 1.d-8) then
    do j = 0, n
      list_FjT(j) = exp(-T) / dble(2*j+1)
    end do

! Eq. (9) in [2] is used for the intermediate region
  else if (T <= 35.d0) then
  ! for j = 0, use the error function representation
    list_FjT(0) = 0.5d0 * sqrt(pi / T) * erf(sqrt(T))
    if (n >= 1) then
      Factor = 1.d0 / dble(2*n+1)
      FnT_old = Factor
      k = 0
      ! recurrence start
      do
        k = k + 1
        Factor = Factor * 2.d0 * T / dble(2*(n+k)+1)
        list_FjT(n) = FnT_old + Factor
        if (abs(list_FjT(n)-FnT_old) <= EPS) exit
        if (k >= 100) stop "Error: subroutine Calc_FjT: FjT not converged"
        FnT_old = list_FjT(n)
      end do
      list_FjT(n) = list_FjT(n) * exp(-T)
      ! Eq. (6) in [2] is used, downward reccurence
      do j = n-1, 1, -1
        list_FjT(j) = ((2.d0 * T) * list_FjT(j+1) + exp(-T)) / dble(2*j+1)
      end do
    end if

  else ! T > 35.d0
  ! First, complimentary functions (Eq. (26) in [1]) are contained in list_FjT
    list_FjT(0) = 0.5d0 * sqrt(pi / T) * erfc(sqrt(T))
    ! upward recurrence for complimentary functions (Eq. (28) in [1])
    do j = 1, n
      list_FjT(j) = (dble(2*j-1) * list_FjT(j-1) + exp(-T)) / (2.d0 * T)
    end do
    ! now, FjT are obtained by Eq. (27) in [1]
    Factor = 0.5d0 * sqrt(pi / T)
    list_FjT(0) = Factor - list_FjT(0)
    do j = 1, n
      Factor = Factor * (dble(j)-0.5d0) / T
      list_FjT(j) = Factor - list_FjT(j)
    end do
  end if

  return

end subroutine

!========================================================================
subroutine calc_d(alphaP,PAx,PBx,n,nbar,d)
! calculate coefficients for Lambda polynomical expansion
! (p.220 of Eq.(2.17)~(2.22) McMurchie & Davidson)
!========================================================================
  implicit none
  
  real(kind=8), intent(in) :: alphaP,PAx,PBx
  integer, intent(in) :: n,nbar
  real(kind=8), intent(out) :: d(0:n,0:nbar,0:n+nbar)

!  real(kind=8) :: PAx,PBx
  integer :: i,j,k

!  PAx = vecPA(1); PBx = vecPB(1)

  d = 0.d0
  d(0,0,0) = 1.d0 ! (2.22) McMurchie & Davidson

  ! Make table for n=0 (loop over nbar & N) 
  do j=1,nbar 
     do k=0,j  ! loop for N 
        if(k==0) then ! left colum 
           d(0,j,k) = PBx*d(0,j-1,k) +(k+1)*d(0,j-1,k+1)
        elseif(k==j) then ! right side
           d(0,j,k) = d(0,j-1,k-1)/2.d0/alphaP
        else
           d(0,j,k) = d(0,j-1,k-1)/2.d0/alphaP +PBx*d(0,j-1,k) +(k+1)*d(0,j-1,k+1)
        end if
     end do
  end do
  
  ! loop over n & N starting from n=0 table (for every nbar)
  do j=0,nbar
     do i=1,n 
        do k=0,j+i
           if(k==0) then
              d(i,j,k) = PAx*d(i-1,j,k) +(k+1)*d(i-1,j,k+1)
           elseif(k==(j+i)) then
              d(i,j,k) = d(i-1,j,k-1)/2.d0/alphaP
           else
              d(i,j,k) = d(i-1,j,k-1)/2.d0/alphaP +PAx*d(i-1,j,k) +(k+1)*d(i-1,j,k+1)
           end if
        end do
     end do
  end do
      
  return
end subroutine calc_d
  

!========================================================================
subroutine calc_gaussP(posA,alphaA,posB,alphaB,c_PAB,posP,alphaP,vecPA,vecPB)
! modified 120301 by Y. Ikeda
!========================================================================
  implicit none
  
  real(kind=8), intent(in) :: posA(3),alphaA,posB(3),alphaB
  real(kind=8), intent(out) :: c_PAB,posP(3),alphaP,vecPA(3),vecPB(3)
  real(kind=8) :: distAB2
  integer :: i

  distAB2 = 0.d0
  alphaP = alphaA + alphaB
  do i=1,3
!     posP(i) = (alphaA*posA(i)+alphaB*posB(i))/alphaP
     posP(i) = posA(i)+alphaB/alphaP*(posB(i)-posA(i))
     vecPA(i) = posP(i) - posA(i)
     vecPB(i) = posP(i) - posB(i)
     distAB2 = distAB2 +(posA(i)-posB(i))**2
  end do

  c_PAB = exp(-alphaA*alphaB/alphaP *distAB2)
  return
end subroutine calc_gaussP
  

