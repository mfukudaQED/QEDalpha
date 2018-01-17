!=============================================================================
!
!  2011.5.30
! -subroutines to set large (quadruple) matrix (4*NBS x 4*NBS)
!  4 = 2 (+/-) x 2 (Kramers pair)
!  ordering of index may be reconsidered.
!  index_from_Qmat : utility routine for converting indice
!  setQmat_h
!  density_Qmat  : position probability density
!
!  6.3
!  setQmat_calFj0
!  setQmat_twoele
!  
!  2011.8.28 edited by fukuda
!  zeta_Qmat    :zeta force density
!  t_Qmat       :spintorque density
!  tau_Qmat     :stress tensor density
!  divj_Qmat    :div of current density
!  j_Qmat       :current density
!  s_Qmat　　   :spin density
!
!  10.31
!  rho_Qmat  : charge density
!
!  11.4
!  TM_Qmat : kinetic+mass energy
!
! 1 ~ 2*NBS : "+"
! 2*NBS +1 ~ 2*NBS : "-"
! in "+" and "-", 1,1bar,2,2bar, ... respectively

!  2012.5.20
!  modified setQmat_twoele to use faster routines.
!
!  11.1
!  subroutine setQmat_E(posR,E_Qmat)
!
!  11.16
!  subroutine setQmat_intJcA(posR,vecJ,intJcA_Qmat)
!
!  2013.3.13
!  subroutine setQmat_twoele_nz(twoele_Qmat_nz)
!
!=============================================================================


!====================================================================
! densities
!   zeta_Qmat(i,x,y,z,nn,mm)
!   t_Qmat(i,x,y,z,nn,mm)
!   tau_Qmat(k,l,x,y,z,nn,mm)
!   divj_Qmat(x,y,z,nn,mm)
!   j_Qmat(i,x,y,z,nn,mm)
!   s_Qmat(i,x,y,z,nn,mm)
!   rho_Qmat(x,y,z,nn,mm)
!   density_Qmat(x,y,z,nn,mm)
!====================================================================

!============================================================
function zeta_Qmat(i,x,y,z,nn,mm)
! nn,mm : labels for molecular orbitals (KP) and +/-. 1~4*NBS
!============================================================
  use Precision
  use DiracOutput
  implicit none

  complex(kind=dp) :: zeta_Qmat
  integer,intent(in) :: i
  real(kind=dp),intent(in) :: x,y,z
  integer,intent(in) :: nn,mm

  complex(kind=dp) :: zeta_mat
  integer :: n,m ! index for 1~2*NBS (+ or - specified)
  character(LEN=1) :: a,b
  
  call index_from_Qmat(nn,n,a)
  call index_from_Qmat(mm,m,b)

  zeta_Qmat = zeta_mat(i,x,y,z,n,a,m,b)

  return
end function zeta_Qmat

!============================================================
function t_Qmat(i,x,y,z,nn,mm)
! nn,mm : labels for molecular orbitals (KP) and +/-. 1~4*NBS
!============================================================
  use Precision
  use DiracOutput
  implicit none

  complex(kind=dp) :: t_Qmat
  integer,intent(in) :: i
  real(kind=dp),intent(in) :: x,y,z
  integer,intent(in) :: nn,mm

  complex(kind=dp) :: t_mat
  integer :: n,m ! index for 1~2*NBS (+ or - specified)
  character(LEN=1) :: a,b
  
  call index_from_Qmat(nn,n,a)
  call index_from_Qmat(mm,m,b)

  t_Qmat = t_mat(i,x,y,z,n,a,m,b)

  return
end function t_Qmat

!============================================================
function tau_Qmat(k,l,x,y,z,nn,mm)
! nn,mm : labels for molecular orbitals (KP) and +/-. 1~4*NBS
!============================================================
  use Precision
  use DiracOutput
  implicit none

  complex(kind=dp) :: tau_Qmat
  integer,intent(in) :: k,l
  real(kind=dp),intent(in) :: x,y,z
  integer,intent(in) :: nn,mm

  complex(kind=dp) :: tau_mat
  integer :: n,m ! index for 1~2*NBS (+ or - specified)
  character(LEN=1) :: a,b
  
  call index_from_Qmat(nn,n,a)
  call index_from_Qmat(mm,m,b)

  tau_Qmat = tau_mat(k,l,x,y,z,n,a,m,b)

  return
end function tau_Qmat

!============================================================
function divj_Qmat(x,y,z,nn,mm)
! nn,mm : labels for molecular orbitals (KP) and +/-. 1~4*NBS
!============================================================
  use Precision
  use DiracOutput
  implicit none

  complex(kind=dp) :: divj_Qmat
  real(kind=dp),intent(in) :: x,y,z
  integer,intent(in) :: nn,mm

  complex(kind=dp) :: divj_mat
  integer :: n,m ! index for 1~2*NBS (+ or - specified)
  character(LEN=1) :: a,b
  
  call index_from_Qmat(nn,n,a)
  call index_from_Qmat(mm,m,b)

  divj_Qmat = divj_mat(x,y,z,n,a,m,b)

  return
end function divj_Qmat

!============================================================
function j_Qmat(i,x,y,z,nn,mm)
! nn,mm : labels for molecular orbitals (KP) and +/-. 1~4*NBS
!============================================================
  use Precision
  use DiracOutput
  implicit none

  complex(kind=dp) :: j_Qmat
  integer,intent(in) :: i
  real(kind=dp),intent(in) :: x,y,z
  integer,intent(in) :: nn,mm

  complex(kind=dp) :: j_mat
  integer :: n,m ! index for 1~2*NBS (+ or - specified)
  character(LEN=1) :: a,b
  
  call index_from_Qmat(nn,n,a)
  call index_from_Qmat(mm,m,b)

  j_Qmat = j_mat(i,x,y,z,n,a,m,b)

  return
end function j_Qmat

!============================================================
function s_Qmat(i,x,y,z,nn,mm)
! nn,mm : labels for molecular orbitals (KP) and +/-. 1~4*NBS
!============================================================
  use Precision
  use DiracOutput
  implicit none

  complex(kind=dp) :: s_Qmat
  integer,intent(in) :: i
  real(kind=dp),intent(in) :: x,y,z
  integer,intent(in) :: nn,mm

  complex(kind=dp) :: s_mat
  integer :: n,m ! index for 1~2*NBS (+ or - specified)
  character(LEN=1) :: a,b
  
  call index_from_Qmat(nn,n,a)
  call index_from_Qmat(mm,m,b)

  s_Qmat = s_mat(i,x,y,z,n,a,m,b)

  return
end function s_Qmat

!============================================================
function rho_Qmat(x,y,z,nn,mm)
! Nn,mm : labels for molecular orbitals (KP) and +/-. 1~4*NBS
!============================================================
  use Precision
  use DiracOutput
  implicit none

  complex(kind=dp) :: rho_Qmat
  real(kind=dp),intent(in) :: x,y,z
  integer,intent(in) :: nn,mm

  complex(kind=dp) :: rho_mat
  integer :: n,m ! index for 1~2*NBS (+ or - specified)
  character(LEN=1) :: a,b
  
  call index_from_Qmat(nn,n,a)
  call index_from_Qmat(mm,m,b)

  rho_Qmat = rho_mat(x,y,z,n,a,m,b)

  return
end function rho_Qmat

!============================================================
function density_Qmat(x,y,z,nn,mm)
! nn,mm : labels for molecular orbitals (KP) and +/-. 1~4*NBS
!============================================================
  use Precision
  use DiracOutput
  implicit none

  complex(kind=dp) :: density_Qmat
  real(kind=dp),intent(in) :: x,y,z
  integer,intent(in) :: nn,mm

  complex(kind=dp) :: density_mat
  integer :: n,m ! index for 1~2*NBS (+ or - specified)
  character(LEN=1) :: a,b
  
  call index_from_Qmat(nn,n,a)
  call index_from_Qmat(mm,m,b)

  density_Qmat = density_mat(x,y,z,n,a,m,b)

  return
end function density_Qmat


!====================================================================
! integrals
!   setQmat_intret(NTIME,calJi_Qmat,calLi_Qmat)
!   subroutine setQmat_E(posR,E_Qmat)  [121101]
!   setQmat_twoele(twoele_Qmat)
!   setQmat_calFj0(calFj0_Qmat)
!   setQmat_TM(TM_Qmat)
!   setQmat_h(h_Qmat)
!====================================================================





!!$!======================================================
!!$subroutine setQmat_intret(NTIME,calJi_Qmat,calLi_Qmat)
!!$! Calj^i_NMPQ
!!$! calL^i_NMPQ
!!$!======================================================
!!$  use Precision
!!$  use DiracOutput
!!$  use Constants  ! for CCC
!!$!  use IntegralStorage
!!$  implicit none
!!$
!!$  integer,intent(in) :: NTIME
!!$
!!$  complex(kind=dp),intent(out) :: calJi_Qmat(4*NBS,4*NBS,4*NBS,4*NBS,NTIME)
!!$  complex(kind=dp),intent(out) :: calLi_Qmat(4*NBS,4*NBS,4*NBS,4*NBS,NTIME)
!!$  
!!$  integer :: nn,mm,pp,qq  ! index for 1~4*NBS
!!$
!!$  complex(kind=dp) :: intJret_mat, intLret_mat
!!$  integer :: n,m,p,q
!!$  character(LEN=1) :: a,b,c,d
!!$  real(kind=dp) :: R 
!!$  real(kind=dp) :: rmin(3),rmax(3)
!!$  real(kind=dp) :: r_box
!!$
!!$  integer :: i,j,k,l
!!$
!!$  r_box = 2._dp
!!$  rmin(1) = -r_box;   rmin(2) = -r_box;   rmin(3) = -r_box
!!$  rmax(1) =  r_box;   rmax(2) =  r_box;   rmax(3) =  r_box
!!$!  rmin(1) = -3._dp;   rmin(2) = -3._dp;   rmin(3) = -3._dp
!!$!  rmax(1) =  3._dp;   rmax(2) =  3._dp;   rmax(3) =  3._dp
!!$
!!$  do j=1,NTIME
!!$     R = CCC*j*DeltaT
!!$
!!$     do nn=1,4*NBS
!!$        do mm=1,4*NBS
!!$           do pp=1,4*NBS
!!$              do qq=1,4*NBS
!!$                 
!!$                 call index_from_Qmat(nn,n,a)
!!$                 call index_from_Qmat(mm,m,b)
!!$                 call index_from_Qmat(pp,p,c)
!!$                 call index_from_Qmat(qq,q,d)
!!$                 
!!$                 write(13,"(4es16.6)") intJret_mat(n,a,m,b,p,c,q,d,R,rmin,rmax),intLret_mat(n,a,m,b,p,c,q,d,R,rmin,rmax)  ! for storage
!!$!                 write(*,"(1i6,4a3,4i6,4es16.6)") j,a,b,c,d,n,m,p,q, &
!!$!                      & intJret_mat(n,a,m,b,p,c,q,d,R,rmin,rmax),intLret_mat(n,a,m,b,p,c,q,d,R,rmin,rmax)  ! for storage
!!$                 calJi_Qmat(nn,mm,pp,qq,j) = intJret_mat(n,a,m,b,p,c,q,d,R,rmin,rmax)
!!$                 calLi_Qmat(nn,mm,pp,qq,j) = intLret_mat(n,a,m,b,p,c,q,d,R,rmin,rmax)
!!$              end do
!!$           end do
!!$        end do
!!$     end do
!!$
!!$  end do
!!$
!!$  return
!!$end subroutine setQmat_intret

!======================================================
subroutine setQmat_intI2JcA(vecJ,intI2JcA_Qmat)
! 1 ~ 2*NBS : "+"
! 2*NBS +1 ~ 2*NBS : "-"
! in "+" and "-", 1,1bar,2,2bar, ... respectively
!
! Integration of intJcA over M. With a factor of (-Ze e)/c
!
! 121116
!======================================================
  use Precision
  use DiracOutput
  use Constants
  use System_Medium
  implicit none
  
  real(kind=dp),intent(in) :: vecJ(3)
  complex(kind=dp),intent(out) :: intI2JcA_Qmat(4*NBS,4*NBS) ! (-Ze e)/c *int_M dr IJcA_NM(vecR,vecJ)
  
  integer :: nn,mm  ! index for 1~4*NBS
  
  integer,parameter :: Ng=10 ! grid number
  
  integer :: i,j,k
  real(kind=dp) :: dx,dy,dz,posR(3)
  complex(kind=dp) :: sum(4*NBS,4*NBS)
  complex(kind=dp) :: intJcA_Qmat(4*NBS,4*NBS) ! IJcA_NM(vecR,vecJ)

  L_Med = 10._dp
  R_Sys = 1._dp

  dx = L_Med/Ng
  dy = L_Med/Ng
  dz = L_Med/Ng

  sum(:,:) = (0._dp,0._dp)
  do i = -Ng,Ng
     do j = -Ng,Ng
        do k = -Ng,Ng
           posR(1) = i*dx  
           posR(2) = j*dy  
           posR(3) = k*dz  
           call setQmat_intJcA(posR,vecJ,intJcA_Qmat)
           do nn=1,4*NBS
              do mm=1,4*NBS
                 sum(nn,mm) = sum(nn,mm) +intJcA_Qmat(nn,mm)*dx*dy*dz
              end do
           end do
        end do
     end do
  end do
  intI2JcA_Qmat(:,:) = (-Ze/CCC)*sum(:,:)
  
end subroutine setQmat_intI2JcA

!======================================================
subroutine setQmat_intJcA(posR,vecJ,intJcA_Qmat)
! 1 ~ 2*NBS : "+"
! 2*NBS +1 ~ 2*NBS : "-"
! in "+" and "-", 1,1bar,2,2bar, ... respectively
! 121115
!  11.16 modified to include vecJ
!======================================================
  use Precision
  use DiracOutput
  implicit none

  real(kind=dp),intent(in) :: posR(3),vecJ(3)
  complex(kind=dp),intent(out) :: intJcA_Qmat(4*NBS,4*NBS) ! IJcA_NM(vecR,vecJ)

  complex(kind=dp) :: intJcA(-NBS:NBS,2,-NBS:NBS,2)
  Integer :: nn,mm  ! index for 1~4*NBS
  integer :: n,m ! index for 1~2*NBS (+ or - specified)
  character(LEN=1) :: a,b ! + or -
  integer :: aa,bb ! 1 -> +, 2 -> -

  call calc_intJcA_mat(posR,vecJ,intJcA)
  
  do nn=1,4*NBS
     do mm=1,4*NBS
        
        call index_from_Qmat(nn,n,a)
        call index_from_Qmat(mm,m,b)
        
        ! "+" --> 1, "-" --> 
        if(a.eq."+") aa = 1
        if(a.eq."-") aa = 2
        if(b.eq."+") bb = 1
        if(b.eq."-") bb = 2
        
        intJcA_Qmat(nn,mm) = intJcA(n,aa,m,bb) 
        
     end do
  end do

end subroutine setQmat_intJcA

!======================================================
subroutine setQmat_E(posR,E_Qmat)
! 1 ~ 2*NBS : "+"
! 2*NBS +1 ~ 2*NBS : "-"
! in "+" and "-", 1,1bar,2,2bar, ... respectively
! 121031
!======================================================
  use Precision
  use DiracOutput
  implicit none

  real(kind=dp),intent(in) :: posR(3)  
  complex(kind=dp),intent(out) :: E_Qmat(3,4*NBS,4*NBS) ! E^k_NM

  complex(kind=dp) :: intE_mat(3)  ! --> used for output of slow subroutine
  complex(kind=dp) :: intE(3,-NBS:NBS,2,-NBS:NBS,2)
  Integer :: nn,mm  ! index for 1~4*NBS
  integer :: n,m ! index for 1~2*NBS (+ or - specified)
  character(LEN=1) :: a,b ! + or -
  integer :: aa,bb ! 1 -> +, 2 -> -
  integer :: kk ! vector index

  !--- faster routine but need large memory ---
  call calc_intEk_mat(posR,intE) 
  !------------------------------------------

  do kk=1,3
     do nn=1,4*NBS
        do mm=1,4*NBS
           
           call index_from_Qmat(nn,n,a)
           call index_from_Qmat(mm,m,b)
           
           !--- faster routine but need large memory ---
           ! "+" --> 1, "-" --> 2 to be used in modified routines 
           if(a.eq."+") aa = 1
           if(a.eq."-") aa = 2
           if(b.eq."+") bb = 1
           if(b.eq."-") bb = 2
           !------------------------------------------

!!$           !--- slow subroutine previously used ---
!!$           call calc_intE_mat(posR,n,a,m,b,intE_mat)
!!$           E_Qmat(kk,nn,mm) = intE_mat(kk)
!!$           !------------------------------------------
           
           !--- faster routine but need large memory ---
           E_Qmat(kk,nn,mm) = intE(kk,n,aa,m,bb) 
           !------------------------------------------
           
        end do
     end do
  end do
  write(*,*) " Integrals E_Qmat done."

end subroutine setQmat_E


!======================================================
subroutine setQmat_twoele(twoele_Qmat)
! ( NM | PQ )
!
! 2012.5.19
! modified to use faster routines for twoele
!======================================================
  use Precision
  use DiracOutput
  use IntegralStorage
  implicit none
  
  complex(kind=dp),intent(out) :: twoele_Qmat(4*NBS,4*NBS,4*NBS,4*NBS)
  
  complex(kind=dp) :: inttwoele_mat ! --> very slow function inttwoele_mat(n,a,m,b,p,c,q,d)
  complex(kind=dp) :: inttwoele(-NBS:NBS,2,-NBS:NBS,2,-NBS:NBS,2,-NBS:NBS,2)
  Integer :: nn,mm,pp,qq  ! index for 1~4*NBS
  integer :: n,m,p,q ! index for 1~2*NBS (+ or - specified)
  character(LEN=1) :: a,b,c,d ! + or -
  integer :: aa,bb,cc,dd ! 1 -> +, 2 -> -

  character(LEN=80) :: temp
  real(kind=dp) :: int_real,int_complex
  character(LEN=300) :: readfile

  if(there_is_twoele) then
     readfile = trim(FILESFOLDER)//"/"//file_twoele
     write(*,*) " There is twoele file. Read integrals from ", trim(readfile)
!     open(unit=200,file=file_twoele,status='unknown',form='formatted')
     open(unit=200,file=readfile,status='unknown',form='formatted')

     do nn=1,4*NBS
        do mm=1,4*NBS
           do pp=1,4*NBS
              do qq=1,4*NBS

!                 read(200,*) temp,temp,temp,temp, temp,temp,temp,temp, int_real,int_complex            
                 ! use same format as storage (see below)
                 read(200,*) int_real,int_complex            
                 twoele_Qmat(nn,mm,pp,qq) = cmplx(int_real,int_complex,dp)

              end do
           end do
        end do
     end do
     write(*,*) " Reading done."
     
  else ! if there isn't precomputed integrals, compute integrals.
     write(*,*) " There is no twoele file. Compute integrals and store at fort.12."

     !--- faster routine but need large memory ---
     call calc_inttwoele_mat(inttwoele) 
     !------------------------------------------

     do nn=1,4*NBS
        do mm=1,4*NBS
           do pp=1,4*NBS
              do qq=1,4*NBS
                 
                 call index_from_Qmat(nn,n,a)
                 call index_from_Qmat(mm,m,b)
                 call index_from_Qmat(pp,p,c)
                 call index_from_Qmat(qq,q,d)

                 !--- faster routine but need large memory ---
                 ! "+" --> 1, "-" --> 2 to be used in modified routines 
                 if(a.eq."+") aa = 1
                 if(a.eq."-") aa = 2
                 if(b.eq."+") bb = 1
                 if(b.eq."-") bb = 2
                 if(c.eq."+") cc = 1
                 if(c.eq."-") cc = 2
                 if(d.eq."+") dd = 1
                 if(d.eq."-") dd = 2
                 !------------------------------------------

!!$                 !--- very slow function previously used ---
!!$!!                 write(12,"(4a3,4i6,2es16.6)") a,b,c,d,n,m,p,q,inttwoele_mat(n,a,m,b,p,c,q,d) 
!!$                 write(12,"(2es16.6)") inttwoele_mat(n,a,m,b,p,c,q,d)  ! for storage
!!$                 twoele_Qmat(nn,mm,pp,qq) = inttwoele_mat(n,a,m,b,p,c,q,d) 
!!$                 !------------------------------------------

                 !--- faster routine but need large memory ---
                 write(12,"(2es16.6)") inttwoele(n,aa,m,bb,p,cc,q,dd)  ! for storage
                 twoele_Qmat(nn,mm,pp,qq) = inttwoele(n,aa,m,bb,p,cc,q,dd) 
                 !------------------------------------------

              end do
           end do
        end do
     end do
     write(*,*) " Integrals twoele done."
     
  end if

  return
end subroutine setQmat_twoele



!======================================================
subroutine setQmat_twoele_nz(twoele_Qmat_nz)
! ( NM | PQ )
!
! 2013.3.13
! modified setQmat_twoele to store only small values of twoele
!======================================================
  use Precision
  use DefineTypes
  use DiracOutput
  use IntegralStorage
  implicit none
  
  type(nonzeromat4legs),allocatable,intent(out) :: twoele_Qmat_nz(:)
  
  complex(kind=dp) :: inttwoele_mat ! --> very slow function inttwoele_mat(n,a,m,b,p,c,q,d)
  complex(kind=dp) :: inttwoele(-NBS:NBS,2,-NBS:NBS,2,-NBS:NBS,2,-NBS:NBS,2)
  Integer :: nn,mm,pp,qq  ! index for 1~4*NBS
  integer :: n,m,p,q ! index for 1~2*NBS (+ or - specified)
  character(LEN=1) :: a,b,c,d ! + or -
  integer :: aa,bb,cc,dd ! 1 -> +, 2 -> -

  character(LEN=80) :: temp
  real(kind=dp) :: int_real,int_complex
  real(kind=dp) :: th
  character(LEN=300) :: readfile

  integer :: i
  type(nonzeromat4legs),allocatable :: nonzero_twoele(:) ! need large memory but only used when twoele is computed.

  if(there_is_twoele) then
     readfile = trim(FILESFOLDER)//"/"//file_twoele
     write(*,*) " There is twoele file. Read integrals from "
     write(*,"(a)") trim(readfile)
     write(10,*) "# There is twoele file. Read integrals from "
     write(10,"(a)") "# "//trim(readfile)
     open(unit=200,file=readfile,status='unknown',form='formatted')
     
     read(200,*) N_twoele, th
     write(*,*) "Number of twoele we use:  ", N_twoele
     write(10,*) "# Number of twoele we use:  ", N_twoele
     write(*,*) "threshold of twoele:  ", th
     write(10,*) "# threshold of twoele:  ", th
     allocate(twoele_Qmat_nz(N_twoele))
     do i=1,N_twoele
        read(200,*) twoele_Qmat_nz(i)%a,twoele_Qmat_nz(i)%b,twoele_Qmat_nz(i)%c, &
             & twoele_Qmat_nz(i)%d,int_real,int_complex 
        twoele_Qmat_nz(i)%val = cmplx(int_real,int_complex,dp)
     end do
     
     write(*,*) " Reading done."
     
  else ! if there isn't precomputed integrals, compute integrals.
     write(*,*) " There is no twoele file. Compute integrals and store at fort.12."

     !--- faster routine but need large memory ---
     call calc_inttwoele_mat(inttwoele) 
     !------------------------------------------

     allocate(nonzero_twoele((4*NBS)**4))
     i = 0
     do nn=1,4*NBS
        do mm=1,4*NBS
           do pp=1,4*NBS
              do qq=1,4*NBS
                 
                 call index_from_Qmat(nn,n,a)
                 call index_from_Qmat(mm,m,b)
                 call index_from_Qmat(pp,p,c)
                 call index_from_Qmat(qq,q,d)

                 ! "+" --> 1, "-" --> 2 to be used in modified routines 
                 if(a.eq."+") aa = 1
                 if(a.eq."-") aa = 2
                 if(b.eq."+") bb = 1
                 if(b.eq."-") bb = 2
                 if(c.eq."+") cc = 1
                 if(c.eq."-") cc = 2
                 if(d.eq."+") dd = 1
                 if(d.eq."-") dd = 2

                 if(abs(inttwoele(n,aa,m,bb,p,cc,q,dd)).gt.TH_twoele) then
                    i = i+1
                    nonzero_twoele(i)%a = nn
                    nonzero_twoele(i)%b = mm
                    nonzero_twoele(i)%c = pp
                    nonzero_twoele(i)%d = qq
                    nonzero_twoele(i)%val = inttwoele(n,aa,m,bb,p,cc,q,dd)
                 end if

              end do
           end do
        end do
     end do
     write(*,*) " Integrals twoele done."
     N_twoele = i

     write(*,*) "Number of twoele we use:  ", N_twoele
     write(10,*) "# Number of twoele we use:  ", N_twoele
     write(*,*) "TH_twoele:  ", TH_twoele
     write(10,*) "# TH_twoele:  ", TH_twoele

     write(12,"(1i16,1es14.4)") N_twoele,TH_twoele
     allocate(twoele_Qmat_nz(N_twoele))
     do i=1,N_twoele
        write(12,"(4i6,2es16.6)") nonzero_twoele(i)%a,nonzero_twoele(i)%b,nonzero_twoele(i)%c,&
             & nonzero_twoele(i)%d,nonzero_twoele(i)%val
        twoele_Qmat_nz(i) = nonzero_twoele(i)
     end do

     deallocate(nonzero_twoele)
  end if

  return
end subroutine setQmat_twoele_nz

!======================================================
subroutine setQmat_calFj0(calFj0_Qmat)
! Calf_nmj(t=t0)
!======================================================
  use Precision
  use DiracOutput
  use Constants  ! for Nph
  implicit none

  complex(kind=dp),intent(out) :: calFj0_Qmat(4*NBS,4*NBS,Nph)

  complex(kind=dp) :: intcalFj_mat
  integer :: nn,mm  ! index for 1~4*NBS
  integer :: n,m ! index for 1~2*NBS (+ or - specified)
  character(LEN=1) :: a,b
  integer :: j

  do j=1,Nph
     do nn=1,4*NBS
        do mm=1,4*NBS
           
           call index_from_Qmat(nn,n,a)
           call index_from_Qmat(mm,m,b)
           
!           write(*,"(1i6,2a3,2i6,2es16.6)") j,a,b,n,m,intcalFj_mat(0._dp,j,n,a,m,b)
           calFj0_Qmat(nn,mm,j) = intcalFj_mat(0._dp,j,n,a,m,b)  ! calFj at t=0
           
        end do
     end do
  end do

end subroutine setQmat_calFj0


!======================================================
subroutine setQmat_TM(TM_Qmat)
! 1 ~ 2*NBS : "+"
! 2*NBS +1 ~ 2*NBS : "-"
! in "+" and "-", 1,1bar,2,2bar, ... respectively
!======================================================
  use Precision
  use DiracOutput
  implicit none

  complex(kind=dp),intent(out) :: TM_Qmat(4*NBS,4*NBS)

  complex(kind=dp) :: intTM_mat
  integer :: nn,mm  ! index for 1~4*NBS
  integer :: n,m ! index for 1~2*NBS (+ or - specified)
  character(LEN=1) :: a,b

  do nn=1,4*NBS
     do mm=1,4*NBS
        
        call index_from_Qmat(nn,n,a)
        call index_from_Qmat(mm,m,b)

        TM_Qmat(nn,mm) = intTM_mat(n,a,m,b)
        
     end do
  end do

end subroutine setQmat_TM

!======================================================
subroutine setQmat_h(h_Qmat)
! 1 ~ 2*NBS : "+"
! 2*NBS +1 ~ 2*NBS : "-"
! in "+" and "-", 1,1bar,2,2bar, ... respectively
!======================================================
  use Precision
  use DiracOutput
  implicit none

  complex(kind=dp),intent(out) :: h_Qmat(4*NBS,4*NBS)

  complex(kind=dp) :: inth_mat
  integer :: nn,mm  ! index for 1~4*NBS
  integer :: n,m ! index for 1~2*NBS (+ or - specified)
  character(LEN=1) :: a,b

  do nn=1,4*NBS
     do mm=1,4*NBS
        
!!$        if(nn.le.2*NBS) then
!!$           a = "+"
!!$           n = nn
!!$        else
!!$           a = "-"
!!$           n = nn -2*NBS
!!$        end if
!!$
!!$        if(mm.le.2*NBS) then
!!$           b = "+"
!!$           m = mm
!!$        else
!!$           b = "-"
!!$           m = mm -2*NBS
!!$        end if
!!$
!!$        if(mod(n,2).eq.0) then ! even -> Kramers pair -> negative sign
!!$           n = -n/2
!!$        else ! odd -> not KP
!!$           n = (n+1)/2
!!$        end if
!!$
!!$        if(mod(m,2).eq.0) then ! even -> Kramers pair -> negative sign
!!$           m = -m/2
!!$        else ! odd -> not KP
!!$           m = (m+1)/2
!!$        end if

        call index_from_Qmat(nn,n,a)
        call index_from_Qmat(mm,m,b)

!        write(*,"(2a3,2i6,2es16.6)") a,b,n,m,inth_mat(n,a,m,b)
        h_Qmat(nn,mm) = inth_mat(n,a,m,b)
        
     end do
  end do

end subroutine setQmat_h

!====================================================================
! utility routines
!   index_from_Qmat(nn,n,a)
!====================================================================


!======================================================
subroutine index_from_Qmat(nn,n,a)
! utility routine to convert indices to n^a (inc. KP or not)
! from N 
!======================================================
  use Precision
  use DiracOutput
  implicit none
  
  integer,intent(in) :: nn
  integer,intent(out) :: n
  character(LEN=1),intent(out) :: a
  
  if(nn.le.2*NBS) then
     a = "+"
     n = nn
  else
     a = "-"
     n = nn -2*NBS
  end if

  if(mod(n,2).eq.0) then ! even -> Kramers pair -> negative sign
     n = -n/2
  else ! odd -> not KP
     n = (n+1)/2
  end if

  return
end subroutine index_from_Qmat
  
