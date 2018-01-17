! Last Change:22-May-2012.
!============================================================================
!  subroutine read_DiracOutput 2012/04/16  written by Fukuda
!  subroutine read_qedinp 2012/04/16  written by Fukuda
!  subroutine read_spininp 2012/04/16  written by Fukuda
!============================================================================
!========================================================================
subroutine read_DiracOutput
!========================================================================

  use DiracOutput

  implicit none

  character(len=24) today,version
  character(len=6) SETDC2
  character(len=5) ircop
  integer bnum_L(100),bnum_S(100),i,j,k
  integer bnum_LT,bnum_ST
  integer bnum_LTA(100),bnum_STA(100)
  integer,allocatable:: n(:,:)
  integer :: NBS_Eg, NBS_Pg, NBS_Eu, NBS_Pu!number of molecular Electron and Positron orbitals
  integer :: NBS_TOTALg, NBS_TOTALu !NBS_E + NBS_P
  integer :: NBS_Lg, NBS_Sg, NBS_Lu, NBS_Su ! number of primitive gaussians for large and small components
  integer,allocatable:: cont(:,:)
  double precision, allocatable :: cpg_L(:),cpg_S(:)

  version = '    2012/4/16'
!this program is only for contract!
  write(*,*) 'read_DiracOutput ver. ', version
  call FDATE(today)
  write(*,*) 'Calculation date : ',today

!  call read_qedinp
!  call read_spininp

!!$!-------------read input file (qed.inp)---------------------------------
!!$!  read (*,'(a12)') FIFI1 ! basis.txt of Dirac10
!!$!  read (*,'(a12)') FIFI2 ! vectors.txt of Dirac10
!!$  read (*,'(a)') FILESFOLDER
!!$  FIFI1 = trim(FILESFOLDER)//"/"//"basis.txt"
!!$  FIFI2 = trim(FILESFOLDER)//"/"//"vectors.txt"
!!$  write(*,'(a)') FIFI1
!!$  write(*,'(a)') FIFI2
!!$  read (*,'(a)') FIFI3 !Dirac output file of Dirac10
!!$  FIFI3 = trim(FILESFOLDER)//"/"//FIFI3
!!$  write(*,'(a)') FIFI3
!!$  read (*,*) symm ! 0=>atom, 1=>molecule using symmmetry calc, 2=>molecule not using symmetry calc
!!$  read (*,*) kbal
!!$!  stop

!----------search Basis and Orbitals from DIRAC10 output file-----------------------------
  
  open(unit=43,file=FIFI3)
  do
   read (43,'(44x,a6)') SETDC2
    if (SETDC2.eq."SETDC2") exit
  end do
  do i=1,4
   read (43,'(44x,a6)') SETDC2
  end do

   if(SYMM.eq.2) then 
    read(43,'(a5,6x,5x,2i5,5x,3i5)') ircop, NBS_L, NBS_S, NBS_TOTAL, NBS_E, NBS_P 
    ! NBS_L and NBS_S are revised later at countBasis
    NBS = NBS_E !NBS_E = NBS_P if KBAL
    if(kbal.eq.1) NBS_P = NBS_E
   else !symmetry
    read(43,'(a5,6x,5x,2i5,5x,3i5)') ircop, NBS_Lg, NBS_Sg, NBS_TOTALg, NBS_Eg, NBS_Pg 
    read(43,'(a5,6x,5x,2i5,5x,3i5)') ircop, NBS_Lu, NBS_Su, NBS_TOTALu, NBS_Eu, NBS_Pu 
    NBS_L = NBS_Lg + NBS_Lu
    NBS_S = NBS_Sg + NBS_Su
    NBS_TOTAL = NBS_TOTALg + NBS_TOTALu
    NBS_E = NBS_Eg + NBS_Eu
    NBS_P = NBS_Pg + NBS_Pu
    NBS = NBS_E !NBS_E = NBS_P if KBAL
    if(kbal.eq.1) NBS_P = NBS_E
   end if
!  read (*,*) NBS,NBS_L,NBS_S !NBS is total number of molecular orbitals
!  read (*,*) NEL
  close(unit=43)
  write(*,*)'# NBS :',NBS
  write(*,*)'# NBS_E :',NBS_E
  write(*,*)'# NBS_P :',NBS_P
  write(*,*)'# NBS_TOTAL :',NBS_TOTAL
!-------------finish reading inputfile----------------------------------

!  NBS0 = NBS_L + NBS_S

  write(*,*) 'count primitive gaussian in Basis.txt'
  open(unit=11,file=FIFI1)
  call countBasis(bnum_L,bnum_S,bnum_LT,bnum_ST,bnum_LTA,bnum_STA)
  close(unit=11)

   write(*,*)'NAT=',NAT
   do i=1,NAT
    write(*,*)'bnum_L(i)=',bnum_L(i),'bnum_S(i)=',bnum_S(i)
    write(*,*)'bnum_LTA(i)=',bnum_LTA(i),'bnum_STA(i)=',bnum_STA(i)
   end do
   write(*,*)'bnum_LT=',bnum_LT,'bnum_ST=',bnum_ST

   allocate(n(NBS0,4))
   allocate(cpg_L(NBS_L),cpg_S(NBS_S)) !coefficient of primitive gaussian

   allocate(xc(NAT),yc(NAT),zc(NAT),cn(NAT)) ! position of atom, charge of atom
  !--- primitive gaussian ---
   allocate(aa_L(NBS_L),xx_L(NBS_L),yy_L(NBS_L),zz_L(NBS_L)) ! position of primitive gaussian
   allocate(nx_L(NBS_L),ny_L(NBS_L),nz_L(NBS_L)) ! angular momentum of primitive gaussian
   allocate(aa_S(NBS_S),xx_S(NBS_S),yy_S(NBS_S),zz_S(NBS_S)) ! position of primitive gaussian
   allocate(nx_S(NBS_S),ny_S(NBS_S),nz_S(NBS_S)) ! angular momentum of primitive gaussian
  !--- electron solution coefficients ---
   allocate(c_La(NBS_L,NBS_E), c_Lb(NBS_L,NBS_E), c_Sa(NBS_S,NBS_E), c_Sb(NBS_S,NBS_E)) ! complex coefficients (NBS_L or NBS_S,NBS_E)
  !--- positron solution coefficients ---
   allocate(d_La(NBS_L,NBS_P), d_Lb(NBS_L,NBS_P), d_Sa(NBS_S,NBS_P), d_Sb(NBS_S,NBS_P)) ! complex coefficients (NBS_L or NBS_S,NBS_P) 
   allocate(cont(NBS00,NAT)) ! contract number read from basis.txt
   allocate(e_eig(NBS_E),p_eig(NBS_P))

  open(unit=11,file=FIFI1)
  call readBasis(bnum_L,bnum_S,bnum_LT,bnum_ST,bnum_LTA,bnum_STA,cpg_L,cpg_S,n,cont)
!  set aa_L,xx_L,yy_L,zz_L,aa_S,xx_S,yy_S,zz_S,nx_L,ny_L,nz_L,nx_S,ny_S,nz_S
!  set xc,yc,zc
  close(unit=11)

  open(unit=42,file=FIFI2)
  call readVectors(bnum_L,bnum_S,bnum_LT,bnum_ST,bnum_LTA,bnum_STA,cpg_L,cpg_S,n,cont)
!  set c_La,c_Lb,c_Sa,c_Sb,d_La,d_Lb,d_Sa,d_Sb
  close(unit=42)

!!$ do j=1,NBS
!!$  do i=1,NBS_L
!!$  write(*,*) c_La(i,j),c_Lb(i,j)
!!$  end do
!!$  do k=1,NBS_S
!!$  write(*,*) c_Sa(k,j),c_Sb(k,j)
!!$  end do
!!$ end do
end subroutine read_DiracOutput

!============================================================================
Subroutine read_qedinp
!read input file (qed.inp)
! modified by KI (130313)
!============================================================================
  use DiracOutput

  implicit none

  character(len=80) INPUT

!  INPUT = 'qed.inp'
!  open(unit=11,file=INPUT)
!KI     read (*,'(a)') FILESFOLDER
     FIFI1 = trim(FILESFOLDER)//"/"//"basis.txt"
     FIFI2 = trim(FILESFOLDER)//"/"//"vectors.txt"
!     read (*,'(a12)') FIFI1 ! basis.txt of Dirac10
!     read (*,'(a12)') FIFI2 ! vectors.txt of Dirac10
     write(*,'(a)') FIFI1
     write(*,'(a)') FIFI2
!KI     read (*,'(a)') FIFI3 !Dirac output file of Dirac10
     Fifi3 = trim(FILESFOLDER)//"/"//FIFI3
     write(*,'(a)') FIFI3
!KI     read (*,*) SYMM ! 0=>atom, 1=>molecule using symmmetry calc, 2=>molecule not using symmetry calc
!Ki     read (*,*) KBAL
!     stop
!  close(unit=11)

end subroutine read_qedinp

!============================================================================
subroutine read_spininp
!============================================================================
  use DiracOutput

  implicit none

  character(len=80) INPUT

!  INPUT = 'spin.inp'
!  open(unit=11,file=INPUT)
     read (*,'(a)') FILESFOLDER
     FIFI1 = trim(FILESFOLDER)//"/"//"basis.txt"
     FIFI2 = trim(FILESFOLDER)//"/"//"vectors.txt"
!     read (*,'(a12)') FIFI1 ! basis.txt of Dirac10
!     read (*,'(a12)') FIFI2 ! vectors.txt of Dirac10
     write(*,'(a)') FIFI1
     write(*,'(a)') FIFI2
     read (*,'(a)') FIFI3 !Dirac output file of Dirac10
     FIFI3 = trim(FILESFOLDER)//"/"//FIFI3
     write(*,'(a)') FIFI3
     read (*,*) symm ! 0=>atom, 1=>molecule using symmmetry calc, 2=>molecule not using symmetry calc
     read (*,*) kbal
     read (*,*) NROOT
     read (*,*) u1,u2,u3,uu,u4,u5,u6
!     stop
!  close(unit=11)

end subroutine read_spininp
