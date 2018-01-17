!==================================================================
program SimQED
!==================================================================
  use Precision
  implicit none

  write(*,*) "# Calling time_evolution ... "
  call time_evolution

end program SimQED

!==================================================================
subroutine time_evolution
!================================================================== 
  use Precision
  use DefineTypes
  use DiracOutput
  use Constants
  use IntegralStorage
  use NucBasis
  use Interface_Mod

  implicit none

  character(LEN=80) :: inifile  ! input file which is given in the command line.

  real(kind=dp) :: x,y,z,dx,dy,dz
  integer :: i,j,k,l
  real(kind=dp) :: density
  real(kind=dp) :: density_ele

  real(kind=dp) :: intN_ele
  complex(kind=dp) :: intN_mat

  complex(kind=dp),allocatable :: calE(:,:)  ! expectation value of (e^dagger e)
  complex(kind=dp),allocatable :: calE_past(:,:)  ! expectation value of (e^dagger e) at one time step before.
  complex(kind=dp),allocatable :: dcalE(:,:) 
  complex(kind=dp),allocatable :: h_Qmat(:,:)  ! large matrix (4*NBS x 4*NBS) of one-electron part (used in BO approx. only)
  complex(kind=dp),allocatable :: TM_Qmat(:,:)  ! large matrix (4*NBS x 4*NBS) of kinetic+mass energy
  complex(kind=dp),allocatable :: calFj0_Qmat(:,:,:)  ! large matrix (4*NBS x 4*NBS) of radiation part (at t=0)
  complex(kind=dp),allocatable :: alpha(:)   ! <coherent|hat(a)|coherent>
  complex(kind=dp),allocatable :: twoele_Qmat(:,:,:,:)  ! large matrix (4*NBS)^4 of two-electron integral
!  complex(kind=dp),allocatable :: calJi_Qmat(:,:,:,:,:)  ! large matrix (4*NBS)^4 of calJ^i 
!  complex(kind=dp),allocatable :: calLi_Qmat(:,:,:,:,:)  ! large matrix (4*NBS)^4 of calL^i
  complex(kind=dp),allocatable :: intIJJ_Qmat(:,:,:,:)  ! large matrix (4*NBS)^4 of IJJ_NMPQ

  complex(kind=dp),allocatable :: E_Qmat(:,:,:)  ! large matrix (4*NBS x 4*NBS) of electric field integral (E^k_NM)
  complex(kind=dp),allocatable :: E_Qmatlist(:,:,:,:)  ! list of E_Qmat

!  complex(kind=dp),allocatable :: intJcA_Qmat(:,:,:)  ! large matrix (4*NBS x 4*NBS) of integral IJcA^k_NM
  complex(kind=dp),allocatable :: intJcA_Qmat(:,:)  ! large matrix (4*NBS x 4*NBS) of integral IJcA_NM
  complex(kind=dp),allocatable :: intI2JcA_Qmat(:,:)  ! large matrix (4*NBS x 4*NBS) of integral (-Ze/CCC)*int_M IJcA_NM
  
  complex(kind=dp),allocatable :: calC(:,:)  ! expectation value of (c^dagger c) , assuming one species.
  complex(kind=dp),allocatable :: dcalC(:,:) 
  complex(kind=dp),allocatable :: nucele_Qmat(:,:,:,:)  ! large matrix (NBS_N)^2 x (4*NBS)^2 of nucleus-electron integral
  complex(kind=dp),allocatable :: twonuc_mat(:,:,:,:)  ! (NBS_N)^4 of two-nucleus integral
  complex(kind=dp),allocatable :: Tnuc_mat(:,:)  ! nucleus kinetic energy (NBS_N)^2

  complex(kind=dp),allocatable :: L_Qmat(:,:)  ! to check divergence (4*NBS x 4*NBS)
  complex(kind=dp),allocatable :: L2_Qmat(:,:)  ! to check divergence (4*NBS x 4*NBS)
  complex(kind=dp),allocatable :: calEa(:,:,:)  ! expectation value of (e^dagger a e)
  complex(kind=dp),allocatable :: dcalEa(:,:,:) 

  type(nonzeromat4legs),allocatable :: twoele_Qmat_nz(:)  ! two-electron integral which are larger than threshold
  type(nonzeromat4legs),allocatable :: intIJJ_Qmat_nz(:)  ! IJJ integral which are larger than threshold

  integer :: it ! time counter
  real(kind=dp) :: time 

  !----------------------------------------------------
  !   parameters to be read from inifile
  !----------------------------------------------------
  character(LEN=80) :: outfile 
  integer :: NTIME
  integer :: nprint ! print out every nprint step (it)
  integer :: NEL
  !----------------------------------------------------
  
!  complex(kind=dp) :: dcalEdt,dcalCdt
  integer :: nn,mm,pp,qq,rr
  complex(kind=dp) :: sum, sum_d
  complex(kind=dp) :: sum1,sum2,sum3
 ! complex(kind=dp) :: pos_e,pos_etot,pos_n,pos_ntot
 ! complex(kind=dp) :: density_Qmat, rho_Qmat
 ! complex(kind=dp) :: pos_e0,pos_n0
  complex(kind=dp) :: diag_e,diag_p,diag_e0,diag_p0

  !-------------------------------------
!  integer,parameter :: n_point = 101
  integer,parameter :: n_point = 1
  !-------------------------------------
  real(kind=dp) :: vec_e0(3,n_point),vec_n0(3,n_point)
!  complex(kind=dp) :: poslist_e(n_point),poslist_n(n_point),poslist_e0(n_point),poslist_n0(n_point)

  real(kind=dp) :: tauSeig(3), tauSeiglist(3,n_point), tauSeiglist0(3,n_point) ! eigenvalues of stress tensor (Hermitian part)
  real(kind=dp) :: tauSvec(3,3), tauSveclist(3,3,n_point), tauSveclist0(3,3,n_point) ! eigenvectors of stress tensor (Hermitian part)
  real(kind=dp) :: rhoe, rhoelist(n_point),rhoelist0(n_point) ! electron charge density 
  real(kind=dp) :: rhoe_tot ! total (integration over whole space) electron charge 
  real(kind=dp) :: rhonuc, rhonuclist(n_point),rhonuclist0(n_point) ! nuclear charge density 
  real(kind=dp) :: rhonuc_tot ! total (integration over whole space) nuclear charge 

  real(kind=dp) :: vecR(3)
  real(kind=dp) :: pol(3), pollist(3,n_point), pollist0(3,n_point) ! polarization density 
  real(kind=dp) :: pol_tot(3)
  real(kind=dp) :: ef(3), eflist(3,n_point), eflist0(3,n_point) ! electric field
  real(kind=dp) :: ef_tot(3)

  real(kind=dp) :: vecJ(3)


!  character(LEN=80),parameter :: outfile = 'out_norad_Tm1_3.dat' ! 10  

  complex(kind=dp) :: intSnuc_mat,intTnuc_mat,intF_mat,intcalFj_mat
  complex(kind=dp) :: intnucele_mat,intV_mat,inttwonuc_mat
!  complex(kind=dp) :: rho_nuc_mat!, density_nuc_mat
  real(kind=dp) :: func_pg  ! primitive gaussian function (normalized)
  real(kind=dp) :: vec(3)
  real(kind=dp) :: vec1(3),vec2(3)

  ! for calculation of commutation relation
  complex(kind=dp),allocatable :: calEV(:,:)  ! vacuum expectation value of (e^dagger e)
  complex(kind=dp),allocatable :: dcalEV(:,:) 
  complex(kind=dp),allocatable :: calDV(:,:)  ! vacuum expectation value of (e e^dagger)
  complex(kind=dp),allocatable :: dcalDV(:,:) 
  complex(kind=dp),allocatable :: calCV(:,:)  ! vacuum expectation value of (c^dagger c)
  complex(kind=dp),allocatable :: dcalCV(:,:) 
  complex(kind=dp),allocatable :: calBV(:,:)  ! vacuum expectation value of (c c^dagger)
  complex(kind=dp),allocatable :: dcalBV(:,:) 


  logical :: use_BO_approx    ! .true. --> BO_approximation, .false. --> include nucleus degree of freedom
  logical :: use_exchange    ! .true. --> include exchange terms in diff. eq. , .false. --> without exchange terms
  logical :: use_coherent     ! .true. --> use coherent state for photon, .false. --> no radiation in the initial state.
  integer :: n_coherent ! number of coherent states which are not zero.
  integer :: i_coh ! i_coh-th coherent state is not zero.
  real(kind=dp) :: p0,th,phi ! photon momentum vector in spherical coord.
  character(LEN=1) :: sig ! polarization (+ or -)  
  integer :: ip0,ith,iphi
  real(kind=dp) :: dp0,dth,dphi,dp0_th_phi ! mesh width

  real(kind=dp) :: dAraddt_4E

  logical :: use_Sch  ! .true. --> compute in Schwarzschild spacetime, .false. --> flat spacetime
  real(kind=dp) :: R_Sch ! radius at which computation is performed (normalized by Schwarzshild radius)
  real(kind=dp) :: Mass ! mass of the central objects in units of solar mass

  !-------------------------------
  ! set fundamental constants
  !-------------------------------
  call set_GammaMatrix

  !---------------------------------------------------------------------------------------------
  ! set points to be computed
  !---------------------------------------------------------------------------------------------
  vec_e0(:,:) = 0._dp
  vec_n0(:,:) = 0._dp
  do i=1,n_point
     vec_e0(1,i) = -1._dp +0.02_dp*(i-1)  ! define x-coord
     vec_n0(1,i) = 0._dp +1.e-2_dp*(i-1)
  end do


  !================================================================
  ! set computation parameters
  ! set expansion functions for electrons & nucleus
  ! set photon parameters
  !================================================================
  !---------------------------------------------------------------------------------------------
  ! read computation parameters from the file which is specified in the command line.
  !---------------------------------------------------------------------------------------------
  call getarg(1,inifile)
  write(*,*) "# Read ", trim(inifile)
  open(unit=1000,file=trim(inifile),status='unknown',form='formatted')

  !--------------------------------------------------------------
  ! paramters which are previously given in qed.inp
  read(1000,*)
  read (1000,'(a)') FILESFOLDER
  read(1000,*)
  read (1000,'(a)') FIFI3 !Dirac output file of Dirac10
  read (1000,*) SYMM ! 0=>atom, 1=>molecule using symmmetry calc, 2=>molecule not using symmetry calc
  read (1000,*) KBAL
  !--------------------------------------------------------------

  read(1000,*)
  read(1000,*) outfile
  open(unit=10,file=outfile,status='unknown',form='formatted')  ! output file
  read(1000,*) NTIME
  read(1000,*) DeltaT   ! defined as global parameter
  read(1000,*) nprint
  read(1000,*) NEL  !  NEL = 1 ! H,H2+
                    !  NEL = 2 ! He,H2
                    !  NEL = 6 ! Li2
  read(1000,*) 
  read(1000,*) there_is_twoele
  !  there_is_twoele = .True. ! true : there is already two electron integrals stored in a file whose name is specified below.
  !  there_is_twoele = .false. ! false : compute two electron integrals and store as fort.12
  if(there_is_twoele) then
     read(1000,*) file_twoele ! defined in module IntegralStorage
  end if
  read(1000,*) there_is_intIJJ
  if(there_is_intIJJ) then
     read(1000,*) file_intIJJ ! defined in module IntegralStorage
  end if

  read(1000,*) there_is_nucele
  !  there_is_nucele = .True. ! true : there is already nucleus-electron 4c integrals and have a name defined in module IntegralStorage
  !  there_is_nucele = .false. ! false : compute nucleus-electron 4c integrals and store as fort.13 (--> nucele.dat)
  read(1000,*) there_is_twonuc
  !  there_is_twonuc = .True. ! true : there is already two nucleus integrals and have a name defined in module IntegralStorage
  !  there_is_twonuc = .false. ! false : compute two nucleus integrals and store as fort.14 (--> twonuc.dat)

  read(1000,*) 
  read(1000,*) use_BO_approx
  ! .true. --> BO_approximation, .false. --> include nucleus degree of freedom

  read(1000,*) 
  read(1000,*) use_exchange
  ! .true. --> include exchange terms in diff. eq. , .false. --> without exchange terms

  read(1000,*) 
  read(1000,*) use_Sch
  ! .true. --> compute in Schwarzschild spacetime, .false. --> flat spacetime
  if((use_Sch .eq. .true.).and.(use_BO_approx .eq. .false.)) then
     write(*,*) "Computation in curved spacetime is only supported for BO approximation."
     stop
  end if
  read(1000,*) R_Sch
  read(1000,*) Mass
  if((use_Sch .eq. .true.).and.(R_Sch.le.1._dp)) then
     write(*,*) "R_Sch should be larger than 1."
     stop
  end if

  read(1000,*) 
  read(1000,*) use_coherent
  ! .true. --> use coherent state for photon, .false. --> no radiation in the initial state.
  if((use_coherent .eq. .true.).and.(use_BO_approx .eq. .false.)) then
     write(*,*) "Coherent state is only supported for BO approximation."
     stop
  end if
  if((use_coherent .eq. .true.).and.(use_Sch .eq. .true.)) then
     write(*,*) "Coherent state is only supported for flat spacetime."
     stop
  end if
  if(use_coherent) then
     read(1000,*) P0MAX  ! global
     read(1000,*) NP0    ! global 
     read(1000,*) NPTH   ! global
     read(1000,*) NPPHI  ! global
     Nph = NP0*NPTH*NPPHI*2 ! j=1~Nph, global
     allocate(alpha(Nph))
     alpha(:) = (0._dp,0._dp)
     read(1000,*) n_coherent
     if(n_coherent.gt.Nph) then
        write(*,*) "n_coherent should be smaller than Nph."
        stop
     end if
     do i=1,n_coherent
        read(1000,*) i_coh
        read(1000,*) alpha(i_coh)
     end do
  end if


 
  write(*,*) "#############################################"
  if(use_BO_approx) then
     write(*,*) "# Use BO approximation."
     write(10,*) "# Use BO approximation."
  else  
     write(*,*) "# No BO approximation."
     write(10,*) "# No BO approximation."
  end if
  if(use_exchange) then
     write(*,*) "# Use exchange terms."
     write(10,*) "# Use exchange terms."
  else  
     write(*,*) "# Without exchange terms."
     write(10,*) "# Without exchange terms."
  end if
  if(use_Sch) then
     write(*,*) "# Compute in Schwarzschild spacetime."
     write(10,*) "# Compute in Schwarzschild spacetime."
  else  
     write(*,*) "# Compute in flat spacetime."
     write(10,*) "# Compute in flat spacetime."
  end if
  write(*,*) "#############################################"

  write(*,*) "# Write outputs to ", trim(outfile)
  write(10,*) "# Write outputs to ", trim(outfile)
  write(*,"(1a15,1i16)") "# NTIME : ", NTIME
  write(10,"(1a15,1i16)") "# NTIME : ", NTIME
  write(*,"(1a15,1es16.6)") "# DeltaT : ", DeltaT
  write(10,"(1a15,1es16.6)") "# DeltaT : ", DeltaT
  write(*,"(1a15,1es16.6)") "# t_end : ", NTIME*DeltaT
  write(10,"(1a15,1es16.6)") "# t_end : ", NTIME*DeltaT
  write(*,"(1a15,1i16)") "# nprint : ", nprint
  write(10,"(1a15,1i16)") "# nprint : ", nprint

  !---------------------------------------------------------------------------------------------
  ! set electron related parameters
  !---------------------------------------------------------------------------------------------
  write(*,*) "#############################################"
  write(*,*) "Set expansion functions for electrons."
  write(10,*) "# Set expansion functions for electrons."
  write(*,*) "#############################################"

  call read_qedinp ! read qed.inp
  call read_DiracOutput  ! read and set global variables -> defined in sub_readDiracOutput.f90

  if(mod(NEL,2).eq.0) then
     NOCC = NEL/2  ! closed shell (KR)
  else
     NOCC = (NEL+1)/2 ! open shell (KR)
     ! The occupation number of the NOCCth orbital is 1.
  end if

  write(*,"(1a15,1i6)") "# NEL : ",NEL
  write(10,"(1a15,1i6)") "# NEL : ",NEL
  write(*,"(1a15,1i6)") "# NBS : ",NBS
  write(10,"(1a15,1i6)") "# NBS : ",NBS
  write(*,"(1a15,1i6)") "# NOCC : ",NOCC
  write(10,"(1a15,1i6)") "# NOCC : ",NOCC

  write(*,"(1a15,1i6)") "# NAT : ",NAT
  write(10,"(1a15,1i6)") "# NAT : ",NAT
  do i=1,NAT
     write(*,"(1a15,1i6,3es16.6)") "# cn & xyzc: ",cn(i),xc(i),yc(i),zc(i)
     write(10,"(1a15,1i6,3es16.6)") "# : ",cn(i),xc(i),yc(i),zc(i)
  end do
!  stop


  ! check total charge and orthogonality of molecular orbitals.
  write(*,"(1a15,1es16.8)") "# intN_ele: ", intN_ele(NEL)
  do i=1,NBS
     do j=1,NBS
        write(*,'(2i6,8es16.6)')  i, j,intN_mat(i,"+",j,"+")
        write(*,'(2i6,8es16.6)')  i, j,intN_mat(i,"+",j,"-")
        write(*,'(2i6,8es16.6)')  i, j,intN_mat(i,"-",j,"+")
        write(*,'(2i6,8es16.6)')  i, j,intN_mat(i,"-",j,"-")
     end do
  end do


  if(.not.use_BO_approx) then
!!$     !---------------------------------------------------------------------------------------------
!!$     ! set nucleus related parameters
!!$     !---------------------------------------------------------------------------------------------
!!$     write(*,*) "#############################################"
!!$     write(*,*) "Set expansion functions for atomic nuclei."
!!$     write(10,*) "# Set expansion functions for atomic nuclei."
!!$     write(*,*) "#############################################"
!!$     !call set_NucBasis ! for H
!!$       call set_NucBasis_H2 ! for H2
!!$     if(NUCTYPE.eq.0) then  ! boson
!!$        NBS_N = NBS_PHI
!!$     elseif(NUCTYPE.eq.1) then  ! fermion
!!$        NBS_N = 2*NBS_PHI
!!$     else
!!$        write(*,*) "NUCTYPE should be 0 or 1."
!!$        stop
!!$     end if
!!$     
!!$     write(*,"(1a15,1i6)") "# NBS_PHI : ",NBS_PHI
!!$     write(10,"(1a15,1i6)") "# NBS_PHI : ",NBS_PHI
!!$     write(*,"(1a15,1i6)") "# N_PGN : ",N_PGN
!!$     write(10,"(1a15,1i6)") "# N_PGN : ",N_PGN
!!$     write(*,"(1a15,1es16.6)") "# L_well : ",L_well
!!$     write(10,"(1a15,1es16.6)") "#L_well : ",L_well
!!$     write(*,"(1a15,1es16.6)") "# L_pg : ",L_pg
!!$     write(10,"(1a15,1es16.6)") "#L_pg : ",L_pg
!!$     
!!$     write(*,"(1a15,1es16.6)") "# alpha_N : ",alpha_N
!!$     write(10,"(1a15,1es16.6)") "# alpha_N : ",alpha_N
!!$     write(*,"(1a15,1es16.6)") "# m_N : ",m_N
!!$     write(10,"(1a15,1es16.6)") "# m_N : ",m_N
!!$     write(*,"(1a15,1i16)") "# NUCTYPE : ",NUCTYPE
!!$     write(10,"(1a15,1i16)") "# NUCTYPE : ",NUCTYPE
!!$     write(*,"(1a15,1i16)") "# NBS_N : ",NBS_N
!!$     write(10,"(1a15,1i16)") "# NBS_N : ",NBS_N
!!$     
!!$     
!!$     !primordial gaussian function positions
!!$     do i=1,N_PGN
!!$        write(*,"(1a15,1i6,3es16.6)") "# PGN pos",i,vecR_N(1,i),vecR_N(2,i),vecR_N(3,i)
!!$        write(10,"(1a15,1i6,3es16.6)") "# PGN pos",i,vecR_N(1,i),vecR_N(2,i),vecR_N(3,i)
!!$     end do
!!$     
!!$     !check overlap
!!$     do i=1,NBS_N
!!$        do j=1,NBS_N
!!$           write(98,"(2i6,3es16.6)") i,j,intSnuc_mat(i,j)
!!$        end do
!!$     end do
!!$     
  end if

  if(use_Sch) then
     write(*,*) "#############################################"
     write(*,*) "Parameters for computation in curved spacetime."
     write(10,*) "# Parameters for computation in curved spacetime."
     write(*,*) "#############################################"
     write(*,"(1a15,1es16.6)") "# R_Sch : ",R_Sch
     write(10,"(1a15,1es16.6)") "# R_Sch : ",R_Sch
     write(*,"(1a15,1es16.6)") "# Mass : ",Mass
     write(10,"(1a15,1es16.6)") "# Mass : ",Mass
  end if

  if(use_coherent) then
     !---------------------------------------------------------------------------------------------
     ! set photon parameters & specify coherent states
     !---------------------------------------------------------------------------------------------
     write(*,*) "#############################################"
     write(*,*) "Set photon parameters."
     write(10,*) "# Set photon parameters."
     write(*,*) "#############################################"
     allocate(ALPHA_COH(Nph))
     ALPHA_COH(:) = alpha(:) 

     write(*,"(1a15,1es16.6)") "# P0MAX: ",P0MAX
     write(10,"(1a15,1es16.6)") "# P0MAX: ",P0MAX
     write(*,"(1a15,1i6)") "# NP0: ",NP0
     write(10,"(1a15,1i6)") "# NP0: ",NP0
     write(*,"(1a15,1i6)") "# NPTH: ",NPTH
     write(10,"(1a15,1i6)") "# NPTH: ",NPTH
     write(*,"(1a15,1i6)") "# NPPHI: ",NPPHI
     write(10,"(1a15,1i6)") "# NPPHI: ",NPPHI
     write(*,"(1a15,1i6)") "# Nph: ",Nph
     write(10,"(1a15,1i6)") "# Nph: ",Nph
     
     if(n_coherent.ge.1) then
        write(*,*) "# Photon momentum discretization"
        write(10,*) "# Photon momentum discretization"
        write(*,*) "# Polar coordinates: th measured from z-axis, include both 0 and PI."
        write(*,*) "#                  : phi measured from x-axis, include include 0 but not 2PI."
        write(10,*) "# Polar coordinates: th measured from z-axis, include both 0 and PI."
        write(10,*) "#                  : phi measured from x-axis, include include 0 but not 2PI."
        
        dp0  = P0MAX/NP0 ! P0MAX specified in GammaMatrix. We do not use p0=0
        dth  = PI/(NPTH-1) ! NPTH specified in GammaMatrix. th include both 0 and PI
        dphi = 2._dp*PI/NPPHI ! NPPHI specified in GammaMatrix. phi include 0 but not 2PI
        
        write(*,"(1a3,1a6,1a5,1a16,2a16)") "#  "," j ","pol"," p0      "," th/PI ", " phi/PI "
        write(10,"(1a3,1a6,1a5,1a16,2a16)") "#  "," j ","pol"," p0      "," th/PI ", " phi/PI "
        do j=1,Nph
           call convert_label_2(j,NP0,NPTH,NPPHI,ip0,ith,iphi,sig)  ! defined in sub_int.f90
           p0 = dp0*ip0
           th = dth*(ith-1)  ! start from 0
           phi= dphi*(iphi-1) ! start from 0
           write(*,"(1a3,1i6,1a5,1es16.6,2f16.4)") "# ", j, sig, p0, th/PI, phi/PI
           write(10,"(1a3,1i6,1a5,1es16.6,2f16.4)") "# ", j, sig, p0, th/PI, phi/PI
        end do
        
        write(*,*) "# Print non-zero alpha"
        write(10,*) "# Print non-zero alpha"
        j=0
        do i=1,Nph
           if(ALPHA_COH(i).ne.(0._dp,0._dp)) then
              write(*,"(1a8,1i4,2es16.6)") "# alpha",i,ALPHA_COH(i)
              write(10,"(1a8,1i4,2es16.6)") "# alpha",i,ALPHA_COH(i)
              j=j+1
           end if
        end do
        if(j.ne.n_coherent) then
           write(*,*) "n_coherent does not match. Check input file."
           stop
        end if
     end if
     
  end if
!  stop


  


  !================================================================
  ! allocate density matrices.
  ! set or compute initial condition for them.
  ! compute physical quantities at initial time.
  !================================================================
  allocate(calE(4*NBS,4*NBS),dcalE(4*NBS,4*NBS))
  allocate(calE_past(4*NBS,4*NBS))
  ! for VEV of commutation relation
  allocate(calEV(4*NBS,4*NBS),dcalEV(4*NBS,4*NBS))
  allocate(calDV(4*NBS,4*NBS),dcalDV(4*NBS,4*NBS))

  if(.not.use_BO_approx) then
     allocate(calC(NBS_N,NBS_N),dcalC(NBS_N,NBS_N))
     ! for VEV of commutation relation
     allocate(calCV(NBS_N,NBS_N),dcalCV(NBS_N,NBS_N))
     allocate(calBV(NBS_N,NBS_N),dcalBV(NBS_N,NBS_N))
  end if

  if(use_coherent .and. (n_coherent.eq.0)) then
     allocate(calEa(4*NBS,4*NBS,Nph),dcalEa(4*NBS,4*NBS,Nph))
  end if



  !----------------------------------------------------------
  ! set initial condition for density matrices
  !----------------------------------------------------------

  ! Initial condition for calE
  ! Kronecker delta for electron occupied.
  calE(:,:) = (0._dp,0._dp)
  do i=1,NEL
     calE(i,i) = (1._dp,0._dp)
  end do
!  calE(3,3) = (1._dp,0._dp)

  calE_past(:,:) = (0._dp,0._dp)  ! calE=0 for t<0

  ! Initial condition for calEV
  ! Kronecker delta for "-" (positron)
  calEV(:,:) = (0._dp,0._dp)
  do i=2*NBS+1,4*NBS
     calEV(i,i) = (1._dp,0._dp)
  end do
  ! Initial condition for calDV
  ! Kronecker delta for "+" (electron)
  calDV(:,:) = (0._dp,0._dp)
  do i=1,2*NBS
     calDV(i,i) = (1._dp,0._dp)
  end do


  if(.not.use_BO_approx) then
     
     ! Initial condition for calC
     calC(:,:) = (0._dp,0._dp)
     
     ! H atom, ground state
     calC(1,1) = (1._dp,0._dp)
     
     
!!$  ! ortho H2 (triplet, +1)
!!$  calC(1,1) = (1._dp,0._dp)
!!$  calC(3,3) = (1._dp,0._dp)
!!$
!!$  ! ortho H2 (triplet, 0)
!!$  calC(1,1) = (0.5_dp,0._dp)
!!$  calC(2,2) = (0.5_dp,0._dp)
!!$  calC(3,3) = (0.5_dp,0._dp)
!!$  calC(4,4) = (0.5_dp,0._dp)
!!$
!!$     ! ortho H2 (triplet, -1)
!!$     calC(2,2) = (1._dp,0._dp)
!!$     calC(4,4) = (1._dp,0._dp)
     
!!$  ! para H2 (singlet)
!!$  calC(1,1) = (1._dp,0._dp)
!!$  calC(2,2) = (1._dp,0._dp)
!!$
     
     ! Initial condition for calCV
     ! zero matrix
     calCV(:,:) = (0._dp,0._dp)
     ! Initial condition for calBV
     ! Kronecker delta
     calBV(:,:) = (0._dp,0._dp)
     do i=1,NBS_N
        calBV(i,i) = (1._dp,0._dp)
     end do
     
  end if
  
  if(use_coherent .and. (n_coherent.eq.0)) then
     calEa(:,:,:) = (0._dp,0._dp)
  end if
  

  write(*,*) "#############################################"
  write(*,*) "Compute initial physical quantities."
  write(10,*) "# Compute initial physical quantities."
  write(*,*) "#############################################"

  ! for calculating polarization
  allocate(E_Qmat(3,4*NBS,4*NBS)) 
  allocate(E_Qmatlist(3,4*NBS,4*NBS,n_point)) 
  do i=1,n_point
     vecR(:) = vec_e0(:,i)
     call setQmat_E(vecR,E_Qmat)
     E_Qmatlist(:,:,:,i) = E_Qmat(:,:,:)
  end do

  !----------------------------------------------------------
  ! compute initial value of physical quantities
  !----------------------------------------------------------

  do i=1,n_point
     x = vec_e0(1,i)
     y = vec_e0(2,i)
     z = vec_e0(3,i)

     ! initial electron charge density 
     call calc_rhoe(x,y,z,calE,rhoe)
     rhoelist0(i) = rhoe
     ! output initial electron charge density 
     write(*,"(1a15,1i6,4es16.6)") "# rhoe0: ",i,x,y,z,rhoelist0(i)
     write(10,"(1a15,1i6,4es16.6)") "# rhoe0: ",i,x,y,z,rhoelist0(i)

     ! initial electronic stress tensor density 
     call calc_tauS_eigvec(x,y,z,calE,tauSeig,tauSvec)
     do k=1,3
        tauSeiglist0(k,i) = tauSeig(k)
        do l=1,3
           tauSveclist0(l,k,i) = tauSvec(l,k)
        end do
     end do

     vecR(:) = vec_e0(:,i) 

     if(use_BO_approx) then
        ! initial polarization density
        E_Qmat(:,:,:) = E_Qmatlist(:,:,:,i) 
        call calc_pol_BO(vecR,calE,E_Qmat,pol)
        pollist0(:,i) = pol(:)
     end if

     ! initial (-1/c)dArad^k/dt 
     do k=1,3
        eflist0(k,i) = dAraddt_4E(k,0._dp,vecR,alpha)
     end do
  end do

  ! output initial polarization density
  do i=1,n_point
     x = vec_e0(1,i);y = vec_e0(2,i);z = vec_e0(3,i)
     write(*,"(1a15,1i6,6es16.6)") "# pol0: ",i,x,y,z,pollist0(1,i),pollist0(2,i),pollist0(3,i)
     write(10,"(1a15,1i6,6es16.6)") "# pol0: ",i,x,y,z,pollist0(1,i),pollist0(2,i),pollist0(3,i)
  end do

  ! output initial electric field
  do i=1,n_point
     x = vec_e0(1,i);y = vec_e0(2,i);z = vec_e0(3,i)
     write(*,"(1a15,1i6,6es16.6)") "# ef0: ",i,x,y,z,eflist0(1,i),eflist0(2,i),eflist0(3,i)
     write(10,"(1a15,1i6,6es16.6)") "# ef0: ",i,x,y,z,eflist0(1,i),eflist0(2,i),eflist0(3,i)
  end do
!  stop

  ! output initial electronic stress tensor (eigenvalues and eigenvectors)
  do k=1,3
     do i=1,n_point
        x = vec_e0(1,i);y = vec_e0(2,i);z = vec_e0(3,i)
        write(*,"(1a15,1i1,1a3,1i6,7es16.6)") "# tauS0 eig(",k,"): ",i,x,y,z,tauSeiglist0(k,i),tauSveclist0(1,k,i),tauSveclist0(2,k,i),tauSveclist0(3,k,i)
        write(10,"(1a15,1i1,1a3,1i6,7es16.6)") "# tauS0 eig(",k,"): ",i,x,y,z,tauSeiglist0(k,i),tauSveclist0(1,k,i),tauSveclist0(2,k,i),tauSveclist0(3,k,i)
     end do
  end do
  
!  stop

  if(.not.use_BO_approx) then
!!$     do i=1,n_point
!!$        x = vec_n0(1,i)
!!$        y = vec_n0(2,i)
!!$        z = vec_n0(3,i)
!!$        ! initial nucleus charge density etc
!!$        call calc_rhonuc(x,y,z,calC,rhonuc)
!!$        rhonuclist0(i)=rhonuc
!!$        write(*,"(1a15,1i6,5es16.6)") "# rho_n0: ",i,x,y,z,rhonuclist0(i)
!!$        write(10,"(1a15,1i6,5es16.6)") "# rho_n0: ",i,x,y,x,rhonuclist0(i)
!!$     end do
  end if

  !================================================================
  ! set up (compute or read) integrals
  !================================================================
  if(use_BO_approx) then
     allocate(h_Qmat(4*NBS,4*NBS))  
     allocate(calFj0_Qmat(4*NBS,4*NBS,Nph))  ! Arad is only supported for BO at present
     if(use_coherent .and. (n_coherent.eq.0)) then
!        allocate(L_Qmat(4*NBS,4*NBS)) 
!        allocate(L2_Qmat(4*NBS,4*NBS)) 
     end if
  else
!     allocate(TM_Qmat(4*NBS,4*NBS)) 
!     allocate(nucele_Qmat(NBS_N,NBS_N,4*NBS,4*NBS))
!     allocate(twonuc_mat(NBS_N,NBS_N,NBS_N,NBS_N))
!     allocate(Tnuc_mat(NBS_N,NBS_N))
  end if
  
  allocate(twoele_Qmat(4*NBS,4*NBS,4*NBS,4*NBS))

  write(*,*) "#############################################"
  write(*,*) "Setting matrix of integrals ... "
  write(10,*) "# Setting matrix of integrals ... "
  write(*,*) "#############################################"

  ! need T+M for non-BO (V is not needed)

  if(use_BO_approx) then
     if(use_Sch) then
!        call setQmat_hsch(R_Sch,Mass,h_Qmat)
     else
        call setQmat_h(h_Qmat)
     end if
     call setQmat_calFj0(calFj0_Qmat)
  else
!     call setQmat_TM(TM_Qmat)
!     call setQmat_nucele(nucele_Qmat)
!     call setmat_twonuc(twonuc_mat)
!     call setmat_Tnuc(Tnuc_mat)
  end if
 
!!$  call setQmat_twoele(twoele_Qmat)
!!$
  call setQmat_twoele_nz(twoele_Qmat_nz)
  !unpack
  twoele_Qmat = (0._dp,0._dp)
  do i=1,N_twoele
     nn = twoele_Qmat_nz(i)%a
     mm = twoele_Qmat_nz(i)%b
     pp = twoele_Qmat_nz(i)%c
     qq = twoele_Qmat_nz(i)%d
     twoele_Qmat(nn,mm,pp,qq) = twoele_Qmat_nz(i)%val
  end do

!  call setQmat_twoele_bin(twoele_Qmat)


  !================================================================
  ! compute time evolution (solve diff. eq.)
  !================================================================
  write(*,*) "#############################################"
  write(*,*) "Start time evolution loop ... "
  write(10,*) "# Start time evolution loop ... "
  write(*,*) "#############################################"

  ! loop for time evolution
  do it = 0,NTIME  
     ! it starts from 0 -> t starts from 0
     time = DeltaT*it  ! DeltaT is defined in GammaMatrix and set in the main routine.

     !print out every nprint steps.
     if(mod(it,nprint).eq.0) then

        ! compute expectation values of physical quantities
        do i=1,n_point
           x = vec_e0(1,i)
           y = vec_e0(2,i)
           z = vec_e0(3,i)

           ! charge density
           call calc_rhoe(x,y,z,calE,rhoe)
           rhoelist(i) = rhoe

           
!!$           ! stress tensor density 
!!$           call calc_tauS_eigvec(x,y,z,calE,tauSeig,tauSvec)
!!$           do k=1,3
!!$              tauSeiglist(k,i) = tauSeig(k)
!!$              do l=1,3
!!$                 tauSveclist(l,k,i) = tauSvec(l,k)
!!$              end do
!!$           end do

           vecR(:) = vec_e0(:,i)

           if(use_BO_approx) then
              ! polarization density
              E_Qmat(:,:,:) = E_Qmatlist(:,:,:,i) 
              call calc_pol_BO(vecR,calE,E_Qmat,pol)
              pollist(:,i) = pol(:)
           end if
           
           ! (-1/c)dArad^k/dt 
           do k=1,3
              eflist(k,i) = dAraddt_4E(k,time,vecR,alpha)
           end do

        end do
        ! compute total electron charge
        call calc_rhoe_tot(calE,rhoe_tot)
        
        if(.not.use_BO_approx) then
!!$           ! compute nucleus charge density
!!$           do i=1,n_point
!!$              x = vec_n0(1,i)
!!$              y = vec_n0(2,i)
!!$              z = vec_n0(3,i)
!!$              call calc_rhonuc(x,y,z,calC,rhonuc)
!!$              rhonuclist(i) = rhonuc
!!$           end do
!!$           ! compute total nucleus charge (trace of calC)
!!$           call calc_rhonuc_tot(calC,rhonuc_tot)
        end if


        !output time evolution of physical quantities
        write(10,"(1i6,1es16.6)",advance="no") it/nprint, time
        write(10,"(1es16.6)",advance="no") rhoe_tot -(-NEL)
        do i=1,n_point
           write(10,"(1es16.6)",advance="no") rhoelist(i)-rhoelist0(i)
           write(10,"(3es16.6)",advance="no") pollist(1,i)-pollist0(1,i),pollist(2,i)-pollist0(2,i),pollist(3,i)-pollist0(3,i)
           write(10,"(3es16.6)",advance="no") eflist(1,i),eflist(2,i),eflist(3,i)
!!$           do k=1,3
!!$              write(10,"(1es16.6)",advance="no")  tauSeiglist(k,i)-tauSeiglist0(k,i)
!!$           end do
        end do
        
        if(.not.use_BO_approx) then
           write(10,"(1es16.6)",advance="no") rhonuc_tot -NNUC
           do i=1,n_point
              write(10,"(1es16.6)",advance="no")  rhonuclist(i)-rhonuclist0(i)
!              write(10,"(2es16.6)",advance="no")  poslist_n(i)-poslist_n0(i)
           end do
        end if

        write(10,*)
     end if

     if(use_BO_approx) then
        if(use_coherent .and. (n_coherent.eq.0)) then
!           call calc_derivatives_BO_N(time,calE,calEa,h_Qmat,calFj0_Qmat,dcalE,dcalEa)
        else
!           call calc_derivatives_BOwithEx(use_exchange,time,calE,calEV,calDV, &
!                & h_Qmat,twoele_Qmat,calFj0_Qmat,alpha,dcalE,dcalEV,dcalDV)
!           call calc_derivatives_BO_Phi_test(time,calE,h_Qmat,calFj0_Qmat,twoele_Qmat_nz,dcalE)
           call calc_derivatives_BO_Phi_test(time,calE,h_Qmat,calFj0_Qmat,twoele_Qmat,dcalE)
!           call calc_derivatives_BO_Phi(time,calE,calE_past,h_Qmat,calFj0_Qmat,twoele_Qmat,intIJJ_Qmat,dcalE)
        end if
     else
!        call calc_derivatives_NonBOwithEx(use_exchange,time,calE,calC, & 
!             & TM_Qmat,Tnuc_mat,twoele_Qmat,nucele_Qmat,twonuc_mat,calFj0_Qmat,alpha, &
!             & dcalE,dcalC)
        ! VEV computation is not supported for with exchange case
     end if

     ! update electron
     do pp=1,4*NBS
        do qq=1,4*NBS
           calE_past(pp,qq)  = calE(pp,qq) 
           calE(pp,qq)  = calE(pp,qq)  + dcalE(pp,qq)
!           calEV(pp,qq) = calEV(pp,qq) + dcalEV(pp,qq)
!           calDV(pp,qq) = calDV(pp,qq) + dcalDV(pp,qq)
!           write(10,"(2i6,4es16.6)") pp,qq,dcalE(pp,qq),L2_Qmat(pp,qq)*(-1._dp/(2._dp*PI**2*CCC))*(DeltaT**2)
        end do
     end do

     if(.not.use_BO_approx) then
        ! update nucleus
        do i=1,NBS_N
           do j=1,NBS_N
              calC(i,j)  = calC(i,j) + dcalC(i,j)
              !           calCV(i,j) = calCV(i,j) + dcalCV(i,j)
              !           calBV(i,j) = calBV(i,j) + dcalBV(i,j)
           end do
        end do
     end if

     if(use_coherent .and. (n_coherent.eq.0)) then
        do pp=1,4*NBS
           do qq=1,4*NBS
              do j=1,Nph
                 calEa(pp,qq,j)  = calEa(pp,qq,j)  + dcalEa(pp,qq,j)
              end do
           end do
        end do
     end if
     
  end do
  
  write(*,*) "Time loop done."
  
  
  deallocate(cn,xc,yc,zc)  ! global
  deallocate(aa_L,xx_L,yy_L,zz_L,nx_L,ny_L,nz_L)
  deallocate(aa_S,xx_S,yy_S,zz_S,nx_S,ny_S,nz_S)
  deallocate(c_La,c_Lb,c_Sa,c_Sb)
  deallocate(d_La,d_Lb,d_Sa,d_Sb)

  deallocate(twoele_Qmat)
  deallocate(twoele_Qmat_nz)
  deallocate(E_Qmat) 
  deallocate(calE,dcalE)
  deallocate(calEV,dcalEV,calDV,dcalDV)
  
  if(use_BO_approx) then
     deallocate(h_Qmat)
     deallocate(calFj0_Qmat)
  else
     deallocate(TM_Qmat,Tnuc_mat)
     deallocate(nucele_Qmat,twonuc_mat)
     deallocate(calC,dcalC)
     deallocate(calCV,dcalCV,calBV,dcalBV)
  end if

  if(use_coherent) then
     deallocate(alpha)
  end if

  if(use_coherent .and. (n_coherent.eq.0)) then
     deallocate(calEa,dcalEa)
  end if
  
end subroutine time_evolution

!==================================================================
subroutine set_GammaMatrix
! "Standard representation"
!  set PI 
! 121115 added gamma0*gammak = alpha
!==================================================================
  use Precision
  use Constants
  implicit none

  integer :: a,b,c,d,i,j,k,l
  complex(kind=dp) :: sum

! In this comments, i-> imaginary unit, I -> 2x2 unit matrix 
!
! Pauli matrices
! sig^1 = (0 1)    sig^2 = (0 -i)   sig^3 = (1  0)
!         (1 0)            (i  0)           (0 -1)

!  complex(kind=dp) :: Gam0(4,4),Gam(3,4,4) ! 3 denotes gamma1,gamma2,gamma3
!  gamma^0 = beta = (I  0)      gamma^k = (   0    sig^k)
!                   (0 -I)                (-sig^k    0  )    
!  complex(kind=dp) :: Gam5(4,4)  ! gamma_5 = -gamma^5
!  gamma_5 = ( 0 | I )
!            ( I | 0 )
!  complex(kind=dp) :: Sigma(3,4,4) ! 4x4 Sigma matrix
!  Sigma^k = (sig^k    0  )
!            (  0    sig^k)
!  complex(kind=dp) :: Gam0k(3,4,4) ! gamma^0 gamma^k = alpha
!  gamma^0 gamma^k = alpha = (  0    sig^k)
!                            (sig^k    0  )

  PI = 4*atan(1._dp)
  
  Gam0 = (0._dp,0._dp); Gam = (0._dp,0._dp); Gam5 = (0._dp,0._dp)
  Sigma = (0._dp,0._dp); Gam0k = (0._dp,0._dp)

  ! gamma^0 = beta
  Gam0(1,1)=( 1._dp,0._dp)
  Gam0(2,2)=( 1._dp,0._dp)
  Gam0(3,3)=(-1._dp,0._dp)
  Gam0(4,4)=(-1._dp,0._dp)

  ! gamma^k
  Gam(1,1,4)=  (1._dp,0._dp)
  Gam(1,2,3)=  (1._dp,0._dp)
  Gam(1,3,2)= (-1._dp,0._dp)
  Gam(1,4,1)= (-1._dp,0._dp)

  Gam(2,1,4)=  (0._dp,-1._dp)
  Gam(2,2,3)=  (0._dp, 1._dp)
  Gam(2,3,2)=  (0._dp, 1._dp)
  Gam(2,4,1)=  (0._dp,-1._dp)
  
  Gam(3,1,3)=( 1._dp,0._dp)
  Gam(3,2,4)=(-1._dp,0._dp)
  Gam(3,3,1)=(-1._dp,0._dp)
  Gam(3,4,2)=( 1._dp,0._dp)

  ! gamma_5 ( = -gamma^5)
  Gam5(1,3)=( 1._dp,0._dp)
  Gam5(2,4)=( 1._dp,0._dp)
  Gam5(3,1)=( 1._dp,0._dp)
  Gam5(4,2)=( 1._dp,0._dp)

  ! Sigma^k
  Sigma(1,1,2)=( 1._dp,0._dp)
  Sigma(1,2,1)=( 1._dp,0._dp)
  Sigma(1,3,4)=( 1._dp,0._dp)
  Sigma(1,4,3)=( 1._dp,0._dp)

  Sigma(2,1,2)=( 0._dp,-1._dp)
  Sigma(2,2,1)=( 0._dp, 1._dp)
  Sigma(2,3,4)=( 0._dp,-1._dp)
  Sigma(2,4,3)=( 0._dp, 1._dp)

  Sigma(3,1,1)=( 1._dp,0._dp)
  Sigma(3,2,2)=(-1._dp,0._dp)
  Sigma(3,3,3)=( 1._dp,0._dp)
  Sigma(3,4,4)=(-1._dp,0._dp)

  ! gamma0 gamma^k = alpha
  Gam0k(1,1,4)=  (1._dp,0._dp)
  Gam0k(1,2,3)=  (1._dp,0._dp)
  Gam0k(1,3,2)=  (1._dp,0._dp)
  Gam0k(1,4,1)=  (1._dp,0._dp)

  Gam0k(2,1,4)=  (0._dp,-1._dp)
  Gam0k(2,2,3)=  (0._dp, 1._dp)
  Gam0k(2,3,2)=  (0._dp,-1._dp)
  Gam0k(2,4,1)=  (0._dp, 1._dp)
  
  Gam0k(3,1,3)=( 1._dp,0._dp)
  Gam0k(3,2,4)=(-1._dp,0._dp)
  Gam0k(3,3,1)=( 1._dp,0._dp)
  Gam0k(3,4,2)=(-1._dp,0._dp)

  ! sum_{k=1}^3 [gamma^0 gamma^k]_{alpha beta} [gamma^0 gamma^k]_{gamma delta}
  i=1
  do a=1,4
     do b=1,4
        do c=1,4
           do d=1,4
              sum = (0._dp,0._dp)
              do k=1,3
                 sum = sum + Gam0k(k,a,b)*Gam0k(k,c,d)
              end do
              !              GamJJ(a,b,c,d) = sum
              if(sum.ne.(0._dp,0._dp)) then
                 GamJJ(i)%a = a
                 GamJJ(i)%b = b
                 GamJJ(i)%c = c
                 GamJJ(i)%d = d
                 GamJJ(i)%val = sum
                 i = i+1
              !                 write(*,"(4i6,2es16.6)") a,b,c,d,GamJJ(a,b,c,d)
              end if
           end do
        end do
     end do
  end do
  
  return
end subroutine set_GammaMatrix



