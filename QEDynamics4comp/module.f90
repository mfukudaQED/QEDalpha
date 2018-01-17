!=============================================================================
! 2012.7.10
! - gathered modules (and subroutine SetGammMatrix) previously defined in simQED.f90
! - added module Precision --> used in almost every routines.
!
! 11.16
! - added module System_Medium for parameters of System (A) and Medium (M)
!
!
!=============================================================================



!===================================================================
module Precision
!===================================================================
  implicit none

!  integer, parameter :: qp = kind(1.q0)
  integer, parameter :: dp = kind(1.d0)
!  integer, parameter :: sp = kind(1.0)

end module Precision

!===================================================================
module DiracOutput
!===================================================================
  use Precision
  implicit none
  
  character(LEN=200) :: FILESFOLDER
  Character(len=120) FIFI1,FIFI2,FIFI3
  integer SYMM
  integer KBAL !KBAL->1  ,URKBAL->0
  integer NROOT
  real(kind=dp):: U1,U2,U3,U4,U5,U6,UU

  integer :: NBS  ! total number of molecular orbitals
  integer :: NBS_E, NBS_P !number of molecular Electron and Positron orbitals
  integer :: NBS_TOTAL !NBS_E + NBS_P
  integer :: NBS00  ! total number of primitive gaussians (NBS_L+NBS_S)
  integer :: NBS0  !  number of primitive gaussians (max(NBS_L,NBS_S))

  integer :: NAT ! number of atoms
  real(kind=dp), allocatable :: xc(:),yc(:),zc(:) ! position of atom
  integer, allocatable :: cn(:) ! charge of atom

  !--- primitive gaussian ---
  integer :: NBS_L, NBS_S  ! number of primitive gaussians for large and small components
  real(kind=dp), allocatable :: aa_L(:),xx_L(:),yy_L(:),zz_L(:) ! position of primitive gaussian
  integer, allocatable :: nx_L(:),ny_L(:),nz_L(:) ! angular momentum of primitive gaussian
  real(kind=dp), allocatable :: aa_S(:),xx_S(:),yy_S(:),zz_S(:) ! position of primitive gaussian
  integer, allocatable :: nx_S(:),ny_S(:),nz_S(:) ! angular momentum of primitive gaussian

  !--- electron solution coefficients ---
  complex(kind=dp), allocatable :: c_La(:,:), c_Lb(:,:), c_Sa(:,:), c_Sb(:,:) ! complex coefficients (NBS_L or NBS_S,NBS)

  !--- positron solution coefficients ---
  complex(kind=dp), allocatable :: d_La(:,:), d_Lb(:,:), d_Sa(:,:), d_Sb(:,:) ! complex coefficients (NBS_L or NBS_S,NBS)

  !---orbital energy---
  double precision, allocatable :: e_eig(:),p_eig(:)

end module DiracOutput

!===================================================================
module DefineTypes
!===================================================================
  use Precision
  implicit none

  integer,parameter :: NMAX_PG = 2000  ! Assume NMAX_PG primitive gaussian functions at maximum for both large and small components.
  
  type primitive_gaussian
     ! Large (1,2) components
     real(kind=dp) :: aL(NMAX_PG)                  ! exponents 
     real(kind=dp) :: xL(NMAX_PG),yL(NMAX_PG),zL(NMAX_PG)      ! positions of p.g.
     integer      :: nxL(NMAX_PG),nyL(NMAX_PG),nzL(NMAX_PG)   ! angular momentum of p.g.
     ! Small (3,4) components
     real(kind=dp) :: aS(NMAX_PG)
     real(kind=dp) :: xS(NMAX_PG),yS(NMAX_PG),zS(NMAX_PG)
     integer      :: nxS(NMAX_PG),nyS(NMAX_PG),nzS(NMAX_PG)
  end type primitive_gaussian

  type nonzeromat4legs
     integer :: a,b,c,d ! specify non-zero components
     ! assume that there 4 spinor legs a,b,c,d
     complex(kind=dp) :: val
  end type nonzeromat4legs

end module DefineTypes


!===================================================================
module Constants
! 121115 added gamma0*gammak = alpha
!===================================================================
  use Precision
  use DefineTypes
  implicit none
  real(kind=dp) :: PI
  real(kind=dp),parameter :: CCC= 137.035999679_dp  !light speed
  real(kind=dp),parameter :: Ze = -1._dp  ! factor for electron charge
  complex(kind=dp),parameter :: IU = (0._dp,1._dp)    ! imaginary unit

  complex(kind=dp) :: Gam0(4,4),Gam(3,4,4) ! 3 denotes gamma1,gamma2,gamma3
  complex(kind=dp) :: Gam5(4,4)   ! gamma_5 = -gamma^5
  complex(kind=dp) :: Sigma(3,4,4) ! 4x4 Sigma matrix
  complex(kind=dp) :: Gam0k(3,4,4) ! gamma^0 gamma^k = alpha
!  complex(kind=dp) :: GamJJ(4,4,4,4) ! sum_{k=1}^3 [gamma^0 gamma^k]_{alpha beta} [gamma^0 gamma^k]_{gamma delta}
!                                     ! used in calc_intthetajj_mat
  integer,parameter :: Ngamjj = 24 ! number of non-zero components in GamJJ
  type(nonzeromat4legs) :: GamJJ(Ngamjj)
  
  !-----------------------------------------
  ! for discretization of photon momentum
  !-----------------------------------------
!  real(kind=dp),parameter :: P0MAX = 40._dp
  real(kind=dp) :: P0MAX
!!$  integer,parameter :: NP0   = 2  ! # of norm of photon momentum (photon energy)
!!$  integer,parameter :: NPTH  = 3  ! # in theta direction (0 <= theta <= pi)
!!$  integer,parameter :: NPPHI = 4  ! # in phi direction (0<= phi < 2pi)
!!$  integer,parameter :: Nph = NP0*NPTH*NPPHI*2 ! j=1~Nph
  integer :: NP0     ! # of norm of photon momentum (photon energy)
  integer :: NPTH    ! # in theta direction (0 <= theta <= pi)
  integer :: NPPHI   ! # in phi direction (0<= phi < 2pi)
  integer :: Nph  ! j=1~Nph

  integer :: NOCC  ! 1~nocc electron orbitals are occupied.

  !--------------------------------------------
  ! for numerical integration of calJ and calL
  !--------------------------------------------
  integer,parameter :: NB_Lret = 2  ! 3D integration (dr)
  integer,parameter :: NB_Jret = 2  ! 3D integration (dr)
  integer,parameter :: NB_L2   = 2  ! 2D integration (ds)
  integer,parameter :: NB_J2   = 2  ! 2D integration (ds)

  !-----------------------------------------
  ! for time evolution
  !-----------------------------------------
  real(kind=dp) :: DeltaT ! time step 

  !  complex(kind=dp) :: ALPHA_COH(Nph)   ! <coherent|hat(a)|coherent>
  complex(kind=dp), allocatable :: ALPHA_COH(:)   ! <coherent|hat(a)|coherent>
  
end module Constants

!===================================================================
module IntegralStorage
!===================================================================
  use Precision
  implicit none

  logical :: there_is_twoele
!  character(LEN=80),parameter :: file_twoele = 'twoele.dat'  ! 200
  character(LEN=80) :: file_twoele  ! 200
  real(kind=dp),parameter :: TH_twoele = 1.e-16_dp  ! do not store two electron integral smaller than this value.
  integer :: N_twoele ! number of two-electron integral larger than TH_twoele

  character(LEN=80),parameter :: file_twoele_bin = 'twoele.bin'  ! 300

  logical :: there_is_intIJJ
  character(LEN=80) :: file_intIJJ  ! 204
  real(kind=dp),parameter :: TH_intIJJ = 1.e-16_dp  ! do not store IJJ integral smaller than this value.
  integer :: N_intIJJ ! number of IJJ integral larger than TH_intIJJ

!!$  logical :: there_is_intret
!!$  character(LEN=80),parameter :: file_intret = 'intret.dat'  ! 201
!!$  integer,parameter :: NRET = 7  ! number of time points to do integration for retardation

  logical :: there_is_nucele
  character(LEN=80),parameter :: file_nucele = 'nucele.dat'  ! 202

  logical :: there_is_twonuc
  character(LEN=80),parameter :: file_twonuc = 'twonuc.dat'  ! 203


end module IntegralStorage

!===================================================================
module NucBasis
!===================================================================
  use Precision
  implicit none

  integer,parameter :: Z_N = 1 ! nuclear charge for H
  integer,parameter :: N_PGN = 15 ! number of basis set function for expansion function for nuclear field.
  real(kind=dp) :: vecR_N(3,N_PGN)  ! center position of 1s gaussian
  integer,parameter :: NBS_PHI = 7 ! number of expansion functions for nuclear field
  integer :: NBS_N  ! =NBS_PHI for boson, =2xNBS_PHI for fermion
  real(kind=dp) :: L_well ! size of the single well potential
  real(kind=dp) :: L_pg ! range for p.g. (larger than L_well)

  ! parameters below are set in set_NucBasis in sub_int_nuc.f90
  integer :: NNUC ! number of nucleus
  real(kind=dp) :: m_N  ! nuclear mass  
  real(kind=dp) :: alpha_N  ! exponent for nucleus
!  real(kind=dp) :: R_N 
  integer :: NUCTYPE ! boson(0) of fermion(1)
  complex(kind=dp) :: c_nuc(NBS_PHI,N_PGN) ! coefficients to build plane waves from primitive gaussians

end module NucBasis


!===================================================================
module System_Medium
!===================================================================
  use Precision
  implicit none

  real(kind=dp) :: L_Med ! size of medium (M) (assumed to be cube)
  real(kind=dp) :: R_Sys ! size of system (A) (assumed to be ball)

end module System_Medium



!===================================================================
module Interface_Mod
!===================================================================
  Interface

     subroutine setQmat_twoele_nz(twoele_Qmat_nz)
       use DefineTypes
       type(nonzeromat4legs),allocatable,intent(out) :: twoele_Qmat_nz(:)
     end subroutine setQmat_twoele_nz

     subroutine setQmat_intIJJ_nz(intIJJ_Qmat_nz)
       use DefineTypes
       type(nonzeromat4legs),allocatable :: intIJJ_Qmat_nz(:)
     end subroutine setQmat_intIJJ_nz

  end interface
end module Interface_mod
