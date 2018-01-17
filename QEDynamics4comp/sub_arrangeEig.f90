!===============================================================
!  subroutine arrangeEigenvalueE  version 2012/3/31
!  subroutine arrangeEigenvalueP  version 2012/3/31
!  subroutine arrangeEigenvalue   version 2010/12/17
!===============================================================
!===============================================================
subroutine arrangeEigenvalueE
!===============================================================

  use DiracOutput 

  implicit none

  double precision tmpe
  integer k,kk,i,j,orbnew(NBS_E),tmporb
  complex(kind(0d0)) tmpc_La(NBS_L,NBS_E),tmpc_Lb(NBS_L,NBS_E)
  complex(kind(0d0)) tmpc_Sa(NBS_S,NBS_E),tmpc_Sb(NBS_S,NBS_E)
!  complex(kind(0d0)) tmpc_La(NBS_L),tmpc_Lb(NBS_L)
!  complex(kind(0d0)) tmpc_Sa(NBS_S),tmpc_Sb(NBS_S)

  do kk=1,NBS_E
   orbnew(kk)=kk
   do i=1,NBS_L
    tmpc_La(i,kk)=c_La(i,kk)
    tmpc_Lb(i,kk)=c_Lb(i,kk)
   end do
   do j=1,NBS_S
    tmpc_Sa(j,kk)=c_Sa(j,kk)
    tmpc_Sb(j,kk)=c_Sb(j,kk)
   end do
  end do

  do k=1,NBS_E-1
   do kk=k+1,NBS_E
    if(e_eig(k).gt.e_eig(kk)) then
     tmpe=e_eig(k)
     e_eig(k)=e_eig(kk)
     e_eig(kk)=tmpe
     tmporb=orbnew(k)
     orbnew(k)=orbnew(kk)
     orbnew(kk)=tmporb
!!$     do i=1,NBS_L
!!$      tmpc_La(i)=c_La(i,k)
!!$      tmpc_Lb(i)=c_Lb(i,k)
!!$      c_La(i,k)=c_La(i,kk)
!!$      c_Lb(i,k)=c_Lb(i,kk)
!!$      c_La(i,kk)=tmpc_La(i)
!!$      c_Lb(i,kk)=tmpc_Lb(i)
!!$     end do
!!$     do i=1,NBS_S
!!$      tmpc_Sa(i)=c_Sa(i,k)
!!$      tmpc_Sb(i)=c_Sb(i,k)
!!$      c_Sa(i,k)=c_Sa(i,kk)
!!$      c_Sb(i,k)=c_Sb(i,kk)
!!$      c_Sa(i,kk)=tmpc_Sa(i)
!!$      c_Sb(i,kk)=tmpc_Sb(i)
!!$     end do
    end if
   end do
  end do

   write(*,*) 'k,e_eig(k)'
  do k=1,NBS_E
   write(*,*) k,e_eig(k)
  end do

  do k=1,NBS_E
   do kk=1,NBS_E
    if(k.eq.orbnew(kk)) then
     do i=1,NBS_L
      c_La(i,kk)=tmpc_La(i,k)
      c_Lb(i,kk)=tmpc_Lb(i,k)
     end do
     do j=1,NBS_S
      c_Sa(j,kk)=tmpc_Sa(j,k)
      c_Sb(j,kk)=tmpc_Sb(j,k)
     end do
    end if
   end do
  end do

end subroutine arrangeEigenvalueE

!===============================================================
subroutine arrangeEigenvalueP
!===============================================================

  use DiracOutput 

  implicit none

  double precision tmpp
  integer k,kk,i,j,orbnew(NBS_P),tmporb
  complex(kind(0d0)) tmpd_La(NBS_L,NBS_P),tmpd_Lb(NBS_L,NBS_P)
  complex(kind(0d0)) tmpd_Sa(NBS_S,NBS_P),tmpd_Sb(NBS_S,NBS_P)
!  complex(kind(0d0)) tmpd_La(NBS_L),tmpd_Lb(NBS_L)
!  complex(kind(0d0)) tmpd_Sa(NBS_S),tmpd_Sb(NBS_S)

  do kk=1,NBS_P
   orbnew(kk)=kk
   do i=1,NBS_L
    tmpd_La(i,kk)=d_La(i,kk)
    tmpd_Lb(i,kk)=d_Lb(i,kk)
   end do
   do j=1,NBS_S
    tmpd_Sa(j,kk)=d_Sa(j,kk)
    tmpd_Sb(j,kk)=d_Sb(j,kk)
   end do
  end do

  do k=1,NBS_P-1
   do kk=k+1,NBS_P
    if(p_eig(k).gt.p_eig(kk)) then
     tmpp=p_eig(k)
     p_eig(k)=p_eig(kk)
     p_eig(kk)=tmpp
     tmporb=orbnew(k)
     orbnew(k)=orbnew(kk)
     orbnew(kk)=tmporb
!!$     do i=1,NBS_L
!!$      tmpd_La(i)=d_La(i,k)
!!$      tmpd_Lb(i)=d_Lb(i,k)
!!$      d_La(i,k)=d_La(i,kk)
!!$      d_Lb(i,k)=d_Lb(i,kk)
!!$      d_La(i,kk)=tmpd_La(i)
!!$      d_Lb(i,kk)=tmpd_Lb(i)
!!$     end do
!!$     do i=1,NBS_S
!!$      tmpd_Sa(i)=d_Sa(i,k)
!!$      tmpd_Sb(i)=d_Sb(i,k)
!!$      d_Sa(i,k)=d_Sa(i,kk)
!!$      d_Sb(i,k)=d_Sb(i,kk)
!!$      d_Sa(i,kk)=tmpd_Sa(i)
!!$      d_Sb(i,kk)=tmpd_Sb(i)
!!$     end do
    end if
   end do
  end do
 
   write(*,*) 'k,p_eig(k)'
  do k=1,NBS_P
   write(*,*) k,p_eig(k)
  end do

  do k=1,NBS_P
   do kk=1,NBS_P
    if(k.eq.orbnew(kk)) then
     do i=1,NBS_L
      d_La(i,kk)=tmpd_La(i,k)
      d_Lb(i,kk)=tmpd_Lb(i,k)
     end do
     do j=1,NBS_S
      d_Sa(j,kk)=tmpd_Sa(j,k)
      d_Sb(j,kk)=tmpd_Sb(j,k)
     end do
    end if
   end do
  end do

end subroutine arrangeEigenvalueP


!!$!===============================================================
!!$subroutine arrangeEigenvalue(e_eig,p_eig)
!!$!===============================================================
!!$
!!$  use DiracOutput 
!!$
!!$  implicit none
!!$
!!$  double precision,intent(inout):: e_eig(NBS)
!!$  double precision tmpe
!!$  double precision,intent(inout):: p_eig(NBS)
!!$  double precision tmpp
!!$  integer k,kk,i,j,orbnew(NBS),tmporb
!!$!  complex(kind(0d0)) tmpc_La(NBS_L,NBS),tmpc_Lb(NBS_L,NBS)
!!$!  complex(kind(0d0)) tmpc_Sa(NBS_S,NBS),tmpc_Sb(NBS_S,NBS)
!!$  complex(kind(0d0)) tmpc_La(NBS_L),tmpc_Lb(NBS_L)
!!$  complex(kind(0d0)) tmpc_Sa(NBS_S),tmpc_Sb(NBS_S)
!!$  complex(kind(0d0)) tmpd_La(NBS_L),tmpd_Lb(NBS_L)
!!$  complex(kind(0d0)) tmpd_Sa(NBS_S),tmpd_Sb(NBS_S)
!!$
!!$  do kk=1,NBS
!!$   orbnew(kk)=kk
!!$!!$   do i=1,NBS_L
!!$!!$    tmpc_La(i,kk)=c_La(i,kk)
!!$!!$    tmpc_Lb(i,kk)=c_Lb(i,kk)
!!$!!$   end do
!!$!!$   do j=1,NBS_S
!!$!!$    tmpc_Sa(j,kk)=c_Sa(j,kk)
!!$!!$    tmpc_Sb(j,kk)=c_Sb(j,kk)
!!$!!$   end do
!!$  end do
!!$
!!$  do k=1,NBS-1
!!$   do kk=k+1,NBS
!!$    if(e_eig(k).gt.e_eig(kk)) then
!!$     tmpe=e_eig(k)
!!$     e_eig(k)=e_eig(kk)
!!$     e_eig(kk)=tmpe
!!$     tmpp=p_eig(k)
!!$     p_eig(k)=p_eig(kk)
!!$     p_eig(kk)=tmpp
!!$     tmporb=orbnew(k)
!!$     orbnew(k)=orbnew(kk)
!!$     orbnew(kk)=tmporb
!!$     do i=1,NBS_L
!!$      tmpc_La(i)=c_La(i,k)
!!$      tmpc_Lb(i)=c_Lb(i,k)
!!$      c_La(i,k)=c_La(i,kk)
!!$      c_Lb(i,k)=c_Lb(i,kk)
!!$      c_La(i,kk)=tmpc_La(i)
!!$      c_Lb(i,kk)=tmpc_Lb(i)
!!$      tmpd_La(i)=d_La(i,k)
!!$      tmpd_Lb(i)=d_Lb(i,k)
!!$      d_La(i,k)=d_La(i,kk)
!!$      d_Lb(i,k)=d_Lb(i,kk)
!!$      d_La(i,kk)=tmpd_La(i)
!!$      d_Lb(i,kk)=tmpd_Lb(i)
!!$     end do
!!$     do i=1,NBS_S
!!$      tmpc_Sa(i)=c_Sa(i,k)
!!$      tmpc_Sb(i)=c_Sb(i,k)
!!$      c_Sa(i,k)=c_Sa(i,kk)
!!$      c_Sb(i,k)=c_Sb(i,kk)
!!$      c_Sa(i,kk)=tmpc_Sa(i)
!!$      c_Sb(i,kk)=tmpc_Sb(i)
!!$      tmpd_Sa(i)=d_Sa(i,k)
!!$      tmpd_Sb(i)=d_Sb(i,k)
!!$      d_Sa(i,k)=d_Sa(i,kk)
!!$      d_Sb(i,k)=d_Sb(i,kk)
!!$      d_Sa(i,kk)=tmpd_Sa(i)
!!$      d_Sb(i,kk)=tmpd_Sb(i)
!!$     end do
!!$    end if
!!$   end do
!!$  end do
!!$
!!$   write(*,*) 'k,e_eig(k),p_eig(k)'
!!$  do k=1,NBS
!!$   write(*,*) k,e_eig(k),p_eig(k)
!!$  end do
!!$!!$  do k=1,NBS
!!$!!$   do kk=1,NBS
!!$!!$    if(k.eq.orbnew(kk)) then
!!$!!$     do i=1,NBS_L
!!$!!$      c_La(i,k)=tmpc_La(i,orbnew(kk))
!!$!!$      c_Lb(i,k)=tmpc_Lb(i,orbnew(kk))
!!$!!$     end do
!!$!!$     do j=1,NBS_S
!!$!!$      c_Sa(j,k)=tmpc_Sa(j,orbnew(kk))
!!$!!$      c_Sb(j,k)=tmpc_Sb(j,orbnew(kk))
!!$!!$     end do
!!$!!$    end if
!!$!!$   end do
!!$!!$  end do
!!$
!!$end subroutine arrangeEigenvalue
