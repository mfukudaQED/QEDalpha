!==============================================================================
! subroutines for eigenvalues and eigenvectors real symmetric matrix.
!
! 2012.9.13 
!  -copied from H2p_5.f90 and diag.f90 in 02_H2p (other project)
!  -modified so that eig(1)<eig(2)<eig(3) 
! 
!
! f90 version of tqli, tred2, pythag
!==============================================================================



!======================================================================
!subroutine calc_eigen(Ts,d,a)
subroutine calc_eigen(Ts,eig,vec)
!
! For given 3x3 matrix "Ts", calculate eigenvalues "d" and eigenvectors "a".
! d's are sorted in order of d(1)>d(2)>d(3)
! 
! [a(1,i), a(2,i), a(3,i)] is the eigenvector for i-th eigenvalue d(i)
!
!======================================================================
  implicit none

  real(kind=8), intent(in) :: Ts(3,3)
  REAL(kind=8), intent(out) :: eig(3),vec(3,3)
  REAL(kind=8) :: e(3),d(3),a(3,3)
  integer :: i,j,ii

  do i=1,3
     do j=1,3
        a(i,j) = Ts(i,j)
!        write(*,*) a(i,j)
     end do
  end do

  call tred2(a,3,3,d,e)
  call tqli(d,e,3,3,a)
  !do ii=1,3
  !write(*,'(4es16.6)') d(ii),a(1,ii),a(2,ii),a(3,ii) 
  !end do
  
  ! sort w.r.t. eigenvalue (In the order of d(1) > d(2) > d(3) )
  call eigsrt(d,a,3,3)  

 !            do ii=1,3
 !               write(*,'(4es16.6)') d(ii),a(1,ii),a(2,ii),a(3,ii) 
 !            end do
  !do ii=1,3
  !     write(*,'(4es16.6)',advance='no') d(ii),a(1,ii),a(2,ii),a(3,ii) 
  !end do
  !  write(*,*)

  ! re-sort
  eig(3) = d(1)
  eig(2) = d(2)
  eig(1) = d(3)
  do j=1,3
     vec(j,3) = a(j,1)
     vec(j,2) = a(j,2)
     vec(j,1) = a(j,3)
  end do

  return
end subroutine calc_eigen

!!$PROGRAM xtqli
!!$  !C     driver for routine tqli
!!$  INTEGER,parameter :: NP=10
!!$  REAL(kind=8),parameter :: TINY=1.d-6
!!$  !PARAMETER(NP=10,TINY=1.0e-6)
!!$  INTEGER :: i,j,k
!!$  REAL(kind=8) :: a(NP,NP),c(NP,NP),d(NP),e(NP),f(NP)
!!$  DATA c/5.0,4.3,3.0,2.0,1.0,0.0,-1.0,-2.0,-3.0,-4.0, &
!!$       &     4.3,5.1,4.0,3.0,2.0,1.0,0.0,-1.0,-2.0,-3.0, &
!!$       &     3.0,4.0,5.0,4.0,3.0,2.0,1.0,0.0,-1.0,-2.0, &
!!$       &     2.0,3.0,4.0,5.0,4.0,3.0,2.0,1.0,0.0,-1.0, &
!!$       &     1.0,2.0,3.0,4.0,5.0,4.0,3.0,2.0,1.0,0.0, &
!!$       &     0.0,1.0,2.0,3.0,4.0,5.0,4.0,3.0,2.0,1.0, &
!!$       &     -1.0,0.0,1.0,2.0,3.0,4.0,5.0,4.0,3.0,2.0, &
!!$       &     -2.0,-1.0,0.0,1.0,2.0,3.0,4.0,5.0,4.0,3.0, &
!!$       &     -3.0,-2.0,-1.0,0.0,1.0,2.0,3.0,4.0,5.0,4.0, &
!!$       &     -4.0,-3.0,-2.0,-1.0,0.0,1.0,2.0,3.0,4.0,5.0/
!!$  do i=1,NP
!!$     do j=1,NP
!!$        a(i,j)=c(i,j)
!!$     end do
!!$  end do
!!$  call tred2(a,NP,NP,d,e)
!!$  call tqli(d,e,NP,NP,a)
!!$  write(*,'(/1x,a)') 'Eigenvectors for a real symmetric matrix'
!!$  do i=1,NP
!!$     do j=1,NP
!!$        f(j)=0.0
!!$        do k=1,NP
!!$           f(j)=f(j)+c(j,k)*a(k,i)
!!$        end do
!!$     end do
!!$     write(*,'(/1x,a,i3,a,f10.6)') 'Eigenvalue',i,' =',d(i)
!!$     write(*,'(/1x,t7,a,t17,a,t31,a)') 'Vector','Mtrx*Vect.','Ratio'
!!$     do j=1,NP
!!$        if (abs(a(j,i)).lt.TINY) then
!!$           write(*,'(1x,2f12.6,a12)') a(j,i),f(j),'div. by 0'
!!$        else
!!$           write(*,'(1x,2f12.6,e14.6)') a(j,i),f(j),f(j)/a(j,i)
!!$        endif
!!$     end do
!!$     write(*,'(/1x,a)') 'press ENTER to continue...'
!!$     read(*,*)
!!$  end do
!!$END PROGRAM xtqli
!!$!C  (C) Copr. 1986-92 Numerical Recipes Software z!0(0.

!=====================================
SUBROUTINE tqli(d,e,n,np,z)
!=====================================
  INTEGER :: n,np
  REAL(kind=8) ::  d(np),e(np),z(np,np)
!CU    USES pythag
  INTEGER :: i,iter,k,l,m
  REAL(kind=8) :: b,c,dd,f,g,p,r,s,pythag
  do i=2,n
     e(i-1)=e(i)
  end do
  e(n)=0.d0
  do l=1,n
     iter=0
1    do m=l,n-1
        dd=abs(d(m))+abs(d(m+1))
        if (abs(e(m))+dd.eq.dd) goto 2
     end do
     m=n
2    if(m.ne.l)then
!        if(iter.eq.30)pause 'too many iterations in tqli'
        if(iter.eq.30)stop 'too many iterations in tqli'
        iter=iter+1
        g=(d(l+1)-d(l))/(2.*e(l))
        r=pythag(g,1.d0)
        g=d(m)-d(l)+e(l)/(g+sign(r,g))
        s=1.d0
        c=1.d0
        p=0.d0
        do i=m-1,l,-1
           f=s*e(i)
           b=c*e(i)
           r=pythag(f,g)
           e(i+1)=r
           if(r.eq.0.d0)then
              d(i+1)=d(i+1)-p
              e(m)=0.d0
              goto 1
           endif
           s=f/r
           c=g/r
           g=d(i+1)-p
           r=(d(i)-g)*s+2.d0*c*b
           p=s*r
           d(i+1)=g+p
           g=c*r-b
!C     Omit lines from here ...
           do k=1,n
              f=z(k,i+1)
              z(k,i+1)=s*z(k,i)+c*f
              z(k,i)=c*z(k,i)-s*f
           end do
!C     ... to here when finding only eigenvalues.
        end do
        d(l)=d(l)-p
        e(l)=g
        e(m)=0.d0
        goto 1
     endif
  end do
  return
END SUBROUTINE tqli
!C  (C) Copr. 1986-92 Numerical Recipes Software z!0(0.

!=====================================
SUBROUTINE tred2(a,n,np,d,e)
!=====================================
  INTEGER :: n,np
  REAL(kind=8) :: a(np,np),d(np),e(np)
  INTEGER :: i,j,k,l
  REAL(kind=8) :: f,g,h,hh,scale
  do i=n,2,-1
     l=i-1
     h=0.d0
     scale=0.d0
     if(l.gt.1)then
        do k=1,l
           scale=scale+abs(a(i,k))
        end do
        if(scale.eq.0.d0)then
           e(i)=a(i,l)
        else
           do k=1,l
              a(i,k)=a(i,k)/scale
              h=h+a(i,k)**2
           end do
           f=a(i,l)
           g=-sign(sqrt(h),f)
           e(i)=scale*g
           h=h-f*g
           a(i,l)=f-g
           f=0.d0
           do j=1,l
!C     Omit following line if finding only eigenvalues
              a(j,i)=a(i,j)/h
              g=0.d0
              do k=1,j
                 g=g+a(j,k)*a(i,k)
              end do
              do k=j+1,l
                 g=g+a(k,j)*a(i,k)
              end do
              e(j)=g/h
              f=f+e(j)*a(i,j)
           end do
           hh=f/(h+h)
           do j=1,l
              f=a(i,j)
              g=e(j)-hh*f
              e(j)=g
              do k=1,j
                 a(j,k)=a(j,k)-f*e(k)-g*a(i,k)
              end do
           end do
        endif
     else
        e(i)=a(i,l)
     endif
     d(i)=h
  end do
!C     Omit following line if finding only eigenvalues.
  d(1)=0.d0
  e(1)=0.d0
  do i=1,n
!C     Delete lines from here ...
     l=i-1
     if(d(i).ne.0.d0)then
        do j=1,l
           g=0.d0
           do k=1,l
              g=g+a(i,k)*a(k,j)
           end do
           do k=1,l
              a(k,j)=a(k,j)-g*a(k,i)
           end do
        end do
     endif
!C     ... to here when finding only eigenvalues.
     d(i)=a(i,i)
!C     Also delete lines from here ...
     a(i,i)=1.d0
     do j=1,l
        a(i,j)=0.d0
        a(j,i)=0.d0
     end do
!C     ... to here when finding only eigenvalues.
  end do
  return
END SUBROUTINE tred2
!C  (C) Copr. 1986-92 Numerical Recipes Software z!0(0.

!=====================================
FUNCTION pythag(a,b)
!=====================================
  REAL(kind=8) :: a,b,pythag
  REAL(kind=8) :: absa,absb
  absa=abs(a)
  absb=abs(b)
  if(absa.gt.absb)then
     pythag=absa*sqrt(1.d0+(absb/absa)**2)
  else
     if(absb.eq.0.d0)then
        pythag=0.d0
     else
        pythag=absb*sqrt(1.d0+(absa/absb)**2)
     endif
  endif
  return
END FUNCTION pythag
!C  (C) Copr. 1986-92 Numerical Recipes Software z!0(0.

!===============================
SUBROUTINE eigsrt(d,v,n,np)
!===============================
  INTEGER n,np
  REAL(kind=8) :: d(np),v(np,np)
  INTEGER :: i,j,k
  REAL(kind=8) :: p
  do i=1,n-1
     k=i
     p=d(i)
     do j=i+1,n
        if(d(j).ge.p)then
           k=j
           p=d(j)
        endif
     end do
     if(k.ne.i)then
        d(k)=d(i)
        d(i)=p
        do j=1,n
           p=v(j,i)
           v(j,i)=v(j,k)
           v(j,k)=p
        end do
     endif
  end do
  return
END SUBROUTINE eigsrt
!C  (C) Copr. 1986-92 Numerical Recipes Software z!0(0.
