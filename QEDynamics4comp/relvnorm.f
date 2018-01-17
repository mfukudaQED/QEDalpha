      subroutine normal(aa,bb,n,NBS0,NBS)
C
      implicit real*8 (a-h,o-z)
C
      real*8 aa(NBS0,4), bb(NBS0,4)
C
      integer n(NBS0,4)
C
C
      parameter ( two  = 2.0D0 )
      parameter ( thre = 3.0D0 )
      parameter ( four = 4.0D0 )
C
      parameter ( e34 = 0.75D0 )
C
      pi   = 3.14159265358979323D0
C
C***************************************************************
C
      do k0=1,4
      do j=1,NBS
        if (n(j,k0) .eq. 0) then
         b = ( two*aa(j,k0) / pi )**e34
        else if (n(j,k0) .le. 3) then
         b = two*dsqrt(aa(j,k0)) * ( two*aa(j,k0) / pi )**e34
        else if ((n(j,k0) .eq. 4) .or.
     &           (n(j,k0) .eq. 7) .or.
     &           (n(j,k0) .eq. 9)) then
c         b = four*aa(j,k0) * ( two*aa(j,k0) / pi )**e34 / dsqrt(thre)
c editted by hara/ ref. DIRAC10/dft/basis_info.F90 line:106-116
         b = four*aa(j,k0) * ( two*aa(j,k0) / pi )**e34
c
        else if ((n(j,k0) .eq. 5) .or.
     &           (n(j,k0) .eq. 6) .or.
     &           (n(j,k0) .eq. 8)) then
         b = four*aa(j,k0) * ( two*aa(j,k0) / pi )**e34
        else if ((n(j,k0) .eq. 10) .or.
     &           (n(j,k0) .eq. 16) .or.
     &           (n(j,k0) .eq. 19)) then
c         b = 8.D0*aa(j,k0)*dsqrt(aa(j,k0)/15.D0)*(two*aa(j,k0)/pi)**e34
         b = 8.D0*aa(j,k0)*dsqrt(aa(j,k0))*(two*aa(j,k0)/pi)**e34
        else if ((n(j,k0) .eq. 11) .or.
     &           (n(j,k0) .eq. 12) .or.
     &           (n(j,k0) .eq. 13) .or.
     &           (n(j,k0) .eq. 15) .or.
     &           (n(j,k0) .eq. 17) .or.
     &           (n(j,k0) .eq. 18)) then
c         b = 8.D0*aa(j,k0)*dsqrt(aa(j,k0)/thre)*(two*aa(j,k0)/pi)**e34
         b = 8.D0*aa(j,k0)*dsqrt(aa(j,k0))*(two*aa(j,k0)/pi)**e34
        else if (n(j,k0) .eq. 14) then
         b = 8.0D0*aa(j,k0)*dsqrt(aa(j,k0))*(two*aa(j,k0)/pi)**e34
        else if ((n(j,k0) .eq. 20) .or.
     &           (n(j,k0) .eq. 30) .or.
     &           (n(j,k0) .eq. 34)) then
c         b = 16.D0*aa(j,k0)*aa(j,k0)/dsqrt(105.D0)
c     &            *(two*aa(j,k0)/pi)**e34
         b = 16.D0*aa(j,k0)*aa(j,k0)
     &            *(two*aa(j,k0)/pi)**e34
        else if ((n(j,k0) .eq. 21) .or.
     &           (n(j,k0) .eq. 22) .or.
     &           (n(j,k0) .eq. 26) .or.
     &           (n(j,k0) .eq. 29) .or.
     &           (n(j,k0) .eq. 31) .or.
     &           (n(j,k0) .eq. 33)) then
c         b = 16.D0*aa(j,k0)*aa(j,k0)/dsqrt(15.D0)
c     &            *(two*aa(j,k0)/pi)**e34
         b = 16.D0*aa(j,k0)*aa(j,k0)
     &            *(two*aa(j,k0)/pi)**e34
        else if ((n(j,k0) .eq. 23) .or.
     &           (n(j,k0) .eq. 25) .or.
     &           (n(j,k0) .eq. 32)) then
c         b = 16.D0*aa(j,k0)*aa(j,k0)/thre*(two*aa(j,k0)/pi)**e34
         b = 16.D0*aa(j,k0)*aa(j,k0)*(two*aa(j,k0)/pi)**e34
        else if ((n(j,k0) .eq. 24) .or.
     &           (n(j,k0) .eq. 27) .or.
     &           (n(j,k0) .eq. 28)) then
c         b = 16.D0*aa(j,k0)*aa(j,k0)/dsqrt(thre)
c     &            *(two*aa(j,k0)/pi)**e34
         b = 16.D0*aa(j,k0)*aa(j,k0)
     &            *(two*aa(j,k0)/pi)**e34
        else if ((n(j,k0) .eq. 35) .or.
     &           (n(j,k0) .eq. 50) .or.
     &           (n(j,k0) .eq. 55)) then
c         b = 32.D0*aa(j,k0)*aa(j,k0)*dsqrt(aa(j,k0)/105.D0)
c     &            *(two*aa(j,k0)/pi)**e34/thre
         b = 32.D0*aa(j,k0)*aa(j,k0)*dsqrt(aa(j,k0))
     &            *(two*aa(j,k0)/pi)**e34
        else if ((n(j,k0) .eq. 36) .or.
     &           (n(j,k0) .eq. 37) .or.
     &           (n(j,k0) .eq. 45) .or.
     &           (n(j,k0) .eq. 49) .or.
     &           (n(j,k0) .eq. 51) .or.
     &           (n(j,k0) .eq. 54)) then
c         b = 32.D0*aa(j,k0)*aa(j,k0)*dsqrt(aa(j,k0)/105.D0)
c     &            *(two*aa(j,k0)/pi)**e34
         b = 32.D0*aa(j,k0)*aa(j,k0)*dsqrt(aa(j,k0))
     &            *(two*aa(j,k0)/pi)**e34
        else if ((n(j,k0) .eq. 38) .or.
     &           (n(j,k0) .eq. 40) .or.
     &           (n(j,k0) .eq. 41) .or.
     &           (n(j,k0) .eq. 44) .or.
     &           (n(j,k0) .eq. 52) .or.
     &           (n(j,k0) .eq. 53)) then
c         b = 32.D0*aa(j,k0)*aa(j,k0)*dsqrt(aa(j,k0)/5.D0)
c     &            *(two*aa(j,k0)/pi)**e34/thre
         b = 32.D0*aa(j,k0)*aa(j,k0)*dsqrt(aa(j,k0))
     &            *(two*aa(j,k0)/pi)**e34
        else if ((n(j,k0) .eq. 39) .or.
     &           (n(j,k0) .eq. 46) .or.
     &           (n(j,k0) .eq. 48)) then
c         b = 32.D0*aa(j,k0)*aa(j,k0)*dsqrt(aa(j,k0)/15.D0)
c     &            *(two*aa(j,k0)/pi)**e34
         b = 32.D0*aa(j,k0)*aa(j,k0)*dsqrt(aa(j,k0))
     &            *(two*aa(j,k0)/pi)**e34
        else if ((n(j,k0) .eq. 42) .or.
     &           (n(j,k0) .eq. 43) .or.
     &           (n(j,k0) .eq. 47)) then
c         b = 32.D0*aa(j,k0)*aa(j,k0)*dsqrt(aa(j,k0))
c     &            *(two*aa(j,k0)/pi)**e34/thre
         b = 32.D0*aa(j,k0)*aa(j,k0)*dsqrt(aa(j,k0))
     &            *(two*aa(j,k0)/pi)**e34
        else
         b = 0.0D0
        end if
        bb(j,k0) = b * bb(j,k0)
      end do
      end do
C
      return
      end


ccccccccccccccccccccccccccccccccccccc
      subroutine normals(tmpa,tmpb,n)
C
      implicit real*8 (a-h,o-z)
C
      real*8 tmpa,tmpb
C
      integer n
C
C
      parameter ( two  = 2.0D0 )
      parameter ( thre = 3.0D0 )
      parameter ( four = 4.0D0 )
C
      parameter ( e34 = 0.75D0 )
C
      pi   = 3.14159265358979323D0
C
C***************************************************************
C
        if (n.eq. 0) then
          bnor = ( two*tmpa / pi )**e34
        else if (n.le.3) then
          bnor = two*dsqrt(tmpa) * ( two*tmpa / pi )**e34
        else if ((n.eq. 4) .or. (n.eq. 7) .or. (n.eq. 9)) then
          bnor = four*tmpa * ( two*tmpa / pi )**e34 / dsqrt(thre)
        else if ((n.eq. 5) .or. (n.eq. 6) .or. (n.eq. 8)) then
          bnor = four*tmpa * ( two*tmpa / pi )**e34
        else if ((n.eq. 10) .or. (n.eq. 16) .or. (n.eq. 19)) then
          bnor = 8.D0*tmpa*dsqrt(tmpa/15.D0)*(two*tmpa/pi)**e34
        else if ((n.eq. 11) .or. (n.eq. 12) .or. (n.eq. 13) .or.
     &           (n.eq. 15) .or. (n.eq. 17) .or. (n.eq. 18)) then
          bnor = 8.D0*tmpa*dsqrt(tmpa/thre)*(two*tmpa/pi)**e34
        else if (n.eq. 14) then
          bnor = 8.0D0*tmpa*dsqrt(tmpa)*(two*tmpa/pi)**e34
        else if ((n.eq. 20) .or. (n.eq. 30) .or. (n.eq. 34)) then
          bnor = 16.D0*tmpa*tmpa/dsqrt(105.D0) *(two*tmpa/pi)**e34
        else if ((n.eq. 21) .or. (n.eq. 22) .or. (n.eq. 26) .or.
     &           (n.eq. 29) .or. (n.eq. 31) .or. (n.eq. 33)) then
          bnor = 16.D0*tmpa*tmpa/dsqrt(15.D0) *(two*tmpa/pi)**e34
        else if ((n.eq. 23) .or. (n.eq. 25) .or. (n.eq. 32)) then
          bnor = 16.D0*tmpa*tmpa/thre*(two*tmpa/pi)**e34
        else if ((n.eq. 24) .or. (n.eq. 27) .or. (n.eq. 28)) then
          bnor = 16.D0*tmpa*tmpa/dsqrt(thre) *(two*tmpa/pi)**e34
        else if ((n.eq. 35) .or. (n.eq. 50) .or. (n.eq. 55)) then
          bnor = 32.D0*tmpa*tmpa*dsqrt(tmpa/105.D0) 
     &            *(two*tmpa/pi)**e34/thre
        else if ((n.eq. 36) .or. (n.eq. 37) .or. (n.eq. 45) .or.
     &           (n.eq. 49) .or. (n.eq. 51) .or. (n.eq. 54)) then
          bnor = 32.D0*tmpa*tmpa*dsqrt(tmpa/105.D0)
     &            *(two*tmpa/pi)**e34
        else if ((n.eq. 38) .or. (n.eq. 40) .or. (n.eq. 41) .or.
     &           (n.eq. 44) .or. (n.eq. 52) .or. (n.eq. 53)) then
          bnor = 32.D0*tmpa*tmpa*dsqrt(tmpa/5.D0)
     &            *(two*tmpa/pi)**e34/thre
        else if ((n.eq. 39) .or. (n.eq. 46) .or. (n.eq. 48)) then
          bnor = 32.D0*tmpa*tmpa*dsqrt(tmpa/15.D0)
     &            *(two*tmpa/pi)**e34
        else if ((n.eq. 42) .or. (n.eq. 43) .or. (n.eq. 47)) then
          bnor = 32.D0*tmpa*tmpa*dsqrt(tmpa)
     &            *(two*tmpa/pi)**e34/thre
        else
          bnor = 0.0D0
        end if
        tmpb = bnor * tmpb
C
      return
      end


      subroutine normal2(nrmca,nrmcb,n,ncon,resfac)
C
      implicit real*8 (a-h,o-z)
C
      real*8 sumres, nrmca(100), nrmcb(100)
C
      integer n
C
c      parameter ( e34 = 0.75D0 )
C
      pi   = 3.14159265358979323D0
C
C***************************************************************
C
      sumres = 0.0d0
      if (n.eq.0) then
        do i=1,ncon
c        write(*,*) 'nrmca i',i, nrmca(i), 'nrmcb i',i, nrmcb(i) 
          do j=1,ncon
            sumres = sumres + nrmcb(i) * nrmcb(j) * pi**1.5
     &               / ( nrmca(i) + nrmca(j) )**1.5
          end do
        end do
      else if (n.eq.1) then 
        do i=1,ncon
          do j=1,ncon
            sumres = sumres + 0.5d0 * nrmcb(i) * nrmcb(j) * pi**1.5
     &               / ( nrmca(i) + nrmca(j) )**2.5 
          end do
        end do
      else if (n.eq.2) then
        do i=1,ncon
          do j=1,ncon
            sumres = sumres + 0.75d0 * nrmcb(i) * nrmcb(j) * pi**1.5
     &               / ( nrmca(i) + nrmca(j) )**3.5 
          end do
        end do
      else
        write (*,*) "ERROR : normal2, n = ", n
        stop
      end if
c      write(*,*) 'sumres = ',sumres
      resfac = 1/dsqrt(sumres)
      write(*,*) 'resfac = ', resfac

      return
      end


