! Last Change:28-May-2012.
!============================================================================
subroutine MRDFTCI
!============================================================================
  use DiracOutput
  use CI 
!      implicit real*8 (a-h,o-z)
      implicit none

      double precision, allocatable :: aa(:,:)
      double precision, allocatable :: xx(:,:),yy(:,:),zz(:,:),fracn(:)
      double precision, allocatable :: enorb(:),occdet(:)
      complex(kind(0d0)), allocatable :: a(:,:,:),f(:,:)
      complex(kind(0d0)), allocatable :: f00(:,:),fx(:,:),fy(:,:),fz(:,:)
      complex(kind(0d0)), allocatable :: fx00(:,:),fy00(:,:),fz00(:,:),g00(:,:)
      complex(kind(0d0)), allocatable :: fxx00(:,:),fyy00(:,:),fzz00(:,:)
      complex(kind(0d0)), allocatable :: fxx(:,:),fyy(:,:),fzz(:,:)
!      complex*16, allocatable :: fxy(:,:),fxz(:,:),fyz(:,:)
      integer, allocatable :: ras(:,:),n(:,:)
      double precision stress(9)
      double precision ccc,w1,w2,w3,x,y,z,detdns,cofpdt
      double precision xxx,yyy,zzz
      double precision rr
      double precision detldns,detsdns,ldnsty,sdnsty
      double precision zetax,zetay,zetaz,etaux,etauy,etauz,spina,spinb,spinc
      double precision spinx,spiny,spinz,spin2,srspin
      double precision chiral,dnsty,tmpdns,tmpldns,tmpsdns,dtchrl,dtdkin,dkin
      double precision avgzeta,avgtorque,abschiral,point
      double precision avgspina,avgspinb,avgspinc,avgspinx,avgspiny,avgspinz
      double precision avgsrspin,avgspin2
      character(len=12) version
      character(len=14) MOJI
      character(len=24) today
      character(len=1), allocatable :: KRAM(:)
      complex(kind(0d0)) dstr0(9),detstr(9)
      complex(kind(0d0)) fs,fp,fd0,fd2,ff0,ff6,ff4,fg0,fg1,fg2,fg3
      complex(kind(0d0)) fh0,fh1,fh2,fh3,fh4
      complex(kind(0d0)) zetax0,zetay0,zetaz0,detztx,detzty,detztz
      complex(kind(0d0)) spina0,spinb0,spinc0,gg
      complex(kind(0d0)) spinx0,spiny0
      complex(kind(0d0)) uI,tpchrl,dtspna,dtspnb
      complex(kind(0d0)) dtspnx,dtspny
      integer NWF,i,j,k,k0,kk,mix,miy,miz,ix,iy,iz,itemp
      integer norb,neact,detcnt
      integer ndet1,ndet2,ne1,ne2,ne3,ne4,tmpck
      complex(kind=8) :: intN_mat

      version ='31th March 2012'
      call FDATE(today)
!      memmax = 10000 ! 4*520000


      write(*,*) '*****************************************************'
      write(*,*) 'MRDFT CI Ver. ',version
      write(*,*) '*****************************************************'
      write(*,*) ' '

      write(*,*) '# Calling read_spin.inp'
      call read_spininp
      write(*,*) '# Calling read_DiracOutput'
      call read_DiracOutput

      write(*,*) '       NBS0=',NBS0

      open (unit=41,file=FIFI3) !FIFI3 = 'DIRAC.out'
      do 
         read (41,'(A14)') MOJI
         if (MOJI.eq.' Core Energy :') exit
      end do
      read (41,1300) norb
 1300 FORMAT(47X,I4)
 
      !allocate(ras(20,100)) ! ,iplst(40,memmax)) ! norb >> detcnt is assumed
      allocate(enorb(norb),occdet(1000))
      allocate(KRAM(norb))

      call CIread(norb,NROOT,enorb,ras,KRAM,neact,detcnt,occdet)

!       do j=1, detcnt
!       write(*,*) 'Determinant :',j,(ras(i,j),i=1,neact)
!       write(*,*) 'Coefficient :',occdet(j)
!       end do 

      allocate(aa(NBS0,4),n(NBS0,4),a(NBS0,4,norb))
      allocate(xx(NBS0,4),yy(NBS0,4),zz(NBS0,4))
      allocate(f(NBS0,4),fx(NBS0,4),fy(NBS0,4),fz(NBS0,4))
      allocate(f00(4,norb),fx00(4,norb),fy00(4,norb),fz00(4,norb))
      allocate(fxx00(4,norb),fyy00(4,norb),fzz00(4,norb))
      allocate(g00(4,norb))
      allocate(fxx(NBS0,4),fyy(NBS0,4),fzz(NBS0,4))
!      allocate(fxy(NBS0,4),fxz(NBS0,4),fyz(NBS0,4))
!      allocate(fracn(NWF))
!      do j=1,norb
!         do i=1,NBS0
!            do k=1,4
!               a(i,k,j)=(0.d0,0.d0)
!            end do
!         end do
!      end do

      !convert QED coef to MRDFT coef
      call convert_pg_QEDtoMRDFT(aa,n,xx,yy,zz)
      call convert_coef_MRDFTCI(a,KRAM,norb,enorb)

      deallocate(cn,xc,yc,zc)  ! global
      deallocate(aa_L,xx_L,yy_L,zz_L,nx_L,ny_L,nz_L)
      deallocate(aa_S,xx_S,yy_S,zz_S,nx_S,ny_S,nz_S)
      deallocate(c_La,c_Lb,c_Sa,c_Sb)
      deallocate(d_La,d_Lb,d_Sa,d_Sb)
      deallocate(e_eig,p_eig)

      write(*,*) 'NBS0=',NBS0,'NBS=',NBS

      write(*,*) ' '
      write(*,*) 'Finish : read wave functions '
      write(*,*) ' '
      write(*,*) 'Start Calculations '
      write(*,*) ' '

!      do i=1,NBS0
!         do k=1,4
!            write(50,*) aa(i,k), xx(i,k), yy(i,k), zz(i,k)
!         end do
!      end do
!      do j=1,norb
!      do j=1,norb
!         do i=1,NBS0
!            do k=1,4
!               write(50,*) a(i,k,j)
!            end do
!         end do
!      end do
!      stop

      uI = cmplx(0.0D0,1.0D0)
      ccc = 137.035999809442D0


!      do j=1,NBS0
!        write(*,'(4e12.4,i4)') a(j,3,norb),a(j,4,norb),j
!      end do
      write(*,*) 'Calculation start date : ',today

      avgzeta = 0.0D0
      avgtorque = 0.0D0
      abschiral = 0.0D0

      avgspina = 0.0D0
      avgspinb = 0.0D0
      avgspinc = 0.0D0
      avgspinx = 0.0D0
      avgspiny = 0.0D0
      avgspinz = 0.0D0
      avgspin2 = 0.0D0
      avgsrspin = 0.0D0

      mix = int(u1/uu)
      miy = int(u2/uu)
      miz = int(u3/uu)
      do ix=-mix,mix
        w1 = uu*dble(ix) + u4
        do iy=-miy,miy
          w2 = uu*dble(iy) + u5
          do iz=-miz,miz
           w3 = uu*dble(iz) + u6

          x = w1
          y = w2
          z = w3

          dnsty  = 0.0D0
          dkin   = 0.0D0
          ldnsty = 0.0d0
          sdnsty = 0.0d0
          spina0 = 0.0D0
          spinb0 = 0.0D0
          spinx0 = 0.0D0
          spiny0 = 0.0D0


          do kk=1,9
            stress(kk) = 0.0D0
          end do
          chiral = 0.0D0
          etaux = 0.0D0
          etauy = 0.0D0
          etauz = 0.0D0
          zetax = 0.0D0
          zetay = 0.0D0
          zetaz = 0.0D0

        do k=1, norb
!          write(*,*) ' norb=',k,norb
          do k0=1,4
            do itemp=1,NBS0 !NBS->NBS0 edited by fukuda
            f(itemp,k0) = cmplx(0.0D0,0.0D0)
            fx(itemp,k0) = cmplx(0.0D0,0.0D0)
            fy(itemp,k0) = cmplx(0.0D0,0.0D0)
            fz(itemp,k0) = cmplx(0.0D0,0.0D0)
            end do
            f00(k0,k) = cmplx(0.0D0,0.0D0)
            fx00(k0,k) = cmplx(0.0D0,0.0D0)
            fy00(k0,k) = cmplx(0.0D0,0.0D0)
            fz00(k0,k) = cmplx(0.0D0,0.0D0)
            fxx00(k0,k) = cmplx(0.0D0,0.0D0)
            fyy00(k0,k) = cmplx(0.0D0,0.0D0)
            fzz00(k0,k) = cmplx(0.0D0,0.0D0)
            g00(k0,k) = cmplx(0.0D0,0.0D0)
          end do

          do k0=1, 4
          do j=1, NBS0 !NBS->NBS0 edited by fukuda

            xxx = x-xx(j,k0)
            yyy = y-yy(j,k0)
            zzz = z-zz(j,k0)
            rr  = xxx*xxx + yyy*yyy + zzz*zzz
            gg  = a(j,k0,k)*dexp(-aa(j,k0)*rr)
            if (n(j,k0) .eq. 0) then
             f(j,k0) = fs(aa(j,k0),a(j,k0,k),rr)
             fx(j,k0) = -2.0D0*aa(j,k0)*xxx*f(j,k0)
             fy(j,k0) = -2.0D0*aa(j,k0)*yyy*f(j,k0)
             fz(j,k0) = -2.0D0*aa(j,k0)*zzz*f(j,k0)
             fxx(j,k0) = f(j,k0)*aa(j,k0)*(4.0D0*aa(j,k0)*xxx*xxx-2.0D0)
             fyy(j,k0) = f(j,k0)*aa(j,k0)*(4.0D0*aa(j,k0)*yyy*yyy-2.0D0)
             fzz(j,k0) = f(j,k0)*aa(j,k0)*(4.0D0*aa(j,k0)*zzz*zzz-2.0D0)

!             write(*,*) 'j= ', j, 'k0= ', k0
!             write(*,*) 'xx= ',xx(j,k0),'yy=',yy(j,k0),'zz= ',zz(j,k0)

            else if (n(j,k0) .eq. 1) then
             f(j,k0) = fp(aa(j,k0),xxx,a(j,k0,k),rr)
             fx(j,k0) = -2.0D0*aa(j,k0)*xxx*f(j,k0) + gg
             fy(j,k0) = -2.0D0*aa(j,k0)*yyy*f(j,k0)
             fz(j,k0) = -2.0D0*aa(j,k0)*zzz*f(j,k0)
             fxx(j,k0) = f(j,k0)*aa(j,k0)*(4.0D0*aa(j,k0)*xxx*xxx-6.0D0)
             fyy(j,k0) = f(j,k0)*aa(j,k0)*(4.0D0*aa(j,k0)*yyy*yyy-2.0D0)
             fzz(j,k0) = f(j,k0)*aa(j,k0)*(4.0D0*aa(j,k0)*zzz*zzz-2.0D0)
            else if (n(j,k0) .eq. 2) then
             f(j,k0) = fp(aa(j,k0),yyy,a(j,k0,k),rr)
             fx(j,k0) = -2.0D0*aa(j,k0)*xxx*f(j,k0)
             fy(j,k0) = -2.0D0*aa(j,k0)*yyy*f(j,k0) + gg
             fz(j,k0) = -2.0D0*aa(j,k0)*zzz*f(j,k0)
             fxx(j,k0) = f(j,k0)*aa(j,k0)*(4.0D0*aa(j,k0)*xxx*xxx-2.0D0)
             fyy(j,k0) = f(j,k0)*aa(j,k0)*(4.0D0*aa(j,k0)*yyy*yyy-6.0D0)
             fzz(j,k0) = f(j,k0)*aa(j,k0)*(4.0D0*aa(j,k0)*zzz*zzz-2.0D0)
            else if (n(j,k0) .eq. 3) then
             f(j,k0) = fp(aa(j,k0),zzz,a(j,k0,k),rr)
             fx(j,k0) = -2.0D0*aa(j,k0)*xxx*f(j,k0)
             fy(j,k0) = -2.0D0*aa(j,k0)*yyy*f(j,k0)
             fz(j,k0) = -2.0D0*aa(j,k0)*zzz*f(j,k0) + gg
             fxx(j,k0) = f(j,k0)*aa(j,k0)*(4.0D0*aa(j,k0)*xxx*xxx-2.0D0)
             fyy(j,k0) = f(j,k0)*aa(j,k0)*(4.0D0*aa(j,k0)*yyy*yyy-2.0D0)
             fzz(j,k0) = f(j,k0)*aa(j,k0)*(4.0D0*aa(j,k0)*zzz*zzz-6.0D0)
            else if (n(j,k0) .eq. 4) then
             f(j,k0) = fd0(aa(j,k0),xxx,a(j,k0,k),rr)
             fx(j,k0) = -2.0D0*aa(j,k0)*xxx*f(j,k0) + 2.0D0*xxx*gg
             fy(j,k0) = -2.0D0*aa(j,k0)*yyy*f(j,k0)
             fz(j,k0) = -2.0D0*aa(j,k0)*zzz*f(j,k0)
             fxx(j,k0) = f(j,k0)*aa(j,k0)*(4.0D0*aa(j,k0)*xxx*xxx-10.D0) +2.0D0*gg
             fyy(j,k0) = f(j,k0)*aa(j,k0)*(4.0D0*aa(j,k0)*yyy*yyy-2.0D0)
             fzz(j,k0) = f(j,k0)*aa(j,k0)*(4.0D0*aa(j,k0)*zzz*zzz-2.0D0)
            else if (n(j,k0) .eq. 5) then
             f(j,k0) = fd2(aa(j,k0),xxx,yyy,a(j,k0,k),rr)
             fx(j,k0) = -2.0D0*aa(j,k0)*xxx*f(j,k0) + yyy*gg
             fy(j,k0) = -2.0D0*aa(j,k0)*yyy*f(j,k0) + xxx*gg
             fz(j,k0) = -2.0D0*aa(j,k0)*zzz*f(j,k0)
             fxx(j,k0) = f(j,k0)*aa(j,k0)*(4.0D0*aa(j,k0)*xxx*xxx-6.0D0)
             fyy(j,k0) = f(j,k0)*aa(j,k0)*(4.0D0*aa(j,k0)*yyy*yyy-6.0D0)
             fzz(j,k0) = f(j,k0)*aa(j,k0)*(4.0D0*aa(j,k0)*zzz*zzz-2.0D0)
            else if (n(j,k0) .eq. 6) then
             f(j,k0) = fd2(aa(j,k0),xxx,zzz,a(j,k0,k),rr)
             fx(j,k0) = -2.0D0*aa(j,k0)*xxx*f(j,k0) + zzz*gg
             fy(j,k0) = -2.0D0*aa(j,k0)*yyy*f(j,k0)
             fz(j,k0) = -2.0D0*aa(j,k0)*zzz*f(j,k0) + xxx*gg
             fxx(j,k0) = f(j,k0)*aa(j,k0)*(4.0D0*aa(j,k0)*xxx*xxx-6.0D0)
             fyy(j,k0) = f(j,k0)*aa(j,k0)*(4.0D0*aa(j,k0)*yyy*yyy-2.0D0)
             fzz(j,k0) = f(j,k0)*aa(j,k0)*(4.0D0*aa(j,k0)*zzz*zzz-6.0D0)
            else if (n(j,k0) .eq. 7) then
             f(j,k0) = fd0(aa(j,k0),yyy,a(j,k0,k),rr)
             fx(j,k0) = -2.0D0*aa(j,k0)*xxx*f(j,k0)
             fy(j,k0) = -2.0D0*aa(j,k0)*yyy*f(j,k0) + 2.0D0*yyy*gg
             fz(j,k0) = -2.0D0*aa(j,k0)*zzz*f(j,k0)
             fxx(j,k0) = f(j,k0)*aa(j,k0)*(4.0D0*aa(j,k0)*xxx*xxx-2.0D0)
             fyy(j,k0) = f(j,k0)*aa(j,k0)*(4.0D0*aa(j,k0)*yyy*yyy-10.D0) +2.0D0*gg
             fzz(j,k0) = f(j,k0)*aa(j,k0)*(4.0D0*aa(j,k0)*zzz*zzz-2.0D0)
            else if (n(j,k0) .eq. 8) then
             f(j,k0) = fd2(aa(j,k0),yyy,zzz,a(j,k0,k),rr)
             fx(j,k0) = -2.0D0*aa(j,k0)*xxx*f(j,k0)
             fy(j,k0) = -2.0D0*aa(j,k0)*yyy*f(j,k0) + zzz*gg
             fz(j,k0) = -2.0D0*aa(j,k0)*zzz*f(j,k0) + yyy*gg
             fxx(j,k0) = f(j,k0)*aa(j,k0)*(4.0D0*aa(j,k0)*xxx*xxx-2.0D0)
             fyy(j,k0) = f(j,k0)*aa(j,k0)*(4.0D0*aa(j,k0)*yyy*yyy-6.0D0)
             fzz(j,k0) = f(j,k0)*aa(j,k0)*(4.0D0*aa(j,k0)*zzz*zzz-6.0D0)
            else if (n(j,k0) .eq. 9) then
             f(j,k0) = fd0(aa(j,k0),zzz,a(j,k0,k),rr)
             fx(j,k0) = -2.0D0*aa(j,k0)*xxx*f(j,k0)
             fy(j,k0) = -2.0D0*aa(j,k0)*yyy*f(j,k0)
             fz(j,k0) = -2.0D0*aa(j,k0)*zzz*f(j,k0) + 2.0D0*zzz*gg !fukuda delete "+2.0D0*gg"
             fxx(j,k0) = f(j,k0)*aa(j,k0)*(4.0D0*aa(j,k0)*xxx*xxx-2.0D0)
             fyy(j,k0) = f(j,k0)*aa(j,k0)*(4.0D0*aa(j,k0)*yyy*yyy-2.0D0)
             fzz(j,k0) = f(j,k0)*aa(j,k0)*(4.0D0*aa(j,k0)*zzz*zzz-10.D0) +2.0D0*gg
            else if (n(j,k0) .eq. 10) then
             f(j,k0) = ff0(aa(j,k0),xxx,a(j,k0,k),rr)
             fx(j,k0) = -2.0D0*aa(j,k0)*xxx*f(j,k0) + 3.0D0*xxx*xxx*gg
             fy(j,k0) = -2.0D0*aa(j,k0)*yyy*f(j,k0)
             fz(j,k0) = -2.0D0*aa(j,k0)*zzz*f(j,k0)
             fxx(j,k0) = f(j,k0)*aa(j,k0)*(4.0D0*aa(j,k0)*xxx*xxx-14.D0) +6.0D0*xxx*gg
             fyy(j,k0) = f(j,k0)*aa(j,k0)*(4.0D0*aa(j,k0)*yyy*yyy-2.0D0)
             fzz(j,k0) = f(j,k0)*aa(j,k0)*(4.0D0*aa(j,k0)*zzz*zzz-2.0D0)
            else if (n(j,k0) .eq. 16) then
             f(j,k0) = ff0(aa(j,k0),yyy,a(j,k0,k),rr)
             fx(j,k0) = -2.0D0*aa(j,k0)*xxx*f(j,k0)
             fy(j,k0) = -2.0D0*aa(j,k0)*yyy*f(j,k0) + 3.0D0*yyy*yyy*gg
             fz(j,k0) = -2.0D0*aa(j,k0)*zzz*f(j,k0)
             fxx(j,k0) = f(j,k0)*aa(j,k0)*(4.0D0*aa(j,k0)*xxx*xxx-2.0D0)
             fyy(j,k0) = f(j,k0)*aa(j,k0)*(4.0D0*aa(j,k0)*yyy*yyy-14.D0) +6.0D0*yyy*gg
             fzz(j,k0) = f(j,k0)*aa(j,k0)*(4.0D0*aa(j,k0)*zzz*zzz-2.0D0)
            else if (n(j,k0) .eq. 19) then
             f(j,k0) = ff0(aa(j,k0),zzz,a(j,k0,k),rr)
             fx(j,k0) = -2.0D0*aa(j,k0)*xxx*f(j,k0)
             fy(j,k0) = -2.0D0*aa(j,k0)*yyy*f(j,k0)
             fz(j,k0) = -2.0D0*aa(j,k0)*zzz*f(j,k0) + 3.0D0*zzz*zzz*gg
             fxx(j,k0) = f(j,k0)*aa(j,k0)*(4.0D0*aa(j,k0)*xxx*xxx-2.0D0)
             fyy(j,k0) = f(j,k0)*aa(j,k0)*(4.0D0*aa(j,k0)*yyy*yyy-2.0D0)
             fzz(j,k0) = f(j,k0)*aa(j,k0)*(4.0D0*aa(j,k0)*zzz*zzz-14.D0) +6.0D0*zzz*gg
            else if (n(j,k0) .eq. 11) then
             f(j,k0) = ff6(aa(j,k0),xxx,yyy,a(j,k0,k),rr)
             fx(j,k0) = -2.0D0*aa(j,k0)*xxx*f(j,k0) + 2.0D0*xxx*yyy*gg
             fy(j,k0) = -2.0D0*aa(j,k0)*yyy*f(j,k0) + xxx*xxx*gg
             fz(j,k0) = -2.0D0*aa(j,k0)*zzz*f(j,k0)
             fxx(j,k0) = f(j,k0)*aa(j,k0)*(4.0D0*aa(j,k0)*xxx*xxx-10.D0) +2.0D0*yyy*gg
             fyy(j,k0) = f(j,k0)*aa(j,k0)*(4.0D0*aa(j,k0)*yyy*yyy-6.0D0)
             fzz(j,k0) = f(j,k0)*aa(j,k0)*(4.0D0*aa(j,k0)*zzz*zzz-2.0D0)
            else if (n(j,k0) .eq. 12) then
             f(j,k0) = ff6(aa(j,k0),xxx,zzz,a(j,k0,k),rr)
             fx(j,k0) = -2.0D0*aa(j,k0)*xxx*f(j,k0) + 2.0D0*xxx*zzz*gg
             fy(j,k0) = -2.0D0*aa(j,k0)*yyy*f(j,k0)
             fz(j,k0) = -2.0D0*aa(j,k0)*zzz*f(j,k0) + xxx*xxx*gg
             fxx(j,k0) = f(j,k0)*aa(j,k0)*(4.0D0*aa(j,k0)*xxx*xxx-10.D0) +2.0D0*zzz*gg
             fyy(j,k0) = f(j,k0)*aa(j,k0)*(4.0D0*aa(j,k0)*yyy*yyy-2.0D0)
             fzz(j,k0) = f(j,k0)*aa(j,k0)*(4.0D0*aa(j,k0)*zzz*zzz-6.0D0)
            else if (n(j,k0) .eq. 13) then
             f(j,k0) = ff6(aa(j,k0),yyy,xxx,a(j,k0,k),rr)
             fx(j,k0) = -2.0D0*aa(j,k0)*xxx*f(j,k0) + yyy*yyy*gg
             fy(j,k0) = -2.0D0*aa(j,k0)*yyy*f(j,k0) + 2.0D0*xxx*yyy*gg
             fz(j,k0) = -2.0D0*aa(j,k0)*zzz*f(j,k0)
             fxx(j,k0) = f(j,k0)*aa(j,k0)*(4.0D0*aa(j,k0)*xxx*xxx-6.0D0)
             fyy(j,k0) = f(j,k0)*aa(j,k0)*(4.0D0*aa(j,k0)*yyy*yyy-10.D0) +2.0D0*xxx*gg
             fzz(j,k0) = f(j,k0)*aa(j,k0)*(4.0D0*aa(j,k0)*zzz*zzz-2.0D0)
            else if (n(j,k0) .eq. 15) then
             f(j,k0) = ff6(aa(j,k0),zzz,xxx,a(j,k0,k),rr)
             fx(j,k0) = -2.0D0*aa(j,k0)*xxx*f(j,k0) + zzz*zzz*gg
             fy(j,k0) = -2.0D0*aa(j,k0)*yyy*f(j,k0)
             fz(j,k0) = -2.0D0*aa(j,k0)*zzz*f(j,k0) + 2.0D0*xxx*zzz*gg
             fxx(j,k0) = f(j,k0)*aa(j,k0)*(4.0D0*aa(j,k0)*xxx*xxx-6.0D0)
             fyy(j,k0) = f(j,k0)*aa(j,k0)*(4.0D0*aa(j,k0)*yyy*yyy-2.0D0)
             fzz(j,k0) = f(j,k0)*aa(j,k0)*(4.0D0*aa(j,k0)*zzz*zzz-10.D0) +2.0D0*xxx*gg
            else if (n(j,k0) .eq. 17) then
             f(j,k0) = ff6(aa(j,k0),yyy,zzz,a(j,k0,k),rr)
             fx(j,k0) = -2.0D0*aa(j,k0)*xxx*f(j,k0)
             fy(j,k0) = -2.0D0*aa(j,k0)*yyy*f(j,k0) + 2.0D0*yyy*zzz*gg
             fz(j,k0) = -2.0D0*aa(j,k0)*zzz*f(j,k0) + yyy*yyy*gg
             fxx(j,k0) = f(j,k0)*aa(j,k0)*(4.0D0*aa(j,k0)*xxx*xxx-2.0D0)
             fyy(j,k0) = f(j,k0)*aa(j,k0)*(4.0D0*aa(j,k0)*yyy*yyy-10.D0) +2.0D0*zzz*gg
             fzz(j,k0) = f(j,k0)*aa(j,k0)*(4.0D0*aa(j,k0)*zzz*zzz-6.0D0)
            else if (n(j,k0) .eq. 18) then
             f(j,k0) = ff6(aa(j,k0),zzz,yyy,a(j,k0,k),rr)
             fx(j,k0) = -2.0D0*aa(j,k0)*xxx*f(j,k0)
             fy(j,k0) = -2.0D0*aa(j,k0)*yyy*f(j,k0) + zzz*zzz*gg
             fz(j,k0) = -2.0D0*aa(j,k0)*zzz*f(j,k0) + 2.0D0*yyy*zzz*gg
             fxx(j,k0) = f(j,k0)*aa(j,k0)*(4.0D0*aa(j,k0)*xxx*xxx-2.0D0)
             fyy(j,k0) = f(j,k0)*aa(j,k0)*(4.0D0*aa(j,k0)*yyy*yyy-6.0D0)
             fzz(j,k0) = f(j,k0)*aa(j,k0)*(4.0D0*aa(j,k0)*zzz*zzz-10.D0) +2.0D0*yyy*gg
            else if (n(j,k0) .eq. 14) then
             f(j,k0) = ff4(aa(j,k0),xxx,yyy,zzz,a(j,k0,k),rr)
             fx(j,k0) = -2.0D0*aa(j,k0)*xxx*f(j,k0) + yyy*zzz*gg
             fy(j,k0) = -2.0D0*aa(j,k0)*yyy*f(j,k0) + xxx*zzz*gg
             fz(j,k0) = -2.0D0*aa(j,k0)*zzz*f(j,k0) + xxx*yyy*gg
             fxx(j,k0) = f(j,k0)*aa(j,k0)*(4.0D0*aa(j,k0)*xxx*xxx-6.0D0)
             fyy(j,k0) = f(j,k0)*aa(j,k0)*(4.0D0*aa(j,k0)*yyy*yyy-6.0D0)
             fzz(j,k0) = f(j,k0)*aa(j,k0)*(4.0D0*aa(j,k0)*zzz*zzz-6.0D0)
            else if (n(j,k0) .eq. 20) then
             f(j,k0) = fg0(aa(j,k0),xxx,a(j,k0,k),rr)
             fx(j,k0) = -2.0D0*aa(j,k0)*xxx*f(j,k0) + 4.0D0*xxx*xxx*xxx*gg
             fy(j,k0) = -2.0D0*aa(j,k0)*yyy*f(j,k0)
             fz(j,k0) = -2.0D0*aa(j,k0)*zzz*f(j,k0)
             fxx(j,k0) = f(j,k0)*aa(j,k0)*(4.0D0*aa(j,k0)*xxx*xxx-18.D0) +12.0D0*xxx*xxx*gg
             fyy(j,k0) = f(j,k0)*aa(j,k0)*(4.0D0*aa(j,k0)*yyy*yyy-2.0D0)
             fzz(j,k0) = f(j,k0)*aa(j,k0)*(4.0D0*aa(j,k0)*zzz*zzz-2.0D0)
            else if (n(j,k0) .eq. 30) then
             f(j,k0) = fg0(aa(j,k0),yyy,a(j,k0,k),rr)
             fx(j,k0) = -2.0D0*aa(j,k0)*xxx*f(j,k0)
             fy(j,k0) = -2.0D0*aa(j,k0)*yyy*f(j,k0) + 4.0D0*yyy*yyy*yyy*gg
             fz(j,k0) = -2.0D0*aa(j,k0)*zzz*f(j,k0)
             fxx(j,k0) = f(j,k0)*aa(j,k0)*(4.0D0*aa(j,k0)*xxx*xxx-2.0D0)
             fyy(j,k0) = f(j,k0)*aa(j,k0)*(4.0D0*aa(j,k0)*yyy*yyy-18.D0) +12.0D0*yyy*yyy*gg
             fzz(j,k0) = f(j,k0)*aa(j,k0)*(4.0D0*aa(j,k0)*zzz*zzz-2.0D0)
            else if (n(j,k0) .eq. 34) then
             f(j,k0) = fg0(aa(j,k0),zzz,a(j,k0,k),rr)
             fx(j,k0) = -2.0D0*aa(j,k0)*xxx*f(j,k0)
             fy(j,k0) = -2.0D0*aa(j,k0)*yyy*f(j,k0)
             fz(j,k0) = -2.0D0*aa(j,k0)*zzz*f(j,k0) + 4.0D0*zzz*zzz*zzz*gg
             fxx(j,k0) = f(j,k0)*aa(j,k0)*(4.0D0*aa(j,k0)*xxx*xxx-2.0D0)
             fyy(j,k0) = f(j,k0)*aa(j,k0)*(4.0D0*aa(j,k0)*yyy*yyy-2.0D0)
             fzz(j,k0) = f(j,k0)*aa(j,k0)*(4.0D0*aa(j,k0)*zzz*zzz-18.D0) +12.0D0*zzz*zzz*gg
            else if (n(j,k0) .eq. 21) then
             f(j,k0) = fg1(aa(j,k0),xxx,yyy,a(j,k0,k),rr)
             fx(j,k0) = -2.0D0*aa(j,k0)*xxx*f(j,k0) + 3.0D0*xxx*xxx*yyy*gg
             fy(j,k0) = -2.0D0*aa(j,k0)*yyy*f(j,k0) + xxx*xxx*xxx*gg
             fz(j,k0) = -2.0D0*aa(j,k0)*zzz*f(j,k0)
             fxx(j,k0) = f(j,k0)*aa(j,k0)*(4.0D0*aa(j,k0)*xxx*xxx-14.D0) +6.0D0*xxx*yyy*gg
             fyy(j,k0) = f(j,k0)*aa(j,k0)*(4.0D0*aa(j,k0)*yyy*yyy-6.0D0)
             fzz(j,k0) = f(j,k0)*aa(j,k0)*(4.0D0*aa(j,k0)*zzz*zzz-2.0D0)
            else if (n(j,k0) .eq. 22) then
             f(j,k0) = fg1(aa(j,k0),xxx,zzz,a(j,k0,k),rr)
             fx(j,k0) = -2.0D0*aa(j,k0)*xxx*f(j,k0) + 3.0D0*xxx*xxx*zzz*gg
             fy(j,k0) = -2.0D0*aa(j,k0)*yyy*f(j,k0)
             fz(j,k0) = -2.0D0*aa(j,k0)*zzz*f(j,k0) + xxx*xxx*xxx*gg
             fxx(j,k0) = f(j,k0)*aa(j,k0)*(4.0D0*aa(j,k0)*xxx*xxx-14.D0) +6.0D0*xxx*zzz*gg
             fyy(j,k0) = f(j,k0)*aa(j,k0)*(4.0D0*aa(j,k0)*yyy*yyy-2.0D0)
             fzz(j,k0) = f(j,k0)*aa(j,k0)*(4.0D0*aa(j,k0)*zzz*zzz-6.0D0)
            else if (n(j,k0) .eq. 26) then
             f(j,k0) = fg1(aa(j,k0),yyy,xxx,a(j,k0,k),rr)
             fx(j,k0) = -2.0D0*aa(j,k0)*xxx*f(j,k0) + yyy*yyy*yyy*gg
             fy(j,k0) = -2.0D0*aa(j,k0)*yyy*f(j,k0) + 3.0D0*xxx*yyy*yyy*gg
             fz(j,k0) = -2.0D0*aa(j,k0)*zzz*f(j,k0)
             fxx(j,k0) = f(j,k0)*aa(j,k0)*(4.0D0*aa(j,k0)*xxx*xxx-6.0D0)
             fyy(j,k0) = f(j,k0)*aa(j,k0)*(4.0D0*aa(j,k0)*yyy*yyy-14.D0) +6.0D0*xxx*yyy*gg
             fzz(j,k0) = f(j,k0)*aa(j,k0)*(4.0D0*aa(j,k0)*zzz*zzz-2.0D0)
            else if (n(j,k0) .eq. 31) then
             f(j,k0) = fg1(aa(j,k0),yyy,zzz,a(j,k0,k),rr)
             fx(j,k0) = -2.0D0*aa(j,k0)*xxx*f(j,k0)
             fy(j,k0) = -2.0D0*aa(j,k0)*yyy*f(j,k0) + 3.0D0*yyy*yyy*zzz*gg
             fz(j,k0) = -2.0D0*aa(j,k0)*zzz*f(j,k0) + yyy*yyy*yyy*gg
             fxx(j,k0) = f(j,k0)*aa(j,k0)*(4.0D0*aa(j,k0)*xxx*xxx-2.0D0)
             fyy(j,k0) = f(j,k0)*aa(j,k0)*(4.0D0*aa(j,k0)*yyy*yyy-14.D0) +6.0D0*yyy*zzz*gg
             fzz(j,k0) = f(j,k0)*aa(j,k0)*(4.0D0*aa(j,k0)*zzz*zzz-6.0D0)
            else if (n(j,k0) .eq. 29) then
             f(j,k0) = fg1(aa(j,k0),zzz,xxx,a(j,k0,k),rr)
             fx(j,k0) = -2.0D0*aa(j,k0)*xxx*f(j,k0) + zzz*zzz*zzz*gg
             fy(j,k0) = -2.0D0*aa(j,k0)*yyy*f(j,k0)
             fz(j,k0) = -2.0D0*aa(j,k0)*zzz*f(j,k0) + 3.0D0*xxx*zzz*zzz*gg
             fxx(j,k0) = f(j,k0)*aa(j,k0)*(4.0D0*aa(j,k0)*xxx*xxx-6.0D0)
             fyy(j,k0) = f(j,k0)*aa(j,k0)*(4.0D0*aa(j,k0)*yyy*yyy-2.0D0)
             fzz(j,k0) = f(j,k0)*aa(j,k0)*(4.0D0*aa(j,k0)*zzz*zzz-14.D0) +6.0D0*xxx*zzz*gg
            else if (n(j,k0) .eq. 33) then
             f(j,k0) = fg1(aa(j,k0),zzz,yyy,a(j,k0,k),rr)
             fx(j,k0) = -2.0D0*aa(j,k0)*xxx*f(j,k0)
             fy(j,k0) = -2.0D0*aa(j,k0)*yyy*f(j,k0) + zzz*zzz*zzz*gg
             fz(j,k0) = -2.0D0*aa(j,k0)*zzz*f(j,k0) + 3.0D0*yyy*zzz*zzz*gg
             fxx(j,k0) = f(j,k0)*aa(j,k0)*(4.0D0*aa(j,k0)*xxx*xxx-2.0D0)
             fyy(j,k0) = f(j,k0)*aa(j,k0)*(4.0D0*aa(j,k0)*yyy*yyy-6.0D0)
             fzz(j,k0) = f(j,k0)*aa(j,k0)*(4.0D0*aa(j,k0)*zzz*zzz-14.D0) +6.0D0*yyy*zzz*gg
            else if (n(j,k0) .eq. 23) then
             f(j,k0) = fg2(aa(j,k0),xxx,yyy,a(j,k0,k),rr)
             fx(j,k0) = -2.0D0*aa(j,k0)*xxx*f(j,k0) + 2.0D0*xxx*yyy*yyy*gg
             fy(j,k0) = -2.0D0*aa(j,k0)*yyy*f(j,k0) + 2.0D0*xxx*xxx*yyy*gg
             fz(j,k0) = -2.0D0*aa(j,k0)*zzz*f(j,k0)
             fxx(j,k0) = f(j,k0)*aa(j,k0)*(4.0D0*aa(j,k0)*xxx*xxx-10.D0) +2.0D0*yyy*yyy*gg
             fyy(j,k0) = f(j,k0)*aa(j,k0)*(4.0D0*aa(j,k0)*yyy*yyy-10.D0) +2.0D0*xxx*xxx*gg
             fzz(j,k0) = f(j,k0)*aa(j,k0)*(4.0D0*aa(j,k0)*zzz*zzz-2.0D0)
            else if (n(j,k0) .eq. 25) then
             f(j,k0) = fg2(aa(j,k0),xxx,zzz,a(j,k0,k),rr)
             fx(j,k0) = -2.0D0*aa(j,k0)*xxx*f(j,k0) + 2.0D0*xxx*zzz*zzz*gg
             fy(j,k0) = -2.0D0*aa(j,k0)*yyy*f(j,k0)
             fz(j,k0) = -2.0D0*aa(j,k0)*zzz*f(j,k0) + 2.0D0*xxx*xxx*zzz*gg
             fxx(j,k0) = f(j,k0)*aa(j,k0)*(4.0D0*aa(j,k0)*xxx*xxx-10.D0) +2.0D0*zzz*zzz*gg
             fyy(j,k0) = f(j,k0)*aa(j,k0)*(4.0D0*aa(j,k0)*yyy*yyy-2.0D0)
             fzz(j,k0) = f(j,k0)*aa(j,k0)*(4.0D0*aa(j,k0)*zzz*zzz-10.D0) +2.0D0*xxx*xxx*gg
            else if (n(j,k0) .eq. 32) then
             f(j,k0) = fg2(aa(j,k0),yyy,zzz,a(j,k0,k),rr)
             fx(j,k0) = -2.0D0*aa(j,k0)*xxx*f(j,k0)
             fy(j,k0) = -2.0D0*aa(j,k0)*yyy*f(j,k0) + 2.0D0*yyy*zzz*zzz*gg
             fz(j,k0) = -2.0D0*aa(j,k0)*zzz*f(j,k0) + 2.0D0*yyy*yyy*zzz*gg
             fxx(j,k0) = f(j,k0)*aa(j,k0)*(4.0D0*aa(j,k0)*xxx*xxx-2.0D0)
             fyy(j,k0) = f(j,k0)*aa(j,k0)*(4.0D0*aa(j,k0)*yyy*yyy-10.D0) +2.0D0*zzz*zzz*gg
             fzz(j,k0) = f(j,k0)*aa(j,k0)*(4.0D0*aa(j,k0)*zzz*zzz-10.D0) +2.0D0*yyy*yyy*gg
            else if (n(j,k0) .eq. 24) then !g211
             f(j,k0) = fg3(aa(j,k0),xxx,yyy,zzz,a(j,k0,k),rr)
             fx(j,k0) = -2.0D0*aa(j,k0)*xxx*f(j,k0) + 2.0D0*xxx*yyy*zzz*gg
!             fy(j,k0) = -2.0D0*aa(j,k0)*yyy*f(j,k0)
!     &                  + 2.0D0*xxx*xxx*zzz*gg
             fy(j,k0) = -2.0D0*aa(j,k0)*yyy*f(j,k0) + xxx*xxx*zzz*gg
!             fz(j,k0) = -2.0D0*aa(j,k0)*zzz*f(j,k0)
!     &                  + 2.0D0*xxx*xxx*yyy*gg
             fz(j,k0) = -2.0D0*aa(j,k0)*zzz*f(j,k0) + xxx*xxx*yyy*gg
             fxx(j,k0) = f(j,k0)*aa(j,k0)*(4.0D0*aa(j,k0)*xxx*xxx-10.D0) +2.0D0*yyy*zzz*gg
             fyy(j,k0) = f(j,k0)*aa(j,k0)*(4.0D0*aa(j,k0)*yyy*yyy-6.0D0)
             fzz(j,k0) = f(j,k0)*aa(j,k0)*(4.0D0*aa(j,k0)*zzz*zzz-6.0D0)
            else if (n(j,k0) .eq. 27) then !g121
             f(j,k0) = fg3(aa(j,k0),yyy,xxx,zzz,a(j,k0,k),rr)
!             fx(j,k0) = -2.0D0*aa(j,k0)*xxx*f(j,k0)
!     &                  + 2.0D0*yyy*yyy*zzz*gg
             fx(j,k0) = -2.0D0*aa(j,k0)*xxx*f(j,k0) + yyy*yyy*zzz*gg
             fy(j,k0) = -2.0D0*aa(j,k0)*yyy*f(j,k0) + 2.0D0*xxx*yyy*zzz*gg
!             fz(j,k0) = -2.0D0*aa(j,k0)*zzz*f(j,k0)
!     &                  + 2.0D0*xxx*yyy*yyy*gg
             fz(j,k0) = -2.0D0*aa(j,k0)*zzz*f(j,k0) + xxx*yyy*yyy*gg
             fxx(j,k0) = f(j,k0)*aa(j,k0)*(4.0D0*aa(j,k0)*xxx*xxx-6.0D0)
             fyy(j,k0) = f(j,k0)*aa(j,k0)*(4.0D0*aa(j,k0)*yyy*yyy-10.D0) +2.0D0*xxx*zzz*gg
             fzz(j,k0) = f(j,k0)*aa(j,k0)*(4.0D0*aa(j,k0)*zzz*zzz-6.0D0)
            else if (n(j,k0) .eq. 28) then !g112
             f(j,k0) = fg3(aa(j,k0),zzz,xxx,yyy,a(j,k0,k),rr)
!             fx(j,k0) = -2.0D0*aa(j,k0)*xxx*f(j,k0)
!     &                  + 2.0D0*yyy*zzz*zzz*gg
             fx(j,k0) = -2.0D0*aa(j,k0)*xxx*f(j,k0) + yyy*zzz*zzz*gg
!             fy(j,k0) = -2.0D0*aa(j,k0)*yyy*f(j,k0)
!     &                  + 2.0D0*xxx*zzz*zzz*gg
             fy(j,k0) = -2.0D0*aa(j,k0)*yyy*f(j,k0) + xxx*zzz*zzz*gg
             fz(j,k0) = -2.0D0*aa(j,k0)*zzz*f(j,k0) + 2.0D0*xxx*yyy*zzz*gg
             fxx(j,k0) = f(j,k0)*aa(j,k0)*(4.0D0*aa(j,k0)*xxx*xxx-6.0D0)
             fyy(j,k0) = f(j,k0)*aa(j,k0)*(4.0D0*aa(j,k0)*yyy*yyy-6.0D0)
             fzz(j,k0) = f(j,k0)*aa(j,k0)*(4.0D0*aa(j,k0)*zzz*zzz-10.D0) +2.0D0*xxx*yyy*gg
            else if (n(j,k0) .eq. 35) then
             f(j,k0) = fh0(aa(j,k0),xxx,a(j,k0,k),rr)
             fx(j,k0) = -2.0D0*aa(j,k0)*xxx*f(j,k0) + 5.0D0*xxx*xxx*xxx*xxx*gg
             fy(j,k0) = -2.0D0*aa(j,k0)*yyy*f(j,k0)
             fz(j,k0) = -2.0D0*aa(j,k0)*zzz*f(j,k0)
             fxx(j,k0) = f(j,k0)*aa(j,k0)*(4.0D0*aa(j,k0)*xxx*xxx-22.D0) +20.0D0*xxx*xxx*xxx*gg
             fyy(j,k0) = f(j,k0)*aa(j,k0)*(4.0D0*aa(j,k0)*yyy*yyy-2.0D0)
             fzz(j,k0) = f(j,k0)*aa(j,k0)*(4.0D0*aa(j,k0)*zzz*zzz-2.0D0)
            else if (n(j,k0) .eq. 50) then
             f(j,k0) = fh0(aa(j,k0),yyy,a(j,k0,k),rr)
             fx(j,k0) = -2.0D0*aa(j,k0)*xxx*f(j,k0)
             fy(j,k0) = -2.0D0*aa(j,k0)*yyy*f(j,k0) + 5.0D0*yyy*yyy*yyy*yyy*gg
             fz(j,k0) = -2.0D0*aa(j,k0)*zzz*f(j,k0)
             fxx(j,k0) = f(j,k0)*aa(j,k0)*(4.0D0*aa(j,k0)*xxx*xxx-2.0D0)
             fyy(j,k0) = f(j,k0)*aa(j,k0)*(4.0D0*aa(j,k0)*yyy*yyy-22.D0) +20.0D0*yyy*yyy*yyy*gg
             fzz(j,k0) = f(j,k0)*aa(j,k0)*(4.0D0*aa(j,k0)*zzz*zzz-2.0D0)
            else if (n(j,k0) .eq. 55) then
             f(j,k0) = fh0(aa(j,k0),zzz,a(j,k0,k),rr)
             fx(j,k0) = -2.0D0*aa(j,k0)*xxx*f(j,k0)
             fy(j,k0) = -2.0D0*aa(j,k0)*yyy*f(j,k0)
             fz(j,k0) = -2.0D0*aa(j,k0)*zzz*f(j,k0) + 5.0D0*zzz*zzz*zzz*zzz*gg
             fxx(j,k0) = f(j,k0)*aa(j,k0)*(4.0D0*aa(j,k0)*xxx*xxx-2.0D0)
             fyy(j,k0) = f(j,k0)*aa(j,k0)*(4.0D0*aa(j,k0)*yyy*yyy-2.D0)
             fzz(j,k0) = f(j,k0)*aa(j,k0)*(4.0D0*aa(j,k0)*zzz*zzz-22.D0) +20.0D0*zzz*zzz*zzz*gg
            else if (n(j,k0) .eq. 36) then
             f(j,k0) = fh1(aa(j,k0),xxx,yyy,a(j,k0,k),rr)
             fx(j,k0) = -2.0D0*aa(j,k0)*xxx*f(j,k0) + 4.0D0*xxx*xxx*xxx*yyy*gg
             fy(j,k0) = -2.0D0*aa(j,k0)*yyy*f(j,k0) + xxx*xxx*xxx*xxx*gg
             fz(j,k0) = -2.0D0*aa(j,k0)*zzz*f(j,k0)
             fxx(j,k0) = f(j,k0)*aa(j,k0)*(4.0D0*aa(j,k0)*xxx*xxx-18.D0) +12.0D0*xxx*xxx*yyy*gg
             fyy(j,k0) = f(j,k0)*aa(j,k0)*(4.0D0*aa(j,k0)*yyy*yyy-6.0D0)
             fzz(j,k0) = f(j,k0)*aa(j,k0)*(4.0D0*aa(j,k0)*zzz*zzz-2.0D0)
            else if (n(j,k0) .eq. 37) then
             f(j,k0) = fh1(aa(j,k0),xxx,zzz,a(j,k0,k),rr)
             fx(j,k0) = -2.0D0*aa(j,k0)*xxx*f(j,k0) + 4.0D0*xxx*xxx*xxx*zzz*gg
             fy(j,k0) = -2.0D0*aa(j,k0)*yyy*f(j,k0)
             fz(j,k0) = -2.0D0*aa(j,k0)*zzz*f(j,k0) + xxx*xxx*xxx*xxx*gg
             fxx(j,k0) = f(j,k0)*aa(j,k0)*(4.0D0*aa(j,k0)*xxx*xxx-18.D0) +12.0D0*xxx*xxx*zzz*gg
             fyy(j,k0) = f(j,k0)*aa(j,k0)*(4.0D0*aa(j,k0)*yyy*yyy-2.0D0)
             fzz(j,k0) = f(j,k0)*aa(j,k0)*(4.0D0*aa(j,k0)*zzz*zzz-6.0D0)
            else if (n(j,k0) .eq. 45) then
             f(j,k0) = fh1(aa(j,k0),yyy,xxx,a(j,k0,k),rr)
             fx(j,k0) = -2.0D0*aa(j,k0)*xxx*f(j,k0) + yyy*yyy*yyy*yyy*gg
             fy(j,k0) = -2.0D0*aa(j,k0)*yyy*f(j,k0) + 4.0D0*xxx*yyy*yyy*yyy*gg
             fz(j,k0) = -2.0D0*aa(j,k0)*zzz*f(j,k0)
             fxx(j,k0) = f(j,k0)*aa(j,k0)*(4.0D0*aa(j,k0)*xxx*xxx-6.D0)
             fyy(j,k0) = f(j,k0)*aa(j,k0)*(4.0D0*aa(j,k0)*yyy*yyy-18.D0) +12.0D0*xxx*yyy*yyy*gg
             fzz(j,k0) = f(j,k0)*aa(j,k0)*(4.0D0*aa(j,k0)*zzz*zzz-2.0D0)
            else if (n(j,k0) .eq. 49) then
             f(j,k0) = fh1(aa(j,k0),zzz,xxx,a(j,k0,k),rr)
             fx(j,k0) = -2.0D0*aa(j,k0)*xxx*f(j,k0) + zzz*zzz*zzz*zzz*gg
             fy(j,k0) = -2.0D0*aa(j,k0)*yyy*f(j,k0)
             fz(j,k0) = -2.0D0*aa(j,k0)*zzz*f(j,k0) + 4.0D0*xxx*zzz*zzz*zzz*gg
             fxx(j,k0) = f(j,k0)*aa(j,k0)*(4.0D0*aa(j,k0)*xxx*xxx-6.0D0)
             fyy(j,k0) = f(j,k0)*aa(j,k0)*(4.0D0*aa(j,k0)*yyy*yyy-2.0D0)
             fzz(j,k0) = f(j,k0)*aa(j,k0)*(4.0D0*aa(j,k0)*zzz*zzz-18.D0) +12.0D0*xxx*zzz*zzz*gg
            else if (n(j,k0) .eq. 51) then
             f(j,k0) = fh1(aa(j,k0),yyy,zzz,a(j,k0,k),rr)
             fx(j,k0) = -2.0D0*aa(j,k0)*xxx*f(j,k0)
             fy(j,k0) = -2.0D0*aa(j,k0)*yyy*f(j,k0) + 4.0D0*yyy*yyy*yyy*zzz*gg
             fz(j,k0) = -2.0D0*aa(j,k0)*zzz*f(j,k0) + yyy*yyy*yyy*yyy*gg
             fxx(j,k0) = f(j,k0)*aa(j,k0)*(4.0D0*aa(j,k0)*xxx*xxx-2.0D0)
             fyy(j,k0) = f(j,k0)*aa(j,k0)*(4.0D0*aa(j,k0)*yyy*yyy-18.D0) +12.0D0*yyy*yyy*zzz*gg
             fzz(j,k0) = f(j,k0)*aa(j,k0)*(4.0D0*aa(j,k0)*zzz*zzz-6.0D0)
            else if (n(j,k0) .eq. 54) then
             f(j,k0) = fh1(aa(j,k0),zzz,yyy,a(j,k0,k),rr)
             fx(j,k0) = -2.0D0*aa(j,k0)*xxx*f(j,k0)
             fy(j,k0) = -2.0D0*aa(j,k0)*yyy*f(j,k0) + zzz*zzz*zzz*zzz*gg
             fz(j,k0) = -2.0D0*aa(j,k0)*zzz*f(j,k0) + 4.0D0*yyy*zzz*zzz*zzz*gg
             fxx(j,k0) = f(j,k0)*aa(j,k0)*(4.0D0*aa(j,k0)*xxx*xxx-2.0D0)
             fyy(j,k0) = f(j,k0)*aa(j,k0)*(4.0D0*aa(j,k0)*yyy*yyy-6.0D0)
             fzz(j,k0) = f(j,k0)*aa(j,k0)*(4.0D0*aa(j,k0)*zzz*zzz-18.D0) +12.0D0*yyy*zzz*zzz*gg
            else if (n(j,k0) .eq. 38) then
             f(j,k0) = fh2(aa(j,k0),xxx,yyy,a(j,k0,k),rr)
             fx(j,k0) = -2.0D0*aa(j,k0)*xxx*f(j,k0) + 3.0D0*xxx*xxx*yyy*yyy*gg
             fy(j,k0) = -2.0D0*aa(j,k0)*yyy*f(j,k0) + 2.0D0*xxx*xxx*xxx*yyy*gg
             fz(j,k0) = -2.0D0*aa(j,k0)*zzz*f(j,k0)
             fxx(j,k0) = f(j,k0)*aa(j,k0)*(4.0D0*aa(j,k0)*xxx*xxx-14.D0) +6.0D0*xxx*yyy*yyy*gg
             fyy(j,k0) = f(j,k0)*aa(j,k0)*(4.0D0*aa(j,k0)*yyy*yyy-10.D0) +2.0D0*xxx*xxx*xxx*gg
             fzz(j,k0) = f(j,k0)*aa(j,k0)*(4.0D0*aa(j,k0)*zzz*zzz-2.0D0)
            else if (n(j,k0) .eq. 40) then
             f(j,k0) = fh2(aa(j,k0),xxx,zzz,a(j,k0,k),rr)
             fx(j,k0) = -2.0D0*aa(j,k0)*xxx*f(j,k0) + 3.0D0*xxx*xxx*zzz*zzz*gg
             fy(j,k0) = -2.0D0*aa(j,k0)*yyy*f(j,k0)
             fz(j,k0) = -2.0D0*aa(j,k0)*zzz*f(j,k0) + 2.0D0*xxx*xxx*xxx*zzz*gg
             fxx(j,k0) = f(j,k0)*aa(j,k0)*(4.0D0*aa(j,k0)*xxx*xxx-14.D0) +6.0D0*xxx*zzz*zzz*gg
             fyy(j,k0) = f(j,k0)*aa(j,k0)*(4.0D0*aa(j,k0)*yyy*yyy-2.0D0)
             fzz(j,k0) = f(j,k0)*aa(j,k0)*(4.0D0*aa(j,k0)*zzz*zzz-10.D0) +2.0D0*xxx*xxx*xxx*gg
            else if (n(j,k0) .eq. 41) then
             f(j,k0) = fh2(aa(j,k0),yyy,xxx,a(j,k0,k),rr)
             fx(j,k0) = -2.0D0*aa(j,k0)*xxx*f(j,k0) + 2.0D0*xxx*yyy*yyy*yyy*gg
             fy(j,k0) = -2.0D0*aa(j,k0)*yyy*f(j,k0) + 3.0D0*xxx*xxx*yyy*yyy*gg
             fz(j,k0) = -2.0D0*aa(j,k0)*zzz*f(j,k0)
             fxx(j,k0) = f(j,k0)*aa(j,k0)*(4.0D0*aa(j,k0)*xxx*xxx-10.D0) +2.0D0*yyy*yyy*yyy*gg
             fyy(j,k0) = f(j,k0)*aa(j,k0)*(4.0D0*aa(j,k0)*yyy*yyy-14.D0) +6.0D0*xxx*xxx*yyy*gg
             fzz(j,k0) = f(j,k0)*aa(j,k0)*(4.0D0*aa(j,k0)*zzz*zzz-2.0D0)
            else if (n(j,k0) .eq. 44) then
             f(j,k0) = fh2(aa(j,k0),zzz,xxx,a(j,k0,k),rr)
             fx(j,k0) = -2.0D0*aa(j,k0)*xxx*f(j,k0) + 2.0D0*xxx*zzz*zzz*zzz*gg
             fy(j,k0) = -2.0D0*aa(j,k0)*yyy*f(j,k0)
             fz(j,k0) = -2.0D0*aa(j,k0)*zzz*f(j,k0) + 3.0D0*xxx*xxx*zzz*zzz*gg
             fxx(j,k0) = f(j,k0)*aa(j,k0)*(4.0D0*aa(j,k0)*xxx*xxx-10.D0) +2.0D0*zzz*zzz*zzz*gg
             fyy(j,k0) = f(j,k0)*aa(j,k0)*(4.0D0*aa(j,k0)*yyy*yyy-2.0D0)
             fzz(j,k0) = f(j,k0)*aa(j,k0)*(4.0D0*aa(j,k0)*zzz*zzz-14.D0) +6.0D0*xxx*xxx*zzz*gg
            else if (n(j,k0) .eq. 52) then
             f(j,k0) = fh2(aa(j,k0),yyy,zzz,a(j,k0,k),rr)
             fx(j,k0) = -2.0D0*aa(j,k0)*xxx*f(j,k0)
             fy(j,k0) = -2.0D0*aa(j,k0)*yyy*f(j,k0) + 3.0D0*yyy*yyy*zzz*zzz*gg
             fz(j,k0) = -2.0D0*aa(j,k0)*zzz*f(j,k0) + 2.0D0*yyy*yyy*yyy*zzz*gg
             fxx(j,k0) = f(j,k0)*aa(j,k0)*(4.0D0*aa(j,k0)*xxx*xxx-2.0D0)
             fyy(j,k0) = f(j,k0)*aa(j,k0)*(4.0D0*aa(j,k0)*yyy*yyy-14.D0) +6.0D0*yyy*zzz*zzz*gg
             fzz(j,k0) = f(j,k0)*aa(j,k0)*(4.0D0*aa(j,k0)*zzz*zzz-10.D0) +2.0D0*yyy*yyy*yyy*gg
            else if (n(j,k0) .eq. 53) then
             f(j,k0) = fh2(aa(j,k0),zzz,yyy,a(j,k0,k),rr)
             fx(j,k0) = -2.0D0*aa(j,k0)*xxx*f(j,k0)
             fy(j,k0) = -2.0D0*aa(j,k0)*yyy*f(j,k0) + 2.0D0*yyy*zzz*zzz*zzz*gg
             fz(j,k0) = -2.0D0*aa(j,k0)*zzz*f(j,k0) + 3.0D0*yyy*yyy*zzz*zzz*gg
             fxx(j,k0) = f(j,k0)*aa(j,k0)*(4.0D0*aa(j,k0)*xxx*xxx-2.0D0)
             fyy(j,k0) = f(j,k0)*aa(j,k0)*(4.0D0*aa(j,k0)*yyy*yyy-10.D0) +2.0D0*zzz*zzz*zzz*gg
             fzz(j,k0) = f(j,k0)*aa(j,k0)*(4.0D0*aa(j,k0)*zzz*zzz-14.D0) +6.0D0*yyy*yyy*zzz*gg
            else if (n(j,k0) .eq. 39) then
             f(j,k0) = fh3(aa(j,k0),xxx,yyy,zzz,a(j,k0,k),rr)
             fx(j,k0) = -2.0D0*aa(j,k0)*xxx*f(j,k0) + 3.0D0*xxx*xxx*yyy*zzz*gg
             fy(j,k0) = -2.0D0*aa(j,k0)*yyy*f(j,k0) + xxx*xxx*xxx*zzz*gg
             fz(j,k0) = -2.0D0*aa(j,k0)*zzz*f(j,k0) + xxx*xxx*xxx*yyy*gg
             fxx(j,k0) = f(j,k0)*aa(j,k0)*(4.0D0*aa(j,k0)*xxx*xxx-14.D0) +6.0D0*xxx*yyy*zzz*gg
             fyy(j,k0) = f(j,k0)*aa(j,k0)*(4.0D0*aa(j,k0)*yyy*yyy-6.0D0)
             fzz(j,k0) = f(j,k0)*aa(j,k0)*(4.0D0*aa(j,k0)*zzz*zzz-6.0D0)
            else if (n(j,k0) .eq. 46) then
             f(j,k0) = fh3(aa(j,k0),yyy,xxx,zzz,a(j,k0,k),rr)
             fx(j,k0) = -2.0D0*aa(j,k0)*xxx*f(j,k0) + yyy*yyy*yyy*zzz*gg
             fy(j,k0) = -2.0D0*aa(j,k0)*yyy*f(j,k0) + 3.0D0*xxx*yyy*yyy*zzz*gg
             fz(j,k0) = -2.0D0*aa(j,k0)*zzz*f(j,k0) + xxx*yyy*yyy*yyy*gg
             fxx(j,k0) = f(j,k0)*aa(j,k0)*(4.0D0*aa(j,k0)*xxx*xxx-6.0D0)
             fyy(j,k0) = f(j,k0)*aa(j,k0)*(4.0D0*aa(j,k0)*yyy*yyy-14.D0) +6.0D0*xxx*yyy*zzz*gg
             fzz(j,k0) = f(j,k0)*aa(j,k0)*(4.0D0*aa(j,k0)*zzz*zzz-6.0D0)
            else if (n(j,k0) .eq. 48) then
             f(j,k0) = fh3(aa(j,k0),zzz,xxx,yyy,a(j,k0,k),rr)
             fx(j,k0) = -2.0D0*aa(j,k0)*xxx*f(j,k0) + yyy*zzz*zzz*zzz*gg
             fy(j,k0) = -2.0D0*aa(j,k0)*yyy*f(j,k0) + xxx*zzz*zzz*zzz*gg
             fz(j,k0) = -2.0D0*aa(j,k0)*zzz*f(j,k0) + 3.0D0*xxx*yyy*zzz*zzz*gg
             fxx(j,k0) = f(j,k0)*aa(j,k0)*(4.0D0*aa(j,k0)*xxx*xxx-6.0D0)
             fyy(j,k0) = f(j,k0)*aa(j,k0)*(4.0D0*aa(j,k0)*yyy*yyy-6.0D0)
             fzz(j,k0) = f(j,k0)*aa(j,k0)*(4.0D0*aa(j,k0)*zzz*zzz-14.D0) +6.0D0*xxx*yyy*zzz*gg
            else if (n(j,k0) .eq. 42) then
             f(j,k0) = fh4(aa(j,k0),xxx,yyy,zzz,a(j,k0,k),rr)
             fx(j,k0) = -2.0D0*aa(j,k0)*xxx*f(j,k0) + 2.0D0*xxx*yyy*yyy*zzz*gg
             fy(j,k0) = -2.0D0*aa(j,k0)*yyy*f(j,k0) + 2.0D0*xxx*xxx*yyy*zzz*gg
             fz(j,k0) = -2.0D0*aa(j,k0)*zzz*f(j,k0) + xxx*xxx*yyy*yyy*gg
             fxx(j,k0) = f(j,k0)*aa(j,k0)*(4.0D0*aa(j,k0)*xxx*xxx-10.D0) +2.0D0*yyy*yyy*zzz*gg
             fyy(j,k0) = f(j,k0)*aa(j,k0)*(4.0D0*aa(j,k0)*yyy*yyy-10.D0) +2.0D0*xxx*xxx*zzz*gg
             fzz(j,k0) = f(j,k0)*aa(j,k0)*(4.0D0*aa(j,k0)*zzz*zzz-6.0D0)
            else if (n(j,k0) .eq. 43) then
             f(j,k0) = fh4(aa(j,k0),xxx,zzz,yyy,a(j,k0,k),rr)
             fx(j,k0) = -2.0D0*aa(j,k0)*xxx*f(j,k0) + 2.0D0*xxx*yyy*zzz*zzz*gg
             fy(j,k0) = -2.0D0*aa(j,k0)*yyy*f(j,k0) + xxx*xxx*zzz*zzz*gg
             fz(j,k0) = -2.0D0*aa(j,k0)*zzz*f(j,k0) + 2.0D0*xxx*xxx*yyy*zzz*gg
             fxx(j,k0) = f(j,k0)*aa(j,k0)*(4.0D0*aa(j,k0)*xxx*xxx-10.D0) +2.0D0*yyy*zzz*zzz*gg
             fyy(j,k0) = f(j,k0)*aa(j,k0)*(4.0D0*aa(j,k0)*yyy*yyy-6.0D0)
             fzz(j,k0) = f(j,k0)*aa(j,k0)*(4.0D0*aa(j,k0)*zzz*zzz-10.D0) +2.0D0*xxx*xxx*yyy*gg
            else if (n(j,k0) .eq. 47) then
             f(j,k0) = fh4(aa(j,k0),yyy,zzz,xxx,a(j,k0,k),rr)
             fx(j,k0) = -2.0D0*aa(j,k0)*xxx*f(j,k0) + yyy*yyy*zzz*zzz*gg
             fy(j,k0) = -2.0D0*aa(j,k0)*yyy*f(j,k0) + 2.0D0*xxx*yyy*zzz*zzz*gg
             fz(j,k0) = -2.0D0*aa(j,k0)*zzz*f(j,k0) + 2.0D0*xxx*yyy*yyy*zzz*gg
             fxx(j,k0) = f(j,k0)*aa(j,k0)*(4.0D0*aa(j,k0)*xxx*xxx-6.0D0)
             fyy(j,k0) = f(j,k0)*aa(j,k0)*(4.0D0*aa(j,k0)*yyy*yyy-10.D0) +2.0D0*xxx*zzz*zzz*gg
             fzz(j,k0) = f(j,k0)*aa(j,k0)*(4.0D0*aa(j,k0)*zzz*zzz-10.D0) +2.0D0*xxx*yyy*yyy*gg
            end if

          end do ! j
          end do ! k0

          do k0=1,4
          do j=1,NBS0 !NBS->NBS0 edited by fukuda
             f00(k0,k) = f00(k0,k) + f(j,k0)
            fx00(k0,k) = fx00(k0,k) + fx(j,k0)
            fy00(k0,k) = fy00(k0,k) + fy(j,k0)
            fz00(k0,k) = fz00(k0,k) + fz(j,k0)
            fxx00(k0,k) = fxx00(k0,k) + fxx(j,k0)
            fyy00(k0,k) = fyy00(k0,k) + fyy(j,k0)
            fzz00(k0,k) = fzz00(k0,k) + fzz(j,k0)
          end do
!           write(*,*) f00(k0,k),k,k0
            g00(k0,k) = fxx00(k0,k)+fyy00(k0,k)+fzz00(k0,k)
          end do

!        write(*,*) ' Complete : derivatives of all orbital'

        end do ! k : norb

!
!  each determinant operation
!
        do ndet1 = 1, detcnt
        do ndet2 = 1, detcnt
           detdns = 0.0d0
           dtdkin = 0.0d0
           detldns = 0.0d0
           detsdns = 0.0d0
           dtspna = cmplx(0.0D0,0.0D0)
           dtspnb = cmplx(0.0D0,0.0D0)
           dtspnx = cmplx(0.0D0,0.0D0)
           dtspny = cmplx(0.0D0,0.0D0)
           dtchrl = 0.0d0
           detztx = cmplx(0.0D0,0.0D0)
           detzty = cmplx(0.0D0,0.0D0)
           detztz = cmplx(0.0D0,0.0D0)
          do kk=1,9
            detstr(kk) = cmplx(0.0D0,0.0D0)
          end do
!
!  Inner product patterns
!
!        write(*,*) ndet1,ndet2
!           call patdet(norb,neact,ndet1,ndet2,ras,
!     &                 iplst,nwnl,numele,memmax)

!         write(*,*) numele,nwnl,neact

!        do npat=1,nwnl
         do ne1 = 1, neact
         do ne2 = 1, neact
           tmpck = 0
!        check the other states
           do ne3=1,neact
           do ne4=1,neact
             if (ne3.eq.ne1) goto 200
             if (ne4.eq.ne2) goto 200
             if (ras(ne3,ndet1).eq.ras(ne4,ndet2)) tmpck=tmpck+1
  200        continue
           end do !ne3
           end do !ne4
           if (tmpck.ne.neact-1) goto 300
          
          tmpdns = 0.0d0
          tmpldns = 0.0d0
          tmpsdns = 0.0d0
          tpchrl = 0.0d0
          zetax0 = cmplx(0.0D0,0.0D0)
          zetay0 = cmplx(0.0D0,0.0D0)
          zetaz0 = cmplx(0.0D0,0.0D0)
          do kk=1,9
            dstr0(kk) = cmplx(0.0D0,0.0D0)
          end do ! kk

          tmpdns = tmpdns + dble(                               & 
     &        conjg(f00(1,ras(ne1,ndet1)))*f00(1,ras(ne2,ndet2))&
     &      + conjg(f00(2,ras(ne1,ndet1)))*f00(2,ras(ne2,ndet2))&
     &      + conjg(f00(3,ras(ne1,ndet1)))*f00(3,ras(ne2,ndet2))&
     &      + conjg(f00(4,ras(ne1,ndet1)))*f00(4,ras(ne2,ndet2)))

          dtdkin = dtdkin + dble(                               & 
     &        conjg(f00(1,ras(ne1,ndet1)))*g00(1,ras(ne2,ndet2))&
     &      + conjg(f00(2,ras(ne1,ndet1)))*g00(2,ras(ne2,ndet2))&
     &      + conjg(f00(3,ras(ne1,ndet1)))*g00(3,ras(ne2,ndet2))&
     &      + conjg(f00(4,ras(ne1,ndet1)))*g00(4,ras(ne2,ndet2)))

          tmpldns = tmpldns + dble(                             & 
     &        conjg(f00(1,ras(ne1,ndet1)))*f00(1,ras(ne2,ndet2))& 
     &      + conjg(f00(2,ras(ne1,ndet1)))*f00(2,ras(ne2,ndet2)))

          tmpsdns = tmpsdns + dble(                             & 
     &        conjg(f00(3,ras(ne1,ndet1)))*f00(3,ras(ne2,ndet2))&
     &      + conjg(f00(4,ras(ne1,ndet1)))*f00(4,ras(ne2,ndet2)))

          tpchrl = -1.0d0 * (                                   &
     &        conjg(f00(1,ras(ne1,ndet1)))*f00(3,ras(ne2,ndet2))&
     &      + conjg(f00(2,ras(ne1,ndet1)))*f00(4,ras(ne2,ndet2))&
     &      + conjg(f00(3,ras(ne1,ndet1)))*f00(1,ras(ne2,ndet2))&
     &      + conjg(f00(4,ras(ne1,ndet1)))*f00(2,ras(ne2,ndet2)))

          zetax0 = zetax0 +                                      &
     &       (conjg(fx00(1,ras(ne1,ndet1)))*f00(3,ras(ne2,ndet2))&
     &      + conjg(fx00(2,ras(ne1,ndet1)))*f00(4,ras(ne2,ndet2))&
     &      + conjg(fx00(3,ras(ne1,ndet1)))*f00(1,ras(ne2,ndet2))&
     &      + conjg(fx00(4,ras(ne1,ndet1)))*f00(2,ras(ne2,ndet2))&
     &      + conjg(f00(1,ras(ne1,ndet1)))*fx00(3,ras(ne2,ndet2))&
     &      + conjg(f00(2,ras(ne1,ndet1)))*fx00(4,ras(ne2,ndet2))&
     &      + conjg(f00(3,ras(ne1,ndet1)))*fx00(1,ras(ne2,ndet2))&
     &      + conjg(f00(4,ras(ne1,ndet1)))*fx00(2,ras(ne2,ndet2)))

          zetay0 = zetay0 +                                      &
     &       (conjg(fy00(1,ras(ne1,ndet1)))*f00(3,ras(ne2,ndet2))&
     &      + conjg(fy00(2,ras(ne1,ndet1)))*f00(4,ras(ne2,ndet2))&
     &      + conjg(fy00(3,ras(ne1,ndet1)))*f00(1,ras(ne2,ndet2))&
     &      + conjg(fy00(4,ras(ne1,ndet1)))*f00(2,ras(ne2,ndet2))&
     &      + conjg(f00(1,ras(ne1,ndet1)))*fy00(3,ras(ne2,ndet2))&
     &      + conjg(f00(2,ras(ne1,ndet1)))*fy00(4,ras(ne2,ndet2))&
     &      + conjg(f00(3,ras(ne1,ndet1)))*fy00(1,ras(ne2,ndet2))&
     &      + conjg(f00(4,ras(ne1,ndet1)))*fy00(2,ras(ne2,ndet2)))

          zetaz0 = zetaz0 +                                      &
     &       (conjg(fz00(1,ras(ne1,ndet1)))*f00(3,ras(ne2,ndet2))&
     &      + conjg(fz00(2,ras(ne1,ndet1)))*f00(4,ras(ne2,ndet2))&
     &      + conjg(fz00(3,ras(ne1,ndet1)))*f00(1,ras(ne2,ndet2))&
     &      + conjg(fz00(4,ras(ne1,ndet1)))*f00(2,ras(ne2,ndet2))&
     &      + conjg(f00(1,ras(ne1,ndet1)))*fz00(3,ras(ne2,ndet2))&
     &      + conjg(f00(2,ras(ne1,ndet1)))*fz00(4,ras(ne2,ndet2))&
     &      + conjg(f00(3,ras(ne1,ndet1)))*fz00(1,ras(ne2,ndet2))&
     &      + conjg(f00(4,ras(ne1,ndet1)))*fz00(2,ras(ne2,ndet2)))

          dstr0(1) = dstr0(1) + uI * (                            &
     &         conjg(f00(1,ras(ne1,ndet1)))*fx00(4,ras(ne2,ndet2))&
     &       + conjg(f00(2,ras(ne1,ndet1)))*fx00(3,ras(ne2,ndet2))&
     &       + conjg(f00(3,ras(ne1,ndet1)))*fx00(2,ras(ne2,ndet2))&
     &       + conjg(f00(4,ras(ne1,ndet1)))*fx00(1,ras(ne2,ndet2))&
     &       - conjg(fx00(1,ras(ne1,ndet1)))*f00(4,ras(ne2,ndet2))&
     &       - conjg(fx00(2,ras(ne1,ndet1)))*f00(3,ras(ne2,ndet2))&
     &       - conjg(fx00(3,ras(ne1,ndet1)))*f00(2,ras(ne2,ndet2))&
     &       - conjg(fx00(4,ras(ne1,ndet1)))*f00(1,ras(ne2,ndet2)))
          dstr0(2) = dstr0(2) + uI * (                            & 
     &         conjg(f00(1,ras(ne1,ndet1)))*fy00(4,ras(ne2,ndet2))&
     &       + conjg(f00(2,ras(ne1,ndet1)))*fy00(3,ras(ne2,ndet2))&
     &       + conjg(f00(3,ras(ne1,ndet1)))*fy00(2,ras(ne2,ndet2))&
     &       + conjg(f00(4,ras(ne1,ndet1)))*fy00(1,ras(ne2,ndet2))&
     &       - conjg(fy00(1,ras(ne1,ndet1)))*f00(4,ras(ne2,ndet2))&
     &       - conjg(fy00(2,ras(ne1,ndet1)))*f00(3,ras(ne2,ndet2))&
     &       - conjg(fy00(3,ras(ne1,ndet1)))*f00(2,ras(ne2,ndet2))&
     &       - conjg(fy00(4,ras(ne1,ndet1)))*f00(1,ras(ne2,ndet2)))
          dstr0(3) = dstr0(3) + uI * (                            &
     &         conjg(f00(1,ras(ne1,ndet1)))*fz00(4,ras(ne2,ndet2))&
     &       + conjg(f00(2,ras(ne1,ndet1)))*fz00(3,ras(ne2,ndet2))&
     &       + conjg(f00(3,ras(ne1,ndet1)))*fz00(2,ras(ne2,ndet2))&
     &       + conjg(f00(4,ras(ne1,ndet1)))*fz00(1,ras(ne2,ndet2))&
     &       - conjg(fz00(1,ras(ne1,ndet1)))*f00(4,ras(ne2,ndet2))&
     &       - conjg(fz00(2,ras(ne1,ndet1)))*f00(3,ras(ne2,ndet2))&
     &       - conjg(fz00(3,ras(ne1,ndet1)))*f00(2,ras(ne2,ndet2))&
     &       - conjg(fz00(4,ras(ne1,ndet1)))*f00(1,ras(ne2,ndet2)))
          dstr0(4) = dstr0(4) + (                                 &
     &         conjg(f00(1,ras(ne1,ndet1)))*fx00(4,ras(ne2,ndet2))& 
     &       - conjg(f00(2,ras(ne1,ndet1)))*fx00(3,ras(ne2,ndet2))&
     &       + conjg(f00(3,ras(ne1,ndet1)))*fx00(2,ras(ne2,ndet2))&
     &       - conjg(f00(4,ras(ne1,ndet1)))*fx00(1,ras(ne2,ndet2))&
     &       + conjg(fx00(4,ras(ne1,ndet1)))*f00(1,ras(ne2,ndet2))&
     &       - conjg(fx00(3,ras(ne1,ndet1)))*f00(2,ras(ne2,ndet2))&
     &       + conjg(fx00(2,ras(ne1,ndet1)))*f00(3,ras(ne2,ndet2))&
     &       - conjg(fx00(1,ras(ne1,ndet1)))*f00(4,ras(ne2,ndet2)))
          dstr0(5) = dstr0(5) + (                                 &
     &         conjg(f00(1,ras(ne1,ndet1)))*fy00(4,ras(ne2,ndet2))& 
     &       - conjg(f00(2,ras(ne1,ndet1)))*fy00(3,ras(ne2,ndet2))&
     &       + conjg(f00(3,ras(ne1,ndet1)))*fy00(2,ras(ne2,ndet2))&
     &       - conjg(f00(4,ras(ne1,ndet1)))*fy00(1,ras(ne2,ndet2))&
     &       + conjg(fy00(4,ras(ne1,ndet1)))*f00(1,ras(ne2,ndet2))&
     &       - conjg(fy00(3,ras(ne1,ndet1)))*f00(2,ras(ne2,ndet2))&
     &       + conjg(fy00(2,ras(ne1,ndet1)))*f00(3,ras(ne2,ndet2))&
     &       - conjg(fy00(1,ras(ne1,ndet1)))*f00(4,ras(ne2,ndet2)))
          dstr0(6) = dstr0(6) + (                                 &
     &         conjg(f00(1,ras(ne1,ndet1)))*fz00(4,ras(ne2,ndet2))& 
     &       - conjg(f00(2,ras(ne1,ndet1)))*fz00(3,ras(ne2,ndet2))&
     &       + conjg(f00(3,ras(ne1,ndet1)))*fz00(2,ras(ne2,ndet2))&
     &       - conjg(f00(4,ras(ne1,ndet1)))*fz00(1,ras(ne2,ndet2))&
     &       + conjg(fz00(4,ras(ne1,ndet1)))*f00(1,ras(ne2,ndet2))&
     &       - conjg(fz00(3,ras(ne1,ndet1)))*f00(2,ras(ne2,ndet2))&
     &       + conjg(fz00(2,ras(ne1,ndet1)))*f00(3,ras(ne2,ndet2))&
     &       - conjg(fz00(1,ras(ne1,ndet1)))*f00(4,ras(ne2,ndet2)))
          dstr0(7) = dstr0(7) + uI * (                            &
     &         conjg(f00(1,ras(ne1,ndet1)))*fx00(3,ras(ne2,ndet2))& 
     &       - conjg(f00(2,ras(ne1,ndet1)))*fx00(4,ras(ne2,ndet2))& 
     &       + conjg(f00(3,ras(ne1,ndet1)))*fx00(1,ras(ne2,ndet2))& 
     &       - conjg(f00(4,ras(ne1,ndet1)))*fx00(2,ras(ne2,ndet2))&
     &       - conjg(fx00(3,ras(ne1,ndet1)))*f00(1,ras(ne2,ndet2))&
     &       + conjg(fx00(4,ras(ne1,ndet1)))*f00(2,ras(ne2,ndet2))&
     &       - conjg(fx00(1,ras(ne1,ndet1)))*f00(3,ras(ne2,ndet2))&
     &       + conjg(fx00(2,ras(ne1,ndet1)))*f00(4,ras(ne2,ndet2)))
          dstr0(8) = dstr0(8) + uI * (                            &
     &         conjg(f00(1,ras(ne1,ndet1)))*fy00(3,ras(ne2,ndet2))& 
     &       - conjg(f00(2,ras(ne1,ndet1)))*fy00(4,ras(ne2,ndet2))& 
     &       + conjg(f00(3,ras(ne1,ndet1)))*fy00(1,ras(ne2,ndet2))& 
     &       - conjg(f00(4,ras(ne1,ndet1)))*fy00(2,ras(ne2,ndet2))&
     &       - conjg(fy00(3,ras(ne1,ndet1)))*f00(1,ras(ne2,ndet2))&
     &       + conjg(fy00(4,ras(ne1,ndet1)))*f00(2,ras(ne2,ndet2))&
     &       - conjg(fy00(1,ras(ne1,ndet1)))*f00(3,ras(ne2,ndet2))&
     &       + conjg(fy00(2,ras(ne1,ndet1)))*f00(4,ras(ne2,ndet2)))
          dstr0(9) = dstr0(9) + uI * (                            &
     &         conjg(f00(1,ras(ne1,ndet1)))*fz00(3,ras(ne2,ndet2))& 
     &       - conjg(f00(2,ras(ne1,ndet1)))*fz00(4,ras(ne2,ndet2))& 
     &       + conjg(f00(3,ras(ne1,ndet1)))*fz00(1,ras(ne2,ndet2))& 
     &       - conjg(f00(4,ras(ne1,ndet1)))*fz00(2,ras(ne2,ndet2))&
     &       - conjg(fz00(3,ras(ne1,ndet1)))*f00(1,ras(ne2,ndet2))&
     &       + conjg(fz00(4,ras(ne1,ndet1)))*f00(2,ras(ne2,ndet2))&
     &       - conjg(fz00(1,ras(ne1,ndet1)))*f00(3,ras(ne2,ndet2))&
     &       + conjg(fz00(2,ras(ne1,ndet1)))*f00(4,ras(ne2,ndet2)))
     
            dtspna = dtspna + (                                 &
     &        conjg(f00(1,ras(ne1,ndet1)))*f00(1,ras(ne2,ndet2))&
     &      + conjg(f00(3,ras(ne1,ndet1)))*f00(3,ras(ne2,ndet2)) )
            dtspnb = dtspnb + (                                 &
     &      + conjg(f00(2,ras(ne1,ndet1)))*f00(2,ras(ne2,ndet2))&
     &      + conjg(f00(4,ras(ne1,ndet1)))*f00(4,ras(ne2,ndet2)))

! added by hara

            dtspnx = dtspnx + (                                 &
     &        conjg(f00(1,ras(ne1,ndet1)))*f00(2,ras(ne2,ndet2))&
     &      + conjg(f00(2,ras(ne1,ndet1)))*f00(1,ras(ne2,ndet2))& 
     &      + conjg(f00(3,ras(ne1,ndet1)))*f00(4,ras(ne2,ndet2))&
     &      + conjg(f00(4,ras(ne1,ndet1)))*f00(3,ras(ne2,ndet2)) )

            dtspny = dtspny + uI * (                            &
     &        conjg(f00(2,ras(ne1,ndet1)))*f00(1,ras(ne2,ndet2))&
     &      - conjg(f00(1,ras(ne1,ndet1)))*f00(2,ras(ne2,ndet2))&
     &      + conjg(f00(4,ras(ne1,ndet1)))*f00(3,ras(ne2,ndet2))&
     &      - conjg(f00(3,ras(ne1,ndet1)))*f00(4,ras(ne2,ndet2)) )

!!!

         detdns = detdns + tmpdns
         detldns = detldns + tmpldns
         detsdns = detsdns + tmpsdns
         detztx = detztx + zetax0
         detzty = detzty + zetay0
         detztz = detztz + zetaz0
         dtchrl = dtchrl + dble(tpchrl)
         do kk = 1,9
           detstr(kk) = detstr(kk) + dstr0(kk)
         end do ! kk
  300      continue
         end do !ne1
         end do !ne2

!          write(*,*) sqrt(occdet(ndet1))*sqrt(occdet(ndet2))*detdns
          cofpdt = sqrt(occdet(ndet1))*sqrt(occdet(ndet2))
          dnsty = dnsty + cofpdt * detdns
          dkin = dkin + (-0.5d0)* cofpdt * dtdkin
          ldnsty = ldnsty + cofpdt * detldns
          sdnsty = sdnsty + cofpdt * detsdns
          spina0 = spina0 + cofpdt * dtspna
          spinb0 = spinb0 + cofpdt * dtspnb
          spinx0 = spinx0 + cofpdt * dtspnx
          spiny0 = spiny0 + cofpdt * dtspny
          chiral = chiral + ccc * cofpdt * dtchrl    ! ccc is added as in Ref. of THEOCHEM
          zetax = zetax - 0.5d0 * ccc * cofpdt * dble( detztx )
          zetay = zetay - 0.5d0 * ccc * cofpdt * dble( detzty )
          zetaz = zetaz - 0.5d0 * ccc * cofpdt * dble( detztz )
          do kk=1,9
            stress(kk) = stress(kk) + 0.5d0 * ccc * cofpdt * dble(detstr(kk))
          end do ! kk
        end do !ndet1
        end do !ndet2

          etaux = stress(6) - stress(8)
          etauy = stress(7) - stress(3)
          etauz = stress(2) - stress(4)

          spinc0 = spina0 - spinb0
          spina = dble(spina0)
          spinb = dble(spinb0)
          spinc = dble(spinc0)

          spinx = 0.5D0 * dble(spinx0)
          spiny = 0.5D0 * dble(spiny0)
          spinz = 0.5D0 * dble(spinc0)

          spin2 = spinx**2 + spiny**2 + spinz**2
          srspin = dsqrt(spin2)          

          write (11,1200) (stress(kk),kk=1,9)
          write (11,*) 'check-2'
          write (12,1000) w1,w2,w3,dnsty,chiral
          write (13,1000) w1,w2,w3,ldnsty,sdnsty
          write (16,1000) w1,w2,w3,dkin
          write (21,1000) w1,w2,w3,spina,spinb,spinc
          write (22,1005) w1,w2,w3,spinx,spiny,spinz,spin2,srspin
          write (31,1010) w1,w2,w3,zetax,zetay,zetaz
          write (32,1010) w1,w2,w3,etaux,etauy,etauz

          avgzeta = avgzeta + dsqrt(zetax**2 + zetay**2 + zetaz**2)
          avgtorque = avgtorque + dsqrt(etaux**2 + etauy**2 + etauz**2)
          abschiral = abschiral + dabs(chiral)
          avgspina = avgspina + spina
          avgspinb = avgspinb + spinb
          avgspinc = avgspinc + spinc

          avgspinx = avgspinx + spinx
          avgspiny = avgspiny + spiny
          avgspinz = avgspinz + spinz
          avgspin2 = avgspin2 + spin2
          avgsrspin = avgsrspin + srspin

           end do
         end do
         write (*,*) 1+ix+mix,'/',2*mix+1
       end do

      point = dble((2*mix+1)*(2*miy+1)*(2*miz+1))
      write (88,*) "avg. of zeta           : ", avgzeta/point
      write (88,*) "avg. of abs. of chiral : ", abschiral/point
      write (88,*) "avg. of torque         : ", avgtorque/point
      write (88,*) "avg. of a-spin dens    : ", avgspina/point
      write (88,*) "avg. of b-spin dens    : ", avgspinb/point
      write (88,*) "avg. of spin density   : ", avgspinc/point
      write (88,*) "avg. of x-spin dens    : ", avgspinx/point
      write (88,*) "avg. of y-spin dens    : ", avgspiny/point
      write (88,*) "avg. of z-spin dens    : ", avgspinz/point
      write (88,*) "avg. of spin2 density  : ", avgspin2/point
      write (88,*) "avg. of sroot-spin dens: ", avgsrspin/point

      call FDATE(today)
      write(*,*) 'Finish Calculation'
      write(*,*) ' '
      write(*,*) 'Outputs'
      write(*,*) 'fort.11 : stress (nine components)'
      write(*,*) 'fort.12 : x, y, z, rho_e, J_5^0'
      write(*,*) 'fort.13 : x, y, z, large_density, small_denstiy'
      write(*,*) 'fort.16 : x, y, z, dkin'
      write(*,*) 'fort.21 : x, y, z, alpha, beta, alpha - beta'
      write(*,*) 'fort.22 : x, y, z, spinx, spiny, spinz, spin^2, spin'
      write(*,*) 'fort.31 : x, y, z, zeta_x, zeta_y, zeta_z'
      write(*,*) 'fort.32 : x, y, z, torque_x, torque_y, torque_z'
      write(*,*) ' '
      call FDATE(today)
      write(*,*) 'Calculation end date   : ',today

      deallocate(enorb,occdet)
      deallocate(KRAM)
      deallocate(aa,n)
      deallocate(a)
      deallocate(xx,yy,zz)
      deallocate(f,fx,fy,fz)
      deallocate(f00,fx00,fy00,fz00)
      deallocate(ras)

      stop

      !write(*,*) 'fort.12 : x, y, z, rho_e, eigenvalue of  stress'
      
 1000     format (3f9.3,4e14.6)
 1005     format (3f9.3,5e14.6)
 1010     format (3f9.3,3e14.6)
 1020     format (3f9.3,4e14.6)
 1100     format (a12)
 1111     format (3f9.3,8f14.10)
 1200     format (9e10.3)
 1210     format (6e11.4)
 4100     format(a12)

end subroutine MRDFTCI

