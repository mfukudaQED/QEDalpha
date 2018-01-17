!============================================================================
! Last Change:28-May-2012.
!  module CI
!  subroutine CIread                      2012/04/16 written by Fukuda
!  subroutine convert_coef_MRDFTCI        2012/04/16 written by Fukuda
!  subroutine convert_coef_QEDtoMRDFT     2012/04/16 written by Fukuda
!  subroutine convert_pg_QEDtoMRDFT       2012/04/16 written by Fukuda
!============================================================================
module CI
contains
subroutine CIread(norb,NROOT,enorb,ras,KRAM,neact,detcnt,occdet)
!============================================================================
       implicit none

       integer norb,i,neact,NROOT,nrcnt,detcnt,j,k,l
       integer cnt
!       integer !,ras(20,100),ras2(20,100),ras3(20,100)
       integer,allocatable :: ras(:,:),ras2(:,:),ras3(:,:),ras1(:,:)
       double precision enorb(norb),occdet(1000)
       character(len=1) REPCH(14),KRAM(norb)
       character(len=14) MOJI
       character(len=4) cras
       integer num,nras1,nras2,nras3,tmpnras1,tmpnras3
       character(len=4) dummy

       detcnt = 1

       do
          read (41,4100) MOJI
          if (MOJI.eq.' Orbital Index') exit
       end do

       do i=1,norb
       KRAM(i) = 'F'
       read(41,1100) (REPCH(j),j=1,14),enorb(i)
       write(*,1100) (REPCH(j),j=1,14),enorb(i)
         do j=1,14
           if (REPCH(j).eq.'2') KRAM(i) ='T'
         end do
!         write(*,*) KRAM(i)
       end do

       do
          read (41,4100) MOJI
          if (MOJI.eq.' Generate natu') exit
       end do

       read (41,4100) MOJI
       read (41,1200) neact
       write(*,*)'Number of active electron :', neact

       do
          read (41,4100) MOJI
          if(MOJI.eq.' Allowed RAS c') exit
       end do

       read (41,4100) MOJI
       read (41,4100) MOJI
       read (41,"(3i6)") nras1,nras2,nras3
       write (*,*) 'NRAS1= ' , nras1
       write (*,*) 'NRAS2= ' , nras2
       write (*,*) 'NRAS3= ' , nras3
       if(neact.ne.(nras1+nras2+nras3)) then
          write(*,*) 'error nras123'
       end if

       nrcnt = 0
       do
          read (41,4100) MOJI
          if (MOJI.eq.'Determinant(s)') nrcnt = nrcnt + 1
          if (nrcnt.eq.NROOT) exit
       end do
       write(*,*) 'Determinant of Root :', nrcnt
       
       read (41,4100) MOJI
       read (41,"(20X,I3)") num
       write(*,*) 'num= ',num
!       allocate(ras(neact,num),ras2(neact,num),ras3(neact,num)
!     &          ,ras1(neact,num))
!       allocate(ras(100,num))
!       allocate(ras1(100,num),ras2(100,num),ras3(100,num))

       if(nras1.eq.0) then
          tmpnras1=1 !to avoid allocate error
       else
          tmpnras1=nras1
       end if
       if(nras3.eq.0) then
          tmpnras3=1 !to avoid allocate error
       else
          tmpnras3=nras3
       end if

       allocate(ras(neact,num))
       allocate(ras1(tmpnras1,num),ras2(nras2,num),ras3(tmpnras3,num))

       ras(:,:) =0
       ras1(:,:)=0
       ras2(:,:)=0
       ras3(:,:)=0

       do !read ras start
          read (41,1400) occdet(detcnt)
          write(*,*)'occdet', occdet(detcnt)
!          if (MOJI.eq.'Determinant   ') then
          if (occdet(detcnt).ne.0) then 
            read (41,"(9x,a4)") cras
            backspace(41)
            if (cras.eq."RAS1") then
              if (nras1.eq.0) then
                read(41,'()')
              else if (nras1.le.6) then
                read (41,1300) (ras1(j,detcnt),j=1,nras1)
              else
                read (41,1300) (ras1(j,detcnt),j=1,6)
                l = nras1 / 6
                if ( mod(nras1,6) .eq. 0) then 
                  do i=1,l-1
                    read (41,1310) (ras1(j,detcnt),j=6*i+1,6*i+6)
                  end do
                else if ( mod(nras1,6) .ne. 0) then
                  if (l .eq. 1) then
                    read (41,1310) (ras1(j,detcnt),j=7,nras1)
                  else
                    do i=1,l-1
                      read (41,1310) (ras1(j,detcnt),j=6*i+1,6*i+6)
                    end do
                  read (41,1310) (ras1(j,detcnt),j=6*i+7,nras1)
                  end if
                end if
              end if
            end if
            write (*,1350) (ras1(j,detcnt),j=1,tmpnras1)

            if(nras2.le.6) then
               read (41,1300) (ras2(j,detcnt),j=1,nras2)
            else
               read (41,1300) (ras2(j,detcnt),j=1,6)
            end if
            read (41,'(a4)') dummy
            backspace(41)
            if (dummy .ne. 'Elec') then
!              backspace(41)
              read (41,1310) (ras2(j,detcnt),j=7,nras2)
!              read (41,"(I4,6X,I4,6X,I4,6X,I4)")(ras2(j,detcnt),j=13,neact)
            end if 
            write (*,1350) (ras2(j,detcnt),j=1,nras2)          
            read (41,1350) (ras3(j,detcnt),j=1,tmpnras3)
            write (*,1350) (ras3(j,detcnt),j=1,tmpnras3)
            detcnt = detcnt + 1

          else 
            detcnt = detcnt -1
            write (*,*) 'Number of determinant :', detcnt
            exit
          end if

       end do !read ras end
       
       do j=1, detcnt
         cnt=1
         do i=1,tmpnras1
           if(ras1(i,j).ne.0) then
             ras(cnt,j) = ras1(i,j)
             cnt=cnt+1
           end if
         end do
         do i=1,nras2
           if(ras2(i,j).ne.0) then 
             ras(cnt,j) = ras2(i,j)
             cnt=cnt+1
           end if
         end do
         do i=1,tmpnras3
           if(ras3(i,j).ne.0) then
             ras(cnt,j) = ras3(i,j)
             cnt=cnt+1
           end if
         end do
       end do 

       do j=1, detcnt
       write(*,*) 'Determinant :',j,(ras(i,j),i=1,neact)
       write(*,*) 'Coefficient :',occdet(j)
       end do 

       deallocate(ras1,ras2,ras3)


!       do i=1,norb
!       write(*,1100) (REPCH(j),j=1,14),enorb(i)
!       end do

 1100 FORMAT(20X,14A,6X,G20.10)
 1200 format(36X,i4)
 1300 FORMAT(15X,6(I4,6X))
 1310 FORMAT(6(I4,6X))
 1350 FORMAT(15X,16(I4,6X)) 
 1400 FORMAT(43X,F16.8)
 4100 format(a14)

end subroutine CIread
end module CI

!============================================================================
subroutine convert_coef_MRDFTCI(a,KRAM,norb,enorb)
!============================================================================
  use DiracOutput
  implicit none

  integer i,j,k
  integer ck
  integer,intent(in) :: norb
  double precision,intent(in) :: enorb(norb)
  character(LEN=1),intent(in) :: KRAM(norb)
  character(LEN=1) KRAMF(NBS_E),KRAMT(NBS_E)
  complex(kind(0d0)),intent(inout) :: a(NBS0,4,norb)
  
  write(*,*)'# start convert coef'
  KRAMF(:) = 'T'
  KRAMT(:) = 'T'

  do i=1,norb
     ck=0
     do j=1,NBS_E
        if(ck.eq.0) then
           if(enorb(i).eq.e_eig(j)) then
              if(KRAM(i).eq.'F') then
                 if(KRAMF(j).eq.'T') then
                    call convert_coef_QEDtoMRDFT(a,KRAM(i),norb,i,j)
!!$                    do k=1,NBS_L
!!$                       a(k,1,i) = c_La(k,j)
!!$                       a(k,2,i) = c_Lb(k,j)
!!$                    end do
!!$                    do k=1,NBS_S
!!$                       a(k,3,i) = c_Sa(k,j)
!!$                       a(k,4,i) = c_Sb(k,j)
!!$                    end do
!!$                    if(NBS_L.gt.NBS_S) then
!!$                       do k=NBS_S+1,NBS0
!!$                             a(k,3,i) = (0.d0,0.d0)
!!$                             a(k,4,i) = (0.d0,0.d0)
!!$                       end do
!!$                    else
!!$                       do k=NBS_L+1,NBS0
!!$                             a(k,1,i) = (0.d0,0.d0)
!!$                             a(k,2,i) = (0.d0,0.d0)
!!$                       end do
!!$                    end if
                    KRAMF(j)='O'
                    ck=1
                 end if
              else if(KRAM(i).eq.'T') then !Kramers pair
                 if(KRAMT(j).eq.'T') then
                    call convert_coef_QEDtoMRDFT(a,KRAM(i),norb,i,j)
!!$                    do k=1,NBS_L
!!$                       a(k,1,i) = -conjg(c_Lb(k,j))
!!$                       a(k,2,i) =  conjg(c_La(k,j))
!!$                    end do
!!$                    do k=1,NBS_S
!!$                       a(k,3,i) = -conjg(c_Sb(k,j))
!!$                       a(k,4,i) =  conjg(c_Sa(k,j))
!!$                    end do
!!$                    if(NBS_L.gt.NBS_S) then
!!$                       do k=NBS_S+1,NBS0
!!$                             a(k,3,i) = (0.d0,0.d0)
!!$                             a(k,4,i) = (0.d0,0.d0)
!!$                       end do
!!$                    else
!!$                       do k=NBS_L+1,NBS0
!!$                             a(k,1,i) = (0.d0,0.d0)
!!$                             a(k,2,i) = (0.d0,0.d0)
!!$                       end do
!!$                    end if
                    KRAMT(j)='O'
                    ck=1
                 end if
              else
                 write(*,*)'error KRAM'
              end if
           end if
        end if
     end do
  end do

  do i=1,norb
     write(*,'(a5,i3,a1,a1)') 'KRAM(',i,')',KRAM(i)
  end do
  do j=1,NBS_E
     write(*,'(a6,i3,a1,a1)') 'KRAMF(',j,')',KRAMF(j)
     write(*,'(a6,i3,a1,a1)') 'KRAMT(',j,')',KRAMT(j)
  end do

  write(*,*)'end convert coef'

end subroutine convert_coef_MRDFTCI

!============================================================================
subroutine convert_coef_QEDtoMRDFT(a,kramers,norb,inumnorb,jnumnorb)
!============================================================================
  use DiracOutput
  implicit none

  integer i,k
  integer,intent(in) :: norb,inumnorb, jnumnorb
  character(LEN=1),intent(in) :: kramers
  complex(kind(0d0)),intent(inout) :: a(NBS0,4,norb)

!  write(*,*)'norb,inumnorb,jnumnorb',norb,inumnorb,jnumnorb
  if(kramers.eq.'F') then
     do i=1,NBS_L
        a(i,1,inumnorb) = c_La(i,jnumnorb)
        a(i,2,inumnorb) = c_Lb(i,jnumnorb)
        write(*,*)i,a(i,1,inumnorb),c_La(i,jnumnorb)
     end do
     do i=1,NBS_S
        a(i,3,inumnorb) = c_Sa(i,jnumnorb)
        a(i,4,inumnorb) = c_Sb(i,jnumnorb)
     end do

  else if(kramers.eq.'T') then
     do i=1,NBS_L
        a(i,1,inumnorb) = -conjg(c_Lb(i,jnumnorb))
        a(i,2,inumnorb) =  conjg(c_La(i,jnumnorb))
     end do
     do i=1,NBS_S
        a(i,3,inumnorb) = -conjg(c_Sb(i,jnumnorb))
        a(i,4,inumnorb) =  conjg(c_Sa(i,jnumnorb))
     end do
  else
     write(*,*)'error kramers in convert_QEDtoMRDFT'
  end if

  if(NBS_L.gt.NBS_S) then
     do i=NBS_S+1,NBS0
        do k=3,4
           a(i,k,inumnorb) = (0.d0,0.d0)
        end do
     end do
  else
     do i=NBS_L+1,NBS0
        do k=1,2
           a(i,k,inumnorb) = (0.d0,0.d0)
        end do
     end do
  end if

end subroutine convert_coef_QEDtoMRDFT

!============================================================================
subroutine convert_pg_QEDtoMRDFT(aa,n,xx,yy,zz)
!============================================================================
  use DiracOutput
  implicit none

  integer i,j,k
  integer,intent(out) :: n(NBS0,4)
  double precision,intent(out) :: aa(NBS0,4)
  double precision,intent(out) :: xx(NBS0,4),yy(NBS0,4),zz(NBS0,4)
  double precision :: bb(NBS0,4)

  do i=1,NBS_L
     do k=1,2
        aa(i,k) = aa_L(i)
        xx(i,k) = xx_L(i)
        yy(i,k) = yy_L(i)
        zz(i,k) = zz_L(i)
        call translaten2(nx_L(i),ny_L(i),nz_L(i),n(i,k))
     end do
  end do
  do i=1,NBS_S
     do k=3,4
        aa(i,k) = aa_S(i)
        xx(i,k) = xx_S(i)
        yy(i,k) = yy_S(i)
        zz(i,k) = zz_S(i)
        call translaten2(nx_S(i),ny_S(i),nz_S(i),n(i,k))
     end do
  end do

  if(NBS_L.gt.NBS_S) then
     do i=NBS_S+1,NBS0
        do k=3,4
           xx(i,k) = 0.0D0
           yy(i,k) = 0.0D0
           zz(i,k) = 0.0D0
           aa(i,k) = 1.D0
           n(i,k)  = -1.D0
        end do
     end do
  else
     do i=NBS_L+1,NBS0
        do k=1,2
           xx(i,k) = 0.0D0
           yy(i,k) = 0.0D0
           zz(i,k) = 0.0D0
           aa(i,k) = 1.D0
           n(i,k)  = -1.D0
        end do
     end do
  end if

  ! normalization for MRDFT
  bb(:,:)=1.d0
  call normal(aa,bb,n,NBS0)

!      do j=1,NBS0
!         do k=1,4
!            write(*,*)j,k,bb(j,k)
!         end do
!      end do
!      stop

  do j=1,NBS_E
     do k=1,NBS_L
        c_La(k,j) = c_La(k,j)*bb(k,1)
        c_Lb(k,j) = c_Lb(k,j)*bb(k,2)
     end do
     do k=1,NBS_S
        c_Sa(k,j) = c_Sa(k,j)*bb(k,3)
        c_Sb(k,j) = c_Sb(k,j)*bb(k,4)
     end do
  end do

end subroutine convert_pg_QEDtoMRDFT
