! Last Change:18-Nov-2012.
!==============================================================================
!   subroutine readVectors version 2011/05/30 written by Fukuda
!      2011.8.25  edited by fukuda to calculate contract basis ver.
!   subroutine angular
!   subroutine translaten 2012/04/16 written by Fukuda
!   subroutine translaten2 2012/04/16 written by Fukuda 
!   subroutine translaten3 2012/04/16 written by Fukuda
!   function renorm 2012/04/16 written by Fukuda
!   function renorm1 2012/04/16 written by Fukuda
!==============================================================================
!==============================================================================
subroutine readVectors(bnum_L,bnum_S,bnum_LT,bnum_ST,bnum_LTA,bnum_STA,cpg_L,cpg_S,n,cont)
!==============================================================================
  use DiracOutput

  implicit none

  character(len=3) MOJI,L,S
  character(len=4) Angu
  character(len=5) DUMMY
  character(len=6) VECTOR
  integer,intent(in):: bnum_L(100),bnum_S(100)
  integer,intent(in):: bnum_LT,bnum_ST
  integer,intent(in):: bnum_LTA(NAT),bnum_STA(NAT)
  integer i,j,ip,ie,kk,k0,bnum,num,inat,m(NAT),m_total
  integer nt,Pns(NBS0,4),Ens(NBS0,4)  !---------------checker of n(NBS0,4) 
  integer tmpPns(NAT,NBS0,4),tmpEns(NAT,NBS0,4)
  integer tmpcont
  integer mcont(NAT)
  integer,intent(in):: cont(NBS00,NAT)
  integer,intent(in):: n(NBS0,4)
  double precision eig
  double precision,intent(in):: cpg_L(NBS_L),cpg_S(NBS_S)
  double precision renorm
  complex(kind(0d0)) tmpA, tmpB
  complex(kind(0d0)) tmpd_La(NAT,NBS_L), tmpd_Lb(NAT,NBS_L),tmpd_Sa(NAT,NBS_S), tmpd_Sb(NAT,NBS_S)
  complex(kind(0d0)) tmpc_La(NAT,NBS_L), tmpc_Lb(NAT,NBS_L),tmpc_Sa(NAT,NBS_S), tmpc_Sb(NAT,NBS_S)
  integer tmpNBS


  read(42,'(a6)') VECTOR
  if (VECTOR.eq.'VECTOR') then
   write(*,*) 'start reading Vectors.txt'
  end if

 !----------------single atom--------------------------------------------------
  if (SYMM.eq.0) then
    write(*,*)"SYMM = 0 is not available."
    write(*,*)"SYMM = 0 is symmetry calculation for atom"
    stop
!!$
!!$   ip=0 ! counter of positronic orbitals
!!$   ie=0 ! counter of electronic orbitals
!!$
!!$   do kk=1,NBS*2
!!$    read (42,'(A3,5x,E18.10)') MOJI,eig
!!$     if (MOJI.eq.'PIV') then
!!$      ip=ip+1
!!$      p_eig(ip) = eig
!!$
!!$      read(42,'(A5)')DUMMY
!!$      do i=1,NBS_L
!!$        read (42,'(i4,a3,7x,a4,2x,4e22.14)') bnum, L, Angu, tmpA, tmpB
!!$!        write (*,'(i4,a3,7x,a4,2x,4e18.10)') bnum, L, Angu, tmpA, tmpB
!!$         if(bnum.ne.i) then
!!$          write(*,*)'error read vectors PIV Large'
!!$          write(*,*)'i=',i,'bnum=',bnum
!!$          stop
!!$         end if
!!$        call angular(nt, Angu)
!!$        Pns(i,1) = nt
!!$        Pns(i,2) = nt
!!$!        d_La(i,ip) = tmpA*cpg_L(i)*renorm(nx_L(i),ny_L(i),nz_L(i))   !---------coefficient of La
!!$        d_La(i,ip) = tmpA*cpg_L(i)   !---------coefficient of La
!!$!        d_Lb(i,ip) = tmpB*cpg_L(i)*renorm(nx_L(i),ny_L(i),nz_L(i))   !---------coefficient of Lb
!!$        d_Lb(i,ip) = tmpB*cpg_L(i)   !---------coefficient of Lb
!!$!!$!test start
!!$!!$        if(n(i,1).eq.4) then
!!$!!$         d_La(i,ip)=d_La(i,ip)*(3**(0.5))
!!$!!$         d_Lb(i,ip)=d_Lb(i,ip)*(3**(0.5))
!!$!!$        else if(n(i,1).eq.7) then
!!$!!$         d_La(i,ip)=d_La(i,ip)*(3**(0.5))
!!$!!$         d_Lb(i,ip)=d_Lb(i,ip)*(3**(0.5))
!!$!!$        else if(n(i,1).eq.9) then
!!$!!$         d_La(i,ip)=d_La(i,ip)*(3**(0.5))
!!$!!$         d_Lb(i,ip)=d_Lb(i,ip)*(3**(0.5))
!!$!!$        end if
!!$!!$!test end
!!$      end do
!!$
!!$      do i=1,NBS_S
!!$        read (42,'(i4,a3,7x,a4,2x,4e22.14)') bnum, S, Angu, tmpA, tmpB
!!$!        write (*,'(i4,a3,7x,a4,2x,4e18.10)') bnum, S, Angu, tmpA, tmpB
!!$         if(bnum /= i+NBS_L) then
!!$          write(*,*)'error read vectors PIV Small'
!!$          write(*,*)'i+NBS_L=',i+NBS_L,'bnum=',bnum
!!$          stop
!!$         end if
!!$        call angular(nt, Angu)
!!$        Pns(i,3) = nt
!!$        Pns(i,4) = nt
!!$!        d_Sa(i,ip) = tmpA*cpg_S(i)*renorm(nx_S(i),ny_S(i),nz_S(i))    !---------coefficient of Sa
!!$        d_Sa(i,ip) = tmpA*cpg_S(i)    !---------coefficient of Sa
!!$!        d_Sb(i,ip) = tmpB*cpg_S(i)*renorm(nx_S(i),ny_S(i),nz_S(i))    !---------coefficient of Sb
!!$        d_Sb(i,ip) = tmpB*cpg_S(i)    !---------coefficient of Sb
!!$!!$!test start
!!$!!$        if(n(i,3).eq.4) then
!!$!!$         d_Sa(i,ip)=d_Sa(i,ip)*(3**(0.5))
!!$!!$         d_Sb(i,ip)=d_Sb(i,ip)*(3**(0.5))
!!$!!$        else if(n(i,3).eq.7) then
!!$!!$         d_Sa(i,ip)=d_Sa(i,ip)*(3**(0.5))
!!$!!$         d_Sb(i,ip)=d_Sb(i,ip)*(3**(0.5))
!!$!!$        else if(n(i,3).eq.9) then
!!$!!$         d_Sa(i,ip)=d_Sa(i,ip)*(3**(0.5))
!!$!!$         d_Sb(i,ip)=d_Sb(i,ip)*(3**(0.5))
!!$!!$        end if
!!$!!$!test end
!!$      end do
!!$
!!$
!!$     else if (MOJI.eq.'EIV') then
!!$      ie=ie+1
!!$      e_eig(ie) = eig
!!$
!!$      read(42,'(A5)')DUMMY
!!$      do i=1,NBS_L
!!$        read (42,'(i4,a3,7x,a4,2x,4e22.14)') bnum, L, Angu, tmpA, tmpB
!!$!        write (*,'(i4,a3,7x,a4,2x,4e22.14)') bnum, L, Angu, tmpA, tmpB
!!$         if(bnum.ne.i) then
!!$          write(*,*)'error read vectors EIV Large'
!!$          write(*,*)'i=',i,'bnum=',bnum
!!$          stop
!!$         end if
!!$        call angular(nt, Angu)
!!$        Ens(i,1) = nt
!!$        Ens(i,2) = nt
!!$!        c_La(i,ie) = tmpA*cpg_L(i)*renorm(nx_L(i),ny_L(i),nz_L(i))   !---------coefficient of La
!!$        c_La(i,ie) = tmpA*cpg_L(i)   !---------coefficient of La
!!$!        c_Lb(i,ie) = tmpB*cpg_L(i)*renorm(nx_L(i),ny_L(i),nz_L(i))   !---------coefficient of Lb
!!$        c_Lb(i,ie) = tmpB*cpg_L(i)   !---------coefficient of Lb
!!$!!$!test start
!!$!!$        if(n(i,1).eq.4) then
!!$!!$         c_La(i,ie)=c_La(i,ie)*(3**(0.5))
!!$!!$         c_Lb(i,ie)=c_Lb(i,ie)*(3**(0.5))
!!$!!$        else if(n(i,1).eq.7) then
!!$!!$         c_La(i,ie)=c_La(i,ie)*(3**(0.5))
!!$!!$         c_Lb(i,ie)=c_Lb(i,ie)*(3**(0.5))
!!$!!$        else if(n(i,1).eq.9) then
!!$!!$         c_La(i,ie)=c_La(i,ie)*(3**(0.5))
!!$!!$         c_Lb(i,ie)=c_Lb(i,ie)*(3**(0.5))
!!$!!$        end if
!!$!!$!test end
!!$      end do
!!$
!!$      do i=1,NBS_S
!!$        read (42,'(i4,a3,7x,a4,2x,4e22.14)') bnum, S, Angu, tmpA, tmpB
!!$!        write (*,'(i4,a3,7x,a4,2x,4e22.14)') bnum, S, Angu, tmpA, tmpB
!!$         if(bnum /= i+NBS_L) then
!!$          write(*,*)'error read vectors EIV Small'
!!$          write(*,*)'i+NBS_L=',i+NBS_L,'bnum=',bnum
!!$          stop
!!$         end if
!!$        call angular(nt, Angu)
!!$        Ens(i,3) = nt
!!$        Ens(i,4) = nt
!!$!        c_Sa(i,ie) = tmpA*cpg_S(i)*renorm(nx_S(i),ny_S(i),nz_S(i))    !---------coefficient of Sa
!!$        c_Sa(i,ie) = tmpA*cpg_S(i)    !---------coefficient of Sa
!!$!        c_Sb(i,ie) = tmpB*cpg_S(i)*renorm(nx_S(i),ny_S(i),nz_S(i))    !---------coefficient of Sb
!!$        c_Sb(i,ie) = tmpB*cpg_S(i)    !---------coefficient of Sb
!!$!!$!test start
!!$!!$        if(n(i,3).eq.4) then
!!$!!$         c_Sa(i,ie)=c_Sa(i,ie)*(3**(0.5))
!!$!!$         c_Sb(i,ie)=c_Sb(i,ie)*(3**(0.5))
!!$!!$        else if(n(i,3).eq.7) then
!!$!!$         c_Sa(i,ie)=c_Sa(i,ie)*(3**(0.5))
!!$!!$         c_Sb(i,ie)=c_Sb(i,ie)*(3**(0.5))
!!$!!$        else if(n(i,3).eq.9) then
!!$!!$         c_Sa(i,ie)=c_Sa(i,ie)*(3**(0.5))
!!$!!$         c_Sb(i,ie)=c_Sb(i,ie)*(3**(0.5))
!!$!!$        end if
!!$!!$!test end
!!$      end do
!!$ 
!!$     end if
!!$   end do
  !-------------------------molecule-------------------------------------------
  else    !--------SYMM=1 => molecule using symmmetry calculation, SYMM=2 => molecule using nonsymmmetry calculation
   ip=0 ! counter of positronic orbitals
   ie=0 ! counter of electronic orbitals

   if(KBAL.eq.0) tmpNBS = NBS_TOTAL
   if(KBAL.eq.1) tmpNBS = NBS*2
   write(*,*)'tmpNBS=',tmpNBS

   do kk=1,tmpNBS
!    read (42,'(A3,5x,E18.10)') MOJI,eig
    read (42,'(A3,5x,E22.14)') MOJI,eig
 !-------read PIV---------------------------------------------------------------------------
     if (MOJI.eq.'PIV') then
      ip=ip+1
      p_eig(ip) = eig

      do inat=1,NAT
        m(inat)=0  ! counter of (inat)th atom's primitive gaussian number
        mcont(inat)=0 ! counter of (inat)th atom's contracted primitive gaussian number
      end do
      read(42,'(A5)')DUMMY
 !-------------read large part of PIV--------------------------------------------------------
      do i=1,bnum_LT
      
        if(SYMM.eq.1) then
           read (42,'(i4,a3,4x,i3,a4,2x,4e22.14)') bnum, L, num, Angu, tmpA, tmpB
        else if(SYMM.eq.2) then
           read (42,'(i4,a3,3x,i1,3x,a4,2x,4e22.14)') bnum, L, num, Angu, tmpA, tmpB
!           write (*,'(i4,a3,3x,i1,3x,a4,2x,4e22.14)') bnum, L, num, Angu, tmpA, tmpB
        else if(SYMM.eq.0) then
           read (42,'(i4,a3,3x,i1,3x,a4,2x,4e22.14)') bnum, L, num, Angu, tmpA, tmpB
        end if

         if(bnum.ne.i) then
            write(*,*)'error read vectors PIV Large'
            write(*,*)'i=',i,'bnum=',bnum
            stop
         end if

        call angular(nt, Angu)

        do inat=1,NAT
          if(num.eq.inat) then
            mcont(inat) = mcont(inat) + 1
            tmpcont = cont(mcont(inat),inat)
            do j=1,tmpcont
              m(inat) = m(inat)+1 
              tmpPns(inat,m(inat),1) = nt
              tmpPns(inat,m(inat),2) = nt
              tmpd_La(inat,m(inat)) = tmpA
              tmpd_Lb(inat,m(inat)) = tmpB
	    end do
          end if
        end do

      end do

         m_total=0
       do inat=1,NAT
          do i=1,m(inat)
           Pns(i+m_total,1) = tmpPns(inat,i,1)
           Pns(i+m_total,2) = tmpPns(inat,i,2)
           d_La(i+m_total,ip) = tmpd_La(inat,i)*cpg_L(i+m_total)*renorm(nx_L(i+m_total),ny_L(i+m_total),nz_L(i+m_total))   !---------coefficient of La
!           d_La(i+m_total,ip) = tmpd_La(inat,i)*cpg_L(i+m_total)   !---------coefficient of La
           d_Lb(i+m_total,ip) = tmpd_Lb(inat,i)*cpg_L(i+m_total)*renorm(nx_L(i+m_total),ny_L(i+m_total),nz_L(i+m_total))   !---------coefficient of La
!           d_Lb(i+m_total,ip) = tmpd_Lb(inat,i)*cpg_L(i+m_total)   !---------coefficient of Lb
          end do
         m_total=m_total+m(inat)
       end do
 !--------------read small part of PIV---------------------------------------------------------
      do inat=1,NAT
        m(inat)=0  ! counter of (inat)th atom's primitive gaussian number
        mcont(inat)=0 ! counter of (inat)th atom's contracted primitive gaussian number
      end do

      do i=1,bnum_ST
        if(SYMM.eq.1) then
          read (42,'(i4,a3,4x,i3,a4,2x,4e22.14)') bnum, S, num, Angu, tmpA, tmpB
!          write(*,'(i4,a3,4x,i3,a4,2x,4e22.14)') bnum, S, num, Angu, tmpA, tmpB
        else if(SYMM.eq.2) then
          read (42,'(i4,a3,3x,i1,3x,a4,2x,4e22.14)') bnum, S, num, Angu, tmpA, tmpB
!          write (*,'(i4,a3,3x,i1,3x,a4,2x,4e22.14)') bnum, S, num, Angu, tmpA, tmpB
        else if(SYMM.eq.0) then
          read (42,'(i4,a3,3x,i1,3x,a4,2x,4e22.14)') bnum, S, num, Angu, tmpA, tmpB
        end if
         if(bnum.ne.i+bnum_LT) then
          write(*,*)'error read vectors PIV Small'
          write(*,*)'i+bnum_LT=',i+bnum_LT,'bnum=',bnum
          stop
         end if
        call angular(nt, Angu)
        do inat=1,NAT
          if(num.eq.inat) then
            mcont(inat) = mcont(inat) + 1
            tmpcont = cont(mcont(inat)+bnum_LTA(inat),inat)
            do j=1,tmpcont
              m(inat) = m(inat)+1 
              tmpPns(inat,m(inat),3) = nt
              tmpPns(inat,m(inat),4) = nt
              tmpd_Sa(inat,m(inat)) = tmpA
              tmpd_Sb(inat,m(inat)) = tmpB
	    end do
          end if
        end do
      end do

         m_total=0
       do inat=1,NAT
          do i=1,m(inat)
           Pns(i+m_total,3) = tmpPns(inat,i,3)
           Pns(i+m_total,4) = tmpPns(inat,i,4)
           d_Sa(i+m_total,ip) = tmpd_Sa(inat,i)*cpg_S(i+m_total)*renorm(nx_S(i+m_total),ny_S(i+m_total),nz_S(i+m_total))   !---------coefficient of La
!           d_Sa(i+m_total,ip) = tmpd_Sa(inat,i)*cpg_S(i+m_total)   !---------coefficient of La
           d_Sb(i+m_total,ip) = tmpd_Sb(inat,i)*cpg_S(i+m_total)*renorm(nx_S(i+m_total),ny_S(i+m_total),nz_S(i+m_total))   !---------coefficient of Lb
!           d_Sb(i+m_total,ip) = tmpd_Sb(inat,i)*cpg_S(i+m_total)   !---------coefficient of Lb
          end do
         m_total=m_total+m(inat)
       end do
 !-------read EIV---------------------------------------------------------------------------
     else if (MOJI.eq.'EIV') then
      ie=ie+1
      e_eig(ie) = eig

      do inat=1,NAT
        m(inat)=0  ! counter of (inat)th atom's primitive gaussian number
        mcont(inat)=0 ! counter of (inat)th atom's contracted primitive gaussian number
      end do

      read(42,'(A5)')DUMMY
 !---------------read large part of EIV--------------------------------------------

      do i=1,bnum_LT
        if(SYMM.eq.1) then
          read (42,'(i4,a3,4x,i3,a4,2x,4e22.14)') bnum, L, num, Angu, tmpA, tmpB
        else if(SYMM.eq.2) then
          read (42,'(i4,a3,3x,i1,3x,a4,2x,4e22.14)') bnum, L, num, Angu, tmpA, tmpB
!        write (*,'(i4,a3,3x,i1,3x,a4,2x,4e22.14)') bnum, L, num, Angu, tmpA, tmpB
        else if(SYMM.eq.0) then
          read (42,'(i4,a3,3x,i1,3x,a4,2x,4e22.14)') bnum, L, num, Angu, tmpA, tmpB
        end if
         if(bnum.ne.i) then
          write(*,*)'error read vectors EIV Large'
          write(*,*)'i=',i,'bnum=',bnum
          stop
         end if
        call angular(nt, Angu)

        do inat=1,NAT
          if(num.eq.inat) then
            mcont(inat) = mcont(inat) + 1
            tmpcont = cont(mcont(inat),inat)
            do j=1,tmpcont
              m(inat) = m(inat)+1 
              tmpEns(inat,m(inat),1) = nt
              tmpEns(inat,m(inat),2) = nt
              tmpc_La(inat,m(inat)) = tmpA
              tmpc_Lb(inat,m(inat)) = tmpB
	    end do
          end if
        end do
      end do

         m_total=0
       do inat=1,NAT
          do i=1,m(inat)
           Ens(i+m_total,1) = tmpEns(inat,i,1)
           Ens(i+m_total,2) = tmpEns(inat,i,2)
           c_La(i+m_total,ie) = tmpc_La(inat,i)*cpg_L(i+m_total)*renorm(nx_L(i+m_total),ny_L(i+m_total),nz_L(i+m_total))   !---------coefficient of La
!           c_La(i+m_total,ie) = tmpc_La(inat,i)*cpg_L(i+m_total)   !---------coefficient of La
           c_Lb(i+m_total,ie) = tmpc_Lb(inat,i)*cpg_L(i+m_total)*renorm(nx_L(i+m_total),ny_L(i+m_total),nz_L(i+m_total))   !---------coefficient of Lb
!           c_Lb(i+m_total,ie) = tmpc_Lb(inat,i)*cpg_L(i+m_total)   !---------coefficient of Lb
          end do
         m_total=m_total+m(inat)
       end do

 !--------------read small part of EIV---------------------------------------------------------
      do inat=1,NAT
        m(inat)=0  ! counter of (inat)th atom's primitive gaussian number
        mcont(inat)=0 ! counter of (inat)th atom's contracted primitive gaussian number
      end do

      do i=1,bnum_ST
        if(SYMM.eq.1) then
        read (42,'(i4,a3,4x,i3,a4,2x,4e22.14)') bnum, S, num, Angu, tmpA, tmpB
        else if(SYMM.eq.2) then
        read (42,'(i4,a3,3x,i1,3x,a4,2x,4e22.14)') bnum, S, num, Angu, tmpA, tmpB
        else if(SYMM.eq.0) then
        read (42,'(i4,a3,3x,i1,3x,a4,2x,4e22.14)') bnum, S, num, Angu, tmpA, tmpB
        end if
         if(bnum.ne.i+bnum_LT) then
          write(*,*)'error read vectors EIV Small'
          write(*,*)'i+bnum_LT=',i+bnum_LT,'bnum=',bnum
          stop
         end if
        call angular(nt, Angu)

        do inat=1,NAT
          if(num.eq.inat) then
            mcont(inat) = mcont(inat) + 1
            tmpcont = cont(mcont(inat)+bnum_LTA(inat),inat)
            do j=1,tmpcont
              m(inat) = m(inat)+1 
              tmpEns(inat,m(inat),3) = nt
              tmpEns(inat,m(inat),4) = nt
              tmpc_Sa(inat,m(inat)) = tmpA
              tmpc_Sb(inat,m(inat)) = tmpB
	    end do
          end if
        end do !inat
      end do

         m_total=0
       do inat=1,NAT
          do i=1,m(inat)
           Ens(i+m_total,3) = tmpEns(inat,i,3)
           Ens(i+m_total,4) = tmpEns(inat,i,4)
           c_Sa(i+m_total,ie) = tmpc_Sa(inat,i)*cpg_S(i+m_total)*renorm(nx_S(i+m_total),ny_S(i+m_total),nz_S(i+m_total))   !---------coefficient of La
!           c_Sa(i+m_total,ie) = tmpc_Sa(inat,i)*cpg_S(i+m_total)   !---------coefficient of La
           c_Sb(i+m_total,ie) = tmpc_Sb(inat,i)*cpg_S(i+m_total)*renorm(nx_S(i+m_total),ny_S(i+m_total),nz_S(i+m_total))   !---------coefficient of Lb
!           c_Sb(i+m_total,ie) = tmpc_Sb(inat,i)*cpg_S(i+m_total)   !---------coefficient of Lb
          end do
         m_total=m_total+m(inat)
       end do

     end if
   end do
  end if  !SYMM

!--------------error check-------------------------
  do i=1,NBS_L
     do k0=1,2
        if(Pns(i,k0).ne.n(i,k0)) then
           write(*,*)'error Pns(',i,',',k0,')=',Pns(i,k0),'n=',n(i,k0)
           stop
        end if
        if(Ens(i,k0).ne.n(i,k0)) then
           write(*,*)'error Ens(',i,',',k0,')=',Ens(i,k0),'n=',n(i,k0)
           stop
        end if
     end do
  end do

  do i=1,NBS_S
     do k0=3,4
        if(Pns(i,k0).ne.n(i,k0)) then
           write(*,*)'error Pns(',i,',',k0,')=',Pns(i,k0),'n=',n(i,k0)
           stop
        end if
        if(Ens(i,k0).ne.n(i,k0)) then
           write(*,*)'error Ens(',i,',',k0,')=',Ens(i,k0),'n=',n(i,k0)
           stop
        end if
     end do
  end do
!-------------arrange orbitals to eigenvalue--------
  call arrangeEigenvalueE(e_eig)
  call arrangeEigenvalueP(p_eig)
!  call arrangeEigenvalue(e_eig,p_eig)
  write(*,*)'finish reading vectors.txt'

end subroutine readVectors

!==============================================================================
subroutine angular(ns, Angu)
!==============================================================================
      implicit none
      integer,intent(out):: ns
      character(len=4),intent(in):: Angu


      if (Angu .eq. 's   ') then
         ns =  0
      else if (Angu .eq. 'px  ') then
         ns =  1
      else if (Angu .eq. 'py  ') then
         ns =  2
      else if (Angu .eq. 'pz  ') then
         ns =  3
      else if (Angu .eq. 'dxx ') then
         ns =  4
      else if (Angu .eq. 'dxy ') then
         ns =  5
      else if (Angu .eq. 'dxz ') then
         ns =  6
      else if (Angu .eq. 'dyy ') then
         ns =  7
      else if (Angu .eq. 'dyz ') then
         ns =  8
      else if (Angu .eq. 'dzz ') then
         ns =  9
      else if (Angu .eq. 'fxxx') then
         ns = 10
      else if (Angu .eq. 'fxxy') then
         ns = 11
      else if (Angu .eq. 'fxxz') then
         ns = 12
      else if (Angu .eq. 'fxyy') then
         ns = 13
      else if (Angu .eq. 'fxyz') then
         ns = 14
      else if (Angu .eq. 'fxzz') then
         ns = 15
      else if (Angu .eq. 'fyyy') then
         ns = 16
      else if (Angu .eq. 'fyyz') then
         ns = 17
      else if (Angu .eq. 'fyzz') then
         ns = 18
      else if (Angu .eq. 'fzzz') then
         ns = 19
      else if (Angu .eq. 'g400') then
         ns = 20
      else if (Angu .eq. 'g310') then
         ns = 21
      else if (Angu .eq. 'g301') then
         ns = 22
      else if (Angu .eq. 'g220') then
         ns = 23
      else if (Angu .eq. 'g211') then
         ns = 24
      else if (Angu .eq. 'g202') then
         ns = 25
      else if (Angu .eq. 'g130') then
         ns = 26
      else if (Angu .eq. 'g121') then
         ns = 27
      else if (Angu .eq. 'g112') then
         ns = 28
      else if (Angu .eq. 'g103') then
         ns = 29
      else if (Angu .eq. 'g040') then
         ns = 30
      else if (Angu .eq. 'g031') then
         ns = 31
      else if (Angu .eq. 'g022') then
         ns = 32
      else if (Angu .eq. 'g013') then
         ns = 33
      else if (Angu .eq. 'g004') then
         ns = 34
      else if (Angu .eq. 'h500') then
         ns = 35
      else if (Angu .eq. 'h410') then
         ns = 36
      else if (Angu .eq. 'h401') then
         ns = 37
      else if (Angu .eq. 'h320') then
         ns = 38
      else if (Angu .eq. 'h311') then
         ns = 39
      else if (Angu .eq. 'h302') then
         ns = 40
      else if (Angu .eq. 'h230') then
         ns = 41
      else if (Angu .eq. 'h221') then
         ns = 42
      else if (Angu .eq. 'h212') then
         ns = 43
      else if (Angu .eq. 'h203') then
         ns = 44
      else if (Angu .eq. 'h140') then
         ns = 45
      else if (Angu .eq. 'h131') then
         ns = 46
      else if (Angu .eq. 'h122') then
         ns = 47
      else if (Angu .eq. 'h113') then
         ns = 48
      else if (Angu .eq. 'h104') then
         ns = 49
      else if (Angu .eq. 'h050') then
         ns = 50
      else if (Angu .eq. 'h041') then
         ns = 51
      else if (Angu .eq. 'h032') then
         ns = 52
      else if (Angu .eq. 'h023') then
         ns = 53
      else if (Angu .eq. 'h014') then
         ns = 54
      else if (Angu .eq. 'h005') then
         ns = 55
      else if (Angu .eq. 'i600') then
         ns = 56
      else if (Angu .eq. 'i510') then
         ns = 57
      else if (Angu .eq. 'i501') then
         ns = 58
      else if (Angu .eq. 'i420') then
         ns = 59
      else if (Angu .eq. 'i411') then
         ns = 60
      else if (Angu .eq. 'i402') then
         ns = 61
      else if (Angu .eq. 'i330') then
         ns = 62
      else if (Angu .eq. 'i321') then
         ns = 63
      else if (Angu .eq. 'i312') then
         ns = 64
      else if (Angu .eq. 'i303') then
         ns = 65
      else if (Angu .eq. 'i240') then
         ns = 66
      else if (Angu .eq. 'i231') then
         ns = 67
      else if (Angu .eq. 'i222') then
         ns = 68
      else if (Angu .eq. 'i213') then
         ns = 69
      else if (Angu .eq. 'i204') then
         ns = 70
      else if (Angu .eq. 'i150') then
         ns = 71
      else if (Angu .eq. 'i141') then
         ns = 72
      else if (Angu .eq. 'i132') then
         ns = 73
      else if (Angu .eq. 'i123') then
         ns = 74
      else if (Angu .eq. 'i114') then
         ns = 75
      else if (Angu .eq. 'i105') then
         ns = 76
      else if (Angu .eq. 'i060') then
         ns = 77
      else if (Angu .eq. 'i051') then
         ns = 78
      else if (Angu .eq. 'i042') then
         ns = 79
      else if (Angu .eq. 'i033') then
         ns = 80
      else if (Angu .eq. 'i024') then
         ns = 81
      else if (Angu .eq. 'i015') then
         ns = 82
      else if (Angu .eq. 'i006') then
         ns = 83
      else 
         write (*,*) 'Error subroutine angular :',Angu
         stop 
      end if
      return
end subroutine angular

!==============================================================================
subroutine translaten(n)
!==============================================================================
  use DiracOutput
      implicit none
      integer j,k0,ntmp
      integer,intent(in):: n(NBS0,4)
      integer nx(NBS0,4),ny(NBS0,4),nz(NBS0,4)

      do k0=1,3,2
      do j=1,NBS0
         if(n(j,k0).ge.0) then
            call translaten3(nx(j,k0),ny(j,k0),nz(j,k0),n(j,k0))
!            write(*,*)j,nx(j,k0),ny(j,k0),nz(j,k0),n(j,k0)
            nx(j,k0+1)=nx(j,k0)
            ny(j,k0+1)=ny(j,k0)
            nz(j,k0+1)=nz(j,k0)
         end if
      end do
      end do

!!$      do k0=1,3,2
!!$        do j=1,NBS0
!!$          call translaten2(nx(j,k0),ny(j,k0),nz(j,k0),ntmp)
!!$          write(*,'(a,6I3)')'j,k0,n,nx,ny,nz',j,k0,ntmp,nx(j,k0),ny(j,k0),nz(j,k0)
!!$        end do
!!$      end do
!!$      stop

!     write(*,*) 'j,nx_L(j),ny_L(j),nz_L(j)'
    do j=1,NBS_L
     nx_L(j)=nx(j,1)
     ny_L(j)=ny(j,1)
     nz_L(j)=nz(j,1)
!     write(*,*) j,nx_L(j),ny_L(j),nz_L(j)
    end do
!     write(*,*) 'j,nx_S(j),ny_S(j),nz_S(j)'
    do j=1,NBS_S
     nx_S(j)=nx(j,3)
     ny_S(j)=ny(j,3)
     nz_S(j)=nz(j,3)
!     write(*,*) j,nx_S(j),ny_S(j),nz_S(j)
    end do

      write(*,*) 'check translaten'
end subroutine translaten

!==============================================================================
subroutine translaten2(nx,ny,nz,n)
!==============================================================================
  implicit none
  integer i,j,k
  integer dnx,nsum
  integer,intent(out):: n
  integer,intent(in):: nx,ny,nz

  nsum = nx+ny+nz
  dnx = nsum-nx

  n = dnx*(dnx+1)/2 +dnx -ny +nsum*(nsum+1)*(nsum+2)/6

end subroutine translaten2

!==============================================================================
subroutine translaten3(nx,ny,nz,n)
!==============================================================================
  implicit none
  integer i,j,k
  integer dnx,nsum,ncnt,ntmp
  integer,intent(in):: n
  integer,intent(out):: nx,ny,nz
  
  nsum=0
  do
     ntmp = (nsum+1)*(nsum+2)*(nsum+3)/6
     if(ntmp.gt.n) exit
     nsum = nsum+1
  end do

  ncnt=0

  do i=1,nsum+1
     dnx=0
     do j=1,i
        do k=1,j
           ncnt = ncnt+1
           if(ncnt.gt.n) then
              nz = k-1
              exit
           end if
        end do
        if(ncnt.gt.n) exit
        dnx = dnx+1
     end do
     if(ncnt.gt.n) exit
  end do

  nx = nsum - dnx
  ny = dnx - nz


end subroutine translaten3

!==============================================================================
function renorm(nx,ny,nz)
! written by fukuda 2012/4/15
!
!According to DIRAC10/dft/basis_info.F90
!
!" a note on normalization:
!  AOs are normalized based on angular momentum and exponent
!  such that each shell has a common normalization
!  with this < AO | AO > =   1 for s, px, py, pz, dxy, dxz, dyz, fxyz
!                           3 for dxx, dyy, dzz, fxxy, ...
!                          15     fxxx, fyyy, fzzz, gxxxy, ...
!                         105     gxxxx, ...
!                         ...     ...
!                           9 for gxxyy, ...
!                         etc     ...
!  normalization is hidden in contraction_coef  "
!------------------------------------------------------------------------------
   implicit none
   double precision :: renorm
   double precision :: renorm1
   integer,intent(in) :: nx,ny,nz

   renorm = renorm1(nx)*renorm1(ny)*renorm1(nz)
end function renorm

!============================================================================
function renorm1(n)
!============================================================================
  implicit none
  double precision :: renorm1
  integer,intent(in) :: n
  integer :: i

  renorm1 = 1.d0

  if(n.ge.1) then
     do i=1,n
        renorm1 = renorm1*(2.d0*i - 1.d0)**(1.d0/2.d0)
     end do
  end if

end function renorm1

!!$!==============================================================================
!!$  function renorm(nx,ny,nz)
!!$! written by Takahashi 2010/12/21
!!$!
!!$!According to DIRAC10/dft/basis_info.F90
!!$!
!!$!" a note on normalization:
!!$!  AOs are normalized based on angular momentum and exponent
!!$!  such that each shell has a common normalization
!!$!  with this < AO | AO > =   1 for s, px, py, pz, dxy, dxz, dyz, fxyz
!!$!                           3 for dxx, dyy, dzz, fxxy, ...
!!$!                          15     fxxx, fyyy, fzzz, gxxxy, ...
!!$!                         105     gxxxx, ...
!!$!                         ...     ...
!!$!                           9 for gxxyy, ...
!!$!                         etc     ...
!!$!  normalization is hidden in contraction_coef  "
!!$!==============================================================================
!!$      implicit none
!!$      integer nx,ny,nz
!!$      double precision renorm,dorb,forb,gorb,horb,iorb
!!$      parameter ( dorb = dsqrt(3.d0) )
!!$      parameter ( forb = dsqrt(15.d0) )
!!$      parameter ( gorb = dsqrt(105.d0) )
!!$      parameter ( horb = dsqrt(945.d0) )
!!$      parameter ( iorb = dsqrt(10395.d0) )
!!$
!!$      renorm = 1.d0
!!$
!!$      if (nx .eq. 2) then
!!$      renorm = dorb*renorm
!!$      else if (nx .eq. 3) then
!!$      renorm = forb*renorm
!!$      else if (nx .eq. 4) then
!!$      renorm = gorb*renorm
!!$      else if (nx .eq. 5) then
!!$      renorm = horb*renorm
!!$      else if (nx .eq. 6) then
!!$      renorm = iorb*renorm
!!$      end if
!!$
!!$      if (ny .eq. 2) then
!!$      renorm = dorb*renorm
!!$      else if (ny .eq. 3) then
!!$      renorm = forb*renorm
!!$      else if (ny .eq. 4) then
!!$      renorm = gorb*renorm
!!$      else if (ny .eq. 5) then
!!$      renorm = horb*renorm
!!$      else if (ny .eq. 6) then
!!$      renorm = iorb*renorm
!!$      end if
!!$
!!$      if (nz .eq. 2) then
!!$      renorm = dorb*renorm
!!$      else if (nz .eq. 3) then
!!$      renorm = forb*renorm
!!$      else if (nz .eq. 4) then
!!$      renorm = gorb*renorm
!!$      else if (nz .eq. 5) then
!!$      renorm = horb*renorm
!!$      else if (nz .eq. 6) then
!!$      renorm = iorb*renorm
!!$      end if
!!$
!!$      return
!!$      end function renorm
!!$!==============================================================================
!!$!==============================================================================
!!$subroutine translaten(n)
!!$!==============================================================================
!!$  use DiracOutput
!!$      implicit none
!!$      integer j,k0,ntmp
!!$      integer,intent(in):: n(NBS0,4)
!!$      integer nx(NBS0,4),ny(NBS0,4),nz(NBS0,4)
!!$
!!$      do k0=1,3,2
!!$      do j=1,NBS0
!!$      if (n(j,k0) .eq. 0) then
!!$         nx(j,k0)=0
!!$         ny(j,k0)=0
!!$         nz(j,k0)=0
!!$      else if (n(j,k0) .eq. 1) then
!!$         nx(j,k0)=1
!!$         ny(j,k0)=0
!!$         nz(j,k0)=0
!!$      else if (n(j,k0) .eq. 2) then
!!$         nx(j,k0)=0
!!$         ny(j,k0)=1
!!$         nz(j,k0)=0
!!$      else if (n(j,k0) .eq. 3) then
!!$         nx(j,k0)=0
!!$         ny(j,k0)=0
!!$         nz(j,k0)=1
!!$      else if (n(j,k0) .eq. 4) then
!!$         nx(j,k0)=2
!!$         ny(j,k0)=0
!!$         nz(j,k0)=0
!!$      else if (n(j,k0) .eq. 5) then
!!$         nx(j,k0)=1
!!$         ny(j,k0)=1
!!$         nz(j,k0)=0
!!$      else if (n(j,k0) .eq. 6) then
!!$         nx(j,k0)=1
!!$         ny(j,k0)=0
!!$         nz(j,k0)=1
!!$      else if (n(j,k0) .eq. 7) then
!!$         nx(j,k0)=0 
!!$         ny(j,k0)=2
!!$         nz(j,k0)=0
!!$      else if (n(j,k0) .eq. 8) then
!!$         nx(j,k0)=0
!!$         ny(j,k0)=1
!!$         nz(j,k0)=1
!!$      else if (n(j,k0) .eq. 9) then
!!$         nx(j,k0)=0
!!$         ny(j,k0)=0
!!$         nz(j,k0)=2
!!$      else if (n(j,k0) .eq. 10) then
!!$         nx(j,k0)=3
!!$         ny(j,k0)=0
!!$         nz(j,k0)=0
!!$      else if (n(j,k0) .eq. 11) then
!!$         nx(j,k0)=2
!!$         ny(j,k0)=1
!!$         nz(j,k0)=0
!!$      else if (n(j,k0) .eq. 12) then
!!$         nx(j,k0)=2
!!$         ny(j,k0)=0
!!$         nz(j,k0)=1
!!$      else if (n(j,k0) .eq. 13) then
!!$         nx(j,k0)=1
!!$         ny(j,k0)=2
!!$         nz(j,k0)=0
!!$      else if (n(j,k0) .eq. 14) then
!!$         nx(j,k0)=1
!!$         ny(j,k0)=1
!!$         nz(j,k0)=1
!!$      else if (n(j,k0) .eq. 15) then
!!$         nx(j,k0)=1
!!$         ny(j,k0)=0
!!$         nz(j,k0)=2
!!$      else if (n(j,k0) .eq. 16) then
!!$         nx(j,k0)=0
!!$         ny(j,k0)=3
!!$         nz(j,k0)=0
!!$      else if (n(j,k0) .eq. 17) then
!!$         nx(j,k0)=0
!!$         ny(j,k0)=2
!!$         nz(j,k0)=1
!!$      else if (n(j,k0) .eq. 18) then
!!$         nx(j,k0)=0
!!$         ny(j,k0)=1
!!$         nz(j,k0)=2
!!$      else if (n(j,k0) .eq. 19) then
!!$         nx(j,k0)=0
!!$         ny(j,k0)=0
!!$         nz(j,k0)=3
!!$      else if (n(j,k0) .eq. 20) then
!!$         nx(j,k0)=4
!!$         ny(j,k0)=0
!!$         nz(j,k0)=0
!!$      else if (n(j,k0) .eq. 21) then
!!$         nx(j,k0)=3
!!$         ny(j,k0)=1
!!$         nz(j,k0)=0
!!$      else if (n(j,k0) .eq. 22) then
!!$         nx(j,k0)=3
!!$         ny(j,k0)=0
!!$         nz(j,k0)=1
!!$      else if (n(j,k0) .eq. 23) then
!!$         nx(j,k0)=2
!!$         ny(j,k0)=2
!!$         nz(j,k0)=0
!!$      else if (n(j,k0) .eq. 24) then
!!$         nx(j,k0)=2
!!$         ny(j,k0)=1
!!$         nz(j,k0)=1
!!$      else if (n(j,k0) .eq. 25) then
!!$         nx(j,k0)=2
!!$         ny(j,k0)=0
!!$         nz(j,k0)=2
!!$      else if (n(j,k0) .eq. 26) then
!!$         nx(j,k0)=1
!!$         ny(j,k0)=3
!!$         nz(j,k0)=0
!!$      else if (n(j,k0) .eq. 27) then
!!$         nx(j,k0)=1
!!$         ny(j,k0)=2
!!$         nz(j,k0)=1
!!$      else if (n(j,k0) .eq. 28) then
!!$         nx(j,k0)=1
!!$         ny(j,k0)=1
!!$         nz(j,k0)=2
!!$      else if (n(j,k0) .eq. 29) then
!!$         nx(j,k0)=1
!!$         ny(j,k0)=0
!!$         nz(j,k0)=3
!!$      else if (n(j,k0) .eq. 30) then
!!$         nx(j,k0)=0
!!$         ny(j,k0)=4
!!$         nz(j,k0)=0
!!$      else if (n(j,k0) .eq. 31) then
!!$         nx(j,k0)=0
!!$         ny(j,k0)=3
!!$         nz(j,k0)=1
!!$      else if (n(j,k0) .eq. 32) then
!!$         nx(j,k0)=0
!!$         ny(j,k0)=2
!!$         nz(j,k0)=2
!!$      else if (n(j,k0) .eq. 33) then
!!$         nx(j,k0)=0
!!$         ny(j,k0)=1
!!$         nz(j,k0)=3
!!$      else if (n(j,k0) .eq. 34) then
!!$         nx(j,k0)=0
!!$         ny(j,k0)=0
!!$         nz(j,k0)=4
!!$      else if (n(j,k0) .eq. 35) then
!!$         nx(j,k0)=5
!!$         ny(j,k0)=0
!!$         nz(j,k0)=0
!!$      else if (n(j,k0) .eq. 36) then
!!$         nx(j,k0)=4
!!$         ny(j,k0)=1
!!$         nz(j,k0)=0
!!$      else if (n(j,k0) .eq. 37) then
!!$         nx(j,k0)=4
!!$         ny(j,k0)=0
!!$         nz(j,k0)=1
!!$      else if (n(j,k0) .eq. 38) then
!!$         nx(j,k0)=3
!!$         ny(j,k0)=2
!!$         nz(j,k0)=0
!!$      else if (n(j,k0) .eq. 39) then
!!$         nx(j,k0)=3
!!$         ny(j,k0)=1
!!$         nz(j,k0)=1
!!$      else if (n(j,k0) .eq. 40) then
!!$         nx(j,k0)=3
!!$         ny(j,k0)=0
!!$         nz(j,k0)=2
!!$      else if (n(j,k0) .eq. 41) then
!!$         nx(j,k0)=2
!!$         ny(j,k0)=3
!!$         nz(j,k0)=0
!!$      else if (n(j,k0) .eq. 42) then
!!$         nx(j,k0)=2
!!$         ny(j,k0)=2
!!$         nz(j,k0)=1
!!$      else if (n(j,k0) .eq. 43) then
!!$         nx(j,k0)=2
!!$         ny(j,k0)=1
!!$         nz(j,k0)=2
!!$      else if (n(j,k0) .eq. 44) then
!!$         nx(j,k0)=2
!!$         ny(j,k0)=0
!!$         nz(j,k0)=3
!!$      else if (n(j,k0) .eq. 45) then
!!$         nx(j,k0)=1
!!$         ny(j,k0)=4
!!$         nz(j,k0)=0
!!$      else if (n(j,k0) .eq. 46) then
!!$         nx(j,k0)=1
!!$         ny(j,k0)=3
!!$         nz(j,k0)=1
!!$      else if (n(j,k0) .eq. 47) then
!!$         nx(j,k0)=1
!!$         ny(j,k0)=2
!!$         nz(j,k0)=2
!!$      else if (n(j,k0) .eq. 48) then
!!$         nx(j,k0)=1
!!$         ny(j,k0)=1
!!$         nz(j,k0)=3
!!$      else if (n(j,k0) .eq. 49) then
!!$         nx(j,k0)=1
!!$         ny(j,k0)=0
!!$         nz(j,k0)=4
!!$      else if (n(j,k0) .eq. 50) then
!!$         nx(j,k0)=0
!!$         ny(j,k0)=5
!!$         nz(j,k0)=0
!!$      else if (n(j,k0) .eq. 51) then
!!$         nx(j,k0)=0
!!$         ny(j,k0)=4
!!$         nz(j,k0)=1
!!$      else if (n(j,k0) .eq. 52) then
!!$         nx(j,k0)=0
!!$         ny(j,k0)=3
!!$         nz(j,k0)=2
!!$      else if (n(j,k0) .eq. 53) then
!!$         nx(j,k0)=0
!!$         ny(j,k0)=2
!!$         nz(j,k0)=3
!!$      else if (n(j,k0) .eq. 54) then
!!$         nx(j,k0)=0
!!$         ny(j,k0)=1
!!$         nz(j,k0)=4
!!$      else if (n(j,k0) .eq. 55) then
!!$         nx(j,k0)=0
!!$         ny(j,k0)=0
!!$         nz(j,k0)=5
!!$      else if (n(j,k0) .eq.  56) then
!!$         nx(j,k0) =  6
!!$         ny(j,k0) =  0
!!$         nz(j,k0) =  0
!!$      else if (n(j,k0) .eq.  57) then
!!$         nx(j,k0) =  5
!!$         ny(j,k0) =  1
!!$         nz(j,k0) =  0
!!$      else if (n(j,k0) .eq.  58) then
!!$         nx(j,k0) =  5
!!$         ny(j,k0) =  0
!!$         nz(j,k0) =  1
!!$      else if (n(j,k0) .eq.  59) then
!!$         nx(j,k0) =  4
!!$         ny(j,k0) =  2
!!$         nz(j,k0) =  0
!!$      else if (n(j,k0) .eq.  60) then
!!$         nx(j,k0) =  4
!!$         ny(j,k0) =  1
!!$         nz(j,k0) =  1
!!$      else if (n(j,k0) .eq.  61) then
!!$         nx(j,k0) =  4
!!$         ny(j,k0) =  0
!!$         nz(j,k0) =  2
!!$      else if (n(j,k0) .eq.  62) then
!!$         nx(j,k0) =  3
!!$         ny(j,k0) =  3
!!$         nz(j,k0) =  0
!!$      else if (n(j,k0) .eq.  63) then
!!$         nx(j,k0) =  3
!!$         ny(j,k0) =  2
!!$         nz(j,k0) =  1
!!$      else if (n(j,k0) .eq.  64) then
!!$         nx(j,k0) =  3
!!$         ny(j,k0) =  1
!!$         nz(j,k0) =  2
!!$      else if (n(j,k0) .eq.  65) then
!!$         nx(j,k0) =  3
!!$         ny(j,k0) =  0
!!$         nz(j,k0) =  3
!!$      else if (n(j,k0) .eq.  66) then
!!$         nx(j,k0) =  2
!!$         ny(j,k0) =  4
!!$         nz(j,k0) =  0
!!$      else if (n(j,k0) .eq.  67) then
!!$         nx(j,k0) =  2
!!$         ny(j,k0) =  3
!!$         nz(j,k0) =  1
!!$      else if (n(j,k0) .eq.  68) then
!!$         nx(j,k0) =  2
!!$         ny(j,k0) =  2
!!$         nz(j,k0) =  2
!!$      else if (n(j,k0) .eq.  69) then
!!$         nx(j,k0) =  2
!!$         ny(j,k0) =  1
!!$         nz(j,k0) =  3
!!$      else if (n(j,k0) .eq.  70) then
!!$         nx(j,k0) =  2
!!$         ny(j,k0) =  0
!!$         nz(j,k0) =  4
!!$      else if (n(j,k0) .eq.  71) then
!!$         nx(j,k0) =  1
!!$         ny(j,k0) =  5
!!$         nz(j,k0) =  0
!!$      else if (n(j,k0) .eq.  72) then
!!$         nx(j,k0) =  1
!!$         ny(j,k0) =  4
!!$         nz(j,k0) =  1
!!$      else if (n(j,k0) .eq.  73) then
!!$         nx(j,k0) =  1
!!$         ny(j,k0) =  3
!!$         nz(j,k0) =  2
!!$      else if (n(j,k0) .eq.  74) then
!!$         nx(j,k0) =  1
!!$         ny(j,k0) =  2
!!$         nz(j,k0) =  3
!!$      else if (n(j,k0) .eq.  75) then
!!$         nx(j,k0) =  1
!!$         ny(j,k0) =  1
!!$         nz(j,k0) =  4
!!$      else if (n(j,k0) .eq.  76) then
!!$         nx(j,k0) =  1
!!$         ny(j,k0) =  0
!!$         nz(j,k0) =  5
!!$      else if (n(j,k0) .eq.  77) then
!!$         nx(j,k0) =  0
!!$         ny(j,k0) =  6
!!$         nz(j,k0) =  0
!!$      else if (n(j,k0) .eq.  78) then
!!$         nx(j,k0) =  0
!!$         ny(j,k0) =  5
!!$         nz(j,k0) =  1
!!$      else if (n(j,k0) .eq.  79) then
!!$         nx(j,k0) =  0
!!$         ny(j,k0) =  4
!!$         nz(j,k0) =  2
!!$      else if (n(j,k0) .eq.  80) then
!!$         nx(j,k0) =  0
!!$         ny(j,k0) =  3
!!$         nz(j,k0) =  3
!!$      else if (n(j,k0) .eq.  81) then
!!$         nx(j,k0) =  0
!!$         ny(j,k0) =  2
!!$         nz(j,k0) =  4
!!$      else if (n(j,k0) .eq.  82) then
!!$         nx(j,k0) =  0
!!$         ny(j,k0) =  1
!!$         nz(j,k0) =  5
!!$      else if (n(j,k0) .eq.  83) then
!!$         nx(j,k0) =  0
!!$         ny(j,k0) =  0
!!$         nz(j,k0) =  6
!!$      end if
!!$         nx(j,k0+1)=nx(j,k0)
!!$         ny(j,k0+1)=ny(j,k0)
!!$         nz(j,k0+1)=nz(j,k0)
!!$      end do
!!$      end do
!!$
!!$!!$      do k0=1,3,2
!!$!!$        do j=1,NBS0
!!$!!$          call translaten2(nx(j,k0),ny(j,k0),nz(j,k0),ntmp)
!!$!!$          write(*,'(a,6I3)')'j,k0,n,nx,ny,nz',j,k0,ntmp,nx(j,k0),ny(j,k0),nz(j,k0)
!!$!!$        end do
!!$!!$      end do
!!$!!$      stop
!!$
!!$!     write(*,*) 'j,nx_L(j),ny_L(j),nz_L(j)'
!!$    do j=1,NBS_L
!!$     nx_L(j)=nx(j,1)
!!$     ny_L(j)=ny(j,1)
!!$     nz_L(j)=nz(j,1)
!!$!     write(*,*) j,nx_L(j),ny_L(j),nz_L(j)
!!$    end do
!!$!     write(*,*) 'j,nx_S(j),ny_S(j),nz_S(j)'
!!$    do j=1,NBS_S
!!$     nx_S(j)=nx(j,3)
!!$     ny_S(j)=ny(j,3)
!!$     nz_S(j)=nz(j,3)
!!$!     write(*,*) j,nx_S(j),ny_S(j),nz_S(j)
!!$    end do
!!$
!!$      write(*,*) 'check translaten'
!!$end subroutine translaten
!!$
!==============================================================================
