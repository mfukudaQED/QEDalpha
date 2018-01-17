! Last Change:18-Jun-2012.
!==============================================================================
! subroutine readBasis version 2012/4/16 by Fukuda
! subroutine countBasis version 2011/5/22 by Fukuda
! subroutine ShellNumber 2012/04/16 by Fukuda
!==============================================================================
!==============================================================================
subroutine readBasis(bnum_L,bnum_S,bnum_LT,bnum_ST,bnum_LTA,bnum_STA,cpg_L,cpg_S,n,cont)
!  read basis.txt
!  set cpg_L,cpg_S,aa_L,aa_S
!==============================================================================
  use DiracOutput

  implicit none
 
  character(len=4) Atom
  character(len=5) pg_orbital
  integer atomicnum
  integer ic,itmp,i,j,k,l,inat
  integer numShell, sumShell, numNs, numN2
  integer,intent(out):: n(NBS0,4)  ! number in 'subroutine translaten'
  integer,intent(in):: bnum_L(100),bnum_S(100)
  integer,intent(in):: bnum_LT,bnum_ST
  integer,intent(in):: bnum_LTA(100),bnum_STA(100)
  double precision tmpx,tmpy,tmpz,tmpa,tmpb
  double precision,intent(out):: cpg_L(NBS_L),cpg_S(NBS_S)
  double precision nrmca(100), nrmcb(100), nrmcb2(100), resfac
  integer,intent(out):: cont(NBS00,NAT)
  integer tmpcont
  integer mcontL(NAT)
  integer mcontS(NAT)

  write(*,*)'start reading basis.txt'

  !------------Large component--------------------------
   ic=0 !count NBS_L
   itmp=0
   mcontL(:)=0 ! counter of (inat)th atom's contracted primitive gaussian number
   cont(:,:)=0

 do inat=1,NAT
   read(11,'(a4)') Atom

  if(Atom.ne.'Atom') then
   write(*,*)'error Basis.txt is not correct L1'
   stop
  end if

   read(11,'(a4,6x,i3,x,3f20.10)') Atom, atomicnum, tmpx, tmpy, tmpz
    xc(inat)=tmpx ! position of atom
    yc(inat)=tmpy
    zc(inat)=tmpz
    cn(inat)=atomicnum
   write(*,*)'xc(',inat,'),yc(',inat,'),zc(',inat,')',xc(inat),yc(inat),zc(inat)


   do j=1,bnum_L(inat)
     read(11,'(a5,i5)') pg_orbital, tmpcont
!     write(*,*) 'pg_orbital, tmpcont',pg_orbital, tmpcont
     call ShellNumber(pg_orbital, numShell, sumShell, numNs, numN2)
     do k=1,numShell
        mcontL(inat) = mcontL(inat) + 1
        cont(mcontL(inat),inat) = tmpcont
!        write(*,*)'contb',mcontL(inat),inat,cont(mcontL(inat),inat)
     end do

     if (tmpcont == 1) then
        do l=1,tmpcont
          read(11,'(5x,e15.9,e20.10)') tmpa, tmpb
            do k=1,numShell
               ic=ic+1
               aa_L(ic) = tmpa
               cpg_L(ic) = tmpb
               n(ic,1) = k + sumShell -1
            end do
        end do
     
     else !cont /= 1

        do l=1,tmpcont
          read(11,'(5x,e15.9,e20.10)') tmpa, tmpb
!          write(*,*) 'tmpa, tmpb',tmpa, tmpb
              nrmcb2(l) = tmpb
          call normals(tmpa,tmpb,numNs)
!          write(*,*) 'normals tmpa, tmpb',tmpa, tmpb
              nrmca(l) = tmpa
              nrmcb(l) = tmpb
        end do

        call normal2(nrmca,nrmcb,numN2,tmpcont,resfac)
        do k=1,numShell
           do l=1,tmpcont
              ic=ic+1
              aa_L(ic) = nrmca(l)
              cpg_L(ic) = nrmcb2(l) * resfac
              n(ic,1) = k + sumShell -1
!	      write(*,*)'cpg_L(',ic,')=',cpg_L(ic)
!	      write(*,*)'aa_L(',ic,')=',aa_L(ic)
           end do
        end do

     end if
   end do
!============================================================================

   do j=itmp+1,ic
    xx_L(j)=tmpx    ! position of primitive gaussian
    yy_L(j)=tmpy    ! position of primitive gaussian
    zz_L(j)=tmpz    ! position of primitive gaussian
!    write(*,*)'xx_L(j),yy_L(j),zz_L(j)',xx_L(j),yy_L(j),zz_L(j)
   end do

   itmp=ic
 end do
   !-----------error check large basis----------------
    read(11,'(a4)') Atom
    if(Atom.eq.'DIV') then
      write(*,*)'check pg Large'
    end if
 
    if(ic.ne.NBS_L) then
      write(*,*)'error basis large'
      write(*,*)'ic=',ic,'NBS_L=',NBS_L
      stop
    end if

    do inat=1,NAT
       if(mcontL(inat).ne.bnum_LTA(inat)) then
          write(*,*)'error mcontL',mcontL(inat),bnum_LTA(inat)
          stop
       end if
    end do
 
    do j=1,NBS_L
     n(j,2) = n(j,1)
    end do
    write(*,*)'ic=',ic,'NBS_L=',NBS_L
 
!     write(*,*)'aa_L, cpg_L'
!    do k=1,NBS_L
!     write(*,*)aa_L(k),cpg_L(k)
!    end do

   !-------------Small component----------------------------------
   ic=0 !count NBS_S
   itmp=0
   do inat=1,NAT
     mcontS(inat)=0 ! counter of (inat)th atom's contracted primitive gaussian number
   end do

 do inat=1,NAT
  read(11,'(a4)') Atom

  if(Atom.ne.'Atom') then
   write(*,*)'error Basis.txt is not correct S2'
   stop
  end if

   read(11,'(a4,6x,i3,x,3f20.10)') Atom, atomicnum, tmpx, tmpy, tmpz

   do j=1,bnum_S(inat)
     read(11,'(a5,i5)') pg_orbital, tmpcont
     call ShellNumber(pg_orbital, numShell, sumShell, numNs, numN2)
     do k=1,numShell
        mcontS(inat) = mcontS(inat) + 1
        cont(mcontS(inat)+bnum_LTA(inat),inat) = tmpcont
     end do

     if (tmpcont == 1) then
        do l=1,tmpcont
          read(11,'(5x,e15.9,e20.10)') tmpa, tmpb
            do k=1,numShell
               ic=ic+1
               aa_S(ic) = tmpa
               cpg_S(ic) = tmpb
               n(ic,3) = k + sumShell -1
!	      write(*,*)'cpg_S(',ic,')=',cpg_S(ic)
!	      write(*,*)'aa_S(',ic,')=',aa_S(ic)
            end do
        end do
     
     else !cont /= 1

        do l=1,tmpcont
          read(11,'(5x,e15.9,e20.10)') tmpa, tmpb
              nrmcb2(l) = tmpb
          call normals(tmpa,tmpb,numNs)
              nrmca(l) = tmpa
              nrmcb(l) = tmpb
        end do
   
        call normal2(nrmca,nrmcb,numN2,tmpcont,resfac)
        do k=1,numShell
           do l=1,tmpcont
              ic=ic+1
              aa_S(ic) = nrmca(l)
              cpg_S(ic) = nrmcb2(l) * resfac
              n(ic,3) = k + sumShell -1
           end do
        end do

     end if
   end do
!============================================================================

   do j=itmp+1,ic
    xx_S(j)=tmpx    ! position of primitive gaussian
    yy_S(j)=tmpy    ! position of primitive gaussian
    zz_S(j)=tmpz    ! position of primitive gaussian
   end do

   itmp=ic

 end do

  !-----------error check small basis--------------
   read(11,'(a4)') Atom
    if(Atom.eq.'DIV') then
     write(*,*)'check pg Small'
    end if

    if(ic.ne.NBS_S) then
     write(*,*)'error basis small'
     write(*,*)'ic=',ic,'NBS_S=',NBS_S
     stop
    end if

   write(*,*)'ic=',ic,'NBS_S=',NBS_S

   do j=1,NBS_S
    n(j,4) = n(j,3)
   end do

   if(NBS_L.gt.NBS_S) then
      do i=NBS_S+1,NBS0
         do k=3,4
            n(i,k)  = -1.D0
         end do
      end do
   else
      do i=NBS_L+1,NBS0
         do k=1,2
            n(i,k)  = -1.D0
         end do
      end do
   end if

   call translaten(n) ! set nx,ny,nz

!    write(*,*)'aa_S, cpg_S'
!   do k=1,NBS_S
!    write(*,*)aa_S(k),cpg_S(k)
!   end do
!-------------error check basis.txt--------------------
  read(11,'(a4)') Atom
   if(Atom.eq.'End ') then
    write(*,*)'check readBasis'
   end if

end subroutine readBasis
 
!==============================================================================
subroutine countBasis(bnum_L,bnum_S,bnum_LT,bnum_ST,bnum_LTA,bnum_STA)
! set bnum_L,bnum_S,bnum_LT,bnum_ST,NBS_L,NBS_S,NBS0,NBS00
!==============================================================================
  use DiracOutput

  implicit none

  integer,intent(out):: bnum_L(100),bnum_S(100)
  integer i,ic,jc,countNBS_L,countNBS_S !counter
  integer,intent(out):: bnum_LT,bnum_ST
  integer,intent(out):: bnum_LTA(100),bnum_STA(100)
  integer numShell, sumShell, numNs, numN2
  character(len=5) MOJI

  write(*,*)'start countBasis'
   do i=1,100 !This calculation assumed less than 100 atoms
    bnum_L(i)=0
    bnum_S(i)=0
    bnum_LTA(i)=0
    bnum_STA(i)=0
   end do

   bnum_LT = 0
   bnum_ST = 0

!----------count the number of large basis-----------------------------------
   ic=0 !atom counter
   countNBS_L = 0
    read(11,'(a5)') MOJI
    if (MOJI.ne."Atom ") then
     write(*,*)'error read basis.txt in countprimL'
    end if

   do
      ic = ic+1
      read(11,'(a5)') MOJI !kind of atom
      do ! loopA
        read(11,'(a5)') MOJI
        if(MOJI.eq.'Atom ') then
           exit
        else if(MOJI.eq.'DIV  ') then 
           exit
        else 
          do i=1,5
             if(MOJI(i:i)>='A' .and. MOJI(i:i)<='Z') then
                 bnum_L(ic) = bnum_L(ic) + 1
                 call ShellNumber(MOJI, numShell, sumShell, numNs, numN2)
                 bnum_LT = bnum_LT + numShell
                 bnum_LTA(ic) = bnum_LTA(ic) + numShell
                 countNBS_L = countNBS_L - numShell
             end if
          end do
          countNBS_L = countNBS_L + numShell 
        end if
      end do !back to loopA
    if(MOJI.eq.'DIV  ') exit
   end do

   NAT = ic !the number of atoms
   NBS_L = countNBS_L

!----------count the number of small basis-----------------------------------
   jc=0 !atom counter
   countNBS_S = 0
    read(11,'(a5)') MOJI
    if (MOJI.ne."Atom ") then
     write(*,*)'error read basis.txt in countprimS'
    end if

   do
      jc = jc+1
      read(11,'(a5)') MOJI !kind of atom
     do !loopB
        read(11,'(a5)') MOJI
        if(MOJI.eq.'Atom ') then
           exit
        else if(MOJI.eq.'DIV  ') then
           exit
        else
           do i=1,5
              if(MOJI(i:i)>='A' .and. MOJI(i:i)<='Z') then
                 bnum_S(jc) = bnum_S(jc) + 1
                 call ShellNumber(MOJI, numShell, sumShell, numNs, numN2)
                 bnum_ST = bnum_ST + numShell
                 bnum_STA(jc) = bnum_STA(jc) + numShell
                 countNBS_S = countNBS_S - numShell
              end if
           end do
           countNBS_S = countNBS_S + numShell
        end if
      end do !back to loopB
    if(MOJI.eq.'DIV  ') exit
   end do

   if(jc.ne.NAT) then
     write(*,*)'error read basis.txt in countprimNAT'
   end if

   NBS_S = countNBS_S
   NBS00 = NBS_L + NBS_S
   NBS0 = max(NBS_L,NBS_S) 
 
   read(11,'(a5)') MOJI
    if(MOJI.eq.'End  ') then
     write(*,*)'end countBasis'
    end if
 
end subroutine countBasis

!============================================================================== 
subroutine ShellNumber(pg_orbital, numShell, sumShell, numNs, numN2)
! numbers used in this subroutine are from relativistic MRDFT(CI) program
!   numShell : 1 3 6 10 15 21 28
!              S P D  F  G  H  I
!   sumShell : 0 1 4 10 20 35 56
!              S P D  F  G  H  I
!=============================================================================
  implicit none

  integer num, numShell, sumShell, numNs, numN2
  character(LEN=5) pg_orbital

 
   if (pg_orbital.eq.'    S') then
      num = 0
   else if (pg_orbital.eq.'    P') then
      num = 1
   else if (pg_orbital.eq.'    D') then
      num = 2
   else if (pg_orbital.eq.'    F') then
      num = 3
   else if (pg_orbital.eq.'    G') then
      num = 4
   else if (pg_orbital.eq.'    H') then
      num = 5
   else if (pg_orbital.eq.'    I') then
      num = 6
   else
      write(*,*)'error pg_orbital in subroutine ShellNumber'
   end if

   numShell = (num+1)*(num+2)/2  ! 1 3 6 10 15 21 28
                             ! S P D  F  G  H  I
   sumShell = num*(num+1)*(num+2)/6  ! 0 1 4 10 20 35 56
                                     ! S P D  F  G  H  I
   numNs = sumShell
   numN2 = num

end subroutine ShellNumber
!============================================================================
