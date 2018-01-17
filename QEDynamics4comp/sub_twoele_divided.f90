!======================================================
!  subroutine calctwoele_divided1
!  subroutine calctwoele_divided2  <- calculate twoele
!  subroutine connect_twoele       <- make file twoele.dat
!  subroutine convert_twoele       <- convert twoele.dat to twoele.bin
!  subroutine setQmat_twoele_bin(twoele_Qmat)
!   ( NM | PQ )
!======================================================
!======================================================
subroutine calctwoele_divided1
!======================================================
  use DiracOutput
  implicit none

  call read_DiracOutput  ! read and set global variables -> defined in sub_readDiracOutput.f90

  call calctwoele_divided2
!  call connect_twoele
!  call convert_twoele

end subroutine calctwoele_divided1
!======================================================
subroutine calctwoele_divided2
! ( NM | PQ )
!======================================================
  use DiracOutput
  implicit none
  
!  complex(kind=8),intent(out) :: twoele_Qmat(4*NBS,4*NBS,4*NBS,4*NBS)
  
  complex(kind=8) :: inttwoele_mat
  integer :: nn,mm,pp,qq  ! index for 1~4*NBS
  integer :: n,m,p,q ! index for 1~2*NBS (+ or - specified)
  character(LEN=1) :: a,b,c,d

  character(LEN=80) :: temp, filename
  real(kind=8) :: int_real,int_complex

  character(LEN=10) :: str1,str2,str3,str4


  call getarg(1,str1)
  call getarg(2,str2)
  call getarg(3,str3)
  call getarg(4,str4)

  read(str1,'(i4.4)') nn  ! 
  read(str2,'(i4.4)') mm  ! 
  read(str3,'(i4.4)') pp  ! 
  read(str4,'(i4.4)') qq  ! 
     write(filename,'(a,i4.4,i4.4,i4.4,i4.4,a)') 'twoele_d',nn,mm,pp,qq,'.dat'
!     write(filename,'(a,i4.4,i4.4,i4.4,a)') 'twoele_d',nn,mm,pp,'.dat'
     open(unit=250,file=filename)

!           do pp=1,4*NBS
!              do qq=1,4*NBS
                 
                 call index_from_Qmat(nn,n,a)
                 call index_from_Qmat(mm,m,b)
                 call index_from_Qmat(pp,p,c)
                 call index_from_Qmat(qq,q,d)
                 
!                 write(12,"(4a3,4i6,2es16.6)") a,b,c,d,n,m,p,q,inttwoele_mat(n,a,m,b,p,c,q,d) 
                 write(250,"(2es16.6)") inttwoele_mat(n,a,m,b,p,c,q,d)  ! for storage
!                 twoele_Qmat(nn,mm,pp,qq) = inttwoele_mat(n,a,m,b,p,c,q,d) 
!              end do
!           end do
     write(*,*) " Integrals twoele done."

     close(unit=250)

  return
end subroutine calctwoele_divided2

!======================================================
subroutine connect_twoele
!======================================================
  use DiracOutput
  implicit none

  complex(kind=8) :: twoele_multi
  integer :: nn,mm,pp,qq  ! index for 1~4*NBS
  character(LEN=80) :: filename,filename0

     
    write(filename0,'(a)') 'twoeleall.dat'
    open(unit=254,file=filename0)

     do nn=1,4*NBS
        do mm=1,4*NBS
           do pp=1,4*NBS
              do qq=1,4*NBS
                 write(filename,'(a,i4.4,i4.4,i4.4,i4.4,a)') 'twoele_d',nn,mm,pp,qq,'.dat'
                 open(unit=253,file=filename)

                 read(253,"(2es16.6)") twoele_multi  ! for storage
                 write(254,"(2es16.6)") twoele_multi  ! for storage
              end do
           end do

          close(unit=253)
        end do
     end do

    close(unit=254)

  return
end subroutine connect_twoele
!======================================================
subroutine convert_twoele
!======================================================
  use DiracOutput
  use IntegralStorage
  implicit none

  complex(kind=8) :: inttwoele_mat
  integer :: nn,mm,pp,qq  ! index for 1~4*NBS
  integer :: n,m,p,q ! index for 1~2*NBS (+ or - specified)
  character(LEN=1) :: a,b,c,d

  character(LEN=80) :: temp
  real(kind=8) :: int_real,int_complex

  write(*,*) "Convert twoele.dat to twoele.bin"
  open(unit=200,file=file_twoele,status='unknown',form='formatted')
  open(unit=300,file=file_twoele_bin,status='unknown',form='unformatted')

  do nn=1,4*NBS
     do mm=1,4*NBS
        do pp=1,4*NBS
           do qq=1,4*NBS

              ! use same format as storage (see below)
              read(200,*) int_real,int_complex            
              write(300) int_real,int_complex            

           end do
        end do
     end do
  end do

  close(unit=200)
  close(unit=300)

  return
end subroutine convert_twoele
!======================================================
subroutine setQmat_twoele_bin(twoele_Qmat)
! ( NM | PQ )
!======================================================
  use DiracOutput
  use IntegralStorage
  implicit none
  
  complex(kind=8),intent(out) :: twoele_Qmat(4*NBS,4*NBS,4*NBS,4*NBS)
  
  complex(kind=8) :: inttwoele_mat
  integer :: nn,mm,pp,qq  ! index for 1~4*NBS
  integer :: n,m,p,q ! index for 1~2*NBS (+ or - specified)
  character(LEN=1) :: a,b,c,d

  character(LEN=80) :: temp
  real(kind=8) :: int_real,int_complex
  character(LEN=200) :: file_twoele_bin2

  if(there_is_twoele) then
     write(*,*) " There is twoele file. Read integrals from ",file_twoele_bin
     file_twoele_bin2 = trim(FILESFOLDER)//"/"//file_twoele_bin
     open(unit=300,file=file_twoele_bin2,status='unknown',form='unformatted')

     do nn=1,4*NBS
        do mm=1,4*NBS
           do pp=1,4*NBS
              do qq=1,4*NBS

                 ! use same format as storage (see below)
                 read(300) int_real,int_complex            
                 twoele_Qmat(nn,mm,pp,qq) = cmplx(int_real,int_complex)
!                 write(206,*)int_real, int_complex

              end do
           end do
        end do
     end do
     write(*,*) " Reading done."
     
  else ! if there isn't precomputed integrals, compute integrals.
     write(*,*) " There is not a twoele file. Compute integrals."
     do nn=1,4*NBS
        do mm=1,4*NBS
           do pp=1,4*NBS
              do qq=1,4*NBS
                 
                 call index_from_Qmat(nn,n,a)
                 call index_from_Qmat(mm,m,b)
                 call index_from_Qmat(pp,p,c)
                 call index_from_Qmat(qq,q,d)
                 
!                 write(12,"(4a3,4i6,2es16.6)") a,b,c,d,n,m,p,q,inttwoele_mat(n,a,m,b,p,c,q,d) 
                 write(12,"(2es16.6)") inttwoele_mat(n,a,m,b,p,c,q,d)  ! for storage
                 twoele_Qmat(nn,mm,pp,qq) = inttwoele_mat(n,a,m,b,p,c,q,d) 
              end do
           end do
        end do
     end do
     write(*,*) " Integrals twoele done."

  end if

  close(unit=300)

  return
end subroutine setQmat_twoele_bin

