!============================================================================
function fs(aa,a,r2)
!============================================================================
  implicit none
  complex(kind(0d0)) :: fs
  complex(kind(0d0)),intent(in) :: a
  double precision,intent(in) :: aa
  double precision,intent(in) :: r2

      fs = a * dexp(-aa*r2)

end function fs
!============================================================================
function fp(aa,z,a,r2)
!============================================================================
  implicit none
  complex(kind(0d0)) :: fp
  complex(kind(0d0)),intent(in) :: a
  double precision,intent(in) :: aa
  double precision,intent(in) :: r2
  double precision,intent(in) :: z

      fp = a * z*dexp(-aa*r2)

end function fp
!============================================================================
function fd0(aa,z,a,r2)
!============================================================================
  implicit none
  complex(kind(0d0)) :: fd0
  complex(kind(0d0)),intent(in) :: a
  double precision,intent(in) :: aa
  double precision,intent(in) :: r2
  double precision,intent(in) :: z

      fd0 = a * z*z*dexp(-aa*r2)

end function fd0
!============================================================================
function fd2(aa,x,y,a,r2)
!============================================================================
  implicit none
  complex(kind(0d0)) :: fd2
  complex(kind(0d0)),intent(in) :: a
  double precision,intent(in) :: aa
  double precision,intent(in) :: r2
  double precision,intent(in) :: x,y

      fd2 = a * x*y*dexp(-aa*r2)

endfunction fd2
!============================================================================
function ff0(aa,z,a,r2)
!============================================================================
  implicit none
  complex(kind(0d0)) :: ff0
  complex(kind(0d0)),intent(in) :: a
  double precision,intent(in) :: aa
  double precision,intent(in) :: r2
  double precision,intent(in) :: z

      ff0 = a * z*z*z*dexp(-aa*r2)

end function ff0
!============================================================================
function ff6(aa,z,u,a,r2)
!============================================================================
  implicit none
  complex(kind(0d0)) :: ff6
  complex(kind(0d0)),intent(in) :: a
  double precision,intent(in) :: aa
  double precision,intent(in) :: r2
  double precision,intent(in) :: z,u

      ff6 = a * u*z*z*dexp(-aa*r2)

end function ff6
!============================================================================
function ff4(aa,x,y,z,a,r2)
!============================================================================
  implicit none
  complex(kind(0d0)) :: ff4
  complex(kind(0d0)),intent(in) :: a
  double precision,intent(in) :: aa
  double precision,intent(in) :: r2
  double precision,intent(in) :: x,y,z

      ff4 = a * x*y*z*dexp(-aa*r2)

end function ff4
!============================================================================
function fg0(aa,z,a,r2)
!============================================================================
  implicit none
  complex(kind(0d0)) :: fg0
  complex(kind(0d0)),intent(in) :: a
  double precision,intent(in) :: aa
  double precision,intent(in) :: r2
  double precision,intent(in) :: z

      fg0 = a * z*z*z*z*dexp(-aa*r2)

end function fg0
!============================================================================
function fg1(aa,z,u,a,r2)
!============================================================================
  implicit none
  complex(kind(0d0)) :: fg1
  complex(kind(0d0)),intent(in) :: a
  double precision,intent(in) :: aa
  double precision,intent(in) :: r2
  double precision,intent(in) :: z,u

      fg1 = a * z*z*z*u*dexp(-aa*r2)

end function fg1
!============================================================================
function fg2(aa,z,u,a,r2)
!============================================================================
  implicit none
  complex(kind(0d0)) :: fg2
  complex(kind(0d0)),intent(in) :: a
  double precision,intent(in) :: aa
  double precision,intent(in) :: r2
  double precision,intent(in) :: z,u

      fg2 = a * z*z*u*u*dexp(-aa*r2)

end function fg2
!============================================================================
function fg3(aa,z,u,v,a,r2)
!============================================================================
  implicit none
  complex(kind(0d0)) :: fg3
  complex(kind(0d0)),intent(in) :: a
  double precision,intent(in) :: aa
  double precision,intent(in) :: r2
  double precision,intent(in) :: z,u,v

      fg3 = a * z*z*u*v*dexp(-aa*r2)

end function fg3
!============================================================================
function fh0(aa,z,a,r2)
!============================================================================
  implicit none
  complex(kind(0d0)) :: fh0
  complex(kind(0d0)),intent(in) :: a
  double precision,intent(in) :: aa
  double precision,intent(in) :: r2
  double precision,intent(in) :: z

      fh0 = a * z*z*z*z*z*dexp(-aa*r2)

end function fh0
!============================================================================
function fh1(aa,z,u,a,r2)
!============================================================================
  implicit none
  complex(kind(0d0)) :: fh1
  complex(kind(0d0)),intent(in) :: a
  double precision,intent(in) :: aa
  double precision,intent(in) :: r2
  double precision,intent(in) :: z,u

      fh1 = a * z*z*z*z*u*dexp(-aa*r2)

end function fh1
!============================================================================
function fh2(aa,z,u,a,r2)
!============================================================================
  implicit none
  complex(kind(0d0)) :: fh2
  complex(kind(0d0)),intent(in) :: a
  double precision,intent(in) :: aa
  double precision,intent(in) :: r2
  double precision,intent(in) :: z,u

      fh2 = a * z*z*z*u*u*dexp(-aa*r2)

end function fh2
!============================================================================
function fh3(aa,z,u,v,a,r2)
!============================================================================
  implicit none
  complex(kind(0d0)) :: fh3
  complex(kind(0d0)),intent(in) :: a
  double precision,intent(in) :: aa
  double precision,intent(in) :: r2
  double precision,intent(in) :: z,u,v

      fh3 = a * z*z*z*u*v*dexp(-aa*r2)

end function fh3
!============================================================================ 
function fh4(aa,z,u,v,a,r2)
!============================================================================
  implicit none
  complex(kind(0d0)) :: fh4
  complex(kind(0d0)),intent(in) :: a
  double precision,intent(in) :: aa
  double precision,intent(in) :: r2
  double precision,intent(in) :: z,u,v

      fh4 = a * z*z*u*u*v*dexp(-aa*r2)

end function fh4

