module clusterintmod
  implicit none

  interface 
    subroutine initialvalues(Fedens,intdens,Imprate,error,prevt,y,t)
      use nuclearvars
      implicit none
      real(kind=dp) :: Fedens, intdens,Imprate,corr,prevt,error,t
      real(kind=dp), dimension(neq) :: y
    end subroutine

    subroutine loop1ngrp(Fedens,intdens,Imprate,error,prevt,y,t)
      use nuclearvars
      implicit none
      external funcallngrp,funcallngrpinit,jac
      real(kind=dp), dimension(neq) :: y
      real(kind=dp) :: Fedens,intdens,Imprate,prevt,error,t
    end subroutine

    subroutine loop1(Fedens,intdens,Imprate,error,prevt,y,t)
      use nuclearvars
      implicit none
      external funcallngrp,funcallngrpinit,jac
      real(kind=dp), dimension(neq):: y
      real(kind=dp) :: Fedens,intdens,Imprate,prevt,error,t
    end subroutine

    subroutine loop2ngrp(Fedens,intdens,Imprate,error,prevt,y,t)
      use nuclearvars
      implicit none
      external funcallngrp,funcallngrpinit,jac
      real(kind=dp), dimension(neq) :: y
      real(kind=dp) :: Fedens,intdens,Imprate,prevt,error,t
    end subroutine

    subroutine loop3(Fedens,intdens,Imprate,error,prevt,y,t)
      use nuclearvars
      implicit none
      external funcallngrp,funcallngrpinit,jac
      real(kind=dp), dimension(neq):: y
      real(kind=dp) :: Fedens,intdens,Imprate,prevt,error,t
    end subroutine

    subroutine loop2(Fedens,intdens,Imprate,error,prevt,y,t)
      use nuclearvars
      implicit none
      external funcallngrp,funcallngrpinit,jac
      real(kind=dp), dimension(neq) :: y
      real(kind=dp) :: Fedens,intdens,Imprate,prevt,error,t
    end subroutine

    subroutine dof()
    end subroutine

    subroutine initialsetup1d
    end subroutine

    subroutine initialsetup
    end subroutine

    subroutine speciesgeneration
    end subroutine

    subroutine printinputs
    end subroutine

    subroutine printabsemissrates2d
    end subroutine

    subroutine printabsemissrates1d
    end subroutine

    subroutine initialvalues1d(Fedens,intdens,Imprate,error,prevt,y,t)
      use nuclearvars
      use onedimensionalvars
      implicit none
      real(kind=dp) :: Fedens, intdens,Imprate,prevt,error,t
      real(kind=dp), dimension(neq) :: y
    end subroutine
  end interface
end module clusterintmod
