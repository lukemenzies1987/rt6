module loopngrpintmod
  implicit none

  interface
		subroutine DLSODE(F, neq, y, t, tout, itol, rtol, atol, itask,&
 & istate, iopt, rwork, lrw, iwork, liw, jac, mf)
      !EXTERNAL F, JAC
      !INTEGER, intent(in) :: neq, ITOL, ITASK, IOPT, LRW, LIW, MF
      !INTEGER, intent(inout) :: IWORK
      !integer, intent(inout) :: ISTATE
      !DOUBLE PRECISION, intent(out) :: Y, T, TOUT, RWORK
      !DOUBLE PRECISION, intent(in) :: RTOL, ATOL
      !DIMENSION NEQ(*), Y(*), RTOL(*), ATOL(*), RWORK(LRW), IWORK(LIW)
     ! SUBROUTINE DLSODE (F, NEQ, Y, T, TOUT, ITOL, RTOL, ATOL, ITASK,
     !1                  ISTATE, IOPT, RWORK, LRW, IWORK, LIW, JAC, MF)
      EXTERNAL F, JAC
      INTEGER NEQ, ITOL, ITASK, ISTATE, IOPT, LRW, IWORK, LIW, MF
      DOUBLE PRECISION Y, T, TOUT, RTOL, ATOL, RWORK
      DIMENSION NEQ(*), Y(*), RTOL(*), ATOL(*), RWORK(LRW), IWORK(LIW)
! must use assumed-size specification for f77 array
		end subroutine DLSODE

      subroutine DINTDY (T, K, YH, NYH, DKY, IFLAG)
        INTEGER K, NYH, IFLAG
        DOUBLE PRECISION T, YH, DKY
        DIMENSION YH(NYH,*), DKY(*)
      end subroutine

     subroutine DSTODE (NEQ, Y, YH, NYH, YH1, EWT, SAVF, ACOR, &
     &   WM, IWM, F, JAC, PJAC, SLVS)

      EXTERNAL F, JAC, PJAC, SLVS
      INTEGER NEQ, NYH, IWM
      DOUBLE PRECISION Y, YH, YH1, EWT, SAVF, ACOR, WM
      DIMENSION NEQ(*), Y(*), YH(NYH,*), YH1(*), EWT(*), SAVF(*), &
     &   ACOR(*), WM(*), IWM(*)
    end subroutine

    SUBROUTINE DCFODE (METH, ELCO, TESCO)
      INTEGER METH
      DOUBLE PRECISION ELCO, TESCO
      DIMENSION ELCO(13,12), TESCO(3,12)
      DIMENSION PC(12)
    end subroutine

    SUBROUTINE DPREPJ (NEQ, Y, YH, NYH, EWT, FTEM, SAVF, WM, IWM, &
     &   F, JAC)

      EXTERNAL F, JAC
      INTEGER NEQ, NYH, IWM
      DOUBLE PRECISION Y, YH, EWT, FTEM, SAVF, WM
      DIMENSION NEQ(*), Y(*), YH(NYH,*), EWT(*), FTEM(*), SAVF(*), &
     &   WM(*), IWM(*)
    end subroutine

    SUBROUTINE DSOLSY (WM, IWM, X, TEM)
      INTEGER IWM
      DOUBLE PRECISION WM, X, TEM
      DIMENSION WM(*), IWM(*), X(*), TEM(*)
    end subroutine

    SUBROUTINE DEWSET (N, ITOL, RTOL, ATOL, YCUR, EWT)
      INTEGER N, ITOL
      INTEGER I
      DOUBLE PRECISION RTOL, ATOL, YCUR, EWT
      DIMENSION RTOL(*), ATOL(*), YCUR(N), EWT(N)
    end subroutine 

    SUBROUTINE DSRCOM (RSAV, ISAV, JOB)
      INTEGER ISAV, JOB
      DOUBLE PRECISION RSAV
      DIMENSION RSAV(*), ISAV(*)
    end subroutine

    SUBROUTINE DGEFA (A, LDA, N, IPVT, INFO)
      INTEGER LDA,N,IPVT(*),INFO
      DOUBLE PRECISION A(LDA,*)
    end subroutine

    SUBROUTINE DGESL (A, LDA, N, IPVT, B, JOB)
      INTEGER LDA,N,IPVT(*),JOB
      DOUBLE PRECISION A(LDA,*),B(*)
    end subroutine


    SUBROUTINE DGBFA (ABD, LDA, N, ML, MU, IPVT, INFO)
      INTEGER LDA,N,ML,MU,IPVT(*),INFO
      DOUBLE PRECISION ABD(LDA,*)
    end subroutine

    SUBROUTINE DGBSL (ABD, LDA, N, ML, MU, IPVT, B, JOB)
      INTEGER LDA,N,ML,MU,IPVT(*),JOB
      DOUBLE PRECISION ABD(LDA,*),B(*)
    end subroutine

    SUBROUTINE DAXPY (N, DA, DX, INCX, DY, INCY)
      DOUBLE PRECISION DX(*), DY(*), DA
    end subroutine

    SUBROUTINE DSCAL (N, DA, DX, INCX)
      DOUBLE PRECISION DA, DX(*)
      INTEGER INCX, N
    end subroutine

    SUBROUTINE XERRWD (MSG, NMES, NERR, LEVEL, NI, I1, I2, NR, R1, R2)
      DOUBLE PRECISION R1, R2
      INTEGER NMES, NERR, LEVEL, NI, I1, I2, NR
      CHARACTER*(*) MSG
    end subroutine

    SUBROUTINE XSETF (MFLAG)
      INTEGER MFLAG
    end subroutine

    SUBROUTINE XSETUN (LUN)
      INTEGER LUN
    end subroutine

    SUBROUTINE DUMSUM(A,B,C)
      DOUBLE PRECISION A, B, C
    end subroutine

    INTEGER FUNCTION IXSAV (IPAR, IVALUE, ISET)
      LOGICAL ISET
      INTEGER IPAR, IVALUE
    end function

    DOUBLE PRECISION FUNCTION DDOT (N, DX, INCX, DY, INCY)
      DOUBLE PRECISION DX(*), DY(*)
    end function

    subroutine tolchangengrptest(y)
      use nuclearvars
      implicit none
      real(kind=dp) , dimension(neq), intent(in) :: y
    end subroutine

    subroutine clusterdensityngrp(y,nucdens)
      use nuclearvars
      implicit none
      real(kind=dp), intent(out):: nucdens
      real(kind=dp), dimension(neq), intent(in) :: y
    end subroutine

    subroutine clusterss(y,clss_V,clss_i,clss_He)
      use nuclearvars
      use twodimensionalvars
      implicit none
      real(kind=dp), intent(out) :: clss_V,clss_i,clss_He
      real(kind=dp), dimension(neq), intent(in) :: y
    end subroutine

    subroutine printintvalues(t,clss_V,clss_i,clss_He)
      use nuclearvars
      implicit none
      real(kind=dp), intent(in) :: t
      real(kind=dp), intent(in) :: clss_V,clss_i,clss_He
    end subroutine

    subroutine sinkstrcorrngrp(y)
      use nuclearvars
      use twodimensionalvars
      implicit none
      real(kind=dp), dimension(neq), intent(in) :: y
    end subroutine
	end interface
  

end module
