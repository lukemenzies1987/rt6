subroutine loop3(Fedens,intdens,Imprate,error,prevt,y,t)
!------------------------------------------------------------------
! Loop using the 1d method to solve rate equations
!------------------------------------------------------------------  
  use nuclearvars
  implicit none
  external funcall1d,jac
  real(kind=dp), dimension(neq) :: y
  real(kind=dp) :: Fedens, intdens,Imprate,invdt,corr,prevt, &
& invtime, error,t,tout,tstage
  integer :: i
  integer,dimension(1) :: neq11
  real(kind=dp), dimension(1) :: rtol11

  interface
		subroutine DLSODE(F, neq, y, t, tout, itol, rtol, atol, itask,&
 & istate, iopt, rwork, lrw, iwork, liw, jac, mf)
      EXTERNAL F, JAC
      INTEGER , intent(in) :: neq, ITOL, ITASK, IOPT, LRW, LIW, MF
      INTEGER, intent(inout) :: IWORK
      integer, intent(inout) :: ISTATE
      DOUBLE PRECISION, intent(inout) :: Y, T, TOUT, RWORK
      DOUBLE PRECISION, intent(in) :: RTOL, ATOL
      DIMENSION NEQ(*), Y(*), RTOL(*), ATOL(*), RWORK(LRW), IWORK(LIW)
! must use assumed-size specification for f77 array
		end subroutine DLSODE

    subroutine sinkstrcorr1d(y)
      use nuclearvars
      use onedimensionalvars
      implicit none      
      real(kind=dp), dimension(neq),intent(in) :: y
    end subroutine

    subroutine tolchange1dtest(y)
      use nuclearvars
      implicit none
      real(kind=dp) , dimension(neq), intent(in) :: y
    end subroutine

    subroutine clusterdensity1d(y,nucdens)
      use nuclearvars
      implicit none
      real(kind=dp) , intent(out) :: nucdens
      real(kind=dp), dimension(neq), intent(in) :: y
    end subroutine


    subroutine printintvalues1d(t)
      use nuclearvars
      implicit none
      real(kind=dp) :: t
    end subroutine
	end interface
  neq11(1) =neq
  rtol11(1)=rtol

  do while (t <= tmax)
    invtime=1._dp/(t+1.E-30)
    dt = dtmax*exp(-tmax*dtconvrt*invtime)
    tout =t+dt0+dt
    !if (istate .lt. 0) then
       !call wtf1d(y)
    !endif
    call DLSODE(funcall1d, neq11, y, t, tout, itol, rtol11, atol, itask,&
   & istate, iopt, rwork, lrw, iwork, liw, jac, mf)
    Cv=y(1)
    Ci=y(neq2+1)
    Cg=y(neq2+2)
    call sinkstrcorr1d(y)
    tstage=nint(t*10000./tmax)
    if (tstage.ge.prevt+1 .and. t .lt. tmax .or. t .lt. 2.) then
      call clusterdensity1d(y,cldens)
      call printintvalues1d(t)
      write(7,102) t*Gdpa, Cv,Ci,Cg,cldens!,intdens/vol
      print 22, t/tmax*100_dp 
      do i =1,neq
        error = error + abs(rwork(lrw-neq+i))
       ! write(321,*) y(i),ydot(i)
      end do
      print*, error
    end if
    call tolchange1dtest(y)

    prevt=tstage
    counter=counter+1
  end do 

  22 format(f7.2,'% Complete ')
  100 format(5x,'t',15x,'Cv',15x,'Ci',15x,'CHe',15x,'CGb',15x,'CGc',15x,'Cluster Density') 
  102 format(6(1pe16.7))

  return
end subroutine loop3
