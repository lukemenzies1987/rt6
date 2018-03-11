subroutine loop1(Fedens,intdens,Imprate,error,prevt,y,t)
!------------------------------------------------------------------
! Loop using the 2d grouping method to solve rate equations
!------------------------------------------------------------------  
  use nag_library, Only: d02nbf, d02nby, d02nbz, d02nsf, d02nvf, d02nyf,   &
                             nag_wp, x04abf
  use nuclearvars
  implicit none
  external funcall,funcallinit,jac,monitr
  real(kind=dp), dimension(neq) :: y
  real(kind=dp) :: Fedens, intdens,Imprate,invdt,corr,prevt, &
& invtime, error,t,tout,tstage,clss_V,clss_i,clss_He
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

    subroutine sinkstrcorr(y)
      use nuclearvars
      use twodimensionalvars
      implicit none      
      real(kind=dp), dimension(neq),intent(in) :: y
    end subroutine

    subroutine tolchange(y)
      use nuclearvars
      implicit none
      real(kind=dp) , dimension(neq), intent(in) :: y
    end subroutine

    subroutine clusterdensity(y,nucdens)
      use nuclearvars
      implicit none
      real(kind=dp) , intent(out) :: nucdens
      real(kind=dp), dimension(neq), intent(in) :: y
    end subroutine

    subroutine clusterss(y,clss_V,clss_i,clss_He)
      use nuclearvars
      use twodimensionalvars
      implicit none
      real(kind=dp), intent(out) :: clss_V,clss_i,clss_He
      real(kind=dp), dimension(neq), intent(in) :: y
    end subroutine

    subroutine printvalues(t,clss_V,clss_i,clss_He)
      use nuclearvars
      use twodimensionalvars
      implicit none
      real(kind=dp) :: clss_V,clss_i,clss_He
      real(kind=dp), dimension(neq) :: t
    end subroutine
	end interface
  !Nag library ode solver setup
  neq11(1) =neq
  rtol11(1)=rtol
  call d02nsf(neq,neq,'Numerical',nwkjac,rwork,ifail)


  call d02nvf(neq,sdysav,maxord,'Newton',petzld,con,tcrit,hmin,hmax,h0, &
          maxstp,mxhnil,'Average-L2',rwork,ifail)
  !Initial preliminary iteration
  !-----------------------------------------------
  invtime=1._dp/(t+1.E-30)

  tout =t+dt0+dt
  y(1:neq) = yinit(1:neq)
  call d02nbf(neq,ldysav,t,tout,y,ydot,rwork,rtol11,atol,itol,iwork, &
      funcallinit,ysav,sdysav,d02nbz,wkjac,nwkjac,monitr,itask,itrace,ifail)

  tgdi=tgdi+disl*(dt0+dt)
  tgenHe=tgenHe + GHe*(dt0+dt)
  Imprate=GHe*tout
  tGBHe=0.0
  tstage=nint(t*10000./tmax)
  !call tolchangengrp(y)
  !----------------------------------------------
  prevt=tstage
  counter=counter+1
  !Main loop
  !-----------------------------------
  do while (t <= tmax)

    invtime=1._dp/(t+1.E-30)
    !dt=0.0001
    tout =t+dt0+dt

   !odepack solver
   ! call DLSODE(funcallngrpss, neq11, y, t, tout, itol, rtol11, atol, itask,&
 !  & istate, iopt, rwork, lrw, iwork, liw, jac, mf)

    !Nag Fortran library solver. 
    call d02nbf(neq,ldysav,t,tout,y,ydot,rwork,rtol11,atol,itol,iwork, &
      funcall,ysav,sdysav,d02nbz,wkjac,nwkjac,monitr,itask,itrace,ifail)
    if (ifail .gt. 0) then 
      call d02nyf(neq,neq,hu,h,tcur,tolsf,rwork,nst,nre,nje,nqu,nq,niter, &
            imxer,algequ,iwork,ifail)
      print*, IMXER,'error comp'
    end if

    !very basic solver. Only used for debugging purposes.
    !call euler(neq,y,ydot,t )

    
    tgdi=tgdi+disl*(dt0+dt) !Dislocation sink strength increase
    tgenHe=tgenHe+GHe*(dt0+dt) !He accumilation concentration
    Imprate=GHe*tout
    tGBHe=tGBHe+ DHe*CHe*sHe_gb*(dt+dt0) 
     
    tstage=nint(t*10000./tmax)
    if (tstage.ge.prevt+1 .and. t .lt. tmax) then
      call clusterdensityngrp(y,cldens) !Determines the cluster density. 
      Cv=y(1) !vacancy concentration
      Ci=y(2) !SIA concentration
      CHe=y(3) !He concentration
      CGb=y(4) !He grain boundary concentration
      Cg=y(5) ! SIA glissile concentration
      write(7,102) t, Cv,Ci,CHe,CGb*1.E+06,Cg,cldens,Imprate*1.E+06
      print 22, t/tmax*100_dp 
      !call wtf2(y)

      !Determines sink strength of clusters
      call clusterss(y,clss_V,clss_i,clss_He)
      !Calculates the estimates numerical error
      do i =1,neq
        error = error + abs(rwork(lrw-neq+i))
      end do
      
      dt = dtmax*exp(-tmax*dtconvrt*invtime) !Evolving timestep
      call printintvalues(t,clss_V,clss_i,clss_He) !print sink strength values
      call tolchangengrptest(y)!test(y)  !Change the tolerance 
      call sinkstrcorrngrp(y)  !High density sink strength correction. 
      !print*, error,' Magnitude of error in vector array'
      
    end if

    call tolchangengrptest(y)
    prevt=tstage
    counter=counter+1

  end do

  !call wtf2(y)
  22 format(f7.2,'% Complete ')
  100 format(5x,'t',15x,'Cv',15x,'Ci',15x,'CHe',15x,'CGb',15x,'CGc',15x,'Cluster Density') 
  102 format(8(1pe16.7))

  return

end subroutine loop1
