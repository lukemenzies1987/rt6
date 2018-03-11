subroutine loop2ngrp(Fedens,intdens,Imprate,error,prevt,y,t)
!------------------------------------------------------------------
! Loop using the 2d grouping method to solve rate equations 
!------------------------------------------------------------------  
  use nag_library, Only: d02nbf, d02nby, d02nbz, d02nsf, d02nvf, d02nyf,   &
                             nag_wp, x04abf
  use nuclearvars
  use loopngrpintmod
  implicit none
  external funcallngrp,funcallngrpinit,jac,monitr
  real(kind=dp), dimension(neq) :: y
  real(kind=dp) :: Fedens, intdens,Imprate,invdt,corr,prevt, &
& invtime, error,t,tout,tstage,clss_V,clss_i,clss_He
  integer :: i
  integer,dimension(1) :: neq11
  double precision, dimension(1) :: rtol11
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
      funcallngrpinit,ysav,sdysav,d02nbz,wkjac,nwkjac,monitr,itask,itrace,ifail)
  !call DLSODE(funcallngrpinit, neq11, y, t, tout, itol, rtol11, atol, itask,&
 !& istate, iopt, rwork, lrw, iwork, liw, jac, mf)

  tgdi=tgdi+disl*(dt0+dt)
  Imprate=GHe*tout
  
  tstage=nint(t*10000./tmax)

  prevt=tstage
  counter=counter+1
  !-----------------------------------
  ifail = -1
  !Main loop
  !-----------------------------------
  do while (t <= tmax)


    tout =t+dt0+dt
   !odepack solver
   ! call DLSODE(funcallngrp, neq11, y, t, tout, itol, rtol11, atol, itask,&
   !& istate, iopt, rwork, lrw, iwork, liw, jac, mf)
    call d02nbf(neq,ldysav,t,tout,y,ydot,rwork,rtol11,atol,itol,iwork, &
      funcallngrp,ysav,sdysav,d02nbz,wkjac,nwkjac,monitr,itask,itrace,ifail)
    if (ifail .gt. 0) then 
      call d02nyf(neq,neq,hu,h,tcur,tolsf,rwork,nst,nre,nje,nqu,nq,niter, &
            imxer,algequ,iwork,ifail)
      print*, IMXER,'error comp'
    end if
    tgdi=tgdi+disl*(dt0+dt)
    tgenHe=tgenHe + GHe*(dt0+dt)
    Imprate=GHe*t
  
     
    tstage=nint(t*10000./tmax)
    if (tstage.ge.prevt+1 .and. t .lt. tmax) then
      call clusterdensityngrp(y,cldens)
      Cv=y(1)
      Ci=y(2)
      CHe=y(3)
      CGb=y(4)
      Cg=y(5)
      write(7,102) t, Cv,Ci,CHe,CGb*1.E+06,Cg,cldens,Imprate*1.E+06
      print 22, t/tmax*100_dp 
      call wtf2(y)
      !do i =1,neq
      !  error = error + abs(rwork(lrw-neq+i))
      !end do
      invtime=1._dp/(t+1.E-30)
      dt = dtmax*dexp(-tmax*dtconvrt*invtime) !Evolving timestep

      call clusterss(y,clss_V,clss_i,clss_He) !Determines sink strength of clusters
      call printintvalues(t,clss_V,clss_i,clss_He) !print sink strength values
      call tolchangengrptest(y) !Change the tolerance 

    end if

    prevt=tstage
    counter=counter+1

  end do
  !call wtf2(y)
  22 format(f7.2,'% Complete ')
  100 format(5x,'t',15x,'Cv',15x,'Ci',15x,'CHe',15x,'CGb',15x,'CGc',15x,'Cluster Density') 
  102 format(8(1pe16.7))

  return
end subroutine loop2ngrp
