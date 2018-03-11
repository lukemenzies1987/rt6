subroutine clusters
!------------------------------------------------------------------
! The main subroutine containing the setup of the parameter and the loop used 
! for the calculations. 
!------------------------------------------------------------------ 
  use clusterintmod
  use nuclearvars
  implicit none

  real(kind=dp) :: efinish, etime, finish, &
& invdt,t,tout,invtime,tstage,prevt, error,Fedens,corr, &
& Imprate,intdens
  real(kind=dp), allocatable, dimension(:) :: y
  integer(kind=ip) :: i,j ,test,checkst,lerrind,rlrw,k



  !Deletes old files
  call dof()

  if (flagvoids .eq. 1) then 
    call initialsetup1d !1 dimensional setup
  else
    call initialsetup !2 dimensional setup
  end if 
  !This determines the (fixed) species generation rate. The user may wish to modify this. 
  call speciesgeneration 
  !neq=neq+1000
  lrw=50+4*neq!(22+9*neq+neq*neq)
  liw=(20+neq) 
  nwkjac=neq*(neq+1)
  ldysav = neq
  maxord=5
  sdysav = maxord + 1

  allocate(rwork(lrw),iwork(23),con(6),y(neq),yinit(neq),ydot(neq), & 
&  ysav(ldysav,sdysav),wkjac(nwkjac),algequ(neq),stat=checkst)
  !neq=neq-1000
  if(checkst.ne.0) then
     print*, 'Problem allocating arrays'
     stop
  end if
  ! Sets up the initial values for the one dimensional and two dimensional case. 
  if (flagClusters .eq. 1) then
    call initialvalues(Fedens,intdens,Imprate,error,prevt,y,t)
  else 
    call initialvalues1d(Fedens,intdens,Imprate,error,prevt,y,t)
  end if
  !The output folder that all the output files are stored in
  foldername="Output"
  call system( "mkdir -p " //foldername)
  ! This subroutine prints the appropriate simulation information into a chosen output folder. 
  call printinputs
  !Prints out the appropriate emission rates for the 1 dimensional and 2 dimensional case. 
  if(flagClusters .eq. 1) call printabsemissrates2d
  if(flagVoids .eq. 1) call printabsemissrates1d
  !Intermediatevalues are printed during the simulation where they contain the evolving sink
  !strength terms. 
  open (unit=printint, file=foldername//'/'//'Intermediatevalues.dat')
  !Clusters is the main output file that contains the concentrations for different species
  !as well as the cluster density. 
  open (unit=7, file=foldername//'/'//'Clusters.dat')
  
  print*, 'Starting...'
  !* print the header and initial conditions
  write (7,*) '  Nucleation of Helium within alpha-Fe within DEMO '
  write (7,*) '     Method: d02nbf   ' !chande to LSODE if using ODE pack

  if (flagClusters .eq. 1) then 
    write (printint,11)
    write (7,100)
    write (7,102) t,Cv,Ci,CHe,CGb,Cg,Cv*Fedens,Imprate !, xi(3), xi(4)
  end if
  if (flagVoids .eq. 1) then 
    write (printint, 12)
    write (7,101) 
    write (7,103) t,Cv,Ci,Cg,Cv*Fedens !, xi(3), xi(4)
  end if  
  !The main loops
  if (flagClusters .eq. 1 ) then
    if (flaggrouping .eq. 1) then 
      if (flagsinkstrengthcorrection .eq. 1)then 
        call loop1(Fedens,intdens,Imprate,error,prevt,y,t) !Grouping method
      else if (flagsinkstrengthcorrection .eq. 0) then
        call loop2(Fedens,intdens,Imprate,error,prevt,y,t) !Grouping method with
  !high density sink strength corrections.
      end if
    else
      if (flagsinkstrengthcorrection .eq. 1)then 
        call loop1ngrp(Fedens,intdens,Imprate,error,prevt,y,t) !Non-grouping method
      else if (flagsinkstrengthcorrection .eq. 0) then
        call loop2ngrp(Fedens,intdens,Imprate,error,prevt,y,t) !Non-grouping method
  !with high density sink strength corrections. 
      end if
    end if
  end if 


  if (flagVoids .eq. 1) then 
    call loop3(Fedens,intdens,Imprate,error,prevt,y,t) !1d grouping method (to use the 
  !non-grouping method, change the srtgrpx parameter to the size xnum-1.  
  end if 

  !Used to print out the size distribution function (comment or uncomment). This can also be 
  !placed in the loops but be warned, this slows the simulation time dramatically as well as
  !producing a large file. 

  !call wtffv(y)
  !call wtf2(y)

  call cpu_time(finish)
  print*, 'Finished!'
  print '("The simulation completed in ",f20.2," seconds.")',finish-start
  close(7)
  close(printint)
  11 format('  time(s)     ','K^2_GB Vac      ','K^2_GB SIA      ','K^2_GB He      ', & 
& ' K^2_dis Vac     ','K^2_dis SIA     ','K^2_dis He   ')
  12 format(' time(s)     ','K^2_GB Vac      ','K^2_GB SIA      ', & 
& ' K^2_dis Vac     ','K^2_dis SIA     ')
  22 format(f7.2,'% Complete ')
  100 format(5x,'t',15x,'Cv',15x,'Ci',15x,'CHe',15x,'CGb',15x,'CGc',15x,'Cluster Density') 
  101 format(5x,'t',15x,'Cv',15x,'Ci',15x,'CGc',15x,'Cluster Density') 
  102 format(8(1pe16.7))
  103 format(5(1pe16.7))
  deallocate(rwork,iwork,wkjac,y,ydot,con,yinit,ysav,algequ)
  return
end subroutine clusters

