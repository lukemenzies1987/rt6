subroutine getparameters

  use nuclearvars
  implicit none
  integer :: nunit,ios

  namelist/parameters/a,temp,seng,Ev,mnum,xnum,m0,D0He,D0i,D0v, &
  &  vol,nui,rhod,ksqGb,Emv,Emi,EmHe,AHe,epsilonr,epsilonv,epsiloni,Gdpa,nuHe,ZHe, &
  &  Zv,Zi,initt,dt,tmax,Rg,l,disabs,rcap,maxgv,g00,g11
  namelist/init/filename,Cv0,Ci0,GHe,Geni,Gv
  namelist/grouping/srtgrpx,srtgrpm,gxw,gmw,dr,drm
  namelist/timestep/dt0,dtmax,dtconvrt
  namelist/tolerance/rtol,atolCv,atolCi,atolCHe,atolCGb,atolCg,atolClusters
  namelist/settings/flaggrouping,flagEbtable,flagstatictimestep,& 
& flagvaryingtimestep,flagVoids,flagClusters,flagsinkstrengthcorrection
  ! Read in input data file
  nunit=18
  open(unit=nunit,file='param.in',status='old',iostat=ios)
    if(ios.ne.0) then
       write(6,*) 'problem opening file parameter file'
       stop
    endif
    read(nunit,parameters)
    read(nunit,init)
    read(nunit,grouping)
    read(nunit,timestep)
    read(nunit,tolerance)
    read(nunit,settings)
    !write(6,struct)
    close(nunit)
  !print*, diff(temp,Emv,D0v)

  if (ios .eq. 0) print*, 'Read in parameters from file'
  if (flagVoids .eq. 1 .and. flagClusters .eq. 1) then
    print*, 'Both options for a void distribution and' ,& 
& ' cluster distribution were chosen... choose one'
    stop
  end if
  if (flagVoids .eq. 0 .and. flagClusters .eq. 0) then
    print*, 'Choose either a void distribution or ',&
& 'cluster distribution'
    stop
  end if 
  return
end subroutine getparameters
