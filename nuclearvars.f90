module nuclearvars
  use kindmod
  implicit none
  !integer, parameter :: dp=selected_real_kind(15,300)
  real(kind=nag_wp),parameter :: Pi = 3.1415927_dp, kb =8.6173324*10.0**(-5.0)
  real(kind=nag_wp) :: a,temp,seng,Ev,D0He, D0i,D0v,vol,nui,rhod,ksqGb,start
  real(kind=nag_wp) :: Cv0,Ci0,Cv,Ci,Cg,CHe,CGb,dCv,dCi,dCg,dCHe,dCGb,fCi,fCv,fCHe
  real(kind=nag_wp) :: AHe,Emv,Emi,EmHe,epsilonr,epsiloni,epsilonv
  real(kind=nag_wp) :: Geni,Gv,Gg,GHe,nuHe,ZHe,Zi,Zv,Gdpa,cldens,disl0,CHe0
  real(kind=nag_wp) :: Di,Dv,DHe,Dg,P,rcap,sv_dis,sHe_dis,si_dis,sv_gb,si_gb,sHe_gb
  integer :: maxgv,printint
  real(kind=nag_wp), allocatable, dimension(:) :: genv
  integer :: mnum,xnum,m0,neq,counter

  real(kind=nag_wp):: initt,dt,tmax,tgdi,tgdid,disl,tgenHe,tGBHe

  real(kind=nag_wp)::rtol,g00,g11
  real(kind=nag_wp), allocatable, dimension(:) :: rwork, atol,rgx,rgm
  integer, allocatable, dimension(:) :: iwork
  integer :: istate,itol,mf,liw,lrw,iopt,itask,nsteps

  integer :: noofv,noofa
  real(kind=nag_wp), allocatable, dimension(:,:) :: gotEvals
  character*30:: filename
  character*6 :: foldername

  real(kind=nag_wp):: atolCv,atolCi,atolCHe,atolCGb,atolCg,atolClusters

  integer :: flaggrouping,flagnongrouping,flagEbtable,flagstatictimestep,& 
& flagvaryingtimestep,flagVoids,flagClusters,flagsinkstrengthcorrection

  integer :: xnumg,mnumg,srtgrpx, srtgrpm, gxw,gmw,unevengx, &
& unevengm,neq1,neq2, neqend,minnum
  integer, allocatable, dimension(:) :: gx, gm , g1xi, &
& g1mi
  real(kind=nag_wp) :: dr,drm,kgCg
  real(kind=nag_wp), allocatable, dimension(:) :: itermx,itermm, &
& meanxarr,meanmarr,rg1xi,rg1mi,r
  integer,allocatable,dimension(:) :: idarrayi,idarrayj
  real(kind=nag_wp) :: dt0, dtmax, dtconvrt

  real(kind=nag_wp) :: sigmavNv, fstterm,secndterm,kgsq,Rg,l,&
& disabs,Qgc

  logical :: petzld
  integer :: maxord,maxstp,mxhnil,nwkjac,ifail,sdysav,ldysav, &
& itrace,iset,outchn
  real(kind=nag_wp) :: hmin,hmax,h0,tcrit,aamin
  real(kind=nag_wp), allocatable, dimension(:) :: con,ydot,wkjac, &
& yinit
  real(kind=nag_wp), allocatable, dimension(:,:) :: ysav
  real(kind=nag_wp) :: hu,h,tcur,tolsf
  integer :: nst,nre,nje,nqu,nq,niter,imxer
  logical, allocatable, dimension(:) :: algequ


end module nuclearvars
