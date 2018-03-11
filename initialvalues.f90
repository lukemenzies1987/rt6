subroutine initialvalues(Fedens,intdens,Imprate,error,prevt,y,t)
  use nuclearvars
  use twodimensionalvars
  implicit none
  real(kind=dp) :: Fedens, intdens,Imprate,invdt,corr,prevt,error,t

  real(kind=dp), dimension(neq) :: y

  
  Fedens=1._dp/vol
  corr=1.!1.E+11
  t =    initt              ! initial value for variable t
  
  !*** the itegration limit and step                ! step size for integration
  !dt  = 0.0002_dp
  !tmax = 1000_dp               ! integrate till tmax 
  invdt= 1._dp/dt

  nsteps = tmax*invdt
  !*** end of initial data

  istate=1
  itol=2
  
  Imprate=0.0
  itask=1
  istate=1
  iopt=1

  itrace = 0
  maxstp=10000
  mxhnil=5
  petzld=.true.
  hmin=0.0!1.0E-10
  hmax=0.0!10.0
  h0=0.0
  con(1:6) = 0.0_nag_wp
  iset=1

  aamin=1.E-35
  mf=22
  tgenHe=0.0
  !Qg=0.0
  CHe0=1.E-12
  disl=0.0_dp
  counter=0  
  y=0.0_dp
  Ci=aamin
  CHe=aamin
  Cg=aamin
  yinit(1)=Cv
  yinit(2)=Ci
  yinit(3)=CHe
  yinit(4)=CGb
  yinit(5)=Cg
  yinit(6)=Cv  
  yinit(minnum)=aamin
  rwork=0.0_dp
  iwork=0
  tgdi=0.0
  disl0=DHe*CHe0*ZHe*rhod
  atol=atolClusters
  atol(1)= atolCv 
  atol(2)= atolCi
  atol(3)= atolCHe
  atol(4)= atolCGb
  atol(5)= atolCg
  prevt=0.0_dp
  error=0.0_dp
  !rwork(6)=1.E+10
  !iwork(6)=10000
  !iwork(7)=12
  !rwork(6)=1.E10
  ! integration of ODEs


  return
end subroutine initialvalues

subroutine initialvalues1d(Fedens,intdens,Imprate,error,prevt,y,t)
  use nuclearvars
  use onedimensionalvars
  implicit none
  real(kind=dp) :: Fedens, intdens,Imprate,invdt,corr,prevt,error,t

  real(kind=dp), dimension(neq) :: y

  
  Fedens=1._dp/vol
  corr=1.!1.E+11
  t =    initt              ! initial value for variable t
  
  !*** the itegration limit and step                ! step size for integration
  !dt  = 0.0002_dp
  !tmax = 1000_dp               ! integrate till tmax 
  invdt= 1._dp/dt

  !nsteps = tmax*invdt
  !*** end of initial data

  istate=1
  itol=2
  
  Imprate=0.0
  itask=1
  istate=1
  iopt=1

  mf=22


  disl=0.0
  counter=0  
  y=0.0_dp

  yinit(1)=Cv
  yinit(neq2+1)=Ci
  yinit(neq2+2)=Cg  
  rwork=0.0_dp
  iwork=0
  tgdi=0.0
  atol=atolClusters
  atol(1)= atolCv 
  atol(neq2+1)= atolCi
  atol(neq2+2)= atolCGb

  prevt=0.0
  error=0.0

  !IWORK(6)=1000
  !iwork(7)=12
  !rwork(6)=1.E10
  ! integration of ODEs


  return
end subroutine initialvalues1d
