module onedimensionalvars 
  use kindmod
  implicit none
  !integer, parameter :: dp=selected_real_kind(15,300)

  real(kind=dp), dimension(:), allocatable :: Qx,Qvarr,fx,fi,Pvarr,Qm,Qiarr,PHearr,&
& Qi,Piarr,rgi,ril,gi,bindingV,Qg
  real(kind=dp), allocatable, dimension(:) :: L0x,L0i,L1x,L1i,g1ii,rg1ii,meaniarr, &
& itermi,itermg1
  integer :: Xg,inumg,srtgrpi,b,inum
  real(kind=dp) :: Zil,Zvl,dri
  
  integer, allocatable, dimension(:) :: g1,fcnt,g1p1,fcntp1,g1m1,fcntm1
  real(kind=dp),dimension(:), allocatable :: Geniloop
  integer :: maxgi,Xgtrt,grpnum



end module onedimensionalvars
