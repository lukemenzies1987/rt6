module twodimensionalvars 
  use kindmod
  implicit none
  !integer, parameter :: dp=selected_real_kind(15,300)

  real(kind=nag_wp), dimension(:,:), allocatable :: Qx,Qvarr,f,df,Pvarr,Qm,Qiarr,PHearr,Qg
  real(kind=nag_wp), dimension(:,:), allocatable :: fPvarr,fQvarr,fQm,fQiarr,fPHearr,fQg
  real(kind=nag_wp), allocatable, dimension(:,:) :: L0,L1x,L1m
  real(kind=nag_wp), allocatable, dimension(:,:) :: bindengV, bindengHe
  real(kind=nag_wp), allocatable, dimension(:,:) :: Pvarr0,PHearr0,Qiarr0,& 
& fPvarr0,fPHearr0,fQiarr0
  real(kind=nag_wp), allocatable, dimension(:) :: zerox,zerom


end module twodimensionalvars
