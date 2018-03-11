subroutine speciesgenerationHe(GHe,t)
!------------------------------------------------------------------
! Helium generation rate 
!------------------------------------------------------------------ 
  use kindmod
  implicit none
  integer, parameter :: clen=19!10 !change depending on if one wishes to use FW or 
!blanket region
  integer :: i
  real(kind=dp), intent(out) :: GHe
  real(kind=dp), intent(in) :: t 
  real(kind=dp), dimension(clen) :: c

  !Helium generation rate for the FW region
  c = (/4.06377395e-014_dp,   4.69588506e-006_dp,   1.87574177e-015_dp,&
      &   5.87396012e-022_dp,  -1.42713308e-028_dp,   1.66184024e-035_dp,&
      &  -1.05651484e-042_dp,   3.81345635e-050_dp,  -7.71461252e-058_dp,&
      &   7.84716062e-066_dp,  -2.11507413e-074_dp,  -2.00414380e-082_dp,&
      &   3.84188260e-091_dp,   7.58895690e-099_dp,   2.28883851e-107_dp,&
      &  -1.33645508e-115_dp,  -1.70262506e-123_dp,  -5.51091497e-132_dp,&
      &   6.52356104e-140_dp/)

 ! Helium generation rate for the blanket region
 ! c = (/ -1.21516554e-04_dp,   7.68530646e-07_dp,   1.81956618e-16_dp, &
 !    &    7.71157133e-24_dp,  -4.02200540e-31_dp,   1.03367381e-38_dp, &
 !    &   -1.50734596e-46_dp,   1.24917198e-54_dp,  -5.45378195e-63_dp, &
 !    &    9.69469897e-72_dp/)



  GHe=0.0_dp

  do i=2,clen
    GHe = GHe + (real(i,kind=dp)-1.0_dp)*t**(real(i,kind=dp)-2.0_dp)*c(i)
  end do

  !Fixed value. Change manually and comment expression above if one wishes.
  !GHe=1.2E-03_dp! !When T=623K, dpa
  
  GHe=GHe*1.E-06_dp

  return
end subroutine speciesgenerationHe

