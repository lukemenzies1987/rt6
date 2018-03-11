subroutine clusterdensityngrp(y,nucdens)
!------------------------------------------------------------------
! Calculates the cluster density non-grouping
!------------------------------------------------------------------ 
  use nuclearvars
  implicit none
  integer :: i,j,k
  real(kind=dp) :: nucdens
  real(kind=dp), dimension(neq) :: y
  real(kind=dp), dimension(1:xnumg,0:mnumg) :: sdf
  
  

  nucdens=0.0_dp
  sdf=0.0
  !first half--------------------------------------
  do k=5+2,neq
    i=idarrayi(k)
    j=idarrayj(k)
    sdf(i,j)=abs(y(k))     
    nucdens=nucdens+sdf(i,j)
  end do


!first half--------------------------------------
  !do j=0,mnumg
  !  do i=1,xnumg
  !    nucdens=nucdens+sdf(i,j)
  !  end do
  !end do 


  !nucdens=nucdens-sdf(1,0)
  !do i=1,srtgrpx
   ! nucdens=nucdens+abs(sdf(i,0))*Fedens
  !end do 

  !do i=srtgrpx+1,xnumg
   ! nucdens=nucdens+abs(sdf(i,0))*gx(i)*Fedens
  !end do
  nucdens=nucdens/vol
  return
end subroutine clusterdensityngrp
