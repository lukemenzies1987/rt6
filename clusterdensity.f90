subroutine clusterdensity(y,nucdens)
!------------------------------------------------------------------
! Calculates the cluster density
!------------------------------------------------------------------ 
  use nuclearvars
  implicit none
  integer :: i,j,k
  real(kind=dp) :: nucdens,Fedens
  real(kind=dp), dimension(neq), intent(in) :: y
  real(kind=dp), dimension(1:xnumg,0:mnumg) :: sdf
  
  k=6
  Fedens=8.5E+28/1.E06
  nucdens=0.0_dp
  sdf=0.0
  !first half--------------------------------------
  do k=5+1,neq1
    i=idarrayi(k)
    j=idarrayj(k)
    sdf(i,j)=abs(y(k))     
    !nucdens=nucdens+abs(y(k))*Fedens
  end do
  !second half--------------------------------------
  do k=neq1+1,neq,3   
    i=idarrayi(k)
    j=idarrayj(k) 
    sdf(i,j)=abs(y(k))
    !nucdens=nucdens+abs(y(k))*gx(i)*gm(j)*Fedens
  end do 


!first half--------------------------------------
  do j=0,srtgrpm
    do i=1,srtgrpx
      nucdens=nucdens+sdf(i,j)
    end do
  end do 

!print*, neq1,'neq1'
  do j=0,srtgrpm
    do i=srtgrpx+1,xnumg


      nucdens=nucdens+sdf(i,j)*rgx(i)*rgm(j)
    end do
  end do 

!second half--------------------------------------
  do j=srtgrpm+1,mnumg
    do i=1,xnumg
      nucdens=nucdens+sdf(i,j)*rgx(i)*rgm(j)
    end do    
  end do 
  nucdens=nucdens-sdf(1,0)
  !do i=1,srtgrpx
   ! nucdens=nucdens+abs(sdf(i,0))*Fedens
  !end do 

  !do i=srtgrpx+1,xnumg
   ! nucdens=nucdens+abs(sdf(i,0))*gx(i)*Fedens
  !end do
  nucdens=nucdens/vol
  return
end subroutine clusterdensity

subroutine clusterdensity1d(y,nucdens)
!------------------------------------------------------------------
! Calculates the cluster density in 1d
!------------------------------------------------------------------ 
  use nuclearvars
  implicit none
  integer :: i,j,k
  real(kind=dp) :: nucdens,Fedens
  real(kind=dp), dimension(neq) :: y
  real(kind=dp), dimension(1:xnumg) :: sdfv
  
  
  Fedens=8.5E+28/1.E06
  nucdens=0.0_dp
  sdfv=0.0
  !first half--------------------------------------
  do k=1,neq1
    i=idarrayi(k)
    sdfv(i)=abs(y(k))
!print*,k,i
  end do 

  do k=neq1+1,neq2,2
    i=idarrayi(k)
    sdfv(i)=abs(y(k))
!print*,k,i
  end do 


  do i=2,srtgrpx
    nucdens=nucdens+abs(sdfv(i))
  end do 

  do i=srtgrpx+1,xnumg
    nucdens=nucdens+abs(sdfv(i))*gx(i)
  end do

  !do i=2,srtgrpi
  !  intdensity=intdensity+abs(sdfi(i))
  !end do 

  !do i=srtgrpi+1,inumg
  !  intdensity=intdensity+ abs(sdfi(i))*gi(i)
  !end do

  nucdens=nucdens/vol
  return
end subroutine clusterdensity1d
