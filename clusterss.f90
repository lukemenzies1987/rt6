subroutine clusterss(y,clss_V,clss_i,clss_He)
  use nuclearvars
  use twodimensionalvars
  implicit none
  real(kind=dp) :: clss_V,clss_i,clss_He
  real(kind=dp), dimension(neq) :: y
  real(kind=dp), dimension(1:xnumg,0:mnumg) :: sdf
  integer ::i,j,k

  clss_V=0.0
  clss_i=0.0
  clss_He=0.0
  sdf=0.0
  !tgge = tgge + rkHet *dt                      ! gas generated
  !tggb = tggb + rgrain*dt                      ! gas lost at GB
  !tgdi = tgdi + (rdisl + rdisl_p)*dt/2.0d0     ! gas lost at disl'

  do k=5+2,neq
    i=idarrayi(k)
    j=idarrayj(k)
    sdf(i,j)=abs(y(k))
    clss_V=clss_V+Pvarr(i,j)*sdf(i,j)
    clss_i=clss_i+Qiarr(i,j)*sdf(i,j)
    clss_He=clss_He+PHearr(i,j)*sdf(i,j)
  end do
  clss_V=clss_V/Dv
  clss_i=clss_i/Di
  clss_He=clss_He/DHe


  return
end subroutine
