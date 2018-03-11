subroutine wtffv(y)
  use nuclearvars
  implicit none

  integer :: n,m, nunit, i,j,k,g
  real(kind=dp), dimension(neq) :: y
  real(kind=dp) :: linf,Fedenscorr
  real(kind=dp), dimension(1:xnumg,0:mnumg) :: sdf

  nunit=23
  Fedenscorr=8.5E+28!/1.E06
  k=6
  
  !first half--------------------------------------
  do k=5+1,neq1

    i=idarrayi(k)
    j=idarrayj(k)
    sdf(i,j)=y(k)
  end do 

  do k=neq1+1,neq,3
    i=idarrayi(k)
    j=idarrayj(k)
    sdf(i,j)=y(k)

  end do  
  sdf(1,0)=y(1)
  !open (unit=nunit, access='append',file='arraytest.dat')
  open (unit=nunit,file='sdfforvoids.dat')
  j=1
  write(nunit, *) 'SDF'  ! this gives you the line break
  do i=1,srtgrpx
    write(nunit, '(2(EN20.10,X))') real(j), abs(sdf(i,0))!*Fedenscorr
    !print*,sdf(i,0)
    j=j+1
  end do 
  g=srtgrpx
  do i=srtgrpx,xnumg
    do k=1,gx(i)
      !print*, k,i,L0(i,0)+ L1x(i,0)*(k-meanxarr(i))
      write(nunit, '(2(EN20.10E3,X))') real(j),abs( abs(L0(i,0))+ L1x(i,0)*(g+k-meanxarr(i)))!*Fedenscorr
      !print*, L0(i,0),L1x(i,0)
      !print*, L0(i,0),L1x(i,0)*(g+k-meanxarr(i)),i,g+k
      j=j+1
    end do 
    g=g+gx(i)
  end do 

  
  close(nunit)


!4.0d0*pi*bvv(i)*(rwh(i))**2*1.0d-27/omega**2



end subroutine wtffv

subroutine wtffvalt(y)
  use nuclearvars
  implicit none

  integer :: n,m, nunit, i,j,k,g
  real(kind=dp), dimension(neq) :: y
  real(kind=dp) :: linf,Fedenscorr
  !real(kind=dp), dimension(1:xnumg,0:mnumg) :: sdf

  nunit=23
  Fedenscorr=8.5E+28/1.E06
  
  do j=0,mnumg
    do i = 1,xnumg
      L0(i,j)=y(3*i-2+5+3*j*xnumg)
      L1x(i,j)=y(3*i-1+5+3*j*xnumg)
      L1m(i,j)=y(3*i+5+3*j*xnumg) 
    end do 
  end do 
  !open (unit=nunit, access='append',file='arraytest.dat')
  open (unit=nunit,file='sdfforvoids.dat')
  j=1
  write(nunit, *) 'SDF'  ! this gives you the line break

  do i=1,srtgrpx
    write(nunit, '(2(EN20.10,X))') real(j), abs( L0(i,0)+ L1x(i,0)*(i-meanxarr(i)))!*Fedenscorr
    !print*,sdf(i,0)
  end do 

  g=srtgrpx
  
  do i=srtgrpx+1,xnumg
    do k=1,gx(i)
      !print*, k,i,L0(i,0)+ L1x(i,0)*(k-meanxarr(i))
      write(nunit, '(2(EN20.10,X))') real(j),( abs(L0(i,0))+ (L1x(i,0))*(g+k-meanxarr(i)))!*Fedenscorr
      j=j+1
      !print*, g+k,'FG'
print*, L1x(i,0)*(g+k-meanxarr(i)),i,k,(L1x(i,0)),'lx',abs(L0(i,0)),'l0',(L1x(i,0))*(g+k-meanxarr(i)),(g+k-meanxarr(i))
print*,( abs(L0(i,0))+ (L1x(i,0))*(g+k-meanxarr(i))),'test'
    end do 
    g=g+gx(i)
  end do 


  
  close(nunit)






end subroutine wtffvalt
