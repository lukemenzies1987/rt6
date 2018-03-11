subroutine wtf2(y)
  use nuclearvars
  use twodimensionalvars
  implicit none

  integer :: n,m, nunit, i,j,k
  real(kind=dp), dimension(neq) :: y
  real(kind=dp) :: linf,Fedens
  real(kind=dp), dimension(1:xnumg,0:mnumg) :: sdf,sdf2,sdf3
  !integer, dimension(1:xnumg,0:mnumg) :: sdf,sdf2,sdf3

  nunit=23
  Fedens=8.5E+28
  k=6
  sdf=0!.0
  if (flaggrouping .eq. 1) then 
    do k=5+1,neq1

      i=idarrayi(k)
      j=idarrayj(k)
  !print*, j
      sdf(i,j)=y(k)
    end do 
    
    do k=neq1+1,neq,3
      i=idarrayi(k)
      j=idarrayj(k)
      sdf(i,j)=y(k)!*rgx(i)*rgm(j)
    end do 

  else 
    do k=5+1,neq

      i=idarrayi(k)
      j=idarrayj(k)
  !print*, j
      !if (k .eq. 6) print*, i,k,'ggg'
      sdf(i,j)=y(k)
      !sdf2(i,j)=j!
      !sdf3(i,j)=k
    end do 

  end if
  !sdf(1,0)=y(1)
  !open (unit=nunit, access='append',file='arraytest.dat')
  open (unit=nunit, access='append',file='arraytest.dat')

  write(nunit, *) 'SDF i'  ! this gives you the line break

  do j=0,mnumg
    do i=1,xnumg
      !print*,mnumg
      !if (j .eq. 0 .and. i .eq. 2) print*, Qvarr(i,j),sdf(i,j),'tt'
      !write(nunit, '(I4,X,I4,X,I4,3X)', advance='no') sdf(i,j),sdf2(i,j),sdf3(i,j)
      write(nunit, '(ES20.10E3,X)', advance='no') sdf(i,j)!
    end do
    write(nunit, *) ''  ! this gives you the line break
  end do
  write(nunit, *) ''  ! this gives you the line break


  write(nunit, *) ' '  ! this gives you the line break

  close(nunit)






end subroutine wtf2

subroutine wtf1d(y)
  use nuclearvars
  use onedimensionalvars
  implicit none

  integer :: n,m, nunit, i,j,k
  real(kind=dp), dimension(neq) :: y
  real(kind=dp) :: linf,Fedens
  real(kind=dp), dimension(1:xnumg) :: sdf

  nunit=23
  Fedens=8.5E+28

  sdf=0.0
  if (flaggrouping .eq. 1) then 
    do k=2,neq1

      i=idarrayi(k)
  !print*, j
      sdf(i)=y(k)
    end do 
    
    do k=neq1+1,neq2,2
      i=idarrayi(k)
      sdf(i)=y(k)*rgx(i)
    end do 
  else 
    do k=5+1,neq

      i=idarrayi(k)
  !print*, j
      sdf(i)=y(k)
    end do 

  end if
  sdf(1)=y(1)
  !open (unit=nunit, access='append',file='arraytest.dat')
  open (unit=nunit, access='append',file='arraytest.dat')

  write(nunit, *) 'SDF i'  ! this gives you the line break

  do i=1,xnumg
      !print*,mnumg
      !if (j .eq. 0 .and. i .eq. 2) print*, Qvarr(i,j),sdf(i,j),'tt'
    write(nunit, '(EN20.10,X,I3,X)') sdf(i)
  end do
  write(nunit, *) ''  ! this gives you the line break


  write(nunit, *) ' '  ! this gives you the line break

  close(nunit)






end subroutine wtf1d

