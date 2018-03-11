subroutine wtf(array,arrays,n,m,sf)
  implicit none
  integer, parameter :: dp=selected_real_kind(15,300)
  integer :: n,m, nunit, i,j 
  real(kind=dp), dimension(1:n,0:m) :: array
  !real(kind=dp), dimension(n*m+n+3) :: arrays
  real(kind=dp), dimension(0:n+1) :: arrays
  character :: sf
  nunit=23
  !open (unit=nunit, access='append',file='arraytest.dat')
  open (unit=nunit, access='append',file='arraytest.dat')
  if (sf=='d') then
    write(nunit, *) 'f'  ! this gives you the line break

    do j=0,m
      do i=1,n
        write(nunit, '(EN20.10,X)', advance='no') array(i,j)
      end do
      write(nunit, *) ''  ! this gives you the line break
    end do
    write(nunit, *) ''  ! this gives you the line break
  end if

  if (sf=='s') then
    write(nunit, *) 'ydot'  ! this gives you the line break
    !do i=1,n*(m+1)+3
    do i=1,n
      write(nunit, '(EN20.10,X)', advance='no') arrays(i)
    end do
  end if
  write(nunit, *) ' '  ! this gives you the line break
  close(nunit)

  

end subroutine wtf



