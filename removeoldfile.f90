subroutine dof()
  implicit none

  integer :: nunit,stat

  nunit=23

  open(unit=nunit, iostat=stat, file='arraytest.dat', status='old')
  if (stat == 0) close(nunit, status='delete')
  
  close(nunit)
end subroutine dof
