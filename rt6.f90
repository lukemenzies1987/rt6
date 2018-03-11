program rt6
!---------------------------------------------------------------------
!Rate Theory program, developed by Luke Menzies, University of 
!Manchester, United Kingdom.
!----------------------------------------------------------------------

  use nuclearvars
  implicit none
  
  interface 
    subroutine getparameters
    end subroutine

    subroutine clusters
    end subroutine

    subroutine deallocatearrays1d
    end subroutine

    subroutine deallocatearrays
    end subroutine

  end interface 
  call cpu_time(start)
  !Reads in parameter from param.in file. 
  call getparameters 
  !The main subroutine. 
  call clusters
  !Deallocates arrays
  if (flagVoids .eq. 1) call deallocatearrays1d
  if (flagClusters .eq. 1 ) call deallocatearrays
end program rt6



