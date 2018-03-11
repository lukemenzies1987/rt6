subroutine compressibility (x,Z) 
!------------------------------------------------------------------
! Compressibility (Z) of a He for diven He/vacancy ratio (x)
!------------------------------------------------------------------ 


  use nuclearvars
  implicit real*8 (a-h,o-z)

  Z = 1.0d0


  ! Carnahan-Starling corrected by R.E. Stoller

     ag    = 33.96855525d0
     bg    = 1.87d0
     cg11  = 0.42d0
     dhenm = 0.3135d0 &
&          *(0.8542d0-0.03996d0*dlog((temp)/9.16d0))
     sigma = dhenm*1.0d-9
     y  = x*pi*sigma**3.0d0/6.0d0/vol
     if (y .gt. 1.0d+100) y=1.0D+100
     !print*, y,x,sigma,omega,(1.0d0+y+y**2.0-y**3.0)/(1.0d0-y)**3.0
     
     Z  = (1.0d0+y+y**2.0d0-y**3.0d0)/(1.0d0-y)**3.0d0
     if (y.gt.0.5) Z = ag*dlog(x*bg-cg11)


  return
  end subroutine

