subroutine euler ( n,y,ydot,t )
!------------------------------------------------------------------
! Eulers method, simple ODE solver (very basic).
!------------------------------------------------------------------ 

      use nuclearvars

      implicit none
      external funcallngrp
      integer :: i, ii, iii,n,nmean,nround,iflagt
      real(kind=dp), dimension(n) :: y,ydot,ystar
      real(kind=dp) :: dtt,xmean,dely,t,utol,ytol,atbub
!-----------------------------------------------------------------------

      nround=10
      atbub=1.0E-12
      do i = 1,neq
         ystar(i) = y(i)
      enddo
      nsteps=1000!nint(tmax/dt)
      nmean = 0
      dtt   = dt/nsteps

      do iii = 1,nsteps

         t = t + dtt

         do ii = 1,nround

            call funcallngrp(neq,t,ystar,ydot )
            iflagt = 0

            do i = 1,neq
               dely     = y(i) + ydot(i)*dtt
               utol     = abs(ystar(i)-dely)
               ystar(i) = dely
               ytol     = rtol*dabs(ystar(i)) + atbub
               if (utol.gt.ytol) iflagt = 1
            enddo

            if (iflagt.eq.0) exit         

         enddo

         nmean = nmean + ii   

         do i = 1,neq
            y(i) = ystar(i)
         enddo

      enddo

      xmean = real(nmean)/real(nsteps)


      return
      end

!***********************************************************************

