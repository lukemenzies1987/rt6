subroutine tolchangengrp(y)
  use nuclearvars
  implicit none
  integer :: srt, arrcount, i,j,k,meani,meanj
  real(kind=dp) :: xmeanmax,rtot,rtotD,tot,ntol,meanr
  real(kind=dp) , dimension(neq) :: y
  real(kind=dp), dimension(1:xnumg,0:mnumg) :: ftol
  atol(1)=Cv*rtol
  atol(2)=Ci*rtol
  atol(3)=CHe*rtol
  atol(4)=CGb*rtol
  atol(5)=Cg*rtol
  ntol=10.
  srt=6
  rtot=0.0
  tot=0.0
  ftol=0.0
  !first half--------------------------------------
  do k=5+1,neq
    i=idarrayi(k)
    j=idarrayj(k)
    ftol(i,j)=abs(y(k))
    rtot=rtot+r(i)*ftol(i,j)*1.E+09
    tot=tot+ftol(i,j)
    !if (y(k) .eq. 0.0) print*, k,i,j
    !nucdens=nucdens+abs(y(k))*Fedens
  end do


  xmeanmax=0.0
  meanr = rtot/tot
  meani=1
  do i=1,xnumg
    if (meanr .lt. r(i)*1.E+09) then 
      meani=i 
      exit
    end if 
!print*, r(i)*1.E+09,'ri',meani,meanr
  end do  

  do j=0,mnumg
    if (xmeanmax .lt. abs(ftol(meani,j)) ) then 
      xmeanmax = abs(ftol(meani,j))
      meanj=j
    end if
   ! print*, abs(ftol(meani,j)),'fdsd',j
  end do 
 ! print*, rtot,tot,meanr,xmeanmax,meani,meanj,'lk'
  atol(srt:neq)=xmeanmax*rtol

  !print*, atol(srt:neq)

  return
end subroutine 

subroutine tolchangengrptest(y)
  use nuclearvars

  implicit none
  integer :: srt, arrcount, i,j,k,meani,meanj
  real(kind=dp) :: xmeanmax,rtot,rtotD,tot,ntol,meanr,limit
  real(kind=dp) , dimension(neq) :: y
  real(kind=dp), dimension(1:xnumg,0:mnumg) :: ftol
  atol(1)=Cv*rtol
  atol(2)=Ci*rtol
  atol(3)=CHe*rtol
  atol(4)=CGb*rtol
  atol(5)=Cg*rtol
  ntol=10.
  srt=6
  rtot=0.0
  tot=0.0
  ftol=0.0
  limit=vol!1.0E-40_dp
  !first half--------------------------------------
  do k=5+1,neq
    i=idarrayi(k)
    j=idarrayj(k)
    ftol(i,j)=abs(y(k))
    !rtot=rtot+r(i)*ftol(i,j)*1.E+09
    !tot=tot+ftol(i,j)
    !if (y(k) .eq. 0.0) print*, k,i,j
    !nucdens=nucdens+abs(y(k))*Fedens
    if(ftol(i,j)*rtol .le. limit) then 
      atol(k)=limit
    else
      atol(k)=ftol(i,j)*rtol
    end if
  end do



  !xmeanmax=0.0
  !meanr = rtot/tot
  !meani=1
  !do i=1,xnumg
  !  if (meanr .lt. r(i)*1.E+09) then 
  !    meani=i 
  !    exit
  !  end if 
!print*, r(i)*1.E+09,'ri',meani,meanr
  !end do  

  !do j=0,mnumg
  !  if (xmeanmax .lt. abs(ftol(meani,j)) ) then 
  !    xmeanmax = abs(ftol(meani,j))
  !    meanj=j
  !  end if
   ! print*, abs(ftol(meani,j)),'fdsd',j
  !end do 
 ! print*, rtot,tot,meanr,xmeanmax,meani,meanj,'lk'
  !atol(srt:neq)=xmeanmax*rtol

  !print*, atol(srt:neq)

  return
end subroutine 
