subroutine tolchange(y)
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
  do k=5+1,neq1
    i=idarrayi(k)
    j=idarrayj(k)
    ftol(i,j)=abs(y(k))
    rtot=rtot+r(i)*ftol(i,j)
    tot=tot+ftol(i,j)
    !nucdens=nucdens+abs(y(k))*Fedens
  end do
  !second half--------------------------------------
  do k=neq1+1,neq,3   
    i=idarrayi(k)
    j=idarrayj(k) 
    ftol(i,j)=abs(y(k))
    rtot=rtot+r(i)*ftol(i,j)*rgx(i)*rgm(j)
    tot=tot+ftol(i,j)*rgx(i)*rgm(j)
    !nucdens=nucdens+abs(y(k))*gx(i)*gm(j)*Fedens
  end do 

  xmeanmax=0.0
  meanr = rtot/tot
  

  do i=1,xnumg
    if (meanr .lt. r(i)) then 
      meani=i 
      exit
    end if 
  end do  

  do j=0,mnumg
    if (xmeanmax .lt. abs(ftol(meani,j)) ) then 
      xmeanmax = abs(ftol(meani,j))
      meanj=j
    end if
  end do 
  ! xmeanmax= rtot/tot
!first half--------------------------------------
  atol(srt:neq1)=xmeanmax*rtol
  !do j=0,srtgrpm
  !  do i=1,srtgrpx
  !    atol(arrcount)=atolreq/ntol
      
!print*, arrcount, i,j
  !    arrcount=arrcount+1
  !  end do
  !end do 
  !neq1=arrcount-1
  !arrcount = neq1+1
!print*, neq1,'neq1'
  !do j=0,srtgrpm
  !  do i=srtgrpx+1,xnumg
  !    idarrayi(arrcount)=i
  !    
  !    idarrayi(arrcount+1)=i
  !    
  !    idarrayi(arrcount+2)=i
  !    
  !    arrcount=arrcount+3
  !  end do
  !end do 
  do arrcount=neq1+1,neq,3
    atol(arrcount)=xmeanmax*rtol
    atol(arrcount+1)=xmeanmax*rtol/ntol
    atol(arrcount+2)=xmeanmax*rtol/ntol
  end do  
!second half--------------------------------------
  !do j=srtgrpm+1,mnumg
  !  do i=1,xnumg
  !    idarrayi(arrcount)=i
  !    idarrayj(arrcount)=j
  !    idarrayi(arrcount+1)=i
  !    idarrayj(arrcount+1)=j
  !    idarrayi(arrcount+2)=i
  !    idarrayj(arrcount+2)=j
  !    arrcount=arrcount+3
  !  end do    
  !end do 
  !print*, xmeanmax*rtol,rtot,tot
  

  return
end subroutine 

subroutine tolchange1d(y)
  use nuclearvars
  implicit none
  integer :: srt, arrcount, i,j,k,meani,meanj
  real(kind=dp) :: xmeanmax,rtot,rtotD,tot,ntol,meanr
  real(kind=dp) , dimension(neq) :: y
  real(kind=dp), dimension(1:xnumg) :: ftol
  atol(1)=Cv*rtol
  atol(neq2+1)=Ci*rtol
  atol(neq2+2)=Cg*rtol
  ntol=10.
  srt=6
  rtot=0.0
  tot=0.0
  ftol=0.0
  !first half--------------------------------------
  do k=1,neq1
    i=idarrayi(k)

    ftol(i)=abs(y(k))
    rtot=rtot+r(i)*ftol(i)*1.E+09
    tot=tot+ftol(i)
    !nucdens=nucdens+abs(y(k))*Fedens
  end do
  !second half--------------------------------------
  do k=neq1+1,neq2,2   
    i=idarrayi(k)

    ftol(i)=abs(y(k))
    rtot=rtot+r(i)*ftol(i)*rgx(i)*1.E+09
    tot=tot+ftol(i)*rgx(i)
    !nucdens=nucdens+abs(y(k))*gx(i)*gm(j)*Fedens
  end do 

  xmeanmax=0.0
  meanr = rtot/tot
  

  do i=1,xnumg

    if (meanr .lt. r(i)*1.E+09) then 
      meani=i 
      exit
    end if 
  end do  


  xmeanmax = abs(ftol(meani))
print*, xmeanmax,'xmeanmax',meani
  ! xmeanmax= rtot/tot
!first half--------------------------------------
  atol(srt:neq1)=xmeanmax*rtol
  !do j=0,srtgrpm
  !  do i=1,srtgrpx
  !    atol(arrcount)=atolreq/ntol
      
!print*, arrcount, i,j
  !    arrcount=arrcount+1
  !  end do
  !end do 
  !neq1=arrcount-1
  !arrcount = neq1+1
!print*, neq1,'neq1'
  !do j=0,srtgrpm
  !  do i=srtgrpx+1,xnumg
  !    idarrayi(arrcount)=i
  !    
  !    idarrayi(arrcount+1)=i
  !    
  !    idarrayi(arrcount+2)=i
  !    
  !    arrcount=arrcount+3
  !  end do
  !end do 
  do arrcount=neq1+1,neq2,2
    atol(arrcount)=xmeanmax*rtol
    atol(arrcount+1)=xmeanmax*rtol/ntol
  end do  
!second half--------------------------------------
  !do j=srtgrpm+1,mnumg
  !  do i=1,xnumg
  !    idarrayi(arrcount)=i
  !    idarrayj(arrcount)=j
  !    idarrayi(arrcount+1)=i
  !    idarrayj(arrcount+1)=j
  !    idarrayi(arrcount+2)=i
  !    idarrayj(arrcount+2)=j
  !    arrcount=arrcount+3
  !  end do    
  !end do 
  !print*, xmeanmax*rtol,rtot,tot
  

  return
end subroutine 

subroutine tolchangetest(y)
  use nuclearvars
  implicit none
  integer :: srt, arrcount, i,j,k,meani,meanj
  real(kind=dp) :: xmeanmax,rtot,rtotD,tot,ntol,meanr,limit
  real(kind=dp) , dimension(neq) :: y
  real(kind=dp), dimension(1:xnumg,0:mnumg) :: ftol,L0tol,Lxtol,Lmtol
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
  limit=1.0E-40_dp
  do k=5+1,neq1
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


  do k=neq1+1,neq,3
    !i=idarrayi(k)
    !j=idarrayj(k)
    L0tol(i,j)=abs(y(k))
    Lxtol(i,j)=abs(y(k+1))
    Lmtol(i,j)=abs(y(k+2))

    if (L0tol(i,j)*rtol .lt. limit) then 
      atol(k)=limit
    else 
      atol(k)=abs(L0tol(i,j))*rtol
    end if 

    if (Lxtol(i,j)*rtol .lt. limit) then 
      atol(k)=limit
    else 
      atol(k)=abs(Lxtol(i,j))*rtol/ntol
    end if 

    if (Lmtol(i,j)*rtol .lt. limit) then 
      atol(k)=limit
    else 
      atol(k)=abs(Lmtol(i,j))*rtol/ntol
    end if 

  end do

  

  return
end subroutine



subroutine tolchange1dtest(y)
  use nuclearvars
  implicit none
  integer :: srt, arrcount, i,j,k,meani,meanj
  real(kind=dp) :: xmeanmax,rtot,rtotD,tot,ntol,meanr,limit
  real(kind=dp) , dimension(neq) :: y
  real(kind=dp), dimension(1:xnumg) :: ftol,ftolmom
  atol(1)=Cv*rtol
  atol(neq2+1)=Ci*rtol
  atol(neq2+2)=Cg*rtol
  ntol=10.
  srt=2
  rtot=0.0
  tot=0.0
  ftol=0.0
  limit=1.0E-40_dp
  !first half--------------------------------------
 ! do k=1,neq1
  !  i=idarrayi(k)

  !  ftol(i)=abs(y(k))
    !rtot=rtot+r(i)*ftol(i)*1.E+09
    !tot=tot+ftol(i)
    !nucdens=nucdens+abs(y(k))*Fedens
 ! end do
  !second half--------------------------------------
 ! do k=neq1+1,neq2,2   
 !   i=idarrayi(k)

  !  ftol(i)=abs(y(k))
  !  ftolmom(i)=abs(y(k+1))
    !rtot=rtot+r(i)*ftol(i)*rgx(i)*1.E+09
    !tot=tot+ftol(i)*rgx(i)
    !nucdens=nucdens+abs(y(k))*gx(i)*gm(j)*Fedens
  !end do 


  do arrcount=srt,neq1
    i=idarrayi(arrcount)
    ftol(i)=abs(y(arrcount))
    if (ftol(i)*rtol .lt. limit) then 
      atol(arrcount)=limit
    else 
      atol(arrcount)=abs(ftol(i))*rtol
    end if 
    !print*, atol(arrcount) ,'atol'
  end do  

  do arrcount=neq1+1,neq2,2
    i=idarrayi(arrcount)
    ftol(i)=abs(y(arrcount))
    ftolmom(i)=abs(y(arrcount+1))
    if (ftol(i)*rtol .lt. limit) then 
      atol(arrcount)=limit
      atol(arrcount+1)=limit
    else 
      atol(arrcount)=abs(ftol(i))*rtol
      atol(arrcount+1)=abs(ftolmom(i))*rtol/ntol
    !  print*, atol(arrcount) ,'atol'
    end if 
  end do  

  

  return
end subroutine 
