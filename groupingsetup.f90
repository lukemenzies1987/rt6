subroutine groupingsetup
!------------------------------------------------------------------
! The grouping setup for the 2 dimensional case. 
!------------------------------------------------------------------    
  use nuclearvars
  use twodimensionalvars
  implicit none

  integer :: i,j,k,nox,nom,arrcount
  real(kind=dp):: x,m,dx,dm

  interface 
    function mean(xi,dx)
      use kindmod
      implicit none
      real(kind=dp) :: xi,dx
      real(kind=dp) :: mean
    end function

    function sigmasq(dx)
      use kindmod
      implicit none
      real(kind=dp) :: dx
      real(kind=dp) :: sigmasq
    end function
  end interface

  if ((gxw .gt. xnum-srtgrpx) .or. (gmw .gt. mnum-srtgrpm)) then 
    print*,"Group widths are too large for maximum size"
    stop
  end if 
!  if (srtgrpm .le. m0) then 
!    print*,"Initial grouping starts before m0"
!    stop
!  end if 


  !noofg = noofg + unevengx*(nint((xnum-srtgrpx)/gx)-xnum)+unevengm*(nint((mnum-srtgrpm)/gm)-mnum -1)+ srtgrpx*unevengm


 

  k=0
  x=srtgrpx
  do while(x .lt. xnum)
    k=k+1
    
    dx= nint( ((36.*Pi/vol)**(1._dp/3_dp)*(x-srtgrpx)**(2._dp/3._dp)*dr)) +1
    !if (k .lt. 30) dx=1
    if (x+dx .gt. xnum) then 
      x=xnum
    else 
      x=x+dx
      
    end if
   ! print*, x ,'x',k
  end do 
  xnumg=k+srtgrpx

  k=0
  m=srtgrpm
  do while(m .lt. mnum)
    k=k+1
    dm= nint( (36.*Pi/vol)**(1._dp/3_dp)*(m-srtgrpm)**(2._dp/3._dp)*drm ) +1
    if (m+dm .gt. mnum) then 
      m=mnum
    else 
      m=m+dm
    end if
    !print*, m ,'m',k
  end do 
  mnumg=k+srtgrpm
  
  allocate(gx(0:xnumg+2),gm(-1:mnumg+1),rgx(0:xnumg+1),rgm(-1:mnumg+1),g1mi(-1:mnumg+1),& 
& g1xi(0:xnumg+1),rg1mi(-1:mnumg+1),rg1xi(0:xnumg+1))

  gx=0
  gm=0
  rgx=0.0
  rgm=0.0
  g1mi=0.0_dp
  g1xi=0.0_dp
  rg1xi=0.0_dp
  rg1mi=0.0_dp

  gx=gxw
  gm=gmw
  do i=1,srtgrpx
    gx(i)=1
    rgx(i)=real(gx(i),kind=dp)
  end do

  do j=0,srtgrpm
    gm(j)=1
    rgm(j)=real(gm(j),kind=dp)
  end do  
  x=srtgrpx
  do i=srtgrpx+1,xnumg
    
    dx= nint( (36.*Pi/vol)**(1._dp/3_dp)*(x-srtgrpx)**(2._dp/3._dp)*dr ) +1
    !if (i .lt. srtgrpx+1+ 30) dx=1
    !if (i .eq. srtgrpx+1) dx=1
    !if (i .eq. srtgrpx+2) dx=2
    if (x+dx .gt. xnum) then
      gx(i)= xnum-x
      rgx(i)=real(gx(i),kind=dp)
    else 
      gx(i)=dx
      rgx(i)=real(gx(i),kind=dp)

    end if
    x=x+gx(i)

  end do 
  m=srtgrpm
  do j=srtgrpm+1,mnumg
    dm= nint( (36.*Pi/vol)**(1._dp/3_dp)*(m-srtgrpm)**(2._dp/3._dp)*drm ) +1
    if (m+dm .gt. mnum) then
      gm(j)= mnum-m
      rgm(j)=real(gm(j),kind=dp)
    else 
      gm(j)=dm
      rgm(j)=real(gm(j),kind=dp)
    end if
    m=m+gm(j)

  end do 


  nox=0
  do i=1,xnumg
    g1xi(i)= nox+gx(i)
    !g1xf(i)= nox+gx(i)
    nox=nox+gx(i)
    rg1xi(i)=real(g1xi(i),kind=dp)

  end do
  nom=-1
  do j=0,mnumg
    g1mi(j)= nom+gm(j)

    nom=nom+gm(j)
    rg1mi(j)=real(g1mi(j),kind=dp)

  end do
  
  arrcount=5+1
!first half--------------------------------------
  do j=0,srtgrpm
    do i=1,srtgrpx     
      if (rg1mi(j) .gt. g00 + g11*(rg1xi(i)-1._dp)) then 
        cycle
      end if
      !write(77,*) arrcount,'arrcount',i,j,1

      arrcount=arrcount+1

    end do
  end do 
  neq1=arrcount-1

  do j=0,srtgrpm
    do i=srtgrpx+1,xnumg
      if (rg1mi(j) .gt. g00 + g11*(rg1xi(i)-1._dp)) then 
        cycle
      end if
      !write(77,*) arrcount,'arrcount',i,j,2
      arrcount=arrcount+3

    end do
  end do 
  neq2= arrcount-1
!second half--------------------------------------
  do j=srtgrpm+1,mnumg
    do i=1,xnumg
      if (rg1mi(j) .gt. g00 + g11*(rg1xi(i)-1._dp)) then 
        cycle
      end if
      !write(77,*) arrcount,'arrcount',i,j,3
      arrcount=arrcount+3

    end do    
  end do 

  neq= arrcount-1


  allocate(itermx(1:xnumg),itermm(0:mnumg),zerox(1:xnumg),zerom(0:mnumg),& 
& meanxarr(0:xnumg+1),meanmarr(-1:mnumg+1), &
& L0(0:xnumg+2,-1:mnumg+1),L1x(0:xnumg+1,-1:mnumg+1),L1m(0:xnumg+1,-1:mnumg+1), &
 idarrayi(neq),idarrayj(neq))



  meanxarr=0.0_dp
  meanmarr=0.0_dp


  do i=1,xnumg
    meanxarr(i)= mean(rg1xi(i),rgx(i))
  end do

  do j=0,mnumg
    meanmarr(j)= mean(rg1mi(j),rgm(j))
  end do

  !do i=1,xnumg
  !  rg1xi(i)=meanxarr(i)!= mean(rg1xi(i),rgx(i))
  !end do

  !do j=0,mnumg
  !  rg1mi(j)=meanmarr(j)!= mean(rg1mi(j),rgm(j))
  !end do

  do i=1,xnumg
    if (gx(i)==1) then 
      itermx(i) =0.0_dp
      zerox(i) =0.0_dp
      cycle
    end if
    itermx(i)=(rgx(i)-1._dp)/(2._dp*sigmasq(rgx(i))*rgx(i))
    zerox(i) =1.0_dp

  end do

  do j=0,mnumg
    if (gm(j)==1) then 
      itermm(j) =0.0_dp
      zerom(j) =0.0_dp
      cycle
    end if
    itermm(j)=(rgm(j)-1._dp)/(2._dp*sigmasq(rgm(j))*rgm(j))
    zerom(j) =1.0_dp
  end do
  arrcount=5+1


!first half--------------------------------------
  do j=0,srtgrpm
    do i=1,srtgrpx
      if (rg1mi(j) .gt. g00 + g11*(rg1xi(i)-1._dp)) then 
        !arrcount=arrcount+1
        cycle
      end if   
      if (i .eq. 1 .and. j.eq. 1) minnum=arrcount       
      idarrayi(arrcount)=i
      idarrayj(arrcount)=j
      arrcount=arrcount+1

    end do
  end do 
  neq1=arrcount-1
  do j=0,srtgrpm
    do i=srtgrpx+1,xnumg
      if (rg1mi(j) .gt. g00 + g11*(rg1xi(i)-1._dp)) then 
        cycle
      end if
      
      idarrayi(arrcount)=i
      idarrayj(arrcount)=j
      idarrayi(arrcount+1)=i
      idarrayj(arrcount+1)=j
      idarrayi(arrcount+2)=i
      idarrayj(arrcount+2)=j
      arrcount=arrcount+3

    end do
  end do 

!second half--------------------------------------
  do j=srtgrpm+1,mnumg
    do i=1,xnumg
      if (rg1mi(j) .gt. g00 + g11*(rg1xi(i)-1._dp)) then 
        cycle
      end if
     
      idarrayi(arrcount)=i
      idarrayj(arrcount)=j
      idarrayi(arrcount+1)=i
      idarrayj(arrcount+1)=j
      idarrayi(arrcount+2)=i
      idarrayj(arrcount+2)=j
      arrcount=arrcount+3

    end do    
  end do 



  neqend= arrcount-1

  L0=0.0_dp
  L1x=0.0_dp
  L1m=0.0_dp
  itermm(0)=0.0
  itermx(1)=0.0
  return
end subroutine groupingsetup
