subroutine groupingsetup1d
!------------------------------------------------------------------
! The grouping setup for the 1 dimensional case. 
!------------------------------------------------------------------    
  use nuclearvars
  use onedimensionalvars
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
    print*, x ,'x',k
  end do 
  xnumg=k+srtgrpx

  k=0
  m=srtgrpi
  do while(m .lt. inum)
    k=k+1
    dm= nint( (4.*Pi*b/vol)**(1._dp/2_dp)*(m-srtgrpi)**(1._dp/2._dp)*dri ) +1
    if (m+dm .gt. inum) then 
      m=inum
    else 
      m=m+dm
    end if
    print*, m ,'m',k
  end do 


  inumg=k+srtgrpi
  print*,xnumg,'xnumg',inumg

  arrcount=1
!first half--------------------------------------
  do i=1,srtgrpx
    arrcount=arrcount+1
  end do
  neq1=arrcount-1
!print*, neq1,'neq1
  do i=srtgrpx+1,xnumg
    arrcount=arrcount+2
  end do
  neq2=arrcount-1

!second half--------------------------------------
  !do i=1,srtgrpi
  !  arrcount=arrcount+1
  !end do

  !neq3=arrcount-1
!print*, neq1,'neq1
 ! do i=srtgrpi+1,inumg
 !   arrcount=arrcount+2
  !end do

  !neq4=arrcount-1
  neq=neq2+2
  allocate(gx(0:xnumg+1),gi(0:inumg+1),g1ii(0:inumg+1),g1xi(0:xnumg+1),&!g1mf(-1:inumg),g1xf(0:xnumg),&
& itermx(1:xnumg),itermi(0:inumg),meanxarr(0:xnumg+1),meaniarr(0:inumg+1), &
& L0x(0:xnumg+1),L0i(0:inumg+1),L1x(0:xnumg+1),L1i(0:inumg+1), &
rg1ii(0:inumg+1),rg1xi(0:xnumg+1),rgx(0:xnumg+1),rgi(0:inumg+1),itermg1(0:xnumg),idarrayi(neq2))!,rg1mf(-1:inumg),rg1xf(0:xnumg))

  
  Xg=1

  gx=0
  gi=0
  g1ii=0.0_dp
  g1xi=0.0_dp
  rg1xi=0.0_dp
  rg1ii=0.0_dp
  meanxarr=0.0_dp
  meaniarr=0.0_dp
  !g1mf=0.0_dp
  !g1xf=0.0_dp
  

    
   

  !gx=gxw
  !gi=giw
  do i=1,srtgrpx
    gx(i)=1
  end do

  do j=1,srtgrpi
    gi(j)=1
  end do  
  x=srtgrpx
  do i=srtgrpx+1,xnumg
    
    dx= nint( (36.*Pi/vol)**(1._dp/3_dp)*(x-srtgrpx)**(2._dp/3._dp)*dr ) +1

    if (x+dx .gt. xnum) then
      gx(i)= xnum-x
    else 
      gx(i)=dx
      !print*, gx(i),i
    end if
    x=x+gx(i)
    print*,i,x,dx,'dx',gx(i)
  end do 
  m=srtgrpi
  do j=srtgrpi+1,inumg
    dm= nint( (4.*Pi*b/vol)**(1._dp/2_dp)*(m-srtgrpi)**(1._dp/2._dp)*dri ) +1
    if (m+dm .gt. inum) then
      gi(j)= inum-m
    else 
      gi(j)=dm
    end if

    m=m+gi(j)
    !print*, m,'m',dm,inum
  end do 
  
  rgi=real(gi)
  rgx=real(gx)
  nox=0
  do i=1,xnumg
    g1xi(i)= nox+gx(i)
    nox=nox+gx(i)
    rg1xi(i)=real(g1xi(i))
    print*,gx(i),rgx(i),'g'
    
  end do
  !do i=srtgrpx,xnumg
  !  write(32,*) rg1xi(i),rgx(i),i
  !end do
  nom=0
  do j=1,inumg
    g1ii(j)= nom+gi(j)
    nom=nom+gi(j)
    rg1ii(j)=real(g1ii(j))
    print*, rg1ii(j),'rg'
  end do
  
  do i=1,xnumg
    meanxarr(i)= mean(rg1xi(i),rgx(i))
!print*,meanxarr(i),i,xnumg,'ptd'
  end do

  do j=1,inumg
    meaniarr(j)= mean(g1ii(j),gi(j))
print*,meaniarr(j),gi(j),g1ii(j),'ptc'
  end do

  do i=1,xnumg
    if (gx(i)==1) then 
      itermx(i) =0.0_dp
    !  itermg1(i)=0.0_dp
      cycle
    end if
    itermx(i)=(rgx(i)-1._dp)/(2._dp*sigmasq(rgx(i))*rgx(i))
    itermg1(i)=1._dp/(sigmasq(rgx(i))*rgx(i))
  end do

  !do j=1,inumg
 !   if (gi(j)==1) then 
 !     itermi(j) =0.0_dp
 !     cycle
 !   end if
 !   itermi(j)=(rgi(j)-1._dp)/(2._dp*sigmasq(rgi(j))*rgi(j))
 ! end do
  arrcount=1

  !call glissgrouping

  do i=1,srtgrpx      
    idarrayi(arrcount)=i
    arrcount=arrcount+1
  end do 

  do i=srtgrpx+1,xnumg
    idarrayi(arrcount)=i
    idarrayi(arrcount+1)=i
    arrcount=arrcount+2
  end do

  !do i=1,srtgrpi      
  !  idarrayi(arrcount)=i
  !  arrcount=arrcount+1
  !end do 


  !do i=srtgrpi+1,inumg
  !  idarrayi(arrcount)=i
  !  idarrayi(arrcount+1)=i
  !  arrcount=arrcount+2
  !end do

  
  L0x=0.0_dp
  L0i=0.0_dp
  L1x=0.0_dp
  L1i=0.0_dp
  !itermi(1)=0.0
  itermx(1)=0.0
  return
end subroutine groupingsetup1d
