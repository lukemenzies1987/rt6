subroutine funcallinit(n, t, y, fd, ires)
!------------------------------------------------------------------
! Initial rate equation function
!------------------------------------------------------------------  
  use nuclearvars
  use twodimensionalvars
  implicit none
  external speciesgenerationHe
  
  real(kind=dp) :: siaabsorption, intabsorption, vacabsorption, Heemission, &
&  Heabsorption, vacemission, Hemonomeremission, Hevacabsorption,ss3v,ss3i, &
& sgbv,sgbi,sumf
  real(kind=nag_wp), intent(in) :: t
  integer , intent(in) :: n
  integer, intent(inout) :: ires
  integer :: i,j,k
  real(kind=nag_wp), intent(in) :: y(n)
  real(kind=nag_wp), intent(out) :: fd(n)
  real(kind=dp) :: invgx,of,f1,invgm,jxi,jxf,jmi,jmf,jx,jm, &
& f10,f11,f20,Gv2,Gv3,Gv4,Jaltxs,Jaltms,Jaltxb,Jaltmb,dL0,dLx,dLm
  real(kind=nag_wp) :: edge, rdislc,rdisle,recombine,P0,ss3He,sgbHe

  interface
    function gbss( sss,rgr )
      use kindmod
      implicit none
      real(kind=dp) :: sss,rgr
      real(kind=dp) :: gbss
    end function


    function Jaltm(i,j,x,m,ax,am)
      use nuclearvars
      use twodimensionalvars
      implicit none
      integer :: i,j
      real(kind=dp), dimension (0:xnumg+1) :: x
      real(kind=dp), dimension (-1:mnumg+1) :: m
      real(kind=dp) :: Jaltm,ax,am
    end function

    function Jaltx(i,j,x,m,ax,am)
      use nuclearvars
      use twodimensionalvars
      implicit none
      integer :: i,j
      real(kind=dp), dimension (0:xnumg+1) :: x
      real(kind=dp), dimension (-1:mnumg+1) :: m
      real(kind=dp) :: Jaltx,ax,am
    end function

  end interface

  siaabsorption=0.0_dp
  intabsorption=0.0_dp
  vacabsorption=0.0_dp
  Heemission=0.0_dp
  Heabsorption=0.0_dp
  vacemission=0.0_dp
  Hemonomeremission=0.0_dp
  Hevacabsorption=0.0_dp
  sigmavNv=0.0_dp
  f=0.0_dp
  L0=0.0
  L1x=0.0
  L1m=0.0

  do k=5+1,neq1
    i=idarrayi(k)
    j=idarrayj(k)
    f(i,j)=abs(y(k))
  end do 

  do k=neq1+1,neq,3
    i=idarrayi(k)
    j=idarrayj(k)
    L0(i,j)=abs(y(k))
    L1x(i,j)=y(k+1)
    L1m(i,j)=y(k+2)
  end do 

  Cv=y(1)
  Ci=y(2)
  CHe=y(3)
  CGb=y(4)
  Cg=y(5)
  f(1,0)=Cv
  do j=0,srtgrpm
    f(srtgrpx+1,j)=L0(srtgrpx+1,j)+L1x(srtgrpx+1,j)*(rg1xi(srtgrpx)+1.-meanxarr(srtgrpx+1))
    L0(srtgrpx,j)=f(srtgrpx,j)
    L1x(srtgrpx,j)=0.0_dp
    L1m(srtgrpx,j)=0.0_dp
  end do 

  do i=1,srtgrpx
    f(i,srtgrpm+1)=L0(i,srtgrpm+1)+L1m(i,srtgrpm+1)*(rg1mi(srtgrpm)+1.-meanmarr(srtgrpm+1))
    L0(i,srtgrpm)=f(i,srtgrpm)
    L1x(i,srtgrpm)=0.0_dp
    L1m(i,srtgrpm)=0.0_dp
  end do 

  !si_dis=rhod*Zi
  !sv_dis=rhod*Zv
  !sv_gb=ksqgb
  !si_gb=ksqgb
  !Gv2=4.39118365E-08
  !Gv3=3.11014284636941E-09
  !Gv4=1.49513852784706E-10

  call speciesgenerationHe(GHe,t)
  kgCg=2._dp*Cg*kgsq
  !k=6

!first half--------------------------------------
  do k=5+1,neq1

    i=idarrayi(k)
    j=idarrayj(k)
!print*, k,i,j
    dL0= ( Pvarr(i-1,j)*Cv*f(i-1,j)-(Qiarr(i,j)*Ci+Qvarr(i,j)+Qg(i,j)**kgcg)*f(i,j)-&
& Pvarr(i,j)*Cv*f(i,j)+(Qiarr(i+1,j)*Ci+Qvarr(i+1,j)+Qg(i+1,j)**kgcg)*f(i+1,j)) &
& +(PHearr(i,j-1)*CHe*f(i,j-1)-Qm(i,j)*f(i,j)-PHearr(i,j)*CHe*f(i,j)+Qm(i,j+1)*f(i,j+1))


    siaabsorption=siaabsorption+Qiarr(i,j)*Ci*f(i,j)
    vacabsorption=vacabsorption+Pvarr(i,j)*Cv*f(i,j)
    Heemission=Heemission+Qm(i,j)*f(i,j)
    Heabsorption=Heabsorption+PHearr(i,j)*CHe*f(i,j)
    vacemission=vacemission+Qvarr(i,j)*f(i,j)
		sigmavNv= sigmavNv + Qg(i,j)*f(i,j)


    fd(k)=dL0
  end do
  
!second half--------------------------------------

  do k=neq1+1,neq,3
    i=idarrayi(k)
    j=idarrayj(k)

    f1=L0(1,j)*rgm(j)
    invgm=(1._dp/rgm(j))
    invgx=(1._dp/rgx(i))

    jxi=Jaltx(i-1,j,rg1xi,meanmarr,0._dp,0._dp)
    jxf=Jaltx(i,j,rg1xi,meanmarr,0._dp,0._dp)
    jmi=Jaltm(i,j-1,meanxarr,rg1mi,0._dp,0._dp)
    jmf=Jaltm(i,j,meanxarr,rg1mi,0._dp,0._dp)



    of= L0(i,j)*rgx(i)*rgm(j)

    dL0 = invgx*(Jxi-Jxf) + invgm*(Jmi-Jmf)
    
    dLx = -itermx(i)*(Jxi+Jxf- 2._dp*Jaltx(i,j,meanxarr,meanmarr,-0.5_dp,0._dp)) + &
& invgm*( Jaltm(i,j-1,meanxarr,rg1mi,1._dp,0._dp)-Jaltm(i,j-1,meanxarr,rg1mi,0._dp,0._dp) ) - &
& invgm*(Jaltm(i,j,meanxarr,rg1mi,1._dp,0._dp)-jmf  )

    dLm =-itermm(j)*(Jmi+Jmf- 2._dp*Jaltm(i,j,meanxarr,meanmarr,0._dp,-0.5_dp) ) + &
& invgx*( Jaltx(i-1,j,rg1xi,meanmarr,0._dp,1._dp)-Jaltx(i-1,j,rg1xi,meanmarr,0._dp,0._dp) ) - &
& invgx*(Jaltx(i,j,rg1xi,meanmarr,0._dp,1._dp)-jxf  )

    

    siaabsorption=siaabsorption+Qiarr(i,j)*Ci*of
    vacabsorption=vacabsorption+Pvarr(i,j)*Cv*of
    Heemission=Heemission+Qm(i,j)*of
    Heabsorption=Heabsorption+PHearr(i,j)*CHe*of
    vacemission=vacemission+Qvarr(i,j)*of
		sigmavNv= sigmavNv + Qg(i,j)*of
      
    fd(k)=dL0
    fd(k+1)=dLx
    fd(k+2)=dLm

  end do    
  sigmavNv=sigmavNV/Dg   


  f10=f(1,0)
  f11=f(1,1)
  f20=f(2,0)
 
  ss3v =sv_dis
  ss3i=si_dis 
  disl=DHe*CHe*ZHe*rhod
  

 
  recombine=nui*Di*Cv*Ci 
  P0=2.*Pvarr(1,0)*f10*Cv
  vacemission=vacemission+Qm(1,1)*f11+ (Qiarr(2,0)*Ci+Qvarr(2,0)+Qg(2,0)*kgcg)*f20 !&
!& + Qiarr(1,1)*f11*Ci
  Heabsorption=Heabsorption+PHearr(1,0)*f10*CHe
  Heemission=Heemission+Qiarr(1,1)*f11*Ci+ Qg(1,1)*kgCg*f11

  if(Cv .ne. 0.0_dp) then 
    ss3v =Zv*rhod +(vacabsorption+P0)/(Dv*Cv)+PHearr(1,0)*f10/DHe!sv_dis +(vacabsorption+P0)/(Dv*Cv)+fPHearr(1,0)*f10/DHe!*CHe 
  end if 

  P0=P0+PHearr(1,0)*CHe*f10
  if (Ci .ne. 0.0) then
    ss3i=Zi*rhod + siaabsorption/(Di*Ci)!si_dis + siaabsorption/(Di*Ci)
  end if
  !ss3v =ss3v+(vacabsorption+PHearr(1,0)*f10*CHe+Pvarr(1,0)*f10*Cv)/(Dv*Cv)
  !ss3i=ss3i + siaabsorption/(Di*Ci)
  tgdid=tgdi!/(rhod*vol**(2._dp/3._dp))

  !rdisle = DHe*CHe*rhod*tgdid*dexp(-edge/(kB*temp)) 
  !disl=ZHe*DHe*CHe*rhod*(disl0/(1.E+04*tgdid+disl0))!(1._dp-tgdid)

  !if (disl .le. 0.0) disl=0.0
  if (CHe .ne. 0.0) then 
    ss3He=(disl+Heabsorption)/(DHe*CHe)
  end if
  

  sgbv = gbss(ss3v,l)
  sv_gb=sgbv


  sgbi = gbss(ss3i,l)
  si_gb= sgbi

  sgbHe=gbss(ss3He,l)
  sHe_gb=sgbHe
 



  dCv = Gv  -recombine - Dv*(Cv-Cv0)*sv_dis  -Dv*(Cv-Cv0)*sv_gb & 
& - P0 - vacabsorption + vacemission

  dCi = Geni - recombine -Di*Ci*si_dis -siaabsorption  - Di*Ci*si_gb

  dCHe = GHe  + Heemission - disl &! disl 
& -Heabsorption  -DHe*CHe*sHe_gb

  dCGb= DHe*CHe*sHe_gb

  kgsq=(fstterm+secndterm+sigmavNv)
  dCg= Gg-Dg*Cg*2._dp*kgsq*kgsq

  do i=2, maxgv
    fd(5+i)=fd(5+i)+genv(i)
  end do 
  
  fd(1)=dCv
  fd(2)=dCi
  fd(3)=dCHe
  fd(4)=dCGb
  fd(5)=dCg
  fd(6)=dCv

  return
end subroutine


subroutine funcall(n, t, y, fd, ires)
!------------------------------------------------------------------
! Rate equation function, grouping method
!------------------------------------------------------------------  
  use nuclearvars
  use twodimensionalvars
  implicit none
  external speciesgenerationHe
  
  real(kind=dp) :: siaabsorption, intabsorption, vacabsorption, Heemission, &
&  Heabsorption, vacemission, Hemonomeremission, Hevacabsorption,ss3v,ss3i, &
& sgbv,sgbi,sumf
  real(kind=nag_wp), intent(in) :: t
  integer , intent(in) :: n
  integer, intent(inout) :: ires
  integer :: i,j,k
  real(kind=nag_wp), intent(in) :: y(n)
  real(kind=nag_wp), intent(out) :: fd(n)
  real(kind=dp) :: invgx,of,f1,invgm,jxi,jxf,jmi,jmf,jx,jm, &
& f10,f11,f20,Gv2,Gv3,Gv4,Jaltxs,Jaltms,Jaltxb,Jaltmb,dL0,dLx,dLm
  real(kind=dp) :: fm,fm1,fp,fp1,Um,Um1,Up,Up1,Pxm,Qex,Px,Qex1,Px1,& 
& Pmm,Qem,Pm,Qem1,fs,fs1,Us,Us1,jxs,jms,Ax1m,Ax1x,Dx1pe,Dx1xe,fy1m,fy1m1, &
& fy1p,fy1p1,Qy1pi,Py1m,Py1x,f1y1x,fy1x1,Qy1pe,Qy1pg,Qy1xe,Qy1xg,Qy1xi,&
& rjx11y,rjx1y,rjy11x,fy1x,rjy1x,Ux1x,Ux1x1,Ux1m,Ux1m1
  real(kind=nag_wp) :: edge, rdislc,rdisle,recombine,P0,ss3He,sgbHe,AAx,QQxi,PPx


  integer :: iUpper,ileft,iRight,jbound!, m, k2
  !real(kind=nag_wp) :: javx,javm,jx1,jm1,jx2,jm2,linf,v1,v2,jx11,jx12,sigmasq
  interface
    function gbss( sss,rgr )
      use kindmod
      implicit none
      real(kind=dp) :: sss,rgr
      real(kind=dp) :: gbss
    end function


    function Jaltm(i,j,x,m,ax,am)
      use nuclearvars
      use twodimensionalvars
      implicit none
      integer :: i,j
      real(kind=dp), dimension (0:xnumg+1) :: x
      real(kind=dp), dimension (-1:mnumg+1) :: m
      real(kind=dp) :: Jaltm,ax,am
    end function

    function Jaltx(i,j,x,m,ax,am)
      use nuclearvars
      use twodimensionalvars
      implicit none
      integer :: i,j
      real(kind=dp), dimension (0:xnumg+1) :: x
      real(kind=dp), dimension (-1:mnumg+1) :: m
      real(kind=dp) :: Jaltx,ax,am
    end function

  end interface

  siaabsorption=0.0_dp
  intabsorption=0.0_dp
  vacabsorption=0.0_dp
  Heemission=0.0_dp
  Heabsorption=0.0_dp
  vacemission=0.0_dp
  Hemonomeremission=0.0_dp
  Hevacabsorption=0.0_dp
  sigmavNv=0.0_dp
  f=0.0_dp
  L0=0.0
  L1x=0.0
  L1m=0.0

  do k=5+1,neq1
    i=idarrayi(k)
    j=idarrayj(k)
    f(i,j)=abs(y(k))
   ! if (i .eq. 1 .and. j .eq. 1) print*, k, 'k'
 ! print*, y(k),'y',k,i,j
  end do 

  do k=neq1+1,neq,3
    i=idarrayi(k)
    j=idarrayj(k)
    L0(i,j)=abs(y(k))
    L1x(i,j)=y(k+1)
    L1m(i,j)=y(k+2)
  end do 

  Cv=y(1) !Vacancy concentration
  Ci=y(2) !SIA concentration
  CHe=y(3) !He concentration
  CGb=y(4) !He GB concentration
  Cg=y(5) !SIA glissile concentration
  f(1,0)=Cv
  do j=0,srtgrpm
    f(srtgrpx+1,j)=L0(srtgrpx+1,j)+L1x(srtgrpx+1,j)*(rg1xi(srtgrpx)+1._dp-meanxarr(srtgrpx+1))
    L0(srtgrpx,j)=f(srtgrpx,j)-L1x(srtgrpx,j)*(rg1xi(srtgrpx)-meanxarr(srtgrpx))
    L1x(srtgrpx,j)=0.0_dp
    L1m(srtgrpx,j)=0.0_dp
  end do 

  do i=1,srtgrpx
    f(i,srtgrpm+1)=L0(i,srtgrpm+1)+L1m(i,srtgrpm+1)*(rg1mi(srtgrpm)+1._dp-meanmarr(srtgrpm+1))
    L0(i,srtgrpm)=f(i,srtgrpm)
    L1x(i,srtgrpm)=0.0_dp
    L1m(i,srtgrpm)=0.0_dp
  end do 

  kgCg=2._dp*Cg*kgsq
  !k=6
  
!first half--------------------------------------
  do k=5+1,neq1

    i=idarrayi(k)
    j=idarrayj(k)

    if (rg1mi(j).gt.g00+g11*(rg1xi(i)-1.0_dp)) then
      jbound = j-1                   ! define boundary layer
     cycle
    else
      if (j.eq.mnumg) jbound = mnumg     ! define boundary layer
    endif
    if (j.eq.mnumg) jbound = mnumg 
                 iRight = 1
    if (i.eq.xnumg) iRight = 0             ! He + vac
                                                 ! He + i
                                                     iLeft = 1
    if (i .gt. 1) then
    if (rg1mi(j).gt.g00+g11*(rg1xi(i)-1.0_dp)) iLeft = 0
    end if

                                          iUpper = 1
    if (j.eq.mnumg                         ) iUpper = 0 ! B + He at
    if (rg1mi(j).gt.g00+g11*(rg1xi(i)-1.0_dp))  iUpper = 0 ! B + He at   


!print*, k,i,j
    fm=f(i-1,j)!+L1x(i-1,j)*(rg1xi(i-1)-meanxarr(i-1))
    fm1=f(i,j)!+L1x(i,j)*(rg1xi(i-1)+1._dp-meanxarr(i))
    fp=f(i,j)!+L1x(i,j)*(rg1xi(i)-meanxarr(i))
    fp1=f(i+1,j)!+L1x(i+1,j)*(rg1xi(i)+1._dp-meanxarr(i+1))

    Um=f(i,j-1)!+L1m(i,j-1)*(rg1mi(j-1)-meanmarr(j-1))
    Um1=f(i,j)!+L1m(i,j)*(rg1mi(j-1)+1._dp-meanmarr(j))
    !if (j .eq. 0) Um1=f(i,j)
    Up=f(i,j)!+L1m(i,j)*(rg1mi(j)-meanmarr(j))
    Up1=f(i,j+1)!+L1m(i,j+1)*(rg1mi(j)+1._dp-meanmarr(j+1))

    jxi=Pvarr(i-1,j)*Cv*fm-(Qiarr(i,j)*Ci+Qvarr(i,j)+Qg(i,j)*kgcg)*fm1!*ileft
    !jxf=Pvarr(i,j)*Cv*fp-(Qiarr(i+1,j)*Ci+Qvarr(i+1,j)+Qg(i+1,j)*kgcg)*fp1
    jxf=Pvarr(i,j)*Cv*fp-(Qiarr(i+1,j)*Ci+Qvarr(i+1,j)+Qg(i+1,j)*kgcg)*fp1
    jmi=PHearr(i,j-1)*CHe*Um-Qm(i,j)*Um1
    jmf=PHearr(i,j)*CHe*Up-Qm(i,j+1)*Up1
    !j1(i,j)=jxi!f(i-1,j)
    !j2(i,j)=jxf!f(i,j)
    !j3(i,j)=jmi
    !j4(i,j)=jmf
    dL0=jxi-jxf+jmi-jmf
!if (i .eq. 2 .and. j .eq. 0) print*, dL0,'l0', &
!& PHearr(i,j-1)*CHe*f(i,j-1),Pvarr(i-1,j)*Cv*f(i-1,j),(Qiarr(i,j)*Ci+Qvarr(i,j)+Qg(i))*f(i,j), &
!& Pvarr(i,j)*Cv*f(i,j),(Qiarr(i+1,j)*Ci+Qvarr(i+1,j)+Qg(i+1))*f(i+1,j),'\n',k

    siaabsorption=siaabsorption+Qiarr(i,j)*Ci*f(i,j)!*(ileft)
    vacabsorption=vacabsorption+Pvarr(i,j)*Cv*f(i,j)!*(1-iright)
    Heemission=Heemission+Qm(i,j)*f(i,j)
    Heabsorption=Heabsorption+PHearr(i,j)*CHe*f(i,j)!*(1-iUpper)
    vacemission=vacemission+Qvarr(i,j)*f(i,j)!*(ileft)
		sigmavNv= sigmavNv + Qg(i,j)*f(i,j)!*(ileft)

    fd(k)=dL0
  end do
  
!second half--------------------------------------

  do k=neq1+1,neq,3
    i=idarrayi(k)
    j=idarrayj(k)
!print*, k,i,j

    if (rg1mi(j).gt.g00+g11*(rg1xi(i)-1.0_dp)) then
      jbound = j-1                   ! define boundary layer
     cycle
    else
      if (j.eq.mnumg) jbound = mnumg     ! define boundary layer
    endif
    if (j.eq.mnumg) jbound = mnumg 
                 iRight = 1
    if (i.eq.xnumg) iRight = 0             ! He + vac
                                                 ! He + i
                                                     iLeft = 1
    if (i .gt. 1) then
    if (rg1mi(j).gt.g00+g11*(rg1xi(i)-1.0_dp)) iLeft = 0
    end if

                                          iUpper = 1
    if (j.eq.mnumg                         ) iUpper = 0 
    if (rg1mi(j).gt.g00+g11*(rg1xi(i)-1.0_dp))  iUpper = 0 


    f1=L0(1,j)*rgm(j)
    invgm=(1._dp/rgm(j))
    invgx=(1._dp/rgx(i))
    fm=L0(i-1,j)+L1x(i-1,j)*(rg1xi(i-1)-meanxarr(i-1))
    fm1=L0(i,j)+L1x(i,j)*(rg1xi(i-1)+1._dp-meanxarr(i))
    fp=L0(i,j)+L1x(i,j)*(rg1xi(i)-meanxarr(i))
    fp1=L0(i+1,j)+L1x(i+1,j)*(rg1xi(i)+1._dp-meanxarr(i+1))

    Pxm=Pvarr(i-1,j)*Cv*fm
    Qex=(Qiarr(i,j)*Ci+Qvarr(i,j)+Qg(i,j)*kgcg)*fm1
    Px=Pvarr(i,j)*Cv*fp
    Qex1=(Qiarr(i+1,j)*Ci+Qvarr(i+1,j)+Qg(i+1,j)*kgcg)*fp1
    
    Um=L0(i,j-1)+L1m(i,j-1)*(rg1mi(j-1)-meanmarr(j-1))
    Um1=L0(i,j)+L1m(i,j)*(rg1mi(j-1)+1._dp-meanmarr(j))
    if (j .eq. 0) Um1=L0(i,j)
    Up=L0(i,j)+L1m(i,j)*(rg1mi(j-1)-meanmarr(j))
    Up1=L0(i,j+1)+L1m(i,j+1)*(rg1mi(j)+1._dp-meanmarr(j+1))

    Pmm=PHearr(i,j-1)*CHe*Um
    Qem=Qm(i,j)*Um1
    Pm=PHearr(i,j)*CHe*Up
    Qem1=Qm(i,j+1)*Up1

    fs  = L1x(i-1,j-1)*(-0.5_dp) + L0(i-1,j-1)
    fs1 = L1x(i-1,j-1)*( 0.5_dp) + L0(i-1,j-1)
    
    Us  = L1m(i-1,j-1)*(-0.5_dp) + L0(i,j-1) 
    Us1 = L1m(i-1,j-1)*(0.5_dp) + L0(i,j-1)  

    jxs=Pvarr(i-1,j-1)*Cv*fs-(Qiarr(i-1,j-1)*Ci+Qvarr(i-1,j-1)+Qg(i-1,j-1)*kgcg)*fs1 !Jaltx(i-1,j-1,meanxarr,meanmarr,-0.5_dp,0._dp)!
    jms=fPHearr(i-1,j-1)*CHe*Us-Qm(i-1,j-1)*Us1!Jaltm(i-1,j-1,meanxarr,meanmarr,0._dp,-0.5_dp)!

    
    fy1m  = fm  + L1m(i-1,j)
    fy1x  = fp + L1m(i,j)
    fy1x1 = fp1 + L1m(i+1,j)
    fy1m1 = fm1 + L1m(i,j)

    
    Ux1m = Um + L1x(i,j-1)
    Ux1x  = Up  + L1x(i,j)
    Ux1m1 = Um1 + L1x(i,j)
    Ux1x1 = Up1 + L1x(i,j+1)


    
    Ax1m = PHearr(i,j-1)*Ux1m*CHe          
    Ax1x  = PHearr(i,j)*Ux1x*CHe                      

    Dx1xe = Qm(i,j)*Ux1m1                      
    Dx1pe = Qm(i,j+1)*Ux1x1                       

    Ax1x = Ax1x*iUpper                             

    rJx11y = ( Ax1m - Dx1xe )
    rJx1y  = ( Ax1x - Dx1pe )



    Py1m = Pvarr(i-1,j)*fy1m*Cv         
    Py1x  = Pvarr(i,j)*fy1x *Cv                 
    Qy1xi = Qiarr(i,j)*fy1m1*Ci                    
    Qy1pi = Qiarr(i+1,j)*fy1x1*Ci                    


    Qy1xg = Qg(i,j)*fy1m1*kgcg                
    Qy1pg = Qg(i+1,j)*fy1x1*kgcg               


    Qy1xe = Qvarr(i,j)*fy1m1                      
    Qy1pe = Qvarr(i+1,j)*fy1x1                       

    AAx = Pm*(1-iUpper)                               ! Boundaries
    Pm  = Pm*iUpper
    QQxi = (Qex)*(1-iLeft )                 ! Boundary
    PPx  = Px*(1-iRight)                       ! Boundary
    !QQxg =  Qxg     *(1-iLeft )


    jxi=(Pxm-Qex)!*ileft)!Jaltx(i-1,j,rg1xi,meanmarr,0._dp,0._dp)!
    jxf=(Px-Qex1)!Jaltx(i,j,rg1xi,meanmarr,0._dp,0._dp)!
    jmi=(Pmm-Qem)!Jaltm(i,j-1,meanxarr,rg1mi,0._dp,0._dp)!Pmm-Qem
    jmf=(Pm-Qem1)

    Qy1xi = Qy1xi*iLeft                              
    Qy1xe = Qy1xe*iLeft                             
    Qy1xg = Qy1xg*iLeft                             
    Py1x  = Py1x *iRight                                        

    rJy11x = ( Py1m - (Qy1xi + Qy1xe + Qy1xg) )
    rJy1x  = ( Py1x - (Qy1pi + Qy1pe + Qy1pg) )

    of= L0(i,j)*rgx(i)*rgm(j)






    !(jm1-jm2)/(2.*rgx(i)*sigmasq(rgx(i)))
    !(jx11-jx12)/(2.*rgm(j)*sigmasq(rgm(j)))
    dL0 = invgx*(Jxi-Jxf) + invgm*(Jmi-Jmf)
    
    dLx = -itermx(i)*(Jxi+Jxf- 2._dp*jxs) + &
& invgm*(rJx11y-Jmi)*zerox(i) - invgm*(rJx1y-jmf)*zerox(i) 

!    dLx = -itermx(i)*(Jxi+Jxf- 2._dp*jxs) + &
!& invgm*(L1x(i,j-1)*(PHearr(i,j-1)*CHe-Qm(i,j-1)))*zerox(i) -  &
!& invgm*(L1x(i,j)*(PHearr(i,j)*CHe-Qm(i,j)))*zerox(i) 

!    dLx = -itermx(i)*(Jxi+Jxf- 2._dp*jxs) + &
!& invgm*(L1x(i,j-1)*PHearr(i,j-1)*CHe-Qm(i,j)*L1x(i,j))*zerox(i) -  &
!& invgm*(L1x(i,j)*PHearr(i,j)*CHe-Qm(i,j+1)*L1x(i,j+1))*zerox(i) 
    !dLx = -itermx(i)*(Jxi+Jxf- 2._dp*Jaltx(i,j,meanxarr,meanmarr,-0.5_dp,0._dp)) + &
!& invgm*( Jaltm(i,j-1,meanxarr,rg1mi,1._dp,0._dp)-Jaltm(i,j-1,meanxarr,rg1mi,0._dp,0._dp) ) - &
!& invgm*(Jaltm(i,j,meanxarr,rg1mi,1._dp,0._dp)-jmf  )

    dLm =-itermm(j)*(Jmi+Jmf- 2._dp*jms) + &
& invgx*(rJy11x-Jxi)*zerom(j) - invgx*(rJy1x-jxf)*zerom(j)

!    dLm =-itermm(j)*(Jmi+Jmf- 2._dp*jms) + &
!& invgx*(L1m(i-1,j)*Pvarr(i-1,j)*Cv-L1m(i,j)*(Qiarr(i,j)*Ci+Qvarr(i,j)+Qg(i,j)*kgcg))*zerom(j) &
!&  - invgx*(L1m(i,j)*Pvarr(i,j)*Cv-L1m(i+1,j)*(Qiarr(i+1,j)*Ci+Qvarr(i+1,j)+Qg(i+1,j)*kgcg))*zerom(j)

!    dLm =-itermm(j)*(Jmi+Jmf- 2._dp*jms) + &
!& invgx*(L1m(i-1,j)*(Pvarr(i-1,j)*Cv-(Qiarr(i-1,j)*Ci+Qvarr(i-1,j)+Qg(i-1,j)*kgcg)))*zerom(j) &
!&  - invgx*(L1m(i,j)*(Pvarr(i,j)*Cv-(Qiarr(i,j)*Ci+Qvarr(i,j)+Qg(i,j)*kgcg)))*zerom(j)
    

    siaabsorption=siaabsorption+Qiarr(i,j)*Ci*of!*(ileft)
    vacabsorption=vacabsorption+fPvarr(i,j)*Cv*of!*(1-iright)
    Heemission=Heemission+Qm(i,j)*of
    Heabsorption=Heabsorption+PHearr(i,j)*CHe*of!*(1-iUpper)
    vacemission=vacemission+Qvarr(i,j)*of!*(ileft)
		sigmavNv= sigmavNv + Qg(i,j)*of!*(ileft)
      
    fd(k)=dL0
    fd(k+1)=dLx
    fd(k+2)=dLm
!if (dlx .gt. 1.E-30) print*, i,j,dL0,dLx,dLm
  end do    
  sigmavNv=sigmavNV/Dg   

  call speciesgenerationHe(GHe,t)
  f10=f(1,0)
  f11=f(1,1)
  f20=f(2,0)

  recombine=nui*Di*Cv*Ci 
  P0=2._dp*Pvarr(1,0)*f10*Cv
  vacemission=vacemission+Qm(1,1)*f11+ (Qiarr(2,0)*Ci+Qvarr(2,0)+Qg(2,0)*kgCg)*f20 !&
!& + Qiarr(1,1)*f11*Ci
  Heabsorption=Heabsorption+PHearr(1,0)*f10*CHe
  Heemission=Heemission+ Qg(1,1)*kgCg*f11+Qiarr(1,1)*f11*Ci


  ss3v =Zv*rhod +(vacabsorption+P0)/(Dv*Cv)+PHearr(1,0)*f10/DHe!sv_dis +(vacabsorption+P0)/(Dv*Cv)+fPHearr(1,0)*f10/DHe!*CHe 
  P0=P0+PHearr(1,0)*CHe*f10
  ss3i=Zi*rhod + siaabsorption/(Di*Ci)!si_dis + siaabsorption/(Di*Ci)
  !ss3v =ss3v+(vacabsorption+PHearr(1,0)*f10*CHe+Pvarr(1,0)*f10*Cv)/(Dv*Cv)
  !ss3i=ss3i + siaabsorption/(Di*Ci)
  tgdid=tgdi!/(rhod*vol**(2._dp/3._dp))
!print*, ss3v,(vacabsorption+P0)/(Dv*Cv)+fPHearr(1,0)*f10/DHe
  !rdisle = DHe*CHe*rhod*tgdid*dexp(-edge/(kB*temp)) 
  disl=ZHe*DHe*CHe*rhod*(disl0/(1.E+04*tgdid+disl0))!(1._dp-tgdid)

  !if (disl .le. 0.0) disl=0.0
  ss3He=(disl+Heabsorption)/(DHe*CHe)


  sgbv = gbss(ss3v,l)
  sv_gb=sgbv


  sgbi = gbss(ss3i,l)
  si_gb= sgbi

  sgbHe=gbss(ss3He,l)
  sHe_gb=sgbHe





  dCv = Gv  -recombine - Dv*(Cv-Cv0)*sv_dis  -Dv*(Cv-Cv0)*sv_gb & 
& - P0 - vacabsorption + vacemission -QQxi*rgm(j)

  dCi = Geni - recombine -Di*Ci*si_dis -siaabsorption  - Di*Ci*si_gb -PPx*rgm(j)

  dCHe = GHe  + Heemission  -disl &! disl 
& -Heabsorption  -DHe*CHe*sHe_gb -AAx*rgx(i)


  dCGb= DHe*CHe*sHe_gb
  kgsq=(fstterm+secndterm+sigmavNv)
  dCg= Gg-Dg*Cg*2_dp*kgsq*kgsq-Qg(1,1)*kgcg*f11



  do i=2, maxgv
    fd(5+i)=fd(5+i)+genv(i)
  end do 

  fd(1)=dCv
  fd(2)=dCi
  fd(3)=dCHe
  fd(4)=dCGb
  fd(5)=dCg
  fd(6)=dCv
  return
end subroutine
