function diff(temp,Em,D0)
!------------------------------------------------------------------
! Diffusion coefficient function
!------------------------------------------------------------------  
  use kindmod
  implicit none
  real(kind=dp),parameter :: kb =8.6173324*10.0**(-5.0)!1.3806488_dp*10.0**(-23.0_dp)
  real(kind=dp),intent(in):: temp,Em,D0
  real(kind=dp):: diff
  if (temp==0.0_dp) then
    diff =D0
  else
    diff= D0*exp(-Em/(kb*temp))
  end if
!print*, D0,Em,kB*temp
  return
end function diff

function Qvf5(xr,mr,i,j)
  use nuclearvars
  use twodimensionalvars
  implicit none
  !real(kind=dp),parameter :: Pi = 3.1415927_dp, kb =8.6173324*10.0**(-5.0)
  real(kind=dp):: Qvf5,w,sengc,e2v,gp,wolf,wolfa,wolfn
  real(kind=dp):: x,xr, mr, Z,Eb,r1,mjoulestoev
  integer :: i,j
  interface
      subroutine compressibility (x,Z) 
      use nuclearvars
      implicit real*8 (a-h,o-z)
      end subroutine
  end interface

  E2v=0.80
  wolfn=1.
  mjoulestoev=6.242E18*1.E-03
  sengc=seng*mjoulestoev
  x  = ( mr)/(xr)
  wolfa=0.0
  r1=(3.*xr*vol/(4.*Pi))**(1._dp/3._dp)
  call compressibility (x,Z)              ! EOS for Gas atoms

  gp  = x*Z*(kB*temp)/vol                   ! gas pressure 

                                          ! Wolfer corr
  if (i .eq. 2 .and. j.eq. 0) then
    wolfa = 4.0d0*(1.0d0-(Ev-e2v)*r1/(2.0d0*sengc*vol)) 
  end if
!                                                   ! v-bubble binding E
  if (i .gt. 10000 .or. i.eq. 1) then
     bindengV(i,j) = Ev
  else
    wolf = 1.0d0-wolfa/(4.0d0+(xr-2.0d0)/(2.0d0+wolfn))
    bindengV(i,j) = Ev - (2.0d0*sengc*wolf/r1-gp)*vol
  endif
  w= ((48.0_dp*Pi*Pi)/(vol*vol))**(1.0_dp/3.0_dp)
  Qvf5 = w*xr**(1._dp/3._dp)*Dv*Cv0*dexp(-bindengV(i,j)/(kB*temp))
  return
end function

function Qvf4(x,m,i,j)
!------------------------------------------------------------------
! Vacancy emission rate function 
!------------------------------------------------------------------ 
  use nuclearvars
  use twodimensionalvars
  implicit none
  !real(kind=dp),parameter :: Pi = 3.1415927_dp, kb =8.6173324*10.0**(-5.0)
  real(kind=dp):: Qvf4,w,sengc,e2v,gp,wolf,wolfa,wolfn
  real(kind=dp):: x, m, Z,Eb,r1,mjoulestoev
  integer :: i,j
  interface
      subroutine compressibility (x,Z) 
      use nuclearvars
      implicit real*8 (a-h,o-z)
      end subroutine
  end interface

  wolfn=1.
  E2v=0.20
  
  mjoulestoev=6.242E18*1.E-03
  sengc=seng*mjoulestoev
  !m=m+1.
  call compressibility(m/x,Z)              ! EOS for Gas atoms
  !print*, Z, m/x, m,x


  gp  = m*Z*kB*temp/vol/x                  ! gas pressure 
  wolfa=0.0
  wolf=0.0

  !wolfa = 4._dp*(1._dp-(Ev-E2v)*r1/(2._dp*sengc*vol))
  if (i .eq. 1 .and. j .ge. 0)  then 
    wolfa=0.0
  else 
    r1=(3._dp*2._dp*vol/(4._dp*Pi))**(1._dp/3._dp)
    wolfa = 4._dp*(1._dp-(Ev-E2v)*r1/(2._dp*sengc*vol))
  end if
  r1=(3._dp*x*vol/(4._dp*Pi))**(1._dp/3._dp)
  
                                          ! Wolfer corr
  !if (x .eq. 2 .and. m .eq. 0) then
  !   wolfa = 4._dp*(1._dp-(Ev-E2v)*r1/(2._dp*sengc*vol)) 
  !endif
!                                                   ! v-bubble binding E
  if (x .gt. 1000 .or. x .eq. 1) then
     Eb = Ev
  else
     wolf = 1._dp-wolfa/(4._dp+(x-2._dp)/(2._dp+wolfn))
     Eb = Ev - (2._dp*sengc*wolf/r1-gp)*vol
  endif
  !write(32,45) i,j,Eb,gp,wolfa,Z,wolf,r1
  !45 format(2(I3,X),6(1pd20.8,X))
  !write(32,*) i,j,Eb,gp,kB*temp,m/x,sengc,wolfa,Z,wolf
  bindengV(i,j)=Eb
  w= ((48.0_dp*Pi*Pi)/(vol*vol))**(1.0_dp/3.0_dp)
  Qvf4 = w*x**(1._dp/3._dp)*Dv*Cv0*dexp(-Eb/(kB*temp)) ! RC v-emiss'
  !write(223,*) Eb,x,m,'t1'
  return
end function Qvf4


function Qvf3(x,m,Dv,vol,seng,Ev,temp,i,j)
!------------------------------------------------------------------
! Vacancy emission rate function 
!------------------------------------------------------------------ 
  use twodimensionalvars
  implicit none
  real(kind=dp),parameter :: Pi = 3.1415927_dp, kb =8.6173324*10.0**(-5.0)
  real(kind=dp):: D, Qvf3, w, Dv, vol,  &
& seng, sengc, temp, Ev,e2v
  real(kind=dp):: x, m, &
& ratio,bindeng
  integer :: i,j
  E2v=0.2
  ratio= (m/x)
  if (m .eq. 0.0) then 
    bindeng = Ev
  else
!print*, ratio, m,x,log10(ratio)
    bindeng=1.99+3.08*log10(ratio)+2.7*log10(ratio)*log10(ratio)
  end if 
  
  if (bindeng .lt. 0.5) bindeng =0.5
  if (bindeng .gt. 5.5) bindeng =5.5
  if (m .eq. 0) bindeng=Ev+(E2v-Ev)*&
& ((x**(2._dp/3._dp)-(x-1._dp)**(2._dp/3._dp))/(2._dp**(2._dp/3._dp)-1._dp))
  if (m .eq. 0 .and. x .eq. 2) bindeng=0.80
  if (m .eq. 0 .and. x .eq. 3) bindeng=0.92
  if (m .eq. 0 .and. x .eq. 4) bindeng=1.64
  bindengV(i,j)=bindeng
  !write(88,*) ratio,'ratio',bindeng,'bind',x,m
  !write(88,'(EN20.10,X,EN20.10)') ratio,bindeng
  w= ((48.0_dp*Pi*Pi)/(vol*vol))**(1.0_dp/3.0_dp)
  Qvf3=w*x**(1._dp/3._dp)*Dv*dexp(-bindeng/(kB*temp))
  !write(555,*) x,m, Qvf3
  return
end function Qvf3


function Qvf2(x,m,Dv,vol,seng,Ev,temp,i,j)
  use twodimensionalvars
  implicit none

  real(kind=dp),parameter :: Pi = 3.1415927_dp, kb =8.6173324*10.0**(-5.0)
  real(kind=dp):: D, Qvf2, w, Dv, vol, E,y,alpha,Ev,seng,sengc,temp,Z
  real(kind=dp):: x,m,mjoulestoev
  integer :: i,j 
  mjoulestoev=6.242E18*1.E-03
  !xr=real(x)
  !mr=real(m)
  sengc=seng*mjoulestoev
  w= ((48.0_dp*Pi*Pi)/(vol*vol))**(1.0_dp/3.0_dp)
  alpha = 2.*sengc*(4.*Pi*vol*vol/3.0_dp)**(1.0_dp/3.0_dp)
  d= 0.3135*(0.8542-0.03996*log(temp/9.16))
  y= ((Pi*d**(3._dp))/(6._dp*vol))*(m/x)
  Z=(1.+y+y*y-y**(3._dp))/((1.-y)**(3._dp))
  E= Ev- alpha*x**(-1.0_dp/3.0_dp)+(m/x)*Z*kb*temp
  bindengV(i,j)=E
  Qvf2 = w*Dv*x**(1.0_dp/3.0_dp)*exp(-E/(kb*temp))
  !print*,alpha,'alpha',d,'d',w,'w',y,'y',z,'z',e,'E',qvf,'qvf'
  !print*,xr,'x',mr,'m',vol,'vol',sengc,'seng',Ev,'Ev',temp,'temp',dv,'dv'
  return
end function Qvf2

function Qvf(x,m,Dv,vol,seng,Ev,temp,i,j)
  use twodimensionalvars
  implicit none

  real(kind=dp),parameter :: Pi = 3.1415927_dp, kb =8.6173324*10.0**(-5.0)
  real(kind=dp):: D, Qvf, w, Dv, vol, E,y,alpha,Ev,seng,sengc,temp,Z
  real(kind=dp):: x,m,mjoulestoev,rho,Vm,Zm,b,T
  integer :: i,j
  mjoulestoev=6.242E18*1.E-03
  !xr=real(x)
  !mr=real(m)
  sengc=seng*mjoulestoev
  T=kB*temp
  w= ((48.0_dp*Pi*Pi)/(vol*vol))**(1.0_dp/3.0_dp)
  alpha = 2._dp*sengc*(4._dp*Pi*vol*vol/3.0_dp)**(1.0_dp/3.0_dp)
  Vm=56._dp*temp**(-0.25_dp)*exp(-0.145_dp*temp**(0.25_dp))
  Zm=0.1225*Vm*temp**(0.555_dp)
  b=170.*temp**(-1._dp/3._dp)-1750.*temp**(-1._dp)
  rho=(Vm/vol)*(m/x)*1.E-30
  Z=(1._dp-rho)*(1._dp+rho-52._dp*rho*rho)+(b/Vm)*rho*(1._dp-rho)*(1._dp-rho) &
& +Zm*rho*rho*(3._dp-2._dp*rho)
  E= Ev- alpha*x**(-1.0_dp/3.0_dp)+(m/x)*Z*T
  if (x .eq. 1.) E=Ev
  bindengV(i,j)=E
  Qvf = w*Dv*x**(1.0_dp/3.0_dp)*exp(-E/(T))
  !print*,alpha,'alpha',d,'d',w,'w',y,'y',z,'z',e,'E',qvf,'qvf'
  !print*,xr,'x',mr,'m',vol,'vol',sengc,'seng',Ev,'Ev',temp,'temp',dv,'dv'
   !write(555,*) x,m, Qvf,E
  return
end function Qvf

  
function Qif(x,D,vol,rcap)
!------------------------------------------------------------------
! Vacancy emission rate function 
!------------------------------------------------------------------ 
  implicit none
  integer, parameter :: dp=selected_real_kind(15,300)
  real(kind=dp),parameter :: Pi = 3.1415927_dp
  real(kind=dp) :: Qif, w, D, vol,rcap
  real(kind=dp) :: x, r, Zi_void
  r=(3.*x*vol/(4.*Pi))**(1._dp/3._dp)
  Zi_void= 1.+rcap/r
  Qif=4.*Pi*r*Zi_void*D/vol
  !print*, D,'d',w,'w'
  return
end function Qif
  
function Pf(x,D,vol)
  use kindmod
  implicit none
  real(kind=dp),parameter :: Pi = 3.1415927_dp
  real(kind=dp) Pf, w, D, vol,xr
  real(kind=dp) :: x
  
  w= ((48.0_dp*Pi*Pi)/(vol*vol))**(1.0_dp/3.0_dp)
  Pf=w*D*x**(1.0_dp/3.0_dp)
  !print*,D,'D',w*x**(1.0_dp/3.0_dp),'w',pf
  return
end function Pf


function Jaltx(i,j,x,m,ax,am)
  use nuclearvars
  use twodimensionalvars
  implicit none  
  integer :: i,j
  real(kind=dp), dimension(0:xnumg+1) :: x
  real(kind=dp), dimension(-1:mnumg+1) :: m
  real(kind=dp) :: Jaltx,ax,am


  interface
    function linf(L0,Lx,Lm,x,meanx,m,meanm)
      use kindmod
      implicit none 
      real(kind=dp) :: linf,meanx,meanm,x,m,L0,Lx,Lm
    end function
  end interface

  Jaltx=Pvarr(i,j)*Cv*linf( L0(i,j),L1x(i,j),L1m(i,j),x(i)+ax,meanxarr(i),m(j)+am,meanmarr(j) ) - &
& (Qiarr(i+1,j)*Ci+Qvarr(i+1,j) +Qg(i+1,j)*kgCg)*linf( L0(i+1,j),L1x(i+1,j),L1m(i+1,j),&
& x(i)+1._dp+ax,meanxarr(i+1),m(j)+am,meanmarr(j) )  
!print*, linf( L0(i,j),L1x(i,j),L1m(i,j),x(i)+ax,meanxarr(i),m(j)+am,meanmarr(j) ),'f'
!write(8,*) L0(i,j),L1x(i,j),L1m(i,j),'Ls',i,j
!write(8,*) Pvarr(int(x(i)))*Cv,Qiarr(int(x(i)))*Ci+Qvarr(int(x(i)),int(m(j))), 'P & Q',i,j
!print*, linf( L0(i+1,j),L1x(i+1,j),L1m(i+1,j),&
!& x(i+1)+ax,meanxarr(i+1),m(j)+am,meanmarr(j) ) ,L0(i+1,j),L1x(i+1,j),L1m(i+1,j)
!print*, x(i+1), int(x(i+1)),i,j

  return
end function  



function Jaltxs(i,j)
  use nuclearvars
  use twodimensionalvars
  implicit none  
  integer :: i,j
  real(kind=dp) :: Jaltxs


  interface
    function linf(L0,Lx,Lm,x,meanx,m,meanm)
      use kindmod
      implicit none 
      real(kind=dp) :: linf,meanx,meanm,x,m,L0,Lx,Lm
    end function
  end interface

  Jaltxs=Pvarr(i,j)*Cv*linf( L0(i,j),L1x(i,j),L1m(i,j),meanxarr(i)-0.5_dp,meanxarr(i),meanmarr(j),meanmarr(j) ) - &
& (Qiarr(i,j)*Ci+Qvarr(i,j) +Qg(i,j))*linf( L0(i,j),L1x(i,j),L1m(i,j),&
& meanxarr(i)+0.5_dp,meanxarr(i),meanmarr(j),meanmarr(j) )
  return
end function 


function Jaltm(i,j,x,m,ax,am)
  use nuclearvars
  use twodimensionalvars
  implicit none
  integer :: i,j
  real(kind=dp), dimension (0:xnumg+1) :: x
  real(kind=dp), dimension (-1:mnumg+1) :: m
  real(kind=dp) :: Jaltm,ax,am

  interface
    function linf(L0,Lx,Lm,x,meanx,m,meanm)
      use kindmod
      implicit none 
      real(kind=dp) :: linf,meanx,meanm,x,m,L0,Lx,Lm
    end function
  end interface

  Jaltm=PHearr(i,j)*CHe*linf( L0(i,j),L1x(i,j),L1m(i,j),x(i)+ax,meanxarr(i),m(j)+am,meanmarr(j) ) - &
& Qm(i,j+1)*linf( L0(i,j+1),L1x(i,j+1),L1m(i,j+1),x(i)+ax,meanxarr(i),m(j)+1._dp+am,meanmarr(j+1) )  
!write(8,*) PHearr(int(x(i)))*CHe,Qm(int(m(j))),'P & Q',i,j
!      if (j.eq.mnumg) then
!        write(13,*) PHearrg(i)*CHe, Qmg(j+1),linf( L0(i,j),L1x(i,j),L1m(i,j),x(i)+ax,meanxarr(i),m(j)+am,meanmarr(j) )&
!&,1./gm(j), Jaltm/gm(j),j,i,Che,'Jx,Jm'
!      end if 
  return
end function  

function Jaltms(i,j)
  use nuclearvars
  use twodimensionalvars
  implicit none
  integer :: i,j
  real(kind=dp) :: Jaltms

  interface
    function linf(L0,Lx,Lm,x,meanx,m,meanm)
      use kindmod
      implicit none 
      real(kind=dp) :: linf,meanx,meanm,x,m,L0,Lx,Lm
    end function
  end interface

  Jaltms=fPHearr(i,j)*CHe*linf( L0(i,j),L1x(i,j),L1m(i,j),meanxarr(i),meanxarr(i),meanmarr(j)-0.5_dp,meanmarr(j) ) - &
& Qm(i,j)*linf( L0(i,j),L1x(i,j),L1m(i,j),meanxarr(i),meanxarr(i),meanmarr(j)+0.5_dp,meanmarr(j) ) 
  return
end function


function Jaltxb(i,j,x,m,ax,am)
  use nuclearvars
  use twodimensionalvars
  implicit none  
  integer :: i,j
  real(kind=dp), dimension(0:xnumg+1) :: x
  real(kind=dp), dimension(-1:mnumg+1) :: m
  real(kind=dp) :: Jaltxb,ax,am

  interface
    function linf(L0,Lx,Lm,x,meanx,m,meanm)
      use kindmod
      implicit none 
      real(kind=dp) :: linf,meanx,meanm,x,m,L0,Lx,Lm
    end function
  end interface 

  Jaltxb= -(Qiarr(i+1,j)*Ci+Qvarr(i+1,j) +Qg(i+1,j))*linf( L0(i+1,j),L1x(i+1,j),L1m(i+1,j),&
& x(i)+1.+ax,meanxarr(i+1),m(j)+am,meanmarr(j) )  
  return
end function  

function Jaltmb(i,j,x,m,ax,am)
  use nuclearvars
  use twodimensionalvars
  implicit none
  integer :: i,j
  real(kind=dp), dimension (0:xnumg+1) :: x
  real(kind=dp), dimension (-1:mnumg+1) :: m
  real(kind=dp) :: Jaltmb,ax,am

  interface
    function linf(L0,Lx,Lm,x,meanx,m,meanm)
      use kindmod
      implicit none 
      real(kind=dp) :: linf,meanx,meanm,x,m,L0,Lx,Lm
    end function
  end interface

  Jaltmb= - Qm(i,j+1)*linf( L0(i,j+1),L1x(i,j+1),L1m(i,j+1),x(i)+ax,meanxarr(i),m(j)+1.+am,meanmarr(j+1) )   
  return
end function  



function sigmasq(dx)
  use kindmod
  implicit none
  real(kind=dp) :: dx
  real(kind=dp) :: sigmasq
  sigmasq=(dx*dx-1._dp)/12._dp !- xi*xi+xi
  return
end function

function mean(xi,dx)
  use kindmod
  implicit none
  real(kind=dp) :: xi,dx
  real(kind=dp) :: mean
  mean=xi-0.5_dp*(dx-1._dp)
  return
end function

function linf(L0,Lx,Lm,x,meanx,m,meanm)
  use kindmod
  implicit none 
  real(kind=dp) :: linf,meanx,meanm,x,m,L0,Lx,Lm
  linf = L0+Lx*(x-meanx)+Lm*(m-meanm)
!if (x .ne. meanx .and. lx .gt. 0.0) print*, x,meanx,Lx*(x-meanx),lx
  return
end function

function ksqst(l,Rg)
  use kindmod
  implicit none
  real(kind=dp) :: l,Rg,ksqst
  ksqst=dsqrt(2._dp/(l*(2._dp*Rg-l)))
  return
end function ksqst

function ksqft(rhod,disabs)
  use kindmod
  implicit none
  !integer, parameter :: dp=selected_real_kind(15,300)
  real(kind=dp),parameter :: Pi = 3.1415927_dp
  real(kind=dp) :: rhod,disabs,ksqft
  ksqft=Pi*rhod*disabs*0.5_dp
  return
end function ksqft

function Qgcoeff(Dg,vol) 
  use kindmod
  implicit none
  real(kind=dp),parameter :: Pi = 3.1415927_dp  
  real(kind=dp) :: Qgcoeff,Dg,vol,r
  r=(0.75_dp*vol/Pi)**(1._dp/3._dp)
  Qgcoeff= Pi*r*r*Dg/vol
  return
end function Qgcoeff

function gbss( sss,rgr )
  use kindmod
  implicit none
  real(kind=dp) ::rkr,sss,fun,rgr,gbss
  rkr  = rgr*dsqrt(abs(sss))
  fun  = 3._dp*(rkr/dtanh(rkr)-1._dp)
  gbss = rgr**(-2._dp)*rkr**2._dp*fun/(rkr**2._dp-fun)
  return
end

function Qemiss(x,m,D,vol,EbHe,temp,i,j)
  use twodimensionalvars
  implicit none
  real(kind=dp),parameter :: Pi = 3.1415927_dp
  real(kind=dp),parameter :: kb =8.6173324*10.0**(-5.0)
  real(kind=dp) :: x, m,D,vol,EbHe,temp,Qemiss,w,Z
  integer :: i,j
  interface
      subroutine compressibility (x,Z) 
      use nuclearvars
      implicit real*8 (a-h,o-z)
      end subroutine
  end interface

  bindengHe(i,j)=EbHe
  !m=m+1.
  call compressibility(m/x,Z)
  w= ((48.0_dp*Pi*Pi)/(vol*vol))**(1.0_dp/3.0_dp)
  Qemiss=w*x**(1._dp/3._dp)*D*dexp(-EbHe/(kB*temp))*Z*m/x
  !write(32,45) i,j,EbHe,Z,Qemiss
  !45 format(2(I3,X),3(1pd20.8,X))
  return
end function Qemiss

function Qemissalt(x,m,D,vol,temp,i,j)
  use kindmod
  use twodimensionalvars
  implicit none
  real(kind=dp),parameter :: Pi = 3.1415927_dp
  real(kind=dp),parameter :: kb =8.6173324*10.0**(-5.0)
  real(kind=dp) :: x,m, D,vol,temp,Qemissalt,w, &
& ratio,bindeng,atmpfreq,Z
  integer :: i,j
  interface
      subroutine compressibility (x,Z) 
      use nuclearvars
      implicit real*8 (a-h,o-z)
      end subroutine
  end interface


  atmpfreq=1.E-0
  !print*, m,x
  ratio= (m/x)
  if (m.eq. 0.0) then
    bindeng=0.0
  else
    bindeng=2.20_dp-1.55_dp*log10(ratio)-0.53_dp*log10(ratio)*log10(ratio)
  end if
  if (bindeng .lt. 0.5) bindeng =0.5
  if (bindeng .gt. 3.5) bindeng =3.5
  if (m .eq. 1 .and. x .eq. 1) bindeng=2.4
  if (m .eq. 2 .and. x .eq. 1) bindeng=1.17
  if (m .eq. 3 .and. x .eq. 1) bindeng=1.39
  if (m .eq. 4 .and. x .eq. 1) bindeng=1.03
  if (m .eq. 4 .and. x .eq. 2) bindeng=1.04

  bindengHe(i,j)=bindeng
  !m=m+1.
  call compressibility(m/x,Z)
  !write(88,*) ratio,'ratio',bindeng,'bind',x,m
  !write(88,'(EN20.10,X,EN20.10)') ratio,bindeng
  w= ((48.0_dp*Pi*Pi)/(vol*vol))**(1.0_dp/3.0_dp)
  !bindeng=2.3
  Qemissalt=w*x**(1._dp/3._dp)*D*atmpfreq*dexp(-bindeng/(kB*temp))*m/x*Z
  if (m .eq. 0.0) Qemissalt=0.0
  !write(88,'(EN20.10,X,EN20.10)') ratio,bindeng,Qemissalt
  return
end function Qemissalt

function poisson(mean)
  use kindmod
	implicit none
	real(kind=dp) :: t,mean,rand
	integer(kind=ip) :: poisson
	t=0.0
	poisson=0
	do while (t.lt. mean) 
		call random_number(rand)
		t=t -log(1.-rand)
		poisson=poisson+1
		!print*, t, 't',mean
	end do 
	poisson= poisson-1
	return
end function poisson



function binomial(n,p)
  use kindmod
	implicit none
	integer(kind=ip) :: binomial,n,i
	real(kind=dp) :: p,rand
	binomial=0
	do i=1,n 
	  call random_number(rand)
	  if (rand .lt. 1.-p) then 
			binomial=binomial
	  else 
			binomial=binomial+1
	  end if 
	end do
	return
end function binomial

function find_nearest(array, value, arraydim)
  use kindmod
	implicit none
	integer(kind=ip) ::find_nearest,arraydim,i
	integer(kind=ip),dimension(1) :: x
	real(kind=dp) :: value
	real(kind=dp),dimension(2,1:arraydim) :: array
  real(kind=dp), dimension(1:arraydim) :: test
	do i=1,arraydim
!print*, array(2,i),array(1,i),value,abs(array(2,i)-value),i
		test(i)= abs(array(2,i)-value)
	end do   

	x= minloc(test(1:arraydim))

	find_nearest=x(1)

	return 
end function

