

function Qvf1d(x,Dv,vol,E2v,Ev,temp)
  use kindmod
  implicit none
  real(kind=dp),parameter :: Pi = 3.1415927_dp, kb =8.6173324*10.0**(-5.0)
  real(kind=dp):: D, Qvf1d, w, Dv, vol, E,Ev,temp,x,E2v
  w= ((48.0_dp*Pi*Pi)/(vol*vol))**(1.0_dp/3.0_dp)
  E=Ev+(E2v-Ev)*((x**(2._dp/3._dp)-(x-1._dp)**(2._dp/3._dp))/(2._dp**(2._dp/3._dp)-1._dp))
  Qvf1d = w*Dv*x**(1.0_dp/3.0_dp)*exp(-abs(E)/(kb*temp))
  !print*, E,'e',Qvf,w,Dv,x
  return
end function Qvf1d

function Qvf21d(x,m,Dv,vol,seng,Ev,temp)
  use onedimensionalvars
  implicit none
  real(kind=dp),parameter :: Pi = 3.1415927_dp, kb =8.6173324*10.0**(-5.0)
  real(kind=dp):: D, Qvf21d, w, Dv, vol, E,y,alpha,Ev,seng,sengc,temp,Z
  real(kind=dp):: x,m,mjoulestoev
  integer :: i
  mjoulestoev=6.242E18*1.E-03
  !xr=real(x)
  !mr=real(m)
  sengc=seng*mjoulestoev
  w= ((48.0_dp*Pi*Pi)/(vol*vol))**(1.0_dp/3.0_dp)
  alpha = 2.*sengc*(4.*Pi*vol*vol/3.0_dp)**(1.0_dp/3.0_dp)
  d= 0.3135*(0.8542-0.03996*log(temp/9.16))
  y= ((Pi*d**(3._dp))/(6._dp*vol))*m/x
  Z=(1.+y+y*y-y**(3._dp))/((1.-y)**(3._dp))
  E= Ev- alpha*x**(-1.0_dp/3.0_dp)+(m/x)*Z*kb*temp
  i=nint(x)
  bindingV(i)=E
  Qvf21d = w*Dv*x**(1.0_dp/3.0_dp)*exp(-E/(kb*temp))
  !print*,alpha,'alpha',d,'d',w,'w',y,'y',z,'z',e,'E',qvf,'qvf'
  !print*,xr,'x',mr,'m',vol,'vol',sengc,'seng',Ev,'Ev',temp,'temp',dv,'dv'
  return
end function Qvf21d

  
function Qif1d(x,D,vol,rcap)
  use kindmod
  implicit none
  real(kind=dp),parameter :: Pi = 3.1415927_dp
  real(kind=dp) :: Qif1d, w, D, vol,rcap
  real(kind=dp) :: x, r, Zi_void
  r=(3._dp*x*vol/(4._dp*Pi))**(1._dp/3._dp)
  Zi_void= 1._dp+rcap/r
  Qif1d=4.*Pi*r*Zi_void*D/vol
  !print*, D,'d',w,'w'
  return
end function Qif1d
  

function Pif(x,b,D,Z,vol)
  use kindmod
  implicit none
  real(kind=dp),parameter :: Pi = 3.1415927_dp
  real(kind=dp) Pif, w, D, Z, vol,b
  real(kind=dp) :: x  
  w= (4.0_dp*Pi/(vol*b))**(0.5_dp)
  Pif=Z*w*D*x**(0.5_dp)
  return
end function Pif

function Jx(i)
  use nuclearvars
  use onedimensionalvars
  implicit none  
  integer :: i
  real(kind=dp) :: Jx
  interface 
    function linf1d(L0,Lx,x,meanx)
      use kindmod
      implicit none 
      real(kind=dp) :: linf1d,meanx,x,L0,Lx
    end function
  end interface 
  Jx=Pvarr(i)*Cv*linf1d( L0x(i),L1x(i),rg1xi(i),meanxarr(i)) - &
& (Qiarr(i+1)*Ci+Qvarr(i+1)+Qg(i+1)*kgCg)*linf1d(L0x(i+1),L1x(i+1),&
& rg1xi(i)+1._dp,meanxarr(i+1))  
  return
end function 

function Ji(i)
  use nuclearvars
  use onedimensionalvars
  implicit none  
  integer :: i
  real(kind=dp) :: Ji
  interface 
    function linf1d(L0,Lx,x,meanx)
      use kindmod
      implicit none 
      real(kind=dp) :: linf1d,meanx,x,L0,Lx
    end function
  end interface 
  Ji=Piarr(i)*Ci*linf1d( L0i(i),L1i(i),rg1ii(i),meaniarr(i)) - &
& Qi(i+1)*Cv*linf1d(L0i(i+1),L1i(i+1),rg1ii(i)+1._dp,meaniarr(i+1))  
  !if (i.eq. 5) then
  !  print*, Ji,Cv,Ci,L0i(i),linf1d( L0i(i),L1i(i),rg1ii(i),meaniarr(i))
  !  print*,L0i(i+1),L1i(i+1),linf1d(L0i(i+1),L1i(i+1),rg1ii(i)+1._dp,meaniarr(i+1)) 
  !end if 
  return
end function 


function JgL0(i)
  use nuclearvars
  use onedimensionalvars
  implicit none
  integer :: i
  real(kind=dp) :: JgL0,rXg,meanXg,meanXgp1,iterm,fgi,fgip1, &
& fgimeanx
  interface 
    function linf1d(L0,Lx,x,meanx)
      use kindmod
      implicit none 
      real(kind=dp) :: linf1d,meanx,x,L0,Lx
    end function
  end interface 
  kgCg=Cg*sqrt(kgsq) 
  rXg=real(Xg)
  meanXg = rg1xi(i) + 0.5_dp*(rXg-1._dp)
  meanXgp1= rg1xi(i+1) + 0.5_dp*(rXg-1._dp)

  fgi=linf1d(L0x(i),L1x(i),meanXg,meanxarr(i))
  fgip1=linf1d(L0x(i+1),L1x(i+1),meanXgp1,meanxarr(i+1))
  JgL0= rXg*(Qg(i)*fgi-Qg(i+1)*fgip1)*kgCg
  return
end function

function JgL1(i)
  use nuclearvars
  use onedimensionalvars
  implicit none
  integer :: i
  real(kind=dp) :: JgL1,rXg,meanXg,meanXgp1,iterm,fgi,fgip1, &
& fgimeanx
  interface 
    function linf1d(L0,Lx,x,meanx)
      use kindmod
      implicit none 
      real(kind=dp) :: linf1d,meanx,x,L0,Lx
    end function
    function sigmasq(dx)
      use kindmod
      implicit none
      real(kind=dp) :: dx
      real(kind=dp) :: sigmasq
    end function
  end interface 
  kgCg=Cg*sqrt(kgsq) 
  rXg=real(Xg)
  meanXg = rg1xi(i) + 0.5_dp*(rXg-1._dp)
  meanXgp1= rg1xi(i+1) + 0.5_dp*(rXg-1._dp)
  iterm=(meanXg-1_dp)/(rgx(i)*sigmasq(rgx(i)) )
  if (rgx(i) .eq. 1.0) iterm = 0.0
  fgi=linf1d(L0x(i),L1x(i),meanXg,meanxarr(i))
  fgip1=linf1d(L0x(i+1),L1x(i+1),meanXgp1,meanxarr(i+1))
  fgimeanx= linf1d(L0x(i),L1x(i),meanxarr(i)+rXg,meanxarr(i))
!print*, iterm,meanxg,sigmasq(rgx(i)),rgx(i),fgi,fgip1,fgimeanx
  JgL1= iterm*rXg*((Qg(i)*fgi-Qg(i+1)*fgip1)*kgCg + rgx(i) &
& * (fgimeanx-Qg(i+1)*fgip1)*kgCg)

  return
end function

function Jis(i)
  use nuclearvars
  use onedimensionalvars
  implicit none  
  integer :: i
  real(kind=dp) :: Jis,ax

  interface 
    function linf1d(L0,Lx,x,meanx)
      use kindmod
      implicit none 
      real(kind=dp) :: linf1d,meanx,x,L0,Lx
    end function
  end interface 
  !Jis=Piarr(i)*Ci*linf1d(L0i(i),L1i(i),real(gi(i)),0.5_dp) - &
!& Qi(i+1)*Cv*linf1d(L0i(i+1),L1i(i+1),0.5_dp*(1._dp+real(gi(i))),0.5_dp*real(gi(i+1)))  
  Jis=Piarr(i)*Ci*(L0i(i)+L1i(i)*(-0.5_dp)) &
& -Qi(i+1)*Cv*(L0i(i+1)+L1i(i+1)*(0.5_dp*(1._dp-rgi(i)-rgi(i+1))))
  return
end function 

!function Jis(i)
!  use nuclearvars
!  implicit none  
!  integer :: i
!  real(kind=dp) :: Jis,linf1d,ax
!  Jis=Piarr(i)*Ci*linf1d(L0i(i),L1i(i),0.0_dp,0.5_dp) - &
!& Qi(i+1)*Cv*linf1d(L0i(i+1),L1i(i+1),0.5_dp,0.0_dp)  
!  return
!end function 

function Jxs(i)
  use nuclearvars
  use onedimensionalvars
  implicit none  
  integer :: i
  real(kind=dp) :: Jxs,linf1d,ax,rgw,rgwp1
  Jxs=Pvarr(i)*Cv*(L0x(i)+L1x(i)*(-0.5_dp)) - &
& (Qiarr(i+1)*Ci+Qvarr(i+1)+Qg(i+1)*kgCg)*(L0x(i+1)+L1x(i+1)*&
& (0.5_dp*(1._dp-rgx(i)-rgx(i+1))))   
  return
end function 

function linf1d(L0,Lx,x,meanx)
  use kindmod
  implicit none 
  real(kind=dp) :: linf1d,meanx,x,L0,Lx
  linf1d = L0+Lx*(x-meanx)
  return
end function


