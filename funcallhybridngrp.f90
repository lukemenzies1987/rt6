subroutine funcallngrpinit(n, t, y, fd, ires)
!------------------------------------------------------------------
! Initial rate equation function, non-grouping method
!------------------------------------------------------------------ 
  use nuclearvars
  use twodimensionalvars
  implicit none

  

  real(kind=nag_wp) :: siaabsorption, intabsorption, vacabsorption, Heemission, &
&  Heabsorption, vacemission, Hemonomeremission, Hevacabsorption,ss3v,ss3i, &
& sgbv,sgbi,sumf
  real(kind=nag_wp), intent(in) :: t
  integer , intent(in) :: n
  integer, intent(inout) :: ires
  integer :: i,j,k
  real(kind=nag_wp), intent(in) :: y(n)
  real(kind=nag_wp), intent(out) :: fd(n)
  real(kind=nag_wp) :: invgx,of,f1,invgm,jxi,jxf,jmi,jmf,jx,jm,linf,Jaltx,Jaltm, &
& f10,f11,f20,Gv2,Gv3,Gv4,Jaltxs,Jaltms,Jaltxb,Jaltmb,dL0,dLx,dLm
  real(kind=nag_wp) :: edge, rdislc,rdisle,recombine,P0,ss3He,sgbHe


  interface
    function gbss( sss,rgr )
      use kindmod
      implicit none
      real(kind=dp),intent(in) :: sss,rgr
      real(kind=dp) :: gbss
    end function

    subroutine speciesgenerationHe(GHe,t)
      use kindmod
      implicit none
      real(kind=dp) :: GHe,t
    end subroutine

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


  do k=5+2,neq
    i=idarrayi(k)
    j=idarrayj(k)
    f(i,j)=abs(y(k))
  end do 

  Cv=y(1)
  Ci=y(2)
  CHe=y(3)
  CGb=y(4)
  Cg=y(5)
  f(1,0)=Cv
  kgsq=(fstterm+secndterm)
  kgCg=2._dp*Cg*kgsq
  !Qm(1,1)=Qm(1,1)+ nuHe*Di*Ci*f(1,1)
  call speciesgenerationHe(GHe,t) !He generation rate
  !Qg=Qg*2._dp*Cg*sqrt(kgsq/2._dp)
  !k=6
!first half--------------------------------------
  do k=5+2,neq

    i=idarrayi(k)
    j=idarrayj(k)
!print*, k,i,j
    jxi=Pvarr(i-1,j)*Cv*f(i-1,j)-(Qiarr(i,j)*Ci+Qvarr(i,j)+Qg(i,j)*kgcg)*f(i,j)
    jxf=Pvarr(i,j)*Cv*f(i,j)-(Qiarr(i+1,j)*Ci+Qvarr(i+1,j)+Qg(i+1,j)*kgcg)*f(i+1,j)
    jmi=PHearr(i,j-1)*CHe*f(i,j-1)-Qm(i,j)*f(i,j)
    jmf=PHearr(i,j)*CHe*f(i,j)-Qm(i,j+1)*f(i,j+1)
    dL0=jxi-jxf+jmi-jmf

    siaabsorption=siaabsorption+Qiarr(i,j)*Ci*f(i,j)
    vacabsorption=vacabsorption+Pvarr(i,j)*Cv*f(i,j)
    Heemission=Heemission+Qm(i,j)*f(i,j)
    Heabsorption=Heabsorption+PHearr(i,j)*CHe*f(i,j)
    vacemission=vacemission+Qvarr(i,j)*f(i,j)
		sigmavNv= sigmavNv + Qg(i,j)*f(i,j)

    fd(k)=dL0
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

subroutine funcallngrp(n, t, y, fd, ires)
!------------------------------------------------------------------
! Rate equation function, non-grouping method
!------------------------------------------------------------------ 
  use nuclearvars
  use twodimensionalvars
  implicit none
  
  real(kind=nag_wp) :: siaabsorption, intabsorption, vacabsorption, Heemission, &
&  Heabsorption, vacemission, Hemonomeremission, Hevacabsorption,ss3v,ss3i, &
& sgbv,sgbi,sumf
  real(kind=nag_wp), intent(in) :: t
  integer , intent(in) :: n
  integer, intent(inout) :: ires
  integer :: i,j,k,concgas
  real(kind=nag_wp), intent(in) :: y(n)
  real(kind=nag_wp), intent(out) :: fd(n)
  real(kind=nag_wp) :: invgx,of,f1,invgm,jxi,jxf,jmi,jmf,jx,jm,linf,Jaltx,Jaltm, &
& f10,f11,f20,Gv2,Gv3,Gv4,Jaltxs,Jaltms,Jaltxb,Jaltmb,dL0,dLx,dLm
  real(kind=nag_wp) :: edge, rdislc,rdisle,recombine,P0,ss3He,sgbHe



  interface
    function gbss( sss,rgr )
      use kindmod
      implicit none
      real(kind=dp) :: sss,rgr
      real(kind=dp) :: gbss
    end function

    subroutine speciesgenerationHe(GHe,t)
      use kindmod
      implicit none
      real(kind=dp) :: GHe,t
    end subroutine
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
  fd=0.0_dp
  !call sinkstrcorrngrp(y)
  
  do k=5+2,neq
    i=idarrayi(k)
    j=idarrayj(k)
    f(i,j)=abs(y(k))
  end do 

  Cv=y(1) !Vacancy concentration
  Ci=y(2) !SIA concentration
  CHe=y(3) !He concentration
  CGb=y(4) !He GB concentration
  Cg=y(5) !SIA glissile concentration

  f(1,0)=Cv
  !Qm(1,1)=Qm(1,1)+ nuHe*Di*Ci*f(1,1)
  !Heemission=Heemission+fQiarr(1,1)*f(1,1)*Ci

  call speciesgenerationHe(GHe,t) !He generation rate
  !Qg=Qg*2._dp*Cg*sqrt(kgsq/2._dp)
  kgCg=2._dp*Cg*kgsq

  fd(6)=0.0
!first half--------------------------------------
  do k=5+2,neq
    i=idarrayi(k)
    j=idarrayj(k)
   !write(274,21), k,i,j,f(i,j),f(i+1,j),i+1,f(i-1,j),i-1,f(i,j-1),j-1,f(i,j+1) ,j+1  
!print*, k,i,j
    jxi=Pvarr(i-1,j)*Cv*f(i-1,j)-(Qiarr(i,j)*Ci+Qvarr(i,j)+Qg(i,j)*kgcg)*f(i,j)
    jxf=Pvarr(i,j)*Cv*f(i,j)-(Qiarr(i+1,j)*Ci+Qvarr(i+1,j)+Qg(i+1,j)*kgcg)*f(i+1,j)
    jmi=PHearr(i,j-1)*CHe*f(i,j-1)-Qm(i,j)*f(i,j)
    jmf=PHearr(i,j)*CHe*f(i,j)-Qm(i,j+1)*f(i,j+1)
    !j1(i,j)=jxi!f(i-1,j)
    !j2(i,j)=jxf!f(i,j)
    !j3(i,j)=jmi
    !j4(i,j)=jmf
    dL0=jxi-jxf+jmi-jmf

    !dL0= ( Pvarr(i-1,j)*Cv*f(i-1,j)-(Qiarr(i,j)*Ci+Qvarr(i,j)+Qg(i,j))*f(i,j)-&
!& Pvarr(i,j)*Cv*f(i,j)+(Qiarr(i+1,j)*Ci+Qvarr(i+1,j)+Qg(i+1,j))*f(i+1,j)) &
!& +(PHearr(i,j-1)*CHe*f(i,j-1)-Qm(i,j)*f(i,j)-PHearr(i,j)*CHe*f(i,j)+Qm(i,j+1)*f(i,j+1))

    siaabsorption=siaabsorption+Qiarr(i,j)*Ci*f(i,j)
    vacabsorption=vacabsorption+Pvarr(i,j)*Cv*f(i,j)
    Heemission=Heemission+Qm(i,j)*f(i,j)
    Heabsorption=Heabsorption+PHearr(i,j)*CHe*f(i,j)
    vacemission=vacemission+Qvarr(i,j)*f(i,j)
		sigmavNv= sigmavNv + Qg(i,j)*f(i,j)
    !sdf(i,j)=dL0
    fd(k)=dL0

  end do
  sigmavNv=sigmavNV/Dg   


  f10=f(1,0)
  f11=f(1,1)
  f20=f(2,0)

  recombine=nui*Di*Cv*Ci 
  P0=Pvarr(1,0)*f10*Cv
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
& - P0 - vacabsorption + vacemission

  dCi = Geni - recombine -Di*Ci*si_dis -siaabsorption  - Di*Ci*si_gb

  dCHe = GHe  + Heemission  -disl &! disl 
& -Heabsorption  -DHe*CHe*sHe_gb


  dCGb= DHe*CHe*sHe_gb
  kgsq=(fstterm+secndterm+sigmavNv)
  dCg= Gg-Dg*Cg*2_dp*kgsq*kgsq-Qg(1,1)*kgcg*f11


  do i=2, maxgv
    fd(5+i)=fd(5+i)+genv(i)
  end do 
  !print*, Cv,dCv,DCi,dCHe
  fd(1)=dCv
  fd(2)=dCi
  fd(3)=dCHe
  fd(4)=dCGb
  fd(5)=dCg
  fd(6)=dCv



  
  return
end subroutine
