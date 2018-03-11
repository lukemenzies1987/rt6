subroutine funcallngrpss(n, t, y, fd, ires)
  use nuclearvars
  use twodimensionalvars
  implicit none
  external speciesgenerationHe
  

  real(kind=dp) :: siaabsorption, intabsorption, vacabsorption, Heemission, &
&  Heabsorption, vacemission, Hemonomeremission, Hevacabsorption,ss3v,ss3i,ss3He, &
& sgbv,sgbi,sgbHe,sumf
  real(kind=nag_wp), intent(in) :: t
  integer , intent(in) :: n
  integer, intent(inout) :: ires
  integer :: i,j,k
  real(kind=nag_wp), intent(in) :: y(n)
  real(kind=nag_wp), intent(out) :: fd(n)
  real(kind=dp) :: invgx,of,f1,invgm,jxi,jxf,jmi,jmf,jx,jm,linf,Jaltx,Jaltm, &
& f10,f11,f20,Gv2,Gv3,Gv4,Jaltxs,Jaltms,Jaltxb,Jaltmb,dL0,dLx,dLm
  real(kind=dp) :: edge, rdislc,rdisle,recombine,P0
  real(kind=dp) :: sv_dis0,si_dis0,tot_sv,tot_si,sv_void,si_void,ss01v,ss01i, &
& sv_void0,si_void0,tot_sv0,tot_si0


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
  !call sinkstrcorrngrp(y)
  
  do k=5+2,neq
    i=idarrayi(k)
    j=idarrayj(k)
    f(i,j)=abs(y(k))
  end do 
  sv_dis0=Zv*rhod
  si_dis0=Zi*rhod
  Cv=y(1)
  Ci=y(2)
  CHe=y(3)
  CGb=y(4)
  Cg=y(5)
  f(1,0)=Cv
  !Qm(1,1)=Qm(1,1)+ nuHe*Di*Ci*f(1,1)
  !Heemission=Heemission+fQiarr(1,1)*f(1,1)*Ci


  call speciesgenerationHe(GHe,t)
  !Qg=Qg*2._dp*Cg*sqrt(kgsq/2._dp)
  kgCg=2.*Cg*kgsq

  !k=6
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
    dL0=jxi-jxf+jmi-jmf
    !dL0= ( Pvarr(i-1,j)*Cv*f(i-1,j)-(Qiarr(i,j)*Ci+Qvarr(i,j)+Qg(i,j))*f(i,j)-&
!& Pvarr(i,j)*Cv*f(i,j)+(Qiarr(i+1,j)*Ci+Qvarr(i+1,j)+Qg(i+1,j))*f(i+1,j)) &
!& +(PHearr(i,j-1)*CHe*f(i,j-1)-Qm(i,j)*f(i,j)-PHearr(i,j)*CHe*f(i,j)+Qm(i,j+1)*f(i,j+1))

    siaabsorption=siaabsorption+fQiarr(i,j)*Ci*f(i,j)
    vacabsorption=vacabsorption+fPvarr(i,j)*Cv*f(i,j)
    Heemission=Heemission+fQm(i,j)*f(i,j)
    Heabsorption=Heabsorption+fPHearr(i,j)*CHe*f(i,j)
    vacemission=vacemission+fQvarr(i,j)*f(i,j)
		sigmavNv= sigmavNv + fQg(i,j)*f(i,j)

    fd(k)=dL0

  end do
  sigmavNv=sigmavNV/(Dg) 
  !do j=1,m0
  !  f1=f(1,j)    
    !Hevacabsorption=Hevacabsorption+PHearr(1,j)*CHe*f1
    !intabsorption=intabsorption+Qiarr(1,j)*Ci*f1
    !Hemonomeremission=Hemonomeremission+Qm(1,j)*f1
	!  sigmavNv= sigmavNv + Qg(1)*f(1,j)
    !print*, f1,j,Hemonomeremission,'f1',y(11),counter
    !if (y(3) .gt. 0.0) print*, y(3), counter,'che' 
  !end do

 ! do j=0,mnumg
    !siaabsorption=siaabsorption-Qiarr(1,j)*Ci*f(1,j)
    !Heabsorption=Heabsorption-PHearr(1,j)*CHe*f(1,j)
    !vacemission=vacemission-Qvarr(1,j)*f(1,j)
    !vacemission=vacemission-Qvarr(2,j)*f(2,j)
	!	sigmavNv= sigmavNv - Qg(1)*f(1,j)
  !end do

  !21 format (3(I3),'check',1pe16.7, 4(1pe16.7,I3))

  !do i=2,xnumg
    !Heemission=Heemission-Qm(i,0)*f(i,0)
  !end do 


  !Hevacabsorption=Hevacabsorption-PHearr(1,m0)*CHe*f(1,m0)

  f10=f(1,0)
  f11=f(1,1)
  f20=f(2,0)

  recombine=nui*Di*Cv*Ci 
  P0=Pvarr(1,0)*f10*Cv+PHearr(1,0)*f10*CHe
  vacemission=vacemission+Qm(1,1)*f11 + (Qiarr(2,0)*Ci+Qvarr(2,0)+Qg(2,0)*kgcg)*f20
  Heabsorption=Heabsorption+fPHearr(1,0)*f10*CHe
  Heemission=Heemission+ Qg(1,1)*kgcg*f11+Qiarr(1,1)*f11*Ci
  tgdid=tgdi!/(rhod*vol**(2._dp/3._dp))

  disl=ZHe*DHe*CHe*rhod*(disl0/(1.E+04*tgdid+disl0))!(1._dp-tgdid)!

  sHe_dis=disl/(DHe*CHe)



  !dCv = Gv  - nui*Di*Cv*Ci - Dv*(Cv-Cv0)*sv_dis  -Dv*(Cv-Cv0)*sv_gb &
!& - 2.*Pvarr(1,0)*f10*Cv- vacabsorption - PHearr(1,0)*f10*CHe &
!& + Qiarr(1,1)*f11*Ci + (Qiarr(2,0)*Ci+Qvarr(2,0)+Qg(2,0)*kgcg)*f20 + vacemission


!  dCi = Geni - nui*Di*Ci*Cv-Di*Ci*si_dis -siaabsorption-intabsorption  - Di*Ci*si_gb
  !dCi = Gi - nui*Di*Ci*Cv-Dv*Ci*si_dis -siaabsorption-intabsorption

!  dCHe = GHe  + Heemission + Hemonomeremission  - &
!&  disl-PHearr(1,0)*f10*CHe+Qg(1,0)*kgcg*f10 -Heabsorption - Hevacabsorption -DHe*CHe*sHe_gb



  dCv = Gv  -recombine - Dv*(Cv-Cv0)*sv_dis  -Dv*(Cv-Cv0)*sv_gb & 
& - P0 - vacabsorption + vacemission

  dCi = Geni - recombine -Di*Ci*si_dis -siaabsorption  - Di*Ci*si_gb

  dCHe = GHe  + Heemission    -disl &! disl 
& -Heabsorption  -DHe*CHe*sHe_gb

  dCGb= DHe*CHe*sHe_gb

  dCg= Gg-Dg*Cg*kgsq
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
