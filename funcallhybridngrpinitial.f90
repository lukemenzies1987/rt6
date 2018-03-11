subroutine funcallngrpinit(n, t, y, ydot)
  use nuclearvars
  use twodimensionalvars
  implicit none
  external speciesgenerationHe
  

  real(kind=dp) :: t, siaabsorption, intabsorption, vacabsorption, Heemission, &
&  Heabsorption, vacemission, Hemonomeremission, Hevacabsorption,ss3v,ss3i, &
& sgbv,sgbi,gbss,sumf
  integer :: i,j,k,n,ii,ji
  real(kind=dp),dimension(n) :: y, ydot
  real(kind=dp) :: invgx,of,f1,invgm,jxi,jxf,jmi,jmf,jx,jm,linf,Jaltx,Jaltm, &
& f10,f11,f20,Gv2,Gv3,Gv4,Jaltxs,Jaltms,Jaltxb,Jaltmb,dL0,dLx,dLm
  real(kind=dp) :: edge, rdislc,rdisle

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


  do k=5+1,neq
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
  Qm(1,1)=Qm(1,1)+ nuHe*Di*Ci*f(1,1)
  call speciesgenerationHe(GHe,t)
  Qg=Qg*Cg*sqrt(kgsq/2._dp)
  !k=6
!first half--------------------------------------
  do k=5+1,neq

    i=idarrayi(k)
    j=idarrayj(k)
!print*, k,i,j
    dL0= ( Pvarr(i-1,j)*Cv*f(i-1,j)-(Qiarr(i,j)*Ci+Qvarr(i,j)+Qg(i))*f(i,j)-&
& Pvarr(i,j)*Cv*f(i,j)+(Qiarr(i+1,j)*Ci+Qvarr(i+1,j)+Qg(i+1))*f(i+1,j)) &
& +(PHearr(i,j-1)*CHe*f(i,j-1)-Qm(i,j)*f(i,j)-PHearr(i,j)*CHe*f(i,j)+Qm(i,j+1)*f(i,j+1))

    siaabsorption=siaabsorption+Qiarr(i,j)*Ci*f(i,j)
    vacabsorption=vacabsorption+Pvarr(i,j)*Cv*f(i,j)
    Heemission=Heemission+Qm(i,j)*f(i,j)
    Heabsorption=Heabsorption+PHearr(i,j)*CHe*f(i,j)
    vacemission=vacemission+Qvarr(i,j)*f(i,j)
		sigmavNv= sigmavNv + Qg(i)*f(i,j)

    ydot(k)=dL0
  end do
   
  do j=1,m0
    f1=f(1,j)    
    Hevacabsorption=Hevacabsorption+PHearr(1,j)*CHe*f1
    intabsorption=intabsorption+Qiarr(1,j)*Ci*f1
    Hemonomeremission=Hemonomeremission+(Qm(1,j))*f1
	  sigmavNv= sigmavNv + Qg(1)*f(1,j)
    !print*, f1,j,Hemonomeremission,'f1',y(11),counter
    !if (y(3) .gt. 0.0) print*, y(3), counter,'che' 
  end do

  do j=0,mnumg
    siaabsorption=siaabsorption-Qiarr(1,j)*Ci*f(1,j)
    Heabsorption=Heabsorption-PHearr(1,j)*CHe*f(1,j)
    vacemission=vacemission-Qvarr(1,j)*f(1,j)
    vacemission=vacemission-Qvarr(2,j)*f(2,j)
		sigmavNv= sigmavNv - Qg(1)*f(1,j)
  end do



  do i=2,xnumg
    Heemission=Heemission-Qm(i,0)*f(i,0)
  end do 


  Hevacabsorption=Hevacabsorption-PHearr(1,m0)*CHe*f(1,m0)

  f10=f(1,0)
  f11=f(1,1)
  f20=f(2,0)
 
  ss3v =sv_dis
  ss3i=si_dis 
  disl=DHe*CHe*ZHe*rhod
  


  sgbv = gbss(ss3v,l)
  sv_gb=sgbv


  sgbi = gbss(ss3i,l)
  si_gb= sgbi


 
  kgsq=2._dp*(fstterm+secndterm)*(fstterm+secndterm)


  dCv = Gv - nui*Di*Cv*Ci - nuHe*DHe*CHe*Cv - Dv*(Cv-Cv0)*sv_dis  -Dv*(Cv-Cv0)*sv_gb &
& - 2.*Pvarr(1,0)*f10*Cv- vacabsorption - PHearr(1,0)*f10*CHe &
& + Qm(1,1)*f11 + (Qiarr(2,0)*Ci+Qvarr(2,0)+Qg(2))*f20 + vacemission


  dCi = Geni - nui*Di*Ci*Cv-Dv*Ci*si_dis -siaabsorption-intabsorption  - Di*Ci*si_gb


  dCHe = GHe + nuHe*Di*Ci*f11 + Heemission + Hemonomeremission - nuHe*DHe*CHe*Cv - &
& disl  - Heabsorption - Hevacabsorption -DHe*CHe*ksqGb


  dCGb= DHe*CHe*ksqGb

  dCg= Gg-Dg*Cg*kgsq

  do i=2, maxgv
    ydot(5+i)=ydot(5+i)+genv(i)
  end do 
  
  ydot(1)=dCv
  ydot(2)=dCi
  ydot(3)=dCHe
  ydot(4)=dCGb
  ydot(5)=dCg
  ydot(6)=dCv
  Qm(1,1)=Qm(1,1)- nuHe*Di*Ci*f(1,1)
  return
end subroutine
