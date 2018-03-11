subroutine funcall1d(n, t, y, fd, ires)
!------------------------------------------------------------------
! Rate equation function, 1 dimensional
!------------------------------------------------------------------ 
  use nuclearvars
  use onedimensionalvars
  implicit none
  

  real(kind=dp)::siaabsorption, intabsorption, vacabsorption,vacemission,& 
& intvoidabs, intemission, glissemiss,JgL0,JgL1
  real(kind=nag_wp), intent(in) :: t
  integer , intent(in) :: n
  integer, intent(inout) :: ires
  integer :: i,j,k
  real(kind=nag_wp), intent(in) :: y(n)
  real(kind=nag_wp), intent(out) :: fd(n)
  real(kind=dp) :: invgx,of,f1,invgi,jxi,jxf,ji,linf,Jis, &
& Gv2,Gv3,Gv4,jii,jif,jp,sgbv,sgbi,ss3v,ss3i
  
  interface
    function gbss( sss,rgr )
      use kindmod
      implicit none
      real(kind=dp) :: sss,rgr
      real(kind=dp) :: gbss
    end function


    function Jx(i)
      use nuclearvars
      use onedimensionalvars
      implicit none  
      integer :: i
      real(kind=dp) :: Jx
    end function

    function Jxs(i)
      use nuclearvars
      use onedimensionalvars
      implicit none  
      integer :: i
      real(kind=dp) :: Jxs
    end function

  end interface

  intvoidabs=0.0
  intabsorption=0.0
  intemission=0.0
  vacabsorption=0.0
  vacemission=0.0
  glissemiss=0.0
  sigmavNv=0.0_dp

  do k=1,neq1
    i=idarrayi(k)
    fx(i)=abs(y(k))
!print*,k,i
  end do 

  do k=neq1+1,neq2,2
    i=idarrayi(k)
    L0x(i)=abs(y(k))
    L1x(i)=y(k+1)
!print*,k,i,neq2,k+1
  end do 


  !do k=neq2+1,neq3
  !  i=idarrayi(k)
  !  fi(i)=abs(y(k))
!print*,k,i
  !end do 

  !do k=neq3+1,neq4,2
  !  i=idarrayi(k)
  !  L0i(i)=abs(y(k))
  !  L1i(i)=y(k+1)
!print*,k,i
  !end do 
!print*, neq,'neq',xnumg,inumg,neq2+1
  fx(srtgrpx+1)=L0x(srtgrpx+1)+L1x(srtgrpx+1)*(rg1xi(srtgrpx)+1.-meanxarr(srtgrpx+1))
  L0x(srtgrpx)=fx(srtgrpx)
  L1x(srtgrpx)=0.0

  !fi(srtgrpi+1)=L0i(srtgrpi+1)+L1i(srtgrpi+1)*(rg1ii(srtgrpi)+1.-meaniarr(srtgrpi+1))
  !L0i(srtgrpi)=fi(srtgrpi)
  !L1i(srtgrpi)=0.0
  !do i=1, Xg
    !if Xg .gt. gx(srtgrpx+1) then
    ! fx(srtgrpx+i)=L0x(srtgrpx+1)+L1x(srtgrpx+1)*(rg1xi(srtgrpx)+i-meanxarr(srtgrpx+1))
    !end if 
  !  fx(srtgrpx+i)=L0x(srtgrpx+i)+L1x(srtgrpx+i)*(rg1xi(srtgrpx+i)-meanxarr(srtgrpx+i))
  !end do
  Cv=y(1) !Vacancy concentration 
  Ci=y(neq2+1) !SIA concentration
  Cg=y(neq2+2) !SIA glissile concentration
  !fx(1)=Cv
  !fi(1)=Ci
  !y(1)=Cv
  !Qg=Qg!*Cg*sqrt(kgsq)
  kgCg=Cg*dsqrt(kgsq/2._dp) !Glissile sink strength
  k=1
!sdf for x--------------------------------------
  k=k+1
  do i=2,srtgrpx
    jxi=Pvarr(i-1)*Cv*fx(i-1)-(Qiarr(i)*Ci+Qvarr(i)+Qg(i)*kgCg)*fx(i)
    jxf=Pvarr(i)*Cv*fx(i)-(Qiarr(i+1)*Ci+Qvarr(i+1))*fx(i+1)-Qg(i+1)*kgCg*fx(i+1)
    ydot(k)= Genv(i)+jxi-jxf 

    vacabsorption=vacabsorption+Pvarr(i)*Cv*fx(i)
    vacemission=vacemission+Qvarr(i)*fx(i)  
    intvoidabs=intvoidabs+Qiarr(i)*Ci*fx(i)
		sigmavNv= sigmavNv + Qg(i)*fx(i)*kgCg

    k=k+1
  end do
  
  vacabsorption=vacabsorption+Pvarr(1)*Cv*fx(1)
  vacemission=vacemission-Qvarr(2)*fx(2)
  
  do i=srtgrpx+1,xnumg   
    invgx=(1._dp/rgx(i))
    jxi=Jx(i-1)
    jxf=Jx(i)
    ydot(k) = invgx*(Jxi-Jxf)
    ydot(k+1) = -itermx(i)*(Jxi+Jxf- 2._dp*Jxs(i))

    vacabsorption=vacabsorption+Pvarr(i)*Cv*L0x(i)*rgx(i)
    vacemission=vacemission+Qvarr(i)*L0x(i)*rgx(i)
    intvoidabs=intvoidabs+Qiarr(i)*Ci*L0x(i)*rgx(i)
		sigmavNv= sigmavNv + Qg(i)*L0x(i)*rgx(i)*kgCg
    !print*,k,i,'t'
    !print*,k+1,i,'t'
    k=k+2
  end do 
  !print*,neq2,k,k+Xg-1,'ks'
  !print*, k,'ci',neq2+1
!sdf for i--------------------------------------
  !k=k+Xg-1
!print*, k,'k1'
  
 !   jii=-Qi(Xg)*Cv*fi(Xg)
 !   jif=Piarr(Xg)*Ci*fi(Xg)-Qi(Xg+1)*Cv*fi(Xg+1)
 !   ydot(k)=Geni(Xg)+ jii-jif
!print*,ydot(k),Geni(i),jii,jif,'sia',i
 !   intabsorption=intabsorption+Piarr(Xg)*Ci*fi(Xg)
 !   intemission=intemission+Qi(Xg)*Cv*fi(Xg)
 ! k=k+1

 ! do i=Xg+1,srtgrpi
!print*, k,'k'
 !   jii=Piarr(i-1)*Ci*fi(i-1)-Qi(i)*Cv*fi(i)
 !   jp=0.0
 !   if (rg1xi(i) .gt. Xg+1 .and. rg1xi(i) .le. 2.*Xg) jp=Qg(i)*kgCg*fi(i)
 !   jif=Piarr(i)*Ci*fi(i)-Qi(i+1)*Cv*fi(i+1)
 !   ydot(k)=Geni(i)+ jii-jif-jp
!print*,ydot(k),Geni(i),jii,jif,'sia',i
 !   intabsorption=intabsorption+Piarr(i)*Ci*fi(i)
 !   intemission=intemission+Qi(i)*Cv*fi(i)
    !intemission=intemission+Qi(i)*Cv*fi(i)
   ! print*,k,i,'t'
 !   k=k+1
 ! end do
  !intabsorption=intabsorption+Piarr(1)*Ci*fi(1)
 ! do i=srtgrpi+1,inumg   
 !   invgi=(1._dp/rgi(i))
 !   jii=Ji(i-1)
 !   jif=Ji(i)      
    !if (i-1 .eq.srtgrpi) print*, jii,'juu' 
 !   ydot(k) = invgi*(Jii-Jif) 
 !   ydot(k+1) = -itermi(i)*(Jii+Jif- 2._dp*Jis(i))
 !   intabsorption=intabsorption+Piarr(i)*Ci*L0i(i)*gi(i)
 !   intemission=intemission+Qi(i)*Cv*L0i(i)*gi(i)
    !print*,k,i,'t'
    !print*,k+1,i,'t'
 !   k=k+2
 ! end do 

  !do i=2,Xg-2
  !  glissemiss=glissemiss+((real(Xg)-rg1xi(i))/real(Xg))*Qg(i)*kgCg*fx(i)
  !end do 
!  dCv = Genv(1)-nui*Di*Cv*Ci -Dv*(Cv-Cv0)*sv_dis &!-Dv*(Cv-Cv0)*sv_gb &
!& -(2._dp*Pvarr(1)*Cv*fx(1)-2._dp*(Qiarr(2)*Ci+Qvarr(2))*fx(2)) &
!& -(vacabsorption-vacemission)+Qg(Xg+1)*kgCg*fx(Xg+1)!-intemission

  if (Ci .eq. 0.0) then 
    ss3v =sv_dis
    ss3i=si_dis 
    sigmavNv= sigmavNv
  else
    ss3v =sv_dis +(vacabsorption+2._dp*Pvarr(1)*fx(1)*Cv)/(Dv*Cv)
    ss3i=si_dis + intvoidabs/(Di*Ci)
    sigmavNv= sigmavNv/(Dg*Cg)
  end if 

!print*, ss3v,l
  sgbv = gbss(ss3v,l)
  sv_gb=sgbv


  sgbi = gbss(ss3i,l)
  si_gb= sgbi

  dCv = Genv(1)-nui*Di*Cv*Ci  -Dv*(Cv-Cv0)*sv_dis  -Dv*(Cv-Cv0)*sv_gb &
& -(Pvarr(1)*Cv*fx(1)-(Qiarr(2)*Ci+Qvarr(2)+Qg(2)*kgCg)*fx(2)) &
& -(vacabsorption-vacemission)!+Qg(Xg+1)*kgCg*fx(Xg+1)!-intabsorption


!  dCi = Geni(1)-nui*Di*Ci*Cv-Dv*Ci*si_dis &! -Di*Ci*si_gb &
!& -(2._dp*Piarr(1)*Ci*fi(1))-(intabsorption) &
!& -intvoidabs

  dCi = Geni-nui*Di*Ci*Cv  -Di*Ci*si_dis -Di*Ci*si_gb &
& &! -(Piarr(1)*Ci*fi(1)-Qi(2)*Cv*fi(2))&! -(intabsorption) &
& -intvoidabs !- Piarr(1)*Ci*fi(1)+Qg(Xg-1)*kgCg*fx(Xg-1)

  kgsq=2._dp*(fstterm+secndterm+sigmavNv)*(fstterm+secndterm+sigmavNv)

  dCg= Gg-Dg*Cg*kgsq!+glissemiss


  ydot(1)=dCv
  ydot(neq2+1)=dCi
  ydot(neq2+2)=dCg

  return
end subroutine
