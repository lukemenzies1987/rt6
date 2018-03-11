subroutine initialsetup1d
!------------------------------------------------------------------
! The initial setup subroutine for the 2 dimensional case.
!------------------------------------------------------------------    
  use nuclearvars
  use onedimensionalvars
  implicit none    
  integer :: nunit,ios,i,j,checkst,arrcount
  real(kind=dp) :: Fedens,E2v,mr

  interface

    function Qvf1d(x,Dv,vol,E2v,Ev,temp)
      use kindmod
      implicit none
      real(kind=dp):: Qvf1d, Dv, vol,Ev,temp,x,E2v
    end function Qvf1d

    function Qvf21d(x,m,Dv,vol,seng,Ev,temp)
      use onedimensionalvars
      implicit none
      real(kind=dp):: D, Qvf21d, w, Dv, vol,Ev,seng,temp,Z
      real(kind=dp):: x,m
    end function Qvf21d

  
    function Qif1d(x,D,vol,rcap)
      use kindmod
      implicit none
      real(kind=dp) :: Qif1d, D, vol,rcap
      real(kind=dp) :: x
    end function Qif1d
      

    function Pif(x,b,D,Z,vol)
      use kindmod
      implicit none
      real(kind=dp) Pif, D, Z, vol,b
      real(kind=dp) :: x  
    end function Pif

    function diff(temp,Em,D0)
      use kindmod
      implicit none
      real(kind=dp),intent(in):: temp,Em,D0
      real(kind=dp) :: diff
    end function

    function Qgcoeff(Dg,vol) 
      use kindmod
      implicit none
      real(kind=dp),intent(in) :: Dg,vol
      real(kind=dp) :: Qgcoeff    
    end function Qgcoeff

    function ksqst(l,Rg)
      use kindmod
      implicit none
      real(kind=dp),intent (in) :: l,Rg
      real(kind=dp) :: ksqst 
    end function ksqst

    function Pf(x,D,vol)
      use kindmod
      implicit none
      real(kind=dp),intent(in) :: x, D, vol
      real(kind=dp) :: Pf
    end function Pf

    function Qif(x,D,vol,rcap)
      use kindmod
      implicit none
      real(kind=dp), intent(in) :: x, D, vol,rcap
      real(kind=dp) :: Qif
    end function Qif

    function Qvf4(x,m,i,j)
      use nuclearvars
      use twodimensionalvars
      implicit none
      real(kind=dp):: Qvf4
      real(kind=dp), intent(in) :: x, m
      integer, intent(in) :: i,j
    end function

    function ksqft(rhod,disabs)
      use kindmod
      implicit none
      real(kind=dp), intent(in) :: rhod,disabs
      real(kind=dp) :: ksqft
    end function ksqft


    function Qvf3(x,m,Dv,vol,seng,Ev,temp,i,j)
      use twodimensionalvars
      implicit none
      real(kind=dp), intent(in) :: x, m,Dv,vol,seng,Ev,temp
      real(kind=dp) :: Qvf3
      integer, intent(in) :: i,j
    end function

    subroutine bindingengtest2
    end subroutine

    subroutine groupingsetup1d
    end subroutine

  end interface

  b= a*1.E-10*dsqrt(3._dp)/2._dp !burgers vector

  !Atomic volume determined by the lattice paramter of bcc iron. 
  !if one wishes to set this manually in the input parameter file, just comment out vol
  vol=((a*1.E-10)**3._dp)/2._dp

  !Controls whether to partition the sdf grid using the 
  !grouping method or not. 
  if( flaggrouping .eq. 1) then
    print*, 'Setting up grouping partitions'
    call groupingsetup1d
    print*,'done!'
  else
    xnumg=xnum
    srtgrpx=xnumg
    arrcount=1
    do i=1,xnumg
      arrcount=arrcount+1
    end do
    neq2=arrcount-1
    neq=neq2+2
    allocate(idarrayi(neq2),stat=checkst)
    if(checkst.ne.0) then
      write(6,*) 'Problem allocating arrays'
      stop
    end if
    arrcount=1
    do i=1,xnumg      
      idarrayi(arrcount)=i
      arrcount=arrcount+1
    end do 

  end if 
  !Allocation of arrays in 1d.
  !neq=xnum*mnum+xnum
  allocate(Qiarr(0:xnumg+1),Pvarr(0:xnumg+1),Qvarr(0:xnumg+1),fx(0:srtgrpx+Xg), &
  & fi(0:srtgrpi+1),Piarr(0:inumg+1),Qi(0:inumg+1),Qg(0:xnumg+1),bindingV(xnumg),& 
  & atol(neq),r(xnumg),ril(inumg),Genv(xnumg),Geniloop(inumg),stat=checkst)
  if(checkst.ne.0) then
    write(6,*) 'Problem allocating arrays'
    stop
  end if

  Dv = diff(temp,Emv,D0v) ! Vacancy diffusion coefficient
  Di = diff(temp,Emi,D0i) ! SIA diffusion coefficient
  Dg = diff(temp,Emi,D0i) ! Glissile SIA diffusion coefficient

  Qgc=Qgcoeff(Dg,vol) ! Glissile SIA capture coefficient term

  Cv0=Cv0*dexp(-Ev/(kB*temp)) !Initial vacancy concentration (equilibrium concentration).
  Qiarr=0.0_dp
  Pvarr=0.0_dp
  Qvarr=0.0_dp
  Piarr=0.0_dp
  Qi=0.0_dp
  E2v=0.2_dp !di-vacancy binding energy

  !Setup of capture and emission rates
  !--------------------------------------
  do i=1,xnumg-1
    mr=rg1xi(i)-0.5_dp*rgx(i)
    if (i .le. srtgrpx) mr=rg1xi(i)
    Pvarr(i)= Pf(mr,Dv,vol) ! Vacancy capture rates
  end do

  do i=1,xnumg
    mr=rg1xi(i)-0.5_dp*rgx(i)
    if (i .le. srtgrpx) mr=rg1xi(i)
    Qvarr(i)= Qvf1d(mr,Dv,vol,E2v,Ev,temp) ! Vacancy emission rates (1d). 
    !Qvarr(i)=Qvf21d(mr,0.0_dp,Dv,vol,seng,Ev,temp)
    Qiarr(i)= Qif1d(mr,Di,vol,rcap)!Qif(mr,Di,vol,0.0_dp) !SIA capture rates (1d). 
    !print*, Qiarr(i),Pvarr(i),Qvarr(i),i
  end do

  !if (flagEbtable.eq.1) call setPQ
  !poisson=0.33
  !shear= 3.24E29
  do i=1,inumg-1
    !Piarr(i)=Pif(rg1ii(i),b,Dv,Di,Zvl,Zil,shear,poisson,seng,Ev,temp,vol)
    mr=rg1ii(i)-0.5_dp*rgi(i)
    if (i .le. srtgrpi) mr=rg1ii(i)
    !Piarr(i)=Pif(mr,b,Di,Zil,vol)
  end do

  do i=1,inumg
    mr=rg1ii(i)-0.5_dp*rgi(i)
    if (i .le. srtgrpi) mr=rg1ii(i)
    !Qi(i)=Pif(mr,b,Dv,Zvl,vol)
  !print*, Piarr(i),Qi(i),rg1ii(i)
  end do 

  !Sink strength terms for SIA glissile clusters. Set to zero for
  ! no SIA glissile contribution. 
  fstterm=ksqft(rhod,disabs) !dislocations
  secndterm=ksqst(l,Rg) !Grain boundaries

  ! Small SIA glissile capture rates
  do i=1,xnumg
    mr=rg1xi(i)-0.5_dp*rgx(i)
    if (i .le. srtgrpx) mr=rg1xi(i)
	  Qg(i)=mr**(2._dp/3._dp)*Qgc
  end do 

  
  ! Void radius
  do i=1,xnumg
    r(i) = (3._dp*meanxarr(i)*vol/(4._dp*Pi))**(1._dp/3._dp)        
  end do 

  !Void radius interstitial loop (not used)
  do i=1,inumg
    ril(i) = (meaniarr(i)*vol/(b*Pi))**(1._dp/2._dp)        
  !print*, r(i)
  end do 

  !Binding energy test (not used)
  !call bindingengtest(Dv,vol,E2v,Ev,temp)  
  !call bindingengtest2

  !Initial Conditions
  Cv= Cv0 !Initial vacancy concentration (equilibrium concentration).
  Ci= Ci0 !Initial SIA concentration.
  Cg=0.0
  fx=0.0_dp
  fi=0.0_dp
  fx(1)=Cv
  fi(1)=Ci
  sv_dis=rhod*Zv  !sink strength for vacancies to dislocaitons. 
  si_dis=rhod*Zi  !sink strength for SIA to dislocaitons. 
  return
end subroutine initialsetup1d

