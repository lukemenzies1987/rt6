  subroutine initialsetup
!------------------------------------------------------------------
! The initial setup subroutine for the 2 dimensional case.
!------------------------------------------------------------------ 
  use nuclearvars
  use twodimensionalvars
  implicit none    
  integer :: nunit,ios,i,j,checkst,arrcount
  real(kind=dp) :: Qvf,Qvf2,Qvf5, &
& Fedens,EbHe,mr,mm,test

  interface
    function diff(temp,Em,D0)
      use kindmod
      implicit none
      real(kind=dp):: temp,Em,D0
      real(kind=dp) :: diff
    end function

    function Qgcoeff(Dg,vol) 
      use kindmod
      implicit none
      real(kind=dp) :: Dg,vol
      real(kind=dp) :: Qgcoeff    
    end function Qgcoeff

    function ksqst(l,Rg)
      use kindmod
      implicit none
      real(kind=dp):: l,Rg
      real(kind=dp) :: ksqst 
    end function ksqst

    function Pf(x,D,vol)
      use kindmod
      implicit none
      real(kind=dp):: x, D, vol
      real(kind=dp) :: Pf
    end function Pf

    function Qif(x,D,vol,rcap)
      use kindmod
      implicit none
      real(kind=dp):: x, D, vol,rcap
      real(kind=dp) :: Qif
    end function Qif

    function Qvf4(x,m,i,j)
      use nuclearvars
      use twodimensionalvars
      implicit none
      real(kind=dp):: Qvf4
      real(kind=dp):: x, m
      integer :: i,j
    end function

    function ksqft(rhod,disabs)
      use kindmod
      implicit none
      real(kind=dp):: rhod,disabs
      real(kind=dp) :: ksqft
    end function ksqft

    function Qemiss(x,m,DHe,vol,EbHe,temp,i,j)
      use twodimensionalvars
      implicit none
      real(kind=dp) :: x, m,DHe,vol,EbHe,temp
      real(kind=dp) :: Qemiss
      integer :: i,j
    end function

    function Qemissalt(x,m,DHe,vol,temp,i,j)
      use twodimensionalvars
      implicit none
      real(kind=dp) :: x, m,DHe,vol,temp
      real(kind=dp) :: Qemissalt
      integer :: i,j
    end function

    function Qvf3(x,m,Dv,vol,seng,Ev,temp,i,j)
      use twodimensionalvars
      implicit none
      real(kind=dp):: x, m,Dv,vol,seng,Ev,temp
      real(kind=dp) :: Qvf3
      integer :: i,j
    end function

    subroutine groupingsetup
    end subroutine

    subroutine setPQ
    end subroutine
  end interface

  print*, 'Setting up variables and constants for run'

  !Controls whether to partition the sdf grid using the 
  !grouping method or not. 
  if( flaggrouping .eq. 1) then
    print*, 'Setting up grouping partitions'
    call groupingsetup
    allocate(f(0:srtgrpx+1,-1:srtgrpm+1),stat=checkst)
    if(checkst.ne.0) then
      write(6,*) 'Problem allocating arrays'
      stop
    end if
    print*,'done!'
  
  else
    xnumg=xnum
    mnumg=mnum
    srtgrpx=xnum
    srtgrpm=mnum
    arrcount=5+1
    ! sets the cutoff of the grid.
    do j=0,mnumg
      do i=1,xnumg     
        if (real(j) .gt. g00 + g11*(real(i)-1._dp)) then 
          cycle
        end if        
        arrcount=arrcount+1
      end do
    end do 

    neq= arrcount-1 !number of equations. 
    allocate(idarrayi(neq),idarrayj(neq),stat=checkst)
    if(checkst.ne.0) then
      write(6,*) 'Problem allocating arrays'
      stop
    end if
    arrcount=5+1
    
    do j=0,mnumg
      do i=1,xnumg     
        if (real(j) .gt. g00 + g11*(real(i)-1._dp)) then 
          cycle
        end if
        if (i .eq. 1 .and. j.eq. 1) minnum=arrcount
        idarrayi(arrcount)=i
        idarrayj(arrcount)=j
        arrcount=arrcount+1
      end do
    end do 

    allocate(f(0:xnumg+1,-1:mnumg+1),stat=checkst)
    if(checkst.ne.0) then
      write(6,*) 'Problem allocating arrays'
      stop
    end if
  end if 

  !Allocation of arrays in 2d. 
  !neq=xnum*mnum+xnum
  allocate(Qiarr(0:xnumg+1,-1:mnumg+1),Pvarr(0:xnumg+1,-1:mnumg+1),Qvarr(0:xnumg+1,-1:mnumg+1), &
  & Qm(0:xnumg+1,-1:mnumg+1),Qg(0:xnumg+1,-1:mnumg+1),bindengHe(1:xnumg,0:mnumg) ,bindengV(1:xnumg,0:mnumg),& 
  & atol(neq), PHearr(0:xnumg+1,-1:mnumg+1),r(1:xnumg),genv(maxgv),stat=checkst)
  if(checkst.ne.0) then
    write(6,*) 'Problem allocating arrays'
    stop
  end if

  allocate(fPvarr(0:xnumg+1,-1:mnumg+1),fQiarr(0:xnumg+1,-1:mnumg+1),fQvarr(0:xnumg+1,-1:mnumg+1),&
& fQm(0:xnumg+1,-1:mnumg+1),fQg(0:xnumg+1,-1:mnumg+1),fPHearr(0:xnumg+1,-1:mnumg+1),stat=checkst)
  if(checkst.ne.0) then
    write(6,*) 'Problem allocating arrays'
    stop
  end if

  allocate(fPvarr0(0:xnumg+1,-1:mnumg+1),fQiarr0(0:xnumg+1,-1:mnumg+1),fPHearr0(0:xnumg+1,-1:mnumg+1),&
& Pvarr0(0:xnumg+1,-1:mnumg+1),Qiarr0(0:xnumg+1,-1:mnumg+1),PHearr0(0:xnumg+1,-1:mnumg+1),stat=checkst)
  if(checkst.ne.0) then
    write(6,*) 'Problem allocating arrays'
    stop
  end if


  !Atomic volume determined by the lattice paramter of bcc iron. 
  !if one wishes to set this manually in the input parameter file, just comment out vol
  vol=((a*1.E-10)**3._dp)/2._dp
  !Binging energy of helium to clusters (only used for certain helium emission options).
  EbHe=2.3
  ksqGB=15._dp/(l*l) !Static grain boundary sink strength (not used anymore).

  !print*,size(Qiarr), size(Pvarr) , size(Qvarr)
  Dv = diff(temp,Emv,D0v) ! Vacancy diffusion coefficient
  Di = diff(temp,Emi,D0i) ! SIA diffusion coefficient
  DHe = diff(temp,EmHe,D0He) ! He diffusion coefficient
  Dg = diff(temp,Emi,D0i) ! Glissile SIA diffusion coefficient
  Qgc=Qgcoeff(Dg,vol) ! Glissile SIA capture coefficient term

  !Setup of capture and emission rates
  !--------------------------------------
  Qiarr=0.0_dp
  Pvarr=0.0_dp
  Qvarr=0.0_dp
  PHearr=0.0_dp
  Qm=0.0_dp
  Qg=0.0_dp
  do j=0,mnumg
    do i=1,xnumg
      if (flaggrouping .eq. 1) then
        mr=meanxarr(i)!rg1xi(i)-0.5_dp*rgx(i)
        if (i .le. srtgrpx) mr=rg1xi(i)
      else 
        mr = real(i,kind=dp) 
      end if
      Qiarr(i,j)= Qif(mr,Di,vol,rcap) !SIA capture rate
      fQiarr(i,j)=Qiarr(i,j)
      if (i .gt. 1) then 
        if (flaggrouping .eq. 1) then
          if( rg1mi(j) .gt. g00 + g11*(rg1xi(i-1)-1._dp)) then 
            Qiarr(i,j)=0.0
            fQiarr(i,j)=0.0
          end if
        else 
          if ( real(j) .gt. g00 + g11*(real(i-1)-1._dp)) then
            Qiarr(i,j)=0.0
            fQiarr(i,j)=0.0
          end if
        end if
      end if
    end do 
  end do


  do j=0,mnumg
    do i=1,xnumg
      if (flaggrouping .eq. 1) then
        mr=meanxarr(i)!rg1xi(i)-0.5_dp*rgx(i)
        if (i .le. srtgrpx) mr=rg1xi(i)
      else 
        mr=real(i,kind=dp)
      end if 
      Pvarr(i,j)= Pf(mr,Dv,vol) !Vacancy capture rate
      fPvarr(i,j)=Pvarr(i,j)
      if (i .eq. xnumg) Pvarr(i,j) =0.0



      if (flaggrouping .eq. 1) then
        if (rg1mi(j) .gt. g00 + g11*(rg1xi(i)-1.)) Pvarr(i,j)=0.0
      else 
        if (real(j) .gt. g00 + g11*(real(i)-1.)) Pvarr(i,j)=0.0
      end if
    end do
  end do 

  do i=1,xnumg
    if (flaggrouping .eq. 1) then
      r(i) = (3._dp*rg1xi(i)*vol/(4._dp*Pi))**(1._dp/3._dp)   !Cluster radius
    else 
      r(i) = (3._dp*real(i)*vol/(4._dp*Pi))**(1._dp/3._dp) 
    end if
  end do 

  do j=0,mnumg
    do i=1,xnumg
      if (flaggrouping .eq. 1) then
        mr=meanxarr(i)!rg1xi(i)-0.5_dp*rgx(i)
        if (i .le. srtgrpx) mr=rg1xi(i)
      else
        mr=real(i,kind=dp)
      end if
      PHearr(i,j)=Pf(mr,DHe,vol) !Helium capture rates
      fPHearr(i,j)=PHearr(i,j)
      if (j .eq. mnumg) PHearr(i,j)=0.0

      if (flaggrouping .eq. 1) then
        if (rg1mi(j+1) .gt. g00 + g11*(rg1xi(i)-1._dp)) PHearr(i,j)=0.0
        if (rg1mi(j) .gt. g00 + g11*(rg1xi(i)-1._dp)) fPHearr(i,j)=0.0
      else 
        if (real(j+1) .gt. g00 + g11*(real(i)-1._dp)) PHearr(i,j)=0.0
      end if
    end do
  end do


  
  
  do j=0,mnumg
    do i=1,xnumg
      if (flaggrouping .eq. 1) then
        mr=meanxarr(i)!rg1xi(i)-0.5_dp*rgx(i)
        mm=meanmarr(j)!rg1mi(j)-0.5_dp*rgm(j)
        if (i .le. srtgrpx) mr=rg1xi(i)
        if (j .le. srtgrpm) mm=rg1mi(j)
      else
        mr=real(i,kind=dp)
        mm=real(j,kind=dp)
      end if
      Qg(i,j)=mr**(2._dp/3._dp)*Qgc !Glissile SIA capture rates
      fQg(i,j)=Qg(i,j)
      !Vacancy emission rates (currently uses Qvf3 and Qvf4
      Qvarr(i,j)= Qvf3(mr,mm,Dv,vol,seng,Ev,temp,i,j)!Qvf4(mr,mm,i,j)!
      fQvarr(i,j)=Qvarr(i,j)
      if (i .eq. 1) Qvarr(i,j)=0.0
      if (i .gt. 1) then
        if (flaggrouping .eq. 1) then 
          if (rg1mi(j) .gt. g00 + g11*(rg1xi(i-1)-1._dp)) then 
            Qvarr(i,j)=0.0
            fQvarr(i,j)=0.0
            Qg(i,j)=0.0
            fQg(i,j)=0.0
          end if 
        else 
          if (real(j) .gt. g00 + g11*(real(i-1)-1._dp)) then
            Qvarr(i,j)=0.0
            fQvarr(i,j)=0.0
            Qg(i,j)=0.0
            fQg(i,j)=0.0
          end if
        end if
      end if 
    end do
  end do


  do j=0,mnumg
    do i=1,xnumg
      if (flaggrouping .eq. 1) then
        mr=meanxarr(i)!rg1xi(i)-0.5_dp*rgx(i)
        mm=meanmarr(j)!rg1mi(j)-0.5_dp*rgm(j)
        if (i .le. srtgrpx) mr=rg1xi(i)
        if (j .le. srtgrpm) mm=rg1mi(j)
      else 
        mr=real(i,kind=dp)
        mm=real(j,kind=dp)
      end if
      !Helium emission rates (currently uses Qemiss and Qemissalt).
      Qm(i,j)=Qemiss(mr,mm,DHe,vol,EbHe,temp,i,j)!Qemissalt(mr,mm,DHe,vol,temp,i,j)!
      fQm(i,j)=Qm(i,j)
      if (flaggrouping .eq. 1) then
        if (i.eq. 1 .and. rg1mi(j) .gt. g00) then 
          Qm(i,j)=0.0
          Qg(i,j)=0.0
        end if
      else 
        if (i.eq. 1 .and. real(j) .gt. g00) then
          Qm(i,j)=0.0
          Qg(i,j)=0.0
        end if
      end if
      if (i.gt.1) then
        if (flaggrouping .eq. 1) then 
          if (rg1mi(j) .gt. g00 + g11*(rg1xi(i)-1._dp)) then 
            Qm(i,j)=0.0
            Qg(i,j)=0.0
          end if
        else
          if (real(j) .gt. g00 + g11*(real(i)-1._dp)) then 
            Qm(i,j)=0.0
            Qg(i,j)=0.0
          
          end if
        end if
      end if
      !write(24,*) Qm(i,j),i,j
    end do 
  end do 
  if (flagEbtable.eq.1) call setPQ !sets tabled option. See MDDFTterms

  !Sink strength terms for SIA glissile clusters. Set to zero for
  ! no SIA glissile contribution. 
  fstterm=ksqft(rhod,disabs) !dislocation   
  secndterm=ksqst(l,Rg) !Grain boundaries

  Qvarr(1,1)=0.0
  fQvarr(1,1)=0.0
  Qiarr(1,0)=0.0
  fQiarr(1,0)=0.0

  do j=0,mnumg
    do i=1,xnumg
      if (i .eq. 1 .and. j .gt. 0) then 
        fQiarr(i,j)=0.0
        fQvarr(i,j)=0.0
        Qiarr(i,j)=0.0
        Qvarr(i,j)=0.0
        Qg(i,j)=0.0
        fQg(i,j)=0.0
      end if 
      if (flaggrouping .eq. 1) then
        if (rg1mi(j) .gt. g00 + g11*(rg1xi(i)-1._dp)) then 
          fPvarr(i,j)=0.0
          fQiarr(i,j)=0.0
          fQvarr(i,j)=0.0
          !fPHearr(i,j)=0.0
          fQm(i,j)=0.0
          fQg(i,j)=0.0
        end if 
      else 
        if (real(j) .gt. g00 + g11*(real(i)-1._dp)) then 
          fPvarr(i,j)=0.0
          fQiarr(i,j)=0.0
          fQvarr(i,j)=0.0
          !fPHearr(i,j)=0.0
          fQm(i,j)=0.0
          fQg(i,j)=0.0
        end if 
      end if
    end do
  end do 
  !-----------------------------------------
  !This subroutine determines the species generation rates
  !specific for DEMO. The results are then fed into speciesgerneration subroutine.

  !call DEMOspeciesgen

  !-----------------------------------------
  fPvarr(1,0)=2._dp*fPvarr(1,0)

  Pvarr(1,0)=2._dp*Pvarr(1,0)

  fPvarr0=fPvarr
  fPHearr0=fPHearr
  fQiarr0=fQiarr
  Pvarr0=Pvarr
  PHearr0=PHearr
  Qiarr0=Qiarr
  !Initial Conditions

  Cv=Cv0*dexp(-Ev/(kB*temp)) !Initial vacancy concentration (equilibrium concentration).
  Cv0=Cv0*dexp(-Ev/(kB*temp)) ! equilibrium vacancy concentration.
  Ci= Ci0!*1.E6 !Initial SIA concentration.
  CHe=0.0 !Initial He concentration
  CGb=0.0 !Initial He concentration at GB.
  Cg=0.0 !Initial small SIA glissile conc. 
  f=0.0_dp !sets size distribution funcion to zero. 
  f(1,0)=Cv
  sv_dis=rhod*Zv !sink strength for vacancies to dislocaitons. 
  si_dis=rhod*Zi !sink strength for SIA to dislocaitons. 
  print*,'done!'
  return
end subroutine initialsetup

