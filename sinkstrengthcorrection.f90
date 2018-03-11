subroutine sinkstrcorr(y)
!------------------------------------------------------------------
! High density sink strength correction (grouping)
!------------------------------------------------------------------ 

  use nuclearvars
  use twodimensionalvars
  implicit none
  real(kind=dp) :: sv_dis0,si_dis0,dddd,ss01,tat0,cr3d,&
& si_void,si_void0,ss01i,ss01v,sv_void,sv_void0,tot_si,&
& tot_si0,tot_sv,tot_sv0,zi_void
  real(kind=dp), dimension(neq) :: y
  integer :: i,j,k

  interface

    function gbss( sss,rgr )
      use kindmod
      implicit none
      real(kind=dp),intent(in) :: sss,rgr
      real(kind=dp) :: gbss
    end function

    function sigmasq(dx)
      use kindmod
      implicit none
      real(kind=dp) :: dx
      real(kind=dp) :: sigmasq
    end function

  end interface
  !-----------------------------------------------------------------------
  ! Dislocations

  sv_dis0 = Zv*rhod ! SS for vac
  si_dis0 =Zi*rhod  ! SS for SIA  

  !-----------------------------------------------------------------------


  sv_void0 = 0.0
  si_void0 = 0.0


  do k=5+1,neq1

    i=idarrayi(k)
    j=idarrayj(k)
    L0(i,j)=abs(y(k))
  end do 

  do k=neq1+1,neq,3
    i=idarrayi(k)
    j=idarrayj(k)
    L0(i,j)=abs(y(k))
    L1x(i,j)=y(k+1)
  end do 



  do j = 1,mnumg
    do i = 2,xnumg

      !if (i.eq.1.and.j.eq.0) cycle

      !r       = rg1xi(i)
      zi_void = 1. + rcap/r(i)

      dddd = gx(i)*gm(j)
      ss01 = 4.*pi*r(i)*L0(i,j)/vol*dddd
      tat0 = (L1x(i,j)*sigmasq(rgx(i))+L0(i,j)*meanxarr(i))*dddd


      sv_void0 = sv_void0 + ss01
      si_void0 = si_void0 + ss01*zi_void

    enddo
  enddo



  tot_sv0 = sv_dis0 + sv_void0     ! total SS vac 
  tot_si0 = si_dis0 + si_void0      ! total SS SIAs

  !-----------------------------------------------------------------------
  !                              High-density corrections to sink strength
  !-----------------------------------------------------------------------

  sv_void = 0.0d0 
  si_void = 0.0d0 


                 
  ! Bubbles
  do j=0,mnumg
    do i = 1,xnumg

      !r       = rg1xi(i) 
      cr3d    = 4.*pi*r(i)/vol
      zi_void = 1. + rcap/r(i)



      Pvarr(i,j) = cr3d*Dv *(1.+r(i)*dsqrt(tot_sv0))
      Qiarr(i,j) = cr3d*Di *(1.+r(i)*dsqrt(tot_si0))
      PHearr(i,j) = cr3d*DHe*(1.+r(i)*dsqrt(tot_sv0))     ! ???

      dddd = gm(j)*gx(i)
      ss01v = Pvarr(i,j)*L0(i,j)*dddd/Dv
      ss01i = Qiarr(i,j)*L0(i,j)*dddd/Di

      sv_void = sv_void + ss01v
      si_void = si_void + ss01i*zi_void

    end do
  end do






  !-----------------------------------------------------------------------
  !                                                          Dislocations

  sv_dis = sv_dis0/(1.-Zv/(4.*pi)*dlog(tot_sv0/sv_dis0))
  si_dis = si_dis0/(1.-Zi/(4.*pi)*dlog(tot_si0/si_dis0))

  !-----------------------------------------------------------------------

  tot_sv = sv_dis + sv_void ! total SS vac 
  tot_si = si_dis + si_void  ! total SS SIAs 

  !-----------------------------------------------------------------------
  !                                              ! GB sink strength for 3D
  sv_gb = 0.0d0
  si_gb = 0.0d0

  !if (key.eq.0) then
  sv_gb = Zv*gbss(tot_sv,l)
  si_gb = Zi*gbss(tot_si,l)
  !endif

  !-----------------------------------------------------------------------


  return
end subroutine

subroutine sinkstrcorr1d(y)
!------------------------------------------------------------------
! High density sink strength correction (grouping) 1 dimensional
!------------------------------------------------------------------ 

  use nuclearvars
  use onedimensionalvars
  implicit none
  real(kind=dp) :: sv_dis0,si_dis0,dddd,ss01,tat0,cr3d,&
& si_void,si_void0,ss01i,ss01v,sv_void,sv_void0,tot_si,&
& tot_si0,tot_sv,tot_sv0,zi_void,si_il,sv_il, &
& si_il0,sv_il0,sv_gb0,si_gb0
 

  real(kind=dp):: besk0,besk1,ktb,n
  real(kind=dp), dimension(neq) :: y
  integer :: i,j,k

  interface

    function gbss( sss,rgr )
      use kindmod
      implicit none
      real(kind=dp),intent(in) :: sss,rgr
      real(kind=dp) :: gbss
    end function

    function sigmasq(dx)
      use kindmod
      implicit none
      real(kind=dp) :: dx
      real(kind=dp) :: sigmasq
    end function

  end interface

  !-----------------------------------------------------------------------
  !                                                          Dislocations

  sv_dis0 = Zv*rhod  ! SS for vac
  si_dis0 =Zi*rhod   ! SS for SIAs  

  !-----------------------------------------------------------------------

  !sa_void  = 0.0d0
  sv_void0 = 0.0
  si_void0 = 0.0

  sv_gb0=ksqgb
  si_gb0=ksqgb

  do k=2+1,neq1
    i=idarrayi(k)
    L0x(i)=abs(y(k))
  end do 

  do k=neq1+1,neq2,2
    i=idarrayi(k)
    L0x(i)=abs(y(k))
    L1x(i)=y(k+1)
  end do 


 ! do k=neq2+1,neq3
 !   i=idarrayi(k)
 !   fi(i)=abs(y(k))
!  end do 

 ! do k=neq3+1,neq,2
 !   i=idarrayi(k)
 !   L0i(i)=abs(y(k))
 !   L1i(i)=y(k+1)
 ! end do 


  do i = 2,xnumg

    !if (i.eq.1.and.j.eq.0) cycle

    !r       = rg1xi(i)
    zi_void = 1. + rcap/r(i)

    dddd = rgx(i)
    ss01 = 4.*pi*r(i)*L0x(i)/vol*dddd
    tat0 = (L1x(i)*sigmasq(rgx(i))+L0x(i)*meanxarr(i))*dddd

    !sa_void  = sa_void  + tat0 
    sv_void0 = sv_void0 + ss01
    si_void0 = si_void0 + ss01*zi_void


  end do

  sv_il0 = 0.0d0 
  si_il0 = 0.0d0
                           

  tot_sv0 = sv_dis0 + sv_void0   ! total SS vac 
  tot_si0 = si_dis0 + si_void0    ! total SS SIAs

  !-----------------------------------------------------------------------
  ! High-density corrections to sink strength
  !-----------------------------------------------------------------------

  sv_void = 0.0d0 
  si_void = 0.0d0 


  !-----------------------------------------------------------------------
                 
   ! Voids
  
  do i = 2,xnumg

    !r       = rg1xi(i) 
    cr3d    = 4.*pi*r(i)/vol
    zi_void = 1. + rcap/r(i)



    Pvarr(i) = cr3d*Dv *(1.+r(i)*dsqrt(tot_sv0))
    Qiarr(i) = cr3d*Di *(1.+r(i)*dsqrt(tot_si0))


    dddd = rgx(i)
    ss01v = Pvarr(i)*L0x(i)*dddd/Dv
    ss01i = Qiarr(i)*L0x(i)*dddd/Di

    sv_void = sv_void + ss01v
    si_void = si_void + ss01i*zi_void

  end do
  

  !-----------------------------------------------------------------------
  !  Dislocations

  sv_dis = sv_dis0/(1.-Zv/(4.*pi)*dlog(tot_sv0/sv_dis0))
  si_dis = si_dis0/(1.-Zi/(4.*pi)*dlog(tot_si0/si_dis0))


  !-----------------------------------------------------------------------

  tot_sv = sv_dis + sv_void ! total SS vac 
  tot_si = si_dis + si_void ! total SS SIAs 

  !-----------------------------------------------------------------------
  !                                              ! GB sink strength for 3D
  sv_gb = 0.0d0
  si_gb = 0.0d0

  ! Using modified bessels function to determine dislocaiton ss (not used but written in).
  !Functions contained in modifiedbessels.f90 

  !n=1.0+0.4
  !ktb=dsqrt(tot_sv)*disabs
  !sv_dis=2.*Pi*ktb*besk1(ktb)/(ktb*besk1(ktb)+besk0(ktb))*rhod
  !si_dis=2.*Pi*n*ktb*besk1(ktb)/(ktb*besk1(ktb)+n*besk0(ktb))*rhod

  sv_gb = gbss(tot_sv,l) !GB sink strength
  si_gb = gbss(tot_si,l)


  !-----------------------------------------------------------------------


  return
end subroutine



