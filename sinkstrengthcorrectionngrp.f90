subroutine sinkstrcorrngrp(y)
!------------------------------------------------------------------
! High density sink strength correction (non-grouping)
!------------------------------------------------------------------ 

  use nuclearvars
  use twodimensionalvars
  implicit none
  real(kind=dp) :: sv_dis0,si_dis0,dddd,ss01,tat0,cr3d,&
& si_void,si_void0,ss01i,ss01v,sv_void,sv_void0,tot_si,&
& tot_si0,tot_sv,tot_sv0,zi_void,sigmasq,n,ktbv,ktbi, &
& besk0,besk1
  real(kind=dp) :: tot_sHe,sHe_void0,ss01He,sHe_void,tot_sHe0
  real(kind=dp), dimension(neq), intent(in) :: y
  real(kind=dp), dimension(1:xnumg,0:mnumg) :: sdf
  integer :: i,j,k

  interface
    function gbss( sss,rgr )
      use kindmod
      implicit none
      real(kind=dp),intent(in) :: sss,rgr
      real(kind=dp) :: gbss
    end function

  end interface

  !-----------------------------------------------------------------------
  !  Dislocations

  sv_dis0 = Zv*rhod  ! SS for vac
  si_dis0 =Zi*rhod   ! SS for SIAs  

  !-----------------------------------------------------------------------
  sv_void=0.0
  si_void=0.0
  sHe_void=0.0
  sv_void0=0.0
  si_void0=0.0
  sHe_void0=0.0
  !do k=5+1,neq
    
  !  i=idarrayi(k)
  !  j=idarrayj(k)
           
  !  sdf(i,j)=abs(y(k))
  !end do
  sdf=0.0
  do k=5+2,neq
    
    i=idarrayi(k)
    j=idarrayj(k)
    sdf(i,j)=abs(y(k))       
    ss01v = fPvarr0(i,j)*sdf(i,j)
    ss01i = fQiarr0(i,j)*sdf(i,j)
    ss01He = fPHearr0(i,j)*sdf(i,j)

    sHe_void0 = sHe_void0 + ss01He

    sv_void0 = sv_void0 + ss01v
    si_void0 = si_void0 + ss01i!*zi_void
  end do
    sv_void0=sv_void0/Dv
    si_void0=si_void0/Di
    sHe_void0=sHe_void0/DHe

    tot_sHe0=sHe_dis+sHe_void0 
    tot_sv0=sv_dis+sv_void0 
    tot_si0=si_dis+si_void0 

  do k=5+1,neq
    
    i=idarrayi(k)
    j=idarrayj(k)
    !Capture rate increase.
    Pvarr(i,j) = Pvarr0(i,j)*(1._dp+r(i)*dsqrt(tot_sv0))
    Qiarr(i,j) = Qiarr0(i,j)*(1._dp+r(i)*dsqrt(tot_si0))
    PHearr(i,j) = PHearr0(i,j)*(1._dp+r(i)*dsqrt(tot_sHe0))     
    fPvarr(i,j) = fPvarr0(i,j)*(1._dp+r(i)*dsqrt(tot_sv0))
    fQiarr(i,j) = fQiarr0(i,j)*(1._dp+r(i)*dsqrt(tot_si0))
    fPHearr(i,j) = fPHearr0(i,j)*(1._dp+r(i)*dsqrt(tot_sHe0))        

    ss01v = fPvarr(i,j)*sdf(i,j)
    ss01i = fQiarr(i,j)*sdf(i,j)
    ss01He = fPHearr(i,j)*sdf(i,j)


    sv_void = sv_void + ss01v
    si_void = si_void + ss01i
    sHe_void = sHe_void + ss01He
    
  end do
  
  sv_void = sv_void-fPvarr(1,0)*sdf(1,0)
  si_void = si_void-fQiarr(1,0)*sdf(1,0)
  sHe_void = sHe_void-fPHearr(1,0)*sdf(1,0)
  sv_void=sv_void/Dv
  si_void=si_void/Di
  sHe_void=sHe_void/DHe

  !Dislocations
  sv_dis = sv_dis0/(1._dp-Zv/(2._dp*pi)*dlog(dsqrt(tot_sv0/sv_dis0)))
  si_dis = si_dis0/(1._dp-Zi/(2._dp*pi)*dlog(dsqrt(tot_si0/si_dis0)))

  tot_sv=sv_void+sv_dis
  tot_si=si_void+si_dis
  tot_SHe=sHe_void+sHe_dis


  ! Using modified bessels function to determine dislocaiton ss(not used but written in).

  !Functions contained in modifiedbessels.f90 
  !n=1.0+0.4
  !ktbv=dsqrt(tot_sv0)*disabs
  !ktbi=dsqrt(tot_si0)*disabs
  !sv_dis=2.*Pi*ktbv*besk1(ktbv)/(ktbv*besk1(ktbv)+besk0(ktbv))*rhod
  !si_dis=2.*Pi*n*ktbi*besk1(ktbi)/(ktbi*besk1(ktbi)+n*besk0(ktbi))*rhod
  ! GB sink strength for 3D

  !--------------------------------------------------------
  sv_gb = 0.0d0
  si_gb = 0.0d0
  sHe_gb=0.0


  sv_gb = Zv*gbss(tot_sv0,l) !GB sink strength
  si_gb = Zi*gbss(tot_si0,l)
  sHe_gb = ZHe*gbss(tot_sHe0,l)


  !-----------------------------------------------------------------------


  return
end subroutine

