subroutine speciesgeneration
!------------------------------------------------------------------
! Species generation rate 
!------------------------------------------------------------------
  use nuclearvars
  implicit none
  real(kind=dp) :: Fedens, &
& s,sak,Gv1,Gi1,ab,cctomc
  integer :: i
  s=0.0
  sak=0.0
  !Fedens=8.5E+28
  cctomc=1.E6

  !All values are detemined using the subroutine DEMOspeciesgen

  !Fw-Region-------------
	!Gv=3.08242802129412E-07
	!Geni=1.96123306347059E-07  
	!Gg=78.1311246547E-09
  !----------------------

  !Blanket-Rear----------
  !Gv1= 10.8367697618E+15     
  !Geni= 6.9521544356E+15     
  !Gg= 2.6916942970E+15  
  !----------------------
  
  !Blanket-Backplate-----
  !Gv=  9.3674752997E+15     
  !Gi=  6.0165379688E+15     
  !Gg=2.3249559427E+15   
  !----------------------

  !Blanket-Middle--------
  !Gv1=18.3952852575E+15
  !Geni=11.7343271720E+15
  !Gg=4.6059384014E+15
  !----------------------

  !These are written in the parameter file and can be used instead of DEMO
  !generation rates (NRT model).
  !Testing model--------
  Gv1=Gdpa*(1.-epsilonr)*(1.-epsilonv)
  Geni=Gdpa*(1.-epsilonr)*(1.-epsiloni)
  Gg=Gdpa*(1.-epsilonr)*epsiloni
  !----------------------

  !Uncomment and comment the above if one wishes to use DEMO generation rates.
  !Gv1=Gv1*vol*cctomc
  !Geni=Geni*vol*cctomc
  !Gg=Gg*vol*cctomc
  Gv=Gv1
  if( flaggrouping .eq. 0) then 
    if (maxgv .gt. xnumg) maxgv=xnumg
  end if 
  if( flaggrouping .eq. 1) then 
    if (maxgv .gt. srtgrpx) maxgv=srtgrpx
  end if 
  do i = 2,maxgv
    genv(i) = (1._dp*i)**(-1._dp)
    s      = s   + genv(i)
    sak    = sak + genv(i)*real(i)
  end do

  ab = Gdpa*(1.-epsilonr)*epsilonv/sak

  do i = 2,maxgv
    genv(i) = ab*genv(i)
  end do
  !Blanket-Middle-----------
  !genv(2)=2.6024292909E+15
  !genv(3)=183.8567852863E+12   
  !genv(4)=8.6014413507E+12
  !genv(5)=299.6619894310E+09
  !genv(6)=8.7206202745E+09
  !genv(7)=177.5300055872E+06
  !-------------------------

  !Blanket-Rear-----------
  !genv(2)=1.5147025028E+15   
  !genv(3)=105.8783205010E+12    
  !genv(4)=4.8867299527E+12 
  !genv(5)=168.2638069782E+09    
  !genv(6)=4.7586701974E+09  
  !genv(7)=86.6000035915E+06 
  !-------------------------

  !Blanket-Backplate-----------
  !genv(2)=1.3041686582E+15    
  !genv(3)=90.3577682673E+12    
  !genv(4)=4.1270141012E+12    
  !genv(5)=139.9932358058E+09     
  !genv(6)=3.9706101647E+09   
  !genv(7)=103.9200043098E+06 
  !-------------------------

  !Fw region---------------
  !genv(2)=44.1986822230E-09 
  !genv(3)=3.1342947307E-09   
  !genv(4)=147.4041853243E-12
  !genv(5)=5.0266058096E-12
  !genv(6)=128.7552717623E-15
  !genv(7)=5.1502108705E-15
  !------------------------

  !Switch between the two for DEMO generation rates or 
  ! NRT rates.

  !genv=genv*vol*cctomc
  genv(1)=Gv1

  return
end subroutine speciesgeneration

