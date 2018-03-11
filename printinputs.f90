!------------------------------------------------------------------
! The file contains the printing subroutines 
!------------------------------------------------------------------    

subroutine printinputs
!------------------------------------------------------------------
! The prints out parameters and setting into the file Initialvaluesandparameters
!------------------------------------------------------------------    
  use nuclearvars
  implicit none
  integer :: nunit
  real(kind=dp) :: daysconv,yearsconv
  nunit=45
  open (unit=nunit,file=foldername//'/'//'Initialvaluesandparameters.dat')

  write(nunit,*)'Initial values and Parameters used in simulation run...'
  if (flagClusters .eq. 1) then 
    write(nunit,*)'for irradiation with the inclusion of helium (2 dimensional SDF)'
  end if 
  if (flagVoids .eq. 1) then 
    write(nunit,*)'for irradiation that does not consider helium (1 dimensional SDF)'
  end if 
  write(nunit,*)
  write(nunit,*)
  write(nunit,10) temp
  write(nunit,20) Dv
  write(nunit,21) Di
  write(nunit,22) DHe
  write(nunit,23) Dg
  write(nunit,*)
  write(nunit,24) a,vol,1._dp/vol
  write(nunit,*)
  write(nunit,25) Cv0
  write(nunit,*)
  if (flaggrouping .eq. 1) then 
    if (flagVoids .eq. 1) then 
      write(nunit,30) xnumg-srtgrpx
      write(nunit,*)
    else
      write(nunit,26) xnumg-srtgrpx,mnumg-srtgrpm
      write(nunit,*)
    end if
  end if 
  daysconv=1._dp/86400._dp
  yearsconv=1._dp/3.154E+07
  write(nunit,27) neq
  write(nunit,*)
  write(nunit,28) rhod
  write(nunit,*)
  write(nunit,29) tmax,tmax*daysconv,tmax*yearsconv
  write(nunit,*)
  if (flagEbtable .eq. 1) then 
    write(nunit,*)'Note: tabled values were used for binding energies',&
& ' for helium and/or vacancies to clusters'
  end if 
  if (flaggrouping .eq. 1) then 
    write(nunit,*)'Note: the grouping method was invoked in this simulation run'
  end if 
  10 format('Diffusion coefficitents for a temperature of ',f6.2,' K:')
  20 format('Dv (vacancies)',1pe16.7)
  21 format('Di (SIA)',1pe16.7)
  22 format('DHe (Helium)',1pe16.7)
  23 format('Dcrow (crowdions)',1pe16.7)
  24 format('Atomic volume for a lattice parameter of ',f5.2,' Angstroms:',& 
& EN15.5, ' and metal density: ' ,1pe16.7)
  25 format('Vacancy concentration prior to irradiation',1pe16.7,' (per atom)')
  26 format('Number of groups not including group widths of 1 for vacancies: ',&
& I4, ' and helium atoms: ',I4)
  27 format('Number of equations the solver has to solve: ' , I5)
  28 format('Dislocation density (m^-3)',1pe16.7)
  29 format('Simulated irradiation time: ',1pe16.7,X,' s',0pf16.4,X,' days ',0pf16.4,X,' years')
  30 format('Number of groups not including group widths of 1 for vacancies: ', I4)
  close(nunit)



end subroutine printinputs

subroutine printabsemissrates2d
!------------------------------------------------------------------
! The prints out the 2 dimensional binding energies, capture and emission rates
! for clusters
!------------------------------------------------------------------  
  use nuclearvars
  use twodimensionalvars
  implicit none
  integer :: vabs,vemiss,heabs,heemiss,intabs,bindHe, &
& bindV,crowd,i,j

  vabs=45
  vemiss=46
  heabs=47
  heemiss=48
  intabs=49
  bindHe=50
  bindV=51
  crowd=52

  open (unit=vabs,file=foldername//'/'//'Rates-VacancyAbsorptionRates.dat')
  open (unit=vemiss,file=foldername//'/'//'Rates-VacancyEmissionRates.dat')
  open (unit=heabs,file=foldername//'/'//'Rates-HeliumAbsorptionRates.dat')
  open (unit=heemiss,file=foldername//'/'//'Rates-HeliumEmissionRates.dat')
  open (unit=intabs,file=foldername//'/'//'Rates-InterstitialAbsorptionRates.dat')
  open (unit=bindHe,file=foldername//'/'//'BindingEnergyHe-Cluster.dat')
  open (unit=bindV,file=foldername//'/'//'BindingEnergyV-Cluster.dat')
  open (unit=crowd,file=foldername//'/'//'Rates-VacancyAbsorptionRatesCrowdions.dat')

  do j=0,mnumg
    do i=1,xnumg
      write(vabs, '(EN11.2,X)', advance='no') Pvarr(i,j)
      write(vemiss, '(EN11.2,X)', advance='no') Qvarr(i,j)
      write(intabs, '(EN11.2,X)', advance='no') Qiarr(i,j)
      write(heabs, '(EN11.2,X)', advance='no') PHearr(i,j)
      write(heemiss, '(EN11.2,X)', advance='no') Qm(i,j)
      write(bindHe, '(f6.2,X)', advance='no') bindengHe(i,j)
      write(bindV, '(f6.2,X)', advance='no') bindengV(i,j)
      write(crowd, '(EN11.2,X)', advance='no') Qg(i,j)
    end do
    write(vabs, *) ''  ! this gives you the line break
    write(vemiss, *) ''  ! this gives you the line break
    write(heabs, *) ''  ! this gives you the line break
    write(heemiss, *) ''  ! this gives you the line break
    write(intabs, *) ''  ! this gives you the line break
    write(bindHe, *) ''  ! this gives you the line break
    write(bindV, *) ''  ! this gives you the line break
    write(crowd, *) ''  ! this gives you the line break
  end do

  close(vabs)
  close(vemiss)
  close(heabs)
  close(heemiss)
  close(intabs)
  close(bindHe)
  close(bindV)
  close(crowd)

end subroutine printabsemissrates2d

subroutine printabsemissrates1d
!------------------------------------------------------------------
! The prints out the 1 dimensional binding energies, capture and emission rates
! for clusters
!------------------------------------------------------------------  
  use nuclearvars
  use onedimensionalvars
  implicit none
  integer :: vabs,vemiss,intabs,i

  vabs=45
  vemiss=46
  intabs=49

  open (unit=vabs,file=foldername//'/'//'Rates-VacancyAbsorptionRates.dat')
  open (unit=vemiss,file=foldername//'/'//'Rates-VacancyEmissionRates.dat')
  open (unit=intabs,file=foldername//'/'//'Rates-InterstitialAbsorptionRates.dat')
    do i=1,xnumg
      write(vabs, '(EN11.2,X)') Pvarr(i)
      write(vemiss, '(EN11.2,X)') Qvarr(i)
      write(intabs, '(EN11.2,X)') Qiarr(i)
  end do

  close(vabs)
  close(vemiss)
  close(intabs)

end subroutine printabsemissrates1d


subroutine printintvalues(t,clss_V,clss_i,clss_He)
!------------------------------------------------------------------
! The prints evolving sink strength values for 2d 
!------------------------------------------------------------------  
  use nuclearvars
  implicit none
  real(kind=dp) :: t,clss_V,clss_i,clss_He


  write(printint,11) t,sv_gb,si_gb,sHe_gb,sv_dis,si_dis,disl/(DHe*CHe),clss_V,clss_i,clss_He

  11 format(1pe10.2,9(EN15.5,X))
end subroutine printintvalues

subroutine printintvalues1d(t)
!------------------------------------------------------------------
! The prints evolving sink strength values for 1d 
!------------------------------------------------------------------  
  use nuclearvars
  real(kind=dp) :: t

  write(printint,11) t,sv_gb,si_gb,sv_dis,si_dis

  11 format(1pe10.2,4(EN15.5,X))
end subroutine printintvalues1d
