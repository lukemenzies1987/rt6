subroutine deallocatearrays
  use nuclearvars
  use twodimensionalvars
  implicit none
deallocate (Qiarr,Pvarr,Qvarr,f,fQg,fQiarr,fPvarr,fPHearr,fQvarr,fQm,&
& Qm,Qg,PHearr,idarrayi,idarrayj,atol,r,genv, & 
& bindengV,bindengHe,Pvarr0,PHearr0,Qiarr0,& 
& fPvarr0,fPHearr0,fQiarr0)


  if (flagEbtable .eq. 1) deallocate(gotEvals)
  if (flaggrouping .eq. 1)  deallocate(g1mi,g1xi,&
& meanxarr,meanmarr,gx,gm,itermx,itermm,&
& L0,L1x,L1m,rg1xi,rg1mi,rgx,rgm)

end subroutine deallocatearrays
