subroutine deallocatearrays1d
  use nuclearvars
  use onedimensionalvars
  implicit none
  deallocate(Qiarr,Pvarr,Qvarr,fx,fi, &
& Qg,g1xi,bindingV, &
& meanxarr,gx,itermx,itermi,&
& L0x,L0i,L1x,L1i,rg1xi, &
& idarrayi,atol,r,ril,genv,rgx, &
& Qi,gi,rgi,meaniarr,Geniloop,rg1ii,g1ii,itermg1)
  if (flagEbtable .eq. 1) deallocate(gotEvals)


end subroutine deallocatearrays1d
