subroutine setPQ
!------------------------------------------------------------------
! Tabled option for vacancy and He emission terms.
! The data is stored in the file bindingengdat.dat,
! keep this formatting.
!------------------------------------------------------------------ 
  use nuclearvars
  use twodimensionalvars
  implicit none
  
  integer :: dunit,ios,i,j,k,length,tag
  character*30 :: nov,noa

  real(kind=dp) ::Evaltemp1,Evaltemp2,Evaltemp3,Evaltemp4,& 
& Evaltemp5,Evaltemp6, xr,w,Eb,EbHE,expont,expontHe,qvf
  dunit = 15

  open (unit=dunit, file=filename, status='old',    &
             access='sequential', form='formatted', action='read' , iostat=ios)
  if(ios.ne.0) then
    print*, 'problem opening file: ',filename
    stop
  endif 
  read(dunit, '(i5)') length
  read(dunit, *) 

  print*, length!,noofa

  allocate(gotEvals(6,1:length))

  
  do i=1,length

    read(dunit, *, iostat=ios) Evaltemp1,Evaltemp2,&
& Evaltemp3,Evaltemp4,Evaltemp5,Evaltemp6
    if(ios.ne.0) print*, 'problem reading file: ',filename
    gotEvals(1,i)=Evaltemp1
    gotEvals(2,i)=Evaltemp2
    gotEvals(3,i)=Evaltemp3
    gotEvals(4,i)=Evaltemp4
    gotEvals(5,i)=Evaltemp5
    gotEvals(6,i)=Evaltemp6
   
    110 format( 6(f6.4,2X) )

  end do


  close(dunit)
  
  !do j=0,noofa
  !  do i=2,noofv 
  !    write(9,*) gotEvals(i,j),gotEvals(i-1,j)+Ev-gotEvals(i,j),'MD',i,j
      !print*, Qvf(i,j,Dv,vol,seng,Ev,temp),'Ana',i,j
  !  end do
  !end do  
  tag=0
  if (flagGrouping .eq. 1) then 
    do j=0,mnumg
      do i=1,xnumg 
        tag=0
        do k=1,length
          if(gotEvals(2,k) .eq. rg1xi(i) .and. gotEvals(1,k) .eq. rg1mi(j)) then 
            Eb=gotEvals(6,k)
            EbHe=gotEvals(5,k)
            tag=1
            exit
          else
            cycle
          end if 
        end do

        if (tag .eq. 0) then 
          tag=0
          cycle
        end if
        w= ((48.0_dp*Pi*Pi)/(vol*vol))**(1.0_dp/3.0_dp)
        xr=rg1xi(i)
        !if (Eb .lt. 0.6) Eb=0.6
        expont=dexp(-Eb/(kb*temp))
        expontHe=dexp(-EbHe/(kb*temp))
        if (Qvarr(i,j) .ne. 0.0) Qvarr(i,j)=w*Dv*xr**(1.0_dp/3.0_dp)*Cv0*expont
        if (Qm(i,j) .ne. 0.0) Qm(i,j)=w*DHe*xr**(1.0_dp/3.0_dp)*expontHe
        !Qvarr(i,j)=w*Dv*xr**(1.0_dp/3.0_dp)*expont

      end do
    end do  
  else
    do j=0,mnumg
      do i=1,xnumg 
        tag=0
        do k=1,length
          if(nint(gotEvals(2,k)) .eq. i .and. nint(gotEvals(1,k)) .eq. j) then 
            Eb=gotEvals(6,k)
            EbHe=gotEvals(5,k)
            bindengHe(i,j)=EbHe
            bindengV(i,j)=Eb
            tag=1
            exit
          else
            cycle
          end if 
        end do
        
        !if (g1xi(i) .gt. noofv .or. g1mi(j) .gt. noofa) exit
        if (tag .eq. 0) then 
          tag=0
          cycle
        end if
        w= ((48.0_dp*Pi*Pi)/(vol*vol))**(1.0_dp/3.0_dp)
        xr=real(i)
        !if (Eb .lt. 0.6) Eb=0.6
        expont=dexp(-Eb/(kb*temp))
        expontHe=dexp(-EbHe/(kb*temp))
        if (Qvarr(i,j) .ne. 0.0) Qvarr(i,j)=w*Dv*xr**(1.0_dp/3.0_dp)*Cv0*expont
        if (Qm(i,j) .ne. 0.0) Qm(i,j)=w*DHe*xr**(1.0_dp/3.0_dp)*expontHe
        !Qvarr(i,j)=w*Dv*xr**(1.0_dp/3.0_dp)*expont

      end do
    end do
  end if
  return
end subroutine
