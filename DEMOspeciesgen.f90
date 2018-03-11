subroutine DEMOspeciesgen
  use nuclearvars
  implicit none
  integer(kind=ip) :: dunit,ios,i,j,types,pkacount,tleft,totalpka, &
& ind,n,nv,k,nsia,ndeftypes,nunit,tabledim
  character*30 :: filenamepka

  real(kind=dp),dimension(2,322) :: distfunc
  real(kind=dp), allocatable,dimension(:,:) :: damage
  real(kind=dp), allocatable,dimension(:,:) :: averagedamage
  real(kind=dp) :: disteng,distno,Eth,fmd,fclv,fcli, &
& Np,randomno,Epka,psia,pv,xi,corr


  interface
	  function poisson(mean)
      use kindmod
		  implicit none
		  real(kind=dp) :: mean
		  integer(kind=ip) :: poisson
	  end function poisson

	  function binomial(n,p)
      use kindmod
		  implicit none
		  integer(kind=ip) :: binomial,n
		  real(kind=dp) :: p
	  end function binomial

	  function find_nearest(array, value, arraydim)
      use kindmod
		  implicit none
		  integer(kind=ip) ::find_nearest,arraydim
		  real(kind=dp) :: value
		  real(kind=dp),dimension(2,1:arraydim) :: array
	  end function
  end interface 

  dunit = 15
  filenamepka='CDF_results_Fe.dat'
  open (unit=dunit, file=filenamepka, status='old',    &
             access='sequential', form='formatted', action='read' , iostat=ios)
  if(ios.ne.0) then
    print*, 'problem opening file '
    stop
  endif

  read(dunit, *) 
  read(dunit, *)


  do j=1,322
    

    read(dunit, 110, iostat=ios) disteng, distno
    if(ios.ne.0) print*, 'problem reading file '
    distfunc(1,j)=disteng
		distfunc(2,j)=distno
    110 format( 2(EN13.8,x) )
    !print*, gotEvals(i,j),i,j
    !120 format ( I18, 4 (2X, F12.3) )

    
  end do 
    


  close(dunit)

  
!  dunit=30
!  filename='sampleddist.dat'
!  open (unit=dunit, file=filename,   &
!              iostat=ios)
 ! if(ios.ne.0) then
!    print*, 'problem opening file '
!    stop
!  endif

	!print*,distfunc(2,1:100)




	corr=4.33E14
	types=2
	Eth= 30._dp 
  n=0
	fmd=0.65_dp

	fcli = 0.55_dp
	fclv = 0.25_dp


	pkacount=1
	tleft=0
	totalpka=1E6
  ndeftypes=10
  tabledim=322
	!totalpka=int(totalpka)
  allocate(damage(2,ndeftypes),averagedamage(2,ndeftypes))
	damage=0.0
  averagedamage=0.0
  write(dunit,9)  totalpka
  write(dunit,*) ' '
	do while (pkacount.lt.totalpka)

  !do i=1,322
  !  print*, distfunc(2,i),i,'test'
  !end do 
		tleft=nint(pkacount*10000./real(totalpka))

    n=0


		do while (n .lt. 1)
			call random_number(randomno)
      
			ind= find_nearest(distfunc,randomno,tabledim)
      
			Epka=distfunc(1,ind)
      !print*, ind , Epka, randomno
			Np=Epka*(dexp(-3.56_dp*Epka)+0.3_dp)/(2._dp*Eth)
			n =poisson(Np)
		end do

    write(dunit,10) Epka,randomno

		if (n .gt. 1) then
      
		  n =(n*fmd)
		  
!print*,n,'n'
		  psia=1.-(1.-fcli)**(1./(2.*n-1.))
		  !Parameter for the binomial sampling of SIA-clusters.
		  pv = 1.-(1.-fclv)**(1./(2.*n - 1.))
		  !Parameter for the binomial sampling of V-clusters.
		   
		  !Randomly sample SIA-clusters from the Binomial distribution:
		  nv = 0
		  do while (nv .lt. n)
		    k = binomial(n-1, pv) + 1
		    !assert(k!=0)
		    damage(1,k)=damage(1,k)+1._dp
		    nv = nv + k
			end do

		  nsia = 0
		  do while (nsia .lt. n)
		    k = binomial(n-1,psia)+1
		    !assert(k!=0)
		    damage(2,k)=damage(2,k)+1._dp
		    nsia = nsia + k
			end do

			if (nv .lt. nsia)then 
		    damage(1,1)=damage(1,1)+real(nsia-nv)
		  else if (nsia .lt. nv) then 
		    damage(2,1)=damage(2,1)+real(nv-nsia)
		  end if 

		else 
		  call random_number(xi)
		  if (xi.lt.fmd) then
		    damage(1,1) = damage(1,1)+1._dp
		    damage(2,1) = damage(2,1)+1._dp
			end if
    end if 
!print*, damage(1,:,pkacount),pkacount
		pkacount = pkacount+1
		
		if ((tleft .lt. nint(pkacount*10000./real(totalpka))) .and. (pkacount .lt. totalpka) ) then
      print 22,(pkacount*100./real(totalpka))
		end if
		
	end do 
  pkacount=pkacount-1
	print*, '100.00 % Complete'
  22 format(f7.2,' % Complete ')
  9 format('Sampled distribution from ' , I20,X, ' PKAs')
  10 format(2(EN20.10,X))
  close(dunit)

  do j=1,ndeftypes
    averagedamage(1,j)=damage(1,j)/real(pkacount)
    averagedamage(2,j)=damage(2,j)/real(pkacount)
  end do 

  !print*, pkacount,averagedamage(1,:)
  !print*, pkacount,averagedamage(2,:)
  do i=1,ndeftypes
    averagedamage(1,i)=averagedamage(1,i)*corr*1.E06*vol
    averagedamage(2,i)=averagedamage(2,i)*corr*1.E06*vol
  end do 
  nunit=23
  !open (unit=nunit, access='append',file='arraytest.dat')

  open (unit=nunit, access='append',status='unknown',file='generationarray2.dat')
  write(nunit, *) PKAcount, 'Number of PKAs Sampled'  
  write(nunit, *) 'Vacancy Generation'  ! this gives you the line break
  do i=1,ndeftypes
    write(nunit, '(EN20.10,X)', advance='no') averagedamage(1,i)
  
  print*, averagedamage(1,i)
  end do
  write(nunit, *) ''  ! this gives you the line break
  write(nunit, *) 'SIA Generation'  ! this gives you the line break
  do i=1,ndeftypes
    write(nunit, '(EN20.10,X)', advance='no') averagedamage(2,i)
!print*,damage(2,i)
  end do
  write(nunit, *) ''  ! this gives you the line break
  write(nunit, *) ''  ! this gives you the line break
  close(nunit)

	deallocate(damage,averagedamage)

end subroutine
