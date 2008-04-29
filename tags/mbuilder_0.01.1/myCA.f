c	INITIAL MICROSTRUCTURE
c       
c       This version reads in a list of ellipsoid centres, and their shapes
c       and uses these as nucleation points in a CA to generate a microstructure
!       The ellipsoid nuclei (centers) are created using recursivesampler
!       and bell_curve.pl which are refined by ellipticalFoam
!       The active.list and elliptical.cells files contain ID's and aspect
!       ratios of the ellipsoids
!
!  To compile: ${F90} -c myCA.f -o myCA.o
!              ${CC} -c mt19937.c -o mt19937.o
!              ${F90} myCA.o mt19937.o -o myCA 
!
!  written by ADR 2002-2004
!  modified by CGR 11/05 
!  modified by CGR 3/06  asks whether making columnar or equiaxed grains
!  modified by CGR 1/07 --> reads elliptical.cells file format
!  modified by CGR 4/07 --> currently reads active.list and elliptical.cells
!                           for grain ID and aspect ratios
!  CGR 10/13/07 --> fixed the scalability problem on non-cubic grids.

!___________________________________________________________________
	module var

	implicit none
        save

! Defines boundary conditions
	logical::periodic

! Multiplication factor which increases the search radius after each step	
	double precision::delta_t

!RVE dimensions & scaling factor
	real::rve_x,rve_y,rve_z
	double precision::mag

	type lattice
	real,dimension(3):: center
	real,dimension(3):: size
	real,dimension(3):: shape
	end type lattice

	integer,dimension(:,:,:),allocatable::spins
	integer::xmax,ymax,zmax,m3
	integer,parameter::nbors=26
	integer,dimension(nbors,3)::nnbors

	type(lattice),dimension(:),allocatable::mesh
	
	integer,dimension(1:1000)::ticktock
	real,dimension(3)::axes
	integer::ngrains,num_grown
	integer,dimension(:),allocatable:: ell_volume
	
	end module var
!___________________________________________________________________

	program structure_builder

	use var
	implicit none
	character:: name*30
	integer:: inbr,jnbr,knbr
	integer:: i,j,k,l,ii,jj,kk,m
	integer::irandom,switch,timebomb
	integer:: num_ells
	integer::mcs_count  !MCS counters
	integer:: nnspin,neigh
	logical:: all_0_nn,inside
	integer:: try_it,checked(18)  !  flags for neighbors
	integer:: index,num
	real:: distance,distance_nn(18),x,y,z
	integer:: state=1
	real,parameter:: pi=3.14159265

	i=1
	mcs_count=1
	switch=1
! Obtain DATA from file or user input
	call input(num_ells)

	call prep() !set up nbor tables

	goto 923
	open(1,file='test.output')

	do i=1,ngrains
	   write(1,*)mesh(i)%center
	   write(1,*)mesh(i)%shape
	   write(1,*)mesh(i)%size
	enddo
	close(1)
 923	continue

	do while(num_grown.lt.m3) !1MCS step after each loop

	   ticktock(mcs_count)=num_grown

! increases size of each cell in the x,y, and z directions
	   do j=1,ngrains	
	      mesh(j)%size(1)=mesh(j)%shape(1)*real(mcs_count)*delta_t
	      mesh(j)%size(2)=mesh(j)%shape(2)*real(mcs_count)*delta_t
	      mesh(j)%size(3)=mesh(j)%shape(3)*real(mcs_count)*delta_t
	   enddo	

	   if(switch.eq.1)switch=timebomb(mcs_count) 
	   write(*,*)num_grown,'MCS time=',mcs_count,' switch ',switch

	   do i=1,m3  !m3 random site selections = 1 MCS

	      select case(switch)
	      case(1)
		 ii=irandom(xmax)
		 jj=irandom(ymax)
		 kk=irandom(zmax)
		 !num=irandom(m3)
		 !kk=mod(num,zmax)
		 !jj=mod(1+((num-1)/zmax),ymax)
		 !ii=mod(1+((num-1)/(zmax*ymax)),xmax)
	      case(2)
		 do j=1,xmax
		    do k=1,ymax
		       do l=1,zmax
			  if(spins(j,k,l) == 0)then
			     ii=j
			     jj=k
			     kk=l
			  endif
		       enddo
		    enddo
		 enddo
	      end select

 456	      if(spins(ii,jj,kk) == 0)then
		 all_0_nn = .true.
		 index=0
		 
		 do m = 1,18
		    checked(m) = -1
		    distance_nn(m)=99999.
		 enddo
! Search 1st and 2nd NN for a non-cero spin

		 do try_it = 1,18 !  go around in order
		    call neighs3d(ii,jj,kk,try_it,inbr,jnbr,knbr)
		    nnspin=spins(inbr,jnbr,knbr)
		    if(nnspin /= 0) then
		       all_0_nn=.false.
		       x=real(ii)
		       y=real(jj)
		       z=real(kk) 
		       call test_inside(x,y,z,nnspin,inside,distance)
		       if(inside) then
			  checked(try_it) = nnspin
			  distance_nn(try_it) = distance
		       endif	!  if(inside)
		    endif	!  if(nnspin

		 enddo		!  do try_it=1,18

		 if(.not.all_0_nn)then
		       call growth(distance_nn,index)
		 endif
		    
		 if(index.gt.0) then
		    neigh=checked(index)
		    spins(ii,jj,kk)=neigh
		    ell_volume(neigh) = ell_volume(neigh) + 1
		    num_grown = num_grown+1
		 endif
		 
	      endif !  if(spin == 0)
	   enddo !  do i=1,m3
	   mcs_count=mcs_count+1
	enddo			!  do i=1,99999
	
 	write(*,*)'Done with growth'
	name='ellipsoid'
	call write_voxels(name,ngrains,mcs_count)

	end
!________________________________________________________
	subroutine write_voxels(name,num_ells,mcs_count)

	use var
	implicit none
	character:: name*30,name2*30,keyword*6
	integer::j,mcs_count
	integer:: isite,jsite,ksite,ilength,num_ells
	integer::unit,tens
	real::rbar=1.0,temp=1.0,t
	character::a,b

	keyword='BaseCA'
	write(*,*)'time=',mcs_count

	unit=int(mcs_count/10.0**tens)
	tens=int(log10(mcs_count*1.0))
	a=char(48+unit)
	b=char(48+tens)
        t=mcs_count*1.0
        write(*,*)'t=',t

	name2 = name(1:len_trim(name))//'_'//keyword//'.ph'
	write(*,*)'  name2 = ',name2, ' is the output filename'

	open(17,file=name2,status='unknown')
	write(17,989) xmax,ymax,zmax
	write(17,990) keyword,t,rbar,temp,num_ells
	write(17,*) 3,0.0,0.0,0,0
	do ksite=1,zmax
	   do jsite=1,ymax
	      write(17,992) (spins(isite,jsite,ksite),isite=1,xmax)
 29	   enddo
 30	enddo
 989	format(3i8)
 990	format("'",a6,"'",f19.2,f7.3,f5.1,i10)
 991	format(f6.3,f6.3,f6.3,(21x,i10))
 992	format(20i6)
	close(17)

	return
	end
!       __________________________________________________________

	subroutine test_inside(ix,iy,iz,nnspin,inside,distance)
!
! Including periodic B.C. if periodic=.true. in module var
! If distance between center and point is > 0.5
! then 'temporarily' move the point, calculate distance, and 
! if <1.0, update checked() and distance() arrays.
!	
	use var
	implicit none
	integer::nnspin
	logical:: inside
	real::temp1,temp2,temp3,distance
	real,save::size_x,size_y,size_z
	real::ix,iy,iz
	real::xcenter,ycenter,zcenter

! Symmetry line along each box axis is 1/2 of the length
	size_x=0.5*xmax
	size_y=0.5*ymax
	size_z=0.5*zmax
	
! Convert RVE ellipsoid centers to the
! voxel coordinates by multiplying by the scale factor
	xcenter=mesh(nnspin)%center(1)*mag
	ycenter=mesh(nnspin)%center(2)*mag
	zcenter=mesh(nnspin)%center(3)*mag
	

	if(periodic)then
	
	   if(abs(ix - xcenter).gt.size_x) then
	      if((ix - xcenter).gt.size_x) then
		 ix = ix - 1.0*xmax
	      endif
	      if((xcenter - ix).gt.size_x) then
		 ix = ix + 1.0*xmax
	      endif
	   endif
	   
	   if(abs(iy - ycenter).gt.size_y) then
	      if((iy -ycenter).gt.size_y) then
		 iy = iy - 1.0*ymax
	      endif
	      if((ycenter - iy).gt.size_y) then
		 iy = iy + 1.0*ymax
	      endif
	   endif
	   
	   if(abs(iz - zcenter).gt.size_z) then
	      if((iz -zcenter).gt.size_z) then
		 iz = iz - 1.0*zmax
	      endif
	      if((zcenter - iz).gt.size_z) then
		 iz = iz + 1.0*zmax
	      endif
	   endif
	
	endif

!	temp1=((ix-mesh(nnspin)%center(1))/mesh(nnspin)%size(1))**2
!	temp2=((iy-mesh(nnspin)%center(2))/mesh(nnspin)%size(2))**2
!	temp3=((iz-mesh(nnspin)%center(3))/mesh(nnspin)%size(3))**2

! New way
! ix iy and iz are the voxel coordinates
! xcenter ycenter and zcenter are the voxel coordinates of the
! center of the ellipsoid
! ix - xcenter = diff in x
! iy - ycenter = diff in y
! iz - zcenter = diff in z
! 
! The size of the ellipsoid must be converted to the voxel length
! scale by multiplying by the scale factor, mag.
 
! I could have just easily divided the voxel coordinates
! ix, iy, and iz by Mag and left everything else the same.

	temp1=((ix-xcenter)/(mesh(nnspin)%size(1)*mag))**2
	temp2=((iy-ycenter)/(mesh(nnspin)%size(2)*mag))**2
	temp3=((iz-zcenter)/(mesh(nnspin)%size(3)*mag))**2

	distance=sqrt(temp1+temp2+temp3)

	inside =(distance.le.1.0)
	return
	end

   !__________________________________________________________

	subroutine neighs3d(isite,jsite,ksite,nbor,inbr,jnbr,knbr)

      ! (isite,jsite,ksite)=current
      ! (inbr,jnbr,knbr)=neighbor
      ! periodic = periodic BC in ALL directiosn ( x and  y and z)

	use var
	
	implicit none
	integer:: isite,jsite,ksite,nbor
	integer,intent(out)::inbr,jnbr,knbr
	!logical::periodic
	integer,dimension(3)::low,high
	
	high(1)=xmax
	high(2)=ymax
	high(3)=zmax
	low=1
	
	
	inbr=isite+nnbors(nbor,1)
	jnbr=jsite+nnbors(nbor,2)
	knbr=ksite+nnbors(nbor,3)
	
	
	if(periodic)then
	   inbr=mod((inbr+high(1)-1),high(1))+1
	   jnbr=mod((jnbr+high(2)-1),high(2))+1
	   knbr=mod((knbr+high(3)-1),high(3))+1
	else
	   if(inbr>high(1))inbr=high(1)
	   if(jnbr>high(2))jnbr=high(2)
	   if(knbr>high(3))knbr=high(3)
	   if(inbr<low(1))inbr=low(1)
	   if(jnbr<low(2))jnbr=low(2)
	   if(knbr<low(3))knbr=low(3)
	endif
	
	return
	end subroutine neighs3d
	
!__________________________

!	 PREP :  Prepares the array and parameters for the run
!	-----------------------------------------------------------------
!     initializes nearest neighbor array
! taken from rex3d code

	subroutine prep()
	use var
	implicit none
	integer:: i,j

        nnbors=0

	nnbors(1,1)=1
	nnbors(2,2)=1
	nnbors(3,3)=1
	nnbors(4,1)=-1
	nnbors(5,2)=-1
	nnbors(6,3)=-1
	nnbors(7,1)=1
	nnbors(7,2)=1
	nnbors(8,1)=1
	nnbors(8,3)=1
	nnbors(9,2)=1
	nnbors(9,3)=1
	nnbors(10,1)=-1
	nnbors(10,2)=1
	nnbors(11,1)=-1
	nnbors(11,3)=1
	nnbors(12,1)=-1
	nnbors(12,2)=-1
	nnbors(13,1)=-1
	nnbors(13,3)=-1
	nnbors(14,2)=-1
	nnbors(14,1)=1
	nnbors(15,2)=-1
	nnbors(15,3)=1
	nnbors(16,2)=-1
	nnbors(16,3)=-1
	nnbors(17,3)=-1
	nnbors(17,1)=1
	nnbors(18,3)=-1
	nnbors(18,2)=1
	nnbors(19,1)=1
	nnbors(19,2)=1
	nnbors(19,3)=1
	nnbors(20,1)=-1
	nnbors(20,2)=1
	nnbors(20,3)=1
	nnbors(21,1)=1
	nnbors(21,2)=-1
	nnbors(21,3)=1
	nnbors(22,1)=1
	nnbors(22,2)=1
	nnbors(22,3)=-1
	nnbors(23,1)=-1
	nnbors(23,2)=-1
	nnbors(23,3)=1
	nnbors(24,1)=-1
	nnbors(24,2)=-1
	nnbors(24,3)=-1
	nnbors(25,1)=1
	nnbors(25,2)=-1
	nnbors(25,3)=-1
	nnbors(26,1)=-1
	nnbors(26,2)=1
	nnbors(26,3)=-1

	return
	end
!__________________________
	integer function timebomb(time)

! if system cannot find an empty site(i,j,k) with spin=0
! at late times when using a RNG, it will permit the switch to
! an ordered search (i.e. do i=1,mx   do j=1,my   do k=1,mz)
!
! 1/15/07 -- modified it so that delta_t= 2 x delta_t
!            once and if a problem occurs again, switch to
!            ordered searches
! 
!
	use var
	implicit none
	logical:: loop_problem
	integer::time,i
	integer,save::onetry=1

	if(time.gt.5.and.onetry==1)then
	   i=time-5
	   loop_problem=.true.
	   do while(i<time)
	      if(ticktock(time-5).ne.ticktock(i))then
		 loop_problem=.false.
		 exit
	      endif
	      i=i+1
	   enddo
	   if(loop_problem)then
	      timebomb=1
	      delta_t=2.0D0*delta_t
	      write(*,*)'Increasing delta_t by factor of 2'
	      onetry=onetry+1
	      return
	   endif
	endif

	if(time.gt.5.and.onetry/=1)then
	   i=time-5
	   loop_problem=.true.
	   do while(i<time)
	      if(ticktock(time-5).ne.ticktock(i))then
		 loop_problem=.false.
		 exit
	      endif
	      i=i+1
	   enddo
	   if(loop_problem)then
	      timebomb=2
	      write(*,*)'Switching to ordered searches'
	      return
	   endif
	endif
	
	timebomb=1
	return
	end
!______________________________
	subroutine input(num_ells)

	use var
	implicit none
	integer::state,num_ells,seed,irandom
	integer::i,j,k,m,ii,column,io=6
	real::tempnum
        double precision::drandom
	real:: sh_min,sh_max,sh_av
        character(10)::str
        type(lattice),dimension(:),allocatable::tempmesh
        integer::numpoints
        integer,allocatable,dimension(:)::alist

	real::dummy1, dummy2, dummy3


! Read active.list to grab # of nuclei and corresponding nuclei IDs
! A grain may contain a spin ID=0 but we renumber from 1..Q
	open(1,file='active.list',iostat=io,status='old',err=943)
        if(io/=0)stop 'error opening active.list'

	read(1,*)ngrains
        write(*,*)'System will have ',ngrains,'grains'

        allocate(alist(ngrains))

        i=1
        do while(io==0)
           read(1,*,iostat=io,end=943)alist(i)
	   read(1,*)  ! number of inactive cells
           read(1,*)  ! inactive cell IDs
           i=i+1
        enddo
 943    close(1)

        call system("wc -l elliptical.cells | awk '{print $1}'>bugs")
        open(3,file="bugs")
        read(3,*)numpoints ! Six lines per entry
        write(*,*)'There are ',numpoints/6,'cells in elliptical.cells'
        close(3,status='delete')        
       
        allocate(mesh(ngrains),tempmesh(0:numpoints/6),stat=state)

! READ data into temporary structure since elliptical cells
! contains centers for ALL points not just grains
! Start counter at zero since spin=0 is the first entry

        open(2,file='elliptical.cells',iostat=io,status='old')
        if(io/=0)stop 'error opening elliptical.cells'
 
        i=0
        do while(io==0)
	   read(2,*,iostat=io,end=100) (tempmesh(i)%center(j),j=1,3)
           read(2,*,end=100) (tempmesh(i)%shape(j),j=1,3)
           read(2,*)! Sample Axis 1 vector
           read(2,*)! Sample Axis 2 vector
           read(2,*)! Sample Axis 3 vector
           read(2,*)! Random #'s
           i=i+1
        enddo          
 100    close(2)

! Define the true spin() and mesh() arrays by renumbering
! spin ids to be consecutive
	sh_max = 0.0
	sh_min = 99999.0
	sh_av = 0.0
        do i=1,ngrains
          mesh(i)%center=tempmesh(alist(i))%center
          mesh(i)%shape=tempmesh(alist(i))%shape
	  sh_max=amax1(sh_max,mesh(i)%shape(1),mesh(i)%shape(2),
     &    mesh(i)%shape(3))
	  sh_min=amin1(sh_min,mesh(i)%shape(1),mesh(i)%shape(2),
     &    mesh(i)%shape(3))
	  sh_av = sh_av + mesh(i)%shape(1)+mesh(i)%shape(2)+
     &    mesh(i)%shape(3)
        enddo

        deallocate(alist,tempmesh)

	if(iargc()/=7)stop 'USAGE: ./myCA Rx Ry Rz Nx Ny Nz periodic'
	                   !'periodic'
c       Read the data for the RVE box size from the command line
        call getarg(1,str)
        read(str,*)rve_x
        call getarg(2,str)
        read(str,*)rve_y
        call getarg(3,str)
        read(str,*)rve_z

!       Read the data for the Voxel definition from the command line
	call getarg(4,str)
        read(str,*)xmax
        call getarg(5,str)
        read(str,*)ymax
        call getarg(6,str)
        read(str,*)zmax

	call getarg(7,str)
	read(str,*)periodic
c
	write(*,*) 'RVE size: ', rve_x, rve_y, rve_z
	write(*,*) 'Voxel size: ', xmax, ymax, zmax
	write(*,*) 'Periodic boundary conditions are: ', periodic
	

	! this loop will continue until the user enters the correct 
	! box dimensions which maintain the RVE box ratios
!	do 
	
c$$$	   write(*,*)'RVE box dimensions were',rve_x,rve_y,rve_z
c$$$	   write(*,*)'The CA grid must maintain the same ratio'
c$$$	   write(*,*)'as the RVE box such that'
c$$$	   write(*,*)'CA_x      CA_y     CA_z  '
c$$$	   write(*,*)'----   = -----  = ------ '
c$$$	   write(*,*)'RVE_x     RVE_y    RVE_z '
c$$$	   write(*,*)
c$$$	   write(*,*)'What box dimensions would you like (x,y,z)?'
c$$$	   read(*,*)xmax,ymax,zmax
c$$$	   write(*,*)'You have chosen dimensions',xmax,ymax,zmax
	   
	   mag=1.0D0*xmax/rve_x
	   
!	   if((1.0D0*xmax/rve_x) == (1.0D0* ymax/rve_y) .and. 
!     &   (1.0D0*ymax/rve_y) == (1.0D0*zmax/rve_z)) exit

	   dummy1 = 1.0D0*xmax/rve_x
	   dummy2 = 1.0D0* ymax/rve_y
	   dummy3 = 1.0D0*zmax/rve_z

!	   write(*,*)dummy1 ,dummy2, dummy3
	   
!	   write(*,'(a)')'Scaling does not match between' 
!	   write(*,*)'RVE box and VOXEL domain'
!	   write(*,'(a)')'Please select a box size the '
!           write(*,*)'maintains the ratios X:Y:Z'
	   
!	enddo
	

	state=1
	allocate(spins(xmax,ymax,zmax),ell_volume(ngrains),stat=state)
	if(state/=0)stop 'error allocating spins array'

c$$$	write(*,*)'Do you want to enforce periodic Boundary conditions?'
c$$$	write(*,*)'1=Yes 2=No'
c$$$	read(*,*)state
c$$$	if(state==1)then
c$$$	   periodic=.true.
c$$$	elseif(state==2)then
c$$$	   periodic=.false.
c$$$	else
c$$$	   write(*,*)'Invalid response, we will automatically',
c$$$     &'set periodic=.true.'
c$$$	   periodic=.true.
c$$$	endif

! Initialize variables and RNG
	seed=ngrains ! RNG seed set equal to # of nuclei
	call irandom_seed(seed)

	num_ells=ngrains
        m3=xmax*ymax*zmax
        num_grown=0
	spins=0


! Convert ellipsoid centers from REAL to INTEGER values
! Using anint() instead of int() to round to nearest whole #
	do m=1,ngrains
	   i=anint(mag*mesh(m)%center(1))	   
	   if(i.lt.1) i=1
	   if(i.gt.xmax) i=xmax
	   j=anint(mag*mesh(m)%center(2))
	   if(j.lt.1) j=1
	   if(j.gt.ymax) j=ymax
	   k=anint(mag*mesh(m)%center(3))
	   if(k.lt.1) k=1
	   if(k.gt.zmax) k=zmax
	   if(spins(i,j,k)/=0)then
	      write(*,*)'site ',i,j,k, 'already filled with spin=',
     &	   spins(i,j,k)
	   else
	      spins(i,j,k)=m
	   endif
	   !write(*,*)'spin ID=',m,'with nuclei center=',i,j,k
	   num_grown=num_grown+1
 	enddo		!  do m=1,ngrains

! Determine delt_t
! The Shape Average=sum of all shapes along 3 axes divided by (3*#ellipsoids)
	sh_av = sh_av / (3.*num_ells)
	write(*,*)'average shape = ',sh_av
	write(*,*)'sh_min, sh_max = ',sh_min, sh_max
!	delta_t = 0.01 / sh_min  !  much smaller <= checking all 18 NN
!	delta_t = 0.05  / sh_max  !  make the increase in size slightly smaller than 1 pixel
!	delta_t=0.1
!  sh_max is the largest increment that any grain will undergo according to
!  the line below that reads
!            ell_size(j,k) = ell_shape(j,k) * float(i)*delta_t
!  therefore it is logical to limit the time/growth increment 
!  to a few voxels in the largest dimension

	delta_t = amin1((1./sh_max/real(xmax)) , (0.01 / sh_min)) 

	write(*,*)'this is a predetermined delta_t = ',delta_t

	write(*,*)'Would you like to modify delta_t yourself?'
	write(*,*)'1=Yes, 2=No'
	read(*,*)i

	if(i==1)then
	   write(*,*)'assign a value from 0.05 to 0.5'
	write(*,*)'Best to assign a small value and work towards larger'
	   read(*,*)delta_t
	endif

	return

 	end subroutine input

!__________________________
	subroutine growth(distance_nn,index)
!
! Looks for smallest distance to center from the availalble list

	use var
	implicit none
	
	real:: dist_min,distance_nn(18)
	integer::index,m

	dist_min = 999999.
	index = 0
	do m = 1,18
	   if(distance_nn(m).lt.dist_min.and.
     &		distance_nn(m).lt.1.0) then
	      dist_min = distance_nn(m)
	      index = m
	   endif
	enddo

	return
	end subroutine growth
