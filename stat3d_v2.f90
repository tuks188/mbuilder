! AUTHORS:  Chris Roberts  12/28/06
!           Anthony Rollett

! modified 7/18/07 because a 500^3 volume caused a seg fault
! when allocating the size_x size_y and size_z arrays
! Does not execute size() or ellipsoid() routines

! OVERVIEW:
! Program will extract ellipsoid information from
! discrete digital microstructures using moments of inertia.
! The grain ID's begin at ZERO not ONE to maintain consistency with
! other programs within the CABuilder package.

! INPUT: discrete structure on a simple cubic grid
! OUTPUT: stats file with ellipsoid ID and corresponding geometric properties
!         XML file format for annealing sequence

! KNOWN ISSUES
! If the libraries, g2c and lapack, are not in standard locations
! (i.e. /usr/lib/), then you will need to explicitly define their locations 
! by adding the -L/DIR flag during compilation

!Ex: 
! To find the location of libg2c.so use a syntax such as
! $ find /usr -name "libg2c*" -print
! /usr/lib/gcc-lib/i386-redhat-linux/3.2.2/libg2c.a
! /usr/lib/gcc-lib/i386-redhat-linux/3.2.2/libg2c.so
! /usr/lib/libg2c.so.0.0.0
! /usr/lib/libg2c.so.0

! Then, add an -L/DIR/ compiler flag to your command line
! -L/usr/lib/gcc-lib/i386-redhat-linux/3.2.2

! Must have a Fortran90 compiler and LAPACK libraries installed
! To compile:  (F90 compiler) -o stat3d stat3d.f90  -llapack -lg2c -L/usr/lib/gcc-lib/i386-redhat-linux/3.2.2
!
! To run:  ./stat3d INPUTFILE
!_________________________________________________________
      module global
      implicit none

! Ellipsoid Structure  
      type ellipsoid
      integer:: count
      real:: xc,yc,zc
      real:: euler(3)
      real:: F,G,H,A,B,C,A0,B0,C0
      real:: a01,b01,c01,a1,b1,c1
      real:: volume
      end type ellipsoid

      type(ellipsoid),allocatable,dimension(:)::shape


!Does structure have periodic boundary conditions
      logical,parameter::periodic=.false.

      integer,parameter::nset=200
! NEIGHBOR STORAGE ARRAYS
      type gstruct
      integer,dimension(nset)::nid
      integer,dimension(nset)::element
      integer::totalneighs
      end type

      type(gstruct),dimension(:),allocatable::grain
! Grain Boundary Area Array
  ! MPP=microns per pixel -- relationship between digital
  !     microstructure length scale and experimental length scale
      double precision::mpp

      integer,allocatable,dimension(:,:)::gbarea
      integer::total_gbarea
! # of Neighbors to search
! Prefer to search to 1st NN      --> nbors=6
! May be modified to 1st&2nd      --> nbors=18
!                    1st,2nd,&3rd --> nbors=26
      integer, parameter:: nbors=26


!part_id --> particle id, not used
!minsize --> cutoff value for a grain to be counted, not used
      integer,parameter::part_id=0,minsize=1

! CONSTANTS
      real, parameter::TOL=5.000,ca=0.0873
      double precision::PI

! File Characteristics
! Q = largest spin ID
! grains = total # of unique grains
      integer:: grains,Q

      logical,dimension(:,:,:),allocatable::spinsburn
      integer,dimension(:,:,:),allocatable::spins
      integer,dimension(:,:,:),allocatable::id
      integer,dimension(26,3)::nnbors
      double precision:: t
      character:: keyword*6
      real:: ave_rad,temp
      real:: radius,linevol,nucvol

! Dimensions of Structure
      integer:: mx,my,mz,m2,m3
      
      end module global
!_______________________________________________________________
      program stat3d

      use global
      implicit none
      character(len=30) :: fname,aspectfile
      integer::ngrain,i,alloc=10
      integer::rate,t1,t2
      integer::j,k,flag

      PI=dacos(-1.0D0)

! General Questions
      write(*,*)'Periodic Boundary Conditions=',periodic
      write(*,*)"If you need to change the boundary conditions edit the logical parameter: periodic and recompile"

      write(*,*)'What is micron/pixel ratio?'
      write(*,*)'If unknown, enter 1.0'
      read(*,*) mpp
!      mpp=1.0D0
      mpp=mpp*1.0D0 ! ensures double precision is assigned and stored


! Defining local neighborhood
      call neighbor_config()
      
! Obtaining filename from command line
      if(iargc()/=1)stop 'No filename given on command line'
      call getarg(1,fname)
! Read Data from File
      call system_clock(count=t1)
      call readDATA(fname)
      call system_clock(count=t2,count_rate=rate)
      write(*,*)'readData=',(t2-t1)/(1.0*rate),'seconds'

! Assign unique spin ID from 0 to ????
      write(*,*)'Do you need to uniquely ID each grain within the structure? 1=yes 2=no'
      read(*,*)flag
      if(flag==1)call unique_id()

      allocate(shape(0:Q),stat=alloc)
      if(alloc/=0)stop 'Error allocating shape array'
      do i=0,Q
         shape(i)%count=0
         shape(i)%xc=0.0
         shape(i)%yc=0.0
         shape(i)%zc=0.0
         shape(i)%F=0.0
         shape(i)%G=0.0          
         shape(i)%H=0.0
         shape(i)%A=0.0
         shape(i)%B=0.0
         shape(i)%C=0.0
      enddo

!Tally up the volume of each grain
      do i=1,mx
         do j=1,my
            do k=1,mz
               shape(spins(i,j,k))%count=shape(spins(i,j,k))%count+1
            enddo
         enddo
      enddo

      open(1,file='tempvol.txt',status='replace')
      write(1,'(a)')'GID   VOLUME'
      do i=0,Q
         write(1,*)i,shape(i)%count
      enddo
      close(1)
     
! Rewrite Data file with updated spin array
      if(flag==1)then
         call system_clock(count=t1)
         call rewrite(fname)
         call system_clock(count=t2,count_rate=rate)
         write(*,*)'rewrite=',(t2-t1)/(1.0*rate) ,'seconds'
      endif

! Determines gb area among all grains without using any arrays()
      call unique_id2()

! Calculate center of mass, volume, and translated coordinates
      !call system_clock(count=t1)
      !call size(aspectfile,ngrain)
      !call system_clock(count=t2,count_rate=rate)
      !write(*,*)'size=',(t2-t1)/(1.0*rate) ,'seconds'

      !if(Q/=ngrain-1)stop 'Ngrain not equal to Q'

! Calculate semi-axes and orientation
      !call system_clock(count=t1)      
      !call FitEllipsoids(aspectfile)
      !call system_clock(count=t2,count_rate=rate)
      !write(*,*)'FitEllipsoids=',(t2-t1)/(1.0*rate) ,'seconds'

! Determine GB area shared between (Si,Sj) spin pairs
      !call system_clock(count=t1)
      !call check_neighborhood()
      !call neighs_with_faces()
      !deallocate(spins,gbarea)
      !call system_clock(count=t2,count_rate=rate)
      !write(*,*)'Neighborhood routines=',(t2-t1)/(1.0*rate) ,'seconds'

      deallocate(spins,spinsburn)

! Output in XML file format which will be used in annealing procedure     
!      call xml_file_format()
      call xml_file_format_supersized()


      deallocate(shape)
      end
!______________________________________________________
      subroutine readDATA(fname)
!
! Reads the standard MC format data file & stores values in the spins array
!      
      use global
      implicit none
      integer::good=0
      integer::isite,jsite,ksite
      character(len=30)::fname

      open(1,file=fname,status='old',iostat=good)
      if(good/=0)stop 'error opening file given on commandline'
      write(*,*)'Reading data from ',fname
      read(1,*) mx,my,mz                  !mx,my,mz
      read(1,*) keyword,t,ave_rad,temp,q  !header
      read(1,*) radius,linevol,nucvol     !header

      m3=mx*my*mz

      write(*,*)'Initial Q = ',q
! Allocate memory and initialize to zero
      good=10
      allocate(spins(mx,my,mz),spinsburn(my,my,mz),stat=good)
      if(good.gt.0)stop 'allocation issue'

      spins=0
      spinsburn=.false.

! READ in spin data
      do ksite=1,mz
         do  jsite=1,my
            read(1,654,end=18)(spins(isite,jsite,ksite),isite=1,mx)
         enddo
      enddo
 18   close(1)
 654  format(20i6)

      return
      end
!______________________________________________________
      subroutine rewrite(fname)
!
! Outputs new spins array in standard MC format
!      
      use global
      implicit none
      integer::good=0
      integer::isite,jsite,ksite,samax=0
      character(len=30)::fname

      open(1,file=fname,status='replace',iostat=good)
      if(good/=0)stop 'error opening file'
      write(*,*)'Writing data to ',fname

      write(1,989) mx,my,mz
      write(1,990) keyword,t,ave_rad,temp,q
      write(1,991) radius,linevol,nucvol,samax,q

! WRITE spin data to file
      do ksite=1,mz
         do  jsite=1,my
            write(1,655)(spins(isite,jsite,ksite),isite=1,mx)
         enddo
      enddo
      close(1)
 655  format(20i6)
 989  format(3i8)
 990  format("'",a6,"'",f19.2,f7.3,f5.1,i10)
 991  format(f6.3,f6.3,f6.3,2(1x,i10))

      return
      end
! ___________________________________

      subroutine FitEllipsoids(aspectfile)
!
!   Program opens aspectXY.'keyword'file
!   Data is in TRANSLATED coordinates to handle periodic B.C.
!   This implies the true center of mass coordinates may be given
!   at points that exceed the box dimensions. 
!   --> Subroutine untranslate() will correct these values at the end
!       of the ellipsoid calculations
!   Fits ellipsoids to each grain
!   Provides volume, semi-axes, center, rotation angle
!
      use global
      implicit none

      integer:: x,y,z,num,i,j,k,filestat=1
      character:: aspectfile*20,line*30
      integer:: nPart=0,tally,gid
      real::a2,b2,c2
      real::ri,rj,rk
      double precision,dimension(3)::vec

! LAPACK VARIABLES
      integer n,order,info,lwork
      parameter(n=3,order=3)
      parameter(lwork=8)
      REAL*8 MAT(3,3),MATV(3,3),TL,eig(3)
      real*8 work(lwork)
! END LAPACK

      write(*,*)'Entering FitEllipsoids'

! Open Data File and extract info
      open(1,file=aspectfile,iostat=filestat)  
      if(filestat/=0)stop 'Error opening file'

 990  format(a20)

! Determine center of mass for each particle
!                msubi * xsubi
! xc= summation  --------------   where M=total mass of particle
!                      M
!
!   COMPLETED in size.f to reduce redundant calculations

! Determine product moments[F,G,H] and 2nd moments of inertia [A,B,C]
! where C>B>A
! using center of mass as origin of coordinate system
! and applying Parallel Axis Theorem
! 
! A = summation   msubi*(ysubi - yc)^2 + (zsubi - zc)^2
! B = summation   msubi*(zsubi - zc)^2 + (xsubi - xc)^2
! C = summation   msubi*(xsubi - xc)^2 + (ysubi - yc)^2

! F = summation   msubi*(ysubi-yc) * (zsubi-zc)
! G = summation   msubi*(zsubi-zc) * (xsubi-xc)
! H = summation   msubi*(xsubi-xc) * (ysubi-yc)

      do 302 while(filestat==0)
         read(1,990,iostat=filestat,end=981)line

         if(line(2:2).ne.'#')then
! Extract integer values from character string
! changed shape(npart) to shape(num) since each array index
! corresponds exactly to the spin ID
            read(line,*)num,x,y,z
            ri=real(x)          !Prevent any arithmetic mistakes
            rj=real(y)
            rk=real(z)
                                ! 2nd Moments
            shape(num)%A=shape(num)%A+((rj-shape(num)%yc)**2 + (rk-shape(num)%zc)**2)
            shape(num)%B=shape(num)%B+((rk-shape(num)%zc)**2 + (ri-shape(num)%xc)**2)
            shape(num)%C=shape(num)%C+((ri-shape(num)%xc)**2 + (rj-shape(num)%yc)**2)
                                ! PRODUCT MOMENTS
            shape(num)%F=shape(num)%F+(rj-shape(num)%yc)*(rk-shape(num)%zc)
            shape(num)%G=shape(num)%G+(rk-shape(num)%zc)*(ri-shape(num)%xc)
            shape(num)%H=shape(num)%H+(ri-shape(num)%xc)*(rj-shape(num)%yc)
         endif
 302  enddo

 981  close(1)

      do i=0,Q

! Create Inertia Matrix
!      |   A    -H   -G  |
!  I = |  -H     B   -F  |
!      |  -G    -F    C  |
!
         mat(1,1)=shape(i)%A
         mat(1,2)=-1.0*shape(i)%H
         mat(1,3)=-1.0*shape(i)%G
         mat(2,1)=-1.0*shape(i)%H
         mat(2,2)=shape(i)%B
         mat(2,3)=-1.0*shape(i)%F
         mat(3,1)=-1.0*shape(i)%G
         mat(3,2)=-1.0*shape(i)%F
         mat(3,3)=shape(i)%C

! Determines eigenvalues of matrix
! also known as principal moments (A0,B0,C0)
! eig(1)=A0   eig(2)=B0  and eig(3)=C0 where C0>B0>A0

! matrix mat is returned and contains the direction cosines
! for each of the principal axes (i.e. Rotation Matrix)   

         call dsyev('V','U',order,mat,n,eig,work,lwork,info)
!         call la_syev(mat,eig,'V','U',info)
!
!        M (a^2 + b^2)
! C0 =   ------------  where M=total mass of ellipsoid
!            5.0             a & b are semi-axes 
!        M (c^2 + a^2)
! B0 =   ------------  where M=total mass of ellipsoid
!            5.0             c & a are semi-axes
!        M (b^2 + c^2)
! A0 =   ------------  where M=total mass of ellipsoid
!            5.0             b & c are semi-axes
!         
! A,B,C or A0,B0,C0 can be used to calculate semi-axes
! but reference frame will be different
!
! Safety Check:       A+B+C=A0+B0+C0
         
         if((shape(i)%A+shape(i)%B+shape(i)%C)-(eig(1)+eig(2)+eig(3)) > TOL)then
            write(*,*)'Error in A0B0C0 calculation with spin=',i
            stop
         endif

         if(i<10)then
            write(*,*)'A0 B0 C0 ',eig
            write(*,*)'A  B   C ',shape(i)%A,shape(i)%B,shape(i)%C
            write(*,*)shape(i)%xc,shape(i)%yc,shape(i)%zc
         endif

! Calculate Semi-axes along principal moments (A0,B0,C0)

         shape(i)%A0=eig(1)
         shape(i)%B0=eig(2)
         shape(i)%C0=eig(3)
         a2=5.0/(2.0*shape(i)%count)*(eig(3)+eig(2)-eig(1))
         b2=5.0/(2.0*shape(i)%count)*(eig(1)-eig(2)+eig(3))
         c2=5.0/(2.0*shape(i)%count)*(eig(2)-eig(3)+eig(1))
         shape(i)%a01 = sqrt(a2)
         shape(i)%b01 = sqrt(b2)
         shape(i)%c01 = sqrt(c2)
         if(i<10)write(*,*)'a1 b1 and c1 using (A0B0C0)',shape(i)%a01,shape(i)%b01,shape(i)%c01 

! Calculate semi-axes w.r.t. sample axes
         a2=5.0/(2.0*shape(i)%count)*(shape(i)%C+shape(i)%B-shape(i)%A)
         b2=5.0/(2.0*shape(i)%count)*(shape(i)%A+shape(i)%C-shape(i)%B)
         c2=5.0/(2.0*shape(i)%count)*(shape(i)%A+shape(i)%B-shape(i)%C)
         shape(i)%a1 = sqrt(a2)
         shape(i)%b1 = sqrt(b2)
         shape(i)%c1 = sqrt(c2)
         if(i<10)write(*,*)'a1 b1 and c1 using(ABC)',sqrt(a2),sqrt(b2),sqrt(c2)

! Calculate Ellopsoid Volume=4/3*PI*a*b*c
! Discrete Voxel Volume is contained within shape(i)%count
         shape(i)%volume=4.0/3.0*PI*shape(i)%a1*shape(i)%b1*shape(i)%c1

! Convert rotaion matrix(eigenvector matrix) to Euler angles (units=radians)
         call mat2Euler(mat,vec)
         do k=1,3
            shape(i)%euler(k)=vec(k)
         enddo
      
      enddo

      call untranslate() ! Corrects any center of mass coordinates outside box

      return
      end subroutine FitEllipsoids
!_____________________________________
      subroutine untranslate()
! If center of mass is less than lower bound, add length of box
! If center of mass is greater than upper bound, subtract length of box
      use global
      implicit none
      
      integer::i

      do i=0,Q

         if(shape(i)%xc>mx)shape(i)%xc=shape(i)%xc-real(mx)
         if(shape(i)%xc<1)shape(i)%xc=shape(i)%xc+real(mx)

         if(shape(i)%yc>my)shape(i)%yc=shape(i)%yc-real(my)
         if(shape(i)%yc<1)shape(i)%yc=shape(i)%yc+real(my)

         if(shape(i)%zc>mz)shape(i)%zc=shape(i)%zc-real(mz)
         if(shape(i)%zc<1)shape(i)%zc=shape(i)%zc+real(mz)

      enddo
      return
      end

!____________________________________
      subroutine neighbor_config()
! Loads values for out to 3rd NN
! Only uses nnbors list out to value specified by parameter, nbors
! Taken from rex3d code

      use global
      implicit none
      integer::i,j

! Initialize to zero
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
!__________________________________________________________
! nnbors are global values initialized in neighbor_config
! written by A.D. Rollett
!
      subroutine neighs(isite,jsite,ksite,nbor,inbr,jnbr,knbr)
      use global
      implicit none
      integer:: isite,jsite,ksite,nbor,inbr,jnbr,knbr

      inbr=isite+nnbors(nbor,1)
      jnbr=jsite+nnbors(nbor,2)
      knbr=ksite+nnbors(nbor,3)

      if(periodic)then
         inbr=mod((inbr+mx-1),mx)+1
         jnbr=mod((jnbr+my-1),my)+1
         knbr=mod((knbr+mz-1),mz)+1
      else
         if(inbr>mx)inbr=mx
         if(jnbr>my)jnbr=my
         if(knbr>mz)knbr=mz
         if(inbr<1)inbr=1
         if(jnbr<1)jnbr=1
         if(knbr<1)knbr=1
      endif

      return
      end subroutine neighs
!___________________________________________________
!        SIZE.F: Measures grain volumes via a burn algorithm.
!       ---------------------------------------------------------------
!  written by A.D. Rollett
!  modified by C.G. Roberts  12/06
        subroutine size(fname2,ngrain)
        use global
        implicit none

        integer::ngrain
        integer tail,in,jn,kn,is,js,ks
        integer nbor,head,spin,area,neigh,ot,ot2
        real:: sumrbar,rbar,abar
        integer::i,j,k
        integer*8 lsite,nucl,list(m3),site
        integer:: count,expon,small
        character ctime*10,number*3
        character fname*20,fname2*20

        integer samax,spinamax
!       samax is the maximum area, SPINAMAX is its spin
! NEW TRANSLATION VARIABLES

        integer::size_x,size_y,size_z,x_curr,y_curr,z_curr
        integer,allocatable,dimension(:)::list_x,list_y,list_z
        integer::alloc=10

! Initialize variables to zero
        write(*,*) 'entering size.f'
        samax=0
        spinamax=-1
        ot=1000
        ot2=ot**2
        abar=0.
        rbar=0.
        sumrbar=0.
        tail=0
        head=0
        ngrain=0
        small=0        

        size_x = mx/2
        size_y = my/2
        size_z = mz/2
        allocate(list_x(m3),list_y(m3),list_z(m3),stat=alloc)
        if(alloc/=0)stop 'Error allocating list_x array'

! convert time to character
        if(t<0.0D0)t=0.0D0 !safety check
        write(ctime,'(i10)')int(t)
        do i=10,1,-1
          if(ctime(i:i)==" ")then
             goto 567
          endif
        enddo
!        adjustl(ctime)

! Open Output files
567     fname='size'//ctime(i+1:10)//'.'//keyword
        fname2='aspect'//ctime(i+1:10)//'.'//keyword

        write(*,*)fname,fname2,size_x,size_y,size_z
        open(4,file=fname2,err=599)       
        open(10,file=fname,err=599)
        write(10,*)'GID   VOLUME'

        do 100 i=1,mx
           do 90 j=1,my
              do 80 k=1,mz
                 nucl=ot2*i+ot*j+k

!        if the nucleus has been burned previously, or is a pore
!        or second phase, go to next nucleus
! modified to include spin ID=0

                 if(spins(i,j,k)<0) go to 80

!        find the spin of the new nucleus, then burn it
!        and set the area of the new cluster equal to one

                 spin=spins(i,j,k)
                 write(4,*)spin,i,j,k
                 !write(*,*)'Grain = ',spin, 'found at ',i,j,k
                 spins(i,j,k)=-1*spins(i,j,k)-q
                 area=1

!        set the current tail and head position and let 'nucl' be
!        the tail

                 tail=tail+1
                 head=head+1
                 list(tail)=nucl
                 
                 list_x(tail) = i
                 list_y(tail) = j
                 list_z(tail) = k

!        let the tail of the list be the current site, 'site'
!        (first pass, site=nucl; later passes, site=next site
!        in the current cluster)

 50              site=list(tail)
                 is=int(site/ot2)
                 js=int(site/ot)-is*ot
                 ks=int(site-is*ot2-js*ot)
                 

                 x_curr = list_x(tail)
                 y_curr = list_y(tail)
                 z_curr = list_z(tail)
               
!        check all the neighbors 'neigh' of the site:  if 'neigh'
!        belongs to the cluster, increment the cluster area,
!        make 'neigh' the head of the list, and burn 'neigh'

                 do 60 nbor=1,nbors
                    call neighs(is,js,ks,nbor,in,jn,kn)
                    if(spins(in,jn,kn).eq.spin) then
                       neigh=in*ot2+jn*ot+kn
                       head=head+1
                       area=area+1
                       list(head)=neigh
                       spins(in,jn,kn)=-1*spins(in,jn,kn)-q
                       
!     now for the coordinate translation
                       list_x(head) = in
                       list_y(head) = jn
                       list_z(head) = kn

                       if(iabs(x_curr - in).gt.size_x) then
                          if((x_curr - in).gt.size_x) then
                             list_x(head) = in + mx
                          endif
                          if((in - x_curr).gt.size_x) then
                             list_x(head) = in - mx
                          endif
                       endif

                       if(iabs(y_curr - jn).gt.size_y) then
                          if((y_curr - jn).gt.size_y) then
                             list_y(head) = jn + my
                          endif
                          if((jn - y_curr).gt.size_y) then
                             list_y(head) = jn - my
                          endif
                       endif
                       

                       if(iabs(z_curr - kn).gt.size_z) then
                          if((z_curr - kn).gt.size_z) then
                             list_z(head) = kn + mz
                          endif
                          if((kn - z_curr).gt.size_z) then
                             list_z(head) = kn - mz
                          endif
                       endif
 543                   continue

                       write(4,*)spin,list_x(head),list_y(head),list_z(head)

! Store center of mass for later use in fit_ellipsoids

                       shape(spin)%xc=shape(spin)%xc+real(list_x(head))
                       shape(spin)%yc=shape(spin)%yc+real(list_y(head))
                       shape(spin)%zc=shape(spin)%zc+real(list_z(head))

                    endif       ! if (spins(in,jn,kn) == spin)
                    

 60              enddo ! do nbor=1,nbors
     
!        If the tail of the list does not equal the head of the list, then:
!        sites were added to the cluster in the last pass, so all the sites
!        in the cluster have not been checked for neighbors.  So, check the
!        unchecked site closest to the tail of the list.

                 if(tail.ne.head) then
                    tail=tail+1
                    go to 50
                 else
                    write(4,*)'#' !END of each grain
                    
                    shape(spin)%xc=shape(spin)%xc/area
                    shape(spin)%yc=shape(spin)%yc/area                    
                    shape(spin)%zc=shape(spin)%zc/area
                    shape(spin)%count=area
     
!        If the tail = the head, the cluster has been completely burned.
!        Record the area, and go to the next cluster.

                    if(area.ge.minsize)then
                       ngrain=ngrain+1
                       sumrbar=sumrbar+(float(area))**(.33333)
                       write(10,*) spin,area
                    else
                       small=small+1
                    end if

                    if(area.gt.samax) then
                       samax=area
                       spinamax=spin
                    endif
                 endif
                 
 80           enddo
 90        enddo
 100    enddo

!       tabulate results
!
        abar=(1.0*m3-lsite)/(1.0*ngrain)
        rbar=sumrbar/(1.0*ngrain)
        write(*,*)'# of grains=  ',ngrain
        write(*,*)'# of small grains= ',small

! Return Spins Array to original state
        do i=1,mx
           do j=1,my
              do k=1,mz
                 spins(i,j,k)=-1*(spins(i,j,k)+Q)
              enddo
           enddo
        enddo


 599    continue
!        write(4,573)'END'
        close(4)
        close(10)
        deallocate(list_x,list_y,list_z)
        return
 573    format(a3)
        end

!__________________________________________

        subroutine unique_id()
          !  written by A.D. Rollett to ensure each grain has a unique spin ID
          !  modified by C.G. Roberts  12/06

          use global
          implicit none

          integer::in,jn,kn,is,js,ks
          integer::i,j,k,myid
          integer*8 lsite,nucl,list(m3),site
          integer:: nbor,head,tail,spin,neigh,ot,ot2
          integer::numneighs,ii
          type(gstruct)::neighbors
          logical::exist

          !  samax is the maximum area, SPINAMAX is its spin
          ! Initialize variables to zero
          write(*,*) 'entering unique_id.f of stat3d_v2.f90'

          open(2,file='neighlist.txt',status='replace')

          ot=1000
          ot2=ot**2
          tail=0
          head=0
          !Lowest spin value given to a grain
          ! (For C/C++, may need to set to 0, so that 1st Grain ID=0)
          ! (For Fortran, set myid=1, so that 1st Grain ID=1 ) <-- not working dont use
          myid=0

          do 100 i=1,mx
             do 90 j=1,my
                do 80 k=1,mz
                   nucl=ot2*i+ot*j+k
                   
                   neighbors%nid=0
                   neighbors%element=0
                   neighbors%totalneighs=0

                   !      if the nucleus has been burned previously, or is a pore
                   !        or second phase, go to next nucleus
                   !     Uses spin=0 as a spin ID
                   
                   if(spins(i,j,k)<0) go to 80
                   myid=myid+1
                   
                   !        find the spin of the new nucleus, then burn it
                   !        and set the area of the new cluster equal to one
                   !
                   spin=spins(i,j,k)
                   spins(i,j,k)=-1*myid
                   !
                   !        set the current tail and head position and let 'nucl' be
                   !        the tail
                   !
                   tail=tail+1
                   head=head+1
                   list(tail)=nucl
                   
                   !        let the tail of the list be the current site, 'site'
                   !        (first pass, site=nucl; later passes, site=next site
                   !        in the current cluster)
                   
50                 site=list(tail)
                   is=int(site/ot2)
                   js=int(site/ot)-is*ot
                   ks=int(site-is*ot2-js*ot)
                   
                   
                   !        check all the neighbors 'neigh' of the site:  if 'neigh'
                   !        belongs to the cluster, increment the cluster area,
                   !        make 'neigh' the head of the list, and burn 'neigh'
                   
                   do 60 nbor=1,nbors
                      call neighs(is,js,ks,nbor,in,jn,kn)
                      if(spins(in,jn,kn).eq.spin) then
                         neigh=in*ot2+jn*ot+kn
                         head=head+1
                         list(head)=neigh
                         spins(in,jn,kn)=-1*myid
                      endif
60                 enddo

                      !        If the tail of the list does not equal the head of the list, then:
                      !        sites were added to the cluster in the last pass, so all the sites
                      !        in the cluster have not been checked for neighbors.  So, check the
                      !        unchecked site closest to the tail of the list.

                   if(tail.ne.head) then
                      tail=tail+1
                      go to 50

                   endif
                   
                   !        If the tail = the head, the cluster has been completely burned.
                   !        Record the area, and go to the next cluster.
                   
80              enddo
90           enddo
100       enddo

          write(*,*)'Actual # of grains=  ',myid
          write(*,*)'while Q=', myid-1, 'since spin=0 is accepted'
          ! Set Q to myid since myid = max(spin value)
          Q=myid-1
          
          ! Renumbering Spin Array
          ! The subtraction of 1 causes the lowest GID to be set to 0
          ! as desired.
          do i=1,mx
             do j=1,my
                do k=1,mz
                   spins(i,j,k)=abs(spins(i,j,k))-1
                enddo
             enddo
          enddo
          
          
          
          return
        end subroutine unique_id



!__________________________________________

        subroutine unique_id2()
          !  written by A.D. Rollett to ensure each grain has a unique spin ID
          !  modified by C.G. Roberts  7/07
! Due to memory constraints the gbarea among grains is measured and outputted to file
! which is subsuquently opened and used in the XML file output subroutine


          use global
          implicit none

          integer::in,jn,kn,is,js,ks
          integer::i,j,k,myid
          integer*8 lsite,nucl,list(m3),site
          integer:: nbor,head,tail,spin,neigh,ot,ot2
          integer::numneighs,ii,area
          type(gstruct)::neighbors
          logical::exist

          !  samax is the maximum area, SPINAMAX is its spin
          ! Initialize variables to zero
          write(*,*) 'entering unique_id.f part dos of stat3d_v2.f90'

          open(2,file='neighlist.txt')
          open(3,file='tempvol_unique.txt')

          ot=1000
          ot2=ot**2
          tail=0
          head=0

          !Lowest spin value given to a grain
          ! (For C/C++, may need to set to 0, so that 1st Grain ID=0)
          ! (For Fortran, set myid=1, so that 1st Grain ID=1 ) <-- not working dont use
          myid=0

          do 100 i=1,mx
             do 90 j=1,my
                do 80 k=1,mz
                   nucl=ot2*i+ot*j+k
                   
                   do ii=1,nset
                      neighbors%nid(ii)=0
                      neighbors%element(ii)=0
                   enddo
                   neighbors%totalneighs=0

                   !      if the nucleus has been burned previously, or is a pore
                   !        or second phase, go to next nucleus
                   !     Uses spin=0 as a spin ID
                   
                   if(spinsburn(i,j,k)) go to 80
                   
                   !        find the spin of the new nucleus, then burn it
                   !        and set the area of the new cluster equal to one
                   !
                   spin=spins(i,j,k)
                   !write(*,*)'Working on spin ID=',spin,'and myid=',myid
                   myid=myid+1
                   spinsburn(i,j,k)=.true.
                   !
                   !        set the current tail and head position and let 'nucl' be
                   !        the tail
                   !
                   tail=tail+1
                   head=head+1
                   area=1
                   list(tail)=nucl
                   
                   !        let the tail of the list be the current site, 'site'
                   !        (first pass, site=nucl; later passes, site=next site
                   !        in the current cluster)
                   
50                 site=list(tail)
                   is=int(site/ot2)
                   js=int(site/ot)-is*ot
                   ks=int(site-is*ot2-js*ot)
                   
                   
                   !        check all the neighbors 'neigh' of the site:  if 'neigh'
                   !        belongs to the cluster, increment the cluster area,
                   !        make 'neigh' the head of the list, and burn 'neigh'
                   
                   do 60 nbor=1,nbors
                      call neighs(is,js,ks,nbor,in,jn,kn)
                      if(.not.spinsburn(in,jn,kn).and.spins(in,jn,kn).eq.spin) then
                         neigh=in*ot2+jn*ot+kn
                         head=head+1
                         list(head)=neigh
                         spinsburn(in,jn,kn)=.true.
                         area=area+1

                         ! If spin is different and one of the 1st NNs, then add to neighlist
                      elseif(spins(in,jn,kn)/=spin.and.nbor<=6)then

                         exist=.false.
                         if(neighbors%totalneighs>0)then
                            do ii=1,neighbors%totalneighs
                               if(neighbors%nid(ii)==spins(in,jn,kn))then
                                  neighbors%element(ii)= neighbors%element(ii)+1
                                  exist=.true.
                               endif
                            enddo
                         endif

                         if(.not.exist)then
                            neighbors%totalneighs=neighbors%totalneighs+1
                            neighbors%nid(neighbors%totalneighs)=spins(in,jn,kn)
                            neighbors%element(neighbors%totalneighs)=neighbors%element(neighbors%totalneighs)+1
                         endif
                      endif
60                 enddo

                      !        If the tail of the list does not equal the head of the list, then:
                      !        sites were added to the cluster in the last pass, so all the sites
                      !        in the cluster have not been checked for neighbors.  So, check the
                      !        unchecked site closest to the tail of the list.

                   if(tail.ne.head) then
                      tail=tail+1
                      go to 50
                   else
                      ! Output neigh list to file
                      write(2,*)spin
                      do ii=1,neighbors%totalneighs
                         write(2,*)neighbors%nid(ii),neighbors%element(ii)
                      enddo
                      write(2,'(a)')'#'

                      write(3,*)spin,area                      

                   endif
                   
                   !        If the tail = the head, the cluster has been completely burned.
                   !        Record the area, and go to the next cluster.
                   
80              enddo
90           enddo
100       enddo

          write(*,*)'Actual # of grains=  ',myid
          write(*,*)'while Q=', myid-1, 'since spin=0 is accepted'
          ! Set Q to myid since myid = max(spin value)
          Q=myid-1
          
          
          close(2)
          close(3)
          
          return
        end subroutine unique_id2

          !__________________________________________

          subroutine mat2Euler(mat,v)
            ! Converts a passive rotation matrix to Euler Angles (units=radians)
            ! As PHI approaches zero (cos(phi)-->1.0), phi1 and phi2 are dependent
            ! The TOL defines the threshold below which PHI is taken as 1.0

            implicit none
            integer,parameter::r8=selected_real_kind(10,200)
            integer,parameter::TOL=0.005
            real(kind=r8),dimension(3,3),intent(in)::mat
            real(kind=r8),dimension(3),intent(out)::v
            real(kind=r8)::sp
            real(kind=r8),parameter::dp=1.0D0

            v(2)=acos(mat(3,3))
            if(1.0D0-mat(3,3)<TOL)then
               v(2)=0.0D0
               v(1)=atan2(mat(1,2),mat(1,1))/(2*dp)
               v(3)=-1.0D0*v(1)
            else
               v(2)=acos(mat(3,3))
               sp=sin(v(2))
               v(3)=atan2(mat(1,3)/sp,mat(2,3)/sp)
               v(1)=atan2(mat(3,1)/sp,-dp*mat(3,2)/sp)
            endif

            return
          end subroutine mat2Euler
          !    ____________________________________________________

          subroutine check_neighborhood()

            ! Accummulates a listing of grain boundary area shared
            ! among the Q grains in the microstucture

            ! The area between (i,j)=(j,i)
            ! Total gbarea is 2X(actual value), but is needed to define
            ! cellidealization file.
            !
            use global
            implicit none
            integer::good=20
            integer::nx,ny,nz,ga,gb
            integer::i,j,k,in,jn,kn,n
            integer::nbr

            write(*,*)'Entering check_neighborhood'
            ! Allocate and initialize to zero
            !      allocate(gbarea(0:Q,0:Q),burnid(mx,my,mz),stat=good)
            allocate(gbarea(0:Q,0:Q),stat=good)

            if(good/=0)stop 'error allocating area array'
            gbarea=0
            total_gbarea=0

            !      do id=0,Q

            ! Searching 1st NN around each lattice site
            do i=1,mx
               do j=1,my
                  do k=1,mz
                     ga=spins(i,j,k)
                     if(ga>=0)then
                        ! Only search 1st NN
                        do n=1,6                   
                           call neighs(i,j,k,n,in,jn,kn)
                           gb=spins(in,jn,kn)
                           if(ga/=gb.and.gb>=0)then
                              total_gbarea=total_gbarea+1
                              gbarea(ga,gb)= gbarea(ga,gb)+1
                           endif
                        enddo

                     endif ! if(ga>=0)

                  enddo ! k=1,mz
               enddo ! j=1,my
            enddo !i=1,mx


            write(*,*)'total_gbarea=',total_gbarea
            write(*,*)'REMEMBER, this is 2X the actual value'

            return
            open(7,file='gbarea.txt')
            do i=0,Q
               do j=0,Q
                  write(7,*)i,j,gbarea(i,j)
               enddo
            enddo
            close(7)
            return
          end subroutine check_neighborhood

          !___________________________________________________________
          subroutine neighs_with_faces()

            ! written by C.G. Roberts
            ! Subroutine creates a struct for each grain and its corresponding
            ! neighbors & shared boundary areas
            use global
            implicit none

            integer::i,j
            integer::good=20
            integer::numneighs,maxnumneighs

            write(*,*)'Entering neighs_with_faces'
            allocate(grain(0:Q),stat=good)  
            if(good/=0)stop 'error allocating grain array'

            maxnumneighs=0
            do i=0,Q
               numneighs=0
               do j=0,Q
                  if(gbarea(i,j)/=0)then
                     numneighs=numneighs+1
                     if(numneighs>nset)then
                        write(*,*)'Exceeded grain array limits!! Increase the parameter value nset in the global module'
                        write(*,*)'at array indices ',i,j
                        stop
                     endif
                     grain(i)%nid(numneighs)=j
                     grain(i)%element(numneighs)=gbarea(i,j)
                  endif
               enddo
               if(numneighs>maxnumneighs)maxnumneighs=numneighs
               grain(i)%totalneighs=numneighs
            enddo

            return
          end subroutine neighs_with_faces
          !____________________________________________________________

          subroutine xml_file_format()
            ! written by Joe Fridy
            ! modified by C.G. Roberts
            ! Subroutine will create the xml file format for the final_anneal.c program
            !
            use global
            implicit none
            integer::i,j,tneighs,io=10

            write(*,*)'Entering xml_file_format'
            open(1,file="cellIdealization2.xml",status='replace',iostat=io)
            if(io/=0)stop 'Error opening cellIdealization file'

            write(1,*)"<nregions> ",Q+1 ," </nregions>\n" !need Q+1 since 0 spin is included

            do i=0,Q
               write(1,*)'<region>'
               write(1,*)'   <id>',i,' </id>'
               ! Volume = #Voxels * (microns/pixels)^3
               write(1,*)'   <volume> ',shape(i)%count*mpp**3,' </volume>'
               write(1,*)'   <center> ',shape(i)%xc,shape(i)%yc,shape(i)%zc,' </center>'
               write(1,*)'   <approximating_ellipsoid>'
               write(1,*)'      <ecenter> ',shape(i)%xc,shape(i)%yc,shape(i)%zc,' </ecenter>'

               write(1,*)'      <semi_axis_lengths> ',shape(i)%a01,shape(i)%b01,shape(i)%c01,' </semi_axis_lengths>'

               ! Axis 1 is < 1.0  0.0  0.0 >
               write(1,*)'      <axis1> 1.0  0.0 0.0 </axis1>'
               ! Axis 2 is < 0.0  1.0  0.0 >
               write(1,*)'      <axis2> 0.0 1.0 0.0 </axis2>'

               ! Axis 3 is < 0.0  0.0  1.0 >
               write(1,*)'      <axis3> 0.0 0.0 1.0 </axis3>'
               ! Euler angles given in RADIANS
               write(1,*)'      <orientation> ',shape(i)%euler(1),shape(i)%euler(2),shape(i)%euler(3),' </orientation>'

               write(1,*)'   </approximating_ellipsoid>'
               write(1,*)'   <npatches> ',grain(i)%totalneighs,' </npatches>'

               tneighs=grain(i)%totalneighs

               do j=1,tneighs
                  write(1,*)'   <patch>'
                  write(1,*)'      <neighbor> ',grain(i)%nid(j),' </neighbor>'

                  ! Just output at 1.0 1.0 1.0
                  write(1,*)'      <average_outward_normal> 1.0 1.0 1.0',' </average_outward_normal>'

                  ! Area given as Real # (Si,Sj) area * (microns/pixels)^2
                  write(1,*)'      <area> ',grain(i)%element(j)*mpp**2,' </area>'
                  ! Area given as # of pixels shared between i and j
                  write(1,*)'      <nfacets> ',grain(i)%element(j),' </nfacets>'
                  write(1,*)'   </patch>'
               enddo

               write(1,*)'</region>'
            enddo

            close(1)
            return
          end subroutine xml_file_format


!__________________________________________________________________

          subroutine xml_file_format_supersized()
            ! written by Joe Fridy
            ! modified by C.G. Roberts
            ! Subroutine will create the xml file format for the final_anneal.c program
            !
            use global
            implicit none
            integer::i,j,tneighs,io=10

            character::line*40
            type(gstruct)::neighs
            logical::flag,lineflag
            integer::gid,num,io2

            write(*,*)'Entering xml_file_format'
            open(1,file="cellIdealization2.xml",status='replace',iostat=io)
            if(io/=0)stop 'Error opening cellIdealization file'

            open(2,file='neighlist.txt',iostat=io2)
            if(io2/=0)stop 'error opening neighlist.txt'
            flag=.true.

            write(1,*)"<nregions> ",Q+1 ," </nregions>\n" !need Q+1 since 0 spin is included

            do i=0,Q

               lineflag=.true.
               num=0

               do while(lineflag)

                  read(2,567,iostat=io2,end=468)line
567               format(a40)
                  !write(*,*)line
                  
                  if(line(1:1)=='#')then
                     flag=.true.
                     lineflag=.false.
                  elseif(flag)then
                     read(line,*)gid
                     flag=.false.
                  elseif(.not.flag)then
                     num=num+1
                     read(line,*)neighs%nid(num),neighs%element(num)
                  endif
                  
               enddo
               

               write(1,*)'<region>'
               write(1,*)'   <id>',i,' </id>'
               ! Volume = #Voxels * (microns/pixels)^3
               write(1,*)'   <volume> ',shape(i)%count*mpp**3,' </volume>'
               write(1,*)'   <center> ',shape(i)%xc,shape(i)%yc,shape(i)%zc,' </center>'
               write(1,*)'   <approximating_ellipsoid>'
               write(1,*)'      <ecenter> ',shape(i)%xc,shape(i)%yc,shape(i)%zc,' </ecenter>'

               write(1,*)'      <semi_axis_lengths> ',shape(i)%a01,shape(i)%b01,shape(i)%c01,' </semi_axis_lengths>'

               ! Axis 1 is < 1.0  0.0  0.0 >
               write(1,*)'      <axis1> 1.0  0.0 0.0 </axis1>'
               ! Axis 2 is < 0.0  1.0  0.0 >
               write(1,*)'      <axis2> 0.0 1.0 0.0 </axis2>'

               ! Axis 3 is < 0.0  0.0  1.0 >
               write(1,*)'      <axis3> 0.0 0.0 1.0 </axis3>'
               ! Euler angles given in RADIANS
               write(1,*)'      <orientation> ',shape(i)%euler(1),shape(i)%euler(2),shape(i)%euler(3),' </orientation>'

               write(1,*)'   </approximating_ellipsoid>'
               write(1,*)'   <npatches> ',num,' </npatches>'

               

               do j=1,num
                  write(1,*)'   <patch>'
                  write(1,*)'      <neighbor> ',neighs%nid(j),' </neighbor>'

                  ! Just output at 1.0 1.0 1.0
                  write(1,*)'      <average_outward_normal> 1.0 1.0 1.0',' </average_outward_normal>'

                  ! Area given as Real # (Si,Sj) area * (microns/pixels)^2
                  write(1,*)'      <area> ',neighs%element(j)*mpp**2,' </area>'
                  ! Area given as # of pixels shared between i and j
                  write(1,*)'      <nfacets> ',neighs%element(j),' </nfacets>'
                  write(1,*)'   </patch>'
               enddo

               write(1,*)'</region>'
            enddo

468            close(1)
            close(2)

            return
          end subroutine xml_file_format_supersized
