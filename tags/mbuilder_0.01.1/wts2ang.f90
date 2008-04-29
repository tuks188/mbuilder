      module master

      implicit none
      integer::num_orients
      real,allocatable,dimension(:)::phi1,Phi,phi2
      integer,allocatable,dimension(:)::area

      end module
!_______________________________________________

      program wts2ang
!
!     AUTHOR:  Chris Roberts
!
!     Updates:
!        CGR v2 1/07  corrected weighting of orientations in ANG file
!        CGR v2 2/07  use dynamic memory allocation
!        CGR v2 6/07  calculates power of 10 based on minval in file'

!     Purpose: Adapts a wts file to the typical ang format for compatability with
!     the TEXSEM OIM 3.5 and newer software

!     Compile:  $(f90) wts2ang_v2.f90  -o exe
!     Run:  ./exe inputfile.wts output.ang (requires 2 command line arguments)

      use master

      implicit none
      integer:: bunny,i
      real:: pi,p1,p,p2
      parameter(pi=3.14159265359)
      integer:: ng,iargc,exponent
      real:: tau,wt,minval,threshold
      character::fname*40,outputname*40,cmd*140
      character:: input*40,output*40,line*100
      external:: getarg

      ng=1
      bunny=1

      i=iargc()

      if(i.ne.2)then
         write(*,'(a70)')'you do not have the input and output file names on the command line'
         write(*,*) ' enter syntax like   ./exe  inputfile outputfile'
         stop
      endif

      call getarg(1,input)
      call getarg(2,output)

      cmd="wc -l "//trim(input)//" > cgrcmu"
      call system(cmd)

      cmd="awk -v min=10 '{if(NF==5){if($4<min)min=$4}} END {print min} ' "//trim(input)//" >> cgrcmu " 

      call system(cmd)

      open(1,file='cgrcmu')
      read(1,*)num_orients
      read(1,*)minval
      write(*,*)'smallest VF = ',minval
      close(1,status='delete')

      allocate(phi1(num_orients),phi(num_orients),phi2(num_orients),area(num_orients),stat=i)
      if(i.ne.0)stop 'Failure allocating memory'

      exponent=6
      !do while(minval*10**exponent<10)
      !  exponent=exponent+1
      !enddo
      write(*,*)'Using a power exponent=',exponent

      open(2,file=input,status='old')

      do i=1,4
         read(2,'(a)') line !dummy line read
         !write(*,*)line
      enddo

      i=1
      threshold=0.01

 101  read(2,*,end=901)p1,p,p2,wt,tau

      !write(*,*)p1,p,p2,wt,tau,i,threshold

      if(p1>threshold)threshold=p1
      if(p>threshold)threshold=p
      if(p2>threshold)threshold=p2

      phi1(i)=p1
      phi(i)=p
      phi2(i)=p2

      area(i)=int(wt*10.0**exponent)
      
      i=i+1

      goto 101
     
 901  close(2)

!determines whether input is in degrees or radians
! and modifies according to the result
      if(abs(threshold)>10.0)then
        do ng=1,i
          phi1(ng)=phi1(ng)*pi/180.0      
          phi(ng)=phi(ng)*pi/180.0      
          phi2(ng)=phi2(ng)*pi/180.0      
        enddo
      endif

      call angfile2(i,output)    ! outputs texture to ang file
     
      cmd="du -BMB "//trim(output)
      write(*,*)'Outputfile=',trim(output),'will be of filesize='
      call system(cmd)

      deallocate(phi1,phi,phi2,area)

      end program
!________________________________
      subroutine angfile2(arraysize,fname)

! Adds standard comment lines to top of new ANG file

        use master

        implicit none

        character:: tempword*90,fname*40
        integer:: i,j,k,fileio
        integer:: arraysize,grn_area
        real,parameter::pi=3.14159265359

        fileio=1

        open(5,file=fname,iostat=fileio)
        if(fileio==1)then
           write(*,*)'error opening ',fname
           stop 'error opening file'
        endif

 900    format(a)
        write(5,900)'# TEM_PIXperUM          1.000000'
        write(5,900)'# x-star                231'
        write(5,900)'# y-star                244'
        write(5,900)'# z-star                302'
        write(5,900)'# WorkingDistance       17'
        write(5,900)'#'
        write(5,900)'# Phase 1'
        write(5,900)'# MaterialName  	Nickel'
        write(5,900)'# Formula     	'
        write(5,900)'# Info		'
        write(5,900)'# Symmetry              43'
        
        tempword='# LatticeConstants      3.520 3.520 3.520 90.000'
        tempword=  tempword(1:48) // '  90.000  90.000'
        write(5,900)tempword(1:70)
        write(5,900)'# NumberFamilies        4'
        write(5,900)'# hklFamilies   	 1  1  1 1 0.000000'
        write(5,900)'# hklFamilies   	 2  0  0 1 0.000000'
        write(5,900)'# hklFamilies   	 2  2  0 1 0.000000'
        write(5,900)'# hklFamilies   	 3  1  1 1 0.000000'
        write(5,900)'# Categories 1 1 1 1 1 '
        write(5,900)'#'
        write(5,900)'# GRID: SqrGrid#'

! (X,Y) coordinates follow C syntax (0:x-1) (0:y-1)
        i=0
        j=0

        do k=1,arraysize
! Write same orientation X times depending on area.
           do grn_area=1,area(k) 
              write(5,901)phi1(k),Phi(k),phi2(k),(1.0*i),(1.0*j),74.0,0.771,0,-700        
              i=i+1
              if(i > 200)then
                 i=0
                 j=j+1
              endif
           enddo
        enddo

 901    format(f8.3,f9.3,f9.3,f10.3,f10.3,f7.1,f7.3,i3,i8)
 500    close(5)

        return
        end subroutine angfile2
