!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     This is a program to                                                   !
!     1. generate micro.input format file from the file 'output'             !
!     2. find the largest regionID number of J. Fridy's 'findPoints.c"       !
!        program.                                                            !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 
	Program microFormat
	implicit none

	integer i,j,k,m,x,y,z,ii,maxID
	integer regionID2(125000000)

        character fname*20,keyfile*6
        character(10)::str

c       The A,B,C must be the same as that in "voxelpoints.for" This
c       program needs the total number of voxels

c       Obtaining filename from command line
	if(iargc()/=4)stop 'USAGE: ./microFormat X Y Z input_filename'
        call getarg(1,str)
        read(str,'(i4)')x
        call getarg(2,str)
        read(str,'(i4)')y
	call getarg(3,str)
        read(str,'(i4)')z
	call getarg(4,str)
	read(str,*)fname

c        write(*,*) 'What is the input filename (from findPoints)?'
c        read(*,'(a)') fname
c        write(*,*) 'What keyword (6 chars) to use?'
c        read(*,'(a)') keyfile
	open(23, file='ellipsoid_Growth.ph', status='unknown')
	m=x*y*z

	
	write(*,*) 'generating "ellipsoid_Growth.ph"..'
	
	open(21, file=fname, status='old')
	maxID=0
	do ii=1, m
		read(21,*) regionID2(ii)
		if (regionID2(ii).gt.maxID) then
		maxID=regionID2(ii)
		endif
	enddo
        close(21)

	write(23,989) x,y,z
	write(23,*) "'grwXXX'             170.00  1.000  1.0", maxID
	write(23,*) '1.000 0.000 0.000'	
	write(23,1000) (regionID2(ii),ii=1,m)

	write(*,*) maxID
 989	format(3i8)
1000  format(20i6)
	end
