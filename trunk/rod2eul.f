	program rod2eul
	
c       converts Rodrigues vector to Bunge Euler angles 
c       d1== phi1,  d2== Phi,  d3== phi2
c
c       input files : orts.txt (obtained from D. Saylor's annealfinal.c)
c                     evodf.txt (contains Vf of each orientation)
c       output file : fridyTexin (Texin1/2 format for rex3d)

        implicit none
	integer:: i,j,id,index,io=10,io2=10
        integer::ngrains
        double precision,allocatable,dimension(:):: vol
	real:: d1,d2,d3,rodr(3)
        real:: sum,diff,phi1,Phi,phi2	
        real,parameter:: pi=3.1415926535
        character::cmd*270
	open(21,file='orts.txt', iostat=io,status='old')

        cmd="awk ' { if( $1 ~ /<id>/ ) x=$2 } { if( $1 ~ 
     &  /<nregions>/ ) xmax=$2 } { if( $1 ~ /<volume>/ ) arr[x]=$2 } 
     &  END { print xmax ; for (i=0; i<xmax ; i++) print i,arr[i] } ' 
     &   cellIdealization.xml > tempfreq.txt "
 !       write(*,*)cmd
        call system(cmd)

	open(25,file='tempfreq.txt',iostat=io2,status='old')
	open(22,file='fridyTexin', status='unknown')
        if(io/=0)stop 'orts.txt file does not exist'
        if(io2/=0)stop 'tempfreq.txt file does not exist'
	
        read(25,*)ngrains
	allocate(vol(0:ngrains))

	write(22,*) 'tex.comp.out'  
	write(22,*) 'Evm    F11    F12    F13    F21    F22    F23    F31
     &    F32    F33    nstate'
	
	write(22,*) '0.000  1.000  0.000  0.000  0.000  1.000  0.000  0.00
     &0  0.000  1.000  1'
      
        write(22,601) 
601	format('Bunge:phi1   PHI   phi2  ,gr.wt.,  tau,    taus;taumo
     &des/tau;       XYZ= 1 2 3')

	do while(io==0)
          read(21,*,iostat=io,end=676)index,rodr(1),rodr(2),rodr(3)
                read(25,*)id,vol(id)
		sum=atan(rodr(3))
		diff=atan(rodr(2)/rodr(1))
		d1=sum+diff
		d2=2.*atan(rodr(1)*cos(sum)/cos(diff))
		d3=sum-diff
	
		phi1=(180/pi)*d1
		Phi=(180/pi)*d2
		phi2=(180/pi)*d3
		write(22,1000) phi1,Phi,phi2,1.00,vol(id)
	enddo
676     close(21)
        close(22)
        close(25,status='delete')

        deallocate(vol)
1000  format(2x,3(g10.4),2g8.3)
        write(*,*) "fridyTexin should be moved to texin1 and texin2."

        end
