
      program sumitup

! Tests whether evodf.txt and evmdf.txt sum up to 1.0
! since each cell represents a volume fraction of each component
!

      implicit none
      integer::cell,iargc,io
      double precision::sum,frac
      character(len=40)::fname

      sum=0.0D0
      if(iargc().ne.1)stop 'no filename on command line'
      call getarg(1,fname)

      open(1,file=fname,iostat=io)

      do while(io==0)
       read(1,*,iostat=io,end=999)cell,frac
       sum=sum+frac
      enddo
     
999      close(1)
      write(*,*)'sum=',sum

      end
     
