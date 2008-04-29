! Author: Chris Roberts
! Modified: 1/31/07
! Converts a traditional MC photo file to a binary format (.raw)
! for Paraview
! To compile:  g95 MCconvert2paraview.f -o mc2pv
! To run program on a discrete digital photo  file:  ./mc2pv file.ph

      program convertMCfile

      implicit none
      integer:: isite,jsite,ksite
      integer:: xx,yy,zz,i,j,k
      character:: fname*40,binfile*40
      integer,allocatable,dimension(:,:,:)::spins
      real::r1,r2,r3
      character::word*10
      integer::Q,iargc
      integer,allocatable,dimension(:)::rngspins
      integer::tempspin,tempspin2
      integer::ir1,ir2,iran

! RANDOMIZE the spin ID orders since 2 neighbors will only 
! differ by 10 units which means the color shade will be similar
! when rendered in Paraview.
! just swap two values into different array indices
      
      if(iargc().ne.1)then
         write(*,*)'Syntax is as follows:  ./mc2pv filename'
         stop 'No filename given as commandline argument'
      endif
      call getarg(1,fname)

      open(17,file=fname,status='old')
      read(17,*) xx,yy,zz
      read(17,*)word,r1,r2,r3,Q 
      read(17,*) ! Line read

      i=10
      allocate(spins(xx,yy,zz),rngspins(0:Q),stat=i)
      if(i/=0)stop 'error allocating spins array'

      do i=0,Q
         rngspins(i)=i
      enddo

      do i=0,Q
         ir1=iran(Q)
         ir2=iran(Q)
         tempspin=rngspins(ir1)
         tempspin2=rngspins(ir2)
         rngspins(ir1)=tempspin2
         rngspins(ir2)=tempspin
      enddo

      do  ksite=1,zz
         do  jsite=1,yy
            read(17,992) (spins(isite,jsite,ksite),isite=1,xx)
 992        format(20i6)
         enddo
      enddo
      close(17)

      i=len_trim(fname)
      binfile=fname(1:i)//'.raw'
      write(*,*)'Name of paraview file=',binfile
      
      open(3,file=binfile,form='unformatted',status='replace')
! ROW MAJOR
!	write(3)(((rngspins(spins(i,j,k)),k=1,zz),j=1,yy),i=1,xx)
! COLUMN MAJOR
      write(3)(((rngspins(spins(i,j,k)),i=1,xx),j=1,yy),k=1,zz)
      close(3)
      
      end
!______________________________________

      integer function iran(Q)

! Obtain an integer on the scale 1<=iran<=Q

      implicit none
      integer::Q
      real::rng
      logical,save::flag=.false.
      integer,dimension(8)::array
      integer::n,k

      if(.not.flag)then
        call random_seed(size=n)
!        call random_seed(get=array(1:n))
!        write(*,*)(array(k),k=1,n)
        call random_seed()  ! Reseed the RNG
!        call random_seed(get=array(1:n))
!        write(*,*)'RNG seeds=',(array(k),k=1,n)
        flag=.true.
      endif

      call random_number(rng)
      iran=int(1.0D0*Q*rng)

      return
      end
