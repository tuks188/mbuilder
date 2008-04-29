! AUTHORS:  Chris Roberts  12/28/06
!           Anthony Rollett
!
! Debugged on 7/6/07

!OVERVIEW: renumbers spin ID's to be a consecutive list from
!          1 to Q  or 0 to Q-1 depending on user responses.

! OPTIONAL: input the orientation(texture) file to have the
!           (phi1,phi,phi2)--> spin # relationship updated.

! INPUT: discrete structure on a simple cubic grid
! OUTPUT:  structure with updated spin IDs

!_________________________________________________________
module global
  implicit none
  
  
  !Does structure have periodic boundary conditions
  logical,parameter::periodic=.true.
  
  ! # of Neighbors to search
  ! Prefer to search to 1st NN      --> nbors=6
  ! May be modified to 1st&2nd      --> nbors=18
  !                    1st,2nd,&3rd --> nbors=26
  integer, parameter:: nbors=26
  
  !part_id --> particle id
  !lowgid --> lowest grain ID allowed (either 0 or 1)
  integer::part_id,lowgid
  
  ! CONSTANTS
  double precision::PI
  
  ! File Characteristics
  integer:: ngrains,Q

  integer,dimension(:,:,:),allocatable::spins
  integer,dimension(:),allocatable::old2newmap
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
program spinchange 
  
  use global
  implicit none
  character(len=30) :: fname
  integer::iargc,answer
  integer::rate,t1,t2
  integer::currentmin
  
  PI=dacos(-1.0D0)
  
  write(*,*)'Predefined parameters'
  write(*,*)'periodic=',periodic
  write(*,*)'nbors=',nbors

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
  write(*,*)'What is the spin ID for the inert particles?'
  write(*,*)'If no inert particles, enter -10'
  part_id=-10
!  read(*,*)part_id
  write(*,*)'Do you want the first available spin ID to be'
  write(*,*)'0 or 1 ?'
  lowgid=0
!  read(*,*)lowgid
  if(lowgid/=0.and.lowgid/=1)stop 'Lowgid not equal to 1 or 0'

  currentmin=minval(spins,mask= spins>part_id)
  write(*,*)'currentmin=',currentmin
  

  write(*,*)'Was the microstructure generated using CA or MS Builder?'
  write(*,*)'Or are the spins continuous (1,2,3,4,5,..) not (1,2,5,6,10...)'
  write(*,*)'Enter 1 for continuous spins or CABuilder or MSBuilder'
  write(*,*)'Enter 2 otherwise'
!  read(*,*)answer
  answer=2

  if(answer==1)then
     call switcharoo(lowgid,currentmin,part_id)
  else
     call unique_id()
  endif

  
  ! Rewrite Data file with updated spin array
  call system_clock(count=t1)
  call rewrite(fname)
  call system_clock(count=t2,count_rate=rate)
  write(*,*)'rewrite=',(t2-t1)/(1.0*rate) ,'seconds'
  
  call orient_reassign()

  deallocate(spins,old2newmap)

end program spinchange
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
  
  ! Allocate memory and initialize to zero
  good=10
  allocate(spins(mx,my,mz),old2newmap(0:q+50),stat=good)
  if(good.gt.0)stop 'allocation issue'
  
  spins=0
  
  ! READ in spin data
  do ksite=1,mz
     do  jsite=1,my
        read(1,654,end=18)(spins(isite,jsite,ksite),isite=1,mx)
     enddo
  enddo
18 close(1)
654 format(20i6)
  
  return
end subroutine readDATA
!______________________________________________________
subroutine rewrite(fname)
  !
  ! Outputs new spins array in standard MC format
  !      
  use global
  implicit none
  integer::good=0
  integer::isite,jsite,ksite,samax=0
  character::fname*(*),fname2*40
  
  fname2=trim(fname)//".renumber"
  
  open(1,file=fname2,status='replace',iostat=good)
  if(good/=0)stop 'error opening file'
  write(*,*)'Writing data to ',fname2
  write(*,*)'Q value=',q
  
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
655 format(20i6)
989 format(3i8)
990 format("'",a6,"'",f19.2,f7.3,f5.1,i10)
991 format(f6.3,f6.3,f6.3,2(1x,i10))
  
  return
end subroutine rewrite

!_______________________________________________

subroutine switcharoo(futuremin,currentmin,particle)
  
  use global
  implicit none
  integer::futuremin,currentmin,particle
  integer::x,y,z,adjustment
  
  adjustment=futuremin-currentmin
  old2newmap=-100

  write(*,*)'Old max GID =',maxval(spins),q
  write(*,*)'Old min GId =',currentmin
  
  do z=1,mz
     do y=1,my
        do x=1,mx
           if(spins(x,y,z)>particle)then
              old2newmap(spins(x,y,z))=spins(x,y,z)+adjustment
              spins(x,y,z)=spins(x,y,z)+adjustment
           endif
        enddo
     enddo
  enddo
  
  q=maxval(spins)
  write(*,*)'New max GID =',q
  write(*,*)'New min GID =',futuremin,minval(spins,mask=spins>particle)
  
  return
end subroutine switcharoo
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
end subroutine neighbor_config
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
  !  samax is the maximum area, SPINAMAX is its spin
  ! Initialize variables to zero
  write(*,*) 'entering unique_id.f'
  
  ot=1000
  ot2=ot**2
  tail=0
  head=0
  !Lowest spin value given to a grain
  ! (For C/C++, may need to set to 0, so that 1st Grain ID=0)
  ! (For Fortran, set myid=1, so that 1st Grain ID=1 )
  myid=lowgid
  old2newmap=-100

  do 100 i=1,mx
     do 90 j=1,my
        do 80 k=1,mz
           nucl=ot2*i+ot*j+k
           
           
           !      if the nucleus has been burned previously, or is a pore
           !        or second phase, go to next nucleus
           !     Uses spin=0 as a spin ID
           
           if(spins(i,j,k)<=part_id) go to 80
           
           !        find the spin of the new nucleus, then burn it
           !        and set the area of the new cluster equal to one
           !
           spin=spins(i,j,k)
           old2newmap(spin)=myid
           spins(i,j,k)=-1*myid-ot
           
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
           
50         site=list(tail)
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
                 spins(in,jn,kn)=-1*myid-ot
              endif
60         enddo
           
           !If the tail of the list does not equal the head of the list, then:
           !sites were added to the cluster in the last pass, so all the sites
           !in the cluster have not been checked for neighbors.  So, check the
           !unchecked site closest to the tail of the list.
           
           if(tail.ne.head) then
              tail=tail+1
              go to 50
           else
              myid=myid+1 ! Increment spin ID
           endif
           
           !If the tail = the head, the cluster has been completely burned.
           !Record the area, and go to the next cluster.
           
80      enddo
90   enddo
100 enddo
  
  ! Renumbering Spin Array
  
  do k=1,mz
     do j=1,my
        do i=1,mx
           if(spins(i,j,k)<part_id)spins(i,j,k)=abs(spins(i,j,k)+ot)
        enddo
     enddo
  enddo

!  do k=1,mz
!     do j=1,my
!        do i=1,mx
!           if(spins(i,j,k)==1)write(*,*)spins(i,j,k),i,j,k
!        enddo
!     enddo
!  enddo

  i=minval(spins,mask=spins>part_id)
  j=maxval(spins)
  write(*,*)i,j
  myid=j-i+1
  write(*,*)'# unique grains =',myid
  q=j
  write(*,*)'while the max gid is Q=',q

  return
end subroutine unique_id
!__________________________________________

subroutine orient_reassign()
  
  ! Orientation reasignment must occur since it each grain
  ! may not be simply bumped up or down in spin ID's.
  ! (i.e. domain with a non-consecutive list of spin numbers)
  ! Expects a simple list file with GID and corresponding 
  ! Euler angles on the same line
  !   GID   phi1  phi   phi2
  !    1   2.34  1.23  0.79
  
  use global
  implicit none
  integer::header,r2d
  character::fname*40,fname2*40,line*80
  integer::i,io,gid
  double precision::phi1,phi,phi2
  type angles
     double precision::phi1,phi,phi2
  end type angles
  type(angles),dimension(:),allocatable::orient
  
  if(lowgid==1)then
     allocate(orient(1:Q))
  elseif(lowgid==0)then
     allocate(orient(0:Q))
  else
     write(*,*)'lowgid does not equal 1 or 0'
     stop
  endif
  
  write(*,*)'What is the name of the orientation file?'
  write(*,*)'enter none if no file exist'
  !read(*,*)fname
  fname='none'
  if(trim(fname)=="none")return
  
  fname2=trim(fname)//".renumber"
  
  write(*,*)'how many header lines?'
  read(*,*)header
  
  open(1,file=fname,status='old',iostat=io)
  open(2,file=fname2,status='replace')
  
  if(header>0)then
     do i=1,header
        read(1,*)line
        write(2,*)line
     enddo
  endif
  
  do while(io==0)

     read(1,*,iostat=io,end=605)gid,phi1,phi,phi2
     orient(old2newmap(gid))%phi1=phi1
     orient(old2newmap(gid))%phi=phi
     orient(old2newmap(gid))%phi2=phi2

  enddo
  
605 close(1)
  
  if(lowgid==1)then
     do i=1,Q
        write(2,702)i,orient(i)
     enddo
  elseif(lowgid==0)then
     do i=0,Q
        write(2,702)i,orient(i)
     enddo
  endif
  
702 format(i10,1x,3(f8.4,1x))
  close(2)
  deallocate(orient)
  
  return
end subroutine orient_reassign
