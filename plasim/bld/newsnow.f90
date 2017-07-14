

      program newsnow
      
      integer, parameter :: NLAT = 32
      integer, parameter :: NLON = 64
      integer, parameter :: NUGP = NLAT*NLON
      
      real :: dsnow1(NUGP) = 0.
      real :: dsnow2(NUGP) = 0.
      real :: dsnow3(NUGP) = 0.
      real :: dsnowz(NUGP) = 0.
      real :: xsnowz(NUGP) = 0.
      real :: lsm(NUGP) = 0.
      
      real :: deltyrs = 0.
      
      character (len=32) :: deltat
      character (len=32) :: sfile0
      character (len=32) :: sfile1
      character (len=32) :: sfile2
      
      call getarg(1,sfile0)
      call getarg(2,sfile1)
      call getarg(3,sfile2)
      call getarg(4,deltat)
      
      read(deltat,*)deltyrs
      
      open(35,file=trim(sfile0),form='unformatted')
      read(35) dsnow1(:)
      close(35)
      
      open(36,file=trim(sfile1),form='unformatted')
      read(36) dsnow2(:)
      close(36)
      
      open(37,file=trim(sfile2),form='unformatted')
      read(37) dsnow3(:)
      
      open(38,file='newdsnow',form='unformatted')
      read(38) dsnowz(:)
      close(38)
      
      open(39,file='newxsnow',form='unformatted')
      read(39) xsnowz(:)
      close(39)
      
      open(40,file='lsm',form='unformatted')
      read(40) lsm(:)
      close(40)
      
      
      where (lsm .gt. 0.5)
         dsnowz(:) = min(max(dsnowz(:) + deltyrs*(dsnow1(:)+dsnow2(:)+dsnow3(:))/3.0,0.0),3.0e3)
      endwhere
      
      where (lsm .le. 0.5)
         xsnowz(:) = min(max(xsnowz(:) + deltyrs*(dsnow1(:)+dsnow2(:)+dsnow3(:))/3.0,0.0),10.0)
      endwhere
      
      open(17,file='newdsnow',form='unformatted')
      open(18,file='newxsnow',form='unformatted')
      
      write(17) dsnowz(:)
      write(18) xsnowz(:)
      
      close(17)
      close(18)
      
      stop
      end
      