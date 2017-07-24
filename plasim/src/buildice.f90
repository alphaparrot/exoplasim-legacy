

      program buildice
      
      integer, parameter :: NLAT = 32
      integer, parameter :: NLON = 64
      integer, parameter :: NUGP = NLAT*NLON

      real :: dsnowz(NUGP) = 0.
      real :: lsm(NUGP) = 0.
      
      open(40,file='newdsnow',form='unformatted')
      read(40) dsnowz(:)
      close(40)
      
      open(42,file='lsm',form='unformatted')
      read(42) lsm(:)
      close(42)
      
      where (lsm .gt. 0.5)
         dsnowz(:) = min(max(dsnowz(:) + 400.0,0.0),3.0e3)
      endwhere
      
      open(17,file='newdsnow',form='unformatted')  
      write(17) dsnowz(:)
      close(17)
      
      stop
      end
      