

      program newsnow
      
      integer, parameter :: NLAT = 32
      integer, parameter :: NLON = 64
      integer, parameter :: NUGP = NLAT*NLON
      real, parameter :: small = 1.0e-10
      
      real :: dsnow1(NUGP) = 0.
      real :: dsnow2(NUGP) = 0.
      real :: dsnow3(NUGP) = 0.
      real :: dsnow4(NUGP) = 0.
      real :: dsnow5(NUGP) = 0.
      real :: dsnowz(NUGP) = 0.
      real :: xsnowz(NUGP) = 0.
      real :: lsm(NUGP) = 0.
      real :: persist(NUGP) = 0.
      
      real :: deltyrs = 0.
      real :: delsnow = 0.
      
      character (len=32) :: deltat
      character (len=32) :: sfile0
      character (len=32) :: sfile1
      character (len=32) :: sfile2
      character (len=32) :: sfile3
      character (len=32) :: sfile4
      
      integer nflag
      integer iflag
      
      call getarg(1,sfile0)
      call getarg(2,sfile1)
      call getarg(3,sfile2)
      call getarg(4,sfile3)
      call getarg(5,sfile4)
      call getarg(6,deltat)
      
      read(deltat,*)deltyrs
      
      open(35,file=trim(sfile0),form='unformatted')
      read(35) dsnow1(:)
      close(35)
      
      open(36,file=trim(sfile1),form='unformatted')
      read(36) dsnow2(:)
      close(36)
      
      open(37,file=trim(sfile2),form='unformatted')
      read(37) dsnow3(:)
      close(37)
      
      open(38,file=trim(sfile3),form='unformatted')
      read(38) dsnow4(:)
      close(38)
      
      open(39,file=trim(sfile4),form='unformatted')
      read(39) dsnow5(:)
      close(39)
      
      open(40,file='newdsnow',form='unformatted')
      read(40) dsnowz(:)
      close(40)
      
      open(41,file='newxsnow',form='unformatted')
      read(41) xsnowz(:)
      close(41)
      
      open(42,file='lsm',form='unformatted')
      read(42) lsm(:)
      close(42)
      
      open(42,file='persist',form='unformatted')
      read(42) persist(:)
      close(42)
      
! On land, only change the snowpack if
!      a) there's already snow there, and either
!        b1) the snowpack is decreasing, or
!        b2) the snow has persisted there for an entire year.
     
      do j=1,NUGP
         iflag=1
         if (lsm(j) .lt. 0.5) iflag = 0
         if (dsnowz(j) .lt. small) iflag = 0
         delsnow = deltyrs * (dsnow1(j)+dsnow2(j)+dsnow3(j)+dsnow4(j)+dsnow5(j)) * 0.2
         nflag = 0
         if (delsnow .lt. 0) nflag = nflag + 1
         if (persist(j) .gt. 0.5) nflag = nflag + 1
         nflag = nflag*iflag
         if (nflag .gt. 0)  dsnowz(j) = min(max(dsnowz(j)+delsnow, 0.0),3.0e3)
      enddo
      
      do j=1,NUGP
         if (lsm(j) .le. 0.5) then
           delsnow = deltyrs * (dsnow1(j)+dsnow2(j)+dsnow3(j)+dsnow4(j)+dsnow5(j)) * 0.2
           xsnowz(j) = min(max(xsnowz(j) + delsnow,0.0),10.0)
         endif
      enddo
      
      open(17,file='newdsnow',form='unformatted')
      open(18,file='newxsnow',form='unformatted')
      
      write(17) dsnowz(:)
      write(18) xsnowz(:)
      
      close(17)
      close(18)
      
      stop
      end
      