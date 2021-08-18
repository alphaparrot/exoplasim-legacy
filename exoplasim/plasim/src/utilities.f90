
!==============Utilities useful for global stuff========================
!-----------------------------
      subroutine makeareas
      use carbonmod
      
      real :: zdeglat(NLAT)
      real :: lat1(NLAT)
      real :: lat2(NLAT)
      
      call mpgarn(zdeglat,deglat,NLPP)
      
      if (mypid == NROOT) then
      
      lat1(1) = 0.5*PI
      lat1(NLAT) = 0.5*(zdeglat(NLAT-1)+zdeglat(NLAT))*PI/180.0
      lat2(1) = 0.5*(zdeglat(1)+zdeglat(2))*PI/180.0
      lat2(NLAT) = -0.5*PI
      do jlat=2,NLAT-1
         lat1(jlat) = 0.5*(zdeglat(jlat-1)+zdeglat(jlat))*PI/180.0
         lat2(jlat) = 0.5*(zdeglat(jlat+1)+zdeglat(jlat))*PI/180.0
      enddo
      
      endif
      
      call mpscrn(lat1,NLPP)
      call mpscrn(lat2,NLPP)
      
      do jlat=1,NLPP
         do jlon=1,NLON !0.5 for 2pi / 4pi, since we want the sum to be 1.
            dglobe((jlat-1)*NLON+jlon) = 0.5/NLON*(sin(lat1(jlat))-sin(lat2(jlat)))
         enddo
      enddo    
            
      return
      end subroutine makeareas

!============ADDITIONAL MPIMOD UTILITIES===========================   
      
!     =================
!     SUBROUTINE MPGADN
!     =================

      subroutine mpgadn(pf,pp,n) ! gather double-precision values
      use mpimod
      
      real (kind = 8) :: pf(*)
      real (kind = 8) :: pp(*)
      
      call mpi_gather(pp,n,MPI_REAL8,&
     &                pf,n,MPI_REAL8,NROOT,myworld,mpinfo)
      
      return
      end subroutine mpgadn
      
!     =================
!     SUBROUTINE MPGADN
!     =================

      subroutine mpgarn(pf,pp,n) ! gather single-precision values
      use mpimod
      
      real :: pf(*)
      real :: pp(*)
      
      call mpi_gather(pp,n,mpi_rtype,&
     &                pf,n,mpi_rtype,NROOT,myworld,mpinfo)
      
      return
      end subroutine mpgarn  
!      
!============================Output and Summation utilities==========72   
!
!        Compute a weighted global average using cell areas
!
      subroutine write_short(pf,doutvar)
      use carbonmod
      
      real :: pf(NHOR)
      real :: zf(NUGP)
      real :: zzf(NUGP)
      real :: da(NUGP)
      

      doutvar = 0.
      
      call mpgagp(zf,pf,1)
      call mpgagp(da,dglobe,1)
      
      zzf(:) = zf(:)
      
      if (mypid == NROOT) then
      
      do i=1,NUGP
         doutvar = doutvar + zzf(i)*da(i)

      enddo
      
      endif
      
      return
      end subroutine write_short
      
!--------------------------------------------------------------------72
!              Data output to a text file

      subroutine writegtextarray(dd,nn,fname)
      
      real :: dd(nn)
      
      integer nn,i
      character fname*13
      
      
      open(unit=23,file=trim(adjustl(fname)),status='unknown')
      do i=1,nn
         write(23,*) dd(i)
      enddo
      close(23)
      
      
      return
      end subroutine writegtextarray     
      
     
!--------------------------------------------------------------------72
!
!     Write a gridpoint array to an unformatted file

      subroutine finishup(dd,fname)
      use pumamod
      
      character (len=*) :: fname
      real :: dd(NHOR)
      real :: ddn(NUGP)
      
      
      call mpgagp(ddn,dd,1)
      
      if (mypid==NROOT) then
         open(93,file=fname,form='unformatted')
         write(93) ddn(:)
         close(93)
      endif
      
      return
      end subroutine finishup

     
!--------------------------------------------------------------------72
!
!     Write a gridpoint array to a text file

      subroutine finishuptext(dd,fname)
      use pumamod
      
      character (len=*) :: fname
      real :: dd(NHOR)
      real :: ddn(NUGP)
      
      
      call mpgagp(ddn,dd,1)
      
      if (mypid==NROOT) then
         open(93,file=fname,status='unknown')
         do i=1,nn
           write(93) ddn(:)
         enddo
         close(93)
      endif
      
      return
      end subroutine finishuptext      
      
!--------------------------------------------------------------------72
!
!     Write a spectral array to an unformatted file

      subroutine finishsp(dd,fname)
      use pumamod
      
      character (len=*) :: fname
      real :: dd(NSPP)
      real :: ddn(NESP)
      
      
      call mpgasp(ddn,dd,1)
      
      if (mypid==NROOT) then
         open(93,file=fname,form='unformatted')
         write(93) ddn(:)
         close(93)
      endif
      
      return
      end subroutine finishsp
       
!--------------------------------------------------------------------72
!
!  In single-thread mode, write a spectral array to an unformatted file

      subroutine finishfsp(ddn,fname)
      use pumamod
      
      character (len=*) :: fname
      real ddn(NESP)
      
      if (mypid==NROOT) then
         open(93,file=fname,form='unformatted')
         write(93) ddn(:)
         close(93)
      endif
      
      return
      end subroutine finishfsp
           
!--------------------------------------------------------------------72
!
!       Write a latitude array to an unformatted file

      subroutine finishlat(dd,fname)
      use pumamod
      
      character (len=*) :: fname
      real :: dd(NLPP)
      real :: ddn(NLAT)
      
      call mpgarn(ddn,dd,NLPP)
      
      if (mypid==NROOT) then
         open(93,file=fname,form='unformatted')
         write(93) ddn(:)
         close(93)
      endif
      
      return
      end subroutine finishlat

      
!--------------------------------------------------------------------72
!
!      Read a gridpoint array from an unformatted file

      subroutine readarray(dd,fname)
      use pumamod
      
      character (len=*) :: fname
      real :: dd(NHOR)
      real :: ddn(NUGP)
      
      if (mypid==NROOT) then
         open(93,file=fname,form='unformatted')
         read(93) ddn(:)
         close(93)
      endif
      
      call mpscgp(ddn,dd,1)
      
      return
      end subroutine readarray      

!--------------------------------------------------------------------72
!
! Read a 2D array from a formatted text file with single-line header
      
      subroutine readdat(filename,ndim,nitems,kdata)
      
      character (len=*), intent(in) :: filename
      character (len=128) :: header
      integer, intent(in) :: ndim
      integer, intent(in) :: nitems
      real, intent(out) :: kdata(nitems,ndim)
      
      
      open(87,file=filename,form='formatted')
      read(87,'(a)') header
      
      do nl = 1,nitems
         read(87,*) kdata(nl,:)
      enddo
      
      close(87)
      
      
      return
      end subroutine readdat      
      
! >>>--ayp      
            
      
                               