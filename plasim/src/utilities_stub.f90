
!==============Utilities useful for global stuff========================
!-----------------------------
      subroutine makeareas
      use carbonmod
      
      real :: lat1(NLAT)
      real :: lat2(NLAT)
      
      lat1(1) = 0.5*PI
      lat1(NLAT) = 0.5*(deglat(NLAT-1)+deglat(NLAT))*PI/180.0
      lat2(1) = 0.5*(deglat(1)+deglat(2))*PI/180.0
      lat2(NLAT) = -0.5*PI
      do jlat=2,NLAT-1
         lat1(jlat) = 0.5*(deglat(jlat-1)+deglat(jlat))*PI/180.0
         lat2(jlat) = 0.5*(deglat(jlat+1)+deglat(jlat))*PI/180.0
      enddo
      
      do jlat=1,NLPP
         do jlon=1,NLON !0.5 for 2pi / 4pi, since we want the sum to be 1.
            dglobe((jlat-1)*NLON+jlon) = 0.5/NLON*(sin(lat1(jlat))-sin(lat2(jlat)))
         enddo
      enddo    
            
      end subroutine makeareas

!============ADDITIONAL MPIMOD UTILITIES===========================   
      
!     =================
!     SUBROUTINE MPGADN
!     =================

      subroutine mpgadn(pf,pp,n) ! gather double-precision values

      return
      end subroutine mpgadn
      
!     =================
!     SUBROUTINE MPGADN
!     =================

      subroutine mpgarn(pf,pp,n) ! gather single-precision values

      
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
     
      doutvar = 0.
      
      do i=1,NUGP
         doutvar = doutvar + pf(i)*dglobe(i)

      enddo
      
      
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
!     Write a gridpoint array to a text file

      subroutine finishuptext(dd,fname)
      use pumamod
      
      character (len=*) :: fname
      real :: dd(NHOR)
 
      
      open(93,file=fname,status='unknown')
      do i=1,nn
        write(93) dd(:)
      enddo
      close(93)
      
      return
      end subroutine finishuptext      
      
      
!--------------------------------------------------------------------72
!
!     Write a gridpoint array to an unformatted file

      subroutine finishup(dd,fname)
      use pumamod
      
      character (len=*) :: fname
      real :: dd(NHOR)
      
     

      open(93,file=fname,form='unformatted')
      write(93) dd(:)
      close(93)

      
      return
      end subroutine finishup
      
!--------------------------------------------------------------------72
!
!     Write a spectral array to an unformatted file

      subroutine finishsp(dd,fname)
      use pumamod
      
      character (len=*) :: fname
      real :: dd(NSPP)
      
      
      open(93,file=fname,form='unformatted')
      write(93) dd(:)
      close(93)

      
      return
      end subroutine finishsp
       
!--------------------------------------------------------------------72
!
!  In single-thread mode, write a spectral array to an unformatted file

      subroutine finishfsp(ddn,fname)
      use pumamod
      
      character (len=*) :: fname
      real :: dd(NESP)
      
      open(93,file=fname,form='unformatted')
      write(93) dd(:)
      close(93)
     
      return
      end subroutine finishfsp
           
!--------------------------------------------------------------------72
!
!       Write a latitude array to an unformatted file

      subroutine finishlat(dd,fname)
      use pumamod
      
      character (len=*) :: fname
      real :: dd(NLPP)


      
      open(93,file=fname,form='unformatted')
      write(93) dd(:)
      close(93)

      
      return
      end subroutine finishlat

      
!--------------------------------------------------------------------72
!
!      Read a gridpoint array from an unformatted file

      subroutine readarray(dd,fname)
      use pumamod
      
      character (len=*) :: fname
      real :: dd(NHOR)
      
      open(93,file=fname,form='unformatted')
      read(93) dd(:)
      close(93)
      
      
      return
      end subroutine readarray      
      
! >>>--ayp      
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
      
                               