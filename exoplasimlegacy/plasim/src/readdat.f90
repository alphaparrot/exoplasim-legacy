
       program readfile
       
       real :: filedata(2048,2)
       
       call readdat('specs/trappist1.dat',2,2048,filedata)
       
       do nl = 1,8
          print *,filedata(nl,:)
       enddo

       end

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