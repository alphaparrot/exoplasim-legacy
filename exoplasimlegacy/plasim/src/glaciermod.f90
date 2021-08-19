!	This is a rudimentary glacier model for capturing the effects of large ice sheets
!	on atmospheric circulation and thermodynamics. It's not intended to be used as a 
!	real-time coupling, merely as a coupling between an external gridpoint glacier
!	model that may run on timescales of years to decades to centuries, and PlaSim, running
!	on timescales of minutes to hours. The external model sets the glacier/snowpack height,
!	and this translates that into a change in orography through the surface geopotential
!	height. This does this by keeping track of two orography fields; the lithographic orography
!	and the glacier orography. The surface orography is the sum of the two. Every timestep, the
!	model checks to see if existing snow on a gridpoint has melted away (default reaching a
!       depth below 2.0 meters liquid water equivalent. If a gridpoint retains some snow cover 
!	continuously for an entire year, its glacier flag is turned on for the subsequent year.
!	If its snow depth in subsequent years falls below 2 meters, the glacier flag is turned off.
!
!	Due to the ground orography implementation being introduced, this module also introduces
!	the possibility of coupling to an orogeny model.

      module glaciermod
      use landmod
!
!     version identifier (date)
!
      character(len=80) :: gversion = '07.07.2017 by Adiv'
!
!     Parameter   
!

!
!     namelist parameters
!
      character (256) :: glacier_namelist = "glacier_namelist"
      
      integer :: nglacier = 1 ! 1 = implement glaciation, 0 = ignore
      real :: glacelim = 2.0 ! Minimum snow depth in meters liquid water equivalent that has to be maintained
                                 ! year-round to convert the gridpoint to a glacier.
      real :: icesheeth = -1.0 !Initial snow depth
!
!     global arrays
!
      real :: groundoro(NHOR)  = 0.0
      real :: glacieroro(NHOR) = 0.0
      real :: netoro(NHOR) = 0.0
      logical :: persistflag(NHOR) = .TRUE.
      
!
!     global scalars
!
      real :: rhoglac  = 850.    ! glacial ice density (kg/m**3)

!
      end module glaciermod
          
!      
!     =============================
!

      subroutine glacierprep
      use glaciermod

      namelist/glacier_nl/nglacier,glacelim,icesheeth
      
      if (mypid==NROOT) then
         open(23,file=glacier_namelist)
         read(23,glacier_nl)
         close(23)
         write(nud,'(/," *********************************************")')
         write(nud,'(" * GLACIERMOD ",a34," *")') trim(gversion)
         write(nud,'(" *********************************************")')
         write(nud,'(" * Namelist GLACIER_NL from <glacier_namelist> *")')
         write(nud,'(" *********************************************")')
         write(nud,glacier_nl)
      endif
      
      call mpbci(nglacier)
      call mpbcr(glacelim)
      call mpbcr(icesheeth)
      
      if (mypid==NROOT) nglspec = nglacier
      
      call mpbci(nglspec)
      
      end subroutine glacierprep
      
!      
!     =============================
!

      subroutine glacierini(noromax)
      use glaciermod
      
      logical :: ldsnow
      real :: zoro(NUGP) = 0.0
      real :: foro(NHOR) = 0.0
      real :: zsnow(NUGP) = 0.0
      
      integer noromax
      
      if (nrestart > 0.) then
        call mpgetgp('dglac',dglac,NHOR,1)
        call mpgetgp('groundsg',groundoro ,NHOR,1)
        call mpgetgp('dglacsg' ,glacieroro,NHOR,1)
        call mpgetgp('doro'    ,doro      ,NHOR,1)
      else
        call mpsurfgp('doro'    ,doro    ,NHOR,1)
        groundoro(:) = doro(:)
        glacieroro(:) = 0.0
        if (icesheeth .ge. 0.0) then
          dglac(:) = 1.0
          dsnowz(:) = icesheeth
          do jhor=1,NHOR
            if (groundoro(jhor) == 0.0) then
              dsnowz(jhor) = 0.0
              dglac(jhor) = 0.0
            endif
          enddo
        endif
          
        call mpgagp(zsnow,dsnowz,1)
        
        if (mypid==NROOT) then
          write(nud,'(/,"Initial Icesheet Heights")')
          write(nud,'("Maximum: ",f10.2," [m]")') (maxval(zsnow))
          write(nud,'("Minimum: ",f10.2," [m]")') (minval(zsnow))
          write(nud,'("Mean:    ",f10.2," [m]")') (sum(zsnow) / (NUGP))
        endif 
      endif
      
      if (nglacier .eq. 1) then
      
        if (nrestart > 0.) then     
          call mpgagp(zoro,doro,1)
        
          if (mypid==NROOT) then
            write(nud,'(/,"Topography before glaciers")')
            write(nud,'("Maximum: ",f10.2," [m]")') (maxval(zoro) / ga)
            write(nud,'("Minimum: ",f10.2," [m]")') (minval(zoro) / ga)
            write(nud,'("Mean:    ",f10.2," [m]")') (sum(zoro) / (ga*NUGP))
          endif        
        endif
        
        if (mypid==NROOT) inquire(file='restart_dsnow',exist=ldsnow)
        call mpbcl(ldsnow)
        
        if (ldsnow .and. (nrestart > 0.)) then
        
           call readarray(dsnowz,'restart_dsnow')
           
        endif   
        
        call oroini       
           
!          Compute spectral orography
        so(:) = 0.
        doro(:) = doro(:) * oroscale  ! Scale orography
        foro(:) = doro(:)
        
        
        call gp2fc(foro,NLON,NLPP)
        call fc2sp(foro,so)
        call mpsum(so,1)
        
        if (npro == 1) then ! print only in single core runs
           call sp2fc(so,doro)
           call fc2gp(doro,nlon,nlpp)
           write(nud,'(/,"Topography after spectral fitting")')
           write(nud,'("Maximum: ",f10.2," [m]")') maxval(doro) / ga
           write(nud,'("Minimum: ",f10.2," [m]")') minval(doro) / ga
           write(nud,'("Mean:    ",f10.2," [m]")') sum(doro) / (ga * NUGP)
        endif
        
        if (mypid == NROOT) then
           so(:) = so(:) / (cv*cv)
           if (noromax < NTRU) then
            jr=-1
            do jm=0,NTRU
             do jn=jm,NTRU
              jr=jr+2
              ji=jr+1
              if(jn > noromax) then
               so(jr)=0.
               so(ji)=0.
              endif
             enddo
            enddo
           endif ! (noromax < NTRU)
        
!          Initialize surface pressure
        
          if (nspinit > 0) then
             sp(:) = -so(:)*cv*cv / (gascon * tgr)
          endif
        endif ! (mypid == NROOT)
        call mpscsp(sp,spm,1)
        
        where (dsnowz(:) > 30.0) dglac = 1.0 !If we have more than 30 m of lq H2O equivalent in snow/ice, it's a glacier
                                             !30 meters of ice is the minimum thickness for an ice sheet to flow
        where (dls(:) < 0.5) persistflag = .FALSE.
        
        call mpgagp(zoro,doro,1)
        
        if (mypid==NROOT) then
          write(nud,'(/,"New Topography after glaciers")')
          write(nud,'("Maximum: ",f10.2," [m]")') (maxval(zoro) / ga)
          write(nud,'("Minimum: ",f10.2," [m]")') (minval(zoro) / ga)
          write(nud,'("Mean:    ",f10.2," [m]")') (sum(zoro) / (ga*NUGP))
        endif
      
      else    !nglacier == 0
       
        if (nrestart > 0.) then     
          call mpgagp(zoro,doro,1)
        
          if (mypid==NROOT) then
            write(nud,'(/,"Topography before smoothing")')
            write(nud,'("Maximum: ",f10.2," [m]")') (maxval(zoro) / ga)
            write(nud,'("Minimum: ",f10.2," [m]")') (minval(zoro) / ga)
            write(nud,'("Mean:    ",f10.2," [m]")') (sum(zoro) / (ga*NUGP))
          endif    
        endif
          
        call oroini 
                
!        Compute spectral orography
        so(:) = 0.
        doro(:) = doro(:) * oroscale  ! Scale orography
        foro(:) = doro(:)
        
        
        call gp2fc(foro,NLON,NLPP)
        call fc2sp(foro,so)
        call mpsum(so,1)
        
        if (npro == 1) then ! print only in single core runs
           call sp2fc(so,doro)
           call fc2gp(doro,nlon,nlpp)
           write(nud,'(/,"Topography after spectral fitting")')
           write(nud,'("Maximum: ",f10.2," [m]")') maxval(doro) / ga
           write(nud,'("Minimum: ",f10.2," [m]")') minval(doro) / ga
           write(nud,'("Mean:    ",f10.2," [m]")') sum(doro) / (ga * NUGP)
        endif
        
        if (mypid == NROOT) then
           so(:) = so(:) / (cv*cv)
           if (noromax < NTRU) then
            jr=-1
            do jm=0,NTRU
             do jn=jm,NTRU
              jr=jr+2
              ji=jr+1
              if(jn > noromax) then
               so(jr)=0.
               so(ji)=0.
              endif
             enddo
            enddo
           endif ! (noromax < NTRU)
        
!          Initialize surface pressure
        
          if (nspinit > 0) then
             sp(:) = -so(:)*cv*cv / (gascon * tgr)
          endif
        endif ! (mypid == NROOT)
        call mpscsp(sp,spm,1)
        
        call mpgagp(zoro,doro,1)
        
        if (mypid==NROOT) then
          write(nud,'(/,"New Topography after glaciers")')
          write(nud,'("Maximum: ",f10.2," [m]")') (maxval(zoro) / ga)
          write(nud,'("Minimum: ",f10.2," [m]")') (minval(zoro) / ga)
          write(nud,'("Mean:    ",f10.2," [m]")') (sum(zoro) / (ga*NUGP))
        endif
      
      endif !nglacier switch
      
      end subroutine glacierini

     
!     ==================
!     SUBROUTINE oroini
!     ==================

      subroutine oroini
      use glaciermod
      
      parameter(zcvel=4.2)
      parameter(zcexp=0.18)
      parameter(gconst=6.67408e-11)
      parameter(p_mass=5.9722e24)
      parameter(rad_p1=6378137.0)
      parameter(rad_p2=6356752.3)
!
      real zuroff(NLON,NLAT)
      real zvroff(NLON,NLAT)
      real zoro(NLON,NLAT)
      real zoron(NLON,NLAT)
      real zlsm(NLON,NLAT)
      real zsi(NLON,NLAT)
      real zsir(NLON,NLPP)
      
      real thing
      real dhsnow
      real radius
      
!
      ilat = NLAT ! using ilat suppresses compiler warnings for T1
!
      do jlat=1,NLPP
       do jlon=1,NLON
        jhor=(jlat-1)*NLON+jlon
        darea(jhor)=gwd(jlat)
        zsir(jlon,jlat)=sid(jlat)
       enddo
      enddo 
      
      call mpgagp(zoro,doro,1)

      if(mypid==NROOT) then

       write(nud,'(/,"Topography before glaciers and before smoothing")')
       write(nud,'("Maximum: ",f10.2," [m]")') (maxval(zoro) / ga)
       write(nud,'("Minimum: ",f10.2," [m]")') (minval(zoro) / ga)
       write(nud,'("Mean:    ",f10.2," [m]")') (sum(zoro) / (ga*NUGP))
      
      endif
      
      if (nglacier .gt. 0.5) then
      
      groundoro(:) = groundoro(:)*oroscale
      
      do i=1,NHOR ! Add elevation of snowpack/ice sheet
         jlat = i/NLON + 1
!          radius = sqrt(((rad_p1**2*cola(jlat))**2+(rad_p2**2*sid(jlat))**2)/&
!      &                 ((rad_p1*cola(jlat))**2+(rad_p2*sid(jlat))**2))
         radius=6.371e6
         dz = radius**2*groundoro(i)/(gconst*p_mass-radius*groundoro(i))
         dhsnow = dsnowz(i)/(rhoglac/1000.0) !Convert from meters of lq H20 equivalent to actual snow depth
         grav = gconst*p_mass/(radius+dz)**2
         glacieroro(i) = grav*dhsnow - grav/(dz+radius)*dhsnow**2
         doro(i) = groundoro(i) + glacieroro(i)
!          if(mypid==NROOT) then
!          write(nud,'("Gr,Sn: ",f10.2," ",f10.2," ",f10.2," [m]")') (groundoro(i) / ga,glacieroro(i)/ga)
!          write(nud,'("Net: ",f10.2," [m]")') (doro(i)/ga)
!          endif 

      enddo
      
      call mpgagp(zoro,glacieroro,1)
      
      if(mypid==NROOT) then

       write(nud,'(/,"Glacial Topography")')
       write(nud,'("Maximum: ",f10.2," [m]")') (maxval(zoro) / ga)
       write(nud,'("Minimum: ",f10.2," [m]")') (minval(zoro) / ga)
       write(nud,'("Mean:    ",f10.2," [m]")') (sum(zoro) / (ga*NUGP))
      
      endif
      
      call mpgagp(zoro,groundoro,1)
      if(mypid==NROOT) then

       write(nud,'(/,"Ground Topography")')
       write(nud,'("Maximum: ",f10.2," [m]")') (maxval(zoro) / ga)
       write(nud,'("Minimum: ",f10.2," [m]")') (minval(zoro) / ga)
       write(nud,'("Mean:    ",f10.2," [m]")') (sum(zoro) / (ga*NUGP))
      
      endif
      call mpgagp(zoro,doro-groundoro,1)
      if(mypid==NROOT) then

       write(nud,'(/,"Net-Ground Topography")')
       write(nud,'("Maximum: ",f10.2," [m]")') (maxval(zoro) / ga)
       write(nud,'("Minimum: ",f10.2," [m]")') (minval(zoro) / ga)
       write(nud,'("Mean:    ",f10.2," [m]")') (sum(zoro) / (ga*NUGP))
      
      endif
      call mpgagp(zoro,doro-glacieroro,1)
      if(mypid==NROOT) then

       write(nud,'(/,"Net-Glacier Topography")')
       write(nud,'("Maximum: ",f10.2," [m]")') (maxval(zoro) / ga)
       write(nud,'("Minimum: ",f10.2," [m]")') (minval(zoro) / ga)
       write(nud,'("Mean:    ",f10.2," [m]")') (sum(zoro) / (ga*NUGP))
      
      endif
      call mpgagp(zoro,glacieroro+groundoro,1)
      if(mypid==NROOT) then

       write(nud,'(/,"Ground + Glacier Topography")')
       write(nud,'("Maximum: ",f10.2," [m]")') (maxval(zoro) / ga)
       write(nud,'("Minimum: ",f10.2," [m]")') (minval(zoro) / ga)
       write(nud,'("Mean:    ",f10.2," [m]")') (sum(zoro) / (ga*NUGP))
      
      endif
      
      endif
      
      call mpgagp(zsi,zsir,1)
      call mpgagp(zoro,doro,1)
      call mpgagp(zlsm,dls,1)

      if(mypid==NROOT) then

       write(nud,'(/,"Topography after glaciers and before smoothing")')
       write(nud,'("Maximum: ",f10.2," [m]")') (maxval(zoro) / ga)
       write(nud,'("Minimum: ",f10.2," [m]")') (minval(zoro) / ga)
       write(nud,'("Mean:    ",f10.2," [m]")') (sum(zoro) / (ga*NUGP))
      
       zoro(:,:)=MAX(zoro(:,:),0.)
       where(zlsm(:,:) < 1.) zoro(:,:)=zlsm(:,:)-1.

!
!     iterate to remove local minima, but only if we're not already flat.
!
       if (ndesert < 0.5 .and. naqua < 0.5 .and. mars < 0.5 .and. ((maxval(zoro)-minval(zoro))/ga>0.01)) then

 1000   continue
        jconv=0
        do jlat=2,ilat-1
         do jlon=2,NLON-1
          if(zlsm(jlon,jlat) > 0.                                       &
     &      .and. zoro(jlon,jlat) <= zoro(jlon,jlat-1)                  &
     &      .and. zoro(jlon,jlat) <= zoro(jlon,jlat+1)                  &
     &      .and. zoro(jlon,jlat) <= zoro(jlon+1,jlat)                  &
     &      .and. zoro(jlon,jlat) <= zoro(jlon-1,jlat)) then
           zoron(jlon,jlat)=1.+MIN(zoro(jlon+1,jlat),zoro(jlon-1,jlat)  &
     &                            ,zoro(jlon,jlat+1),zoro(jlon,jlat-1))
           jconv=jconv+1
          else
           zoron(jlon,jlat)=zoro(jlon,jlat)
          endif
         enddo
         if(zlsm(1,jlat) > 0.                                           &
     &     .and. zoro(1,jlat) <= zoro(1,jlat-1)                         &
     &     .and. zoro(1,jlat) <= zoro(1,jlat+1)                         &
     &     .and. zoro(1,jlat) <= zoro(2,jlat)                           &
     &     .and. zoro(1,jlat) <= zoro(NLON,jlat)) then
          zoron(1,jlat)=1.+MIN(zoro(2,jlat),zoro(NLON,jlat)             &
     &                        ,zoro(1,jlat+1),zoro(1,jlat-1))
          jconv=jconv+1
         else
          zoron(1,jlat)=zoro(1,jlat)
         endif
         if(zlsm(NLON,jlat) > 0.                                        &
     &     .and. zoro(NLON,jlat) <= zoro(NLON,jlat-1)                   &
     &     .and. zoro(NLON,jlat) <= zoro(NLON,jlat+1)                   &
     &     .and. zoro(NLON,jlat) <= zoro(1,jlat)                        &
     &     .and. zoro(NLON,jlat) <= zoro(NLON-1,jlat)) then
          zoron(NLON,jlat)=1.+MIN(zoro(1,jlat),zoro(NLON-1,jlat)        &
     &                           ,zoro(NLON,jlat+1),zoro(NLON,jlat-1))
          jconv=jconv+1
         else
          zoron(NLON,jlat)=zoro(NLON,jlat)
         endif
        enddo
        do jlon=2,NLON-1
         if(zlsm(jlon,1) > 0.                                           &
     &    .and. zoro(jlon,1) <= zoro(jlon,2)                            &
     &    .and. zoro(jlon,1) <= zoro(jlon+1,1)                          &
     &    .and. zoro(jlon,1) <= zoro(jlon-1,1)) then
          zoron(jlon,1)=1.+MIN(zoro(jlon+1,1),zoro(jlon-1,1)            &
     &                        ,zoro(jlon,2))
          jconv=jconv+1
         else
          zoron(jlon,1)=zoro(jlon,1)
         endif
         if(zlsm(jlon,NLAT) > 0.                                        &
     &    .and. zoro(jlon,NLAT) <= zoro(jlon,NLAT-1)                    &
     &    .and. zoro(jlon,NLAT) <= zoro(jlon+1,NLAT)                    &
     &    .and. zoro(jlon,NLAT) <= zoro(jlon-1,NLAT)) then
          zoron(jlon,NLAT)=1.+MIN(zoro(jlon+1,NLAT),zoro(jlon-1,NLAT)   &
     &                           ,zoro(jlon,NLAT-1))
          jconv=jconv+1
         else
          zoron(jlon,NLAT)=zoro(jlon,NLAT)
         endif
        enddo
        if(zlsm(1,1) > 0.                                               &
     &    .and. zoro(1,1) <= zoro(1,2)                                  &
     &    .and. zoro(1,1) <= zoro(2,1)                                  &
     &    .and. zoro(1,1) <= zoro(NLON,1)) then
         zoron(1,1)=1.+MIN(zoro(2,1),zoro(NLON,1)                       &
     &                    ,zoro(1,2))
         jconv=jconv+1
        else
         zoron(1,1)=zoro(1,1)
        endif
        if(zlsm(NLON,NLAT) > 0.                                         &
     &    .and. zoro(NLON,NLAT) <= zoro(NLON,NLAT-1)                    &
     &    .and. zoro(NLON,NLAT) <= zoro(1,NLAT)                         &
     &    .and. zoro(NLON,NLAT) <= zoro(NLON-1,NLAT)) then
         zoron(NLON,NLAT)=1.+MIN(zoro(1,NLAT),zoro(NLON-1,NLAT)         &
     &                          ,zoro(NLON,NLAT-1))
         jconv=jconv+1
        else
         zoron(NLON,NLAT)=zoro(NLON,NLAT)
        endif
        if(zlsm(NLON,1) > 0.                                            &
     &    .and. zoro(NLON,1) <= zoro(NLON,2)                            &
     &    .and. zoro(NLON,1) <= zoro(1,1)                               &
     &    .and. zoro(NLON,1) <= zoro(NLON-1,1)) then
         zoron(NLON,1)=1.+MIN(zoro(1,1),zoro(NLON-1,1)                  &
     &                       ,zoro(NLON,2))
         jconv=jconv+1
        else
         zoron(NLON,1)=zoro(NLON,1)
        endif
        if(zlsm(1,NLAT) > 0.                                            &
     &    .and. zoro(1,NLAT) <= zoro(1,NLAT-1)                          &
     &    .and. zoro(1,NLAT) <= zoro(2,NLAT)                            &
     &    .and. zoro(1,NLAT) <= zoro(NLON,NLAT)) then
         zoron(1,NLAT)=1.+MIN(zoro(NLON,NLAT),zoro(2,NLAT)              &
     &                       ,zoro(1,NLAT-1))
         jconv=jconv+1
        else
         zoron(1,NLAT)=zoro(1,NLAT)
        endif
        zoro(:,:)=zoron(:,:)
        if(jconv > 0 ) goto 1000

       do jlat=1,NLAT
        do jlon=1,NLON-1
         zdx=TWOPI*cos(ASIN(zsi(jlon,jlat)))*plarad/real(NLON)
         zdh=(zoro(jlon+1,jlat)-zoro(jlon,jlat))/zdx
         zfac=1.
         if(zdh > 0.) zfac=-1.
         zuroff(jlon,jlat)=zfac*zcvel/zdx*ABS(zdh)**zcexp
        enddo
        zdx=TWOPI*cos(ASIN(zsi(NLON,jlat)))*plarad/real(NLON)
        zdh=(zoro(1,jlat)-zoro(NLON,jlat))/zdx
        zfac=1.
        if(zdh > 0.) zfac=-1.
        zuroff(NLON,jlat)=zfac*zcvel/zdx*ABS(zdh)**zcexp
       enddo

       do jlat=1,NLAT-1
        do jlon=1,NLON
         zdy=(ASIN(zsi(jlon,jlat+1))-ASIN(zsi(jlon,jlat)))*plarad
         zdh=(zoro(jlon,jlat+1)-zoro(jlon,jlat))/zdy
         zfac=1.
         if(zdh < 0.) zfac=-1.
         zvroff(jlon,jlat)=zfac*zcvel/zdy*ABS(zdh)**zcexp
        enddo
       enddo
       zvroff(1:NLON,NLAT)=0.
       
       endif

      endif

      call mpscgp(zuroff,duroff,1)
      call mpscgp(zvroff,dvroff,1)

      driver(:)=0.
      drunoff(:)=0.

      return
      end subroutine oroini           
      
!      
!     =============================
!
      subroutine glacierstep
      use glaciermod
      
      if (nglacier .eq. 1) then
      
      do jhor = 1,NHOR
         if (dls(jhor) > 0.5) then
            if (dsnowz(jhor) < glacelim) persistflag(jhor) = .FALSE. !Snowpack below minimum persistent threshhold
            if (dsnowz(jhor) < 0.1) dglac(jhor) = 0.0 !If the snow/ice is basically gone, so is the glacier
         endif
      enddo
      
      endif
      
      return
      end subroutine glacierstep
      
!     
!     ============================
!

      subroutine glacierstop
      use glaciermod

      real :: rpersist(NHOR) = 0.0
      
      if (nglacier .eq. 1) then
      
      call finishup(dsnowz,'newdsnow')
      call finishup(asndch,'restart_snow')
      
      where (persistflag(:)) dglac = 1.0
      where (persistflag(:)) rpersist = 1.0
      
      endif
      
      call finishup(dls,'lsm')
      call finishup(rpersist,'persist')
      
      call mpputgp('groundsg',groundoro ,NHOR,1)
      call mpputgp('dglacsg' ,glacieroro,NHOR,1)
      call mpputgp('doro'    ,doro      ,NHOR,1)
      
!       endif
      
      call mpputgp('dglac'   ,dglac     ,NHOR,1)

      return
      end subroutine glacierstop
      
!
!     ===============================


