#define kasting 0
#define realv 1
#define keeprecord 0
#define outputrain 0
#define precipweath 1

      module carbonmod
      use pumamod
!
!     version identifier (date)
!
      character(len=80) :: version = '07.14.2017 by Adiv'
!
!     Parameter   
!                          !Outgassing rates given in ubars/year
! #if realv == 1
!       parameter(VEARTH=5.0e-2) ! Rate based on volcanic outgassing estimates from Gerlach 2011
! #else
! #if kasting==1
!       parameter(VEARTH=6.3e-2)
! #else
!       parameter(VEARTH=7.0e-3) ! Earth volcanic CO2 outgassing rate (ubars), assumed to be equal to Earth
!                                ! weathering rate.
! #endif
! #endif
      parameter(CO2EARTH=330.0)  ! Earth CO2 level in ubars
      parameter(RAD_EARTH=6371220.0) !Earth radius
      parameter(RAD_EARTHSQ=6371220.0*6371220.0) !Earth radius squared
!
!     namelist parameters
!

      character (256) :: carbon_namelist     = "carbonmod_namelist"
      
      integer :: ncarbon = 1 ! 1 = do weathering, 0 = constant pCO2
      integer :: nsupply = 0 ! Should it be supply-limited weathering?
      integer :: nco2evolve = 0 !Should CO2 and pressure evolve each year? (0/1)
      real :: volcanco2 = 0.02 ! Volcanic outgassing rate in units of VEARTH
      real :: kact = 0.09 ! activation energy factor
      real :: krun = 0.045 ! runoff efficiency factor
      real :: beta = 0.5 ! pCO2 dependence
      real :: frequency = 4.0 ! Number of times to compute weathering per day. Can be a float, i.e.
                              ! for once every 2 days, use frequency=0.5.
      real :: VEARTH = 5.0e-2 !ubars/year
      real :: PEARTH = 79.0 ! Annual precipitation on modern Earth that is relevant for weathering.
                           ! 79 cm/yr (Chen 2002 & Schneider 2014)
      real :: WMAX = 1.0 ! Maximum weathering rate for the supply-limited case, in ubar/yr
      real :: zeta = 0.0 !Dependence of max weathering on precipitation (not used)
      
      real :: psurf0 = 101100.0 !Mean sea-level pressure
!
!     global arrays
!
      real :: dglobe(NHOR) = 0. ! Normalized cell area

      real :: localweathering(NHOR) = 0.0 ! Local instantaneous weathering rate relative to Earth standard
      real :: localavgweather(NHOR) = 0.0 ! Annual average weathering rate relative to Earth standard
      real :: localavgtemps(NHOR) = 0.0 ! Annual average surface temperature
      real :: localprecip(NHOR) = 0.0 !Local precipitation
      real :: localavgprecip(NHOR) = 0.0!Local net precipitation
      
      real :: aweathering(NHOR) = 0.0 !Monthly average weathering
      
#if outputrain==1      
      real :: cprecip(NHOR) = 0.0
#endif

!
!     global scalars
!
      integer :: interval = 8   ! Number of timesteps to go between weathering updates
      integer :: cstep = 0 ! Current timestep number in the cycle. When cstep = interval, do weathering
      integer :: istep = 0
!       integer :: nco2evolve = 0 ! Do we allow CO2 to evolve year-to-year?
      real :: timeweight = 0.0 ! Weight for computing annual averages
      real :: avgweathering = 0.0 ! Global annual average weathering rate
      real :: dpco2dt = 0.0 ! Change in pCO2 with respect to time. = (volcanco2 - avgweathering)*VEARTH
      real :: tune1 = 5.41 ! Tuning adjustment to make global average match for non-precip model.
      real :: tune2 = 2.20 ! Tuning adjustment to make precip model match non-precip model.
!
      end module carbonmod
      
!      
!     =============================
!

      subroutine carbonini
      use carbonmod
      
      namelist/carbonmod_nl/ncarbon,volcanco2,kact,krun,beta,frequency,VEARTH,PEARTH,nsupply,WMAX,zeta,nco2evolve
      
      if (mypid==NROOT) then
         open(23,file=carbon_namelist)
         read(23,carbonmod_nl)
         close(23)
         write(nud,'(/," *********************************************")')
         write(nud,'(" * CARBONMOD ",a34," *")') trim(version)
         write(nud,'(" *********************************************")')
         write(nud,'(" * Namelist CARBON_NL from <carbon_namelist> *")')
         write(nud,'(" *********************************************")')
         write(nud,carbonmod_nl)
      endif
      
      frequency = min(frequency,real(mtspd))
      interval = int(mtspd/frequency)
      timeweight = (1.0 / (frequency*m_days_per_year))
      
      psurf0 = psurf
      
      call makeareas
      
      call mpbci(interval)
      call mpbci(ncarbon)
      call mpbcr(volcanco2)
      call mpbcr(kact)
      call mpbcr(krun)
      call mpbcr(beta)
      call mpbcr(frequency)
      call mpbcr(timeweight)
      call mpbcr(PEARTH)
      call mpbcr(VEARTH)
      call mpbcr(WMAX)
      call mpbcr(nsupply)
      call mpbcr(psurf0)
      
      end subroutine carbonini
      
!      
!     =============================
!
      subroutine carbonstep
      use carbonmod
      use rainmod
      
      real tsurf(NHOR)
      real pco2
      real cweathering
      real landfraction
      real globalavgt

#if outputrain==1      
      character srainstep*5   
      real :: fullprecip(NUGP) = 0.0
      real :: fullweather(NUGP) = 0.0
      real :: fullevap(NUGP) = 0.0
      real :: fullqvi(NUGP) = 0.0
#endif   
      
      integer i
      
      if (ntime > 0) call mksecond(zsec,0.)
      
      if (mypid == NROOT) cstep = cstep + 1
      if (mypid == NROOT) istep = istep + 1
      
      call mpbci(cstep)
      
      localprecip(:) = (dprc(:) + dprl(:))*3.154e9 !cm/yr
      localavgprecip(:) = localavgprecip(:) + localprecip(:)*timeweight
      
      
      if (cstep .eq. interval) then
      
             
         if (mypid == NROOT) cstep = 0
         
         tsurf(:) = dt(:,NLEP)
         
         if (ncarbon .gt. 0.5) then !we're doing this
         
         localweathering(:) = 0.
         
         do i=1,NHOR
            if (dls(i) .gt. 0.5) then !land
               if (tsurf(i) .ge. 273.15) then !not frozen
                  pco2 = co2*dp(i)*1e-5 ! Convert ppmv to ubars (1e-6 for ppmv->ppv, and 10 for Pa->ubars)
#if precipweath == 1
                  localweathering(i) = ((pco2/CO2EARTH)**beta) * exp(kact*(tsurf(i)-288.0)) * &
     &                                 (localprecip(i)/PEARTH)**0.65 !* tune1*tune2 
#else
                  localweathering(i) = ((pco2/CO2EARTH)**beta) * exp(kact*(tsurf(i)-288.0)) * &
                                       (1+krun*(tsurf(i)-288.0))**0.65 * tune1
#endif
                  if (nsupply .eq. 1) then !Apply a weathering supply limit following Foley 2015
                     wmaxp = max((WMAX/VEARTH),2.0e-15)
                     localweathering(i) = wmaxp*(1-exp(-localweathering(i)/wmaxp))
                  endif
                  
               endif
            endif
#if outputrain==1
            cprecip(i) = (dprc(i) + dprl(i))
#endif
         enddo
         
         
         localweathering(:) = localweathering(:)*plarad**2 / (RAD_EARTHSQ*0.29) !Account for land fraction
         
!          endif
         
#if outputrain==1     
         call mpgagp(fullprecip,cprecip,1)
         call mpgagp(fullweather,localweathering*timeweight,1)
         call mpgagp(fullevap,devap,1)
         call mpgagp(fullqvi,dqvi,1)
         if (mypid == NROOT) then
            write(srainstep,'(I5.5)') istep
            call writegtextarray(fullprecip,NUGP,'precipmp'//trim(adjustl(srainstep)))
            call writegtextarray(fullweather,NUGP,'weathrmp'//trim(adjustl(srainstep)))
            call writegtextarray(fullevap,NUGP,'evaprtmp'//trim(adjustl(srainstep)))
            call writegtextarray(fullqvi,NUGP,'humdtymp'//trim(adjustl(srainstep)))
         endif
#endif        
        
         endif
         
         
#if keeprecord==1
         call write_short(localweathering,cweathering)
         call write_short(dls,landfraction)
         if (mypid==NROOT) then
            write(nud,'(" Land Fraction        = ",f10.2," [...]")') landfraction
            write(nud,'(" Average Weathering   = ",f10.2," [earth]")') cweathering
            cweathering = cweathering / (landfraction * PEARTH)
            open(unit=44,file="weathertrack.pso",position="append",status="unknown")
            write(44,'(1p8e13.5)') cweathering*ncarbon*VEARTH
            close(44)
         endif
#endif
         
          
         localavgtemps(:) = localavgtemps(:) + tsurf(:)*timeweight
         localavgweather(:) = localavgweather(:) + localweathering(:)*timeweight
         call mpbci(cstep)
         
         
      endif
      
      if (ntime > 0) then
         call mksecond(zsec,zsec)
         time4co2 = time4co2 + zsec
         write(nud,*)
         write(nud,*)'Carbon cycle routines: ',zsec,'s'
      endif
      
      return
      end subroutine carbonstep
      
!     
!     ============================
!

      subroutine carbonstop
      use carbonmod
      
      real landfraction
      real newco2
      real pco2
      real newpco2
      real globalavgt
      real :: globalweath(NUGP) = 0.0
      
      
      call mpgagp(globalweath,localavgweather,1)
      
      call write_short(localavgweather,avgweathering)
!       avgweathering = avgweathering / PEARTH
      call write_short(localavgtemps,globalavgt)
      call write_short(dls,landfraction)
      
      
      !
      ! We only computed weathering on land, but our global average included the ocean, so if we 
      ! used Earth conditions, we would so far have a weathering rate of W_earth * (land fraction). 
      ! Therefore, we have to divide by the land fraction. If we were to include seafloor weathering,
      ! then we would divide by the present-day (landfraction+seafraction), where the fraction that is
      ! sea ice is excluded. This way, the advance of sea ice would actually reduce the global weathering
      ! regardless of temperature--the fraction of the Earth undergoing weathering would go down from 
      ! present-day.
      !
      !
      ! Nevermind. We seem to get the right answer without dividing by the land fraction...
      ! Maybe the equation is already parameterized to account for land fraction...?
      !
      ! From Abbot, et al.... this isn't the global weathering rate; it's the continental weathering rate.
      ! This can be modulated by reduction in land due to glacial advance or sea level rise/fall, but the
      ! overall weathering rate is this PLUS the seafloor weathering rate, modulated by changes in sea surface
      ! area.
      !
      ! No, we actually should divide by the land fraction. Because we're summing area-weighted and time-weighted 
      ! estimates of the global weathering rate due to continental weathering--each cell's value before area-weighting
      ! is an estimate of what the global weathering rate would be if every cell on the planet had that weathering
      ! rate. So when we do the global average, if we don't divide by the land fraction, then we're actually giving
      ! the continental weathering rate multiplied by the land fraction. We have to divide to get the global overall
      ! weathering rate (factor of 3.33 increase, seems okay based on 1 AU tests) due to continental weathering. If 
      ! we were to add in seafloor weathering, then we'd add factors to account for changing land and sea area, but
      ! there would also be weights to determine the relative contributions of each, and those are parameterized (poorly).
!       if (mypid==NROOT) then 
!          if (landfraction > 0) then
!            avgweathering = avgweathering / landfraction
!          else
!            avgweathering = 0.0
!          endif
!       endif

      ! We should divide by Earth's land fraction as a constant pre-factor, and then upon integrating
      ! over the whole planet, we'll pick up a multiplicative factor of the actual land fraction.
      
      if (mypid==NROOT) then
         
         globalweath(:) = globalweath(:)*VEARTH*1000.0
         call writegtextarray(globalweath,NUGP,'annualweather')
         dpco2dt = ncarbon*VEARTH*(volcanco2 - avgweathering)
         pco2 = co2*1e-6*psurf0 ! Pa
         newco2 = (pco2 + dpco2dt*0.1)/(psurf0+dpco2dt*0.1) ! ppv [Pa/Pa]
         newpco2 = newco2*(psurf0+dpco2dt*0.1)*1.e-5 ! Convert to bars (from Pa)
         open(unit=43,file="weathering.pso",position="append",status="unknown")
         write(43,'(1p8e13.5)') pco2*1e-5,globalavgt,avgweathering,volcanco2,dpco2dt*1e-6,newpco2
         close(43)
         co2 = newco2*1e6 ! ppmv
         psurf = psurf0 + dpco2dt*0.1 ! Pa
      endif
      call mpbcr(avgweathering)
      call mpbcr(dpco2dt)
      call mpbcr(co2)
      call mpbcr(psurf)
      
      n_run_months = 0
      call mpbci(n_run_months)
      if (nco2evolve > 0.5) then
        call co2update !Change CO2 for following year
        call psurfupdate !Change surface pressure for following year
      endif
      
      return
      end subroutine carbonstop
      
!
!     ===============================
!
      
      subroutine co2update
      use radmod
      
      namelist/radmod_nl/ndcycle,ncstsol,solclat,solcdec,no3,co2        &
     &               ,iyrbp,nswr,nlwr,nfixed,fixedlon                   &
     &               ,a0o3,a1o3,aco3,bo3,co3,toffo3,o3scale             &
     &               ,nsol,nswrcl,nrscat,rcl1,rcl2,acl2,clgray,tpofmt   &
     &               ,acllwr,tswr1,tswr2,tswr3,th2oc,dawn
      if (mypid == NROOT) then
      open(23,file=radmod_namelist)
      write(23,radmod_nl)
      close(23)
      endif
      return
      end subroutine co2update
      
!
!     ===============================
!
   
      subroutine psurfupdate
      use pumamod
      
      namelist /plasim_nl/ &
     &               noutput,ngui,n_start_year,       &
     &               n_days_per_year,n_run_years,n_run_months,n_run_days,   &
     &               kick,mpstep,naqua,ndiag,nguidbg,nqspec,  &  
     &               nveg,nwpd,nprint,nsync,syncstr,psurf
      if (mypid == NROOT) then
      open(23,file=plasim_namelist,form='formatted')
      write(23,plasim_nl)
      close(23)
      endif
      return
      end subroutine psurfupdate
                   

      
