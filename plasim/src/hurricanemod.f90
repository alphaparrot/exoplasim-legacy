
      module hurricanemod
      use resmod
      
      character (256) :: hc_indlog = "hurricane_indicators"
      character (256) :: hc_output = "hurricane_output"
      
      real :: capes(NHOR) = 0.  ! convective available potential energy at sea level at sat.
      real :: capeb(NHOR) = 0.  ! convective available potential energy of boundary layer
      real :: toutflow(NHOR) = 0. !Outflow temperature [K]
      real :: absvort(NHOR,NLEV) = 0. !Absolute vorticity (dimensionalized)
      real :: entropy(NHOR,NLEP) = 0. !Entropy
      real :: satentropy(NHOR,NLEP) = 0. !Saturation entropy
      real :: environtropy(NHOR,NLEV) = 0. !Environmental entropy
      real :: chim(NHOR) = 0. !Entropy deficit
      real :: shear(NHOR) = 0. !Change in velocity between sigma=0.1 and 0.85
      
      
! Thermodynamic Constants
      real :: CPD=1005.7       ! [J/kg.K] Specific heat of dry air at constant pressure
      real :: CPV=1870.0       ! [J/kg.K] Specific heat of water vapor at constant pressure
      real :: !CL=4190.0       ! [J/kg.K] Specific heat of liquid water
      real :: CL=2500.0        ! [J/kg.K] Modified specific heat of liquid water
      real :: CPVMCL=-630.0    ! CPV-CL
      real :: RVp=461.5         ! [J/kg.K] gas constant of water vapor
      real :: RD=287.04        ! [J/kg.K] gas constant of dry air
      real :: EPS=0.62197183  !RD/Rvp      ! [unitless] epsilon, the ratio of gas constants
      real :: ALV0=2.501e6     ! [J/kg] latent heat of vaporization at 0 degrees C

! pLCL Empirical Parameters
      real :: AA=1669.0
      real :: BB=122.0

! PI Auxiliaries
      real :: baz=2.0        ! Exponent for estimating azimuthal velocity in the eye, V=V_m(r/r_m)**b (Emanuel 1995, EQN. 25)
      real :: top=0.05      ! Pressure below which sounding is ignored (hPa)
      
      
      end module hurricanemod
      
!
!     ============================
!     SUBROUTINE HURRICANEINI
!     ============================
!
!     This routine computes the various hurricane formation probability metrics
!

      subroutine hurricaneini
      use hurricanemod
      
      
      
      end subroutine hurricaneini     
      
!
!     ============================
!     SUBROUTINE HURRICANESTEP
!     ============================
!
!     This routine computes the various hurricane formation probability metrics
!

      subroutine hurricanestep
      use hurricanemod
      
      
      
      end subroutine hurricanestep   
      
!
!     ============================
!     SUBROUTINE HURRICANESTOP
!     ============================
!
!     This routine computes the various hurricane formation probability metrics
!

      subroutine hurricanestop
      use hurricanemod
      
      
      
      end subroutine hurricanestop
      
!
! =================================================    
!    Convenience Functions
! ================================================= 
!
! ---------------------- Thermodynamic Calculations ---------------------- %

! Saturated water vapor pressure
! from Clausius-Clapeyron relation/August-Roche-Magnus formula
! Input temperature (TC) in Celsius
! Output saturation vapor pressure in hPa
      subroutine es_cc(TC,svp)
      real, intent(in ) TC
      real, intent(out) svp
      
      svp = 6.112*exp(17.67*TC/(243.5+TC))
      
      return
      end subroutine es_cc

! Latent heat of vaporization, as a function of temperature
! Input temperature (TC) in Celsius
! Output Latent heat of vaporization in J/kg
      subroutine Lv(TC,lvp)
      real, intent(in ) TC
      real, intent(out) lvp
      
      lvp=ALV0+CPVMCL*TC
      
      return 
      end subroutine Lv

! Parcel vapor pressure (hPa)
! Input mixing ratio (R) in gram/gram
! Input pressure (P) in hPA
      subroutine ev(R,P,pvp)
      use hurricanemod
      real, intent(in ) R
      real, intent(in ) P
      real, intent(out) pvp
      
      pvp=R*P/(EPS+R)
      
      return
      end subroutine ev

! Parcel mixing ratio (gram/gram)
! Input vapor pressure (E) in hPa
! Input pressure (P) in hPA
      subroutine rv(E,P,rvap)
      use hurricanemod
      real, intent(in ) E
      real, intent(in ) P
      real, intent(out) rvap
      
      rvap = EPS*E/(P-E)
      
      return
      end subroutine rv

! Total specific Entropy per unit mass of dry air (E94, EQN. 4.5.9)
! Input temperature (T) in kelvin
! Input mixing ratio (R) in gram/gram
! Input pressure (P) in hPA
      subroutine entropy_S(T,R,P,entropy)
      use hurricanemod
      real, intent(in ) T
      real, intent(in ) R
      real, intent(in ) P
      real, intent(out) entropy
      
      real evp,ES,RH
      
      ev(R,P,evp)
      es_cc(T-273.15,ES)
      RH=min(evp/ES,1.0)
      Lv(T-273.15,ALV)
      entropy=(CPD+R*CL)*log(T)-RD*log(P-evp)+ALV*R/T-R*RVp*log(RH)
      
      return
      end subroutine entropy_S

! Density temperature in K
! Input temperature (T) in kelvin
! Input total water content mixing ratio (RT) in gram/gram
! Input parcel water vapor mixing ratio (R) in gram/gram
      subroutine Trho(T,RT,R,rhoT)
      use hurricanemod
      real, intent(in ) T
      real, intent(in ) RT
      real, intent(in ) R
      real, intent(out) rhoT
      
      rhoT=T*(1.+R/EPS)/(1.+RT)

! Empirical pLCL
! Input parcel pressure (PP), temperature (TP), and relative humidity
! Output lifting condensation level pressure
      subroutine e_pLCL(TP,RH,PP,pLCL)
      use hurricanemod
      real TP,RH,PP,pLCL
      
      pLCL=PP*(RH**(TP/(AA-BB*RH-TP)))
      
      return
      end subroutine e_pLCL
      
      
!==============================================================================
!      Main analysis routines
!==============================================================================
! Compute convective available Potential Energy
      
      subroutine CAPE(TP,RP,PP,T,R,P,MSL,CAPED,TOB,LNB,IFLAG)
      use hurricanemod
      
      real,     intent(in ) :: TP
      real,     intent(in ) :: RP
      real,     intent(in ) :: PP
      real,     intent(in ) :: T(NLEV)
      real,     intent(in ) :: R(NLEV)
      real,     intent(in ) :: P(NLEV)
      real,     intent(in ) :: MSL
                
      real,     intent(out) :: CAPED
      real,     intent(out) :: TOB
      real,     intent(out) :: LNB
      integer, integer(out) :: IFLAG
      
      
      ! Internal variables
      
      real :: TVRDIF(NLEV) = 0.0
      
      real pdiff, pdiffnu, TPC, ESP, EVP, RH, S, PLCL, SL, SG
      real TG, RG, TLVR, TVENV, TGNEW, TJC, ES, ENEW, ALV, EM
      real AP, RMEAN, PFAC, PA, NA, PAT, PINB
      
      integer N, jlev, NCMAS, jmin, NC, INB

!     function [CAPED,TOB,LNB,IFLAG]= cape(TP,RP,PP,T,R,P,ascent_flag=0,ptop=50,miss_handle=1)
!
!       This function calculates the CAPE of a parcel given parcel pressure PP (hPa), 
!       temperature TP (K) and mixing ratio RP (gram/gram) and given a sounding
!       of temperature (T in K) and mixing ratio (R in gram/gram) as a function
!       of pressure (P in hPa). CAPED is the calculated value of CAPE following
!       Emanuel 1994 (E94) Equation 6.3.6 and TOB is the temperature at the
!       level of neutral buoyancy ("LNB") for the displaced parcel. IFLAG is a flag
!       integer. If IFLAG = 1, routine is successful; if it is 0, routine did
!       not run owing to improper sounding (e.g. no water vapor at parcel level).
!       IFLAG=2 indicates that the routine did not converge, IFLAG=3 indicates that
!       the input profile had missing values.         
!
!  INPUT:   TP,RP,PP: floating point numbers of Parcel pressure (hPa), 
!             temperature (K), and mixing ratio (gram/gram)
!
!           T,R,P: One-dimensional arrays 
!             containing environmental pressure (hPa), temperature (K),
!             and mixing ratio (gram/gram) profiles. The arrays MUST be
!             arranged so that the lowest index corresponds
!             to the lowest model level, with increasing index
!             corresponding to decreasing pressure.
!
!           ascent_flag: Adjustable constant fraction for buoyancy of displaced  
!             parcels, where 0=Reversible ascent;  1=Pseudo-adiabatic ascent
!
!           ptop: Pressure below which sounding is ignored (hPa)
!
!           miss_handle: Flag that determines how missing (NaN) values are handled.
!             If = 0 (BE02 default), NaN values in profile are ignored and PI is still calcuated
!             If = 1 (pyPI default), given NaN values PI will be set to missing (with IFLAG=3)
!             NOTE: If any missing values are between the lowest valid level and ptop
!             then PI will automatically be set to missing (with IFLAG=3)
!
!
!  OUTPUT:  CAPED (J/kg) is Convective Available Potential Energy of an air parcel
!             consistent with its parcel and environmental properties.
!
!           TOB is the Temperature (K) at the level of neutral bouyancy 
!             for the displaced air parcel
!
!           LNB is the pressure level of neutral bouyancy (hPa) for the 
!             displaced air parcel
!
!           IFLAG is a flag where the value of 1 means OK; a value of 0
!             indicates an improper sounding or parcel; a value of 2
!             means that the routine failed to converge
!
! Translated to Fortran for ExoPlaSim from Daniel Gilford's Python implementation of 
! Kerry Emanuel's Fortran code, by Adiv Paradise (I know this seems convoluted to 
! go from Fortran to Python to Fortran, but this is the best way to ensure correct
! indexing since PlaSim arrays are reversed relative to what the Emanuel code expects)
!


    ! Populate new environmental profiles removing values above ptop and
    ! find new number, N, of profile levels with which to calculate CAPE
    
      N = 1
      
      ptop = top*MSL
      
      pdiff = 1.0e4
      do jlev=1,NLEV
         pdiffnu = sigmah(jlev)-top
         if (pdiffnu .lt. pdiff) then
            pdiff = pdiffnu
            N = jlev
         endif
      enddo
!       P=P[N:]
!       T=T[N:]
!       R=R[N:]
!       nlvl=len(P)
!       TVRDIF = np.zeros((nlvl,))
      
      !
      !   ***  Run checks   ***
      !
      
      ! CHECK: Is the input parcel suitable? If not, return missing CAPE
      if ((RP < 1e-6) or (TP < 200)) then
          CAPED=0
          TOB=-1
          LNB=-1
          IFLAG=0
          ! Return the unsuitable values
      else
    
    !
    !  ***  Define various parcel quantities, including reversible   ***
    !  ***                       entropy, S                          ***
    !                         
       TPC=TP-273.15                 ! Parcel temperature in Celsius
       es_cc(TPC,ESP)                ! Parcel's saturated vapor pressure
       ev(RP,PP,EVP)                 ! Parcel's partial vapor pressure
       RH=EVP/ESP                              ! Parcel's relative humidity
       RH=min(RH,1.0)                        ! ensure that the relatively humidity does not exceed 1.0
       ! calculate reversible total specific entropy per unit mass of dry air (E94, EQN. 4.5.9)
       entropy_S(TP,RP,PP,S)
       
       
       !
       !   ***  Estimate lifted condensation level pressure, PLCL   ***
       !     Based on E94 "calcsound.f" code at http://texmex.mit.edu/pub/emanuel/BOOK/
       !     see also https://psl.noaa.gov/data/composites/day/calculation.html
       !
       !   NOTE: Modern PLCL calculations are made following the exact expressions of Romps (2017),
       !   see https://journals.ametsoc.org/doi/pdf/10.1175/JAS-D-17-0102.1
       !   and Python PLCL code at http://romps.berkeley.edu/papers/pubdata/2016/lcl/lcl.py
       !
       e_pLCL(TP,RH,PP,PLCL)
       
       ! Initial default values before loop
       CAPED=0
       TOB=T(NLEV)
       IFLAG=1
       ! Values to help loop
       NCMAX=0
       jmin=NLEV!int(1e6)
       
       !
       !   ***  Begin updraft loop   ***
       !
       
       ! loop over each level in the profile
       do jlev=NLEV,N,-1
           
           ! jmin is the index of the lowest pressure level evaluated in the loop
           !jmin=int(min([jmin,j]))
       
           !
           !   *** Calculate Parcel quantities BELOW lifted condensation level   ***
           !
           if (P(jlev) .ge. PLCL) then
               ! Parcel temperature at this pressure
               TG=TP*(P(jlev)/PP)**(RD/CPD)
               ! Parcel Mixing ratio
               RG=RP
               ! Parcel and Environmental Density Temperatures at this pressure (E94, EQN. 4.3.1 and 6.3.7)
               Trho(TG,RG,RG,TLVR)
               Trho(T(jlev),R(jlev),R(jlev),TVENV)
               ! Bouyancy of the parcel in the environment (Proxy of E94, EQN. 6.1.5)
               TVRDIF(jlev)=TLVR-TVENV
               
           !
           !   *** Calculate Parcel quantities ABOVE lifted condensation level   ***
           ! 
           else
               
               ! Initial default values before loop
               TGNEW=T(jlev)
               TJC=T(jlev)-273.15
               es_cc(TJC,ES)
               rv(ES,P(jlev),RG)
               
               !
               !   ***  Iteratively calculate lifted parcel temperature and mixing   ***
               !   ***                ratio for reversible ascent                    ***
               !
               
               ! set loop counter and initial condition
               NC=0
               TG=0
       
               ! loop until loop converges or bails out
               do while (abs(TGNEW-TG) .gt. 0.001)
               
                   ! Parcel temperature and mixing ratio during this iteration
                   TG=TGNEW
                   TC=TG-273.15
                   es_cc(TC,ENEW)
                   rv(ENEW,P(jlev),RG)
                   
                   ! increase iteration count in the loop
                   NC = NC+1
                   
                   !
                   !   ***  Calculate estimates of the rates of change of the entropy    ***
                   !   ***           with temperature at constant pressure               ***
                   !
       
                   Lv(TC,ALV)
                   ! calculate the rate of change of entropy with temperature, s_ell
                   SL=(CPD+RP*CL+ALV*ALV*RG/(RVp*TG*TG))/TG
                   ev(RG,P(jlev),EM)
                   ! calculate the saturated entropy, s_k, noting r_T=RP and
                   ! the last term vanishes with saturation, i.e. RH=1
                   SG=(CPD+RP*CL)*log(TG)-RD*log(P(jlev)-EM)+ALV*RG/TG
                   ! convergence speed (AP, step in entropy fraction) varies as a function of 
                   ! number of iterations
                   if (NC < 3) then
                       ! converge slowly with a smaller step
                       AP=0.3
                   else
                       ! speed the process with a larger step when nearing convergence
                       AP=1.0
                   endif
                   ! find the new temperature in the iteration
                   TGNEW=TG+AP*(S-SG)/SL
                   
                   !
                   !   ***   If the routine does not converge, set IFLAG=2 and bail out   ***
                   !
                   if ((NC > 500) .OR. (ENEW > (P(jlev)-1))) then
                       CAPED=0
                       TOB=T(NLEV)
                       LNB=P(NLEV)
                       IFLAG=2
                       ! Return the uncoverged values
                       RETURN
                   endif
                   
                   ! store the number of iterations
                   NCMAX=NC
               enddo
                   
               !
               !   *** Calculate buoyancy   ***
               !
               RMEAN=RP
               ! Parcel and Environmental Density Temperatures at this pressure (E94, EQN. 4.3.1 and 6.3.7)
               Trho(TG,RMEAN,RG,TLVR)
               Trho(T(jlev),R(jlev),R(jlev),TVENV)
               ! Bouyancy of the parcel in the environment (Proxy of E94, EQN. 6.1.5)
               TVRDIF(jlev)=TLVR-TVENV
           endif
       enddo
               
       
       !
       !  ***  Begin loop to find Positive areas (PA) and Negative areas (NA) ***
       !                  ***  and CAPE from reversible ascent ***
       NA=0.0
       PA=0.0
       
       !
       !   ***  Find maximum level (minimum pressure) of positive buoyancy, INB    ***
       !
       INB=NLEV
       do jlev=N,NLEV
           if (TVRDIF(jlev) > 0) INB=min(INB,jlev)
       enddo
           !once we find a pressure where this is true, that becomes what we keep 
       !print jmin,nlvl,INB        
       ! CHECK: Is the LNB higher than the surface? If not, return zero CAPE  
       if (INB .eq. NLEV) then
           CAPED=0
           TOB=T(NLEV)
           LNB=P(INB)
!            TOB=np.nan
           LNB=0
           ! Return the unconverged values
           RETURN
       
       ! if check is passed, continue with the CAPE calculation
       else
       
       !
       !   ***  Find positive and negative areas and CAPE  ***
       !                  via E94, EQN. 6.3.6)
       !
           do jlev=NLEV-1,INB,-1
               PFAC=RD*(TVRDIF(jlev)+TVRDIF(jlev+1))*(P(jlev+1)-P(jlev))/(P(jlev)+P(jlev+1))
               PA=PA+max(PFAC,0.0)
               NA=NA-min(PFAC,0.0)
           enddo
       
       !
       !   ***   Find area between parcel pressure and first level above it ***
       !
           PMA=(PP+P(NLEV))
           PFAC=RD*(PP-P(NLEV))/PMA
           PA=PA+PFAC*max(TVRDIF(NLEV),0.0)
           NA=NA-PFAC*min(TVRDIF(NLEV),0.0)
           
       !
       !   ***   Find residual positive area above INB and TO  ***
       !         and finalize estimate of LNB and its temperature
       !
           PAT=0.0
           TOB=T(INB)
           LNB=P(INB)
           if (INB > 1) then
               PINB=(P(INB-1)*TVRDIF(INB)-P(INB)*TVRDIF(INB-1))/(TVRDIF(INB)-TVRDIF(INB-1))
               LNB=PINB
               PAT=RD*TVRDIF(INB)*(P(INB)-PINB)/(P(INB)+PINB)
               TOB=(T(INB)*(PINB-P(INB-1))+T(INB-1)*(P(INB)-PINB))/(P(INB)-P(INB-1))
           endif
       
       !
       !   ***   Find CAPE  ***
       !            
           CAPED=PA+PAT-NA
           CAPED=max(CAPED,0.0)
           ! set the flag to OK if procedure reached this point
           IFLAG=1
           ! Return the calculated outputs to the above program level 
           RETURN
       endif
      endif
      
      RETURN
      end subroutine CAPE
      
      
 !===========================================================================
 !    Potential Intensity Calculation
 !===========================================================================
 
 ! Takes as input SST in C, mean sea-level pressure, air temperature in C, Specific
 ! humidity in g/g, and returns potential intensity (m/s), convective available potential
 ! energy (J/kg), and the level of neutral buoyancy (hPa)
 
      subroutine POTI(SSTC,MSL,P,TC,R,VMAX,CAPEA,LNBA)
      use hurricanemod
      real, intent(in ) :: SSTC
      real, intent(in ) :: MSL
      real, intent(in ) :: P(NLEV)
      real, intent(in ) :: TC(NLEV)
      real, intent(in ) :: R(NLEV)
      real, intent(out) :: VMAX
      real, intent(out) :: CAPEA
      real, intent(out) :: LNBA
      
      !Internal variables
      real SSTK, PMIN, IFL, TOF, OTL, RAT, RS0
      real tmin, ES0, CAPEA, CAPEM, TNB,LNB, TNBMS, LNBS
      real TVAV, CAT, PNEW, PMOLD, PM, PNEW, CATFAC,FAC
      real T(NLEV)
      
      integer NK, IFLAG, IFL, NP
      
 !     function [VMAX,PMIN,IFL,TO,OTL] = pi(SSTC,MSL,P,TC,R,CKCD=0.9,ascent_flag=0,diss_flag=1,V_reduc=0.8,ptop=50,miss_handle=0)
 !   
 !   ***    This function calculates the maximum wind speed         ***
 !   ***             and mimimum central pressure                   ***
 !   ***    achievable in tropical cyclones, given a sounding       ***
 !   ***             and a sea surface temperature.                 ***
 !
 !   Thermodynamic and dynamic technical backgrounds (and calculations) are found in Bister 
 !   and Emanuel (2002; BE02) and Emanuel's "Atmospheric Convection" (E94; 1994; ISBN: 978-0195066302)
 !
 !  INPUT:   SSTC: Sea surface temperature (C)
 !
 !           MSL: Mean Sea level pressure (hPa)
 !
 !           P,TC,R: One-dimensional arrays 
 !             containing pressure (hPa), temperature (C),
 !             and mixing ratio (g/kg). The arrays MUST be
 !             arranged so that the lowest index corresponds
 !             to the lowest model level, with increasing index
 !             corresponding to decreasing pressure. The temperature
 !             sounding should extend to at least the tropopause and 
 !             preferably to the lower stratosphere, however the
 !             mixing ratios are not important above the boundary
 !             layer. Missing mixing ratios can be replaced by zeros
 !
 !           CKCD: Ratio of C_k to C_D (unitless number), i.e. the ratio
 !             of the exchange coefficients of enthalpy and momentum flux
 !             (e.g. see Bister and Emanuel 1998, EQN. 17-18). More discussion
 !             on CK/CD is found in Emanuel (2003). Default is 0.9 based
 !             on e.g. Wing et al. (2015)
 !
 !           ascent_flag: Adjustable constant fraction (unitless fraction) 
 !             for buoyancy of displaced parcels, where 
 !             0=Reversible ascent (default) and 1=Pseudo-adiabatic ascent
 !
 !           diss_flag: Adjustable switch integer (flag integer; 0 or 1)
 !             for whether dissipative heating is 1=allowed (default) or 0=disallowed.
 !             See Bister and Emanuel (1998) for inclusion of dissipative heating.
 !
 !           V_reduc: Adjustable constant fraction (unitless fraction) 
 !             for reduction of gradient winds to 10-m winds see 
 !             Emanuel (2000) and Powell (1980). Default is 0.8
 !
 !           ptop: Pressure below which sounding is ignored (hPa)
 !
 !           miss_handle: Flag that determines how missing (NaN) values are handled in CAPE calculation
 !             If = 0 (BE02 default), NaN values in profile are ignored and PI is still calcuated
 !             If = 1, given NaN values PI will be set to missing (with IFLAG=3)
 !             NOTE: If any missing values are between the lowest valid level and ptop
 !             then PI will automatically be set to missing (with IFLAG=3)
 !
 !  OUTPUT:  VMAX is the maximum surface wind speed (m/s)
 !             reduced to reflect surface drag via V_reduc
 !
 !           PMIN is the minimum central pressure (hPa)
 !
 !           IFL is a flag: A value of 1 means OK; a value of 0
 !             indicates no convergence; a value of 2
 !             means that the CAPE routine failed to converge;
 !             a value of 3  means the CAPE routine failed due to
 !             missing data in the inputs
 !
 !           TO is the outflow temperature (K)
 !
 !           OTL is the outflow temperature level (hPa), defined as the level of neutral bouyancy 
 !             where the outflow temperature is found, i.e. where buoyancy is actually equal 
 !             to zero under the condition of an air parcel that is saturated at sea level pressure
 !

     !if P[1]>P[0]:
     !    P = P[::-1]
     !    TC = TC[::-1]
     !    R = R[::-1]

     ! convert units
     SSTK=SSTC+273.15 ! SST in kelvin
     T(:)=TC(:)+273.15      ! Temperature profile in kelvin
     !R=R*0.001                   ! Mixing ratio profile in g/g

     ! CHECK 1: do SSTs exceed 5C? If not, set IFL=0 and return missing PI
     if (SSTC .le. 5.0) then
         VMAX=-1
         PMIN=-1
         IFL=0
         TOF=-1
         OTL=-1
         RETURN
     endif

     ! CHECK 2: do Temperature profiles exceed 100K? If not, set IFL=0 and return missing PI
     
     tmin = 1.0e4
     do jlev=1,NLEV
        tmin = min(tmin,T(jlev))
     enddo
     
     if (tmin .le. 100) then
         VMAX=-1
         PMIN=-1
         IFL=0
         TOF=-1
         OTL=-1
         RETURN
     endif
     
     ! Saturated water vapor pressure
     ! from Clausius-Clapeyron relation/August-Roche-Magnus formula
     es_cc(SSTC,ES0)

     ! define the level from which parcels lifted (first pressure level)
     NK=NLEV
     
     !
     !   ***   Find environmental CAPE *** 
     !
     TP=T(NK)
     RP=R(NK)
     PP=P(NK)
     CAPE(TP,RP,PP,T,R,P,MSL,CAPEA,TNB,LNBA,IFLAG)
     LNB = LNBA
     
     ! if the CAPE function tripped a flag, set the output IFL to it
     if (IFLAG .ne. 1) IFL=IFLAG
     
     !
     !   ***   Begin iteration to find mimimum pressure   ***
     !
     
     ! set loop counter and initial condition
     NP=0         ! loop counter
     PM=970.0
     PMOLD=PM     ! initial condition from minimum pressure
     PNEW=0.0     ! initial condition from minimum pressure
     IFL=1        ! Default flag for CAPE calculation

     ! loop until convergence or bail out
     do while (abs(PNEW-PMOLD) .gt. 0.5)
         
         !
         !   ***  Find CAPE at radius of maximum winds   ***
         !
         TP=T(NK)
         PP=min(PM,1000.0)
         ! find the mixing ratio with the average of the lowest level pressure and MSL
         RP=EPS*R(NK)*MSL/(PP*(EPS+R(NK))-R(NK)*MSL)
         
         CAPE(TP,RP,PP,T,R,P,MSL,CAPEM,TNB,LNB,IFLAG)
         ! if the CAPE function tripped a different flag, set the output IFL to it
         if (IFLAG .ne. 1) IFL=IFLAG
         
         !
         !  ***  Find saturation CAPE at radius of maximum winds    ***
         !  *** Note that TO and OTL are found with this assumption ***
         !
         TP=SSTK
         PP=min(PM,1000.0)
         rv(ES0,PP,RP)
         
         CAPE(TP,RP,PP,T,R,P,MSL,CAPEMS,TNBMS,LNBS,IFLAG)
         ! if the CAPE function tripped a flag, set the output IFL to it
         if (IFLAG .ne. 1) IFL=IFLAG
         ! Store the outflow temperature and level of neutral bouyancy at the outflow level (OTL)
         TOF=TNBMS   
         OTL=LNBS
         ! Calculate the proxy for TC efficiency (BE02, EQN. 1-3)
         RAT=SSTK/TOF
         
         !
         !  ***  Initial estimate of pressure at the radius of maximum winds  ***
         !
         RS0=RP
         ! Lowest level and Sea-surface Density Temperature (E94, EQN. 4.3.1 and 6.3.7)
         Trho(T(NK),R(NK),R(NK),TV0)
         Trho(SSTK,RS0,RS0,TVSST)
         ! Average Surface Density Temperature, e.g. 1/2*[Tv(Tsfc)+Tv(sst)]
         TVAV=0.5*(TV0+TVSST)
         ! Converge toward CAPE*-CAPEM (BE02, EQN 3-4)
         CAT=(CAPEM-CAPEA)+0.5*CKCD*RAT*(CAPEMS-CAPEM)
         CAT=max(CAT,0.0)
         ! Iterate on pressure
         PNEW=MSL*exp(-CAT/(RD*TVAV))
         
         !
         !   ***  Test for convergence (setup for possible next while iteration)  ***
         !
         ! store the previous step's pressure       
         PMOLD=PM
         ! store the current step's pressure
         PM=PNEW
         ! increase iteration count in the loop
         NP = NP+1
         
         !
         !   ***   If the routine does not converge, set IFL=0 and return missing PI   ***
         !
         if ((NP > 200)  .OR. (PM < 400)) then
             VMAX=-1
             PMIN=-1
             IFL=0
             TOF=-1
             OTL=-1
             RETURN
         endif
     enddo
     
     ! Once converged, set potential intensity at the radius of maximum winds
     CATFAC=0.5*(1.+1/baz)
     CAT=(CAPEM-CAPEA)+CKCD*RAT*CATFAC*(CAPEMS-CAPEM)
     CAT=max(CAT,0.0)
     
     ! Calculate the minimum pressure at the eye of the storm
     ! BE02 EQN. 4
     PMIN=MSL*exp(-CAT/(RD*TVAV))
                  
     ! Calculate the potential intensity at the radius of maximum winds
     ! BE02 EQN. 3, reduced by some fraction (default 20%) to account for the reduction 
     ! of 10-m winds from gradient wind speeds (Emanuel 2000, Powell 1980)
     FAC=max(0.0,(CAPEMS-CAPEM))
     VMAX=sqrt(CKCD*RAT*FAC)
         
     ! Return the calculated outputs to the above program level
     RETURN
     end subroutine POTI
