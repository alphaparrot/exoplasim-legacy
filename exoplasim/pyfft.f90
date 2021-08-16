

      
      ! =================
      ! SUBROUTINE LEGINI
      ! =================
      
      subroutine legini(NLAT,NLON,NTRU,NLEV,qi,qj,qc,qe,qm,qq,qu,qv,&
    &                   sfilt,nfs)
      implicit none
      
      integer :: jlat ! Latitude
      integer :: lm
      integer :: m
      integer :: n
      
      integer :: nfs
      
      integer :: NLAT, NLON, NTRU, NLEV
      
      real (kind=8) :: EZ
      real (kind=8) :: PI, TWOPI
      real (kind=8) :: amsq
      real (kind=8) :: z1
      real (kind=8) :: z2
      real (kind=8) :: z3
      real (kind=8) :: f1m
      real (kind=8) :: f2m
      real (kind=8) :: znn1
      real (kind=8) :: zsin    ! sin
      real (kind=8) :: zcsq    ! cos2
      real (kind=8) :: zgwd    ! gw
      real (kind=8) :: zgwdcsq ! gw / cos2
      
      real (kind=8) :: zpli((NTRU+1)*(NTRU+2)/2)
      real (kind=8) :: zpld((NTRU+1)*(NTRU+2)/2)
      
      real (kind=8) :: sid(NLAT)     ! sin(phi)
      real (kind=8) :: gwd(NLAT)     ! Gaussian weights
      real (kind=8) :: csq(NLAT)     ! cos(phi)**2
      real (kind=8) :: cola(NLAT)    ! cos(phi)
      real (kind=8) :: rcs(NLAT)     ! 1 / cos(phi)
      real (kind=8) :: deglat(NLat)  ! latitude in degrees
            
      real (kind=8):: qi((NTRU+1)*(NTRU+2)/2,NLAT) ! P(m,n) = Associated Legendre Polynomials
      real (kind=8):: qj((NTRU+1)*(NTRU+2)/2,NLAT) ! Q(m,n) = Used for d/d(mu)
      real (kind=8):: qc((NTRU+1)*(NTRU+2)/2,NLAT) ! P(m,n) * gwd              used in fc2sp
      real (kind=8):: qe((NTRU+1)*(NTRU+2)/2,NLAT) ! Q(mn,) * gwd / cos2       used in mktend
      real (kind=8):: qm((NTRU+1)*(NTRU+2)/2,NLAT) ! P(m,n) * gwd / cos2 * m   used in mktend
      real (kind=8):: qq((NTRU+1)*(NTRU+2)/2,NLAT) ! P(m,n) * gwd / cos2 * n * (n+1) / 2  "
      real (kind=8):: qu((NTRU+1)*(NTRU+2)/2,NLAT) ! P(m,n) / (n*(n+1)) * m    used in dv2uv
      real (kind=8):: qv((NTRU+1)*(NTRU+2)/2,NLAT) ! Q(m,n) / (n*(n+1))        used in dv2uv
          
      real (kind=8):: sfilt(NTRU+1) ! Physics filter
      
      integer NLPP, NHOR, NUGP, NPGP, NLEM, NLEP, NLSQ, NTP1
      integer NRSP, NCSP, NSPP, NESP, NVCT
          
      EZ     = 1.63299310207D0
      PI     = 3.14159265359D0
      TWOPI  = PI * PI
      
      NLPP = NLAT         ! Latitudes per process
      NHOR = NLON * NLPP        ! Horizontal part
      NUGP = NLON * NLAT        ! Number of gridpoints
      NPGP = NLON * NLAT / 2    ! Dimension of packed fields
      NLEM = NLEV - 1           ! Levels - 1
      NLEP = NLEV + 1           ! Levels + 1
      NLSQ = NLEV * NLEV        ! Levels squared
      NTP1 = NTRU + 1           ! Truncation + 1
      NRSP =(NTRU+1)*(NTRU+2)   ! No of real global    modes
      NCSP = NRSP / 2           ! No of complex global modes
      NSPP = (NRSP+1-1)/1 ! Modes per process
      NESP = NSPP * 1        ! Dim of spectral fields
      NVCT = 2 * (NLEV+1)       ! Dim of Vert. Coord. Tab
      
      call inigau(NLAT,sid,gwd)
      
      do jlat = 1 , NLPP
      
      ! set p(0,0) and p(0,1)
      
         zgwd    = gwd(jlat)            ! gaussian weight - from inigau
         zsin    = sid(jlat)            ! sin(phi) - from inigau
         zcsq    = 1.0_8 - zsin * zsin  ! cos(phi) squared
         zgwdcsq = zgwd / zcsq          ! weight / cos squared
         f1m     = sqrt(1.5_8)
         zpli(1) = sqrt(0.5_8)
         zpli(2) = f1m * zsin
         zpld(1) = 0.0
         lm      = 2
      
      ! loop over wavenumbers
      
         do m = 0 , NTRU
            if (m > 0) then
               lm  = lm + 1
               f2m = -f1m * sqrt(zcsq / (m+m))
               f1m =  f2m * sqrt(m+m + 3.0_8)
               zpli(lm) = f2m
               if (lm < NCSP) then
                  lm = lm + 1
                  zpli(lm  ) =       f1m * zsin
                  zpld(lm-1) =  -m * f2m * zsin
               endif ! (lm < NCSP)
            endif ! (m > 0)
      
            amsq = m * m
      
            do n = m+2 , NTRU
               lm = lm + 1
               z1 = sqrt(((n-1)*(n-1) - amsq) / (4*(n-1)*(n-1)-1))
               z2 = zsin * zpli(lm-1) - z1 * zpli(lm-2)
               zpli(lm  ) = z2 * sqrt((4*n*n-1) / (n*n-amsq))
               zpld(lm-1) = (1-n) * z2 + n * z1 * zpli(lm-2)
            enddo ! n
      
            if (lm < NCSP) then ! mode (m,NTRU)
               z3 = sqrt((NTRU*NTRU-amsq) / (4*NTRU*NTRU-1))
               zpld(lm)=-NTRU*zsin*zpli(lm) +  &
     &                  (NTRU+NTRU+1)*zpli(lm-1)*z3
            else                ! mode (NTRU,NTRU)
               zpld(lm)=-NTRU*zsin*zpli(lm)
            endif
         enddo ! m
      
         lm = 0
         do m = 0 , NTRU
            do n = m , NTRU
                 lm = lm + 1
                 znn1 = 0.0
                 if (n > 0) znn1 = 1.0_8 / (n*(n+1))
                 qi(lm,jlat) = zpli(lm)
                 qj(lm,jlat) = zpld(lm)
                 qc(lm,jlat) = zpli(lm) * zgwd
                 qu(lm,jlat) = zpli(lm) * znn1 * m
                 qv(lm,jlat) = zpld(lm) * znn1
                 qe(lm,jlat) = zpld(lm) * zgwdcsq
                 qq(lm,jlat) = zpli(lm) * zgwdcsq * n * (n+1) * 0.5_8
                 qm(lm,jlat) = zpli(lm) * zgwdcsq * m
            enddo ! n
         enddo ! m
         
      enddo! jlat
      
      sfilt(:) = 1.0
         
      do n=1,NTP1
         sfilt(n) = (1-nfs)*sfilt(n)+nfs*exp(-8*(real(n)/NTRU)**8)
      enddo
      
      return
      end
            
      
      ! ================
      ! SUBROUTINE FC2SP
      ! ================
      
      subroutine fc2sp(fc,sp,NLAT,NLON,NTRU,NLEV,nfs)
      implicit none
      
      integer, intent(in ) :: NLAT
      integer, intent(in ) :: NLON
      integer, intent(in ) :: NLEV
      integer, intent(in ) :: NTRU
      integer, intent(in ) :: nfs
      
      real (kind=8), intent(in ) :: fc(2,NLON/2,NLAT)
!f2py intent(in) :: fc
      real (kind=8), intent(out) :: sp(2,(NTRU+1)*(NTRU+2)/2)
!f2py intent(out) :: sp
      
      integer :: l ! Index for latitude
      integer :: m ! Index for zonal wavenumber
      integer :: n ! Index for total wavenumber
      integer :: w ! Index for spherical harmonic
       
      real (kind=8):: qi((NTRU+1)*(NTRU+2)/2,NLAT) ! P(m,n) = Associated Legendre Polynomials
      real (kind=8):: qj((NTRU+1)*(NTRU+2)/2,NLAT) ! Q(m,n) = Used for d/d(mu)
      real (kind=8):: qc((NTRU+1)*(NTRU+2)/2,NLAT) ! P(m,n) * gwd              used in fc2sp
      real (kind=8):: qe((NTRU+1)*(NTRU+2)/2,NLAT) ! Q(mn,) * gwd / cos2       used in mktend
      real (kind=8):: qm((NTRU+1)*(NTRU+2)/2,NLAT) ! P(m,n) * gwd / cos2 * m   used in mktend
      real (kind=8):: qq((NTRU+1)*(NTRU+2)/2,NLAT) ! P(m,n) * gwd / cos2 * n * (n+1) / 2  "
      real (kind=8):: qu((NTRU+1)*(NTRU+2)/2,NLAT) ! P(m,n) / (n*(n+1)) * m    used in dv2uv
      real (kind=8):: qv((NTRU+1)*(NTRU+2)/2,NLAT) ! Q(m,n) / (n*(n+1))        used in dv2uv
           
      integer NLPP, NHOR, NUGP, NPGP, NLEM, NLEP, NLSQ, NTP1
      integer NRSP, NCSP, NSPP, NESP, NVCT
      
      real (kind=8) :: EZ
      real (kind=8) :: PI
      real (kind=8) :: TWOPI
      
      real (kind=8) :: sfilt(NTRU+1)
      
      EZ     = 1.63299310207D0
      PI     = 3.14159265359D0
      TWOPI  = PI * PI
      
      NLPP = NLAT         ! Latitudes per process
      NHOR = NLON * NLPP        ! Horizontal part
      NUGP = NLON * NLAT        ! Number of gridpoints
      NPGP = NLON * NLAT / 2    ! Dimension of packed fields
      NLEM = NLEV - 1           ! Levels - 1
      NLEP = NLEV + 1           ! Levels + 1
      NLSQ = NLEV * NLEV        ! Levels squared
      NTP1 = NTRU + 1           ! Truncation + 1
      NRSP =(NTRU+1)*(NTRU+2)   ! No of real global    modes
      NCSP = NRSP / 2           ! No of complex global modes
      NSPP = (NRSP+1-1)/1 ! Modes per process
      NESP = NSPP * 1        ! Dim of spectral fields
      NVCT = 2 * (NLEV+1)       ! Dim of Vert. Coord. Tab
      
      call legini(NLAT,NLON,NTRU,NLEV,qi,qj,qc,qe,qm,qq,qu,qv,&
     &            sfilt,nfs)
      
           
      sp(:,:) = 0.0
      
      if (NLPP < NLAT) then  ! Universal (parallel executable) version
      !----------------------------------------------------------------------
        do l = 1 , NLPP
          w = 1
          do m = 1 , NTP1
            do n = m , NTP1
              sp(1,w) = sp(1,w) + qc(w,l) * fc(1,m,l)*sfilt(n)
              sp(2,w) = sp(2,w) + qc(w,l) * fc(2,m,l)*sfilt(n)
              w = w + 1
            enddo ! n
          enddo ! m
        enddo ! l
      else                   ! Single CPU version (symmetry conserving)
      !----------------------------------------------------------------------
        do l = 1 , NLAT/2
          w = 1
          do m = 1 , NTP1
            do n = m , NTP1
              if (mod(m+n,2) == 0) then ! Symmetric modes
                sp(1,w) = sp(1,w) + qc(w,l) * (fc(1,m,l) + &
     &                                    fc(1,m,NLAT+1-l))*sfilt(n)
                sp(2,w) = sp(2,w) + qc(w,l) * (fc(2,m,l) + &
     &                                    fc(2,m,NLAT+1-l))*sfilt(n)
              else                      ! Antisymmetric modes
                sp(1,w) = sp(1,w) + qc(w,l) * (fc(1,m,l) - &
     &                                    fc(1,m,NLAT+1-l))*sfilt(n)
                sp(2,w) = sp(2,w) + qc(w,l) * (fc(2,m,l) - &
     &                                    fc(2,m,NLAT+1-l))*sfilt(n)
              endif
              w = w + 1
            enddo ! n
          enddo ! m
        enddo ! l
      !----------------------------------------------------------------------
      endif ! parallel ?
      return
      end
      
      
      ! ================
      ! SUBROUTINE SP2FC
      ! ================
      
      subroutine sp2fc(sp,fc,NLAT,NLON,NTRU,NLEV,nfs) ! Spectral to Fourier
      implicit none
      
      integer, intent(in ) :: NLAT
      integer, intent(in ) :: NLON
      integer, intent(in ) :: NLEV
      integer, intent(in ) :: NTRU
      integer, intent(in ) :: nfs
      
      real (kind=8), intent(in ) :: sp(2,(NTRU+1)*(NTRU+2)/2) ! Coefficients of spherical harmonics
!f2py intent(in) :: sp
      real (kind=8), intent(out) :: fc(2,NLON/2,NLAT) ! Fourier coefficients
!f2py intent(out) :: fc
      
      integer :: l ! Loop index for latitude
      integer :: m ! Loop index for zonal wavenumber m
      integer :: n ! Loop index for total wavenumber n
      integer :: w ! Loop index for spectral mode
      
       
      real (kind=8):: qi((NTRU+1)*(NTRU+2)/2,NLAT) ! P(m,n) = Associated Legendre Polynomials
      real (kind=8):: qj((NTRU+1)*(NTRU+2)/2,NLAT) ! Q(m,n) = Used for d/d(mu)
      real (kind=8):: qc((NTRU+1)*(NTRU+2)/2,NLAT) ! P(m,n) * gwd              used in fc2sp
      real (kind=8):: qe((NTRU+1)*(NTRU+2)/2,NLAT) ! Q(mn,) * gwd / cos2       used in mktend
      real (kind=8):: qm((NTRU+1)*(NTRU+2)/2,NLAT) ! P(m,n) * gwd / cos2 * m   used in mktend
      real (kind=8):: qq((NTRU+1)*(NTRU+2)/2,NLAT) ! P(m,n) * gwd / cos2 * n * (n+1) / 2  "
      real (kind=8):: qu((NTRU+1)*(NTRU+2)/2,NLAT) ! P(m,n) / (n*(n+1)) * m    used in dv2uv
      real (kind=8):: qv((NTRU+1)*(NTRU+2)/2,NLAT) ! Q(m,n) / (n*(n+1))        used in dv2uv
                      
      integer NLPP, NHOR, NUGP, NPGP, NLEM, NLEP, NLSQ, NTP1
      integer NRSP, NCSP, NSPP, NESP, NVCT
      
      real (kind=8) :: EZ
      real (kind=8) :: PI
      real (kind=8) :: TWOPI
      
      real (kind=8) :: sfilt(NTRU+1)
      
      EZ     = 1.63299310207D0
      PI     = 3.14159265359D0
      TWOPI  = PI * PI
      
      NLPP = NLAT         ! Latitudes per process
      NHOR = NLON * NLPP        ! Horizontal part
      NUGP = NLON * NLAT        ! Number of gridpoints
      NPGP = NLON * NLAT / 2    ! Dimension of packed fields
      NLEM = NLEV - 1           ! Levels - 1
      NLEP = NLEV + 1           ! Levels + 1
      NLSQ = NLEV * NLEV        ! Levels squared
      NTP1 = NTRU + 1           ! Truncation + 1
      NRSP =(NTRU+1)*(NTRU+2)   ! No of real global    modes
      NCSP = NRSP / 2           ! No of complex global modes
      NSPP = (NRSP+1-1)/1 ! Modes per process
      NESP = NSPP * 1        ! Dim of spectral fields
      NVCT = 2 * (NLEV+1)       ! Dim of Vert. Coord. Tab
      
      call legini(NLAT,NLON,NTRU,NLEV,qi,qj,qc,qe,qm,qq,qu,qv, &
     &            sfilt,nfs)
      
      fc(:,:,:) = 0.0
      
      do l = 1 , NLPP
         w = 1  
         do m = 1 , NTP1
            do n = m , NTP1
               fc(1,m,l) = fc(1,m,l) + qi(w,l) * sp(1,w)*sfilt(n) !72
               fc(2,m,l) = fc(2,m,l) + qi(w,l) * sp(2,w)*sfilt(n)
               w = w + 1
            enddo ! n
         enddo ! m
      enddo ! l
      return
      end      
      
      
      ! ================
      ! SUBROUTINE SP2FC
      ! ================
      
      
      subroutine sp3fc(spp,fc,NLAT,NLON,NTRU,NLEV,nfs)
      implicit none
      integer :: v, k ! Loop index for level
      integer, intent(in) :: NLEV
      integer, intent(in) :: NLON
      integer, intent(in) :: NTRU
      integer, intent(in) :: NLAT
      integer, intent(in) :: nfs
      
      real (kind=8), intent(in ) :: spp((NTRU+1)*(NTRU+2), NLEV) 
!f2py intent(in ) :: spp
      real (kind=8):: sp(2,(NTRU+1)*(NTRU+2)/2)
      real (kind=8):: fcc(2,NLON/2,NLAT)
      real (kind=8), intent(out) :: fc(2,NLON/2,NLAT, NLEV) ! Fourier coefficients
!f2py intent(out) :: fc
      
      do v = 1 , NLEV
         do k = 1, (NTRU+1)*(NTRU+2)/2
           sp(1,k) = spp(2*k-1,v)
           sp(2,k) = spp(2*k  ,v)
         enddo
         call sp2fc(sp,fcc,NLAT,NLON,NTRU,NLEV,nfs)
         fc(:,:,:,v) = fcc(:,:,:)
      enddo
      return
      end      
            

      ! ===================
      ! SUBROUTINE SP2FCDMU
      ! ===================
      
      subroutine sp2fcdmu(sp,fc,NLAT,NLON,NTRU,NLEV,nfs) ! Spectral to Fourier d/dmu
      implicit none
      
      integer, intent(in) :: NLEV
      integer, intent(in) :: NLON
      integer, intent(in) :: NTRU
      integer, intent(in) :: NLAT
      integer, intent(in) :: nfs
      
      real(kind=8), intent(in) :: sp(2,(NTRU+1)*(NTRU+2)/2)! Coefficients of spherical harmonics
!f2py intent(in ) :: sp
      real(kind=8), intent(out):: fc(2,NLON/2,NLAT) ! Fourier coefficients
!f2py intent(out) :: fc      
      
      integer :: l ! Loop index for latitude
      integer :: m ! Loop index for zonal wavenumber m
      integer :: n ! Loop index for total wavenumber n
      integer :: w ! Loop index for spectral mode
        
      real (kind=8):: qi((NTRU+1)*(NTRU+2)/2,NLAT) ! P(m,n) = Associated Legendre Polynomials
      real (kind=8):: qj((NTRU+1)*(NTRU+2)/2,NLAT) ! Q(m,n) = Used for d/d(mu)
      real (kind=8):: qc((NTRU+1)*(NTRU+2)/2,NLAT) ! P(m,n) * gwd              used in fc2sp
      real (kind=8):: qe((NTRU+1)*(NTRU+2)/2,NLAT) ! Q(mn,) * gwd / cos2       used in mktend
      real (kind=8):: qm((NTRU+1)*(NTRU+2)/2,NLAT) ! P(m,n) * gwd / cos2 * m   used in mktend
      real (kind=8):: qq((NTRU+1)*(NTRU+2)/2,NLAT) ! P(m,n) * gwd / cos2 * n * (n+1) / 2  "
      real (kind=8):: qu((NTRU+1)*(NTRU+2)/2,NLAT) ! P(m,n) / (n*(n+1)) * m    used in dv2uv
      real (kind=8):: qv((NTRU+1)*(NTRU+2)/2,NLAT) ! Q(m,n) / (n*(n+1))        used in dv2uv
                      
      integer NLPP, NHOR, NUGP, NPGP, NLEM, NLEP, NLSQ, NTP1
      integer NRSP, NCSP, NSPP, NESP, NVCT
      
      real (kind=8) :: EZ
      real (kind=8) :: PI
      real (kind=8) :: TWOPI
      
      real (kind=8) :: sfilt(NTRU+1)
      
      EZ     = 1.63299310207D0
      PI     = 3.14159265359D0
      TWOPI  = PI * PI
      
      NLPP = NLAT         ! Latitudes per process
      NHOR = NLON * NLPP        ! Horizontal part
      NUGP = NLON * NLAT        ! Number of gridpoints
      NPGP = NLON * NLAT / 2    ! Dimension of packed fields
      NLEM = NLEV - 1           ! Levels - 1
      NLEP = NLEV + 1           ! Levels + 1
      NLSQ = NLEV * NLEV        ! Levels squared
      NTP1 = NTRU + 1           ! Truncation + 1
      NRSP =(NTRU+1)*(NTRU+2)   ! No of real global    modes
      NCSP = NRSP / 2           ! No of complex global modes
      NSPP = (NRSP+1-1)/1 ! Modes per process
      NESP = NSPP * 1        ! Dim of spectral fields
      NVCT = 2 * (NLEV+1)       ! Dim of Vert. Coord. Tab
      
      call legini(NLAT,NLON,NTRU,NLEV,qi,qj,qc,qe,qm,qq,qu,qv, &
     &            sfilt,nfs)
     
      fc(:,:,:) = 0.0
      
      do l = 1 , NLPP
         w = 1  
         do m = 1 , NTP1
            do n = m , NTP1
               fc(1,m,l) = fc(1,m,l) + qj(w,l) * sp(1,w)*sfilt(n)
               fc(2,m,l) = fc(2,m,l) + qj(w,l) * sp(2,w)*sfilt(n)
               w = w + 1
            enddo ! n
         enddo ! m
      enddo ! l
      return
      end


      ! ================
      ! SUBROUTINE DV2UV        !SP->GP
      ! ================
      
      subroutine dv2uv(sd,sz,pu,pv,NLAT,NLON,NTRU,NLEV,nfs)
      implicit none
      
      integer, intent(in) :: NLEV
      integer, intent(in) :: NLON
      integer, intent(in) :: NTRU
      integer, intent(in) :: NLAT
      integer, intent(in) :: nfs
      
      real(kind=8), intent(in)  :: sd((NTRU+1)*(NTRU+2),NLEV)
!f2py intent(in) :: sd
      real(kind=8), intent(in)  :: sz((NTRU+1)*(NTRU+2),NLEV)
!f2py intent(in) :: sz
      real(kind=8)  :: pd(2,(NTRU+1)*(NTRU+2)/2,NLEV)
      real(kind=8)  :: pz(2,(NTRU+1)*(NTRU+2)/2,NLEV)
      real(kind=8), intent(out) :: pu(2,NLON/2,NLAT,NLEV)
!f2py intent(out) :: pu
      real(kind=8), intent(out) :: pv(2,NLON/2,NLAT,NLEV)
!f2py intent(out) :: pv
      integer :: l ! Loop index for latitude
      integer :: m ! Loop index for zonal wavenumber m
      integer :: n ! Loop index for total wavenumber n
      integer :: v ! Loop index for level
      integer :: w ! Loop index for spectral mode
      integer :: k
         
      real (kind=8):: qi((NTRU+1)*(NTRU+2)/2,NLAT) ! P(m,n) = Associated Legendre Polynomials
      real (kind=8):: qj((NTRU+1)*(NTRU+2)/2,NLAT) ! Q(m,n) = Used for d/d(mu)
      real (kind=8):: qc((NTRU+1)*(NTRU+2)/2,NLAT) ! P(m,n) * gwd              used in fc2sp
      real (kind=8):: qe((NTRU+1)*(NTRU+2)/2,NLAT) ! Q(mn,) * gwd / cos2       used in mktend
      real (kind=8):: qm((NTRU+1)*(NTRU+2)/2,NLAT) ! P(m,n) * gwd / cos2 * m   used in mktend
      real (kind=8):: qq((NTRU+1)*(NTRU+2)/2,NLAT) ! P(m,n) * gwd / cos2 * n * (n+1) / 2  "
      real (kind=8):: qu((NTRU+1)*(NTRU+2)/2,NLAT) ! P(m,n) / (n*(n+1)) * m    used in dv2uv
      real (kind=8):: qv((NTRU+1)*(NTRU+2)/2,NLAT) ! Q(m,n) / (n*(n+1))        used in dv2uv
                      
      integer NLPP, NHOR, NUGP, NPGP, NLEM, NLEP, NLSQ, NTP1
      integer NRSP, NCSP, NSPP, NESP, NVCT
      
      real (kind=8) :: EZ
      real (kind=8) :: PI
      real (kind=8) :: TWOPI
      
      real (kind=8) :: sfilt(NTRU+1)
      
      EZ     = 1.63299310207D0
      PI     = 3.14159265359D0
      TWOPI  = PI * PI
      
      NLPP = NLAT         ! Latitudes per process
      NHOR = NLON * NLPP        ! Horizontal part
      NUGP = NLON * NLAT        ! Number of gridpoints
      NPGP = NLON * NLAT / 2    ! Dimension of packed fields
      NLEM = NLEV - 1           ! Levels - 1
      NLEP = NLEV + 1           ! Levels + 1
      NLSQ = NLEV * NLEV        ! Levels squared
      NTP1 = NTRU + 1           ! Truncation + 1
      NRSP =(NTRU+1)*(NTRU+2)   ! No of real global    modes
      NCSP = NRSP / 2           ! No of complex global modes
      NSPP = (NRSP+1-1)/1 ! Modes per process
      NESP = NSPP * 1        ! Dim of spectral fields
      NVCT = 2 * (NLEV+1)       ! Dim of Vert. Coord. Tab
      
      call legini(NLAT,NLON,NTRU,NLEV,qi,qj,qc,qe,qm,qq,qu,qv, &
     &            sfilt,nfs)
     
      do v=1,NLEV
        do k=1,(NTRU+1)*(NTRU+2)/2
          pd(1,k,v) = sd(2*k-1,v)
          pd(2,k,v) = sd(2*K  ,v)
          pz(1,k,v) = sz(2*k-1,v)
          pz(2,k,v) = sz(2*k  ,v)
        enddo
      enddo
     
      pu(:,:,:,:) = 0.0
      pv(:,:,:,:) = 0.0
      
      do v = 1 , NLEV
!         zsave = pz(1,2,v)
!         pz(1,2,v) = zsave - plavor
        do l = 1 , NLPP
          w = 1
          do m = 1 , NTP1
            do n = m , NTP1
              pu(1,m,l,v)=pu(1,m,l,v)+qv(w,l)*pz(1,w,v)*sfilt(n) &
     &                               +qu(w,l)*pd(2,w,v)*sfilt(n)
              pu(2,m,l,v)=pu(2,m,l,v)+qv(w,l)*pz(2,w,v)*sfilt(n) &
     &                               -qu(w,l)*pd(1,w,v)*sfilt(n)
              pv(1,m,l,v)=pv(1,m,l,v)+qu(w,l)*pz(2,w,v)*sfilt(n) &
     &                               -qv(w,l)*pd(1,w,v)*sfilt(n)
              pv(2,m,l,v)=pv(2,m,l,v)-qu(w,l)*pz(1,w,v)*sfilt(n) &
     &                               -qv(w,l)*pd(2,w,v)*sfilt(n)
              w = w + 1
            enddo ! n
          enddo ! m
        enddo ! l
!         pz(1,2,v) = zsave
      enddo ! jv
      return
      end


      ! ================
      ! SUBROUTINE UV2DV        !GP->SP
      ! ================
      
      subroutine uv2dv(pu,pv,pd,pz,NLAT,NLON,NTRU,NLEV,nfs)
      implicit none
      
      integer, intent(in) :: NLEV
      integer, intent(in) :: NLON
      integer, intent(in) :: NTRU
      integer, intent(in) :: NLAT
      integer, intent(in) :: nfs
      
      real(kind=8), intent(out) :: pd(2,(NTRU+1)*(NTRU+2)/2,NLEV)
!f2py intent(out) :: pd
      real(kind=8), intent(out) :: pz(2,(NTRU+1)*(NTRU+2)/2,NLEV)
!f2py intent(out) :: pz
      real(kind=8), intent(in) :: pu(2,NLON/2,NLAT,NLEV)
!f2py intent(in) :: pu
      real(kind=8), intent(in) :: pv(2,NLON/2,NLAT,NLEV)
!f2py intent(in) :: pv
      
      integer :: k ! Loop index for southern latitude
      integer :: l ! Loop index for latitude
      integer :: m ! Loop index for zonal wavenumber m
      integer :: n ! Loop index for total wavenumber n
      integer :: v ! Loop index for level
      integer :: w ! Loop index for spectral mode
         
      real (kind=8):: qi((NTRU+1)*(NTRU+2)/2,NLAT) ! P(m,n) = Associated Legendre Polynomials
      real (kind=8):: qj((NTRU+1)*(NTRU+2)/2,NLAT) ! Q(m,n) = Used for d/d(mu)
      real (kind=8):: qc((NTRU+1)*(NTRU+2)/2,NLAT) ! P(m,n) * gwd              used in fc2sp
      real (kind=8):: qe((NTRU+1)*(NTRU+2)/2,NLAT) ! Q(mn,) * gwd / cos2       used in mktend
      real (kind=8):: qm((NTRU+1)*(NTRU+2)/2,NLAT) ! P(m,n) * gwd / cos2 * m   used in mktend
      real (kind=8):: qq((NTRU+1)*(NTRU+2)/2,NLAT) ! P(m,n) * gwd / cos2 * n * (n+1) / 2  "
      real (kind=8):: qu((NTRU+1)*(NTRU+2)/2,NLAT) ! P(m,n) / (n*(n+1)) * m    used in dv2uv
      real (kind=8):: qv((NTRU+1)*(NTRU+2)/2,NLAT) ! Q(m,n) / (n*(n+1))        used in dv2uv
                      
      integer NLPP, NHOR, NUGP, NPGP, NLEM, NLEP, NLSQ, NTP1
      integer NRSP, NCSP, NSPP, NESP, NVCT
      
      real (kind=8) :: EZ
      real (kind=8) :: PI
      real (kind=8) :: TWOPI
      
      real (kind=8) :: sfilt(NTRU+1)
      
      EZ     = 1.63299310207D0
      PI     = 3.14159265359D0
      TWOPI  = PI * PI
      
      NLPP = NLAT         ! Latitudes per process
      NHOR = NLON * NLPP        ! Horizontal part
      NUGP = NLON * NLAT        ! Number of gridpoints
      NPGP = NLON * NLAT / 2    ! Dimension of packed fields
      NLEM = NLEV - 1           ! Levels - 1
      NLEP = NLEV + 1           ! Levels + 1
      NLSQ = NLEV * NLEV        ! Levels squared
      NTP1 = NTRU + 1           ! Truncation + 1
      NRSP =(NTRU+1)*(NTRU+2)   ! No of real global    modes
      NCSP = NRSP / 2           ! No of complex global modes
      NSPP = (NRSP+1-1)/1 ! Modes per process
      NESP = NSPP * 1        ! Dim of spectral fields
      NVCT = 2 * (NLEV+1)       ! Dim of Vert. Coord. Tab
      
      call legini(NLAT,NLON,NTRU,NLEV,qi,qj,qc,qe,qm,qq,qu,qv, &
     &            sfilt,nfs)
     
      pd(:,:,:) = 0.0
      pz(:,:,:) = 0.0
      
      if (NLPP < NLAT) then  ! Universal (parallel executable) version
      !----------------------------------------------------------------------
      do v = 1 , NLEV
        do l = 1 , NLPP
          w = 1
          do m = 1 , NTP1
            do n = m , NTP1
              pz(1,w,v) = pz(1,w,v)+qe(w,l)*pu(1,m,l,v)*sfilt(n)&
     &                             -qm(w,l)*pv(2,m,l,v)*sfilt(n)
              pz(2,w,v) = pz(2,w,v)+qe(w,l)*pu(2,m,l,v)*sfilt(n)&
     &                             +qm(w,l)*pv(1,m,l,v)*sfilt(n)
              pd(1,w,v) = pd(1,w,v)-qe(w,l)*pv(1,m,l,v)*sfilt(n)&
     &                             -qm(w,l)*pu(2,m,l,v)*sfilt(n)
              pd(2,w,v) = pd(2,w,v)-qe(w,l)*pv(2,m,l,v)*sfilt(n)&
     &                             +qm(w,l)*pu(1,m,l,v)*sfilt(n)
              w = w + 1
            enddo ! n
          enddo ! m
        enddo ! l
      enddo ! v
      else                   ! Single CPU version (symmetry conserving)
      !----------------------------------------------------------------------
      do v = 1 , NLEV
        do l = 1 , NLAT/2
          k = NLAT+1-l
          w = 1
          do m = 1 , NTP1
            do n = m , NTP1
              if (mod(m+n,2) == 0) then ! symmetric -----------------
                pz(1,w,v) = pz(1,w,v) + qe(w,l) * &
     &                           (pu(1,m,l,v)-pu(1,m,k,v))*sfilt(n) &
     &               - qm(w,l) * (pv(2,m,l,v)+pv(2,m,k,v))*sfilt(n)
                pz(2,w,v) = pz(2,w,v) + qe(w,l) *  &
     &                           (pu(2,m,l,v)-pu(2,m,k,v))*sfilt(n) &
     &               + qm(w,l) * (pv(1,m,l,v)+pv(1,m,k,v))*sfilt(n)
                pd(1,w,v) = pd(1,w,v) - qe(w,l) *  &
     &                           (pv(1,m,l,v)-pv(1,m,k,v))*sfilt(n) &
     &               - qm(w,l) * (pu(2,m,l,v)+pu(2,m,k,v))*sfilt(n)
                pd(2,w,v) = pd(2,w,v) - qe(w,l) *  &
     &                           (pv(2,m,l,v)-pv(2,m,k,v))*sfilt(n) &
     &               + qm(w,l) * (pu(1,m,l,v)+pu(1,m,k,v))*sfilt(n)
              else ! ---------------- antisymmetric -----------------
                pz(1,w,v) = pz(1,w,v) + qe(w,l) *  &
     &                           (pu(1,m,l,v)+pu(1,m,k,v))*sfilt(n) &
     &               - qm(w,l) * (pv(2,m,l,v)-pv(2,m,k,v))*sfilt(n)
                pz(2,w,v) = pz(2,w,v) + qe(w,l) *  &
     &                           (pu(2,m,l,v)+pu(2,m,k,v))*sfilt(n) &
     &               + qm(w,l) * (pv(1,m,l,v)-pv(1,m,k,v))*sfilt(n)
                pd(1,w,v) = pd(1,w,v) - qe(w,l) *  &
     &                           (pv(1,m,l,v)+pv(1,m,k,v))*sfilt(n) &
     &               - qm(w,l) * (pu(2,m,l,v)-pu(2,m,k,v))*sfilt(n)
                pd(2,w,v) = pd(2,w,v) - qe(w,l) *  &
     &                           (pv(2,m,l,v)+pv(2,m,k,v))*sfilt(n) &
     &               + qm(w,l) * (pu(1,m,l,v)-pu(1,m,k,v))*sfilt(n)
              endif
              w = w + 1
            enddo ! n
          enddo ! m
        enddo ! l
      enddo ! v
      !----------------------------------------------------------------------
      endif ! symmetric?
      return
      end 

            
      ! =================
      ! SUBROUTINE INIGAU
      ! =================
      
      subroutine inigau(klat,pz0,pzw)        ! pz0 & pzw are (kind=8) reals !!!
      implicit none
      integer ,intent(IN)        :: klat       ! Number of Gaussian latitudes
      real (kind=8), intent(out) :: pz0(klat)  ! Gaussian abscissas
      real (kind=8), intent(out) :: pzw(klat)  ! Gaussian weights
      integer                  :: jlat       ! Latitudinal loop index
      integer                  :: jiter      ! Iteration loop index
      integer      , parameter :: NITER = 50 ! Maximum # of iterations
      real (kind=8), parameter :: PI    =  3.14159265358979_8
      real (kind=8), parameter :: ZEPS  =  1.0e-16 ! Convergence criterion
      real (kind=8) :: z0,z1,z2,z3,z4,z5
      real (kind=8) :: ql,qld
      
      ! Compute Gaussian abscissas & weights
      
      z0 = PI / (2*klat+1)
      z1 = 1.0_8 / (klat*klat*8)
      z4 = 2.0_8 / (klat*klat)
      
      do jlat = 1 , klat/2
         z2 = z0 * (2*jlat - 0.5_8)
         z2 = cos(z2 + z1 / tan(z2))
         do jiter = 1 , NITER
            z3 = ql(klat,z2) * qld(klat,z2)
            z2 = z2 - z3
            if (abs(z3) < ZEPS) exit ! converged
         enddo ! jiter
         z5 = ql(klat-1,z2) / sqrt(klat - 0.5_8)
         pz0(jlat) = z2
         pzw(jlat) = z4 * (1.0_8 - z2 * z2) / (z5 * z5)
         pz0(klat-jlat+1) = -z2
         pzw(klat-jlat+1) = pzw(jlat)
      enddo ! jlat
      
      return
      end subroutine inigau
      
      ! ===========
      ! FUNCTION QL
      ! ===========
      
      real (kind=8) function ql(k,p)
      implicit none
      integer      , intent(IN) :: k
      real (kind=8), intent(IN) :: p
      real (kind=8) :: z0,z1,z2,z3,z4
      integer :: j
      z0 = acos(p)
      z1 = 1.0
      z2 = 0.0
      do j = k , 0 , -2
         z3 = z1 * cos(z0 * j)
         z2 = z2 + z3
         z4 = (k-j+1) * (k+j) * 0.5_8
         z1 = z1 * z4 / (z4 + (j-1))
      enddo ! j
      if (mod(k,2) == 0) z2 = z2 - 0.5_8 * z3
      
      z0 = sqrt(2.0_8)
      do j = 1 ,k
         z0 = z0 * sqrt(1.0_8 - 0.25_8 / (j*j))
      enddo ! j
      ql = z0 * z2
      return
      end function ql
      
      ! ============
      ! FUNCTION QLD
      ! ============
      
      real (kind=8) function qld(k,p)
      implicit none
      integer      , intent(IN) :: k
      real (kind=8), intent(IN) :: p
      real (kind=8) :: z
      real (kind=8) :: ql
      
      z = p * ql(k,p) - sqrt((k + k + 1.0_8) / &
     &            (k + k - 1.0_8)) * ql(k-1,p)
      qld = (p * p - 1.0_8) / (k * z)
      
      return
      end function qld
       
       



!     =============
!     MODULE FFTMOD
!     =============

      module fftmod
      parameter(NRES = 13)
      integer :: nallowed(NRES)=(/8,16,32,48,64,96,128,256,384,512, &
     &                           1024,2048,4096/)
!     T5    - N16   : 8-2
!     T10   - N32   : 8-2-2
!     T15   - N48   : 8-3-2
!     T21   - N64   : 8-4-2
!     T31   - N96   : 8-4-3
!     T42   - N128  : 8-4-4
!     T85   - N256  : 8-4-4-2
!     T127  - N384  : 8-4-4-3
!     T170  - N512  : 8-4-4-4
!     T341  - N1024 : 8-4-4-4-2
!     T682  - N2048 : 8-4-4-4-4
!     T1365 - N4096 : 8-4-4-4-4-2

      integer :: lastn = 0
      real (kind=8),allocatable :: trigs(:)
      end module fftmod

!     ================
!     SUBROUTINE GP2FC
!     ================

      subroutine gp2fc(a,c,n,lot)
      use fftmod
      real (kind=8), intent(in) :: a
      real (kind=8), intent(out):: c
      integer, intent(in) :: n
      integer, intent(in) :: lot
      
      dimension a(n,lot)
      dimension c(n,lot)

      if (n /= lastn) then
         if (allocated(trigs)) deallocate(trigs)
         allocate(trigs(n))
         lastn = n
         call fftini(n)
      endif

      call dfft8(a,c,n,lot)
      la = n / 8
      do while (la >= 4)
         call dfft4(c,trigs,n,lot,la)
      enddo

      if (la == 3) then
         do l = 1 , lot
            call dfft3(c(1,l),trigs,n)
         enddo
      endif

      if (la == 2) then
         do l = 1 , lot
            call dfft2(c(1,l),trigs,n)
         enddo
      endif
      return
      end subroutine gp2fc

!     ================
!     SUBROUTINE FC2GP
!     ================

      subroutine fc2gp(a,c,n,lot)
      use fftmod
      real (kind=8), intent(in) :: a
!f2py intent(in) :: a
      real (kind=8), intent(out) :: c
!f2py intent(out) :: c
      integer, intent(in) :: n
!f2py intent(in) :: n
      integer, intent(in) :: lot
!f2py intent(in) :: lot
      real (kind=8) b
      
      dimension a(n,lot)
      dimension b(n,lot)
      dimension c(n,lot)
      
!       write(*,*) n
!       write(*,*) lot
      
      if (n /= lastn) then
         if (allocated(trigs)) deallocate(trigs)
         allocate(trigs(n))
         lastn = n
         call fftini(n)
      endif

      b(:,:) = a(:,:)
      
      nf = n/8
      do while (nf >= 4)
         nf = nf/4
      enddo
      la = 1
      if (nf == 2) call ifft2(b,trigs,n,lot,la)
      if (nf == 3) call ifft3(b,trigs,n,lot,la)
      do while (la < n/8)
         call ifft4(b,trigs,n,lot,la)
      enddo
      call ifft8(b,c,n,lot)
      
      c(:,:) = c(:,:)*1.4142135623730951 ! *sqrt(2)
      
      return
      end subroutine fc2gp
      
!     ================
!     SUBROUTINE FC3GP
!     ================
      
      
      subroutine fc3gp(a,c,n,lot,NLEV)
      implicit none
      real (kind=8), intent(in) :: a
!f2py intent(in) :: a
      real (kind=8), intent(out) :: c
!f2py intent(out) :: c
      integer, intent(in) :: n
!f2py intent(in) :: n
      integer, intent(in) :: lot
!f2py intent(in) :: lot
      integer, intent(in) :: NLEV
!f2py intent(in) :: NLEV

      integer jlev
      real (kind=8) :: dd
      
      dimension dd(n,lot)
      dimension a(n,lot,NLEV)
      dimension c(n,lot,NLEV)
      
      do jlev=1,NLEV
        call fc2gp(a(:,:,jlev),dd,n,lot)
        c(:,:,jlev) = dd(:,:)
      enddo
      
      return
      end subroutine
      
!     ================
!     SUBROUTINE GP3FC
!     ================
      
      
      subroutine gp3fc(a,c,n,lot,NLEV)
      implicit none
      real (kind=8), intent(in) :: a
!f2py intent(in) :: a
      real (kind=8), intent(out) :: c
!f2py intent(out) :: c
      integer, intent(in) :: n
!f2py intent(in) :: n
      integer, intent(in) :: lot
!f2py intent(in) :: lot
      integer, intent(in) :: NLEV
!f2py intent(in) :: NLEV

      integer jlev
      real (kind=8) :: dd
      
      dimension dd(n,lot)
      dimension a(n,lot,NLEV)
      dimension c(n,lot,NLEV)
      
      do jlev=1,NLEV
        call gp2fc(a(:,:,jlev),dd,n,lot)
        c(:,:,jlev) = dd(:,:)
      enddo
      
      return
      end subroutine
      
      
!     ================
!     SUBROUTINE SP2GP
!     ================
      
      subroutine sp2gp(sp,gp,NLAT,NLON,NTRU,NLEV,nfs)
      implicit none
      integer :: v, k ! Loop index for level
      integer, intent(in) :: NLEV
      integer, intent(in) :: NLON
      integer, intent(in) :: NTRU
      integer, intent(in) :: NLAT
      integer, intent(in) :: nfs !1 or 0, whether use physics filter
      
      real (kind=8), intent(in ) :: sp((NTRU+1)*(NTRU+2), NLEV) 
!f2py intent(in ) :: sp
      real (kind=8):: spp(2,(NTRU+1)*(NTRU+2)/2)
      real (kind=8):: fcc(2,NLON/2,NLAT)
      real (kind=8):: fcl(NLON,NLAT)
      real (kind=8):: gpp(NLON,NLAT)
      real (kind=8), intent(out) :: gp(NLON, NLAT, NLEV) ! Fourier coefficients
!f2py intent(out) :: gp
      
      do v = 1 , NLEV
         do k = 1, (NTRU+1)*(NTRU+2)/2
           spp(1,k) = sp(2*k-1,v)
           spp(2,k) = sp(2*k  ,v)
         enddo
         call sp2fc(spp,fcc,NLAT,NLON,NTRU,NLEV,nfs)
         do k = 1, NLON/2
           fcl(2*k-1,:) = fcc(1,k,:)
           fcl(2*k  ,:) = fcc(2,k,:)
         enddo
         call fc2gp(fcl,gpp,NLON,NLAT)
         gp(:,:,v) = gpp(:,:)
!          write(*,*) "Computed layer",v
      enddo
      
      return
      end subroutine
      
      
!     ================
!     SUBROUTINE SPVGP
!     ================
      
      subroutine spvgp(sd,sz,rdcostheta,gu,gv,NLAT,NLON,NTRU,NLEV,nfs)
      implicit none
      integer :: v, k, l ! Loop index for level
      integer, intent(in ) :: NLEV
      integer, intent(in ) :: NLON
      integer, intent(in ) :: NTRU
      integer, intent(in ) :: NLAT
      integer, intent(in ) :: nfs
      real (kind=8), intent(in ) :: sd((NTRU+1)*(NTRU+2),NLEV)
!f2py intent(in ) :: sd
      real (kind=8), intent(in ) :: sz((NTRU+1)*(NTRU+2),NLEV)
!f2py intent(in ) :: sz
      real (kind=8), intent(in ) :: rdcostheta(NLAT) !radius / cos(latitude)
!f2py intent(in ) :: rdcostheta
      real (kind=8), intent(out) :: gu(NLON, NLAT, NLEV)
!f2py intent(out) :: gu
      real (kind=8), intent(out) :: gv(NLON, NLAT, NLEV)
!f2py intent(out) :: gv

!       real (kind=8) sdd(2, (NTRU+1)*(NTRU+2)/2, NLEV)
!       real (kind=8) szz(2, (NTRU+1)*(NTRU+2)/2, NLEV)
      real (kind=8) fuu(2, NLON/2, NLAT, NLEV)
      real (kind=8) fvv(2, NLON/2, NLAT, nLEV)
      real (kind=8) fup(NLON, NLAT)
      real (kind=8) fvp(NLON, NLAT)
      real (kind=8) gup(NLON, NLAT)
      real (kind=8) gvp(NLON, NLAT)
      
!       do v = 1, NLEV
!         do k = 1, (NTRU+1)*(NTRU+2)/2
!           sdd(1,k,v) = sd(2*k-1,v)
!           sdd(2,k,v) = sd(2*k  ,v)
!           szz(1,k,v) = sz(2*k-1,v)
!           szz(2,k,v) = sz(2*k  ,v)
!         enddo
!       enddo
      
      call dv2uv(sd,sz,fuu,fvv,NLAT,NLON,NTRU,NLEV,nfs)
      
      do v = 1, NLEV
        do l = 1, NLAT
          do k = 1, NLON/2
            fup(2*k-1,l) = fuu(1,k,l,v)
            fup(2*k  ,l) = fuu(2,k,l,v)
            fvp(2*k-1,l) = fvv(1,k,l,v)
            fvp(2*k  ,l) = fvv(2,k,l,v)
          enddo
        enddo
        call fc2gp(fup,gup,NLON,NLAT)
        call fc2gp(fvp,gvp,NLON,NLAT)
        do l = 1, NLAT
          gu(:,l,v) = gup(:,l)*rdcostheta(l)
          gv(:,l,v) = gvp(:,l)*rdcostheta(l)
        enddo
      enddo
      
      return
      end subroutine
      
      
!     ================
!     SUBROUTINE GPVSP
!     ================
      
      subroutine gpvsp(gu,gv,costhetadr,sd,sz,NLAT,NLON,NTRU,NLEV,nfs)
      implicit none
      integer :: v, k, l ! Loop index for level
      integer, intent(in ) :: NLEV
      integer, intent(in ) :: NLON
      integer, intent(in ) :: NTRU
      integer, intent(in ) :: NLAT
      integer, intent(in ) :: nfs  
      real (kind=8), intent(out) :: sd((NTRU+1)*(NTRU+2),NLEV)
!f2py intent(out) :: sd
      real (kind=8), intent(out) :: sz((NTRU+1)*(NTRU+2),NLEV)
!f2py intent(out) :: sz
      real (kind=8), intent(in ) :: gu(NLON, NLAT, NLEV)
!f2py intent(in ) :: gu
      real (kind=8), intent(in ) :: gv(NLON, NLAT, NLEV)
!f2py intent(in ) :: gv
      real (kind=8), intent(in ) :: costhetadr(NLAT) !cos(latitude) / radius
!f2py intent(in ) :: costhetadr    

      real (kind=8) sdd(2, (NTRU+1)*(NTRU+2)/2, NLEV)
      real (kind=8) szz(2, (NTRU+1)*(NTRU+2)/2, NLEV)
      real (kind=8) guu(2, NLON/2, NLAT, NLEV)
      real (kind=8) gvv(2, NLON/2, NLAT, nLEV)
      real (kind=8) gur(NLON, NLAT)
      real (kind=8) gvr(NLON, NLAT)
      real (kind=8) fuu(NLON, NLAT)
      real (kind=8) fvv(NLON, NLAT)
      
      do v = 1, NLEV
        do l = 1, NLAT
          gur(:,l) = gu(:,l,v)*costhetadr(l)/1.4142135623730951
          gvr(:,l) = gv(:,l,v)*costhetadr(l)/1.4142135623730951
        enddo
        call gp2fc(gur,fuu,NLON,NLAT)
        call gp2fc(gvr,fvv,NLON,NLAT)
        do l = 1, NLAT
          do k = 1, NLON/2
            guu(1,k,l,v) = fuu(2*k-1,l)
            guu(2,k,l,v) = fuu(2*k  ,l)
            gvv(1,k,l,v) = fvv(2*k-1,l)
            gvv(2,k,l,v) = fvv(2*k  ,l)
          enddo
        enddo
      enddo
      
      call uv2dv(guu,gvv,sdd,szz,NLAT,NLON,NTRU,NLEV,nfs)
      
      do v = 1, NLEV
        do k = 1, (NTRU+1)*(NTRU+2)/2
          sd(2*k-1,v) = sdd(1,k,v)
          sd(2*k  ,v) = sdd(2,k,v)
          sz(2*k-1,v) = szz(1,k,v)
          sz(2*k  ,v) = szz(2,k,v)
        enddo
      enddo
      
      return
      end subroutine
            

      
!     ================
!     SUBROUTINE GP2SP
!     ================
      
      subroutine gp2sp(gp,sp,NLAT,NLON,NTRU,NLEV,nfs)
      implicit none
      integer :: v, k ! Loop index for level
      integer, intent(in) :: NLEV
      integer, intent(in) :: NLON
      integer, intent(in) :: NTRU
      integer, intent(in) :: NLAT
      integer, intent(in) :: nfs !1 or 0, whether use physics filter
      
      real (kind=8), intent(out ) :: sp((NTRU+1)*(NTRU+2), NLEV) 
!f2py intent(out) :: sp
      real (kind=8):: spp(2,(NTRU+1)*(NTRU+2)/2)
      real (kind=8):: fcc(2,NLON/2,NLAT)
      real (kind=8):: fcl(NLON,NLAT)
      real (kind=8):: gpp(NLON,NLAT)
      real (kind=8), intent(in) :: gp(NLON, NLAT, NLEV) ! Gridpoint variable
!f2py intent(in ) :: gp
      
      do v = 1 , NLEV
         call gp2fc(gp(:,:,v)/1.4142135623730951,fcl,NLON,NLAT)
         do k = 1, NLON/2
            fcc(1,k,:) = fcl(2*k-1,:) 
            fcc(2,k,:) = fcl(2*k  ,:) 
         enddo
         call fc2sp(fcc,spp,NLAT,NLON,NTRU,NLEV,nfs)
         do k = 1, (NTRU+1)*(NTRU+2)/2
            sp(2*k-1,v) = spp(1,k) 
            sp(2*k  ,v) = spp(2,k) 
         enddo
!          write(*,*) "Computed layer",v
      enddo
      
      return
      end subroutine

!     =================
!     SUBROUTINE FFTINI
!     =================

      subroutine fftini(n)
      use fftmod
      integer, intent(in) :: n
      logical labort

!     check for allowed values of n

      labort = .true.
      do j = 1 , NRES
         if (n == nallowed(j)) labort = .false.
      enddo

      if (labort) then
         write (*,*) '*** FFT does not support n = ',n,' ***'
         write (*,*) 'Following resolutions may be used:'
         write (*,*) '----------------------------------'
         do j = 1 , NRES
            write (*,1000) nallowed(j), nallowed(j)/2, nallowed(j)/3
         enddo
         stop
      endif
 1000 format(' NLON=',I5,'  NLAT=',I5,'  NTRU=',I5)

      del = 4.0 * asin(1.0) / n
      do k=0,n/2-1
        angle = k * del
        trigs(2*k+1) = cos(angle)
        trigs(2*k+2) = sin(angle)
      enddo
      return
      end subroutine fftini

!     ================
!     SUBROUTINE DFFT2
!     ================

      subroutine dfft2(a,trigs,n)
      real (kind=8), intent(inout) :: a(n)
      real (kind=8) c(n)
      real (kind=8), intent(in) :: trigs(n)
      integer, intent(in) :: n

      c(1) = a(1) + a(2)
      c(2) = 0.0

      ja = 3
      jb = n - 1

      do i=3,n-5,4
         c1 = trigs(ja  )
         s1 = trigs(ja+1)
         a1p3 = c1 * a(i+1) + s1 * a(i+3)
         a3m1 = c1 * a(i+3) - s1 * a(i+1)
         c(ja  ) = a(i) + a1p3
         c(jb  ) = a(i) - a1p3
         c(ja+1) = a3m1 + a(i+2)
         c(jb+1) = a3m1 - a(i+2)
         ja = ja + 2
         jb = jb - 2
      enddo

      c(ja  ) =  a(n-1)
      c(ja+1) = -a(n  )

      a = c
      return
      end subroutine dfft2

!     ================
!     SUBROUTINE DFFT3
!     ================

      subroutine dfft3(a,trigs,n)
      real(kind=8) SIN60
      parameter(SIN60 = 0.866025403784438D0)
      integer, intent(in) :: n
      real(kind=8), intent(inout) :: a(n)
!f2py intent(in,out) :: a
      real(kind=8) c(n)
      real(kind=8), intent(in) :: trigs(n)

      ja = 1              !  1
      jb = 2 * (n/3)  + 1 ! 65
      jc = jb             ! 65

      c(ja  ) = a(1) + a(2) + a(3)
      c(ja+1) = 0.0
      c(jb  ) = a(1) - 0.5 * (a(2) + a(3))
      c(jb+1) =      SIN60 * (a(3) - a(2))

      ja = 3         !  3, 5, 7, ... ,31
      jb = jb + 2    ! 67,69,71, ... ,95
      jc = jc - 2    ! 63,61,59, ... ,35

      do i = 4 , n-8 , 6 ! 88
         c1 = trigs(ja  )
         s1 = trigs(ja+1)
         c2 = trigs(ja+ja-1)
         s2 = trigs(ja+ja  )
         a1 = (c1*a(i+1)+s1*a(i+4))+(c2*a(i+2)+s2*a(i+5))
         b1 = (c1*a(i+4)-s1*a(i+1))+(c2*a(i+5)-s2*a(i+2))
         a2 = a(i  ) - 0.5 * a1
         b2 = a(i+3) - 0.5 * b1
         a3 = SIN60*((c1*a(i+1)+s1*a(i+4))-(c2*a(i+2)+s2*a(i+5)))
         b3 = SIN60*((c1*a(i+4)-s1*A(i+1))-(c2*a(i+5)-s2*a(i+2)))
         c(ja  ) = a(i  ) + a1
         c(ja+1) = a(i+3) + b1
         c(jb  ) = a2 + b3
         c(jb+1) = b2 - a3
         c(jc  ) = a2 - b3
         c(jc+1) =-b2 - a3
         ja = ja + 2
         jb = jb + 2
         jc = jc - 2
      enddo

      if (ja <= jc) then ! ja=33  jc=33
         c(ja  ) = a(n-2) + 0.5 * (a(n-1) - a(n)) ! 33
         c(ja+1) =       -SIN60 * (a(n-1) + a(n)) ! 34
      endif
      a(:) = c(:)
      return
      end subroutine dfft3

!     ================
!     SUBROUTINE DFFT4
!     ================

      subroutine dfft4(a,trigs,n,lot,la)
      real(kind=8), intent(inout) :: a(n,lot)
!f2py intent(in,out) :: a
      real(kind=8) c(n,lot)
      real(kind=8), intent(in) :: trigs(n)
      integer, intent(in) :: n
      integer, intent(in) :: lot
      integer, intent(inout) :: la
!f2py intent(in,out) :: la
      la = la / 4

      i1 = la
      i2 = la + i1
      i3 = la + i2
      i4 = la + i3
      i5 = la + i4
      i6 = la + i5
      i7 = la + i6

      j1 = n/2 - la
      j2 = n - la
      j3 = j1
      j5 = j1 + la

      do i=1,la
         do l=1,lot
         a0p2 = a(i   ,l) + a(i2+i,l)
         a1p3 = a(i1+i,l) + a(i3+i,l)
         c(   i,l) = a0p2 + a1p3
         c(j2+i,l) = a0p2 - a1p3
         c(j1+i,l) = a(   i,l) - a(i2+i,l)
         c(j5+i,l) = a(i3+i,l) - a(i1+i,l)
         enddo
      enddo

      jink = 2 * la
      j0 = la
      j1 = j1 + jink
      j2 = j2 - jink
      j3 = j3 - jink
      j4 = j0 + la
      j5 = j1 + la
      j6 = j2 + la
      j7 = j3 + la

      ibase=4*la

      do 450 k=la,(n-4)/8,la
         kb=k+k
         kc=kb+kb
         kd=kc+kb
         c1=trigs(kb+1)
         s1=trigs(kb+2)
         c2=trigs(kc+1)
         s2=trigs(kc+2)
         c3=trigs(kd+1)
         s3=trigs(kd+2)

         i=ibase+1
         do j=1,la
            do l=1,lot
            a1p5 = c1 * a(i1+i,l) + s1 * a(i5+i,l)
            a2p6 = c2 * a(i2+i,l) + s2 * a(i6+i,l)
            a3p7 = c3 * a(i3+i,l) + s3 * a(i7+i,l)
            a5m1 = c1 * a(i5+i,l) - s1 * a(i1+i,l)
            a6m2 = c2 * a(i6+i,l) - s2 * a(i2+i,l)
            a7m3 = c3 * a(i7+i,l) - s3 * a(i3+i,l)
            a0 = a(i,l) + a2p6
            a2 = a(i,l) - a2p6
            a1 = a1p5 + a3p7
            a3 = a3p7 - a1p5
            b0 = a(i4+i,l) + a6m2
            b2 = a(i4+i,l) - a6m2
            b1 = a5m1 + a7m3
            b3 = a5m1 - a7m3
            c(j0+j,l) = a0+a1
            c(j2+j,l) = a0-a1
            c(j4+j,l) = b0+b1
            c(j6+j,l) = b1-b0
            c(j1+j,l) = a2+b3
            c(j3+j,l) = a2-b3
            c(j5+j,l) = a3+b2
            c(j7+j,l) = a3-b2
            enddo
            i=i+1
         enddo

        ibase=ibase+8*la
        j0 = j0 + jink
        j1 = j1 + jink
        j2 = j2 - jink
        j3 = j3 - jink
        j4 = j0 + la
        j5 = j1 + la
        j6 = j2 + la
        j7 = j3 + la
450   continue
      if (j1 <= j2) then
         sin45=sqrt(0.5)
         i=ibase+1
         do j=1,la
            do l=1,lot
            a1p3 = sin45 * (a(i1+i,l) + a(i3+i,l))
            a1m3 = sin45 * (a(i1+i,l) - a(i3+i,l))
            c(j0+j,l) =  a(   i,l) + a1m3
            c(j1+j,l) =  a(   i,l) - a1m3
            c(j4+j,l) = -a(i2+i,l) - a1p3
            c(j5+j,l) =  a(i2+i,l) - a1p3
            enddo
            i=i+1
         enddo
      endif
      if (la == 1) then
         do l=1,lot
         a(1,l) = c(1,l)
         a(2,l) = 0.0
         a(3:n,l) = c(2:n-1,l)
         enddo
      else
         a = c
      endif
      return
      end subroutine dfft4

!     ================
!     SUBROUTINE DFFT8
!     ================

      subroutine dfft8(a,c,n,lot)
      real(kind=8), intent(in ) :: a(n*lot)
      real(kind=8), intent(out) :: c(n*lot)
      integer, intent(in) :: n
      integer, intent(in) :: lot
      
      la = n / 8
      z  = 1.0 / n
      zsin45 = z * sqrt(0.5)

      do i=0,la*lot-1
         i0 = (i/la) * n + mod(i,la) + 1
         i1 = i0 + la
         i2 = i1 + la
         i3 = i2 + la
         i4 = i3 + la
         i5 = i4 + la
         i6 = i5 + la
         i7 = i6 + la

         a0p4 =  a(i0) + a(i4)
         a1p5 =  a(i1) + a(i5)
         a2p6 =  a(i2) + a(i6)
         a3p7 =  a(i3) + a(i7)
         a5m1 =  a(i5) - a(i1)
         a7m3 =  a(i7) - a(i3)
         a0m4 = (a(i0) - a(i4)) * z
         a6m2 = (a(i6) - a(i2)) * z

         a0p4p2p6 = a0p4 + a2p6
         a1p5p3p7 = a1p5 + a3p7
         a7m3p5m1 = (a7m3 + a5m1) * zsin45
         a7m3m5m1 = (a7m3 - a5m1) * zsin45

         c(i0) = z * (a0p4p2p6 + a1p5p3p7)
         c(i7) = z * (a0p4p2p6 - a1p5p3p7)
         c(i3) = z * (a0p4 - a2p6)
         c(i4) = z * (a3p7 - a1p5)
         c(i1) = a0m4 + a7m3m5m1
         c(i5) = a0m4 - a7m3m5m1
         c(i2) = a7m3p5m1 + a6m2
         c(i6) = a7m3p5m1 - a6m2
      enddo
      return
      end subroutine dfft8

!     ================
!     SUBROUTINE IFFT4
!     ================

      subroutine ifft4(c,trigs,n,lot,la)
      real(kind=8) a(n,lot)
      real(kind=8), intent(inout) :: c(n,lot)
!f2py intent(in,out) :: c
      real(kind=8), intent(in) :: trigs(n)
      integer, intent(in) :: n
      integer, intent(in) :: lot
      integer, intent(inout) :: la
!f2py intent(in,out) :: c

      if (la == 1) then
         a(1,:) = 0.5 * c(1,:)
         a(n,:) = 0.0
         a(2:n-1,:) = c(3:n,:)
      else
         a = c
      endif

      kstop=(n-4)/8

      i1 = n/2 - la
      i2 = n   - la
      i5 = i1  + la

      j1 = la
      j2 = la+j1
      j3 = la+j2
      j4 = la+j3
      j5 = la+j4
      j6 = la+j5
      j7 = la+j6

      do i=1,la
      do l=1,lot
         c(   i,l) = a(i,l) + a(i2+i,l) + a(i1+i,l)
         c(j1+i,l) = a(i,l) - a(i2+i,l) - a(i5+i,l)
         c(j2+i,l) = a(i,l) + a(i2+i,l) - a(i1+i,l)
         c(j3+i,l) = a(i,l) - a(i2+i,l) + a(i5+i,l)
      enddo
      enddo

      iink  = 2 * la
      jbase = 4 * la + 1
      i0    = la
      i1    = i0 + n/2
      i2    = n - 3 * la
      i3    = i2 - n/2
      i4    = i0 + la
      i5    = i1 + la
      i6    = i2 + la
      i7    = i3 + la

      do 450 k=la,kstop,la
        kb=k+k
        kc=kb+kb
        kd=kc+kb
        c1=trigs(kb+1)
        s1=trigs(kb+2)
        c2=trigs(kc+1)
        s2=trigs(kc+2)
        c3=trigs(kd+1)
        s3=trigs(kd+2)
        do i = 1 , la
           j = jbase
           do l=1,lot
           a0p2 = a(i0+i,l) + a(i2+i,l)
           a0m2 = a(i0+i,l) - a(i2+i,l)
           a1p3 = a(i1+i,l) + a(i3+i,l)
           a1m3 = a(i1+i,l) - a(i3+i,l)
           a4p6 = a(i4+i,l) + a(i6+i,l)
           a4m6 = a(i4+i,l) - a(i6+i,l)
           a5p7 = a(i5+i,l) + a(i7+i,l)
           a5m7 = a(i5+i,l) - a(i7+i,l)

           a0p2m1p3 = a0p2 - a1p3
           a4m6m5m7 = a4m6 - a5m7

           c(   j,l) = a0p2 + a1p3
           c(j4+j,l) = a4m6 + a5m7
           c(j2+j,l) = c2 * a0p2m1p3 - s2 * a4m6m5m7
           c(j6+j,l) = s2 * a0p2m1p3 + c2 * a4m6m5m7
           c(j1+j,l) = c1*(a0m2-a5p7)-s1*(a4p6+a1m3)
           c(j5+j,l) = s1*(a0m2-a5p7)+c1*(a4p6+a1m3)
           c(j3+j,l) = c3*(a0m2+a5p7)-s3*(a4p6-a1m3)
           c(j7+j,l) = s3*(a0m2+a5p7)+c3*(a4p6-a1m3)
           enddo
           jbase=jbase+1
        enddo
        i0 = i0 + iink
        i1 = i1 + iink
        i2 = i2 - iink
        i3 = i3 - iink
        i4 = i4 + iink
        i5 = i5 + iink
        i6 = i6 - iink
        i7 = i7 - iink
        jbase=jbase+7*la
450     continue

      if (i1 <= i2) then
         sin45=sqrt(0.5)
         do i=1,la
            j=jbase
            do l=1,lot
            c(   j,l)=a(i0+i,l)+a(i1+i,l)
            c(j1+j,l)=sin45*((a(i0+i,l)-a(i1+i,l))-(a(la+i0+i,l)+&
     &                a(la+i1+i,l)))
            c(j2+j,l)=a(la+i1+i,l)-a(la+i0+i,l)
            c(j3+j,l)=-sin45*((a(i0+i,l)-a(i1+i,l))+(a(la+i0+i,l)+&
     &                a(la+i1+i,l)))
            enddo
            jbase=jbase+1
         enddo
      endif
      la = la * 4
      return
      end subroutine ifft4

!     ================
!     SUBROUTINE IFFT2
!     ================

      subroutine ifft2(a,trigs,n,lot,la)
      real(kind=8), intent(inout):: a(n,lot)
!f2py intent(in,out) :: a
      real(kind=8) c(n,lot)
      real(kind=8), intent(in) :: trigs(n)
      integer, intent(in) :: n
      integer, intent(in) :: lot
      integer, intent(inout) :: la
!f2py intent(in,out) :: la
      integer j,ia,ib
    
      c(1,:) = 0.5 * a(1,:)
      c(2,:) = c(1,:)

      ia    =   3
      ib    = n-1

      do j = 3 , n-5 , 4
         c1 = trigs(ia  )
         s1 = trigs(ia+1)
         do l=1,lot
         amb = a(ia  ,l) - a(ib  ,l)
         apb = a(ia+1,l) + a(ib+1,l)
         c(j  ,l) = a(ia  ,l) + a(ib  ,l)
         c(j+2,l) = a(ia+1,l) - a(ib+1,l)
         c(j+1,l) = c1 * amb - s1 * apb
         c(j+3,l) = s1 * amb + c1 * apb
         enddo
         ia = ia + 2
         ib = ib - 2
      enddo
      c(n-1,:) =  a(ia  ,:)
      c(n  ,:) = -a(ia+1,:)

      a(:,:) = c(:,:)
      la = 2
      return
      end subroutine ifft2

!     ================
!     SUBROUTINE IFFT3
!     ================

      subroutine ifft3(a,trigs,n,lot,la)
      real(kind=8), intent(inout) :: a(n,lot)
!f2py intent(in,out) :: a
      real(kind=8) c(n,lot)
      real(kind=8), intent(in) :: trigs(n)
      integer, intent(in) :: n
      integer, intent(in) :: lot
      integer, intent(inout) :: la
!f2py intent(in,out) :: la
      real(kind=8) SIN60
      parameter(SIN60 = 0.866025403784438D0)

      ib = 2 * (n/3) + 1

      c(1,:) = 0.5 * a(1,:) + a(ib,:)
      c(2,:) = 0.5 * a(1,:) - 0.5 * a(ib,:) - SIN60 * a(ib+1,:)
      c(3,:) = 0.5 * a(1,:) - 0.5 * a(ib,:) + SIN60 * a(ib+1,:)

      ia = 3
      ic = ib - 2
      ib = ib + 2

      do j = 4 , n-8 , 6
         c1 = trigs(ia  )
         s1 = trigs(ia+1)
         c2 = trigs(ia+ia-1)
         s2 = trigs(ia+ia  )

         do l = 1 , lot
            hbpc = a(ia  ,l) - 0.5 * (a(ib  ,l) + a(ic  ,l))
            hbmc = a(ia+1,l) - 0.5 * (a(ib+1,l) - a(ic+1,l))
            sbmc = SIN60 * (a(ib  ,l) - a(ic  ,l))
            sbpc = SIN60 * (a(ib+1,l) + a(ic+1,l))

            c(j  ,l) = a(ia  ,l) + a(ib  ,l) + a(ic  ,l)
            c(j+3,l) = a(ia+1,l) + a(ib+1,l) - a(ic+1,l)
            c(j+1,l) = c1 * (hbpc-sbpc) - s1 * (hbmc+sbmc)
            c(j+4,l) = s1 * (hbpc-sbpc) + c1 * (hbmc+sbmc)
            c(j+2,l) = c2 * (hbpc+sbpc) - s2 * (hbmc-sbmc)
            c(j+5,l) = s2 * (hbpc+sbpc) + c2 * (hbmc-sbmc)
         enddo
         ia = ia + 2
         ib = ib + 2
         ic = ic - 2
      enddo

      c(n-2,:) = a(ia,:)
      c(n-1,:) =   0.5 * a(ia,:) - SIN60 * a(ia+1,:)
      c(n  ,:) = - 0.5 * a(ia,:) - SIN60 * a(ia+1,:)

      a(:,:)  = c(:,:)
      la = 3
      return
      end subroutine ifft3

!     ================
!     SUBROUTINE IFFT8
!     ================

      subroutine ifft8(a,c,n,lot)
      real(kind=8) SQRT2
      parameter(SQRT2 = 1.414213562373095D0)
      integer, intent(in ) :: n
      integer, intent(in ) :: lot
      real(kind=8), intent(in ) :: a(n*lot)
      real(kind=8), intent(out) :: c(n*lot)
      la = n / 8

      do i=0,la*lot-1
         i0 = (i/la) * n + mod(i,la) + 1
         i1 = i0 + la
         i2 = i1 + la
         i3 = i2 + la
         i4 = i3 + la
         i5 = i4 + la
         i6 = i5 + la
         i7 = i6 + la

         a0p7 = a(i0) + a(i7)
         a0m7 = a(i0) - a(i7)
         a1p5 = a(i1) + a(i5)
         a1m5 = a(i1) - a(i5)
         a2p6 = a(i2) + a(i6)
         a2m6 = a(i2) - a(i6)

         a0p7p3   = a0p7 + a(i3)
         a0p7m3   = a0p7 - a(i3)
         a0m7p4   = 2.0 * (a0m7 + a(i4))
         a0m7m4   = 2.0 * (a0m7 - a(i4))
         a1m5p2p6 = SQRT2 * (a1m5 + a2p6)
         a1m5m2p6 = SQRT2 * (a1m5 - a2p6)

         c(i0)  = 2.0 * (a0p7p3 + a1p5)
         c(i2)  = 2.0 * (a0p7m3 - a2m6)
         c(i4)  = 2.0 * (a0p7p3 - a1p5)
         c(i6)  = 2.0 * (a0p7m3 + a2m6)

         c(i1)  = a0m7m4 + a1m5m2p6
         c(i3)  = a0m7p4 - a1m5p2p6
         c(i5)  = a0m7m4 - a1m5m2p6
         c(i7)  = a0m7p4 + a1m5p2p6
      enddo
      return
      end
