module module_chem_utilities
   use module_domain
   use module_model_constants
   use module_state_description
   use module_configure

! Add the calculation of pi_phy in chem_prep
contains
   subroutine chem_prep(config_flags, &
                        u, v, p, pb, alt, ph, &
                        phb, t, moist, n_moist, &
                        rho, p_phy, pi_phy, &
                        u_phy, v_phy, p8w, t_phy, t8w, &
                        z, z_at_w, dz8w, rh, &
                        fzm, fzp, &
                        ids, ide, jds, jde, kds, kde, &
                        ims, ime, jms, jme, kms, kme, &
                        its, ite, jts, jte, kts, kte)
      implicit none

      type(grid_config_rec_type), intent(in) :: config_flags
      integer, intent(in) :: ids, ide, jds, jde, kds, kde, &
                             ims, ime, jms, jme, kms, kme, &
                             its, ite, jts, jte, kts, kte
      integer, intent(in) :: n_moist

      real, dimension(ims:ime, kms:kme, jms:jme, n_moist), intent(in) :: moist

      real, dimension(ims:ime, kms:kme, jms:jme), &
         intent(out) :: u_phy, v_phy, p_phy, p8w, &
                        t_phy, t8w, rho, z, dz8w, &
                        rh, z_at_w, pi_phy

      real, dimension(ims:ime, kms:kme, jms:jme), &
         intent(in) :: pb, &
                       p, &
                       u, &
                       v, &
                       alt, &
                       ph, &
                       phb, &
                       t

      real, dimension(kms:kme), intent(in) :: fzm, &
                                              fzp

      integer :: i_start, i_end, j_start, j_end, k_start, k_end
      integer :: i, j, k
      real    :: w1, w2, z0, z1, z2

!-----------------------------------------------------------------------
!  set up loop bounds for this grid's boundary conditions

      i_start = its
      i_end = min(ite, ide - 1)
      j_start = jts
      j_end = min(jte, jde - 1)

      k_start = kts
      k_end = min(kte, kde - 1)

!  compute thermodynamics and velocities at pressure points
      do j = j_start, j_end
      do k = k_start, k_end
      do i = i_start, i_end

         p_phy(i, k, j) = p(i, k, j) + pb(i, k, j)
         t_phy(i, k, j) = (t(i, k, j) + t0)*(p_phy(i, k, j)/p1000mb)**rcp
         rho(i, k, j) = 1./alt(i, k, j)*(1.+moist(i, k, j, P_QV))
         u_phy(i, k, j) = 0.5*(u(i, k, j) + u(i + 1, k, j))
         v_phy(i, k, j) = 0.5*(v(i, k, j) + v(i, k, j + 1))
         pi_phy(i, k, j) = (p_phy(i, k, j)/p1000mb)**rcp

      enddo
      enddo
      enddo

!  wig: added to make sure there is no junk in the top level even
!       though it should not be used
      do j = j_start, j_end
      do i = i_start, i_end
         p_phy(i, kte, j) = p_phy(i, k_end, j)
         t_phy(i, kte, j) = t_phy(i, k_end, j)
         rho(i, kte, j) = rho(i, k_end, j)
         u_phy(i, kte, j) = u_phy(i, k_end, j)
         v_phy(i, kte, j) = v_phy(i, k_end, j)
         pi_phy(i, kte, j) = pi_phy(i, k_end, j)         

      enddo
      enddo

!  compute z at w points
      do j = j_start, j_end
      do k = k_start, kte
      do i = i_start, i_end
         z_at_w(i, k, j) = (phb(i, k, j) + ph(i, k, j))/g
      enddo
      enddo
      enddo

      do j = j_start, j_end
      do k = k_start, kte - 1
      do i = i_start, i_end
         dz8w(i, k, j) = z_at_w(i, k + 1, j) - z_at_w(i, k, j)
      enddo
      enddo
      enddo

      do j = j_start, j_end
      do i = i_start, i_end
         dz8w(i, kte, j) = 0.
      enddo
      enddo

!  compute z at p points (average of z at w points)
      do j = j_start, j_end
      do k = k_start, k_end
      do i = i_start, i_end
         z(i, k, j) = 0.5*(z_at_w(i, k, j) + z_at_w(i, k + 1, j))
         rh(i, k, j) = max(.1, MIN(.95, moist(i, k, j, p_qv)/ &
                                   (3.80*exp(17.27*(t_phy(i, k, j) - 273.)/ &
                                             (t_phy(i, k, j) - 36.))/(.01*p_phy(i, k, j)))))
      enddo
      enddo
      enddo

!  interp t and p at w points

      do j = j_start, j_end
      do k = 2, k_end
      do i = i_start, i_end
         p8w(i, k, j) = fzm(k)*p_phy(i, k, j) + fzp(k)*p_phy(i, k - 1, j)
         t8w(i, k, j) = fzm(k)*t_phy(i, k, j) + fzp(k)*t_phy(i, k - 1, j)
      enddo
      enddo
      enddo

!  extrapolate p and t to surface and top.
!  we'll use an extrapolation in z for now

      do j = j_start, j_end
      do i = i_start, i_end

! bottom
         z0 = z_at_w(i, 1, j)
         z1 = z(i, 1, j)
         z2 = z(i, 2, j)
         w1 = (z0 - z2)/(z1 - z2)
         w2 = 1.-w1
         p8w(i, 1, j) = w1*p_phy(i, 1, j) + w2*p_phy(i, 2, j)
         t8w(i, 1, j) = w1*t_phy(i, 1, j) + w2*t_phy(i, 2, j)

! top
         z0 = z_at_w(i, kte, j)
         z1 = z(i, k_end, j)
         z2 = z(i, k_end - 1, j)
         w1 = (z0 - z2)/(z1 - z2)
         w2 = 1.-w1

!    p8w(i,kde,j) = w1*p_phy(i,kde-1,j)+w2*p_phy(i,kde-2,j)
!!!  bug fix      extrapolate ln(p) so p is positive definite
         p8w(i, kde, j) = exp(w1*log(p_phy(i, kde - 1, j)) + w2*log(p_phy(i, kde - 2, j)))
         t8w(i, kde, j) = w1*t_phy(i, kde - 1, j) + w2*t_phy(i, kde - 2, j)

      enddo
      enddo
   END SUBROUTINE chem_prep

   SUBROUTINE calc_zenithb(lat, long, ijd, dtstep, gmt, curr_secs1, zenith)

      ! input 
      ! lat - latitude in decimal degrees
      ! long - longitude in decimal degrees
      ! NOTE: Nonstandard convention for long: >0 for W, <0 for E!!
      ! ijd - grid%julday current julian day
      ! dtstep - grid%dt
      ! gmt - grid%gmt
      ! curr_secs1 - curr_secs number of seconds into the simulation
      ! output
      ! zenith - in radians

      IMPLICIT NONE

      REAL, INTENT (IN) :: lat, long, dtstep, gmt

      real(KIND=8) :: curr_secs1

      INTEGER, INTENT (IN) :: ijd

      REAL, INTENT (OUT) :: zenith

      ! .. Scalar Arguments ..
      REAL :: csz, cw, d, decl, dr, ec, epsi, eqt, eyt, feqt, feqt1, &
         feqt2, feqt3, feqt4, feqt5, feqt6, feqt7, lbgmt, lzgmt, ml, pepsi, &
         pi, ra, rdecl, reqt, rlt, rml, rphi, rra, ssw, sw, tab, w, wr, &
         yr, zpt, zr, tmidh, gmtp, xmin,yt

      REAL(KIND=8) :: xtime, xhour

      INTEGER :: ixhour, jd

        xtime = curr_secs1/60._8 + real(dtstep/120.,8)
        ixhour = int(gmt + 0.01) + int(xtime/60._8)
        xhour=real(ixhour,8)
        xmin = 60.*gmt +real(xtime-xhour*60._8,8)
        gmtp=MOD(xhour,24._8)
        tmidh = gmtp + xmin/60.

      ! .. Intrinsic Functions ..
      ! INTRINSIC acos, atan, cos, min, sin, tan
      ! convert to radians
        pi = 3.1415926535590
        dr = pi/180.
        rlt = lat*dr
        rphi = long*dr

        ! ???? + (yr - yref)

        jd = ijd

        d = jd + tmidh/24.0
        ! calc geom mean longitude
        ml = 279.2801988 + .9856473354*d + 2.267E-13*d*d
        rml = ml*dr

        ! calc equation of time in sec
        ! w = mean long of perigee
        ! e = eccentricity
        ! epsi = mean obliquity of ecliptic
        w = 282.4932328 + 4.70684E-5*d + 3.39E-13*d*d
        wr = w*dr
        ec = 1.6720041E-2 - 1.1444E-9*d - 9.4E-17*d*d
        epsi = 23.44266511 - 3.5626E-7*d - 1.23E-15*d*d
        pepsi = epsi*dr
        yt = (tan(pepsi/2.0))**2
        cw = cos(wr)
        sw = sin(wr)
        ssw = sin(2.0*wr)
        eyt = 2.*ec*yt
        feqt1 = sin(rml)*(-eyt*cw-2.*ec*cw)
        feqt2 = cos(rml)*(2.*ec*sw-eyt*sw)
        feqt3 = sin(2.*rml)*(yt-(5.*ec**2/4.)*(cw**2-sw**2))
        feqt4 = cos(2.*rml)*(5.*ec**2*ssw/4.)
        feqt5 = sin(3.*rml)*(eyt*cw)
        feqt6 = cos(3.*rml)*(-eyt*sw)
        feqt7 = -sin(4.*rml)*(.5*yt**2)
        feqt = feqt1 + feqt2 + feqt3 + feqt4 + feqt5 + feqt6 + feqt7
        eqt = feqt*13751.0

        ! convert eq of time from sec to deg
        reqt = eqt/240.
        ! calc right ascension in rads
        ra = ml - reqt
        rra = ra*dr
        ! calc declination in rads, deg
        tab = 0.43360*sin(rra)
        rdecl = atan(tab)
        decl = rdecl/dr
        ! calc local hour angle
        lbgmt = 12.0 - eqt/3600. + long*24./360.
        lzgmt = 15.0*(tmidh-lbgmt)
        zpt = lzgmt*dr
        csz = sin(rlt)*sin(rdecl) + cos(rlt)*cos(rdecl)*cos(zpt)
        csz = min(1.,csz)
        zr = acos(csz)
      ! keep zenith angle in radians for later use (GJF 6/2004)
        zenith = zr

        RETURN

   END SUBROUTINE calc_zenithb

   SUBROUTINE getpar( tsolar, pres, zen, pardb, pardif )

!***********************************************************************
!  subroutine body starts at line  
!
!  DESCRIPTION:
!  
!        Based on code from Bart Brashers (10/2000), which was based on
!        code from Weiss and Norman (1985).  
!
!
!  PRECONDITIONS REQUIRED:
!     Solar radiation (W/m2) and pressure (mb)
!
!  SUBROUTINES AND FUNCTIONS CALLED:
!
!  REVISION  HISTORY:
!    3/01 : Prototype by JMV
! 
!***********************************************************************
!
! Project Title: Sparse Matrix Operator Kernel Emissions (SMOKE) Modeling
!                System
! File: @(#)Id: getpar.f,v 1.1.1.1 2001/03/27 19:08:49 smith_w Exp 
!
! COPYRIGHT (C) 2001, MCNC--North Carolina Supercomputing Center
! All Rights Reserved
!
! See file COPYRIGHT for conditions of use.
!
! MCNC-Environmental Programs Group
! P.O. Box 12889
! Research Triangle Park, NC  27709-2889
!
! env_progs@mcnc.org
!
! Pathname: Source: /env/proj/archive/cvs/jmv/beis3v0.9/getpar.f,v 
! Last updated: Date: 2001/03/27 19:08:49  
!
!***********************************************************************
      IMPLICIT NONE

!........ Inputs

        REAL , INTENT  (IN) :: tsolar   ! modeled or observed total radiation (W/m2)
        REAL , INTENT  (IN) :: pres     ! atmospheric pressure (mb)
        REAL , INTENT  (IN) :: zen      ! solar zenith angle (radians)

!........ Outputs

        REAL, INTENT (OUT) :: pardb     ! direct beam PAR( umol/m2-s) now (W/m^2)
        REAL, INTENT (OUT) :: pardif    ! diffuse PAR ( umol/m2-s) now (W/m^2)

!...........   PARAMETERS and their descriptions:

!        REAL, PARAMETER :: watt2umol = 4.6  ! convert W/m^2 to umol/m^2-s (4.6)

!      
        REAL ratio              ! transmission fraction for total radiation
        REAL ot                 ! optical thickness
        REAL rdvis              ! possible direct visible beam (W/m^2)
        REAL rfvis              ! possible visible diffuse (W/m^2)
        REAL wa                 ! water absorption in near-IR (W/m^2)
        REAL rdir               ! direct beam in near-IR (W/m^2)
        REAL rfir               ! diffuse near-IR (W/m^2)
        REAL rvt                ! total possible visible radiation (W/m^2)
        REAL rirt               ! total possible near-IR radiation (W/m^2)
        REAL fvis               ! fraction of visible to total 
        REAL fvb                ! fraction of visible that is direct beam
        REAL fvd                ! fraction of visible that is diffuse

!***************************************
!   begin body of subroutine  

!............ Assume that PAR = 0 if zenith angle is greater than 87 degrees
!............ or if solar radiation is zero

        IF (zen .GE. 1.51844 .OR. tsolar .LE. 0.) THEN
             pardb  = 0.
             pardif = 0.
             RETURN
        ENDIF

!............ Compute clear sky (aka potential) radiation terms

        ot    = pres / 1013.25 / COS(zen)              !Atmospheric Optical thickness
        rdvis = 600. * EXP(-.185 * ot) * COS(zen)      !Direct visible beam, eqn (1)
        rfvis = 0.42 * (600 - rdvis) * COS(zen)        !Visible Diffuse, eqn (3)
        wa    = 1320 * .077 * (2. * ot)**0.3           !water absorption in near-IR, eqn (6)
        rdir  = (720. * EXP(-0.06 * ot)-wa) * COS(zen) !Direct beam near-IR, eqn (4)
        rfir  = 0.65 * (720. - wa - rdir) * COS(zen)   !Diffuse near-IR, eqn (5)
        rvt   = rdvis + rfvis                    !Total visible radiation, eqn (9)
        rirt  = rdir + rfir                      !Total near-IR radiation, eqn (10) 
        fvis  = rvt/(rirt + rvt)                 !Fraction of visible to total radiation, eqn 7
        ratio = tsolar /(rirt + rvt)             !Ratio of "actual" to clear sky solar radiation

!............ Compute fraction of visible that is direct beam

        IF (ratio .GE. 0.89) THEN
           fvb = rdvis/rvt * 0.941124
        ELSE IF (ratio .LE. 0.21) THEN
           fvb = rdvis/rvt * 9.55E-3
        ELSE
           fvb = rdvis/rvt * (1.-((0.9 - ratio)/0.7)**0.666667)
        ENDIF
        fvd = 1. - fvb

!............ Compute PAR (direct beam and diffuse) in umol/m2-sec

        pardb  = tsolar * fvis * fvb
        pardif = tsolar * fvis * fvd

        RETURN

   END SUBROUTINE getpar
      
   SUBROUTINE get_cloud_optical_depth(t3d, p8w3d, qc3d, qi3d, qndrop3d,                                   &
                                      taucldc, taucldi, optd,                                             &
                                      f_qc, f_qi, f_qndrop, warm_rain,                                    &
                                      ims, ime, jms, jme, kms, kme,                                       &
                                      its, ite, jts, jte, kts, kte                                        )

   ! calculating the ice/liquid cloud optical depth (visible broadband) baed on
   ! the GSFC radiation parameterization
     
         integer, parameter                                                :: fp_kind = selected_real_kind(15)   !Real precision (6-single 15-double 20-quad)
         real, dimension(ims:ime, kms:kme, jms:jme), intent(in)            :: p8w3d,                         &   ! pressure at full levels (Pa)
                                                                              t3d                                ! temperature (K)              
         real, dimension(ims:ime, kms:kme, jms:jme), optional, intent(in)  :: qc3d,                          &   ! cloud water mixint ratio (kg/kg)
                                                                              qi3d,                          &   ! cloud ice mixint ratio (kg/kg)
                                                                              qndrop3d                             !                                                                                                                              
         real, dimension(ims:ime, kms:kme, jms:jme), intent(inout)         :: taucldc,                       &   ! liquid cloud optical depth
                                                                              taucldi,                       &   ! ice cloud optical depth
                                                                              optd                               ! cloud optical depth
         logical, optional,                          intent(in)            :: f_qc, f_qi, f_qndrop, warm_rain
         integer,                                    intent(in)            :: ims, ime, jms, jme, kms, kme,  &   
                                                                              its, ite, jts, jte, kts, kte    
         real(kind=fp_kind), parameter                                     :: re = 10.                           ! cloud droplet effective radius (micron)
         
   !-----------------------------local variable-----------------------------
         integer :: i, j, k, nk
         real(kind = fp_kind), dimension(its:ite, kts-1:kte, 3)  :: reff,  &  ! 1-ice cloud, 2-liquid cloud, 3-rain particale effective size (micron)
                                                                    cwc       ! hydrometer mixing ratio (kg/kg) or (g/g)
         real(kind = fp_kind), dimension(its:ite, kts-1:kte+1)   :: p8w2d     ! pressure in full level (mb)
         real(kind = fp_kind), dimension(its:ite, kts-1:kte)     :: t2d       ! temperature (K)
         real(kind = fp_kind), dimension(its:ite, kts-1:kte)     :: qndrop2d
         real(kind = fp_kind)                                    :: x, lwpmin, pi, third, relconst, rhoh2o
         logical                                                 :: predicate                                                                                     
         
     j_loop: do j = jts,jte
   !------------------------------------------------------------------------
        do k = kts, kte+1
           do i = its, ite
              taucldc(i, k, j) = 0._fp_kind
              taucldi(i, k, j) = 0._fp_kind
              optd(i, k, j) = 0._fp_kind
           enddo
        enddo
            
        do k = kts, kte
           do i = its, ite
              cwc(i, k, 1) = 0.
              cwc(i, k, 2) = 0.
           enddo
        enddo
                
   ! reverse vars
        do k = kts, kte+1
           do i = its, ite
              nk = kme-k+kms
              p8w2d(i, k) = p8w3d(i, nk, j)*0.01 ! p8w2d in hPa
           enddo
        enddo
                
        do i = its, ite
           p8w2d(i, 0) = .0
        enddo

        do k = kts, kte
           do i = its, ite
              nk = kme-1-k+kms
              t2d(i, k) = t3d(i, nk, j)
              cwc(i, k, 2) = qc3d(i, nk, j)
              cwc(i, k, 2) = max(0., cwc(i, k, 2))
           enddo
        enddo

        if ( present( f_qi ) ) then
           predicate = f_qi
        else
           predicate = .false.
        endif
   
        if (.not. warm_rain .and. .not. predicate ) then
           do k = kts, kte
              do i = its, ite
                 if (t2d(i, k) .lt. 273.15) then
                    cwc(i, k, 1) = cwc(i, k, 2)
                    cwc(i, k, 2) = 0.
                 endif
              enddo
           enddo
        endif

        if ( present( f_qndrop ) ) then
           if ( f_qndrop ) then
              do k = kts, kte
                 do i = its, ite
                    nk = kme-1-k+kms
                    qndrop2d(i,k) = qndrop3d(i, nk, j)
                 enddo
              enddo
              qndrop2d(:,kts-1) = 0.
           endif
        endif

        do i = its, ite
           t2d(i, 0) = t2d(i, 1)
           cwc(i, 0, 2) = 0.
           cwc(i, 0, 1) = 0.
        enddo

        if ( present( f_qi ) .and. present( qi3d ) ) then
           if ( f_qi ) then
              do k = kts, kte
                 do i = its, ite
                    nk = kme-1-k+kms
                    cwc(i, k, 1) = qi3d(i, nk, j)
                    cwc(i, k, 1) = max(0., cwc(i, k, 1))
                 enddo
              enddo
           endif
        endif
                  
   ! cloud drop effective radius based on GSFCSW

        lwpmin = 3.e-5
        pi = 4.*atan(1.0)
        third=1./3.
        rhoh2o=1.e3
        relconst=3/(4.*pi*rhoh2o)

        do k = kts-1, kte
           do i = its, ite
              reff(i, k, 2) = re
              if ( present( f_qndrop ) )then
                if( f_qndrop ) then
                  if( cwc(i, k, 2) * (p8w2d(i, k+1)-p8w2d(i, k)) .gt. lwpmin .and. &
                      qndrop2d(i, k) .gt. 1000. ) then
                      reff(i, k, 2) = (relconst*cwc(i, k, 2)/qndrop2d(i, k))**third ! effective radius in m 
                      ! apply scaling from Martin et al., JAS 51, 1830.
                      reff(i, k, 2) = 1.1*reff(i, k, 2)
                      reff(i, k, 2) = reff(i, k, 2)*1.e6 ! convert from m to microns
                      reff(i, k, 2) = max(reff(i, k, 2), 4.)
                      reff(i, k, 2) = min(reff(i, k, 2), 20.)
                  endif
                endif
              endif
              reff(i, k, 1) = 80.
           enddo
        enddo

   ! cloud water mixint ratio   

        if ( present( f_qc ) .and. present( qc3d ) )then
           do k = kts,kte
              do i = its, ite
                 nk = kme-1-k+kms
                    x=1.02*10000.*( p8w2d(i, k+1) - p8w2d(i, k) )
                    taucldc(i, nk, j) = x * cwc(i, k, 2) * ( -6.59e-3 + 1.65/reff(i, k, 2) )
              enddo
           enddo
        endif

   ! cloud ice mixing ratio

        if ( present( f_qi ) .and. present( qi3d ) ) then
           do k = kts, kte
              do i = its, ite
                 nk = kme-1-k+kms
                      reff(i, k, 1) = 125. + (t2d(i, k) - 243.16)*5.   ! ice effective radius depends on temperature
                      reff(i, k, 1) = min( 125._fp_kind, max(25._fp_kind, reff(i, k, 1)) )
                      x = 1.02*10000. * (p8w2d(i, k+1) - p8w2d(i, k) )
                      taucldi(i, nk, j) = x * cwc(i, k, 1) * ( 3.33e-4 + 2.52/reff(i, k, 1) )
               enddo
            enddo
        endif
                
        do k = kts, kte
           do i = its, ite
              optd(i, k, j) = taucldc(i, k, j) + taucldi(i, k, j)
           enddo
        enddo
                                                    
     enddo j_loop
   
   END SUBROUTINE get_cloud_optical_depth
   
   SUBROUTINE  calc_slp(t, p, pb, qv, ph, phb,                                   &
                        slp,                                                     &
                        ids, ide, jds, jde, kds, kde,                            &
                        ims, ime, jms, jme, kms, kme,                            &
                        its, ite, jts, jte, kts, kte                             )

! calculating the sea level pressure

            integer, parameter                                                :: fp_kind = selected_real_kind(15)     !Real precision (6-single 15-double 20-quad)
            real, dimension(ims:ime, kms:kme, jms:jme), intent(in)            :: p,                               &   ! base state pressure (Pa)
                                                                                 pb,                              &   ! perturbation pressure (Pa) 
                                                                                 t,                               &   ! perturbation potential temperature (K)
                                                                                 ph,                              &   ! perturbation geopotential (m^2 s^-2)
                                                                                 phb                                  ! base state geopotential (m^2 s^-2)                                                                                                                                                               
            real, dimension(ims:ime, kms:kme, jms:jme), optional, intent(in)  :: qv                                   !  Water vapor mixing ratio(kg/kg)

            real, dimension(ims:ime, jms:jme), intent(out)                    :: slp                                  ! sea level pressure (hPa)

            integer,                           intent(in)                     :: ids, ide, jds, jde, kds, kde,      &
                                                                                 ims, ime, jms, jme, kms, kme,      &
                                                                                 its, ite, jts, jte, kts, kte
! local variables
            real, dimension(its:ite, kts:kte, jts:jte)                        :: p_phy, t_phy, hgt
            real, dimension(its:ite, jts:jte)                                 :: level, plo, phi, tlo, thi, zlo, zhi
            integer                                                           :: i, j, k,                           &
                                                                               i_start, i_end, j_start, j_end,    &
                                                                               k_start, k_end, klo, khi

            real                                                              :: TC, PCONST, Rd, GAMMA1,            &
                                                                                 p_at_pconst, t_at_pconst,          &
                                                                                 z_at_pconst, t_surf, t_sea_lev, z_half_lowest

            logical                                                         :: traditional_comp, l1, l2, l3

!  set up loop bounds for this grid's boundary conditions

            i_start = its
            i_end = min(ite, ide - 1)
            j_start = jts
            j_end = min(jte, jde - 1)
            k_start = kts
            k_end = min(kte, kde - 1)

            TC = 273.16+17.5
            PCONST = 10000.
            Rd = 287.04
            GAMMA1 = .0065
            traditional_comp = .true.
            
!  compute pressure (Pa), temperature (K) and height (m) at pressure points
           do j = j_start, j_end
           do k = k_start, k_end
           do i = i_start, i_end
              p_phy(i, k, j) = p(i, k, j) + pb(i, k, j)
              t_phy(i, k, j) = (t(i, k, j) + t0)*(p_phy(i, k, j)/p1000mb)**rcp
           enddo
           enddo
           enddo

!  wig: added to make sure there is no junk in the top level even
!       though it should not be used      
           do j = j_start, j_end
           do i = i_start, i_end
              p_phy(i, kte, j) = p_phy(i, k_end, j)
              t_phy(i, kte, j) = t_phy(i, k_end, j)
           enddo
           enddo
           
! compute height (m) at pressure points   
           do j = j_start, j_end
           do k = k_start, kte
           do i = i_start, i_end
              hgt(i, k, j) = (phb(i, k, j) + ph(i, k, j))/g
           enddo
           enddo
           enddo

!  Calculation of SLP        
           do j = jts, jte
              do i = its, ite
                 level(i, j) = -1
                 do k = kts, kte
                    if( p_phy(i, k, j) .lt. p_phy(i, 1, j) - PCONST) then
                        level(i, j) = k
                    exit
                    endif
                 enddo
                    if( level(i, j) .eq. -1) then
                        call wrf_message("Error in finding 100hPa up")
                    endif
              enddo
           enddo

           do j = jts, jte
           do i = its, ite
              klo = max( level(i, j)-1, 1.0)
              khi = min( klo+1, kte-1)
              if( klo .eq. khi) then
                  call wrf_message("Trapping levels are weird and they should not be equal")
              endif
              plo(i, j) = p_phy(i, klo, j)
              phi(i, j) = p_phy(i, khi, j)
              tlo(i, j) = t_phy(i, klo, j)*(1.+.608*qv(i, klo, j))
              thi(i, j) = t_phy(i, khi, j)*(1.+.608*qv(i, klo, j))
              zlo(i, j) = hgt(i, klo, j)
              zhi(i, j) = hgt(i, khi, j)
           enddo
           enddo

           do j = jts, jte
           do i = its, ite
              p_at_pconst = p_phy(i, 1, j) - PCONST;
              t_at_pconst = thi(i, j) - (thi(i, j) - tlo(i, j)) * log(p_at_pconst/phi(i, j)) * log(plo(i, j)/phi(i, j))
              z_at_pconst = zhi(i, j) - (zhi(i, j) - zlo(i, j)) * log(p_at_pconst/phi(i, j)) * log(plo(i, j)/phi(i, j))
              t_surf = t_at_pconst*(p_phy(i, 1, j)/p_at_pconst)**(GAMMA1*Rd/g)
              t_sea_lev = t_at_pconst + GAMMA1*z_at_pconst

              if(traditional_comp) then
                 if( t_sea_lev .lt. TC) then
                     l1 = .true.
                 endif
                 if( t_surf .lt. TC) then
                     l2 = .true.
                 endif

                  l3 = .not.l1
                  if( l2 .and. l3) then
                      t_sea_lev = TC
                  else
                      t_sea_lev = TC - 0.005*(t_surf - TC)**2
                  endif

              endif
              ! The grand finale
                z_half_lowest = hgt(i, 1, j);
                slp(i, j) = .01 * p_phy(i, 1, j)*exp((2*g*z_half_lowest)/(Rd*(t_sea_lev + t_surf)))
      enddo
      enddo

END SUBROUTINE calc_slp

   subroutine UPCASE( lstring )
!----------------------------------------------------------------------
!       ... Convert character string lstring to upper case
!----------------------------------------------------------------------
   implicit none

!-----------------------------------------------------------------------
!  ... Dummy args
!-----------------------------------------------------------------------
   character(len=*), intent(inout) ::  lstring

!-----------------------------------------------------------------------
!  ... Local variables
!-----------------------------------------------------------------------
   integer :: i

   do i = 1,LEN_TRIM( lstring )
     if( ICHAR(lstring(i:i)) >= 97 .and.  ICHAR(lstring(i:i)) <= 122 ) then
       lstring(i:i) = CHAR(ICHAR(lstring(i:i)) - 32)
     end if
   end do

   end subroutine UPCASE


END MODULE module_chem_utilities
