!------------------------------------------------------------------------------
!                   Harmonized Emissions Component (HEMCO)                    !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: hcox_lightnox_mod.F90
!
! !DESCRIPTION: Module HCOX\_LightNOx\_Mod contains routines to
!  compute NO lightning emissions, according to the GEOS-Chem lightning
!  algorithms.
!\\
!\\
! This is a HEMCO extension module that uses many of the HEMCO core
! utilities. In particular, the LIS-OTD local redistribution factors are
! now read through the HEMCO framework, and the corresponding netCDF
! input file is specified in the HEMCO configuration file. The table of
! cumulative distribution functions used to vertically distribute lightning
! NOx emissions is specified in the extension switch section of the
! configuration file.
!\\
!\\
! References:
! \begin{itemize}
! \item Murray, L. T., Jacob, D. J., Logan, J. A., Hudman, R. C., and
!       Koshak, W. J.: \emph{Optimized regional and interannual variability
!       of lightning in a global chemical transport model con- strained
!       by LIS/OTD satellite data}, \underline{J. Geophys. Res.},
!       Atmospheres, 117, 2012.
! \item Ott, L. E., K. E. Pickering, G. L. Stenchikov, D. J. Allen,
!       A. J. DeCaria, B. Ridley, R.-F. Lin, S. Lang, and W.-K. Tao,
!       \emph{Production of lightning NOx and its vertical distribution
!       calculated  from three-dimensional cloud-scale chemical transport
!       model simulations}, \underline{J. Geophys. Res.}, 115, D04301, 2010.
! \end{itemize}
!
! !INTERFACE:
!
MODULE HCOX_LightNOx_Mod
!
! !USES:
!
  USE HCO_Error_Mod
  USE HCO_Diagn_Mod
  USE HCOX_TOOLS_MOD
  USE HCO_State_Mod,  ONLY : HCO_State
  USE HCOX_State_MOD, ONLY : Ext_State

  IMPLICIT NONE
  PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!
  PUBLIC  :: HCOX_LightNOX_Run
  PUBLIC  :: HCOX_LightNOX_Final
  PUBLIC  :: HCOX_LightNOX_Init
!
! !PRIVATE MEMBER FUNCTIONS:
!
  PRIVATE :: LIGHTNOX
  PRIVATE :: LIGHTDIST
!
! !PUBLIC DATA MEMBERS:
!
!
! !REMARKS:
!  %%% NOTE: MFLUX and PRECON methods are now deprecated (ltm, bmy, 7/9/09)
!                                                                             .
!  References:
!  ============================================================================
!  (1 ) Price & Rind (1992), JGR, vol. 97, 9919-9933.
!  (2 ) Price & Rind (1994), M. Weather Rev, vol. 122, 1930-1939.
!  (3 ) Allen & Pickering (2002), JGR, 107, D23, 4711, doi:10.1029/2002JD002066
!  (4 ) Hudman et al (2007), JGR, 112, D12S05, doi:10.1029/2006JD007912
!  (5 ) Sauvage et al, 2007, ACP,
!        http://www.atmos-chem-phys.net/7/815/2007/acp-7-815-2007.pdf
!  (6 ) Ott et al., (2010), JGR
!  (7 ) Allen et al., (2010), JGR
!  (8 ) Murray et al., (2011), in prep.
!
! !REVISION HISTORY:
!  14 Apr 2004 - L. Murray, R. Hudman - Initial version
!  See https://github.com/geoschem/hemco for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !DEFINED PARAMETERS:
!
  INTEGER, PARAMETER            :: NLTYPE        = 4
  INTEGER, PARAMETER            :: NNLIGHT       = 3200
  REAL*8,  PARAMETER            :: RFLASH_MIDLAT = 3.011d26   ! 500 mol/flash applied in the N extratropics
  REAL*8,  PARAMETER            :: RFLASH_TROPIC = 1.566d26   ! 260 mol/flash applied in the tropics / S extratropics
  REAL*8,  PARAMETER            :: EAST_WEST_DIV = -30d0
  REAL*8,  PARAMETER            :: WEST_NS_DIV   =  35d0
  REAL*8,  PARAMETER            :: EAST_NS_DIV   =  35d0
!
! !PRIVATE TYPES:
!
  ! Scalars
  TYPE :: MyInst
   INTEGER                       :: Instance
   INTEGER                       :: IDTNO     ! NO tracer ID
   INTEGER                       :: ExtNr     ! HEMCO Extension ID
   LOGICAL                       :: LCNVFRC   ! Use convective fractions?
   LOGICAL                       :: LLFR      ! Use GEOS-5 flash rates

   ! Arrays
   REAL(dp), POINTER             :: PROFILE(:,:)
   REAL(hp), POINTER             :: SLBASE(:,:,:)
   REAL(sp), POINTER             :: FLASH_DENS_TOT(:,:)
   REAL(sp), POINTER             :: FLASH_DENS_IC(:,:)
   REAL(sp), POINTER             :: FLASH_DENS_CG(:,:)
   REAL(sp), POINTER             :: CONV_DEPTH(:,:)

   ! Overall scale factor to be applied to lightning NOx emissions. Must
   ! be defined in the HEMCO configuration file as extension attribute
   ! 'Scaling_NO'.
   ! SpcScalFldNme is the name of the gridded scale factor. Must be provided
   ! in the HEMCO configuration file as extension attribute 'ScaleField_NO'.
   REAL(sp), ALLOCATABLE          :: SpcScalVal(:)
   CHARACTER(LEN=61), ALLOCATABLE :: SpcScalFldNme(:)

   TYPE(MyInst), POINTER           :: NextInst => NULL()
  END TYPE MyInst

  ! Pointer to all instances
  TYPE(MyInst), POINTER            :: AllInst => NULL()

CONTAINS
!EOC
!------------------------------------------------------------------------------
!                   Harmonized Emissions Component (HEMCO)                    !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: HCOX_LightNOx_Run
!
! !DESCRIPTION: Subroutine HCOX\_LIGHTNOX\_RUN is the driver routine
! to calculate lightning NOx emissions and return them to the HEMCO
! driver routine.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE HCOX_LightNOx_Run( ExtState, HcoState, RC )
!
! !USES:
!
    USE HCO_FluxArr_Mod,  ONLY : HCO_EmisAdd
!
! !INPUT PARAMETERS:
!
    TYPE(Ext_State), POINTER        :: ExtState   ! Module options
    TYPE(HCO_State), POINTER        :: HcoState   ! HEMCO options
!
! !INPUT/OUTPUT PARAMETERS:
!
    INTEGER,         INTENT(INOUT)  :: RC
!
! !REVISION HISTORY:
!  09 Oct 1997 - R. Yantosca - Initial version
!  See https://github.com/geoschem/hemco for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    TYPE(MyInst), POINTER :: Inst
    INTEGER               :: Yr, Mt
    LOGICAL               :: FOUND
    CHARACTER(LEN=255)    :: MSG, LOC

    !=================================================================
    ! HCOX_LIGHTNOX_RUN begins here!
    !=================================================================
    LOC = 'HCOX_LIGHTNOX_RUN (HCOX_LIGHTNOX_MOD.F90)'

    ! Enter
    CALL HCO_ENTER( HcoState%Config%Err, LOC, RC )
    IF ( RC /= HCO_SUCCESS ) THEN
        CALL HCO_ERROR( 'ERROR 0', RC, THISLOC=LOC )
        RETURN
    ENDIF

    ! Return if extension disabled
    IF ( ExtState%LightNOx <= 0 ) THEN
       CALL HCO_LEAVE( HcoState%Config%Err,RC )
       RETURN
    ENDIF

    ! Get pointer to this instance. Varible Inst contains all module
    ! variables for the current instance. The instance number is
    ! ExtState%<yourname>.
    Inst => NULL()
    CALL InstGet ( ExtState%LightNOx, Inst, RC )
    IF ( RC /= HCO_SUCCESS ) THEN
       WRITE(MSG,*) 'Cannot find lightning NOx instance Nr. ', ExtState%LightNOx
       CALL HCO_ERROR(MSG,RC)
       RETURN
    ENDIF

    ! Update lightnox NOx emissions (fill SLBASE)
    CALL LIGHTNOX( HcoState, ExtState, Inst, RC )
    IF ( RC /= HCO_SUCCESS ) THEN
        CALL HCO_ERROR( 'ERROR 1', RC, THISLOC=LOC )
        RETURN
    ENDIF

    !=================================================================
    ! Pass to HEMCO State and update diagnostics
    !=================================================================
    IF ( Inst%IDTNO > 0 ) THEN

       ! Add flux to emission array
       CALL HCO_EmisAdd( HcoState, Inst%SLBASE, Inst%IDTNO, &
                         RC, ExtNr=Inst%ExtNr)
       IF ( RC /= HCO_SUCCESS ) THEN
          CALL HCO_ERROR( 'HCO_EmisAdd error: SLBASE', RC )
          RETURN
       ENDIF

    ENDIF

    ! Return w/ success
    Inst => NULL()
    CALL HCO_LEAVE( HcoState%Config%Err,RC )

  END SUBROUTINE HCOX_LightNOx_Run
!EOC
!------------------------------------------------------------------------------
!                   Harmonized Emissions Component (HEMCO)                    !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: LightNOx
!
! !DESCRIPTION: Subroutine LIGHTNOX uses Price \& Rind's formulation for
!  computing NOx emission from lightning (with various updates).
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE LightNOx( HcoState, ExtState, Inst, RC )
!
! !USES:
!
    USE HCO_Calc_Mod,     ONLY : HCO_EvalFld
    USE HCO_EmisList_Mod, ONLY : HCO_GetPtr
    USE HCO_GeoTools_Mod, ONLY : HCO_LANDTYPE
    USE HCO_Clock_Mod,    ONLY : HcoClock_Get
    USE HCO_Clock_Mod,    ONLY : HcoClock_First
    USE HCO_ExtList_Mod,  ONLY : GetExtOpt
!
! !INPUT PARAMETERS:
!
    TYPE(HCO_State), POINTER        :: HcoState  ! Output obj
    TYPE(Ext_State), POINTER        :: ExtState    ! Module options
    TYPE(MyInst   ), POINTER        :: Inst
!
! !INPUT/OUTPUT PARAMETERS:
!
    INTEGER,         INTENT(INOUT)  :: RC
!
! !REMARKS:
!
! !REVISION HISTORY:
!  10 May 2006 - L. Murray - Initial version
!  See https://github.com/geoschem/hemco for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER             :: I, J, L
    INTEGER             :: LTOP
    INTEGER             :: MONTH
    INTEGER             :: MTYPE
    REAL*8              :: A_M2
    REAL*8              :: A_KM2
    REAL*8              :: H0
    REAL*8              :: IC_CG_RATIO
    REAL*8              :: RATE
    REAL*8              :: RATE_SAVE
    REAL*8              :: TOTAL
    REAL*8              :: TOTAL_CG
    REAL*8              :: TOTAL_IC
    REAL*8              :: X
    REAL*8              :: YMID
    REAL*8              :: XMID
    REAL*8              :: VERTPROF(HcoState%NZ)
    INTEGER             :: LMAX
    INTEGER             :: LNDTYPE
    INTEGER             :: SFCTYPE
    REAL(hp)            :: TROPP
    REAL(dp)            :: TmpScale
    CHARACTER(LEN=255)  :: LOC

    !=================================================================
    ! LIGHTNOX begins here!
    !=================================================================
    LOC = 'LIGHTNOX (HCOX_LIGHTNOX_MOD.F90)'

    ! Enter
    CALL HCO_ENTER( HcoState%Config%Err, LOC, RC )
    IF ( RC /= HCO_SUCCESS ) THEN
        CALL HCO_ERROR( 'ERROR 2', RC, THISLOC=LOC )
        RETURN
    ENDIF

    ! Reset arrays
    Inst%SLBASE         = 0.0_hp
    Inst%FLASH_DENS_TOT = 0.0_sp
    !Inst%FLASH_DENS_IC  = 0.0_sp
    !Inst%FLASH_DENS_CG  = 0.0_sp
    Inst%CONV_DEPTH     = 0.0_sp

    ! LMAX: the highest L-level to look for lightning NOx (usually LLPAR-1)
    LMAX   = HcoState%NZ - 1

    ! Get current month (to be passed to LIGHTDIST)
    CALL HcoClock_Get( HcoState%Clock, cMM=MONTH, RC=RC)
    IF ( RC /= HCO_SUCCESS ) THEN
        CALL HCO_ERROR( 'ERROR 3', RC, THISLOC=LOC )
        RETURN
    ENDIF

    !=================================================================
    ! Compute lightning NOx emissions for each (I,J) column
    !=================================================================

!$OMP PARALLEL DO                                                     &
!$OMP DEFAULT( SHARED )                                               &
!$OMP PRIVATE( I,         J,           L,        A_M2,        A_KM2  ) &
!$OMP PRIVATE( YMID,      XMID,        LTOP,     MTYPE               ) &
!$OMP PRIVATE( LNDTYPE,   SFCTYPE,     TROPP                         ) &
!$OMP PRIVATE( RATE,      RATE_SAVE,   H0,       IC_CG_RATIO         ) &
!$OMP PRIVATE( TOTAL,     TOTAL_IC,    TOTAL_CG, VERTPROF,    X      ) &
!$OMP SCHEDULE( DYNAMIC )

    ! Loop over surface boxes
    DO J = 1, HcoState%NY
    DO I = 1, HcoState%NX

       ! Grid box surface areas in [m2] and [km2]
       A_M2     = HcoState%Grid%AREA_M2%Val( I, J )
       A_KM2    = A_M2 / 1d6

       ! Grid box latitude and longitude [degrees]
       YMID     = HcoState%Grid%YMID%Val( I, J )
       XMID     = HcoState%Grid%XMID%Val( I, J )

       ! Make sure xmid is between -180 and +180
       IF ( XMID >= 180.0d0 ) XMID = XMID - 360.0d0

       ! Get surface type. Note that these types are different than
       ! the types used elsewhere in this module. HCO_LANDTYPE returns
       ! 0=ocean, 1=land, 2=ice; elsewhere we use 0=land, 1=ocean, 2=ice!
       LNDTYPE = HCO_LANDTYPE( ExtState%FRLAND%Arr%Val(I,J),   &
                               ExtState%FRLANDIC%Arr%Val(I,J), &
                               ExtState%FROCEAN%Arr%Val(I,J),  &
                               ExtState%FRSEAICE%Arr%Val(I,J), &
                               ExtState%FRLAKE%Arr%Val(I,J))

       ! Set surface type (0=land, 1=ocean, 2=ice) based on land type
       IF ( LNDTYPE == 2 ) THEN
          SFCTYPE = 2    ! Ice
       ELSEIF ( LNDTYPE == 1 ) THEN
          SFCTYPE = 0    ! Land
       ELSE
          SFCTYPE = 1    ! Ocean (default)
       ENDIF

       ! Tropopause pressure. Convert to Pa
       TROPP = ExtState%TROPP%Arr%Val(I,J) !* 100.0_hp

       !===========================================================
       ! Initialize
       !===========================================================
       RATE          = 0.0
       RATE_SAVE     = 0.0
       H0            = 0.0
       IC_CG_RATIO   = 1.0

       TOTAL         = 0d0
       TOTAL_IC      = 0d0
       TOTAL_CG      = 0d0

       !===========================================================
       ! (1) Get flash density [#/km2/s] from meteorology
       !===========================================================

       !-----------------------------------------------------------
       ! (1a) Prescribed in HEMCO_Config.rc
       !-----------------------------------------------------------
       RATE        = ExtState%FLASH_DENS%Arr%Val(I,J)

       !-----------------------------------------------------------
       ! (1b) From GEOS-5
       !-----------------------------------------------------------
       IF ( Inst%LLFR                                  &
            .AND. ASSOCIATED( ExtState%LFR%Arr%Val   ) &
            .AND. ASSOCIATED( ExtState%BYNCY%Arr%Val )  ) THEN
          RATE     = ExtState%LFR%Arr%Val(I,J)
       ENDIF

       ! Error check: do not continue if flash rate is zero
       IF ( RATE <= 0.0 ) CYCLE

       !===========================================================
       ! (2) Get depth of convection [m] and find associated LTOP
       !===========================================================

       !-----------------------------------------------------------
       ! (2a) Prescribed in HEMCO_Config.rc
       !-----------------------------------------------------------
       H0          = ExtState%CONV_DEPTH%Arr%Val(I,J)
       LTOP = 1
       DO L = 1, HcoState%NZ
          IF ( SUM(HcoState%Grid%BXHEIGHT_M%Val(I,J,1:L)) > H0 ) THEN
             LTOP = L
             EXIT
          ENDIF
       ENDDO
       ! Reset H0 to be the height of that layer in the model,
       ! to avoid negative values in the partitioning
       H0          = SUM(HcoState%Grid%BXHEIGHT_M%Val(I,J,1:LTOP))

       !-----------------------------------------------------------
       ! (2b) From GEOS-5
       !-----------------------------------------------------------
       IF ( Inst%LLFR                                  &
            .AND. ASSOCIATED( ExtState%LFR%Arr%Val   ) &
            .AND. ASSOCIATED( ExtState%BYNCY%Arr%Val )  ) THEN

          ! Set LTOP to top of buoyancy
          DO L = HcoState%NZ, 1, -1
             IF ( ExtState%BYNCY%Arr%Val(I,J,L) >= 0.0_sp ) THEN
                LTOP = L + 1
                EXIT
             ENDIF
          ENDDO
          !LTOP = MAX( LTOP, LMAX )
          ! H0 is the convective cloud top height [m].  This is the
          ! distance from the surface to the top edge of box (I,J,LTOP).
          H0 = SUM(HcoState%Grid%BXHEIGHT_M%Val(I,J,1:LTOP))
       ENDIF

       ! Save out convective cloud depth
       Inst%CONV_DEPTH(I,J) = LTOP

       !===========================================================
       ! (3) Compute ratio of CG vs total flashes
       !===========================================================

       ! Ratio of cloud-to-ground flashes to total # of flashes
       X    = 1d0 / ( 1d0 + IC_CG_RATIO )

       !-----------------------------------------------------------
       ! Store flash rates [flashes/km2/min]
       !-----------------------------------------------------------
       IF ( RATE > 0d0 ) THEN

          ! Flashes per km2 per minute
          RATE_SAVE   = RATE / 60d0

          ! Store total, IC, and CG flash rates
          Inst%FLASH_DENS_TOT(I,J) = RATE_SAVE
          !Inst%FLASH_DENS_IC(I,J)  = RATE_SAVE * X
          !Inst%FLASH_DENS_CG(I,J)  = RATE_SAVE * ( 1d0 - X )

       ENDIF

       !===========================================================
       ! (4) Compute LNOx yield for IC and CG flashes (molec/km2/s)
       !===========================================================

       ! Compute LNOx emissions for tropics or midlats
       IF ( XMID > EAST_WEST_DIV ) THEN

          !--------------------------------------------------------
          ! (4a) We are in EURASIA
          !--------------------------------------------------------
          IF ( YMID > EAST_NS_DIV ) THEN

             ! Eurasian Mid-Latitudes
             TOTAL_IC = RFLASH_MIDLAT * RATE * ( 1d0 - X )
             TOTAL_CG = RFLASH_MIDLAT * RATE * X

          ELSE

             ! Eurasian Tropics
             TOTAL_IC = RFLASH_TROPIC * RATE * ( 1d0 - X )
             TOTAL_CG = RFLASH_TROPIC * RATE * X

          ENDIF

       ELSE

          !--------------------------------------------------------
          ! (4b) We are in the AMERICAS
          !--------------------------------------------------------
          IF ( YMID > WEST_NS_DIV ) THEN

             ! American Mid-Latitudes
             TOTAL_IC = RFLASH_MIDLAT * RATE * ( 1d0 - X )
             TOTAL_CG = RFLASH_MIDLAT * RATE * X

          ELSE

             ! American Tropics
             TOTAL_IC = RFLASH_TROPIC * RATE * ( 1d0 - X )
             TOTAL_CG = RFLASH_TROPIC * RATE * X

          ENDIF
       ENDIF

       !===========================================================
       ! (5) Compute column total lightning NOx yield (kg/km2/s)
       !===========================================================

       ! Sum of IC + CG
       TOTAL = TOTAL_IC + TOTAL_CG

       ! Convert from molec km-2 s-1 to kg(NO) m-2 s-1
       TOTAL = TOTAL * ( HcoState%Spc(Inst%IDTNO)%MW_g / 1000.0_hp ) / &
                         HcoState%Phys%Avgdr / 1000000.0_hp

       !===========================================================
       ! (6) Distribute column LNOx vertically from surface to LTOP
       !===========================================================

       !-----------------------------------------------------------
       ! LIGHTDIST computes the lightning NOx distribution from
       ! the ground to the convective cloud top using cumulative
       ! distribution functions for ocean flashes, tropical land
       ! flashes, and non-tropical land flashes, as specified by
       ! Lesley Ott [JGR, 2010]
       !-----------------------------------------------------------

       ! If there's lightning NOx w/in the column ...
       IF ( TOTAL > 0d0 ) THEN

          ! Partition the column total NOx [kg/m2/s] from lightning
          ! into the vertical using Ott et al. [2010] PDF functions
          CALL LIGHTDIST( I, J, LTOP, H0, YMID, TOTAL, VERTPROF, &
                          ExtState, HcoState, SFCTYPE, MONTH, MTYPE, Inst )

          ! Add vertically partitioned NOx into SLBASE array
          DO L = 1, HcoState%NZ
             Inst%SLBASE(I,J,L) = VERTPROF(L)

             ! No lightning NOx emissions in the stratosphere (cdh, 4/25/2013)
             IF ( HcoState%Grid%PEDGE%Val(I,J,L) < TROPP ) THEN
                Inst%SLBASE(I,J,L) = 0.0_hp
             ENDIF

          ENDDO
       ENDIF

    ENDDO
    ENDDO
!$OMP END PARALLEL DO

    !-----------------------------------------------------------------
    ! Eventually add scale factors
    !-----------------------------------------------------------------

    ! Eventually apply species specific scale factor
    IF ( Inst%SpcScalVal(1) /= 1.0_sp ) THEN
       Inst%SLBASE = Inst%SLBASE * Inst%SpcScalVal(1)
    ENDIF

    ! Eventually apply spatiotemporal scale factors
    CALL HCOX_SCALE( HcoState, Inst%SLBASE, TRIM(Inst%SpcScalFldNme(1)), RC )
    IF ( RC /= HCO_SUCCESS ) THEN
        CALL HCO_ERROR( 'ERROR 4', RC, THISLOC=LOC )
        RETURN
    ENDIF

    ! Return w/ success
    CALL HCO_LEAVE( HcoState%Config%Err,RC )

  END SUBROUTINE LightNOx
!EOC
!------------------------------------------------------------------------------
!                   Harmonized Emissions Component (HEMCO)                    !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: LightDist
!
! !DESCRIPTION: Subroutine LightDist reads in the CDF used to partition the
!  column lightning NOx into the GEOS-Chem vertical layers.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE LightDist( I, J, LTOP, H0, XLAT, TOTAL, VERTPROF, &
                        ExtState, HcoState, SFCTYPE, MONTH, MTYPE, Inst )
!
! !INPUT PARAMETERS:
!
    INTEGER,         INTENT(IN)  :: I          ! Longitude index
    INTEGER,         INTENT(IN)  :: J          ! Latitude index
    INTEGER,         INTENT(IN)  :: LTOP       ! Level of conv cloud top
    REAL*8,          INTENT(IN)  :: H0         ! Conv cloud top height [m]
    REAL*8,          INTENT(IN)  :: XLAT       ! Latitude value [degrees]
    REAL*8,          INTENT(IN)  :: TOTAL      ! Column Total # of LNOx molec
    TYPE(Ext_State), POINTER     :: ExtState   ! Module options
    TYPE(HCO_State), POINTER     :: HcoState   ! Hemco state object
    INTEGER,         INTENT(IN)  :: SFCTYPE    ! Surface type
    INTEGER,         INTENT(IN)  :: MONTH      ! Current month
    TYPE(MyInst),    POINTER     :: Inst       ! Hemco state object
!
! !OUTPUT PARAMETERS:
!
    REAL*8,          INTENT(OUT) :: VERTPROF(HcoState%NZ) ! Vertical profile
    INTEGER,         INTENT(OUT) :: MTYPE                 ! lightning type
!
! !REMARKS:
!  References:
!  ============================================================================
!  (1 ) Pickering et al., JGR 103, 31,203 - 31,316, 1998.
!  (2 ) Ott et al., JGR, 2010
!  (3 ) Allen et al., JGR, 2010
!
! !REVISION HISTORY:
!  18 Sep 2002 - M. Evans - Initial version (based on Yuhang Wang's code)
!  See https://github.com/geoschem/hemco for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER :: L
    REAL*8  :: ZHEIGHT, YMID
    REAL*8  :: FRAC(HcoState%NZ)

      !=================================================================
      ! LIGHTDIST begins here!
      !=================================================================

    ! Initialize
    MTYPE    = 0
    VERTPROF = 0d0

    !%%% NOTE: Use L=1 for GRID_MOD functions.  This is OK for the
    !%%% existing GEOS-Chem with a pure cartesian grid, but may be an
    !%%% issue when interfaced with a GCM with a non-regular grid
    !%%% (bmy, 3/1/12)
    YMID     = HcoState%Grid%YMID%Val( I, J )

    !=================================================================
    ! Test whether location (I,J) is continental, marine, or snow/ice
    !
    ! Depending on the combination of land/water and latitude,
    ! assign a flag describing the type of lightning:
    !
    !   MTYPE = 1: ocean lightning
    !   MTYPE = 2: tropical continental lightning
    !   MTYPE = 3: midlatitude continental lightning
    !   MTYPE = 4: subtropical lightning
    !
    ! (ltm, bmy, 1/25/11)
    !=================================================================

    ! Assign profile kind to grid box, following Allen et al.
    ! [JGR, 2010] (ltm, 1/25,11)

    SELECT CASE ( MONTH )

       ! Southern Hemisphere Summer
       CASE ( 1,2,3,12 )

           IF ( ABS(YMID) .le. 15 ) THEN
              IF ( SFCTYPE == 0 ) THEN
                 MTYPE = 2        ! Tropical continental
              ELSE
                 MTYPE = 1        ! Tropical marine
              ENDIF
           ELSE IF ( ( YMID .gt. 15. ) .and. ( YMID .le. 30. ) ) THEN
              MTYPE = 4           ! N. Subtropics
           ELSE IF ( ( YMID .ge. -40. ) .and. ( YMID .lt. -15. ) ) THEN
              MTYPE = 4           ! S. Subtropics
           ELSE
              MTYPE = 3           ! Midlatitude
           ENDIF

        ! Equinox months
        CASE ( 4,5,10,11 )

           IF ( ABS(YMID) .le. 15 ) THEN
              IF ( SFCTYPE == 0 ) THEN
                 MTYPE = 2        ! Tropical continental
              ELSE
                 MTYPE = 1        ! Tropical marine
              ENDIF
           ELSE IF ( ABS(YMID) .le. 30 ) THEN
              MTYPE = 4           ! Subtropics
           ELSE
              MTYPE = 3           ! Midlatitude
           ENDIF

        ! Northern Hemisphere Summer
        CASE ( 6,7,8,9 )

           IF ( ABS(YMID) .le. 15 ) THEN
              IF ( SFCTYPE == 0 ) THEN
                 MTYPE = 2        ! Tropical continental
              ELSE
                 MTYPE = 1        ! Tropical marine
              ENDIF
           ELSE IF ( ( YMID .gt. 15. ) .and. ( YMID .le. 40. ) ) THEN
              MTYPE = 4           ! N. Subtropics
           ELSE IF ( ( YMID .ge. -30. ) .and. ( YMID .lt. -15. ) ) THEN
              MTYPE = 4           ! S. Subtropics
           ELSE
              MTYPE = 3           ! Midlatitude
           ENDIF

    END SELECT

    ! Extra safety check for pathological grid boxes (bmy, 11/29/06)
    IF ( MTYPE == 0 ) RETURN

    !=================================================================
    ! Use the CDF for this type of lightning to partition the total
    ! column lightning NOx into the layers
    !=================================================================
    ZHEIGHT = 0.0

    ! Compute the height [km] at the top of each vertical level.
    ! Look up the cumulative fraction of NOx for each vertical level
    DO L = 1, LTOP
       ZHEIGHT = ZHEIGHT + HcoState%Grid%BXHEIGHT_M%Val(I,J,L)
       FRAC(L) = Inst%PROFILE( NINT( ( ZHEIGHT/H0 )*3200. ), MTYPE ) *0.01
    ENDDO

    ! Convert from cumulative fraction to fraction for each level
    DO L = LTOP, 2, - 1
       FRAC(L) = FRAC(L) - FRAC(L-1)
    ENDDO

    ! Partition lightning NOx by layer into VERTPROF
    DO L = 1, LTOP
       VERTPROF(L) = ( FRAC(L) * TOTAL )
    ENDDO

  END SUBROUTINE LightDist
!EOC
!------------------------------------------------------------------------------
!                   Harmonized Emissions Component (HEMCO)                    !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: HCOX_LightNOx_Init
!
! !DESCRIPTION: Subroutine HCOX\_LIGHTNOX\_INIT allocates all module arrays.
!  It also reads the lightning CDF data from disk before the first lightning
!  timestep.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE HCOX_LightNOx_Init( HcoState, ExtName, ExtState, RC )
!
! !USES:
!
    USE HCO_Chartools_Mod, ONLY : HCO_CharParse
    USE HCO_ExtList_Mod,   ONLY : GetExtNr
    USE HCO_ExtList_Mod,   ONLY : GetExtOpt
    USE HCO_ExtList_Mod,   ONLY : GetExtSpcVal
    USE HCO_State_Mod,     ONLY : HCO_GetHcoID
    USE HCO_State_Mod,     ONLY : HCO_GetExtHcoID
    USE HCO_ReadList_Mod,  ONLY : ReadList_Remove
    USE HCO_inquireMod,    ONLY : findfreeLUN
!
! !INPUT PARAMETERS:
!
    TYPE(HCO_State),  POINTER        :: HcoState   ! Hemco options
    CHARACTER(LEN=*), INTENT(IN   )  :: ExtName    ! Extension name
    TYPE(Ext_State),  POINTER        :: ExtState     ! Module options
!
! !OUTPUT PARAMETERS:
!
    INTEGER,          INTENT(  OUT)  :: RC
!
! !REVISION HISTORY:
!  14 Apr 2004 - R. Yantosca - Initial version
!  See https://github.com/geoschem/hemco for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER                        :: AS, III, IOS, JJJ, IU_FILE, nSpc
    INTEGER                        :: ExtNr
    LOGICAL                        :: FOUND, FileExists
    INTEGER, ALLOCATABLE           :: HcoIDs(:)
    CHARACTER(LEN=31), ALLOCATABLE :: SpcNames(:)
    CHARACTER(LEN=255)             :: MSG, LOC, FILENAME, FileMsg
    TYPE(MyInst), POINTER          :: Inst

    !=======================================================================
    ! HCOX_LightNOX_Init begins here!
    !=======================================================================
    LOC = 'HCOX_LightNOX_Init (HCOX_LIGHTNOX_MOD.F90)'

    ! Extension Nr.
    ExtNr = GetExtNr( HcoState%Config%ExtList, TRIM(ExtName) )
    IF ( ExtNr <= 0 ) RETURN

    ! Enter
    CALL HCO_ENTER( HcoState%Config%Err, LOC, RC )
    IF ( RC /= HCO_SUCCESS ) THEN
        CALL HCO_ERROR( 'ERROR 5', RC, THISLOC=LOC )
        RETURN
    ENDIF

    ! Create LightNOx instance for this simulation
    Inst => NULL()
    CALL InstCreate ( ExtNr, ExtState%LightNOx, Inst, RC )
    IF ( RC /= HCO_SUCCESS ) THEN
       CALL HCO_ERROR ( 'Cannot create LightNOx instance', RC )
       RETURN
    ENDIF

    !=======================================================================
    ! Obtain lightning CDF's from Ott et al [JGR, 2010].
    !
    ! PART 1 --- Move the file name check to the front of this routine to
    ! facilitate the GEOS-Chem dry-run and HEMCO-standalone dry-run.
    !=======================================================================

    ! Get filename from configuration file
    CALL GetExtOpt( HcoState%Config, ExtNr, 'CDF table',                     &
                    OptValChar=FILENAME, RC=RC                              )
    IF ( RC /= HCO_SUCCESS ) THEN
        CALL HCO_ERROR( 'ERROR 6', RC, THISLOC=LOC )
        RETURN
    ENDIF

    ! Call HEMCO parser to replace tokens such as $ROOT, $MET, or $RES.
    ! There shouldn't be any date token in there ($YYYY, etc.), so just
    ! provide some dummy variables here
    CALL HCO_CharParse( HcoState%Config, FILENAME, -999, -1, -1, -1, -1, RC )
    IF ( RC /= HCO_SUCCESS ) THEN
        CALL HCO_ERROR( 'ERROR 7', RC, THISLOC=LOC )
        RETURN
    ENDIF

    !-----------------------------------------------------------------------
    ! In dry-run mode, print file path to dryrun log and exit.
    ! Otherwise, print file path to the HEMCO log file and continue.
    !-----------------------------------------------------------------------

    ! Test if the file exists
    INQUIRE( FILE=TRIM( FileName ), EXIST=FileExists )

    ! Create a display string based on whether or not the file is found
    IF ( FileExists ) THEN
       FileMsg = 'HEMCO (LIGHTNOX): Opening'
    ELSE
       FileMsg = 'HEMCO (LIGHTNOX): REQUIRED FILE NOT FOUND'
    ENDIF

    ! Write file status to stdout and the HEMCO log
    IF ( HcoState%amIRoot ) THEN
       WRITE( 6,   300 ) TRIM( FileMsg ), TRIM( FileName )
       WRITE( MSG, 300 ) TRIM( FileMsg ), TRIM( FileName )
       CALL HCO_MSG( HcoState%Config%Err, MSG )
 300   FORMAT( a, ' ', a )
    ENDIF

    ! For dry-run simulation, return to calling program.
    ! For regular simulations, throw an error if we can't find the file.
    IF ( HcoState%Options%IsDryRun ) THEN
       RETURN
    ELSE
       IF ( .not. FileExists ) THEN
          WRITE( MSG, 300 ) TRIM( FileMsg ), TRIM( FileName )
          CALL HCO_ERROR(MSG, RC )
          RETURN
       ENDIF
    ENDIF

    !=======================================================================
    ! Exit if this is a GEOS-Chem or HEMCO-standalone dry-run
    !=======================================================================
    IF ( HcoState%Options%IsDryRun ) THEN
       Inst => NULL()
       CALL HCO_LEAVE( HcoState%Config%Err,RC )
       RETURN
    ENDIF

    !=======================================================================
    ! Continue for regular simulations ...
    !=======================================================================

    ! Check for usage of convective fractions. This becomes only active
    ! if both the convective fraction and the buoyancy field are available.
    CALL GetExtOpt( HcoState%Config, ExtNr, 'Use CNV_FRC', &
                     OptValBool=Inst%LCNVFRC, FOUND=FOUND, RC=RC )
    IF ( RC /= HCO_SUCCESS ) THEN
        CALL HCO_ERROR( 'ERROR 8', RC, THISLOC=LOC )
        RETURN
    ENDIF
    IF ( .NOT. FOUND ) Inst%LCNVFRC = .FALSE.

    ! Check for usage of GEOS-5 lightning flash rates. If on, the GEOS-5
    ! flash rates (where available) are used instead of the computed flash
    ! rates. This is off by default.
    CALL GetExtOpt( HcoState%Config, ExtNr, 'GEOS-5 flash rates', &
                     OptValBool=Inst%LLFR, FOUND=FOUND, RC=RC )
    IF ( RC /= HCO_SUCCESS ) THEN
        CALL HCO_ERROR( 'ERROR 9', RC, THISLOC=LOC )
        RETURN
    ENDIF
    IF ( .NOT. FOUND ) Inst%LLFR = .FALSE.

    ! Get species ID
    CALL HCO_GetExtHcoID( HcoState, ExtNr, HcoIDs, SpcNames, nSpc, RC )
    IF ( RC /= HCO_SUCCESS ) THEN
        CALL HCO_ERROR( 'ERROR 10', RC, THISLOC=LOC )
        RETURN
    ENDIF
    IF ( nSpc /= 1 ) THEN
       MSG = 'Lightning NOx module must have exactly one species!'
       CALL HCO_ERROR(MSG, RC )
       RETURN
    ENDIF
    Inst%IDTNO = HcoIDs(1)

    ! Get species scale factor
    CALL GetExtSpcVal( HcoState%Config, ExtNr, nSpc, &
                       SpcNames, 'Scaling', 1.0_sp, Inst%SpcScalVal, RC )
    IF ( RC /= HCO_SUCCESS ) THEN
        CALL HCO_ERROR( 'ERROR 11', RC, THISLOC=LOC )
        RETURN
    ENDIF

    CALL GetExtSpcVal( HcoState%Config, ExtNr, nSpc, &
                       SpcNames, 'ScaleField', HCOX_NOSCALE, Inst%SpcScalFldNme, RC )
    IF ( RC /= HCO_SUCCESS ) THEN
        CALL HCO_ERROR( 'ERROR 12', RC, THISLOC=LOC )
        RETURN
    ENDIF

    ! Echo info about this extension
    IF ( HcoState%amIRoot ) THEN
       MSG = 'Use lightning NOx emissions (extension module)'
       CALL HCO_MSG(HcoState%Config%Err,MSG, SEP1='-' )
       WRITE(MSG,*) ' - Use species ', TRIM(SpcNames(1)), '->', Inst%IDTNO
       CALL HCO_MSG(HcoState%Config%Err,MSG)
       WRITE(MSG,*) ' - Use GEOS-5 flash rates: ', Inst%LLFR
       CALL HCO_MSG(HcoState%Config%Err,MSG)
       WRITE(MSG,*) ' - Use scalar scale factor: ', Inst%SpcScalVal(1)
       CALL HCO_MSG(HcoState%Config%Err,MSG)
       WRITE(MSG,*) ' - Use gridded scale field: ', TRIM(Inst%SpcScalFldNme(1))
       CALL HCO_MSG(HcoState%Config%Err,MSG)
    ENDIF

    !=======================================================================
    ! Allocate arrays
    !=======================================================================

    ALLOCATE( Inst%PROFILE( NNLIGHT, NLTYPE ), STAT=AS )
    IF( AS /= 0 ) THEN
       CALL HCO_ERROR ( 'PROFILE', RC )
       RETURN
    ENDIF
    Inst%PROFILE = 0.0_hp

    ALLOCATE( Inst%SLBASE(HcoState%NX,HcoState%NY,HcoState%NZ), STAT=AS )
    IF( AS /= 0 ) THEN
       CALL HCO_ERROR ( 'SLBASE', RC )
       RETURN
    ENDIF
    Inst%SLBASE = 0.0_hp

    ALLOCATE ( Inst%FLASH_DENS_TOT( HcoState%NX, HcoState%NY), STAT=AS )
    IF ( AS/=0 ) THEN
       CALL HCO_ERROR( 'FLASH_DENS_TOT', RC )
       RETURN
    ENDIF
    Inst%FLASH_DENS_TOT = 0.0_sp

    !ALLOCATE ( Inst%FLASH_DENS_IC( HcoState%NX, HcoState%NY), STAT=AS )
    !IF ( AS/=0 ) THEN
    !   CALL HCO_ERROR( 'FLASH_DENS_IC', RC )
    !   RETURN
    !ENDIF
    !Inst%FLASH_DENS_IC = 0.0_sp

    !ALLOCATE ( Inst%FLASH_DENS_CG( HcoState%NX, HcoState%NY), STAT=AS )
    !IF ( AS/=0 ) THEN
    !   CALL HCO_ERROR( 'FLASH_DENS_CG', RC )
    !   RETURN
    !ENDIF
    !Inst%FLASH_DENS_CG = 0.0_sp

    ALLOCATE ( Inst%CONV_DEPTH( HcoState%NX, HcoState%NY), STAT=AS )
    IF ( AS/=0 ) THEN
       CALL HCO_ERROR( 'CONV_DEPTH', RC )
       RETURN
    ENDIF
    Inst%CONV_DEPTH = 0.0_sp

    !=======================================================================
    ! Obtain lightning CDF's from Ott et al [JGR, 2010].
    !
    ! PART 2 --- Read the data!
    !=======================================================================

    ! Find a free file LUN
    IU_FILE = findFreeLUN()

    ! Open file containing lightning NOx PDF data
    OPEN( IU_FILE, FILE=TRIM( FILENAME ), STATUS='OLD', IOSTAT=IOS )
    IF ( IOS /= 0 ) THEN
       MSG = 'IOERROR: LightDist: 1'
       CALL HCO_ERROR(MSG, RC )
       RETURN
    ENDIF

    ! Read 12 header lines
    DO III = 1, 12
       READ( IU_FILE, '(a)', IOSTAT=IOS )
       IF ( IOS /= 0 ) THEN
          MSG = 'IOERROR: LightDist: 2'
          CALL HCO_ERROR(MSG, RC )
          RETURN
       ENDIF
    ENDDO

    ! Read NNLIGHT types of lightning profiles
    DO III = 1, NNLIGHT
       READ( IU_FILE,*,IOSTAT=IOS) (Inst%PROFILE(III,JJJ),JJJ=1,NLTYPE)
       IF ( IOS /= 0 ) THEN
          MSG = 'IOERROR: LightDist: 3'
          CALL HCO_ERROR(MSG, RC )
          RETURN
       ENDIF
    ENDDO

    ! Close file
    CLOSE( IU_FILE )

    !=======================================================================
    ! Create diagnostics for lightning flash rates and convective cloud height
    !=======================================================================
    CALL Diagn_Create( HcoState  = HcoState,                        &
                       cName     = 'HcoLightningFlashRate_Total',   &
                       ExtNr     = ExtNr,                           &
                       Cat       = -1,                              &
                       Hier      = -1,                              &
                       HcoID     = -1,                              &
                       SpaceDim  = 2,                               &
                       OutUnit   = 'flashes/min/km2',               &
                       AutoFill  = 0,                               &
                       Trgt2D    = Inst%FLASH_DENS_TOT,             &
                       RC        = RC )
    IF ( RC /= HCO_SUCCESS ) THEN
        CALL HCO_ERROR( 'ERROR 13', RC, THISLOC=LOC )
        RETURN
    ENDIF

    !CALL Diagn_Create( HcoState  = HcoState,                        &
    !                   cName     = 'HcoLightningFlashRate_IntraCld', &
    !                   ExtNr     = ExtNr,                           &
    !                   Cat       = -1,                              &
    !                   Hier      = -1,                              &
    !                   HcoID     = -1,                              &
    !                   SpaceDim  = 2,                               &
    !                   OutUnit   = 'flashes/min/km2',               &
    !                   AutoFill  = 0,                               &
    !                   Trgt2D    = Inst%FLASH_DENS_IC,              &
    !                   RC        = RC )
    !IF ( RC /= HCO_SUCCESS ) RETURN
    !    
    !CALL Diagn_Create( HcoState  = HcoState,                         &
    !                   cName     = 'HcoLightningFlashRate_CldGround', &
    !                   ExtNr     = ExtNr,                            &
    !                   Cat       = -1,                               &
    !                   Hier      = -1,                               &
    !                   HcoID     = -1,                               &
    !                   SpaceDim  = 2,                                &
    !                   OutUnit   = 'flashes/min/km2',                &
    !                   AutoFill  = 0,                                &
    !                   Trgt2D    = Inst%FLASH_DENS_CG,               &
    !                   RC        = RC )
    !IF ( RC /= HCO_SUCCESS ) RETURN

    CALL Diagn_Create( HcoState  = HcoState,                         &
                       cName     = 'HcoConvectiveCloudTopHeight',       &
                       ExtNr     = ExtNr,                            &
                       Cat       = -1,                               &
                       Hier      = -1,                               &
                       HcoID     = -1,                               &
                       SpaceDim  = 2,                                &
                       OutUnit   = '1',                              &
                       AutoFill  = 0,                                &
                       Trgt2D    = Inst%CONV_DEPTH,                  &
                       RC        = RC )
    IF ( RC /= HCO_SUCCESS ) THEN
        CALL HCO_ERROR( 'ERROR 14', RC, THISLOC=LOC )
        RETURN
    ENDIF

    !=======================================================================
    ! Activate met fields required by this module
    !=======================================================================
    ExtState%TK%DoUse         = .TRUE.
    ExtState%TROPP%DoUse      = .TRUE.
    ExtState%CNV_MFC%DoUse    = .TRUE.
    ExtState%CNV_FRC%DoUse    = .TRUE.
    ExtState%LFR%DoUse        = .TRUE.
    ExtState%FLASH_DENS%DoUse = .TRUE.
    ExtState%CONV_DEPTH%DoUse = .TRUE.
    ExtState%FRLAND%DoUse     = .TRUE.
    ExtState%FRLANDIC%DoUse   = .TRUE.
    ExtState%FROCEAN%DoUse    = .TRUE.
    ExtState%FRSEAICE%DoUse   = .TRUE.
    ExtState%FRLAKE%DoUse     = .TRUE.

    ! Only activate BYNCY and LFR if they are needed
    IF ( Inst%LCNVFRC .OR. Inst%LLFR ) ExtState%BYNCY%DoUse = .TRUE.
    IF ( Inst%LLFR ) ExtState%LFR%DoUse = .TRUE.

    ! Cleanup
    Inst => NULL()

    ! Leave w/ success
    IF ( ALLOCATED(HcoIDs  ) ) DEALLOCATE(HcoIDs  )
    IF ( ALLOCATED(SpcNames) ) DEALLOCATE(SpcNames)
    CALL HCO_LEAVE( HcoState%Config%Err,RC )

  END SUBROUTINE HCOX_LightNOx_Init
!EOC
!------------------------------------------------------------------------------
!                   Harmonized Emissions Component (HEMCO)                    !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: hcox_lightnox_final
!
! !DESCRIPTION: Subroutine HCOX\_LIGHTNOX\_FINAL deallocates all module
!  arrays.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE HCOX_LightNOx_Final( ExtState )
!
! !INPUT PARAMETERS:
!
    TYPE(Ext_State),  POINTER       :: ExtState   ! Module options
!
! !REVISION HISTORY:
!  14 Apr 2004 - R. Yantosca - Initial version
!  See https://github.com/geoschem/hemco for complete history
!EOP
!------------------------------------------------------------------------------
!BOC

    !=================================================================
    ! Cleanup module arrays
    !=================================================================
    CALL InstRemove ( ExtState%LightNOx )

  END SUBROUTINE HCOX_LightNOx_Final
!EOC
!------------------------------------------------------------------------------
!                   Harmonized Emissions Component (HEMCO)                    !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: InstGet
!
! !DESCRIPTION: Subroutine InstGet returns a pointer to the desired instance.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE InstGet ( Instance, Inst, RC, PrevInst )
!
! !INPUT PARAMETERS:
!
    INTEGER                             :: Instance
    TYPE(MyInst),     POINTER           :: Inst
    INTEGER                             :: RC
    TYPE(MyInst),     POINTER, OPTIONAL :: PrevInst
!
! !REVISION HISTORY:
!  18 Feb 2016 - C. Keller   - Initial version
!  See https://github.com/geoschem/hemco for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
    TYPE(MyInst),     POINTER    :: PrvInst

    !=================================================================
    ! InstGet begins here!
    !=================================================================

    ! Get instance. Also archive previous instance.
    PrvInst => NULL()
    Inst    => AllInst
    DO WHILE ( ASSOCIATED(Inst) )
       IF ( Inst%Instance == Instance ) EXIT
       PrvInst => Inst
       Inst    => Inst%NextInst
    END DO
    IF ( .NOT. ASSOCIATED( Inst ) ) THEN
       RC = HCO_FAIL
       RETURN
    ENDIF

    ! Pass output arguments
    IF ( PRESENT(PrevInst) ) PrevInst => PrvInst

    ! Cleanup & Return
    PrvInst => NULL()
    RC = HCO_SUCCESS

  END SUBROUTINE InstGet
!EOC
!------------------------------------------------------------------------------
!                   Harmonized Emissions Component (HEMCO)                    !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: InstCreate
!
! !DESCRIPTION: Subroutine InstCreate adds a new instance to the list of
!  instances, assigns a unique instance number to this new instance, and
!  archives this instance number to output argument Instance.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE InstCreate ( ExtNr, Instance, Inst, RC )
!
! !INPUT PARAMETERS:
!
    INTEGER,       INTENT(IN)       :: ExtNr
!
! !OUTPUT PARAMETERS:
!
    INTEGER,       INTENT(  OUT)    :: Instance
    TYPE(MyInst),  POINTER          :: Inst
!
! !INPUT/OUTPUT PARAMETERS:
!
    INTEGER,       INTENT(INOUT)    :: RC
!
! !REVISION HISTORY:
!  18 Feb 2016 - C. Keller   - Initial version
!  See https://github.com/geoschem/hemco for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
    TYPE(MyInst), POINTER          :: TmpInst
    INTEGER                        :: nnInst

    !=================================================================
    ! InstCreate begins here!
    !=================================================================

    ! ----------------------------------------------------------------
    ! Generic instance initialization
    ! ----------------------------------------------------------------

    ! Initialize
    Inst => NULL()

    ! Get number of already existing instances
    TmpInst => AllInst
    nnInst = 0
    DO WHILE ( ASSOCIATED(TmpInst) )
       nnInst  =  nnInst + 1
       TmpInst => TmpInst%NextInst
    END DO

    ! Create new instance
    ALLOCATE(Inst)
    Inst%Instance = nnInst + 1
    Inst%ExtNr    = ExtNr

    ! Attach to instance list
    Inst%NextInst => AllInst
    AllInst       => Inst

    ! Update output instance
    Instance = Inst%Instance

    ! Return w/ success
    RC = HCO_SUCCESS

  END SUBROUTINE InstCreate
!EOC
!------------------------------------------------------------------------------
!                   Harmonized Emissions Component (HEMCO)                    !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: InstRemove
!
! !DESCRIPTION: Subroutine InstRemove removes an instance from the list of
! instances.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE InstRemove ( Instance )
!
! !INPUT PARAMETERS:
!
    INTEGER                         :: Instance
!
! !REVISION HISTORY:
!  18 Feb 2016 - C. Keller   - Initial version
!  See https://github.com/geoschem/hemco for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
    INTEGER                     :: RC
    TYPE(MyInst), POINTER       :: PrevInst => NULL()
    TYPE(MyInst), POINTER       :: Inst     => NULL()

    !=================================================================
    ! InstRemove begins here!
    !=================================================================

    ! Init
    PrevInst => NULL()
    Inst     => NULL()

    ! Get instance. Also archive previous instance.
    CALL InstGet ( Instance, Inst, RC, PrevInst=PrevInst )

    ! Instance-specific deallocation
    IF ( ASSOCIATED(Inst) ) THEN

       !---------------------------------------------------------------------
       ! Deallocate fields of Inst before popping off from the list
       ! in order to avoid memory leaks (Bob Yantosca (17 Aug 2022)
       !---------------------------------------------------------------------   
       IF ( ASSOCIATED( Inst%PROFILE ) ) THEN
          DEALLOCATE( Inst%PROFILE )
       ENDIF
       Inst%PROFILE => NULL()

       IF ( ASSOCIATED( Inst%SLBASE ) ) THEN
          DEALLOCATE( Inst%SLBASE )
       ENDIF
       Inst%SLBASE => NULL()

       IF ( ASSOCIATED( Inst%FLASH_DENS_TOT ) ) THEN
          DEALLOCATE( Inst%FLASH_DENS_TOT )
       ENDIF
       Inst%FLASH_DENS_TOT => NULL()

       !IF ( ASSOCIATED( Inst%FLASH_DENS_IC ) ) THEN
       !   DEALLOCATE( Inst%FLASH_DENS_IC )
       !ENDIF
       Inst%FLASH_DENS_IC => NULL()

       !IF ( ASSOCIATED( Inst%FLASH_DENS_CG ) ) THEN
       !   DEALLOCATE ( Inst%FLASH_DENS_CG )
       !ENDIF
       Inst%FLASH_DENS_CG => NULL()

       IF ( ASSOCIATED( Inst%CONV_DEPTH ) ) THEN
          DEALLOCATE( Inst%CONV_DEPTH )
       ENDIF
       Inst%CONV_DEPTH => NULL()

       IF ( ALLOCATED ( Inst%SpcScalVal ) ) THEN
          DEALLOCATE( Inst%SpcScalVal )
       ENDIF

       IF ( ALLOCATED( Inst%SpcScalFldNme ) ) THEN
          DEALLOCATE( Inst%SpcScalFldNme )
       ENDIF

       !---------------------------------------------------------------------
       ! Pop off instance from list
       !---------------------------------------------------------------------
       IF ( ASSOCIATED(PrevInst) ) THEN
          PrevInst%NextInst => Inst%NextInst
       ELSE
          AllInst => Inst%NextInst
       ENDIF
       DEALLOCATE(Inst)
    ENDIF

    ! Free pointers before exiting
    PrevInst => NULL()
    Inst     => NULL()

   END SUBROUTINE InstRemove
!EOC
END MODULE HCOX_LightNOx_Mod
