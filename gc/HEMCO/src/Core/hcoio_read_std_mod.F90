!BOC
#if defined ( MODEL_GCCLASSIC ) || defined( MODEL_WRF ) || defined( MODEL_CESM ) || defined( HEMCO_STANDALONE )
! The 'standard' HEMCO I/O module is used for:
! - HEMCO Standalone (HEMCO_STANDALONE)
! - GEOS-Chem 'Classic' (MODEL_GCCLASSIC)
! - WRF-GC (MODEL_WRF)
! - CESM-GC and CAM-Chem / HEMCO-CESM (MODEL_CESM)
!EOC
!------------------------------------------------------------------------------
!                   Harmonized Emissions Component (HEMCO)                    !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: hcoio_read_std_mod.F90
!
! !DESCRIPTION: Module HCOIO\_read\_mod controls data processing
! (file reading, unit conversion, regridding) for HEMCO in the
! 'standard' environment (i.e. non-ESMF).
!
! This module implements the 'standard' environment (i.e. non-ESMF).
!\\
!\\
! !INTERFACE:
!
MODULE HCOIO_Read_Mod
!
! !USES:
!
  USE HCO_Types_Mod
  USE HCO_Error_Mod
  USE HCO_CharTools_Mod
  USE HCO_State_Mod,       ONLY : Hco_State
  USE HCOIO_Util_Mod

  IMPLICIT NONE
  PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!
  PUBLIC  :: HCOIO_Read
  PUBLIC  :: HCOIO_CloseAll
!
! !REMARKS:
!  Beginning with HEMCO 3.0.0, all I/O modules use the same module names,
!  and their compilation depends on pre-processor flags defined at the top
!  of the file.
!
!  This is to streamline the implementation of one unified Data Input Layer,
!  that can be switched in and out at compile time, and reduce branching of
!  code paths elsewhere.
!
! !REVISION HISTORY:
!  22 Aug 2013 - C. Keller   - Initial version
!  See https://github.com/geoschem/hemco for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !DEFINED PARAMETERS
!
  ! Parameter used for difference testing of floating points
  REAL(dp), PRIVATE, PARAMETER :: EPSILON = 1.0e-5_dp

#if defined( MODEL_CESM ) || defined( MODEL_WRF )
  REAL(hp), PRIVATE            :: GC_72_EDGE_SIGMA(73) = (/1.000000E+00, 9.849998E-01, 9.699136E-01, 9.548285E-01, 9.397434E-01, 9.246593E-01, 9.095741E-01, 8.944900E-01, 8.794069E-01, 8.643237E-01, 8.492406E-01, 8.341584E-01, 8.190762E-01, 7.989697E-01, 7.738347E-01, 7.487007E-01, 7.235727E-01, 6.984446E-01, 6.733175E-01, 6.356319E-01, 5.979571E-01, 5.602823E-01, 5.226252E-01, 4.849751E-01, 4.473417E-01, 4.097261E-01, 3.721392E-01, 3.345719E-01, 2.851488E-01, 2.420390E-01, 2.055208E-01, 1.746163E-01, 1.484264E-01, 1.261653E-01, 1.072420E-01, 9.115815E-02, 7.748532E-02, 6.573205E-02, 5.565063E-02, 4.702097E-02, 3.964964E-02, 3.336788E-02, 2.799704E-02, 2.341969E-02, 1.953319E-02, 1.624180E-02, 1.346459E-02, 1.112953E-02, 9.171478E-03, 7.520355E-03, 6.135702E-03, 4.981002E-03, 4.023686E-03, 3.233161E-03, 2.585739E-03, 2.057735E-03, 1.629410E-03, 1.283987E-03, 1.005675E-03, 7.846040E-04, 6.089317E-04, 4.697755E-04, 3.602270E-04, 2.753516E-04, 2.082408E-04, 1.569208E-04, 1.184308E-04, 8.783617E-05, 6.513694E-05, 4.737232E-05, 3.256847E-05, 1.973847E-05, 9.869233E-06/)
#endif

CONTAINS
!EOC
!------------------------------------------------------------------------------
!                   Harmonized Emissions Component (HEMCO)                    !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: HCOIO_Read
!
! !DESCRIPTION: Reads a netCDF file and returns the regridded array in proper
! units. This routine uses the HEMCO generic data reading and regridding
! routines.
!\\
!\\
! Two different regridding algorithm are used: NCREGRID for 3D data with
! vertical regridding, and map\_a2a for all other data. map\_a2a also
! supports index-based remapping, while this feature is currently not
! possible in combination with NCREGRID.
!\\
!\\
! 3D data is vertically regridded onto the simulation grid on the sigma
! interface levels. In order to calculate these levels correctly, the netCDF
! vertical coordinate description must adhere to the CF - conventions. See
! routine NC\_Get\_Sigma\_Levels in Ncdf\_Mod for more details.
!\\
!\\
! A simpler vertical interpolation scheme is used if (a) the number of
! vertical levels of the input data corresponds to the number of levels
! on the simulation grid (direct mapping, no remapping), (b) the vertical
! level variable name (long\_name) contains the word "GEOS-Chem level". In
! the latter case, the vertical levels of the input data is interpreted as
! GEOS vertical levels and mapped onto the simulation grid using routine
! ModelLev\_Interpolate.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE HCOIO_Read( HcoState, Lct, RC )
!
! !USES:
!
    USE HCO_Ncdf_Mod,       ONLY : NC_Open
    USE HCO_Ncdf_Mod,       ONLY : NC_Close
    USE HCO_Ncdf_Mod,       ONLY : NC_Read_Var
    USE HCO_Ncdf_Mod,       ONLY : NC_Read_Arr
    USE HCO_Ncdf_Mod,       ONLY : NC_Get_Grid_Edges
    USE HCO_Ncdf_Mod,       ONLY : NC_Get_Sigma_Levels
    USE HCO_Ncdf_Mod,       ONLY : NC_IsModelLevel
    USE HCO_Ncdf_Mod,       ONLY : NC_IsSigmaLevel
    USE HCO_CHARPAK_MOD,    ONLY : TRANLC
    USE HCO_Unit_Mod,       ONLY : HCO_Unit_Change
    USE HCO_Unit_Mod,       ONLY : HCO_Unit_ScalCheck
    USE HCO_Unit_Mod,       ONLY : HCO_IsUnitless
    USE HCO_Unit_Mod,       ONLY : HCO_IsIndexData
    USE HCO_Unit_Mod,       ONLY : HCO_UnitTolerance
    USE HCO_GeoTools_Mod,   ONLY : HCO_ValidateLon
    USE HCO_FileData_Mod,   ONLY : FileData_ArrCheck
    USE HCO_FileData_Mod,   ONLY : FileData_ArrInit
    USE HCO_FileData_Mod,   ONLY : FileData_Cleanup
    USE HCOIO_MESSY_MOD,    ONLY : HCO_MESSY_REGRID
    USE HCO_INTERP_MOD,     ONLY : REGRID_MAPA2A
    USE HCO_INTERP_MOD,     ONLY : ModelLev_Check
    USE HCO_CLOCK_MOD,      ONLY : HcoClock_Get
    USE HCO_DIAGN_MOD,      ONLY : Diagn_Update
    USE HCO_EXTLIST_MOD,    ONLY : HCO_GetOpt
    USE HCO_TIDX_MOD,       ONLY : tIDx_IsInRange

    include "netcdf.inc"
!
! !INPUT PARAMETERS:
!
    TYPE(HCO_State),  POINTER        :: HcoState   ! HEMCO state object
    TYPE(ListCont),   POINTER        :: Lct        ! HEMCO list container
!
! !INPUT/OUTPUT PARAMETERS:
!
    INTEGER,          INTENT(INOUT)  :: RC         ! Success or failure?
!
! !REVISION HISTORY:
!  13 Mar 2013 - C. Keller   - Initial version
!  See https://github.com/geoschem/hemco for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    CHARACTER(LEN=255)            :: thisUnit, LevUnit, LevName
    CHARACTER(LEN=1023)           :: MSG, LOC
    CHARACTER(LEN=1023)           :: srcFile, srcFile2
    INTEGER                       :: NX, NY
    INTEGER                       :: NCRC, Flag, AS
    INTEGER                       :: ncLun, ncLun2
    INTEGER                       :: ierr,   v_id
    INTEGER                       :: nlon,   nlat,  nlev, nTime
    INTEGER                       :: lev1,   lev2,  dir
    INTEGER                       :: tidx1,  tidx2,  ncYr,  ncMt
    INTEGER                       :: tidx1b, tidx2b, ncYr2, ncMt2
    INTEGER                       :: HcoID
    INTEGER                       :: ArbIdx
    INTEGER                       :: nlatEdge, nlonEdge
    INTEGER                       :: Direction
    REAL(hp)                      :: MW_g
    REAL(sp)                      :: wgt1,   wgt2
    REAL(sp), POINTER             :: ncArr(:,:,:,:)
    REAL(sp), POINTER             :: ncArr2(:,:,:,:)
    REAL(hp), POINTER             :: SigEdge(:,:,:)
    REAL(hp), POINTER             :: SigLev (:,:,:)
    REAL(hp), POINTER             :: LonMid   (:)
    REAL(hp), POINTER             :: LatMid   (:)
    REAL(hp), POINTER             :: LevMid   (:)
    REAL(hp), POINTER             :: LonEdge  (:)
    REAL(hp), POINTER             :: LatEdge  (:)
    REAL(hp)                      :: UnitFactor
    LOGICAL                       :: KeepSpec
    LOGICAL                       :: FOUND
    LOGICAL                       :: IsModelLevel
    LOGICAL                       :: DoReturn
    INTEGER                       :: UnitTolerance
    INTEGER                       :: AreaFlag, TimeFlag
    REAL(dp)                      :: YMDhma,  YMDhmb, YMDhm1
    REAL(dp)                      :: oYMDhm1, oYMDhm2
    INTEGER                       :: cYr, cMt, cDy, cHr, Yr1, Yr2
    INTEGER                       :: nYears, iYear
    INTEGER                       :: I, J

    ! Use MESSy regridding routines?
    LOGICAL                       :: UseMESSy

    ! SAVEd scalars
    LOGICAL, SAVE                 :: doPrintWarning = .TRUE.

    !=================================================================
    ! HCOIO_READ begins here
    !=================================================================
    LOC = 'HCOIO_READ (HCOIO_READ_STD_MOD.F90)'

    ! Enter
    CALL HCO_ENTER( HcoState%Config%Err, LOC, RC )
    IF ( RC /= HCO_SUCCESS ) THEN
        CALL HCO_ERROR( 'ERROR 0', RC, THISLOC=LOC )
        RETURN
    ENDIF

    ! Initialize pointers
    ncArr   => NULL()
    ncArr2  => NULL()
    SigEdge => NULL()
    SigLev  => NULL()
    LonMid  => NULL()
    LatMid  => NULL()
    LevMid  => NULL()
    LonEdge => NULL()
    LatEdge => NULL()

    ! Zero local variables for safety's sake
    dir     =  0
    lev1    =  0
    lev2    =  0
    ncYr    =  0
    ncMt    =  0
    ncYr2   =  0
    ncMt2   =  0
    nLon    =  0
    nLat    =  0
    nLev    =  0
    nTime   =  0
    tIdx1   =  0
    tIdx2   =  0
    tidx1b  =  0
    tidx2b  =  0
    wgt1    =  0.0_sp
    wgt2    =  0.0_sp

    ! Get unit tolerance set in configuration file
    UnitTolerance = HCO_UnitTolerance( HcoState%Config )

    ! For convenience, copy horizontal grid dimensions from HEMCO
    ! state object
    NX = HcoState%NX
    NY = HcoState%NY

    ! Verbose
    IF ( HCO_IsVerb(HcoState%Config%Err,2) ) THEN
       WRITE(MSG,*) 'Processing container: ', TRIM(Lct%Dct%cName)
       CALL HCO_MSG( HcoState%Config%Err, MSG, SEP1='-' )
    ENDIF

    ! If the file has cycle flag "E" (e.g. it's a restart file), then we will
    ! read it only once and then never again.  If the file has already been
    ! read on a previous call, then don't call HCOIO_READ. (bmy, 10/4/18)
    !
    ! Moved this handling from hcoio_dataread_mod, as it is non-MAPL specific
    ! (hplin, 4/5/21)
    IF ( Lct%Dct%Dta%CycleFlag == HCO_CFLAG_EXACT .and.                      &
         Lct%Dct%Dta%UpdtFlag  == HCO_UFLAG_ONCE  .and.                      &
         Lct%Dct%Dta%isTouched                          ) THEN

       ! Print a warning message only once
       IF ( doPrintWarning ) THEN
          doPrintWarning = .FALSE.
          MSG = 'No further attempts will be made to read file: ' //         &
                TRIM( Lct%Dct%Dta%NcFile )
          CALL HCO_WARNING ( HcoState%Config%Err, MSG, RC, WARNLEV=1 )
       ENDIF

       ! Return without reading
       CALL HCO_LEAVE( HcoState%Config%Err, RC )
       RETURN
    ENDIF

    ! ----------------------------------------------------------------
    ! Parse source file name. This will replace all tokens ($ROOT,
    ! ($YYYY), etc., with valid values.
    ! ----------------------------------------------------------------
    CALL SrcFile_Parse( HcoState, Lct, srcFile, FOUND, RC )
    IF ( RC /= HCO_SUCCESS ) THEN
       MSG = 'Error encountered in routine "SrcFile_Parse", located '     // &
             'module src/Core/hcoio_read_std_mod.F90!'
        CALL HCO_ERROR( MSG, RC )
        RETURN
    ENDIF

    ! Handle found or not in the standard way if HEMCO is in regular run mode.
    IF ( .NOT. HcoState%Options%isDryRun ) THEN

       !====================================================================
       ! HEMCO is in regular simulation mode (not dry-run)!
       !====================================================================

       ! If file not found, return w/ error. No error if cycling attribute is
       ! select to range. In that case, just make sure that array is empty.
       IF ( .NOT. FOUND ) THEN
          IF ( ( Lct%Dct%Dta%CycleFlag == HCO_CFLAG_RANGE ) .OR.      &
               ( Lct%Dct%Dta%CycleFlag == HCO_CFLAG_EXACT )     ) THEN

             ! If MustFind flag is enabled, return with error if field is not
             ! found
             IF ( Lct%Dct%Dta%MustFind ) THEN
                MSG = 'Cannot find file for current simulation time: ' // &
                     TRIM(srcFile) // ' - Cannot get field ' // &
                     TRIM(Lct%Dct%cName) // '. Please check file name ' // &
                     'and time (incl. time range flag) in the config. file'
                CALL HCO_ERROR( MSG, RC )
                RETURN

             ! If MustFind flag is not enabled, ignore this field and return
             ! with a warning.
             ELSE
                CALL FileData_Cleanup( Lct%Dct%Dta, DeepClean=.FALSE. )
                MSG = 'No valid file found for current simulation time - data '// &
                     'will be ignored for time being - ' // TRIM(Lct%Dct%cName)
                CALL HCO_WARNING ( HcoState%Config%Err, MSG, RC, WARNLEV=3 )
                CALL HCO_LEAVE ( HcoState%Config%Err,  RC )
                RETURN
             ENDIF

          ELSE
             MSG = 'Cannot find file for current simulation time: ' // &
                  TRIM(srcFile) // ' - Cannot get field ' // &
                  TRIM(Lct%Dct%cName) // '. Please check file name ' // &
                  'and time (incl. time range flag) in the config. file'
             CALL HCO_ERROR( MSG, RC )
             RETURN
          ENDIF
       ENDIF

    ELSE

       !====================================================================
       ! HEMCO is in a "dry-run" mode!
       !====================================================================

       ! Simulate file read buffer
       IF ( TRIM(HcoState%ReadLists%FileInArchive) == TRIM(srcFile) ) THEN
          CALL HCO_LEAVE ( HcoState%Config%Err,  RC )
          RETURN
       ENDIF

       ! If file exists, print the result. If NOT, then handle accordingly
       ! But NEVER error out (HCO_ERROR), as we want to get a list of all
       ! files. (hplin, 11/2/19)
       IF ( .NOT. FOUND ) THEN
          IF ( ( Lct%Dct%Dta%CycleFlag == HCO_CFLAG_RANGE )       .OR.       &
               ( Lct%Dct%Dta%CycleFlag == HCO_CFLAG_EXACT )     ) THEN

             ! If MustFind flag is enabled, return with error if field is not
             ! found
             IF ( Lct%Dct%Dta%MustFind ) THEN
                MSG = 'Cannot find file for current simulation time: '    // &
                     TRIM(srcFile) // ' - Cannot get field '              // &
                     TRIM(Lct%Dct%cName) // '. Please check file name '   // &
                     'and time (incl. time range flag) in the config. file'
                CALL HCO_Warning( HcoState%Config%Err, MSG, RC, WARNLEV=3 )

                ! Write a msg to stdout (NOT FOUND)
                WRITE( 6, 300 ) TRIM( srcFile )
 300            FORMAT( 'HEMCO: REQUIRED FILE NOT FOUND ', a )

             ! If MustFind flag is not enabled, ignore this field and return
             ! with a warning.
             ELSE
                CALL FileData_Cleanup( Lct%Dct%Dta, DeepClean=.FALSE. )
                MSG = 'No valid file found for current simulation time - '// &
                     'data will be ignored for time being - '             // &
                     TRIM(Lct%Dct%cName)
                CALL HCO_WARNING ( HcoState%Config%Err, MSG, RC, WARNLEV=3 )

                ! Write a msg to stdout (OPTIONAL)
                WRITE( 6, 310 ) TRIM( srcFile )
 310            FORMAT( 'HEMCO: OPTIONAL FILE NOT FOUND ', a )

             ENDIF

          ! Not range or exact
          ELSE
             MSG = 'Cannot find file for current simulation time: '       // &
                  TRIM(srcFile) // ' - Cannot get field '                 // &
                  TRIM(Lct%Dct%cName) // '. Please check file name '      // &
                  'and time (incl. time range flag) in the config. file'
             CALL HCO_WARNING ( HcoState%Config%Err, MSG, RC, WARNLEV=3 )

             ! Write a msg to stdout (NOT FOUND)
             WRITE( 6, 300 ) TRIM(srcFile)

          ENDIF
       ELSE

          ! Write a mesage to stdout (HEMCO: Opening...)
          WRITE( 6, 100 ) TRIM( srcFile )

       ENDIF

       ! It is safe to leave now, we do not need to handle opening the file.
       ! This may be changed in the future if the "dry-run" mode requires
       ! a check of the file contents... this may be beyond the scope for now.

       ! Simulate the "reading" in netCDF to prevent duplicate entries
       ! in the log
       HcoState%ReadLists%FileInArchive = TRIM(srcFile)

       ! Skip further processing
       CALL HCO_LEAVE ( HcoState%Config%Err,  RC )
       RETURN
    ENDIF ! End of dry-run mode else clause

    ! ----------------------------------------------------------------
    ! Open netCDF
    ! ----------------------------------------------------------------

    ! Check if file is already in buffer. In that case use existing
    ! open stream. Otherwise open new file. At any given time there
    ! can only be one file in buffer.
    ncLun = -1
    IF ( HcoState%ReadLists%FileLun > 0 ) THEN
       IF ( TRIM(HcoState%ReadLists%FileInArchive) == TRIM(srcFile) ) THEN
          ncLun = HcoState%ReadLists%FileLun
       ELSE
          CALL NC_CLOSE ( HcoState%ReadLists%FileLun )
          HcoState%ReadLists%FileLun = -1
       ENDIF
    ENDIF

    ! To read from existing stream:
    IF ( ncLun > 0 ) THEN

       ! Verbose mode
       IF ( HCO_IsVerb(HcoState%Config%Err,2) ) THEN
          WRITE(MSG,*) 'Reading from existing stream: ', TRIM(srcFile)
          CALL HCO_MSG( HcoState%Config%Err, MSG )
       ENDIF

    ! To open a new file:
    ELSE
       CALL NC_OPEN ( TRIM(srcFile), ncLun )

       ! Verbose mode
       IF ( HCO_IsVerb(HcoState%Config%Err,1) ) THEN
          WRITE(MSG,*) 'Opening file: ', TRIM(srcFile)
          CALL HCO_MSG( HcoState%Config%Err, MSG )
       ENDIF

       ! Also write to standard output
       WRITE( 6, 100 ) TRIM( srcFile )
 100   FORMAT( 'HEMCO: Opening ', a )

       ! This is now the file in archive
       HcoState%ReadLists%FileInArchive = TRIM(srcFile)
       HcoState%ReadLists%FileLun       = ncLun
    ENDIF

    ! ----------------------------------------------------------------
    ! Extract time slice information
    ! This determines the lower and upper time slice index (tidx1
    ! and tidx2) to be read based upon the time slice information
    ! extracted from the file and the time stamp settings set in the
    ! HEMCO configuration file. Multiple time slices are only selected
    ! for weekdaily data or for 'autodetected' hourly data (using the
    ! wildcard character in the configuration file time attribute) or
    ! if data shall be interpolated between two (consecutive) time
    ! slices. The weights to be assigned to those two time slices is
    ! also calculated in GET_TIMEIDX and returned as variables wgt1
    ! and wgt2, respectively.
    ! ----------------------------------------------------------------
    CALL GET_TIMEIDX ( HcoState,  Lct,                &
                       ncLun,     tidx1,    tidx2,    &
                       wgt1,      wgt2,     oYMDhm1,  &
                       YMDhma,    YMDhm1,   RC        )
    IF ( RC /= HCO_SUCCESS ) THEN
        CALL HCO_ERROR( 'ERROR 1', RC, THISLOC=LOC )
        RETURN
    ENDIF

    !-----------------------------------------------------------------
    ! Check for negative tidx1. tidx1 can still be negative if:
    ! (a) CycleFlag is set to range and the current simulation
    ! time is outside of the data time range. In this case, we
    ! prompt a warning and make sure that there is no data
    ! associated with this FileData container.
    ! (b) CycleFlag is set to exact and none of the data time
    ! stamps matches the current simulation time exactly. Return
    ! with error!
    !-----------------------------------------------------------------
    IF ( tidx1 < 0 ) THEN
       DoReturn = .FALSE.
       IF ( Lct%Dct%Dta%CycleFlag == HCO_CFLAG_CYCLE ) THEN
          MSG = 'Invalid time index in ' // TRIM(srcFile)
          CALL HCO_ERROR( MSG, RC )
          DoReturn = .TRUE.
       ELSEIF ( ( Lct%Dct%Dta%CycleFlag == HCO_CFLAG_RANGE ) .OR.      &
                ( Lct%Dct%Dta%CycleFlag == HCO_CFLAG_EXACT )     ) THEN
          IF ( Lct%Dct%Dta%MustFind ) THEN
             MSG = 'Cannot find field with valid time stamp in ' // &
                   TRIM(srcFile) // ' - Cannot get field ' // &
                   TRIM(Lct%Dct%cName) // '. Please check file name ' // &
                   'and time (incl. time range flag) in the config. file'
             CALL HCO_ERROR( MSG, RC )
             DoReturn = .TRUE.
          ELSE
             CALL FileData_Cleanup( Lct%Dct%Dta, DeepClean=.FALSE.)
             MSG = 'Simulation time is outside of time range provided for '//&
                  TRIM(Lct%Dct%cName) // ' - field is ignored for the time being!'
             CALL HCO_WARNING ( HcoState%Config%Err, MSG, RC, WARNLEV=3 )
             DoReturn = .TRUE.
             CALL HCO_LEAVE ( HcoState%Config%Err,  RC )
          ENDIF
       ENDIF

       ! Eventually return here
       IF ( DoReturn ) THEN
          RETURN
       ENDIF
    ENDIF

    ! ----------------------------------------------------------------
    ! Check if variable is in file
    ! ----------------------------------------------------------------
    ierr = Nf_Inq_Varid( ncLun, Lct%Dct%Dta%ncPara, v_id )
    IF ( ierr /= NF_NOERR ) THEN

       ! If MustFind flag is enabled, return with error if field is not
       ! found
       IF ( Lct%Dct%Dta%MustFind ) THEN
          MSG = 'Cannot find field ' // TRIM(Lct%Dct%cName) // &
                '. Please check variable name in the config. file'
          CALL HCO_ERROR( MSG, RC )
          RETURN

       ! If MustFind flag is not enabled, ignore this field and return
       ! with a warning.
       ELSE
          CALL FileData_Cleanup( Lct%Dct%Dta, DeepClean=.FALSE. )
          MSG = 'Cannot find field ' // TRIM(Lct%Dct%cName) // &
                '. Will be ignored for time being.'
          CALL HCO_WARNING ( HcoState%Config%Err, MSG, RC, WARNLEV=3 )
          CALL HCO_LEAVE ( HcoState%Config%Err,  RC )
          RETURN
       ENDIF
    ENDIF

    ! ----------------------------------------------------------------
    ! Read grid
    ! ----------------------------------------------------------------

    ! Extract longitude midpoints
    CALL NC_READ_VAR ( ncLun, 'lon', nlon, thisUnit, LonMid, NCRC )
    IF ( NCRC /= 0 ) THEN
       CALL HCO_ERROR( 'NC_READ_VAR: lon', RC )
       RETURN
    ENDIF

    IF ( nlon == 0 ) THEN
       CALL NC_READ_VAR ( ncLun, 'longitude', nlon, thisUnit, LonMid, NCRC )
    ENDIF
    IF ( NCRC /= 0 ) THEN
       CALL HCO_ERROR( 'NC_READ_VAR: longitude', RC )
       RETURN
    ENDIF

    IF ( nlon == 0 ) THEN
       CALL NC_READ_VAR ( ncLun, 'Longitude', nlon, thisUnit, LonMid, NCRC )
    ENDIF
    IF ( NCRC /= 0 ) THEN
       CALL HCO_ERROR( 'NC_READ_LON: Longitude', RC )
       RETURN
    ENDIF

    IF ( nlon == 0 ) THEN
       MSG = 'Cannot find longitude variable in ' // TRIM(srcFile) // &
             ' - Must be one of `lon`, `longitude`, `Longitude`'
       CALL HCO_ERROR( MSG, RC )
       RETURN
    ENDIF

    ! Unit must be degrees_east
    CALL TRANLC( thisUnit)
    IF ( INDEX( thisUnit, 'degrees_east' ) == 0 ) THEN
       MSG = 'illegal longitude unit in ' // TRIM(srcFile) // &
             ' - Must be `degrees_east`.'
       CALL HCO_ERROR( MSG, RC )
       RETURN
    ENDIF

    ! Make sure longitude is steadily increasing.
    CALL HCO_ValidateLon( HcoState, nlon, LonMid, RC )
    IF ( RC /= HCO_SUCCESS ) THEN
        CALL HCO_ERROR( 'ERROR 2', RC, THISLOC=LOC )
        RETURN
    ENDIF

    ! Extract latitude midpoints
    CALL NC_READ_VAR ( ncLun, 'lat', nlat, thisUnit, LatMid, NCRC )
    IF ( NCRC /= 0 ) THEN
       CALL HCO_ERROR( 'NC_READ_LON: lat', RC )
       RETURN
    ENDIF

    IF ( nlat == 0 ) THEN
       CALL NC_READ_VAR ( ncLun, 'latitude', nlat, thisUnit, LatMid, NCRC )
    ENDIF
    IF ( NCRC /= 0 ) THEN
       CALL HCO_ERROR( 'NC_READ_LON: latitude', RC )
       RETURN
    ENDIF

    IF ( nlat == 0 ) THEN
       CALL NC_READ_VAR ( ncLun, 'Latitude', nlat, thisUnit, LatMid, NCRC )
    ENDIF
    IF ( NCRC /= 0 ) THEN
       CALL HCO_ERROR( 'NC_READ_LON: Latitude', RC )
       RETURN
    ENDIF

    IF ( nlat == 0 ) THEN
       MSG = 'Cannot find latitude variable in ' // TRIM(srcFile) // &
             ' - Must be one of `lat`, `latitude`, `Latitude`'
       CALL HCO_ERROR( MSG, RC )
       RETURN
    ENDIF

    ! Unit must be degrees_north
    CALL TRANLC( thisUnit)
    IF ( INDEX( thisUnit, 'degrees_north' ) == 0 ) THEN
       MSG = 'illegal latitude unit in ' // TRIM(srcFile) // &
             ' - Must be `degrees_north`.'
       CALL HCO_ERROR( MSG, RC )
       RETURN
    ENDIF

    ! Get level index if we are dealing with 3D data
    IF ( Lct%Dct%Dta%SpaceDim == 3 ) THEN

       ! Try to extract level midpoints
       LevName = 'lev'
       CALL NC_READ_VAR ( ncLun, LevName, nlev, LevUnit, LevMid, NCRC )
       IF ( NCRC /= 0 ) THEN
          CALL HCO_ERROR( 'NC_READ_VAR: lev', RC )
          RETURN
       ENDIF
       IF ( nlev == 0 ) THEN
          LevName = 'height'
          CALL NC_READ_VAR ( ncLun, LevName, nlev, LevUnit, LevMid, NCRC )
          IF ( NCRC /= 0 ) THEN
             CALL HCO_ERROR( 'NC_READ_VAR: height', RC )
             RETURN
          ENDIF
       ENDIF
       IF ( nlev == 0 ) THEN
          LevName = 'level'
          CALL NC_READ_VAR ( ncLun, LevName, nlev, LevUnit, LevMid, NCRC )
          IF ( NCRC /= 0 ) THEN
             CALL HCO_ERROR( 'NC_READ_VAR: level', RC )
             RETURN
          ENDIF
       ENDIF

       ! Error check
       IF ( nlev == 0 ) THEN
          MSG = 'Cannot find vertical coordinate variable in ' // &
                 TRIM(SrcFile) // ' - Must be one of `lev`, `level`, `height`.'
          CALL HCO_ERROR( MSG, RC )
          RETURN
       ENDIF

       ! Are these model levels? This will only return true if the long
       ! name of the level variable contains "GEOS-Chem level".
       ! For now, we assume levels are already on model levels if the
       ! number of levels to be read is explicitly set in the configuration
       ! file (ckeller, 5/20/15).
       IF ( Lct%Dct%Dta%Levels == 0 ) THEN

          ! Check if vertical coordinate is GEOS-Chem levels
          IsModelLevel = NC_IsModelLevel( ncLun, LevName )

          ! Further check if the given number of vertical levels should be
          ! treated as model levels. This is the case if e.g. the nuber of
          ! levels found on the file exactly matches the number of vertical
          ! levels of the grid. Some of these assumptions are rather arbitrary.
          ! IsModelLev will stay True if is was set so in NC_ISMODELLEVEL
          ! above. (ckeller, 9/29/15)
          CALL ModelLev_Check( HcoState, nlev, IsModelLevel, RC )
          IF ( RC /= HCO_SUCCESS ) THEN
              CALL HCO_ERROR( 'ERROR 3', RC, THISLOC=LOC )
              RETURN
          ENDIF

          ! Override IsModelLevel if the long_name contains
          ! "atmospheric_hybrid_sigma_pressure_coordinate"
          IsModelLevel = ( .not. NC_IsSigmaLevel( ncLun, LevName ) )

          ! Set level indeces to be read
          lev1 = 1
          lev2 = nlev

       ! If levels are explicitly given:
       ELSE

          ! If long_name is "atmospheric_hybrid_sigma_pressure_coordinate",
          ! then treat it as sigma levels; otherwise assume model levels.
          IsModelLevel = ( .not. NC_IsSigmaLevel( ncLun, LevName ) )

          ! Number of levels to be read must be smaller or equal to total
          ! number of available levels
          IF ( ABS(Lct%Dct%Dta%Levels) > nlev ) THEN
             WRITE(MSG,*) Lct%Dct%Dta%Levels, ' levels requested but file ', &
                'has only ', nlev, ' levels: ', TRIM(Lct%Dct%cName)
             CALL HCO_ERROR( MSG, RC )
             RETURN
          ENDIF

          ! Set levels to be read
          IF ( Lct%Dct%Dta%Levels > 0 ) THEN
             lev1 = 1
             lev2 = Lct%Dct%Dta%Levels

          ! Reverse axis!
          ELSE
             lev1 = nlev
             lev2 = nlev + Lct%Dct%Dta%Levels + 1
          ENDIF

       ENDIF

       ! Verbose
       IF ( HCO_IsVerb(HcoState%Config%Err,2) ) THEN
          WRITE(MSG,*) 'Will read vertical levels ', lev1, ' to ', lev2
          CALL HCO_MSG(HcoState%Config%Err,MSG)
       ENDIF

    ! For 2D data, set lev1 and lev2 to zero. This will ignore
    ! the level dimension in the netCDF reading call that follows.
    ELSE
       nlev        = 0
       lev1        = 0
       lev2        = 0
       IsModelLevel = .FALSE.
    ENDIF

    ! ----------------------------------------------------------------
    ! Check for arbitrary additional dimension. Will return -1 if not
    ! set.
    ! ----------------------------------------------------------------
    CALL GetArbDimIndex( HcoState, ncLun, Lct, ArbIdx, RC )
    IF ( RC /= HCO_SUCCESS ) THEN
        CALL HCO_ERROR( 'ERROR 4', RC, THISLOC=LOC )
        RETURN
    ENDIF

    ! ----------------------------------------------------------------
    ! Read data
    ! ----------------------------------------------------------------

    ! Verbose mode
    IF ( HCO_IsVerb(HcoState%Config%Err,2) ) THEN
       WRITE(MSG,*) 'Reading variable ', TRIM(Lct%Dct%Dta%ncPara)
       CALL HCO_MSG(HcoState%Config%Err,MSG)
    ENDIF

    CALL NC_READ_ARR( fID     = ncLun,              &
                      ncVar   = Lct%Dct%Dta%ncPara, &
                      lon1    = 1,                  &
                      lon2    = nlon,               &
                      lat1    = 1,                  &
                      lat2    = nlat,               &
                      lev1    = lev1,               &
                      lev2    = lev2,               &
                      time1   = tidx1,              &
                      time2   = tidx2,              &
                      ncArr   = ncArr,              &
                      varUnit = thisUnit,           &
                      wgt1    = wgt1,               &
                      wgt2    = wgt2,               &
                      MissVal = HCO_MISSVAL,        &
                      ArbIdx  = ArbIdx,             &
                      RC      = NCRC                 )

    IF ( NCRC /= 0 ) THEN
       CALL HCO_ERROR( 'NC_READ_ARRAY', RC )
       RETURN
    ENDIF

    ! Check for missing values: set base emissions and masks to 0, and
    ! scale factors to 1. This will make sure that these entries will
    ! be ignored.
!!! CALL CheckMissVal ( Lct, ncArr )

    !-----------------------------------------------------------------
    ! Eventually do interpolation between files. This is a pretty
    ! crude implementation for data interpolation between different
    ! files. It is only applied to data that is marked as interpolated
    ! data and if no appropriate interpolation date could be found in
    ! the first file. This will only be the case if the preferred date-
    ! time is outside the file range.
    !-----------------------------------------------------------------
    IF ( Lct%Dct%Dta%CycleFlag == HCO_CFLAG_INTER .AND. wgt1 < 0.0_sp ) THEN

       ! No need to read another file if the previous file had exactly the
       ! time stamp we were looking for.
       IF ( oYMDhm1 == YMDhma ) THEN
          FOUND = .FALSE.
       ELSE
          ! Check if there exists another file for a future/previous date
          ! oYMDhm1 is the originally preferred date, YMDhma is the date read
          ! from the file. If YMDhma is in the future, we need to look for a
          ! file in the past. Else, need to look for a file in the future.
          IF ( oYMDhm1 < YMDhma ) THEN
             Direction = -1
          ELSE
             Direction = +1
          ENDIF
          CALL SrcFile_Parse ( HcoState,  Lct, srcFile2, &
                               FOUND, RC, Direction = Direction )
          IF ( RC /= HCO_SUCCESS ) THEN
              CALL HCO_ERROR( 'ERROR 5', RC, THISLOC=LOC )
              RETURN
          ENDIF
       ENDIF

       ! If found, read data. Assume that all meta-data is the same.
       IF ( FOUND ) THEN

          ! Open file
          CALL NC_OPEN ( TRIM(srcFile2), ncLun2 )

          ! Define time stamp to be read. Use this call only
          ! to get the datetime of the first time slice (YMDhm1).
          ! All other values will be ignored and reset below.
          CALL GET_TIMEIDX ( HcoState,  Lct,               &
                             ncLun2,    tidx1,    tidx2,   &
                             wgt1,      wgt2,     oYMDhm2, &
                             YMDhmb,    YMDhm1,   RC       )
          IF ( RC /= HCO_SUCCESS ) THEN
              CALL HCO_ERROR( 'ERROR 6', RC, THISLOC=LOC )
              RETURN
          ENDIF

          ! Always read first time slice
          tidx1 = 1
          tidx2 = 1
          wgt1  = -1.0_sp
          wgt2  = -1.0_sp

          ! Read data and write into array ncArr2
          CALL NC_READ_ARR( fID     = ncLun2,             &
                            ncVar   = Lct%Dct%Dta%ncPara, &
                            lon1    = 1,                  &
                            lon2    = nlon,               &
                            lat1    = 1,                  &
                            lat2    = nlat,               &
                            lev1    = lev1,               &
                            lev2    = lev2,               &
                            time1   = tidx1,              &
                            time2   = tidx2,              &
                            ncArr   = ncArr2,             &
                            varUnit = thisUnit,           &
                            wgt1    = wgt1,               &
                            wgt2    = wgt2,               &
                            MissVal = HCO_MISSVAL,        &
                            ArbIdx  = ArbIdx,             &
                            RC      = NCRC                 )
          IF ( NCRC /= 0 ) THEN
             CALL HCO_ERROR( 'NC_READ_ARRAY (2)', RC )
             RETURN
          ENDIF

          ! Eventually fissing values
!!!          CALL CheckMissVal ( Lct, ncArr2 )

          ! Calculate weights to be applied to ncArr2 and ncArr1. These
          ! weights are calculated based on the originally preferred
          ! datetime oYMDh1 and the selected datetime of file 1 (YMDhma)
          ! and file 2 (YMDhm1)
          ! If date on file 1 < date on file 2:
          IF ( YMDhma < YMDhm1 ) THEN
             CALL GetWeights ( YMDhma, YMDhm1, oYMDhm1, wgt1, wgt2 )
          ! If date on file 1 > date on file 2:
          ELSEIF ( YMDhma > YMDhm1 ) THEN
             CALL GetWeights ( YMDhm1, YMDhma, oYMDhm1, wgt2, wgt1 )
          ! If both datetimes are for some reason the same (this should
          ! not happen!)
          ELSE
             wgt1 = 0.5_sp
             wgt2 = 0.5_sp
          ENDIF

          ! Apply weights
          ncArr = (wgt1 * ncArr) + (wgt2 * ncArr2)

          ! Verbose
          IF ( HCO_IsVerb(HcoState%Config%Err,2) ) THEN
             MSG = 'Interpolated data between two files:'
             CALL HCO_MSG(HcoState%Config%Err,MSG)
             MSG = '- File 1: ' // TRIM(srcFile)
             CALL HCO_MSG(HcoState%Config%Err,MSG)
             WRITE(MSG,*) '   Time stamp used: ', YMDhma
             CALL HCO_MSG(HcoState%Config%Err,MSG)
             WRITE(MSG,*) '   Applied weight: ', wgt1
             CALL HCO_MSG(HcoState%Config%Err,MSG)
             MSG = '- File 2: ' // TRIM(srcFile2)
             CALL HCO_MSG(HcoState%Config%Err,MSG)
             WRITE(MSG,*) '   Time stamp used: ', YMDhm1
             CALL HCO_MSG(HcoState%Config%Err,MSG)
             WRITE(MSG,*) '   Applied weight: ', wgt2
             CALL HCO_MSG(HcoState%Config%Err,MSG)
          ENDIF

          ! Cleanup
          IF ( ASSOCIATED(ncArr2) ) DEALLOCATE(ncArr2)

          ! Close file
          CALL NC_CLOSE ( ncLun2 )
       ENDIF !FOUND

    !-----------------------------------------------------------------
    ! Eventually calculate averages. Currently, averages are only
    ! calculated on the year dimension, e.g. over years.
    !-----------------------------------------------------------------
    ELSEIF ( Lct%Dct%Dta%CycleFlag == HCO_CFLAG_AVERG    .OR. &
             Lct%Dct%Dta%CycleFlag == HCO_CFLAG_RANGEAVG       ) THEN

       ! cYr is the current simulation year
       CALL HcoClock_Get( HcoState%Clock, cYYYY=cYr, cMM=cMt, cDD=cDy, &
                          cH=cHr, RC=RC )
       IF ( RC /= HCO_SUCCESS ) THEN
           CALL HCO_ERROR( 'ERROR 7', RC, THISLOC=LOC )
           RETURN
       ENDIF

       ! Determine year range to be read:
       ! By default, we would like to average between the year range given
       ! in the time attribute
       Yr1 = Lct%Dct%Dta%ncYrs(1)
       Yr2 = Lct%Dct%Dta%ncYrs(2)

       ! If averaging shall only be performed if outside the given
       ! range, check if current simulation date is within the range
       ! provied in the configuration file. If so, set year range to
       ! be read to current year only.
       IF ( ( Lct%Dct%Dta%CycleFlag == HCO_CFLAG_RANGEAVG ) ) THEN
          IF ( tIDx_IsInRange(Lct,cYr,cMt,cDy,cHr) ) THEN
             Yr1 = cYr
             Yr2 = cYr
          ENDIF
       ENDIF

       ! Total number of years to be read
       nYears = Yr2 - Yr1 + 1

       ! Read and add annual data if there is more than one year to be
       ! used.
       IF ( nYears > 1 ) THEN

          ! Cleanup ncArr. This is refilled again
          ncArr = 0.0_sp

          DO iYear = Yr1, Yr2

             ! Get file name for this year
             CALL SrcFile_Parse ( HcoState, Lct, srcFile2, &
                                  FOUND, RC, Year=iYear )
             IF ( RC /= HCO_SUCCESS ) THEN
                 CALL HCO_ERROR( 'ERROR 8', RC, THISLOC=LOC )
                 RETURN
             ENDIF

             ! If found, read data. Assume that all meta-data is the same.
             IF ( .NOT. FOUND ) THEN
                WRITE(MSG,*) 'Cannot find file for year ', iYear, ' - needed ', &
                   'to perform time-averaging on file ', TRIM(Lct%Dct%Dta%ncFile)
                CALL HCO_ERROR( MSG, RC )
                RETURN
             ENDIF

             ! Open file
             CALL NC_OPEN ( TRIM(srcFile2), ncLun2 )

             ! Define time stamp to be read.
             CALL GET_TIMEIDX ( HcoState,  Lct,               &
                                ncLun2,    tidx1,    tidx2,   &
                                wgt1,      wgt2,     oYMDhm2, &
                                YMDhmb,    YMDhm1,   RC,      &
                                Year=iYear                    )
             IF ( RC /= HCO_SUCCESS ) THEN
                 CALL HCO_ERROR( 'ERROR 9', RC, THISLOC=LOC )
                 RETURN
             ENDIF

             ! Do not perform weights
             wgt1  = -1.0_sp
             wgt2  = -1.0_sp

             ! Read data and write into array ncArr2
             CALL NC_READ_ARR( fID     = ncLun2,             &
                               ncVar   = Lct%Dct%Dta%ncPara, &
                               lon1    = 1,                  &
                               lon2    = nlon,               &
                               lat1    = 1,                  &
                               lat2    = nlat,               &
                               lev1    = lev1,               &
                               lev2    = lev2,               &
                               time1   = tidx1,              &
                               time2   = tidx2,              &
                               ncArr   = ncArr2,             &
                               varUnit = thisUnit,           &
                               wgt1    = wgt1,               &
                               wgt2    = wgt2,               &
                               MissVal = HCO_MISSVAL,        &
                               ArbIdx  = ArbIdx,             &
                               RC      = NCRC                 )
             IF ( NCRC /= 0 ) THEN
                CALL HCO_ERROR( 'NC_READ_ARRAY (3)', RC )
                RETURN
             ENDIF

             ! Eventually fissing values
!!!             CALL CheckMissVal ( Lct, ncArr2 )

             ! Add all values to ncArr
             ncArr = ncArr + ncArr2

             ! Cleanup
             IF ( ASSOCIATED(ncArr2) ) DEALLOCATE(ncArr2)

             ! Close file
             CALL NC_CLOSE ( ncLun2 )

          ENDDO !iYear

          ! Now calculate average
          ncArr = ncArr / REAL(nYears,sp)

          ! Verbose
          IF ( HcoState%amIRoot .AND. HCO_IsVerb(HcoState%Config%Err,1) ) THEN
             WRITE(MSG,110) TRIM(Lct%Dct%cName), Yr1, Yr2
             CALL HCO_MSG(HcoState%Config%Err,MSG)
          ENDIF
 110      FORMAT( 'Field ', a, ': Average data over years ', I4.4, ' to ', I4.4 )

       ENDIF !nYears>1
    ENDIF !Averaging

    !-----------------------------------------------------------------
    ! Convert to HEMCO units
    ! HEMCO data are all in kg/m2/s for fluxes and kg/m3 for
    ! concentrations. Unit conversion is performed based on the
    ! unit on the input file and the srcUnit attribute given in the
    ! configuration file. By default, HEMCO will attempt to convert
    ! the units found in the input file to the standard quantities
    ! for mass (kg), area (m2 or m3), and time (s). For instance,
    ! g/cm2/hr will be converted to kg/m2/s. The exceptions to this
    ! rule are:
    ! 1. If srcUnit is set to '1', the input data are expected to
    !    be unitless. If the units string on the input file is none
    !    of the units recognized by HEMCO as unitless, an error is
    !    returned if the unit tolerance setting is set to zero, or
    !    a warning is prompted if unit tolerance is greater than zero.
    ! 2. If srcUnit is set to 'count', no unit conversion is performed
    !    and data will be treated as 'index' data, e.g. regridding will
    !    preserve the absolute values.
    !-----------------------------------------------------------------

    ! If OrigUnit is set to wildcard character: use unit from source file
    IF ( TRIM(Lct%Dct%Dta%OrigUnit) == &
         TRIM(HCO_GetOpt(HcoState%Config%ExtList,'Wildcard')) ) THEN
       Lct%Dct%Dta%OrigUnit = TRIM(thisUnit)
    ENDIF

    ! If OrigUnit is set to '1' or to 'count', perform no unit
    ! conversion.
    IF ( HCO_IsUnitLess(Lct%Dct%Dta%OrigUnit)  .OR. &
         HCO_IsIndexData(Lct%Dct%Dta%OrigUnit)       ) THEN

       ! Check if file unit is also unitless. This will return 0 for
       ! unitless, 1 for HEMCO emission unit, 2 for HEMCO conc. unit,
       ! -1 otherwise.
       Flag = HCO_UNIT_SCALCHECK( thisUnit )

       ! Return with error if: (1) thisUnit is recognized as HEMCO unit and
       ! unit tolerance is set to zero; (2) thisUnit is neither unitless nor
       ! a HEMCO unit and unit tolerance is set to zero or one.
       ! The unit tolerance is defined in the configuration file.
       IF ( Flag /= 0 .AND. UnitTolerance == 0 ) THEN
          MSG = 'Illegal unit: ' // TRIM(thisUnit) // '. File: ' // &
                TRIM(srcFile)
          CALL HCO_ERROR( MSG, RC )
          RETURN
       ENDIF

       ! Prompt a warning if thisUnit is not recognized as unitless.
       IF ( Flag /= 0 ) THEN
          MSG = 'Data is treated as unitless, but file attribute suggests ' // &
                'it is not: ' // TRIM(thisUnit) // '. File: ' // TRIM(srcFile)
          CALL HCO_WARNING( HcoState%Config%Err, MSG, RC, WARNLEV=1 )
       ENDIF

       ! Verbose mode
       IF ( HCO_IsVerb(HcoState%Config%Err,2) ) THEN
          WRITE(MSG,*) 'Based on srcUnit attribute (', TRIM(Lct%Dct%Dta%OrigUnit), &
                       '), no unit conversion is performed.'
          CALL HCO_MSG(HcoState%Config%Err,MSG)
       ENDIF

    ! Convert to HEMCO units in all other cases.
    ELSE

       ! For zero unit tolerance, make sure that thisUnit matches
       ! with unit set in configuration file. For higher unit
       ! tolerances, prompt a level 3 warning.
       IF ( TRIM(Lct%Dct%Dta%OrigUnit) /= TRIM(thisUnit) ) THEN
          MSG = 'File units do not match: ' // TRIM(thisUnit) // &
                ' vs. ' // TRIM(Lct%Dct%Dta%OrigUnit)    // &
                '. File: ' // TRIM(srcFile)

          IF ( UnitTolerance == 0 ) THEN
             CALL HCO_ERROR( MSG, RC )
             RETURN
          ELSE
             CALL HCO_WARNING( HcoState%Config%Err, MSG, RC, WARNLEV=3 )
          ENDIF
       ENDIF

       ! Shadow species properties needed for unit conversion
       HcoID = Lct%Dct%HcoID
       IF ( HcoID > 0 ) THEN
          MW_g = HcoState%Spc(HcoID)%MW_g
       ELSE
          MW_g = -999.0_hp
       ENDIF

       ! Now convert to HEMCO units. This attempts to convert mass,
       ! area/volume and time to HEMCO standards (kg, m2/m3, s).
       ncYr  = FLOOR( MOD( oYMDhm1, 1.0e12_dp ) / 1.0e8_dp )
       ncMt  = FLOOR( MOD( oYMDhm1, 1.0e8_dp  ) / 1.0e6_dp )

       IF ( ncYr == 0 ) THEN
          CALL HcoClock_Get( HcoState%Clock, cYYYY = ncYr, RC=RC )
          IF ( RC /= HCO_SUCCESS ) THEN
              CALL HCO_ERROR( 'ERROR 10', RC, THISLOC=LOC )
              RETURN
          ENDIF
       ENDIF
       IF ( ncMt == 0 ) THEN
          CALL HcoClock_Get( HcoState%Clock, cMM   = ncMt, RC=RC )
          IF ( RC /= HCO_SUCCESS ) THEN
              CALL HCO_ERROR( 'ERROR 11', RC, THISLOC=LOC )
              RETURN
          ENDIF
       ENDIF

       ! Verbose mode
       IF ( HcoState%amIRoot .and. HCO_IsVerb(HcoState%Config%Err,3) ) THEN
          WRITE(MSG,*) 'Unit conversion settings: '
          CALL HCO_MSG(HcoState%Config%Err,MSG)
          WRITE(MSG,*) '- Year, month        : ', ncYr, ncMt
          CALL HCO_MSG(HcoState%Config%Err,MSG)
       ENDIF

       CALL HCO_UNIT_CHANGE(                 &
            HcoConfig     = HcoState%Config, &
            Array         = ncArr,           &
            Units         = thisUnit,        &
            MW            = MW_g,            &
            YYYY          = ncYr,            &
            MM            = ncMt,            &
            AreaFlag      = AreaFlag,        &
            TimeFlag      = TimeFlag,        &
            FACT          = UnitFactor,      &
            RC            = RC                )
       IF ( RC /= HCO_SUCCESS ) THEN
          MSG = 'Cannot convert units for ' // TRIM(Lct%Dct%cName)
          CALL HCO_ERROR( MSG , RC )
          RETURN
       ENDIF

       ! Verbose mode
       IF ( UnitFactor /= 1.0_hp ) THEN
          IF ( HCO_IsVerb(HcoState%Config%Err,1) ) THEN
             WRITE(MSG,*) 'Data was in units of ', TRIM(thisUnit), &
                          ' - converted to HEMCO units by applying ', &
                          'scale factor ', UnitFactor
             CALL HCO_MSG(HcoState%Config%Err,MSG)
          ENDIF
       ELSE
          IF ( HCO_IsVerb(HcoState%Config%Err,2) ) THEN
             WRITE(MSG,*) 'Data was in units of ', TRIM(thisUnit), &
                          ' - unit conversion factor is ', UnitFactor
             CALL HCO_MSG(HcoState%Config%Err,MSG)
          ENDIF
       ENDIF

       ! Check for valid unit combinations, i.e. emissions must be kg/m2/s,
       ! concentrations kg/m3. Eventually multiply by emission time step
       ! or divide by area to obtain those values.

       ! Concentration data
       IF ( AreaFlag == 3 .AND. TimeFlag == 0 ) THEN
          Lct%Dct%Dta%IsConc = .TRUE.

       ! If concentration data is per second (kg/m3/s), multiply by emission
       ! time step to get concentration (kg/m3).
       ELSEIF ( AreaFlag == 3 .AND. TimeFlag == 1 ) THEN
          Lct%Dct%Dta%IsConc = .TRUE.

          ncArr = ncArr * HcoState%TS_EMIS
          MSG = 'Data converted from kg/m3/s to kg/m3: ' // &
                TRIM(Lct%Dct%cName) // ': ' // TRIM(thisUnit)
          CALL HCO_WARNING( HcoState%Config%Err, MSG, RC, WARNLEV=1 )

       ! Unitless data
       ELSEIF ( AreaFlag == -1 .AND. TimeFlag == -1 ) THEN
          ! nothing do to

       ! Emission data
       ELSEIF ( AreaFlag == 2 .AND. TimeFlag == 1 ) THEN
          ! nothing do to

       ! Emission data that is not per time (kg/m2): convert to kg/m2/s
       ELSEIF ( AreaFlag == 2 .AND. TimeFlag == 0 ) THEN
          ncArr = ncArr / HcoState%TS_EMIS
          MSG = 'Data converted from kg/m2 to kg/m2/s: ' // &
                TRIM(Lct%Dct%cName) // ': ' // TRIM(thisUnit)
          CALL HCO_WARNING( HcoState%Config%Err, MSG, RC, WARNLEV=1 )

       ! Emission data that is not per area (i.e. kg/s) needs to be converted
       ! to per area manually.
       ELSEIF ( AreaFlag == 0 .AND. TimeFlag == 1 ) THEN

          ! Get lat edges: those are read from file if possible, otherwise
          ! calculated from the lat midpoints.
          ! ==> Sine of lat is needed. Do conversion right here.
          CALL NC_GET_GRID_EDGES ( ncLun, 2, LatMid,   nlat, &
                                   LatEdge,  nlatEdge, NCRC   )
          IF ( NCRC /= 0 ) THEN
             MSG = 'Cannot read lat edge of ' // TRIM(srcFile)
             CALL HCO_ERROR( MSG, RC )
             RETURN
          ENDIF

          ! Now normalize data by area calculated from lat edges.
          CALL NORMALIZE_AREA( HcoState, ncArr,   nlon, &
                               LatEdge,  srcFile, RC     )
          IF ( RC /= HCO_SUCCESS ) THEN
              CALL HCO_ERROR( 'ERROR 12', RC, THISLOC=LOC )
              RETURN
          ENDIF

       ! All other combinations are invalid
       ELSE
          MSG = 'Unit must be unitless, emission or concentration: ' // &
                TRIM(Lct%Dct%cName) // ': ' // TRIM(thisUnit)
          CALL HCO_ERROR( MSG, RC )
          RETURN
       ENDIF
    ENDIF ! Unit conversion

    !-----------------------------------------------------------------
    ! Get horizontal grid edges
    !-----------------------------------------------------------------

    ! Get longitude edges and make sure they are steadily increasing.
    CALL NC_GET_GRID_EDGES ( ncLun, 1, LonMid,   nlon, &
                             LonEdge,  nlonEdge, NCRC   )
    IF ( NCRC /= 0 ) THEN
       MSG = 'Cannot read lon edge of ' // TRIM(srcFile)
       CALL HCO_ERROR( MSG, RC )
       RETURN
    ENDIF
    CALL HCO_ValidateLon( HcoState, nlonEdge, LonEdge, RC )
    IF ( RC /= HCO_SUCCESS ) THEN
        CALL HCO_ERROR( 'ERROR 13', RC, THISLOC=LOC )
        RETURN
    ENDIF

    ! Get latitude edges (only if they have not been read yet
    ! for unit conversion)
    IF ( .NOT. ASSOCIATED( LatEdge ) ) THEN
       CALL NC_GET_GRID_EDGES ( ncLun, 2, LatMid,   nlat, &
                                LatEdge,  nlatEdge, NCRC   )
       IF ( NCRC /= 0 ) THEN
          MSG = 'Cannot read lat edge of ' // TRIM(srcFile)
          CALL HCO_ERROR( MSG, RC )
          RETURN
       ENDIF
    ENDIF

    !-----------------------------------------------------------------
    ! Determine regridding algorithm to be applied: use NCREGRID from
    ! MESSy only if we need to regrid vertical levels. For all other
    ! fields, use the much faster map_a2a.
    ! Perform no vertical regridding if the vertical levels are model
    ! levels. Model levels are assumed to start at the surface, i.e.
    ! the first input level must correspond to the surface level. The
    ! total number of vertical levels must not match the number of
    ! vertical levels on the simulation grid. Data is not extrapolated
    ! beyond the existing levels.
    ! Vertical regridding based on NCREGRID will always map the input
    ! data onto the entire simulation grid (no extrapolation beyond
    ! the vertical input coordinates).
    ! Index-based remapping can currently not be done with the MESSy
    ! routines, i.e. it is not possible to vertically regrid index-
    ! based data.
    !-----------------------------------------------------------------

    UseMESSy = .FALSE.
    IF ( nlev > 1 .AND. .NOT. IsModelLevel ) THEN
       UseMESSy = .TRUE.
    ENDIF

#if defined( MODEL_CESM ) || defined( MODEL_WRF )
    ! If in WRF or the CESM environment, the vertical grid is arbitrary.
    ! MESSy regridding ALWAYS has to be used.
    IF ( nlev > 1 ) THEN
      UseMESSy = .TRUE.

      IF ( HCO_IsVerb(HcoState%Config%Err,2) ) THEN
        WRITE(MSG,*) '  ==> WRF/CESM: Always forcing MESSy regridding for number of verticals', nlev, IsModelLevel
        CALL HCO_MSG(HcoState%Config%Err,MSG)
      ENDIF
    ENDIF
#endif

    IF ( HCO_IsIndexData(Lct%Dct%Dta%OrigUnit) .AND. UseMESSy ) THEN
       MSG = 'Cannot do MESSy regridding for index data: ' // &
             TRIM(srcFile)
       CALL HCO_ERROR( MSG, RC )
       RETURN
    ENDIF

    !-----------------------------------------------------------------
    ! Use MESSy regridding
    !-----------------------------------------------------------------
    IF ( UseMESSy ) THEN
       IF ( HCO_IsVerb(HcoState%Config%Err,2) ) THEN
          WRITE(MSG,*) '  ==> Use MESSy regridding (NCREGRID)'
          CALL HCO_MSG(HcoState%Config%Err,MSG)
       ENDIF

#if !defined( MODEL_CESM ) && !defined( MODEL_WRF )
       ! If we do MESSy regridding, we can only do one time step
       ! at a time at the moment!
       IF ( tidx1 /= tidx2 ) THEN
          MSG = 'Cannot do MESSy regridding for more than one time step; ' &
                // TRIM(srcFile)
          CALL HCO_ERROR( MSG, RC )
          RETURN
       ENDIF

       ! Note: This seems to be a soft restriction - removing this
       ! does not conflict with MESSy regridding. Need to check (hplin, 5/30/20)
       ! This has to be used for WRF-GC and CESM so ifdefd out
#endif

#if defined( MODEL_WRF ) || defined( MODEL_CESM )
       !--------------------------------------------------------------
       ! Eventually get sigma levels
       ! For files that have hardcoded GEOS-Chem "index"-based levels,
       ! translate these levels back into a sigma representation
       ! of the GEOS-Chem levels (sigma = p/ps on INTERFACE)
       !
       ! There are caveats with this. This is essentially a copy of the
       ! hardcoded hPa lists from
       ! http://wiki.seas.harvard.edu/geos-chem/index.php/GEOS-Chem_vertical_grids
       ! hard-coded by hand, and we only assume that the data is either
       ! 47-levels or 72-levels.
       !
       ! Parse the 72 list using regex like so: ^ ?\d{1,2} then remove the lines
       ! Then you have the 73 edges.
       !
       ! psfc = PEDGE(0) = 1013.250 hPa
       !
       ! Ported from the original WRF-GC implementation (hplin, 5/27/20)
       !--------------------------------------------------------------
       IF ( nlev > 1 .AND. IsModelLevel ) THEN
          IF ( HCO_IsVerb(HcoState%Config%Err,2) ) THEN
            WRITE(MSG,*) '  ==> WRF/CESM: Writing in fixed sigma coordinates for GEOS-Chem levels', nlon, nlat
            CALL HCO_MSG(HcoState%Config%Err,MSG)
          ENDIF

          ALLOCATE(SigEdge(nlon, nlat, nlev))

          DO I = 1, nlon
             DO J = 1, nlat
               ! Fill with pre-defined, hard coded sigma levels computed.
               SigEdge(I, J, :) = GC_72_EDGE_SIGMA(1:nlev)
             ENDDO
           ENDDO
       ENDIF
#endif

       !--------------------------------------------------------------
       ! Eventually get sigma levels
       ! Vertical regridding is performed on sigma interface levels:
       ! sigma(i,j,l) = p(i,j,l) / ps(i,j)
       ! NC_Get_Sigma_Levels attempts to create the sigma levels from
       ! the content of the netCDF file.
       ! For now, it is assumed that all input data is on vertical
       ! mid-point levels, and the interface values are calculated
       ! by linear interpolation of the mid-point values in a second
       ! step.
       ! For model levels, the sigma levels don't need to be known
       ! as vertical interpolation will be done based on subroutine
       ! ModelLev_Interpolate (within HCO_MESSY_REGRID).
       !--------------------------------------------------------------
       IF ( nlev > 1 .AND. .NOT. IsModelLevel ) THEN

          ! Get sigma levels
          CALL NC_Get_Sigma_Levels ( fID     = ncLun,   &
                                     ncFile  = srcFile, &
                                     levName = LevName, &
                                     lon1    = 1,       &
                                     lon2    = nlon,    &
                                     lat1    = 1,       &
                                     lat2    = nlat,    &
                                     lev1    = 1,       &
                                     lev2    = nlev,    &
                                     time    = tidx1,   &
                                     SigLev  = SigLev,  &
                                     Dir     = dir,     &
                                     RC      = NCRC      )
          IF ( NCRC /= 0 ) THEN
             CALL HCO_ERROR( 'Cannot read sigma levels of '//TRIM(srcFile), RC )
             RETURN
          ENDIF

          ! Interpolate onto edges
          CALL SigmaMidToEdges ( HcoState, SigLev, SigEdge, RC )
          IF ( RC /= HCO_SUCCESS ) THEN
              CALL HCO_ERROR( 'ERROR 14', RC, THISLOC=LOC )
              RETURN
          ENDIF

          ! Sigma levels are not needed anymore
          IF ( ASSOCIATED(SigLev) ) DEALLOCATE(SigLev)

          !-----------------------------------------------------------
          ! Flip vertical axis if positive axis is 'down', i.e. level
          ! index 1 is the top of the atmosphere
          !-----------------------------------------------------------
          IF ( dir == -1 ) THEN
             SigEdge(:,:,:  ) = SigEdge(:,:,nlev+1:1:-1  )
             NcArr  (:,:,:,:) = NcArr  (:,:,nlev  :1:-1,:)
          ENDIF

       ENDIF ! nlev>1

#if defined( MODEL_WRF ) || defined( MODEL_CESM )
       ! Input data is "never" on model levels because model levels can change! (hplin, 5/29/20)
       IsModelLevel = .false.
#endif

       ! Now do the regridding
       CALL HCO_MESSY_REGRID ( HcoState,  NcArr,                 &
                               LonEdge,   LatEdge,      SigEdge, &
                               Lct,       IsModelLevel, RC        )
       IF ( RC /= HCO_SUCCESS ) THEN
           CALL HCO_ERROR( 'ERROR 15', RC, THISLOC=LOC )
           RETURN
       ENDIF

       ! Cleanup
       IF ( ASSOCIATED(SigEdge) ) DEALLOCATE(SigEdge)

    !-----------------------------------------------------------------
    ! Use map_a2a regridding
    !-----------------------------------------------------------------
    ELSE
       IF ( HCO_IsVerb(HcoState%Config%Err,2) ) THEN
          WRITE(MSG,*) '  ==> Use map_a2a regridding'
          CALL HCO_MSG(HcoState%Config%Err,MSG)
       ENDIF

       CALL REGRID_MAPA2A ( HcoState, NcArr, LonEdge, LatEdge, Lct, RC )
       IF ( RC /= HCO_SUCCESS ) THEN
           CALL HCO_ERROR( 'ERROR 16', RC, THISLOC=LOC )
           RETURN
       ENDIF

    ENDIF

    !-----------------------------------------------------------------
    ! Add to diagnostics (if it exists)
    !-----------------------------------------------------------------
    IF ( HcoState%Options%Field2Diagn ) THEN
       IF ( Lct%Dct%Dta%SpaceDim == 3 .AND. ASSOCIATED(Lct%Dct%Dta%V3) ) THEN
          IF ( ASSOCIATED(Lct%Dct%Dta%V3(1)%Val) ) THEN
             CALL Diagn_Update ( HcoState, cName=TRIM(Lct%Dct%cName), &
                                 Array3D=Lct%Dct%Dta%V3(1)%Val, COL=-1, RC=RC )
             IF ( RC /= HCO_SUCCESS ) THEN
                 CALL HCO_ERROR( 'ERROR 17', RC, THISLOC=LOC )
                 RETURN
             ENDIF
          ENDIF
       ELSEIF ( Lct%Dct%Dta%SpaceDim == 2 .AND. ASSOCIATED(Lct%Dct%Dta%V2) ) THEN
          IF ( ASSOCIATED(Lct%Dct%Dta%V2(1)%Val) ) THEN
             CALL Diagn_Update ( HcoState, cName=TRIM(Lct%Dct%cName), &
                                 Array2D=Lct%Dct%Dta%V2(1)%Val, COL=-1, RC=RC )
             IF ( RC /= HCO_SUCCESS ) THEN
                 CALL HCO_ERROR( 'ERROR 18', RC, THISLOC=LOC )
                 RETURN
             ENDIF
          ENDIF
       ENDIF
    ENDIF

    !-----------------------------------------------------------------
    ! Cleanup and leave
    !-----------------------------------------------------------------
    IF ( ASSOCIATED ( ncArr   ) ) DEALLOCATE ( ncArr   )
    IF ( ASSOCIATED ( LonMid  ) ) DEALLOCATE ( LonMid  )
    IF ( ASSOCIATED ( LatMid  ) ) DEALLOCATE ( LatMid  )
    IF ( ASSOCIATED ( LevMid  ) ) DEALLOCATE ( LevMid  )
    IF ( ASSOCIATED ( LonEdge ) ) DEALLOCATE ( LonEdge )
    IF ( ASSOCIATED ( LatEdge ) ) DEALLOCATE ( LatEdge )

    ! Return w/ success
    CALL HCO_LEAVE ( HcoState%Config%Err,  RC )

  END SUBROUTINE HCOIO_Read
!EOC
!------------------------------------------------------------------------------
!                   Harmonized Emissions Component (HEMCO)                    !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: HCOIO_CloseAll
!
! !DESCRIPTION: Subroutine HCOIO\_CloseAll makes sure that there is no open
! netCDF file left in the stream.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE HCOIO_CloseAll( HcoState, RC )
!
! !USES:
!
    USE HCO_Ncdf_Mod,   ONLY : NC_CLOSE
!
! !INPUT PARAMTERS:
!
    TYPE(HCO_State), POINTER          :: HcoState    ! HEMCO state
!
! !INPUT/OUTPUT PARAMETERS:
!
    INTEGER,          INTENT(INOUT)   :: RC
!
! !REVISION HISTORY:
!  24 Mar 2016 - C. Keller: Initial version
!  See https://github.com/geoschem/hemco for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    !======================================================================
    ! HCOIO_CloseAll begins here
    !======================================================================
    IF ( HcoState%ReadLists%FileLun > 0 ) THEN
       CALL NC_CLOSE( HcoState%ReadLists%FileLun )
       HcoState%ReadLists%FileLun = -1
    ENDIF

    ! Return w/ success
    RC = HCO_SUCCESS

  END SUBROUTINE HCOIO_CloseAll
!EOC
END MODULE HCOIO_Read_Mod
#endif
