!------------------------------------------------------------------------------
!                   Harmonized Emissions Component (HEMCO)                    !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: hcoio_util_mod.F90
!
! !DESCRIPTION: Module HCOIO\_Util\_Mod contains utility functions
! for use in data processing including file reading, unit conversions,
! and regridding.
!\\
!\\
! !INTERFACE:
!
MODULE HCOIO_Util_Mod
!
! !USES:
!
  USE HCO_Types_Mod
  USE HCO_Error_Mod
  USE HCO_CharTools_Mod
  USE HCO_State_Mod,       ONLY : Hco_State

  IMPLICIT NONE
  PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!
#if !defined(ESMF_)
  PUBLIC :: GET_TIMEIDX
  PUBLIC :: Check_AvailYMDhm
  PUBLIC :: prefYMDhm_Adjust
  PUBLIC :: Set_tIdx2
  PUBLIC :: IsClosest
  PUBLIC :: GetIndex2Interp
  PUBLIC :: GetWeights
  PUBLIC :: YMDhm2hrs
  PUBLIC :: Normalize_Area
  PUBLIC :: SrcFile_Parse
  PUBLIC :: SigmaMidToEdges
  PUBLIC :: CheckMissVal
  PUBLIC :: GetArbDimIndex
#endif
  PUBLIC :: HCOIO_ReadOther
  PUBLIC :: HCOIO_ReadCountryValues
  PUBLIC :: HCOIO_ReadFromConfig
  PUBLIC :: GetDataVals
  PUBLIC :: GetSliceIdx
  PUBLIC :: FillMaskBox
  PUBLIC :: ReadMath
!
! !REVISION HISTORY:
!  12 Jun 2020 - E. Lundgren - Initial version, created from subset of
!                              hcoio_util_mod.F90
!  See https://github.com/geoschem/hemco for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !DEFINED PARAMETERS
!
  ! Parameter used for difference testing of floating points
  REAL(dp), PRIVATE, PARAMETER :: EPSILON = 1.0e-5_dp

CONTAINS
!EOC
#if !defined( ESMF_ )
!------------------------------------------------------------------------------
!                   Harmonized Emissions Component (HEMCO)                    !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Get_TimeIdx
!
! !DESCRIPTION: Returns the lower and upper time slice index (tidx1
! and tidx2, respectively) to be read. These values are determined
! based upon the time slice information extracted from the netCDF file,
! the time stamp settings set in the config. file, and the current
! simulation date.
!\\
!\\
! Return arguments wgt1 and wgt2 denote the weights to be given to
! the two time slices. This is only of relevance for data that shall
! be interpolated between two (not necessarily consecutive) time slices.
! In all other cases, the returned weights are negative and will be
! ignored.
!\\
!\\
! Also returns the time slice year and month, as these values may be
! used for unit conversion.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE GET_TIMEIDX( HcoState,  Lct,               &
                          ncLun,     tidx1,    tidx2,   &
                          wgt1,      wgt2,     oYMDhm,  &
                          YMDhm,     YMDhm1,   RC,      &
                          Year )
!
! !USES:
!
    USE HCO_Ncdf_Mod,  ONLY : NC_Read_Time_YYYYMMDDhhmm
    USE HCO_tIdx_Mod,  ONLY : HCO_GetPrefTimeAttr
!
! !INPUT PARAMETERS:
!
    TYPE(HCO_State),  POINTER                  :: HcoState  ! HcoState object
    TYPE(ListCont),   POINTER                  :: Lct       ! List container
    INTEGER,          INTENT(IN   )            :: ncLun     ! open ncLun
    INTEGER,          INTENT(IN   ), OPTIONAL  :: Year      ! year to be used
!
! !OUTPUT PARAMETERS:
!
    INTEGER,          INTENT(  OUT)            :: tidx1  ! lower time idx
    INTEGER,          INTENT(  OUT)            :: tidx2  ! upper time idx
    REAL(sp),         INTENT(  OUT)            :: wgt1   ! weight to tidx1
    REAL(sp),         INTENT(  OUT)            :: wgt2   ! weight to tidx2
    REAL(dp),         INTENT(  OUT)            :: oYMDhm ! preferred time slice
    REAL(dp),         INTENT(  OUT)            :: YMDhm  ! selected time slice
    REAL(dp),         INTENT(  OUT)            :: YMDhm1 ! 1st time slice in file
!
! !INPUT/OUTPUT PARAMETERS:
!
    INTEGER,          INTENT(INOUT)            :: RC
!
! !REVISION HISTORY:
!  13 Mar 2013 - C. Keller - Initial version
!  See https://github.com/geoschem/hemco for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOcAL VARIABLES:
!
    CHARACTER(LEN=255)    :: MSG, LOC
    CHARACTER(LEN=1023)   :: MSG_LONG
    INTEGER               :: tidx1a
    INTEGER               :: nTime,  T, CNT, NCRC
    INTEGER               :: prefYr, prefMt, prefDy, prefHr, prefMn
    INTEGER               :: refYear
    REAL(dp)              :: origYMDhm, prefYMDhm
    REAL(dp),   POINTER   :: availYMDhm(:)
    LOGICAL               :: ExitSearch
    LOGICAL               :: verb

    !=================================================================
    ! GET_TIMEIDX begins here
    !=================================================================

    ! Initialize
    LOC         = 'GET_TIMEIDX (HCOIO_UTIL_MOD.F90)'

    ! Officially enter Get_TimeIdx
    CALL HCO_ENTER( HcoState%Config%Err, LOC, RC )
    IF ( RC /= HCO_SUCCESS ) THEN
        CALL HCO_ERROR( 'ERROR 0', RC, THISLOC=LOC )
        RETURN
    ENDIF
    verb = HCO_IsVerb(HcoState%Config%Err,3)

    ! Initialize local variables for safety's sake
    nTime      =  0
    cnt        =  0
    prefYr     =  0
    prefMt     =  0
    prefDy     =  0
    prefHr     =  0
    prefMn     =  0
    refYear    =  0
    origYMDhm  =  0
    prefYMDhm  =  0
    tidx1      =  0
    tidx2      =  0
    tidx1a     =  0
    wgt1       = -1.0_sp
    wgt2       = -1.0_sp
    oYMDhm     =  0.0_dp
    YMDhm      =  0.0_dp
    YMDhm1     =  0.0_dp
    ExitSearch = .FALSE.
    availYMDhm => NULL()

    ! ----------------------------------------------------------------
    ! Extract netCDF time slices (YYYYMMDDhhmm)
    ! ----------------------------------------------------------------
    CALL NC_READ_TIME_YYYYMMDDhhmm( ncLun, nTime,    availYMDhm,  &
                                    refYear=refYear, RC=NCRC     )
    IF ( NCRC /= 0 ) THEN
       CALL HCO_ERROR( 'NC_READ_TIME_YYYYMMDDhhmm', RC )
       RETURN
    ENDIF

    ! Return warning if netCDF reference year prior to 1901: it seems
    ! like there are some problems with that and the time slices can be
    ! off by one day!
    IF ( (refYear <= 1900) .AND. (nTime > 0) ) THEN
       MSG = 'ncdf reference year is prior to 1901 - ' // &
            'time stamps may be wrong!'
       CALL HCO_WARNING ( HcoState%Config%Err, MSG, RC, WARNLEV=1 )
    ENDIF

    ! verbose mode
    IF ( verb ) THEN
       write(MSG,*) 'Number of time slices found: ', nTime
       CALL HCO_MSG(HcoState%Config%Err,MSG)
       IF ( nTime > 0 ) THEN
          write(MSG,*) 'Time slice range : ', &
                       availYMDhm(1), availYMDhm(nTime)
          CALL HCO_MSG(HcoState%Config%Err,MSG)
       ENDIF
    ENDIF

    ! ----------------------------------------------------------------
    ! Select time slices to read
    ! ----------------------------------------------------------------

    ! ----------------------------------------------------------------
    ! Get preferred time stamp to read based upon the specs set in the
    ! config. file.
    ! This can return value -1 for prefHr, indicating that all
    ! corresponding time slices shall be read.
    ! This call will return -1 for all date attributes if the
    ! simulation date is outside of the data range given in the
    ! configuration file.
    ! ----------------------------------------------------------------
    CALL HCO_GetPrefTimeAttr ( HcoState, Lct, &
                               prefYr, prefMt, prefDy, prefHr, prefMn, RC )
    IF ( RC /= HCO_SUCCESS ) THEN
       MSG = &
         'Error encountered in HCO_GetPrefTimeAttr for ' // TRIM(Lct%Dct%cName)
       CALL HCO_ERROR( MSG, RC )
       IF ( ASSOCIATED(availYMDhm) ) THEN
          DEALLOCATE(availYMDhm)
          availYMDhm => NULL()
       ENDIF
       RETURN
    ENDIF

    ! Eventually force preferred year to passed value
    IF ( PRESENT(Year) ) prefYr = Year

    ! Check if we are outside of provided range
    IF ( prefYr < 0 .OR. prefMt < 0 .OR. prefDy < 0 ) THEN

       ! This should only happen for 'range' data
       IF ( Lct%Dct%Dta%CycleFlag /= HCO_CFLAG_RANGE ) THEN
          MSG = 'Cannot get preferred datetime for ' // TRIM(Lct%Dct%cName)
          CALL HCO_ERROR( MSG, RC )
          IF ( ASSOCIATED(availYMDhm) ) THEN
             DEALLOCATE(availYMDhm)
             availYMDhm => NULL()
          ENDIF
          RETURN
       ENDIF

       ! If this part of the code gets executed, the data associated
       ! with this container shall not be used at the current date.
       ! To do so, set the time indeces to -1 and leave right here.
       tidx1 = -1
       tidx2 = -1

       ! Leave w/ success
       CALL HCO_LEAVE( HcoState%Config%Err,  RC )
       RETURN
    ENDIF

    ! origYMDhm is the preferred datetime. Store into shadow variable
    ! prefYMDhm. prefYMDhm may be adjusted if origYMDhm is outside of the
    ! netCDF datetime range.
    ! Now put origYMDhm, prefYMDhm in YYYYMMDDhhmm format (bmy, 4/10/17)
    origYMDhm = ( DBLE(      prefYr      ) * 1.0e8_dp ) + &
                ( DBLE(      prefMt      ) * 1.0e6_dp ) + &
                ( DBLE(      prefDy      ) * 1.0e4_dp ) + &
                ( DBLE( MAX( prefHr, 0 ) ) * 1.0e2_dp ) + &
                ( DBLE( MAX( prefMn, 0 ) )          )
    prefYMDhm = origYMDhm

    ! verbose mode
    IF ( verb ) THEN
       write(MSG,'(A30,f14.0)') 'preferred datetime: ', prefYMDhm
       CALL HCO_MSG(HcoState%Config%Err,MSG)
    ENDIF

    ! ================================================================
    ! Case 1: Only one time slice available.
    ! ================================================================
    IF ( nTime == 1 ) THEN
       tidx1 = 1
       tidx2 = 1

    ! ================================================================
    ! Case 2: More than one time slice available. Determine lower
    ! and upper time slice index from file & HEMCO settings.
    ! ================================================================
    ELSEIF ( nTime > 1 ) THEN

       ! Init
       tidx1   = -1
       tidx2   = -1

       ! -------------------------------------------------------------
       ! Check if preferred datetime prefYMDhm is within the range
       ! available time slices, e.g. it falls within the interval
       ! of availYMDhm. In this case, set tidx1 to the index of the
       ! closest time slice that is not in the future.
       ! -------------------------------------------------------------
       CALL Check_AvailYMDhm ( Lct, nTime, availYMDhm, prefYMDhm, tidx1a )

       ! -------------------------------------------------------------
       ! Check if we need to continue search. Even if the call above
       ! returned a time slice, it may be possible to continue looking
       ! for a better suited time stamp. This is only the case if
       ! there are discontinuities in the time stamps, e.g. if a file
       ! contains monthly data for 2005 and 2020. In that case, the
       ! call above would return the index for Dec 2005 for any
       ! simulation date between 2005 and 2010 (e.g. July 2010),
       ! whereas it makes more sense to use July 2005 (and eventually
       ! interpolate between the July 2005 and July 2020 data).
       ! The IsClosest command checks if there are any netCDF time
       ! stamps (prior to the selected one) that are closer to each
       ! other than the difference between the preferred time stamp
       ! prefYMDhm and the currently selected time stamp
       ! availYMDhm(tidx1a). In that case, it continues the search by
       ! updating prefYMDhm so that it falls within the range of the
       ! 'high-frequency' interval.
       ! -------------------------------------------------------------
       ExitSearch = .FALSE.
       IF ( Lct%Dct%Dta%CycleFlag == HCO_CFLAG_EXACT ) THEN
          ExitSearch = .TRUE.
       ELSE IF ( tidx1a > 0 ) THEN
          ExitSearch = IsClosest( prefYMDhm, availYMDhm, nTime, tidx1a )
       ENDIF

       ! When using the interpolation flag, use the first or last timestep
       ! when outside of the available date range
       IF ( Lct%Dct%Dta%CycleFlag == HCO_CFLAG_INTER .and. tidx1a < 0 ) THEN
          IF ( prefYMDhm < availYMDhm(1) ) THEN
             tidx1a = 1
          ELSE IF ( prefYMDhm > availYMDhm(nTime) ) THEN
             tidx1a = nTime
          ENDIF
       ENDIF

       ! Do not continue search if data is to be interpolated and is
       ! not discontinuous (mps, 10/23/19)
       IF ( Lct%Dct%Dta%CycleFlag == HCO_CFLAG_INTER .and. &
            .not. Lct%Dct%Dta%Discontinuous ) THEN
          ExitSearch = .TRUE.
       ENDIF

       ! Write to tidx1 if this is the best match.
       IF ( ExitSearch ) THEN
          tidx1 = tidx1a

       ! -------------------------------------------------------------
       ! If search shall be continued, adjust preferred year, then
       ! month, then day to the closest available year (month, day)
       ! in the time slices, and check if this is a better match.
       ! -------------------------------------------------------------
       ELSE

          ! Adjust year, month, and day (in this order).
          CNT  = 0
          DO
             CNT = CNT + 1
             IF ( ExitSearch .OR. CNT > 3 ) EXIT

             ! Adjust prefYMDhm at the given level (1=Y, 2=M, 3=D)
             CALL prefYMDhm_Adjust ( nTime, availYMDhm, prefYMDhm, CNT, tidx1a )

             ! verbose mode
             IF ( verb ) THEN
                write(MSG,'(A30,f14.0)') 'adjusted preferred datetime: ', &
                     prefYMDhm
                CALL HCO_MSG(HcoState%Config%Err,MSG)
             ENDIF

             ! check for time stamp with updated date/time
             CALL Check_AvailYMDhm ( Lct, nTime, availYMDhm, prefYMDhm, tidx1a )

             ! Can we leave now?
             ExitSearch = IsClosest( prefYMDhm, availYMDhm, nTime, tidx1a )
             IF ( ExitSearch ) tidx1 = tidx1a

          ENDDO
       ENDIF

       ! -------------------------------------------------------------
       ! If tidx1 still isn't defined, i.e. prefYMDhm is still
       ! outside the range of availYMDhm, set tidx1 to the closest
       ! available date. This must be 1 or nTime!
       ! -------------------------------------------------------------
       IF ( .NOT. ExitSearch ) THEN
          IF ( prefYMDhm < availYMDhm(1) ) THEN
             tidx1 = 1
          ELSE
             tidx1 = nTime
          ENDIF
       ENDIF

       ! -------------------------------------------------------------
       ! If we are dealing with 3-hourly or hourly data, select all timesteps
       ! -------------------------------------------------------------

       ! Hour flag is -1: wildcard
       IF ( Lct%Dct%Dta%ncHrs(1) == -1 .AND. nTime == 8 ) THEN
          tidx1 = 1
          tidx2 = nTime

          ! verbose mode
          IF ( verb ) THEN
             WRITE(MSG,*) 'Data is 3-hourly. Entire day will be read.'
             CALL HCO_MSG(HcoState%Config%Err,MSG)
          ENDIF
       ENDIF
       IF ( Lct%Dct%Dta%ncHrs(1) == -1 .AND. nTime == 24 ) THEN
          tidx1 = 1
          tidx2 = nTime

          ! verbose mode
          IF ( verb ) THEN
             WRITE(MSG,*) 'Data is hourly. Entire day will be read.'
             CALL HCO_MSG(HcoState%Config%Err,MSG)
          ENDIF
       ENDIF

       ! -------------------------------------------------------------
       ! If we are dealing with weekday data, pick the slice to be
       ! used based on the current day of week.
       ! The ncDys flag has been set in subroutine HCO_ExtractTime
       ! (hco_tidx_mod.F90) based upon the time attributes set in the
       ! configuration file. It can have the following values:
       ! >0  : specific days are given.
       ! -1  : wildcard (autodetect)
       ! -10 : WD (weekday).
       ! -999: determine from current simulation day.
       ! For specific days or if determined from the current datetime
       ! (flags >0 or -999), the weekday is not taken into account.
       ! If auto-detection is enabled, days are treated as weekday if
       ! (and only if) there are exactly 7 time slices. Otherwise, they
       ! are interpreted as 'regular' day data.
       ! If flag is set to -10, e.g. time attribute is 'WD', the current
       ! time index is assumed to hold Sunday data, with the following
       ! six slices being Mon, Tue, ..., Sat. For weekdaily data, all
       ! seven time slices will be read into memory so that at any given
       ! time, the local weekday can be taken (weekdaily data is always
       ! assumed to be in local time).
       ! -------------------------------------------------------------

       ! Day flag is -1: wildcard
       IF ( Lct%Dct%Dta%ncDys(1) == -1 .AND. nTime == 7 ) THEN
          tidx1 = 1
          tidx2 = nTime

          ! Make sure data is treated in local time
          Lct%Dct%Dta%IsLocTime = .TRUE.

       ! Day flag is -10: WD
       ELSEIF ( Lct%Dct%Dta%ncDys(1) == -10 ) THEN

          ! There must be at least 7 time slices
          IF ( nTime < 7 ) THEN
             MSG = 'Data must have exactly 7 time slices '// &
                   'if you set day attribute to WD: '//TRIM(Lct%Dct%cName)
             CALL HCO_ERROR( MSG, RC )
             IF ( ASSOCIATED(availYMDhm) ) THEN
                DEALLOCATE(availYMDhm)
                availYMDhm => NULL()
             ENDIF
             RETURN
          ENDIF

          ! If there are exactly seven time slices, interpret them as
          ! the seven weekdays.
          IF ( nTime == 7 ) THEN
             tidx1 = 1
             tidx2 = 7

          ! If there are more than 7 time slices, interpret the current
          ! selected index as sunday of the current time frame (e.g. sunday
          ! data of current month), and select the time slice index
          ! accordingly. This requires that there are at least 6 more time
          ! slices following the current one.
          ELSE
             IF ( tidx1 < 0 ) THEN
                WRITE(MSG,*) 'Cannot get weekday slices for: ', &
                   TRIM(Lct%Dct%cName), '. Cannot find first time slice.'
                CALL HCO_ERROR( MSG, RC )
                IF ( ASSOCIATED(availYMDhm) ) THEN
                   DEALLOCATE(availYMDhm)
                   availYMDhm => NULL()
                ENDIF
                RETURN
             ENDIF

             IF ( (tidx1+6) > nTime ) THEN
                WRITE(MSG,*) 'Cannot get weekday for: ',TRIM(Lct%Dct%cName), &
                   '. There are less than 6 additional time slices after ',  &
                   'selected start date ', availYMDhm(tidx1)
                CALL HCO_ERROR( MSG, RC )
                IF ( ASSOCIATED(availYMDhm) ) THEN
                   DEALLOCATE(availYMDhm)
                   availYMDhm => NULL()
                ENDIF
                RETURN
             ENDIF
             tidx2 = tidx1 + 6
          ENDIF

          ! Make sure data is treated in local time
          Lct%Dct%Dta%IsLocTime = .TRUE.

       ENDIF

       ! -------------------------------------------------------------
       ! Now need to set upper time slice index tidx2. This index
       ! is only different from tidx1 if:
       ! (1) We interpolate between two time slices, i.e. TimeCycle
       !     attribute is set to 'I'. In this case, we simply pick
       !     the next higher time slice index and calculate the
       !     weights for time1 and time2 based on the current time.
       ! (2) Multiple hourly slices are read (--> prefHr = -1 or -10,
       !     e.g. hour attribute in config. file was set to wildcard
       !     character or data is in local hours). In this case,
       !     check if there are multiple time slices for the selected
       !     date (y/m/d).
       ! tidx2 has already been set to proper value above if it's
       ! weekday data.
       ! -------------------------------------------------------------
       IF ( tidx2 < 0 ) THEN

          ! Interpolate between dates
          IF ( Lct%Dct%Dta%CycleFlag == HCO_CFLAG_INTER ) THEN

             CALL GetIndex2Interp( HcoState,   Lct,       nTime,      &
                                   availYMDhm, prefYMDhm, origYMDhm,  &
                                   tidx1,      tidx2,     wgt1,       &
                                   wgt2,       RC                   )
             IF ( RC /= HCO_SUCCESS ) THEN
                MSG = 'Error encountered in GetIndex2Interp for: '        // &
                     TRIM(Lct%Dct%Cname)
                CALL HCO_ERROR( MSG, RC )
                IF ( ASSOCIATED(availYMDhm) ) THEN
                   DEALLOCATE(availYMDhm)
                   availYMDhm => NULL()
                ENDIF
                RETURN
             ENDIF

          ! Check for multiple hourly data
          ELSEIF ( tidx1 > 0 .AND. prefHr < 0 ) THEN
             CALL SET_TIDX2 ( nTime, availYMDhm, tidx1, tidx2 )

             ! Denote as local time if necessary
             IF ( Lct%Dct%Dta%ncHrs(1) == -10 ) THEN
                Lct%Dct%Dta%IsLocTime = .TRUE.
             ENDIF
          ELSE
             tidx2 = tidx1
          ENDIF
       ENDIF

    ! ================================================================
    ! Case 3: No time slice available. Set both indeces to zero. Data
    ! with no time stamp must have CycleFlag 'Cycling'.
    ! ================================================================
    ELSE
       IF ( Lct%Dct%Dta%CycleFlag /= HCO_CFLAG_CYCLE ) THEN
          MSG = 'Field has no time/date variable - cycle flag must' // &
                'be set to `C` in the HEMCO configuration file:'    // &
                TRIM(Lct%Dct%cName)
          CALL HCO_ERROR( MSG, RC )
          IF ( ASSOCIATED(availYMDhm) ) THEN
             DEALLOCATE(availYMDhm)
             availYMDhm => NULL()
          ENDIF
          RETURN
       ENDIF

       tidx1 = 0
       tidx2 = 0
    ENDIF

    !-----------------------------------------------------------------
    ! Sanity check: if CycleFlag is set to 'Exact', the file time stamp
    ! must exactly match the current time.
    !-----------------------------------------------------------------
    IF ( (Lct%Dct%Dta%CycleFlag == HCO_CFLAG_EXACT) .AND. (tidx1 > 0) ) THEN
       IF ( availYMDhm(tidx1) /= prefYMDhm ) THEN
          tidx1 = -1
          tidx2 = -1
       ENDIF
    ENDIF

    !-----------------------------------------------------------------
    ! If multiple time slices are read, extract time interval between
    ! time slices in memory (in hours). This is to make sure that the
    ! cycling between the slices will be done at the correct rate
    ! (e.g. every hour, every 3 hours, ...).
    !-----------------------------------------------------------------
    IF ( (tidx2>tidx1) .AND. (Lct%Dct%Dta%CycleFlag/=HCO_CFLAG_INTER) ) THEN
       Lct%Dct%Dta%DeltaT = YMDhm2hrs( availYMDhm(tidx1+1) - availYMDhm(tidx1) )
    ELSE
       Lct%Dct%Dta%DeltaT = 0
    ENDIF

    ! verbose mode
    IF ( verb ) THEN
       WRITE(MSG,'(A30,I14)') 'selected tidx1: ', tidx1
       CALL HCO_MSG(HcoState%Config%Err,MSG)
       IF ( tidx1 > 0 ) THEN
          WRITE(MSG,'(A30,f14.0)') 'corresponding datetime 1: ', &
               availYMDhm(tidx1)
          CALL HCO_MSG(HcoState%Config%Err,MSG)
          IF ( wgt1 >= 0.0_sp ) THEN
             WRITE(MSG,*) 'weight1: ', wgt1
             CALL HCO_MSG(HcoState%Config%Err,MSG)
          ENDIF
       ENDIF

       IF ( (tidx2 /= tidx1) ) THEN
          WRITE(MSG,'(A30,I14)') 'selected tidx2: ', tidx2
          CALL HCO_MSG(HcoState%Config%Err,MSG)
          WRITE(MSG,'(A30,f14.0)') 'corresponding datetime 2: ', &
               availYMDhm(tidx2)
          CALL HCO_MSG(HcoState%Config%Err,MSG)
          IF ( wgt1 >= 0.0_sp ) THEN
             WRITE(MSG,*) 'weight2: ', wgt2
             CALL HCO_MSG(HcoState%Config%Err,MSG)
          ENDIF
       ENDIF

       WRITE(MSG,'(A30,I14)') 'assigned delta t [h]: ', Lct%Dct%Dta%DeltaT
       CALL HCO_MSG(HcoState%Config%Err,MSG)
       WRITE(MSG,*) 'local time? ', Lct%Dct%Dta%IsLocTime
       CALL HCO_MSG(HcoState%Config%Err,MSG)
    ENDIF

    ! ----------------------------------------------------------------
    ! TODO: set time brackets
    ! --> In future, we may want to set time brackets denoting the
    ! previous and next time slice available in the netCDF file. This
    ! may become useful for temporal interpolations and more efficient
    ! data update calls (only update if new time slice is available).
    ! ----------------------------------------------------------------

    !-----------------------------------------------------------------
    ! Prepare output, cleanup and leave
    !-----------------------------------------------------------------

    ! ncYr and ncMt are the year and month fo the time slice to be
    ! used. These values may be required to convert units to 'per
    ! seconds'.
    IF ( tidx1 > 0 ) THEN
       YMDhm  = availYMDhm(tidx1)
       YMDhm1 = availYMDhm(1)
       oYMDhm = origYMDhm
    ENDIF

    ! Deallocate and nullify the pointer
    IF ( ASSOCIATED(availYMDhm) ) THEN
       DEALLOCATE(availYMDhm)
       availYMDhm => NULL()
    ENDIF

    ! Return w/ success
    CALL HCO_LEAVE ( HcoState%Config%Err,  RC )

  END SUBROUTINE GET_TIMEIDX
!EOC
!------------------------------------------------------------------------------
!                   Harmonized Emissions Component (HEMCO)                    !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Check_AvailYMDhm
!
! !DESCRIPTION: Checks if prefYMDhm is within the range of availYMDhm
! and returns the location of the closest vector element that is in
! the past (--> tidx1). tidx1 is set to -1 otherwise.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Check_AvailYMDhm( Lct, N, availYMDhm, prefYMDhm, tidx1 )
!
! !INPUT PARAMETERS:
!
    TYPE(ListCont),   POINTER      :: Lct
    INTEGER,          INTENT(IN)   :: N
    REAL(dp),         INTENT(IN)   :: availYMDhm(N)
    REAL(dp),         INTENT(IN)   :: prefYMDhm
!
! !OUTPUT PARAMETERS:
!
    INTEGER,          INTENT(OUT)  :: tidx1
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
    INTEGER :: I, nTime

    !=================================================================
    ! Check_availYMDhm begins here
    !=================================================================

    ! Init
    tidx1 = -1

    ! Return if preferred datetime not within the vector range
    IF ( prefYMDhm < availYMDhm(1) .OR. prefYMDhm > availYMDhm(N) ) RETURN

    ! To avoid out-of-bounds error in the loop below:
    ! (1) For interpolated data, the upper loop limit should be N;
    ! (2) Otherwise, the upper loop limit should be N-1.
    ! (bmy, 4/28/21)
    nTime = N - 1
    IF ( Lct%Dct%Dta%CycleFlag == HCO_CFLAG_INTER ) nTime = N

    ! Get closest index that is not in the future
    DO I = 1, nTime

       ! NOTE: Epsilon test is more robust than an equality test
       ! for double-precision variables (bmy, 4/11/17)
       IF ( ABS( availYMDhm(I) - prefYMDhm ) < EPSILON ) THEN
          tidx1 = I
          EXIT
       ENDIF

       ! Check if next time slice is in the future, in which case the
       ! current slice is selected. Don't do this for a CycleFlag of
       ! 3 (==> exact match).
       IF ( (availYMDhm(I+1)       >  prefYMDhm        ) .AND. &
            (Lct%Dct%Dta%CycleFlag /= HCO_CFLAG_EXACT) ) THEN
          tidx1 = I
          EXIT
       ENDIF
    ENDDO

  END SUBROUTINE Check_AvailYMDhm
!EOC
!------------------------------------------------------------------------------
!                   Harmonized Emissions Component (HEMCO)                    !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: prefYMDhm_Adjust
!
! !DESCRIPTION: Adjusts prefYMDhm to the closest available time attribute. Can
! be adjusted for year (level=1), month (level=2), or day (level=3).
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE prefYMDhm_Adjust( N, availYMDhm, prefYMDhm, level, tidx1 )
!
! !INPUT PARAMETERS:
!
    INTEGER   , INTENT(IN)     :: N
    REAL(dp)  , INTENT(IN)     :: availYMDhm(N)
    INTEGER   , INTENT(IN)     :: level
    INTEGER   , INTENT(IN)     :: tidx1
!
! !INPUT/OUTPUT PARAMETERS:
!
    REAL(dp)  , INTENT(INOUT)  :: prefYMDhm
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
    ! Scalars
    INTEGER          :: I, IMIN, IMAX
    REAL(dp)         :: origYr,  origMt,  origDy, origHr, origMi
    REAL(dp)         :: refAttr, tmpAttr, newAttr
    REAL(dp)         :: iDiff,   minDiff
    REAL(dp)         :: modVal
    REAL(dp)         :: div

    !=================================================================
    ! prefYMDhm_Adjust begins here!
    !=================================================================

    ! Get original Yr, Mt, Day, Hr, Mi
    ! Time values are now in YYYYMMDDhhmm format (bmy, 4/11/17)
    origYr = FLOOR( MOD( prefYMDhm, 1.0e12_dp ) / 1.0e8_dp )
    origMt = FLOOR( MOD( prefYMDhm, 1.0e8_dp  ) / 1.0e6_dp )
    origDy = FLOOR( MOD( prefYMDhm, 1.0e6_dp  ) / 1.0e4_dp )
    origHr = FLOOR( MOD( prefYMDhm, 1.0e4_dp  ) / 1.0e2_dp )
    origMi = FLOOR( MOD( prefYMDhm, 1.0e2_dp  )            )

    ! Extract new attribute from availYMDhm and insert into prefYMDhm. Pick
    ! closest available value.
    SELECT CASE ( level )
       ! --- Year
       CASE ( 1 )
          modVal  = 1.0e12_dp
          div     = 1.0e8_dp
          refAttr = origYr

       ! --- Month
       CASE ( 2 )
          modVal  = 1.0e8_dp
          div     = 1.0e6_dp
          refAttr = origMt

       ! --- Day
       CASE ( 3 )
          modVal  = 1.0e6_dp
          div     = 1.0e4_dp
          refAttr = origMt

       ! --- Hour
       CASE ( 4 )
          modval  = 1.0e4_dp
          div     = 1.0e2_dp
          refAttr = origHr

       ! --- Minute
       CASE ( 5 )
          modVal  = 1.0e2_dp
          div     = 1.0_dp
          refAttr = origMi

       CASE DEFAULT
          RETURN
    END SELECT

    ! Maximum loop number:
    ! If tidx1 is already set, only search values in the past.
    IF ( tidx1 > 0 ) THEN
       IMIN = 1
       IMAX = tidx1

    ! If tidx1 is not yet set, prefYMDhm must be outside the range of
    ! availYMDhm. Pick only the closest available time stamp.
    ELSE
       IF ( prefYMDhm > availYMDhm(1) ) THEN
          IMIN = N
          IMAX = N
       ELSE
          IMIN = 1
          IMAX = 1
       ENDIF
    ENDIF

    ! Select current minimum value
    minDiff = 10000000000000000.0_dp
    newAttr = -1d0
    DO I = IMIN, IMAX
       tmpAttr = FLOOR( MOD(availYMDhm(I),modVal) / div )
       iDiff   = ABS( tmpAttr - refAttr )
       IF ( iDiff < minDiff ) THEN
          newAttr = tmpAttr
          minDiff = iDiff
       ENDIF
    ENDDO

    ! Just reuse current value if no better value could be found
    IF ( newAttr < 0 ) THEN
       newAttr = refAttr
    ENDIF

    ! Update variable
    ! --- Year
    IF ( level == 1 ) THEN
       prefYMDhm = ( newAttr * 1.0e8_dp ) + &
                   ( origMt  * 1.0e6_dp ) + &
                   ( origDy  * 1.0e4_dp ) + &
                   ( origHr  * 1.0e2_dp ) + &
                   ( origMi          )

    ! --- Month
    ELSEIF ( level == 2 ) THEN
       prefYMDhm = ( origYr  * 1.0e8_dp ) + &
                   ( newAttr * 1.0e6_dp ) + &
                   ( origDy  * 1.0e4_dp ) + &
                   ( origHr  * 1.0e2_dp ) + &
                   ( origMi             )

    ! --- Day
    ELSEIF ( level == 3 ) THEN
       prefYMDhm = ( origYr  * 1.0e8_dp  ) + &
                   ( origMt  * 1.0e6_dp  ) + &
                   ( newAttr * 1.0e4_dp  ) + &
                   ( origHr  * 1.0e2_dp  ) + &
                   ( origMi              )

    ! --- Hour
    ELSEIF ( level == 4 ) THEN
       prefYMDhm = ( origYr  * 1.0e8_dp  ) + &
                   ( origMt  * 1.0e6_dp  ) + &
                   ( origDy  * 1.0e4_dp  ) + &
                   ( newAttr * 1.0e2_dp  ) + &
                   ( origMi              )
    ! --- Minute
    ELSEIF ( level == 5 ) THEN
       prefYMDhm = ( origYr  * 1.0e8_dp  ) + &
                   ( origMt  * 1.0e6_dp  ) + &
                   ( origDy  * 1.0e4_dp  ) + &
                   ( origHr  * 1.0e2_dp  ) + &
                   ( newAttr             )

    ENDIF

  END SUBROUTINE prefYMDhm_Adjust
!EOC
!------------------------------------------------------------------------------
!                   Harmonized Emissions Component (HEMCO)                    !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Set_tIdx2
!
! !DESCRIPTION: sets the upper time slice index by selecting the range
! of all elements in availYMDhm with the same date (year,month,day) as
! availYMDh(tidx1).
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Set_tIdx2( N, availYMDhm, tidx1, tidx2 )
!
! !INPUT PARAMETERS:
!
    INTEGER,  INTENT(IN)  :: N               ! Number of times
    REAL(dp), INTENT(IN)  :: availYMDhm(N)   ! Time stamp vector
    INTEGER,  INTENT(IN)  :: tidx1           ! Lower time slice index
!
! !INPUT/OUTPUT PARAMETERS:
!
    INTEGER,  INTENT(OUT) :: tidx2           ! Upper time slice index
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
    INTEGER :: YMD, I, IYMD

    !=================================================================
    ! SET_TIDX2 begins here!
    !=================================================================

    ! Init
    tidx2 = tidx1

    ! Sanity check
    IF ( tidx1 == N ) RETURN

    ! Get wanted YMD
    YMD = floor(availYMDhm(tidx1) / 1.0e4_dp)

    ! See how many more tile slices with the same YMD exist from index
    ! tidx1 onwards.
    DO I = tidx1, N
       iYMD = floor(availYMDhm(I) / 1.0e4_dp)
       IF ( iYMD == YMD ) THEN
          tidx2 = I
       ELSEIF ( iYMD > YMD ) THEN
          EXIT
       ENDIF
    ENDDO

  END SUBROUTINE Set_tIdx2
!EOC
!------------------------------------------------------------------------------
!                   Harmonized Emissions Component (HEMCO)                    !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: IsClosest
!
! !DESCRIPTION: function IsClosest returns true if the selected time index
! is the 'closest' one. It is defined as being closest if:
! (a) the currently selected index exactly matches the preferred one.
! (b) the time gap between the preferred time stamp and the currently selected
! index is at least as small as any other gap of consecutive prior time stamps.
!\\
!\\
! !INTERFACE:
!
  FUNCTION IsClosest ( prefYMDhm, availYMDhm, nTime, ctidx1 ) RESULT ( Closest )
!
! !INPUT PARAMETERS:
!
    REAL(dp),   INTENT(IN)  :: prefYMDhm
    REAL(dp),   INTENT(IN)  :: availYMDhm(nTime)
    INTEGER,    INTENT(IN)  :: nTime
    INTEGER,    INTENT(IN)  :: ctidx1
!
! !OUTPUT PARAMETERS:
!
    LOGICAL              :: Closest
!
! !REVISION HISTORY:
!  03 Mar 2015 - C. Keller   - Initial version
!  See https://github.com/geoschem/hemco for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER :: N
    INTEGER :: diff, idiff

    !=================================================================
    ! IsClosest begins here!
    !=================================================================

    ! Init
    Closest = .TRUE.

    ! It's not closest if index is not defined
    IF ( ctidx1 <= 0 ) THEN
       Closest = .FALSE.
       RETURN
    ENDIF

    ! It's closest if it is the first index
    IF ( ctidx1 == 1 ) RETURN

    ! It's closest if it matches date exactly
    ! NOTE: Epsilon test is more robust than an equality test
    ! for double-precision variables (bmy, 4/11/17)
    IF ( ABS( availYMDhm(ctidx1) - prefYMDhm ) < EPSILON ) RETURN

    ! It's closest if current select one is in the future
    IF ( availYMDhm(ctidx1) > prefYMDhm ) RETURN

    ! Check if any of the time stamps in the past have closer intervals
    ! than the current select time stamp to it's previous one
    diff = prefYMDhm - availYMDhm(ctidx1)
    DO N = 2, ctidx1
       idiff = availYMDhm(N) - availYMDhm(N-1)
       IF ( idiff < diff ) THEN
          Closest = .FALSE.
          RETURN
       ENDIF
    ENDDO

  END FUNCTION IsClosest
!EOC
!------------------------------------------------------------------------------
!                   Harmonized Emissions Component (HEMCO)                    !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: GetIndex2Interp
!
! !DESCRIPTION: GetIndex2Interp
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE GetIndex2Interp ( HcoState,  Lct,                   &
                               nTime,     availYMDhm,            &
                               prefYMDhm, origYMDhm,  tidx1,     &
                               tidx2,     wgt1,       wgt2,  RC   )
!
! !INPUT PARAMETERS:
!
    TYPE(HCO_State),  POINTER       :: HcoState
    TYPE(ListCont),   POINTER       :: Lct
    INTEGER,          INTENT(IN)    :: nTime
    REAL(dp),         INTENT(IN)    :: availYMDhm(nTime)
    REAL(dp),         INTENT(IN)    :: prefYMDhm
    REAL(dp),         INTENT(IN)    :: origYMDhm
    INTEGER,          INTENT(IN)    :: tidx1
!
! !OUTPUT PARAMETERS:
!
    INTEGER,          INTENT(OUT)   :: tidx2
!
! !INPUT/OUTPUT PARAMETERS:
!
    REAL(sp),         INTENT(INOUT) :: wgt1
    REAL(sp),         INTENT(INOUT) :: wgt2
    INTEGER,          INTENT(INOUT) :: RC
!
! !REVISION HISTORY:
!  02 Mar 2015 - C. Keller   - Initial version
!  See https://github.com/geoschem/hemco for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Scalars
    INTEGER             :: I
    REAL(dp)            :: tmpYMDhm
    LOGICAL             :: verb

    ! Strings
    CHARACTER(LEN=255)  :: MSG
    CHARACTER(LEN=255)  :: LOC = 'GetIndex2Interp (hcoio_util_mod.F90)'

    !=================================================================
    ! GetIndex2Interp begins here
    !=================================================================

    ! Verbose mode?
    verb = HCO_IsVerb(HcoState%Config%Err,3)

    ! If the originally wanted datetime was below the available data
    ! range, set all weights to the first index.
    IF ( origYMDhm <= availYMDhm(1) ) THEN
       tidx2 = tidx1
       wgt1  = 1.0_sp
       wgt2  = 0.0_sp

    ! If the originally wanted datetime is beyond the available data
    ! range, set tidx2 to tidx1 but leave weights in their original
    ! values (-1.0). The reason is that we will attempt to interpolate
    ! between a second file, which is only done if the weights are
    ! negative.
    ELSEIF ( origYMDhm >= availYMDhm(nTime) ) THEN
       tidx2 = tidx1

    ! No interpolation needed if there is a time slices that exactly
    ! matches the (originally) preferred datetime.
    ! NOTE: An Epsilon test is more robust than an equality test
    ! for double-precision variables (bmy, 4/11/17)
    ELSEIF ( ABS( origYMDhm - availYMDhm(tidx1) ) < EPSILON ) THEN
       tidx2 = tidx1
       wgt1  = 1.0_sp
       wgt2  = 0.0_sp

    ! If we are inside the data range but none of the time slices
    ! matches the preferred datetime, get the second time slices that
    ! shall be used for data interpolation. This not necessarily needs
    ! to be the consecutive time slice. For instance, imagine a data
    ! set that contains montlhly data for years 2005 and 2010. For
    ! Feb 2007, we would want to interpolate between Feb 2005 and Feb
    ! 2010 data. The index tidx1 already points to Feb 2005, but the
    ! upper index tidx2 needs to be set accordingly.
    ELSE

       ! Init
       tidx2 = -1

       ! Search for a time slice in the future that has the same
       ! month/day/hour as currently selected time slice.
       tmpYMDhm = availYMDhm(tidx1)
       DO
          ! Increase by one year
          tmpYMDhm = tmpYMDhm + 1.0e8_dp

          ! Exit if we are beyond available dates
          IF ( tmpYMDhm > availYMDhm(nTime) ) EXIT

          ! Check if there is a time slice with that date
          DO I = tidx1,nTime
             IF ( tmpYMDhm == availYMDhm(I) ) THEN
                tidx2 = I
                EXIT
             ENDIF
          ENDDO
          IF ( tidx2 > 0 ) EXIT
       ENDDO

       ! Repeat above but now only modify month.
       IF ( tidx2 < 0 ) THEN
          tmpYMDhm = availYMDhm(tidx1)
          DO
             ! Increase by one month
             tmpYMDhm = tmpYMDhm + 1.0e6_dp

             ! Exit if we are beyond available dates
             IF ( tmpYMDhm > availYMDhm(nTime) ) EXIT

             ! Check if there is a time slice with that date
             DO I = tidx1,nTime
                IF ( ABS( tmpYMDhm - availYMDhm(I) ) < EPSILON ) THEN
                   tidx2 = I
                   EXIT
                ENDIF
             ENDDO
             IF ( tidx2 > 0 ) EXIT
          ENDDO
       ENDIF

       ! Repeat above but now only modify day
       IF ( tidx2 < 0 ) THEN
          tmpYMDhm = availYMDhm(tidx1)
          DO
             ! Increase by one day
             tmpYMDhm = tmpYMDhm + 1.0e4_dp

             ! Exit if we are beyond available dates
             IF ( tmpYMDhm > availYMDhm(nTime) ) EXIT

             ! Check if there is a time slice with that date
             DO I = tidx1,nTime
                IF ( tmpYMDhm == availYMDhm(I) ) THEN
                   tidx2 = I
                   EXIT
                ENDIF
             ENDDO
             IF ( tidx2 > 0 ) EXIT
          ENDDO
       ENDIF

       ! If all of those tests failed, simply get the next time
       ! slice.
       IF ( tidx2 < 0 ) THEN
          tidx2 = tidx1 + 1

          ! Make sure that tidx2 does not exceed nTime, which is
          ! the number of time slices in the file. This can cause
          ! an out-of-bounds error. (bmy, 3/7/19)
          IF ( tidx2 > nTime ) tidx2 = nTime

          ! Prompt warning
          WRITE(MSG,*) 'Having problems in finding the next time slice ', &
                'to interpolate from, just take the next available ',     &
                'slice. Interpolation will be performed from ',           &
                availYMDhm(tidx1), ' to ', availYMDhm(tidx2), '. Data ',    &
                'container: ', TRIM(Lct%Dct%cName)
          CALL HCO_WARNING(HcoState%Config%Err, MSG, RC, WARNLEV=1, THISLOC=LOC)
       ENDIF

       ! Calculate weights wgt1 and wgt2 to be given to slice 1 and
       ! slice2, respectively.
       CALL GetWeights ( availYMDhm(tidx1), availYMDhm(tidx2), origYMDhm, &
                         wgt1, wgt2 )

    ENDIF

    ! Return w/ success
    RC = HCO_SUCCESS

  END SUBROUTINE GetIndex2Interp
!EOC
!------------------------------------------------------------------------------
!                   Harmonized Emissions Component (HEMCO)                    !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: GetWeights
!
! !DESCRIPTION: Helper function to get the interpolation weights between
! two datetime intervals (int1, int2) and for a given time cur.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE GetWeights ( int1, int2, cur, wgt1, wgt2 )
!
! !INPUT PARAMETERS:
!
    REAL(dp),         INTENT(IN   )   :: int1, int2, cur
!
! !INPUT/OUTPUT PARAMETERS:
!
    REAL(sp),         INTENT(  OUT)   :: wgt1, wgt2
!
! !REVISION HISTORY:
!  04 Mar 2015 - C. Keller - Initial version
!  See https://github.com/geoschem/hemco for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    REAL(dp)              :: diff1, diff2
    REAL(dp)              :: jdc, jd1, jd2

    !=================================================================
    ! GetWeights begins here!
    !=================================================================

    ! Convert dates to Julian dates
    jdc = YMDhm2jd ( cur  )
    jd1 = YMDhm2jd ( int1 )
    jd2 = YMDhm2jd ( int2 )

    ! Check if outside of range
    IF ( jdc <= jd1 ) THEN
       wgt1 = 1.0_sp
    ELSEIF ( jdc >= jd2 ) THEN
       wgt1 = 0.0_sp
    ELSE
       diff1 = jd2 - jdc
       diff2 = jd2 - jd1
       wgt1  = diff1 / diff2
    ENDIF

    ! second weight is just complement of wgt1
    wgt2  = 1.0_sp - wgt1

  END SUBROUTINE GetWeights
!EOC
!------------------------------------------------------------------------------
!                   Harmonized Emissions Component (HEMCO)                    !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: YMDhm2jd
!
! !DESCRIPTION: returns the julian date of element YMDhm.
!\\
!\\
! !INTERFACE:
!
  FUNCTION YMDhm2jd ( YMDhm ) RESULT ( jd )
!
! !USES:
!
    USE HCO_Julday_Mod
!
! !INPUT PARAMETERS:
!
    REAL(dp), INTENT(IN)  :: YMDhm
!
! !INPUT/OUTPUT PARAMETERS:
!
    REAL(hp) :: jd
!
! !REVISION HISTORY:
!  24 Feb 2019 - C. Keller - Initial version
!  See https://github.com/geoschem/hemco for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER               :: yr, mt, dy, hr, mn
    REAL(dp)              :: utc, day

    !=================================================================
    ! YMDh2jd begins here!
    !=================================================================
    yr  = FLOOR( MOD( YMDhm, 1.0e12_dp ) / 1.0e8_dp )
    mt  = FLOOR( MOD( YMDhm, 1.0e8_dp  ) / 1.0e6_dp )
    dy  = FLOOR( MOD( YMDhm, 1.0e6_dp  ) / 1.0e4_dp )
    hr  = FLOOR( MOD( YMDhm, 1.0e4_dp  ) / 1.0e2_dp )
    mn  = FLOOR( MOD( YMDhm, 1.0e2_dp  ) )
    utc = ( REAL(hr,dp) / 24.0_dp    ) + &
          ( REAL(mn,dp) / 1440.0_dp  ) + &
          ( REAL(0 ,dp) / 86400.0_dp )
    day = REAL(dy,dp) + utc
    jd  = JULDAY( yr, mt, day )

  END FUNCTION YMDhm2jd
!EOC
!------------------------------------------------------------------------------
!                   Harmonized Emissions Component (HEMCO)                    !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: YMDhm2hrs
!
! !DESCRIPTION: returns the hours of element YMDhm. For simplicity, 30 days are
! assigned to every month. At the moment, this routine is only called to
! determine the time interval between two emission time slices (DeltaT) and
! this approximation is good enough.
!\\
!\\
! !INTERFACE:
!
  FUNCTION YMDhm2hrs ( YMDhm ) RESULT ( hrs )
!
! !INPUT PARAMETERS:
!
    REAL(dp), INTENT(IN)  :: YMDhm
!
! !INPUT/OUTPUT PARAMETERS:
!
    INTEGER              :: hrs
!
! !REVISION HISTORY:
!  26 Jan 2015 - C. Keller - Initial version
!  See https://github.com/geoschem/hemco for complete history
!EOP
!------------------------------------------------------------------------------
!BOC

    !=================================================================
    ! YMDh2hrs begins here!
    !=================================================================
    hrs = FLOOR( MOD( YMDhm, 1.0e12_dp ) / 1.0e8_dp ) * 8760 + &
          FLOOR( MOD( YMDhm, 1.0e8_dp  ) / 1.0e6_dp ) * 720  + &
          FLOOR( MOD( YMDhm, 1.0e6_dp  ) / 1.0e4_dp ) * 24   + &
          FLOOR( MOD( YMDhm, 1.0e4_dp  ) / 1.0e2_dp )

  END FUNCTION YMDhm2hrs
!EOC
!------------------------------------------------------------------------------
!                   Harmonized Emissions Component (HEMCO)                    !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Normalize_Area
!
! !DESCRIPTION: Subroutine Normalize\_Area normalizes the given array
! by the surface area calculated from the given netCDF file.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Normalize_Area( HcoState, Array, nlon, LatEdge, FN, RC )
!
! !INPUT PARAMETERS:
!
    TYPE(HCO_State),  POINTER         :: HcoState    ! HEMCO state object
    INTEGER,          INTENT(IN   )   :: nlon        ! # of lon midpoints
    REAL(hp),         POINTER         :: LatEdge(:)  ! lat edges
    CHARACTER(LEN=*), INTENT(IN   )   :: FN          ! filename
!
! !INPUT/OUTPUT PARAMETERS:
!
    REAL(sp),         POINTER         :: Array(:,:,:,:) ! Data
    INTEGER,          INTENT(INOUT)   :: RC             ! Return code
!
! !REVISION HISTORY:
!  13 Mar 2013 - C. Keller - Initial version
!  See https://github.com/geoschem/hemco for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    REAL(hp)              :: DLAT, AREA
    INTEGER               :: NLAT, J
    CHARACTER(LEN=255)    :: MSG, LOC

    !=================================================================
    ! NORNALIZE_AREA begins here!
    !=================================================================

    ! Initialize
    LOC    = 'NORMALIZE_AREA (hcoio_util_mod.F90 )'

    ! Check array size
    NLAT = SIZE(LatEdge,1) - 1

    IF ( SIZE(Array,1) /= nlon ) THEN
       MSG = 'Array size does not agree with nlon: ' // TRIM(FN)
       CALL HCO_ERROR( MSG, RC, THISLOC=LOC )
       RETURN
    ENDIF
    IF ( SIZE(Array,2) /= NLAT ) THEN
       MSG = 'Array size does not agree with nlat: ' // TRIM(FN)
       CALL HCO_ERROR( MSG, RC, THISLOC=LOC )
       RETURN
    ENDIF

    ! Loop over all latitudes
    DO J = 1, NLAT
       ! get grid box area in m2 for grid box with lower and upper latitude
       ! llat/ulat:  Area = 2 * PI * Re^2 * DLAT / nlon,
       ! where DLAT = abs( sin(ulat) - sin(llat) )
       DLAT = ABS( SIN(LatEdge(J+1)*HcoState%Phys%PI_180)  &
                   - SIN(LatEdge(J)*HcoState%Phys%PI_180) )
       AREA = ( 2_hp * HcoState%Phys%PI * DLAT * HcoState%Phys%Re**2 ) &
              / REAL(nlon,hp)

       ! convert array data to m-2
       ARRAY(:,J,:,:) = ARRAY(:,J,:,:) / AREA
    ENDDO

    ! Prompt a warning
    WRITE(MSG,*) 'No area unit found in ' // TRIM(FN) // ' - convert to m-2!'
    CALL HCO_WARNING ( HcoState%Config%Err, MSG, RC, WARNLEV=1, THISLOC=LOC )

    ! Leave w/ success
    RC = HCO_SUCCESS

  END SUBROUTINE Normalize_Area
!EOC
!------------------------------------------------------------------------------
!                   Harmonized Emissions Component (HEMCO)                    !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: SrcFile_Parse
!
! !DESCRIPTION: Routine SrcFile\_Parse parses the source file name ('ncFile')
! of the provided list container Lct. In particular, it searches for tokens
! such as $ROOT, $YYYY, etc., within the file name and replaces those values
! with the intendend characters. The parsed file name is returned in string
! srcFile, while the original file name is retained in Lct.
!\\
!\\
! It now also checks if the file exists. If the file does not exist and the
! file name contains date tokens, it tries to adjust the file name to the
! closest available date in the past. The optional flag FUTURE can be used
! to denote that the next available file in the future shall be selected,
! even if there is a file that exactly matches the preferred date time. This
! is useful for interpolation between fields.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE SrcFile_Parse ( HcoState, Lct, srcFile, FOUND, RC, &
                             Direction, Year )
!
! !USES:
!
    USE HCO_TIDX_MOD,         ONLY : HCO_GetPrefTimeAttr
    USE HCO_TIDX_MOD,         ONLY : tIDx_IsInRange
    USE HCO_CLOCK_MOD,        ONLY : HcoClock_Get
    USE HCO_CLOCK_MOD,        ONLY : Get_LastDayOfMonth
!
! !INPUT PARAMETERS:
!
    TYPE(HCO_State),  POINTER                 :: HcoState   ! HEMCO state object
    TYPE(ListCont),   POINTER                 :: Lct        ! HEMCO list
    INTEGER,          INTENT(IN   ), OPTIONAL :: Direction  ! Look for file in
                                                            ! future (+1) or
                                                            ! past (-1)
    INTEGER,          INTENT(IN   ), OPTIONAL :: Year       ! To use fixed year
!
! !OUTPUT PARAMETERS:
!
    CHARACTER(LEN=*), INTENT(  OUT)           :: srcFile    ! output string
    LOGICAL,          INTENT(  OUT)           :: FOUND      ! Does file exist?
!
! !INPUT/OUTPUT PARAMETERS:
!
    INTEGER,          INTENT(INOUT)           :: RC         ! return code
!
! !REVISION HISTORY:
!  01 Oct 2014 - C. Keller - Initial version
!  See https://github.com/geoschem/hemco for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER :: INC,     CNT,    TYPCNT, TYP,   NEWTYP
    INTEGER :: prefYr,  prefMt, prefDy, prefHr, prefMn
    INTEGER :: origYr,  origMt, origDy, origHr
    LOGICAL :: hasFile, hasYr,  hasMt,  hasDy, hasHr
    LOGICAL :: nextTyp
    CHARACTER(LEN=1023) :: MSG, LOC
    CHARACTER(LEN=1023) :: srcFileOrig

    ! maximum # of iterations for file search
    INTEGER, PARAMETER :: MAXIT = 10000

    !=================================================================
    ! SrcFile_Parse
    !=================================================================

    ! Initialize
    LOC     = 'SrcFile_Parse (HCOIO_UTIL_MOD.F90)'
    RC      = HCO_SUCCESS
    found   = .FALSE.
    srcFile = Lct%Dct%Dta%ncFile

    ! verbose mode
    IF ( HCO_IsVerb(HcoState%Config%Err,3) ) THEN
       WRITE(MSG,*) 'Parsing source file and replacing tokens'
       CALL HCO_MSG(HcoState%Config%Err,MSG)
    ENDIF

    ! Get preferred dates (to be passed to parser)
    CALL HCO_GetPrefTimeAttr ( HcoState, Lct, &
                               prefYr, prefMt, prefDy, prefHr, prefMn, RC )
    IF ( RC /= HCO_SUCCESS ) THEN
        CALL HCO_ERROR( 'ERROR 1', RC, THISLOC=LOC )
        RETURN
    ENDIF

    ! Make sure dates are not negative
    IF ( prefYr <= 0 ) THEN
       CALL HcoClock_Get( HcoState%Clock, cYYYY = prefYr, RC = RC )
       IF ( RC /= HCO_SUCCESS ) THEN
           CALL HCO_ERROR( 'ERROR 2', RC, THISLOC=LOC )
           RETURN
       ENDIF
    ENDIF
    IF ( prefMt <= 0 ) THEN
       CALL HcoClock_Get( HcoState%Clock, cMM   = prefMt, RC = RC )
       IF ( RC /= HCO_SUCCESS ) THEN
           CALL HCO_ERROR( 'ERROR 3', RC, THISLOC=LOC )
           RETURN
       ENDIF
    ENDIF
    IF ( prefDy <= 0 ) THEN
       CALL HcoClock_Get( HcoState%Clock, cDD   = prefDy, RC = RC )
       IF ( RC /= HCO_SUCCESS ) THEN
           CALL HCO_ERROR( 'ERROR 4', RC, THISLOC=LOC )
           RETURN
       ENDIF
    ENDIF
    IF ( prefHr <  0 ) THEN
       CALL HcoClock_Get( HcoState%Clock, cH    = prefHr, RC = RC )
       IF ( RC /= HCO_SUCCESS ) THEN
           CALL HCO_ERROR( 'ERROR 5', RC, THISLOC=LOC )
           RETURN
       ENDIF
    ENDIF

    ! Eventually replace default preferred year with specified one
    IF ( PRESENT(Year) ) prefYr = Year

    ! Call the parser
    CALL HCO_CharParse ( HcoState%Config, srcFile, prefYr, prefMt, &
                         prefDy, prefHr, prefMn, RC )
    IF ( RC /= HCO_SUCCESS ) THEN
        CALL HCO_ERROR( 'ERROR 6', RC, THISLOC=LOC )
        RETURN
    ENDIF
    srcFileOrig = TRIM(srcFile)

    ! Check if file exists
    INQUIRE( FILE=TRIM(srcFile), EXIST=HasFile )

    ! If the direction flag is on, force HasFile to be false.
    IF ( PRESENT(Direction) ) THEN
       IF ( Direction /= 0 ) HasFile = .FALSE.
    ENDIF

    !-----------------------------------------------------------------------
    ! If this is a HEMCO dry-run simulation, then do not enter the loop
    ! where we will attempt to go back in time until a file is found.
    ! For the dry-run we need to report all files, even missing.
    ! This fixes Github issue geoschem/geos-chem #312. (bmy, 6/9/20)
    !-----------------------------------------------------------------------
    IF ( HcoState%Options%isDryRun ) THEN

       ! Make sure that the year is not 1, this indicates that the
       ! preferred year is outside of the years specified in the
       ! time range settings in the configuration file, and will
       ! lead to files with a year of "0001" in the path.
       ! (bmy, 6/9/20)
       IF ( prefyr == 1 ) THEN
          MSG = 'Cannot find file for current simulation time: ' // &
               TRIM(srcFile) // ' - Cannot get field ' // &
               TRIM(Lct%Dct%cName) // '. Please check file name ' // &
               'and time (incl. time range flag) in the config. file'
          CALL HCO_ERROR( MSG, RC )
          RETURN
       ENDIF

       ! Otherwise return with success
       RC    = HCO_SUCCESS
       Found = HasFile
       RETURN
    ENDIF

    ! If file does not exist, check if we can adjust prefYr, prefMt, etc.
    IF ( .NOT. HasFile .AND. Lct%Dct%DctType /= HCO_CFLAG_EXACT ) THEN

       ! Check if any token exist
       HasYr = ( INDEX(TRIM(Lct%Dct%Dta%ncFile),'YYYY') > 0 )
       HasMt = ( INDEX(TRIM(Lct%Dct%Dta%ncFile),'MM'  ) > 0 )
       HasDy = ( INDEX(TRIM(Lct%Dct%Dta%ncFile),'DD'  ) > 0 )
       HasHr = ( INDEX(TRIM(Lct%Dct%Dta%ncFile),'HH'  ) > 0 )

       ! Search for file
       IF ( HasYr .OR. HasMt .OR. HasDy .OR. HasHr ) THEN

          ! Date increments
          INC = -1
          IF ( PRESENT(Direction) ) THEN
             INC = Direction
          ENDIF

          ! Initialize counters
          CNT = 0

          ! Type is the update type (see below)
          TYP = 0

          ! Mirror preferred variables
          origYr = prefYr
          origMt = prefMt
          origDy = prefDy
          origHr = prefHr

          ! Do until file is found or counter exceeds threshold
          DO WHILE ( .NOT. HasFile )

             ! Inrease counter
             CNT = CNT + 1
             IF ( CNT > MAXIT ) EXIT

             ! Increase update type if needed:
             nextTyp = .FALSE.

             ! Type 0: Initialization
             IF ( TYP == 0 ) THEN
                nextTyp = .TRUE.
             ! Type 1: update hour only
             ELSEIF ( TYP == 1 .AND. TYPCNT > 24 ) THEN
                nextTyp = .TRUE.
             ! Type 2: update day only
             ELSEIF ( TYP == 2 .AND. TYPCNT > 31 ) THEN
                nextTyp = .TRUE.
             ! Type 3: update month only
             ELSEIF ( TYP == 3 .AND. TYPCNT > 12 ) THEN
                nextTyp = .TRUE.
             ! Type 4: update year only
             ELSEIF ( TYP == 4 .AND. TYPCNT > 300 ) THEN
                nextTyp = .TRUE.
             ! Type 5: update hour and day
             ELSEIF ( TYP == 5 .AND. TYPCNT > 744 ) THEN
                nextTyp = .TRUE.
             ! Type 6: update day and month
             ELSEIF ( TYP == 6 .AND. TYPCNT > 372 ) THEN
                nextTyp = .TRUE.
             ! Type 7: update month and year
             ELSEIF ( TYP == 7 .AND. TYPCNT > 3600 ) THEN
                EXIT
             ENDIF

             ! Get next type
             IF ( nextTyp ) THEN
                NEWTYP = -1
                IF     ( hasHr .AND. TYP < 1 ) THEN
                   NEWTYP = 1
                ELSEIF ( hasDy .AND. TYP < 2 ) THEN
                   NEWTYP = 2
                ELSEIF ( hasMt .AND. TYP < 3 ) THEN
                   NEWTYP = 3
                ELSEIF ( hasYr .AND. TYP < 4 ) THEN
                   NEWTYP = 4
                ELSEIF ( hasDy .AND. TYP < 2 ) THEN
                   NEWTYP = 5
                ELSEIF ( hasDy .AND. TYP < 2 ) THEN
                   NEWTYP = 6
                ELSEIF ( hasDy .AND. TYP < 2 ) THEN
                   NEWTYP = 7
                ENDIF

                ! Exit if no other type found
                IF ( NEWTYP < 0 ) EXIT

                ! This is the new type, reset type counter
                TYP    = NEWTYP
                TYPCNT = 0

                ! Make sure we reset all values
                prefYr = origYr
                prefMt = origMt
                prefDy = origDy
                prefHr = origHr

             ENDIF

             ! Update preferred datetimes
             SELECT CASE ( TYP )
                ! Adjust hour only
                CASE ( 1 )
                   prefHr = prefHr + INC
                ! Adjust day only
                CASE ( 2 )
                   prefDy = prefDy + INC
                ! Adjust month only
                CASE ( 3 )
                   prefMt = prefMt + INC
                ! Adjust year only
                CASE ( 4 )
                   prefYr = prefYr + INC
                ! Adjust hour and day
                CASE ( 5 )
                   prefHr = prefHr + INC
                   IF ( MOD(TYPCNT,24) == 0 ) prefDy = prefDy + INC
                ! Adjust day and month
                CASE ( 6 )
                   prefDy = prefDy + INC
                   IF ( MOD(TYPCNT,31) == 0 ) prefMt = prefMt + INC
                ! Adjust month and year
                CASE ( 7 )
                   prefMt = prefMt + INC
                   IF ( MOD(TYPCNT,12) == 0 ) prefYr = prefYr + INC
                CASE DEFAULT
                   EXIT
             END SELECT

             ! Check if we need to adjust a year/month/day/hour
             IF ( prefHr < 0 ) THEN
                prefHr = 23
                prefDy = prefDy - 1
             ENDIF
             IF ( prefHr > 23 ) THEN
                prefHr = 0
                prefDy = prefDy + 1
             ENDIF
             IF ( prefDy < 1  ) THEN
                prefDy = 31
                prefMt = prefMt - 1
             ENDIF
             IF ( prefDy > 31 ) THEN
                prefDy = 1
                prefMt = prefMt + 1
             ENDIF
             IF ( prefMt < 1  ) THEN
                prefMt = 12
                prefYr = prefYr - 1
             ENDIF
             IF ( prefMt > 12 ) THEN
                prefMt = 1
                prefYr = prefYr + 1
             ENDIF

             ! Make sure day does not exceed max. number of days in this month
             prefDy = MIN( prefDy, Get_LastDayOfMonth( prefMt, prefYr ) )

             ! Mirror original file
             srcFile = Lct%Dct%Dta%ncFile

             ! Call the parser with adjusted values
             CALL HCO_CharParse ( HcoState%Config, srcFile, prefYr, &
                                  prefMt, prefDy, prefHr, prefMn, RC )
             IF ( RC /= HCO_SUCCESS ) THEN
                 CALL HCO_ERROR( 'ERROR 7', RC, THISLOC=LOC )
                 RETURN
             ENDIF

             ! Check if this file exists
             INQUIRE( FILE=TRIM(srcFile), EXIST=HasFile )

             ! Update counter
             TYPCNT = TYPCNT + 1
          ENDDO
       ENDIF
    ENDIF

    ! Additional check for data with a given range: make sure that the selected
    ! field is not outside of the given range
    IF ( HasFile .AND. ( Lct%Dct%Dta%CycleFlag == HCO_CFLAG_RANGE ) ) THEN
       HasFile = TIDX_IsInRange ( Lct, prefYr, prefMt, prefDy, prefHr )
    ENDIF

    ! Restore original source file name and date to avoid confusion in log file
    IF ( .not. HasFile ) THEN
       srcFile = Trim(srcFileOrig)
    ENDIF

    ! Return variable
    FOUND = HasFile

    ! Return w/ success
    RC = HCO_SUCCESS

  END SUBROUTINE SrcFile_Parse
!EOC
!------------------------------------------------------------------------------
!                   Harmonized Emissions Component (HEMCO)                    !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: SigmaMidToEdges
!
! !DESCRIPTION: Helper routine to interpolate sigma mid point values to edges.
! A simple linear interpolation is performed.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE SigmaMidToEdges ( HcoState, SigMid, SigEdge, RC )
!
! !INPUT PARAMETERS:
!
    TYPE(HCO_State),  POINTER                 :: HcoState        ! HEMCO state
    REAL(hp),         POINTER                 :: SigMid(:,:,:)   ! sigma levels
!
! !OUTPUT PARAMETERS:
!
    REAL(hp),         POINTER                 :: SigEdge(:,:,:)  ! sigma edges
    INTEGER,          INTENT(  OUT)           :: RC              ! return code
!
! !REVISION HISTORY:
!  03 Oct 2013 - C. Keller - Initial version
!  See https://github.com/geoschem/hemco for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER            :: L, AS
    INTEGER            :: nx, ny, nz
    CHARACTER(LEN=255) :: MSG
    CHARACTER(LEN=255) :: LOC = 'SigmaMidToEdges (hcoio_util_mod.F90)'

    !=================================================================
    ! SigmaMidToEdges begins here!
    !=================================================================

    ! Allocate space as required
    nx = SIZE(SigMid,1)
    ny = SIZE(SigMid,2)
    nz = SIZE(SigMid,3)
    IF ( ASSOCIATED(SigEdge) ) DEALLOCATE(SigEdge)
    ALLOCATE(SigEdge(nx,ny,nz+1),STAT=AS)
    IF ( AS/=0 ) THEN
       CALL HCO_ERROR( 'Allocate SigEdge', RC, &
                       THISLOC=LOC )
       RETURN
    ENDIF
    SigEdge = 0.0_hp

    ! Calculate sigma edges by linear interpolation (symmetric mid-points)
    DO L = 1, nz-1
       SigEdge(:,:,L+1) = ( SigMid(:,:,L) + SigMid(:,:,L+1) ) / 2.0_hp
    ENDDO

    ! Get outermost values:
    SigEdge(:,:,1   ) = SigMid(:,:,1 ) - ( SigEdge(:,:,2) - SigMid(:,:,1)   )
    SigEdge(:,:,nz+1) = SigMid(:,:,nz) + ( SigMid(:,:,nz) - SigEdge(:,:,nz) )

    ! Return w/ success
    RC = HCO_SUCCESS

  END SUBROUTINE SigmaMidToEdges
!EOC
!------------------------------------------------------------------------------
!                   Harmonized Emissions Component (HEMCO)                    !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: CheckMissVal
!
! !DESCRIPTION: Checks for missing values in the passed array. Missing values
! of base emissions and masks are set to 0, missing values of scale factors
! are set to 1.
!\\
! !INTERFACE:
!
  SUBROUTINE CheckMissVal ( Lct, Arr )
!
! !INPUT PARAMETERS:
!
    TYPE(ListCont),   POINTER                 :: Lct
    REAL(sp),         POINTER                 :: Arr(:,:,:,:)
!
! !REVISION HISTORY:
!  04 Mar 2015 - C. Keller - Initial version
!  See https://github.com/geoschem/hemco for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    !=================================================================
    ! CheckMissVal begins here!
    !=================================================================

    ! Error trap
    IF ( .NOT. ASSOCIATED(Arr) ) RETURN

    IF ( ANY(Arr == HCO_MISSVAL) ) THEN
       ! Base emissions
       IF ( Lct%Dct%DctType == HCO_DCTTYPE_BASE ) THEN
          WHERE(Arr == HCO_MISSVAL) Arr = 0.0_sp
       ! Scale factor
       ELSEIF ( Lct%Dct%DctType == HCO_DCTTYPE_SCAL ) THEN
          WHERE(Arr == HCO_MISSVAL) Arr = 1.0_sp
       ! Mask
       ELSEIF ( Lct%Dct%DctType == HCO_DCTTYPE_MASK ) THEN
          WHERE(Arr == HCO_MISSVAL) Arr = 0.0_sp
       ENDIF
    ENDIF

  END SUBROUTINE CheckMissVal
!EOC
!------------------------------------------------------------------------------
!                   Harmonized Emissions Component (HEMCO)                    !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: GetArbDimIndex
!
! !DESCRIPTION: Subroutine GetArbDimIndex returns the index of the arbitrary
! file dimension. -1 if no such dimension is defined.
!\\
! !INTERFACE:
!
  SUBROUTINE GetArbDimIndex( HcoState, Lun, Lct, ArbIdx, RC )
!
! !USES:
!
    USE HCO_m_netcdf_io_checks
    USE HCO_m_netcdf_io_get_dimlen
    USE HCO_ExtList_Mod,    ONLY : GetExtOpt
!
! !INPUT PARAMETERS:
!
    TYPE(HCO_State),  POINTER                 :: HcoState
    INTEGER,          INTENT(IN   )           :: Lun
    TYPE(ListCont),   POINTER                 :: Lct
!
! !OUTPUT PARAMETERS:
!
    INTEGER,          INTENT(  OUT)           :: ArbIdx
    INTEGER,          INTENT(  OUT)           :: RC
!
! !REVISION HISTORY:
!  22 Sep 2015 - C. Keller - Initial version
!  See https://github.com/geoschem/hemco for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER             :: TargetVal, nVal
    LOGICAL             :: Found
    CHARACTER(LEN=255)  :: ArbDimVal
    CHARACTER(LEN=511)  :: MSG
    CHARACTER(LEN=255)  :: LOC = 'GetArbDimIndex (hcoio_util_mod.F90)'

    !=================================================================
    ! GetArbDimIndex
    !=================================================================

    ! Assume success until otherwise
    RC = HCO_SUCCESS

    ! Init
    ArbIdx = -1
    IF ( TRIM(Lct%Dct%Dta%ArbDimName) == 'none' ) RETURN

    ! Check if variable exists
    Found = Ncdoes_Dim_Exist ( Lun, TRIM(Lct%Dct%Dta%ArbDimName) )
    IF ( .NOT. Found ) THEN
       MSG = 'Cannot read dimension ' // TRIM(Lct%Dct%Dta%ArbDimName) &
             // ' from file ' // &
             TRIM(Lct%Dct%Dta%ncFile)
       CALL HCO_ERROR( MSG, RC, THISLOC=LOC )
       RETURN
    ENDIF

    ! Get dimension length
    CALL Ncget_Dimlen ( Lun, TRIM(Lct%Dct%Dta%ArbDimName), nVal )

    ! Get value to look for. This is archived in variable ArbDimVal.
    ! Eventually need to extract value from HEMCO settings
    ArbDimVal = TRIM(Lct%Dct%Dta%ArbDimVal)

    ! If string starts with a number, evaluate value directly
    IF ( ArbDimVal(1:1) == '0' .OR. &
         ArbDimVal(1:1) == '1' .OR. &
         ArbDimVal(1:1) == '2' .OR. &
         ArbDimVal(1:1) == '3' .OR. &
         ArbDimVal(1:1) == '4' .OR. &
         ArbDimVal(1:1) == '5' .OR. &
         ArbDimVal(1:1) == '6' .OR. &
         ArbDimVal(1:1) == '7' .OR. &
         ArbDimVal(1:1) == '8' .OR. &
         ArbDimVal(1:1) == '9'       ) THEN
       READ(ArbDimVal,*) TargetVal

    ! Otherwise, assume this is a HEMCO option (including a token)
    ELSE
       IF ( ArbDimVal(1:1) == '$' ) ArbDimVal = ArbDimVal(2:LEN(ArbDimVal))
       CALL GetExtOpt ( HcoState%Config, ExtNr=-999, &
                        OptName=TRIM(ArbDimVal), &
                        OptValInt=TargetVal, FOUND=Found, RC=RC )
       IF ( RC /= HCO_SUCCESS ) THEN
           CALL HCO_ERROR( 'ERROR 8', RC, THISLOC=LOC )
           RETURN
       ENDIF
       IF ( .NOT. Found ) THEN
          WRITE(MSG,*) 'Cannot evaluate additional dimension value ', &
             TRIM(ArbDimVal), '. This does not seem to be a number nor ', &
             'a HEMCO token/setting. This error happened when evaluating ', &
             'dimension ', TRIM(Lct%Dct%Dta%ArbDimName), ' belonging to ', &
             'file ', TRIM(Lct%Dct%Dta%ncFile)
          CALL HCO_ERROR( MSG, RC, THISLOC=LOC )
          RETURN
       ENDIF
    ENDIF

    IF ( TargetVal > nVal ) THEN
       WRITE(MSG,*) 'Desired dimension value ', TargetVal, &
          ' exceeds corresponding dimension length on that file: ', nVal, &
          'This error happened when evaluating ', &
          'dimension ', TRIM(Lct%Dct%Dta%ArbDimName), ' belonging to ', &
          'file ', TRIM(Lct%Dct%Dta%ncFile)
       CALL HCO_ERROR( MSG, RC, THISLOC=LOC )
       RETURN

    ELSE
       ArbIdx = TargetVal
    ENDIF

    ! Verbose
    IF ( HcoState%amIRoot .AND. HCO_IsVerb( HcoState%Config%Err, 2 ) ) THEN
       WRITE(MSG,*) 'Additional dimension ', TRIM(Lct%Dct%Dta%ArbDimName), &
                    ' in ', TRIM(Lct%Dct%Dta%ncFile), ': use index ',      &
                    ArbIdx, ' (set: ', Lct%Dct%Dta%ArbDimVal, ')'
       CALL HCO_MSG(HcoState%Config%Err,MSG)
    ENDIF

    ! Return w/ success
    RC = HCO_SUCCESS

  END SUBROUTINE GetArbDimIndex
!EOC
#endif
!------------------------------------------------------------------------------
!                   Harmonized Emissions Component (HEMCO)                    !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: HCOIO_ReadOther
!
! !DESCRIPTION: Subroutine HCOIO\_ReadOther is a wrapper routine to
! read data from sources other than netCDF.
!\\
!\\
! If a file name is given (ending with '.txt'), the data are assumed
! to hold country-specific values (e.g. diurnal scale factors). In all
! other cases, the data is directly read from the configuration file
! (scalars).
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE HCOIO_ReadOther( HcoState, Lct, RC )
!
! !USES:
!
!
! !INPUT PARAMTERS:
!
    TYPE(HCO_State), POINTER          :: HcoState    ! HEMCO state
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(ListCont),   POINTER         :: Lct
    INTEGER,          INTENT(INOUT)   :: RC
!
! !REVISION HISTORY:
!  22 Dec 2014 - C. Keller: Initial version
!  See https://github.com/geoschem/hemco for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    CHARACTER(LEN=255) :: MSG, LOC

    !======================================================================
    ! HCOIO_ReadOther begins here
    !======================================================================
    LOC = 'HCOIO_ReadOther (HCOIO_UTIL_MOD.F90)'

    ! Error check: data must be in local time
    IF ( .NOT. Lct%Dct%Dta%IsLocTime ) THEN
       MSG = 'Cannot read data from file that is not in local time: ' // &
             TRIM(Lct%Dct%cName)
       CALL HCO_ERROR( MSG, RC, THISLOC='HCOIO_ReadOther (hcoio_dataread_mod.F90)' )
       RETURN
    ENDIF

    ! Read an ASCII file as country values
    IF ( INDEX( TRIM(Lct%Dct%Dta%ncFile), '.txt' ) > 0 ) THEN
       CALL HCOIO_ReadCountryValues( HcoState, Lct, RC )
       IF ( RC /= HCO_SUCCESS ) THEN
           CALL HCO_ERROR( 'ERROR 9', RC, THISLOC=LOC )
           RETURN
       ENDIF

    ! Directly read from configuration file otherwise
    ELSE
       CALL HCOIO_ReadFromConfig( HcoState, Lct, RC )
       IF ( RC /= HCO_SUCCESS ) THEN
           CALL HCO_ERROR( 'ERROR 10', RC, THISLOC=LOC )
           RETURN
       ENDIF
    ENDIF

    ! Return w/ success
    RC = HCO_SUCCESS

  END SUBROUTINE HCOIO_ReadOther
!EOC
!------------------------------------------------------------------------------
!                   Harmonized Emissions Component (HEMCO)                    !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: HCOIO_ReadCountryValues
!
! !DESCRIPTION: Subroutine HCOIO\_ReadCountryValues
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE HCOIO_ReadCountryValues ( HcoState, Lct, RC )
!
! !USES:
!
    USE HCO_inquireMod,     ONLY : findFreeLUN
    USE HCO_CHARTOOLS_MOD,  ONLY : HCO_CMT, HCO_SPC, NextCharPos
    USE HCO_EmisList_Mod,   ONLY : HCO_GetPtr
    USE HCO_FileData_Mod,   ONLY : FileData_ArrCheck
!
! !INPUT PARAMTERS:
!
    TYPE(HCO_State), POINTER          :: HcoState    ! HEMCO state
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(ListCont),   POINTER         :: Lct
    INTEGER,          INTENT(INOUT)   :: RC
!
! !REVISION HISTORY:
!  22 Dec 2014 - C. Keller: Initial version
!  See https://github.com/geoschem/hemco for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER               :: IUFILE, IOS
    INTEGER               :: ID1, ID2, I, NT, CID, NLINE
    REAL(sp), POINTER     :: CNTR(:,:)
    INTEGER,  ALLOCATABLE :: CIDS(:,:)
    REAL(hp), POINTER     :: Vals(:)
    LOGICAL               :: Verb
    CHARACTER(LEN=2047)   :: LINE
    CHARACTER(LEN=255)    :: MSG, DUM, CNT
    CHARACTER(LEN=255)    :: LOC = 'HCOIO_ReadCountryValues (hcoio_util_mod.F90)'

    !======================================================================
    ! HCOIO_ReadCountryValues begins here
    !======================================================================

    ! Init
    CNTR => NULL()
    Vals => NULL()

    ! verbose mode?
    Verb = HCO_IsVerb(HcoState%Config%Err,2)

    ! Verbose
    IF ( Verb ) THEN
       MSG = 'Use country-specific values for ' // TRIM(Lct%Dct%cName)
       CALL HCO_MSG(HcoState%Config%Err,MSG)
       MSG = '- Source file: ' // TRIM(Lct%Dct%Dta%ncFile)
       CALL HCO_MSG(HcoState%Config%Err,MSG)
    ENDIF

    ! Open file
    IUFILE = FindFreeLun()
    OPEN ( IUFILE, FILE=TRIM( Lct%Dct%Dta%ncFile ), STATUS='OLD', IOSTAT=IOS )
    IF ( IOS /= 0 ) THEN
       MSG = 'Cannot open ' // TRIM(Lct%Dct%Dta%ncFile)
       CALL HCO_ERROR( MSG, RC, THISLOC=LOC )
       RETURN
    ENDIF

    ! Repeat for every line
    NLINE = 0
    DO

       ! Read line
       READ( IUFILE, '(a)', IOSTAT=IOS ) LINE

       ! End of file?
       IF ( IOS < 0 ) EXIT

       ! Error?
       IF ( IOS > 0 ) THEN
          MSG = 'Error reading ' // TRIM(Lct%Dct%Dta%ncFile)
          MSG = TRIM(MSG) // ' - last valid line: ' // TRIM(LINE)
          CALL HCO_ERROR( MSG, RC, THISLOC=LOC )
          RETURN
       ENDIF

       ! Skip commented lines and/or empty lines
       IF ( TRIM(LINE) == ''      ) CYCLE
       IF ( LINE(1:1)  == HCO_CMT ) CYCLE

       ! First (valid) line holds the name of the mask container
       IF ( NLINE == 0 ) THEN

          ! Get pointer to mask. Convert to integer
          CALL HCO_GetPtr( HcoState, TRIM(LINE), CNTR, RC )
          IF ( RC /= HCO_SUCCESS ) THEN
              CALL HCO_ERROR( 'ERROR 11', RC, THISLOC=LOC )
              RETURN
          ENDIF
          ALLOCATE( CIDS(HcoState%NX, HcoState%NY), STAT=IOS )
          IF ( IOS /= 0 ) THEN
             CALL HCO_ERROR( 'Cannot allocate CIDS', RC, THISLOC=LOC )
             RETURN
          ENDIF
          CIDS = NINT(CNTR)

          ! Verbose
          IF ( HCO_IsVerb(HcoState%Config%Err,3) ) THEN
             MSG = '- Use ID mask ' // TRIM(LINE)
             CALL HCO_MSG(HcoState%Config%Err,MSG)
          ENDIF

          ! Go to next line
          NLINE = NLINE + 1
          CYCLE
       ENDIF

       ! Get first space character to skip country name.
       ! We assume here that a country name is given right at the
       ! beginning of the line, e.g. 'USA 744 1.05/1.02/...'
       ID1 = NextCharPos( LINE, HCO_SPC )
       CNT = LINE(1:ID1)

       ! Get country ID
       DO I = ID1, LEN(LINE)
          IF ( LINE(I:I) /= HCO_SPC ) EXIT
       ENDDO
       ID1 = I
       ID2 = NextCharPos( LINE, HCO_SPC, START=ID1 )

       IF ( ID2 >= LEN(LINE) .OR. ID2 < 0 ) THEN
          MSG = 'Cannot extract country ID from: ' // TRIM(LINE)
          CALL HCO_ERROR( MSG, RC, THISLOC=LOC )
          RETURN
       ENDIF
       DUM = LINE(ID1:ID2)
       READ( DUM, * ) CID

       ! Extract data values
       ID1  = ID2+1
       ID2  = LEN(LINE)
       LINE = LINE(ID1:ID2)
       CALL GetDataVals( HcoState, Lct, LINE, Vals, RC )
       IF ( RC /= HCO_SUCCESS ) THEN
           CALL HCO_ERROR( 'ERROR 12', RC, THISLOC=LOC )
           RETURN
       ENDIF

       ! Check data / array dimensions
       NT = SIZE(Vals,1)
       CALL FileData_ArrCheck( HcoState%Config, Lct%Dct%Dta, &
                               HcoState%NX, HcoState%NY, NT, RC )
       IF ( RC /= HCO_SUCCESS ) THEN
           CALL HCO_ERROR( 'ERROR 13', RC, THISLOC=LOC )
           RETURN
       ENDIF

       ! Pass to data array. If the country ID is larger than zero, fill
       ! only those grid boxes. Otherwise, fill all grid boxes that have
       ! not yet been filled.
       DO I = 1, NT
          IF ( CID == 0 ) THEN
             WHERE ( Lct%Dct%Dta%V2(I)%Val <= 0.0_sp )
                Lct%Dct%Dta%V2(I)%Val = Vals(I)
             ENDWHERE
          ELSE
             WHERE ( CIDS == CID )
                Lct%Dct%Dta%V2(I)%Val = Vals(I)
             ENDWHERE
          ENDIF
       ENDDO

       ! Verbose
       IF ( HCO_IsVerb(HcoState%Config%Err,3) ) THEN
          WRITE(MSG,*) '- Obtained values for ',TRIM(CNT),' ==> ID:', CID
          CALL HCO_MSG(HcoState%Config%Err,MSG)
       ENDIF

       ! Cleanup
       IF ( ASSOCIATED(Vals) ) DEALLOCATE( Vals )
       Vals => NULL()

       ! Update # of read lines
       NLINE = NLINE + 1
    ENDDO

    ! Close file
    CLOSE ( IUFILE )

    ! Data is 2D
    Lct%Dct%Dta%SpaceDim  = 2

    ! Make sure data is in local time
    IF ( .NOT. Lct%Dct%Dta%IsLocTime ) THEN
       Lct%Dct%Dta%IsLocTime = .TRUE.
       MSG = 'Data assigned to mask regions will be treated in local time: '//&
              TRIM(Lct%Dct%cName)
       CALL HCO_WARNING( HcoState%Config%Err, MSG, RC, WARNLEV=2, THISLOC=LOC )
    ENDIF

    ! Cleanup
    Cntr => NULL()
    IF ( ALLOCATED(CIDS) ) DEALLOCATE ( CIDS )

    ! Return w/ success
    RC = HCO_SUCCESS

  END SUBROUTINE HCOIO_ReadCountryValues
!EOC
!------------------------------------------------------------------------------
!                   Harmonized Emissions Component (HEMCO)                    !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: HCOIO_ReadFromConfig
!
! !DESCRIPTION: Subroutine HCOIO\_ReadFromConfig reads data directly from
! the configuration file (instead of reading it from a netCDF file).
! These data is always assumed to be spatially uniform, but it is possible
! to specify multiple time slices by separating the individual time slice
! values by the HEMCO separator sign ('/' by default). The time dimension
! of these data is either determined from the srcTime attribute or estimated
! from the number of time slices provided. For example, if no srcTime is
! specified and 24 time slices are provided, data is assumed to represent
! hourly data. Similarly, data is assumed to represent weekdaily or monthly
! data for 7 or 12 time slices, respectively.
!\\
!\\
! If the srcTime attribute is defined, the time slices are determined from
! this attribute. Only one time dimension (year, month, day, or hour) can
! be defined for scalar fields!
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE HCOIO_ReadFromConfig( HcoState, Lct, RC )
!
! !USES:
!
    USE HCO_FILEDATA_MOD,   ONLY : FileData_ArrCheck
!
! !INPUT PARAMTERS:
!
    TYPE(HCO_State), POINTER          :: HcoState    ! HEMCO state
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(ListCont),   POINTER         :: Lct
    INTEGER,          INTENT(INOUT)   :: RC
!
! !REVISION HISTORY:
!  24 Jul 2014 - C. Keller: Initial version
!  See https://github.com/geoschem/hemco for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER            :: I, NT
    REAL(hp), POINTER  :: Vals(:)
    CHARACTER(LEN=255) :: MSG
    CHARACTER(LEN=255) :: LOC = 'HCOIO_ReadFromConfig (hcoio_util_mod.F90)'

    !======================================================================
    ! HCOIO_ReadFromConfig begins here
    !======================================================================

    ! Init
    Vals => NULL()

    ! Verbose
    IF ( HCO_IsVerb(HcoState%Config%Err,2) ) THEN
       WRITE(MSG, *) 'Read from config file: ', TRIM(Lct%Dct%cName)
       CALL HCO_MSG(HcoState%Config%Err,MSG)
    ENDIF

    !-------------------------------------------------------------------
    ! Get data values for this time step.
    !-------------------------------------------------------------------
    CALL GetDataVals( HcoState, Lct, Lct%Dct%Dta%ncFile, Vals, RC )
    IF ( RC /= HCO_SUCCESS ) THEN
        CALL HCO_ERROR( 'ERROR 14', RC, THISLOC=LOC )
        RETURN
    ENDIF

    !-------------------------------------------------------------------
    ! Copy data into array.
    !-------------------------------------------------------------------

    ! Number of values
    NT = SIZE(Vals,1)

    ! For masks, interpret data as mask corners (lon1/lat1/lon2/lat2)
    ! with no time dimension
    IF ( Lct%Dct%DctType == HCO_DCTTYPE_MASK ) THEN

       ! Make sure data is allocated
       CALL FileData_ArrCheck( HcoState%Config, Lct%Dct%Dta, &
                               HcoState%NX, HcoState%NY, 1, RC )
       IF ( RC /= HCO_SUCCESS ) THEN
           CALL HCO_ERROR( 'ERROR 15', RC, THISLOC=LOC )
           RETURN
       ENDIF

       ! Fill array: 1.0 within grid box, 0.0 outside.
       CALL FillMaskBox( HcoState, Lct, Vals, RC )
       IF ( RC /= HCO_SUCCESS ) THEN
           CALL HCO_ERROR( 'ERROR 16', RC, THISLOC=LOC )
           RETURN
       ENDIF

       ! Data is 2D
       Lct%Dct%Dta%SpaceDim = 2

    ! For base emissions and scale factors, interpret data as scalar
    ! values with a time dimension.
    ELSE

       CALL FileData_ArrCheck( HcoState%Config, Lct%Dct%Dta, 1, 1, NT, RC )
       IF ( RC /= HCO_SUCCESS ) THEN
           CALL HCO_ERROR( 'ERROR 17', RC, THISLOC=LOC )
           RETURN
       ENDIF
       DO I = 1, NT
          Lct%Dct%Dta%V2(I)%Val(1,1) = Vals(I)
!==============================================================================
! KLUDGE BY BOB YANTOSCA (05 Jan 2016)
!
! This WRITE statement avoids a seg fault in some Intel Fortran Compiler
! versions, such as ifort 12 and ifort 13.  The ADVANCE="no" prevents
! carriage returns from being added to the log file, and the '' character
! will prevent text from creeping across the screen.
!
! NOTE: This section only gets executed during the initialization phase,
! when we save data not read from netCDF files into the HEMCO data structure.
! This type of data includes scale factors and mask data specified as vectors
! in the HEMCO configuration file.  Therefore, this section will only get
! executed at startup, so the WRITE statment should not add significant
! overhead to the simulation.
!
! The root issue seems to be an optimization bug in the compiler.
!==============================================================================
#if defined( LINUX_IFORT )
          WRITE( 6, '(a)', ADVANCE='no' ) ''
#endif

       ENDDO

       ! Data is 1D
       Lct%Dct%Dta%SpaceDim  = 1

       ! Make sure data is in local time
       IF ( .NOT. Lct%Dct%Dta%IsLocTime ) THEN
          Lct%Dct%Dta%IsLocTime = .TRUE.
          MSG = 'Scale factors read from file are treated as local time: '// &
                 TRIM(Lct%Dct%cName)
          CALL HCO_WARNING( HcoState%Config%Err, MSG, RC, WARNLEV=2, &
                            THISLOC=LOC )
       ENDIF

    ENDIF

    ! Cleanup
    IF ( ASSOCIATED(Vals) ) DEALLOCATE(Vals)

    ! Return w/ success
    RC = HCO_SUCCESS

  END SUBROUTINE HCOIO_ReadFromConfig
!EOC
!------------------------------------------------------------------------------
!                   Harmonized Emissions Component (HEMCO)                    !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: GetSliceIdx
!
! !DESCRIPTION: gets the time slice index to be used for data directly
! read from the HEMCO configuration file. prefDt denotes the preferred
! time attribute (year, month, or day). DtType is used to identify the
! time attribute type (1=year, 2=month, 3=day). The time slice index will
! be selected based upon those two variables. IDX is the selected time
! slice index. It will be set to -1 if the current simulation date
! is outside of the specified time range and the time cycle attribute is
! not enabled for this field.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE GetSliceIdx ( HcoState, Lct, DtType, prefDt, IDX, RC )
!
! !INPUT PARAMETERS:
!
    TYPE(HCO_State),  POINTER                 :: HcoState
    TYPE(ListCont),   POINTER                 :: Lct
    INTEGER,          INTENT(IN   )           :: DtType
    INTEGER,          INTENT(IN   )           :: prefDt
!
! !INPUT/OUTPUT PARAMETERS:
!
    INTEGER,          INTENT(INOUT)           :: IDX
    INTEGER,          INTENT(INOUT)           :: RC
!
! !REVISION HISTORY:
!  13 Mar 2013 - C. Keller - Initial version
!  See https://github.com/geoschem/hemco for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER            :: lowDt, uppDt
    CHARACTER(LEN=255) :: MSG
    CHARACTER(LEN=255) :: LOC = 'GetSliceIdx (hcoio_util_mod.F90)'

    !=================================================================
    ! GetSliceIdx begins here!
    !=================================================================

    ! Init
    RC = HCO_SUCCESS

    ! Get upper and lower time range
    IF ( DtType == 1 ) THEN
       lowDt = Lct%Dct%Dta%ncYrs(1)
       uppDt = Lct%Dct%Dta%ncYrs(2)
    ELSEIF ( DtType == 2 ) THEN
       lowDt = Lct%Dct%Dta%ncMts(1)
       uppDt = Lct%Dct%Dta%ncMts(2)
    ELSEIF ( DtType == 3 ) THEN
       lowDt = Lct%Dct%Dta%ncDys(1)
       uppDt = Lct%Dct%Dta%ncDys(2)
    ELSE
       WRITE(MSG,*) "DtType must be one of 1, 2, 3: ", DtType
       CALL HCO_ERROR( MSG, RC, THISLOC=LOC )
       RETURN
    ENDIF

    ! Check for cycle flags:

    ! Data cycle set to range or exact date: in these cases, the
    ! the preferred date will be equal to the current date, so
    ! check if the preferred date is indeed within the available
    ! range (lowDt, uppDt).
    ! For data only to be used within the specified range, set
    ! index to -1. This will force the scale factors to be set to
    ! zero!
    IF ( prefDt < lowDt .OR. prefDt > uppDt ) THEN
       IF ( ( Lct%Dct%Dta%CycleFlag == HCO_CFLAG_EXACT ) .OR.      &
            ( Lct%Dct%Dta%CycleFlag == HCO_CFLAG_RANGE )     ) THEN
          IDX = -1
          RETURN
       ELSE
          ! this here should never happen, since for a cycle flag of 1,
          ! the preferred date should always be restricted to the range
          ! of available time stamps.
          MSG = 'preferred date is outside of range: ' // TRIM(Lct%Dct%cName)
          CALL HCO_ERROR( MSG, RC, THISLOC=LOC )
          RETURN
       ENDIF
    ENDIF

    ! If the code makes it to here, prefDt is within the available data range
    ! and we simply get the wanted index from the current index and the lowest
    ! available index.
    IDX = prefDt - lowDt + 1

  END SUBROUTINE GetSliceIdx
!EOC
!------------------------------------------------------------------------------
!                   Harmonized Emissions Component (HEMCO)                    !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: GetDataVals
!
! !DESCRIPTION: Subroutine GetDataVals extracts the data values from ValStr
! and writes them into vector Vals. ValStr is typically a character string
! read from an external ASCII file or directly from the HEMCO configuration
! file. Depending on the time specifications provided in the configuration
! file, Vals will be filled with only a subset of the values of ValStr.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE GetDataVals ( HcoState, Lct, ValStr, Vals, RC )
!
! !USES:
!
    USE HCO_CHARTOOLS_MOD,  ONLY : HCO_CharSplit
    USE HCO_EXTLIST_MOD,    ONLY : HCO_GetOpt
    USE HCO_UNIT_MOD,       ONLY : HCO_Unit_Change
    USE HCO_tIdx_Mod,       ONLY : HCO_GetPrefTimeAttr
    USE HCO_CLOCK_MOD,      ONLY : HcoClock_Get
!
! !INPUT PARAMTERS:
!
    TYPE(HCO_State),  POINTER         :: HcoState    ! HEMCO state
    CHARACTER(LEN=*), INTENT(IN   )   :: ValStr
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(ListCont),   POINTER         :: Lct
    INTEGER,          INTENT(INOUT)   :: RC
!
! !OUTPUT PARAMETERS:
!
    REAL(hp),         POINTER         :: Vals(:)
!
! !REVISION HISTORY:
!  22 Dec 2014 - C. Keller: Initial version
!  See https://github.com/geoschem/hemco for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER            :: HcoID
    INTEGER            :: I, N, NUSE, AS
    INTEGER            :: IDX1, IDX2
    INTEGER            :: AreaFlag, TimeFlag, Check
    INTEGER            :: prefYr, prefMt, prefDy, prefHr, prefMn
    INTEGER            :: cYr,    cMt,    cDy,    cHr
    REAL(hp)           :: MW_g
    REAL(hp)           :: UnitFactor
    REAL(hp)           :: FileVals(100)
    REAL(hp), POINTER  :: FileArr(:,:,:,:)
    LOGICAL            :: IsPerArea
    LOGICAL            :: IsMath
    CHARACTER(LEN=255) :: MSG
    CHARACTER(LEN=255) :: LOC = 'GetDataVals (hcoio_util_mod.F90)'

    !======================================================================
    ! GetDataVals begins here
    !======================================================================

    ! Initialize
    FileArr => NULL()

    ! Shadow species properties needed for unit conversion
    HcoID = Lct%Dct%HcoID
    IF ( HcoID > 0 ) THEN
       MW_g = HcoState%Spc(HcoID)%MW_g
    ELSE
       MW_g = -999.0_hp
    ENDIF

    ! Is this a math expression?
    IsMath = .FALSE.
    IF ( LEN(ValStr) > 5 ) THEN
       IF ( ValStr(1:5)=='MATH:' ) IsMath = .TRUE.
    ENDIF

    ! Evaluate math expression if string starts with 'MATH:'
    IF ( IsMath ) THEN
       CALL ReadMath ( HcoState, Lct, ValStr, FileVals, N, RC )
       IF ( RC /= HCO_SUCCESS ) THEN
           CALL HCO_ERROR( 'ERROR 18', RC, THISLOC=LOC )
           RETURN
       ENDIF

    ! Use regular string parser otherwise
    ELSE
       CALL HCO_CharSplit ( ValStr, &
                            HCO_GetOpt(HcoState%Config%ExtList,'Separator'), &
                            HCO_GetOpt(HcoState%Config%ExtList,'Wildcard'), &
                            FileVals, N, RC )
       IF ( RC /= HCO_SUCCESS ) THEN
           CALL HCO_ERROR( 'ERROR 19', RC, THISLOC=LOC )
           RETURN
       ENDIF
    ENDIF

    ! Return w/ error if no scale factor defined
    IF ( N == 0 ) THEN
       MSG = 'Cannot read data: ' // TRIM(Lct%Dct%cName) // &
             ': ' // TRIM(ValStr)
       CALL HCO_ERROR( MSG, RC, THISLOC=LOC)
       RETURN
    ENDIF

    ! Get the preferred times, i.e. the preferred year, month, day,
    ! or hour (as specified in the configuration file).
    CALL HCO_GetPrefTimeAttr( HcoState, Lct, &
                              prefYr, prefMt, prefDy, prefHr, prefMn, RC )
    IF ( RC /= HCO_SUCCESS ) THEN
        CALL HCO_ERROR( 'ERROR 20', RC, THISLOC=LOC )
        RETURN
    ENDIF

    ! ----------------------------------------------------------------
    ! For masks, assume that values represent the corners of the mask
    ! box, e.g. there must be four values. Masks are time-independent
    ! and unitless
    ! ----------------------------------------------------------------
    IF ( Lct%Dct%DctType == HCO_DCTTYPE_MASK ) THEN

       ! There must be exactly four values
       IF ( N /= 4 ) THEN
          MSG = 'Mask values are not lon1/lat1/lon2/lat2: ' // &
                TRIM(ValStr) // ' --> ' // TRIM(Lct%Dct%cName)
          CALL HCO_ERROR( MSG, RC, THISLOC=LOC )
          RETURN
       ENDIF

       ! Pass to FileArr array (will be used below)
       NUSE = 4
       ALLOCATE( FileArr(1,1,1,NUSE), STAT=AS )
       IF ( AS /= 0 ) THEN
          MSG = 'Cannot allocate FileArr'
          CALL HCO_ERROR( MSG, RC, THISLOC=LOC )
          RETURN
       ENDIF
       FileArr(1,1,1,:) = FileVals(1:NUSE)

    ! ----------------------------------------------------------------
    ! For non-masks, the data is interpreted as uniform values with
    ! a time dimension. Need to select the time slices to be used at
    ! this time (depending on the provided time attributes), as well
    ! as to ensure that values are in the correct units.
    ! Use all time slices unless a time interval is provided in
    ! attribute srcTime of the configuration file.
    ! ----------------------------------------------------------------
    ELSE

       ! If there is only one value use this one and ignore any time
       ! preferences.
       IF ( N == 1 ) THEN
          NUSE = 1
          IDX1 = 1
          IDX2 = 1

       ! If it's a math expression use all passed values
       ELSEIF ( IsMath ) THEN
          NUSE = N
          IDX1 = 1
          IDX2 = N

       ELSE
          ! Currently, data read directly from the configuration file can only
          ! represent one time dimension, i.e. it can only be yearly, monthly,
          ! daily (or hourly data, but this is read all at the same time).

          ! Annual data
          IF ( Lct%Dct%Dta%ncYrs(1) /= Lct%Dct%Dta%ncYrs(2) ) THEN
             ! Error check
             IF ( Lct%Dct%Dta%ncMts(1) /= Lct%Dct%Dta%ncMts(2) .OR. &
                  Lct%Dct%Dta%ncDys(1) /= Lct%Dct%Dta%ncDys(2) .OR. &
                  Lct%Dct%Dta%ncHrs(1) /= Lct%Dct%Dta%ncHrs(2)       ) THEN
                MSG = 'Data must not have more than one time dimension: ' // &
                       TRIM(Lct%Dct%cName)
                CALL HCO_ERROR( MSG, RC, THISLOC=LOC )
                RETURN
             ENDIF

             CALL GetSliceIdx ( HcoState, Lct, 1, prefYr, IDX1, RC )
             IF ( RC /= HCO_SUCCESS ) THEN
                 CALL HCO_ERROR( 'ERROR 21', RC, THISLOC=LOC )
                 RETURN
             ENDIF
             IDX2 = IDX1
             NUSE = 1

          ! Monthly data
          ELSEIF ( Lct%Dct%Dta%ncMts(1) /= Lct%Dct%Dta%ncMts(2) ) THEN
             ! Error check
             IF ( Lct%Dct%Dta%ncDys(1) /= Lct%Dct%Dta%ncDys(2) .OR. &
                  Lct%Dct%Dta%ncHrs(1) /= Lct%Dct%Dta%ncHrs(2)       ) THEN
                MSG = 'Data must only have one time dimension: ' // &
                      TRIM(Lct%Dct%cName)
                CALL HCO_ERROR( MSG, RC, THISLOC=LOC )
                RETURN
             ENDIF

             CALL GetSliceIdx ( HcoState, Lct, 2, prefMt, IDX1, RC )
             IF ( RC /= HCO_SUCCESS ) THEN
                 CALL HCO_ERROR( 'ERROR 22', RC, THISLOC=LOC )
                 RETURN
             ENDIF
             IDX2 = IDX1
             NUSE = 1

          ! Daily data
          ELSEIF ( Lct%Dct%Dta%ncDys(1) /= Lct%Dct%Dta%ncDys(2) ) THEN
             ! Error check
             IF ( Lct%Dct%Dta%ncHrs(1) /= Lct%Dct%Dta%ncHrs(2) ) THEN
                MSG = 'Data must only have one time dimension: ' // &
                      TRIM(Lct%Dct%cName)
                CALL HCO_ERROR( MSG, RC, THISLOC=LOC )
                RETURN
             ENDIF

             CALL GetSliceIdx ( HcoState, Lct, 3, prefDy, IDX1, RC )
             IF ( RC /= HCO_SUCCESS ) THEN
                 CALL HCO_ERROR( 'ERROR 23', RC, THISLOC=LOC )
                 RETURN
             ENDIF
             IDX2 = IDX1
             NUSE = 1

          ! All other cases (incl. hourly data): read all time slices).
          ELSE
             IDX1 = 1
             IDX2 = N
             NUSE = N
          ENDIF
       ENDIF

       ! ----------------------------------------------------------------
       ! Read selected time slice(s) into data array
       ! ----------------------------------------------------------------
       IF ( IDX2 > N ) THEN
          WRITE(MSG,*) 'Index ', IDX2, ' is larger than number of ', &
                       'values found: ', TRIM(Lct%Dct%cName)
          CALL HCO_ERROR( MSG, RC, THISLOC=LOC )
          RETURN
       ENDIF

       ALLOCATE( FileArr(1,1,1,NUSE), STAT=AS )
       IF ( AS /= 0 ) THEN
          MSG = 'Cannot allocate FileArr'
          CALL HCO_ERROR( MSG, RC, THISLOC=LOC )
          RETURN
       ENDIF

       ! Check for range/exact flag
       ! If range is given, the preferred Yr/Mt/Dy/Hr will be negative
       ! if we are outside the desired range.
       IF ( Lct%Dct%Dta%CycleFlag == HCO_CFLAG_RANGE ) THEN
          IF ( prefYr == -1 .OR. prefMt == -1 .OR. prefDy == -1 ) IDX1 = -1
          IF ( Lct%Dct%Dta%ncHrs(1) >= 0 .AND. prefHr == -1 )     IDX1 = -1

       ! If flag is exact, the preferred date must be equal to the current
       ! simulation date.
       ELSEIF ( Lct%Dct%Dta%CycleFlag == HCO_CFLAG_EXACT ) THEN
          IF ( Lct%Dct%Dta%ncYrs(1) > 0 ) THEN
             IF ( prefYr < Lct%Dct%Dta%ncYrs(1) .OR. &
                  prefYr > Lct%Dct%Dta%ncYrs(2) ) IDX1 = -1
          ENDIF
          IF ( Lct%Dct%Dta%ncMts(1) > 0 ) THEN
             IF ( prefMt < Lct%Dct%Dta%ncMts(1) .OR. &
                  prefMt > Lct%Dct%Dta%ncMts(2) ) IDX1 = -1
          ENDIF
          IF ( Lct%Dct%Dta%ncDys(1) > 0 ) THEN
             IF ( prefDy < Lct%Dct%Dta%ncDys(1) .OR. &
                  prefDy > Lct%Dct%Dta%ncDys(2) ) IDX1 = -1
          ENDIF
          IF ( Lct%Dct%Dta%ncHrs(1) >= 0 ) THEN
             IF ( prefHr < Lct%Dct%Dta%ncHrs(1) .OR. &
                  prefHr > Lct%Dct%Dta%ncHrs(2) ) IDX1 = -1
          ENDIF
       ENDIF

       ! IDX1 becomes -1 for data that is outside of the valid range
       ! (and no time cycling enabled). In this case, make sure that
       ! scale factor is set to zero.
       IF ( IDX1 < 0 ) THEN
          IF ( Lct%Dct%DctType == HCO_DCTTYPE_BASE ) THEN
             FileArr(1,1,1,:) = 0.0_hp
             MSG = 'Base field outside of range - set to zero: ' // &
                   TRIM(Lct%Dct%cName)
             CALL HCO_WARNING ( HcoState%Config%Err, MSG, RC, WARNLEV=1, &
                                THISLOC=LOC )
#if defined( MODEL_GEOS )
          ELSEIF ( Lct%Dct%DctType == HCO_DCTTYPE_MASK ) THEN
             FileArr(1,1,1,:) = 0.0_hp
             MSG = 'Mask outside of range - set to zero: ' // &
                   TRIM(Lct%Dct%cName)
             CALL HCO_WARNING ( HcoState%Config%Err, MSG, RC, WARNLEV=1, &
                                THISLOC=LOC )
#endif
          ELSE
             FileArr(1,1,1,:) = 1.0_hp
             MSG = 'Scale factor outside of range - set to one: ' // &
                   TRIM(Lct%Dct%cName)
             CALL HCO_WARNING ( HcoState%Config%Err, MSG, RC, WARNLEV=1, &
                                THISLOC=LOC )
          ENDIF
       ELSE
          FileArr(1,1,1,:) = FileVals(IDX1:IDX2)
       ENDIF

       ! ----------------------------------------------------------------
       ! Convert data to HEMCO units
       ! ----------------------------------------------------------------
       CALL HCO_UNIT_CHANGE( HcoConfig     = HcoState%Config,            &
                             Array         = FileArr,                    &
                             Units         = TRIM(Lct%Dct%Dta%OrigUnit), &
                             MW            = MW_g,                       &
                             YYYY          = -999,                       &
                             MM            = -999,                       &
                             AreaFlag      = AreaFlag,                   &
                             TimeFlag      = TimeFlag,                   &
                             FACT          = UnitFactor,                 &
                             RC            = RC                           )
       IF ( RC /= HCO_SUCCESS ) THEN
           CALL HCO_ERROR( 'ERROR 24', RC, THISLOC=LOC )
           RETURN
       ENDIF

       ! Verbose mode
       IF ( UnitFactor /= 1.0_hp ) THEN
          IF ( HCO_IsVerb(HcoState%Config%Err,1) ) THEN
             WRITE(MSG,*) 'Data was in units of ', TRIM(Lct%Dct%Dta%OrigUnit), &
                          ' - converted to HEMCO units by applying ', &
                          'scale factor ', UnitFactor
             CALL HCO_MSG(HcoState%Config%Err,MSG)
          ENDIF
       ELSE
          IF ( HCO_IsVerb(HcoState%Config%Err,2) ) THEN
             WRITE(MSG,*) 'Data was in units of ', TRIM(Lct%Dct%Dta%OrigUnit), &
                          ' - unit conversion factor is ', UnitFactor
             CALL HCO_MSG(HcoState%Config%Err,MSG)
          ENDIF
       ENDIF

       ! Data must be ...
       ! ... concentration ...
       IF ( AreaFlag == 3 .AND. TimeFlag == 0 ) THEN
          Lct%Dct%Dta%IsConc = .TRUE.

       ELSEIF ( AreaFlag == 3 .AND. TimeFlag == 1 ) THEN
          Lct%Dct%Dta%IsConc = .TRUE.
          FileArr = FileArr * HcoState%TS_EMIS
          MSG = 'Data converted from kg/m3/s to kg/m3: ' // &
                TRIM(Lct%Dct%cName) // ': ' // TRIM(Lct%Dct%Dta%OrigUnit)
          CALL HCO_WARNING ( HcoState%Config%Err, MSG, RC, WARNLEV=1, &
                             THISLOC=LOC )

       ! ... emissions or unitless ...
       ELSEIF ( (AreaFlag == -1 .AND. TimeFlag == -1) .OR. &
                (AreaFlag ==  2 .AND. TimeFlag ==  1)       ) THEN
          Lct%Dct%Dta%IsConc = .FALSE.

       ! ... invalid otherwise:
       ELSE
          MSG = 'Unit must be unitless, emission or concentration: ' // &
                TRIM(Lct%Dct%cName) // ': ' // TRIM(Lct%Dct%Dta%OrigUnit)
          CALL HCO_ERROR( MSG, RC, THISLOC=LOC )
          RETURN
       ENDIF

       ! Auto-detect delta t [in hours] between time slices.
       ! Scale factors can be:
       ! length 1 : constant
       ! length 7 : weekday factors: Sun, Mon, ..., Sat
       ! length 12: monthly factors: Jan, Feb, ..., Dec
       ! length 24: hourly  factors: 12am, 1am, ... 11pm
       IF ( NUSE == 1 ) THEN
          Lct%Dct%Dta%DeltaT = 0
       ELSEIF ( NUSE == 7 ) THEN
          Lct%Dct%Dta%DeltaT = 24
       ELSEIF ( NUSE == 12 ) THEN
          Lct%Dct%Dta%DeltaT = 720
       ELSEIF ( NUSE == 24 ) THEN
          Lct%Dct%Dta%DeltaT = 1
       ELSE
          MSG = 'Factor must be of length 1, 7, 12, or 24!' // &
                 TRIM(Lct%Dct%cName)
          CALL HCO_ERROR( MSG, RC, THISLOC=LOC)
          RETURN
       ENDIF

    ENDIF ! Masks vs. non-masks

    ! Copy data into output array.
    IF ( ASSOCIATED(Vals) ) DEALLOCATE( Vals )
    ALLOCATE( Vals(NUSE), STAT=AS )
    IF ( AS /= 0 ) THEN
       MSG = 'Cannot allocate Vals'
       CALL HCO_ERROR( MSG, RC, THISLOC=LOC )
       RETURN
    ENDIF
    Vals(:) = FileArr(1,1,1,:)

    ! Cleanup
    IF ( ASSOCIATED(FileArr) ) DEALLOCATE(FileArr)
    FileArr => NULL()

    ! Return w/ success
    RC = HCO_SUCCESS

  END SUBROUTINE GetDataVals
!EOC
!------------------------------------------------------------------------------
!                   Harmonized Emissions Component (HEMCO)                    !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: FillMaskBox
!
! !DESCRIPTION: Subroutine FillMaskBox fills the data array of the passed list
! container Lct according to the mask region provided in Vals. Vals contains
! the mask region of interest, denoted by the lower left and upper right grid
! box corners: lon1, lat1, lon2, lat2. The data array of Lct is filled such
! that all grid boxes are set to 1 whose mid-point is inside of the given box
! range.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE FillMaskBox ( HcoState, Lct, Vals, RC )
!
! !USES:
!
!
! !INPUT PARAMTERS:
!
    TYPE(HCO_State),  POINTER         :: HcoState    ! HEMCO state
    REAL(hp)        , POINTER         :: Vals(:)
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(ListCont),   POINTER         :: Lct
    INTEGER,          INTENT(INOUT)   :: RC
!
! !REVISION HISTORY:
!  29 Dec 2014 - C. Keller - Initial version
!  See https://github.com/geoschem/hemco for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    LOGICAL            :: GridPoint
    INTEGER            :: I, J
    REAL(hp)           :: LON1, LON2, LAT1, LAT2
    REAL(hp)           :: XDG1, XDG2, YDG1, YDG2
    REAL(hp)           :: ILON, ILAT
    CHARACTER(LEN=255) :: MSG
    CHARACTER(LEN=255) :: LOC = 'FillMaskBox (hcoio_util_mod.F90)'

    !=================================================================
    ! FillMaskBox begins here!
    !=================================================================

    ! Extract lon1, lon2, lat1, lat2
    LON1 = VALS(1)
    LAT1 = VALS(2)
    LON2 = VALS(3)
    LAT2 = VALS(4)

    ! Check if this is mask is a point. In this case, we need the grid
    ! box edges being defined.
    GridPoint = .FALSE.
    IF ( ( LON1 == LON2 ) .AND. ( LAT1 == LAT2 ) ) THEN
       IF ( .NOT. ASSOCIATED(HcoState%Grid%XEDGE%Val) .OR. &
            .NOT. ASSOCIATED(HcoState%Grid%YEDGE%Val)       ) THEN
          MSG = 'Cannot evaluate grid point mask - need grid box '   // &
                'edges for this. This error occurs if a mask covers '// &
                'a fixed grid point (e.g. lon1=lon2 and lat1=lat2) ' // &
                'but HEMCO grid edges are not defined.'
          CALL HCO_ERROR( MSG, RC, THISLOC=LOC )
          RETURN
       ENDIF
       GridPoint = .TRUE.
    ENDIF

    ! Check for every grid box if mid point is within mask region.
    ! Set to 1.0 if this is the case.
!$OMP PARALLEL DO                        &
!$OMP DEFAULT( SHARED                 )  &
!$OMP PRIVATE( I, J, ILON, ILAT       )  &
!$OMP PRIVATE( XDG1, XDG2, YDG1, YDG2 )  &
!$OMP SCHEDULE( DYNAMIC               )
    DO J = 1, HcoState%NY
    DO I = 1, HcoState%NX

       ! If it's a grid point, check if it's within this
       ! grid box
       IF ( GridPoint ) THEN
          XDG1 = HcoState%Grid%XEDGE%Val(I  ,J  )
          XDG2 = HcoState%Grid%XEDGE%Val(I+1,J  )
          YDG1 = HcoState%Grid%YEDGE%Val(I  ,J  )
          YDG2 = HcoState%Grid%YEDGE%Val(I  ,J+1)
          IF ( XDG1 >= 180.0_hp ) XDG1 = XDG1 - 360.0_hp
          IF ( XDG2 >= 180.0_hp ) XDG2 = XDG2 - 360.0_hp

          IF ( LON1 >= XDG1 .AND. LON1 <= XDG2 .AND. &
               LAT1 >= YDG1 .AND. LAT1 <= YDG2        ) THEN
             Lct%Dct%Dta%V2(1)%Val(I,J) = 1.0_sp
          ENDIF

       ! Check if mid point is within mask region
       ELSE
          ! Get longitude and latitude at this grid box
          ILON = HcoState%Grid%XMID%Val(I,J)
          ILAT = HcoState%Grid%YMID%Val(I,J)
          IF ( ILON >= 180.0_hp ) ILON = ILON - 360.0_hp

          IF ( ILON >= LON1 .AND. ILON <= LON2 .AND. &
               ILAT >= LAT1 .AND. ILAT <= LAT2        ) THEN
             Lct%Dct%Dta%V2(1)%Val(I,J) = 1.0_sp
          ENDIF
       ENDIF

    ENDDO
    ENDDO
!$OMP END PARALLEL DO

    ! Return w/ success
    RC = HCO_SUCCESS

  END SUBROUTINE FillMaskBox
!EOC
!------------------------------------------------------------------------------
!                   Harmonized Emissions Component (HEMCO)                    !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: ReadMath
!
! !DESCRIPTION: Subroutine ReadMath reads and evaluates a mathematical
! expression. Mathematical expressions can combine time-stamps with
! mathematical functions, e.g. to yield the sine of current simulation hour.
! Mathematical expressions must start with the identifier 'MATH:', followed
! by the actual expression. Each expression must include at least one
! variable (evaluated at runtime). The following variables are currently
! supported: YYYY (year), MM (month), DD (day), HH (hour), LH (local hour),
! NN (minute), SS (second), WD (weekday), LWD (local weekday),
! DOY (day of year), ELH (elapsed hours), ELS (elapsed seconds).
! In addition, the following variables can be used: PI (3.141...), DOM
! (\# of days of current month).
! For example, the following expression would yield a continuous sine
! curve as function of hour of day: 'MATH:sin(HH/24*PI*2)'.
!\\
!\\
! For a full list of valid mathematical expressions, see module interpreter.F90.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE ReadMath( HcoState, Lct, ValStr, Vals, N, RC )
!
! !USES:
!
    USE HCO_CLOCK_MOD,      ONLY : HcoClock_Get
    USE HCO_tIdx_Mod,       ONLY : HCO_GetPrefTimeAttr
    USE INTERPRETER
!
! !INPUT PARAMTERS:
!
    TYPE(HCO_State),  POINTER         :: HcoState    ! HEMCO state
    TYPE(ListCont),   POINTER         :: Lct
    CHARACTER(LEN=*), INTENT(IN   )   :: ValStr
!
! !INPUT/OUTPUT PARAMETERS:
!
    REAL(hp),         INTENT(INOUT)   :: Vals(:)
    INTEGER,          INTENT(INOUT)   :: RC
!
! !OUTPUT PARAMETERS:
!
    INTEGER,          INTENT(  OUT)   :: N
!
! !REVISION HISTORY:
!  11 May 2017 - C. Keller - Initial version
!  See https://github.com/geoschem/hemco for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    LOGICAL            :: EOS
    INTEGER            :: STRL
    INTEGER            :: I, NVAL, LHIDX, LWDIDX
    INTEGER            :: prefYr, prefMt, prefDy, prefHr, prefMn
    INTEGER            :: prefWD, prefDOY, prefS, LMD, cHr
    INTEGER            :: nSteps
    REAL(hp)           :: ELH, ELS
    REAL(hp)           :: Val
    CHARACTER(LEN=255) :: MSG
    CHARACTER(LEN=255) :: LOC = 'ReadMath (hcoio_util_mod.F90)'

    ! Variables used by the evaluator to build and to determine the value
    ! of the expressions
    character(len = 10) :: all_variables(12)
    real(hp)            :: all_variablesvalues(12)

    !String variable that will store the function that the evaluator will build
    character (len = 275)  :: func

    !String variable that will return the building of the expression result
    !If everything was ok then statusflag = 'ok', otherwise statusflag = 'error'
    character (len = 5)  :: statusflag

    !======================================================================
    ! ReadMath begins here
    !======================================================================

    ! Substring (without flag 'MATH:')
    STRL = LEN(ValStr)
    IF ( STRL < 6 ) THEN
       MSG = 'Math expression is too short - expected `MATH:<expr>`: ' &
             //TRIM(ValStr)
       CALL HCO_ERROR( MSG, RC, THISLOC=LOC )
       RETURN
    ENDIF
    func = ValStr(6:STRL)

    ! Get preferred time stamps
    CALL HCO_GetPrefTimeAttr( HcoState, Lct, &
                              prefYr, prefMt, prefDy, prefHr, prefMn, RC )
    IF ( RC /= HCO_SUCCESS ) THEN
        CALL HCO_ERROR( 'ERROR 25', RC, THISLOC=LOC )
        RETURN
    ENDIF

    ! Get some other current time stamps
    CALL HcoClock_Get( HcoState%Clock,  cS=prefS,     cH=cHr, &
                       cWEEKDAY=prefWD, cDOY=prefDOY, LMD=LMD,      &
                       nSteps=nSteps,   RC=RC )
    IF ( RC /= HCO_SUCCESS ) THEN
        CALL HCO_ERROR( 'ERROR 26', RC, THISLOC=LOC )
        RETURN
    ENDIF

    ! GetPrefTimeAttr can return -999 for hour. In this case set to current
    ! simulation hour
    IF ( prefHr < 0 ) prefHr = cHr

    ! Parse function. This will replace any tokens in the function with the
    ! actual token values. (ckeller, 7/7/17)
    CALL HCO_CharParse ( HcoState%Config, func, &
                         prefYr, prefMt, prefDy, prefHr, prefMn, RC )
    IF ( RC /= HCO_SUCCESS ) THEN
        CALL HCO_ERROR( 'ERROR 27', RC, THISLOC=LOC )
        RETURN
    ENDIF

    ! Elapsed hours and seconds since start time
    ELS = HcoState%TS_DYN * nSteps
    ELH = ELS / 3600.0_hp

    ! Check which variables are in string.
    ! Possible variables are YYYY, MM, DD, WD, HH, NN, SS, DOY, ELH, ELS
    NVAL   = 0
    LHIDX  = -1
    LWDIDX = -1

    IF ( INDEX(func,'YYYY') > 0 ) THEN
       NVAL                      = NVAL + 1
       all_variables(NVAL)       = 'yyyy'
       all_variablesvalues(NVAL) = prefYr
    ENDIF
    IF ( INDEX(func,'MM') > 0 ) THEN
       NVAL                      = NVAL + 1
       all_variables(NVAL)       = 'mm'
       all_variablesvalues(NVAL) = prefMt
    ENDIF
    IF ( INDEX(func,'DD') > 0 ) THEN
       NVAL                      = NVAL + 1
       all_variables(NVAL)       = 'dd'
       all_variablesvalues(NVAL) = prefDy
    ENDIF
    IF ( INDEX(func,'WD') > 0 ) THEN
       NVAL                      = NVAL + 1
       all_variables(NVAL)       = 'wd'
       all_variablesvalues(NVAL) = prefWD
    ENDIF
    IF ( INDEX(func,'LWD') > 0 ) THEN
       NVAL                      = NVAL + 1
       all_variables(NVAL)       = 'lwd'
       all_variablesvalues(NVAL) = prefWD
       LWDIDX                    = NVAL
    ENDIF
    IF ( INDEX(func,'HH') > 0 ) THEN
       NVAL                      = NVAL + 1
       all_variables(NVAL)       = 'hh'
       all_variablesvalues(NVAL) = prefHr
    ENDIF
    IF ( INDEX(func,'LH') > 0 ) THEN
       NVAL                      = NVAL + 1
       all_variables(NVAL)       = 'lh'
       all_variablesvalues(NVAL) = prefHr
       LHIDX                     = NVAL
    ENDIF
    IF ( INDEX(func,'NN') > 0 ) THEN
       NVAL                      = NVAL + 1
       all_variables(NVAL)       = 'nn'
       all_variablesvalues(NVAL) = prefMn
    ENDIF
    IF ( INDEX(func,'SS') > 0 ) THEN
       NVAL                      = NVAL + 1
       all_variables(NVAL)       = 'ss'
       all_variablesvalues(NVAL) = prefS
    ENDIF
    IF ( INDEX(func,'DOY') > 0 ) THEN
       NVAL                      = NVAL + 1
       all_variables(NVAL)       = 'doy'
       all_variablesvalues(NVAL) = prefDOY
    ENDIF
    IF ( INDEX(func,'PI') > 0 ) THEN
       NVAL                      = NVAL + 1
       all_variables(NVAL)       = 'pi'
       all_variablesvalues(NVAL) = HcoState%Phys%PI
    ENDIF
    IF ( INDEX(func,'DOM') > 0 ) THEN
       NVAL                      = NVAL + 1
       all_variables(NVAL)       = 'dom'
       all_variablesvalues(NVAL) = LMD
    ENDIF
    IF ( INDEX(func,'ELH') > 0 ) THEN
       NVAL                      = NVAL + 1
       all_variables(NVAL)       = 'elh'
       all_variablesvalues(NVAL) = ELH
    ENDIF
    IF ( INDEX(func,'ELS') > 0 ) THEN
       NVAL                      = NVAL + 1
       all_variables(NVAL)       = 'els'
       all_variablesvalues(NVAL) = ELS
    ENDIF

    ! Error trap: cannot have local hour and local weekday in
    ! same expression
    IF ( LHIDX > 0 .AND. LWDIDX > 0 ) THEN
       MSG = 'Cannot have local hour and local weekday in '//&
             'same expression: '//TRIM(func)
       CALL HCO_ERROR( MSG, RC, THISLOC=LOC )
       RETURN
    ENDIF

    ! N is the number of expressions.
    Vals(:) = -999.0_hp
    IF ( LHIDX > 0 ) THEN
       N = 24
    ELSEIF ( LWDIDX > 0 ) THEN
       N = 7
    ELSE
       N = 1
    ENDIF

    ! Evaluate expression
    !Initialize function
    call init (func, all_variables(1:NVAL), statusflag)
    IF(statusflag == 'ok') THEN
       DO I=1,N
          IF ( LHIDX  > 0 ) all_variablesvalues(LHIDX)  = I-1
          IF ( LWDIDX > 0 ) all_variablesvalues(LWDIDX) = I-1
          Val = evaluate( all_variablesvalues(1:NVAL) )
          Vals(I) = Val

          ! Verbose
          IF ( HCO_IsVerb(HcoState%Config%Err,2) ) THEN
             WRITE(MSG,*) 'Evaluated function: ',TRIM(func),' --> ', Val
             CALL HCO_MSG(HcoState%Config%Err,MSG)
          ENDIF
       ENDDO
    ELSE
       MSG = 'Error evaluation function: '//TRIM(func)
       CALL HCO_ERROR( MSG, RC, THISLOC=LOC )
       RETURN
    ENDIF
    call destroyfunc()

    ! Return w/ success
    RC = HCO_SUCCESS

  END SUBROUTINE ReadMath
!EOC
END MODULE HCOIO_Util_Mod
