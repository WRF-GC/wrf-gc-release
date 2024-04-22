!------------------------------------------------------------------------------ 
!                   Harmonized Emissions Component (HEMCO)                    ! 
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: hco_interface_common
!
! !DESCRIPTION: Module HCO\_Interface\_Common defines common utilities for interfacing
!  with HEMCO from other models. The common toolbox should be present for all
!  models interacting with HEMCO.
!\\
!\\
! !INTERFACE:
!
MODULE HCO_Interface_Common
!
! !USES:
!
  USE HCO_Error_Mod
  USE HCO_State_Mod,  ONLY: HCO_State
  USE HCOX_State_Mod, ONLY: Ext_State

  IMPLICIT NONE
  PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!
  PUBLIC  :: SetHcoTime
  PUBLIC  :: GetHcoVal
  PUBLIC  :: GetHcoDiagn
!
! !REMARKS:
!  These utilities were mostly migrated from GEOS-Chem HCO_Interface_Mod.
!  All functions now accept as input the HEMCO state (HcoState) as there may
!  be multiple instances and all variables should not be inferred.
!
! !REVISION HISTORY:
!  12 Mar 2020 - H.P. Lin    - Initial version
!  See https://github.com/geoschem/hemco for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
CONTAINS
!EOC
!------------------------------------------------------------------------------
!                    Harmonized Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: SetHcoTime
!
! !DESCRIPTION: SUBROUTINE SetHcoTime sets the current simulation
! datetime in HcoState.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE SetHcoTime( HcoState,  ExtState,   year,   month,  &
                         day,       dayOfYr,    hour,   minute, &
                         second,    IsEmisTime, RC             )
!
! !USES:
!
    USE HCO_CLOCK_MOD, ONLY : HcoClock_Set
!
! !INPUT PARAMETERS:
!
    TYPE(HCO_State), POINTER       :: HcoState    ! HEMCO state object
    TYPE(Ext_State), POINTER       :: ExtState    ! HEMCO extensions state object
    INTEGER,         INTENT(IN)    :: year        ! UTC year 
    INTEGER,         INTENT(IN)    :: month       ! UTC month
    INTEGER,         INTENT(IN)    :: day         ! UTC day
    INTEGER,         INTENT(IN)    :: dayOfYr     ! UTC day of year
    INTEGER,         INTENT(IN)    :: hour        ! UTC hour
    INTEGER,         INTENT(IN)    :: minute      ! UTC minute
    INTEGER,         INTENT(IN)    :: second      ! UTC second
    LOGICAL,         INTENT(IN)    :: IsEmisTime  ! Time for emissions?
!
! !INPUT/OUTPUT PARAMETERS:
!
    INTEGER,         INTENT(INOUT) :: RC
!
! !REVISION HISTORY:
!  23 Oct 2012 - C. Keller - Initial Version
!  See https://github.com/geoschem/hemco for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! LOCAL VARIABLES:
!
    !=================================================================
    ! SetHcoTime begins here
    !=================================================================

    CALL HcoClock_Set ( HcoState, year, month, day, hour, minute, &
                        second, dayOfYr, IsEmisTime=IsEmisTime, RC=RC )

  END SUBROUTINE SetHcoTime
!EOC
!------------------------------------------------------------------------------
!                    Harmonized Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: GetHcoVal
!
! !DESCRIPTION: Subroutine GetHcoVal is a wrapper routine to return an
! emission (kg/m2/s) or deposition (1/s) value from the HEMCO state object
! for a given species at position I, J, L.
! A value of zero is returned if no HEMCO species is defined for the given
! tracer, and the output parameter Found is set to false.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE GetHcoVal ( HcoState, ExtState, &
                         HcoID, I, J, L, Found, Emis, Dep )
!
! !USES:
!
!
! !INPUT ARGUMENTS:
!
    TYPE(HCO_State),    POINTER        :: HcoState    ! HEMCO state object
    TYPE(Ext_State),    POINTER        :: ExtState    ! HEMCO extensions state object
    INTEGER,            INTENT(IN   )  :: HcoID       ! HEMCO tracer ID
    INTEGER,            INTENT(IN   )  :: I, J, L     ! Position
!
! !OUTPUT ARGUMENTS:
!
    LOGICAL,            INTENT(  OUT)  :: Found       ! Was this tracer ID found?
    REAL(hp), OPTIONAL, INTENT(  OUT)  :: Emis        ! Emissions  [kg/m2/s]
    REAL(hp), OPTIONAL, INTENT(  OUT)  :: Dep         ! Deposition [1/s]
!
! !REMARKS:
!  Tracer ID passed is now always HcoID. This is because GEOS-Chem uses the same
!  HcoID = TrcID. If your model is not, please do a mapping internally in the
!  interface.
!
!  Note that HEMCO expects the grid to be 3-D in IJL indices. If your model
!  (e.g. CAM) stores data in 2-D columns (K, I) where I is a chunked set of
!  columns, one dummy dimension needs to be added. Refer to ESCOMP/HEMCO_CESM 
!  (hplin, 3/12/20)
!
! !REVISION HISTORY:
!  20 Oct 2014 - C. Keller - Initial Version
!  See https://github.com/geoschem/hemco for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
    !=================================================================
    ! GetHcoVal begins here
    !=================================================================

    ! Init
    FOUND = .FALSE.
    IF ( PRESENT(Emis) ) Emis = 0.0_hp
    IF ( PRESENT(Dep ) ) Dep  = 0.0_hp

    ! If HEMCO species exists, get value from HEMCO state
    IF ( HcoID > 0 ) THEN
       IF ( PRESENT(Emis) ) THEN
          IF ( ASSOCIATED(HcoState%Spc(HcoID)%Emis%Val) ) THEN
             Emis  = HcoState%Spc(HcoID)%Emis%Val(I,J,L)
             FOUND = .TRUE.
          ENDIF
       ENDIF
       IF ( PRESENT(Dep) ) THEN
          IF ( ASSOCIATED(HcoState%Spc(HcoID)%Depv%Val) ) THEN
             Dep   = HcoState%Spc(HcoID)%Depv%Val(I,J)
             FOUND = .TRUE.
          ENDIF
       ENDIF
    ENDIF

  END SUBROUTINE GetHcoVal
!EOC
!------------------------------------------------------------------------------
!                    Harmonized Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: GetHcoDiagn
!
! !DESCRIPTION: Subroutine GetHcoDiagn is a convenience wrapper routine to
!  get a HEMCO diagnostics from an external model.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE GetHcoDiagn ( HcoState,       ExtState,  DiagnName, &
                           StopIfNotFound, RC,        Ptr2D,     &
                           Ptr3D,          COL,       AutoFill  )
!
! !USES:
!
    USE HCO_TYPES_MOD,      ONLY : DiagnCont
    USE HCO_DIAGN_MOD
!
! !INPUT PARAMETERS:
!
    TYPE(HCO_State), POINTER               :: HcoState       ! HEMCO state object
    TYPE(Ext_State), POINTER               :: ExtState       ! HEMCO extensions state object
    CHARACTER(LEN=*), INTENT(IN)           :: DiagnName      ! Name of diagnostics
    LOGICAL,          INTENT(IN)           :: StopIfNotFound ! Stop if diagnostics
                                                             ! does not exist?
    INTEGER,          INTENT(IN), OPTIONAL :: COL            ! Collection Nr.
    INTEGER,          INTENT(IN), OPTIONAL :: AutoFill       ! AutoFill diagnostics only?
!
! !INPUT/OUTPUT PARAMETERS:
!
    INTEGER,          INTENT(INOUT)        :: RC             ! Error return code
!
! !OUTPUT PARAMETERS:
!
    REAL(sp),         POINTER, OPTIONAL    :: Ptr2D(:,:)     ! Pointer to 2D data
    REAL(sp),         POINTER, OPTIONAL    :: Ptr3D(:,:,:)   ! Pointer to 3D data
!
! !REVISION HISTORY:
!  24 Sep 2014 - C. Keller   - Initial version
!  See https://github.com/geoschem/hemco for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER                   :: FLAG, LevIDx, PS, AF
    TYPE(DiagnCont), POINTER  :: DgnCont  => NULL()

    ! Strings
    CHARACTER(LEN=255) :: ErrMsg
    CHARACTER(LEN=255) :: ThisLoc

    !=======================================================================
    ! GetHcoDiagn begins here
    !=======================================================================

    ! Initialize
    RC      = HCO_SUCCESS

    ! For error handling
    ErrMsg  = ''

    ! Set collection number
    PS = HcoState%Diagn%HcoDiagnIDManual
    IF ( PRESENT(COL) ) PS = COL

    ! Set AutoFill flag
    AF = -1
    IF ( PRESENT(AutoFill) ) AF = AutoFill

    ! Get diagnostics by name. Search all diagnostics, i.e. both AutoFill
    ! and manually filled diagnostics. Also include those with a manual
    ! output interval.
    CALL Diagn_Get( HcoState, .FALSE., DgnCont, FLAG, RC, &
                    cName=TRIM(DiagnName), AutoFill=AF, COL=PS )

    ! Trap potential errors
    IF ( RC /= HCO_SUCCESS ) THEN
       ErrMsg = 'Error in getting diagnostics: ' // TRIM(DiagnName)
       CALL HCO_Error( ErrMsg, RC )
       RETURN
    ENDIF

    IF ( (FLAG /= HCO_SUCCESS) .AND. StopIfNotFound ) THEN
       ErrMsg = 'Cannot get diagnostics for this time stamp: ' //    &
                 TRIM(DiagnName)
       CALL HCO_Error( ErrMsg, RC )
       RETURN
    ENDIF

    ! Pass data to output pointer (only if diagnostics defined):
    IF ( FLAG == HCO_SUCCESS ) THEN

       ! 2D pointer
       IF ( PRESENT(Ptr2D) ) THEN

          ! Pass 2D data
          IF ( ASSOCIATED(DgnCont%Arr2D%Val) ) THEN
             Ptr2D => DgnCont%Arr2D%Val

          ! Pass 3D data. Get level index from diagnostics (if set)
          ELSEIF ( ASSOCIATED(DgnCont%Arr3D%Val) ) THEN
             LevIDx = DgnCont%LevIdx
             IF ( LevIdx < 1 ) LevIdx = 1
             Ptr2D => DgnCont%Arr3D%Val(:,:,LevIDx)

          ! Error if no 2D or 3D data available
          ELSE
             ErrMsg = 'no data defined: '// TRIM(DiagnName)
             CALL HCO_Error( ErrMsg, RC )
             RETURN
          ENDIF

       ! 3D pointer: must point to 3D data
       ELSEIF ( PRESENT(Ptr3D) ) THEN
          IF ( ASSOCIATED(DgnCont%Arr3D%Val) ) THEN
             Ptr3D => DgnCont%Arr3D%Val
          ELSE
             ErrMsg = 'no 3D data defined: '// TRIM(DiagnName)
             CALL HCO_Error( ErrMsg, RC )
             RETURN
          ENDIF

       ! Error otherwise
       ELSE
          ErrMsg = 'Please define output data pointer: ' // TRIM(DiagnName)
          CALL HCO_Error( ErrMsg, RC )
          RETURN
       ENDIF
    ENDIF

    ! Free pointer
    DgnCont  => NULL()

    ! Leave with success
    RC = HCO_SUCCESS

  END SUBROUTINE GetHcoDiagn

!EOC
END MODULE HCO_Interface_Common
