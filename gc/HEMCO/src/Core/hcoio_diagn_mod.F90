!EOC
!------------------------------------------------------------------------------
!                   Harmonized Emissions Component (HEMCO)                    !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: hcoio_diagn_mod.F90
!
! !DESCRIPTION: Module HCOIO\_Diagn\_Mod.F90 is the data interface module
! for the HEMCO diagnostics. It contains routines to write out diagnostics
! into a netCDF file.
!\\
!\\
! In an ESMF/MAPL environment, the HEMCO diagnostics are not directly
! written to disk but passed to the gridded component export state, where
! they can be picked up by the MAPL HISTORY component.
!\\
!\\
! !INTERFACE:
!
MODULE HCOIO_DIAGN_MOD
!
! !USES:
!
  USE HCO_ERROR_MOD

  IMPLICIT NONE
  PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!
  PUBLIC :: HcoDiagn_Write
  PUBLIC :: HCOIO_Diagn_WriteOut
!
! !REMARKS:
!  HEMCO diagnostics are still in testing mode. We will fully activate them
!  at a later time.  They will be turned on when debugging & unit testing.
!
! !REVISION HISTORY:
!  04 May 2014 - C. Keller   - Initial version.
!  See https://github.com/geoschem/hemco for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !DEFINED PARAMETERS:
!
  ! Fill value used in HEMCO diagnostics netCDF files.
!  REAL(hp), PARAMETER :: FillValue = 1.e-31_hp
  REAL(sp), PARAMETER :: FillValue = HCO_MISSVAL

CONTAINS
!EOC
!------------------------------------------------------------------------------
!                   Harmonized Emissions Component (HEMCO)                    !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: HcoDiagn_Write
!
! !DESCRIPTION: Subroutine HcoDiagn\_Write is the wrapper routine to write out
! the content of the built-in HEMCO diagnostics. If input argument Restart is
! set to TRUE, only the restart collection will be written out. Otherwise,
! the default collection
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE HcoDiagn_Write( HcoState, Restart, RC )
!
! !USES:
!
    USE HCO_State_Mod,        ONLY : HCO_State
    USE HCO_Clock_Mod,        ONLY : HcoClock_SetLast
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(HCO_State), POINTER          :: HcoState     ! HEMCO state object
    LOGICAL,         INTENT(IN   )    :: Restart      ! write restart (enforced)?
!
! !INPUT/OUTPUT PARAMETERS:
!
    INTEGER,         INTENT(INOUT)    :: RC           ! Return code
!
! !REVISION HISTORY:
!  03 Apr 2015 - C. Keller   - Initial version
!  See https://github.com/geoschem/hemco for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER            :: I, COL
    CHARACTER(LEN=255) :: MSG, LOC
#ifdef ADJOINT
    INTEGER            :: MaxIdx
#endif

    !=================================================================
    ! HcoDiagn_Write begins here!
    !=================================================================

    ! Init
    LOC = 'HcoDiagn_Write (hcoio_diagn_mod.F90)'

    ! To write restart (enforced)
    IF ( RESTART ) THEN
       CALL HCOIO_DIAGN_WRITEOUT ( HcoState,                               &
                                   ForceWrite  = .TRUE.,                   &
                                   UsePrevTime = .FALSE.,                  &
                                   COL = HcoState%Diagn%HcoDiagnIDRestart, &
                                   RC          = RC                         )
       IF( RC /= HCO_SUCCESS) RETURN

       ! Set last flag for use below
       CALL HcoClock_SetLast ( HcoState%Clock, .TRUE., RC )
       IF( RC /= HCO_SUCCESS) RETURN

       CALL HCOIO_DIAGN_WRITEOUT ( HcoState,                               &
                                   ForceWrite  = .FALSE.,                  &
                                   UsePrevTime = .FALSE.,                  &
                                   COL = HcoState%Diagn%HcoDiagnIDDefault, &
                                   RC          = RC                         )
       IF( RC /= HCO_SUCCESS) RETURN

#ifdef ADJOINT
       IF (HcoState%isAdjoint) THEN
       CALL HCOIO_DIAGN_WRITEOUT ( HcoState,                               &
                                   ForceWrite  = .FALSE.,                  &
                                   UsePrevTime = .FALSE.,                  &
                                   COL = HcoState%Diagn%HcoDiagnIDAdjoint, &
                                   RC          = RC                         )
       IF( RC /= HCO_SUCCESS) RETURN 
       ENDIF
#endif

       ! Reset IsLast flag. This is to ensure that the last flag is not
       ! carried over (ckeller, 11/1/16).
       CALL HcoClock_SetLast ( HcoState%Clock, .FALSE., RC )
       IF( RC /= HCO_SUCCESS) RETURN

    ELSE

       ! Loop over all collections that shall be written out.
       ! HCOIO_DIAGN_WRITEOUT will determine whether it is time to
       ! write a collection or not.
#ifndef ADJOINT
       DO I = 1, 3
#else
       MaxIdx = 3
       IF (HcoState%isAdjoint) MaxIdx = 4
       DO I = 1, MaxIdx
#endif

          ! Define collection ID
          SELECT CASE ( I )
             CASE ( 1 )
                COL = HcoState%Diagn%HcoDiagnIDDefault
             CASE ( 2 )
                COL = HcoState%Diagn%HcoDiagnIDRestart
             CASE ( 3 )
                COL = HcoState%Diagn%HcoDiagnIDManual
#ifdef ADJOINT
             CASE ( 4 )
                COL = HcoState%Diagn%HcoDiagnIDAdjoint
#endif
          END SELECT

#if       !defined ( ESMF_ )
          ! If not ESMF environment, never write the manual diagnostics
          ! to disk. Instead, the content of the manual diagnostics needs
          ! to be fetched explicitly.
          IF ( I == 3 ) CYCLE
#else
          ! Don't write restart variables to EXPORT in an ESMF environment.
          ! They are already passed to the INTERNAL state when calling
          ! HCO_RestartWrite. (ckeller, 10/9/17)
          IF ( I == 2 ) CYCLE
#endif

          ! Restart file
          CALL HCOIO_DIAGN_WRITEOUT ( HcoState,                        &
                                      ForceWrite  = .FALSE.,           &
                                      UsePrevTime = .FALSE.,           &
                                      COL         = COL,               &
                                      RC          = RC                  )
          IF(RC /= HCO_SUCCESS) RETURN
       ENDDO
    ENDIF

  END SUBROUTINE HcoDiagn_Write
!EOC
!------------------------------------------------------------------------------
!                   Harmonized Emissions Component (HEMCO)                    !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: HCOIO_Diagn_WriteOut
!
! !DESCRIPTION: Subroutine HCOIO\_Diagn\_WriteOut writes diagnostics to
! output. Depending on the model environment, different subroutines will
! be invoked.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE HCOIO_Diagn_WriteOut( HcoState, ForceWrite,  &
                                   RC,          PREFIX,   UsePrevTime, &
                                   OnlyIfFirst, COL                     )
!
! !USES:
!
    USE HCO_State_Mod,        ONLY : HCO_State
    USE HCOIO_Write_Mod,      ONLY : HCOIO_Write
!
! !INPUT PARAMETERS:
!
    TYPE(HCO_State),  POINTER                 :: HcoState    ! HEMCO state object
    LOGICAL,                    INTENT(IN   ) :: ForceWrite  ! Write all diagnostics?
    CHARACTER(LEN=*), OPTIONAL, INTENT(IN   ) :: PREFIX      ! File prefix
    LOGICAL,          OPTIONAL, INTENT(IN   ) :: UsePrevTime ! Use previous time
    LOGICAL,          OPTIONAL, INTENT(IN   ) :: OnlyIfFirst ! Only write if nnDiagn is 1
    INTEGER,          OPTIONAL, INTENT(IN   ) :: COL         ! Collection Nr.
!
! !INPUT/OUTPUT PARAMETERS:
!
    INTEGER,          INTENT(INOUT) :: RC          ! Failure or success
!
! !REVISION HISTORY:
!  12 Sep 2013 - C. Keller   - Initial version
!  See https://github.com/geoschem/hemco for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    CHARACTER(LEN=255), PARAMETER :: LOC = 'HCOIO_DIAGN_WRITEOUT (hcoio_diagn_mod.F90)'

    !=================================================================
    ! HCOIO_DIAGN_WRITEOUT begins here!
    !=================================================================

#if defined(ESMF_)
    !-----------------------------------------------------------------
    ! ESMF environment: call ESMF output routines
    !-----------------------------------------------------------------
    CALL HCOIO_Write     ( HcoState,                &
                           RC,                      &
                           OnlyIfFirst=OnlyIfFirst, &
                           COL=COL                   )
    IF ( RC /= HCO_SUCCESS ) THEN
        CALL HCO_ERROR( 'ERROR 0', RC, THISLOC=LOC )
        RETURN
    ENDIF

#else
    !-----------------------------------------------------------------
    ! Standard environment: call default output routines
    !-----------------------------------------------------------------
    CALL HCOIO_Write    ( HcoState,                 &
                          ForceWrite,               &
                          RC,                       &
                          PREFIX      =PREFIX,      &
                          UsePrevTime =UsePrevTime, &
                          OnlyIfFirst =OnlyIfFirst, &
                          COL         = COL          )
    IF ( RC /= HCO_SUCCESS ) THEN
        CALL HCO_ERROR( 'ERROR 1', RC, THISLOC=LOC )
        RETURN
    ENDIF

#endif

    ! Return
    RC = HCO_SUCCESS

  END SUBROUTINE HCOIO_DIAGN_WRITEOUT
!EOC
END MODULE HCOIO_Diagn_Mod

