!BOC
#if defined ( ESMF_ )
! The 'standard' HEMCO I/O module is used for:
! - GEOS-Chem High Performance / GCHP and GEOS (ESMF_)
!EOC
!------------------------------------------------------------------------------
!                   Harmonized Emissions Component (HEMCO)                    !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: hcoio_write_esmf_mod.F90
!
! !DESCRIPTION: Module HCOIO\_Write\_Mod.F90 is the HEMCO output
! interface for the ESMF environment.
! In an ESMF/MAPL environment, the HEMCO diagnostics are not directly
! written to disk but passed to the gridded component export state, where
! they can be picked up by the MAPL HISTORY component.
!\\
!\\
! !INTERFACE:
!
MODULE HCOIO_Write_Mod
!
! !USES:
!
  USE HCO_ERROR_MOD
  USE HCO_DIAGN_MOD

  IMPLICIT NONE
  PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!
  PUBLIC :: HCOIO_Write
!
! !REMARKS:
!  HEMCO diagnostics are still in testing mode. We will fully activate them
!  at a later time.  They will be turned on when debugging & unit testing.
!
! !REVISION HISTORY:
!  04 May 2014 - C. Keller   - Initial version
!  See https://github.com/geoschem/hemco for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !DEFINED PARAMETERS:
!
CONTAINS
!EOC
!------------------------------------------------------------------------------
!                   Harmonized Emissions Component (HEMCO)                    !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: HCOIO_Write
!
! !DESCRIPTION: Subroutine HCOIO\_Diagn\_WriteOut is the interface routine to
! link the HEMCO diagnostics arrays to the corresponding data pointers of the
! MAPL/ESMF history component.
!\\
!\\
! Since the history component internally organizes many diagnostics tasks such
! as output scheduling, file writing, and data averaging, all HEMCO diagnostics
! are made available to the history component on every time step, e.g. the
! entire content of the HEMCO diagnostics list is 'flushed' every time this
! subroutine is called.
!\\
!\\
! For now, all diagnostics data is copied to the corresponding MAPL data
! pointer so that this routine works for cases where the HEMCO precision is
! not equal to the ESMF precision.
!\\
!\\
! Once the HEMCO precision is pegged to the ESMF precision, we can just
! establish pointers between the export arrays and the diagnostics the first
! time this routine is called.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE HCOIO_WRITE ( HcoState, RC, OnlyIfFirst, COL )
!
! !USES:
!
    USE ESMF
    USE MAPLBase_MOD
    USE HCO_Types_Mod, ONLY : DiagnCont
    USE HCO_State_Mod, ONLY : HCO_State

# include "MAPL_Generic.h"
!
! !INPUT PARAMETERS:
!
    TYPE(HCO_State),  POINTER                 :: HcoState    ! HEMCO state object
    LOGICAL,          OPTIONAL, INTENT(IN   ) :: OnlyIfFirst !
    INTEGER,          OPTIONAL, INTENT(IN   ) :: COL         ! Collection Nr.
!
! !INPUT/OUTPUT PARAMETERS:
!

    INTEGER,                    INTENT(INOUT) :: RC          ! Failure or success
!
! !REVISION HISTORY:
!  05 Aug 2014 - C. Keller    - Initial version
!  See https://github.com/geoschem/hemco for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    TYPE(DiagnCont), POINTER  :: ThisDiagn
    INTEGER                   :: PS, FLAG, STAT
    CHARACTER(LEN=255)        :: MSG
    LOGICAL                   :: EOI
    REAL, POINTER             :: Ptr2D(:,:)
    REAL, POINTER             :: Ptr3D(:,:,:)

    CHARACTER(LEN=255), PARAMETER :: LOC = 'HCOIO_WRITE_ESMF (hcoio_write_esmf_mod.F90)'

    !=================================================================
    ! HCOIO_WRITE_ESMF begins here!
    !=================================================================

    ! Assume success until otherwise
    RC  = HCO_SUCCESS

    ! Init
    ThisDiagn => NULL()
    Ptr2D     => NULL()
    Ptr3D     => NULL()

    ! Collection number
    PS = HcoState%Diagn%HcoDiagnIDDefault
    IF ( PRESENT(COL) ) PS = COL

    ! In an ESMF environment, always get all diagnostics since output
    ! is scheduled through MAPL History!
    EOI = .FALSE.

    !-----------------------------------------------------------------
    ! Connect diagnostics to export state.
    !-----------------------------------------------------------------

    ! Loop over all diagnostics in diagnostics list
    ThisDiagn => NULL()
    DO WHILE ( .TRUE. )

       ! IF (MAPL_am_I_Root()) WRITE(*,*) "Getting next diagnostic in list"

       ! Get next diagnostics in list. This will return the next
       ! diagnostics container that contains content to be written
       ! out on this time step.
       CALL Diagn_Get ( HcoState, EOI, ThisDiagn, FLAG, RC, COL=PS )

       ! IF (MAPL_am_I_Root()) THEN
       !    IF ( RC == HCO_SUCCESS .and. FLAG == HCO_SUCCESS ) WRITE(*,*) "Got it! Name: ", TRIM(ThisDiagn%cName)
       !    IF ( RC /= HCO_SUCCESS) WRITE(*,*) "Fail! RC: ", RC, "   Flag: ", FLAG
       ! ENDIF

       IF ( RC /= HCO_SUCCESS ) THEN
           CALL HCO_ERROR( 'ERROR 0', RC, THISLOC=LOC )
           RETURN
       ENDIF
       IF ( FLAG /= HCO_SUCCESS ) EXIT

       ! Only write diagnostics if this is the first Diagn_Get call for
       ! this container and time step.
       IF ( PRESENT(OnlyIfFirst) ) THEN
          IF ( OnlyIfFirst .AND. ThisDiagn%nnGetCalls > 1 ) CYCLE
       ENDIF

       ! Get pointer to ESMF EXPORT field and pass data to it (if found):

       ! 2D...
       IF ( ThisDiagn%SpaceDim == 2 ) THEN
          CALL MAPL_GetPointer ( HcoState%EXPORT, Ptr2D, &
             TRIM(ThisDiagn%cName), NotFoundOk=.TRUE., RC=STAT )
          IF ( ASSOCIATED(Ptr2D) ) THEN
             IF ( ASSOCIATED(ThisDiagn%Arr2D) ) THEN
                Ptr2D = ThisDiagn%Arr2D%Val
                !Ptr2D => ThisDiagn%Arr2D%Val
             ENDIF
          ENDIF

       ! ... or 3D
       ELSEIF ( ThisDiagn%SpaceDim == 3 ) THEN
          CALL MAPL_GetPointer ( HcoState%EXPORT, Ptr3D, &
             TRIM(ThisDiagn%cName), NotFoundOk=.TRUE., RC=STAT )
          IF ( ASSOCIATED(Ptr3D) ) THEN
             IF ( ASSOCIATED(ThisDiagn%Arr3D) ) THEN
                Ptr3D(:,:,:) = ThisDiagn%Arr3D%Val(:,:,HcoState%NZ:1:-1)
                !Ptr3D => ThisDiagn%Arr3D%Val
             ENDIF
          ENDIF
       ENDIF

       ! Free pointer
       Ptr2D => NULL()
       Ptr3D => NULL()
    ENDDO

    ! Cleanup
    ThisDiagn => NULL()

    ! Return
    RC = HCO_SUCCESS

  END SUBROUTINE HCOIO_Write
!EOC
END MODULE HCOIO_Write_Mod
#endif
