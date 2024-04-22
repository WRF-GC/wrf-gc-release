!------------------------------------------------------------------------------
!                   Harmonized Emissions Component (HEMCO)                    !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: hcox_custom_mod.F90
!
! !DESCRIPTION: Customizable HEMCO emission extension.
!\\
!\\
! !INTERFACE:
!
MODULE HCOX_Custom_Mod
!
! !USES:
!
  USE HCO_Error_MOD
  USE HCO_Diagn_MOD
  USE HCOX_State_MOD, ONLY : Ext_State
  USE HCO_State_MOD,  ONLY : HCO_State

  IMPLICIT NONE
  PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!
  PUBLIC :: HCOX_Custom_Run
  PUBLIC :: HCOX_Custom_Init
  PUBLIC :: HCOX_Custom_Final
!
! !REVISION HISTORY:
!  13 Dec 2013 - C. Keller   - Initial version
!  See https://github.com/geoschem/hemco for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !MODULE VARIABLES:
!
  TYPE :: MyInst
   INTEGER                         :: Instance
   INTEGER                         :: ExtNr   = -1
   INTEGER                         :: nOcWind = -1
   INTEGER                         :: nIceSrc = -1
   INTEGER,      POINTER           :: OcWindIDs(:)
   INTEGER,      POINTER           :: IceSrcIDs(:)
   TYPE(MyInst), POINTER           :: NextInst => NULL()
  END TYPE MyInst

  ! Pointer to instances
  TYPE(MyInst), POINTER            :: AllInst => NULL()

CONTAINS
!EOC
!------------------------------------------------------------------------------
!                   Harmonized Emissions Component (HEMCO)                    !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: HCOX_Custom_Run
!
! !DESCRIPTION: Subroutine HCOX\_Custom\_Run is the driver routine
! for the customizable HEMCO extension.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE HCOX_Custom_Run( ExtState, HcoState, RC )
!
! !USES:
!
    USE HCO_FluxArr_Mod,  ONLY : HCO_EmisAdd
    USE HCO_GeoTools_Mod, ONLY : HCO_LANDTYPE
!
! !INPUT PARAMETERS:
!
    TYPE(Ext_State), POINTER       :: ExtState    ! Module options
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(HCO_State), POINTER       :: HcoState    ! Hemco state
    INTEGER,         INTENT(INOUT) :: RC          ! Success or failure
!
! !REMARKS:
!
!
! !REVISION HISTORY:
!  13 Dec 2013 - C. Keller   - Initial version
!  See https://github.com/geoschem/hemco for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER               :: I, J, N, AS, LANDTYPE
    INTEGER               :: tmpID
    REAL*8                :: W10M
    REAL(hp), ALLOCATABLE :: FLUXICE(:,:)
    REAL(hp), ALLOCATABLE :: FLUXWIND(:,:)
    LOGICAL               :: ERR
    CHARACTER(LEN=255)    :: MSG, LOC

    TYPE(MyInst), POINTER :: Inst
!
! !DEFINED PARAMETERS:
!
    REAL*8,   PARAMETER :: SCALICE  = 1.0d-14
    REAL*8,   PARAMETER :: SCALWIND = 1.0d-14

    !=================================================================
    ! HCOX_CUSTOM_RUN begins here!
    !=================================================================
    LOC = 'HCOX_CUSTOM_RUN (HCOX_CUSTOM_MOD.F90)'

    ! Enter
    CALL HCO_ENTER( HcoState%Config%Err, LOC, RC )
    IF ( RC /= HCO_SUCCESS ) THEN
        CALL HCO_ERROR( 'ERROR 0', RC, THISLOC=LOC )
        RETURN
    ENDIF

    ! Set error flag
    ERR = .FALSE.

    ! Sanity check: return if extension not turned on
    IF ( ExtState%Custom <= 0 ) RETURN

    ! Get instance
    Inst => NULL()
    CALL InstGet ( ExtState%Custom, Inst, RC )
    IF ( RC /= HCO_SUCCESS ) THEN
       WRITE(MSG,*) 'Cannot find custom instance Nr. ', ExtState%Custom
       CALL HCO_ERROR(MSG,RC)
       RETURN
    ENDIF

    ! Initialize flux arrays
    ALLOCATE ( FLUXICE( HcoState%NX,HcoState%NY),        &
               FLUXWIND(HcoState%NX,HcoState%NY), STAT=AS )
    IF ( AS/= 0 ) THEN
       CALL HCO_ERROR( 'ALLOCATION ERROR', RC )
       RETURN
    ENDIF
    FLUXICE  = 0.0_hp
    FLUXWIND = 0.0_hp

!$OMP PARALLEL DO                                            &
!$OMP DEFAULT( SHARED )                                      &
!$OMP PRIVATE( I, J, W10M, LANDTYPE                        ) &
!$OMP SCHEDULE( DYNAMIC )
    ! Loop over surface grid boxes
    DO J = 1, HcoState%NY
    DO I = 1, HcoState%NX

       ! Get the land type for grid box (I,J)
       LANDTYPE = HCO_LANDTYPE( ExtState%FRLAND%Arr%Val(I,J),   &
                                ExtState%FRLANDIC%Arr%Val(I,J), &
                                ExtState%FROCEAN%Arr%Val(I,J),  &
                                ExtState%FRSEAICE%Arr%Val(I,J), &
                                ExtState%FRLAKE%Arr%Val(I,J)   )

       ! Check surface type
       ! Ocean:
       IF ( LANDTYPE == 0 ) THEN

          ! 10m wind speed [m/s]
          W10M = ExtState%U10M%Arr%Val(I,J)**2 + &
                 ExtState%V10M%Arr%Val(I,J)**2
          W10M = SQRT(W10M)

          ! Set flux to wind speed
          FLUXWIND(I,J) = W10M * SCALWIND

       ! Ice:
       ELSE IF ( LANDTYPE == 2 ) THEN

          ! Set uniform flux
          FLUXICE(I,J) = SCALICE
       ENDIF

    ENDDO !I
    ENDDO !J
!$OMP END PARALLEL DO

    ! Check exit status
    IF ( ERR ) THEN
       RC = HCO_FAIL
       RETURN
    ENDIF

    ! Add wind fluxes to emission arrays & diagnostics
    DO N = 1, Inst%nOcWind

       ! Emissions array
       CALL HCO_EmisAdd( HcoState, FLUXWIND, Inst%OcWindIDs(N), &
                         RC,       ExtNr=Inst%ExtNr )
       IF ( RC /= HCO_SUCCESS ) THEN
           CALL HCO_ERROR( 'ERROR 1', RC, THISLOC=LOC )
           RETURN
       ENDIF
    ENDDO !N

    ! Add ice fluxes to emission arrays & diagnostics
    DO N = 1, Inst%nIceSrc

       ! Emissions array
       CALL HCO_EmisAdd( HcoState, FLUXICE, Inst%IceSrcIDs(N), &
                         RC,       ExtNr=Inst%ExtNr )
       IF ( RC /= HCO_SUCCESS ) THEN
           CALL HCO_ERROR( 'ERROR 2', RC, THISLOC=LOC )
           RETURN
       ENDIF
    ENDDO !N

    ! Return w/ success
    Inst => NULL()
    CALL HCO_LEAVE( HcoState%Config%Err,RC )

  END SUBROUTINE HCOX_Custom_Run
!EOC
!------------------------------------------------------------------------------
!                   Harmonized Emissions Component (HEMCO)                    !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: HCOX_Custom_Init
!
! !DESCRIPTION: Subroutine HCOX\_Custom\_Init initializes the HEMCO
! CUSTOM extension.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE HCOX_Custom_Init( HcoState, ExtName, ExtState, RC )
!
! !USES:
!
    USE HCO_ExtList_Mod,    ONLY : GetExtNr
    USE HCO_STATE_MOD,      ONLY : HCO_GetExtHcoID
!
! !INPUT PARAMETERS:
!
    CHARACTER(LEN=*), INTENT(IN   ) :: ExtName    ! Extension name
    TYPE(Ext_State),  POINTER       :: ExtState   ! Module options
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(HCO_State),  POINTER       :: HcoState   ! Hemco state
    INTEGER,          INTENT(INOUT) :: RC

! !REVISION HISTORY:
!  13 Dec 2013 - C. Keller   - Now a HEMCO extension
!  See https://github.com/geoschem/hemco for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER                        :: ExtNr, N, nSpc, AS
    INTEGER,           ALLOCATABLE :: HcoIDs(:)
    LOGICAL                        :: verb
    CHARACTER(LEN=31), ALLOCATABLE :: SpcNames(:)
    CHARACTER(LEN=255)             :: MSG, LOC
    TYPE(MyInst), POINTER          :: Inst

    !=================================================================
    ! HCOX_CUSTOM_INIT begins here!
    !=================================================================
    LOC = 'HCOX_CUSTOM_INIT (HCOX_CUSTOM_MOD.F90)'

    ! Extension Nr.
    ExtNr = GetExtNr( HcoState%Config%ExtList, TRIM(ExtName) )
    IF ( ExtNr <= 0 ) RETURN

    ! Enter
    CALL HCO_ENTER( HcoState%Config%Err, LOC, RC )
    IF ( RC /= HCO_SUCCESS ) THEN
        CALL HCO_ERROR( 'ERROR 3', RC, THISLOC=LOC )
        RETURN
    ENDIF
    verb = HCO_IsVerb(HcoState%Config%Err,1)

    Inst => NULL()
    CALL InstCreate ( ExtNr, ExtState%Custom, Inst, RC )
    IF ( RC /= HCO_SUCCESS ) THEN
       CALL HCO_ERROR (  'Cannot create custom instance', RC )
       RETURN
    ENDIF

    ! Set species IDs
    CALL HCO_GetExtHcoID( HcoState, Inst%ExtNr, HcoIDs, SpcNames, nSpc, RC )
    IF ( RC /= HCO_SUCCESS ) THEN
        CALL HCO_ERROR( 'ERROR 4', RC, THISLOC=LOC )
        RETURN
    ENDIF

    ! Assume first half are 'wind species', second half are ice.
    IF ( MOD(nSpc,2) /= 0 ) THEN
       MSG = 'Cannot set species IDs for custom emission module!'
       CALL HCO_ERROR(MSG, RC )
       RETURN
    ENDIF

    ! Pass # of sources
    Inst%nOcWind = nSpc / 2
    Inst%nIceSrc = nSpc / 2

    ! Allocate vector w/ the species IDs
    ALLOCATE ( Inst%OcWindIDs(Inst%nOcWind) )
    ALLOCATE ( Inst%IceSrcIDs(Inst%nIceSrc) )
    Inst%OcWindIDs(:) = HcoIDs(1:Inst%nOcWind)
    N = Inst%nOcWind + 1
    Inst%IceSrcIDs(:) = HcoIDs(N:nSpc)

    ! Verbose mode
    IF ( verb ) THEN
       MSG = 'Use custom emissions module (extension module)'
       CALL HCO_MSG(HcoState%Config%Err,MSG )

       MSG = 'Use the following species (Name: HcoID):'
       CALL HCO_MSG(HcoState%Config%Err,MSG)
       DO N = 1, nSpc
          WRITE(MSG,*) TRIM(SpcNames(N)), ':', HcoIDs(N)
          CALL HCO_MSG(HcoState%Config%Err,MSG)
       ENDDO
    ENDIF

    ! Activate met fields required by this extension
    ExtState%U10M%DoUse = .TRUE.
    ExtState%V10M%DoUse = .TRUE.
    ExtState%FRLAND%DoUse   = .TRUE.
    ExtState%FRLANDIC%DoUse = .TRUE.
    ExtState%FROCEAN%DoUse  = .TRUE.
    ExtState%FRSEAICE%DoUse = .TRUE.
    ExtState%FRLAKE%DoUse   = .TRUE.

    ! Activate this extension
    !ExtState%Custom = .TRUE.

    ! Leave w/ success
    IF ( ALLOCATED(HcoIDs  ) ) DEALLOCATE(HcoIDs  )
    IF ( ALLOCATED(SpcNames) ) DEALLOCATE(SpcNames)

    CALL HCO_LEAVE( HcoState%Config%Err,RC )

  END SUBROUTINE HCOX_Custom_Init
!EOC
!------------------------------------------------------------------------------
!                   Harmonized Emissions Component (HEMCO)                    !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: HCOX_Custom_Final
!
! !DESCRIPTION: Subroutine HCOX\_Custom\_Final finalizes the HEMCO
! CUSTOM extension.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE HCOX_Custom_Final ( ExtState )
!
! !INPUT PARAMETERS:
!
    TYPE(Ext_State),  POINTER       :: ExtState   ! Module options
!
! !REVISION HISTORY:
!  13 Dec 2013 - C. Keller   - Now a HEMCO extension
!  See https://github.com/geoschem/hemco for complete history
!EOP
!------------------------------------------------------------------------------
!BOC

    !=================================================================
    ! HCOX_CUSTOM_FINAL begins here!
    !=================================================================
    CALL InstRemove ( ExtState%Custom )

  END SUBROUTINE HCOX_Custom_Final
!EOC
!------------------------------------------------------------------------------
!                   Harmonized Emissions Component (HEMCO)                    !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: InstGet
!
! !DESCRIPTION: Subroutine InstGet returns a poiner to the desired instance.
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
! !DESCRIPTION: Subroutine InstCreate creates a new instance.
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

    ! ----------------------------------------------------------------
    ! Type specific initialization statements follow below
    ! ----------------------------------------------------------------

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
! !DESCRIPTION: Subroutine InstRemove creates a new instance.
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
    TYPE(MyInst), POINTER       :: PrevInst
    TYPE(MyInst), POINTER       :: Inst

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
       IF ( ASSOCIATED( Inst%OcWindIDs ) ) THEN
          DEALLOCATE ( Inst%OcWindIDs )
       ENDIF
       Inst%OcWindIDs => NULL()

       IF ( ASSOCIATED( Inst%IceSrcIDs ) ) THEN
          DEALLOCATE ( Inst%IceSrcIDs )
       ENDIF
       Inst%IceSrcIDs => NULL()

       !---------------------------------------------------------------------
       ! Pop off instance from list
       !---------------------------------------------------------------------
       IF ( ASSOCIATED(PrevInst) ) THEN
          PrevInst%NextInst => Inst%NextInst
       ELSE
          AllInst => Inst%NextInst
       ENDIF
       DEALLOCATE(Inst)
       Inst => NULL()
    ENDIF

    ! Free pointers before exiting
    PrevInst => NULL()
    Inst     => NULL()

   END SUBROUTINE InstRemove
!EOC
END MODULE HCOX_Custom_Mod
