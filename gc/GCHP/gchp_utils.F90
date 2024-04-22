!------------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 910.1 and      !
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: GCHP_Utils
!
! !DESCRIPTION: Utility module for the ESMF interface to GEOS-Chem.
!\\
!\\
! !INTERFACE:
!
MODULE GCHP_Utils
!
! !USES:
!
  IMPLICIT NONE

  PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!
  PUBLIC :: Set_Background_Conc
!
! !REVISION HISTORY:
!  09 Oct 2012 - M. Long     - Initial version
!  09 Oct 2012 - R. Yantosca - Added ProTeX headers
!  09 Oct 2012 - R. Yantosca - Use F90 free-format indenting (Emacs F90 mode)
!  22 Oct 2012 - R. Yantosca - Renamed to gigc_chem_utils.F90

!EOP
!------------------------------------------------------------------------------
!BOC
  CONTAINS
!EOC
  SUBROUTINE SET_BACKGROUND_CONC( am_I_Root, SpcInfo, State_Chm, State_Met, State_Grid, Input_Opt, IND, RC)

    USE Species_Mod,      ONLY : Species
    USE State_Chm_Mod,    ONLY : ChmState
    USE State_Met_Mod,    ONLY : MetState
    USE State_Grid_Mod,   ONLY : GrdState
    USE Input_Opt_Mod,    ONLY : OptInput
    USE ErrCode_Mod
    USE Precision_Mod

    TYPE(Species),    POINTER :: SpcInfo
    TYPE(ChmState)            :: State_Chm
    TYPE(MetState)            :: State_Met
    TYPE(GrdState)            :: State_Grid
    TYPE(OptInput)            :: Input_Opt
    INTEGER, INTENT(IN)       :: IND
    INTEGER, INTENT(OUT)      :: RC
    LOGICAL, INTENT(IN)       :: am_I_Root

    INTEGER                   :: I,J,L,N

    ! Assume success
    RC        = GC_SUCCESS

    DO L = 1, State_Grid%NZ 
    DO J = 1, State_Grid%NY
    DO I = 1, State_Grid%NX
       ! Special handling for MOH (mimicking GEOS-Chem Classic)
       IF ( TRIM( SpcInfo%Name ) == 'MOH' ) THEN
          ! Test for altitude (L < 9 is always in the trop)
          IF ( L <= 9 ) THEN
             ! Test for ocean/land boxes
             IF ( State_Met%FRCLND(I,J) >= 0.5 ) THEN
                ! Continental boundary layer: 2 ppbv MOH
                State_Chm%Species(IND)%Conc(I,J,L) = 2.000e-9_fp
             ELSE
                ! Marine boundary layer: 0.9 ppbv MOH
                State_Chm%Species(IND)%Conc(I,J,L) = 0.900e-9_fp
             ENDIF
          ELSE
             ! Test for troposphere
             IF ( State_Met%InTroposphere(I,J,L) ) THEN
                ! Free troposphere: 0.6 ppbv MOH
                State_Chm%Species(IND)%Conc(I,J,L) = 0.600e-9_fp
             ELSE
                ! Strat/mesosphere:
                State_Chm%Species(IND)%Conc(I,J,L) = 1.0E-30_FP
             ENDIF
          ENDIF
       ELSEIF ( L > State_Grid%MaxChemLev .AND. &
                ( .NOT. SpcInfo%Is_Advected ) ) THEN
          ! For non-advected spc at L > State_Grid%MaxChemLev, use small number
          State_Chm%Species(IND)%Conc(I,J,L) = 1.0E-30_fp
       ELSE
          ! For all other cases, use the background value in spc db
          State_Chm%Species(IND)%Conc(I,J,L) = SpcInfo%BackgroundVV
       ENDIF
    ENDDO
    ENDDO
    ENDDO

  END SUBROUTINE SET_BACKGROUND_CONC

END MODULE GCHP_Utils
