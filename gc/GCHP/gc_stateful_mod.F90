!------------------------------------------------------------------------------
!                             The WRF-GCHP Project                            !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: GC_Stateful_Mod.F
!
! !DESCRIPTION: Module GIGC\_Stateful\_Mod stores derived type objects related
!  to Grid-Independent GEOS-Chem (GIGC)'s State structures for a coupled multi-domain
!  model. (hplin, 6/4/18)
!\\
!\\
! !INTERFACE:
module GC_Stateful_Mod
!
! !USES:
!
    use ErrCode_Mod
    use Input_Opt_Mod
    use State_Met_Mod
    use State_Chm_Mod
    use State_Diag_Mod
    use State_Grid_Mod
    use DiagList_Mod, only: DgnList
    use HCO_TYPES_MOD, only: ConfigObj
    use HCO_State_Mod, only: HCO_State
    use HCOX_State_Mod, only: Ext_State

    implicit none
    private

!
! !PUBLIC MEMBER FUNCTIONS:
!
    public :: GIGC_State_Boot
!
! !PUBLIC TYPES:
!
    ! Derived type for domain properties
    type, public :: GIGC_Stateful_Object
        integer                                        :: ID = -999
        logical                                        :: Init = .false.
        type(MetState)                                 :: State_Met
        type(ChmState)                                 :: State_Chm
        type(DgnState)                                 :: State_Diag
        type(GrdState)                                 :: State_Grid
        type(HCO_State), pointer                       :: HcoState
        type(Ext_State), pointer                       :: ExtState
    end type GIGC_Stateful_Object

    ! Global options objects
    type(OptInput), public                             :: Global_Input_Opt
    type(ConfigObj), pointer, public                   :: Global_HcoConfig => NULL()
    type(DgnList), public                              :: Global_DiagList

    ! Stateful objects
#if defined ( EXTERNAL_GRID ) || defined( MODEL_ )
    !-----------------------------------------------------------------
    !         %%%%%%% GEOS-Chem HP (with ESMF & MPI) %%%%%%%
    !
    ! The max number of domain configurations can be flexible when coupling
    ! with an external model. (hplin, 6/7/18)
    !-----------------------------------------------------------------
    integer                                            :: EXTERNAL_MAX_DOM = 8
    type(GIGC_Stateful_Object), dimension(1:8), public :: GIGC_States
#else
    !-----------------------------------------------------------------
    !         %%%%%%% GEOS-Chem CLASSIC (with OpenMP) %%%%%%%
    !
    ! For GEOS-Chem "Classic", the world is seen as a contiguous set of
    ! columns with no domain extensions possible.
    !-----------------------------------------------------------------
    integer                                            :: EXTERNAL_MAX_DOM = 1
    type(GIGC_Stateful_Object), dimension(1:1), public :: GIGC_States
#endif
    logical                                            :: Init = .false.

!
! !REMARKS:
!  This module is intended to be the ONLY stateful module in WRF-GC, meaning that all other
!  modules should not be statefully programmed (in an ideal world, of course) & only rely on
!  this module to retrieve/set derived type objects with storage within CPUs.
!
!  Note: Fortran passes derived type objects by intrinsic copy, not by reference.
!  This means that there are some maneuverings required in your coding, specifically relating to
!  how to work with Input_Opt, a global variable. Namely, the workflow in external models would be
!  (EM - External Model, CH - GIGC_CHUNK_MOD, ST - GC_Stateful_Mod)
!
!  > On Initialization <
!   EM - check if initialized GIGC_State, if not
!         ST - call GIGC_State_Boot, initializes Global_Input_Opt, Global_HcoConfig
!      - call GIGC_Chunk_Init, do all the chunk initialization routines
!         CH - operate on datas and return
!      - call GIGC_State_Init, ... to store data into appropriate containers by ID - opens a new container
!      - if needing other operations, call GIGC_State_Get_... to get initialized datas or _Set_ to set.
!
!  > On Column-Code Running <
!   EM - call GIGC_State_Get_... to get initialized datas to pass into GIGC_CHUNK_MOD
!         CH - operate on datas and return
!   EM - call GIGC_State_Set_Opt, Met, Chm ... to store data into appropriate containers by ID
!      - operate as necessary
!      - call GIGC_State_Set_Opt, Met, Chm ... to store data if operated on, otherwise skip
!
!  Note the drastic workflow change in GIGC_Chunk_Init: This means that the part where options are read
!  is NO LONGER in GIGC_Chunk_Init but part of an upper model. This facilitates future I/O changes, as
!  I/O should be really handled by a parent model or a gridded component anyway.
!
! !REVISION HISTORY:
!  04 Jun 2018 - H.P. Lin  - First crack.
!  07 Jun 2018 - H.P. Lin  - Added workflow description and clarified code.
!  28 Dec 2019 - H.P. Lin  - Update to GEOS-Chem 12.6.3
!  19 May 2020 - H.P. Lin  - Booyah!
!------------------------------------------------------------------------------
!BOC
contains
!EOC
!------------------------------------------------------------------------------
!                             The WRF-GCHP Project                            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: GIGC_State_Boot
!
! !DESCRIPTION: Subroutine GIGC\_State\_Boot initializes the stateful module and
!  read configuration variables to store in Input_Opt, HcoConfig. (hplin, 6/7/18)
!\\
!\\
! !INTERFACE:
!
    subroutine GIGC_State_Boot(am_I_Root,  MPI_COMM, &
                               State_Grid, RC)
!
! !USES:
! 
        USE GC_Environment_Mod, only: GC_Allocate_All
        USE INPUT_MOD,          only: Read_Input_File
        USE HCO_CONFIG_MOD,     only: Config_Readfile
        USE HCO_DIAGN_MOD,      only: DiagnFileOpen
        USE LINOZ_MOD,          only: Linoz_Read
        USE DiagList_Mod,       only: Init_DiagList, Print_DiagList

!
! !INPUT PARAMETERS:
!
        logical, intent(in)           :: am_I_Root          ! Are we on the root CPU?
        integer, intent(in)           :: MPI_COMM           ! MPI Communicator #
!
! !INPUT/OUTPUT PARAMETERS:
!
        type(GrdState), intent(inout) :: State_Grid         ! Grid State object
        integer, intent(inout)        :: RC                 ! Success or failure?
!
! !REMARKS:
!  This code was written with the intent to supercede GIGC_Get_Options, with
!  the role of reading configuration files & initializing global vars such as
!  Input_Opt, HcoConfig, ... into Stateful_Mod. Note that in this design, you
!  cannot change "input.geos"/"HEMCO_Config.rc" across domains, so they are still
!  not independent.
!
!  It is OK to call this method more than once - it will not crash and instead work
!  gracefully.
!
! !REVISION HISTORY:
!  07 Jun 2018 - H.P. Lin   - Initial version.
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
        integer :: HCO_DIAGN_LUN

        ! Assume success
        RC = GC_SUCCESS

        ! If we are already initialized, return
#if defined ( LINUX_GFORTRAN )
        if(Init .eqv. .true.) then
#else
        if(Init .eq. .true.) then
#endif
            return
        endif

        ! Initialize Input_Opt fields to zeros or equivalent (v12)
        call Set_Input_Opt( am_I_Root, Global_Input_Opt, RC )

        ! Some necessary set-up for GEOS-Chem HP
#if defined ( EXTERNAL_GRID ) || defined( MODEL_ )
        !-----------------------------------------------------------------
        !         %%%%%%% GEOS-Chem HP (with ESMF & MPI) %%%%%%%
        !
        ! We have to set up some HPC and Root CPU parameters into Input_Opt
        ! for GEOS-Chem HP / GIGC. This must be done after GC_Allocate_All
        ! so it is not overwritten by GC_Allocate_All. (hplin, 6/13/18)
        !-----------------------------------------------------------------
        Global_Input_Opt%isMPI   = .true.
        Global_Input_Opt%AmIRoot = am_I_Root
        Global_Input_Opt%MPIComm = MPI_COMM

        ! Set some DEFAULT time-steps which will be overwritten by GIGC_Chunk_Mod later on
        Global_Input_Opt%TS_CHEM = 10   ! Chemistry timestep [min]
        Global_Input_Opt%TS_EMIS = 10   ! Chemistry timestep [min]
        Global_Input_Opt%TS_DYN  = 20   ! Dynamic   timestep [min]
        Global_Input_Opt%TS_CONV = 20   ! Dynamic   timestep [min]
#endif

        ! Read input.geos, now done on all threads since GC(HP) 12.2.0
        ! to remove GIGC MPI wrapper dependency
        ! Read input.geos at very beginning of simulation on every thread
        CALL Read_Input_File( Global_Input_Opt, State_Grid, RC )
        if(RC /= GC_SUCCESS) then 
            write(6, *) "STOP GC_Stateful_Mod :: Return Code /= GC_SUCCESS (Read_Input_File)"
            return
        endif
        ! Echo info
        write(6, *) '### GC_Stateful_Mod: Root CPU, after READ_INPUT_FILE'

        ! In the ESMF/MPI environment, we can get the total overhead ozone
        ! either from the met fields (GIGCsa) or from the Import State (GEOS-5)
        Global_Input_Opt%USE_O3_FROM_MET = .TRUE.

        ! Read LINOZ climatology
        ! This is here because it is read into Input_Opt (hplin, 12/28/19)
        IF ( Global_Input_Opt%LLINOZ ) THEN
            CALL Linoz_Read( Global_Input_Opt, RC ) 
            if(RC /= GC_SUCCESS) then 
                write(6, *) "STOP GC_Stateful_Mod :: Return Code /= GC_SUCCESS (Linoz_Read)"
                return
            endif

            ! Echo info
            write(6, *) '### GC_Stateful_Mod: Root CPU, after LINOZ_READ'
        ENDIF

        ! Read HEMCO Configuration File
        ! No longer required as EMISSIONS_INIT will handle it (for now) and pass back later.
        ! (hplin, 8/6/18)
        ! call Config_ReadFile(am_I_Root, Global_HcoConfig, Global_Input_Opt%HcoConfigFile, 0, RC)

        ! Read HEMCO Diagnostic File
        ! call DiagnFileOpen(am_I_Root, Global_HcoConfig, HCO_DIAGN_LUN, RC)

        ! Read HISTORY.rc & initialize diagnostics list object
        ! Ported from gigc_historyexports_mod, v12
        CALL Init_DiagList(am_I_Root, "HISTORY.rc", Global_DiagList, RC)
        if(RC /= GC_SUCCESS) then 
            write(6, *) "STOP GC_Stateful_Mod :: Return Code /= GC_SUCCESS (Init_DiagList)"
            return
        endif

        CALL Print_DiagList(am_I_Root, Global_DiagList, RC)

        !
        ! Allocate GEOS-Chem module arrays,
        !
        ! GC_Allocate_All's task is to initialize the CMN_ module arrays.
        ! Input_Opt is no longer initialized as of GEOS-Chem v12, so this should be called
        ! AFTER Input_Opt is read, as some options depend on Input_Opt%LUCX, for example.
        ! (hplin, 8/7/18)
        !
        ! Generic sizings are passed in -- they will be switched by GIGC_Switch_Dims at runtime.
        !
        ! Allocate all lat/lon arrays
        CALL GC_Allocate_All( Global_Input_Opt, State_Grid, RC )
        if(RC /= GC_SUCCESS) then 
            write(6, *) "STOP GC_Stateful_Mod :: Return Code /= GC_SUCCESS (GC_Allocate_All)"
            return
        endif

        ! We are finished!
        if(RC /= GC_SUCCESS) then 
            write(6, *) "STOP GC_Stateful_Mod :: Return Code /= GC_SUCCESS"
        endif
        
        ! We are now initialized!
        Init = .true.
    
    end subroutine GIGC_State_Boot
!EOC
end module GC_Stateful_Mod