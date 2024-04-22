!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: gchp_chunk_mod
!
! !DESCRIPTION: Module GC\_CHUNK\_MOD is the module that contains the init,
!  and run methods for the ESMF interface to GEOS-Chem.
!\\
!\\
! !INTERFACE:
!
MODULE GIGC_Chunk_Mod
!
! !USES:
!      
  USE GCHP_Utils

  USE ErrCode_Mod
  USE Precision_Mod

  IMPLICIT NONE
  PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!
  PUBLIC :: GIGC_Chunk_Init
  PUBLIC :: GIGC_Chunk_Run
!
! !PRIVATE MEMBER FUNCTIONS:
!
  PRIVATE  ::  GIGC_Switch_Dims
  INTEGER  ::  MemDebugLevel
!
! !PUBLIC TYPES: (WRF-GC)
!
  ! Derived type for chunk operator options (for GIGC_Chunk_Run), replacing previous
  ! phases.
  TYPE GIGC_Chunk_Operators
    logical                     :: Conv
    logical                     :: DryDep
    logical                     :: Emis
    logical                     :: Tend
    logical                     :: Turb
    logical                     :: Chem
    logical                     :: WetDep
    logical                     :: Rad

    logical                     :: GCDiagn      ! Use GEOS-Chem History-based GCC/NetCDF Diagns?
  END TYPE GIGC_Chunk_Operators
  PUBLIC :: GIGC_Chunk_Operators

  INTEGER, PARAMETER :: KIND_R4 = selected_real_kind(3,25)

!
! !REVISION HISTORY:
!  22 Jun 2009 - R. Yantosca & P. Le Sager - Chunkized & cleaned up.
!  28 Dec 2019 - H.P. Lin    - Build version for WRF-GC v1.1 (12.6.3)
!  22 May 2020 - H.P. Lin    - Build version for WRF-GC Tech 2020 (12.8.1)
!  18 Apr 2022 - H.P. Lin    - Build version for WRF-GC Tech 2022 (13.4.0)
!EOP
!------------------------------------------------------------------------------
!BOC

CONTAINS
!EOC
!------------------------------------------------------------------------------
!                             The WRF-GCHP Project                            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: GIGC_Switch_Dims
!
! !DESCRIPTION: Subroutine GIGC\_Switch\_Dims contains a set of very, very hackish
!  operations to attempt to "manage" the in-module "state" variables that should've
!  been stored inside a derived type object. This is a temporary shim until we can
!  coordinate with GCST regarding moving everything into derived type objects.
!  (This is also a good reference of which variables are stateful) (hplin, 6/4/18)
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE GIGC_Switch_Dims( ID,                                          &
                               lonCtr,     latCtr,      lonEdge,    latEdge,&
                               Input_Opt,                                   &
                               State_Met,  State_Chm,   State_Diag,         &
                               State_Grid,                                  &
                               CMN_Alloc,  GC_Alloc,                        &
                               RC )
!
! !USES:
!
    USE ErrCode_Mod

    ! USE HCOI_GC_Main_Mod, ONLY : HCOI_GC_Init, HCOI_GC_Final

    ! GEOS-Chem Modules: COMMON
    ! Manually re-allocated
    USE CMN_FJX_MOD,        ONLY : IRHARR, ZPJ, ODMDUST, ODAER, ISOPOD
    USE CMN_FJX_MOD,        ONLY : JVN_, NWVAA, NDUST, NAER

    ! GEOS-Chem Modules: To fix soon into State_Chm
    ! Auto re-allocated
    USE GC_Environment_Mod, ONLY : GC_Init_Grid
    USE AEROSOL_MOD,        ONLY : INIT_AEROSOL, CLEANUP_AEROSOL
    USE CARBON_MOD,         ONLY : INIT_CARBON, CLEANUP_CARBON
    USE PRESSURE_MOD,       ONLY : CLEANUP_PRESSURE, INIT_PRESSURE
    USE PRESSURE_MOD,       ONLY : GET_AP, GET_BP, Accept_External_ApBp
    USE Linear_CHEM_MOD,     ONLY : CLEANUP_Linear_CHEM, INIT_Linear_CHEM
    USE SULFATE_MOD,        ONLY : INIT_SULFATE, CLEANUP_SULFATE

    USE CO2_MOD,            ONLY : INIT_CO2, CLEANUP_CO2

    USE Vdiff_Mod,          ONLY : Cleanup_Vdiff, Init_Vdiff

    USE GC_Grid_Mod,        ONLY : SetGridFromCtrEdges

    ! HEMCO Modules
    USE HCO_Error_Mod,      ONLY : hp
    USE HCO_State_GC_Mod,   ONLY : HcoState, ExtState

    USE Input_Opt_Mod,      ONLY : OptInput
    USE State_Met_Mod,      ONLY : MetState
    USE State_Chm_Mod,      ONLY : ChmState
    USE State_Diag_Mod,     ONLY : DgnState
    USE State_Grid_Mod,     ONLY : GrdState, Cleanup_State_Grid

    USE GC_Stateful_Mod
!
! !INPUT PARAMETERS:
!
    INTEGER,            INTENT(IN)    :: ID           ! Domain identifier within this CPU
    REAL(KIND_R4),      INTENT(IN)    :: lonCtr(:,:)  ! Lon centers [radians]
    REAL(KIND_R4),      INTENT(IN)    :: latCtr(:,:)  ! Lat centers [radians]
    REAL(KIND_R4),      INTENT(IN)    :: lonEdge(:,:) ! Lat centers [radians]
    REAL(KIND_R4),      INTENT(IN)    :: latEdge(:,:) ! Lat centers [radians]
    LOGICAL,            INTENT(IN)    :: CMN_Alloc    ! Do CMN allocations? If you skip GC_Allocate_All, use this
    LOGICAL,            INTENT(IN)    :: GC_Alloc     ! Do other allocations? Not in init

!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(OptInput),     INTENT(INOUT) :: Input_Opt
    TYPE(MetState),     INTENT(INOUT) :: State_Met
    TYPE(ChmState),     INTENT(INOUT) :: State_Chm
    TYPE(DgnState),     INTENT(INOUT) :: State_Diag
    TYPE(GrdState),     INTENT(INOUT) :: State_Grid
!
! !OUTPUT PARAMETERS:
!
    INTEGER,            INTENT(OUT)   :: RC
!
! !REMARKS:
!  This code is basically a set of hacks originally called GIGC\_Do\_Nested\_Hacks.
!  The name says it all - it attempts to "patch" the lack of derived type objects at a 
!  higher level than you're supposed to, but it enables in-CPU switching of array dimensions
!  without resorting to changing core GEOS-CHEM code. This should be later coordinated with
!  GCST to eventually remove all this using the "hacks" here as a reference to where should be
!  changed.
!
!  CMN_Alloc: Only allocate bare essentials and deallocate all, as there will be
!  a pass of GC_Allocate_All and GC_Init_Extra later.
!
!  GC_Alloc: Deallocate AND reallocate for domain switching.
!
! !REVISION HISTORY:
!  04 Jun 2018 - H.P. Lin   - First crack.
!  19 May 2020 - H.P. Lin   - Massively improved implementation.
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER, SAVE                     :: PREVIOUS_ID = -1
    INTEGER, SAVE                     :: PREVIOUS_GC_ID = -1
    LOGICAL                           :: HCO_ERROR
    LOGICAL                           :: INIT_ID

    REAL(fp)                          :: Ap(State_Grid%NZ+1), Bp(State_Grid%NZ+1)
    REAL(hp)                          :: Ap_Pa(State_Grid%NZ+1)
    INTEGER                           :: I, J, L
    INTEGER                           :: IM, JM, LM

    INIT_ID = GIGC_States(ID)%Init

    ! All other state objects received here are already in the correct configuration,
    ! because they have been retrieved from the GEOS-Chem state container.
    !
    ! Only State_Grid needs to be allocated in this case - so check INIT_ID (initialization status)
    ! to see if we need to re-allocate, later in the code.

    ! Get IM, JM, LM. State_Grid has already been modified at this point
    ! from WRFGC_Get_WRF.
    IM = State_Grid%NX
    JM = State_Grid%NY
    LM = State_Grid%NZ

    ! If we are already in this domain we do not need to perform switching...
    if(PREVIOUS_GC_ID .eq. ID) then
        write(6, *) "Switch_Dims: Switching is not necessary as we are already"
        write(6, *) "Switch_Dims: in the right domain!"

        return
    endif

    write(6, *) "%%%%%%%%%% WRF-GC Domain Switcher (Debug) %%%%%%%%%%"
    write(6, *) "    (Not very) Proudly brought to you by hplin      "
    write(6, *) "rv: 20220810 1825"
    write(6, *) "Some debug information:"
    write(6, *) "ID:", ID
    write(6, *) "IM, JM, LM:", IM, JM, LM
    write(6, *) "Allocs: CMN", CMN_Alloc, "GC", GC_Alloc
    write(6, *) "PREVIOUS_ID: ", PREVIOUS_ID, "PREVIOUS_GC_ID", PREVIOUS_GC_ID
    write(6, *) "RC at entry-point:", RC
    write(6, *) "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"

    ! Retrieve and put the HEMCO state in the GEOS-Chem state container
    if(PREVIOUS_ID .ne. -1) then
        ! Now also update HEMCO post-initialization (hplin, 2/26/19)
        GIGC_States(PREVIOUS_ID)%HcoState => HcoState
        GIGC_States(PREVIOUS_ID)%ExtState => ExtState
        write(6,*) "Switch_Dims: Saving HcoState from ", PREVIOUS_ID
    endif

    if(PREVIOUS_ID .ne. -1 .and. PREVIOUS_ID .ne. ID) then
        ! Switch out HEMCO for re-initialization, if not initialized previously
        if(.not. INIT_ID) then
            nullify(HcoState)
            nullify(ExtState)
            write(6,*) "Switch_Dims: HcoState for ID not previously init'd, assigning", ID
        else
            ! Or just retrieve it from the Stateful container...
            HcoState => GIGC_States(ID)%HcoState
            ExtState => GIGC_States(ID)%ExtState

            write(6,*) "Switch_Dims: Switch HcoState from PREV->ID", PREVIOUS_ID, ID
            write(6,*) "Switch_Dims: Sanity check assoc", associated(HcoState), associated(ExtState)
            write(6,*) "Switch_Dims: Locations of Hco, Ext", loc(HcoState), loc(ExtState)
        endif
    else
      write(6,*) "Switch_Dims: No need to switch HcoState for ID", ID
    endif

    if(associated(HcoState)) then
      write(6,*) "Switch_Dims: HcoState nSpc, NX, NY, NZ", HcoState%nSpc, HcoState%NX, HcoState%NY, HcoState%NZ
    endif

    ! (re)Allocate State_Grid arrays
    ! Initialize GEOS-Chem horizontal grid structure
    ! Allocate State_Grid arrays is also done in GC_Init_Grid so 
    !
    ! No need to re-allocate if previously performed in chemics_init, I think
    if(.not. INIT_ID) then
        write(6,*) "Switch_Dims: Debug ID", ID, "not initialized, initializing State_Grid"

        ! Initialize GEOS-Chem horizontal grid structure
        call GC_Init_Grid( Input_Opt, State_Grid, RC )

        ! Now write necessary information back...
        State_Grid%ID = ID
        State_Grid%NX = IM
        State_Grid%NY = JM
        State_Grid%NZ = LM
        State_Grid%GlobalNX = IM
        State_Grid%GlobalNY = JM
        State_Grid%NativeNZ = LM

        write(6,*) "Switch_Dims: After GC_Init_Grid", ID, IM, JM, LM, RC
    endif

    if(.not. associated( State_Grid%Area_M2 )) then
        write(6,*) "Switch_Dims: FATAL: Area_M2 is not allocated in State_Grid!!!"
    endif

    ! This is necessary to avoid crashing the whole Area_M2 debacle
    ! probably have to merge it into a core GC code soon
    ! (Only necessary in GIGC_Chunk_Run, in Init we run before things are initialized)
    if(associated( State_Met%Area_M2 )) then
        State_Grid%Area_M2 = State_Met%Area_M2
    else
        write(6,*) "Switch_Dims: WARNING: Area_M2 is not allocated in State_Met."
    endif

    ! Set lat/lon centers and calculate grid box parameters for the GEOS-Chem Grid.
    ! This is to satisfy HEMCO geographical coordinate requirements
    ! while using GIGC_Chunk_Run as a pure 1-D model that
    ! is technically unaware of any surrounding grid boxes.
    !
    ! As a note, you must ensure that IM, LM, JM are consistent across your calls to
    ! GIGC_Chunk_Init and GIGC_Chunk_Run, as the G-C grid is initialized by
    ! GIGC_Chunk_Init -> GIGC_Init_Simulation and arrays are opened to the extent of
    ! the original IM, JM, LM (_WORLD) parameters.
    ! Failure to do so will result in array overflows and almost certainly SIGSEGV.
    !
    ! Note that we've directly inferred NX, NY (IM, JM in GCHP-speak) from lonCtr, latCtr.
    ! This serves as a consistency check by making use of GC_GRID_MOD::SetGridFromCtr's
    ! checks against SIZE(XMID, 1) & SIZE(XMID, 2) respectively. (hplin, 4/24/18)
    IF (SIZE(lonCtr, 1) .ne. SIZE(latCtr, 1)) THEN
      WRITE(6, *) "GIGC_Chunk_Mod Fatal: Size of lonCtr, latCtr arrays on I do not match"
      RC = GC_FAILURE
      RETURN
    ENDIF

    IF (SIZE(lonCtr, 2) .ne. SIZE(latCtr, 2)) THEN
      WRITE(6, *) "GIGC_Chunk_Mod Fatal: Size of lonCtr, latCtr arrays on J do not match"
      RC = GC_FAILURE
      RETURN
    ENDIF

    CALL SetGridFromCtrEdges( Input_Opt, State_Grid, lonCtr, latCtr, lonEdge, latEdge, RC )
    IF ( RC /= GC_SUCCESS ) THEN
      WRITE(6, *) "GIGC_Chunk_Mod Fatal: Could not set GC grid parameters from lat/lon ctr/edges"
      RETURN
    ENDIF

    !--------------------------------------------------------------
    ! WRF-GC 2.1 Two-domain allocations
    ! Coordinate with GCST in the future to move these to other appropriate places.
    ! (hplin, 3/19/20)
    !--------------------------------------------------------------

    if(Input_Opt%ITS_A_CO2_SIM) then
        call Cleanup_CO2()
        if(GC_Alloc) then
            call Init_CO2(Input_Opt, State_Grid, RC)
        endif
    endif

    ! CMN_FJX_MOD
    ! Only allocate those who are grid-dependent.
    if(allocated(ZPJ)) deallocate(ZPJ)
    if(allocated(ODMDUST)) deallocate(ODMDUST)
    if(allocated(ODAER)) deallocate(ODAER)
    if(allocated(ISOPOD)) deallocate(ISOPOD)
    if(allocated(IRHARR)) deallocate(IRHARR)

    allocate(ZPJ(State_Grid%NZ, JVN_, State_Grid%NX, State_Grid%NY), stat=RC)
    allocate(ODMDUST(State_Grid%NX, State_Grid%NY, State_Grid%NZ, NWVAA, NDUST), stat=RC)
    allocate(ODAER(State_Grid%NX, State_Grid%NY, State_Grid%NZ, NWVAA, NAER), stat=RC)
    allocate(ISOPOD(State_Grid%NX, State_Grid%NY, State_Grid%NZ, NWVAA), stat=RC)
    allocate(IRHARR(State_Grid%NX, State_Grid%NY, State_Grid%NZ), stat=RC)
    write(6,*) "Switch_Dims: After CMN_FJX_MOD stuff", RC
    ZPJ = 0e+0_fp
    ODMDUST = 0e+0_fp
    ODAER = 0e+0_fp
    ISOPOD = 0e+0_fp
    IRHARR = 0d0

    ! VDIFF_MOD.F90 does vertical mixing (hplin, 8/20/22)
    call Cleanup_Vdiff(RC)
    write(6,*) "Switch_Dims: After Cleanup_Vdiff", RC
    if(GC_Alloc .and. (Input_Opt%LTURB .and. Input_Opt%LNLPBL)) then
        call Init_Vdiff(Input_Opt, State_Chm, State_Grid, RC)
        write(6,*) "Switch_Dims: After Init_Vdiff", RC
    endif

    ! CARBON_MOD.F mostly contains diagnostics. Safe to re-allocate.
    ! Will be removed in dev/12.8.0 (hplin, 3/19/20)
    call Cleanup_Carbon()
    if(GC_Alloc) then
        call Init_Carbon( Input_Opt, State_Chm, State_Diag, State_Grid, RC )
        write(6,*) "Switch_Dims: After Init_Carbon", RC
    endif

    ! PRESSURE_MOD.F holds Ap, Bp which need to be backed up.
    if(GC_Alloc) then
        ! Saved as [hPa] from Ap, Bp
        ! Remember to convert to [Pa] when passing to HEMCO later
        do L = 1, LM+1
            Ap(L) = GET_AP(L)
            Ap_Pa(L) = GET_AP(L) * 100_hp
            Bp(L) = GET_BP(L)
        enddo
    endif
    call Cleanup_Pressure()

    if(GC_Alloc) then
        ! Write to HcoState not necessary, as stored in HcoState and not pointed
        ! through pointers. Note that the horizontal coordinates are assigned
        ! through pointers to State_Grid, so they should NEVER be re-allocated
        ! hence the check above.

        ! Rellocate the pressure module
        call Init_Pressure( Input_Opt, State_Grid, RC )
        write(6,*) "Switch_Dims: After Init_Pressure", RC

        ! Reassign the Ap, Bp coordinates
        call Accept_External_ApBp( State_Grid, Ap, Bp, RC )

        ! PEDGEs are set later in GIGC_Chunk_Run.
    endif

    ! strat_chem_mod.F90, now linear_chem_mod.F90 in v13.3+
    ! Only used in Linoz. Probably not high priority.
    call Cleanup_Linear_Chem()
    if(GC_Alloc .and. Input_Opt%LINEAR_CHEM) then
        call Init_Linear_Chem( Input_Opt, State_Chm, State_Met, State_Grid, RC )
        write(6,*) "Switch_Dims: After Init_Linear_Chem", RC
    endif

    ! SULFATE_MOD.F
    call CLEANUP_SULFATE()
    if(GC_Alloc .and. Input_Opt%LSULF) then
        call INIT_SULFATE( Input_Opt, State_Chm, State_Diag, State_Grid, RC )
        write(6,*) "Switch_Dims: After Init_Sulfate", RC
    endif

    ! AEROSOL_MOD.F
    call CLEANUP_AEROSOL()
    if(GC_Alloc .and. &
        (Input_Opt%LSULF .or. Input_Opt%LCARB    .or. &
         Input_Opt%LDUST .or. Input_Opt%LSSALT)) then
        call INIT_AEROSOL( Input_Opt, State_Chm, State_Diag, State_Grid, RC )
        write(6,*) "Switch_Dims: After INIT_AEROSOL", RC
    endif

    !
    ! Oh my, here it goes -
    !  We found love on an empty page
    !  Kill the stars above, trying to fight the fade
    !                 - Oh Wonder "Bigger than Love"
    !
    PREVIOUS_ID = ID

    if(GC_Alloc) then
        PREVIOUS_GC_ID = ID
    endif

    write(6, *) "Exit RC: ", RC
    write(6, *) "%%%%%%%%%    FINISHED GIGC_Switch_Dims    %%%%%%%%%"
    
  END SUBROUTINE GIGC_Switch_Dims
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: gchp_chunk_init
!
! !DESCRIPTION: Subroutine GCHP\_CHUNK\_INIT is the ESMF init method for
!  GEOS-Chem.  This routine calls routines within core GEOS-Chem to allocate
!  arrays and read input files.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE GIGC_Chunk_Init( nymdB,         nhmsB,      nymdE,           &
                              nhmsE,         tsChem,     tsDyn,           &
                              tsRad,         lonCtr,     latCtr,          &
                              lonEdge,       latEdge,                     &
                              Input_Opt,     State_Chm,  State_Diag,      &
                              State_Grid,    State_Met,  HcoConfig,       &
                              RC,            MPI_COMM,                    &
                              Ap,            Bp,         ID)
!
! !USES:
!
    USE Chemistry_Mod,           ONLY : Init_Chemistry
    USE Emissions_Mod,           ONLY : Emissions_Init
    USE GC_Environment_Mod

    USE HCO_Types_Mod,           ONLY : ConfigObj
    USE Input_Mod,               ONLY : Read_Input_File
    USE Input_Opt_Mod,           ONLY : OptInput, Set_Input_Opt
    USE Linear_Chem_Mod,         ONLY : Init_Linear_Chem
    USE Linoz_Mod,               ONLY : Linoz_Read
    USE PhysConstants,           ONLY : PI_180
    USE Pressure_Mod,            ONLY : Init_Pressure
    USE Roundoff_Mod,            ONLY : RoundOff
    USE State_Chm_Mod,           ONLY : ChmState, Ind_
    USE State_Diag_Mod,          ONLY : DgnState
    USE State_Grid_Mod,          ONLY : GrdState, Init_State_Grid
    USE State_Met_Mod,           ONLY : MetState
    USE Time_Mod,                ONLY : Set_Timesteps
    USE UCX_MOD,                 ONLY : INIT_UCX
    USE UnitConv_Mod,            ONLY : Convert_Spc_Units

    ! Diagnostics list replacement (WRF-GC)
    USE DiagList_Mod,            ONLY : DgnList, Init_DiagList, Print_DiagList

    ! For conc initialization
    USE Species_Mod,             ONLY : Species

    ! For HEMCO state boot (Don't remember why this is needed, hplin, 12/28/19)
    USE HCO_State_GC_Mod,       ONLY : HcoState, ExtState

    USE GC_Stateful_Mod,         ONLY : GIGC_State_Boot, GIGC_States

    ! Get the Diagnostics List from GIGC_State_Boot
    ! Stored in Global_DiagList
    USE GC_Stateful_Mod,         ONLY : DiagList => Global_DiagList
    USE GC_Stateful_Mod,         ONLY : TaggedDiag_List => Global_TaggedDiag_List

    USE Pressure_Mod,            ONLY : Accept_External_ApBp
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(OptInput),      INTENT(INOUT) :: Input_Opt      ! Input Options object
    TYPE(ChmState),      INTENT(INOUT) :: State_Chm      ! Chem State object
    TYPE(DgnState),      INTENT(INOUT) :: State_Diag     ! Diag State object
    TYPE(GrdState),      INTENT(INOUT) :: State_Grid     ! Grid State object
    TYPE(MetState),      INTENT(INOUT) :: State_Met      ! Met State object
    TYPE(ConfigObj),     POINTER       :: HcoConfig      ! HEMCO config obj
!
! !INPUT PARAMETERS:
! Do not move to above - requires State_Grid to be defined first
!
    INTEGER,            INTENT(IN)    :: nymdB       ! YYYYMMDD @ start of run
    INTEGER,            INTENT(IN)    :: nhmsB       ! hhmmss   @ start of run
    INTEGER,            INTENT(IN)    :: nymdE       ! YYYYMMDD @ end of run
    INTEGER,            INTENT(IN)    :: nhmsE       ! hhmmss   @ end of run
    REAL,               INTENT(IN)    :: tsChem      ! Chemistry timestep [s]
    REAL,               INTENT(IN)    :: tsDyn       ! Chemistry timestep [s]
    REAL,               INTENT(IN)    :: tsRad       ! Chemistry timestep [s]
    REAL(KIND_R4),      INTENT(IN)    :: lonCtr(:,:) ! Lon centers [radians]
    REAL(KIND_R4),      INTENT(IN)    :: latCtr(:,:) ! Lat centers [radians]

    REAL(KIND_R4),      INTENT(IN)    :: lonEdge(:,:)! Lon edges   [radians]
    REAL(KIND_R4),      INTENT(IN)    :: latEdge(:,:)! Lat edges   [radians]
    INTEGER,            INTENT(IN)    :: MPI_COMM    ! MPI Communicator #
    REAL(fp),           INTENT(IN)    :: Ap(State_Grid%NZ+1)    ! Ap coord
    REAL(fp),           INTENT(IN)    :: Bp(State_Grid%NZ+1)    ! Bp coord
    INTEGER,            INTENT(IN)    :: ID          ! Domain identifier within this CPU
!
! !OUTPUT PARAMETERS:
!
    INTEGER,             INTENT(OUT)   :: RC             ! Success or failure?
!
! !REMARKS:
!  Need to add better error checking
!
! !REVISION HISTORY:
!  18 Jul 2011 - M. Long     - Initial Version
!  See https://github.com/geoschem/geos-chem for history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER                        :: I, J, L, STATUS
    CHARACTER(LEN=255)             :: Iam


    TYPE(Species),      POINTER    :: ThisSpc     ! Species pointer for init bg values
    INTEGER                        :: IND         ! Current species index

    !=======================================================================
    ! GIGC_CHUNK_INIT begins here 
    !=======================================================================

    ! Error trap
    Iam = 'GIGC_CHUNK_INIT (gigc_chunk_mod.F90)'

    ! Assume success
    RC = GC_SUCCESS

    ! Update Input_Opt with timing fields
    ! We will skip defining these in READ_INPUT_FILE
    Input_Opt%NYMDb   = nymdB           ! YYYYMMDD @ start of simulation
    Input_Opt%NHMSb   = nhmsB           ! hhmmss   @ end   of simulation
    Input_Opt%NYMDe   = nymdE           ! YYYYMMDD @ start of simulation
    Input_Opt%NHMSe   = nhmsE           ! hhmmss   @ end   of simulation
    Input_Opt%TS_CHEM = INT( tsChem )   ! Chemistry timestep [sec]
    Input_Opt%TS_EMIS = INT( tsChem )   ! Chemistry timestep [sec]
    Input_Opt%TS_DYN  = INT( tsDyn  )   ! Dynamic   timestep [sec]
    Input_Opt%TS_CONV = INT( tsDyn  )   ! Dynamic   timestep [sec]
    Input_Opt%TS_RAD  = INT( tsRad  )   ! RRTMG     timestep [sec]

    ! Initialize GEOS-Chem
    ! GIGC_State_Boot does: Read_Input_File, Linoz_Read, and GC_Allocate_All (only CMN)
    call GIGC_State_Boot(am_I_Root = Input_Opt%AmIRoot,& ! Are we on the root PET?
                         MPI_COMM  = MPI_COMM,         & ! MPI Communicator
                         State_Grid= State_Grid,       & ! Grid State object
                         RC        = RC)
    ! If input options are not read successfully, abort the simulation
    IF ( RC /= GC_SUCCESS ) THEN
      WRITE(6, *) "### GIGC_CHUNK_INIT/INIT_SIMULATION: fatal error in GIGC_State_Boot. stop."
      stop
    ENDIF

    !=======================================================================
    ! Temporary support for in-PET grid switching: verify if we need to
    ! pre-update CMN_SIZE_MOD variables and de/reallocate ahead of time.
    !=======================================================================
    ! Now reallocations are required every time, because GIGC_State_Boot initializes a dummy
    ! sized I/J/L. (hplin, 6/12/18)
    !
    ! chemics_init will pass in %NX, %NY, and %NZ.
    ! GIGC_State_Boot will run SetGridFromCtrEdges
    call GIGC_Switch_Dims( ID,                               &
                           lonCtr, latCtr, lonEdge, latEdge, &
                           Input_Opt, State_Met, State_Chm, State_Diag, State_Grid, &
                           .true., .false., RC )
    write(6, *) "GIGC_Chunk_Init: Received new grid IM, JM, LM =", State_Grid%NX, State_Grid%NY, State_Grid%NZ

    ! GIGC_Switch_Dims will call Read_Input_File (hplin, 4/18/22)
    ! GIGC_Switch_Dims will call GC_Init_Grid (hplin, 12/28/19)

    ! Set maximum number of levels in the chemistry grid
    State_Grid%MaxChemLev  = State_Grid%MaxStratLev

    ! Set GEOS-Chem timesteps on all CPUs
    CALL Set_Timesteps( Input_Opt  = Input_Opt,                              &
                        Chemistry  = Input_Opt%TS_CHEM,                      &
                        Convection = Input_Opt%TS_CONV,                      &
                        Dynamics   = Input_Opt%TS_DYN,                       &
                        Emission   = Input_Opt%TS_EMIS,                      &
                        Radiation  = Input_Opt%TS_RAD,                       &
                        Unit_Conv  = MAX( Input_Opt%TS_DYN,                  &
                                          Input_Opt%TS_CONV ),               &
                        Diagnos    = Input_Opt%TS_DIAG         )

    ! Initialize derived-type objects for met, chem, and diag
    CALL GC_Init_StateObj( DiagList,  TaggedDiag_List, Input_Opt, &
                           State_Chm, State_Diag, State_Grid, State_Met, RC )

    ! Initialize other GEOS-Chem modules
    CALL GC_Init_Extra( DiagList, Input_Opt,    &
                        State_Chm, State_Diag, State_Grid, RC ) 

    ! Feed correct vertical coordinates into pressure_mod... needed for emissions_mod initialization
    ! This needs to be *after* Switch_Dims (which deallocates all) AND GC_init_extra which inits them again (hplin, 4/25/22)
    call ACCEPT_EXTERNAL_APBP(State_Grid, Ap, Bp, RC)

    ! Set initial species units to internal state units, the same
    ! units as the restart file values. Note that species concentrations
    ! are all still zero at this point since internal state values are not
    ! copied to State_Chm%Species%Conc until Run (post-initialization).
# if defined( MODEL_GEOS )
    State_Chm%Spc_Units = 'kg/kg total'
#else
    State_Chm%Spc_Units = 'v/v dry'
#endif

    ! Initialize chemistry mechanism
    IF ( Input_Opt%ITS_A_FULLCHEM_SIM .OR. Input_Opt%ITS_AN_AEROSOL_SIM ) THEN
       CALL INIT_CHEMISTRY ( Input_Opt,  State_Chm, State_Diag, &
                             State_Grid, RC )
    ENDIF

    ! Initialize HEMCO
    CALL EMISSIONS_INIT( Input_Opt, State_Chm, State_Grid, State_Met, RC, &
                         HcoConfig=HcoConfig )

    write(6,*) "Chunk_Init: ID", ID, "State_Grid%ID", State_Grid%ID, "RC", RC
    write(6,*) "Chunk_Init: ID", ID, "%NX, %NY, %NZ", State_Grid%NX, State_Grid%NY, State_Grid%NZ
    write(6,*) "Chunk_Init: ID", ID, "Status XMid",  associated( State_Grid%XMid )
    write(6,*) "Chunk_Init: ID", ID, "Status YEdge", associated( State_Grid%YEdge )
    write(6,*) "Chunk_Init: ID", ID, "Size YEdge", size( State_Grid%YEdge, 1 ), size( State_Grid%YEdge, 2 )
    write(6,*) "State_Grid%YEdge(1,:)", State_Grid%YEdge(1,:)

    ! Stratosphere - can't be initialized without HEMCO because of STATE_PSC
    ! Initialize UCX routines
    CALL INIT_UCX( Input_Opt, State_Chm, State_Diag, State_Grid )

    IF ( Input_Opt%LINEAR_CHEM ) THEN
       CALL INIT_LINEAR_CHEM( Input_Opt, State_Chm, State_Met, State_Grid, RC )
    ENDIF

    !-------------------------------------------------------------------------
    ! Diagnostics and tendencies 
    !-------------------------------------------------------------------------

#if defined( MODEL_WRF )
    !=======================================================================
    ! Initialize simulation with background values.
    !=======================================================================
    ! These values are eventually replaced in State_Chm% by whatever external model
    ! you use to handle chemistry data (in case of WRF, it's the chem array)
    ! every time you run GIGC_Chunk_Run. GIGC in particular is NOT aware of
    ! any location-specific background values.
    !
    ! Ported from Includes_Before_Run.H (GCHP) (hplin, 4/26/2018)
    do I = 1, State_Chm%nSpecies
      ThisSpc => State_Chm%SpcData(I)%Info

      if(trim(ThisSpc%Name) == '') cycle
      IND = IND_(trim(ThisSpc%Name))
      if(IND < 0) cycle

      ! Initialize using background values from species database.
      call SET_BACKGROUND_CONC( Input_Opt%AmIRoot, ThisSpc, &
                               State_Chm, State_Met, State_Grid, &
                               Input_Opt, IND, RC)

      ThisSpc => NULL()
    enddo
    write(6,*) "Chunk_Init: Set background concentrations"
#endif

    !=======================================================================
    ! %%% Replicate GEOS-Chem Classic (non-ESMF) functionality in GIGC %%%
    !
    ! Initializes the History component in GEOS-Chem.
    !=======================================================================
    ! For now, just hardwire the input file for the History component
    Input_Opt%HistoryInputFile = './HISTORY.rc'

    !=======================================================================
    ! Save stateful information into GC_Stateful_Mod (hplin, 5/19/20)
    !=======================================================================
    GIGC_States(ID)%ID   = ID
    GIGC_States(ID)%Init = .true.

    ! Return success
    RC = GC_Success

  END SUBROUTINE GIGC_Chunk_Init
!EOC

!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: gigc_chunk_run
!
! !DESCRIPTION: Subroutine GIGC\_CHUNK\_RUN is the ESMF run method for
!  GEOS-Chem.
!
! !INTERFACE:
!
  SUBROUTINE GIGC_Chunk_Run( ID,                                             &
                             nymd,       nhms,       year,       month,      &
                             day,        dayOfYr,    hour,       minute,     &
                             second,     utc,        hElapsed,   Input_Opt,  &
                             State_Chm,  State_Diag, State_Grid, State_Met,  &
                             lonCtr,     latCtr,     lonEdge,    latEdge,    &
                             Operators,  IsChemTime, IsRadTime,              &
                             RC )
!
! !USES:
!
    ! GEOS-Chem state objects
    USE Input_Opt_Mod,      ONLY : OptInput
    USE State_Chm_Mod,      ONLY : ChmState
    USE State_Diag_Mod
    USE State_Grid_Mod,     ONLY : GrdState
    USE State_Met_Mod,      ONLY : MetState

    ! GEOS-Chem components
    USE Chemistry_Mod,      ONLY : Do_Chemistry, Recompute_OD
    USE Convection_Mod,     ONLY : Do_Convection
    USE DryDep_Mod,         ONLY : Do_DryDep
    USE Emissions_Mod,      ONLY : Emissions_Run
    USE Mixing_Mod,         ONLY : Do_Tend, Do_Mixing
    USE WetScav_Mod,        ONLY : Setup_WetScav, Do_WetDep

    ! HEMCO components (eventually moved to a separate GridComp?)
    USE HCO_State_GC_Mod,   ONLY : HcoState, ExtState
    USE HCO_Interface_Common, ONLY : SetHcoTime
    USE HCO_Interface_GC_Mod, ONLY : Compute_Sflx_For_Vdiff

    ! Specialized subroutines
    USE Calc_Met_Mod,       ONLY : AirQnt
    USE Calc_Met_Mod,       ONLY : Set_Dry_Surface_Pressure
    USE Calc_Met_Mod,       ONLY : Set_Clock_Tracer
    USE Calc_Met_Mod,       ONLY : GCHP_Cap_Tropopause_Prs
    USE Set_Global_CH4_Mod, ONLY : Set_CH4
    USE MODIS_LAI_Mod,      ONLY : Compute_XLAI
    USE PBL_Mix_Mod,        ONLY : Compute_PBL_Height
    USE Pressure_Mod,       ONLY : Set_Floating_Pressures
    USE TOMS_Mod,           ONLY : Compute_Overhead_O3
    USE UCX_Mod,            ONLY : Set_H2O_Trac
    USE Vdiff_Mod,          ONLY : Max_PblHt_for_Vdiff

    ! Utilities
    USE ErrCode_Mod
    USE HCO_Error_Mod
    USE Pressure_Mod,       ONLY : Accept_External_Pedge
    USE State_Chm_Mod,      ONLY : IND_
    USE Time_Mod,           ONLY : Accept_External_Date_Time
    USE UnitConv_Mod,       ONLY : Convert_Spc_Units, Print_Global_Species_Kg

    ! Diagnostics
    USE Diagnostics_Mod,    ONLY : Zero_Diagnostics_StartofTimestep
    USE Diagnostics_Mod,    ONLY : Set_Diagnostics_EndofTimestep
    USE Aerosol_Mod,        ONLY : Set_AerMass_Diagnostic

    USE Species_Mod,   ONLY : Species

    ! Added hplin 12/1/18
    USE Olson_Landmap_Mod,  ONLY : Compute_Olson_Landmap

!
! !INPUT PARAMETERS:
!
    INTEGER,        INTENT(IN)    :: ID          ! Domain identifier within this CPU

    INTEGER,        INTENT(IN)    :: nymd        ! YYYY/MM/DD @ current time
    INTEGER,        INTENT(IN)    :: nhms        ! hh:mm:ss   @ current time
    INTEGER,        INTENT(IN)    :: year        ! UTC year
    INTEGER,        INTENT(IN)    :: month       ! UTC month
    INTEGER,        INTENT(IN)    :: day         ! UTC day
    INTEGER,        INTENT(IN)    :: dayOfYr     ! UTC day of year
    INTEGER,        INTENT(IN)    :: hour        ! UTC hour
    INTEGER,        INTENT(IN)    :: minute      ! UTC minute
    INTEGER,        INTENT(IN)    :: second      ! UTC second
    REAL*4,         INTENT(IN)    :: utc         ! UTC time [hrs]
    REAL*4,         INTENT(IN)    :: hElapsed    ! Elapsed hours
    TYPE(GIGC_Chunk_Operators),  INTENT(IN)    :: Operators   ! Operators to run (derived type)
    LOGICAL,        INTENT(IN)    :: IsChemTime  ! Time for chemistry?
    LOGICAL,        INTENT(IN)    :: IsRadTime   ! Time for RRTMG?

    ! WRF-GC
    REAL(KIND_R4),  INTENT(IN)    :: lonCtr(:,:) ! Lon centers [radians]
    REAL(KIND_R4),  INTENT(IN)    :: latCtr(:,:) ! Lat centers [radians]
    REAL(KIND_R4), INTENT(IN)     :: lonEdge(:,:)! Lon edges   [radians]
    REAL(KIND_R4), INTENT(IN)     :: latEdge(:,:)! Lat edges   [radians]
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(OptInput),      INTENT(INOUT) :: Input_Opt   ! Input Options obj
    TYPE(ChmState),      INTENT(INOUT) :: State_Chm   ! Chemistry State obj
    TYPE(DgnState),      INTENT(INOUT) :: State_Diag  ! Diagnostics State obj
    TYPE(GrdState),      INTENT(INOUT) :: State_Grid  ! Grid State obj
    TYPE(MetState),      INTENT(INOUT) :: State_Met   ! Meteorology State obj
!
! !OUTPUT PARAMETERS:
!
    INTEGER,             INTENT(OUT)   :: RC          ! Return code
!
! !REMARKS:
!
! !REVISION HISTORY:
!  18 Jul 2011 - M. Long     - Initial Version
!  See https://github.com/geoschem/geos-chem for history
!EOP
!------------------------------------------------------------------------------
!BOC
    REAL*8                         :: DT
    CHARACTER(LEN=255)             :: Iam, OrigUnit
    INTEGER                        :: STATUS, HCO_PHASE, RST
#if defined( MODEL_GEOS ) || defined( MODEL_WRF )
    INTEGER                        :: I, J, L
#endif

    ! Local logicals to turn on/off individual components
    ! The parts to be executed are based on the input options,
    ! the time step and the phase.
    LOGICAL                        :: DoConv
    LOGICAL                        :: DoDryDep
    LOGICAL                        :: DoEmis
    LOGICAL                        :: DoTend
    LOGICAL                        :: DoTurb
    LOGICAL                        :: DoChem
    LOGICAL                        :: DoWetDep
    LOGICAL                        :: DoRad

    ! First call?
    LOGICAL, SAVE                  :: FIRST    = .TRUE.
    LOGICAL, SAVE                  :: FIRST_RT = .TRUE. ! RRTMG

    ! # of times this routine has been called. Only temporary for printing 
    ! processes on the first 10 calls.
    INTEGER, SAVE                  :: NCALLS = 0

    ! Strat. H2O settings
    LOGICAL                        :: SetStratH2O 
#if defined( MODEL_GEOS )
    LOGICAL, SAVE                  :: LSETH2O_orig
#endif

    ! For RRTMG
    INTEGER                        :: N

    ! Whether to scale mixing ratio with meteorology update in AirQnt
    LOGICAL, SAVE                  :: scaleMR = .FALSE.

    ! Debug variables
    INTEGER, parameter             :: I_DBG = 6, J_DBG = 5, L_DBG=1

    !=======================================================================
    ! GIGC_CHUNK_RUN begins here 
    !=======================================================================

    ! Error trap
    Iam = 'GIGC_CHUNK_RUN (gigc_chunk_mod.F90)'

    ! Assume success
    RC = GC_SUCCESS

    !=======================================================================
    ! Define processes to be covered in this phase
    !
    ! In the standard GEOS-Chem, the following operator sequence is used:
    ! 1. DryDep (kg)
    ! 2. Emissions (kg)
    ! 3. Turbulence (v/v)
    ! 4. Convection (v/v)
    ! 5. Chemistry (kg)
    ! 6. Wetdep (kg)
    !
    ! The GEOS-5 operator sequence is:
    ! 1. Gravity wave drag
    ! 2. Moist (convection)
    ! 3. Chemistry 1 (drydep and emissions)
    ! 4. Surface 1
    ! 5. Turbulence 1
    ! 6. Surface 2
    ! 7. Turbulence 2
    ! 8. Chemistry 2 (chemistry and wet deposition)
    ! 9. Radiation
    !
    ! Here, we use the following operator sequence:
    !
    ! 1.  Convection (v/v) --> Phase 1
    ! 2.  DryDep (kg)      --> Phase 1
    ! 3.  Emissions (kg)   --> Phase 1
    ! 4a. Tendencies (v/v) --> Phase 1
    ! -------------------------------
    ! 4b. Turbulence (v/v) --> Phase 2
    ! 5.  Chemistry (kg)   --> Phase 2
    ! 6.  WetDep (kg)      --> Phase 2
    !
    ! Any of the listed processes is only executed if the corresponding switch
    ! in the geoschem_config.yml file is enabled. If the physics component
    ! already covers convection or turbulence, they should not be applied here!
    ! The tendencies are only applied if turbulence is not done within
    ! GEOS-Chem (ckeller, 10/14/14).
    !
    ! The standard number of phases in GCHP is 1, set in GCHP.rc, which
    ! results in Phase -1 in gchp_chunk_run. This results in executing
    ! all GEOS-Chem components in a single run rather than splitting up
    ! across two runs as is done in GEOS-5. (ewl, 10/26/18)
    !=======================================================================

    ! By default, do processes as defined in geoschem_config.yml. DoTend
    ! defined below.
    DoConv   = Input_Opt%LCONV                    ! dynamic time step
    DoDryDep = Input_Opt%LDRYD .AND. IsChemTime   ! chemistry time step
    DoEmis   = IsChemTime                         ! chemistry time step
#if defined( MODEL_GEOS )
    DoTurb   = Input_Opt%LTURB .AND. IsChemTime   ! dynamic time step
#else
    DoTurb   = Input_Opt%LTURB                    ! dynamic time step
#endif
    DoChem   = Input_Opt%LCHEM .AND. IsChemTime   ! chemistry time step
    DoWetDep = Input_Opt%LWETD                    ! dynamic time step 
    DoRad    = Input_Opt%LRAD  .AND. IsRadTime    ! radiation time step

    ! Only do selected processes for given operator options.
    DoConv   = DoConv .AND. Operators%Conv
    DoDryDep = DoDryDep .AND. Operators%DryDep
    DoEmis   = DoEmis .AND. Operators%Emis
    DoTurb   = DoTurb .AND. Operators%Turb
    DoChem   = DoChem .AND. Operators%Chem
    DoWetDep = DoWetDep .AND. Operators%WetDep
    DoRad    = DoRad .AND. Operators%Rad

    ! Check if tendencies need be applied. The drydep and emission calls
    ! only calculates the emission / drydep rates, but do not apply the
    ! tendencies to the tracer array yet. If turbulence is done as part of
    ! GEOS-5, we need to make sure that these tendencies are applied to the
    ! tracer array. If turbulence is explicitly covered by GEOS-Chem,
    ! however, the tendencies become automatically applied within the PBL
    ! mixing routines (DO_MIXING), so we should never apply the tendencies
    ! in this case.
#if defined( MODEL_WRF )
    DoTend = ( DoEmis .OR. DoDryDep ) .AND. ((.NOT. Input_Opt%LTURB) .OR. (.NOT. Operators%Turb))
#else
    DoTend = ( DoEmis .OR. DoDryDep ) .AND. .NOT. Input_Opt%LTURB
#endif

    ! testing only
    IF ( Input_Opt%AmIRoot .and. NCALLS < 10 ) THEN
       write(*,*) 'GEOS-Chem Column Code Operators:'
       write(*,*) 'DoConv   : ', DoConv
       write(*,*) 'DoDryDep : ', DoDryDep
       write(*,*) 'DoEmis   : ', DoEmis
       write(*,*) 'DoTend   : ', DoTend
       write(*,*) 'DoTurb   : ', DoTurb
       write(*,*) 'DoChem   : ', DoChem
       write(*,*) 'DoWetDep : ', DoWetDep
       write(*,*) ' '
       write(*,*) 'Write G-C Diagns : ', Operators%GCDiagn
    ENDIF

    !-------------------------------------------------------------------------
    ! Pre-Run assignments
    !-------------------------------------------------------------------------

    ! Zero out certain State_Diag arrays. This should not be done in a phase 2
    ! call since this can erase diagnostics filled during phase 1 (e.g., drydep)
    ! (ckeller, 1/21/2022).
    CALL Zero_Diagnostics_StartOfTimestep( Input_Opt, State_Diag, RC )

    !=======================================================================
    ! Temporary support for in-PET grid switching: verify if we need to
    ! pre-update CMN_SIZE_MOD variables and de/reallocate ahead of time.
    !=======================================================================
    call GIGC_Switch_Dims( ID,                               &
                           lonCtr, latCtr, lonEdge, latEdge, &
                           Input_Opt, State_Met, State_Chm, State_Diag, State_Grid, .true., .true., RC )
    write(6, *) "GIGC_Chunk_Run: Running for in-PET grid #", ID
    write(6, *) "GIGC_Chunk_Run: Received new grid IM, JM, LM =", State_Grid%NX, State_Grid%NY, State_Grid%NZ

    ! Update trop, strat heights
    ! using climatology data from tropopause (WRF) (hplin, 19/12/28)
    ! example: 47 = 38, 44; 72 = 40, 59
    State_Grid%MaxTropLev  = State_Grid%NZ ! # trop. levels below
    State_Grid%MaxStratLev = State_Grid%NZ ! # strat. levels below

    DO L = 1, State_Grid%NZ
        IF(State_Met%TROPP(1, 1) .ge. State_Met%PEDGE(1, 1, L)) THEN
            State_Grid%MaxTropLev = State_Grid%MaxTropLev - 1
        ENDIF
    ENDDO

    ! Compute Olson Landmap State_Met values (IREG, ILAND...)
    ! from State_Met%LandTypeFrac(I,J,N). Copied from chem_gridcompmod, hplin 12/1/18
    ! Compute State_Met variables IREG, ILAND, IUSE, and FRCLND
    CALL Compute_Olson_Landmap( Input_Opt, State_Grid, State_Met, RC )

    ! Pass time values obtained from the ESMF environment to GEOS-Chem
    CALL Accept_External_Date_Time( value_NYMD     = nymd,       &
                                    value_NHMS     = nhms,       &
                                    value_YEAR     = year,       &
                                    value_MONTH    = month,      &
                                    value_DAY      = day,        &
                                    value_DAYOFYR  = dayOfYr,    &
                                    value_HOUR     = hour,       &
                                    value_MINUTE   = minute,     &
                                    value_HELAPSED = hElapsed,   &
                                    value_UTC      = utc,        &
                                    RC             = RC         )

    ! Pass time values obtained from the ESMF environment to HEMCO
#if !defined( MODEL_GEOS )
    CALL SetHcoTime ( HcoState,   ExtState,   year,    month,   day,   &
                      dayOfYr,    hour,       minute,  second,  DoEmis,  RC )
#endif

    ! Calculate MODIS leaf area indexes needed for dry deposition
    ! WRF-GC only: Avoid re-computing the State_Met%XLAI and State_Met%MODISLAI
    ! State_Met%XLAI and State_Met%MODISLAI are transferred from WRF in wrfgc_convert_mod.
#if !defined( MODEL_WRF )
    CALL Compute_XLAI( Input_Opt, State_Grid, State_Met, RC )
#endif

    ! Set the pressure at level edges [hPa] from the ESMF environment
    CALL Accept_External_Pedge( State_Met  = State_Met,   &
                                State_Grid = State_Grid,  &
                                RC         = RC          )

    ! Set dry surface pressure (PS1_DRY) from State_Met%PS1_WET
    CALL SET_DRY_SURFACE_PRESSURE( State_Grid, State_Met, 1 )

    ! Set dry surface pressure (PS2_DRY) from State_Met%PS2_WET
    CALL SET_DRY_SURFACE_PRESSURE( State_Grid, State_Met, 2 )

    ! Initialize surface pressures to match the post-advection pressures
    State_Met%PSC2_WET = State_Met%PS1_WET
    State_Met%PSC2_DRY = State_Met%PS1_DRY
    CALL SET_FLOATING_PRESSURES( State_Grid, State_Met, RC )
    IF ( RC /= GC_SUCCESS ) RETURN

    ! Define airmass and related quantities
    ! Using MODEL_GEOS approach (hplin, 5/22/20)
    CALL AirQnt( Input_Opt, State_Chm, State_Grid, State_Met, RC, .FALSE. )

    ! Initialize/reset wetdep after air quantities computed
    IF ( DoConv .OR. DoChem .OR. DoWetDep ) THEN
       CALL SETUP_WETSCAV( Input_Opt, State_Chm, State_Grid, State_Met, RC )
    ENDIF

    ! Cap the polar tropopause pressures at 200 hPa, in order to avoid
    ! tropospheric chemistry from happening too high up (cf. J. Logan)
    CALL GCHP_Cap_Tropopause_Prs( Input_Opt      = Input_Opt,  &
                                  State_Grid     = State_Grid, &
                                  State_Met      = State_Met,  &
                                  RC             = RC         )

    ! Update clock tracer if relevant
    IF (  IND_('CLOCK','A') > 0 ) THEN
       CALL Set_Clock_Tracer( State_Chm, State_Grid )
    ENDIF

    ! Call PBL quantities. Those are always needed
    CALL Compute_Pbl_Height( Input_Opt, State_Grid, State_Met, RC )
    
    ! Convert to dry mixing ratio
    CALL Convert_Spc_Units ( Input_Opt, State_Chm, State_Grid, State_Met, &
                             'kg/kg dry', RC, OrigUnit=OrigUnit )

    ! SDE 05/28/13: Set H2O to STT if relevant
    IF ( IND_('H2O','A') > 0 ) THEN
       SetStratH2O = .FALSE.
       IF ( Input_Opt%LSETH2O ) THEN
          SetStratH2O = .TRUE.
       ENDIF
       CALL SET_H2O_TRAC( SetStratH2O, Input_Opt, State_Chm, & 
                          State_Grid,  State_Met, RC )

      ! Only force strat once if using UCX
      ! FIXME: hplin this needs to be handled in another way for nested domains --
       IF (Input_Opt%LSETH2O) Input_Opt%LSETH2O = .FALSE.
    ENDIF

    !=======================================================================
    ! EMISSIONS. Pass HEMCO Phase 1 which only updates the HEMCO clock
    ! and the HEMCO data list. Should be called every time to make sure
    ! that the HEMCO clock and the HEMCO data list are up to date.
    !=======================================================================
    HCO_PHASE = 1
    CALL EMISSIONS_RUN( Input_Opt, State_Chm, State_Diag, &
                        State_Grid, State_Met, DoEmis, HCO_PHASE, RC  )

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!!                                PHASE 1 or -1                           !!!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !=======================================================================
    ! 1. Convection
    !
    ! Call GEOS-Chem internal convection routines if convection is enabled
    ! in geoschem_config.yml. This should only be done if convection is not
    ! covered by another gridded component and/or the GC species are not made
    ! friendly to this component!!
    !=======================================================================
    IF ( DoConv ) THEN
       if(Input_Opt%AmIRoot.and.NCALLS<10) write(*,*) ' --- Do convection now'

       CALL DO_CONVECTION ( Input_Opt, State_Chm, State_Diag, &
                            State_Grid, State_Met, RC )

       if(Input_Opt%AmIRoot.and.NCALLS<10) write(*,*) ' --- Convection done!'
    ENDIF

    !=======================================================================
    ! 2. Dry deposition
    !
    ! Calculates the deposition rates in [s-1].
    !=======================================================================
    IF ( DoDryDep ) THEN
       if(Input_Opt%AmIRoot.and.NCALLS<10) THEN
          write(*,*) ' --- Do drydep now'
          write(*,*) '     Use FULL PBL: ', Input_Opt%PBL_DRYDEP
       endif
    
       ! Do dry deposition
       CALL Do_DryDep ( Input_Opt, State_Chm, State_Diag, &
                        State_Grid, State_Met, RC )

       if(Input_Opt%AmIRoot.and.NCALLS<10) write(*,*) ' --- Drydep done!'
    ENDIF

    !=======================================================================
    ! 3. Emissions (HEMCO)
    !
    ! HEMCO must be called on first time step to make sure that the HEMCO
    ! data lists are all properly set up.
    !=======================================================================
    IF ( DoEmis ) THEN
       if(Input_Opt%AmIRoot.and.NCALLS<10) write(*,*) ' --- Do emissions now', DoEmis

       ! Do emissions. Pass HEMCO Phase 2 which performs the emissions
       ! calculations.
       HCO_PHASE = 2
       CALL EMISSIONS_RUN( Input_Opt, State_Chm, State_Diag, &
                           State_Grid, State_Met, DoEmis, HCO_PHASE, RC )

       if(Input_Opt%AmIRoot.and.NCALLS<10) write(*,*) ' --- Emissions done!'
    ENDIF

    !=======================================================================
    ! If physics covers turbulence, simply add the emission and dry
    ! deposition fluxes calculated above to the tracer array, without caring
    ! about the vertical distribution. The tracer tendencies are only added
    ! to the tracers array after emissions, drydep. So we need to use the
    ! emissions time step here.
    !=======================================================================
    IF ( DoTend ) THEN
       if(Input_Opt%AmIRoot.and.NCALLS<10) write(*,*)   &
                           ' --- Add emissions and drydep to tracers'
       DT = HcoState%TS_EMIS 

       ! Apply tendencies over entire PBL. Use emission time step.
       CALL DO_TEND( Input_Opt, State_Chm, State_Diag, &
                     State_Grid, State_Met, .FALSE., RC, DT=DT )

       ! testing only
       if(Input_Opt%AmIRoot.and.NCALLS<10) write(*,*)   &
                                 '     Tendency time step [s]: ', DT 

       if(Input_Opt%AmIRoot.and.NCALLS<10) write(*,*)   &
                                 ' --- Fluxes applied to tracers!' 
    ENDIF ! Tendencies

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!!                              PHASE 2 or -1                             !!!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !=======================================================================
    ! 4. Turbulence
    !
    ! Call GEOS-Chem internal turbulence routines if turbulence is enabled
    ! in geoschem_config.yml. This should only be done if turbulence is not
    ! covered by another gridded component and/or the GC species are not made
    ! friendly to this component!!
    !=======================================================================

    IF ( DoTurb ) THEN
       if(Input_Opt%AmIRoot.and.NCALLS<10) write(*,*) ' --- Do turbulence now'

       ! Only do the following for the non-local PBL mixing
       IF ( Input_Opt%LNLPBL ) THEN

          ! Once the initial met fields have been read in, we need to find
          ! the maximum PBL level for the non-local mixing algorithm.
          ! This only has to be done once. (bmy, 5/28/20)
          !
          ! hplin: this might need change for nested domains
          ! but it appears to be just a 0-D quantity, so leave it here.
          ! do not check for FIRST, because the first call due to off-centered clock
          ! only runs emissions. Check for NCALLS instead as a hack (hplin, 4/25/22)
          !
          ! Note: NCALLS needs to account for MAX_DOM * 2 at the very least,
          ! because GIGC_Chunk_Run will be called empty for the first MAX_DOM times.
          ! Some hours of debugging were wasted here. Bug reported by axzhang from SUSTech.
          ! (hplin, 2/16/23)
          IF ( NCALLS < 17 ) THEN
             CALL Max_PblHt_For_Vdiff( Input_Opt, State_Grid, State_Met, RC )
          ENDIF

          ! Compute the surface flux for the non-local mixing,
          ! (which means getting emissions & drydep from HEMCO)
          ! and store it in State_Chm%Surface_Flux
          CALL Compute_Sflx_For_Vdiff( Input_Opt,  State_Chm, State_Diag,    &
                                       State_Grid, State_Met, RC            )
       ENDIF

       ! Do mixing and apply tendencies. This will use the dynamic time step,
       ! which is fine since this call will be executed on every time step.
       CALL DO_MIXING ( Input_Opt, State_Chm, State_Diag,                    &
                        State_Grid, State_Met, RC                           )

       if(Input_Opt%AmIRoot.and.NCALLS<10) write(*,*) ' --- Turbulence done!'
    ENDIF

    ! Set tropospheric CH4 concentrations and fill species array with
    ! current values.
#if defined( MODEL_GEOS ) || defined ( MODEL_WRF )
    IF ( DoTurb .OR. DoTend ) THEN
#else
    IF ( Phase /= 2 .AND. Input_Opt%ITS_A_FULLCHEM_SIM  &
         .AND. IND_('CH4','A') > 0 ) THEN
#endif
       CALL SET_CH4 ( Input_Opt, State_Chm, State_Diag, &
                      State_Grid, State_Met, RC )
    ENDIF

    !=======================================================================
    ! 5. Chemistry
    !=======================================================================
    IF ( DoChem ) THEN
       if(Input_Opt%AmIRoot.and.NCALLS<10) write(*,*) ' --- Do chemistry now'

       IF ( Input_Opt%ITS_A_FULLCHEM_SIM ) THEN
          ! Calculate TOMS O3 overhead. For now, always use it from the
          ! Met field. State_Met%TO3 is imported from PCHEM (ckeller, 10/21/2014).
          CALL COMPUTE_OVERHEAD_O3( Input_Opt, State_Grid, State_Chm, DAY, &
                                    .TRUE., State_Met%TO3, RC )
       ENDIF

#if !defined( MODEL_GEOS )
       ! Set H2O to species value if H2O is advected
       IF ( IND_('H2O','A') > 0 ) THEN
          CALL SET_H2O_TRAC( .FALSE., Input_Opt, &
                             State_Chm, State_Grid, State_Met, RC )
       ENDIF
#endif

       ! Do chemistry
       CALL Do_Chemistry( Input_Opt, State_Chm, State_Diag, &
                          State_Grid, State_Met, RC ) 

       if(Input_Opt%AmIRoot.and.NCALLS<10) write(*,*) ' --- Chemistry done!'
    ENDIF

    !=======================================================================
    ! 6. Wet deposition
    !=======================================================================
    IF ( DoWetDep ) THEN
       if(Input_Opt%AmIRoot.and.NCALLS<10) write(*,*) ' --- Do wetdep now'

       ! Fix negatives which absolutely break wet deposition.
       ! Note there is a negative flip sign somewhere in univconv_mod. It should not affect us here
       ! but this might be very worth checking. Thanks to xlu for the tip. hplin, 5/24/20
       ! where(State_Chm%Species < 0.0e0)
       !     State_Chm%Species = 1.0e-36
       ! endwhere

       ! Do wet deposition
       CALL DO_WETDEP( Input_Opt, State_Chm, State_Diag, &
                       State_Grid, State_Met, RC )

       if(Input_Opt%AmIRoot.and.NCALLS<10) write(*,*) ' --- Wetdep done!'
    ENDIF

    !=======================================================================
    ! Diagnostics
    !=======================================================================

    !==============================================================
    !      ***** U P D A T E  O P T I C A L  D E P T H *****
    !==============================================================
    ! Recalculate the optical depth at the wavelength(s) specified
    ! in the Radiation Menu. This must be done before the call to any
    ! diagnostic and only on a chemistry timestep.
    ! (skim, 02/05/11)
    IF ( DoChem ) THEN
       CALL RECOMPUTE_OD ( Input_Opt, State_Chm, State_Diag, &
                           State_Grid, State_Met, RC )
    ENDIF

    if(Input_Opt%AmIRoot.and.NCALLS<10) write(*,*) ' --- Do diagnostics now'

    ! Set certain diagnostics dependent on state at end of step. This
    ! includes species concentration and dry deposition flux.
    CALL Set_Diagnostics_EndofTimestep( Input_Opt,  State_Chm, State_Diag, &
                                        State_Grid, State_Met, RC )

    ! Archive aerosol mass and PM2.5 diagnostics
    IF ( State_Diag%Archive_AerMass ) THEN
       CALL Set_AerMass_Diagnostic( Input_Opt,  State_Chm, State_Diag, &
                                    State_Grid, State_Met, RC )
    ENDIF

    if(Input_Opt%AmIRoot.and.NCALLS<10) write(*,*) ' --- Diagnostics done!'

    !=======================================================================
    ! Convert State_Chm%Species units
    !=======================================================================
    CALL Convert_Spc_Units ( Input_Opt, State_Chm, State_Grid, State_Met, &
                             OrigUnit, RC )

    !=======================================================================
    ! Clean up
    !=======================================================================

    ! testing only
    IF ( NCALLS < 10 ) NCALLS = NCALLS + 1 

    ! First call is done
    FIRST = .FALSE.

    ! Return success
    RC = GC_SUCCESS

  END SUBROUTINE GIGC_Chunk_Run

!BOP
  SUBROUTINE GCHP_PRINT_MET(I, J, L,         &
       Input_Opt, State_Grid, State_Met, LOC, RC )

    !
    ! !USES:
    !
    USE State_Met_Mod,        ONLY : MetState
    USE Input_Opt_Mod,        ONLY : OptInput
    USE State_Grid_Mod,       ONLY : GrdState

    !
    ! !INPUT PARAMETERS:
    !
    INTEGER,          INTENT(IN)    :: I         ! Grid cell lat index
    INTEGER,          INTENT(IN)    :: J         ! Grid cell lon index
    INTEGER,          INTENT(IN)    :: L         ! Grid cell lev index
    CHARACTER(LEN=*), INTENT(IN)    :: LOC       ! Call location string
    TYPE(OptInput),   INTENT(IN)    :: Input_Opt ! Input Options object
    TYPE(GrdState),   INTENT(IN)    :: State_Grid! Grid State object
    TYPE(MetState),   INTENT(IN)    :: State_Met ! Meteorology State object
    !
    ! !INPUT/OUTPUT PARAMETERS:
    !

    !
    ! !OUTPUT PARAMETERS:
    !
    INTEGER,          INTENT(OUT)   :: RC        ! Success or failure?!
    ! !REMARKS:
    !
    ! !REVISION HISTORY:
    !EOP
    !------------------------------------------------------------------------------
    !BOC
    !
    ! !LOCAL VARIABLES:
    !
    CHARACTER(LEN=255) :: ErrorMsg, ThisLoc


    !=========================================================================
    ! GCHP_PRINT_MET begins here!
    !=========================================================================

    ErrorMsg  = ''
    ThisLoc   = ' -> at GCHP_Print_Met (in module ' // &
         'Interfaces/GCHP/gchp_chunk_mod.F)'

    ! Assume success
    RC = GC_SUCCESS

    ! Echo info
    IF ( Input_Opt%amIRoot ) THEN
       WRITE( 6, 100 ) TRIM( LOC )
       WRITE( 6, 113 ) State_Grid%YMid(I,J), State_Grid%XMid(I,J)
    ENDIF
100 FORMAT( /, '%%%%% GCHP_PRINT_MET at ', a )
113 FORMAT( 'Lat: ', f5.1, '   Lon: ', f5.1 )

    ! Write formatted output
    IF ( Input_Opt%amIRoot ) THEN
       ! 2-D Fields
       WRITE( 6, 114 ) 'PBLH',     State_Met%PBLH(I,J),     I, J
       WRITE( 6, 114 ) 'PSC2_WET', State_Met%PSC2_WET(I,J), I, J
       WRITE( 6, 114 ) 'PSC2_DRY', State_Met%PSC2_DRY(I,J), I, J
       WRITE( 6, 114 ) 'PS1_WET',  State_Met%PS1_WET(I,J), I, J
       WRITE( 6, 114 ) 'PS1_DRY',  State_Met%PS1_DRY(I,J), I, J
       WRITE( 6, 114 ) 'PS2_WET',  State_Met%PS2_WET(I,J), I, J
       WRITE( 6, 114 ) 'PS2_DRY',  State_Met%PS2_DRY(I,J), I, J
       WRITE( 6, 114 ) 'TS',       State_Met%TS(I,J),       I, J
       WRITE( 6, 114 ) 'U10M',     State_Met%U10M(I,J),     I, J
       ! 3-D Fields
       WRITE( 6, 115 ) 'CLDF',     State_Met%CLDF(I,J,L),      I, J, L
       WRITE( 6, 115 ) 'OMEGA',    State_Met%OMEGA(I,J,L),     I, J, L
       WRITE( 6, 115 ) 'PEDGE',    State_Met%PEDGE(I,J,L),     I, J, L
       WRITE( 6, 115 ) 'T',        State_Met%T(I,J,L),         I, J, L
       WRITE( 6, 115 ) 'U',        State_Met%U(I,J,L),         I, J, L
       WRITE( 6, 115 ) 'V',        State_Met%V(I,J,L),         I, J, L
       WRITE( 6, 115 ) 'AD',       State_Met%AD(I,J,L),        I, J, L
       WRITE( 6, 115 ) 'PREVSPHU', State_Met%SPHU_PREV(I,J,L), I, J, L
       WRITE( 6, 115 ) 'SPHU',     State_Met%SPHU(I,J,L),      I, J, L
       ! terminator
       WRITE( 6, 120 )
    ENDIF
114 FORMAT( 'Grid cell  for ', a8, ' = ', es24.16, ', I,J  = ',2I4 )
115 FORMAT( 'Grid cell  for ', a8, ' = ', es24.16, ', I,J,L= ',3I4 )
120 FORMAT( / )


  END SUBROUTINE GCHP_PRINT_MET
!EOC
END MODULE GIGC_Chunk_Mod
