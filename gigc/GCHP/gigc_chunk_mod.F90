!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!               __          _______  ______       _____  _____                 !
!               \ \        / /  __ \|  ____|     / ____|/ ____|                !
!                \ \  /\  / /| |__) | |__ ______| |  __| |                     !
!                 \ \/  \/ / |  _  /|  __|______| | |_ | |                     !
!                  \  /\  /  | | \ \| |         | |__| | |____                 !
!                   \/  \/   |_|  \_\_|          \_____|\_____|                !
!                                                                              !
!----------------------- ALPHA VERSION, v0.9 (20190802) -----------------------!
!
! WRF-GC: GEOS-Chem High Performance-powered Chemistry Add-On for WRF Model
! Developed by Haipeng Lin <hplin@g.harvard.edu>, Xu Feng <fengx7@pku.edu.cn>
!    January 2018, Peking University, Dept of Atmospheric and Oceanic Sciences
!    Correspondence to: Tzung-May Fu <fuzm@sustech.edu.cn>
!
! ALPHA INFORMATION:
!    WRF-GC Alpha (version 0.9) is experimental. Please notify the authors of
!    any bugs, suggestions and feature requests through email or the GitHub
!    repository.
!
! COPYRIGHT STATEMENT:
!    Permission is hereby granted, free of charge, to any person obtaining a copy
!   of this software and associated documentation files (the "Software"), to 
!   use, copy, modify the Software, and to permit persons to whom the Software is
!   furnished to do so, subject to the following conditions:
!
!   - The above copyright notice and this permission notice shall be included in all
!   copies or substantial portions of the Software.
! 
!   - The Software, modified in part or in full may not be redistributed without
!   express permission from the copyright holder(s).
! 
!   Except as contained in this notice or in attribution, the name of the WRF-GC model
!   shall not be used as an endorsement for distributing modified copies of the
!   Software without prior written permission from the copyright holder(s).
! 
!   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
!   IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
!   FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
!   AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
!   LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
!   OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
!   SOFTWARE.
! 
!  WRF and the GEOS-Chem model, GCHP are (c) their original authors.
!
!------------------------------------------------------------------------------- 

MODULE GIGC_Chunk_Mod
  USE GCHP_Utils
  USE Input_Opt_Mod,           ONLY : OptInput
  USE State_Chm_Mod,           ONLY : ChmState
  USE State_Diag_Mod,          ONLY : DgnState
  USE State_Met_Mod,           ONLY : MetState
  USE HCO_TYPES_MOD,           ONLY : ConfigObj

  USE GIGC_Stateful_Mod

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: GIGC_Chunk_Init
  PUBLIC :: GIGC_Chunk_Run
  PUBLIC :: GIGC_Chunk_Final

  PRIVATE :: GIGC_Switch_Dims

  TYPE GC_DIAG
     LOGICAL                    :: DO_PRINT     ! Should we print out?
     INTEGER                    :: N_DIAG       ! # of diag quantities
     INTEGER                    :: COUNT        ! Counter for averaging
     CHARACTER(LEN=10), POINTER :: NAME(:)      ! Tracer names
     REAL*8,            POINTER :: TRACER(:,:)  ! Tracer concentrations
     CHARACTER(LEN=40)          :: FILENAME     ! File name for output
     INTEGER                    :: LUN          ! File unit # for output
  END TYPE GC_DIAG

  TYPE(GC_DIAG)                 :: DIAG_COL

  TYPE GIGC_Chunk_Operators
    logical                     :: Conv
    logical                     :: DryDep
    logical                     :: Emis
    logical                     :: Tend
    logical                     :: Turb
    logical                     :: Chem
    logical                     :: WetDep

    logical                     :: GCDiagn      ! Use GEOS-Chem History-based GCC/NetCDF Diagns?
  END TYPE GIGC_Chunk_Operators
  PUBLIC :: GIGC_Chunk_Operators

  INTEGER, PARAMETER :: KIND_R4 = selected_real_kind(3,25)
CONTAINS
  SUBROUTINE GIGC_Switch_Dims( am_I_Root,                                   &
                               ID,         IM,          JM,         LM,     &
                               lonCtr,     latCtr,      lonEdge,    latEdge,&
                               Input_Opt,                                   &
                               State_Met,  State_Chm,   State_Diag,         &
                               CMN_Alloc,  GC_Alloc,                        &
                               RC )
    USE ErrCode_Mod

    USE CMN_SIZE_MOD
    USE CMN_O3_MOD
    USE CMN_FJX_MOD,    ONLY : ZPJ, ODMDUST, ODAER, ISOPOD, IRHARR, JVN_, NWVAA, NDUST, NAER
    USE CMN_FJX_MOD,    ONLY : L_, L1_, L2_, JVL_, JXL_, JXL1_, JXL2_, JTAUMX, N_

    USE VDIFF_PRE_Mod,  ONLY : Init_Vdiff_Pre
    USE GC_GRID_MOD
    USE GC_Environment_Mod, ONLY : GC_Init_Grid, GC_Init_Regridding
    USE PRESSURE_MOD,   ONLY : Cleanup_Pressure, Init_Pressure, GET_AP, GET_BP
    USE GRID_REGISTRY_MOD, ONLY : Init_Grid_Registry, Cleanup_Grid_Registry
    USE CARBON_MOD,     ONLY : INIT_CARBON, CLEANUP_CARBON
    USE AEROSOL_MOD,    ONLY : INIT_AEROSOL, CLEANUP_AEROSOL
    USE DIAG_OH_MOD,    ONLY : INIT_DIAG_OH, CLEANUP_DIAG_OH
    USE UCX_MOD,        ONLY : INIT_UCX, CLEANUP_UCX
    USE SULFATE_MOD,    ONLY : INIT_SULFATE, CLEANUP_SULFATE
    USE C2H6_MOD,       ONLY : INIT_C2H6, CLEANUP_C2H6
    USE HCOI_GC_Main_Mod, ONLY : HCOI_GC_Init, HCOI_GC_Final
    USE PBL_MIX_MOD,    ONLY : INIT_PBL_MIX, CLEANUP_PBL_MIX
    USE Regrid_A2A_Mod, ONLY : Cleanup_Map_A2A
    USE TOMS_MOD,       ONLY : Init_Toms, Cleanup_Toms

    USE HCO_Error_Mod,  ONLY : hp ! For hemco precision
    USE HCO_VertGrid_Mod, ONLY : HCO_VertGrid_Define, HCO_VertGrid_Cleanup

    USE HCO_INTERFACE_MOD, ONLY : HcoState, ExtState

    USE Input_Opt_Mod,  ONLY : OptInput
    USE State_Met_Mod,  ONLY : MetState
    USE State_Chm_Mod,  ONLY : ChmState
    USE State_Diag_Mod, ONLY : DgnState

    USE GIGC_Stateful_Mod
    LOGICAL,            INTENT(IN)    :: am_I_Root    ! Are we on the root CPU?
    INTEGER,            INTENT(IN)    :: ID           ! Domain identifier within this CPU
    INTEGER,            INTENT(IN)    :: IM           ! # lons, this CPU
    INTEGER,            INTENT(IN)    :: JM           ! # lats, this CPU
    INTEGER,            INTENT(IN)    :: LM           ! # levs, this CPU
    REAL(KIND_R4),      INTENT(IN)    :: lonCtr(:,:)  ! Lon centers [radians]
    REAL(KIND_R4),      INTENT(IN)    :: latCtr(:,:)  ! Lat centers [radians]
    REAL(KIND_R4),      INTENT(IN)    :: lonEdge(:,:) ! Lat centers [radians]
    REAL(KIND_R4),      INTENT(IN)    :: latEdge(:,:) ! Lat centers [radians]
    LOGICAL,            INTENT(IN)    :: CMN_Alloc    ! Do CMN allocations? If you skip GC_Allocate_All, use this
    LOGICAL,            INTENT(IN)    :: GC_Alloc     ! Do other allocations? Not in init

    TYPE(OptInput),     INTENT(INOUT) :: Input_Opt
    TYPE(MetState),     INTENT(INOUT) :: State_Met
    TYPE(ChmState),     INTENT(INOUT) :: State_Chm
    TYPE(DgnState),     INTENT(INOUT) :: State_Diag

    INTEGER,            INTENT(OUT)   :: RC
    INTEGER, SAVE                     :: PREVIOUS_ID = -1
    INTEGER, SAVE                     :: PREVIOUS_GC_ID = -1
    LOGICAL                           :: HCO_ERROR
    LOGICAL, SAVE                     :: FIRST = .TRUE.
    LOGICAL                           :: INIT_ID

    REAL(hp)                          :: Ap(LM+1), Bp(LM+1)
    INTEGER                           :: L

    CALL GIGC_State_Get_Status(am_I_Root, ID, INIT_ID)

    if(PREVIOUS_GC_ID .eq. ID) then
        write(6, *) "Switching is not necessary as we are already"
        write(6, *) "in the right domain!"

        return
    endif

    write(6, *) "%%%%%%%%% GIGC IN-PET GRID SWITCHING HACKS %%%%%%%%%"
    write(6, *) "    (Not very) Proudly brought to you by hplin      "
    write(6, *) "Some debug information:"
    write(6, *) "ID:", ID
    write(6, *) "IM, JM, LM:", IM, JM, LM, " will be set as PAR"
    write(6, *) "Allocs: CMN", CMN_Alloc, "GC", GC_Alloc
    write(6, *) "RC at entry-point:", RC
    write(6, *) "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"



    call Cleanup_CMN_SIZE(am_I_Root, RC)

    IIPAR = IM
    JJPAR = JM
    LLPAR = LM

    IGLOB = IIPAR
    JGLOB = JJPAR
    LGLOB = LLPAR

    LLTROP  = LM
    LLSTRAT = LM

    IF ( Input_Opt%LUCX ) THEN
      LLCHEM     = LLSTRAT
      LLCHEM_FIX = LLSTRAT
    ELSE
      LLCHEM     = LLTROP
      LLCHEM_FIX = LLTROP_FIX
    ENDIF

    IM_WORLD = IM
    JM_WORLD = JM
    LM_WORLD = LM

    if(CMN_Alloc) then
      ALLOCATE( DLON( IIPAR, JJPAR, LLPAR ), STAT=RC )
      DLON = 0e+0_fpp

      ALLOCATE( DLAT( IIPAR, JJPAR, LLPAR ), STAT=RC )
      DLAT = 0e+0_fpp
    endif

    call Cleanup_CMN_O3(am_I_Root, RC)
    if(CMN_Alloc) call Init_CMN_O3(am_I_Root, RC)

    if(allocated(ZPJ)) deallocate(ZPJ)
    if(allocated(ODMDUST)) deallocate(ODMDUST)
    if(allocated(ODAER)) deallocate(ODAER)
    if(allocated(ISOPOD)) deallocate(ISOPOD)
    if(allocated(IRHARR)) deallocate(IRHARR)
    if(CMN_Alloc) then
      allocate(ZPJ(LLPAR, JVN_, IIPAR, JJPAR))
      allocate(ODMDUST(IIPAR, JJPAR, LLPAR, NWVAA, NDUST))
      allocate(ODAER(IIPAR, JJPAR, LLPAR, NWVAA, NAER))
      allocate(ISOPOD(IIPAR, JJPAR, LLPAR, NWVAA))
      allocate(IRHARR(IIPAR, JJPAR, LLPAR))

      L_     = LLPAR    ! Number of CTM layers
      L1_    = L_+1     ! Number of CTM layer edges
      L2_    = L1_*2    ! Number of levels in FJX grid that
      JVL_   = LLPAR    ! Vertical levels for J-values

      JXL_   = LLPAR    ! Vertical levels for J-values computed within Fast-JX
      JXL1_  = JXL_+1   ! Vertical levels edges for J-values
      JXL2_  = 2*JXL_+2 ! Max # levels in the basic Fast-JX grid (mid-level)

      JTAUMX = ( N_ - 4*JXL_ ) / 2  ! Maximum number of divisions ( i.e., may
    endif

    if(PREVIOUS_ID .ne. -1) then
      call GIGC_State_Set_HCO(am_I_Root, PREVIOUS_ID, HcoState, RC)
      call GIGC_State_Set_HCOX(am_I_Root, PREVIOUS_ID, ExtState, RC)
    endif

    if(PREVIOUS_ID .ne. -1 .and. PREVIOUS_ID .ne. ID) then
      if(.not. INIT_ID) then
        nullify(HcoState)
        nullify(ExtState)
        write(6,*) "HcoState for ID not previously init'd, assigning", ID
      else
        CALL GIGC_State_Get_HCO(am_I_Root, ID, HcoState, RC)
        CALL GIGC_State_Get_HCOX(am_I_Root, ID, ExtState, RC)
        write(6,*) "Switch HcoState from PREV->ID", PREVIOUS_ID, ID
        write(6,*) "Sanity check assoc", associated(HcoState), associated(ExtState)
        write(6,*) "Locations of Hco, Ext", loc(HcoState), loc(ExtState)
      endif
    else
      write(6,*) "No need to switch HcoState for ID", ID
    endif
    if(associated(HcoState)) then
      write(6,*) "Debug HcoState nSpc, NX, NY, NZ", HcoState%nSpc, HcoState%NX, HcoState%NY, HcoState%NZ
    endif

    call Cleanup_Grid(am_I_Root, RC)
    if(CMN_Alloc) then
      call GC_Init_Grid( am_I_Root, Input_Opt, RC )

      call SetGridFromCtrEdges(am_I_Root, IM, JM, lonCtr, latCtr, lonEdge, latEdge, RC)
      if(RC /= GC_SUCCESS) then
        write(6, *) "GIGC_Switch_Dims: Failed to SetGridFromCtrEdges! Abort"
        write(6, *) "lonCtr(s):"
        write(6, *) lonCtr
        write(6, *) "latCtr(s):"
        write(6, *) latCtr
        write(6, *) "lonEdge(s):"
        write(6, *) lonEdge
        write(6, *) "latEdge(s):"
        write(6, *) latEdge
        return
      endif

      if(CMN_Alloc .and. GC_Alloc .and. associated(State_Met%AREA_M2)) then
        AREA_M2 = State_Met%AREA_M2
      endif

      if(GC_Alloc) then
        DO L = 1, LLPAR + 1
          Ap(L) = GET_AP(L) * 100_hp
          Bp(L) = GET_BP(L)
        ENDDO

        call HCO_VertGrid_Cleanup(HcoState%Grid%zGrid)
        call HCO_VertGrid_Define(am_I_Root, HcoState%Config, zGrid = HcoState%Grid%zGrid, &
                                 nz = LM, Ap = Ap, Bp = Bp, RC = RC)

        HcoState%Grid%XMID%Val       => XMID   (:,:,1)
        HcoState%Grid%YMID%Val       => YMID   (:,:,1)
        HcoState%Grid%XEDGE%Val      => XEDGE  (:,:,1)
        HcoState%Grid%YEDGE%Val      => YEDGE  (:,:,1)
        HcoState%Grid%YSIN%Val       => YSIN   (:,:,1)
        HcoState%Grid%AREA_M2%Val    => AREA_M2(:,:,1)
      endif
    endif

    if(CMN_Alloc) then
      call Cleanup_Map_A2A()
      call GC_Init_Regridding(am_I_Root, Input_Opt, RC)
    endif


    call CLEANUP_PBL_MIX()
    if(GC_Alloc) call INIT_PBL_MIX(am_I_Root, RC)

    call Cleanup_Pressure()
    if(GC_Alloc) call INIT_PRESSURE(am_I_Root)


    call Cleanup_Toms(am_I_Root, RC)
    if(GC_Alloc) call Init_Toms( am_I_Root, Input_Opt, State_Chm, State_Diag, RC )

    call CLEANUP_CARBON()
    if(GC_Alloc) call INIT_CARBON( am_I_Root, Input_Opt, State_Chm, State_Diag, RC )

    call CLEANUP_AEROSOL()
    if(GC_Alloc) call INIT_AEROSOL( am_I_Root, Input_Opt, State_Chm, State_Diag, RC )

    call CLEANUP_DIAG_OH()
    if(GC_Alloc) call INIT_DIAG_OH(am_I_Root, Input_Opt, RC)

    call CLEANUP_UCX(am_I_Root)
    if(GC_Alloc) call INIT_UCX(am_I_Root, Input_Opt, State_Chm, State_Diag)

    call CLEANUP_SULFATE()
    if(GC_Alloc) call INIT_SULFATE(am_I_Root, Input_Opt, State_Chm, State_Diag, RC)

    call CLEANUP_C2H6()
    if(GC_Alloc) call INIT_C2H6(am_I_Root, Input_Opt, RC)

    FIRST = .FALSE.
    PREVIOUS_ID = ID

    if(GC_Alloc) then
        PREVIOUS_GC_ID = ID
    endif

    write(6, *) "Exit RC: ", RC
    write(6, *) "%%%%%%%%%    FINISHED GIGC_Switch_Dims    %%%%%%%%%"
    
  END SUBROUTINE GIGC_Switch_Dims
  SUBROUTINE GIGC_Chunk_Init( am_I_Root, I_LO,      J_LO,       I_HI,      &
                              J_HI,      IM,        JM,         LM,        &
                              ID,                                          &
                              IM_WORLD,  JM_WORLD,  LM_WORLD,   nymdB,     &
                              nhmsB,     nymdE,     nhmsE,                 &
                              tsChem,    tsDyn,                            &
                              lonCtr,    latCtr,    lonEdge,    latEdge,   & 
                              myPET,                                       &
                              Input_Opt, State_Chm, State_Diag, State_Met, &
                              HcoConfig, RC,        MPI_COMM)
    USE Chemistry_Mod,           ONLY : Init_Chemistry
    USE CMN_Size_Mod,            ONLY : IIPAR, JJPAR, LLPAR, dLon, dLat
    USE DiagList_Mod
    USE Emissions_Mod,           ONLY : Emissions_Init
    USE Fast_JX_Mod,             ONLY : Init_FJX
    USE GC_Environment_Mod
    USE GC_Grid_Mod,             ONLY : SetGridFromCtrEdges
    USE HCO_Types_Mod,           ONLY : ConfigObj
    USE Input_Mod,               ONLY : Read_Input_File
    USE Input_Opt_Mod,           ONLY : OptInput, Set_Input_Opt
    USE Linoz_Mod,               ONLY : Linoz_Read
    USE PBL_Mix_Mod,             ONLY : Init_PBL_Mix
    USE PhysConstants,           ONLY : PI_180
    USE Pressure_Mod,            ONLY : Init_Pressure
    USE Roundoff_Mod,            ONLY : RoundOff
    USE State_Chm_Mod,           ONLY : ChmState
    USE State_Diag_Mod,          ONLY : DgnState
    USE State_Met_Mod,           ONLY : MetState
    USE Time_Mod,                ONLY : Set_Timesteps
    USE UCX_MOD,                 ONLY : INIT_UCX
    USE UnitConv_Mod,            ONLY : Convert_Spc_Units

    USE Species_Mod,             ONLY : Species
    USE State_Chm_Mod,           ONLY : IND_
    USE ErrCode_Mod
    USE Error_Mod,               ONLY : Debug_Msg
    USE HCO_INTERFACE_MOD,       ONLY : HcoState, ExtState

    USE History_Mod,             ONLY : History_Init

    USE GIGC_Stateful_Mod
    LOGICAL,            INTENT(IN)    :: am_I_Root   ! Are we on the root CPU?
    INTEGER,            INTENT(IN)    :: I_LO        ! Min lon index, this CPU
    INTEGER,            INTENT(IN)    :: J_LO        ! Min lat index, this CPU
    INTEGER,            INTENT(IN)    :: I_HI        ! Max lon index, this CPU
    INTEGER,            INTENT(IN)    :: J_HI        ! Max lat index, this CPU
    INTEGER,            INTENT(IN)    :: IM          ! # lons, this CPU
    INTEGER,            INTENT(IN)    :: JM          ! # lats, this CPU
    INTEGER,            INTENT(IN)    :: LM          ! # levs, this CPU
    INTEGER,            INTENT(IN)    :: ID          ! Domain identifier within this CPU
    INTEGER,            INTENT(IN)    :: IM_WORLD    ! # lons, global grid
    INTEGER,            INTENT(IN)    :: JM_WORLD    ! # lats, global grid
    INTEGER,            INTENT(IN)    :: LM_WORLD    ! # levs, global grid
    INTEGER,            INTENT(IN)    :: myPET       ! Local PET
    INTEGER,            INTENT(IN)    :: nymdB       ! YYYYMMDD @ start of run
    INTEGER,            INTENT(IN)    :: nhmsB       ! hhmmss   @ start of run
    INTEGER,            INTENT(IN)    :: nymdE       ! YYYYMMDD @ end of run
    INTEGER,            INTENT(IN)    :: nhmsE       ! hhmmss   @ end of run
    REAL,               INTENT(IN)    :: tsChem      ! Chemistry timestep [s]
    REAL,               INTENT(IN)    :: tsDyn       ! Chemistry timestep [s]
    REAL(KIND_R4), INTENT(IN)    :: lonCtr(:,:)      ! Lon centers [radians]
    REAL(KIND_R4), INTENT(IN)    :: latCtr(:,:)      ! Lat centers [radians]
    REAL(KIND_R4), INTENT(IN)    :: lonEdge(:,:)     ! Lon edges   [radians]
    REAL(KIND_R4), INTENT(IN)    :: latEdge(:,:)     ! Lat edges   [radians]

    INTEGER, INTENT(IN)  :: MPI_COMM                 ! MPI Communicator #
    TYPE(OptInput),     INTENT(INOUT) :: Input_Opt   ! Input Options object
    TYPE(ChmState),     INTENT(INOUT) :: State_Chm   ! Chemistry State object
    TYPE(DgnState),     INTENT(INOUT) :: State_Diag  ! Diagnostics State object
    TYPE(MetState),     INTENT(INOUT) :: State_Met   ! Meteorology State object
    TYPE(ConfigObj),    POINTER       :: HcoConfig   ! HEMCO config obj 
    INTEGER,            INTENT(OUT)   :: RC          ! Success or failure?
    INTEGER                           :: STATUS
    TYPE(Species),      POINTER       :: ThisSpc     ! Species pointer for init bg values
    INTEGER                           :: IND         ! Current species index
    INTEGER                           :: I           ! Loop idx...
    INTEGER                           :: II, JJ, LL
    LOGICAL                           :: prtDebug
    LOGICAL,            SAVE          :: FIRST = .TRUE.

    TYPE(DgnList)                     :: DiagList

    RC = GC_SUCCESS

    IF ( am_I_Root ) THEN 
       WRITE(6, *) "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"
       WRITE(6, *) "%               GEOS-CHEM HIGH PERFORMANCE               %"
       WRITE(6, *) "%                                                        %"
       WRITE(6, *) "% This GEOS-Chem column code was based on work by the    %"
       WRITE(6, *) "% original GCHP authors.                                 %"
       WRITE(6, *) "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"
       WRITE(6, *) "Simulation Start Date/Time:", nymdB, nhmsB
       WRITE(6, *) "Simulation End Date/Time:", nymdE, nhmsE
       WRITE(6, *) "Chemical Timestep (s):", tsChem
       WRITE(6, *) "Dynamics Timestep (s):", tsDyn
       WRITE(6, *) "%%% CONFIGURATION FOR DOMAIN", ID, "%%%"
       WRITE(6, *) "In-PET LON IDX:", I_LO, I_HI
       WRITE(6, *) "In-PET LAT IDX:", J_LO, J_HI
       WRITE(6, *) "IM, JM, LM:", IM, JM, LM
       WRITE(6, *) "%%% WORLD CONFIGURATION %%%"
       WRITE(6, *) "World IM, JM, LM:", IM_WORLD, JM_WORLD, LM_WORLD
    ENDIF


    call GIGC_State_Boot(am_I_Root = Am_I_Root,        & ! Are we on the root PET?
                         MPI_COMM  = MPI_COMM,         & ! MPI Communicator
                         RC        = RC)

    IF ( RC /= GC_SUCCESS ) THEN
      WRITE(6, *) "### GIGC_CHUNK_INIT/INIT_SIMULATION: fatal error in GIGC_State_Boot. stop."
      stop
    ENDIF

    IF ( prtDebug ) THEN
      CALL DEBUG_MSG( '### GIGC_INIT_SIMULATION: after GIGC_State_Boot' )
    ENDIF

    call GIGC_State_Get_Opt(am_I_Root, Input_Opt)

    call GIGC_State_Get_DiagList(am_I_Root, DiagList)

    call GIGC_Switch_Dims( am_I_Root, ID, IM, JM, LM,        &
                           lonCtr, latCtr, lonEdge, latEdge, &
                           Input_Opt, State_Met, State_Chm, State_Diag, .true., .false., RC )
    write(6, *) "GIGC_Chunk_Init: Received new grid IM, JM, LM =", IM, JM, LM

    prtDebug = .true.

    Input_Opt%NYMDb   = nymdB
    Input_Opt%NHMSb   = nhmsB
    Input_Opt%NYMDe   = nymdE
    Input_Opt%NHMSe   = nhmsE
    Input_Opt%TS_CHEM = INT( tsChem )   ! Chemistry timestep [sec]
    Input_Opt%TS_EMIS = INT( tsChem )   ! Chemistry timestep [sec]
    Input_Opt%TS_DYN  = INT( tsDyn  )   ! Dynamic   timestep [sec]
    Input_Opt%TS_CONV = INT( tsDyn  )   ! Dynamic   timestep [sec]

    Input_Opt%myCPU = myPET

    IF ( prtDebug ) THEN
      CALL DEBUG_MSG( '### GIGC_INIT_SIMULATION: after Input_Opt% setups' )
    ENDIF

    DO LL = 1, LLPAR
    DO JJ = 1, JJPAR
    DO II = 1, IIPAR
       dLon(II,JJ,LL)    = RoundOff( (lonEdge(II+1, JJ) - lonEdge(II, JJ)) / PI_180, 4 )

       dLat(II,JJ,LL)    = RoundOff( (latEdge(II, JJ+1) - latEdge(II, JJ)) / PI_180, 4 )
    ENDDO
    ENDDO
    ENDDO



    CALL Set_Timesteps( am_I_Root  = am_I_Root,                          &
                        Chemistry  = Input_Opt%TS_CHEM,                  &
                        Convection = Input_Opt%TS_CONV,                  &
                        Dynamics   = Input_Opt%TS_DYN,                   &
                        Emission   = Input_Opt%TS_EMIS,                  &
                        Radiation  = Input_Opt%TS_RAD,                   &
                        Unit_Conv  = MAX( Input_Opt%TS_DYN,              &
                                          Input_Opt%TS_CONV ),           &
                        Diagnos    = Input_Opt%TS_DIAG         )

    CALL GC_Init_StateObj( am_I_Root, DiagList,   Input_Opt, &
                           State_Chm, State_Diag, State_Met, RC )
    IF ( RC /= GC_SUCCESS ) RETURN

    IF ( prtDebug ) THEN
      CALL DEBUG_MSG( '### GIGC_INIT_SIMULATION: after GC_Init_All' )
    ENDIF

    CALL GC_Init_Extra( am_I_Root, DiagList,   Input_Opt,    &
                        State_Chm, State_Diag, RC ) 
    IF ( RC /= GC_SUCCESS ) RETURN

    IF ( prtDebug ) THEN
      CALL DEBUG_MSG( '### GIGC_INIT_SIMULATION: after GC_Init_Extra' )
    ENDIF

    IF ( Input_Opt%ITS_A_FULLCHEM_SIM .OR.                     &
         Input_Opt%ITS_AN_AEROSOL_SIM ) THEN
       CALL Init_FJX( am_I_Root, Input_Opt, State_Chm, State_Diag, RC ) 
       IF ( RC /= GC_SUCCESS ) RETURN
          
          IF ( prtDebug ) THEN
             CALL DEBUG_MSG( '### GIGC_INIT_SIMULATION: after INIT_FJX' )        
          ENDIF
    ENDIF

    State_Chm%Spc_Units = 'kg/kg dry'

    CALL Init_Pressure( am_I_Root )

    if(am_I_Root) then
        write(*,*) '# GIGC_CHUNK_MOD: after Init_Pressure'
    endif



    CALL Init_PBL_Mix( am_I_Root, RC )
    IF ( RC /= GC_SUCCESS ) RETURN

    if(am_I_Root) then
        write(*,*) '# GIGC_CHUNK_MOD: after Init_PBL_Mix'
    endif
    
    IF ( Input_Opt%ITS_A_FULLCHEM_SIM .OR. Input_Opt%ITS_AN_AEROSOL_SIM ) THEN
       CALL Init_Chemistry( am_I_Root, Input_Opt, State_Chm, State_Diag, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

    if(am_I_Root) then
        write(*,*) '# GIGC_CHUNK_MOD: after Init_Chemistry'
    endif

    CALL EMISSIONS_INIT ( am_I_Root, Input_Opt, State_Met, State_Chm, RC, &
                          HcoConfig=HcoConfig )

    if(am_I_Root) then
        write(*,*) '# GIGC_CHUNK_MOD: after EMISSIONS_INIT'
    endif

    IF ( Input_Opt%LUCX ) THEN

       CALL INIT_UCX( am_I_Root, Input_Opt, State_Chm, State_Diag )

    ENDIF

    CALL Convert_Spc_Units( am_I_Root, Input_Opt, State_Met, &
                            State_Chm, 'v/v dry', RC )

    do I = 1, State_Chm%nSpecies
      ThisSpc => State_Chm%SpcData(I)%Info

      if(trim(ThisSpc%Name) == '') cycle
      IND = IND_(trim(ThisSpc%Name))
      if(IND < 0) cycle

      call SET_BACKGROUND_CONC(am_I_Root, ThisSpc, State_Chm, State_Met, Input_Opt, IND, RC)

      where(State_Chm%Species < 0.0e0)
        State_Chm%Species = 1.0e-36
      endwhere

      ThisSpc => NULL()
    enddo

    Input_Opt%HistoryInputFile = './HISTORY.rc'


    CALL GIGC_State_Set_Opt(am_I_Root, Input_Opt, HcoConfig)
    CALL GIGC_State_Init(am_I_Root    = am_I_Root,   &
                         ID           = ID,          &
                         State_Met    = State_Met,   &
                         State_Chm    = State_Chm,   &
                         State_Diag   = State_Diag,  &
                         HcoState     = HcoState,    &
                         ExtState     = ExtState,    &
                         RC           = RC           )

    write(6, *) '# GIGC_CHUNK_MOD: after SV save, ID =', ID

    FIRST = .FALSE.

  END SUBROUTINE GIGC_Chunk_Init

  SUBROUTINE GIGC_Chunk_Run( am_I_Root,                                   &
                             IM,        JM,        LM,        ID,         &
                             nymd,      nhms,      year,      month,      &
                             day,       dayOfYr,   hour,      minute,     &
                             second,    hElapsed,                         &
                             lonCtr,    latCtr,    lonEdge,   latEdge,    &
                             Input_Opt, State_Chm, State_Met, State_Diag, &
                             Operators, IsChemTime,                       &
                             RC                                           )
    USE HCO_Interface_Mod,  ONLY : HcoState
    USE Input_Opt_Mod,      ONLY : OptInput
    USE State_Chm_Mod,      ONLY : ChmState
    USE State_Diag_Mod
    USE State_Met_Mod,      ONLY : MetState

    USE Aerosol_Mod,        ONLY : Set_AerMass_Diagnostic
    USE Chemistry_Mod,      ONLY : Do_Chemistry, Recompute_OD
    USE Convection_Mod,     ONLY : Do_Convection
    USE DryDep_Mod,         ONLY : Do_DryDep
    USE Emissions_Mod,      ONLY : Emissions_Run
    USE Mixing_Mod,         ONLY : Do_Tend, Do_Mixing
    USE Strat_Chem_Mod,     ONLY : Init_Strat_Chem, Minit_is_Set
    USE WetScav_Mod,        ONLY : Setup_WetScav, Do_WetDep

    USE Dao_Mod,            ONLY : AirQnt, Set_Dry_Surface_Pressure
    USE Dao_Mod,            ONLY : GIGC_Cap_Tropopause_Prs
    USE Set_Global_CH4_Mod, ONLY : Set_CH4
    USE PBL_Mix_Mod,        ONLY : Compute_PBL_Height
    USE Pressure_Mod,       ONLY : Set_Floating_Pressures
    USE TOMS_Mod,           ONLY : Compute_Overhead_O3
    USE UCX_Mod,            ONLY : Set_H2O_Trac

    USE ErrCode_Mod
    USE GC_Grid_Mod,        ONLY : AREA_M2
    USE HCO_Error_Mod
    USE HCO_Interface_Mod,  ONLY : SetHcoTime
    USE Pressure_Mod,       ONLY : Accept_External_Pedge
    USE State_Chm_Mod,      ONLY : IND_
    USE Time_Mod,           ONLY : Accept_External_Date_Time
    Use UnitConv_Mod,       ONLY : Convert_Spc_Units

    USE Diagnostics_Mod,    ONLY : Set_Diagnostics_EndofTimestep

    USE CMN_Size_Mod
    USE Precision_Mod
    USE Error_Mod,          ONLY : Debug_Msg
    USE GC_Grid_Mod,        ONLY : SetGridFromCtrEdges
    USE History_Mod,        ONLY : History_Write, History_Update, History_SetTime

    USE Olson_Landmap_Mod,  ONLY : Compute_Olson_Landmap_GCHP

    LOGICAL,        INTENT(IN)                 :: am_I_Root   ! Are we on root CPU?
    INTEGER,        INTENT(IN)                 :: IM          ! # of lons on this CPU
    INTEGER,        INTENT(IN)                 :: JM          ! # of lats on this CPU
    INTEGER,        INTENT(IN)                 :: LM          ! # of levs on this CPU
    INTEGER,        INTENT(IN)                 :: ID          ! Domain identifier within this CPU
    INTEGER,        INTENT(IN)                 :: nymd        ! YYYY/MM/DD @ current time
    INTEGER,        INTENT(IN)                 :: nhms        ! hh:mm:ss   @ current time
    INTEGER,        INTENT(IN)                 :: year        ! UTC year 
    INTEGER,        INTENT(IN)                 :: month       ! UTC month
    INTEGER,        INTENT(IN)                 :: day         ! UTC day
    INTEGER,        INTENT(IN)                 :: dayOfYr     ! UTC day of year
    INTEGER,        INTENT(IN)                 :: hour        ! UTC hour
    INTEGER,        INTENT(IN)                 :: minute      ! UTC minute
    INTEGER,        INTENT(IN)                 :: second      ! UTC second
    REAL(KIND_R4),  INTENT(IN)                 :: hElapsed    ! Hours elapsed in run
    TYPE(GIGC_Chunk_Operators),  INTENT(IN)    :: Operators   ! Operators to run (derived type)
    LOGICAL,        INTENT(IN)                 :: IsChemTime  ! Time for chemistry? 
    REAL(KIND_R4),  INTENT(IN)                 :: lonCtr(:,:) ! Lon centers [radians]
    REAL(KIND_R4),  INTENT(IN)                 :: latCtr(:,:) ! Lat centers [radians]
    REAL(KIND_R4), INTENT(IN)                  :: lonEdge(:,:)! Lon edges   [radians]
    REAL(KIND_R4), INTENT(IN)                  :: latEdge(:,:)! Lat edges   [radians]
    TYPE(OptInput),      INTENT(INOUT) :: Input_Opt   ! Input Options obj
    TYPE(ChmState),      INTENT(INOUT) :: State_Chm   ! Chemistry State obj
    TYPE(MetState),      INTENT(INOUT) :: State_Met   ! Meteorology State obj
    TYPE(DgnState),      INTENT(INOUT) :: State_Diag  ! Diagnostics State obj
    INTEGER,             INTENT(OUT)   :: RC          ! Return code
    REAL*8                         :: DT
    INTEGER                        :: STATUS

    LOGICAL                        :: DoConv 
    LOGICAL                        :: DoDryDep
    LOGICAL                        :: DoEmis
    LOGICAL                        :: DoTend 
    LOGICAL                        :: DoTurb 
    LOGICAL                        :: DoChem
    LOGICAL                        :: DoWetDep

    LOGICAL, SAVE                  :: FIRST = .TRUE.

    LOGICAL                        :: pUpdate

    INTEGER, SAVE                  :: NCALLS = 0

    INTEGER                        :: L = 1
    INTEGER                        :: N

    REAL*4                         :: utc


    RC = GC_SUCCESS

    AREA_M2 = State_Met%AREA_M2


    DoConv   = Input_Opt%LCONV                    ! dynamic time step
    DoDryDep = Input_Opt%LDRYD .AND. IsChemTime   ! chemistry time step
    DoEmis   = Input_Opt%LEMIS .AND. IsChemTime   ! chemistry time step
    DoTurb   = Input_Opt%LTURB                    ! dynamic time step
    DoChem   = Input_Opt%LCHEM .AND. IsChemTime   ! chemistry time step
    DoWetDep = Input_Opt%LWETD                    ! dynamic time step 

    DoConv   = DoConv .AND. Operators%Conv
    DoDryDep = DoDryDep .AND. Operators%DryDep
    DoEmis   = DoEmis .AND. Operators%Emis
    DoTurb   = DoTurb .AND. Operators%Turb
    DoChem   = DoChem .AND. Operators%Chem
    DoWetDep = DoWetDep .AND. Operators%WetDep

    DoTend = ((DoEmis .OR. DoDryDep) .AND. Operators%Tend) .AND. .NOT. DoTurb

    IF ( am_I_Root .and. NCALLS < 10 ) THEN 
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


    call GIGC_Switch_Dims( am_I_Root, ID, IM, JM, LM, &
                           lonCtr, latCtr, lonEdge, latEdge, &
                           Input_Opt, State_Met, State_Chm, State_Diag, .true., .true., RC )
    write(6, *) "GIGC_Chunk_Run: Running for in-PET grid #", ID
    write(6, *) "GIGC_Chunk_Run: Received new grid IM, JM, LM =", IM, JM, LM

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

    CALL SetGridFromCtrEdges( am_I_Root, SIZE(lonCtr, 1), SIZE(lonCtr, 2), lonCtr, latCtr, lonEdge, latEdge, RC )
    IF ( RC /= GC_SUCCESS ) THEN
      WRITE(6, *) "GIGC_Chunk_Mod Fatal: Could not set GC grid parameters from lat/lon ctr/edges"
      RETURN
    ENDIF

    where(State_Chm%Species < 0.0e0)
      State_Chm%Species = 1.0e-36
    endwhere

    LLTROP  = LM
    LLSTRAT = LM
    DO L = 1, LM
      IF(State_Met%TROPP(1, 1) .ge. State_Met%PEDGE(1, 1, L)) THEN
        LLTROP  = LLTROP - 1
      ENDIF
    ENDDO

    if(Input_Opt%LUCX) then
      LLCHEM     = LLSTRAT
      LLCHEM_FIX = LLSTRAT
    else
      LLCHEM     = LLTROP
      LLCHEM_FIX = LLTROP_FIX
    endif

    if(am_I_Root.and.NCALLS<10) write(*,*) 'GIGC_CHUNK_RUN: Diagnosed LLTROP, LLSTRAT =', LLTROP, LLSTRAT 

    if(am_I_Root.and.NCALLS<10) write(*,*) "- Running GEOS-Chem Column Code, Call", NCALLS

    IF ( DoConv .OR. DoChem .OR. DoWetDep ) THEN
       CALL SETUP_WETSCAV( am_I_Root, Input_Opt, State_Met, State_Chm, RC )
    ENDIF

    CALL Compute_Olson_Landmap_GCHP( am_I_Root, State_Met, RC )

    utc =           ( DBLE( hour ) ) + ( DBLE( minute ) / 60e+0_f8 ) + ( DBLE( second ) / 3600e+0_f8 )
    CALL Accept_External_Date_Time( am_I_Root      = am_I_Root,  &
                                    value_NYMD     = nymd,       &  
                                    value_NHMS     = nhms,       &  
                                    value_YEAR     = year,       &  
                                    value_MONTH    = month,      &  
                                    value_DAY      = day,        &  
                                    value_DAYOFYR  = dayOfYr,    &  
                                    value_HOUR     = hour,       &  
                                    value_MINUTE   = minute,     &  
                                    value_SECOND   = second,     &
                                    value_HELAPSED = hElapsed,   &
                                    value_UTC      = utc,        &
                                    RC             = RC         )

    CALL SetHcoTime ( am_I_Root, DoEmis, RC )

    CALL Accept_External_Pedge    ( am_I_Root      = am_I_Root,  &
                                    State_Met      = State_Met,  &
                                    RC             = RC         )

    State_Met%PSC2_WET = State_Met%PS2_WET
    State_Met%PSC2_DRY = State_Met%PS2_DRY
    CALL SET_FLOATING_PRESSURES( am_I_Root, State_Met, RC )
    IF ( RC /= GC_SUCCESS ) RETURN

    pUpdate = ((.not.FIRST).and.(.not.Input_Opt%LTRAN))
    CALL AirQnt( am_I_Root, Input_Opt, State_Met, State_Chm, RC, pUpdate )

    IF ( FIRST .and. Input_Opt%LSCHEM ) THEN
       CALL INIT_STRAT_CHEM( am_I_Root, Input_Opt, State_Chm, State_Met, RC )
       Minit_is_set = .true.
    ENDIF

    CALL GIGC_Cap_Tropopause_Prs  ( am_I_Root      = am_I_Root,  &
                                    IM             = IM,         &
                                    JM             = JM,         &
                                    Input_Opt      = Input_Opt,  &
                                    State_Met      = State_Met,  &
                                    RC             = RC         )

    CALL COMPUTE_PBL_HEIGHT( am_I_Root, State_Met, RC )

    CALL Convert_Spc_Units( am_I_Root, Input_Opt, State_Met, State_Chm, &
                            'kg/kg dry', RC )
    
    IF ( IND_('H2O','A') > 0 ) THEN
       CALL SET_H2O_TRAC( am_I_Root, ((.NOT. Input_Opt%LUCX) .OR.    &
                          Input_Opt%LSETH2O ), Input_Opt, State_Met, &
                          State_Chm, RC )
       IF (Input_Opt%LSETH2O) Input_Opt%LSETH2O = .FALSE.
    ENDIF

    CALL EMISSIONS_RUN( am_I_Root, Input_Opt,  State_Met,         &
                        State_Chm, State_Diag, DoEmis, 1, RC       )



    IF ( DoConv ) THEN
       if(am_I_Root.and.NCALLS<10) write(*,*) ' --- Do convection now' 

       CALL DO_CONVECTION ( am_I_Root, Input_Opt, State_Met, State_Chm, &
                            State_Diag, RC )

       if(am_I_Root.and.NCALLS<10) write(*,*) ' --- Convection done!'
    ENDIF

    IF ( DoDryDep ) THEN
       if(am_I_Root.and.NCALLS<10) THEN
          write(*,*) ' --- Do drydep now'
          write(*,*) '     Use FULL PBL: ', Input_Opt%PBL_DRYDEP
       endif

       CALL Do_DryDep( am_I_Root, Input_Opt=Input_Opt, State_Chm=State_Chm, &
                       State_Met=State_Met, State_Diag=State_Diag, RC=RC ) 

       if(am_I_Root.and.NCALLS<10) write(*,*) ' --- Drydep done!'
    ENDIF ! Do drydep

    IF ( DoEmis ) THEN
       if(am_I_Root.and.NCALLS<10) write(*,*) ' --- Do emissions now'

       CALL EMISSIONS_RUN ( am_I_Root,  Input_Opt, State_Met, State_Chm, &
                            State_Diag, DoEmis, 2, RC )

       if(am_I_Root.and.NCALLS<10) write(*,*) ' --- Emissions done!'
    ENDIF

    IF ( DoTend ) THEN
       if(am_I_Root.and.NCALLS<10) write(*,*)   &
                           ' --- Add emissions and drydep to tracers'

       DT = HcoState%TS_EMIS 

       CALL DO_TEND ( am_I_Root, Input_Opt, State_Met, State_Chm,  &
                      State_Diag, .FALSE., RC, DT=DT )

       if(am_I_Root.and.NCALLS<10) write(*,*)   &
                                 ' --- Fluxes applied to tracers!' 
    ENDIF ! Tendencies


    IF ( DoTurb ) THEN
       if(am_I_Root.and.NCALLS<10) write(*,*) ' --- Do turbulence now'

       CALL DO_MIXING ( am_I_Root, Input_Opt, State_Met, State_Chm, &
                        State_Diag, RC )

       if(am_I_Root.and.NCALLS<10) write(*,*) ' --- Turbulence done!'
    ENDIF

    IF ( .NOT. DoTurb .AND. Input_Opt%ITS_A_FULLCHEM_SIM  &
         .AND. IND_('CH4','A') > 0 ) THEN
       CALL SET_CH4 ( am_I_Root,  Input_Opt, State_Met, State_Chm, &
                      State_Diag, RC )
    ENDIF

    IF ( DoChem ) THEN
       if(am_I_Root.and.NCALLS<10) write(*,*) ' --- Do chemistry now'

       CALL COMPUTE_OVERHEAD_O3( am_I_Root, DAY, .TRUE., State_Met%TO3 )

       IF ( IND_('H2O','A') > 0 ) THEN
          CALL SET_H2O_TRAC( am_I_Root, (.not. Input_Opt%LUCX), Input_Opt, &
                             State_Met, State_Chm, RC )
       ENDIF

       CALL Do_Chemistry( am_I_Root  = am_I_Root,            & ! Root CPU?
                          Input_Opt  = Input_Opt,            & ! Input Options
                          State_Chm  = State_Chm,            & ! Chemistry State
                          State_Met  = State_Met,            & ! Met State
                          State_Diag = State_Diag,           & ! Diagn State
                          RC         = RC                   )  ! Success?

       if(am_I_Root.and.NCALLS<10) write(*,*) ' --- Chemistry done!'
    ENDIF

    IF ( DoWetDep ) THEN
       if(am_I_Root.and.NCALLS<10) write(*,*) ' --- Do wetdep now'

       CALL DO_WETDEP( am_I_Root, Input_Opt, State_Met, State_Chm,  &
                       State_Diag, RC )

       if(am_I_Root.and.NCALLS<10) write(*,*) ' --- Wetdep done!'
    ENDIF


    IF ( DoChem ) THEN
       CALL Recompute_OD( am_I_Root, Input_Opt,  State_Met,  &
                          State_Chm, State_Diag, RC         )
    ENDIF
    CALL Set_Diagnostics_EndofTimestep( am_I_Root,  Input_Opt, &
                                        State_Met,  State_Chm, &
                                        State_Diag, RC )

    IF ( State_Diag%Archive_AerMass ) THEN
       CALL Set_AerMass_Diagnostic( am_I_Root, Input_Opt,  State_Met, &
                                    State_Chm, State_Diag, RC         )
    ENDIF

    IF ( Operators%GCDiagn .and. .false. ) THEN
      IF (am_I_Root .and. NCALLS < 10) write(*,*) ' --- Do history now'

      CALL History_Update( am_I_Root, RC )

      CALL History_SetTime( am_I_Root, RC )

      CALL History_Write( am_I_Root, State_Chm%Spc_Units, RC )
      IF (am_I_Root .and. NCALLS < 10) write(*,*) ' --- History done!'
    ENDIF


    IF ( NCALLS < 10 ) NCALLS = NCALLS + 1 

    FIRST = .FALSE.

    CALL Convert_Spc_Units( am_I_Root, Input_Opt, State_Met, State_Chm, &
                            'v/v dry', RC )

    RC = GC_SUCCESS

  END SUBROUTINE GIGC_Chunk_Run
  SUBROUTINE GIGC_Chunk_Final( am_I_Root, Input_Opt,  State_Chm,             &
                               State_Met, State_Diag, RC                    )
    USE Input_Opt_Mod,    ONLY : OptInput, Cleanup_Input_Opt
    USE State_Chm_Mod,    ONLY : ChmState, Cleanup_State_Chm
    USE State_Met_Mod,    ONLY : MetState, Cleanup_State_Met
    USE State_Diag_Mod,   ONLY : DgnState, Cleanup_State_Diag
    USE HCOI_GC_MAIN_MOD, ONLY : HCOI_GC_FINAL
    USE ErrCode_Mod
    LOGICAL,        INTENT(IN)    :: am_I_Root     ! Are we on the root CPU?
    TYPE(OptInput), INTENT(INOUT) :: Input_Opt     ! Input Options object
    TYPE(ChmState), INTENT(INOUT) :: State_Chm     ! Chemistry State object
    TYPE(MetState), INTENT(INOUT) :: State_Met     ! Meteorology State object
    TYPE(DgnState), INTENT(INOUT) :: State_Diag    ! Diagnostics State object
    INTEGER,        INTENT(OUT)   :: RC            ! Success or failure

    RC = GC_SUCCESS

    CALL HCOI_GC_FINAL( am_I_Root, .FALSE., RC )
    IF ( am_I_Root ) THEN
       IF ( RC == GC_SUCCESS ) THEN
          write(*,'(a)') 'HEMCO::Finalize... OK.'
       ELSE
          write(*,'(a)') 'HEMCO::Finalize... FAILURE.'
       ENDIF
    ENDIF

    CALL Cleanup_State_Diag( am_I_Root, State_Diag, RC )
    IF ( am_I_Root ) THEN
       IF ( RC == GC_SUCCESS ) THEN
          write(*,'(a)') 'Chem::State_Diag Finalize... OK.'
       ELSE
          write(*,'(a)') 'Chem::State_Diag Finalize... FAILURE.'
       ENDIF
    ENDIF

    CALL Cleanup_State_Chm( am_I_Root, State_Chm, RC )
    IF ( am_I_Root ) THEN
       IF ( RC == GC_SUCCESS ) THEN
          write(*,'(a)') 'Chem::State_Chm Finalize... OK.'
       ELSE
          write(*,'(a)') 'Chem::State_Chm Finalize... FAILURE.'
       ENDIF
    ENDIF

    CALL Cleanup_State_Met( am_I_Root, State_Met, RC )
    IF ( am_I_Root ) THEN
       IF ( RC == GC_SUCCESS ) THEN
          write(*,'(a)') 'Chem::State_Met Finalize... OK.'
       ELSE
          write(*,'(a)') 'Chem::State_Met Finalize... FAILURE.'
       ENDIF
    ENDIF
  END SUBROUTINE GIGC_Chunk_Final
END MODULE GIGC_Chunk_Mod
