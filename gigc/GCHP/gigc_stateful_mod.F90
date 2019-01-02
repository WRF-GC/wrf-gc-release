!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!               __          _______  ______       _____  _____                 !
!               \ \        / /  __ \|  ____|     / ____|/ ____|                !
!                \ \  /\  / /| |__) | |__ ______| |  __| |                     !
!                 \ \/  \/ / |  _  /|  __|______| | |_ | |                     !
!                  \  /\  /  | | \ \| |         | |__| | |____                 !
!                   \/  \/   |_|  \_\_|          \_____|\_____|                !
!                                                                              !
!----------------------- ALPHA VERSION, v0.1 (20190101) -----------------------!
!
! WRF-GC: GEOS-Chem High Performance-powered Chemistry Add-On for WRF Model
! Developed by Haipeng Lin <linhaipeng@pku.edu.cn> (GEOS-Chem Stateful Module)
!    January 2018, Peking University, Dept of Atmospheric and Oceanic Sciences
!    Correspondence to: Tzung-May Fu <tmfu@pku.edu.cn>
!
! ALPHA INFORMATION:
!    WRF-GC Alpha (version 0.1) is experimental. Please notify the authors of
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

module GIGC_Stateful_Mod
    use ErrCode_Mod
    use Input_Opt_Mod
    use State_Met_Mod
    use State_Chm_Mod
    use State_Diag_Mod
    use DiagList_Mod, only: DgnList
    use HCO_TYPES_MOD, only: ConfigObj
    use HCO_State_Mod, only: HCO_State
    use HCOX_State_Mod, only: Ext_State
    implicit none
    private
    public :: GIGC_State_Boot
    public :: GIGC_State_Init
    public :: GIGC_State_Final
    public :: GIGC_State_Get_Status
    public :: GIGC_State_Get_Opt
    public :: GIGC_State_Get_DiagList
    public :: GIGC_State_Get_Met
    public :: GIGC_State_Get_Chm
    public :: GIGC_State_Get_Diag
    public :: GIGC_State_Get_HCO
    public :: GIGC_State_Get_HCOX
    public :: GIGC_State_Set_Opt
    public :: GIGC_State_Set_Met
    public :: GIGC_State_Set_Chm
    public :: GIGC_State_Set_Diag
    type GIGC_Stateful_Object
        integer                                :: ID = -999
        logical                                :: Init = .false.
        type(MetState)                         :: State_Met
        type(ChmState)                         :: State_Chm
        type(DgnState)                         :: State_Diag
        type(HCO_State), pointer               :: HcoState
        type(Ext_State), pointer               :: ExtState
    end type GIGC_Stateful_Object
    type(OptInput)                             :: Global_Input_Opt
    type(ConfigObj), pointer                   :: Global_HcoConfig => NULL()
    type(DgnList)                              :: Global_DiagList
#if defined ( EXTERNAL_GRID ) || defined( MODEL_ )
    integer                                    :: EXTERNAL_MAX_DOM = 8
    type(GIGC_Stateful_Object), dimension(1:8) :: States
#else
    integer                                    :: EXTERNAL_MAX_DOM = 1
    type(GIGC_Stateful_Object), dimension(1:1) :: States
#endif
    logical                                    :: Init = .false.
contains
    subroutine GIGC_State_Boot(am_I_Root, MPI_COMM, RC)
        USE GC_Environment_Mod, only: GC_Allocate_All
        USE INPUT_MOD,          only: Read_Input_File
        USE HCO_CONFIG_MOD,     only: Config_Readfile
        USE HCO_DIAGN_MOD,      only: DiagnFileOpen
        USE LINOZ_MOD,          only: Linoz_Read
        USE DiagList_Mod,       only: Init_DiagList, Print_DiagList
        USE GIGC_Mpi_Wrap,      only: GIGC_Input_Bcast, GIGC_IDT_Bcast
        logical, intent(in)           :: am_I_Root          ! Are we on the root CPU?
        integer, intent(in)           :: MPI_COMM           ! MPI Communicator #
        integer, intent(inout)        :: RC                 ! Success or failure?
        integer :: HCO_DIAGN_LUN
        RC = GC_SUCCESS
        if(Init .eq. .true.) then
            return
        endif
        call Set_Input_Opt( am_I_Root, Global_Input_Opt, RC )
#if defined ( EXTERNAL_GRID ) || defined( MODEL_ )
        Global_Input_Opt%HPC     = .true.
        Global_Input_Opt%RootCPU = am_I_Root
        Global_Input_Opt%TS_CHEM = 10   ! Chemistry timestep [min]
        Global_Input_Opt%TS_EMIS = 10   ! Chemistry timestep [min]
        Global_Input_Opt%TS_DYN  = 20   ! Dynamic   timestep [min]
        Global_Input_Opt%TS_CONV = 20   ! Dynamic   timestep [min]
#endif
        if(am_I_Root) then
            call Read_Input_File(am_I_Root, Global_Input_Opt, RC)
            if(RC /= GC_SUCCESS) return
            Global_Input_Opt%USE_O3_FROM_MET = .TRUE.
            write(6, *) '### GIGC_Stateful_Mod: Root CPU, after READ_INPUT_FILE'
            if(Global_Input_Opt%LLINOZ) then
              call Linoz_Read(am_I_Root, Global_Input_Opt, RC) 
              if(RC /= GC_SUCCESS) return
              write(6, *) '### GIGC_Stateful_Mod: Root CPU, after LINOZ_READ'
            endif
        endif
        CALL Init_DiagList(am_I_Root, "HISTORY.rc", Global_DiagList, RC)
        if(RC /= GC_SUCCESS) return
        CALL Print_DiagList(am_I_Root, Global_DiagList, RC)
        CALL GIGC_Input_Bcast(am_I_Root = am_I_Root,                           &
                              Input_Opt = Global_Input_Opt,                    &
                              RC        = RC,                                  &
                              MPI_COMM_ = MPI_COMM)
        if(RC /= GC_SUCCESS) return
        write(6, *) '### GIGC_Stateful_Mod: after GIGC_Input_Bcast'
        CALL GIGC_IDT_Bcast(am_I_Root = am_I_Root,                             &  
                            Input_Opt = Global_Input_Opt,                      &  
                            RC        = RC)
        if(RC /= GC_SUCCESS) return
        write(6, *) '### GIGC_Stateful_Mod: after GIGC_IDT_Bcast'
        call GC_Allocate_All( am_I_Root      = am_I_Root,               &
                              Input_Opt      = Global_Input_Opt,        &
                              value_I_LO     = 1,                       &
                              value_J_LO     = 1,                       &
                              value_I_HI     = 1,                       &
                              value_J_HI     = 1,                       &
                              value_IM       = 1,                       &
                              value_JM       = 1,                       &
                              value_LM       = 1,                       &
                              value_IM_WORLD = 1,                       &
                              value_JM_WORLD = 1,                       &
                              value_LM_WORLD = 1,                       &
                              RC             = RC                       )            
        if(RC /= GC_SUCCESS) return
        if(RC /= GC_SUCCESS) then 
            write(6, *) "STOP GIGC_Stateful_Mod :: Return Code /= GC_SUCCESS"
        endif
        Init = .true.
    
    end subroutine GIGC_State_Boot
    subroutine GIGC_State_Init(am_I_Root, &
                               ID, State_Met, State_Chm, State_Diag, &
                               HcoState, ExtState, &
                               RC)
        logical, intent(in)           :: am_I_Root          ! Are we on the root CPU?
        integer, intent(in)           :: ID                 ! Domain identifier ID
        type(MetState), intent(in)    :: State_Met
        type(ChmState), intent(in)    :: State_Chm
        type(DgnState), intent(in)    :: State_Diag
        type(HCO_State), pointer      :: HcoState
        type(Ext_State), pointer      :: ExtState
        integer, intent(inout)        :: RC                 ! Success or failure?
        integer                       :: D
        integer                       :: Found_Index = -1
        RC = GC_SUCCESS
        Found_Index = -1
        do D = 1, EXTERNAL_MAX_DOM
            if(States(D)%Init .and. States(D)%ID .eq. ID) then
                Found_Index = D
                exit
            endif
        enddo
        if(Found_Index .ne. -1) then
            write(6, *) "%%% GIGC_Stateful_Mod: Detected duplicate initialization of ID", ID, "%%%"
            write(6, *) "Check your code for GIGC_State_Init duplicate calls on the same ID."
            write(6, *) "You are attempting to initialize ID", ID
            write(6, *) "But domains already stored are at States%ID", States%ID
            write(6, *) "Of which Found_Index =", Found_Index, "matches ID..."
            stop
        endif
        Found_Index = -1
        do D = 1, EXTERNAL_MAX_DOM
            if(.not. States(D)%Init) then
                Found_Index = D
                exit
            endif
        enddo
        if(Found_Index .eq. -1) then
            write(6, *) "%%% GIGC_Stateful_Mod: We have no more space (max number of ext. domains is", EXTERNAL_MAX_DOM
            write(6, *) "You might need to change the headers in gigc_stateful_mod.f90"
            write(6, *) "You are attempting to initialize ID", ID
            write(6, *) "But domains already stored are at States%ID", States%ID
            stop
        endif
        States(D)%Init       = .true.
        States(D)%ID         = ID
        States(D)%State_Met  = State_Met
        States(D)%State_Chm  = State_Chm
        States(D)%State_Diag = State_Diag
        States(D)%HcoState   => HcoState
        States(D)%ExtState   => ExtState
    end subroutine GIGC_State_Init
    subroutine GIGC_State_Final(am_I_Root, &
                                ID, RC)
        logical, intent(in)           :: am_I_Root          ! Are we on the root CPU?
        integer, intent(in)           :: ID                 ! Domain identifier ID
        integer, intent(inout)        :: RC                 ! Success or failure?
        integer                       :: D
        integer                       :: Found_Index = -1
        RC = GC_SUCCESS
        Found_Index = -1
        do D = 1, EXTERNAL_MAX_DOM
            if(States(D)%Init .and. States(D)%ID .eq. ID) then
                Found_Index = D
                exit
            endif
        enddo
        if(Found_Index .eq. -1) then
            RC = GC_FAILURE
            write(6, *) "%%% GIGC_Stateful_Mod: The ID requested to be finalized in GIGC_State_Final could NOT be found"
            write(6, *) "ID requested:", ID
            return
        endif
        CALL Cleanup_State_Chm(am_I_Root, States(Found_Index)%State_Chm, RC)
        IF (am_I_Root) write(6, *) 'GIGC_Stateful_Mod State_Chm Finalize ID =', ID
        CALL Cleanup_State_Met(am_I_Root, States(Found_Index)%State_Met, RC)
        IF (am_I_Root) write(6, *) 'GIGC_Stateful_Mod State_Met Finalize ID =', ID
        CALL Cleanup_State_Diag(am_I_Root, States(Found_Index)%State_Diag, RC)
        IF (am_I_Root) write(6, *) 'GIGC_Stateful_Mod State_Diag Finalize ID =', ID
    end subroutine GIGC_State_Final
    subroutine GIGC_State_Get_Status(am_I_Root, ID, Status)
        logical, intent(in)           :: am_I_Root          ! Are we on the root CPU?
        integer, intent(in)           :: ID                 ! Domain identifier ID
        logical, intent(inout)        :: Status             ! Have we found it?
        integer                       :: D
        integer                       :: Found_Index = -1
        Found_Index = -1
        do D = 1, EXTERNAL_MAX_DOM
            if(States(D)%Init .and. States(D)%ID .eq. ID) then
                Found_Index = D
                exit
            endif
        enddo
        Status = Found_Index .ne. -1
    end subroutine GIGC_State_Get_Status
    subroutine GIGC_State_Get_Opt(am_I_Root, Input_Opt, HcoConfig)
        logical, intent(in)                :: am_I_Root          ! Are we on the root CPU?
        type(OptInput), intent(inout)      :: Input_Opt          ! Input_Opt object.
        type(ConfigObj), pointer, optional :: HcoConfig          ! HEMCO Configuration object (optional)
        if(.not. Init) then
            write(6, *) "%%% GIGC_Stateful_Mod: Not initialized but requested GIGC_State_Get_Opt. Stop."
            stop
        endif
        Input_Opt = Global_Input_Opt
        if(present(HcoConfig)) then
            HcoConfig => Global_HcoConfig
        endif
    end subroutine GIGC_State_Get_Opt
    subroutine GIGC_State_Get_DiagList(am_I_Root, DiagList)
        logical, intent(in)                :: am_I_Root          ! Are we on the root CPU?
        type(DgnList), intent(inout)      :: DiagList           ! DiagList object.
        if(.not. Init) then
            write(6, *) "%%% GIGC_Stateful_Mod: Not initialized but requested GIGC_State_Get_DiagList. Stop."
            stop
        endif
        DiagList = Global_DiagList
    end subroutine GIGC_State_Get_DiagList
    subroutine GIGC_State_Get_Met(am_I_Root, ID, State_Met, RC)
        logical, intent(in)           :: am_I_Root          ! Are we on the root CPU?
        integer, intent(in)           :: ID                 ! Domain identifier ID
        type(MetState), intent(inout) :: State_Met          ! Met state.
        integer, intent(inout)        :: RC                 ! Success or failure?
        integer                       :: D
        integer                       :: Found_Index = -1
        RC = GC_SUCCESS
        Found_Index = -1
        do D = 1, EXTERNAL_MAX_DOM
            if(States(D)%Init .and. States(D)%ID .eq. ID) then
                Found_Index = D
                exit
            endif
        enddo
        if(Found_Index .eq. -1) then
            RC = GC_FAILURE
            write(6, *) "%%% GIGC_Stateful_Mod: The ID requested could NOT be found"
            write(6, *) "ID requested:", ID
            return
        endif
        State_Met = States(Found_Index)%State_Met
    end subroutine GIGC_State_Get_Met
    subroutine GIGC_State_Get_Chm(am_I_Root, ID, State_Chm, RC)
        logical, intent(in)           :: am_I_Root          ! Are we on the root CPU?
        integer, intent(in)           :: ID                 ! Domain identifier ID
        type(ChmState), intent(inout) :: State_Chm          ! Chemistry state.
        integer, intent(inout)        :: RC                 ! Success or failure?
        integer                       :: D
        integer                       :: Found_Index = -1
        RC = GC_SUCCESS
        Found_Index = -1
        do D = 1, EXTERNAL_MAX_DOM
            if(States(D)%Init .and. States(D)%ID .eq. ID) then
                Found_Index = D
                exit
            endif
        enddo
        if(Found_Index .eq. -1) then
            RC = GC_FAILURE
            write(6, *) "%%% GIGC_Stateful_Mod: The ID requested could NOT be found"
            write(6, *) "ID requested:", ID
            return
        endif
        State_Chm = States(Found_Index)%State_Chm
    end subroutine GIGC_State_Get_Chm
    subroutine GIGC_State_Get_Diag(am_I_Root, ID, State_Diag, RC)
        logical, intent(in)           :: am_I_Root          ! Are we on the root CPU?
        integer, intent(in)           :: ID                 ! Domain identifier ID
        type(DgnState), intent(inout) :: State_Diag         ! Diagnostic state.
        integer, intent(inout)        :: RC                 ! Success or failure?
        integer                       :: D
        integer                       :: Found_Index = -1
        RC = GC_SUCCESS
        Found_Index = -1
        do D = 1, EXTERNAL_MAX_DOM
            if(States(D)%Init .and. States(D)%ID .eq. ID) then
                Found_Index = D
                exit
            endif
        enddo
        if(Found_Index .eq. -1) then
            RC = GC_FAILURE
            write(6, *) "%%% GIGC_Stateful_Mod: The ID requested could NOT be found"
            write(6, *) "ID requested:", ID
            return
        endif
        State_Diag = States(Found_Index)%State_Diag
    end subroutine GIGC_State_Get_Diag
    subroutine GIGC_State_Get_HCO(am_I_Root, ID, HcoState, RC)
        logical, intent(in)           :: am_I_Root          ! Are we on the root CPU?
        integer, intent(in)           :: ID                 ! Domain identifier ID
        type(HCO_State), pointer      :: HcoState           ! HEMCO state.
        integer, intent(inout)        :: RC                 ! Success or failure?
        integer                       :: D
        integer                       :: Found_Index = -1
        RC = GC_SUCCESS
        Found_Index = -1
        do D = 1, EXTERNAL_MAX_DOM
            if(States(D)%Init .and. States(D)%ID .eq. ID) then
                Found_Index = D
                exit
            endif
        enddo
        if(Found_Index .eq. -1) then
            RC = GC_FAILURE
            write(6, *) "%%% GIGC_Stateful_Mod: The ID requested could NOT be found"
            write(6, *) "ID requested:", ID
            return
        endif
        HcoState => States(Found_Index)%HcoState
    end subroutine GIGC_State_Get_HCO
    subroutine GIGC_State_Get_HCOX(am_I_Root, ID, ExtState, RC)
        logical, intent(in)           :: am_I_Root          ! Are we on the root CPU?
        integer, intent(in)           :: ID                 ! Domain identifier ID
        type(Ext_State), pointer      :: ExtState           ! HEMCO extensions state.
        integer, intent(inout)        :: RC                 ! Success or failure?
        integer                       :: D
        integer                       :: Found_Index = -1
        RC = GC_SUCCESS
        Found_Index = -1
        do D = 1, EXTERNAL_MAX_DOM
            if(States(D)%Init .and. States(D)%ID .eq. ID) then
                Found_Index = D
                exit
            endif
        enddo
        if(Found_Index .eq. -1) then
            RC = GC_FAILURE
            write(6, *) "%%% GIGC_Stateful_Mod: The ID requested could NOT be found"
            write(6, *) "ID requested:", ID
            return
        endif
        ExtState => States(Found_Index)%ExtState
    end subroutine GIGC_State_Get_HCOX
    subroutine GIGC_State_Set_Opt(am_I_Root, Input_Opt, HcoConfig)
        logical, intent(in)                :: am_I_Root          ! Are we on the root CPU?
        type(OptInput), intent(inout)      :: Input_Opt          ! Input_Opt object.
        type(ConfigObj), pointer, optional :: HcoConfig          ! HEMCO Configuration object (optional)
        if(.not. Init) then
            write(6, *) "%%% GIGC_Stateful_Mod: Not initialized but requested GIGC_State_Set_Opt. Stop."
            write(6, *) "(Developer Note: Are you trying to set Input_Opt manually instead? This option"
            write(6, *) " is currently unsupported, but you can edit GIGC_STATEFUL_MOD manually.)"
            stop
        endif
        Global_Input_Opt = Input_Opt
        if(present(HcoConfig)) then
            Global_HcoConfig => HcoConfig
        endif
    end subroutine GIGC_State_Set_Opt
    subroutine GIGC_State_Set_Met(am_I_Root, ID, State_Met, RC)
        logical, intent(in)           :: am_I_Root          ! Are we on the root CPU?
        integer, intent(in)           :: ID                 ! Domain identifier ID
        type(MetState), intent(inout) :: State_Met          ! Met state.
        integer, intent(inout)        :: RC                 ! Success or failure?
        integer                       :: D
        integer                       :: Found_Index = -1
        RC = GC_SUCCESS
        Found_Index = -1
        do D = 1, EXTERNAL_MAX_DOM
            if(States(D)%Init .and. States(D)%ID .eq. ID) then
                Found_Index = D
                exit
            endif
        enddo
        if(Found_Index .eq. -1) then
            RC = GC_FAILURE
            write(6, *) "%%% GIGC_Stateful_Mod: The ID requested could NOT be found"
            write(6, *) "ID requested:", ID
            return
        endif
        States(Found_Index)%State_Met = State_Met
    end subroutine GIGC_State_Set_Met
    subroutine GIGC_State_Set_Chm(am_I_Root, ID, State_Chm, RC)
        logical, intent(in)           :: am_I_Root          ! Are we on the root CPU?
        integer, intent(in)           :: ID                 ! Domain identifier ID
        type(ChmState), intent(inout) :: State_Chm          ! Chemistry state.
        integer, intent(inout)        :: RC                 ! Success or failure?
        integer                       :: D
        integer                       :: Found_Index = -1
        RC = GC_SUCCESS
        Found_Index = -1
        do D = 1, EXTERNAL_MAX_DOM
            if(States(D)%Init .and. States(D)%ID .eq. ID) then
                Found_Index = D
                exit
            endif
        enddo
        if(Found_Index .eq. -1) then
            RC = GC_FAILURE
            write(6, *) "%%% GIGC_Stateful_Mod: The ID requested could NOT be found"
            write(6, *) "ID requested:", ID
            return
        endif
        States(Found_Index)%State_Chm = State_Chm
    end subroutine GIGC_State_Set_Chm
    subroutine GIGC_State_Set_Diag(am_I_Root, ID, State_Diag, RC)
        logical, intent(in)           :: am_I_Root          ! Are we on the root CPU?
        integer, intent(in)           :: ID                 ! Domain identifier ID
        type(DgnState), intent(inout) :: State_Diag         ! Diagnostics state.
        integer, intent(inout)        :: RC                 ! Success or failure?
        integer                       :: D
        integer                       :: Found_Index = -1
        RC = GC_SUCCESS
        Found_Index = -1
        do D = 1, EXTERNAL_MAX_DOM
            if(States(D)%Init .and. States(D)%ID .eq. ID) then
                Found_Index = D
                exit
            endif
        enddo
        if(Found_Index .eq. -1) then
            RC = GC_FAILURE
            write(6, *) "%%% GIGC_Stateful_Mod: The ID requested could NOT be found"
            write(6, *) "ID requested:", ID
            return
        endif
        States(Found_Index)%State_Diag = State_Diag
    end subroutine GIGC_State_Set_Diag
end module GIGC_Stateful_Mod