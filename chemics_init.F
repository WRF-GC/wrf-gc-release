!WRF:MODEL_LAYER:INITIALIZATION
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
! Developed by Haipeng Lin <linhaipeng@pku.edu.cn>, Xu Feng <fengx7@pku.edu.cn>
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

subroutine chem_init(id, chem, emis_ant, scalar, dt, bioemdt, photdt, chemdt, stepbioe, &
                    stepphot, stepchem, stepfirepl, plumerisefire_frq, z_at_w, xlat, xlong, &
                    g, aerwrf, config_flags, grid, alt, t, p, CONVFAC, ttday, tcosz, julday, gmt, &
                    tauaer1, tauaer2, tauaer3, tauaer4, &
                    gaer1, gaer2, gaer3, gaer4, &
                    waer1, waer2, waer3, waer4, &
                    l2aer, l3aer, l4aer, l5aer, l6aer, l7aer, &
                    extaerlw1, extaerlw2, extaerlw3, extaerlw4, &
                    extaerlw5, extaerlw6, extaerlw7, extaerlw8, &
                    extaerlw9, extaerlw10, extaerlw11, extaerlw12, &
                    extaerlw13, extaerlw14, extaerlw15, extaerlw16, &
                    tauaerlw1, tauaerlw2, tauaerlw3, tauaerlw4, &
                    tauaerlw5, tauaerlw6, tauaerlw7, tauaerlw8, &
                    tauaerlw9, tauaerlw10, tauaerlw11, tauaerlw12, &
                    tauaerlw13, tauaerlw14, tauaerlw15, tauaerlw16, &
                    dgnum4d, dgnumwet4d, dgnum_a1, dgnum_a2, dgnum_a3, &
                    dgnumwet_a1, dgnumwet_a2, dgnumwet_a3, &
                    pm2_5_dry, pm2_5_water, pm2_5_dry_ec, &
                    last_chem_time_year, last_chem_time_month, &
                    last_chem_time_day, last_chem_time_hour, &
                    last_chem_time_minute, last_chem_time_second, &
                    chem_in_opt, kemit, num_vert_mix, &
                    oldids, oldide, oldjds, oldjde, oldkds, oldkde, &
                    oldims, oldime, oldjms, oldjme, oldkms, oldkme, &
                    oldits, oldite, oldjts, oldjte, oldkts, oldkte)

      use module_domain
      use module_configure
      use module_state_description
      use module_model_constants, only: DEGRAD
      use module_dm
      use module_cam_support, only: gas_pcnst_modal_aero, gas_pcnst_modal_aero_pos
      use module_tropopause, only: tropopause_init
      use module_input_chem_data, only: get_last_gas, last_chem_time, setup_gasprofile_maps
      use GIGC_Chunk_Mod
      use HCO_TYPES_MOD, only: ConfigObj
      use Input_Opt_Mod, only: OptInput
      use State_Chm_Mod, only: ChmState
      use State_Met_Mod, only: MetState
      use State_Diag_Mod, only: DgnState
      use GIGC_Convert_State_Mod, only: GIGC_Set_WRF, GIGC_IdxSetup
      use GIGC_Stateful_Mod, only: GIGC_State_Get_Status
      implicit none
      logical, external :: wrf_dm_on_monitor

      real, intent(in) :: bioemdt, photdt, chemdt, dt, gmt
      integer, intent(in) :: plumerisefire_frq
      integer, intent(in) :: chem_in_opt
      integer, intent(inout) :: num_vert_mix
      integer, intent(in) :: id, julday, kemit, &
                             last_chem_time_year, &
                             last_chem_time_month, &
                             last_chem_time_day, &
                             last_chem_time_hour, &
                             last_chem_time_minute, &
                             last_chem_time_second, &
                             oldids, oldide, oldjds, oldjde, oldkds, oldkde, &
                             oldims, oldime, oldjms, oldjme, oldkms, oldkme, &
                             oldits, oldite, oldjts, oldjte, oldkts, oldkte

      integer :: ids, ide, jds, jde, kds, kde, &
                 ims, ime, jms, jme, kms, kme, &
                 its, ite, jts, jte, kts, kte
                 
      integer, intent(out) :: stepbioe, stepphot, stepchem, stepfirepl
      type(grid_config_rec_type), intent(in) :: config_flags
      type(domain), intent(inout) :: grid

      real, dimension(grid%sm31:grid%em31, grid%sm32:grid%em32, grid%sm33:grid%em33), &
         intent(INOUT) :: pm2_5_dry, pm2_5_water, pm2_5_dry_ec, &
         tauaer1, tauaer2, tauaer3, tauaer4, &
         extaerlw1, extaerlw2, extaerlw3, extaerlw4, &
         extaerlw5, extaerlw6, extaerlw7, extaerlw8, &
         extaerlw9, extaerlw10, extaerlw11, extaerlw12, &
         extaerlw13, extaerlw14, extaerlw15, extaerlw16, &
         tauaerlw1, tauaerlw2, tauaerlw3, tauaerlw4, &
         tauaerlw5, tauaerlw6, tauaerlw7, tauaerlw8, &
         tauaerlw9, tauaerlw10, tauaerlw11, tauaerlw12, &
         tauaerlw13, tauaerlw14, tauaerlw15, tauaerlw16, &
         gaer1, gaer2, gaer3, gaer4, &
         waer1, waer2, waer3, waer4

      real, dimension(grid%sm31:grid%em31, grid%sm32:grid%em32, grid%sm33:grid%em33, 3), intent(INOUT) :: dgnum4d, dgnumwet4d
      real, dimension(grid%sm31:grid%em31, grid%sm32:grid%em32, grid%sm33:grid%em33), intent(INOUT) :: dgnum_a1, dgnum_a2, dgnum_a3, dgnumwet_a1, dgnumwet_a2, dgnumwet_a3


      real, dimension(grid%sm31:grid%em31, grid%sm32:grid%em32, grid%sm33:grid%em33, 1:4), &
         intent(INOUT) :: l2aer, l3aer, l4aer, l5aer, l6aer, l7aer
      real, dimension(grid%sm31:grid%em31, grid%sm32:grid%em32, grid%sm33:grid%em33), intent(IN) :: z_at_w, t, p, alt, convfac
      real, dimension(grid%sm31:grid%em31, grid%sm32:grid%em32, grid%sm33:grid%em33, num_chem), intent(INOUT) :: chem
      real, dimension(grid%sm31:grid%em31, grid%sm32:grid%em32, grid%sm33:grid%em33, 1) :: moist_dummy
      real, dimension(grid%sm31:grid%em31, grid%sm32:grid%em32, grid%sm33:grid%em33, num_emis_ant), intent(INOUT) :: emis_ant
      real, dimension(grid%sm31:grid%em31, grid%sm32:grid%em32, grid%sm33:grid%em33, num_scalar), intent(INOUT) :: scalar
      real, dimension(grid%sm31:grid%em31, grid%sm32:grid%em32, grid%sm33:grid%em33), intent(INOUT) :: aerwrf
      real, dimension(grid%sm31:grid%em31, grid%sm33:grid%em33), intent(INOUT) :: ttday, tcosz, xlat, xlong
      real, intent(IN) :: g

      character*256 :: mminlu_loc
      character*256 :: message_txt
      type(WRFU_TimeInterval) :: tmpTimeInterval
      integer :: i, j, k, l, numgas, n, kk, nv, gigc_ptr
      integer :: GEOS_CHEM_RC
      logical :: IsInit
      type(OptInput)       :: Input_Opt  
      type(MetState)       :: State_Met  
      type(ChmState)       :: State_Chm  
      type(DgnState)       :: State_Diag 
      type(ConfigObj), pointer :: HcoConfig
      real(4), allocatable :: lonCtr(:, :) 
      real(4), allocatable :: latCtr(:, :) 
      real(4), allocatable :: lonEdge(:, :)
      real(4), allocatable :: latEdge(:, :)
      integer              :: IM, JM
      integer              :: II, JJ
      integer              :: nymdB, nhmsB, nymdE, nhmsE
      logical, allocatable :: is_aerosol(:) ! true if field is aerosol (any phase)
      real                 :: epsilc
      logical              :: Am_I_Root
      integer              :: WRF_DM_MyProc, WRF_DM_NProc, WRF_DM_Comm

      call get_ijk_from_grid(grid, &
                             ids, ide, jds, jde, kds, kde, &
                             ims, ime, jms, jme, kms, kme, &
                             its, ite, jts, jte, kts, kte)

      ite = min(ite, ide - 1)
      jte = min(jte, jde - 1)

      IM = ite - its + 1
      JM = jte - jts + 1

      GEOS_CHEM_RC = 0

      if(wrf_dm_on_monitor()) then
         Am_I_Root = .true.
      else
         Am_I_Root = .false.
      endif

      call wrf_get_nproc(WRF_DM_NProc)
      call wrf_get_myproc(WRF_DM_MyProc)
      call wrf_get_dm_communicator(WRF_DM_Comm)

      Input_Opt%HPC     = .true.
      Input_Opt%RootCPU = Am_I_Root
      Input_Opt%myCPU   = WRF_DM_MyProc
      Input_Opt%NPES    = WRF_DM_NProc
      Input_Opt%MPICOMM = WRF_DM_Comm

      Input_Opt%LPRT    = .true. ! Print debug information? (LD70)
      Input_Opt%MAX_FAM = 250

      call GIGC_State_Get_Status(Am_I_Root, grid%id, IsInit)
      if(IsInit) then
        call wrf_message("GIGC chemics_init: Domain already initialized, skipping.")
        return
      endif
      allocate(lonCtr (IM,   JM))
      allocate(latCtr (IM,   JM))
      allocate(lonEdge(IM+1, JM))
      allocate(latEdge(IM,   JM+1))
      do j = jts, jte
         do i = its, ite
            II = i - its + 1
            JJ = j - jts + 1
            lonCtr(II, JJ) = (grid%XLONG(i, j)) * DEGRAD
            latCtr(II, JJ) = (grid%XLAT (i, j)) * DEGRAD
            lonEdge(II, JJ) = grid%XLONG_U(i, j) * DEGRAD
            if(i .eq. ite) then
              lonEdge(II+1, JJ) = grid%XLONG_U(i+1, j) * DEGRAD
            endif
            latEdge(II, JJ) = grid%XLAT_V(i, j) * DEGRAD
            if(j .eq. jte) then
              latEdge(II, JJ+1) = grid%XLAT_V(i, j+1) * DEGRAD
            endif
         enddo
      enddo
      nymdB = config_flags%start_year * 10000 + config_flags%start_month * 100 + config_flags%start_day
      nhmsB = config_flags%start_hour * 10000 + config_flags%start_minute * 100 + config_flags%start_second
      nymdE = config_flags%end_year * 10000 + config_flags%end_month * 100 + config_flags%end_day
      nhmsE = config_flags%end_hour * 10000 + config_flags%end_minute * 100 + config_flags%end_second

#ifdef CHEM_DBG_I
      call print_chem_species_index(config_flags%chem_opt)
#endif

      write(message_txt, *) "   This is domain ", grid%id
      call wrf_message("*****************************************************************")
      call wrf_message("  __          _______  ______       _____  _____ _    _ _____    ")
      call wrf_message("  \ \        / /  __ \|  ____|     / ____|/ ____| |  | |  __ \   ")
      call wrf_message("   \ \  /\  / /| |__) | |__ ______| |  __| |    | |__| | |__) |  ")
      call wrf_message("    \ \/  \/ / |  _  /|  __|______| | |_ | |    |  __  |  ___/   ")
      call wrf_message("     \  /\  /  | | \ \| |         | |__| | |____| |  | | |       ")
      call wrf_message("      \/  \/   |_|  \_\_|          \_____|\_____|_|  |_|_|       ")
      call wrf_message("*****************************************************************")
      call wrf_message("         WRF-GC MODEL, VERSION 1901.01 ALPHA - WRF BINDINGS      ")
      call wrf_message(" For errors, suggestions & feedback email linhaipeng@pku.edu.cn  ")
      call wrf_message("*****************************************************************")
      call wrf_message(message_txt)
      numgas = get_last_gas(config_flags%chem_opt)
      IF (config_flags%chem_opt == 0 .AND. config_flags%aer_ra_feedback .NE. 0) THEN
         call wrf_error_fatal("chemics_init: If no chemistry, aer_ra_feedback must be 0")
      ENDIF
      CALL nl_get_mminlu(1, mminlu_loc)
      IF (trim(mminlu_loc) /= 'USGS' .and. trim(mminlu_loc) /= 'MODIFIED_IGBP_MODIS_NOAH') THEN
         print *, mminlu_loc
         call wrf_error_fatal("chemics_init: Chemistry routines require USGS or MODIS_NOAH land use maps.")
      ELSE
         IF (trim(mminlu_loc) == 'USGS' .and. grid%num_land_cat <= 23) THEN
            call wrf_error_fatal("chemics_init: USGS land use map should have 24 or more categories.")
         ELSEIF (trim(mminlu_loc) == 'MODIFIED_IGBP_MODIS_NOAH' .and. grid%num_land_cat <= 19) THEN
            call wrf_error_fatal("chemics_init: MODIS_NOAH land use map should have 20 or more categories.")
         ENDIF
      ENDIF
      if (config_flags%progn > 0) then
         call wrf_error_fatal("Pumpkin chemics_init: progn > 0 option is not supported by this Chemistry abstraction layer.")
      endif
      if (config_flags%have_bcs_upper) then
         call wrf_error_fatal("Pumpkin chemics_init: have_bcs_upper = .true. is not supported by this Chemistry abstraction layer.")
      endif
      if(config_flags%chem_opt == 0) then
        CALL wrf_debug(15, 'Pumpkin chemics_init: **** NO CHEMISTRY IS USED IN THIS WRF RUN ****')
        return
      endif
      num_vert_mix = 0
      IF (config_flags%bl_pbl_physics == ACMPBLSCHEME) THEN
         mix_select: select case(config_flags%chem_opt)
            case (radm2sorg_aq, radm2sorg_aqchem, racmsorg_aq)

               num_vert_mix = numgas
            case default
               num_vert_mix = num_chem
         end select mix_select

         if (num_vert_mix .gt. config_flags%ndepvel) then
            write (message_txt, '(A30,2(I8,2x))') 'chem_init: num_vert_mix and ndepvel ', num_vert_mix, config_flags%ndepvel
            call wrf_message(trim(message_txt))
            call wrf_error_fatal("Pumpkin chemics_init: FATAL - num_vert_mix > ndepvel ")
         endif
      endif

      stepbioe = nint(bioemdt*60./dt)
      stepphot = nint(photdt*60./dt)
      stepchem = nint(chemdt*60./dt)
      stepfirepl = nint(plumerisefire_frq*60/dt)
      stepbioe = max(stepbioe, 1)
      stepphot = max(stepphot, 1)
      stepchem = max(stepchem, 1)
      stepfirepl = max(stepfirepl, 1)
      call wrf_debug(15, 'Pumpkin chemics_init: Starting chemistry init')
      call setup_gasprofile_maps(config_flags%chem_opt, numgas)

      if (.not. allocated(is_aerosol)) then
         allocate (is_aerosol(num_chem))
      else
         if (size(is_aerosol) /= num_chem) &
            call wrf_error_fatal("The number of chemistry species has changed between nests. Use same chem_opt for all domains.")
      end if
      kte = kte - 1
      write(message_txt, *) "Pumpkin chemics_init: its, ite, jts, jte, kts, kte: ", its, " ", ite, " ", jts, " ", jte, " ", kts, " ", kte, "domain", grid%id
      call wrf_debug(15, message_txt)

      write(message_txt, *) "Pumpkin chemics_init: ims, ime, jms, jme, kms, kme: ", ims, " ", ime, " ", jms, " ", jme, " ", kms, " ", kme, "domain", grid%id
      call wrf_debug(15, message_txt)
     
      write(message_txt, *) "Pumpkin chemics_init: ids, ide, jds, jde, kds, kde: ", ids, " ", ide, " ", jds, " ", jde, " ", kds, " ", kde, "domain", grid%id
      call wrf_debug(15, message_txt)
      call GIGC_Chunk_Init(am_I_Root = Am_I_Root,           &
                           I_LO      = 1,                   &
                           J_LO      = 1,                   &
                           I_HI      = ite-its+1,           &
                           J_HI      = jte-jts+1,           &
                           IM        = ite-its+1,           &
                           JM        = jte-jts+1,           &
                           LM        = kte-kts+1,           &
                           ID        = grid%id,             &
                           IM_WORLD  = ite-its+1,           &
                           JM_WORLD  = jte-jts+1,           &
                           LM_WORLD  = kte-kts+1,           &
                           nymdB     = nymdB,               &
                           nhmsB     = nhmsB,               &
                           nymdE     = nymdE,               &
                           nhmsE     = nhmsE,               &
                           tsChem    = chemdt*60.,          &
                           tsDyn     = chemdt*60.,          &
                           lonCtr    = lonCtr,              &
                           latCtr    = latCtr,              &
                           lonEdge   = lonEdge,             &
                           latEdge   = latEdge,             &
                           myPET     = WRF_DM_MyProc,       &
                           Input_Opt = Input_Opt,           &
                           State_Chm = State_Chm,           &
                           State_Diag= State_Diag,          &
                           State_Met = State_Met,           &
                           HcoConfig = HcoConfig,           &
                           RC        = GEOS_CHEM_RC,        &
                           MPI_COMM  = WRF_DM_Comm)          

      write(message_txt, *) "GIGC chemics_init: Return from GIGC_Chunk_Init return code ", GEOS_CHEM_RC
      call wrf_debug(1, message_txt)

      write(message_txt, *) "GIGC chemics_init: GIGC_Chunk_Init completed for ID #", grid%id
      call wrf_debug(1, message_txt)
      write(6, *) "GIGC chemics_init: Following information regarding Input_Opt% (for all doms)"
      write(6, *) "After domain", grid%id
      write(6, *) "HPC: ", Input_Opt%HPC
      write(6, *) "LPRT: ", Input_Opt%LPRT
      write(6, *) "DoConv: ", Input_Opt%LCONV
      write(6, *) "DoDryDep: ", Input_Opt%LDRYD
      write(6, *) "DoEmis: ", Input_Opt%LEMIS
      write(6, *) "DoChem: ", Input_Opt%LCHEM
      write(6, *) "DoTurb: ", Input_Opt%LTURB
      write(6, *) "DoWetDep: ", Input_Opt%LWETD
      call GIGC_IdxSetup(Am_I_Root)
      if(config_flags%have_bcs_chem .eq. .false.) then
        call GIGC_Set_WRF(Am_I_Root, config_flags, grid, num_chem, chem, num_scalar, scalar, 1, moist_dummy, & 
                          its, ite, jts, jte, &
                          ide, jde, &
                          kts, kte, &
                          Input_Opt, State_Met, State_Chm, State_Diag)
      endif
      do nv = 1, num_chem
         do j = jts, jte
            do i = its, ite
               chem(i, kde, j, nv) = chem(i, kde - 1, j, nv)
            enddo
         enddo
      enddo
      if (config_flags%restart) then
         call wrf_debug(15, "Setting last_chem_time from restart file")
         call WRFU_TimeSet(last_chem_time(id), &
                           YY=last_chem_time_year, &
                           MM=last_chem_time_month, &
                           DD=last_chem_time_day, &
                           H=last_chem_time_hour, &
                           M=last_chem_time_minute, &
                           S=last_chem_time_second)
      else
         call wrf_debug(15, "Setting last_chem_time to model start time-dt")
         call WRFU_TimeIntervalSet(tmpTimeInterval, s_=real(dt, 8))
         last_chem_time(id) = domain_get_current_time(grid) - tmpTimeInterval
      end if

      call tropopause_init( id, xlat, xlong, config_flags,   &
                            ids, ide, jds,jde, kds,kde,      &
                            ims, ime, jms,jme, kms,kme,      &
                            its, ite, jts,jte, kts,kte       )

      gas_pcnst_modal_aero_pos = max(1, gas_pcnst_modal_aero) ! BSINGH - Fix for non-MAM simulations
END SUBROUTINE chem_init