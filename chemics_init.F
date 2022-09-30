!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!               __          _______  ______       _____  _____                 !
!               \ \        / /  __ \|  ____|     / ____|/ ____|                !
!                \ \  /\  / /| |__) | |__ ______| |  __| |                     !
!                 \ \/  \/ / |  _  /|  __|______| | |_ | |                     !
!                  \  /\  /  | | \ \| |         | |__| | |____                 !
!                   \/  \/   |_|  \_\_|          \_____|\_____|                !
!                                                                              !
!-------------------     v2.0 RELEASE VERSION (20210223)    -------------------!
!
! WRF-GC: GEOS-Chem High Performance-powered Chemistry Add-On for WRF Model
! Developed by Haipeng Lin <hplin@g.harvard.edu>, Xu Feng <fengx7@pku.edu.cn>
!    February 2021, SUSTech, Atmospheric Chemistry and Climate Group
!    Correspondence to: Tzung-May Fu <fuzm@sustech.edu.cn>
!
! VERSION INFORMATION:
!    Please visit http://wrf.geos-chem.org for information.
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
      use module_mixactivate_wrappers, only: wrfgc_mixactivate_init
      use module_diag_aero_size_info, only: diag_aero_wrfgc_aercld_ptrs
      use GIGC_Chunk_Mod
      use HCO_TYPES_MOD, only: ConfigObj
      use Input_Opt_Mod, only: OptInput
      use State_Chm_Mod, only: ChmState
      use State_Met_Mod, only: MetState
      use State_Diag_Mod, only: DgnState
      use State_Grid_Mod
      use WRFGC_Convert_State_Mod, only: WRFGC_Set_WRF, WRFGC_IdxSetup
      use GC_Stateful_Mod, only: GIGC_States, Global_Input_Opt, Global_HcoConfig, Global_DiagList

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

      real(4), allocatable :: lonCtr(:, :)  ! Lon Centers (rad)
      real(4), allocatable :: latCtr(:, :)  ! Lat Centers (rad)
      real(4), allocatable :: lonEdge(:, :) ! Lon Edges (rad)
      real(4), allocatable :: latEdge(:, :) ! Lat Edges (rad)
      integer              :: IM, JM
      integer              :: II, JJ       ! Loop variable
      integer              :: nymdB, nhmsB, nymdE, nhmsE
      logical, allocatable :: is_aerosol(:) ! true if field is aerosol (any phase)
      real                 :: epsilc
      logical              :: Am_I_Root
      integer              :: WRF_DM_MyProc, WRF_DM_NProc, WRF_DM_Comm
      integer              :: debug_level
      call get_ijk_from_grid(grid, &
                             ids, ide, jds, jde, kds, kde, &
                             ims, ime, jms, jme, kms, kme, &
                             its, ite, jts, jte, kts, kte)

      call nl_get_debug_level(1, debug_level)
      call set_wrf_debug_level(debug_level)
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
      Global_Input_Opt%isMPI   = .true.
      Global_Input_Opt%amIRoot = Am_I_Root
      Global_Input_Opt%thisCPU = WRF_DM_MyProc
      Global_Input_Opt%numCPUs = WRF_DM_NProc
      Global_Input_Opt%MPIComm = WRF_DM_Comm

      if(debug_level .ge. 2) then
        Global_Input_Opt%LPRT    = .true. ! Print debug information? (LD70)
      endif

      if(Global_Input_Opt%KppStop .eq. .true.) then
        if(grid%id .ge. 2) then
          Global_Input_Opt%KppStop = .false.
        endif
      endif
      if(Global_Input_Opt%KppStop .eq. .false.) then
        write(6,*) "*** WRF-GC WARNING ***"
        write(6,*) "YOU ARE RUNNING WRF-GC IN MULTI-DOMAIN NESTED MODE, WHICH DISABLES"
        write(6,*) "CRITICAL INTEGRITY FAILSAFES WITHIN GEOS-CHEM KPP INTEGRATOR."
        write(6,*) "Your results may NOT be scientifically valid until you verify the run"
        write(6,*) "has completed without KPP errors. Refer to the log and State_Diag%KppError"
      endif
      if(grid%id .gt. 8 .or. grid%id .lt. 1) then
        call wrf_error_fatal("WRFGC chemics_init: Cannot run with domains < 1 or > 8. This is an arbitrary code limitation.")
      endif
      if(GIGC_States(grid%id)%Init .eq. .true.) then
        write(6, *) "WRFGC chemics_init: Domain", grid%id, "already initialized, skipping."
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

      write(message_txt, *) "   This is domain ", grid%id
      call wrf_message("*****************************************************************")
      call wrf_message("        __          _______  ______       _____  _____           ")
      call wrf_message("        \ \        / /  __ \|  ____|     / ____|/ ____|          ")
      call wrf_message("         \ \  /\  / /| |__) | |__ ______| |  __| |               ")
      call wrf_message("          \ \/  \/ / |  _  /|  __|______| | |_ | |               ")
      call wrf_message("           \  /\  /  | | \ \| |         | |__| | |____           ")
      call wrf_message("            \/  \/   |_|  \_\_|          \_____|\_____|          ")
      call wrf_message("*****************************************************************")
      call wrf_message("           WRF-GC MODEL, BUILD 2209.3X - WRF BINDINGS            ")
      call wrf_message("                 WRF-GC RELEASE 2.02    20220930                 ")
      call wrf_message(" For errors, suggestions & feedback email hplin@seas.harvard.edu ")
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

      if(.NOT. config_flags%restart) then
       do j=jts,jte
          do k=kts,kte
             do i=its,ite
                tauaer1(i,k,j) = 0.
                tauaer2(i,k,j) = 0.
                tauaer3(i,k,j) = 0.
                tauaer4(i,k,j) = 0.
                gaer1(i,k,j) = 0.
                gaer2(i,k,j) = 0.
                gaer3(i,k,j) = 0.
                gaer4(i,k,j) = 0.
                waer1(i,k,j) = 0.
                waer2(i,k,j) = 0.
                waer3(i,k,j) = 0.
                waer4(i,k,j) = 0.
                l2aer(i,k,j,1) = 0.
                l2aer(i,k,j,2) = 0.
                l2aer(i,k,j,3) = 0.
                l2aer(i,k,j,4) = 0.
                l3aer(i,k,j,1) = 0.
                l3aer(i,k,j,2) = 0.
                l3aer(i,k,j,3) = 0.
                l3aer(i,k,j,4) = 0.
                l4aer(i,k,j,1) = 0.
                l4aer(i,k,j,2) = 0.
                l4aer(i,k,j,3) = 0.
                l4aer(i,k,j,4) = 0.
                l5aer(i,k,j,1) = 0.
                l5aer(i,k,j,2) = 0.
                l5aer(i,k,j,3) = 0.
                l5aer(i,k,j,4) = 0.
                l6aer(i,k,j,1) = 0.
                l6aer(i,k,j,2) = 0.
                l6aer(i,k,j,3) = 0.
                l6aer(i,k,j,4) = 0.
                l7aer(i,k,j,1) = 0.
                l7aer(i,k,j,2) = 0.
                l7aer(i,k,j,3) = 0.
                l7aer(i,k,j,4) = 0.
                extaerlw1(i,k,j) = 0.
                extaerlw2(i,k,j) = 0.
                extaerlw3(i,k,j) = 0.
                extaerlw4(i,k,j) = 0.
                extaerlw5(i,k,j) = 0.
                extaerlw6(i,k,j) = 0.
                extaerlw7(i,k,j) = 0.
                extaerlw8(i,k,j) = 0.
                extaerlw9(i,k,j) = 0.
                extaerlw10(i,k,j) = 0.
                extaerlw11(i,k,j) = 0.
                extaerlw12(i,k,j) = 0.
                extaerlw13(i,k,j) = 0.
                extaerlw14(i,k,j) = 0.
                extaerlw15(i,k,j) = 0.
                extaerlw16(i,k,j) = 0.
                tauaerlw1(i,k,j) = 0.
                tauaerlw2(i,k,j) = 0.
                tauaerlw3(i,k,j) = 0.
                tauaerlw4(i,k,j) = 0.
                tauaerlw5(i,k,j) = 0.
                tauaerlw6(i,k,j) = 0.
                tauaerlw7(i,k,j) = 0.
                tauaerlw8(i,k,j) = 0.
                tauaerlw9(i,k,j) = 0.
                tauaerlw10(i,k,j) = 0.
                tauaerlw11(i,k,j) = 0.
                tauaerlw12(i,k,j) = 0.
                tauaerlw13(i,k,j) = 0.
                tauaerlw14(i,k,j) = 0.
                tauaerlw15(i,k,j) = 0.
                tauaerlw16(i,k,j) = 0.
             end do
          end do
       end do
      end if

      if (config_flags%have_bcs_upper) then
         call wrf_error_fatal("Pumpkin chemics_init: have_bcs_upper = .true. is not supported by this Chemistry abstraction layer.")
      endif
      if(config_flags%chem_opt == 0) then
        CALL wrf_debug(1, 'Pumpkin chemics_init: **** NO CHEMISTRY IS USED IN THIS WRF RUN ****')
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
      call diag_aero_wrfgc_aercld_ptrs(num_chem, is_aerosol, config_flags)
      kte = kte - 1

      write(message_txt, *) "Pumpkin chemics_init: its, ite, jts, jte, kts, kte: ", its, " ", ite, " ", jts, " ", jte, " ", kts, " ", kte, "domain", grid%id
      call wrf_debug(1, message_txt)

      write(message_txt, *) "Pumpkin chemics_init: ims, ime, jms, jme, kms, kme: ", ims, " ", ime, " ", jms, " ", jme, " ", kms, " ", kme, "domain", grid%id
      call wrf_debug(1, message_txt)
     
      write(message_txt, *) "Pumpkin chemics_init: ids, ide, jds, jde, kds, kde: ", ids, " ", ide, " ", jds, " ", jde, " ", kds, " ", kde, "domain", grid%id
      call wrf_debug(1, message_txt)
      if(debug_level .ge. 2) then
        write(6,*) "Pumpkin chemics_init, old indices: its, ite, jts, jte, kts, kte", oldits, oldite, oldjts, oldjte, oldkts, oldkte
        write(6,*) "Pumpkin chemics_init, old indices: ids, ide, jds, jde, kds, kde", oldids, oldide, oldjds, oldjde, oldkds, oldkde
        write(6,*) "Pumpkin chemics_init, old indices: ims, ime, jms, jme, kms, kme", oldims, oldime, oldjms, oldjme, oldkms, oldkme
      endif
      call Init_State_Grid(Global_Input_Opt, GIGC_States(grid%id)%State_Grid, GEOS_CHEM_RC)

      GIGC_States(grid%id)%State_Grid%ID          = grid%id       ! ID of this grid
      GIGC_States(grid%id)%State_Grid%NX          = ite-its+1     ! # lons   on this PET
      GIGC_States(grid%id)%State_Grid%NY          = jte-jts+1     ! # lats   on this PET
      GIGC_States(grid%id)%State_Grid%NZ          = kte-kts+1     ! # levels on this PET
      GIGC_States(grid%id)%State_Grid%GlobalNX    = ite-its+1     ! # lons   in global grid
      GIGC_States(grid%id)%State_Grid%GlobalNY    = jte-jts+1     ! # lats   in global grid
      GIGC_States(grid%id)%State_Grid%NativeNZ    = kte-kts+1     ! # levels in global grid
      GIGC_States(grid%id)%State_Grid%XMinOffset  = 1             ! X offset from global grid
      GIGC_States(grid%id)%State_Grid%XMaxOffset  = GIGC_States(grid%id)%State_Grid%NX ! X offset from global grid
      GIGC_States(grid%id)%State_Grid%YMinOffset  = 1             ! Y offset from global grid
      GIGC_States(grid%id)%State_Grid%YMaxOffset  = GIGC_States(grid%id)%State_Grid%NY ! Y offset from global grid
      GIGC_States(grid%id)%State_Grid%MaxTropLev  = GIGC_States(grid%id)%State_Grid%NZ ! # trop. levels below
      GIGC_States(grid%id)%State_Grid%MaxStratLev = GIGC_States(grid%id)%State_Grid%NZ ! # strat. levels below
      CALL GIGC_Chunk_Init( nymdB     = nymdB,                           & ! YYYYMMDD @ start of run
                            nhmsB     = nhmsB,                           & ! hhmmss   @ start of run
                            nymdE     = nymdE,                           & ! YYYYMMDD @ end of run
                            nhmsE     = nhmsE,                           & ! hhmmss   @ end of run
                            tsChem    = chemdt*60.,                      & ! Chemical timestep [s]
                            tsDyn     = chemdt*60.,                      & ! Dynamic  timestep [s]
                            lonCtr    = lonCtr,                          & ! Lon centers [radians]
                            latCtr    = latCtr,                          & ! Lat centers [radians]
                            lonEdge   = lonEdge,                         & ! Lon edges [radians]
                            latEdge   = latEdge,                         & ! Lat edges [radians]
                            Input_Opt = Global_Input_Opt,                & ! Input Options obj
                            State_Chm = GIGC_States(grid%id)%State_Chm,  & ! Chemistry State obj
                            State_Diag= GIGC_States(grid%id)%State_Diag, & ! Diagnostics State obj
                            State_Grid= GIGC_States(grid%id)%State_Grid, & ! Grid State obj
                            State_Met = GIGC_States(grid%id)%State_Met,  & ! Meteorology State obj
                            HcoConfig = Global_HcoConfig,                & ! HEMCO Configuration Object
                            RC        = GEOS_CHEM_RC,                    & ! GC Status Code (Success?)
                            MPI_COMM  = WRF_DM_Comm,                     & ! MPI Communicator
                            ID        = grid%id)                           ! Domain identifier, this PET

      write(message_txt, *) "GIGC chemics_init: Return from GIGC_Chunk_Init return code ", GEOS_CHEM_RC
      call wrf_debug(1, message_txt)

      write(message_txt, *) "GIGC chemics_init: GIGC_Chunk_Init completed for ID #", grid%id
      call wrf_debug(1, message_txt)
      Global_Input_Opt%isMPI   = .true.
      Global_Input_Opt%amIRoot = Am_I_Root
      Global_Input_Opt%thisCPU = WRF_DM_MyProc
      Global_Input_Opt%numCPUs = WRF_DM_NProc
      Global_Input_Opt%MPIComm = WRF_DM_Comm

      if(debug_level .ge. 2) then
        Global_Input_Opt%LPRT    = .true. ! Print debug information? (LD70)
      endif
      write(6, *) "GIGC chemics_init: Following information regarding Input_Opt% (for all doms)"
      write(6, *) "After domain", grid%id
      write(6, *) "HPC: ", Global_Input_Opt%isMPI
      write(6, *) "LPRT: ", Global_Input_Opt%LPRT
      write(6, *) "DoConv: ", Global_Input_Opt%LCONV
      write(6, *) "DoDryDep: ", Global_Input_Opt%LDRYD
      write(6, *) "DoEmis: ", Global_Input_Opt%LEMIS
      write(6, *) "DoChem: ", Global_Input_Opt%LCHEM
      write(6, *) "DoTurb: ", Global_Input_Opt%LTURB
      write(6, *) "DoWetDep: ", Global_Input_Opt%LWETD
      call WRFGC_IdxSetup(Am_I_Root)
      if(config_flags%chem_in_opt .eq. 0) then
        call WRFGC_Set_WRF(Am_I_Root, config_flags, grid, num_chem, chem, num_scalar, scalar, 1, moist_dummy, & 
                           its, ite, jts, jte, &
                           ide, jde, &
                           kts, kte, &
                           Global_Input_Opt, &
                           GIGC_States(grid%id)%State_Met, &
                           GIGC_States(grid%id)%State_Chm, &
                           GIGC_States(grid%id)%State_Diag)
      endif
      do nv = 1, num_chem
         do j = jts, jte
            do i = its, ite
               chem(i, kde, j, nv) = chem(i, kde - 1, j, nv)
            enddo
         enddo
      enddo
      if ( (.not.config_flags%restart) .and. (config_flags%progn > 0) ) then
         call wrfgc_mixactivate_init(                   &
          config_flags, chem, scalar,                 &
         chem_in_opt,                                &
         ims, ime, jms, jme, kms, kme,               &
         its, ite, jts, jte, kts, kte)
      endif
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
      deallocate(lonCtr )
      deallocate(latCtr )
      deallocate(lonEdge)
      deallocate(latEdge)
END SUBROUTINE chem_init