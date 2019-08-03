!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!               __          _______  ______       _____  _____!               \ \        / /  __ \|  ____|     / ____|/ ____|!                \ \  /\  / /| |__) | |__ ______| |  __| |!                 \ \/  \/ / |  _  /|  __|______| | |_ | |!                  \  /\  /  | | \ \| |         | |__| | |____!                   \/  \/   |_|  \_\_|          \_____|\_____|!!----------------------- ALPHA VERSION, v0.9 (20190802) -----------------------!
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

subroutine chem_driver(grid, config_flags &
#include "dummy_new_args.inc"
)

   use module_domain, only: domain
   use module_configure
   use module_driver_constants
   use module_machine
   use module_tiles
   use module_dm
   use module_model_constants
   use module_state_description
   use module_chem_utilities
   use module_input_chem_data, only: last_chem_time, &
#if (defined(CHEM_DBG_I) && defined(CHEM_DBG_J) && defined(CHEM_DBG_K))
                                    chem_dbg, &
#endif
      get_last_gas
   use module_upper_bc_driver, only: upper_bc_driver
   use module_tropopause, only: tropopause_driver
   use GIGC_Chunk_Mod
   use PRECISION_MOD
   use Input_Opt_Mod, only: OptInput
   use State_Chm_Mod, only: ChmState
   use State_Met_Mod, only: MetState
   use State_Diag_Mod, only: DgnState
   use GIGC_Convert_State_Mod
   use GIGC_Stateful_Mod
   implicit none
   logical, external :: wrf_dm_on_monitor

   type(domain), target :: grid

# include <dummy_new_decl.inc>
# define NO_I1_OLD
   TYPE(grid_config_rec_type), INTENT(IN) :: config_flags

   integer :: ids, ide, jds, jde, kds, kde, &
              ims, ime, jms, jme, kms, kme, &
              ips, ipe, jps, jpe, kps, kpe, &
              its, ite, jts, jte, kts, kte

   integer :: GEOS_CHEM_RC

   TYPE(OptInput) :: Input_Opt   TYPE(MetState) :: State_Met   TYPE(ChmState) :: State_Chm   TYPE(DgnState) :: State_Diag
   TYPE(GIGC_Chunk_Operators) :: GIGC_Ops

   real(4), allocatable :: lonCtr(:, :)   real(4), allocatable :: latCtr(:, :)   real(4), allocatable :: lonEdge(:, :)   real(4), allocatable :: latEdge(:, :)
   integer :: IM, JM, LM

   integer :: II, JJ

   integer :: i, j, k, l, numgas, nv, n, nr, ktau
   integer :: ijulian, nymd, nhms

   logical, save :: FIRST = .TRUE.

   real, dimension(grid%sm31:grid%em31, grid%sm32:grid%em32, grid%sm33:grid%em33) :: &
      p_phy, u_phy, v_phy, t_phy, dz8w, t8w, p8w, rho, z_at_w, vvel, zmid, rh, pi_phy

   real, dimension(grid%sm31:grid%em31, grid%sm33:grid%em33) :: pbl_h

   TYPE(WRFU_TimeInterval) :: tmpTimeInterval
   real(KIND=8) :: curr_secs
   real(KIND=8) :: real_time_r8   logical      :: do_chemstep

   integer :: debug_level
   integer :: ij
   integer :: num_3d_m, ic, num_3d_c, num_3d_s

   real :: epsilc
   parameter(epsilc=1.e-30)

   real :: chem_minval

   character(len=256) :: current_date_char
   integer :: WRF_dateY, WRF_dateM, WRF_dateD, WRF_dateH, WRF_dateI, WRF_dateS

   real     :: WRF_minutesElapsed
   real(f4) :: WRF_hoursElapsed
   
   character*256 :: message_txt
   character*256 :: debug_format

   real(KIND=8), dimension(8), save ::        &
                         WRFGC_Overhead_Time, &                         WRFGC_GC_Time,       &                         WRFGC_Phys_Time,     &                         WRFGC_Diag_Time,     &                         WRF_Total_Time   real(KIND=8), save :: WRF_Time_Last_Call
   real(KIND=8)       :: WRFGC_Time_Temp_Start, WRFGC_Time_Temp_End

   intrinsic max, min

   logical :: am_I_Root
   integer :: WRF_DM_MyProc, WRF_DM_NProc, WRF_DM_Comm


   if(FIRST) then
   else
      WRF_Total_Time(grid%id) = WRF_Total_Time(grid%id) + (MPI_Wtime() - WRF_Time_Last_Call)
   endif

   if(wrf_dm_on_monitor()) then
      am_I_Root = .true.
   else
      am_I_Root = .false.
   endif

   call wrf_get_nproc(WRF_DM_NProc)
   call wrf_get_myproc(WRF_DM_MyProc)
   call wrf_get_dm_communicator(WRF_DM_Comm)

   call nl_get_debug_level(1, debug_level)
   call set_wrf_debug_level(debug_level)

   GIGC_Ops%Conv   = config_flags%gc_do_convection
   GIGC_Ops%Emis   = config_flags%gc_do_hemco
   GIGC_Ops%Tend   = .true.
   GIGC_Ops%Turb   = config_flags%gc_do_pblmix
   GIGC_Ops%Chem   = config_flags%gc_do_chemistry
   GIGC_Ops%DryDep = config_flags%gc_do_drydep
   GIGC_Ops%WetDep = config_flags%gc_do_wetdep

   GIGC_Ops%GCDiagn = config_flags%gc_do_gcdiagn

   ktau = grid%itimestep
   tmpTimeInterval = domain_get_time_since_sim_start(grid)
   curr_secs = real_time_r8(tmpTimeInterval)
   ijulian = ifix(grid%julian)

   if (ktau == 1) then
      grid%conv_ct(:, :, :, :) = 0.
      grid%chem_ct(:, :, :, :) = 0.
      grid%vmix_ct(:, :, :, :) = 0.
   endif

   do_chemstep = .false.
   if (ktau == 1) then
      do_chemstep = .true.
      grid%ktauc = 1
   else
      if ((grid%chemdt <= 0) .or. &
          (curr_secs + real(grid%dt, 8) + 0.01 >= &
           (int(curr_secs/real(grid%chemdt*60., 8) + 1, 8)*real(grid%chemdt*60., 8))) &
          ) then
         do_chemstep = .true.
         grid%ktauc = grid%ktauc + 1
         last_chem_time(grid%id) = domain_get_current_time(grid)

         call WRFU_TimeGet(last_chem_time(grid%id), &
                           YY=grid%last_chem_time_year, &
                           MM=grid%last_chem_time_month, &
                           DD=grid%last_chem_time_day, &
                           H=grid%last_chem_time_hour, &
                           M=grid%last_chem_time_minute, &
                           S=grid%last_chem_time_second)
      endif
   endif

   call get_ijk_from_grid(grid, &
                          ids, ide, jds, jde, kds, kde, &
                          ims, ime, jms, jme, kms, kme, &
                          ips, ipe, jps, jpe, kps, kpe)

   call domain_clock_get(grid, current_timestr=current_date_char)

   call domain_clock_get(grid, minutesSinceSimulationStart=WRF_minutesElapsed)
   WRF_hoursElapsed = WRF_minutesElapsed / 60.0_f4


   read (current_date_char(1:4), FMT='(I4)') WRF_dateY
   read (current_date_char(6:7), FMT='(I2)') WRF_dateM
   read (current_date_char(9:10), FMT='(I2)') WRF_dateD
   read (current_date_char(12:13), FMT='(I2)') WRF_dateH
   read (current_date_char(15:16), FMT='(I2)') WRF_dateI
   read (current_date_char(18:19), FMT='(I2)') WRF_dateS

   nymd = WRF_dateY * 10000 + WRF_dateM * 100 + WRF_dateD
   nhms = WRF_dateH * 10000 + WRF_dateI * 100 + WRF_dateS

   grid%raincv_b(:, :) = grid%raincv(:, :)

   num_3d_m = num_moist
   num_3d_c = num_chem
   num_3d_s = num_scalar
   numgas = get_last_gas(config_flags%chem_opt)

   call set_tiles(grid, ids, ide, jds, jde, ips, ipe, jps, jpe)

   chem_minval = epsilc   chem_select: select case(config_flags%chem_opt)
      case (233)
         call wrf_debug(15, 'GIGC chem_driver: chem_opt = 233. Running GEOS-Chem HP chemistry option.')
      case default
         call wrf_error_fatal("Pumpkin chem_driver: Unrecognized chem_opt. WRF-GCHP is chem_opt=233.")
   end select chem_select


   do j = jps, min(jde - 1, jpe)
      do k = kps, kpe
         do i = ips, min(ide - 1, ipe)
            vvel(i, k, j) = grid%w_2(i, k, j)
            zmid(i, k, j) = grid%z(i, k, j)
         enddo
      enddo
   enddo
   do j = jps, min(jde - 1, jpe)
      do i = ips, min(ide - 1, ipe)
         pbl_h(i, j) = grid%pblh(i, j)
      enddo
   enddo

!------------------------------------------------------------------------
! Main chemistry tile loop
!------------------------------------------------------------------------
   chem_tile_loop_1: DO ij = 1, grid%num_tiles

      its = grid%i_start(ij)
      ite = min(grid%i_end(ij), ide - 1)
      jts = grid%j_start(ij)
      jte = min(grid%j_end(ij), jde - 1)

      kts = kps
      kte = min(kpe, kde - 1)

      IM = ite - its + 1
      JM = jte - jts + 1
      LM = kte - kts + 1

      WRFGC_Time_Temp_Start = MPI_Wtime()

      CALL wrf_debug(15, 'Pumpkin chem_driver: Calling module_chem_utilities::chem_prep')

      CALL chem_prep(config_flags, &
                     grid%u_2, grid%v_2, grid%p, grid%pb, &
                     grid%alt, grid%ph_2, grid%phb, grid%t_2, &
                     moist, num_3d_m, rho, &
                     p_phy, pi_phy, u_phy, v_phy, &
                     p8w, t_phy, t8w, grid%z, z_at_w, &
                     dz8w, rh, grid%fnm, grid%fnp, &
                     ids, ide, jds, jde, kds, kde, &
                     ims, ime, jms, jme, kms, kme, &
                     its, ite, jts, jte, &
                     kps, kpe)

      WRFGC_Time_Temp_End = MPI_Wtime()
      WRFGC_Phys_Time(grid%id) = WRFGC_Phys_Time(grid%id) + (WRFGC_Time_Temp_End - WRFGC_Time_Temp_Start)

#if (defined(CHEM_DBG_I) && defined(CHEM_DBG_J) && defined(CHEM_DBG_K))
      if ((its <= CHEM_DBG_I .and. ite >= CHEM_DBG_I) .and. &
          (jts <= CHEM_DBG_J .and. jte >= CHEM_DBG_J) .and. &
          (kps <= CHEM_DBG_K .and. kpe >= CHEM_DBG_K)) then
         call wrf_debug(15, 'Pumpkin chem_driver: Calling chem_dbg before doing chemistry')
         call chem_dbg(CHEM_DBG_I, CHEM_DBG_J, CHEM_DBG_K, grid%dt, ktau, &
                       dz8w, t_phy, p_phy, rho, chem, emis_ant, &
                       ids, ide, jds, jde, kds, kde, &
                       ims, ime, jms, jme, kms, kme, &
                       its, ite, jts, jte, kts, kte, &
                       config_flags%kemit)
      end if
#endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!  CHEMISTRY DRIVER/ MECHANISM (do_timestep)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   if (do_chemstep) then
      if(am_I_Root) then
         call wrf_debug(15, 'Pumpkin chem_driver: THIS IS THE MASTER PROCESS - CALLING CHEMISTRY')
      else
         write(message_txt, *) "Pumpkin chem_driver: Calling chemistry on subprocess ", WRF_DM_MyProc, "(+1) of ", WRF_DM_NProc
         call wrf_debug(15, message_txt)
      endif

      write(message_txt, *) 'GIGC chem_driver: Before State_Met conversion IM = ', IM, ' JM = ', JM
      call wrf_debug(15, message_txt)

      write(message_txt, *) "Pumpkin chem_driver: its, ite, jts, jte, kts, kte: ", its, " ", ite, " ", jts, " ", jte, " ", kts, " ", kte
      call wrf_debug(15, message_txt)

      write(message_txt, *) "Pumpkin chem_driver: ims, ime, jms, jme, kms, kme: ", ims, " ", ime, " ", jms, " ", jme, " ", kms, " ", kme
      call wrf_debug(15, message_txt)
     
      write(message_txt, *) "Pumpkin chem_driver: ids, ide, jds, jde, kds, kde: ", ids, " ", ide, " ", jds, " ", jde, " ", kds, " ", kde
      call wrf_debug(15, message_txt)

      WRFGC_Time_Temp_Start = MPI_Wtime()

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







      call GIGC_State_Get_Opt(am_I_Root, Input_Opt)
      call GIGC_State_Get_Met(am_I_Root, grid%id, State_Met, GEOS_CHEM_RC)
      call GIGC_State_Get_Chm(am_I_Root, grid%id, State_Chm, GEOS_CHEM_RC)
      call GIGC_State_Get_Diag(am_I_Root, grid%id, State_Diag, GEOS_CHEM_RC)

      Input_Opt%HPC = .true.
      Input_Opt%RootCPU = am_I_Root
      Input_Opt%myCPU = WRF_DM_MyProc
      Input_Opt%NPES = WRF_DM_NProc
      Input_Opt%haveImpRst = .true.

      if(WRF_dateM .eq. 2 .and. WRF_dateD .eq. 15 .and. WRF_dateH .eq. 3 .and. WRF_dateI .eq. 10) then
            call wrf_debug(1, "        ___iiiii___            ")
            call wrf_debug(1, "       |           |           ")
            call wrf_debug(1, "     __|___________|__         ")
            call wrf_debug(1, "    |^^^^^^^^^^^^^^^^^|        ")
            call wrf_debug(1, "    |                 |        ")
            call wrf_debug(1, "    |                 |        ")
            call wrf_debug(1, "    ~~~~~~~~~~~~~~~~~~~        ")
            call wrf_debug(1, " hplin 02/15/1997 wrf|gchp     ")
      endif

      call wrf_debug(15, 'Pumpkin chem_driver: Diagnosing tropopause location')
      call tropopause_driver(grid%id, grid%dt, current_date_char,           &
                               t_phy, p_phy, p8w, zmid, z_at_w,             &
                               grid%tropo_lev, grid%tropo_p,  grid%tropo_z, &
                               ids, ide, jds, jde, kds, kde,                &
                               ims, ime, jms, jme, kms, kme,                &
                               its, min(grid%i_end(ij), ide - 1), jts, min(grid%j_end(ij), jde - 1), kts, kte)
      call wrf_debug(15, 'Pumpkin chem_driver: After tropopause_driver')


      call GIGC_Get_WRF(am_I_Root, config_flags, grid, num_chem, chem, num_scalar, scalar, num_moist, moist, &
                        dz8w, p8w, t8w, pi_phy, &
                        f_qc, f_qi, f_qs,   &
                        its, ite, jts, jte, &
                        ide, jde, &
                        kts, kte, &
                        ids, jds, kds, &
                        curr_secs, &
                        Input_Opt, State_Met, State_Chm, &
                        WRFGC_Phys_Time(grid%id))

      call wrf_debug(15, 'GIGC chem_driver: After GIGC_Get_WRF conversion')

      WRFGC_Time_Temp_End = MPI_Wtime()
      WRFGC_Overhead_Time(grid%id) = WRFGC_Overhead_Time(grid%id) + (WRFGC_Time_Temp_End - WRFGC_Time_Temp_Start)

      WRFGC_Time_Temp_Start = MPI_Wtime()

      if(.not. FIRST) then

         call wrf_debug(15, 'GIGC chem_driver: Before GIGC_Chunk_Run')
         call GIGC_Chunk_Run( am_I_Root  = am_I_Root,                      &                              IM         = IM,                             &                              JM         = JM,                             &                              LM         = LM,                             &                              ID         = grid%id,                        &                              nymd       = nymd,                           &                              nhms       = nhms,                           &                              year       = WRF_dateY,                      &                              month      = WRF_dateM,                      &                              day        = WRF_dateD,                      &                              dayOfYr    = ijulian,                        &                              hour       = WRF_dateH,                      &                              minute     = WRF_dateI,                      &                              second     = WRF_dateS,                      &                              hElapsed   = WRF_hoursElapsed,               &                              lonCtr     = lonCtr,                         &                              latCtr     = latCtr,                         &                              lonEdge    = lonEdge,                        &                              latEdge    = latEdge,                        &                              Input_Opt  = Input_Opt,                      &                              State_Chm  = State_Chm,                      &                              State_Met  = State_Met,                      &                              State_Diag = State_Diag,                     &                              Operators  = GIGC_Ops,                       &                              IsChemTime = .true.,                         &                              RC         = GEOS_CHEM_RC)         call wrf_debug(15, 'GIGC chem_driver: After GIGC_Chunk_Run')

      endif

      WRFGC_Time_Temp_End = MPI_Wtime()
      WRFGC_GC_Time(grid%id) = WRFGC_GC_Time(grid%id) + (WRFGC_Time_Temp_End - WRFGC_Time_Temp_Start)

      WRFGC_Time_Temp_Start = MPI_Wtime()

      call GIGC_Set_WRF(am_I_Root, config_flags, grid, num_chem, chem, num_scalar, scalar, num_moist, moist, &
                        its, ite, jts, jte, &
                        ide, jde, &
                        kts, kte, &
                        Input_Opt, State_Met, State_Chm, State_Diag)

      call wrf_debug(15, 'GIGC chem_driver: After GIGC_Set_WRF conversion')

      call GIGC_Chunk_Cleanup(config_flags, kds, kde, Input_Opt, State_Met, State_Chm)

      call GIGC_State_Set_Met(am_I_Root, grid%id, State_Met, GEOS_CHEM_RC)
      call GIGC_State_Set_Chm(am_I_Root, grid%id, State_Chm, GEOS_CHEM_RC)
      call GIGC_State_Set_Diag(am_I_Root, grid%id, State_Diag, GEOS_CHEM_RC)

      WRFGC_Time_Temp_End = MPI_Wtime()
      WRFGC_Overhead_Time(grid%id) = WRFGC_Overhead_Time(grid%id) + (WRFGC_Time_Temp_End - WRFGC_Time_Temp_Start)


      if(debug_level .ge. 1) then
         debug_format = "(2x, a30, 4x, f14.3)"
         write(6, *) "================================================="
         write(6, *) "W R F - G C   T I M E R S"
         write(6, *) "   Domain", grid%id, "(Units: s)"
         write(6, *) "================================================="
         write(6, debug_format) "=> Coupler Component          ", WRFGC_Overhead_Time(grid%id)
         write(6, debug_format) " -- Physics & Parameterization", WRFGC_Phys_Time(grid%id)
         write(6, debug_format) " -- Diagnostics Out           ", WRFGC_Diag_Time(grid%id)
         write(6, debug_format) "=> GEOS-Chem Column           ", WRFGC_GC_Time(grid%id)
         write(6, debug_format) "===> Total Chemistry          ", WRFGC_Overhead_Time(grid%id)+WRFGC_GC_Time(grid%id)
         write(6, debug_format) "=> WRF                        ", WRF_Total_Time(grid%id)
         write(6, debug_format) "===> Total WRF-GC             ", WRFGC_Overhead_Time(grid%id)+WRFGC_GC_Time(grid%id)+WRF_Total_Time(grid%id)
      endif
   end if
   do nv = 1, num_chem
      do j = jts, jte
         do i = its, ite
            chem(i, kpe, j, nv) = chem(i, kte, j, nv)
         enddo
      enddo
   enddo

!-- set upper boundary condition
   if (config_flags%have_bcs_upper) then
      call wrf_debug(15, 'Pumpkin chem_driver: set upper boundary condition')
      call upper_bc_driver(grid%id, grid%dt, current_date_char, &
                           chem, p_phy, p8w, grid%tropo_lev, &
                           ids, ide, jds, jde, kds, kde, &
                           ims, ime, jms, jme, kms, kme, &
                           its, ite, jts, jte, kts, kte)
   endif

   call wrf_debug(15, 'Pumpkin chem_driver: Exiting chem_driver (chem_tile_loop_1)')
end do chem_tile_loop_1
   FIRST = .FALSE.
   WRF_Time_Last_Call = MPI_Wtime()

end subroutine chem_driver