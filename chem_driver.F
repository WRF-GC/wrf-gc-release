!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!               __          _______  ______       _____  _____                 !
!               \ \        / /  __ \|  ____|     / ____|/ ____|                !
!                \ \  /\  / /| |__) | |__ ______| |  __| |                     !
!                 \ \/  \/ / |  _  /|  __|______| | |_ | |                     !
!                  \  /\  /  | | \ \| |         | |__| | |____                 !
!                   \/  \/   |_|  \_\_|          \_____|\_____|                !
!                                                                              !
!------------------------------------------------------------------------------!
!         WRF-GC - two-way, nested coupled meteorology-chemistry model         !
!                             v3.0.0 (April 2024)                              !
!------------------------------------------------------------------------------!
!  Formerly known as WRF-GCHP, WRF-GIGC. GCHP-powered chemistry add-on for WRF !
!  WRF & GCHP are (c) their original authors.                                  !
!  WRF-GC coupling layer (c) Haipeng Lin, Xu Feng, Tzung-May Fu                !
!  Developed at Atmospheric Chemistry and Climate Group, Peking U./SUSTech     !
!    by Haipeng Lin <hplin@seas.harvard.edu>, Xu Feng <xfeng@g.harvard.edu>    !
!    April 2018                                                                !
!                                                                              !
!  Atmospheric Chemistry and Climate Group, SUSTech: https://atmoschem.org.cn/ !
!  WRF-GC Documentation: https://wrfgc.readthedocs.io/                         !
!------------------------------------------------------------------------------!
!  Chemistry driver developed based on Project "Pumpkin":                      !
!    Abstracted Bindings for Chemistry-to-WRF.                                 !
!  (c) 2018 Haipeng Lin, Peking U.                                             !
!  This chemical interface (chem/) is written after comprehensive study of the !
!  original chem_driver.F from WRF-Chem v3.6.1, (c) respective authors.        !
!------------------------------------------------------------------------------!
!  Refer to the chem_driver.F chemistry driver for release history             !
!------------------------------------------------------------------------------!

!------------------------------------------------------------------------------
!                               The WRF-GC Project                            !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: chem_driver.F
!
! !DESCRIPTION: Subroutine chem_driver is the entry point for WRF to call chem
!  routines and adapted to enter (GI)GC entry point.
!\\
!\\
! !INTERFACE:
subroutine chem_driver(grid, config_flags &
#include "dummy_new_args.inc"
)

!
! !USES:
!
   !! WRF Interfaces
   use module_domain, only: domain
   use module_configure
   use module_driver_constants
   use module_machine
   use module_tiles
   use module_dm
   use module_model_constants
   use module_state_description

   !! Chemistry stub interfaces
   ! module_chem_utilities: Provides chem_prep, which passes grid meteo information to chemistry
   ! to local variables below.
   use module_chem_utilities

   ! module_input_chem_data: get chemistry information, namely last_chem_time (last runtime of chemistry),
   ! get_last_gas, last gas species index #
   use module_input_chem_data, only: last_chem_time, get_last_gas

   ! upper_bc_driver: provides upper boundary conditions
   use module_upper_bc_driver, only: upper_bc_driver

   ! module_tropopause: diagnoses tropopause location
   use module_tropopause, only: tropopause_driver
   
   ! module_diag_aero_size_info: diagnoses aerosol size distributions
   use module_diag_aero_size_info, only: diag_aero_size_info
   
   ! wrfgc_mixactivate: aerosol activation
   use module_mixactivate_wrappers

   !! GIGC: Grid Independent GEOS-Chem
   ! GIGC_Chunk_Mod provides the GEOS-Chem Column Model Code
   use GIGC_Chunk_Mod
   use PRECISION_MOD ! Variable Precision code (fp) parameter

   ! Input_Opt_Mod provides the Input_Opt% type (OptInput) struct
   use Input_Opt_Mod, only: OptInput
   use State_Chm_Mod, only: ChmState
   use State_Met_Mod, only: MetState
   use State_Diag_Mod, only: DgnState
   use State_Grid_Mod, only: GrdState

   ! Use our conversion module for GEOS-Chem code
   use WRFGC_Convert_State_Mod

   ! Use our stateful module for GEOS-Chem code
   use GC_Stateful_Mod

   ! End of module use declarations
   implicit none

!
! !INPUT/OUTPUT PARAMETERS:
!
   ! provides information on WRF DMPAR parallel (MPI) process information
   ! namely whether this is the master process, using wrf_dm_on_monitor()
   logical, external :: wrf_dm_on_monitor

   !! Input data.
   type(domain), target :: grid

   !! Definitions of dummy arguments to solve
# include <dummy_new_decl.inc>
# define NO_I1_OLD

   ! namestep configuration flags (config_flags%chem_opt, for example)
   TYPE(grid_config_rec_type), INTENT(IN) :: config_flags
!
! !REMARKS:
!  This code is based upon the original interface for WRF-Chem and adapted to become
!  an abstraction layer for generic chemistry in WRF. See the "Pumpkin" project.
!
! !REVISION HISTORY:
!  30 May 2018 - H.P. Lin   - Added ProTeX headers
!  14 Apr 2022 - H.P. Lin   - Add fix for fractional timesteps breaking HEMCO daily update
!  25 Sep 2022 - H.P. Lin   - Update for specialty simulations CO2 and CH4
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
   !-- ids           start index for i in domain
   !-- ide           end index for i in domain
   !-- jds           start index for j in domain
   !-- jde           end index for j in domain
   !-- kds           start index for k in domain
   !-- kde           end index for k in domain

   !-- ims           start index for i in memory
   !-- ime           end index for i in memory
   !-- jms           start index for j in memory
   !-- jme           end index for j in memory
   !-- kms           start index for k in memory
   !-- kme           end index for k in memory

   !-- its           start index for i in tile
   !-- ite           end index for i in tile
   !-- jts           start index for j in tile
   !-- jte           end index for j in tile
   !-- kts           start index for k in tile
   !-- kte           end index for k in tile

   !-- ips           start index for i in patch
   !-- ipe           end index for i in patch
   !-- jps           start index for j in patch
   !-- jpe           end index for j in patch
   integer :: ids, ide, jds, jde, kds, kde, &
              ims, ime, jms, jme, kms, kme, &
              ips, ipe, jps, jpe, kps, kpe, &
              its, ite, jts, jte, kts, kte

   !! GEOS-Chem Integration Variables
   ! GEOS_CHEM_RC (RC in GEOS-Chem speak) is the return result code for all GIGC Subroutines.
   ! 0 = Success, > 0 is error
   integer :: GEOS_CHEM_RC

   ! Operator Options
   TYPE(GIGC_Chunk_Operators) :: GIGC_Ops

   ! For passing grid location (latlon) information
   real(4), allocatable :: lonCtr(:, :)  ! Lon Centers (rad)
   real(4), allocatable :: latCtr(:, :)  ! Lat Centers (rad)
   real(4), allocatable :: lonEdge(:, :) ! Lon Edges (rad)
   real(4), allocatable :: latEdge(:, :) ! Lat Edges (rad)

   ! For tiling information
   ! # of lons, lats, levs (I, J, L) in-PET
   integer :: IM, JM, LM

   ! Loop variables
   integer :: II, JJ

   ! Local Scalars
   integer :: i, j, k, l, numgas, nv, n, nr, ktau
   integer :: ijulian, nymd, nhms

   ! Check if first time running - if yes, skip because PBL Height is not computed in the first WRF timestep
   ! Note: Arbitrary hard-coded 8-domain limit here (hplin, 5/21/20)
   logical, save, dimension(1:8) :: FIRST = .TRUE.

   ! Memorder XZY: %s-em31 is I, %s-em32 is K, %s-em33 is J
   ! Local Variables computed for chemistry usage.
   !
   ! To pass physics (meteo) variable to CHEMISTRY, edit this to create new local variables
   ! (sized as ims:ime, kms:kme, jms:jme)
   !
   ! Modify the call to module_chem_utilities::chem_prep,
   ! and modify the calls to external_driver as necessary.
   real, dimension(grid%sm31:grid%em31, grid%sm32:grid%em32, grid%sm33:grid%em33) :: &
      p_phy, u_phy, v_phy, t_phy, dz8w, t8w, p8w, rho, z_at_w, vvel, zmid, rh, pi_phy, rri

   ! PBL Height, dimension is (I x J) obviously...
   real, dimension(grid%sm31:grid%em31, grid%sm33:grid%em33) :: pbl_h

   ! Variables for calculating adaptive timesteps.
   TYPE(WRFU_TimeInterval) :: tmpTimeInterval
   real(KIND=8) :: curr_secs
   real(KIND=8) :: real_time_r8 ! ext. function in adapt_timestep_em.F
   logical      :: do_chemstep

   integer :: debug_level
   integer :: ij
   integer :: num_3d_m, ic, num_3d_c, num_3d_s

   ! From module_ctrans_grell - used for setting a minimum range on chemistry species.
   real :: epsilc
   parameter(epsilc=1.e-30)

   real :: chem_minval

   character(len=256) :: current_date_char
   integer :: WRF_dateY, WRF_dateM, WRF_dateD, WRF_dateH, WRF_dateI, WRF_dateS
   ! Y-M-D H:I:S

   ! Calculate last time chem_driver was called, in order to
   ! do an extra hemco advance if this is a new day
   integer, save, dimension(1:8) :: WRF_lastNewDay
   logical  :: WRFGC_isAdvanceNewDayMidnight

   real     :: WRF_minutesElapsed, WRF_dateUTC
   real(f4) :: WRF_hoursElapsed
   
   ! For optical driver
   logical  :: haveaer
   ! For activate driver
   logical  :: is_gc
   integer  :: nbin_o
   
   ! For debug
   character*256 :: message_txt
   character*256 :: debug_format

   ! For timing and scalability tests
   ! Note: Cannot use config_flags%max_dom in a SAVE statement, so we just save the outermost 8
   !       domains' information.
   real(KIND=8), dimension(8), save ::        &
                         WRFGC_Overhead_Time, &          ! WRF-GC Coupling Overhead Time, Generic [s]
                         WRFGC_GC_Time,       &          ! GEOS-Chem Column Code Time [s]
                         WRFGC_Phys_Time,     &          ! WRF-GC Coupling, Physics/Parameterization [s]
                         WRFGC_Diag_Time,     &          ! WRF-GC Coupling, Diagnostics [s]
                         WRF_Total_Time                  ! Total WRF Time (Estimate from MPI) [s]
   real(KIND=8), save :: WRF_Time_Last_Call
   real(KIND=8)       :: WRFGC_Time_Temp_Start, WRFGC_Time_Temp_End

   intrinsic max, min

   ! WRF DM (MPI) Parallel Information - is master process?
   logical :: am_I_Root
   integer :: WRF_DM_MyProc, WRF_DM_NProc, WRF_DM_Comm

   ! ---------------- Code begins here ---------------- !

   ! Measure time WRF between chem_driver calls
   if(FIRST(1)) then
      ! Do not measure first timestep (initialization time is not taken into account)
   else
      ! The time spent in WRF is essentially the time between chem_driver calls.
      ! Though we cannot distinguish between domain times...
      !
      ! ...wrf does stuff...
      ! in chem_driver
      ! out chem_driver --- measure time 1
      ! ...wrf does stuf...
      ! in chem_driver  --- measure time 2 => time in WRF is 2 - 1
      ! out chem_driver --- measure time 1
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

   ! Get debug level
   call nl_get_debug_level(1, debug_level)
   call set_wrf_debug_level(debug_level)

   ! Setup the adaptive timestep for the chem routines. Most of this follows
   ! what is in solve_em, except for the call to adjust time_step.
   !
   ! The necessary variables exist for the EM core and using the adaptive
   ! techniques will work even with a constant time step. In fact, they
   ! prevent issues with restarts and changed time steps. So, we will
   ! always use them with the EM core.
   ktau = grid%itimestep
   tmpTimeInterval = domain_get_time_since_sim_start(grid)
   curr_secs = real_time_r8(tmpTimeInterval)
   ijulian = ifix(grid%julian)

   ! initializing diagnostics and macros
   if (ktau == 1) then
      grid%conv_ct(:, :, :, :) = 0.
      grid%chem_ct(:, :, :, :) = 0.
      grid%vmix_ct(:, :, :, :) = 0.
   endif

   ! Check whether to do chemistry steps (do_chemstep logical)
   do_chemstep = .false.
   if (ktau == 1) then
      do_chemstep = .true.
      grid%ktauc = 1
   else if ( config_flags%restart .and. FIRST(grid%id) ) then
      do_chemstep = .true.
   else
      if ((grid%chemdt <= 0) .or. &
          (curr_secs + real(grid%dt, 8) + 0.01 >= &
           (int(curr_secs/real(grid%chemdt*60., 8) + 1, 8)*real(grid%chemdt*60., 8))) &
          ) then
         do_chemstep = .true.
         grid%ktauc = grid%ktauc + 1
         last_chem_time(grid%id) = domain_get_current_time(grid)

         ! Don't be fooled, this is for ESMF, it's to Get data into WRFU's Ancient ESMF bindings,
         ! not get from grid! hplin 4/25/18
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

   ! Note this gets MINUTES first then convert it to a seconds number...
   call domain_clock_get(grid, minutesSinceSimulationStart=WRF_minutesElapsed)
   WRF_hoursElapsed = WRF_minutesElapsed / 60.0_f4

   ! write(6, *) "Pumpkin chem_driver: current_date_char = ", current_date_char
   ! 2018-02-22_00:00:00

   read (current_date_char(1:4), FMT='(I4)') WRF_dateY
   read (current_date_char(6:7), FMT='(I2)') WRF_dateM
   read (current_date_char(9:10), FMT='(I2)') WRF_dateD
   read (current_date_char(12:13), FMT='(I2)') WRF_dateH
   read (current_date_char(15:16), FMT='(I2)') WRF_dateI
   read (current_date_char(18:19), FMT='(I2)') WRF_dateS

   ! Compute nymd (YYYYMMDD) & nhms (hhmmss) for GEOS-Chem Column Code
   nymd = WRF_dateY * 10000 + WRF_dateM * 100 + WRF_dateD
   nhms = WRF_dateH * 10000 + WRF_dateI * 100 + WRF_dateS

   ! UTC hours
   WRF_dateUTC = WRF_dateH + WRF_dateI/60.0 + WRF_dateS/3600.0

   grid%raincv_b(:, :) = grid%raincv(:, :)

   ! Number of moist variables
   num_3d_m = num_moist
   num_3d_c = num_chem
   num_3d_s = num_scalar
   numgas = get_last_gas(config_flags%chem_opt)

   ! Compute these starting and stopping locations for each tile and number of tiles.
   call set_tiles(grid, ids, ide, jds, jde, ips, ipe, jps, jpe)

   chem_minval = epsilc ! chem_minval can be case dependant and set below...
   chem_select: select case(config_flags%chem_opt)
      case (233)
         call wrf_debug(15, 'GEOS-Chem chem_driver: chem_opt = 233. Running GEOS-Chem HP chemistry option.')
         haveaer = .false.
      case default
         call wrf_error_fatal("Pumpkin chem_driver: Unrecognized chem_opt. WRF-GC is chem_opt=233.")
   end select chem_select

   ! Some leftovers from original WRF chem_select. Don't touch unless you know what you are doing (I don't)
   ! $ hplin 20180327 $
   do j = jps, min(jde - 1, jpe)
      do k = kps, kpe
         do i = ips, min(ide - 1, ipe)
            vvel(i, k, j) = grid%w_2(i, k, j)
            zmid(i, k, j) = grid%z(i, k, j)
         enddo
      enddo
   enddo
   do j = jps, min(jde - 1, jpe)
      do k = kps, min(kde - 1, kpe)
         do i = ips, min(ide - 1, ipe) 
            rri(i, k, j) = grid%alt(i, k, j)
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
      ! hplin 5/7/2018: I have no idea why it loops only until ide-1, but
      ! I have removed this so allocs to XMID (GC) in chemics_init are not
      ! broken afterwards.
      !
      ! hplin 5/15/2018: Learn everyday. WRF does everything for a reason.
      ! Do NOT touch IDE, JDE -- ALWAYS allocate ONE less when you encounter
      ! the end boundaries. You have been warned.

      its = grid%i_start(ij)
      ite = min(grid%i_end(ij), ide - 1)
      ! ite = grid%i_end(ij)
      jts = grid%j_start(ij)
      ! jte = grid%j_end(ij)
      jte = min(grid%j_end(ij), jde - 1)

      kts = kps
      kte = min(kpe, kde - 1)

      ! Compute IM, JM, LM (lon, lat, lev) on this particular run instance
      IM = ite - its + 1
      JM = jte - jts + 1
      LM = kte - kts + 1

      ! Start WRF-GC Coupling Physics Timing (Timer-ID: chem_driver_3)
      ! Timer: WRFGC_Phys_Time
      WRFGC_Time_Temp_Start = MPI_Wtime()

      CALL wrf_debug(15, 'Pumpkin chem_driver: Calling module_chem_utilities::chem_prep')

      ! chem_prep @ module_chem_utilities.f
      ! The chem_prep routine reads in the grid & configuration flags (grid%...)
      ! And does conversions onto memorder XZY local variables, so they can be used by chemistry.
      ! This is considered for timer purposes part of the physics/parameterization process.
      CALL chem_prep(config_flags, &
                     grid%u_2, grid%v_2, grid%p, grid%pb, grid%alt, grid%ph_2, &
                     grid%phb, grid%t_2, moist, num_3d_m, &
                     rho, p_phy, pi_phy, &
                     u_phy, v_phy, p8w, t_phy, t8w, &
                     grid%z, z_at_w, dz8w, rh, &
                     grid%fnm, grid%fnp, &
                     ids, ide, jds, jde, kds, kde, &
                     ims, ime, jms, jme, kms, kme, &
                     its, ite, jts, jte, &
                     kps, kpe)

      ! End WRF-GC Coupling Physics Timing (Timer-ID: chem_driver_3)
      ! Timer: WRFGC_Phys_Time
      WRFGC_Time_Temp_End = MPI_Wtime()
      WRFGC_Phys_Time(grid%id) = WRFGC_Phys_Time(grid%id) + (WRFGC_Time_Temp_End - WRFGC_Time_Temp_Start)
      ! End Timer chem_driver_3

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

      ! Start WRF-GC Coupling Overhead Timing (Timer-ID: chem_driver_0)
      ! Timer: WRFGC_Overhead_Time
      WRFGC_Time_Temp_Start = MPI_Wtime()

      ! Allocate array for lon/lat centers for passing lat/lon geo. information to GIGC. (hplin, 4/24/18)
      ! Added grid box edge information (xfeng, 11/11/18)
      allocate(lonCtr (IM,   JM))
      allocate(latCtr (IM,   JM))
      allocate(lonEdge(IM+1, JM))
      allocate(latEdge(IM,   JM+1))

      ! For every grid box populate with relevant info...
      do j = jts, jte
         do i = its, ite
            ! Loop indexes in GC
            II = i - its + 1
            JJ = j - jts + 1

            ! Prepare lat/lon grid box center information in [rad]
            ! Grab lat lon directly from grid%XLAT, grid%XLONG [deg] (DEGRAD is piconst/180. in model_constants)
            lonCtr(II, JJ) = (grid%XLONG(i, j)) * DEGRAD
            latCtr(II, JJ) = (grid%XLAT (i, j)) * DEGRAD

            ! Grab grid%XLONG_U, grid%XLAT_V for lonEdge, latEdge
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

      ! Debugging only -- output longitude, latitudes for grid boxes, first in WRF format (so edges)
      ! then in GEOS-Chem format as diagnosed by WRF-GC coupling chem_driver
      ! write(6, *) "=============================================================="
      ! write(6, *) "D I A G N O S T I C   W R F   G R I D   C O O R D I N A T E S"
      ! write(6, *) "=============================================================="

      ! write(6, *) "its, ite, jts, jte", its, ite, jts, jte
      ! write(6, *) "IM, JM", IM, JM

      ! do i = its, ite+1
      !    do j = jts, jte+1
      !       write(6, *) "(", grid%XLONG(i, j), ",", grid%XLAT(i, j), ")"
      !    enddo
      !    write(6, *) " "
      ! enddo

      ! write(6, *) "===  G E O S - C H E M    L O N / L A T    C E N T E R S  ==="

      ! do i = 1, IM+1
      !    do j = 1, JM+1
      !       write(6, *) "(", lonCtr(i, j), ",", latCtr(i, j), ")"
      !    enddo
      !    write(6, *) " "
      ! enddo

      ! write(6, *) "GIGC chem_driver lonCtr debug:", lonCtr
      ! write(6, *) "GIGC chem_driver latCtr debug:", latCtr
      ! End debugging

      ! Retrieve stateful information from GIGC_Stateful_Mod (passing in grid#)
      ! No longer required -- directly tap into Global_Input_Opt, or the appropriate
      ! GIGC_States(grid%id)%

      ! Pass some HPC Information to Input_Opt...
      Global_Input_Opt%isMPI   = .true.
      Global_Input_Opt%amIRoot = Am_I_Root
      Global_Input_Opt%thisCPU = WRF_DM_MyProc
      Global_Input_Opt%numCPUs = WRF_DM_NProc
      Global_Input_Opt%MPIComm = WRF_DM_Comm

      ! write(6, *) 'GIGC chem_driver: NYMD, NHMS debug: nymd =', nymd, 'nhms =', nhms, 'ijulian =', ijulian
      if(WRF_dateM .eq. 2 .and. WRF_dateD .eq. 15 .and. WRF_dateH .eq. 3 .and. WRF_dateI .eq. 10) then
            call wrf_debug(1, "        ___iiiii___            ")
            call wrf_debug(1, "       |           |           ")
            call wrf_debug(1, "     __|___________|__         ")
            call wrf_debug(1, "    |^^^^^^^^^^^^^^^^^|        ")
            call wrf_debug(1, "    |                 |        ")
            call wrf_debug(1, "    |                 |        ")
            call wrf_debug(1, "    ~~~~~~~~~~~~~~~~~~~        ")
            call wrf_debug(1, "   hplin 02/15/1997 wrf|gc     ")
      endif

      ! Diagnose the tropopause location using WRF-Chem's tropopause_driver.
      call wrf_debug(15, 'Pumpkin chem_driver: Diagnosing tropopause location')
      call tropopause_driver(grid%id, grid%dt, current_date_char,           &
                               t_phy, p_phy, p8w, zmid, z_at_w,             &
                               grid%tropo_lev, grid%tropo_p,  grid%tropo_z, &
                               ids, ide, jds, jde, kds, kde,                &
                               ims, ime, jms, jme, kms, kme,                &
                               its, min(grid%i_end(ij), ide - 1), jts, min(grid%j_end(ij), jde - 1), kts, kte)
      call wrf_debug(15, 'Pumpkin chem_driver: After tropopause_driver')

      ! WRITE(6, *) "%%%% REGISTERED STATE_CHM SPECIES DUMP BEFORE CHUNK-CNV %%%%"
      ! do N = 1, MIN(State_Chm%nSpecies, 8)
      !     WRITE(6, *) "N:", N, "Name:", State_Chm%SpcData(N)%Info%Name
      !     ! WRITE(6, *) "Full Name:", State_Chm%SpcData(N)%Info%FullName
      !     WRITE(6, *) "Value at PET (1,1,1):", State_Chm%Species(1, 1, 1, N)
      !     WRITE(6, *) "Value at PET (2,1,1):", State_Chm%Species(2, 1, 1, N)
      !     WRITE(6, *) "Value at PET (7,1,1):", State_Chm%Species(7, 1, 1, N)
      ! enddo

      ! Run the conversion routines for GEOS-Chem column chemistry code
      ! Now passes full tile information for conversion (hplin, 5/7/2018)
      call WRFGC_Get_WRF(am_I_Root, config_flags, grid, num_chem, chem, num_scalar, scalar, num_moist, moist, &
                        dz8w, p8w, pi_phy, &
                        f_qc, f_qi, f_qndrop, grid%warm_rain, &
                        its, ite, jts, jte, &
                        ide, jde, &
                        kts, kte, &
                        ids, jds, kds, &
                        curr_secs, &
                        Global_Input_Opt, &
                        GIGC_States(grid%id)%State_Met,  &
                        GIGC_States(grid%id)%State_Chm,  &
                        GIGC_States(grid%id)%State_Grid, &
                        WRFGC_Phys_Time(grid%id))

      call wrf_debug(15, 'GIGC chem_driver: After WRFGC_Get_WRF conversion')

      ! End WRF-GC Coupling Overhead Timing (Timer-ID: chem_driver_0)
      ! Timer: WRFGC_Overhead_Time
      WRFGC_Time_Temp_End = MPI_Wtime()
      WRFGC_Overhead_Time(grid%id) = WRFGC_Overhead_Time(grid%id) + (WRFGC_Time_Temp_End - WRFGC_Time_Temp_Start)
      ! End Timer chem_driver_0

      ! Start WRF-GC GC Column Timing (Timer-ID: chem_driver_1)
      ! Timer: WRFGC_Overhead_Time
      WRFGC_Time_Temp_Start = MPI_Wtime()

      if((.not. FIRST(grid%id)) .and. (WRF_lastNewDay(grid%id) .ne. WRF_dateD)) then
         ! 10/30 00:00 -- no ops. special case (FIRST(grid%id)). update WRF_lastNewDay = 30
         ! 10/30 00:08 -- ...
         ! 10/30 23:58 -- ...
         ! 10/31 00:08 -- WRF_lastNewDay was 30, but now WRF_dateD is 31. additional call

         ! fix for hemco day list not updating correctly at fractional time steps.
         ! reported by xfeng, April 2022. hplin 4/14/22
         
         GIGC_Ops%Conv   = .false.
         GIGC_Ops%Emis   = .false.
         GIGC_Ops%Tend   = .false.
         GIGC_Ops%Turb   = .false.
         GIGC_Ops%Chem   = .false.
         GIGC_Ops%DryDep = .false.
         GIGC_Ops%WetDep = .false.
         GIGC_Ops%Rad    = .false.
         GIGC_Ops%GCDiagn = .false.

         ! do a dry-run
         ! warning: this nhms and WRF_dateUTC are inconsistent with the minute, second = 0, 0
         ! input here. this is good enough (checked) for SetHcoTime, but do not try
         ! to run anything else within here or things will be screwed up
         CALL GIGC_Chunk_Run(                                         &
                      ID         = grid%id,                           & ! Domain identifier, this PET
                      nymd       = nymd,                              & ! Current YYYYMMDD
                      nhms       = nhms,                              & ! Current hhmmss
                      year       = WRF_dateY,                         & ! Current year
                      month      = WRF_dateM,                         & ! Current month
                      day        = WRF_dateD,                         & ! Current day
                      dayOfYr    = ijulian,                           & ! Current doy
                      hour       = WRF_dateH,                         & ! Current hour
                      minute     = 0,                                 & ! Current minute
                      second     = 0,                                 & ! Current second
                      utc        = WRF_dateUTC,                       & ! Current UTC [hrs]
                      hElapsed   = WRF_hoursElapsed,                  & ! Elapsed hours
                      Input_Opt  = Global_Input_Opt,                  & ! Input Options
                      State_Chm  = GIGC_States(grid%id)%State_Chm,    & ! Chemistry State
                      State_Diag = GIGC_States(grid%id)%State_Diag,   & ! Diagnostics State
                      State_Grid = GIGC_States(grid%id)%State_Grid,   & ! Grid State
                      State_Met  = GIGC_States(grid%id)%State_Met,    & ! Meteorology State
                      lonCtr     = lonCtr,                            & ! Column Lon Centers [rad]
                      latCtr     = latCtr,                            & ! Column Lat Centers [rad]
                      lonEdge    = lonEdge,                           & ! Column Lon Edges [rad]
                      latEdge    = latEdge,                           & ! Column Lat Edges [rad]
                      Operators  = GIGC_Ops,                          & ! Run operators
                      IsChemTime = .false.,                           & ! Time for chem?
                      IsRadTime  = .false.,                           & ! Time for rad? (stub)
                      RC         = GEOS_CHEM_RC     )                   ! Success or fail?
      endif

      ! If first timestep, some meteorology data isn't quite ready and needs some computing
      ! before running GEOS-Chem. Hence skip first timestep.

      ! First timestep needs to run GIGC_Chunk_Run in order to compute emission quantities.
      ! Otherwise HEMCO will fail as it receives a off-centered clock.
      if(.not. FIRST(grid%id)) then
         ! Options to run operators in GIGC_Chunk_Run.
         ! It now accepts configuration in the WRF namelist for each individual process,
         ! so you can turn off turbulence if it causes spurious results
         ! (hplin, 8/15/18)
         GIGC_Ops%Conv   = config_flags%gc_do_convection
         GIGC_Ops%Emis   = config_flags%gc_do_hemco
         GIGC_Ops%Tend   = .not. config_flags%gc_do_pblmix
         GIGC_Ops%Turb   = config_flags%gc_do_pblmix
         GIGC_Ops%Chem   = config_flags%gc_do_chemistry
         GIGC_Ops%DryDep = config_flags%gc_do_drydep
         GIGC_Ops%WetDep = config_flags%gc_do_wetdep
         GIGC_Ops%Rad    = .false.

         GIGC_Ops%GCDiagn = config_flags%gc_do_gcdiagn
      else
         GIGC_Ops%Conv   = .false.
         GIGC_Ops%Emis   = .false.
         GIGC_Ops%Tend   = .false.
         GIGC_Ops%Turb   = .false.
         GIGC_Ops%Chem   = .false.
         GIGC_Ops%DryDep = .false.
         GIGC_Ops%WetDep = .false.
         GIGC_Ops%Rad    = .false.

         GIGC_Ops%GCDiagn = .false.
      endif

      call wrf_debug(15, 'GIGC chem_driver: Before GIGC_Chunk_Run')

      ! Run the GEOS-Chem column chemistry code for the given phase
      ! $ hplin 4/3/2018 1419 $
      ! $ hplin 28/12/2019 1759 $
      ! $ hplin 17/3/2020 1004 $
      ! $ hplin 14/4/2022 1401 $
      CALL GIGC_Chunk_Run(  ID         = grid%id,                           & ! Domain identifier, this PET
                            nymd       = nymd,                              & ! Current YYYYMMDD
                            nhms       = nhms,                              & ! Current hhmmss
                            year       = WRF_dateY,                         & ! Current year
                            month      = WRF_dateM,                         & ! Current month
                            day        = WRF_dateD,                         & ! Current day
                            dayOfYr    = ijulian,                           & ! Current doy
                            hour       = WRF_dateH,                         & ! Current hour
                            minute     = WRF_dateI,                         & ! Current minute
                            second     = WRF_dateS,                         & ! Current second
                            utc        = WRF_dateUTC,                       & ! Current UTC [hrs]
                            hElapsed   = WRF_hoursElapsed,                  & ! Elapsed hours
                            Input_Opt  = Global_Input_Opt,                  & ! Input Options
                            State_Chm  = GIGC_States(grid%id)%State_Chm,    & ! Chemistry State
                            State_Diag = GIGC_States(grid%id)%State_Diag,   & ! Diagnostics State
                            State_Grid = GIGC_States(grid%id)%State_Grid,   & ! Grid State
                            State_Met  = GIGC_States(grid%id)%State_Met,    & ! Meteorology State
                            lonCtr     = lonCtr,                            & ! Column Lon Centers [rad]
                            latCtr     = latCtr,                            & ! Column Lat Centers [rad]
                            lonEdge    = lonEdge,                           & ! Column Lon Edges [rad]
                            latEdge    = latEdge,                           & ! Column Lat Edges [rad]
                            Operators  = GIGC_Ops,                          & ! Run operators
                            IsChemTime = .true.,                            & ! Time for chem?
                            IsRadTime  = .false.,                           & ! Time for rad? (stub)
                            RC         = GEOS_CHEM_RC     )                   ! Success or fail?

      call wrf_debug(15, 'GIGC chem_driver: After GIGC_Chunk_Run')

      ! End WRF-GC GC Column Timing (Timer-ID: chem_driver_1)
      ! Timer: WRFGC_Overhead_Time
      WRFGC_Time_Temp_End = MPI_Wtime()
      WRFGC_GC_Time(grid%id) = WRFGC_GC_Time(grid%id) + (WRFGC_Time_Temp_End - WRFGC_Time_Temp_Start)
      ! End Timer chem_driver_1

      ! Start WRF-GC Coupling Overhead Timing (Timer-ID: chem_driver_2)
      ! Timer: WRFGC_Overhead_Time
      WRFGC_Time_Temp_Start = MPI_Wtime()

      ! Update the last day WRF-GC chem_driver was called (hplin, 4/14/22)
      WRF_lastNewDay = WRF_dateD

      ! Run the conversion routines from GEOS-Chem column chemistry code
      call WRFGC_Set_WRF(am_I_Root, config_flags, grid, num_chem, chem, num_scalar, scalar, num_moist, moist, &
                        its, ite, jts, jte, &
                        ide, jde, &
                        kts, kte, &
                        Global_Input_Opt, &
                        GIGC_States(grid%id)%State_Grid, &
                        GIGC_States(grid%id)%State_Met,  &
                        GIGC_States(grid%id)%State_Chm,  &
                        GIGC_States(grid%id)%State_Diag)

      call wrf_debug(15, 'GIGC chem_driver: After WRFGC_Set_WRF conversion')
     
      
      !-- Calculate the aerosol optical properties
      if(.not. FIRST(grid%id)) then
      if(config_flags%aer_ra_feedback .eq. 1) then
         call optical_driver(grid%id,curr_secs,grid%dt,config_flags,haveaer, &
                             is_gc, chem,dz8w,rri,rh, &
                             grid%tauaer1,grid%tauaer2,grid%tauaer3,grid%tauaer4, &
                             grid%extaer1,grid%extaer2,grid%extaer3,grid%extaer4, &
                             grid%gaer1,grid%gaer2,grid%gaer3,grid%gaer4, &
                             grid%waer1,grid%waer2,grid%waer3,grid%waer4, &
                             grid%bscoef1,grid%bscoef2,grid%bscoef3,grid%bscoef4, &
                             grid%l2aer,grid%l3aer,grid%l4aer,grid%l5aer,grid%l6aer,grid%l7aer, &
                             grid%extaerlw1,grid%extaerlw2,grid%extaerlw3,grid%extaerlw4,grid%extaerlw5, &
                             grid%extaerlw6,grid%extaerlw7,grid%extaerlw8,grid%extaerlw9,grid%extaerlw10, &
                             grid%extaerlw11,grid%extaerlw12,grid%extaerlw13,grid%extaerlw14,grid%extaerlw15, &
                             grid%extaerlw16, &
                             grid%tauaerlw1,grid%tauaerlw2,grid%tauaerlw3,grid%tauaerlw4,grid%tauaerlw5, &
                             grid%tauaerlw6,grid%tauaerlw7,grid%tauaerlw8,grid%tauaerlw9,grid%tauaerlw10, &
                             grid%tauaerlw11,grid%tauaerlw12,grid%tauaerlw13,grid%tauaerlw14,grid%tauaerlw15, &
                             grid%tauaerlw16, &
                             ids,ide, jds,jde, kds,kde, &
                             ims,ime, jms,jme, kms,kme, &
                             its,ite, jts,jte, kts,kte)
      endif
      endif


      !--- Diagnose the aerosol size distributions and calculate the aerosol activation
      if(.not. FIRST(grid%id)) then
      if(config_flags%aer_cu_feedback .eq. 1) then
         nbin_o = 4
         is_gc = .TRUE.
         call diag_aero_size_info(nbin_o, chem, num_chem, rh, is_gc,       &
                                     ids, ide, jds, jde, kds, kde,         &
                                     ims, ime, jms, jme, kms, kme,         &
                                     its, ite, jts, jte, kts, kte          )

         call mixactivate_driver(grid%id, ktau, grid%dt, config_flags, &
                                 rho, t_phy, vvel, grid%cldfra, grid%cldfra_old,  &
                                 zmid, dz8w, p8w, t8w, grid%exch_h, &
                                 moist, scalar, chem, &
                                 grid%ccn1, grid%ccn2, grid%ccn3, grid%ccn4, &
                                 grid%ccn5, grid%ccn6, grid%qndropsource, &
                                 ids, ide, jds, jde, kds, kde, &
                                 ims, ime, jms, jme, kms, kme, &
                                 its, ite, jts, jte, kts, kte  )
      endif
      endif        

      ! End WRF-GC Coupling Overhead Timing (Timer-ID: chem_driver_2)
      ! Timer: WRFGC_Overhead_Time
      WRFGC_Time_Temp_End = MPI_Wtime()
      WRFGC_Overhead_Time(grid%id) = WRFGC_Overhead_Time(grid%id) + (WRFGC_Time_Temp_End - WRFGC_Time_Temp_Start)
      ! End Timer chem_driver_2


      ! Output diagnostic debug output for timers (do_chemstep)
      ! Debug output for timers --
      if(debug_level .ge. 1) then
         ! Output timer information
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
   end if ! Chemistry time step check

   ! Fill top level to prevent spurious interpolation results (no extrapolation)
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

! No longer first call
   FIRST(grid%id) = .FALSE.

   ! Cleanup
   ! deallocate(lonCtr )
   ! deallocate(latCtr )
   ! deallocate(lonEdge)
   ! deallocate(latEdge)

   ! For WRF time measurement
   WRF_Time_Last_Call = MPI_Wtime()

end subroutine chem_driver

!
! While I am waiting for unit tests to complete, I would like to bury a message here at
! this lonely corner of the code. (hplin, 8/11/18)
!
! I would like to take this opportunity to thank those along the long and tedious,
! no doubt frustrating journey of the WRF-GC (formerly "WRF-GCHP") project.
!
! This project was envisioned by my group's leader, Prof. Tzung-May Fu, and assigned to
! me and my fellow collaborator Feng Xu, in the winter of 2017. Without Prof. Fu's
! tremendous support and trust that we could make it, this project would never come into
! fruition as it has today. Thank you.
!
! I'd also like to thank Feng Xu, who has provided the day-to-day scientific support
! that was crucial as I looked at WRF & GEOS-Chem from a programmer's perspective only.
! Thanks for bearing with me as I broke, built, and broke things again and again,
! gradually building the core architecture that makes WRF-GC today.
!
! Thanks to the Harvard Atmospheric Chemistry Modeling Group & the GEOS-Chem support team,
! without whose help, especially integrating the changes (or disruptions) WRF-GC 
! has "constributed" to the GEOS-Chem code base, WRF-GC would not be in this position today.
! Special thanks to Prof. Daniel Jacob, who provided me with an amazing opportunity to
! be in Harvard in Summer '18 to work closer with the GEOS-Chem community.
!
! I am deeply indebted to my family, who although hear nothing specific about my work,
! have no doubt had more than enough of my hair-pulling mumbling while I worked late
! weekend nights. Their unconditional support of my endeavors throughout my life have
! been vital to make me who I am today.
!
! Last but not least, thanks to my fellow friend and companion "Angela" Huang, whom
! I had met mid-April at the peak of my frustrations with the WRF-to-GCHP interface.
! Long walks and rants at night in the PKU campus have helped provide the emotional
! support necessary to push this project further.
!
! Thank you everyone.
!                                                             - Haipeng Lin, 2018.8.11
! 
!
