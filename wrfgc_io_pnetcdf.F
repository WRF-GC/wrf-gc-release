module wrfgc_io_pnetcdf
   use mpi
   use pnetcdf
   ! WRF
   use module_domain_type, only: domain
   use module_utility
   use module_date_time
   ! GC
   use GC_Stateful_Mod, only: Global_DiagList, Global_Input_Opt,Global_HcoConfig
   use PRECISION_MOD
   use Input_Opt_Mod, only: OptInput
   use State_Chm_Mod, only: ChmState
   use State_Met_Mod, only: MetState
   use State_Diag_Mod, only: DgnState, Get_TagInfo, Get_NameInfo
   use State_Grid_Mod, only: GrdState
   USE Grid_Registry_Mod,     ONLY : Lookup_Grid
   use Registry_Mod
   use registry_params_mod
   use Charpak_Mod
   use ERROR_MOD
   use MetaHistContainer_Mod, ONLY: MetaHistContainer
   USE HistContainer_Mod,     ONLY : HistContainer
   USE MetaHistItem_Mod
   USE HistItem_Mod,         ONLY : HistItem
   ! HEMCO
   USE HCO_TYPES_MOD,     ONLY : DiagnCollection, DiagnCont, DiagnBundle, ConfigObj
   use HCO_Interface_Common, only : GetHcoDiagn
   use HCO_State_GC_Mod, only : HcoState, ExtState

   ! WRF-GC
   use WRFGC_History_Mod, only: CollectionList, CollectionCount, IndexVariables


   ! use DiagList_Mod,      ONLY : CollList, ColItem, DgnList, DgnItem
   implicit none

   INTEGER, PUBLIC, PARAMETER :: GC_SUCCESS =  0   ! Routine returns success
   INTEGER, PUBLIC, PARAMETER :: GC_FAILURE = -1   ! Routine returns failure
   CHARACTER(LEN=9), PARAMETER, PUBLIC :: UNDEFINED_STR      = 'not found'


   public :: writeDiag
contains
   subroutine check(err, message)
      use mpi
      use pnetcdf
      implicit none
      integer :: err
      character(len=*) :: message

      ! It is a good idea to check returned value for possible error
      if (err .NE. NF90_NOERR) then
         write(*,*)'error: something wrong ', trim(message)
         write(6,*) trim(message), " ",trim(nf90mpi_strerror(err))
         call MPI_Abort(MPI_COMM_WORLD, -1, err)
      end if
   end subroutine check

   subroutine construct_filename(time_str, filenametemplate1, filename, domainName)
      character(len=*), intent(in) :: time_str, filenametemplate1, domainName
      character(len=*), intent(out) :: filename


      character :: c
      character(len=5) :: cs
      integer :: i, start, end, ifilename
      character(len=255) :: filenametemplate

      filename = 'not'

      filenametemplate = filenametemplate1
      start = 0
      end = 0
      ifilename = 1


      do i = 1, 5
         start = index(filenametemplate, "%")
         if (start > 0) then
            ! copy the string before the %
            filename(ifilename:ifilename+start-end-1) = filenametemplate(end+1:start-1)
            ifilename = ifilename + start - end -1
            ! copy the string after the %
            if (filenametemplate(start+2:start+2) == "2") then
               end = start + 2
               select case(filenametemplate(start+1:start+1))
                case('m')
                  filename(ifilename:ifilename+1) = time_str(6:7)
                case('d')
                  filename(ifilename: ifilename+1) = time_str(9:10)
                case('h')
                  filename(ifilename:ifilename+1) = time_str(12:13)
                case('n')
                  filename(ifilename:ifilename+1) = time_str(15:16)
               end select
               ifilename = ifilename + 2
            else
               filename(ifilename: ifilename+3) = time_str(1:4)
               end = start + 2
               ifilename = ifilename + 4
            end if
            filenametemplate(start:start) = '_'

         end if

      end do
      filename(ifilename:ifilename+len_trim(filenametemplate)-end) = filenametemplate(end+1:len_trim(filenametemplate))
      filename = trim(domainName)//filename
   endsubroutine construct_filename

   subroutine checkgc(RC, message)
      use mpi
      implicit none
      integer :: RC
      character(len=*) :: message

      ! It is a good idea to check returned value for possible error
      if (RC .NE. GC_SUCCESS) then
         call MPI_Abort(MPI_COMM_WORLD, -1, RC)
      end if
   end subroutine checkgc

   subroutine writeDiag(debug_level, Input_Opt, State_Met, State_Chm, grid, State_Grid, State_Diag, its, ite, jts, jte, ide, jde, kte)
      USE GC_Stateful_Mod,   ONLY : GIGC_States


      type(OptInput), intent(inout) :: Input_Opt
      type(MetState), intent(inout) :: State_Met
      type(ChmState), intent(inout) :: State_Chm
      type(GrdState), intent(inout) :: State_Grid
      type(DgnState), intent(in) :: State_Diag
      type(domain), target :: grid
      integer(kind=MPI_OFFSET_KIND):: i, time_index, time_count
      integer, intent(in)::debug_level
      INTEGER:: RC, N           ! Success or failure?
      integer :: time_unixs(1), dimensions(3)
      integer(kind=MPI_OFFSET_KIND), intent(in) :: ide, jde, kte, its, ite, jts, jte

      integer :: ncid, ierr, dimids(8), xtimei, varid, time_unix, nTags, year, month, day, hour, minute, second, time_interval
      character(len=50) :: filename, time_unix_start
      real :: xtime
      type(WRFU_Clock) :: clock
      TYPE(WRFU_Time) :: currTime, startTime, stopTime, referenceTimem, time
      TYPE(WRFU_TimeInterval) :: timeStep
      CHARACTER (LEN=64) :: currTime_str, startTime_str, stopTime_str, filename_time_str, time_str
      CHARACTER (LEN=19) :: referenceTime_str, new_referenceTime_str
      CHARACTER (LEN=64) :: diag_name, state_name
      integer, dimension(2) :: s
      type(MetaRegItem), pointer:: current_reg
      ! type(DgnItem), pointer :: current_diag
      character(len=255) :: ItemTemplate, ItemName, ItemTemplateUC, ItemNameUC, ItemPrefix,SubStrs(255), tagId, &
         tagName, ErrMsg, OutputName, Description, Units
      character(len=80) :: ErrorLine
      character(len=5) :: StateMetUC, StateChmUC, domainName
      logical :: Found
      type(HistContainer) , POINTER:: Current_container
      type(MetaHistContainer), POINTER :: Current_metahistcontainer
      type(HistItem), POINTER :: Current_Item
      type(MetaHistItem), POINTER :: Current_metahistitem




      REAL(f8),            POINTER :: Grid_Lat (:    )
      REAL(f8),            POINTER :: Grid_LatE(:    )
      REAL(f8),            POINTER :: Grid_Lon (:    )
      REAL(f8),            POINTER :: Grid_LonE(:    )

      REAL(f8), POINTER  :: Source_0d_8           ! Ptr to 0D 8-byte    data
      REAL(f4), pointer  :: Source_0d_4           ! Ptr to 0D 4-byte    data
      INTEGER,  POINTER  :: Source_0d_I           ! Ptr to 0D integer   data

      REAL(f8), POINTER  :: Source_1d_8(:    )    ! Ptr to 1D 8-byte    data
      REAL(f4), POINTER  :: Source_1d_4(:    )    ! Ptr to 1D 4-byte    data
      INTEGER,  POINTER  :: Source_1d_I(:    )    ! Ptr to 1D integer   data

      REAL(f8), POINTER  :: Source_2d_8(:,:  )    ! Ptr to 2D 8-byte    data
      REAL(f4), POINTER  :: Source_2d_4(:,:  )    ! Ptr to 2D 4-byte    data
      INTEGER,  POINTER  :: Source_2d_I(:,:  )    ! Ptr to 2D integer   data

      REAL(f8), POINTER  :: Source_3d_8(:,:,:)    ! Ptr to 3D 8-byte    data
      REAL(f4), POINTER  :: Source_3d_4(:,:,:)    ! Ptr to 3D 4-byte    data
      INTEGER,  POINTER  :: Source_3d_I(:,:,:)    ! Ptr to 3D integer   data

      REAL(f8), POINTER  :: Data_0d               ! 0D scalar
      REAL(f8), POINTER  :: Data_1d(:    )        ! 1D vector
      REAL(f8), POINTER  :: Data_2d(:,:  )        ! 2D array
      REAL(f8), POINTER  :: Data_3d(:,:,:)        ! 3D array

      REAL :: T1, T2 ! for measuring exctuion time

      

      ! for hemco
      type(DiagnCollection), pointer :: hco_collection
      type(DiagnCont), pointer :: hco_diag
      logical :: do_diagn_now 
      integer :: delta_s

      Source_0d_8 => NULL()
      Source_1d_8 => NULL()
      Source_1d_4 => NULL()
      Source_1d_I => NULL()
      Source_2d_8 => NULL()
      Source_2d_4 => NULL()
      Source_2d_I => NULL()
      Source_3d_8 => NULL()
      Source_3d_4 => NULL()
      Source_3d_I => NULL()

      Data_0d => NULL()
      Data_1d => NULL()
      Data_2d => NULL()
      Data_3d => NULL()

      RC = GC_SUCCESS


      ! do nothing if function called from chemics_init.F
      if (.not. associated(GIGC_States(grid%id)%HcoState))then
         return
      endif

      if (debug_level .gt. 200)then
         write(*,*) 'The current output domain: ',State_Grid%ID, grid%id;
      endif
      write(domainName, '(A,I2.2)'), 'd', grid%id;
      ! get time
      clock = grid%domain_clock
      CALL WRFU_ClockGet( clock, CurrTime=currTime, StartTime=startTime, &
         StopTime=stopTime, TimeStep=timeStep, rc=rc )

      call wrf_timetoa( startTime, startTime_str )
      call wrf_timetoa( currTime, currTime_str )
      call wrf_timetoa( stopTime, stopTime_str )

      xtime = grid%xtime
      ! get filename
      ! use wrf history interval as output interval
      xtimei = int(xtime/grid%history_interval)

      write(*,*)'before get collection, line 241'
      ! write hemco diagnostic
      hco_collection => GIGC_States(State_Grid%id)%HcoState%Diagn%Collections

      write(*,*)'after get collection, line 245'

        do while(associated(hco_collection))
         do_diagn_now = .False.
         if (debug_level .ge. 500) then
            write(*,*)'write hemco diagn: ', trim(hco_collection%PREFIX)
            write(*,*)'frequency: ',hco_collection%deltaYMD, hco_collection%deltaHMS
            
         endif
         if ((hco_collection%deltaHMS == 1 .and. hco_collection%deltaYMD == 0) .or. hco_collection%deltaYMD == -1)then
            ! always output
            do_diagn_now = .True.
            write(filename, *)trim(hco_collection%PREFIX),"_",trim(domainName), "_", currTime_str(1:4),'_',currTime_str(6:7),"_", currTime_str(9:10),currTime_str(12:13),currTime_str(15:16), currTime_str(18:19), '.nc'
         else if (hco_collection%lastYMD == -1)then 
            ! first time output, using start time and update lastymd
            do_diagn_now = .True.
            read(startTime_str(1:4), *) year
            read(startTime_str(6:7), *) month
            read(startTime_str(9:10), *) day
            read(startTime_str(12:13), *) hour
            read(startTime_str(15:16), *) minute
            read(startTime_str(18:19), *) second
            write(filename, *)trim(hco_collection%PREFIX),"_",trim(domainName), "_", startTime_str(1:4),'_',currTime_str(6:7),"_", currTime_str(9:10),currTime_str(12:13),currTime_str(15:16), currTime_str(18:19), '.nc'

            hco_collection%lastYMD = year*10000 + month*100 + day
            hco_collection%lastHms = hour*10000 + minute*100 + second
         else
            ! if interval > time interval, output and update lastymd
            year = hco_collection%lastYMD / 10000
            month = (hco_collection%lastYMD - year*10000) / 100
            day = hco_collection%lastYMD - year*10000 - month*100
            hour = hco_collection%lastHms / 10000
            minute = (hco_collection%lastHms - hour*10000) / 100
            second = hco_collection%lastHms - hour*10000 - minute*100
            write(referenceTime_str,'(I4.4, "-", I2.2, "-", I2.2, " ", I2.2, ":", I2.2, ":", I2.2)') year, month, day, hour, minute, second
            call geth_idts(currTime_str(1:19), referenceTime_str(1:19), time_interval)
            year = hco_collection%deltaYMD / 10000
            month = (hco_collection%deltaYMD - year*10000) /100
            day = hco_collection%deltaYMD - year*10000 - month*100
            hour = hco_collection%deltaHms / 10000
            minute = (hco_collection%deltaHms - hour*10000) / 100
            second = hco_collection%deltaHms - hour*10000 - minute*100
            delta_s = year*365*24*60*60 + month*30*24*60*60 + day*24*60*60 + hour*60*60 + minute*60 + second
            if (time_interval > delta_s) then
               do_diagn_now = .True.
               call geth_newdate(new_referenceTime_str, referenceTime_str, delta_s)
               read(new_referenceTime_str(1:4), *) year
               read(new_referenceTime_str(6:7), *) month
               read(new_referenceTime_str(9:10), *) day
               read(new_referenceTime_str(12:13), *) hour
               read(new_referenceTime_str(15:16), *) minute
               read(new_referenceTime_str(18:19), *) second
               write(filename, *)trim(hco_collection%PREFIX),"_",trim(domainName), "_", new_referenceTime_str(1:4),'_',new_referenceTime_str(6:7),"_", new_referenceTime_str(9:10),new_referenceTime_str(12:13),new_referenceTime_str(15:16), new_referenceTime_str(18:19), '.nc'
               hco_collection%lastYMD = year*10000 + month*100 + day
               hco_collection%lastHms = hour*10000 + minute*100 + second
            else if (time_interval<11*60)then
               ! if interval < 11 minutes, output
               do_diagn_now = .True.
               write(filename, *)trim(hco_collection%PREFIX),"_",trim(domainName), "_", stopTime_str(1:4),'_',stopTime_str(6:7),"_", stopTime_str(9:10),stopTime_str(12:13),stopTime_str(15:16), stopTime_str(18:19), '.nc'
            ENDIF
         endif
         ! do_diagn_now = .False.
         ! ! last time = -1
         ! if (hco_collection%lastYMD == -1) do_diagn_now = .True.
         ! if (hco_collection%deltaYMD == -1) do_diagn_now = .True.
         ! call geth_idts(stopTime_str(1:19), currTime_str(1:19), time_interval)
         ! ! if reach end of simulation, do diagn
         ! if (time_interval<11*60) do_diagn_now = .True.
         ! ! if current time - last time > delta time
         ! if ( .not. do_diagn_now)then
         !    year = hco_collection%lastYMD / 10000
         !    month = (hco_collection%lastYMD - year*10000) / 100
         !    day = hco_collection%lastYMD - year*10000 - month*100
         !    hour = hco_collection%lastHms / 10000
         !    minute = (hco_collection%lastHms - hour*10000) / 100
         !    second = hco_collection%lastHms - hour*10000 - minute*100
         !    write(referenceTime_str,&
         !    '(I4.4, "-", I2.2, "-", I2.2, " ", I2.2, ":", I2.2, ":", I2.2)'&
         !    ) year, month, day, hour, minute, second
         !    call geth_idts(currTime_str(1:19), referenceTime_str(1:19), time_interval)
         !    ! convert deltaymd and deltahms to second
         !    year = hco_collection%deltaYMD / 10000
         !    month = (hco_collection%deltaYMD - year*10000) /100
         !    day = hco_collection%deltaYMD - year*10000 - month*100
         !    hour = hco_collection%deltaHms / 10000
         !    minute = (hco_collection%deltaHms - hour*10000) / 100
         !    second = hco_collection%deltaHms - hour*10000 - minute*100
         !    delta_s = year*365*24*60*60 + month*30*24*60*60 + day*24*60*60 + hour*60*60 + minute*60 + second
         !    if (time_interval > delta_s) do_diagn_now = .True.
         ! endif

         T1 = MPI_Wtime()
            if ( do_diagn_now )then
         !       !  write(*,* )'wrfgc_io_pnetcdf: write hemco diagn.'
         !        ! update lastymd and last hms
         !       if (hco_collection%deltaYMD == -1)then
         !          ! always output
         !          write(filename, *)trim(hco_collection%PREFIX), "_", currTime_str(1:4),'_',currTime_str(6:7),"_", currTime_str(9:10),currTime_str(12:13),currTime_str(15:16), currTime_str(18:19), '.nc'
         !       else if (hco_collection%lastYMD == -1)then
         !          read(startTime_str(1:4), *) year
         !          read(startTime_str(6:7), *) month
         !          read(startTime_str(9:10), *) day
         !          read(startTime_str(12:13), *) hour
         !          read(startTime_str(15:16), *) minute
         !          read(startTime_str(18:19), *) second
         !          write(filename, *)trim(hco_collection%PREFIX), "_", startTime_str(1:4),'_',currTime_str(6:7),"_", currTime_str(9:10),currTime_str(12:13),currTime_str(15:16), currTime_str(18:19), '.nc'
         !       else if (time_interval<11*60)then
         !          write(filename, *)trim(hco_collection%PREFIX), "_", stopTime_str(1:4),'_',stopTime_str(6:7),"_", stopTime_str(9:10),stopTime_str(12:13),stopTime_str(15:16), stopTime_str(18:19), '.nc'
         !       else
         !          ! convert deltaymd and deltahms to second
         !          year = hco_collection%deltaYMD / 10000
         !          month = (hco_collection%deltaYMD - year*10000) /100
         !          day = hco_collection%deltaYMD - year*10000 - month*100
         !          hour = hco_collection%deltaHms / 10000
         !          minute = (hco_collection%deltaHms - hour*10000) / 100
         !          second = hco_collection%deltaHms - hour*10000 - minute*100
         !          delta_s = year*365*24*60*60 + month*30*24*60*60 + day*24*60*60 + hour*60*60 + minute*60 + second
                  
         !          year = hco_collection%lastYMD / 10000
         !          month = (hco_collection%lastYMD - year*10000) / 100
         !          day = hco_collection%lastYMD - year*10000 - month*100
         !          hour = hco_collection%lastHms / 10000
         !          minute = (hco_collection%lastHms - hour*10000) / 100
         !          second = hco_collection%lastHms - hour*10000 - minute*100
         !          write(referenceTime_str,'(I4.4, "-", I2.2, "-", I2.2, " ", I2.2, ":", I2.2, ":", I2.2)') year, month, day, hour, minute, second
         !          call geth_newdate(new_referenceTime_str, referenceTime_str, delta_s)
         !          read(new_referenceTime_str(1:4), *) year
         !          read(new_referenceTime_str(6:7), *) month
         !          read(new_referenceTime_str(9:10), *) day
         !          read(new_referenceTime_str(12:13), *) hour
         !          read(new_referenceTime_str(15:16), *) minute
         !          read(new_referenceTime_str(18:19), *) second
         !          write(filename, *)trim(hco_collection%PREFIX), "_", new_referenceTime_str(1:4),'_',new_referenceTime_str(6:7),"_", new_referenceTime_str(9:10),new_referenceTime_str(12:13),new_referenceTime_str(15:16), new_referenceTime_str(18:19), '.nc'
         !       endif
         !       hco_collection%lastYMD = year*10000 + month*100 + day
         !       hco_collection%lastHms = hour*10000 + minute*100 + second

               filename = trim(ADJUSTL(filename))


               call create_file(filename, ncid, ide, jde, kte, dimids, State_Grid, Input_Opt)
               ierr = nf90mpi_redef(ncid)
               call check(ierr, 'redef')
               ! def var all at once
               hco_diag => hco_collection%DiagnList
               do while(associated(hco_diag))
                  if ( hco_diag%SpaceDim == 2 .and. associated(hco_diag%Arr2D))then
                     if (debug_level .ge. 500) then
                        write(*,*)'def hemco diagn 2d: name: ', hco_diag%cName
                     endif
                     if (associated(hco_diag%Arr2D%Val))then
                        ierr = nf90mpi_def_var(ncid, trim(hco_diag%cName), NF90_FLOAT, (/dimids(1), dimids(2), dimids(4)/), varid)
                        call check(ierr, 'def var')
                     else
                        write(*,*)'val not allocated'
                     endif
                  endif
                  if (hco_diag%SpaceDim == 3 .and. associated(hco_diag%Arr3D))then
                     if (debug_level .ge. 500) then
                        write(*,*)'def hemco diagn 3d : name: ', hco_diag%cName
                     endif
                     if (associated(hco_diag%Arr3D%Val))then
                        ierr = nf90mpi_def_var(ncid, trim(hco_diag%cName), NF90_FLOAT, (/dimids(1), dimids(2), dimids(3), dimids(4)/), varid)
                        call check(ierr, 'def var')
                     else 
                        write(*,*)'val not allocated'
                     endif
                  endif
                  ierr = nf90mpi_put_att(ncid, varid, 'units', trim(hco_diag%OutUnit))
                  ierr = nf90mpi_put_att(ncid, varid, 'long_name', trim(hco_diag%long_name))
                  hco_diag => hco_diag%NextCont
               enddo
               ierr = nf90mpi_enddef(ncid)
               call check(ierr, 'enddef')

                ! write dimension
                allocate(Grid_Lat(State_Grid%NY), STAT=RC)
                call checkgc(RC, 'allocate Grid_lat')
                do i = 1,State_Grid%NY
                Grid_Lat(i) = State_Grid%YMid(1, i)
                enddo
                call writeDimension('lat', Grid_Lat, ncid, jts, jte)

                allocate(Grid_Lon(State_Grid%NX), STAT=RC)
                call checkgc(RC, 'allocate Grid_lon')
                do i = 1,State_Grid%NX
                Grid_Lon(i) = State_Grid%XMid(i, 1)
                enddo
                call writeDimension('lon', Grid_Lon, ncid, its, ite)

                allocate(Grid_LatE(State_Grid%NY+1), STAT=RC)
                call checkgc(RC, 'allocate Grid_lat')
                do i = 1,State_Grid%NY+1
                Grid_LatE(i) = State_Grid%YEdge(1, i)
                enddo
                call writeDimension('lat_e', Grid_LatE, ncid, jts, jte+1)

                allocate(Grid_LonE(State_Grid%NX+1), STAT=RC)
                call checkgc(RC, 'allocate Grid_lon')
                do i = 1,State_Grid%NX+1
                Grid_LonE(i) = State_Grid%XEdge(i, 1)
                enddo
                call writeDimension('lon_e', Grid_LonE, ncid, its, ite+1)

                ! write time
                ierr = nf90mpi_inq_varid(ncid, 'time', varid)
                call check(ierr, 'inq varid time')
                time_index = 1
                time_count = 1
                time_unix_start = "1970-01-01_00:00:00"
                call geth_idts(currTime_str(1:19), time_unix_start(1:19), time_unix)
                time_unixs(1) = time_unix
                ierr = nf90mpi_put_var_all(ncid, varid, time_unixs, (/time_index/), (/time_count/))
                call check(ierr, 'put var time')

                ! write data array
                hco_diag => hco_collection%DiagnList
                do while(associated(hco_diag))
                  Source_2d_4 => NULL()
                  Source_3d_4 => NULL()
                  if(debug_level .ge. 500)then
                     write(*,* )'write hemco diagn: ', hco_diag%cID
                  endif
                    if ( hco_diag%SpaceDim == 2 .and. associated(hco_diag%Arr2D))then
                        CALL GetHcoDiagn( HcoState, ExtState,hco_diag%cName, .false., RC, Ptr2D=Source_2d_4, &
                                          COL=HcoState%Diagn%HcoDiagnIDDefault)
                        if (debug_level .ge. 500) then
                            write(*,*)'write hemco diagn 2d: name: ', hco_diag%cName
                        endif
                        if (associated(Source_2d_4))then
                            call write2file(trim(hco_diag%cName), Source_2d_4, &
                                            Data_2d, its, ite, jts, jte, ide, jde, &
                                            time_index, ncid, (/dimids(1:2), dimids(4)/),&
                                            trim(hco_diag%OutUnit), trim(hco_diag%long_name))
                        else
                            write(*,*)'val not allocated'
                        endif
                    endif
                    if (hco_diag%SpaceDim == 3 .and. associated(hco_diag%Arr3D))then
                        CALL GetHcoDiagn( HcoState, ExtState,hco_diag%cName, .false., RC, Ptr3D=Source_3d_4, &
                                          COL=HcoState%Diagn%HcoDiagnIDDefault)
                        if (debug_level .ge. 500) then
                            write(*,*)'write hemco diagn 3d : name: ', hco_diag%cName
                        endif
                        if (associated(Source_3d_4))then
                            call write2file3d(trim(hco_diag%cName), Source_3d_4, &
                                            Data_3d, its, ite, jts, jte, ide, jde, &
                                            kte, time_index, ncid, dimids, &
                                            trim(hco_diag%OutUnit), trim(hco_diag%long_name))
                        else 
                            write(*,*)'val not allocated'
                        endif
                    endif
                    Source_2d_4 => NULL()
                     Source_3d_4 => NULL()
                    hco_diag => hco_diag%NextCont
                enddo
                call closefile(ncid)
            endif
            T2 = MPI_Wtime()
            write(*,*)'hemco write time: ', T2-T1
            hco_collection => hco_collection%NextCollection
        enddo

      ! clean up 
      hco_collection => null()
      Grid_Lat => null()
      Grid_Lon => null()
      Grid_LatE => null()
      Grid_LonE => null()
      Source_2d_4 => NULL()
      Source_3d_4 => NULL()
    
      ! write wrfchem diagnostic

      Current_metahistcontainer => CollectionList(State_Grid%id)%p
      do while(associated(Current_metahistcontainer))
         ! open file according to the reference time
         Current_container => Current_metahistcontainer%Container
         year = Current_container%ReferenceYmd / 10000
         month = (Current_container%ReferenceYmd - year*10000) / 100
         day = Current_container%ReferenceYmd - year*10000 - month*100
         hour = Current_container%ReferenceHms / 10000
         minute = (Current_container%ReferenceHms - hour*10000) / 100
         second = Current_container%ReferenceHms - hour*10000 - minute*100
         write(referenceTime_str, '(I4.4, "-", I2.2, "-", I2.2, " ", I2.2, ":", I2.2, ":", I2.2)') year, month, day, hour, minute, second
         ! if referencetime < starttime, create new file
         ! referencetime  = starttime + n*fileclosetime

         call geth_idts(referenceTime_str(1:19), startTime_str(1:19), time_interval)
         ! write(*,*)'file close time_interval: ', time_interval
         ! write(*,*)'reference time: ', trim(referenceTime_str)
         ! write(*,*)'start time: ', trim(startTime_str)
         if (time_interval<0) then
            ! At the begining of simulation
            call construct_filename(startTime_str, Current_container%Filename, filename, domainName)
            call create_file(filename, Current_container%FileId, ide, jde, kte, dimids, State_Grid, Input_Opt)
            read(currTime_str(1:4), *) year
            read(currTime_str(6:7), *) month
            read(currTime_str(9:10), *) day
            read(currTime_str(12:13), *) hour
            read(currTime_str(15:16), *) minute
            read(currTime_str(18:19), *) second
            ! update current container
            Current_container%ReferenceYmd = year*10000 + month*100 + day
            Current_container%ReferenceHms = hour*10000 + minute*100 + second
            Current_container%xDimId = dimids(1)
            Current_container%yDimId = dimids(2)
            Current_container%zDimId = dimids(3)
            Current_container%tDimId = dimids(4)
            Current_container%IsFileOpen = .true.
            ! write(*,*)'after create file: ', Current_container%FileId
            ! write(*,*)'file dimids: ', Current_container%xDimId, Current_container%yDimId, Current_container%zDimId, Current_container%tDimId
            ! write(*,*)'file is open: ', Current_container%IsFileOpen
         else
            ! if currenttime > referencetime + fileclosetime, close file then create new file
            call geth_idts(currTime_str(1:19), referenceTime_str(1:19), time_interval)
            ! write(*,*)'time interval: ', time_interval
            if ( time_interval >= Current_container%FileCloseIvalSec)then
               ! write(*,*)'close file', Current_container%FileId
               call closefile(Current_container%FileId)
               call geth_newdate(new_referenceTime_str(1:19), referenceTime_str(1:19), int(Current_container%FileCloseIvalSec))
               call construct_filename(new_referenceTime_str, Current_container%Filename, filename, domainName)
               filename = trim(filename)
               call create_file(filename, Current_container%FileId, ide, jde, kte, dimids, State_Grid, Input_Opt)
               read(new_referenceTime_str(1:4), *) year
               read(new_referenceTime_str(6:7), *) month
               read(new_referenceTime_str(9:10), *) day
               read(new_referenceTime_str(12:13), *) hour
               read(new_referenceTime_str(15:16), *) minute
               read(new_referenceTime_str(18:19), *) second
               ! update current container
               Current_container%ReferenceYmd = year*10000 + month*100 + day
               Current_container%ReferenceHms = hour*10000 + minute*100 + second
               Current_container%xDimId = dimids(1)
               Current_container%yDimId = dimids(2)
               Current_container%zDimId = dimids(3)
               Current_container%tDimId = dimids(4)
               Current_container%IsFileOpen = .true.
            else
               ! write(*,*)'is file opened?',Current_container%IsFileOpen
               if (Current_container%IsFileOpen) then
                  ! write(*,*)'file is open: ', Current_container%FileId
                  ! if currenttime < referencetime + fileclosetime, use current file
                  dimids(1) = Current_container%xDimId
                  dimids(2) = Current_container%yDimId
                  dimids(3) = Current_container%zDimId
                  dimids(4) = Current_container%tDimId
               else
                  ! write(*,*)'file is not open'
                  ! currenttime < referencetime + fileclosetime, but file is not open, open file
                  call construct_filename(referenceTime_str, Current_container%Filename, filename, domainName)
                  filename = trim(filename)
                  call create_file(filename, Current_container%FileId, ide, jde, kte, dimids, State_Grid, Input_Opt)
                  Current_container%xDimId = dimids(1)
                  Current_container%yDimId = dimids(2)
                  Current_container%zDimId = dimids(3)
                  Current_container%tDimId = dimids(4)
                  Current_container%IsFileOpen = .true.
               endif
            endif
         endif

         ! write file according to currentcontainer%currentYmd
         !
         year = Current_container%CurrentYmd / 10000
         month = (Current_container%CurrentYmd - year*10000) / 100
         day = Current_container%CurrentYmd - year*10000 - month*100
         hour = Current_container%CurrentHms / 10000
         minute = (Current_container%CurrentHms - hour*10000) / 100
         second = Current_container%CurrentHms - hour*10000 - minute*100
         write(time_str, '(I4.4, "-", I2.2, "-", I2.2, " ", I2.2, ":", I2.2, ":", I2.2)') year, month, day, hour, minute, second
         call geth_idts( currTime_str(1:19), time_str(1:19), time_interval)
         ! 1. interval < 0, write to current file, update container%currentYmd
         ! 2. interval > 0, interval < fileupdateinterval, skip
         ! 3. interval > 0, interval > fileupdateinterval, write to current file and update container%currentYmd
         ! write(*,*)'current diag: ', trim(Current_container%Name)
         ! write(*,*)'time_interval: ', time_interval
         ! write(*,*)'updateinterval: ', Current_container%UpdateIvalSec
         ! write(*,*)'current time:', trim(currTime_str)
         ! write(*,*)'last write time: ', trim(time_str)
         if (time_interval <= 0 .or. time_interval >= Current_container%UpdateIvalSec)then
            ! write dimensions to variables
            ierr = nf90mpi_inquire_dimension(Current_container%FileId, dimids(4), len=time_index)
            call check(ierr, 'inq dim time')
            ! write(*,*)'time length: ', time_index

            if (time_index == 0) then
               if (debug_level .ge. 500) then
                  write(*,*)'def vars'
               endif
               Current_metahistitem => Current_container%HistItems
               ierr = nf90mpi_redef(Current_container%FileId)
               call check(ierr, 'redef')
               ! def var all at once
               do while(associated(Current_metahistitem))
                  Current_Item => Current_metahistitem%Item
                  SELECT CASE(Current_Item%SpaceDim)
                     CASE( 3 )
                        ierr = nf90mpi_def_var(  Current_container%FileId, \
                                                trim(Current_Item%Name), \
                                                NF90_FLOAT, \
                                                (/dimids(1),dimids(2),dimids(3),dimids(4)/), \
                                                varid)
                        call check(ierr, 'def var')
                        Current_Item%NcVarId = varid
                     CASE( 2 )
                        ierr = nf90mpi_def_var(  Current_container%FileId, \
                                                trim(Current_Item%Name), \
                                                NF90_FLOAT, \
                                                (/dimids(1),dimids(2),dimids(4)/), \
                                                varid)
                        call check(ierr, 'def var')
                        Current_Item%NcVarId = varid
                  end select
                  ierr = nf90mpi_put_att(Current_container%FileId, \
                                          varid,\
                                          'units', \
                                          trim(Current_Item%Units))
                  ierr = nf90mpi_put_att(Current_container%FileId, \
                                          varid,\
                                          'long_name', \
                                          trim(Current_Item%LongName))
                  Current_metahistitem => Current_metahistitem%Next
               enddo
               ierr = nf90mpi_enddef(Current_container%FileId)
               call check(ierr, 'enddef, may caused by too many variables in the collection')

               if (debug_level .ge. 500) then
                  write(*,*)'write dimensions'
               endif

               allocate(Grid_Lat(State_Grid%NY), STAT=RC)
               call checkgc(RC, 'allocate Grid_lat')
               do i = 1,State_Grid%NY
                  Grid_Lat(i) = State_Grid%YMid(1, i)
               enddo
               call writeDimension('lat', Grid_Lat, Current_container%FileId, jts, jte)

               allocate(Grid_Lon(State_Grid%NX), STAT=RC)
               call checkgc(RC, 'allocate Grid_lon')
               do i = 1,State_Grid%NX
                  Grid_Lon(i) = State_Grid%XMid(i, 1)
               enddo
               call writeDimension('lon', Grid_Lon, Current_container%FileId, its, ite)

               allocate(Grid_LatE(State_Grid%NY+1), STAT=RC)
               call checkgc(RC, 'allocate Grid_lat')
               do i = 1,State_Grid%NY+1
                  Grid_LatE(i) = State_Grid%YEdge(1, i)
               enddo
               call writeDimension('lat_e', Grid_LatE, Current_container%FileId, jts, jte+1)

               allocate(Grid_LonE(State_Grid%NX+1), STAT=RC)
               call checkgc(RC, 'allocate Grid_lon')
               do i = 1,State_Grid%NX+1
                  Grid_LonE(i) = State_Grid%XEdge(i, 1)
               enddo
               call writeDimension('lon_e', Grid_LonE, Current_container%FileId, its, ite+1)

               ! write index variables only for domain 1
               if (State_Grid%ID == 1)then
                  write(*,*)'index variables: ', trim(IndexVariables(State_Grid%ID)%p%Item%Name)
                  Current_metahistitem => IndexVariables(State_Grid%ID)%p
               
                  do while(associated(Current_metahistitem))
                     Current_Item => Current_metahistitem%Item
                     SELECT CASE(Current_Item%SpaceDim)
                     CASE( 3 )
                        if (Current_Item%Source_KindVal == KINDVAL_F8) then
                           Current_Item%Data_3d = Current_Item%Source_3d_8
                        else if (Current_Item%Source_KindVal == KINDVAL_F4) then
                           Current_Item%Data_3d = Current_Item%Source_3d_4
                        else if (Current_Item%Source_KindVal == KINDVAL_I4) then
                           Current_Item%Data_3d = Current_Item%Source_3d_I
                        endif
                     CASE( 2 )
                        if (Current_Item%Source_KindVal == KINDVAL_F8) then
                           Current_Item%Data_2d = Current_Item%Source_2d_8
                        else if (Current_Item%Source_KindVal == KINDVAL_F4) then
                           Current_Item%Data_2d = Current_Item%Source_2d_4
                        else if (Current_Item%Source_KindVal == KINDVAL_I4) then
                           Current_Item%Data_2d = Current_Item%Source_2d_I
                        endif
                     CASE( 1 )
                        if (Current_Item%Source_KindVal == KINDVAL_F8) then
                           Current_Item%Data_1d = Current_Item%Source_1d_8
                        else if (Current_Item%Source_KindVal == KINDVAL_F4) then
                           Current_Item%Data_1d = Current_Item%Source_1d_4
                        else if (Current_Item%Source_KindVal == KINDVAL_I4) then
                           Current_Item%Data_1d = Current_Item%Source_1d_I
                        endif
                     CASE( 0 )
                        
                        Current_Item%Data_0d = Current_Item%Source_0d_8
                        
                     end select

                     SELECT CASE(Current_Item%Name)
                     ! CASE('lat_bnds')
                     !    ierr = nf90mpi_redef(Current_container%FileId)
                     !    call check(ierr, 'redef')
                     !    ierr = nf90mpi_def_var(  Current_container%FileId, \
                     !                            'lat_bnds', \
                     !                            NF90_FLOAT, \
                     !                            (/dimids(8),dimids(2)/), \
                     !                            varid)
                     !    call check(ierr, 'def var lat_bnds')
                     !    ierr = nf90mpi_enddef(Current_container%FileId)
                     !    call check(ierr, 'enddef')
                        
                     !    ierr = nf90mpi_put_var_all(   Current_container%FileId,\
                     !                                  varid,\
                     !                                  Current_Item%Data_2d,\
                     !                                  (/int(1,kind=MPI_OFFSET_KIND),jts/),\
                     !                                  (/int(2,kind=MPI_OFFSET_KIND),jte-jts+1/) )
                     !    call check(ierr, 'put var latbnds')
                     ! CASE('lon_bnds')
                     !    ierr = nf90mpi_redef(Current_container%FileId)
                     !    call check(ierr, 'redef')
                     !    ierr = nf90mpi_def_var(  Current_container%FileId, \
                     !                            'lon_bnds', \
                     !                            NF90_FLOAT, \
                     !                            (/dimids(8),dimids(1)/), \
                     !                            varid)
                     !    call check(ierr, 'def var lon_bnds')
                     !    ierr = nf90mpi_enddef(Current_container%FileId)
                     !    call check(ierr, 'enddef')

                     !    ierr = nf90mpi_put_var_all(   Current_container%FileId,\
                     !                                  varid,\
                     !                                  Current_Item%Data_2d,\
                     !                                  (/int(1,kind=MPI_OFFSET_KIND),its/),\
                     !                                  (/int(2,kind=MPI_OFFSET_KIND),ite-its+1/) )
                     !    call check(ierr, 'put var lonbnds')
                     CASE('lev')
                        ierr = nf90mpi_inq_varid(Current_container%FileId, 'lev', varid)
                        call check(ierr, 'inq lev')

                        ierr = nf90mpi_put_var_all(   Current_container%FileId,\
                                                      varid,\
                                                      Current_Item%Data_1d,\
                                                      (/int(1, MPI_OFFSET_KIND)/),\
                                                      (/kte/) )
                        call check(ierr, 'put var lev')
                     CASE('ilev')
                        ierr = nf90mpi_inq_varid(Current_container%FileId, 'ilev', varid)
                        call check(ierr, 'inq ilev')

                        ierr = nf90mpi_put_var_all(   Current_container%FileId,\
                                                      varid,\
                                                      Current_Item%Data_1d,\
                                                      (/int(1, MPI_OFFSET_KIND)/),\
                                                      (/kte+1/) )
                        call check(ierr, 'put var ilev')
                     CASE('hyam')
                        ierr = nf90mpi_redef(Current_container%FileId)
                        call check(ierr, 'redef')
                        ierr = nf90mpi_def_var(  Current_container%FileId, \
                                                'hyam', \
                                                NF90_FLOAT, \
                                                dimids(3), \
                                                varid)
                        call check(ierr, 'def var hyam')
                        ierr = nf90mpi_put_att(   Current_container%FileId,\
                                                   varid,\
                                                   'units',\
                                                   Current_Item%Units)
                        ierr = nf90mpi_put_att(   Current_container%FileId,\
                                                   varid,\
                                                   'long_name',\
                                                   Current_Item%LongName)  
                        ierr = nf90mpi_enddef(Current_container%FileId)
                        call check(ierr, 'enddef')

                        ierr = nf90mpi_put_var_all(   Current_container%FileId,\
                                                      varid,\
                                                      Current_Item%Data_1d,\
                                                      (/int(1, MPI_OFFSET_KIND)/),\
                                                      (/kte/) )
                        call check(ierr, 'put var hyam')
                     CASE('hybm')
                        ierr = nf90mpi_redef(Current_container%FileId)
                        call check(ierr, 'redef')
                        ierr = nf90mpi_def_var(  Current_container%FileId, \
                                                'hybm', \
                                                NF90_FLOAT, \
                                                dimids(3), \
                                                varid)
                        call check(ierr, 'def var hybm')
                        ierr = nf90mpi_put_att(   Current_container%FileId,\
                                                   varid,\
                                                   'units',\
                                                   Current_Item%Units)
                        ierr = nf90mpi_put_att(   Current_container%FileId,\
                                                   varid,\
                                                   'long_name',\
                                                   Current_Item%LongName)  
                        ierr = nf90mpi_enddef(Current_container%FileId)
                        call check(ierr, 'enddef')

                        ierr = nf90mpi_put_var_all(   Current_container%FileId,\
                                                      varid,\
                                                      Current_Item%Data_1d,\
                                                      (/int(1, MPI_OFFSET_KIND)/),\
                                                      (/kte/) )
                        call check(ierr, 'put var hybm')
                     CASE('hyai')
                        ierr = nf90mpi_redef(Current_container%FileId)
                        call check(ierr, 'redef')
                        ierr = nf90mpi_def_var(  Current_container%FileId, \
                                                'hyai', \
                                                NF90_FLOAT, \
                                                dimids(7), \
                                                varid)
                        call check(ierr, 'def var hyai')
                        ierr = nf90mpi_put_att(   Current_container%FileId,\
                                                   varid,\
                                                   'units',\
                                                   Current_Item%Units)
                        ierr = nf90mpi_put_att(   Current_container%FileId,\
                                                   varid,\
                                                   'long_name',\
                                                   Current_Item%LongName)  
                        ierr = nf90mpi_enddef(Current_container%FileId)
                        call check(ierr, 'enddef')

                        ierr = nf90mpi_put_var_all(   Current_container%FileId,\
                                                      varid,\
                                                      Current_Item%Data_1d,\
                                                      (/int(1, MPI_OFFSET_KIND)/),\
                                                      (/kte/) )
                        call check(ierr, 'put var hyai')
                     CASE('hybi')
                        ierr = nf90mpi_redef(Current_container%FileId)
                        call check(ierr, 'redef')
                        ierr = nf90mpi_def_var(  Current_container%FileId, \
                                                'hybi', \
                                                NF90_FLOAT, \
                                                dimids(7), \
                                                varid)
                        call check(ierr, 'def var hybi')
                        ierr = nf90mpi_put_att(   Current_container%FileId,\
                                                   varid,\
                                                   'units',\
                                                   Current_Item%Units)
                        ierr = nf90mpi_put_att(   Current_container%FileId,\
                                                   varid,\
                                                   'long_name',\
                                                   Current_Item%LongName)  
                        ierr = nf90mpi_enddef(Current_container%FileId)
                        call check(ierr, 'enddef')

                        ierr = nf90mpi_put_var_all(   Current_container%FileId,\
                                                      varid,\
                                                      Current_Item%Data_1d,\
                                                      (/int(1, MPI_OFFSET_KIND)/),\
                                                      (/kte/) )
                        call check(ierr, 'put var hybi')
                     ! ======= AREA is not working for now =======
                     ! CASE('AREA')
                     !    ierr = nf90mpi_redef(Current_container%FileId)
                     !    call check(ierr, 'redef')
                     !    ierr = nf90mpi_def_var(  Current_container%FileId, \
                     !                            'AREA', \
                     !                            NF90_FLOAT, \
                     !                            (/dimids(1), dimids(2)/), \
                     !                            varid)
                     !    call check(ierr, 'def var AREA')
                     !    ierr = nf90mpi_put_att(   Current_container%FileId,\
                     !                               varid,\
                     !                               'units',\
                     !                               Current_Item%Units)
                     !    ierr = nf90mpi_put_att(   Current_container%FileId,\
                     !                               varid,\
                     !                               'long_name',\
                     !                               Current_Item%LongName)  
                     !    ierr = nf90mpi_enddef(Current_container%FileId)
                     !    call check(ierr, 'enddef')

                     !    ierr = nf90mpi_put_var_all(   Current_container%FileId,\
                     !                                  varid,\
                     !                                  Current_Item%Data_2d,\
                     !                                  (/its, jts/),\
                     !                                  (/ite-its+1, jte-jts+1/) )
                     !    call check(ierr, 'put var AREA')
                     CASE('P0')
                        ierr = nf90mpi_redef(Current_container%FileId)
                        call check(ierr, 'redef')
                        ierr = nf90mpi_def_var(  Current_container%FileId, \
                                                'P0', \
                                                NF90_FLOAT, \
                                                varid)
                        call check(ierr, 'def var P0')
                        
                        ierr = nf90mpi_put_att(   Current_container%FileId,\
                                                   varid,\
                                                   'units',\
                                                   Current_Item%Units)
                        ierr = nf90mpi_put_att(   Current_container%FileId,\
                                                   varid,\
                                                   'long_name',\
                                                   Current_Item%LongName)                                                
                        call check(ierr, 'put att P0')

                        ierr = nf90mpi_enddef(Current_container%FileId)
                        call check(ierr, 'enddef')

                        ierr = nf90mpi_put_var_all(   Current_container%FileId,\
                                                      varid,\
                                                      Current_Item%Data_0d)
                        call check(ierr, 'put var P0')

                     END SELECT
                     Current_metahistitem => Current_metahistitem%Next
                  end do
               endif

            endif

            ! write time
            ierr = nf90mpi_inq_varid(Current_container%FileId, 'time', varid)
            call check(ierr, 'inq varid time')
            time_index = time_index + 1
            time_count = 1
            time_unix_start = "1970-01-01_00:00:00"
            call geth_idts(currTime_str(1:19), time_unix_start(1:19), time_unix)
            time_unixs(1) = time_unix
            ierr = nf90mpi_put_var_all(Current_container%FileId, varid, time_unixs, (/time_index/), (/time_count/))
            call check(ierr, 'put var time')

            ! write the data
            ! if instantaneous, write the source data
            dimensions = (/Current_container%nX, Current_container%nY, Current_container%nZ/)
            if (To_UpperCase(trim(Current_container%UpdateMode)) == To_UpperCase('instantaneous'))then
               Current_metahistitem => Current_container%HistItems
               do while (associated(Current_metahistitem))
                  ! write(*,*) 'in while loop: ', trim(Current_metahistitem%Item%Name)
                  Current_Item => Current_metahistitem%Item
                  ! write(*,*) 'current item ', Current_Item%Name
                  if (wrf_debug_level .ge. 500) then
                     write(*,*)'write item: ', Current_Item%Name
                     write(*,*)'write item: ', Current_Item%Units
                     write(*,*)'write item: ', Current_Item%LongName
                  endif
                  write(*,*)'write item: ', Current_Item%LongName
                  call writeItemName( Current_Item%Name, Input_Opt, State_Diag, State_Chm, State_Met,&
                     its, ite, jts, jte, ide, jde, kte, time_index, Current_container%FileId, dimids(1:4), &
                     Current_Item%SpaceDim, dimensions, &
                     Current_Item%Source_0d_8, Current_Item%Source_1d_8, Current_Item%Source_2d_8, Current_Item%Source_3d_8, &
                     Source_0d_4, Current_Item%Source_1d_4, Current_Item%Source_2d_4, Current_Item%Source_3d_4, &
                     Source_0d_I, Current_Item%Source_1d_I, Current_Item%Source_2d_I, Current_Item%Source_3d_I, &
                     Current_Item%Units, Current_Item%LongName)
                  Current_metahistitem => Current_metahistitem%Next
               enddo
               ! if timeaveraged, average the data_2d/data_3d, and write
            else if (To_UpperCase(trim(Current_container%UpdateMode)) == To_UpperCase('TIMEAVERAGED'))then
               Current_metahistitem => Current_container%HistItems
               do while(associated(Current_metahistitem))
                  Current_Item => Current_metahistitem%Item
                  if (Current_Item%SpaceDim == 2) Then
                     if (.not. associated(Data_2d))then
                        allocate(Data_2d(Current_container%nX, Current_container%nY))
                     endif
                     ! nUpdates was used to store the # of time averaged data
                     if ( associated(Current_Item%Source_2d_4))then
                        Data_2d = (Current_Item%Data_2d + Current_Item%Source_2d_4)/(Current_Item%nUpdates + 1)
                     else if ( associated(Current_Item%Source_2d_8))then
                        Data_2d = (Current_Item%Data_2d + Current_Item%Source_2d_8)/(Current_Item%nUpdates + 1)
                     endif
                     Current_Item%Data_2d = 0.0_f8

                  else if (Current_Item%SpaceDim == 3) Then
                     if (.not. associated(Data_3d))then
                        allocate(Data_3d(Current_container%nX, Current_container%nY, Current_container%nZ))
                     endif
                     ! write(*,*)'in average: 3d_4', Current_Item%Source_3d_4
                     ! write(*,*)Current_Item%Data_3d + Current_Item%Source_3d_4
                     ! write(*,*)Current_Item%nUpdates+1
                     ! write(*,*)(Current_Item%Data_3d + Current_Item%Source_3d_4)/(Current_Item%nUpdates + 1)
                     ! nUpdates was used to store the # of time averaged data
                     if (associated(Current_Item%Source_3d_8))then
                        Data_3d = (Current_Item%Data_3d + Current_Item%Source_3d_8)/(Current_Item%nUpdates + 1)
                     else if (associated(Current_Item%Source_3d_4))then
                        Data_3d = (Current_Item%Data_3d + Current_Item%Source_3d_4)/(Current_Item%nUpdates + 1)
                     endif
                     Current_Item%Data_3d = 0.0_f8
                  endif
                  call writeItemName( Current_Item%Name, Input_Opt, State_Diag, State_Chm, State_Met, &
                     its, ite, jts, jte, ide, jde, kte, time_index, Current_container%FileId, dimids(1:4), &
                     Current_Item%SpaceDim, dimensions, &
                     Data_0d, Data_1d, Data_2d, Data_3d, &
                     Source_0d_4, Source_1d_4, Source_2d_4, Source_3d_4, &
                     Source_0d_I, Source_1d_I, Source_2d_I, Source_3d_I, &
                     Current_Item%Units, Current_Item%LongName)
                  Data_2d => NULL()
                  Data_3d => NULL()
                  Current_Item%nUpdates = 0
                  Current_metahistitem => Current_metahistitem%Next
               end do

            endif
            ! update currentcontainer%currentYmd
            read(currTime_str(1:4), *) year
            read(currTime_str(6:7), *) month
            read(currTime_str(9:10), *) day
            read(currTime_str(12:13), *) hour
            read(currTime_str(15:16), *) minute
            read(currTime_str(18:19), *) second
            ! write(*,*)'update currentymd: ', year, month, day, hour, minute, second
            ! write(*,*)'currentymd before: ', Current_container%CurrentYmd
            Current_container%CurrentYmd = year*10000 + month*100 + day
            Current_container%CurrentHms = hour*10000 + minute*100 + second
            ! write(*,*)'currentymd after: ', Current_container%CurrentYmd
         else
            ! skip output, update data_2d/data_3d and nUpdates
            Current_metahistitem => Current_container%HistItems
            do while(associated(Current_metahistitem))
               Current_Item => Current_metahistitem%Item
               if (Current_Item%SpaceDim == 2) Then
                  if ( associated(Current_Item%Source_2d_4))then
                     Current_Item%Data_2d = Current_Item%Data_2d + Current_Item%Source_2d_4
                  else if ( associated(Current_Item%Source_2d_8))then
                     Current_Item%Data_2d = Current_Item%Data_2d + Current_Item%Source_2d_8
                  endif
               else if (Current_Item%SpaceDim == 3) Then
                  if ( associated(Current_Item%Source_3d_4))then
                     Current_Item%Data_3d = Current_Item%Data_3d + Current_Item%Source_3d_4
                  else if ( associated(Current_Item%Source_3d_8))then
                     Current_Item%Data_3d = Current_Item%Data_3d + Current_Item%Source_3d_8
                  endif
               endif
               if ( Current_Item%nUpdates < 0)then
                  Current_Item%nUpdates = 0
               endif
               Current_Item%nUpdates = Current_Item%nUpdates + 1

               Current_metahistitem => Current_metahistitem%Next
            enddo
         endif
         ! todo: close file if the currenttime + 10*60 > endtime
         call geth_idts(stopTime_str(1:19), currTime_str(1:19),  time_interval)
         if (time_interval < 10*60) then
            ! write(*,*)'close file: ', Current_container%FileId
            ! write(*,*)'time interval: ', time_interval
            call closefile(Current_container%FileId)
            Current_container%IsFileOpen = .false.
         endif
         Current_metahistcontainer => Current_metahistcontainer%Next
      enddo

   end subroutine writeDiag


   subroutine writeItemName(ItemName, Input_Opt, State_Diag, State_Chm, &
      State_Met, its, ite, jts, jte, ide, jde, kte, &
      time_index, ncid, dimids, rank, dimensions, &
      Ptr0d_8, Ptr1d_8, Ptr2d_8, Ptr3d_8, &
      Ptr0d_4, Ptr1d_4, Ptr2d_4, Ptr3d_4, &
      Ptr0d_I, Ptr1d_I, Ptr2d_I, Ptr3d_I,&
      Units, long_name)
      USE Charpak_Mod,           ONLY : To_UpperCase

      integer(kind=MPI_OFFSET_KIND), intent(in) :: ide, jde, kte, its, ite, jts, jte, time_index
      integer, intent(in) :: ncid, dimids(4), rank, Dimensions(3)
      character(len=255) :: ItemName, StateChmUC, StateMetUC,Description, ItemNameUC
      integer :: RC
      type(OptInput), intent(inout) :: Input_Opt
      type(MetState), intent(inout) :: State_Met
      type(ChmState), intent(inout) :: State_Chm
      type(DgnState), intent(in) :: State_Diag
      ! Pointer arrays
      REAL(f8),            POINTER , intent(inout):: Ptr0d_8
      REAL(f4),            POINTER , intent(inout):: Ptr0d_4
      INTEGER,             POINTER , intent(inout):: Ptr0d_I
      REAL(f8),            POINTER , intent(inout):: Ptr1d_8(:    )
      REAL(f4),            POINTER , intent(inout):: Ptr1d_4(:    )
      INTEGER,             POINTER , intent(inout):: Ptr1d_I(:    )
      REAL(f8),            POINTER , intent(inout):: Ptr2d_8(:,:  )
      REAL(f4),            POINTER , intent(inout):: Ptr2d_4(:,:  )
      INTEGER,             POINTER , intent(inout):: Ptr2d_I(:,:  )
      REAL(f8),            POINTER , intent(inout):: Ptr3d_8(:,:,:)
      REAL(f4),            POINTER , intent(inout):: Ptr3d_4(:,:,:)
      INTEGER,             POINTER , intent(inout):: Ptr3d_I(:,:,:)

      LOGICAL                      :: OnLevelEdges
      INTEGER                      :: Source_KindVal
      INTEGER                      :: Output_KindVal
      INTEGER                      :: ItemDims(3)
      CHARACTER(LEN=3)             :: DimNames
      character(len=*), intent(in) :: Units, long_name

      ItemNameUC     = To_UpperCase( ItemName )
      StateMetUC     = State_Met%State // '_'   ! State_Met%State is uppercase
      StateChmUC     = State_Chm%State // '_'   ! State_Chm%State is uppercase

      RC             =  GC_SUCCESS
      Description    =  ''
      Source_KindVal =  0
      Output_KindVal =  0


      ! get the associated prt and output
      ! write(*,*)'get variable:',trim(ItemName)
      ! write(*,*)'dimenssion: ', Dimensions(1), dimensions(2), dimensions(3)
      ! write(*,*)'rank: ', rank
      select case (rank)
       case(2)
         if (Dimensions(1) == ite-its+1 .and. Dimensions(2) == jte-jts+1) then
            if (associated(Ptr2d_4) .or. associated(Ptr2d_8)) then
               call write2file(  ItemName, Ptr2d_4, Ptr2d_8, its, ite, jts, jte, ide, jde, &
                                 time_index, ncid, (/dimids(1:2), dimids(4)/),&
                                 Units, long_name)
            else if ( associated(Ptr2d_I) ) then
               write(*,*)'error: integer type not supported'
            end if
         else
            write(*,*)'error: dimension not match'
            write(*,*)'ide, jde: ', ite-its+1, jte-jts+1
            write(*,*)'Dimensions: ', Dimensions(1), Dimensions(2)
         endif
       case(3)
         if ( Dimensions(1) == ite-its+1 .and. Dimensions(2) == jte-jts+1  ) then
            if ( associated(Ptr3d_4) .or. associated(Ptr3d_8)) then
               call write2file3d(ItemName, Ptr3d_4, Ptr3d_8, its, ite, jts, jte, ide, jde, kte, &
                                 time_index, ncid, dimids(1:4), &
                                 Units, long_name)
            else if ( associated(Ptr3d_I) ) then
               write(*,*)'error: integer type not supported'
            end if
         else
            write(*,*)'error: dimension not match'
            write(*,*)'ide, jde, kte: ', ite-its+1, jte-jts+1, kte
            write(*,*)'Dimensions: ', Dimensions(1), Dimensions(2), Dimensions(3)
         end if
      end select

   end subroutine writeItemName

   subroutine writeDimension(dimname, data, ncid, ts, te)
      ! Writes data to a netcdf file
      character(len=*), intent(in) :: dimname
      integer(kind=MPI_OFFSET_KIND) :: starts(1), counts(1)
      integer(kind=MPI_OFFSET_KIND), intent(in) :: ts, te
      integer, intent(in) :: ncid
      REAL(KIND( REAL( 0.0, 8 ) )), intent(inout), dimension(:) :: data
      integer:: ierr, i, varid

      ierr = nf90mpi_inq_varid(ncid, dimname, varid)
      ! write data
      starts = (/ ts /)
      counts = (/ te-ts+1 /)
      ierr = nf90mpi_put_var_all(ncid, varid, data, starts, counts)
      call check(ierr, 'put var domain')
   end subroutine writeDimension

   ! Create a new file and define dimensions
   subroutine create_file(filename, ncid, ide, jde, kte, dimids, State_Grid, Input_Opt)
      use mpi
      use pnetcdf
      implicit none
      character(len=*), intent(in) :: filename
      integer(kind=MPI_OFFSET_KIND), intent(in) :: ide, jde, kte
      type(GrdState), intent(in) :: State_Grid
      type(OptInput), intent(in) :: Input_Opt
      integer, intent(out) :: ncid, dimids(8)
      integer :: ierr, varid
      logical :: exists
      REAL(f8),            POINTER :: Grid_Lat (:    )
      REAL(f8),            POINTER :: Grid_LatE(:    )
      REAL(f8),            POINTER :: Grid_Lon (:    )
      REAL(f8),            POINTER :: Grid_LonE(:    )
      integer :: RC, I, cmode


      Grid_Lat       => NULL()
      Grid_LatE      => NULL()
      Grid_Lon       => NULL()
      Grid_LonE      => NULL()

      RC = GC_SUCCESS

      ! initialize file&define dimensions
      cmode = IOR(NF90_CLOBBER, NF90_64BIT_DATA)
      ierr = nf90mpi_create(MPI_COMM_WORLD, trim(filename), cmode, &
         MPI_INFO_NULL, ncid)
      call check(ierr, 'create file '//filename)
      ierr = nf90mpi_def_dim(ncid, 'lon', ide-1, dimids(1))
      call check(ierr, 'def dim lon')
      ierr = nf90mpi_def_dim(ncid, 'lat', jde-1, dimids(2))
      call check(ierr, 'def dim lat')
      ierr = nf90mpi_def_dim(ncid, 'lev', kte, dimids(3))
      call check(ierr, 'def dim lev')
      ierr = nf90mpi_def_dim(ncid, 'time', 0, dimids(4))
      call check(ierr, 'def dim lev')
      ierr = nf90mpi_def_dim(ncid, 'lon_e', ide, dimids(5))
      call check(ierr, 'def dim lon_e')
      ierr = nf90mpi_def_dim(ncid, 'lat_e', jde, dimids(6))
      call check(ierr, 'def dim lat_e')
      ierr = nf90mpi_def_dim(ncid, 'ilev', kte+1, dimids(7))
      call check(ierr, 'def dim ilev')
      ierr = nf90mpi_def_dim(ncid, 'nb', 2, dimids(8))

      ierr = nf90mpi_def_var(ncid, 'lon', NF90_FLOAT, dimids(1), varid)
      call check(ierr, 'def var lon')
      ierr = nf90mpi_def_var(ncid, 'lat', NF90_FLOAT, dimids(2), varid)
      call check(ierr, 'def var lat')
      ierr = nf90mpi_def_var(ncid, 'lev', NF90_FLOAT, dimids(3), varid)
      call check(ierr, 'def var lev')
      ierr = nf90mpi_def_var(ncid, 'time', NF90_INT, dimids(4), varid)
      call check(ierr, 'def var time')

      ierr = nf90mpi_put_att(ncid, varid, 'units', 's')
      ierr = nf90mpi_put_att(ncid, varid, 'description', 'Unix time stamp')

      ierr = nf90mpi_def_var(ncid, 'lon_e', NF90_FLOAT, dimids(5), varid=varid)
      call check(ierr, 'def var lon_e')
      ierr = nf90mpi_def_var(ncid, 'lat_e', NF90_FLOAT, dimids(6), varid=varid)
      call check(ierr, 'def var lat_e')
      ierr = nf90mpi_def_var(ncid, 'ilev', NF90_FLOAT, dimids(7), varid=varid)
      call check(ierr, 'def var ilev')

      ierr = nf90mpi_enddef(ncid)
      call check(ierr, 'enddef create')


   end subroutine create_file

   subroutine closefile(ncid)
      use mpi
      use pnetcdf
      implicit none
      integer, intent(in) :: ncid
      integer :: ierr

      ierr = nf90mpi_close(ncid)
      call check(ierr, 'close file')
   end subroutine closefile

   ! the index in wrf is its, ite, jts, jte, kts, kte
   subroutine write2file3d(varname, data_4, data_8, its, ite, jts, jte, ide, jde, kte,time_s, ncid, dimids, units, longname)
      ! Writes data to a netcdf file
      character(len=*), intent(in) :: varname
      integer(kind=MPI_OFFSET_KIND) :: starts(4), counts(4)
      integer(kind=MPI_OFFSET_KIND), intent(in) :: its, ite, jts, jte, ide, jde, kte, time_s
      integer(kind=MPI_OFFSET_KIND) :: kts
      integer, intent(in) :: ncid, dimids(4)
      integer :: varid, ierr, i
      character(len=26) :: filename
      character(len=*), intent(in) :: units, longname
      ! two dimensions data
      REAL(KIND( REAL( 0.0, 8 ) )), intent(inout), dimension(:,:,:), pointer :: data_8
      REAL(KIND( REAL( 0.0, 4 ) )), intent(inout), dimension(:,:,:), pointer :: data_4

      kts = 1

      ! skip coordinates
      ierr = nf90mpi_inq_varid(ncid, varname, varid)
      ! call check(ierr, 'inq varid')
      if (ierr .eq. NF90_ENOTVAR) then
         ! define variable
         ierr = nf90mpi_redef(ncid)
         call check(ierr, 'redef')
         if (associated(data_4))then
            ierr = nf90mpi_def_var(ncid, varname, NF90_FLOAT, dimids, varid)
         else if (associated(data_8))then
            ierr = nf90mpi_def_var(ncid, varname, NF90_DOUBLE, dimids, varid)
         else
            write(*,*)"In write2file3d, no data attached to ", varname
         endif
         call check(ierr, 'def var')
         ierr = nf90mpi_put_att( ncid,\
                                 varid,\
                                 'units',\
                                 Units)
         ierr = nf90mpi_put_att( ncid,\
                                 varid,\
                                 'long_name',\
                                 LongName)  
         ierr = nf90mpi_enddef(ncid)
         call check(ierr, 'enddef')
      end if

      ! write data
      starts = (/ its, jts, kts, time_s /)
      counts = (/ ite-its+1, jte-jts+1, kte , int(1,kind=MPI_OFFSET_KIND)/)
      if (associated(data_4))then
         ierr = nf90mpi_put_var_all(ncid, varid, data_4, starts, counts)
      else if (associated(data_8)) then
         ierr = nf90mpi_put_var_all(ncid, varid, data_8, starts, counts)
      end if
      call check(ierr, 'put var')

   endsubroutine write2file3d


   subroutine write2file(varname, data_4, data_8, its, ite, jts, jte, ide, jde, time_s, ncid, dimids, units, longname)
      ! Writes data to a netcdf file
      character(len=*), intent(in) :: varname
      integer(kind=MPI_OFFSET_KIND) :: starts(3), counts(3)
      integer(kind=MPI_OFFSET_KIND), intent(in) :: its, ite, jts, jte, ide, jde, time_s
      integer, intent(in) :: ncid, dimids(3)
      integer :: varid, ierr, i
      character(len=26) :: filename
      character(len=*), intent(in) :: units, longname
      ! two dimensions data
      REAL(KIND( REAL( 0.0, 8 ) )), intent(inout), dimension(:,:), pointer :: data_8
      REAL(KIND( REAL( 0.0, 4 ) )), intent(inout), dimension(:,:), pointer :: data_4
      logical :: exists


      ! skip coordinates
      ierr = nf90mpi_inq_varid(ncid, varname, varid)
    !   call check(ierr, 'inq varid')
      if (ierr .eq. NF90_ENOTVAR) then
         ! define variable
         ierr = nf90mpi_redef(ncid)
         call check(ierr, 'redef')
         if (associated(data_4))then
            ierr = nf90mpi_def_var(ncid, varname, NF90_FLOAT, dimids, varid)
         else if (associated(data_8))then
            ierr = nf90mpi_def_var(ncid, varname, NF90_DOUBLE, dimids, varid)
         else
            write(*,*)'In write2file, no data attached to ', varname
         endif
         call check(ierr, 'def var')
         ierr = nf90mpi_put_att( ncid,\
                                 varid,\
                                 'units',\
                                 Units)
         ierr = nf90mpi_put_att( ncid,\
                                 varid,\
                                 'long_name',\
                                 LongName)
         ierr = nf90mpi_enddef(ncid)
         call check(ierr, 'enddef')
      end if

      ! write data
      starts = (/ its, jts, time_s /)
    !   write(*,*)starts
      counts = (/ ite-its+1, jte-jts+1, int(1,kind=MPI_OFFSET_KIND)/)
    !   write(*,*)counts
    !   write(*,*)'data_4: ', associated(data_4)
    !   write(*,*)'data_8: ', associated(data_8)
      if (associated(data_4))then
         ierr = nf90mpi_put_var_all(ncid, varid, data_4, starts, counts)
      else if (associated(data_8)) then
         ierr = nf90mpi_put_var_all(ncid, varid, data_8, starts, counts)
      end if
      call check(ierr, 'put var')

   endsubroutine write2file



end module wrfgc_io_pnetcdf

