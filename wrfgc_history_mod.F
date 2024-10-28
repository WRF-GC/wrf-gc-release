!  DESCRIPTION: Implement history.rc defined diagnostics output
!               for wrf-gc like history_mod.F90 in GEOS-Chem


MODULE WRFGC_History_Mod
   USE Precision_Mod
   USE HistContainer_Mod,     ONLY : HistContainer, HistContainer_Create
   USE MetaHistContainer_Mod, ONLY : MetaHistContainer
   USE Input_Opt_Mod, ONLY : OptInput
   USE MetaHistItem_Mod
   USE State_Chm_Mod ,     ONLY : ChmState
   USE State_Diag_Mod,     ONLY : DgnState
   USE State_Grid_Mod,     ONLY : GrdState
   USE State_Met_Mod,      ONLY : MetState
   USE HistItem_Mod

   USE GC_Stateful_Mod,   ONLY : GIGC_States

   IMPLICIT NONE

   type MetaHistContainerptr
      type(MetaHistContainer), pointer :: p;
   end type MetaHistContainerptr

   type MetaHistItemptr
      type(MetaHistItem), pointer :: p;
   end type MetaHistItemptr
   ! !LOCAL VARIABLES:
   !
   public :: wrfgc_history_init

   TYPE(MetaHistContainerptr), dimension(1:8), PUBLIC     :: CollectionList
   TYPE(MetaHistItemptr), dimension(1:8), PUBLIC         :: IndexVariables

   ! Scalars
   INTEGER, PUBLIC                              :: CollectionCount

   ! Strings
   CHARACTER(LEN=255),      ALLOCATABLE :: CollectionName       (:  )
   CHARACTER(LEN=255),      ALLOCATABLE :: CollectionFileName   (:  )
   CHARACTER(LEN=255),      ALLOCATABLE :: CollectionTemplate   (:  )
   CHARACTER(LEN=255),      ALLOCATABLE :: CollectionFormat     (:  )
   CHARACTER(LEN=255),      ALLOCATABLE :: CollectionFrequency  (:  )
   CHARACTER(LEN=255),      ALLOCATABLE :: CollectionAccInterval(:  )
   CHARACTER(LEN=255),      ALLOCATABLE :: CollectionDuration   (:  )
   CHARACTER(LEN=255),      ALLOCATABLE :: CollectionMode       (:  )
   CHARACTER(LEN=255),      ALLOCATABLE :: CollectionLonRange   (:  )
   CHARACTER(LEN=255),      ALLOCATABLE :: CollectionLatRange   (:  )
   INTEGER,                 ALLOCATABLE :: CollectionSubsetInd  (:,:)
   CHARACTER(LEN=255),      ALLOCATABLE :: CollectionLevels     (:  )
   INTEGER,                 ALLOCATABLE :: CollectionLevelInd   (:,:)
   CHARACTER(LEN=255),      ALLOCATABLE :: CollectionHrRange    (:  )

   ! Objects
   !
   !   !DEFINED PARAMETERS:
   !
   ! Maximum number of collections (set to a ridiculously big number)
   INTEGER,                 PARAMETER   :: MAX_COLLECTIONS = 500
   INTEGER :: N

   INTEGER,          PARAMETER, PUBLIC :: UNDEFINED_INT      = -999
   REAL(f4),         PARAMETER, PUBLIC :: UNDEFINED          = -1.0e+31_f4
   REAL(f8),         PARAMETER, PUBLIC :: UNDEFINED_DBL      = -1.0e+31_f8
   CHARACTER(LEN=9), PARAMETER, PUBLIC :: UNDEFINED_STR      = 'not found'

   INTEGER, PUBLIC, PARAMETER :: GC_SUCCESS =  0   ! Routine returns success
   INTEGER, PUBLIC, PARAMETER :: GC_FAILURE = -1   ! Routine returns failure

  


   CONTAINS
   subroutine wrfgc_history_init(Input_Opt, State_Chm, State_Diag, State_Met, State_Grid, comm)
      USE InquireMod,        ONLY : FindFreeLun
      ! USE ErrorMod,          ONLY : GC_Error
      USE Charpak_Mod
      USE History_Util_Mod
      USE MetaHistContainer_Mod
      USE state_diag_mod
      USe DiagList_Mod,      ONLY : CollList, ColItem, Search_CollList
      USE HCO_TYPES_MOD,     ONLY : DiagnCollection, DiagnCont, DiagnBundle, ConfigObj
      USE HCO_State_Mod,     ONLY : HCO_State
      USE Grid_Registry_Mod, ONLY : Lookup_Grid, Init_Grid_Registry
      
      IMPLICIT NONE

      type(OptInput),   intent(in)                 :: Input_Opt
      TYPE(ChmState),   intent(in)                 :: State_Chm  ! GEOS-Chem chemical state
      TYPE(DgnState),   intent(in)                 :: State_Diag ! GEOS-Chem diagnostic state
      TYPE(MetState),   intent(in)                 :: State_Met  ! GEOS-Chem meteorological state
      TYPE(GrdState),   intent(in)                 :: State_Grid ! GEOS-Chem grid state
      integer,          intent(in)                 :: comm

      TYPE(HistItem),     POINTER     :: Item       ! HISTORY ITEM object
      type(HistContainer), POINTER     :: Container  ! HISTORY CONTAINER object
      TYPE(ColItem), POINTER :: Current
      INTEGER     :: RC        
      !
      ! !LOCAL VARIABLES:
      !
      ! Scalars
      LOGICAL                      :: EOF,            Found
      LOGICAL                      :: FileExists
      INTEGER                      :: yyyymmdd,       hhmmss
      INTEGER                      :: yyyymmdd_end,   hhmmss_end
      INTEGER                      :: DeltaYMD,       DeltaHMS
      INTEGER                      :: X,              Y
      INTEGER                      :: C,              N,             W
      INTEGER                      :: nX,             nY,            nZ
      INTEGER                      :: fId,            IOS,           LineNum
      INTEGER                      :: nSubs1,         nSubs2
      INTEGER                      :: Ind1,           Ind2
      INTEGER                      :: UpdateYmd,      UpdateHms
      INTEGER                      :: FileCloseYmd,   FileCloseHms
      INTEGER                      :: FileWriteYmd,   FileWriteHms
      INTEGER                      :: ItemCount,      SpaceDim,      Operation
      INTEGER                      :: Ind_All,        Ind_Adv,       Ind_Aer
      INTEGER                      :: Ind_Dry,        Ind_Fix,       Ind_Gas
      INTEGER                      :: Ind_Kpp,        Ind_Pho,       Ind_Rst
      INTEGER                      :: Ind_Var,        Ind_Wet,       Ind
      INTEGER                      :: HbHrs,          HbMin,         HbSec
      INTEGER                      :: HeartBeatHms,   nTags
      integer                      :: charsize, totalsize, offset
      REAL(f8)                     :: UpdateAlarm,    HeartBeatDtSec
      REAL(f8)                     :: FileWriteAlarm, FileCloseAlarm
      REAL(f8)                     :: JulianDate,     JulianDateEnd
      REAL(f8)                     :: UpdateCheck,    FileWriteCheck
      REAL(f8)                     :: SimLengthSec

      ! Strings
      CHARACTER(LEN=6  )           :: TStr
      CHARACTER(LEN=8  )           :: DStr
      CHARACTER(LEN=20 )           :: StartTimeStamp, EndTimeStamp
      CHARACTER(LEN=63 )           :: CName
      CHARACTER(LEN=80 )           :: ErrorLine
      CHARACTER(LEN=255)           :: FileExpId
      CHARACTER(LEN=255)           :: Line,           FileName
      CHARACTER(LEN=255)           :: OutputName,     ThisLoc
      CHARACTER(LEN=255)           :: MetaData,       Reference
      CHARACTER(LEN=255)           :: Title,          Units
      CHARACTER(LEN=255)           :: ItemTemplate,   ItemTemplateUC
      CHARACTER(LEN=255)           :: ItemName,       Description
      CHARACTER(LEN=255)           :: TmpMode,        Contact
      CHARACTER(LEN=255)           :: Pattern,        ItemPrefix
      CHARACTER(LEN=255)           :: tagId,          tagName
      CHARACTER(LEN=512)           :: ErrMsg,         FileMsg

      REAL(f8)                     :: Subset(2)
      INTEGER                      :: Levels(200)
      CHARACTER(LEN=255)           :: Subs1(255)
      CHARACTER(LEN=255)           :: Subs2(255)
      CHARACTER(LEN=255)           :: SubStrs(255)


      ! Strings
      CHARACTER(LEN=20)        :: ItemDimNames(13)
      CHARACTER(LEN=20)        :: ItemNames(13)
      CHARACTER(LEN=20)        :: RegistryNames(13)

      ! Pointer arrays
      REAL(f8),        POINTER :: Ptr0d_8
      REAL(f8),        POINTER :: Ptr1d_8(:  )
      REAL(f8),        POINTER :: Ptr2d_8(:,:)
      REAL(f4),        POINTER :: Ptr2d_4(:,:)

      LOGICAL                  :: OnLevelEdges
      INTEGER                  :: Output_KindVal
      INTEGER                  :: Source_KindVal
      INTEGER                  :: Rank
      INTEGER                  :: nILev
      INTEGER                  :: nLev

      INTEGER                  :: Dimensions(3)
      INTEGER                  :: Subset_X(2), Subset_Xc(2)
      INTEGER                  :: Subset_Y(2), Subset_Yc(2)
      INTEGER                  :: Subset_Z(2), Subset_Zc(2), Subset_Ze(2)

      
      

      ! initialize variables
      EOF            =  .FALSE.
      offset         =  0
      IOS            =  0
      UpdateYmd      =  0
      UpdateHms      =  0
      FileCloseYmd   =  0
      FileCloseHms   =  0
      FileWriteYmd   =  0
      FileWriteHms   =  0
      LineNum        =  0
      SpaceDim       =  0
      HeartBeatDtSec =  DBLE( Input_Opt%TS_DYN )
      yyyymmdd       =  Input_Opt%NymdB
      hhmmss         =  Input_Opt%NhmsB
      yyyymmdd_end   =  Input_Opt%NymdE
      hhmmss_end     =  Input_Opt%NhmsE
      Container => NULL()
      Item      => NULL()
      C              =  0
      Ptr0d_8 => NULL()
      Ptr1d_8 => NULL()
      Ptr2d_8 => NULL()
      Ptr2d_4 => NULL()
      
      ! do N = 1, 8
      !    allocate(IndexVariables(State_Grid%ID));
      !    allocate(CollectionList(State_Grid%ID));
      ! enddo

      ThisLoc = " -> at wrfgc_history_init (in module wfrcg_history_mod.F)"

      if(State_Grid%id == 1)then
         call Init_Grid_Registry(Input_Opt, State_Grid, RC)
         call checkgc(RC, 'Error initializing grid registry'//ThisLoc )
      endif

      CALL Compute_DeltaYmdHms_For_End( yyyymmdd,     hhmmss,               &
         yyyymmdd_end, hhmmss_end,           &
         deltaYMD,     deltaHMS             )

      ! Initialize the collection list
      ! read collection : 1. collection name,
      !                       1). template, 2). frequency, 3). duration, 4). mode, 5). fields
      call get_collection_name()

      fId = FindFreeLun()
      OPEN( fId, FILE=TRIM(Input_Opt%HistoryInputFile), STATUS='OLD', IOSTAT=RC )
      ErrMsg = 'Error opening "' //TRIM(Input_Opt%HistoryInputFile) // '"!'
      CALL checkgc( RC, ThisLoc//ErrMsg)

      do
500      continue
         Line = ReadOneLine(fId, EOF, IOS, Squeeze = .TRUE.)
         LineNum = LineNum + 1

         if (EOF) exit

         IF ( IOS > 0 ) THEN
            ErrMsg = 'Unexpected end-of-file in "'                          // &
               TRIM( Input_Opt%HistoryInputFile )
            WRITE( ErrorLine, 250 ) LineNum
250         FORMAT( ' -> ERROR occurred at (or near) line ', i6,               &
               ' of the HISTORY.rc file' )
            CALL checkgc( RC, ErrMsg//ThisLoc//ErrorLine )
            RETURN
         ENDIF

         ! skip commented line
         if(Line(1:1) == '#') cycle

         !====================================================================
         ! Get the EXPID string.  This is the "front part" of the netCDF
         ! file path for each collection.  In other words, if EXPID is
         ! "OutputDir/GEOSChem", then the default SpeciesConc collection file
         ! names will be "OutputDir/GEOSChem.SpeciesConc_YYYYMMDD_hhmmz.nc4"
         !====================================================================
         IF ( INDEX( Line, 'EXPID' ) > 0  ) THEN

            ! Split the line on the colon
            CALL StrSplit( Line, ":", Subs1, nSubs1 )

            ! Stop with error if there are more than 2 substrings
            IF ( nSubs1 /= 2 ) THEN
               ErrMsg = 'Error in extracting the EXPID value from the '     // &
                  'HISTORY.rc file.  This forms the start of the '    // &
                  'netCDF file name for each collection.  Please '    // &
                  'check the HISTORY.rc file for typos.'
               WRITE( ErrorLine, 250 ) LineNum
               CALL checkgc(  RC, ThisLoc//ErrorLine //ErrMsg)
            ENDIF

            ! Save the EXPID parameter
            FileExpId = Subs1(2)
            CALL CStrip( FileExpId )
         ENDIF

         ! "template": Specifies the year/month/day/hr/min/sec in filenames
         ! Can be omitted if "filename" is specified
         Pattern = 'template'
         IF ( INDEX( TRIM( Line ), TRIM( Pattern ) ) > 0 ) THEN
            CALL GetCollectionMetaData( Input_Opt, Line, Pattern, MetaData, C )
            IF ( C > 0 ) then
               ! write(*,*)"get template of collection: ", TRIM(CollectionName(C))
               CollectionTemplate(C) = Metadata
            endif
         ENDIF

         ! "format": Specifies the file output format (e.g. netCDF-4, CFIO)
         Pattern = 'format'
         IF ( INDEX( TRIM( Line ), TRIM( Pattern ) ) > 0 ) THEN
            CALL GetCollectionMetaData( Input_Opt, Line, Pattern, MetaData, C )
            IF ( C > 0 )then
               ! write(*,*)'get format of collection: ', TRIM(CollectionName(C))
               CollectionFormat(C) = Metadata
            endif
         ENDIF

         ! "frequency": Specifies how often diagnostics are updated,
         ! Must be either in "YYYYMMDD hhmmss" or "hhmmss" format.
         Pattern = 'frequency'
         IF ( INDEX( TRIM( Line ), TRIM( Pattern ) ) > 0 ) THEN
            CALL GetCollectionMetaData( Input_Opt, Line, Pattern, MetaData, C )
            IF ( C > 0 ) THEN
               IF ( LEN_TRIM( MetaData ) == 6     .or.                         &
                  LEN_TRIM( MetaData ) == 14    .or.                         &
                  TRIM(     MetaData ) == 'End' .or.                         &
                  TRIM(     MetaData ) == 'end' ) THEN
                  CollectionFrequency(C) = Metadata
                  ! write(*,*)'get frequency of collection: ', TRIM(CollectionName(C))
               ELSE
                  ErrMsg = 'Error in defining "frequency" for collection "' // &
                     TRIM( CollectionName(C) ) // '"!  This field '   // &
                     'must either be of the format "YYYYMMDD '        // &
                     'hhmmss", "hhmmss", or "End".  Please check the '// &
                     '"frequency" setting in the HISTORY.rc file.'
                  WRITE( ErrorLine, 250 ) LineNum
                  CALL checkgc(  RC, ThisLoc//ErrorLine //ErrMsg)
                  RETURN
               ENDIF
            ENDIF
         ENDIF

         Pattern = 'duration'
         IF ( INDEX( TRIM( Line ), TRIM( Pattern ) ) > 0 ) THEN
            CALL GetCollectionMetaData( Input_Opt, Line, Pattern, MetaData, C )
            IF ( C > 0 ) THEN
               IF ( LEN_TRIM( MetaData ) == 6     .or.                         &
                  LEN_TRIM( MetaData ) == 14    .or.                         &
                  TRIM(     MetaData ) == 'End' .or.                         &
                  TRIM(     MetaData ) == 'end' ) THEN
                  CollectionDuration(C) = Metadata
                  ! write(*,*)'get duration of collection: ', TRIM(CollectionName(C))
               ELSE
                  ErrMsg = 'Error in defining "duration" for collection "'  // &
                     TRIM( CollectionName(C) ) // '"!  This field '   // &
                     'must either be of the format "YYYYMMDD '        // &
                     'hhmmss", "hhmmss", or "End".  Please check the '// &
                     '"duration" setting in the HISTORY.rc file.'
                  WRITE( ErrorLine, 250 ) LineNum
                  CALL checkgc(  RC, ThisLoc//ErrorLine //ErrMsg)
                  RETURN
               ENDIF
            ENDIF
         ENDIF

         Pattern = 'mode'
         IF ( INDEX( TRIM( Line ), TRIM( Pattern ) ) > 0 ) THEN
            CALL GetCollectionMetaData( Input_Opt, Line, Pattern, MetaData, C )
            IF ( C > 0 ) THEN
               TmpMode = Metadata
               ! write(*,*)'get mode of collection: ', TRIM(CollectionName(C))
               CALL TranUc( TmpMode )
               SELECT CASE( TmpMode )
                CASE( 'INSTANTANEOUS', 'TIME-AVERAGED', 'TIMEAVERAGED' )
                  CollectionMode(C) = Metadata
                CASE DEFAULT
                  ErrMsg = 'Error in defining "mode" for collection "'   // &
                     TRIM( CollectionName(C) ) // '"!  The mode ' // &
                     'value can either be "instantaneous" or '     // &
                     '"time-averaged".  Please check the "mode" '  // &
                     'setting in the HISTORY.rc file.'
                  WRITE( ErrorLine, 250 ) LineNum
                  CALL checkgc(  RC, ThisLoc//ErrorLine //ErrMsg)
                  RETURN
               END SELECT
            ENDIF
         ENDIF


         if (INDEX(trim(Line), 'fields') > 0) then
            if (C>0)then
               ! write(*,*)'get fields of collection: ', TRIM(CollectionName(C))
               ! write(*,*)'Line: ', trim(Line),                           &
               !    'CollectionName: ', TRIM(CollectionName(C)),           &
               !    'CollectionFileName: ', TRIM(CollectionFileName(C)),   &
               !    'CollectionTemplate: ', TRIM(CollectionTemplate(C)),   &
               !    'CollectionFormat: ', TRIM(CollectionFormat(C)),       &
               !    'CollectionFrequency: ', TRIM(CollectionFrequency(C)), &
               !    'CollectionAccInterval: ', TRIM(CollectionAccInterval(C)), &
               !    'CollectionDuration: ', TRIM(CollectionDuration(C)),   &
               !    'CollectionMode: ', TRIM(CollectionMode(C))           


               ! Frequecy indicates how often the data will be writen to file
               IF ( LEN_TRIM( CollectionFrequency(C) ) == 6 ) THEN
                  READ( CollectionFrequency(C), '(i6.6)'  ) FileWriteHms
               ELSE IF ( LEN_TRIM( CollectionFrequency(C) ) == 14 ) THEN
                  READ( CollectionFrequency(C), '(i8,i6)' ) FileWriteYmd,         &
                     FileWriteHms
               ELSE IF ( TRIM( CollectionFrequency(C) ) == 'End'   .or.         &
                  TRIM( CollectionFrequency(C) ) == 'end' ) THEN
                  FileWriteYmd = DeltaYMD
                  FileWriteHms = DeltaHMS
               ENDIF

               IF ( FileWriteHms == 240000 ) THEN
                  FileWriteYmd = 00000001
                  FileWriteHms = 000000
               ENDIF

               UpdateYmd = FileWriteYmd
               UpdateHms = FileWriteHms

               ! Duration indicates how often the new file will be created
               IF ( TRIM( CollectionDuration(C) ) == UNDEFINED_STR ) THEN
                  FileCloseYmd = FileWriteYmd
                  FileCloseHms = FileWriteHms
               ELSE IF ( LEN_TRIM( CollectionDuration(C) ) == 6 ) THEN
                  READ( CollectionDuration(C), '(i6.6)'  ) FileCloseHms
               ELSE IF ( LEN_TRIM( CollectionDuration(C) ) == 14 ) THEN
                  READ( CollectionDuration(C), '(i8,i6)' ) FileCloseYmd,          &
                     FileCloseHms
               ELSE IF ( TRIM( CollectionDuration(C) ) == 'End'   .or.          &
                  TRIM( CollectionDuration(C) ) == 'end' ) THEN
                  FileCloseYmd = DeltaYMD
                  FileCloseHms = DeltaHMS
               ENDIF

               call HistContainer_Create( Input_Opt = Input_Opt,           &
                  Container  = Container,          &
                  id         = C,                  &
                  name       = CollectionName(C),  &
                  CurrentYmd = yyyymmdd,           &
                  CurrentHms = hhmmss,             &
                  UpdateMode = CollectionMode(C),  &
                  UpdateYmd  = UpdateYmd,          &
                  UpdateHms  = UpdateHms,          &
                  FileWriteYmd = FileWriteYmd,     &
                  FileWriteHms = FileWriteHms,     &
                  FileCloseYmd = FileCloseYmd,     &
                  FileCloseHms = FileCloseHms,     &
                  FileName   = CollectionFileName(C),&
                  FileTemplate = CollectionTemplate(C),&
                  FileExpId  = FileExpId,          &
                  RC         = RC                  )
            endif



            ! get collection name
            N = INDEX( Trim(Line), '.' )
            CName = Line(1:N-1)
            ! skip all fileds if collection name is not in the list
            CALL Search_CollList( Input_Opt%amIRoot, CollList, CName, Found, RC )
            IF ( .not. Found ) then
               do
                  Line = ReadOneLine(fId, EOF, IOS, Squeeze = .TRUE.)
                  LineNum = LineNum + 1
                  call StrSqueeze(Line)
                  if (TRIM(Line) == '::') GOTO 500
               end do
            END IF

            ! Zero the counter of items
            ItemCount = 0

            ! Create title string for collection
            Title = 'GEOS-Chem diagnostic collection: ' //                     &
               TRIM( CollectionName(C) )
            ! read item name
            do
               if (ItemCount == 0)then
                  call GetCollectionMetaData(Input_Opt, Line, 'fields', MetaData, C)
                  call StrSplit(MetaData, " ", Subs1, nSubs1)
                  ItemName = Subs1(1)
               else
                  Line = ReadOneLine(fId, EOF, IOS, Squeeze = .TRUE.)
                  LineNum = LineNum + 1
                  Line = CleanText(Line)
                  if(Line(1:1) == "#") cycle
                  if (Line(1:2)== '::') exit
                  if (EOF) exit
                  if (IOS > 0) then
                     ErrMsg = 'Unexpected end-of-file in "'                          // &
                        TRIM( Input_Opt%HistoryInputFile )
                     WRITE( ErrorLine, 250 ) LineNum
                     CALL checkgc( RC, ErrMsg//ThisLoc//ErrorLine )
                     RETURN
                  ENDIF

                  call StrSplit(Line, " ", Subs1, nSubs1)
                  ItemName = Subs1(1)
               end if

               ItemTemplate = ItemName
               ItemTemplateUC = To_UpperCase(ItemTemplate)
               if (index(ItemTemplate, "?")>0) then
                  call StrSplit(ItemTemplate, '?', SubStrs, N)
                  tagId = SubStrs(N-1)
                  ItemPrefix = SubStrs(1)

                  ! Get number of tags for this wildcard
                  CALL Get_TagInfo( Input_Opt, tagId, State_Chm, Found, RC, &
                     nTags=nTags )
                  ErrMsg = 'Error retrieving # of tags for' //              &
                     ' wildcard ' // TRIM(tagId)
                  call checkgc(RC, ErrMsg)

                  ! add tag to the item
                  do N = 1, nTags
                     call Get_TagInfo(Input_Opt, tagId, State_Chm, Found, RC, &
                        N=N, tagName=tagName)
                     ErrMsg = 'Error retreving tag name for wildcard '//trim(tagId)
                     call checkgc(RC, ErrMsg)
                     ItemName = trim(ItemPrefix)//trim(tagName)

                     call Get_NameInfo(Input_Opt, ItemName, OutputName, RC)

                     ItemCount = ItemCount + 1
                     call History_AddItemToCollection(Input_Opt   = Input_Opt,         &
                        State_Chm   = State_Chm,         &
                        State_Diag  = State_Diag,        &
                        State_Met   = State_Met,         &
                        ItemName    = OutputName,        &
                        ItemCount   = ItemCount,         &
                        Collection  = Container,         &
                        CollectionId = C,                &
                        RC=RC                            &
                        )
                     ErrMsg = 'Could not add diagnostic "'               // &
                        TRIM( OutputName ) // '" to collection: '  // &
                        TRIM( CollectionName(C) )
                  enddo
               else
                  ! the item has no wildcard
                  call Get_NameInfo(Input_Opt, ItemName, OutputName, RC)

                  ItemCount = ItemCount + 1

                  call History_AddItemToCollection(Input_Opt   = Input_Opt,         &
                     State_Chm   = State_Chm,         &
                     State_Diag  = State_Diag,        &
                     State_Met   = State_Met,         &
                     ItemName    = OutputName,        &
                     ItemCount   = ItemCount,         &
                     Collection  = Container,         &
                     CollectionId = C,                &
                     RC=RC                            &
                     )
                  ErrMsg = 'Could not add diagnostic "'               // &
                     TRIM( OutputName ) // '" to collection: '  // &
                     TRIM( CollectionName(C) )
               endif

            enddo

            
            
            call MetaHistContainer_AddNew(   Input_Opt   = Input_Opt,      &
            Node        = CollectionList(State_Grid%ID)%p, &
            Container   = Container,      &
            RC          = RC              )
            ErrMsg = 'Could not add Container' //                           &
            TRIM( CollectionName(C) ) //                           &
            ' to the list of collections!'
            call checkgc(RC, ErrMsg)
            
            C = 0
         endif
         
      enddo
      
      ! add HISTORY ITEM for each netCDF index variable (e.g.
      !  lon, lat, lev, time, area, hyam, hybm, p0)
      ! Fields saved in the Registry object in GeosUtil/grid_registry_mod.F90
      ! The registry object current only works for domain 1
      if (State_Grid%ID == 1) then
         RegistryNames(1 ) = 'GRID_AREA'
         RegistryNames(2 ) = 'GRID_P0'
         RegistryNames(3 ) = 'GRID_HYBI'
         RegistryNames(4 ) = 'GRID_HYAI'
         RegistryNames(5 ) = 'GRID_HYBM'
         RegistryNames(6 ) = 'GRID_HYAM'
         RegistryNames(7 ) = 'GRID_LON'
         RegistryNames(8 ) = 'GRID_LONBND'
         RegistryNames(9 ) = 'GRID_LAT'
         RegistryNames(10) = 'GRID_LATBND'
         RegistryNames(11) = 'GRID_ILEV'
         RegistryNames(12) = 'GRID_LEV'
         RegistryNames(13) = 'GRID_TIME'

         ! Name for each HISTORY ITEM
         ItemNames(1 )     = 'AREA'
         ItemNames(2 )     = 'P0'
         ItemNames(3 )     = 'hybi'
         ItemNames(4 )     = 'hyai'
         ItemNames(5 )     = 'hybm'
         ItemNames(6 )     = 'hyam'
         ItemNames(7 )     = 'lon'
         ItemNames(8 )     = 'lon_bnds'
         ItemNames(9 )     = 'lat'
         ItemNames(10)     = 'lat_bnds'
         ItemNames(11)     = 'ilev'
         ItemNames(12)     = 'lev'
         ItemNames(13)     = 'time'

         ! Dimensions for each HISTORY ITEM
         ItemDimNames(1 )  = 'xy'
         ItemDimNames(2 )  = '-'
         ItemDimNames(3 )  = 'z'
         ItemDimNames(4 )  = 'z'
         ItemDimNames(5 )  = 'z'
         ItemDimNames(6 )  = 'z'
         ItemDimNames(7 )  = 'x'
         ItemDimNames(8 )  = 'bx'
         ItemDimNames(9 )  = 'y'
         ItemDimNames(10)  = 'by'
         ItemDimNames(11)  = 'z'
         ItemDimNames(12)  = 'z'
         ItemDimNames(13)  = 't'

         DO N = 2, SIZE(RegistryNames)
            !---------------------------------------------------------------------
            ! Look up one of the index fields from grid_registry_mod.F90
            !---------------------------------------------------------------------
            CALL Lookup_Grid( Input_Opt      = Input_Opt,                         &
                              Variable       = RegistryNames(N),                   &
                              Description    = Description,                       &
                              Dimensions     = Dimensions,                        &
                              Source_KindVal = Source_KindVal,                    &
                              Output_KindVal = Output_KindVal,                    &
                              Rank           = Rank,                              &
                              Units          = Units,                             &
                              OnLevelEdges   = OnLevelEdges,                      &
                              Ptr0d_8        = Ptr0d_8,                           &
                              Ptr1d_8        = Ptr1d_8,                           &
                              Ptr2d_8        = Ptr2d_8,                           &
                              Ptr2d_4        = Ptr2d_4,                           &
                              RC             = RC                                )
            ErrMsg = 'Error looking up grid variable "' // TRIM( RegistryNames(N) ) // '"'
            CALL checkgc( RC, ErrMsg )

            CALL Get_Number_Of_Levels( Container, nLev, nIlev )

            ! Subset indices
            Subset_Xc = (/CollectionList(State_Grid%ID)%p%Container%X0,CollectionList(State_Grid%ID)%p%Container%X1/)
            Subset_Yc = (/ CollectionList(State_Grid%ID)%p%Container%Y0,CollectionList(State_Grid%ID)%p%Container%Y1 /)
            Subset_Zc = (/ CollectionList(State_Grid%ID)%p%Container%Z0, nLev         /)
            Subset_Ze = (/ CollectionList(State_Grid%ID)%p%Container%Z0, nILev        /)

            SELECT CASE( N )
               CASE( 2 )                   ! P0
                  Subset_X = (/ 0, 0 /)
                  Subset_Y = (/ 0, 0 /)
                  Subset_Z = (/ 0, 0 /)
               CASE( 3, 4, 11 )            ! hybi, hyai, ilev
                  Subset_X = Subset_Xc
                  Subset_Y = Subset_Yc
                  Subset_Z = Subset_Ze
               CASE( 8  )                  ! lon_bnds
                  Subset_X = (/ 1, 2 /)
                  Subset_Y = Subset_Xc
                  Subset_Z = (/ 0, 0 /)
               CASE( 10 )                  ! lat_bnds
                  Subset_X = (/ 1, 2 /)
                  Subset_Y = Subset_Yc
                  Subset_Z = (/ 0, 0 /)
               CASE DEFAULT                ! everything else
                  Subset_X = Subset_Xc
                  Subset_Y = Subset_Yc
                  Subset_Z = Subset_Zc
            END SELECT

            write(*,*)'dingao debug: ', ItemNames(N)      
            !---------------------------------------------------------------------
            ! Create a HISTORY ITEM for this index field
            !---------------------------------------------------------------------
         if(N==7)then 
               write(*,*)N
               write(*,*)Ptr1d_8
         endif
            CALL HistItem_Create( Input_Opt      = Input_Opt,                     &
            Item           = Item,                          &
            Id             = N,                             &
            ContainerId    = 0,                             &
            Name           = ItemNames(N),                   &
            LongName       = Description,                   &
            Units          = Units,                         &
            SpaceDim       = Rank,                          &
            OnLevelEdges   = OnLevelEdges,                  &
            DimNames       = ItemDimNames(N),                &
            Operation      = 0,                             &
            Subset_X       = Subset_X,                      &
            Subset_Y       = Subset_Y,                      &
            Subset_Z       = Subset_Z,                      &
            Source_KindVal = Source_KindVal,                &
            Output_KindVal = Output_KindVal,                &
            Source_0d_8    = Ptr0d_8,                       &
            Source_1d_8    = Ptr1d_8,                       &
            Source_2d_4    = Ptr2d_4,                       &
            Source_2d_8    = Ptr2d_8,                       &
            RC             = RC                            )
            ErrMsg = 'Error creating HISTORY ITEM for grid variable "' // &
            TRIM( RegistryNames(N) ) // '"'
            CALL checkgc( RC, ErrMsg )

            

            call MetaHistItem_AddNew( Input_Opt   = Input_Opt,             &
                                    Item        = Item,                    &
                                    Node        = IndexVariables(State_Grid%ID)%p,     &
                                    RC          = RC                        )
            ErrMsg = 'Could not add HISTORY ITEM for grid variable "' // &
            TRIM( RegistryNames(N) ) // '" to the list of HISTORY ITEMS!'
            call checkgc(RC, ErrMsg)
         enddo 
      endif
      !============================
      ! clean up and quit
      !============================
      close(fId)
      ! Write output
      IF ( Input_Opt%amIRoot ) THEN
         WRITE( 6, '(/,a)' ) REPEAT( '=', 79 )
         WRITE( 6, '(a  )' ) 'DEFINED DIAGNOSTIC COLLECTIONS:'
         WRITE( 6, '(a  )' ) REPEAT( '=', 79 )

         DO C = 1, CollectionCount
            print*, 'Collection        ', TRIM( CollectionName       (C) )
            print*, '  -> FileName     ', TRIM( CollectionFileName   (C) )
            print*, '  -> Format       ', TRIM( CollectionFormat     (C) )
            print*, '  -> Frequency    ', TRIM( CollectionFrequency  (C) )
            IF ( CollectionAccInterval(C) /= UNDEFINED_STR ) THEN
               print*, '  -> Acc_Interval ', TRIM( CollectionAccInterval(C) )
            ENDIF
            IF ( CollectionDuration(C) /= UNDEFINED_STR ) THEN
               print*, '  -> Duration     ', TRIM( CollectionDuration   (C) )
            ENDIF
            print*, '  -> Mode         ', TRIM( CollectionMode       (C) )
            IF ( CollectionLonRange(C) /= UNDEFINED_STR ) THEN
               print*, '  -> LON_RANGE    ',                                   &
                  TRIM(ADJUSTL(ADJUSTR( CollectionLonRange(C) )))
               print*, '     -> X0 X1  ', ((CollectionSubsetInd(N,C)), N=1,2)
            ENDIF
            IF ( CollectionLatRange(C) /= UNDEFINED_STR ) THEN
               print*, '  -> LAT_RANGE    ',                                   &
                  TRIM(ADJUSTL(ADJUSTR( CollectionLatRange(C) )))
               print*, '     -> Y0 Y1  ', ((CollectionSubsetInd(N,C)), N=3,4)
            ENDIF
            IF ( CollectionLevels(C) /= UNDEFINED_STR ) THEN
               print*, '  -> Levels    ' , TRIM( CollectionLevels(C) )
               print*, '     -> Z0 Z1  ', ((CollectionLevelInd(N,C)), N=1,2)
            ENDIF
            IF ( CollectionHrRange(C) /= UNDEFINED_STR ) THEN
               print*, '  -> hrrange  ' , TRIM( CollectionHrRange(C) )
            ENDIF

            ! Trap error if the collection frequency is undefined
            ! This indicates an error in parsing the file
            IF ( TRIM( CollectionFrequency(C) ) == UNDEFINED_STR ) THEN
               ErrMsg = 'Collection: ' // TRIM( CollectionName(C) ) //         &
                  ' is undefined!'
               WRITE( ErrorLine, 250 ) LineNum
               call checkgc(  RC, ThisLoc//ErrorLine //ErrMsg)
               RETURN
            ENDIF
         ENDDO
      ENDIF

      IF (Input_Opt%amIRoot ) THEN
         ! Print information about each diagnostic collection
         CALL MetaHistContainer_Print( Input_Opt, CollectionList(State_Grid%ID)%p, RC )
      ENDIF

      ! Write spacer
      WRITE( 6, '(a,/)' ) REPEAT( '=', 79 )

   end subroutine wrfgc_history_init

   ! !IROUTINE: GetCollectionMetaData
!
! !DESCRIPTION: Parses a line of the HISTORY.rc file and returns metadata
!  for a given attribute (e.g. "frequency", "template", etc.)
!\\
!\\
! !INTERFACE:
!
   SUBROUTINE GetCollectionMetaData( Input_Opt, Line, Pattern, MetaData, &
      nCollection )
!
! !USES:
!
      USE Charpak_Mod,      ONLY : CleanText, StrSplit
      USE DiagList_Mod,     ONLY : CollList,  Search_CollList
      USE History_Util_Mod
      USE Input_Opt_Mod,    ONLY : OptInput
!
! !INPUT PARAMETERS:
!
      TYPE(OptInput),     INTENT(IN)  :: Input_Opt     ! Input Options object
      CHARACTER(LEN=*),   INTENT(IN)  :: Line          ! Line to be searched
      CHARACTER(LEN=*),   INTENT(IN)  :: Pattern       ! Search pattern
!
! !OUTPUT PARAMETERS:
!
      CHARACTER(LEN=255), INTENT(OUT) :: MetaData      ! Metadata value
      INTEGER,            INTENT(OUT) :: nCollection   ! Collection Id
!
! !REVISION HISTORY:
!  16 Aug 2017 - R. Yantosca - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      ! Scalars
      LOGICAL                  :: Found
      INTEGER                  :: C, Ind, nSubStr, N, P, RC

      ! Strings
      CHARACTER(LEN=255)       :: Name
      CHARACTER(LEN=255)       :: SubStr(255)

      !=======================================================================
      ! Initialize
      !=======================================================================
      nCollection = UNDEFINED_INT
      MetaData    = UNDEFINED_STR

      !=======================================================================
      ! Find the metadata for the given collection
      !=======================================================================

      ! The collection name is between column 1 and the first "." character
      Ind  = INDEX( TRIM( Line ), '.' )
      Name = Line(1:Ind-1)

      ! Exit if the collection name is not in the list of active collections
      CALL Search_CollList( Input_Opt%amIRoot, CollList, Name, Found, RC )
      IF ( .not. Found ) RETURN

      ! Non-white-space lengths of the collection cname and search pattern
      N = LEN_TRIM( Name    )
      P = LEN_TRIM( Pattern )

      ! Loop over all collection names
      ! NOTE: This algorithm may not be the most efficient, as it will
      ! not skip collections that we have already encountered.  But it
      ! only gets done during the init phase, so it might not be a huge
      ! expenditure of time anyway.  Worry about this later. (bmy, 1/18/18)
      DO C = 1, CollectionCount

         ! Check to see if the current line matches the collection name
         ! Then check to see which collection this is in
         IF ( Name(1:30) == CollectionName(C)(1:30) ) THEN
            Ind = 1
         ELSE
            Ind = 0
         ENDIF

         ! If the we match the current collection, then ...
         IF ( Ind > 0 ) THEN

            ! Split the line on the the colon
            CALL StrSplit( Line, ':', SubStr, nSubStr )

            ! If there are 2 substrings ...
            IF ( nSubStr == 2 ) THEN

               ! Make sure the first substring matches the name
               ! of the metadata field we would like to obtain.
               ! if it does, then we have found a match, and so return
               IF ( SubStr(1)(N+2:P+N+1) == Pattern(1:P) ) THEN
                  nCollection = C
                  MetaData    = CleanText( SubStr(2) )
                  EXIT
               ENDIF
            ENDIF
         ENDIF
      ENDDO

   END SUBROUTINE GetCollectionMetaData


   ! read collection name into CollectionName and CollectionCount
   subroutine get_collection_name()
      USe DiagList_Mod,      ONLY : CollList, ColItem
      TYPE(ColItem), POINTER :: Current

      CHARACTER(LEN=255)     :: TmpCollectionName(MAX_COLLECTIONS)
      CHARACTER(LEN=255)     :: ErrMsg, ThisLoc, Line,  Line2
      INTEGER                :: RC

      ! Zero local variables
      RC      = GC_SUCCESS
      ErrMsg  = 'No error'
      ThisLoc = &
         ' -> at get_collection_name (in module wfrcg_history_mod.F)'


      CollectionCount = 0
      Current => CollList%head
      do while(associated(Current))
         CollectionCount = CollectionCount + 1
         TmpCollectionName(CollectionCount) = Current%cname
         Current => Current%next
      end do
      Current => NULL()
      if ( .not. allocated(CollectionName) ) then
         allocate(CollectionName(CollectionCount), STAT=RC)
         call checkgc(RC, 'Could not allocate CollectionName!'//ThisLoc)
      end if
      ! Copy the collection names from the temporary array
      ! NOTE: The linked-list stors in reverse order, so revert here
      DO N = 1, CollectionCount
         CollectionName(N) = TmpCollectionName(CollectionCount-N+1)
      ENDDO
      write(*,*)'collection name:', ThisLoc
      DO N = 1, CollectionCount
         write(*,*)CollectionName(N)
      enddo
      ! Allocate CollectionFileName
      IF ( .not. ALLOCATED( CollectionFileName ) ) THEN
         ALLOCATE( CollectionFileName( CollectionCount ), STAT=RC )
         IF ( RC /= GC_SUCCESS ) THEN
            ErrMsg = 'Could not allocate CollectionFileName!'
            CALL checkgc( RC, ErrMsg//ThisLoc )
            RETURN
         ENDIF
         CollectionFileName = UNDEFINED_STR
      ENDIF

      ! Allocate CollectionTemplate
      IF ( .not. ALLOCATED( CollectionTemplate ) ) THEN
         ALLOCATE( CollectionTemplate( CollectionCount ), STAT=RC )
         IF ( RC /= GC_SUCCESS ) THEN
            ErrMsg = 'Could not allocate CollectionTemplate!'
            CALL checkgc( RC, ErrMsg//ThisLoc )
            RETURN
         ENDIF
         CollectionTemplate = UNDEFINED_STR
      ENDIF

      ! Allocate CollectionFormat
      IF ( .not. ALLOCATED( CollectionFormat ) ) THEN
         ALLOCATE( CollectionFormat( CollectionCount ), STAT=RC )
         IF ( RC /= GC_SUCCESS ) THEN
            ErrMsg = 'Could not allocate CollectionFormat!'
            CALL checkgc( RC, ErrMsg//ThisLoc )
            RETURN
         ENDIF
         CollectionFormat = UNDEFINED_STR
      ENDIF

      ! Allocate CollectionFrequency
      IF ( .not. ALLOCATED( CollectionFrequency ) ) THEN
         ALLOCATE( CollectionFrequency( CollectionCount ), STAT=RC )
         IF ( RC /= GC_SUCCESS ) THEN
            ErrMsg = 'Could not allocate CollectionFrequency!'
            CALL checkgc( RC, ErrMsg//ThisLoc )
            RETURN
         ENDIF
         CollectionFrequency = UNDEFINED_STR
      ENDIF

      ! Allocate CollectionAccInterval
      IF ( .not. ALLOCATED( CollectionAccInterval ) ) THEN
         ALLOCATE( CollectionAccInterval( CollectionCount ), STAT=RC )
         IF ( RC /= GC_SUCCESS ) THEN
            ErrMsg = 'Could not allocate CollectionAccInterval!'
            CALL checkgc( RC, ErrMsg//ThisLoc )
            RETURN
         ENDIF
         CollectionAccInterval = UNDEFINED_STR
      ENDIF

      ! Allocate CollectionDuration
      IF ( .not. ALLOCATED( CollectionDuration ) ) THEN
         ALLOCATE( CollectionDuration( CollectionCount ), STAT=RC )
         IF ( RC /= GC_SUCCESS ) THEN
            ErrMsg = 'Could not allocate CollectionDuration!'
            CALL checkgc( RC, ErrMsg//ThisLoc )
            RETURN
         ENDIF
         CollectionDuration = UNDEFINED_STR
      ENDIF

      ! Allocate CollectionMode
      IF ( .not. ALLOCATED( CollectionMode ) ) THEN
         ALLOCATE( CollectionMode( CollectionCount ), STAT=RC )
         IF ( RC /= GC_SUCCESS ) THEN
            ErrMsg = 'Could not allocate CollectionMode!'
            CALL checkgc( RC, ErrMsg//ThisLoc )
            RETURN
         ENDIF
         CollectionMode = UNDEFINED_STR
      ENDIF

      ! Allocate CollectionLonRange
      IF ( .not. ALLOCATED( CollectionLonRange ) ) THEN
         ALLOCATE( CollectionLonRange( CollectionCount ), STAT=RC )
         IF ( RC /= GC_SUCCESS ) THEN
            ErrMsg = 'Could not allocate CollectionLonRange!'
            CALL checkgc( RC, ErrMsg//ThisLoc )
            RETURN
         ENDIF
         CollectionLonRange = UNDEFINED_STR
      ENDIF

      ! Allocate CollectionLatRange
      IF ( .not. ALLOCATED( CollectionLatRange ) ) THEN
         ALLOCATE( CollectionLatRange( CollectionCount ), STAT=RC )
         IF ( RC /= GC_SUCCESS ) THEN
            ErrMsg = 'Could not allocate CollectionLatRange!'
            CALL checkgc( RC, ErrMsg//ThisLoc )
            RETURN
         ENDIF
         CollectionLatRange = UNDEFINED_STR
      ENDIF

      ! Allocate CollectionSubsetInd
      IF ( .not. ALLOCATED( CollectionSubsetInd ) ) THEN
         ALLOCATE( CollectionSubsetInd( 4, CollectionCount ), STAT=RC )
         IF ( RC /= GC_SUCCESS ) THEN
            ErrMsg = 'Could not allocate CollectionSubsetInd!'
            CALL checkgc( RC, ErrMsg//ThisLoc )
            RETURN
         ENDIF
         CollectionSubsetInd = UNDEFINED_INT
      ENDIF

      ! Allocate CollectionLevels
      IF ( .not. ALLOCATED( CollectionLevels ) ) THEN
         ALLOCATE( CollectionLevels( CollectionCount ), STAT=RC )
         IF ( RC /= GC_SUCCESS ) THEN
            ErrMsg = 'Could not allocate CollectionLevels!'
            CALL checkgc( RC, ErrMsg//ThisLoc )
            RETURN
         ENDIF
         CollectionLevels = UNDEFINED_STR
      ENDIF

      ! Allocate CollectionLevelInd
      IF ( .not. ALLOCATED( CollectionLevelInd ) ) THEN
         ALLOCATE( CollectionLevelInd( 2, CollectionCount ), STAT=RC )
         IF ( RC /= GC_SUCCESS ) THEN
            ErrMsg = 'Could not allocate CollectionLevelInt!'
            CALL checkgc( RC, ErrMsg//ThisLoc )
            RETURN
         ENDIF
         CollectionLevelInd = UNDEFINED_INT
      ENDIF

      ! Allocate CollectionHrRange
      IF ( .not. ALLOCATED( CollectionHrRange ) ) THEN
         ALLOCATE( CollectionHrRange( CollectionCount ), STAT=RC )
         IF ( RC /= GC_SUCCESS ) THEN
            ErrMsg = 'Could not allocate CollectionHrRange!'
            CALL checkgc( RC, ErrMsg//ThisLoc )
            RETURN
         ENDIF
         CollectionHrRange = UNDEFINED_STR
      ENDIF

   end subroutine get_collection_name

   subroutine checkgc(RC, message)
      implicit none
      integer :: RC
      character(len=*) :: message

      ! It is a good idea to check returned value for possible error
      if (RC .NE. GC_SUCCESS) then
         write(*,*) message
         return
      end if
   end subroutine checkgc

   !------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: History_AddItemToCollection
!
! !DESCRIPTION: Creates a HISTORY ITEM object for a given diagnostic quantity,
!  and then attaches it to a given diagnostic collection.  Given the name
!  of the diagnostic quantity, it will obtain metadata (and pointers to the
!  data array) via the appropriate state registry.
!\\
!\\
! !INTERFACE:
!
   SUBROUTINE History_AddItemToCollection( Input_Opt,                         &
      State_Chm,    State_Diag,          &
      State_Met,    Collection,          &
      CollectionId, ItemName,            &
      ItemCount,    SubsetInd,           &
      LevelInd,     RC                  )
!
! !USES:
!
      USE Charpak_Mod,           ONLY : To_UpperCase
      USE ErrCode_Mod
      USE HistContainer_Mod
      USE HistItem_Mod
      USE History_Util_Mod,      ONLY : UNDEFINED_INT
      USE Input_Opt_Mod,         ONLY : OptInput
      USE MetaHistContainer_Mod
      USE MetaHistItem_Mod
      USE Registry_Mod,          ONLY : Registry_Lookup
      USE State_Chm_Mod
      USE State_Diag_Mod
      USE State_Met_Mod
!
! !INPUT PARAMETERS:
!
! Required arguments
      TYPE(OptInput),      INTENT(IN)  :: Input_Opt      ! Input Options State
      TYPE(ChmState),      INTENT(IN)  :: State_Chm      ! Chemistry State
      TYPE(DgnState),      INTENT(IN)  :: State_Diag     ! Diagnostic State
      TYPE(MetState),      INTENT(IN)  :: State_Met      ! Meteorology State
      INTEGER,             INTENT(IN)  :: CollectionID   ! Collection ID number
      CHARACTER(LEN=255),  INTENT(IN)  :: ItemName       ! Name of HISTORY ITEM
      INTEGER,             INTENT(IN)  :: ItemCount      ! Index of HISTORY ITEM

! Optional arguments
      INTEGER,             OPTIONAL    :: SubsetInd(4)    ! X0,X1,Y0,Y1 indices
      INTEGER,             OPTIONAL    :: LevelInd(2)     ! Z0,Z1 indices
!
! !INPUT/OUTPUT PARAMETERS:
!
      TYPE(HistContainer), POINTER     :: Collection     ! Diagnostic Collection
!
! !OUTPUT PARAMETERS:
!
      INTEGER,             INTENT(OUT) :: RC             ! Success or failure?
!
! !REMARKS:
!  Private routine, called from History_Init.
!
! !REVISION HISTORY:
!  06 Jan 2015 - R. Yantosca - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
! Scalars
      LOGICAL                      :: OnLevelEdges
      INTEGER                      :: Source_KindVal
      INTEGER                      :: Output_KindVal
      INTEGER                      :: Rank
      INTEGER                      :: NX, X0, X1
      INTEGER                      :: NY, Y0, Y1
      INTEGER                      :: NZ, Z0, Z1

! Arrays
      INTEGER                      :: Dimensions(3)
      INTEGER                      :: ItemDims(3)
      INTEGER                      :: Subset_X(2)
      INTEGER                      :: Subset_Y(2)
      INTEGER                      :: Subset_Z(2)

! Strings
      CHARACTER(LEN=4  )           :: StateMetUC
      CHARACTER(LEN=5  )           :: StateChmUC
      CHARACTER(LEN=255)           :: ItemNameUC
      CHARACTER(LEN=255)           :: Description
      CHARACTER(LEN=255)           :: ThisLoc
      CHARACTER(LEN=255)           :: Units
      CHARACTER(LEN=512)           :: ErrMsg

! Objects
      TYPE(HistItem),      POINTER :: Item

! Pointer arrays
      REAL(fp),            POINTER :: Ptr0d
      REAL(f8),            POINTER :: Ptr0d_8
      REAL(f4),            POINTER :: Ptr0d_4
      INTEGER,             POINTER :: Ptr0d_I
      REAL(fp),            POINTER :: Ptr1d  (:    )
      REAL(f8),            POINTER :: Ptr1d_8(:    )
      REAL(f4),            POINTER :: Ptr1d_4(:    )
      INTEGER,             POINTER :: Ptr1d_I(:    )
      REAL(fp),            POINTER :: Ptr2d  (:,:  )
      REAL(f8),            POINTER :: Ptr2d_8(:,:  )
      REAL(f4),            POINTER :: Ptr2d_4(:,:  )
      INTEGER,             POINTER :: Ptr2d_I(:,:  )
      REAL(fp),            POINTER :: Ptr3d  (:,:,:)
      REAL(f8),            POINTER :: Ptr3d_8(:,:,:)
      REAL(f4),            POINTER :: Ptr3d_4(:,:,:)
      INTEGER,             POINTER :: Ptr3d_I(:,:,:)

!=======================================================================
! Initialize
!=======================================================================
      RC             =  GC_SUCCESS
      Description    =  ''
      Dimensions     =  0
      Source_KindVal =  0
      Output_KindVal =  0
      Rank           =  0
      Units          =  ''
      ErrMsg         =  ''
      ThisLoc        =  &
         ' -> History_AddItemToCollection (in History/history_mod.F90)'
      ItemNameUC     = To_UpperCase( ItemName )
      StateMetUC     = State_Met%State // '_'   ! State_Met%State is uppercase
      StateChmUC     = State_Chm%State // '_'   ! State_Chm%State is uppercase

! Free pointers
      Ptr0d    => NULL()
      Ptr0d_8  => NULL()
      Ptr0d_4  => NULL()
      Ptr0d_I  => NULL()
      Ptr1d    => NULL()
      Ptr1d_8  => NULL()
      Ptr1d_4  => NULL()
      Ptr1d_I  => NULL()
      Ptr2d    => NULL()
      Ptr2d_8  => NULL()
      Ptr2d_4  => NULL()
      Ptr2d_I  => NULL()
      Ptr3d    => NULL()
      Ptr3d_8  => NULL()
      Ptr3d_4  => NULL()
      Ptr3d_I  => NULL()

!=======================================================================
! For each HISTORY ITEM, find the matching entry in the relevant
! registry (in State_Chm, State_Diag, State_Met) and get a pointer
! to the data source
!=======================================================================
      IF ( ItemNameUC(1:5) == StateChmUC ) THEN

!--------------------------------------------------------------------
! Chemistry State
!--------------------------------------------------------------------
         CALL Registry_Lookup( am_I_Root      = Input_Opt%amIRoot,             &
            Registry       = State_Chm%Registry,            &
            RegDict        = State_Chm%RegDict,             &
            State          = State_Chm%State,               &
            Variable       = ItemName,                      &
            Description    = Description,                   &
            Dimensions     = Dimensions,                    &
            Source_KindVal = Source_KindVal,                &
            Output_KindVal = Output_KindVal,                &
            Rank           = Rank,                          &
            Units          = Units,                         &
            OnLevelEdges   = OnLevelEdges,                  &
            Ptr0d_8        = Ptr0d_8,                       &
            Ptr1d_8        = Ptr1d_8,                       &
            Ptr2d_8        = Ptr2d_8,                       &
            Ptr3d_8        = Ptr3d_8,                       &
            Ptr0d_4        = Ptr0d_4,                       &
            Ptr1d_4        = Ptr1d_4,                       &
            Ptr2d_4        = Ptr2d_4,                       &
            Ptr3d_4        = Ptr3d_4,                       &
            Ptr0d_I        = Ptr0d_I,                       &
            Ptr1d_I        = Ptr1d_I,                       &
            Ptr2d_I        = Ptr2d_I,                       &
            Ptr3d_I        = Ptr3d_I,                       &
            RC             = RC                            )

! Trap potential not found error
         IF ( RC /= GC_SUCCESS ) THEN
            ErrMsg = 'Could not locate ' // TRIM( ItemName )  // &
               ' chemistry state registry.'
            CALL GC_Error( ErrMsg, RC, ThisLoc )
            RETURN
         ENDIF

      ELSE IF ( ItemNameUC(1:4) == StateMetUC ) THEN

!--------------------------------------------------------------------
! Meteorology State
!--------------------------------------------------------------------
         CALL Registry_Lookup( am_I_Root      = Input_Opt%amIRoot,             &
            Registry       = State_Met%Registry,            &
            RegDict        = State_Met%RegDict,             &
            State          = State_Met%State,               &
            Variable       = ItemName,                      &
            Description    = Description,                   &
            Dimensions     = Dimensions,                    &
            Source_KindVal = Source_KindVal,                &
            Output_KindVal = Output_KindVal,                &
            Rank           = Rank,                          &
            Units          = Units,                         &
            OnLevelEdges   = OnLevelEdges,                  &
            Ptr0d_8        = Ptr0d_8,                       &
            Ptr1d_8        = Ptr1d_8,                       &
            Ptr2d_8        = Ptr2d_8,                       &
            Ptr3d_8        = Ptr3d_8,                       &
            Ptr0d_4        = Ptr0d_4,                       &
            Ptr1d_4        = Ptr1d_4,                       &
            Ptr2d_4        = Ptr2d_4,                       &
            Ptr3d_4        = Ptr3d_4,                       &
            Ptr0d_I        = Ptr0d_I,                       &
            Ptr1d_I        = Ptr1d_I,                       &
            Ptr2d_I        = Ptr2d_I,                       &
            Ptr3d_I        = Ptr3d_I,                       &
            RC             = RC                            )

! Trap potential not found error
         IF ( RC /= GC_SUCCESS ) THEN
            ErrMsg = 'Could not locate ' // TRIM( ItemName )  // &
               ' meteorology state registry.'
            CALL GC_Error( ErrMsg, RC, ThisLoc )
            RETURN
         ENDIF

      ELSE

!--------------------------------------------------------------------
! Diagnostic State
!--------------------------------------------------------------------
         CALL Registry_Lookup( am_I_Root      = Input_Opt%amIRoot,             &
            Registry       = State_Diag%Registry,           &
            RegDict        = State_Diag%RegDict,            &
            State          = State_Diag%State,              &
            Variable       = ItemName,                      &
            Description    = Description,                   &
            Dimensions     = Dimensions,                    &
            Source_KindVal = Source_KindVal,                &
            Output_KindVal = Output_KindVal,                &
            Rank           = Rank,                          &
            Units          = Units,                         &
            OnLevelEdges   = OnLevelEdges,                  &
            Ptr0d_8        = Ptr0d_8,                       &
            Ptr1d_8        = Ptr1d_8,                       &
            Ptr2d_8        = Ptr2d_8,                       &
            Ptr3d_8        = Ptr3d_8,                       &
            Ptr0d_4        = Ptr0d_4,                       &
            Ptr1d_4        = Ptr1d_4,                       &
            Ptr2d_4        = Ptr2d_4,                       &
            Ptr3d_4        = Ptr3d_4,                       &
            Ptr0d_I        = Ptr0d_I,                       &
            Ptr1d_I        = Ptr1d_I,                       &
            Ptr2d_I        = Ptr2d_I,                       &
            Ptr3d_I        = Ptr3d_I,                       &
            RC             = RC                            )

! Trap potential not found error
         IF ( RC /= GC_SUCCESS ) THEN
            ErrMsg = 'Could not locate ' // TRIM( ItemName )  // &
               ' diagnostics state registry.'
            CALL GC_Error( ErrMsg, RC, ThisLoc )
            RETURN
         ENDIF
      ENDIF

!=======================================================================
! If the optional SUBSETS and/or LEVELS arguments are passed, then use
! these to size the data arrays.  Otherwise assume the data arrays will
! be the same size as the pointer to the data source (as will be true
! in most cases.)
!=======================================================================

!-------------------------
! Default values
!-------------------------

! By default, use the size of the data array to define the
! X0, Y0, X1, Y1, Z0, and Z1 indices for the subset region.
      X0 = 1
      X1 = MAX( Dimensions(1), 1 )
      Y0 = 1
      Y1 = MAX( Dimensions(2), 1 )
      Z0 = 1
      Z1 = MAX( Dimensions(3), 1 )

!-------------------------
! Horizontal subsetting
!-------------------------

! If SubsetInd has valid values, use them to redefine X0, Y0, X1, and Y1.
      IF ( PRESENT( SubsetInd ) ) THEN
         IF ( SubsetInd(1) /= UNDEFINED_INT ) X0 = SubsetInd(1)
         IF ( SubsetInd(2) /= UNDEFINED_INT ) X1 = SubsetInd(2)
         IF ( SubsetInd(3) /= UNDEFINED_INT ) Y0 = SubsetInd(3)
         IF ( SubsetInd(4) /= UNDEFINED_INT ) Y1 = SubsetInd(4)
      ENDIF

!-------------------------
! Vertical subsetting
!-------------------------

! If LevelInd has valid values, use them to redefine Z0 and Z1.
      IF ( PRESENT( LevelInd ) ) THEN
         IF ( LevelInd(1) /= UNDEFINED_INT ) Z0 = LevelInd(1)
         IF ( LevelInd(2) /= UNDEFINED_INT ) Z1 = LevelInd(2)
      ENDIF

! Error check X-dimension indices
      IF ( X1 < X0 ) THEN
         WRITE( ErrMsg, 100 ) X0, X1, TRIM( Collection%Name )
100      FORMAT(  'Invalid X-dimension indices: ', 2i6, ' for collection', a )
         CALL GC_Error( ErrMsg, RC, ThisLoc )
         RETURN
      ENDIF

! Error-check Y-dimension indices
      IF ( Y1 < Y0 ) THEN
         WRITE( ErrMsg, 110 ) Y0, Y1, TRIM( Collection%Name )
110      FORMAT(  'Invalid Y-dimension indices: ', 2i6, ' for collection', a )
         CALL GC_Error( ErrMsg, RC, ThisLoc )
         RETURN
      ENDIF

! Error-check Z-dimension indices
      IF ( Z1 < Z0 ) THEN
         WRITE( ErrMsg, 120 ) Z0, Z1, TRIM( Collection%Name )
120      FORMAT(  'Invalid Y-dimension indices: ', 2i6, ' for collection', a )
         CALL GC_Error( ErrMsg, RC, ThisLoc )
         RETURN
      ENDIF

! Compute dimension extent
      NX = X1 - X0 + 1
      NY = Y1 - Y0 + 1
      NZ = Z1 - Z0 + 1

! Indices for subsetting the data
      Subset_X = (/ X0, X1 /)
      Subset_Y = (/ Y0, Y1 /)
      Subset_Z = (/ Z0, Z1 /)

! Save the subsets to the collection
      Collection%X0 = X0
      Collection%X1 = X1
      Collection%Y0 = Y0
      Collection%Y1 = Y1
      Collection%Z0 = Z0
!NOTE: Z1 is not needed, we compute that later!

!=======================================================================
! Now that we have obtained information (and pointers to the data)
! corresponding to the given diagnostic quantity, use that to create
! a HISTORY ITEM object for that diagnostic quantity.
!=======================================================================
      CALL HistItem_Create( Input_Opt      = Input_Opt,                        &
         Item           = Item,                             &
         Id             = ItemCount,                        &
         ContainerId    = CollectionId,                     &
         Name           = ItemName,                         &
         LongName       = Description,                      &
         Units          = Units,                            &
         OnLevelEdges   = OnLevelEdges,                     &
         SpaceDim       = Rank,                             &
         Operation      = Collection%Operation,             &
         Subset_X       = Subset_X,                         &
         Subset_Y       = Subset_Y,                         &
         Subset_Z       = Subset_Z,                         &
         Source_KindVal = Source_KindVal,                   &
         Output_KindVal = Output_KindVal,                   &
         Source_0d_8    = Ptr0d_8,                          &
         Source_1d_8    = Ptr1d_8,                          &
         Source_1d_4    = Ptr1d_4,                          &
         Source_1d_I    = Ptr1d_I,                          &
         Source_2d_8    = Ptr2d_8,                          &
         Source_2d_4    = Ptr2d_4,                          &
         Source_2d_I    = Ptr2d_I,                          &
         Source_3d_8    = Ptr3d_8,                          &
         Source_3d_4    = Ptr3d_4,                          &
         Source_3d_I    = Ptr3d_I,                          &
         Dimensions     = ItemDims,                         &
         RC             = RC                               )

! Trap potential error
      IF ( RC /= GC_SUCCESS ) THEN
         ErrMsg = 'Could not create Item: "' // TRIM( ItemName ) // '"!'
         CALL GC_Error( ErrMsg, RC, ThisLoc )
         RETURN
      ENDIF

!=======================================================================
! Attach this HISTORY ITEM to the METAHISTORY ITEM (aka list of HISTORY
! ITEMS) belonging to the HISTORY CONTAINER object for the given
! diagnostic collection.
!
! In other words, we are adding this diagnostic quantity to the list
! of diagnostic quantities that belong to this diagnostic collection.
! These quantities will be written to the netCDF file described by
! the collection, with the specified update frequency.
!=======================================================================
      CALL MetaHistItem_AddNew( Input_Opt = Input_Opt,                         &
         Node      = Collection%HistItems,              &
         Item      = Item,                              &
         RC        = RC                                )

! Trap potential error
      IF ( RC /= GC_SUCCESS ) THEN
         ErrMsg = 'Could not add Item "' //                                    &
            TRIM( ItemName )       //  '" to '   //                      &
            TRIM( CollectionName(CollectionId) ) // '%HistItems!'
         CALL GC_Error( ErrMsg, RC, ThisLoc )
         RETURN
      ENDIF

!=======================================================================
! Define the dimensions of the collection from the HISTORY ITEMS,
! and also whether vertical data is on the level centers or edges
!=======================================================================

! Define the X dimension of the collection from the
! first HISTORY ITEM that has an X dimension
      IF ( Collection%NX == UNDEFINED_INT ) THEN
         SELECT CASE( Item%DimNames )
          CASE( 'xyz', 'xz', 'xy', 'x' )
            Collection%NX = ItemDims(1)
          CASE DEFAULT
! Nothing
         END SELECT
      ENDIF

! Define the Y dimension of the collection from the
! first HISTORY ITEM that has a Y dimension
      IF ( Collection%NY == UNDEFINED_INT ) THEN
         SELECT CASE( Item%DimNames )
          CASE( 'xyz', 'xy' )
            Collection%NY = ItemDims(2)
          CASE( 'yz', 'y' )
            Collection%NY = ItemDims(1)
          CASE DEFAULT
! Nothing
         END SELECT
      ENDIF

! Define the Z dimension of the collection from the first HISTORY ITEM
! that has a Z dimension. Also define whether the collection will
! contain data that is centered or edged on vertical levels.
      IF ( Collection%NZ == UNDEFINED_INT ) THEN
         SELECT CASE( Item%DimNames )
          CASE( 'xyz' )
            Collection%NZ = ItemDims(3)
          CASE( 'xz', 'yz' )
            Collection%NZ = ItemDims(2)
          CASE( 'z' )
            Collection%NZ = ItemDims(1)
          CASE DEFAULT
! Nothing
         END SELECT

         Collection%OnLevelEdges = Item%OnLevelEdges
      ENDIF

! NB is always two (lon0, lon1) or (lat0, lat1)
      IF ( Collection%NB == UNDEFINED_INT ) THEN
         Collection%NB = 2
      ENDIF

!=======================================================================
! Make sure that all the HISTORY ITEMS in this collection are
! placed on the level centers or edges, but not both.  The netCDF
! COARDS/CF conventions do not allow for data on more than one
! vertical dimension per file.
!=======================================================================
      IF ( Item%SpaceDim == 3 ) THEN
         IF ( Collection%OnLevelEdges .neqv. Item%OnLevelEdges ) THEN
            ErrMsg = TRIM( Item%Name )                                      // &
               ' has the wrong vertical alignment for collection: "'  // &
               TRIM( Collection%Name )  // '".  Please check your '   // &
               'HISTORY.rc file to make sure that this collection '   // &
               'only contains 3-D diagnostics with the same vertical '// &
               'alignment.  You cannot add diagnostics that are '     // &
               'defined on level centers and diagnostics that are '   // &
               'defined on level edges in the same collection, as '   // &
               'per netCDF conventions.'
            CALL GC_Error( ErrMsg, RC, ThisLoc )
            RETURN
         ENDIF
      ENDIF

!=======================================================================
! Cleanup and quit
!=======================================================================

! Increment number of HISTORY ITEMS in this collection
      Collection%nHistItems = Collection%nHistItems + 1

! Free Item (we have added it to the container and don't need it anymore)
      DEALLOCATE( Item )
      Item => NULL()

! Free pointers
      Ptr0d   => NULL()
      Ptr0d_8 => NULL()
      Ptr0d_4 => NULL()
      Ptr0d_I => NULL()
      Ptr1d   => NULL()
      Ptr1d_8 => NULL()
      Ptr1d_4 => NULL()
      Ptr1d_I => NULL()
      Ptr2d   => NULL()
      Ptr2d_8 => NULL()
      Ptr2d_4 => NULL()
      Ptr2d_I => NULL()
      Ptr3d   => NULL()
      Ptr3d_8 => NULL()
      Ptr3d_4 => NULL()
      Ptr3d_I => NULL()

   END SUBROUTINE History_AddItemToCollection

   SUBROUTINE Get_Number_Of_Levels( Container, nLev, nILev )
      !
      ! !USES:
      !
          USE HistContainer_Mod, ONLY : HistContainer
      !
      ! !INPUT PARAMETERS:
      !
          TYPE(HistContainer), POINTER     :: Container ! Diagnostic collection obj
      !
      ! !OUTPUT PARAMETERS:
      !
          INTEGER,             INTENT(OUT) :: nLev      ! Number of levels
          INTEGER,             INTENT(OUT) :: nIlev     ! Number of level interfaces
      !
      ! !REVISION HISTORY:
      !  05 Jun 2019 - R. Yantosca - Initial version
      !  See https://github.com/geoschem/geos-chem for complete history
      !EOP
      !------------------------------------------------------------------------------
      !BOC
      
          !========================================================================
          ! Pick the dimensions of the lev and ilev variables properly
          ! so that we can use that for writing to then netCDF files.
          !
          ! If the vertical dimension (Container%NZ) is undefined, then
          ! this indicates that there is only 2-D data in the collection.
          ! Thus, there will be 1 level (the surface) and 2 level edges.
          !========================================================================
          IF ( Container%OnLevelEdges ) THEN
             nILev = MAX( Container%NZ, 2 )
             nLev  = nILev - 1
          ELSE
             nLev  = MAX( Container%NZ, 1 )
             nILev = nLev  + 1
          ENDIF
      
        END SUBROUTINE Get_Number_Of_Levels


end module WRFGC_History_Mod
