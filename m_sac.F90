!! ----------------------------------------------------------------------------------------------------------------------------- !!
!>
!! Read/Write SAC-formatted seismograms: Fortran 2003 Version
!<
!! ----------------------------------------------------------------------------------------------------------------------------- !!
module m_sac

  !! -- Dependency
  use m_std
  use m_endian
  use m_daytim

  !! -- Declarations
  implicit none
  private

  !! -- Public Procedures
  public :: sac__hdr     ! sac data type
  public :: sac__read    ! read  sac datafile
  public :: sac__write   ! write sac datafile
  public :: sac__init    ! initialize sac data type
  public :: sac__rhdr    ! read header
  public :: sac__whdr    ! read header


  !! --------------------------------------------------------------------------------------------------------------------------- !!
  !>
  !! Sac Header Type Definition
  !<
  !! --
  type sac__hdr
     
     !!               var name          description                   record#
     real(DP)      :: delta           ! sampling interval             (001)
     real(DP)      :: depmin          ! minimum value                 (002)
     real(DP)      :: depmax          ! maximum value                 (003)
     real(DP)      :: b               ! begenning independent value   (006)
     real(DP)      :: e               ! ending independent value      (007)
     real(DP)      :: o               ! event origin time             (008)
     real(DP)      :: a               ! first arrival time            (009)
     real(DP)      :: t0              ! time picks                    (011)
     real(DP)      :: t1              ! time picks                    (012)
     real(DP)      :: t2              ! time picks                    (013)
     real(DP)      :: t3              ! time picks                    (014)
     real(DP)      :: t4              ! time picks                    (015)
     real(DP)      :: t5              ! time picks                    (016)
     real(DP)      :: t6              ! time picks                    (017)
     real(DP)      :: t7              ! time picks                    (018)
     real(DP)      :: t8              ! time picks                    (019)
     real(DP)      :: t9              ! time picks                    (020)
     real(DP)      :: stla            ! station latitude              (032)
     real(DP)      :: stlo            ! station longitude             (033)
     real(DP)      :: stel            ! station elevation (m)         (034)
     real(DP)      :: stdp            ! station depth (m)             (035)
     real(DP)      :: evla            ! event latitude                (036)
     real(DP)      :: evlo            ! event longitude               (037)
     real(DP)      :: evel            ! event elevation (m)           (038)
     real(DP)      :: evdp            ! event depth (m)               (039)
     real(DP)      :: mag             ! event magnitude               (040)
     real(DP)      :: user0           ! user header                   (041)
     real(DP)      :: user1           ! user header                   (042)
     real(DP)      :: user2           ! user header                   (043)
     real(DP)      :: user3           ! user header                   (044)
     real(DP)      :: user4           ! user header                   (045)
     real(DP)      :: user5           ! user header                   (046)
     real(DP)      :: user6           ! user header                   (047)
     real(DP)      :: user7           ! user header                   (048)
     real(DP)      :: user8           ! user header                   (049)
     real(DP)      :: user9           ! user header                   (050)
     real(DP)      :: dist            ! distance (km)                 (051)
     real(DP)      :: az              ! azimuth (deg)                 (052)
     real(DP)      :: baz             ! back azimuth (deg)            (053)
     real(DP)      :: gcarc           ! angular distance (deg)        (054)
     real(DP)      :: depmen          ! mean value                    (057)
     real(DP)      :: cmpaz           ! component azimuth             (058)
     real(DP)      :: cmpinc          ! component incident angle      (059)     
     integer       :: nzyear          ! reference time, year          (071)
     integer       :: nzjday          ! reference time, julian day    (072)
     integer       :: nzhour          ! reference time, hour          (073)
     integer       :: nzmin           ! reference time, minute        (074)
     integer       :: nzsec           ! reference time, second        (075)
     integer       :: nzmsec          ! reference time, millisecond   (076)
     integer       :: nvhdr           ! header version                (077)
     integer       :: npts            ! number of data points         (080)
     integer       :: iftype          ! type of file                  (086)
     integer       :: idep            ! type of dependent var.        (087)
     integer       :: ievtyp          ! event type                    (093)
     logical       :: leven           ! is evenly spaced file         (106)
     logical       :: lpspol          ! is positive polarity          (107)
     logical       :: lovrok          ! is overwrite ok?              (108)
     logical       :: lcalda          ! is calc distance azimuth      (109)
     character(8)  :: kstnm           ! station name                  (111)
     character(16) :: kevnm           ! event name                    (113)
     character(8)  :: ka              ! time pick name                (123)
     character(8)  :: kt0             ! time pick name                (123)
     character(8)  :: kt1             ! time pick name                (125)
     character(8)  :: kt2             ! time pick name                (127)
     character(8)  :: kt3             ! time pick name                (129)
     character(8)  :: kt4             ! time pick name                (131)
     character(8)  :: kt5             ! time pick name                (133)
     character(8)  :: kt6             ! time pick name                (135)
     character(8)  :: kt7             ! time pick name                (137)
     character(8)  :: kt8             ! time pick name                (139)
     character(8)  :: kt9             ! time pick name                (141)
     character(8)  :: kuser0          ! user area                     (145)
     character(8)  :: kuser1          ! user area                     (147)
     character(8)  :: kuser2          ! user area                     (149)
     character(8)  :: kcmpnm          ! component name                (151)
     character(8)  :: knetwk          ! network name                  (153)
     character(8)  :: kdatrd          ! date data onto comp.          (155)
     character(8)  :: kinst           ! instrument                    (157)
     
     !! associated information from sac header
     integer :: nzmonth ! month of begin time from nzjday
     integer :: nzday   ! day   of begin time from nzjday
     
     integer :: tim     ! absolute begin time from 1970/1/1 0:0:0 in second
     logical :: is_same_endian
     
  end type sac__hdr
  !! --------------------------------------------------------------------------------------------------------------------------- !!

  
  !! --------------------------------------------------------------------------------------------------------------------------- !!
  !>
  !! Read SAC file
  !!
  !!
  !! @par Usage
  !! call sac__read( char filename, sac__hdr sh, real(*) dat(:) )
  !!
  !! Data array can be double or single precisions. 
  !! If the data array dat(:) are not allocated, the routine automatically 
  !! allocate memory according to the file length defined in the sac header
  !! If the routine is called without data array, it reads only header part.
  !<
  !! --
  interface sac__read
     
     module procedure rsac_d, rsac_s
     
  end interface
  !! --------------------------------------------------------------------------------------------------------------------------- !!
  
  !! --------------------------------------------------------------------------------------------------------------------------- !!
  !>
  !! Write SAC file
  !!
  !! @par Usage
  !! call sac__write( char filename, type__header, real data(:), logical sw )
  !! data can be single or double precisions
  !! if sw = true, the existing file is automatically replaced. 
  !<
  !! --
  interface sac__write
     
     module procedure wsac_d, wsac_s
     
  end interface
  !! --------------------------------------------------------------------------------------------------------------------------- !!
  
contains
  
  !! --------------------------------------------------------------------------------------------------------------------------- !!
  !>
  !! Read the sac file fn_sac
  !! If optional argument header_only = .true., rsac will not read data part.
  !<
  !! --
  subroutine rsac_d ( fn_sac, ss, dat )

    !! --Arguments
    character(*),          intent(in)    :: fn_sac  !< sac filename
    type(sac__hdr),        intent(out)   :: ss      !< header info
    real(DP), allocatable, intent(inout) :: dat(:)  !< waveform data
    real(SP), allocatable                :: fdat(:)

    !! ----
    
    call rsac_s( fn_sac, ss, fdat )

    if( .not. allocated( dat ) )  allocate( dat(1:ss%npts) )
    dat = dble(fdat)
    
  end subroutine rsac_d
  !! --------------------------------------------------------------------------------------------------------------------------- !!
  
  !! --------------------------------------------------------------------------------------------------------------------------- !!
  !>
  !! Read the sac file fn_sac
  !! If optional argument header_only = .true., rsac will not read data part.
  !<
  !! --
  subroutine rsac_s( fn_sac, ss, dat )

    !! --Arguments
    character(*),                    intent(in)    :: fn_sac !< sac filename
    type(sac__hdr),                  intent(out)   :: ss     !< header info
    real(SP), allocatable, optional, intent(inout) :: dat(:) !< waveform data

    integer :: io
    integer :: i
    integer :: nmax

    !! ----
    
    call std__getio( io )
    open( io, file=fn_sac, action='read', access='stream', form='unformatted' )
    
    call sac__rhdr(io,ss)
    
    ! store month, day, localtime info
    call daytim__jul2md   ( ss%nzjday, ss%nzyear, ss%nzmonth, ss%nzday )
    call daytim__timelocal( ss%nzyear, ss%nzmonth, ss%nzday,  ss%nzhour, ss%nzmin, ss%nzsec, ss%tim )
    
    if( present( dat ) ) then       
       
       if( .not. allocated( dat ) ) then
          allocate( dat(1:ss%npts))
          nmax = ss%npts
       else
          nmax = max( ss%npts, size(dat) )
       end if
       
       if( ss%npts > nmax ) then
          write(STDERR,*) 'rsac: data array does not have enough size' 
       end if
       
       read(io) dat
       if( .not. ss%is_same_endian ) then
          do i=1, ss%npts
             call endian__change( dat(i) )
          end do
       end if
       
    end if
    
    close( io )
    
  end subroutine rsac_s
  !! --------------------------------------------------------------------------------------------------------------------------- !!
  
  !! --------------------------------------------------------------------------------------------------------------------------- !!
  !>
  !! Write SAC file
  !<
  !! --
  subroutine wsac_d( fn_sac, ss, dat, overwrite )

    !! -- Arguments
    character(*),   intent(in)           :: fn_sac
    type(sac__hdr), intent(in)           :: ss
    real(DP),       intent(in)           :: dat(1:ss%npts)
    logical,        intent(in), optional :: overwrite

    !! ----
    
    if( present( overwrite) ) then
       call wsac_s( fn_sac, ss, real(dat), overwrite) 
    else
       call wsac_s( fn_sac, ss, real(dat) )
    end if
    
  end subroutine wsac_d
  !! --------------------------------------------------------------------------------------------------------------------------- !!

  !! --------------------------------------------------------------------------------------------------------------------------- !!
  !>
  !! Write SAC file
  !<
  subroutine wsac_s( fn_sac, ss, dat, overwrite )

    !! -- Arguments
    character(*),   intent(in)           :: fn_sac
    type(sac__hdr), intent(in)           :: ss
    real(SP),       intent(in)           :: dat(1:ss%npts)
    logical,        intent(in), optional :: overwrite

    logical        :: isexist
    integer        :: io
    character(1)   :: yn
    !! ----
    
    !! overwrite check
    inquire( file = fn_sac, exist=isexist )
    if( isexist ) then
       if( present( overwrite) ) then
          if( .not. overwrite ) then
             write(STDERR,*) 'wsac: file '//trim(fn_sac)//' exists.' 
             write(STDERR,*) 'wsac: could not overwrite the file.'
             write(STDERR,*) 'wsac: return without success'
             write(STDERR,*)
             return
          end if
       else
          write(STDERR,*) 'wsac: file '//trim(fn_sac)//' exists.' 
          write(STDERR,*) 'wsac: Overwrite ? (y/n)' 
          read(STDIN,'(A)') yn
          if( yn /= 'y' .and. yn /='Y' ) then
             write(STDERR,*) 'wsac: could not overwrite the file.'
             write(STDERR,*) 'wsac: return without success'
             write(STDERR,*)
             return
          end if
       end if
    end if

    
    call std__getio( io )
    open( io, file=trim(fn_sac), action='write', access='stream', form='unformatted')

    call sac__whdr(io, ss) 

    write( io ) dat(1:ss%npts)
    close( io )
    
  end subroutine wsac_s
  !! --------------------------------------------------------------------------------------------------------------------------- !!

  !! --------------------------------------------------------------------------------------------------------------------------- !!
  !>
  !! Write SAC data header from pre-opened file io
  !! No endian conversion will be made. Always write in machine-endian.
  !<
  !! --
  subroutine sac__whdr(io, ss)

    !! -- Arguments
    integer,        intent(in) :: io
    type(sac__hdr), intent(in) :: ss

    real(SP)       :: fheader(70)
    integer        :: iheader(71:105)
    logical        :: lheader(106:110)
    character(4)   :: aheader(111:158)
    integer        :: i
    !! ----

    !! header initialize
    fheader(1:70) = -12345.0
    iheader(71:105) = -12345
    lheader(106:110) = .false. 
    do i=111, 157, 2
       aheader( i ) = '-123'
       aheader( i+1 ) = '45'
    end do


    ! Copy header data to temprary arrays
    fheader(  1) = ss % delta
    fheader(  2) = ss % depmin
    fheader(  3) = ss % depmax
    fheader(  6) = ss % b     
    fheader(  7) = ss % e     
    fheader(  8) = ss % o     
    fheader(  9) = ss % a     
    fheader( 11) = ss % t0    
    fheader( 12) = ss % t1    
    fheader( 13) = ss % t2    
    fheader( 14) = ss % t3    
    fheader( 15) = ss % t4    
    fheader( 16) = ss % t5    
    fheader( 17) = ss % t6    
    fheader( 18) = ss % t7    
    fheader( 19) = ss % t8    
    fheader( 20) = ss % t9    
    fheader( 32) = ss % stla  
    fheader( 33) = ss % stlo  
    fheader( 34) = ss % stel  
    fheader( 35) = ss % stdp  
    fheader( 36) = ss % evla  
    fheader( 37) = ss % evlo  
    fheader( 38) = ss % evel  
    fheader( 39) = ss % evdp  
    fheader( 40) = ss % mag   
    fheader( 41) = ss % user0 
    fheader( 42) = ss % user1 
    fheader( 43) = ss % user2 
    fheader( 44) = ss % user3 
    fheader( 45) = ss % user4 
    fheader( 46) = ss % user5 
    fheader( 47) = ss % user6 
    fheader( 48) = ss % user7 
    fheader( 49) = ss % user8 
    fheader( 50) = ss % user9 
    fheader( 51) = ss % dist  
    fheader( 52) = ss % az    
    fheader( 53) = ss % baz   
    fheader( 54) = ss % gcarc 
    fheader( 57) = ss % depmen
    fheader( 58) = ss % cmpaz 
    fheader( 59) = ss % cmpinc
    
    iheader( 71) = ss % nzyear
    iheader( 72) = ss % nzjday
    iheader( 73) = ss % nzhour
    iheader( 74) = ss % nzmin
    iheader( 75) = ss % nzsec
    iheader( 76) = ss % nzmsec
    iheader( 77) = ss % nvhdr
    iheader( 80) = ss % npts
    iheader( 86) = ss % iftype
    iheader( 87) = ss % idep
    iheader( 93) = ss % ievtyp
    
    lheader(106) = ss % leven
    lheader(107) = ss % lpspol
    lheader(108) = ss % lovrok
    lheader(109) = ss % lcalda
    
    aheader(111) = ss%kstnm(1:4);  aheader(112) = ss%kstnm(5:8)
    aheader(113) = ss%kevnm(1:4);  aheader(114) = ss%kevnm(5:8)
    aheader(115) = ss%kevnm(9:12); aheader(116) = ss%kevnm(13:16)
    aheader(123) = ss%kt0(1:4);    aheader(124) = ss%kt0(5:8)
    aheader(125) = ss%kt1(1:4);    aheader(126) = ss%kt1(5:8)
    aheader(127) = ss%kt2(1:4);    aheader(128) = ss%kt2(5:8)
    aheader(129) = ss%kt3(1:4);    aheader(130) = ss%kt3(5:8)
    aheader(131) = ss%kt4(1:4);    aheader(132) = ss%kt4(5:8)
    aheader(133) = ss%kt5(1:4);    aheader(134) = ss%kt5(5:8)
    aheader(135) = ss%kt6(1:4);    aheader(136) = ss%kt6(5:8)
    aheader(137) = ss%kt7(1:4);    aheader(138) = ss%kt7(5:8)
    aheader(139) = ss%kt8(1:4);    aheader(140) = ss%kt8(5:8)
    aheader(141) = ss%kt9(1:4);    aheader(142) = ss%kt9(5:8)
    aheader(145) = ss%kuser0(1:4); aheader(146) = ss%kuser0(5:8)
    aheader(147) = ss%kuser1(1:4); aheader(148) = ss%kuser1(5:8)
    aheader(149) = ss%kuser2(1:4); aheader(150) = ss%kuser2(5:8)
    aheader(151) = ss%kcmpnm(1:4); aheader(152) = ss%kcmpnm(5:8)
    aheader(153) = ss%knetwk(1:4); aheader(154) = ss%knetwk(5:8)
    aheader(155) = ss%kdatrd(1:4); aheader(156) = ss%kdatrd(5:8)
    aheader(157) = ss%kinst(1:4);  aheader(158) = ss%kinst(5:8)

    !! write
    write( io ) fheader(1:70)
    write( io ) iheader(71:105)
    write( io ) lheader(106:110)
    write( io ) aheader(111:158)
    
  end subroutine sac__whdr
  !! --------------------------------------------------------------------------------------------------------------------------- !!

  !! --------------------------------------------------------------------------------------------------------------------------- !!
  !>
  !! Initialize SAC header
  !<
  !! --
  subroutine sac__init( ss )

    !! -- Arguments
    type(sac__hdr), intent(inout) :: ss

    real(SP)     :: ferr = -12345.0_SP
    integer      :: ierr = -12345
    character(6) :: cerr = '-12345'
    !! ----
    
    ss%delta   = ferr
    ss%depmin  = ferr
    ss%depmax  = ferr
    ss%b       = ferr
    ss%e       = ferr
    ss%o       = ferr
    ss%a       = ferr
    ss%t0      = ferr
    ss%t1      = ferr
    ss%t2      = ferr
    ss%t3      = ferr
    ss%t4      = ferr
    ss%t5      = ferr
    ss%t6      = ferr
    ss%t7      = ferr
    ss%t8      = ferr
    ss%t9      = ferr
    ss%stla    = ferr
    ss%stlo    = ferr
    ss%stel    = ferr
    ss%stdp    = ferr
    ss%evla    = ferr
    ss%evlo    = ferr
    ss%evel    = ferr
    ss%evdp    = ferr
    ss%mag     = ferr
    ss%user0   = ferr
    ss%user1   = ferr
    ss%user2   = ferr
    ss%user3   = ferr
    ss%user4   = ferr
    ss%user5   = ferr
    ss%user6   = ferr
    ss%user7   = ferr
    ss%user8   = ferr
    ss%user9   = ferr
    ss%dist    = ferr
    ss%az      = ferr
    ss%baz     = ferr
    ss%gcarc   = ferr
    ss%depmen  = ferr
    ss%cmpaz   = ferr
    ss%cmpinc  = ferr
    ss%nzyear  = ierr
    ss%nzjday  = ierr
    ss%nzhour  = ierr
    ss%nzmin   = ierr
    ss%nzsec   = ierr
    ss%nzmsec  = ierr
    ss%nvhdr   = 6
    ss%npts    = ierr
    ss%iftype  = 01
    ss%ievtyp  = ierr
    ss%idep    = ierr
    ss%leven   = .true.
    ss%lpspol  = .false.
    ss%lovrok  = .true.
    ss%lcalda  = .true.
    ss%kstnm   = cerr
    ss%kcmpnm  = cerr
    ss%kevnm   = cerr
    ss%kt0     = cerr
    ss%kt1     = cerr
    ss%kt2     = cerr
    ss%kt3     = cerr
    ss%kt4     = cerr
    ss%kt5     = cerr
    ss%kt6     = cerr
    ss%kt7     = cerr
    ss%kt8     = cerr
    ss%kt9     = cerr
    ss%kuser0  = cerr
    ss%kuser1  = cerr
    ss%kuser2  = cerr
    ss%knetwk  = cerr
    ss%kdatrd  = cerr
    ss%kinst   = cerr
    
    ss%nzmonth = ierr
    ss%nzday   = ierr
    ss%tim     = ierr
    
  end subroutine sac__init
  !! --------------------------------------------------------------------------------------------------------------------------- !!

  !! --------------------------------------------------------------------------------------------------------------------------- !!
  !>
  !! Read SAC data header from pre-opened file io
  !<
  !! --
  subroutine sac__rhdr( io, ss )

    integer,        intent(in) :: io
    type(sac__hdr), intent(out) :: ss

    integer :: i
    real(SP)     :: fheader(1:70)
    integer      :: iheader(71:105)
    logical      :: lheader(106:110)
    integer      :: nvhdr
    character(4) :: aheader(111:158)

    !! ----
    
    ! endian check
    read(io,pos=77*4-3) nvhdr
    if( 1<= nvhdr .and. nvhdr <= 6 ) then
       ss%is_same_endian = .true.
    else
       ss%is_same_endian = .false.
    end if

    read(io,pos=1) fheader
    read(io)       iheader
    read(io)       lheader
    read(io)       aheader
    
    if( .not. ss%is_same_endian ) then
       do i=  1, 70
          call endian__change(fheader(i))
       end do
       do i= 71,105
          call endian__change(iheader(i))
       end do
       do i=106,110
          call endian__change(lheader(i))
       end do
    end if
    
    ss % delta  = dble( int( fheader(  1) * 1d7 + 0.5 ) ) / 1d7
    ss % depmin = dble( fheader(  2)  )
    ss % depmax = dble( fheader(  3)  ) 
    ss % b      = dble( fheader(  6)  )
    ss % e      = dble( fheader(  7)  )
    ss % o      = dble( fheader(  8)  )
    ss % a      = dble( fheader(  9)  )
    ss % t0     = dble( fheader( 11)  )
    ss % t1     = dble( fheader( 12)  )
    ss % t2     = dble( fheader( 13)  )
    ss % t3     = dble( fheader( 14)  )
    ss % t4     = dble( fheader( 15)  )
    ss % t5     = dble( fheader( 16)  )
    ss % t6     = dble( fheader( 17)  )
    ss % t7     = dble( fheader( 18)  )
    ss % t8     = dble( fheader( 19)  )
    ss % t9     = dble( fheader( 20)  )
    ss % stla   = dble( fheader( 32)  )
    ss % stlo   = dble( fheader( 33)  )
    ss % stel   = dble( fheader( 34)  )
    ss % stdp   = dble( fheader( 35)  )
    ss % evla   = dble( fheader( 36)  )
    ss % evlo   = dble( fheader( 37)  )
    ss % evel   = dble( fheader( 38)  )
    ss % evdp   = dble( fheader( 39)  )
    ss % mag    = dble( fheader( 40)  )
    ss % user0  = dble( fheader( 41)  )
    ss % user1  = dble( fheader( 42)  )
    ss % user2  = dble( fheader( 43)  )
    ss % user3  = dble( fheader( 44)  )
    ss % user4  = dble( fheader( 45)  )
    ss % user5  = dble( fheader( 46)  )
    ss % user6  = dble( fheader( 47)  )
    ss % user7  = dble( fheader( 48)  )
    ss % user8  = dble( fheader( 49)  )
    ss % user9  = dble( fheader( 50)  )
    ss % dist   = dble( fheader( 51)  )
    ss % az     = dble( fheader( 52)  )
    ss % baz    = dble( fheader( 53)  )
    ss % gcarc  = dble( fheader( 54)  )
    ss % depmen = dble( fheader( 57)  )
    ss % cmpaz  = dble( fheader( 58)  )
    ss % cmpinc = dble( fheader( 59)  )
    
    ss % nzyear = iheader( 71) 
    ss % nzjday = iheader( 72) 
    ss % nzhour = iheader( 73) 
    ss % nzmin  = iheader( 74) 
    ss % nzsec  = iheader( 75) 
    ss % nzmsec = iheader( 76) 
    ss % nvhdr  = iheader( 77) 
    ss % npts   = iheader( 80) 
    ss % iftype = iheader( 86) 
    ss % idep   = iheader( 87) 
    ss % ievtyp = iheader( 93) 
    
    ss % leven  = lheader(106)
    ss % lpspol = lheader(107)
    ss % lovrok = lheader(108)
    ss % lcalda = lheader(109)
    
    ss % kstnm  = aheader(111) // aheader(112)
    ss % kevnm  = aheader(113) // aheader(114) // aheader(115) // aheader(116)
    ss % kt0    = aheader(123) // aheader(124)
    ss % kt1    = aheader(125) // aheader(126)
    ss % kt2    = aheader(127) // aheader(128)
    ss % kt3    = aheader(129) // aheader(130)
    ss % kt4    = aheader(131) // aheader(132)
    ss % kt5    = aheader(133) // aheader(134)
    ss % kt6    = aheader(135) // aheader(136)
    ss % kt7    = aheader(137) // aheader(138)
    ss % kt8    = aheader(139) // aheader(140)
    ss % kt9    = aheader(141) // aheader(142)
    ss % kuser0 = aheader(145) // aheader(146)
    ss % kuser1 = aheader(147) // aheader(148)
    ss % kuser2 = aheader(149) // aheader(150)
    ss % kcmpnm = aheader(151) // aheader(152)
    ss % knetwk = aheader(153) // aheader(154)
    ss % kdatrd = aheader(155) // aheader(156)
    ss % kinst  = aheader(157) // aheader(158)
    
    call char_zeropad( ss%kstnm  )
    call char_zeropad( ss%kevnm  )
    call char_zeropad( ss%kt0    )
    call char_zeropad( ss%kt1    )
    call char_zeropad( ss%kt2    )
    call char_zeropad( ss%kt3    )
    call char_zeropad( ss%kt4    )
    call char_zeropad( ss%kt5    )
    call char_zeropad( ss%kt6    )
    call char_zeropad( ss%kt7    )
    call char_zeropad( ss%kt8    )
    call char_zeropad( ss%kt9    )
    call char_zeropad( ss%kuser0 )
    call char_zeropad( ss%kuser1 )
    call char_zeropad( ss%kuser2 )
    call char_zeropad( ss%kcmpnm )
    call char_zeropad( ss%knetwk )
    
    
  end subroutine sac__rhdr
  !! --------------------------------------------------------------------------------------------------------------------------- !!
  
  !! --------------------------------------------------------------------------------------------------------------------------- !!
  !>
  !! Substitute blanks in character ch by padding '0'
  !<
  !! --
  subroutine char_zeropad( ch )

    !! -- Arguments
    character(*), intent(inout) :: ch

    integer :: i
    !! ----
    
    do i=1, len(ch)
       if( ichar(ch(i:i))== 0 )then
          ch(i:i) = ' ' 
       end if
    end do
    
  end subroutine char_zeropad
  !! --------------------------------------------------------------------------------------------------------------------------- !!
  
end module m_sac
!! ----------------------------------------------------------------------------------------------------------------------------- !!
