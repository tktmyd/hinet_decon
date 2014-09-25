!! ----------------------------------------------------------------------------------------------------------------------------- !!
!>
!! Fortran2003 codes to deconvolve Hi-net velocity record by its seismometer response by using inverse filtering technique. 
!!
!! @see
!!    Maeda, T., K. Obara, T. Furumura, and T. Saito,
!!       Interference of long-period seismic wavefield observed by dense Hi-net array in Japan,
!!       J. Geophys. Res., 116, B10303, doi:10.1029/2011JB008464, 2011.
!!
!! @par Usage
!!
!!   hinet_decon.x sacfile [-o <outfile>] [-f0 <freq>] [-h0 <damp>] [-f1 <freq>]
!!                         [-h1 <damp>] [-decon] [-int] [-fint <freq>] [-hint <freq>]
!!
!! @par Options <default>
!!
!!   -  '-o  <outfile>' ( indicate output sac filename <sacfile.out>          )
!!   -  '-f0 <freq>'    ( eigen frequency of input signal <1.0>               )
!!   -  '-h0 <damp>'    ( damping constant of input signal <0.7>              )
!!   -  '-f1 <freq>'    ( eigen frequency of output signal <0.00833333>       )
!!   -  '-h1 <damp>'    ( damping constant of output signal <0.707>           )
!!   -  '-decon'        ( deconvolve without simulation seismometer           )
!!   -  '-int'          ( output displacement record by numerical integration )
!!   -  '-fint <freq>'  ( corner frequency of low-cut filter for integration  )
!!   -  '-hint <damp>'  ( damping constant of low-cut filter for integration  )
!!  
!! @par Hint
!!  ./hinet_decon (input.sac) -o (output.sac)
!!      This default option will convert Hi-net type velocity seismometer response to that of (approximated) STS-2. 
!!      By adding "-int" option, deconvolved displacement waveform will be obtained.
!<
!! --
program hinet_decon

  !! -- Dependencices
  use m_std
  use m_sac
  use m_tdecon
  use m_getopt
  use m_system
  use m_rtrend

  !! -- Declarations
  implicit none

  !! -- Parameters
  real(DP), parameter :: D_F0   = 1.0_DP             !< Default input eigenfreq
  real(DP), parameter :: D_H0   = 0.7_DP             !< Default input damping constant
  real(DP), parameter :: D_F1   = 1.0_DP / 120.0_DP  !< Default output eigenfreq
  real(DP), parameter :: D_H1   = 0.707_DP           !< Default output damping constant
  real(DP), parameter :: D_FINT = 1.0_DP / 80.0_DP   !< Default cut-off freq for integration
  real(DP), parameter :: D_HINT = 0.6321_DP          !< Default damping constant for integration
  integer,  parameter :: NC     = 256                !< character length

  !! --

  !=Variables
  character(NC)         :: fn_in
  character(NC)         :: fn_sac
  character(NC)         :: fn_out
  type(sac__hdr)        :: sh
  logical               :: is_exist
  logical               :: is_decon
  logical               :: is_integ
  logical               :: is_list
  real(DP)              :: f0
  real(DP)              :: f1
  real(DP)              :: h0
  real(DP)              :: h1
  real(DP)              :: hint
  real(DP)              :: fint
  real(DP)              :: a(0:2)
  real(DP)              :: b(0:2)
  real(DP)              :: ai(0:2)
  real(DP)              :: bi(0:2)
  real(DP), allocatable :: wav(:)
  integer               :: io
  integer               :: ierr
  integer               :: i
  integer               :: nfile
  !! ----

  if( system__iargc() == 0 ) call usage_exit()

  call system__getarg(1, fn_in)

  !! Option processing
  call getopt( 'o'    , is_exist, fn_out, trim(fn_in)//'.out' )
  call getopt( 'f0'   , is_exist, f0,     D_F0                )
  call getopt( 'h0'   , is_exist, h0,     D_H0                )
  call getopt( 'f1'   , is_exist, f1,     D_F1                )
  call getopt( 'h1'   , is_exist, h1,     D_H1                )
  call getopt( 'decon', is_decon                              )
  call getopt( 'int'  , is_integ                              )
  call getopt( 'fint' , is_exist, fint,   D_FINT              )
  call getopt( 'hint' , is_exist, hint,   D_HINT              )
  call getopt( 'l'    , is_list ) 


  if( is_list ) then

     !! List mode
     !! 
     !! In this mode, input file fn_in is assumed to be an ascii-file that list-ups sac files.
     !! All sac files are assumed to have common sampling interval.

     call std__getio( io )
     open( io, file=trim(fn_in), action='read', status='old', iostat=ierr )
     
     if( ierr /= 0 ) then
        write(STDERR,*) "File "//trim(fn_in)//" not found."
        call usage_exit()
     end if

     !! in this mode, output filename is automatically determined. -o option does not work.
     fn_out = trim(fn_in) // '.dec'
     call std__countline( io, nfile )
     read(io,'(A)') fn_sac

  else

     !! Single file mode.
     !!
     !! Command-line orgument is SAC filename
     
     fn_sac = fn_in
     
  end if

  !! Read (first) SAC file with header
  call sac__read( fn_sac, sh, wav )
  call rtrend( sh%npts, wav )

  !!
  !! Filter design
  !!
  if( is_decon ) then
     call tdecon__tf_decon  ( f0, h0,         sh%delta, a, b ) !< Only deconvlolution
  else
     call tdecon__tf_simfilt( f0, h0, f1, h1, sh%delta, a, b ) !< Simluation Filter: Deconvolution & Convolution
  end if
  
  !! Integration with respect to time (if specified)
  if( is_integ ) then
     call tdecon__tf_integ( fint, hint, sh%delta, ai, bi )
     sh%idep = sh%idep - 1 !< Change physical unit from vel->disp or acc->vel
  end if

  !!
  !! Filter application
  !!

  !! First file
  call tdecon__apply( sh%npts, a,  b,  wav )
  if( is_integ ) call tdecon__apply( sh%npts, ai, bi, wav )
  call sac__write( fn_out, sh, wav, .true. )

  !! List mode
  !!
  !! read additional files from the given file list
  if( is_list ) then
     
     do i=2, nfile
        read(io,'(A)') fn_sac
        fn_out = trim(fn_sac) // '.dec'
        call sac__read( fn_sac, sh, wav )
        call rtrend( sh%npts, wav )

        call tdecon__apply( sh%npts, a,  b,  wav )
        if( is_integ ) call tdecon__apply( sh%npts, ai, bi, wav )
        call sac__write( fn_out, sh, wav, .true. )

     end do

  end if

contains
  
  !! --------------------------------------------------------------------------------------------------------------------------- !!  
  subroutine usage_exit()
    
  write(STDERR,'(A)')
  write(STDERR,'(A)')  ' hinet_decon.x sacfile [-o <outfile>] [-f0 <freq>] [-h0 <damp>] [-f1 <freq>] [-h1 <damp>]'
  write(STDERR,'(A)')  '                       [-h1 <damp>] [-decon] [-int] [-fint <freq>] [-hint <freq>]'
  write(STDERR,'(A)') 
  write(STDERR,'(A)')  " options: "
  write(STDERR,'(A)')  "   [-o  <outfile>]  indicate output sac filename <sacfile.out> "
  write(STDERR,'(A)')  "   [-f0 <freq>   ]  eigen frequency of input signal <1.0> "
  write(STDERR,'(A)')  "   [-h0 <damp>   ]  damping constant of input signal <0.7> "
  write(STDERR,'(A)')  "   [-f1 <freq>   ]  eigen frequency of output signal <0.00833333> "
  write(STDERR,'(A)')  "   [-h1 <damp>   ]  damping constant of output signal <0.707> "
  write(STDERR,'(A)')  "   [-decon       ]  deconvolve without simulation seismometer "
  write(STDERR,'(A)')  "   [-int         ]  output displacement record by numerical integration "
  write(STDERR,'(A)')  "   [-fint <freq> ]  ( corner frequency of low-cut filter for integration "
  write(STDERR,'(A)')  "   [-hint <damp> ]  ( damping constant of low-cut filter for integration "
  write(STDERR,'(A)')

  stop
  
  end subroutine usage_exit
  !! --------------------------------------------------------------------------------------------------------------------------- !!  
  
  
end program hinet_decon
!! ----------------------------------------------------------------------------------------------------------------------------- !!
