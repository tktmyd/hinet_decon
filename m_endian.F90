!! ----------------------------------------------------------------------------------------------------------------------------- !!
!>
!! Detection and conversion of machine endian of binary data
!<
!! ----------------------------------------------------------------------------------------------------------------------------- !!
module m_endian

  !! -- Dependency
  use m_std

  !! -- Declarations
  implicit none
  private

  !! -- Public Procedures
  public :: endian__check
  public :: endian__change
  public :: endian__byteswap

  !! -- Internal Constants
  integer, parameter, public :: BIG_ENDIAN     =  0
  integer, parameter, public :: LITTLE_ENDIAN  =  1
  integer, parameter, public :: UNKNOWN_ENDIAN = -1

  
  !! --------------------------------------------------------------------------------------------------------------------------- !!
  !>
  !! Check machine endian
  !!
  !! call endian__check( endian )
  !! endian: integer or string
  !!
  !!   return integer value has: 
  !!   1: little_endian ,  0: big_endian , -1: unknown_endian as defined parameters in this module
  !<
  !! --
  interface endian__check
     
     module procedure checkendian_i, checkendian_a

  end interface
  !! --------------------------------------------------------------------------------------------------------------------------- !!
  
  !! --------------------------------------------------------------------------------------------------------------------------- !!
  !!>
  !! Change between big and little endians
  !!
  !! Example call chgEndian( var )
  !<
  !! --
  interface endian__change
     
     module procedure chgEndian_r,  &
                      chgEndian_i,  &
                      chgEndian_l,  &
                      chgEndian_a4, &
                      chgEndian_d,  &
                      chgEndian_c,  &
                      chgEndian_cd
  end interface
  !! --------------------------------------------------------------------------------------------------------------------------- !!
  
contains
  
  !! --------------------------------------------------------------------------------------------------------------------------- !!
  !>
  !!  Check machine endian and returns integer
  !< 
  subroutine checkendian_i( endian )

    !! -- Arguments
    integer, intent(out) :: endian

    !! -- Local parameters
    integer, parameter :: ascii_0 = 48 ! ascii codes of '0'
    integer, parameter :: ascii_1 = 49 ! ascii codes of '1'
    integer, parameter :: ascii_2 = 50 ! ascii codes of '2'
    integer, parameter :: ascii_3 = 51 ! ascii codes of '3'

    integer :: i, io
    character(4) :: wk

    !! ----
    
    ! test variable
    i = ascii_0 + ascii_1*(256**1) + ascii_2*(256**2) + ascii_3*(256**3)

    ! create temporary file
    call std__getio( io )

    ! old direct access mode
    !open( io, status='scratch', access='direct', recl = 4) ! for old machine, need to care on recl dependency
    !write( io, rec = 1 ) i
    !read ( io, rec = 1 ) wk

    ! stream i/o version
    open( io, status='scratch', access='stream' )
    write( io ) i
    rewind( io )
    read ( io ) wk
    
    select case ( wk )
    case( '0123' )
       endian =  little_endian
    case( '3210' )
       endian =  big_endian
    case default
       endian = unknown_endian
    end select

    close( io ) ! automatically delete file in scratch mode
    
  end subroutine checkendian_i
  !! --------------------------------------------------------------------------------------------------------------------------- !!
  
  !! --------------------------------------------------------------------------------------------------------------------------- !!
  !>
  !!  Check machine endian and returns characters "Little" "BIG" or "Unknown"
  !<
  !! --
  subroutine checkendian_a( endian )

    !! -- Arguments
    character(*), intent(out) :: endian
    
    integer :: iendian
    !! ----
    
    call checkendian_i( iendian )
    if( iendian == LITTLE_ENDIAN ) then
       endian = 'Little'
    else if ( iendian == BIG_ENDIAN ) then
       endian = 'Big'
    else
       endian = 'Unknown'
    end if
    
  end subroutine checkendian_a
  !! --------------------------------------------------------------------------------------------------------------------------- !!
  
  !! --------------------------------------------------------------------------------------------------------------------------- !!
  !>
  !! Swap byte order of character foo. Used by endian change routines
  !<
  !! --
  subroutine endian__byteSwap( foo, nbyte )

    !! -- Arguments
    integer,          intent(in)    :: nbyte ! must be even
    character(nbyte), intent(inout) :: foo

    integer  :: i, j
    character(1) :: wk
    !! ----
    
    do i=1, nbyte/2
       j = nbyte - i + 1
       wk       = foo(i:i)
       foo(i:i) = foo(j:j)
       foo(j:j) = wk
    end do
    
  end subroutine endian__byteSwap
  !! --------------------------------------------------------------------------------------------------------------------------- !!
  

  
  !! --------------------------------------------------------------------------------------------------------------------------- !!
  !>
  !! Change endian of real-value
  !<
  !! --
  subroutine chgEndian_r( var )

    !! -- Arguments
    real, intent(inout) :: var

    character(4) :: c
    !! ----
    
    c = transfer( var, c )
    call endian__byteswap( c, 4 )
    var = transfer( c, var )
    
  end subroutine chgEndian_r
  !! --------------------------------------------------------------------------------------------------------------------------- !!
  
  !! --------------------------------------------------------------------------------------------------------------------------- !!
  !>
  !! Change endian of integer value
  !<
  !! --
  subroutine chgEndian_i( var )

    !! -- Arguments
    integer, intent(inout) :: var

    character(4) :: c
    !! ----
    
    c = transfer( var, c )
    call endian__byteswap( c, 4 )
    var = transfer( c, var )
    
  end subroutine chgEndian_i
  !! --------------------------------------------------------------------------------------------------------------------------- !!
  
  !! --------------------------------------------------------------------------------------------------------------------------- !!
  !>
  !! Change endian of logical value
  !<
  !! --
  subroutine chgEndian_l( var )

    !! -- Arguments
    logical, intent(inout) :: var

    character(4) :: c
    !! ----

    c = transfer( var, c )
    call endian__byteswap( c, 4 )
    var = transfer( c, var )
    
  end subroutine chgEndian_l
  !! --------------------------------------------------------------------------------------------------------------------------- !!
  
  !! --------------------------------------------------------------------------------------------------------------------------- !!
  !>
  !! Change endian of 4-character
  !<
  !! --
  subroutine chgEndian_a4( var )

    !! -- Arguments
    character(4), intent(inout) :: var
    !! ----
    
    call endian__byteswap( var, 4 )
    
  end subroutine chgEndian_a4
  !! --------------------------------------------------------------------------------------------------------------------------- !!
  
  !! --------------------------------------------------------------------------------------------------------------------------- !!
  !>
  !! Change endian of double precision value
  !<
  !! --
  subroutine chgEndian_d( var )

    !! -- Arguments
    real(DP), intent(inout) :: var

    character(8) :: c
    !! ----

    c = transfer( var, c )
    call endian__byteswap( c, 8 )
    var = transfer( c, var )
    
  end subroutine chgEndian_d
  !! --------------------------------------------------------------------------------------------------------------------------- !!
  
  !! --------------------------------------------------------------------------------------------------------------------------- !!
  !>
  !! Change endian of complex value
  !<
  !! --
  subroutine chgEndian_c( var )

    !! -- Arguments
    complex, intent(inout) :: var

    character(8) :: c
    !! ----

    c = transfer( var, c )
    call endian__byteswap( c, 8 )
    var = transfer( c, var )
    
  end subroutine chgEndian_c
  !! --------------------------------------------------------------------------------------------------------------------------- !!
  
  
  !! --------------------------------------------------------------------------------------------------------------------------- !!
  !>
  !! Change endian of double precision complex value
  !<
  !! --
  subroutine chgEndian_cd( var )

    !! -- Arguments
    complex(DP), intent(inout) :: var

    character(16) :: c
    !! ----

    c = transfer( var, c )
    call endian__byteswap( c, 16 )
    var = transfer( c, var )
    
  end subroutine chgEndian_cd
  !! --------------------------------------------------------------------------------------------------------------------------- !!
  
end module m_endian
!! ----------------------------------------------------------------------------------------------------------------------------- !!
