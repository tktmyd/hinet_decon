!+-----------------------------------------------------------------------------!
module m_rtrend
  !
  !=Description
  ! 時系列のトレンド除去
  !
  !=Version
  ! $Id: m_rtrend.F90 4abbb83a8c5f 2014/09/07 07:18:37 maeda $
  !
  !=Declarations
  use m_std
  implicit none
  public :: rtrend
  !+
  !--

  interface rtrend
     module procedure  rtrend_s, rtrend_d
  end interface rtrend
  

contains
  
  !+---------------------------------------------------------------------------!
  subroutine rtrend_d( n, y )
    !
    !=Description
    ! Remove trendline 
    !
    !=Arguments
    integer,  intent(in)    :: n
    real(DP), intent(inout) :: y(n)
    !+
    !=Variables
    integer ::  i
    real(DP) :: sum_y,sum_xy, sum_xx, sum_x, delta, a, b
    real(DP) :: x(n)
    !--
    ! 1d fitting
    do i=1,n
       x(i) = dble(i) / dble(n)
    end do
    
    sum_x  = sum(x(1:n))
    sum_xx = sum(x(1:n)*x(1:n))
    sum_y  = sum(y(1:n))
    sum_xy = sum(x(1:n)*y(1:n))
    
    delta = n * sum_xx - sum_x * sum_x
    a = ( n * sum_xy - sum_y * sum_x ) / delta
    b = ( sum_y * sum_xx - sum_x * sum_xy ) / delta
    
    ! remove trend
    do i=1, n
       y(i) = y(i) - a*x(i)-b
    end do
    
  end subroutine rtrend_d
  !----------------------------------------------------------------------------!
  
  
  !+---------------------------------------------------------------------------!
  subroutine rtrend_s( n, y )
    !
    !=Description
    ! Remove trendline 
    !
    !=Arguments
    integer,  intent(in)    :: n
    real(SP), intent(inout) :: y(n)
    !+
    !=Variables
    integer ::  i
    real(SP) :: sum_y,sum_xy, sum_xx, sum_x, delta, a, b
    real(SP) :: x(n)
    !--
    ! 1d fitting
    do i=1,n
       x(i) = dble(i) / dble(n)
    end do
    
    sum_x  = sum(x(1:n))
    sum_xx = sum(x(1:n)*x(1:n))
    sum_y  = sum(y(1:n))
    sum_xy = sum(x(1:n)*y(1:n))
    
    delta = n * sum_xx - sum_x * sum_x
    a = ( n * sum_xy - sum_y * sum_x ) / delta
    b = ( sum_y * sum_xx - sum_x * sum_xy ) / delta
    
    ! remove trend
    do i=1, n
       y(i) = y(i) - a*x(i)-b
    end do
    
  end subroutine rtrend_s
  !----------------------------------------------------------------------------!

end module m_rtrend
!------------------------------------------------------------------------------!
