!! ----------------------------------------------------------------------------------------------------------------------------- !!
!>
!! Time-domain deconvolution filter
!!
!! Designs transfer function in Z-domain for seismometer instrumental response (and optionally integration with respect to time),
!! and applys it as time-domain digital filter
!!
!! Tranfer function notation in this routine is written as follows.
!!
!!           a(0) + a(1)*(1/z) + a(2)*(1/z)**2
!!   T(z) = -----------------------------------
!!           b(0) + b(1)*(1/z) + b(2)*(1/z)**2
!!
!! @see
!!   Maeda, T. et al. (2011), J. Geophys. Res. doi:10.1029/2011JB008464
!<
!! ----------------------------------------------------------------------------------------------------------------------------- !!
module m_tdecon

  !! -- Dependency
  use m_std

  !! -- Declarations
  implicit none
  private

  !! -- Public Procedures
  public :: tdecon__tf_simfilt
  public :: tdecon__tf_conv
  public :: tdecon__tf_decon
  public :: tdecon__tf_integ
  public :: tdecon__apply

  !! --------------------------------------------------------------------------------------------------------------------------- !!
  !>
  !! Filter design of seismometer response simulation filter
  !<
  !! --
  interface tdecon__tf_simfilt
     module procedure tf_simfilt_d, tf_simfilt_s
  end interface
  !! --------------------------------------------------------------------------------------------------------------------------- !!

  !! --------------------------------------------------------------------------------------------------------------------------- !!
  !>
  !! Filter design of seismometer response
  !<
  !! --
  interface tdecon__tf_conv
     module procedure tf_conv_d, tf_conv_s
  end interface
  !! --------------------------------------------------------------------------------------------------------------------------- !!

  !! --------------------------------------------------------------------------------------------------------------------------- !!
  !>
  !! Filter design of inverse seismometer response
  !<
  !! --
  interface tdecon__tf_decon
     module procedure tf_decon_d, tf_decon_s
  end interface
  !! --------------------------------------------------------------------------------------------------------------------------- !!

  !! --------------------------------------------------------------------------------------------------------------------------- !!
  !>
  !! Filter design for integration with time as a form of filter
  !<
  !!--
  interface tdecon__tf_integ
     module procedure tf_integ_d, tf_integ_s
  end interface
  !! --------------------------------------------------------------------------------------------------------------------------- !!

  !! --------------------------------------------------------------------------------------------------------------------------- !!
  !>
  !! Application of time-domain digital filter with coefficients of transfer function a(:) and b(:)
  !<
  !! --
  interface tdecon__apply
     module procedure apply_d, apply_s
  end interface
  !! --------------------------------------------------------------------------------------------------------------------------- !!
  

contains

  !! --------------------------------------------------------------------------------------------------------------------------- !!
  !>
  !! Filter design of seismometer response simulation filter
  !!
  !! Returns coefficients of transfer function of 
  !! Convolution/Deconvolution simulation filter to change instrument response from one to another by time-domain filtering
  !<
  !! --
  subroutine tf_simfilt_d ( f0, h0, f1, h1, dt, a, b )

    !! --Arguments
    real(DP), intent(in)  :: f0     !< eigenfrequency of input instrument
    real(DP), intent(in)  :: h0     !< damping factor of input instrument
    real(DP), intent(in)  :: f1     !< eigenfrequency of output instrument
    real(DP), intent(in)  :: h1     !< damping factor of output instrument
    real(DP), intent(in)  :: dt     !< sampling interval
    real(DP), intent(out) :: a(0:2) !< numerator coefficient
    real(DP), intent(out) :: b(0:2) !< denominator coefficient

    real(DP) :: tw0
    real(DP) :: tw1

    !--
    
    tw0 = tan( PI_D * f0 * dt )
    tw1 = tan( PI_D * f1 * dt )

    a(0) =  1.0_DP + 2*h0*tw0 +   tw0*tw0
    a(1) = -2.0_DP            + 2*tw0*tw0 
    a(2) =  1.0_DP - 2*h0*tw0 +   tw0*tw0

    b(0) =  1.0_DP + 2*h1*tw1 +   tw1*tw1
    b(1) = -2.0_DP            + 2*tw1*tw1 
    b(2) =  1.0_DP - 2*h1*tw1 +   tw1*tw1

  end subroutine tf_simfilt_d
  !! --------------------------------------------------------------------------------------------------------------------------- !!

  !! --------------------------------------------------------------------------------------------------------------------------- !!
  !>
  !! Filter design of seismometer response simulation filter
  !!
  !! Returns coefficients of transfer function of 
  !! Convolution/Deconvolution simulation filter to change instrument response from one to another by time-domain filtering
  !<
  !! --
  subroutine tf_simfilt_s ( f0, h0, f1, h1, dt, a, b )

    !! --Arguments
    real(SP), intent(in)  :: f0     !< eigenfrequency of input instrument
    real(SP), intent(in)  :: h0     !< damping factor of input instrument
    real(SP), intent(in)  :: f1     !< eigenfrequency of output instrument
    real(SP), intent(in)  :: h1     !< damping factor of output instrument
    real(SP), intent(in)  :: dt     !< sampling interval
    real(SP), intent(out) :: a(0:2) !< numerator coefficient
    real(SP), intent(out) :: b(0:2) !< denominator coefficient

    real(DP) :: tw0
    real(DP) :: tw1

    !--
    
    tw0 = tan( PI_D * f0 * dt )
    tw1 = tan( PI_D * f1 * dt )

    a(0) =  1.0 + 2*h0*tw0 +   tw0*tw0
    a(1) = -2.0            + 2*tw0*tw0 
    a(2) =  1.0 - 2*h0*tw0 +   tw0*tw0

    b(0) =  1.0 + 2*h1*tw1 +   tw1*tw1
    b(1) = -2.0            + 2*tw1*tw1 
    b(2) =  1.0 - 2*h1*tw1 +   tw1*tw1

  end subroutine tf_simfilt_s
  !! --------------------------------------------------------------------------------------------------------------------------- !!
  

  !! --------------------------------------------------------------------------------------------------------------------------- !!
  !>
  !! Filter design of seismometer response
  !!
  !! Returns coefficients of transfer function of seismometer for time-domain filtering
  !<
  !! --
  subroutine tf_conv_d ( f0, h, dt, a, b )

    !! -- Arguments 
    real(DP), intent(in)  :: f0     !< eigenfrequency of input instrument
    real(DP), intent(in)  :: h      !< damping factor of input instrument
    real(DP), intent(in)  :: dt     !< sampling interval
    real(DP), intent(out) :: a(0:2) !< numerator coefficient
    real(DP), intent(out) :: b(0:2) !< denominator coefficient

    real(DP) :: w0
    real(DP) :: tw
    !! ----
    
    w0 = 2 * PI_D * f0
    tw = tan( PI_D * f0 * dt )

    a(0) =  1.0_DP
    a(1) = -2.0_DP
    a(2) =  1.0_DP
    b(0) =  1.0_DP + 2*h*tw  +   tw*tw
    b(1) = -2.0_DP           + 2*tw*tw
    b(2) =  1.0_DP - 2*h*tw  +   tw*tw

  end subroutine tf_conv_d
  !! --------------------------------------------------------------------------------------------------------------------------- !!

  !! --------------------------------------------------------------------------------------------------------------------------- !!
  !>
  !! Filter design of seismometer response
  !!
  !! Returns coefficients of transfer function of seismometer for time-domain filtering
  !<
  !! --
  subroutine tf_conv_s( f0, h, dt, a, b )

    !! -- Arguments 
    real(SP), intent(in)  :: f0     !< eigenfrequency of input instrument
    real(SP), intent(in)  :: h      !< damping factor of input instrument
    real(SP), intent(in)  :: dt     !< sampling interval
    real(SP), intent(out) :: a(0:2) !< numerator coefficient
    real(SP), intent(out) :: b(0:2) !< denominator coefficient

    real(DP) :: w0
    real(DP) :: tw
    !! ----
    
    w0 = 2 * PI_D * f0
    tw = tan( PI_D * f0 * dt )

    a(0) =  1.0_DP
    a(1) = -2.0_DP
    a(2) =  1.0_DP
    b(0) =  1.0_DP + 2*h*tw  +   tw*tw
    b(1) = -2.0_DP           + 2*tw*tw
    b(2) =  1.0_DP - 2*h*tw  +   tw*tw

  end subroutine tf_conv_s
  !! --------------------------------------------------------------------------------------------------------------------------- !!
  
  !! --------------------------------------------------------------------------------------------------------------------------- !!
  !>
  !! Filter design of inverse seismometer response
  !!
  !! Returns coefficients of transfer function of inverse seismometer response
  !<
  !! --
  subroutine tf_decon_d( f0, h, dt, a, b )

    !! -- Arguments
    real(DP), intent(in)  :: f0     !< eigenfrequency of input instrument
    real(DP), intent(in)  :: h      !< damping factor of input instrument
    real(DP), intent(in)  :: dt     !< sampling interval
    real(DP), intent(out) :: a(0:2) !< numerator coefficient
    real(DP), intent(out) :: b(0:2) !< denominator coefficient

    real(DP) :: w0
    real(DP) :: tw
    !! ----

    w0 = 2 * PI_D * f0
    tw = tan( PI_D * f0 * dt )

    a(0) =  1.0_DP + 2*h*tw  +   tw*tw
    a(1) = -2.0_DP           + 2*tw*tw
    a(2) =  1.0_DP - 2*h*tw  +   tw*tw
    b(0) =  1.0_DP
    b(1) = -2.0_DP
    b(2) =  1.0_DP

  end subroutine tf_decon_d
  !! --------------------------------------------------------------------------------------------------------------------------- !!

  !! --------------------------------------------------------------------------------------------------------------------------- !!
  !>
  !! Filter design of inverse seismometer response
  !!
  !! Returns coefficients of transfer function of inverse seismometer response
  !<
  !! --
  subroutine tf_decon_s( f0, h, dt, a, b )

    !! -- Arguments
    real(SP), intent(in)  :: f0     !< eigenfrequency of input instrument
    real(SP), intent(in)  :: h      !< damping factor of input instrument
    real(SP), intent(in)  :: dt     !< sampling interval
    real(SP), intent(out) :: a(0:2) !< numerator coefficient
    real(SP), intent(out) :: b(0:2) !< denominator coefficient

    real(DP) :: w0
    real(DP) :: tw
    !! ----

    w0 = 2 * PI_D * f0
    tw = tan( PI_D * f0 * dt )

    a(0) =  1.0_DP + 2*h*tw  +   tw*tw
    a(1) = -2.0_DP           + 2*tw*tw
    a(2) =  1.0_DP - 2*h*tw  +   tw*tw
    b(0) =  1.0_DP
    b(1) = -2.0_DP
    b(2) =  1.0_DP

  end subroutine tf_decon_s
  !! --------------------------------------------------------------------------------------------------------------------------- !!

  
  !! --------------------------------------------------------------------------------------------------------------------------- !!
  !>
  !! Filter design for integration with time as a form of filter
  !!
  !! Returns coefficients of transfer function of time-domain integration with cut-off corner frequency f0 and damping h
  !<
  !! --
  subroutine tf_integ_d ( f0, h, dt, a, b )

    !! -- Arguments
    real(DP), intent(in)  :: f0     !< eigenfrequency of input instrument
    real(DP), intent(in)  :: h      !< damping factor of input instrument
    real(DP), intent(in)  :: dt     !< sampling interval
    real(DP), intent(out) :: a(0:2) !< numerator coefficient
    real(DP), intent(out) :: b(0:2) !< denominator coefficient

    real(DP) :: w0
    real(DP) :: tw
    !! ----

    w0 = 2 * PI_D * f0
    tw = tan( PI_D * f0 * dt )

    a(0) =  dt/2
    a(1) =  0.0_DP
    a(2) = -dt/2

    b(0) =  1.0_DP + 2*h*tw  +   tw*tw
    b(1) = -2.0_DP           + 2*tw*tw
    b(2) =  1.0_DP - 2*h*tw  +   tw*tw

  end subroutine tf_integ_d
  !! --------------------------------------------------------------------------------------------------------------------------- !!

  !! --------------------------------------------------------------------------------------------------------------------------- !!
  !>
  !! Filter design for integration with time as a form of filter
  !!
  !! Returns coefficients of transfer function of time-domain integration with cut-off corner frequency f0 and damping h
  !<
  !! --
  subroutine tf_integ_s( f0, h, dt, a, b )

    !! -- Arguments
    real(SP), intent(in)  :: f0     !< eigenfrequency of input instrument
    real(SP), intent(in)  :: h      !< damping factor of input instrument
    real(SP), intent(in)  :: dt     !< sampling interval
    real(SP), intent(out) :: a(0:2) !< numerator coefficient
    real(SP), intent(out) :: b(0:2) !< denominator coefficient

    real(DP) :: w0
    real(DP) :: tw
    !! ----

    w0 = 2 * PI_D * f0
    tw = tan( PI_D * f0 * dt )

    a(0) =  dt/2
    a(1) =  0.0_DP
    a(2) = -dt/2

    b(0) =  1.0_DP + 2*h*tw  +   tw*tw
    b(1) = -2.0_DP           + 2*tw*tw
    b(2) =  1.0_DP - 2*h*tw  +   tw*tw

  end subroutine tf_integ_s
  !! --------------------------------------------------------------------------------------------------------------------------- !!
  

  !! --------------------------------------------------------------------------------------------------------------------------- !!
  !>
  !! Application of time-domain digital filter with coefficients of transfer function a(:) and b(:)
  !<
  !! --
  subroutine apply_d( n, a, b, y )

    integer,  intent(in)    :: n      !< number of the data sample
    real(DP), intent(in)    :: a(0:2) !< transfer function numeretor coefficients
    real(DP), intent(in)    :: b(0:2) !< transfer funciton denominator coefficients
    real(DP), intent(inout) :: y(n)   !< time series: input will be overwritten by filtered trace

    real(DP) :: a0, a1, a2
    real(DP) :: b1, b2
    real(DP) :: y1, y2, x1, x2
    real(DP) :: output
    integer  :: i
    
    !! ----

    a0 = a(0) / b(0)
    a1 = a(1) / b(0)
    a2 = a(2) / b(0)
    
    b1 = b(1) / b(0)
    b2 = b(2) / b(0)

    y1 = 0.0_DP
    y2 = 0.0_DP
    x1 = 0.0_DP
    x2 = 0.0_DP
    
    do i=1, n

       output = a0*y(i) + a1*x1 + a2*x2 - b1*y1 - b2*y2
       y2 = y1
       y1 = output
       x2 = x1
       x1 = y(i)
       y(i) = output
       
    end do
    
  end subroutine apply_d
  !! --------------------------------------------------------------------------------------------------------------------------- !!
  
  !! --------------------------------------------------------------------------------------------------------------------------- !!
  !>
  !! Application of time-domain digital filter with coefficients of transfer function a(:) and b(:)
  !!
  !! This routine works as the interface to double-precision routine apply_d.
  !<
  !! --
  subroutine apply_s( n, a, b, y )

    integer,  intent(in)    :: n      !< number of the data sample
    real(SP), intent(in)    :: a(0:2) !< transfer function numeretor coefficients
    real(SP), intent(in)    :: b(0:2) !< transfer funciton denominator coefficients
    real(SP), intent(inout) :: y(n)   !< time series: input will be overwritten by filtered trace
    real(DP) :: aa(0:2), bb(0:2)
    real(DP) :: yy(n)

    aa(0:2) = a(0:2)
    bb(0:2) = b(0:2)
    
    call apply_d( n, aa, bb, yy )

    y(1:n) = real(yy(1:n))

  end subroutine apply_s
  !! --------------------------------------------------------------------------------------------------------------------------- !!
  
end module m_tdecon
!! ----------------------------------------------------------------------------------------------------------------------------- !!
