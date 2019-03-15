MODULE constants
USE SECOND_PRECISION,  ONLY : dp, ik
IMPLICIT NONE
PUBLIC

real(dp), parameter   :: Na = 6.022140857d23  ! 1/mol Avogadro constant
real(dp), parameter   :: R  = 8.3144598       ! [J/K/mol] Universal gas constant
real(dp), parameter   :: kb = R/Na            ! [J/K] Boltzmann constant
real(dp), parameter   :: pi = ACOS(-1d0)      ! pi
real(dp), parameter   :: K0 = 273.15          ! [K] Zero degree celcius in K
integer(ik), parameter:: min_s = 60           ! [s] seconds in minute
integer(ik), parameter:: hour_s  = 3600       ! [s] seconds in hour
integer(ik), parameter:: day_s = 24*hour_s

! Container for input parameters for creating simulated datapoints
! Uses gaussian function as generator
type parametered_input
    real(dp) :: sigma       ! Standard deviation for the Gaussian
    real(dp) :: min_c       ! Minimum value for the parametrized concentration
                            ! OR constant value if max_c <= min_c
    real(dp) :: max_c       ! Peak value
    real(dp) :: peaktime    ! Time of peak value
    real(dp) :: omega       ! Angular frequency [hours] of modifying sine function
    real(dp) :: amplitude   ! Amplitude of modificaion
    LOGICAL  :: LOGSCALE    ! Is concetration scale logaritmically or linearily.
                            ! True for logaritmically scaled
end type parametered_input

CONTAINS


end MODULE constants
