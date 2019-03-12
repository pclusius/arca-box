MODULE constants
USE SECOND_PRECISION,  ONLY : dp, ik
IMPLICIT NONE
PUBLIC

real(dp), parameter   :: Na = 6.022140857d23  ! 1/mol Avogadro constant
real(dp), parameter   :: R  = 8.3144598       ! [J/K/mol] Universal gas constant
real(dp), parameter   :: kb = R/Na            ! [J/K] Boltzmann constant
real(dp), parameter   :: pi = ACOS(-1d0)         ! pi
integer(ik), parameter:: min_s = 60           ! [s] seconds in minute
integer(ik), parameter:: hour_s  = 3600       ! [s] seconds in hour
integer(ik), parameter:: day_s = 24*hour_s

CONTAINS
  CHARACTER(8) function time_hms(s)
    real(dp) :: s
    write(time_hms, '(i2.2, ":" i2.2, ":" i2.2)') int(s+0.5_dp)/3600_ik, &
      int(MODULO(int(s+0.5_dp),3600_ik)/60_ik), MODULO(MODULO(int(s+0.5_dp),3600_ik), 60_ik)
  end function time_hms

end MODULE constants
