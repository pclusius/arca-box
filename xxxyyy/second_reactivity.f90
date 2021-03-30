module second_reactivity

  implicit none

  private

  integer, parameter :: NREACTIVITY=5

  character(len=20), parameter :: reactivity_name(NREACTIVITY) = (/ &
    'rOH_inorg           ', &
    'rOH_inorg2          ', &
    'rOH_monoter         ', &
    'rOH_other           ', &
    'rNO3_monoter        ' &
 /)

  public :: calculate_reactivities, NREACTIVITY, reactivity_name

contains

  subroutine calculate_reactivities(CONC, reactivity)

    use second_Precision   ! dp
    use second_Parameters  ! NSPEC, ind_*
    use second_Global      ! TEMP, M, O2, N2, H2O, RO2, K values

    real(dp), intent(in) :: CONC(NSPEC)
    real(dp), intent(out) :: reactivity(NREACTIVITY)

    ! rOH_inorg
    reactivity(1) = &
      CONC(ind_O3)*(1.70D-12*EXP(-940/TEMP)) + CONC(ind_NO)*(KMT07) + CONC(ind_NO2)*(KMT08) + CONC(ind_SO2)*(KMT12) + CONC(ind_CO)*(KMT05)

    ! rOH_inorg2
    reactivity(2) = &
      CONC(ind_H2)*(7.7D-12*EXP(-2100/TEMP)) + CONC(ind_H2O2)*(2.9D-12*EXP(-160/TEMP)) + CONC(ind_HONO)*(2.5D-12*EXP(260/TEMP)) + CONC(ind_HO2NO2)*(3.2D-13*EXP(690/TEMP)*1.0) + CONC(ind_HNO3)*(KMT11) &
      + CONC(ind_NO3)*(2.0D-11)

