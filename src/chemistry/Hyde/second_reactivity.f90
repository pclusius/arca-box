module second_reactivity

  implicit none

  private

  integer, parameter :: NREACTIVITY=0

  character(len=20), parameter :: reactivity_name(NREACTIVITY) = (/ &
REAL:: &
 /)

  public :: calculate_reactivities, NREACTIVITY, reactivity_name

contains

  subroutine calculate_reactivities(CONC, reactivity)

    use second_Precision   ! dp
    use second_Parameters  ! NSPEC, ind_*
    use second_Global      ! TEMP, M, O2, N2, H2O, RO2, K values

    real(dp), intent(in) :: CONC(NSPEC)
    real(dp), intent(out) :: reactivity(NREACTIVITY)

  end subroutine calculate_reactivities

end module second_reactivity