module custom_functions

USE input                               ! input
USE second_MAIN                         ! Main second file
USE second_PARAMETERS                   ! CH_NSPEC (originally NSPEC) and chemical indices, ind_xxxx, come from here
USE second_Precision,  ONLY : dp        ! KPP Numerical type
USE second_Monitor,    ONLY : SPC_NAMES ! Names of chemicals from KPP
USE constants
USE AUXILLARIES
implicit NONE

contains
subroutine AFTER_CHEM(TSTEP_CONC,CH_GAS,J_TOTAL_M3)
    implicit none
    real(dp) :: TSTEP_CONC(:)
    real(dp) :: CH_GAS(:)
    real(dp) :: J_TOTAL_M3
    ! -------------------------------
end subroutine AFTER_CHEM

subroutine AFTER_NUCL(TSTEP_CONC,CH_GAS,J_TOTAL_M3)
    implicit none
    real(dp) :: TSTEP_CONC(:)
    real(dp) :: CH_GAS(:)
    real(dp) :: J_TOTAL_M3
    ! -------------------------------
end subroutine AFTER_NUCL

end module custom_functions
