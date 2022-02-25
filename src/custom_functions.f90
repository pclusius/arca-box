! ===================================================================================
! Atmospherically Relevant Chemistry and Aerosol box model
! Copyright (C) 2021  Multi-Scale Modelling group
! Institute for Atmospheric and Earth System Research (INAR), University of Helsinki
! Contact information arca@helsinki.fi
!
! This program is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with this program.  If not, see <https://www.gnu.org/licenses/>.
! ===================================================================================


module custom_functions

USE input                               ! input
USE second_MAIN                         ! Main second file
USE second_PARAMETERS                   ! CH_NSPEC and chemical indices, ind_xxxx, come from here
USE second_Precision,  ONLY : dp        ! KPP Numerical type
USE second_Monitor,    ONLY : SPC_NAMES ! Names of chemicals from KPP
USE constants
USE AUXILLARIES
implicit NONE

! ================================================================================================
! -----------------------------------  INDICES IN TSTEP_CONC -------------------------------------
! ------------------------------------------------------------------------------------------------
! 'TEMPK'         = inm_TempK
! 'PRESSURE'      = inm_pres
! 'REL_HUMIDITY'  = inm_RH
! 'CONDENS_SINK'  = inm_CS
! 'CON_SIN_NITR'  = inm_CS_NA
! 'SW_RADIATION'  = inm_swr
! 'ION_PROD_RATE' = inm_IPR
! 'H2SO4'         = inm_H2SO4
! 'NH3'           = inm_NH3
! 'DMA'           = inm_DMA
! 'SO2'           = inm_SO2
! 'NO'            = inm_NO
! 'NO2'           = inm_NO2
! 'CO'            = inm_CO
! 'H2'            = inm_H2
! 'O3'            = inm_O3
! 'NUC_RATE_IN '  = inm_JIN
!
! For other variables in TSTEP_CON, use the function IndexFromName(NAME) where NAME is
! the name as a string as it appears in the GUI (e.g. "APINENE").
! ------------------------------------------------------------------------------------------------


! ================================================================================================
! --------------------------------------  INDICES IN CH_GAS --------------------------------------
! ------------------------------------------------------------------------------------------------
! In KPP the indexing logic is so that any COMPOUNDS index in ind_COMPOUND. You can also get the
! index by using the function IndexFromName(NAME, list_of_names) where NAME is the name as a string,
! e.g. "APINENE", and list_of_names is the list of names whose index is searched for. In the case
! of Chemistry module, this list is called SPC_NAMES. So for example to find the index of PROPACID,
! use IndexFromName('PROPACID', SPC_NAMES). This integer would correspond to ind_PROPACID.
! ------------------------------------------------------------------------------------------------

contains

subroutine AFTER_CHEM(TSTEP_CONC,CH_GAS,CH_GAS_old,CH_RO2,CH_Beta,CH_H2O,CH_RO2_old,J_TOTAL_M3)
    implicit none
    real(dp), INTENT(INOUT) :: TSTEP_CONC(:)
    real(dp), INTENT(INOUT) :: CH_GAS(:)
    real(dp), INTENT(INOUT) :: CH_GAS_old(:)
    real(dp), INTENT(INOUT) :: CH_RO2,CH_Beta,CH_H2O
    real(dp), INTENT(INOUT) :: CH_RO2_old
    real(dp), INTENT(INOUT) :: J_TOTAL_M3

    ! BEGIN CUSTOM CODE -------------------------------
    ! CH_GAS(IND_O3) = TSTEP_CONC(inm_O3)
    ! IF (CH_GAS(IND_O3)>1d0) THEN
    !     IF (CH_GAS(IND_O3)>250*1d-9*GC_AIR_NOW) THEN
    !         FLOAT_EMIS_AFTER_HRS = 0d0
    !         CH_GAS(ind_emiO3) = 0d0
    !     ELSE IF (CH_GAS(IND_O3)<245*1d-9*GC_AIR_NOW.and.GTIME%hrs>FLOAT_EMIS_AFTER_HRS) THEN
    !          CH_GAS(ind_emiO3) = 25*(1d-12*GC_AIR_NOW)
    !     END IF
    ! END IF
    ! if (gtime%printnow) print'(es12.3,es12.3,es12.3)', FLOAT_EMIS_AFTER_HRS, CH_GAS(IND_O3)/GC_AIR_NOW*1d9,CH_GAS(ind_emiO3)

    ! END CUSTOM CODE ---------------------------------

end subroutine AFTER_CHEM

subroutine AFTER_NUCL(TSTEP_CONC,CH_GAS,J_TOTAL_M3)
    implicit none
    real(dp), INTENT(INOUT) :: TSTEP_CONC(:)
    real(dp), INTENT(INOUT) :: CH_GAS(:)
    real(dp), INTENT(INOUT) :: J_TOTAL_M3

    ! BEGIN CUSTOM CODE -------------------------------

    ! END CUSTOM CODE ---------------------------------

end subroutine AFTER_NUCL

end module custom_functions
