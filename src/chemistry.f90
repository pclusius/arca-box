!ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
!
!
!
!   CHEMISTRY MODULE
!
!
!   For calculations of all J-values
!
!   This version is for MCM version 3.3.1
!
!ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

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



MODULE CHEMISTRY


    !Include the following second_*.f90 files that are specific for different chemistries:
    USE SECOND_MAIN
    USE SECOND_PRECISION,  ONLY : DP           ! KPP Numerical type
    USE SECOND_PARAMETERS, ONLY : NSPEC        ! Number of chemical species
    USE SECOND_GLOBAL,     ONLY : NPHOT, RO2   ! NPHOT = Number of photochemical reactions used in KPP
                                               ! RO2 = organic peroxy radicals
    USE SECOND_REACTIVITY
    USE AUXILLARIES,  ONLY : i2chr, COLCOUNT   ! helper

    IMPLICIT NONE

    PRIVATE

    ! Public subroutines which will be called in the main program to start the chemistry
    PUBLIC :: CHEMCALC

    ! Making the J_values and number of complex reaction reates public so they are available in the main program
    PUBLIC :: J_values!, NKVALUES

    ! Global variables:
    REAL(DP), DIMENSION(NPHOT) :: J_values

!    INTEGER, PARAMETER :: NKVALUES = 42   ! Number of rate coefficients used in KPP. Hand-copied from second_Main.f90

    ! File path for the input values:
    CHARACTER(LEN=*), PARAMETER :: data_root = 'ModelLib'

CONTAINS



!oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
!
!   CHEMCALC
!
!   Calculates the chemically related properties for KPP
!
!   This subroutine includes three subroutines:
!
!        - RATECONS
!        - PHOTOLYSIS_KPP
!        - KPP_PROCEED
!
!ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo



SUBROUTINE CHEMCALC(CONS, TIME_kpp, END_kpp, Tp, SWDWN, Beta, H2O, AIR, RES1, RES2, Albedo, RO2_out,reactivities,swr_dis,swr_is_af)

    IMPLICIT NONE

    ! Global variables for chemistry module:

    INTEGER                        :: J1, I
    LOGICAL, SAVE                  :: first_call = .TRUE.

    LOGICAL, INTENT(in)            :: swr_is_af

    REAL(dp), INTENT(inout)        :: CONS(NSPEC)                ! Concentrations of all gas compounds

    REAL(dp), INTENT(out)          :: RO2_out                    ! Concentration of peroxy radical

    REAL(dp), DIMENSION(75)        :: I_ACF                      ! Calculated actinic flux from 280-650nm
    ! REAL(dp), DIMENSION(84), SAVE  :: swr_dis                    ! Distribution of solar radiation measured at SMEAR II

    REAL(dp), INTENT(in)           :: TIME_kpp, END_kpp,       & ! Model time in sec and same time plus chemical timestep
                                      H2O,                     & ! Water vapour concentration
                                      AIR,                     & ! air molecules - calculated in the main file
                                      Beta,                    & ! Solar elevation angle
                                      Tp,                      & ! Temperature [K]
                                      SWDWN,                   & ! Global short wave radiation
                                      RES1, RES2,              & ! Condensation sink for H2SO4 and HNO3
                                      Albedo,                  & ! Albedo
                                      swr_dis(:)                 ! Spectral data - either distribution or absolute values

    REAL(dp)                       :: O2,                      & ! O2-concentration in molecule / cm3
                                      H2,                      & ! H2-concentration in molecule / cm3
                                      N2,                      & ! N2-concentration in molecule / cm3
                                      MO2N2,                   & ! MO2N2: Sum of O2 and N2 concentration
                                      ZSD                        ! Solar zenith angle
    REAL(DP)                       :: reactivities(:)

    ! Calculation of O2, N2, H2 and MO2N2 in molecule / cm3:
    O2        = 0.209460 * Air
    N2        = 0.780840 * Air
    H2        = 0.5 * Air / 1E6
    MO2N2     = O2 + N2
    ZSD       = 90. - MAX(Beta, 0d0)

    ! Actinic flux and J are zero by default
    I_ACF = 0.
    J_values = 0.

    ! Calculated actinic flux based on supplied global shortwave radiation
    IF (swr_is_af) THEN
        I_ACF(:) = SWDWN * swr_dis(1:75)
    ELSE
        I_ACF = 1d-4/(6.62607004D-34*2.99792458D17) * [(i,i=280,650,5)] ! conversion from W/m² to photons/cm²/s
        I_ACF = I_ACF * SWDWN * swr_dis(1:75) * (1 + 2 * Albedo * COS(ZSD*3.1416/180) + Albedo)
    END IF

    ! Calculation of the J-values(1:56) based on MCM
    CALL PHOTOLYSIS_KPP(I_ACF, TP, J_values)

    ! In the new chemistry from MCM no J-values above 56 are used
    ! But this is handled in PHOTOLYSIS_KPP
    ! J_values(57:) = 0

    ! P.C. commented out as unnecessary. SWR is measured and should not be put to zero.
    ! AF will not go to zero after sunset (though small it is)
    ! IF (ZSD .GE. 90. .AND. .NOT. swr_is_af) THEN
    !     DO J1 = 1,NPHOT
    !        J_values(J1) = 0.
    !     ENDDO
    ! ENDIF

    ! Some safety check by Sampo Smolander
    IF (MINVAL(J_values) < -1d9 .OR. MAXVAL(J_values) > 1d12 ) &
        STOP 'J-values are unphysical, check PHOTOLYSIS_KPP'

    ! Calling KPP to calculate the new concentrations of the chemical compounds
    CALL KPP_Proceed(CONS, TIME_kpp, END_kpp, Tp, O2, N2, MO2N2, H2O, RES1, RES2, J_values)

    ! RO2 redicals are summed and provided for output - not used anywhere else
    RO2_out = RO2

    if (NREACTIVITY>0) CALL calculate_reactivities(CONS, reactivities)

END SUBROUTINE CHEMCALC


!ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
!
!   Subroutine for the photodissociation reactions
!
!ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
SUBROUTINE PHOTOLYSIS_KPP(SAFN, T, J_pc)
    IMPLICIT NONE

    INTEGER :: I, ioi, JJ, cols, R

    INTEGER,  DIMENSION(75), SAVE      :: WL        ! Wavelengths from the input files
    REAL(dp), DIMENSION(56,75,7), SAVE :: CS_QY=0d0 ! cross section parameters. First rank indices refer to reaction rate
    REAL(dp), INTENT(IN)               :: SAFN(:)   ! spectral actinic flux
    REAL(dp), INTENT(INOUT)            :: J_pc(:)   ! Photochemistry rates (some are empty), second rank to parameters (some are empty)
    REAL(dp), INTENT(IN)               :: T         ! Temperature [K]
    REAL(dp), PARAMETER                :: DL=5d0    ! Delta lambda [nm]

    LOGICAL, SAVE :: first_call = .TRUE.

    ! The following is run only once, only on the first time CHEMCALC is called
    ! Input of absorption cross spectrum and quantum yields for different molecules:
    IF (first_call) THEN ! do only once
        first_call = .FALSE.

        CS_QY = 0d0
        DO I=1,56
            OPEN(900, FILE = data_root//'/Photolyse/MCM3_CS_QY/J'//TRIM(i2chr(I))//'.txt',STATUS = 'OLD', IOSTAT=ioi)
            if (ioi==0) THEN
                cols = COLCOUNT(900)
                DO JJ=1,75
                    READ(900,*)  WL(JJ), CS_QY(I,JJ,1:cols-1)
                END DO
                CLOSE(900)
            END IF
        END DO

    ENDIF ! do only once

    J_pc = 0d0
    ! Calculate J-values:
    DO I=1,75
        R=1;  J_pc(R) = J_pc(R) &
                        + CS_QY(R,I,1) * EXP(CS_QY(R,I,2)/T) &
                        * CS_QY(R,I,3) * SAFN(I) * DL

        R=2;  J_pc(R) = J_pc(R) &
                        + CS_QY(R,I,1) * EXP(CS_QY(R,I,2)/T) &
                        * CS_QY(R,I,3) * SAFN(I) * DL
        R=3;  J_pc(R) = J_pc(R) + CS_QY(R,I,1) * CS_QY(R,I,2) * SAFN(I) * DL
        R=4;  J_pc(R) = J_pc(R) + CS_QY(R,I,1) * CS_QY(R,I,2) * SAFN(I) * DL
        R=5;  J_pc(R) = J_pc(R) + CS_QY(R,I,1) * CS_QY(R,I,2) * SAFN(I) * DL
        R=6;  J_pc(R) = J_pc(R) + CS_QY(R,I,1) * CS_QY(R,I,2) * SAFN(I) * DL
        R=7;  J_pc(R) = J_pc(R) + CS_QY(R,I,1) * CS_QY(R,I,2) * SAFN(I) * DL
        ! (T-independent) ACS, QY=1
        R=8;  J_pc(R) = J_pc(R) + CS_QY(R,I,1) * CS_QY(R,I,4) * SAFN(I) * DL

        R=11; J_pc(R) = J_pc(R) &
                       + (CS_QY(R,I,1) + CS_QY(R,I,2) * (T-298)) &
                       * CS_QY(R,I,3) * SAFN(I) * DL

        R=12; J_pc(R) = J_pc(R) &
                        + (CS_QY(R,I,1) + CS_QY(R,I,2) * (T-298)) &
                        * CS_QY(R,I,3) * SAFN(I) * DL

        R=13; J_pc(R) = J_pc(R) + CS_QY(R,I,1) * CS_QY(R,I,2) * SAFN(I) * DL
        R=14; J_pc(R) = J_pc(R) + CS_QY(R,I,1) * CS_QY(R,I,2) * SAFN(I) * DL
        R=15; J_pc(R) = J_pc(R) + CS_QY(R,I,1) * CS_QY(R,I,2) * SAFN(I) * DL
        R=16; J_pc(R) = J_pc(R) + CS_QY(R,I,1) * CS_QY(R,I,2) * SAFN(I) * DL
        R=17; J_pc(R) = J_pc(R) + CS_QY(R,I,1) * CS_QY(R,I,2) * SAFN(I) * DL
        R=18; J_pc(R) = J_pc(R) + CS_QY(R,I,1) * CS_QY(R,I,2) * SAFN(I) * DL
        R=19; J_pc(R) = J_pc(R) + CS_QY(R,I,1) * CS_QY(R,I,2) * SAFN(I) * DL
        ! Reaction J20 was previously missing
        R=20; J_pc(R) = J_pc(R) + CS_QY(R,I,1) * CS_QY(R,I,2) * SAFN(I) * DL
        R=21; J_pc(R) = J_pc(R) &
                   + CS_QY(R,I,1) * (1 + (CS_QY(R,I,2) * T) + (CS_QY(R,I,3) * (T**2)) + (CS_QY(R,I,4) * (T**3))) &
                   ! Using summed quantum yields (from column 7)
                   * CS_QY(R,I,7) * SAFN(I) * DL

        R=22; J_pc(R) = J_pc(R) + CS_QY(R,I,1) * CS_QY(R,I,2) * SAFN(I) * DL
        R=23; J_pc(R) = J_pc(R) + CS_QY(R,I,1) * CS_QY(R,I,2) * SAFN(I) * DL
        R=24; J_pc(R) = J_pc(R) + CS_QY(R,I,1) * CS_QY(R,I,2) * SAFN(I) * DL
        R=31; J_pc(R) = J_pc(R) + CS_QY(R,I,1) * CS_QY(R,I,2) * SAFN(I) * DL
        R=32; J_pc(R) = J_pc(R) + CS_QY(R,I,1) * CS_QY(R,I,2) * SAFN(I) * DL
        R=33; J_pc(R) = J_pc(R) + CS_QY(R,I,1) * CS_QY(R,I,2) * SAFN(I) * DL
        R=34; J_pc(R) = J_pc(R) + CS_QY(R,I,1) * CS_QY(R,I,2) * SAFN(I) * DL
        R=35; J_pc(R) = J_pc(R) + CS_QY(R,I,1) * CS_QY(R,I,2) * SAFN(I) * DL
        R=41; J_pc(R) = J_pc(R) + CS_QY(R,I,1) * CS_QY(R,I,2) * SAFN(I) * DL

        R=51; J_pc(R) = J_pc(R) &
                   + CS_QY(R,I,1) * EXP(CS_QY(R,I,2) * (T-298)) &
                   * CS_QY(R,I,3) * SAFN(I) * DL

        R=52; J_pc(R) = J_pc(R) &
                   + CS_QY(R,I,1) * EXP(CS_QY(R,I,2) * (T-298)) &
                   * CS_QY(R,I,3) * SAFN(I) * DL

        R=53; J_pc(R) = J_pc(R) + CS_QY(R,I,1) * CS_QY(R,I,2) * SAFN(I) * DL

        R=54; J_pc(R) = J_pc(R) &
                   + CS_QY(R,I,1) * EXP(CS_QY(R,I,2) * (T-298)) &
                   * CS_QY(R,I,3) * SAFN(I) * DL

        R=55; J_pc(R) = J_pc(R) + CS_QY(R,I,1) * CS_QY(R,I,2) * SAFN(I) * DL
        R=56; J_pc(R) = J_pc(R) + CS_QY(R,I,1) * CS_QY(R,I,2) * SAFN(I) * DL

    ENDDO

END SUBROUTINE PHOTOLYSIS_KPP


END MODULE


! !ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
! !
! !
! !
! !   CHEMISTRY MODULE
! !
! !
! !   For calculations of all J-values
! !
! !   This version is for MCM version 3.3.1
! !
! !ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
!
! ! ===================================================================================
! ! Atmospherically Relevant Chemistry and Aerosol box model
! ! Copyright (C) 2021  Multi-Scale Modelling group
! ! Institute for Atmospheric and Earth System Research (INAR), University of Helsinki
! ! Contact information arca@helsinki.fi
! !
! ! This program is free software: you can redistribute it and/or modify
! ! it under the terms of the GNU General Public License as published by
! ! the Free Software Foundation, either version 3 of the License, or
! ! (at your option) any later version.
! !
! ! This program is distributed in the hope that it will be useful,
! ! but WITHOUT ANY WARRANTY; without even the implied warranty of
! ! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! ! GNU General Public License for more details.
! !
! ! You should have received a copy of the GNU General Public License
! ! along with this program.  If not, see <https://www.gnu.org/licenses/>.
! ! ===================================================================================
!
!
!
!
!
! MODULE CHEMISTRY
!
!
!     !Include the following second_*.f90 files that are specific for different chemistries:
!     USE SECOND_MAIN
!     USE SECOND_PRECISION,  ONLY : DP           ! KPP Numerical type
!     USE SECOND_PARAMETERS, ONLY : NSPEC        ! Number of chemical species
!     USE SECOND_GLOBAL,     ONLY : NPHOT, RO2   ! NPHOT = Number of photochemical reactions used in KPP
!                                                ! RO2 = organic peroxy radicals
!     USE SECOND_REACTIVITY
!
!     IMPLICIT NONE
!
!     PRIVATE
!
!     ! Public subroutines which will be called in the main program to start the chemistry
!     PUBLIC :: CHEMCALC
!
!     ! Making the J_values and number of complex reaction reates public so they are available in the main program
!     PUBLIC :: J_values!, NKVALUES
!
!     ! Global variables:
!     REAL(DP), DIMENSION(NPHOT) :: J_values
!
! !    INTEGER, PARAMETER :: NKVALUES = 42   ! Number of rate coefficients used in KPP. Hand-copied from second_Main.f90
!
!     ! File path for the input values:
!     CHARACTER(LEN=*), PARAMETER :: filename1 = 'ModelLib'
!
! CONTAINS
!
!
!
!     !oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
!     !
!     !   CHEMCALC
!     !
!     !   Calculates the chemically related properties for KPP
!     !
!     !   This subroutine includes three subroutines:
!     !
!     !        - RATECONS
!     !        - PHOTOLYSIS_KPP
!     !        - KPP_PROCEED
!     !
!     !ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
!
!
!
!     SUBROUTINE CHEMCALC(CONS, TIME_kpp, END_kpp, Tp, SWDWN, Beta, H2O, AIR, RES1, RES2, Albedo, RO2_out,reactivities,swr_dis,swr_is_af)
!
!         IMPLICIT NONE
!
!         ! Global variables for chemistry module:
!
!         INTEGER                        :: J1, I
!         LOGICAL, SAVE                  :: first_call = .TRUE.
!
!         LOGICAL, INTENT(in)            :: swr_is_af
!
!         REAL(dp), INTENT(inout)        :: CONS(NSPEC)                ! Concentrations of all gas compounds
!
!         REAL(dp), INTENT(out)          :: RO2_out                    ! Concentration of peroxy radical
!
!         REAL(dp), DIMENSION(75)        :: I_ACF                      ! Calculated actinic flux from 280-650nm
!         ! REAL(dp), DIMENSION(84), SAVE  :: swr_dis                    ! Distribution of solar radiation measured at SMEAR II
!
!         REAL(dp), INTENT(in)           :: TIME_kpp, END_kpp,       & ! Model time in sec and same time plus chemical timestep
!                                           H2O,                     & ! Water vapour concentration
!                                           AIR,                     & ! air molecules - calculated in the main file
!                                           Beta,                    & ! Solar elevation angle
!                                           Tp,                      & ! Temperature at this level of atmosphere. So in this box: temperature at ground level
!                                           SWDWN,                   & ! Global short wave radiation
!                                           RES1, RES2,              & ! Condensation sink for H2SO4 and HNO3
!                                           Albedo,                  & ! Albedo
!                                           swr_dis(:)                 ! Spectral data - either distribution or absolute values
!
!         REAL(dp)                       :: O2,                      & ! O2-concentration in molecule / cm3
!                                           H2,                      & ! H2-concentration in molecule / cm3
!                                           N2,                      & ! N2-concentration in molecule / cm3
!                                           MO2N2,                   & ! MO2N2: Sum of O2 and N2 concentration
!                                           ZSD                        ! Solar zenith angle
!         REAL(DP)                       :: reactivities(:)
!
!         ! CHARACTER(len=*), INTENT(in)   :: spectrumfile               ! Spectral function from 280-695 nm (min 84 rows)
!
! !        write(*,*) 'Gases: ',    CONS
! !        write(*,*) 'Time_kpp: ', TIME_kpp
! !        write(*,*) 'End_kpp: ',  END_kpp
! !        write(*,*) 'Temp: ',     Tp
! !        write(*,*) 'Global: ',   SWDWN
! !        write(*,*) 'Beta: ',     Beta
! !        write(*,*) 'H2O: ',      H2O
! !        write(*,*) 'Air: ',      AIR
! !        write(*,*) 'Res1: ',     RES1
! !        write(*,*) 'Res2: ',     RES2
! !        write(*,*) 'Albedo: ',   Albedo
! !        write(*,*) 'RO2: ',      RO2_out
!
!
!         ! Calculation of O2, N2, H2 and MO2N2 in molecule / cm3:
!         O2        = 0.209460 * Air
!         N2        = 0.780840 * Air
!         H2        = 0.5 * Air / 1E6
!         MO2N2     = O2 + N2
!         ZSD       = 90. - Beta
!
!         ! IF (first_call) THEN ! do only once
!         !     first_call = .FALSE.
!         !     ! For solar radiation distribution: measured data is read in for 280-500 nm but value at 500 nm is used for 505 - 650 nm
!         !     OPEN(960, file = TRIM(spectrumfile), status = "old")
!         !     DO J1 = 1,84
!         !         READ(960,*) swr_dis(J1)
!         !     ENDDO
!         !     CLOSE(960)
!         ! ENDIF
!
!         I_ACF = 0.       ! all values for actinic flux are set to zero at the first time step
!
!         ! Calculated actinic flux based on supplied global shortwave radiation
!         IF (swr_is_af) THEN
!             I_ACF(:) = SWDWN * swr_dis(1:75)
!         ELSE
!             I_ACF(:) = SWDWN * swr_dis(1:75) * (1 + 2 * Albedo * COS(ZSD*3.1416/180) + Albedo)
!         END IF
!
!         J_values = 0.    ! all J-values are set to zero at the first time step
!
!         ! Calculation of the J-values based on MCM; J-value 70 (HO2NO2) and 71 (N2O5) are added by Ditte
!         CALL PHOTOLYSIS_KPP     (I_ACF,        TP,           J_values(1),  J_values(2),  J_values(3),  J_values(4),  J_values(5),  J_values(6),      &
!                                  J_values(7),  J_values(8),  J_values(11), J_values(12), J_values(13), J_values(14), J_values(15), J_values(16),     &
!                                  J_values(17), J_values(18), J_values(19), J_values(21), J_values(22), J_values(23), J_values(24), J_values(31),     &
!                                  J_values(32), J_values(33), J_values(34), J_values(35), J_values(41), J_values(51), J_values(52), J_values(53),     &
!                                  J_values(54), J_values(55), J_values(56), J_values(57), J_values(67), J_values(68), J_values(70), J_values(71))
!
!         ! In the new chemistry from MCM no J-values above 56 are used
!         J_values(57:size(J_values,1)) = 0
!
!         IF (ZSD .GE. 90. .AND. .NOT. swr_is_af) THEN
!             DO J1 = 1,NPHOT
!                J_values(J1) = 0.
!             ENDDO
!         ENDIF
!
!         ! Some safety check by Sampo Smolander
!         DO I = 1,NPHOT
!             IF (J_values(I) < -1d9 .OR. J_values(I) > 1d12 ) then
!                 WRITE(*,*) 'Note by Sampo: J_value(I) has bad value.'
!                 WRITE(*,*) 'I, J_value(I) = ', I, J_values(I)
!                 STOP
!             ENDIF
!         ENDDO
!
!         ! Calling KPP to calculate the new concentrations of the chemical compounds
!         CALL KPP_Proceed(CONS, TIME_kpp, END_kpp, Tp, O2, N2, MO2N2, H2O, RES1, RES2, J_values)
!
!         ! RO2 redicals are summed and provided for output - not used anywhere else
!         RO2_out = RO2
!
!         if (NREACTIVITY>0) CALL calculate_reactivities(CONS, reactivities)
!
!     END SUBROUTINE
!
!
!
! !ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
! !
! !   Subroutine for the photodissociation reactions
! !
! !ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
!
!
!
!     SUBROUTINE PHOTOLYSIS_KPP(SAFN, T,    PR1,  PR2,  PR3,  PR4,  PR5,  PR6,  PR7,  PR8,    &
!                               PR11, PR12, PR13, PR14, PR15, PR16, PR17, PR18, PR19, PR21,   &
!                               PR22, PR23, PR24, PR31, PR32, PR33, PR34, PR35, PR41, PR51,   &
!                               PR52, PR53, PR54, PR55, PR56, PR57, PR67, PR68, PR70, PR71)
!
!     IMPLICIT NONE
!
!     INTEGER, PARAMETER :: real_x = SELECTED_REAL_KIND(p=14, r=300)
!
!     INTEGER I
!
!     INTEGER, DIMENSION(75), SAVE :: WL
!
!     REAL(real_x), DIMENSION(75) :: SAFN ! spectral actinic flux - based on hyytiala data
!
!     REAL(real_x) :: T,                                                                         & ! Temperature in K
!                     DL,                                                                        & ! Delta lambda
!                     VA, VAA, VAB, VAC, VAD,                                                    & ! Dimension parameter
!                     PR1, PR2, PR3, PR4, PR5, PR6, PR7, PR8, PR11, PR12, PR13, PR14, PR15,      & ! PRX : photolysis rate
!                     PR16, PR17, PR18, PR19, PR21, PR22, PR23, PR24, PR31, PR32, PR33, PR34,    &
!                     PR35, PR41, PR51, PR52, PR53, PR54, PR55, PR56, PR57, PR67, PR68, PR70,    &
!                     PR71,                                                                      &
!                     ACS_o3, ACS_hcho, ACS_ch3coch3, ACS_ch3no3, ACS_c2h5no3, ACS_i_c3h7no3
!
!     REAL(real_x), DIMENSION(75), SAVE :: ACS_mek, ACS_h2o2, ACSA_o3, ACSB_o3, ACS_no2,         & ! ACS_X : absorption cross section
!                                          QY_mek, QY_h2o2, QY_o1d, QY_o3p, QY_no2, QY_hono,     & ! QY_x : quantum yield for photolysis
!                                          ACS_hono, QY_hno3, ACS_hno3, ACS_no3, ACSA_hcho,      &
!                                          QY_no3_no, QY_no3_no2, ACSB_hcho, QY_hcho_h_hco,      &
!                                          QY_hcho_h2_co, ACS_ch3cho, QY_ch3cho, ACS_c2h5cho,    &
!                                          QY_c2h5cho, ACS_ic3h7cho, QY_ic3h7cho, ACS_nc3h7cho,  &
!                                          QY1_nc3h7cho, QY2_nc3h7cho, ACS_macr, QY1_macr,       &
!                                          QY2_macr, ACSA_ch3coch3, ACSB_ch3coch3,               &
!                                          ACSC_ch3coch3, ACSD_ch3coch3, QY_ch3coch3, ACS_mvk,   &
!                                          QY1_mvk, QY2_mvk, ACS_glyox, QY1_glyox, QY3_glyox,    &
!                                          QY2_glyox, ACS_mglyox, QY_mglyox, ACS_biacet,         &
!                                          QY_biacet, ACS_ch3ooh, QY_ch3ooh, ACSA_ch3no3,        &
!                                          ACSB_ch3no3, QY_ch3no3, ACSA_c2h5no3, ACSB_c2h5no3,   &
!                                          QY_c2h5no3, ACS_n_c3h7no3, QY_n_c3h7no3,              &
!                                          QY_i_c3h7no3, ACSA_i_c3h7no3, ACSB_i_c3h7no3,         &
!                                          ACS_t_c4h9no3, QY_t_c4h9no3, ACS_noa, QY1_noa, QY2_noa, &
!                                          ACS_ho2no2, QY1_ho2no2, QY2_ho2no2, ACS_n2o5,         &
!                                          ACSA_n2o5, ACSB_n2o5, QY1_n2o5, QY2_n2o5
!
!     LOGICAL, SAVE :: first_call = .TRUE.
!
!     ! The following is run only once, only on the first time CHEMCALC is called
!     ! Input of absorption cross spectrum and quantum yields for different molecules:
!     IF (first_call) THEN ! do only once
!         first_call = .FALSE.
!
!          OPEN(900, FILE = ''//filename1//'/Photolyse/QY_CS/o3_mcm_v2.dat',           STATUS = 'OLD')
!          OPEN(901, FILE = ''//filename1//'/Photolyse/QY_CS/o1d_mcm_v2.dat',          STATUS = 'OLD')
!          OPEN(902, FILE = ''//filename1//'/Photolyse/QY_CS/o3p_mcm_v2.dat',          STATUS = 'OLD')
!          OPEN(903, FILE = ''//filename1//'/Photolyse/QY_CS/h2o2_mcm_v2.dat',         STATUS = 'OLD')
!          OPEN(904, FILE = ''//filename1//'/Photolyse/QY_CS/no2_mcm_v2.dat',          STATUS = 'OLD')
!          OPEN(905, FILE = ''//filename1//'/Photolyse/QY_CS/no2_qy_mcm_v2.dat',       STATUS = 'OLD')
!          OPEN(906, FILE = ''//filename1//'/Photolyse/QY_CS/no3_mcm_v2.dat',          STATUS = 'OLD')
!          OPEN(907, FILE = ''//filename1//'/Photolyse/QY_CS/no3_no_mcm_v2.dat',       STATUS = 'OLD')
!          OPEN(908, FILE = ''//filename1//'/Photolyse/QY_CS/no3_no2_mcm_v2.dat',      STATUS = 'OLD')
!          OPEN(909, FILE = ''//filename1//'/Photolyse/QY_CS/hono_mcm_v2.dat',         STATUS = 'OLD')
!          OPEN(910, FILE = ''//filename1//'/Photolyse/QY_CS/hno3_mcm_v2_298.dat',     STATUS = 'OLD')
!          OPEN(911, FILE = ''//filename1//'/Photolyse/QY_CS/hcho_mcm_v2.dat',         STATUS = 'OLD')
!          OPEN(912, FILE = ''//filename1//'/Photolyse/QY_CS/hcho_h_hco_mcm_v2.dat',   STATUS = 'OLD')
!          OPEN(913, FILE = ''//filename1//'/Photolyse/QY_CS/hcho_h2_co_mcm_v2.dat',   STATUS = 'OLD')
!          OPEN(914, FILE = ''//filename1//'/Photolyse/QY_CS/ch3cho_mcm_v2.dat',       STATUS = 'OLD')
!          OPEN(915, FILE = ''//filename1//'/Photolyse/QY_CS/ch3cho_qy_mcm_v2.dat',    STATUS = 'OLD')
!          OPEN(916, FILE = ''//filename1//'/Photolyse/QY_CS/c2h5cho_mcm_v2.dat',      STATUS = 'OLD')
!          OPEN(917, FILE = ''//filename1//'/Photolyse/QY_CS/c2h5cho_qy_mcm_v2.dat',   STATUS = 'OLD')
!          OPEN(918, FILE = ''//filename1//'/Photolyse/QY_CS/nc3h7cho_mcm_v2.dat',     STATUS = 'OLD')
!          OPEN(919, FILE = ''//filename1//'/Photolyse/QY_CS/ic3h7cho_mcm_v2.dat',     STATUS = 'OLD')
!          OPEN(920, FILE = ''//filename1//'/Photolyse/QY_CS/ic3h7cho_qy_mcm_v2.dat',  STATUS = 'OLD')
!          OPEN(921, FILE = ''//filename1//'/Photolyse/QY_CS/macr_mcm_v2.dat',         STATUS = 'OLD')
!          OPEN(922, FILE = ''//filename1//'/Photolyse/QY_CS/ch3coch3_mcm_v2.dat',     STATUS = 'OLD')
!          OPEN(923, FILE = ''//filename1//'/Photolyse/QY_CS/ch3coch3_qy_mcm_v2.dat',  STATUS = 'OLD')
!          OPEN(924, FILE = ''//filename1//'/Photolyse/QY_CS/mek_mcm_v2.dat',          STATUS = 'OLD')
!          OPEN(925, FILE = ''//filename1//'/Photolyse/QY_CS/mvk_mcm_v2.dat',          STATUS = 'OLD')
!          OPEN(926, FILE = ''//filename1//'/Photolyse/QY_CS/mvk_qy_mcm_v2.dat',       STATUS = 'OLD')
!          OPEN(927, FILE = ''//filename1//'/Photolyse/QY_CS/glyox_mcm_v2.dat',        STATUS = 'OLD')
!          OPEN(928, FILE = ''//filename1//'/Photolyse/QY_CS/glyox_qy_mcm_v2.dat',     STATUS = 'OLD')
!          OPEN(929, FILE = ''//filename1//'/Photolyse/QY_CS/mglyox_mcm_v2.dat',       STATUS = 'OLD')
!          OPEN(930, FILE = ''//filename1//'/Photolyse/QY_CS/mglyox_qy_mcm_v2.dat',    STATUS = 'OLD')
!          OPEN(931, FILE = ''//filename1//'/Photolyse/QY_CS/biacet_mcm_v2.dat',       STATUS = 'OLD')
!          OPEN(932, FILE = ''//filename1//'/Photolyse/QY_CS/ch3ooh_mcm_v2.dat',       STATUS = 'OLD')
!          OPEN(933, FILE = ''//filename1//'/Photolyse/QY_CS/ch3no3_mcm_v2.dat',       STATUS = 'OLD')
!          OPEN(934, FILE = ''//filename1//'/Photolyse/QY_CS/c2h5no3_mcm_v2.dat',      STATUS = 'OLD')
!          OPEN(935, FILE = ''//filename1//'/Photolyse/QY_CS/n_c3h7no3_mcm_v2.dat',    STATUS = 'OLD')
!          OPEN(936, FILE = ''//filename1//'/Photolyse/QY_CS/i_c3h7no3_mcm_v2.dat',    STATUS = 'OLD')
!          OPEN(937, FILE = ''//filename1//'/Photolyse/QY_CS/tc4h9no3_mcm_v2.dat',     STATUS = 'OLD')
!          OPEN(938, FILE = ''//filename1//'/Photolyse/QY_CS/noa_mcm_v2.dat',          STATUS = 'OLD')
!          OPEN(939, FILE = ''//filename1//'/Photolyse/QY_CS/ho2no2_aitkinson.dat',    STATUS = 'OLD')
!          OPEN(940, FILE = ''//filename1//'/Photolyse/QY_CS/n2o5_aitkinson.dat',      STATUS = 'OLD')
!
!         DO I=1,75
!            READ(900,*)  WL(I), ACSA_o3(I), ACSB_o3(I)
!            READ(901,*)  WL(I), QY_o1d(I)
!            READ(902,*)  WL(I), QY_o3p(I)
!            READ(903,*)  WL(I), ACS_h2o2(I), QY_h2o2(I)
!            READ(904,*)  WL(I), ACS_no2(I)
!            READ(905,*)  WL(I), QY_no2(I)
!            READ(906,*)  WL(I), ACS_no3(I)
!            READ(907,*)  WL(I), QY_no3_no(I)
!            READ(908,*)  WL(I), QY_no3_no2(I)
!            READ(909,*)  WL(I), ACS_hono(I), QY_hono(I)
!            READ(910,*)  WL(I), ACS_hno3(I), QY_hno3(I)
!            READ(911,*)  WL(I), ACSA_hcho(I), ACSB_hcho(I)
!            READ(912,*)  WL(I), QY_hcho_h_hco(I)
!            READ(913,*)  WL(I), QY_hcho_h2_co(I)
!            READ(914,*)  WL(I), ACS_ch3cho(I)
!            READ(915,*)  WL(I), QY_ch3cho(I)
!            READ(916,*)  WL(I), ACS_c2h5cho(I)
!            READ(917,*)  WL(I), QY_c2h5cho(I)
!            READ(918,*)  WL(I), ACS_nc3h7cho(I), QY1_nc3h7cho(I), QY2_nc3h7cho(I)
!            READ(919,*)  WL(I), ACS_ic3h7cho(I)
!            READ(920,*)  WL(I), QY_ic3h7cho(I)
!            READ(921,*)  WL(I), ACS_macr(I), QY1_macr(I), QY2_macr(I)
!            READ(922,*)  WL(I), ACSA_ch3coch3(I), ACSB_ch3coch3(I), ACSC_ch3coch3(I), ACSD_ch3coch3(I)
!            READ(923,*)  WL(I), QY_ch3coch3(I)
!            READ(924,*)  WL(I), ACS_mek(I), QY_mek(I)
!            READ(925,*)  WL(I), ACS_mvk(I)
!            READ(926,*)  WL(I), QY1_mvk(I), QY2_mvk(I)
!            READ(927,*)  WL(I), ACS_glyox(I)
!            READ(928,*)  WL(I), QY1_glyox(I), QY3_glyox(I), QY2_glyox(I)
!            READ(929,*)  WL(I), ACS_mglyox(I)
!            READ(930,*)  WL(I), QY_mglyox(I)
!            READ(931,*)  WL(I), ACS_biacet(I), QY_biacet(I)
!            READ(932,*)  WL(I), ACS_ch3ooh(I), QY_ch3ooh(I)
!            READ(933,*)  WL(I), ACSA_ch3no3(I), ACSB_ch3no3(I), QY_ch3no3(I)
!            READ(934,*)  WL(I), ACSA_c2h5no3(I), ACSB_c2h5no3(I), QY_c2h5no3(I)
!            READ(935,*)  WL(I), ACS_n_c3h7no3(I), QY_n_c3h7no3(I)
!            READ(936,*)  WL(I), ACSA_i_c3h7no3(I), ACSB_i_c3h7no3(I), QY_i_c3h7no3(I)
!            READ(937,*)  WL(I), ACS_t_c4h9no3(I), QY_t_c4h9no3(I)
!            READ(938,*)  WL(I), ACS_noa(I), QY1_noa(I), QY2_noa(I)
!            READ(939,*)  WL(I), ACS_ho2no2(I), QY1_ho2no2(I), QY2_ho2no2(I)
!            READ(940,*)  WL(I), ACSA_n2o5(I), ACSB_n2o5(I), QY1_n2o5(I), QY2_n2o5(I)
!         ENDDO
!
!         CLOSE(900)
!         CLOSE(901)
!         CLOSE(902)
!         CLOSE(903)
!         CLOSE(904)
!         CLOSE(905)
!         CLOSE(906)
!         CLOSE(907)
!         CLOSE(908)
!         CLOSE(909)
!         CLOSE(910)
!         CLOSE(911)
!         CLOSE(912)
!         CLOSE(913)
!         CLOSE(914)
!         CLOSE(915)
!         CLOSE(916)
!         CLOSE(917)
!         CLOSE(918)
!         CLOSE(919)
!         CLOSE(920)
!         CLOSE(921)
!         CLOSE(922)
!         CLOSE(923)
!         CLOSE(924)
!         CLOSE(925)
!         CLOSE(926)
!         CLOSE(927)
!         CLOSE(928)
!         CLOSE(929)
!         CLOSE(930)
!         CLOSE(931)
!         CLOSE(932)
!         CLOSE(933)
!         CLOSE(934)
!         CLOSE(935)
!         CLOSE(936)
!         CLOSE(937)
!         CLOSE(938)
!         CLOSE(939)
!         CLOSE(940)
!      ENDIF ! do only once
!
!     ! Put all J-values to 0 when new time step begins:
!     PR1 = 0.
!     PR2 = 0.
!     PR3 = 0.
!     PR4 = 0.
!     PR5 = 0.
!     PR6 = 0.
!     PR7 = 0.
!     PR8 = 0.
!     PR11 = 0.
!     PR12 = 0.
!     PR13 = 0.
!     PR14 = 0.
!     PR15 = 0.
!     PR16 = 0.
!     PR17 = 0.
!     PR18 = 0.
!     PR19 = 0.
!     PR21 = 0.
!     PR22 = 0.
!     PR23 = 0.
!     PR24 = 0.
!     PR31 = 0.
!     PR32 = 0.
!     PR33 = 0.
!     PR34 = 0.
!     PR35 = 0.
!     PR41 = 0.
!     PR51 = 0.
!     PR52 = 0.
!     PR53 = 0.
!     PR54 = 0.
!     PR55 = 0.
!     PR56 = 0.
!     PR57 = 0.
!     PR67 = 0.
!     PR68 = 0.
!     PR70 = 0.
!     PR71 = 0.
!
!     ! Delta lambda and dimension parameters are given:
!     DL  = 5.
!     VA  = 1.E14
!     VAA = 1.E-20
!     VAB = 1.E-19
!     VAC = 1.E-21
!     VAD = 1.E-24
!
!     ! Calculate J-values:
!       DO I=1,75
!          ACS_o3 = ACSA_o3(I) * EXP(ACSB_o3(I)/T)
!          PR1 = PR1 + ACS_o3 * QY_o1d(I) * SAFN((WL(I)-280)/5+1) * DL * VA
!          PR2 = PR2 + ACS_o3 * QY_o3p(I) * SAFN((WL(I)-280)/5+1) * DL * VA
!          PR3 = PR3 + ACS_h2o2(I) * QY_h2o2(I) * SAFN((WL(I)-280)/5+1) * DL * VA
!          PR4 = PR4 + ACS_no2(I) * QY_no2(I) * SAFN((WL(I)-280)/5+1) * DL * VA * VAA
!          PR5 = PR5 + ACS_no3(I) * VAB * QY_no3_no(I) *SAFN((WL(I)-280)/5+1) * DL * VA
!          PR6 = PR6 + ACS_no3(I) * VAB * QY_no3_no2(I) *SAFN((WL(I)-280)/5+1) * DL * VA
!          PR7 = PR7 + ACS_hono(I) * QY_hono(I) *SAFN((WL(I)-280)/5+1) * DL * VA
!          PR8 = PR8 + ACS_hno3(I) * VAA * QY_hno3(I) * SAFN((WL(I)-280)/5+1) * DL * VA
!          ACS_hcho = (ACSA_hcho(I) * VAC) + (ACSB_hcho(I) * VAD * (T-298))
!          PR11 = PR11 + ACS_hcho * QY_hcho_h_hco(I) * SAFN((WL(I)-280)/5+1) * DL * VA
!          PR12 = PR12 + ACS_hcho * QY_hcho_h2_co(I) * SAFN((WL(I)-280)/5+1) * DL * VA
!          PR13 = PR13 + ACS_ch3cho(I) * QY_ch3cho(I) * SAFN((WL(I)-280)/5+1) * DL * VA
!          PR14 = PR14 + ACS_c2h5cho(I) * QY_c2h5cho(I) * SAFN((WL(I)-280)/5+1) * DL * VA
!          PR15 = PR15 + ACS_nc3h7cho(I) * VAC * QY1_nc3h7cho(I) * SAFN((WL(I)-280)/5+1) * DL * VA
!          PR16 = PR16 + ACS_nc3h7cho(I) * VAC * QY2_nc3h7cho(I) * SAFN((WL(I)-280)/5+1) * DL * VA
!          PR17 = PR17 + ACS_ic3h7cho(I) * VAC * QY_ic3h7cho(I) * SAFN((WL(I)-280)/5+1) * DL * VA
!          PR18 = PR18 + ACS_macr(I) * QY1_macr(I) * SAFN((WL(I)-280)/5+1) * DL * VA
!          PR19 = PR19 + ACS_macr(I) * QY2_macr(I) * SAFN((WL(I)-280)/5+1) * DL * VA
!          ACS_ch3coch3 = ACSA_ch3coch3(I) * (1 + (ACSB_ch3coch3(I) * T) + (ACSC_ch3coch3(I) * (T**2)) + (ACSD_ch3coch3(I) * (T**3)))
!          PR21 = PR21 + ACS_ch3coch3 * QY_ch3coch3(I) * SAFN((WL(I)-280)/5+1) * DL * VA
!          PR22 = PR22 + ACS_mek(I) * QY_mek(I) * SAFN((WL(I)-280)/5+1) * DL * VA
!          PR23 = PR23 + ACS_mvk(I) * QY1_mvk(I) * SAFN((WL(I)-280)/5+1) * DL * VA
!          PR24 = PR24 + ACS_mvk(I) * QY2_mvk(I) * SAFN((WL(I)-280)/5+1) * DL * VA
!          PR31 = PR31 + ACS_glyox(I) * VAA * QY1_glyox(I) * SAFN((WL(I)-280)/5+1) * DL * VA
!          PR32 = PR32 + ACS_glyox(I) * VAA * QY2_glyox(I) * SAFN((WL(I)-280)/5+1) * DL * VA
!          PR33 = PR33 + ACS_glyox(I) * VAA * QY3_glyox(I) * SAFN((WL(I)-280)/5+1) * DL * VA
!          PR34 = PR34 + ACS_mglyox(I) * VAA * QY_mglyox(I) * SAFN((WL(I)-280)/5+1) * DL * VA
!          PR35 = PR35 + ACS_biacet(I) * VAA * QY_biacet(I) * SAFN((WL(I)-280)/5+1) * DL * VA
!          PR41 = PR41 + ACS_ch3ooh(I) * VAA * QY_ch3ooh(I) * SAFN((WL(I)-280)/5+1) * DL * VA
!          ACS_ch3no3 = ACSA_ch3no3(I) * VAA * EXP(ACSB_ch3no3(I) * 1.E-3 * (T-298))
!          PR51 = PR51 + ACS_ch3no3 * QY_ch3no3(I) * SAFN((WL(I)-280)/5+1) * DL * VA
!          ACS_c2h5no3 = ACSA_c2h5no3(I) * EXP(ACSB_c2h5no3(I) * (T-298))
!          PR52 = PR52 + ACS_c2h5no3 * QY_c2h5no3(I) * SAFN((WL(I)-280)/5+1) * DL * VA
!          PR53 = PR53 + ACS_n_c3h7no3(I) * QY_n_c3h7no3(I) * SAFN((WL(I)-280)/5+1) * DL * VA
!          ACS_i_c3h7no3 = ACSA_i_c3h7no3(I) * EXP(ACSB_i_c3h7no3(I) * (T-298))
!          PR54 = PR54 + ACS_i_c3h7no3 * QY_i_c3h7no3(I) * SAFN((WL(I)-280)/5+1) * DL * VA
!          PR55 = PR55 + ACS_t_c4h9no3(I) * QY_t_c4h9no3(I) * SAFN((WL(I)-280)/5+1) * DL * VA
!          PR56 = PR56 + ACS_noa(I) * QY1_noa(I) * SAFN((WL(I)-280)/5+1) * DL * VA
!          PR57 = PR57 + ACS_noa(I) * QY2_noa(I) * SAFN((WL(I)-280)/5+1) * DL * VA
!          PR67 = PR67 + ACS_ho2no2(I) * QY1_ho2no2(I) * SAFN((WL(I)-280)/5+1) * DL * VA * VAA
!          PR68 = PR68 + ACS_ho2no2(I) * QY2_ho2no2(I) * SAFN((WL(I)-280)/5+1) * DL * VA * VAA
!          ACS_n2o5 = ACSA_n2o5(I) * 10**(1000*ACSB_n2o5(I)*((1./T)-(1./298.)))
!          PR70 = PR70 + ACS_n2o5(I) * QY1_n2o5(I)  * SAFN((WL(I)-280)/5+1) * DL * VA * VAA
!          PR71 = PR71 + ACS_n2o5(I) * QY2_n2o5(I)  * SAFN((WL(I)-280)/5+1) * DL * VA * VAA
!       ENDDO
!
!
!     END SUBROUTINE
!
!     END MODULE
