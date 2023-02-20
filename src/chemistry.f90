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

    LOGICAL,  INTENT(in)    :: swr_is_af
    REAL(dp), INTENT(out)   :: CONS(NSPEC)               ! Concentrations of all gas compounds
    REAL(dp), INTENT(out)   :: RO2_out                   ! Concentration of peroxy radical
    REAL(dp), INTENT(in)    :: TIME_kpp, END_kpp         ! Model time in sec and same time plus chemical timestep
    REAL(dp), INTENT(in)    :: H2O                       ! Water vapour concentration
    REAL(dp), INTENT(in)    :: AIR                       ! air molecules - calculated in the main file
    REAL(dp), INTENT(in)    :: Beta                      ! Solar elevation angle
    REAL(dp), INTENT(in)    :: Tp                        ! Temperature [K]
    REAL(dp), INTENT(in)    :: SWDWN                     ! Global short wave radiation
    REAL(dp), INTENT(in)    :: RES1, RES2                ! Condensation sink for H2SO4 and HNO3
    REAL(dp), INTENT(in)    :: Albedo                    ! Albedo
    REAL(dp), INTENT(in)    :: swr_dis(:)                ! Spectral data - either distribution or absolute values
    REAL(DP), INTENT(out)   :: reactivities(:)
    REAL(dp)                :: I_ACF(75)                 ! Calculated actinic flux from 280-650nm
    REAL(dp)                :: O2                        ! O2-concentration in molecule / cm3
    REAL(dp)                :: H2                        ! H2-concentration in molecule / cm3
    REAL(dp)                :: N2                        ! N2-concentration in molecule / cm3
    REAL(dp)                :: MO2N2                     ! MO2N2: Sum of O2 and N2 concentration
    REAL(dp)                :: ZSD                       ! Solar zenith angle
    INTEGER                 :: I

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
    REAL(dp), DIMENSION(56,75,7), SAVE :: CS_QY=0d0 ! cross section parameters. First rank indices refer to reaction (J1, J2 ...)
    REAL(dp), INTENT(IN)               :: SAFN(:)   ! spectral actinic flux
    REAL(dp), INTENT(INOUT)            :: J_pc(:)   ! Photochemistry rates (some are empty)
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
