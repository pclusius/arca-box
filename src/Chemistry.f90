!ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
!
!   CHEMISTRY MODULE
!
!   Calculations of all J-values and link to KPP-files
!
!   This version is for MCM version 3.3.1 coupled to PRAM
!
!ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

MODULE CHEMISTRY

    USE CONSTANTS
    USE SECOND_PRECISION,  ONLY : DP           ! KPP Numerical type
    USE SECOND_PARAMETERS, ONLY : NSPEC        ! Number of chemical species
    USE SECOND_GLOBAL,     ONLY : NPHOT        ! Number of photochemical reactions used in KPP

    IMPLICIT NONE

    PRIVATE

    PUBLIC :: NKVALUES
    PUBLIC :: current_j_values

    REAL(real_x), DIMENSION(57) :: current_j_values

    INTEGER, PARAMETER :: NKVALUES = 42

CONTAINS

    !ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
    !
    !   Subroutine for the photodissociation reactions
    !
    !ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

    SUBROUTINE PHOTOLYSIS_KPP(TS,   Beta, SAFN, T,                             &
        PR1,  PR2,  PR3,  PR4,  PR5,  PR6,  PR7,  PR8,  PR11, PR12, PR13,      &
        PR14, PR15, PR16, PR17, PR18, PR19, PR21, PR22, PR23, PR24, PR31,      &
        PR32, PR33, PR34, PR35, PR41, PR51, PR52, PR53, PR54, PR55, PR56,      &
        PR57, PR67, PR68, PR70, PR71)

        IMPLICIT NONE

        INTEGER I
  
        INTEGER, DIMENSION(75), SAVE :: WL

        REAL(DP), DIMENSION(75) :: SAFN ! spectral actinic flux - based on hyytiala data

        REAL(DP) :: TS, T,                                                             &
            Beta,                                                                      & ! Zenit angle
            DL,                                                                        & ! Delta lambda
            VA, VAA, VAB, VAC, VAD,                                                    & ! Dimension parameter
            PR1, PR2, PR3, PR4, PR5, PR6, PR7, PR8, PR11, PR12, PR13, PR14, PR15,      & ! PRX: photolysis rate
            PR16, PR17, PR18, PR19, PR21, PR22, PR23, PR24, PR31, PR32, PR33, PR34,    &
            PR35, PR41, PR51, PR52, PR53, PR54, PR55, PR56, PR57, PR67, PR68, PR70,    &
            PR71,                                                                      &
            ACS_o3, ACS_hcho, ACS_ch3coch3, ACS_ch3no3, ACS_c2h5no3, ACS_i_c3h7no3
 
        REAL(DP), DIMENSION(75), SAVE :: ACS_mek, ACS_h2o2,       & ! ACS_X : absorption cross sections
            ACSA_o3, ACSB_o3, ACS_no2,                            & ! QY_x : quantum yields for photolysis
            QY_mek, QY_h2o2, QY_o1d, QY_o3p, QY_no2, QY_hono,     &
            ACS_hono, QY_hno3, ACS_hno3, ACS_no3, ACSA_hcho,      &
            QY_no3_no, QY_no3_no2, ACSB_hcho, QY_hcho_h_hco,      &
            QY_hcho_h2_co, ACS_ch3cho, QY_ch3cho, ACS_c2h5cho,    &
            QY_c2h5cho, ACS_ic3h7cho, QY_ic3h7cho, ACS_nc3h7cho,  &
            QY1_nc3h7cho, QY2_nc3h7cho, ACS_macr, QY1_macr,       &
            QY2_macr, ACSA_ch3coch3, ACSB_ch3coch3,               &
            ACSC_ch3coch3, ACSD_ch3coch3, QY_ch3coch3, ACS_mvk,   &
            QY1_mvk, QY2_mvk, ACS_glyox, QY1_glyox, QY3_glyox,    &
            QY2_glyox, ACS_mglyox, QY_mglyox, ACS_biacet,         &
            QY_biacet, ACS_ch3ooh, QY_ch3ooh, ACSA_ch3no3,        &
            ACSB_ch3no3, QY_ch3no3, ACSA_c2h5no3, ACSB_c2h5no3,   &
            QY_c2h5no3, ACS_n_c3h7no3, QY_n_c3h7no3,              &
            QY_i_c3h7no3, ACSA_i_c3h7no3, ACSB_i_c3h7no3,         &
            ACS_t_c4h9no3, QY_t_c4h9no3, ACS_noa, QY1_noa,        &
            QY2_noa, ACS_ho2no2, QY1_ho2no2, QY2_ho2no2,          &
            ACS_n2o5, ACSA_n2o5, ACSB_n2o5, QY1_n2o5, QY2_n2o5

        LOGICAL, SAVE :: first_call = .TRUE.

        ! The following is run only once, only on the first time CHEMCALC is called
        ! Input of absorption cross spectrum and quantum yields for different molecules:
        IF (first_call) THEN ! do only once
            first_call = .FALSE.

            OPEN(900,  FILE = ''//filename1//'/General/Photolyse/o3_mcm_v2.dat',       STATUS = 'OLD')
            OPEN(901,  FILE = ''//filename1//'/General/Photolyse/o1d_mcm_v2.dat',       STATUS = 'OLD')
            OPEN(902,  FILE = ''//filename1//'/General/Photolyse/o3p_mcm_v2.dat',       STATUS = 'OLD')
            OPEN(903,  FILE = ''//filename1//'/General/Photolyse/h2o2_mcm_v2.dat',       STATUS = 'OLD')
            OPEN(904,  FILE = ''//filename1//'/General/Photolyse/no2_mcm_v2.dat',       STATUS = 'OLD')
            OPEN(905,  FILE = ''//filename1//'/General/Photolyse/no2_qy_mcm_v2.dat',       STATUS = 'OLD')
            OPEN(906,  FILE = ''//filename1//'/General/Photolyse/no3_mcm_v2.dat',       STATUS = 'OLD')
            OPEN(907,  FILE = ''//filename1//'/General/Photolyse/no3_no_mcm_v2.dat',       STATUS = 'OLD')
            OPEN(908,  FILE = ''//filename1//'/General/Photolyse/no3_no2_mcm_v2.dat',       STATUS = 'OLD')
            OPEN(909,  FILE = ''//filename1//'/General/Photolyse/hono_mcm_v2.dat',       STATUS = 'OLD')
            OPEN(910,  FILE = ''//filename1//'/General/Photolyse/hno3_mcm_v2_298.dat',       STATUS = 'OLD')
            OPEN(911,  FILE = ''//filename1//'/General/Photolyse/hcho_mcm_v2.dat',       STATUS = 'OLD')
            OPEN(912,  FILE = ''//filename1//'/General/Photolyse/hcho_h_hco_mcm_v2.dat',       STATUS = 'OLD')
            OPEN(913,  FILE = ''//filename1//'/General/Photolyse/hcho_h2_co_mcm_v2.dat',       STATUS = 'OLD')
            OPEN(914,  FILE = ''//filename1//'/General/Photolyse/ch3cho_mcm_v2.dat',       STATUS = 'OLD')
            OPEN(915,  FILE = ''//filename1//'/General/Photolyse/ch3cho_qy_mcm_v2.dat',       STATUS = 'OLD')
            OPEN(916,  FILE = ''//filename1//'/General/Photolyse/c2h5cho_mcm_v2.dat',       STATUS = 'OLD')
            OPEN(917,  FILE = ''//filename1//'/General/Photolyse/c2h5cho_qy_mcm_v2.dat',       STATUS = 'OLD')
            OPEN(918,  FILE = ''//filename1//'/General/Photolyse/nc3h7cho_mcm_v2.dat',       STATUS = 'OLD')
            OPEN(919,  FILE = ''//filename1//'/General/Photolyse/ic3h7cho_mcm_v2.dat',       STATUS = 'OLD')
            OPEN(920,  FILE = ''//filename1//'/General/Photolyse/ic3h7cho_qy_mcm_v2.dat',       STATUS = 'OLD')
            OPEN(921,  FILE = ''//filename1//'/General/Photolyse/macr_mcm_v2.dat',       STATUS = 'OLD')
            OPEN(922,  FILE = ''//filename1//'/General/Photolyse/ch3coch3_mcm_v2.dat',       STATUS = 'OLD')
            OPEN(923,  FILE = ''//filename1//'/General/Photolyse/ch3coch3_qy_mcm_v2.dat',       STATUS = 'OLD')
            OPEN(924,  FILE = ''//filename1//'/General/Photolyse/mek_mcm_v2.dat',       STATUS = 'OLD')
            OPEN(925,  FILE = ''//filename1//'/General/Photolyse/mvk_mcm_v2.dat',       STATUS = 'OLD')
            OPEN(926,  FILE = ''//filename1//'/General/Photolyse/mvk_qy_mcm_v2.dat',       STATUS = 'OLD')
            OPEN(927,  FILE = ''//filename1//'/General/Photolyse/glyox_mcm_v2.dat',       STATUS = 'OLD')
            OPEN(928,  FILE = ''//filename1//'/General/Photolyse/glyox_qy_mcm_v2.dat',       STATUS = 'OLD')
            OPEN(929,  FILE = ''//filename1//'/General/Photolyse/mglyox_mcm_v2.dat',       STATUS = 'OLD')
            OPEN(930,  FILE = ''//filename1//'/General/Photolyse/mglyox_qy_mcm_v2.dat',       STATUS = 'OLD')
            OPEN(931,  FILE = ''//filename1//'/General/Photolyse/biacet_mcm_v2.dat',       STATUS = 'OLD')
            OPEN(932,  FILE = ''//filename1//'/General/Photolyse/ch3ooh_mcm_v2.dat',       STATUS = 'OLD')
            OPEN(933,  FILE = ''//filename1//'/General/Photolyse/ch3no3_mcm_v2.dat',       STATUS = 'OLD')
            OPEN(934,  FILE = ''//filename1//'/General/Photolyse/c2h5no3_mcm_v2.dat',       STATUS = 'OLD')
            OPEN(935,  FILE = ''//filename1//'/General/Photolyse/n_c3h7no3_mcm_v2.dat',       STATUS = 'OLD')
            OPEN(936,  FILE = ''//filename1//'/General/Photolyse/i_c3h7no3_mcm_v2.dat',       STATUS = 'OLD')
            OPEN(937,  FILE = ''//filename1//'/General/Photolyse/tc4h9no3_mcm_v2.dat',       STATUS = 'OLD')
            OPEN(938,  FILE = ''//filename1//'/General/Photolyse/noa_mcm_v2.dat',       STATUS = 'OLD')
            OPEN(939,  FILE = ''//filename1//'/General/Photolyse/ho2no2_aitkinson.dat',       STATUS = 'OLD')
            OPEN(940,  FILE = ''//filename1//'/General/Photolyse/n2o5_aitkinson.dat',       STATUS = 'OLD')

            DO I=1,75
                READ(900,*)  WL(I), ACSA_o3(I), ACSB_o3(I)
                READ(901,*)  WL(I), QY_o1d(I)
                READ(902,*)  WL(I), QY_o3p(I)
                READ(903,*)  WL(I), ACS_h2o2(I), QY_h2o2(I)
                READ(904,*)  WL(I), ACS_no2(I)
                READ(905,*)  WL(I), QY_no2(I)
                READ(906,*)  WL(I), ACS_no3(I)
                READ(907,*)  WL(I), QY_no3_no(I)
                READ(908,*)  WL(I), QY_no3_no2(I)
                READ(909,*)  WL(I), ACS_hono(I), QY_hono(I)
                READ(910,*)  WL(I), ACS_hno3(I), QY_hno3(I)
                READ(911,*)  WL(I), ACSA_hcho(I), ACSB_hcho(I)
                READ(912,*)  WL(I), QY_hcho_h_hco(I)
                READ(913,*)  WL(I), QY_hcho_h2_co(I)
                READ(914,*)  WL(I), ACS_ch3cho(I)
                READ(915,*)  WL(I), QY_ch3cho(I)
                READ(916,*)  WL(I), ACS_c2h5cho(I)
                READ(917,*)  WL(I), QY_c2h5cho(I)
                READ(918,*)  WL(I), ACS_nc3h7cho(I), QY1_nc3h7cho(I), QY2_nc3h7cho(I)
                READ(919,*)  WL(I), ACS_ic3h7cho(I)
                READ(920,*)  WL(I), QY_ic3h7cho(I)
                READ(921,*)  WL(I), ACS_macr(I), QY1_macr(I), QY2_macr(I)
                READ(922,*)  WL(I), ACSA_ch3coch3(I), ACSB_ch3coch3(I), ACSC_ch3coch3(I), ACSD_ch3coch3(I)
                READ(923,*)  WL(I), QY_ch3coch3(I)
                READ(924,*)  WL(I), ACS_mek(I), QY_mek(I)
                READ(925,*)  WL(I), ACS_mvk(I)
                READ(926,*)  WL(I), QY1_mvk(I), QY2_mvk(I)
                READ(927,*)  WL(I), ACS_glyox(I)
                READ(928,*)  WL(I), QY1_glyox(I), QY3_glyox(I), QY2_glyox(I)
                READ(929,*)  WL(I), ACS_mglyox(I)
                READ(930,*)  WL(I), QY_mglyox(I)
                READ(931,*)  WL(I), ACS_biacet(I), QY_biacet(I)
                READ(932,*)  WL(I), ACS_ch3ooh(I), QY_ch3ooh(I)
                READ(933,*)  WL(I), ACSA_ch3no3(I), ACSB_ch3no3(I), QY_ch3no3(I)
                READ(934,*)  WL(I), ACSA_c2h5no3(I), ACSB_c2h5no3(I), QY_c2h5no3(I)
                READ(935,*)  WL(I), ACS_n_c3h7no3(I), QY_n_c3h7no3(I)
                READ(936,*)  WL(I), ACSA_i_c3h7no3(I), ACSB_i_c3h7no3(I), QY_i_c3h7no3(I)
                READ(937,*)  WL(I), ACS_t_c4h9no3(I), QY_t_c4h9no3(I)
                READ(938,*)  WL(I), ACS_noa(I), QY1_noa(I), QY2_noa(I)
                READ(939,*)  WL(I), ACS_ho2no2(I), QY1_ho2no2(I), QY2_ho2no2(I)
                READ(940,*)  WL(I), ACSA_n2o5(I), ACSB_n2o5(I), QY1_n2o5(I), QY2_n2o5(I)
            ENDDO
  
            CLOSE(900)
            CLOSE(901)
            CLOSE(902)
            CLOSE(903)
            CLOSE(904)
            CLOSE(905)
            CLOSE(906)
            CLOSE(907)
            CLOSE(908)
            CLOSE(909)
            CLOSE(910)
            CLOSE(911)
            CLOSE(912)
            CLOSE(913)
            CLOSE(914)
            CLOSE(915)
            CLOSE(916)
            CLOSE(917)
            CLOSE(918)
            CLOSE(919)
            CLOSE(920)
            CLOSE(921)
            CLOSE(922)
            CLOSE(923)
            CLOSE(924)
            CLOSE(925)
            CLOSE(926)
            CLOSE(927)
            CLOSE(928)
            CLOSE(929)
            CLOSE(930)
            CLOSE(931)
            CLOSE(932)
            CLOSE(933)
            CLOSE(934)
            CLOSE(935)
            CLOSE(936)
            CLOSE(937)
            CLOSE(938)
            CLOSE(939)
            CLOSE(940)
        ENDIF ! do only once

        ! Put all J-values to 0 when new time step begins:
        PR1 = 0.
        PR2 = 0.
        PR3 = 0.
        PR4 = 0.
        PR5 = 0.
        PR6 = 0.
        PR7 = 0.
        PR8 = 0.
        PR11 = 0.
        PR12 = 0.
        PR13 = 0.
        PR14 = 0.
        PR15 = 0.
        PR16 = 0.
        PR17 = 0.
        PR18 = 0.
        PR19 = 0.
        PR21 = 0.
        PR22 = 0.
        PR23 = 0.
        PR24 = 0.
        PR31 = 0.
        PR32 = 0.
        PR33 = 0.
        PR34 = 0.
        PR35 = 0.
        PR41 = 0.
        PR51 = 0.
        PR52 = 0.
        PR53 = 0.
        PR54 = 0.
        PR55 = 0.
        PR56 = 0.
        PR57 = 0.
        PR67 = 0.
        PR68 = 0.
        PR70 = 0.
        PR71 = 0.
    
        ! Delta lambda and dimension parameters are given:
        DL  = 5.
        VA  = 1.E14
        VAA = 1.E-20
        VAB = 1.E-19
        VAC = 1.E-21
        VAD = 1.E-24
   
        ! Calculate J-values:
        DO I=1,75
            ACS_o3 = ACSA_o3(I) * EXP(ACSB_o3(I)/T)
            PR1 = PR1 + ACS_o3 * QY_o1d(I) * SAFN((WL(I)-280)/5+1) * DL * VA
            PR2 = PR2 + ACS_o3 * QY_o3p(I) * SAFN((WL(I)-280)/5+1) * DL * VA
            PR3 = PR3 + ACS_h2o2(I) * QY_h2o2(I) * SAFN((WL(I)-280)/5+1) * DL * VA
            PR4 = PR4 + ACS_no2(I) * QY_no2(I) * SAFN((WL(I)-280)/5+1) * DL * VA * VAA
            PR5 = PR5 + ACS_no3(I) * VAB * QY_no3_no(I) *SAFN((WL(I)-280)/5+1) * DL * VA
            PR6 = PR6 + ACS_no3(I) * VAB * QY_no3_no2(I) *SAFN((WL(I)-280)/5+1) * DL * VA
            PR7 = PR7 + ACS_hono(I) * QY_hono(I) *SAFN((WL(I)-280)/5+1) * DL * VA
            PR8 = PR8 + ACS_hno3(I) * VAA * QY_hno3(I) * SAFN((WL(I)-280)/5+1) * DL * VA
            ACS_hcho = (ACSA_hcho(I) * VAC) + (ACSB_hcho(I) * VAD * (T-298))
            PR11 = PR11 + ACS_hcho * QY_hcho_h_hco(I) * SAFN((WL(I)-280)/5+1) * DL * VA
            PR12 = PR12 + ACS_hcho * QY_hcho_h2_co(I) * SAFN((WL(I)-280)/5+1) * DL * VA
            PR13 = PR13 + ACS_ch3cho(I) * QY_ch3cho(I) * SAFN((WL(I)-280)/5+1) * DL * VA
            PR14 = PR14 + ACS_c2h5cho(I) * QY_c2h5cho(I) * SAFN((WL(I)-280)/5+1) * DL * VA
            PR15 = PR15 + ACS_nc3h7cho(I) * VAC * QY1_nc3h7cho(I) * SAFN((WL(I)-280)/5+1) * DL * VA
            PR16 = PR16 + ACS_nc3h7cho(I) * VAC * QY2_nc3h7cho(I) * SAFN((WL(I)-280)/5+1) * DL * VA
            PR17 = PR17 + ACS_ic3h7cho(I) * VAC * QY_ic3h7cho(I) * SAFN((WL(I)-280)/5+1) * DL * VA
            PR18 = PR18 + ACS_macr(I) * QY1_macr(I) * SAFN((WL(I)-280)/5+1) * DL * VA
            PR19 = PR19 + ACS_macr(I) * QY2_macr(I) * SAFN((WL(I)-280)/5+1) * DL * VA
            ACS_ch3coch3 = ACSA_ch3coch3(I) * (1 + (ACSB_ch3coch3(I) * T) + (ACSC_ch3coch3(I) * (T**2)) + (ACSD_ch3coch3(I) * (T**3)))
            PR21 = PR21 + ACS_ch3coch3 * QY_ch3coch3(I) * SAFN((WL(I)-280)/5+1) * DL * VA
            PR22 = PR22 + ACS_mek(I) * QY_mek(I) * SAFN((WL(I)-280)/5+1) * DL * VA
            PR23 = PR23 + ACS_mvk(I) * QY1_mvk(I) * SAFN((WL(I)-280)/5+1) * DL * VA
            PR24 = PR24 + ACS_mvk(I) * QY2_mvk(I) * SAFN((WL(I)-280)/5+1) * DL * VA
            PR31 = PR31 + ACS_glyox(I) * VAA * QY1_glyox(I) * SAFN((WL(I)-280)/5+1) * DL * VA
            PR32 = PR32 + ACS_glyox(I) * VAA * QY2_glyox(I) * SAFN((WL(I)-280)/5+1) * DL * VA
            PR33 = PR33 + ACS_glyox(I) * VAA * QY3_glyox(I) * SAFN((WL(I)-280)/5+1) * DL * VA
            PR34 = PR34 + ACS_mglyox(I) * VAA * QY_mglyox(I) * SAFN((WL(I)-280)/5+1) * DL * VA
            PR35 = PR35 + ACS_biacet(I) * VAA * QY_biacet(I) * SAFN((WL(I)-280)/5+1) * DL * VA
            PR41 = PR41 + ACS_ch3ooh(I) * VAA * QY_ch3ooh(I) * SAFN((WL(I)-280)/5+1) * DL * VA
            ACS_ch3no3 = ACSA_ch3no3(I) * VAA * EXP(ACSB_ch3no3(I) * 1.E-3 * (T-298))
            PR51 = PR51 + ACS_ch3no3 * QY_ch3no3(I) * SAFN((WL(I)-280)/5+1) * DL * VA
            ACS_c2h5no3 = ACSA_c2h5no3(I) * EXP(ACSB_c2h5no3(I) * (T-298))
            PR52 = PR52 + ACS_c2h5no3 * QY_c2h5no3(I) * SAFN((WL(I)-280)/5+1) * DL * VA
            PR53 = PR53 + ACS_n_c3h7no3(I) * QY_n_c3h7no3(I) * SAFN((WL(I)-280)/5+1) * DL * VA
            ACS_i_c3h7no3 = ACSA_i_c3h7no3(I) * EXP(ACSB_i_c3h7no3(I) * (T-298))
            PR54 = PR54 + ACS_i_c3h7no3 * QY_i_c3h7no3(I) * SAFN((WL(I)-280)/5+1) * DL * VA
            PR55 = PR55 + ACS_t_c4h9no3(I) * QY_t_c4h9no3(I) * SAFN((WL(I)-280)/5+1) * DL * VA
            PR56 = PR56 + ACS_noa(I) * QY1_noa(I) * SAFN((WL(I)-280)/5+1) * DL * VA
            PR57 = PR57 + ACS_noa(I) * QY2_noa(I) * SAFN((WL(I)-280)/5+1) * DL * VA
            PR67 = PR67 + ACS_ho2no2(I) * QY1_ho2no2(I) * SAFN((WL(I)-280)/5+1) * DL * VA * VAA
            PR68 = PR68 + ACS_ho2no2(I) * QY2_ho2no2(I) * SAFN((WL(I)-280)/5+1) * DL * VA * VAA
            ACS_n2o5 = ACSA_n2o5(I) * 10**(1000*ACSB_n2o5(I)*((1/T)-(1/298)))
            PR70 = PR70 + ACS_n2o5(I) * QY1_n2o5(I)  * SAFN((WL(I)-280)/5+1) * DL * VA * VAA
            PR71 = PR71 + ACS_n2o5(I) * QY2_n2o5(I)  * SAFN((WL(I)-280)/5+1) * DL * VA * VAA
        ENDDO

    END SUBROUTINE

    !ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
    !
    !   RATECONS
    !
    !   Calculation of simple and complex rate coefficients
    !
    !ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

    SUBROUTINE RATECONS(Temp, O2, N2, MO2N2, H2O, NKVALUES, RCON, DT_kpp)

        IMPLICIT NONE

        INTEGER, PARAMETER :: real_x = SELECTED_REAL_KIND(p=14, r=300)

        INTEGER :: NKVALUES

        REAL(real_x), DIMENSION(NKVALUES) :: RCON

        REAL(real_x) :: O2, N2, MO2N2, H2O, K0, KI, Fc, K1, K2, K3, K4, NC, FCC, avo, m_h2so4, rhov, amu, Kb, pil

        REAL(real_x) :: Temp, DT_kpp

        REAL(real_x) :: KRO2NO, KRO2HO2, KAPHO2, KAPNO, KRO2NO3, KNO3AL, KDEC, KROPRIM, KROSEC, KCH3O2, K298CH3O2, K_h2so4


        ! O(1P) + NO -> NO2,   KMT01
        K0        =  1.0E-31 * (Temp / 300.)**(-1.6) * MO2N2
        KI        =  3.0E-11 * (Temp / 300.)**0.3
        Fc        =  0.85
        NC        =  0.75 - 1.27 * (LOG10(Fc))
        FCC       =  10**(LOG10(Fc)/(1+(LOG10(K0/KI)/NC)**2))
        RCON(1)   = (K0*KI) * FCC / (K0+KI)

        ! O(1P) + NO2 -> NO3,  KMT02
        K0        =  1.3E-31 * (Temp / 300.)**(-1.5) * MO2N2
        KI        =  2.3E-11 * (Temp / 300.)**0.24
        Fc        =  0.6
        NC        =  0.75 - 1.27 * (LOG10(Fc))
        FCC       =  10**(LOG10(Fc)/(1+(LOG10(K0/KI)/NC)**2))
        RCON(2)   = (K0*KI) * FCC / (K0+KI)

        ! NO2 + NO3 -> N2O5,   KMT03
        K0        =  3.6E-30 * (Temp / 300.)**(-4.1) * MO2N2
        KI        =  1.9E-12 * (Temp / 300.)**0.2
        Fc        =  0.35
        NC        =  0.75 - 1.27 * (LOG10(Fc))
        FCC       =  10**(LOG10(Fc)/(1+(LOG10(K0/KI)/NC)**2))
        RCON(3)   = (K0*KI) * FCC / (K0+KI)

        ! N2O5 -> NO2 + NO3,    KMT04
        K0        =  1.3E-3 * (Temp / 300.)**(-3.5) * EXP(-11000./Temp) * MO2N2
        KI        =  9.7E+14 * (Temp / 300.)**0.1 * EXP(-11080./Temp)
        Fc        =  0.35
        NC        =  0.75 - 1.27 * (LOG10(Fc))
        FCC       =  10**(LOG10(Fc)/(1+(LOG10(K0/KI)/NC)**2))
        RCON(4)   = (K0*KI) * FCC / (K0+KI)

        ! OH + CO -> HO2 + CO2,  KMT05
        RCON(5)   =  1.44D-13*(1+(MO2N2/4.2D+19))

        !HO2 + HO2        ->  H2O2 + O2,   KMT06
        RCON(6)   =  1 + (1.4E-21 * EXP(2200./Temp) * H2O)

        ! OH + NO -> HONO,    KMT07
        K0        =  7.40E-31 * (Temp / 300.)**(-2.4) * MO2N2
        KI        =  3.30E-11 * (Temp / 300.)**(-0.3)
        FC        =  EXP(-Temp / 1420.)
        NC        =  0.75 - 1.27 * (LOG10(Fc))
        FCC       =  10**(LOG10(Fc)/(1+(LOG10(K0/KI)/NC)**2))
        RCON(7)   = (K0*KI) * FCC / (K0+KI)

        ! OH + NO2 ->  HNO3,    KMT08
        K0        =  3.3E-30 * (Temp / 300.)**(-3.) * MO2N2
        KI        =  4.1E-11
        Fc        =  0.4
        NC        =  0.75 - 1.27 * (LOG10(Fc))
        FCC       =  10**(LOG10(Fc)/(1+(LOG10(K0/KI)/NC)**2))
        RCON(8)   = (K0*KI) * FCC / (K0+KI)

        ! HO2 + NO2 -> HO2NO2,    KMT09
        K0        =  1.8E-31 * (Temp / 300.)**(-3.2) * MO2N2
        KI        =  4.7E-12
        Fc        =  0.6
        NC        =  0.75 - 1.27 * (LOG10(Fc))
        FCC       =  10**(LOG10(Fc)/(1+(LOG10(K0/KI)/NC)**2))
        RCON(9)   = (K0*KI) * FCC / (K0+KI)

        ! HO2NO2 -> HO2 + NO2,    KMT10
        K0        =  4.1E-5 * EXP(-10650. / Temp) * MO2N2
        KI        =  4.8E+15 * EXP(-11170. / Temp)
        Fc        = 0.6
        NC        =  0.75 - 1.27 * (LOG10(Fc))
        FCC       =  10**(LOG10(Fc)/(1+(LOG10(K0/KI)/NC)**2))
        RCON(10)   = (K0*KI) * FCC / (K0+KI)

        ! OH + HNO3 -> NO3 + H2O,    KMT11
        K1        =  2.4E-14 * EXP(460./Temp)
        K3        =  6.5E-34 * EXP(1335./Temp)
        K4        =  2.7E-17 * EXP(2199./Temp)
        K2        =  (K3 * MO2N2) / (1 + (K3*MO2N2/K4))
        RCON(11)  =  K1 + K2

        ! OH + SO2 -> HSO3,    KMT12

        !Aitkinson 2004
        K0        =  4.5E-31 * (Temp/300.)**(-3.9) * MO2N2
        KI        =  1.3E-12 * (Temp/300.)**(-0.7)
        Fc        =  0.525
        NC        =  0.75 - 1.27 * (LOG10(Fc))
        FCC       =  10**(LOG10(Fc)/(1+(LOG10(K0/KI)/NC)**2))
        RCON(12)   = (K0*KI* FCC) / (K0+KI)

        !       MCM version 3.2
        !       K0        =  4.0E-31 * (Temp/300.)**(-3.3) * MO2N2
        !       KI        =  2.0E-12
        !       Fc        =  0.45
        !       RCON(12)  =  (K0 / (1. + K0/KI)) * Fc**(1 / (1 + (LOG10(K0 / KI))**2))

        ! CH3O2 + NO2 -> CH3O2NO2,   KMT13
        K0        =  2.50E-30 * (Temp/300)**(-5.5) * MO2N2
        KI        =  1.8E-11
        Fc        =  0.36
        NC        =  0.75 - 1.27 * (LOG10(Fc))
        FCC       =  10**(LOG10(Fc)/(1+(LOG10(K0/KI)/NC)**2))
        RCON(13)   = (K0*KI) * FCC / (K0+KI)

        ! CH3O2NO2 -> CH3O2 + NO2,    KMT14
        K0        =  9.00E-5 * EXP(-9690/Temp) * MO2N2
        KI        =  1.10E+16 * EXP(-10560/Temp)
        Fc        =  0.4
        NC        =  0.75 - 1.27 * (LOG10(Fc))
        FCC       =  10**(LOG10(Fc)/(1+(LOG10(K0/KI)/NC)**2))
        RCON(14)   = (K0*KI) * FCC / (K0+KI)

        !  ->     KMT15
        K0        =  8.6D-29 * (Temp/300)**(-3.1) * MO2N2
        KI        =  9.00D-12*(Temp/300)**(-0.85)
        Fc        =  0.48
        NC        =  0.75 - 1.27 * (LOG10(Fc))
        FCC       =  10**(LOG10(Fc)/(1+(LOG10(K0/KI)/NC)**2))
        RCON(15)   = (K0*KI) * FCC / (K0+KI)


        ! ->     KMT16
        K0        =  8.00E-27 * (Temp/300)**(-3.5) * MO2N2
        KI        =  3.00E-11*(Temp/300.)**(-1)
        Fc        =  0.5
        NC        =  0.75 - 1.27 * (LOG10(Fc))
        FCC       =  10**(LOG10(Fc)/(1+(LOG10(K0/KI)/NC)**2))
        RCON(16)   = (K0*KI) * FCC / (K0+KI)


        ! ->     KMT17
        K0        =  5.0E-30 * MO2N2 * (Temp/300.)**(-1.5)
        KI        =  1.0E-12
        Fc        =  0.17*EXP(-51/Temp) + EXP(-Temp/204)
        NC        =  0.75 - 1.27 * (LOG10(Fc))
        FCC       =  10**(LOG10(Fc)/(1+(LOG10(K0/KI)/NC)**2))
        RCON(17)   = (K0*KI*FCC) / (K0+KI)
        write(*,*) 'in old', RCON(17)

        !HO + CH3SCH3 ->    KMT18
        RCON(18)  =  9.5E-39 * O2 * EXP(5270/Temp) / (1 + 7.5E-29 * O2 * EXP(5610/Temp))

        ! KFPAN
        K0       = 3.28E-28*(TEMP/300)**(-6.87) * MO2N2
        KI       = 1.125E-11*(TEMP/300)**(-1.105)
        FCC       = 0.30
        NC        = 0.75-1.27*(LOG10(FCC))
        FC        = 10**(LOG10(FCC)/(1+(LOG10(K0/KI)/NC)**2))
        RCON(19)  =   (K0*KI)*FC/(K0+KI) ! KFPAN


        ! KBPAN
        K0        =  4.90E-03 * EXP(-12100/Temp) * MO2N2
        KI        =  5.40E+16 * EXP(-13830/Temp)   !Previously there was a mistake, so that it was 5.4E-D16
        Fc        =  0.3
        NC        =  0.75 - 1.27 * (LOG10(Fc))
        FCC       =  10**(LOG10(Fc)/(1+(LOG10(K0/KI)/NC)**2))
        RCON(20)  =  (K0*KI) * FCC / (K0+KI)

        !Here comes the simple rate coefficients:
        ! KRO2NO
        RCON(21) = 2.7D-12*EXP(360/Temp)

        ! KRO2HO2
        RCON(22) = 2.91D-13*EXP(1300/Temp)

        ! KAPHO2
        RCON(23) = 5.2D-13*EXP(980/Temp)

        ! KAPNO
        RCON(24) = 7.5D-12*EXP(290/Temp)

        ! KRO2NO3
        RCON(25) = 2.3D-12

        ! KNO3AL
        RCON(26) = 1.4D-12*EXP(-1860/Temp)

        ! KDEC
        RCON(27) = 1.00D+06

        ! KROPRIM
        RCON(28) = 2.50D-14*EXP(-300/Temp)

        ! KROSEC
        RCON(29) = 2.50D-14*EXP(-300/Temp)

        !KCH3O2
        RCON(30) = 1.03D-13*EXP(365/Temp)

        !K298CH3O2
        RCON(31) = 3.5D-13

        !K_h2so4, collision rate between two H2SO4 molecules
        avo  = 6.0221D23 ! avogadro no.,   /mol
        m_h2so4    = 98.0775 ! molar mass,  g/mol
        rhov = 1830. ! vapor molecule density,  kg/m^3
        amu  = 1.6605D-27 ! kg
        Kb   = 1.380658D-23 ! J/K
        pil = 3.14
        Fc = amu * (m_h2so4 / 2.)
        FCC = ((m_h2so4 / avo) * (1D-03 / rhov) * (6. / pil))**(1./3.)
        NC = SQRT(8. * Kb * Temp / (pil * Fc))
        RCON(32) = pil * FCC**2. * NC * 10.**6.
        RCON(32) = 0.

        !Criegee intermediate reaction rate: K_CI
        RCON(33) = 7.00D-14

        !!!Now come the rate coefficients from Aitkinson 2004
        ! 2OH-> H2O2,   KMT34
        K0        =  6.9D-31*(Temp/300.)**(-0.8)*MO2N2
        KI        =  2.6D-11
        Fc        =  0.50
        NC        =  0.75 - 1.27 * (LOG10(Fc))
        FCC       =  10**(LOG10(Fc)/(1+(LOG10(K0/KI)/NC)**2))
        RCON(34)   = (K0*KI) * FCC / (K0+KI)

        ! NO + NO2 -> N2O3,   KMT35
        K0        =  3.1D-34*(Temp/300.)**(-7.7)*MO2N2
        KI        =  7.9D-12*(Temp/300)**1.4
        Fc        =  0.6
        NC        =  0.75 - 1.27 * (LOG10(Fc))
        FCC       =  10**(LOG10(Fc)/(1+(LOG10(K0/KI)/NC)**2))
        RCON(35)   = (K0*KI) * FCC / (K0+KI)

        ! N2O3 + M -> NO + NO2 + M,   KMT36
        K0        =  1.9D-7*(Temp/300.)**(-8.7)*EXP(-4880/Temp)*MO2N2
        KI        =  4.7D15*(Temp/300)**(0.4)*EXP(-4880/Temp)
        Fc        =  0.6
        NC        =  0.75 - 1.27 * (LOG10(Fc))
        FCC       =  10**(LOG10(Fc)/(1+(LOG10(K0/KI)/NC)**2))
        RCON(36)   = (K0*KI) * FCC / (K0+KI)

        ! 2NO2 + M -> N2O4 + M,   KMT37
        K0        =  1.4D-33*(Temp/300.)**(-3.8)*MO2N2
        KI        =  1.D-12
        Fc        =  0.4
        NC        =  0.75 - 1.27 * (LOG10(Fc))
        FCC       =  10**(LOG10(Fc)/(1+(LOG10(K0/KI)/NC)**2))
        RCON(37)   = (K0*KI) * FCC / (K0+KI)

        ! N2O4 + M -> 2NO2 + M,   KMT38
        K0        =  1.3D-5*(Temp/300.)**(-3.8)*EXP(-6400/Temp)*MO2N2
        KI        =  1.15D16*EXP(-6400/Temp)
        Fc        =  0.4
        NC        =  0.75 - 1.27 * (LOG10(Fc))
        FCC       =  10**(LOG10(Fc)/(1+(LOG10(K0/KI)/NC)**2))
        RCON(38)   = (K0*KI) * FCC / (K0+KI)

        ! OH + CS2 + M -> HOCS2 + M,   KMT39
        K0        =  8.D-31*MO2N2
        KI        =  8.D-12
        Fc        =  0.8
        NC        =  0.75 - 1.27 * (LOG10(Fc))
        FCC       =  10**(LOG10(Fc)/(1+(LOG10(K0/KI)/NC)**2))
        RCON(39)   = (K0*KI) * FCC / (K0+KI)

        ! HOCS2 + M -> OH + CS2 + M,   KMT40
        K0        =  1.6D-6*EXP(-5160/Temp)*MO2N2
        KI        =  1.6D13*EXP(-5160/Temp)
        Fc        =  0.8
        NC        =  0.75 - 1.27 * (LOG10(Fc))
        FCC       =  10**(LOG10(Fc)/(1+(LOG10(K0/KI)/NC)**2))
        RCON(40)   = (K0*KI) * FCC / (K0+KI)

        ! HS + NO + M -> HSNO + M,   KMT41
        K0        =  2.4D-31*(Temp/300)**(-2.5)*MO2N2
        KI        =  2.7D-11
        Fc        =  0.6
        NC        =  0.75 - 1.27 * (LOG10(Fc))
        FCC       =  10**(LOG10(Fc)/(1+(LOG10(K0/KI)/NC)**2))
        RCON(41)   = (K0*KI) * FCC / (K0+KI)

        ! CH3S + NO + M -> CH3SNO + M,   KMT42
        K0        =  3.3D-29*(Temp/300)**(-4)*MO2N2
        KI        =  4.D-11
        Fc        =  0.54
        NC        =  0.75 - 1.27 * (LOG10(Fc))
        FCC       =  10**(LOG10(Fc)/(1+(LOG10(K0/KI)/NC)**2))
        RCON(42)   = (K0*KI) * FCC / (K0+KI)

    END SUBROUTINE

END MODULE ! Chemsitry
