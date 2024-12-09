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


MODULE OUTPUT
use netcdf
use CONSTANTS
use INPUT
use second_Monitor
USE SECOND_REACTIVITY, ONLY : reactivity_name, NREACTIVITY
use AUXILLARIES
USE PSD_scheme, ONLY: old_PSD,get_composition, get_dp, get_volume, get_mass, get_conc,G_COAG_SINK

IMPLICIT NONE
private

type parsave
  integer :: i
  integer :: d
  integer :: type
  character(10) :: u
  character(30) :: name
end type parsave

! Storing of indices of files
INTEGER, PARAMETER :: N_FILES = 3
INTEGER, PARAMETER :: shuff=1, compress=1, compression=9

INTEGER       :: ncfile_ids(N_FILES)
CHARACTER(200) :: ncfile_names(N_FILES) = (['General  ', 'Chemistry', 'Particles'])

real,allocatable     :: TIMESERIES(:,:)
integer              :: i_TIMESERIES = 1
integer              :: reclen

INTEGER, allocatable :: shifter_ind(:)
INTEGER, allocatable :: multipl_ind(:)
INTEGER, allocatable :: mods_ind(:)
INTEGER, allocatable :: chem_ind(:)
INTEGER, allocatable :: acdc_ind(:)
INTEGER, allocatable :: par_ind(:)

INTEGER :: dtime_id
INTEGER :: dbins_id
INTEGER :: dcond_id
INTEGER :: dorg_id
INTEGER :: dvbs_id
INTEGER :: dchrg_id
INTEGER :: dstring_id
INTEGER :: dconstant_id
INTEGER :: timearr_id
INTEGER :: hrsarr_id
INTEGER :: gJ_out_SUM_id
INTEGER :: gJ_out_TOT_id
INTEGER :: gCS_calc_id

public :: OPEN_FILES, SAVE_GASES,CLOSE_FILES,INITIALIZE_WITH_LAST,WRITE_GAS_BINARY_DUMP,TIMESERIES,i_TIMESERIES

type(parsave) :: parbuf(25)
type(parsave), ALLOCATABLE :: savepar(:)
REAL(dp), ALLOCATABLE,PUBLIC :: Depos_composition(:) ! cumulative sum of aerosol mass (composition) lost to walls [kg/m³]
REAL(dp), ALLOCATABLE,PUBLIC :: c_org_wall_old(:)    ! total concentrations of organic vapours on the walls in the previous time step [molecules/chamber]
REAL(dp), ALLOCATABLE,PUBLIC :: c_org_wall(:)        ! total concentrations of organic vapours on the walls [molecules/chamber]
                                                     !   :.. divide with chamber surface area to get [molecules/m²], dim=n_cond_org
REAL(dp), ALLOCATABLE,PUBLIC :: CONDENSED_MASS(:)    ! average mass transfer. Positive value = Mass goes to particles [kg/m3/s]
REAL(dp), ALLOCATABLE,PUBLIC :: CONDENSED_MASS_VBS(:,:)  ! size resolved average mass transfer from VBS bins. Positive value = Mass goes to particles [kg/m3/s]
REAL(dp), ALLOCATABLE,PUBLIC :: VAP_DEP_MASS_WALLS(:)    ! average mass transfer to walls. Positive value = Mass goes to walls [kg/m3/s]

CONTAINS

! --------------------------------------------------------------------------------------------------------------------
! CREATE (OR OVERWRITE OLD) THREE TEMPORARY NETCDF-FILES TO STORE OUTPUT. AFTER FINISHED RUN TEMPORARY FILES ARE RENAMED
! TO THE FINAL FILES. THIS WAY KILLED RUN DOESN'T SPOIL PREVIOUS SUCCESFUL OUTPUT.
! --------------------------------------------------------------------------------------------------------------------
SUBROUTINE OPEN_FILES(filename, Description, CurrentChem,CurrentVers,SHA, MODS, CH_GAS,reactivities, vapours,V_chamber)
  IMPLICIT NONE

  CHARACTER(LEN=*), INTENT(IN)    :: filename
  CHARACTER(LEN=*), INTENT(IN)    :: CurrentChem
  CHARACTER(LEN=*), INTENT(IN)    :: CurrentVers
  CHARACTER(LEN=*), INTENT(IN)    :: SHA
  CHARACTER(*), INTENT(IN)        :: Description
  TYPE(input_mod),INTENT(INOUT)   :: MODS(:)
  TYPE(vapour_ambient),INTENT(IN) :: vapours
  REAL(dp), INTENT(IN)            :: CH_GAS(:)
  REAL(dp), INTENT(IN)            :: reactivities(:)
  REAL(dp), INTENT(IN)            :: V_chamber
  CHARACTER(255)                  :: PROGRAM_NAME
  INTEGER                         :: i,j,jj,ioi,lenD
  INTEGER, PARAMETER              :: textdim = len(SPC_NAMES(1))
  CHARACTER(textdim), ALLOCATABLE :: COND_NAMES(:)
  CHARACTER(len=256)              :: realdate, realtime
  CHARACTER(len=16)               :: compounds,cache

  call date_and_time(realdate,realtime)


  ! ---------------------------
  ! ------dimension bits-------
  ! time | bins | condensibles | name  | VBS | n_cond_org
  !   1  |   2  |      4       |   8   | 16  |    32
  ! eg. number concentration:
  ! time | bins | condensibles
  !  1      1          0         => 011 => parbuf(i)%d = 3
  i=1
  parbuf(i)%name = 'VAPOURS'               ; parbuf(i)%u = '[]'         ; parbuf(i)%d = -12 ; parbuf(i)%type = NF90_CHAR   ; i=i+1
  parbuf(i)%name = 'NUMBER_CONCENTRATION'  ; parbuf(i)%u = '[/cm^3]'    ; parbuf(i)%d = 3   ; parbuf(i)%type = NF90_DOUBLE ; i=i+1
  parbuf(i)%name = 'INPUT_CONCENTRATION'   ; parbuf(i)%u = '[/cm^3]'    ; parbuf(i)%d = 3   ; parbuf(i)%type = NF90_DOUBLE ; i=i+1
  parbuf(i)%name = 'DIAMETER'              ; parbuf(i)%u = '[m]'        ; parbuf(i)%d = 3   ; parbuf(i)%type = NF90_DOUBLE ; i=i+1
  parbuf(i)%name = 'GROWTH_RATE'           ; parbuf(i)%u = '[nm/h]'     ; parbuf(i)%d = 3   ; parbuf(i)%type = NF90_DOUBLE ; i=i+1
  parbuf(i)%name = 'COAG_SINK'             ; parbuf(i)%u = '[/s]'       ; parbuf(i)%d = 3   ; parbuf(i)%type = NF90_DOUBLE ; i=i+1
  parbuf(i)%name = 'MASS'                  ; parbuf(i)%u = '[kg]'       ; parbuf(i)%d = 3   ; parbuf(i)%type = NF90_DOUBLE ; i=i+1
  parbuf(i)%name = 'PARTICLE_COMPOSITION'  ; parbuf(i)%u = '[kg/particle]'; parbuf(i)%d = 7 ; parbuf(i)%type = NF90_DOUBLE ; i=i+1
  parbuf(i)%name = 'MASS_FLUX_ON_PAR'      ; parbuf(i)%u = '[kg/m^3/s]' ; parbuf(i)%d = 5   ; parbuf(i)%type = NF90_DOUBLE ; i=i+1
  if (vapours%n_cond_org>1) THEN
    parbuf(i)%name = 'VBS_MATRIX'            ; parbuf(i)%u = '[]'         ; parbuf(i)%d = 49  ; parbuf(i)%type = NF90_INT    ; i=i+1
    parbuf(i)%name = 'MASS_FLUX_ON_PAR_VBS'  ; parbuf(i)%u = '[µg/m^3/s]' ; parbuf(i)%d = 19  ; parbuf(i)%type = NF90_DOUBLE ; i=i+1
  END IF
  if (Deposition) THEN
      parbuf(i)%name = 'DEPOSITED_PAR_COMP' ; parbuf(i)%u = '[kg/m^3]'  ; parbuf(i)%d = 5   ; parbuf(i)%type = NF90_DOUBLE ; i=i+1
      parbuf(i)%name = 'PARTICLE_LOSS_RATE' ; parbuf(i)%u = '[/s]'      ; parbuf(i)%d = 3   ; parbuf(i)%type = NF90_DOUBLE ; i=i+1
  END IF
  if (Chem_Deposition) THEN
      parbuf(i)%name = 'MASS_FLUX_ON_WALLS' ; parbuf(i)%u = '[kg/s]'    ; parbuf(i)%d = 5   ; parbuf(i)%type = NF90_DOUBLE ; i=i+1
  END IF

  allocate(savepar(i-1))
  savepar = parbuf(1:i-1)

  if (Aerosol_flag) THEN
      ALLOCATE(COND_NAMES(vapours%n_cond_tot))
      COND_NAMES = vapours%vapour_names
  END IF

  ! Print run description
  lenD = LEN(TRIM(Description))
  print FMT_SUB, 'Description for run:'
  i=0
  do while (i< 1+lenD/90)
    print FMT_MSG, '    '//TRIM(Description(((i*90)+1):90*(i+1)))
    i = i+1
  end do
  ! end print run description

  print FMT_HDR, 'PREPARING OUTPUT FILES'
  print FMT_SUB, 'NetCDF version: '//trim(nf90_inq_libvers())
  print FMT_SUB, 'Create files to: '//TRIM(filename)

  DO I=1, N_FILES
    if (I == 3.and..not.Aerosol_flag) cycle

    ncfile_names(I) = trim(filename)//'/'//TRIM(ncfile_names(I))
    ! Clearing file; Opening file. Overwrites
    open(720+I, FILE=TRIM(ncfile_names(I))//'.tmp', iostat = ioi)
    CALL handle_file_io(ioi, TRIM(ncfile_names(I))//'.tmp', 'Terminating when trying to open netCDF-file, does the CASE and RUN directory exist')
    close(720+I)

    ! Added compression for particles.nc, so we need to use netCDF4-file. Here used in classic mode
    call handler(__LINE__, nf90_create(TRIM(ncfile_names(I))//'.tmp', IOR(NF90_NETCDF4, NF90_CLASSIC_MODEL), ncfile_ids(I)) )

    ! Defining dimensions: time(unlimited), size sections, vapor_species
    call handler(__LINE__, nf90_def_dim(ncfile_ids(I), "time",NF90_UNLIMITED, dtime_id) )
    IF (I==1) call handler(__LINE__, nf90_def_dim(ncfile_ids(I), "formation_rates",4, dchrg_id) )
    IF (I==3) call handler(__LINE__, nf90_def_dim(ncfile_ids(I), "string",textdim, dstring_id) )
    IF (I>1)  call handler(__LINE__, nf90_def_dim(ncfile_ids(I), "bins",n_bins_par, dbins_id) )
    IF ((I == 3) .and. Aerosol_flag) call handler(__LINE__, nf90_def_dim(ncfile_ids(I), "vapours",vapours%n_cond_tot, dcond_id) )

    IF ((I == 3) .and. vapours%n_cond_org>1) call handler(__LINE__, nf90_def_dim(ncfile_ids(I), "org vapours",vapours%n_cond_org-1, dorg_id) )

    IF ((I == 3) .and. Aerosol_flag) call handler(__LINE__, nf90_def_dim(ncfile_ids(I), "vbs_bins",NVBS, dvbs_id) )
    call handler(__LINE__, nf90_def_dim(ncfile_ids(I), "Constant",1, dconstant_id) )

    ! Create attributes for general stuff
    CALL get_command_argument(0, PROGRAM_NAME)
    call handler(__LINE__, nf90_put_att(ncfile_ids(I), NF90_GLOBAL, 'Information', '(c) Multiscale Modelling Group and (c) Computational Aerosol Physics Group (ACDC)'))
    call handler(__LINE__, nf90_put_att(ncfile_ids(I), NF90_GLOBAL, 'Contact', 'arca@helsinki.fi (ARCA box), tinja.olenius@alumni.helsinki.fi (ACDC)'))
    call handler(__LINE__, nf90_put_att(ncfile_ids(I), NF90_GLOBAL, 'Software', 'ARCA box '//TRIM(CurrentVers)))
    call handler(__LINE__, nf90_put_att(ncfile_ids(I), NF90_GLOBAL, 'git_hash', TRIM(SHA)))
    call handler(__LINE__, nf90_put_att(ncfile_ids(I), NF90_GLOBAL, 'Package_Name', TRIM(PROGRAM_NAME(3:))))
    call handler(__LINE__, nf90_put_att(ncfile_ids(I), NF90_GLOBAL, 'Chemistry_module', TRIM(CurrentChem)))
    call handler(__LINE__, nf90_put_att(ncfile_ids(I), NF90_GLOBAL, 'Notes', TRIM(Description)))
    call handler(__LINE__, nf90_put_att(ncfile_ids(I), NF90_GLOBAL, 'experiment', filename))
    call handler(__LINE__, nf90_put_att(ncfile_ids(I), NF90_GLOBAL, 'INITFILE', Fname_init))
    call handler(__LINE__, nf90_put_att(ncfile_ids(I), NF90_GLOBAL, 'Real_start_date', TRIM(realdate(1:4))//'/'//TRIM(realdate(5:6))//'/'//TRIM(realdate(7:8)) ))
    call handler(__LINE__, nf90_put_att(ncfile_ids(I), NF90_GLOBAL, 'Real_start_time', TRIM(realtime(1:2))//':'//TRIM(realtime(3:4))//':'//TRIM(realtime(5:6)) ))
    call handler(__LINE__, nf90_put_att(ncfile_ids(I), NF90_GLOBAL, 'Nominal_save_interval_s', GTIME%FSAVE_INTERVAL ))
    call handler(__LINE__, nf90_put_att(ncfile_ids(I), NF90_GLOBAL, 'missing_value', -9999d0 ))
    if (Chem_Deposition) &
      call handler(__LINE__, nf90_put_att(ncfile_ids(I), NF90_GLOBAL, 'Chamber_volume', V_chamber ))
    call handler(__LINE__, nf90_def_var(ncfile_ids(I), "TIME_IN_SEC", NF90_DOUBLE, dtime_id, timearr_id))
    call handler(__LINE__, nf90_def_var(ncfile_ids(I), "TIME_IN_HRS", NF90_DOUBLE, dtime_id, hrsarr_id))
    ! COMPRESSION
    call handler(__LINE__, nf90_def_var_deflate(ncfile_ids(I), timearr_id, shuff, compress, compression) )
    call handler(__LINE__, nf90_def_var_deflate(ncfile_ids(I), hrsarr_id, shuff, compress, compression) )
    call handler(__LINE__, nf90_put_att(ncfile_ids(I), timearr_id, 'unit' , '[s]'))
    call handler(__LINE__, nf90_put_att(ncfile_ids(I), hrsarr_id, 'unit' , '[hrs]'))
    call handler(__LINE__, nf90_put_att(ncfile_ids(I), timearr_id, 'type' , 'time'))
    call handler(__LINE__, nf90_put_att(ncfile_ids(I), hrsarr_id, 'type' , 'time'))

  END DO

  ALLOCATE(multipl_ind(size(MODS)))
  ALLOCATE(shifter_ind(size(MODS)))
  ALLOCATE(mods_ind(size(MODS)))
  ALLOCATE(chem_ind(size(CH_GAS)+NREACTIVITY))
  if (Aerosol_flag) ALLOCATE(par_ind(size(vapours%vapour_names)))
  ALLOCATE(acdc_ind(size(G_ACDC)))
  if (size(G_ACDC)>0) THEN
    do i=1,size(G_ACDC)
      if (G_ACDC(i)%inuse) THEN
        WRITE(cache, '(a,i0.2)') 'ACDC_',i
        call handler(__LINE__, nf90_put_att(ncfile_ids(1), NF90_GLOBAL, TRIM(cache)//'_System_file', TRIM(G_ACDC(i)%SYSTEM_FILE) ))
        call handler(__LINE__, nf90_put_att(ncfile_ids(1), NF90_GLOBAL, TRIM(cache)//'_Energy_file', TRIM(G_ACDC(i)%ENERGY_FILE) ))
        call handler(__LINE__, nf90_put_att(ncfile_ids(1), NF90_GLOBAL, TRIM(cache)//'_Dipole_file', TRIM(G_ACDC(i)%DIPOLE_FILE) ))
        call handler(__LINE__, nf90_put_att(ncfile_ids(1), NF90_GLOBAL, TRIM(cache)//'_Name'       , TRIM(G_ACDC(i)%NICKNAME   ) ))
      END IF
    end do
  END IF


  I=1 ! GENERAL FILE
  do j = 1,size(MODS)
    IF (TRIM(MODS(j)%name) /= '#') THEN
      IF (j<LENV .or. INDRELAY_CH(j)>0 .or. MODS(j)%ISPROVIDED) THEN
        call handler(__LINE__, nf90_def_var(ncfile_ids(I), TRIM(MODS(j)%name), NF90_DOUBLE, dtime_id, mods_ind(j)) )
        call handler(__LINE__, nf90_def_var_deflate(ncfile_ids(I), mods_ind(j), shuff, compress, compression) )
        call handler(__LINE__, nf90_put_att(ncfile_ids(I), mods_ind(j), 'unit' , TRIM(UNITS(MODS(I)%UNIT))))
      END IF

      IF ((ABS(MODS(j)%shift) > 1d-100) .or. (ABS(MODS(j)%multi-1d0) > 1d-100)) THEN
        call handler(__LINE__, nf90_def_var(ncfile_ids(I), TRIM(MODS(j)%name)//'_Multipl', NF90_DOUBLE, dtime_id, multipl_ind(j)))
        call handler(__LINE__, nf90_def_var(ncfile_ids(I), TRIM(MODS(j)%name)//'_Shifter', NF90_DOUBLE, dtime_id, shifter_ind(j)))

        call handler(__LINE__, nf90_def_var_deflate(ncfile_ids(I), multipl_ind(j), shuff, compress, compression) )
        call handler(__LINE__, nf90_def_var_deflate(ncfile_ids(I), shifter_ind(j), shuff, compress, compression) )

        call handler(__LINE__, nf90_put_att(ncfile_ids(I), multipl_ind(j), 'unit' , '[]'))
        call handler(__LINE__, nf90_put_att(ncfile_ids(I), shifter_ind(j), 'unit' , '[same]'))

      END IF
    END IF
  end do

  do jj=1,size(G_ACDC)
    if (G_ACDC(jj)%inuse) THEN
      call handler(__LINE__, nf90_def_var(ncfile_ids(I), 'J_ACDC_'//TRIM(i2chr(jj))//'_CM3', NF90_DOUBLE, [dchrg_id,dtime_id], acdc_ind(jj)))
      call handler(__LINE__, nf90_def_var_deflate(ncfile_ids(I), acdc_ind(jj),shuff, compress, compression) )
      call handler(__LINE__, nf90_put_att(ncfile_ids(I), acdc_ind(jj), 'unit' , '[1/s/cm^3]'))
      call handler(__LINE__, nf90_put_att(ncfile_ids(I), acdc_ind(jj), 'Entries_in_J_vector' , 'J_tot/J_neut/J_pos/J_neg'))
      call handler(__LINE__, nf90_put_att(ncfile_ids(I), acdc_ind(jj), 'Missing_value_in_J_vector' , '-1'))
      call handler(__LINE__, nf90_put_att(ncfile_ids(I), acdc_ind(jj), 'Missing_value_means' , 'Charge_not_in_ACDC_system'))
      do j=1,size(G_ACDC(jj)%ACDC_MONOMER_NAMES)
        if (G_ACDC(jj)%ACDC_LINK_IND(j)>0) compounds = MODS( (G_ACDC(jj)%ACDC_LINK_IND(j)) )%NAME
        if (G_ACDC(jj)%ACDC_LINK_IND(j)<0) compounds = SPC_NAMES(-1*G_ACDC(jj)%ACDC_LINK_IND(j))
        call handler(__LINE__, nf90_put_att(ncfile_ids(I), acdc_ind(jj), 'Component_links'//TRIM(i2chr(j)) , &
                                          TRIM( G_ACDC(jj)%ACDC_MONOMER_NAMES(j) )//':'//TRIM(compounds) ))
      end do
    end if
  end do

  call handler(__LINE__, nf90_def_var(ncfile_ids(I), 'J_ACDC_SUM_CM3', NF90_DOUBLE, dtime_id, gJ_out_SUM_id))
  call handler(__LINE__, nf90_def_var(ncfile_ids(I), 'J_TOTAL_CM3', NF90_DOUBLE, dtime_id, gJ_out_TOT_id))
  call handler(__LINE__, nf90_def_var(ncfile_ids(I), 'CS_CALC', NF90_DOUBLE, dtime_id, gCS_calc_id))
  call handler(__LINE__, nf90_def_var_deflate(ncfile_ids(I), gJ_out_SUM_id,shuff, compress, compression) )
  call handler(__LINE__, nf90_def_var_deflate(ncfile_ids(I), gJ_out_TOT_id,shuff, compress, compression) )
  call handler(__LINE__, nf90_def_var_deflate(ncfile_ids(I), gCS_calc_id,shuff, compress, compression) )
  call handler(__LINE__, nf90_put_att(ncfile_ids(I), gJ_out_SUM_id, 'unit' , '[1/s/cm^3]'))
  call handler(__LINE__, nf90_put_att(ncfile_ids(I), gJ_out_TOT_id, 'unit' , '[1/s/cm^3]'))
  call handler(__LINE__, nf90_put_att(ncfile_ids(I), gCS_calc_id, 'unit' , '[1/s]'))

  I=2 ! Chemical file
  do j = 1,size(CH_GAS)
    call handler(__LINE__, nf90_def_var(ncfile_ids(I), TRIM(SPC_NAMES(j)), NF90_DOUBLE, dtime_id, chem_ind(j)) )
    call handler(__LINE__, nf90_def_var_deflate(ncfile_ids(I), chem_ind(j), shuff, compress, compression) )
    call handler(__LINE__, nf90_put_att(ncfile_ids(I), chem_ind(j), 'unit' , '1/cm^3'))
    call handler(__LINE__, nf90_put_att(ncfile_ids(I), chem_ind(j), 'unit' , 'gas'))
  end do

  if (NREACTIVITY>0) THEN
    do j = 1,NREACTIVITY
      call handler(__LINE__, nf90_def_var(ncfile_ids(I), TRIM(reactivity_name(j)), NF90_DOUBLE, dtime_id, chem_ind(size(CH_GAS)+j)) )
      call handler(__LINE__, nf90_def_var_deflate(ncfile_ids(I), chem_ind(size(CH_GAS)+j), shuff, compress, compression) )
      call handler(__LINE__, nf90_put_att(ncfile_ids(I), chem_ind(size(CH_GAS)+j), 'unit' , '1/s'))
      call handler(__LINE__, nf90_put_att(ncfile_ids(I), chem_ind(size(CH_GAS)+j), 'type' , 'reactivity'))
    end do
  end if


  I=3 ! Particle file.
  if (Aerosol_flag) THEN
    do j = 1,size(savepar)
      if (savepar(J)%d == 0) call handler(__LINE__, nf90_def_var(ncfile_ids(I), TRIM(  savepar(J)%name  ), savepar(J)%type, dconstant_id  , savepar(J)%i) )
      if (savepar(J)%d == 3) call handler(__LINE__, nf90_def_var(ncfile_ids(I), TRIM(  savepar(J)%name  ), savepar(J)%type, ([dbins_id, dtime_id])  , savepar(J)%i) )
      if (savepar(J)%d == 4) call handler(__LINE__, nf90_def_var(ncfile_ids(I), TRIM(  savepar(J)%name  ), savepar(J)%type, dcond_id  , savepar(J)%i) )
      if (savepar(J)%d == 5) call handler(__LINE__, nf90_def_var(ncfile_ids(I), TRIM(  savepar(J)%name  ), savepar(J)%type, ([dcond_id,dtime_id])  , savepar(J)%i) )
      if (savepar(J)%d == 7) call handler(__LINE__, nf90_def_var(ncfile_ids(I), TRIM(  savepar(J)%name  ), savepar(J)%type, ([dcond_id,dbins_id,dtime_id])  , savepar(J)%i) )
      if (savepar(J)%d == 19) call handler(__LINE__, nf90_def_var(ncfile_ids(I), TRIM(  savepar(J)%name  ), savepar(J)%type, ([dvbs_id,dbins_id,dtime_id])  , savepar(J)%i) )
      if (savepar(J)%d == 49) call handler(__LINE__, nf90_def_var(ncfile_ids(I), TRIM(  savepar(J)%name  ), savepar(J)%type, ([dvbs_id,dorg_id,dtime_id])  , savepar(J)%i) )
      if (savepar(J)%d == -12) call handler(__LINE__, nf90_def_var(ncfile_ids(I), TRIM(  savepar(J)%name  ), savepar(J)%type, ([dstring_id,dcond_id])  , savepar(J)%i) )
      call handler(__LINE__, nf90_def_var_fill(ncfile_ids(I), savepar(J)%i, 1, -9999d0) )
      call handler(__LINE__, nf90_def_var_deflate(ncfile_ids(I), savepar(J)%i, shuff, compress, compression) )
      call handler(__LINE__, nf90_put_att(ncfile_ids(I), savepar(J)%i, 'unit' , savepar(J)%u))
      call handler(__LINE__, nf90_put_att(ncfile_ids(I), savepar(J)%i, 'type' , 'misc'))
    end do

    do j = 1,vapours%n_cond_tot
      call handler(__LINE__, nf90_def_var(ncfile_ids(I), TRIM(  vapours%vapour_names(j)  ), NF90_DOUBLE, dtime_id, par_ind(j)) )
      call handler(__LINE__, nf90_def_var_deflate(ncfile_ids(I), par_ind(j), shuff, compress, compression) )
      call handler(__LINE__, nf90_put_att(ncfile_ids(I), par_ind(j), 'unit' , '1/cm^3'))
      call handler(__LINE__, nf90_put_att(ncfile_ids(I), par_ind(j), 'type' , TRIM(i2chr(vapours%cond_type(j))) ))
      call handler(__LINE__, nf90_put_att(ncfile_ids(I), par_ind(j), 'Molar_mass_kg' , vapours%molar_mass(j) ))
      call handler(__LINE__, nf90_put_att(ncfile_ids(I), par_ind(j), 'bulk_density' , vapours%density(j) ))
      if (VBS_CSAT) THEN
        call handler(__LINE__, nf90_put_att(ncfile_ids(I), par_ind(j), 'Csat_300' , vapours%psat_a(j) ))
        call handler(__LINE__, nf90_put_att(ncfile_ids(I), par_ind(j), 'Enthalpy' , vapours%psat_b(j) ))
      ELSE
        call handler(__LINE__, nf90_put_att(ncfile_ids(I), par_ind(j), 'Psat_A' , vapours%psat_a(j) ))
        call handler(__LINE__, nf90_put_att(ncfile_ids(I), par_ind(j), 'Psat_B' , vapours%psat_b(j) ))
      END IF
      if (vapours%cond_type(j)<2) &
        call handler(__LINE__, nf90_put_att(ncfile_ids(I), par_ind(j), 'surface_tension' , vapours%surf_tension(j) ))
      if (vapours%cond_type(j)<2) &
        call handler(__LINE__, nf90_put_att(ncfile_ids(I), par_ind(j), 'alpha_wall' , vapours%alphawall(j) ))
    end do

  end if

  ! End definition mode
  DO I=1,N_FILES
    if (I == 3.and..not.Aerosol_flag) cycle
    call handler(__LINE__, nf90_enddef(ncfile_ids(I)))
  END DO

  if (Aerosol_flag) THEN
    do j = 1,size(vapours%vapour_names)
      call handler(__LINE__, nf90_put_var(ncfile_ids(3), savepar(1)%i, COND_NAMES))
    end do
  END IF

  print FMT_LEND,

END SUBROUTINE OPEN_FILES


! ====================================================================================================================
! Here the input is written to netcdf-files.
! --------------------------------------------------------------------------------------------------------------------
SUBROUTINE SAVE_GASES(TSTEP_CONC,MODS,CH_GAS,reactivities,conc_vapours, VAPOURS, save_measured, GR, Lossrate)
  IMPLICIT NONE
  type(input_mod), INTENT(in)     :: MODS(:)
  real(dp), INTENT(in)            :: TSTEP_CONC(:)
  real(dp), INTENT(in)            :: CH_GAS(:)
  real(dp), INTENT(in)            :: reactivities(:)
  real(dp), INTENT(in)            :: conc_vapours(:)
  real(dp), INTENT(in)            :: Lossrate(:)
  TYPE(vapour_ambient),INTENT(IN) :: vapours
  real(dp), INTENT(in)            :: save_measured(:)
  real(dp), INTENT(in)            :: GR(:)
  INTEGER                         :: i,j,jj

  DO I = 1,N_FILES
    if (I == 3.and..not.Aerosol_flag) cycle
    call handler(__LINE__, nf90_put_var(ncfile_ids(I), timearr_id, GTIME%sec, (/GTIME%ind_netcdf/) ))
    call handler(__LINE__, nf90_put_var(ncfile_ids(I), hrsarr_id, GTIME%hrs, (/GTIME%ind_netcdf/) ))
  END DO

  I=1 ! General file

  do j = 1,size(MODS)
    IF ((TRIM(MODS(j)%name) /= '#')) THEN
      IF (j<LENV .or. INDRELAY_CH(j)>0 .or. MODS(j)%ISPROVIDED) THEN
        call handler(__LINE__, nf90_put_var(ncfile_ids(I), mods_ind(j), TSTEP_CONC(j), (/GTIME%ind_netcdf/)) )
      END IF
      IF ((ABS(MODS(j)%shift-0d0) > 1d-100) .or. (ABS(MODS(j)%multi-1d0) > 1d-100)) THEN
        call handler(__LINE__, nf90_put_var(ncfile_ids(I), multipl_ind(j), MODS(j)%multi, (/GTIME%ind_netcdf/) ) )
        call handler(__LINE__, nf90_put_var(ncfile_ids(I), shifter_ind(j), MODS(j)%shift,  (/GTIME%ind_netcdf/) ) )
      END IF
    END IF
  end do

  do jj=1,size(G_ACDC)
      if (G_ACDC(jj)%inuse) THEN
          call handler(__LINE__, nf90_put_var(ncfile_ids(I), acdc_ind(jj), G_ACDC(jj)%J_OUT_CM3 , start=(/1,GTIME%ind_netcdf/), count=(/4/)  ) )
      end if
  end do

  call handler(__LINE__, nf90_put_var(ncfile_ids(I), gJ_out_SUM_id, SUM(G_ACDC%J_OUT_CM3(1)), (/GTIME%ind_netcdf/)) )
  call handler(__LINE__, nf90_put_var(ncfile_ids(I), gJ_out_TOT_id, 1d-6*(J_TOTAL_M3), (/GTIME%ind_netcdf/)) )
  call handler(__LINE__, nf90_put_var(ncfile_ids(I), gCS_calc_id, GCS, (/GTIME%ind_netcdf/)) )

  I=2 ! Chemical file
  do j = 1,size(CH_GAS)
    call handler(__LINE__, nf90_put_var(ncfile_ids(I), chem_ind(j), CH_GAS(j), (/GTIME%ind_netcdf/)) )
  end do

  if (NREACTIVITY>0) THEN
    do j = 1,NREACTIVITY
      call handler(__LINE__, nf90_put_var(ncfile_ids(I), chem_ind(size(CH_GAS)+j), reactivities(j), (/GTIME%ind_netcdf/)) )
    end do
  end if


  I=3 ! Particle file

  if (Aerosol_flag) THEN

    do j = 1,vapours%n_cond_tot
      call handler(__LINE__, nf90_put_var(ncfile_ids(I), par_ind(j), conc_vapours(j), (/GTIME%ind_netcdf/)) )
    end do

    do j = 1,size(savepar)
      if (savepar(j)%name == 'NUMBER_CONCENTRATION'  ) call handler(__LINE__, nf90_put_var(ncfile_ids(I), savepar(J)%i, &
          1d-6*(get_conc(old_PSD)), start=(/1,GTIME%ind_netcdf/), count=(/n_bins_par/)))

      if (savepar(j)%name == 'INPUT_CONCENTRATION'   ) call handler(__LINE__, nf90_put_var(ncfile_ids(I), savepar(J)%i, &
          save_measured, start=(/1,GTIME%ind_netcdf/), count=(/n_bins_par/)))

      if (savepar(j)%name == 'DIAMETER'              ) call handler(__LINE__, nf90_put_var(ncfile_ids(I), savepar(J)%i, &
          get_dp(old_PSD), start=(/1,GTIME%ind_netcdf/), count=(/n_bins_par/)))

      if (savepar(j)%name == 'GROWTH_RATE'           ) call handler(__LINE__, nf90_put_var(ncfile_ids(I), savepar(J)%i, &
          GR, start=(/1,GTIME%ind_netcdf/), count=(/n_bins_par/)))

      if (savepar(j)%name == 'COAG_SINK'             ) call handler(__LINE__, nf90_put_var(ncfile_ids(I), savepar(J)%i, &
          G_COAG_SINK, start=(/1,GTIME%ind_netcdf/), count=(/n_bins_par/)))

      if (savepar(j)%name == 'CORE_VOLUME'           ) call handler(__LINE__, nf90_put_var(ncfile_ids(I), savepar(J)%i, &
          get_volume(old_PSD), start=(/1,GTIME%ind_netcdf/), count=(/n_bins_par/)))

      if (savepar(j)%name == 'MASS'                  ) call handler(__LINE__, nf90_put_var(ncfile_ids(I), savepar(J)%i, &
          get_mass(old_PSD), start=(/1,GTIME%ind_netcdf/), count=(/n_bins_par/)))

      if (savepar(j)%name == 'PARTICLE_COMPOSITION'  ) call handler(__LINE__, nf90_put_var(ncfile_ids(I), savepar(J)%i, &
          TRANSPOSE(get_composition(old_PSD)), (/1,1,GTIME%ind_netcdf/), (/vapours%n_cond_tot,n_bins_par,1/)))

      if (savepar(j)%name == 'VBS_MATRIX'            ) call handler(__LINE__, nf90_put_var(ncfile_ids(I), savepar(J)%i, &
          TRANSPOSE(VAPOUR_PROP%VBS_BINS), (/1,1,GTIME%ind_netcdf/), (/NVBS,vapours%n_cond_org-1,1/)))

      if (savepar(j)%name == 'MASS_FLUX_ON_PAR'      ) call handler(__LINE__, nf90_put_var(ncfile_ids(I), savepar(J)%i, &
          CONDENSED_MASS/Gtime%FSAVE_INTERVAL, start=(/1,GTIME%ind_netcdf/), count=(/vapours%n_cond_tot/)))

      if (savepar(j)%name == 'MASS_FLUX_ON_PAR_VBS'  ) call handler(__LINE__, nf90_put_var(ncfile_ids(I), savepar(J)%i, &
          TRANSPOSE(CONDENSED_MASS_VBS/Gtime%FSAVE_INTERVAL), start=(/1,1,GTIME%ind_netcdf/), count=(/NVBS,n_bins_par,1/)))

      if (savepar(j)%name == 'DEPOSITED_PAR_COMP'    ) call handler(__LINE__, nf90_put_var(ncfile_ids(I), savepar(J)%i, &
          Depos_composition, start=(/1,GTIME%ind_netcdf/), count=(/vapours%n_cond_tot/)))

      if (savepar(j)%name == 'PARTICLE_LOSS_RATE'    ) call handler(__LINE__, nf90_put_var(ncfile_ids(I), savepar(J)%i, &
          Lossrate, start=(/1,GTIME%ind_netcdf/), count=(/n_bins_par/)))

      if (savepar(j)%name == 'MASS_FLUX_ON_WALLS'    ) call handler(__LINE__, nf90_put_var(ncfile_ids(I), savepar(J)%i, &
          VAP_DEP_MASS_WALLS/Gtime%FSAVE_INTERVAL, start=(/1,GTIME%ind_netcdf/), count=(/vapours%n_cond_org/)))
          ! vapours%molar_mass(1:vapours%n_cond_org)*c_org_wall, start=(/1,GTIME%ind_netcdf/), count=(/vapours%n_cond_org/)))

    end do
    ! CONDENSED_MASS is zeroed after the mean flux is written to file
    if (Condensation) CONDENSED_MASS = 0d0
    if (Condensation) CONDENSED_MASS_VBS = 0d0
    if (Chem_Deposition) VAP_DEP_MASS_WALLS = 0d0
  end if

  ! Advance netcdf index by one
  GTIME%ind_netcdf = GTIME%ind_netcdf + 1

END SUBROUTINE SAVE_GASES

! Closes netcdf-files at the end of the run
subroutine CLOSE_FILES(filename)
  IMPLICIT NONE
  CHARACTER(LEN=*), INTENT(IN) :: filename
  INTEGER :: I,ioi
  DO I=1,N_FILES
    if (I == 3.and..not.Aerosol_flag) cycle
    call handler(__LINE__, nf90_close(ncfile_ids(I)))
  END DO
  DO I=1,N_FILES
    if (I == 3.and..not.Aerosol_flag) cycle
    ioi = RENAME(TRIM(ncfile_names(I))//'.tmp', TRIM(ncfile_names(I))//'.nc')
    if (ioi /= 0) print*, 'Error while copying final file.'
  END DO
  Write(*,FMT_INTRMT) 'Outputfiles in '//TRIM(filename)//' closed. '
  ! Write(*,FMT_LEND)

end subroutine CLOSE_FILES

! handler for NETCDF-calls, monitors errors
subroutine HANDLER(line, status)
  IMPLICIT NONE
  INTEGER, intent(in) :: line, status
  if (status /= nf90_noerr) then
    write(*,'(a, i0, a, a)') 'On line ', line, ': Error writing netcdf, code: ',NF90_STRERROR(status)
    stop
  end if
end subroutine HANDLER

! Basically a dictionary for unit names in MODS. Relates input units to output units
PURE CHARACTER(10) FUNCTION UNITS(unit)
  CHARACTER(*), intent(in) :: unit
  IF (UNIT == '#'   ) UNITS = '[1/cm^3]'
  IF (UNIT == 'ppm' ) UNITS = '[1/cm^3]'
  IF (UNIT == 'ppb' ) UNITS = '[1/cm^3]'
  IF (UNIT == 'ppt' ) UNITS = '[1/cm^3]'
  IF (UNIT == 'ppq' ) UNITS = '[1/cm^3]'
  IF (UNIT == 'Pa'  ) UNITS = '[Pa]'
  IF (UNIT == 'hPa' ) UNITS = '[Pa]'
  IF (UNIT == 'kPa' ) UNITS = '[Pa]'
  IF (UNIT == 'mbar') UNITS = '[Pa]'
  IF (UNIT == 'atm' ) UNITS = '[Pa]'
  IF (UNIT == 'K'   ) UNITS = '[K]'
  IF (UNIT == 'C'   ) UNITS = '[K]'

END FUNCTION UNITS

SUBROUTINE INITIALIZE_WITH_LAST(c_pp, n_c, chemconc)
    IMPLICIT NONE
    real(dp) :: c_pp(:,:)
    real(dp), ALLOCATABLE :: c_pp_tp(:,:)
    real(dp) :: n_c(:)
    real(dp) :: chemconc(:)
    integer :: file_id,i,varid, len_time, dimid, datarow = 0
    ALLOCATE(c_pp_tp(size(c_pp, dim=2), size(c_pp, dim=1) ))
    if (Chemistry_flag) THEN
        print FMT_MSG, 'Initializing chemistry with '//TRIM(INITIALIZE_WITH)//'. This takes a while.'
        call handler(__LINE__, nf90_open(TRIM(INITIALIZE_WITH)//'/Chemistry.nc', NF90_NOWRITE, file_id) )

    if (INITIALIZE_FROM>0) THEN
        datarow = INITIALIZE_FROM
    ELSE
        call handler(__LINE__, nf90_inq_dimid(file_id, 'time', dimid) )
        call handler(__LINE__, nf90_inquire_dimension(file_id, dimid, len=len_time) )
        datarow = len_time
    END if
    do i = 1,size(chemconc)
        call handler(__LINE__, nf90_inq_varid(file_id, TRIM(SPC_NAMES(i)), varid))
        call handler(__LINE__, nf90_get_var(file_id, varid, chemconc(i), start = ([( datarow )]) ))
    end do

    call handler(__LINE__, nf90_close(file_id) ) ! close netcdf dataset
    END IF

if (Aerosol_flag) THEN
    print FMT_MSG, 'Initializing aerosols with '//TRIM(INITIALIZE_WITH)//'. This takes a while.'
    call handler(__LINE__, nf90_open(TRIM(INITIALIZE_WITH)//'/Particles.nc', NF90_NOWRITE, file_id) )

    ! Get the last index if no index was given
    if (INITIALIZE_FROM>0) THEN
        datarow = INITIALIZE_FROM
    ELSE
        call handler(__LINE__, nf90_inq_dimid(file_id, 'time', dimid) )
        call handler(__LINE__, nf90_inquire_dimension(file_id, dimid, len=len_time) )
        datarow = len_time
    END if

    n_c = 0d0
    call handler(__LINE__, nf90_inq_varid(file_id, 'NUMBER_CONCENTRATION', varid))
    call handler(__LINE__, nf90_get_var(file_id, varid, n_c,  start=(/1,datarow/), count=(/n_bins_par/)))
    n_c = n_c*1d6

    call handler(__LINE__, nf90_inq_varid(file_id, 'PARTICLE_COMPOSITION', varid))
    ! TRANSPOSE(get_composition()), (/1,1,GTIME%ind_netcdf/), (/vapours%n_cond_tot,n_bins_par,1/))


    call handler(__LINE__, nf90_get_var(file_id, varid, c_pp_tp, start=(/1,1,datarow/), count=(/size(c_pp, 2), n_bins_par/)))
    c_pp = TRANSPOSE(c_pp_tp)
    call handler(__LINE__, nf90_close(file_id) ) ! close netcdf dataset
END IF

END SUBROUTINE INITIALIZE_WITH_LAST

SUBROUTINE WRITE_GAS_BINARY_DUMP(dir)
    CHARACTER(LEN=*) :: dir
    INQUIRE(iolength=reclen) TIMESERIES

    call system('rm -f '//trim(dir)//'/TIMESERIES.r4')
    open (unit=1,file=trim(dir)//'/TIMESERIES.r4',form='unformatted',access='direct',recl=reclen,STATUS='new')
    write (1,rec=1) TIMESERIES
    close(1)
END SUBROUTINE WRITE_GAS_BINARY_DUMP

END MODULE OUTPUT
