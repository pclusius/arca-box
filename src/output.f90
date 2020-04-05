MODULE OUTPUT
use netcdf
use CONSTANTS
use INPUT
use second_Monitor
use AUXILLARIES
USE Aerosol_auxillaries, only: vapour_ambient
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
CHARACTER(200) :: ncfile_names(N_FILES) = (['General  ', 'Chemistry', 'Particle '])

INTEGER, allocatable        :: shifter_ind(:)
INTEGER, allocatable        :: multipl_ind(:)
INTEGER, allocatable        :: mods_ind(:)
INTEGER, allocatable        :: chem_ind(:)
INTEGER, allocatable        :: par_ind(:)

INTEGER :: dtime_id
INTEGER :: dbins_id
INTEGER :: dcond_id
INTEGER :: dstring_id
INTEGER :: dconstant_id
INTEGER :: timearr_id
INTEGER :: hrsarr_id
INTEGER :: gJ_out_NH3_id
INTEGER :: gJ_out_DMA_id
INTEGER :: gJ_out_SUM_id
INTEGER :: gRES_BASE
INTEGER :: gRES_J
INTEGER :: n_condensables=0
public :: OPEN_FILES, SAVE_GASES,CLOSE_FILES

type(parsave) :: parbuf(20)
type(parsave), ALLOCATABLE :: savepar(:)

CONTAINS

  ! --------------------------------------------------------------------------------------------------------------------
  ! CREATE (OR OVERWRITE OLD) THREE NETCDF-FILES TO STORE OUTPUT. PARTICLES IS STILL VERY ROUGH SINCE WE DON'T HAVE THEM
  ! YET, BUT GENERAL AND CHEMISTRY ARE THERE ALREADY.
  ! --------------------------------------------------------------------------------------------------------------------
  SUBROUTINE OPEN_FILES(filename, Description, MODS, CH_GAS, vapours, current_psd)
    IMPLICIT NONE

    CHARACTER(LEN=*), INTENT(IN)    :: filename
    CHARACTER(*), INTENT(IN)        :: Description
    TYPE(input_mod),INTENT(INOUT)   :: MODS(:)
    TYPE(vapour_ambient),INTENT(IN) :: vapours
    TYPE(psd),INTENT(IN)            :: current_PSD
    REAL(dp), INTENT(IN)            :: CH_GAS(:)
    CHARACTER(255)                  :: PROGRAM_NAME
    INTEGER                         :: i,j,k,ioi,lenD
    INTEGER, PARAMETER              :: textdim = len(SPC_NAMES(1))
    CHARACTER(textdim), ALLOCATABLE :: COND_NAMES(:)

    ! ---------------------------
    ! ------dimension bits-------
    ! time | bins | condensables | name
    !  1      2          4           8
    ! eg. number concentration:
    ! time | bins | condensables
    !  1      1          0         => 011 => parbuf(i)%d = 3
    i=1
    parbuf(i)%name = 'CONDENSABLES'          ; parbuf(i)%u = '[]'         ; parbuf(i)%d = -12; parbuf(i)%type = NF90_CHAR  ; i=i+1
    parbuf(i)%name = 'NUMBER_CONCENTRATION'  ; parbuf(i)%u = '[1/cm^3]'   ; parbuf(i)%d = 3 ; parbuf(i)%type = NF90_DOUBLE ; i=i+1
    parbuf(i)%name = 'DIAMETER'              ; parbuf(i)%u = '[m]'        ; parbuf(i)%d = 3 ; parbuf(i)%type = NF90_DOUBLE ; i=i+1
    parbuf(i)%name = 'DRY_DIAMETER'          ; parbuf(i)%u = '[m]'        ; parbuf(i)%d = 3 ; parbuf(i)%type = NF90_DOUBLE ; i=i+1
    parbuf(i)%name = 'ORIGINAL_DRY_DIAMETER' ; parbuf(i)%u = '[m]'        ; parbuf(i)%d = 3 ; parbuf(i)%type = NF90_DOUBLE ; i=i+1
    parbuf(i)%name = 'GROWTH_RATE'           ; parbuf(i)%u = '[m/s]'      ; parbuf(i)%d = 3 ; parbuf(i)%type = NF90_DOUBLE ; i=i+1
    parbuf(i)%name = 'CORE_VOLUME'           ; parbuf(i)%u = '[m^3]'      ; parbuf(i)%d = 3 ; parbuf(i)%type = NF90_DOUBLE ; i=i+1
    parbuf(i)%name = 'MASS'                  ; parbuf(i)%u = '[kg]'       ; parbuf(i)%d = 3 ; parbuf(i)%type = NF90_DOUBLE ; i=i+1
    parbuf(i)%name = 'SURFACE_TENSION'       ; parbuf(i)%u = '[N/m]'      ; parbuf(i)%d = 4 ; parbuf(i)%type = NF90_DOUBLE ; i=i+1
    parbuf(i)%name = 'VAPOR_CONCENTRATION'   ; parbuf(i)%u = '[#/m^3]'    ; parbuf(i)%d = 5 ; parbuf(i)%type = NF90_DOUBLE ; i=i+1
    parbuf(i)%name = 'PARTICLE_COMPOSITION'  ; parbuf(i)%u = '[kg/m^3]'   ; parbuf(i)%d = 7 ; parbuf(i)%type = NF90_DOUBLE ; i=i+1
!    parbuf(i)%name = 'VOLUME_CONCENTRATION'  ; parbuf(i)%u = '[um^3/m^3]' ; parbuf(i)%d = 7 ; parbuf(i)%type = NF90_DOUBLE ; i=i+1

    allocate(savepar(i-1))
    savepar = parbuf(1:i-1)

    do j = 1,size(vapours%vapour_names)
      k = IndexFromName( vapours%vapour_names(j), SPC_NAMES )
      if (k>0) n_condensables=n_condensables+1
    end do
    ALLOCATE(COND_NAMES(n_condensables))
    i=1
    do j = 1,size(vapours%vapour_names)
      k = IndexFromName( vapours%vapour_names(j), SPC_NAMES )
      if (k>0) THEN
        COND_NAMES(i) = vapours%vapour_names(j)
        i = i + 1
      END IF
    end do

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

      ncfile_names(I) = trim(filename)//'/'//TRIM(ncfile_names(I))//'.nc'

      !Clearing file; Opening file. Overwrites
      open(720+I, FILE=ncfile_names(I), iostat = ioi)
      CALL handle_file_io(ioi, ncfile_names(I), 'Terminating when trying to open netCDF-file, does the CASE and RUN directory exist')
      close(720+I)

      ! Added compression for particle.nc, so we need to use netCDF4-file. Here used in classic mode
      ! call handler( nf90_create(ncfile_names(I), NF90_NETCDF4, ncfile_ids(I)) )
      call handler( nf90_create(ncfile_names(I), IOR(NF90_NETCDF4, NF90_CLASSIC_MODEL), ncfile_ids(I)) )

      ! Defining dimensions: time(unlimited), size sections, vapor_species
      call handler(nf90_def_dim(ncfile_ids(I), "time",NINT(MODELTIME%SIM_TIME_S/MODELTIME%FSAVE_INTERVAL+1), dtime_id) )
      IF (I==3) call handler(nf90_def_dim(ncfile_ids(I), "string",textdim, dstring_id) )
      call handler(nf90_def_dim(ncfile_ids(I), "bins",n_bins_particle, dbins_id) )
      call handler(nf90_def_dim(ncfile_ids(I), "condensables",n_condensables, dcond_id) )
      call handler(nf90_def_dim(ncfile_ids(I), "Constant",1, dconstant_id) )

      !Create attributes for general stuff
      CALL get_command_argument(0, PROGRAM_NAME)
      call handler(nf90_put_att(ncfile_ids(I), NF90_GLOBAL, 'Information', '(c) Atmospheric modelling group 2019 and (c) Simugroup 2019 (ACDC)'))
      call handler(nf90_put_att(ncfile_ids(I), NF90_GLOBAL, 'Contact', 'michael.boy@helsinki.fi (Superbox), tinja.olenius@alumni.helsinki.fi (ACDC)'))
      call handler(nf90_put_att(ncfile_ids(I), NF90_GLOBAL, 'Software', 'Superbox 0.0.1'))
      call handler(nf90_put_att(ncfile_ids(I), NF90_GLOBAL, 'Package_Name:', TRIM(PROGRAM_NAME(3:))))
      call handler(nf90_put_att(ncfile_ids(I), NF90_GLOBAL, 'Notes', TRIM(Description)))
      call handler(nf90_put_att(ncfile_ids(I), NF90_GLOBAL, 'experiment', Fname_init))
      call handler(nf90_def_var(ncfile_ids(I), "Time_in_sec", NF90_DOUBLE, dtime_id, timearr_id))
      call handler(nf90_def_var(ncfile_ids(I), "time_in_hrs", NF90_DOUBLE, dtime_id, hrsarr_id))
      ! COMPRESSION
      call handler(nf90_def_var_deflate(ncfile_ids(I), timearr_id, shuff, compress, compression) )
      call handler(nf90_def_var_deflate(ncfile_ids(I), hrsarr_id, shuff, compress, compression) )
      call handler(nf90_put_att(ncfile_ids(I), timearr_id, 'unit' , '[s]'))
      call handler(nf90_put_att(ncfile_ids(I), hrsarr_id, 'unit' , '[hrs]'))

    END DO

  ALLOCATE(multipl_ind(size(MODS)))
  ALLOCATE(shifter_ind(size(MODS)))
  ALLOCATE(mods_ind(size(MODS)))
  ALLOCATE(chem_ind(size(CH_GAS)))
  ALLOCATE(par_ind(size(vapours%vapour_names)))


  I=1 ! GENERAL FILE
  do j = 1,size(MODS)
    IF ((TRIM(MODS(j)%name) /= '#')) THEN
      IF (j<LENV .or. INDRELAY_CH(j)>0) THEN
        call handler(nf90_def_var(ncfile_ids(I), TRIM(MODS(j)%name), NF90_DOUBLE, dtime_id, mods_ind(j)) )
        call handler(nf90_def_var_deflate(ncfile_ids(I), mods_ind(j), shuff, compress, compression) )
        call handler(nf90_put_att(ncfile_ids(I), mods_ind(j), 'unit' , TRIM(UNITS(MODS(I)%UNIT))))
      END IF

      IF ((ABS(MODS(j)%shift-0d0) > 1d-100) .or. (ABS(MODS(j)%multi-1d0) > 1d-100)) THEN
        call handler(nf90_def_var(ncfile_ids(I), TRIM(MODS(j)%name)//'_Multipl', NF90_DOUBLE, dtime_id, multipl_ind(j)))
        call handler(nf90_def_var(ncfile_ids(I), TRIM(MODS(j)%name)//'_Shifter', NF90_DOUBLE, dtime_id, shifter_ind(j)))

        call handler(nf90_def_var_deflate(ncfile_ids(I), multipl_ind(j), shuff, compress, compression) )
        call handler(nf90_def_var_deflate(ncfile_ids(I), shifter_ind(j), shuff, compress, compression) )

        call handler(nf90_put_att(ncfile_ids(I), multipl_ind(j), 'unit' , '[]'))
        call handler(nf90_put_att(ncfile_ids(I), shifter_ind(j), 'unit' , '[same]'))

      END IF
    END IF
  end do


  call handler(nf90_def_var(ncfile_ids(I), 'J_ACDC_NH3', NF90_DOUBLE, dtime_id, gJ_out_NH3_id))
  call handler(nf90_def_var(ncfile_ids(I), 'J_ACDC_DMA', NF90_DOUBLE, dtime_id, gJ_out_DMA_id))
  call handler(nf90_def_var(ncfile_ids(I), 'J_ACDC_SUM', NF90_DOUBLE, dtime_id, gJ_out_SUM_id))

  call handler(nf90_def_var_deflate(ncfile_ids(I), gJ_out_NH3_id,         shuff, compress, compression) )
  call handler(nf90_def_var_deflate(ncfile_ids(I), gJ_out_DMA_id,         shuff, compress, compression) )
  call handler(nf90_def_var_deflate(ncfile_ids(I), gJ_out_SUM_id,         shuff, compress, compression) )

  call handler(nf90_put_att(ncfile_ids(I), gJ_out_NH3_id, 'unit' , '[1/s/cm^3]'))
  call handler(nf90_put_att(ncfile_ids(I), gJ_out_DMA_id, 'unit' , '[1/s/cm^3]'))
  call handler(nf90_put_att(ncfile_ids(I), gJ_out_SUM_id, 'unit' , '[1/s/cm^3]'))

  IF (RESOLVE_BASE) THEN
    call handler(nf90_def_var(ncfile_ids(I), 'RESOLVED_BASE', NF90_DOUBLE, dtime_id, gRES_BASE))
    call handler(nf90_def_var_deflate(ncfile_ids(I), gRES_BASE, shuff, compress, compression) )
    call handler(nf90_put_att(ncfile_ids(I), gRES_BASE, 'unit' , '[1/cm^3]'))

    call handler(nf90_def_var(ncfile_ids(I), 'RESOLVED_J', NF90_DOUBLE, dtime_id, gRES_J))
    call handler(nf90_def_var_deflate(ncfile_ids(I), gRES_J, shuff, compress, compression) )
    call handler(nf90_put_att(ncfile_ids(I), gRES_J, 'unit' , '[1/s/cm^3]'))
  end if



  I=2 ! Chemical file
  do j = 1,size(CH_GAS)
      call handler(nf90_def_var(ncfile_ids(I), TRIM(SPC_NAMES(j)), NF90_DOUBLE, dtime_id, chem_ind(j)) )
      call handler(nf90_def_var_deflate(ncfile_ids(I), chem_ind(j), shuff, compress, compression) )
      call handler(nf90_put_att(ncfile_ids(I), chem_ind(j), 'unit' , '1/cm^3'))
  end do



  I=3 ! Particle file. Currently only condensibles are stored here. Particles are added when we get them.
  do j = 1,size(savepar)

    if (savepar(J)%d == 0) call handler(nf90_def_var(ncfile_ids(I), TRIM(  savepar(J)%name  ), savepar(J)%type, dconstant_id  , savepar(J)%i) )
    if (savepar(J)%d == 3) call handler(nf90_def_var(ncfile_ids(I), TRIM(  savepar(J)%name  ), savepar(J)%type, ([dbins_id, dtime_id])  , savepar(J)%i) )
    if (savepar(J)%d == 4) call handler(nf90_def_var(ncfile_ids(I), TRIM(  savepar(J)%name  ), savepar(J)%type, dcond_id  , savepar(J)%i) )
    if (savepar(J)%d == 5) call handler(nf90_def_var(ncfile_ids(I), TRIM(  savepar(J)%name  ), savepar(J)%type, ([dcond_id,dtime_id])  , savepar(J)%i) )
    if (savepar(J)%d == 7) call handler(nf90_def_var(ncfile_ids(I), TRIM(  savepar(J)%name  ), savepar(J)%type, ([dcond_id,dbins_id,dtime_id])  , savepar(J)%i) )
    if (savepar(J)%d == -12) call handler(nf90_def_var(ncfile_ids(I), TRIM(  savepar(J)%name  ), savepar(J)%type, ([dstring_id,dcond_id])  , savepar(J)%i) )
    call handler(nf90_def_var_deflate(ncfile_ids(I), savepar(J)%i, shuff, compress, compression) )
    call handler(nf90_put_att(ncfile_ids(I), savepar(J)%i, 'unit' , savepar(J)%u))
  end do

  do j = 1,size(vapours%vapour_names)
    k = IndexFromName( vapours%vapour_names(j), SPC_NAMES )
    if (k>0) THEN
      call handler(nf90_def_var(ncfile_ids(I), TRIM(  vapours%vapour_names(j)  ), NF90_DOUBLE, dtime_id, par_ind(j)) )
      call handler(nf90_def_var_deflate(ncfile_ids(I), par_ind(j), shuff, compress, compression) )
      call handler(nf90_put_att(ncfile_ids(I), par_ind(j), 'unit' , '1/cm^3'))
    end if
  end do


  ! Ending definition mode
  DO I=1,N_FILES
    call handler( nf90_enddef(ncfile_ids(I)))
  END DO

  I=1

  do j = 1,size(vapours%vapour_names)
    k = IndexFromName( vapours%vapour_names(j), SPC_NAMES )
    if (k>0) THEN
      call handler( nf90_put_var(ncfile_ids(3), savepar(1)%i, COND_NAMES))
      i=i+1
    END IF

  end do


  print FMT_LEND,

  ORIGINAL_TEMP(2)  = MODS(inm_TempK)
  ORIGINAL_press(2) = MODS(inm_pres)
  MODS(inm_TempK) = ORIGINAL_TEMP(1)
  MODS(inm_pres)  = ORIGINAL_press(1)

  ! Also save all settings to initfile. Use this file to rerun if necessary
  open(9129, file=filename//'/RUN_INI.conf', action='WRITE')
  write(9129,NML=NML_Path)! directories and test cases
  write(9129,NML=NML_Flag)! flags
  write(9129,NML=NML_TIME)! time related stuff
  write(9129,NML=NML_PARTICLE)! dmps_file information
  write(9129,NML=NML_ENV)! environmental information
  write(9129,NML=NML_MCM)! MCM_file information
  write(9129,NML=NML_MODS)! modification parameters
  write(9129,NML=NML_MISC)! misc input
  write(9129,NML=NML_VAP)! vapour input
  close(9129)

  MODS(inm_TempK) = ORIGINAL_TEMP(2)
  MODS(inm_pres)  = ORIGINAL_press(2)


  RETURN

END SUBROUTINE OPEN_FILES



! --------------------------------------------------------------------------------------------------------------------
! Here the input is written to netcdf-files. Again, particles still rudimentary
! --------------------------------------------------------------------------------------------------------------------
SUBROUTINE SAVE_GASES(TSTEP_CONC,MODS,CH_GAS,J_ACDC_NH3, J_ACDC_DMA, VAPOURS, current_PSD)
  IMPLICIT NONE
  type(input_mod), INTENT(in)     :: MODS(:)
  real(dp), INTENT(in)            :: TSTEP_CONC(:)
  real(dp), INTENT(in)            :: CH_GAS(:)
  real(dp), INTENT(in)            :: J_ACDC_NH3
  real(dp), INTENT(in)            :: J_ACDC_DMA
  TYPE(vapour_ambient),INTENT(IN) :: vapours
  TYPE(psd),INTENT(IN) :: current_PSD
  INTEGER                         :: i,j,k

  DO I = 1,N_FILES
    call handler( nf90_put_var(ncfile_ids(I), timearr_id, MODELTIME%sec, (/MODELTIME%ind_netcdf/) ))
    call handler( nf90_put_var(ncfile_ids(I), hrsarr_id, MODELTIME%hrs, (/MODELTIME%ind_netcdf/) ))
  END DO

  I=1

  do j = 1,size(MODS)
    IF ((TRIM(MODS(j)%name) /= '#')) THEN
      IF (j<LENV .or. INDRELAY_CH(j)>0) THEN
        call handler( nf90_put_var(ncfile_ids(I), mods_ind(j), TSTEP_CONC(j), (/MODELTIME%ind_netcdf/)) )
      END IF
      IF ((ABS(MODS(j)%shift-0d0) > 1d-100) .or. (ABS(MODS(j)%multi-1d0) > 1d-100)) THEN
        call handler(nf90_put_var(ncfile_ids(I), multipl_ind(j), MODS(j)%multi, (/MODELTIME%ind_netcdf/) ) )
        call handler(nf90_put_var(ncfile_ids(I), shifter_ind(j), MODS(j)%shift,  (/MODELTIME%ind_netcdf/) ) )
      END IF
    END IF
  end do

  call handler( nf90_put_var(ncfile_ids(I), gJ_out_NH3_id, J_ACDC_NH3, (/MODELTIME%ind_netcdf/)) )
  call handler( nf90_put_var(ncfile_ids(I), gJ_out_DMA_id, J_ACDC_DMA, (/MODELTIME%ind_netcdf/)) )
  call handler( nf90_put_var(ncfile_ids(I), gJ_out_SUM_id, J_ACDC_DMA+J_ACDC_NH3, (/MODELTIME%ind_netcdf/)) )
  IF (RESOLVE_BASE) THEN
    call handler( nf90_put_var(ncfile_ids(I), gRES_BASE, RESOLVED_BASE, (/MODELTIME%ind_netcdf/)) )
    call handler( nf90_put_var(ncfile_ids(I), gRES_J,    RESOLVED_J,     (/MODELTIME%ind_netcdf/)) )
  END IF

  I=2 ! Chemical file
  do j = 1,size(CH_GAS)
    call handler( nf90_put_var(ncfile_ids(I), chem_ind(j), CH_GAS(j), (/MODELTIME%ind_netcdf/)) )
  end do



  I=3 ! Particle file
  if (Aerosol_flag) THEN
    do j = 1,size(vapours%vapour_names)
      k = IndexFromName( vapours%vapour_names(j),   SPC_NAMES )
      if (k>0) then
        call handler(nf90_put_var(ncfile_ids(I), par_ind(j), CH_GAS(k), (/MODELTIME%ind_netcdf/)) )
      end if
    end do
    do j = 1,size(savepar)
      ! if (parbuf(i)%name == 'CONDENSABLES'          ) call handler(nf90_put_var(ncfile_ids(I), savepar(J)%i, current_PSD%conc_fs, start=(/1,MODELTIME%ind_netcdf/), count=(/n_bins_particle/)))
      if (savepar(j)%name == 'NUMBER_CONCENTRATION'  ) call handler(nf90_put_var(ncfile_ids(I), savepar(J)%i, current_PSD%conc_fs, start=(/1,MODELTIME%ind_netcdf/), count=(/n_bins_particle/)))
      if (savepar(j)%name == 'DIAMETER'              ) call handler(nf90_put_var(ncfile_ids(I), savepar(J)%i, current_PSD%diameter_fs, start=(/1,MODELTIME%ind_netcdf/), count=(/n_bins_particle/)))
      ! if (savepar(j)%name == 'DRY_DIAMETER'          ) call handler(nf90_put_var(ncfile_ids(I), savepar(J)%i, current_PSD%dp_dry_fs, start=(/1,MODELTIME%ind_netcdf/), count=(/n_bins_particle/)))
      ! if (parbuf(i)%name == 'ORIGINAL_DRY_DIAMETER' ) call handler(nf90_put_var(ncfile_ids(I), savepar(J)%i, current_PSD%original_dry_radius, start=(/1,MODELTIME%ind_netcdf/), count=(/n_bins_particle/)))
      ! if (parbuf(i)%name == 'GROWTH_RATE'           ) call handler(nf90_put_var(ncfile_ids(I), savepar(J)%i, current_PSD%conc_fs, start=(/1,MODELTIME%ind_netcdf/), count=(/n_bins_particle/)))
      if (savepar(j)%name == 'CORE_VOLUME'           ) call handler(nf90_put_var(ncfile_ids(I), savepar(J)%i, current_PSD%volume_fs, start=(/1,MODELTIME%ind_netcdf/), count=(/n_bins_particle/)))
      if (savepar(j)%name == 'MASS'                  ) call handler(nf90_put_var(ncfile_ids(I), savepar(J)%i, current_PSD%particle_mass_fs, start=(/1,MODELTIME%ind_netcdf/), count=(/n_bins_particle/)))
      if (savepar(j)%name == 'PARTICLE_COMPOSITION'  ) call handler(nf90_put_var(ncfile_ids(I), savepar(J)%i, current_PSD%composition_fs, start=(/1,1,MODELTIME%ind_netcdf/), count=(/n_condensables,n_bins_particle/)))
      ! if (parbuf(i)%name == 'SURFACE_TENSION'       ) call handler(nf90_put_var(ncfile_ids(I), savepar(J)%i, current_PSD%conc_fs, start=(/1,MODELTIME%ind_netcdf/), count=(/n_bins_particle/)))
      ! if (parbuf(i)%name == 'VAPOR_CONCENTRATION'   ) call handler(nf90_put_var(ncfile_ids(I), savepar(J)%i, current_PSD%conc_fs, start=(/1,MODELTIME%ind_netcdf/), count=(/n_bins_particle/)))

      ! if (savepar(J)%d == 0) call handler(nf90_put_var(ncfile_ids(I), TRIM(  savepar(J)%name  ), savepar(J)%type, dconstant_id  , savepar(J)%i) )
      ! if (savepar(J)%d == 3) call handler(nf90_put_var(ncfile_ids(I), savepar(J)%i, current_PSD%conc_fs, start=(/1,MODELTIME%ind_netcdf/), count=(/n_bins_particle/)))
      ! if (savepar(J)%d == 4) call handler(nf90_put_var(ncfile_ids(I), TRIM(  savepar(J)%name  ), savepar(J)%type, dcond_id  , savepar(J)%i) )
      ! if (savepar(J)%d == 5) call handler(nf90_put_var(ncfile_ids(I), TRIM(  savepar(J)%name  ), savepar(J)%type, ([dcond_id,dtime_id])  , savepar(J)%i) )
      ! if (savepar(J)%d == 7) call handler(nf90_put_var(ncfile_ids(I), TRIM(  savepar(J)%name  ), savepar(J)%type, ([dcond_id,dbins_id,dtime_id])  , savepar(J)%i) )
      ! if (savepar(J)%d == -12) call handler(nf90_put_var(ncfile_ids(I), TRIM(  savepar(J)%name  ), savepar(J)%type, ([dstring_id,dcond_id])  , savepar(J)%i) )
    end do
  end if


END SUBROUTINE SAVE_GASES


! Closes netcdf-files at the end of the run
subroutine CLOSE_FILES()
  IMPLICIT NONE
  INTEGER :: I
  DO I=1,N_FILES
    call handler( nf90_close(ncfile_ids(I)))
  END DO

  Write(*,FMT_MSG) 'Outputfiles closed.'
  Write(*,FMT_LEND)

end subroutine CLOSE_FILES

! handler for NETCDF-calls, monitors errors
subroutine HANDLER(status)
  IMPLICIT NONE
  INTEGER, intent(in) :: status
  if (status /= nf90_noerr) then
    write(*,*) 'Error writing netcdf, code: ',NF90_STRERROR(status)
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



END MODULE OUTPUT
