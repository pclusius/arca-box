MODULE OUTPUT
use netcdf
use CONSTANTS
use INPUT
use second_Monitor
use AUXILLARIES
USE PSD_scheme, ONLY: get_composition, get_dp, get_volume, get_mass, get_conc

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
INTEGER :: gJ_out_TOT_id
INTEGER :: gRES_BASE
INTEGER :: gRES_J
INTEGER :: gRES_F

public :: OPEN_FILES, SAVE_GASES,CLOSE_FILES,INITIALIZE_WITH_LAST

type(parsave) :: parbuf(20)
type(parsave), ALLOCATABLE :: savepar(:)

CONTAINS

  ! --------------------------------------------------------------------------------------------------------------------
  ! CREATE (OR OVERWRITE OLD) THREE NETCDF-FILES TO STORE OUTPUT. PARTICLES IS STILL VERY ROUGH SINCE WE DON'T HAVE THEM
  ! YET, BUT GENERAL AND CHEMISTRY ARE THERE ALREADY.
  ! --------------------------------------------------------------------------------------------------------------------
  SUBROUTINE OPEN_FILES(filename, Description, currentChemistry, MODS, CH_GAS, vapours)
    IMPLICIT NONE

    CHARACTER(LEN=*), INTENT(IN)    :: filename
    CHARACTER(LEN=*), INTENT(IN)    :: currentChemistry
    CHARACTER(*), INTENT(IN)        :: Description
    TYPE(input_mod),INTENT(INOUT)   :: MODS(:)
    TYPE(vapour_ambient),INTENT(IN) :: vapours
    ! TYPE(psd),INTENT(IN)            :: current_PSD
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
    parbuf(i)%name = 'INPUT_CONCENTRATION'   ; parbuf(i)%u = '[1/cm^3]'   ; parbuf(i)%d = 3 ; parbuf(i)%type = NF90_DOUBLE ; i=i+1
    parbuf(i)%name = 'DIAMETER'              ; parbuf(i)%u = '[m]'        ; parbuf(i)%d = 3 ; parbuf(i)%type = NF90_DOUBLE ; i=i+1
    ! parbuf(i)%name = 'DRY_DIAMETER'          ; parbuf(i)%u = '[m]'        ; parbuf(i)%d = 3 ; parbuf(i)%type = NF90_DOUBLE ; i=i+1
    parbuf(i)%name = 'GROWTH_RATE'           ; parbuf(i)%u = '[nm/h]'     ; parbuf(i)%d = 3 ; parbuf(i)%type = NF90_DOUBLE ; i=i+1
    ! parbuf(i)%name = 'CORE_VOLUME'           ; parbuf(i)%u = '[m^3]'      ; parbuf(i)%d = 3 ; parbuf(i)%type = NF90_DOUBLE ; i=i+1
    parbuf(i)%name = 'MASS'                  ; parbuf(i)%u = '[kg]'       ; parbuf(i)%d = 3 ; parbuf(i)%type = NF90_DOUBLE ; i=i+1
    ! parbuf(i)%name = 'SURFACE_TENSION'       ; parbuf(i)%u = '[N/m]'      ; parbuf(i)%d = 4 ; parbuf(i)%type = NF90_DOUBLE ; i=i+1
    ! parbuf(i)%name = 'VAPOR_CONCENTRATION'   ; parbuf(i)%u = '[#/m^3]'    ; parbuf(i)%d = 5 ; parbuf(i)%type = NF90_DOUBLE ; i=i+1
    parbuf(i)%name = 'PARTICLE_COMPOSITION'  ; parbuf(i)%u = '[kg/particle]'   ; parbuf(i)%d = 7 ; parbuf(i)%type = NF90_DOUBLE ; i=i+1
!    parbuf(i)%name = 'VOLUME_CONCENTRATION'  ; parbuf(i)%u = '[um^3/m^3]' ; parbuf(i)%d = 7 ; parbuf(i)%type = NF90_DOUBLE ; i=i+1

    allocate(savepar(i-1))
    savepar = parbuf(1:i-1)

    if (Aerosol_flag) THEN
        ALLOCATE(COND_NAMES(vapours%n_condtot))
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

      ncfile_names(I) = trim(filename)//'/'//TRIM(ncfile_names(I))

      !Clearing file; Opening file. Overwrites
      open(720+I, FILE=TRIM(ncfile_names(I))//'.tmp', iostat = ioi)
      CALL handle_file_io(ioi, TRIM(ncfile_names(I))//'.tmp', 'Terminating when trying to open netCDF-file, does the CASE and RUN directory exist')
      close(720+I)

      ! Added compression for particles.nc, so we need to use netCDF4-file. Here used in classic mode
      ! call handler(__LINE__, nf90_create(TRIM(ncfile_names(I))//'.tmp', NF90_NETCDF4, ncfile_ids(I)) )
      call handler(__LINE__, nf90_create(TRIM(ncfile_names(I))//'.tmp', IOR(NF90_NETCDF4, NF90_CLASSIC_MODEL), ncfile_ids(I)) )

      ! Defining dimensions: time(unlimited), size sections, vapor_species
      ! call handler(__LINE__, nf90_def_dim(ncfile_ids(I), "time",NINT(GTIME%SIM_TIME_S/GTIME%FSAVE_INTERVAL+1), dtime_id) )
      call handler(__LINE__, nf90_def_dim(ncfile_ids(I), "time",NF90_UNLIMITED, dtime_id) )
      IF (I==3) call handler(__LINE__, nf90_def_dim(ncfile_ids(I), "string",textdim, dstring_id) )
      IF (I>1)  call handler(__LINE__, nf90_def_dim(ncfile_ids(I), "bins",n_bins_par, dbins_id) )
      IF ((I == 3) .and. Aerosol_flag) call handler(__LINE__, nf90_def_dim(ncfile_ids(I), "condensables",vapours%n_condtot, dcond_id) )
      call handler(__LINE__, nf90_def_dim(ncfile_ids(I), "Constant",1, dconstant_id) )

      !Create attributes for general stuff
      CALL get_command_argument(0, PROGRAM_NAME)
      call handler(__LINE__, nf90_put_att(ncfile_ids(I), NF90_GLOBAL, 'Information', '(c) Multiscale Modelling Group and (c) Computational Aerosol Physics Group (ACDC)'))
      call handler(__LINE__, nf90_put_att(ncfile_ids(I), NF90_GLOBAL, 'Contact', 'arca@helsinki.fi (ARCA box), tinja.olenius@alumni.helsinki.fi (ACDC)'))
      call handler(__LINE__, nf90_put_att(ncfile_ids(I), NF90_GLOBAL, 'Software', 'ARCA box 1.0.1'))
      call handler(__LINE__, nf90_put_att(ncfile_ids(I), NF90_GLOBAL, 'Package_Name', TRIM(PROGRAM_NAME(3:))))
      call handler(__LINE__, nf90_put_att(ncfile_ids(I), NF90_GLOBAL, 'Chemistry_module', TRIM(CurrentChemistry)))
      call handler(__LINE__, nf90_put_att(ncfile_ids(I), NF90_GLOBAL, 'Notes', TRIM(Description)))
      call handler(__LINE__, nf90_put_att(ncfile_ids(I), NF90_GLOBAL, 'experiment', filename))
      call handler(__LINE__, nf90_put_att(ncfile_ids(I), NF90_GLOBAL, 'INITFILE', Fname_init))
      call handler(__LINE__, nf90_def_var(ncfile_ids(I), "TIME_IN_SEC", NF90_DOUBLE, dtime_id, timearr_id))
      call handler(__LINE__, nf90_def_var(ncfile_ids(I), "TIME_IN_HRS", NF90_DOUBLE, dtime_id, hrsarr_id))
      ! COMPRESSION
      call handler(__LINE__, nf90_def_var_deflate(ncfile_ids(I), timearr_id, shuff, compress, compression) )
      call handler(__LINE__, nf90_def_var_deflate(ncfile_ids(I), hrsarr_id, shuff, compress, compression) )
      call handler(__LINE__, nf90_put_att(ncfile_ids(I), timearr_id, 'unit' , '[s]'))
      call handler(__LINE__, nf90_put_att(ncfile_ids(I), hrsarr_id, 'unit' , '[hrs]'))

    END DO

  ALLOCATE(multipl_ind(size(MODS)))
  ALLOCATE(shifter_ind(size(MODS)))
  ALLOCATE(mods_ind(size(MODS)))
  ALLOCATE(chem_ind(size(CH_GAS)))
  ALLOCATE(par_ind(size(vapours%vapour_names)))


  I=1 ! GENERAL FILE
  do j = 1,size(MODS)
    IF (TRIM(MODS(j)%name) /= '#') THEN
      IF (j<LENV .or. INDRELAY_CH(j)>0) THEN
        call handler(__LINE__, nf90_def_var(ncfile_ids(I), TRIM(MODS(j)%name), NF90_DOUBLE, dtime_id, mods_ind(j)) )
        call handler(__LINE__, nf90_def_var_deflate(ncfile_ids(I), mods_ind(j), shuff, compress, compression) )
        call handler(__LINE__, nf90_put_att(ncfile_ids(I), mods_ind(j), 'unit' , TRIM(UNITS(MODS(I)%UNIT))))
      END IF

      IF ((ABS(MODS(j)%shift-0d0) > 1d-100) .or. (ABS(MODS(j)%multi-1d0) > 1d-100)) THEN
        call handler(__LINE__, nf90_def_var(ncfile_ids(I), TRIM(MODS(j)%name)//'_Multipl', NF90_DOUBLE, dtime_id, multipl_ind(j)))
        call handler(__LINE__, nf90_def_var(ncfile_ids(I), TRIM(MODS(j)%name)//'_Shifter', NF90_DOUBLE, dtime_id, shifter_ind(j)))

        call handler(__LINE__, nf90_def_var_deflate(ncfile_ids(I), multipl_ind(j), shuff, compress, compression) )
        call handler(__LINE__, nf90_def_var_deflate(ncfile_ids(I), shifter_ind(j), shuff, compress, compression) )

        call handler(__LINE__, nf90_put_att(ncfile_ids(I), multipl_ind(j), 'unit' , '[]'))
        call handler(__LINE__, nf90_put_att(ncfile_ids(I), shifter_ind(j), 'unit' , '[same]'))

      END IF
    END IF
  end do


  call handler(__LINE__, nf90_def_var(ncfile_ids(I), 'J_ACDC_NH3_CM3', NF90_DOUBLE, dtime_id, gJ_out_NH3_id))
  call handler(__LINE__, nf90_def_var(ncfile_ids(I), 'J_ACDC_DMA_CM3', NF90_DOUBLE, dtime_id, gJ_out_DMA_id))
  call handler(__LINE__, nf90_def_var(ncfile_ids(I), 'J_ACDC_SUM_CM3', NF90_DOUBLE, dtime_id, gJ_out_SUM_id))
  call handler(__LINE__, nf90_def_var(ncfile_ids(I), 'J_TOTAL_CM3', NF90_DOUBLE, dtime_id, gJ_out_TOT_id))

  call handler(__LINE__, nf90_def_var_deflate(ncfile_ids(I), gJ_out_NH3_id,shuff, compress, compression) )
  call handler(__LINE__, nf90_def_var_deflate(ncfile_ids(I), gJ_out_DMA_id,shuff, compress, compression) )
  call handler(__LINE__, nf90_def_var_deflate(ncfile_ids(I), gJ_out_SUM_id,shuff, compress, compression) )
  call handler(__LINE__, nf90_def_var_deflate(ncfile_ids(I), gJ_out_TOT_id,shuff, compress, compression) )

  call handler(__LINE__, nf90_put_att(ncfile_ids(I), gJ_out_NH3_id, 'unit' , '[1/s/cm^3]'))
  call handler(__LINE__, nf90_put_att(ncfile_ids(I), gJ_out_DMA_id, 'unit' , '[1/s/cm^3]'))
  call handler(__LINE__, nf90_put_att(ncfile_ids(I), gJ_out_SUM_id, 'unit' , '[1/s/cm^3]'))
  call handler(__LINE__, nf90_put_att(ncfile_ids(I), gJ_out_TOT_id, 'unit' , '[1/s/cm^3]'))

  IF (RESOLVE_BASE) THEN
    call handler(__LINE__, nf90_def_var(ncfile_ids(I), 'RESOLVED_BASE', NF90_DOUBLE, dtime_id, gRES_BASE))
    call handler(__LINE__, nf90_def_var_deflate(ncfile_ids(I), gRES_BASE, shuff, compress, compression) )
    call handler(__LINE__, nf90_put_att(ncfile_ids(I), gRES_BASE, 'unit' , '[1/cm^3]'))

    call handler(__LINE__, nf90_def_var(ncfile_ids(I), 'RESOLVED_J', NF90_DOUBLE, dtime_id, gRES_J))
    call handler(__LINE__, nf90_def_var_deflate(ncfile_ids(I), gRES_J, shuff, compress, compression) )
    call handler(__LINE__, nf90_put_att(ncfile_ids(I), gRES_J, 'unit' , '[1/s/cm^3]'))

    call handler(__LINE__, nf90_def_var(ncfile_ids(I), 'RESOLVED_J_FACTR', NF90_DOUBLE, dtime_id, gRES_F))
    call handler(__LINE__, nf90_def_var_deflate(ncfile_ids(I), gRES_F, shuff, compress, compression) )
    call handler(__LINE__, nf90_put_att(ncfile_ids(I), gRES_F, 'unit' , '[]'))
  end if



  I=2 ! Chemical file
  do j = 1,size(CH_GAS)
      call handler(__LINE__, nf90_def_var(ncfile_ids(I), TRIM(SPC_NAMES(j)), NF90_DOUBLE, dtime_id, chem_ind(j)) )
      call handler(__LINE__, nf90_def_var_deflate(ncfile_ids(I), chem_ind(j), shuff, compress, compression) )
      call handler(__LINE__, nf90_put_att(ncfile_ids(I), chem_ind(j), 'unit' , '1/cm^3'))
      call handler(__LINE__, nf90_put_att(ncfile_ids(I), chem_ind(j), 'condensing' , 1))
  end do



  I=3 ! Particle file.
  if (Aerosol_flag) THEN

      do j = 1,size(savepar)

        if (savepar(J)%d == 0) call handler(__LINE__, nf90_def_var(ncfile_ids(I), TRIM(  savepar(J)%name  ), savepar(J)%type, dconstant_id  , savepar(J)%i) )
        if (savepar(J)%d == 3) call handler(__LINE__, nf90_def_var(ncfile_ids(I), TRIM(  savepar(J)%name  ), savepar(J)%type, ([dbins_id, dtime_id])  , savepar(J)%i) )
        if (savepar(J)%d == 4) call handler(__LINE__, nf90_def_var(ncfile_ids(I), TRIM(  savepar(J)%name  ), savepar(J)%type, dcond_id  , savepar(J)%i) )
        if (savepar(J)%d == 5) call handler(__LINE__, nf90_def_var(ncfile_ids(I), TRIM(  savepar(J)%name  ), savepar(J)%type, ([dcond_id,dtime_id])  , savepar(J)%i) )
        if (savepar(J)%d == 7) call handler(__LINE__, nf90_def_var(ncfile_ids(I), TRIM(  savepar(J)%name  ), savepar(J)%type, ([dcond_id,dbins_id,dtime_id])  , savepar(J)%i) )
        if (savepar(J)%d == -12) call handler(__LINE__, nf90_def_var(ncfile_ids(I), TRIM(  savepar(J)%name  ), savepar(J)%type, ([dstring_id,dcond_id])  , savepar(J)%i) )
        call handler(__LINE__, nf90_def_var_deflate(ncfile_ids(I), savepar(J)%i, shuff, compress, compression) )
        call handler(__LINE__, nf90_put_att(ncfile_ids(I), savepar(J)%i, 'unit' , savepar(J)%u))
      end do

      ! if (DONT_SAVE_CONDENSIBLES .eqv. .false.) THEN
        do j = 1,vapours%n_condorg-1
          k = IndexFromName( vapours%vapour_names(j), SPC_NAMES )
          if (k>0) THEN
            call handler(__LINE__, nf90_def_var(ncfile_ids(I), TRIM(  vapours%vapour_names(j)  ), NF90_DOUBLE, dtime_id, par_ind(j)) )
            call handler(__LINE__, nf90_def_var_deflate(ncfile_ids(I), par_ind(j), shuff, compress, compression) )
            call handler(__LINE__, nf90_put_att(ncfile_ids(I), par_ind(j), 'unit' , '1/cm^3'))
          end if
        end do
        ! Generic vapour - these values are only for consistency, they will be zero for gas concetrations
        call handler(__LINE__, nf90_def_var(ncfile_ids(I), TRIM(  vapours%vapour_names(vapours%ind_GENERIC)  ), NF90_DOUBLE, dtime_id, par_ind(j)) )
        call handler(__LINE__, nf90_def_var_deflate(ncfile_ids(I), par_ind(j), shuff, compress, compression) )
        call handler(__LINE__, nf90_put_att(ncfile_ids(I), par_ind(j), 'unit' , '1/cm^3'))
        ! Sulfuric acid
        call handler(__LINE__, nf90_def_var(ncfile_ids(I), TRIM(  vapours%vapour_names(vapours%ind_H2SO4)  ), NF90_DOUBLE, dtime_id, par_ind(j+1)) )
        call handler(__LINE__, nf90_def_var_deflate(ncfile_ids(I), par_ind(j+1), shuff, compress, compression) )
        call handler(__LINE__, nf90_put_att(ncfile_ids(I), par_ind(j+1), 'unit' , '1/cm^3'))
      ! end if
  end if


  ! Ending definition mode
  DO I=1,N_FILES
    call handler(__LINE__, nf90_enddef(ncfile_ids(I)))
  END DO

  I=1

  ! if (Aerosol_flag) call handler(__LINE__, nf90_put_var(ncfile_ids(3), savepar(1)%i, vapours%vapour_names, count=([vapours%n_condtot])))
  if (Aerosol_flag) THEN
    do j = 1,size(vapours%vapour_names)
      k = IndexFromName( vapours%vapour_names(j), SPC_NAMES )
      if (k>0) THEN
        call handler(__LINE__, nf90_put_var(ncfile_ids(3), savepar(1)%i, COND_NAMES))
        i=i+1
      END IF
    end do
  END IF

  print FMT_LEND,

  RETURN

END SUBROUTINE OPEN_FILES


! ====================================================================================================================
! Here the input is written to netcdf-files. Again, particles still rudimentary
! --------------------------------------------------------------------------------------------------------------------
SUBROUTINE SAVE_GASES(TSTEP_CONC,MODS,CH_GAS,conc_vapours,J_ACDC_NH3_M3, J_ACDC_DMA_M3, VAPOURS, save_measured, GR)
  IMPLICIT NONE
  type(input_mod), INTENT(in)     :: MODS(:)
  real(dp), INTENT(in)            :: TSTEP_CONC(:)
  real(dp), INTENT(in)            :: CH_GAS(:)
  real(dp), INTENT(in)            :: conc_vapours(:)
  real(dp), INTENT(in)            :: J_ACDC_NH3_M3
  real(dp), INTENT(in)            :: J_ACDC_DMA_M3
  TYPE(vapour_ambient),INTENT(IN) :: vapours
  real(dp), INTENT(in)            :: save_measured(:)
  real(dp), INTENT(in)            :: GR(:)
  INTEGER                         :: i,j

  DO I = 1,N_FILES
    call handler(__LINE__, nf90_put_var(ncfile_ids(I), timearr_id, GTIME%sec, (/GTIME%ind_netcdf/) ))
    call handler(__LINE__, nf90_put_var(ncfile_ids(I), hrsarr_id, GTIME%hrs, (/GTIME%ind_netcdf/) ))
  END DO

  I=1

  do j = 1,size(MODS)
    IF ((TRIM(MODS(j)%name) /= '#')) THEN
      IF (j<LENV .or. INDRELAY_CH(j)>0) THEN
        call handler(__LINE__, nf90_put_var(ncfile_ids(I), mods_ind(j), TSTEP_CONC(j), (/GTIME%ind_netcdf/)) )
      END IF
      IF ((ABS(MODS(j)%shift-0d0) > 1d-100) .or. (ABS(MODS(j)%multi-1d0) > 1d-100)) THEN
        call handler(__LINE__, nf90_put_var(ncfile_ids(I), multipl_ind(j), MODS(j)%multi, (/GTIME%ind_netcdf/) ) )
        call handler(__LINE__, nf90_put_var(ncfile_ids(I), shifter_ind(j), MODS(j)%shift,  (/GTIME%ind_netcdf/) ) )
      END IF
    END IF
  end do

  call handler(__LINE__, nf90_put_var(ncfile_ids(I), gJ_out_NH3_id, 1d-6*(J_ACDC_NH3_M3), (/GTIME%ind_netcdf/)) )
  call handler(__LINE__, nf90_put_var(ncfile_ids(I), gJ_out_DMA_id, 1d-6*(J_ACDC_DMA_M3), (/GTIME%ind_netcdf/)) )
  call handler(__LINE__, nf90_put_var(ncfile_ids(I), gJ_out_SUM_id, 1d-6*(J_ACDC_DMA_M3+J_ACDC_NH3_M3), (/GTIME%ind_netcdf/)) )
  call handler(__LINE__, nf90_put_var(ncfile_ids(I), gJ_out_TOT_id, 1d-6*(J_TOTAL_M3), (/GTIME%ind_netcdf/)) )
  IF (RESOLVE_BASE) THEN
    call handler(__LINE__, nf90_put_var(ncfile_ids(I), gRES_BASE,RESOLVED_BASE,   (/GTIME%ind_netcdf/)) )
    call handler(__LINE__, nf90_put_var(ncfile_ids(I), gRES_J,   RESOLVED_J,      (/GTIME%ind_netcdf/)) )
    call handler(__LINE__, nf90_put_var(ncfile_ids(I), gRES_F,   RESOLVED_J_FACTR,(/GTIME%ind_netcdf/)) )
  END IF

  I=2 ! Chemical file
  do j = 1,size(CH_GAS)
    call handler(__LINE__, nf90_put_var(ncfile_ids(I), chem_ind(j), CH_GAS(j), (/GTIME%ind_netcdf/)) )
  end do



  I=3 ! Particle file

  if (Aerosol_flag) THEN
    ! if (DONT_SAVE_CONDENSIBLES .eqv. .false.) THEN
        do j = 1,vapours%n_condtot
            call handler(__LINE__, nf90_put_var(ncfile_ids(I), par_ind(j), conc_vapours(j), (/GTIME%ind_netcdf/)) )
        end do
        ! call handler(__LINE__, nf90_put_var(ncfile_ids(I), par_ind(vapours%ind_GENERIC), 0d0, (/GTIME%ind_netcdf/)) )
        ! call handler(__LINE__, nf90_put_var(ncfile_ids(I), par_ind(vapours%ind_H2SO4), conc_vapours(vapours%ind_H2SO4), (/GTIME%ind_netcdf/)) )
    ! end if

    do j = 1,size(savepar)
      if (savepar(j)%name == 'NUMBER_CONCENTRATION'  ) call handler(__LINE__, nf90_put_var(ncfile_ids(I), savepar(J)%i, 1d-6*(get_conc()), start=(/1,GTIME%ind_netcdf/), count=(/n_bins_par/)))
      if (savepar(j)%name == 'INPUT_CONCENTRATION'   ) call handler(__LINE__, nf90_put_var(ncfile_ids(I), savepar(J)%i, save_measured, start=(/1,GTIME%ind_netcdf/), count=(/n_bins_par/)))
      if (savepar(j)%name == 'DIAMETER'              ) call handler(__LINE__, nf90_put_var(ncfile_ids(I), savepar(J)%i, get_dp(), start=(/1,GTIME%ind_netcdf/), count=(/n_bins_par/)))
      ! if (savepar(j)%name == 'DRY_DIAMETER'          ) call handler(__LINE__, nf90_put_var(ncfile_ids(I), savepar(J)%i, current_PSD%dp_dry_fs, start=(/1,GTIME%ind_netcdf/), count=(/n_bins_par/)))
      ! if (savepar(j)%name == 'ORIGINAL_DRY_DIAMETER' ) call handler(__LINE__, nf90_put_var(ncfile_ids(I), savepar(J)%i, current_PSD%original_dry_radius, start=(/1,GTIME%ind_netcdf/), count=(/n_bins_par/)))
      if (savepar(j)%name == 'GROWTH_RATE'           ) call handler(__LINE__, nf90_put_var(ncfile_ids(I), savepar(J)%i, GR, start=(/1,GTIME%ind_netcdf/), count=(/n_bins_par/)))
      if (savepar(j)%name == 'CORE_VOLUME'           ) call handler(__LINE__, nf90_put_var(ncfile_ids(I), savepar(J)%i, get_volume(), start=(/1,GTIME%ind_netcdf/), count=(/n_bins_par/)))
      if (savepar(j)%name == 'MASS'                  ) call handler(__LINE__, nf90_put_var(ncfile_ids(I), savepar(J)%i, get_mass(), start=(/1,GTIME%ind_netcdf/), count=(/n_bins_par/)))
      if (savepar(j)%name == 'PARTICLE_COMPOSITION'  ) call handler(__LINE__, nf90_put_var(ncfile_ids(I), savepar(J)%i, TRANSPOSE(get_composition()), (/1,1,GTIME%ind_netcdf/), (/vapours%n_condtot,n_bins_par,1/)))
      ! if (savepar(j)%name == 'SURFACE_TENSION'       ) call handler(__LINE__, nf90_put_var(ncfile_ids(I), savepar(J)%i, current_PSD%conc_fs, start=(/1,GTIME%ind_netcdf/), count=(/n_bins_par/)))
      ! if (savepar(j)%name == 'VAPOR_CONCENTRATION'   ) call handler(__LINE__, nf90_put_var(ncfile_ids(I), savepar(J)%i, current_PSD%conc_fs, start=(/1,GTIME%ind_netcdf/), count=(/n_bins_par/)))

      ! if (savepar(J)%d == 0) call handler(__LINE__, nf90_put_var(ncfile_ids(I), TRIM(  savepar(J)%name  ), savepar(J)%type, dconstant_id  , savepar(J)%i) )
      ! if (savepar(J)%d == 3) call handler(__LINE__, nf90_put_var(ncfile_ids(I), savepar(J)%i, current_PSD%conc_fs, start=(/1,GTIME%ind_netcdf/), count=(/n_bins_par/)))
      ! if (savepar(J)%d == 4) call handler(__LINE__, nf90_put_var(ncfile_ids(I), TRIM(  savepar(J)%name  ), savepar(J)%type, dcond_id  , savepar(J)%i) )
      ! if (savepar(J)%d == 5) call handler(__LINE__, nf90_put_var(ncfile_ids(I), TRIM(  savepar(J)%name  ), savepar(J)%type, ([dcond_id,dtime_id])  , savepar(J)%i) )
      ! if (savepar(J)%d == 7) call handler(__LINE__, nf90_put_var(ncfile_ids(I), TRIM(  savepar(J)%name  ), savepar(J)%type, ([dcond_id,dbins_id,dtime_id])  , savepar(J)%i) )
      ! if (savepar(J)%d == -12) call handler(__LINE__, nf90_put_var(ncfile_ids(I), TRIM(  savepar(J)%name  ), savepar(J)%type, ([dstring_id,dcond_id])  , savepar(J)%i) )
    end do
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
    call handler(__LINE__, nf90_close(ncfile_ids(I)))
  END DO
  DO I=1,N_FILES
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
    ! TRANSPOSE(get_composition()), (/1,1,GTIME%ind_netcdf/), (/vapours%n_condtot,n_bins_par,1/))


    call handler(__LINE__, nf90_get_var(file_id, varid, c_pp_tp, start=(/1,1,datarow/), count=(/size(c_pp, 2), n_bins_par/)))
    c_pp = TRANSPOSE(c_pp_tp)
    call handler(__LINE__, nf90_close(file_id) ) ! close netcdf dataset
END IF

END SUBROUTINE INITIALIZE_WITH_LAST

END MODULE OUTPUT
