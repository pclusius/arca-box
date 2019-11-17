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
CHARACTER(200) :: ncfile_names(N_FILES) = (['general ', 'chemical', 'particle'])

INTEGER, allocatable        :: shifter_ind(:)
INTEGER, allocatable        :: multipl_ind(:)
INTEGER, allocatable        :: mods_ind(:)
INTEGER, allocatable        :: chem_ind(:)
INTEGER, allocatable        :: par_ind(:)

INTEGER :: dtime_id
INTEGER :: dbins_id
INTEGER :: dcons_id
INTEGER :: dconstant_id
INTEGER :: timearr_id
INTEGER :: hrsarr_id
INTEGER :: gJ_out_NH3_id
INTEGER :: gJ_out_DMA_id

public :: OPEN_FILES, SAVE_GASES,CLOSE_FILES

type(parsave) :: parbuf(20)
type(parsave), ALLOCATABLE :: savepar(:)

CONTAINS


  SUBROUTINE OPEN_FILES(filename, Description, MODS, CH_GAS, vapours)
    IMPLICIT NONE

    CHARACTER(LEN=*), INTENT(IN)    :: filename
    CHARACTER(*), INTENT(IN)        :: Description
    TYPE(input_mod),INTENT(IN)      :: MODS(:)
    TYPE(vapour_ambient),INTENT(IN) :: vapours
    REAL(dp), INTENT(IN)            :: CH_GAS(:)
    CHARACTER(255)                  :: PROGRAM_NAME
    INTEGER                         :: i,j,k,lenD, n_condensables=0

    ! ---------------------------
    ! ------dimension bits-------
    ! time | bins | condensables
    !  1      2          4
    ! eg. number concentration:
    ! time | bins | condensables
    !  1      1          0         => 011 => parbuf(i)%d = 3
    i=1
    parbuf(i)%name = 'NUMBER_CONCENTRATION'  ; parbuf(i)%u = '[1/m^3]'    ; parbuf(i)%d = 3 ; parbuf(i)%type = NF90_DOUBLE ; i=i+1
    parbuf(i)%name = 'DIAMETER'              ; parbuf(i)%u = '[m]'        ; parbuf(i)%d = 3 ; parbuf(i)%type = NF90_DOUBLE ; i=i+1
    parbuf(i)%name = 'DRY_DIAMETER'          ; parbuf(i)%u = '[m]'        ; parbuf(i)%d = 3 ; parbuf(i)%type = NF90_DOUBLE ; i=i+1
    parbuf(i)%name = 'ORIGINAL_DRY_DIAMETER' ; parbuf(i)%u = '[m]'        ; parbuf(i)%d = 3 ; parbuf(i)%type = NF90_DOUBLE ; i=i+1
    parbuf(i)%name = 'GROWTH_RATE'           ; parbuf(i)%u = '[m/s]'      ; parbuf(i)%d = 3 ; parbuf(i)%type = NF90_DOUBLE ; i=i+1
    parbuf(i)%name = 'CORE_VOLUME'           ; parbuf(i)%u = '[m^3]'      ; parbuf(i)%d = 3 ; parbuf(i)%type = NF90_DOUBLE ; i=i+1
    parbuf(i)%name = 'MASS'                  ; parbuf(i)%u = '[kg]'       ; parbuf(i)%d = 3 ; parbuf(i)%type = NF90_DOUBLE ; i=i+1
    parbuf(i)%name = 'VOLUME_CONCENTRATION'  ; parbuf(i)%u = '[um^3/m^3]' ; parbuf(i)%d = 7 ; parbuf(i)%type = NF90_DOUBLE ; i=i+1
    parbuf(i)%name = 'VAPOR_CONCENTRATION'   ; parbuf(i)%u = '[#/m^3]'    ; parbuf(i)%d = 5 ; parbuf(i)%type = NF90_DOUBLE ; i=i+1
    parbuf(i)%name = 'SURFACE_TENSION'       ; parbuf(i)%u = '[N/m]'      ; parbuf(i)%d = 4 ; parbuf(i)%type = NF90_DOUBLE ; i=i+1

    allocate(savepar(i-1))
    savepar = parbuf(1:i-1)

    do j = 1,size(vapours%vapour_names)
      k = IndexFromName( vapours%vapour_names(j), SPC_NAMES )
      if (k>0) n_condensables=n_condensables+1
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
    print FMT_SUB, 'Create chemfiles: '//TRIM(filename)//'_*.nc'

    DO I=1, N_FILES

      ncfile_names(I) = trim(filename)//'_'//TRIM(ncfile_names(I))//'.nc'

      !Clearing file; Opening file. Overwrites
      open(720+I, FILE=ncfile_names(I), ERR = 100)
      close(720+I)
      ! Added compression for particle.nc, so we need to use netCDF4-file. Here used in classic mode
      ! call handler( nf90_create(ncfile_names(I), IOR(NF90_NETCDF4, NF90_CLASSIC_MODEL), ncfile_ids(I)) )
      call handler( nf90_create(ncfile_names(I), NF90_NETCDF4, ncfile_ids(I)) )

      ! Defining dimensions: time(unlimited), size sections, vapor_species
      call handler(nf90_def_dim(ncfile_ids(I), "time",NINT(MODELTIME%SIM_TIME_S/MODELTIME%FSAVE_INTERVAL+1), dtime_id) )
      call handler(nf90_def_dim(ncfile_ids(I), "bins",n_bins_particle, dbins_id) )
      call handler(nf90_def_dim(ncfile_ids(I), "condensables",n_condensables, dcons_id) )
      call handler(nf90_def_dim(ncfile_ids(I), "Constant",1, dconstant_id) )

      !call handler(nf90_def_dim(ncfile_ids(I), "StringL",16, gstrlen_id) )

      !Create attributes for general stuff
      CALL get_command_argument(0, PROGRAM_NAME)
      call handler(nf90_put_att(ncfile_ids(I), NF90_GLOBAL, 'Information', '(c) Atmospheric modelling group 2019 and (c) Simugroup 2019 (ACDC)'))
      call handler(nf90_put_att(ncfile_ids(I), NF90_GLOBAL, 'Contact', 'michael.boy@helsinki.fi (Superbox), tinja.olenius@alumni.helsinki.fi (ACDC)'))
      call handler(nf90_put_att(ncfile_ids(I), NF90_GLOBAL, 'Software', 'Superbox 0.0.1'))
      call handler(nf90_put_att(ncfile_ids(I), NF90_GLOBAL, 'Package_Name:', TRIM(PROGRAM_NAME(3:))))
      call handler(nf90_put_att(ncfile_ids(I), NF90_GLOBAL, 'Notes', TRIM(Description)))
      call handler(nf90_put_att(ncfile_ids(I), NF90_GLOBAL, 'experiment', 'Experiment set here'))

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

  call handler(nf90_def_var_deflate(ncfile_ids(I), gJ_out_NH3_id,         shuff, compress, compression) )
  call handler(nf90_def_var_deflate(ncfile_ids(I), gJ_out_DMA_id,         shuff, compress, compression) )

  call handler(nf90_put_att(ncfile_ids(I), gJ_out_NH3_id, 'unit' , '[1/s/m^3]'))
  call handler(nf90_put_att(ncfile_ids(I), gJ_out_DMA_id, 'unit' , '[1/s/m^3]'))


  I=2 ! Chemical file

  do j = 1,size(CH_GAS)

      call handler(nf90_def_var(ncfile_ids(I), TRIM(SPC_NAMES(j)), NF90_DOUBLE, dtime_id, chem_ind(j)) )
      call handler(nf90_def_var_deflate(ncfile_ids(I), chem_ind(j), shuff, compress, compression) )
      call handler(nf90_put_att(ncfile_ids(I), chem_ind(j), 'unit' , '1/cm^3'))

  end do


  I=3 ! Paricle file

  do j = 1,size(savepar)
    if (savepar(J)%d == 0) call handler(nf90_def_var(ncfile_ids(I), TRIM(  savepar(J)%name  ), savepar(J)%type, dconstant_id  , savepar(J)%i) )
    if (savepar(J)%d == 3) call handler(nf90_def_var(ncfile_ids(I), TRIM(  savepar(J)%name  ), savepar(J)%type, ([dbins_id, dtime_id])  , savepar(J)%i) )
    if (savepar(J)%d == 4) call handler(nf90_def_var(ncfile_ids(I), TRIM(  savepar(J)%name  ), savepar(J)%type, dcons_id  , savepar(J)%i) )
    if (savepar(J)%d == 5) call handler(nf90_def_var(ncfile_ids(I), TRIM(  savepar(J)%name  ), savepar(J)%type, ([dcons_id,dtime_id])  , savepar(J)%i) )
    if (savepar(J)%d == 7) call handler(nf90_def_var(ncfile_ids(I), TRIM(  savepar(J)%name  ), savepar(J)%type, ([dcons_id,dbins_id,dtime_id])  , savepar(J)%i) )

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



  !Ending definition
  DO I=1,N_FILES
    call handler( nf90_enddef(ncfile_ids(I)))
  END DO
  print FMT_LEND,


  open(9129, file=filename//'_INITFILE.txt', action='WRITE')
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





  RETURN
  100 continue
  print *, "Error in opening NetCDF file: "//filename//""
  stop


END SUBROUTINE OPEN_FILES














SUBROUTINE SAVE_GASES(TSTEP_CONC,MODS,CH_GAS,J_ACDC_NH3, J_ACDC_DMA, vapours)
  IMPLICIT NONE
  type(input_mod), INTENT(in)     :: MODS(:)
  real(dp), INTENT(in)            :: TSTEP_CONC(:)
  real(dp), INTENT(in)            :: CH_GAS(:)
  TYPE(vapour_ambient),INTENT(IN) :: vapours
  real(dp), INTENT(in)            :: J_ACDC_NH3
  real(dp), INTENT(in)            :: J_ACDC_DMA
  INTEGER                         :: i,j, k

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


  I=2 ! Chemical file

  do j = 1,size(CH_GAS)
    call handler( nf90_put_var(ncfile_ids(I), chem_ind(j), CH_GAS(j), (/MODELTIME%ind_netcdf/)) )
  end do

  I=3 ! Chemical file

  do j = 1,size(vapours%vapour_names)
    k = IndexFromName( vapours%vapour_names(j),   SPC_NAMES )
    if (k>0) call handler(nf90_put_var(ncfile_ids(I), par_ind(j), CH_GAS(k), (/MODELTIME%ind_netcdf/)) )

  end do


END SUBROUTINE SAVE_GASES



subroutine CLOSE_FILES()
  IMPLICIT NONE
  INTEGER :: I
  DO I=1,N_FILES
    call handler( nf90_close(ncfile_ids(I)))
  END DO

  Write(*,FMT_MSG) 'Outputfiles closed.'
  Write(*,FMT_LEND)

end subroutine CLOSE_FILES

subroutine HANDLER(status)
  IMPLICIT NONE
  INTEGER, intent(in) :: status
  if (status /= nf90_noerr) then
    write(*,*) 'Error writing netcdf, code: ',NF90_STRERROR(status)
    stop
  end if
end subroutine HANDLER

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
