MODULE OUTPUT
use netcdf
use CONSTANTS
use AUXILLARIES
IMPLICIT NONE
private
! Storing of indices of gases.nc
character(16), allocatable  :: gasnames(:)
INTEGER, allocatable        :: shifter_ind(:)
INTEGER, allocatable        :: multipl_ind(:)

integer :: gtime_id = 0
integer :: gas_ncfile_id = 0
integer :: gcompounds_id = 0
integer :: gas_concentrations_id
integer :: gas_names_id
integer :: gtemperature_id
integer :: gpressure_id
integer :: gJ_out_NH3_id
integer :: gJ_out_DMA_id
integer :: gc_sink_id
integer :: gnotes_id
integer :: gconstant_id = 0
integer :: gtimearr_id
integer :: gstrlen_id

! Storing of indices of particle.nc
integer :: par_ncfile_id = 0
integer :: variable_id = 0
integer :: time_id = 0
integer :: section_id = 0
integer :: composition_id = 0
integer :: nconc_id = 0
integer :: radius_id = 0
integer :: mass_id = 0
integer :: Vol_conc_id = 0
integer :: timevar_id = 0
integer :: dry_radius_id = 0
integer :: orig_dry_radius_id = 0
integer :: growth_rate_id = 0
integer :: condensable_id = 0
integer :: vap_conc_id = 0
integer :: constant_id = 0
integer :: temperature_id = 0
integer :: pressure_id = 0
integer :: saturation_ratio_id = 0
integer :: nucl_coef_id = 0
integer :: boundaryheight_id = 0
integer :: vapordensity_id = 0
integer :: molarmass_id = 0
integer :: saturation_id = 0
integer :: surfacetension_id = 0
integer :: diffusionvolume_id = 0
integer :: massaccomondation_id = 0
integer :: core_volume_id = 0
integer :: order_id = 0

integer, dimension(2) :: particles_nconc_dimension_ids, particles_radius_dimension_ids, &
particles_mass_dimension_ids, ambient_vap_conc_ids
integer, dimension(3) ::  part_volconc_dim_id

public :: output_particles,output_ambient, initialize_output,CLOSE_FILES, &
 output_constants, output_time, OPEN_GASFILE, SAVE_GASES, OPEN_NETCDF

CONTAINS

SUBROUTINE OPEN_NETCDF(filename, aerosol_options)
  type(type_options), intent(in)  :: aerosol_options
  character(len=*), intent(in)    :: filename
  !Initialize netcdf-output
  write(*,*) TRIM(filename)
  call initialize_output(TRIM(filename), 'explanation', 'experiment_set', aerosol_options)
end SUBROUTINE OPEN_NETCDF


SUBROUTINE OPEN_GASFILE(filename, modifier_names, cons_multipliers, cons_shifters)
  implicit none
  CHARACTER(255)     :: PROGRAM_NAME
  real(dp)      :: cons_multipliers(:)
  real(dp)      :: cons_shifters(:)
  CHARACTER(*)  :: modifier_names(:)
  integer :: i
  integer                         :: n_compounds = 5 ! needed for develempoment, later will be calculated from KPP input and stored in e.g. Constants.f90
  character(len=*), intent(in)    :: filename
  ! settings for netCDF4-file compression. shuff=1 might improve compression but is slower
  integer:: shuff=1, compress=1, compression=9
  integer::gambient_gases_ids(2)

  ALLOCATE(multipl_ind(size(modifier_names)))
  ALLOCATE(shifter_ind(size(modifier_names)))
  print *, 'NetCDF version: ', trim(nf90_inq_libvers())
  write(*,*) 'Create chemfile', TRIM(filename)


  !Clearing file; Opening file. Overwrites
  open(999, FILE=filename, ERR = 100)
  close(999)
  ! Added compression for particle.nc, so we need to use netCDF4-file. Here used in classic mode
  call handler( nf90_create(trim(filename), IOR(NF90_NETCDF4, NF90_CLASSIC_MODEL), gas_ncfile_id) )

  ! Defining dimensions: time(unlimited), size sections, vapor_species
  call handler(nf90_def_dim(gas_ncfile_id, "time",1442, gtime_id) )
  call handler(nf90_def_dim(gas_ncfile_id, "Compound",n_compounds, gcompounds_id) )
  call handler(nf90_def_dim(gas_ncfile_id, "Constant",1, gconstant_id) )
  call handler(nf90_def_dim(gas_ncfile_id, "StringL",16, gstrlen_id) )

  !Identifying different shapes for arrays
  !Ambient:
  gambient_gases_ids = (/gcompounds_id, gtime_id/)

  !Ambient:
  call handler(nf90_def_var(gas_ncfile_id, "time_in_sec", NF90_DOUBLE, gtime_id, gtimearr_id))
  call handler(nf90_def_var(gas_ncfile_id, "gas_names", NF90_CHAR, ([gstrlen_id, gcompounds_id]), gas_names_id) )
  call handler(nf90_def_var(gas_ncfile_id, "gas_concentrations", NF90_DOUBLE, gambient_gases_ids, gas_concentrations_id ) )
  call handler(nf90_def_var(gas_ncfile_id, "temperature", NF90_DOUBLE, gtime_id, gtemperature_id))
  call handler(nf90_def_var(gas_ncfile_id, "pressure", NF90_DOUBLE, gtime_id, gpressure_id))
  call handler(nf90_def_var(gas_ncfile_id, "J_out_NH3", NF90_DOUBLE, gtime_id, gJ_out_NH3_id))
  call handler(nf90_def_var(gas_ncfile_id, "J_out_DMA", NF90_DOUBLE, gtime_id, gJ_out_DMA_id))
  call handler(nf90_def_var(gas_ncfile_id, "c_sink", NF90_DOUBLE, gtime_id, gc_sink_id))
  ! COMPRESSION
  call handler(nf90_def_var_deflate(gas_ncfile_id, gtime_id,              shuff, compress, compression) )
  call handler(nf90_def_var_deflate(gas_ncfile_id, gas_names_id,          shuff, compress, compression) )
  call handler(nf90_def_var_deflate(gas_ncfile_id, gas_concentrations_id, shuff, compress, compression) )
  call handler(nf90_def_var_deflate(gas_ncfile_id, gtemperature_id,       shuff, compress, compression) )
  call handler(nf90_def_var_deflate(gas_ncfile_id, gpressure_id,          shuff, compress, compression) )
  call handler(nf90_def_var_deflate(gas_ncfile_id, gJ_out_NH3_id,         shuff, compress, compression) )
  call handler(nf90_def_var_deflate(gas_ncfile_id, gJ_out_DMA_id,         shuff, compress, compression) )
  call handler(nf90_def_var_deflate(gas_ncfile_id, gc_sink_id,            shuff, compress, compression) )
  ! END COMPRESSION

  !Constants:

  !Create attributes for general stuff
  CALL get_command_argument(0, PROGRAM_NAME)
  call handler(nf90_put_att(gas_ncfile_id, NF90_GLOBAL, 'Information', '(c) Atmospheric modelling group 2019 and (c) Simugroup 2019 (ACDC)'))
  call handler(nf90_put_att(gas_ncfile_id, NF90_GLOBAL, 'Contact', 'michael.boy@helsinki.fi (Superbox), tinja.olenius@alumni.helsinki.fi (ACDC)'))
  call handler(nf90_put_att(gas_ncfile_id, NF90_GLOBAL, 'Software', 'Superbox 0.1'))
  call handler(nf90_put_att(gas_ncfile_id, NF90_GLOBAL, 'Package_Name:', TRIM(PROGRAM_NAME(3:))))
  call handler(nf90_put_att(gas_ncfile_id, NF90_GLOBAL, 'Notes', 'e.g. Sulfuric acid concentration multiplied by 0.1'))
  call handler(nf90_put_att(gas_ncfile_id, NF90_GLOBAL, 'experiment', 'Experiment set here'))
  ! call output_options(options)

  do i = 1,size(modifier_names)
    call handler(nf90_def_var(gas_ncfile_id, TRIM(modifier_names(i))//'_Multipl', NF90_DOUBLE, gtime_id, multipl_ind(i)))
    call handler(nf90_def_var(gas_ncfile_id, TRIM(modifier_names(i))//'_Shifter', NF90_DOUBLE, gtime_id, shifter_ind(i)))
  end do

  !defining 'units' as attributes.

  !Particles
  call handler(nf90_put_att(gas_ncfile_id, gas_concentrations_id, 'units' , '1/m^3'))
  call handler(nf90_put_att(gas_ncfile_id, gas_names_id, 'units', '[]'))
  call handler(nf90_put_att(gas_ncfile_id, gtemperature_id, 'units', 'K'))
  call handler(nf90_put_att(gas_ncfile_id, gpressure_id, 'units', 'Pa'))
  call handler(nf90_put_att(gas_ncfile_id, gJ_out_NH3_id, 'units', '1/s/m3'))
  call handler(nf90_put_att(gas_ncfile_id, gJ_out_DMA_id, 'units', '1/s/m3'))
  call handler(nf90_put_att(gas_ncfile_id, gc_sink_id, 'units', '1/s'))


  !Ending definition
  call handler( nf90_enddef(gas_ncfile_id))

ALLOCATE(gasnames(5))
gasnames(1) = "Sulf_acid"
gasnames(2) = "Ammonia"
gasnames(3) = "DMA"
gasnames(4) = "undefined"
gasnames(5) = "undefined"

call handler(nf90_put_var(gas_ncfile_id, gas_names_id, (gasnames), start=([1]), count=([16, 5])))

! call handler(nf90_put_var(gas_ncfile_id, gas_names_id, (['time']), start=([1]), count=([4])))
  return
  ! error message if file open failed
  100 continue
  print *, "Error in  opening gasfile: "//filename//""
  stop

end SUBROUTINE OPEN_GASFILE

subroutine output_constants(ambient)
  IMPLICIT NONE

  type(type_ambient), intent(in) :: ambient

  call handler( nf90_put_var(par_ncfile_id, nucl_coef_id, ambient%nuc_coeff, (/1/)))
  call handler( nf90_put_var(par_ncfile_id, boundaryheight_id, ambient%boundaryheight, (/1/)))
  call handler( nf90_put_var(par_ncfile_id, vapordensity_id, ambient%density))
  call handler( nf90_put_var(par_ncfile_id, molarmass_id, ambient%molarmass))
  call handler( nf90_put_var(par_ncfile_id, surfacetension_id, ambient%surf_ten))
  call handler( nf90_put_var(par_ncfile_id, diffusionvolume_id, ambient%diff_vol))
  call handler( nf90_put_var(par_ncfile_id, massaccomondation_id, ambient%alpha))
  !    call handler( nf90_put_var(par_ncfile_id, temperature_id, ambient%temp, (/1/)))
  !    call handler( nf90_put_var(par_ncfile_id, pressure_id, ambient%pres, (/1/)))
  !    call handler( nf90_put_var(par_ncfile_id, saturation_ratio_id, ambient%rh, (/1/)))
  !    call handler( nf90_put_var(par_ncfile_id, saturation_id, ambient%c_sat))

end subroutine output_constants

subroutine output_options(options)
  IMPLICIT NONE

  type(type_options), intent(in) :: options

  call handler(nf90_put_att(par_ncfile_id, NF90_GLOBAL, 'solver', options%solver))
  call handler(nf90_put_att(par_ncfile_id, NF90_GLOBAL, 'approach', options%dist_approach))
  call handler(nf90_put_att(par_ncfile_id, NF90_GLOBAL, 'nucleation', options%nuc_approach))

end subroutine output_options

subroutine output_time(time)
  IMPLICIT NONE
  type(timetype), intent(in) :: time
  call handler(nf90_put_var(par_ncfile_id, timevar_id, time%sec, (/time%ind_netcdf/) ))
end subroutine output_time

subroutine output_particles(particles, timestep)
  IMPLICIT NONE
  type(type_particles), intent(in) :: particles
  integer, intent(in) :: timestep

  integer :: bin, compound
  bin=0
  compound=0

  call handler( nf90_put_var(par_ncfile_id, nconc_id, particles%n_conc, start=(/1, timestep/), count=(/uhma_sections, 1/)))
  call handler( nf90_put_var(par_ncfile_id, radius_id, particles%radius, start=(/1, timestep/), count=(/uhma_sections, 1/)))
  call handler( nf90_put_var(par_ncfile_id, dry_radius_id, particles%rdry, start=(/1, timestep/),count=(/uhma_sections, 1/) ))
  call handler( nf90_put_var(par_ncfile_id, orig_dry_radius_id, particles%rdry_orig, start=(/1, timestep/), count=(/uhma_sections, 1/)))
  call handler( nf90_put_var(par_ncfile_id, growth_rate_id, particles%gr, start=(/1, timestep/), count=(/uhma_sections,1/)))
  call handler( nf90_put_var(par_ncfile_id, mass_id, particles%mass, start=(/1, timestep/), count=(/uhma_sections,1/) ))
  call handler( nf90_put_var(par_ncfile_id, core_volume_id, particles%core, start=(/1, timestep/), count=(/uhma_sections,1/)))
  call handler( nf90_put_var(par_ncfile_id, order_id, particles%order, start=(/1, timestep/), count=(/uhma_sections, 1/)))
  call handler(nf90_put_var(par_ncfile_id, vol_conc_id, particles%vol_conc, (/1, 1, timestep/), (/uhma_sections, uhma_compo, 1/) ))

  return

end subroutine output_particles

subroutine output_ambient(ambient, timestep)
  IMPLICIT NONE
  type(type_ambient), intent(in) :: ambient
  integer, intent(in) :: timestep
  integer :: compound, condensable

  compound=0
  condensable= 0

  call handler( nf90_put_var(par_ncfile_id, temperature_id, ambient%temp, (/timestep/)))
  call handler( nf90_put_var(par_ncfile_id, pressure_id, ambient%pres, (/timestep/)))
  call handler( nf90_put_var(par_ncfile_id, saturation_id,ambient%c_sat, (/1, timestep/), (/uhma_cond, 1/)))
  call handler( nf90_put_var(par_ncfile_id, vap_conc_id,ambient%vap_conc, (/1, timestep/), (/uhma_cond,1/)))

end subroutine output_ambient

subroutine initialize_output(filename, explanation, experiment_set, options)
  IMPLICIT NONE
  ! settings for netCDF4-file compression. shuff=1 might improve compression but is slower
  integer:: shuff=0, compress=1, compression=9

  !Extend as needed. Define stuff here.
  character(len=*), intent(in) :: filename, explanation, experiment_set
  type(type_options) :: options

  !Clearing file; Opening file. Overwrites
  open(999, FILE=filename, ERR = 100)
  close(999)

  ! Added compression for particle.nc, so we need to use netCDF4-file. Here used in classic mode
  call handler( nf90_create(trim(filename), IOR(NF90_NETCDF4, NF90_CLASSIC_MODEL), par_ncfile_id) )


  ! Defining dimensions: time(unlimited), size sections, vapor_species
  call handler(nf90_def_dim(par_ncfile_id, "time",NF90_UNLIMITED, time_id) )
  call handler(nf90_def_dim(par_ncfile_id, "Sections",uhma_sections, section_id) )
  call handler(nf90_def_dim(par_ncfile_id, "Composition",uhma_compo, composition_id) )
  call handler(nf90_def_dim(par_ncfile_id, "Condensable",uhma_cond, condensable_id) )
  call handler(nf90_def_dim(par_ncfile_id, "Constant", 1, constant_id))

  ! Identifying different shapes for arrays

  !Particles
  particles_nconc_dimension_ids = (/section_id, time_id/)
  particles_radius_dimension_ids = (/section_id, time_id/)
  particles_mass_dimension_ids = (/section_id, time_id/)
  part_volconc_dim_id = (/section_id, composition_id, time_id/)
  ambient_vap_conc_ids = (/condensable_id, time_id/)

  !Ambient:
  ! none


  ! Defining variables to write

  !Particles:
  call handler(nf90_def_var(par_ncfile_id, "number_concentration", NF90_DOUBLE, particles_nconc_dimension_ids, nconc_id) )
  call handler(nf90_def_var(par_ncfile_id, "radius", NF90_DOUBLE, particles_radius_dimension_ids, radius_id) )
  call handler(nf90_def_var(par_ncfile_id, "dry_radius", NF90_DOUBLE, particles_radius_dimension_ids, dry_radius_id) )
  call handler(nf90_def_var(par_ncfile_id, "original_dry_radius", NF90_DOUBLE, particles_radius_dimension_ids, orig_dry_radius_id) )
  call handler(nf90_def_var(par_ncfile_id, "growth_rate", NF90_DOUBLE, particles_radius_dimension_ids, growth_rate_id) )
  call handler(nf90_def_var(par_ncfile_id, "core_volume", NF90_DOUBLE, particles_radius_dimension_ids, core_volume_id) )
  call handler(nf90_def_var(par_ncfile_id, "mass", NF90_DOUBLE, particles_mass_dimension_ids, mass_id) )
  call handler(nf90_def_var(par_ncfile_id, "time_in_'units'", NF90_DOUBLE, time_id, timevar_id) )
  call handler(nf90_def_var(par_ncfile_id, "volume_concentration", NF90_DOUBLE, part_volconc_dim_id, vol_conc_id) )
  call handler(nf90_def_var(par_ncfile_id, "Bin order", NF90_INT, particles_radius_dimension_ids, order_id))

  !Ambient:
  call handler(nf90_def_var(par_ncfile_id, "vapor_concentration", NF90_DOUBLE, ambient_vap_conc_ids, vap_conc_id) )
  call handler(nf90_def_var(par_ncfile_id, "temperature", NF90_DOUBLE, time_id, temperature_id))
  call handler(nf90_def_var(par_ncfile_id, "pressure", NF90_DOUBLE, time_id, pressure_id))

  !Constants:
  call handler(nf90_def_var(par_ncfile_id, "nucleation_coefficent", NF90_DOUBLE, constant_id, nucl_coef_id))
  call handler(nf90_def_var(par_ncfile_id, "boundary_layer_height", NF90_DOUBLE, constant_id, boundaryheight_id))

  !vapor parameters
  call handler(nf90_def_var(par_ncfile_id, "vapor_density", NF90_DOUBLE, condensable_id, vapordensity_id))
  call handler(nf90_def_var(par_ncfile_id, "molarmass", NF90_DOUBLE, condensable_id, molarmass_id))
  call handler(nf90_def_var(par_ncfile_id, "saturation_concentration", NF90_DOUBLE, ambient_vap_conc_ids, saturation_id))
  ! call handler(nf90_def_var(par_ncfile_id, "saturation_concentration", NF90_DOUBLE, condensable_id, saturation_id))
  call handler(nf90_def_var(par_ncfile_id, "surface_tension", NF90_DOUBLE, condensable_id, surfacetension_id))
  call handler(nf90_def_var(par_ncfile_id, "diffusion_volume", NF90_DOUBLE, condensable_id, diffusionvolume_id))
  call handler(nf90_def_var(par_ncfile_id, "mass_accomondation_coefficient", NF90_DOUBLE, condensable_id, massaccomondation_id))

  !Making info variable(well attribute) for general stuff
  call handler(nf90_put_att(par_ncfile_id, NF90_GLOBAL, 'info', explanation))
  call handler(nf90_put_att(par_ncfile_id, NF90_GLOBAL, 'experiment', experiment_set))
  call output_options(options)

!defining 'units' as attributes.

  !Particles
  call handler(nf90_put_att(par_ncfile_id, nconc_id, 'units' , '1/m^3'))
  call handler(nf90_put_att(par_ncfile_id, radius_id, 'units', 'm'))
  call handler(nf90_put_att(par_ncfile_id, dry_radius_id, 'units', 'm'))
  call handler(nf90_put_att(par_ncfile_id, orig_dry_radius_id, 'units', 'm'))
  call handler(nf90_put_att(par_ncfile_id, growth_rate_id, 'units', 'm/s'))
  call handler(nf90_put_att(par_ncfile_id, core_volume_id, 'units', 'm^3'))
  call handler(nf90_put_att(par_ncfile_id, mass_id, 'units', 'kg'))
  call handler(nf90_put_att(par_ncfile_id, timevar_id, 'units', 's'))
  call handler(nf90_put_att(par_ncfile_id, vol_conc_id, 'units', 'um^3/m^3'))
  call handler(nf90_put_att(par_ncfile_id, order_id, 'units', ''))

  !Ambient
  call handler(nf90_put_att(par_ncfile_id, vap_conc_id, 'units', '1/m^3'))
  call handler(nf90_put_att(par_ncfile_id, temperature_id, 'units', 'K'))
  call handler(nf90_put_att(par_ncfile_id, pressure_id, 'units', 'Pa'))

  !constants:
  ! call handler(nf90_put_att(par_ncfile_id, saturation_ratio_id, 'units', '1/1'))
  call handler(nf90_put_att(par_ncfile_id, nucl_coef_id, 'units', 'Si units'))
  call handler(nf90_put_att(par_ncfile_id, boundaryheight_id, 'units', 'm?'))

  ! gas parameters
  call handler(nf90_put_att(par_ncfile_id, vapordensity_id, 'units', 'kg/m^3'))
  call handler(nf90_put_att(par_ncfile_id, molarmass_id, 'units', 'g/mol'))
  call handler(nf90_put_att(par_ncfile_id, saturation_id, 'units', '1/m^3'))
  call handler(nf90_put_att(par_ncfile_id, surfacetension_id, 'units', 'N/m'))
  call handler(nf90_put_att(par_ncfile_id, diffusionvolume_id, 'units', '????'))
  call handler(nf90_put_att(par_ncfile_id, massaccomondation_id, 'units', '1/1'))

!Ending definition

! COMPRESSION (added by PC)
  call handler(nf90_def_var_deflate(par_ncfile_id, nconc_id             , shuff, compress, compression) )
  call handler(nf90_def_var_deflate(par_ncfile_id, radius_id            , shuff, compress, compression) )
  call handler(nf90_def_var_deflate(par_ncfile_id, dry_radius_id        , shuff, compress, compression) )
  call handler(nf90_def_var_deflate(par_ncfile_id, orig_dry_radius_id   , shuff, compress, compression) )
  call handler(nf90_def_var_deflate(par_ncfile_id, growth_rate_id       , shuff, compress, compression) )
  call handler(nf90_def_var_deflate(par_ncfile_id, core_volume_id       , shuff, compress, compression) )
  call handler(nf90_def_var_deflate(par_ncfile_id, mass_id              , shuff, compress, compression) )
  call handler(nf90_def_var_deflate(par_ncfile_id, timevar_id           , shuff, compress, compression) )
  call handler(nf90_def_var_deflate(par_ncfile_id, vol_conc_id          , shuff, compress, compression) )
  call handler(nf90_def_var_deflate(par_ncfile_id, order_id             , shuff, compress, compression) )
  call handler(nf90_def_var_deflate(par_ncfile_id, vap_conc_id          , shuff, compress, compression) )
  call handler(nf90_def_var_deflate(par_ncfile_id, temperature_id       , shuff, compress, compression) )
  call handler(nf90_def_var_deflate(par_ncfile_id, pressure_id          , shuff, compress, compression) )
  call handler(nf90_def_var_deflate(par_ncfile_id, nucl_coef_id         , shuff, compress, compression) )
  call handler(nf90_def_var_deflate(par_ncfile_id, boundaryheight_id    , shuff, compress, compression) )
  call handler(nf90_def_var_deflate(par_ncfile_id, vapordensity_id      , shuff, compress, compression) )
  call handler(nf90_def_var_deflate(par_ncfile_id, molarmass_id         , shuff, compress, compression) )
  call handler(nf90_def_var_deflate(par_ncfile_id, saturation_id        , shuff, compress, compression) )
  call handler(nf90_def_var_deflate(par_ncfile_id, surfacetension_id    , shuff, compress, compression) )
  call handler(nf90_def_var_deflate(par_ncfile_id, diffusionvolume_id   , shuff, compress, compression) )
  call handler(nf90_def_var_deflate(par_ncfile_id, massaccomondation_id , shuff, compress, compression) )
! END COMPRESSION

  call handler( nf90_enddef(par_ncfile_id))

  return
  ! error message if file open failed
  100 continue
  print *, "Error in output, in opening: "//filename//""
  stop

end subroutine initialize_output

SUBROUTINE SAVE_GASES(time, temperature, C_H2SO4, C_NH3, C_DMA, J_NH3, J_DMA, CS, Gases, cons_multipliers, cons_shifters)
  IMPLICIT NONE
  type(type_ambient), INTENT(in)  :: Gases
  type(timetype), INTENT(in)      :: time
  real(dp), INTENT(in)            :: temperature
  real(dp), INTENT(in)            :: C_H2SO4
  real(dp), INTENT(in)            :: C_NH3
  real(dp), INTENT(in)            :: C_DMA
  real(dp), INTENT(in)            :: J_NH3
  real(dp), INTENT(in)            :: J_DMA
  real(dp), INTENT(in)            :: CS
  real(dp)      :: cons_multipliers(:)
  real(dp)      :: cons_shifters(:)
  integer       :: i

  do i = 1,size(cons_shifters)
    call handler(nf90_put_var(gas_ncfile_id, multipl_ind(i), cons_multipliers(i), (/time%ind_netcdf/) ) )
    call handler(nf90_put_var(gas_ncfile_id, shifter_ind(i), cons_shifters(i),  (/time%ind_netcdf/) ) )
  end do


  call handler( nf90_put_var(gas_ncfile_id, gtimearr_id, time%sec, (/time%ind_netcdf/) ))
  call handler( nf90_put_var(gas_ncfile_id, gtemperature_id, temperature, (/time%ind_netcdf/)) )!, count=(/uhma_sections, 1/)))
  call handler( nf90_put_var(gas_ncfile_id, gas_concentrations_id, ([C_H2SO4, C_NH3, C_DMA]), start=(/1, time%ind_netcdf/), count=(/3/)))
  call handler( nf90_put_var(gas_ncfile_id, gJ_out_NH3_id, J_NH3, (/time%ind_netcdf/)) )
  call handler( nf90_put_var(gas_ncfile_id, gJ_out_DMA_id, J_DMA, (/time%ind_netcdf/)) )
  call handler( nf90_put_var(gas_ncfile_id, gc_sink_id, CS, (/time%ind_netcdf/)) )
  ! call handler( nf90_put_var(gas_ncfile_id, gas_concentrations_id, ([ C_H2SO4 ]), start=(/1, time%ind_netcdf/), count=(/1/)))
  ! call handler( nf90_put_var(gas_ncfile_id, gcompounds_id, C_NH3, (/time%ind_netcdf/), count=(/2, 1/)))
  ! call handler( nf90_put_var(gas_ncfile_id, gcompounds_id, C_DMA, (/time%ind_netcdf/), count=(/3, 1/)))
  ! call handler( nf90_put_var(gas_ncfile_id, gtemperature_id, temperature, (/time%ind_netcdf/)) )!, count=(/uhma_sections, 1/)))
  ! call handler( nf90_put_var(gas_ncfile_id, gtemperature_id, temperature, (/time%ind_netcdf/)) )!, count=(/uhma_sections, 1/)))
  ! call handler( nf90_put_var(gas_ncfile_id, nconc_id, particles%n_conc, start=(/1, timestep/), count=(/uhma_sections, 1/)))
END SUBROUTINE SAVE_GASES



subroutine CLOSE_FILES()
  IMPLICIT NONE
  ! call handler( nf90_close(par_ncfile_id))
  call handler( nf90_close(gas_ncfile_id))
  Write(*,*) 'Outputfiles closed.'
end subroutine CLOSE_FILES

subroutine HANDLER(status)
  IMPLICIT NONE
  integer, intent(in) :: status
  if (status /= nf90_noerr) then
    write(*,*) 'Error writing netcdf, code: ',NF90_STRERROR(status)
    stop
  end if
end subroutine HANDLER

END MODULE OUTPUT
