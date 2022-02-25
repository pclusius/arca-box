#!/bin/bash

ACDC_EXECUTABLE=acdc_2020_04_28.pl

echo
echo "---------------------------------------------------------------------"
echo "Parsing cluster system with ACDC version \"${ACDC_EXECUTABLE}"\"
echo "---------------------------------------------------------------------"
echo

cd $2

DIPOLES=''
NAME=''

if [ "$1" = "" ];then
  echo Using settings from settings.conf
  source settings.conf
else
  echo Using settings from $1
  source $1
fi

loc=`pwd`
i=$((${#loc}-1))
no_system="${loc:$i:1}"
suffix=_0x$no_system

# Create the Perl option string

perl_opt=""
# cluster_file="input_"
# vapor_suffix="_"

if [ -n "$TEMPERATURE" ]; then
    perl_opt+=" --TEMPERATURE $TEMPERATURE"
else
    perl_opt+=" --variable_temp"
fi

if [ -n "$RH" ]; then
    perl_opt+=" --rh $RH"
fi

for vapor in "${VAPORS[@]}"; do
    perl_opt+=" --cs_only 1$vapor,0"
    [ $L_CONST_VAPOR -eq 1 ] && perl_opt+=" --no_eq 1$vapor"

    # cluster_file+="$vapor"
    # vapor_suffix+="$vapor"
done

# cluster_file+="narrow_neutral"

if [ $L_INCL_IONS -eq 1 ]; then
    perl_opt+=" --variable_ion_source"
    DIPOLES="Perl_input/"$DIPOLES
    perl_opt+=" --dip "
    perl_opt+=$DIPOLES
    # cluster_file+="_neg_pos"
    # vapor_suffix+="_ions"
# else
    # vapor_suffix+="_noions"
fi

#cluster_file+=".inp"
ENERGIES="Perl_input/"$ENERGIES
SYSCONF="Perl_input/"$SYSCONF

# Generate the equations
echo Perl call:
echo --fortran --save_outgoing --variable_cs --cs exp_loss --exp_loss_exponent -1.6 --e $ENERGIES `echo "$perl_opt --i $SYSCONF --append $suffix"`
echo ''

perl ${ACDC_EXECUTABLE} --fortran --save_outgoing --variable_cs --cs exp_loss --exp_loss_exponent -1.6 --e $ENERGIES `echo "$perl_opt --i $SYSCONF --append $suffix"`

echo ''
echo Modifying ACDC files a little to make them compatible with ARCA...

# Modify acdc_equations.f90 file to make it unique for ARCA:
echo ! Input files used to create the current system. Saved in Generic.nc: > ACDC_RECORD_NML
echo \&ACDC_RECORD_NML >> ACDC_RECORD_NML
echo "  System   = \"$SYSCONF\"" >> ACDC_RECORD_NML
echo "  Energies = \"$ENERGIES\"" >> ACDC_RECORD_NML
echo "  Dipoles  = \"$DIPOLES\"" >> ACDC_RECORD_NML
echo "  Name     = \"$NAME\"" >> ACDC_RECORD_NML
echo '/' >> ACDC_RECORD_NML
echo ! Input files used to create this file: > temp
echo ! System name: $NAME >> temp
echo ! System file: $SYSCONF >> temp
echo ! Energy file: $ENERGIES >> temp
echo ! Dipole file: $DIPOLES >> temp
echo module acdc_equations$suffix >> temp
echo implicit none >> temp && echo private >> temp
echo public ::  feval, jeval, formation$suffix >> temp && echo contains  >> temp
cat acdc_equations$suffix.f90 >> temp && echo end module acdc_equations$suffix >> temp && mv temp acdc_equations$suffix.f90
sed -i "s/feval/feval$suffix/g" acdc_equations$suffix.f90
sed -i "s/jeval/jeval$suffix/g" acdc_equations$suffix.f90
sed -i "s/subroutine formation/subroutine formation$suffix/g" acdc_equations$suffix.f90
sed -i "s/subroutine large_formation/subroutine large_formation$suffix/g" acdc_equations$suffix.f90
sed -i "s/call formation/call formation$suffix/g" acdc_equations$suffix.f90
sed -i "s/use acdc_simulation_setup/use acdc_simulation_setup$suffix/g" acdc_equations$suffix.f90

# Modify acdc_system.f90 file to make it unique for ARCA:
echo ! Input files used to create this file: > temp
echo ! System name: $NAME >> temp
echo ! System file: $SYSCONF >> temp
echo ! Energy file: $ENERGIES >> temp
echo ! Dipole file: $DIPOLES >> temp
cat acdc_system$suffix.f90 >> temp && mv temp acdc_system$suffix.f90

sed -i  "s/module acdc_system/module acdc_system$suffix/g" acdc_system$suffix.f90

echo Finished. From ACDC perl script output above, verify that the file writing was succesful.








#
