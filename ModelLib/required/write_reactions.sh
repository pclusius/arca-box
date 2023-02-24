#/user/bin/bash
head -7 tmp.f90>save_reactions.f90
grep -E 'RCONST\([0-9]*\)' second_Initialize.f90>>save_reactions.f90
grep -E 'RCONST\([0-9]*\)' second_Rates.f90|grep -v !>text_rates.txt
sed -i 's|RCONST|TEXT|g' text_rates.txt
sed -i 's| = | = "|g' text_rates.txt
sed -i 's|)$|)"|g' text_rates.txt

cat text_rates.txt>>save_reactions.f90
tail -7 tmp.f90>>save_reactions.f90

gfortran second_Monitor.f90 save_reactions.f90 -o write_reactions.exe
./write_reactions.exe>reactions.txt
rm second_monitor.mod write_reactions.exe tmp.f90 save_reactions.f90
