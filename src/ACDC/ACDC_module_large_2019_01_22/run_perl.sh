#!/bin/bash

./acdc_2018_08_30.pl --fortran --variable_temp --variable_cs --cs exp_loss --exp_loss_exponent -1.6 --variable_ion_source --save_outgoing --e ./Perl_input/HS298.15K_426clusters2016Apr25.txt --dip ./Perl_input/dip_pol_298.15K_426clusters2016Apr25.txt --i ./Perl_input/input_ANnarrow_neutral_neg_pos.inp --append _AN_ions --no_eq 1A --no_eq 1N
