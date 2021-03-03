#!/bin/bash

perl acdc_2018_05_17.pl --fortran --variable_temp --variable_cs --cs exp_loss --exp_loss_exponent -1.6 --save_outgoing --no_eq 1A --no_eq 1D --e ./Perl_input/dH_dS.txt --i ./Perl_input/input_AD.inp --append _AD_new
