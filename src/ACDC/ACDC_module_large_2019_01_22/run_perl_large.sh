#!/bin/bash

# Extended set (up to ~100's-1000's of clusters; loop mode)
./acdc_2019_01_22.pl --fortran --save_outgoing --variable_temp --variable_cs --cs exp_loss --exp_loss_exponent -1.6 --loop --boundary --hs_in --i ./Perl_input/input_ANO_loop.inp --append _large --cs_only 1A,0 --cs_only 1N,0 --cs_only 1O,0