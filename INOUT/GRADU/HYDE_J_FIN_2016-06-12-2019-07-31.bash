#!/bin/bash
# This file should reside in INOUT/GRADU
cd ../../
./arcabox.exe INOUT/GRADU/HYDE_2016-06-12/input/HYDE_2016-06-12_J_FIN.conf |tee INOUT/GRADU/HYDE_2016-06-12/J_FIN/runReport.txt
./arcabox.exe INOUT/GRADU/HYDE_2016-10-02/input/HYDE_2016-10-02_J_FIN.conf |tee INOUT/GRADU/HYDE_2016-10-02/J_FIN/runReport.txt
./arcabox.exe INOUT/GRADU/HYDE_2017-05-28/input/HYDE_2017-05-28_J_FIN.conf |tee INOUT/GRADU/HYDE_2017-05-28/J_FIN/runReport.txt
./arcabox.exe INOUT/GRADU/HYDE_2017-06-03/input/HYDE_2017-06-03_J_FIN.conf |tee INOUT/GRADU/HYDE_2017-06-03/J_FIN/runReport.txt
./arcabox.exe INOUT/GRADU/HYDE_2017-07-15/input/HYDE_2017-07-15_J_FIN.conf |tee INOUT/GRADU/HYDE_2017-07-15/J_FIN/runReport.txt
./arcabox.exe INOUT/GRADU/HYDE_2018-05-27/input/HYDE_2018-05-27_J_FIN.conf |tee INOUT/GRADU/HYDE_2018-05-27/J_FIN/runReport.txt
./arcabox.exe INOUT/GRADU/HYDE_2019-07-31/input/HYDE_2019-07-31_J_FIN.conf |tee INOUT/GRADU/HYDE_2019-07-31/J_FIN/runReport.txt
