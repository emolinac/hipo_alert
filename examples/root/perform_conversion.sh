#!/bin/bash

# This directory should contain the hipo files that need conversion
dir="/lustre24/expphy/volatile/clas12/rg-l/production/p0v6_calib/calib/recon/022994"

file_name="files.txt"
ls -1 ${dir} > ${file_name}

Temporary folder that will contain the hipo file being converted
temp_hipo_file_directory="../../../service-work-2025/022994/"

for i in $(cat ${file_name}); do
        cp ${dir}/${i} ${temp_hipo_file_directory}

        ./converter_alertbanks.exe ${temp_hipo_file_directory}${i} 0

        rm ${temp_hipo_file_directory}*.hipo
done