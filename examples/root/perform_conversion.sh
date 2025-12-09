#!/bin/bash

dir="/lustre24/expphy/volatile/clas12/rg-l/production/p0v6_calib/calib/recon/022994"

# ls -1 ${dir} > files.txt

for i in $(cat files.txt); do
        cp ${dir}/${i} ../../../service-work-2025/022994/

        ./converter_alertbanks.exe ../../../service-work-2025/022994/${i} 0

        rm ../../../service-work-2025/022994/*.hipo
done