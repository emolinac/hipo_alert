#!/bin/bash

dir="/lustre24/expphy/volatile/clas12/rg-l/production/p0v6_calib/calib/recon/022992"

# ls -1 ${dir} > files.txt

for i in $(cat files.txt); do
        cp ${dir}/${i} ../../../service-work-2025/022992/

        ./converter_alertbanks.exe ../../../service-work-2025/022992/${i} 0

        rm ../../../service-work-2025/022992/*.hipo
done