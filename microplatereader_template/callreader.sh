#!/bin/bash

ml numpy
python process_plate6.py $1

if [ -f "20200530_erika_taga_00001_reader_plate_1_kinetic_supp_3_high_F.csv" ]; then

mkdir -p Fluorescence 
cp 20200530_erika_taga_00001_reader_plate_1_kinetic_supp_abs.csv Fluorescence/.
mv 20200530_erika_taga_00001_reader_plate_1_kinetic_supp_3_high_F.csv Fluorescence/20200530_erika_taga_00001_reader_plate_1_kinetic_supp_3_high.csv
cp callanalysis.sh analysis.py user_friendly_beta.py helper_functions.py Manifest.csv Fluorescence/.
fi

if [ -f "20200530_erika_taga_00001_reader_plate_1_kinetic_supp_3_high_L.csv" ]; then

mkdir -p Luminescence 
mv 20200530_erika_taga_00001_reader_plate_1_kinetic_supp_abs.csv Luminescence/.
mv 20200530_erika_taga_00001_reader_plate_1_kinetic_supp_3_high_L.csv Luminescence/20200530_erika_taga_00001_reader_plate_1_kinetic_supp_3_high.csv
cp callanalysis.sh analysis.py user_friendly_beta.py helper_functions.py Manifest.csv Luminescence/.
fi
