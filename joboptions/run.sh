#!/bin/bash

# first MC request
python makeOfficialJobOptionGrid.py ../grids/grid_10_mcrequest00.csv

cd 100xxx
tar -cvf ../joboptions.tar.gz .
cd ..
