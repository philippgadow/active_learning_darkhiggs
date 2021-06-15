#!/bin/bash

DIR=data/xampp_311378_simpleana_100000_monoSbb_zp500_dm200_dh50/
DSID=311378
python compare_cutflows.py $DIR/user.changqia.*.root $DIR/histograms.root $DSID

DIR=data/xampp_311382_simpleana_100024_monoSbb_zp2500_dm200_dh50/
DSID=311382
python compare_cutflows.py $DIR/user.changqia.*.root $DIR/histograms.root $DSID

DIR=data/xampp_312360_simpleana_100004_monoSbb_zp500_dm200_dh130/
DSID=312360
python compare_cutflows.py $DIR/user.changqia.*.root $DIR/histograms.root $DSID

DIR=data/xampp_312364_simpleana_100028_monoSbb_zp2500_dm200_dh130/
DSID=312364
python compare_cutflows.py $DIR/user.changqia.*.root $DIR/histograms.root $DSID
