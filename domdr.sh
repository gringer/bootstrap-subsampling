#!/bin/sh

fileName=$1;
/itsshared/phd/scripts/makemdrcfg.pl ${fileName} MAXLOCIVALUE 2 MODELBUILDINTERVAL 1000 COMBOSTART 1 COMBOEND 3 CROSSVALINTERVAL 5 NOTRECOGNIZEDRESPONSE -1 POWERMODEL ON RANDOMSHUFFLE ON RANDOMSEED 2 > ./$(basename ${fileName} .mdr).cfg;
/home/gringer/install/pMDR/pMDR/bin/pMDR ./$(basename ${fileName} .mdr).cfg;
