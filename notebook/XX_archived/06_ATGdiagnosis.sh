#!/bin/bash

#diagnosis of variable 'ATG'
mkdir output/ATGProfile
find /workdir/sh2246/p_evolBNI/output/CDSMSAPerOG_v3/ -type f| parallel -j 20 "python ./src/04_ATGdiagnosis.py {} output/ATGProfile/{/.}.atgProfile.txt"

