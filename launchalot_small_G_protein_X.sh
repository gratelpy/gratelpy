#!/bin/bash
N=30
CLIENT='/usr/users/cbu/waltherg/JIC/Dropbox/GratelPy/src/subclient.py mechanisms/small_G_protein_protein_X_mechanism.txt 15'
for i in `seq $N`; do
    bsub -q normal -o ~/tmp/out-%J.txt $CLIENT
done
