#!/bin/bash
N=80
CLIENT='/usr/users/cbu/waltherg/JIC/Dropbox/GratelPy/src/subclient.py mechanisms/small_G_protein_protein_X_mechanism.txt 15'
for i in `seq $N`; do
    bsub -q normal -o ~/tmp/out-%J.txt -e ~/tmp/err-%J.txt $CLIENT
done
