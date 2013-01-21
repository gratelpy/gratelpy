#!/bin/bash
N=80
CLIENT='/usr/users/cbu/waltherg/JIC/Dropbox/GratelPy/src/subclient.py mechanisms/double_layer_mapk_mechanism.txt 12'
for i in `seq $N`; do
    bsub -q normal -o ~/tmp/out-%J.txt $CLIENT
done

