#!/bin/bash
N=20
CLIENT='/usr/users/cbu/waltherg/JIC/Projects/my_stuff/networkx/LPA/src/subclient.py mechanisms/small_G_protein_mechanism.txt 11'
for i in `seq $N`; do
    bsub -q normal -o ~/tmp/out-%J.txt $CLIENT
done
