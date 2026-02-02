#!/bin/bash

# for loop copying parameter files for each chain
for i in $(seq 2 20); do
    cp params_noZ_chain1 params_noZ_chain${i}
    sed -i "s/chain1/chain${i}/g" params_noZ_chain${i}
    cp params_Z_chain1 params_Z_chain${i}
    sed -i "s/chain1/chain${i}/g" params_Z_chain${i}
done
