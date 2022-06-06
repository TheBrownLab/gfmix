#!/bin/bash
set -eu -o pipefail
mkdir -p ${PREFIX}/bin
make CC=$CC CPP=$CXX F77=$F77

cp gfmix ${PREFIX}/bin
cp treecns ${PREFIX}/bin
cp rert ${PREFIX}/bin
cp alpha_est_mix_rt ${PREFIX}/bin
cp *dat ${PREFIX}/bin
