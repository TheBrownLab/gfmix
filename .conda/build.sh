#!/bin/bash
set -eu -o pipefail
mkdir -p ${PREFIX}/bin
make CC=$CC CPP=$CXX F77=$F77

cp gfmix ${PREFIX}/bin
cp treecns ${PREFIX}/bin
cp rert ${PREFIX}/bin
cp alpha_est_mix_rt ${PREFIX}/bin
cp *dat ${PREFIX}/bin

mkdir -p ${PREFIX}/etc/conda/activate.d

cat > ${PREFIX}/etc/conda/activate.d/${PKG_NAME}_activate.sh <<EOF
#!/bin/bash

GFMIX_DIR=$HOME/.gfmix

if [ -d "$GFMIX_DIR" ]; then rm -Rf $GFMIX_DIR; fi
mkdir $GFMIX_DIR

cp ${PREFIX}/bin/*.dat $GFMIX_DIR
EOF
