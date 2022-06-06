#!/bin/bash
set -eu -o pipefail
mkdir -p ${PREFIX}/bin
make CC=$CC CPP=$CXX F77=$F77

cp gfmix ${PREFIX}/bin
cp treecns ${PREFIX}/bin
cp rert ${PREFIX}/bin
cp alpha_est_mix_rt ${PREFIX}/bin
cp *dat ${PREFIX}/bin

# mkdir -p ${PREFIX}/etc/conda/activate.d
# cat > ${PREFIX}/etc/conda/activate.d/${PKG_NAME}_activate.sh <<EOF
# #!/bin/bash

# # Remove gfmix dir if exits
# if [ -d "$HOME/.gfmix" ]; then rm -Rf $HOME/.gfmix; fi

# # Make gfmix dir and copy .dat file there
# mkdir $HOME/.gfmix
# cp ${PREFIX}/bin/*.dat $HOME/.gfmix

# EOF

# mkdir -p ${PREFIX}/etc/conda/deactivate.d
# cat > ${PREFIX}/etc/conda/deactivate.d/${PKG_NAME}_deactivate.sh <<EOF
# #!/bin/bash

# # Clean up gfmix dir
# if [ -d "$HOME/.gfmix" ]; then rm -Rf $HOME/.gfmix; fi

# EOF