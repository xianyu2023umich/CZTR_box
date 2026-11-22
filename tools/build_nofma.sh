#!/bin/sh
# build the no-FMA variant of CZTR_box (bit-identity checks): build-nofma/, cztr_box-nofma.x
cd /Users/xianyu/Documents/FALCON/CZTR_box
F=$(grep '^ARBOR_FFLAGS' $(grep ARBOR_MK config.mk | awk '{print $3}') | sed 's/ARBOR_FFLAGS *:= *//')
make -s BUILD=build-nofma EXE=cztr_box-nofma.x ARBOR_FFLAGS="$F -ffp-contract=off" 2>&1 | grep -E "Error" -A3
exit 0
