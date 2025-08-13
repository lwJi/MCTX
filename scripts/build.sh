#!/bin/bash

set -ex

export MCTXSPACE="$PWD"
export WORKSPACE="$PWD/../workspace"
mkdir -p "$WORKSPACE"
cd "$WORKSPACE"

cd Cactus

# Set up SimFactory
cp "$MCTXSPACE/scripts/actions-$ACCELERATOR-$REAL_PRECISION.cfg" ./simfactory/mdb/optionlists
cp "$MCTXSPACE/scripts/actions-$ACCELERATOR-$REAL_PRECISION.ini" ./simfactory/mdb/machines
cp "$MCTXSPACE/scripts/actions-$ACCELERATOR-$REAL_PRECISION.run" ./simfactory/mdb/runscripts
cp "$MCTXSPACE/scripts/actions-$ACCELERATOR-$REAL_PRECISION.sub" ./simfactory/mdb/submitscripts
cp "$MCTXSPACE/scripts/mctx.th" .
cp "$MCTXSPACE/scripts/defs.local.ini" ./simfactory/etc

# For Formaline
git config --global user.email "mctx@einsteintoolkit.org"
git config --global user.name "Github Actions"

case "$MODE" in
    debug) mode='--debug';;
    optimize) mode='--optimise';;
    *) exit 1;;
esac

# Build
# The build log needs to be stored for later.
time ./simfactory/bin/sim --machine="actions-$ACCELERATOR-$REAL_PRECISION" build "$mode" --jobs $(nproc) sim 2>&1 |
    tee build.log

# Check whether the executable exists and is executable
test -x exe/cactus_sim
