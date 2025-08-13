#!/bin/bash

set -ex

export MCTXSPACE="$PWD"
export WORKSPACE="$PWD/../workspace"
mkdir -p "$WORKSPACE"
cd "$WORKSPACE"

# Check out Cactus
wget https://raw.githubusercontent.com/gridaphobe/CRL/master/GetComponents
chmod a+x GetComponents
./GetComponents --no-parallel --shallow "$MCTXSPACE/scripts/mctx.th"

cd Cactus

# Create a link to the MCTX repository
ln -s "$MCTXSPACE" repos
# Create links for the MCTX thorns
mkdir -p arrangements/MCTX
pushd arrangements/MCTX
ln -s ../../repos/MCTX/* .
popd
