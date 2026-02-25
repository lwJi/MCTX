#!/bin/bash
docker exec $CONTAINERARMLOCAL zsh -c '
  cd /home/joe/Meitner/Cactus &&
  ./simfactory/bin/sim build sim -j8 \
    --machine actions-arm-real64 \
    --thornlist=thornlists/mctx.th &&
  make sim-testsuite
'
