#!/bin/bash
log=$(mktemp)
trap 'rm -f "$log"' EXIT
if docker exec $CONTAINERARMLOCAL zsh -c '
  cd /home/joe/Meitner/Cactus &&
  ./simfactory/bin/sim build sim -j8 \
    --machine actions-arm-real64 \
    --thornlist=thornlists/mctx.th
' > "$log" 2>&1; then
  echo "✓ build"
else
  cat "$log"
  echo "✗ build failed" >&2
  exit 1
fi
