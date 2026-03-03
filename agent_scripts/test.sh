#!/bin/bash
log=$(mktemp)
trap 'rm -f "$log"' EXIT
if docker exec $CONTAINERARMLOCAL zsh -c '
  cd /home/joe/Meitner/Cactus &&
  make sim-testsuite
' > "$log" 2>&1; then
  if grep -q 'Number failed            -> 0' "$log"; then
    echo "✓ test"
  else
    cat "$log"
    echo "✗ test failed" >&2
    exit 1
  fi
else
  cat "$log"
  echo "✗ test failed" >&2
  exit 1
fi
