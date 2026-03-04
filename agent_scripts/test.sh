#!/bin/bash

if [ -z "${CONTAINERARMLOCAL:-}" ]; then
  echo "✗ missing required env var: CONTAINERARMLOCAL" >&2
  exit 1
fi

if [ -z "${CONTAINERLOCALCACTUS:-}" ]; then
  echo "✗ missing required env var: CONTAINERLOCALCACTUS" >&2
  exit 1
fi

log=$(mktemp)
trap 'rm -f "$log"' EXIT
if docker exec \
  -e CONTAINERLOCALCACTUS="$CONTAINERLOCALCACTUS" \
  "$CONTAINERARMLOCAL" zsh -c '
  cd "$CONTAINERLOCALCACTUS" &&
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
