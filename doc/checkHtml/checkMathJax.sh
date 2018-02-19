#!/bin/bash

DIR=$(dirname $0)

if ! which node &> /dev/null; then
  echo "Warning: Skipping MathJax check. 'node' program not found."
  exit 0
fi

node $DIR/checkMathJax.js "$@"
