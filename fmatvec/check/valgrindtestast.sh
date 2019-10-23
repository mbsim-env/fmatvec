#! /bin/bash

valgrind --help >& /dev/null
if [ $? -ne 0 ]; then
  echo "WARNING!"
  echo "valgrind not found!"
  echo "Skipping valgrind test of testdump"
  exit 0
else
  if file .libs/lt-testast | grep ELF > /dev/null; then
    valgrind --error-exitcode=200 --num-callers=150 --leak-check=full --show-reachable=yes .libs/lt-testast || exit
  elif file ./testast | grep ELF > /dev/null; then
    valgrind --error-exitcode=200 --num-callers=150 --leak-check=full --show-reachable=yes ./testast || exit
  elif file .libs/testast | grep ELF > /dev/null; then
    valgrind --error-exitcode=200 --num-callers=150 --leak-check=full --show-reachable=yes .libs/testast || exit
  else
    echo "Unknown install/build stage"
    exit 1
  fi
fi
