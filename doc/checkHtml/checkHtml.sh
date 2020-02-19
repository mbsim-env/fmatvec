#!/bin/bash

DIR=$(dirname $0)

if ! which java &> /dev/null; then
  echo "Warning: Skipping HTML validation check. 'java' program not found."
  exit 0
fi

if ! java nu.validator.client.SimpleCommandLineValidator --version &> /dev/null; then
  echo "Warning: Skipping HTML validation check. nu-html validator jar-file not found in CLASSPATH."
  exit 0
fi

echo "HTML validation using nu-html validator of files:"
for F in "$@"; do
  echo $F
done

java nu.validator.client.SimpleCommandLineValidator --no-langdetect --skip-non-html "$@"
