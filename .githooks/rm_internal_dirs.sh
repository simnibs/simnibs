#!/bin/sh

set -e


repo_root=$(git rev-parse --show-toplevel)

# For all internal files
for d in $(find $repo_root/**/internal/ -type d)
do
    rm -r "$d"
done
