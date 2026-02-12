#!/bin/bash
DIR=$(cd "$(dirname "$0")"; pwd)
export DYLD_LIBRARY_PATH="$DIR/lib"
"$DIR/bin/gcta64" "$@"
