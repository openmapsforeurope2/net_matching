#!/bin/bash
# /usr/local/bin/env_setup.sh
export LD_LIBRARY_PATH=/usr/lib/x86_64-linux-gnu:$LD_LIBRARY_PATH
exec "$@"