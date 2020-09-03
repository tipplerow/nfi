#!/bin/sh
########################################################################
# Usage: allele-footprint-driver.sh [JVM OPTIONS] PROP_FILE1 [PROP_FILE2 ...]
########################################################################

if [ $# -lt 1 ]
then
    echo "Usage:" `basename $0` "[JVM OPTIONS] PROP_FILE1 [PROP_FILE2 ...]"
    exit 1
fi

if [ -z "${NFI_HOME}" ]
then
    echo "Environment variable NFI_HOME is not set; exiting."
    exit 1
fi

${NFI_HOME}/bin/nfi-run.sh nfi.model.AlleleFootprintDriver "$@"
